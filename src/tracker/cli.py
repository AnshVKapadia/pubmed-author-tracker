"""Compose collection, reporting, backups, and explicit publication."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import json
import logging
import os
from pathlib import Path
import sys

from .config import collection_window, normalize, parse_config, previous_quarter, quarter_window
from .matching import HISTORY_COLUMNS, match_article, merge_history, record_key
from .pubmed import PubMedClient, parse_articles
from .reports import build_reports, summary_table, write_reports
from .sheets import SheetsStore, build_publish_request, history_from_snapshot
from .themes import extract_themes

LOG = logging.getLogger(__name__)
ROOT = Path(__file__).resolve().parents[2]


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description="Track PubMed faculty publications. Dry run unless --publish is supplied.")
    source = result.add_mutually_exclusive_group()
    source.add_argument("--config", type=Path, help="Local Config TSV; no Google credentials needed")
    source.add_argument("--snapshot", type=Path, help="Saved sheet snapshot JSON; read-only")
    result.add_argument("--history", type=Path, help="Previous history CSV for a local run")
    result.add_argument("--pubmed-xml", type=Path, help="Offline PubMed fixture; cannot be published")
    result.add_argument("--start", help="Manual collection start (YYYY-MM-DD); does not advance checkpoint")
    result.add_argument("--end", help="Manual collection end (YYYY-MM-DD)")
    result.add_argument("--quarter", help="Summary quarter (YYYY-Qn); default is the last completed quarter")
    result.add_argument("--output", type=Path, default=ROOT / "out")
    result.add_argument("--check-config", action="store_true", help="Validate configuration and history without querying PubMed")
    mode = result.add_mutually_exclusive_group()
    mode.add_argument("--publish", action="store_true", help="Publish all reports to Sheets atomically")
    mode.add_argument("--dry-run", action="store_true", help="Explicit read-only mode (the default)")
    result.add_argument("--themes", action="store_true", help="Enable optional AI themes for this run")
    return result


def run(args) -> dict:
    if args.publish and (args.snapshot or args.pubmed_xml or args.history):
        raise ValueError("Publishing saved snapshots, fixtures, or local history is disabled. Read live Sheets for publication.")
    if args.publish and args.config:
        raise ValueError("Publishing requires live Config from the target spreadsheet; omit --config.")
    if args.pubmed_xml and not (args.config or args.snapshot):
        raise ValueError("--pubmed-xml requires --config or --snapshot.")
    if args.history and not args.config:
        raise ValueError("--history is supported with --config only; live Sheets supplies its own history.")
    store, snapshot = None, None
    if args.config:
        with args.config.open(encoding="utf-8-sig", newline="") as handle:
            values = list(csv.reader(handle, delimiter="\t"))
    elif args.snapshot:
        snapshot = json.loads(args.snapshot.read_text(encoding="utf-8"))
        values = snapshot["tables"]["Config"]
    else:
        store = SheetsStore(os.getenv("SPREADSHEET_ID", ""), os.getenv("GOOGLE_APPLICATION_CREDENTIALS_JSON", ""))
        snapshot = store.read()
        values = snapshot["tables"]["Config"]
    settings, authors = parse_config(values)
    now = datetime.now(timezone.utc)
    first, last, checkpoint = collection_window(settings, now.date(), args.start, args.end)
    quarter = args.quarter or previous_quarter(now.date())
    quarter_window(quarter)
    history = history_from_snapshot(snapshot, authors) if snapshot else []
    if args.history:
        with args.history.open(encoding="utf-8-sig", newline="") as handle:
            history = list(csv.DictReader(handle))
    history = merge_history(history, [], now.isoformat())
    before_history = history
    LOG.info("Validated %d faculty; %d historical matches. Collection %s to %s; report %s.",
             len(authors), len(history), first, last, quarter)
    if args.check_config:
        return {"status": "configuration_valid", "faculty_count": len(authors), "history_matches": len(history),
                "collection_start": first.isoformat(), "collection_end": last.isoformat(), "report_quarter": quarter}
    # Back up before network collection so even a failed first migration is recoverable.
    args.output.mkdir(parents=True, exist_ok=True)
    if snapshot is not None:
        backup = args.output / "backups" / now.strftime("%Y%m%dT%H%M%S%fZ")
        backup.mkdir(parents=True, exist_ok=True)
        (backup / "sheets-backup.json").write_text(json.dumps(snapshot, ensure_ascii=False), encoding="utf-8")
    if args.pubmed_xml:
        articles = parse_articles(args.pubmed_xml.read_text(encoding="utf-8"))
        candidates = {a.full_name: set(articles) for a in authors}
    else:
        client = PubMedClient(settings.get("ncbi_email", ""), settings.get("ncbi_tool", ""), os.getenv("NCBI_API_KEY", ""))
        candidates = {}
        for index, author in enumerate(authors, 1):
            LOG.info("Searching faculty %d/%d: %s", index, len(authors), author.full_name)
            candidates[author.full_name] = client.search(author, first, last)
        # Refresh missing legacy metadata and revalidate all known historical matches.
        pmids = set().union(*candidates.values()) | {r["pmid"] for r in history}
        LOG.info("Fetching %d unique PMIDs once each.", len(pmids))
        articles = client.fetch(pmids)
    updates = []
    previous = {record_key(row): row for row in history}
    # Examine each fetched record against the entire roster, recovering coauthor
    # relationships that old global PMID deduplication discarded.
    for article in articles.values():
        for author in authors:
            row = match_article(article, author)
            if row:
                updates.append(row)
            elif (article["pmid"], normalize(author.full_name)) in previous:
                row = {**previous[(article["pmid"], normalize(author.full_name))],
                       "status": "NEEDS REVIEW (identity not found in refreshed metadata)", "pub_date": article["pub_date"]}
                updates.append(row)
    history = merge_history(history, updates, now.isoformat())
    report = build_reports(history, authors, quarter, now.isoformat())
    report["tables"]["History Backup"] = (HISTORY_COLUMNS, before_history)
    if args.themes:
        settings["themes_enabled"] = "true"
    report["analytics"]["themes"] = extract_themes(report["papers_for_themes"], settings,
                                                os.getenv("GEMINI_API_KEY", "")) if not args.pubmed_xml else {"status": "offline", "themes": []}
    report["analytics"].update(collection_start=first.isoformat(), collection_end=last.isoformat(),
                               mode="publish" if args.publish else "dry_run",
                               date_policy="Journal issue date; electronic date fallback when journal month is unavailable")
    report["tables"]["Summary"] = summary_table(report["analytics"])
    write_reports(report, args.output, snapshot)
    if snapshot is not None:
        # Validate the complete proposed mutation even during a read-only dry run.
        build_publish_request(snapshot, report["tables"], now.isoformat(), checkpoint)
    if args.publish:
        store.publish(snapshot, report["tables"], now.isoformat(), checkpoint)
        LOG.info("Google Sheets publication verified. %s", "Checkpoint advanced." if checkpoint else "Manual collection: checkpoint preserved.")
    else:
        LOG.info("Dry run complete. Reports: %s. Google Sheets and its checkpoint were not modified.", args.output)
    summary = os.getenv("GITHUB_STEP_SUMMARY")
    if summary:
        analytics = report["analytics"]
        with Path(summary).open("a", encoding="utf-8") as handle:
            handle.write(f"### PubMed {quarter}\n\nMode: {analytics['mode']}\n\n"
                         f"- Quarter publications: {analytics['total_articles']}\n"
                         f"- Cumulative publications: {analytics['cumulative_total_articles']}\n"
                         f"- Faculty: {len(authors)}\n- Pending review: {analytics['review_matches']}\n\n"
                         "Reports are in Google Sheets; History Backup retains the previous history. "
                         "Downloadable artifacts are available only when explicitly enabled.\n")
    return report["analytics"]


def main(argv: list[str] | None = None) -> int:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    from dotenv import load_dotenv
    load_dotenv(ROOT / ".env")
    args = parser().parse_args(argv)
    log_dir = ROOT / "logs"
    log_dir.mkdir(exist_ok=True)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s",
                        handlers=[logging.StreamHandler(), logging.FileHandler(log_dir / "tracker.log", encoding="utf-8")], force=True)
    try:
        analytics = run(args)
        print(json.dumps(analytics, indent=2, ensure_ascii=False, allow_nan=False))
        return 0
    except (ValueError, RuntimeError, OSError, KeyError, TypeError) as error:
        # Only our own errors are intentionally detailed; never dump provider bodies/credentials.
        message = str(error) if isinstance(error, (ValueError, RuntimeError)) else type(error).__name__
        LOG.error("Tracker stopped: %s", message)
        return 1
    except Exception as error:
        LOG.error("Tracker stopped (%s). Check credentials and API availability; reports may be available in out/.", type(error).__name__)
        return 1
