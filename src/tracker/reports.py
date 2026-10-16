"""Derive publication-level and faculty-level reports from cumulative history."""

from __future__ import annotations

from collections import Counter
from datetime import date
import csv
import json
from pathlib import Path
from statistics import median

from .config import Author, normalize, quarter_window
from .matching import BASE_COLUMNS, HISTORY_COLUMNS, accepted


def publication_rows(matches: list[dict]) -> list[dict]:
    grouped = {}
    for row in matches:
        grouped.setdefault(str(row["pmid"]), []).append(row)
    rows = []
    for pmid in sorted(grouped, key=int, reverse=True):
        members = grouped[pmid]
        item = {column: members[0].get(column, "") for column in BASE_COLUMNS + ["pub_date"]}
        item["tracked_author"] = "; ".join(sorted({m["tracked_author"] for m in members}))
        item["other_institutions"] = int(any(str(m.get("other_institutions", "0")) == "1" for m in members))
        rows.append(item)
    return rows


def quarter_for(row: dict) -> str:
    value = str(row.get("pub_date", ""))
    if len(value) < 7:
        return ""
    try:
        parsed = date.fromisoformat(value if len(value) == 10 else value + "-01")
    except ValueError:
        raise ValueError("History pub_date must be YYYY, YYYY-MM, or YYYY-MM-DD.") from None
    return f"{parsed.year}-Q{(parsed.month - 1) // 3 + 1}"


def analyze(matches: list[dict], authors: list[Author]) -> dict:
    publications = publication_rows(matches)
    counts = {a.full_name: len({str(r["pmid"]) for r in matches if normalize(r["tracked_author"]) == normalize(a.full_name)})
              for a in authors}
    values = list(counts.values())
    types = Counter(t for p in publications for t in set(p.get("article_types", "").split("; ")) if t)
    return {"total_articles": len(publications), "faculty_count": len(authors),
            "per_faculty": {"median": median(values) if values else 0,
                            "range": {"min": min(values, default=0), "max": max(values, default=0)}},
            "faculty_productivity_percentages": {
                f">{n}_publication" + ("s" if n != 1 else ""): round(sum(v > n for v in values) / len(values) * 100, 2) if values else 0
                for n in (1, 3, 5)},
            "authorship": {"first_author_publications": len({r["pmid"] for r in matches if str(r.get("first_author")) == "1"}),
                           "last_author_publications": len({r["pmid"] for r in matches if str(r.get("last_author")) == "1"})},
            "collaboration": {"articles_with_other_institutions": sum(p["other_institutions"] for p in publications)},
            "article_types": dict(types.most_common()), "faculty_counts": counts}


def build_reports(history: list[dict], authors: list[Author], report_quarter: str, timestamp: str) -> dict:
    start, end = quarter_window(report_quarter)
    roster = {normalize(a.full_name) for a in authors}
    selected = [r for r in history if normalize(r["tracked_author"]) in roster]
    approved = [r for r in selected if accepted(r)]
    current = [r for r in approved if quarter_for(r) == report_quarter]
    analytics = analyze(current, authors)
    analytics.update(report_quarter=report_quarter, report_start=start.isoformat(), report_end=end.isoformat(),
                     generated_at_utc=timestamp, period_months=3,
                     cumulative_total_articles=len(publication_rows(approved)),
                     unassigned_date_articles=len({r["pmid"] for r in approved if not quarter_for(r)}),
                     review_matches=sum(not accepted(r) and r.get("decision") != "reject" for r in selected))
    quarters = sorted({quarter_for(r) for r in approved if quarter_for(r)} | {report_quarter})
    quarterly = []
    for label in quarters:
        summary = analyze([r for r in approved if quarter_for(r) == label], authors)
        quarterly.append({"quarter": label, "total_articles": summary["total_articles"],
                          "faculty_with_publications": sum(v > 0 for v in summary["faculty_counts"].values()),
                          "faculty_count": len(authors), "median_per_faculty": summary["per_faculty"]["median"],
                          **summary["authorship"]})
    tables = {"History": (HISTORY_COLUMNS, history),
              "Master": (BASE_COLUMNS + ["pub_date"], publication_rows(approved)),
              "Low Confidence Articles": (HISTORY_COLUMNS, [r for r in selected if not accepted(r)]),
              "Quarterly Summary": (["quarter", "total_articles", "faculty_with_publications", "faculty_count",
                                     "median_per_faculty", "first_author_publications", "last_author_publications"], quarterly),
              "Faculty Summary": (["tracked_author", "report_quarter", "publications", "cumulative_publications"], [
                  {"tracked_author": a.full_name, "report_quarter": report_quarter,
                   "publications": analytics["faculty_counts"][a.full_name],
                   "cumulative_publications": len({r["pmid"] for r in approved if normalize(r["tracked_author"]) == normalize(a.full_name)})}
                  for a in authors])}
    for author in authors:
        tables[author.full_name] = (HISTORY_COLUMNS, [r for r in approved if normalize(r["tracked_author"]) == normalize(author.full_name)])
    return {"analytics": analytics, "tables": tables, "papers_for_themes": publication_rows(current)}


def summary_table(analytics: dict) -> tuple[list[str], list[dict]]:
    keys = ["schema_version", "run_status", "report_quarter", "report_start", "report_end", "generated_at_utc",
            "total_articles", "cumulative_total_articles", "review_matches", "analytics_json"]
    values = {**analytics, "schema_version": 1, "run_status": "complete",
              "analytics_json": json.dumps(analytics, ensure_ascii=False, allow_nan=False, separators=(",", ":"))}
    return ["key", "value"], [{"key": key, "value": values[key]} for key in keys]


def write_reports(report: dict, output: Path, snapshot: dict | None = None) -> None:
    output.mkdir(parents=True, exist_ok=True)
    names = {"History": "history", "Master": "master", "Low Confidence Articles": "review",
             "Quarterly Summary": "quarterly", "Faculty Summary": "faculty", "Summary": "summary",
             "History Backup": "history-before"}
    for title, filename in names.items():
        columns, rows = report["tables"][title]
        with (output / f"{filename}.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    (output / "analytics.json").write_text(json.dumps(report["analytics"], indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    if snapshot is not None:
        (output / "sheets-backup.json").write_text(json.dumps(snapshot, ensure_ascii=False), encoding="utf-8")
