"""Build a fresh Master for one publication-date window."""

from __future__ import annotations

from collections import Counter
from datetime import date
import calendar
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


def date_membership(row: dict, start: date, end: date) -> str:
    value = str(row.get("pub_date", ""))
    try:
        if len(value) == 10:
            first = last = date.fromisoformat(value)
        elif len(value) == 7:
            first = date.fromisoformat(value + "-01")
            last = date(first.year, first.month, calendar.monthrange(first.year, first.month)[1])
        elif len(value) == 4 and value.isdigit():
            first, last = date(int(value), 1, 1), date(int(value), 12, 31)
        else:
            return "uncertain"
    except ValueError:
        raise ValueError("Publication dates must use YYYY, YYYY-MM, or YYYY-MM-DD.") from None
    if last < start or first > end:
        return "outside"
    return "inside" if start <= first and last <= end else "uncertain"


def build_reports(matches: list[dict], authors: list[Author], period: str, timestamp: str,
                  end_date: date | None = None) -> dict:
    # --quarter remains a shortcut; all outputs use this same selected window.
    start, end = quarter_window(period) if end_date is None else (date.fromisoformat(period), end_date)
    roster = {normalize(a.full_name) for a in authors}
    selected = []
    uncertain_dates = set()
    for original in matches:
        if normalize(original["tracked_author"]) not in roster:
            continue
        membership = date_membership(original, start, end)
        if membership == "outside":
            continue
        row = dict(original)
        if membership == "uncertain":
            uncertain_dates.add(row["pmid"])
            row["status"] += " + NEEDS REVIEW (publication date not precise enough)"
        selected.append(row)
    approved = [r for r in selected if accepted(r) and date_membership(r, start, end) == "inside"]
    review = [r for r in selected if not accepted(r) or date_membership(r, start, end) != "inside"]
    analytics = analyze(approved, authors)
    analytics.update(report_start=start.isoformat(), report_end=end.isoformat(),
                     period_key=f"{start.isoformat()}/{end.isoformat()}", generated_at_utc=timestamp,
                     unassigned_date_articles=len(uncertain_dates),
                     review_matches=len(review))
    tables = {"Master": (BASE_COLUMNS + ["pub_date"], publication_rows(approved))}
    return {"analytics": analytics, "tables": tables,
            "review": review,
            "papers_for_themes": publication_rows(approved)}


def write_reports(report: dict, output: Path, snapshot: dict | None = None) -> None:
    output.mkdir(parents=True, exist_ok=True)
    columns, rows = report["tables"]["Master"]
    for filename, fields, data in (("master", columns, rows), ("review", HISTORY_COLUMNS, report["review"])):
        with (output / f"{filename}.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(data)
    (output / "analytics.json").write_text(json.dumps(report["analytics"], indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    if snapshot is not None:
        (output / "sheets-backup.json").write_text(json.dumps(snapshot, ensure_ascii=False), encoding="utf-8")
