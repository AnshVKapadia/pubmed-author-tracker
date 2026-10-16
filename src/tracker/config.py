"""Validated configuration and separate collection/reporting date windows."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
import re
import unicodedata

RESERVED_TABS = {"config", "master", "low confidence articles", "history", "summary",
                 "quarterly summary", "faculty summary", "history backup"}


def normalize(value: str) -> str:
    value = unicodedata.normalize("NFKD", value).casefold()
    value = "".join(c for c in value if not unicodedata.combining(c))
    return " ".join(re.sub(r"[^a-z0-9 ]", " ", value).split())


@dataclass(frozen=True)
class Author:
    full_name: str
    other_names: tuple[str, ...] = ()
    affiliations: tuple[str, ...] = ()
    orcid: str = ""

    @property
    def names(self) -> tuple[str, ...]:
        return tuple(dict.fromkeys((self.full_name, *self.other_names)))


def pipe_list(value: str) -> tuple[str, ...]:
    return tuple(x.strip() for x in value.split("|") if x.strip())


def parse_config(values: list[list[str]]) -> tuple[dict[str, str], list[Author]]:
    settings, authors = {}, []
    mode = None
    for row in values:
        cells = [str(c).strip() for c in row] + [""] * 4
        key, value = cells[:2]
        if not key or key.startswith("#"):
            continue
        if key in ("[SETTINGS]", "[AUTHORS]"):
            mode = key
            continue
        if mode == "[SETTINGS]":
            if key.casefold() == "key":
                continue
            if key in settings:
                raise ValueError(f"Duplicate Config setting: {key}")
            settings[key] = value
        elif mode == "[AUTHORS]":
            if key.casefold() == "full_name":
                continue
            if key in {"last_run_utc", "starting_date", "ending_date", "ncbi_email"}:
                raise ValueError(f"Move {key} into Config's [SETTINGS] section.")
            authors.append(Author(key, pipe_list(value), pipe_list(cells[2]), cells[3]))
        else:
            raise ValueError("Config must begin with [SETTINGS] or [AUTHORS].")
    validate_authors(authors)
    for key in ("starting_date", "ending_date", "history_start_date"):
        if settings.get(key):
            parse_date(settings[key])
    if settings.get("last_run_utc"):
        parse_timestamp(settings["last_run_utc"])
    theme_limit = settings.get("theme_max_papers") or "100"
    if not theme_limit.isdigit() or not 1 <= int(theme_limit) <= 500:
        raise ValueError("theme_max_papers must be blank or an integer between 1 and 500.")
    if settings.get("themes_enabled", "false").casefold() not in {"true", "false"}:
        raise ValueError("themes_enabled must be true or false.")
    return settings, authors


def validate_authors(authors: list[Author]) -> None:
    if not authors:
        raise ValueError("No faculty found under Config's [AUTHORS] section.")
    seen = set()
    for author in authors:
        name = author.full_name
        if not normalize(name) or name.casefold() in RESERVED_TABS:
            raise ValueError(f"Invalid or reserved faculty tab name: {name}")
        if len(name) > 100 or re.search(r"[\[\]:*?/\\]", name) or name.startswith("'") or name.endswith("'"):
            raise ValueError(f"Faculty name cannot be used as a worksheet title: {name}")
        if normalize(name) in seen:
            raise ValueError(f"Duplicate faculty name: {name}")
        seen.add(normalize(name))
        if any('"' in n or '[' in n or ']' in n for n in author.names):
            raise ValueError(f"Name aliases must contain names only: {name}")
        if any('"' in n or '[' in n or ']' in n for n in author.affiliations):
            raise ValueError(f"Institution keywords must be plain text: {name}")
        if author.orcid and not re.fullmatch(r"(?:https://orcid.org/)?\d{4}-\d{4}-\d{4}-\d{3}[\dX]", author.orcid):
            raise ValueError(f"Invalid ORCID for {name}")


def parse_date(value: str) -> date:
    try:
        return date.fromisoformat(value.replace("/", "-"))
    except ValueError:
        raise ValueError("Dates must use YYYY-MM-DD or YYYY/MM/DD.") from None


def parse_timestamp(value: str) -> datetime:
    try:
        result = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError:
        raise ValueError("last_run_utc must be an ISO timestamp with a time zone.") from None
    if result.tzinfo is None:
        raise ValueError("last_run_utc must include its time zone (for example +00:00).")
    return result.astimezone(timezone.utc)


def quarter_window(label: str) -> tuple[date, date]:
    match = re.fullmatch(r"(\d{4})-Q([1-4])", label)
    if not match:
        raise ValueError("Report quarter must use YYYY-Q1, YYYY-Q2, YYYY-Q3, or YYYY-Q4.")
    year, quarter = map(int, match.groups())
    start = date(year, quarter * 3 - 2, 1)
    next_start = date(year + 1, 1, 1) if quarter == 4 else date(year, quarter * 3 + 1, 1)
    return start, next_start - timedelta(days=1)


def previous_quarter(today: date) -> str:
    last_day = date(today.year, ((today.month - 1) // 3) * 3 + 1, 1) - timedelta(days=1)
    return f"{last_day.year}-Q{(last_day.month - 1) // 3 + 1}"


def collection_window(settings: dict, today: date, start: str | None = None,
                      end: str | None = None) -> tuple[date, date, bool]:
    manual_start = start or settings.get("starting_date")
    manual_end = end or settings.get("ending_date")
    if manual_end and not manual_start:
        raise ValueError("ending_date/--end requires starting_date/--start.")
    if manual_start:
        first, last = parse_date(manual_start), parse_date(manual_end) if manual_end else today
        checkpoint = False
    else:
        last = today
        if settings.get("last_run_utc"):
            checkpoint_day = parse_timestamp(settings["last_run_utc"]).date()
            if checkpoint_day > today:
                raise ValueError("last_run_utc is in the future; correct Config before running.")
            first = checkpoint_day - timedelta(days=91)
        else:
            first = quarter_window(previous_quarter(today))[0]
        if settings.get("history_start_date"):
            first = min(first, parse_date(settings["history_start_date"]))
        checkpoint = True
    if first > last or last > today:
        raise ValueError("Collection start must be on/before end, and end cannot be in the future.")
    return first, last, checkpoint
