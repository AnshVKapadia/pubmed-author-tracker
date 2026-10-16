"""Snapshot reads and atomic replacement of the current Master report."""

from __future__ import annotations

import copy
import json
import time

from .config import parse_config

OUTPUT_TITLES = {"History", "History Backup", "Master", "Low Confidence Articles", "Quarterly Summary", "Faculty Summary", "Summary"}


def records(values: list[list], title: str) -> list[dict]:
    if not values:
        return []
    headers = [str(v).strip() for v in values[0]]
    if len(headers) != len(set(headers)):
        raise ValueError(f"Duplicate column headers in {title}; fix them before publishing.")
    result = []
    for row in values[1:]:
        if any(str(v).strip() for v in row):
            if len(row) > len(headers) or any(not h and i < len(row) and str(row[i]).strip() for i, h in enumerate(headers)):
                raise ValueError(f"Unlabelled data in {title}; no data will be overwritten.")
            result.append({h: str(row[i]) if i < len(row) else "" for i, h in enumerate(headers) if h})
    return result


def quote_title(title: str) -> str:
    return "'" + title.replace("'", "''") + "'"


def checkpoint_request(values: list[list], timestamp: str, sheet_id: int) -> dict:
    in_settings = False
    fallback = None
    row_index = None
    for index, row in enumerate(values):
        key = str(row[0]).strip() if row else ""
        if key == "[SETTINGS]":
            in_settings = True
            continue
        if key == "[AUTHORS]":
            break
        if in_settings and key == "last_run_utc":
            row_index = index
            break
        if in_settings and not key and not any(str(v).strip() for v in row) and fallback is None:
            fallback = index
    if row_index is None:
        row_index = fallback
    if row_index is None:
        raise ValueError("Add a last_run_utc row (blank value is fine) inside Config's [SETTINGS] block.")
    return {"updateCells": {"range": {"sheetId": sheet_id, "startRowIndex": row_index,
                                      "endRowIndex": row_index + 1, "startColumnIndex": 0, "endColumnIndex": 2},
                            "rows": [{"values": [cell("last_run_utc"), cell(timestamp)]}],
                            "fields": "userEnteredValue"}}


def cell(value) -> dict:
    # stringValue never interprets titles, abstracts, or review notes as formulas.
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return {"userEnteredValue": {"numberValue": value}}
    value = str(value) if value is not None else ""
    if len(value) > 50000:
        raise ValueError("A report cell exceeds Google's 50,000-character limit; reports are saved locally.")
    result = {"userEnteredValue": {"stringValue": value}} if value else {}
    if value.startswith("https://pubmed.ncbi.nlm.nih.gov/"):
        result["userEnteredFormat"] = {"textFormat": {"link": {"uri": value}}}
    return result


def report_setting_requests(values: list[list], serialized: str, sheet_id: int, setting_key: str = "last_report_json") -> list[dict]:
    in_settings = False
    empty = None
    insertion = None
    target = None
    for index, row in enumerate(values):
        key = str(row[0]).strip() if row else ""
        if key == "[SETTINGS]":
            in_settings = True
        elif key == "[AUTHORS]":
            insertion = index
            break
        elif in_settings and key == setting_key:
            target = index
            break
        elif in_settings and not any(str(v).strip() for v in row) and empty is None:
            empty = index
    requests = []
    if target is None:
        target = empty
    if target is None:
        if insertion is None:
            raise ValueError("Config needs [SETTINGS] and [AUTHORS] sections.")
        target = insertion
        values.insert(target, [])
        requests.append({"insertDimension": {"range": {"sheetId": sheet_id, "dimension": "ROWS",
                                                      "startIndex": target, "endIndex": target + 1},
                                              "inheritFromBefore": True}})
    requests.append({"updateCells": {"range": {"sheetId": sheet_id, "startRowIndex": target,
                                                "endRowIndex": target + 1, "startColumnIndex": 0, "endColumnIndex": 2},
                                      "rows": [{"values": [cell(setting_key), cell(serialized)]}],
                                      "fields": "userEnteredValue"}})
    values[target] = [setting_key, serialized, *values[target][2:]]
    return requests


def build_publish_request(snapshot: dict, tables: dict, timestamp: str, advance_checkpoint: bool, report: dict | None = None) -> dict:
    metadata = {s["properties"]["title"]: s["properties"] for s in snapshot["metadata"]["sheets"]}
    next_id = max((v["sheetId"] for v in metadata.values()), default=0) + 1
    # Google sheet IDs are non-negative signed 32-bit integers.
    if next_id + len(tables) >= 2 ** 31:
        next_id = 1
    allocated_ids = {v["sheetId"] for v in metadata.values()}
    requests = []
    config_values = copy.deepcopy(snapshot["tables"]["Config"])
    for title, (columns, rows) in tables.items():
        previous = snapshot["tables"].get(title, [])
        if previous:
            records(previous, title)
            extra_headers = set(str(h) for h in previous[0] if h) - set(columns)
            if extra_headers:
                raise ValueError(f"{title} contains custom columns. Move them to a separate tab or a notes column outside the managed report before publishing.")
            if any(str(value).startswith("=") for row in previous for value in row):
                raise ValueError(f"{title} contains formulas. Move them to a separate tab before publishing.")
        needed_rows, needed_cols = len(rows) + 1, len(columns)
        if title not in metadata:
            while next_id in allocated_ids:
                next_id += 1
            props = {"sheetId": next_id, "title": title,
                     "gridProperties": {"rowCount": max(2, needed_rows), "columnCount": needed_cols, "frozenRowCount": 1}}
            metadata[title] = props
            allocated_ids.add(next_id)
            next_id += 1
            requests.append({"addSheet": {"properties": props}})
        else:
            props = metadata[title]
            grid = props["gridProperties"]
            requests.append({"updateSheetProperties": {"properties": {"sheetId": props["sheetId"], "gridProperties": {
                "rowCount": max(grid["rowCount"], needed_rows, 2), "columnCount": max(grid["columnCount"], needed_cols),
                "frozenRowCount": 1}}, "fields": "gridProperties(rowCount,columnCount,frozenRowCount)"}})
        sheet_id = metadata[title]["sheetId"]
        values = [columns] + [[row.get(column, "") for column in columns] for row in rows]
        values += [[""] * needed_cols for _ in range(max(0, len(previous) - len(values)))]
        # Clear trailing userEnteredValues in the same atomic operation. Clear
        # obsolete explicit links too, so a reordered PMID cannot keep an old URL.
        requests.append({"updateCells": {"range": {"sheetId": sheet_id, "startRowIndex": 0,
                                                    "endRowIndex": max(needed_rows, len(previous)),
                                                    "startColumnIndex": 0, "endColumnIndex": needed_cols},
                                         "rows": [{"values": [cell(v) for v in row]} for row in values],
                                         "fields": "userEnteredValue,userEnteredFormat.textFormat.link"}})
        requests.append({"setBasicFilter": {"filter": {"range": {"sheetId": sheet_id, "startRowIndex": 0,
                                                                  "endRowIndex": needed_rows, "startColumnIndex": 0,
                                                                  "endColumnIndex": needed_cols}}}})
    if report is not None:
        _, authors = parse_config(snapshot["tables"]["Config"])
        obsolete = (OUTPUT_TITLES - {"Master"}) | {a.full_name for a in authors}
        known_headers = {"quarter", "total_articles", "faculty_with_publications", "faculty_count",
                         "median_per_faculty", "first_author_publications", "last_author_publications",
                         "key", "value", "report_quarter", "publications", "cumulative_publications"}
        from .matching import HISTORY_COLUMNS
        known_headers.update(HISTORY_COLUMNS)
        for title in sorted(obsolete & metadata.keys()):
            old_values = snapshot["tables"].get(title, [])
            records(old_values, title)  # Reject unlabelled data before deleting a generated legacy tab.
            if old_values and (set(str(c) for c in old_values[0] if c) - known_headers or
                               any(str(c).startswith("=") for row in old_values for c in row)):
                raise ValueError(f"{title} has custom content. Move it to a separate tab before removing legacy reports.")
            requests.append({"deleteSheet": {"sheetId": metadata[title]["sheetId"]}})
        serialized = json.dumps(report, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
        requests.extend(report_setting_requests(config_values, serialized, metadata["Config"]["sheetId"]))
        settings, _ = parse_config(config_values)
        if not settings.get("report_interval_months"):
            requests.extend(report_setting_requests(config_values, "3", metadata["Config"]["sheetId"], "report_interval_months"))
    if advance_checkpoint:
        requests.append(checkpoint_request(config_values, timestamp, metadata["Config"]["sheetId"]))
    body = {"requests": requests}
    if len(json.dumps(body).encode("utf-8")) > 8 * 1024 * 1024:
        raise ValueError("Atomic publication exceeds 8 MiB. Local reports are saved; reduce output before publishing.")
    return body


class SheetsStore:
    def __init__(self, spreadsheet_id: str, credentials_json: str):
        if not spreadsheet_id:
            raise ValueError("Set SPREADSHEET_ID to the target Google spreadsheet ID.")
        if not credentials_json:
            raise ValueError("Set GOOGLE_APPLICATION_CREDENTIALS_JSON in .env or GitHub Actions secrets.")
        import gspread
        from google.oauth2.service_account import Credentials
        try:
            credentials = Credentials.from_service_account_info(json.loads(credentials_json),
                scopes=["https://www.googleapis.com/auth/spreadsheets"])
        except (ValueError, KeyError, TypeError):
            raise ValueError("GOOGLE_APPLICATION_CREDENTIALS_JSON is not valid service-account JSON.") from None
        self.client = gspread.authorize(credentials)
        self.sh = self.read_retry(lambda: self.client.open_by_key(spreadsheet_id))

    @staticmethod
    def read_retry(operation):
        import gspread
        import requests
        for attempt in range(5):
            try:
                return operation()
            except gspread.exceptions.APIError as error:
                status = error.response.status_code
                if status not in {429, 500, 502, 503, 504} or attempt == 4:
                    raise RuntimeError(f"Google Sheets request failed (HTTP {status}). Check API access and spreadsheet sharing.") from None
            except (requests.ConnectionError, requests.Timeout):
                if attempt == 4:
                    raise RuntimeError("Google Sheets connection failed after five attempts.") from None
            except (PermissionError, gspread.exceptions.SpreadsheetNotFound):
                raise RuntimeError("Cannot access the spreadsheet. Enable Sheets API and share the spreadsheet with the service account as Editor.") from None
            time.sleep(2 ** attempt)

    def read(self) -> dict:
        metadata = self.read_retry(self.sh.fetch_sheet_metadata)
        properties = {s["properties"]["title"]: s["properties"] for s in metadata["sheets"]}
        if "Config" not in properties:
            raise ValueError("The spreadsheet has no Config tab. See examples/Config.tsv.")
        config = self.read_retry(lambda: self.sh.values_get(quote_title("Config"))).get("values", [])
        _, authors = parse_config(config)
        owned = OUTPUT_TITLES | {a.full_name for a in authors}
        titles = sorted(owned & properties.keys())
        # One batch read avoids dozens of per-faculty requests and quota bursts.
        data = self.read_retry(lambda: self.sh.values_batch_get([quote_title(t) for t in titles],
                              params={"valueRenderOption": "FORMULA"})) if titles else {"valueRanges": []}
        tables = {t: item.get("values", []) for t, item in zip(titles, data["valueRanges"])}
        tables["Config"] = config
        return {"metadata": metadata, "tables": tables}

    def publish(self, snapshot: dict, tables: dict, timestamp: str, advance_checkpoint: bool, report: dict) -> None:
        import gspread
        import requests
        body = build_publish_request(copy.deepcopy(snapshot), tables, timestamp, advance_checkpoint, report)
        # Recheck human-editable input immediately before the atomic mutation.
        for title in ("Config", "Master"):
            if title in snapshot["tables"]:
                current = self.read_retry(lambda title=title: self.sh.values_get(quote_title(title),
                                          params={"valueRenderOption": "FORMULA"})).get("values", [])
                if current != snapshot["tables"][title]:
                    raise RuntimeError(f"{title} changed during collection. Rerun to preserve the latest edits.")
        expected = json.dumps(report, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
        for attempt in range(5):
            try:
                self.sh.batch_update(body)
                break
            except gspread.exceptions.APIError as error:
                status = error.response.status_code
                if status == 429 and attempt < 4:
                    time.sleep(2 ** attempt)
                    continue
                if status < 500:
                    raise RuntimeError(f"Sheets rejected the atomic publish (HTTP {status}); no partial batch was applied.") from None
                # A server error may follow an applied request; verify instead of
                # retrying addSheet operations and risking duplicate-sheet errors.
                if not self._published(expected):
                    raise RuntimeError("Sheets publication could not be confirmed. Inspect Config.last_report_json before retrying; local reports are saved.") from None
                return
            except (requests.ConnectionError, requests.Timeout):
                if not self._published(expected):
                    raise RuntimeError("Sheets connection interrupted during publication. Inspect Config.last_report_json before retrying.") from None
                return
        if not self._published(expected):
            raise RuntimeError("Sheets publish verification failed; inspect Config.last_report_json and rerun if necessary.")

    def _published(self, expected: str) -> bool:
        try:
            values = self.read_retry(lambda: self.sh.values_get(quote_title("Config"))).get("values", [])
            return any(len(row) > 1 and row[0] == "last_report_json" and row[1] == expected for row in values)
        except RuntimeError:
            return False
