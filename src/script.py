import json
import os
from pathlib import Path
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List, Optional, TextIO
import requests
import pandas as pd
import yaml
import xml.etree.ElementTree as ET
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

import gspread
from google.oauth2.service_account import Credentials

from dotenv import load_dotenv
load_dotenv()

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

api_key = os.getenv("NCBI_API_KEY")
google_json = os.getenv("GOOGLE_APPLICATION_CREDENTIALS_JSON")

if not api_key:
    print("Warning: NCBI_API_KEY not set; PubMed requests will use lower rate limits.")
if not google_json:
    raise RuntimeError(
        "Missing GOOGLE_APPLICATION_CREDENTIALS_JSON. "
        "Set it in your local .env or as a GitHub Actions secret."
    )

# ----------------- LOGGING HELPERS --------------------

def ensure_logs_dir() -> None:
    Path("logs").mkdir(parents=True, exist_ok=True)

def open_debug_log() -> TextIO:
    ensure_logs_dir()
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%SZ")
    return open(f"logs/debug_{ts}.txt", "w", encoding="utf-8")

def log_write(dbg: TextIO, msg: str) -> None:
    dbg.write(msg + "\n")

# -------------------- TIME HELPERS --------------------

def iso_to_utc_dt(iso_str: str) -> datetime:
    return datetime.fromisoformat(iso_str).astimezone(timezone.utc)

def ymd(dt_utc: datetime) -> str:
    return dt_utc.strftime("%Y/%m/%d")

# -------------------- FILE LOADERS --------------------

def load_yaml(path: str) -> Any:
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def load_state(path: str) -> Dict[str, Any]:
    if not os.path.exists(path):
        default_last = (datetime.now(timezone.utc) - timedelta(days=30)).replace(microsecond=0)
        return {"last_run_utc": default_last.isoformat(), "seen_pmids": []}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def save_state(path: str, state: Dict[str, Any]) -> None:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2, sort_keys=True)

# -------------------- PUBMED SEARCH --------------------

def build_query(author_name: str, mindate: str, maxdate: str) -> str:
    return (
        f'"{author_name}"[Author] '
        f'AND ("{mindate}"[PDAT] : "{maxdate}"[PDAT])'
    )

@retry(
    reraise=True,
    stop=stop_after_attempt(4),
    wait=wait_exponential(min=1, max=10),
    retry=retry_if_exception_type(requests.RequestException),
)
def http_get(url: str, params: Dict[str, Any], timeout: int = 30) -> requests.Response:
    r = requests.get(url, params=params, timeout=timeout)
    r.raise_for_status()
    return r

def esearch_pmids(query: str, tool, email, api_key) -> List[str]:
    params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": 500,
        "sort": "pub+date",
    }
    if tool: params["tool"] = tool
    if email: params["email"] = email
    if api_key: params["api_key"] = api_key

    data = http_get(f"{EUTILS_BASE}/esearch.fcgi", params).json()
    return data.get("esearchresult", {}).get("idlist", [])

def _text(node: Optional[ET.Element]) -> str:
    return node.text.strip() if node is not None and node.text else ""

# -------------------- EFETCH + FILTER --------------------

def efetch_details(
    pmids: List[str],
    affiliation_keywords: List[str],
    tracked_author: str,
    tool: Optional[str],
    email: Optional[str],
    api_key: Optional[str],
    dbg: Optional[TextIO] = None,
) -> List[Dict[str, Any]]:

    if not pmids:
        return []

    rows: List[Dict[str, Any]] = []
    url = f"{EUTILS_BASE}/efetch.fcgi"

    for i in range(0, len(pmids), 200):
        chunk = pmids[i:i+200]

        params = {
            "db": "pubmed",
            "id": ",".join(chunk),
            "retmode": "xml",
        }
        if tool: params["tool"] = tool
        if email: params["email"] = email
        if api_key: params["api_key"] = api_key

        xml_text = http_get(url, params).text
        root = ET.fromstring(xml_text)

        for article in root.findall(".//PubmedArticle"):
            pmid = _text(article.find(".//PMID"))

            affiliations = [
                (aff.text or "").strip()
                for aff in article.findall(".//Affiliation")
                if aff.text
            ]

            match = any(
                kw.lower() in aff.lower()
                for kw in affiliation_keywords
                for aff in affiliations
            )

            if not match:
                if dbg:
                    if not affiliations:
                        log_write(dbg, f"[SKIP] PMID {pmid}: NO affiliations | tracked={tracked_author}")
                    else:
                        log_write(dbg, f"[SKIP] PMID {pmid}: affiliation mismatch | tracked={tracked_author}")
                        for aff in affiliations:
                            log_write(dbg, f"   AFF: {aff}")
                continue

            # OPTIONAL DEBUG: investigator but not author
            author_names = [
                f"{_text(a.find('LastName'))} {_text(a.find('Initials'))}".strip()
                for a in article.findall(".//Author")
            ]
            if tracked_author not in author_names and dbg:
                log_write(dbg, f"[INFO] PMID {pmid}: tracked author NOT in AuthorList")

            rows.append({
                "pmid": pmid,
                "title": _text(article.find(".//ArticleTitle")),
                "journal": _text(article.find(".//Journal/Title")),
                "pub_year": _text(article.find(".//PubDate/Year")),
                "doi": next(
                    (_text(a) for a in article.findall(".//ArticleId")
                     if a.attrib.get("IdType", "").lower() == "doi"),
                    ""
                ),
                "authors": "; ".join(author_names),
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            })

    return rows

# -------------------- GOOGLE SHEETS --------------------

def connect_gsheets() -> gspread.Client:
    creds = Credentials.from_service_account_info(
        json.loads(google_json),
        scopes=[
            "https://www.googleapis.com/auth/spreadsheets",
            "https://www.googleapis.com/auth/drive",
        ],
    )
    return gspread.authorize(creds)

def write_df_to_worksheet(sh, title: str, df: pd.DataFrame):
    try:
        ws = sh.worksheet(title)
    except gspread.WorksheetNotFound:
        ws = sh.add_worksheet(title=title, rows=1, cols=1)

    ws.clear()
    ws.update([df.columns.tolist()] + df.fillna("").values.tolist())

    if len(df) > 0:
        ws.freeze(rows=1)
        ws.set_basic_filter()

def write_meta_to_worksheet(sh):
    ts = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    try:
        ws = sh.worksheet("Meta")
    except gspread.WorksheetNotFound:
        ws = sh.add_worksheet(title="Meta", rows=10, cols=5)

    ws.update(range_name="A1:B1", values=[["last_updated_utc", ts]])

# -------------------- MAIN --------------------

def main():
    authors = load_yaml("config/authors.yaml")
    settings = load_yaml("config/settings.yaml")

    state = load_state("data/state.json")
    now = datetime.now(timezone.utc)

    mindate, maxdate = ymd(iso_to_utc_dt(state["last_run_utc"])), ymd(now)

    seen = set(state["seen_pmids"])
    all_rows = []

    dbg = open_debug_log()

    try:
        for a in authors:
            pmids = esearch_pmids(
                build_query(a["name"], mindate, maxdate),
                settings.get("ncbi_tool"),
                settings.get("ncbi_email"),
                api_key,
            )

            new_pmids = [p for p in pmids if p not in seen]

            details = efetch_details(
                new_pmids,
                affiliation_keywords=a.get("affiliation_keywords", []),
                tracked_author=a["name"],
                tool=settings.get("ncbi_tool"),
                email=settings.get("ncbi_email"),
                api_key=api_key,
                dbg=dbg,
            )

            for d in details:
                d["tracked_author"] = a["name"]
                all_rows.append(d)

            seen.update(new_pmids)

    finally:
        dbg.close()
        print(f"Debug log written to {dbg.name}")

        EXPECTED_COLS = [
            "pmid",
            "title",
            "journal",
            "pub_year",
            "doi",
            "authors",
            "pubmed_url",
            "tracked_author",
        ]

        df = pd.DataFrame(all_rows, columns=EXPECTED_COLS)

    os.makedirs("out", exist_ok=True)
    df.to_csv("out/master.csv", index=False)

    gc = connect_gsheets()
    sh = gc.open_by_key(settings["spreadsheet_id"])

    write_df_to_worksheet(sh, "Master", df)

    for a in authors:
        name = a["name"]
        sub = df[df["tracked_author"] == name]

        try:
            ws = sh.worksheet(name)
            ws.clear()
        except gspread.WorksheetNotFound:
            ws = sh.add_worksheet(title=name, rows=5, cols=5)

        if not sub.empty:
            ws.update([sub.columns.tolist()] + sub.fillna("").values.tolist())
            ws.freeze(rows=1)
            ws.set_basic_filter()
        else:
            ws.update(
                range_name="A1",
                values=[["No publications found in this date window."]]
            )


    write_meta_to_worksheet(sh)

    state["last_run_utc"] = now.isoformat()
    state["seen_pmids"] = sorted(seen)
    save_state("data/state.json", state)

    print("Run complete: Google Sheet updated.")

if __name__ == "__main__":
    main()