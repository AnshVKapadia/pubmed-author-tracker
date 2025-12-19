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

# -------------------- CONSTANTS --------------------

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# -------------------- ENV --------------------

load_dotenv()
API_KEY = os.getenv("NCBI_API_KEY")
GOOGLE_JSON = os.getenv("GOOGLE_APPLICATION_CREDENTIALS_JSON")

if not GOOGLE_JSON:
    raise RuntimeError("Missing GOOGLE_APPLICATION_CREDENTIALS_JSON")

# -------------------- LOGGING --------------------

def ensure_logs_dir():
    Path("logs").mkdir(exist_ok=True)

def open_debug_log() -> TextIO:
    ensure_logs_dir()
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%SZ")
    return open(f"logs/debug_{ts}.txt", "w", encoding="utf-8")

def log(dbg: Optional[TextIO], msg: str):
    if dbg:
        dbg.write(msg + "\n")

# -------------------- TIME --------------------

def iso_to_utc(iso: str) -> datetime:
    return datetime.fromisoformat(iso).astimezone(timezone.utc)

def ymd(dt: datetime) -> str:
    return dt.strftime("%Y/%m/%d")

# -------------------- FILE LOADERS --------------------

def load_yaml(path: str):
    with open(path, "r", encoding="utf-8") as f:
        return yaml.safe_load(f)

def load_state(path: str) -> Dict[str, Any]:
    if not os.path.exists(path):
        default = datetime.now(timezone.utc) - timedelta(days=30)
        return {"last_run_utc": default.isoformat(), "seen_pmids": []}
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)

def save_state(path: str, state: Dict[str, Any]):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        json.dump(state, f, indent=2)

# -------------------- NAME HELPERS --------------------

def build_search_name(full_name: str) -> str:
    """
    Vishal Kapadia            -> Kapadia V
    Rachel L Leon             -> Leon RL
    Myra H Wyckoff             -> Wyckoff MH
    """
    parts = full_name.strip().split()
    if len(parts) < 2:
        raise ValueError(f"Invalid full_name: {full_name}")

    first = parts[0]
    middle = parts[1:-1]
    last = parts[-1]

    initials = first[0]
    if middle:
        initials += middle[0][0]

    return f"{last} {initials}"

# -------------------- HTTP --------------------

@retry(
    stop=stop_after_attempt(4),
    wait=wait_exponential(min=1, max=10),
    retry=retry_if_exception_type(requests.RequestException),
    reraise=True,
)
def http_get(url: str, params: Dict[str, Any]) -> requests.Response:
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    return r

# -------------------- PUBMED SEARCH --------------------

def build_query(author_search: str, mindate: str, maxdate: str) -> str:
    return f'"{author_search}"[Author] AND ("{mindate}"[EDAT] : "{maxdate}"[EDAT])'

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
    full_name: str,
    affiliations: List[str],
    tool: Optional[str],
    email: Optional[str],
    api_key: Optional[str],
    dbg: Optional[TextIO] = None,
) -> List[Dict[str, Any]]:

    if not pmids:
        return []

    parts = full_name.lower().split()
    tracked_first = parts[0]
    tracked_last = parts[-1]
    tracked_initial = tracked_first[0]

    rows = []
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

        root = ET.fromstring(http_get(url, params).text)

        for article in root.findall(".//PubmedArticle"):
            pmid = _text(article.find(".//PMID"))

            # ---- AFFILIATION FILTER ----
            affs = [
                (a.text or "").strip()
                for a in article.findall(".//Affiliation")
                if a.text
            ]

            aff_match = any(
                kw.lower() in aff.lower()
                for kw in affiliations
                for aff in affs
            )

            if not aff_match:
                log(dbg, f"[SKIP] PMID {pmid}: affiliation mismatch")
                continue

            # ---- AUTHOR NAME FILTER ----
            author_match = False
            author_list = []

            for a in article.findall(".//Author"):
                last = _text(a.find("LastName")).lower()
                fore = _text(a.find("ForeName")).lower()
                init = _text(a.find("Initials")).lower()

                if last:
                    author_list.append(f"{last} {init}")

                if last == tracked_last:
                    if fore.startswith(tracked_first) or init.startswith(tracked_initial):
                        author_match = True

            if not author_match:
                log(dbg, f"[SKIP] PMID {pmid}: affiliation OK but AUTHOR mismatch")
                log(dbg, f"   authors: {author_list}")
                continue

            # ---- METADATA ----
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
                "authors": "; ".join(author_list),
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
            })

    return rows

# -------------------- GOOGLE SHEETS --------------------

def connect_gsheets():
    creds = Credentials.from_service_account_info(
        json.loads(GOOGLE_JSON),
        scopes=[
            "https://www.googleapis.com/auth/spreadsheets",
            "https://www.googleapis.com/auth/drive",
        ],
    )
    return gspread.authorize(creds)

def write_df(sh, title: str, df: pd.DataFrame):
    try:
        ws = sh.worksheet(title)
    except gspread.WorksheetNotFound:
        ws = sh.add_worksheet(title=title, rows=1, cols=1)

    ws.clear()
    ws.update([df.columns.tolist()] + df.fillna("").values.tolist())

    if len(df) > 0:
        ws.freeze(rows=1)
        ws.set_basic_filter()

# -------------------- MAIN --------------------

def main():
    authors = load_yaml("config/authors.yaml")
    settings = load_yaml("config/settings.yaml")
    state = load_state("data/state.json")

    now = datetime.now(timezone.utc)
    mindate, maxdate = ymd(iso_to_utc(state["last_run_utc"])), ymd(now)

    seen = set(state["seen_pmids"])
    all_rows = []

    dbg = open_debug_log()

    for a in authors:
        full_name = a["full_name"]
        search_name = build_search_name(full_name)

        pmids = esearch_pmids(
            build_query(search_name, mindate, maxdate),
            settings.get("ncbi_tool"),
            settings.get("ncbi_email"),
            API_KEY,
        )

        new_pmids = [p for p in pmids if p not in seen]

        rows = efetch_details(
            new_pmids,
            full_name,
            a["affiliations"],
            settings.get("ncbi_tool"),
            settings.get("ncbi_email"),
            API_KEY,
            dbg,
        )

        for r in rows:
            r["tracked_author"] = full_name
            all_rows.append(r)

        seen.update(new_pmids)

    dbg.close()
    print(f"Debug log written to {dbg.name}")

    df = pd.DataFrame(all_rows)
    os.makedirs("out", exist_ok=True)
    df.to_csv("out/master.csv", index=False)

    gc = connect_gsheets()
    sh = gc.open_by_key(settings["spreadsheet_id"])

    write_df(sh, "Master", df)

    for a in authors:
        name = a["full_name"]
        sub = df[df["tracked_author"] == name]
        if not sub.empty:
            write_df(sh, name, sub)

    state["last_run_utc"] = now.isoformat()
    state["seen_pmids"] = sorted(seen)
    save_state("data/state.json", state)

    print("Run complete.")

if __name__ == "__main__":
    main()