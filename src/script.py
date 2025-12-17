import json
import os
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List, Optional
import requests
import pandas as pd
import yaml
import xml.etree.ElementTree as ET
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type

import gspread
from google.oauth2.service_account import Credentials

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

from dotenv import load_dotenv
load_dotenv()
api_key = os.getenv("NCBI_API_KEY")
google_json = os.getenv("GOOGLE_APPLICATION_CREDENTIALS_JSON")
if not api_key:
    print("Warning: NCBI_API_KEY not set; PubMed requests will use lower rate limits.")
if not google_json:
    raise RuntimeError(
        "Missing GOOGLE_APPLICATION_CREDENTIALS_JSON. "
        "Set it in your local .env (for local runs) or as a GitHub Actions secret."
    )


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


# -------------------- PUBMED QUERY --------------------

def build_query(author_name: str, affiliation: Optional[str], mindate: str, maxdate: str) -> str:
    parts = [f'"{author_name}"[Author]']
    if affiliation and affiliation.strip():
        parts.append(f'"{affiliation.strip()}"[Affiliation]')
    parts.append(f'("{mindate}"[Date - Publication] : "{maxdate}"[Date - Publication])')
    return " AND ".join(parts)


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


def efetch_details(pmids: List[str], tool, email, api_key) -> List[Dict[str, Any]]:
    if not pmids:
        return []

    rows = []
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

        xml_text = http_get(f"{EUTILS_BASE}/efetch.fcgi", params).text
        root = ET.fromstring(xml_text)

        for article in root.findall(".//PubmedArticle"):
            pmid = _text(article.find(".//PMID"))
            title = _text(article.find(".//ArticleTitle"))
            journal = _text(article.find(".//Journal/Title"))
            pub_year = _text(article.find(".//PubDate/Year"))

            authors = []
            for a in article.findall(".//Author"):
                last = _text(a.find("LastName"))
                init = _text(a.find("Initials"))
                if last:
                    authors.append(f"{last} {init}".strip())

            doi = ""
            for aid in article.findall(".//ArticleId"):
                if aid.attrib.get("IdType") == "doi":
                    doi = _text(aid)

            rows.append({
                "pmid": pmid,
                "title": title,
                "journal": journal,
                "pub_year": pub_year,
                "doi": doi,
                "authors": "; ".join(authors),
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
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

    ws.freeze(rows=1)
    ws.set_basic_filter()


# -------------------- MAIN --------------------

def main():
    authors = load_yaml("config/authors.yaml")
    settings = load_yaml("config/settings.yaml")

    state = load_state("data/state.json")
    last_run = iso_to_utc_dt(state["last_run_utc"])
    now = datetime.now(timezone.utc)

    mindate, maxdate = ymd(last_run), ymd(now)

    seen = set(state["seen_pmids"])
    all_rows = []

    for a in authors:
        pmids = esearch_pmids(
            build_query(a["name"], a.get("affiliation"), mindate, maxdate),
            settings.get("ncbi_tool"),
            settings.get("ncbi_email"),
            api_key,
        )

        new_pmids = [p for p in pmids if p not in seen]
        details = efetch_details(new_pmids,
                                 settings.get("ncbi_tool"),
                                 settings.get("ncbi_email"),
                                 api_key)

        for d in details:
            d["tracked_author"] = a["name"]
            all_rows.append(d)

        seen.update(new_pmids)

    df = pd.DataFrame(all_rows)
    os.makedirs("out", exist_ok=True)
    df.to_csv("out/master.csv", index=False)

    # ---- WRITE TO GOOGLE SHEETS ----
    gc = connect_gsheets()
    sh = gc.open_by_key(settings["spreadsheet_id"])

    write_df_to_worksheet(sh, "Master", df)

    for author in df["tracked_author"].unique():
        write_df_to_worksheet(
            sh,
            author,
            df[df["tracked_author"] == author]
        )

    state["last_run_utc"] = now.isoformat()
    state["seen_pmids"] = sorted(seen)
    save_state("data/state.json", state)

    print("Run complete: Google Sheet updated.")


if __name__ == "__main__":
    main()
