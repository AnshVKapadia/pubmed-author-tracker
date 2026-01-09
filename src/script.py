import json
import os
import time
from pathlib import Path
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List, Optional, TextIO
from collections import Counter
import requests
import pandas as pd
import yaml
import xml.etree.ElementTree as ET
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type
import unicodedata
import re

import gspread
from google.oauth2.service_account import Credentials

from dotenv import load_dotenv
load_dotenv()
from threading import Lock

from google import genai
from pydantic import BaseModel, Field

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

api_key = os.getenv("NCBI_API_KEY")
google_json = os.getenv("GOOGLE_APPLICATION_CREDENTIALS_JSON")
llm_key = os.getenv("GEMINI_API_KEY")

if not api_key:
    print("Warning: NCBI_API_KEY not set; PubMed requests will use lower rate limits.")
if not google_json:
    raise RuntimeError(
        "Missing GOOGLE_APPLICATION_CREDENTIALS_JSON. "
        "Set it in your local .env or as a GitHub Actions secret."
    )

# ----------------- GEMINI LOGIC + HELPERS  ---------------------
client = genai.Client(api_key=llm_key)

def build_theme_prompt(papers: list[dict]) -> str:
        return f"""
    You are a research analyst.

    TASK:
    Identify the top 5 research themes represented across the provided articles.

    WHAT IS A THEME:
    A recurring scholarly focus such as:
    - topic or subject area
    - methodology or study design
    - population or sample type
    - intervention, exposure, or technique
    - outcome, objective, or research question

    RULES:
    - Use ONLY the provided titles and abstracts.
    - Do NOT infer information that is not explicitly present.
    - Do NOT assume author intent, institutions, or fields.
    - Group closely related concepts into a single theme.
    - Rank themes by number of articles contributing (rank 1 = most frequent).

    ARTICLES:
    {json.dumps(papers, indent=2)}
    """.strip()

def extract_json_llm(text: str) -> dict:
    start = text.find("{")
    end = text.rfind("}")
    if start == -1 or end == -1:
        raise ValueError("No JSON found")
    return json.loads(text[start:end+1])

class Theme(BaseModel):
    rank: int = Field(ge=1, le=5)
    name: str
    article_count: int = Field(ge=0)
    description: str

class ThemesOut(BaseModel):
    themes: list[Theme] = Field(min_length=5, max_length=5)

def gemini_extract_themes(papers: list[dict]) -> dict:
    prompt = build_theme_prompt(papers)

    resp = client.models.generate_content(
        model="gemini-2.5-flash",
        contents=prompt,
        config={
            "temperature": 0.2,
            "response_mime_type": "application/json",
            "response_json_schema": ThemesOut.model_json_schema(),
        },
    )
    parsed = ThemesOut.model_validate_json(resp.text)
    return parsed.model_dump()


# ----------------- LOGGING HELPERS --------------------
def ensure_logs_dir() -> None:
    Path("logs").mkdir(parents=True, exist_ok=True)

def open_debug_log() -> TextIO:
    ensure_logs_dir()
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%SZ")
    return open(f"logs/debug_{ts}.txt", "w", encoding="utf-8")

def log_write(dbg: TextIO, msg: str) -> None:
    dbg.write(msg + "\n")

def log_section(dbg: TextIO, title: str):
    log_write(dbg, "=" * 60)
    log_write(dbg, title)
    log_write(dbg, "=" * 60)

# -------------------- SHEETS SAFETY --------------------
SHEETS_MIN_INTERVAL = 1.2  # seconds (~50 writes/min)
_last_sheets_write = 0.0

def sheets_guard():
    global _last_sheets_write
    now = time.monotonic()
    delta = now - _last_sheets_write
    if delta < SHEETS_MIN_INTERVAL:
        time.sleep(SHEETS_MIN_INTERVAL - delta)
    _last_sheets_write = time.monotonic()

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
def get_author_search_names(author_entry: Dict[str, Any]) -> List[str]:
    """
    Use ONLY explicitly declared name variants from authors.yaml.
    No auto-generation.
    """

    names = []

    if author_entry.get("full_name"):
        names.append(author_entry["full_name"])

    for n in author_entry.get("other_names", []):
        if n:
            names.append(n)

    # Normalize whitespace + dedupe while preserving order
    seen = set()
    final = []
    for n in names:
        key = normalize_name(n)
        if key not in seen:
            seen.add(key)
            final.append(n.strip())

    return final

def normalize_name(s: str) -> str:
    if not s:
        return ""
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.lower()
    s = re.sub(r"[^a-z\s]", "", s)
    return re.sub(r"\s+", " ", s).strip()

def first_n_words(text: str, n: int = 75) -> str:
    """
    Return the first n words of text.
    Safe for None / empty input.
    """
    if not text:
        return ""
    words = text.split()
    return " ".join(words[:n])

def build_query(author_name: str, mindate: str, maxdate: str) -> str:
    return (
        f'"{author_name}" '
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

def normalize_article_types(pub_types: List[str]) -> Dict[str, bool]:
    pt = {p.lower() for p in pub_types}

    return {
        "is_review": "review" in pt,
        "is_clinical_trial": any("clinical trial" in p for p in pt),
        "is_case_report": "case reports" in pt,
        "is_meta_analysis": "meta-analysis" in pt,
        "is_guideline": "practice guideline" in pt,
        "is_editorial": "editorial" in pt,
        "is_letter": "letter" in pt,
        "is_journal_article": "journal article" in pt,
    }

def extract_abstract(article: ET.Element) -> str:
    parts = []
    for a in article.findall(".//AbstractText"):
        if a.text:
            parts.append(a.text.strip())
    return " ".join(parts)

# -------------------- EFETCH + FILTER --------------------
def efetch_details(
    pmids: List[str],
    author_entry: Dict[str, Any],
    search_names: List[str],
    tool: Optional[str],
    email: Optional[str],
    api_key: Optional[str],
    dbg: Optional[TextIO] = None,
) -> List[Dict[str, Any]]:

    if not pmids:
        return [], 0, 0
    
    first_author_count = 0
    last_author_count = 0
 
    affiliations = author_entry.get("affiliations", [])

    normalized_search_names = {
        normalize_name(n) for n in search_names
    }

    rows: List[Dict[str, Any]] = []
    url = f"{EUTILS_BASE}/efetch.fcgi"

    for i in range(0, len(pmids), 200):
        chunk = pmids[i:i + 200]

        params = {
            "db": "pubmed",
            "id": ",".join(chunk),
            "retmode": "xml",
        }
        if tool:
            params["tool"] = tool
        if email:
            params["email"] = email
        if api_key:
            params["api_key"] = api_key

        xml_text = http_get(url, params).text
        root = ET.fromstring(xml_text)

        for article in root.findall(".//PubmedArticle"):
            pmid = _text(article.find(".//PMID"))

            abstract_full = extract_abstract(article)
            abstract_short = first_n_words(abstract_full, 75)

            # ------------------ ARTICLE TYPE ------------------
            pub_types = [
                _text(pt)
                for pt in article.findall(".//PublicationType")
                if _text(pt)
            ]

            ptypes = normalize_article_types(pub_types)

            # ------------------ AFFILIATIONS ------------------
            affs = [
                (aff.text or "").strip().replace(',','')
                for aff in article.findall(".//Affiliation")
                if aff.text
            ]

            aff_match_tracking = [
                kw.lower() in aff.lower()
                for kw in affiliations
                for aff in affs
            ]

            affiliation_match = any(aff_match_tracking)
            other_affs_exist = not all(aff_match_tracking) and len(affs) > 0

            # ------------------ AUTHOR NAME MATCHING ------------------
            author_elems = article.findall(".//Author")

            author_names = []
            author_match = False
            author_match_index = None
            wrong_author_flagged = False

            for idx, a in enumerate(author_elems):
                last = normalize_name(_text(a.find("LastName")))
                fore = normalize_name(_text(a.find("ForeName")))
                init = normalize_name(_text(a.find("Initials")))    

                together = f"{fore} {last}"
                if together.strip() == "": continue
                author_names.append(together)
                    
                if len(fore) < 2:
                    forelast = f"{last} {fore}"
                elif len(init) > 1:
                    forelast = f"{fore.split(' ')[0]} {" ".join(init[1:])} {last}"
                else:
                    forelast = f"{fore} {last}"

                if forelast in normalized_search_names:
                    #print(f"{forelast} in the list of authors matched with {normalized_search_names}")
                    author_match = True
                    author_match_index = idx
                    wrong_author_flagged = False
                    break
                elif any([last in n for n in normalized_search_names]) and len(last) > 1:
                    wrong_author_flagged = True
                    #print(f"{forelast} wrong flagged the list of authors matched with {normalized_search_names}, last is {last} pmid is {pmid}")
            
            if author_match_index is not None:
                if author_match_index == 0:
                    first_author_count += 1
                elif author_match_index == len(author_elems) - 1:
                    last_author_count += 1

            collab_match = False
            wrong_collab_flagged = False
            if not (author_match or wrong_author_flagged): #If author not found, check investigators
                investigators = article.findall(".//Investigator")
                for a in investigators:
                    last = normalize_name(_text(a.find("LastName")))
                    fore = normalize_name(_text(a.find("ForeName")))
                    init = normalize_name(_text(a.find("Initials")))        

                    together = f"{fore} {last}"
                    if together.strip() == "": continue
                    # fore = normalize_name(together.split(' ')[0]) # Re-split and take fore
                    # last = normalize_name(together.split(' ')[-1]) # Re-split and take last
                        
                    if len(fore) < 2:
                        forelast = f"{last} {fore}"
                    else:
                        forelast = f"{fore} {last}"

                    if normalize_name(forelast) in normalized_search_names:
                        #print(f"{forelast} in the list of collabs matched with {normalized_search_names}")
                        collab_match = True
                        wrong_collab_flagged = False
                        break
                    elif any([normalize_name(last) in n for n in normalized_search_names]) and len(last) > 1:
                        wrong_collab_flagged = True
                        #print(f"{forelast} wrong flagged the list of collabs matched with {normalized_search_names} pmid is {pmid}")

            # --------------- FILTERING LOGIC 2.0 -----------------
            status = []
            if dbg:
                if wrong_author_flagged:
                    status.append("AUTHOR MISMATCH")
                    log_write(dbg, f"[INFO] PMID {pmid}: AUTHOR MISMATCH | tracked={search_names[0]}")
                    log_write(dbg, f"   Authors found: {author_names}")
                    continue
                elif wrong_collab_flagged:
                    status.append("COLLABORATOR MISMATCH")
                    continue
                elif author_match:
                    status.append("AUTHOR MATCH")
                elif collab_match:
                    status.append("COLLABORATOR MATCH")
                else:
                    status.append("NOT FOUND")
                    log_write(dbg, f"[INFO] PMID {pmid}: PERSON NOT FOUND | tracked={search_names[0]}")
                    continue

                if affiliation_match and affs:
                    status.append("AFF MATCH")
                elif affs:
                    status.append("AFF MISMATCH")
                    log_write(dbg, f"[INFO] PMID {pmid}: AFF MISMATCH | tracked={search_names[0]}")
                    log_write(dbg, f"    - {affs[0:3]}")
                    if author_match: continue
                else:
                    status.append("AFF NOT FOUND")
                    log_write(dbg, f"[INFO] PMID {pmid}: AFF NOT FOUND | tracked={search_names[0]}")
                
            status = " + ".join(status)

            # ------------------ METADATA ------------------
            doi = next(
                (
                    _text(a)
                    for a in article.findall(".//ArticleId")
                    if a.attrib.get("IdType", "").lower() == "doi"
                ),
                ""
            )

            rows.append({
                "pmid": pmid,
                "title": _text(article.find(".//ArticleTitle")),
                "journal": _text(article.find(".//Journal/Title")),
                "pub_year": _text(article.find(".//PubDate/Year")),
                "doi": doi,
                "authors": "; ".join(author_names),
                "article_types": "; ".join(pub_types),
                "abstract": abstract_short,
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "status": status,
                "other_institutions": int(other_affs_exist),
            })

    return rows, first_author_count, last_author_count

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

def safe_update(ws, values, range_name=None):
    sheets_guard()
    if range_name:
        ws.update(range_name=range_name, values=values)
    else:
        ws.update(values)

def safe_clear(ws):
    sheets_guard()
    ws.clear()

def safe_add_worksheet(sh, title, rows=1, cols=1):
    sheets_guard()
    return sh.add_worksheet(title=title, rows=rows, cols=cols)

def write_df_to_worksheet(sh, title: str, df: pd.DataFrame):
    try:
        ws = sh.worksheet(title)
    except gspread.WorksheetNotFound:
        ws = safe_add_worksheet(sh, title=title, rows=1, cols=1)

    safe_clear(ws)
    safe_update(ws, [df.columns.tolist()] + df.fillna("").values.tolist())

    if len(df) > 0:
        sheets_guard()
        ws.freeze(rows=1)
        sheets_guard()
        ws.set_basic_filter()

def write_meta_to_worksheet(sh):
    ts = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    try:
        ws = sh.worksheet("Meta")
    except gspread.WorksheetNotFound:
        ws = safe_add_worksheet(sh, title="Meta", rows=10, cols=5)

    safe_update(ws, [["last_updated_utc", ts]], range_name="A1:B1")

    ws.update(range_name="A1:B1", values=[["last_updated_utc", ts]])

# -------------------- MAIN --------------------
def main():
    authors = load_yaml("config/authors.yaml") 
    settings = load_yaml("config/settings.yaml")

    state = load_state("data/state.json")
    now = datetime.now(timezone.utc)

    state = load_state("data/state.json")
    now = datetime.now(timezone.utc)

    # -------------------- DATE WINDOW LOGIC --------------------
    settings_start = (settings.get("starting_date") or "").strip()
    settings_end = (settings.get("ending_date") or "").strip()

    if settings_start:
        # ---- MANUAL OVERRIDE MODE ----
        mindate = settings_start
        maxdate = settings_end if settings_end else ymd(now)
        update_state = False
        seen = set()
    else:
        # ---- INCREMENTAL MODE ----
        last_run = state.get("last_run_utc")
        seen = set(state["seen_pmids"])
        if last_run:
            mindate = ymd(iso_to_utc_dt(last_run))
        else:
            # first-ever run fallback
            mindate = ymd(now - timedelta(days=30))

        maxdate = ymd(now)
        update_state = True
    
    print(f"[DEBUG] Date window: {mindate} → {maxdate}")

    start_dt = datetime.strptime(mindate, "%Y-%m-%d")
    end_dt   = datetime.strptime(maxdate, "%Y/%m/%d")
    period_months = (end_dt.year - start_dt.year) * 12 + (end_dt.month - start_dt.month)

    all_rows = []
    dbg = open_debug_log()

    total_fa, total_la = 0, 0

    for a in authors:
        full_name = a["full_name"]
        search_names = get_author_search_names(a)

        log_section(dbg, f"AUTHOR IDENTITY: {a['full_name']}")
        log_write(dbg, f"Search name variants:")
        for n in search_names:
            log_write(dbg, f"  - {n}")

        pmids = set()

        for name in search_names:
            query = build_query(name, mindate, maxdate)
            print(f"[DEBUG]: Query: {query}")
            result_pmids = esearch_pmids(
                query,
                settings.get("ncbi_tool"),
                settings.get("ncbi_email"),
                api_key=api_key,
            )
            pmids.update(result_pmids)

        new_pmids = list(pmids)

        rows, fa_cnt, la_cnt = efetch_details(
            new_pmids,
            a,                    # ← PASS FULL YAML ENTRY
            search_names,
            settings.get("ncbi_tool"),
            settings.get("ncbi_email"),
            api_key=api_key,
            dbg=dbg,
        )

        total_fa += fa_cnt
        total_la += la_cnt

        for r in rows:
            r["tracked_author"] = full_name
            all_rows.append(r)

        seen.update(new_pmids)

    dbg.close()
    print(f"Debug log written to {dbg.name}")

    EXPECTED_COLS = [
        "pmid",
        "title",
        "journal",
        "pub_year",
        "doi",
        "authors",
        "abstract",
        "article_types",
        "pubmed_url",
        "tracked_author",
        "status",
        "other_institutions",  # ← REQUIRED
    ]

    LOW_CONF_KEYWORDS = ("AFF NOT FOUND", "AFF MISMATCH")

    df = pd.DataFrame(all_rows, columns=EXPECTED_COLS)
    df = df.drop_duplicates(subset=["pmid"])

    low_conf_df = df[df["status"].str.contains("|".join(LOW_CONF_KEYWORDS), na=False)]
    master_df = df[~df.index.isin(low_conf_df.index)]

    total_articles = master_df["pmid"].nunique()

    per_faculty_counts = (
        master_df.groupby("tracked_author")["pmid"]
        .nunique()
    )

    median_per_faculty = per_faculty_counts.median()
    range_per_faculty = (per_faculty_counts.min(), per_faculty_counts.max())

    n_faculty = per_faculty_counts.shape[0]
    pct_gt_1 = (per_faculty_counts > 1).sum() / n_faculty * 100
    pct_gt_3 = (per_faculty_counts > 3).sum() / n_faculty * 100
    pct_gt_5 = (per_faculty_counts > 5).sum() / n_faculty * 100

    article_type_counter = Counter()

    for types in master_df["article_types"]:
        article_type_counter.update(types.split("; "))

    article_types_summary = dict(
        sorted(article_type_counter.items(), key=lambda x: x[1], reverse=True)
    )

    # -------------------- THEME ANALYTICS (LLM) --------------------
    THEME_MAX_PAPERS = 50  # safety cap to limit tokens
    theme_df = df.head(THEME_MAX_PAPERS)
    papers_for_llm = []

    for _, row in theme_df.iterrows():
        papers_for_llm.append({
            "pmid": row["pmid"],
            "title": row["title"],
            "abstract": row["abstract"]
        })

    themes = {}
    if papers_for_llm:
        # try:
            themes = gemini_extract_themes(papers_for_llm)
            print("[DEBUG] Extracted themes via Gemini")
        # except Exception as e:
        #     print("[WARN] Theme extraction failed:", e)


    os.makedirs("out", exist_ok=True)
    df.to_csv("out/master.csv", index=False)

    gc = connect_gsheets()
    sh = gc.open_by_key(os.getenv("SPREADSHEET_ID"))

    write_df_to_worksheet(sh, "Low Confidence Articles", low_conf_df)
    write_df_to_worksheet(sh, "Master", master_df)

    for a in authors:
        name = a["full_name"]
        sub = master_df[master_df["tracked_author"] == name]

        if sub.empty:
            continue  # NO writes for empty authors

        try:
            ws = sh.worksheet(name)
        except gspread.WorksheetNotFound:
            ws = safe_add_worksheet(sh, title=name, rows=1, cols=1)

        safe_clear(ws)
        safe_update(ws, [sub.columns.tolist()] + sub.fillna("").values.tolist())

        sheets_guard()
        ws.freeze(rows=1)
        sheets_guard()
        ws.set_basic_filter()

    write_meta_to_worksheet(sh)

    if update_state:
        state["last_run_utc"] = now.isoformat()
        state["seen_pmids"] = sorted(seen)
        save_state("data/state.json", state)


    print("Run complete: Google Sheet updated.")

    analytics = {
        "period_months": period_months,
        "total_articles": int(total_articles),

        "per_faculty": {
            "median": float(median_per_faculty),
            "range": {
                "min": int(range_per_faculty[0]),
                "max": int(range_per_faculty[1]),
            }
        },

        "faculty_productivity_percentages": {
            ">1_publication": round(pct_gt_1, 2),
            ">3_publications": round(pct_gt_3, 2),
            ">5_publications": round(pct_gt_5, 2),
        },

        "authorship": {
            "first_author_publications": int(total_fa),
            "last_author_publications": int(total_la),
        },

        "collaboration": {
            "articles_with_other_institutions": int(
                master_df["other_institutions"].sum()
            )
        },

        "article_types": article_types_summary,

        "themes": themes
    }

    print("\n===== ANALYTICS SUMMARY =====")
    print(json.dumps(analytics, indent=2))
    print("===== END ANALYTICS =====\n")

if __name__ == "__main__":
    main()