import json
import os
import time
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

import google.generativeai as genai
from dotenv import load_dotenv
load_dotenv()
from threading import Lock

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

api_key = os.getenv("NCBI_API_KEY")
google_json = os.getenv("GOOGLE_APPLICATION_CREDENTIALS_JSON")
genai.configure(api_key=os.getenv("GEMINI_API_KEY"))

GEMINI_MODEL = genai.GenerativeModel(
    model_name="gemma-3-27b-it",
    generation_config={
        "temperature": 0,
        "response_mime_type": "application/json",
    },
)

GEMINI_MIN_INTERVAL = 7.5  # seconds → ~9 requests/min
_last_gemini_ts = 0.0
_gemini_lock = Lock()

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

def log_section(dbg: TextIO, title: str):
    log_write(dbg, "=" * 60)
    log_write(dbg, title)
    log_write(dbg, "=" * 60)


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

def build_search_names_from_author(author_entry: Dict[str, Any]) -> List[str]:
    """
    Build PubMed search names from authors.yaml entry.
    Uses explicit other_names + derived initials.
    """

    base_names = set()

    # Always include full_name
    if author_entry.get("full_name"):
        base_names.add(author_entry["full_name"])

    # Include any manually specified variants
    for n in author_entry.get("other_names", []):
        if n:
            base_names.add(n)

    search_names = set()

    for name in base_names:
        parts = name.strip().split()
        if len(parts) < 2:
            continue

        first = parts[0]
        middle = parts[1:-1]
        last = parts[-1]

        # Full name as-is
        search_names.add(name)

        # First + Last
        search_names.add(f"{first} {last}")

        # Last + First initial
        search_names.add(f"{last} {first[0]}")

        # Last + First+Middle initials
        if middle:
            initials = first[0] + "".join(m[0] for m in middle)
            search_names.add(f"{last} {initials}")

    return sorted(search_names)

import unicodedata
import re

def normalize_name(s: str) -> str:
    if not s:
        return ""
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.lower()
    s = re.sub(r"[^a-z\s]", "", s)
    return re.sub(r"\s+", " ", s).strip()

def first_n_words(text: str, n: int = 100) -> str:
    words = text.split()
    return " ".join(words[:n])

def build_query(author_name: str, mindate: str, maxdate: str) -> str:
    return (
        f'"{author_name}" '
        f'AND ("{mindate}"[PDAT] : "{maxdate}"[PDAT])'
    )

def build_author_profile(author_entry: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert authors.yaml entry into a compact JSON profile for the LLM.
    """
    return {
        "full_name": author_entry.get("full_name"),
        "affiliations": author_entry.get("affiliations", []),
        "notes": author_entry.get("notes", ""),  # optional future field
    }

def build_paper_evidence(article: ET.Element) -> Dict[str, Any]:
    def texts(path):
        return [
            (x.text or "").strip()
            for x in article.findall(path)
            if x.text
        ]
    
    abstract_full = " ".join(texts(".//AbstractText"))

    return {
        "title": _text(article.find(".//ArticleTitle")),
        "journal": _text(article.find(".//Journal/Title")),
        "abstract": first_n_words(abstract_full, 100),
        "authors": [
            {
                "last": _text(a.find("LastName")),
                "fore": _text(a.find("ForeName")),
                "initials": _text(a.find("Initials")),
                "affiliations": texts(".//Affiliation"),
            }
            for a in article.findall(".//Author")
        ],
        "mesh_terms": texts(".//MeshHeading/DescriptorName"),
        "keywords": texts(".//Keyword"),
    }

def gemini_validate_batch(
    author_profile: Dict[str, Any],
    papers: List[Dict[str, Any]],
    dbg
) -> Dict[str, Any]:
    """
    papers = [
      {"pmid": str, "evidence": {...}},
      ...
    ]
    """

    wait_for_gemini_slot()

    prompt = f"""
    You are validating which papers belong to a specific researcher.

    AUTHOR PROFILE:
    {json.dumps(author_profile, indent=2)}

    PAPERS:
    {json.dumps(papers, indent=2)}

    Instructions:
    - Return YES only for papers that clearly belong to this person.
    - Return NO if clearly different.
    - Return UNCERTAIN if unsure.
    - Affiliations matter more than names.
    - Do NOT assume name uniqueness.

    Respond strictly in JSON:
    {{
    "results": {{
        "<PMID>": {{
        "decision": "YES|NO|UNCERTAIN",
        "confidence": 0.0-1.0,
        "reason": "..."
        }}
    }}
    }}
    """

    log_write(dbg, f"\n\n{prompt}\n\n")

    response = GEMINI_MODEL.generate_content(prompt)

    try:
        return json.loads(response.text)
    except Exception:
        return {"results": {}}
    
def wait_for_gemini_slot():
    global _last_gemini_ts
    with _gemini_lock:
        now = time.monotonic()
        delta = now - _last_gemini_ts
        if delta < GEMINI_MIN_INTERVAL:
            time.sleep(GEMINI_MIN_INTERVAL - delta)
        _last_gemini_ts = time.monotonic()

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
    author_entry: Dict[str, Any],   # ← ADD THIS
    search_names: List[str],
    tool: Optional[str],
    email: Optional[str],
    api_key: Optional[str],
    dbg: Optional[TextIO] = None,
) -> List[Dict[str, Any]]:

    if not pmids:
        return []
    
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
            
            # ------------------ AUTHOR NAME MATCHING ------------------
            author_elems = article.findall(".//Author")

            author_names = []
            author_match = False
            wrong_author_flagged = False

            for a in author_elems:
                last = _text(a.find("LastName"))
                fore = _text(a.find("ForeName"))
                init = _text(a.find("Initials"))

                together = f"{fore} {init} {last}"
                if together.strip() == "": continue
                fore = normalize_name(together.split(' ')[0]) # Re-split and take fore
                last = normalize_name(together.split(' ')[-1]) # Re-split and take last
                author_names.append(together)
                    
                forelast = ((fore + " " + last))

                if forelast in normalized_search_names:
                    author_match = True
                    wrong_author_flagged = False
                    break
                elif any([last in n for n in normalized_search_names]):
                    wrong_author_flagged = True

            collab_match = False
            wrong_collab_flagged = False
            if not (author_match or wrong_author_flagged): #If author not found, check investigators
                print(pmid)
                investigators = article.findall(".//Investigator")
                for a in investigators:
                    last = _text(a.find("LastName"))
                    fore = _text(a.find("ForeName"))
                    init = _text(a.find("Initials"))

                    together = f"{fore} {init} {last}"
                    if together.strip() == "": continue
                    fore = normalize_name(together.split(' ')[0]) # Re-split and take fore
                    last = normalize_name(together.split(' ')[-1]) # Re-split and take last
                        
                    forelast = ((fore + " " + last))
                    if pmid == "41364689" and last == "kapadia":
                        print(forelast)
                        print(normalized_search_names)

                    if normalize_name(forelast) in normalized_search_names:
                        collab_match = True
                        wrong_collab_flagged = False
                        break
                    elif any([normalize_name(last) in n for n in normalized_search_names]):
                        wrong_collab_flagged = True

            # ------------------ AFFILIATIONS ------------------
            affs = [
                (aff.text or "").strip().replace(',','')
                for aff in article.findall(".//Affiliation")
                if aff.text
            ]

            affiliation_match = any(
                kw.lower() in aff.lower()
                for kw in affiliations
                for aff in affs
            )

            # --------------- FILTERING LOGIC 2.0 -----------------
            status = []
            if dbg:
                if wrong_author_flagged:
                    status.append("AUTHOR MISMATCH")
                    log_write(dbg, f"[INFO] PMID {pmid}: AUTHOR MISMATCH | tracked={search_names[0]}")
                    log_write(dbg, f"   Authors found: {author_names}")
                    continue
                elif author_match:
                    status.append("AUTHOR MATCH")
                elif wrong_collab_flagged:
                    status.append("COLLABORATOR MISMATCH")
                    continue
                elif collab_match:
                    status.append("COLLABORATOR MATCH")
                else:
                    status.append("AUTHOR NOT FOUND")
                    log_write(dbg, f"[INFO] PMID {pmid}: AUTHOR NOT FOUND | tracked={search_names[0]}")

                if affiliation_match and affs:
                    status.append("AFF MATCH")
                elif affs:
                    status.append("AFF MISMATCH")
                    log_write(dbg, f"[INFO] PMID {pmid}: AFF MISMATCH | tracked={search_names[0]}")
                    log_write(dbg, f"    - {affs[0:3]}")
                else:
                    status.append("AFF NOT FOUND")
                    log_write(dbg, f"[INFO] PMID {pmid}: AFF NOT FOUND | tracked={search_names[0]}")
                
            status = " + ".join(status)

            # ------------------ FILTERING LOGIC ------------------
            # status = "AFFILIATIONS VALID"
            # log_write(dbg, "")
            # if not affiliation_match:
            #     if affs and author_match:
            #         if dbg:
            #             status = "AFFILIATIONS MISMATCH + AUTHOR MISMATCH"
            #             log_write(
            #                 dbg,
            #                 f"[INFO] PMID {pmid}: affiliation mismatch → WRONG person | tracked={search_names[0]}"
            #             )
            #             # for aff in affs:
            #             #     log_write(dbg, f"    - {aff}")
            #             log_write(dbg, f"    - {affs[0]}")
            #             log_write(dbg, f"   Authors found: {author_names}")
            #         # continue
            #     elif affs:
            #         if dbg:
            #             status = "AFFILIATIONS MISMATCH + COLLABORATOR LIKELY"
            #             log_write(
            #                 dbg,
            #                 f"[INFO] PMID {pmid}: collaborator hit (non-author) | tracked={search_names[0]}"
            #             )
            #     else:
            #         status = "AFFILIATIONS NOT FOUND"
            #         if dbg:
            #             log_write(
            #                 dbg,
            #                 f"[INFO] PMID {pmid}: no affiliations found | tracked={search_names[0]}"
            #             )
            #             log_write(dbg, f"   Authors found: {author_names}")
            #             log_write(dbg, f"Do the authors match? {author_match}")

            # ------------- LLM QUERY BUILDING -------------
            # llm_candidates = []

            # use_llm = (
            #     (status == "AFFILIATIONS NOT FOUND" and author_match)
            #     or
            #     (status == "AFFILIATIONS MISMATCH + COLLABORATOR LIKELY")
            # )

            # use_llm = False

            # if use_llm:
            #     log_write(dbg, f"[LLM-CANDIDATE] PMID {pmid}")
            #     log_write(dbg, f"  Status: {status}")
            #     log_write(dbg, f"  Reason: deterministic filters inconclusive")

            #     llm_candidates.append({
            #         "pmid": pmid,
            #         "evidence": build_paper_evidence(article),
            #     })


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
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "status": status,
            })

    # if llm_candidates:
    #     log_section(dbg, f"LLM VALIDATION for {author_entry['full_name']}")

    #     author_profile = build_author_profile(author_entry)
    #     llm_response = gemini_validate_batch(author_profile, llm_candidates, dbg)

    #     for r in rows:
    #         pmid = r["pmid"]
    #         decision = llm_response.get("results", {}).get(pmid)

    #         if decision:
    #             r["status"] = f"LLM_{decision['decision']}"
    #             log_write(dbg, f"PMID {pmid}: {r['status']}")
    #             log_write(dbg, f"  Confidence: {decision['confidence']}")
    #             log_write(dbg, f"  Reason: {decision['reason']}")

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

    
    all_rows = []

    dbg = open_debug_log()

    for a in authors:
        full_name = a["full_name"]
        search_names = build_search_names_from_author(full_name)

        log_section(dbg, f"AUTHOR IDENTITY: {a['full_name']}")
        log_write(dbg, f"Search name variants:")
        for n in search_names:
            log_write(dbg, f"  - {n}")

        pmids = set()

        for name in search_names:
            query = build_query(name, mindate, maxdate)
            #print(f"[DEBUG]: Query: {query}")
            result_pmids = esearch_pmids(
                query,
                settings.get("ncbi_tool"),
                settings.get("ncbi_email"),
                api_key=api_key,
            )
            pmids.update(result_pmids)

        new_pmids = list(pmids)

        rows = efetch_details(
            new_pmids,
            a,                    # ← PASS FULL YAML ENTRY
            search_names,
            settings.get("ncbi_tool"),
            settings.get("ncbi_email"),
            api_key=api_key,
            dbg=dbg,
        )

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
        "pubmed_url",
        "tracked_author",
        "status",
    ]

    df = pd.DataFrame(all_rows, columns=EXPECTED_COLS)

    os.makedirs("out", exist_ok=True)
    df.to_csv("out/master.csv", index=False)

    gc = connect_gsheets()
    sh = gc.open_by_key(os.getenv("SPREADSHEET_ID"))

    write_df_to_worksheet(sh, "Master", df)

    for a in authors:
        name = a["full_name"]
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

    if update_state:
        state["last_run_utc"] = now.isoformat()
        state["seen_pmids"] = sorted(seen)
        save_state("data/state.json", state)


    print("Run complete: Google Sheet updated.")

if __name__ == "__main__":
    main()