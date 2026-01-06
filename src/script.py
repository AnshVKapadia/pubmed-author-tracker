import json
import os
import time
from pathlib import Path
from datetime import datetime, timedelta, timezone
from typing import Any, Dict, List, Optional, TextIO, Tuple
import requests
import pandas as pd
import yaml
import xml.etree.ElementTree as ET
from tenacity import retry, stop_after_attempt, wait_exponential, retry_if_exception_type
import re
import unicodedata

import gspread
from google.oauth2.service_account import Credentials

from dotenv import load_dotenv
load_dotenv()

EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

api_key = os.getenv("NCBI_API_KEY")
google_json = os.getenv("GOOGLE_APPLICATION_CREDENTIALS_JSON")

# -------------------- LOCAL LLM (OLLAMA) --------------------
OLLAMA_HOST = os.getenv("OLLAMA_HOST", "http://localhost:11434")
OLLAMA_MODEL = os.getenv("OLLAMA_MODEL", "qwen2.5:1.5b-instruct")
OLLAMA_TIMEOUT = int(os.getenv("OLLAMA_TIMEOUT", "30"))

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

# -------------------- TEXT / NAME HELPERS --------------------

def normalize_name(s: str) -> str:
    if not s:
        return ""
    s = unicodedata.normalize("NFKD", s)
    s = "".join(c for c in s if not unicodedata.combining(c))
    s = s.lower()
    s = re.sub(r"[^a-z\s]", "", s)
    return re.sub(r"\s+", " ", s).strip()

def safe_last_name(full_name: str) -> str:
    parts = normalize_name(full_name).split()
    return parts[-1] if parts else ""

def _text(node: Optional[ET.Element]) -> str:
    return node.text.strip() if node is not None and node.text else ""

def first_n_words(text: str, n: int = 100) -> str:
    words = text.split()
    return " ".join(words[:n])

# -------------------- FILE LOADERS --------------------

def load_yaml(path: str) -> Any:
    """
    Be tolerant to weird encodings (your earlier UnicodeDecodeError).
    Tries UTF-8, UTF-8-SIG, then latin-1 as a last resort.
    """
    for enc in ("utf-8", "utf-8-sig", "latin-1"):
        try:
            with open(path, "r", encoding=enc) as f:
                return yaml.safe_load(f)
        except UnicodeDecodeError:
            continue
    # final fallback with replacement
    with open(path, "r", encoding="utf-8", errors="replace") as f:
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
    names: List[str] = []

    if author_entry.get("full_name"):
        names.append(author_entry["full_name"])

    for n in author_entry.get("other_names", []):
        if n:
            names.append(n)

    # Normalize whitespace + dedupe while preserving order
    seen = set()
    final: List[str] = []
    for n in names:
        key = normalize_name(n)
        if key not in seen:
            seen.add(key)
            final.append(n.strip())

    return final

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

# -------------------- LOCAL LLM: OLLAMA --------------------

def _extract_first_json_object(text: str) -> Optional[dict]:
    """
    Ollama sometimes returns extra text.
    Extract the first {...} JSON object we can parse.
    """
    if not text:
        return None
    start = text.find("{")
    end = text.rfind("}")
    if start == -1 or end == -1 or end <= start:
        return None
    snippet = text[start:end+1]
    try:
        return json.loads(snippet)
    except Exception:
        return None

@retry(
    reraise=True,
    stop=stop_after_attempt(3),
    wait=wait_exponential(min=1, max=8),
    retry=retry_if_exception_type(requests.RequestException),
)
def ollama_generate(prompt: str) -> str:
    """
    Calls local Ollama server.
    Requires: Ollama running on localhost:11434
    """
    url = f"{OLLAMA_HOST.rstrip('/')}/api/generate"
    payload = {
        "model": OLLAMA_MODEL,
        "prompt": prompt,
        "stream": False,
        "options": {
            "temperature": 0
        }
    }
    r = requests.post(url, json=payload, timeout=OLLAMA_TIMEOUT)
    r.raise_for_status()
    data = r.json()
    return data.get("response", "")

def llm_name_match(
    target_names: List[str],
    candidate_full_names: List[str],
    dbg: Optional[TextIO] = None,
    context_label: str = "",
) -> Dict[str, Any]:
    """
    Ask local LLM if ANY candidate refers to the same person as the target.

    Returns dict:
      {
        "match": bool,
        "best_match": str|None,
        "confidence": float,
        "reason": str
      }
    """

    # Normalize inputs BEFORE LLM
    target_names_norm = [normalize_name(n) for n in target_names]
    candidate_names_norm = [
        normalize_name(c) for c in candidate_full_names if c.strip()
    ]

    # ---------------- HARD DETERMINISTIC MATCH ----------------
    target_names_norm = {normalize_name(t) for t in target_names if t}
    candidate_names_norm = {normalize_name(c) for c in candidate_full_names if c}

    exact_matches = target_names_norm & candidate_names_norm
    if exact_matches:
        match = exact_matches.pop()
        result = {
            "match": True,
            "best_match": match,
            "confidence": 1.0,
            "reason": "Exact normalized name match (deterministic)"
        }
        if dbg:
            log_write(dbg, f"[DETERMINISTIC MATCH] {result}")
        return result

    payload = {
        "target_names": list(target_names),
        "candidate_names": list(candidate_full_names),
    }

    prompt = f"""
    You are an identity matcher for academic author names.

    TASK:
    Determine whether ANY candidate name refers to the SAME PERSON as ANY target name variant.

    IMPORTANT CONTEXT:
    - Target names are ALIASES of the SAME person.
    - Candidates are names found on a paper.
    - A SINGLE good match is sufficient for match=true.

    CRITICAL RULES (READ CAREFULLY):
    1. First name + last name alignment is REQUIRED for any match.
    2. Middle-name handling:
    - A middle INITIAL (e.g., "T") MAY correspond to a full middle name (e.g., "Talha").
    - Treat initial ↔ expanded-name matches as POSSIBLE, not exact.
    - Assign moderate confidence (≈0.5–0.75) for these cases.
    3. Missing middle names are acceptable and should LOWER confidence, not force NO.
    4. Conflicting middle initials (e.g., target "T", candidate "M") → NOT a match.
    5. Initial-based formats (e.g., "Ahmed ST", "Ahmed S") are valid aliases.
    6. Do NOT require exact string matches.
    7. Do NOT assume name uniqueness.
    8. Compare target ↔ candidate ONE PAIR AT A TIME.
    9. If ANY pair reasonably matches → match=true.

    CONFIDENCE GUIDELINES:
    - Exact full-name match → 0.9–1.0
    - Initial ↔ expanded middle name → 0.5–0.75
    - Initial-only vs full name → 0.4–0.6
    - Different first OR last name → 0.0

    INPUT JSON:
    {json.dumps(payload, indent=2)}

    OUTPUT JSON ONLY (no markdown, no text):
    {{
    "match": true|false,
    "best_match": "<target ↔ candidate pair or null>",
    "confidence": 0.0,
    "reason": "brief explanation"
    }}
    """.strip()

    if dbg:
        log_section(dbg, f"LLM NAME MATCH {context_label}".strip())
        #log_write(dbg, prompt)
        log_write(dbg, f"{json.dumps(payload, indent=2)}")

    raw = ollama_generate(prompt)
    out = _extract_first_json_object(raw) or {}

    # Harden outputs
    match = bool(out.get("match", False))
    best_match = out.get("best_match")
    confidence = out.get("confidence", 0.0)
    reason = out.get("reason", "")

    if best_match is not None and not isinstance(best_match, str):
        best_match = None

    try:
        confidence = float(confidence)
    except Exception:
        confidence = 0.0

    result = {
        "match": match,
        "best_match": best_match,
        "confidence": max(0.0, min(1.0, confidence)),
        "reason": reason if isinstance(reason, str) else str(reason)
    }

    if dbg:
        #log_write(dbg, f"\n[LLM RAW]\n{raw}\n")
        log_write(dbg, f"[LLM PARSED]\n{json.dumps(result, indent=2)}\n")

    return result

# -------------------- EFETCH + LLM MATCH --------------------

def _collect_people(article: ET.Element, tag: str) -> List[Dict[str, str]]:
    """
    tag: "Author" or "Investigator"
    Returns list of dicts with fore, last, initials, display_name.
    """
    out = []
    for p in article.findall(f".//{tag}"):
        last = _text(p.find("LastName"))
        fore = _text(p.find("ForeName"))
        init = _text(p.find("Initials"))

        # display name: prefer ForeName + LastName; fallback to Initials + LastName
        if fore and last:
            display = f"{fore} {last}".strip()
        elif init and last:
            display = f"{init} {last}".strip()
        else:
            display = " ".join([x for x in [fore, init, last] if x]).strip()

        if display:
            out.append({
                "fore": fore,
                "last": last,
                "initials": init,
                "display_name": display
            })
    return out

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
        return []

    affiliations_keywords = author_entry.get("affiliations", [])
    target_full_name = author_entry.get("full_name", "")
    target_last = safe_last_name(target_full_name)

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

            # Collect author + investigator candidates
            authors = _collect_people(article, "Author")
            investigators = _collect_people(article, "Investigator")

            author_display_names = [p["display_name"] for p in authors]
            collab_display_names = [p["display_name"] for p in investigators]

            # --- last-name prefilter to keep LLM tiny ---
            author_candidates = [
                n for n in author_display_names
                if safe_last_name(n) == target_last and target_last != ""
            ]
            collab_candidates = [
                n for n in collab_display_names
                if safe_last_name(n) == target_last and target_last != ""
            ]

            # --- LLM match ---
            author_match = False
            collab_match = False
            llm_info: Optional[Dict[str, Any]] = None
            llm_kind = ""

            if author_candidates:
                llm_info = llm_name_match(
                    target_names=search_names,
                    candidate_full_names=author_candidates,
                    dbg=dbg,
                    context_label=f"(AUTHOR) PMID {pmid}"
                )
                author_match = bool(llm_info.get("match", False))
                llm_kind = "AUTHOR"
            elif collab_candidates:
                llm_info = llm_name_match(
                    target_names=search_names,
                    candidate_full_names=collab_candidates,
                    dbg=dbg,
                    context_label=f"(COLLAB) PMID {pmid}"
                )
                collab_match = bool(llm_info.get("match", False))
                llm_kind = "COLLAB"

            # If last name exists in paper but LLM says no match -> treat as mismatch
            # If neither list had candidates -> name not present at all
            status_parts = []

            if author_candidates:
                if author_match:
                    status_parts.append("AUTHOR MATCH (LLM)")
                else:
                    status_parts.append("AUTHOR MISMATCH (LLM)")
                    # skip adding the paper if you only want true matches:
                    #continue
            elif collab_candidates:
                if collab_match:
                    status_parts.append("COLLABORATOR MATCH (LLM)")
                else:
                    status_parts.append("COLLABORATOR MISMATCH (LLM)")
                    #continue
            else:
                status_parts.append("TARGET LAST NAME NOT FOUND")

            # --- affiliations on the article (not per-author reliably) ---
            affs = [
                (aff.text or "").strip().replace(",", "")
                for aff in article.findall(".//Affiliation")
                if aff.text
            ]

            affiliation_match = any(
                kw.lower() in aff.lower()
                for kw in affiliations_keywords
                for aff in affs
            )

            if affs:
                status_parts.append("AFF MATCH" if affiliation_match else "AFF MISMATCH")
            else:
                status_parts.append("AFF NOT FOUND")

            # DOI
            doi = next(
                (
                    _text(a)
                    for a in article.findall(".//ArticleId")
                    if a.attrib.get("IdType", "").lower() == "doi"
                ),
                ""
            )

            # pub year
            pub_year = _text(article.find(".//PubDate/Year"))

            # attach brief LLM info
            llm_conf = ""
            llm_best = ""
            if llm_info:
                llm_conf = str(llm_info.get("confidence", ""))
                llm_best = str(llm_info.get("best_match", ""))

            status = " + ".join(status_parts)

            rows.append({
                "pmid": pmid,
                "title": _text(article.find(".//ArticleTitle")),
                "journal": _text(article.find(".//Journal/Title")),
                "pub_year": pub_year,
                "doi": doi,
                "authors": "; ".join(author_display_names),
                "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                "status": status,
                "llm_kind": llm_kind,
                "llm_best_match": llm_best,
                "llm_confidence": llm_conf,
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
    authors = load_yaml("config/authors2.yaml")
    settings = load_yaml("config/settings.yaml")

    state = load_state("data/state.json")
    now = datetime.now(timezone.utc)

    # -------------------- DATE WINDOW LOGIC --------------------
    settings_start = (settings.get("starting_date") or "").strip()
    settings_end = (settings.get("ending_date") or "").strip()

    if settings_start:
        mindate = settings_start
        maxdate = settings_end if settings_end else ymd(now)
        update_state = False
        seen = set()
    else:
        last_run = state.get("last_run_utc")
        seen = set(state.get("seen_pmids", []))
        if last_run:
            mindate = ymd(iso_to_utc_dt(last_run))
        else:
            mindate = ymd(now - timedelta(days=30))
        maxdate = ymd(now)
        update_state = True

    print(f"[DEBUG] Date window: {mindate} → {maxdate}")
    print(f"[DEBUG] Local LLM: {OLLAMA_MODEL} at {OLLAMA_HOST}")

    all_rows = []
    dbg = open_debug_log()

    for a in authors:
        full_name = a["full_name"]
        search_names = get_author_search_names(a)

        log_section(dbg, f"AUTHOR IDENTITY: {full_name}")
        log_write(dbg, "Search name variants:")
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

        rows = efetch_details(
            new_pmids,
            a,
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
        "llm_kind",
        "llm_best_match",
        "llm_confidence",
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
            ws.update(range_name="A1", values=[["No publications found in this date window."]])

    write_meta_to_worksheet(sh)

    if update_state:
        state["last_run_utc"] = now.isoformat()
        state["seen_pmids"] = sorted(seen)
        save_state("data/state.json", state)

    print("Run complete: Google Sheet updated.")

if __name__ == "__main__":
    main()
