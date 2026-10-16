"""Rate-limited NCBI access and XML parsing, independent of Google services."""

from __future__ import annotations

import calendar
from datetime import date, timedelta
import re
import time
import xml.etree.ElementTree as ET

import requests

from .config import Author

BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


class PubMedError(RuntimeError):
    pass


def text(node: ET.Element | None) -> str:
    return "".join(node.itertext()).strip() if node is not None else ""


def build_query(author: Author, start: date, end: date) -> str:
    full_names, initials = [], []
    for name in author.names:
        last_token = name.split()[-1].replace(".", "")
        is_initials = name != author.full_name and last_token.isalpha() and (len(last_token) == 1 or (last_token.isupper() and len(last_token) <= 5))
        (initials if is_initials else full_names).append(name)
    full_query = " OR ".join(f'"{name}"[au] OR "{name}"[ir]' for name in full_names)
    initial_query = " OR ".join(f'"{name}"[au] OR "{name}"[ir]' for name in initials)
    affiliation_query = " OR ".join(f'"{word}"[ad]' for word in author.affiliations)
    # Broad initials identify thousands of unrelated people. Require an article-
    # level institutional clue only for these searches; full names still retrieve
    # missing-affiliation records. Match each person's own affiliation afterward.
    parts = [f"({full_query})"] if full_query else []
    if initial_query and affiliation_query:
        parts.append(f"(({initial_query}) AND ({affiliation_query}))")
    names = " OR ".join(parts)
    if not names:
        raise ValueError(f"Provide a full name or institutional keywords for {author.full_name}.")
    first, last = start.strftime("%Y/%m/%d"), end.strftime("%Y/%m/%d")
    return f'({names}) AND (("{first}"[dp]:"{last}"[dp]) OR ("{first}"[crdt]:"{last}"[crdt]))'


class PubMedClient:
    def __init__(self, email: str, tool: str = "pubmed-author-tracker", api_key: str = "",
                 session=None, sleep=time.sleep, monotonic=time.monotonic):
        if not email or "@" not in email:
            raise ValueError("Set ncbi_email in Config to your contact email before querying PubMed.")
        self.params = {"tool": tool or "pubmed-author-tracker", "email": email, "db": "pubmed"}
        if api_key:
            self.params["api_key"] = api_key
        self.session = session or requests.Session()
        self.sleep, self.monotonic = sleep, monotonic
        self.interval, self.last_request = (0.11 if api_key else 0.36), float("-inf")

    def request(self, endpoint: str, params: dict, *, post: bool = False):
        for attempt in range(5):
            delay = self.interval - (self.monotonic() - self.last_request)
            if delay > 0:
                self.sleep(delay)
            self.last_request = self.monotonic()
            try:
                arguments = {"data" if post else "params": {**self.params, **params}, "timeout": (10, 60)}
                response = (self.session.post if post else self.session.get)(f"{BASE_URL}/{endpoint}.fcgi", **arguments)
                if response.status_code in {429, 500, 502, 503, 504}:
                    if attempt == 4:
                        raise PubMedError(f"NCBI {endpoint} unavailable after five attempts (HTTP {response.status_code}).")
                    retry_after = response.headers.get("Retry-After", "")
                    self.sleep(min(60, max(2 ** attempt, float(retry_after) if retry_after.isdigit() else 0)))
                    continue
                if not 200 <= response.status_code < 300:
                    raise PubMedError(f"NCBI {endpoint} rejected the request (HTTP {response.status_code}).")
                return response
            except (requests.ConnectionError, requests.Timeout):
                if attempt == 4:
                    raise PubMedError(f"NCBI {endpoint} connection failed after five attempts.") from None
                self.sleep(2 ** attempt)
        raise PubMedError("NCBI request did not complete.")

    def search(self, author: Author, start: date, end: date) -> set[str]:
        query = build_query(author, start, end)
        result: set[str] = set()
        total = None
        while total is None or len(result) < total:
            offset = len(result)
            try:
                data = self.request("esearch", {"term": query, "retmode": "json", "retmax": 1000,
                                                "retstart": offset, "sort": "pub date"}).json()
                found = data["esearchresult"]
                if data.get("error") or found.get("errorlist", {}).get("fieldsnotfound"):
                    raise PubMedError("NCBI returned a search error; collection stopped.")
                total = int(found["count"])
                page = found["idlist"]
                if total < 0 or not isinstance(page, list) or any(not str(p).isdigit() for p in page):
                    raise ValueError
            except (ValueError, KeyError, TypeError):
                raise PubMedError("Invalid NCBI search response; existing history has not been replaced.") from None
            if total > 9999:
                if start == end:
                    raise PubMedError(f"More than 9,999 records on one date for {author.full_name}; narrow the search.")
                middle = start + (end - start) // 2
                return self.search(author, start, middle) | self.search(author, middle + timedelta(days=1), end)
            result.update(map(str, page))
            if len(result) == total:
                return result
            if len(result) <= offset or len(result) > total:
                raise PubMedError("NCBI pagination was incomplete or changed mid-search; retry the run.")
        return result

    def fetch(self, pmids: set[str]) -> dict[str, dict]:
        result = {}
        ordered = sorted(pmids, key=int)
        for offset in range(0, len(ordered), 150):
            chunk = ordered[offset:offset + 150]
            response = self.request("efetch", {"id": ",".join(chunk), "retmode": "xml"}, post=True)
            result.update(parse_articles(response.text))
            if set(chunk) - result.keys():
                raise PubMedError("NCBI omitted requested records. No report was published; retry or investigate removed PMIDs.")
        return result


def parse_publication_date(node: ET.Element | None) -> str:
    if node is None:
        return ""
    year = text(node.find("Year"))
    month = text(node.find("Month"))
    day = text(node.find("Day"))
    if not year:
        match = re.search(r"\b(\d{4})\b", text(node.find("MedlineDate")))
        return match[1] if match else ""
    if not re.fullmatch(r"\d{4}", year):
        return ""
    if not month:
        return year
    if not month.isdigit():
        month = next((str(i) for i in range(1, 13) if calendar.month_abbr[i].casefold() == month[:3].casefold()), "")
    if not month or not 1 <= int(month) <= 12:
        return year
    prefix = f"{year}-{int(month):02d}"
    if day.isdigit():
        try:
            return date(int(year), int(month), int(day)).isoformat()
        except ValueError:
            pass
    return prefix


def parse_person(node: ET.Element) -> dict:
    last, fore, initials = (text(node.find(key)) for key in ("LastName", "ForeName", "Initials"))
    return {"last": last, "fore": fore, "initials": initials,
            "name": f"{fore or initials} {last}".strip() or text(node.find("CollectiveName")),
            "affiliations": [text(aff) for aff in node.findall("./AffiliationInfo/Affiliation") if text(aff)],
            "orcids": [text(i).removeprefix("https://orcid.org/") for i in node.findall("Identifier")
                       if i.get("Source", "").casefold() == "orcid"]}


def parse_articles(xml: str) -> dict[str, dict]:
    try:
        root = ET.fromstring(xml)
    except ET.ParseError:
        raise PubMedError("NCBI returned invalid XML; no report was published.") from None
    if root.find(".//ERROR") is not None:
        raise PubMedError("NCBI returned an XML error; no report was published.")
    result = {}
    for article in root.findall(".//PubmedArticle"):
        citation = article.find("MedlineCitation")
        pmid = text(citation.find("PMID")) if citation is not None else ""
        if not pmid.isdigit():
            raise PubMedError("Publication has no valid PMID.")
        byline = [parse_person(a) for a in article.findall("./MedlineCitation/Article/AuthorList/Author")]
        investigators = [parse_person(a) for a in article.findall("./MedlineCitation/InvestigatorList/Investigator")]
        pub_date = parse_publication_date(article.find("./MedlineCitation/Article/Journal/JournalIssue/PubDate"))
        if len(pub_date) < 7:
            electronic = [parse_publication_date(n) for n in article.findall("./MedlineCitation/Article/ArticleDate")
                          if n.get("DateType") == "Electronic"]
            pub_date = next((value for value in electronic if len(value) >= 7), pub_date)
        result[pmid] = {"pmid": pmid, "title": text(article.find(".//ArticleTitle")),
                        "journal": text(article.find(".//Journal/Title")),
                        "pub_year": pub_date[:4], "pub_date": pub_date,
                        "doi": next((text(n) for n in article.findall("./PubmedData/ArticleIdList/ArticleId")
                                     if n.get("IdType") == "doi"), ""),
                        "authors": "; ".join(p["name"] for p in byline),
                        "abstract": " ".join(text(n) for n in article.findall(".//Abstract/AbstractText")),
                        "article_types": "; ".join(text(n) for n in article.findall(".//PublicationType")),
                        "pubmed_url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        "byline": byline, "investigators": investigators}
    return result
