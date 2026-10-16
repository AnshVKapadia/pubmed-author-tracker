"""Keep publication-faculty relationships and review decisions distinct."""

from __future__ import annotations

from .config import Author, normalize

BASE_COLUMNS = ["pmid", "title", "journal", "pub_year", "doi", "authors", "abstract",
                "article_types", "pubmed_url", "tracked_author", "status", "other_institutions"]
HISTORY_COLUMNS = BASE_COLUMNS + ["pub_date", "author_role", "first_author", "last_author",
                                  "author_affiliations", "decision", "review_notes", "first_seen_utc", "last_seen_utc"]


def internal_affiliation(affiliation: str, keywords: tuple[str, ...]) -> bool:
    haystack = f" {normalize(affiliation)} "
    return any(f" {normalize(word)} " in haystack for word in keywords if normalize(word))


def name_strength(person: dict, author: Author) -> int:
    """2 = full given name/ORCID evidence; 1 = initials only; 0 = no match."""
    if author.orcid and person["orcids"]:
        return 2 if author.orcid.removeprefix("https://orcid.org/") in person["orcids"] else 0
    last, fore = normalize(person["last"]), normalize(person["fore"])
    if not last:
        return 0
    initials = normalize(person["initials"]).replace(" ", "")
    strength = 0
    for raw in author.names:
        alias = normalize(raw)
        if alias.endswith(" " + last):
            given = alias[:-(len(last) + 1)]
        elif alias.startswith(last + " "):
            given = alias[len(last) + 1:]
        else:
            continue
        if not given:
            continue
        words = given.split()
        if fore and len(fore.split()[0]) > 1 and len(words[0]) > 1:
            if words[0] != fore.split()[0]:
                continue
            # Compare middle initials whether they are joined, spaced, or expanded.
            expected_middle = "".join(w if len(w) <= 2 else w[0] for w in words[1:])
            actual_middle = initials[1:] or "".join(w[0] for w in fore.split()[1:])
            if expected_middle and expected_middle != actual_middle:
                continue
            return 2
        expected = "".join(w[0] if len(w) > 2 else w for w in words)
        if initials and (expected == initials or (len(expected) == 1 and initials.startswith(expected))):
            strength = 1
    return strength


def match_article(article: dict, author: Author) -> dict | None:
    candidates = []
    for role, people in (("author", article["byline"]), ("investigator", article["investigators"])):
        for index, person in enumerate(people):
            strength = name_strength(person, author)
            if strength:
                own_match = any(internal_affiliation(a, author.affiliations) for a in person["affiliations"])
                candidates.append((strength, own_match, role == "author", role, index, person))
    if not candidates:
        return None
    # Prioritize supported identity + affiliation; byline breaks equal-confidence ties.
    _, affiliation_match, _, role, index, person = max(candidates, key=lambda item: item[:3])
    strength = name_strength(person, author)
    aff_status = "AFF MATCH" if affiliation_match else ("AFF MISMATCH" if person["affiliations"] else "AFF NOT FOUND")
    identity_status = "AUTHOR MATCH" if role == "author" else "COLLABORATOR MATCH"
    if strength < 2:
        identity_status += " (INITIALS ONLY)"
    all_affiliations = [a for p in article["byline"] + article["investigators"] for a in p["affiliations"]]
    row = {column: article.get(column, "") for column in BASE_COLUMNS}
    row.update(tracked_author=author.full_name, status=f"{identity_status} + {aff_status}",
               pub_date=article["pub_date"], author_role=role,
               first_author=int(role == "author" and index == 0),
               last_author=int(role == "author" and index == len(article["byline"]) - 1),
               author_affiliations="; ".join(person["affiliations"]),
               other_institutions=int(any(not internal_affiliation(a, author.affiliations) for a in all_affiliations)),
               decision="auto", review_notes="", first_seen_utc="", last_seen_utc="")
    return row


def accepted(row: dict) -> bool:
    decision = str(row.get("decision", "auto")).strip().casefold() or "auto"
    if decision not in {"auto", "accept", "reject"}:
        raise ValueError("Review decision must be auto, accept, or reject.")
    return decision == "accept" or (decision == "auto" and "AFF MATCH" in row.get("status", "")
                                    and "INITIALS ONLY" not in row.get("status", "")
                                    and "NEEDS REVIEW" not in row.get("status", ""))


def record_key(row: dict) -> tuple[str, str]:
    return str(row["pmid"]), normalize(row["tracked_author"])


def merge_history(previous: list[dict], updates: list[dict], timestamp: str) -> list[dict]:
    merged = {}
    for row in previous:
        if not str(row.get("pmid", "")).isdigit() or not row.get("tracked_author"):
            raise ValueError("History rows require a numeric PMID and a tracked_author.")
        accepted(row)  # Validate manual decisions even when not in the current roster.
        normalized = {column: row.get(column, "") for column in HISTORY_COLUMNS}
        normalized["decision"] = str(normalized["decision"] or "auto").strip().casefold()
        key = record_key(normalized)
        if key in merged:
            for field in ("decision", "review_notes"):
                left, right = merged[key].get(field), normalized.get(field)
                if left and right and left != right and (field != "decision" or "auto" not in (left, right)):
                    raise ValueError("Conflicting review decisions/notes for one PMID and faculty member.")
            if normalized["decision"] == "auto" and merged[key].get("decision") != "auto":
                normalized["decision"] = merged[key]["decision"]
            normalized["review_notes"] = normalized["review_notes"] or merged[key].get("review_notes", "")
        merged[key] = normalized
    for row in updates:
        key = record_key(row)
        old = merged.get(key, {})
        item = {column: row.get(column, "") for column in HISTORY_COLUMNS}
        item.update(decision=old.get("decision") or "auto", review_notes=old.get("review_notes", ""),
                    first_seen_utc=old.get("first_seen_utc") or timestamp, last_seen_utc=timestamp)
        merged[key] = item
    return sorted(merged.values(), key=lambda r: (str(r.get("pub_date", "")), int(r["pmid"]), r["tracked_author"]), reverse=True)
