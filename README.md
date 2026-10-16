# PubMed Author Tracker

Track a faculty roster in Google Sheets, retain publication history, and generate calendar-quarter summaries. Python collects PubMed records; Google Sheets is the editable roster and reporting interface. Runs are read-only unless `--publish` is explicitly supplied.

## Quick start

Use Python 3.11 or 3.13. From the repository root:

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# macOS/Linux: source .venv/bin/activate
python -m pip install -r requirements.txt -c requirements.lock
python -m unittest discover -s tests -v
```

Try the complete pipeline with synthetic data, without Google credentials or network requests:

```bash
python src/script.py --config examples/Config.tsv --pubmed-xml tests/fixtures/pubmed.xml --quarter 2026-Q3 --dry-run
```

Reports appear in `out/`. The original `python src/script.py` entry point is retained, but now defaults to a dry run.

## Existing Google Sheet setup

Keep the existing `Config` tab. The original three-column author layout is supported. [examples/Config.tsv](examples/Config.tsv) is a copyable template with optional settings and an optional fourth ORCID column.

```text
[SETTINGS]
key                 value
starting_date
ending_date
last_run_utc
history_start_date
ncbi_tool           pubmed-author-tracker
ncbi_email          your-contact-email@example.com
themes_enabled      false
gemini_model
theme_max_papers     100

[AUTHORS]
full_name           other_names                  affiliations            orcid
Jane A Smith        Jane Smith|Smith JA          UT Southwestern|UTSW
```

Use `|` between aliases or institution keywords. Faculty names must be unique and valid worksheet titles. Keep `last_run_utc` inside `[SETTINGS]`, even when its value is blank. Dates use `YYYY-MM-DD` or `YYYY/MM/DD`. The tracker validates these inputs before querying PubMed.

1. Enable **Google Sheets API** in the existing Google Cloud project. The tracker uses only the Sheets scope; Google Drive API is not required for normal runs.
2. Share the spreadsheet with the existing service account's email as **Editor**.
3. Copy `.env.example` to `.env` and set:

   | Environment variable | Purpose |
   | --- | --- |
   | `SPREADSHEET_ID` | ID between `/d/` and `/edit` in the spreadsheet URL |
   | `GOOGLE_APPLICATION_CREDENTIALS_JSON` | Full service-account JSON, stored as one environment variable |
   | `NCBI_API_KEY` | Optional higher NCBI rate limit |
   | `GEMINI_API_KEY` | Optional AI themes |

Do not commit credentials or put them in spreadsheet cells. Existing GitHub Actions secrets with these names are supported; a new Google Cloud project or key is not inherently necessary.

Validate and preview before publishing:

```bash
python src/script.py --check-config
python src/script.py --dry-run
python src/script.py --publish
```

A local TSV can be used with `--config` for collection without Google authentication. Add `--history out/history.csv` to merge a prior local report. Publishing always reads live Config and History from the target spreadsheet; fixtures, saved snapshots, and local history cannot be published.

## History and reports

- **History:** one row per PMID and faculty member, including uncertain matches. This is the durable source of publication-faculty relationships.
- **Master:** accepted publications, one row per PMID. `tracked_author` lists every matched faculty member.
- **Low Confidence Articles:** uncertain or rejected faculty-paper matches, with the evidence status. This is a generated view.
- **Individual faculty tabs:** all accepted historical matches for that faculty member. Empty faculty tabs are refreshed too.
- **Faculty Summary:** quarter counts and cumulative counts for every currently configured faculty member, including zero counts.
- **Quarterly Summary:** unique publication totals and productivity statistics for quarters represented in the retained history.
- **Summary:** the selected quarter and a versioned `analytics_json` value suitable for an email integration.
- **History Backup:** the preceding history, replaced atomically with each successful publication. Use it to recover the previous state if necessary.

Only the current roster contributes to Master and summary metrics; records for removed faculty remain in History. A shared paper counts once in department totals and once for each matching faculty member's productivity. First/last authorship statistics count unique publications with a tracked faculty member in that position; a sole author is both first and last. Investigator credit does not count as byline authorship.

### Review decisions

Edit only the **History** row's `decision` (`auto`, `accept`, or `reject`) and `review_notes` columns. Find the row by both PMID and faculty member. Decisions and notes survive metadata refreshes and reruns. Other report tabs are regenerated from History; edits to their copied review fields do not become source decisions.

`auto` accepts a full-name/ORCID match with that person's own institutional affiliation. Initials-only evidence, absent affiliations, affiliation mismatches, and legacy records awaiting validation remain reviewable. ORCID is used when configured and present in the PubMed record. Rule-based matching still requires human review for ambiguous identities.

Do not rename managed headers or add custom columns to generated tabs. Use a separate tab for formulas and annotations, or History's review notes. Publishing stops on unexpected headers or formulas instead of overwriting them.

### First migration

When History does not exist, the tracker imports records from Master, Low Confidence Articles, and existing faculty tabs. It fetches the associated PubMed metadata again and checks every fetched publication against the current roster, recovering shared-paper relationships discarded by the old global deduplication. Existing records that cannot be revalidated remain in History for review.

A pre-run spreadsheet snapshot is saved locally before collection. Existing retained records are preserved, but papers absent from every old tab cannot be recovered without a backfill. Copy the original workbook before the first production publication; verify a dry run or test copy first.

## Reporting dates and collection dates

The default report is the **last completed calendar quarter**: a run on October 16 reports July 1–September 30. Select another quarter with `--quarter 2026-Q2`.

Collection is independent of reporting:

- An incremental collection searches from 91 days before `last_run_utc` through today. With no checkpoint, it starts at the previous quarter's beginning.
- It searches **publication date OR PubMed create date**, helping catch older papers added to PubMed recently.
- Explicit full-name searches include missing-affiliation records. Initials-only aliases require an article-level institutional clue to avoid thousands of unrelated candidates. The subsequent matcher still checks the tracked person's own affiliation. An initials-only record without any institutional clue may therefore be missed.
- Known historical PMIDs are fetched again and refreshed. Corrections to unknown older records outside the overlap window may require a backfill.
- `history_start_date`, when populated, expands every incremental search back to that date. Leave it blank for routine incremental runs; use a manual backfill for a one-time recovery.

Backfill example:

```bash
python src/script.py --start 2024-01-01 --end 2026-09-30 --quarter 2026-Q3 --dry-run
# After inspection, repeat with --publish.
```

Manual `--start`/`--end` or Config `starting_date`/`ending_date` overrides merge into history but do not advance `last_run_utc`. Clear Config overrides after use. Summary dates use the journal issue date, with electronic publication date as a fallback if no journal month is available. Year-only/seasonal dates remain in cumulative history but are not assigned to a quarter; their count is disclosed.

## Reliable publishing and recovery

The tracker retrieves shared PMIDs once, paginates search results, splits queries that exceed PubMed's 9,999-record retrieval window, paces requests, and retries transient NCBI and Sheets read failures. Incomplete API responses stop publication rather than appearing as an empty successful search.

Local outputs include `master.csv`, `history.csv`, `history-before.csv`, `review.csv`, `faculty.csv`, `quarterly.csv`, `summary.csv`, and `analytics.json`. With Sheets input, snapshots are stored in `out/backups/` and `out/sheets-backup.json`.

Publishing sizes grids, replaces report data, removes obsolete trailing rows, and advances the successful collection checkpoint in **one atomic Sheets batch**. It rechecks Config and History immediately before publishing so edits made during collection cause a safe stop. GitHub concurrency prevents overlapping workflow publishers; avoid running another local publisher at the same time.

A network failure during a write can make the result uncertain. The tracker verifies Summary's unique report JSON before declaring success; it does not blindly replay sheet creation. If verification fails, inspect Summary and the saved outputs before retrying.

For rollback, stop publishers, copy the affected workbook, and restore History from History Backup or the saved JSON/CSV. Restore the previous checkpoint from the pre-run snapshot if necessary, then dry-run before publishing again. Repeated runs merge by PMID plus faculty member, so retrying collection does not duplicate history.

Atomic batches above 8 MiB or cells above Google's 50,000-character limit stop with local reports saved. This implementation favors a verifiable atomic report over partial publication; a larger installation may need staged writes or a database.

## GitHub Actions

`Tracker checks` runs on pushes and pull requests, installing the tested core dependency set and running secret-free regression tests on Python 3.11 and 3.13.

`Quarterly PubMed Update` runs at 06:17 UTC on January 1, April 1, July 1, and October 1. It tests before collection and publishes on scheduled runs. Manual dispatch defaults to a **dry run** and supports publication, backfill dates, a report quarter, and optional themes.

Configure the same environment variables above as repository **Actions secrets**. Scheduled and manual publication runs share a concurrency group, have a 30-minute timeout, and leave a run summary. Update the dependency lock intentionally after testing a fresh environment.

This repository is public. Report artifacts can contain review notes and a snapshot of private spreadsheet configuration. Uploads are disabled by default. Set the repository variable `UPLOAD_REPORT_ARTIFACTS=true` only if that artifact visibility is appropriate; enabled artifacts expire after 30 days. Otherwise reports and the previous-history backup remain in Sheets. Logs and run summaries are visible according to repository visibility.

GitHub schedules can be delayed, and [public-repository workflows may be disabled after 60 days of inactivity](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/disable-and-enable-workflows). Check the Actions page for failures or disabled schedules. Use an external scheduler if unattended operation without this GitHub limitation is required.

## Optional Gemini themes

Core tracking does not depend on Gemini or its SDK. To enable themes:

```bash
python -m pip install -r requirements-ai.txt -c requirements.lock
```

Set `themes_enabled=true`, `gemini_model` to a currently available model in your Gemini account, and optionally `theme_max_papers` (blank defaults to 100; maximum 500). Supply `GEMINI_API_KEY`. `--themes` enables it for one run. Actions attempts the optional SDK installation when a Gemini key is configured, and core reporting continues if it fails.

Themes use only accepted publications in the selected quarter. Sampling is deterministic by descending PMID, abstracts are limited to 6,000 characters, and coverage is disclosed. Generated PMID references are checked against the sample; counts are derived from those references. This checks evidence IDs, not the scientific quality of the thematic interpretation. Missing keys, unavailable models, missing SDKs, and provider failures do not prevent reports from being published.

## Email automation

The previous README suggested Apps Script but did not contain an implementation. API writes do not fire edit triggers. See [apps-script/README.md](apps-script/README.md) for recovery steps and a new disabled-by-default clock-trigger integration. It requires deliberate recipient setup and installation; updating this repository alone does not send email or install a Google trigger.

## Project structure

```text
src/script.py          Existing entry point
src/tracker/           Configuration, PubMed, matching, reports, Sheets, CLI, optional AI
examples/Config.tsv    Synthetic, copyable configuration
tests/                 Offline regression tests and PubMed fixture
apps-script/           Optional email integration and setup guide
.github/workflows/     CI and quarterly publication workflow
```

The tracked legacy `config/*.yaml` and journal lookup CSV are not used by the new or previous sheet-driven entry point. They are retained as historical files; the editable roster lives in Config. Files under ignored `data/`, `out/`, `logs/`, and `.env` should remain local.
