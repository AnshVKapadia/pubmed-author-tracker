# Publication tracker: audit and implementation plan

Audited October 16, 2026. The baseline findings below informed the implementation. The owner selected Google Sheets, retained history, and quarterly summaries.

## Implementation status

The tracker has been split into testable modules; core dependencies are pinned with a tested full lock; Google Sheets publication is atomic; history retains every faculty-paper match and manual decision; AI is optional; and CI plus the quarterly workflow have been rewritten. See README.md for current behavior.

Validation completed: 35 regression tests in a clean Python 3.13 environment, a successful live PubMed dry run against the actual Config roster, and a full atomic publish to a separate copy of the workbook, followed by readback of summary metrics, history count, faculty rows, and PubMed hyperlinks. A known shared paper now credits all five configured faculty involved. The production workbook and checkpoint have not been changed.

The original Apps Script source and triggers were not located in the repository, accessible Drive files, or the connected inbox's Apps Script notifications. apps-script/ contains an explicitly labeled optional replacement and recovery instructions; no recipients or triggers have been enabled. GitHub-hosted workflow verification and production service-account validation remain separate from the connector-based test-copy publication.

## Recommendation

Keep Python, PubMed, Google Sheets, and GitHub Actions for the first repair. Sheets already provides an editable author roster and a useful review interface. Replacing it with a web application would add deployment, authentication, and storage work without fixing the publication-matching problems.

Separate collection and analysis from Sheets publishing. Produce CSV and JSON reports locally before publishing so results can be checked, downloaded, and recovered. Make Gemini an optional final enrichment step whose failure cannot prevent publication tracking.

Proposed default, pending the owner's reporting preference: preserve publication history and produce calendar-quarter summaries. Use one record per publication plus separate publication-to-faculty matches, so a paper counts once department-wide while crediting every relevant faculty member.

## Evidence and current failures

### Workflow startup

- GitHub reports four scheduled runs in 2026: January 1, April 1, July 1, and October 1. All failed. The two earlier manual runs in December 2025 succeeded on older commits.
- The [October 1 run](https://github.com/AnshVKapadia/pubmed-author-tracker/actions/runs/36860472359) installed dependencies successfully, then failed at `src/script.py:351` with `SyntaxError: f-string: expecting '}'`.
- `.github/workflows/quarterly.yml` selects Python 3.11. Line 351 reuses double quotes inside a double-quoted f-string, a feature introduced in Python 3.12. Fix the expression and test on the declared Python versions.
- `requirements.txt:9` installs `google-generativeai`, while `src/script.py:23` imports `google.genai`, which belongs to `google-genai`. This is a separate startup blocker once the syntax is repaired.
- The latest run shows masked, populated values for all four configured secrets. Their presence is established; their validity and Google permissions are not.

### Publication accuracy

| Problem | Evidence | Planned correction |
| --- | --- | --- |
| Shared papers lose faculty credit | `script.py:691` drops duplicate PMIDs before faculty grouping and author tabs. A two-faculty example retains only the first faculty member. | Deduplicate publications by PMID and matches by PMID plus faculty identity. Calculate each statistic at the appropriate level. |
| Another author's affiliation can validate the tracked person | Lines 316-328 compare institution keywords against every affiliation in the article. Reproduced with a tracked person at another university and a coauthor at UTSW. | Match a person's identity and that person's own affiliation together. Keep ambiguous or missing evidence in the review queue. |
| Collaboration counts are incorrect | Line 329 checks all keyword/affiliation combinations. A single UTSW affiliation with two aliases can be labeled external collaboration. | An affiliation is internal if any configured institutional alias matches; classify genuinely external affiliations separately. |
| Author lists are truncated | Lines 339-359 stop building the displayed author list at the first matched author. | Parse the complete author list before classifying matches and positions. |
| Matching depends on logging | The identity and affiliation filtering at line 395 runs only when a debug stream is supplied. | Make matching independent of logging; logging only records decisions. |
| Review candidates can disappear | Line 419 skips regular author matches with affiliation mismatch, instead of retaining them in the documented review tab. | Retain uncertain matches with explicit reasons. |
| Zero-publication faculty are excluded | Lines 698-709 use only faculty present in accepted results. Medians and productivity percentages use the wrong population. | Include every configured faculty member, including zero counts. |
| Empty results can fail after publishing | An isolated empty-result run reaches the write and checkpoint calls, then raises `ValueError: cannot convert float NaN to integer` at the analytics construction. All external calls and writes were stubbed during this check. | Define empty statistics, validate the report before publishing, and advance the checkpoint only after required outputs succeed. |
| Authorship counts do not use the accepted dataset | Lines 363-367 count first/last authors before confidence filtering. A single author is counted as first but never last because of `elif`. | Calculate documented authorship metrics from accepted publication-faculty matches. |
| Some title/abstract text is lost | `_text` and `extract_abstract` read `.text` rather than all nested XML text. A synthetic title containing italics loses the italicized text and the remainder. | Preserve XML text with `itertext()` and handle alternate publication-date formats. |

### Dates, history, and publishing

- The date window starts 91 days before the previous run, not at a calendar-quarter boundary. A run three months after its predecessor can search roughly six months, and the README describes different behavior. Define report periods independently of collection catch-up and overlap.
- Publication-date-only searches can miss older papers added to PubMed later. Establish a bounded recheck policy and explicitly document its limits; validate PubMed indexing/update-date options before choosing a collection strategy.
- Output tabs are cleared and replaced; there is no cumulative history merge. A later run can remove older records or review work.
- Google Sheets writes are separate clear/update operations. A failure after clearing can leave an empty or partly updated report. There is pacing but no Sheets retry/backoff strategy.
- New output worksheets are created with one row and one column, then larger tables are sent without resizing. Allocate the required grid size before writing.
- Faculty tabs with no current results are skipped, leaving old results visible without indicating they are stale.
- If `last_run_utc` is missing, `update_config_value` appends it after the authors block. A parser round trip treats it as a faculty member. Insert settings inside the settings section and reject malformed configuration.
- Gemini is initialized at import time and fails without a key, despite being documented as optional. A failed Gemini call also prevents all subsequent report publishing.
- The theme input is an unordered half of the accepted papers, not a fixed safety limit, and the schema demands exactly five themes even for very small datasets. Use deterministic inputs, a real limit with disclosed coverage, and a schema suitable for the dataset.
- Analytics are printed to stdout but never written to a summary cell or tab. The README's proposed Apps Script email trigger therefore has no implemented analytics-cell update to react to. Any existing Apps Script must be inspected before changing its interface.

### Maintenance and automation

- There are no automated tests, pull-request checks, clean-install checks, or offline dry-run controls.
- Dependencies have minimum versions only; there is no reproducible lock or tested constraints file.
- The scheduled workflow has no overlap protection, explicit job timeout, result artifacts, or run summary.
- The latest job log warns that the checkout and setup-python action versions use a deprecated Node runtime. Update to supported releases as part of workflow maintenance.
- `config/` is ignored but its legacy YAML files and roughly 11 MB journal lookup CSV remain tracked. They are not read by the current pipeline. Inventory them before removing obsolete files; ignoring files does not untrack them.
- The local GitHub CLI reports invalid stored tokens. The connected GitHub tool did allow read-only inspection. Only repair CLI authentication if a later step actually needs it.
- GitHub documents that public-repository schedules can be disabled after 60 days without repository activity. This is an operational limitation to address/document, not the demonstrated cause of the October failure. A quarterly cron is not an exact-time delivery guarantee.

## Implementation sequence

1. **Establish a testable baseline.** Add focused regression tests around the reproduced defects, make imports safe without credentials, fix Python compatibility and the SDK dependency, and introduce a reproducible dependency installation. Keep the existing script entry point.
2. **Correct collection and matching.** Separate configuration, PubMed access/XML parsing, identity matching, analytics, and publishing into small modules. Use author-qualified searches, inspect search counts, retrieve results beyond the current 500-record cap, throttle NCBI calls, retry transient failures, and reject API error payloads rather than treating them as empty results. Fetch shared PMIDs once per run.
3. **Define reporting and persistence.** Apply the agreed report periods, preserve every faculty-paper relationship, retain review decisions, and make repeated runs idempotent. Validate dates, duplicate faculty entries, reserved/invalid worksheet names, and configuration structure before collection. Decide the history/backfill range with the owner.
4. **Publish safely.** Generate validated CSV/JSON reports first. Read existing state where needed, size and batch Sheets writes, clear obsolete trailing data in the same publish operation where practical, and preserve human-maintained review fields. Record the run status and last successful collection only after required publishing succeeds. Expose a Summary tab with an agreed interface for existing Apps Script.
5. **Repair automation.** Add secret-free CI on pushes and pull requests; use supported action versions, dependency caching, minimum permissions, concurrency control, a timeout, downloadable reports/logs, and a useful GitHub run summary. Keep manual execution and quarterly scheduling. Add manual date overrides and a dry run that cannot update Sheets or checkpoints. Surface failures through GitHub's existing workflow status; add email behavior only after reviewing the current setup.
6. **Validate the integration.** Run offline tests and a clean dependency install, then a read-only configuration check and a small live PubMed dry run. Compare known publications against expected faculty credit. Publish to a test spreadsheet/copy, check a repeated run, an empty period, a missing Gemini key, and simulated API failures. Finally validate the workflow on GitHub and document setup, recovery, backfills, and operational limits.

## Information needed from the owner

Start with:

1. **The existing spreadsheet link.** Inspect its Config layout, tabs, existing history, formulas, and any manually edited fields before designing migration or publishing.
2. **The intended report.** Latest quarter, cumulative history with quarterly summaries, or another period; desired historical start date; whether to retain individual faculty tabs and AI themes.
3. **Any attached Apps Script and triggers.** Especially email summaries, recipient selection, summary-cell expectations, and scheduled triggers. The script source is outside this repository.

For final integration validation, confirm or inspect through an authorized route:

- The Google Cloud project/service-account email, whether Sheets API is enabled, and whether the target spreadsheet is shared with that service account as an editor. Do not assume a new project or key is necessary.
- Presence and validity of `SPREADSHEET_ID` and `GOOGLE_APPLICATION_CREDENTIALS_JSON` in GitHub Actions. `NCBI_API_KEY` is optional with a lower request rate; `GEMINI_API_KEY` should only be needed if themes are enabled.
- Whether an institutional match must be the faculty member's own affiliation; how consortium/investigator credit, former affiliations, and manual review decisions should be handled. Default to the person's own institution and keep uncertainty reviewable.
- A few publications that should be included and, if available, known false positives for validation.

Do not paste private keys or API-key values into chat. Update credentials through GitHub Secrets or an ignored local environment file if validation shows that an update is necessary.

## Alternatives

- **Keep Sheets (recommended):** retain the editable roster, review queue, and shared reports while fixing the pipeline.
- **CSV/JSON reports only:** remove Google service-account authentication and sheet publishing entirely. Appropriate if editable/shared Sheets and its email integration are not needed. An automatic runner still needs somewhere durable for cumulative state; expiring Actions artifacts alone are not a permanent database.
- **A dashboard/database:** consider later if review history, permissions, or reporting volume require it. It introduces more infrastructure than this repair needs.
- **Move everything into Apps Script:** potentially reduces platforms, but requires checking execution quotas and workload, rewriting the collector, and providing a testing/recovery strategy. It is not yet established as a simpler solution for this roster.

## Completion criteria

- A clean checkout installs and passes CI without production secrets.
- The documented local and Actions entry points work.
- One shared paper counts once in department totals and once for every matching faculty member.
- Missing affiliation evidence is reviewable, and another person's institution cannot validate the tracked author.
- Empty periods and disabled/unavailable AI do not break the tracker.
- Failed collection does not erase reports or advance the successful-run checkpoint.
- Repeated runs preserve history and human review decisions without adding duplicate matches.
- Outputs disclose reporting dates, collection time, counts, warnings, and AI coverage where relevant.
- A test-sheet integration run and a GitHub Actions run succeed before claiming the production integration is verified.

## References checked

- [Python 3.12 f-string changes](https://docs.python.org/3/whatsnew/3.12.html#pep-701-syntactic-formalization-of-f-strings)
- [Google GenAI SDK migration](https://ai.google.dev/gemini-api/docs/migrate)
- [GitHub scheduled-workflow inactivity behavior](https://docs.github.com/en/actions/how-tos/manage-workflow-runs/disable-and-enable-workflows)
- [gspread worksheet API](https://docs.gspread.org/en/latest/api/models/worksheet.html)

The audit used source inspection, GitHub run metadata/logs, installed-library inspection, and isolated in-memory reproductions with synthetic publication records. It did not authenticate to the Google Sheet, contact Gemini, execute the production tracker, or modify a live report.
