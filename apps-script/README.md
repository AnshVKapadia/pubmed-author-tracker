# Optional email summaries

The old repository contains only a README suggestion about Apps Script, not the original source, recipient list, or trigger configuration. Inspection of the shared sheet does not expose bound Apps Script source or the owner's installed triggers. Searches of accessible Drive files and Apps Script notification emails did not locate the original automation. This folder is a **new optional replacement**, not recovered code. No email has been sent and no trigger has been installed.

Google documents that [script executions and API requests do not activate installable triggers](https://developers.google.com/apps-script/guides/triggers/installable#restrictions). Therefore, an `onEdit`/edit trigger cannot detect Python's Sheets API updates. Use a clock trigger that checks the completed Summary report instead.

## Locate an existing script first

1. Open the original spreadsheet as an editor using the account that owns it.
2. Choose **Extensions → Apps Script**. Review the source and the **Triggers** and **Executions** pages.
3. Also check **My Projects** and **My Triggers** at [Apps Script](https://script.google.com/home) in the original owner's account. Bound projects are not ordinary standalone Drive files, and one account cannot see another account's installed triggers.
4. Preserve a copy of the source and note the current recipients and trigger owner before replacing anything.

The connected Sheets/Drive tools cannot read Apps Script project source or list its triggers. To inspect a specific old script here, provide its project link/source or use an authenticated browser with editor access. Do not create competing email triggers blindly.

## Install the replacement only if needed

The versioned `EmailSummary.gs` checks a completed published Summary once per hour and sends at most one report per new quarter. It ignores dry runs, repeated runs, and backfills for already-sent quarters. It does not collect PubMed data.

1. In the spreadsheet's Apps Script editor, add the file's source after inspecting any existing automation.
2. In **Project Settings → Script Properties**, configure:

   | Property | Value |
   | --- | --- |
   | `SPREADSHEET_ID` | The original report spreadsheet ID |
   | `REPORT_RECIPIENTS` | Explicit comma-separated recipient email addresses |
   | `EMAIL_ENABLED` | `true` only when ready to enable sending |
   | `LAST_SENT_QUARTER` | Optional; e.g. `2026-Q3` to suppress sending the current quarter on first installation |

3. Run `installHourlySummaryCheck` as the intended sending account and authorize it. This is a user action; Codex has not performed it.
4. Inspect **Triggers** to verify there is exactly one intended clock trigger. Existing triggers created by another account must be managed from that account.

Emails are sent as the trigger creator. Delivery remains subject to [Apps Script quotas](https://developers.google.com/apps-script/guides/services/quotas). Disable sending by changing `EMAIL_ENABLED` to `false`. There is a small retry ambiguity if delivery succeeds but saving `LAST_SENT_QUARTER` fails; inspect Executions before manually retrying a failed send.

## Summary interface

The tracker writes `Summary!A:B` with `key`/`value` headers, `schema_version=1`, `run_status=complete`, `report_quarter`, `generated_at_utc`, and `analytics_json`. Read by key rather than assuming a fixed analytics-cell position. The data and successful collection checkpoint are published in the same atomic Sheets batch.
