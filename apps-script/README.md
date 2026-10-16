# Autonomous publication-summary emails

The reporting period belongs to the email automation, not a history worksheet. The tracker keeps one current Master. Scheduled runs use Config's `report_interval_months` (default `3`) to rebuild it for the latest completed calendar period.

## How it runs

1. GitHub checks daily. If a new period is due, it collects publications and atomically publishes Master plus `Config.last_report_json`.
2. The metadata includes an email summary for that automatic report. Manual runs never queue a new email and preserve an existing pending summary.
3. An Apps Script clock trigger checks once per hour and sends the queued summary once. It remembers the most recently emailed period and does not resend reports ending on the same or an earlier date, including after interval changes.

Google documents that [API writes do not activate edit triggers](https://developers.google.com/apps-script/guides/triggers/installable#restrictions), so this must be a clock trigger, not `onEdit`.

## Setup (not enabled automatically)

First inspect **Extensions → Apps Script** in the original spreadsheet using its owner's account. Check Triggers and Executions before adding a replacement. The original source and recipients were not found in the repository or connected searches. This file is a new replacement, not recovered code.

In the Apps Script editor, add `EmailSummary.gs`. In **Project Settings → Script Properties**, configure:

| Property | Value |
| --- | --- |
| `SPREADSHEET_ID` | The original tracker spreadsheet ID |
| `REPORT_RECIPIENTS` | Explicit comma-separated recipient email addresses |
| `EMAIL_ENABLED` | `true` only when ready to send |
| `LAST_SENT_PERIOD` | Optional; e.g. `2026-07-01/2026-09-30` to suppress the current report on first setup |

Run `installHourlySummaryCheck` as the intended sender and authorize it. Verify there is one intended clock trigger. The assistant has not installed a trigger, configured recipients, or sent email. Disable sending with `EMAIL_ENABLED=false`.

Emails are sent as the trigger creator and are subject to [Apps Script quotas](https://developers.google.com/apps-script/guides/services/quotas). Another account's triggers must be managed from that account. If sending succeeds but saving the sent-period marker fails, a later retry may duplicate delivery; inspect Executions before manually retrying.

To change cadence later, edit `report_interval_months` inside Config's `[SETTINGS]`: `1`, `2`, `3`, `4`, `6`, or `12`. No cron or Apps Script edits are needed. GitHub processes the latest completed period; it does not email every historical missed period after a long outage.

## Testing

`node tests/test_email.js` uses mocked Google services to test disabled mode, manual-run behavior, duplicate prevention, interval changes, and new-period delivery. It never sends email or creates a real trigger.
