// Runs only mocked Google services; never sends real email or installs a trigger.
const assert = require('node:assert/strict');
const fs = require('node:fs');
const vm = require('node:vm');
const properties = {};
const sent = [];
let rows = [['[SETTINGS]'], ['[AUTHORS]']];
const context = {
  PropertiesService: { getScriptProperties: () => ({ getProperty: key => properties[key], setProperty: (key, value) => { properties[key] = value; } }) },
  LockService: { getScriptLock: () => ({ tryLock: () => true, releaseLock: () => {} }) },
  SpreadsheetApp: { openById: () => ({ getSheetByName: name => name === 'Config' ? { getDataRange: () => ({ getValues: () => rows }) } : null,
                                    getUrl: () => 'https://example.invalid/test-workbook' }) },
  MailApp: { sendEmail: message => sent.push(message) }
};
vm.createContext(context);
vm.runInContext(fs.readFileSync('apps-script/EmailSummary.gs', 'utf8'), context);
context.checkPublicationSummary();
assert.equal(sent.length, 0, 'disabled integration must do nothing');
Object.assign(properties, { EMAIL_ENABLED: 'true', SPREADSHEET_ID: 'test-only', REPORT_RECIPIENTS: 'nobody@example.invalid' });
function report(start, end) {
  return { period_key: `${start}/${end}`, report_start: start, report_end: end, total_articles: 2, faculty_count: 3,
           per_faculty: { median: 0 }, review_matches: 1 };
}
function queue(value) { rows = [['[SETTINGS]'], ['last_report_json', JSON.stringify(value)], ['[AUTHORS]']]; }
queue({ mode: 'dry_run', email_ready: false });
context.checkPublicationSummary();
assert.equal(sent.length, 0, 'manual previews must not queue email');
queue({ email_report: report('2026-07-01', '2026-09-30') });
context.checkPublicationSummary();
context.checkPublicationSummary();
assert.equal(sent.length, 1, 'repeated clock checks must send once');
assert.match(sent[0].subject, /2026-07-01 to 2026-09-30/);
queue({ email_report: report('2026-09-01', '2026-09-30') });
context.checkPublicationSummary();
assert.equal(sent.length, 1, 'changing interval must not resend a report ending on the same date');
queue({ email_report: report('2026-10-01', '2026-10-31') });
context.checkPublicationSummary();
assert.equal(sent.length, 2, 'a new monthly period must send');
console.log('Email automation checks passed (mocked services; no email sent).');
