/**
 * Optional replacement email integration. Not installed or enabled by default.
 * Run only after setting the Script Properties documented in README.md.
 */
function checkPublicationSummary() {
  const properties = PropertiesService.getScriptProperties();
  if (properties.getProperty('EMAIL_ENABLED') !== 'true') return;
  const lock = LockService.getScriptLock();
  if (!lock.tryLock(10000)) return;
  try {
    const spreadsheetId = properties.getProperty('SPREADSHEET_ID');
    const recipients = properties.getProperty('REPORT_RECIPIENTS');
    if (!spreadsheetId || !recipients) throw new Error('Configure SPREADSHEET_ID and REPORT_RECIPIENTS first.');
    const spreadsheet = SpreadsheetApp.openById(spreadsheetId);
    const sheet = spreadsheet.getSheetByName('Summary');
    if (!sheet) return;
    const summary = Object.fromEntries(sheet.getDataRange().getValues().slice(1));
    if (summary.run_status !== 'complete' || !summary.analytics_json) return;
    const report = JSON.parse(summary.analytics_json);
    if (report.mode !== 'publish' || !/^\d{4}-Q[1-4]$/.test(report.report_quarter)) return;
    // Send once per new quarter; reruns and backfills do not send duplicate mail.
    const lastSent = properties.getProperty('LAST_SENT_QUARTER') || '';
    if (report.report_quarter <= lastSent) return;
    const message = [
      'PubMed publication report: ' + report.report_quarter,
      'Reporting dates: ' + report.report_start + ' to ' + report.report_end,
      'Unique publications: ' + report.total_articles,
      'Faculty tracked: ' + report.faculty_count,
      'Median publications per faculty: ' + report.per_faculty.median,
      'Pending publication matches for review: ' + report.review_matches,
      'Cumulative unique publications: ' + report.cumulative_total_articles,
      '',
      'Full report: ' + spreadsheet.getUrl()
    ].join('\n');
    MailApp.sendEmail({to: recipients, subject: 'PubMed publication report — ' + report.report_quarter, body: message});
    properties.setProperty('LAST_SENT_QUARTER', report.report_quarter);
  } finally {
    lock.releaseLock();
  }
}

function installHourlySummaryCheck() {
  if (PropertiesService.getScriptProperties().getProperty('EMAIL_ENABLED') !== 'true') {
    throw new Error('Review recipients and set EMAIL_ENABLED=true before installing.');
  }
  const exists = ScriptApp.getProjectTriggers().some(trigger => trigger.getHandlerFunction() === 'checkPublicationSummary');
  if (!exists) ScriptApp.newTrigger('checkPublicationSummary').timeBased().everyHours(1).create();
}
