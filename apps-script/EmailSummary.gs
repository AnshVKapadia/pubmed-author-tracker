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
    const sheet = spreadsheet.getSheetByName('Config');
    if (!sheet) return;
    let metadata = null;
    let inSettings = false;
    for (const row of sheet.getDataRange().getValues()) {
      if (row[0] === '[SETTINGS]') inSettings = true;
      if (row[0] === '[AUTHORS]') break;
      if (inSettings && row[0] === 'last_report_json' && row[1]) metadata = JSON.parse(row[1]);
    }
    const report = metadata && metadata.email_report;
    if (!report || !/^\d{4}-\d{2}-\d{2}\/\d{4}-\d{2}-\d{2}$/.test(report.period_key)) return;
    const lastSent = properties.getProperty('LAST_SENT_PERIOD') || '';
    const lastEnd = lastSent.includes('/') ? lastSent.split('/')[1] : '';
    if (report.report_end <= lastEnd) return;
    const message = [
      'PubMed publication report: ' + report.report_start + ' to ' + report.report_end,
      'Reporting dates: ' + report.report_start + ' to ' + report.report_end,
      'Unique publications: ' + report.total_articles,
      'Faculty tracked: ' + report.faculty_count,
      'Median publications per faculty: ' + report.per_faculty.median,
      'Pending publication matches for review: ' + report.review_matches,
      '',
      'Tracker workbook (Master shows the latest run): ' + spreadsheet.getUrl()
    ].join('\n');
    MailApp.sendEmail({to: recipients, subject: 'PubMed publication report — ' + report.report_start + ' to ' + report.report_end, body: message});
    properties.setProperty('LAST_SENT_PERIOD', report.period_key);
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
