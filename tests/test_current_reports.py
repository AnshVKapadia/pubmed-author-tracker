"""Current-window regression coverage, including the reported 2024-paper bug."""
import csv
from datetime import date, datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch

from test_tracker import ROOT, FIXTURE, NOW, configuration, snapshot
from tracker.config import collection_window, parse_config
from tracker.matching import match_article, merge_history
from tracker.pubmed import parse_articles
from tracker.reports import build_reports, date_membership
from tracker.sheets import build_publish_request, report_setting_requests
from tracker.cli import parser, run


class CurrentReportTests(unittest.TestCase):
    def setUp(self):
        _, self.authors = parse_config(configuration())
        articles = parse_articles(FIXTURE.read_text(encoding="utf-8"))
        self.rows = merge_history([], [row for a in articles.values() for person in self.authors
                                      if (row := match_article(a, person))], NOW)

    def test_shared_paper_counts_once_but_credits_both_faculty(self):
        report = build_reports(self.rows, self.authors, "2026-Q3", NOW)
        self.assertEqual({"Master"}, set(report["tables"]))
        self.assertEqual(2, len(report["tables"]["Master"][1]))
        self.assertEqual({"Jane A Smith": 2, "John B Doe": 1, "Ada Zero": 0}, report["analytics"]["faculty_counts"])
        shared = next(r for r in report["tables"]["Master"][1] if r["pmid"] == "1001")
        self.assertIn("Jane A Smith", shared["tracked_author"])
        self.assertIn("John B Doe", shared["tracked_author"])

    def test_old_paper_is_excluded_unless_explicitly_requested(self):
        old = {**self.rows[0], "pmid": "35765233", "pub_date": "2024-06", "status": "AUTHOR MATCH + AFF MATCH"}
        rows = self.rows + [old]
        current = build_reports(rows, self.authors, "2026-Q3", NOW)
        self.assertNotIn("35765233", {r["pmid"] for r in current["tables"]["Master"][1]})
        broad = build_reports(rows, self.authors, "2024-01-01", NOW, date(2026, 9, 30))
        self.assertIn("35765233", {r["pmid"] for r in broad["tables"]["Master"][1]})

    def test_partial_dates_do_not_invent_quarter_membership(self):
        self.assertEqual("uncertain", date_membership({"pub_date": "2026"}, date(2026, 7, 1), date(2026, 9, 30)))
        self.assertEqual("inside", date_membership({"pub_date": "2026"}, date(2026, 1, 1), date(2026, 12, 31)))
        self.assertEqual("uncertain", date_membership({"pub_date": "2026-07"}, date(2026, 7, 15), date(2026, 7, 31)))

    def test_date_window_cannot_be_bypassed_by_an_old_accept_decision(self):
        row = {**self.rows[0], "pub_date": "2026", "decision": "accept"}
        report = build_reports([row], self.authors, "2026-Q3", NOW)
        self.assertEqual([], report["tables"]["Master"][1])
        self.assertEqual(1, len(report["review"]))

    def test_configurable_calendar_periods(self):
        cases = {1: ("2026-09-01", "2026-09-30"), 2: ("2026-07-01", "2026-08-31"),
                 3: ("2026-07-01", "2026-09-30"), 4: ("2026-05-01", "2026-08-31"),
                 6: ("2026-01-01", "2026-06-30"), 12: ("2025-01-01", "2025-12-31")}
        for months, expected in cases.items():
            first, last, _ = collection_window({"report_interval_months": str(months)}, date(2026, 10, 16))
            self.assertEqual(expected, (first.isoformat(), last.isoformat()))
        first, last, _ = collection_window({}, date(2026, 1, 1))
        self.assertEqual(("2025-10-01", "2025-12-31"), (first.isoformat(), last.isoformat()))

    def test_old_checkpoint_does_not_extend_report_dates(self):
        first, last, _ = collection_window({"last_run_utc": "2020-01-01T00:00:00+00:00"}, date(2026, 10, 16))
        self.assertEqual(date(2026, 7, 1), first)

    def test_invalid_interval_is_rejected(self):
        values = configuration()
        values.insert(2, ["report_interval_months", "5"])
        # Remove the example interval to avoid testing duplicate validation instead.
        values = [r for i, r in enumerate(values) if i == 2 or not r or r[0] != "report_interval_months"]
        with self.assertRaises(ValueError):
            parse_config(values)

    def test_report_state_fits_inside_config_even_without_empty_settings_rows(self):
        values = [["[SETTINGS]"], ["last_run_utc", ""], ["[AUTHORS]"], ["full_name", "other_names", "affiliations"], ["Jane Smith", "", "UTSW"]]
        requests = report_setting_requests(values, "{}", 0)
        self.assertEqual(2, len(requests))
        self.assertEqual(2, requests[0]["insertDimension"]["range"]["startIndex"])
        self.assertEqual("last_report_json", requests[1]["updateCells"]["rows"][0]["values"][0]["userEnteredValue"]["stringValue"])

    def test_publish_removes_generated_history_and_summary_tabs_but_not_custom_tabs(self):
        data = snapshot()
        for identifier, title in enumerate(("History", "Quarterly Summary", "Summary", "Notes"), 1):
            data["metadata"]["sheets"].append({"properties": {"sheetId": identifier, "title": title,
                "gridProperties": {"rowCount": 10, "columnCount": 25}}})
        report = build_reports(self.rows, self.authors, "2026-Q3", NOW)
        body = build_publish_request(data, report["tables"], NOW, True, report["analytics"])
        deleted = {r["deleteSheet"]["sheetId"] for r in body["requests"] if "deleteSheet" in r}
        self.assertEqual({1, 2, 3}, deleted)
        self.assertTrue(any("last_report_json" in json.dumps(r) for r in body["requests"]))

    def test_scheduled_run_skips_already_completed_period_without_pubmed(self):
        data = snapshot()
        first, last, _ = collection_window({}, date.today())
        data["tables"]["Config"].insert(2, ["last_report_json", json.dumps({"last_automatic_period": f"{first}/{last}"})])
        store = Mock()
        store.read.return_value = data
        with patch("tracker.cli.SheetsStore", return_value=store), patch("tracker.cli.PubMedClient") as client:
            result = run(parser().parse_args(["--scheduled", "--publish"]))
        self.assertEqual("not_due", result["status"])
        client.assert_not_called()
        store.publish.assert_not_called()

    def test_new_scheduled_period_publishes_and_queues_email_ignoring_manual_dates(self):
        data = snapshot()
        for row in data["tables"]["Config"]:
            if row and row[0] == "starting_date":
                row[:] = ["starting_date", "2020-01-01"]
        store, client = Mock(), Mock()
        store.read.return_value = data
        articles = parse_articles(FIXTURE.read_text(encoding="utf-8"))
        client.search.return_value = set(articles)
        client.fetch.return_value = articles
        with tempfile.TemporaryDirectory(dir=ROOT / "out") as temporary:
            with patch("tracker.cli.SheetsStore", return_value=store), patch("tracker.cli.PubMedClient", return_value=client), patch("tracker.cli.datetime") as clock:
                clock.now.return_value = datetime(2026, 10, 16, tzinfo=timezone.utc)
                result = run(parser().parse_args(["--scheduled", "--publish", "--output", temporary]))
        self.assertEqual("2026-07-01/2026-09-30", result["last_automatic_period"])
        self.assertTrue(result["email_ready"])
        self.assertEqual(result["period_key"], result["email_report"]["period_key"])
        self.assertEqual(2, result["email_report"]["total_articles"])
        store.publish.assert_called_once()

    def test_manual_reports_do_not_queue_email_and_preserve_pending_automatic_email(self):
        with tempfile.TemporaryDirectory(dir=ROOT / "out") as temporary:
            data = snapshot()
            queued = {"period_key": "2026-04-01/2026-06-30", "total_articles": 7}
            data["tables"]["Config"].insert(2, ["last_report_json", json.dumps({"email_report": queued})])
            path = Path(temporary) / "snapshot.json"
            path.write_text(json.dumps(data), encoding="utf-8")
            args = parser().parse_args(["--snapshot", str(path), "--pubmed-xml", str(FIXTURE), "--quarter", "2026-Q3", "--output", temporary])
            result = run(args)
            self.assertFalse(result["email_ready"])
            self.assertEqual(queued, result["email_report"])
            self.assertFalse((Path(temporary) / "history.csv").exists())

    def test_offline_reruns_replace_master_without_adding_history(self):
        with tempfile.TemporaryDirectory(dir=ROOT / "out") as temporary:
            command = [sys.executable, str(ROOT / "src/script.py"), "--config", str(ROOT / "examples/Config.tsv"),
                       "--pubmed-xml", str(FIXTURE), "--quarter", "2026-Q3", "--output", temporary]
            for _ in range(2):
                result = subprocess.run(command, capture_output=True, text=True, encoding="utf-8")
                self.assertEqual(0, result.returncode, result.stderr)
            with (Path(temporary) / "master.csv").open(encoding="utf-8", newline="") as handle:
                self.assertEqual(2, len(list(csv.DictReader(handle))))
            self.assertFalse((Path(temporary) / "history.csv").exists())


if __name__ == "__main__":
    unittest.main()
