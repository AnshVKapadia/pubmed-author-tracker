import copy
import csv
from datetime import date
import importlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import Mock, patch
import xml.etree.ElementTree as ET

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tracker.config import Author, collection_window, parse_config, previous_quarter, quarter_window
from tracker.matching import accepted, match_article, merge_history
from tracker.pubmed import PubMedClient, PubMedError, build_query, parse_articles, parse_publication_date
from tracker.reports import analyze, build_reports, quarter_for, summary_table
from tracker.sheets import SheetsStore, build_publish_request, checkpoint_request, history_from_snapshot
from tracker.themes import extract_themes
from tracker.cli import parser, run

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/pubmed.xml"
NOW = "2026-10-16T12:00:00+00:00"


def configuration():
    with (ROOT / "examples/Config.tsv").open(encoding="utf-8", newline="") as handle:
        return list(csv.reader(handle, delimiter="\t"))


def snapshot():
    return {"tables": {"Config": configuration()}, "metadata": {"sheets": [
        {"properties": {"title": "Config", "sheetId": 0, "gridProperties": {"rowCount": 1000, "columnCount": 26}}}]}}


def responses(count, ids):
    response = Mock(status_code=200, headers={})
    response.json.return_value = {"esearchresult": {"count": str(count), "idlist": ids}}
    return response


class TrackerTests(unittest.TestCase):
    def setUp(self):
        self.settings, self.authors = parse_config(configuration())
        self.articles = parse_articles(FIXTURE.read_text(encoding="utf-8"))
        self.rows = [row for article in self.articles.values() for author in self.authors
                     if (row := match_article(article, author))]
        self.history = merge_history([], self.rows, NOW)

    def test_imports_do_not_need_credentials_or_gemini(self):
        for module in ("config", "matching", "pubmed", "reports", "sheets", "themes", "cli"):
            importlib.import_module("tracker." + module)

    def test_shared_publication_credits_each_faculty_once(self):
        report = build_reports(self.history, self.authors, "2026-Q3", NOW)
        analytics = report["analytics"]
        self.assertEqual(2, analytics["total_articles"])
        self.assertEqual({"Jane A Smith": 2, "John B Doe": 1, "Ada Zero": 0}, analytics["faculty_counts"])
        self.assertEqual(3, analytics["cumulative_total_articles"])
        master = report["tables"]["Master"][1]
        self.assertEqual(3, len(master))
        shared = next(r for r in master if r["pmid"] == "1001")
        self.assertIn("Jane A Smith", shared["tracked_author"])
        self.assertIn("John B Doe", shared["tracked_author"])

    def test_other_person_cannot_validate_affiliation(self):
        row = match_article(self.articles["1002"], self.authors[0])
        self.assertIn("AFF MISMATCH", row["status"])
        self.assertFalse(accepted(row))
        self.assertEqual(1, row["other_institutions"])

    def test_aliases_do_not_invent_external_collaborators(self):
        row = match_article(self.articles["1001"], self.authors[0])
        self.assertEqual(0, row["other_institutions"])

    def test_full_author_list_and_nested_xml_text_are_preserved(self):
        article = self.articles["1001"]
        self.assertEqual("A study with nested title text", article["title"])
        self.assertEqual("Before bold evidence and after.", article["abstract"])
        self.assertEqual("Jane Alice Smith; John B Doe", article["authors"])
        self.assertEqual("10.1234/example", article["doi"])

    def test_sole_author_is_first_and_last(self):
        row = match_article(self.articles["1005"], self.authors[1])
        self.assertEqual((1, 1), (row["first_author"], row["last_author"]))

    def test_investigator_has_own_affiliation_and_no_byline_credit(self):
        row = match_article(self.articles["1003"], self.authors[0])
        self.assertTrue(accepted(row))
        self.assertEqual("investigator", row["author_role"])
        self.assertEqual((0, 0), (row["first_author"], row["last_author"]))

    def test_initials_only_and_missing_affiliation_need_review(self):
        row = match_article(self.articles["1004"], self.authors[0])
        self.assertFalse(accepted(row))
        article = copy.deepcopy(self.articles["1001"])
        article["byline"][0]["affiliations"] = []
        row = match_article(article, self.authors[0])
        self.assertIn("AFF NOT FOUND", row["status"])
        self.assertFalse(accepted(row))

    def test_full_given_name_conflict_never_auto_accepts(self):
        article = copy.deepcopy(self.articles["1001"])
        article["byline"][0]["fore"] = "James A"
        row = match_article(article, self.authors[0])
        self.assertTrue(row is None or not accepted(row))

    def test_matching_orcid_and_conflicting_orcid(self):
        author = Author("Jane A Smith", affiliations=("UTSW", "UT Southwestern"), orcid="0000-0001-2345-6789")
        article = copy.deepcopy(self.articles["1001"])
        article["byline"][0]["orcids"] = [author.orcid]
        self.assertTrue(accepted(match_article(article, author)))
        article["byline"][0]["orcids"] = ["0000-0001-2345-6780"]
        self.assertIsNone(match_article(article, author))

    def test_empty_reports_are_finite_and_include_zero_faculty(self):
        report = build_reports([], self.authors, "2026-Q3", NOW)
        self.assertEqual(0, report["analytics"]["total_articles"])
        self.assertEqual(0, report["analytics"]["per_faculty"]["median"])
        self.assertEqual(3, report["analytics"]["faculty_count"])
        json.dumps(report["analytics"], allow_nan=False)
        self.assertEqual([], report["tables"]["Ada Zero"][1])

    def test_zero_faculty_counts_in_percentages(self):
        report = analyze([r for r in self.history if accepted(r)], self.authors)
        self.assertEqual(2, report["faculty_counts"]["Jane A Smith"])
        self.assertAlmostEqual(66.67, report["faculty_productivity_percentages"][">1_publication"])
        self.assertEqual(0, report["per_faculty"]["range"]["min"])

    def test_year_only_dates_not_assigned_a_quarter(self):
        report = build_reports(self.history, self.authors, "2026-Q3", NOW)
        self.assertEqual(1, report["analytics"]["unassigned_date_articles"])
        self.assertEqual("", quarter_for({"pub_date": "2026"}))
        self.assertEqual("2026", parse_publication_date(ET.fromstring("<PubDate><MedlineDate>2026 Winter</MedlineDate></PubDate>")))

    def test_merge_is_idempotent_and_preserves_decisions_and_notes(self):
        before = copy.deepcopy(self.history)
        before[0].update(decision="reject", review_notes="Checked against CV")
        twice = merge_history(before, self.rows, "2026-10-17T12:00:00+00:00")
        self.assertEqual(len(before), len(twice))
        changed = next(r for r in twice if (r["pmid"], r["tracked_author"]) == (before[0]["pmid"], before[0]["tracked_author"]))
        self.assertEqual("reject", changed["decision"])
        self.assertEqual("Checked against CV", changed["review_notes"])
        self.assertEqual(NOW, changed["first_seen_utc"])

    def test_invalid_or_conflicting_review_decisions_stop_run(self):
        row = {**self.history[0], "decision": "maybe"}
        with self.assertRaises(ValueError):
            merge_history([row], [], NOW)
        with self.assertRaises(ValueError):
            merge_history([{**self.history[0], "decision": "accept"}, {**self.history[0], "decision": "reject"}], [], NOW)

    def test_calendar_quarters_and_date_modes(self):
        self.assertEqual("2025-Q4", previous_quarter(date(2026, 1, 1)))
        self.assertEqual((date(2026, 7, 1), date(2026, 9, 30)), quarter_window("2026-Q3"))
        first, last, checkpoint = collection_window({}, date(2026, 10, 16))
        self.assertEqual(date(2026, 7, 1), first)
        self.assertTrue(checkpoint)
        first, last, checkpoint = collection_window({}, date(2026, 10, 16), "2025-01-01", "2025-03-31")
        self.assertFalse(checkpoint)
        self.assertEqual(date(2025, 3, 31), last)

    def test_invalid_dates_and_configuration_fail_early(self):
        for config in ({"ending_date": "2026-01-01"}, {"starting_date": "2026-12-01"}, {"last_run_utc": "2026-12-01T00:00:00+00:00"}):
            with self.assertRaises(ValueError):
                collection_window(config, date(2026, 10, 16))
        for row in (["Master", "", "UTSW"], ["Ada Zero", "", "UTSW"]):
            with self.assertRaises(ValueError):
                parse_config(configuration() + [row])

    def test_checkpoint_insert_stays_inside_settings(self):
        values = [row for row in configuration() if not row or row[0] != "last_run_utc"]
        request = checkpoint_request(values, NOW, 0)
        index = request["updateCells"]["range"]["startRowIndex"]
        self.assertLess(index, values.index(["[AUTHORS]"]))
        values[index] = ["last_run_utc", NOW]
        settings, authors = parse_config(values)
        self.assertEqual(NOW, settings["last_run_utc"])
        self.assertEqual(3, len(authors))

    def test_checkpoint_missing_without_blank_slot_fails_safely(self):
        with self.assertRaises(ValueError):
            checkpoint_request([["[SETTINGS]"], ["[AUTHORS]"]], NOW, 0)

    def test_legacy_migration_recovers_faculty_tabs(self):
        data = snapshot()
        data["tables"].update({"Master": [["pmid", "tracked_author"], ["1001", "Jane A Smith"]],
                               "John B Doe": [["pmid", "tracked_author"], ["1001", "John B Doe"]]})
        migrated = merge_history(history_from_snapshot(data, self.authors), [], NOW)
        self.assertEqual(2, len(migrated))
        self.assertTrue(all(not accepted(r) for r in migrated))

    def test_atomic_batch_contains_checkpoint_and_clears_empty_tabs(self):
        data = snapshot()
        data["tables"]["Ada Zero"] = [["pmid", "tracked_author"], ["1000", "Ada Zero"]]
        data["metadata"]["sheets"].append({"properties": {"title": "Ada Zero", "sheetId": 10,
                                                          "gridProperties": {"rowCount": 2, "columnCount": 2}}})
        report = build_reports(self.history, self.authors, "2026-Q3", NOW)
        report["tables"]["Summary"] = summary_table(report["analytics"])
        body = build_publish_request(data, report["tables"], NOW, True)
        updates = [r["updateCells"] for r in body["requests"] if "updateCells" in r]
        empty = next(r for r in updates if r["range"]["sheetId"] == 10)
        self.assertEqual(2, len(empty["rows"]))
        self.assertTrue(all(not cell for cell in empty["rows"][1]["values"]))
        self.assertEqual(2, empty["range"]["endRowIndex"])
        self.assertEqual(0, updates[-1]["range"]["sheetId"])
        self.assertGreater(len(body["requests"]), 5)

    def test_custom_columns_and_formulas_block_overwrite(self):
        report = build_reports([], self.authors, "2026-Q3", NOW)
        for existing in ([["pmid", "custom"], ["1", "keep me"]], [["pmid"], ["=SUM(A2:A3)"]]):
            data = snapshot()
            data["tables"]["Master"] = existing
            with self.assertRaises(ValueError):
                build_publish_request(data, report["tables"], NOW, False)

    def test_no_clear_calls_and_formula_like_titles_written_as_text(self):
        data = snapshot()
        body = build_publish_request(data, {"Master": (["title"], [{"title": "=IMPORTXML(unsafe)"}])}, NOW, False)
        update = next(r["updateCells"] for r in body["requests"] if "updateCells" in r)
        self.assertEqual({"stringValue": "=IMPORTXML(unsafe)"}, update["rows"][1]["values"][0]["userEnteredValue"])

    def test_manual_run_does_not_include_checkpoint_update(self):
        body = build_publish_request(snapshot(), {"Master": (["pmid"], [])}, NOW, False)
        self.assertTrue(all(r.get("updateCells", {}).get("range", {}).get("sheetId") != 0 for r in body["requests"]))

    def test_search_paginates_beyond_500_and_uses_correct_fields(self):
        session = Mock()
        session.get.side_effect = [responses(1200, [str(i) for i in range(1, 1001)]), responses(1200, [str(i) for i in range(1001, 1201)])]
        client = PubMedClient("test@example.com", session=session, sleep=lambda _: None)
        result = client.search(self.authors[0], date(2026, 7, 1), date(2026, 9, 30))
        self.assertEqual(1200, len(result))
        self.assertEqual(1000, session.get.call_args_list[1].kwargs["params"]["retstart"])
        query = build_query(self.authors[0], date(2026, 7, 1), date(2026, 9, 30))
        self.assertIn("[au]", query)
        self.assertIn("[ir]", query)
        self.assertIn("[crdt]", query)
        self.assertIn('"Smith JA"[au]', query)
        self.assertIn('"UTSW"[ad]', query)

    def test_collection_failure_never_calls_publisher(self):
        (ROOT / "out").mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT / "out") as temporary:
            store = Mock()
            store.read.return_value = snapshot()
            client = Mock()
            client.search.side_effect = PubMedError("simulated collection failure")
            with patch("tracker.cli.SheetsStore", return_value=store), patch("tracker.cli.PubMedClient", return_value=client):
                with self.assertRaises(PubMedError):
                    run(parser().parse_args(["--publish", "--output", temporary]))
            store.publish.assert_not_called()

    def test_publisher_stops_if_review_input_changed(self):
        store = object.__new__(SheetsStore)
        store.sh = Mock()
        store.sh.values_get.return_value = {"values": [["changed during collection"]]}
        report = build_reports([], self.authors, "2026-Q3", NOW)
        report["tables"]["Summary"] = summary_table(report["analytics"])
        with self.assertRaises(RuntimeError):
            store.publish(snapshot(), report["tables"], NOW, True)
        store.sh.batch_update.assert_not_called()

    def test_successful_publish_uses_one_atomic_batch(self):
        data = snapshot()
        store = object.__new__(SheetsStore)
        store.sh = Mock()
        store.sh.values_get.return_value = {"values": data["tables"]["Config"]}
        store._published = Mock(return_value=True)
        report = build_reports([], self.authors, "2026-Q3", NOW)
        report["tables"]["Summary"] = summary_table(report["analytics"])
        store.publish(data, report["tables"], NOW, True)
        store.sh.batch_update.assert_called_once()
        final_request = store.sh.batch_update.call_args.args[0]["requests"][-1]
        self.assertEqual(0, final_request["updateCells"]["range"]["sheetId"])

    def test_reordered_pubmed_links_receive_the_current_pmid(self):
        body = build_publish_request(snapshot(), {"Master": (["pubmed_url"], [{"pubmed_url": "https://pubmed.ncbi.nlm.nih.gov/1001/"}])}, NOW, False)
        update = next(r["updateCells"] for r in body["requests"] if "updateCells" in r)
        self.assertIn("userEnteredFormat.textFormat.link", update["fields"])
        link = update["rows"][1]["values"][0]["userEnteredFormat"]["textFormat"]["link"]["uri"]
        self.assertEqual("https://pubmed.ncbi.nlm.nih.gov/1001/", link)

    def test_api_errors_and_incomplete_searches_never_become_empty_success(self):
        for response in (Mock(status_code=400), responses(3, []), Mock(status_code=200, json=lambda: {"error": "failure"})):
            client = PubMedClient("test@example.com", session=Mock(get=Mock(return_value=response)), sleep=lambda _: None)
            with self.assertRaises(PubMedError):
                client.search(self.authors[0], date(2026, 7, 1), date(2026, 9, 30))

    def test_retry_transient_errors_without_exposing_api_key(self):
        session = Mock()
        session.get.side_effect = [Mock(status_code=429, headers={"Retry-After": "2"}), responses(0, [])]
        sleeps = []
        client = PubMedClient("test@example.com", api_key="sensitive-key", session=session, sleep=sleeps.append)
        self.assertEqual(set(), client.search(self.authors[0], date(2026, 7, 1), date(2026, 9, 30)))
        self.assertIn(2, sleeps)
        session.get.side_effect = None
        session.get.return_value = Mock(status_code=403)
        with self.assertRaises(PubMedError) as caught:
            client.request("esearch", {})
        self.assertNotIn("sensitive-key", str(caught.exception))

    def test_missing_fetch_and_invalid_xml_stop_collection(self):
        response = Mock(status_code=200, headers={}, text="<PubmedArticleSet/>")
        client = PubMedClient("test@example.com", session=Mock(post=Mock(return_value=response)))
        with self.assertRaises(PubMedError):
            client.fetch({"1001"})
        with self.assertRaises(PubMedError):
            parse_articles("bad xml")

    def test_themes_disabled_missing_key_and_provider_failure_are_nonfatal(self):
        self.assertEqual("disabled", extract_themes([], {}, "")["status"])
        self.assertEqual("no_articles", extract_themes([], {"themes_enabled": "true"}, "")["status"])
        self.assertEqual("skipped", extract_themes(list(self.articles.values()), {"themes_enabled": "true"}, "")["status"])
        with patch.dict(sys.modules, {"google.genai": None}):
            result = extract_themes(list(self.articles.values()), {"themes_enabled": "true", "gemini_model": "example"}, "dummy")
        self.assertEqual("unavailable", result["status"])

    def test_cli_offline_end_to_end_and_rerun(self):
        (ROOT / "out").mkdir(exist_ok=True)
        with tempfile.TemporaryDirectory(dir=ROOT / "out") as temporary:
            args = [sys.executable, str(ROOT / "src/script.py"), "--config", str(ROOT / "examples/Config.tsv"),
                    "--pubmed-xml", str(FIXTURE), "--quarter", "2026-Q3", "--output", temporary]
            first = subprocess.run(args, capture_output=True, text=True, encoding="utf-8")
            self.assertEqual(0, first.returncode, first.stderr)
            output = Path(temporary)
            analytics = json.loads((output / "analytics.json").read_text(encoding="utf-8"))
            self.assertEqual(2, analytics["total_articles"])
            self.assertEqual("dry_run", analytics["mode"])
            with (output / "history.csv").open(encoding="utf-8", newline="") as handle:
                original = list(csv.DictReader(handle))
            second = subprocess.run(args + ["--history", str(output / "history.csv")], capture_output=True, text=True, encoding="utf-8")
            self.assertEqual(0, second.returncode, second.stderr)
            with (output / "history.csv").open(encoding="utf-8", newline="") as handle:
                self.assertEqual(len(original), len(list(csv.DictReader(handle))))

    def test_cli_refuses_publishing_fixtures(self):
        result = subprocess.run([sys.executable, str(ROOT / "src/script.py"), "--config", str(ROOT / "examples/Config.tsv"),
                                 "--pubmed-xml", str(FIXTURE), "--publish"], capture_output=True)
        self.assertEqual(1, result.returncode)


if __name__ == "__main__":
    unittest.main()
