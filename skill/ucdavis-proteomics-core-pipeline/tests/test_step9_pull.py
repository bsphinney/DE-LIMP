#!/usr/bin/env python3
"""The step-9 pull brings what a report built on a laptop reads (2.15 review2 HIGH, H4).

Core staff on a laptop build the report from a pull of the HIVE session (`core_sync.sh --get`,
SKILL.md step 1c.7 and core-submissions.md; the slim-down moved the router's hand-copied loop
there, so its GET_RECORDS is the one pull list). The notes read input/raw_names.tsv when no stage
record is above the copy (submission_report.stage_map), and a Core client work mark counts only
with logs/core_client.staff.json (core_client_mark). The list had neither: a laptop-built report
then said stage_record_not_found and matched the runs by their names as acquired, or left the
survey out, while AGENTS.md (finalize on HIVE) was right.

This test reads the record list from core_sync.sh itself, copies exactly those files from a
HIVE-like session (a stage record above it) into a laptop folder (none above), and compares the
notes and core_run. It fails loudly when the list cannot be found. Invented records only
(PROT_9901 and runs of the fixtures). Stdlib only, no network.
"""
import json
import os
import re
import shutil
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL_DIR = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL_DIR, "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import staff  # noqa: E402
import submission_report as sr  # noqa: E402
from test_sheet_names_as_acquired import (EXTRA_ACQ, EXTRA_DEL, N, acquired_name,  # noqa: E402
                                          plate_session, run)


def pull_list():
    """The session records `core_sync.sh --get` pulls: its GET_RECORDS, the one list (the
    output/ items and the search's files are pulled too, but no record the notes read is
    among them)."""
    with open(os.path.join(SCRIPTS, "core_sync.sh"), encoding="utf-8") as fh:
        text = fh.read()
    found = re.findall(r'^GET_RECORDS="([^"]*)"', text, re.M)
    if len(found) != 1:
        raise AssertionError(f"core_sync.sh must assign GET_RECORDS exactly once, found {len(found)}")
    return found[0].split()


def pull(hive_session, laptop):
    """Copy the records core_sync.sh --get pulls, and nothing else; a file not there is
    skipped (core_sync.sh reports it "missing", harmlessly)."""
    for rel in pull_list():
        src = os.path.join(hive_session, rel)
        if os.path.isfile(src):
            os.makedirs(os.path.dirname(os.path.join(laptop, rel)), exist_ok=True)
            shutil.copy2(src, os.path.join(laptop, rel))
    os.makedirs(os.path.join(laptop, "output"), exist_ok=True)


def notes(s):
    return [(n["id"], n["text"]) for n in sr.quality_notes(sr.load(s), s)]


class TheStep9Pull(unittest.TestCase):
    def test_it_lists_every_record_the_laptop_report_reads(self):
        got = pull_list()
        for rel in ("session.json", "input/submission.json", "input/raw_files.txt",
                    "input/raw_names.tsv", sr.CORE_CLIENT_STAFF_RECORD.replace(os.sep, "/")):
            self.assertIn(rel, got)
        ref = os.path.join(SKILL_DIR, "references", "core-submissions.md")
        with open(ref, encoding="utf-8") as fh:
            part = " ".join(fh.read().split())
        self.assertIn("as `core_sync.sh --get` lists them (the one list)", part)
        # and no second copy of the list: the router names core_sync.sh, never a loop of records
        with open(os.path.join(SKILL_DIR, "SKILL.md"), encoding="utf-8") as fh:
            router = fh.read()
        self.assertIn("core_sync.sh --get", router)
        self.assertNotIn("locate_decisions.json raw_names.tsv", router)
        self.assertNotRegex(router, r'for f in [^;]*submission\.json[^;]*; do')

    def test_a_renamed_plate_pulled_to_a_laptop_says_what_hive_says(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = [acquired_name(i) for i in range(1, N + 1)] + [EXTRA_ACQ]
            ren = {acquired_name(i): run(i) for i in range(1, 11)}
            ren[EXTRA_ACQ] = EXTRA_DEL
            s = plate_session(tmp, raws, renames=ren)
            self.assertEqual(sr.write_name_map(s)[0], "OK")       # as attach does, on HIVE
            laptop = os.path.join(tmp, "laptop", "core", "PROT_9901", "session")
            pull(s, laptop)
            hive, here = notes(s), notes(laptop)
            self.assertIn("raw_names_corrected", dict(here))
            self.assertNotIn("stage_record_not_found", dict(here))
            self.assertEqual(here, hive)

    def test_a_marked_session_pulled_to_a_laptop_carries_the_survey(self):
        import session as session_mod
        from test_core_client_survey import STAFF, mark, with_analysis
        with tempfile.TemporaryDirectory() as tmp:
            s = with_analysis(tmp)
            self.assertEqual(mark(s, STAFF)[0], 0)
            laptop = os.path.join(tmp, "laptop", "session")
            pull(s, laptop)
            self.assertEqual(sr.core_run(laptop), {"prot": None, "core_client": True})
            with mock.patch.object(staff, "core_staff", lambda: STAFF):
                self.assertIsNone(sr.survey_off_note(laptop))
            def sj(x):
                with open(session_mod.paths_for(x)["session_json"], encoding="utf-8") as fh:
                    return json.load(fh)
            self.assertEqual(sj(laptop), sj(s))


if __name__ == "__main__":
    unittest.main()
