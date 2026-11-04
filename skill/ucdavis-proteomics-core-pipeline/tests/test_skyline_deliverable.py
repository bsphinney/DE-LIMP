"""
skyline_peaks.deliverable(output_dir) -> (ok, state, why) is the ONE rule for whether the peak
check ships (architectural rule 3; C1 F3). core_submission.skyline_deliverable is a thin wrapper
returning (ok, why) from it, so the two can never disagree.

The rule does not require the Skyline document: a finished ("ok") record is deliverable on its
own, as core_submission.skyline_deliverable has always decided.

Stdlib only, no network.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))

import core_submission as cs  # noqa: E402
import skyline_peaks as sp  # noqa: E402


def write_record(out, rec=None, raw=None):
    sky = os.path.join(out, "skyline")
    os.makedirs(sky, exist_ok=True)
    with open(os.path.join(sky, sp.RECORD), "w") as f:
        f.write(raw if raw is not None else json.dumps(rec))
    return sky


class Deliverable(unittest.TestCase):
    def setUp(self):
        self._t = tempfile.TemporaryDirectory()
        self.out = self._t.name
        self.addCleanup(self._t.cleanup)

    def test_missing(self):
        ok, state, why = sp.deliverable(self.out)
        self.assertFalse(ok)
        self.assertEqual(state, "missing")
        self.assertIn(sp.RECORD, why)

    def test_unreadable(self):
        write_record(self.out, raw="not json {")
        ok, state, why = sp.deliverable(self.out)
        self.assertEqual((ok, state), (False, "unreadable"))
        self.assertEqual(why, sp.STATE_WORDS["unreadable"])

    def test_skipped_carries_reason(self):
        write_record(self.out, {"schema": sp.SCHEMA, "status": "skipped", "reason": "no Cys fixed mod"})
        ok, state, why = sp.deliverable(self.out)
        self.assertEqual((ok, state), (False, "skipped"))
        self.assertIn("no Cys fixed mod", why)

    def test_skipped_without_reason(self):
        write_record(self.out, {"schema": sp.SCHEMA, "status": "skipped"})
        ok, state, why = sp.deliverable(self.out)
        self.assertFalse(ok)
        self.assertIn("no reason recorded", why)

    def test_planned_and_submitted(self):
        sky = write_record(self.out, {"schema": sp.SCHEMA, "status": "planned"})
        ok, state, why = sp.deliverable(self.out)
        self.assertEqual((ok, state), (False, "planned"))
        self.assertIn(sp.STATE_WORDS["planned"], why)
        os.makedirs(os.path.join(sky, sp.WORK))
        open(os.path.join(sky, sp.WORK, "jobs.txt"), "w").close()
        ok, state, why = sp.deliverable(self.out)
        self.assertEqual((ok, state), (False, "submitted"))
        self.assertIn(sp.STATE_WORDS["submitted"], why)

    def test_stale(self):
        write_record(self.out, {"schema": sp.SCHEMA, "status": "ok",
                                "de_record_sha256": "0" * 64})
        ok, state, why = sp.deliverable(self.out)
        self.assertEqual((ok, state), (False, "stale"))
        self.assertEqual(why, sp.STATE_WORDS["stale"])

    def test_ok_is_deliverable_without_document(self):
        write_record(self.out, {"schema": sp.SCHEMA, "status": "ok"})
        ok, state, why = sp.deliverable(self.out)
        self.assertEqual((ok, state, why), (True, "ok", None))

    def test_ok_with_document(self):
        sky = write_record(self.out, {"schema": sp.SCHEMA, "status": "ok"})
        open(os.path.join(sky, "T" + sp.SKYLINE_DOCUMENT_SUFFIX), "w").close()
        self.assertEqual(sp.deliverable(self.out), (True, "ok", None))

# The reviewer's 16 output folders (L9 review F10). FROZEN is the old rule's (ok, why) on each:
# core_submission.skyline_deliverable as it was at 71f1e143, before b8a04605 made it a wrapper of
# skyline_peaks.deliverable(), computed ONCE from that commit's scripts and written here as
# literals. Comparing the wrapper with the function it calls (the test this replaces) could not
# catch a drift from the old rule; this table does.
FROZEN = [
    ("no folder", (False, "output/skyline/ holds no peak-check record (skyline_provenance.json)")),
    ("empty folder", (False, "output/skyline/ holds no peak-check record (skyline_provenance.json)")),
    ("garbage", (False, "output/skyline/skyline_provenance.json could not be read as a peak-check "
                        "record")),
    ("a list", (False, "output/skyline/skyline_provenance.json could not be read as a peak-check "
                       "record")),
    ("null", (False, "output/skyline/skyline_provenance.json could not be read as a peak-check "
                     "record")),
    ("skipped with a reason", (False, "the peak check was skipped: the container is gone")),
    ("skipped without a reason", (False, "the peak check was skipped: no reason recorded")),
    ("skipped with an empty reason", (False, "the peak check was skipped: no reason recorded")),
    ("planned", (False, "the peak check is not finished: it was planned, but its jobs were never "
                        "submitted")),
    ("submitted", (False, "the peak check is not finished: its jobs were submitted and had not "
                          "written their results; make the report again once they finish")),
    ("ok", (True, None)),
    ("ok and stale", (False, "the differential expression was run again after this peak check, so "
                             "the proteins it shows were chosen from an earlier result; it needs "
                             "to be run again")),
    ("ok and current", (True, None)),
    ("ok made with no DE, a DE now", (False, "the differential expression was run again after this "
                                             "peak check, so the proteins it shows were chosen "
                                             "from an earlier result; it needs to be run again")),
    ("an unknown status", (False, "output/skyline/skyline_provenance.json could not be read as a "
                                  "peak-check record")),
    ("no status", (False, "output/skyline/skyline_provenance.json could not be read as a "
                          "peak-check record")),
]


def build_case(o, name):
    """Output folder `o` as the case `name` of FROZEN has it."""
    sky = os.path.join(o, "skyline")

    def de():
        os.makedirs(os.path.join(o, "tables"), exist_ok=True)
        p = os.path.join(o, "tables", "de_provenance.json")
        with open(p, "w") as f:
            f.write('{"adjp": 0.05}\n')
        return sp.sha256_file(p)
    ok = {"schema": sp.SCHEMA, "status": "ok"}
    if name == "no folder":
        return
    if name == "empty folder":
        os.makedirs(sky)
        return
    raw = {"garbage": "not json {", "a list": "[1, 2]", "null": "null"}.get(name)
    if raw is not None:
        write_record(o, raw=raw)
        return
    rec = {"skipped with a reason": {"schema": sp.SCHEMA, "status": "skipped",
                                     "reason": "the container is gone"},
           "skipped without a reason": {"schema": sp.SCHEMA, "status": "skipped"},
           "skipped with an empty reason": {"schema": sp.SCHEMA, "status": "skipped",
                                            "reason": ""},
           "planned": {"schema": sp.SCHEMA, "status": "planned", "plan_id": "p"},
           "submitted": {"schema": sp.SCHEMA, "status": "planned", "plan_id": "p"},
           "ok": ok,
           "ok and stale": dict(ok, de_record_sha256="0" * 64),
           "ok made with no DE, a DE now": dict(ok, de_record_sha256=None),
           "an unknown status": {"schema": sp.SCHEMA, "status": "weird"},
           "no status": {"schema": sp.SCHEMA}}.get(name)
    if name == "ok and current":
        rec = dict(ok, de_record_sha256=de())
    if rec is None:
        raise KeyError(name)
    write_record(o, rec)
    if name == "submitted":
        os.makedirs(os.path.join(sky, sp.WORK))
        with open(os.path.join(sky, sp.WORK, "jobs.txt"), "w") as f:
            f.write("build\t1\n")
    if name == "ok made with no DE, a DE now":
        de()


class FrozenOldRule(unittest.TestCase):
    """skyline_peaks.deliverable() and its core_submission wrapper give, on every one of the
    16 folders, exactly what the old core_submission rule gave (L9 review F10)."""

    def test_sixteen_cases(self):
        self.assertEqual(len(FROZEN), 16)
        for name, want in FROZEN:
            with self.subTest(case=name), tempfile.TemporaryDirectory() as d:
                build_case(d, name)
                ok, state, why = sp.deliverable(d)
                self.assertEqual((ok, why), want)
                self.assertEqual(cs.skyline_deliverable(d), want)
                self.assertEqual(ok, state == "ok")

if __name__ == "__main__":
    unittest.main()
