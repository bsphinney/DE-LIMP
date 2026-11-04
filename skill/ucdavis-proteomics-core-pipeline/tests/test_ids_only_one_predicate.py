#!/usr/bin/env python3
"""One answer to "was this search for identifications only?" (2.14.2 C2 follow-up; DE-LIMP rule 3).

Two readings existed after the C2 merge:
  * sage_lfq_check.search_quantified(folder) -- THE predicate the DE gate, the README, the email,
    the report page and the analysis brief ask: the folder's sage_lfq_check.json, else the copies
    in search_provenance.json (`quantified: false`, `sage_lfq_check` / `quantification` status
    not_applicable);
  * provenance.py's own `_record_key(record, "quantified") is False` for reproduce.sh's step-5
    notes and REPRODUCE.md.
They agreed on every record run_search.py writes today, but not on a record with only one of the
signs (an older Sage record whose folder's sage_lfq_check.json says not_applicable). Now
provenance.identified_only() asks search_quantified (with the record it read), so there is one
definition. What these tests pin:
  * the two readers agree, and are right, on an identifications-only Sage record, an
    identifications-only FragPipe record, a quantifying DIA-NN record and a quantifying Sage
    record -- and on the older Sage record, where provenance.py now follows the predicate;
  * no record: not identifications-only;
  * the bundle says "identifications only" exactly when the predicate does.
Stdlib only. Fixtures: test_reproduce_replays_what_ran.py (fake setup.json, no conda, no R).
"""
import json
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import provenance                                                 # noqa: E402
import sage_lfq_check as slc                                      # noqa: E402
import test_reproduce_replays_what_ran as rw                      # noqa: E402

write = rw.write
NA = {"status": "not_applicable", "why": "identifications only (quant.lfq is false)"}


class Cases(rw.Base):
    """Each case: a search folder with its record (and, for some, the folder's own
    sage_lfq_check.json) -> (folder, record path, record, identifications only?)."""

    def folder(self, name, record, gate_record=None, files=()):
        d = os.path.join(self.d, "output", name)
        for f in files:
            write(os.path.join(d, f), "x")
        if gate_record is not None:
            write(os.path.join(d, slc.RECORD), json.dumps(gate_record))
        path = write(os.path.join(d, slc.PROVENANCE), json.dumps(record))
        return d, path, record

    def cases(self):
        return {
            "sage identifications only": self.folder("search_sage_ids", {
                "engine": "sage", "version": "0.14.7", "quantified": False,
                "sage_lfq_check": NA, "sage_identifications": {"row_unit": "r", "kept_rule": "k"},
                "result": {"engine": "sage", "quantified": False}},
                gate_record=NA, files=("identifications.tsv", "sage_identifications.json"))
            + (True,),
            "fragpipe identifications only": self.folder("search_fp_ids", {
                "engine": "fragpipe", "version": "24.0", "quantified": False,
                "quantification": dict(NA, why="the workflow does not quantify"),
                "result": {"engine": "fragpipe", "quantified": False}},
                files=("identifications.tsv", "fragpipe_identifications.json")) + (True,),
            "diann quantifying": self.folder("search_diann", {
                "engine": "diann", "version": "2.7.0", "result": {"engine": "diann"}},
                files=("report.parquet",)) + (False,),
            "sage quantifying": self.folder("search_sage_lfq", {
                "engine": "sage", "version": "0.14.7", "quantified": True,
                "sage_lfq_check": {"status": "ok"}, "result": {"engine": "sage"}},
                gate_record={"status": "ok"}, files=("report.parquet",)) + (False,),
            # an older record without `quantified`: only the folder's gate record says it
            "older sage record, the gate's record says so": self.folder("search_sage_old", {
                "engine": "sage", "version": "0.14.7", "result": {"engine": "sage"}},
                gate_record=NA, files=("identifications.tsv",)) + (True,),
        }


class TheReadersAgree(Cases):
    def test_on_every_record(self):
        for name, (folder, path, rec, ids_only) in self.cases().items():
            with self.subTest(case=name):
                self.assertIs(slc.search_quantified(folder), not ids_only)
                self.assertIs(provenance.identified_only(path, rec), ids_only)

    def test_no_record_is_not_identifications_only(self):
        self.assertIs(provenance.identified_only(None, None), False)
        self.assertIs(provenance.identified_only(os.path.join(self.d, "nowhere.json"), None),
                      False)


class TheBundleSaysWhatThePredicateSays(Cases):
    def test_reproduce_sh_and_reproduce_md(self):
        for name, (_folder, path, _rec, ids_only) in self.cases().items():
            with self.subTest(case=name):
                self.out = os.path.join(self.d, "repro_" + name.replace(" ", "_").replace(",", ""))
                sh, md, _man = self.provenance("--search-record", path, "--engine",
                                               _rec["engine"])
                self.assertIs("Identifications only (its record says it did not quantify)"
                              in sh, ids_only)
                self.assertIs("Identifications only: the search did not quantify" in md,
                              ids_only)
                self.assertIs("No differential expression ran: the search's report is the "
                              "result." in md, not ids_only)


if __name__ == "__main__":
    unittest.main()
