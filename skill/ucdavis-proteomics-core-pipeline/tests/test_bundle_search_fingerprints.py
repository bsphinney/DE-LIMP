#!/usr/bin/env python3
"""The reproducibility bundle fingerprints what the search wrote (2.14.2 C2, E2E F3).

finalize gave provenance.py only `--report` -- the DE's input or the record's result.report --
and an identifications-only Sage search has neither. So its bundle said "left out for want of a
record: the search report ... so the bundle fingerprints none", checksums.json held only the raw
inputs and the FASTA, and REPRODUCE.md still claimed checksums/ held the "search report".

Now provenance.py takes --search-output (repeatable) and fingerprints the search record it reads
(--search-record); finalize passes session.search_result_files(engine_tables=True): the
identifications and their record, and Sage's own results.sage.parquet. What these tests pin:
  * checksums.json `search_outputs` and `search_record` carry each file's sha256;
  * REPRODUCE.md says what checksums/ holds, from what was fingerprinted;
  * finalize's command (make_deposit.bundle_command) names those outputs and says nothing is
    missing for an identifications-only search;
  * finalize's build puts them in the bundle;
  * a re-adapted identifications.tsv makes the bundle stale (report_files.bundle_state), as a
    changed report already did;
  * a quantifying DIA-NN session's command is unchanged: its report, no other output (pin).
Stdlib only. Fixtures: test_reproduce_replays_what_ran.py.
"""
import hashlib
import json
import os
import sys
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

from job_env import job_env                                       # noqa: E402
from test_reproduce_sh_runs import arg                            # noqa: E402
import test_reproduce_replays_what_ran as rw                      # noqa: E402

read = rw.read


def sha(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def search_outputs(cmd):
    return [cmd[i + 1] for i, w in enumerate(cmd) if w == "--search-output"]


class ProvenanceFingerprints(rw.Base):
    def test_the_search_outputs_and_record_are_fingerprinted(self):
        cfg, rec = self.sage_ids_only()
        search = os.path.dirname(rec)
        outs = [os.path.join(search, n) for n in ("identifications.tsv",
                                                   "sage_identifications.json",
                                                   "results.sage.parquet")]
        more = []
        for o in outs:
            more += ["--search-output", o]
        _sh, md, man = self.provenance("--workflow-manifest", self.manifest(), "--params", cfg,
                                       "--search-record", rec, "--engine", "sage", *more)
        checks = json.loads(read(os.path.join(self.out, "checksums", "checksums.json")))
        self.assertEqual([(x["path"], x["sha256"]) for x in checks["search_outputs"]],
                         [(o, sha(o)) for o in outs])
        self.assertEqual((checks["search_record"]["path"], checks["search_record"]["sha256"]),
                         (rec, sha(rec)))
        self.assertNotIn("report", checks, "an identifications-only search has no report")
        self.assertIn("[OK]      fingerprinted 3 search output(s) (identifications.tsv, "
                      "sage_identifications.json, results.sage.parquet)", man)
        self.assertIn("[OK]      fingerprinted the search record (search_provenance.json)", man)
        captured = md.split("## What's captured", 1)[1].split("##", 1)[0]
        self.assertIn("the search's outputs (identifications.tsv, sage_identifications.json, "
                      "results.sage.parquet)", captured)
        self.assertIn("the search record", captured)
        self.assertNotIn("search report", captured)
        self.assertNotIn("DE outputs", captured)

    def test_a_missing_output_is_said(self):
        cfg, rec = self.sage_ids_only()
        gone = os.path.join(os.path.dirname(rec), "nothing_here.tsv")
        _sh, _md, man = self.provenance("--workflow-manifest", self.manifest(), "--params", cfg,
                                        "--search-record", rec, "--search-output", gone)
        self.assertIn(f"[SKIPPED] search output -- not found: {gone}", man)


class FinalizePassesThem(rw.Base):
    session = rw.FinalizesCommand.session        # the identifications-only Sage session

    def test_the_command_names_the_search_outputs(self):
        import make_deposit
        p, _cfg, ids = self.session()
        cmd, missing, unnamed = make_deposit.bundle_command(make_deposit.gather(p["session_dir"]))
        self.assertEqual(missing, [])
        self.assertEqual(search_outputs(cmd),
                         [ids["identifications"],
                          os.path.join(p["search_out"], "sage_identifications.json"),
                          ids["source"]])
        self.assertFalse([u for u in unnamed if "fingerprints none" in u], unnamed)

    def test_finalizes_bundle_holds_them_and_goes_stale_when_they_change(self):
        import make_deposit
        import report_files as rf
        p, _cfg, ids = self.session()
        env = job_env(self.d, PROTEOMICS_PIPELINE_HOME=self.home)
        with mock.patch.dict(os.environ, env, clear=True):
            man = make_deposit.Manifest()
            self.assertTrue(make_deposit.ensure_bundle(p["session_dir"], man), man.lines)
        self.assertNotIn("so the bundle fingerprints none", " ".join(man.lines))
        checks = json.loads(read(os.path.join(p["repro_dir"], "checksums", "checksums.json")))
        self.assertEqual([os.path.basename(x["path"]) for x in checks["search_outputs"]],
                         ["identifications.tsv", "sage_identifications.json",
                          "results.sage.parquet"])
        self.assertEqual(checks["search_record"]["sha256"], sha(p["search_prov"]))
        self.assertEqual(rf.bundle_state(p["output_dir"])[0], "present")
        with open(ids["identifications"], "a") as fh:      # --adapt-only run again
            fh.write("runA\tPEPTIDEK\tP1\t2\n")
        state, why = rf.bundle_state(p["output_dir"])
        self.assertEqual(state, "stale")
        self.assertIn("the search output identifications.tsv changed", why)


class QuantifyingSession(rw.Base):
    def test_a_diann_sessions_command_is_unchanged(self):
        """(pin) the report, as --report; no other search output."""
        import make_deposit
        import test_deposit_package as tdp
        root = os.path.join(self.d, "s")
        os.makedirs(root)
        p = tdp.dia_session(root, de=False)
        cmd, missing, unnamed = make_deposit.bundle_command(make_deposit.gather(p["session_dir"]))
        self.assertEqual(missing, [])
        self.assertEqual(arg(cmd, "--report"), os.path.join(p["search_out"], "report.parquet"))
        self.assertEqual(search_outputs(cmd), [])
        self.assertEqual(unnamed, [])


if __name__ == "__main__":
    unittest.main()
