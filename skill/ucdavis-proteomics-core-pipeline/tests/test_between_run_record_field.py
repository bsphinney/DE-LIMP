#!/usr/bin/env python3
"""The between-run step is read from the record's own `between_run` field (C1 review, F4; rule 1).

provenance.py (reproduce.sh's --no-norm notes) and normalization_check.step8_commands branched on
`method == "maxlfq"`: a DE route other than MaxLFQ that applies and records a between-run step
would have needed edits in both. Now they key off the record's `between_run` (de_provenance.json,
build_maxlfq.R maxlfq_between_run; a design_request's requested value) through
normalization_check.between_run_of. The one method test left is there, for what a record from
before 2.14.2 -- which names no step -- meant: quantile for normalised MaxLFQ, none for raw.

Stdlib only: reproduce.sh's steps 5 and 6 run with stub run_search.py / Rscript.
"""
import json
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402
from test_reproduce_sh_runs import arg, provenance, run_steps_5_and_6  # noqa: E402

# a route that is not MaxLFQ but applies, and records, a between-run step of its own
NEW_ROUTE = "another_route"


class BetweenRunOf(unittest.TestCase):
    def test_a_record_that_names_its_step_is_read_whatever_its_route(self):
        rec = {"method": NEW_ROUTE, "quantities": et.NORMALISED,
               "between_run": {"method": "median", "requested": "median", "why": "x"}}
        self.assertEqual(nc.between_run_of(rec, "method"), "median")
        self.assertEqual(nc.between_run_of({"method": NEW_ROUTE, "between_run": "auto"}), "auto")

    def test_a_record_from_before_2_14_2_is_read_by_the_old_rule(self):
        """(pin) the legacy rule, unchanged: the one place that decides what an old record meant."""
        self.assertEqual(nc.between_run_of({"method": "maxlfq", "quantities": et.NORMALISED},
                                           "method"), nc.LEGACY_MAXLFQ_BETWEEN_RUN)
        self.assertEqual(nc.between_run_of({"method": "maxlfq", "quantities": et.RAW}), "none")
        self.assertIsNone(nc.between_run_of({"method": "dpc"}))
        # run_de.R writes a dpc DE's between_run as null
        self.assertIsNone(nc.between_run_of({"method": "dpc", "between_run": None}))


class Step8Commands(unittest.TestCase):
    def test_the_opt_in_is_read_from_the_design_not_the_method(self):
        cmds = nc.step8_commands("/s/output/search/report.parquet", "/s/input/conditions.csv",
                                 "/s", method=NEW_ROUTE,
                                 design={"method": NEW_ROUTE, "between_run": "median"})
        norm, raw, _ = cmds["check"]
        self.assertIn("--between-run median", norm)
        self.assertNotIn("--between-run", raw)
        self.assertIn("--between-run median", cmds["final"])
        dpc = nc.step8_commands("/s/output/search/report.parquet", "/s/input/conditions.csv",
                                "/s", design={"method": "dpc", "between_run": None})
        self.assertNotIn("--between-run", dpc["next"])


class ReproduceSh(unittest.TestCase):
    def replay(self, rec):
        """-> (the replayed --between-run, the replay's stderr)."""
        with tempfile.TemporaryDirectory() as d:
            de = os.path.join(d, "tables")
            os.makedirs(de)
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump(rec, fh)
            out = os.path.join(d, "repro")
            sh = provenance(out, de, "--de-method", rec.get("method") or "dpc")
            calls, err = run_steps_5_and_6(d, out, sh)
        return arg(calls[1]["argv"], "--between-run"), err

    def test_a_no_norm_reports_normalised_de_is_noted_for_any_route_that_records_a_step(self):
        flag, err = self.replay({"method": NEW_ROUTE, "quantities": et.NORMALISED,
                                 "between_run": {"method": "quantile", "requested": "auto",
                                                 "why": "x"},
                                 "dia_nn_normalisation": {"no_norm": True, "source": "x"}})
        self.assertEqual(flag, "quantile")
        self.assertIn("re-quantify with --no-norm", err)

    def test_a_raw_de_of_a_route_that_records_a_step_is_noted(self):
        flag, err = self.replay({"method": NEW_ROUTE, "quantities": et.RAW,
                                 "between_run": {"method": "none", "requested": "auto",
                                                 "why": "--quantities raw"}})
        self.assertIsNone(flag)
        self.assertIn("re-quantify with --no-norm", err)

    def test_a_2_14_1_maxlfq_record_and_a_dpc_record(self):
        """(pin) a 2.14.1 MaxLFQ record replays the old rule's steps and notes; dpc none."""
        self.assertEqual(self.replay({"method": "maxlfq", "quantities": et.NORMALISED})[0],
                         "quantile")
        flag, err = self.replay({"method": "maxlfq", "quantities": et.RAW})
        self.assertIsNone(flag)
        self.assertIn("re-quantify with --no-norm", err)
        flag, err = self.replay({"method": "dpc", "quantities": et.RAW, "between_run": None})
        self.assertIsNone(flag)
        self.assertNotIn("re-quantify", err)


if __name__ == "__main__":
    unittest.main()
