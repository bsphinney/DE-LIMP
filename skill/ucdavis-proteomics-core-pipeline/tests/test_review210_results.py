#!/usr/bin/env python3
"""
Independent review of 2.10 (review210-results): the normalisation decision -- ways a client can
get the wrong quantities without anyone noticing. Each test FAILS on 40c5b9f.

  1. The check is skippable: the chain's recovery step (diann_parallel.py --next, which a resumed
     session reads from RECOVERY.md / checkpoint.py status) runs the final DE with no check, and
     run_de.R without --normalization-check reads normalised quantities whatever the session's
     recorded type says -- an IP gets DIA-NN's normalised quantities, flagged only as INFO.
  2. The gate is not tied to the data: a decision made on one report authorises a DE of another
     report (and the DE then records the other report's check as its own).
  3. The recommendation ignores the experiment type: B1 alone makes it say "raw" for a
     whole-proteome cohort whose raw volcano is the broken one, and "normalised" for an IP; for a
     TurboID with low controls it says raw and, in the same sentence, that normalisation corrects.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402
from test_normalization_check import HAVE_ARROW, NEEDS, make_report, session  # noqa: E402
from test_run_de_contaminants import r_has  # noqa: E402

CHECK = os.path.join(SCRIPTS, "normalization_check.py")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")
GENERATOR = os.path.join(SCRIPTS, "diann_parallel.py")
PINNED = "--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n--cont-quant-exclude Cont_\n"


class TheCheckIsNotSkippable(unittest.TestCase):
    def test_the_recovery_next_step_runs_the_normalisation_check(self):
        """diann_parallel.py:1893 records `--next "Rscript run_de.R ... --method dpc ..."` (no
        check, no --quantities) -- the command RECOVERY.md and `checkpoint.py status`
        (next_commands) hand a session that resumes after the chain finishes."""
        d = os.path.realpath(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, d, True)
        sess = os.path.join(d, "sess")
        os.makedirs(os.path.join(sess, "input"))
        raws = []
        for i in range(4):
            raws.append(os.path.join(d, f"f{i}.d"))
            os.makedirs(raws[-1])
        with open(os.path.join(d, "db.fasta"), "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDER\n")
        with open(os.path.join(d, "params.cfg"), "w") as fh:
            fh.write(PINNED)
        out = os.path.join(sess, "output", "search")
        p = subprocess.run([sys.executable, GENERATOR, "--diann", "/bin/true", "--raw", *raws,
                            "--fasta", os.path.join(d, "db.fasta"), "--out", out,
                            "--cfg", os.path.join(d, "params.cfg")],
                           capture_output=True, text=True, cwd=d, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(out, "submit.sh")) as fh:
            nxt = [ln for ln in fh.read().splitlines() if "--next" in ln]
        self.assertTrue(nxt, "no --next recorded")
        self.assertIn("normalization_check", nxt[0],
                      "the recovery step goes straight to the final DE: " + nxt[0].strip())

    @unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
    def test_a_de_with_no_check_in_an_ip_session_is_not_silently_normalised(self):
        """The session records an IP (default non-normalised). A DE run without
        --normalization-check -- the recovery command above, or any re-run -- reads
        Precursor.Normalised; the record says status not_run, the audit says INFO."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        report, cond, contrast = make_report(os.path.join(d, "ip"), "ip")
        s = session(os.path.join(d, "ip"), "ip", bait="BAITX")
        meta = os.path.join(s, "input", "conditions.csv")
        shutil.copy(cond, meta)                       # where the chain's DE reads it
        out = os.path.join(s, "output", "tables")
        r = subprocess.run(["Rscript", RUN_DE, "--input", report, "--metadata", meta,
                            "--contrasts", contrast, "--outdir", out],
                           capture_output=True, text=True, timeout=900, cwd=d)
        if r.returncode != 0:
            return                                    # refused: not silent
        with open(os.path.join(out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov.get("quantities"), et.RAW,
                         "an IP session's DE read DIA-NN's NORMALISED quantities with no check: "
                         f"{prov.get('normalisation')!r}, check status "
                         f"{(prov.get('normalization_check') or {}).get('status')!r}")


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class TheGateIsTiedToTheData(unittest.TestCase):
    def test_a_check_of_one_report_does_not_authorise_a_de_of_another(self):
        """normalization_check.gate(path, quantities) never looks at the DE's --input or
        --metadata; the record names its own report and conditions (rec["report"],
        rec["conditions"]). A decision on search v1 then passes the DE of search v2 (or of
        corrected conditions), and de_provenance/Methods say "a check of the data agreed"."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        report, cond, contrast = make_report(os.path.join(d, "ip"), "ip")
        s = session(os.path.join(d, "ip"), "ip", bait="BAITX")
        dirs = {}
        for q in (et.NORMALISED, et.RAW):
            dirs[q] = os.path.join(d, "ip", "norm_check", q)
            p = subprocess.run(["Rscript", RUN_DE, "--input", report, "--metadata", cond,
                                "--contrasts", contrast, "--quantities", q, "--outdir", dirs[q]],
                               capture_output=True, text=True, timeout=900, cwd=d)
            self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        p = subprocess.run([sys.executable, CHECK, "run", "--report", report, "--conditions", cond,
                            "--de-normalised", dirs[et.NORMALISED], "--de-raw", dirs[et.RAW],
                            "--session", s], capture_output=True, text=True, timeout=300)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        path = json.loads(p.stdout)["check"]
        # another report (a re-search, another project) with the same decision's quantities
        other, ocond, ocontrast = make_report(os.path.join(d, "other"), "whole", seed=11)
        out = os.path.join(d, "final")
        r = subprocess.run(["Rscript", RUN_DE, "--input", other, "--metadata", ocond,
                            "--contrasts", ocontrast, "--quantities", et.RAW,
                            "--normalization-check", path, "--outdir", out],
                           capture_output=True, text=True, timeout=900, cwd=d)
        if r.returncode == 0:
            with open(os.path.join(out, "de_provenance.json")) as fh:
                n = json.load(fh)["normalization_check"]
            self.fail(f"a check of {n.get('report')} was accepted for a DE of {other} "
                      f"(recorded as status {n.get('status')!r}, decided by "
                      f"{(n.get('decision') or {}).get('by')!r})")


def rec_for(etype, fc, anomalies, carbox=None):
    """A tripped record as normalization_check.run() builds it, for recommend()."""
    d = et.default_for(etype)
    vol = {k: {"contrasts": {"A-B": {"anomalies": [f"x{i}" for i in range(n)]}}}
           for k, n in anomalies.items()}
    return {"experiment_type": {"type": etype}, "default": {"quantities": d["quantities"],
                                                            "why": d["why"]},
            "checks": {"factor_confound": fc}, "volcano": vol, "tripped": True,
            "trips": ["B1 ..."], "carboxylases": carbox}


class TheRecommendationRespectsTheType(unittest.TestCase):
    CONF = {"large_and_confounded": True, "large_unexplained": False, "gap_fold": 2.1}
    LOOSE = {"large_and_confounded": False, "large_unexplained": True, "gap_fold": 1.4}

    def test_a_whole_proteome_is_not_told_raw_when_the_raw_volcano_is_the_broken_one(self):
        """A whole-proteome 2 vs 2 whose loading happens to follow the groups (B1: 2-fold,
        separated -- 4-13% of 2 vs 2 cohorts with 0.5-0.7 log2 loading noise, simulated) and
        whose RAW volcano is shifted and lopsided: the recommendation still says raw."""
        rec = rec_for("whole_proteome", self.CONF, {et.NORMALISED: 0, et.RAW: 3})
        r = nc.recommend(rec)
        self.assertFalse(r.startswith("non-normalised (raw) quantities"),
                         "recommended raw for a whole proteome against its own volcanoes: " + r)

    def test_an_ip_is_not_told_to_normalise_when_the_normalised_volcano_is_the_broken_one(self):
        """An IP with one IgG lane carrying more background (factors not separating): B1
        'unexplained' and the recommendation is to normalise -- the thing the IP default
        exists to prevent -- though the normalised volcano shows depleted background."""
        rec = rec_for("ip", self.LOOSE, {et.NORMALISED: 2, et.RAW: 0})
        r = nc.recommend(rec)
        self.assertFalse(r.startswith("normalised quantities"),
                         "recommended normalising an IP against its own volcanoes: " + r)

    def test_separation_fractions_are_not_steered_to_normalisation_by_the_volcano(self):
        """SEC / organelle fractions (default non-normalised, per DIA-NN's README): their raw
        volcano is shifted and lopsided BY DESIGN (each fraction holds a different proteome),
        and only enrichment contrasts are exempted (volcano_version raw=True) -- so the check
        trips on the raw volcano and the anomaly count recommends normalising."""
        none = {"large_and_confounded": False, "large_unexplained": False}
        rec = rec_for("fractions_separation", none, {et.NORMALISED: 0, et.RAW: 3})
        r = nc.recommend(rec)
        self.assertFalse(r.startswith("normalised quantities"),
                         "recommended normalising separation fractions: " + r)

    def test_a_turboid_recommendation_does_not_contradict_itself(self):
        """test_turboid_whose_controls_came_out_low_still_trips: proximity (default normalised),
        controls low -> 'non-normalised (raw) quantities ...' followed by '... consistent with
        the controls carrying less material, which normalisation corrects'."""
        cx = {et.NORMALISED: {"median_sd_log2": 0.1}, et.RAW: {"median_sd_log2": 0.9}}
        rec = rec_for("proximity", self.CONF, {et.NORMALISED: 1, et.RAW: 1}, cx)
        r = nc.recommend(rec)
        self.assertFalse(r.startswith("non-normalised (raw) quantities")
                         and "which normalisation corrects" in r, r)


if __name__ == "__main__":
    unittest.main()
