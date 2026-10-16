#!/usr/bin/env python3
"""
The SILAC ratio analysis (run_silac.R + silac_ratios.R) on a synthetic plexDIA report
(tests/silac_fixture.py: planted effects, a mixing offset, channel-specific proteins).

The method was validated on a real triple-SILAC experiment (4 replicates, DIA-NN 2.7.0): it
reproduced the hand-made analysis's tables to floating-point precision (1,224 / 1,153 / 1,093
proteins tested; 816 / 799 / 0 significant). These tests pin, on data whose answer is known:
  * the planted changes are found, in the right direction, and nothing else;
  * the ratio maths, recomputed here in Python: the q-value filter (every column the report has,
    or the run-specific two with --run-level-only), the precursor median, the per-run centring,
    the mean over a replicate's re-injections, the limma mean, Benjamini-Hochberg, the minimum
    number of replicates;
  * the statistician's checks: the zero point (median against the density peak), single-precursor
    calls, the realised BH cutoff behind the power reading, detection-pattern wording;
  * every default is tagged in the record and in the Methods;
  * a label swap compares groups, not channels; a label-check run gives the incorporation;
  * a search that declared the label twice is flagged; channel imbalance is flagged;
  * figures.json is the contract make_analysis_html.py / analysis_prompt.py read;
  * de_provenance.json describes itself (pipeline_id silac_ratio, its normalisation check not
    applicable, no set tests, its replay, its Methods paragraph);
  * another pipeline's results are never overwritten; a turnover search is refused;
  * reproducibility_log.R replays the tables exactly;
  * the design's rules (silac.py read_design).
"""
import csv
import json
import math
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import silac  # noqa: E402
import silac_fixture as fx  # noqa: E402

RSCRIPT = shutil.which("Rscript")


def r_has(*pkgs):
    if not RSCRIPT:
        return False
    code = "q(status = if (all(sapply(c(%s), requireNamespace, quietly = TRUE))) 0 else 1)" % \
        ", ".join(f'"{p}"' for p in pkgs)
    return subprocess.run([RSCRIPT, "-e", code], capture_output=True).returncode == 0


def py_has(mod):
    try:
        __import__(mod)
        return True
    except ImportError:
        return False


# R runs the analysis; the synthetic report (silac_fixture.make_report) is written with pyarrow.
# CI's unittest stage has neither, and skips these.
HAVE_R = r_has("limma", "data.table", "arrow", "jsonlite") and py_has("pyarrow")
NEEDS = "needs Rscript with limma, data.table, arrow and jsonlite, and Python's pyarrow"


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def num(x):
    return float("nan") if x in ("", "NA", None) else float(x)


GROUP = {"L": "WT", "M": "KO", "H": "DKO"}
RUN_LEVEL = ("Q.Value", "PG.Q.Value")
EXPERIMENT_WIDE = ("Global.Q.Value", "Global.PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value")


def hand_ratios(report, num_g, den_g, rep_of, run_level_only=False):
    """The ratio maths by hand: -> (per-run protein medians {(run, pg): [log2 ratios]},
    per-run offsets {run: median}, per-replicate centred means {(rep, pg): value})."""
    import pyarrow.parquet as pq
    cols = RUN_LEVEL if run_level_only else RUN_LEVEL + EXPERIMENT_WIDE
    cut = {c: 0.05 if c == "PG.Q.Value" else 0.01 for c in cols}
    q = {}
    for r in pq.read_table(report).to_pylist():
        if (r["Channel"] in GROUP and r["Run"] in rep_of and r["Channel.Q.Value"] <= 0.01
                and all(r[c] <= v for c, v in cut.items())
                and not r["Protein.Ids"].startswith("Cont_")):
            q.setdefault((r["Run"], r["Protein.Group"], r["Precursor.Id"]), {})[
                GROUP[r["Channel"]]] = r["Precursor.Quantity"]
    prot = {}
    for (run, pg, _), v in q.items():
        if num_g in v and den_g in v:
            prot.setdefault((run, pg), []).append(math.log2(v[num_g] / v[den_g]))
    med = {k: statistics.median(v) for k, v in prot.items()}
    offset = {run: statistics.median([m for (r, _), m in med.items() if r == run])
              for run in {r for r, _ in med}}
    per_rep = {}
    for (run, pg), m in med.items():
        per_rep.setdefault((rep_of[run], pg), []).append(m - offset[run])
    return prot, offset, {k: statistics.mean(v) for k, v in per_rep.items()}


def run_silac(d, *more, report=None, design=None, ok=True, outdir=None):
    outdir = outdir or os.path.join(d, "out", "tables")
    p = subprocess.run([RSCRIPT, os.path.join(SCRIPTS, "run_silac.R"),
                        "--input", report or os.path.join(d, "report.parquet"),
                        "--design", design or os.path.join(d, "design.csv"),
                        "--outdir", outdir, *more], capture_output=True, text=True, timeout=600)
    if ok:
        assert p.returncode == 0, p.stdout + p.stderr
    return p, outdir


@unittest.skipUnless(HAVE_R, NEEDS)
class Analysis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.d = cls.tmp.name
        fx.make_report(os.path.join(cls.d, "report.parquet"))
        fx.write_design(os.path.join(cls.d, "design.csv"))
        with open(os.path.join(cls.d, "report.log.txt"), "w") as fh:   # DIA-NN 2.7.0's wording
            fh.write("Registered SILAC channel H with following masses: K[8.0142] R[10.0083] \n"
                     "Registered SILAC channel decoy with following masses: K[12.0033] R[13.9964] \n")
        cls.p, cls.out = run_silac(cls.d)
        cls.prov = json.load(open(os.path.join(cls.out, "de_provenance.json")))
        cls.ko = {r["Protein.Group"]: r for r in read(os.path.join(cls.out, "DE_silac_KO.WT.csv"))}

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_planted_changes_are_found_and_nothing_else(self):
        sig = {k: r for k, r in self.ko.items() if num(r["adj.P.Val"]) < 0.05}
        self.assertEqual(sorted(sig), [f"P{i:02d}" for i in range(1, 21)])
        for k, r in sig.items():
            self.assertEqual(num(r["logFC"]) > 0, int(k[1:]) <= 10, k)
            self.assertAlmostEqual(abs(num(r["logFC"])), 2.0, delta=0.25)
        dko_ko = read(os.path.join(self.out, "DE_silac_DKO.KO.csv"))
        self.assertFalse([r for r in dko_ko if num(r["adj.P.Val"]) < 0.05])
        self.assertEqual(self.prov["contrasts"], ["KO-WT", "DKO-WT", "DKO-KO"])
        self.assertIn("DEFAULT", self.prov["comparisons_source"])

    def test_untested_channel_specific_and_contaminants(self):
        self.assertIn("ratio in 2 of 4 replicates: not tested (minimum 3)", self.ko["P56"]["Evidence"])
        self.assertTrue(math.isnan(num(self.ko["P56"]["adj.P.Val"])))
        for k in ("P51", "P55", "P57", "Cont_P00761"):
            self.assertNotIn(k, self.ko)
        spec = {r["Protein.Group"]: r for r in read(os.path.join(self.out,
                                                                 "SILAC_channel_specific_detections.csv"))}
        self.assertEqual(sorted(spec), [f"P{i}" for i in range(51, 56)])
        self.assertEqual(spec["P51"]["pattern"],
                         "detected in WT (4/4); not detected in KO (0/4), DKO (0/4)")
        # where it is detected: its abundance among that channel's proteins; elsewhere nothing
        self.assertTrue(0 < num(spec["P51"]["percentile_WT"]) <= 100)
        self.assertEqual(spec["P51"]["percentile_KO"], "")
        self.assertEqual(self.prov["filter_counts"]["contaminant_protein_groups"], 1)
        self.assertEqual(self.ko["P01"]["Detected_KO"], "4/4")
        # P57: in WT in 2 replicates -- no ratio, too few for a pattern: only in the sample matrix
        em = {r["Protein.Group"] for r in read(os.path.join(self.out, "SILAC_sample_matrix.csv"))}
        self.assertIn("P57", em)
        self.assertEqual(self.prov["sample_matrix"]["proteins_only_here"], 1)
        chk = [c for c in self.prov["checks"] if c["check"] == "silac_unlisted_proteins"]
        self.assertIn("1 protein groups", chk[0]["message"])

    def test_ratio_maths_by_hand(self):
        rep = dict(zip(fx.RUNS, fx.REPS))
        prot, offset, _ = hand_ratios(os.path.join(self.d, "report.parquet"), "KO", "WT", rep)
        med = {k: statistics.median(v) for k, v in prot.items()}
        for run, rp in rep.items():
            self.assertAlmostEqual(offset[run], self.prov["de_tables"]["KO-WT"]["offset_log2"][rp], places=9)
        for pg in ("P03", "P12", "P30", "P41", "P42", "P44"):
            row = self.ko[pg]
            for run, rp in rep.items():
                self.assertAlmostEqual(num(row[f"uncentred_{rp}"]), med[(run, pg)], places=9)
                self.assertAlmostEqual(num(row[f"centred_{rp}"]), med[(run, pg)] - offset[run], places=9)
                self.assertEqual(int(row[f"n_precursors_{rp}"]), len(prot[(run, pg)]))
            # intercept-only limma: the estimate is the mean of the centred values
            self.assertAlmostEqual(num(row["logFC"]),
                                   statistics.mean(num(row[f"centred_{rp}"]) for rp in fx.REPS), places=9)
        # the experiment-wide columns apply: P40 fails Global.PG.Q.Value, P41 / P42 lose the
        # precursor that fails Global.Q.Value / Lib.Q.Value
        self.assertNotIn("P40", self.ko)
        self.assertEqual(int(self.ko["P41"]["n_precursors_Rep1"]), 2, "3 searched, 1 fails")
        self.assertEqual(int(self.ko["P42"]["n_precursors_Rep1"]), 3, "4 searched, 1 fails")
        fc = self.prov["filter_counts"]
        self.assertEqual((fc["rows_failing_Global.PG.Q.Value"], fc["rows_failing_Global.Q.Value"],
                          fc["rows_failing_Lib.Q.Value"]), (24, 12, 12))
        self.assertEqual(self.prov["q_columns"],
                         ["Q.Value", "PG.Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value", "Global.Q.Value",
                          "Global.PG.Q.Value", "Channel.Q.Value"])
        self.assertEqual(self.prov["q_cutoffs"], [0.01, 0.05, 0.01, 0.01, 0.01, 0.01, 0.01])
        # Benjamini-Hochberg over the tested rows
        tested = sorted((num(r["P.Value"]), k) for k, r in self.ko.items() if r["Evidence"] == "tested")
        m, prev = len(tested), 1.0
        bh = {}
        for i in range(m - 1, -1, -1):
            prev = min(prev, tested[i][0] * m / (i + 1))
            bh[tested[i][1]] = prev
        for k, v in bh.items():
            self.assertAlmostEqual(num(self.ko[k]["adj.P.Val"]), v, places=12)
        self.assertEqual(self.prov["de_tables"]["KO-WT"]["n_tested"], m)

    def test_checks(self):
        by = {}
        for c in self.prov["checks"]:
            by.setdefault(c["check"], []).append(c)
        self.assertEqual(by["silac_channel_balance"][0]["status"], "WARN")
        self.assertEqual(by["silac_label_incorporation"][0]["status"], "INFO")
        self.assertNotIn("silac_label_declaration", by)
        # the zero point: median and density peak agree here (a symmetric design)
        self.assertEqual([c["status"] for c in by["silac_zero_point"]], ["PASS"] * 3)
        # P04's KO ratio rests on one precursor in Rep1: counted, said
        self.assertEqual(self.prov["de_tables"]["KO-WT"]["n_significant_one_precursor"], 1)
        self.assertEqual(num(self.ko["P04"]["min_precursors"]), 1)
        one = by["silac_single_precursor"][0]
        self.assertEqual(one["status"], "INFO")
        self.assertIn("KO-WT: 1 of the 20 significant", one["message"])
        # replicates not confirmed as separate cultures: said, tagged
        self.assertIn("DEFAULT", by["silac_replicate_independence"][0]["message"])
        pat = by["silac_detection_patterns"][0]["message"]
        self.assertIn("not a statistical test", pat)
        self.assertIn("not detected is not absent", pat)
        for c in self.prov["checks"]:
            self.assertEqual(set(c), {"check", "status", "message", "detail"})
            self.assertIn(c["status"], ("PASS", "WARN", "FAIL", "INFO"))

    def test_record_describes_itself(self):
        p = self.prov
        self.assertEqual(p["pipeline_id"], "silac_ratio")
        self.assertEqual(p["normalization_check"]["status"], "not_applicable")
        self.assertTrue(p["normalization_check"]["reason"])
        self.assertIs(p["set_tests"]["available"], False)
        self.assertEqual(p["figures"]["written_by"], "run_silac.R")
        self.assertEqual(p["replay"]["script"], "run_silac.R")
        args = p["replay"]["args"]
        self.assertEqual(args[args.index("--comparisons") + 1], "KO-WT,DKO-WT,DKO-KO")
        tag = "(DEFAULT \u2014 not user-confirmed)"
        for want in ("median detected protein", "KO/WT: 3 of 4", "robust = TRUE",
                     "Benjamini-Hochberg", "Derks", "Phipson", "no fold-change filter",
                     f"the run-specific Q.Value <= 0.01 {tag}",
                     f"the experiment-wide Global.PG.Q.Value <= 0.01 {tag}",
                     f"Channel.Q.Value <= 0.01 {tag}", f"adj.P.Val < 0.05 {tag}",
                     f"max(3, half the replicates) {tag}",
                     "The comparisons were KO/WT, DKO/WT, DKO/KO",
                     f"the group in the heavier channel over the lighter {tag}",
                     f"taken to be a separate culture and mix {tag}",
                     "not detected is not absent",
                     "decoy channel (K[12.0033] R[13.9964], from the search log)"):
            self.assertIn(want, p["methods_paragraph"])
        self.assertNotIn("1% global FDR", p["methods_paragraph"])
        self.assertEqual(p["de_tables"]["KO-WT"]["min_replicates"], 3)
        # every setting nobody gave says so, in its first word
        for k, v in p["settings"].items():
            self.assertIs(v["given"], False, k)
            self.assertTrue(v["source"].startswith("DEFAULT, not user-confirmed"), k)
        self.assertEqual(p["comparisons_sentence"].split(";")[0],
                         "The comparisons were KO/WT, DKO/WT, DKO/KO, each a log2 ratio "
                         "numerator/denominator (a positive value is higher in the first-named group)")
        self.assertEqual(p["decoy_channel"]["masses"], "K[12.0033] R[13.9964]")
        self.assertEqual(p["quantities"], "raw")
        self.assertTrue(any("MEDIAN-CENTRED" in n for n in p["brief_notes"]))
        self.assertEqual([s["title"] for s in p["report_sections"]],
                         ["SILAC channels and mixing", "Channel-specific proteins"])

    def test_figures_json_contract(self):
        from make_analysis_html import read_figures_json
        figdir = os.path.join(os.path.dirname(self.out), "figures")
        fj = read_figures_json(figdir)
        self.assertFalse(fj["error"])
        files = [f["file"] for f in fj["figures"]]
        for want in ("volcano_KO.WT.png", "volcano_DKO.WT.png", "volcano_DKO.KO.png",
                     "silac_channel_composition.png", "silac_channel_detections.png",
                     "qc_pvalue_panel.png"):
            self.assertIn(want, files)
        for f in fj["figures"]:
            self.assertTrue(os.path.getsize(os.path.join(figdir, f["file"])) > 0, f["file"])
            self.assertTrue(f["caption"], f["file"])
        self.assertEqual(fj["failed"], [])
        # the default scatter: two comparisons over the same reference (M/L against H/L), and
        # the correlation the shared channel builds in, with its sign
        self.assertIn("silac_scatter_KO.WT_vs_DKO.WT.png", files)
        scatter = next(f for f in fj["figures"] if f["file"].startswith("silac_scatter"))
        self.assertIn("share WT as the denominator in both", scatter["caption"])
        self.assertIn("a positive correlation", scatter["caption"])

    def test_power_reading_uses_the_realised_bh_cutoff(self):
        ko, dk = self.prov["de_tables"]["KO-WT"], self.prov["de_tables"]["DKO-KO"]
        self.assertAlmostEqual(ko["p_cutoff"], 0.05 * 20 / 52, places=12)
        self.assertIn("20 significant of 52 tested", ko["p_cutoff_basis"])
        self.assertAlmostEqual(dk["p_cutoff"], 0.05 / 52, places=12)
        self.assertIn("none of 52 tested significant", dk["p_cutoff_basis"])
        for t in (ko, dk):
            self.assertLess(t["critical_log2_raw_p05"], t["critical_log2_significant"])
            self.assertLess(t["critical_log2_significant"], t["true_log2_80pct_power"])
        self.assertNotIn("detectable_log2", ko)
        self.assertTrue(any("80% power" in n and "BH's realised cutoff" in n
                            for n in self.prov["brief_notes"]))

    def test_qc_matrix(self):
        em = read(os.path.join(self.out, "SILAC_sample_matrix.csv"))
        samples = read(os.path.join(self.out, "SILAC_samples.csv"))
        self.assertEqual([s["File.Name"] for s in samples][:3], ["Rep1_WT", "Rep1_KO", "Rep1_DKO"])
        col = [num(r["Rep1_KO"]) for r in em if r["Rep1_KO"] not in ("", "NA")]
        self.assertAlmostEqual(statistics.median(col), 0.0, places=9)

    def test_reproducibility_log_replays_exactly(self):
        p = subprocess.run([RSCRIPT, "reproducibility_log.R"], cwd=self.out, capture_output=True,
                           text=True, timeout=300)
        self.assertEqual(p.returncode, 0, p.stderr)
        replay = {r["Protein.Group"]: r for r in read(os.path.join(self.out,
                                                                   "replay_DE_silac_KO.WT.csv"))}
        for k, r in self.ko.items():
            self.assertEqual(r["logFC"], replay[k]["logFC"], k)
            self.assertEqual(r["adj.P.Val"], replay[k]["adj.P.Val"], k)
        for f in os.listdir(self.out):
            if f.startswith("replay_"):
                os.unlink(os.path.join(self.out, f))


@unittest.skipUnless(HAVE_R, NEEDS)
class Variants(unittest.TestCase):
    def test_min_replicates_is_a_recorded_parameter(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet")); fx.write_design(os.path.join(d, "design.csv"))
            _, out = run_silac(d, "--min-replicates", "4", "--scatter", "KO-WT:DKO-KO")
            ko = {r["Protein.Group"]: r for r in read(os.path.join(out, "DE_silac_KO.WT.csv"))}
            self.assertIn("ratio in 2 of 4", ko["P56"]["Evidence"])
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            self.assertEqual(prov["settings"]["min_replicates"],
                             {"value": 4, "given": True, "source": "given (--min-replicates)"})
            self.assertIn("(KO/WT: 4 of 4; DKO/WT: 4 of 4; DKO/KO: 4 of 4)", prov["methods_paragraph"])
            # a scatter of comparisons sharing a group on opposite sides: a negative coupling
            fj = json.load(open(os.path.join(os.path.dirname(out), "figures", "figures.json")))
            sc = next(f for f in fj["figures"] if f["file"].startswith("silac_scatter"))
            self.assertIn("a negative correlation", sc["caption"])
            self.assertTrue(any(c["check"] == "silac_replicates" and c["status"] == "WARN"
                                for c in prov["checks"]))

    def test_run_level_only_is_the_validated_hand_made_rule(self):
        # --run-level-only: the run-specific Q.Value and channel-specific PG.Q.Value alone, the
        # filter of the hand-made analysis the method reproduced; recorded and replayed
        with tempfile.TemporaryDirectory() as d:
            report = fx.make_report(os.path.join(d, "report.parquet"))
            fx.write_design(os.path.join(d, "design.csv"))
            _, out = run_silac(d, "--run-level-only")
            ko = {r["Protein.Group"]: r for r in read(os.path.join(out, "DE_silac_KO.WT.csv"))}
            rep = dict(zip(fx.RUNS, fx.REPS))
            prot, offset, centred = hand_ratios(report, "KO", "WT", rep, run_level_only=True)
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            for run, rp in rep.items():
                self.assertAlmostEqual(offset[run], prov["de_tables"]["KO-WT"]["offset_log2"][rp], places=9)
            for pg in ("P40", "P41", "P42"):
                self.assertEqual(ko[pg]["Evidence"], "tested", pg)
                for rp in fx.REPS:
                    self.assertAlmostEqual(num(ko[pg][f"centred_{rp}"]), centred[(rp, pg)], places=9)
            self.assertEqual(int(ko["P41"]["n_precursors_Rep1"]), 3)
            self.assertEqual(prov["q_columns"], ["Q.Value", "PG.Q.Value", "Channel.Q.Value"])
            self.assertEqual(prov["settings"]["run_level_only"]["given"], True)
            self.assertIn("--run-level-only", prov["replay"]["args"])
            self.assertNotIn("Global", prov["methods_paragraph"])
            with open(os.path.join(out, "reproducibility_log.R")) as fh:
                self.assertIn("cutoffs <- c('Q.Value' = 0.01, 'PG.Q.Value' = 0.05)", fh.read())

    def test_reinjections_are_one_replicate(self):
        # run_A2 re-injects run_A1's mix: they share Rep1, are averaged, and count once
        with tempfile.TemporaryDirectory() as d:
            report = fx.make_report(os.path.join(d, "report.parquet"))
            design = fx.write_design(os.path.join(d, "design.csv"))
            text = open(design).read()
            with open(design, "w") as fh:
                fh.write(text.replace("run_A2,Rep2", "run_A2,Rep1"))
            _, out = run_silac(d, "--replicates", "biological")
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            self.assertEqual(prov["design_record"]["reinjected_replicates"], {"Rep1": ["run_A1", "run_A2"]})
            self.assertEqual(prov["n_replicates"], 3)
            t = prov["de_tables"]["KO-WT"]
            self.assertEqual((t["n_replicates_available"], t["min_replicates"]), (3, 3))
            self.assertEqual(sorted(t["offset_log2"]), ["Rep1_run_A1", "Rep1_run_A2", "Rep3", "Rep4"])
            rep = {"run_A1": "Rep1", "run_A2": "Rep1", "run_B1": "Rep3", "run_B2": "Rep4"}
            _, _, centred = hand_ratios(report, "KO", "WT", rep)
            ko = {r["Protein.Group"]: r for r in read(os.path.join(out, "DE_silac_KO.WT.csv"))}
            self.assertNotIn("centred_Rep2", ko["P03"])
            for pg in ("P03", "P12", "P44"):
                for rp in ("Rep1", "Rep3", "Rep4"):
                    self.assertAlmostEqual(num(ko[pg][f"centred_{rp}"]), centred[(rp, pg)], places=9)
                # distinct precursors: a re-injection does not count one twice
                self.assertEqual(ko[pg]["n_precursors_Rep1"], ko[pg]["n_precursors_Rep3"])
            self.assertIn("averaged into one replicate (Rep1: 2 runs)", prov["methods_paragraph"])
            self.assertIn("Each replicate was a separate culture and mix.", prov["methods_paragraph"])
            self.assertEqual(prov["settings"]["replicates"]["source"],
                             "given (--replicates biological): each Replicate is a separate culture and mix")
            self.assertNotIn("silac_replicate_independence", [c["check"] for c in prov["checks"]])
            samples = [r["File.Name"] for r in read(os.path.join(out, "SILAC_samples.csv"))]
            self.assertIn("Rep1_run_A2_KO", samples)
            self.assertIn("Rep3_KO", samples)
            # different channel maps under one Replicate are two mixes, not a re-injection
            text = open(design).read()
            with open(design, "w") as fh:
                fh.write(text.replace("run_A2,Rep1,L,WT", "run_A2,Rep1,L,XX"))
            p, _ = run_silac(d, ok=False, outdir=os.path.join(d, "o2"))
            self.assertIn("different channel assignments", p.stderr)

    def test_default_minimum_grows_with_the_replicates(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet"), n_runs=8)
            fx.write_design(os.path.join(d, "design.csv"), n_runs=8)
            _, out = run_silac(d, "--comparisons", "KO-WT")
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            self.assertEqual(prov["de_tables"]["KO-WT"]["min_replicates"], 4, "max(3, ceiling(8 / 2))")
            self.assertIn("KO/WT: 4 of 8", prov["methods_paragraph"])

    def test_zero_point_check_reports_the_mode_centred_count(self):
        # one channel's losses fill the lower half: the median (the zero) sits 0.5 log2 below the
        # density peak -- silac_compare says so and recounts with the peak as the zero
        code = r"""
suppressMessages({library(data.table); library(limma)})
source(file.path(Sys.getenv("SCRIPTS"), "silac_ratios.R"))
set.seed(3)
pr <- c(rnorm(90, 0, 0.08), runif(60, -1.2, -0.4), rnorm(50, -2.5, 0.08))
w <- rbindlist(lapply(1:4, function(i) data.table(Run = paste0("r", i), Replicate = paste0("R", i),
  Protein.Group = sprintf("P%03d", 1:200), Precursor.Id = sprintf("pep%03d", 1:200),
  WT = 2^12, KO = 2^(12 - 1.5 + pr + rnorm(200, 0, 0.1)))))
res <- silac_compare(w, "KO", "WT", paste0("R", 1:4), 3)
cat(jsonlite::toJSON(list(mc = res$mode_check, off = res$offsets), auto_unbox = TRUE, digits = NA))
"""
        p = subprocess.run([RSCRIPT, "-e", code], capture_output=True, text=True, timeout=300,
                           env=dict(os.environ, SCRIPTS=SCRIPTS))
        self.assertEqual(p.returncode, 0, p.stderr)
        r = json.loads(p.stdout)
        self.assertGreater(r["mc"]["max_abs_difference"], 0.25)
        self.assertIsInstance(r["mc"]["n_significant_mode_centred"], int)
        for o in r["off"]:
            self.assertAlmostEqual(o["mode"], -1.5, delta=0.15)
            self.assertLess(o["offset"], o["mode"] - 0.25)

    def test_label_swap_compares_groups_not_channels(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet"), swap=True)
            fx.write_design(os.path.join(d, "design.csv"), swap=True)
            _, out = run_silac(d, "--comparisons", "KO-WT")
            ko = {r["Protein.Group"]: r for r in read(os.path.join(out, "DE_silac_KO.WT.csv"))}
            self.assertGreater(num(ko["P01"]["logFC"]), 1.5)
            self.assertLess(num(ko["P12"]["logFC"]), -1.5)
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            self.assertEqual(sorted(prov["design_record"]["label_swapped_groups"]), ["KO", "WT"])
            self.assertEqual(prov["comparisons_source"], "given (--comparisons)")

    def test_label_check_run_gives_incorporation(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet"), label_check_run="run_C1")
            fx.write_design(os.path.join(d, "design.csv"), extra=[
                {"File.Name": "run_C1", "Replicate": "Check1", "Channel": "H", "Group": "DKO"}])
            _, out = run_silac(d)
            inc = read(os.path.join(out, "SILAC_label_incorporation.csv"))
            self.assertEqual(len(inc), 1)
            self.assertAlmostEqual(num(inc[0]["median_percent"]), 100 / 1.02, delta=0.5)
            self.assertEqual(int(inc[0]["precursors_without_light_row"]), 1)
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            self.assertTrue(prov["label_incorporation"]["measured"])
            self.assertEqual(prov["n_replicates"], 4, "a label-check run is not a replicate")
            chk = next(c for c in prov["checks"] if c["check"] == "silac_label_incorporation")
            self.assertEqual(chk["status"], "PASS")

    def test_a_label_declared_twice_is_flagged(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet"), doubled=True)
            fx.write_design(os.path.join(d, "design.csv"))
            _, out = run_silac(d)
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            chk = [c for c in prov["checks"] if c["check"] == "silac_label_declaration"]
            self.assertEqual(chk and chk[0]["status"], "WARN")
            self.assertGreater(prov["label_declaration"]["doubled_fraction"], 0.9)

    def test_labels_and_channel_norm_from_the_params(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet")); fx.write_design(os.path.join(d, "design.csv"))
            cfg = os.path.join(d, "params.cfg")
            open(cfg, "w").write("--fixed-mod SILAC,0.0,KR,label\n")
            rec = silac.labels_record("silac3", "K4R6,K8R10")       # normalisation: the default
            json.dump({"labels": rec}, open(cfg + ".rationale.json", "w"))
            _, out = run_silac(d, "--params", cfg)
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            self.assertEqual(prov["labels"]["channels"][2]["preset"], "K8R10")
            self.assertTrue(any(c["check"] == "silac_channel_norm" for c in prov["checks"]))
            # pulsed SILAC / turnover: refused, nothing written
            rec = silac.labels_record("silac3", "K4R6,K8R10", "run")
            json.dump({"labels": rec}, open(cfg + ".rationale.json", "w"))
            p, out2 = run_silac(d, "--params", cfg, ok=False, outdir=os.path.join(d, "o2"))
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("turnover", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(out2, "de_provenance.json")))

    def test_another_pipelines_results_are_not_overwritten(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet")); fx.write_design(os.path.join(d, "design.csv"))
            out = os.path.join(d, "tables")
            os.makedirs(out)
            json.dump({"pipeline_id": "dpc"}, open(os.path.join(out, "de_provenance.json"), "w"))
            p, _ = run_silac(d, ok=False, outdir=out)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("holds another analysis (dpc)", p.stderr)
            self.assertEqual(json.load(open(os.path.join(out, "de_provenance.json"))),
                             {"pipeline_id": "dpc"})

    def test_a_run_left_out_is_the_users_recorded_decision(self):
        with tempfile.TemporaryDirectory() as d:
            report = fx.make_report(os.path.join(d, "report.parquet"), label_check_run="run_C1")
            design = fx.write_design(os.path.join(d, "design.csv"))
            p, _ = run_silac(d, ok=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("report run(s) with no row: run_C1", p.stderr)
            self.assertIn("--exclude <run> --reason", p.stderr)
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "silac.py"), "design",
                                "--validate", design, "--against", report, "--exclude", "run_C1",
                                "--reason", "failed injection"], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            _, out = run_silac(d)
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            self.assertEqual(prov["runs_left_out"]["runs"],
                             [{"run": "run_C1", "reason": "failed injection"}])
            self.assertEqual(prov["runs_left_out"]["n_analysed"], 4)
            m = os.path.join(d, "methods.md")
            subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"), "--raw",
                            os.path.join(d, "x.d"), "--instrument", "timsTOF HT", "--acquisition",
                            "DIA", "--de-dir", out, "--out", m], check=True, capture_output=True)
            self.assertIn("Of the 5 runs searched, 4 were analysed and 1 was left out",
                          open(m, encoding="utf-8").read())
            # an exclusion without a reason is refused
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "silac.py"), "design",
                                "--validate", design, "--against", report, "--exclude", "run_C1",
                                "--reason", " "], capture_output=True, text=True)
            self.assertNotEqual(r.returncode, 0)

    def test_a_bad_design_stops_with_the_reason(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet"))
            fx.write_design(os.path.join(d, "design.csv"))
            rows = open(os.path.join(d, "design.csv")).read().replace("run_B2", "run_X9")
            open(os.path.join(d, "design.csv"), "w").write(rows)
            p, _ = run_silac(d, ok=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("run_X9", p.stderr)


def ncbi_names(report):
    """The report as DIA-NN 2.7.0 writes it on an NCBI FASTA: Genes and Protein.Names both hold
    the description's first word for P01-P30; P31 on keeps its UniProt-style names."""
    import pyarrow.parquet as pq
    import pyarrow as pa
    t = pq.read_table(report).to_pydict()
    words = ("MULTISPECIES:", "hypothetical", "30S", "fimbrial")
    for i, pg in enumerate(t["Protein.Group"]):
        if pg.startswith("P") and int(pg[1:]) <= 30:
            t["Genes"][i] = t["Protein.Names"][i] = words[int(pg[1:]) % 4]
    pq.write_table(pa.table(t), report)


def write_ncbi_fasta(path):
    """P01-P20 with NCBI headers (P01 a CDS-style one with [gene=...]); P31 a UniProt one."""
    with open(path, "w") as fh:
        fh.write(">P01 [gene=sampA] [locus_tag=X_0001] [protein=outer membrane porin A]\nMPEPTIDEK\n")
        for i in range(2, 21):
            fh.write(f">P{i:02d} MULTISPECIES: hypothetical protein {i} [Strain-X]\nMPEPTIDEK\n")
        fh.write(">sp|P31|P31_TEST Test protein 31 OS=Strain-X OX=9 GN=G31 PE=1 SV=1\nMPEPTIDEK\n")
    return path


@unittest.skipUnless(HAVE_R, NEEDS)
class NonUniprotNames(unittest.TestCase):
    """On an NCBI database DIA-NN's Genes and Protein.Names hold the description's first word
    ("MULTISPECIES:", "hypothetical"): a client read them as gene names in the tables and in
    the report's Top proteins. They are blanked, a [gene=...] token kept, Description used."""

    def check(self, out, with_fasta):
        ko = {r["Protein.Group"]: r for r in read(os.path.join(out, "DE_silac_KO.WT.csv"))}
        for pg in ("P02", "P12", "P21", "P30"):
            self.assertEqual((ko[pg]["Genes"], ko[pg]["Protein.Names"]), ("", ""), pg)
        self.assertEqual(ko["P31"]["Genes"], "G31")
        self.assertEqual(ko["P31"]["Protein.Names"], "P31_TEST")
        em = {r["Protein.Group"]: r for r in read(os.path.join(out, "SILAC_sample_matrix.csv"))}
        self.assertEqual(em["P02"]["Genes"], "")
        prov = json.load(open(os.path.join(out, "de_provenance.json")))
        self.assertEqual(prov["protein_names"]["n_not_uniprot"], 30)
        self.assertTrue(any("not UniProt's" in n for n in prov["brief_notes"]))
        self.assertIn("Genes", prov["de_columns"])
        if with_fasta:
            self.assertEqual(ko["P01"]["Genes"], "sampA", "the [gene=...] token is the gene")
            self.assertEqual(ko["P01"]["Description"], "outer membrane porin A")
            self.assertEqual(ko["P02"]["Description"], "hypothetical protein 2")
            self.assertEqual(prov["protein_names"]["n_gene_from_token"], 1)
        else:
            self.assertEqual(ko["P01"]["Genes"], "")
        return ko

    def test_names_from_an_ncbi_database_are_not_shown_as_genes(self):
        for with_fasta in (True, False):
            with self.subTest(with_fasta=with_fasta), tempfile.TemporaryDirectory() as d:
                fx.make_report(os.path.join(d, "report.parquet"))
                ncbi_names(os.path.join(d, "report.parquet"))
                fx.write_design(os.path.join(d, "design.csv"))
                more = ("--fasta", write_ncbi_fasta(os.path.join(d, "db.fasta"))) if with_fasta else ()
                _, out = run_silac(d, *more)
                self.check(out, with_fasta)
                if not with_fasta:
                    continue
                # the report's Top proteins: a Description column, no first-word "genes"
                html = os.path.join(d, "Analysis_Report.html")
                p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                                    "--tables", out, "--figures", os.path.join(d, "out", "figures"),
                                    "--out", html, "--no-pdf"], capture_output=True, text=True,
                                   timeout=300)
                self.assertEqual(p.returncode, 0, p.stdout[-2000:] + p.stderr[-2000:])
                md = open(os.path.join(os.path.dirname(html), "Analysis_Report.md")).read()
                top = md[md.index("## Top proteins per contrast"):]
                self.assertIn("| Protein | Gene | Description | log2FC |", top)
                self.assertIn("| P01 | sampA | outer membrane porin A |", top)
                self.assertRegex(top, r"\| P02 \| \u2014 \| hypothetical protein 2 \|")
                for w in ("MULTISPECIES:", "| hypothetical |", "| 30S |", "| fimbrial |"):
                    self.assertNotIn(w, top)


@unittest.skipUnless(r_has("limma", "jsonlite"), "needs Rscript with limma and jsonlite")
class RunDeLeavesSilacResultsAlone(unittest.TestCase):
    """run_de.R refuses an --outdir holding SILAC results (as run_silac.R refuses run_de.R's):
    the audit, report and delivery read every DE_*.csv there, and DE_silac_*.csv must never pass
    as a label-free result."""

    def run_de(self, outdir):
        d = os.path.dirname(outdir)
        for f in ("report.parquet", "conditions.csv"):
            open(os.path.join(d, f), "w").write("File.Name,Group\n")
        return subprocess.run([RSCRIPT, os.path.join(SCRIPTS, "run_de.R"), "--input",
                               os.path.join(d, "report.parquet"), "--metadata",
                               os.path.join(d, "conditions.csv"), "--outdir", outdir],
                              capture_output=True, text=True, timeout=300)

    def test_a_silac_record_or_silac_tables_stop_it_before_anything_is_written(self):
        for case in ("record", "tables"):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as d:
                out = os.path.join(d, "tables")
                os.makedirs(out)
                if case == "record":
                    rec = {"pipeline_id": "silac_ratio", "replay": {"script": "run_silac.R"}}
                    json.dump(rec, open(os.path.join(out, "de_provenance.json"), "w"))
                else:
                    open(os.path.join(out, "DE_silac_KO.WT.csv"), "w").write("Protein.Group\n")
                before = sorted(os.listdir(out))
                p = self.run_de(out)
                self.assertNotEqual(p.returncode, 0)
                self.assertIn("holds another analysis", p.stderr)
                self.assertIn("silac", p.stderr.lower())
                self.assertEqual(sorted(os.listdir(out)), before, "nothing was written")


class DesignRules(unittest.TestCase):
    def write(self, d, text):
        p = os.path.join(d, "design.csv")
        open(p, "w").write(text)
        return p

    def test_default_comparisons_heavier_over_lighter(self):
        with tempfile.TemporaryDirectory() as d:
            des = silac.read_design(fx.write_design(os.path.join(d, "x.csv")), ["L", "M", "H"])
            self.assertEqual([c["contrast"] for c in des["comparisons"]],
                             ["KO-WT", "DKO-WT", "DKO-KO"])
            self.assertTrue(des["comparisons_source"].startswith("DEFAULT"))
            self.assertEqual({a["role"] for a in des["assignments"]}, {"ratio"})

    def test_rules(self):
        head = "File.Name,Replicate,Channel,Group\n"
        bad = {
            "missing column": "File.Name,Channel\nr1,L\n",
            "unknown channel": head + "r1,R1,X,WT\nr1,R1,H,KO\n",
            "group twice in a run": head + "r1,R1,L,WT\nr1,R1,H,WT\n",
            "channel twice": head + "r1,R1,L,WT\nr1,R1,L,KO\n",
            "one replicate, two mixes": head + "r1,R1,L,WT\nr1,R1,H,KO\nr2,R1,L,KO\nr2,R1,H,WT\n",
            "bad group name": head + "r1,R1,L,W-T\nr1,R1,H,KO\n",
            "no run with two groups": head + "r1,R1,H,KO\n",
        }
        for why, text in bad.items():
            with self.subTest(why=why), tempfile.TemporaryDirectory() as d:
                with self.assertRaises(silac.SilacError):
                    silac.read_design(self.write(d, text), ["L", "M", "H"])
        with tempfile.TemporaryDirectory() as d:
            p = self.write(d, head + "r1,R1,L,WT\nr1,R1,H,KO\nr2,R2,L,WT\nr2,R2,H,KO\n")
            with self.assertRaises(silac.SilacError) as e:
                silac.read_design(p, ["L", "H"], comparisons="KO-XX")
            self.assertIn("KO-XX", str(e.exception))
            with self.assertRaises(silac.SilacError) as e:
                silac.read_design(p, ["L", "H"], report_runs=["r1", "r2", "r3"])
            self.assertIn("r3", str(e.exception))
            des = silac.read_design(p, ["L", "H"], comparisons="WT-KO")
            self.assertEqual(des["comparisons"][0]["contrast"], "WT-KO")

    def test_reinjected_replicate(self):
        with tempfile.TemporaryDirectory() as d:
            p = self.write(d, "File.Name,Replicate,Channel,Group\nr1,R1,L,WT\nr1,R1,H,KO\n"
                              "r1b,R1,L,WT\nr1b,R1,H,KO\nr2,R2,L,WT\nr2,R2,H,KO\n")
            des = silac.read_design(p, ["L", "H"])
            self.assertEqual(des["reinjected_replicates"], {"R1": ["r1", "r1b"]})
            self.assertEqual(len(des["comparisons"]), 1)

    def test_label_check_and_unused_runs(self):
        with tempfile.TemporaryDirectory() as d:
            p = self.write(d, "File.Name,Replicate,Channel,Group\nr1,R1,L,WT\nr1,R1,H,KO\n"
                              "chk,C,H,KO\nlight,Lt,L,WT\n")
            des = silac.read_design(p, ["L", "H"])
            self.assertEqual(des["label_check_runs"], ["chk"])
            self.assertEqual(des["unused_runs"], ["light"])

    def test_cli_json(self):
        with tempfile.TemporaryDirectory() as d:
            p = fx.write_design(os.path.join(d, "x.csv"))
            runs = os.path.join(d, "runs.txt")
            open(runs, "w").write("\n".join(fx.RUNS) + "\n")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "silac.py"), "design",
                                "--validate", p, "--runs-file", runs, "--json"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertTrue(json.loads(r.stdout)["ok"])
            open(runs, "a").write("run_extra\n")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "silac.py"), "design",
                                "--validate", p, "--runs-file", runs, "--json"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 1)
            self.assertIn("run_extra", json.loads(r.stdout.splitlines()[0])["error"])


if __name__ == "__main__":
    unittest.main()
