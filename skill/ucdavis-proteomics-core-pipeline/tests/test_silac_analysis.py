#!/usr/bin/env python3
"""
The SILAC ratio analysis (run_silac.R + silac_ratios.R) on a synthetic plexDIA report
(tests/silac_fixture.py: planted effects, a mixing offset, channel-specific proteins).

The method was validated on a real triple-SILAC experiment (4 replicates, DIA-NN 2.7.0): it
reproduced the hand-made analysis's tables to floating-point precision (1,224 / 1,153 / 1,093
proteins tested; 816 / 799 / 0 significant). These tests pin, on data whose answer is known:
  * the planted changes are found, in the right direction, and nothing else;
  * the ratio maths, recomputed here in Python: the precursor median, the per-replicate
    centring, the limma mean, Benjamini-Hochberg, the minimum number of replicates;
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


HAVE_R = r_has("limma", "data.table", "arrow", "jsonlite")


def read(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def num(x):
    return float("nan") if x in ("", "NA", None) else float(x)


def run_silac(d, *more, report=None, design=None, ok=True, outdir=None):
    outdir = outdir or os.path.join(d, "out", "tables")
    p = subprocess.run([RSCRIPT, os.path.join(SCRIPTS, "run_silac.R"),
                        "--input", report or os.path.join(d, "report.parquet"),
                        "--design", design or os.path.join(d, "design.csv"),
                        "--outdir", outdir, *more], capture_output=True, text=True, timeout=600)
    if ok:
        assert p.returncode == 0, p.stdout + p.stderr
    return p, outdir


@unittest.skipUnless(HAVE_R, "needs Rscript with limma, data.table, arrow and jsonlite")
class Analysis(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.d = cls.tmp.name
        fx.make_report(os.path.join(cls.d, "report.parquet"))
        fx.write_design(os.path.join(cls.d, "design.csv"))
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
        for k in ("P51", "P55", "Cont_P00761"):
            self.assertNotIn(k, self.ko)
        spec = {r["Protein.Group"]: r for r in read(os.path.join(self.out,
                                                                 "SILAC_channel_specific_detections.csv"))}
        self.assertEqual(sorted(spec), [f"P{i}" for i in range(51, 56)])
        self.assertEqual(spec["P51"]["pattern"], "in WT; not in KO, DKO")
        self.assertEqual(self.prov["filter_counts"]["contaminant_protein_groups"], 1)
        self.assertEqual(self.ko["P01"]["Detected_KO"], "4/4")

    def test_ratio_maths_by_hand(self):
        import pyarrow.parquet as pq
        t = pq.read_table(os.path.join(self.d, "report.parquet")).to_pylist()
        rep = dict(zip(fx.RUNS, fx.REPS))
        grp = {"L": "WT", "M": "KO", "H": "DKO"}
        q = {}
        for r in t:
            if (r["Channel"] in grp and r["Channel.Q.Value"] <= 0.01 and r["Q.Value"] <= 0.01
                    and r["PG.Q.Value"] <= 0.05 and not r["Protein.Ids"].startswith("Cont_")):
                q.setdefault((rep[r["Run"]], r["Protein.Group"], r["Precursor.Id"]), {})[
                    grp[r["Channel"]]] = r["Precursor.Quantity"]
        prot = {}
        for (rp, pg, _), v in q.items():
            if "KO" in v and "WT" in v:
                prot.setdefault((rp, pg), []).append(math.log2(v["KO"] / v["WT"]))
        med = {k: statistics.median(v) for k, v in prot.items()}
        offset = {rp: statistics.median([m for (r, _), m in med.items() if r == rp]) for rp in fx.REPS}
        for rp in fx.REPS:
            self.assertAlmostEqual(offset[rp], self.prov["de_tables"]["KO-WT"]["offset_log2"][rp], places=9)
        for pg in ("P03", "P12", "P30", "P44"):
            row = self.ko[pg]
            for rp in fx.REPS:
                self.assertAlmostEqual(num(row[f"uncentred_{rp}"]), med[(rp, pg)], places=9)
                self.assertAlmostEqual(num(row[f"centred_{rp}"]), med[(rp, pg)] - offset[rp], places=9)
                self.assertEqual(int(row[f"n_precursors_{rp}"]), len(prot[(rp, pg)]))
            # intercept-only limma: the estimate is the mean of the centred values
            self.assertAlmostEqual(num(row["logFC"]),
                                   statistics.mean(num(row[f"centred_{rp}"]) for rp in fx.REPS), places=9)
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
        for want in ("median-centred", "at least 3 of the 4 replicates", "robust = TRUE",
                     "Benjamini-Hochberg", "Derks", "Phipson", "no fold-change filter"):
            self.assertIn(want, p["methods_paragraph"])
        self.assertEqual(p["settings"]["min_replicates"]["value"], 3)
        self.assertIn("default", p["settings"]["min_replicates"]["source"])
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
                     "qc_pvalue_panel.png", "silac_scatter_KO.WT_vs_DKO.KO.png"):
            self.assertIn(want, files)
        for f in fj["figures"]:
            self.assertTrue(os.path.getsize(os.path.join(figdir, f["file"])) > 0, f["file"])
            self.assertTrue(f["caption"], f["file"])
        self.assertEqual(fj["failed"], [])
        scatter = next(f for f in fj["figures"] if f["file"].startswith("silac_scatter"))
        self.assertIn("share the KO channel", scatter["caption"])

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


@unittest.skipUnless(HAVE_R, "needs Rscript with limma, data.table, arrow and jsonlite")
class Variants(unittest.TestCase):
    def test_min_replicates_is_a_recorded_parameter(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet")); fx.write_design(os.path.join(d, "design.csv"))
            _, out = run_silac(d, "--min-replicates", "4")
            ko = {r["Protein.Group"]: r for r in read(os.path.join(out, "DE_silac_KO.WT.csv"))}
            self.assertIn("ratio in 2 of 4", ko["P56"]["Evidence"])
            prov = json.load(open(os.path.join(out, "de_provenance.json")))
            self.assertEqual(prov["settings"]["min_replicates"],
                             {"value": 4, "source": "given (--min-replicates)"})
            self.assertIn("at least 4 of the 4 replicates", prov["methods_paragraph"])
            self.assertTrue(any(c["check"] == "silac_replicates" and c["status"] == "WARN"
                                for c in prov["checks"]))

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

    def test_a_bad_design_stops_with_the_reason(self):
        with tempfile.TemporaryDirectory() as d:
            fx.make_report(os.path.join(d, "report.parquet"))
            fx.write_design(os.path.join(d, "design.csv"))
            rows = open(os.path.join(d, "design.csv")).read().replace("run_B2", "run_X9")
            open(os.path.join(d, "design.csv"), "w").write(rows)
            p, _ = run_silac(d, ok=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("run_X9", p.stderr)


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
            "replicate names two runs": head + "r1,R1,L,WT\nr1,R1,H,KO\nr2,R1,L,WT\nr2,R1,H,KO\n",
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
