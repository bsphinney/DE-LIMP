#!/usr/bin/env python3
"""
2.15.1 audit (2026-11-17), the MaxLFQ route on DDA engines' quantities (build_maxlfq.R, run_de.R).

L2-3  (lane3 #4, lane4 C2) every declared quantity (Sage, FragPipe, AlphaDIA) got a quantile step
      on top. Quantile forces every run to one distribution, so a lopsided experiment (a quarter of
      the proteome up) is pulled back and the unchanged proteins come out "down". Now --between-run
      auto reads the declaration's engine_normalised: "true" (IonQuant / directLFQ normalised)
      -> no second step; "false" (Sage LFQ) -> median_ratio, a location-only shift fitted on the
      proteins quantified in every run; not recorded -> median_ratio, said "[confirm]".
L2-2b a non-normalised (raw) DE of an adapted report: allowed when the engine did not normalise
      (Sage), refused with the way out when it did (FragPipe ionquant.normalization=1).
L2-4  a Sage protein = Tukey median polish of its peptides (one interfered peptide does not move
      it), not the single most intense peptide.
L2-6  Identified.In.Run: a value measured by matching between runs (MS1 signal, no MS2
      identification in that run) is quantified, tested and counts as present (one definition,
      2.15.1 review finding 3); the DE tables flag groups whose values are all matched
      (matched_between_runs_only_in).
#16   an undeclared report that is not DIA-NN-shaped is not described as "DIA-NN PG.MaxLFQ".
L2-10 / lane3 #9: contrast directions printed and recorded; default thresholds tagged.

Synthetic reports only. Needs pyarrow and R with limma/arrow/dplyr/tidyr/jsonlite.
"""
import csv
import json
import os
import random
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

from test_run_de_contaminants import r_has  # noqa: E402

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

NEEDS = ("limma", "arrow", "dplyr", "tidyr", "jsonlite")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")
RUNS = [(f"c{i}", "Control") for i in range(1, 5)] + [(f"t{i}", "Treated") for i in range(1, 5)]
Q_COLS = ("Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value")


def declared(label="Sage", engine_normalised=None):
    m = {"delimp.quantity.level": "peptide" if label == "Sage" else "protein",
         "delimp.quantity.label": label,
         "delimp.quantity.value": f"{label} synthetic intensity",
         "delimp.quantity.identification": "synthetic",
         "delimp.quantity.q_columns": "placeholder",
         "delimp.quantity.citation": "synthetic"}
    if engine_normalised is not None:
        m["delimp.quantity.engine_normalised"] = engine_normalised
    return m


def write_report(out, meta=None, n_prot=400, frac_up=0.25, up=1.5, n_pep=3, seed=5,
                 interfere=None, identified=None, peptide_col=True):
    """A Sage-shaped (peptide x run) report: each run a loading offset, `frac_up` of the proteins
    `up` log2 higher in Treated, no missing values. `interfere`: {protein: (peptide index, run,
    +log2)} spikes one peptide in one run. `identified`: {(protein, run): 0} marks every peptide
    of that protein in that run not identified by MS2 there (Identified.In.Run 0: measured by
    matching between runs), others 1."""
    rnd = random.Random(seed)
    os.makedirs(out, exist_ok=True)
    load = {r: rnd.uniform(-0.6, 0.6) for r, _ in RUNS}
    rows = {k: [] for k in ("Run", "Protein.Group", "Protein.Names", "Peptide.Id", "PG.MaxLFQ")
            + Q_COLS}
    if identified is not None:
        rows["Identified.In.Run"] = []
    n_up = int(n_prot * frac_up)
    for i in range(n_prot):
        acc, base = f"P{i:05d}", rnd.uniform(16, 24)
        for k in range(n_pep):
            off = rnd.gauss(0, 1)
            for r, g in RUNS:
                lv = base + off + load[r] + (up if g == "Treated" and i < n_up else 0.0) \
                    + rnd.gauss(0, 0.15)
                if interfere and acc in interfere and interfere[acc][:2] == (k, r):
                    lv += interfere[acc][2]
                for c, x in (("Run", r), ("Protein.Group", acc), ("Protein.Names", f"N{i}"),
                             ("Peptide.Id", f"PEP{i}K{k}"), ("PG.MaxLFQ", 2 ** lv)):
                    rows[c].append(x)
                for c in Q_COLS:
                    rows[c].append(0.0)
                if identified is not None:
                    rows["Identified.In.Run"].append(identified.get((acc, r), 1))
    if not peptide_col:
        del rows["Peptide.Id"]
    t = pa.table(rows)
    if meta:
        t = t.replace_schema_metadata(meta)
    report = os.path.join(out, "report.parquet")
    pq.write_table(t, report)
    cond = os.path.join(out, "conditions.csv")
    with open(cond, "w") as fh:
        fh.write("File.Name,Group\n" + "".join(f"{r},{g}\n" for r, g in RUNS))
    return report, cond, n_up


R_CALL = r'''
setwd("{scripts}")
source("build_maxlfq.R")
m <- build_maxlfq("{report}", quantities = "{quantities}", between_run = "{between_run}")
write.csv(m$E, file.path("{out}", "E.csv"))
jsonlite::write_json(m$descriptor, file.path("{out}", "descriptor.json"), auto_unbox = TRUE,
                     null = "null", na = "null", digits = NA)
'''


def run_maxlfq(report, out, between_run="auto", quantities="normalised"):
    os.makedirs(out, exist_ok=True)
    p = subprocess.run(["Rscript", "-e", R_CALL.format(scripts=SCRIPTS, report=report, out=out,
                                                       quantities=quantities,
                                                       between_run=between_run)],
                       capture_output=True, text=True, timeout=300)
    if p.returncode != 0:
        return None, p.stderr
    with open(os.path.join(out, "E.csv"), newline="") as fh:
        r = csv.reader(fh)
        runs = next(r)[1:]
        E = {row[0]: {run: float(v) if v not in ("NA", "") else None
                      for run, v in zip(runs, row[1:])} for row in r}
    with open(os.path.join(out, "descriptor.json")) as fh:
        return E, json.load(fh)


def null_logfc(E, n_up):
    """Median Treated - Control log2 difference over the unchanged proteins."""
    d = []
    for p, row in E.items():
        if int(p[1:]) < n_up:
            continue
        t = [v for r, v in row.items() if r.startswith("t")]
        c = [v for r, v in row.items() if r.startswith("c")]
        d.append(sum(t) / len(t) - sum(c) / len(c))
    d.sort()
    return d[len(d) // 2]


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class BetweenRunByDeclaration(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_sage_lopsided_median_ratio_keeps_the_null_at_zero(self):
        report, _, n_up = write_report(os.path.join(self.d, "s"),
                                       declared("Sage", engine_normalised="false"))
        E, desc = run_maxlfq(report, os.path.join(self.d, "auto"))
        self.assertEqual(desc["between_run"]["method"], "median_ratio")
        self.assertNotIn("[confirm]", desc["normalisation"])
        self.assertLess(abs(null_logfc(E, n_up)), 0.15)
        # the step 2.15.0 applied pulls the unchanged quarter of ... the rest of the proteome down
        Eq, dq = run_maxlfq(report, os.path.join(self.d, "q"), between_run="quantile")
        self.assertEqual(dq["between_run"]["method"], "quantile")
        self.assertLess(null_logfc(Eq, n_up), -0.2)

    def test_engine_normalised_true_gets_no_second_step(self):
        report, _, _ = write_report(os.path.join(self.d, "f"),
                                    declared("FragPipe", engine_normalised="true"))
        _, desc = run_maxlfq(report, os.path.join(self.d, "o"))
        self.assertEqual(desc["between_run"]["method"], "none")
        self.assertIs(desc["between_run"]["engine_normalised"]["value"], True)

    def test_not_recorded_is_median_ratio_said_confirm(self):
        report, _, _ = write_report(os.path.join(self.d, "u"), declared("FragPipe"))
        _, desc = run_maxlfq(report, os.path.join(self.d, "o"))
        self.assertEqual(desc["between_run"]["method"], "median_ratio")
        self.assertIn("[confirm]", desc["normalisation"])

    def test_raw_sage_is_allowed_and_raw_normalised_fragpipe_refused(self):
        report, _, _ = write_report(os.path.join(self.d, "s"),
                                    declared("Sage", engine_normalised="false"))
        E, desc = run_maxlfq(report, os.path.join(self.d, "raw"), quantities="raw")
        self.assertIsNotNone(E, desc)
        self.assertEqual(desc["between_run"]["method"], "none")
        report, _, _ = write_report(os.path.join(self.d, "f"),
                                    declared("FragPipe", engine_normalised="true"))
        E, err = run_maxlfq(report, os.path.join(self.d, "rawf"), quantities="raw")
        self.assertIsNone(E)
        self.assertIn("ionquant.normalization", " ".join(err.split()))


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class Rollup(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_one_interfered_peptide_does_not_move_the_protein(self):
        meta = declared("Sage", engine_normalised="false")
        base, _, _ = write_report(os.path.join(self.d, "a"), meta, n_prot=60, n_pep=5)
        spiked, _, _ = write_report(os.path.join(self.d, "b"), meta, n_prot=60, n_pep=5,
                                    interfere={"P00059": (0, "t1", 6.0)})
        E0, d0 = run_maxlfq(base, os.path.join(self.d, "oa"), between_run="none")
        E1, _ = run_maxlfq(spiked, os.path.join(self.d, "ob"), between_run="none")
        self.assertEqual(d0["pipeline_id"], "peptide_medpolish")
        self.assertIn("median polish", d0["rollup_method"])
        self.assertLess(abs(E1["P00059"]["t1"] - E0["P00059"]["t1"]), 0.5)

    def test_undeclared_non_diann_report_is_not_called_diann(self):
        report, _, _ = write_report(os.path.join(self.d, "u"), None, n_prot=40, peptide_col=False)
        _, desc = run_maxlfq(report, os.path.join(self.d, "o"))
        self.assertEqual(desc["pipeline_id"], "undeclared_pg_maxlfq")
        self.assertIn("NOT RECORDED", desc["rollup_method"] + desc.get("caveat", ""))


@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class RunDeRecords(unittest.TestCase):
    """One run_de.R of a Sage-shaped report: matched values, directions, default tags, coverage."""

    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        # P00000: identified by MS2 in Treated only; its Control values are all matched between runs
        ident = {("P00000", r): 0 for r, g in RUNS if g == "Control"}
        cls.report, cls.cond, _ = write_report(cls.d, declared("Sage", engine_normalised="false"),
                                               n_prot=120, identified=ident)
        cls.out = os.path.join(cls.d, "de")
        cls.p = subprocess.run(["Rscript", RUN_DE, "--input", cls.report, "--metadata", cls.cond,
                                "--method", "maxlfq", "--outdir", cls.out],
                               capture_output=True, text=True, cwd=cls.d, timeout=600)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.d, ignore_errors=True)

    def setUp(self):
        self.assertEqual(self.p.returncode, 0, self.p.stderr[-3000:])

    def prov(self):
        with open(os.path.join(self.out, "de_provenance.json")) as fh:
            return json.load(fh)

    def test_matched_values_are_present_tested_and_flagged(self):
        """2.15.1 review finding 3: one definition of "present in a run" -- a value in the matrix
        limma tests, identified in its run or matched between runs. P00000 has values everywhere
        (its Control values are all matched): it is tested, flagged matched_between_runs_only_in
        = Control, and NOT on the on/off list. No protein is on/off AND has a p-value for the
        same contrast."""
        with open(os.path.join(self.out, "Detection_Matrix.csv"), newline="") as fh:
            rows = {r[next(iter(r))]: r for r in csv.DictReader(fh)}
        row = rows["P00000"]
        self.assertEqual({row[r] for r, _ in RUNS}, {"1"})
        de = os.path.join(self.out, "DE_maxlfq_Treated.Control.csv")
        with open(de, newline="") as fh:
            tab = {r["Protein.Group"]: r for r in csv.DictReader(fh)}
        self.assertNotIn(tab["P00000"]["P.Value"], ("", "NA"))
        self.assertEqual(tab["P00000"]["matched_between_runs_only_in"], "Control")
        self.assertEqual(tab["P00001"]["matched_between_runs_only_in"], "")
        self.assertNotIn("quantified_by_transfer_only_in", tab["P00000"])
        with open(os.path.join(self.out, "OnOff_Proteins.csv"), newline="") as fh:
            onoff = list(csv.DictReader(fh))
        self.assertNotIn("P00000", {r["Protein.Group"] for r in onoff})
        for r in onoff:
            if r["Contrast"] == "Treated - Control" and r["Protein.Group"] in tab:
                self.assertIn(tab[r["Protein.Group"]]["P.Value"], ("", "NA"), r)
        dm = self.prov()["detection_matrix"]
        self.assertEqual(dm["zero_means"], "missing")
        self.assertEqual(dm["matched_cells"], 4)
        self.assertIn("matched_between_runs_only_in", dm["matched_rule"])
        self.assertNotIn("transferred_cells", dm)
        # this report declares no definition: run_de.R says so rather than invent one
        self.assertIn("the report does not say how such a value was measured",
                      dm["matched_definition"])
        for text in (dm["matched_rule"], dm["values"]):
            self.assertNotIn("transferred", text)

    def test_contrast_direction_printed_and_recorded(self):
        self.assertIn("positive logFC = higher in Treated than in Control", self.p.stderr
                      + self.p.stdout)
        prov = self.prov()
        self.assertEqual(prov["contrast_directions"]["Treated-Control"],
                         "positive logFC = higher in Treated than in Control")
        self.assertIn("first group alphabetically (Control)", prov["contrasts_source"])
        with open(os.path.join(self.out, "methods.txt")) as fh:
            txt = fh.read()
        self.assertIn("Treated-Control: positive logFC = higher in Treated than in Control", txt)

    def test_default_thresholds_are_tagged(self):
        prov = self.prov()
        self.assertEqual(prov["adjp_source"], "DEFAULT — not user-confirmed")
        with open(os.path.join(self.out, "methods.txt")) as fh:
            txt = fh.read()
        self.assertIn("adj.P.Val < 0.05 (DEFAULT — not user-confirmed)", txt)

    def test_the_methods_paragraph_tags_the_defaults(self):
        import make_methods as mm
        para = mm.de_paragraph(self.prov())
        self.assertIn("adj.P.Val < 0.05 (DEFAULT — not user-confirmed)", para)
        self.assertIn("The contrasts were the analysis default", para)
        self.assertIn("Treated-Control, positive logFC = higher in Treated than in Control", para)
        user = dict(self.prov(), adjp_source="--adjp", contrasts_source="--contrasts")
        para = mm.de_paragraph(user)
        self.assertNotIn("DEFAULT", para.split("called significant")[0].split("adj.P.Val")[-1])
        self.assertNotIn("The contrasts were the analysis default", para)
        q = dict(self.prov(), q_columns_role="real", identification_fdr=None,
                 q_cutoff_source="DEFAULT — not user-confirmed")
        self.assertIn("(DEFAULT — not user-confirmed)", mm.identification_filter(q)[1])

    def test_coverage_record(self):
        cov = self.prov()["coverage"]
        self.assertEqual(cov["min_all"], 4)
        self.assertEqual(cov["min_group"], {"Control": 2, "Treated": 2})
        self.assertIn("of the runs of any one group", cov["text"])


@unittest.skipUnless(r_has("jsonlite"), "needs Rscript")
class CoverageRule(unittest.TestCase):
    def test_group_specific_protein_is_kept(self):
        r = subprocess.run(["Rscript", "-e", f'''setwd("{SCRIPTS}"); source("build_maxlfq.R")
          E <- rbind(all = c(1,1,1,1,1,1,1,1),
                     tonly = c(NA,NA,NA,NA,NA,1,1,1),   # 3 of 8 overall, 3 of 5 Treated
                     one = c(1,NA,NA,NA,NA,NA,NA,NA))
          cr <- coverage_rule(E, c(rep("C",3), rep("T",5)), 0.5)
          cat(cr$keep, cr$min_all, cr$min_group[["C"]], cr$min_group[["T"]],
              cr$n_kept_by_group_only, sep = ",")'''],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr[-1500:])
        self.assertEqual(r.stdout.strip(), "TRUE,TRUE,FALSE,4,2,3,1")


if __name__ == "__main__":
    unittest.main()
