"""make_figures.R: violin_top_<contrast>.png -- the top proteins of every contrast, with each
point marked measured or not.

The mark is the point of the figure. DPC-Quant gives every protein a value in every run, so a
large fold change can rest on values that were never measured; a group with no measured run
makes the fold change a detection event rather than a magnitude. The mark must come from
Detection_Matrix.csv and nowhere else: with no detection record the figure says so and draws
no status legend (never fabricate visualisation data), and under MaxLFQ a zero means "Missing"
(no value exists), not "Inferred".

Needs Rscript with ggplot2; skipped otherwise (CI has no R).
"""
import csv
import json
import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "make_figures.R")

# A dot and an underscore in one group name: make.names() turns the contrast
# "Treat.high_dose-Ctrl" into the file name "Treat.high_dose.Ctrl", which cannot be split
# on "." -- the groups must come from de_provenance.json or a pair search.
GROUPS = ["Ctrl", "Other", "Treat.high_dose"]
CONTRASTS = ["Treat.high_dose-Ctrl", "Other-Ctrl"]
N_PROT = 30
EVENT = 1          # protein P00001: never measured in Ctrl


def make_name(s):
    """R's make.names() for the characters used here."""
    return "".join(c if c.isalnum() or c in "._" else "." for c in s)


def have_r_ggplot():
    if not shutil.which("Rscript"):
        return False
    r = subprocess.run(["Rscript", "-e", "quit(status = !requireNamespace('ggplot2', quietly=TRUE))"],
                       capture_output=True)
    return r.returncode == 0


def write_inputs(d, pipeline="dpc", detection=True, provenance=True):
    runs = [(f"20260924_Exploris480_DIA_60SPD_{g.replace('.', '')}_S{i:02d}_rep{i}", g)
            for g in GROUPS for i in range(1, 4)]
    with open(os.path.join(d, "conditions.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        w.writerows(runs)
    tables = os.path.join(d, "tables")
    os.makedirs(tables)

    def value(p, s, g):
        base = 18 + (p % 7) - (p % 3) * 0.5
        shift = {"Treat.high_dose": 3.0 - 0.1 * p, "Other": 0.4, "Ctrl": 0.0}[g] if p < 12 else 0.0
        return round(base + shift + ((p * 7 + s * 3) % 5) / 10, 3)

    def nobs(p, g, s):
        if p == EVENT and g == "Ctrl":
            return 0
        if p in (3, 5) and s % 3 == 0:            # a few scattered unmeasured runs
            return 0
        return 1 + (p + s) % 4

    expr, det = [], []
    for p in range(N_PROT):
        er, dr = [], []
        for s, (run, g) in enumerate(runs):
            n = nobs(p, g, s)
            v = value(p, s, g)
            if pipeline == "maxlfq" and n == 0:
                v = "NA"                           # MaxLFQ leaves a hole; DPC fills it
            er.append(v)
            dr.append(n if pipeline == "dpc" else int(n > 0))
        expr.append([f"P{p:05d}", f"GENE{p}", f"PROT{p}_HUMAN"] + er)
        det.append([f"P{p:05d}"] + dr)
    with open(os.path.join(tables, "Expression_Matrix.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Protein.Group", "Genes", "Protein.Names"] + [r for r, _ in runs])
        w.writerows(expr)
    if detection:
        with open(os.path.join(tables, "Detection_Matrix.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Protein.Group"] + [r for r, _ in runs])
            w.writerows(det)

    for c in CONTRASTS:
        with open(os.path.join(tables, f"DE_{pipeline}_{make_name(c)}.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Protein.Group", "Genes", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"])
            for p in range(N_PROT):
                lfc = (3.0 - 0.1 * p if c.startswith("Treat") else 0.4) if p < 12 else 0.05
                adj = 10 ** -(8 - p * 0.6) if p < 12 else 0.5 + p / 100
                w.writerow([f"P{p:05d}", f"GENE{p}", lfc, 19, 5, adj / 10, min(adj, 1.0), 2])

    if provenance:
        rollup = ("DPC-Quant (Detection Probability Curve quantification, dpcCN)" if pipeline == "dpc"
                  else "MaxLFQ (DIA-NN PG.MaxLFQ)")
        with open(os.path.join(tables, "de_provenance.json"), "w") as fh:
            json.dump({"pipeline_id": pipeline, "method": pipeline, "rollup_method": rollup,
                       "contrasts": CONTRASTS}, fh, indent=2)
    return tables


@unittest.skipUnless(have_r_ggplot(), "Rscript with ggplot2 not available")
class TopProteinViolins(unittest.TestCase):
    def run_figures(self, **kw):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        tables = write_inputs(d, **kw)
        out = os.path.join(d, "figs")
        r = subprocess.run(["Rscript", SCRIPT, "--de-dir", tables,
                            "--conditions", os.path.join(d, "conditions.csv"), "--outdir", out],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(out, "figures.json")) as fh:
            figs = json.load(fh)
        violins = {f["file"]: f["caption"] for f in figs if f["type"] == "violin"}
        return out, violins, r.stderr + r.stdout

    def line(self, log, contrast, what):
        tag = f"[figures] violin {make_name(contrast)} {what}: "
        hits = [ln[len(tag):] for ln in log.splitlines() if ln.startswith(tag)]
        self.assertEqual(len(hits), 1, f"no '{tag}' line in:\n{log}")
        return hits[0]

    def test_one_violin_per_contrast_and_listed(self):
        out, violins, log = self.run_figures()
        for c in CONTRASTS:
            fn = f"violin_top_{make_name(c)}.png"
            self.assertTrue(os.path.getsize(os.path.join(out, fn)) > 0, fn)
            self.assertIn(fn, violins)
        self.assertEqual(len(violins), len(CONTRASTS))
        self.assertNotIn("skipped", "".join(ln for ln in log.splitlines() if "violin" in ln))
        # the groups came from the formula, not from splitting the file name on "."
        self.assertIn("Treat.high_dose vs Ctrl", violins["violin_top_Treat.high_dose.Ctrl.png"])

    def test_detection_matrix_marks_inferred_values_and_detection_events(self):
        _, violins, log = self.run_figures()
        cap = violins["violin_top_Treat.high_dose.Ctrl.png"]
        sub = self.line(log, "Treat.high_dose-Ctrl", "subtitle")
        for text in (cap, sub):
            self.assertIn("inferred", text)
            self.assertIn("DPC-Quant", text)          # named by the pipeline's own record
        self.assertIn("hollow", sub)
        self.assertIn("GENE1 in Ctrl", sub)           # the never-measured group is named
        self.assertIn("detection event", cap)
        self.assertEqual(self.line(log, "Treat.high_dose-Ctrl", "status legend"), "Detected, Inferred")

    def test_without_detection_matrix_no_status_is_invented(self):
        _, violins, log = self.run_figures(detection=False)
        self.assertEqual(len(violins), len(CONTRASTS))
        for c in CONTRASTS:
            self.assertEqual(self.line(log, c, "status legend"), "none")
            sub = self.line(log, c, "subtitle")
            self.assertIn("not recorded", sub)
            self.assertNotIn("hollow", sub)
            cap = violins[f"violin_top_{make_name(c)}.png"]
            self.assertIn("not recorded", cap)
            self.assertNotIn("hollow", cap)
            self.assertNotIn("detection event", cap)

    def test_maxlfq_says_missing_not_inferred(self):
        _, violins, log = self.run_figures(pipeline="maxlfq")
        cap = violins["violin_top_Treat.high_dose.Ctrl.png"]
        sub = self.line(log, "Treat.high_dose-Ctrl", "subtitle")
        for text in (cap, sub):
            self.assertIn("Missing", text)
            self.assertNotIn("inferred", text.lower())
            self.assertNotIn("DPC", text)

    def test_no_provenance_derives_status_from_the_data_and_says_so(self):
        # A zero-detection cell that still holds a value can only have been inferred; the
        # wording must admit the pipeline was not recorded rather than name one.
        _, violins, log = self.run_figures(provenance=False)
        self.assertEqual(len(violins), len(CONTRASTS))  # groups found without the formula
        cap = violins["violin_top_Treat.high_dose.Ctrl.png"]
        self.assertIn("inferred", cap)
        self.assertIn("pipeline not recorded", cap)
        self.assertNotIn("DPC", cap)


if __name__ == "__main__":
    unittest.main()
