"""make_figures.R: short sample names on the heatmap / PCA / QC plots.

Raw run names ("08132026__60SPD_DIA-LRS-124_S3-F4_1_23659") took ~40% of the heatmap and
buried the PCA (Silva08172026). short_sample_names() is the one definition: a Label/Sample
column in conditions.csv wins; otherwise the tokens every run shares are dropped and the
shortest distinguishing stretch is kept; any empty or colliding label falls back to the full
run name. The mapping is written to sample_labels.csv so no label is ambiguous.

Needs Rscript (+ ggplot2 for the end-to-end case); skipped otherwise (CI has no R).
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

MICHELLE = ["08132026__60SPD_DIA-LRS-96_S3-B1_1_23630", "08132026__60SPD_DIA-LRS-100_S3-F1_1_23648",
            "08132026__60SPD_DIA-LRS-124_S3-F4_1_23659", "08132026__60SPD_DIA-LRS-125_S3-G4_1_23639"]


def have_r(pkg=None):
    if not shutil.which("Rscript"):
        return False
    if pkg is None:
        return True
    r = subprocess.run(["Rscript", "-e", f"quit(status = !requireNamespace('{pkg}', quietly=TRUE))"],
                       capture_output=True)
    return r.returncode == 0


def short_names(runs, meta=None):
    """Run short_sample_names() from make_figures.R (only that definition is evaluated)."""
    d = tempfile.mkdtemp()
    try:
        with open(os.path.join(d, "runs.txt"), "w") as fh:
            fh.write("\n".join(runs) + "\n")
        meta_expr = "NULL"
        if meta is not None:
            with open(os.path.join(d, "meta.csv"), "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(list(meta[0].keys()))
                w.writerows([list(m.values()) for m in meta])
            meta_expr = f"read.csv('{d}/meta.csv', stringsAsFactors = FALSE, check.names = FALSE)"
        code = f"""
          for (e in parse('{SCRIPT}')) if (is.call(e) && identical(e[[1]], as.name('<-')) &&
              identical(e[[2]], as.name('short_sample_names'))) eval(e)
          res <- short_sample_names(readLines('{d}/runs.txt'), {meta_expr})
          write.csv(res, '{d}/out.csv', row.names = FALSE)"""
        r = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, timeout=120)
        if r.returncode != 0:
            raise AssertionError(r.stderr)
        with open(os.path.join(d, "out.csv")) as fh:
            return list(csv.DictReader(fh))
    finally:
        shutil.rmtree(d, True)


@unittest.skipUnless(have_r(), "Rscript not available")
class ShortSampleNames(unittest.TestCase):
    def test_shared_prefix_and_suffix_are_stripped(self):
        rows = short_names(MICHELLE)
        # "08132026__60SPD_DIA-" and "_S3-..._1_<id>" go; the bare number keeps its "LRS-"
        self.assertEqual([r["Label"] for r in rows], ["LRS-96", "LRS-100", "LRS-124", "LRS-125"])
        self.assertTrue(all(r["Source"] == "derived from the run name" for r in rows))

    def test_shortest_distinguishing_stretch_and_shared_tail_dropped(self):
        runs = [f"20260924_Exploris_DIA_{g}_S{i:02d}_Astral" for g in ("Ctrl", "Treat") for i in (1, 2, 3)]
        labels = [r["Label"] for r in short_names(runs)]
        self.assertEqual(labels[:2], ["Ctrl_S01", "Ctrl_S02"])   # "Ctrl" alone would not be unique
        self.assertEqual(len(set(labels)), len(labels))

    def test_labels_are_always_unique(self):
        runs = ["S1", "S1_extra", "a-b", "a_b", "x.1", "x.10", "plain"]   # awkward shapes on purpose
        labels = [r["Label"] for r in short_names(runs)]
        self.assertEqual(len(set(labels)), len(runs), labels)
        self.assertTrue(all(labels))

    def test_conditions_label_column_takes_precedence(self):
        meta = [{"File.Name": r, "Group": "A", "Label": f"mouse{i}"} for i, r in enumerate(MICHELLE, 1)]
        rows = short_names(MICHELLE, meta)
        self.assertEqual([r["Label"] for r in rows], ["mouse1", "mouse2", "mouse3", "mouse4"])
        self.assertTrue(all(r["Source"] == "conditions.csv Label" for r in rows))

    def test_empty_or_colliding_labels_fall_back_to_the_full_run_name(self):
        given = ["dup", "dup", "", "ok"]
        meta = [{"File.Name": r, "Group": "A", "Label": g} for r, g in zip(MICHELLE, given)]
        rows = short_names(MICHELLE, meta)
        self.assertEqual([r["Label"] for r in rows], MICHELLE[:3] + ["ok"])
        self.assertTrue(all("full run name" in r["Source"] for r in rows[:3]))


@unittest.skipUnless(have_r("ggplot2"), "Rscript with ggplot2 not available")
class SampleLabelsFile(unittest.TestCase):
    def test_figures_write_the_mapping_and_say_so(self):
        from test_make_figures_violin import write_inputs
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        tables = write_inputs(d)
        out = os.path.join(d, "figs")
        r = subprocess.run(["Rscript", SCRIPT, "--de-dir", tables,
                            "--conditions", os.path.join(d, "conditions.csv"), "--outdir", out],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(out, "sample_labels.csv")) as fh:
            rows = list(csv.DictReader(fh))
        with open(os.path.join(tables, "Expression_Matrix.csv")) as fh:
            runs = next(csv.reader(fh))[3:]
        self.assertEqual(sorted(r["File.Name"] for r in rows), sorted(runs))
        self.assertEqual(len({r["Label"] for r in rows}), len(runs))
        self.assertTrue(all(len(r["Label"]) < len(r["File.Name"]) for r in rows))
        with open(os.path.join(out, "figures.json")) as fh:
            caps = {f["file"]: f["caption"] for f in json.load(fh)}
        for fn in ("pca.png", "heatmap_top.png"):
            self.assertIn("sample_labels.csv", caps[fn])


if __name__ == "__main__":
    unittest.main()
