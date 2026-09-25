"""make_figures.R: the sample PCA circles every group and names it on the plot.

Groups are enclosed by a rounded hull (a circle around each sample, then the convex hull), which
works for any n -- circle for 1, capsule for 2, smooth hull for 3+ -- where a 95% normal ellipse
on n = 3 is unstable and huge. When the group names cross two factors (Old_JPH3 / Young_IgG =
Age x Bait) one factor becomes colour and the other point shape + outline style, but only when
the "_" split is unambiguous.

Needs Rscript (+ ggplot2 for the end-to-end cases); skipped otherwise (CI has no R).
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
MICHELLE_GROUPS = ["Old_IgG", "Old_JPH3", "Old_JPH4", "Old_Kv21", "Old_RyR",
                   "Young_IgG", "Young_JPH3", "Young_JPH4", "Young_Kv21", "Young_RyR"]


def have_r(pkg=None):
    if not shutil.which("Rscript"):
        return False
    if pkg is None:
        return True
    r = subprocess.run(["Rscript", "-e", f"quit(status = !requireNamespace('{pkg}', quietly=TRUE))"],
                       capture_output=True)
    return r.returncode == 0


def run_r(body):
    """Evaluate make_figures.R's top-level function definitions only, then `body`; return stdout."""
    code = f"""
      for (e in parse('{SCRIPT}')) if (is.call(e) && identical(e[[1]], as.name('<-')) &&
          (is.call(e[[3]]) && identical(e[[3]][[1]], as.name('function')) ||
           identical(e[[2]], as.name('GROUP_PAL')))) eval(e)
      {body}"""
    r = subprocess.run(["Rscript", "-e", code], capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise AssertionError(r.stderr)
    return r.stdout


def inside(px, py, poly):
    """Ray-casting point-in-polygon."""
    hit = False
    for (x1, y1), (x2, y2) in zip(poly, poly[1:]):
        if (y1 > py) != (y2 > py) and px < (x2 - x1) * (py - y1) / (y2 - y1) + x1:
            hit = not hit
    return hit


@unittest.skipUnless(have_r(), "Rscript not available")
class RoundedHull(unittest.TestCase):
    def hull(self, xs, ys, r=1.0):
        out = run_r(f"h <- rounded_hull(c({','.join(map(str, xs))}), c({','.join(map(str, ys))}), {r}); "
                    "cat(sprintf('%.6f %.6f', h$x, h$y), sep = '\\n')")
        return [tuple(map(float, ln.split())) for ln in out.strip().splitlines()]

    def test_closed_and_encloses_every_sample_for_n_1_2_3(self):
        for xs, ys in (([0], [0]), ([0, 5], [0, 1]), ([0, 5, 2], [0, 1, 4])):
            poly = self.hull(xs, ys)
            self.assertGreaterEqual(len(poly), 4, (xs, ys))
            self.assertEqual(poly[0], poly[-1], "polygon must be closed")
            for x, y in zip(xs, ys):
                self.assertTrue(inside(x, y, poly), f"sample ({x},{y}) outside its hull")

    def test_single_sample_is_a_circle_of_the_given_radius(self):
        poly = self.hull([3], [4], r=2.0)
        for x, y in poly:
            self.assertAlmostEqual(((x - 3) ** 2 + (y - 4) ** 2) ** 0.5, 2.0, places=5)


@unittest.skipUnless(have_r(), "Rscript not available")
class TwoFactorSplit(unittest.TestCase):
    def split(self, groups):
        out = run_r(f"tf <- two_factor_split(c({','.join(repr(g) for g in groups)})); "
                    "if (is.null(tf)) cat('NULL') else cat(paste(names(tf$colour), tf$colour, tf$style, sep = '|'), sep = '\\n')")
        if out.strip() == "NULL":
            return None
        return {ln.split("|")[0]: tuple(ln.split("|")[1:]) for ln in out.strip().splitlines()}

    def test_age_by_bait_is_detected(self):
        tf = self.split(MICHELLE_GROUPS)
        self.assertIsNotNone(tf)
        self.assertEqual(tf["Old_JPH3"], ("JPH3", "Old"))       # 5 baits = colour, 2 ages = shape
        self.assertEqual(tf["Young_IgG"], ("IgG", "Young"))

    def test_two_by_two(self):
        tf = self.split(["WT_ctrl", "WT_drug", "KO_ctrl", "KO_drug"])
        self.assertEqual(tf["KO_drug"], ("drug", "KO"))

    def test_not_split_when_it_is_not_a_clean_crossing(self):
        for groups in (["Ctrl", "Treat.high_dose", "Other"],            # no shared "_" structure
                       ["A_x", "B_y", "C_z", "D_w"],                    # 4 x 4 with 4 cells: not crossed
                       ["Old_JPH3", "Young_IgG", "Mid_RyR", "Old_IgG"],  # no 2-level factor
                       ["a_b_c", "d_e"]):                               # different token counts
            self.assertIsNone(self.split(groups), groups)


def write_pca_inputs(d, groups):
    """groups: {name: n_runs}. Expression only -- the PCA needs nothing else."""
    runs = [(f"20260925_Astral_DIA_{g}_run{i:02d}_final", g) for g, n in groups.items() for i in range(1, n + 1)]
    with open(os.path.join(d, "conditions.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        w.writerows(runs)
    tables = os.path.join(d, "tables")
    os.makedirs(tables)
    with open(os.path.join(tables, "Expression_Matrix.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Protein.Group", "Genes"] + [r for r, _ in runs])
        for p in range(40):
            w.writerow([f"P{p:05d}", f"G{p}"] +
                       [round(20 + ((p * 7 + s * 3 + sum(map(ord, g)) % 5) % 11) / 3, 3) for s, (_, g) in enumerate(runs)])
    return tables


@unittest.skipUnless(have_r("ggplot2"), "Rscript with ggplot2 not available")
class PcaFigure(unittest.TestCase):
    def render(self, groups):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        tables = write_pca_inputs(d, groups)
        out = os.path.join(d, "figs")
        r = subprocess.run(["Rscript", SCRIPT, "--de-dir", tables,
                            "--conditions", os.path.join(d, "conditions.csv"), "--outdir", out],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(r.returncode, 0, r.stderr)
        log = r.stderr + r.stdout
        self.assertNotIn("PCA failed", log)
        with open(os.path.join(out, "figures.json")) as fh:
            caps = {f["file"]: f["caption"] for f in json.load(fh)}
        self.assertIn("pca.png", caps)
        self.assertGreater(os.path.getsize(os.path.join(out, "pca.png")), 0)
        return caps["pca.png"], log

    def test_groups_of_one_two_and_three_render_with_a_two_factor_encoding(self):
        cap, log = self.render({"Old_A": 1, "Old_B": 2, "Young_A": 3, "Young_B": 2})
        self.assertIn("two factors", log)
        self.assertIn("cross two factors", cap)
        self.assertIn("Old: circle, solid outline", cap)
        self.assertIn("with a value in every sample", cap)   # how it was computed, from the code
        self.assertIn("scaled to unit variance", cap)

    def test_names_that_do_not_cross_get_one_colour_per_group(self):
        cap, log = self.render({"Ctrl": 3, "Treat": 3, "Other": 2})
        self.assertIn("one colour per group", log)
        self.assertNotIn("cross two factors", cap)
        self.assertIn("rounded hull", cap)


if __name__ == "__main__":
    unittest.main()
