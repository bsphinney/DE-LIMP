#!/usr/bin/env python3
"""make_screen_data.R: the protein screen's core data (values, precursors, report rows, manifest).

End to end on a synthetic finished limpa analysis (screen_fixture.py; needs R with limpa/limma/
arrow/nanoparquet/jsonlite/statmod, skips otherwise). numpy/scipy are used when present.
NOTE: CI runs Python only and skips every R-invoking test -- run them locally.
"""
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)

try:
    import numpy as np
    from scipy import stats as sps
except ImportError:      # the tests that need them skip
    np = sps = None

import screen_fixture as fx  # noqa: E402
from screen_fixture import r_has, rscript  # noqa: E402

SCREEN = os.path.join(SCRIPTS, "make_screen_data.R")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def run_screen_data(session, *extra, env=None, proteins=None):
    """Run make_screen_data.R; the CompletedProcess gains `.out`, the folder it wrote.

    `env` adds environment variables. `proteins` writes a partial folder, into a fresh directory beside
    the session unless `extra` names an --out, so it never replaces the session's full screen_data.
    """
    args = list(extra)
    if proteins is not None:
        args += ["--proteins", proteins]
        if "--out" not in args:
            args += ["--out", tempfile.mkdtemp(prefix="screen_partial_", dir=os.path.dirname(session))]
    out = args[args.index("--out") + 1] if "--out" in args else os.path.join(session, "output", "screen_data")
    p = subprocess.run(["Rscript", SCREEN, "--session", session, *args], capture_output=True, text=True,
                       env=dict(os.environ, **(env or {})))
    p.out = out
    return p


def read_record(out, pg):
    """The record of protein group `pg` in the screen_data folder `out`."""
    with open(os.path.join(out, "index.tsv"), newline="") as fh:
        part = {r["Protein.Group"]: r["part"] for r in csv.DictReader(fh, delimiter="\t")}[pg]
    with gzip.open(os.path.join(out, "parts", part), "rt") as fh:
        return json.load(fh)[pg]


# Shares by nudging, of limpa's posterior mode converged tightly: the independent check of
# make_screen_data.R's analytic shares (design section 9b allows the derivative "provided a nudge check
# on a sample of proteins agrees"). Each fit is BFGS from limpa's own start, then restarted twice from
# its Newton-polished point; the derivative is a central difference of step 1e-3 (on the fixture it
# agrees with the exact derivative to ~4e-11). Priors and DPC from y_protein, sigma from
# dpcQuantHyperparam. Prints {pg: {own: [[precursor][sample]], other, total}}, precursors in raw_data order.
SHARES_REF_R = r'''
suppressMessages({ library(limma); library(limpa) })
s <- readRDS("@TABLES@/DE-LIMP_session.rds"); y <- s$raw_data; yp <- s$y_protein
pid <- as.character(y$genes$Protein.Group); o <- order(pid)
sig <- dpcQuantHyperparam(y[o, ], protein.id = pid[o], dpc.slope = yp$dpc[2])$sigma
conv <- function(E, sg) {
  f <- function(start) peptides2ProteinBFGS(E, sigma = sg, dpc = yp$dpc, prior.mean = yp$prior.mean,
         prior.sd = yp$prior.sd, prior.logFC = yp$prior.logFC, standard.errors = TRUE,
         newton.polish = TRUE, start = start)
  out <- f(NULL)
  for (r in 1:2) out <- f(drop(out$par))
  drop(out$protein.expression)
}
res <- lapply(setNames(nm = c(@PGS@)), function(pg) {
  E <- y$E[pid == pg, , drop = FALSE]; sg <- sig[[pg]]; ns <- ncol(E); np <- nrow(E); d <- 1e-3
  W <- array(0, c(ns, np, ns))
  for (p in seq_len(np)) for (k in seq_len(ns)) if (!is.na(E[p, k])) {
    up <- E; dn <- E; up[p, k] <- up[p, k] + d; dn[p, k] <- dn[p, k] - d
    W[, p, k] <- (conv(up, sg) - conv(dn, sg)) / (2 * d)
  }
  own <- matrix(vapply(seq_len(ns), function(j) W[j, , j], numeric(np)), np, ns); own[is.na(E)] <- NA
  total <- vapply(seq_len(ns), function(j) sum(W[j, , ]), 0)
  list(own = own, other = total - colSums(own, na.rm = TRUE), total = total)
})
cat(jsonlite::toJSON(res, digits = NA, na = "null"))
'''


def first_found(rec):
    """(precursor index, sample index) of the first found value of a record."""
    for p, prec in enumerate(rec["precursors"]):
        for j, v in enumerate(prec["log2"]):
            if v is not None:
                return p, j
    raise AssertionError("record %s has no found value" % rec["pg"])


def approx(x, tol=1e-6):
    class _A:
        def __eq__(self, other):
            return abs(other - x) <= tol * max(1.0, abs(x))
        __hash__ = None
    return _A()


@unittest.skipUnless(r_has(*fx.R_NEEDS), "needs R with " + "/".join(fx.R_NEEDS))
class ScreenData(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        os.makedirs(os.path.join(cls.tmp, "dpc"))
        s = fx.make_session(os.path.join(cls.tmp, "dpc"))
        cls.session, cls.output, cls.tables = s["session"], s["output"], s["tables"]
        cls.proc = run_screen_data(cls.session)
        cls.out = os.path.join(cls.output, "screen_data")
        if cls.proc.returncode != 0:
            raise AssertionError("make_screen_data.R failed:\n" + cls.proc.stderr[-3000:])
        with open(os.path.join(cls.out, "manifest.json")) as fh:
            cls.manifest = json.load(fh)
        tj = rscript(f'''suppressMessages(library(limma)); s <- readRDS("{cls.tables}/DE-LIMP_session.rds")
          E <- s$y_protein$E
          cat(jsonlite::toJSON(list(samples = colnames(E),
              rows = lapply(setNames(rownames(E), rownames(E)), function(p) unname(E[p, ]))),
              auto_unbox = TRUE, digits = NA))''')
        j = json.loads(tj)
        cls.y_samples, cls.y_protein = j["samples"], j["rows"]
        with open(os.path.join(cls.out, "index.tsv"), newline="") as fh:
            cls.index = {r["Protein.Group"]: r for r in csv.DictReader(fh, delimiter="\t")}
        cls.index_groups = set(cls.index)
        cls._maxlfq = None
        # limma's weights as run_de.R kept them (the unblocked dpc fit's vooma weights), and the
        # residual / prior degrees of freedom of the delivered fit; samples in y_protein's order
        wj = json.loads(rscript(f'''suppressMessages(library(limma))
          s <- readRDS("{cls.tables}/DE-LIMP_session.rds"); smp <- colnames(s$y_protein$E)
          w <- readRDS("{cls.tables}/set_test_inputs.rds")$weights$independent
          cat(jsonlite::toJSON(list(
              w = lapply(setNames(rownames(w), rownames(w)), function(p) unname(w[p, smp])),
              df_residual = as.list(setNames(s$fit$df.residual, rownames(s$fit$coefficients))),
              df_prior = s$fit$df.prior), auto_unbox = TRUE, digits = NA))'''))
        cls.vooma_weights, cls.df_residual, cls.df_prior = wj["w"], wj["df_residual"], wj["df_prior"]
        cls.groups = ["A"] * 3 + ["B"] * 3 + ["C"] * 3      # the fixture's design, samples in order

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def record(self, pg):
        return read_record(self.out, pg)

    def de_table(self, contrast):
        with open(os.path.join(self.tables, "DE_dpc_%s.csv" % contrast.replace("-", "."))) as fh:
            return {r["Protein.Group"]: {k: (float(v) if k == "logFC" and v != "NA" else v)
                                         for k, v in r.items()} for r in csv.DictReader(fh)}

    def test_runs_clean(self):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr[-2000:])

    def test_manifest_basics(self):
        m = self.manifest
        self.assertEqual(m["schema"], "screen_data/1")
        self.assertTrue(m["supported"])
        self.assertEqual(m["samples"], self.y_samples)
        self.assertEqual(m["groups"], ["A"] * 3 + ["B"] * 3 + ["C"] * 3)
        self.assertEqual([c["name"] for c in m["contrasts"]], ["B-A", "C-A"])
        self.assertEqual((m["contrasts"][0]["numerator"], m["contrasts"][0]["denominator"]), ("B", "A"))
        self.assertEqual(m["analysis"]["pipeline_id"], "dpc")
        # the page says which q-value checks and cutoffs ran, and what the DE engine was, from the analysis
        with open(os.path.join(self.tables, "de_provenance.json")) as fh:
            prov = json.load(fh)
        for k in ("q_columns", "q_cutoffs", "significance_rule", "de_engine", "rollup_method", "filters_applied"):
            self.assertEqual(m["analysis"][k], prov[k], k)
        self.assertEqual(len(m["dpc"]), 2)
        self.assertEqual(set(m["prior"]), {"prior.mean", "prior.sd", "prior.logFC"})
        self.assertEqual(m["index"], "index.tsv")
        self.assertTrue(m["written_by"].startswith("make_screen_data.R "))
        self.assertIn("tables/de_provenance.json", m["sources"])
        self.assertIn("search/report.stats.tsv", m["sources"])
        for rel, sha in m["sources"].items():
            self.assertEqual(sha256(os.path.join(self.output, rel)), sha)
        self.assertTrue(m["parts"])
        for rel, sha in m["parts"].items():
            self.assertEqual(sha256(os.path.join(self.out, rel)), sha)
        self.assertEqual(sorted(m["parts"]),
                         sorted({"parts/" + r["part"] for r in self.index.values()} | {"index.tsv"}))
        self.assertEqual(m["parts"]["index.tsv"], sha256(os.path.join(self.out, "index.tsv")))
        self.assertEqual({s["block"] for s in m["skipped"]},
                         {"xic", "library_uniqueness", "check_routes", "fran", "report_text"})
        # report.stats.tsv names runs by full raw-file path (the fixture's last is a Windows .d path)
        self.assertEqual(len(m["run_stats"]), 9)
        self.assertEqual(sorted(m["run_stats"]), sorted(m["samples"]))
        self.assertNotIn("run_stats", {s["block"] for s in m["skipped"]})
        self.assertIn("Precursors.Identified", m["run_stats"][m["samples"][-1]])

    @unittest.skipIf(np is None, "needs numpy")
    def test_values_match_delivered(self):
        r = self.record("P00001")
        self.assertEqual(len(r["values"]["protein"]), 9)
        np.testing.assert_allclose(r["values"]["protein"], self.y_protein["P00001"], atol=1e-4)
        self.assertEqual(r["de"]["B-A"]["logFC"], approx(self.de_table("B-A")["P00001"]["logFC"]))
        self.assertEqual(len(r["values"]["se"]), 9)
        self.assertEqual(len(r["values"]["n_obs"]), 9)
        self.assertGreater(len(r["precursors"]), 0)
        self.assertEqual(len(r["precursors"][0]["log2"]), 9)

    def test_ci_from_t_and_df(self):
        d = self.record("P00001")["de"]["B-A"]
        se = d["logFC"] / d["t"]
        if sps is not None:
            q = sps.t.ppf(0.975, d["df_total"])
        else:
            q = float(rscript(f'cat(qt(0.975, {d["df_total"]!r}))'))
        self.assertAlmostEqual(d["CI.L"], d["logFC"] - q * se, places=6)
        self.assertAlmostEqual(d["CI.R"], d["logFC"] + q * se, places=6)
        self.assertEqual(d["detected"], {"B": "3/3", "A": "3/3"})

    def test_failed_check_and_contaminant(self):
        rows = self.record("P00006")["rows"]
        bad = [x for x in rows if x["Run"] == "S-04.a" and x["Q.Value"] > 0.01]
        self.assertEqual(len(bad), 1)
        self.assertEqual(bad[0]["failed"], ["Q.Value"])
        self.assertFalse(bad[0]["used"])
        self.assertEqual(bad[0]["sample"], "S-04.a")
        self.assertFalse(bad[0]["contaminant"])
        self.assertTrue(any(x["used"] for x in rows))
        self.assertNotIn("P00090;Cont_P99999", self.index_groups)

    def test_untested_contrast_row(self):
        d = self.record("P00051")["de"]["B-A"]
        self.assertFalse(d["tested"])
        self.assertIsNone(d["logFC"])
        self.assertIn("not measured in either group", d["evidence"])
        self.assertTrue(self.record("P00051")["de"]["C-A"]["tested"])

    def test_refuses_maxlfq_session(self):
        d = os.path.join(self.tmp, "maxlfq")
        os.makedirs(d)
        s = fx.make_session(d, method="maxlfq")
        p = run_screen_data(s["session"])
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("limpa", p.stderr)
        with open(os.path.join(s["output"], "screen_data", "manifest.json")) as fh:
            m = json.load(fh)
        self.assertFalse(m["supported"])
        self.assertIn("limpa", m["unsupported_reason"])
        # the route is named from its own description (de_provenance.json), never hard-coded
        self.assertIn("MaxLFQ + limma (pipeline_id maxlfq", m["unsupported_reason"])
        self.assertIn("later part", m["unsupported_reason"])
        self.assertEqual(m["analysis"]["pipeline_id"], "maxlfq")
        self.assertEqual(m["analysis"]["display_label"], "MaxLFQ + limma")
        st = self.status_of(s["session"])                    # and the session's status record says refused
        self.assertEqual((st["state"], st["route"]), ("refused", "maxlfq"))
        self.assertIn("limpa", st["reason"])
        # and the page builder says so, naming the route, instead of drawing anything
        b = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_protein_screen.py"), "--screen-data",
                            os.path.join(s["output"], "screen_data"), "--protein", "P00001",
                            "--out", os.path.join(d, "P00001.html"), "--offline"], capture_output=True, text=True)
        self.assertEqual(b.returncode, 2, b.stderr[-1500:])
        self.assertIn("does not cover this analysis", b.stderr)
        self.assertIn("It used MaxLFQ + limma (pipeline_id maxlfq)", b.stderr)
        self.assertFalse(os.path.exists(os.path.join(d, "P00001.html")))

    def copied_session(self, name, repoint=True):
        """A copy of the session. A copy's de_provenance.json still names the ORIGINAL's report; `repoint`
        names the copy's own, as a DE run in the copy would (the sha256 is unchanged: same file)."""
        d = os.path.join(self.tmp, name)
        self.copy_to(d, repoint)
        shutil.rmtree(os.path.join(d, "output", "screen_data"), ignore_errors=True)
        return d

    def copy_to(self, d, repoint=True):
        shutil.copytree(self.session, d)
        if repoint:
            self.set_input(d, os.path.join(d, "output", "search", "report.parquet"))

    @staticmethod
    def set_input(d, path, drop=()):
        pf = os.path.join(d, "output", "tables", "de_provenance.json")
        with open(pf) as fh:
            prov = json.load(fh)
        prov["input"] = path
        for k in drop:
            prov.pop(k, None)
        with open(pf, "w") as fh:
            json.dump(prov, fh)

    def test_reads_the_report_the_de_read_from_its_own_record(self):
        # a session whose search folder is not output/search/ (search_dda/): the DE's record names the report
        # (input path + sha256) and the manifest says where it was read. The screen finds it from the record
        # alone; it does not look for a report by itself.
        d = self.copied_session("moved_search", repoint=False)
        os.rename(os.path.join(d, "output", "search"), os.path.join(d, "output", "search_dia"))
        self.set_input(d, os.path.join(d, "output", "search_dia", "report.parquet"))
        p = run_screen_data(d, proteins="P00001")
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.join(p.out, "manifest.json")) as fh:
            m = json.load(fh)
        self.assertEqual(m["files"]["report"], "search_dia/report.parquet")
        self.assertEqual(m["files"]["stats"], "search_dia/report.stats.tsv")
        self.assertIn("search_dia/report.parquet", m["sources"])
        self.assertIn("search_dia/report.stats.tsv", m["sources"])
        self.assertNotIn("search/report.parquet", m["sources"])

    def test_a_record_that_does_not_name_the_report_is_refused_not_searched(self):
        # an older analysis (no input path): the report is sitting in output/search/ with the right hash,
        # and the screen still refuses -- it will not re-derive which search the session chose
        d = self.copied_session("no_input", repoint=False)
        self.set_input(d, "", drop=("input",))
        p = run_screen_data(d)
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("does not record the report this analysis read", p.stderr)
        with open(os.path.join(d, "output", "tables", "de_provenance.json")) as fh:
            prov = json.load(fh)
        prov["input"] = os.path.join(d, "output", "search", "report.parquet")
        prov.pop("input_sha256")
        with open(os.path.join(d, "output", "tables", "de_provenance.json"), "w") as fh:
            json.dump(prov, fh)
        p = run_screen_data(d)
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("does not record the report this analysis read", p.stderr)

    def test_a_copied_session_is_accepted_by_the_same_relative_path(self):
        # the record still names the original's report (an absolute path outside this copy, which still
        # exists): the report at the same relative path under this session, with the DE's hash, is read
        d = self.copied_session("copied", repoint=False)
        p = run_screen_data(d, proteins="P00001")
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.join(p.out, "manifest.json")) as fh:
            m = json.load(fh)
        self.assertEqual(m["files"]["report"], "search/report.parquet")
        self.assertIn("search/report.parquet", m["sources"])
        self.assertFalse(any(os.path.isabs(x) for x in m["sources"]))

    def test_a_copied_session_whose_report_changed_is_refused(self):
        d = self.copied_session("copied_changed", repoint=False)
        with open(os.path.join(d, "output", "search", "report.parquet"), "ab") as fh:
            fh.write(b"x")
        p = run_screen_data(d)
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("the report changed since the DE ran", p.stderr)

    def test_a_copied_session_without_the_report_at_that_path_is_refused(self):
        # the report is in the session, but under another folder name: the screen does not go looking
        d = self.copied_session("copied_renamed", repoint=False)
        os.rename(os.path.join(d, "output", "search"), os.path.join(d, "output", "search_other"))
        p = run_screen_data(d)
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("outside this session's output/ folder", p.stderr)
        self.assertIn("same path under this session's output/ is not there", p.stderr)

    def test_a_record_naming_a_missing_report_is_refused(self):
        d = self.copied_session("gone", repoint=False)
        self.set_input(d, os.path.join(d, "output", "search_gone", "report.parquet"))
        p = run_screen_data(d)
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("is not at", p.stderr)

    def test_refuses_limpa_over_another_engines_report(self):
        # run_de.R --method dpc does not refuse an adapter's report; every block of part 1 reads DIA-NN's
        # precursor rows, so the screen refuses one that declares its own quantity, naming what it declares
        d = self.copied_session("declared")
        rp = os.path.join(d, "output", "search", "report.parquet")
        rscript(f'''suppressMessages(library(arrow)); tb <- read_parquet("{rp}", as_data_frame = FALSE)
          tb$metadata[["delimp.quantity.level"]] <- "protein"
          tb$metadata[["delimp.quantity.label"]] <- "FragPipe (IonQuant MaxLFQ)"
          tb$metadata[["delimp.quantity.value"]] <- "IonQuant MaxLFQ Intensity"
          write_parquet(tb, "{rp}")''')
        pf = os.path.join(d, "output", "tables", "de_provenance.json")
        with open(pf) as fh:
            prov = json.load(fh)
        prov["input"], prov["input_sha256"] = rp, sha256(rp)
        with open(pf, "w") as fh:
            json.dump(prov, fh)
        p = run_screen_data(d)
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("written by FragPipe (IonQuant MaxLFQ)", p.stderr)
        self.assertIn("not by DIA-NN", p.stderr)
        with open(os.path.join(d, "output", "screen_data", "manifest.json")) as fh:
            m = json.load(fh)
        self.assertFalse(m["supported"])
        self.assertEqual(m["analysis"]["pipeline_id"], "dpc")

    def test_changed_report_refused(self):
        d = self.copied_session("changed")
        with open(os.path.join(d, "output", "search", "report.parquet"), "ab") as fh:
            fh.write(b"x")
        p = run_screen_data(d)
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("is not the report this analysis read", p.stderr)
        self.assertIn("changed since the DE ran", p.stderr)

    def test_index_has_contrast_columns(self):
        """The index carries P.Value / adj.P.Val / logFC per contrast, equal to the record's, empty if untested."""
        with open(os.path.join(self.out, "index.tsv"), newline="") as fh:
            header = next(csv.reader(fh, delimiter="\t"))
        for c in self.manifest["contrasts"]:
            for col in ("P.Value", "adj.P.Val", "logFC"):
                self.assertIn("%s__%s" % (col, c["name"]), header)
        n_tested = {c["name"]: 0 for c in self.manifest["contrasts"]}
        for pg, row in self.index.items():
            rec = self.record(pg)
            for c in self.manifest["contrasts"]:
                de, n = rec["de"][c["name"]], c["name"]
                for col, key in (("P.Value", "P.Value"), ("adj.P.Val", "adj.P.Val"), ("logFC", "logFC")):
                    got = row["%s__%s" % (col, n)]
                    if de["tested"]:
                        self.assertAlmostEqual(float(got), de[key], places=9, msg=(pg, n, col))
                    else:
                        self.assertEqual(got, "", (pg, n, col))
                n_tested[n] += bool(de["tested"])
        self.assertTrue(all(n_tested.values()))

    def test_index_calls_match_manifest_n_sig(self):
        """Counting adj.P.Val < adjp in the index gives the manifest's n_sig (the rule is the table's)."""
        for c in self.manifest["contrasts"]:
            called = [r for r in self.index.values() if r["adj.P.Val__" + c["name"]] != ""
                      and float(r["adj.P.Val__" + c["name"]]) < c["adjp"]]
            self.assertEqual(len(called), c["n_sig"], c["name"])

    def test_partial(self):
        out = os.path.join(self.tmp, "partial")
        p = run_screen_data(self.session, "--out", out, "--proteins", "P00001,P00002")
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.join(out, "manifest.json")) as fh:
            self.assertTrue(json.load(fh)["partial"])
        with open(os.path.join(out, "index.tsv")) as fh:
            self.assertEqual(len(fh.read().split("\n")) - 2, 2)

    def test_proteins_refuses_the_default_out(self):
        # --proteins writes a partial folder: into the session's own screen_data it would replace the full one
        d = os.path.join(self.tmp, "default_out")
        self.copy_to(d)
        out = os.path.join(d, "output", "screen_data")

        def state():
            return {os.path.relpath(os.path.join(r, f), out): sha256(os.path.join(r, f))
                    for r, _, fs in os.walk(out) for f in fs}
        before = state()
        self.assertIn("manifest.json", before)
        for extra in ([], ["--out", out], ["--out", os.path.join(d, "output", "..", "output", "screen_data")],
                      ["--out", os.path.join(out, ".")]):
            p = run_screen_data(d, "--proteins", "P00001", *extra)
            self.assertEqual(p.returncode, 2, (extra, p.stderr[-1500:]))
            self.assertIn("--proteins", p.stderr)
            self.assertIn("--out", p.stderr)
            self.assertEqual(state(), before, extra)          # the full screen data is untouched
        # the path check alone (no screen data there yet): <S>/output/screen_data/. is the session's folder
        shutil.rmtree(out)
        os.makedirs(out)
        p = run_screen_data(d, "--proteins", "P00001", "--out", os.path.join(out, "."))
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("give it its own --out", p.stderr)
        self.assertEqual(os.listdir(out), [])

    def test_proteins_refuses_a_folder_of_full_screen_data(self):
        # wherever full screen data lies (another spelling, another letter case, a copy), --proteins must not
        # replace it: its manifest says supported and not partial
        full = os.path.join(self.tmp, "full_copy", "screen_data")
        shutil.copytree(self.out, full)
        with open(os.path.join(full, "manifest.json")) as fh:
            before = fh.read()
        p = run_screen_data(self.session, "--out", full, proteins="P00001")
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("holds full screen data", p.stderr)
        with open(os.path.join(full, "manifest.json")) as fh:
            self.assertEqual(fh.read(), before)
        self.assertEqual(sorted(os.listdir(os.path.join(full, "parts"))), sorted(os.listdir(os.path.join(self.out, "parts"))))

    def test_proteins_refuses_an_old_or_unreadable_manifest(self):
        # full screen data written before manifests carried `partial` has no such key, and a manifest that cannot
        # be read says nothing about what the folder holds: --proteins replaces neither
        for name, edit in (("no_partial_key", lambda m: {k: v for k, v in json.loads(m).items() if k != "partial"}),
                           ("unreadable", None)):
            full = os.path.join(self.tmp, name, "screen_data")
            shutil.copytree(self.out, full)
            mf = os.path.join(full, "manifest.json")
            with open(mf) as fh:
                text = fh.read()
            text = json.dumps(edit(text)) if edit else text[: len(text) // 2]
            with open(mf, "w") as fh:
                fh.write(text)
            p = run_screen_data(self.session, "--out", full, proteins="P00001")
            self.assertEqual(p.returncode, 2, (name, p.stderr[-1500:]))
            self.assertIn("own --out", p.stderr, name)
            with open(mf) as fh:
                self.assertEqual(fh.read(), text, name)
            self.assertEqual(sorted(os.listdir(os.path.join(full, "parts"))),
                             sorted(os.listdir(os.path.join(self.out, "parts"))), name)

    def test_an_argument_error_keeps_the_manifest(self):
        # a typo in --proteins refuses, and leaves the screen data already in --out (supported) as it was
        out = os.path.join(self.tmp, "partial_keep")
        p = run_screen_data(self.session, "--out", out, proteins="P00001")
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.join(out, "manifest.json")) as fh:
            before = fh.read()
        self.assertTrue(json.loads(before)["supported"])
        p = run_screen_data(self.session, "--out", out, proteins="P0000X")
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("P0000X", p.stderr)
        self.assertIn("manifest.json left as it was", p.stderr)
        with open(os.path.join(out, "manifest.json")) as fh:
            self.assertEqual(fh.read(), before)
        # into a fresh folder an argument error writes nothing that could pass for screen data, and says so
        fresh = os.path.join(self.tmp, "partial_typo")
        p = run_screen_data(self.session, "--out", fresh, proteins="P0000X")
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("(nothing written)", p.stderr)
        self.assertNotIn("left as it was", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(fresh, "manifest.json")))

    def test_shares_shape_and_gaps(self):
        r = self.record("P00002")
        own = r["shares"]["own"]
        self.assertEqual(len(own), len(r["precursors"]))
        for p, prec in enumerate(r["precursors"]):
            for j, v in enumerate(prec["log2"]):
                self.assertEqual(own[p][j] is None, v is None)      # a gap carries no own share
        for j in range(9):
            s = sum(x[j] for x in own if x[j] is not None) + r["shares"]["other"][j]
            self.assertAlmostEqual(s, r["shares"]["total"][j], places=3)

    def test_refit_reproduces_delivered(self):
        self.assertEqual(self.manifest["shares"]["withheld"], [])
        self.assertLessEqual(self.record("P00002")["shares"]["refit_max_diff"], 1e-6)

    def test_refit_mismatch_withholds(self):
        # SCREEN_DATA_SIGMA_SCALE=1.01 is a test-only env hook that perturbs sigma before refitting
        p = run_screen_data(self.session, env={"SCREEN_DATA_SIGMA_SCALE": "1.01"}, proteins="P00002")
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        rec = read_record(p.out, "P00002")
        self.assertIn("withheld", rec["shares"])
        self.assertIn("refit", rec["shares"]["withheld"])
        # ... and the manifest lists it with the same reason, and says the hook was set
        with open(os.path.join(p.out, "manifest.json")) as fh:
            sh = json.load(fh)["shares"]
        self.assertEqual(sh["withheld"], [{"pg": "P00002", "reason": rec["shares"]["withheld"]}])
        self.assertEqual(sh["test_hook_sigma_scale"], 1.01)
        # the detection chances come from the same fit, so they are withheld with the same reason
        self.assertEqual(rec["detection"], {"withheld": rec["shares"]["withheld"]})
        with open(os.path.join(p.out, "manifest.json")) as fh:
            self.assertEqual(json.load(fh)["detection"]["withheld"], sh["withheld"])

    def test_hyperparameter_mismatch_withholds_every_protein(self):
        # The cross-CPU-family path: limpa's hyperparameters recomputed here differ from the delivered ones (as
        # dpcQuant's do between zen2 and zen4), so no refit describes the analysis and every protein's shares and
        # detection chances are withheld, the reason naming both families. The delivered prior.sd is nudged by
        # 1e-6; the DE is recorded as run on zen2, and this "node" is zen4 (a stand-in scontrol on PATH, read the
        # way compute_record.R reads SLURM's).
        d = os.path.join(self.tmp, "hyper")
        self.copy_to(d)
        tables = os.path.join(d, "output", "tables")
        f = os.path.join(tables, "DE-LIMP_session.rds")
        rscript(f'''suppressMessages(library(limma)); s <- readRDS("{f}")
          s$y_protein$prior.sd <- s$y_protein$prior.sd + 1e-6; saveRDS(s, "{f}")''')
        pf = os.path.join(tables, "de_provenance.json")
        with open(pf) as fh:
            prov = json.load(fh)
        prov["compute"] = dict(prov.get("compute") or {}, cpu_family="zen2")
        with open(pf, "w") as fh:
            json.dump(prov, fh)
        fake_bin = tempfile.mkdtemp(dir=self.tmp)
        with open(os.path.join(fake_bin, "scontrol"), "w") as fh:
            fh.write('#!/bin/sh\necho "NodeName=fake Arch=x86_64 AvailableFeatures=amd,zen4 ActiveFeatures=amd,zen4"\n')
        os.chmod(os.path.join(fake_bin, "scontrol"), 0o755)
        pgs = ["P00001", "P00002", "P00051"]
        p = run_screen_data(d, proteins=",".join(pgs),
                            env={"PATH": fake_bin + os.pathsep + os.environ["PATH"], "SLURMD_NODENAME": "fake"})
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        with open(os.path.join(p.out, "manifest.json")) as fh:
            m = json.load(fh)
        self.assertEqual(m["compute"], {"de": "zen2", "here": "zen4"})
        why = {w["pg"]: w["reason"] for w in m["shares"]["withheld"]}
        self.assertEqual(sorted(why), pgs)                                  # every protein
        self.assertEqual(len(set(why.values())), 1)
        reason = why["P00001"]
        self.assertIn("hyperparameters", reason)
        self.assertIn("compute.de zen2", reason)
        self.assertIn("compute.here zen4", reason)
        self.assertEqual(m["detection"]["withheld"], m["shares"]["withheld"])
        self.assertIn(reason, m["shares"]["nudge_check"]["reason"])
        for pg in pgs:
            rec = read_record(p.out, pg)
            self.assertEqual(rec["shares"], {"withheld": reason})
            self.assertEqual(rec["detection"], {"withheld": reason})

    def test_shares_match_converged_reference(self):
        # the page shows shares to 4 decimals: they must be right to that, not carry BFGS's stopping
        # point (plain-BFGS nudges at 0.01 put P00031's and P00037's shares ~0.02 off). P00052/P00055
        # are C-only (6 of 9 samples gaps), P00006/P00021 have a few gaps.
        self.assertEqual(self.manifest["shares"]["method"], "derivative")
        pgs = ["P00002", "P00006", "P00021", "P00031", "P00037", "P00052", "P00055"]
        ref = json.loads(rscript(SHARES_REF_R.replace("@TABLES@", self.tables)
                                 .replace("@PGS@", ", ".join('"%s"' % p for p in pgs))))
        for pg in pgs:
            sh, rf = self.record(pg)["shares"], ref[pg]
            for p, row in enumerate(sh["own"]):
                for j, v in enumerate(row):
                    if v is None:
                        self.assertIsNone(rf["own"][p][j])
                    else:
                        self.assertLessEqual(abs(v - rf["own"][p][j]), 1e-4, (pg, "own", p, j))
            for key in ("other", "total"):
                for j, v in enumerate(sh[key]):
                    self.assertLessEqual(abs(v - rf[key][j]), 1e-4, (pg, key, j))

    def test_refuses_row_by_row_analysis(self):
        # every protein one precursor: limpa used dpcQuantByRow, which the screen does not cover
        d = os.path.join(self.tmp, "byrow")
        os.makedirs(d)
        s = fx.make_session(d, precursors=1)
        # the premise: limpa did quantify row by row, so y_protein's rows are precursors, not proteins
        rows = rscript(f'''suppressMessages(library(limma))
          cat(rownames(readRDS("{s["tables"]}/DE-LIMP_session.rds")$y_protein$E), sep = "\\n")''').split()
        self.assertTrue(rows)
        self.assertEqual([r for r in rows if not r.startswith(("PEP", "CONT"))], [])
        p = run_screen_data(s["session"])
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("dpcQuantByRow", p.stderr)
        with open(os.path.join(s["output"], "screen_data", "manifest.json")) as fh:
            m = json.load(fh)
        self.assertFalse(m["supported"])
        self.assertIn("dpcQuantByRow", m["unsupported_reason"])
        self.assertIn("row by row", m["unsupported_reason"])
        self.assertFalse(os.path.exists(os.path.join(s["output"], "screen_data", "parts")))  # no records

    def test_precursor_count_is_limpas(self):
        # limpa >= 1.4 names the per-protein count genes$NPrec (older limpa: NPeptides). It counts the protein's
        # precursors (its rows in raw_data), so the record calls it n_precursors, never n_peptides.
        cols = rscript(f'''suppressMessages(library(limma))
          cat(names(readRDS("{self.tables}/DE-LIMP_session.rds")$y_protein$genes), sep = " ")''').split()
        self.assertTrue({"NPrec", "NPeptides"} & set(cols), cols)
        r = self.record("P00001")
        self.assertNotIn("n_peptides", r)
        self.assertEqual(r["n_precursors"], len(r["precursors"]))
        self.assertNotIn("n_precursors", {s["block"] for s in self.manifest["skipped"]})

    def test_precursor_count_of_older_limpa_or_none(self):
        d = os.path.join(self.tmp, "nprec")
        self.copy_to(d)
        f = os.path.join(d, "output", "tables", "DE-LIMP_session.rds")

        def rename(to):
            rscript(f'''suppressMessages(library(limma)); s <- readRDS("{f}"); g <- s$y_protein$genes
              names(g)[names(g) %in% c("NPrec", "NPeptides")] <- "{to}"; s$y_protein$genes <- g; saveRDS(s, "{f}")''')
        rename("NPeptides")                                    # an analysis made with limpa < 1.4
        p = run_screen_data(d, proteins="P00001")
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        rec = read_record(p.out, "P00001")
        self.assertEqual(rec["n_precursors"], len(rec["precursors"]))
        rename("SomethingElse")                                # neither: null, and the manifest says why
        p = run_screen_data(d, proteins="P00001")
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertIsNone(read_record(p.out, "P00001")["n_precursors"])
        with open(os.path.join(p.out, "manifest.json")) as fh:
            sk = {s["block"]: s["reason"] for s in json.load(fh)["skipped"]}
        self.assertIn("NPrec", sk["n_precursors"])
        self.assertIn("NPeptides", sk["n_precursors"])

    def test_sigma_is_limpas(self):
        # the protein's own sigma from limpa's dpcQuantHyperparam, as dpcQuant computes it
        sg = float(rscript(f'''suppressMessages({{library(limma); library(limpa)}})
          s <- readRDS("{self.tables}/DE-LIMP_session.rds"); y <- s$raw_data
          pid <- as.character(y$genes$Protein.Group); o <- order(pid)
          h <- dpcQuantHyperparam(y[o, ], protein.id = pid[o], dpc.slope = s$y_protein$dpc[2])
          cat(sprintf("%.15g", h$sigma[["P00002"]]))'''))
        r = self.record("P00002")
        self.assertIsInstance(r["sigma"], float)
        self.assertAlmostEqual(r["sigma"], sg, places=12)
        self.assertNotIn("test_hook_sigma_scale", self.manifest["shares"])
        self.assertNotIn("test_hook_derivative_scale", self.manifest["shares"])
        self.assertIn("shares", self.manifest["elapsed_s"])

    def test_nudge_check_recorded(self):
        # every run checks 3 proteins' derivatives against nudges of fully converged fits, and records
        # the limpa that ran. On the fixture the gappy, several-precursor proteins are the C-only ones.
        sh = self.manifest["shares"]
        ver = rscript('cat(as.character(packageVersion("limpa")))').strip()
        self.assertEqual(sh["limpa_version"], ver)
        nc = sh["nudge_check"]
        self.assertEqual(nc["proteins"], ["P00051", "P00052", "P00053"])
        self.assertTrue(nc["passed"])
        self.assertLessEqual(nc["max_diff"], 1e-4)
        self.assertEqual(nc["tolerance"], 1e-4)
        self.assertEqual(nc["cells_checked"], nc["cells_found"])      # every found value was nudged
        self.assertGreater(nc["cells_checked"], 0)

    def test_nudge_check_failure_withholds_everything(self):
        # SCREEN_DATA_DERIVATIVE_SCALE (a test-only hook) makes every derivative 1% wrong, as a change in
        # limpa's objective would: the run-time check must catch it and withhold every protein's shares
        p = run_screen_data(self.session, env={"SCREEN_DATA_DERIVATIVE_SCALE": "1.01"},
                            proteins="P00002,P00006,P00051,P00052,P00053")
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertIn("nudge check failed", p.stderr)
        with open(os.path.join(p.out, "manifest.json")) as fh:
            sh = json.load(fh)["shares"]
        nc = sh["nudge_check"]
        self.assertFalse(nc["passed"])
        self.assertGreater(nc["max_diff"], 1e-4)
        self.assertIn("nudge check failed", nc["reason"])
        self.assertEqual(sh["test_hook_derivative_scale"], 1.01)
        pgs = ["P00002", "P00006", "P00051", "P00052", "P00053"]
        self.assertEqual(sh["withheld"], [{"pg": pg, "reason": nc["reason"]} for pg in pgs])
        for pg in pgs:
            self.assertEqual(read_record(p.out, pg)["shares"], {"withheld": nc["reason"]})
            self.assertEqual(read_record(p.out, pg)["detection"], {"withheld": nc["reason"]})

    def test_a_killed_worker_stops_the_run(self):
        # SCREEN_DATA_DROP_RECORD (a test-only hook) empties the 2nd protein's record after the workers
        # return, as mclapply leaves the record of a worker that was killed (out of memory): the run must
        # stop, never file the other records under the wrong keys and exit 0
        p = run_screen_data(self.session, env={"SCREEN_DATA_DROP_RECORD": "2"}, proteins="P00001,P00002,P00003")
        self.assertNotEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertIn("P00002", p.stderr)
        self.assertIn("returned nothing", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(p.out, "manifest.json")))     # nothing claims to be screen data

    @staticmethod
    def status_of(session):
        with open(os.path.join(session, "output", "screen_data_status.json")) as fh:
            return json.load(fh)

    def test_status_record_built(self):
        st = self.status_of(self.session)
        self.assertEqual((st["schema"], st["state"], st["route"]), ("screen_data_status/1", "built", "dpc"))
        self.assertIsNone(st["reason"])

    def test_status_record_failed_with_one_line(self):
        d = self.copied_session("status_failed")
        p = run_screen_data(d, env={"SCREEN_DATA_DROP_RECORD": "2"})
        self.assertEqual(p.returncode, 1, p.stderr[-1500:])
        st = self.status_of(d)
        self.assertEqual(st["state"], "failed")
        self.assertIn("returned nothing", st["reason"])
        self.assertNotIn("\n", st["reason"])
        self.assertEqual(st["route"], "dpc")
        # and the page builder, on the folder the crashed rebuild left, says what happened
        sd = os.path.join(d, "output", "screen_data")
        self.assertFalse(os.path.exists(os.path.join(sd, "manifest.json")))
        b = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_protein_screen.py"), "--screen-data", sd,
                            "--protein", "P00001", "--offline", "--out", os.path.join(d, "x.html")],
                           capture_output=True, text=True)
        self.assertEqual(b.returncode, 2, b.stderr)
        self.assertIn("the last build of the screen data failed (", b.stderr)
        self.assertIn("returned nothing", b.stderr)
        self.assertNotIn("made before", b.stderr)

    def test_status_record_is_the_sessions_own_only(self):
        # a partial / --out run describes no session screen data: it leaves the record as it was
        before = self.status_of(self.session)
        self.assertEqual(run_screen_data(self.session, proteins="P00001").returncode, 0)
        self.assertEqual(self.status_of(self.session), before)
        d = self.copied_session("status_partial_only")
        os.remove(os.path.join(d, "output", "screen_data_status.json"))
        self.assertEqual(run_screen_data(d, proteins="P00001").returncode, 0)     # into its own --out
        self.assertFalse(os.path.exists(os.path.join(d, "output", "screen_data_status.json")))

    def test_a_killed_rebuild_leaves_no_old_manifest(self):
        # a rebuild that empties parts/ also removes the old manifest.json and index.tsv: a run killed midway
        # then reads "has no manifest.json ... Run make_screen_data.R", not "part-000N is missing"
        out = os.path.join(self.tmp, "rebuild_killed")
        pgs = "P00001,P00002,P00003"
        self.assertEqual(run_screen_data(self.session, "--out", out, proteins=pgs).returncode, 0)
        self.assertTrue(os.path.isfile(os.path.join(out, "manifest.json")))
        p = run_screen_data(self.session, "--out", out, env={"SCREEN_DATA_DROP_RECORD": "2"}, proteins=pgs)
        self.assertNotEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertFalse(os.path.exists(os.path.join(out, "manifest.json")))
        self.assertFalse(os.path.exists(os.path.join(out, "index.tsv")))
        b = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_protein_screen.py"), "--screen-data", out,
                            "--protein", "P00001", "--offline", "--out", os.path.join(self.tmp, "killed.html")],
                           capture_output=True, text=True)
        self.assertEqual(b.returncode, 2, b.stderr)
        self.assertIn("has no manifest.json", b.stderr)
        self.assertIn("Run make_screen_data.R", b.stderr)

    def test_detection_ranges(self):
        d = self.record("P00002")["detection"]
        for p, prec in enumerate(self.record("P00002")["precursors"]):
            for j, v in enumerate(prec["log2"]):
                self.assertTrue(0 <= d["pf_fit"][p][j] <= 1)
                self.assertEqual(d["pf_obs"][p][j] is None, v is None)
        det = self.manifest["detection"]
        self.assertEqual(det["withheld"], [])
        self.assertEqual(det["quadrature_points"], 16)
        self.assertLessEqual(det["objective_max_diff"], det["objective_tolerance"])   # limpa's own objective
        self.assertEqual(det["dpc"], self.manifest["dpc"])     # the curve limpa used is the analysis's

    def test_found_chance_is_dpc_at_value(self):
        r = self.record("P00002"); b0, b1 = self.manifest["dpc"]
        p, j = first_found(r)
        self.assertAlmostEqual(r["detection"]["pf_obs"][p][j],
                               1 / (1 + math.exp(-(b0 + b1 * r["precursors"][p]["log2"][j]))), places=3)

    @unittest.skipIf(np is None, "needs numpy")
    def test_weights_from_the_fit(self):
        w = self.record("P00001")["weights"]
        np.testing.assert_allclose(w["w"], self.vooma_weights["P00001"], rtol=1e-6)   # set_test_inputs.rds weights$independent
        for g in ("A", "B", "C"):
            self.assertAlmostEqual(sum(s for s, gg in zip(w["share_in_group"], self.groups) if gg == g), 1.0, places=6)
        self.assertEqual(w["df_residual"], self.df_residual["P00001"])
        self.assertEqual(w["df_prior"], approx(self.df_prior, 1e-12))
        # a protein with a group of inferred values: limma leaves that group out of its residual df
        self.assertLess(self.record("P00051")["weights"]["df_residual"], w["df_residual"])
        self.assertEqual(self.record("P00051")["weights"]["df_residual"], self.df_residual["P00051"])

    def test_weights_absent_is_skipped(self):
        # a set_test_inputs.rds without the fit's weights: weights are null, with the reason in
        # manifest.skipped -- never refitted
        d = os.path.join(self.tmp, "noweights")
        self.copy_to(d)
        f = os.path.join(d, "output", "tables", "set_test_inputs.rds")
        rscript(f'x <- readRDS("{f}"); x$weights$independent <- NULL; saveRDS(x, "{f}")')
        p = run_screen_data(d, proteins="P00001")
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertIsNone(read_record(p.out, "P00001")["weights"])
        with open(os.path.join(p.out, "manifest.json")) as fh:
            sk = {s["block"]: s["reason"] for s in json.load(fh)["skipped"]}
        self.assertIn("weights", sk)
        self.assertIn("weights$independent", sk["weights"])


# The derivative at a realistic detection curve. The fixture's DPC slope (~1e-4) leaves the detection
# term almost inert, so this checks screen_shares.R on one synthetic protein with DPC (-9, 0.7), gaps
# in several samples and one sample with a single found value, against a central difference of fully
# converged limpa fits written here independently. Prints max |analytic - nudge| over the derivatives
# and over own/other/total, the largest derivative, and the seconds the analytic part took.
DERIV_R = r'''
suppressMessages(library(limpa))
source("@SCRIPTS@/screen_shares.R")
set.seed(7)
ns <- 6; np <- 5
E <- matrix(rnorm(np * ns, 16, 1), np, ns) + rep(c(0, 0, 0, 1, 1, 1), each = np)
E[1, c(2, 5)] <- NA; E[3, 1:3] <- NA; E[5, 6] <- NA; E[, 4] <- c(NA, NA, NA, NA, 15)
h <- list(dpc = c(-9, 0.7), prior.mean = 15, prior.sd = 2, prior.logFC = 1.2); sg <- 0.4
t0 <- proc.time()[["elapsed"]]
d <- shares_derivative(E, sg, h, shares_fit(E, sg, h))
a <- shares_split(d$Wm, E)
secs <- proc.time()[["elapsed"]] - t0
conv <- function(E) {
  f <- peptides2ProteinBFGS(E, sigma = sg, dpc = h$dpc, prior.mean = h$prior.mean, prior.sd = h$prior.sd,
                            prior.logFC = h$prior.logFC, standard.errors = TRUE, newton.polish = TRUE)
  for (r in 1:4) f <- peptides2ProteinBFGS(E, sigma = sg, dpc = h$dpc, prior.mean = h$prior.mean,
    prior.sd = h$prior.sd, prior.logFC = h$prior.logFC, standard.errors = TRUE, newton.polish = TRUE,
    start = drop(f$par))
  drop(f$par[1:ns])
}
Wn <- matrix(0, ns, np * ns); st <- 1e-4
for (i in which(!is.na(as.vector(E)))) {
  up <- E; dn <- E; up[i] <- up[i] + st; dn[i] <- dn[i] - st
  Wn[, i] <- (conv(up) - conv(dn)) / (2 * st)
}
b <- shares_split(Wn, E)
cat(max(abs(d$Wm - Wn)), max(abs(a$own - b$own), abs(a$other - b$other), abs(a$total - b$total), na.rm = TRUE),
    max(abs(d$Wm)), secs, "\n")
'''


@unittest.skipUnless(r_has("limpa", "statmod"), "needs R with limpa/statmod")
class SharesDerivative(unittest.TestCase):
    def test_derivative_at_a_realistic_detection_curve(self):
        w, s, wmax, secs = map(float, rscript(DERIV_R.replace("@SCRIPTS@", SCRIPTS)).split())
        self.assertLessEqual(w, 1e-4)        # every d value_j / d E[p, k]
        self.assertLessEqual(s, 1e-4)        # own / other / total
        self.assertGreater(wmax, 0.1)        # the derivatives are not trivially small
        self.assertLess(secs, 1.0)           # the analytic part, for this protein


# The detection chances at a realistic detection curve (the fixture's DPC slope ~1e-4 makes every
# chance ~0.93), on DERIV_R's protein. Written here independently of screen_shares.R: the fitted
# precursor level from the contr.sum offsets, and the chance of finding as a numerical integral of
# the curve over N(level, sigma). Then limpa's own objective at its fit, rebuilt from the chances,
# must equal limpa's value there; given a wrong value, detection_at_fit must withhold. Prints
# max |pf_fit - integral|, the relative objective difference, max |pf_obs - curve at value|, whether
# pf_obs is NA exactly at the gaps, the gaps' smallest and largest chance, and whether the wrong
# value was withheld.
DETECT_R = r'''
suppressMessages(library(limpa))
source("@SCRIPTS@/screen_shares.R")
set.seed(7)
ns <- 6; np <- 5
E <- matrix(rnorm(np * ns, 16, 1), np, ns) + rep(c(0, 0, 0, 1, 1, 1), each = np)
E[1, c(2, 5)] <- NA; E[3, 1:3] <- NA; E[5, 6] <- NA; E[, 4] <- c(NA, NA, NA, NA, 15)
h <- list(dpc = c(-9, 0.7), prior.mean = 15, prior.sd = 2, prior.logFC = 1.2); sg <- 0.4
f <- shares_fit(E, sg, h)
d <- detection_at_fit(E, f, sg, h)
b <- drop(f$par.bgfs); off <- b[-seq_len(ns)]
lvl <- outer(c(off, -sum(off)), b[seq_len(ns)], "+")
pint <- matrix(vapply(lvl, function(m) integrate(function(x) plogis(h$dpc[1] + h$dpc[2] * x) * dnorm(x, m, sg),
                                               m - 12 * sg, m + 12 * sg, rel.tol = 1e-12)$value, 0), np, ns)
f_bad <- f; f_bad$value.bfgs <- f$value.bfgs * (1 + 1e-6)
bad <- detection_at_fit(E, f_bad, sg, h)
cat(max(abs(d$pf_fit - pint)), d$objective_diff, max(abs(d$pf_obs - plogis(-9 + 0.7 * E)), na.rm = TRUE),
    identical(is.na(d$pf_obs), is.na(E)), range(d$pf_fit[is.na(E)]), is.character(bad$withheld), "\n")
'''


@unittest.skipUnless(r_has("limpa", "statmod"), "needs R with limpa/statmod")
class DetectionChances(unittest.TestCase):
    def test_detection_at_a_realistic_curve(self):
        out = rscript(DETECT_R.replace("@SCRIPTS@", SCRIPTS)).split()
        fit_diff, obj_diff, obs_diff = map(float, out[:3])
        lo, hi = map(float, out[4:6])
        self.assertLessEqual(fit_diff, 1e-8)       # 16-point quadrature of the logistic-normal integral
        self.assertLessEqual(obj_diff, 1e-10)      # limpa's own objective, rebuilt from the chances
        self.assertLessEqual(obs_diff, 1e-12)
        self.assertEqual(out[3], "TRUE")           # pf_obs only at found values
        self.assertGreater(hi - lo, 0.1)           # the gaps' chances are not all alike
        self.assertEqual(out[6], "TRUE")           # a fit whose objective is not limpa's is withheld


@unittest.skipUnless(r_has(*fx.R_NEEDS), "needs R with " + "/".join(fx.R_NEEDS))   # make_session builds its fixture in R
class DesignRefusals(unittest.TestCase):
    """screen_support() refuses what the page cannot describe, with its reason (M1, M2)."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def refused(self, name, **kw):
        d = os.path.join(self.tmp, name)
        os.makedirs(d)
        s = fx.make_session(d, **kw)
        p = run_screen_data(s["session"])
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertFalse(os.path.exists(os.path.join(s["output"], "screen_data", "parts")))   # no records
        with open(os.path.join(s["output"], "screen_data", "manifest.json")) as fh:
            m = json.load(fh)
        self.assertFalse(m["supported"])
        self.assertEqual(m["analysis"]["pipeline_id"], "dpc")        # named from the analysis, like a route
        return p, m

    def test_refuses_a_contrast_that_is_not_one_group_against_one(self):
        # limma takes (B+C)/2-A; the page needs a numerator and a denominator group
        p, m = self.refused("combo", contrasts=["(B+C)/2-A"])
        self.assertIn("one group with one other", m["unsupported_reason"])
        self.assertIn("'(B+C)/2-A'", m["unsupported_reason"])
        self.assertIn("later part", m["unsupported_reason"])

    def test_refuses_when_any_one_contrast_is_not_pairwise(self):
        p, m = self.refused("mixed", contrasts=["B-A", "(B+C)/2-A"])
        self.assertIn("'(B+C)/2-A'", m["unsupported_reason"])
        self.assertNotIn("'B-A'", m["unsupported_reason"])

    def test_refuses_a_blocked_design(self):
        # a --block fit's numbers are not the unblocked fit's; the page would show the unblocked fit's weights
        # beside them. The refusal reads de_provenance.json block.applied
        p, m = self.refused("blocked", block_col="Mouse", de_args=["--block", "Mouse"])
        self.assertIn("blocked or with covariates", m["unsupported_reason"])
        self.assertIn("later part", m["unsupported_reason"])

    def test_refuses_a_design_with_covariates(self):
        # a Batch column the DE models as a covariate (covariates.terms in de_provenance.json, written by run_de.R)
        d = os.path.join(self.tmp, "covar")
        os.makedirs(d)
        s = fx.make_session(d, block_col="Batch")
        with open(os.path.join(s["tables"], "de_provenance.json")) as fh:
            self.assertTrue(json.load(fh)["covariates"]["terms"])      # the premise: the DE kept a term
        p = run_screen_data(s["session"])
        self.assertEqual(p.returncode, 2, p.stderr[-1500:])
        self.assertIn("blocked or with covariates", p.stderr)


if __name__ == "__main__":
    unittest.main(verbosity=2)
