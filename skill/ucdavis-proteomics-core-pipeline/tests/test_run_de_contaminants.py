#!/usr/bin/env python3
"""run_de.R removes contaminants before quantification, records it, and the Methods follow.

msalemi, 2026-09-24 (Silva08172026, mouse brain IPs): run_de.R re-quantified the DIA-NN
report with limpa and had no contaminant filter. DIA-NN's --cont-quant-exclude Cont_ only
shapes DIA-NN's own quantities, so all 121 Cont_ protein groups (bovine serum proteins from
the antibody prep) entered the DE and BH and came out as hits (bovine HBB +10.5 log2 in the
Kv2.1 IPs) -- while make_methods.py wrote that contaminants "were excluded from
quantification and normalisation".

Guards:
  * the R mirrors (contaminants.R) agree with fetch_fasta.py: the tag, and which sidecars
    are "legacy" (built before target-identical contaminants were removed);
  * end to end on a synthetic DIA-NN report (needs R with limpa/arrow/dplyr/tidyr/jsonlite;
    skips otherwise): a Cont_ group and a precursor SHARED with a Cont_ entry are gone from
    the matrix, for dpc and maxlfq; the counts are recorded; --keep-contaminants keeps them;
    a legacy sidecar is named; the emitted reproducibility script removes the same rows;
  * make_methods.py words both the DIA-NN step and the DE step from their records.
"""
import csv
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fetch_fasta as ff    # noqa: E402
import make_methods as mm   # noqa: E402

CONT_R = os.path.join(SCRIPTS, "contaminants.R")


def r_has(*pkgs):
    if not shutil.which("Rscript"):
        return False
    expr = "quit(status = if (all(vapply(c(%s), requireNamespace, logical(1), quietly = TRUE))) 0 else 1)" \
        % ", ".join(f'"{p}"' for p in pkgs)
    return subprocess.run(["Rscript", "-e", expr], capture_output=True).returncode == 0


def rscript(expr):
    p = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True)
    if p.returncode != 0:
        raise AssertionError(p.stderr[-1500:])
    return p.stdout


# A small DIA-NN-shaped report: 80 sample proteins x 4 precursors over 6 runs (A x3, B x3),
# low-intensity precursors missing more often (limpa's detection model needs that), a bovine
# serum contaminant strongly up in B -- the shape of the hit that went out -- and one sample
# precursor whose Protein.Ids also names a Cont_ entry (DIA-NN's rule excludes it too).
SYNTH_R = r'''
set.seed(1)
runs <- sprintf("run%02d", 1:6); grp <- rep(c("A", "B"), each = 3)
rows <- list()
add <- function(pg, ids, gene, prec, base, effB = 0) {
  for (r in seq_along(runs)) {
    lv <- base + rnorm(1, 0, 0.2) + (grp[r] == "B") * effB
    int <- 2^(lv + rnorm(1, 0, 0.3))
    if (runif(1) < plogis(-(lv - 13) * 2)) next
    rows[[length(rows) + 1]] <<- data.frame(Run = runs[r], Precursor.Id = prec,
      Protein.Group = pg, Protein.Ids = ids, Protein.Names = paste0(gene, "_X"), Genes = gene,
      Proteotypic = 1L, Precursor.Normalised = int, Precursor.Quantity = int,
      Q.Value = 0.001, Lib.Q.Value = 0.001, Lib.PG.Q.Value = 0.001, PG.Q.Value = 0.001,
      Global.Q.Value = 0.001, Global.PG.Q.Value = 0.001, stringsAsFactors = FALSE)
  }
}
for (i in 1:80) {
  pg <- sprintf("P%05d", i); base <- runif(1, 11, 20)
  for (k in 1:4) add(pg, pg, sprintf("G%d", i), sprintf("PEP%dK%d2", i, k), base + rnorm(1, 0, 1),
                     effB = if (i <= 5) 1.5 else 0)
}
for (k in 1:3) add("Cont_P02070", "Cont_P02070", "HBB", sprintf("HBBPEP%dK2", k), 18, effB = 4)
add("P00001", "P00001;Cont_P00761", "G1", "SHAREDPEPK2", 17)
d <- do.call(rbind, rows)
pgm <- aggregate(Precursor.Normalised ~ Protein.Group + Run, d, sum)
names(pgm)[3] <- "PG.MaxLFQ"
d <- merge(d, pgm, by = c("Protein.Group", "Run"))
arrow::write_parquet(d, file.path(OUT, "report.parquet"))
write.csv(data.frame(File.Name = runs, Group = grp), file.path(OUT, "conditions.csv"),
          row.names = FALSE)
'''

# What an old sidecar looks like: contaminants appended, no contaminant_target_rule
# (Silva08172026's own search.fasta.meta.json has exactly these keys).
LEGACY_SIDECAR = {"organism": "Mus musculus", "taxid": 10090, "n_contaminants_appended": 381,
                  "contaminant_set": "universal", "diann_cont_quant_exclude": "Cont_"}
CLEAN_SIDECAR = dict(LEGACY_SIDECAR, contaminant_target_rule=ff.CONTAMINANT_TARGET_RULE,
                     min_unique_peptides=ff.MIN_UNIQUE_PEPTIDES,
                     contaminants_dropped_as_target=[], contaminants_identical_to_target_kept=[])
# fetch_fasta.py before 2.8.0: the identity rule alone -- the rule, no min_unique_peptides.
IDENTITY_SIDECAR = dict(LEGACY_SIDECAR, contaminant_target_rule=ff.contaminant_target_rule(0),
                        contaminants_dropped_as_target=[], contaminants_identical_to_target_kept=[])


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


class RMirrorsAgreeWithFetchFasta(unittest.TestCase):
    def test_tag_is_fetch_fastas_tag(self):
        with open(CONT_R) as fh:
            m = re.search(r'^CONTAMINANT_TAG\s*<-\s*"([^"]+)"', fh.read(), re.M)
        self.assertIsNotNone(m, "CONTAMINANT_TAG not found in contaminants.R")
        self.assertEqual(m.group(1), ff.CONT_TAG)

    def test_constants_are_fetch_fastas(self):
        with open(CONT_R) as fh:
            r = fh.read()
        keep = re.search(r'^KEEP_TARGET_CONTAMINANTS_RULE\s*<-\s*"([^"]+)"', r, re.M)
        self.assertEqual(keep.group(1), ff.KEEP_TARGET_CONTAMINANTS_RULE)
        adv = re.search(r'^REBUILD_ADVICE\s*<-\s*paste0\(((?:\s*"[^"]*",?)+)\)', r, re.M)
        self.assertEqual("".join(re.findall(r'"([^"]*)"', adv.group(1))), ff.REBUILD_ADVICE)

    def test_sidecar_state_matches(self):
        """contaminants.R's sidecar_state() is fetch_fasta.sidecar_state(), case for case."""
        if not r_has("jsonlite"):
            self.skipTest("needs Rscript + jsonlite")
        cases = [LEGACY_SIDECAR, CLEAN_SIDECAR, IDENTITY_SIDECAR, {}, {"contaminant_set": "none"},
                 {"contaminant_set": None, "n_contaminants_appended": 0},
                 {"n_contaminants_already_present": 12},
                 {"contaminant_target_rule": None, "n_contaminants_appended": 5},
                 {"contaminant_set": "cell_culture"},
                 dict(IDENTITY_SIDECAR, min_unique_peptides=0),
                 dict(IDENTITY_SIDECAR, min_unique_peptides=None),
                 dict(IDENTITY_SIDECAR, min_unique_peptides=3),
                 dict(LEGACY_SIDECAR, contaminant_target_rule=ff.KEEP_TARGET_CONTAMINANTS_RULE)]
        self.assertEqual({ff.sidecar_state(c) for c in cases}, {"legacy", "identity_only", "current"})
        with tempfile.TemporaryDirectory() as tmp:
            paths = []
            for i, c in enumerate(cases):
                paths.append(os.path.join(tmp, f"{i}.json"))
                with open(paths[-1], "w") as fh:
                    json.dump(c, fh)
            out = rscript(f'source("{CONT_R}"); for (p in c({", ".join(repr(p) for p in paths)})) '
                          'cat(sidecar_state(jsonlite::fromJSON(p, simplifyVector = FALSE)), "\\n")'
                          .replace("'", '"'))
        self.assertEqual(out.split(), [ff.sidecar_state(c) for c in cases])

    def test_legacy_rule_matches(self):
        if not r_has("jsonlite"):
            self.skipTest("needs Rscript + jsonlite")
        cases = [LEGACY_SIDECAR, CLEAN_SIDECAR, IDENTITY_SIDECAR, {}, {"contaminant_set": "none"},
                 {"contaminant_set": None, "n_contaminants_appended": 0},
                 {"n_contaminants_already_present": 12},
                 {"contaminant_target_rule": None, "n_contaminants_appended": 5},
                 {"contaminant_set": "cell_culture"}]
        with tempfile.TemporaryDirectory() as tmp:
            paths = []
            for i, c in enumerate(cases):
                paths.append(os.path.join(tmp, f"{i}.json"))
                with open(paths[-1], "w") as fh:
                    json.dump(c, fh)
            out = rscript(f'source("{CONT_R}"); for (p in c({", ".join(repr(p) for p in paths)})) '
                          'cat(sidecar_is_legacy(jsonlite::fromJSON(p, simplifyVector = FALSE)), "\\n")'
                          .replace("'", '"'))
        self.assertEqual([x == "TRUE" for x in out.split()],
                         [ff._is_legacy_sidecar(c) for c in cases])

    def test_rule_is_any_accession_not_a_substring(self):
        if not shutil.which("Rscript"):
            self.skipTest("Rscript not available")
        ids = ["P1", "Cont_P2", "P3;Cont_P4", "Cont_P5;P6", "XCont_P7", "P8;YCont_P9", "NA"]
        out = rscript(f'source("{CONT_R}"); x <- c({", ".join(repr(i) for i in ids)}); '
                      'x[x == "NA"] <- NA; cat(is_contaminant(x))'.replace("'", '"'))
        self.assertEqual(out.split(), ["FALSE", "TRUE", "TRUE", "TRUE", "FALSE", "FALSE", "FALSE"])


@unittest.skipUnless(r_has("limpa", "limma", "arrow", "dplyr", "tidyr", "jsonlite"),
                     "needs R with limpa/limma/arrow/dplyr/tidyr/jsonlite")
class RunDeRemovesContaminants(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        rscript(f'OUT <- "{cls.tmp}"\n' + SYNTH_R)
        for name, meta in (("legacy.json", LEGACY_SIDECAR), ("clean.json", CLEAN_SIDECAR),
                           ("identity.json", IDENTITY_SIDECAR)):
            with open(os.path.join(cls.tmp, name), "w") as fh:
                json.dump(meta, fh)
        cls.runs = {}
        for key, args in (("dpc", ["--method", "dpc"]),
                          ("keep", ["--method", "dpc", "--keep-contaminants"]),
                          ("maxlfq", ["--method", "maxlfq"]),
                          ("legacy", ["--method", "dpc", "--fasta-meta",
                                      os.path.join(cls.tmp, "legacy.json")]),
                          ("clean", ["--method", "dpc", "--fasta-meta",
                                     os.path.join(cls.tmp, "clean.json")]),
                          ("identity", ["--method", "dpc", "--fasta-meta",
                                        os.path.join(cls.tmp, "identity.json")])):
            out = os.path.join(cls.tmp, key)
            p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"),
                                "--input", os.path.join(cls.tmp, "report.parquet"),
                                "--metadata", os.path.join(cls.tmp, "conditions.csv"),
                                "--outdir", out, *args],
                               capture_output=True, text=True, cwd=cls.tmp)
            cls.runs[key] = (p, out)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def out(self, key):
        p, out = self.runs[key]
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        return out

    def prov(self, key):
        with open(os.path.join(self.out(key), "de_provenance.json")) as fh:
            return json.load(fh)

    def methods(self, key):
        with open(os.path.join(self.out(key), "methods.txt")) as fh:
            return fh.read()

    def groups(self, key, name="Expression_Matrix.csv"):
        return {r["Protein.Group"] for r in read_csv(os.path.join(self.out(key), name))}

    def test_contaminant_group_is_gone_from_matrix_and_de(self):
        for key, de in (("dpc", "DE_dpc_B.A.csv"), ("maxlfq", "DE_maxlfq_B.A.csv")):
            with self.subTest(method=key):
                self.assertNotIn("Cont_P02070", self.groups(key))
                self.assertNotIn("Cont_P02070", self.groups(key, de))
                self.assertIn("P00001", self.groups(key))   # the sample protein stays

    def test_counts_are_recorded(self):
        for key in ("dpc", "maxlfq"):
            with self.subTest(method=key):
                c = self.prov(key)["contaminants"]
                self.assertEqual(c["policy"], "removed")
                self.assertTrue(c["removed"])
                self.assertEqual(c["id_column"], "Protein.Ids")
                self.assertEqual(c["n_precursors"], 4)        # 3 HBB + the shared one
                self.assertEqual(c["n_protein_groups"], 1)
                self.assertEqual(c["n_sample_groups_sharing"], 1)
                self.assertEqual(c["n_sample_groups_all_shared"], 0)
                self.assertTrue(any(f.startswith("contaminants removed: 4 precursors")
                                    for f in self.prov(key)["filters_applied"]))
                rem = {r["Protein.Group"]: r for r in
                       read_csv(os.path.join(self.out(key), c["removed_table"]))}
                self.assertEqual(rem["Cont_P02070"]["Removed.Entirely"], "TRUE")
                self.assertEqual(rem["P00001"]["Contaminant.Precursors"], "1")
                self.assertEqual(rem["P00001"]["Removed.Entirely"], "FALSE")
                self.assertIn("REMOVED before quantification", self.methods(key))

    def test_shared_precursor_is_removed_from_the_limpa_input(self):
        rds = os.path.join(self.out("dpc"), "DE-LIMP_session.rds")
        out = rscript(f's <- readRDS("{rds}"); cat(nrow(s$raw_data$E), '
                      '"SHAREDPEPK2" %in% rownames(s$raw_data$E), '
                      'any(grepl("^HBBPEP", rownames(s$raw_data$E))), '
                      '"Protein.Ids" %in% names(s$raw_data$genes))')
        n, shared, hbb, ids_col = out.split()
        self.assertEqual((shared, hbb, ids_col), ("FALSE", "FALSE", "FALSE"))
        self.assertEqual(int(n), self.prov("dpc")["contaminants"]["n_precursors_total"] - 4)

    def test_share_is_a_qc_table(self):
        rows = read_csv(os.path.join(self.out("dpc"), "QC_contaminant_share.csv"))
        self.assertEqual(len(rows), 6)
        pct = {r["Run"]: float(r["Contaminant.Pct"]) for r in rows}
        a = [pct[r] for r in ("run01", "run02", "run03")]
        b = [pct[r] for r in ("run04", "run05", "run06")]
        self.assertGreater(min(b), max(a))            # the contaminant is 16x up in B
        self.assertEqual({r["Group"] for r in rows}, {"A", "B"})
        # written with --keep-contaminants too: it describes the sample, not the filter
        self.assertTrue(os.path.exists(os.path.join(self.out("keep"), "QC_contaminant_share.csv")))

    def test_keep_contaminants_keeps_them_and_says_so(self):
        self.assertIn("Cont_P02070", self.groups("keep"))
        c = self.prov("keep")["contaminants"]
        self.assertEqual(c["policy"], "kept")
        self.assertFalse(c["removed"])
        self.assertIsNone(c.get("removed_table"))
        self.assertFalse(os.path.exists(os.path.join(self.out("keep"), "contaminants_removed.csv")))
        self.assertIn("KEPT (--keep-contaminants)", self.methods("keep"))
        self.assertNotIn("REMOVED", self.methods("keep"))

    def test_legacy_database_risk_is_named(self):
        c = self.prov("legacy")["contaminants"]
        self.assertTrue(c["database_checked"])
        self.assertTrue(c["database_risk"])
        self.assertIn("built before fetch_fasta.py removed", c["database_note"])
        self.assertIn("real Mus musculus proteins", c["database_note"])
        self.assertIn("CAUTION:", self.methods("legacy"))
        self.assertIn("contaminant filter:", self.runs["legacy"][0].stderr)   # R warning

    def test_identity_only_database_risk_names_the_near_identical_set(self):
        """A FASTA built by the identity rule alone (fetch_fasta.py before 2.8.0) kept bovine
        EEF1A1 / YWHAZ as Cont_ entries: for mouse, DIA-NN reports Eef1a1 / Ywhaz only as those,
        and removing the Cont_ groups takes them out of the DE -- a CAUTION, with the fix."""
        c = self.prov("identity")["contaminants"]
        self.assertTrue(c["database_checked"])
        self.assertTrue(c["database_risk"])
        for part in ("identity rule alone", "NEAR-identical to Mus musculus proteins",
                     "with the universal set, 10 of them: bovine EEF1A1, YWHAZ and TUBA1D",
                     "now missing from the DE", "audit_results.py --fasta-meta",
                     ff.REBUILD_ADVICE):
            self.assertIn(part, c["database_note"])
        self.assertIn("CAUTION:", self.methods("identity"))

    def test_clean_sidecar_carries_no_risk(self):
        c = self.prov("clean")["contaminants"]
        self.assertTrue(c["database_checked"])
        self.assertFalse(c["database_risk"])
        self.assertNotIn("CAUTION", self.methods("clean"))

    def test_without_a_sidecar_the_check_is_said_not_run(self):
        c = self.prov("dpc")["contaminants"]
        self.assertFalse(c["database_checked"])
        self.assertIsNone(c.get("database_risk"))
        self.assertIn("not run: no FASTA sidecar", c["database_note"])

    def test_reproducibility_script_removes_the_same_rows(self):
        for key, meth in (("dpc", "dpc"), ("maxlfq", "maxlfq")):
            with self.subTest(method=key):
                with open(os.path.join(self.out(key), "reproducibility_log.R")) as fh:
                    src = fh.read()
                self.assertIn("Remove contaminants", src)
                # every contaminant tag: the skill's Cont_ and FragPipe's contam_
                self.assertIn("grepl('(^|;)(Cont_|contam_)'", src)
                rerun = os.path.join(self.tmp, f"rerun_{key}")
                os.makedirs(rerun, exist_ok=True)
                shutil.copy(os.path.join(self.out(key), "reproducibility_log.R"), rerun)
                p = subprocess.run(["Rscript", "reproducibility_log.R"], cwd=rerun,
                                   capture_output=True, text=True)
                self.assertEqual(p.returncode, 0, p.stderr[-1500:])
                de = f"DE_{meth}_B.A.csv"
                orig = {r["Protein.Group"]: float(r["logFC"])
                        for r in read_csv(os.path.join(self.out(key), de))}
                again = {r["Protein.Group"]: float(r["logFC"])
                         for r in read_csv(os.path.join(rerun, "de_results_rerun", de))}
                self.assertEqual(set(orig), set(again))
                self.assertLess(max(abs(orig[k] - again[k]) for k in orig), 1e-9)

    def test_detection_matrix_matches_the_expression_matrix_and_the_qc_totals(self):
        for key in ("dpc", "maxlfq"):
            with self.subTest(method=key):
                out = self.out(key)
                with open(os.path.join(out, "Expression_Matrix.csv"), newline="") as fh:
                    em = list(csv.reader(fh))
                with open(os.path.join(out, "Detection_Matrix.csv"), newline="") as fh:
                    dm = list(csv.reader(fh))
                samples = [c for c in em[0] if c not in ("Protein.Group", "Genes", "Protein.Names")]
                self.assertEqual(dm[0], ["Protein.Group"] + samples)            # same columns
                self.assertEqual([r[0] for r in dm[1:]], [r[0] for r in em[1:]])   # same rows
                det = self.prov(key)["detection_matrix"]
                self.assertEqual(det["file"], "Detection_Matrix.csv")
                self.assertEqual(det["zero_means"], "inferred" if key == "dpc" else "missing")
                if key == "maxlfq":      # 0 exactly where the matrix is NA
                    ecol = {c: i for i, c in enumerate(em[0])}
                    for er, dr in zip(em[1:], dm[1:]):
                        for j, smp in enumerate(samples, 1):
                            self.assertEqual(dr[j] == "0", er[ecol[smp]] in ("", "NA"))
        rows = read_csv(os.path.join(self.out("dpc"), "Detection_Matrix.csv"))
        qc = {r["Sample"]: int(r["Detected"])
              for r in read_csv(os.path.join(self.out("dpc"), "QC_detected_vs_inferred.csv"))}
        for smp, n in qc.items():
            self.assertEqual(sum(1 for r in rows if r[smp] not in ("", "NA") and int(r[smp]) > 0), n)

    def test_keep_run_emits_no_filter(self):
        with open(os.path.join(self.out("keep"), "reproducibility_log.R")) as fh:
            self.assertNotIn("Remove contaminants", fh.read())


def de_prov(**contaminants):
    p = {"display_label": "DPC-Quant + limma (limpa)", "q_columns": ["Q.Value"],
         "q_cutoffs": [0.01], "design": "~ 0 + groups", "de_engine": "limpa::dpcDE",
         "adjp": 0.05, "logfc": 1, "logfc_role": "reference_line_only"}
    if contaminants:
        p["contaminants"] = contaminants
    return p


REMOVED = dict(policy="removed", removed=True, tag="Cont_", id_column="Protein.Ids",
               n_precursors=2196, n_protein_groups=121, n_sample_groups_sharing=47,
               n_sample_groups_all_shared=0)


class MethodsFollowTheRecord(unittest.TestCase):
    def test_removed(self):
        s = mm.de_paragraph(de_prov(**REMOVED))
        self.assertIn("2,196 precursors mapping to a Cont_-tagged contaminant entry", s)
        self.assertIn("taking out 121 contaminant protein groups", s)
        self.assertIn("47 lost some", s)
        self.assertNotIn("kept", s)

    def test_kept_never_says_removed_or_excluded(self):
        s = mm.de_paragraph(de_prov(policy="kept", removed=False, tag="Cont_",
                                    n_precursors=2196, n_protein_groups=121))
        self.assertIn("were kept in the differential-expression analysis (--keep-contaminants)", s)
        self.assertNotIn("removed", s)
        self.assertNotIn("excluded", s)

    def test_a_record_without_the_step_is_tagged_not_assumed(self):
        s = mm.de_paragraph(de_prov())
        self.assertIn(f"Contaminant handling in the differential-expression step: "
                      f"{mm.NOT_RECORDED}", s)
        self.assertNotIn("excluded", s)

    def test_diann_step_comes_from_the_parameters_not_the_sidecar(self):
        self.assertIn("excluded from normalisation",
                      mm.diann_contaminant_sentence({"engine": "diann", "cont_quant_exclude":
                                                     {"value": "Cont_", "source": "p.cfg"}}))
        self.assertIn("was not set", mm.diann_contaminant_sentence(
            {"engine": "diann", "cont_quant_exclude": {"value": None, "source": "p.cfg"}}))
        self.assertIn(mm.NOT_RECORDED, mm.diann_contaminant_sentence({"engine": None}))
        self.assertEqual(mm.diann_contaminant_sentence({"engine": "sage"}), "")

    def test_search_record_reads_the_flag(self):
        with tempfile.TemporaryDirectory() as tmp:
            cfg = os.path.join(tmp, "params.cfg")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n--cont-quant-exclude Cont_\n--cut K*,R*\n")
            self.assertEqual(mm.search_record(params=cfg)["cont_quant_exclude"]["value"], "Cont_")
            with open(cfg, "w") as fh:
                fh.write("--qvalue 0.01\n--cut K*,R*\n")
            self.assertIsNone(mm.search_record(params=cfg)["cont_quant_exclude"]["value"])

    def test_methods_md_no_longer_claims_exclusion_from_the_sidecar(self):
        """The Silva08172026 shape: the sidecar recommends --cont-quant-exclude, the DE kept
        contaminants. Nothing in the Methods may say they were excluded."""
        with tempfile.TemporaryDirectory() as tmp:
            meta = os.path.join(tmp, "search.fasta.meta.json")
            with open(meta, "w") as fh:
                json.dump(dict(LEGACY_SIDECAR, proteome="UP000000589", uniprot_release="2026_03",
                               content_used="one_per_gene", n_proteome=21860), fh)
            de = os.path.join(tmp, "de")
            os.makedirs(de)
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump(de_prov(policy="kept", removed=False, tag="Cont_",
                                  n_precursors=2196, n_protein_groups=121), fh)
            out = os.path.join(tmp, "methods.md")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                                "--raw", os.path.join(tmp, "Ex_sample.raw"), "--fasta-meta", meta,
                                "--de-dir", de, "--out", out], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(out) as fh:
                text = fh.read()
        self.assertNotIn("excluded from quantification and normalisation", text)
        self.assertIn("was appended.", text)
        self.assertIn("--keep-contaminants", text)
        self.assertIn("Whether DIA-NN's --cont-quant-exclude was set", text)

    def test_caveat_when_the_filter_probably_removed_real_proteins(self):
        with tempfile.TemporaryDirectory() as tmp:
            de = os.path.join(tmp, "de")
            os.makedirs(de)
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump(de_prov(**REMOVED, database_checked=True, database_risk=True,
                                  database_note="the search database was built before ..."), fh)
            out = os.path.join(tmp, "methods.md")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                                "--raw", os.path.join(tmp, "Ex_sample.raw"), "--de-dir", de,
                                "--out", out], capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(out) as fh:
                text = fh.read()
        self.assertIn("> Contaminant filter caveat (resolve before publication): the search "
                      "database was built before", text)


class AuditKeepsTheEvidence(unittest.TestCase):
    """The filter takes Cont_ groups out of Expression_Matrix.csv, so audit_results.py's
    "real proteins present only as Cont_" evidence now also comes from the removed table."""

    def test_removed_real_protein_is_still_named(self):
        with tempfile.TemporaryDirectory() as tmp:
            de = os.path.join(tmp, "de")
            os.makedirs(de)
            with open(os.path.join(de, "Expression_Matrix.csv"), "w") as fh:
                fh.write("Protein.Group,Genes,Protein.Names,S1,S2\nP04406,GAPDH,G3P_HUMAN,18,18\n")
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump(de_prov(**dict(REMOVED, removed_table="contaminants_removed.csv")), fh)
            with open(os.path.join(de, "contaminants_removed.csv"), "w") as fh:
                fh.write("Protein.Group,Genes,Contaminant.Group,Precursors,"
                         "Contaminant.Precursors,Removed.Entirely\n"
                         "Cont_P60712,ACTB,TRUE,9,9,TRUE\n")
            kept = [{"cont_acc": "Cont_P60712", "gene": "ACTB", "target_acc": "P60709"}]
            meta = os.path.join(tmp, "m.json")
            with open(meta, "w") as fh:
                json.dump({"organism": "Homo sapiens", "contaminant_target_rule": "x",
                           "contaminants_identical_to_target_kept": kept}, fh)
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "audit_results.py"),
                                "--out", os.path.join(tmp, "AUDIT.md"), "--de-dir", de,
                                "--fasta-meta", meta], capture_output=True, text=True, cwd=tmp)
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(os.path.join(tmp, "AUDIT.json")) as fh:
                f = next(x for x in json.load(fh)["findings"]
                         if x["check"] == "contaminant_overlap")
        self.assertEqual(f["detail"]["removed_by_de_filter"], ["ACTB (Cont_P60712)"])
        self.assertIn("run_de.R's contaminant filter then removed every Cont_ group", f["message"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
