#!/usr/bin/env python3
"""No sidecar must not mean no caveat.

Release review 2026-09-28: with no fetch_fasta.py sidecar, run_de.R's contaminant filter said only
"not run" -- so a search on the Core's superseded MRS human FASTA (the Siegel entries staged
2026-09-25) removed ACTB, EEF1A1, KRT8 and ~150 other real proteins from the DE with no caveat.
The check now reads the FASTA the search itself names (its provenance, or --fasta in its log),
re-checks it with fetch_fasta.py's rule, and -- when that file cannot be read here -- recognises a
known superseded database by md5 or name (superseded_databases.json). fetch_fasta.py defines it
once (database_without_sidecar, `check-db`); contaminants.R only words the answer.

Also here: run_search.py warns when the search digests differently from the digest the FASTA's
contaminants were judged on (the sidecar's contaminant_digest).
"""
import contextlib
import io
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

import fetch_fasta as ff    # noqa: E402
import run_search as rs     # noqa: E402
import estimate_params as ep  # noqa: E402

LEGACY_NAME = "UP000005640_9606_plus_universal_contam.fasta"
LEGACY_MD5 = "8de1d9bd0a052b175f88f66f82500d92"      # verified on HIVE, 2026-09-25

ACTB = "MDDDIAALVVDNGSGMCKAGFAGDDAPRAVFPSIVGRPRHQGVMVGMGQKDSYVGDEAQSKRGILTLKYPIEHGIVTNWDDMEK"
GAPDH = ("MGKVKVGVNGFGRIGRLVTRAAFNSGKVDIVAINDPFIDLNYMVYMFQYDSTHGKFHGTVKAENGKLVINGNPITIFQERDPSK"
         "IKWGDAGAEYVVESTGVFTTMEKAGAHLQGGAKRVIISAPSADAPMFVMGVNHEKYDNSLKIISNASCTTNCLAPLAKVIHDNF")
# Enough human entries that OX=9606 holds the >= 95% majority fran_deposit.organism_from_headers
# needs even counting the Cont_ entries (it does count `>sp|Cont_...` ones -- reported upstream).
FILLER = "".join(f">sp|Q{i:05d}|F{i}_HUMAN Filler {i} OS=Homo sapiens OX=9606 GN=F{i} PE=1 SV=1\n"
                 + "".join("ACDEFGHILMNPQSTVWY"[(i * 7 + j) % 18] for j in range(15)) + "K\n"
                 for i in range(40))
HUMAN = (f">sp|P60709|ACTB_HUMAN Actin OS=Homo sapiens OX=9606 GN=ACTB PE=1 SV=1\n{ACTB}\n"
         f">sp|P04406|G3P_HUMAN GAPDH OS=Homo sapiens OX=9606 GN=GAPDH PE=1 SV=3\n{GAPDH}\n"
         + FILLER)
TRYPSIN = ">sp|Cont_P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1\nFPTDDDDKIVGGYTCAANSIPYQVSLNSG\n"
ACTB_BOVIN = f">sp|Cont_P60712|ACTB_BOVIN Actin OS=Bos taurus OX=9913 GN=ACTB PE=1 SV=1\n{ACTB}\n"


def search_dir(root, fasta_arg, name="search"):
    """A search folder whose DIA-NN log names `fasta_arg` -- what fran_deposit.search_fastas reads."""
    d = os.path.join(root, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "report.log.txt"), "w") as fh:
        fh.write(f"DIA-NN 2.7.0 Academia\ndiann-linux --f a.raw --fasta {fasta_arg} --out report.parquet\n")
    return d


def write(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    return path


class DatabaseWithoutSidecar(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name

    def tearDown(self):
        self._td.cleanup()

    def test_named_fasta_is_rechecked(self):
        fa = write(os.path.join(self.root, "human_plus.fasta"), HUMAN + ACTB_BOVIN + TRYPSIN)
        r = ff.database_without_sidecar(search_dir(self.root, fa))
        self.assertTrue(r["checked"])
        self.assertTrue(r["risk"])
        self.assertEqual(r["organism"], "Homo sapiens")          # from the FASTA's own OX=/OS=
        self.assertEqual((r["n_lost"], r["genes"]), (1, ["ACTB"]))
        self.assertIn(f"re-checked {fa}: 1 of its Cont_ entries", r["why"])

    def test_clean_fasta_is_checked_and_passes(self):
        fa = write(os.path.join(self.root, "clean.fasta"), HUMAN + TRYPSIN)
        r = ff.database_without_sidecar(search_dir(self.root, fa))
        self.assertEqual((r["checked"], r["risk"], r["genes"]), (True, False, []))
        self.assertIn("none of its Cont_ entries is a target protein", r["why"])

    def test_unreadable_superseded_mrs_file_is_recognised_by_name(self):
        """A Core search run on Windows logs Y:\\MRS\\...; HIVE cannot read that path."""
        r = ff.database_without_sidecar(search_dir(self.root, "Y:\\MRS\\" + LEGACY_NAME))
        self.assertTrue(r["checked"])
        self.assertTrue(r["risk"])
        self.assertEqual((r["organism"], r["n_lost"]), ("Homo sapiens", 161))
        for g in ("ACTB", "EEF1A1", "YWHAZ", "TUBB", "KRT8"):
            self.assertIn(g, r["genes"])
        self.assertLess(r["genes"].index("ACTB"), r["genes"].index("KRT8"))   # non-keratins first
        self.assertIn("by its name it is the Core's Sep-2025 MRS human", r["why"])

    def test_a_readable_file_is_rechecked_even_with_the_legacy_name(self):
        """The name is only a fallback: a readable file answers for itself."""
        fa = write(os.path.join(self.root, LEGACY_NAME), HUMAN + TRYPSIN)
        r = ff.database_without_sidecar(search_dir(self.root, fa))
        self.assertEqual((r["risk"], r["n_lost"]), (False, 0))
        self.assertNotIn("Sep-2025", r["why"])

    def test_unknown_unreadable_fasta_is_not_guessed(self):
        r = ff.database_without_sidecar(search_dir(self.root, "/nowhere/other.fasta"))
        self.assertFalse(r["checked"])
        self.assertIsNone(r["risk"])
        self.assertIn("/nowhere/other.fasta", r["why"])

    def test_a_search_that_names_no_fasta(self):
        d = os.path.join(self.root, "bare")
        os.makedirs(d)
        r = ff.database_without_sidecar(d)
        self.assertFalse(r["checked"])
        self.assertIn("names no FASTA", r["why"])

    def test_check_db_cli_prints_the_answer(self):
        fa = write(os.path.join(self.root, "human_plus.fasta"), HUMAN + ACTB_BOVIN)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "fetch_fasta.py"), "check-db",
                            "--search-dir", search_dir(self.root, fa)],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(json.loads(p.stdout)["genes"], ["ACTB"])


class SupersededData(unittest.TestCase):
    def test_the_legacy_record_is_recognised_by_md5_and_by_name(self):
        by_md5 = ff.superseded_database(md5=LEGACY_MD5)
        self.assertEqual(by_md5["basename"], LEGACY_NAME)
        self.assertEqual(ff.superseded_database(path="/quobyte/proteomics-grp/MRS/" + LEGACY_NAME)["md5"],
                         LEGACY_MD5)
        self.assertEqual(ff.superseded_database(path="Y:\\MRS\\" + LEGACY_NAME)["md5"], LEGACY_MD5)
        # A read file with another md5 is NOT the known database, whatever its name.
        self.assertIsNone(ff.superseded_database(path=LEGACY_NAME, md5="0" * 32))

    def test_the_record_is_internally_consistent(self):
        with open(ff.SUPERSEDED_FILE) as fh:
            d = json.load(fh)["databases"][0]
        self.assertEqual(sum(d["by_reason"].values()), d["n_lost"])
        self.assertEqual(len(d["genes"]), d["n_lost"])
        self.assertEqual(d["by_reason"], {"identical": 152, "substring": 1, "shared_peptides": 8})


def r_ok():
    if not shutil.which("Rscript"):
        return False
    return subprocess.run(["Rscript", "-e", 'quit(status = !requireNamespace("jsonlite", quietly = TRUE))'],
                          capture_output=True).returncode == 0


@unittest.skipUnless(r_ok(), "needs Rscript + jsonlite")
class RWordsTheAnswer(unittest.TestCase):
    """contaminants.R asks fetch_fasta.py (check-db) and words the answer; it keeps no copy."""

    def risk(self, search_dir_expr):
        expr = (f'.script_dir <- "{SCRIPTS}"; source(file.path(.script_dir, "contaminants.R")); '
                f'r <- contaminant_database_risk(NULL, 140, search_dir = {search_dir_expr}); '
                'cat(jsonlite::toJSON(list(checked = r$checked, risk = r$risk, note = r$note), '
                'auto_unbox = TRUE, na = "null"))')
        p = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        return json.loads(p.stdout)

    def test_rechecked_database_risk_names_the_proteins(self):
        with tempfile.TemporaryDirectory() as root:
            fa = write(os.path.join(root, "human_plus.fasta"), HUMAN + ACTB_BOVIN + TRYPSIN)
            r = self.risk(f'"{search_dir(root, fa)}"')
        self.assertTrue(r["checked"])
        self.assertTrue(r["risk"])
        for part in ("no fetch_fasta.py sidecar", f"re-checked {fa}", "-- ACTB.",
                     "140 Cont_ protein groups removed here are probably real Homo sapiens proteins",
                     ff.REBUILD_ADVICE):
            self.assertIn(part, r["note"])

    def test_superseded_mrs_file_is_named_even_when_unreadable(self):
        with tempfile.TemporaryDirectory() as root:
            r = self.risk(f'"{search_dir(root, "/quobyte/proteomics-grp/MRS/" + LEGACY_NAME)}"')
        self.assertTrue(r["risk"])
        self.assertIn("Sep-2025 MRS human", r["note"])
        self.assertIn("ACTB", r["note"])

    def test_clean_database_passes(self):
        with tempfile.TemporaryDirectory() as root:
            fa = write(os.path.join(root, "clean.fasta"), HUMAN + TRYPSIN)
            r = self.risk(f'"{search_dir(root, fa)}"')
        self.assertEqual((r["checked"], r["risk"]), (True, False))
        self.assertTrue(r["note"].startswith("passed: "))

    def test_without_a_search_folder_it_is_said_not_run(self):
        r = self.risk("NULL")
        self.assertFalse(r["checked"])
        self.assertIsNone(r["risk"])
        self.assertIn("not run: no FASTA sidecar", r["note"])


def r_has_de():
    if not shutil.which("Rscript"):
        return False
    expr = ('quit(status = if (all(vapply(c("limpa", "arrow", "dplyr", "tidyr", "jsonlite"), '
            'requireNamespace, logical(1), quietly = TRUE))) 0 else 1)')
    return subprocess.run(["Rscript", "-e", expr], capture_output=True).returncode == 0


@unittest.skipUnless(r_has_de(), "needs R with limpa/arrow/dplyr/tidyr/jsonlite")
class RunDeWithoutSidecar(unittest.TestCase):
    """End to end: run_de.R with no --fasta-meta, beside a report whose log names a database
    holding a target-identical Cont_ entry, writes a CAUTION -- it used to say "not run"."""

    def test_run_de_cautions_from_the_searchs_own_fasta(self):
        from test_run_de_contaminants import SYNTH_R, rscript
        with tempfile.TemporaryDirectory() as tmp:
            rscript(f'OUT <- "{tmp}"\n' + SYNTH_R)
            fa = write(os.path.join(tmp, "human_plus.fasta"), HUMAN + ACTB_BOVIN + TRYPSIN)
            search_dir(os.path.dirname(tmp), fa, name=os.path.basename(tmp))   # log beside the report
            out = os.path.join(tmp, "de")
            p = subprocess.run(["Rscript", os.path.join(SCRIPTS, "run_de.R"),
                                "--input", os.path.join(tmp, "report.parquet"),
                                "--metadata", os.path.join(tmp, "conditions.csv"),
                                "--outdir", out, "--method", "dpc"],
                               capture_output=True, text=True, cwd=tmp)
            self.assertEqual(p.returncode, 0, p.stderr[-2000:])
            with open(os.path.join(out, "de_provenance.json")) as fh:
                c = json.load(fh)["contaminants"]
            with open(os.path.join(out, "methods.txt")) as fh:
                methods = fh.read()
        self.assertTrue(c["database_checked"])
        self.assertTrue(c["database_risk"])
        self.assertIn("no fetch_fasta.py sidecar", c["database_note"])
        self.assertIn("ACTB", c["database_note"])
        self.assertIn("CAUTION:", methods)


class SearchDigestWarning(unittest.TestCase):
    """run_search.py: the search's digest vs the one the FASTA's contaminants were judged on."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.root = self._td.name
        self.fasta = write(os.path.join(self.root, "db.fasta"), HUMAN + TRYPSIN)
        self.meta = {"contaminant_digest": dict(ff.DIANN_DIGEST),
                     "min_unique_peptides": ff.MIN_UNIQUE_PEPTIDES}

    def tearDown(self):
        self._td.cleanup()

    def sidecar(self, meta):
        write(self.fasta + ".meta.json", json.dumps(meta))

    def diann_cfg(self, missed=1):
        return write(os.path.join(self.root, "diann.cfg"),
                     f"--qvalue 0.01\n--cut K*,R*\n--missed-cleavages {missed}\n--min-pep-len 7\n"
                     "--max-pep-len 30\n--met-excision\n--unimod4\n")

    def warn(self, engine, params):
        with contextlib.redirect_stderr(io.StringIO()):
            return rs.warn_contaminant_digest(engine, params, self.fasta)

    def test_matching_diann_cfg_is_silent(self):
        self.sidecar(self.meta)
        self.assertIsNone(self.warn("diann", self.diann_cfg()))

    def test_overridden_missed_cleavages_warns(self):
        self.sidecar(self.meta)
        msg = self.warn("diann", self.diann_cfg(missed=2))
        self.assertIn("missed_cleavages 2 (FASTA: 1)", msg)
        self.assertIn("contaminants_dropped_as_target", msg)

    def sage_cfg(self, **change):
        enz = dict(ep.SAGE_ENZYME, **change)       # what estimate_params.build_sage writes
        return write(os.path.join(self.root, "sage.json"), json.dumps({"database": {"enzyme": enz}}))

    def test_sage_default_enzyme_is_one_info_line_not_a_warning(self):
        """Every Sage search differs from the DIA-NN digest by Sage's own default (2 missed
        cleavages, no cleavage before P). That is KNOWN: a WARNING on every Sage search would
        teach people to ignore the warning, so it is one INFO line that says why."""
        self.sidecar(self.meta)
        msg = self.warn("sage", self.sage_cfg())
        self.assertTrue(msg.startswith("[run_search] INFO: Sage's default digest"), msg)
        self.assertIn("cut 'K*,R*,!*P' (FASTA: 'K*,R*')", msg)
        self.assertIn("missed_cleavages 2 (FASTA: 1)", msg)
        self.assertIn("expected for every Sage search", msg)
        self.assertNotIn("WARNING", msg)

    def test_sage_override_beyond_its_default_warns(self):
        """--overrides changing Sage's missed cleavages or length is NOT the known default."""
        self.sidecar(self.meta)
        msg = self.warn("sage", self.sage_cfg(missed_cleavages=3, min_len=6))
        self.assertTrue(msg.startswith("[run_search] WARNING:"), msg)
        self.assertIn("missed_cleavages 3 (FASTA: 1)", msg)
        self.assertIn("min_pep_len 6 (FASTA: 7)", msg)
        self.assertIn("(besides Sage's known default: cut 'K*,R*,!*P' (FASTA: 'K*,R*'))", msg)

    def test_build_sage_writes_the_one_sage_enzyme(self):
        text, _rationale = ep.build_sage("DDA", "orbitrap_generic", "", {})
        self.assertEqual(json.loads(text)["database"]["enzyme"], ep.SAGE_ENZYME)

    def test_no_peptide_rule_or_no_sidecar_is_silent(self):
        self.assertIsNone(self.warn("diann", self.diann_cfg(missed=2)))      # no sidecar
        self.sidecar({"contaminant_digest": dict(ff.DIANN_DIGEST), "min_unique_peptides": 0})
        self.assertIsNone(self.warn("diann", self.diann_cfg(missed=2)))      # identity rule only


if __name__ == "__main__":
    unittest.main()
