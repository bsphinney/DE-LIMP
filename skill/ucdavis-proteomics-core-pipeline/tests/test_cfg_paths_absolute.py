#!/usr/bin/env python3
"""
A relative --cfg can be read from the wrong file (SKILL_OPEN_DEFECTS, sage-review of 6d49ec3).

DIA-NN logs `--cfg <file>` as it was given, and does not log its working directory. So
make_methods' reader (the one reader of DIA-NN's --cont-quant-exclude, which the DE descriptors,
the Methods and the run record ask) tried the log's folder, then its parent -- and with a cfg of
that name in both, read the log's folder's even when DIA-NN ran in the parent. Now:

  * every place the skill takes a DIA-NN cfg path makes it absolute at entry and refuses a
    missing one: run_search.run_diann (which writes `--cfg` into the job), diann_parallel.py
    (which records the path and hands it to later jobs) -- run_search.main() already did;
  * the reader never picks between two different files: a relative --cfg found in more than one
    candidate folder with different contents is NOT RECORDED, naming both, and the reason
    reaches contaminants.R instead of "no DIA-NN log beside the report".

Hermetic: fixture cfgs and logs, a fake DIA-NN, no SLURM, no network.
"""
import glob
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import make_methods  # noqa: E402
import run_search  # noqa: E402
from test_cfg_reader_quoting import _estimate, _fake_diann  # noqa: E402
from test_run_de_contaminants import r_has  # noqa: E402

GENERATOR = os.path.join(SCRIPTS, "diann_parallel.py")
PINNED = "--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n--window 7\n--cont-quant-exclude Cont_\n"


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    return path


def diann_log(path, line):
    write(path, "\nDIA-NN 2.7.0 Academia  (Data-Independent Acquisition by Neural Networks)\n"
                "Compiled on Sep 15 2026 15:30:01\nLogical CPU cores: 64\n" + line + "\n\n")


class TheReaderNeverPicksBetweenTwoFiles(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.out = os.path.join(self.d, "run", "search")
        self.report = os.path.join(self.out, "report.parquet")
        write(self.report, "")
        diann_log(os.path.join(self.out, "report.log.txt"),
                  f"/x/diann-linux --cfg params.cfg --f /d/a.raw --out {self.report}")

    def test_the_same_name_in_both_folders_with_different_contents_is_not_recorded(self):
        here = write(os.path.join(self.out, "params.cfg"), "--cont-quant-exclude Cont_\n")
        parent = write(os.path.join(self.d, "run", "params.cfg"), "--qvalue 0.01\n")
        rec, why = make_methods.diann_cont_quant_exclude_why(self.report)
        self.assertIsNone(rec)
        self.assertIsNone(make_methods.diann_cont_quant_exclude(self.report))
        self.assertIn("relative path", why)
        self.assertIn(here, why)
        self.assertIn(parent, why)
        self.assertIn("which one DIA-NN read is not recorded", why)

    def test_the_same_file_contents_in_both_is_read(self):
        write(os.path.join(self.out, "params.cfg"), "--cont-quant-exclude Cont_\n")
        write(os.path.join(self.d, "run", "params.cfg"), "--cont-quant-exclude Cont_\n")
        rec = make_methods.diann_cont_quant_exclude(self.report)
        self.assertEqual(rec["value"], "Cont_")

    def test_one_candidate_is_read_as_before(self):
        write(os.path.join(self.d, "run", "params.cfg"), "--cont-quant-exclude Cont_\n")
        self.assertEqual(make_methods.diann_cont_quant_exclude(self.report)["value"], "Cont_")

    def test_the_parameters_file_the_search_ran_with_still_answers(self):
        """An ambiguous --cfg is unread, not absent: the search's own record is asked next."""
        write(os.path.join(self.out, "params.cfg"), "--cont-quant-exclude Cont_\n")
        write(os.path.join(self.d, "run", "params.cfg"), "--qvalue 0.01\n")
        pf = write(os.path.join(self.d, "wf", "params.resolved.cfg"),
                   "--cont-quant-exclude Cont_\n")
        with open(os.path.join(self.out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "diann", "resolved_params_file": pf}, fh)
        self.assertEqual(make_methods.diann_cont_quant_exclude(self.report),
                         {"value": "Cont_", "source": "params.resolved.cfg"})

    def test_the_reason_reaches_the_methods_and_the_cli(self):
        write(os.path.join(self.out, "params.cfg"), "--cont-quant-exclude Cont_\n")
        write(os.path.join(self.d, "run", "params.cfg"), "--qvalue 0.01\n")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_methods.py"),
                            "cont-quant-exclude", self.report],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        j = json.loads(p.stdout)
        self.assertIs(j["recorded"], False)
        self.assertIn("exist and differ", j["why"])
        sentence = make_methods.diann_contaminant_sentence(
            {"engine": "diann", "cont_quant_exclude": None,
             "cont_quant_exclude_why": j["why"]})
        self.assertIn("NOT RECORDED", sentence.upper())
        self.assertIn("exist and differ", sentence)

    def test_an_unreadable_cfg_says_so_not_no_log(self):
        diann_log(os.path.join(self.out, "report.log.txt"),
                  f"/x/diann-linux --cfg /gone/params.cfg --out {self.report}")
        rec, why = make_methods.diann_cont_quant_exclude_why(self.report)
        self.assertIsNone(rec)
        self.assertIn("the log names --cfg /gone/params.cfg, which could not be read", why)
        self.assertNotIn("no DIA-NN log", why)

    @unittest.skipUnless(r_has("jsonlite"), "needs R with jsonlite")
    def test_contaminants_r_quotes_the_reason(self):
        write(os.path.join(self.out, "params.cfg"), "--cont-quant-exclude Cont_\n")
        write(os.path.join(self.d, "run", "params.cfg"), "--qvalue 0.01\n")
        expr = (f'.script_dir <- "{SCRIPTS}"; source(file.path(.script_dir, "contaminants.R")); '
                f'k <- diann_kept_quant(diann_cont_quant_exclude("{self.report}"), FALSE); '
                f'cat(k$text)')
        p = subprocess.run(["Rscript", "-e", expr], capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        self.assertIn("[not recorded -- confirm]", p.stdout)
        self.assertIn("exist and differ", p.stdout)
        self.assertNotIn("no DIA-NN log beside the report", p.stdout)


class EveryEntryMakesTheCfgAbsolute(unittest.TestCase):
    def setUp(self):
        self.d = os.path.realpath(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.d, True)
        self.cwd = os.getcwd()
        self.addCleanup(os.chdir, self.cwd)

    def inputs(self):
        raws = []
        for i in range(6):
            p = os.path.join(self.d, f"f{i}.d")
            os.makedirs(p, exist_ok=True)
            raws.append(p)
        return raws, write(os.path.join(self.d, "db.fasta"), ">sp|P1|X\nPEPTIDER\n")

    def test_the_chain_generator_records_an_absolute_cfg(self):
        write(os.path.join(self.d, "wf", "params.cfg"), PINNED)
        raws, fasta = self.inputs()
        p = subprocess.run([sys.executable, GENERATOR, "--diann", "/bin/true", "--raw", *raws,
                            "--fasta", fasta, "--out", "out", "--cfg", "wf/params.cfg"],
                           capture_output=True, text=True, cwd=self.d, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        resolved = json.loads(p.stdout)["resolved_params"]["file"]
        self.assertEqual(resolved, os.path.join(self.d, "wf", "params.cfg"))

    def test_the_chain_generator_refuses_a_missing_cfg_naming_what_was_given(self):
        raws, fasta = self.inputs()
        p = subprocess.run([sys.executable, GENERATOR, "--diann", "/bin/true", "--raw", *raws,
                            "--fasta", fasta, "--out", "out", "--cfg", "wf/nope.cfg"],
                           capture_output=True, text=True, cwd=self.d, timeout=120)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn(f"cfg not found: {os.path.join(self.d, 'wf', 'nope.cfg')} "
                      "(given as wf/nope.cfg)", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.d, "out", "submit.sh")))

    def test_the_single_shot_job_passes_an_absolute_cfg(self):
        cfg = _estimate(self.d, "astral", "--acquisition", "DIA", "--instrument",
                        "Orbitrap Astral")
        fake = _fake_diann(self.d)
        raws, fasta = self.inputs()
        os.chdir(self.d)
        run_search.run_diann(fake, os.path.basename(cfg), raws, fasta,
                             os.path.join(self.d, "out"), 8, os.path.join(self.d, "job.sh"),
                             acquisition="DIA")
        cfg_args = []
        for script in glob.glob(os.path.join(self.d, "job*.sh")):
            with open(script) as fh:
                words = shlex.split(fh.read(), comments=True)
            cfg_args += [words[i + 1] for i, w in enumerate(words[:-1]) if w == "--cfg"]
        self.assertTrue(cfg_args, "no --cfg on any job's command line")
        self.assertEqual(set(cfg_args), {cfg})

    def test_the_single_shot_route_refuses_a_missing_cfg_before_writing_anything(self):
        raws, fasta = self.inputs()
        os.chdir(self.d)
        with self.assertRaises(SystemExit) as cm:
            run_search.run_diann("/bin/true", "wf/nope.cfg", raws, fasta,
                                 os.path.join(self.d, "out"), 8, os.path.join(self.d, "job.sh"),
                                 acquisition="DIA")
        self.assertIn(f"cfg not found: {os.path.join(self.d, 'wf', 'nope.cfg')}",
                      str(cm.exception))
        self.assertFalse(os.path.exists(os.path.join(self.d, "out")))


if __name__ == "__main__":
    unittest.main()
