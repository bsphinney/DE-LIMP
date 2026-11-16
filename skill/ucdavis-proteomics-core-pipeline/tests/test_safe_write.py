#!/usr/bin/env python3
"""Job scripts and the files jobs read are never written through a planted link.

Security review of 2.12 (reproduced on make_podcast.py): Core session folders are group-writable
(umask 0002, no sticky bit), and the generators wrote with open(path, "w") + os.chmod(path), both
of which follow a symlink -- a group member who planted `step2_firstpass.sbatch -> ~victim/.bashrc`
had the victim's file overwritten with a job script and made executable. safe_write.write_file()
creates a new file beside the target (O_CREAT|O_EXCL|O_NOFOLLOW), sets its mode, and os.replace()s
it over the link. Every generator writes its scripts, lists and cfgs through it.

Also: run_search.py's `#SBATCH --output=` is double-quoted, because sbatch refuses an unquoted
path with a space ("Invalid directive found in batch script"; HIVE, 2026-10-15) and accepts and
writes to a double-quoted one (a test job). diatracer_parallel.py's -o/-e likewise.
"""
import json
import os
import stat
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import safe_write  # noqa: E402
from test_single_shot_sbatch import HIGH, LIBFREE_CFG, _Workspace, _read, _write  # noqa: E402

VICTIM_TEXT = "the victim's own file\n"


def plant(link, d):
    """A victim file elsewhere and a symlink to it at `link`; returns the victim's path."""
    os.makedirs(d, exist_ok=True)
    victim = _write(os.path.join(d, "victim_" + os.path.basename(link)), VICTIM_TEXT)
    os.makedirs(os.path.dirname(link), exist_ok=True)
    os.symlink(victim, link)
    return victim


class WriteFile(unittest.TestCase):
    def test_a_planted_link_is_replaced_never_followed(self):
        with tempfile.TemporaryDirectory() as d:
            target = os.path.join(d, "job.sh")
            victim = plant(target, d)
            safe_write.write_file(target, "#!/bin/bash\necho hi\n", 0o755)
            self.assertEqual(_read(victim), VICTIM_TEXT)
            self.assertFalse(os.path.islink(target))
            self.assertEqual(_read(target), "#!/bin/bash\necho hi\n")
            self.assertEqual(stat.S_IMODE(os.stat(target).st_mode), 0o755)
            self.assertEqual(stat.S_IMODE(os.stat(victim).st_mode) & 0o111, 0,
                             "the victim was made executable")
            self.assertEqual(sorted(os.listdir(d)), sorted(["job.sh", os.path.basename(victim)]),
                             "a temporary file was left behind")

    def test_mode_is_exact_whatever_the_umask(self):
        with tempfile.TemporaryDirectory() as d:
            old = os.umask(0o002)
            try:
                p = safe_write.write_file(os.path.join(d, "list.txt"), "a\n")
            finally:
                os.umask(old)
            self.assertEqual(stat.S_IMODE(os.stat(p).st_mode), 0o644)

    def test_a_failure_leaves_nothing_and_the_old_file_stands(self):
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "adir"))
            with self.assertRaises(OSError):
                safe_write.write_file(os.path.join(d, "adir"), "x")
            self.assertEqual(os.listdir(d), ["adir"])
            p = _write(os.path.join(d, "f.json"), "old")
            with self.assertRaises(ValueError):
                with safe_write.open_text(p) as fh:
                    fh.write("half")
                    raise ValueError("boom")
            self.assertEqual(_read(p), "old")
            with safe_write.open_text(p) as fh:
                json.dump({"a": 1}, fh)
            self.assertEqual(json.loads(_read(p)), {"a": 1})


class Generators(unittest.TestCase):
    def test_run_search_single_shot_jobs_submit_and_list(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            os.makedirs(w.out)
            victims = {p: plant(p, os.path.join(d, "elsewhere"))
                       for p in (w.lib_job, w.search_job, os.path.join(w.out, "submit.sh"),
                                 os.path.join(w.out, "diann_job_input_files.txt"),
                                 os.path.join(w.out, "search_provenance.json"))}
            w.generate(*HIGH)
            for p, victim in victims.items():
                self.assertEqual(_read(victim), VICTIM_TEXT, p)
                self.assertFalse(os.path.islink(p), p)
            self.assertTrue(os.access(os.path.join(w.out, "submit.sh"), os.X_OK))
            self.assertIn("#!/bin/bash", _read(w.search_job))

    def test_diann_parallel_steps_and_file_list(self):
        with tempfile.TemporaryDirectory() as d:
            # --window left out: step 1b measures it, so params.base.cfg is written too
            w = _Workspace(d, cfg_lines=LIBFREE_CFG, n_runs=3, ext=".d")
            out = os.path.join(d, "chain_out")
            os.makedirs(out)
            names = ("step2_firstpass.sbatch", "step5_report.sbatch", "submit.sh",
                     "file_list.txt", "params.base.cfg")
            victims = {n: plant(os.path.join(out, n), os.path.join(d, "elsewhere"))
                       for n in names}
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                                "--diann", w.diann, "--fasta", w.fasta, "--out", out,
                                "--cfg", w.cfg, "--partition", "high", "--account",
                                "genome-center-grp", "--raw", *w.runs],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            for n, victim in victims.items():
                self.assertEqual(_read(victim), VICTIM_TEXT, n)
                self.assertFalse(os.path.islink(os.path.join(out, n)), n)
                self.assertTrue(os.path.isfile(os.path.join(out, n)), n)
            self.assertEqual(_read(os.path.join(out, "file_list.txt")).split(), w.runs)

    def test_radiant_parallel_step_scripts(self):
        with tempfile.TemporaryDirectory() as d:
            runs = [_write(os.path.join(d, f"r{i}.mzML"), "") for i in range(2)]
            listing = _write(os.path.join(d, "list.txt"), "\n".join(runs) + "\n")
            fasta = _write(os.path.join(d, "db.fasta"), ">sp|P1|X\nPEPTIDEK\n")
            config = _write(os.path.join(d, "x.radiantConfig"), "")
            out = os.path.join(d, "out")
            D = os.path.join(out, "radiant_parallel")
            victim = plant(os.path.join(D, "step2_search.sbatch"), os.path.join(d, "elsewhere"))
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "radiant_parallel.py"),
                                "--runtime", "apptainer", "--image", "/x/radiant.sif",
                                "--raw-list", listing, "--fasta", fasta, "--config", config,
                                "--out", out, "--diann", "/opt/diann/diann-linux",
                                "--partition", "high", "--account", "genome-center-grp"],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(_read(victim), VICTIM_TEXT)
            self.assertFalse(os.path.islink(os.path.join(D, "step2_search.sbatch")))

    def test_diatracer_parallel_script_and_quoted_logs(self):
        with tempfile.TemporaryDirectory() as d:
            stage = os.path.join(d, "stage")
            os.makedirs(os.path.join(stage, "run1.d"))
            out = os.path.join(d, "out with space")
            victim = plant(os.path.join(out, "diatracer_array.sbatch"),
                           os.path.join(d, "elsewhere"))
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diatracer_parallel.py"),
                                "--stage", stage, "--out", out, "--jar", "/x/diatracer.jar",
                                "--partition", "high", "--account", "genome-center-grp"],
                               capture_output=True, text=True)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertEqual(_read(victim), VICTIM_TEXT)
            text = _read(os.path.join(out, "diatracer_array.sbatch"))
            self.assertIn(f'#SBATCH -o "{out}/diatracer_%A_%a.log"', text)
            self.assertIn(f'#SBATCH -e "{out}/diatracer_%A_%a.log"', text)
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diatracer_parallel.py"),
                                "--stage", stage, "--out", os.path.join(d, 'o"x'),
                                "--jar", "/x/diatracer.jar"], capture_output=True, text=True)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("REFUSING --out", p.stderr)


class SbatchOutputQuoted(unittest.TestCase):
    def test_a_session_folder_with_a_space_gets_a_quoted_log_path(self):
        with tempfile.TemporaryDirectory() as d:
            w = _Workspace(d)
            w.out = os.path.join(d, "Service run", "search out")
            w.generate(*HIGH)
            line = next(ln for ln in _read(w.search_job).splitlines()
                        if ln.startswith("#SBATCH --output="))
            self.assertEqual(line, f'#SBATCH --output="{w.out}/diann_search_%j.log"')


if __name__ == "__main__":
    unittest.main()
