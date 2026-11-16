"""One definition of "this is a SLURM cluster login node" and of where SLURM lives (cleanup audit E2).

on_compute.py holds the definition (SLURM_BIN_DEFAULT, CLUSTER_CONF_DEFAULT, on_login_node). Bash
cannot import it, so setup.sh, acquire_tools.sh and ensure_dotnet8.sh carry the same three-way test
for the case with no python3, and record_run.py is piped to a remote python with no siblings, so it
keeps a literal path too. This ties every copy to the one definition: change the path or the marker
file in on_compute.py and this fails until all of them follow.
"""
import os
import re
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
import on_compute  # noqa: E402


def read(name):
    with open(os.path.join(SCRIPTS, name), encoding="utf-8") as fh:
        return fh.read()


class OneDefinition(unittest.TestCase):
    BASH_TWINS = ("setup.sh", "acquire_tools.sh", "ensure_dotnet8.sh")

    def test_every_bash_twin_uses_the_same_path_and_marker_file(self):
        for name in self.BASH_TWINS:
            text = read(name)
            self.assertIn("${HPCCF_SLURM_BIN:-%s}" % on_compute.SLURM_BIN_DEFAULT, text, name)
            self.assertIn("${UCD_HPC_CLUSTER_CONF:-%s}" % on_compute.CLUSTER_CONF_DEFAULT, text, name)

    def test_every_bash_twin_tests_the_same_three_things(self):
        # sbatch on PATH, sbatch at the fixed path, the marker file: no more, no fewer
        for name in self.BASH_TWINS:
            text = read(name)
            self.assertRegex(text, r"(have|command -v) sbatch", name)
            self.assertRegex(text, r'-x "[^"]*/sbatch"|-x "\$UCD_SLURM_BIN/sbatch"', name)
            self.assertRegex(text, r"-r \"[^\"]*\"", name)

    def test_the_standalone_scripts_name_the_same_slurm_directory(self):
        self.assertIn(on_compute.SLURM_BIN_DEFAULT + "/sacct", read("record_run.py"))
        # run_search takes it from on_compute: it must not carry its own copy
        self.assertNotIn(on_compute.SLURM_BIN_DEFAULT, read("run_search.py"))

    def test_the_venv_builder_honours_the_same_override_as_the_others(self):
        text = read("build_directlfq_venv.sh")
        self.assertIn("${UCD_HPC_CLUSTER_CONF:-%s}" % on_compute.CLUSTER_CONF_DEFAULT, text)

    def test_python_says_login_node_for_each_of_the_three_signs_and_only_those(self):
        with tempfile.TemporaryDirectory() as d:
            bindir = os.path.join(d, "bin")
            os.makedirs(bindir)
            conf = os.path.join(d, "cluster_name.conf")
            base = {"PATH": os.path.join(d, "empty"), "HPCCF_SLURM_BIN": bindir,
                    "UCD_HPC_CLUSTER_CONF": conf}
            self.assertFalse(on_compute.on_login_node(dict(base)))            # none of the signs
            with open(conf, "w") as fh:
                fh.write("x")
            self.assertTrue(on_compute.on_login_node(dict(base)))             # the marker file
            os.unlink(conf)
            sb = os.path.join(bindir, "sbatch")
            with open(sb, "w") as fh:
                fh.write("#!/bin/sh\n")
            os.chmod(sb, 0o755)
            self.assertTrue(on_compute.on_login_node(dict(base)))             # sbatch at the path
            self.assertFalse(on_compute.on_login_node(dict(base, SLURM_JOB_ID="1")))   # in a job


if __name__ == "__main__":
    unittest.main()
