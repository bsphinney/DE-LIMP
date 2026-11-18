#!/usr/bin/env python3
"""
2.15.1 review, minor 10: "is this HIVE?" had two readings of HPC@UCD's cluster file -- the shell
(hive_marker.sh) took the first 40 characters of the whole file, Python (run_search.cluster_is_hive)
the first 200 -- so a longer file could make them disagree, and a comment that named "hive" made
any cluster HIVE. HPCCF's own ucd-hpc skill reads the file as shell (`. cluster_name.conf; echo
$CLUSTER_NAME`). Now one rule, written twice (bash cannot import Python): a CLUSTER_NAME=<name>
line, else the bare name on the first non-comment line. This runs both on the same files and
requires the same name and the same HIVE verdict.
"""
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import run_search  # noqa: E402

CASES = {
    "hive\n": "hive",
    "Hive": "hive",
    "  farm  \n": "farm",
    "CLUSTER_NAME=hive\n": "hive",
    "export CLUSTER_NAME=\"franklin\"\n": "franklin",
    "CLUSTER_NAME='hive' # set by HPCCF\n": "hive",
    "# HPCCF cluster name file: one of hive, farm, franklin -- maintained by HPCCF, do not edit\n"
    "CLUSTER_NAME=farm\n": "farm",
    "# a long header that names hive in a comment, longer than forty characters by far\n"
    "farm\n": "farm",
    "# only a comment naming hive\n": "",
    "": "",
    "CLUSTER_NAME=\n": "",
}


def shell(conf, what):
    env = dict(os.environ, UCD_HPC_CLUSTER_CONF=conf, HIVE_GROUP_DIR="/nonexistent/grp")
    script = f'. "{SCRIPTS}/hive_marker.sh"; ' + (
        "cluster_name" if what == "name" else "if is_hive; then echo yes; else echo no; fi")
    return subprocess.run(["bash", "-c", script], capture_output=True, text=True,
                          env=env).stdout.strip()


class OneRule(unittest.TestCase):
    def test_shell_and_python_read_the_same_name_and_verdict(self):
        with tempfile.TemporaryDirectory() as d:
            conf = os.path.join(d, "cluster_name.conf")
            for text, want in CASES.items():
                with self.subTest(text=text):
                    with open(conf, "w") as fh:
                        fh.write(text)
                    self.assertEqual(run_search.cluster_name_from_conf(text), want)
                    self.assertEqual(shell(conf, "name"), want)
                    env = {"UCD_HPC_CLUSTER_CONF": conf, "HIVE_GROUP_DIR": "/nonexistent/grp"}
                    py = run_search.cluster_is_hive(env)
                    self.assertEqual(py, "hive" in want)
                    self.assertEqual(shell(conf, "hive"), "yes" if py else "no")

    def test_no_file_falls_back_to_the_group_folder_in_both(self):
        with tempfile.TemporaryDirectory() as d:
            env = {"UCD_HPC_CLUSTER_CONF": os.path.join(d, "none"), "HIVE_GROUP_DIR": d}
            self.assertTrue(run_search.cluster_is_hive(env))
            r = subprocess.run(["bash", "-c", f'. "{SCRIPTS}/hive_marker.sh"; is_hive && echo yes'],
                               capture_output=True, text=True, env=dict(os.environ, **env))
            self.assertEqual(r.stdout.strip(), "yes")


if __name__ == "__main__":
    unittest.main()
