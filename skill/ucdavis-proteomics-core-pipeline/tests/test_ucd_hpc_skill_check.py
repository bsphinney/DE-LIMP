#!/usr/bin/env python3
"""check_access.sh says whether HPC@UCD's own agent skill (ucd-hpc) is installed, and when it
is required.

HPC@UCD made its `ucd-hpc` skill mandatory for agents on Farm, Franklin and Hive on 2026-10-16
("agents that do not follow it put the user's account at risk of suspension"); on cluster policy
it supersedes this skill. Claude reads it only when it is installed where Claude Code looks for
skills -- on the computer Claude Code runs on: the laptop in hive_remote, HIVE in hive_local.
So step 0a's check_access.sh reports `ucd_hpc_skill`: installed or not, required or not (anyone
working on a UC Davis cluster), and the install line from the ucd-hpc README.

No test contacts a host: HOME is a temp folder, no HIVE login is given (so no ssh), and "on a
cluster" is a stub sbatch on PATH plus a stand-in for /opt/hpccf/etc/cluster_name.conf.
"""
import json
import os
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
CHECK_ACCESS = os.path.join(os.path.dirname(HERE), "scripts", "check_access.sh")
CLONE = "git clone https://github.com/ucdavis/ucdhpc-cluster-skill"


class UcdHpcSkillCheck(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.home = self._tmp.name
        self.bin = os.path.join(self.home, "fakebin")
        os.makedirs(self.bin)
        self.conf = os.path.join(self.home, "cluster_name.conf")

    def tearDown(self):
        self._tmp.cleanup()

    def on_cluster(self, name="hive"):
        """A stub sbatch on PATH (check_access.sh's on_hive) and HPCCF's cluster file."""
        p = os.path.join(self.bin, "sbatch")
        with open(p, "w") as fh:
            fh.write("#!/bin/sh\nexit 0\n")
        os.chmod(p, 0o755)
        with open(self.conf, "w") as fh:
            fh.write(f'export CLUSTER_NAME="{name}"\n')

    def install(self, via_clone=True):
        """The README's install: a clone, and a link to its ucd-hpc/ in ~/.claude/skills."""
        repo = os.path.join(self.home, "ucdhpc-cluster-skill")
        skill = os.path.join(repo, "ucd-hpc")
        os.makedirs(skill)
        with open(os.path.join(skill, "SKILL.md"), "w") as fh:
            fh.write("---\nname: ucd-hpc\ndescription: HPCCF guidelines\n---\n# UC Davis HPC\n")
        if via_clone:
            os.makedirs(os.path.join(self.home, ".claude", "skills"))
            os.symlink(skill, os.path.join(self.home, ".claude", "skills", "ucd-hpc"))
        return skill

    def run_check(self, **extra):
        env = {"HOME": self.home, "PATH": self.bin + os.pathsep + os.environ["PATH"],
               "UCD_HPC_CLUSTER_CONF": self.conf}
        for k in ("SLURM_JOB_ID", "HIVE_USER", "HIVE_KEY", "CLAUDE_CONFIG_DIR"):
            env.pop(k, None)
        env.update(extra)
        r = subprocess.run(["bash", CHECK_ACCESS], capture_output=True, text=True, env=env,
                           timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        return json.loads(r.stdout)["ucd_hpc_skill"]

    @unittest.skipIf(shutil.which("sbatch"), "this machine has SLURM: it would read as a cluster")
    def test_a_laptop_without_a_cluster_does_not_need_it(self):
        u = self.run_check()
        self.assertIs(u["installed"], False)
        self.assertIs(u["required"], False)
        self.assertEqual(u["path"], os.path.join(self.home, ".claude", "skills", "ucd-hpc"))

    def test_on_a_uc_davis_cluster_without_it_it_is_required_with_the_readme_install(self):
        self.on_cluster()
        u = self.run_check()
        self.assertIs(u["installed"], False)
        self.assertIs(u["required"], True)
        self.assertIn(f"{CLONE} ~/ucdhpc-cluster-skill", u["install"])
        self.assertIn("mkdir -p ~/.claude/skills", u["install"])
        self.assertIn("ln -s ~/ucdhpc-cluster-skill/ucd-hpc ~/.claude/skills/ucd-hpc", u["install"])
        self.assertIn("cp -r", u["note"])                  # Windows: a copy when ln -s cannot link
        self.assertIn("supersedes", u["note"])

    def test_another_slurm_cluster_does_not_require_it(self):
        self.on_cluster()
        os.remove(self.conf)                               # sbatch, but no HPCCF cluster file
        self.assertIs(self.run_check()["required"], False)

    def test_an_installed_skill_is_found_through_its_link(self):
        self.on_cluster()
        self.install()
        u = self.run_check()
        self.assertIs(u["installed"], True)
        self.assertIs(u["required"], True)
        self.assertIsNone(u["install"])
        self.assertIn("git -C ~/ucdhpc-cluster-skill pull --ff-only", u["note"])

    def test_a_clone_without_the_link_needs_only_the_link(self):
        self.on_cluster()
        self.install(via_clone=False)
        u = self.run_check()
        self.assertIs(u["installed"], False)
        self.assertNotIn(CLONE, u["install"])
        self.assertIn("ln -s ~/ucdhpc-cluster-skill/ucd-hpc ~/.claude/skills/ucd-hpc", u["install"])

    def test_a_dangling_link_is_not_installed(self):
        self.on_cluster()
        skill = self.install()
        shutil.rmtree(os.path.dirname(skill))              # the clone deleted, the link left
        u = self.run_check()
        self.assertIs(u["installed"], False)
        self.assertIn(CLONE, u["install"])

    def test_a_dangling_link_is_removed_first_by_the_install_line(self):
        self.on_cluster()
        shutil.rmtree(os.path.dirname(self.install()))
        self.assertIn("rm -f ~/.claude/skills/ucd-hpc && mkdir -p ~/.claude/skills && ln -s ",
                      self.run_check()["install"])

    def test_a_folder_of_the_users_in_the_way_is_never_deleted(self):
        self.on_cluster()
        os.makedirs(os.path.join(self.home, ".claude", "skills", "ucd-hpc", "notes"))
        u = self.run_check()
        self.assertIs(u["installed"], False)
        self.assertIsNone(u["install"])
        self.assertIn("exists but is not the ucd-hpc skill", u["note"])

    def test_a_skill_of_another_name_is_not_it(self):
        self.on_cluster()
        skill = self.install()
        with open(os.path.join(skill, "SKILL.md"), "w") as fh:
            fh.write("---\nname: something-else\n---\n")
        self.assertIs(self.run_check()["installed"], False)

    def test_claude_config_dir_moves_the_skills_folder(self):
        self.on_cluster()
        cfg = os.path.join(self.home, "claude-config")
        u = self.run_check(CLAUDE_CONFIG_DIR=cfg)
        self.assertEqual(u["path"], os.path.join(cfg, "skills", "ucd-hpc"))
        self.assertIn(f"mkdir -p {cfg}/skills", u["install"])


if __name__ == "__main__":
    unittest.main()
