#!/usr/bin/env python3
"""`session.py finalize` MAKES the reproducibility bundle (2.14.2, L5b review H1; Ruling 59).

Brett (2026-10-31): "I want to make sure an agents.md file and all the reproducibility and
methods are always done." Until this, finalize only *required* output/reproducibility/: a session
whose agent never ran provenance.py (step 10) -- the no-DE sessions above all -- was refused with
a pointer to "the flags SKILL.md step 10 lists", and the documented way out (--override-finish)
then left a session with no bundle that `session.py status` called "Complete". What these tests
pin:
  * a missing or stale bundle is built by finalize, running provenance.py with every flag read
    from the session's own records (workflow manifest, params, conditions, FASTA and its sidecar,
    the raw list, the search record and its report, commands.log, the DE record and its
    contrasts), for a session with a DE and one without;
  * a current bundle is kept as it is;
  * when a record it needs is absent, finalize refuses and names the exact provenance.py command;
  * the command the gates print is the one finalize runs (make_deposit.bundle_command, one
    definition);
  * a build that fails leaves the session not finalized (never stamped without a bundle).
provenance.py runs for real against a fake pipeline home (setup.json naming a fake python and
Rscript), so no conda env is listed and no R is started. Stdlib only, no network.
"""
import json
import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import finish_fixture                          # noqa: E402
import report_files as rf                      # noqa: E402
import test_deposit_package as tdp             # noqa: E402
import test_finish_gates as tfg                # noqa: E402  (FinalizeBase: a no-DE DIA session)

read, write = tfg.read, tfg.write
PY = sys.executable
BUNDLE_PART = "Reproducibility bundle (output/reproducibility)"


def line(manifest, part):
    """MANIFEST.txt's line for `part` (its level and note), or None."""
    for ln in manifest.splitlines():
        if ln.startswith("[") and part in ln and " -- " in ln:
            return ln
    return None


class Base(tfg.FinalizeBase):
    def run_manifest(self):
        return json.loads(read(os.path.join(self.out, "reproducibility", "run_manifest.json")))

    def status(self):
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "status", "--session",
                            self.sd], capture_output=True, text=True, env=self.env)
        return r.returncode, r.stdout

    def de(self):
        """A DE in this session, as run_de.R leaves its record."""
        write(os.path.join(self.p["de_dir"], "de_provenance.json"),
              json.dumps(dict(tdp.DE_PROV, input=os.path.join(self.p["search_out"],
                                                              "report.parquet"))))
        write(os.path.join(self.p["de_dir"], "DE_Treated_vs_Control.csv"), "Protein,logFC\nP1,1\n")


class FinalizeMakesTheBundle(Base):
    def test_a_no_de_session_without_a_bundle_gets_one(self):
        tfg.skill_report(self.out)                       # step 9 ran; step 10 never did
        r = self.finalize("--zip")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        rm = self.run_manifest()
        self.assertEqual(rm["engine"], "diann")
        self.assertEqual(rm["inputs"]["raw"], [x.rstrip("/") for x in
                                               finish_fixture_raws(self.sd)])
        self.assertEqual(rm["inputs"]["conditions"], self.p["conditions"])
        self.assertEqual(rm["inputs"]["fasta"], self.p["fasta"])
        self.assertEqual(rm["query"]["acquisition"], "DIA")
        self.assertEqual(rm["query"]["organism_taxid"], "9606")
        self.assertIsNone(rm["de"]["method"], "no DE: no --de-method")
        for name in ("MANIFEST.txt", "REPRODUCE.md", "reproduce.sh", "checksums/checksums.json"):
            self.assertTrue(os.path.isfile(os.path.join(self.out, "reproducibility", name)), name)
        checks = json.loads(read(os.path.join(self.out, "reproducibility", "checksums",
                                              "checksums.json")))
        self.assertEqual(checks["report"]["path"], os.path.join(self.p["search_out"],
                                                                "report.parquet"))
        self.assertNotIn("de_outputs", checks)
        ln = line(self.manifest(), BUNDLE_PART)
        self.assertTrue(ln and ln.startswith("[OK]"), self.manifest())
        self.assertIn("built by finalize", ln)
        st = json.loads(read(self.stamp_path()))
        self.assertEqual(st["outputs"]["output/reproducibility/run_manifest.json"],
                         rf._sha256(os.path.join(self.out, "reproducibility", "run_manifest.json")))
        rc, out = self.status()
        self.assertEqual(rc, 0, out)
        self.assertIn("Complete: every item is done.", out)

    def test_a_de_session_bundle_replays_the_de_record(self):
        self.de()
        tfg.skill_report(self.out)
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        rm = self.run_manifest()
        self.assertEqual(rm["de"]["method"], "dpc")
        self.assertEqual(rm["de"]["contrasts"], "Treated-Control")
        self.assertEqual((rm["de"]["q_cutoff"], rm["de"]["logfc"], rm["de"]["adjp"]),
                         ("0.01", "1", "0.05"))
        checks = json.loads(read(os.path.join(self.out, "reproducibility", "checksums",
                                              "checksums.json")))
        self.assertIn("de_provenance.json",
                      [os.path.basename(x["path"]) for x in checks["de_outputs"]])
        self.assertIn('--contrasts "Treated-Control"',
                      read(os.path.join(self.out, "reproducibility", "reproduce.sh")))

    def test_a_current_bundle_is_kept(self):
        tfg.skill_report(self.out)
        self.assertEqual(self.finalize().returncode, 0)
        before = read(os.path.join(self.out, "reproducibility", "run_manifest.json"))
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(read(os.path.join(self.out, "reproducibility", "run_manifest.json")),
                         before, "a current bundle is kept as it is")
        ln = line(self.manifest(), BUNDLE_PART)
        self.assertTrue(ln.startswith("[OK]"), ln)
        self.assertIn("present and current", ln)

    def test_a_bundle_made_before_the_de_was_run_again_is_built_again(self):
        self.de()
        tfg.skill_report(self.out)
        self.assertEqual(self.finalize().returncode, 0)
        rec = dict(tdp.DE_PROV, contrasts=["Control-Treated"], written="later")
        write(os.path.join(self.p["de_dir"], "de_provenance.json"), json.dumps(rec))
        state, why = rf.bundle_state(self.out)
        self.assertEqual(state, "stale")
        self.assertIn("de_provenance.json", why)
        self.assertEqual(rf.finish_check(self.out, ("reproducibility",))[0]["state"], "stale")
        rc, out = self.status()
        self.assertEqual(rc, 3, out)
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(self.run_manifest()["de"]["contrasts"], "Control-Treated")
        self.assertEqual(rf.bundle_state(self.out)[0], "present")
        self.assertIn("built by finalize", line(self.manifest(), BUNDLE_PART))

    def test_a_search_only_in_search_dda_is_found(self):
        """A Core member's submission's kind: the search is in output/search_dda/, and the raw list is
        only in its record."""
        os.remove(self.p["raw_list"])
        import shutil
        shutil.rmtree(self.p["search_out"])
        dda = os.path.join(self.out, "search_dda")
        raws = [os.path.join(self.root, "raw", r + ".d") for r in tdp.RUNS]
        write(os.path.join(dda, "report.parquet"), "y" * 50)
        write(os.path.join(dda, "search_provenance.json"), json.dumps({
            "engine": "diann", "version": "2.7.0", "files": raws,
            "result": {"engine": "diann", "report": os.path.join(dda, "report.parquet")}}))
        tfg.skill_report(self.out)
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        rm = self.run_manifest()
        self.assertEqual(rm["inputs"]["raw"], raws)
        checks = json.loads(read(os.path.join(self.out, "reproducibility", "checksums",
                                              "checksums.json")))
        self.assertEqual(checks["report"]["path"], os.path.join(dda, "report.parquet"))


class RecordsAbsent(Base):
    def test_no_raw_record_refuses_with_the_exact_command(self):
        os.remove(self.p["raw_list"])
        tfg.skill_report(self.out)
        r = self.finalize()
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("the reproducibility bundle (output/reproducibility/run_manifest.json)",
                      r.stderr)
        self.assertIn("finalize cannot build it without a record of the raw files", r.stderr)
        cmd = r.stderr.split("the reproducibility bundle", 1)[1].split(" -- make it with: ", 1)[1]
        self.assertTrue(cmd.startswith("python3 scripts/provenance.py --outdir "
                                       + self.p["repro_dir"]), cmd)
        self.assertIn("--raw <the raw files>", cmd)
        self.assertIn("--workflow-manifest " + self.p["workflow_manifest"], cmd)
        self.assertIn("--conditions " + self.p["conditions"], cmd)
        self.assertNotIn("...", cmd)
        self.assertNotIn("SKILL.md step 10 lists", r.stderr)
        self.assertFalse(os.path.exists(self.p["manifest_txt"]), "nothing written")

    def test_an_override_without_the_bundle_is_never_complete(self):
        os.remove(self.p["raw_list"])
        tfg.skill_report(self.out)
        write(os.path.join(self.out, "methods.md"), "# Methods\nWritten at step 9d.\n")
        r = self.finalize("--override-finish", "the raw files were deleted before step 10")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "reproducibility",
                                                     "run_manifest.json")))
        ln = line(self.manifest(), BUNDLE_PART)
        self.assertTrue(ln and ln.startswith("[SKIPPED]"), self.manifest())
        rc, out = self.status()
        self.assertEqual(rc, 0, out)
        self.assertNotIn("every item is done", out)
        self.assertIn("The reproducibility bundle (output/reproducibility/run_manifest.json) "
                      "(overridden: the raw files were deleted before step 10)", out)

    def test_a_failed_build_leaves_the_session_not_finalized(self):
        write(self.p["fasta_meta"], "{ not json")       # provenance.py refuses a broken sidecar
        tfg.skill_report(self.out)
        r = self.finalize()
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        ln = line(self.manifest(), BUNDLE_PART)
        self.assertTrue(ln and ln.startswith("[SKIPPED]"), self.manifest())
        self.assertIn("provenance.py failed", ln)
        self.assertFalse(os.path.exists(self.stamp_path()))
        self.assertIn("NOT FINALIZED", r.stderr)


class OneCommand(Base):
    def test_the_command_the_gates_print_is_the_one_finalize_runs(self):
        import make_deposit
        gap = rf.finish_gaps(self.out, ("reproducibility",))[0]
        cmd = gap.split(" -- make it with: ", 1)[1]
        self.assertTrue(cmd.startswith("python3 scripts/provenance.py --outdir "), cmd)
        self.assertIn("(step 10; session.py finalize also makes it)", cmd)
        want, missing, _ = make_deposit.bundle_command(make_deposit.gather(self.sd))
        self.assertEqual(missing, [])
        for arg in want[2:]:
            self.assertIn(arg, cmd)

    def test_a_long_raw_list_is_read_from_the_session(self):
        raws = [f"/data/runs/S{i:02d}.d" for i in range(1, 9)]
        write(self.p["raw_list"], "# raw\n" + "\n".join(raws) + "\n")
        gap = rf.finish_gaps(self.out, ("reproducibility",))[0]
        self.assertIn("--raw $(grep -v '^#' " + self.p["raw_list"] + ")", gap)
        self.assertNotIn("S08.d", gap)


def finish_fixture_raws(sd):
    """The raw paths the session's input/raw_files.txt holds."""
    import session
    return session.read_raw_list(sd)


if __name__ == "__main__":
    unittest.main()
