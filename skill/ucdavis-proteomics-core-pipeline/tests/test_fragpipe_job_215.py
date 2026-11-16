#!/usr/bin/env python3
"""
The FragPipe search job, 2.15 (lane F3; review H3, H6, M2, M4, M10, L6, L7, the L11 notes;
Ruling 90).

  * Each FragPipe start runs with its OWN settings folder (XDG_CONFIG_HOME under
    <out>/fragpipe_settings), never the account's shared ~/.config/FragPipe -- which every start
    rewrote (H5) and which filled in a workflow's missing database with the LAST search's (H5
    finding 11). It is given the tools folder and bundled DIA-NN from tools.json fragpipe_tools,
    the one carrier (M4): the FRAGPIPE_TOOLS_FOLDER variable no longer reaches the job.
  * The job records `java -version` of the java FragPipe's launcher runs and stops below 11.
  * A target-only --fasta gets its decoys in the job, from the install's Philosopher (never
    --contam), and the workflow FragPipe runs names that .fas; a half decoy set is refused; the
    database is hashed in the job (or taken from its fetch_fasta sidecar), not on the login node.
    reproduce.sh's search, which rebuilds a target-only FASTA, therefore runs.
  * Refused before anything is written: a missing --params (L6), a DDA workflow that does not
    say whether it quantifies, two runs FragPipe would name as one experiment (it maps every
    character outside [A-Za-z0-9_] to '_': A-1 and A_1 were pooled silently), a tools.json without
    fragpipe_tools, a settings folder that is a link.
  * The decoy tag's source is recorded (the workflow's, or FragPipe's default rev_, said so).
  * The contaminant digest warning covers FragPipe (L7); a FragPipe job's time follows its runs
    (M10), as a default said so.
"""
import datetime
import hashlib
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
sys.path.insert(0, HERE)

try:
    import fragpipe_db  # noqa: E402
except ImportError:         # a fail-first base without fragpipe_db.py: each test fails on its own
    fragpipe_db = None
import run_search  # noqa: E402
import scratch_files  # noqa: E402
from fake_fragpipe_install import fake_install, tools_json  # noqa: E402
from test_fragpipe_dda_stage import FP_WORKFLOW, StageHarness, _write, listing  # noqa: E402

TARGETS = ">sp|P1|ONE_HUMAN one\nPEPTIDER\n>sp|P2|TWO_HUMAN two\nSAMPLERK\n"
COMPLETE = TARGETS + ">rev_sp|P1|ONE_HUMAN one\nREDITPEP\n>rev_sp|P2|TWO_HUMAN two\nKRELPMAS\n"
HALF = TARGETS + ">rev_sp|P1|ONE_HUMAN one\nREDITPEP\n"


def lfq(text=None):
    """The fixture workflow with label-free quantification on (a quantifying DDA workflow)."""
    text = text or open(FP_WORKFLOW, newline="").read()
    return text.replace("quantitation.run-label-free-quant=false",
                        "quantitation.run-label-free-quant=true")


def without(text, key):
    return "".join(ln for ln in text.splitlines(keepends=True) if not ln.startswith(key + "="))


def tryptic(text=None):
    """The fixture workflow searching the skill's shared DDA digest (KR, not before P, 2 missed,
    7-30): FragPipe's known default for the contaminant-digest check."""
    text = text or open(FP_WORKFLOW, newline="").read()
    for k, v in (("msfragger.num_enzyme_termini", "2"), ("msfragger.search_enzyme_cut_1", "KR"),
                 ("msfragger.search_enzyme_nocut_1", "P"),
                 ("msfragger.search_enzyme_name_1", "trypsin"),
                 ("msfragger.digest_max_length", "30")):
        text = re.sub(rf"(?m)^{re.escape(k)}=.*$", f"{k}={v}\r", text)
    # an explicit rule is written without FragPipe's enzyme dropdown, which it would load over
    # the rule (C0 fact 1) -- as make_presets.py writes one
    return without(text, "msfragger.misc.fragger.enzyme-dropdown-1")


class Experiments(unittest.TestCase):
    def test_fragpipes_rule(self):
        """InputLcmsFile.java (24.0): trim, then every UTF-16 unit outside [A-Za-z0-9_] -> '_'"""
        for name, want in (("10282026__60SPD_Dda-LW1CT_S3-A2_1_25599",
                            "10282026__60SPD_Dda_LW1CT_S3_A2_1_25599"),
                           ("A-1", "A_1"), ("A_1", "A_1"), (" A.1 \t", "A_1"),
                           ("s 1+x", "s_1_x"), ("Größe", "Gr__e"), ("a\U0001F600b", "a__b")):
            self.assertEqual(run_search.fragpipe_experiment(name), want, name)


class Refusals(StageHarness):
    """Each refused before anything is written (the search folder is not even made)."""

    def assertRefused(self, p, *words):
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        for w in ("REFUSED", "nothing was written") + words:
            self.assertIn(w, p.stderr)
        self.assertFalse(os.path.exists(self.out), "the search folder was made")
        self.assertFalse(os.path.exists(self.job))
        self.assertEqual(listing(self.raw), self.before)

    def test_a_missing_params_file(self):
        """L6: it crashed after the stage and the manifest were written"""
        params = os.path.join(self.d, "no_such.workflow")
        tools = tools_json(self.d, self.fp, self.fp_install)
        bundle = _write(os.path.join(self.d, "bundle.json"), json.dumps(
            {"engine": {"name": "fragpipe"}, "acquisition": "DDA"}))
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"), "--tools",
                            tools, "--bundle", bundle, "--params", params, "--fasta",
                            _write(os.path.join(self.d, "db.fasta"), COMPLETE), "--out",
                            self.out, "--files", *self.files, "--engine", "fragpipe",
                            "--allow-damaged-tdf", "--sbatch", self.job],
                           capture_output=True, text=True, timeout=120)
        self.assertRefused(p, "no such file")

    def test_a_dda_workflow_that_does_not_say_whether_it_quantifies(self):
        p = self.generate(workflow=without(open(FP_WORKFLOW, newline="").read(),
                                           "quantitation.run-label-free-quant"))
        self.assertRefused(p, "quantitation.run-label-free-quant")

    def test_two_runs_fragpipe_would_name_alike(self):
        a = _write(os.path.join(self.raw, "A-1.mzML"), "")
        b = _write(os.path.join(self.raw, "A_1.mzML"), "")
        self.before = listing(self.raw)
        p = self.generate(files=[a, b], workflow=lfq())
        self.assertRefused(p, "A-1, A_1", "experiment 'A_1'", "ONE experiment")

    def test_alike_names_are_fine_when_the_search_pools_its_runs(self):
        """identifications only: one pooled experiment, nothing to collide"""
        a = _write(os.path.join(self.raw, "A-1.mzML"), "")
        b = _write(os.path.join(self.raw, "A_1.mzML"), "")
        p = self.generate(files=[a, b])
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_a_tools_json_without_fragpipe_tools(self):
        tools = _write(os.path.join(self.d, "old_tools.json"), json.dumps({"fragpipe": self.fp}))
        p = self.generate(tools=tools)
        self.assertRefused(p, "fragpipe_tools", "acquire_tools.sh")
        # Ruling 104 (4): it says which step re-runs it, on both paths
        self.assertIn("SKILL.md step 5", p.stderr)
        self.assertIn("reproduce.sh step 3", p.stderr)

    def test_fragpipe_tools_without_a_tools_folder(self):
        fp = fake_install(os.path.join(self.d, "unlicensed"), licensed=False)
        p = self.generate(tools=tools_json(self.d, self.fp, fp, name="unlicensed.json"))
        self.assertRefused(p, "no tools folder", "no MSFragger/IonQuant jars")

    def test_half_a_decoy_set(self):
        p = self.generate(fasta=_write(os.path.join(self.d, "half.fasta"), HALF))
        self.assertRefused(p, "1 decoys", "2 targets", "sp|P2|TWO_HUMAN", "half a decoy set")

    def test_a_target_only_fasta_without_a_philosopher(self):
        fp = fake_install(os.path.join(self.d, "nophi"))
        fp["fragpipe_tools"]["philosopher"] = None
        p = self.generate(tools=tools_json(self.d, self.fp, fp, name="nophi.json"))
        self.assertRefused(p, "no decoys", "fragpipe_tools.philosopher")

    def test_a_settings_folder_that_is_a_link(self):
        os.makedirs(self.out)
        os.symlink(os.path.join(self.d, "home"), os.path.join(self.out, "fragpipe_settings"))
        p = self.generate()
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("is a link", p.stderr)
        self.assertFalse(os.path.exists(self.job))


class SettingsFolder(StageHarness):
    """Ruling 90: the per-start settings folder is on the job line (with --out holding a space),
    the tools come from tools.json, the shared ~/.config/FragPipe is never touched, the folder is
    scratch and goes after a finished search."""

    def setUp(self):
        super().setUp()
        self.out = os.path.join(self.d, "my out")                       # one word with a space

    def test_the_job_line(self):
        p = self.generate(env_extra={"FRAGPIPE_TOOLS_FOLDER": "/elsewhere/tools"})
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        text = open(self.job).read()
        word = f"XDG_CONFIG_HOME='{self.out}/fragpipe_settings'/job_\"${{SLURM_JOB_ID:-inline}}\""
        self.assertIn(word + " " + self.fp + " --headless", text)
        tools = self.fp_install["tools_folder"]
        self.assertIn(f"--config-tools-folder {tools}", text)
        self.assertIn("--config-diann " + self.fp_install["fragpipe_tools"]["diann"]["path"],
                      text)
        self.assertNotIn("/elsewhere/tools", text)             # M4: tools.json, never the env
        prov = json.load(open(os.path.join(self.out, "search_provenance.json")))
        self.assertEqual(prov["fragpipe_config"]["tools_folder"], tools)
        self.assertEqual(prov["fragpipe_settings"]["env"], "XDG_CONFIG_HOME")

    def test_the_job_runs_fragpipe_in_its_own_folder_and_removes_it(self):
        self.assertEqual(self.generate().returncode, 0)
        home = os.path.join(self.d, "home")
        shared = os.path.join(home, ".config", "FragPipe")
        r = self.run_job()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(self.out, ".fake_fragpipe_ran")) as fh:
            seen = fh.read()
        self.assertIn(f"XDG_CONFIG_HOME={self.out}/fragpipe_settings/job_inline\n", seen)
        self.assertFalse(os.path.exists(shared), "the shared settings folder was written")
        # FragPipe wrote its settings there (the fake does, only under XDG_CONFIG_HOME) -- and
        # the job removed them after the finished search
        self.assertFalse(os.path.exists(os.path.join(self.out, "fragpipe_settings")))

    def test_a_search_that_failed_keeps_its_settings(self):
        self.assertEqual(self.generate().returncode, 0)
        r = self.run_job(FAKE_FP_RC="2")
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertTrue(os.path.isdir(os.path.join(self.out, "fragpipe_settings", "job_inline")))

    def test_the_diatracer_route_keeps_its_settings_without_all_jobs_done(self):
        """review L7: on the diaTracer route (DIA on timsTOF .d, no stage) the settings folder
        went after any exit 0; a FragPipe that exits 0 after a crash prints no ALL JOBS DONE,
        and its settings folder is now kept, said so. A finished search still removes it."""
        dotd = [os.path.join(self.raw, "S3.d")]
        settings = os.path.join(self.out, "fragpipe_settings", "job_inline")
        self.assertEqual(self.generate(files=dotd, acquisition="DIA").returncode, 0)
        self.assertTrue(os.path.isdir(os.path.join(self.out, "diatracer_stage")))
        r = self.run_job(FAKE_FP_DONE="0")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertTrue(os.path.isdir(settings), "removed without ALL JOBS DONE")
        self.assertIn("printed no ALL JOBS DONE", r.stderr)
        r = self.run_job()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "fragpipe_settings")))
        self.assertIn("ALL JOBS DONE", open(os.path.join(
            self.out, run_search.FRAGPIPE_CONSOLE_LOG)).read())

    def test_the_folders_are_scratch(self):
        for name in ("fragpipe_settings", "fragpipe_db"):
            self.assertTrue(scratch_files.is_scratch_dir(name), name)
            self.assertIn(name, scratch_files.LABEL)
            self.assertIn("excluded", scratch_files.dir_reason(name, "default"))
        self.assertTrue(scratch_files.is_scratch_file("2026-11-04-decoys-x.fasta.fas.1.pepindex"))


class SpectralLibraryPython(StageHarness):
    """A spectral-library workflow (DIA_SpecLib_*: EasyPQP) gets FragPipe's --config-python from
    tools.json fragpipe_tools.python: its own settings folder no longer supplies the account's
    bin-python, and EasyPQP on the node's bare python3 fails."""

    WF = ("workflow.misc.save-sdrf=false\ndiann.run-dia-nn=true\nspeclibgen.run-speclibgen=true\n"
          "quantitation.run-label-free-quant=false\ndatabase.decoy-tag=rev_\n")

    def test_it_is_passed(self):
        mz = _write(os.path.join(self.raw, "S4.mzML"), "")
        p = self.generate(files=[mz], acquisition="DIA", workflow=self.WF)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("--config-python " + self.fp_install["fragpipe_tools"]["python"]["path"],
                      open(self.job).read())

    def test_a_workflow_without_a_library_is_not_given_one(self):
        self.assertEqual(self.generate().returncode, 0)
        self.assertNotIn("--config-python", open(self.job).read())

    def test_a_python_without_easypqp_is_refused(self):
        """F2: python recorded but no EasyPQP found in it"""
        fp = fake_install(os.path.join(self.d, "noeasypqp"))
        fp["fragpipe_tools"]["python"]["easypqp"] = None
        mz = _write(os.path.join(self.raw, "S4.mzML"), "")
        p = self.generate(files=[mz], acquisition="DIA", workflow=self.WF,
                          tools=tools_json(self.d, self.fp, fp, name="noeasypqp.json"))
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("has no EasyPQP", p.stderr)
        self.assertFalse(os.path.exists(self.out))

    def test_none_recorded_is_refused(self):
        fp = fake_install(os.path.join(self.d, "nopy"))
        fp["fragpipe_tools"]["python"] = None
        mz = _write(os.path.join(self.raw, "S4.mzML"), "")
        p = self.generate(files=[mz], acquisition="DIA", workflow=self.WF,
                          tools=tools_json(self.d, self.fp, fp, name="nopy.json"))
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("speclibgen.run-speclibgen=true", p.stderr)
        self.assertIn("EasyPQP", p.stderr)
        self.assertFalse(os.path.exists(self.out))


JAVA_LINES = (('openjdk version "25.0.4.1" 2026-08-18', 25), ('java version "1.8.0_292"', 8),
              ('openjdk version "11.0.2"', 11), ('openjdk version "21" 2023-09-19', 21),
              ('openjdk version "17-ea" 2021-03-16', 17), ("bash: java: not found", None))


class Java(StageHarness):
    def test_java_major(self):
        for text, want in JAVA_LINES:
            self.assertEqual(fragpipe_db.java_major(text), want, text)

    def test_one_java_floor_and_one_reading(self):
        """the floor and the parser are engine_choice's (detect_env.sh, resolve_defaults and the
        job record read them); the job's shell line gives the same answers"""
        import engine_choice
        self.assertIs(fragpipe_db.java_major, engine_choice.java_major)
        self.assertEqual(fragpipe_db.MIN_JAVA, engine_choice.FRAGPIPE_JAVA_MIN)
        self.assertIn(f'-lt {engine_choice.FRAGPIPE_JAVA_MIN} ]',
                      "\n".join(run_search.fragpipe_java_lines()))
        for text, want in JAVA_LINES:
            sh = subprocess.run(["bash", "-c", run_search.FRAGPIPE_JAVA_AWK
                                 + '; printf "%s" "$_fp_java_major"'],
                                env=dict(os.environ, _fp_java_text=text),
                                capture_output=True, text=True)
            self.assertEqual(sh.stdout, "" if want is None else str(want), text)

    def test_below_java_11_nothing_is_copied_or_searched(self):
        self.assertEqual(self.generate().returncode, 0)
        r = self.run_job(FAKE_JAVA_VERSION="1.8.0_292")
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertIn("needs Java 11 or newer", r.stderr)
        self.assertIn('1.8.0_292', r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, ".fake_fragpipe_ran")))
        self.assertNotIn("S1-a.raw", os.listdir(self.stage()))      # the copy never started

    def test_no_java_at_all_is_refused_too(self):
        self.assertEqual(self.generate().returncode, 0)
        r = self.run_job(JAVA_HOME=os.path.join(self.d, "no_jdk"))
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertIn("needs Java 11 or newer", r.stderr)

    def test_the_java_is_recorded_and_merged_by_adapt_only(self):
        self.assertEqual(self.generate().returncode, 0)
        r = self.run_job(FAKE_JAVA_VERSION="17.0.9")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("[fragpipe] java", r.stdout)
        job = json.load(open(os.path.join(self.out, fragpipe_db.JOB_RECORD)))
        self.assertEqual(job["java"]["major"], 17)
        tools = tools_json(self.d, self.fp, self.fp_install)
        bundle = os.path.join(self.d, "bundle.json")
        a = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"), "--tools",
                            tools, "--bundle", bundle, "--params",
                            os.path.join(self.d, "search.workflow"), "--fasta",
                            os.path.join(self.d, "db.fasta"), "--out", self.out, "--files",
                            *self.files, "--engine", "fragpipe", "--adapt-only"],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(a.returncode, 0, a.stdout + a.stderr)
        prov = json.load(open(os.path.join(self.out, "search_provenance.json")))
        self.assertEqual(prov["fragpipe_java"]["major"], 17)
        db = prov["fragpipe_database"]
        self.assertEqual(db["sha256"], job["database"]["sha256"])
        self.assertEqual(db["sha256_source"], "measured in the job, before FragPipe started")


class Decoys(StageHarness):
    """H3: a target-only FASTA is searched; its decoys are built in the job."""

    def test_the_job_builds_them_and_fragpipe_searches_that_fas(self):
        p = self.generate()                                    # db.fasta: targets only
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        prov = json.load(open(os.path.join(self.out, "search_provenance.json")))
        db = prov["fragpipe_database"]
        import fetch_fasta
        name = os.path.basename(db["path"])
        self.assertTrue(fetch_fasta.PHILOSOPHER_DECOY_FAS.match(name), name)  # its own name
        self.assertEqual(name, f"{datetime.date.today().isoformat()}-decoys-db.fasta.fas")
        self.assertEqual(os.path.dirname(db["path"]), os.path.join(self.out, "fragpipe_db"))
        self.assertTrue(db["decoys_built_in_job"])
        self.assertFalse(db["decoy_build"]["contaminants_added"])
        # review L1: before the job has run nothing was measured, and the record says so -- a job
        # that dies before measuring leaves a null sha256 that never claims a measurement
        self.assertIsNone(db["sha256"])
        self.assertFalse(db["sha256_source"].startswith("measured"), db["sha256_source"])
        self.assertTrue(db["sha256_source"].startswith("not yet measured"), db["sha256_source"])
        with open(prov["fragpipe_workflow"]["run"]) as fh:
            self.assertIn("database.db-path=" + db["path"].replace(":", "\\:"), fh.read())
        text = open(self.job).read()
        build = [ln for ln in text.splitlines() if "database --custom" in ln]
        self.assertEqual(len(build), 1, build)
        self.assertNotIn("--contam", build[0])
        r = self.run_job()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertTrue(os.path.isfile(db["path"]))
        job = json.load(open(os.path.join(self.out, fragpipe_db.JOB_RECORD)))
        self.assertEqual(job["database"]["state"], "target_decoy")
        self.assertEqual(job["database"]["n_decoys"], job["database"]["n_targets"])
        self.assertEqual(job["database"]["sha256"],
                         hashlib.sha256(open(db["path"], "rb").read()).hexdigest())
        self.assertFalse(os.path.exists(os.path.join(self.out, "fragpipe_db", ".build")))

    def test_half_a_set_from_the_build_stops_the_job_before_fragpipe(self):
        self.assertEqual(self.generate(
            fasta=_write(os.path.join(self.d, "two.fasta"), TARGETS)).returncode, 0)
        r = self.run_job(FAKE_PHI_HALF="1")
        self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("not a complete target+decoy database", r.stdout + r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, ".fake_fragpipe_ran")))

    def test_a_complete_set_is_searched_as_it_is(self):
        fasta = _write(os.path.join(self.d, "td.fasta"), COMPLETE)
        p = self.generate(fasta=fasta)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        db = json.load(open(os.path.join(self.out, "search_provenance.json")))["fragpipe_database"]
        self.assertEqual(db["path"], fasta)
        self.assertFalse(db["decoys_built_in_job"])
        self.assertNotIn("database --custom", open(self.job).read())

    def test_the_decoy_tag_source(self):
        """L11 note: rev_ when the workflow names none is FragPipe's default, said so"""
        wf = os.path.join(self.d, "x.workflow")
        _write(wf, "database.decoy-tag=decoy_\n")
        self.assertEqual(fragpipe_db.decoy_tag(wf), ("decoy_", fragpipe_db.TAG_FROM_WORKFLOW))
        _write(wf, "quantitation.run-label-free-quant=false\n")
        tag, src = fragpipe_db.decoy_tag(wf)
        self.assertEqual(tag, "rev_")
        self.assertTrue(src.startswith("default"), src)

    def test_the_sidecar_sha256_is_used_not_a_login_node_hash(self):
        fasta = _write(os.path.join(self.d, "sc.fasta"), COMPLETE)
        _write(fasta + ".meta.json", json.dumps({"fasta": fasta, "sha256": "ab" * 32}))
        p = self.generate(fasta=fasta)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        db = json.load(open(os.path.join(self.out, "search_provenance.json")))["fragpipe_database"]
        self.assertEqual(db["sha256"], "ab" * 32)
        self.assertIn("sidecar", db["sha256_source"])
        # the job measures the real file; a sidecar that disagrees is said, the measured kept
        self.assertEqual(self.run_job().returncode, 0)
        merged, _java, note = fragpipe_db.merge_job_record(self.out, db)
        self.assertEqual(merged["sha256"], hashlib.sha256(open(fasta, "rb").read()).hexdigest())
        self.assertIn("changed since its sidecar", note)

    def test_the_replay_searches_from_the_target_fasta(self):
        """H3: reproduce.sh's step 5 searches the target FASTA step 4 rebuilt, with the workflow
        that ran minus its database line, from the bundle folder (relative paths) -- it runs"""
        bundle_dir = os.path.join(self.d, "bundle")
        _write(os.path.join(bundle_dir, "search.fasta"), TARGETS)
        wf = without(lfq(), "database.db-path")
        _write(os.path.join(bundle_dir, "inputs", "fragpipe_input.workflow"), wf)
        tools = tools_json(self.d, self.fp, self.fp_install)
        bundle = _write(os.path.join(bundle_dir, "wf", "workflow.manifest.json"), json.dumps(
            {"engine": {"name": "fragpipe"}, "acquisition": "DDA"}))
        from job_env import job_env
        env = job_env(self.d, FRAGPIPE_MONO_BIN=self.mono)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"), "--tools",
                            tools, "--bundle", "./wf/workflow.manifest.json", "--params",
                            "./inputs/fragpipe_input.workflow", "--fasta", "./search.fasta",
                            "--out", "./search_out", "--files", *self.files[:2], "--engine",
                            "fragpipe", "--sbatch", "./search_job.sh"],
                           capture_output=True, text=True, cwd=bundle_dir, env=env, timeout=120)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.job, self.out = (os.path.join(bundle_dir, "search_job.sh"),
                              os.path.join(bundle_dir, "search_out"))
        r = self.run_job(FAKE_FP_LFQ="1")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("ALL JOBS DONE", r.stdout)
        self.assertEqual(listing(self.raw), self.before)


class DiaStaged(StageHarness):
    """LAST item: FragPipe's DIA route on .raw/mzML stages through raw_stage; a '--' in a staged
    FILE name still reaches DIA-NN and is refused; the job waits for DIA-NN's report."""

    DIA_WF = ("workflow.misc.save-sdrf=false\ndiann.run-dia-nn=true\n"
              "quantitation.run-label-free-quant=false\ndatabase.decoy-tag=rev_\n")

    def test_a_dashed_raw_file_name_is_refused(self):
        dashed = _write(os.path.join(self.raw, "x--y.raw"), "thermo\n")
        self.before = listing(self.raw)
        p = self.generate(files=[dashed, self.files[0]], acquisition="DIA", workflow=self.DIA_WF)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("REFUSING: this FragPipe workflow runs DIA-NN", p.stderr)
        self.assertFalse(os.path.exists(self.out))

    def test_the_job_waits_for_diann_and_then_empties_the_stage(self):
        mz = _write(os.path.join(self.raw, "S4.mzML"), "")
        p = self.generate(files=[self.files[0], mz], acquisition="DIA", workflow=self.DIA_WF)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        r = self.run_job()                                  # the fake wrote no DIA-NN report
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("dia-quant-output/report.parquet or dia-quant-output/report.tsv",
                      r.stderr)
        self.assertIn("S1-a.raw", os.listdir(self.stage()))
        _write(os.path.join(self.out, "dia-quant-output", "report.tsv"), "Run\n")
        r = self.run_job()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(sorted(os.listdir(self.stage())), ["S4.mzML"])   # the link stays
        self.assertEqual(listing(self.raw), self.before + ["S4.mzML"])


class DiaRunNames(unittest.TestCase):
    """A staged DIA .raw: FragPipe hands DIA-NN MSFragger's <run>_uncalibrated.mzML, so its
    report's runs are '<run>_uncalibrated' (HIVE, 2026-11-04); the adapter names them as the
    inputs, from the search's own manifest"""

    def setUp(self):
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            self.skipTest("pyarrow")
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        _write(os.path.join(self.d, "fragpipe.fp-manifest"),
               f"{self.d}/raw_stage/S1-a.raw\t\t\tDIA\n{self.d}/raw_stage/S2.mzML\t\t\tDIA\n")

    def test_the_tsv_report(self):
        import pyarrow.parquet as pq
        cols = ["Run", "Protein.Group", "PG.MaxLFQ"] + list(run_search.FDR_REQUIRED)
        rows = [["S1-a_uncalibrated", "P1", "10"] + ["0.001"] * len(run_search.FDR_REQUIRED),
                ["S2", "P1", "12"] + ["0.001"] * len(run_search.FDR_REQUIRED)]
        _write(os.path.join(self.d, "dia-quant-output", "report.tsv"),
               "\t".join(cols) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
        rep = run_search.adapt_fragpipe_dia(self.d)
        self.assertEqual(sorted(set(pq.read_table(rep).column("Run").to_pylist())),
                         ["S1-a", "S2"])

    def test_the_parquet_report_is_copied_with_the_names(self):
        import pyarrow as pa
        import pyarrow.parquet as pq
        os.makedirs(os.path.join(self.d, "dia-quant-output"))
        src = os.path.join(self.d, "dia-quant-output", "report.parquet")
        pq.write_table(pa.table({"Run": ["S1-a_uncalibrated", "S2"], "Protein.Group": ["P1", "P1"],
                                 "PG.MaxLFQ": [1.0, 2.0]}), src)
        rep = run_search.adapt_fragpipe_dia(self.d)
        self.assertEqual(rep, os.path.join(self.d, "report.parquet"))
        self.assertEqual(pq.read_table(rep).column("Run").to_pylist(), ["S1-a", "S2"])
        self.assertEqual(pq.read_table(src).column("Run").to_pylist()[0], "S1-a_uncalibrated")


class ContaminantDigest(unittest.TestCase):
    """L7: warn_contaminant_digest skipped FragPipe; it reads the shared digest
    (cleavage_rules.fragpipe_digest)"""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        import estimate_params
        self.fasta = _write(os.path.join(self.d, "db.fasta"), TARGETS)
        _write(self.fasta + ".meta.json", json.dumps(
            {"contaminant_digest": estimate_params.DIANN_DIGEST, "min_unique_peptides": 2}))

    def test_a_non_specific_fragpipe_search_is_warned(self):
        msg = run_search.warn_contaminant_digest("fragpipe", FP_WORKFLOW, self.fasta)
        self.assertIsNotNone(msg)
        self.assertIn("WARNING", msg)
        self.assertIn("(none: a non-specific digest)", msg)

    def test_a_tryptic_rule_in_the_shared_syntax(self):
        wf = _write(os.path.join(self.d, "t.workflow"), tryptic())
        self.assertEqual(run_search.search_digest("fragpipe", wf)["cut"], "K*,R*,!*P")

    def test_the_shared_dda_digest_is_one_info_line(self):
        """Ruling 57.1: FragPipe's DDA preset digests as Sage does (the one shared digest), so
        that known difference is an INFO, not a WARNING"""
        wf = _write(os.path.join(self.d, "t.workflow"), tryptic())
        msg = run_search.warn_contaminant_digest("fragpipe", wf, self.fasta)
        self.assertIn("INFO: The FragPipe DDA preset's digest (the shared DDA digest", msg)
        self.assertIn("expected for every FragPipe DDA preset search", msg)

    def test_a_make_presets_default_preset_is_the_known_digest(self):
        """the preset make_presets.py writes by default (no --enzyme) states that digest"""
        import estimate_params
        wf = os.path.join(self.d, "mp.workflow")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_presets.py"), "--engine",
                            "fragpipe", "--instrument", "timsTOF Pro 2", "--acquisition", "DDA",
                            "--out", wf], capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        got = run_search.search_digest("fragpipe", wf)
        got.pop("met_excision", None)
        self.assertEqual(got, run_search.sage_digest(estimate_params.DDA_DIGEST))
        msg = run_search.warn_contaminant_digest("fragpipe", wf, self.fasta)
        self.assertNotIn("WARNING", msg)

    def test_an_n_terminal_sage_rule_in_the_one_cut_syntax(self):
        """sage_digest writes cleavage_rules.diann_cut's syntax: an N-terminal rule (Sage aspn,
        c_terminal false) cuts BEFORE D (*D), not after it"""
        self.assertEqual(run_search.sage_digest({"cleave_at": "D", "restrict": None,
                                                 "c_terminal": False})["cut"], "*D")
        self.assertEqual(run_search.sage_digest({"cleave_at": "KR", "restrict": "P"})["cut"],
                         "K*,R*,!*P")


class EnzymeDropdown(StageHarness):
    """C0 fact 1: a dropdown FragPipe would load over the workflow's rule. Ruling 104: a user's or
    vendor's workflow is never refused for it -- it is said (a NOTE), and the search runs as
    FragPipe loads it; only the skill's own presets are refused on it (make_presets.py)"""

    def test_one_that_changes_the_search_is_said_not_refused(self):
        wf = open(FP_WORKFLOW, newline="").read()
        wf = (tryptic(wf).replace("msfragger.search_enzyme_cut_1=KR",
                                  "msfragger.search_enzyme_cut_1=D")
              .replace("msfragger.search_enzyme_sense_1=C", "msfragger.search_enzyme_sense_1=N")
              .replace("msfragger.search_enzyme_nocut_1=P", "msfragger.search_enzyme_nocut_1="))
        wf += "msfragger.misc.fragger.enzyme-dropdown-1=trypsin\r\n"
        p = self.generate(workflow=wf)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertNotIn("REFUSED", p.stderr)
        self.assertIn("NOTE", p.stderr)
        self.assertIn("enzyme-dropdown-1=trypsin", p.stderr)
        self.assertIn("fragger.params", p.stderr)

    def test_h5s_hand_made_workflow_is_not_refused(self):
        """H5's dropdown `custom` beside its FLMWY / no-cut P rule (it ran as argc): said, and
        the search is generated"""
        wf = open(FP_WORKFLOW, newline="").read()
        wf = (tryptic(wf).replace("msfragger.search_enzyme_cut_1=KR",
                                  "msfragger.search_enzyme_cut_1=FLMWY"))
        wf += "msfragger.misc.fragger.enzyme-dropdown-1=custom\r\n"
        p = self.generate(workflow=wf)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertNotIn("REFUSED", p.stderr)
        self.assertIn("argc", p.stderr)

    def test_no_dropdown_no_note(self):
        p = self.generate(workflow=tryptic())
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertNotIn("enzyme-dropdown", p.stderr)


class InlineHeap(unittest.TestCase):
    """an inline search inside an allocation (reproduce.sh's step 5 in a job) sizes FragPipe's
    heap from the allocation, not the node (H5 gap 4)"""

    def test_from_the_allocation(self):
        from unittest import mock
        with mock.patch.dict(os.environ, {"SLURM_MEM_PER_NODE": str(128 * 1024)}):
            self.assertEqual(run_search.fragpipe_inline_ram_gb(), run_search.fragpipe_ram_gb(128))
        with mock.patch.dict(os.environ, {}, clear=False):
            for k in ("SLURM_MEM_PER_NODE", "SLURM_MEM_PER_CPU", "SLURM_CPUS_ON_NODE"):
                os.environ.pop(k, None)
            self.assertIsNone(run_search.fragpipe_inline_ram_gb())

    def test_from_a_mem_per_cpu_allocation(self):
        """review L4: under --mem-per-cpu SLURM sets no SLURM_MEM_PER_NODE; the allocation is the
        per-CPU memory times the CPUs on this node -- not the node's, which FragPipe would take"""
        from unittest import mock
        with mock.patch.dict(os.environ, {"SLURM_MEM_PER_CPU": "4096", "SLURM_CPUS_ON_NODE": "16"}):
            os.environ.pop("SLURM_MEM_PER_NODE", None)
            self.assertEqual(run_search.fragpipe_inline_ram_gb(), run_search.fragpipe_ram_gb(64))
        with mock.patch.dict(os.environ, {"SLURM_MEM_PER_CPU": "4096"}):
            for k in ("SLURM_MEM_PER_NODE", "SLURM_CPUS_ON_NODE"):
                os.environ.pop(k, None)
            self.assertIsNone(run_search.fragpipe_inline_ram_gb())     # no CPU count: unknown
        # --mem wins when both are present
        with mock.patch.dict(os.environ, {"SLURM_MEM_PER_NODE": str(32 * 1024),
                                          "SLURM_MEM_PER_CPU": "4096", "SLURM_CPUS_ON_NODE": "16"}):
            self.assertEqual(run_search.fragpipe_inline_ram_gb(), run_search.fragpipe_ram_gb(32))


class RunCountTime(unittest.TestCase):
    """M10: SEARCH_WALL_HOURS' 12 h would cut a 330-run FragPipe search"""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.wf = _write(os.path.join(self.d, "t.workflow"), tryptic())

    def test_the_time_follows_the_runs_as_a_default_said_so(self):
        z = run_search.single_job_sizing("fragpipe", self.wf, n_runs=330, threads=32)
        self.assertEqual(z["hours"], 44)                        # 330 x 4 min x 2 / 60
        self.assertTrue(z["hours_source"].startswith("default for 330 runs"), z)
        run_search.EMITTED_JOBS[:] = [{"path": "j", "mem": "64G", "hours": 44}]
        try:
            rec = run_search.job_size_record(z, run_search.EMITTED_JOBS)
        finally:
            del run_search.EMITTED_JOBS[:]
        self.assertEqual(rec["job_time"], "44h (default for 330 runs)")

    def test_few_runs_keep_the_route_default(self):
        z = run_search.single_job_sizing("fragpipe", self.wf, n_runs=2, threads=32)
        self.assertIsNone(z["hours"])
        self.assertEqual(z["hours_source"], "default")

    def test_fewer_cpus_take_longer(self):
        self.assertEqual(run_search.fragpipe_run_hours(330, 8), 88)

    def test_the_users_time_wins(self):
        z = run_search.single_job_sizing("fragpipe", self.wf, hours=6, n_runs=330, threads=32)
        self.assertEqual((z["hours"], z["hours_source"]), (6, "user (--time)"))

    def test_a_broad_digest_takes_the_larger(self):
        z = run_search.single_job_sizing("fragpipe", FP_WORKFLOW, n_runs=330, threads=32)
        self.assertEqual(z["hours"], 44)
        z = run_search.single_job_sizing("fragpipe", FP_WORKFLOW, n_runs=20, threads=32)
        self.assertEqual(z["hours"], run_search.BROAD_DIGEST_HOURS)

    def test_a_generated_job_asks_for_it(self):
        h = StageHarness("setUp")
        h.setUp()
        self.addCleanup(shutil.rmtree, h.d, True)
        files = [_write(os.path.join(h.raw, f"r{i:03d}.mzML"), "") for i in range(100)]
        p = h.generate(files=files, workflow=tryptic(lfq()))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        # 100 runs x 4 min x 2 headroom, x4 at the harness's 4 CPUs (measured at 16): 54 h
        self.assertEqual(run_search.fragpipe_run_hours(100, 4), 54)
        self.assertIn("#SBATCH --time=54:00:00", open(h.job).read())
        prov = json.load(open(os.path.join(h.out, "search_provenance.json")))
        self.assertEqual(prov["job_time"], "54h (default for 100 runs)")


class EngineChoiceRecord(StageHarness):
    def test_the_selection_is_copied_verbatim(self):
        sel = {"engine": "fragpipe", "default": "fragpipe", "how": "default", "sentence": None,
               "decided_by": "engine_choice.py"}
        p = self.generate(bundle={"engine": {"name": "fragpipe", "selection": sel},
                                  "acquisition": "DDA"})
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        prov = json.load(open(os.path.join(self.out, "search_provenance.json")))
        self.assertEqual(prov["engine_selection"], sel)

    def test_fragpipe_not_ready_here_is_refused(self):
        fp = fake_install(os.path.join(self.d, "notready"))
        fp["fragpipe_tools"].update(ready=False, why="no Java here")
        fp["fragpipe_tools"]["java"] = {"version": None, "major": None, "ok": False}
        p = self.generate(tools=tools_json(self.d, self.fp, fp, name="notready.json"))
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("FragPipe cannot run here", p.stderr)
        self.assertFalse(os.path.exists(self.out))

    def test_no_engine_named_asks_the_one_rule(self):
        class A:
            engine = None
        got = run_search.pick_engine(A(), {"acquisition": "DDA", "instruments": ["timsTOF HT"],
                                           "instrument_class": "timstof"},
                                     {"fragpipe_tools": self.fp_install["fragpipe_tools"]})
        import engine_choice
        want = engine_choice.choose("DDA", "timstof",
                                    fp_tools=self.fp_install["fragpipe_tools"],
                                    cls="timstof")[0]["engine"]
        self.assertEqual(got, want)
        self.assertEqual(got, "fragpipe")       # timsTOF DDA, FragPipe ready: the default

    def test_no_instrument_class_reads_the_route_tables_unclassified_row(self):
        """a bundle with no instrument class: ROUTES[(acq, "other")] -- DIA-NN / Sage, as 2.14.2"""
        class A:
            engine = None
        got = {acq: run_search.pick_engine(A(), {"acquisition": acq}) for acq in ("DIA", "DDA", "")}
        self.assertEqual(got, {"DIA": "diann", "DDA": "sage", "": "sage"})


if __name__ == "__main__":
    unittest.main()
