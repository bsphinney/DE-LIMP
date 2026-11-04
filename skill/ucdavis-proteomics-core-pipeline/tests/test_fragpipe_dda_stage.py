#!/usr/bin/env python3
"""
A FragPipe DDA search never writes into the input raw folder (2.14.2, L11 N1, N4).

A 38-run hydrolysate search on Thermo .raw: run_search.py --engine fragpipe wrote the raw paths
into fragpipe.fp-manifest as given, so MSFragger wrote <run>_uncalibrated.mzML and .mzBIN (~1 GB a
run) next to each .raw -- unasked writes into a collaborator's raw folder, and a read-only share
would fail the search. The job also had no mono on PATH (FragPipe reads .raw through
tools/ext/thermo/BatmassIoThermoServer.exe, a .NET Framework program), and the manifest's
experiment column was empty, so IonQuant wrote no combined_protein.tsv (H5 gap 1 too). Staff
copied the raws into <search>/raw_stage by hand (a symlink to a .raw fails in MSFragger:
"Instrument index not available for requested device"), and later had to move that 50 GB folder
out of the session before finalize --zip.

Now: the inputs are staged in <out>/raw_stage (a .raw copied in the job, a .d or mzML linked as
diatracer_stage.py links a .d), the manifest names the staged paths with one experiment per run,
the job puts mono on PATH (or the route refuses .raw before anything is written), the job removes
its staged copies and MSFragger's uncalibrated files after a search that finished and keeps them
after one that did not -- and raw_stage is scratch (scratch_files.py): never zipped, catalogued
or delivered.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import run_search  # noqa: E402
import scratch_files  # noqa: E402
from job_env import job_env  # noqa: E402

FP_WORKFLOW = os.path.join(HERE, "fixtures", "fragpipe24_nonspecific_peptidome",
                           "Nonspecific-peptidome.workflow")

# A fake FragPipe: reads --manifest, writes what MSFragger writes beside each input, each
# experiment's identification tables under FragPipe's own experiment name ('-' -> '_', and its
# fragpipe-files.fp-manifest, as FragPipe 24.0 does), IonQuant's combined_protein.tsv only with
# FAKE_FP_LFQ=1, and what FragPipe prints at the end. FAKE_FP_RC / FAKE_FP_DONE / FAKE_FP_IDS=0
# shape a failure.
FAKE_FRAGPIPE = r"""#!/bin/bash
man=""; wd=""
while [ $# -gt 0 ]; do case "$1" in --manifest) man="$2"; shift;; --workdir) wd="$2"; shift;;
  esac; shift; done
: > "$wd/fragpipe-files.fp-manifest"; : > "$wd/.fake_fragpipe_ran"
while IFS=$'\t' read -r path exp rep kind; do
  stem="${path%.*}"
  : > "${stem}_uncalibrated.mzML"; : > "${stem}_uncalibrated.mzBIN"; : > "${stem}.mzBIN"
  fexp="${exp//-/_}"
  printf '%s\t%s\t%s\t%s\n' "$path" "$fexp" "$rep" "$kind" >> "$wd/fragpipe-files.fp-manifest"
  if [ "${FAKE_FP_IDS:-1}" = 1 ]; then
    # one table set per experiment; with the column empty, ONE set for all runs (FragPipe's
    # pooled filter), so rows are appended
    t="$wd/$fexp"; mkdir -p "$t"
    run="$(basename "${path%.*}")"
    [ -f "$t/psm.tsv" ] || printf 'Spectrum\tPeptide\tCharge\tProbability\tAssigned Modifications\tProtein ID\tMapped Proteins\n' > "$t/psm.tsv"
    printf '%s.00001.00001.2\tPEPTIDEK\t2\t0.99\t\tP02768\t\n' "$run" >> "$t/psm.tsv"
    echo "protein" > "$t/protein.tsv"
    [ -f "$t/peptide.tsv" ] || printf 'Peptide\tCharges\tProbability\tAssigned Modifications\tProtein ID\tMapped Proteins\nPEPTIDEK\t2\t0.99\t\tP02768\t\n' > "$t/peptide.tsv"
  fi
done < "$man"
[ "${FAKE_FP_LFQ:-0}" = 1 ] && printf 'Protein\tProtein ID\n' > "$wd/combined_protein.tsv"
[ "${FAKE_FP_DONE:-1}" = 1 ] && echo "=====ALL JOBS DONE IN 0.1 MINUTES====="
exit "${FAKE_FP_RC:-0}"
"""


def _write(path, text, mode=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    if mode:
        os.chmod(path, mode)
    return path


def listing(d):
    return sorted(os.path.relpath(os.path.join(dp, f), d)
                  for dp, dns, fns in os.walk(d) for f in fns + dns)


class StageHarness(unittest.TestCase):
    """The fixture and helpers, no tests (test_fragpipe_ids_only.py builds on it)."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.raw = os.path.join(self.d, "raw data")              # a folder we must never write
        self.files = [_write(os.path.join(self.raw, "S1-a.raw"), "thermo one\n"),
                      _write(os.path.join(self.raw, "S2-b.raw"), "thermo two\n")]
        os.makedirs(os.path.join(self.raw, "S3.d"))
        _write(os.path.join(self.raw, "S3.d", "analysis.tdf"), "")
        self.files.append(os.path.join(self.raw, "S3.d"))
        self.before = listing(self.raw)
        self.mono = os.path.join(self.d, "mono_bin")
        _write(os.path.join(self.mono, "mono"), "#!/bin/sh\n", 0o755)
        self.fp = _write(os.path.join(self.d, "bin", "fragpipe"), FAKE_FRAGPIPE, 0o755)
        self.out = os.path.join(self.d, "out")
        self.job = os.path.join(self.d, "fp_job.sh")

    def generate(self, *extra, mono=True, files=None, workflow=None, acquisition="DDA",
                 fasta=None):
        """`workflow`: the .workflow text to search with (default: the shipped
        Nonspecific-peptidome.workflow, which runs no quantification). `fasta`: the database
        given as --fasta (default: a one-entry db.fasta written here)."""
        fasta = fasta or _write(os.path.join(self.d, "db.fasta"), ">sp|P1|X\nPEPTIDER\n")
        params = os.path.join(self.d, "search.workflow")
        if workflow is None:
            shutil.copy(FP_WORKFLOW, params)
        else:
            with open(params, "w", newline="") as fh:
                fh.write(workflow)
        tools = _write(os.path.join(self.d, "tools.json"), json.dumps({"fragpipe": self.fp}))
        bundle = _write(os.path.join(self.d, "bundle.json"), json.dumps(
            {"engine": {"name": "fragpipe"}, "acquisition": acquisition}))
        env = job_env(self.d, PATH=os.pathsep.join([os.path.dirname(sys.executable),
                                                    "/usr/bin", "/bin"]),
                      FRAGPIPE_MONO_BIN=self.mono if mono else os.path.join(self.d, "no_mono"))
        env.pop("FRAGPIPE_TOOLS_FOLDER", None)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"),
                            "--tools", tools, "--bundle", bundle, "--params", params,
                            "--fasta", fasta, "--out", self.out,
                            "--files", *(files or self.files),
                            "--engine", "fragpipe", "--threads", "4", "--no-notify", "--no-fran",
                            # the fixture's .d is a folder, not a readable timsTOF run
                            "--allow-damaged-tdf", "--sbatch", self.job, *extra],
                           capture_output=True, text=True, env=env, cwd=self.d, timeout=240)
        return p

    def run_job(self, **extra):
        env = job_env(self.d, PATH=os.environ["PATH"], **extra)
        return subprocess.run(["bash", self.job], cwd=self.out, capture_output=True, text=True,
                              env=env, timeout=240)

    def stage(self):
        return os.path.join(self.out, scratch_files.FRAGPIPE_STAGE_DIR)

    def manifest(self):
        with open(os.path.join(self.out, "fragpipe.fp-manifest")) as fh:
            return [ln.rstrip("\n").split("\t") for ln in fh]



class Stage(StageHarness):
    def test_the_manifest_points_at_the_stage_with_one_experiment_per_run(self):
        """one experiment per run when the workflow quantifies (IonQuant needs them); the
        shipped IDs-only workflow keeps FragPipe's one pooled experiment (follow-up 3)"""
        lfq = open(FP_WORKFLOW, newline="").read().replace(
            "quantitation.run-label-free-quant=false", "quantitation.run-label-free-quant=true")
        p = self.generate(workflow=lfq)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        rows = self.manifest()
        self.assertEqual([r[1] for r in rows],
                         [run_search.sage_run_name(f) for f in self.files])     # S1-a, S2-b, S3
        p = self.generate()                                       # the shipped, IDs only
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        rows = self.manifest()
        self.assertEqual([os.path.dirname(r[0]) for r in rows], [self.stage()] * 3)
        self.assertEqual([os.path.basename(r[0]) for r in rows], ["S1-a.raw", "S2-b.raw", "S3.d"])
        self.assertEqual([r[1] for r in rows], ["", "", ""])
        self.assertIn("pooled", json.load(open(os.path.join(
            self.out, "search_provenance.json")))["fragpipe_stage"]["experiment_rule"])
        self.assertEqual({r[3] for r in rows}, {"DDA"})
        self.assertTrue(os.path.islink(os.path.join(self.stage(), "S3.d")))  # a .d is linked
        self.assertFalse(os.path.exists(os.path.join(self.stage(), "S1-a.raw")))   # copied in the job
        self.assertEqual(listing(self.raw), self.before)
        prov = json.load(open(os.path.join(self.out, "search_provenance.json")))
        st = prov["fragpipe_stage"]
        self.assertEqual(st["dir"], self.stage())
        self.assertEqual(sorted(os.path.basename(x["src"]) for x in st["copied"]),
                         ["S1-a.raw", "S2-b.raw"])
        self.assertEqual([os.path.basename(x["src"]) for x in st["linked"]], ["S3.d"])
        self.assertEqual(prov["fragpipe_mono"]["bin"], self.mono)

    def test_the_job_copies_then_removes_after_a_search_that_finished(self):
        self.assertEqual(self.generate().returncode, 0)
        text = open(self.job).read()
        self.assertIn(f'export PATH={self.mono}:"$PATH"', text)          # mono for the .raw
        r = self.run_job()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        left = sorted(os.listdir(self.stage()))
        self.assertEqual(left, ["S3.d"], left)          # copies and MSFragger's files removed
        self.assertIn("removed the staged copies", r.stdout + r.stderr)
        self.assertEqual(listing(self.raw), self.before)                  # nothing beside the raws

    def test_a_search_that_failed_keeps_the_staged_copies(self):
        self.assertEqual(self.generate().returncode, 0)
        r = self.run_job(FAKE_FP_RC="3")
        self.assertEqual(r.returncode, 3, r.stdout + r.stderr)
        self.assertIn("S1-a.raw", os.listdir(self.stage()))
        self.assertIn("KEPT for a re-run", r.stdout + r.stderr)
        self.assertEqual(listing(self.raw), self.before)

    def test_exit_0_without_all_jobs_done_is_not_a_finished_search(self):
        self.assertEqual(self.generate().returncode, 0)
        r = self.run_job(FAKE_FP_DONE="0")
        self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("S1-a.raw", os.listdir(self.stage()))

    def test_a_rerun_reuses_a_copy_already_staged(self):
        self.assertEqual(self.generate().returncode, 0)
        self.assertEqual(self.run_job(FAKE_FP_RC="1").returncode, 1)
        r = self.run_job()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("already staged", r.stdout + r.stderr)

    def test_without_mono_a_raw_is_refused_before_anything_is_written(self):
        p = self.generate(mono=False)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("mono", p.stderr)
        self.assertIn("mzML", p.stderr)
        self.assertIn("ThermoRawFileParser", p.stderr)
        self.assertFalse(os.path.exists(self.job))
        self.assertFalse(os.path.exists(self.out))
        self.assertEqual(listing(self.raw), self.before)

    def test_mzml_and_d_need_no_mono(self):
        mz = _write(os.path.join(self.raw, "S4.mzML"), "")
        p = self.generate(mono=False, files=[mz, os.path.join(self.raw, "S3.d")])
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        text = open(self.job).read()
        self.assertNotIn("no_mono", text)
        self.assertNotIn("[fragpipe] mono", text)

    def test_two_inputs_with_one_file_name_are_refused(self):
        other = _write(os.path.join(self.d, "elsewhere", "S1-a.raw"), "another\n")
        p = self.generate(files=[self.files[0], other])
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("S1-a.raw", p.stderr)
        self.assertFalse(os.path.exists(self.job))


class StageSafety(StageHarness):
    """Review F2, F10, F17: free space is checked before the copies (at generation, and again in
    the job), a failed copy removes its partial file and says what is kept; a link left in the
    stage is replaced by a copy; an input given from inside the stage is refused (the job empties
    the stage); after a finished search every regular file in the stage goes, links stay."""

    def test_not_enough_space_is_refused_before_anything_is_written(self):
        import run_search
        from unittest import mock
        out = os.path.join(self.d, "out_nospace")
        with mock.patch("shutil.disk_usage", return_value=(10, 9, 1)), \
                self.assertRaises(SystemExit) as cm:
            run_search.refuse_fragpipe_stage({"acquisition": "DDA"}, self.files, out)
        msg = str(cm.exception.code)
        for words in ("REFUSED", "GB free", "needs", "--out", "mzML"):
            self.assertIn(words, msg)
        self.assertFalse(os.path.exists(out))

    def test_the_job_checks_space_again_before_copying(self):
        self.assertEqual(self.generate().returncode, 0)
        text = open(self.job).read()
        self.assertIn("df -Pk", text)
        self.assertLess(text.index("df -Pk"), text.index("_stage_copy '"))

    def test_a_failed_copy_removes_its_partial_file_and_says_so(self):
        self.assertEqual(self.generate().returncode, 0)
        os.chmod(self.files[1], 0)                               # S2-b.raw cannot be read
        self.addCleanup(os.chmod, self.files[1], 0o644)
        r = self.run_job()
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("could not be copied", r.stdout + r.stderr)
        self.assertNotIn("S2-b.raw.part", os.listdir(self.stage()))
        self.assertFalse(os.path.exists(os.path.join(self.out, ".fake_fragpipe_ran")))

    def test_a_link_left_in_the_stage_is_replaced_by_a_copy(self):
        self.assertEqual(self.generate().returncode, 0)
        link = os.path.join(self.stage(), "S1-a.raw")
        os.symlink(self.files[0], link)                          # a hand attempt's link
        r = self.run_job(FAKE_FP_RC="1")                         # keep the stage to look at
        self.assertEqual(r.returncode, 1, r.stdout + r.stderr)
        self.assertFalse(os.path.islink(link))
        self.assertTrue(os.path.isfile(link))

    def test_an_input_inside_the_stage_is_refused(self):
        inside = _write(os.path.join(self.out, "raw_stage", "S9.raw"), "copied by hand\n")
        p = self.generate(files=[inside])
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("raw_stage", p.stderr)
        self.assertTrue(os.path.isfile(inside))
        self.assertFalse(os.path.exists(self.job))

    def test_after_a_finished_search_every_regular_file_in_the_stage_goes(self):
        self.assertEqual(self.generate().returncode, 0)
        for extra in ("S1-a_calibrated.mzML", "S1-a.quantindex"):
            _write(os.path.join(self.stage(), extra), "side file\n")
        r = self.run_job()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(sorted(os.listdir(self.stage())), ["S3.d"])     # the link stays


class ExpectedOutputsUnderPipefail(unittest.TestCase):
    """C2 errors_ux M3: the job checked each expected table with `find ... | grep -q .` under
    `set -euo pipefail`. With ~30 or more runs (one table set per run since 2.14.2) find prints
    more than a pipe block, grep -q closes the pipe after the first line, find gets SIGPIPE (141),
    pipefail fails the test, and a FINISHED search was reported unfinished with its stage kept.
    Now `[ -n "$(find ... -print -quit)" ]`: find stops by itself."""

    def job(self, out):
        stage = {"dir": os.path.join(out, "raw_stage"), "copied": []}
        return run_search.fragpipe_stage_job(stage, "true", os.path.join(out, "c.log"),
                                             ["psm.tsv"], out)

    def test_the_check_uses_print_quit_and_no_pipe_to_grep(self):
        with tempfile.TemporaryDirectory() as out:
            text = self.job(out)
        lines = [ln for ln in text.splitlines() if "-name psm.tsv" in ln]
        self.assertEqual(len(lines), 1, text)
        self.assertIn("-print -quit", lines[0])
        self.assertNotIn("grep -q", lines[0])

    def test_many_runs_tables_are_found_under_pipefail(self):
        with tempfile.TemporaryDirectory() as out:
            for i in range(1500):                      # ~300 KB of paths: far past a pipe block
                _write(os.path.join(out, f"experiment_{'x' * 150}_{i:04d}", "psm.tsv"), "t\n")
            line = [ln for ln in self.job(out).splitlines() if "-name psm.tsv" in ln][0]
            r = subprocess.run(["bash", "-c", "set -euo pipefail\n_fp_missing=''\n" + line +
                                '\nprintf "missing=[%s]" "$_fp_missing"'],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(r.stdout, "missing=[]", r.stderr)

    def test_a_table_never_written_is_still_missing(self):
        with tempfile.TemporaryDirectory() as out:
            _write(os.path.join(out, "S1", "psm.tsv"), "")             # empty: not written
            line = [ln for ln in self.job(out).splitlines() if "-name psm.tsv" in ln][0]
            r = subprocess.run(["bash", "-c", "set -euo pipefail\n_fp_missing=''\n" + line +
                                '\nprintf "missing=[%s]" "$_fp_missing"'],
                               capture_output=True, text=True, timeout=60)
            self.assertEqual(r.stdout, "missing=[ psm.tsv]", r.stderr)


class ManifestRowsStayWhole(StageHarness):
    """C2 security n3: FragPipe's manifest is one tab-separated line a run; a tab or a line break
    in an input's path split its row (only --out was checked). Refused before anything is
    written."""

    def test_a_tab_or_a_line_break_in_an_input_path_is_refused(self):
        for name in ("S9\tx.raw", "S9\nx.raw"):
            with self.subTest(name=repr(name)):
                bad = _write(os.path.join(self.d, "odd", name), "thermo\n")
                p = self.generate(files=[bad])
                self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
                self.assertIn("tab or a line break", p.stderr)
                self.assertFalse(os.path.exists(self.job))


class TheDecoyFasFindsItsTargetSidecar(StageHarness):
    """H7b-1 (HIVE, 2026-11-03): the FragPipe database fix asks for --fasta to be Philosopher's
    target+decoy .fas (MSFragger refuses a target-only FASTA), which has no .meta.json. The search
    record then said keratin_sample_source null, and nothing that reads the sidecar from the
    record's `fasta` found one; the Methods were right only because --fasta-meta was given by
    hand. Now the target FASTA's sidecar is found (fetch_fasta.sidecar_for: <name>.meta.json
    beside the .fas or in the session's input/, only when it describes the .fas's target
    entries) and recorded (search_provenance.json `fasta_meta`), or the search is refused with
    the --fasta-meta hint -- never null facts recorded silently."""

    FAS = "2026-11-03-decoys-search.fasta.fas"

    def database(self, n_entries=1, sidecar=True):
        db = os.path.join(self.d, "db")
        _write(os.path.join(db, "search.fasta"), ">sp|P1|X\nPEPTIDER\n")
        side = os.path.join(db, "search.fasta.meta.json")
        if sidecar:
            _write(side, json.dumps({"fasta": os.path.join(db, "search.fasta"),
                                     "n_entries": n_entries, "keratin_sample": False,
                                     "keratin_sample_source": "user"}))
        fas = _write(os.path.join(db, self.FAS),
                     ">sp|P1|X\nPEPTIDER\n>rev_sp|P1|X\nREDITPEP\n")
        return fas, side

    def prov(self):
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            return json.load(fh)

    def test_the_targets_sidecar_beside_the_fas_is_found_and_recorded(self):
        fas, side = self.database()
        p = self.generate(fasta=fas)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        prov = self.prov()
        self.assertEqual(prov["fasta"], fas)                        # the database searched
        self.assertEqual(prov["fasta_meta"]["path"], side)
        self.assertIn("Philosopher's target+decoy database of search.fasta",
                      prov["fasta_meta"]["source"])
        self.assertEqual(prov["keratin_sample_source"], "user")     # not null (H7b-1)
        self.assertIs(prov["keratin_sample_check"]["sidecar_keratin_sample"], False)

    def test_no_sidecar_for_the_fas_is_refused_with_the_hint(self):
        fas, _ = self.database(sidecar=False)
        p = self.generate(fasta=fas)
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        for words in ("REFUSED", "nothing was converted, written or submitted",
                      "no sidecar describing its 1 target entries", "--fasta-meta",
                      "--fasta-meta none"):
            self.assertIn(words, p.stderr)
        self.assertFalse(os.path.exists(self.job))

    def test_a_sidecar_that_does_not_describe_its_targets_is_not_taken(self):
        fas, _ = self.database(n_entries=5)
        p = self.generate(fasta=fas)
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("describes 5 entries, the .fas has 1 targets", p.stderr)
        self.assertFalse(os.path.exists(self.job))

    def test_fasta_meta_names_it_or_says_there_is_none(self):
        fas, side = self.database(sidecar=False)
        other = _write(os.path.join(self.d, "elsewhere", "target.fasta.meta.json"),
                       json.dumps({"n_entries": 1, "keratin_sample": False,
                                   "keratin_sample_source": "default"}))
        p = self.generate("--fasta-meta", other, fasta=fas)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual((self.prov()["fasta_meta"]["path"], self.prov()["keratin_sample_source"]),
                         (other, "default"))
        p = self.generate("--fasta-meta", "none", fasta=fas)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue(self.prov()["fasta_meta"]["source"].startswith("--fasta-meta none"))

    def test_a_plain_fasta_without_a_sidecar_is_unchanged(self):
        """(pin) a database that is not a Philosopher .fas: no refusal, as before"""
        p = self.generate()
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIsNone(self.prov()["fasta_meta"]["path"])

    def test_a_renamed_copys_fas_takes_the_session_sidecar_by_exact_count(self):
        """Controller ruling: the reference doc's own recipe (`cp search.fasta target.fasta`,
        before 2.14.2) gives <date>-decoys-target.fasta.fas, with no target.fasta.meta.json. The
        session's input/ sidecar is taken when its n_entries equals the .fas's target entries
        EXACTLY (rule "count"), end to end through run_search.py."""
        sess = os.path.join(self.d, "S")
        side = _write(os.path.join(sess, "input", "search.fasta.meta.json"),
                      json.dumps({"n_entries": 1, "keratin_sample": False,
                                  "keratin_sample_source": "user"}))
        os.makedirs(os.path.join(sess, "output"))
        fas = _write(os.path.join(self.d, "philo", "2026-11-03-decoys-target.fasta.fas"),
                     ">sp|P1|X\nPEPTIDER\n>rev_sp|P1|X\nREDITPEP\n")
        self.out = os.path.join(sess, "output", "search_fp")
        p = self.generate(fasta=fas)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        rec = self.prov()["fasta_meta"]
        self.assertEqual((rec["path"], rec["rule"]), (side, "count"))
        self.assertIn("by rule 'count'", rec["source"])
        self.assertEqual(self.prov()["keratin_sample_source"], "user")

    def test_the_count_rule_takes_only_an_exact_and_single_match(self):
        import fetch_fasta as ff
        inp = os.path.join(self.d, "S", "input")
        fas = _write(os.path.join(self.d, "philo", "2026-11-03-decoys-target.fasta.fas"),
                     ">sp|P1|X\nA\n>sp|P2|Y\nC\n>rev_sp|P1|X\nA\n>rev_sp|P2|Y\nC\n")
        _write(os.path.join(inp, "search.fasta.meta.json"), json.dumps({"n_entries": 3}))
        rec = ff.sidecar_for(fas, look_in=[inp])
        self.assertIsNone(rec["path"])                          # 3 is not 2: not taken
        self.assertIn("describes 3 entries", rec["problem"])
        a = _write(os.path.join(inp, "mouse.fasta.meta.json"), json.dumps({"n_entries": 2}))
        self.assertEqual(ff.sidecar_for(fas, look_in=[inp])["rule"], "count")
        b = _write(os.path.join(inp, "rat.fasta.meta.json"), json.dumps({"n_entries": 2}))
        rec = ff.sidecar_for(fas, look_in=[inp])                 # two match: none is taken
        self.assertIsNone(rec["path"])
        self.assertIn("2 sidecars describe its 2 target entries, so none is taken", rec["problem"])
        for path in (a, b):
            self.assertIn(path, rec["problem"])
        # the named rule wins over the count rule
        named = _write(os.path.join(inp, "target.fasta.meta.json"), json.dumps({"n_entries": 2}))
        rec = ff.sidecar_for(fas, look_in=[inp])
        self.assertEqual((rec["path"], rec["rule"]), (named, "named"))

    def test_the_reference_recipe_keeps_the_sidecars_name(self):
        with open(os.path.join(os.path.dirname(HERE), "references", "fragpipe-diatracer.md"),
                  encoding="utf-8") as fh:
            text = fh.read()
        self.assertNotIn("cp search.fasta target.fasta", text)
        self.assertIn("cp search.fasta search.fasta.meta.json .", text)
        self.assertIn("--custom ./search.fasta --prefix rev_", text)

    def test_the_session_input_is_looked_in_and_readers_use_the_record(self):
        import contaminant_competition as cc
        import fetch_fasta as ff
        sess_input = os.path.join(self.d, "S", "input")
        side = _write(os.path.join(sess_input, "search.fasta.meta.json"),
                      json.dumps({"n_entries": 1, "keratin_sample": True}))
        fas = _write(os.path.join(self.d, "fasdir", self.FAS),
                     ">sp|P1|X\nPEPTIDER\n>rev_sp|P1|X\nREDITPEP\n")
        rec = ff.sidecar_for(fas, look_in=[sess_input])
        self.assertEqual(rec["path"], side)
        out = os.path.join(self.d, "search_out")
        _write(os.path.join(out, "search_provenance.json"),
               json.dumps({"fasta": fas, "fasta_meta": rec}))
        _report, fasta, meta, _why = cc.resolve_inputs(search_out=out)
        self.assertEqual(fasta, fas)
        self.assertIs(meta["keratin_sample"], True)


class StageNeverFollowsALink(StageHarness):
    """C2 security H2: a raw_stage made by hand as a link to the raw folder (`ln -s <raws>
    <out>/raw_stage` -- the refusal text tells staff to `ln -s` into the stage) made each .raw
    "already staged" (the dest IS the original), had MSFragger write into the raw folder, and the
    inline clean-up's os.scandir followed the link and deleted EVERY regular file in the raw
    folder. Now raw_stage must be a real folder inside --out (lstat), or it is refused before
    anything is written; the clean-up never follows a link (the job's and the inline one)."""

    def link_the_stage_to(self, target):
        os.makedirs(self.out, exist_ok=True)
        if os.path.isdir(self.stage()) and not os.path.islink(self.stage()):
            shutil.rmtree(self.stage())
        os.symlink(target, self.stage())

    def test_a_raw_stage_linked_to_the_raw_folder_is_refused_and_the_raws_survive(self):
        self.link_the_stage_to(self.raw)
        p = self.generate()
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        for words in ("REFUSED", "raw_stage", "is a link", "nothing was written"):
            self.assertIn(words, p.stderr)
        self.assertFalse(os.path.exists(self.job))
        self.assertEqual(listing(self.raw), self.before)          # every fake raw is still there
        self.assertTrue(os.path.islink(self.stage()))              # the link itself left alone

    def test_a_raw_stage_that_is_a_file_is_refused(self):
        os.makedirs(self.out)
        _write(self.stage(), "not a folder\n")
        p = self.generate()
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("not a folder", p.stderr)
        self.assertFalse(os.path.exists(self.job))

    def test_an_input_reached_through_an_alias_of_the_stage_is_refused(self):
        inside = _write(os.path.join(self.stage(), "S9.raw"), "copied by hand\n")
        alias = os.path.join(self.d, "alias_of_stage")
        os.symlink(self.stage(), alias)
        p = self.generate(files=[os.path.join(alias, "S9.raw")])
        self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertIn("raw_stage", p.stderr)
        self.assertTrue(os.path.isfile(inside))
        self.assertFalse(os.path.exists(self.job))

    def test_the_job_refuses_a_stage_made_a_link_after_it_was_generated(self):
        self.assertEqual(self.generate().returncode, 0)
        self.link_the_stage_to(self.raw)
        r = self.run_job()
        self.assertNotEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertIn("REFUSED", r.stdout + r.stderr)
        self.assertIn("nothing was copied or searched", r.stdout + r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, ".fake_fragpipe_ran")))
        self.assertEqual(listing(self.raw), self.before)

    def test_the_job_checks_the_stage_before_copying_and_before_the_clean_up(self):
        self.assertEqual(self.generate().returncode, 0)
        text = open(self.job).read()
        q = run_search.shlex.quote(self.stage())
        guard = run_search.stage_dir_is_safe_sh(q)
        self.assertEqual(text.count(guard), 2, text)
        self.assertLess(text.index(guard), text.index(f"mkdir -p {q}"))
        self.assertLess(text.index(guard), text.index("_stage_copy '"))
        self.assertIn(f"find -P {q} -maxdepth 1 -type f", text)
        self.assertLess(text.rindex(guard), text.index("find -P"))

    def test_the_jobs_stage_test_never_follows_a_link(self):
        """stage_dir_is_safe_sh: a real folder or nothing yet is safe; a link (to a folder or
        dangling) or a file is not"""
        real = os.path.join(self.d, "real_stage")
        os.makedirs(real)
        link = os.path.join(self.d, "linked_stage")
        os.symlink(self.raw, link)
        dangling = os.path.join(self.d, "dangling_stage")
        os.symlink(os.path.join(self.d, "nowhere"), dangling)
        a_file = _write(os.path.join(self.d, "file_stage"), "x\n")
        for path, safe in ((real, True), (os.path.join(self.d, "absent"), True),
                           (link, False), (dangling, False), (a_file, False)):
            test = run_search.stage_dir_is_safe_sh(run_search.shlex.quote(path))
            r = subprocess.run(["bash", "-c", f"if {test}; then echo safe; else echo unsafe; fi"],
                               capture_output=True, text=True)
            self.assertEqual(r.stdout.strip(), "safe" if safe else "unsafe", path)

    def test_the_inline_search_refuses_a_linked_stage_and_the_raws_survive(self):
        out = os.path.join(self.d, "inline_out")
        os.makedirs(out)
        os.symlink(self.raw, os.path.join(out, scratch_files.FRAGPIPE_STAGE_DIR))
        params = os.path.join(self.d, "search.workflow")
        shutil.copy(FP_WORKFLOW, params)
        fasta = _write(os.path.join(self.d, "db.fasta"),
                       ">sp|P1|X\nPEPTIDER\n>rev_sp|P1|X\nREDITPEP\n")
        with self.assertRaises(SystemExit) as cm:
            run_search.run_fragpipe(self.fp, {"acquisition": "DDA"}, params, self.files[:2],
                                    fasta, out, 2, None,
                                    mono={"bin": self.mono, "source": "test"})
        self.assertIn("is a link", str(cm.exception.code))
        self.assertFalse(os.path.exists(os.path.join(out, ".fake_fragpipe_ran")))
        self.assertEqual(listing(self.raw), self.before)

    def test_the_inline_stage_refuses_a_link_before_copying(self):
        self.link_the_stage_to(self.raw)
        stage = {"dir": self.stage(), "copied": [{"src": self.files[0],
                                                  "dest": os.path.join(self.stage(), "S1-a.raw")}]}
        log = os.path.join(self.d, "console.log")
        with self.assertRaises(SystemExit) as cm:
            run_search.fragpipe_stage_inline(stage, "echo ALL JOBS DONE", log, (), self.out)
        self.assertIn("nothing was copied or searched", str(cm.exception.code))
        self.assertFalse(os.path.exists(log))
        self.assertEqual(listing(self.raw), self.before)

    def test_the_clean_up_never_follows_a_stage_swapped_for_a_link(self):
        """a stage that became a link after the copies: NOTHING is removed, and that is said"""
        self.link_the_stage_to(self.raw)
        res = run_search.empty_fragpipe_stage(self.stage(), self.out)
        self.assertIsInstance(res, str)
        self.assertIn("is a link", res)
        self.assertEqual(listing(self.raw), self.before)

    def test_the_clean_up_removes_regular_files_and_keeps_links_and_their_targets(self):
        os.makedirs(self.stage())
        _write(os.path.join(self.stage(), "S1-a.raw"), "copy\n")
        _write(os.path.join(self.stage(), "S1-a.mzBIN"), "side\n")
        os.symlink(self.files[1], os.path.join(self.stage(), "S2-b.raw"))   # a link to a raw
        os.symlink(os.path.join(self.raw, "S3.d"), os.path.join(self.stage(), "S3.d"))
        res = run_search.empty_fragpipe_stage(self.stage(), self.out)
        self.assertEqual(sorted(res), ["S1-a.mzBIN", "S1-a.raw"])
        self.assertEqual(sorted(os.listdir(self.stage())), ["S2-b.raw", "S3.d"])
        self.assertEqual(listing(self.raw), self.before)


class MonoAndDia(StageHarness):
    """Review F12, F18: the no-mono refusal tells someone outside the Core how to proceed; the
    FragPipe DIA route on .raw/mzML -- which would write beside its inputs, with no mono check --
    is refused before anything is written, with the way out (staged properly in 2.15). DIA on .d
    keeps diatracer_stage.py (pin)."""

    def test_the_no_mono_refusal_says_how_outside_the_core(self):
        p = self.generate(mono=False)
        self.assertNotEqual(p.returncode, 0)
        for words in ("conda-forge", "FRAGPIPE_MONO_BIN", "outside the Core"):
            self.assertIn(words, p.stderr)

    def test_fragpipe_dia_on_raw_or_mzml_is_refused(self):
        mz = _write(os.path.join(self.raw, "S4.mzML"), "")
        for files in ([self.files[0]], [mz]):
            p = self.generate(files=files, acquisition="DIA")
            self.assertNotEqual(p.returncode, 0, files)
            for words in ("REFUSED", "nothing was written", "--engine diann", "mkdir -p"):
                self.assertIn(words, p.stderr)
            self.assertFalse(os.path.exists(self.job))
            self.assertFalse(os.path.exists(self.out))
        self.assertEqual(listing(self.raw), self.before + ["S4.mzML"])

    def test_fragpipe_dia_on_d_keeps_its_stage(self):
        """(pin) diatracer_stage.py links the .d"""
        p = self.generate(files=[os.path.join(self.raw, "S3.d")], acquisition="DIA")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertTrue(os.path.islink(os.path.join(self.out, "diatracer_stage", "S3.d")))


class TheDatabaseSearched(StageHarness):
    """Follow-up 3 (H7 bug 1): the shipped workflows have no database.db-path, and headless
    FragPipe then fills it from the account's per-user cache -- the last search's database. The
    route writes its own copy of the workflow with database.db-path = --fasta, refuses a workflow
    that names another database, and records what it searched (path + sha256)."""

    def test_the_job_searches_the_fasta_it_was_given(self):
        import hashlib
        p = self.generate()
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        fasta = os.path.join(self.d, "db.fasta")
        run_wf = os.path.join(self.out, run_search.FRAGPIPE_RUN_WORKFLOW)
        self.assertIn(f"--workflow {run_wf}", open(self.job).read())
        with open(run_wf, newline="") as fh:
            lines = fh.read().splitlines()
        self.assertEqual([ln for ln in lines if ln.startswith("database.db-path=")],
                         ["database.db-path=" + fasta.replace(":", "\\:")])
        prov = json.load(open(os.path.join(self.out, "search_provenance.json")))
        db = prov["fragpipe_database"]
        self.assertEqual(db["path"], fasta)
        self.assertEqual(db["sha256"], hashlib.sha256(open(fasta, "rb").read()).hexdigest())
        self.assertIsNone(db["workflow_db_path"])          # the shipped workflow names none

    def test_a_workflow_naming_another_database_is_refused(self):
        wf = open(FP_WORKFLOW, newline="").read() + "database.db-path=/elsewhere/other.fas\r\n"
        p = self.generate(workflow=wf)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("REFUSED", p.stderr)
        self.assertIn("/elsewhere/other.fas", p.stderr)
        self.assertFalse(os.path.exists(self.job))

    def test_a_workflow_naming_the_same_database_is_searched(self):
        fasta = os.path.join(self.d, "db.fasta")
        wf = open(FP_WORKFLOW, newline="").read() + f"database.db-path={fasta}\r\n"
        p = self.generate(workflow=wf)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        prov = json.load(open(os.path.join(self.out, "search_provenance.json")))
        self.assertEqual(prov["fragpipe_database"]["workflow_db_path"], fasta)

    def test_a_target_only_fasta_is_said(self):
        p = self.generate()                              # db.fasta has no rev_ entries
        self.assertIn("no decoys", p.stdout + p.stderr)


class InlineRoute(StageHarness):
    """F1 (review of 9abc8390): inline, on Linux, a .raw search puts mono on PATH before FragPipe.
    The mono line ended in `echo "...$(command -v mono)"` and was joined to the FragPipe command
    with a space, so echo printed the command and FragPipe never ran."""

    def test_the_inline_search_runs_fragpipe_with_mono(self):
        import run_search
        out = os.path.join(self.d, "inline_out")
        params = os.path.join(self.d, "search.workflow")
        shutil.copy(FP_WORKFLOW, params)
        # a real (target + decoy) FASTA: the route reads --fasta for the database record
        # (follow-up 3, item 1) and refuses one it cannot read
        fasta = _write(os.path.join(self.d, "db.fasta"),
                       ">sp|P1|X\nPEPTIDER\n>rev_sp|P1|X\nREDITPEP\n")
        res = run_search.run_fragpipe(self.fp, {"acquisition": "DDA"}, params, self.files[:2],
                                      fasta, out, 2, None,
                                      mono={"bin": self.mono, "source": "test"})
        self.assertTrue(os.path.exists(os.path.join(out, ".fake_fragpipe_ran")))
        self.assertIsNone(res["report"])                      # an IDs-only workflow
        self.assertTrue(os.path.isfile(res["identifications"]))
        self.assertEqual(listing(self.raw), self.before)


class ExperimentNamesBackToRuns(unittest.TestCase):
    """FragPipe turns '-' into '_' in an experiment name (its own fragpipe-files.fp-manifest:
    S-1 -> S_1; H5 gap 8), so combined_protein.tsv's columns are not the run names
    conditions.csv uses. The adapter maps them back through the two manifests."""

    def write(self, d, ours, theirs):
        cols = ["Protein", "Protein ID", "Entry Name", "Gene"]
        exps = [e for _, e in theirs]
        with open(os.path.join(d, "combined_protein.tsv"), "w") as fh:
            fh.write("\t".join(cols + [f"{e} Intensity" for e in exps]
                               + [f"{e} MaxLFQ Intensity" for e in exps]) + "\n")
            fh.write("\t".join(["sp|P02768|ALBU_HUMAN", "P02768", "ALBU_HUMAN", "ALB"]
                               + ["100"] * (2 * len(exps))) + "\n")
        for name, rows in (("fragpipe.fp-manifest", ours), ("fragpipe-files.fp-manifest", theirs)):
            if rows is not None:
                with open(os.path.join(d, name), "w") as fh:
                    fh.write("".join(f"{p}\t{e}\t\tDDA\n" for p, e in rows))

    def runs(self, d):
        import pyarrow.parquet as pq
        return sorted(set(pq.read_table(run_search.adapt_fragpipe_dda(d)).column("Run")
                          .to_pylist()))

    def test_fragpipe_s_names_go_back_to_the_run_names(self):
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            self.skipTest("pyarrow")
        with tempfile.TemporaryDirectory() as d:
            st = os.path.join(d, "raw_stage")
            self.write(d, [(f"{st}/S1-a.raw", "S1-a"), (f"{st}/S2-b.raw", "S2-b"),
                           (f"{st}/S3.d", "S3")],
                       [(f"{st}/S1-a.raw", "S1_a"), (f"{st}/S2-b.raw", "S2_b"), (f"{st}/S3.d", "S3")])
            self.assertEqual(self.runs(d), ["S1-a", "S2-b", "S3"])

    def test_without_fragpipe_s_manifest_the_columns_are_kept(self):
        """(pin) an output from before, or from outside the skill"""
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            self.skipTest("pyarrow")
        with tempfile.TemporaryDirectory() as d:
            self.write(d, None, [("x/S1_a.raw", "S1_a"), ("x/S2.raw", "S2")])
            os.remove(os.path.join(d, "fragpipe-files.fp-manifest"))
            self.assertEqual(self.runs(d), ["S1_a", "S2"])


class StageIsScratch(unittest.TestCase):
    def test_the_one_rule(self):
        self.assertTrue(scratch_files.is_scratch_dir("raw_stage"))
        for fn in ("S1_uncalibrated.mzML", "S1_uncalibrated.mzBIN", "S1.mzBIN"):
            self.assertTrue(scratch_files.is_scratch_file(fn), fn)
        self.assertFalse(scratch_files.is_scratch_file("S1.mzML"))
        self.assertIn("raw_stage", scratch_files.LABEL)

    def test_an_interim_delivery_skips_it_with_its_reason(self):
        import core_submission
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "output")
            _write(os.path.join(out, "search", "raw_stage", "S1.raw"), "raw\n")
            _write(os.path.join(out, "search", "psm.tsv"), "x\n")
            items = core_submission._plan_tree(out, "search")
            by = {i["name"]: i for i in items}
            st = by[os.path.join("search", "raw_stage") + "/"]
            self.assertEqual(st["status"], "SKIPPED")
            self.assertIn("staged copies of the raw files", st["reason"])
            self.assertEqual(by[os.path.join("search", "psm.tsv")]["status"], "OK")

    def test_finalize_zip_leaves_it_out(self):
        from test_finalize_zip_logs import finish_fixture
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(job_env(tmp), CLAUDE_CONFIG_DIR=os.path.join(tmp, "cc"))
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "init",
                                "--name", "Stage zip", "--base", os.path.join(tmp, "base")],
                               capture_output=True, text=True, env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            sdir = json.loads(r.stdout)["paths"]["session_dir"]
            finish_fixture(sdir)
            for rel in ("output/search/raw_stage/S1.raw", "output/search/raw_stage/S1.mzBIN",
                        "output/search/raw_stage/S1_uncalibrated.mzML", "output/search/psm.tsv"):
                _write(os.path.join(sdir, rel), "x\n")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "finalize",
                                "--dir", sdir, "--zip", "--no-deposit", "--no-notify"],
                               capture_output=True, text=True, env=env, timeout=600)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            res = json.loads(r.stdout)
            names = zipfile.ZipFile(res["zip"]).namelist()
            self.assertFalse([n for n in names if "raw_stage" in n], names)
            self.assertTrue([n for n in names if n.endswith("output/search/psm.tsv")], names)
            self.assertGreaterEqual(res["zip_excluded"][scratch_files.LABEL], 3)


if __name__ == "__main__":
    unittest.main()
