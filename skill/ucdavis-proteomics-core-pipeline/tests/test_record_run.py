#!/usr/bin/env python3
"""
record_run.py is the Core's run registry: one session folder per search on HIVE (the layout of
Brett's DataAnalysis sessions), an append-only master log and an activity log. Nobody is waiting
on it, so every failure here is SILENT -- a record that is missing, duplicated, misattributed,
wrong, or that leaked a credential into a folder the whole Core group reads would go unnoticed
for months. The assertions are about those:

  * the folder is sessions/<YYYY-MM-DD_Short-Description>, found again by the search's identity
    (not its name), so a re-record UPDATES it; a different search wanting the name gets _2;
  * the job-end hook's contract: ONE JSON object on stdout, exit 0 whatever happened;
  * data_analysis.md and activity_log.csv are append-only, rows are single write() calls with the
    exact columns, and the same event is never logged twice in the master log;
  * the CoreOmics submission is recorded when known and says "not recorded" (plus a Data Quality
    Note) when not -- never guessed from a folder name;
  * SEARCH_LOG.md always has a Data Quality Notes section;
  * big files, the zip over its cap, per-run .quant and anything secret are never copied;
  * a collaborator (not the Core) is skipped; RECORD_RUN=off touches nothing and sends nothing;
  * the SSH route (a laptop in hive_remote mode) works through hive_exec.sh.

Nothing contacts HIVE or SLURM, and nothing can reach the real /quobyte registry: SKILL_RUNS_DIR
always points into a temp dir, HIVE_EXEC is a stub over another temp dir, sacct is disabled or a
stub, and SLURM_* / HIVE_* are removed from the environment.
"""
import csv
import getpass
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time
import unittest
from unittest import mock
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
SCRIPT = os.path.join(SCRIPTS, "record_run.py")
sys.path.insert(0, SCRIPTS)
import record_run  # noqa: E402

FAKE_ROOT = "/nonexistent_hive_root_for_record_run_tests"    # never exists -> never "direct"
SECRET = "ghp_" + "A1b2C3d4E5f6G7h8I9j0K1l2M3n4"
USER = record_run.clean(getpass.getuser())
DAY = time.mktime((2026, 9, 23, 17, 31, 0, 0, 0, -1))
TS_RE = re.compile(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}[+-]\d{2}:\d{2}$")
SESSION = "2026-09-23_chkLUppm_HeLa50"

RUNS = ["FL030226_HeL50_35m_3", "FL050326_HeL50_35m", "FL190326_HeL50_35m"]
STATS = ("File.Name\tPrecursors.Identified\tProteins.Identified\tMS1.Signal\tMS2.Signal\n"
         + "".join(f"/nfs/lssc0/flinders/proteomics/Data/raw_data/Lumos1/x/{r}.raw\t{p}\t{q}\t1e11\t1e10\n"
                   for r, p, q in zip(RUNS, (28474, 28563, 29552), (3593, 3730, 3723))))
CFG = """--qvalue 0.01
--matrices
--xic 10
--fasta-search
--predictor
--reanalyse
--cut K*,R*
--missed-cleavages 1
--min-pep-len 7
--max-pep-len 30
--min-pr-mz 357
--max-pr-mz 1105
--min-pr-charge 2
--max-pr-charge 4
--unimod4
--cont-quant-exclude Cont_
"""


def write(path, text="", mode="w"):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, mode) as fh:
        fh.write(text)
    return path


def make_search(root, name="search_out", report=True, runs=RUNS, secrets=False, log_extra="",
                clean_db=False, resolution=False, run_dir="run1"):
    """A search folder shaped like the real HIVE e2e search (DIA-NN 2.7.0, 3 Lumos .raw, single
    shot, library predicted in a first job), with its wf/ and FASTA sidecar beside it."""
    run = os.path.join(root, run_dir)
    out = os.path.join(run, name)
    files = [f"/nfs/lssc0/flinders/proteomics/Data/raw_data/Lumos1/x/{r}.raw" for r in runs]
    write(os.path.join(run, "wf", "params.cfg"), CFG)
    rationale = {
        "engine": "diann", "acquisition": "DIA", "instrument": "Orbitrap Fusion Lumos",
        "instrument_class": "orbitrap_generic", "class_label": "Orbitrap (resolution unknown)",
        "mass_accuracy_source": "no documented DIA-NN value (resolution unknown)",
        "mass_accuracy_plan": "auto",
        "rationale": {"--qvalue": {"value": 0.01, "source": "standard 1% precursor FDR"},
                      "--min-pr-mz": {"value": 357, "source": "measured from the acquired "
                                                               "isolation windows"},
                      "--cut": {"value": "K*,R*", "source": "trypsin/P"}}}
    if resolution:
        rationale.update(instrument_class="orbitrap_120k", class_label="Orbitrap 120k/30k",
                         resolution={"ms1": 120000, "ms2": 30000, "source": "user",
                                     "source_label": "the user said so"})
    write(os.path.join(run, "wf", "params.cfg.rationale.json"), json.dumps(rationale))
    write(os.path.join(run, "wf", "workflow.manifest.json"), json.dumps({
        "acquisition": "DIA", "instruments": ["Orbitrap Fusion Lumos"],
        "engine": {"name": "diann", "version": "2.7.0"}}))
    fasta = os.path.join(run, "search.fasta")
    write(fasta, ">sp|P1|X\nPEPTIDE\n")
    meta = {"fasta": fasta, "md5": "e697e1f6", "n_entries": 21044, "n_sequences": 21044,
            "organism": "Homo sapiens", "taxid": 9606, "organism_source": "uniprot_api:UP000005640",
            "proteome": "UP000005640", "proteome_type": "Reference proteome", "n_proteome": 20663,
            "source": "hive:/quobyte/proteomics-grp/MRS/UP000005640_9606.fasta",
            "content_requested": "one_per_gene", "content_used": "as_staged",
            "staged_file": {"mtime_utc": "2025-04-25T23:25:36Z"},
            "n_contaminants_appended": 381, "n_contaminants_already_present": 0,
            "contaminant_set": "universal", "diann_cont_quant_exclude": "Cont_",
            "digestion_enzymes_used": ["trypsin", "lysc"]}
    if clean_db:
        meta["contaminant_target_rule"] = "drop contaminants identical to a target protein"
    write(fasta + ".meta.json", json.dumps(meta))
    write(os.path.join(out, "search_provenance.json"), json.dumps({
        "engine": "diann", "version": "2.7.0",
        "engine_version": {"value": "2.7.0", "source": "tools.json versions.diann"},
        "resolved_command": "/quobyte/proteomics-grp/dia-nn/build_270/diann-2.7.0/diann-linux",
        "params_file": os.path.join(run, "wf", "params.cfg"),
        "resolved_params_file": os.path.join(run, "wf", "params.cfg"),
        "scan_window": {"source": "not in the cfg -- DIA-NN chooses the radius itself",
                        "value": None},
        "fasta": fasta, "threads": 16, "n_files": len(files), "files": files,
        "search_mode": "single_shot",
        "parallel_routing_reason": "3 file(s), at or below the threshold of 5",
        "result": {"engine": "diann", "report": os.path.join(out, "report.parquet"),
                   "mode": "two_job_libfree"}}))
    os.utime(os.path.join(out, "search_provenance.json"), (DAY, DAY))
    write(os.path.join(out, "report.log.txt"),
          "\nDIA-NN 2.7.0 Academia  (Data-Independent Acquisition by Neural Networks)\n"
          "[1:45] Optimised mass accuracy: 24 ppm\n"
          "[0:57] Recommended MS1 mass accuracy setting: 8 ppm\n"
          "[4:20] Recommended MS1 mass accuracy setting: 7 ppm\n" + log_extra)
    write(os.path.join(out, "diann_libpred_23978018.log"), "library predicted\n")
    write(os.path.join(out, "diann_search_23978112.log"), "searching\n" + log_extra)
    write(os.path.join(out, "job_input_files.txt"), "\n".join(files) + "\n")
    for r in runs:
        write(os.path.join(out, "quant", f"{r}.quant"), "Q" * 2048)
        write(os.path.join(out, "report_xic", f"{r}.xic.parquet"), "X" * 64)
    if report:
        write(os.path.join(out, "report.parquet"), "P" * 4096)
        write(os.path.join(out, "report.stats.tsv"), STATS)
        write(os.path.join(out, "fran_deposit.json"), json.dumps({
            "status": "staged", "entry": "/quobyte/proteomics-grp/fran/incoming/x__abcd1234",
            "organism": "Homo sapiens", "taxon": 9606, "staged_by": USER,
            "linked": ["report.parquet", "report.log.txt"]}))
    if secrets:
        write(os.path.join(out, "hive.env"), "HIVE_USER=someone\nHIVE_KEY=~/.ssh/id_ed25519\n")
        write(os.path.join(out, ".pgfarm_token"), "eyJhbGciOi.secret.token\n")
        write(os.path.join(out, "slack_webhook.cfg"), "https://hooks.slack.com/services/T/B/X\n")
        write(os.path.join(out, "notes.log.txt"), f"pushed with {SECRET}\n")
    return out


def make_session(root, out=None, quant_in_zip=True, zip_secret=False, big=0, extras=True,
                 name=SESSION, fasta_in_zip=False, manifest_in_zip=False):
    """A finalized session (session.py layout) and its zip, the way finalize --zip builds it:
    everything under the session folder, including a search run into output/search."""
    sess = os.path.join(root, name)
    write(os.path.join(sess, "input", "conditions.csv"), "Run,Group\nA1,A\nB1,B\n")
    write(os.path.join(sess, "input", "raw_files.txt"), "# raw\n/nfs/x/A1.raw\n")
    write(os.path.join(sess, "output", "tables", "de_provenance.json"), json.dumps({
        "method": "dpc", "contrasts": ["B-A"], "significant_per_contrast": {"B-A": 42},
        "q_cutoff": 0.01, "logfc": 1.0, "adjp": 0.05}))
    write(os.path.join(sess, "output", "tables", "methods.txt"), "DE methods\n")
    write(os.path.join(sess, "output", "tables", "DE_B-A.csv"), "protein,logFC\nP1,2\n")
    write(os.path.join(sess, "output", "methods.md"), "# Methods\n")
    write(os.path.join(sess, "output", "HeLa50_Report.docx"), "PK-docx-report")
    write(os.path.join(sess, "output", "METHODS.docx"), "PK-docx-methods")
    write(os.path.join(sess, "output", "DATA_SUBMISSION", "HOW_TO_SUBMIT.md"), "# How\n")
    write(os.path.join(sess, "output", "reproducibility", "REPRODUCE.md"), "# Reproduce\n")
    write(os.path.join(sess, "output", "reproducibility", "reproduce.sh"), "#!/bin/bash\n")
    write(os.path.join(sess, "logs", "commands.log"), "python3 run_search.py ...\n")
    write(os.path.join(sess, "MANIFEST.txt"),
          "Session export manifest\n[OK] Publication methods\n[SKIPPED] Word copy -- no pandoc\n")
    readme = "# session\n"
    if extras:
        readme += ("\n## Expert Review Notes\n\n- **warning** -- n = 3 per group; "
                   "underpowered for 1.5-fold changes.\n\n## Next\n")
        write(os.path.join(sess, "output", "AUDIT.json"), json.dumps({
            "overall": "WARN", "findings": [
                {"check": "replication", "status": "PASS", "message": "ok"},
                {"check": "de_signal", "status": "WARN",
                 "message": "B-A: 0/3000 proteins significant"}]}))
        write(os.path.join(sess, "output", "SAMPLE_QUALITY.json"), json.dumps({
            "flags": ["**HEMOLYSIS is CONFOUNDED WITH GROUP** (B high)."]}))
        write(os.path.join(sess, "input", "detect.json"), json.dumps({
            "overall": "DIA", "instrument": "Orbitrap Fusion Lumos",
            "orbitrap_resolution_unknown": {"files": ["/x/A1.raw", "/x/B1.raw"]},
            "files": [{"file": "/x/A1.raw", "warnings": ["instrument clock skew"]},
                      {"file": "/x/B1.raw", "warnings": ["instrument clock skew"]}]}))
    write(os.path.join(sess, "README.md"), readme)
    if out:
        # a hive_remote session holds a COPY of the provenance, pulled from where it ran
        write(os.path.join(sess, "output", "search", "search_provenance.json"), json.dumps(
            {"engine": "diann", "version": "2.7.0",
             "result": {"report": os.path.join(out, "report.parquet")}}))
    zpath = sess + ".zip"
    with zipfile.ZipFile(zpath, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("s/README.md", "# session\n")
        z.writestr("s/output/tables/DE_B-A.csv", "protein,logFC\nP1,2.0\n" * 50)
        z.writestr("s/output/search/report.parquet", "P" * 4096)
        z.writestr(zipfile.ZipInfo("s/output/search/report.stats.tsv"), STATS)  # STORED
        if big:
            z.writestr("s/output/figures/big.bin", os.urandom(big))
        if quant_in_zip:
            for r in RUNS:
                z.writestr(f"s/output/search/quant/{r}.quant", os.urandom(4096))
            z.writestr("s/output/search/quant_step2_orig/x.quant", os.urandom(512))
            z.writestr("s/output/search/step1.predicted.speclib", os.urandom(6000))
            z.writestr("s/output/search/report-lib.parquet.skyline.speclib", "kept\n")
        if zip_secret:
            z.writestr("s/input/hive.env", "HIVE_USER=x\n")
        if fasta_in_zip:
            z.writestr("s/input/search.fasta", ">sp|P1|X\nPEPTIDE\n" * 400)
            z.writestr("s/input/extra.fa.gz", os.urandom(800))
        if manifest_in_zip:
            z.writestr("s/MANIFEST.txt", "zipped manifest\n")
    return sess, zpath


class Base(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = os.path.realpath(self._tmp.name)
        self.runs = os.path.join(self.d, "grp", "skill_runs")
        self.issues = os.path.join(self.d, "grp", "skill_issues")
        os.makedirs(self.runs)
        os.makedirs(self.issues)
        self.remote = os.path.join(self.d, "remote")                  # FAKE_ROOT on "HIVE"
        self.key = write(os.path.join(self.d, "id_test"))
        self.relay_log = os.path.join(self.d, "relay_calls.log")
        self.fake_exec = os.path.join(self.d, "fake_hive_exec.sh")
        write(self.fake_exec, f"""#!/usr/bin/env bash
echo "$*" >> '{self.relay_log}'
FAKE='{FAKE_ROOT}'; REMOTE='{self.remote}'
if [ "$1" = "--put" ]; then
  src="$2"; dst="$3"; case "$dst" in /*) ;; *) dst="$HOME/$dst" ;; esac
  dst="${{dst//$FAKE/$REMOTE}}"
  mkdir -p "$dst" && cp -R "$src" "$dst/" || exit 1
  b="$dst/$(basename "$src")/bundle.json"
  # the paths inside the upload name HIVE locations; on the stand-in they live under REMOTE
  [ -f "$b" ] && python3 -c 'import sys; p,a,b=sys.argv[1:]; s=open(p).read().replace(a,b); open(p,"w").write(s)' "$b" "$FAKE" "$REMOTE"
  exit 0
fi
cmd="${{1//$FAKE/$REMOTE}}"
exec bash -c "$cmd"
""")
        os.chmod(self.fake_exec, 0o755)

    def tearDown(self):
        for root, dirs, _ in os.walk(self.d):
            for x in dirs:
                try:
                    os.chmod(os.path.join(root, x), 0o755)
                except OSError:
                    pass
        self._tmp.cleanup()

    def env(self, route="direct", **extra):
        e = {k: v for k, v in os.environ.items()
             if not k.startswith(("HIVE_", "SKILL_", "RECORD_RUN", "SLURM_", "FRAN_"))}
        e.update(HOME=self.d, HIVE_ENV_FILE=os.path.join(self.d, "no-hive.env"),
                 RECORD_RUN_SACCT="none", SKILL_ISSUES_DIR=self.issues, TMPDIR=self.d)
        if route == "direct":
            e["SKILL_RUNS_DIR"] = self.runs
        else:
            e["SKILL_RUNS_DIR"] = FAKE_ROOT + "/skill_runs"
            e["SKILL_ISSUES_DIR"] = FAKE_ROOT + "/skill_issues"
        if route == "ssh":
            e.update(HIVE_USER="gabrig", HIVE_KEY=self.key, HIVE_EXEC=self.fake_exec)
        e.update({k: str(v) for k, v in extra.items()})
        return e

    def run_raw(self, *args, route="direct", **extra):
        return subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True,
                              env=self.env(route, **extra), timeout=120)

    def run_it(self, *args, route="direct", **extra):
        """The hook contract: exit 0 and ONE JSON object on stdout."""
        p = self.run_raw(*args, route=route, **extra)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(len(p.stdout.strip().splitlines()), 1, p.stdout)
        res = json.loads(p.stdout)
        self.assertIsInstance(res.get("recorded"), bool, res)
        return res

    def folders(self, root=None):
        sd = os.path.join(root or self.runs, "sessions")
        return sorted(os.path.join(sd, f) for f in record_run.listdir(sd))

    def only_folder(self, root=None):
        f = self.folders(root)
        self.assertEqual(len(f), 1, f)
        return f[0]

    def read(self, folder, name="SEARCH_LOG.md"):
        with open(os.path.join(folder, name)) as fh:
            return json.load(fh) if name.endswith(".json") else fh.read()

    def all_files(self, folder):
        return [os.path.relpath(os.path.join(r, f), folder)
                for r, _, fs in os.walk(folder) for f in fs]

    def activity(self, root=None):
        with open(os.path.join(root or self.runs, "activity_log.csv")) as fh:
            text = fh.read()
        return text, list(csv.reader(io.StringIO(text)))

    def master(self, root=None):
        with open(os.path.join(root or self.runs, "data_analysis.md")) as fh:
            return fh.read()


class Layout(Base):
    def test_search_done_writes_a_session_folder(self):
        out = make_search(self.d)
        res = self.run_it("search-done", "--out", out, "--status", "completed",
                          "--exit-code", "0")
        self.assertTrue(res["recorded"], res)
        folder = self.only_folder()
        self.assertEqual(res["path"], folder)
        # no session: <date of submission>_<search dir name>; "search_out" climbs to its parent
        self.assertEqual(os.path.basename(folder), "2026-09-23_run1")
        rec = self.read(folder, "run_record.json")
        self.assertEqual(rec["schema_version"], record_run.SCHEMA_VERSION)
        self.assertEqual(rec["identity"]["out"], os.path.realpath(out))
        s = rec["search"]
        self.assertEqual((s["status"], s["engine"], s["engine_version"], s["exit_code"]),
                         ("completed", "diann", "2.7.0", 0))
        self.assertEqual(s["results"]["median_precursors"], 28563)
        log = self.read(folder)
        for want in ("**COMPLETED** (exit 0)", "**CoreOmics submission:** not recorded",
                     "## Data Quality Notes", "DIA-NN 2.7.0", "tools.json versions.diann",
                     "Orbitrap Fusion Lumos", "357–1105",
                     "measured from the acquired isolation windows",
                     "not fixed -- the engine optimised it per run", "optimised MS2 24 ppm",
                     "0.01 (precursor)", "Trypsin/P", "Homo sapiens (taxid 9606)",
                     "381 (universal)", "trypsin, lysc", "28,563", "3,723", "staged",
                     "left in place, never copied", "legacy database",
                     "Orbitrap resolution not recorded"):
            self.assertIn(want, log, want)
        files = set(self.all_files(folder))
        for want in ("input/params.cfg", "input/params.cfg.rationale.json",
                     "input/workflow.manifest.json", "input/search.fasta.meta.json",
                     "input/raw_files.txt", "output/search/search_provenance.json",
                     "output/search/report.stats.tsv", "output/search/report.log.txt",
                     "output/search/fran_deposit.json", "output/search/diann_search_23978112.log",
                     "SEARCH_LOG.md", "run_record.json"):
            self.assertIn(want, files)
        with open(os.path.join(folder, "input", "raw_files.txt")) as fh:
            self.assertIn("FL190326_HeL50_35m.raw", fh.read())
        self.assertTrue(os.path.islink(os.path.join(folder, "output", "search", "report.parquet")))
        self.assertEqual(os.readlink(os.path.join(folder, "output", "search", "search_out")), out)
        real = [f for f in files if not os.path.islink(os.path.join(folder, f))]
        self.assertFalse([f for f in real if f.endswith((".quant", ".xic.parquet", ".parquet"))])

    def test_session_name_is_the_folder_name(self):
        sess, _ = make_session(self.d, quant_in_zip=False)
        out = make_search(os.path.join(sess, "output"), name="search", run_dir=".")
        self.run_it("search-done", "--out", out, "--session", sess)
        self.assertEqual(os.path.basename(self.only_folder()), SESSION)

    def test_a_different_search_with_the_same_name_gets_a_suffix(self):
        a = make_search(os.path.join(self.d, "a"))
        b = make_search(os.path.join(self.d, "b"))      # also "run1", also 2026-09-23
        self.run_it("search-done", "--out", a)
        self.run_it("search-done", "--out", b)
        self.run_it("search-done", "--out", a)           # still the first folder
        names = [os.path.basename(f) for f in self.folders()]
        self.assertEqual(names, ["2026-09-23_run1", "2026-09-23_run1_2"])
        self.assertEqual(self.read(self.folders()[0], "run_record.json")["identity"]["out"],
                         os.path.realpath(a))
        self.assertEqual(len(self.read(self.folders()[0], "run_record.json")["history"]), 2)

    def test_re_record_updates_in_place_and_logs_append(self):
        out = make_search(self.d)
        self.run_it("search-done", "--out", out, "--status", "completed", "--exit-code", "0")
        master1 = self.master()
        self.run_it("search-done", "--out", out, "--status", "completed", "--exit-code", "0")
        sess, _ = make_session(self.d, quant_in_zip=False)
        self.run_it("analysis-done", "--session", sess, "--out", out)
        folder = self.only_folder()
        rec = self.read(folder, "run_record.json")
        self.assertEqual([h["event"] for h in rec["history"]],
                         ["search-done", "search-done", "analysis-done"])
        master = self.master()
        self.assertTrue(master.startswith(master1), "the master log was rewritten, not appended")
        self.assertEqual(master.count("\n## "), 1, master)                # one entry per run
        self.assertEqual(master.count("search completed -->"), 1)         # never twice
        self.assertEqual(master.count(": analysis complete"), 1, master)
        text, rows = self.activity()
        actions = [r[2] for r in rows[1:]]
        self.assertEqual(actions.count("search_completed"), 2)            # every call appends
        self.assertIn("re-recorded", text)
        self.assertEqual(actions.count("analysis_completed"), 1)
        self.assertEqual(actions.count("fran_staged"), 1)                 # state events once
        # finalize re-reads the search without --status: the job's own exit code must survive
        self.assertEqual((rec["search"]["status"], rec["search"]["exit_code"]), ("completed", 0))
        log = self.read(folder)
        self.assertIn("**COMPLETED** (exit 0)", log)
        self.assertIn("B-A = 42", log)
        self.assertIn("DATA_SUBMISSION", log)

    def test_intermediate_failure_then_success_is_one_record(self):
        out = make_search(self.d, report=False,
                          log_extra="ERROR: DIA-NN tried but failed to load the following files\n")
        write(os.path.join(out, "s2_firstpass_23977910_4.log"), "ERROR: task 4 died\n")
        res = self.run_it("search-done", "--out", out, "--status", "failed", "--exit-code", "1",
                          SLURM_JOB_ID="23977911", SLURM_JOB_NAME="s3_assembly")
        self.assertTrue(res["recorded"], res)
        folder = self.only_folder()
        log = self.read(folder)
        self.assertIn("**FAILED** (exit 1) at step `s3_assembly`", log)
        self.assertIn("job 23977911 `s3_assembly`", log)
        self.assertIn("failed to load the following files", log)
        self.assertIn("output/search/s2_firstpass_23977910_4.log", self.all_files(folder))
        # the fixed search completes: same folder, a dated update line, a new activity row
        write(os.path.join(out, "report.parquet"), "P" * 4096)
        write(os.path.join(out, "report.stats.tsv"), STATS)
        self.run_it("search-done", "--out", out, "--status", "completed", "--exit-code", "0")
        self.assertEqual(self.only_folder(), folder)
        master = self.master()
        self.assertIn("FAILED (exit 1, step s3_assembly)", master)
        self.assertEqual(master.count("\n## "), 1)
        update = master.split(": search completed", 1)[1]
        self.assertNotIn("Run by", update)              # an update is a dated line, not an entry
        _, rows = self.activity()
        self.assertEqual([r[2] for r in rows[1:] if r[2].startswith("search_")],
                         ["search_failed", "search_completed"])

    def test_dry_run_writes_nothing(self):
        out = make_search(self.d)
        p = self.run_raw("search-done", "--out", out, "--dry-run")
        self.assertEqual(p.returncode, 0)
        res = json.loads(p.stdout)
        self.assertEqual((res["recorded"], res["reason"]), (False, "dry_run"))
        self.assertIn("# Search log -- 2026-09-23_run1", p.stderr)
        self.assertIn("## Data Quality Notes", p.stderr)
        self.assertEqual(os.listdir(self.runs), [])


class Logs(Base):
    def test_activity_log_columns_and_timestamps(self):
        out = make_search(self.d)
        self.run_it("search-done", "--out", out, "--status", "completed", "--exit-code", "0")
        text, rows = self.activity()
        self.assertTrue(text.startswith("timestamp,session,action,tool,target,status,notes\n"))
        self.assertEqual(text.count("timestamp,session"), 1)
        for r in rows[1:]:
            self.assertEqual(len(r), 7, r)
            self.assertRegex(r[0], TS_RE)
            self.assertEqual(r[1], "2026-09-23_run1")
        done = [r for r in rows if r[2] == "search_completed"][0]
        self.assertEqual((done[3], done[4], done[5]), ("DIA-NN 2.7.0", out, "completed"))
        self.assertIn("median 28,563 precursors", done[6])

    def test_each_row_is_one_write_call(self):
        path = os.path.join(self.d, "a.csv")
        record_run.append_locked(path, ["x\n"], record_run.Deadline(10), header="h\n")
        calls = []
        real = record_run.os.write
        record_run.os.write = lambda fd, b: calls.append(b) or real(fd, b)
        rows = [record_run.csv_row(["t", "s", "a", "b", "c", "d", f"note {i}, with a comma"])
                for i in range(3)]
        try:
            record_run.append_locked(path, rows, record_run.Deadline(10), header="h\n")
        finally:
            record_run.os.write = real
        self.assertEqual(len(calls), 3)
        self.assertTrue(all(c.count(b"\n") == 1 for c in calls))

    def test_two_writers_200_rows_each_nothing_lost_or_interleaved(self):
        """The lead's check: two processes, 200 appends each, each append taking the lock."""
        path = os.path.join(self.d, "act.csv")
        prog = ("import sys; sys.path.insert(0, %r); import record_run as r\n"
                "for i in range(200):\n"
                "    r.append_locked(%r, [r.csv_row([r.ts_min(), 'S', 'act', 'tool', 't', 'ok', "
                "'writer ' + sys.argv[1] + ' row ' + str(i) + ' ' + 'x' * 180])], "
                "r.Deadline(60), header='timestamp,session,action,tool,target,status,notes\\n')\n"
                "print(len(r.LOCK_TIMEOUTS))\n" % (SCRIPTS, path))
        procs = [subprocess.Popen([sys.executable, "-c", prog, str(n)], stdout=subprocess.PIPE,
                                  text=True) for n in range(2)]
        for p in procs:
            out, _ = p.communicate(timeout=120)
            self.assertEqual((p.returncode, out.strip()), (0, "0"))      # no lock timeouts
        with open(path) as fh:
            lines = fh.read().splitlines()
        self.assertEqual(lines[0], ",".join(record_run.ACTIVITY_COLUMNS))
        rows = list(csv.reader(lines[1:]))
        self.assertTrue(all(len(r) == 7 for r in rows))                 # nothing torn
        seen = sorted(r[6].split(" x")[0] for r in rows)
        want = sorted(f"writer {w} row {i}" for w in range(2) for i in range(200))
        self.assertEqual(seen, want)                                      # nothing lost, no dupes
        self.assertFalse(os.path.exists(path + ".lock.d"))

    def test_a_stale_lock_is_broken_and_a_live_one_times_out(self):
        path = os.path.join(self.d, "m.csv")
        os.mkdir(path + ".lock.d")
        old = time.time() - 120
        os.utime(path + ".lock.d", (old, old))          # a dead holder's lock
        del record_run.LOCK_TIMEOUTS[:]
        t0 = time.monotonic()
        record_run.append_locked(path, ["a\n"], record_run.Deadline(30))
        self.assertLess(time.monotonic() - t0, 5)
        self.assertEqual(record_run.LOCK_TIMEOUTS, [])
        os.mkdir(path + ".lock.d")                       # a live holder: wait, then write anyway
        os.environ["RECORD_RUN_LOCK_WAIT"] = "0.5"
        try:
            record_run.append_locked(path, ["b\n"], record_run.Deadline(30))
        finally:
            del os.environ["RECORD_RUN_LOCK_WAIT"]
        self.assertEqual(record_run.LOCK_TIMEOUTS, ["m.csv"])
        with open(path) as fh:
            self.assertEqual(fh.read(), "a\nb\n")
        self.assertTrue(os.path.isdir(path + ".lock.d"))  # someone else's lock is left alone

    def test_lock_race_abc_a_live_lock_is_never_broken(self):
        """A died holding the lock. B judges A's lock stale -- and, before B renames it, C breaks
        A's lock and takes a FRESH one. B's rename then moves C's LIVE lock: B must notice (the
        moved directory's owner is C, not the A it judged), put it back, and wait like anyone
        else. (FRAN ad81863; plain mkdir-and-break let two writers hold the lock at once.)"""
        import contextlib
        path = os.path.join(self.d, "abc.csv")
        lkd = path + ".lock.d"
        os.mkdir(lkd)
        write(os.path.join(lkd, "owner"), "A:dead:1\n")
        os.utime(lkd, (time.time() - 3600,) * 2)
        b = record_run.DirLock(path, record_run.Deadline(30), wait=0.5)

        def c_interleaves(judged):
            self.assertEqual(judged, "A:dead:1")
            b.break_hook = None                          # once
            os.rename(lkd, lkd + ".c_grave")             # C breaks A's stale lock...
            shutil.rmtree(lkd + ".c_grave")
            os.mkdir(lkd)                                # ...and takes a fresh one
            write(os.path.join(lkd, "owner"), "C:live:2\n")
        b.break_hook = c_interleaves
        del record_run.LOCK_TIMEOUTS[:]
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            with b:
                self.assertFalse(b.held)                 # B waited, then went on unlocked
        self.assertEqual(record_run.DirLock.owner_of(lkd), "C:live:2")   # C's lock survives
        self.assertIn("restored it to its live holder", err.getvalue())
        self.assertEqual(record_run.LOCK_TIMEOUTS, ["abc.csv"])
        self.assertEqual([f for f in os.listdir(self.d) if f.startswith("abc.csv.lock.d.")], [])

    def test_a_break_that_keeps_failing_times_out_instead_of_spinning(self):
        """FRAN auto_ingest_state._lock spun forever at 100% CPU when the stale lock's rename kept
        failing (EACCES): every failed break `continue`d past the deadline check. Here the rename
        always fails; the lock must give up at its wait and write unlocked."""
        import contextlib
        path = os.path.join(self.d, "stuck.csv")
        lkd = path + ".lock.d"
        os.mkdir(lkd)
        write(os.path.join(lkd, "owner"), "A:dead:1\n")
        os.utime(lkd, (time.time() - 3600,) * 2)
        real_rename = os.rename

        def failing_rename(src, dst):
            if src == lkd:
                raise PermissionError(13, "Permission denied", src)
            return real_rename(src, dst)
        del record_run.LOCK_TIMEOUTS[:]
        t0 = time.monotonic()
        with mock.patch.object(record_run.os, "rename", failing_rename), \
                contextlib.redirect_stderr(io.StringIO()):
            with record_run.DirLock(path, record_run.Deadline(30), wait=0.5) as lk:
                self.assertFalse(lk.held)
        self.assertLess(time.monotonic() - t0, 5)
        self.assertEqual(record_run.LOCK_TIMEOUTS, ["stuck.csv"])
        self.assertTrue(os.path.isdir(lkd))              # never removed by a failed break

    def test_an_ownerless_lock_is_broken_only_if_a_second_look_agrees(self):
        """No owner file can mean a dead holder OR a new holder between mkdir and writing its
        token. The breaker looks again; if an owner appeared meanwhile, the lock is left alone."""
        import contextlib
        path = os.path.join(self.d, "fresh.csv")
        lkd = path + ".lock.d"
        os.mkdir(lkd)
        os.utime(lkd, (time.time() - 3600,) * 2)       # stale, owner-less

        def new_holder_writes_token(_seconds):
            write(os.path.join(lkd, "owner"), "C:live:2\n")
            os.utime(lkd, None)
        with mock.patch.object(record_run, "OWNERLESS_SECOND_LOOK_S", 0.0), \
                mock.patch.object(record_run.time, "sleep", new_holder_writes_token), \
                contextlib.redirect_stderr(io.StringIO()):
            b = record_run.DirLock(path, record_run.Deadline(30), wait=0.0)
            b._try_break("B:x", 3600)
        self.assertEqual(record_run.DirLock.owner_of(lkd), "C:live:2")
        # and a truly dead owner-less lock IS broken
        shutil.rmtree(lkd)
        os.mkdir(lkd)
        os.utime(lkd, (time.time() - 3600,) * 2)
        with mock.patch.object(record_run, "OWNERLESS_SECOND_LOOK_S", 0.0), \
                contextlib.redirect_stderr(io.StringIO()):
            record_run.DirLock(path, record_run.Deadline(30), wait=0.0)._try_break("B:x", 3600)
        self.assertFalse(os.path.exists(lkd))

    def test_release_removes_only_our_own_lock(self):
        import contextlib
        path = os.path.join(self.d, "own.csv")
        lkd = path + ".lock.d"
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            with record_run.DirLock(path, record_run.Deadline(10)) as lk:
                self.assertTrue(lk.held)
                write(os.path.join(lkd, "owner"), "someone-else\n")  # broken + re-taken meanwhile
        self.assertTrue(os.path.isdir(lkd))
        self.assertIn("left in place", err.getvalue())
        shutil.rmtree(lkd)
        with record_run.DirLock(path, record_run.Deadline(10)) as lk:
            self.assertTrue(lk.held)
            self.assertEqual(record_run.DirLock.owner_of(lkd), lk.token)
        self.assertFalse(os.path.exists(lkd))

    def test_lock_timeout_is_reported_in_the_result(self):
        out = make_search(self.d)
        os.mkdir(os.path.join(self.runs, "activity_log.csv.lock.d"))
        res = self.run_it("search-done", "--out", out, RECORD_RUN_LOCK_WAIT="0.3")
        self.assertTrue(res["recorded"], res)
        self.assertEqual(res["lock_timeout"], ["activity_log.csv"])
        _, rows = self.activity()
        self.assertIn("search_completed", [r[2] for r in rows])        # still written

    def test_master_marker_is_checked_under_the_lock(self):
        path = os.path.join(self.d, "m.md")
        d = record_run.Deadline(10)
        self.assertEqual(record_run.append_locked(path, ["\n## e\n<!-- k -->\n"], d,
                                                  header="H\n", marker="<!-- k -->"), 1)
        self.assertEqual(record_run.append_locked(path, ["\n## e\n<!-- k -->\n"], d,
                                                  header="H\n", marker="<!-- k -->"), 0)
        with open(path) as fh:
            self.assertEqual(fh.read(), "H\n\n## e\n<!-- k -->\n")

    def test_an_unreadable_record_is_moved_aside_never_overwritten(self):
        out = make_search(self.d)
        self.run_it("search-done", "--out", out)
        folder = self.only_folder()
        with open(os.path.join(folder, "run_record.json"), "w") as fh:
            fh.write('{"half": ')                       # torn: another writer, or a crash
        res = self.run_it("search-done", "--out", out)
        self.assertTrue(res["recorded"], res)
        self.assertEqual(self.only_folder(), folder)     # the index still finds it
        aside = [f for f in os.listdir(folder) if f.startswith("run_record.json.unreadable-")]
        self.assertEqual(len(aside), 1)
        with open(os.path.join(folder, aside[0])) as fh:
            self.assertEqual(fh.read(), '{"half": ')
        self.assertEqual(self.read(folder, "run_record.json")["search"]["status"], "completed")
        self.assertIn("could not be read and was kept as", self.read(folder))


class CoreOmics(Base):
    def test_prot_forms_are_normalised(self):
        p = record_run.parse_prot
        self.assertEqual(p(["807"]), {"prot": "PROT_0807", "id": None})
        self.assertEqual(p(["PROT_0807"])["prot"], "PROT_0807")
        self.assertEqual(p(["prot-0807"])["prot"], "PROT_0807")
        self.assertEqual(p(["99922F5337F8"]), {"prot": None, "id": "99922f5337f8"})
        self.assertEqual(p(["PROT_0807 / 99922f5337f8"]),
                         {"prot": "PROT_0807", "id": "99922f5337f8"})
        self.assertIsNone(p(["hela"]))

    def test_prot_given_is_recorded_everywhere(self):
        out = make_search(self.d)
        self.run_it("search-done", "--out", out, "--prot", "807", "--prot", "99922f5337f8")
        folder = self.only_folder()
        self.assertIn("**CoreOmics submission:** PROT_0807 / `99922f5337f8`", self.read(folder))
        self.assertEqual(self.read(folder, "run_record.json")["prot"]["prot"], "PROT_0807")
        self.assertIn("PROT_0807", self.master())
        self.assertIn("PROT_0807", self.activity()[0])
        self.assertNotIn("CoreOmics submission: not recorded", self.read(folder))

    def test_prot_from_session_metadata_and_core_submission_receipt(self):
        sess, _ = make_session(self.d, quant_in_zip=False)
        write(os.path.join(sess, "input", "submission.json"),
              json.dumps({"coreomics": {"prot": "PROT_0744", "id": "3cba9a067ee4"}}))
        out = make_search(self.d)
        self.run_it("analysis-done", "--session", sess, "--out", out)
        self.assertIn("PROT_0744 / `3cba9a067ee4`", self.read(self.only_folder()))

        other = os.path.join(self.d, "work")
        sess2, _ = make_session(other, quant_in_zip=False, name="2026-09-24_Other_Study")
        write(os.path.join(other, ".core_submission.json"), json.dumps(
            {"schema": "core_submission/1", "internal_id": "PROT_0807", "id": "99922f5337f8",
             "session": sess2}))
        self.run_it("analysis-done", "--session", sess2)
        f2 = [f for f in self.folders() if f.endswith("Other_Study")][0]
        self.assertIn("PROT_0807 / `99922f5337f8`", self.read(f2))

    def test_a_parent_receipt_for_another_run_is_ignored(self):
        """P1: a .core_submission.json ABOVE the session belongs to whatever lives under that
        folder; it is used only when it names this session or search out dir."""
        work = os.path.join(self.d, "work")
        sess, _ = make_session(work, quant_in_zip=False, name="2026-09-24_Mine")
        write(os.path.join(work, ".core_submission.json"), json.dumps(
            {"schema": "core_submission/1", "internal_id": "PROT_0999", "id": "aaaaaaaaaaaa",
             "session": os.path.join(work, "2026-09-24_Someone_Else")}))
        self.run_it("analysis-done", "--session", sess)
        log = self.read(self.only_folder())
        self.assertIn("**CoreOmics submission:** not recorded", log)
        self.assertNotIn("PROT_0999", log)

    def test_prot_is_never_guessed_from_a_folder_name(self):
        sess, _ = make_session(self.d, quant_in_zip=False, name="2026-09-24_PROT_0999_Smith")
        self.run_it("analysis-done", "--session", sess)
        log = self.read(self.only_folder())
        self.assertIn("**CoreOmics submission:** not recorded", log)
        self.assertIn("CoreOmics submission: not recorded", log)          # the DQ note
        self.assertNotIn("PROT_0999 /", log)


class DataQuality(Base):
    def test_section_is_always_there_even_when_clean(self):
        out = make_search(self.d, clean_db=True, resolution=True)
        self.run_it("search-done", "--out", out, "--prot", "PROT_0807")
        log = self.read(self.only_folder())
        self.assertIn("## Data Quality Notes", log)
        self.assertIn("Nothing anomalous observed", log)

    def test_notes_collect_what_the_skill_found(self):
        out = make_search(self.d)
        sess, _ = make_session(self.d, quant_in_zip=True)
        write(os.path.join(self.issues, f"2026-09-22_{USER}_chkLUppm_HeLa50.md"),
              "# Skill issues\n<!-- end of header -->\n## ThermoRawFileParser not found  [bug]\n")
        write(os.path.join(self.issues, "2026-09-23_someoneelse_chkLUppm_HeLa50.md"),
              "## not ours\n")
        self.run_it("analysis-done", "--session", sess, "--out", out)
        folder = self.only_folder()
        log = self.read(folder)
        dq = log.split("## Data Quality Notes", 1)[1].split("\n## ", 1)[0]
        for want in ("B-A: 0/3000 proteins significant",           # AUDIT.json WARN
                     "CRITICAL** -- HEMOLYSIS is CONFOUNDED",       # SAMPLE_QUALITY flag
                     "Orbitrap resolution not recorded",            # detection + rationale
                     "instrument clock skew",                       # per-file detection warning
                     "legacy database",                             # FASTA sidecar
                     "per-run .quant files",                        # the zip finding
                     "ThermoRawFileParser not found",               # skill issue
                     "CoreOmics submission: not recorded",
                     "*Why it matters:*", "*Suggested fix:*"):
            self.assertIn(want, dq, want)
        self.assertNotIn("not ours", log)
        self.assertIn("## Expert Review Notes", log)
        self.assertIn("underpowered for 1.5-fold", log)
        _, rows = self.activity()
        self.assertIn("issue_recorded", [r[2] for r in rows])


class Readme(Base):
    def test_readme_is_written_when_missing_and_replaced_when_old(self):
        out = make_search(self.d)
        self.run_it("search-done", "--out", out)
        with open(os.path.join(self.runs, "README.md")) as fh:
            self.assertEqual(fh.read(), record_run.README_TEXT)
        old = "# Skill run registry\n\nOne folder per search: <YYYY>/<date>_<user>_<session>/\n"
        write(os.path.join(self.runs, "README.md"), old)          # the first, unversioned one
        self.run_it("search-done", "--out", out)
        with open(os.path.join(self.runs, "README.md")) as fh:
            self.assertEqual(fh.read(), record_run.README_TEXT)

    def test_a_current_readme_is_left_alone(self):
        path = os.path.join(self.runs, "README.md")
        mine = record_run.README_TEXT + "\nlocal note\n"
        write(path, mine)
        self.assertFalse(record_run.ensure_readme(self.runs, record_run.Deadline(10)))
        with open(path) as fh:
            self.assertEqual(fh.read(), mine)

    def test_the_reference_doc_quotes_the_readme_verbatim(self):
        with open(os.path.join(os.path.dirname(HERE), "references", "run-registry.md")) as fh:
            doc = fh.read()
        m = re.search(r"<!-- README:BEGIN -->\n```markdown\n(.*?)```\n<!-- README:END -->",
                      doc, re.S)
        self.assertIsNotNone(m)
        self.assertEqual(m.group(1), record_run.README_TEXT,
                         "references/run-registry.md and README_TEXT differ -- re-paste it")
        self.assertIn(f"README v{record_run.README_VERSION}", m.group(1))


class WhatIsCopied(Base):
    def test_the_session_is_mirrored_and_the_docx_listed_first(self):
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=False)
        self.run_it("analysis-done", "--session", sess, "--out", out)
        folder = self.only_folder()
        files = set(self.all_files(folder))
        for want in ("README.md", "MANIFEST.txt", os.path.basename(zpath),
                     "input/conditions.csv", "input/raw_files.txt",
                     "output/HeLa50_Report.docx", "output/METHODS.docx", "output/methods.md",
                     "output/tables/DE_B-A.csv", "output/tables/methods.txt",
                     "output/DATA_SUBMISSION/HOW_TO_SUBMIT.md", "scripts/commands.log",
                     "scripts/REPRODUCE.md", "scripts/reproduce.sh"):
            self.assertIn(want, files)
        with open(os.path.join(folder, "input", "raw_files.txt")) as fh:
            self.assertIn("/nfs/x/A1.raw", fh.read())                # the session's own list
        entry = self.master().split(": analysis complete", 1)[1]
        first = [ln for ln in entry.splitlines() if ln.startswith("- ")][0]
        self.assertIn("Report (Word)", first)
        self.assertIn(f"sessions/{SESSION}/output/HeLa50_Report.docx", first)
        self.assertIn(f"Reproducibility / zip:** `sessions/{SESSION}/{os.path.basename(zpath)}`",
                      entry)

    def test_big_files_are_skipped_with_a_note(self):
        out = make_search(self.d)
        write(os.path.join(out, "report.log.txt"), "DIA-NN 2.7.0\n" + "x" * 5000)
        self.run_it("search-done", "--out", out, "--file-cap-mb", "0.002")
        folder = self.only_folder()
        self.assertNotIn("output/search/report.log.txt", self.all_files(folder))
        log = self.read(folder)
        self.assertIn("## Not copied", log)
        self.assertIn("per-file cap", log)

    def test_zip_over_the_cap_is_not_copied(self):
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=False, big=200_000)
        self.run_it("analysis-done", "--session", sess, "--out", out, "--zip-cap-gb", "0.0001")
        folder = self.only_folder()
        self.assertFalse(any(f.endswith(".zip") for f in self.all_files(folder)))
        self.assertFalse(self.read(folder, "run_record.json")["zip_copy"]["copied"])
        self.assertIn("cap", self.read(folder))
        self.assertIn("not copied", self.master())

    def test_quant_and_predicted_library_are_never_copied_and_are_reported(self):
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=True)
        res = self.run_it("analysis-done", "--session", sess, "--out", out)
        self.assertEqual(res["findings"], ["session_zip_trimmed"])
        folder = self.only_folder()
        self.assertFalse([f for f in self.all_files(folder) if f.endswith(".quant")])
        with zipfile.ZipFile(os.path.join(folder, os.path.basename(zpath))) as z, \
                zipfile.ZipFile(zpath) as orig:
            self.assertIsNone(z.testzip())            # every kept member's CRC checks out
            self.assertEqual(sorted(z.namelist()), sorted(
                [n for n in orig.namelist() if not n.endswith((".quant", ".predicted.speclib"))]
                + ["s/MANIFEST.txt"]))                       # C4: added from the session
            self.assertIn("s/output/search/report-lib.parquet.skyline.speclib", z.namelist())
            for n in z.namelist():
                if n != "s/MANIFEST.txt":
                    self.assertEqual(z.read(n), orig.read(n))
        f = self.read(folder, "run_record.json")["findings"][0]
        self.assertEqual((f["n_quant"], f["n_predicted_speclib"]), (4, 1))
        self.assertIn("4 per-run .quant files and 1 predicted spectral library", f["detail"])
        self.assertIn("re-run session.py finalize --zip with skill", f["detail"])
        dq = self.read(folder).split("## Data Quality Notes", 1)[1].split("\n## ", 1)[0]
        self.assertIn("4 per-run .quant files and 1 predicted spectral library", dq)

    def test_the_zip_cap_counts_only_what_is_copied(self):
        """A zip over the cap only because of its .quant and predicted library is still copied."""
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=True)
        with zipfile.ZipFile(zpath) as z:
            kept = sum(i.compress_size for i in z.infolist()
                       if not i.filename.endswith((".quant", ".predicted.speclib")))
            spec = sum(i.compress_size for i in z.infolist()
                       if i.filename.endswith(".predicted.speclib"))
        cap = kept + spec / 2              # fits only if the predicted library is not counted
        self.assertLess(kept, cap)
        self.assertLess(cap, kept + spec)
        self.run_it("analysis-done", "--session", sess, "--out", out,
                    "--zip-cap-gb", str(cap / (1 << 30)))
        rec = self.read(self.only_folder(), "run_record.json")
        self.assertTrue(rec["zip_copy"]["copied"], rec["zip_copy"])

    def test_secrets_are_never_copied(self):
        out = make_search(self.d, secrets=True)
        sess, zpath = make_session(self.d, quant_in_zip=False, zip_secret=True)
        write(os.path.join(sess, "input", "hive.env"), "HIVE_USER=x\n")
        self.run_it("analysis-done", "--session", sess, "--out", out)
        folder = self.only_folder()
        files = self.all_files(folder)
        for bad in ("hive.env", ".pgfarm_token", "slack_webhook.cfg", "notes.log.txt"):
            self.assertFalse([f for f in files if f.endswith(bad)], bad)
        for rel in files:
            p = os.path.join(folder, rel)
            if os.path.islink(p):
                continue
            with open(p, "rb") as fh:
                data = fh.read()
            self.assertNotIn(SECRET.encode(), data, rel)
            self.assertNotIn(b"hooks.slack.com/services", data, rel)
        for shared in ("data_analysis.md", "activity_log.csv"):
            with open(os.path.join(self.runs, shared), "rb") as fh:
                self.assertNotIn(SECRET.encode(), fh.read())
        with zipfile.ZipFile(os.path.join(folder, os.path.basename(zpath))) as z:
            self.assertNotIn("s/input/hive.env", z.namelist())
        log = self.read(folder)
        self.assertIn("credential", log)
        self.assertIn("look like a key, token or webhook", log)

    def test_copy_zip_without_moves_compressed_bytes_intact(self):
        src = os.path.join(self.d, "a.zip")
        with zipfile.ZipFile(src, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("keep/a.txt", "a" * 10000)
            z.writestr(zipfile.ZipInfo("keep/stored.bin"), os.urandom(3000))
            z.writestr("drop/x.quant", os.urandom(5000))
        dst = os.path.join(self.d, "b", "a.zip")
        kept = record_run.copy_zip_without(src, dst, ["drop/x.quant"], record_run.Deadline(60))
        self.assertEqual(kept, 2)
        with zipfile.ZipFile(dst) as z, zipfile.ZipFile(src) as o:
            self.assertIsNone(z.testzip())
            self.assertEqual(z.namelist(), ["keep/a.txt", "keep/stored.bin"])
            for n in z.namelist():
                self.assertEqual(z.read(n), o.read(n))
                self.assertEqual(z.getinfo(n).compress_type, o.getinfo(n).compress_type)


class Gate(Base):
    def test_non_writable_destination_is_skipped(self):
        if os.geteuid() == 0:
            self.skipTest("root can write anywhere")
        out = make_search(self.d)
        os.chmod(self.runs, 0o550)            # the group folder of a group you are not in
        res = self.run_it("search-done", "--out", out)
        os.chmod(self.runs, 0o755)
        self.assertEqual((res["recorded"], res["reason"]), (False, "not_core_member"))
        self.assertIn("Proteomics Core", res["detail"])
        self.assertEqual(os.listdir(self.runs), [])

    def test_no_hive_and_no_login(self):
        res = self.run_it("search-done", "--out", make_search(self.d), route="none")
        self.assertEqual((res["recorded"], res["reason"]), (False, "not_on_hive"))

    def test_off_switches_touch_nothing_and_send_nothing(self):
        out = make_search(self.d)
        for extra in ({"RECORD_RUN": "off"}, {"SKILL_RUNS_DIR": "off"}):
            for route in ("direct", "ssh"):
                env = dict(extra)
                res = self.run_it("search-done", "--out", out, route=route, **env)
                self.assertEqual((res["recorded"], res["reason"]), (False, "disabled"))
        res = self.run_it("analysis-done", "--session", self.d, route="ssh", RECORD_RUN="off")
        self.assertEqual(res["reason"], "disabled")
        self.assertEqual(os.listdir(self.runs), [])
        self.assertFalse(os.path.exists(self.relay_log), "the SSH relay was invoked")

    def test_never_fatal(self):
        f = write(os.path.join(self.d, "not_a_dir"), "x")
        for args, reason in ((["search-done", "--out", f], "bad_input"),
                             (["search-done", "--out", "/no/such/dir"], "out_not_found"),
                             (["analysis-done"], "bad_input"),
                             (["analysis-done", "--session", "/no/such"], "session_not_found")):
            res = self.run_it(*args)
            self.assertEqual((res["recorded"], res["reason"]), (False, reason), args)
        self.assertEqual(self.run_raw("search-done", "--bogus").returncode, 2)   # usage only
        self.assertEqual(os.listdir(self.runs), [])

    def test_stays_under_the_hooks_60_s(self):
        a = record_run.build_parser().parse_args(["search-done", "--out", "x"])
        self.assertLessEqual(a.timeout + 5, 55)

    def test_where(self):
        self.assertTrue(self.run_raw("--where").stdout.startswith("direct:"))
        self.assertTrue(self.run_raw("--where", route="ssh").stdout.startswith("ssh: gabrig@hive"))
        self.assertTrue(self.run_raw("--where", route="none").stdout.startswith("not recorded"))


class Sacct(Base):
    def fake_sacct(self, rows):
        p = write(os.path.join(self.d, "bin", "sacct"),
                  "#!/usr/bin/env bash\ncat <<'EOF'\n" + "\n".join(rows) + "\nEOF\n")
        os.chmod(p, 0o755)
        return p

    def test_times_and_state_come_from_sacct(self):
        out = make_search(self.d, report=False)
        exe = self.fake_sacct([
            "23978018|diann_libpred|COMPLETED|0:0|2026-09-23T17:22:43|2026-09-23T17:22:45|"
            "2026-09-23T17:28:38|00:05:53|hive-as-11-3-71|high|genome-center-grp|u",
            "23978112|diann_search|OUT_OF_MEMORY|0:125|2026-09-23T17:31:34|2026-09-23T17:31:35|"
            "2026-09-23T17:42:01|00:10:26|hive-dc-7-4-18|high|genome-center-grp|u"])
        self.run_it("search-done", "--out", out, RECORD_RUN_SACCT=exe)
        folder = self.only_folder()
        s = self.read(folder, "run_record.json")["search"]
        self.assertEqual((s["status"], s["exit_code"], s["failing_step"]),
                         ("failed", "0:125", "diann_search"))
        self.assertIn("OUT_OF_MEMORY", s["status_source"])
        self.assertEqual((s["times"]["submitted"], s["times"]["finished"]),
                         ("2026-09-23T17:22:43", "2026-09-23T17:42:01"))
        _, rows = self.activity()
        sub = [r for r in rows if r[2] == "search_submitted"][0]
        self.assertTrue(sub[0].startswith("2026-09-23T17:22"), sub)


class ListRuns(Base):
    def put(self, root, user, date, status, name):
        write(os.path.join(root, "sessions", name, "run_record.json"), json.dumps({
            "schema_version": 2, "date": date, "user": user, "name": name,
            "search": {"status": status, "engine": "diann", "engine_version": "2.7.0",
                       "data": {"n_files": 3}, "results": {"median_precursors": 28563}}}))

    def test_list_aggregates_and_filters(self):
        self.put(self.runs, "alice", "2026-08-30", "completed", "2026-08-30_old")
        self.put(self.runs, "alice", "2026-09-23", "failed", "2026-09-23_broken")
        self.put(self.runs, "bob", "2026-09-24", "completed", "2026-09-24_good")
        rows = json.loads(self.run_raw("list", "--json").stdout)
        self.assertEqual(len(rows), 3)
        self.assertEqual({r["session"] for r in json.loads(
            self.run_raw("list", "--json", "--user", "alice").stdout)},
            {"2026-08-30_old", "2026-09-23_broken"})
        self.assertEqual([r["session"] for r in json.loads(
            self.run_raw("list", "--json", "--status", "failed").stdout)], ["2026-09-23_broken"])
        self.assertEqual({r["session"] for r in json.loads(
            self.run_raw("list", "--json", "--since", "2026-09-01").stdout)},
            {"2026-09-23_broken", "2026-09-24_good"})
        tsv = self.run_raw("list", "--tsv", "--user", "bob").stdout.splitlines()
        self.assertEqual(tsv[0].split("\t")[:3], ["date", "session", "user"])
        self.assertEqual(tsv[1].split("\t")[:3], ["2026-09-24", "2026-09-24_good", "bob"])
        self.assertIn("(3 run(s)", self.run_raw("list").stdout)

    def test_list_over_ssh(self):
        self.put(os.path.join(self.remote, "skill_runs"), "gabrig", "2026-09-23", "completed",
                 "2026-09-23_remote")
        rows = json.loads(self.run_raw("list", "--json", route="ssh").stdout)
        self.assertEqual([r["session"] for r in rows], ["2026-09-23_remote"])


class SshRoute(Base):
    """A laptop in hive_remote mode: the record is written ON HIVE through hive_exec.sh."""

    def setUp(self):
        super().setUp()
        self.rruns = os.path.join(self.remote, "skill_runs")
        os.makedirs(self.rruns)
        os.makedirs(os.path.join(self.remote, "skill_issues"))

    def test_search_that_lives_on_hive(self):
        out = make_search(os.path.join(self.remote, "searches"))
        res = self.run_it("search-done", "--out", out.replace(self.remote, FAKE_ROOT),
                          route="ssh")
        self.assertTrue(res["recorded"], res)
        self.assertTrue(res["path"].startswith("gabrig@hive:"))
        folder = self.only_folder(self.rruns)
        rec = self.read(folder, "run_record.json")
        self.assertEqual(rec["search"]["engine_version"], "2.7.0")
        self.assertEqual(rec["history"][-1]["route"], "ssh")
        self.assertIn("output/search/search_provenance.json", self.all_files(folder))
        self.assertTrue(os.path.isfile(os.path.join(self.rruns, "data_analysis.md")))
        self.assertEqual(os.listdir(os.path.join(self.d, record_run.UPLOAD_DIR)), [],
                         "the upload must be removed after use")

    def test_login_from_hive_env_file_only(self):
        out = make_search(os.path.join(self.remote, "searches"))
        envfile = write(os.path.join(self.d, "hive.env"),
                        f"HIVE_USER=gabrig\nHIVE_KEY={self.key}\n")
        e = self.env("ssh")
        e.pop("HIVE_USER")
        e.pop("HIVE_KEY")
        e["HIVE_ENV_FILE"] = envfile
        p = subprocess.run([sys.executable, SCRIPT, "search-done", "--out",
                            out.replace(self.remote, FAKE_ROOT)], capture_output=True,
                           text=True, env=e, timeout=120)
        self.assertTrue(json.loads(p.stdout)["recorded"], p.stdout + p.stderr)

    def test_local_session_with_the_search_on_hive(self):
        out = make_search(os.path.join(self.remote, "searches"))
        sess, zpath = make_session(os.path.join(self.d, "laptop"),
                                   out=out.replace(self.remote, FAKE_ROOT), quant_in_zip=True)
        res = self.run_it("analysis-done", "--session", sess, "--prot", "807", route="ssh")
        self.assertTrue(res["recorded"], res)
        folder = self.only_folder(self.rruns)
        self.assertEqual(os.path.basename(folder), SESSION)
        rec = self.read(folder, "run_record.json")
        # the search was read on HIVE and keyed on its HIVE path; the session came from the laptop
        self.assertEqual(rec["identity"]["out"], os.path.realpath(out))
        self.assertEqual(rec["search"]["status"], "completed")
        self.assertEqual(rec["analysis"]["de"]["significant_per_contrast"], {"B-A": 42})
        self.assertEqual(rec["prot"]["prot"], "PROT_0807")
        self.assertEqual(rec["user"], USER)      # owner of the search folder, not the laptop login
        files = self.all_files(folder)
        self.assertIn("output/HeLa50_Report.docx", files)
        with zipfile.ZipFile(os.path.join(folder, os.path.basename(zpath))) as z:
            self.assertFalse([n for n in z.namelist() if n.endswith(".quant")])
            self.assertIn("s/MANIFEST.txt", z.namelist())            # C4 on the SSH route

    def test_a_non_core_account_over_ssh_is_refused(self):
        if os.geteuid() == 0:
            self.skipTest("root can write anywhere")
        out = make_search(os.path.join(self.remote, "searches"))
        os.chmod(self.rruns, 0o550)
        res = self.run_it("search-done", "--out", out.replace(self.remote, FAKE_ROOT),
                          route="ssh")
        os.chmod(self.rruns, 0o755)
        self.assertEqual((res["recorded"], res["reason"]), (False, "not_core_member"))
        self.assertEqual(self.folders(self.rruns), [])

    def test_ssh_failure_is_a_reason_not_an_error(self):
        broken = write(os.path.join(self.d, "broken_exec.sh"),
                       "#!/usr/bin/env bash\necho 'ssh: connect to host hive: timed out' >&2\n"
                       "exit 255\n")
        os.chmod(broken, 0o755)
        res = self.run_it("search-done", "--out", FAKE_ROOT + "/x", route="ssh",
                          HIVE_EXEC=broken)
        self.assertEqual((res["recorded"], res["reason"]), (False, "ssh_failed"))
        self.assertIn("timed out", res["detail"])


class ReviewFixes(Base):
    """The independent review's fix list (C1-C6, P1-P4), one test per item."""

    def ssh_setup(self):
        self.rruns = os.path.join(self.remote, "skill_runs")
        os.makedirs(self.rruns, exist_ok=True)
        os.makedirs(os.path.join(self.remote, "skill_issues"), exist_ok=True)

    def relay_calls(self):
        return open(self.relay_log).read().splitlines() if os.path.exists(self.relay_log) else []

    def upload_dir(self):
        return os.path.join(self.d, record_run.UPLOAD_DIR)

    # ---- C1 (a): the gate is probed before anything is uploaded
    def test_ssh_non_core_uploads_nothing(self):
        if os.geteuid() == 0:
            self.skipTest("root can write anywhere")
        self.ssh_setup()
        os.chmod(self.rruns, 0o550)
        sess, _ = make_session(os.path.join(self.d, "laptop"), quant_in_zip=True)
        res = self.run_it("analysis-done", "--session", sess, route="ssh")
        os.chmod(self.rruns, 0o755)
        self.assertEqual((res["recorded"], res["reason"]), (False, "not_core_member"))
        calls = self.relay_calls()
        self.assertEqual(len(calls), 1, calls)                  # the probe, and nothing else
        self.assertNotIn("--put", calls[0])
        self.assertFalse(os.path.exists(self.upload_dir()))

    # ---- C1 (f): a dry run over SSH uploads nothing
    def test_ssh_dry_run_uploads_nothing(self):
        self.ssh_setup()
        sess, _ = make_session(os.path.join(self.d, "laptop"), quant_in_zip=True)
        p = self.run_raw("analysis-done", "--session", sess, "--dry-run", route="ssh")
        res = json.loads(p.stdout)
        self.assertEqual((res["recorded"], res["reason"]), (False, "dry_run"))
        self.assertIn("would upload", p.stderr)
        self.assertEqual(len(self.relay_calls()), 1)
        self.assertFalse(os.path.exists(self.upload_dir()))
        self.assertEqual(os.listdir(self.rruns), [])

    # ---- C1 (b): a zip that cannot be staged is a note, and the record is still written
    def test_ssh_zip_staging_failure_still_records(self):
        self.ssh_setup()
        sess, zpath = make_session(os.path.join(self.d, "laptop"), quant_in_zip=True)
        raw = bytearray(open(zpath, "rb").read())
        with zipfile.ZipFile(zpath) as z:
            kept = [i for i in z.infolist() if not i.filename.endswith(".quant")][0]
        raw[kept.header_offset:kept.header_offset + 4] = b"XXXX"   # central directory still valid
        open(zpath, "wb").write(bytes(raw))
        res = self.run_it("analysis-done", "--session", sess, "--prot", "807", route="ssh")
        self.assertTrue(res["recorded"], res)
        self.assertFalse(res["zip_copied"])
        rec = self.read(self.only_folder(self.rruns), "run_record.json")
        self.assertIn("copy failed", rec["zip_copy"]["reason"])
        self.assertIn("the session zip is on the user's machine", self.read(
            self.only_folder(self.rruns)))

    # ---- C1 (c)+(d): the upload is removed even when the merge runs out of time
    def test_ssh_upload_removed_when_the_merge_times_out(self):
        self.ssh_setup()
        slow = os.path.join(self.d, "slow_exec.sh")
        body = open(self.fake_exec).read().replace(
            'fi\ncmd=', 'fi\ncase "$1" in *record_run.py\\ merge*) sleep 30 ;; esac\ncmd=', 1)
        self.assertIn("sleep 30", body)
        write(slow, body)
        os.chmod(slow, 0o755)
        sess, _ = make_session(os.path.join(self.d, "laptop"), quant_in_zip=True)
        t0 = time.monotonic()
        res = self.run_it("analysis-done", "--session", sess, "--timeout", "8", route="ssh",
                          HIVE_EXEC=slow)
        # The HIVE call is cut at its own timeout (about 2 s here) because the whole process group
        # is killed; without that, the remote `sleep` holds the pipes open and only the SIGALRM
        # backstop (timeout + 5 = 13 s) ends it.
        self.assertLess(time.monotonic() - t0, 8)
        self.assertEqual((res["recorded"], res["reason"]), (False, "timeout"))
        self.assertEqual(os.listdir(self.upload_dir()), [])

    def test_merge_sweeps_uploads_older_than_a_day(self):
        self.ssh_setup()
        old = os.path.join(self.upload_dir(), "record_run_old")
        fresh = os.path.join(self.upload_dir(), "record_run_fresh")
        for d in (old, fresh):
            os.makedirs(os.path.join(d, "files"))
        t = time.time() - 2 * 86400
        os.utime(old, (t, t))
        out = make_search(os.path.join(self.remote, "searches"))
        res = self.run_it("search-done", "--out", out.replace(self.remote, FAKE_ROOT), route="ssh")
        self.assertTrue(res["recorded"], res)
        self.assertEqual(sorted(os.listdir(self.upload_dir())), ["record_run_fresh"])

    # ---- C1 (e): over the upload limit the zip stays on the laptop, identified
    def test_ssh_zip_over_the_upload_limit_stays_with_its_checksum(self):
        import hashlib
        import socket
        self.ssh_setup()
        sess, zpath = make_session(os.path.join(self.d, "laptop"), quant_in_zip=True)
        res = self.run_it("analysis-done", "--session", sess, route="ssh",
                          RECORD_RUN_SSH_ZIP_CAP_MB="0.0001")
        self.assertTrue(res["recorded"], res)
        folder = self.only_folder(self.rruns)
        zc = self.read(folder, "run_record.json")["zip_copy"]
        self.assertEqual((zc["copied"], zc["on_host"], zc["local_path"], zc["local_bytes"]),
                         (False, socket.gethostname(), zpath, os.path.getsize(zpath)))
        self.assertEqual(zc["sha256"], hashlib.sha256(open(zpath, "rb").read()).hexdigest())
        self.assertFalse([f for f in os.listdir(folder) if f.endswith(".zip")])
        log = self.read(folder)
        self.assertIn("the session zip is on the user's machine", log)
        self.assertIn(zc["sha256"], log)

    def test_the_sigalrm_backstop_follows_the_timeout(self):
        import contextlib
        calls = []
        real_alarm = record_run.signal.alarm
        saved = dict(os.environ)
        record_run.signal.alarm = lambda n: calls.append(n) or 0
        try:
            os.environ.update(self.env("direct"))
            os.environ.pop("RECORD_RUN", None)
            with contextlib.redirect_stdout(io.StringIO()):
                record_run.main(["analysis-done", "--session", "/no/such", "--timeout", "300"])
        finally:
            record_run.signal.alarm = real_alarm
            os.environ.clear()
            os.environ.update(saved)
        self.assertEqual((calls[0], calls[-1]), (305, 0))

    # ---- C3: FASTAs never reach the registry's copy of the zip
    def test_fastas_are_left_out_and_counted(self):
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=True, fasta_in_zip=True)
        self.run_it("analysis-done", "--session", sess, "--out", out)
        folder = self.only_folder()
        with zipfile.ZipFile(os.path.join(folder, os.path.basename(zpath))) as z:
            self.assertFalse([n for n in z.namelist() if record_run.FASTA_RE.search(n)])
        f = self.read(folder, "run_record.json")["findings"][0]
        self.assertEqual((f["id"], f["n_quant"], f["n_predicted_speclib"], f["n_fasta"]),
                         ("session_zip_trimmed", 4, 1, 2))
        self.assertIn("4 per-run .quant files, 1 predicted spectral library and 2 FASTA files",
                      f["detail"])
        self.assertIn("sidecar", f["detail"])

    def test_a_fasta_only_zip_is_trimmed_without_a_note(self):
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=False, fasta_in_zip=True)
        res = self.run_it("analysis-done", "--session", sess, "--out", out)
        self.assertEqual(res["findings"], [])                  # a FASTA is not an anomaly
        folder = self.only_folder()
        with zipfile.ZipFile(os.path.join(folder, os.path.basename(zpath))) as z:
            self.assertFalse([n for n in z.namelist() if record_run.FASTA_RE.search(n)])
        self.assertIn("2 FASTA", self.read(folder))

    # ---- C4: the copy carries a MANIFEST.txt
    def test_manifest_is_added_on_the_plain_copy_path_too(self):
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=False)        # nothing to drop
        self.run_it("analysis-done", "--session", sess, "--out", out)
        folder = self.only_folder()
        with zipfile.ZipFile(os.path.join(folder, os.path.basename(zpath))) as z:
            self.assertEqual(z.read("s/MANIFEST.txt").decode(),
                             open(os.path.join(sess, "MANIFEST.txt")).read())
            self.assertIsNone(z.testzip())
        self.assertTrue(self.read(folder, "run_record.json")["zip_copy"]["manifest_added"])
        with zipfile.ZipFile(zpath) as z:                               # the original untouched
            self.assertNotIn("s/MANIFEST.txt", z.namelist())

    def test_a_zip_that_has_its_manifest_keeps_it(self):
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=True, manifest_in_zip=True)
        self.run_it("analysis-done", "--session", sess, "--out", out)
        with zipfile.ZipFile(os.path.join(self.only_folder(), os.path.basename(zpath))) as z:
            self.assertEqual(z.namelist().count("s/MANIFEST.txt"), 1)
            self.assertEqual(z.read("s/MANIFEST.txt"), b"zipped manifest\n")

    # ---- C5: a second, different failure gets its own dated line
    def test_a_different_failure_is_logged_again(self):
        out = make_search(self.d, report=False)
        for st, code, step in (("failed", 1, "libpred"), ("failed", 137, "search"),
                               ("failed", 137, "search"), ("completed", 0, "search")):
            self.run_it("search-done", "--out", out, "--status", st, "--exit-code", str(code),
                        "--step", step)
        md = self.master()
        self.assertIn("FAILED (exit 1, step libpred)", md)
        self.assertEqual(md.count("search failed (exit 137, step search)"), 1)   # not repeated
        self.assertEqual(md.count(": search completed"), 1)
        self.assertEqual(md.count("\n## "), 1)
        notes = [r[6] for r in self.activity()[1] if r[2] == "search_failed"]
        self.assertEqual([("re-recorded" in n) for n in notes], [False, False, True])

    # ---- C6: copies cannot eat phase C's lock wait; the message states the real wait
    def test_copies_leave_phase_c_time_for_its_locks(self):
        """A zip copy that takes all the time it is given must still leave phase C's locks a real
        wait (not a single try that falls straight through to an unlocked write)."""
        import contextlib
        sess, _ = make_session(os.path.join(self.d, "laptop"), quant_in_zip=False)
        lk = os.path.join(self.runs, "data_analysis.md.lock.d")
        os.makedirs(lk)
        write(os.path.join(lk, "owner"), "othernode:123:1\n")         # held throughout
        real_copy, saved = record_run.copy_file, dict(os.environ)

        def slow_copy(src, dst, deadline):
            if src.endswith(".zip"):                                    # uses its whole budget
                while deadline.left() > 0.2:
                    time.sleep(0.05)
                raise TimeoutError("the --timeout was reached")
            return real_copy(src, dst, deadline)
        record_run.copy_file = slow_copy
        del record_run.LOCK_TIMEOUTS[:]
        out, err = io.StringIO(), io.StringIO()
        try:
            os.environ.update(self.env("direct"), RECORD_RUN_LOCK_WAIT="4")
            os.environ.pop("RECORD_RUN", None)
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
                record_run.main(["analysis-done", "--session", sess, "--timeout", "14"])
        finally:
            record_run.copy_file = real_copy
            os.environ.clear()
            os.environ.update(saved)
        res = json.loads(out.getvalue())
        self.assertTrue(res["recorded"], res)
        self.assertEqual(res["lock_timeout"], ["data_analysis.md"])
        waited = re.search(r"data_analysis\.md\.lock\.d still held after (\d+\.\d) s",
                           err.getvalue())
        self.assertIsNotNone(waited, err.getvalue())
        self.assertGreaterEqual(float(waited.group(1)), 3.0)          # a real wait, not one try

    def test_the_timeout_message_states_the_real_wait(self):
        import contextlib
        path = os.path.join(self.d, "w.csv")
        os.mkdir(path + ".lock.d")
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            with record_run.DirLock(path, record_run.Deadline(4)):       # wait 10, but 1 s left
                pass
        m = re.search(r"still held after (\d+\.\d) s", err.getvalue())
        self.assertIsNotNone(m, err.getvalue())
        self.assertLess(float(m.group(1)), 3)

    # ---- P2: a lock that shows life while being broken is put back
    def test_a_lock_that_shows_life_during_the_break_is_restored(self):
        import contextlib
        path = os.path.join(self.d, "p2.csv")
        lkd = path + ".lock.d"
        os.mkdir(lkd)
        write(os.path.join(lkd, "owner"), "A:slow:1\n")
        os.utime(lkd, (time.time() - 3600,) * 2)
        b = record_run.DirLock(path, record_run.Deadline(30), wait=0.5)

        def a_shows_life(judged):
            b.break_hook = None
            os.utime(lkd, None)                          # A is alive after all: fresh again
        b.break_hook = a_shows_life
        err = io.StringIO()
        with contextlib.redirect_stderr(err):
            with b:
                self.assertFalse(b.held)
        self.assertEqual(record_run.DirLock.owner_of(lkd), "A:slow:1")
        self.assertIn("restored", err.getvalue())
        self.assertEqual([f for f in os.listdir(self.d) if f.startswith("p2.csv.lock.d.")], [])

    # ---- P3: search-done reads the FRAN receipt as it is now
    def test_search_done_reads_the_fran_receipt_fresh(self):
        out = make_search(self.d)
        receipt = os.path.join(out, "fran_deposit.json")
        body = open(receipt).read()
        os.remove(receipt)
        self.run_it("search-done", "--out", out, "--status", "completed", "--exit-code", "0")
        self.assertIn("no fran_deposit.json", self.read(self.only_folder()))
        write(receipt, body)                              # the hook staged it after all
        self.run_it("search-done", "--out", out, "--status", "completed", "--exit-code", "0")
        log = self.read(self.only_folder())
        self.assertIn("- **staged** ->", log)
        self.assertIn("fran_staged", [r[2] for r in self.activity()[1]])

    # ---- P4: an ingested receipt is logged as fran_ingested
    def test_an_ingested_receipt_is_logged_as_ingested(self):
        out = make_search(self.d)
        r = json.load(open(os.path.join(out, "fran_deposit.json")))
        r["status"] = "ingested"
        write(os.path.join(out, "fran_deposit.json"), json.dumps(r))
        self.run_it("search-done", "--out", out)
        rows = self.activity()[1]
        self.assertIn("fran_ingested", [x[2] for x in rows])
        self.assertNotIn("fran_staged", [x[2] for x in rows])

    # ---- C7: findings are replaced per part, and only for the parts a call re-evaluated
    def test_search_done_keeps_the_zip_findings(self):
        out = make_search(self.d)
        sess, _ = make_session(self.d, quant_in_zip=True)
        self.run_it("analysis-done", "--session", sess, "--out", out)
        res = self.run_it("search-done", "--out", out, "--status", "completed", "--exit-code", "0")
        self.assertEqual(res["findings"], ["session_zip_trimmed"])   # search-done: no zip look
        rec = self.read(self.only_folder(), "run_record.json")
        self.assertEqual([f["part"] for f in rec["findings"]], ["zip"])
        self.assertIn(".quant files", self.read(self.only_folder()))

    def test_a_legacy_untagged_finding_is_replaced_by_its_part(self):
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=False)             # a clean zip
        self.run_it("analysis-done", "--session", sess, "--out", out)
        folder = self.only_folder()
        path = os.path.join(folder, "run_record.json")
        rec = json.load(open(path))
        rec["findings"] = [{"id": "session_zip_contains_quant", "detail": "legacy: .quant held"}]
        write(path, json.dumps(rec))                     # written before findings had a part
        self.run_it("search-done", "--out", out)         # does not look at the zip: kept
        self.assertEqual([f["id"] for f in self.read(folder, "run_record.json")["findings"]],
                         ["session_zip_contains_quant"])
        self.run_it("analysis-done", "--session", sess, "--out", out)    # re-evaluates the zip
        self.assertEqual(self.read(folder, "run_record.json")["findings"], [])
        self.assertNotIn("legacy: .quant held", self.read(folder))

    def test_an_untagged_finding_of_unknown_part_is_replaced_by_any_evaluating_call(self):
        out = make_search(self.d)
        self.run_it("search-done", "--out", out)
        path = os.path.join(self.only_folder(), "run_record.json")
        rec = json.load(open(path))
        rec["findings"] = [{"id": "mystery", "detail": "who raised this?"}]
        write(path, json.dumps(rec))
        self.run_it("search-done", "--out", out)
        self.assertEqual(self.read(self.only_folder(), "run_record.json")["findings"], [])

    # ---- a re-finalized, clean zip clears the old zip's note
    def test_a_clean_refinalize_clears_the_zip_note(self):
        out = make_search(self.d)
        sess, zpath = make_session(self.d, quant_in_zip=True)
        self.assertEqual(self.run_it("analysis-done", "--session", sess, "--out", out)
                         ["findings"], ["session_zip_trimmed"])
        with zipfile.ZipFile(zpath) as z:
            keep = [(i, z.read(i)) for i in z.infolist()
                    if not i.filename.endswith((".quant", ".predicted.speclib"))]
        os.remove(zpath)
        with zipfile.ZipFile(zpath, "w") as z:
            for i, b in keep:
                z.writestr(i, b)
        self.assertEqual(self.run_it("analysis-done", "--session", sess, "--out", out)
                         ["findings"], [])
        dq = self.read(self.only_folder()).split("## Data Quality Notes", 1)[1].split("\n## ")[0]
        self.assertNotIn(".quant files", dq)




class FindingPartsAreDeclared(unittest.TestCase):
    """merge() replaces findings by the PART that raised them. A finding with no part is replaced
    by ANY call that re-evaluates something, and a misspelled part is never replaced -- so every
    finding the code raises must carry a part from FINDING_PARTS."""

    def test_every_finding_literal_names_a_declared_part(self):
        src = open(record_run.__file__, encoding="utf-8").read()
        parts = re.findall(r'"part":\s*"([^"]+)"', src)
        self.assertTrue(parts, "no finding declares a part")
        self.assertEqual(sorted(set(parts) - set(record_run.FINDING_PARTS)), [])
        # every finding dict (an "id" that ends in a finding) carries a "part" beside it
        for m in re.finditer(r'\{\s*"id":\s*"(session_[a-z_]+)"(.{0,200})', src, re.S):
            self.assertIn('"part":', m.group(2), m.group(1))


class RunBoundedOnWindows(unittest.TestCase):
    """gabrig's laptop is Windows: no os.killpg, no signal.SIGKILL. A timed-out call must still be
    killed (tree and all) and still raise TimeoutExpired -- not AttributeError, which the SSH
    route reported as "error" while the upload kept running."""

    def test_timeout_without_killpg_kills_the_child_and_raises_timeout(self):
        calls = []
        real_run = subprocess.run

        def fake_run(argv, *a, **k):
            if argv and argv[0] == "taskkill":
                calls.append(argv)
                return subprocess.CompletedProcess(argv, 0)
            return real_run(argv, *a, **k)
        saved = record_run.os.__dict__.pop("killpg", None)
        try:
            with mock.patch.object(record_run.subprocess, "run", fake_run):
                t0 = time.monotonic()
                with self.assertRaises(subprocess.TimeoutExpired):
                    record_run.run_bounded([sys.executable, "-c", "import time; time.sleep(30)"], 1)
                self.assertLess(time.monotonic() - t0, 15)
        finally:
            if saved is not None:
                record_run.os.killpg = saved
        self.assertEqual(len(calls), 1)
        self.assertEqual(calls[0][:4], ["taskkill", "/T", "/F", "/PID"])

    def test_timeout_with_killpg_still_raises_timeout(self):
        t0 = time.monotonic()
        with self.assertRaises(subprocess.TimeoutExpired):
            record_run.run_bounded(["bash", "-c", "sleep 30 & sleep 30"], 1)
        self.assertLess(time.monotonic() - t0, 15)

if __name__ == "__main__":
    unittest.main()
