#!/usr/bin/env python3
"""
run_search.py --engine metamorpheus: MetaMorpheus as a DDA engine for identifications only.

What the route promises, each pinned here:
  * resolve_defaults.py --engine metamorpheus ships the VALIDATED task TOML (its md5 checked)
    as wf/metamorpheus.toml, with a manifest naming engine metamorpheus 1.1.12 and the validation
    it came from; DIA, an ion-trap MS2 and --ms1-ppm/--ms2-ppm are refused, not silently applied;
  * generation (--sbatch) writes a job that checks the pinned image's sha256, copies the image's
    app AND the database to node-local $TMPDIR (MetaMorpheus writes into its own folder and a
    DatabaseIndex/ beside its database), installs the custom modification files, runs
    `dotnet CMD.dll --acceptThermoLicence -t <task copy> -d <db copy> -s <files> -o <out>/run`,
    checks AllPSMs.psmtsv was written, and builds identifications.tsv;
  * the task copy the job runs has MaxThreadsToUsePerFile = --threads (the given file says 16);
  * search_provenance.json records engine metamorpheus, version 1.1.12 from tools.json, the
    given TOML's md5 and validated preset, the custom mods, quantified false, and a job-end hook
    that never stages for FRAN; sage_lfq_check.json is the not_applicable record (DE refused);
  * refusals before anything is written: DIA, a .d, no image, not a .toml, a Custom mod no
    installed file defines, an earlier run/ in --out;
  * the adapter keeps target PSMs with QValue AND QValue Notch <= 0.01 (MetaMorpheus' own 1%
    rule: on a real search it reproduces results.txt's PSM and peptide totals exactly) whose
    peptide passes the same rule in AllPeptides, one row per (run, Full Sequence, charge);
  * the job, run with a stub apptainer, ends with identifications.tsv and its record.
"""
import csv
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import metamorpheus_params as mmp  # noqa: E402
import resolve_defaults  # noqa: E402
from job_env import job_env  # noqa: E402

PRESET = os.path.join(SKILL, "presets", "metamorpheus", "NonSpecific_Peptidomics_HCD.toml")
PIN = resolve_defaults.METAMORPHEUS_PIN

COLS = ["File Name", "Scan Number", "Precursor Charge", "Base Sequence", "Full Sequence",
        "Accession", "Decoy/Contaminant/Target", "QValue", "QValue Notch", "PEP_QValue"]
OX = "M[Common Variable:Oxidation on M]ASTER"
PSMS = [  # File Name, scan, z, base, full, accession, DCT, q, q notch, pep_q
    ("A", 1, 2, "PEPTIDEK", "PEPTIDEK", "P1", "T", 0.001, 0.001, 0.002),
    ("A", 2, 2, "PEPTIDEK", "PEPTIDEK", "P1", "T", 0.005, 0.005, 0.002),   # same precursor
    ("A", 3, 3, "PEPTIDEK", "PEPTIDEK", "P1", "T", 0.002, 0.001, 0.002),   # another charge
    ("B", 4, 2, "MASTER", OX, "P2", "T", 0.003, 0.003, 0.004),
    ("B", 5, 2, "DECOYPEP", "DECOYPEP", "DECOY_P3", "D", 0.001, 0.001, 0.001),
    ("B", 6, 2, "LOWQ", "LOWQ", "P4", "T", 0.05, 0.05, 0.05),
    ("B", 9, 2, "NOTCHQ", "NOTCHQ", "P7", "T", 0.004, 0.03, 0.004),       # fails in its notch
    ("A", 7, 2, "NOPEPQ", "NOPEPQ", "P5", "T", 0.004, 0.004, 0.004),      # peptide fails q
    ("A", 8, 2, "CONTAM", "CONTAM", "Cont_P6", "C", 0.001, 0.001, 0.001),
]
PEPTIDES = [
    ("A", 1, 2, "PEPTIDEK", "PEPTIDEK", "P1", "T", 0.001, 0.001, 0.002),
    ("B", 4, 2, "MASTER", OX, "P2", "T", 0.003, 0.003, 0.004),
    ("B", 9, 2, "NOTCHQ", "NOTCHQ", "P7", "T", 0.004, 0.004, 0.004),
    ("A", 7, 2, "NOPEPQ", "NOPEPQ", "P5", "T", 0.002, 0.02, 0.02),        # fails in its notch
    ("A", 8, 2, "CONTAM", "CONTAM", "Cont_P6", "C", 0.001, 0.001, 0.001),
    ("B", 5, 2, "DECOYPEP", "DECOYPEP", "DECOY_P3", "D", 0.001, 0.001, 0.001),
]

# A stand-in for apptainer: the app copy makes CMD.dll + Mods/; the search writes the fixture's
# tables under its -o, records the -d it was given, and fails if that is not a readable file.
STUB_APPTAINER = r"""#!/bin/bash
[ "$1" = exec ] || exit 2; shift
while [ $# -gt 0 ]; do case "$1" in --cleanenv) shift ;; --bind) shift 2 ;; *) break ;; esac; done
shift   # the image
if [ "$1" = bash ]; then app="${@: -1}"; mkdir -p "$app/Mods"; : > "$app/CMD.dll"; exit 0; fi
out=""; db=""; app=""
while [ $# -gt 0 ]; do
  case "$1" in -o) out="$2"; shift ;; -d) db="$2"; shift ;; *CMD.dll) app="${1%/CMD.dll}" ;; esac
  shift
done
[ -r "$db" ] || { echo "no database $db" >&2; exit 3; }
ls "$app/Mods" > "$MM_STUB_LOG.mods"; echo "$db" > "$MM_STUB_LOG.db"
mkdir -p "$out/Task1SearchTask" "$(dirname "$db")/DatabaseIndex"
cp "$MM_FIXTURE"/All*.psmtsv "$out/Task1SearchTask/"
printf 'MetaMorpheus: version Not a release version\nAll target PSMs with q-value <= 0.01: 6\n' \
  > "$out/allResults.txt"
"""


def _write(path, text, mode=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        fh.write(text)
    if mode:
        os.chmod(path, mode)
    return path


def write_table(path, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(COLS)
        for r in rows:
            w.writerow(r)


class Harness(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.raw = os.path.join(self.d, "raw")
        self.files = [_write(os.path.join(self.raw, "A.raw"), "one\n"),
                      _write(os.path.join(self.raw, "B.mzML"), "two\n")]
        self.dbdir = os.path.join(self.d, "db")
        self.fasta = _write(os.path.join(self.dbdir, "search.fasta"), ">sp|P1|X\nPEPTIDEK\n")
        self.out = os.path.join(self.d, "out")
        self.job = os.path.join(self.d, "mm_job.sh")
        self.wf = os.path.join(self.d, "wf")

    def sif(self, name="metamorpheus_1.1.12.sif"):
        return _write(os.path.join(self.d, "img", name), "image\n")

    def resolve(self, *extra, acquisition="DDA"):
        return subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
             "--acquisition", acquisition, "--instrument", "Orbitrap Exploris 480",
             "--ms1-resolution", "120000", "--ms2-resolution", "15000",
             "--engine", "metamorpheus", "--dest", self.wf, *extra],
            capture_output=True, text=True, timeout=120)

    def generate(self, *extra, files=None, params=None, bundle=None, sif=None, tools=None):
        if bundle is None:
            r = self.resolve()
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            bundle = os.path.join(self.wf, "workflow.manifest.json")
        sif = sif or self.sif()
        tools = tools or _write(os.path.join(self.d, "tools.json"), json.dumps(
            {"metamorpheus": sif, "versions": {"metamorpheus": "1.1.12"}}))
        env = job_env(self.d, PATH=os.pathsep.join([os.path.dirname(sys.executable),
                                                    "/usr/bin", "/bin"]))
        return subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "run_search.py"),
             "--tools", tools, "--bundle", bundle,
             "--params", params or os.path.join(self.wf, "metamorpheus.toml"),
             "--fasta", self.fasta, "--out", self.out, "--files", *(files or self.files),
             "--engine", "metamorpheus", "--threads", "8", "--no-notify", "--not-qc",
             "--sbatch", self.job, *extra],
            capture_output=True, text=True, env=env, cwd=self.d, timeout=240)

    def prov(self):
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            return json.load(fh)

    def script(self):
        with open(self.job) as fh:
            return fh.read()


class ResolveDefaults(Harness):
    def test_the_validated_toml_is_shipped_and_recorded(self):
        r = self.resolve()
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(self.wf, "workflow.manifest.json")) as fh:
            m = json.load(fh)
        self.assertEqual(m["engine"], {"name": "metamorpheus", "version": "1.1.12"})
        self.assertIsNone(m["de"]["method"])
        self.assertFalse(m["search"]["estimate_params"])
        toml = os.path.join(self.wf, "metamorpheus.toml")
        self.assertEqual(m["search"]["params_file"], toml)
        self.assertEqual(mmp.md5(toml), mmp.md5(PRESET))
        self.assertEqual(mmp.md5(PRESET),
                         resolve_defaults.METAMORPHEUS_PRESETS["nonspecific_peptidomics"]["md5"])
        self.assertEqual(m["validated"]["date"], "2026-11-04")
        self.assertEqual((m["search"]["ms1_ppm"], m["search"]["ms2_ppm"]), (10.0, 20.0))
        self.assertEqual(m["search"]["metamorpheus_preset"]["settings"]["search_type"],
                         "NonSpecific")

    def test_dia_is_refused(self):
        r = self.resolve(acquisition="DIA")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("DDA only", r.stderr)
        self.assertFalse(os.path.exists(self.wf))

    def test_a_ppm_override_is_refused_not_ignored(self):
        r = self.resolve("--ms1-ppm", "5")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("not applied to a MetaMorpheus task TOML", r.stderr)

    def test_an_ion_trap_ms2_is_refused(self):
        r = subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
             "--acquisition", "DDA", "--instrument", "Orbitrap Fusion Lumos",
             "--ms1-resolution", "120000", "--ms2-analyzer", "ITMS",
             "--engine", "metamorpheus", "--dest", self.wf],
            capture_output=True, text=True, timeout=120)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("ion trap", r.stderr)


class TaskToml(unittest.TestCase):
    def test_the_shipped_toml_reads_as_validated(self):
        st = mmp.summary(mmp.read_task(PRESET))
        self.assertEqual(st["search_type"], "NonSpecific")
        self.assertEqual(st["precursor_tolerance"], "±10.0000 PPM")
        self.assertEqual(st["q_value_threshold"], "0.01")
        self.assertEqual(st["label_free_quantification"], "false")
        self.assertIn("Custom: Gln to PyroGlu peptide N-term on Q", st["variable_mods"])

    def test_its_custom_mods_are_defined_by_the_shipped_file(self):
        mods = [os.path.join(resolve_defaults.METAMORPHEUS_PRESET_DIR, m) for m in
                resolve_defaults.METAMORPHEUS_PRESETS["nonspecific_peptidomics"]["mods"]]
        self.assertEqual(mmp.problems(mmp.read_task(PRESET), mods), [])
        self.assertEqual(len(mmp.problems(mmp.read_task(PRESET), [])), 2)

    def test_the_run_copy_sets_threads_only(self):
        with open(PRESET, encoding="utf-8", newline="") as fh:
            text = fh.read()
        new, old = mmp.run_copy_text(text, 8)
        self.assertEqual(old, "16")
        diff = [(a, b) for a, b in zip(text.splitlines(), new.splitlines()) if a != b]
        self.assertEqual(diff, [("MaxThreadsToUsePerFile = 16", "MaxThreadsToUsePerFile = 8")])


class Generation(Harness):
    def test_the_job_and_the_record(self):
        p = self.generate()
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        s = self.script()
        self.assertIn(f"echo '{PIN['sif_sha256']}  ", s)                  # sha256 checked first
        self.assertIn('MM_TMP="$(mktemp -d', s)
        self.assertIn('cp ' + self.fasta + ' "$MM_TMP/db/search.fasta"', s)
        self.assertIn("for f in /MetaMorpheus/*", s)
        self.assertIn('= publish ]', s)
        self.assertIn("dotnet \"$MM_APP/CMD.dll\" --acceptThermoLicence -t ", s)
        self.assertIn('-d "$MM_TMP/db/search.fasta"', s)
        self.assertIn(f"-o {os.path.join(self.out, 'run')} -v normal", s)
        self.assertIn(f"--bind {self.raw}:{self.raw}:ro", s)
        self.assertNotIn(f"--bind {self.dbdir}", s)          # the database is copied, not bound
        self.assertIn("AllPSMs.psmtsv", s)
        self.assertIn("run_search.adapt_metamorpheus_ids", s)
        self.assertIn("--no-fran", s)                                    # the hook never stages
        with open(os.path.join(self.out, "metamorpheus_task.toml")) as fh:
            self.assertIn("MaxThreadsToUsePerFile = 8\n", fh.read())
        self.assertTrue(os.path.isfile(os.path.join(self.out, "metamorpheus_mods",
                                                    "CustomPyroGluPeptideNterm.txt")))
        prov = self.prov()
        self.assertEqual(prov["engine"], "metamorpheus")
        self.assertEqual(prov["version"], "1.1.12")
        self.assertIs(prov["engine_version"]["mismatch"], False)
        self.assertIs(prov["quantified"], False)
        self.assertEqual(prov["job_end_hook"]["fran"], "off")
        mm = prov["metamorpheus"]
        self.assertEqual(mm["task_given"]["validated_preset"], "nonspecific_peptidomics")
        self.assertEqual(mm["task_given"]["md5"], mmp.md5(PRESET))
        self.assertEqual(mm["task_run"]["changed"][mmp.THREADS_KEY]["run"], "8")
        self.assertEqual(mm["container"]["sif_sha256"], PIN["sif_sha256"])
        self.assertEqual(mm["database"]["md5"], mmp.md5(self.fasta))
        self.assertEqual([m["installed_as"] for m in mm["custom_mods"]],
                         ["Mods/CustomPyroGluPeptideNterm.txt"])
        with open(os.path.join(self.out, "sage_lfq_check.json")) as fh:
            q = json.load(fh)
        self.assertEqual((q["status"], q["engine"]), ("not_applicable", "metamorpheus"))
        self.assertIn("no DE", q["de_refusal"])

    def test_fran_is_switched_off_even_when_not_asked(self):
        p = self.generate()
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertIn("never staged for FRAN", p.stderr)

    def test_a_toml_that_is_not_the_validated_one_says_so(self):
        r = self.resolve()
        self.assertEqual(r.returncode, 0, r.stderr)
        mine = os.path.join(self.d, "mine.toml")
        with open(PRESET, encoding="utf-8") as fh:
            _write(mine, fh.read().replace("±10.0000 PPM", "±5.0000 PPM"))
        p = self.generate(params=mine, bundle=os.path.join(self.wf, "workflow.manifest.json"))
        self.assertEqual(p.returncode, 0, p.stderr)
        tg = self.prov()["metamorpheus"]["task_given"]
        self.assertIsNone(tg["validated_preset"])
        self.assertIn("not one of the skill's validated", tg["note"])


class Refusals(Harness):
    def refused(self, p, words):
        self.assertNotEqual(p.returncode, 0, p.stdout)
        self.assertIn(words, p.stdout + p.stderr)
        self.assertFalse(os.path.exists(self.job))

    def test_dia(self):
        b = _write(os.path.join(self.d, "dia.json"), json.dumps(
            {"engine": {"name": "metamorpheus"}, "acquisition": "DIA"}))
        self.refused(self.generate(bundle=b), "DDA identifications only")

    def test_a_bruker_d(self):
        d = os.path.join(self.raw, "C.d")
        os.makedirs(d)
        # (the fixture's .d is a bare folder, not a readable timsTOF run)
        self.refused(self.generate("--allow-damaged-tdf", files=[d]), "converted to mzML")

    def test_no_image(self):
        tools = _write(os.path.join(self.d, "t2.json"), json.dumps(
            {"metamorpheus": os.path.join(self.d, "nowhere", "metamorpheus_1.1.12.sif")}))
        self.refused(self.generate(tools=tools), "not an Apptainer image on disk")

    def test_params_not_a_toml(self):
        cfg = _write(os.path.join(self.d, "sage.json"), "{}")
        self.refused(self.generate(params=cfg), "search-task .toml")

    def test_an_undefined_custom_mod(self):
        r = self.resolve()
        self.assertEqual(r.returncode, 0, r.stderr)
        bad = os.path.join(self.d, "bad.toml")
        with open(PRESET, encoding="utf-8") as fh:
            text = fh.read().replace(
                'ListOfModsVariable = "', 'ListOfModsVariable = "Custom\tMy mod on K\t\t')
        _write(bad, text)
        self.refused(self.generate(params=bad,
                                   bundle=os.path.join(self.wf, "workflow.manifest.json")),
                     "'My mod on K'")

    def test_an_earlier_run_folder(self):
        os.makedirs(os.path.join(self.out, "run"))
        self.refused(self.generate(), "from an earlier MetaMorpheus search")


class Methods(unittest.TestCase):
    """make_methods.search_record reads a MetaMorpheus search from its task TOML and record:
    non-specific cleavage, length, tolerances, the record's FDR sentence; its modifications are
    MetaMorpheus' own names, so they go to the writer as a warning, never mapped to UNIMOD."""

    def test_the_search_paragraph(self):
        import make_methods as mm
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        p = _write(os.path.join(d, "search_provenance.json"), json.dumps(
            {"engine": "metamorpheus", "version": "1.1.12", "params_file": PRESET,
             "metamorpheus_identifications": {"methods_sentence": "Identifications were Y."}}))
        r = mm.search_record(search_prov=p)
        self.assertEqual(r["engine_label"], "MetaMorpheus")
        self.assertEqual(r["cleavage"]["specificity"], "none")
        self.assertEqual(r["pep_len"]["value"], (7, 25))
        self.assertEqual((r["ms1_tol"]["value"], r["ms2_tol"]["value"]), (10.0, 20.0))
        self.assertFalse(r["mods_read"])
        self.assertTrue(any("PyroGlu" in w for w in r["warnings"]))
        text = mm.search_paragraph(r)
        self.assertIn("MetaMorpheus 1.1.12", text)
        self.assertIn("identifications only", text)
        self.assertIn("non-specific", text)
        self.assertIn("Identifications were Y.", text)


class Adapter(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        task = os.path.join(self.d, "run", "Task1SearchTask")
        write_table(os.path.join(task, "AllPSMs.psmtsv"), PSMS)
        write_table(os.path.join(task, "AllPeptides.psmtsv"), PEPTIDES)

    def test_the_kept_rows(self):
        import run_search
        path = run_search.adapt_metamorpheus_ids(self.d)
        with open(path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        self.assertEqual(list(rows[0]), list(run_search.METAMORPHEUS_IDS_COLUMNS))
        got = [(r["run"], r["peptide"], r["charge"], r["psm_q"], r["target_type"]) for r in rows]
        self.assertEqual(got, [("A", "CONTAM", "2", "0.001", "C"),
                               ("A", "PEPTIDEK", "2", "0.001", "T"),
                               ("A", "PEPTIDEK", "3", "0.002", "T"),
                               ("B", OX, "2", "0.003", "T")])
        with open(os.path.join(self.d, "metamorpheus_identifications.json")) as fh:
            rec = json.load(fh)
        t = rec["totals"]
        self.assertEqual((t["psms"], t["decoy_dropped"], t["psm_q_dropped"],
                          t["peptide_q_dropped"], t["psms_kept"], t["rows"], t["peptides"]),
                         (9, 1, 2, 1, 5, 4, 3))
        self.assertEqual(rec["per_run"]["A"]["precursors"], 3)
        self.assertIn("1% FDR", rec["methods_sentence"])

    def test_no_peptide_table_writes_nothing(self):
        import run_search
        os.remove(os.path.join(self.d, "run", "Task1SearchTask", "AllPeptides.psmtsv"))
        with self.assertRaises(SystemExit) as e:
            run_search.adapt_metamorpheus_ids(self.d)
        self.assertIn("AllPeptides.psmtsv", str(e.exception))
        self.assertFalse(os.path.exists(os.path.join(self.d, "identifications.tsv")))


class JobRuns(Harness):
    """The generated job, run with a stub apptainer and an image that is not the pinned one (its
    sha256 is then not checked): it ends with identifications.tsv, and writes nothing beside the
    database."""

    def test_the_job_ends_with_the_identifications(self):
        sif = self.sif("metamorpheus_test.sif")
        p = self.generate(sif=sif, tools=_write(os.path.join(self.d, "tools.json"), json.dumps(
            {"metamorpheus": sif, "versions": {"metamorpheus": ""}})))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertNotIn("sha256sum", self.script())
        fix = os.path.join(self.d, "fixture")
        write_table(os.path.join(fix, "AllPSMs.psmtsv"), PSMS)
        write_table(os.path.join(fix, "AllPeptides.psmtsv"), PEPTIDES)
        stub = os.path.join(self.d, "stubbin")
        _write(os.path.join(stub, "apptainer"), STUB_APPTAINER, 0o755)
        log = os.path.join(self.d, "stub")
        env = job_env(self.d, PATH=os.pathsep.join([stub, os.path.dirname(sys.executable),
                                                    "/usr/bin", "/bin"]),
                      MM_FIXTURE=fix, MM_STUB_LOG=log, TMPDIR=os.path.join(self.d, "nodetmp"))
        os.makedirs(env["TMPDIR"])
        r = subprocess.run(["bash", self.job], cwd=self.out, capture_output=True, text=True,
                           env=env, timeout=240)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        with open(os.path.join(self.out, "identifications.tsv")) as fh:
            self.assertEqual(len(fh.read().splitlines()), 5)
        with open(os.path.join(self.out, "metamorpheus_identifications.json")) as fh:
            self.assertEqual(json.load(fh)["metamorpheus_reported"]["All"]["psms"], 6)
        with open(log + ".mods") as fh:
            self.assertIn("CustomPyroGluPeptideNterm.txt", fh.read())
        with open(log + ".db") as fh:
            self.assertTrue(fh.read().startswith(env["TMPDIR"]))
        self.assertEqual(os.listdir(self.dbdir), ["search.fasta"])   # no DatabaseIndex here
        self.assertEqual(os.listdir(env["TMPDIR"]), [])              # the copies are removed


if __name__ == "__main__":
    unittest.main()
