#!/usr/bin/env python3
"""
`run_search.py --engine sage --sbatch` converted .raw -> mzML BEFORE writing the job: run_sage()
called ensure_mzml() first, so msconvert ran on the HIVE login node for every file (gabrig
2026-09-29, 4 Fusion Lumos .raw, ~2.4 GB, login2). It failed at once only because bioconda's
Linux msconvert has no vendor readers; had it worked it would have broken golden rule 3.

Now the conversion is planned at generation and RUN in the job, on the compute node, with
ThermoRawFileParser for .raw (the parser step 2 uses, found the same way) -- nothing heavy runs
before submission. Radiant's single-job --sbatch route had the same pattern, plus the DIA-NN
predicted library built inline; both move into its job too.

The converters are stand-ins on PATH that LOG every call: generating must log none; running the
generated job (tests/job_env.py) must log one per file and leave complete mzML where Sage reads
them, then run the LFQ check (test_sage_lfq_check.py's fixture plays Sage).
"""
import json
import os
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
RUN_SEARCH = os.path.join(SCRIPTS, "run_search.py")
sys.path.insert(0, HERE)

from job_env import job_env  # noqa: E402
import test_sage_lfq_check as lfqfx  # noqa: E402  (Sage output fixture)

FAKE_TRFP = r'''#!{py}
import json, os, sys
with open(os.environ["FAKE_CONV_LOG"], "a") as fh:
    fh.write(json.dumps({{"tool": "trfp", "argv": sys.argv[1:],
                         "DOTNET_ROOT": os.environ.get("DOTNET_ROOT")}}) + "\n")
if os.environ.get("FAKE_TRFP_EXIT"):
    print("ThermoRawFileParser: pretend failure", file=sys.stderr)
    sys.exit(int(os.environ["FAKE_TRFP_EXIT"]))
out = next((a[3:] for a in sys.argv[1:] if a.startswith("-b=")), None)
if out and not os.environ.get("FAKE_TRFP_WRITE_NOTHING"):
    with open(out, "w") as fh:
        fh.write('<?xml version="1.0" encoding="utf-8"?>\n<indexedmzML><mzML/>\n')
        if not os.environ.get("FAKE_TRFP_TRUNCATE"):
            fh.write("<indexOffset>0</indexOffset>\n</indexedmzML>\n")
'''

FAKE_MSCONVERT = r'''#!{py}
import json, os, sys
with open(os.environ["FAKE_CONV_LOG"], "a") as fh:
    fh.write(json.dumps({{"tool": "msconvert", "argv": sys.argv[1:]}}) + "\n")
'''

# Plays Sage: logs its argv, then writes Sage 0.14.7-shaped output at the offsets it is told.
FAKE_SAGE = r'''#!{py}
import json, os, sys
sys.path.insert(0, {here!r})
import test_sage_lfq_check as fx
args = sys.argv[1:]
out = args[args.index("-o") + 1]
files = [a for a in args if a.endswith(".mzML")]
with open(os.environ["FAKE_CONV_LOG"], "a") as fh:
    fh.write(json.dumps({{"tool": "sage", "argv": args,
                         "missing": [f for f in files if not os.path.isfile(f)]}}) + "\n")
off = float(os.environ.get("FAKE_SAGE_OFFSET", "7.2"))
fx.write_sage_outputs(out, {{os.path.basename(f): off for f in files}}, log=False)
print("[INFO sage] discovered %s target MS1 peaks at 5%% FDR"
      % os.environ.get("FAKE_SAGE_PEAKS", "0"), file=sys.stderr)
'''


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, os.stat(path).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


@unittest.skipUnless(os.name == "posix", "the converter stand-ins are POSIX shims")
class _Harness(unittest.TestCase):
    def setUp(self):
        self.d = os.path.realpath(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.d, True)
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        py = sys.executable
        _exe(os.path.join(self.bin, "ThermoRawFileParser"), FAKE_TRFP.format(py=py))
        _exe(os.path.join(self.bin, "msconvert"), FAKE_MSCONVERT.format(py=py))
        self.sage = os.path.join(self.bin, "fake_sage")
        _exe(self.sage, FAKE_SAGE.format(py=py, here=HERE))
        self.log = os.path.join(self.d, "calls.jsonl")
        raw_dir = os.path.join(self.d, "raw")
        os.makedirs(raw_dir)
        self.raws = []
        for n in ("HeLa_1.raw", "HeLa_2.raw"):
            p = os.path.join(raw_dir, n)
            with open(p, "wb") as fh:
                fh.write(b"\x01\xa1F\x00i\x00n\x00n\x00i\x00g\x00a\x00n\x00")   # not read
            self.raws.append(p)
        self.fasta = os.path.join(self.d, "search.fasta")
        with open(self.fasta, "w") as fh:
            fh.write(">sp|P1|X_HUMAN\nPEPTIDEK\n")
        self.cfg = os.path.join(self.d, "sage_config.json")
        with open(self.cfg, "w") as fh:
            json.dump({"database": {"fasta": "REPLACED_AT_RUNTIME.fasta"},
                       "precursor_tol": {"ppm": [-10.0, 10.0]},
                       "fragment_tol": {"ppm": [-10.0, 10.0]}, "quant": {"lfq": True}}, fh)
        self.bundle = os.path.join(self.d, "workflow.manifest.json")
        with open(self.bundle, "w") as fh:
            json.dump({"acquisition": "DDA", "engine": {"name": "sage", "version": "0.14.7"}},
                      fh)
        self.tools = os.path.join(self.d, "tools.json")
        with open(self.tools, "w") as fh:
            json.dump({"sage": self.sage, "versions": {"sage": "v0.14.7"}}, fh)
        self.out = os.path.join(self.d, "search")

    def env(self, **extra):
        # No sbatch on PATH (the login-node guard stays out of it), the stand-ins first, and no
        # real parser reachable: not by $THERMORAWFILEPARSER, a shared copy or the pipeline env.
        env = job_env(self.d, PATH=self.bin + os.pathsep + "/usr/bin:/bin",
                      FAKE_CONV_LOG=self.log, THERMORAWFILEPARSER_SHARED="",
                      PROTEOMICS_PIPELINE_HOME=os.path.join(self.d, "home"), **extra)
        env.pop("THERMORAWFILEPARSER", None)
        return env

    def calls(self):
        if not os.path.exists(self.log):
            return []
        with open(self.log) as fh:
            return [json.loads(ln) for ln in fh if ln.strip()]

    def run_search(self, *extra, files=None, env=None, tools=None, engine="sage"):
        return subprocess.run(
            [sys.executable, RUN_SEARCH, "--tools", tools or self.tools, "--bundle", self.bundle,
             "--params", self.cfg, "--fasta", self.fasta, "--out", self.out,
             "--files", *(files or self.raws), "--engine", engine, "--threads", "4", *extra],
            capture_output=True, text=True, timeout=180, env=env or self.env(), cwd=self.d)


class SageSbatchConvertsInTheJob(_Harness):
    def generate(self):
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        return job, open(job).read()

    def test_nothing_is_converted_before_submission(self):
        job, text = self.generate()
        self.assertEqual(self.calls(), [], "a converter ran at generation, on the login node")
        mzdir = os.path.join(self.out, "mzml")
        self.assertFalse(os.path.isdir(mzdir) and os.listdir(mzdir),
                         "mzML appeared before the job ran")

    def test_the_job_converts_with_thermorawfileparser_before_sage(self):
        job, text = self.generate()
        trfp = os.path.join(self.bin, "ThermoRawFileParser")
        for raw in self.raws:
            self.assertIn(f"-i={raw}", text)
        self.assertIn(trfp, text)
        self.assertIn("-f=2", text)
        self.assertNotIn("msconvert", text)
        conv = text.index(trfp)
        sage = text.index(self.sage)
        self.assertLess(conv, sage, "the job must convert before it searches")
        self.assertIn("| tee ", text[sage:])                  # Sage's log kept for the check
        self.assertIn("sage_lfq_check.py", text[sage:])
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            prov = json.load(fh)
        conv_rec = prov["result"]["mzml_conversion"]
        self.assertEqual(conv_rec["where"], "in the job, before Sage")
        self.assertEqual([f["input"] for f in conv_rec["files"]], self.raws)

    def test_running_the_job_converts_searches_and_warns(self):
        job, _ = self.generate()
        r = subprocess.run(["bash", job], capture_output=True, text=True, timeout=180,
                           env=self.env(FAKE_SAGE_OFFSET="7.4"), cwd=self.d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        calls = self.calls()
        self.assertEqual([c["tool"] for c in calls], ["trfp", "trfp", "sage"], calls)
        self.assertEqual(calls[-1]["missing"], [], "Sage was handed mzML that did not exist")
        for raw in self.raws:
            mz = os.path.join(self.out, "mzml",
                              os.path.splitext(os.path.basename(raw))[0] + ".mzML")
            self.assertTrue(os.path.isfile(mz), mz)
            self.assertIn(mz, calls[-1]["argv"])
        self.assertIn("discovered 0 target MS1 peaks", open(os.path.join(self.out,
                                                                         "sage.log")).read())
        self.assertIn("[sage_lfq_check] WARNING:", r.stderr)
        with open(os.path.join(self.out, "sage_lfq_check.json")) as fh:
            rec = json.load(fh)
        self.assertEqual(rec["status"], "warn")
        self.assertEqual(rec["target_ms1_peaks_5pct_fdr"], 0)
        self.assertIn("sage.log", rec["target_ms1_peaks_source"])
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            self.assertEqual(json.load(fh)["sage_lfq_check"]["status"], "warn")

    def test_a_failing_converter_says_so_before_the_job_ends(self):
        """Under set -euo pipefail a converter that exits non-zero used to end the job before
        its FAILED line was printed, leaving the watcher nothing to classify."""
        job, _ = self.generate()
        r = subprocess.run(["bash", job], capture_output=True, text=True, timeout=180,
                           env=self.env(FAKE_TRFP_EXIT="3"), cwd=self.d)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("FAILED: no mzML from ThermoRawFileParser", r.stderr)
        self.assertNotIn("sage", [c["tool"] for c in self.calls()])

    def test_a_truncated_mzml_stops_the_job_before_sage(self):
        """Exit 0 is not proof: a parser that stops part way leaves an mzML without its
        closing tag, and the job must not search it."""
        job, _ = self.generate()
        r = subprocess.run(["bash", job], capture_output=True, text=True, timeout=180,
                           env=self.env(FAKE_TRFP_TRUNCATE="1"), cwd=self.d)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("FAILED: no mzML from ThermoRawFileParser", r.stderr)
        self.assertNotIn("sage", [c["tool"] for c in self.calls()])
        self.assertFalse(os.path.exists(os.path.join(self.out, "mzml", "HeLa_1.mzML")))


class SageInlineStillConverts(_Harness):
    def test_inline_converts_here_then_checks(self):
        p = self.run_search(env=self.env(FAKE_SAGE_OFFSET="0.6", FAKE_SAGE_PEAKS="350"))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual([c["tool"] for c in self.calls()], ["trfp", "trfp", "sage"])
        self.assertIn("[sage_lfq_check] INFO:", p.stderr)
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertEqual(prov["sage_lfq_check"]["status"], "ok")
        self.assertEqual(prov["result"]["mzml_conversion"]["where"], "inline, before Sage")
        self.assertTrue(os.path.isfile(os.path.join(self.out, "report.parquet")))


class ConversionRefusals(_Harness):
    def test_no_converter_is_refused_before_anything_is_written(self):
        for n in ("ThermoRawFileParser", "msconvert"):
            os.remove(os.path.join(self.bin, n))
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("ThermoRawFileParser", p.stderr)
        self.assertIn("Nothing was converted, written or submitted", p.stderr)
        self.assertFalse(os.path.exists(job))

    def test_msconvert_is_the_fallback_for_raw_and_warns_on_linux(self):
        os.remove(os.path.join(self.bin, "ThermoRawFileParser"))
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(self.calls(), [])
        text = open(job).read()
        self.assertIn(os.path.join(self.bin, "msconvert"), text)
        # centroided MS1+MS2, and the same indexed-mzML completeness check as the parser
        self.assertIn("'peakPicking vendor msLevel=1-'", text)
        self.assertEqual(text.count("*'</indexedmzML>'*"), len(self.raws))   # one per file
        if sys.platform.startswith("linux"):
            self.assertIn("no vendor", p.stderr)

    def test_bruker_d_keeps_msconvert_in_the_job(self):
        d_dir = os.path.join(self.d, "raw", "S1.d")
        os.makedirs(d_dir)
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job, files=[d_dir])
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(self.calls(), [])
        text = open(job).read()
        self.assertIn(os.path.join(self.bin, "msconvert"), text)
        self.assertIn(os.path.join(self.out, "mzml", "S1.mzML"), text)

    def test_two_inputs_with_one_name_are_refused(self):
        other = os.path.join(self.d, "raw2")
        os.makedirs(other)
        twin = os.path.join(other, "HeLa_1.raw")
        shutil.copy(self.raws[0], twin)
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job, files=[self.raws[0], twin])
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("share an mzML name", p.stderr)
        self.assertFalse(os.path.exists(job))

    def test_mzml_input_needs_no_converter(self):
        mz = os.path.join(self.d, "raw", "ready.mzML")
        open(mz, "w").close()
        for n in ("ThermoRawFileParser", "msconvert"):
            os.remove(os.path.join(self.bin, n))
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job, files=[mz])
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertNotIn("converting", open(job).read())


class RadiantSbatchBuildsInTheJob(_Harness):
    """Radiant's single-job --sbatch route: .raw conversion AND the DIA-NN predicted library
    both used to run here, before the job was written."""

    def test_conversion_and_library_are_in_the_job(self):
        diann = os.path.join(self.bin, "diann-linux")
        _exe(diann, "#!/bin/sh\necho DIANN-RAN >> \"$FAKE_CONV_LOG.diann\"\n")
        tools = os.path.join(self.d, "tools_radiant.json")
        with open(tools, "w") as fh:
            json.dump({"radiant": "apptainer exec", "radiant_runtime": "apptainer",
                       "radiant_image": "/img/radiant.sif", "diann": diann,
                       "versions": {"radiant": "2.3.3"}}, fh)
        cfg = os.path.join(self.d, "default.radiantConfig")
        open(cfg, "w").close()
        job = os.path.join(self.d, "radiant_job.sh")
        p = subprocess.run(
            [sys.executable, RUN_SEARCH, "--tools", tools, "--bundle", self.bundle,
             "--params", cfg, "--fasta", self.fasta, "--out", self.out,
             "--files", self.raws[0], "--engine", "radiant", "--threads", "4",
             "--sbatch", job],
            capture_output=True, text=True, timeout=180, env=self.env(), cwd=self.d)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual(self.calls(), [])
        self.assertFalse(os.path.exists(self.log + ".diann"), "DIA-NN ran at generation")
        text = open(job).read()
        conv = text.index(os.path.join(self.bin, "ThermoRawFileParser"))
        lib = text.index("make_radiant_library.py")
        search = text.index("radiant_fulcrum")
        self.assertLess(conv, search)
        self.assertLess(lib, search)
        # the library the job builds is the one the container is given
        self.assertIn(f"--bind {os.path.join(self.out, 'radiant_lib')}:/mnt/in0", text)
        self.assertIn("--library /mnt/in0/radiant_library.tsv", text[search:])
        self.assertIn(os.path.join(self.out, "mzml") + ":/mnt/in2", text)


if __name__ == "__main__":
    unittest.main()
