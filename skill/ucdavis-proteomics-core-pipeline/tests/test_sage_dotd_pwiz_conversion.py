#!/usr/bin/env python3
"""
Sage quantification converts a timsTOF .d in the job with the pinned vendor-enabled ProteoWizard
container (I45; T21 of the 2.14.2 plan).

Sage 0.14.7 has no MS1 quantification on a .d, so a quantifying search (quant.lfq true, the
default) needs mzML. mzml_plan() sent every .d to msconvert from PATH, and bioconda's Linux
msconvert has no Bruker reader: on HIVE the job could only fail. The Core already pins a
vendor-enabled ProteoWizard image for Skyline (resolve_defaults.SKYLINE_PIN), with msconvert in it
under Wine, so the .d is converted with that, in the job:

  * the job checks the pinned container's sha256 before the first conversion (a sif other than
    the Core's is used as given, as skyline_peaks.verify_container treats one);
  * every path stays one word, a service folder with a space or a quote included (Review Focus 2);
  * with no container there, nothing is written, and the refusal names both ways out: an mzML
    from a vendor-enabled msconvert, or a search for identifications only (--ids-only);
  * a folder apptainer cannot bind (':' or ',') is refused before anything is written;
  * search_provenance.json names the container and its version (only the pinned one's is known).

msconvert under Wine is given plain Linux paths, as ProteoWizard's own container README runs it
(`wine msconvert /data/file.raw`); the HIVE run decides whether Windows paths (Z:\\...) are needed.

apptainer, Sage and ThermoRawFileParser are stand-ins; the sif is a small file. Hermetic.
"""
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
WATCH = os.path.join(SCRIPTS, "watch_run.sh")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import run_search  # noqa: E402
from resolve_defaults import SKYLINE_PIN  # noqa: E402
from test_sage_sbatch_conversion import _Harness, _exe  # noqa: E402
# self.env() is job_env() (see _Harness.env); imported so the guard can see that.
from job_env import job_env  # noqa: E402,F401
import test_sage_lfq_check as lfqfx  # noqa: E402

# The Core's Windows msconvert settings for timsTOF DDA-PASEF (3.0.26058, from the <dataProcessing>
# of the mzML it made for a Core submission; H2-F1), written as MSConvertGUI writes them. Literal here, so
# a change to run_search's constant has to change this test too.
SCAN_SUMMING = "scanSumming precursorTol=0.05 scanTimeTol=10 ionMobilityTol=0.1 sumMs1=0"
MSCONVERT_TIMS = ["--mzML", "--zlib", "--combineIonMobilitySpectra",
                  "--filter", "peakPicking vendor msLevel=1-", "--filter", SCAN_SUMMING]

# Plays `apptainer exec ... <sif> wine msconvert <d> ... -o <dir>`: logs its argv, then writes the
# indexed mzML msconvert would (<dir>/<d name>.mzML), complete unless told otherwise.
FAKE_APPTAINER = r'''#!{py}
import json, os, sys
args = sys.argv[1:]
with open(os.environ["FAKE_CONV_LOG"], "a") as fh:
    fh.write(json.dumps({{"tool": "apptainer", "argv": args}}) + "\n")
i = args.index("msconvert")
d, outdir = args[i + 1], args[args.index("-o") + 1]
if not os.path.isdir(d):
    print("msconvert: no such run " + d, file=sys.stderr)
    sys.exit(2)
with open(os.path.join(outdir, os.path.splitext(os.path.basename(d))[0] + ".mzML"), "w") as fh:
    fh.write('<?xml version="1.0" encoding="utf-8"?>\n<indexedmzML><mzML/>\n')
    if not os.environ.get("FAKE_PWIZ_TRUNCATE"):
        fh.write("<indexOffset>0</indexOffset>\n</indexedmzML>\n")
'''


class _PwizHarness(_Harness):
    def setUp(self):
        super().setUp()
        _exe(os.path.join(self.bin, "apptainer"), FAKE_APPTAINER.format(py=sys.executable))
        self.tims = os.path.join(self.d, "tims")
        self.dots = []
        for n in ("runA_S1-A1_1_1001.d", "runB_S1-B1_1_1002.d"):
            p = os.path.join(self.tims, n)
            os.makedirs(p)
            self.dots.append(p)
        # a container that is not the Core's (any name but the pinned one's)
        self.sif = os.path.join(self.d, "containers", "my-pwiz.sif")
        os.makedirs(os.path.dirname(self.sif))
        with open(self.sif, "wb") as fh:
            fh.write(b"not really a sif")
        # one with the pinned name, whose bytes do not match the pinned sha256
        self.pinned = os.path.join(self.d, "apptainers", os.path.basename(SKYLINE_PIN["sif"]))
        os.makedirs(os.path.dirname(self.pinned))
        with open(self.pinned, "wb") as fh:
            fh.write(b"a different image under the pinned name")

    def env(self, **extra):
        # the real sha256sum (coreutils on HIVE, /sbin on macOS) beside the stand-ins
        env = super().env(**extra)
        env["PATH"] += os.pathsep + "/sbin"
        return env

    def generate(self, files=None, sif=None, expect=0):
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job, "--pwiz-sif", sif or self.sif,
                            files=files or self.dots)
        if expect is None:
            return job, p
        self.assertEqual(p.returncode, expect, p.stdout + p.stderr)
        with open(job) as fh:
            return job, fh.read()

    def provenance(self):
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            return json.load(fh)

    def apptainer_lines(self, text):
        return [ln for ln in text.splitlines() if ln.startswith("apptainer ")]


class TheJobConvertsWithTheContainer(_PwizHarness):
    def test_dotd_quant_converts_with_the_container(self):
        _job, text = self.generate(sif=self.pinned)
        self.assertIn("apptainer exec", text)
        self.assertIn("wine msconvert", text)
        self.assertNotIn(os.path.join(self.bin, "msconvert"), text, "bioconda's msconvert")
        sha = text.index("sha256sum -c")
        self.assertLess(sha, text.index("apptainer exec"), "the container is checked first")
        line = text[text.rindex("\n", 0, sha) + 1:text.index("\n", sha)]
        self.assertIn(shlex.quote(f"{SKYLINE_PIN['sif_sha256']}  {self.pinned}"), line)
        self.assertLess(text.index("apptainer exec"), text.index(self.sage))
        self.assertIn("sage_lfq_check.py", text[text.index(self.sage):])
        self.assertEqual(self.calls(), [], "nothing runs at generation")

    def test_the_container_argv(self):
        """apptainer exec --cleanenv --env WINEDEBUG=-all, the .d's folder read-only and the
        partial mzML folder read-write, then wine msconvert as the old step ran it."""
        _job, text = self.generate()
        lines = self.apptainer_lines(text)
        self.assertEqual(len(lines), 2, text)
        argv = shlex.split(lines[0].split(" || ")[0])
        partial = os.path.join(self.out, "mzml", ".partial")
        self.assertEqual(argv, [
            "apptainer", "exec", "--cleanenv", "--env", "WINEDEBUG=-all",
            "--bind", f"{self.tims}:{self.tims}:ro", "--bind", f"{partial}:{partial}",
            self.sif, "wine", "msconvert", self.dots[0], *MSCONVERT_TIMS, "-o", partial])

    def test_the_job_combines_ion_mobility_and_sums_scans(self):
        """H2-F1: without --combineIonMobilitySpectra and scanSumming, msconvert writes one
        spectrum per mobility scan (HIVE, one timsTOF HT run: 11,139,896 spectra, past 9.9 GB,
        against 82,905 from the Core's Windows msconvert). Vendor peakPicking stays the FIRST
        filter."""
        _job, text = self.generate()
        line = self.apptainer_lines(text)[0]
        self.assertIn(" --combineIonMobilitySpectra ", line)
        self.assertIn("--filter " + shlex.quote(SCAN_SUMMING), line)
        argv = shlex.split(line.split(" || ")[0])
        filters = [argv[i + 1] for i, x in enumerate(argv) if x == "--filter"]
        self.assertEqual(filters, ["peakPicking vendor msLevel=1-", SCAN_SUMMING])

    def test_path_with_a_space_is_one_word(self):
        """Review Focus 2: a service folder with a space or a quote stays one word."""
        svc = os.path.join(self.d, "UC Davis", "PROT_0000 'Lab A'")
        d = os.path.join(svc, "run 1.d")
        os.makedirs(d)
        _job, text = self.generate(files=[d])
        self.assertIn(shlex.quote(d), text)
        argv = shlex.split(self.apptainer_lines(text)[0].split(" || ")[0])
        self.assertIn(d, argv)
        self.assertIn(f"{svc}:{svc}:ro", argv)

    def test_non_core_sif_skips_the_sha_line(self):
        _job, text = self.generate()
        self.assertNotIn("sha256sum", text)
        self.assertIn(shlex.quote(self.sif), text)

    def test_mixed_quant_raw_and_dotd(self):
        """A quantifying search with both: the .raw by ThermoRawFileParser, the .d by the
        container, both before Sage."""
        _job, text = self.generate(files=[self.raws[0], self.dots[0]])
        self.assertIn("-i=" + self.raws[0], text)
        self.assertEqual(len(self.apptainer_lines(text)), 1)
        sage = text.index(self.sage)
        self.assertLess(text.index("-i=" + self.raws[0]), sage)
        self.assertLess(text.index("apptainer exec"), sage)


class Refusals(_PwizHarness):
    def test_missing_sif_is_refused_with_both_ways_out(self):
        job, p = self.generate(sif=os.path.join(self.d, "nowhere", "pwiz.sif"), expect=None)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("--ids-only", p.stderr)
        self.assertIn("vendor-enabled msconvert", p.stderr)
        self.assertIn(os.path.join(self.d, "nowhere", "pwiz.sif"), p.stderr)
        self.assertIn("Nothing was converted, written or submitted", p.stderr)
        self.assertFalse(os.path.exists(job))
        self.assertFalse(os.path.exists(os.path.join(self.out, "search_provenance.json")))
        # L7R-F1: the command it hands over is the one the job would run (H2-F1's options)
        self.assertIn(run_search.pwiz_msconvert_options(windows=True), p.stderr)
        self.assertIn("--combineIonMobilitySpectra", p.stderr)

    def test_ids_only_needs_no_container(self):
        cfg = os.path.join(self.d, "ids.json")
        with open(cfg, "w") as fh:
            json.dump({"quant": {"lfq": False}}, fh)
        self.cfg = cfg
        _job, text = self.generate(sif=os.path.join(self.d, "nowhere", "pwiz.sif"))
        self.assertNotIn("apptainer", text)

    def test_a_folder_apptainer_cannot_bind_is_refused(self):
        d = os.path.join(self.d, "plate:3", "r1.d")
        os.makedirs(d)
        job, p = self.generate(files=[d], expect=None)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("apptainer", p.stderr)
        self.assertIn("plate:3", p.stderr)
        self.assertFalse(os.path.exists(job))


class TheRecord(_PwizHarness):
    def test_provenance_names_the_container_and_version(self):
        self.generate(sif=self.pinned)
        prov = self.provenance()
        conv = prov["result"]["mzml_conversion"]
        self.assertEqual([f["input"] for f in conv["files"]], self.dots)
        self.assertEqual(conv["tools"], [f"ProteoWizard msconvert (pinned Skyline container "
                                         f"{SKYLINE_PIN['version']})"])
        c = conv["containers"]
        self.assertEqual(len(c), 1)
        self.assertEqual(c[0]["sif"], self.pinned)
        self.assertEqual(c[0]["version"], SKYLINE_PIN["version"])
        self.assertEqual(c[0]["sif_sha256"], SKYLINE_PIN["sif_sha256"])
        self.assertEqual(prov["sage_input"]["dotd"], "converted")
        self.assertIn(SKYLINE_PIN["version"], prov["sage_input"]["why"])

    def test_the_record_names_the_combine_and_the_filter(self):
        """H2-F1: the provenance text says how the .d was converted, from the one constant."""
        self.generate(sif=self.pinned)
        prov = self.provenance()
        c = prov["result"]["mzml_conversion"]["containers"][0]
        options = " ".join(shlex.quote(x) for x in MSCONVERT_TIMS)
        self.assertEqual(c["msconvert_options"], options)
        self.assertIn("dataProcessing", c["msconvert_options_source"])
        self.assertIn(options, prov["sage_input"]["why"])
        for f in prov["result"]["mzml_conversion"]["files"]:
            self.assertIn(options, f["command"])

    def test_no_submission_number_in_any_record(self):
        """L7R-F10: search_provenance.json is delivered to every collaborator."""
        self.generate(sif=self.pinned)
        with open(os.path.join(self.out, "search_provenance.json")) as fh:
            self.assertNotRegex(fh.read(), r"PROT_?\d{4}")

    def test_the_check_is_said_as_the_jobs_not_as_done(self):
        """L7R-F14: written at generation, before anything was checked."""
        self.generate(sif=self.pinned)
        c = self.provenance()["result"]["mzml_conversion"]["containers"][0]
        self.assertIn("checked by the job before the first conversion", c["verified"])
        self.assertIn("stops on a mismatch", c["verified"])

    def test_msconvert_version_is_the_images_own(self):
        """L7R-F11: the pinned image's own msconvert, read on HIVE (`wine msconvert --help`:
        "ProteoWizard release: 3.0.26266 (578b175)"; its mzML's pwiz_Reader_Bruker agrees) --
        not Skyline's 26.1.0.266, nor 3.0.26058, the Windows build the settings came from."""
        self.assertEqual(SKYLINE_PIN["msconvert_version"], "3.0.26266")
        self.generate(sif=self.pinned)
        c = self.provenance()["result"]["mzml_conversion"]["containers"][0]
        self.assertEqual(c["msconvert_version"], "3.0.26266")
        self.assertNotIn("not recorded", c["msconvert_version_source"])

    def test_the_record_says_what_summs1_does_in_the_pinned_image(self):
        """H2-7: the pinned 3.0.26266 takes no sumMs1 and writes sumMS1="true" into its mzML,
        yet MS1 was not summed in effect (spectrum counts, IDs, MS1 peaks and LFQ all equal the
        Windows 3.0.26058 output). The argv keeps sumMs1=0; the record says the truth beside it,
        and the Methods sentence stays about what was done."""
        self.generate(sif=self.pinned)
        prov = self.provenance()
        c = prov["result"]["mzml_conversion"]["containers"][0]
        self.assertIn("sumMs1=0", c["msconvert_options"])
        note = c["msconvert_options_note"]
        for part in ("3.0.26266", "sumMS1=true", "not summed in effect"):
            self.assertIn(part, note)
        self.assertNotIn("sumMS1=true", prov["sage_input"]["methods_sentence"])

    def test_another_sifs_msconvert_version_is_not_recorded(self):
        self.generate()
        c = self.provenance()["result"]["mzml_conversion"]["containers"][0]
        self.assertIsNone(c["msconvert_version"])
        self.assertIn("not recorded", c["msconvert_version_source"])

    def test_a_rerun_keeps_a_container_that_is_not_the_pins(self):
        """L7R-F11: the LFQ gate's re-run command and reproduce.sh carry --pwiz-sif."""
        import sage_lfq_check as slc
        self.generate()
        prov = self.provenance()
        self.assertEqual(prov["pwiz_sif_arg"], ["--pwiz-sif", self.sif])
        cmd = slc.rerun_plan(self.out, os.path.join(self.out, "corrected.json"), 10)["command"]
        self.assertIn(f"--pwiz-sif {shlex.quote(self.sif)}", cmd)
        de = os.path.join(self.d, "de")
        os.makedirs(de)
        repro = os.path.join(self.d, "repro")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                            "--outdir", repro, "--de-dir", de,
                            "--report", os.path.join(self.out, "report.parquet"),
                            "--engine", "sage", "--acquisition", "DDA",
                            "--instrument", "timsTOF HT", "--params", self.cfg],
                           capture_output=True, text=True, timeout=300)
        self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        with open(os.path.join(repro, "reproduce.sh")) as fh:
            sh = fh.read()
        line = sh[sh.index('python3 "$SKILL/scripts/run_search.py"'):]
        # the search command ends at the first blank line, whatever follows it: a DE step 6, or
        # none when no DE ran (provenance.py writes no step 6 then -- 2.14.2 C2 fixer D, F2)
        end = line.find("\n\n")
        line = line if end < 0 else line[:end]
        self.assertNotIn("run_de.R", line)
        self.assertIn(f"--pwiz-sif {shlex.quote(self.sif)}", line)

    def test_a_rerun_keeps_the_searchs_memory(self):
        """L7R-mem: --mem <GB> is recorded as queue.mem "<GB>G"; the gate's re-run command asks
        for the same, and for nothing when the search used the default."""
        import sage_lfq_check as slc
        job = os.path.join(self.d, "sage_job.sh")
        p = self.run_search("--sbatch", job, "--pwiz-sif", self.sif, "--mem", "96",
                            files=self.dots)
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(self.provenance()["queue"]["mem"], "96G")
        cmd = slc.rerun_plan(self.out, os.path.join(self.out, "corrected.json"), 10)["command"]
        self.assertIn(" --mem 96 ", cmd + " ")
        self.generate()                                  # the default memory
        self.assertNotIn("mem", self.provenance()["queue"])
        cmd = slc.rerun_plan(self.out, os.path.join(self.out, "corrected.json"), 10)["command"]
        self.assertNotIn("--mem", cmd)

    def test_the_pins_own_path_needs_no_flag(self):
        self.assertEqual(run_search.pwiz_rerun_args(SKYLINE_PIN["sif"]), [])
        self.assertEqual(run_search.pwiz_rerun_args(self.sif), ["--pwiz-sif", self.sif])

    def test_the_methods_describe_the_conversion(self):
        """L7R-F3: the .d -> mzML step, its options (the one constant) and the container, in
        the record's own sentence, quoted by make_methods."""
        import make_methods
        self.generate(sif=self.pinned)
        prov = self.provenance()
        sent = prov["sage_input"]["methods_sentence"]
        self.assertIn("converted to mzML with ProteoWizard msconvert", sent)
        self.assertIn(run_search.pwiz_msconvert_options(), sent)
        self.assertIn(SKYLINE_PIN["tag"], sent)
        self.assertIn("ProteoWizard msconvert 3.0.26266", sent)
        self.assertNotIn("not recorded", sent)
        rec = make_methods.search_record(params=self.cfg, search_prov=os.path.join(
            self.out, "search_provenance.json"))
        self.assertIn(sent, make_methods.search_paragraph(rec))

    def test_another_sif_has_no_version_it_cannot_know(self):
        self.generate()
        c = self.provenance()["result"]["mzml_conversion"]["containers"][0]
        self.assertEqual(c["sif"], self.sif)
        self.assertIsNone(c["version"])
        self.assertIn("not the Core's pinned image", c["verified"])
        self.assertNotIn(SKYLINE_PIN["version"],
                         self.provenance()["result"]["mzml_conversion"]["tools"][0])


@unittest.skipUnless(os.path.exists("/sbin/sha256sum") or shutil.which("sha256sum"),
                     "the job's container check runs sha256sum")
class RunningTheJob(_PwizHarness):
    @unittest.skipUnless(lfqfx.HAVE_ARROW, "pyarrow is needed to write the Sage fixtures")
    def test_the_job_converts_in_the_container_then_searches(self):
        job, _ = self.generate()
        r = subprocess.run(["bash", job], capture_output=True, text=True, timeout=180,
                           env=self.env(), cwd=self.d)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        calls = self.calls()
        self.assertEqual([c["tool"] for c in calls], ["apptainer", "apptainer", "sage"], calls)
        self.assertEqual(calls[-1]["missing"], [])
        for d in self.dots:
            mz = os.path.join(self.out, "mzml", os.path.basename(d)[:-2] + ".mzML")
            self.assertTrue(os.path.isfile(mz), mz)
            self.assertIn(mz, calls[-1]["argv"])

    def test_a_pinned_name_that_does_not_match_stops_the_job_first(self):
        job, _ = self.generate(sif=self.pinned)
        r = subprocess.run(["bash", job], capture_output=True, text=True, timeout=180,
                           env=self.env(), cwd=self.d)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn(f"FAILED: the ProteoWizard container {self.pinned} does not match its "
                      f"pinned sha256", r.stderr)
        self.assertEqual(self.calls(), [], "nothing was converted or searched")

    def test_a_truncated_mzml_stops_the_job_before_sage(self):
        job, _ = self.generate()
        r = subprocess.run(["bash", job], capture_output=True, text=True, timeout=180,
                           env=self.env(FAKE_PWIZ_TRUNCATE="1"), cwd=self.d)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("FAILED: no mzML from ProteoWizard msconvert", r.stderr)
        self.assertNotIn("sage", [c["tool"] for c in self.calls()])


class Inline(_PwizHarness):
    def test_inline_checks_the_pinned_sha_before_converting(self):
        p = self.run_search("--pwiz-sif", self.pinned, files=self.dots)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("does not match its pinned sha256", p.stderr)
        self.assertEqual(self.calls(), [])

    @unittest.skipUnless(lfqfx.HAVE_ARROW, "pyarrow is needed to write the Sage fixtures")
    def test_inline_converts_with_the_same_step(self):
        p = self.run_search("--pwiz-sif", self.sif, files=self.dots,
                            env=self.env(FAKE_SAGE_OFFSET="0.6", FAKE_SAGE_PEAKS="350"))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual([c["tool"] for c in self.calls()], ["apptainer", "apptainer", "sage"])
        self.assertEqual(self.provenance()["result"]["mzml_conversion"]["where"],
                         "inline, before Sage")


class WatcherNamesTheContainer(unittest.TestCase):
    def watch(self, log_text):
        with tempfile.TemporaryDirectory() as d:
            log = os.path.join(d, "sage_search_1.log")
            with open(log, "w") as fh:
                fh.write(log_text)
            for name, body in (("sacct", 'echo "  FAILED"'), ("squeue", "")):
                p = os.path.join(d, name)
                with open(p, "w") as fh:
                    fh.write(f"#!/bin/bash\n{body}\n")
                os.chmod(p, 0o755)
            env = dict(os.environ, PATH=d + os.pathsep + os.environ.get("PATH", ""))
            r = subprocess.run(["bash", WATCH, "--slurm", "999", "--log", log],
                               capture_output=True, text=True, env=env, timeout=60)
            self.assertEqual(r.returncode, 0, r.stderr)
            return json.loads(r.stdout)

    def test_a_failed_container_conversion_is_classified_with_its_fix(self):
        out = self.watch("FAILED: no mzML from ProteoWizard msconvert (pinned Skyline container "
                         "26.1.0.266) for /q/r1.d (expected a complete /q/x.mzML). Nothing was "
                         "searched.\n")
        self.assertEqual(out["error_class"], "sage_no_mzml")
        self.assertIn("ProteoWizard container", out["fix"])
        self.assertIn("--ids-only", out["fix"])

    def test_a_container_that_does_not_match_its_pin_is_classified(self):
        out = self.watch("FAILED: the ProteoWizard container /q/pwiz.sif does not match its "
                         "pinned sha256\n")
        self.assertEqual(out["error_class"], "sage_no_mzml")
        self.assertIn("--pwiz-sif", out["fix"])


if __name__ == "__main__":
    unittest.main()
