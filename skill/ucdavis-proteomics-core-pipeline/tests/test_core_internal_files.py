#!/usr/bin/env python3
"""C2 security M1 (2.14.2): Core-internal files never leave the Core.

`deliver --mode interim` walks whole search folders. Its JOB_FILES rule matched `slurm-*.out`, a
name the skill never writes: the skill's own scheduler logs are <job>_<jobid>.log
(run_search.py, diann_parallel.py, radiant_parallel.py, on_compute.py ...), so they were
delivered -- node names, paths with the HIVE login, the job-end hook's FRAN, Slack and run-log
lines. So were fran_deposit.json (the login that staged the search into the Core's corpus) and
node_faults.json, and sage_lfq_check.json / search_provenance.json named who accepted an LFQ
window by login, as did AUDIT.md. The session zip, which travels to collaborators too, kept the
same files outside logs/.

Now one rule, staff.is_core_internal (by name, wherever the file is), keeps them out of every
delivery mode and the zip, and staff.client_copy hands a client the role, never the login, of who
accepted (the session keeps its own copy, as commands.log keeps --by).
Stdlib only, no network.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import staff                                   # noqa: E402
import test_core_submission as tcs             # noqa: E402  (the deliver fixture, the fake API)
import test_deliver_interim as tdi             # noqa: E402  (a searched, unfinished session)
from job_env import job_env                    # noqa: E402

LOGIN = "someanalyst"                          # a made-up login


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def accepted_record(by=LOGIN):
    return {"status": "warn", "accepted": {"reason": "the user accepts it", "by": by,
                                           "by_source": "the login that ran it", "at": "2026-11-01"}}


class TheOneRule(unittest.TestCase):
    def test_the_skills_own_job_logs_and_records_are_core_internal(self):
        for name in ("diann_search_25404434.log", "fragpipe_search_1.log", "sage_search_7.log",
                     "diann_step3_123456_4.log", "diatracer_88_2.log", "deliver_999.log",
                     "s5_report_12.log", "slurm-4242.out", "slurm-4242.err",
                     "fran_deposit.json", "node_faults.json"):
            self.assertTrue(staff.is_core_internal(name), name)
            self.assertTrue(staff.is_core_internal(os.path.join("output", "search", name)), name)

    def test_deliverables_are_not(self):
        for name in ("report.log.txt", "report.parquet", "search_provenance.json",
                     "sage_lfq_check.json", "run_search.sbatch", "submit.sh", "commands.log",
                     "identifications.tsv", "AUDIT.md", "S1_01.raw", "params.resolved.cfg"):
            self.assertFalse(staff.is_core_internal(name), name)

    def test_the_names_are_the_ones_the_writers_use(self):
        import fran_deposit
        import node_fault
        self.assertTrue(staff.is_core_internal(fran_deposit.RECEIPT))
        self.assertTrue(staff.is_core_internal(node_fault.RECORD))
        # the generators' own patterns (#SBATCH --output ..._%j.log, _%A_%a.log, %x_%j.log)
        for script, pattern in (("run_search.py", '_%j.log"'), ("diann_parallel.py", "_%A_%a.log"),
                                ("radiant_parallel.py", "%x_%j.log")):
            with open(os.path.join(SCRIPTS, script), encoding="utf-8") as fh:
                self.assertIn(pattern, fh.read(), script)
        for name in ("diann_search_%j.log", "diann_step2_%A_%a.log", "radiant_%j.log"):
            filled = name.replace("%j", "25404434").replace("%A", "77").replace("%a", "3")
            self.assertTrue(staff.is_core_internal(filled), filled)

    def test_client_copy_gives_the_role_never_the_login(self):
        with tempfile.TemporaryDirectory() as d:
            lfq = write(os.path.join(d, "sage_lfq_check.json"), json.dumps(accepted_record()))
            prov = write(os.path.join(d, "search_provenance.json"),
                         json.dumps({"engine": "sage", "sage_lfq_check": accepted_record()}))
            for path, keys in ((lfq, ("accepted",)), (prov, ("sage_lfq_check", "accepted"))):
                text = staff.client_copy(path)
                self.assertIsNotNone(text, path)
                self.assertNotIn(LOGIN, text)
                block = json.loads(text)
                for k in keys:
                    block = block[k]
                self.assertIn(block["by"], (staff.ROLE, staff.USER_ROLE))
                self.assertEqual(block["reason"], "the user accepts it")
                self.assertIn(LOGIN, open(path).read())        # the session's copy keeps it
            # nothing to redact: copied as it is
            plain = write(os.path.join(d, "x", "search_provenance.json"),
                          json.dumps({"engine": "diann"}))
            self.assertIsNone(staff.client_copy(plain))
            role = write(os.path.join(d, "y", "sage_lfq_check.json"),
                         json.dumps(accepted_record(by=staff.USER_ROLE)))
            self.assertIsNone(staff.client_copy(role))
            self.assertIsNone(staff.client_copy(write(os.path.join(d, "report.log.txt"), LOGIN)))

    def test_the_audit_names_the_role_not_the_login(self):
        import sage_lfq_check as slc
        rec = dict(accepted_record(), message="LFQ window warning.")
        status, msg, detail, _plain = slc.audit_finding(rec)
        self.assertNotIn(LOGIN, msg)
        self.assertNotIn(LOGIN, json.dumps(detail))
        self.assertEqual(rec["accepted"]["by"], LOGIN)             # the record itself untouched


class InterimLeavesThemOut(tdi.InterimBase):
    """deliver --mode interim (and every mode: _plan_tree, _plan_file) on a search folder holding
    the skill's real job-log names, the Core's records, and an acceptance by login."""

    def test_job_logs_and_core_records_stay_out_and_the_login_never_travels(self):
        c = self
        dda = c.dda
        for rel, text in (("diann_search_25404434.log", "node hpc-node-1 /home/x"),
                          ("diann_step3_123_4.log", "array task"),
                          ("fran_deposit.json", json.dumps({"staged_by": LOGIN})),
                          ("node_faults.json", json.dumps({"nodes": ["hpc-node-1"]}))):
            tcs.write(dda, rel, text)
        tcs.write(dda, "sage_lfq_check.json", json.dumps(accepted_record()))
        tcs.write(dda, "search_provenance.json",
                  json.dumps({"engine": "diann", "version": "2.7.0",
                              "sage_lfq_check": accepted_record()}))
        rc, out, p = c.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertTrue(out["verified"], out)
        delivered = []
        for dp, _, fs in os.walk(c.delivery):
            delivered += [os.path.relpath(os.path.join(dp, f), c.delivery) for f in fs]
        for rel in ("diann_search_25404434.log", "diann_step3_123_4.log", "fran_deposit.json",
                    "node_faults.json"):
            self.assertNotIn(f"search_dda/{rel}", delivered)
        man = tcs.read(os.path.join(c.delivery, "MANIFEST.txt"))
        self.assertIn("[SKIPPED] search_dda/diann_search_25404434.log -- Core-internal: a "
                      "scheduler job log", man)
        self.assertIn("[SKIPPED] search_dda/fran_deposit.json -- Core-internal", man)
        for rel in ("search_dda/sage_lfq_check.json", "search_dda/search_provenance.json"):
            self.assertIn(rel, delivered)
            self.assertNotIn(LOGIN, tcs.read(os.path.join(c.delivery, rel)), rel)
        for dp, _, fs in os.walk(c.delivery):
            for f in fs:
                with open(os.path.join(dp, f), encoding="utf-8", errors="replace") as fh:
                    self.assertNotIn(LOGIN, fh.read(), os.path.join(dp, f))
        self.assertIn(LOGIN, tcs.read(os.path.join(dda, "sage_lfq_check.json")))  # kept at home


class InputLinksAreNotDeliveredTwice(tdi.InterimBase):
    """C2 security L7: an interim delivery walks a search folder following links, so a DIA-NN
    search's input_links/ (links to whole .d runs, diann_parallel.INPUT_LINKS_DIR) was copied --
    each run's analysis.tdf and .tdf_bin a second time (the raw files are linked in raw/).
    input_links is scratch now (scratch_files, the one rule)."""

    def test_input_links_stay_out_with_their_reason(self):
        import diann_parallel
        import scratch_files
        self.assertEqual(scratch_files.INPUT_LINKS_DIR, diann_parallel.INPUT_LINKS_DIR)
        run = os.path.join(self.tmp, "elsewhere", "S1.d")
        os.makedirs(run)
        with open(os.path.join(run, "analysis.tdf"), "w") as fh:
            fh.write("raw data")
        os.makedirs(os.path.join(self.dda, "input_links"))
        os.symlink(run, os.path.join(self.dda, "input_links", "S1.d"))
        rc, out, p = self.interim("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertFalse(os.path.lexists(os.path.join(self.delivery, "search_dda", "input_links")))
        man = tcs.read(os.path.join(self.delivery, "MANIFEST.txt"))
        self.assertIn("[SKIPPED] search_dda/input_links/ -- excluded: the search's links to the "
                      "raw files", man)


class TheAcceptanceAsksOnlyWhy(unittest.TestCase):
    """Verify residual M1(a): the gate asked for --accept-lfq-window "<who accepted, and why>", so
    a typed name stayed verbatim in accepted.reason, AUDIT.md and commands.log. It asks for the
    reason only: who is the login that runs it, named to clients by role."""

    def test_the_refusal_and_the_help_ask_for_the_reason_only(self):
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            self.skipTest("pyarrow is needed for the Sage fixtures")
        import sage_lfq_check as slc
        import test_sage_lfq_check as fx
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "search")
            fx.write_sage_outputs(out, {"HeLa_1.mzML": 7.1, "HeLa_2.mzML": 7.7}, peaks=0)
            msg = slc.refusal(out, slc.check(out))
        self.assertIn('--accept-lfq-window "<why, in the user\'s words>"', msg)
        self.assertIn("put no name in the reason", msg)
        self.assertNotIn("<who", msg)
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "run_search.py"), "--help"],
                           capture_output=True, text=True, timeout=60)
        flat = " ".join(p.stdout.split())
        self.assertNotIn("give who and why", flat)
        self.assertIn("give WHY (no name", flat)


class TheDepositCarriesNoLogin(unittest.TestCase):
    """Verify residual M1(b): the PRIDE deposit's staging copied search_provenance.json verbatim
    (prepare_upload.sbatch `copy`), so an accepted-by login would reach a PUBLIC repository. The
    upload plan now points at a redacted copy (make_deposit.client_copies, staff.client_copy)."""

    def test_the_upload_plan_and_the_staging_use_the_redacted_copy(self):
        import csv
        import make_deposit as md
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as tmp:
            p = tdp.dia_session(tmp)
            with open(p["search_prov"]) as fh:
                prov = json.load(fh)
            prov["sage_lfq_check"] = accepted_record()
            write(p["search_prov"], json.dumps(prov))
            md.build(p["session_dir"], md.Manifest())
            with open(os.path.join(p["deposit_dir"], "files_to_upload.tsv"), newline="",
                      encoding="utf-8") as fh:
                rows = {r["upload_name"]: r for r in csv.DictReader(fh, delimiter="\t")}
            row = rows["search_provenance.json"]
            red = os.path.join(p["deposit_dir"], md.CLIENT_COPY_DIR, "search_provenance.json")
            self.assertEqual(row["source_path"], red)
            with open(red, encoding="utf-8") as fh:
                text = fh.read()
            self.assertNotIn(LOGIN, text)
            self.assertIn("redacted copy", row["notes"])
            with open(os.path.join(p["deposit_dir"], "prepare_upload.sbatch"),
                      encoding="utf-8") as fh:
                script = fh.read()
            self.assertIn(f"copy\t{red}\tsearch_provenance.json", script)
            self.assertNotIn(f"copy\t{p['search_prov']}\t", script)
            with open(p["search_prov"]) as fh:
                self.assertIn(LOGIN, fh.read())                 # the session keeps its copy


class TheZipLeavesThemOut(unittest.TestCase):
    def test_finalize_zip(self):
        from test_finalize_zip_logs import finish_fixture
        with tempfile.TemporaryDirectory() as tmp:
            env = dict(job_env(tmp), CLAUDE_CONFIG_DIR=os.path.join(tmp, "cc"))
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "init",
                                "--name", "Core internal zip", "--base",
                                os.path.join(tmp, "base")],
                               capture_output=True, text=True, env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            sdir = json.loads(r.stdout)["paths"]["session_dir"]
            finish_fixture(sdir)
            search = os.path.join(sdir, "output", "search")
            write(os.path.join(search, "diann_search_25404434.log"), "node hpc-node-1\n")
            write(os.path.join(search, "fran_deposit.json"), json.dumps({"staged_by": LOGIN}))
            write(os.path.join(search, "sage_lfq_check.json"), json.dumps(accepted_record()))
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "finalize",
                                "--dir", sdir, "--zip", "--no-deposit", "--no-notify"],
                               capture_output=True, text=True, env=env, timeout=600)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            res = json.loads(r.stdout)
            z = zipfile.ZipFile(res["zip"])
            names = z.namelist()
            self.assertFalse([n for n in names if n.endswith(("diann_search_25404434.log",
                                                              "fran_deposit.json"))], names)
            lfq = [n for n in names if n.endswith("output/search/sage_lfq_check.json")]
            self.assertEqual(len(lfq), 1, names)
            self.assertNotIn(LOGIN, z.read(lfq[0]).decode())
            self.assertEqual(sum(v for k, v in res["zip_excluded"].items()
                                 if k.startswith("Core-internal files outside logs/")), 2)


if __name__ == "__main__":
    unittest.main()
