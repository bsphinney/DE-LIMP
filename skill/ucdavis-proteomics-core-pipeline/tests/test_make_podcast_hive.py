#!/usr/bin/env python3
"""The podcast on HIVE and for Core runs (2.12.0).

Brett, 2026-10-15: "It should also create the podcast by default" -- for a CoreOmics submission;
internal Core work is offered it (round-2 decision). --core-run records the Core director's
standing approval, make_podcast.CORE_CLOUD_OK, the one place that text lives. Windows staff have
no usable local Python, so their podcast runs on HIVE: check through hive_exec.sh, then a SLURM
job (make_podcast.py job) that renders in Gemini voices, verifies, links and shares -- minutes of
network calls never run on the login node. On HIVE the key is the Core's shared one
(/quobyte/proteomics-grp/.config/gemini_key: owner brettsp, group proteomics-grp, mode 640),
read last, for a Core run only, and only while the Core's trust rule holds and accounts outside
the group cannot read it. setup.sh puts google-genai >= GENAI_MIN in the pipeline env.

Round 2 (2.12.0 reviews): every member of proteomics-grp can read that key, but only Core staff
may use it, so --core-run needs the person on the Core staff list (CORE_STAFF_FILE, trusted by
the Core's rule) -- checked by uid on HIVE, asked of HIVE from a laptop -- and `internal` needs
session.json to say internal Core work. Their login goes to the staff-only
logs/podcast.staff.json; podcast.json says "Core staff". On HIVE a personal key others can read
is refused. The job script and its log live in ~/.proteomics-pipeline/jobs/, a folder only its
owner can change (Core session folders are group-writable with no sticky bit: another member
could swap the script before sbatch reads it), written without following a link, and its
--output path is quoted.

Offline: google.genai is test_make_podcast's stand-in; the Core key and the staff list are
temporary files (CORE_KEY_FILE, the constant -- nothing in the environment can move it -- and
CORE_STAFF_FILE); the Core's admins are the account running the tests; CORE_ON_HIVE says where we
are; SLURM's queue choice is patched; the job's python is a stub that logs what it is asked to
run, in job_env's environment.
"""
import json
import os
import pwd
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import make_podcast as mp                             # noqa: E402
import slack_collab                                   # noqa: E402
import staff                                          # noqa: E402
import submission_report                              # noqa: E402
import test_make_podcast as T                         # noqa: E402  helpers only (no classes run here)
import test_setup_limpa as setup_limpa                # noqa: E402
from job_env import job_env                           # noqa: E402

ME = pwd.getpwuid(os.getuid()).pw_name
CORE_KEY = "AIza" + "SyCOREfake0123456789abcdefghijklmn"     # shaped like a key; not one
JOB_PY = r'''#!/bin/bash
# the job's python: answers job_python's google-genai probe, logs every other call
if [ "$1" = "-c" ]; then
  [ "${FAKE_GENAI_VERSION:-2.31.0}" = none ] && { echo "ModuleNotFoundError: google" >&2; exit 1; }
  echo "${FAKE_GENAI_VERSION:-2.31.0}"; exit 0
fi
( IFS='|'; echo "$2|${*:3}" ) >> "$FAKE_PY_LOG"      # one field per argument
case "$2" in
  render) exit "${FAKE_RENDER_RC:-0}" ;;
  link)   exit "${FAKE_LINK_RC:-0}" ;;
  share)  exit "${FAKE_SHARE_RC:-0}" ;;
esac
'''


def write(path, text, mode=None):
    T.write(path, text)
    if mode is not None:
        os.chmod(path, mode)
    return path


class CoreKeyFile(object):
    """HIVE as the Core sets it up, in a temporary folder: the Core's .config folder (mode 750;
    2750 on HIVE) with the shared key (640) and the Core staff list (644) naming the test runner,
    all the test runner's -- a Core admin here. CORE_ON_HIVE=1: this is HIVE."""

    def core_key_setup(self, staff_logins=None):
        self.cfg = os.path.join(self.d, "proteomics-grp", ".config")
        os.makedirs(self.cfg)
        os.chmod(self.cfg, 0o750)
        self.core_file = write(os.path.join(self.cfg, "gemini_key"), CORE_KEY + "\n", 0o640)
        self.staff_file = write(os.path.join(self.cfg, "core_staff.txt"),
                                "\n".join(staff_logins or ["someone_else", ME]) + "\n", 0o644)
        self.home = os.path.join(self.d, "home")
        os.makedirs(self.home)
        env = {k: v for k, v in os.environ.items() if k != "GEMINI_API_KEY"}
        env.update(HOME=self.home, PROTEOMICS_PIPELINE_HOME=os.path.join(self.home, "pp"),
                   CORE_ON_HIVE="1", **{staff.STAFF_FILE_ENV: self.staff_file})
        for p in (mock.patch.dict(os.environ, env, clear=True),
                  mock.patch.object(mp, "CORE_KEY_FILE", self.core_file),
                  mock.patch.object(slack_collab, "core_admins", lambda: (ME,))):
            p.start()
            self.addCleanup(p.stop)

    def user_key(self, key, mode=0o600):
        return write(os.path.join(self.home, ".config", "ucdavis-proteomics", "gemini_key"),
                     key + "\n", mode)

    def internal_session(self, title="Astral HeLa gradient test"):
        """session.json for internal Core work, and no submission record."""
        write(os.path.join(self.d, "session.json"), json.dumps({"analysis_for": {"internal": title}}))
        rec = os.path.join(self.d, "input", "submission.json")
        if os.path.exists(rec):
            os.remove(rec)


# ------------------------------------------------------------------------------- the key
class CoreKey(CoreKeyFile, unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.addCleanup(self._td.cleanup)
        self.d = self._td.name
        self.core_key_setup()

    def refused(self, core=True):
        with self.assertRaises(mp.RenderError) as cm:
            mp.read_key(core=core)
        msg = str(cm.exception)
        self.assertNotIn(CORE_KEY, msg)
        return msg

    def test_a_core_run_reads_the_core_key_last(self):
        self.assertEqual(mp.find_key(core=True)[:2], (CORE_KEY, self.core_file))
        self.user_key("user-key")
        self.assertEqual(mp.find_key(core=True)[:2], ("user-key", mp.KEY_FILE))
        with mock.patch.dict(os.environ, {"GEMINI_API_KEY": "env-key"}):
            self.assertEqual(mp.find_key(core=True)[:2], ("env-key", "GEMINI_API_KEY"))

    def test_never_outside_a_core_run(self):
        msg = self.refused(core=False)
        self.assertIn(f"{self.core_file} is used only for a Core run (--core-run)", msg)

    def test_group_readable_is_the_point_so_no_warning(self):
        rc, out, err = T.run("key", "--core-run")
        self.assertEqual(rc, 0, err)
        self.assertNotIn("[WARN]", err)
        got = json.loads(out)
        self.assertEqual((got["found"], got["source"]), (True, self.core_file))
        self.assertNotIn(CORE_KEY, out + err)

    def test_not_used_when_others_could_read_or_replace_it(self):
        cases = [
            (lambda: os.chmod(self.core_file, 0o644), "accounts outside its group can read it"),
            (lambda: os.chmod(self.core_file, 0o660), "people other than its owner can write it"),
            (lambda: os.chmod(self.cfg, 0o770), "others can write in its folder"),
            (lambda: setattr(slack_collab, "core_admins", lambda: ("someone",)),
             "not to a Core admin"),
        ]
        for spoil, why in cases:
            with self.subTest(why=why), \
                    mock.patch.object(slack_collab, "core_admins", slack_collab.core_admins):
                spoil()
                self.assertIn(why, self.refused())
                os.chmod(self.core_file, 0o640)
                os.chmod(self.cfg, 0o750)
        self.assertEqual(mp.read_key(core=True), CORE_KEY)              # restored: used again

    def test_a_link_or_an_empty_file_is_not_a_key(self):
        real = write(os.path.join(self.d, "elsewhere"), CORE_KEY, 0o640)
        os.remove(self.core_file)
        os.symlink(real, self.core_file)
        self.assertIn("is not readable", self.refused())
        os.remove(self.core_file)
        write(self.core_file, "\n", 0o640)
        self.assertIn("is empty", self.refused())

    def test_no_file_is_no_core_key(self):
        os.remove(self.core_file)
        self.assertEqual(mp.core_key(), (None, None))
        self.assertNotIn("Core", self.refused())

    def test_the_core_key_is_scrubbed_once_read(self):
        self.assertEqual(mp.read_key(core=True), CORE_KEY)
        self.assertNotIn(CORE_KEY, mp.scrub(f"400 Bad Request for key {CORE_KEY} (x)"))

    def test_key_without_a_core_run_says_why(self):
        rc, out, err = T.run("key")
        self.assertEqual(rc, 1)
        got = json.loads(out)
        self.assertEqual(got["found"], False)
        self.assertIn("used only for a Core run", got["note"])
        self.assertEqual(got["core_key"], "usable, for a Core run (--core-run)")

    # review L4: on HIVE every home folder is open to every account
    def test_on_hive_a_personal_key_others_can_read_is_refused(self):
        for mode in (0o644, 0o640, 0o604):
            with self.subTest(mode=oct(mode)):
                self.user_key("user-key", mode)
                key, src, why = mp.find_key(core=False)
                self.assertEqual((key, src), (None, None))
                self.assertIn(f"{mp.KEY_FILE} is not used", why)
                self.assertIn("other HIVE accounts can read it", why)
                self.assertIn("chmod 600", why)
                msg = self.refused(core=True)                # never quietly the Core's key instead
                self.assertIn("is not used", msg)
                self.assertNotIn("user-key", msg)
        self.user_key("user-key", 0o600)
        self.assertEqual(mp.find_key(core=False)[:2], ("user-key", mp.KEY_FILE))

    def test_on_a_laptop_it_is_still_only_a_warning(self):
        self.user_key("user-key", 0o644)
        with mock.patch.dict(os.environ, {"CORE_ON_HIVE": "0"}):
            rc, out, err = T.run("key")
        self.assertEqual((rc, json.loads(out)["source"]), (0, mp.KEY_FILE))
        self.assertIn("[WARN]", err)


# ------------------------------------------------------------------------------- consent
class CoreRunConsent(CoreKeyFile, T.Workspace):
    def setUp(self):
        super().setUp()
        self.core_key_setup()
        for p in (mock.patch.object(mp, "_sleep", lambda s: None),
                  mock.patch.object(mp, "encode_aac", return_value=(None, "test"))):
            p.start()
            self.addCleanup(p.stop)
        self.assertEqual(self.check()[0], 0)

    def core_session(self, prot="PROT_0807"):
        os.makedirs(os.path.join(self.d, "input"), exist_ok=True)
        submission_report.attach(self.d, {"internal_id": prot})

    def render(self, *extra):
        mods, log = T.fake_genai(lambda m, c, k: T.pcm_for(c))
        keys = []

        class Client(mods["google.genai"].Client):
            def __init__(self, api_key):
                keys.append(api_key)
                super().__init__(api_key)
        mods["google.genai"].Client = Client
        with mock.patch.dict(sys.modules, mods):
            rc, out, err = T.run("render", self.script, "--tts", "gemini", "--no-verify",
                                 "--no-preflight", *extra)
        return rc, out, err, log, keys

    def manifest(self):
        return json.loads(T.read(os.path.join(self.pod, "podcast.json")))

    def staff_record(self):
        return staff.read(os.path.join(self.d, "logs", "podcast.staff.json"))

    def assert_no_key_in_files(self):
        for root, _, files in os.walk(self.d):
            for f in files:
                if os.path.join(root, f) != self.core_file:
                    with open(os.path.join(root, f), "rb") as fh:
                        self.assertNotIn(CORE_KEY.encode(), fh.read(), f)

    def test_a_submission_records_the_standing_approval_and_uses_the_core_key(self):
        self.core_session()
        rc, out, err, log, keys = self.render("--core-run")
        self.assertEqual(rc, 0, err)
        self.assertEqual(keys, [CORE_KEY])
        man = self.manifest()
        self.assertEqual(man["cloud_tts_consent"], mp.CORE_CLOUD_OK)
        self.assertEqual(man["core_run"], {"kind": "submission", "prot": "PROT_0807",
                                           "by": staff.ROLE, "staff_check":
                                           f"on the Core staff list ({self.staff_file})"})
        self.assertNotIn(ME, json.dumps(man), "a login never goes in podcast.json")
        rec = self.staff_record()
        self.assertEqual([(e["event"], e["by"], e["kind"], e["prot"]) for e in rec],
                         [("render", ME, "submission", "PROT_0807")])
        self.assertTrue(staff.is_staff_only(os.path.join(self.d, "logs", "podcast.staff.json")))
        self.assert_no_key_in_files()

    def test_not_a_core_run_sends_nothing(self):
        rc, out, err, log, keys = self.render("--core-run")
        self.assertEqual((rc, log["clients"], keys), (2, 0, []))
        self.assertIn("answers no CoreOmics submission", err)
        self.assertIn("--core-run internal", err)
        self.assertIn("--cloud-ok", err)
        self.assertFalse(os.path.exists(os.path.join(self.pod, "podcast.json")))

    # review M3: the key is readable by all of proteomics-grp; the approval is for Core staff
    def test_only_core_staff_may_claim_a_core_run(self):
        self.core_session()
        write(self.staff_file, "someone_else\nanother\n", 0o644)
        for env in ({}, {"USER": "someone_else", "LOGNAME": "someone_else"}):   # $USER is no proof
            with self.subTest(env=env), mock.patch.dict(os.environ, env):
                rc, out, err, log, keys = self.render("--core-run")
                self.assertEqual((rc, log["clients"], keys), (2, 0, []))
                self.assertIn(f"--core-run is for Core staff: {ME} is not on the Core staff "
                              f"list ({self.staff_file})", err)
        self.assertFalse(os.path.exists(os.path.join(self.pod, "podcast.json")))

    def test_an_untrusted_or_missing_staff_list_is_a_no(self):
        self.core_session()
        cases = [(lambda: os.chmod(self.staff_file, 0o666), "is not trusted"),
                 (lambda: os.remove(self.staff_file), "does not exist")]
        for spoil, why in cases:
            with self.subTest(why=why):
                spoil()
                rc, out, err, log, keys = self.render("--core-run")
                self.assertEqual((rc, log["clients"]), (2, 0))
                self.assertIn("the Core staff list cannot be used", err)
                self.assertIn(why, err)

    def test_internal_work_needs_session_json_to_say_so(self):
        rc, out, err, log, keys = self.render("--core-run", "internal")
        self.assertEqual((rc, log["clients"]), (2, 0))
        self.assertIn("does not say it is internal Core work", err)
        write(os.path.join(self.d, "session.json"),
              json.dumps({"analysis_for": {"submissions": ["PROT_0807"]}}))
        self.assertEqual(self.render("--core-run", "internal")[0], 2)
        self.internal_session()
        rc, out, err, log, keys = self.render("--core-run", "internal")
        self.assertEqual(rc, 0, err)
        man = self.manifest()
        self.assertEqual((man["cloud_tts_consent"], man["core_run"]["kind"],
                          man["core_run"]["title"], man["core_run"]["by"]),
                         (mp.CORE_CLOUD_OK, "internal", "Astral HeLa gradient test", staff.ROLE))
        self.assertEqual(self.staff_record()[-1]["by"], ME)

    # review E3: client data must never be labelled the Core's own
    def test_internal_is_refused_on_a_submissions_session(self):
        self.core_session()                                    # attached: PROT_0807
        write(os.path.join(self.d, "session.json"),
              json.dumps(dict(json.loads(T.read(os.path.join(self.d, "session.json"))),
                              analysis_for={"internal": "looks internal"})))
        rc, out, err, log, keys = self.render("--core-run", "internal")
        self.assertEqual((rc, log["clients"]), (2, 0))
        self.assertIn("this session is for PROT_0807, a CoreOmics submission", err)
        self.assertIn("pass --core-run.", err)

    def test_a_named_but_unattached_submission_says_attach(self):
        write(os.path.join(self.d, "session.json"),
              json.dumps({"analysis_for": {"submissions": ["PROT_0807"]}}))
        rc, out, err, log, keys = self.render("--core-run", "internal")
        self.assertEqual(rc, 2)
        self.assertIn("this session is for PROT_0807", err)
        rc, out, err, log, keys = self.render("--core-run")
        self.assertEqual((rc, log["clients"]), (2, 0))
        self.assertIn("no submission record is attached: run submission_report.py attach", err)
        self.assertNotIn("--core-run internal", err)

    def test_a_laptop_asks_hive_about_the_account_its_key_reaches(self):
        self.core_session()
        calls = []

        def fake_run(argv, **kw):
            calls.append(argv)
            return subprocess.CompletedProcess(argv, 0, stdout=self.hive_answer, stderr="")
        answers = [
            (json.dumps({"ok": True, "by": "brettsp", "checked": "on the Core staff list (/q)"}),
             0, None),
            (json.dumps({"ok": True, "by": "brettsp", "checked": "NOT VERIFIED: no list"}), 2,
             "NOT VERIFIED"),
            (json.dumps({"ok": False, "error": "--by 'x' is not on the Core staff list"}), 2,
             "is not on the Core staff list"),
            ("ssh: connect to host hive.hpc.ucdavis.edu port 22: timed out", 2,
             "could not be asked from here"),
        ]
        for self.hive_answer, want_rc, said in answers:
            with self.subTest(said=said), mock.patch.dict(os.environ, {"CORE_ON_HIVE": "0"}), \
                    mock.patch.object(mp.subprocess, "run", fake_run), \
                    mock.patch.dict(os.environ, {"GEMINI_API_KEY": "laptop-key"}):
                mods, log = T.fake_genai(lambda m, c, k: T.pcm_for(c))
                with mock.patch.dict(sys.modules, mods):
                    rc, out, err = T.run("render", self.script, "--tts", "gemini", "--no-verify",
                                         "--no-preflight", "--core-run")
                self.assertEqual(rc, want_rc, err)
                if said:
                    self.assertIn(said, err)
        self.assertEqual(calls[0][-1], 'python3 ~/proteomics-pipeline/scripts/staff.py check '
                                       '--by "$(id -un)"')
        self.assertTrue(calls[0][1].endswith("hive_exec.sh"))
        self.assertEqual(self.staff_record()[0]["by"], "brettsp")

    def test_cloud_ok_never_reaches_the_core_key(self):
        self.core_session()
        rc, out, err, log, keys = self.render("--cloud-ok", "the lab's PI, 2026-10-15")
        self.assertEqual((rc, log["clients"]), (1, 0))
        self.assertIn("used only for a Core run (--core-run)", err)

    def test_one_consent_or_the_other(self):
        with self.assertRaises(SystemExit):
            T.run("render", self.script, "--tts", "gemini", "--core-run", "--cloud-ok", "x")

    def test_verify_under_a_core_run(self):
        self.core_session()
        mp.BACKENDS["fake"] = T.FakeTTS
        self.addCleanup(mp.BACKENDS.pop, "fake", None)
        self.assertEqual(T.run("render", self.script, "--tts", "fake")[0], 0)
        heard = "\n".join(t.text for t in mp.parse_script(self.script).turns())
        mods, log = T.fake_genai(lambda m, c, k: T.pcm_for(c), transcribe=lambda m, c: heard)
        with mock.patch.dict(sys.modules, mods), \
                mock.patch.object(mp, "_verify_audio", return_value=(b"AAC", "audio/aac", "t")):
            rc, out, err = T.run("verify", self.script, "--core-run")
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.manifest()["verify"]["cloud_consent"], mp.CORE_CLOUD_OK)
        self.assertEqual([e["event"] for e in self.staff_record()], ["verify"])

    def test_the_approval_text_lives_in_one_place(self):
        hits = []
        for root, _, files in os.walk(SKILL):
            for f in files:
                p = os.path.join(root, f)
                if f.endswith((".md", ".py", ".sh", ".R")) and os.path.abspath(p) not in (
                        os.path.abspath(mp.__file__), os.path.abspath(__file__)):
                    with open(p, encoding="utf-8", errors="replace") as fh:
                        if "standing approval for Core runs" in fh.read():
                            hits.append(os.path.relpath(p, SKILL))
        self.assertEqual(hits, [], "point at make_podcast.CORE_CLOUD_OK; never retype it")


# ------------------------------------------------------------------------------- verify
class StaleVerify(T.Workspace):
    """Review E2: a re-render whose verify failed (ASR 429) left the previous render's verify.txt
    saying OK beside the new audio, rc 0, and record_run copies verify.txt. Now render clears
    both verify files when it writes podcast.json, and a verify that does not run says so."""

    def setUp(self):
        super().setUp()
        for p in (mock.patch.dict(os.environ, {"GEMINI_API_KEY": T.FAKE_KEY}),
                  mock.patch.object(mp, "_sleep", lambda s: None),
                  mock.patch.object(mp, "encode_aac", return_value=(None, "test")),
                  mock.patch.object(mp, "_verify_audio", return_value=(b"AAC", "audio/aac", "t"))):
            p.start()
            self.addCleanup(p.stop)
        self.assertEqual(self.check()[0], 0)
        self.heard = "\n".join(t.text for t in mp.parse_script(self.script).turns())

    def render(self, transcribe, *extra):
        mods, log = T.fake_genai(lambda m, c, k: T.pcm_for(c), transcribe=transcribe)
        with mock.patch.dict(sys.modules, mods):
            return T.run("render", self.script, "--tts", "gemini", "--cloud-ok", "PI, today",
                         "--no-preflight", *extra)

    def files(self):
        verify = os.path.join(self.pod, "verify.txt")
        man = json.loads(T.read(os.path.join(self.pod, "podcast.json")))
        return (T.read(verify) if os.path.exists(verify) else None,
                os.path.exists(os.path.join(self.pod, "verify_transcript.txt")), man.get("verify"))

    def test_a_failed_verify_never_leaves_the_last_ok(self):
        self.assertEqual(self.render(lambda m, c: self.heard)[0], 0)
        txt, transcript, block = self.files()
        self.assertIn("Podcast audio check (ASR round trip): OK", txt)
        self.assertTrue(transcript)

        def quota(m, c):
            raise T.ApiError(429, "RESOURCE_EXHAUSTED GenerateRequestsPerDayPerProjectPerModel")
        rc, out, err = self.render(quota)
        self.assertEqual(rc, 0, err)                          # the render itself is fine
        txt, transcript, block = self.files()
        self.assertIn("Podcast audio check (ASR round trip): NOT RUN", txt)
        self.assertIn("reason: verify failed:", txt)
        self.assertNotIn(": OK", txt)
        self.assertFalse(transcript, "the last render's transcript is gone")
        self.assertEqual(block["status"], "NOT RUN")
        self.assertNotIn(T.FAKE_KEY, txt + json.dumps(block))

    def test_no_verify_and_offline_renders_say_it_was_not_run(self):
        self.assertEqual(self.render(lambda m, c: self.heard)[0], 0)
        self.assertEqual(self.render(lambda m, c: self.heard, "--no-verify")[0], 0)
        txt, transcript, block = self.files()
        self.assertIn("NOT RUN\nreason: --no-verify was given", txt)
        self.assertEqual((transcript, block["status"]), (False, "NOT RUN"))
        mp.BACKENDS["fake"] = T.FakeTTS
        self.addCleanup(mp.BACKENDS.pop, "fake", None)
        self.assertEqual(T.run("render", self.script, "--tts", "fake")[0], 0)
        txt, transcript, block = self.files()
        self.assertIn("reason: an offline render (fake) sends nothing", txt)
        self.assertEqual(block["status"], "NOT RUN")


# ------------------------------------------------------------------------------- the job
class Job(CoreKeyFile, T.Workspace):
    QUEUE = ("high", "genome-center-grp", "genome-center-grp-high-qos")

    def setUp(self):
        super().setUp()
        self.core_key_setup()
        os.makedirs(os.path.join(self.d, "logs"))
        os.makedirs(os.path.join(self.d, "input"))
        submission_report.attach(self.d, {"internal_id": "PROT_0807"})
        self.py = setup_limpa._exe(os.path.join(self.d, "bin", "py"), JOB_PY)
        self.py_log = os.path.join(self.d, "py.log")
        import run_search
        self.queue = mock.patch.object(run_search, "slurm_queue", return_value=self.QUEUE)
        self.queue.start()
        self.addCleanup(self.queue.stop)
        self.assertEqual(self.check()[0], 0)
        self.rd = os.path.realpath(self.d)             # the job names the output folder's real path
        self.jobs = os.path.join(os.environ["PROTEOMICS_PIPELINE_HOME"], "jobs")
        sha = mp.sha256_bytes(self.rd.encode("utf-8"))[:8]
        self.job_path = os.path.join(self.jobs, f"podcast_PROT_0807_{sha}.sh")
        self.log_path = os.path.join(self.jobs, f"podcast_PROT_0807_{sha}_%j.log")

    def job(self, *extra, python=True):
        args = ["job", self.script] + list(extra) + (["--python", self.py] if python else [])
        return T.run(*args)

    def run_job(self, **env):
        e = job_env(self.d, FAKE_PY_LOG=self.py_log, **env)
        r = subprocess.run(["bash", self.job_path], capture_output=True, text=True, env=e,
                           timeout=60)
        self.job_stderr = r.stderr
        return r.returncode, [ln.split("|")[0] for ln in setup_limpa._lines(self.py_log)]

    def job_args(self):
        """What the job passed to each step, one list of arguments per call."""
        return [ln.split("|") for ln in setup_limpa._lines(self.py_log)]

    def test_the_job_script(self):
        rc, out, err = self.job("--core-run")
        self.assertEqual(rc, 0, err)
        got = json.loads(out)
        self.assertEqual(got["job"], self.job_path)
        self.assertEqual(got["submit"], f"sbatch {self.job_path}")
        self.assertEqual(got["log"], self.log_path.replace("%j", "<jobid>"))
        self.assertEqual((got["python"], got["google_genai"], got["key"]),
                         (self.py, "2.31.0", self.core_file))
        text = T.read(self.job_path)
        lines = text.splitlines()
        self.assertEqual(lines[0], "#!/bin/bash -l", "a login shell: HIVE's `module` command")
        # sized for the local voice, which takes over in this job when Gemini cannot be used
        for want in ("#SBATCH --open-mode=append",              # a requeued job keeps its log
                     "#SBATCH --ntasks=1", "#SBATCH --cpus-per-task=16", "#SBATCH --mem=24G",
                     "#SBATCH --time=2:00:00",
                     "#SBATCH --partition=high", "#SBATCH --account=genome-center-grp",
                     "#SBATCH --qos=genome-center-grp-high-qos",
                     f'#SBATCH --output="{self.log_path}"', "module load ffmpeg/7.0.2",
                     "module load apptainer"):
            self.assertIn(want, text)
        # ucd-hpc: modules.sh is sourced before the first `module`
        self.assertLess(text.index(". /etc/profile.d/modules.sh"), text.index("module load"))
        self.assertNotIn("--requeue", text)
        steps = [ln for ln in lines if ln.startswith('"$PY" "$MP"')]
        self.assertEqual([s.split()[2] for s in steps], ["render", "link", "share"])
        self.assertIn(f"render {self.script} --tts gemini --core-run || exit $?", steps[0])
        self.assertEqual(steps[1], f'"$PY" "$MP" link {os.path.join(self.rd, "output")} || exit $?')
        self.assertIn(f"MP={os.path.join(SCRIPTS, 'make_podcast.py')}", text)
        self.assertNotIn(mp.CORE_CLOUD_OK, text, "render records it; the job names --core-run")
        self.assertNotIn(CORE_KEY, text)
        self.assertEqual(subprocess.run(["bash", "-n", self.job_path]).returncode, 0)
        self.assertFalse(os.listdir(os.path.join(self.d, "logs")), "nothing in the shared session")

    # review M2: a group-writable session folder let another member swap the script
    def test_the_job_is_private_and_never_written_through_a_link(self):
        self.assertEqual(self.job("--core-run")[0], 0)
        self.assertEqual(stat.S_IMODE(os.stat(self.jobs).st_mode), 0o700)
        self.assertEqual(stat.S_IMODE(os.stat(self.job_path).st_mode), 0o700)
        self.assertEqual([f for f in os.listdir(self.jobs) if f.endswith(".tmp")], [])
        victim = write(os.path.join(self.d, "victim.txt"), "keep me\n", 0o600)
        os.remove(self.job_path)
        os.symlink(victim, self.job_path)
        rc, out, err = self.job("--core-run")
        self.assertEqual(rc, 2)
        self.assertIn("is a link: not writing through it", err)
        self.assertEqual(T.read(victim), "keep me\n")
        self.assertTrue(os.path.islink(self.job_path))

    def test_a_jobs_folder_others_could_change_is_refused_or_made_private(self):
        os.makedirs(self.jobs)
        os.chmod(self.jobs, 0o777)
        self.assertEqual(self.job("--core-run")[0], 0)
        self.assertEqual(stat.S_IMODE(os.stat(self.jobs).st_mode), 0o700)      # made private
        shutil.rmtree(self.jobs)
        os.symlink(self.d, self.jobs)                                       # a link: refused
        rc, out, err = self.job("--core-run")
        self.assertEqual(rc, 2)
        self.assertIn("is not a folder (a link?)", err)
        os.remove(self.jobs)
        pp = os.environ["PROTEOMICS_PIPELINE_HOME"]
        os.chmod(pp, 0o777)                                                 # its parent: refused
        try:
            rc, out, err = self.job("--core-run")
        finally:
            os.chmod(pp, 0o755)
        self.assertEqual(rc, 2)
        self.assertIn("can be changed by other accounts", err)

    # review: an unquoted --output with a space ends the directive (sbatch --test-only, HIVE)
    def test_the_log_path_is_quoted_for_sbatch(self):
        home = os.path.join(self.d, "a home with spaces")
        with mock.patch.dict(os.environ, {"PROTEOMICS_PIPELINE_HOME": home}):
            rc, out, err = self.job("--core-run")
        self.assertEqual(rc, 0, err)
        got = json.loads(out)
        text = T.read(got["job"])
        self.assertIn(f'#SBATCH --output="{os.path.join(home, "jobs")}/', text)
        self.assertIn("'", got["submit"], "the submit line quotes the path too")

    # review E1: 2,557 Core service folders have a space in their names
    def test_a_session_path_with_a_space(self):
        sd = os.path.join(self.d, "Core service", "PROT 0807 session")
        out = os.path.join(sd, "output")
        T.write(os.path.join(out, "AI_Analysis_Report.md"), T.REPORT)
        T.write(os.path.join(out, "Analysis_Report.html"), T.REPORT_HTML)
        script = os.path.join(out, "podcast", "podcast_script.md")
        T.write(script, T.script_text())
        rc, o, e = T.run("check", script, "--source", os.path.join(out, "AI_Analysis_Report.md"))
        self.assertEqual(rc, 0, e)
        T.meaning_pass(script)
        os.makedirs(os.path.join(sd, "input"))
        submission_report.attach(sd, {"internal_id": "PROT_0807"})
        home = os.path.join(self.d, "pipeline home")
        with mock.patch.dict(os.environ, {"PROTEOMICS_PIPELINE_HOME": home}):
            rc, o, err = T.run("job", script, "--core-run", "--python", self.py)
        self.assertEqual(rc, 0, err)
        got = json.loads(o)
        self.assertEqual(got["submit"], "sbatch " + shlex.quote(got["job"]))
        self.assertIn(f'#SBATCH --output="{os.path.join(home, "jobs")}/', T.read(got["job"]))
        self.job_path = got["job"]
        self.assertEqual(self.run_job(), (0, ["render", "link", "share"]))
        render, link, share = self.job_args()
        self.assertEqual(render[1:5], [script, "--tts", "gemini", "--core-run"])
        self.assertEqual(link[1:], [os.path.realpath(out)])
        self.assertEqual(share[1:], [os.path.realpath(out)])

    # review E4: a failed module load used to ship a 58 MB podcast.wav without a word
    def test_no_ffmpeg_is_said_in_the_log(self):
        self.assertEqual(self.job("--core-run")[0], 0)
        only_bash = os.path.join(self.d, "only_bash")                   # no ffmpeg, no afconvert
        os.makedirs(only_bash)
        os.symlink(shutil.which("bash"), os.path.join(only_bash, "bash"))
        self.assertEqual(self.run_job(PATH=only_bash)[0], 0)
        self.assertIn("[WARN] no ffmpeg here (module load ffmpeg/7.0.2 failed)", self.job_stderr)
        self.assertIn("podcast.wav", self.job_stderr)

    def test_a_preemptible_queue_requeues(self):
        self.queue.stop()
        import run_search
        with mock.patch.object(run_search, "slurm_queue",
                               return_value=("low", "publicgrp", "publicgrp-low-qos")):
            self.assertEqual(self.job("--core-run")[0], 0)
        self.queue.start()
        self.assertIn("#SBATCH --requeue", T.read(self.job_path))

    def test_it_runs_render_link_share_and_stops_at_a_failure(self):
        self.internal_session()
        rc, out, err = self.job("--core-run", "internal")
        self.assertEqual(rc, 0, err)
        self.job_path = json.loads(out)["job"]
        self.assertIn("--core-run internal ||", T.read(self.job_path))
        self.assertEqual(self.run_job(), (0, ["render", "link", "share"]))
        os.remove(self.py_log)
        self.assertEqual(self.run_job(FAKE_RENDER_RC="1"), (1, ["render"]))
        os.remove(self.py_log)
        self.assertEqual(self.run_job(FAKE_LINK_RC="3"), (3, ["render", "link"]))

    def test_refuses_someone_not_on_the_staff_list(self):
        write(self.staff_file, "someone_else\n", 0o644)
        rc, out, err = self.job("--core-run")
        self.assertEqual(rc, 2)
        self.assertIn("--core-run is for Core staff", err)
        self.assertFalse(os.path.exists(self.job_path))

    def test_refuses_a_script_whose_check_does_not_hold(self):
        T.write(self.script, T.read(self.script).replace("30 IPs", "31 IPs"))
        rc, out, err = self.job("--core-run")
        self.assertEqual(rc, 2)
        self.assertIn("the script changed after check.txt was written", err)
        self.assertFalse(os.path.exists(self.job_path))

    def test_refuses_without_consent(self):
        rc, out, err = self.job()
        self.assertEqual(rc, 2)
        self.assertIn("--core-run for a Core run", err)
        rc, out, err = self.job("--cloud-ok", "declined")
        self.assertEqual(rc, 2)
        self.assertIn("is a refusal", err)
        self.assertFalse(os.path.exists(self.job_path))

    def test_refuses_without_a_new_enough_google_genai(self):
        for version, said in (("2.20.0", "google-genai 2.20.0"), ("none", "no google-genai")):
            with self.subTest(version=version), \
                    mock.patch.dict(os.environ, {"FAKE_GENAI_VERSION": version}):
                rc, out, err = self.job("--core-run")
                self.assertEqual(rc, 2)
                self.assertIn(f"no python here has google-genai >= {mp.GENAI_MIN}", err)
                self.assertIn(said, err)
                self.assertIn("setup.sh", err)
        self.assertFalse(os.path.exists(self.job_path))

    def test_the_pipeline_env_is_the_default_python(self):
        home = os.environ["PROTEOMICS_PIPELINE_HOME"]
        write(os.path.join(home, "setup.json"), json.dumps({"python": self.py}))
        rc, out, err = self.job("--core-run", python=False)
        self.assertEqual(rc, 0, err)
        self.assertEqual(json.loads(out)["python"], self.py)

    def test_refuses_without_a_key(self):
        os.remove(self.core_file)
        rc, out, err = self.job("--core-run")
        self.assertEqual(rc, 2)
        self.assertIn("no Gemini API key for the job", err)
        self.assertFalse(os.path.exists(self.job_path))

    def test_a_consent_note_is_quoted_for_the_shell(self):
        note = "Dr O'Neil's lab; agreed 2026-10-15 $(id)"
        with mock.patch.dict(os.environ, {"GEMINI_API_KEY": "env-key"}):
            rc, out, err = self.job("--cloud-ok", note, "--model", "gemini-3.8-flash-tts")
        self.assertEqual(rc, 0, err)
        self.assertEqual(json.loads(out)["key"], "GEMINI_API_KEY")
        self.job_path = json.loads(out)["job"]          # not a Core run: named for the session
        self.assertEqual(self.run_job()[0], 0)
        render = self.job_args()[0]
        self.assertEqual(render[-4:], ["--cloud-ok", note, "--model", "gemini-3.8-flash-tts"])

    def test_a_missing_script_is_one_line(self):
        rc, out, err = T.run("job", os.path.join(self.pod, "nope.md"), "--core-run")
        self.assertEqual(rc, 2)
        self.assertIn("is not there", err)
        self.assertNotIn("Traceback", err)

    def test_the_script_must_sit_in_output_podcast(self):
        other = os.path.join(self.d, "elsewhere", "podcast_script.md")
        T.write(other, T.read(self.script))
        rc, out, err = T.run("job", other, "--core-run", "--python", self.py)
        self.assertEqual(rc, 2)
        self.assertIn("must be output/podcast/podcast_script.md", err)


# ------------------------------------------------------------------------------- setup.sh
class GenaiFloor(unittest.TestCase):
    def test_setup_sh_installs_the_floor_make_podcast_needs(self):
        with open(os.path.join(SCRIPTS, "setup.sh")) as fh:
            got = re.search(r'(?m)^GENAI_MIN="([^"]+)"$', fh.read())
        self.assertEqual(got and got.group(1), mp.GENAI_MIN)


class SetupInstallsGenai(setup_limpa.setup_t.SetupCheckHarness):
    """setup.sh step 2d with test_setup_limpa's stubs (no conda, no network)."""
    run_setup = setup_limpa.SetupLimpaGate.run_setup
    prebuilt_env = setup_limpa.SetupLimpaGate.prebuilt_env

    def setUp(self):            # SetupLimpaGate.setUp, whose super() needs that class
        super().setUp()
        for tool in ("chmod", "ln", "rm"):
            os.symlink(shutil.which(tool), os.path.join(self.sys, tool))
        setup_limpa._exe(os.path.join(self.sys, "micromamba"), setup_limpa.FAKE_MICROMAMBA)
        self.ensure = setup_limpa._exe(os.path.join(self.d, "ensure_dotnet8_stub.sh"),
                                       setup_limpa.FAKE_ENSURE)
        self.dotnet_root = setup_limpa.fake_dotnet_root(
            os.path.join(self.d, ".proteomics-pipeline", "dotnet8"))
        self.mm_log = os.path.join(self.d, "mm.log")
        self.r_log = os.path.join(self.d, "r.log")
        self.state = os.path.join(self.d, "limpa_version")
        self.prefix = os.path.join(self.home, "micromamba", "envs", "proteomics-pipeline")

    def genai_record(self, version=None):
        meta = os.path.join(self.prefix, "conda-meta")
        for f in os.listdir(meta):
            if f.startswith("google-genai-"):
                os.remove(os.path.join(meta, f))
        if version:
            write(os.path.join(meta, f"google-genai-{version}-pyh0000000_0.json"), "{}")

    def genai_installs(self):
        return [c for c in setup_limpa._lines(self.mm_log)
                if c.startswith("install") and "google-genai" in c]

    def test_an_env_without_it_gets_it_from_conda_forge(self):
        self.prebuilt_env("1.4.1")
        self.genai_record(None)
        r, s = self.run_setup()
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        got = self.genai_installs()
        self.assertEqual(len(got), 1, setup_limpa._lines(self.mm_log))
        # --freeze-installed: the env is shared, and the podcast's package moves nothing else
        self.assertIn(f"-p {self.prefix} --freeze-installed --override-channels -c conda-forge "
                      f"-c bioconda google-genai>={mp.GENAI_MIN}", got[0])
        # setup.json says which skill's setup.sh wrote it, and when (skill_version.sh
        # --check-hive reads them on HIVE)
        # (a copy of setup.sh without plugin.json beside it: skill_version.sh's tagged unknown)
        self.assertRegex(s["skill_version"], r"^(\d+\.\d+\.\d+|\(unknown)")
        self.assertIn("written", s)               # (this harness's PATH has no `date`)
        self.assertEqual(s["podcast"], {"google_genai": "2.31.0",
                                        "required": f">= {mp.GENAI_MIN}", "ok": True})
        self.assertFalse([n for n in s["notes"] if "google-genai" in n], s["notes"])
        r, s = self.run_setup()                                   # and never again
        self.assertEqual(len(self.genai_installs()), 1)

    def test_an_older_one_is_upgraded(self):
        self.prebuilt_env("1.4.1")
        self.genai_record("2.20.0")
        r, s = self.run_setup()
        self.assertEqual(len(self.genai_installs()), 1)
        self.assertIs(s["podcast"]["ok"], True)

    def test_check_reports_it_and_installs_nothing(self):
        self.prebuilt_env("1.4.1")
        self.genai_record("2.20.0")
        r, s = self.run_setup("--check")
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertEqual(self.genai_installs(), [])
        self.assertEqual(s["podcast"], {"google_genai": "2.20.0",
                                        "required": f">= {mp.GENAI_MIN}", "ok": False})
        note = [n for n in s["notes"] if "google-genai" in n]
        self.assertEqual(len(note), 1, s["notes"])
        self.assertIn(f"google-genai 2.20.0 is older than {mp.GENAI_MIN}", note[0])
        self.assertIn(f"'google-genai>={mp.GENAI_MIN}'", note[0])
        self.assertIn("--freeze-installed", note[0])
        self.assertIs(s["check_only"], True)

    def test_a_failed_install_is_a_note_not_a_failure(self):
        self.prebuilt_env("1.4.1")
        self.genai_record(None)
        r, s = self.run_setup(FAKE_MM_INSTALL_FAIL="1")
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        self.assertIs(s["podcast"]["ok"], False)
        self.assertIn("only the podcast needs it", r.stderr)
        self.assertTrue([n for n in s["notes"] if "google-genai is not installed" in n])

    def test_a_current_one_is_left_alone(self):
        self.prebuilt_env("1.4.1")
        r, s = self.run_setup()
        self.assertEqual(self.genai_installs(), [])
        self.assertIs(s["podcast"]["ok"], True)


if __name__ == "__main__":
    unittest.main()
