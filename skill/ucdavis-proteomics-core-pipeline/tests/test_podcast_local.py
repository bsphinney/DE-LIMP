#!/usr/bin/env python3
"""
The local voice (make_podcast.py --tts local, and the fallback when Gemini cannot be used):
Kokoro-82M speaks the episode and Whisper hears it back, in the Core's pinned containers on HIVE
(podcast_local.py). Brett Phinney (Core director) decided it on 2026-10-20: Gemini stays the
default; Kokoro takes over -- making EVERY chunk again, so one episode never mixes voices --
only when Gemini fails in a way waiting cannot fix (a 402 "prepay credits depleted", the key
refused, no quota on any model, every model stalled). The voices are bf_emma (Maya) and af_heart
(the second host); both are women, so a Kokoro render says "Lena" where the script says "Leo".

No network and no HIVE: the containers are a fake `apptainer` (and `sha256sum`) on PATH, run by
podcast_local's real job scripts through the Here runner, sized as a 16-CPU, 24 GB job.
"""
import json
import os
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import test_make_podcast as T  # noqa: E402
from test_make_podcast import (FAKE_KEY, ApiError, fake_genai, pcm_for, read, rpc_error,  # noqa: E402
                               run, mp, write)

import podcast_local as pl  # noqa: E402

STALLS = "gemini-3.8-flash-tts"
ANSWERS = "gemini-2.5-flash-preview-tts"
JOB_ENV = {"SLURM_JOB_ID": "1", "SLURM_CPUS_PER_TASK": "16", "SLURM_MEM_PER_NODE": "24576"}
PREPAY = {"error": {"code": 402, "status": "RESOURCE_EXHAUSTED",
                    "message": "Your prepayment credits are depleted. Please go to AI Studio "
                               "to manage your project and billing."}}

FAKE_APPTAINER = r'''#!PYTHON
"""A stand-in for `apptainer exec ... IMAGE python3 /opt/kokoro/synth.py|/opt/whisper/transcribe.py`."""
import array, json, os, sys, wave
argv = sys.argv[1:]
if os.environ.get("FAKE_APPTAINER_LOG"):
    with open(os.environ["FAKE_APPTAINER_LOG"], "a") as fh:
        fh.write(json.dumps(argv) + "\n")
    if os.environ.get("GEMINI_API_KEY"):
        open(os.environ["FAKE_APPTAINER_LOG"] + ".key", "w").write("the key reached a container")
i = argv.index("python3")
image, prog, args = argv[i - 1], argv[i + 1], argv[i + 2:]
def opt(name):
    return args[args.index(name) + 1] if name in args else None
if prog.endswith("synth.py"):
    rc = int(os.environ.get("FAKE_SYNTH_RC") or 0)
    if rc:
        sys.stderr.write(os.environ.get("FAKE_SYNTH_ERR", "") + "\n")
        sys.exit(rc)
    job = json.load(open(args[0]))
    d = opt("--per-turn")
    os.makedirs(d, exist_ok=True)
    total = 0.0
    for k, t in enumerate(job["turns"], 1):
        n = int(24000 * 0.05 * max(1, len(t["text"].split())))
        a = array.array("h", [1500, -1500]) * (n // 2)
        with wave.open(os.path.join(d, "turn-%04d.wav" % k), "wb") as w:
            w.setnchannels(1); w.setsampwidth(2); w.setframerate(24000); w.writeframes(a.tobytes())
        total += n / 24000.0
    th, wk = int(opt("--threads")), int(opt("--workers"))
    json.dump({"turns": len(job["turns"]), "turns_audio_s": round(total, 2), "threads": th,
               "workers": wk, "threads_per_worker": th // wk, "synth": "1.2",
               "startup_in_process_s": 1.0, "render_s": 2.0, "memory_gb": 24.0,
               "voices": job["voices"], "per_turn": [{"turn": 1}]}, open(opt("--timing"), "w"))
elif prog.endswith("transcribe.py"):
    audio, out, model = args[0], args[1], opt("--model")
    if not os.path.isfile(audio):
        sys.exit("no audio at " + audio)
    src = os.environ.get("FAKE_HEARD_" + model.replace(".", "_").replace("-", "_")) \
        or os.environ.get("FAKE_HEARD")
    text = open(src).read() if src else ""
    words, t = [], 0.0
    for w in text.split():
        words.append({"w": " " + w, "start": round(t, 2), "end": round(t + 0.3, 2), "p": 0.9})
        t += 0.35
    json.dump({"text": text, "words": words, "device": "cpu", "compute_type": "int8",
               "transcribe_s": 1.0, "asr_model": "whisper " + model}, open(out, "w"))
else:
    sys.exit("unknown program " + prog)
'''

FAKE_SHA = """#!/bin/sh
case "$1" in
  *kokoro*) s="${FAKE_KOKORO_SHA:-@KOKORO_PIN@}" ;;
  *whisper*) s="${FAKE_WHISPER_SHA:-@WHISPER_PIN@}" ;;
  *) s=0000 ;;
esac
echo "$s  $1"
"""


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)


def prepay_error():
    return ApiError(402, "RESOURCE_EXHAUSTED", PREPAY)


class FakeContainers(object):
    """A bin folder with the fake apptainer and sha256sum, first on PATH."""

    def start_containers(self):
        self._bin = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self._bin, True)
        _exe(os.path.join(self._bin, "apptainer"), FAKE_APPTAINER.replace("PYTHON", sys.executable))
        _exe(os.path.join(self._bin, "sha256sum"),
             FAKE_SHA.replace("@KOKORO_PIN@", pl.KOKORO_SHA256)
             .replace("@WHISPER_PIN@", pl.WHISPER_SHA256))
        self.calls = os.path.join(self._bin, "calls.jsonl")
        env = mock.patch.dict(os.environ, {"PATH": self._bin + os.pathsep + os.environ["PATH"],
                                           "FAKE_APPTAINER_LOG": self.calls})
        env.start()
        self.addCleanup(env.stop)

    def container_calls(self):
        if not os.path.isfile(self.calls):
            return []
        with open(self.calls) as f:
            return [json.loads(ln) for ln in f if ln.strip()]

    def kokoro_calls(self):
        return [c for c in self.container_calls() if any(a.endswith("synth.py") for a in c)]

    def heard(self, text, model=None):
        path = os.path.join(self._bin, f"heard_{model or 'all'}.txt")
        write(path, text)
        key = "FAKE_HEARD" + ("_" + model.replace(".", "_").replace("-", "_") if model else "")
        p = mock.patch.dict(os.environ, {key: path})
        p.start()
        self.addCleanup(p.stop)


class Local(FakeContainers, T.Workspace):
    """A checked script, Gemini available as a fake, and the local voice on a fake HIVE job."""

    def setUp(self):
        super().setUp()
        self.start_containers()
        self.env = mock.patch.dict(os.environ, {"GEMINI_API_KEY": FAKE_KEY})
        self.env.start()
        self.addCleanup(self.env.stop)
        for name, value in (("_sleep", lambda s: None),
                            ("encode_aac", lambda w, m: (None, "test"))):
            p = mock.patch.object(mp, name, value)
            p.start()
            self.addCleanup(p.stop)
        self.runner = pl.Here(JOB_ENV)
        self.reach = mock.patch.object(pl, "find_runner", return_value=self.runner)
        self.reach.start()
        self.addCleanup(self.reach.stop)
        self.assertEqual(self.check()[0], 0)
        self.heard(self.said())                       # Whisper hears what Kokoro was given

    def said(self, voices=None):
        s = mp.parse_script(self.script)
        speak = mp.local_speaker(s, voices)
        return " ".join(speak(t.text) for t in s.turns())

    def render(self, behaviour=None, *extra, preflight=False, tts="gemini"):
        mods, log = fake_genai(behaviour or (lambda m, c, k: pcm_for(c)))
        flags = () if preflight or tts != "gemini" else ("--no-preflight",)
        consent = ("--cloud-ok", "PI, 2026-10-20") if tts == "gemini" else ()
        with mock.patch.dict(sys.modules, mods):
            rc, out, err = run("render", self.script, "--tts", tts, *consent, *flags, *extra)
        return rc, out, err, log

    def manifest(self):
        return json.loads(read(os.path.join(self.pod, "podcast.json")))

    def cached(self):
        out = {}
        cache = os.path.join(self.pod, ".cache")
        for f in os.listdir(cache):
            if f.endswith(".json") and len(f) == 69:
                out[f[:-5]] = json.loads(read(os.path.join(cache, f))).get("model")
        return out


def stuck(release):
    def behaviour(model, contents, config):
        release.wait(30)
        return pcm_for(contents)
    return behaviour


class WhenItFallsBack(Local):
    """The four failures waiting cannot fix hand the whole episode to the local voice."""

    def assert_fell_back(self, rc, err, why):
        self.assertEqual(rc, 0, err)
        man = self.manifest()
        self.assertEqual(man["tts"]["backend"], "kokoro")
        self.assertEqual(man["tts"]["fallback"]["why"], why)
        self.assertEqual(man["tts"]["models_used"], [pl.KOKORO_LABEL])
        self.assertIn("making the WHOLE episode with the local voice", err)
        self.assertEqual(len(self.kokoro_calls()), 1, "all lines in ONE container call")
        return man

    def test_a_402_in_a_chunk(self):
        seen = {"n": 0}

        def behaviour(model, contents, config):
            seen["n"] += 1
            if seen["n"] > 1:
                raise prepay_error()
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", ANSWERS)
        man = self.assert_fell_back(rc, err, "billing")
        self.assertIn("prepay credits are depleted (402)", man["tts"]["fallback"]["reason"])
        self.assertEqual(man["tts"]["fallback"]["stopped_at"], "chunk 2/2 (segment 2)")
        self.assertEqual(man["key_tier"], "prepaid")              # what the 402 showed

    def test_a_402_in_the_preflight(self):
        rc, out, err, log = self.render(lambda m, c, k: (_ for _ in ()).throw(prepay_error()),
                                        "--model", ANSWERS, preflight=True)
        man = self.assert_fell_back(rc, err, "billing")
        self.assertEqual(man["tts"]["fallback"]["stopped_at"], "preflight")
        fb = man["tts"]["fallback"]
        self.assertEqual((fb["chunks_sent_to_gemini"], fb["gemini_chunks_sent"],
                          fb["preflight_requests"]), (0, 0, 1))
        # exactly what Google had: its sound check, none of this script
        self.assertIn("was sent 1 two-line sound check(s), generic sentences that are not this "
                      "script (1 in this render)", man["tts"]["what_was_sent"])
        self.assertNotIn("chunk(s) of this script", man["tts"]["what_was_sent"])
        self.assertIn("only a sound check, none of this transcript",
                      read(os.path.join(self.pod, "transcript.html")))

    def test_an_auth_failure(self):
        for where, pre in (("chunk", False), ("preflight", True)):
            with self.subTest(where=where):
                shutil.rmtree(os.path.join(self.pod, ".cache"), ignore_errors=True)
                os.remove(self.calls) if os.path.exists(self.calls) else None

                def behaviour(model, contents, config):
                    raise ApiError(403, "PERMISSION_DENIED", {"error": {
                        "code": 403, "status": "PERMISSION_DENIED",
                        "message": "API key not valid. Please pass a valid API key."}})
                rc, out, err, log = self.render(behaviour, "--model", ANSWERS, preflight=pre)
                self.assert_fell_back(rc, err, "auth")

    def test_no_quota_on_any_model(self):
        def behaviour(model, contents, config):
            raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(
                "GenerateRequestsPerDayPerProjectPerModel-FreeTier", value="0"))
        rc, out, err, log = self.render(behaviour)
        man = self.assert_fell_back(rc, err, "quota")
        self.assertGreater(len(man["tts"]["fallback"]["gemini_models"]), 1, "every model tried")

    def test_every_model_stalled(self):
        release = threading.Event()
        self.addCleanup(release.set)
        for name, value in (("HEARTBEAT", 0.05), ("PREFLIGHT_TIMEOUT", 0.3)):
            p = mock.patch.object(mp, name, value)
            p.start()
            self.addCleanup(p.stop)
        rc, out, err, log = self.render(stuck(release), "--model", STALLS, "--model", ANSWERS,
                                        "--request-timeout", "0.3")
        man = self.assert_fell_back(rc, err, "stalled")
        self.assertIn("stopped answering", man["tts"]["fallback"]["reason"])


class TimeoutBudget(Local):
    """Review: with the defaults (600 s, twice per model, a 30 min budget) the budget ran out on
    the second model's first timeout, so "every model stalled" mid-episode never handed over.
    The defaults' own ratio here (budget / timeout = 3), at 0.3 s a timeout."""
    TIMEOUT = 0.3
    BUDGET_MIN = TIMEOUT * (mp.REQUEST_TIMEOUT_BUDGET_MIN * 60 / mp.REQUEST_TIMEOUT) / 60

    def setUp(self):
        super().setUp()
        self.release = threading.Event()
        self.addCleanup(self.release.set)
        p = mock.patch.object(mp, "HEARTBEAT", 0.05)
        p.start()
        self.addCleanup(p.stop)

    def test_the_budget_after_a_stall_hand_over_is_every_model_stalled(self):
        self.assertEqual(mp.REQUEST_TIMEOUT_BUDGET_MIN * 60 / mp.REQUEST_TIMEOUT, 3.0)
        rc, out, err, log = self.render(stuck(self.release), "--model", STALLS, "--model",
                                        ANSWERS, "--request-timeout", str(self.TIMEOUT),
                                        "--request-timeout-budget", str(self.BUDGET_MIN))
        self.assertEqual(rc, 0, err)
        fb = self.manifest()["tts"]["fallback"]
        self.assertEqual(fb["why"], "stalled")
        self.assertIn("unanswered requests have cost", fb["message"])

    def test_the_budget_after_only_a_preflight_hand_over_is_the_resume_message(self):
        """Third review: a hand-over at the preflight is not a stall mid-episode."""
        self.render(None, "--model", STALLS)                       # cached by STALLS
        cache = os.path.join(self.pod, ".cache")
        os.remove(os.path.join(cache, sorted(f for f in os.listdir(cache) if f.endswith(".wav"))[0]))
        p = mock.patch.object(mp, "PREFLIGHT_TIMEOUT", 0.3)
        p.start()
        self.addCleanup(p.stop)

        def behaviour(model, contents, config):
            if model == STALLS or "short sound check" not in json.dumps(contents):
                self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", STALLS, "--model", ANSWERS,
                                        "--request-timeout", "0.3", "--request-timeout-budget",
                                        "0.004", preflight=True)
        self.assertEqual(rc, 1, err)
        self.assertIn("did not answer the preflight", err)
        self.assertIn("unanswered requests have cost", err)
        self.assertEqual(self.container_calls(), [])

    def test_the_budget_without_a_stall_is_the_resume_message(self):
        seen = {"n": 0}

        def slow_now_and_then(model, contents, config):
            seen["n"] += 1
            if seen["n"] % 2:                       # every other request goes unanswered once
                self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(slow_now_and_then, "--model", ANSWERS,
                                        "--request-timeout", "0.3",
                                        "--request-timeout-budget", "0.01")
        self.assertEqual(rc, 1)
        self.assertIn("unanswered requests have cost", err)
        self.assertEqual(self.container_calls(), [])


class WhenItWaits(Local):
    """Errors that waiting or another Gemini model clears never reach the local voice."""

    def test_a_rate_limit_with_a_retry_is_waited_out(self):
        seen = {"n": 0}

        def behaviour(model, contents, config):
            seen["n"] += 1
            if seen["n"] == 1:
                raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(
                    "GenerateRequestsPerMinutePerProjectPerModel-FreeTier", delay="3s"))
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        man = self.manifest()
        self.assertEqual(man["tts"]["backend"], "gemini")
        self.assertNotIn("fallback", man["tts"])
        self.assertEqual(self.container_calls(), [], "the local voice was never asked")

    def test_one_stalled_model_hands_over_to_the_next_gemini_model(self):
        release = threading.Event()
        self.addCleanup(release.set)
        p = mock.patch.object(mp, "HEARTBEAT", 0.05)
        p.start()
        self.addCleanup(p.stop)

        def behaviour(model, contents, config):
            if model == STALLS:
                release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", STALLS, "--model", ANSWERS,
                                        "--request-timeout", "0.3")
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.manifest()["tts"]["models_used"], [ANSWERS])
        self.assertEqual(self.container_calls(), [])

    def test_a_transient_error_stops_with_the_resume_message(self):
        def behaviour(model, contents, config):
            raise ApiError(503, "UNAVAILABLE", {"error": {"code": 503, "status": "UNAVAILABLE",
                                                          "message": "overloaded"}})
        rc, out, err, log = self.render(behaviour, "--model", ANSWERS)
        self.assertEqual(rc, 1)
        self.assertIn("re-run the same command to resume", err)
        self.assertEqual(self.container_calls(), [])

    def test_a_transient_error_on_every_model_at_the_preflight_waits(self):
        """Review: a 503 from every model at the two-line test fell back as "no quota"."""
        def behaviour(model, contents, config):
            raise ApiError(503, "UNAVAILABLE", {"error": {"code": 503, "status": "UNAVAILABLE",
                                                          "message": "The model is overloaded."}})
        rc, out, err, log = self.render(behaviour, preflight=True)
        self.assertEqual(rc, 1)
        self.assertIn("no Gemini TTS model answered a two-line test", err)
        self.assertIn("re-run later", err)
        self.assertNotIn("local voice", err)
        self.assertEqual(self.container_calls(), [])

    def assert_waits(self, rc, err):
        self.assertEqual(rc, 1, err)
        self.assertIn("no Gemini TTS model answered a two-line test", err)
        self.assertNotIn("local voice", err)
        self.assertEqual(self.container_calls(), [])

    def test_an_empty_answer_at_the_preflight_waits(self):
        rc, out, err, log = self.render(lambda m, c, k: b"", "--model", ANSWERS, preflight=True)
        self.assertIn("held no audio", err)
        self.assert_waits(rc, err)

    def test_a_reset_connection_at_the_preflight_waits(self):
        def reset(model, contents, config):
            raise ConnectionResetError("Connection reset by peer")
        rc, out, err, log = self.render(reset, preflight=True)
        self.assert_waits(rc, err)

    def test_one_transient_model_among_ones_with_no_quota_waits_too(self):
        def behaviour(model, contents, config):
            if model == ANSWERS:
                raise ApiError(503, "UNAVAILABLE", {"error": {"code": 503, "message": "busy"}})
            raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(
                "GenerateRequestsPerDayPerProjectPerModel-FreeTier", value="0"))
        rc, out, err, log = self.render(behaviour, preflight=True)
        self.assertEqual(rc, 1)
        self.assertEqual(self.container_calls(), [])

    def test_no_fallback_stops_as_before(self):
        rc, out, err, log = self.render(lambda m, c, k: (_ for _ in ()).throw(prepay_error()),
                                        "--model", ANSWERS, "--no-fallback")
        self.assertEqual(rc, 1)
        self.assertIn("prepay credits are depleted", err)
        self.assertNotIn("local voice", err)
        self.assertEqual(self.container_calls(), [])

    def test_error_info_kinds_behind_the_rule(self):
        self.assertEqual(mp.error_info(prepay_error())["kind"], "billing")
        self.assertEqual(mp.error_info(ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(
            "GenerateRequestsPerMinutePerProjectPerModel", delay="7s")))["kind"], "rate")
        self.assertEqual(mp.error_info(Exception("API key not valid"))["kind"], "fatal")


class EveryChunkAgain(Local):
    def test_gemini_chunks_are_never_in_the_kokoro_episode_and_stay_cached(self):
        seen = {"n": 0}

        def behaviour(model, contents, config):
            seen["n"] += 1
            if seen["n"] > 1:
                raise prepay_error()
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        man = self.manifest()
        self.assertEqual((man["chunks"], man["chunks_cached"], man["chunks_synthesized"]),
                         (2, 0, 2), "the chunk Gemini made is not reused")
        self.assertEqual(man["tts"]["fallback"]["chunks_sent_to_gemini"], 2)
        made_by = self.cached()
        self.assertEqual(sorted(set(made_by.values())), sorted({ANSWERS, pl.KOKORO_LABEL}),
                         "Gemini's chunk is kept for a later Gemini render")
        self.assertEqual(sum(1 for v in made_by.values() if v == pl.KOKORO_LABEL), 2)
        # every line of the script went to Kokoro, each once
        job = self.kokoro_calls()[0]
        self.assertIn("--per-turn", job)

    def test_what_was_sent_counts_chunks_from_an_earlier_render(self):
        """Review: a resume that met a 402 said Gemini "was sent no text of this script"."""
        seen = {"n": 0}

        def first(model, contents, config):
            seen["n"] += 1
            if seen["n"] > 1:
                raise ApiError(503, "UNAVAILABLE", {"error": {"code": 503, "message": "busy"}})
            return pcm_for(contents)
        rc, out, err, log = self.render(first, "--model", ANSWERS)
        self.assertEqual(rc, 1, "a transient error stops with the resume message")
        rc, out, err, log = self.render(lambda m, c, k: (_ for _ in ()).throw(prepay_error()),
                                        "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        fb = self.manifest()["tts"]["fallback"]
        # render 1 sent chunks 1 and 2 (chunk 2 failed); render 2 sent chunk 2 again
        self.assertEqual((fb["gemini_chunks_sent"], fb["chunks"], fb["chunks_sent_to_gemini"],
                          fb["gemini_chunks_cached"]), (2, 2, 1, 1))
        self.assertIn("the text of 2 of this script's 2 chunk(s), after pronunciation "
                      "substitutions (1 of them in this render)",
                      self.manifest()["tts"]["what_was_sent"])

    def test_chunks_whose_text_is_no_longer_in_the_script_are_counted(self):
        """Third review: a line edited after Gemini had it was not counted."""
        self.render(None, "--model", ANSWERS)                   # Gemini has both chunks
        self.assertEqual(self.check(T.read(self.script).replace("30 IPs, 4 baits", "4 baits, 30 IPs"))[0],
                         0)
        rc, out, err, log = self.render(lambda m, c, k: (_ for _ in ()).throw(prepay_error()),
                                        "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        fb = self.manifest()["tts"]["fallback"]
        self.assertEqual((fb["gemini_chunks_sent"], fb["gemini_chunks_not_in_script"]), (2, 1))
        self.assertIn("the text of 1 chunk(s) that are not in this script as it is now",
                      self.manifest()["tts"]["what_was_sent"])
        self.assertIn("it had been sent the text of 2 of the 2 parts of this transcript and the "
                      "text of 1 part(s) not in this transcript",
                      read(os.path.join(self.pod, "transcript.html")))

    def test_a_broken_record_never_stops_a_gemini_render(self):
        os.makedirs(os.path.join(self.pod, ".cache"), exist_ok=True)
        write(os.path.join(self.pod, ".cache", mp.SENT_FILE), json.dumps({"chunks": []}))
        rc, out, err, log = self.render(None, "--model", ANSWERS)
        self.assertEqual(rc, 0, err)

    def test_a_broken_record_never_stops_a_fallback_and_unknown_keys_stay(self):
        """Third review, round 3: a bad value inside the record broke the fallback itself."""
        os.makedirs(os.path.join(self.pod, ".cache"), exist_ok=True)
        path = os.path.join(self.pod, ".cache", mp.SENT_FILE)
        write(path, json.dumps({"preflight": {"requests": "x"}, "chunks": {"k": [1]},
                                "from_a_later_version": {"kept": True}}))
        rc, out, err, log = self.render(lambda m, c, k: (_ for _ in ()).throw(prepay_error()),
                                        "--model", ANSWERS, preflight=True)
        self.assertEqual(rc, 0, err)
        fb = self.manifest()["tts"]["fallback"]
        self.assertEqual((fb["preflight_requests"], fb["preflight_requests_before"]), (1, 0))
        got = json.loads(read(path))
        self.assertEqual(got["from_a_later_version"], {"kept": True})
        self.assertEqual(got["preflight"]["requests"], 1)

    def test_no_value_in_the_record_can_break_a_render(self):
        """Review of 4664ebd7: 1e400 read as inf (OverflowError); deep nesting (RecursionError)."""
        g = mp.GeminiTTS.__new__(mp.GeminiTTS)
        g.cache_dir = os.path.join(self.pod, ".cache")
        os.makedirs(g.cache_dir, exist_ok=True)
        path = os.path.join(g.cache_dir, mp.SENT_FILE)
        for raw in ('{"preflight": {"requests": 1e400}}', '{"preflight": {"requests": Infinity}}',
                    '{"preflight": {"requests": NaN}}', '{"preflight": {"requests": -3}}',
                    '{"preflight": {"requests": %s}}' % ("9" * 4000), "[" * 100000, "[1, 2]",
                    '{"chunks": {"k": null}, "preflight": []}'):
            with self.subTest(raw=raw[:40]):
                write(path, raw)
                rec = g.sent_record()
                self.assertTrue(0 <= rec["preflight"]["requests"] <= 10 ** 6)
                self.assertTrue(all(isinstance(v, dict) for v in rec["chunks"].values()))

    def test_the_page_names_everything_google_had(self):
        cases = {(2, 0, 0): "the text of 2 of the 2 parts of this transcript",
                 (2, 1, 0): "the text of 2 of the 2 parts of this transcript and the text of 1 "
                            "part(s) not in this transcript",
                 (2, 0, 1): "the text of 2 of the 2 parts of this transcript and a sound check",
                 (2, 1, 1): "the text of 2 of the 2 parts of this transcript, the text of 1 "
                            "part(s) not in this transcript and a sound check",
                 (0, 0, 1): "only a sound check, none of this transcript",
                 (0, 1, 0): "the text of 1 part(s) not in this transcript, none of this transcript",
                 (0, 1, 1): "the text of 1 part(s) not in this transcript and a sound check, none "
                            "of this transcript",
                 (0, 0, 0): "none of this transcript"}
        for (k, other, tests), want in cases.items():
            with self.subTest(k=k, other=other, tests=tests):
                fb = {"gemini_chunks_sent": k, "chunks": 2, "gemini_chunks_not_in_script": other,
                      "preflight_requests": tests, "preflight_requests_before": 0}
                self.assertEqual(mp.gemini_had(fb, page=True), want)

    def test_the_record_survives_the_cache_pruning(self):
        self.test_what_was_sent_counts_chunks_from_an_earlier_render()
        self.assertTrue(os.path.isfile(os.path.join(self.pod, ".cache", mp.SENT_FILE)))

    def test_a_later_gemini_render_resumes_from_its_own_chunk(self):
        self.test_gemini_chunks_are_never_in_the_kokoro_episode_and_stay_cached()
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        man = self.manifest()
        self.assertEqual(man["tts"]["backend"], "gemini")
        self.assertEqual((man["chunks_cached"], man["chunks_synthesized"]), (1, 1))


class Record(Local):
    def fall_back(self):
        rc, out, err, log = self.render(lambda m, c, k: (_ for _ in ()).throw(prepay_error()),
                                        "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        return out, err

    def test_podcast_json_says_what_made_it_and_why(self):
        out, err = self.fall_back()
        man = self.manifest()
        tts = man["tts"]
        self.assertEqual((tts["backend"], tts["model"], tts["sent_to_cloud"]),
                         ("kokoro", pl.KOKORO_LABEL, False))
        self.assertEqual(tts["consent"], "not needed for this audio, which was made offline on "
                                         "HIVE; Gemini, asked first, was sent the text of 1 of "
                                         "this script's 2 chunk(s), after pronunciation "
                                         "substitutions (1 of them in this render) under "
                                         "cloud_tts_consent")
        eng = tts["engine"]
        self.assertEqual((eng["image"], eng["image_sha256"], eng["synth"]),
                         (pl.KOKORO_IMAGE, pl.KOKORO_SHA256, "1.2"))
        self.assertEqual(eng["voices"], {"Maya": "bf_emma", "Leo": "af_heart"})
        self.assertEqual(eng["spoken_names"], {"Leo": "Lena"})
        self.assertEqual((eng["container_call"]["workers"], eng["container_call"]["threads"]),
                         (8, 16))
        fb = tts["fallback"]
        self.assertEqual((fb["from"], fb["to"], fb["why"]), ("gemini", "kokoro", "billing"))
        self.assertIn("Your prepayment credits are depleted", fb["google_said"])
        self.assertTrue(fb["at"])
        self.assertIn("the whole episode was made with Kokoro on HIVE instead", fb["note"])
        self.assertEqual(man["warnings"][0], fb["note"])
        # what went to Google before it failed, and the consent it went under
        self.assertEqual(man["cloud_tts_consent"], "PI, 2026-10-20")
        self.assertIn("was sent the text of 1 of this script's 2 chunk(s)", tts["what_was_sent"])
        self.assertNotIn("sound check", tts["what_was_sent"], "--no-preflight: none was sent")
        self.assertEqual([h["name"] for h in man["hosts"]], ["Maya", "Leo"])      # labels kept
        self.assertEqual([h["voice"] for h in man["hosts"]], ["bf_emma", "af_heart"])
        self.assertNotIn(FAKE_KEY, read(os.path.join(self.pod, "podcast.json")) + err)
        self.assertFalse(os.path.exists(self.calls + ".key"), "the Gemini key reached a container")
        self.assertIn("Nothing was sent to Google for the audio or its check.", err)
        self.assertEqual(json.loads(out)["fallback"], fb["note"])

    def test_the_check_is_whisper_on_hive(self):
        self.fall_back()
        v = self.manifest()["verify"]
        self.assertEqual(v["status"], "OK", v)
        self.assertEqual(v["audio_sent"], "nothing sent: transcribed on HIVE")
        self.assertFalse(v["cloud_consent"])
        self.assertEqual((v["image"], v["image_sha256"]), (pl.WHISPER_IMAGE, pl.WHISPER_SHA256))
        self.assertEqual(v["agreed_departures"], [])
        txt = read(os.path.join(self.pod, "verify.txt"))
        self.assertIn("Podcast audio check (ASR round trip, on HIVE): OK", txt)
        whisper = [c for c in self.container_calls() if any(a.endswith("transcribe.py") for a in c)]
        self.assertEqual(sorted(c[c.index("--model") + 1] for c in whisper),
                         sorted(pl.WHISPER_MODELS), "both transcriptions")

    def test_the_transcript_page_keeps_leo_and_names_the_spoken_name(self):
        self.fall_back()
        page = read(os.path.join(self.pod, "transcript.html"))
        self.assertIn("Leo: af_heart (said as Lena)", page)
        self.assertIn("Gemini, tried first, failed; it had been sent the text of 1 of the 2 "
                      "parts of this transcript", page)
        self.assertIn("<b>Leo</b>", page)
        self.assertNotIn("<b>Lena</b>", page)


class SpokenName(Local):
    LINES = (("MAYA", "Welcome to Signal to Noise. Quick note first: this is an AI-generated "
                      "discussion, and our voices are synthetic. Leo, you start."),
             ("LEO", "The report is the record. Leo's rule: 30 IPs, 4 baits."),
             ("MAYA", "Thanks, Leo. That is 6,112 proteins across 12 comparisons."))

    def setUp(self):
        super().setUp()
        self.assertEqual(self.check(T.script_text(segs=(self.LINES, T.SEG2)))[0], 0)

    def test_whole_word_in_the_spoken_text_only(self):
        s = mp.parse_script(self.script)
        speak = mp.local_speaker(s)
        lines = [speak(t.text) for t in s.turns()][:3]
        self.assertTrue(lines[0].endswith("Lena, you start."), lines[0])
        self.assertIn("Lena's rule", lines[1])
        self.assertIn("Thanks, Lena. That is", lines[2])
        self.assertEqual(speak("Not Leonardo, not LEO, not Leon: Leo."),
                         "Not Leonardo, not LEO, not Leon: Lena.")
        self.assertEqual([t.speaker for t in s.turns()][:3], ["Maya", "Leo", "Maya"])

    def test_gemini_says_leo(self):
        rc, out, err, log = self.render(None, "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        sent = json.dumps([c for _, c, _ in log["calls"]])
        self.assertIn("Leo, you start", sent)
        self.assertNotIn("Lena", sent)
        self.assertEqual(self.manifest()["hosts"][1]["voice"], "Charon")

    def test_a_mans_voice_keeps_the_name(self):
        s = mp.parse_script(self.script)
        self.assertEqual(pl.spoken_names({"Maya": "bf_emma", "Leo": "am_michael"}), {})
        speak = mp.local_speaker(s, {"Maya": "bf_emma", "Leo": "am_michael"})
        self.assertTrue(speak(s.turns()[0].text).endswith("Leo, you start."))

    def test_the_rendered_lines_say_lena(self):
        rc, out, err, log = self.render(None, tts="local")
        self.assertEqual(rc, 0, err)
        job = self.kokoro_calls()[0]
        d = job[job.index("--per-turn") + 1]
        self.assertTrue(d.endswith("turns"))
        self.assertIn("Thanks, Lena.", self.said())
        self.assertNotIn("Thanks, Leo.", self.said())


class Voices(unittest.TestCase):
    def test_the_voices_brett_chose(self):
        self.assertEqual(pl.VOICES, ("bf_emma", "af_heart"))
        self.assertEqual(pl.voices_for(["Maya", "Leo"]), {"Maya": "bf_emma", "Leo": "af_heart"})

    def test_a_header_can_name_others_and_an_unknown_voice_is_refused(self):
        self.assertEqual(pl.voices_for(["Maya", "Leo"], "Maya=af_bella, Leo=bm_george"),
                         {"Maya": "af_bella", "Leo": "bm_george"})
        with self.assertRaises(ValueError) as e:
            pl.voices_for(["Maya", "Leo"], "Leo=nobody")
        self.assertIn("nobody", str(e.exception))

    def test_gemini_voices_are_unchanged(self):
        self.assertEqual(mp.GEMINI_VOICES, ("Kore", "Charon"))


class TtsLocal(Local):
    def test_renders_with_no_consent_and_sends_nothing(self):
        mods = {"google": None, "google.genai": None, "google.genai.types": None}
        with mock.patch.dict(sys.modules, mods):          # importing google.genai would fail
            rc, out, err = run("render", self.script, "--tts", "local")
        self.assertEqual(rc, 0, err)
        self.assertNotIn("GEMINI_API_KEY", os.environ, "render --tts local passes no key on")
        man = self.manifest()
        self.assertEqual(man["tts"]["backend"], "kokoro")
        self.assertIsNone(man["tts"]["fallback"])
        self.assertEqual(man["tts"]["consent"], "not needed: made offline on HIVE, nothing was "
                                                "sent to any service")
        self.assertFalse(man["cloud_tts_consent"])
        self.assertEqual(man["tts"]["what_was_sent"], "nothing for this audio: made on HIVE with "
                                                      "Kokoro, and heard back there with Whisper")
        self.assertEqual(man["verify"]["status"], "OK")
        self.assertNotIn("Nothing was sent to Google", err)

    def test_a_second_render_makes_nothing_again(self):
        self.assertEqual(self.render(None, tts="local")[0], 0)
        self.assertEqual(len(self.kokoro_calls()), 1)
        rc, out, err, log = self.render(None, "--no-verify", tts="local")
        self.assertEqual(rc, 0, err)
        self.assertEqual(len(self.kokoro_calls()), 1, "every chunk was cached")
        self.assertEqual(self.manifest()["chunks_cached"], 2)

    def test_the_preview_shows_what_kokoro_reads(self):
        rc, out, err = run("render", self.script, "--tts", "local", "--dry-run")
        self.assertEqual(rc, 0, err)
        self.assertIn("J P H three", out)
        self.assertEqual(self.container_calls(), [])

    def test_verify_of_a_kokoro_episode_is_local_whatever_the_flags(self):
        self.assertEqual(self.render(None, "--no-verify", tts="local")[0], 0)
        rc, out, err = run("verify", self.script)            # no consent: none needed
        self.assertEqual(rc, 0, err)
        self.assertEqual(json.loads(out)["status"], "OK")

    def test_a_number_both_transcriptions_misheard_warns(self):
        self.heard(self.said().replace("6,112", "6,412"))
        rc, out, err, log = self.render(None, tts="local")
        self.assertEqual(rc, 0, err)
        v = self.manifest()["verify"]
        self.assertEqual(v["status"], "WARN")
        self.assertEqual([d["kind"] for d in v["agreed_departures"]], ["number"])
        self.assertIn("6112", v["agreed_departures"][0]["script"])

    def test_one_transcriber_alone_is_a_place_to_listen_not_a_warn(self):
        self.heard(self.said().replace("6,112", "6,412"), model="small.en")
        rc, out, err, log = self.render(None, tts="local")
        self.assertEqual(rc, 0, err)
        v = self.manifest()["verify"]
        self.assertEqual(v["status"], "OK")
        self.assertEqual([(d["kind"], d["by"]) for d in v["one_transcriber_only"]],
                         [("number", "small.en")])


class CheckFails(Local):
    def test_the_delivered_files_name_no_hive_folder(self):
        def fail(*a, **k):
            raise pl.LocalError("could not copy episode.m4a to HIVE: rsync: mkdir "
                                "\"/home/staff1/.proteomics-pipeline/podcast-local/ab12\" failed")
        with mock.patch.object(pl, "hear", fail):
            rc, out, err, log = self.render(None, tts="local")
        self.assertEqual(rc, 0, err)
        self.assertIn("/home/staff1/", err)                      # the log says why
        delivered = read(os.path.join(self.pod, "verify.txt")) + read(
            os.path.join(self.pod, "podcast.json"))
        self.assertNotIn("/home/", delivered)
        self.assertNotIn("staff1", delivered)
        self.assertEqual(self.manifest()["verify"]["status"], "NOT RUN")


class NoKey(Local):
    """Q2 (Brett, 2026-10-20): no Gemini key falls back like a refused key, said plainly."""

    def setUp(self):
        super().setUp()
        self.env.stop()                                          # no GEMINI_API_KEY
        for k in ("GEMINI_API_KEY",):
            p = mock.patch.dict(os.environ)
            p.start()
            self.addCleanup(p.stop)
            os.environ.pop(k, None)
        p = mock.patch.object(mp, "KEY_FILE", os.path.join(self.d, "no", "gemini_key"))
        p.start()
        self.addCleanup(p.stop)

    def render_no_key(self, *extra):
        mods = {"google": None, "google.genai": None, "google.genai.types": None}
        with mock.patch.dict(sys.modules, mods):              # never imported without a key
            return run("render", self.script, "--tts", "gemini", "--cloud-ok", "PI, 2026-10-20",
                       *extra)

    def test_no_key_makes_the_episode_with_the_local_voice(self):
        rc, out, err = self.render_no_key()
        self.assertEqual(rc, 0, err)
        self.assertIn("No Gemini key was set, which waiting cannot fix: making the WHOLE episode "
                      "with the local voice", err)
        self.assertNotIn("Gemini's chunks stay in the cache", err)
        man = self.manifest()
        fb = man["tts"]["fallback"]
        self.assertEqual((fb["why"], fb["reason"]), ("nokey", "No Gemini key was set"))
        self.assertEqual(man["tts"]["what_was_sent"],
                         "nothing for this audio: made on HIVE with Kokoro, and heard back there "
                         "with Whisper. Gemini was not asked: no Gemini key was set")
        self.assertEqual(man["tts"]["consent"], "not needed: made offline on HIVE, nothing was "
                                                "sent to any service")
        self.assertFalse(man["cloud_tts_consent"], "nothing was sent")
        self.assertNotIn("Gemini, tried first", read(os.path.join(self.pod, "transcript.html")))
        self.assertIn("No Gemini key was set; the whole episode was made with Kokoro", fb["note"])

    def test_a_key_file_found_but_not_used_says_why(self):
        """Review of 86671ee1: a key file refused for its mode was reported as no key."""
        why = f"{mp.KEY_FILE} is not used: other accounts can read it"   # find_key's own words
        write(mp.KEY_FILE, "a key that may not be used\n")              # it is there
        with mock.patch.object(mp, "find_key", return_value=(None, None, why)):
            rc, out, err = self.render_no_key()
        self.assertEqual(rc, 0, err)
        self.assertIn(f"No usable Gemini key ({why})", err)       # the log says why, in full
        fb = self.manifest()["tts"]["fallback"]
        # podcast.json is delivered: no key-file paths in it
        self.assertEqual(fb["reason"], "No usable Gemini key (a key file was found but not used; "
                                       "the render log says why)")
        self.assertNotIn(".config", json.dumps(self.manifest()))
        self.assertIn("Gemini was not asked: no usable Gemini key", self.manifest()["tts"]
                      ["what_was_sent"])

    def test_the_cores_key_for_another_run_is_a_note_not_a_refusal(self):
        """HIVE, 2026-10-20: a run that is not a Core run read "No usable Gemini key"."""
        note = "/quobyte/proteomics-grp/.config/gemini_key is used only for a Core run (--core-run)"
        with mock.patch.object(mp, "find_key", return_value=(None, None, note)):
            rc, out, err = self.render_no_key()
        self.assertEqual(rc, 0, err)
        self.assertIn(f"No Gemini key was set ({note})", err)
        write(mp.KEY_FILE, "")                       # an empty personal key file changes nothing
        with mock.patch.object(mp, "find_key", return_value=(None, None, note)):
            rc, out, err = self.render_no_key()
        self.assertIn(f"No Gemini key was set ({note})", err)
        self.assertEqual(self.manifest()["tts"]["fallback"]["reason"], "No Gemini key was set")
        self.assertNotIn("gemini_key", json.dumps(self.manifest()), "no key-file path delivered")

    def test_no_fallback_still_stops(self):
        rc, out, err = self.render_no_key("--no-fallback")
        self.assertEqual(rc, 1)
        self.assertIn("no Gemini key was set", err)
        self.assertEqual(self.container_calls(), [])

    def test_no_key_and_no_local_voice_says_both(self):
        self.reach.stop()
        self.reach = mock.patch.object(pl, "find_runner", side_effect=pl.Unavailable(
            "this computer has no HIVE login saved"))
        self.reach.start()
        rc, out, err = self.render_no_key()
        self.assertEqual(rc, 1)
        self.assertIn("No Gemini key was set, and the local voice could not make the episode "
                      "instead: this computer has no HIVE login saved", err)


class QuotaHandOver(Local):
    """Q3 (Brett, 2026-10-20): one model out of its daily quota mid-episode hands over to the
    next Gemini model, which makes every chunk again; only when every model is out of quota or
    stalled does the local voice take over."""
    DAILY = "GenerateRequestsPerDayPerProjectPerModel-FreeTier"

    def test_the_next_gemini_model_makes_every_chunk_again(self):
        seen = {"n": 0}

        def behaviour(model, contents, config):
            if model == STALLS:
                seen["n"] += 1
                if seen["n"] > 1:                                # chunk 2: out for the day
                    raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(self.DAILY))
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", STALLS, "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        man = self.manifest()
        self.assertEqual(man["tts"]["backend"], "gemini")
        self.assertEqual(man["tts"]["models_used"], [ANSWERS])
        self.assertEqual(man["chunks_synthesized"], 2)
        self.assertEqual([(x["from"], x["to"], x["kind"]) for x in man["tts"]["switched"]],
                         [(STALLS, ANSWERS, "quota")])
        self.assertIn("ran out of its daily quota", man["warnings"][0])
        self.assertEqual(self.container_calls(), [], "Gemini still had a model")

    def test_a_spent_budget_after_a_quota_hand_over_is_not_a_stall(self):
        release = threading.Event()
        self.addCleanup(release.set)
        p = mock.patch.object(mp, "HEARTBEAT", 0.05)
        p.start()
        self.addCleanup(p.stop)
        seen = {"n": 0}

        def behaviour(model, contents, config):
            if model == STALLS:
                seen["n"] += 1
                if seen["n"] > 1:
                    raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(self.DAILY))
                return pcm_for(contents)
            release.wait(30)                                     # the next model hangs once
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", STALLS, "--model", ANSWERS,
                                        "--request-timeout", "0.3", "--request-timeout-budget",
                                        "0.004")
        self.assertEqual(rc, 1, err)
        self.assertIn("unanswered requests have cost", err)
        self.assertEqual(self.container_calls(), [])

    def test_one_model_out_of_quota_falls_back_with_a_true_message(self):
        seen = {"n": 0}

        def behaviour(model, contents, config):
            seen["n"] += 1
            if seen["n"] > 1:
                raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(self.DAILY))
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        fb = self.manifest()["tts"]["fallback"]
        self.assertIn("no other Gemini model can be used", fb["message"])
        self.assertNotIn(f"resume with {ANSWERS}", fb["message"])

    def test_every_model_out_of_quota_falls_back(self):
        seen = {"n": 0}

        def behaviour(model, contents, config):
            seen["n"] += 1
            if seen["n"] > 1:
                raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(self.DAILY))
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", STALLS, "--model", ANSWERS)
        self.assertEqual(rc, 0, err)
        man = self.manifest()
        self.assertEqual((man["tts"]["backend"], man["tts"]["fallback"]["why"]), ("kokoro", "quota"))
        self.assertIn(f"{STALLS} ran out of its daily quota",
                      " ".join(man["tts"]["fallback"]["gemini_switched"]))


class ShortOutage(Local):
    """Q3 (Brett, 2026-10-20): a one-minute outage must not make a Kokoro episode. When every
    model misses the two-line test, the render waits PREFLIGHT_RETRY_WAIT (5 min) and asks them
    once more; only a second miss counts as "every model stalled". The wait is _sleep's, so
    these tests record it instead of waiting."""

    def setUp(self):
        super().setUp()
        self.release = threading.Event()
        self.addCleanup(self.release.set)                      # the hung worker threads end
        self.slept = []
        for name, value in (("HEARTBEAT", 0.05), ("PREFLIGHT_TIMEOUT", 0.3),
                            ("_sleep", self.slept.append)):
            p = mock.patch.object(mp, name, value)
            p.start()
            self.addCleanup(p.stop)

    def probes(self, log):
        return [m for m, c, _ in log["calls"] if mp.PREFLIGHT_LINES[0] in json.dumps(c)]

    def waits(self):
        return [x for x in self.slept if x == mp.PREFLIGHT_RETRY_WAIT]

    def test_an_outage_over_within_the_wait_renders_in_gemini(self):
        def behaviour(model, contents, config):
            if not self.waits():                               # Google is out until the wait
                self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", STALLS, "--model", ANSWERS,
                                        preflight=True)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.waits(), [300.0], "one 5-minute wait")
        self.assertIn(f"Waiting 5 min, then asking {STALLS}, {ANSWERS} once more, in case this "
                      "is a short outage at Google", err)
        self.assertEqual(self.probes(log), [STALLS, ANSWERS, STALLS])
        self.assertIn(f"asked again 5 min later: {STALLS} answered; using {STALLS}", err)
        man = self.manifest()
        self.assertEqual(man["tts"]["backend"], "gemini")
        self.assertNotIn("fallback", man["tts"])
        self.assertEqual(man["tts"]["models_used"], [STALLS])
        self.assertEqual((man["tts"]["preflight"]["model"], man["tts"]["preflight"]["waited_s"]),
                         (STALLS, 300.0))
        self.assertEqual(man["tts"]["preflight"]["skipped"][-1],
                         f"asked again 5 min later: {STALLS} answered")
        self.assertEqual(self.container_calls(), [])

    def test_a_second_miss_by_every_model_is_every_model_stalled(self):
        rc, out, err, log = self.render(stuck(self.release), "--model", STALLS, "--model",
                                        ANSWERS, preflight=True)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.waits(), [300.0], "one wait, then no more")
        self.assertEqual(self.probes(log), [STALLS, ANSWERS, STALLS, ANSWERS])
        man = self.manifest()
        fb = man["tts"]["fallback"]
        self.assertEqual((man["tts"]["backend"], fb["why"]), ("kokoro", "stalled"))
        self.assertIn("asked again 5 min later", fb["message"])
        self.assertEqual((fb["preflight_requests"], fb["preflight_waited_s"]), (4, 300.0))
        self.assertLess(err.index("Waiting 5 min"), err.index("making the WHOLE episode"))

    def test_the_wait_is_the_constant(self):
        p = mock.patch.object(mp, "PREFLIGHT_RETRY_WAIT", 7.0)
        p.start()
        self.addCleanup(p.stop)
        rc, out, err, log = self.render(stuck(self.release), "--model", ANSWERS, "--no-fallback",
                                        preflight=True)
        self.assertEqual(rc, 1, err)
        self.assertEqual(self.waits(), [7.0])
        self.assertIn("Waiting 7 s, then asking", err)
        self.assertIn("no Gemini TTS model answered a two-line test", err)
        self.assertEqual(self.container_calls(), [], "--no-fallback still stops")

    def test_no_wait_when_no_model_missed_the_test(self):
        """A 503 or a missing model is not an unanswered test: no wait, and the same outcome as
        before (the resume message; the local voice)."""
        for status, code, want in (("UNAVAILABLE", 503, 1), ("NOT_FOUND", 404, 0)):
            with self.subTest(status=status):
                del self.slept[:]
                shutil.rmtree(os.path.join(self.pod, ".cache"), ignore_errors=True)

                def behaviour(model, contents, config):
                    raise ApiError(code, status, {"error": {"code": code, "status": status,
                                                            "message": "no"}})
                rc, out, err, log = self.render(behaviour, "--model", ANSWERS, preflight=True)
                self.assertEqual(rc, want, err)
                self.assertEqual(self.waits(), [])
                self.assertNotIn("Waiting 5 min", err)
                self.assertNotIn("Traceback", err)
                if want:
                    self.assertIn("re-run later", err)
                    self.assertEqual(self.container_calls(), [])
                else:
                    self.assertEqual(self.manifest()["tts"]["fallback"]["why"], "quota")

    def test_a_free_tier_key_waits_too(self):
        """The shape a free-tier key meets: the pro model out of quota, both flash models
        missing the test. That fell back at once ("none"); now the flash models get the wait."""
        pro = "gemini-9-pro-tts"

        def behaviour(model, contents, config):
            if model == pro:
                raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(
                    "GenerateRequestsPerDayPerProjectPerModel-FreeTier", value="0"))
            if model == STALLS or not self.waits():
                self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", pro, "--model", STALLS, "--model",
                                        ANSWERS, preflight=True)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.waits(), [300.0])
        self.assertIn(f"then asking {STALLS}, {ANSWERS} once more", err)
        self.assertEqual(self.probes(log), [pro, STALLS, ANSWERS, STALLS, ANSWERS])
        man = self.manifest()
        self.assertEqual((man["tts"]["backend"], man["tts"]["models_used"]), ("gemini", [ANSWERS]))
        self.assertEqual(self.container_calls(), [])

    def test_a_model_used_after_the_wait_can_still_hand_over(self):
        """Every model that missed the first test has a fresh chance after the wait, asked again
        or not: one that answers then and stops mid-episode hands over to the next."""
        seen = {"n": 0}

        def behaviour(model, contents, config):
            if not self.waits():
                self.release.wait(30)
            elif model == STALLS and mp.PREFLIGHT_LINES[0] not in json.dumps(contents):
                seen["n"] += 1
                if seen["n"] >= 2:                               # chunk 2 never answers
                    self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", STALLS, "--model", ANSWERS,
                                        "--request-timeout", "0.3", preflight=True)
        self.assertEqual(rc, 0, err)
        man = self.manifest()
        self.assertEqual((man["tts"]["backend"], man["tts"]["models_used"]), ("gemini", [ANSWERS]))
        self.assertEqual([(x["from"], x["to"], x["kind"]) for x in man["tts"]["switched"]],
                         [(STALLS, ANSWERS, "stalled")])

    def test_a_remembered_model_that_answers_after_the_wait_is_forgotten(self):
        cache = os.path.join(self.pod, ".cache")
        os.makedirs(cache, exist_ok=True)
        with open(os.path.join(cache, mp.STALL_FILE), "w") as fh:
            json.dump({m: {"at": mp.now_iso(), "where": "chunk 2/2 (segment 2)",
                           "why": "2 requests of 600 s unanswered"} for m in (STALLS, ANSWERS)},
                      fh)

        def behaviour(model, contents, config):
            if "pro" in model:
                raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(
                    "GenerateRequestsPerDayPerProjectPerModel-FreeTier", value="0"))
            if model != ANSWERS or not self.waits():
                self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, preflight=True)     # no --model: remembered
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.waits(), [300.0])
        self.assertIn(f"preflight: {ANSWERS} stopped answering earlier; it answers now", err)
        self.assertEqual(err.count(f"preflight: {ANSWERS} stopped answering earlier ("), 1)
        with open(os.path.join(cache, mp.STALL_FILE)) as fh:
            self.assertEqual(sorted(json.load(fh)), [STALLS])
        self.assertEqual(self.manifest()["tts"]["models_used"], [ANSWERS])

    def remember(self, *models):
        cache = os.path.join(self.pod, ".cache")
        os.makedirs(cache, exist_ok=True)
        with open(os.path.join(cache, mp.STALL_FILE), "w") as fh:
            json.dump({m: {"at": mp.now_iso(), "where": "chunk 2/2 (segment 2)",
                           "why": "2 requests of 600 s unanswered"} for m in models}, fh)

    @staticmethod
    def pro_exhausted(model):
        if "pro" in model:
            raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(
                "GenerateRequestsPerDayPerProjectPerModel-FreeTier", value="0"))

    def test_a_remembered_model_not_asked_again_stays_out(self):
        """Second review: a model an earlier render saw stall, which missed this test and was
        not asked again (another answered first), is not handed the episode mid-way."""
        self.remember(STALLS, ANSWERS)
        state = {"used": None, "n": 0}

        def behaviour(model, contents, config):
            self.pro_exhausted(model)
            probe = mp.PREFLIGHT_LINES[0] in json.dumps(contents)
            if not self.waits() or model not in (STALLS, ANSWERS):
                self.release.wait(30)                            # the lite model never answers
            elif probe and state["used"] in (None, model):
                state["used"] = model
            elif model == state["used"]:
                state["n"] += 1
                if state["n"] >= 2:                              # chunk 2 never answers
                    self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--request-timeout", "0.3", preflight=True)
        other = ANSWERS if state["used"] == STALLS else STALLS
        self.assertEqual(rc, 0, err)
        self.assertIsNotNone(state["used"], err)
        self.assertEqual(self.probes(log).count(other), 1, "asked once, before the wait")
        self.assertNotIn(other, [m for m, c, _ in log["calls"]
                                 if mp.PREFLIGHT_LINES[0] not in json.dumps(c)])
        self.assertEqual(self.manifest()["tts"]["backend"], "kokoro")

    def test_the_fallback_message_keeps_its_advice(self):
        """Second review: the longest realistic shape (pro models out of free-tier quota, the
        flash models missing twice, remembered stalls) still ends with the advice."""
        self.remember(STALLS, ANSWERS)

        def behaviour(model, contents, config):
            self.pro_exhausted(model)
            self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, preflight=True)
        self.assertEqual(rc, 0, err)
        fb = self.manifest()["tts"]["fallback"]
        self.assertIn("asked again 5 min later", fb["message"])
        self.assertTrue(fb["message"].endswith("--tts say"), fb["message"][-200:])

    def test_a_resume_whose_cached_model_misses_twice_makes_every_chunk_again(self):
        rc, out, err, log = self.render(None, "--model", STALLS)
        self.assertEqual(rc, 0, err)
        cache = os.path.join(self.pod, ".cache")
        os.remove(os.path.join(cache, sorted(f for f in os.listdir(cache)
                                             if f.endswith(".wav"))[0]))

        def behaviour(model, contents, config):
            if model == STALLS or not self.waits():
                self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--model", STALLS, "--model", ANSWERS,
                                        preflight=True)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.waits(), [300.0])
        man = self.manifest()
        self.assertEqual(man["tts"]["models_used"], [ANSWERS])
        self.assertEqual((man["chunks_cached"], man["chunks_synthesized"]), (0, 2))
        self.assertEqual([(x["from"], x["to"], x["kind"]) for x in man["tts"]["switched"]],
                         [(STALLS, ANSWERS, "preflight")])
        self.assertIn("asked again 5 min later", man["tts"]["switched"][0]["why"])


class NotAvailable(Local):
    WHY = "this computer has no HIVE login saved (check_access.sh saves it)"

    def setUp(self):
        super().setUp()
        self.reach.stop()
        self.reach = mock.patch.object(pl, "find_runner", side_effect=pl.Unavailable(self.WHY))
        self.reach.start()

    def test_gemini_failed_and_the_local_voice_cannot_take_over(self):
        rc, out, err, log = self.render(lambda m, c, k: (_ for _ in ()).throw(prepay_error()),
                                        "--model", ANSWERS)
        self.assertEqual(rc, 1)
        self.assertIn("prepay credits are depleted", err)
        self.assertIn("and the local voice could not make the episode instead: " + self.WHY, err)
        self.assertIn("No episode this time", err)
        self.assertFalse(os.path.exists(os.path.join(self.pod, "podcast.json")))

    def test_tts_local_says_why(self):
        rc, out, err, log = self.render(None, tts="local")
        self.assertEqual(rc, 1)
        self.assertIn("the local voice (Kokoro on HIVE) is not available here: " + self.WHY, err)


class WhereItRuns(unittest.TestCase):
    def test_off_hive_without_a_hive_login(self):
        with mock.patch.object(pl, "_on_hive", return_value=False), \
                mock.patch("record_run.hive_login", return_value=None):
            with self.assertRaises(pl.Unavailable) as e:
                pl.find_runner()
        self.assertIn("no HIVE login saved", str(e.exception))

    def test_on_hive_the_images_must_be_readable(self):
        with tempfile.TemporaryDirectory() as d:
            k = os.path.join(d, "kokoro-1.2.sif")
            with mock.patch.object(pl, "_on_hive", return_value=True), \
                    mock.patch.object(pl, "KOKORO_IMAGE", k):
                with self.assertRaises(pl.Unavailable) as e:
                    pl.find_runner()
                self.assertIn("is not here", str(e.exception))
                write(k, "x")
                with mock.patch.object(pl.os, "access", return_value=False):
                    with self.assertRaises(pl.Unavailable) as e:
                        pl.find_runner()
                self.assertIn("members of proteomics-grp", str(e.exception))

    def test_on_a_hive_login_node_each_call_is_an_srun(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ("k.sif", "w.sif"):
                write(os.path.join(d, name), "x")
            with mock.patch.object(pl, "_on_hive", return_value=True), \
                    mock.patch.object(pl, "KOKORO_IMAGE", os.path.join(d, "k.sif")), \
                    mock.patch.object(pl, "WHISPER_IMAGE", os.path.join(d, "w.sif")), \
                    mock.patch("on_compute.on_login_node", return_value=True):
                r = pl.find_runner()
        self.assertIsInstance(r, pl.Srun)
        self.assertEqual((r.resources("speak"), r.resources("hear")), ((16, 24.0), (16, 12.0)))
        # never this login node's own /tmp, which the compute node cannot see (HIVE, 2026-10-20)
        with tempfile.TemporaryDirectory() as home:
            with mock.patch.dict(os.environ, {"PROTEOMICS_PIPELINE_HOME": home, "TMPDIR": "/tmp"}):
                d, _ = r.open()
                self.addCleanup(r.close, d)
                self.assertEqual(os.path.dirname(d), os.path.join(home, "podcast-local"))
                r.scratch = os.path.join(home, "cache")
                os.makedirs(r.scratch)
                d2, _ = r.open()
                self.assertEqual(os.path.dirname(d2), r.scratch)
                r.close(d2)
        with mock.patch("on_compute.srun_argv", return_value=["true"]) as srun:
            self.assertEqual(r.run("echo hi", "speak")[0], 0)
        argv, minutes, mem, cpus = srun.call_args[0][:4]
        self.assertEqual((argv[:2], minutes, mem, cpus), (["bash", "-c"], 30, "24G", 16))
        with mock.patch("on_compute.srun_argv", return_value=["true"]) as srun:
            r.run("echo hi", "hear")
        self.assertEqual(srun.call_args[0][1:4], (45, "12G", 16))

    def images(self, d):
        for name in ("k.sif", "w.sif"):
            write(os.path.join(d, name), "x")
        return (mock.patch.object(pl, "_on_hive", return_value=True),
                mock.patch.object(pl, "KOKORO_IMAGE", os.path.join(d, "k.sif")),
                mock.patch.object(pl, "WHISPER_IMAGE", os.path.join(d, "w.sif")),
                mock.patch("on_compute.on_login_node", return_value=False))

    def test_a_job_with_room_runs_it_itself_and_a_small_one_asks_for_its_own(self):
        """Review (ucd-hpc): the Gemini podcast job (1 CPU, 2 GB) must not hold Kokoro's 16 CPUs
        just in case; the local voice asks for them only when it takes over."""
        small = {"SLURM_JOB_ID": "7", "SLURM_CPUS_PER_TASK": "1", "SLURM_MEM_PER_NODE": "2048"}
        with tempfile.TemporaryDirectory() as d:
            a, b, c, e = self.images(d)
            with a, b, c, e:
                with mock.patch.dict(os.environ, JOB_ENV):
                    self.assertIsInstance(pl.find_runner(), pl.Here)
                with mock.patch.dict(os.environ, small):
                    r = pl.find_runner()
        self.assertIsInstance(r, pl.Srun)
        self.assertTrue(r.inside_job)
        self.assertIn("allocation of its own", r.where)

    def test_an_srun_from_a_job_asks_for_a_new_allocation_and_no_key(self):
        """srun inside a job makes a step of that job unless the job's SLURM variables are left
        out; and neither srun nor the container gets the Gemini key."""
        env = dict(JOB_ENV, SLURM_JOBID="1", SLURM_NODELIST="n1", PROTEOMICS_ON_COMPUTE="1",
                   GEMINI_API_KEY=FAKE_KEY, SLURM_CONF="/etc/slurm/slurm.conf")
        with mock.patch.dict(os.environ, env), \
                mock.patch("on_compute.srun_argv", return_value=["/usr/bin/env"]):
            rc, out, err = pl.Srun(inside_job=True).run("true", "speak")
        self.assertEqual(rc, 0)
        self.assertEqual([ln for ln in out.splitlines() if ln.startswith("SLURM")],
                         ["SLURM_CONF=/etc/slurm/slurm.conf"], "where SLURM's configuration is")
        self.assertNotIn("PROTEOMICS_ON_COMPUTE", out)
        self.assertNotIn(FAKE_KEY, out)
        self.assertIn("PATH=", out)

    def test_from_another_computer_hive_is_asked(self):
        def answer(said):
            return lambda *a, **k: (0, said, "")
        with mock.patch.object(pl, "_on_hive", return_value=False), \
                mock.patch("record_run.hive_login", return_value="staff1"):
            with mock.patch.object(pl.Remote, "_hx", answer("PODCAST_LOCAL_SKILL_OK\n")):
                with self.assertRaises(pl.Unavailable) as e:
                    pl.find_runner()
            self.assertIn("staff1 cannot read the Core's voice containers", str(e.exception))
            with mock.patch.object(pl.Remote, "_hx", answer("PODCAST_LOCAL_IMAGES_OK\n")):
                with self.assertRaises(pl.Unavailable) as e:
                    pl.find_runner()
            self.assertIn("the skill is not on HIVE", str(e.exception))
            self.assertIn("--put-skill", str(e.exception))
            with mock.patch.object(pl.Remote, "_hx", answer(        # 2.13.1: no on_compute.py
                    "PODCAST_LOCAL_IMAGES_OK\nPODCAST_LOCAL_SKILL_THERE\n")):
                with self.assertRaises(pl.Unavailable) as e:
                    pl.find_runner()
            self.assertIn("the skill on HIVE is older than this one", str(e.exception))
            self.assertIn("needs skill 2.13.2 or later", str(e.exception))
            with mock.patch.object(pl.Remote, "_hx", answer(
                    "PODCAST_LOCAL_IMAGES_OK\nPODCAST_LOCAL_SKILL_OK\n")):
                self.assertIsInstance(pl.find_runner(), pl.Remote)
            with mock.patch.object(pl.Remote, "_hx", lambda *a, **k: (255, "", "Connection refused")):
                with self.assertRaises(pl.Unavailable) as e:
                    pl.find_runner()
            self.assertIn("HIVE could not be asked (Connection refused)", str(e.exception))

            def hangs(*a, **k):                         # review: a traceback lost both reasons
                raise pl.LocalError("no answer in 2 min: stopped")
            with mock.patch.object(pl.Remote, "_hx", hangs):
                with self.assertRaises(pl.Unavailable) as e:
                    pl.find_runner()
            self.assertIn("HIVE did not answer (no answer in 2 min", str(e.exception))


class Plumbing(unittest.TestCase):
    def test_no_apptainer_is_said_plainly(self):
        empty = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, empty, True)
        r = subprocess.run(["/bin/bash", "-c", "set -eo pipefail\n" + pl._apptainer_setup()],
                           capture_output=True, text=True, env={"PATH": empty})
        self.assertEqual(r.returncode, 5)
        self.assertIn("no apptainer here", r.stderr)

    def test_a_timeout_stops_the_whole_process_group(self):
        class Slow(object):
            pid = 4242
            returncode = None
            calls = 0

            def communicate(self, timeout=None):
                Slow.calls += 1
                if Slow.calls == 1:
                    raise subprocess.TimeoutExpired("bash", timeout)
                return "", ""
        with mock.patch.object(pl.subprocess, "Popen", return_value=Slow()) as popen, \
                mock.patch.object(pl.os, "killpg") as killpg:
            with self.assertRaises(pl.LocalError) as e:
                pl._run(["bash", "-c", "sleep 99"], 1)
        self.assertTrue(popen.call_args.kwargs["start_new_session"])
        # SIGTERM first (an srun then cancels its allocation); SIGKILL to what is left
        self.assertEqual([c.args for c in killpg.call_args_list],
                         [(4242, pl.signal.SIGTERM), (4242, pl.signal.SIGKILL)])
        self.assertIn("no answer in 1 min", str(e.exception))

    def test_a_group_that_ignores_sigterm_is_killed(self):
        class Stubborn(object):
            pid, returncode, calls = 4343, None, 0

            def communicate(self, timeout=None):
                Stubborn.calls += 1
                if Stubborn.calls <= 2:                     # the run, then the SIGTERM wait
                    raise subprocess.TimeoutExpired("bash", timeout)
                return "", ""
        with mock.patch.object(pl.subprocess, "Popen", return_value=Stubborn()), \
                mock.patch.object(pl.os, "killpg") as killpg:
            with self.assertRaises(pl.LocalError):
                pl._run(["bash", "-c", "trap '' TERM; sleep 99"], 1)
        self.assertEqual([c.args[1] for c in killpg.call_args_list],
                         [pl.signal.SIGTERM, pl.signal.SIGKILL])

    def test_a_ctrl_c_during_the_sigterm_wait_still_kills_the_group(self):
        class Interrupted(object):
            pid, returncode, calls = 4444, None, 0

            def communicate(self, timeout=None):
                Interrupted.calls += 1
                if Interrupted.calls == 1:
                    raise subprocess.TimeoutExpired("bash", timeout)
                if Interrupted.calls == 2:                  # the SIGTERM wait
                    raise KeyboardInterrupt
                return "", ""
        with mock.patch.object(pl.subprocess, "Popen", return_value=Interrupted()), \
                mock.patch.object(pl.os, "killpg") as killpg:
            with self.assertRaises(KeyboardInterrupt):
                pl._run(["bash", "-c", "sleep 99"], 1)
        self.assertIn(pl.signal.SIGKILL, [c.args[1] for c in killpg.call_args_list])

    def test_the_gemini_key_is_not_handed_on(self):
        with mock.patch.dict(os.environ, {"GEMINI_API_KEY": FAKE_KEY}):
            rc, out, err = pl._run(["/bin/sh", "-c", 'echo "[${GEMINI_API_KEY:-none}]"'], 1)
        self.assertEqual(out.strip(), "[none]")


class RemoteCleanup(unittest.TestCase):
    def test_closing_never_hides_the_real_error_and_names_what_is_left(self):
        r = pl.Remote("staff1")
        with mock.patch.object(pl.Remote, "_hx", return_value=(0, "", "")):
            local, word = r.open()
        self.assertIn("$HOME/.proteomics-pipeline/podcast-local/", word)

        def gone(*a, **k):
            raise pl.LocalError("no answer in 2 min: stopped")
        err = __import__("io").StringIO()
        with mock.patch.object(pl.Remote, "_hx", gone), mock.patch.object(pl.sys, "stderr", err):
            r.close(local)                                       # no exception
        self.assertFalse(os.path.exists(local))
        self.assertIn("could not remove ~/.proteomics-pipeline/podcast-local/", err.getvalue())
        self.assertIn("remove it by hand", err.getvalue())

    def test_a_folder_that_cannot_be_made_leaves_nothing_here(self):
        made = []
        real = tempfile.mkdtemp

        def mkdtemp(*a, **k):
            made.append(real(*a, **k))
            return made[-1]

        def hangs(*a, **k):
            raise pl.LocalError("no answer in 2 min: stopped")
        with mock.patch.object(pl.tempfile, "mkdtemp", mkdtemp), \
                mock.patch.object(pl.Remote, "_hx", hangs), \
                mock.patch.object(pl.sys, "stderr", __import__("io").StringIO()):
            with self.assertRaises(pl.LocalError):
                pl.Remote("staff1").open()
        self.assertFalse(os.path.exists(made[0]))


class Sizing(unittest.TestCase):
    def test_workers_from_the_cpus_and_the_memory(self):
        for (cpus, mem), want in (((16, 24), (8, 16)), ((16, 32), (8, 16)), ((8, 24), (4, 8)),
                                  ((16, 10), (2, 4)), ((32, 64), (16, 32)), ((1, 2), (1, 1)),
                                  ((16, None), (8, 16)), ((2, 8), (1, 2))):
            with self.subTest(cpus=cpus, mem=mem):
                self.assertEqual(pl.layout(cpus, mem), want)

    def test_about_two_and_a_half_gb_per_worker(self):
        self.assertEqual((pl.GB_PER_WORKER, pl.THREADS_PER_WORKER), (2.5, 2))
        workers, _ = pl.layout(64, 100)
        self.assertLessEqual(workers * pl.GB_PER_WORKER + pl.RESERVE_GB, 100)

    def test_the_gemini_job_has_time_to_wait_for_a_fallback(self):
        """Gemini's own waits (the two-line test twice, 5 min apart, then the budgets), then
        each local-voice step with up to 15 min in the queue."""
        preflight = (2 * len(mp.KNOWN_TTS_MODELS) * mp.PREFLIGHT_TIMEOUT
                     + mp.PREFLIGHT_RETRY_WAIT) / 60
        need = (preflight + mp.RATE_WAIT_BUDGET_MIN + mp.REQUEST_TIMEOUT_BUDGET_MIN
                + sum(m + 15 for _, _, m in pl.STEPS.values()))
        self.assertLessEqual(need, mp.JOB_HOURS * 60)
        self.assertEqual((mp.JOB_CPUS, mp.JOB_MEM), (1, "2G"), "never Kokoro's CPUs just in case")

    def test_the_jobs_resources(self):
        self.assertEqual(pl.job_resources(JOB_ENV), (16, 24.0))
        self.assertEqual(pl.job_resources({"SLURM_JOB_ID": "1", "SLURM_CPUS_PER_TASK": "8",
                                           "SLURM_MEM_PER_CPU": "4096"}), (8, 32.0))
        self.assertEqual(pl.layout(pl.JOB_CPUS, pl._gb(pl.JOB_MEM)), (8, 16),
                         "the podcast job holds 8 workers")

    def test_an_srun_step_says_its_job_id_and_cleans_up_the_note(self):
        """Q1: the podcast job (1 CPU) asks for the step's own job; its id is said in the log
        and kept beside the episode's cache while it runs, for whoever cancels the podcast."""
        with tempfile.TemporaryDirectory() as d:
            r = pl.Srun(inside_job=True)
            r.scratch = d
            said, seen = [], {}

            def srun(cmd, *a, **k):
                return ["env", "SLURM_JOB_ID=4242", "bash", "-c",
                        cmd[2] + "\nsleep 3\ncat " + shlex.quote(os.path.join(d, pl.JOB_FILE))]
            with mock.patch("on_compute.srun_argv", side_effect=srun):
                rc, out, err = r.run("true", "speak", said.append)
            self.assertEqual(rc, 0, err)
            self.assertTrue(out.startswith("4242 speak "), out)
            self.assertTrue(any("runs as HIVE job 4242" in m and "`scancel 4242`" in m
                                for m in said), said)
            self.assertFalse(os.path.exists(os.path.join(d, pl.JOB_FILE)), "removed after")

    def test_a_note_left_by_a_cancelled_podcast_job_is_named_then_replaced(self):
        """Review of 86671ee1: a cancelled podcast job leaves its step's note behind."""
        with tempfile.TemporaryDirectory() as d:
            write(os.path.join(d, pl.JOB_FILE), "9999 speak 2026-10-20T20:00:00Z\n")
            r = pl.Srun(inside_job=True)
            r.scratch = d
            said = []

            def srun(cmd, *a, **k):
                return ["env", "SLURM_JOB_ID=4242", "bash", "-c", cmd[2] + "\nsleep 3"]
            with mock.patch("on_compute.srun_argv", side_effect=srun):
                rc, out, err = r.run("true", "hear", said.append)
            self.assertEqual(rc, 0, err)
            self.assertTrue(any("an earlier local-voice job, 9999" in m and "`scancel 9999`" in m
                                for m in said), said)
            self.assertTrue(any(m.startswith("[verify]") and "HIVE job 4242" in m for m in said),
                            said)
            self.assertFalse(os.path.exists(os.path.join(d, pl.JOB_FILE)))

    def test_srun_refusals_are_read_as_refusals(self):
        for said in ("srun: error: Unable to allocate resources: Invalid qos specification",
                     "srun: error: Job submit/allocate failed: Requested time limit is invalid",
                     "srun: Job allocation 4242 has been revoked",
                     "srun: error: Unable to allocate resources: QOSMaxCpuPerJobLimit"):
            with self.subTest(said):
                self.assertTrue(pl._REFUSED.search(said))
        self.assertFalse(pl._REFUSED.search("srun: error: n1: task 0: Exited with exit code 4"))

    def test_a_job_that_started_and_failed_is_not_called_refused(self):
        with tempfile.TemporaryDirectory() as d:
            r = pl.Srun(inside_job=True)
            r.scratch = d
            for name, argv in (
                    ("wrote its note", lambda cmd: ["env", "SLURM_JOB_ID=4242", "bash", "-c",
                                                    cmd[2] + "\nexit 5"]),
                    ("a task's exit", lambda cmd: ["bash", "-c", "echo 'srun: error: n1: task 0: "
                                                   "Exited with exit code 4' >&2; exit 4"])):
                with self.subTest(name), \
                        mock.patch("on_compute.srun_argv", side_effect=lambda c, *a, **k: argv(c)):
                    rc, out, err = r.run("true", "speak")          # no LocalError
                    self.assertNotEqual(rc, 0)

    def test_a_job_hive_refuses_is_said_plainly(self):
        with tempfile.TemporaryDirectory() as d:
            r = pl.Srun(inside_job=True)
            r.scratch = d
            refused = ["bash", "-c", "echo 'srun: error: Unable to allocate resources: Invalid "
                                     "qos specification' >&2; exit 1"]
            with mock.patch("on_compute.srun_argv", return_value=refused):
                with self.assertRaises(pl.LocalError) as e:
                    r.run("true", "speak")
        self.assertIn("HIVE did not give the local voice a job (srun: error: Unable to allocate "
                      "resources: Invalid qos specification)", str(e.exception))

    def test_work_folders_apptainer_can_bind(self):
        """Review: a , or : in a path breaks apptainer's --bind (it splits on them)."""
        self.assertFalse(pl._bindable("/quobyte/x/PROT 0807, v2/output/podcast/.cache"))
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "plain"))
            os.symlink(os.path.join(d, "plain"), os.path.join(d, "a:b"))
            self.assertFalse(pl._bindable(os.path.join(d, "a:b")), "the path as given counts")
        self.assertTrue(pl._bindable("/quobyte/x/Core service/PROT 0807 (v2)/.cache"))
        with tempfile.TemporaryDirectory() as home:
            bad = os.path.join(home, "a,b")
            os.makedirs(bad)
            r = pl.Srun()
            r.scratch = bad
            with mock.patch.dict(os.environ, {"PROTEOMICS_PIPELINE_HOME": home}):
                d, word = r.open()
            self.assertEqual(os.path.dirname(d), os.path.join(home, "podcast-local"))
            r.close(d)
            with mock.patch.dict(os.environ, {"TMPDIR": bad}):
                d, _ = pl.Here(JOB_ENV).open()
            self.assertEqual(os.path.dirname(d), "/tmp")
            pl.Here(JOB_ENV).close(d)
            r.scratch = bad
            with mock.patch.dict(os.environ, {"PROTEOMICS_PIPELINE_HOME": bad}):
                with self.assertRaises(pl.LocalError) as e:
                    r.open()
            self.assertIn("hold a , or :", str(e.exception))


class ContainerCall(FakeContainers, unittest.TestCase):
    """podcast_local.speak/hear: the real job scripts, a fake container."""

    def setUp(self):
        self.start_containers()
        self.runner = pl.Here(JOB_ENV)

    def test_speak_sizes_the_call_and_returns_each_line(self):
        d, paths, timing = pl.speak(self.runner, {"A": "bf_emma", "B": "af_heart"},
                                    [("A", "one two three"), ("B", "four")])
        self.addCleanup(self.runner.close, d)
        self.assertEqual([os.path.basename(p) for p in paths], ["turn-0001.wav", "turn-0002.wav"])
        call = self.kokoro_calls()[0]
        self.assertIn("--cleanenv", call)
        self.assertEqual(call[call.index("--no-mount") + 1], "home,cwd", "the key file is in home")
        self.assertEqual(call[call.index("--workers") + 1], "8")
        self.assertEqual(call[call.index("--threads") + 1], "16")
        self.assertIn(pl.KOKORO_IMAGE, call)
        self.assertIn("wall_s", timing)
        self.assertNotIn("per_turn", timing)

    def test_an_image_that_is_not_the_pinned_one_is_refused(self):
        with mock.patch.dict(os.environ, {"FAKE_KOKORO_SHA": "f" * 64}):
            with self.assertRaises(pl.LocalError) as e:
                pl.speak(self.runner, {"A": "bf_emma"}, [("A", "hello")])
        self.assertIn("is not the pinned image", str(e.exception))
        self.assertIn(pl.KOKORO_SHA256, str(e.exception))
        self.assertEqual(self.kokoro_calls(), [], "never run")

    def test_a_dead_worker_fails_at_once_with_its_reason(self):
        said = ("[synth 1.2] a Kokoro worker process ended abruptly (A process in the process "
                "pool was terminated abruptly). Most often it was killed for using more memory")
        with mock.patch.dict(os.environ, {"FAKE_SYNTH_RC": "3", "FAKE_SYNTH_ERR": said}):
            with self.assertRaises(pl.LocalError) as e:
                pl.speak(self.runner, {"A": "bf_emma"}, [("A", "hello")])
        self.assertIn("a worker process died", str(e.exception))
        self.assertIn("killed for using more memory", str(e.exception))
        with mock.patch.dict(os.environ, {"FAKE_SYNTH_RC": "137"}):
            with self.assertRaises(pl.LocalError) as e:
                pl.speak(self.runner, {"A": "bf_emma"}, [("A", "hello")])
        self.assertIn("most often for using more memory than the job has", str(e.exception))

    def test_hear_asks_both_models_in_one_call(self):
        self.heard("hello there")
        with tempfile.TemporaryDirectory() as d:
            audio = os.path.join(d, "podcast.wav")
            write(audio, "RIFF")
            got = pl.hear(self.runner, audio)
        self.assertEqual(sorted(got), sorted(pl.WHISPER_MODELS))
        self.assertEqual(got["small.en"]["text"], "hello there")
        calls = [c for c in self.container_calls() if any(a.endswith("transcribe.py") for a in c)]
        self.assertEqual([c[c.index("--threads") + 1] for c in calls], ["8", "8"])
        for c in calls:          # review: never the episode's folder, read-only, over the work one
            self.assertIn("--cleanenv", c)
            self.assertEqual(c[c.index("--no-mount") + 1], "home,cwd")
            binds = [c[i + 1] for i, a in enumerate(c) if a == "--bind"]
            self.assertEqual(len(binds), 1, c)
            audio = c[c.index("/opt/whisper/transcribe.py") + 1]
            self.assertEqual(os.path.dirname(audio), binds[0])
            self.assertEqual(os.path.basename(audio), "episode.wav")

    def test_cleanup_on_any_exception(self):
        class Stop(pl.Here):
            def run(self, script, step, log=None):
                raise KeyboardInterrupt
        r = Stop(JOB_ENV)
        made = []
        real_open = r.open

        def opened():
            d = real_open()
            made.append(d[0])
            return d
        r.open = opened
        with self.assertRaises(KeyboardInterrupt):
            pl.speak(r, {"A": "bf_emma"}, [("A", "hello")])
        self.assertFalse(os.path.exists(made[0]), "the folder with the script's lines is gone")


class Agreement(unittest.TestCase):
    def asr(self, text):
        return {"text": text, "words": [{"w": " " + w, "start": float(i)}
                                        for i, w in enumerate(text.split())]}

    def test_a_number_both_misheard_is_agreed(self):
        said = "We found 6,112 proteins in 12 comparisons."
        heard = self.asr("We found 6,412 proteins in 12 comparisons.")
        agreed, one = pl.compare_two(mp, said, heard, heard)
        self.assertEqual([d["kind"] for d in agreed], ["number"])
        self.assertEqual(one, [])

    def test_numbers_do_not_join_across_a_full_stop(self):
        said = "It went up to thirty. Three cultures grew."
        heard = self.asr("It went up to 30. 3 cultures grew.")
        agreed, one = pl.compare_two(mp, said, heard, heard)
        self.assertEqual((agreed, one), ([], []))

    def test_a_leading_zero_still_counts(self):
        self.assertEqual(pl.same_value("4.10"), pl.same_value("4.1"))
        self.assertNotEqual(pl.same_value("0807"), pl.same_value("807"))


class Job(unittest.TestCase):
    """make_podcast.py job --tts local: no consent, key or google-genai; render --tts local."""

    def test_the_local_job(self):
        import test_make_podcast_hive as H
        t = H.Job("test_the_job_script")
        t.setUp()
        self.addCleanup(t.doCleanups)
        self.addCleanup(t.tearDown)
        os.remove(t.core_file)                                  # no key at all
        with mock.patch.dict(os.environ, {"FAKE_GENAI_VERSION": "none"}):
            rc, out, err = t.job("--tts", "local")
        self.assertEqual(rc, 0, err)
        got = json.loads(out)
        self.assertEqual(got["consent"], "not needed: the local voice sends nothing")
        self.assertIsNone(got["key"])
        self.assertEqual(got["resources"], {"cpus": 16, "mem": "24G", "hours": 3})
        text = read(got["job"])
        self.assertIn(f"render {t.script} --tts local || exit $?", text)
        self.assertNotIn("--core-run", text)
        self.assertNotIn("--cloud-ok", text)
        self.assertIn("#SBATCH --cpus-per-task=16", text)
        self.assertEqual(subprocess.run(["bash", "-n", got["job"]]).returncode, 0)

    def test_a_gemini_job_without_a_key_is_written_for_the_local_voice(self):
        import test_make_podcast_hive as H
        t = H.Job("test_the_job_script")
        t.setUp()
        self.addCleanup(t.doCleanups)
        self.addCleanup(t.tearDown)
        os.remove(t.core_file)
        rc, out, err = t.job("--core-run")
        self.assertEqual(rc, 0, err)
        self.assertIn("the job's render will make the episode with the local voice", err)
        # a python with google-genai is still preferred: a key may be saved before it runs
        self.assertEqual(json.loads(out)["google_genai"], "2.31.0")


class Documented(unittest.TestCase):
    def read(self, *parts):
        with open(os.path.join(os.path.dirname(HERE), *parts), encoding="utf-8") as fh:
            return " ".join(fh.read().split())

    def test_the_reference_describes_the_local_voice(self):
        t = self.read("references", "podcast.md")
        missing = [s for s in ("--tts local", "--no-fallback", "bf_emma", "af_heart", "Lena",
                               "tts.fallback", "kokoro-1.2.sif", pl.KOKORO_SHA256,
                               "whisper-1.0.sif", "verify --local", "every chunk", "16 CPUs")
                   if s not in t]
        self.assertEqual(missing, [], "references/podcast.md")

    def test_skill_md_points_to_it(self):
        self.assertTrue("--tts local" in self.read("SKILL.md"), "SKILL.md names --tts local")


if __name__ == "__main__":
    unittest.main(verbosity=2)
