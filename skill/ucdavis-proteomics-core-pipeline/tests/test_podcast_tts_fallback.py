#!/usr/bin/env python3
"""
A podcast render is never pinned to a Gemini TTS model that stopped answering
(make_podcast.py: GeminiTTS.stall_switch / preflight / gemini_client / note_tier).

A staff render (PROT_0000, skill 2.11.2, 2026-10-16) sat 28 minutes on a stalled
gemini-3.8-flash-tts while a two-line two-speaker request to gemini-2.5-flash-preview-tts answered
in 2.9 s. The episode's first chunks were cached from the stalled model, and a model was handed
over only before the first chunk and only when missing / out of quota / rejected, so the render
could never fall back. The staff member probed each model by hand and re-rendered every chunk with
--model gemini-2.5-flash-preview-tts. Its consent note said "paid"; the key was free tier.

No network: google.genai is test_make_podcast.fake_genai (the HTTP-options case adds the two
option types and a Client that takes http_options, as google-genai 2.31.0's do).
"""
import json
import os
import shutil
import sys
import threading
import types
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from test_make_podcast import (FAKE_KEY, ApiError, Workspace, fake_genai, pcm_for,  # noqa: E402
                               read, rpc_error, run, mp)

STALLS = "gemini-3.8-flash-tts"
ANSWERS = "gemini-2.5-flash-preview-tts"
FREE_DAY = "GenerateRequestsPerDayPerProjectPerModel-FreeTier"


class Gemini(Workspace):
    """test_make_podcast.Gemini's setting, without its tests (a subclass would run them again)."""

    def setUp(self):
        super().setUp()
        self.env = mock.patch.dict(os.environ, {"GEMINI_API_KEY": FAKE_KEY})
        self.env.start()
        self.nosleep = mock.patch.object(mp, "_sleep", lambda s: None)
        self.nosleep.start()
        self.enc = mock.patch.object(mp, "encode_aac", return_value=(None, "test"))
        self.enc.start()
        self.assertEqual(self.check()[0], 0)

    def tearDown(self):
        for p in (self.env, self.nosleep, self.enc):
            p.stop()
        super().tearDown()

    def render(self, behaviour, *extra, preflight=False, **kw):
        mods, log = fake_genai(behaviour, **kw)
        flags = ("--no-verify",) + (() if preflight else ("--no-preflight",))
        with mock.patch.dict(sys.modules, mods):
            rc, out, err = run("render", self.script, "--tts", "gemini", *flags, *extra)
        return rc, out, err, log


def is_probe(contents):
    text = json.dumps(contents) if not isinstance(contents, str) else contents
    return mp.PREFLIGHT_LINES[0] in text


class StallBase(Gemini):
    def setUp(self):
        super().setUp()
        self.release = threading.Event()
        self.addCleanup(self.release.set)                      # the hung worker threads end
        for name, value in (("HEARTBEAT", 0.05), ("PREFLIGHT_TIMEOUT", 0.3)):
            p = mock.patch.object(mp, name, value)
            p.start()
            self.addCleanup(p.stop)

    def manifest(self):
        return json.loads(read(os.path.join(self.pod, "podcast.json")))

    def wavs(self):
        return [f for f in os.listdir(os.path.join(self.pod, ".cache")) if f.endswith(".wav")]


class Fallback(StallBase):

    def test_a_model_that_stops_answering_hands_over_and_every_chunk_is_made_again(self):
        seen = {"n": 0}

        def behaviour(model, contents, config):
            if model == STALLS:
                seen["n"] += 1
                if seen["n"] >= 2:                          # chunk 2 never answers
                    self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--cloud-ok", "PI, 2026-10-16",
                                        "--model", STALLS, "--model", ANSWERS,
                                        "--request-timeout", "0.4")
        self.assertEqual(rc, 0, err)
        self.assertIn(f"{STALLS} stopped answering (chunk 2/2", err)
        self.assertIn(f"Switching to {ANSWERS} and making EVERY chunk again with it", err)
        man = self.manifest()
        self.assertEqual(man["tts"]["models_used"], [ANSWERS])             # one model's voices
        self.assertEqual(man["chunks_synthesized"], 2)
        self.assertEqual([(x["from"], x["to"]) for x in man["tts"]["switched"]],
                         [(STALLS, ANSWERS)])
        self.assertTrue(any("stopped answering" in w for w in man["warnings"]))
        # the requests: chunk 1 on the stalled model, chunk 2 twice unanswered, then both on the next
        self.assertEqual([m for m, _, _ in log["calls"]], [STALLS] * 3 + [ANSWERS] * 2)

    def test_with_no_other_model_it_stops_with_the_resume_message(self):
        def behaviour(model, contents, config):
            self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--cloud-ok", "--model", STALLS,
                                        "--request-timeout", "0.3")
        self.assertEqual(rc, 1)
        self.assertIn("did not answer within 0 s, 2 times", err)
        self.assertIn("re-run the same command to resume", err)

    def test_preflight_passes_over_a_model_that_does_not_answer(self):
        def behaviour(model, contents, config):
            if model == STALLS:
                self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--cloud-ok", "--model", STALLS,
                                        "--model", ANSWERS, preflight=True)
        self.assertEqual(rc, 0, err)
        self.assertIn(f"preflight: {STALLS}: no answer in 0 s; using {ANSWERS}", err)
        probes = [m for m, c, _ in log["calls"] if is_probe(c)]
        self.assertEqual(probes, [STALLS, ANSWERS])
        self.assertEqual([m for m, c, _ in log["calls"] if not is_probe(c)], [ANSWERS, ANSWERS])
        man = self.manifest()
        self.assertEqual(man["tts"]["models_used"], [ANSWERS])
        self.assertTrue(man["tts"]["preflight"]["run"])
        self.assertEqual(man["tts"]["preflight"]["model"], ANSWERS)

    def test_a_resume_whose_cached_model_stalls_makes_every_chunk_again(self):
        """The staff render's case: chunks cached from a model that then stalls."""
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--cloud-ok", "--model",
                                        STALLS)
        self.assertEqual(rc, 0, err)
        cache = os.path.join(self.pod, ".cache")
        os.remove(os.path.join(cache, sorted(self.wavs())[0]))       # one chunk left to make

        def behaviour(model, contents, config):
            if model == STALLS:
                self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--cloud-ok", "--model", STALLS, "--model",
                                        ANSWERS, preflight=True)
        self.assertEqual(rc, 0, err)
        self.assertIn(f"{STALLS} made this episode's cached chunks but did not answer the "
                      f"preflight; {ANSWERS} makes every chunk again", err)
        man = self.manifest()
        self.assertEqual(man["tts"]["models_used"], [ANSWERS])
        self.assertEqual((man["chunks_cached"], man["chunks_synthesized"]), (0, 2))

    def test_no_model_answering_the_preflight_stops_before_any_chunk(self):
        def behaviour(model, contents, config):
            self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--cloud-ok", "--model", STALLS, "--model",
                                        ANSWERS, preflight=True)
        self.assertEqual(rc, 1)
        self.assertIn("no Gemini TTS model answered a two-line test", err)
        self.assertEqual([c for _, c, _ in log["calls"] if not is_probe(c)], [])

    def test_a_rate_limited_model_counts_as_answering(self):
        state = {"n": 0}

        def behaviour(model, contents, config):
            state["n"] += 1
            if state["n"] == 1:
                raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(
                    "GenerateRequestsPerMinutePerProjectPerModel", delay="2s"))
            return pcm_for(contents)
        rc, out, err, log = self.render(behaviour, "--cloud-ok", "--model", ANSWERS,
                                        preflight=True)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.manifest()["tts"]["preflight"]["model"], ANSWERS)

    def test_no_preflight_is_recorded(self):
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--cloud-ok", "--model",
                                        ANSWERS)
        self.assertEqual(rc, 0, err)
        self.assertEqual(self.manifest()["tts"]["preflight"],
                         {"run": False, "why": "--no-preflight"})
        self.assertFalse(any(is_probe(c) for _, c, _ in log["calls"]))


class ResumeAfterAStall(StallBase):
    """Review of 2.13.2 (H3), the staff incident's own shape, WITHOUT --model on the resume (so
    prepare pins the cached model): chunks cached from a model that then stalls, pro models with
    no free-tier quota, and the 2.5 flash model that answers."""
    PRO = "gemini-9-pro-tts"

    def first_render_stops_on_a_stall(self):
        n = {"k": 0}

        def first(model, contents, config):
            n["k"] += 1
            if n["k"] >= 2:
                self.release.wait(30)
            return pcm_for(contents)
        rc, out, err, log = self.render(first, "--cloud-ok", "--model", STALLS,
                                        "--request-timeout", "0.3")
        self.assertEqual(rc, 1, err)
        self.assertEqual(len(self.wavs()), 1)

    def second(self, model, contents, config):
        if "pro" in model:
            raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(FREE_DAY, value="0"))
        if model == STALLS and not is_probe(contents):
            self.release.wait(30)
        return pcm_for(contents)

    def resume(self):
        return self.render(self.second, "--cloud-ok", "--request-timeout", "0.3", preflight=True)

    def test_a_stall_is_remembered_so_the_resume_never_pins_that_model(self):
        self.first_render_stops_on_a_stall()
        self.assertTrue(os.path.isfile(os.path.join(self.pod, ".cache", "stalled.json")))
        rc, out, err, log = self.resume()
        self.assertEqual(rc, 0, err)
        self.assertIn(f"not resuming with {STALLS}", err)
        man = self.manifest()
        self.assertEqual(man["tts"]["models_used"], [ANSWERS])
        self.assertNotIn(STALLS, [m for m, c, _ in log["calls"] if not is_probe(c)])
        self.assertNotIn("made this episode's cached chunks", err)

    def test_without_the_memory_the_hand_over_skips_models_that_cannot_be_used(self):
        self.first_render_stops_on_a_stall()
        os.remove(os.path.join(self.pod, ".cache", "stalled.json"))       # an older cache
        rc, out, err, log = self.resume()
        self.assertEqual(rc, 0, err)
        self.assertNotIn(f"{self.PRO} made this episode's cached chunks", err)
        self.assertEqual(self.manifest()["tts"]["models_used"], [ANSWERS])
        # and the next resume does not wait on the stalled model again
        rc, out, err, log = self.resume()
        self.assertEqual(rc, 0, err)
        self.assertNotIn("still waiting for", err)


class KeyTier(Gemini):
    def free_tier_pro(self, model, contents, config):
        if "pro" in model:
            raise ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(FREE_DAY, value="0"))
        return pcm_for(contents)

    def test_a_free_tier_quota_is_recorded_and_a_paid_note_gets_one_warning(self):
        rc, out, err, log = self.render(self.free_tier_pro, "--cloud-ok",
                                        "PI, 2026-10-16, paid key")
        self.assertEqual(rc, 0, err)                                     # never a stop
        self.assertIn("the Gemini key's Google project is on the FREE tier", err)
        self.assertIn("generate_content_free_tier_requests", err)
        self.assertEqual(err.count("does not match the key"), 1)
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertEqual(man["key_tier"], "free")
        self.assertIn("free-tier quota", man["key_tier_basis"])
        self.assertTrue(any("does not match the key" in w for w in man["warnings"]))

    def test_a_note_that_says_free_or_nothing_gets_no_warning(self):
        for note in ("PI, 2026-10-16, free key", "PI, 2026-10-16", "PI, not paid"):
            with self.subTest(note=note):
                shutil.rmtree(os.path.join(self.pod, ".cache"), ignore_errors=True)   # pro first
                rc, out, err, log = self.render(self.free_tier_pro, "--cloud-ok", note)
                self.assertEqual(rc, 0, err)
                self.assertNotIn("does not match the key", err)
                man = json.loads(read(os.path.join(self.pod, "podcast.json")))
                self.assertEqual(man["key_tier"], "free")

    def test_no_tier_named_is_never_guessed(self):
        rc, out, err, log = self.render(lambda m, c, k: pcm_for(c), "--cloud-ok", "PI, paid key")
        self.assertEqual(rc, 0, err)
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertIsNone(man["key_tier"])
        self.assertIn("no answer from Google named", man["key_tier_basis"])
        self.assertNotIn("does not match the key", err)

    def test_error_info_reads_the_tier(self):
        e = ApiError(429, "RESOURCE_EXHAUSTED", rpc_error(FREE_DAY, value="0"))
        self.assertIn("generate_content_free_tier_requests", mp.error_info(e)["free_tier"])
        e = ApiError(429, "RESOURCE_EXHAUSTED", {"error": {"details": [{"violations": [
            {"quotaMetric": "x/generate_content_paid_requests", "quotaId": "PerMinute"}]}]}})
        self.assertIsNone(mp.error_info(e)["free_tier"])


class SdkOptions(Gemini):
    def sdk(self, with_options=True):
        mods, log = fake_genai(lambda m, c, k: pcm_for(c))
        gt, made = mods["google.genai.types"], {}
        if with_options:
            gt.HttpRetryOptions = lambda **kw: types.SimpleNamespace(**kw)
            gt.HttpOptions = lambda **kw: types.SimpleNamespace(**kw)
        base = mods["google.genai"].Client

        class Client(base):
            def __init__(self, api_key, http_options=None):
                made["http_options"] = http_options
                super().__init__(api_key)
        if with_options:
            mods["google.genai"].Client = Client
        return mods, made

    def test_a_request_timeout_and_no_sdk_retries_when_the_sdk_has_them(self):
        mods, made = self.sdk()
        with mock.patch.dict(sys.modules, mods):
            rc, out, err = run("render", self.script, "--tts", "gemini", "--no-verify",
                               "--no-preflight", "--cloud-ok", "--request-timeout", "240")
        self.assertEqual(rc, 0, err)
        self.assertEqual(made["http_options"].timeout, 240000)            # milliseconds
        self.assertEqual(made["http_options"].retry_options.attempts, 1)
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertEqual(man["tts"]["http"], {"timeout_ms": 240000, "retry_attempts": 1,
                                              "set": True})

    def test_an_sdk_without_them_still_renders_and_says_so(self):
        mods, _ = self.sdk(with_options=False)
        with mock.patch.dict(sys.modules, mods):
            rc, out, err = run("render", self.script, "--tts", "gemini", "--no-verify",
                               "--no-preflight", "--cloud-ok")
        self.assertEqual(rc, 0, err)
        self.assertIn("has no HttpOptions / HttpRetryOptions", err)
        man = json.loads(read(os.path.join(self.pod, "podcast.json")))
        self.assertFalse(man["tts"]["http"]["set"])


class OfflineRecord(unittest.TestCase):
    def test_an_offline_voice_says_why_there_is_no_test(self):
        """Review of 2.13.2 (L2): --tts say recorded "every chunk is cached" with none cached."""
        from test_make_podcast import FakeTTS, meaning_pass, script_text, write
        import tempfile
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "output")
            pod = os.path.join(out, "podcast")
            script = os.path.join(pod, "podcast_script.md")
            write(os.path.join(out, "AI_Analysis_Report.md"), __import__("test_make_podcast").REPORT)
            write(os.path.join(out, "Analysis_Report.html"), "<html></html>")
            write(script, script_text())
            self.assertEqual(run("check", script, "--source",
                                 os.path.join(out, "AI_Analysis_Report.md"))[0], 0)
            meaning_pass(script)
            mp.BACKENDS["fake"] = FakeTTS
            self.addCleanup(mp.BACKENDS.pop, "fake", None)
            with mock.patch.object(mp, "encode_aac", return_value=(None, "test")):
                self.assertEqual(run("render", script, "--tts", "fake")[0], 0)
            man = json.loads(read(os.path.join(pod, "podcast.json")))
        self.assertFalse(man["tts"]["preflight"]["run"])
        self.assertIn("no two-line test", man["tts"]["preflight"]["why"])


class Documented(unittest.TestCase):
    def test_the_brief_describes_the_fallback(self):
        with open(os.path.join(os.path.dirname(HERE), "references", "podcast.md"),
                  encoding="utf-8") as fh:
            t = " ".join(fh.read().split())
        for s in ("two-line test", "every chunk again", "--no-preflight", "key_tier"):
            self.assertIn(s, t)


if __name__ == "__main__":
    unittest.main(verbosity=2)
