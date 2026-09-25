#!/usr/bin/env python3
"""
FRAN hand-off, the half after staging: is the cron actually taking what the skill staged, is the
ingest code on HIVE the code on FRAN's main, and were any past searches never staged at all?

Every failure guarded here was silent in production. After 2026-09-17 the cron ran every 4 h and
ingested nothing in 41 consecutive runs (the recent ones "0 ingested, 3 duplicate-skipped, 2 failed,
~181 still queued") while none of the skill's drop entries appeared in one log; `verify` said
"staged_pending_cron", which was true and useless. A mouse search was staged with the HUMAN
database because three FASTA sidecars sat in one folder and the first one sorted wins. And the
skill's searches from before staging was automatic are simply not in FRAN.

No network and no HIVE: logs are synthetic in FRAN's exact format (checked against real logs and
auto_ingest.py on 2026-09-24), GitHub is a fake fetcher, and the service tree is a temp dir.
"""
import contextlib
import io
import json
import os
import sys
import tempfile
import time
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import fran_deposit as fd  # noqa: E402


def setUpModule():
    # verify() asks the LIVE corpus when this account can read a PG Farm token -- which a suite run
    # on HIVE can. fran_deposit.corpus_query() returns None before touching anything with this set.
    global _SAVED_CORPUS_QUERY
    _SAVED_CORPUS_QUERY = os.environ.get("FRAN_CORPUS_QUERY")
    os.environ["FRAN_CORPUS_QUERY"] = "off"


def tearDownModule():
    if _SAVED_CORPUS_QUERY is None:
        os.environ.pop("FRAN_CORPUS_QUERY", None)
    else:
        os.environ["FRAN_CORPUS_QUERY"] = _SAVED_CORPUS_QUERY


_SAVED_CORPUS_QUERY = None


def _load_json(path):
    """json.load with the file closed (no ResourceWarning)."""
    with open(path) as fh:
        return json.load(fh)


def _read(path, mode="r"):
    with open(path, mode) as fh:
        return fh.read()

# The QC name rule reads the last three components of a search dir's path, which in these tests
# include a random temp name. tempfile draws from [a-z0-9_], so "tmpx_qc_1" would make a test's
# search a "QC run" about once in 20,000 names. Without "_" in the alphabet the rule cannot fire on
# a temp name at all (it needs a non-alphanumeric character before "qc").
tempfile._RandomNameSequence.characters = "abcdefghijklmnopqrstuvwxyz0123456789"

NOW = time.mktime((2026, 9, 24, 12, 30, 0, 0, 0, -1))
H = 3600


def ts(e):
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(e))


def run_log(start, ingested=0, dup=0, failed=0, queued=0, items=(), skips=(), cands=(),
            done=True, abort=None):
    """One auto_ingest_<jobid>.out, in the exact shape FRAN's auto_ingest.py prints."""
    L = ["  preflight: PG Farm:5432 reachable",
         f"===== fran auto-ingest {ts(start)} on hive-as-11-2-70 limit=5 =====",
         f"===== auto_ingest {ts(start)} on hive-as-11-2-70 =====",
         "scanning: /envs/python /ingest/find_uningested.py --json-out /tmp/uningested_1.json",
         "4_142944_SpN_cut_mid_line/20260826_212930_cut",        # the scan tail is cut mid-line
         *[f"  diann        {c}" for c in cands],
         "  ... and 204 more", "", "wrote /tmp/uningested_1.json", "",
         "244 dirs -> 186 distinct searches (1 skipped)",
         *[f"  SKIP {n}  ({why})" for n, why in skips],
         'WARNING:  database "delimp" has a collation version mismatch', ""]
    if abort:
        L.append(abort)
    for i, (engine, d, outcome) in enumerate(items, 1):
        L += ["", f"[{i}/{len(items)}] {engine} {os.path.basename(d)[:52]}", f"      {d}"]
        if outcome == "ok":
            L.append("      OK in 120s")
        elif outcome == "dup":
            L += ["      SKIPPED-DUPLICATE in 245s (guard refused — not a failure)",
                  "      exists: R:\\Data\\lab\\service\\x.sne"]
        elif outcome == "fail":
            L += ["      FAILED rc=1 in 27s", "      --- last output ---", "      | --- stderr ---",
                  '      | WARNING:  database "delimp" has a collation version mismatch',
                  "      | HINT:  Rebuild all objects in this database",
                  "      | No precursor records parsed (check the report / --engine)."]
    if done:
        L.append(f"===== done: {ingested} ingested, {dup} duplicate-skipped, {failed} failed, "
                 f"{queued} still queued — {ts(start + 600)} =====")
        L.append(f"===== auto-ingest exit rc=0 {ts(start + 601)} =====")
    return "\n".join(L) + "\n"


class Env:
    """A fake HIVE: log dir, submit log, drop dir, all under one temp root, wired in by env."""
    def __init__(self, root, now=NOW):
        # Every mtime and the submit line are relative to `now`: NOW for tests that hand
        # fran_deposit the same `now`; time.time() for tests of the CLI entry points (health,
        # verify), which read the real clock -- a fixed date there goes stale the next day.
        self.root = root
        self.now = now
        self.logs = os.path.join(root, "logs")
        self.drop = os.path.join(root, "incoming")
        os.makedirs(self.logs)
        os.makedirs(self.drop)
        self.job = 1000

    def log(self, text, age_h=None):
        self.job += 1
        p = os.path.join(self.logs, f"auto_ingest_{self.job}.out")
        with open(p, "w") as fh:
            fh.write(text)
        if age_h is not None:
            os.utime(p, (self.now - age_h * H, self.now - age_h * H))
        return p

    def submit(self, age_h=0.2, line="Submitted batch job 23990862"):
        p = os.path.join(self.logs, fd.SUBMIT_LOG)
        with open(p, "w") as fh:
            fh.write(f"{ts(self.now - age_h * H - 420)} {line}\n")
        os.utime(p, (self.now - age_h * H, self.now - age_h * H))

    def entry(self, name, age_h=100, output_dir=None):
        d = os.path.join(self.drop, name)
        os.makedirs(d)
        with open(os.path.join(d, fd.MANIFEST), "w") as fh:
            json.dump({"output_dir": output_dir or f"/real/{name}", "engine": "diann",
                       "staged_by": "brettsp"}, fh)
        os.utime(d, (self.now - age_h * H, self.now - age_h * H))
        return d


@contextlib.contextmanager
def env_vars(**kv):
    old = {k: os.environ.get(k) for k in kv}
    for k, v in kv.items():
        if v is None:
            os.environ.pop(k, None)
        else:
            os.environ[k] = v
    try:
        yield
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def runs_of(env):
    runs, info = fd.load_ingest_logs(env.logs)
    return runs, info


# ------------------------------------------------------------------------------ progress --
class ProgressTests(unittest.TestCase):
    def test_the_real_stall_is_called_stuck(self):
        """The 2026-09-17..24 pattern: every run 0 ingested with ~181 queued, the last real ingest
        a week earlier. Twelve healthy-looking COMPLETED jobs; one verdict: stuck."""
        with tempfile.TemporaryDirectory() as d:
            e = Env(d)
            e.log(run_log(NOW - 170 * H, ingested=4, dup=1, queued=198))
            for k in range(12, 0, -1):
                e.log(run_log(NOW - k * 4 * H, dup=3, failed=2, queued=181))
            e.submit()
            runs, info = runs_of(e)
            p = fd.progress_health(runs, NOW, fd.submit_log_info(e.logs, NOW), info)
            self.assertEqual(p["verdict"], "stuck", p)
            self.assertEqual(p["consecutive_runs_without_ingest"], 12)
            self.assertEqual(p["queued"], 181)
            self.assertIn("181 searches still queued", p["detail"])
            self.assertEqual(p["last_ingest"]["ingested"], 4)
            self.assertTrue(p["cron"]["submitting"])

    def test_a_run_that_ingests_is_healthy(self):
        with tempfile.TemporaryDirectory() as d:
            e = Env(d)
            for k in (3, 2):
                e.log(run_log(NOW - k * 4 * H, dup=3, queued=50))
            e.log(run_log(NOW - 3 * H, ingested=2, queued=48))
            e.submit()
            runs, info = runs_of(e)
            p = fd.progress_health(runs, NOW, fd.submit_log_info(e.logs, NOW), info)
            self.assertEqual(p["verdict"], "healthy", p)
            self.assertEqual(p["consecutive_runs_without_ingest"], 0)

    def test_an_idle_cron_with_nothing_queued_is_not_stuck(self):
        """Zero ingested because there is nothing to do is the GOOD state."""
        with tempfile.TemporaryDirectory() as d:
            e = Env(d)
            for k in (4, 3, 2, 1):
                e.log(run_log(NOW - k * 4 * H, queued=0))
            e.submit()
            runs, info = runs_of(e)
            self.assertEqual(fd.progress_health(runs, NOW, fd.submit_log_info(e.logs, NOW),
                                                info)["verdict"], "healthy")

    def test_no_run_for_half_a_day_is_not_running(self):
        with tempfile.TemporaryDirectory() as d:
            e = Env(d)
            e.log(run_log(NOW - 20 * H, ingested=3, queued=0))
            e.submit(age_h=20)
            runs, info = runs_of(e)
            p = fd.progress_health(runs, NOW, fd.submit_log_info(e.logs, NOW), info)
            self.assertEqual(p["verdict"], "not_running", p)
            self.assertFalse(p["cron"]["submitting"])

    def test_a_cron_that_stopped_submitting_is_not_running_even_with_a_fresh_log(self):
        with tempfile.TemporaryDirectory() as d:
            e = Env(d)
            e.log(run_log(NOW - 2 * H, ingested=3, queued=0))
            e.submit(age_h=30)
            runs, info = runs_of(e)
            self.assertEqual(fd.progress_health(runs, NOW, fd.submit_log_info(e.logs, NOW),
                                                info)["verdict"], "not_running")

    def test_a_run_in_progress_is_not_counted_and_dead_runs_are(self):
        """The newest log has no summary yet -- it is running, not failing. Runs that died (PG
        Farm unreachable) ingest nothing and DO count, or a cron that aborts every time would look
        fine forever."""
        with tempfile.TemporaryDirectory() as d:
            e = Env(d)
            e.log(run_log(NOW - 40 * H, ingested=1, queued=20))
            for k in (3, 2, 1):
                e.log("  preflight FAILED: timed out\nABORT: PG Farm unreachable from hive-x\n",
                      age_h=k * 4)
            e.log(run_log(NOW - 0.1 * H, done=False))
            e.submit()
            runs, info = runs_of(e)
            p = fd.progress_health(runs, NOW, fd.submit_log_info(e.logs, NOW), info)
            self.assertEqual(p["last_run"]["state"], "running")
            self.assertEqual(p["consecutive_runs_without_ingest"], 3)
            self.assertEqual(p["of_which_aborted"], 3)
            self.assertEqual(p["verdict"], "stuck")

    def test_unreadable_logs_are_unknown_not_an_exception(self):
        runs, info = fd.load_ingest_logs("/definitely/not/a/log/dir")
        self.assertIsNone(runs)
        self.assertEqual(fd.progress_health(runs, NOW, None, info)["verdict"], "unknown")


# ------------------------------------------------------------------------ one entry's fate --
class EntryStateTests(unittest.TestCase):
    ENTRY = "/quobyte/proteomics-grp/fran/incoming/search_out__e14aac29"

    def _runs(self, *texts):
        return [dict(fd.parse_ingest_log(t), log=f"L{i}", mtime=NOW, jobid=i)
                for i, t in enumerate(texts)]

    def test_ingested(self):
        runs = self._runs(run_log(NOW - H, ingested=1, items=[("diann", self.ENTRY, "ok")]))
        s = fd.entry_log_state(runs, self.ENTRY)
        self.assertEqual((s["state"], s["outcome"]), ("ingested", "ok"))

    def test_a_duplicate_counts_as_in_fran(self):
        """The guard refused it because the corpus already holds the same search."""
        runs = self._runs(run_log(NOW - H, dup=1, items=[("diann", self.ENTRY, "dup")]))
        s = fd.entry_log_state(runs, self.ENTRY)
        self.assertEqual((s["state"], s["outcome"]), ("ingested", "duplicate"))

    def test_failed_reports_the_real_reason_not_the_collation_noise(self):
        runs = self._runs(run_log(NOW - H, failed=1, items=[("diann", self.ENTRY, "fail")]))
        s = fd.entry_log_state(runs, self.ENTRY)
        self.assertEqual(s["state"], "failed")
        self.assertIn("No precursor records parsed", s["detail"])
        self.assertNotIn("collation", s["detail"])

    def test_a_later_success_beats_an_earlier_failure(self):
        runs = self._runs(run_log(NOW - H, ingested=1, items=[("diann", self.ENTRY, "ok")]),
                          run_log(NOW - 5 * H, failed=1, items=[("diann", self.ENTRY, "fail")]))
        self.assertEqual(fd.entry_log_state(runs, self.ENTRY)["state"], "ingested")

    def test_skipped_by_name_is_a_failure_with_the_reason(self):
        runs = self._runs(run_log(NOW - H, skips=[("search_out__e14aac29",
                                                   "no usable report in 1 export(s) — empty/failed")]))
        s = fd.entry_log_state(runs, self.ENTRY)
        self.assertEqual((s["state"], s["outcome"]), ("failed", "skipped"))
        self.assertIn("no usable report", s["detail"])

    def test_never_reached_even_when_listed_as_a_candidate(self):
        """The head-of-line case: the scan finds the entry every run, the limit=5 slice never
        reaches it. Being LISTED is not being attempted."""
        runs = self._runs(*[run_log(NOW - k * H, dup=3, failed=2, queued=181, cands=[self.ENTRY],
                                    items=[("spectronaut", f"/FRAN_reports/x{j}", "dup")
                                           for j in range(3)])
                            for k in range(1, 4)])
        s = fd.entry_log_state(runs, self.ENTRY)
        self.assertEqual(s["state"], "never_reached")
        self.assertEqual(s["listed_as_candidate"], 3)
        self.assertEqual(s["attempts"], 0)

    def test_a_plain_folder_name_never_matches_another_search(self):
        """Every parallel-chain search is called `search_out`; only entry names are unique."""
        runs = self._runs(run_log(NOW - H, ingested=1,
                                  items=[("diann", "/other/lab/search_out", "ok")],
                                  skips=[("search_out", "named fail*")]))
        s = fd.entry_log_state(runs, "/quobyte/SERVICE/mine/search_out")
        self.assertEqual(s["state"], "never_reached")

    def test_matched_by_the_real_search_dir_too(self):
        """A queue row or a manifest-aware cron prints the search's REAL dir as the identity."""
        real = "/quobyte/proteomics-grp/SERVICE/on_campus/X/search_out"
        text = run_log(NOW - H, ingested=1, items=[("diann", "/elsewhere/q", "ok")]).replace(
            "      /elsewhere/q\n", f"      /elsewhere/q\n      -> {real}\n")
        runs = self._runs(text)
        self.assertEqual(fd.entry_log_state(runs, self.ENTRY, real)["state"], "ingested")


class IncomingTests(unittest.TestCase):
    def test_entries_never_reached_for_days_are_starved(self):
        """The cron can be ingesting OTHER searches while the skill's entries wait forever at the
        back of an alphabetical queue. That has to surface even when progress looks fine."""
        with tempfile.TemporaryDirectory() as d:
            e = Env(d)
            old = e.entry("GallPlasCer__5b11a0d9", age_h=29 * 24)
            e.entry("search__9ff203cf", age_h=20)
            e.log(run_log(NOW - H, ingested=2, queued=100,
                          items=[("spectronaut", "/FRAN_reports/a", "ok")]))
            runs, _ = runs_of(e)
            inc = fd.incoming_health(runs, NOW, e.drop)
            self.assertEqual(inc["verdict"], "starved", inc)
            self.assertEqual(inc["n_entries"], 2)
            self.assertEqual(inc["by_state"], {"never_reached": 2})
            self.assertAlmostEqual(inc["oldest_never_reached_days"], 29, delta=0.2)
            self.assertTrue(os.path.isdir(old))

    def test_broken_links_are_listed(self):
        with tempfile.TemporaryDirectory() as d:
            e = Env(d)
            ent = e.entry("x__00000000", age_h=1)
            os.symlink(os.path.join(d, "gone.parquet"), os.path.join(ent, "report.parquet"))
            inc = fd.incoming_health([], NOW, e.drop)
            self.assertEqual(inc["entries"][0]["broken_links"], ["report.parquet"])


class HealthReportTests(unittest.TestCase):
    def test_overall_verdict_and_one_line_summary(self):
        with tempfile.TemporaryDirectory() as d:
            e = Env(d)
            e.entry("search_out__e14aac29", age_h=72)
            e.log(run_log(NOW - 170 * H, ingested=4, queued=198))
            for k in (3, 2, 1):
                e.log(run_log(NOW - k * 4 * H, dup=3, failed=2, queued=181))
            e.submit()
            with env_vars(FRAN_INGEST_LOG_DIR=e.logs, FRAN_DROP_DIR=e.drop):
                h = fd.health_report(now=NOW, check_code=False)
            self.assertEqual(h["verdict"], "stuck")
            self.assertFalse(h["healthy"])
            self.assertNotIn("\n", h["summary"])
            self.assertIn("STUCK", h["summary"])
            self.assertIn("never reached", h["summary"])

    def test_nothing_readable_is_unknown(self):
        with env_vars(FRAN_INGEST_LOG_DIR="/nope/logs", FRAN_DROP_DIR="/nope/drop"):
            h = fd.health_report(now=NOW, check_code=False)
        self.assertEqual(h["verdict"], "unknown")


# ---------------------------------------------------------------------------- ingest code --
class FakeGitHub:
    """raw.githubusercontent.com + the commits API, from dicts. Counts every call."""
    def __init__(self, files_at, commits, down=False, rate_limited=False):
        self.files_at, self.commits = files_at, commits   # {ref: {name: bytes}}, {name: [sha..]}
        self.down, self.rate_limited, self.calls = down, rate_limited, []

    def __call__(self, url, timeout):
        self.calls.append(url)
        if self.down:
            return None, "URLError: no route to host", {}
        if "api.github.com" in url:
            if self.rate_limited:
                return 403, b'{"message":"API rate limit exceeded"}', {
                    "x-ratelimit-remaining": "0", "x-ratelimit-reset": str(int(NOW + 1800))}
            name = url.split("path=ingest/")[1].split("&")[0]
            return 200, json.dumps([{"sha": s, "commit": {"committer": {"date": f"2026-09-{10 + i:02d}T10:00:00Z"}}}
                                    for i, s in enumerate(reversed(self.commits.get(name, [])))][::-1]
                                   ).encode(), {}
        ref, name = url.split("/FRAN/")[1].split("/ingest/")
        body = self.files_at.get(ref, {}).get(name)
        return (200, body, {}) if body is not None else (404, b"404: Not Found", {})


class IngestCodeTests(unittest.TestCase):
    V1, V2, LOCAL = b"version = 1\n", b"version = 2\n", b"version = 2  # edited on HIVE\n"

    def _run(self, hive, gh, files=("corpus_ingest.py", "raw_metadata.py"), cache=None):
        with tempfile.TemporaryDirectory() as d:
            ing = os.path.join(d, "fran_ingest")
            os.makedirs(ing)
            for n, body in hive.items():
                with open(os.path.join(ing, n), "wb") as fh:
                    fh.write(body)
            return fd.ingest_code_health(ing, fetch=gh, now=NOW, files=files,
                                         cache_path=cache or os.path.join(d, "cache.json"))

    def _gh(self, **kw):
        return FakeGitHub({"main": {"corpus_ingest.py": self.V2, "raw_metadata.py": self.V2},
                           "c2": {"corpus_ingest.py": self.V2, "raw_metadata.py": self.V2},
                           "c1": {"corpus_ingest.py": self.V1, "raw_metadata.py": self.V1}},
                          {"corpus_ingest.py": ["c2", "c1"], "raw_metadata.py": ["c2", "c1"]}, **kw)

    def test_current(self):
        r = self._run({"corpus_ingest.py": self.V2, "raw_metadata.py": self.V2}, self._gh())
        self.assertEqual(r["verdict"], "current")
        self.assertEqual(r["by_status"], {"current": ["corpus_ingest.py", "raw_metadata.py"]})

    def test_stale_names_the_older_commit_it_matches(self):
        r = self._run({"corpus_ingest.py": self.V1, "raw_metadata.py": self.V2}, self._gh())
        f = next(x for x in r["files"] if x["file"] == "corpus_ingest.py")
        self.assertEqual(f["status"], "stale")
        self.assertEqual(f["matches"]["sha"], "c1")
        self.assertIn("matches c1 from 2026-09-", f["detail"])
        self.assertEqual(r["verdict"], "stale")

    def test_an_edit_nobody_pushed_is_a_local_modification(self):
        """What raw_metadata.py looked like on HIVE on 2026-09-24: in-progress work, matching no
        commit. Outside the refuse gate that is a warning (`modified`), not a stale verdict."""
        r = self._run({"corpus_ingest.py": self.V2, "raw_metadata.py": self.LOCAL}, self._gh())
        f = next(x for x in r["files"] if x["file"] == "raw_metadata.py")
        self.assertEqual(f["status"], "local_modification")
        self.assertEqual(r["verdict"], "modified")

    def test_a_local_edit_to_a_refuse_gated_file_is_stale(self):
        """corpus_ingest.py is on FRAN's own refuse list: rows written by code nobody committed
        are rows nobody can reproduce."""
        r = self._run({"corpus_ingest.py": self.LOCAL, "raw_metadata.py": self.V2}, self._gh())
        self.assertEqual(r["verdict"], "stale")

    def test_no_network_is_unknown_never_an_exception(self):
        r = self._run({"corpus_ingest.py": self.V2, "raw_metadata.py": self.V2}, self._gh(down=True))
        self.assertEqual(r["verdict"], "unknown")
        self.assertTrue(all(f["status"] == "unknown" for f in r["files"]))

    def test_rate_limit_is_reported_and_remembered(self):
        """60 API requests/h is shared by everyone behind a login node's address. Hitting the
        limit leaves the file 'differs', and the cache stops the next run asking again."""
        with tempfile.TemporaryDirectory() as d:
            cache = os.path.join(d, "c.json")
            gh = self._gh(rate_limited=True)
            r = self._run({"corpus_ingest.py": self.V2, "raw_metadata.py": self.LOCAL}, gh, cache=cache)
            f = next(x for x in r["files"] if x["file"] == "raw_metadata.py")
            self.assertEqual(f["status"], "differs")
            self.assertIn("rate limit", f["detail"])
            api = sum("api.github.com" in u for u in gh.calls)
            gh2 = self._gh()
            self._run({"corpus_ingest.py": self.V2, "raw_metadata.py": self.LOCAL}, gh2, cache=cache)
            self.assertEqual(api, 1)
            self.assertEqual(sum("api.github.com" in u for u in gh2.calls), 0)
            self.assertEqual(sum("/main/" in u for u in gh2.calls), 0, "main should be cached")

    def test_missing_on_hive_and_not_on_main(self):
        gh = self._gh()
        gh.files_at["main"]["engine_fasta.py"] = self.V1
        r = self._run({"corpus_ingest.py": self.V2}, gh,
                      files=("corpus_ingest.py", "engine_fasta.py", "diann_to_corpus.py"))
        st = {f["file"]: f["status"] for f in r["files"]}
        self.assertEqual(st, {"corpus_ingest.py": "current", "engine_fasta.py": "missing",
                              "diann_to_corpus.py": "not_on_main"})
        self.assertEqual(r["verdict"], "stale")

    def test_a_github_that_never_answers_costs_at_most_the_budget(self):
        """urlopen's timeout does not cover a hung DNS lookup; the check must return anyway."""
        import threading
        never = threading.Event()

        def hang(url, timeout):
            never.wait()                          # a lookup that never comes back
            return None, "released at test end", {}
        with tempfile.TemporaryDirectory() as d:
            ing = os.path.join(d, "ing")
            os.makedirs(ing)
            with open(os.path.join(ing, "corpus_ingest.py"), "wb") as fh:
                fh.write(self.V2)
            t0 = time.monotonic()
            r = fd.ingest_code_health(ing, fetch=hang, now=NOW, budget_s=1,
                                      files=("corpus_ingest.py", "versions.py"),
                                      cache_path=os.path.join(d, "c.json"))
            self.assertLess(time.monotonic() - t0, 3.0)
        self.assertEqual(r["verdict"], "unknown")
        self.assertIn("no answer from GitHub within 1 s", r["files"][0]["detail"])
        never.set()

    def test_the_default_budget_is_ten_seconds_and_no_credential_is_sent(self):
        import inspect
        self.assertEqual(inspect.signature(fd.ingest_code_health).parameters["budget_s"].default, 10)
        self.assertNotIn("Authorization", inspect.getsource(fd._default_fetch))

    def test_no_ingest_dir_asks_nothing(self):
        gh = self._gh()
        with env_vars(FRAN_INGEST_DIR=None):
            saved = fd.INGEST_DIRS
            fd.INGEST_DIRS = ["/no/such/ingest"]
            try:
                r = fd.ingest_code_health(fetch=gh, now=NOW)
            finally:
                fd.INGEST_DIRS = saved
        self.assertEqual(r["verdict"], "unknown")
        self.assertEqual(gh.calls, [])


# --------------------------------------------------------------------- stage + verify use it --
def search_dir(root, name="search_out", engine="diann", empty=False, prov=True, extra=()):
    d = os.path.join(root, name)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "report.parquet"), "w") as fh:
        fh.write("" if empty else "x" * 64)
    if prov:
        with open(os.path.join(d, "search_provenance.json"), "w") as fh:
            json.dump({"engine": engine, "version": "2.7.0"}, fh)
    for f in extra:
        with open(os.path.join(d, f), "w") as fh:
            fh.write("x")
    return d


class Args:
    def __init__(self, out, **kw):
        self.out = out
        for k, v in dict(skip=False, force=False, dry_run=False, organism=None, taxon=None,
                         fasta_meta=None, name=None, qc=False, not_qc=False).items():
            setattr(self, k, kw.get(k, v))


def run_quiet(fn, *a):
    out, err = io.StringIO(), io.StringIO()
    with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
        try:
            fn(*a)
        except SystemExit as e:
            code = e.code
        else:
            code = None
    return code, json.loads(out.getvalue()) if out.getvalue().strip() else None, err.getvalue()


class StageHealthTests(unittest.TestCase):
    """stage runs inside every search job. It READS the status file `health` writes; it never
    reaches GitHub, the database, or the cron's logs."""

    def _status(self, drop, verdict, age_h=0.5, summary="FRAN ingest STUCK: 41 runs"):
        os.makedirs(drop, exist_ok=True)
        t = time.time() - age_h * H
        with open(fd.health_file(drop), "w") as fh:
            json.dump({"checked_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t)),
                       "verdict": verdict, "summary": summary}, fh)

    def _stage(self, d, **kw):
        out = search_dir(d)
        with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"), FRAN_HEALTH=None):
            return run_quiet(fd.stage, Args(out, **kw))

    def test_an_unhealthy_verdict_warns_but_never_fails_the_stage(self):
        with tempfile.TemporaryDirectory() as d:
            self._status(os.path.join(d, "incoming"), "stuck")
            code, res, err = self._stage(d)
            self.assertEqual(code, 0)
            self.assertTrue(res["staged"])
            self.assertEqual(res["health_warning"], "FRAN ingest STUCK: 41 runs")
            self.assertEqual(res["fran_health"]["verdict"], "stuck")
            self.assertIn("WARNING: FRAN ingest STUCK: 41 runs", err)

    def test_a_healthy_verdict_is_carried_quietly(self):
        with tempfile.TemporaryDirectory() as d:
            self._status(os.path.join(d, "incoming"), "healthy", summary="ok")
            code, res, err = self._stage(d)
            self.assertEqual(res["fran_health"]["verdict"], "healthy")
            self.assertNotIn("health_warning", res)
            self.assertNotIn("WARNING", err)

    def test_a_stale_verdict_is_unknown_and_silent(self):
        """Old news is not bad news: a >12 h verdict must not make the job's Slack post call FRAN
        unhealthy. Same as a missing file, plus the date it was last checked."""
        for verdict in ("healthy", "stuck"):
            with tempfile.TemporaryDirectory() as d:
                self._status(os.path.join(d, "incoming"), verdict, age_h=13)
                code, res, err = self._stage(d)
                self.assertTrue(res["staged"])
                self.assertEqual(res["fran_health"]["verdict"], "unknown")
                self.assertIn("last checked", res["fran_health"]["summary"])
                self.assertNotIn("health_warning", res)
                self.assertNotIn("WARNING", err)

    def test_no_status_file_or_a_broken_one_says_nothing(self):
        with tempfile.TemporaryDirectory() as d:
            code, res, err = self._stage(d)
            self.assertTrue(res["staged"])
            self.assertNotIn("fran_health", res)
            self.assertNotIn("health_warning", res)
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(os.path.join(d, "incoming"))
            with open(fd.health_file(os.path.join(d, "incoming")), "w") as fh:
                fh.write("{not json")
            code, res, err = self._stage(d)
            self.assertTrue(res["staged"])
            self.assertNotIn("fran_health", res)

    def test_stage_makes_no_network_call_and_reads_no_logs(self):
        """The job-end hook runs stage on a compute node inside the job's time limit. Any network,
        database or log-scan call would now raise -- and stage must still stage."""
        import socket
        import urllib.request

        def banned(*a, **kw):
            raise AssertionError("stage reached for the network / a full health check")
        patches = [(urllib.request, "urlopen"), (socket, "socket"), (socket, "create_connection"),
                   (socket, "getaddrinfo"), (fd, "load_ingest_logs"), (fd, "health_report"),
                   (fd, "ingest_code_health"), (fd, "corpus_query"), (fd, "_default_fetch")]
        saved = [(m, n, getattr(m, n)) for m, n in patches]
        with tempfile.TemporaryDirectory() as d:
            self._status(os.path.join(d, "incoming"), "stuck")
            try:
                for m, n in patches:
                    setattr(m, n, banned)
                code, res, err = self._stage(d)
            finally:
                for m, n, f in saved:
                    setattr(m, n, f)
            self.assertEqual(code, 0, err)
            self.assertTrue(res["staged"])
            self.assertEqual(res["health_warning"], "FRAN ingest STUCK: 41 runs")

    def test_a_hung_mount_cannot_hold_the_job(self):
        """A status file that never answers (a FIFO with no writer blocks open() like a hung NFS
        mount) costs at most the read timeout, then stage says nothing."""
        with tempfile.TemporaryDirectory() as d:
            drop = os.path.join(d, "incoming")
            os.makedirs(drop)
            fifo = fd.health_file(drop)
            os.mkfifo(fifo)
            t0 = time.monotonic()
            got = fd.read_health_status(drop=drop, timeout_s=0.5)
            self.assertLess(time.monotonic() - t0, 2.0)
            self.assertEqual(got, (None, None))
            try:                                  # release the reader thread
                os.close(os.open(fifo, os.O_WRONLY | os.O_NONBLOCK))
            except OSError:
                pass

    def test_health_writes_the_file_stage_reads(self):
        with tempfile.TemporaryDirectory() as d:
            e = Env(d, now=time.time())          # fd.health reads the real clock
            e.log(run_log(e.now - 170 * H, ingested=4, queued=198))
            for k in (3, 2, 1):
                e.log(run_log(e.now - k * 4 * H, dup=3, failed=2, queued=181))
            e.submit(age_h=0)

            class A:
                no_code, alert = True, False
            with env_vars(FRAN_INGEST_LOG_DIR=e.logs, FRAN_DROP_DIR=e.drop):
                code, h, _ = run_quiet(fd.health, A())
                fh, warn = fd.read_health_status()
            self.assertEqual(h["status_file"]["written"], os.path.join(d, "ingest_health.json"))
            self.assertEqual(oct(os.stat(h["status_file"]["written"]).st_mode & 0o777), "0o664")
            # beside the drop dir, never in it: incoming/ holds drop entries and nothing else
            self.assertEqual([f for f in os.listdir(e.drop) if not os.path.isdir(
                os.path.join(e.drop, f))], [])
            self.assertEqual(fh["verdict"], "stuck")
            self.assertEqual(warn, h["summary"])

    def test_an_opt_out_is_recorded_so_backfill_honours_it(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming")):
                code, res, _ = run_quiet(fd.stage, Args(out, skip=True))
                self.assertEqual(res["reason"], "opted_out")
                self.assertEqual(fd.read_receipt(out)["status"], "opted_out")
                # ...which does not block an explicit stage later
                self.assertTrue(fd.check(Args(out))["eligible"])

    def test_an_opt_out_never_overwrites_a_staged_receipt(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"), FRAN_HEALTH="off"):
                run_quiet(fd.stage, Args(out))
                run_quiet(fd.stage, Args(out, skip=True))
            self.assertEqual(fd.read_receipt(out)["status"], "staged")


class StagedAtTests(unittest.TestCase):
    """FRAN's runner takes drop entries oldest-first by `staged_at`; its fallback is the entry
    directory's mtime, which every re-stage resets."""
    ISO = r"^\d{4}-\d\d-\d\dT\d\d:\d\d:\d\dZ$"

    def _stage(self, d, out, **kw):
        with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"), FRAN_HEALTH="off"):
            code, res, _ = run_quiet(fd.stage, Args(out, **kw))
        return res, _load_json(os.path.join(res["entry"], fd.MANIFEST))

    def test_staged_at_is_iso_utc_and_every_old_key_is_still_there(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            res, man = self._stage(d, out, organism="Homo sapiens", taxon=9606, name="Run1")
            self.assertRegex(man["staged_at"], self.ISO)
            self.assertEqual(man["staged_by"], __import__("getpass").getuser())
            for k in ("output_dir", "search_name", "organism", "taxon", "engine", "report",
                      "suggested_ingest", "suggested_ingest_shell", "xic", "linked"):
                self.assertIn(k, man)
            self.assertEqual(fd.read_receipt(out)["staged_at"], man["staged_at"])

    def test_a_restage_keeps_the_original_time(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            _, first = self._stage(d, out)
            mp = os.path.join(os.path.join(d, "incoming"), fd.entry_name(out), fd.MANIFEST)
            first["staged_at"] = "2026-08-26T17:40:50Z"          # as if staged a month ago
            with open(mp, "w") as fh:
                json.dump(first, fh)
            _, again = self._stage(d, out, force=True)
            self.assertEqual(again["staged_at"], "2026-08-26T17:40:50Z")
            self.assertRegex(again["restaged_at"], self.ISO)

    def test_an_entry_from_before_staged_at_keeps_its_manifest_time(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            _, first = self._stage(d, out)
            mp = os.path.join(os.path.join(d, "incoming"), fd.entry_name(out), fd.MANIFEST)
            del first["staged_at"]
            with open(mp, "w") as fh:
                json.dump(first, fh)
            when = time.mktime((2026, 8, 26, 10, 40, 50, 0, 0, -1))
            os.utime(mp, (when, when))
            _, again = self._stage(d, out, force=True)
            self.assertEqual(again["staged_at"],
                             time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(when)))


# --------------------------------------------------------------- FRAN's manifest contract --
# A COPY of the validation in FRAN's ingest/find_uningested.py `read_manifest` / `_staged_epoch`
# (branch fix/auto-ingest-starvation, 2026-09-24). It cannot be imported across repos, so it is
# copied, and must be updated when FRAN's changes. What it enforces: fran_manifest_version 1 or
# absent; an absolute output_dir; a known engine; search_name / organism a non-empty string or
# absent; taxon numeric; qc / exclude JSON booleans ("anything else makes the manifest
# malformed"); staged_at ISO 8601 or epoch.
_FRAN_ENGINES = ("diann", "spectronaut", "fragpipe", "radiant")


def _fran_staged_epoch(v):
    import datetime as _dt
    if v is None or isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v) if v > 0 else None
    try:
        return float(str(v).strip())
    except ValueError:
        pass
    try:
        return _dt.datetime.fromisoformat(str(v).strip().replace("Z", "+00:00")).timestamp()
    except (ValueError, OverflowError, OSError):
        return None


def fran_read_manifest(d, detected_engine=None):
    try:
        with open(os.path.join(d, "fran_manifest.json"), encoding="utf-8") as fh:
            m = json.load(fh)
    except FileNotFoundError:
        return None, "no fran_manifest.json"
    except (OSError, ValueError) as e:
        return None, f"unreadable ({type(e).__name__})"
    if not isinstance(m, dict):
        return None, "not a JSON object"
    if m.get("fran_manifest_version", 1) not in (1,):
        return None, "version"
    od = m.get("output_dir")
    if not isinstance(od, str) or not od.strip() or not os.path.isabs(od.strip()):
        return None, "no absolute output_dir"
    eng = m.get("engine")
    if eng is not None:
        if not isinstance(eng, str) or eng.strip().lower() not in _FRAN_ENGINES:
            return None, f"unknown engine {eng!r}"
        eng = eng.strip().lower()
        if detected_engine and eng != detected_engine:
            return None, "engine mismatch"
    out = {"output_dir": od.strip().rstrip("/"), "engine": eng or detected_engine}
    for k in ("search_name", "organism"):
        v = m.get(k)
        if v is not None and (not isinstance(v, str) or not v.strip()):
            return None, f"field {k!r} is not a non-empty string"
        out[k] = v.strip() if v is not None else None
    tx = m.get("taxon")
    if tx is None or tx == "":
        out["taxon"] = None
    elif isinstance(tx, bool) or not str(tx).strip().isdigit():
        return None, f"taxon {tx!r}"
    else:
        out["taxon"] = str(tx).strip()
    for k in ("qc", "exclude"):
        v = m.get(k)
        if v is not None and not isinstance(v, bool):
            return None, f"{k} {v!r} is not true/false"
        out[k] = v
    for k in ("fasta_path", "staged_by"):
        v = m.get(k)
        out[k] = v.strip() if isinstance(v, str) and v.strip() else None
    out["staged_at"] = _fran_staged_epoch(m.get("staged_at"))
    return out, None


class FranManifestContractTests(unittest.TestCase):
    """Every manifest stage can write must pass FRAN's reader, or FRAN skips the entry as
    malformed -- silently, from the skill's side."""

    def _stage(self, d, out, **kw):
        with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"), FRAN_HEALTH="off"):
            return run_quiet(fd.stage, Args(out, **kw))[1]

    def _read(self, entry):
        m, why = fran_read_manifest(entry, fd.detect_engine(entry)[0])
        self.assertIsNone(why, why)
        return m

    def test_a_full_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            res = self._stage(d, out, organism="Homo sapiens", taxon=9606, name="Run 1")
            m = self._read(res["entry"])
            self.assertEqual((m["organism"], m["taxon"], m["search_name"]),
                             ("Homo sapiens", "9606", "Run 1"))
            self.assertIs(m["qc"], False)
            self.assertAlmostEqual(m["staged_at"], time.time(), delta=120)
            self.assertEqual(m["output_dir"], os.path.realpath(out))

    def test_blank_name_and_organism_are_absent_not_empty(self):
        """FRAN rejects "" for search_name / organism. An empty --name must not produce one."""
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            m = self._read(self._stage(d, out, name="  ", organism="")["entry"])
            self.assertIsNone(m["search_name"])
            self.assertIsNone(m["organism"])

    def test_not_qc_and_withdrawn_manifests(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(os.path.join(d, "lab", "proj"), "hela_qc_2")
            m = self._read(self._stage(d, out, not_qc=True)["entry"])
            self.assertIs(m["qc"], False)
            res = self._stage(d, out, qc=True)                   # withdraws it
            m = self._read(res["withdrawn"])
            self.assertEqual((m["qc"], m["exclude"]), (True, True))


class QcRuleTests(unittest.TestCase):
    """QC runs are never handed to FRAN (Brett, 2026-09-24). ONE rule, is_qc_run(), mirrored by
    FRAN's ingest/find_uningested.py `qc_reason` / `QC_NAME_RE`; FRAN pins the same vectors."""
    EXCLUDED = ["chkLUppm_HeLa50_2026 Lumos QC", "QC_run_01", "hela_qc_2", "Exploris QC2"]
    KEPT = ["HeLa_digest_timecourse", "aqc_buffer_study", "QCM_study", "Plasma_liver2"]

    def test_the_pinned_vectors(self):
        for n in self.EXCLUDED:
            self.assertTrue(fd.QC_NAME_RE.search(n), n)
        for n in self.KEPT:
            self.assertFalse(fd.QC_NAME_RE.search(n), n)
        self.assertEqual(fd.QC_NAME_RE.pattern, r"(?i)(?<![a-z0-9])qc(?![a-z])")

    def _session(self, d, title, qc_marker=None, conditions=False):
        """session.py's layout: <session>/README.md, input/, output/search."""
        sess = os.path.join(d, "lab", title.replace(" ", "_"))
        os.makedirs(os.path.join(sess, "input"))
        with open(os.path.join(sess, "README.md"), "w") as fh:
            fh.write(f"# {title}\n")
        if conditions:
            with open(os.path.join(sess, "input", "conditions.csv"), "w") as fh:
                fh.write("run,condition\na,ctrl\nb,treated\n")
        if qc_marker is not None:
            with open(os.path.join(sess, "session.json"), "w") as fh:
                json.dump({"qc": qc_marker}, fh)
        return search_dir(os.path.join(sess, "output"), "search")

    def _stage(self, d, out, **kw):
        with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"), FRAN_HEALTH="off"):
            return run_quiet(fd.stage, Args(out, **kw))[1]

    def test_gabrigs_run_is_caught_once_its_name_is_known_and_withdrawn(self):
        """The real case: session '2026-09-23_chkLUppm_HeLa50_2026' (no QC token, no conditions),
        staged by the job-end hook, then named '... Lumos QC' by the agent at step 7c."""
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "2026-09-23_chkLUppm_HeLa50_2026")
            first = self._stage(d, out)                         # the job-end hook: no --name
            self.assertTrue(first["staged"])
            self.assertEqual(first["qc_rule"][:6], "not QC")
            res = self._stage(d, out, name="chkLUppm_HeLa50_2026 Lumos QC")
            self.assertFalse(res["staged"])
            self.assertEqual(res["reason"], "qc_run")
            self.assertIn("--not-qc", res["detail"])
            self.assertEqual(res["withdrawn"], first["entry"])
            man = _load_json(os.path.join(first["entry"], fd.MANIFEST))
            self.assertEqual((man["qc"], man["exclude"]), (True, True))
            self.assertIn("Lumos QC", man["qc_rule"])
            self.assertEqual(man["output_dir"], os.path.realpath(out))    # nothing else lost
            rec = fd.read_receipt(out)
            self.assertEqual(rec["status"], "qc_run")
            self.assertIn("Lumos QC", rec["qc_rule"])
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"),
                          FRAN_INGEST_LOG_DIR=os.path.join(d, "nologs")):
                v = run_quiet(fd.verify, Args(out))[1]
            self.assertEqual(v["state"], "qc_excluded")

    def test_a_qc_run_named_up_front_is_never_staged(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "Exploris QC2")
            res = self._stage(d, out)
            self.assertEqual(res["reason"], "qc_run")
            self.assertFalse(os.path.exists(os.path.join(d, "incoming", fd.entry_name(out))))

    def test_a_hela_experiment_with_conditions_is_not_qc(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "HeLa_digest_timecourse", conditions=True)
            res = self._stage(d, out, name="HeLa digest timecourse")
            self.assertTrue(res["staged"], res.get("detail"))
            man = _load_json(os.path.join(res["entry"], fd.MANIFEST))
            self.assertIs(man["qc"], False)
            self.assertTrue(man["qc_rule"].startswith("not QC"))

    def test_an_explicit_marker_beats_the_name_both_ways(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "Plasma_liver2", qc_marker=True)
            q, why = fd.is_qc_run(out)
            self.assertTrue(q)
            self.assertIn("session metadata qc: true", why)
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "QC_run_01", qc_marker=False)
            self.assertFalse(fd.is_qc_run(out)[0])
            self.assertTrue(fd.is_qc_run(out, override=True)[0])      # the flag beats metadata

    def test_not_qc_stages_a_false_positive_with_user_override(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "hela_qc_2")
            self.assertEqual(self._stage(d, out)["reason"], "qc_run")
            res = self._stage(d, out, not_qc=True)
            self.assertTrue(res["staged"], res.get("detail"))
            man = _load_json(os.path.join(res["entry"], fd.MANIFEST))
            self.assertEqual((man["qc"], man["qc_rule"]), (False, "user override"))

    def test_qc_flag_excludes_a_normal_name(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "Plasma_liver2")
            res = self._stage(d, out, qc=True)
            self.assertEqual(res["reason"], "qc_run")
            self.assertEqual(res["qc_rule"], "QC run: excluded by policy (qc: true, user override)")
            self.assertEqual(fd.read_receipt(out)["status"], "qc_run")
            # ...and backfill honours the recorded --qc, whatever the name says
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"), FRAN_HEALTH="off"):
                rows = fd.plan_backfill([out], prefixes=(os.path.realpath(d) + "/",))
            self.assertEqual((rows[0]["decision"], rows[0]["reason"]), ("excluded", "qc_run"))

    def test_precedence_matches_frans_policy_exclusion(self):
        """qc:true > DEFAULT_EXCLUDES tree (beats --not-qc) > qc:false > name rule."""
        glendon = "/quobyte/proteomics-grp/brett/glendon/sweep/search_out"
        self.assertTrue(fd.is_qc_run(glendon)[0])
        q, why = fd.is_qc_run(glendon, override=False)            # --not-qc cannot rescue it
        self.assertTrue(q)
        self.assertEqual(why, "QC run: excluded by policy (output_dir is under "
                              "/quobyte/proteomics-grp/brett/glendon/ (DEFAULT_EXCLUDES))")
        self.assertTrue(fd.is_qc_run("/quobyte/proteomics-grp/brett/v1_smoke/x")[0])
        self.assertEqual(fd.is_qc_run("/data/lab/QC_run_01/search", override=False),
                         (False, "user override"))
        self.assertEqual(fd.is_qc_run("/data/lab/x/search", names=["hela_qc_2"])[1],
                         "QC run: excluded by policy (search_name 'hela_qc_2' matches QC_NAME_RE)")
        self.assertEqual(fd.FRAN_DEFAULT_EXCLUDES, (
            "/quobyte/proteomics-grp/STAN/", "/quobyte/proteomics-grp/hela_qcs/",
            "/quobyte/proteomics-grp/brett/v1_smoke", "/quobyte/proteomics-grp/brett/glendon/",
            "/Data/lab/ToFEvoQC/"))

    def test_session_metadata_qc_false_does_not_beat_an_excluded_tree(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "Plasma_liver2", qc_marker=False)
            self.assertFalse(fd.is_qc_run(out)[0])
            real = fd.FRAN_DEFAULT_EXCLUDES
            fd.FRAN_DEFAULT_EXCLUDES = real + (os.path.realpath(d) + "/",)
            try:
                self.assertTrue(fd.is_qc_run(out)[0])
            finally:
                fd.FRAN_DEFAULT_EXCLUDES = real

    def test_frans_qc_trees_and_the_out_dir_path(self):
        self.assertTrue(fd.is_qc_run("/quobyte/proteomics-grp/hela_qcs/2026/x/search_out")[0])
        self.assertTrue(fd.is_qc_run("/nfs/lssc0/flinders/proteomics/Data/lab/ToFEvoQC/a/b")[0])
        self.assertTrue(fd.is_qc_run("/data/lab/QC_run_01/search")[0])
        self.assertFalse(fd.is_qc_run("/data/QC_run_01/lab/sub/search")[0])   # beyond 3 levels
        self.assertFalse(fd.is_qc_run("/data/lab/aqc_buffer_study/search")[0])

    def test_the_hook_decides_at_generation_so_nothing_is_ever_staged(self):
        """The race: the job-end hook stages at search end, the cron can ingest within 4 h, and
        the agent's later --name used to be the first time "QC" was seen. With the name baked
        into the hook's argv, the QC run never reaches the drop dir at all."""
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "2026-09-23_chkLUppm_HeLa50_2026")
            argv = fd.stage_argv(out, name="chkLUppm_HeLa50_2026 Lumos QC")
            self.assertEqual(argv[2:], ["stage", "--out", out, "--name",
                                        "chkLUppm_HeLa50_2026 Lumos QC"])
            a = Args(out, name=argv[argv.index("--name") + 1])
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"), FRAN_HEALTH="off"):
                res = run_quiet(fd.stage, a)[1]
            self.assertEqual(res["reason"], "qc_run")
            self.assertFalse(os.path.exists(os.path.join(d, "incoming", fd.entry_name(out))))
            self.assertEqual(fd.stage_argv(out, qc=True)[-1], "--qc")
            self.assertEqual(fd.stage_argv(out, qc=False)[-1], "--not-qc")
            self.assertNotIn("--name", fd.stage_argv(out, name="  "))

    def test_a_recorded_decision_outlives_a_later_flagless_stage(self):
        """The hook staged with --not-qc; the agent's later stage names it "... QC" and carries no
        flag. The user's word stands: no withdrawal. And a recorded --qc stays QC."""
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "Plasma_liver2")
            self.assertTrue(self._stage(d, out, not_qc=True)["staged"])
            res = self._stage(d, out, name="Plasma liver QC2", force=True)
            self.assertTrue(res["staged"], res.get("detail"))
            self.assertEqual(res["qc_rule"], "user override")
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "Plasma_liver2")
            self.assertEqual(self._stage(d, out, qc=True)["reason"], "qc_run")
            self.assertEqual(self._stage(d, out, name="Plasma liver 2")["reason"], "qc_run")

    def test_backfill_lists_qc_runs_separately_and_withdraws_on_apply(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            qc_out = search_dir(os.path.join(t.svc_root, "on_campus", "LabQ", "QC_run_01"),
                                "search_out")
            staged_qc = self._session(os.path.join(t.svc_root, "on_campus"), "2026-09-23_HeLa50")
            with env_vars(FRAN_DROP_DIR=t.drop, FRAN_HEALTH="off"):
                run_quiet(fd.stage, Args(staged_qc, name="chkLUppm_HeLa50_2026 Lumos QC",
                                         not_qc=True))          # staged, override recorded
                entry_qc2 = os.path.join(t.drop, fd.entry_name(qc_out))
                run_quiet(fd.stage, Args(qc_out, not_qc=True))
                # now take the override back out of the record, as if it were a plain stage
                rec = fd.read_receipt(qc_out)
                rec["qc_rule"] = "not QC: legacy"
                with open(os.path.join(qc_out, fd.RECEIPT), "w") as fh:
                    json.dump(rec, fh)
                man = _load_json(os.path.join(entry_qc2, fd.MANIFEST))
                man.pop("qc"), man.pop("qc_rule")
                with open(os.path.join(entry_qc2, fd.MANIFEST), "w") as fh:
                    json.dump(man, fh)
                found, _, _ = fd.discover([t.svc_root])
                pre = os.path.realpath(t.root) + "/"
                dry = {r["out"]: r for r in fd.plan_backfill(found, prefixes=(pre,))}
                self.assertEqual(dry[qc_out]["decision"], "excluded")
                self.assertEqual(dry[qc_out]["reason"], "qc_run")
                self.assertEqual(dry[qc_out]["would_withdraw"], entry_qc2)
                # the recorded --not-qc is honoured: not excluded, just already staged
                self.assertEqual(dry[staged_qc]["reason"], "already_staged")
                self.assertEqual(dry[t.skill]["decision"], "would_stage")
                fd.plan_backfill(found, apply=True, prefixes=(pre,))
            man = _load_json(os.path.join(entry_qc2, fd.MANIFEST))
            self.assertEqual((man["qc"], man["exclude"]), (True, True))
            self.assertEqual(fd.read_receipt(qc_out)["status"], "qc_run")


class VerifyFromLogsTests(unittest.TestCase):
    def _staged(self, d):
        out = search_dir(d)
        e = Env(os.path.join(d, "hive"), now=time.time())    # fd.verify reads the real clock
        with env_vars(FRAN_DROP_DIR=e.drop, FRAN_HEALTH="off"):
            run_quiet(fd.stage, Args(out))
        return out, e, os.path.join(e.drop, fd.entry_name(out))

    def _verify(self, out, e):
        with env_vars(FRAN_DROP_DIR=e.drop, FRAN_INGEST_LOG_DIR=e.logs):
            return run_quiet(fd.verify, Args(out))[1]

    def test_ingested_per_the_cron_log(self):
        with tempfile.TemporaryDirectory() as d:
            out, e, entry = self._staged(d)
            e.log(run_log(e.now - H, ingested=1, items=[("diann", entry, "ok")]))
            r = self._verify(out, e)
            self.assertEqual(r["state"], "ingested", r)
            self.assertTrue(r["ingested"])
            self.assertEqual(fd.read_receipt(out)["status"], "ingested")

    def test_failed_per_the_cron_log(self):
        with tempfile.TemporaryDirectory() as d:
            out, e, entry = self._staged(d)
            e.log(run_log(e.now - H, failed=1, items=[("diann", entry, "fail")]))
            r = self._verify(out, e)
            self.assertEqual(r["state"], "ingest_failed", r)
            self.assertIn("No precursor records parsed", r["detail"])
            self.assertEqual(fd.read_receipt(out)["status"], "staged")   # still handed over

    def test_pending_says_when_the_cron_is_stuck(self):
        with tempfile.TemporaryDirectory() as d:
            out, e, entry = self._staged(d)
            for k in (3, 2, 1):
                e.log(run_log(e.now - k * 4 * H, dup=3, failed=2, queued=181))
            e.submit(age_h=0)
            r = self._verify(out, e)
            self.assertEqual(r["state"], "staged_pending_cron")
            self.assertEqual(r["cron"]["verdict"], "stuck")
            self.assertIn("BUT FRAN's cron is stuck", r["detail"])


class AlertTests(unittest.TestCase):
    """health --alert goes through notify_slack.send_alert, the notifier's API for other scripts;
    the webhook never passes through this file."""
    STUCK = {"verdict": "stale_code", "summary": "FRAN ingest STUCK: 41 runs; ingest code STALE",
             "ingest_code": {"verdict": "stale",
                             "detail": "corpus_ingest.py stale: matches abc1234 (41 runs)"}}

    def _fake(self, sent=True):
        mod = type(sys)("notify_slack")
        mod.calls = []

        def send_alert(text, *, title=None, with_status=False):
            mod.calls.append((text, title))
            return (sent, "sent" if sent else "not sent: no webhook configured") \
                if with_status else sent
        mod.send_alert = send_alert
        return mod

    def test_a_stuck_cron_alone_never_pages(self):
        """FRAN's runner alerts on its own progress; a second page from here is a duplicate."""
        mod = self._fake()
        sys.modules["notify_slack"] = mod
        try:
            r = fd.send_alert({"verdict": "stuck", "summary": "FRAN ingest STUCK: 41 runs",
                               "ingest_code": {"verdict": "current"}}, cache_path="/nonexistent/c")
        finally:
            del sys.modules["notify_slack"]
        self.assertFalse(r["sent"])
        self.assertEqual(mod.calls, [])
        self.assertIsNone(fd.alert_reason({"verdict": "not_running", "ingest_code": {}}))

    def test_stale_or_modified_code_pages_with_the_code_detail(self):
        for verdict in ("stale", "modified"):
            h = {"verdict": "healthy", "summary": "x",
                 "ingest_code": {"verdict": verdict, "detail": "raw_metadata.py local modification"}}
            self.assertIn("raw_metadata.py", fd.alert_reason(h))
            self.assertIn(verdict, fd.alert_reason(h))

    def test_sends_once_then_holds_the_same_alert_for_a_day(self):
        with tempfile.TemporaryDirectory() as d:
            mod = self._fake()
            sys.modules["notify_slack"] = mod
            try:
                cache = os.path.join(d, "c.json")
                r1 = fd.send_alert(self.STUCK, cache_path=cache, now=NOW)
                r2 = fd.send_alert(self.STUCK, cache_path=cache, now=NOW + 3600)
                r3 = fd.send_alert(self.STUCK, cache_path=cache, now=NOW + 25 * 3600)
            finally:
                del sys.modules["notify_slack"]
            self.assertEqual((r1["sent"], r2["sent"], r3["sent"]), (True, False, True))
            self.assertEqual(len(mod.calls), 2)
            self.assertIn("41 runs", mod.calls[0][0])

    def test_healthy_is_never_posted(self):
        mod = self._fake()
        sys.modules["notify_slack"] = mod
        try:
            r = fd.send_alert({"verdict": "healthy", "summary": "ok",
                               "ingest_code": {"verdict": "current"}}, cache_path="/nonexistent/c")
        finally:
            del sys.modules["notify_slack"]
        self.assertFalse(r["sent"])
        self.assertEqual(mod.calls, [])

    def test_no_notifier_installed_is_a_no_op(self):
        from unittest import mock
        with mock.patch.dict(sys.modules, {"notify_slack": None}):     # import now raises
            r = fd.send_alert(self.STUCK, cache_path="/nonexistent/c")
        self.assertEqual(r, {"sent": False, "why": "notify_slack.py is not installed"})

    def test_a_webhook_that_is_not_configured_is_reported_not_retried_forever(self):
        with tempfile.TemporaryDirectory() as d:
            sys.modules["notify_slack"] = self._fake(sent=False)
            try:
                r = fd.send_alert(self.STUCK, cache_path=os.path.join(d, "c.json"), now=NOW)
            finally:
                del sys.modules["notify_slack"]
            self.assertFalse(r["sent"])
            self.assertIn("no webhook", r["why"])


# ------------------------------------------------------------------------------ meta fix --
class MetaSelectionTests(unittest.TestCase):
    def _meta(self, root, stem, organism, taxid):
        fa = os.path.join(root, f"{stem}.fasta")
        with open(fa, "w") as fh:
            fh.write(">sp|P1|A\nPEPTIDEK\n")
        with open(fa + ".meta.json", "w") as fh:
            json.dump({"fasta": fa, "organism": organism, "taxid": taxid, "n_sequences": 1}, fh)
        return fa

    def test_the_search_fasta_picks_the_sidecar_not_the_alphabet(self):
        """PROT_0793: human_, mouse_ and mouse_mousecont sidecars side by side; sorted() picks
        human for the mouse search. The search's own --fasta decides."""
        with tempfile.TemporaryDirectory() as d:
            self._meta(d, "human_UP000005640", "Homo sapiens", 9606)
            mfa = self._meta(d, "mouse_UP000000589_mousecont", "Mus musculus", 10090)
            out = search_dir(d, "search_mouse_mousecont", prov=False)
            with open(os.path.join(out, "report.log.txt"), "w") as fh:
                fh.write(f"diann-linux --f a.d --f b.d --lib x.parquet --fasta {mfa} --threads 8\n")
            self.assertEqual(fd.organism_from_meta(out)[:2], ("Mus musculus", 10090))
            self.assertEqual(fd.fasta_from_meta(out)[0], mfa)

    def test_several_sidecars_and_no_known_fasta_is_absent_not_guessed(self):
        with tempfile.TemporaryDirectory() as d:
            self._meta(d, "human", "Homo sapiens", 9606)
            self._meta(d, "mouse", "Mus musculus", 10090)
            out = search_dir(d, "s", prov=False)
            self.assertEqual(fd.organism_from_meta(out), (None, None, None))

    def test_a_lone_sidecar_is_still_used(self):
        with tempfile.TemporaryDirectory() as d:
            self._meta(d, "dog", "Canis lupus familiaris", 9615)
            out = search_dir(d, "search_out", prov=False)
            self.assertEqual(fd.organism_from_meta(out)[0], "Canis lupus familiaris")

    def test_session_layout_is_found_through_provenance(self):
        """<session>/output/search with the sidecar in <session>/input -- the old globs looked in
        <session>/output/input and missed it, so the gabrig search reached FRAN with no database."""
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "sess", "input")
            os.makedirs(inp)
            fa = self._meta(inp, "search", "Homo sapiens", 9606)
            out = search_dir(os.path.join(d, "sess", "output"), "search", prov=False)
            with open(os.path.join(out, "search_provenance.json"), "w") as fh:
                json.dump({"engine": "diann", "fasta": fa}, fh)
            self.assertEqual(fd.organism_from_meta(out)[0], "Homo sapiens")
            self.assertEqual(fd.fasta_from_meta(out)[0], fa)

    def test_a_staging_path_in_the_sidecar_becomes_the_path_the_search_read(self):
        """Silva / gabrig on 2026-09-24: the sidecar named ~/proteomics-pipeline/staging/search.fasta
        (it EXISTS, and the next session overwrites it); the search read <session>/input/search.fasta."""
        with tempfile.TemporaryDirectory() as d:
            staging = os.path.join(d, "home", "staging")
            os.makedirs(staging)
            with open(os.path.join(staging, "search.fasta"), "w") as fh:
                fh.write(">a\nPEPTIDEK\n")
            inp = os.path.join(d, "sess", "input")
            os.makedirs(inp)
            real = os.path.join(inp, "search.fasta")
            with open(real, "w") as fh:
                fh.write(">a\nPEPTIDEK\n")
            with open(real + ".meta.json", "w") as fh:
                json.dump({"fasta": os.path.join(staging, "search.fasta"), "organism": "Mus musculus",
                           "taxid": 10090, "md5": "abc", "n_sequences": 1}, fh)
            out = search_dir(os.path.join(d, "sess", "output"), "search", prov=False)
            with open(os.path.join(out, "search_provenance.json"), "w") as fh:
                json.dump({"engine": "diann", "fasta": real}, fh)
            self.assertEqual(fd.fasta_from_meta(out), (real, "abc", 1))

    def test_a_laptop_path_in_the_sidecar_becomes_the_path_the_search_read(self):
        """hive_remote writes the sidecar on the laptop, so its `fasta` is /Users/...; the search
        on HIVE read the same file name from its own input dir."""
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "input")
            os.makedirs(inp)
            real = os.path.join(inp, "dog_contam.fasta")
            with open(real, "w") as fh:
                fh.write(">a\nPEPTIDEK\n")
            with open(real + ".meta.json", "w") as fh:
                json.dump({"fasta": "/Users/someone/sessions/x/input/dog_contam.fasta",
                           "organism": "Canis lupus familiaris", "taxid": 9615,
                           "md5": "abc", "n_sequences": 1}, fh)
            out = search_dir(d, "search_out", prov=False)
            with open(os.path.join(out, "report.log.txt"), "w") as fh:
                fh.write(f"diann-linux --fasta {real} --out report.parquet\n")
            self.assertEqual(fd.fasta_from_meta(out), (real, "abc", 1))


class ReviewFixTests(unittest.TestCase):
    """The independent review's FIX-FIRST list (items A-I). Each test is the reproduction."""

    def _session(self, d, title, **kw):
        return QcRuleTests._session(self, d, title, **kw)

    def _stage(self, d, out, **kw):
        with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"), FRAN_HEALTH="off"):
            return run_quiet(fd.stage, Args(out, **kw))

    # B ----------------------------------------------------------------------------------------
    def test_B_a_withdrawal_survives_a_later_plain_stage(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "2026-09-23_chkLUppm_HeLa50_2026")
            entry = self._stage(d, out)[1]["entry"]                       # the hook: no name
            w = self._stage(d, out, name="chkLUppm_HeLa50_2026 Lumos QC")[1]
            self.assertEqual(w["withdrawn"], entry)
            self.assertEqual(fd.read_receipt(out)["search_name"], "chkLUppm_HeLa50_2026 Lumos QC")
            for kw in ({}, {"force": True}):                             # plain stage, twice
                r = self._stage(d, out, **kw)[1]
                self.assertEqual(r["reason"], "qc_run", r.get("detail"))
                man = _load_json(os.path.join(entry, fd.MANIFEST))
                self.assertEqual((man["qc"], man["exclude"]), (True, True))
            # only an EXPLICIT --not-qc turns it back
            r = self._stage(d, out, not_qc=True)[1]
            self.assertTrue(r["staged"], r.get("detail"))
            self.assertEqual(r["qc_rule"], "user override")

    # A ----------------------------------------------------------------------------------------
    def test_A_verify_never_invents_a_staged_receipt(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"),
                          FRAN_INGEST_LOG_DIR=os.path.join(d, "nologs")):
                v = run_quiet(fd.verify, Args(out))[1]
                self.assertEqual(v["state"], "not_staged")
                self.assertIsNone(v["receipt"])
                self.assertIsNone(fd.read_receipt(out))
                self.assertTrue(fd.check(Args(out))["eligible"])       # not already_staged

    def test_A_verify_of_a_staged_entry_without_a_receipt_records_staged(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            self._stage(d, out)
            os.unlink(os.path.join(out, fd.RECEIPT))
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"),
                          FRAN_INGEST_LOG_DIR=os.path.join(d, "nologs")):
                run_quiet(fd.verify, Args(out))
            self.assertEqual(fd.read_receipt(out)["status"], "staged")

    # C ----------------------------------------------------------------------------------------
    def test_C_the_name_rule_reads_the_real_path_like_fran(self):
        """FRAN judges output_dir = realpath(out). A symlinked search dir whose real home is a
        QC path must be QC here too."""
        with tempfile.TemporaryDirectory() as d:
            real = search_dir(os.path.join(d, "instr", "QC_run_01"), "search_out")
            link_parent = os.path.join(d, "lab", "proj")
            os.makedirs(link_parent)
            link = os.path.join(link_parent, "search_out")
            os.symlink(real, link)
            q, why = fd.is_qc_run(link)
            self.assertTrue(q, why)
            self.assertIn("QC_run_01", why)

    def test_C_the_smb_spelling_is_mapped_to_hive(self):
        q, why = fd.is_qc_run("/Volumes/proteomics-grp/brett/glendon/sweep/search_out",
                              override=False)
        self.assertTrue(q)
        self.assertIn("/quobyte/proteomics-grp/brett/glendon/", why)
        self.assertEqual(fd._excludes_path("/Volumes/proteomics-grp/brett/glendon"),
                         "/quobyte/proteomics-grp/brett/glendon/")

    # D ----------------------------------------------------------------------------------------
    def test_D_a_022_umask_still_leaves_the_entry_group_writable(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            os.makedirs(os.path.join(out, "report_xic"))
            with open(os.path.join(out, "report_xic", "a.xic.parquet"), "w") as fh:
                fh.write("x")
            old = os.umask(0o022)
            try:
                entry = self._stage(d, out)[1]["entry"]
            finally:
                os.umask(old)
            mode = lambda p: os.stat(p).st_mode & 0o777               # noqa: E731
            self.assertEqual(mode(entry), 0o775)
            self.assertEqual(mode(os.path.join(entry, "report_xic")), 0o775)
            self.assertEqual(mode(os.path.join(entry, fd.MANIFEST)), 0o664)
            self.assertEqual(mode(os.path.join(out, fd.RECEIPT)), 0o664)

    def test_D_a_failed_withdrawal_says_so_and_never_claims_success(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "2026-09-23_HeLa50")
            entry = self._stage(d, out)[1]["entry"]
            os.chmod(entry, 0o555)                                        # not ours to rewrite
            try:
                r = self._stage(d, out, qc=True)[1]
            finally:
                os.chmod(entry, 0o775)
            self.assertEqual(r["reason"], "qc_run")
            self.assertTrue(r["detail"].startswith("withdraw FAILED:"), r["detail"])
            self.assertNotIn("now marked", r["detail"])
            self.assertIsNone(r.get("withdrawn"))
            self.assertIn("withdraw_failed", fd.read_receipt(out))
            self.assertIs(_load_json(os.path.join(entry, fd.MANIFEST))["qc"], False)

    def test_D_restaging_an_unwritable_entry_is_reported_not_a_crash(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            entry = self._stage(d, out)[1]["entry"]
            os.chmod(entry, 0o555)
            try:
                code, r, err = self._stage(d, out, force=True)
            finally:
                os.chmod(entry, 0o775)
            self.assertEqual(code, 0, err)
            self.assertEqual(r["reason"], "entry_not_writable")

    def test_D_a_permission_error_is_never_recorded_as_not_core(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            drop = os.path.join(d, "incoming")
            os.makedirs(drop)
            os.chmod(drop, 0o555)
            saved = fd._in_core_group
            try:
                for member, reason, recorded in ((None, "not_core_facility", False),
                                                 (True, "drop_dir_not_writable", False),
                                                 (False, "not_core_facility", True)):
                    fd._in_core_group = lambda m=member: m
                    r = self._stage(d, out)[1]
                    self.assertEqual(r["reason"], reason, member)
                    rec = fd.read_receipt(out)
                    self.assertEqual(bool(rec and rec.get("status") == "not_core_facility"),
                                     recorded, member)
            finally:
                fd._in_core_group = saved
                os.chmod(drop, 0o755)

    def test_D_an_entry_that_is_a_symlink_is_replaced_never_followed(self):
        """Re-review, DESTRUCTIVE: incoming/<name>__<hash8> as a symlink to the search dir made
        the relink loop unlink report.parquet / report.log.txt INSIDE the real search dir, and
        chmod it to 2775."""
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            with open(os.path.join(out, "report.log.txt"), "w") as fh:
                fh.write("diann log\n")
            os.chmod(out, 0o750)
            drop = os.path.join(d, "incoming")
            os.makedirs(drop)
            entry = os.path.join(drop, fd.entry_name(out))
            os.symlink(out, entry)                                    # a legacy bare-link entry
            before = {f: (os.stat(os.path.join(out, f)).st_ino,
                          _read(os.path.join(out, f), "rb"))
                      for f in ("report.parquet", "report.log.txt")}
            mode = os.stat(out).st_mode & 0o7777
            r = self._stage(d, out)[1]
            self.assertTrue(r["staged"], r.get("detail"))
            for f, (ino, body) in before.items():
                p = os.path.join(out, f)
                self.assertFalse(os.path.islink(p), f)
                self.assertEqual((os.stat(p).st_ino, _read(p, "rb")), (ino, body), f)
            self.assertEqual(os.stat(out).st_mode & 0o7777, mode)
            self.assertFalse(os.path.islink(entry))
            self.assertTrue(os.path.isdir(entry))
            self.assertEqual(os.path.realpath(os.path.join(entry, "report.parquet")),
                             os.path.realpath(os.path.join(out, "report.parquet")))
            # and a withdrawal refuses to write through a symlink entry
            os.rename(entry, entry + ".real")
            os.symlink(out, entry)
            done, why = fd._withdraw_qc_entry(entry, "x", "me")
            self.assertIsNone(done)
            self.assertIn("symlink", why)
            self.assertFalse(os.path.exists(os.path.join(out, fd.MANIFEST)))

    def test_D_entry_not_writable_is_never_listed_as_would_stage(self):
        """do_stage's entry_not_writable result carried eligible: True, so backfill --apply listed
        a search it had just failed to stage as would_stage."""
        with tempfile.TemporaryDirectory() as d:
            out = self._session(d, "Plasma_liver2")
            drop = os.path.join(d, "incoming")
            os.makedirs(drop)
            with open(os.path.join(drop, fd.entry_name(out)), "w") as fh:
                fh.write("not a directory")                      # the entry path is occupied
            with env_vars(FRAN_DROP_DIR=drop, FRAN_HEALTH="off"):
                res = fd.do_stage(Args(out))
                self.assertEqual((res["eligible"], res["staged"], res["reason"]),
                                 (False, False, "entry_not_writable"))
                row = fd.plan_backfill([out], apply=True, prefixes=(os.path.realpath(d) + "/",))[0]
            self.assertEqual((row["decision"], row["reason"]), ("skip", "entry_not_writable"))

    # E ----------------------------------------------------------------------------------------
    def test_E_an_explicit_meta_for_another_database_is_ignored(self):
        with tempfile.TemporaryDirectory() as d:
            mouse = os.path.join(d, "mouse.fasta")
            human = os.path.join(d, "human.fasta")
            for fa, org, tx in ((mouse, "Mus musculus", 10090), (human, "Homo sapiens", 9606)):
                with open(fa, "w") as fh:
                    fh.write(">a\nPEPTIDEK\n")
                with open(fa + ".meta.json", "w") as fh:
                    json.dump({"fasta": fa, "organism": org, "taxid": tx, "n_sequences": 1}, fh)
            out = search_dir(os.path.join(d, "s"), "search_out", prov=False)
            with open(os.path.join(out, "report.log.txt"), "w") as fh:
                fh.write(f"diann-linux --f a.d --fasta {mouse} --out report.parquet\n")
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming")):
                err = io.StringIO()
                with contextlib.redirect_stderr(err):
                    c = fd.check(Args(out, fasta_meta=human + ".meta.json"))
            self.assertIn("--fasta-meta ignored", err.getvalue())
            self.assertIn("fasta_meta_ignored", c)
            # the wrong meta is dropped; the sidecar tied to the search's own --fasta decides
            self.assertEqual((c["organism"], c["fasta_path"]), ("Mus musculus", mouse))
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming")):
                c = fd.check(Args(out, fasta_meta=mouse + ".meta.json"))
            self.assertEqual((c["organism"], c["fasta_path"]), ("Mus musculus", mouse))
            self.assertNotIn("fasta_meta_ignored", c)
            # ...and with no sidecar of the search's own FASTA, blank -- never the wrong one
            os.unlink(mouse + ".meta.json")
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming")):
                c = fd.check(Args(out, fasta_meta=human + ".meta.json"))
            self.assertEqual((c["organism"], c["fasta_path"]), (None, None))

    # F ----------------------------------------------------------------------------------------
    def test_F_coursework_is_never_staged_by_stage_or_backfill(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d)
            saved = fd.getpass.getuser
            fd.getpass.getuser = lambda: "proteomics-class-07"
            try:
                r = self._stage(d, out)[1]
            finally:
                fd.getpass.getuser = saved
            self.assertEqual(r["reason"], "not_core_facility")
            self.assertIn("teaching account proteomics-class-07", r["detail"])
            self.assertEqual(fd.read_receipt(out)["status"], "not_core_facility")
            saved_owner = fd._owner
            fd._owner = lambda p: "proteomics-class-03"
            try:
                core, why = fd.core_search(out, {"brettsp"}, prefixes=(os.path.realpath(d) + "/",))
            finally:
                fd._owner = saved_owner
            self.assertFalse(core)
            self.assertIn("coursework", why)

    # G ----------------------------------------------------------------------------------------
    def _fragpipe(self, d, log_tail):
        out = os.path.join(d, "fp_out")
        os.makedirs(os.path.join(out, "dia-quant-output"))
        with open(os.path.join(out, "dia-quant-output", "report.tsv"), "w") as fh:
            fh.write("x" * 64)                  # present in every FragPipe run, finished or not
        with open(os.path.join(out, "search_provenance.json"), "w") as fh:
            json.dump({"engine": "fragpipe"}, fh)
        if log_tail is not None:
            with open(os.path.join(out, "log_2026-06-10_22-23-44.txt"), "w") as fh:
                fh.write("DIA-Quant run DIA-NN: 4.91 minutes\n" + log_tail)
        return out

    def test_G_a_cancelled_fragpipe_run_is_incomplete_despite_its_report(self):
        """Real on HIVE: glendon/fragpipe_bigdog_dia_FULL215 was cancelled and still has
        dia-quant-output/report.tsv."""
        with tempfile.TemporaryDirectory() as d:
            out = self._fragpipe(d, "~~~~~~~~~~~~~~~~~~~~\nCancelling 5 remaining tasks\n")
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming")):
                c = fd.check(Args(out))
            self.assertEqual(c["reason"], "search_incomplete")
            self.assertIn("cancelled", c["detail"])
        with tempfile.TemporaryDirectory() as d:
            out = self._fragpipe(d, "=====ALL JOBS DONE IN 9.5 MINUTES=====\n")
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming")):
                self.assertTrue(fd.check(Args(out))["eligible"])

    def test_G_fulcrum_without_success_is_incomplete(self):
        with tempfile.TemporaryDirectory() as d:
            out = search_dir(d, "rad", engine="radiant")
            fr = os.path.join(out, "radiant_results", "fulcrum-results")
            os.makedirs(fr)
            with open(os.path.join(fr, "part-00000.snappy.parquet"), "w") as fh:
                fh.write("x" * 64)
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming")):
                self.assertEqual(fd.check(Args(out))["reason"], "search_incomplete")
                open(os.path.join(fr, "_SUCCESS"), "w").close()
                self.assertTrue(fd.check(Args(out))["eligible"])

    def test_G_no_marker_the_agent_may_stage_but_backfill_never_does(self):
        with tempfile.TemporaryDirectory() as d:
            out = self._fragpipe(d, None)                                 # no FragPipe log
            with env_vars(FRAN_DROP_DIR=os.path.join(d, "incoming"), FRAN_HEALTH="off"):
                c = fd.check(Args(out))
                self.assertTrue(c["eligible"])                            # the agent watched it
                self.assertIsNone(c["completion"]["done"])
                for apply in (False, True):
                    row = fd.plan_backfill([out], apply=apply,
                                           prefixes=(os.path.realpath(d) + "/",))[0]
                    self.assertEqual((row["decision"], row["reason"]),
                                     ("skip", "needs_agent_check"), apply)
            self.assertFalse(os.path.exists(os.path.join(d, "incoming", fd.entry_name(out))))

    # I ----------------------------------------------------------------------------------------
    def test_I_detection_markers_match_frans_engine_markers(self):
        self.assertEqual(fd.DETECT_MARKERS, [
            ("fragpipe", ("dia-quant-output/report.tsv", "fragpipe.fp-manifest")),
            ("radiant", ("radiant_results/fulcrum-results", "fulcrum-results/_SUCCESS")),
            ("diann", ("report.parquet", "report.tsv"))])
        with tempfile.TemporaryDirectory() as d:
            with open(os.path.join(d, "delimp_report.parquet"), "w") as fh:
                fh.write("x")
            self.assertIsNone(fd.detect_engine(d)[0])


class CorpusQueryGuardTests(unittest.TestCase):
    """No test may ever query the live FRAN database (FRAN's own suite once did)."""

    def test_the_switch_returns_before_the_token_is_even_looked_at(self):
        import builtins
        import subprocess
        tokens = {os.path.expanduser(t) for t in fd.TOKEN_CANDIDATES if t}

        def guard(real):
            def f(path, *a, **kw):
                if os.path.expanduser(str(path)) in tokens:
                    raise AssertionError(f"touched the PG Farm token path {path}")
                return real(path, *a, **kw)
            return f

        def no_subprocess(*a, **kw):
            raise AssertionError("corpus_query started a subprocess")
        patches = [(builtins, "open"), (os.path, "exists"), (os.path, "isfile"), (os, "access"),
                   (subprocess, "run")]
        saved = [(m, n, getattr(m, n)) for m, n in patches]
        try:
            builtins.open, os.path.exists = guard(saved[0][2]), guard(saved[1][2])
            os.path.isfile, os.access = guard(saved[2][2]), guard(saved[3][2])
            subprocess.run = no_subprocess
            with env_vars(FRAN_CORPUS_QUERY="off"):
                self.assertIsNone(fd.corpus_query("/some/search_out", "/drop/entry"))
            # the guard is live: without the switch, the lookup reaches the token path
            with env_vars(FRAN_CORPUS_QUERY=None, FRAN_INGEST_DIR=tempfile.gettempdir(),
                          FRAN_INGEST_PYTHON=sys.executable):
                saved_dirs, saved_py = fd.INGEST_DIRS, fd.PY_CANDIDATES
                fd.INGEST_DIRS, fd.PY_CANDIDATES = [tempfile.gettempdir()], [sys.executable]
                try:
                    with self.assertRaises(AssertionError):
                        fd.corpus_query("/some/search_out", "/drop/entry")
                finally:
                    fd.INGEST_DIRS, fd.PY_CANDIDATES = saved_dirs, saved_py
        finally:
            for m, n, f in saved:
                setattr(m, n, f)

    def test_every_fran_test_module_switches_it_off(self):
        for mod in ("test_fran_deposit_gate.py", "test_fran_manifest_fasta.py",
                    "test_fran_health_backfill.py"):
            src = _read(os.path.join(HERE, mod))
            self.assertIn('os.environ["FRAN_CORPUS_QUERY"] = "off"', src, mod)
            self.assertIn("def setUpModule", src, mod)


class ReviewerReproTests(unittest.TestCase):
    """The independent reviewer's reproduction script (scratchpad/review_fd_repro/repro.sh, A-E),
    as unit tests: same steps, same inputs -- umask 022, a 'PAR1xxxxPAR1' report stub, sidecars in
    fetch_fasta.py's nested {"selected": {...}} form. Each reproduced on af9e46d and must not on
    the fixed code."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.T = self._tmp.name
        self.drop = os.path.join(self.T, "drop", "incoming")
        os.makedirs(self.drop)
        self._umask = os.umask(0o022)
        self._env = env_vars(FRAN_DROP_DIR=self.drop, FRAN_HEALTH="off", FRAN_DEPOSIT=None,
                             FRAN_INGEST_LOG_DIR=os.path.join(self.T, "logs"))
        self._env.__enter__()

    def tearDown(self):
        self._env.__exit__(None, None, None)
        os.umask(self._umask)
        self._tmp.cleanup()

    def mk(self, *parts):
        d = os.path.join(self.T, *parts)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "report.parquet"), "w") as fh:
            fh.write("PAR1xxxxPAR1")
        return d

    def cli(self, fn, out, **kw):
        return run_quiet(fn, Args(out, **kw))[1]

    def test_A_verify_then_stage(self):
        out = self.mk("A", "search_out")
        self.assertEqual(self.cli(fd.verify, out)["state"], "not_staged")
        self.assertFalse(os.path.exists(os.path.join(out, fd.RECEIPT)))
        r = self.cli(fd.stage, out)
        self.assertEqual((r["staged"], r["reason"]), (True, "ok"))

    def test_B_withdrawal_then_name_less_stage(self):
        out = self.mk("B", "job1", "search_out")
        entry = self.cli(fd.stage, out)["entry"]
        self.assertEqual(self.cli(fd.stage, out, name="Lumos QC")["withdrawn"], entry)
        r = self.cli(fd.stage, out)
        self.assertEqual((r["staged"], r["reason"]), (False, "qc_run"))
        m = _load_json(os.path.join(entry, fd.MANIFEST))
        self.assertEqual((m["qc"], m["exclude"]), (True, True))

    def test_C_symlinked_session_search_into_a_qc_path(self):
        real = self.mk("C", "real", "Lumos_QC_2026", "run1", "search_out")
        os.makedirs(os.path.join(self.T, "C", "sess", "output"))
        link = os.path.join(self.T, "C", "sess", "output", "search")
        os.symlink(real, link)
        r = self.cli(fd.stage, link)
        self.assertEqual((r["staged"], r["reason"]), (False, "qc_run"))
        self.assertIn("Lumos_QC_2026", r["qc_rule"])
        self.assertEqual(os.listdir(self.drop), [])

    def test_D_modes_under_umask_022(self):
        entry = self.cli(fd.stage, self.mk("D", "search_out"))["entry"]
        self.assertEqual(os.stat(entry).st_mode & 0o777, 0o775)
        self.assertEqual(os.stat(os.path.join(entry, fd.MANIFEST)).st_mode & 0o777, 0o664)

    def test_E_explicit_meta_of_another_search(self):
        P = os.path.join(self.T, "E", "P")
        sess_in = os.path.join(self.T, "E", "sess", "input")
        os.makedirs(sess_in)
        out = self.mk("E", "P", "search_mouse")
        for fa, org, tx in ((os.path.join(P, "mouse.fasta"), "Mus musculus", 10090),
                            (os.path.join(sess_in, "search.fasta"), "Homo sapiens", 9606)):
            with open(fa, "w") as fh:
                fh.write(">a\nPEPTIDEK\n")
            with open(fa + ".meta.json", "w") as fh:
                json.dump({"selected": {"fasta": fa, "organism": org, "taxid": tx,
                                        "n_sequences": 1}}, fh)
        with open(os.path.join(out, "report.log.txt"), "w") as fh:
            fh.write(f"diann-linux --f a.d --fasta {os.path.join(P, 'mouse.fasta')} --out x\n")
        code, r, err = run_quiet(fd.stage, Args(out, fasta_meta=os.path.join(
            sess_in, "search.fasta.meta.json")))
        self.assertTrue(r["staged"])
        mouse = os.path.join(P, "mouse.fasta")
        self.assertEqual((r["organism"], r["taxon"], r["fasta_path"]), ("Mus musculus", 10090, mouse))
        self.assertIn("the search read", r["fasta_meta_ignored"])
        self.assertIn("--fasta-meta ignored", err)
        m = _load_json(os.path.join(r["entry"], fd.MANIFEST))
        self.assertEqual((m["organism"], m["fasta_path"]), ("Mus musculus", mouse))


# ------------------------------------------------------------------------------ backfill --
class Tree:
    """A fake service tree: skill searches, an app search, a failed one, raw data."""
    def __init__(self, root):
        self.root = root
        svc = os.path.join(root, "service", "on_campus")
        self.skill = search_dir(os.path.join(svc, "LabA", "2026-09-21_Dog"), "search_out")
        self.session = search_dir(os.path.join(svc, "LabB", "2026", "sess", "output"), "search")
        # a 5-step chain run through hive_remote: no provenance on HIVE, only the chain's logs
        self.chain = search_dir(os.path.join(svc, "LabC", "run"), "search_out", prov=False,
                                extra=("step5_report.sbatch", "s5_report_23945276.log"))
        # the DE-LIMP app writes the same step names, but its jobs are diann_<name>_<step>
        self.app = search_dir(os.path.join(svc, "LabD", "app"), "out", prov=False,
                              extra=("step5_report.sbatch", "diann_proj_step5_1.log"))
        self.failed = search_dir(os.path.join(svc, "LabE"), "search_out", empty=True)
        self.sage = search_dir(os.path.join(svc, "LabF"), "search_out", engine="sage")
        raw = os.path.join(svc, "LabA", "raw", "run1.d")
        os.makedirs(raw)
        search_dir(raw, "inside_a_d")                  # must never be reached: .d is pruned
        self.drop = os.path.join(root, "incoming")
        os.makedirs(self.drop)
        self.svc_root = os.path.join(root, "service")


class DiscoveryTests(unittest.TestCase):
    def test_finds_search_dirs_and_prunes_raw_containers(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            found, _, stats = fd.discover([t.svc_root])
            self.assertEqual(set(found), {t.skill, t.session, t.chain, t.app, t.failed, t.sage})
            self.assertFalse(stats["truncated"])

    def test_depth_limit(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            found, _, _ = fd.discover([t.svc_root], max_depth=4)
            self.assertIn(t.skill, found)             # depth 4
            self.assertNotIn(t.session, found)        # depth 6

    def test_time_budget_truncates_instead_of_running_on(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            tick = iter(range(10 ** 6))
            found, _, stats = fd.discover([t.svc_root], budget_s=3, clock=lambda: next(tick))
            self.assertTrue(stats["truncated"])
            self.assertLess(stats["dirs_walked"], 6)

    def test_time_small_roots_do_not_use_goes_to_the_big_one(self):
        """Seven empty roots and one real tree, with a budget the tree needs almost all of: an
        equal one-shot split gave the tree 1/8 and truncated it (seen on HIVE)."""
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            empties = []
            for k in range(7):
                empties.append(os.path.join(d, f"home{k}"))
                os.makedirs(empties[-1])
            # the fake clock advances one tick per call; measure what a full walk of the tree costs
            tick = iter(range(10 ** 6))
            _, _, full = fd.discover([t.svc_root], budget_s=10 ** 5, clock=lambda: next(tick))
            cost = full["elapsed_s"]
            budget = cost + 100              # an equal one-shot 1/8 split would truncate this
            self.assertLess(budget / 8, cost)
            tick = iter(range(10 ** 6))
            _, _, stats = fd.discover(empties + [t.svc_root], budget_s=budget,
                                      clock=lambda: next(tick))
            self.assertFalse(stats["truncated"], stats)
            self.assertEqual(stats["per_root"][t.svc_root]["dirs_walked"], full["dirs_walked"])

    def test_symlinks_are_never_followed(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            os.symlink(os.path.dirname(t.skill), os.path.join(t.svc_root, "link_to_labA"))
            found, _, _ = fd.discover([t.svc_root])
            self.assertEqual(sum(1 for f in found if "link_to_labA" in f), 0)

    def test_a_session_record_points_at_a_search_outside_the_tree(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            sess = os.path.join(t.svc_root, "on_campus", "LabG", "sess")
            os.makedirs(sess)
            elsewhere = search_dir(os.path.join(d, "hive_out"), "search_out")
            with open(os.path.join(sess, ".recovery.json"), "w") as fh:
                json.dump({"report": os.path.join(elsewhere, "report.parquet")}, fh)
            _, sessions, _ = fd.discover([t.svc_root])
            self.assertEqual(sessions, [sess])
            self.assertEqual(fd.session_search_dirs(sess), [elsewhere])


class PlanTests(unittest.TestCase):
    def _plan(self, t, apply=False, prefixes=None, members=None):
        with env_vars(FRAN_DROP_DIR=t.drop, FRAN_HEALTH="off"):
            found, _, _ = fd.discover([t.svc_root])
            return {r["out"]: r for r in fd.plan_backfill(
                found, apply=apply, members=members,
                prefixes=(os.path.realpath(t.root) + "/",) if prefixes is None else prefixes)}

    def test_dry_run_classifies_with_the_existing_reason_codes(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            rows = self._plan(t)
            self.assertEqual(rows[t.skill]["decision"], "would_stage")
            self.assertEqual(rows[t.session]["decision"], "would_stage")
            self.assertEqual(rows[t.chain]["decision"], "would_stage")
            self.assertEqual(rows[t.chain]["marker"], "s5_report_23945276.log")
            self.assertEqual(rows[t.app]["reason"], "not_a_skill_search")
            self.assertEqual(rows[t.failed]["reason"], "search_incomplete")
            self.assertEqual(rows[t.sage]["reason"], "engine_unsupported")
            self.assertEqual(os.listdir(t.drop), [], "a dry run staged something")

    def test_a_dry_run_does_not_create_a_missing_drop_dir(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            os.rmdir(t.drop)
            rows = self._plan(t)
            self.assertEqual(rows[t.skill]["decision"], "would_stage")
            self.assertFalse(os.path.exists(t.drop))

    def test_already_staged_and_opted_out_are_skipped(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            os.makedirs(os.path.join(t.drop, fd.entry_name(t.skill)))            # entry, no receipt
            with open(os.path.join(t.session, fd.RECEIPT), "w") as fh:
                json.dump({"status": "opted_out", "decided_by": "gabrig", "at": "2026-09-01"}, fh)
            with open(os.path.join(t.chain, fd.RECEIPT), "w") as fh:
                json.dump({"status": "staged"}, fh)
            rows = self._plan(t)
            self.assertEqual(rows[t.skill]["reason"], "already_staged")
            self.assertEqual(rows[t.session]["reason"], "opted_out")
            self.assertEqual(rows[t.chain]["reason"], "already_staged")

    def test_a_search_outside_the_core_trees_owned_by_a_non_member_is_refused(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            rows = self._plan(t, prefixes=("/quobyte/proteomics-grp/",), members=set())
            self.assertEqual(rows[t.skill]["reason"], "not_core_facility")
            rows = self._plan(t, prefixes=("/quobyte/proteomics-grp/",),
                              members={__import__("getpass").getuser()})
            self.assertEqual(rows[t.skill]["decision"], "would_stage")

    def test_a_search_fran_already_ingested_by_another_route_is_not_restaged(self):
        """PROT_0793/search_mouse reached FRAN through its queue on 2026-09-08, by its real path,
        with no receipt. Staging it again only feeds the duplicate guard."""
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            log = run_log(NOW - H, ingested=1, items=[("diann", os.path.realpath(t.skill), "ok")])
            runs = [dict(fd.parse_ingest_log(log), log="L", mtime=NOW, jobid=1)]
            with env_vars(FRAN_DROP_DIR=t.drop, FRAN_HEALTH="off"):
                rows = {r["out"]: r for r in fd.plan_backfill(
                    [t.skill, t.session], members=None, runs=runs,
                    prefixes=(os.path.realpath(t.root) + "/",))}
            self.assertEqual(rows[t.skill]["reason"], "already_ingested")
            self.assertEqual(rows[t.session]["decision"], "would_stage")

    def test_apply_stages_through_the_same_path_with_a_labelled_name(self):
        with tempfile.TemporaryDirectory() as d:
            t = Tree(d)
            rows = self._plan(t, apply=True)
            self.assertEqual(rows[t.session]["decision"], "staged")
            entry = rows[t.session]["entry"]
            man = _load_json(os.path.join(entry, fd.MANIFEST))
            self.assertEqual(man["search_name"], "sess")          # not "search"
            self.assertIn("backfill", man["search_name_source"])
            self.assertEqual(man["output_dir"], os.path.realpath(t.session))
            self.assertEqual(fd.read_receipt(t.session)["status"], "staged")
            # a second pass converges: nothing new
            again = self._plan(t, apply=True)
            self.assertEqual(again[t.session]["reason"], "already_staged")


class LoginGuardTests(unittest.TestCase):
    class A:
        def __init__(self, roots=None, list=None, allow=False):
            self.roots, self.list, self.allow_login_node = roots, list, allow

    def test_a_walk_outside_slurm_is_refused(self):
        r = fd.login_guard(self.A(), 0, env={})
        self.assertEqual(r["refused"], "login_node")
        self.assertIn("--sbatch", r["detail"])

    def test_a_short_list_is_fine_on_a_login_node_and_a_long_one_is_not(self):
        self.assertIsNone(fd.login_guard(self.A(list="x"), 3, env={}))
        self.assertIsNotNone(fd.login_guard(self.A(list="x"), fd.LOGIN_LIST_MAX + 1, env={}))

    def test_inside_slurm_or_forced_it_runs(self):
        self.assertIsNone(fd.login_guard(self.A(), 0, env={"SLURM_JOB_ID": "1"}))
        self.assertIsNone(fd.login_guard(self.A(allow=True), 0, env={}))


class SbatchTests(unittest.TestCase):
    def test_writes_a_dry_run_job_with_the_given_queue(self):
        with tempfile.TemporaryDirectory() as d:
            class A:
                sbatch_dir, time_budget, max_depth = d, 1200, 9
                roots, list, no_homes, apply = None, None, False, False
                partition, account, qos = "high", "genome-center-grp", "genome-center-grp-high-qos"
                sbatch_minutes = 30
            _, res, _ = run_quiet(lambda: fd.jout(fd.write_backfill_sbatch(A())))
            body = _read(res["sbatch"])
            self.assertIn("#SBATCH --partition=high", body)
            self.assertIn("#SBATCH --account=genome-center-grp", body)
            self.assertIn("#SBATCH --time=30", body)
            self.assertIn("#SBATCH --mem=2G", body)
            self.assertIn("fran_deposit.py backfill", body)
            self.assertIn('/"fran_backfill_${SLURM_JOB_ID}.json"', body)   # bash expands the id
            self.assertNotIn("--apply", body)
            self.assertIn("DRY RUN", body)
            self.assertEqual(res["mode"], "dry_run")


if __name__ == "__main__":
    unittest.main(verbosity=2)
