#!/usr/bin/env python3
"""save_transcript.py, transcript_hook.py and log_decision.py -- the record of how an analysis
was done, for an AI to review later (Brett, 2026-09-28: "save a transcript of the data analysis
session for an AI to look over later. Just as another reproducibility check").

Guards (stdlib; synthetic transcripts in the shapes Claude Code 2.1.x writes -- never a real
one, which would carry other projects and internal paths):
  * the transcript is found by its session id under $CLAUDE_CONFIG_DIR or ~/.claude/projects;
  * every kind of secret is redacted before anything is written: a planted fake of each is in
    no saved file, and the JSONL stays JSON;
  * no Claude Code, or no transcript: a clean skip, exit 0, nothing written;
  * an entry of an unknown shape is one "[unrecognised entry]" line, never a crash;
  * two conversations: two files and one index; saving one again changes nothing;
  * the hook is a no-op for a session nobody recorded, never fails, and saves a recorded one
    detached;
  * finalize's MANIFEST line; the session zip and the delivery leave the conversation out;
  * the decisions log, AGENTS.md's reviewer checklist, the registry pointer.
"""
import base64
import contextlib
import io
import json
import os
import subprocess
import sys
import tempfile
import time
import unittest
import zipfile
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import save_transcript as st  # noqa: E402
from job_env import job_env  # noqa: E402

PY = sys.executable
SID = "11111111-2222-4333-8444-555555555555"
SID2 = "66666666-7777-4888-9999-aaaaaaaaaaaa"

# One fake of each kind of secret; none may survive into any saved file.
SECRETS = {
    "google": "AIzaSyFAKEfake0123456789abcdefghijklmn",
    "slack_webhook": "https://hooks.slack.com/services/TFAKE0000/BFAKE0000/fakefakefakefake1234",
    "password": "FakePassw0rdXYZ987",
    "env_assign": "faketok_1234567890abcdef",
    "json_key": "fakeapikey1234567890",
    "anthropic": "sk-ant-api03-FAKEfakeFAKEfakeFAKEfake0123456789",
    "bearer": "fakebearer1234567890abcdefXYZ",
    "dict_key": "plantedDictTokenValue12345",
    "env_value": "plantedEnvSecret123456789",             # MY_SERVICE_TOKEN's value, printed bare
    "file_value": "plantedFileSecretABCDEFG123",          # ~/.coreomics_token, printed bare
}


# The skill's text injected when it is loaded: where an analysis's part of a conversation starts.
LOAD = {"type": "user", "isMeta": True, "uuid": "load", "timestamp": "2026-09-28T09:59:00.000Z",
        "message": {"role": "user", "content": "Base directory for this skill: /x/skills/"
                                               "ucdavis-proteomics-core-pipeline\n\n# Proteomics"
                                               " Pipeline ..."}}


def entries(sid=SID, day="2026-09-28", hour=10):
    t = lambda m: f"{day}T{hour:02d}:{m:02d}:00.000Z"      # noqa: E731
    base = {"sessionId": sid, "version": "2.1.284", "cwd": "/x", "userType": "external"}
    s = SECRETS
    return [
        dict(LOAD, timestamp=t(0)),
        {"type": "permission-mode", "permissionMode": "default", "sessionId": sid},
        dict(base, type="user", timestamp=t(0), uuid="u1",
             message={"role": "user", "content": "Analyse PROT_0756 please. My Gemini key is "
                                                 f"{s['google']} and the webhook {s['slack_webhook']}"}),
        dict(base, type="user", timestamp=t(1), isMeta=True,
             message={"role": "user", "content": "<the skill's SKILL.md text>"}),
        dict(base, type="assistant", timestamp=t(2),
             message={"role": "assistant", "content": [
                 {"type": "thinking", "thinking": "plan", "signature": "sig"},
                 {"type": "text", "text": "I will run the search."},
                 {"type": "tool_use", "id": "toolu_1", "name": "Bash", "input": {
                     "command": f"export COREOMICS_TOKEN={s['env_assign']}\n"
                                "python3 scripts/run_search.py --out o",
                     "description": "Run the search"}}]}),
        dict(base, type="user", timestamp=t(3),
             message={"role": "user", "content": [{"type": "tool_result", "tool_use_id": "toolu_1",
                                                   "content": "\n".join(f"search line {i}"
                                                                        for i in range(100))}]}),
        dict(base, type="assistant", timestamp=t(4),
             message={"role": "assistant", "content": [
                 {"type": "tool_use", "id": "toolu_2", "name": "Read",
                  "input": {"file_path": "/x/config.json", "token": s["dict_key"]}}]}),
        dict(base, type="user", timestamp=t(5),
             message={"role": "user", "content": [{"type": "tool_result", "tool_use_id": "toolu_2",
                                                   "content": [{"type": "text", "text":
                                                       f'{{"api_key": "{s["json_key"]}"}}\n'
                                                       f"password={s['password']}\n"
                                                       f"Authorization: Bearer {s['bearer']}\n"
                                                       f"{s['anthropic']}\n{s['env_value']}\n"
                                                       f"{s['file_value']}"}]}]}),
        dict(base, type="system", subtype="compact_boundary", content="Conversation compacted",
             timestamp=t(6)),
        dict(base, type="user", timestamp=t(7), isCompactSummary=True,
             message={"role": "user", "content": "Summary of the earlier conversation."}),
        dict(base, type="attachment", timestamp=t(8), attachment={"type": "hook"}),
        {"type": "brand-new-kind", "foo": 1, "sessionId": sid},
        dict(base, type="assistant", timestamp=t(9), message="not a dict"),
        "a bare string entry",
        dict(base, type="assistant", timestamp=t(10),
             message={"role": "assistant", "content": [{"type": "server_tool_use", "id": "s1"},
                                                       {"type": "text", "text": "Done: 6,112 proteins."}]}),
    ]


def write_transcript(root, sid=SID, project="-Users-someone-project", **kw):
    path = os.path.join(root, "projects", project, sid + ".jsonl")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for e in entries(sid, **kw):
            fh.write(json.dumps(e) + "\n")
        fh.write("{ not json\n")
    return path


class Env(unittest.TestCase):
    """HOME, the Claude Code config dir and the skill's config dir all inside a temp folder."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.d = self._td.name
        self.home = os.path.join(self.d, "home")
        self.claude = os.path.join(self.d, "claude_cfg")
        self.cfg = os.path.join(self.d, "skill_cfg")
        os.makedirs(self.home)
        with open(os.path.join(self.home, ".coreomics_token"), "w") as fh:
            fh.write(SECRETS["file_value"] + "\n")
        self.session = os.path.join(self.d, "PROT_0756", "session")    # PROT_0756 is its own
        os.makedirs(os.path.join(self.session, "logs"))
        self.env = {"HOME": self.home, "CLAUDE_CONFIG_DIR": self.claude,
                    "SKILL_CONFIG_DIR": self.cfg, "CLAUDE_CODE_SESSION_ID": SID,
                    "MY_SERVICE_TOKEN": SECRETS["env_value"]}
        self._patch = mock.patch.dict(os.environ, self.env)
        self._patch.start()

    def tearDown(self):
        self._patch.stop()
        self._td.cleanup()

    def conv(self, *parts):
        return os.path.join(self.session, "logs", "conversation", *parts)

    def read(self, *parts):
        with open(self.conv(*parts), encoding="utf-8") as fh:
            return fh.read()

    def cli(self, *argv):
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = st.main(list(argv))
        return rc, json.loads(out.getvalue())


class Finding(Env):
    def test_found_by_its_id_wherever_the_project_folder_is(self):
        want = write_transcript(self.claude)
        write_transcript(self.claude, SID2, project="-other")                  # another one
        self.assertEqual(st.find_transcript(SID), want)
        older = write_transcript(os.path.join(self.home, ".claude"), SID)      # ~/.claude too
        os.utime(older, (1, 1))
        self.assertEqual(st.find_transcript(SID), want)                         # the newest
        os.remove(want)
        self.assertEqual(st.find_transcript(SID), older)
        for bad in ("../../etc/passwd", "*", "", None, "a"):
            self.assertIsNone(st.find_transcript(bad), bad)

    def test_no_transcript_or_no_claude_code_is_a_clean_skip(self):
        rc, out = self.cli(self.session)
        self.assertEqual(rc, 0)
        self.assertIn(f"no transcript for Claude Code session {SID}", out["skipped"])
        self.assertFalse(os.path.exists(self.conv()))                          # nothing written
        with mock.patch.dict(os.environ, {"CLAUDE_CODE_SESSION_ID": ""}):
            rc, out = self.cli(self.session)
        self.assertEqual(rc, 0)
        self.assertIn("record decisions with log_decision.py", out["skipped"])
        rc, out = self.cli(os.path.join(self.d, "no-such-session"))
        self.assertEqual(rc, 0)
        self.assertIn("no session folder", out["skipped"])


class Saving(Env):
    def test_every_kind_of_secret_is_redacted_from_every_file(self):
        write_transcript(self.claude)
        rc, out = self.cli(self.session)
        self.assertEqual(rc, 0, out)
        self.assertEqual(out["session_id"], SID)
        self.assertEqual(out["redacted_strings"], 4)   # the prompt, a command, a field, a result
        files = os.listdir(self.conv())
        self.assertEqual(sorted(files), sorted([SID + ".jsonl", "index.json", "conversation.md"]))
        for f in files:
            text = self.read(f)
            for kind, secret in SECRETS.items():
                self.assertNotIn(secret, text, f"{kind} in {f}")
        lines = self.read(SID + ".jsonl").splitlines()
        self.assertEqual(len(lines), len(entries()) + 1)
        for ln in lines[:-1]:
            json.loads(ln)                                   # the structure survives redaction
        raw = self.read(SID + ".jsonl")
        self.assertIn("[redacted]", raw)
        self.assertIn("python3 scripts/run_search.py --out o", raw)   # what is not secret stays
        leftovers = [f for f in os.listdir(os.path.dirname(self.conv())) if "part" in f]
        self.assertEqual(leftovers, [])
        self.assertFalse([f for f in os.listdir(self.conv()) if f.endswith(".part")])

    def test_the_readable_file_never_crashes_on_an_unknown_shape(self):
        write_transcript(self.claude)
        self.cli(self.session)
        md = self.read("conversation.md")
        self.assertIn("### 2026-09-28 10:00:00 UTC -- User", md)
        self.assertIn("Analyse PROT_0756 please.", md)
        self.assertIn("### 2026-09-28 10:02:00 UTC -- Assistant\n\nI will run the search.", md)
        self.assertIn("**Tool: Bash** -- Run the search", md)
        self.assertIn("python3 scripts/run_search.py --out o", md)
        self.assertIn("search line 0", md)                     # the result, cut to ~40 lines
        self.assertIn("... 60 lines not shown ...", md)
        self.assertIn("search line 99", md)
        self.assertNotIn("search line 50", md)
        self.assertIn("**Tool: Read**", md)
        self.assertIn("the conversation was compacted here", md)
        self.assertIn("compaction summary", md)
        self.assertIn("[unrecognised entry: type='brand-new-kind']", md)
        self.assertIn("[unrecognised entry: type='assistant']", md)   # message not a dict
        self.assertIn("[unrecognised entry: not JSON]", md)
        self.assertIn("[unrecognised entry: block type='server_tool_use']", md)
        self.assertIn("Done: 6,112 proteins.", md)             # after the odd block, still read
        self.assertIn("*Not shown: 2 metadata entries, 2 meta messages (skill text, reminders), "
                      "1 thinking blocks.*", md)
        self.assertNotIn("SKILL.md text", md)

    def test_two_conversations_two_files_one_index_and_resaving_changes_nothing(self):
        write_transcript(self.claude, SID, hour=14)
        write_transcript(self.claude, SID2, project="-p2", hour=9)           # started earlier
        self.cli(self.session)
        rc, out = self.cli(self.session, "--session-id", SID2)
        self.assertEqual(out["conversations"], 2)
        idx = json.loads(self.read("index.json"))
        self.assertEqual([c["id"] for c in idx["conversations"]], [SID2, SID])  # in order
        c = idx["conversations"][1]
        self.assertEqual((c["first"], c["last"], c["claude_code_version"], c["entries"]),
                         ("2026-09-28T14:00:00.000Z", "2026-09-28T14:10:00.000Z", "2.1.284",
                          len(entries()) + 1))
        self.assertEqual(c["bytes"], os.path.getsize(self.conv(SID + ".jsonl")))
        md = self.read("conversation.md")
        self.assertLess(md.index(f"## Conversation 1: `{SID2}`"),
                        md.index(f"## Conversation 2: `{SID}`"))
        before = {f: self.read(f) for f in (SID + ".jsonl", SID2 + ".jsonl")}
        sha = c["sha256"]
        self.cli(self.session)                                   # the same conversation again
        idx = json.loads(self.read("index.json"))
        self.assertEqual(len(idx["conversations"]), 2)
        self.assertEqual(idx["conversations"][1]["sha256"], sha)
        self.assertEqual({f: self.read(f) for f in before}, before)

    def test_remember_and_lookup(self):
        self.assertEqual(st.lookup(SID), [])
        rc, out = self.cli(self.session, "--remember")
        self.assertTrue(out["remembered"])
        self.assertEqual(st.lookup(SID)[0]["session_dir"], os.path.abspath(self.session))
        self.assertEqual(st._record(st._read_map(), SID)["segments"][0]["start"],
                         {"line": 0, "after_uuid": None})                # at the skill's load
        with mock.patch.dict(os.environ, {"CLAUDE_CODE_SESSION_ID": ""}):
            self.assertEqual(st.remember(self.session), (None, None))
        with open(st.map_path(), "w") as fh:                        # a damaged map: no crash
            fh.write("{")
        self.assertEqual(st.lookup(SID), [])


def write_lines(path, objs, tail=None):
    """A transcript of these entries; `tail`: a last line with no newline (still being written)."""
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for o in objs:
            fh.write(json.dumps(o) + "\n")
        if tail is not None:
            fh.write(tail)
    return path


def msg(i, text, sid=SID):
    return {"type": "user", "uuid": f"u{i}", "sessionId": sid, "version": "2.1.284",
            "timestamp": f"2026-09-28T10:{i:02d}:00.000Z",
            "message": {"role": "user", "content": text}}


class ReviewFixes(Env):
    """podcast-reviewer on 32a0780: redaction gaps, a partial last line, one conversation that
    runs two analyses, permissions, overlapping saves, the HIVE fetch, the delivered AGENTS.md."""

    def path(self, sid=SID):
        return os.path.join(self.claude, "projects", "-p", sid + ".jsonl")

    def saved(self, sid=SID):
        with open(self.conv(sid + ".jsonl"), encoding="utf-8") as fh:
            return fh.read()

    def test_the_redaction_gaps(self):
        v = SECRETS["env_value"]                            # held here (MY_SERVICE_TOKEN)
        b64 = base64.b64encode(b"xy" + v.encode()).decode()
        url = base64.urlsafe_b64encode(b"z" + v.encode()).decode()
        gho = "gho_" + "A1b2C3d4E5f6G7h8I9j0K1l2"
        camel = "camelCaseSecretValue999"
        keyed = "AIza" + "SyKEYASDICTKEY0123456789abcdefgh"
        png = "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAADUlEQVR42mNk"
        objs = [LOAD, msg(0, f"split copy: {v[:9]}\n  {v[9:]} ; base64: {b64} ; urlsafe: {url}"),
                msg(1, f"gh: {gho} ; js: accessToken: '{camel}' ; refreshToken={camel}"),
                {"type": "assistant", "uuid": "a2", "timestamp": "2026-09-28T10:02:00.000Z",
                 "message": {"role": "assistant", "content": [
                     {"type": "tool_use", "id": "t", "name": "X",
                      "input": {keyed: 1, "clientSecret": camel}}]}},
                {"type": "user", "uuid": "u3", "timestamp": "2026-09-28T10:03:00.000Z",
                 "message": {"role": "user", "content": [
                     {"type": "tool_result", "tool_use_id": "t", "content": [
                         {"type": "image", "source": {"type": "base64", "media_type":
                                                      "image/png", "data": png}}]}]},
                 "toolUseResult": {"type": "image", "file": {"base64": png}}}]
        write_lines(self.path(), objs)
        rc, out = self.cli(self.session)
        for f in os.listdir(self.conv()):
            text = self.read(f)
            for bad in (v, v[9:], b64[4:-4], url[4:-4], gho, camel, keyed, png):
                self.assertNotIn(bad, text, f"{bad[:12]}... in {f}")
        self.assertIn('"text": "[image]"', self.saved())
        self.assertIn("[image]", self.read("conversation.md"))

    def test_every_place_core_submission_reads_the_coreomics_key_from(self):
        """A key in a UTF-16 file (Windows PowerShell's `>`) is the key api_token() sends, so the
        redaction knows it too; so does a key in a file check calls misplaced (Notepad's .txt)."""
        k16, ktxt = "utf16CoreOmicsKey0123456789", "notepadCoreOmicsKey0123456789"
        with open(os.path.join(self.home, ".coreomics_token"), "wb") as fh:
            fh.write((k16 + "\r\n").encode("utf-16"))
        with open(os.path.join(self.home, ".coreomics_token.txt"), "w") as fh:
            fh.write(ktxt + "\n")
        vals = st.secret_values()
        self.assertIn(k16, vals)
        self.assertIn(ktxt, vals)

    def test_a_failure_finding_the_coreomics_key_never_stops_the_redaction(self):
        import core_submission as cs
        with mock.patch.object(cs, "key_files", side_effect=RuntimeError("boom")):
            vals = st.secret_values()
        self.assertIn(SECRETS["env_value"], vals)
        self.assertIn(SECRETS["file_value"], vals)          # ~/.coreomics_token, the fallback

    def test_a_partial_last_line_is_left_for_the_next_save(self):
        write_lines(self.path(), [LOAD, msg(0, "complete"), msg(1, "also complete")],
                    tail='{"type": "user", "message": {"content": "HALF-WRITT')
        rc, out = self.cli(self.session)
        self.assertEqual(out["entries"], 3)
        self.assertNotIn("HALF-WRITT", self.saved())
        self.assertEqual(len(self.saved().splitlines()), 3)
        for ln in self.saved().splitlines():
            json.loads(ln)

    def test_permissions_are_the_owner_and_group(self):
        write_transcript(self.claude)
        self.cli(self.session)
        import log_decision
        log_decision.append(self.session, "x", "y")
        mode = lambda p: oct(os.stat(p).st_mode & 0o777)           # noqa: E731
        self.assertEqual(mode(self.conv()), "0o750")
        for f in os.listdir(self.conv()):
            self.assertEqual(mode(self.conv(f)), "0o640", f)
        self.assertEqual(mode(os.path.join(self.session, "logs", "decisions.md")), "0o640")

    def test_the_index_is_rebuilt_from_the_files(self):
        write_lines(self.path(), [LOAD] + [msg(i, f"one {i}") for i in range(4)])
        write_lines(self.path(SID2), [LOAD] + [msg(i, "two", SID2) for i in range(2)])
        self.cli(self.session)
        os.remove(self.conv("index.json"))                          # lost, or a stale list
        self.cli(self.session, "--session-id", SID2)
        idx = json.loads(self.read("index.json"))
        self.assertEqual(sorted(c["id"] for c in idx["conversations"]), sorted([SID, SID2]))

    def test_an_older_snapshot_finishing_later_is_refused(self):
        # the race guard: a save of the transcript as it was, finishing after a newer save
        full = [LOAD] + [msg(i, f"one {i}") for i in range(4)]
        older = os.path.join(self.d, "older.jsonl")
        write_lines(older, full[:3])                                # a snapshot, 3 lines in
        write_lines(self.path(), full)
        self.cli(self.session)                                      # the newer save: 5 lines
        idx = json.loads(self.read("index.json"))["conversations"][0]
        self.assertEqual((idx["reach"], idx["reach_uuid"]), (5, "u3"))
        rc, out = self.cli(self.session, "--transcript", older)
        self.assertIn("reaches further into the conversation", out["kept"])
        self.assertEqual(len(self.saved().splitlines()), 5)
        self.assertIn("one 3", self.saved())

    def test_a_newer_save_finishing_during_the_redaction_wins(self):
        # the same race, closer: the newer save lands while the older one is redacting
        full = [LOAD] + [msg(i, f"one {i}") for i in range(4)]
        older, newer = (os.path.join(self.d, n) for n in ("older.jsonl", "newer.jsonl"))
        write_lines(older, full[:3])
        write_lines(newer, full)
        conv = os.path.join(self.d, "conv")
        real = st._select

        def select(lines, ranges):
            st._select = real
            self.assertTrue(st.save_into(conv, newer, SID)["written"])
            return real(lines, ranges)
        with mock.patch.object(st, "_select", select):
            got = st.save_into(conv, older, SID)
        self.assertFalse(got["written"])
        self.assertIn("reaches further into the conversation", got["why"])
        with open(os.path.join(conv, SID + ".jsonl")) as fh:
            self.assertIn("one 3", fh.read())
        self.assertEqual((got["entry"]["reach"], got["entry"]["reach_uuid"]), (5, "u3"))

    def test_the_delivered_agents_md_leaves_out_the_review_records(self):
        import session_docs
        import test_deposit_package as tdp
        import log_decision
        p = tdp.dia_session(os.path.join(self.d, "s"))
        log_decision.append(p["session_dir"], "x", "y")
        os.makedirs(os.path.join(p["session_dir"], "logs", "conversation"))
        with open(os.path.join(p["session_dir"], "logs", "conversation", "conversation.md"),
                  "w") as fh:
            fh.write("# c\n")
        f = session_docs.gather(p["session_dir"])
        self.assertIn("## Reviewing this analysis", session_docs.agents_md(f))
        delivered = session_docs.agents_md(f, for_delivery=True)
        for word in ("Reviewing this analysis", "logs/conversation", "conversation.md",
                     "decisions.md", "transcript"):
            self.assertNotIn(word, delivered)


def say(i, who, text, sid=SID):
    """One turn: the user's words, or the assistant's (with a tool call for a command)."""
    if who == "user":
        return msg(i, text, sid)
    return {"type": "assistant", "uuid": f"a{i}", "sessionId": sid,
            "timestamp": f"2026-09-28T11:{i:02d}:00.000Z",
            "message": {"role": "assistant", "content": [{"type": "text", "text": text}]}}


class Segments(Env):
    """One conversation, several analyses: each analysis's copy holds its own part only
    (review of bc7845c: client A's lines reached client B's copy; a status peek moved lines)."""

    def setUp(self):
        super().setUp()
        self.lines, self.t = [], os.path.join(self.claude, "projects", "-p", SID + ".jsonl")

    def add(self, *turns):
        for who, text in turns:
            self.lines.append(say(len(self.lines), who, text))
        write_lines(self.t, self.lines)

    def sess(self, name):
        d = os.path.join(self.d, name)
        os.makedirs(d, exist_ok=True)
        return d

    def copy(self, d):
        self.cli(d)
        path = os.path.join(d, "logs", "conversation", SID + ".jsonl")
        with open(path, encoding="utf-8") as fh:
            return fh.read()

    def status(self, d, *extra):
        r = subprocess.run([PY, os.path.join(SCRIPTS, "checkpoint.py"), "status", "--session", d,
                            *extra], capture_output=True, text=True, env=dict(os.environ),
                           timeout=60)
        return r                                   # no checkpoint record: exit 1 after remember

    def test_a_client_named_before_init_stays_out_of_the_next_clients_copy(self):
        # the reviewer's (a): client A (PROT_0801) discussed, then `init` for client B
        self.lines.append(LOAD)
        self.add(("user", "For PROT_0801, drop sample 7 -- the PI asked."),
                 ("assistant", "Noted: PROT_0801 without sample 7."),
                 ("user", "Now the other one, PROT_0807."))
        b = self.sess(os.path.join("PROT_0807", "session"))
        entry, warn = st.remember(b)                              # session.py init for B
        self.assertIn("names PROT_0801, not this session's submission", warn)
        self.add(("assistant", "PROT_0807: 12 runs, two groups."))
        got = self.copy(b)
        for other in ("PROT_0801", "drop sample 7", "the PI asked"):
            self.assertNotIn(other, got)
        self.assertIn("PROT_0807: 12 runs, two groups.", got)

    def test_the_skill_load_starts_the_first_copy_and_is_found(self):
        self.add(("user", "What's the weather like?"), ("assistant", "No idea."),
                 ("user", "Please analyse my DIA runs."))           # the request that loads it
        self.lines.append({"type": "assistant", "uuid": "sk", "timestamp": "2026-09-28T11:05:00Z",
                           "message": {"role": "assistant", "content": [
                               {"type": "tool_use", "id": "s", "name": "Skill",
                                "input": {"skill": "ucdavis-proteomics-core-pipeline"}}]}})
        self.add(("user", "Groups: Old vs Young, 3 each."), ("assistant", "Confirmed."))
        x = self.sess("X")
        self.assertEqual(st.remember(x)[0]["key"], "local:" + os.path.abspath(x))
        got = self.copy(x)
        self.assertNotIn("weather", got)
        self.assertIn("Please analyse my DIA runs.", got)                 # the opening request
        self.assertIn("Groups: Old vs Young, 3 each.", got)              # the design, pre-init
        y = self.sess("Y")                                            # no skill load at all
        other = os.path.join(self.claude, "projects", "-q", SID2 + ".jsonl")
        write_lines(other, [msg(0, "hello", SID2)])
        entry, warn = st.remember(y, session_id=SID2)
        self.assertIn("the skill's first load is not in this conversation", warn)

    def skill_call(self):
        self.lines.append({"type": "assistant", "uuid": f"sk{len(self.lines)}",
                           "timestamp": "2026-09-28T11:30:00Z",
                           "message": {"role": "assistant", "content": [
                               {"type": "tool_use", "id": "s", "name": "Skill",
                                "input": {"skill": "ucdavis-proteomics-core-pipeline"}}]}})
        write_lines(self.t, self.lines)

    def test_the_opening_request_of_an_auto_invoked_skill_is_kept(self):
        self.add(("user", "Analyse PROT_0802, groups Old vs Young."))
        self.skill_call()
        self.add(("assistant", "PROT_0802: 12 runs."))
        b = self.sess(os.path.join("PROT_0802", "session"))
        entry, warn = st.remember(b)
        self.assertIsNone(warn)
        self.assertIn("Analyse PROT_0802, groups Old vs Young.", self.copy(b))

    def test_an_opening_request_naming_another_client_is_left_out(self):
        self.add(("user", "Compare with PROT_0801 later; now analyse PROT_0802."))
        self.skill_call()
        self.add(("assistant", "PROT_0802: 12 runs."))
        b = self.sess(os.path.join("PROT_0802", "session"))
        entry, warn = st.remember(b)
        self.assertIsNone(warn)                                    # the load: after the request
        got = self.copy(b)
        self.assertNotIn("PROT_0801", got)
        self.assertIn("PROT_0802: 12 runs.", got)

    def test_a_compaction_after_a_switch_carries_no_other_analysis(self):
        # re-verification of 85eddbe: a compaction summary sums up the whole conversation, so
        # one made in Y after X's segment put "CLIENT-X PROT_0801 (Dr. Smith, dropped sample
        # 7)" in Y's copy
        self.lines.append(LOAD)
        x, y = self.sess(os.path.join("PROT_0801", "X")), self.sess(os.path.join("PROT_0802", "Y"))
        self.add(("user", "CLIENT-X PROT_0801: Dr. Smith wants sample 7 dropped."))
        st.remember(x)
        self.add(("assistant", "Dropped sample 7 for Dr. Smith."))
        st.remember(y)
        self.add(("assistant", "Y: search submitted."))
        summary = "Summary: worked on CLIENT-X PROT_0801 (Dr. Smith, dropped sample 7), then Y."
        self.lines += [
            {"type": "system", "subtype": "compact_boundary", "content": "Conversation compacted",
             "timestamp": "2026-09-28T12:00:00Z", "uuid": "cb"},
            {"type": "user", "isCompactSummary": True, "isVisibleInTranscriptOnly": True,
             "uuid": "cs", "timestamp": "2026-09-28T12:00:01Z",
             "message": {"role": "user", "content": summary}},
            {"type": "summary", "summary": summary, "leafUuid": "cs"},
            {"type": "ai-title", "aiTitle": "CLIENT-X Dr. Smith sample 7 and Y", "sessionId": SID}]
        self.add(("assistant", "Y: DE done."))
        cy, md = self.copy(y), self.read_md(y)
        for other in ("CLIENT-X", "Dr. Smith", "sample 7", "PROT_0801"):
            self.assertNotIn(other, cy)
            self.assertNotIn(other, md)
        self.assertEqual(cy.count(st.SUMMARY_OMITTED), 1)          # the compaction summary
        self.assertIn("[2 other entries omitted]", cy)             # "summary", "ai-title"
        self.assertIn(f"*{st.SUMMARY_OMITTED}*", md)
        self.assertIn("Y: DE done.", cy)
        cx = self.copy(x)                                          # X's copy holds no Y either
        self.assertNotIn("Y: search submitted", cx)

    def compaction(self, text):
        self.lines += [
            {"type": "system", "subtype": "compact_boundary", "content": "Conversation compacted",
             "timestamp": "2026-09-28T12:00:00Z", "uuid": f"cb{len(self.lines)}"},
            {"type": "user", "isCompactSummary": True, "isVisibleInTranscriptOnly": True,
             "uuid": f"cs{len(self.lines)}", "timestamp": "2026-09-28T12:00:01Z",
             "message": {"role": "user", "content": text}}]
        write_lines(self.t, self.lines)

    def test_from_now_then_a_compaction_carries_no_earlier_client(self):
        # the reviewer's A: an unmapped client before --from-now survived inside a later summary
        self.lines.append(LOAD)
        self.add(("user", "About CLIENT-A PROT_0801: Dr. Jones wants the outlier kept."),
                 ("assistant", "Noted for CLIENT-A."))
        x = self.sess("X")
        self.cli(x, "--remember", "--from-now")
        self.add(("assistant", "X: 12 runs."))
        self.compaction("Summary: CLIENT-A PROT_0801 (Dr. Jones, keep the outlier); then X.")
        self.add(("assistant", "X: DE done."))
        got, md = self.copy(x), self.read_md(x)
        for other in ("CLIENT-A", "PROT_0801", "Dr. Jones", "outlier"):
            self.assertNotIn(other, got)
            self.assertNotIn(other, md)
        self.assertIn(st.SUMMARY_OMITTED, got)
        self.assertIn("Conversation compacted", got)               # the boundary itself stays
        self.assertIn("X: DE done.", got)

    def test_a_single_analysis_from_the_skill_load_drops_a_summary_of_what_came_before(self):
        self.add(("user", "Unrelated: how did CLIENT-A's tau blots look?"),
                 ("assistant", "CLIENT-A's blots were fine."))
        self.lines.append(dict(LOAD, uuid="load2"))                   # the skill loads here
        self.add(("user", "Groups: Old vs Young."))
        x = self.sess("X")
        st.remember(x)
        self.compaction("Summary: CLIENT-A tau blots; then X, Old vs Young.")
        self.add(("assistant", "X: done."))
        got = self.copy(x)
        self.assertNotIn("CLIENT-A", got)
        self.assertIn(st.SUMMARY_OMITTED, got)
        self.assertIn("Groups: Old vs Young.", got)

    def test_a_tool_result_before_the_skill_is_not_the_opening_request(self):
        # the reviewer's B: the request -> a Bash call naming PROT_0801 -> the Skill call
        self.add(("user", "Check Dr. Smith's tau results, then set up the new analysis."))
        self.lines += [
            {"type": "assistant", "uuid": "b1", "timestamp": "2026-09-28T11:20:00Z",
             "message": {"role": "assistant", "content": [
                 {"type": "tool_use", "id": "bash1", "name": "Bash",
                  "input": {"command": "cat ~/core/PROT_0801/tau_summary.txt"}}]}},
            {"type": "user", "uuid": "b2", "timestamp": "2026-09-28T11:20:05Z",
             "message": {"role": "user", "content": [
                 {"type": "tool_result", "tool_use_id": "bash1",
                  "content": "PROT_0801 tau: Dr. Smith's knockout lowers pTau 40%."}]}}]
        self.skill_call()
        self.add(("user", "Groups: WT vs KO."))
        b = self.sess(os.path.join("PROT_0802", "session"))
        st.remember(b)
        got = self.copy(b)
        for other in ("Dr. Smith", "PROT_0801", "tau", "pTau"):
            self.assertNotIn(other, got)
        self.assertIn("Groups: WT vs KO.", got)

    def other_entries(self, who):
        """Each whole-conversation entry type, naming another client."""
        return [
            {"type": "custom-title", "customTitle": f"{who} and more", "sessionId": SID},
            {"type": "last-prompt", "lastPrompt": f"tell me about {who}", "sessionId": SID},
            {"type": "queue-operation", "operation": "enqueue", "content": f"then {who}",
             "timestamp": "2026-09-28T11:40:00Z", "sessionId": SID},
            {"type": "attachment", "attachment": {"type": "file", "content": f"{who} notes"},
             "timestamp": "2026-09-28T11:40:01Z", "uuid": "att"},
            {"type": "system", "subtype": "away_summary", "content": f"While away: {who} ...",
             "timestamp": "2026-09-28T11:40:02Z", "uuid": "away"},
            {"type": "user", "isMeta": True, "uuid": "meta2", "timestamp": "2026-09-28T11:40:03Z",
             "message": {"role": "user", "content": f"<system-reminder>{who}</system-reminder>"}}]

    def test_a_partial_copy_holds_no_other_entry_type(self):
        self.lines.append(LOAD)
        x, y = self.sess("X"), self.sess("Y")
        self.add(("assistant", "X first."))
        st.remember(x)
        st.remember(y)
        self.add(("assistant", "Y first."))
        self.lines += self.other_entries("CLIENT-Q Dr. Lee")
        self.add(("assistant", "Y second."))
        got = self.copy(y)
        self.assertNotIn("CLIENT-Q", got)
        self.assertNotIn("Dr. Lee", got)
        self.assertIn("[6 other entries omitted]", got)            # one line for the run
        self.assertEqual(got.count("other entries omitted"), 1)
        self.assertIn("Y second.", got)
        self.assertIn("[6 other entries omitted]", self.read_md(y))

    def test_a_whole_conversation_copy_keeps_everything(self):
        self.lines.append(LOAD)
        x = self.sess("X")
        self.add(("assistant", "X first."))
        st.remember(x)
        self.lines += self.other_entries("CLIENT-Q Dr. Lee")
        self.compaction("Summary: X so far.")
        self.add(("assistant", "X second."))
        got = self.copy(x)
        for kept in ("customTitle", "lastPrompt", "queue-operation", "away_summary",
                     "attachment", "Summary: X so far.", "system-reminder"):
            self.assertIn(kept, got)
        self.assertNotIn("omitted", got)

    def test_metadata_before_the_first_message_does_not_make_a_copy_partial(self):
        # Claude Code writes permission-mode, mode ... before the first message: nothing to leak
        self.lines += [{"type": "permission-mode", "permissionMode": "default"},
                       {"type": "mode", "mode": "normal"}]
        self.add(("user", "Analyse PROT_0802, groups Old vs Young."))
        self.skill_call()
        b = self.sess(os.path.join("PROT_0802", "session"))
        st.remember(b)
        self.compaction("Summary: PROT_0802, Old vs Young.")
        self.lines.append({"type": "custom-title", "customTitle": "PROT_0802 run"})
        self.add(("assistant", "done."))
        got = self.copy(b)
        self.assertIn("Summary: PROT_0802, Old vs Young.", got)
        self.assertIn("customTitle", got)                          # whole: everything kept
        self.assertNotIn("omitted", got)

    def test_from_now_cleans_a_copy_already_saved(self):
        # the reviewer's 4: another client by name after the load -> init -> a save (the copy
        # holds Jones) -> --remember --from-now -> a save. The narrower copy has fewer entries
        # and was "kept"; now a changed selection always rewrites it.
        self.lines.append(LOAD)
        self.add(("user", "Dr. Jones's samples first -- the outlier stays in."),
                 ("assistant", "Noted for Dr. Jones."))
        b = self.sess("B")
        st.remember(b)                                             # init B, default start
        self.add(("assistant", "B: 12 runs."))
        self.assertIn("Dr. Jones", self.copy(b))                   # step 7's save
        rc, out = self.cli(b, "--remember", "--from-now")
        self.add(("assistant", "B: DE done."))
        got, md = self.copy(b), self.read_md(b)
        for text in (got, md):
            self.assertNotIn("Dr. Jones", text)
            self.assertNotIn("outlier", text)
        self.assertIn("B: DE done.", got)

    def test_a_switch_rewrites_a_copy_even_when_the_allowlist_shortens_it(self):
        self.lines.append(LOAD)
        x, y = self.sess("X"), self.sess("Y")
        self.add(("assistant", "X first."))
        st.remember(x)
        self.lines += self.other_entries("CLIENT-Q Dr. Lee")         # whole copy: kept
        self.add(("assistant", "X second."))
        self.assertIn("CLIENT-Q", self.copy(x))
        st.remember(y)                                              # X's segment closes here
        self.add(("assistant", "Y first."))
        got = self.copy(x)                                          # fewer entries, and newer
        self.assertNotIn("CLIENT-Q", got)
        self.assertIn("[6 other entries omitted]", got)
        self.assertIn("X second.", got)

    def test_a_save_that_started_before_the_segments_changed_writes_nothing(self):
        self.lines.append(LOAD)
        x = self.sess("X")
        self.add(("assistant", "X first."))
        st.remember(x)
        self.copy(x)
        conv = os.path.join(x, "logs", "conversation")
        with open(os.path.join(conv, SID + ".jsonl")) as fh:
            before = fh.read()
        got = st.save_into(conv, self.t, SID, parts=[({"line": 0, "after_uuid": None}, None)],
                           current=lambda: "a newer selection")
        self.assertFalse(got["written"])
        self.assertIn("changed while this save ran", got["why"])
        with open(os.path.join(conv, SID + ".jsonl")) as fh:
            self.assertEqual(fh.read(), before)

    def test_a_single_analysis_keeps_its_summaries(self):
        self.lines.append(LOAD)
        x = self.sess("X")
        self.add(("assistant", "X: conditions confirmed."))
        st.remember(x)
        self.lines.append({"type": "user", "isCompactSummary": True, "uuid": "cs",
                           "timestamp": "2026-09-28T12:00:01Z",
                           "message": {"role": "user", "content": "Summary: X, all of it."}})
        self.add(("assistant", "X: more."))
        got = self.copy(x)
        self.assertIn("Summary: X, all of it.", got)
        self.assertNotIn(st.SUMMARY_OMITTED, got)
        self.assertIn("compaction summary", self.read_md(x))

    def test_a_status_peek_moves_nothing(self):
        # the reviewer's (b): working on Y, a status on X (a peek), back to Y, more Y
        self.lines.append(LOAD)
        x, y = self.sess("X"), self.sess("Y")
        self.add(("user", "Start X."), ("assistant", "X: conditions confirmed."))
        st.remember(x)
        self.add(("assistant", "X: done for now."))
        st.remember(y)                                             # init Y: a switch
        self.add(("assistant", "Y: search submitted."))
        self.status(x)                                             # a peek at X
        self.status(y)                                             # Y again: nothing moves
        self.add(("user", "Y: the PI's private note -- do not share."),
                 ("assistant", "Y: noted."))
        cx, cy = self.copy(x), self.copy(y)
        for line in ("search submitted", "private note", "Y: noted"):
            self.assertNotIn(line, cx)
            self.assertIn(line, cy)
        for line in ("Start X.", "conditions confirmed", "X: done for now"):
            self.assertIn(line, cx)
            self.assertNotIn(line, cy)
        self.assertEqual([s["key"] for s in st.lookup(SID)],
                         ["local:" + os.path.abspath(x), "local:" + os.path.abspath(y)])

    def test_status_resume_switches_back(self):
        self.lines.append(LOAD)
        x, y = self.sess("X"), self.sess("Y")
        self.add(("assistant", "X one."))
        st.remember(x)
        self.add(("assistant", "X two."))
        st.remember(y)
        self.add(("assistant", "Y one."))
        self.status(x, "--resume")                                 # really back on X
        self.add(("assistant", "X three."))
        cx, cy = self.copy(x), self.copy(y)
        self.assertIn("X three.", cx)
        self.assertNotIn("X three.", cy)

    def test_x_y_x_y_interleaving(self):
        self.lines.append(LOAD)
        x, y = self.sess("X"), self.sess("Y")
        self.add(("assistant", "X part 1."))
        st.remember(x)
        self.add(("assistant", "X part 1b."))
        st.remember(y)
        self.add(("assistant", "Y part 1."))
        rc, out = self.cli(x, "--remember")                        # an explicit switch
        self.assertIn("moves here from", out["warning"])
        self.add(("assistant", "X part 2."))
        st.remember(y)
        self.add(("assistant", "Y part 2."))
        cx, cy = self.copy(x), self.copy(y)
        self.assertEqual([ln for ln in ("X part 1.", "X part 1b.", "X part 2.") if ln in cx],
                         ["X part 1.", "X part 1b.", "X part 2."])
        self.assertEqual(cx.count(st.OMITTED), 1)
        self.assertEqual(cy.count(st.OMITTED), 1)
        for ln in ("Y part 1.", "Y part 2."):
            self.assertNotIn(ln, cx)
            self.assertIn(ln, cy)
        for ln in ("X part 1.", "X part 2."):
            self.assertNotIn(ln, cy)
        self.assertLess(cx.index("X part 1b."), cx.index(st.OMITTED))
        self.assertLess(cx.index(st.OMITTED), cx.index("X part 2."))
        md = self.read_md(x)
        self.assertIn(f"*{st.OMITTED}*", md)

    def read_md(self, d):
        with open(os.path.join(d, "logs", "conversation", "conversation.md"), encoding="utf-8") as fh:
            return fh.read()

    def test_from_now(self):
        self.lines.append(LOAD)
        self.add(("user", "About client A, PROT_0801."), ("assistant", "A: fine."))
        x = self.sess("X")
        rc, out = self.cli(x, "--remember", "--from-now")          # an orchestrator's init
        self.assertIsNone(out.get("warning"))
        self.add(("assistant", "X from here."))
        got = self.copy(x)
        self.assertNotIn("client A", got)
        self.assertIn("X from here.", got)
        self.add(("user", "scratch that"), ("assistant", "X restarts."))
        rc, out = self.cli(x, "--remember", "--from-now")          # the active one, restarted
        self.assertIn("its copy now starts here", out["warning"])
        self.add(("assistant", "X after the restart."))
        got = self.copy(x)
        self.assertNotIn("X from here.", got)
        self.assertIn("X after the restart.", got)

    def test_session_init_takes_from_now(self):
        self.lines.append(LOAD)
        self.add(("user", "About client A, PROT_0801."))
        import test_deposit_package as tdp
        raw = os.path.join(self.d, "raw")
        os.makedirs(raw)
        tdp.make_d(os.path.join(raw, "r1.d"))
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "init", "--name", "b",
                            "--date", "2026-09-28", "--raw", os.path.join(raw, "*.d"),
                            "--base", self.d, "--transcript-from-now"], capture_output=True,
                           text=True, env=dict(os.environ), check=True)
        sd = json.loads(r.stdout)["paths"]["session_dir"]
        self.assertEqual(st.lookup(SID)[0]["key"], "local:" + sd)
        self.add(("assistant", "B from init."))
        got = self.copy(sd)
        self.assertNotIn("client A", got)
        self.assertIn("B from init.", got)


class Hook(Env):
    def hook(self, payload, env=None):
        raw = payload if isinstance(payload, str) else json.dumps(payload)
        t0 = time.monotonic()
        r = subprocess.run([PY, os.path.join(SCRIPTS, "transcript_hook.py")], input=raw,
                           capture_output=True, text=True, env=dict(os.environ, **(env or {})),
                           timeout=30)
        return r, time.monotonic() - t0

    def test_a_no_op_for_a_session_nobody_recorded(self):
        path = write_transcript(self.claude)
        r, _ = self.hook({"session_id": SID, "transcript_path": path,
                          "hook_event_name": "SessionEnd", "reason": "other"})
        self.assertEqual((r.returncode, r.stdout, r.stderr), (0, "", ""))
        self.assertFalse(os.path.exists(st.hook_log_path()))
        self.assertFalse(os.path.exists(self.conv()))

    def test_it_never_fails_the_session(self):
        for bad in ("{ not json", "[1, 2]", ""):
            r, _ = self.hook(bad)
            self.assertEqual(r.returncode, 0, bad)
        with open(st.hook_log_path()) as fh:
            self.assertIn("[error] unreadable hook input", fh.read())

    def test_a_recorded_session_is_saved_detached_and_fast(self):
        path = write_transcript(self.claude)
        st.remember(self.session, session_id=SID)
        r, took = self.hook({"session_id": SID, "transcript_path": path,
                             "hook_event_name": "PreCompact", "trigger": "auto"})
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertLess(took, 1.5)                    # SessionEnd hooks share a 1.5 s budget
        target = self.conv(SID + ".jsonl")
        for _ in range(100):                          # the save runs on after the hook exits
            if os.path.isfile(target) and os.path.isfile(self.conv("conversation.md")):
                break
            time.sleep(0.1)
        self.assertTrue(os.path.isfile(target))
        with open(st.hook_log_path()) as fh:
            log = fh.read()
        self.assertIn(f"PreCompact auto {SID}: saving to {os.path.abspath(self.session)}", log)
        for _ in range(50):
            with open(st.hook_log_path()) as fh:
                if "save_transcript {" in fh.read():
                    break
            time.sleep(0.1)
        with open(st.hook_log_path()) as fh:
            self.assertIn('save_transcript {"saved": ', fh.read())
        for secret in SECRETS.values():
            self.assertNotIn(secret, self.read(SID + ".jsonl"))


class Hive(Env):
    """A HIVE session driven from this computer: redacted here, put there -- with a stand-in for
    hive_exec.sh that maps the 'remote' onto a local folder."""

    FAKE = r'''#!/usr/bin/env bash
case "$1" in
  --get) [ -n "$FAIL_GET" ] && { echo "scp: refused" >&2; exit 2; }
         cp -R "$2" "$3" ;;
  --put) cp "$2" "$3" ;;
  *) [ -n "$FAIL_SSH" ] && { echo "ssh: connect: timed out" >&2; exit 255; }
     bash -c "$1" ;;
esac
'''

    def test_saved_on_hive_merging_what_is_there(self):
        fake = os.path.join(self.d, "hive_exec.sh")
        with open(fake, "w") as fh:
            fh.write(self.FAKE)
        remote = os.path.join(self.d, "PROT_0756", "hive session's dir")   # quoting, for real
        os.makedirs(remote)
        write_transcript(self.claude, SID2, project="-p2", hour=9)
        write_transcript(self.claude)
        with mock.patch.dict(os.environ, {"SKILL_HIVE_EXEC": fake}):
            rc, out = self.cli("--hive", remote, "--session-id", SID2)
            self.assertEqual(out["conversations"], 1, out)
            rc, out = self.cli("--hive", remote)
        self.assertTrue(out["hive"])
        self.assertEqual(out["conversations"], 2)                # fetched the first, added this
        conv = os.path.join(remote, "logs", "conversation")
        self.assertEqual(sorted(os.listdir(conv)), sorted([SID + ".jsonl", SID2 + ".jsonl",
                                                           "index.json", "conversation.md"]))
        self.assertEqual(st.lookup(SID)[0]["hive"], remote)
        for f in os.listdir(conv):
            with open(os.path.join(conv, f), encoding="utf-8") as fh:
                text = fh.read()
            for secret in SECRETS.values():
                self.assertNotIn(secret, text, f)


class HiveFailures(Env):
    def setUp(self):
        super().setUp()
        self.fake = os.path.join(self.d, "hive_exec.sh")
        with open(self.fake, "w") as fh:
            fh.write(Hive.FAKE)
        self.remote = os.path.join(self.d, "PROT_0756", "hive")
        self.rconv = os.path.join(self.remote, "logs", "conversation")
        write_transcript(self.claude, SID2, project="-p2", hour=9)
        write_transcript(self.claude)
        with mock.patch.dict(os.environ, {"SKILL_HIVE_EXEC": self.fake}):
            self.cli("--hive", self.remote, "--session-id", SID2)
        with open(os.path.join(self.rconv, "index.json")) as fh:
            self.before = fh.read()

    def test_a_fetch_that_fails_changes_nothing_on_hive(self):
        # review of 32a0780: scp's `--get <dir>/` always exits 2, and every Windows save rebuilt
        # index.json from one conversation
        with mock.patch.dict(os.environ, {"SKILL_HIVE_EXEC": self.fake, "FAIL_GET": "1"}):
            rc, out = self.cli("--hive", self.remote)
        self.assertEqual(rc, 0)
        self.assertIn("could not fetch", out["skipped"])
        self.assertIn("nothing was changed there", out["skipped"])
        self.assertEqual(sorted(os.listdir(self.rconv)),
                         sorted([SID2 + ".jsonl", "index.json", "conversation.md"]))
        with open(os.path.join(self.rconv, "index.json")) as fh:
            self.assertEqual(fh.read(), self.before)

    def test_an_unreachable_hive_changes_nothing(self):
        with mock.patch.dict(os.environ, {"SKILL_HIVE_EXEC": self.fake, "FAIL_SSH": "1"}):
            rc, out = self.cli("--hive", self.remote)
        self.assertIn("cannot reach HIVE", out["skipped"])
        with open(os.path.join(self.rconv, "index.json")) as fh:
            self.assertEqual(fh.read(), self.before)


class Decisions(unittest.TestCase):
    def test_appended_with_one_header_and_redacted(self):
        import log_decision
        with tempfile.TemporaryDirectory() as d:
            with mock.patch.dict(os.environ, {"CLAUDE_CODE_SESSION_ID": SID}):
                log_decision.append(d, "Groups: Old vs Young", "the user confirmed the design "
                                    "table", step="3. design", when="2026-09-28 10:00 UTC")
                log_decision.append(d, "Contrasts: Old - Young per bait",
                                    f"the user said so; token={SECRETS['env_assign']}",
                                    when="2026-09-28 11:00 UTC")
                # review of 32a0780: notify_slack.redact alone let these through
                with mock.patch.dict(os.environ, {"MY_SERVICE_TOKEN": SECRETS["env_value"]}):
                    log_decision._REDACTOR = None
                    log_decision.append(d, "Overrides", f"export COREOMICS_TOKEN="
                                        f"{SECRETS['env_assign']}; {{\"api_key\": "
                                        f"\"{SECRETS['json_key']}\"}}; {SECRETS['env_value']}")
                    log_decision._REDACTOR = None
            with open(os.path.join(d, "logs", "decisions.md"), encoding="utf-8") as fh:
                text = fh.read()
        self.assertEqual(text.count("# Decisions log"), 1)
        self.assertIn("## 2026-09-28 10:00 UTC -- Groups: Old vs Young\n- **Why:** the user "
                      "confirmed the design table\n- **Step:** 3. design\n"
                      f"- **Conversation:** `{SID}`", text)
        self.assertIn("## 2026-09-28 11:00 UTC -- Contrasts: Old - Young per bait", text)
        for s in ("env_assign", "json_key", "env_value"):
            self.assertNotIn(SECRETS[s], text, s)


class SessionIntegration(unittest.TestCase):
    """init records the conversation; finalize saves it with a MANIFEST line and keeps it out of
    the zip; deliver never ships it; AGENTS.md has the reviewer checklist; the registry points."""

    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.d = self._td.name
        self.claude = os.path.join(self.d, "claude_cfg")
        self.env = job_env(self.d, CLAUDE_CONFIG_DIR=self.claude, CLAUDE_CODE_SESSION_ID=SID,
                           HOME=os.path.join(self.d, "home"))

    def tearDown(self):
        self._td.cleanup()

    def session(self):
        import test_deposit_package as tdp
        os.makedirs(os.path.join(self.d, "home"), exist_ok=True)
        raw = os.path.join(self.d, "raw")
        os.makedirs(raw)
        for r in tdp.RUNS:
            tdp.make_d(os.path.join(raw, r + ".d"))
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "init", "--name", "demo",
                            "--date", "2026-09-28", "--raw", os.path.join(raw, "*.d"),
                            "--base", self.d], capture_output=True, text=True, env=self.env,
                           check=True)
        return json.loads(r.stdout)["paths"]

    def test_init_records_it_and_finalize_saves_it_but_never_zips_it(self):
        p = self.session()
        with mock.patch.dict(os.environ, {"SKILL_CONFIG_DIR": self.env["SKILL_CONFIG_DIR"]}):
            self.assertEqual(st.lookup(SID)[0]["session_dir"], p["session_dir"])
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "finalize", "--dir",
                            p["session_dir"], "--zip"], capture_output=True, text=True,
                           env=self.env)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(p["manifest_txt"], encoding="utf-8") as fh:
            man = fh.read()
        # Claude Code, and nothing saved: a failure
        self.assertIn("[SKIPPED] Analysis conversation (logs/conversation)", man)
        self.assertIn(f"no transcript for Claude Code session {SID}", man)
        env = dict(self.env)                                   # not Claude Code, and no log:
        env.pop("CLAUDE_CODE_SESSION_ID")                      # nothing records how it was done
        subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "finalize", "--dir",
                        p["session_dir"]], capture_output=True, text=True, env=env, check=True)
        with open(p["manifest_txt"], encoding="utf-8") as fh:
            self.assertIn("no conversation and no decisions log (log_decision.py)", fh.read())
        import log_decision
        with mock.patch.dict(os.environ, {"SKILL_CONFIG_DIR": self.env["SKILL_CONFIG_DIR"]}):
            log_decision.append(p["session_dir"], "Groups: ctrl vs trt", "the user confirmed")
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "finalize", "--dir",
                            p["session_dir"]], capture_output=True, text=True, env=self.env)
        with open(p["manifest_txt"], encoding="utf-8") as fh:  # Claude Code, nothing saved:
            self.assertIn("[SKIPPED] Analysis conversation", fh.read())   # a failure, even so
        env = dict(self.env)                                   # another agent + its log: INFO
        env.pop("CLAUDE_CODE_SESSION_ID")
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "finalize", "--dir",
                            p["session_dir"]], capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(p["manifest_txt"], encoding="utf-8") as fh:
            man = fh.read()
        self.assertIn("[INFO]    Analysis conversation (logs/conversation)", man)
        self.assertIn("the decisions log (logs/decisions.md) is the record", man)
        write_transcript(self.claude)
        r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "finalize", "--dir",
                            p["session_dir"], "--zip"], capture_output=True, text=True,
                           env=self.env)
        self.assertEqual(r.returncode, 0, r.stderr)
        res = json.loads(r.stdout)
        with open(p["manifest_txt"], encoding="utf-8") as fh:
            man = fh.read()
        self.assertRegex(man, r"\[OK\] +Analysis conversation \(logs/conversation\) +-- 1 "
                              r"conversation\(s\), redacted -- Core-internal: never delivered, "
                              r"not in the zip")
        conv = os.path.join(p["session_dir"], "logs", "conversation")
        self.assertTrue(os.path.isfile(os.path.join(conv, SID + ".jsonl")))
        names = zipfile.ZipFile(res["zip"]).namelist()
        self.assertFalse([n for n in names if "/logs/conversation" in n], names)
        self.assertEqual(res["zip_excluded"]["logs/conversation (the analysis conversation: "
                                             "Core-internal, kept on disk)"], 3)
        self.assertFalse([n for n in names if n.endswith("logs/decisions.md")], names)
        self.assertEqual(res["zip_excluded"]["logs/decisions.md (the decisions log: "
                                             "Core-internal, kept on disk)"], 1)
        with open(os.path.join(p["session_dir"], "AGENTS.md"), encoding="utf-8") as fh:
            agents = fh.read()
        self.assertIn("## Reviewing this analysis", agents)
        self.assertIn("`logs/conversation/conversation.md` (the conversation, readable", agents)
        self.assertLess(agents.index("## Reviewing this analysis"), agents.index("## Do not"))

    def test_the_delivery_never_ships_it(self):
        import core_submission as cs
        out = os.path.join(self.d, "output")
        for rel, text in (("Analysis_Report.html", "<html>r</html>"),
                          ("reproducibility/conversation/" + SID + ".jsonl", "{}"),
                          ("reproducibility/conversation/conversation.md", "# c"),
                          ("reproducibility/" + SID2 + ".jsonl", "{}"),
                          ("tables/decisions.md", "# d"),
                          ("tables/DE_a.csv", "x")):
            path = os.path.join(out, rel)
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as fh:
                fh.write(text)
        items = cs.plan_delivery(out)
        by = {i["name"]: i for i in items}
        for name in ("reproducibility/conversation/", "reproducibility/" + SID2 + ".jsonl",
                     "tables/decisions.md"):
            self.assertEqual((by[name]["status"], by[name]["reason"]),
                             ("SKIPPED", cs.INTERNAL_REASON), name)
        self.assertFalse([i for i in items if i["status"] == "OK" and
                          ("conversation" in i["name"] or i["name"].endswith(".jsonl"))])
        self.assertEqual(by["tables/DE_a.csv"]["status"], "OK")

    def test_the_reviewer_checklist(self):
        import session_docs
        import test_deposit_package as tdp
        p = tdp.dia_session(os.path.join(self.d, "s"))
        text = "\n".join(session_docs.reviewing_lines(session_docs.gather(p["session_dir"])))
        for want in ("## Reviewing this analysis",
                     "- **Groups and contrasts:** were the groups and contrasts analysed the ones "
                     "the user confirmed?",
                     "- **Every number traces to a command:**",
                     "- **Nothing unrecorded:** was any step skipped, changed or re-run without "
                     "being recorded",
                     "- **No warning ignored:** were any warnings ignored",
                     "a record that is missing is a finding, not a pass"):
            self.assertIn(want, text)
        import log_decision
        log_decision.append(p["session_dir"], "x", "y")
        text = "\n".join(session_docs.reviewing_lines(session_docs.gather(p["session_dir"])))
        self.assertIn("`logs/decisions.md` (the decisions, and why)", text)

    def test_the_registry_points_to_it_and_copies_nothing(self):
        import record_run
        sess = os.path.join(self.d, "s2")
        conv = os.path.join(sess, "logs", "conversation")
        os.makedirs(conv)
        for f in (SID + ".jsonl", SID2 + ".jsonl", "conversation.md"):
            with open(os.path.join(conv, f), "w") as fh:
                fh.write("x")
        got = record_run._conversation(sess)
        self.assertEqual((got["dir"], got["n"]), (conv, 2))
        lines = "\n".join(record_run.render_analysis({}, {
            "conversation": got, "decisions": os.path.join(sess, "logs", "decisions.md")}))
        self.assertIn(f"**The analysis conversation (Core-internal, not copied):** `{conv}` -- "
                      "2 conversation(s), redacted", lines)
        self.assertIn("**Decisions log:**", lines)
        self.assertIsNone(record_run._conversation(os.path.join(self.d, "nowhere")))

    def test_the_catalog_describes_it(self):
        import make_report
        for name in ("conversation.md", SID + ".jsonl", "decisions.md"):
            cat, desc = make_report.describe(name)
            self.assertNotIn("unrecognized", desc, name)
        self.assertIn("CORE-INTERNAL", make_report.describe("conversation.md")[1])


if __name__ == "__main__":
    unittest.main()
