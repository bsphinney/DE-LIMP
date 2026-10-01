#!/usr/bin/env python3
"""
save_transcript.py -- keep the conversation that ran an analysis, redacted, in the session, for
an AI (or a person) to review later: another reproducibility check (Brett, 2026-09-28: "save a
transcript of the data analysis session for an AI to look over later").

  python3 save_transcript.py <session dir>            this Claude Code conversation, now
  python3 save_transcript.py --hive '<HIVE session>'  a session on HIVE driven from this computer:
                                                     redacted here, then hive_exec.sh --put
  python3 save_transcript.py <session dir> --remember
                                                     only record this conversation for the hook
  python3 save_transcript.py <session dir> --transcript <file.jsonl> --session-id <id>
                                                     what transcript_hook.py runs

What it writes, in <session>/logs/conversation/ -- CORE-INTERNAL: conversations hold internal
paths, other projects and offhand remarks, so deliver never ships it and the session zip leaves
it out (it stays on disk; the Core's run registry points to it):
  <session-id>.jsonl  the Claude Code transcript, REDACTED; one file per conversation (a resume or
                      a new chat is another), replaced when the same one is saved again
  index.json          per conversation: id, first/last timestamp, entries, bytes, sha256, when
                      saved, the Claude Code version, how many strings were redacted
  conversation.md     the conversations in order, readable: the user's messages, the
                      assistant's text, each tool call (the command, its result cut to ~40
                      lines), with timestamps. Best effort: the raw file is the record

Finding it: Claude Code sets CLAUDE_CODE_SESSION_ID in every Bash subprocess, and keeps the
transcript at <config>/projects/<cwd with non-alphanumerics as ->/<session id>.jsonl, where
<config> is $CLAUDE_CONFIG_DIR or ~/.claude -- found by the id, whatever the project folder is.
Claude Code deletes transcripts after cleanupPeriodDays (30 by default): this copy stays. The
format is internal to Claude Code and changes between versions, so nothing here trusts it: an
entry not understood is one "[unrecognised entry]" line in conversation.md, never a crash.

Redaction (Redactor: the one redaction for the analysis record, which log_decision.py uses
too), before anything is written, of every string -- and every key -- in every entry: the
skill's one list of secret patterns (notify_slack.redact: Google/Gemini keys, GitHub and Slack
tokens and webhooks, key=, token=, password=, bearer, JWT, private keys ...); then the secret
values this computer holds (the CoreOmics token, the Gemini key and Slack webhook files, and any
environment variable named like a token, key, secret, password, webhook or cookie), also when a
copy is wrapped across lines or base64-encoded; then assignments by such names in env or JSON
form (FOO_TOKEN=..., "api_key": "...", accessToken: ...), and any field so named. These last are
here, not in the one list, because they are broad: record_run.py REFUSES whatever that list
matches, and a transcript only needs its copy scrubbed. Images are replaced by "[image]" (a
screenshot can show anything). No raw copy is ever written. What it CANNOT catch: a secret in
ordinary prose ("my password is ..."), or one this computer does not hold and no pattern
matches -- which is why the conversation is Core-internal, never a deliverable.

Files are 0640 in a 0750 folder (the Core's group reads them on HIVE). One conversation that
runs several analyses is split: see remember().

Other agents (not Claude Code): there is no transcript to find, so it says so and exits 0 --
the decisions log (log_decision.py -> logs/decisions.md) is the record there, and a curated
one is easier to review on Claude Code too.

Never fails a step: no Claude Code, no transcript, nothing writable -> one JSON line, exit 0.
"""
import argparse
import base64
import datetime
import glob
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

CONV_DIR = os.path.join("logs", "conversation")
SKILL_NAME = "ucdavis-proteomics-core-pipeline"
OMITTED = "[other work in this conversation omitted]"
INDEX = "index.json"
READABLE = "conversation.md"
RESULT_LINES = 40                        # a tool result in conversation.md: first 30, last 10
DIR_MODE, FILE_MODE = 0o750, 0o640        # the Core's group reads it on HIVE; nobody else
MAP_KEEP_DAYS = 180                      # the hook's map forgets older sessions
_SAFE_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_\-]{7,79}$")   # a session id, never a glob
# an environment variable named like a secret: its value is one
_SECRET_NAME = re.compile(r"(?i)(token|secret|passw(or)?d|api_?key|webhook|cookie|credential)")
# a field or key that holds one: access_token, api_key, password, webhook_url -- not
# input_tokens or tokenizer
# (camelCase too: accessToken, refreshToken, clientSecret, apiKey)
_KEY = (r"[A-Za-z0-9_\-]*?(?:token|secret|passw(?:or)?d|api[_\-]?key|webhook(?:[_\-]?url)?|"
        r"cookie|credentials?)")
_SECRET_KEY = re.compile(r"(?i)^" + _KEY + r"$")
# FOO_TOKEN=value (env style), api_key: value -- a value that starts like a path, a variable or
# an already redacted mark is a name, not a secret
_VALUE = r"(\s*[=:]\s*)(['\"]?)(?![/~.$\[])([^\s'\"]{8,})"
_ASSIGN = [re.compile(r"\b([A-Z0-9_]*(?:TOKEN|SECRET|PASSWORD|PASSWD|API_?KEY|WEBHOOK|COOKIE)"
                      r"[A-Z0-9_]*)" + _VALUE),
           re.compile(r"(?i)(?<![A-Za-z0-9_])(" + _KEY + r")(?![A-Za-z0-9_])" + _VALUE)]
_JSON_KV = re.compile(r"(?i)(\"" + _KEY + r"\"\s*:\s*\")([^\"]{8,})(\")")

# The entry types a transcript carries that are not conversation: shown as a count, not a line
_METADATA = {"attachment", "last-prompt", "mode", "permission-mode", "atis-latch", "ai-title",
             "queue-operation", "file-history-snapshot", "file-history-delta", "pr-link",
             "frame-link", "bridge-session", "custom-title", "tag", "agent-name", "cost-state"}


def now_iso():
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()


ONELINE = False                          # --oneline: the hook's log gets one line per save


def emit(obj, quiet=False):
    if not quiet:
        print(f"{now_iso()} save_transcript {json.dumps(obj)}" if ONELINE else
              json.dumps(obj, indent=2), flush=True)
    return obj


# --------------------------------------------------------------------- the hook's session map
def config_dir():
    return os.path.expanduser(os.environ.get("SKILL_CONFIG_DIR")
                              or os.path.join("~", ".config", "ucdavis-proteomics"))


def map_path():
    return os.path.join(config_dir(), "transcript_sessions.json")


def hook_log_path():
    return os.path.join(config_dir(), "transcript_hook.log")


def _read_map():
    try:
        with open(map_path(), encoding="utf-8") as fh:
            m = json.load(fh)
        return m if isinstance(m, dict) else {}
    except (OSError, ValueError):
        return {}


def _skey(session_dir=None, hive=None):
    return f"hive:{hive}" if hive else f"local:{os.path.abspath(session_dir or '')}"


def _record(m, sid):
    """A conversation's map record: {"sessions": {key: {session_dir, hive, recorded}},
    "segments": [{"session": key, "start": marker}]} -- also from the earlier shapes."""
    rec = m.get(sid) if isinstance(m.get(sid), dict) else {}
    if isinstance(rec.get("segments"), list) and isinstance(rec.get("sessions"), dict):
        return {"sessions": rec["sessions"],
                "segments": [s for s in rec["segments"] if isinstance(s, dict)]}
    olds = rec.get("sessions") if isinstance(rec.get("sessions"), list) else (
        [dict(rec, start={"line": 0, "after_uuid": None})]
        if rec.get("session_dir") or rec.get("hive") else [])
    sessions, segments = {}, []
    for s in sorted((s for s in olds if isinstance(s, dict)),
                    key=lambda s: int((s.get("start") or {}).get("line") or 0)):
        k = _skey(s.get("session_dir"), s.get("hive"))
        sessions[k] = {"session_dir": s.get("session_dir"), "hive": s.get("hive"),
                       "recorded": s.get("recorded")}
        segments.append({"session": k, "start": s.get("start") or {"line": 0,
                                                                    "after_uuid": None}})
    return {"sessions": sessions, "segments": segments}


def lookup(session_id):
    """The analysis sessions a Claude Code session has run, in the order they first appear
    ([] when none) -- the hook saves each one's own part."""
    if not session_id:
        return []
    rec = _record(_read_map(), session_id)
    keys = list(dict.fromkeys(s["session"] for s in rec["segments"]))
    return [dict(rec["sessions"].get(k) or {}, key=k) for k in keys if k in rec["sessions"]]


def _transcript_end(path):
    """(complete lines, the uuid of the last entry) -- where a segment opened now starts."""
    n, last = 0, None
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for raw in fh:
                if not raw.endswith("\n"):
                    break                                   # an entry still being written
                n += 1
                if '"uuid"' in raw:
                    try:
                        u = json.loads(raw).get("uuid")
                        last = u if isinstance(u, str) else last
                    except (ValueError, AttributeError):
                        pass
    except (OSError, TypeError):
        pass
    return n, last


def _now(path):
    n, last = _transcript_end(path) if path else (0, None)
    return {"line": n, "after_uuid": last}


def _at(lines, i):
    """The marker that starts a segment AT line i (after the entry before it)."""
    u = None
    if i > 0:
        try:
            u = json.loads(lines[i - 1]).get("uuid")
        except (ValueError, AttributeError):
            u = None
    return {"line": i, "after_uuid": u if isinstance(u, str) else None}


def _user_text(e):
    """The user's own words in an entry (not a tool result, the skill's text or a summary)."""
    if not isinstance(e, dict) or e.get("type") != "user" or e.get("isMeta") \
            or e.get("isCompactSummary"):
        return None
    c = (e.get("message") or {}).get("content")
    if isinstance(c, str):
        return c
    texts = [str(b.get("text") or "") for b in c if isinstance(b, dict) and b.get("type") == "text"] \
        if isinstance(c, list) else []
    return "\n".join(texts) if texts else None


def _opening_request(lines, load):
    """When Claude invoked the skill itself (a Skill tool call), the user's request that led to
    it ("analyse PROT_0802, Old vs Young") -- but only when the Skill call DIRECTLY answers it:
    the assistant's first turn after the request, with no tool call and no tool result between
    (a tool's output there -- another client's results -- is not this analysis's). Else None. A
    slash command or an injected load carries the request itself."""
    try:
        e = json.loads(lines[load])
    except (ValueError, IndexError):
        return None
    if not isinstance(e, dict) or e.get("type") != "assistant":
        return None
    for i in range(load - 1, -1, -1):
        try:
            e = json.loads(lines[i])
        except ValueError:
            continue
        if not isinstance(e, dict) or e.get("type") not in ("user", "assistant"):
            continue                                   # metadata between: no content
        c = (e.get("message") or {}).get("content")
        blocks = [b for b in c if isinstance(b, dict)] if isinstance(c, list) else []
        if e.get("type") == "assistant":
            if any(b.get("type") == "tool_use" for b in blocks):
                return None                            # a tool ran first: not a direct answer
            continue                                   # the same turn's text or thinking
        if _user_text(e):
            return i
        return None                                    # a tool result or injected text
    return None


def _skill_load(lines):
    """The first entry that loads this skill -- the Skill tool called for it, its slash command,
    or the skill's text injected ("Base directory for this skill") -- or None."""
    for i, raw in enumerate(lines):
        if SKILL_NAME not in raw:
            continue
        try:
            e = json.loads(raw)
            c = (e.get("message") or {}).get("content")
        except (ValueError, AttributeError):
            continue
        for b in c if isinstance(c, list) else []:
            if (isinstance(b, dict) and b.get("type") == "tool_use" and b.get("name") == "Skill"
                    and SKILL_NAME in json.dumps(b.get("input"))):
                return i
        if e.get("type") == "user":
            text = _content_text(c) if c is not None else ""
            if SKILL_NAME in text and ("<command-name>" in text or "Base directory for this skill"
                                       in text or e.get("isMeta")):
                return i
    return None


def _spoken(lines):
    """What the user and the assistant wrote, and the tools' inputs -- not tools' output, which
    lists other submissions (locate, a folder listing) without being about them."""
    out = []
    for raw in lines:
        try:
            e = json.loads(raw)
            c = (e.get("message") or {}).get("content")
        except (ValueError, AttributeError):
            continue
        if e.get("type") == "user" and not e.get("isMeta"):
            if isinstance(c, str):
                out.append(c)
            for b in c if isinstance(c, list) else []:
                if isinstance(b, dict) and b.get("type") == "text":
                    out.append(str(b.get("text") or ""))
        elif e.get("type") == "assistant":
            for b in c if isinstance(c, list) else []:
                if isinstance(b, dict) and b.get("type") == "text":
                    out.append(str(b.get("text") or ""))
                elif isinstance(b, dict) and b.get("type") == "tool_use":
                    out.append(json.dumps(b.get("input"), ensure_ascii=False))
    return "\n".join(out)


def _own_prots(session_dir, hive):
    """The PROT numbers that are this session's: in its path, or its attached record."""
    import core_submission as cs
    own = {v for k, v in cs.ids_in(hive or session_dir or "") if k == "internal_id"}
    if session_dir:
        try:
            import submission_report
            core = submission_report.core_run(session_dir)
            if core and core.get("prot"):
                own.add(core["prot"])
        except Exception:
            pass
    return own


def _first_start(path, session_dir, hive, from_now):
    """Where the FIRST session of a conversation starts: at the skill's first load -- the design
    is agreed before `init` -- unless the talk between there and now names a submission that is
    not this session's (another client), or the load cannot be found: then here, with a warning.
    `from_now`: here, on purpose (an orchestrator, a conversation about several clients)."""
    if from_now or not path:
        return _now(path), None
    lines = _complete_lines(path)
    load = _skill_load(lines)
    if load is None:
        return _now(path), ("the skill's first load is not in this conversation, so its copy "
                            "starts here")
    import core_submission as cs
    own = _own_prots(session_dir, hive)
    named = lambda part: {v for k, v in cs.ids_in(_spoken(part), typed=True)  # noqa: E731
                          if k == "internal_id"} - own
    start = load
    ask = _opening_request(lines, load)          # an auto-invoked skill: the user's request
    if ask is not None and not named(lines[ask:ask + 1]):
        start = ask
    other = sorted(named(lines[start:]))
    if other:
        return _now(path), (f"the conversation before this session names {', '.join(other)}, "
                            "not this session's submission, so its copy starts here, not at the "
                            "skill's first load")
    return _at(lines, start), None


def _label(info):
    return (info or {}).get("hive") or (info or {}).get("session_dir") or "?"


def remember(session_dir=None, hive=None, session_id=None, transcript=None, switch=True,
             from_now=False):
    """Record that this Claude Code conversation (CLAUDE_CODE_SESSION_ID) is now working on the
    analysis in `session_dir` (or the HIVE session `hive`), so the hook saves it there.

    A conversation is a run of SEGMENTS, one analysis each. The first opens at the skill's first
    load (_first_start). Recording a session that is not the one being worked on closes that
    segment and opens one for it, here -- an earlier session too (X, Y, back to X). Recording the
    active session again changes nothing. Each session's copy is its own segments, joined with
    "[other work in this conversation omitted]". switch=False (checkpoint.py status without
    --resume, a save): never a switch -- only a conversation that runs no session yet is
    recorded. from_now: the session starts HERE (a new segment, or the active one restarted).
    -> (the session's entry or None, a warning or None); (None, None) outside Claude Code."""
    sid = session_id or os.environ.get("CLAUDE_CODE_SESSION_ID")
    if not sid or not _SAFE_ID.match(sid):
        return None, None
    m = _read_map()
    cutoff = (datetime.datetime.now(datetime.timezone.utc)
              - datetime.timedelta(days=MAP_KEEP_DAYS)).isoformat()
    m = {k: v for k, v in m.items() if isinstance(v, dict)
         and str(v.get("updated") or v.get("recorded")) >= cutoff}
    rec = _record(m, sid)
    key = _skey(session_dir, hive)
    me = {"session_dir": os.path.abspath(session_dir) if session_dir else None,
          "hive": hive or None}
    segs, warn = rec["segments"], None
    path = transcript or find_transcript(sid)
    if not segs:
        start, warn = _first_start(path, session_dir, hive, from_now)
        segs.append({"session": key, "start": start})
    elif segs[-1]["session"] == key:
        if from_now:
            segs[-1]["start"] = _now(path)
            warn = "its copy now starts here (--from-now)"
    elif switch or from_now:
        prev = rec["sessions"].get(segs[-1]["session"])
        segs.append({"session": key, "start": _now(path)})
        warn = (f"this conversation moves here from {_label(prev)} to {_label(me)}: each copy "
                "keeps only its own part")
    else:
        prev = rec["sessions"].get(segs[-1]["session"])
        info = rec["sessions"].get(key)
        return (dict(info, key=key) if info else None,
                f"not switched: this conversation is working on {_label(prev)} "
                "(checkpoint.py status --resume, or save_transcript.py --remember, switches)")
    rec["sessions"][key] = dict(me, recorded=now_iso())
    m[sid] = dict(rec, updated=now_iso())
    os.makedirs(config_dir(), exist_ok=True)
    part = f"{map_path()}.{os.getpid()}.part"
    with open(part, "w", encoding="utf-8") as fh:
        json.dump(m, fh, indent=1, sort_keys=True)
    os.replace(part, map_path())
    if warn:
        sys.stderr.write(f"[save_transcript] {warn}\n")
    return dict(rec["sessions"][key], key=key), warn


def _parts(sid, key):
    """This session's segments of the conversation: [(start, end)], end None for the last."""
    segs = _record(_read_map(), sid)["segments"]
    return [(s["start"], segs[i + 1]["start"] if i + 1 < len(segs) else None)
            for i, s in enumerate(segs) if s["session"] == key]


# ------------------------------------------------------------------------------- finding it
def find_transcript(session_id):
    """<config>/projects/*/<session id>.jsonl -- the newest, if several. None when absent."""
    if not session_id or not _SAFE_ID.match(session_id):
        return None
    roots = [os.environ.get("CLAUDE_CONFIG_DIR"), os.path.join("~", ".claude")]
    found = []
    for r in roots:
        if r:
            found += glob.glob(os.path.join(os.path.expanduser(r), "projects", "*",
                                            session_id + ".jsonl"))
    found = [f for f in found if os.path.isfile(f)]
    return max(found, key=os.path.getmtime) if found else None


# -------------------------------------------------------------------------------- redaction
def _looks_like_path(v):
    return v.startswith(("/", "~", "./", "../", "$")) or bool(re.match(r"^[A-Za-z]:[\\/]", v))


def secret_values():
    """The secrets this computer holds, as exact strings, longest first: the files the skill
    keeps them in and every environment variable named like one."""
    vals = set()
    for k, v in os.environ.items():
        if _SECRET_NAME.search(k) and v and len(v) >= 8 and not _looks_like_path(v):
            vals.add(v.strip())
    cfg = config_dir()
    files = [os.path.join("~", ".coreomics_token"), os.path.join(cfg, "gemini_key"),
             os.path.join(cfg, "slack_webhook")]
    files += glob.glob(os.path.join(cfg, "*token*")) + glob.glob(os.path.join(cfg, "*secret*"))
    # The CoreOmics key wherever core_submission may find it (on Windows also Git Bash's HOME),
    # and decoded as it decodes it: a UTF-16 file read line by line below matches nothing. Any
    # failure here only loses these extra places -- it must never stop the redaction.
    try:
        import core_submission as cs
        for f in cs.key_files() + cs.stray_key_files():
            files.append(f)
            try:
                vals.add(cs.read_key_file(f))
            except OSError:
                pass
    except Exception:
        pass
    for f in files:
        try:
            with open(os.path.expanduser(f), encoding="utf-8", errors="replace") as fh:
                for line in fh.read().splitlines():
                    line = line.strip()
                    if len(line) >= 8 and not line.startswith("#"):
                        vals.add(line.split("=", 1)[1].strip().strip("'\"")
                                 if re.match(r"^[A-Z_]+=", line) else line)
        except OSError:
            continue
    return sorted((v for v in vals if len(v) >= 8), key=len, reverse=True)


def _b64_forms(value):
    """The value as it appears inside base64 (standard and URL-safe), at each of the 3 byte
    alignments: only the characters that depend on the value alone."""
    raw, forms = value.encode("utf-8"), set()
    for k in range(3):
        enc = base64.b64encode(b"\0" * k + raw).decode("ascii")
        a, b = -(-8 * k // 6), 8 * (k + len(raw)) // 6
        part = enc[a:b]
        if len(part) >= 12:
            forms.update((part, part.replace("+", "-").replace("/", "_")))
    return forms


def _known_patterns(values):
    """One pattern per secret value: the value with any whitespace between its characters (a
    copy wrapped across lines), or any of its base64 forms."""
    out = []
    for v in values:
        alts = [r"\s*".join(re.escape(c) for c in v)] + sorted(re.escape(f) for f in _b64_forms(v))
        out.append(re.compile("|".join(alts)))
    return out


class Redactor(object):
    """THE redaction for the analysis record -- the saved conversation and the decisions log
    (log_decision.py) both use it."""

    def __init__(self, known=None):
        from notify_slack import REDACTED, redact
        self.mark, self._one_list = REDACTED, redact
        self.known = _known_patterns(secret_values() if known is None else known)
        self.n = 0                                      # strings changed

    def text(self, s):
        t = s
        for rx in self.known:
            t = rx.sub(self.mark, t)
        t = self._one_list(t)
        for rx in _ASSIGN:
            t = rx.sub(lambda m: m.group(1) + m.group(2) + m.group(3) + self.mark, t)
        t = _JSON_KV.sub(lambda m: m.group(1) + self.mark + m.group(3), t)
        if t != s:
            self.n += 1
        return t

    def obj(self, o):
        if isinstance(o, str):
            return self.text(o)
        if isinstance(o, list):
            return [self.obj(x) for x in o]
        if isinstance(o, dict):
            if o.get("type") == "image":                # a screenshot can hold anything
                self.n += 1
                return {"type": "text", "text": "[image]"}
            out = {}
            for k, v in o.items():
                k2 = self.text(k) if isinstance(k, str) else k
                if isinstance(v, str) and len(v) >= 8 and _SECRET_KEY.match(str(k)):
                    self.n += 1
                    out[k2] = self.mark
                else:
                    out[k2] = self.obj(v)
            return out
        return o

    def line(self, raw):
        """One JSONL line: parsed and every string redacted (the structure survives), or, when
        it is not JSON, the line as text."""
        try:
            return json.dumps(self.obj(json.loads(raw)), ensure_ascii=False)
        except ValueError:
            return self.text(raw.rstrip("\n"))


# -------------------------------------------------------------------------------- saving
def _write(path, data, mode=FILE_MODE):
    part = f"{path}.{os.getpid()}.part"         # the hook and an explicit save may overlap
    with open(part, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(data)
    _chmod(part, mode)
    os.replace(part, path)


def _chmod(path, mode):
    try:
        os.chmod(path, mode)
    except OSError:
        pass                                    # a filesystem without modes: nothing to do


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def _stats(path):
    """first/last timestamp, entries, Claude Code version -- read defensively."""
    first = last = version = None
    n = 0
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if not raw.strip():
                continue
            n += 1
            try:
                e = json.loads(raw)
            except ValueError:
                continue
            if not isinstance(e, dict):
                continue
            ts = e.get("timestamp")
            if isinstance(ts, str) and ts:
                first = first or ts
                last = ts
            if isinstance(e.get("version"), str):
                version = e["version"]
    return {"first": first, "last": last, "entries": n, "claude_code_version": version}


def _complete_lines(path):
    """The transcript's lines, less a last one with no newline: an entry still being written
    (a PreCompact save racing Claude Code's append)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    lines = text.split("\n")
    return lines[:-1]                    # after a final "\n" that is "", else the partial line


def _index_of(lines, marker, default):
    if not marker:
        return default
    u = marker.get("after_uuid")
    if u:
        needle = json.dumps(u)
        for i, raw in enumerate(lines):
            if needle in raw:
                try:
                    if json.loads(raw).get("uuid") == u:
                        return i + 1
                except (ValueError, AttributeError):
                    pass
    line = marker.get("line")
    return min(int(line), len(lines)) if isinstance(line, int) else default


def _omitted(lines, i):
    """The line that stands for other work between two of a session's segments."""
    ts = None
    for raw in lines[i:i + 1]:
        try:
            ts = json.loads(raw).get("timestamp")
        except (ValueError, AttributeError):
            pass
    return json.dumps({"type": "system", "subtype": "omitted", "content": OMITTED,
                       "timestamp": ts})


SUMMARY_OMITTED = ("[compaction summary omitted: it covers conversation outside this analysis's "
                   "record]")


def _allowed(e):
    """What a copy that is not the whole conversation may hold -- an ALLOWLIST, since Claude Code's
    format is internal and keeps adding entry types (custom-title, last-prompt, queue-operation,
    attachment, a system recap ...) that repeat or sum up what came before: "keep" (the user's
    and the assistant's messages, tool results, a compaction boundary), "summary" (a compaction
    summary: kept only when the copy holds all that it sums up), or "other"."""
    if not isinstance(e, dict):
        return "other"
    if e.get("isCompactSummary"):
        return "summary"
    t = e.get("type")
    if t in ("user", "assistant") and isinstance(e.get("message"), dict) and not e.get("isMeta"):
        return "keep"
    if t == "system" and e.get("subtype") == "compact_boundary":
        return "keep"
    return "other"


def _marker(subtype, content, ts=None):
    return json.dumps({"type": "system", "subtype": subtype, "content": content, "timestamp": ts})


def _select(lines, ranges):
    """The entries of this session's ranges, as they go into its copy. A copy of the WHOLE
    conversation (one range to the end, with no message before it -- only the metadata Claude
    Code writes first) keeps everything. Any other copy keeps only what _allowed allows, each run
    of other entries as one counted line, and a compaction summary only when the copy holds every
    message before it (it falls in the first range, and no message precedes that) -- a summary
    sums up all that came before: another client before --from-now or before the skill load,
    another analysis."""
    def talk_before(a):                    # a message before line a: the copy misses some
        for raw in lines[:a]:
            try:
                if _allowed(json.loads(raw)) != "other":
                    return True
            except ValueError:
                continue
        return False
    first_all = bool(ranges) and not talk_before(ranges[0][0])   # only metadata before it
    whole = len(ranges) == 1 and ranges[0][1] == len(lines) and first_all
    out = []
    for k, (a, b) in enumerate(ranges):
        if out and a > 0:
            out.append(_omitted(lines, a))
        run = 0
        for raw in lines[a:b]:
            if not raw.strip():
                continue
            if whole:
                out.append(raw)
                continue
            try:
                e = json.loads(raw)
            except ValueError:
                e = None
            kind = _allowed(e)
            if kind == "summary" and k == 0 and first_all:
                kind = "keep"                          # it sums up only what this copy holds
            if kind == "other":
                run += 1
                continue
            if run:
                out.append(_marker("omitted-entries", f"[{run} other entr"
                                                      f"{'y' if run == 1 else 'ies'} omitted]"))
                run = 0
            out.append(raw if kind == "keep" else
                       _marker("summary-omitted", SUMMARY_OMITTED, e.get("timestamp")))
        if run:
            out.append(_marker("omitted-entries",
                               f"[{run} other entr{'y' if run == 1 else 'ies'} omitted]"))
    return out


def _fingerprint(parts):
    """The selection a copy was made from: this session's segments (their start and end markers).
    It changes with --from-now, a switch, or any other change to the session's segments."""
    return hashlib.sha256(json.dumps(parts, sort_keys=True).encode("utf-8")).hexdigest()[:16]


def save_into(conv_dir, transcript, session_id, parts=None, redactor=None, current=None):
    """Redact this session's parts of `transcript` -- [(start, end)] markers, joined with an
    "[other work in this conversation omitted]" line; None: all of it -- into
    conv_dir/<id>.jsonl; then rebuild index.json and conversation.md from every .jsonl there.
    What a partial copy may hold: _select.

    WHICH COPY WINS -- by how far into the transcript a copy REACHES and what it was selected
    from, never by how many entries it holds (a narrower selection, or the allowlist after a
    switch, is shorter and still newer). The index records each copy's `selection` (the
    segment markers' fingerprint) and `reach` (the transcript lines it covers, and the uuid of
    the last). A save whose selection no longer matches the map (`current`: the map's selection
    now) started before the segments changed: it writes nothing. Otherwise a different selection
    from the copy's always rewrites it (--from-now, a switch); the same one rewrites when it
    reaches as far or further, so an older snapshot finishing after a newer save is refused.
    Both checks run again after the redaction, just before the write: a newer save that
    finished while this one redacted wins.
    -> {"entry", "written", "why"}."""
    os.makedirs(conv_dir, exist_ok=True)
    _chmod(conv_dir, DIR_MODE)
    r = redactor or Redactor()
    lines = _complete_lines(transcript)
    ranges = []
    for start, end in parts if parts is not None else [(None, None)]:
        a = _index_of(lines, start, 0)
        ranges.append((a, max(a, _index_of(lines, end, len(lines)))))
    fp = _fingerprint(parts)
    reach = max((b for a, b in ranges), default=0)
    uuid = None
    if reach:
        try:
            uuid = json.loads(lines[reach - 1]).get("uuid")
        except (ValueError, AttributeError):
            uuid = None
    dest = os.path.join(conv_dir, session_id + ".jsonl")

    def refused():
        old = next((x for x in read_index(conv_dir) if x.get("id") == session_id), {})
        if current is not None and current() != fp:
            return "the session's part of the conversation changed while this save ran: kept"
        if (os.path.isfile(dest) and old.get("selection") == fp
                and isinstance(old.get("reach"), int) and reach < old["reach"]):
            return "the copy there reaches further into the conversation (a newer save): kept it"
        return None
    why = refused()
    if why is None:
        out = [r.line(raw) for raw in _select(lines, ranges)]
        why = refused()                     # read again: the redaction is the slow part
        if why is None:
            _write(dest, "\n".join(out) + ("\n" if out else ""))
    written = why is None
    index = rebuild_index(conv_dir, {session_id: {
        "redacted_strings": r.n, "selection": fp,
        "reach": reach, "reach_uuid": uuid if isinstance(uuid, str) else None}}
        if written else {})
    entry = next((x for x in index if x["id"] == session_id), {})
    return {"entry": entry, "written": written, "why": why}


def rebuild_index(conv_dir, fresh=None):
    """index.json and conversation.md from the .jsonl files that are there -- never from a list
    that a failed fetch or an overlapping save could have left short. `fresh`: {id: {strings
    redacted, selection, reach, reach_uuid}} for the file this save wrote; the others keep what
    the index said, while their sha256 still matches."""
    fresh = fresh or {}
    old = {x.get("id"): x for x in read_index(conv_dir)}
    index = []
    for f in sorted(glob.glob(os.path.join(conv_dir, "*.jsonl"))):
        sid = os.path.basename(f)[:-len(".jsonl")]
        sha = _sha256(f)
        prev = old.get(sid) if (old.get(sid) or {}).get("sha256") == sha else None
        keep = fresh.get(sid) or {k: (prev or {}).get(k) for k in
                                  ("redacted_strings", "selection", "reach", "reach_uuid")}
        index.append(dict(_stats(f), id=sid, file=os.path.basename(f),
                          bytes=os.path.getsize(f), sha256=sha,
                          saved_at=(now_iso() if sid in fresh or not prev
                                    else prev.get("saved_at")), **keep))
    index.sort(key=lambda x: (x.get("first") or "", x.get("id") or ""))
    _write(os.path.join(conv_dir, INDEX), json.dumps(
        {"schema": "conversation_index/1",
         "note": "Claude Code transcripts of this analysis, redacted (save_transcript.py). "
                 "Core-internal: never delivered, never in the session zip.",
         "conversations": index}, indent=2) + "\n")
    _write(os.path.join(conv_dir, READABLE), render(conv_dir, index))
    return index


def read_index(conv_dir):
    try:
        with open(os.path.join(conv_dir, INDEX), encoding="utf-8") as fh:
            got = json.load(fh).get("conversations")
        return [x for x in got if isinstance(x, dict)] if isinstance(got, list) else []
    except (OSError, ValueError, AttributeError):
        return []


def saved_count(session_dir):
    d = os.path.join(session_dir, CONV_DIR)
    return len(glob.glob(os.path.join(d, "*.jsonl"))) if os.path.isdir(d) else 0


def save(session_dir=None, transcript=None, session_id=None, hive=None, quiet=True,
         from_now=False):
    """The one entry point (CLI, hook, session.py finalize). -> {"saved": ...} or
    {"skipped": why}; never raises for a missing transcript."""
    sid = session_id or os.environ.get("CLAUDE_CODE_SESSION_ID")
    if sid and not _SAFE_ID.match(sid):
        return emit({"skipped": f"not a session id: {sid[:40]!r}"}, quiet)
    if not hive and not (session_dir and os.path.isdir(session_dir)):
        return emit({"skipped": f"no session folder {session_dir!r}"}, quiet)
    if not transcript:
        if not sid:
            return emit({"skipped": "not a Claude Code session (CLAUDE_CODE_SESSION_ID is not "
                                    "set): record decisions with log_decision.py"}, quiet)
        transcript = find_transcript(sid)
        if not transcript:
            return emit({"skipped": f"no transcript for Claude Code session {sid} under "
                                    "$CLAUDE_CONFIG_DIR or ~/.claude/projects"}, quiet)
    if not os.path.isfile(transcript):
        return emit({"skipped": f"transcript not found: {os.path.basename(transcript)}"}, quiet)
    sid = sid or os.path.splitext(os.path.basename(transcript))[0]
    if not _SAFE_ID.match(sid):
        return emit({"skipped": f"not a session id: {sid[:40]!r}"}, quiet)
    try:                                          # a save never switches (see remember)
        entry, warn = remember(session_dir, hive, sid, transcript, switch=from_now,
                               from_now=from_now)
    except OSError as e:
        return emit({"skipped": f"the hook's map could not be read or written "
                                f"({e.strerror or e}): nothing saved"}, quiet)
    key = _skey(session_dir, hive)
    parts = _parts(sid, key)
    current = lambda: _fingerprint(_parts(sid, key))           # noqa: E731
    if not parts:
        return emit({"skipped": warn or "this session is not recorded for this conversation"},
                    quiet)
    if hive:
        res = _save_hive(hive, transcript, sid, parts, current)
    else:
        conv = os.path.join(os.path.abspath(session_dir), CONV_DIR)
        got = save_into(conv, transcript, sid, parts, current=current)
        res = {"saved": os.path.join(conv, sid + ".jsonl"), "session_id": sid,
               "entries": got["entry"].get("entries"),
               "redacted_strings": got["entry"].get("redacted_strings"),
               "conversations": saved_count(session_dir),
               "readable": os.path.join(conv, READABLE)}
        if not got["written"]:
            res["kept"] = got["why"]
    if warn:
        res["warning"] = warn
    return emit(res, quiet)


def _hive(*args):
    exe = os.environ.get("SKILL_HIVE_EXEC") or os.path.join(HERE, "hive_exec.sh")
    return subprocess.run(["bash", exe] + list(args), capture_output=True, text=True,
                          timeout=600)


def _tail(r):
    return (r.stderr or r.stdout or "").strip()[-200:]


def _save_hive(remote, transcript, sid, parts=None, current=None):
    """A HIVE session driven from this computer: redact HERE (the raw transcript never leaves),
    fetch what is saved there -- telling "nothing there yet" from "could not fetch", which aborts
    before anything on HIVE changes -- rebuild the index and readable file from the files, and
    put them back."""
    rconv = remote.rstrip("/") + "/" + CONV_DIR.replace(os.sep, "/")
    q = shlex.quote(rconv)
    probe = _hive(f"if test -d {q}; then echo __there__; else echo __absent__; fi")
    said = (probe.stdout or "").split()
    if probe.returncode != 0 or not ({"__there__", "__absent__"} & set(said)):
        return {"skipped": f"cannot reach HIVE to check {rconv}: {_tail(probe)}"}
    with tempfile.TemporaryDirectory(prefix="conversation-") as tmp:
        local = os.path.join(tmp, "conversation")
        if "__there__" in said:
            got = _hive("--get", rconv, tmp + "/")                 # lands in tmp/conversation
            if got.returncode != 0 or not os.path.isdir(local):
                return {"skipped": f"could not fetch {rconv} from HIVE, so nothing was changed "
                                   f"there: {_tail(got)}"}
        res = save_into(local, transcript, sid, parts, current=current)
        mk = _hive(f"mkdir -p {q} && chmod 750 {q}")
        if mk.returncode != 0:
            return {"skipped": f"cannot make {rconv} on HIVE: {_tail(mk)}"}
        names = ([sid + ".jsonl"] if res["written"] else []) + [INDEX, READABLE]
        for name in names:
            put = _hive("--put", os.path.join(local, name), rconv + "/")
            if put.returncode != 0:
                return {"skipped": f"could not put {name} on HIVE: {_tail(put)}"}
        _hive(f"chmod 640 {q}/*")
        n = len(glob.glob(os.path.join(local, "*.jsonl")))
    out = {"saved": rconv + "/" + sid + ".jsonl", "session_id": sid, "hive": True,
           "entries": res["entry"].get("entries"),
           "redacted_strings": res["entry"].get("redacted_strings"), "conversations": n}
    if not res["written"]:
        out["kept"] = res["why"]
    return out


# ------------------------------------------------------------------------- the readable file
def _fence(text):
    run = max((len(m) for m in re.findall(r"`+", text)), default=0)
    return "`" * max(3, run + 1)


def _block(text, lang=""):
    f = _fence(text)
    return f"{f}{lang}\n{text.rstrip()}\n{f}"


def _shorten(text, n=RESULT_LINES):
    lines = text.splitlines()
    if len(lines) <= n:
        return text
    head, tail = n - 10, 10
    return "\n".join(lines[:head] + [f"... {len(lines) - n} lines not shown ..."]
                     + lines[-tail:])


def _content_text(c):
    """A message's content -- a string, or a list of blocks -- as text; blocks it does not know
    become a short tag."""
    if isinstance(c, str):
        return c
    out = []
    for b in c if isinstance(c, list) else []:
        if isinstance(b, dict) and b.get("type") == "text":
            out.append(str(b.get("text") or ""))
        elif isinstance(b, dict) and b.get("type") == "image":
            out.append("[image]")
        elif isinstance(b, str):
            out.append(b)
        elif isinstance(b, dict):
            out.append(f"[{b.get('type') or 'block'}]")
    return "\n".join(out)


def _tool_call(b):
    name, inp = b.get("name") or "tool", b.get("input")
    if isinstance(inp, dict) and isinstance(inp.get("command"), str):
        head = f"**Tool: {name}**" + (f" -- {inp['description']}"
                                       if isinstance(inp.get("description"), str) else "")
        return head + "\n\n" + _block(inp["command"], "bash")
    body = json.dumps(inp, ensure_ascii=False, indent=1) if inp is not None else ""
    return f"**Tool: {name}**" + ("\n\n" + _block(_shorten(body, 20), "json") if body else "")


def _ts(e):
    t = e.get("timestamp")
    return t.replace("T", " ")[:19] + " UTC" if isinstance(t, str) and len(t) >= 19 else ""


def _entries(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw in fh:
            if not raw.strip():
                continue
            try:
                yield json.loads(raw)
            except ValueError:
                yield None


def render_one(path):
    """One conversation as Markdown lines. Never raises: an entry it cannot read is one line."""
    L, results, skipped = [], {}, {"metadata": 0, "meta": 0, "thinking": 0}
    entries = list(_entries(path))
    for e in entries:                          # tool results first, to show under each call
        try:
            c = e["message"]["content"]
            for b in c if isinstance(c, list) else []:
                if isinstance(b, dict) and b.get("type") == "tool_result":
                    results[b.get("tool_use_id")] = _content_text(b.get("content"))
        except (KeyError, TypeError):
            continue
    for e in entries:
        try:
            L += _render_entry(e, results, skipped)
        except Exception:                     # the format is internal: never crash on it
            L.append(f"- [unrecognised entry: {_kind(e)}]")
    tail = [f"{v} {k}" for k, v in (("metadata entries", skipped["metadata"]),
                                     ("meta messages (skill text, reminders)", skipped["meta"]),
                                     ("thinking blocks", skipped["thinking"])) if v]
    if tail:
        L += ["", "*Not shown: " + ", ".join(tail) + ".*"]
    return L


def _kind(e):
    return f"type={e.get('type')!r}" if isinstance(e, dict) else "not JSON"


def _render_entry(e, results, skipped):
    if not isinstance(e, dict):
        return [f"- [unrecognised entry: {_kind(e)}]"]
    t = e.get("type")
    if t in _METADATA or str(t).startswith("artifact-"):
        skipped["metadata"] += 1
        return []
    who = " (subagent)" if e.get("isSidechain") else ""
    if t == "summary":
        return ["", f"*[summary] {str(e.get('summary') or '')[:300]}*"]
    if t == "system":
        sub = e.get("subtype") or "system"
        if sub == "compact_boundary":
            return ["", f"### {_ts(e)} -- the conversation was compacted here"]
        if sub == "omitted":                  # between two of this session's segments
            return ["", "---", "", f"*{OMITTED}*", "", "---"]
        if sub in ("summary-omitted", "omitted-entries"):
            return ["", f"*{e.get('content') or SUMMARY_OMITTED}*"]
        text = _content_text(e.get("content")) if e.get("content") is not None else ""
        return ["", f"*[{sub}] {text.strip()[:300]}*" if text.strip() else f"*[{sub}]*"]
    msg = e.get("message")
    if t not in ("user", "assistant") or not isinstance(msg, dict):
        return [f"- [unrecognised entry: {_kind(e)}]"]
    c = msg.get("content")
    if t == "user":
        if e.get("isMeta"):
            skipped["meta"] += 1
            return []
        if e.get("isCompactSummary"):
            text = _content_text(c)
            return ["", f"### {_ts(e)} -- compaction summary{who}", "",
                    _block(_shorten(text, 15))]
        if isinstance(c, list) and c and all(isinstance(b, dict) and b.get("type") == "tool_result"
                                             for b in c):
            return []                          # shown under the call
        text = _content_text(c)
        return ["", f"### {_ts(e)} -- User{who}", "", text.strip() or "*(empty)*"]
    out = []
    blocks = c if isinstance(c, list) else [{"type": "text", "text": c}]
    for b in blocks:
        if not isinstance(b, dict):
            out.append(f"- [unrecognised entry: block {type(b).__name__}]")
        elif b.get("type") == "text" and str(b.get("text") or "").strip():
            out += ["", f"### {_ts(e)} -- Assistant{who}", "", str(b["text"]).strip()]
        elif b.get("type") == "thinking":
            skipped["thinking"] += 1
        elif b.get("type") == "tool_use":
            out += ["", f"#### {_ts(e)} -- tool call{who}", "", _tool_call(b)]
            res = results.get(b.get("id"))
            if res is not None:
                out += ["", "Result:", "", _block(_shorten(res.strip() or "(no output)"))]
        elif b.get("type") != "text":
            out.append(f"- [unrecognised entry: block type={b.get('type')!r}]")
    return out


def render(conv_dir, index):
    L = ["# The conversations of this analysis", "",
         "Claude Code transcripts, **redacted**, in the order they started (save_transcript.py). "
         "**Core-internal**: they hold internal paths and remarks, so they are never delivered "
         "and the session zip leaves them out. This file is a best-effort reading; the "
         "`<session-id>.jsonl` files beside it are the record, and `index.json` lists them. "
         "For the decisions made, read `logs/decisions.md` first.", ""]
    for i, x in enumerate(index, 1):
        L += ["", f"## Conversation {i}: `{x.get('id')}`", "",
              f"{x.get('first') or '?'} to {x.get('last') or '?'} · {x.get('entries', '?')} entries"
              + (f" · Claude Code {x['claude_code_version']}" if x.get("claude_code_version")
                 else "") + f" · saved {x.get('saved_at') or '?'}"]
        path = os.path.join(conv_dir, str(x.get("file") or ""))
        if not os.path.isfile(path):
            L.append("\n*(the file is missing)*")
            continue
        L += render_one(path)
    return "\n".join(L).rstrip() + "\n"


# --------------------------------------------------------------------------------------- CLI
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("session", nargs="?", help="the analysis session folder (on this computer)")
    ap.add_argument("--hive", help="the session folder on HIVE, driven from this computer")
    ap.add_argument("--transcript", help="the transcript file (the hook passes it)")
    ap.add_argument("--session-id", help="the Claude Code session id (default: "
                                         "$CLAUDE_CODE_SESSION_ID)")
    ap.add_argument("--remember", action="store_true",
                    help="only record that this conversation now works on this session (a "
                         "switch, when it was working on another)")
    ap.add_argument("--from-now", action="store_true",
                    help="this session's part of the conversation starts HERE, not at the "
                         "skill's first load: for a conversation about several clients")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--oneline", action="store_true",
                    help="one timestamped line of JSON (the hook's log)")
    a = ap.parse_args(argv)
    global ONELINE
    ONELINE = a.oneline
    if not (a.session or a.hive):
        ap.error("give the session folder, or --hive <HIVE session folder>")
    try:
        if a.remember:
            got, warn = remember(a.session, a.hive, a.session_id, a.transcript,
                                 from_now=a.from_now)
            out = ({"remembered": True, "entry": got} if got else
                   {"skipped": "not a Claude Code session (CLAUDE_CODE_SESSION_ID is not set)"})
            if warn:
                out["warning"] = warn
            emit(out, a.quiet)
            return 0
        save(a.session, a.transcript, a.session_id, a.hive, quiet=a.quiet,
             from_now=a.from_now)
    except Exception as e:                    # never fail a step over this
        emit({"skipped": f"{type(e).__name__}: {e}"}, a.quiet)
    return 0


if __name__ == "__main__":
    sys.exit(main())
