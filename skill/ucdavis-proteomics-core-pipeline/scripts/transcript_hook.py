#!/usr/bin/env python3
"""
transcript_hook.py -- the PreCompact and SessionEnd hook (the plugin's hooks/hooks.json): save
this conversation into the analysis session it ran, before a compaction and when it ends.

A plugin hook runs in EVERY Claude Code session of a user who has the plugin, so it is a NO-OP
outside a skill analysis: it reads the hook's JSON on stdin, looks its session_id up in the map
session.py init / checkpoint.py status / save_transcript.py record
(~/.config/ucdavis-proteomics/transcript_sessions.json), and exits 0 when it is not there.
(Hooks in the SKILL.md frontmatter would be scoped to sessions that used the skill, but in a
live `claude -p` run, 2026-09-28, Claude Code 2.1.284, they registered and SessionEnd never
ran; the plugin's hooks.json ran for SessionEnd and PreCompact both.)

It never blocks and never fails the user's session: it starts save_transcript.py DETACHED (its
own session, stdin closed) and exits 0 at once -- SessionEnd hooks share a 1.5-second budget,
and a PreCompact hook that exits 2 would block the compaction. Whatever goes wrong is one line
in ~/.config/ucdavis-proteomics/transcript_hook.log, never an exit status.
"""
import datetime
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
LOG_MAX = 1 << 20                        # the log is cut back to its last half past 1 MB


def log(line):
    try:
        import save_transcript as st
        path = st.hook_log_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        if os.path.isfile(path) and os.path.getsize(path) > LOG_MAX:
            with open(path, "rb") as fh:
                fh.seek(-LOG_MAX // 2, os.SEEK_END)
                keep = fh.read()
            with open(path, "wb") as fh:
                fh.write(keep)
        stamp = datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()
        with open(path, "a", encoding="utf-8") as fh:
            fh.write(f"{stamp} {line}\n")
    except Exception:
        pass


def main():
    try:
        raw = sys.stdin.read() if not sys.stdin.isatty() else ""
        payload = json.loads(raw) if raw.strip() else {}
        if not isinstance(payload, dict):
            payload = {}
    except Exception as e:
        log(f"[error] unreadable hook input: {type(e).__name__}: {e}")
        return 0
    try:
        import save_transcript as st
        sid = payload.get("session_id") or os.environ.get("CLAUDE_CODE_SESSION_ID")
        sessions = st.lookup(sid)
        if not sessions:
            return 0                          # not a skill analysis: nothing to do
        event = payload.get("hook_event_name") or "?"
        why = payload.get("reason") or payload.get("trigger") or ""
        path = st.hook_log_path()
        os.makedirs(os.path.dirname(path), exist_ok=True)
        for entry in sessions:                # one conversation may run several analyses
            args = [sys.executable, os.path.join(HERE, "save_transcript.py"), "--oneline",
                    "--session-id", sid]
            if payload.get("transcript_path"):
                args += ["--transcript", str(payload["transcript_path"])]
            args += (["--hive", entry["hive"]] if entry.get("hive")
                     else [entry.get("session_dir") or ""])
            with open(path, "a", encoding="utf-8") as out:
                kw = {"stdin": subprocess.DEVNULL, "stdout": out, "stderr": out,
                      "close_fds": True}
                if os.name == "nt":
                    kw["creationflags"] = 0x00000008 | 0x00000200    # DETACHED, NEW_PROCESS_GROUP
                else:
                    kw["start_new_session"] = True
                subprocess.Popen(args, **kw)
            log(f"{event} {why} {sid}: saving to {entry.get('hive') or entry.get('session_dir')}")
    except Exception as e:
        log(f"[error] {type(e).__name__}: {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
