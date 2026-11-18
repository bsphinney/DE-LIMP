#!/usr/bin/env python3
"""
finish_stop_hook.py -- the plugin's Stop hook (hooks/hooks.json): when the agent is about to end
its turn while a skill analysis it ran still has finish items missing (the report, Methods, the
reproducibility bundle, AGENTS.md, README, the saved conversation, the finalize stamp, the
delivery), block the stop ONCE with the checklist, so it does not end silently without them
(2.14.2 field sessions skipped finalize / provenance / report steps and just stopped).

The finish items are session.finish_status() -- the very list `session.py status` prints (one
definition; it already adapts to the session: the delivery only for a CoreOmics submission, the
conversation only for a Claude Code analysis). Nothing here lists items of its own.

A plugin hook runs in EVERY Claude Code session of a user who has the plugin, so it is silent
and instant (exit 0, no output) unless ALL of these hold:
  * stop_hook_active is not true (Claude Code sets it when it is already continuing because a
    Stop hook blocked: never loop);
  * this Claude Code session (session_id) ran a LOCAL analysis session -- a session whose folder is on
    HIVE (every CoreOmics submission, hive_remote) is NOT checked: reading it means an ssh call,
    which a hook with a time budget must not make (a dead VPN would hang the turn's end). hooks.json
    says so, and for those sessions SKILL.md rule 16 puts the check on the agent: `session.py status`
    over hive_exec.sh. Audit 2026-11-17 (lane 4, M3): the hook used to read as a safety net for all
    of them. The sessions it does check are those in the map session.py init / checkpoint.py /
    save_transcript.py record (save_transcript.lookup), as transcript_hook.py reads it;
  * that session has produced results (a search or DE output): a conversation still gathering
    conditions is not "ending without finalizing";
  * no person's hold is recorded (`session.py hold`: finish_hold) -- a stop/pause/abandon is
    respected, never argued with;
  * the checklist has a missing item, and this exact set was not already reminded for this
    conversation (~/.config/ucdavis-proteomics/finish_reminders.json): once per set, so a user
    who decides not to finish is not nagged at every turn.

This is NOT an export path: it writes no deliverable, so it must never fail a user's turn. Any
internal error is one line on stderr and exit 0 -- the opposite of the export paths' "no silent
catch" rule (CLAUDE.md rule 4), on purpose: a reminder that cannot be computed is no reminder.

Stop hook contract (Claude Code hooks reference): stdin JSON {session_id, transcript_path, cwd,
hook_event_name, stop_hook_active}; to make the agent continue, print
{"hookSpecificOutput": {"hookEventName": "Stop", "additionalContext": "..."}} on stdout and exit
0 -- the agent continues with the text, and the user sees "Stop hook feedback" (decision:block
also continues but shows the user "Stop hook error occurred"). Claude Code older than the field
ignores it, so the agent just stops: no reminder, never a failed turn. Python 3.9+, stdlib only.
"""
import hashlib
import json
import os
import sys
import tempfile

TIME_BUDGET_S = 10          # the whole hook; a stalled read (iCloud placeholder, /quobyte) ends here

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)


def _warn(msg):
    try:
        sys.stderr.write(f"[finish_stop_hook] {msg}\n")
    except Exception:
        pass


def _state_path():
    import save_transcript as st
    return os.path.join(st.config_dir(), "finish_reminders.json")


def _read_state():
    try:
        with open(_state_path(), encoding="utf-8") as fh:
            m = json.load(fh)
        return m if isinstance(m, dict) else {}
    except (OSError, ValueError):
        return {}


def _write_state(m):
    path = _state_path()
    os.makedirs(os.path.dirname(path), exist_ok=True)
    keep = dict(list(m.items())[-200:])                 # bounded: the newest conversations
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), prefix=".finish_reminders.", suffix=".tmp")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(keep, fh)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)                            # atomic: a reader sees the old or the new
    except BaseException:
        try:
            os.unlink(tmp)
        except OSError:
            pass
        raise


def reminder(payload):
    """The block reason (str) for this Stop event, or None to let it stop."""
    if payload.get("stop_hook_active"):
        return None
    import save_transcript as st
    import session
    sid = payload.get("session_id") or os.environ.get("CLAUDE_CODE_SESSION_ID")
    lines, sig, outward, inward = [], [], [], 0
    for entry in st.lookup(sid):
        sdir = entry.get("session_dir")
        if entry.get("hive") or not sdir or not os.path.isdir(sdir):
            continue                                    # a HIVE session: NOT checked (see the header)
        if session.finish_hold(sdir) or not session.has_results(sdir):
            continue
        missing = [i for i in session.finish_status(sdir) if not i["done"]]
        if not missing:
            continue
        lines.append(f"Session {sdir}:")
        for i in missing:
            if i["item"] in session.OUTWARD_FINISH_ITEMS:
                # outward-facing (the client delivery): never "do it now" -- ask first (ux M5)
                outward.append(i["label"])
                lines.append(f"  - {i['label']}: {i['detail']} -> ASK the person before doing "
                             "this: it sends results to the client, so never run it unasked"
                             + (f" (when they say so: {i['command']})" if i.get("command")
                                else ""))
            else:
                inward += 1
                lines.append(f"  - {i['label']}: {i['detail']}" +
                             (f" -> {i['command']}" if i.get("command") else ""))
            sig.append(f"{sdir}|{i['item']}|{i['state']}")
    if not lines:
        return None
    digest = hashlib.sha256("\n".join(sorted(sig)).encode()).hexdigest()[:16]
    state = _read_state()
    if state.get(sid) == digest:
        return None                                     # this exact set was already reminded
    state.pop(sid, None)                                # re-insert last: the newest are kept
    state[sid] = digest
    try:
        _write_state(state)
    except Exception as e:
        _warn(f"could not record the reminder ({type(e).__name__}: {e}); it may repeat")
    return ("The analysis you ran has finish items missing (golden rule 16: run "
            "`python3 scripts/session.py status --session <dir>` before telling the user it is "
            "finished).\n" + "\n".join(lines) +
            (f"\nDo {'the others' if outward else 'them'} now, or say to the user which you are leaving out and why (finalize "
             "--override-finish \"<reason>\")." if inward else "") +
            (f"\nAsk the user (Core staff) whether to do this now: {'; '.join(outward)}. Do NOT "
             "run it unasked." if outward else "") +
            " If the user has stopped or abandoned this "
            "analysis, or a person paused or stopped you on the board, run `session.py hold --session "
            "<dir> --why \"<what they said>\"` and stop. "
            "This reminder is shown once.")


class _TimeUp(BaseException):
    pass


def _arm_time_budget():
    """POSIX: SIGALRM after TIME_BUDGET_S ends the hook silently (exit 0). Windows has no
    SIGALRM: no budget there (the hook's own `timeout` in hooks.json still applies)."""
    try:
        import signal

        def _up(signum, frame):
            raise _TimeUp()
        signal.signal(signal.SIGALRM, _up)
        signal.alarm(TIME_BUDGET_S)
    except (AttributeError, ValueError, OSError):
        pass


def main():
    _arm_time_budget()
    try:
        raw = sys.stdin.read() if not sys.stdin.isatty() else ""
        payload = json.loads(raw) if raw.strip() else {}
        if not isinstance(payload, dict):
            return 0
        why = reminder(payload)
    except KeyboardInterrupt:
        raise
    except _TimeUp:
        _warn(f"over the {TIME_BUDGET_S} s budget; no reminder")
        return 0
    except BaseException as e:                          # never fail the user's turn (SystemExit too)
        _warn(f"{type(e).__name__}: {e}")
        return 0
    try:
        import signal
        signal.alarm(0)
    except (AttributeError, ValueError, OSError):
        pass
    if why:
        # additionalContext, not decision:block: the agent still continues (same loop guards),
        # but the user sees "Stop hook feedback", not "Stop hook error occurred" (verified live,
        # Claude Code 2.1.322; hooks reference, "Stop decision control").
        sys.stdout.write(json.dumps({"hookSpecificOutput": {
            "hookEventName": "Stop", "additionalContext": why}}))
        sys.stdout.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
