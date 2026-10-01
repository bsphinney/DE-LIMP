#!/usr/bin/env python3
"""
log_decision.py -- one line of the analysis's decisions log, <session>/logs/decisions.md: what
was decided at a decision point, and why. Any agent can write it (the skill is agent-agnostic;
without Claude Code there is no transcript to save), and on Claude Code it is the curated record
beside the raw conversation (save_transcript.py) -- far quicker for a reviewer to read.

  python3 log_decision.py --session <session dir> --what "Contrasts: Old - Young per bait" \\
      --why "the user confirmed them after seeing the design table" [--step "3. design"]

Log it where SKILL.md says: the conditions confirmed, the contrasts chosen, the defaults
confirmed or changed, any override (--unchecked, --force, a gate accepted), anything the user
declined or asked for that the skill does not do by default. Redacted as the saved
conversation is (save_transcript.Redactor: the one redaction for the analysis record). Like the
conversation it is CORE-INTERNAL: never delivered, not in the session zip, readable by the
owner and group only (0640). Appends; never rewrites what is there.
"""
import argparse
import datetime
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

HEADER = ("# Decisions log\n\n"
          "What was decided during this analysis, and why -- written by the agent at each "
          "decision point (log_decision.py). The Claude Code conversation, when there is one, is "
          "in `logs/conversation/` (Core-internal); every command run is in `logs/commands.log`.\n")


_REDACTOR = None


def one(text):
    global _REDACTOR
    if _REDACTOR is None:
        from save_transcript import Redactor
        _REDACTOR = Redactor()
    return " ".join(_REDACTOR.text(str(text)).split())


def append(session, what, why, step=None, when=None):
    path = os.path.join(os.path.abspath(session), "logs", "decisions.md")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    when = when or datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = ["", f"## {when} -- {one(what)}", f"- **Why:** {one(why)}"]
    if step:
        lines.append(f"- **Step:** {one(step)}")
    if os.environ.get("CLAUDE_CODE_SESSION_ID"):
        lines.append(f"- **Conversation:** `{os.environ['CLAUDE_CODE_SESSION_ID']}`")
    new = not os.path.isfile(path) or os.path.getsize(path) == 0
    with open(path, "a", encoding="utf-8") as fh:
        fh.write((HEADER if new else "") + "\n".join(lines) + "\n")
    try:
        os.chmod(path, 0o640)
    except OSError:
        pass
    return path


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--session", required=True, help="the analysis session folder")
    ap.add_argument("--what", required=True, help="what was decided")
    ap.add_argument("--why", required=True, help="why: the user's words, the evidence")
    ap.add_argument("--step", help="the SKILL.md step, e.g. '3. design'")
    a = ap.parse_args(argv)
    if not os.path.isdir(a.session):
        print(json.dumps({"logged": False, "reason": f"no session folder {a.session!r}"}))
        return 0
    print(json.dumps({"logged": True, "path": append(a.session, a.what, a.why, a.step)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
