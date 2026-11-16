#!/usr/bin/env python3
"""
log_decision.py -- one line of the analysis's decisions log, <session>/logs/decisions.md: what
was decided at a decision point, and why. Any agent can write it (the skill is agent-agnostic;
without Claude Code there is no transcript to save), and on Claude Code it is the curated record
beside the raw conversation (save_transcript.py) -- far quicker for a reviewer to read.

  python3 log_decision.py --session <session dir> --step '3. design' --what-file - --why-file - <<'EOF'
  Contrasts: Old - Young per bait
  ---
  the user confirmed them after seeing the design table
  EOF

The text comes in as DATA: --what-file / --why-file F (a file) or - (standard input; when both
are -, the what, a line `---`, then the why), with a heredoc quoted <<'EOF', in which nothing is
expanded. --what "..." / --why "..." still work for text with none of $ ` \\ " ', but in double
quotes bash runs anything in backticks or $(...) -- the user's own words, an error text -- on
this computer, and again on HIVE through hive_exec.sh's `bash -l -c` (2.12 security review).
A text is read up to 64 KiB.

The engine choice of step 4 is logged from its record instead of typed (--from-json): the
workflow manifest (engine.selection), resolve_defaults.py's printed summary (selection) or a
search_provenance.json (engine_selection). Its client sentence is quoted as written, never
re-worded; the staff record beside it is not read.

  python3 log_decision.py --session <S> --from-json wf/workflow.manifest.json

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


TEXT_MAX = 65536
SEPARATOR = "---"      # between the what and the why when both come from standard input


def _read(src, limit):
    """The text in file `src`, or standard input for "-": UTF-8 (a BOM dropped, a bad byte
    replaced), cut at `limit` bytes with a note. Raises OSError."""
    if src == "-":
        data = sys.stdin.buffer.read(limit + 1)
    else:
        with open(src, "rb") as fh:
            data = fh.read(limit + 1)
    text = data[:limit].decode("utf-8-sig", "replace")
    return text + (" [... cut at 64 KiB]" if len(data) > limit else "")


def texts(a, ap):
    """(what, why) from --what/--why or --what-file/--why-file (see the module text)."""
    what, why = a.what, a.why
    try:
        if a.what_file == "-" and a.why_file == "-":
            lines = _read("-", 2 * TEXT_MAX).splitlines()
            cut = next((i for i, ln in enumerate(lines) if ln.strip() == SEPARATOR), None)
            if cut is None:
                ap.error(f"--what-file - --why-file -: standard input needs the what, a line "
                         f"{SEPARATOR!r}, then the why")
            what, why = "\n".join(lines[:cut]), "\n".join(lines[cut + 1:])
        else:
            if a.what_file:
                what = _read(a.what_file, TEXT_MAX)
            if a.why_file:
                why = _read(a.why_file, TEXT_MAX)
    except OSError as e:
        ap.error(f"could not read the text: {e}")
    for name, v in (("what", what), ("why", why)):
        if not (v or "").strip():
            ap.error(f"the {name} is empty")
    return what, why


def selection_of(doc):
    """The engine-selection record (engine_choice.choose) inside a workflow manifest
    (engine.selection), resolve_defaults.py's summary (selection) or a search record
    (engine_selection); None when the document holds none."""
    if not isinstance(doc, dict):
        return None
    for cand in ((doc.get("engine") or {}).get("selection") if isinstance(doc.get("engine"), dict)
                 else None, doc.get("selection"), doc.get("engine_selection")):
        if isinstance(cand, dict) and cand.get("engine"):
            return cand
    return None


def from_json(src):
    """(what, why) for the engine choice in `src` (a file, or - for standard input), from the
    record's own fields: the client sentence quoted when there is one, never re-worded; the
    staff record is not read. Raises ValueError (no selection in it) or OSError."""
    doc = json.loads(_read(src, 4 * TEXT_MAX))
    sel = selection_of(doc)
    if sel is None:
        raise ValueError("it holds no engine selection (engine.selection, selection or "
                         "engine_selection): write the decision with --what/--why")
    from make_methods import ENGINE_LABEL
    lab = lambda e: ENGINE_LABEL.get(str(e or "").lower(), e or "not recorded")  # noqa: E731
    what = (f"Search engine: {lab(sel['engine'])} ({sel.get('how') or 'how not recorded'})"
            + ("; identifications only" if sel.get("ids_only") else ""))
    why = sel.get("sentence") or (
        f"the default for {sel.get('data_type') or 'this data type'} is "
        f"{lab(sel.get('default'))}; how: {sel.get('how') or 'not recorded'}")
    why += f" (decided by {sel.get('decided_by') or 'not recorded'})"
    if sel.get("replayed"):
        why += f"; {sel['replayed']}"
    return what, why


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--session", required=True, help="the analysis session folder")
    ap.add_argument("--from-json", metavar="F|-",
                    help="log step 4's engine choice from its record (a workflow manifest, "
                         "resolve_defaults.py's summary or a search_provenance.json) instead of "
                         "--what/--why")
    w = ap.add_mutually_exclusive_group()
    w.add_argument("--what", help="what was decided (text with none of $ ` \\ \" ')")
    w.add_argument("--what-file", metavar="F|-", help="what was decided, from a file or - (stdin)")
    y = ap.add_mutually_exclusive_group()
    y.add_argument("--why", help="why: the user's words, the evidence (none of $ ` \\ \" ')")
    y.add_argument("--why-file", metavar="F|-", help="why, from a file or - (stdin; with "
                                                     "--what-file - too: what, a line ---, why)")
    ap.add_argument("--step", help="the SKILL.md step, e.g. '3. design'")
    a = ap.parse_args(argv)
    typed = any(x is not None for x in (a.what, a.what_file, a.why, a.why_file))
    if a.from_json:
        if typed:
            ap.error("--from-json writes the what and the why from the record: leave out "
                     "--what/--why")
        try:
            what, why = from_json(a.from_json)
        except (OSError, ValueError) as e:
            ap.error(f"--from-json {a.from_json}: {e}")
        a.step = a.step or "4. defaults (engine)"
    else:
        for need, have in (("--what --what-file", (a.what, a.what_file)),
                           ("--why --why-file", (a.why, a.why_file))):
            if all(x is None for x in have):
                ap.error(f"one of the arguments {need} is required (or --from-json)")
        what, why = texts(a, ap)
    if not os.path.isdir(a.session):
        print(json.dumps({"logged": False, "reason": f"no session folder {a.session!r}"}))
        return 0
    print(json.dumps({"logged": True, "path": append(a.session, what, why, a.step)}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
