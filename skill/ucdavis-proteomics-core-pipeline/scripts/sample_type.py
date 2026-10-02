#!/usr/bin/env python3
"""
sample_type.py -- what the samples are (tissue, cell line, biofluid, IP ...), AS THE USER SAID
IT: the session's one record of it, <session>/input/sample_type.json.

Why (a staff report, 2.10): an abbreviation in a cohort's file names was read as a tissue, and
the agent wrote that tissue through the draft report and Methods. The expert reviewer found the
data did not support it. SKILL.md asked what the samples are only to decide whether they are
keratin, so nothing else
held the answer and the report took one from the file names. The sample type is now asked for
every analysis (SKILL.md step 3), recorded here in the user's own words, and read by the
analysis brief (analysis_prompt.py) -- which says NOT STATED, and forbids guessing, when it
was not recorded. Nothing in the skill infers it.

  python3 sample_type.py set  --session <session> --stated "mouse leg muscle, flash-frozen"
  python3 sample_type.py set  --session <session> --not-stated     # the user does not know
  python3 sample_type.py show --session <session>                  # JSON; exit 0

The record is written only from what the user said; `--stated` refuses an empty answer.
Stdlib only.
"""
import argparse
import datetime
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

FILE = "sample_type.json"
SCHEMA = "sample_type/1"
# The words every reader uses when nothing was recorded: never a guess in its place (rule 2).
NOT_STATED = "NOT STATED by the user"
NEVER_INFER = ("never infer it from file or folder names, sample IDs or abbreviations (a "
               "file-name abbreviation was once written up as a tissue the data did not "
               "support)")


def path_for(session):
    return os.path.join(os.path.abspath(os.path.expanduser(session)), "input", FILE)


def record(session, stated=None, user_did_not_say=False, when=None):
    """Write the session's sample-type record -> its path. `stated`: the user's words, as given.
    `user_did_not_say`: asked, and the user does not know or did not answer."""
    if bool(stated and stated.strip()) == bool(user_did_not_say):
        raise ValueError("give the user's answer (--stated) or say they did not give one "
                         "(--not-stated), not both or neither")
    p = path_for(session)
    if not os.path.isdir(os.path.dirname(p)):
        raise ValueError(f"not a session directory (no input/): {os.path.dirname(os.path.dirname(p))}")
    rec = {"schema": SCHEMA,
           "stated": " ".join(stated.split()) if stated else None,
           "source": "user" if stated else "asked; the user did not say",
           "recorded_at": when or datetime.datetime.now(datetime.timezone.utc).strftime(
               "%Y-%m-%dT%H:%M:%SZ"),
           "rule": "as the user stated it; " + NEVER_INFER}
    tmp = p + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, p)
    return p


def load(session):
    """The record, or None when the session has none (never asked, or a session from before
    2.10). A record that exists but cannot be read raises: it is not the same as no answer."""
    if not session:
        return None
    p = path_for(session)
    if not os.path.isfile(p):
        return None
    with open(p, encoding="utf-8") as fh:
        rec = json.load(fh)
    if not isinstance(rec, dict) or rec.get("schema") != SCHEMA:
        raise ValueError(f"{p} is not a {SCHEMA} record")
    return rec


def describe(rec):
    """One line for a reader: the user's words, or NOT STATED (and why)."""
    if rec and rec.get("stated"):
        return f"{rec['stated']} (as the user stated it)"
    if rec:
        return f"{NOT_STATED} (asked; the user did not say)"
    return f"{NOT_STATED} (not recorded in this session -- ask the user)"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("set", help="record the user's answer")
    s.add_argument("--session", required=True)
    g = s.add_mutually_exclusive_group(required=True)
    g.add_argument("--stated", help="what the user said the samples are, in their words")
    g.add_argument("--not-stated", action="store_true",
                   help="asked, and the user does not know or did not say")
    sh = sub.add_parser("show", help="print the record (JSON)")
    sh.add_argument("--session", required=True)
    a = ap.parse_args()
    try:
        if a.cmd == "set":
            p = record(a.session, a.stated, a.not_stated)
            print(json.dumps({"written": p, "sample_type": describe(load(a.session))}, indent=2))
        else:
            rec = load(a.session)
            print(json.dumps({"record": rec, "sample_type": describe(rec)}, indent=2))
    except (ValueError, OSError) as e:
        sys.exit(f"[sample_type] {e}")


if __name__ == "__main__":
    main()
