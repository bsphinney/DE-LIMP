#!/usr/bin/env python3
"""
staff.py -- who may sign a Core decision, and where their name may go. The ONE check behind
every `--by` that records a person's decision (normalization_check.py decide and ack-legacy,
qc_bracket.py ack).

A decision recorded with --by is the PERSON's -- a Core staff member's -- never Claude's. --by is
that person's HIVE login (shaped like one, lowercased: a full name such as "Staff Example" is
refused, list or no list), listed in the Core's staff file:

  CORE_STAFF_FILE (environment), default /quobyte/proteomics-grp/.config/core_staff.txt
  one HIVE login per line; `#` starts a comment; anything after the login on a line is ignored

The file counts only when a Core admin owns it and its folder and nobody else can write either
(slack_collab.trust_problem, the rule for every Core file on HIVE: /quobyte/proteomics-grp is
group-writable and not sticky, so any member could put a list of their own in its place). Such a
list is AUTHORITATIVE: a login on it is accepted -- whatever it looks like (people are called
Claude, Talbot, Kaiser) -- and one not on it is refused. Without a usable list -- a laptop, a Core
that has not written one yet, a file that fails that rule -- nothing can be verified: a login
that is an agent's name (claude, assistant, agent, bot, ... as a whole part of it: "my-bot",
"ai", never "talbot") is refused, and any other login passes with a warning saying it was NOT
VERIFIED.

Names stay with the staff. What reaches a client -- methods, AUDIT.md, tables/,
reproducibility/, the session zip -- says ROLE ("Core staff"). The name and the login go into
staff-only records: files whose names end in STAFF_SUFFIX (".staff.json"), kept in a session's
logs/ or beside the file they are about. core_submission.py deliver never ships one and session.py
finalize leaves them out of the zip, wherever they are.

  python3 staff.py check --by <login>     # {"ok", "by", "error", "warning", "staff_file", "checked"}
"""
import argparse
import datetime
import getpass
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

STAFF_FILE_ENV = "CORE_STAFF_FILE"
STAFF_FILE_DEFAULT = "/quobyte/proteomics-grp/.config/core_staff.txt"
#: What a client deliverable calls whoever decided.
ROLE = "Core staff"
#: A staff-only record's file name ends in this: never delivered, never zipped.
STAFF_SUFFIX = ".staff.json"
STAFF_SCHEMA = "staff_record"
FILE_MODE = 0o640                 # as logs/decisions.md: the Core's group reads it on HIVE

# An agent's name, never a person's -- consulted only when there is no usable staff list. Whole
# parts of the name only (split at anything that is not a letter or digit: a login's . _ -), never
# substrings, so "talbot", "kaiser" or "claudette" are people.
AGENT_WORDS = {"claude", "claudecode", "claudeai", "assistant", "agent", "agents", "bot", "ai",
               "gpt", "chatgpt", "llm", "copilot", "gemini", "openai", "anthropic", "opus",
               "sonnet", "haiku", "skill", "model", "automated", "auto", "default", "system",
               "script", "pipeline"}
_LOGIN = re.compile(r"^[a-z_][a-z0-9_.-]{0,31}$")      # a HIVE username (slack_collab's shape)


def staff_file():
    return os.environ.get(STAFF_FILE_ENV) or STAFF_FILE_DEFAULT


def is_agent_like(name):
    """Is a whole part of `name` an agent's name ("my-bot", "ai", "Claude Code"; not "talbot")?"""
    return any(w in AGENT_WORDS for w in re.split(r"[^a-z0-9]+", (name or "").lower()) if w)


def is_staff_only(name):
    """Is this file a staff-only record (by its name)?"""
    return os.path.basename(name or "").endswith(STAFF_SUFFIX)


def load_staff(path=None):
    """(logins, None) from a usable staff file, else (None, why it cannot be used)."""
    path = path or staff_file()
    import slack_collab as sc                  # the Core's trust rule and its admins
    info = sc._look(path, True)
    if info.get("error") == "FileNotFoundError":
        return None, f"{path} does not exist"
    why = sc.trust_problem(info, sc._folder(path), sc.core_admins())
    if why:
        return None, f"{path} is not trusted: {why}"
    logins = []
    for ln in (info.get("text") or "").splitlines():
        f = ln.split("#", 1)[0].split()
        if f and _LOGIN.match(f[0].lower()) and f[0].lower() not in logins:
            logins.append(f[0].lower())
    if not logins:
        return None, f"{path} lists no HIVE login"
    return logins, None


def check_by(by):
    """May `by` sign a Core decision? -> {"ok", "by", "error", "warning", "staff_file",
    "checked"}: `by` as recorded (the login, lowercased), `checked` how ("on the Core staff list
    (<file>)" or "NOT VERIFIED: the Core staff list is not set up (<why>)")."""
    by = (by or "").strip()
    path = staff_file()
    out = {"ok": False, "by": by, "error": None, "warning": None, "staff_file": path,
           "checked": None}
    if not by:
        out["error"] = ("--by: the HIVE login of the Core staff member whose decision this is "
                        "(never Claude's)")
        return out
    agent = (f"--by {by!r} names an agent, not a person. This decision is a Core staff member's, "
             "never Claude's: ask them, and pass their HIVE login.")
    if not _LOGIN.match(by.lower()):              # a full name, a phrase: never, list or no list
        out["error"] = agent if is_agent_like(by) else (
            f"--by {by!r} is not a HIVE login (lowercase letters, digits, _ . -; what `whoami` "
            "prints on HIVE for them), and the decision is recorded by login: pass the Core staff "
            "member's HIVE login, not their name.")
        return out
    by = out["by"] = by.lower()
    logins, why = load_staff(path)
    if logins is None:
        # nothing to verify against: the agent-name rule is all there is
        if is_agent_like(by):
            out["error"] = agent
            return out
        out.update(ok=True,
                   checked=f"NOT VERIFIED: the Core staff list is not set up ({why})",
                   warning=(f"WARNING: NOT VERIFIED -- the Core staff list is not set up ({why}), "
                            f"so --by {by!r} was checked only for being a login and not an "
                            f"agent's name. A Core admin writes the list ({STAFF_FILE_DEFAULT}: one "
                            "HIVE login per line, owned by a Core admin, chmod 644)."))
        return out
    if by not in logins:                          # the list is authoritative, both ways
        out["error"] = (f"--by {by!r} is not on the Core staff list ({path}). The decision is a "
                        "Core staff member's: pass their HIVE login (what `whoami` prints on HIVE "
                        "for them); a Core admin adds staff to the list.")
        return out
    out.update(ok=True, checked=f"on the Core staff list ({path})")
    return out


def require(by):
    """check_by(), raising ValueError with its error; the warning, if any, goes to stderr."""
    c = check_by(by)
    if not c["ok"]:
        raise ValueError(c["error"])
    if c["warning"]:
        print(c["warning"], file=sys.stderr)
    return c


# The flags that carry a person's name or login on a command line (normalization_check.py decide
# and ack-legacy, qc_bracket.py ack, estimate_params.py / resolve_defaults.py --override-by).
# commands.log holds every command verbatim, and it travels: reproducibility/inputs/ and the
# session zip get it through redact().
NAME_FLAGS = ("--by", "--override-by")
_NAME_ARG = re.compile(r"(?<![\w-])(" + "|".join(re.escape(f) for f in NAME_FLAGS) + r")"
                       r"""(\s+|=)("(?:[^"\\]|\\.)*"|'[^']*'|[^\s"']+)""")


def redact(text):
    """`text` (a command log) with the value of every name flag replaced by the role."""
    return _NAME_ARG.sub(lambda m: f'{m.group(1)}{m.group(2)}"<{ROLE}>"', text)


def login():
    """The login running this command (who TYPED it), or None."""
    try:
        return getpass.getuser()
    except (KeyError, OSError):
        return None


def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def read(path):
    """The entries of a staff-only record ([] when there is none)."""
    try:
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
    except FileNotFoundError:
        return []
    if not (isinstance(d, dict) and d.get("schema") == STAFF_SCHEMA
            and isinstance(d.get("entries"), list)):
        raise ValueError(f"{path} is not a {STAFF_SCHEMA} record")
    return [e for e in d["entries"] if isinstance(e, dict)]


def append(path, entry):
    """Add one entry to the staff-only record at `path` (a STAFF_SUFFIX file). -> the entry."""
    if not is_staff_only(path):
        raise ValueError(f"{path}: a staff-only record's name ends in {STAFF_SUFFIX}")
    entries = read(path) + [entry]
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump({"schema": STAFF_SCHEMA, "schema_version": 1, "staff_only": True,
                   "note": "Core-internal: who decided, by name and login. Never delivered, not "
                           "in the session zip; client documents say \"" + ROLE + "\".",
                   "entries": entries}, fh, indent=2)
        fh.write("\n")
    try:
        os.chmod(tmp, FILE_MODE)
    except OSError:
        pass
    os.replace(tmp, path)
    return entry


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check", help="may this login sign a Core decision?")
    c.add_argument("--by", required=True)
    a = ap.parse_args()
    r = check_by(a.by)
    print(json.dumps(r, indent=2))
    return 0 if r["ok"] else 2


if __name__ == "__main__":
    sys.exit(main())
