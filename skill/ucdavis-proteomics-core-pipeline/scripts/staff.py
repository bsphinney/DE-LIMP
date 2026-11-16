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
reproducibility/, the session zip -- says ROLE ("Core staff"); an override (estimate_params.py,
resolve_defaults.py --override-by) says it only for a login on the list, else USER_ROLE ("the
user"), through role_for(). The name and the login go into
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
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

STAFF_FILE_ENV = "CORE_STAFF_FILE"
STAFF_FILE_DEFAULT = "/quobyte/proteomics-grp/.config/core_staff.txt"
#: What a client deliverable calls whoever decided.
ROLE = "Core staff"
#: What it calls whoever set a value when they are not known to be Core staff (role_for).
USER_ROLE = "the user"
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


#: THE rule for files that are Core-internal by NAME, wherever they are (C2 security M1): never
#: delivered (core_submission.py, every mode) and never zipped (session.py finalize). They are
#:   * the scheduler's job logs -- node names, paths with the HIVE login, the job-end hook's FRAN,
#:     Slack and run-log lines. Every one this skill writes is <name>_<jobid>.log, or
#:     <name>_<arrayjob>_<task>.log for an array (#SBATCH --output ..._%j.log / _%A_%a.log /
#:     %x_%j.log), plus SLURM's own default slurm-<id>.out/.err;
#:   * the records of the Core's own steps: fran_deposit.json (fran_deposit.RECEIPT: the login
#:     that staged the search into the Core's corpus, and its entry) and node_faults.json
#:     (node_fault.RECORD: the nodes a job failed on). A test ties both names to their writers.
#: Job SCRIPTS (.sh/.sbatch) are not in it: the zip keeps a session's scripts/ as its record.
CORE_INTERNAL_FILES = re.compile(r"^(?:slurm-.*\.(?:out|err)|.+_\d+(?:_\d+)?\.log|"
                                 r"fran_deposit\.json|node_faults\.json)$", re.I)
CORE_INTERNAL_REASON = ("Core-internal: a scheduler job log or a record of the Core's own steps "
                        "(nodes, logins, paths), never delivered")


def is_core_internal(name):
    """Is this file Core-internal by its name (CORE_INTERNAL_FILES)?"""
    return bool(CORE_INTERNAL_FILES.match(os.path.basename(name or "")))


#: Client-facing records that can carry the LOGIN of who accepted something (sage_lfq_check.accept
#: records the login that ran --accept-lfq-window, in the search's quantification record --
#: either of its names, sage_lfq_check.RECORD_NAMES (a test ties them) -- and again in
#: search_provenance.json `sage_lfq_check`): {file name: the key path to the `accepted` block}.
#: The session keeps the login (the Core's own record, as commands.log keeps --by); a delivery and
#: the zip carry client_copy() instead (C2 security M1).
ACCEPTED_BY_RECORDS = {"quantification_check.json": ("accepted",),
                       "sage_lfq_check.json": ("accepted",),
                       "search_provenance.json": ("sage_lfq_check", "accepted")}
ACCEPTED_BY_SOURCE = ("the role of whoever accepted (staff.role_for); the login stays in the "
                      "Core's own copy")


def client_copy(path):
    """The text a client gets of the record at `path` when it names a login as who accepted
    (ACCEPTED_BY_RECORDS: `by` -> role_for(by)); None when it names none, is not such a record, or
    is not a JSON object -- then it is copied as it is."""
    keys = ACCEPTED_BY_RECORDS.get(os.path.basename(path or ""))
    if not keys:
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            rec = json.load(fh)
    except (OSError, ValueError):
        return None
    block = rec
    for k in keys:
        block = block.get(k) if isinstance(block, dict) else None
    by = block.get("by") if isinstance(block, dict) else None
    if not by or by in (ROLE, USER_ROLE):
        return None
    block["by"] = role_for(by)
    block["by_source"] = ACCEPTED_BY_SOURCE
    return json.dumps(rec, indent=2) + "\n"


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


# On a laptop, HIVE is asked about the account this computer's HIVE login reaches (hive_exec.sh,
# hive.env): staff.py's own check of `id -un` there -- what the ssh key proves, never a name this
# computer claims.
STAFF_CHECK_ON_HIVE = 'python3 ~/proteomics-pipeline/scripts/staff.py check --by "$(id -un)"'


def _on_hive():
    """core_submission.on_hive, the skill's test (its CORE_ON_HIVE switch is for tests)."""
    try:
        from core_submission import on_hive
    except ImportError:
        return os.path.isdir("/quobyte")
    return on_hive()


def _uid_login():
    """This account's login, by its uid: never $USER or $LOGNAME, which anyone can set."""
    try:
        import pwd
        return pwd.getpwuid(os.getuid()).pw_name
    except (ImportError, KeyError, AttributeError):
        return None


def _core_staff_on_hive():
    """From a laptop: ask HIVE whether the account this computer's HIVE login reaches is on the
    Core staff list. -> (login, how, None) or (None, None, why)."""
    try:
        # stdin closed: ssh would read the caller's (a heredoc meant for the caller, a pipe)
        r = subprocess.run(["bash", os.path.join(HERE, "hive_exec.sh"), STAFF_CHECK_ON_HIVE],
                           capture_output=True, text=True, timeout=120,
                           stdin=subprocess.DEVNULL)
        got = json.loads(r.stdout[r.stdout.find("{"):])      # past anything a login shell prints
    except (OSError, subprocess.SubprocessError, ValueError) as e:
        return None, None, ("the Core staff list is on HIVE and could not be asked from here "
                            f"({type(e).__name__}): this computer needs its HIVE login "
                            "(check_access.sh saves it) and the skill on HIVE "
                            "(hive_exec.sh --put-skill)")
    if not isinstance(got, dict):
        return None, None, "HIVE's answer about the Core staff list was not readable"
    checked = str(got.get("checked") or "")
    if got.get("ok") and checked.startswith("on the Core staff list"):
        return got.get("by"), checked + ", asked on HIVE", None
    return None, None, str(got.get("error") or checked or "HIVE did not confirm it")


def core_staff():
    """Is the person running this Core staff? -> (login, how it was checked, None) or (None, None,
    why not). THE one test of "Core staff" (2.15 review M1): the podcast's standing consent
    (make_podcast.py --core-run), the PRM follow-up (skyline_peaks.py), marking a session as Core
    client work (session.py core-client) and the line that says the Core's survey is off
    (submission_report.survey_off_note). On HIVE this account, by uid, must be on the list
    load_staff trusts (STAFF_FILE_ENV, owned by a Core admin) -- the list check_by and role_for
    read; on a laptop HIVE is asked about the account its key logs in as (one hive_exec.sh call).
    Reading /quobyte/proteomics-grp is not enough: not every member of the group is Core staff.
    Fails closed: no usable list is no."""
    if not _on_hive():
        return _core_staff_on_hive()
    login = _uid_login()
    logins, why = load_staff()
    if logins is None:
        return None, None, f"the Core staff list cannot be used ({why})"
    if not login or login not in logins:
        return None, None, (f"{login or 'this account'} is not on the Core staff list "
                            f"({staff_file()}); a Core admin adds staff to it")
    return login, f"on the Core staff list ({staff_file()})", None


def role_for(login):
    """The role a client document names for `login`: ROLE ("Core staff") only when check_by() finds
    it on a usable Core staff list; else USER_ROLE ("the user") -- a login not on the list, a name
    rather than a login, an agent's name, or anyone at all when there is no usable list, since
    then nothing can be verified. Records nothing and prints nothing: an override is the user's
    choice whoever typed it (estimate_params.override_source), and only its label depends on this."""
    if not (login or "").strip():
        return USER_ROLE
    c = check_by(login)
    return ROLE if c["ok"] and str(c.get("checked") or "").startswith(
        "on the Core staff list") else USER_ROLE


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
#: A flag's value as a command log holds it, one shell word: '\''...'\'' (inside a single-quoted
#: command), \"...\" (inside a double-quoted one, as given to hive_exec.sh), or a run of "...",
#: '...', $'...', backslash-escaped characters and plain characters (For\ X, 'a'"b").
_ARG_VALUE = (r"""('\\''(?:(?!'\\'').)*'\\''|\\"(?:(?!\\").)*\\"|"""
              r"""(?:"(?:[^"\\]|\\.)*"|'[^'\n]*'|\$'(?:[^'\\\n]|\\.)*'|\\.|[^\s"'\\|;&<>()])+)""")
_NAME_ARG = re.compile(r"(?<![\w-])(" + "|".join(re.escape(f) for f in NAME_FLAGS) + r")"
                       r"(\s+|=)" + _ARG_VALUE)


#: What replaces the value after --override-by: no role at all. Anyone may set an override, and a
#: setter not on the staff list is "the user" in the sidecar, so "<Core staff>" would mislead.
REDACTED_SETTER = "<redacted>"


#: The reasons staff give session.py core-client (who the work is for, why there is no submission
#: yet) never travel (2.15 reviews H3, NH1; ruling 97). The decisions log and the staff-only
#: record keep the reason; every copy of commands.log that travels -- the bundle's inputs/, the
#: session zip, the run registry -- goes through redact_log / redact_report, which
#:   1. removes by VALUE every reason recorded in logs/core_client.staff.json (its `why_text` and
#:      `why`), in any quoting or escaping, standing as words of its own: on a core-client command
#:      any piece of it; on any other line only a piece of 2+ words and 12+ characters, and each
#:      such line is a note (ruling 98: an unrelated path or title is never changed silently);
#:   2. on each core-client command (its `\`-continued lines joined), replaces the --why /
#:      --why-file value (REDACTED_REASON_ARG: one shell word, so the quoting around it holds) and
#:      the body of a heredoc fed to it (one REDACTED_REASON line) -- the body only when its end
#:      line is found: the lines after an opener with no end are kept as logged;
#:   3. withholds a core-client command whole (WITHHELD_COMMAND) when its reason cannot be cut out
#:      exactly -- piped or here-string input, a quoted value that does not end on its own line
#:      (never joined to the lines after it), a heredoc with no end, a value left after (2) -- and
#:      says so (the notes, a MANIFEST [NOTE]). Never silent, and no line it cannot account for
#:      is dropped.
#: Known limit (references/core-submissions.md): lines kept after a withheld command are cleaned
#: by value only, so a reason of one word or under 12 characters may stay in a log the redactor
#: cannot parse -- give the reason with the documented heredoc, and core-client warns when it is
#: that short.
REASON_FLAGS = ("--why", "--why-file")
REASON_COMMANDS = (re.compile(r"session\.py[\"']?\s+core-client(?![\w-])"),)
REDACTED_REASON = "<reason redacted: logs/decisions.md has it>"
REDACTED_REASON_ARG = "[reason-redacted]"
WITHHELD_COMMAND = "[core-client mark: reason withheld]"
CORE_CLIENT_RECORD = "core_client.staff.json"        # in logs/, beside commands.log
_REASON_ARG = re.compile(r"(?<![\w-])(" + "|".join(re.escape(f) for f in REASON_FLAGS) + r")"
                         r"(\s+|=)" + _ARG_VALUE)
_HEREDOC = re.compile(r"""(?<!<)<<(?!<)(-?)\s*(?:'([^'\n]+)'|"([^"\n]+)"|\\?([A-Za-z_][\w.-]*))""")


def reasons_for(*places):
    """The reasons staff gave session.py core-client, from the staff-only record(s) in `places`
    (a logs/ folder, a commands.log beside it, or a session folder) -> ([reason text], problem):
    each entry's `why_text` (the reason as given) and `why` (its one-line form). problem: why a
    record that is there could not be read (None otherwise): the caller says so."""
    out, problem, seen = [], None, set()
    for place in places:
        if not place:
            continue
        place = os.path.abspath(os.path.expanduser(place))
        cands = ([place] if place.endswith(CORE_CLIENT_RECORD) else
                 [os.path.join(os.path.dirname(place), CORE_CLIENT_RECORD),
                  os.path.join(place, CORE_CLIENT_RECORD),
                  os.path.join(place, "logs", CORE_CLIENT_RECORD)])
        for c in cands:
            if c in seen or not os.path.isfile(c):
                continue
            seen.add(c)
            try:
                entries = read(c)
            except (OSError, ValueError) as e:
                problem = f"{c}: {type(e).__name__}: {e}"
                continue
            for e in entries:
                if str(e.get("event") or "").startswith("core-client"):
                    for k in ("why_text", "why"):
                        v = e.get(k)
                        if isinstance(v, str) and v.strip() and v not in out:
                            out.append(v)
    return out, problem


def _char_pattern(ch):
    """`ch` as a command log may hold it: plain, behind backslashes, or a ' in its quoted forms."""
    if ch == "'":
        return r"""(?:'\\''|'"'"'|\\*')"""
    if ch == "\\":
        return r"\\+"
    return r"\\*" + re.escape(ch)


#: between two words of a reason: spaces, a backslash-escaped space, a line break (\n, \ at a
#: line end), or the quotes that close one quoted piece and open the next
_WORD_SEP = r"""(?:\\n|\\t|\\?\s|["'])+"""


#: a by-value hit must stand on its own: not inside a longer word, path or file name
_BOUND_BEFORE, _BOUND_AFTER = r"(?<![\w./-])", r"(?![\w./-])"
#: outside a core-client command, a reason is replaced only when it cannot be an ordinary phrase
#: or a path part (ruling 98): 2+ words AND 12+ characters, and every such hit is a [NOTE]
OUTSIDE_MIN_WORDS, OUTSIDE_MIN_CHARS = 2, 12


def _value_patterns(reasons, outside=False):
    """Compiled patterns, longest first, for each whole reason and each line of it, with word
    boundaries. `outside` (a line that is not a core-client command): only pieces of 2+ words and
    12+ characters, and words joined by spaces only (never across quotes or lines), so an
    unrelated path or title is never touched; on a core-client command any piece counts."""
    pieces = []
    for r in reasons:
        for piece in [r] + r.splitlines():
            p, words = piece.strip(), piece.split()
            if not words or p in pieces:
                continue
            if outside and not (len(words) >= OUTSIDE_MIN_WORDS and len(p) >= OUTSIDE_MIN_CHARS):
                continue
            pieces.append(p)
    pieces.sort(key=len, reverse=True)
    sep = r"(?:\\?[ \t])+" if outside else _WORD_SEP
    return [re.compile(_BOUND_BEFORE + sep.join("".join(_char_pattern(c) for c in w)
                                                 for w in p.split()) + _BOUND_AFTER)
            for p in pieces]


def _continued(line):
    r"""Does `line` end with a `\` that continues the command on the next line?"""
    m = re.search(r"(\\+)$", line)
    return bool(m) and len(m.group(1)) % 2 == 1


def _tokens(cmd):
    import shlex
    lx = shlex.shlex(cmd, posix=True, punctuation_chars=True)
    lx.whitespace_split = True
    lx.commenters = ""
    return list(lx)


def _is_reason_command(text):
    return any(c.search(text) for c in REASON_COMMANDS)


def _reason_left(cmd, depth=0):
    """Could a reason still be in `cmd` (a core-client command after the substitutions)? Yes when
    it does not parse, when input is piped or here-stringed into core-client -- or into the
    command that carries it as one quoted word (hive_exec.sh "...", bash -c '...'), which is looked
    into too -- or when a --why / --why-file value is anything but standard input or the marker."""
    try:
        toks = _tokens(cmd)
    except ValueError:
        return True
    ends = (";", "&&", "||", "&")
    anchors = []
    for k, t in enumerate(toks):
        if t == "core-client" and k and re.search(r"session\.py$", toks[k - 1]):
            anchors.append((k, False))
        elif re.search(r"\s", t) and _is_reason_command(t):
            if depth >= 3 or _reason_left(t, depth + 1):
                return True
            anchors.append((k, True))
    for k, nested in anchors:
        start = max([i + 1 for i, x in enumerate(toks[:k]) if x in ends] or [0])
        if any(x in ("|", "|&") for x in toks[start:k]):
            return True                               # the reason comes down a pipe
        rest = toks[k + 1:]
        for m, x in enumerate(rest):
            if x in ends or x in ("|", "|&"):
                break
            if x == "<<<":
                return True                           # a here-string is the reason
            if nested:
                continue
            if x in REASON_FLAGS:
                v = rest[m + 1] if m + 1 < len(rest) else ""
                if v not in (REDACTED_REASON_ARG,) + (("-",) if x == "--why-file" else ()):
                    return True
            elif x.startswith(tuple(f + "=" for f in REASON_FLAGS)):
                f, v = x.split("=", 1)
                if v not in (REDACTED_REASON_ARG,) + (("-",) if f == "--why-file" else ()):
                    return True
    return False


def _parses(text):
    try:
        _tokens(text)
        return True
    except ValueError:
        return False


_MAX_BODY_LINES = 30                  # a reason is a few lines: an "end" further down is not its


def _heredoc_end(lines, j, delim):
    """The index of the line that ends a heredoc opened just before lines[j], or None: a line that
    is the delimiter, within _MAX_BODY_LINES, before any line that opens another heredoc (whose own
    end could be taken for this one's)."""
    for k in range(j, min(len(lines), j + _MAX_BODY_LINES + 1)):
        if lines[k].strip() == delim:
            return k
        if _HEREDOC.search(lines[k]):
            return None
    return None


def _redact_reasons(text, reasons=()):
    """(text, notes): `text` with the core-client reasons taken out (see REASON_FLAGS above);
    notes say what was withheld whole, and where a reason was found outside a core-client
    command (ruling 98: never an unrelated line changed silently)."""
    on_cc, outside = _value_patterns(reasons), _value_patterns(reasons, outside=True)
    lines, out, notes, i, hits = text.split("\n"), [], [], 0, 0

    def scrub(line):                      # a line that is no core-client command
        nonlocal hits
        new = line
        for pat in outside:
            new = pat.sub(REDACTED_REASON_ARG, new)
        hits += new != line
        return new

    def withhold(why):
        out.append(WITHHELD_COMMAND)
        notes.append(why)

    def sub_arg(text):
        return _REASON_ARG.sub(lambda m: m.group(0) if m.group(3) == "-" else    # standard input
                               f"{m.group(1)}{m.group(2)}{REDACTED_REASON_ARG}", text)

    while i < len(lines):
        j, cmd = i + 1, lines[i]
        while _continued(cmd) and j < len(lines):         # a `\`-continued command, joined
            cmd, j = cmd[:-1].rstrip() + " " + lines[j].strip(), j + 1
        if not _is_reason_command(cmd):
            out += [scrub(x) for x in lines[i:j]]
            i = j
            continue
        for pat in on_cc:
            cmd = pat.sub(REDACTED_REASON_ARG, cmd)
        if not _parses(cmd):
            # a quote that does not close on its own line is never joined to the lines after
            # it (review2 R1: a later line's quote would close it and swallow that line)
            withhold("a core-client command's quoted value does not end on its line: the command "
                     "was withheld, the lines after it kept as logged")
            i = j
            continue
        body_end, docs = None, list(_HEREDOC.finditer(cmd))
        if len(docs) > 1:
            withhold("a core-client command with more than one heredoc was withheld")
            i = j
            continue
        if docs:
            d = docs[0]
            delim = d.group(2) or d.group(3) or d.group(4)
            body_end = _heredoc_end(lines, j, delim)
            if body_end is None:
                withhold("a core-client command's heredoc has no end line in the log: the "
                         "command was withheld, the lines after it kept as logged")
                i = j
                continue
            cmd_wo = cmd[:d.start()] + cmd[d.end():]
        else:
            cmd_wo = cmd
        cmd, cmd_wo = sub_arg(cmd), sub_arg(cmd_wo)
        if _reason_left(cmd_wo):
            withhold("a core-client command whose reason could not be cut out exactly (piped or "
                     "here-string input, or a quoting form not recognised) was withheld")
            i = (body_end + 1) if body_end is not None else j
            continue
        out.append(cmd)
        if body_end is not None:
            if body_end > j:
                out.append(REDACTED_REASON)
            out.append(lines[body_end])
            i = body_end + 1
        else:
            i = j
    if hits:
        notes.append(f"a reason recorded for a core-client mark was found in {hits} other "
                     f"line(s) and replaced there by {REDACTED_REASON_ARG}")
    return "\n".join(out), notes


def redact_report(text, reasons=()):
    """(text, notes): `text` (a command log) as it may travel -- the value of every name flag
    replaced (--by, staff-only by construction, by the role; --override-by, which anyone may set,
    by REDACTED_SETTER) and the core-client reasons taken out (REASON_FLAGS above; `reasons` from
    reasons_for). notes: one line per command withheld whole, for the writer's MANIFEST [NOTE]."""
    def sub(m):
        word = REDACTED_SETTER if m.group(1) == "--override-by" else f"<{ROLE}>"
        return f'{m.group(1)}{m.group(2)}"{word}"'
    text, notes = _redact_reasons(text, reasons)
    return _NAME_ARG.sub(sub, text), notes


def redact(text, reasons=()):
    """redact_report's text alone. THE one redactor of commands.log, with redact_log: provenance.py
    (the bundle), session.py finalize (the zip) and record_run.py (the run registry)."""
    return redact_report(text, reasons)[0]


def redact_log(path, *places):
    """(text, notes) of the command log at `path` as it may travel: redact_report with the reasons
    recorded beside it (logs/core_client.staff.json) and in `places` (a session folder). A staff
    record that cannot be read is a note: the reasons were then cut out by their form alone."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    reasons, problem = reasons_for(path, *places)
    text, notes = redact_report(text, reasons)
    if problem:
        notes = notes + [f"the core-client staff record could not be read ({problem}): reasons "
                         "were taken out by their command's form only"]
    return text, notes


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
