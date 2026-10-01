#!/usr/bin/env python3
"""
notes.py  --  Notes from the Core: Brett (or another trusted sender) leaves a note for a staff
member, and that person's Claude shows it when they next run the skill, and writes a reply back.

WHY THIS EXISTS
---------------
report_issue.sh is the one channel from staff to the Core (skill_issues/). Nothing went the
other way: on 2026-09-29 Brett wanted Michelle's Claude to read what he had to say about her
runs while she was in the middle of one, without her having to notice an email or a Slack
message first. This is that channel in reverse, beside skill_issues/ and skill_runs/:

  /quobyte/proteomics-grp/skill_notes/        (SKILL_NOTES_DIR overrides it -- for the tests)
    SENDERS                                   who else may send (TRUST, below)
    msalemi/                                  one folder per recipient HIVE user; all/ = everyone
      20260929T231000Z_brettsp_keratin.md                               a note
      20260929T231000Z_brettsp_keratin.read.msalemi                     her read receipt (JSON)
      20260929T231000Z_brettsp_keratin.reply.msalemi.20260930T160500Z.md  her reply

A note is header lines -- From:, To:, Date:, Subject:, and Session: when it concerns one session
folder -- then a blank line, then Markdown.

TRUST: A NOTE IS PROMPT TEXT THAT GOES INTO SOMEONE ELSE'S CLAUDE
-----------------------------------------------------------------
* The sender is the file's OWNER on disk, never its From: line, which anyone can type. A note
  whose From: disagrees with its owner is still shown when the owner is trusted, and flagged.
* The trusted senders are the Core admins, plus the names listed in skill_notes/SENDERS. The
  admins are scripts/core_admins.txt, shipped with the skill, not anything found in the shared
  folder: /quobyte/proteomics-grp is writable by the whole group and not sticky, so any member
  could rename skill_notes/ and make a new one of their own. The folder counts only when it is
  the folder itself (not a link), sticky, and a Core admin's; otherwise nothing in it is shown,
  and `check` says why. SENDERS counts only when it is one plain file a Core admin owns that
  nobody else can write (not group- or other-writable, not a link, not hard-linked). A missing,
  unreadable or empty core_admins.txt means no admins, so nothing is trusted. Changing it is a
  code change, shipped with the skill.
* A note from anyone else is never shown. `check` reports only how many there are, and whose.
* Nor is a note that someone other than its owner could have changed: one that is group- or
  other-writable, or hard-linked (its owner may not be its author); one whose name is not its
  own (<Date:>_<owner>_<slug of Subject:>, so a renamed note is caught -- and the date shown is
  the Date: line's); one in a folder where someone else could delete or rename it -- a folder
  must be sticky AND owned by a Core admin or by the note's sender (a sticky folder's owner may
  still delete what is in it).
* `ack` acts only on a note `check` shows (trusted, unread, in one of my folder and all/), and
  writes nothing for any other.
* A read receipt or a reply is by its file's owner. `replies` names the owner, and flags one
  whose owner is not the user its file name gives.
* A sticky folder should stop anyone else deleting a note; that was not testable with one
  account on /quobyte, so `send` also keeps the sender's own record, outside the shared folder:
  notes_sent.jsonl in the config dir of the computer the note was sent from (id, recipient,
  subject, sha256). `replies` flags a recorded note that has gone missing from the inbox or
  whose text no longer matches -- that way only: a note with no record (older than the record,
  or sent from another computer) is never flagged for it.
* `check` also reports a recipient folder that is not sticky or not an admin's, and `replies`
  anything in skill_notes/ the skill does not put there. `say` is words for the user only.

WRITES ON /quobyte
------------------
flock does not lock across HIVE nodes on /quobyte (measured 2026-09-24: 578 of 800 writes
lost), so no two people ever write one file. Every note, receipt and reply is its own file,
written under a dot-name and os.replace'd into place, so a reader sees all of it or nothing.
Folders are 3770 (setgid + sticky: only a file's owner may delete or rename it) and files 0640,
in the skill_notes folder's group (proteomics-grp). Both modes are set with chmod after the
write: on /quobyte a new folder follows the umask (umask 027 gave drwxr-s---, measured
2026-09-29), so neither the umask nor inheritance gives the sticky or group-write bit.

WHERE IT RUNS (the routing of report_issue.sh and record_run.py)
----------------------------------------------------------------
  on HIVE (the group folder is there)        -> the inbox is read and written directly
  elsewhere, with a HIVE login (hive.env)    -> ONE hive_exec.sh call per command. This file goes
                                                over stdin (`python3 - check ...`), so it works
                                                when the HIVE copy of the skill predates notes.py
  elsewhere, with no HIVE login              -> nothing at all: status no_hive, exit 0, silent

  notes.py check [--session DIR] [--json]       unread notes for me and for all/ (read-only)
  notes.py ack NOTE_ID [--reply TEXT | --reply-file F] [--session DIR] [--no-slack] [--dry-run]
  notes.py send --to USER|all --subject S [--body TEXT | --body-file F | stdin] [--session DIR]
                [--no-slack] [--dry-run]
  notes.py replies [--since YYYY-MM-DD] [--json]  the notes I sent: who read them, and replies
  notes.py list [--json]                          a short status

Every command prints one JSON object with --json. Exit 0: done, with notes or none, and also
no_hive / no_access, which the agent says nothing about (their `say` is null). A missing
skill_notes is NOT quiet: the inbox is live, so for anyone who can see the group folder it is
reported (say: tell the Core). Exit 5: the inbox could
not be reached -- never fatal: say so in one line and carry on. Exit 2: a usage error, or a
send / ack that could not be done. `send` and `ack --reply` post one line to the Core's Slack
channel through notify_slack.py (kind `note`); the body of a note is never posted. --dry-run
writes and sends nothing and shows what would be. Guide: references/notes.md.
"""
import os
import sys

# Run from stdin (`python3 - ...`), Python sets __file__ to "<stdin>" and argv[0] to "-" (3.9,
# 3.13: measured). That is how stdin is told -- never by whether a file of that name exists,
# which anyone could plant in the working folder. Python also puts the working folder first on
# sys.path, and on HIVE that is the user's home; it goes before anything else is imported --
# subprocess alone imports select and signal -- so no file there can stand in for a module (the
# hop also runs -I).
_STDIN = globals().get("__file__") == "<stdin>" or sys.argv[:1] == ["-"]
if _STDIN and sys.path and sys.path[0] in ("", "."):
    del sys.path[0]

import argparse  # noqa: E402
import base64  # noqa: E402
import datetime  # noqa: E402
import getpass  # noqa: E402
import hashlib  # noqa: E402
import json  # noqa: E402
import posixpath  # noqa: E402
import re  # noqa: E402
import shlex  # noqa: E402
import shutil  # noqa: E402
import stat  # noqa: E402
import subprocess  # noqa: E402
import tempfile  # noqa: E402
import threading  # noqa: E402
import zlib  # noqa: E402

# `python3 - check ...` (the HIVE side of a call from a laptop) reads this file from stdin, and
# has no folder of its own: no sibling module, no core_admins.txt beside it (_STDIN, above).
_FILE = None if _STDIN else globals().get("__file__")
HERE = os.path.dirname(os.path.abspath(_FILE)) if _FILE else None
if HERE and HERE not in sys.path:
    sys.path.insert(0, HERE)

GROUP_ROOT = "/quobyte/proteomics-grp"
#: The user's copy of the skill on HIVE (hive_exec.sh --put-skill). Used only from stdin on HIVE,
#: to find notify_slack.py for the Slack line (references/notes.md, "Windows").
HIVE_COPY_SCRIPTS = os.path.join("~", "proteomics-pipeline", "scripts")
ALL = "all"
SENDERS = "SENDERS"
#: THE trust anchor, beside this file: the Core admins, one HIVE user per line with `#` comments.
#: Only they may own skill_notes/, its folders and SENDERS. Read here through skill_version.py's
#: core_admins() (skill_version.sh's skill_core_admins in bash; tests/test_skill_version.py keeps
#: the two equal). Changing it is a code change, shipped with the skill.
CORE_ADMINS_FILE = "core_admins.txt"
DIR_MODE, FILE_MODE = 0o3770, 0o640
#: The sender's own record of what they sent (send appends, replies reads): outside the shared
#: folder, in the skill's config dir -- save_transcript.config_dir()'s convention.
OUTBOX = "notes_sent.jsonl"
MAX_BODY = 16000            # characters: a note is for a person to read in one go
MAX_REPLY = 4000
MAX_SUBJECT = 200
MAX_READ = 65536            # bytes read from any one file
#: ssh to HIVE (hive_exec.sh's ConnectTimeout is 20 s) plus the work there -- notify_slack's relay.
HOP_TIMEOUT_S = 90
#: A stalled /quobyte mount must not hold the session: a direct read gives up after this.
DIRECT_TIMEOUT_S = 60
#: The line the HIVE side answers with. hive_exec.sh runs a login shell, whose own chatter may
#: come first, so the answer is found by this prefix rather than by being all of stdout.
MARK = "NOTES_JSON "

USER_RE = re.compile(r"^[A-Za-z0-9_-]{1,64}$")
TS = r"\d{8}T\d{6}Z"
ID = TS + r"_[A-Za-z0-9_-]{1,160}"
NOTE_ID_RE = re.compile(r"^" + ID + r"$")
NOTE_FILE_RE = re.compile(r"^(" + ID + r")\.md$")
#: A read receipt: <id>.read.<user>.<ts> -- a name of its own, never replacing another's -- or
#: <id>.read.<user> as the first version wrote it (still read).
RECEIPT_RE = re.compile(r"^(" + ID + r")\.read\.([A-Za-z0-9_-]{1,64})(?:\.(" + TS + r"))?$")
REPLY_RE = re.compile(r"^(" + ID + r")\.reply\.([A-Za-z0-9_-]{1,64})\.(" + TS + r")\.md$")
#: A note's Date: line, as send writes it -- the only form its name can be checked against.
DATE_RE = re.compile(r"^(\d{4})-(\d{2})-(\d{2})T(\d{2}):(\d{2}):(\d{2})Z$")

STATUS_EXIT = {"ok": 0, "no_hive": 0, "no_access": 0, "unreachable": 5, "error": 2}
QUIET = ("no_hive", "no_access")


class NotesError(Exception):
    """A send / ack that cannot be done as asked (exit 2), with the reason in plain words."""


# ---------------------------------------------------------------------------- small helpers
def notes_root():
    return os.environ.get("SKILL_NOTES_DIR") or posixpath.join(GROUP_ROOT, "skill_notes")


def outbox_path():
    """The sender's own record, in this home: SKILL_CONFIG_DIR or ~/.config/ucdavis-proteomics
    (save_transcript.config_dir(); from stdin on HIVE there is no sibling to import it from)."""
    d = os.environ.get("SKILL_CONFIG_DIR") or os.path.join("~", ".config", "ucdavis-proteomics")
    return os.path.join(os.path.expanduser(d), OUTBOX)


def sha256_of(data):
    return hashlib.sha256(data).hexdigest()


def record_sent(path, rec):
    """Append one note's record to the outbox at `path`. Never fails the send (the note is
    already there): the reason it was not recorded, or None."""
    try:
        os.makedirs(os.path.dirname(path), mode=0o700, exist_ok=True)
        fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o600)
        try:
            os.write(fd, (json.dumps(rec, sort_keys=True) + "\n").encode("utf-8"))
        finally:
            os.close(fd)
        return None
    except OSError as e:
        return f"the note was left, but not recorded in {path} ({type(e).__name__})"


def read_outbox(path):
    """{(to, id): record} from the outbox at `path`; a line that is not a record is skipped."""
    recs = {}
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                try:
                    r = json.loads(ln)
                except ValueError:
                    continue
                if (isinstance(r, dict) and NOTE_ID_RE.match(str(r.get("id") or ""))
                        and (r.get("to") == ALL or USER_RE.match(str(r.get("to") or "")))):
                    recs[(r["to"], r["id"])] = r
    except OSError:
        pass
    return recs


#: How many of the laptop's newest records a `replies` call carries to HIVE: one ssh argument
#: (Linux caps one at 128 KiB), zlib + base64.
MAX_CARRIED = 400


def pack_outbox(recs, since=None):
    """The laptop's records for HIVE, newest MAX_CARRIED first: (base64 text or None, how many
    older ones were left out)."""
    since8 = (since or "").replace("-", "")
    keep = sorted((r for r in recs.values() if not since8 or r["id"][:8] >= since8),
                  key=lambda r: r["id"], reverse=True)
    rows = [[r["to"], r["id"], r.get("sha256"), (r.get("subject") or "")[:80],
             r.get("session")] for r in keep[:MAX_CARRIED]]
    if not rows:
        return None, 0
    blob = zlib.compress(json.dumps(rows).encode("utf-8"), 9)
    return base64.b64encode(blob).decode("ascii"), max(0, len(keep) - MAX_CARRIED)


def unpack_outbox(b64):
    """pack_outbox's text -> {(to, id): record}; anything malformed is dropped."""
    try:
        rows = json.loads(zlib.decompress(base64.b64decode(b64.encode("ascii"))).decode("utf-8"))
    except (ValueError, zlib.error, TypeError):
        return {}
    recs = {}
    for row in rows if isinstance(rows, list) else []:
        if not (isinstance(row, list) and len(row) == 5):
            continue
        to, nid, sha, subject, session = row
        if NOTE_ID_RE.match(str(nid)) and (to == ALL or USER_RE.match(str(to))):
            recs[(to, nid)] = {"to": to, "id": nid, "sha256": sha, "subject": subject,
                               "session": session}
    return recs


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0)


def stamp(dt):
    return dt.strftime("%Y%m%dT%H%M%SZ")


def iso(dt):
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def stamp_iso(ts):
    """20260929T231000Z -> 2026-09-29T23:10:00Z."""
    return f"{ts[0:4]}-{ts[4:6]}-{ts[6:8]}T{ts[9:11]}:{ts[11:13]}:{ts[13:15]}Z"


def one_line(text, cap):
    s = " ".join(str(text or "").split())
    return s if len(s) <= cap else s[:cap - 1] + "..."


def slug(subject):
    s = re.sub(r"[^a-z0-9]+", "-", str(subject).lower()).strip("-")[:40].strip("-")
    return s or "note"


def owner_name(uid):
    """The login name for a uid. `pwd` exists only on HIVE / Linux / macOS, and only the HIVE
    side calls this -- a Windows laptop (Git Bash) never reads the inbox itself."""
    try:
        import pwd
        return pwd.getpwuid(uid).pw_name
    except (ImportError, KeyError, OverflowError):
        return f"uid {uid}"          # never a valid user name, so never a trusted sender


def admin_lines(text):
    """skill_version.admin_lines(), for the copy of this file sent over stdin, which has no
    sibling to import it from. tests/test_notes.py keeps the two equal."""
    out = []
    for ln in text.split("\n"):
        ln = ln.split("#", 1)[0].rstrip(" \t\n\r\x0b\x0c").lstrip(" \t\n\r\x0b\x0c")
        ln = ln.replace("\r", "")
        if ln:
            out.append(ln)
    return out


def load_core_admins(path):
    """The Core admins in the core_admins.txt at `path`: skill_version.core_admins() when this
    runs as a file, the copy above from stdin. Missing, unreadable or empty is the empty set: no
    admins, so nothing is trusted (fail closed)."""
    if _FILE:
        try:
            import skill_version
            return frozenset(skill_version.core_admins(path))
        except (ImportError, AttributeError):
            pass
    try:
        with open(path, "rb") as fh:
            return frozenset(admin_lines(fh.read(MAX_READ).decode("utf-8", "replace")))
    except OSError:
        return frozenset()


def admins_arg(admins):
    """The admins as the HIVE side of a hop gets them: comma-joined (a user name has no comma)."""
    return ",".join(sorted(admins))


def parse_admins_arg(text):
    return frozenset(t for t in (text or "").split(",") if USER_RE.match(t))


def core_admins(arg=None, remote_hop=False):
    """This run's Core admins. On the HIVE side of a laptop's call: only what the laptop sent
    (--core-admins, from its own core_admins.txt; honoured there and nowhere else). Run as a
    file: the core_admins.txt beside it. From stdin on HIVE without a hop (the Windows path,
    references/notes.md): the one in the user's HIVE copy of the skill. Nothing readable: the
    empty set, and nothing trusted."""
    if remote_hop:
        return parse_admins_arg(arg)
    folder = HERE or os.path.expanduser(HIVE_COPY_SCRIPTS)
    return load_core_admins(os.path.join(folder, CORE_ADMINS_FILE))


def current_user():
    """Who runs this, by uid (not $USER, which a shell can set to anything)."""
    try:
        import pwd
        return pwd.getpwuid(os.getuid()).pw_name
    except Exception:  # noqa: BLE001
        return getpass.getuser()


def known_account(user):
    """False only when HIVE says there is no such account; None when that cannot be asked."""
    try:
        import pwd
    except ImportError:
        return None
    try:
        pwd.getpwnam(user)
        return True
    except KeyError:
        return False


def _open_plain(path):
    """A read-only fd that never follows a symlink and never blocks on a FIFO."""
    return os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0))


def parse_note(text):
    """({lowercased header: value}, body): header lines up to the first blank line."""
    lines = text.replace("\r\n", "\n").split("\n")
    headers, i = {}, 0
    while i < len(lines) and lines[i].strip():
        k, sep, v = lines[i].partition(":")
        if sep:
            headers.setdefault(k.strip().lower(), v.strip())
        i += 1
    return headers, "\n".join(lines[i + 1:]).strip("\n")


def same_session(a, b):
    """Does a note's Session: name this session? The same folder, or the same folder name (a
    sender may write just the name)."""
    if not a or not b:
        return False
    na, nb = (posixpath.normpath(x.strip().replace("\\", "/")) for x in (a, b))
    return na == nb or posixpath.basename(na) == posixpath.basename(nb)


def classify(names):
    """A folder listing -> ({note id}, {(id, user, file name)} receipts, {(id, user, ts)}
    replies). Dot-names (a write in progress) and anything else are left out."""
    notes, receipts, replies = set(), set(), set()
    for n in names:
        if n.startswith("."):
            continue
        m = NOTE_FILE_RE.match(n)
        if m:
            notes.add(m.group(1))
            continue
        m = RECEIPT_RE.match(n)
        if m:
            receipts.add((m.group(1), m.group(2), n))
            continue
        m = REPLY_RE.match(n)
        if m:
            replies.add((m.group(1), m.group(2), m.group(3)))
    return notes, receipts, replies


def safe(text):
    """Someone else's text -- a file name, an owner that is no plain user name -- as it may be
    shown: every control and non-ASCII character escaped, and at most 40 characters of it. A
    newline in a name planted in the shared folder must not become a line of its own in what
    Brett's Claude reads."""
    return ascii(str(text))[:40]


def as_user(name):
    """An owner as shown: a plain user name as it is, anything else through safe()."""
    return name if USER_RE.match(str(name)) else safe(name)


def _listdir(d):
    try:
        return os.listdir(d)
    except FileNotFoundError:
        return []


# ------------------------------------------------------------------------------- the inbox
class Inbox:
    """The inbox on this filesystem, read as `me`. `owner_of(path, stat_result)` names a file's
    owner (default: its uid, through pwd); the tests pass their own to play several people.
    `outbox` is my own record of the notes I sent (default: outbox_path()); `admins` the Core
    admins (default: core_admins())."""

    def __init__(self, root, me, owner_of=None, outbox=None, admins=None):
        self.root = root.rstrip("/") or root
        self.me = me
        self._owner_of = owner_of or (lambda path, st: owner_name(st.st_uid))
        self.outbox = outbox or outbox_path()
        self.admins = frozenset(core_admins() if admins is None else admins)

    def owner(self, path, st):
        return self._owner_of(path, st)

    # ---- trust
    def root_problem(self):
        """Why nothing in skill_notes can be trusted, or None. It must be a folder itself (never
        a link: lstat), sticky, and a Core admin's -- /quobyte/proteomics-grp is writable by the
        whole group, so anyone could rename it and put another of their own in its place."""
        if not self.admins:
            return ("no Core admins are listed (scripts/core_admins.txt is missing, unreadable "
                    "or empty), so no note is shown")
        st = os.lstat(self.root)
        if stat.S_ISLNK(st.st_mode):
            return "skill_notes is a link, not the folder itself, so no note in it is shown"
        if not stat.S_ISDIR(st.st_mode):
            return "skill_notes is not a folder, so no note is shown"
        owner = self.owner(self.root, st)
        if owner not in self.admins:
            return (f"skill_notes belongs to {owner}, who is not a Core admin "
                    f"({', '.join(sorted(self.admins))}: scripts/core_admins.txt), so no note "
                    "in it is shown")
        if not st.st_mode & stat.S_ISVTX:
            return (f"skill_notes is not sticky, so anyone in the group could rename what is in "
                    f"it, and no note is shown until {owner} runs chmod 3770 on it")
        return None

    def trusted(self):
        """(trusted sender names, why SENDERS was ignored or None): the Core admins, plus
        SENDERS when it is one plain file, a Core admin's, that nobody else can write."""
        names = set(self.admins)
        p = os.path.join(self.root, SENDERS)
        try:
            fd = _open_plain(p)
        except FileNotFoundError:
            return names, None
        except OSError as e:
            return names, f"SENDERS ignored: it is a link or cannot be opened ({type(e).__name__})"
        try:
            st = os.fstat(fd)                    # before fdopen, which refuses a folder itself
            if not stat.S_ISREG(st.st_mode):
                os.close(fd)
                return names, "SENDERS ignored: not a plain file"
            with os.fdopen(fd, "rb") as fh:
                text = fh.read(MAX_READ).decode("utf-8", "replace")
        except OSError as e:
            return names, f"SENDERS ignored: it cannot be read ({type(e).__name__})"
        who = self.owner(p, st)
        if who not in self.admins:
            return names, f"SENDERS ignored: it belongs to {who}, who is not a Core admin"
        if st.st_nlink != 1:
            return names, ("SENDERS ignored: it is hard-linked, so it may not be what its owner "
                           "wrote")
        if st.st_mode & 0o022:
            return names, "SENDERS ignored: people other than its owner can write it (chmod 644)"
        for ln in text.splitlines():
            tok = ln.split("#", 1)[0].strip()
            if USER_RE.match(tok):
                names.add(tok)
        return names, None

    def folder_guard(self, folder):
        """(the folder's owner, is it sticky?) -- (None, False) for anything but a real folder
        (a link, a file). FileNotFoundError when it is not there."""
        d = os.path.join(self.root, folder)
        st = os.lstat(d)
        if not stat.S_ISDIR(st.st_mode):
            return None, False
        return self.owner(d, st), bool(st.st_mode & stat.S_ISVTX)

    def folder_problem(self, guard, sender):
        """Why someone other than `sender` could have deleted, renamed or swapped a note of
        theirs in a folder with this guard, or None. It must be sticky -- then only a file's
        owner, and the folder's, may do that -- and owned by a Core admin or by the sender
        themself, so the folder's owner is no one else either. A folder a member made first
        (all/, or the folder of someone new) is theirs, and so is not trusted."""
        owner, sticky = guard
        if owner is None:
            return "its folder is a link or not a folder"
        if not sticky:
            return (f"its folder is not sticky, so others could delete or rename what is in it "
                    f"({owner} can run chmod 3770 on it)")
        if owner not in self.admins and owner != sender:
            return (f"its folder belongs to {owner}, neither a Core admin nor its sender, who "
                    "could delete or rename what is in it")
        return None

    def root_entries(self):
        """What in skill_notes/ is not what the skill puts there -- a dot-name, a link, a file
        (but SENDERS), a name that is no user -- as one line per kind: how many, and a few of
        their names through safe(). Anyone in the group can name a file there, so a name is
        never shown as it is; and nothing is skipped silently."""
        kinds = {}
        for name in sorted(os.listdir(self.root)):
            try:
                st = os.lstat(os.path.join(self.root, name))
            except OSError:
                kinds.setdefault("entries that cannot be read", []).append(name)
                continue
            if name == SENDERS:
                if not stat.S_ISREG(st.st_mode):
                    kinds.setdefault("SENDERS that is not a plain file", []).append(name)
            elif name.startswith("."):
                kinds.setdefault("dot-names, which the notes never use", []).append(name)
            elif stat.S_ISLNK(st.st_mode):
                kinds.setdefault("links", []).append(name)
            elif not stat.S_ISDIR(st.st_mode):
                kinds.setdefault("files where folders belong", []).append(name)
            elif name != ALL and not USER_RE.match(name):
                kinds.setdefault("folders not named for a user", []).append(name)
        return [f"skill_notes/ holds {len(names)} {kind}: "
                + ", ".join(safe(n) for n in names[:5]) + (" and more" if len(names) > 5 else "")
                for kind, names in sorted(kinds.items())]

    # ---- reading
    def _mine(self, d, receipts, note_id):
        """True when a read receipt for `note_id` in folder `d` is a plain file of mine -- in
        either form, <id>.read.<me>.<ts> or the first version's <id>.read.<me>."""
        for rid, user, name in receipts:
            if rid != note_id or user != self.me:
                continue
            p = os.path.join(d, name)
            try:
                st = os.lstat(p)
            except OSError:
                continue
            if stat.S_ISREG(st.st_mode) and self.owner(p, st) == self.me:
                return True
        return False

    def read_note(self, folder, note_id, guard=None):
        """The note as a dict. `problems` (not shown) and `warnings` (shown, flagged) say what
        is wrong with it; `sender` is its owner. `guard`: folder_guard(folder), when known.

        Its name must be its own: <UTC Date:>_<owner>_<slug of Subject:>, with send's -2, -3
        when two share a second. A name that is not -- a renamed note, an old one given a new
        date or words -- is not shown, and the date shown is the Date: line's."""
        path = os.path.join(self.root, folder, note_id + ".md")
        fd = _open_plain(path)
        st = os.fstat(fd)                       # the file read, not whatever is at the path now
        if stat.S_ISREG(st.st_mode):
            with os.fdopen(fd, "rb") as fh:
                raw = fh.read(MAX_READ + 1)
        else:
            os.close(fd)
            raw = b""
        owner = self.owner(path, st)
        problems, warnings = [], []
        if not stat.S_ISREG(st.st_mode):
            problems.append("not a plain file")
        if st.st_mode & 0o022:
            problems.append("people other than its owner can write it, so its owner may not "
                            "be its author")
        if st.st_nlink > 1:
            problems.append("it is hard-linked, so its owner may not be its author")
        why = self.folder_problem(guard or self.folder_guard(folder), owner)
        if why:
            problems.append(why)
        headers, body = parse_note(raw[:MAX_READ].decode("utf-8", "replace"))
        missing = [h for h in ("From", "To", "Date", "Subject") if not headers.get(h.lower())]
        if missing:
            problems.append("not a note (no " + ", ".join(missing) + " line)")
        m = DATE_RE.match(headers.get("date") or "")
        own = f"{''.join(m.groups()[:3])}T{''.join(m.groups()[3:])}Z_{owner}_" if m else None
        own = own and own + slug(headers.get("subject") or "")
        if not (own and (note_id == own or (note_id.startswith(own + "-")
                                             and note_id[len(own) + 1:].isdigit()))):
            problems.append("its file name is not its own (<Date:>_<owner>_<Subject:>), so it "
                            "may have been renamed")
        if headers.get("from") and headers["from"] != owner:
            warnings.append(f"its From: line says {headers['from']}, but the file belongs to "
                            f"{owner}, who is the sender")
        if headers.get("to") and headers["to"] not in (folder, ALL):
            warnings.append(f"its To: line says {headers['to']}, but it is in {folder}/")
        if len(raw) > MAX_READ or len(body) > MAX_BODY:
            body = body[:MAX_BODY]
            warnings.append("longer than this shows: cut")
        return {"id": note_id, "folder": folder, "sender": owner,
                "from": headers.get("from"), "to": headers.get("to"),
                "date": headers.get("date") or None, "subject": headers.get("subject"),
                "session": headers.get("session") or None, "body": body,
                "warnings": warnings, "problems": problems, "path": path,
                "sha256": sha256_of(raw)}

    def _folders(self, errors):
        """(folder, guard) for mine and all/, each named in `errors` when it is not usable
        (not there is fine): a bad entry is reported, never the end of the check."""
        out = []
        for folder in (self.me, ALL):
            try:
                guard = self.folder_guard(folder)
            except FileNotFoundError:
                continue
            except OSError as e:
                errors.append(f"skill_notes/{folder} cannot be read ({type(e).__name__})")
                continue
            if guard[0] is None:
                errors.append(f"skill_notes/{folder} is a link or a file, not a folder: nothing "
                              "in it is shown")
                continue
            if not guard[1]:
                errors.append(f"skill_notes/{folder}/ is not sticky, so nothing in it is shown "
                              f"until its owner, {guard[0]}, runs chmod 3770 on it")
            elif guard[0] not in self.admins:
                errors.append(f"skill_notes/{folder}/ belongs to {guard[0]}, who is not a Core "
                              f"admin: only {guard[0]}'s own notes there can be shown")
            out.append((folder, guard))
        return out

    def check(self, session=None):
        """Unread notes for me and for all/, from trusted senders; read-only."""
        untrusted = self.root_problem()
        if untrusted:
            return {"status": "ok", "user": self.me, "unread": [], "not_shown": [],
                    "trusted_senders": [], "senders_file": None, "inbox_trusted": False,
                    "problems": [untrusted], "say": untrusted}
        trusted, senders_why = self.trusted()
        unread, ignored, errors = [], {}, [senders_why] if senders_why else []
        for folder, guard in self._folders(errors):
            d = os.path.join(self.root, folder)
            try:
                notes, receipts, _ = classify(os.listdir(d))
            except OSError as e:
                errors.append(f"skill_notes/{folder}/ cannot be listed ({type(e).__name__})")
                continue
            for nid in sorted(notes):
                if self._mine(d, receipts, nid):
                    continue
                try:
                    n = self.read_note(folder, nid, guard)
                except FileNotFoundError:
                    continue                       # withdrawn while we looked
                except Exception as e:  # noqa: BLE001 -- one bad file never ends the check
                    key = ("?", f"cannot be read ({type(e).__name__})")
                    ignored[key] = ignored.get(key, 0) + 1
                    continue
                why = self.verdict(n, trusted)
                if why:
                    key = (n["sender"], why)
                    ignored[key] = ignored.get(key, 0) + 1
                    continue
                n["for_this_session"] = same_session(n["session"], session)
                for k in ("problems", "sha256"):
                    del n[k]
                unread.append(n)
        unread.sort(key=lambda n: (not n["for_this_session"], n["id"]))
        not_shown = [{"owner": as_user(o), "count": c, "why": w}
                     for (o, w), c in sorted(ignored.items())]
        return {"status": "ok", "user": self.me, "unread": unread, "not_shown": not_shown,
                "trusted_senders": sorted(trusted), "senders_file": senders_why,
                "inbox_trusted": True, "problems": errors,
                "say": check_say(unread, not_shown, errors)}

    @staticmethod
    def verdict(n, trusted):
        """Why note `n` is not shown, or None: what is wrong with it, or who it is from."""
        if n["problems"]:
            return "; ".join(n["problems"])
        if n["sender"] not in trusted:
            return "not a trusted sender"
        return None

    # ---- writing
    def _finish(self, path, mode):
        """The mode, set explicitly -- on /quobyte a new folder follows the umask, so neither
        the sticky nor the group-write bit comes by itself -- and the skill_notes folder's group.
        Never fails the write: where setgid is refused (not a member of the folder's group) the
        rest of the mode is still set."""
        for m in (mode, mode & ~stat.S_ISGID, mode & 0o777):
            try:
                os.chmod(path, m)
                break
            except OSError:
                continue
        try:
            gid = os.stat(self.root).st_gid
            if os.stat(path).st_gid != gid:
                os.chown(path, -1, gid)
        except (OSError, AttributeError):
            pass

    def _folder(self, name, create=True):
        """The recipient's folder: made 3770 when missing, and made sticky when it is mine and
        is not. One that fails folder_problem() for me would hide every note put in it, so a
        note is refused rather than sent there to be lost -- and the reason says who can fix
        it."""
        d = os.path.join(self.root, name)
        if not os.path.lexists(d):
            if not create:
                return d
            try:
                os.mkdir(d, DIR_MODE)
            except FileExistsError:
                pass                    # someone made it meanwhile: judged below, like any other
            else:
                self._finish(d, DIR_MODE)          # made here, so mine: only the mode to check
                if not os.lstat(d).st_mode & stat.S_ISVTX:
                    raise NotesError(f"skill_notes/{name}/ was made but could not be made sticky "
                                     "(chmod 3770), so a note put there would never be shown")
                return d
        guard = self.folder_guard(name)
        if guard[0] == self.me and not guard[1] and create:
            self._finish(d, DIR_MODE)
            guard = self.folder_guard(name)
        why = self.folder_problem(guard, self.me)
        if why:
            owner, sticky = guard
            fix = ("a Core admin must look at it" if owner is None else
                   f"ask {owner} to run: chmod 3770 {d}" if not sticky else
                   f"{owner} made it: ask them to move it away, so that it can be made again")
            raise NotesError(f"skill_notes/{name}/: {why}, so a note put there would never be "
                             f"shown; {fix}")
        return d

    def _put(self, d, name, text):
        """Write `text` to <d>/<name> through a dot-name and a rename: all of it or none."""
        fd, tmp = tempfile.mkstemp(prefix="." + name + ".", suffix=".part", dir=d)
        try:
            with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(text)
            self._finish(tmp, FILE_MODE)
            os.replace(tmp, os.path.join(d, name))
        except BaseException:
            try:
                os.unlink(tmp)
            except OSError:
                pass
            raise
        return os.path.join(d, name)

    def ack(self, note_id, reply=None, session=None, skill_version=None, dry_run=False):
        """Mark a note read (a receipt), and write the reply when there is one -- only for a
        note `check` shows: trusted, and unread. Anything else is refused and nothing is
        written. An id in both my folder and all/ is the copy that passes the trust checks;
        when both do, neither is acknowledged (send never makes such a pair itself)."""
        if not NOTE_ID_RE.match(note_id or ""):
            raise NotesError(f"'{note_id}' is not a note id "
                             "(e.g. 20260929T231000Z_brettsp_keratin)")
        untrusted = self.root_problem()
        if untrusted:
            raise NotesError(untrusted)
        where = [f for f in (self.me, ALL)
                 if os.path.lexists(os.path.join(self.root, f, note_id + ".md"))]
        if not where:
            raise NotesError(f"no note {note_id} for {self.me} (in skill_notes/{self.me}/ or "
                             "all/)")
        trusted, _ = self.trusted()
        passing, why = [], []
        for f in where:
            try:
                n = self.read_note(f, note_id)
            except Exception as e:  # noqa: BLE001
                why.append(f"{f}/: cannot be read ({type(e).__name__})")
                continue
            w = self.verdict(n, trusted)
            if w:
                why.append(f"{f}/: {w}" if len(where) > 1 else w)
            else:
                passing.append((f, n))
        if len(passing) > 1:
            raise NotesError(f"{note_id} is in both skill_notes/{self.me}/ and all/, and both "
                             "look like notes from the Core: neither is acknowledged. Tell the "
                             "Core (report_issue.sh)")
        if not passing:
            raise NotesError(f"{note_id} is not a note from the Core that is shown "
                             f"({'; '.join(why)}): nothing is written")
        folder, n = passing[0]
        d = os.path.join(self.root, folder)
        if self._mine(d, classify(os.listdir(d))[1], note_id):
            raise NotesError(f"{note_id} is already acknowledged")
        now = utc_now()
        # a name of its own: never replaces a file, anyone's, under another receipt's name
        rpath = os.path.join(d, f"{note_id}.read.{self.me}.{stamp(now)}")
        receipt = {"note": note_id, "folder": folder, "user": self.me, "read_at": iso(now),
                   "session": session or None, "skill_version": skill_version or None,
                   "replied": bool(reply)}
        reply_path = None
        if reply:
            reply_path = os.path.join(d, f"{note_id}.reply.{self.me}.{stamp(now)}.md")
        if not dry_run:
            if reply:
                text = (f"From: {self.me}\nTo: {n['sender']}\nDate: {iso(now)}\n"
                        f"In-Reply-To: {note_id}\n"
                        f"Subject: Re: {one_line(n['subject'] or '', MAX_SUBJECT)}\n"
                        + (f"Session: {one_line(session, 500)}\n" if session else "")
                        + "\n" + reply.rstrip("\n") + "\n")
                self._put(d, os.path.basename(reply_path), text)
            self._put(d, os.path.basename(rpath), json.dumps(receipt, indent=2) + "\n")
        return {"status": "ok", "note": note_id, "folder": folder, "user": self.me,
                "sender": n["sender"], "subject": n["subject"], "trusted": True,
                "receipt": rpath, "reply": reply_path, "dry_run": dry_run, "say": None}

    def send(self, to, subject, body, session=None, dry_run=False, record=True):
        """Leave a note for `to` (a HIVE user, or `all`). Only a trusted sender may: a note
        from anyone else is never shown, so sending one would only look like it worked.
        `record`: add it to my outbox here -- False on the HIVE side of a laptop's call, whose
        own outbox (on the laptop) records it from the `outbox_record` returned."""
        if not (to == ALL or USER_RE.match(to or "")):
            raise NotesError(f"--to must be a HIVE user name or 'all', not '{to}'")
        subject = one_line(subject, MAX_SUBJECT)
        if not subject:
            raise NotesError("a note needs a --subject")
        body = (body or "").strip("\n")
        if not body.strip():
            raise NotesError("a note needs a body (--body, --body-file, or text on stdin)")
        if len(body) > MAX_BODY:
            raise NotesError(f"the body is {len(body)} characters; keep a note under {MAX_BODY}")
        untrusted = self.root_problem()
        if untrusted:
            raise NotesError(untrusted + ": it must be a sticky folder a Core admin owns "
                             "(references/notes.md, setup)")
        trusted, senders_why = self.trusted()
        if self.me not in trusted:
            raise NotesError(
                f"{self.me} is not a trusted sender, so the note would never be shown. A Core "
                f"admin ({', '.join(sorted(self.admins))}) can add you to skill_notes/SENDERS"
                + (f" ({senders_why})" if senders_why else "") + ".")
        warnings = []
        if to != ALL and known_account(to) is False:
            warnings.append(f"there is no HIVE account named {to} here; check the spelling")
        d = self._folder(to, create=not dry_run)
        now = utc_now()
        base = f"{stamp(now)}_{self.me}_{slug(subject)}"
        note_id, k = base, 1
        # not an id already in use here or in all/ (or, for all/, in anyone's folder): the same
        # id in two folders a reader sees could never be acked
        taken = {d, os.path.join(self.root, ALL)}
        if to == ALL:
            taken |= {os.path.join(self.root, x) for x in os.listdir(self.root)
                      if USER_RE.match(x)}
        while any(os.path.lexists(os.path.join(f, note_id + ".md")) for f in taken):
            k += 1
            note_id = f"{base}-{k}"
        text = (f"From: {self.me}\nTo: {to}\nDate: {iso(now)}\nSubject: {subject}\n"
                + (f"Session: {one_line(session, 500)}\n" if session else "")
                + "\n" + body + "\n")
        path = os.path.join(d, note_id + ".md")
        rec = {"id": note_id, "to": to, "subject": subject, "session": session or None,
               "sent_at": iso(now), "path": path, "sha256": sha256_of(text.encode("utf-8"))}
        if not dry_run:
            self._put(d, note_id + ".md", text)
            why = record_sent(self.outbox, rec) if record else None
            if why:
                warnings.append(why)
        return {"status": "ok", "id": note_id, "to": to, "sender": self.me, "subject": subject,
                "session": session or None, "path": path, "warnings": warnings,
                "outbox_record": rec, "dry_run": dry_run, "say": "; ".join(warnings) or None}

    # ---- the sender's view
    def _author(self, path, user, folder):
        """(the file's owner, why it does not count or None) for a receipt or a reply named for
        `user`: it counts when `user` owns it and may answer a note in `folder` (its recipient,
        or anyone for all/)."""
        try:
            st = os.lstat(path)
        except OSError as e:
            return None, f"cannot be read ({type(e).__name__})"
        if not stat.S_ISREG(st.st_mode):
            return None, "not a plain file"
        owner = self.owner(path, st)
        if owner != user:
            return owner, f"its file is named for {user}, but it belongs to {as_user(owner)}"
        if folder not in (ALL, user):
            return owner, f"{user} is not this note's recipient"
        return owner, None

    def replies(self, since=None, extra=None):
        """The notes I sent: who has read each one (their own receipts), what they replied, and
        anything wrong with it. A receipt or a reply counts only when its owner is the user its
        name gives, and, for a note to one person, only from that person; any other is listed
        with its owner and a `flag`. A note in my outbox that is missing from the inbox, or no
        longer the text I sent, is flagged too -- only that way round: a note with no record is
        listed, never flagged for it. `extra`: records from the outbox of the laptop calling.
        `problems`: anything in skill_notes/ that is not what the skill puts there."""
        since8 = (since or "").replace("-", "")
        recorded = read_outbox(self.outbox)
        recorded.update(extra or {})
        problems = self.root_entries()
        sent, seen = [], set()
        for folder in sorted(os.listdir(self.root)):
            d = os.path.join(self.root, folder)
            if folder.startswith(".") or folder == SENDERS or not (
                    folder == ALL or USER_RE.match(folder)):
                continue                           # named in `problems` by root_entries()
            try:
                guard = self.folder_guard(folder)
                if guard[0] is None:
                    continue                       # a link or a file: named in `problems`
                notes, receipts, replies = classify(os.listdir(d))
            except OSError as e:
                problems.append(f"skill_notes/{folder}/ cannot be listed ({type(e).__name__})")
                continue
            for nid in sorted(notes):
                if since8 and nid[:8] < since8:
                    continue
                try:
                    n = self.read_note(folder, nid, guard)
                except OSError:
                    continue
                if n["sender"] != self.me:
                    continue
                seen.add((folder, nid))
                rec = recorded.get((folder, nid))
                flags = list(n["problems"])      # what makes `check` skip it for its reader
                if rec and rec.get("sha256") != n["sha256"]:
                    flags.append("changed since it was sent: its text is not what was sent")
                # Each receipt and reply on its own: a bad one (planted to break this) is named
                # in `problems` -- by its user and time only -- and the rest still read.
                marks, read_by = [], []
                for rid, user, name in sorted(receipts):
                    if rid != nid:
                        continue
                    try:
                        p = os.path.join(d, name)
                        owner, why = self._author(p, user, folder)
                        marks.append({"named": user, "owner": as_user(owner) if owner else None,
                                      "read_at": _receipt_time(p), "flag": why})
                        if not why and user not in {x["user"] for x in read_by}:
                            read_by.append({"user": user, "read_at": marks[-1]["read_at"]})
                    except Exception as e:  # noqa: BLE001
                        problems.append(f"a receipt for {nid} named for {user} cannot be read "
                                        f"({type(e).__name__})")
                answers = []
                for rid, user, ts in sorted(replies):
                    if rid != nid:
                        continue
                    try:
                        p = os.path.join(d, f"{nid}.reply.{user}.{ts}.md")
                        owner, why = self._author(p, user, folder)
                        text = None
                        if owner is not None:
                            with os.fdopen(_open_plain(p), "rb") as fh:
                                _, text = parse_note(fh.read(MAX_READ).decode("utf-8", "replace"))
                            text = text[:MAX_REPLY]
                        answers.append({"from": as_user(owner) if owner else None,
                                        "named": user, "date": stamp_iso(ts), "text": text,
                                        "flag": why})
                    except Exception as e:  # noqa: BLE001
                        problems.append(f"a reply to {nid} named for {user} at {ts} cannot be "
                                        f"read ({type(e).__name__})")
                sent.append({"id": nid, "to": folder, "subject": n["subject"],
                             "date": n["date"], "session": n["session"],
                             "recorded": rec is not None, "flags": flags, "read_by": read_by,
                             "receipts": marks, "replies": answers,
                             "unread": folder != ALL and not read_by})
        for (to, nid), rec in sorted(recorded.items()):
            if (to, nid) in seen or (since8 and nid[:8] < since8):
                continue
            p = os.path.join(self.root, to, nid + ".md")
            try:
                st = os.lstat(p)
            except FileNotFoundError:
                flag = "missing from the inbox: deleted by someone?"
            except OSError as e:
                flag = f"cannot be read ({type(e).__name__})"
            else:
                who = self.owner(p, st)
                flag = (f"replaced: the file there now belongs to {who}" if who != self.me
                        else "cannot be read as a note")
            sent.append({"id": nid, "to": to, "subject": rec.get("subject"),
                         "date": stamp_iso(nid[:16]), "session": rec.get("session"),
                         "recorded": True, "flags": [flag], "read_by": [], "receipts": [],
                         "replies": [], "unread": None})
        sent.sort(key=lambda e: e["id"])
        flagged = [f"{e['id']} to {e['to']}: {x}" for e in sent
                   for x in e["flags"] + [m["flag"] for m in e["receipts"] + e["replies"]
                                          if m["flag"]]]
        untrusted = self.root_problem()          # then nobody is shown any of these notes
        if untrusted:
            problems.insert(0, untrusted)
        return {"status": "ok", "user": self.me, "sent": sent, "flagged": len(flagged),
                "inbox_trusted": not untrusted, "outbox": self.outbox, "problems": problems,
                "say": ("; ".join(problems + flagged)[:1000] or None)}

    def status(self, extra=None):
        c = self.check()
        r = self.replies(extra=extra)
        trusted = c["trusted_senders"]
        problems = c["problems"] + [p for p in r["problems"] if p not in c["problems"]]
        return {"status": "ok", "user": self.me, "unread": len(c["unread"]),
                "not_shown": sum(x["count"] for x in c["not_shown"]),
                "sent": len(r["sent"]), "sent_unread": sum(1 for x in r["sent"] if x["unread"]),
                "replies": sum(len(x["replies"]) for x in r["sent"]),
                "flagged": r["flagged"],
                "can_send": self.me in trusted, "trusted_senders": trusted,
                "problems": problems,
                "say": "; ".join(x for x in (c["say"], r["say"]) if x) or None}


def _receipt_time(path):
    """A receipt's read_at, or None. Whatever is in the file -- 30,000 nested brackets make
    json.loads raise RecursionError -- it never raises."""
    try:
        with os.fdopen(_open_plain(path), "rb") as fh:
            t = json.loads(fh.read(MAX_READ).decode("utf-8", "replace")).get("read_at")
        return t if isinstance(t, str) and DATE_RE.match(t) else None
    except Exception:  # noqa: BLE001
        return None


def check_say(unread, not_shown, errors):
    """One line for the user -- words for them, never instructions to the agent -- or None when
    there is nothing to say."""
    bits = []
    if unread:
        who = ", ".join(sorted({n["sender"] for n in unread}))
        bits.append(f"{len(unread)} unread note(s) from the Core ({who})")
    for x in not_shown:
        bits.append(f"{x['count']} note(s) from {x['owner']} not shown: {x['why']}")
    bits.extend(errors)
    return "; ".join(bits) or None


# ------------------------------------------------------------------------------- routing
def hive_login():
    """The HIVE user when this computer can reach HIVE through hive_exec.sh, else None:
    record_run.hive_login(), the skill's own test for it (the environment, then hive.env)."""
    if not _FILE:
        return None                       # from stdin on HIVE: no second hop from here
    try:
        import record_run
        return record_run.hive_login()
    except Exception:  # noqa: BLE001 -- no sibling (from stdin): no second hop from here
        return None


def route(remote_hop=False):
    """direct (the group folder is here: HIVE), ssh (a HIVE login), or none."""
    if os.path.isdir(posixpath.dirname(notes_root().rstrip("/"))):
        return "direct"
    if not remote_hop and hive_login():
        return "ssh"
    return "none"


def hop(argv):
    """This file, run on HIVE with `argv` through hive_exec.sh -- ONE ssh call. (answer, None),
    or (None, why) when HIVE did not answer."""
    import record_run
    hx = record_run.hive_exec_path()
    if not _FILE or not shutil.which("bash"):
        return None, "bash is not available here to run hive_exec.sh"
    with open(_FILE, "rb") as fh:
        src = fh.read()
    # The inbox override, when set, goes along (the tests); in use it is not, and HIVE's
    # default applies there.
    env = " ".join(f"{k}={shlex.quote(os.environ[k])}" for k in ("SKILL_NOTES_DIR",)
                   if os.environ.get(k))
    # -I: HIVE's python imports nothing from the working folder (the user's home) or PYTHON*.
    cmd = (env + " " if env else "") + "python3 -I - " + " ".join(shlex.quote(a) for a in argv)
    # Its own process group, so a timeout ends ssh as well as bash: on Windows (Git Bash) ssh
    # outlives a killed bash and keeps the pipes open, and communicate() would wait on it.
    kw = ({"creationflags": getattr(subprocess, "CREATE_NEW_PROCESS_GROUP", 0)}
          if os.name == "nt" else {"start_new_session": True})
    try:
        p = subprocess.Popen(["bash", hx, cmd], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE, **kw)
    except OSError as e:
        return None, f"could not run hive_exec.sh ({type(e).__name__})"
    try:
        so, se = p.communicate(src, timeout=HOP_TIMEOUT_S)
    except subprocess.TimeoutExpired:
        _kill_tree(p)
        return None, f"no answer from HIVE within {HOP_TIMEOUT_S} s"
    r = subprocess.CompletedProcess(p.args, p.returncode, so, se)
    out = r.stdout.decode("utf-8", "replace").replace("\r", "")
    for ln in reversed(out.splitlines()):
        if ln.startswith(MARK):
            try:
                ans = json.loads(ln[len(MARK):])
            except ValueError:
                break
            if isinstance(ans, dict):
                return ans, None
    err = [x.strip() for x in r.stderr.decode("utf-8", "replace").splitlines() if x.strip()]
    return None, (err[-1] if err else f"hive_exec.sh exited {r.returncode}")[:200]


def _kill_tree(p):
    """End a hop that ran out of time: its whole process group (bash, ssh, whatever they
    started), then collect it briefly -- a pipe still held open is left, not waited on."""
    try:
        if os.name == "nt":
            subprocess.run(["taskkill", "/F", "/T", "/PID", str(p.pid)],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=15)
        else:
            import signal
            os.killpg(p.pid, signal.SIGKILL)
    except (OSError, subprocess.SubprocessError):
        pass
    try:
        p.kill()
    except OSError:
        pass
    try:
        p.communicate(timeout=5)
    except (subprocess.TimeoutExpired, OSError, ValueError):
        pass
    for f in (p.stdin, p.stdout, p.stderr):
        try:
            f.close()
        except (OSError, ValueError, AttributeError):
            pass


def _bounded(fn, seconds):
    """fn() on a daemon thread, given up after `seconds` (a stalled mount must not hold the
    session): its result, or None on timeout. An exception in fn is raised here."""
    box = {}

    def work():
        try:
            box["r"] = fn()
        except BaseException as e:  # noqa: BLE001 -- handed back to the caller
            box["e"] = e

    t = threading.Thread(target=work, daemon=True)
    t.start()
    t.join(seconds)
    if t.is_alive():
        return None
    if "e" in box:
        raise box["e"]
    return box["r"]


def _in_child(fn, seconds, fallback=True):
    """fn() in a forked child that is killed at the deadline: its result, or None on timeout.

    A stalled /quobyte mount can hold a stat() in the kernel (D state), where no thread can be
    given up on for good; a child can be left behind. Everything that touches this computer's
    files runs in it -- route()'s isdir as well. fn's result must be JSON; an exception in fn
    comes back as {"status": "unreachable"} with its name. Where there is no fork (Windows,
    which never reads the inbox itself), a daemon thread does the same, less surely.

    The child answers only through its pipe: its stdin, stdout and stderr are /dev/null, so a
    child that cannot be killed (D state) holds nothing of the caller's open -- whoever reads
    this process's output is not kept waiting for it. A child that dies without an answer is
    run again here, under the thread's deadline, only when `fallback` (a read: check, list,
    replies); a send or an ack may have written already, so it is never done twice -- None."""
    if not hasattr(os, "fork"):
        try:
            return _bounded(fn, seconds)
        except Exception as e:  # noqa: BLE001
            return {"status": "unreachable", "say": f"{type(e).__name__}: {e}"[:300]}
    try:
        import pwd  # noqa: F401 -- imported before the fork, not inside a child that may hang
    except ImportError:
        pass
    import select
    import signal
    import time
    rfd, wfd = os.pipe()
    pid = os.fork()
    if pid == 0:                                       # the child: do it, write it, go
        code = 0
        try:
            os.close(rfd)
            null = os.open(os.devnull, os.O_RDWR)
            for fd in (0, 1, 2):
                os.dup2(null, fd)
            try:
                res = fn()
            except BaseException as e:  # noqa: BLE001
                res = {"status": "unreachable", "say": f"{type(e).__name__}: {e}"[:300]}
            data = json.dumps(res).encode("utf-8")
            while data:
                data = data[os.write(wfd, data):]
        except BaseException:  # noqa: BLE001
            code = 1
        finally:
            os._exit(code)
    os.close(wfd)
    chunks, deadline = [], time.monotonic() + seconds
    try:
        while True:
            left = deadline - time.monotonic()
            if left <= 0 or not select.select([rfd], [], [], left)[0]:
                try:
                    os.kill(pid, signal.SIGKILL)
                    os.waitpid(pid, os.WNOHANG)
                except OSError:
                    pass
                return None
            chunk = os.read(rfd, 65536)
            if not chunk:
                break
            chunks.append(chunk)
    finally:
        os.close(rfd)
    status = None
    try:
        status = os.waitpid(pid, 0)[1]
    except OSError:
        pass
    try:
        return json.loads(b"".join(chunks).decode("utf-8"))
    except ValueError:
        pass
    # It ended without an answer: a crash, not a hang (macOS aborts a forked child -- status 6
    # -- once the parent's threads have set up its system frameworks). A read is done here
    # instead, still under a deadline, rather than lose the check; a write is not repeated.
    if not fallback:
        return None
    try:
        res = _bounded(fn, seconds)
    except Exception as e:  # noqa: BLE001
        res = {"status": "unreachable", "_route": "direct",
               "say": f"the notes check failed: {type(e).__name__}: {e}"[:300]}
    if isinstance(res, dict):
        res.setdefault("child_status", status)
    return res


def local_work(a, body=None, reply=None):
    """Everything this command does to this computer's files, in one piece for _in_child():
    the route, then -- when the inbox is here -- the command on it. {"_route": ...} added."""
    where = route(a.remote_hop)
    if where != "direct":
        return {"_route": where}
    res = run_direct(a.cmd, a, body, reply)
    res["_route"] = "direct"
    return res


def run_direct(cmd, a, body=None, reply=None):
    """Run one command on the inbox on this filesystem. Always returns a result dict."""
    root = notes_root()
    quiet = {"status": "no_access", "say": None,
             "detail": f"{root} is readable only by the Proteomics Core"}
    try:
        os.stat(root)
    except FileNotFoundError:
        # The inbox is live (since 2026-09-30), so a missing one was moved or replaced -- never
        # "not set up yet". Reaching here, the group folder is there and this account can see
        # into it (a Core account): loud, for check, list and replies alike.
        return missing_root(cmd, f"{root} is missing: it may have been moved or replaced. "
                                 "Tell the Core")
    except PermissionError:
        return quiet
    if not os.access(root, os.R_OK | os.X_OK):
        return quiet
    box = Inbox(root, current_user(), admins=core_admins(a.core_admins, a.remote_hop))
    extra = unpack_outbox(a.expect_b64) if getattr(a, "expect_b64", None) else None
    work = {
        "check": lambda: box.check(a.session),
        "list": lambda: box.status(extra),
        "replies": lambda: box.replies(a.since, extra),
        "ack": lambda: box.ack(a.note_id, reply, a.session, a.skill_version or _skill_version(),
                               a.dry_run),
        # a laptop's note is recorded in the laptop's own outbox, once HIVE has written it
        "send": lambda: box.send(a.to, a.subject, body, a.session, a.dry_run,
                                 record=not a.remote_hop),
    }[cmd]
    try:
        res = work()
    except NotesError as e:
        return {"status": "error", "say": str(e)}
    except PermissionError as e:
        if cmd in ("send", "ack"):
            return {"status": "error", "say": f"permission denied: {e.filename or root}"}
        return {"status": "unreachable", "say": f"the notes folder could not be read ({e})"}
    except Exception as e:  # noqa: BLE001 -- the inbox never takes down the session
        return {"status": "unreachable" if cmd not in ("send", "ack") else "error",
                "say": f"the notes folder could not be used: {type(e).__name__}: {e}"[:300]}
    return res


def missing_root(cmd, msg):
    """The answer when skill_notes is not there, in each command's own shape."""
    if cmd in ("send", "ack"):
        return {"status": "error", "say": msg}
    base = {"status": "ok", "user": current_user(), "inbox_trusted": False, "problems": [msg],
            "trusted_senders": [], "say": msg}
    if cmd == "check":
        return dict(base, unread=[], not_shown=[], senders_file=None)
    if cmd == "replies":
        return dict(base, sent=[], flagged=0, outbox=outbox_path())
    return dict(base, unread=0, not_shown=0, sent=0, sent_unread=0, replies=0, flagged=0,
                can_send=False)


def _skill_version():
    """This copy's version (skill_version.py, the one reader), or None where it cannot be read
    -- from stdin on HIVE the laptop passes its own with --skill-version."""
    if not _FILE:
        return None
    try:
        import skill_version
        return skill_version.skill_version()
    except Exception:  # noqa: BLE001
        return None


_NS = []


def _notify_slack():
    """notify_slack.py, the skill's one Slack client, or None. Run as a file: the sibling. Run
    from stdin on HIVE (a Windows laptop with no usable python3: references/notes.md): that
    file in the user's HIVE copy of the skill, loaded by its path and nothing else."""
    if _NS:
        return _NS[0]
    mod = None
    if _FILE:
        try:
            import notify_slack as mod
        except Exception:  # noqa: BLE001
            mod = None
    else:
        p = os.path.join(os.path.expanduser(HIVE_COPY_SCRIPTS), "notify_slack.py")
        if os.path.isfile(p):
            try:
                import importlib.util
                spec = importlib.util.spec_from_file_location("notify_slack", p)
                mod = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(mod)
            except Exception:  # noqa: BLE001
                mod = None
    _NS.append(mod)
    return mod


def refuse_secrets(*texts):
    """A note or reply lands in a folder the whole Core can read: text shaped like a key, token,
    password or webhook is refused, as report_issue.sh does, with notify_slack's one list."""
    ns = _notify_slack()
    if ns is None or not hasattr(ns, "contains_secret"):
        return None
    if any(ns.contains_secret(t) for t in texts if t):
        return ("the text looks like it contains a key, token, password or webhook -- remove "
                "it and send again")
    return None


def slack_line(event, res, a, where, reply=None):
    """The Slack mirror of a send or a reply, through notify_slack.py. Never fatal: returns the
    status line (or, with --dry-run, the payload that would be posted)."""
    if a.no_slack:
        return "off (--no-slack)"
    ns = _notify_slack()
    if ns is None or not hasattr(ns, "note_facts"):
        return ("not sent: this copy of notify_slack.py has no note posts (put the skill on "
                "HIVE again)")
    try:
        if event == "sent":
            facts = ns.note_facts("sent", sender=res.get("sender"), to=res.get("to"),
                                  subject=res.get("subject"))
        else:
            facts = ns.note_facts("reply", sender=res.get("sender"), to=None,
                                  subject=res.get("subject"), reply=reply,
                                  replier=res.get("user"))
        if a.dry_run:
            return {"dry_run": True, "payload": ns.payload(facts)}
        _, msg = ns.deliver(facts, relay_ok=(where == "ssh"))
        return msg
    except Exception as e:  # noqa: BLE001 -- a Slack problem never fails the note
        return f"not sent: {type(e).__name__}"


# ------------------------------------------------------------------------------------ CLI
def _read_text(path):
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read()


def _b64(s):
    return base64.b64encode(s.encode("utf-8")).decode("ascii")


def _unb64(s):
    return base64.b64decode(s.encode("ascii")).decode("utf-8", "replace")


def parser():
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--json", action="store_true", help="print one JSON object")
    common.add_argument("--remote-hop", action="store_true", help=argparse.SUPPRESS)
    common.add_argument("--core-admins", help=argparse.SUPPRESS)   # the laptop's, over a hop
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd")
    c = sub.add_parser("check", parents=[common], help="unread notes for me (read-only)")
    c.add_argument("--session", help="this session's folder: its notes come first")
    k = sub.add_parser("ack", parents=[common], help="mark a note read, with a reply")
    k.add_argument("note_id")
    g = k.add_mutually_exclusive_group()
    g.add_argument("--reply", help="one line: what the user decided, or what was done")
    g.add_argument("--reply-file")
    g.add_argument("--reply-b64", help="the reply, base64 of its UTF-8 -- for a shell that "
                   "would otherwise read its $, backticks or quotes (references/notes.md)")
    k.add_argument("--session")
    k.add_argument("--skill-version", help=argparse.SUPPRESS)
    k.add_argument("--no-slack", action="store_true", help="post nothing to Slack")
    k.add_argument("--dry-run", action="store_true", help="write and send nothing; show it")
    s = sub.add_parser("send", parents=[common], help="leave a note (trusted senders only)")
    s.add_argument("--to", required=True, help="a HIVE user name, or all")
    sg = s.add_mutually_exclusive_group(required=True)
    sg.add_argument("--subject")
    sg.add_argument("--subject-b64", help=argparse.SUPPRESS)
    bg = s.add_mutually_exclusive_group()
    bg.add_argument("--body")
    bg.add_argument("--body-file")
    bg.add_argument("--body-b64", help=argparse.SUPPRESS)
    s.add_argument("--session", help="the session folder the note is about")
    s.add_argument("--no-slack", action="store_true", help="post nothing to Slack")
    s.add_argument("--dry-run", action="store_true", help="write and send nothing; show it")
    r = sub.add_parser("replies", parents=[common], help="the notes I sent: read by, replies")
    r.add_argument("--since", help="YYYY-MM-DD")
    r.add_argument("--expect-b64", help=argparse.SUPPRESS)       # the laptop's outbox records
    li = sub.add_parser("list", parents=[common], help="a short status")
    li.add_argument("--expect-b64", help=argparse.SUPPRESS)
    return ap


def emit(res, a):
    """Print the result (JSON, or a few lines for a person) and return the exit code."""
    code = STATUS_EXIT.get(res.get("status"), 2)
    if a.remote_hop:
        print(MARK + json.dumps(res))
        return 0
    if a.json:
        print(json.dumps(res, indent=2))
        return code
    if res.get("status") in QUIET and a.cmd in ("check", "list", "replies"):
        return code                                   # nothing to say to anyone
    text = render_text(res, a.cmd)
    if text:
        sys.stdout.write(text.encode(sys.stdout.encoding or "utf-8", "replace")
                         .decode(sys.stdout.encoding or "utf-8") + "\n")
    return code


def render_text(res, cmd):
    if res.get("status") != "ok":
        return f"notes: {res.get('status')}: {res.get('say') or ''}".rstrip(": ")
    out = []
    if cmd == "check":
        for n in res["unread"]:
            flag = " -- about this session" if n.get("for_this_session") else ""
            out += [f"Note from {n['sender']}, {n['date']}{flag}", f"Subject: {n['subject']}",
                    f"id: {n['id']}"] + [f"(!) {w}" for w in n["warnings"]] + ["", n["body"], "--"]
        for x in res["not_shown"]:
            out.append(f"{x['count']} note(s) from {x['owner']} not shown: {x['why']}")
        out += res.get("problems") or []
    elif cmd == "replies":
        for n in res["sent"]:
            who = ", ".join(x["user"] for x in n["read_by"]) or "nobody yet"
            out.append(f"{n['date']}  to {n['to']}: {n['subject']}  (read by {who})  [{n['id']}]")
            out += [f"    (!) {f}" for f in n["flags"]]
            out += [f"    (!) receipt: {m['flag']}" for m in n["receipts"] if m["flag"]]
            for x in n["replies"]:
                out.append(f"    {x['from']}, {x['date']}: {x['text']}"
                           + (f"  (!) {x['flag']}" if x["flag"] else ""))
    elif cmd == "list":
        out.append(f"{res['user']}: {res['unread']} unread; sent {res['sent']} "
                   f"({res['sent_unread']} unread by their recipient, {res['replies']} replies, "
                   f"{res['flagged']} flagged)")
        if res.get("say"):
            out.append(res["say"])
    elif cmd == "send":
        out.append(("would leave" if res.get("dry_run") else "left") +
                   f" note {res['id']} for {res['to']}")
        out += [f"(!) {w}" for w in res.get("warnings") or []]
    elif cmd == "ack":
        out.append(("would mark" if res.get("dry_run") else "marked") + f" {res['note']} read"
                   + (" and reply" if res.get("reply") else ""))
    if isinstance(res.get("slack"), str):
        out.append(f"Slack: {res['slack']}")
    elif isinstance(res.get("slack"), dict):
        out.append("Slack (dry run): " + json.dumps(res["slack"]["payload"]))
    return "\n".join(out)


def main(argv=None):
    ap = parser()
    a = ap.parse_args(argv)
    if not a.cmd:
        ap.print_help(sys.stderr)
        return 2
    for k in ("session", "since", "note_id", "to", "no_slack", "dry_run", "skill_version",
              "expect_b64"):
        if not hasattr(a, k):
            setattr(a, k, None)

    # The free text, read HERE: on the far side of a hop stdin is this file, not the body.
    body = reply = None
    try:
        if a.cmd == "send":
            a.subject = _unb64(a.subject_b64) if a.subject_b64 else a.subject
            if a.body_b64:
                body = _unb64(a.body_b64)
            elif a.body is not None:
                body = a.body
            elif a.body_file:
                body = _read_text(a.body_file)
            elif not a.remote_hop and sys.stdin is not None and not sys.stdin.isatty():
                body = sys.stdin.read()
        elif a.cmd == "ack":
            reply = (_unb64(a.reply_b64) if a.reply_b64 else
                     _read_text(a.reply_file) if a.reply_file else a.reply)
            reply = (reply or "").strip() or None
            if reply and len(reply) > MAX_REPLY:
                return emit({"status": "error",
                             "say": f"keep a reply under {MAX_REPLY} characters"}, a)
    except OSError as e:
        return emit({"status": "error", "say": f"cannot read {e.filename}: {e.strerror}"}, a)
    if a.since and not re.match(r"^\d{4}-\d{2}-\d{2}$", a.since):
        return emit({"status": "error", "say": "--since takes YYYY-MM-DD"}, a)
    if not a.remote_hop:
        why = refuse_secrets(getattr(a, "subject", None), body, reply)
        if why:
            return emit({"status": "error", "say": why}, a)

    res = _in_child(lambda: local_work(a, body, reply), DIRECT_TIMEOUT_S,
                    fallback=a.cmd in ("check", "list", "replies"))
    if res is None:
        res = {"status": "unreachable", "_route": "direct",
               "say": f"the inbox did not answer; whether the {a.cmd} was written is not "
                      "known: check before trying again" if a.cmd in ("send", "ack") else
                      f"the notes folder did not answer within {DIRECT_TIMEOUT_S} s"}
    where = res.pop("_route", "direct")              # direct: `res` is the answer already
    if where == "ssh":
        # HIVE has no copy of this file to read: the admins go along, from this computer's own
        # core_admins.txt, and are all the HIVE side trusts
        tail = [a.cmd, "--json", "--remote-hop", "--core-admins", admins_arg(core_admins())]
        if a.session:
            tail += ["--session", a.session]
        if a.cmd == "replies" and a.since:
            tail += ["--since", a.since]
        if a.cmd == "ack":
            tail += [a.note_id, "--skill-version", _skill_version() or ""]
            if reply:
                tail += ["--reply-b64", _b64(reply)]
        if a.cmd == "send":
            tail += ["--to", a.to, "--subject-b64", _b64(a.subject or "")]
            if body is not None:
                tail += ["--body-b64", _b64(body)]
        if a.cmd in ("send", "ack") and a.dry_run:
            tail.append("--dry-run")
        left_out = 0
        if a.cmd in ("replies", "list"):
            packed, left_out = pack_outbox(read_outbox(outbox_path()),
                                           a.since if a.cmd == "replies" else None)
            if packed:
                tail += ["--expect-b64", packed]
        res, why = hop(tail)
        if res is None:
            res = {"status": "unreachable",
                   "say": f"could not reach the notes on HIVE ({why}); carrying on"}
        elif a.cmd == "send" and res.get("status") == "ok" and not a.dry_run \
                and isinstance(res.get("outbox_record"), dict):
            why = record_sent(outbox_path(), res["outbox_record"])      # this laptop's record
            if why:
                res["warnings"] = (res.get("warnings") or []) + [why]
                res["say"] = "; ".join(res["warnings"])
        if left_out:
            res["outbox_left_out"] = left_out           # older records than one call can carry
    elif where == "none":
        res = {"status": "no_hive", "say": None}
        if a.cmd in ("send", "ack"):
            res = {"status": "error", "say": "no HIVE here and no HIVE login (hive.env): the "
                                             "notes live on HIVE"}
    if a.cmd in ("send", "ack") and res.get("status") in QUIET:
        res["status"] = "error"                 # an action that did not happen is not quiet
        res["say"] = res.get("say") or res.get("detail") or "the notes are not reachable here"
    if a.cmd == "check":
        res.setdefault("unread", [])

    # The Slack line, from this side of any hop: only here is notify_slack.py importable.
    if not a.remote_hop and res.get("status") == "ok":
        if a.cmd == "send":
            res["slack"] = slack_line("sent", res, a, where)
        elif a.cmd == "ack" and reply and res.get("trusted"):
            res["slack"] = slack_line("reply", res, a, where, reply=reply)
    return emit(res, a)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001 -- checking notes must never end a session
        hop_side = "--remote-hop" in sys.argv       # the laptop finds the answer by its mark
        print((MARK if hop_side else "") + json.dumps(
            {"status": "unreachable", "say": f"notes.py failed: {type(e).__name__}: {e}"[:300]}))
        sys.exit(0 if hop_side else 5)
