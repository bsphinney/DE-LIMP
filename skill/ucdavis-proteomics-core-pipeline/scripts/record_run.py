#!/usr/bin/env python3
"""
record_run.py  --  Record every search the UC Davis Proteomics Core runs with this skill in the
Core's run registry on HIVE, laid out the way Brett's DataAnalysis sessions are.

Why this exists: a search leaves its evidence scattered -- search_provenance.json in the search
folder, the FASTA sidecar beside the FASTA, the parameters in wf/, the SLURM logs wherever the job
wrote them, the report and the session zip on the instrument share or a laptop. "What did the
Core search last month, for which submission, with which DIA-NN, and how deep did it go?" meant
opening each of them, if anyone still knew where they were. So every search gets one session
folder in the registry, written the moment it ends (or fails) and updated when it is finalized,
plus a line in an append-only master log and rows in an activity log.

  # when a search ends -- the job-end hook calls this, and so can the agent (ON HIVE)
  python3 record_run.py search-done --out <search out dir> --status completed|failed \
      --exit-code N [--session <session dir>] [--prot PROT_0807]
  # after `session.py finalize` -- the finalize hook calls this
  python3 record_run.py analysis-done --session <session dir> [--out <dir>] [--zip <zip>] [--prot ..]
  # the registry, from every run_record.json (the logs are for people; this is the index)
  python3 record_run.py list [--since 2026-09-01] [--user gabrig] [--status failed] [--tsv|--json]
  python3 record_run.py --where        # which route and destination from here, and nothing else
  ... --dry-run                        # what it would write, written nowhere (preview on stderr)

stdout is ONE JSON object -- {"recorded": true, "path": "<session folder>", ...} or
{"recorded": false, "reason": "<code>", "detail": "..."} -- and the exit status is 0 whatever
happened; only an argparse usage error exits non-zero. Diagnostics go to stderr.

Layout (DataAnalysis "Session Structure"):
  <SKILL_RUNS_DIR>/
    data_analysis.md        MASTER LOG: one entry per run, later events appended as dated lines
    activity_log.csv        timestamp,session,action,tool,target,status,notes -- one row per action
    sessions/<YYYY-MM-DD_Short-Description>/     the skill's session name (else <date>_<out name>)
      SEARCH_LOG.md         the run: CoreOmics submission, Data Quality Notes (always present),
                            who/when/status, data, engine + version, key parameters and their
                            sources, headline results, where everything is, FRAN, skill issues
      run_record.json       the same, machine-readable (schema_version)
      README.md             the session README (after finalize)
      <session>.zip         the session zip, minus per-run .quant files, when under the cap
      input/                conditions.csv, FASTA sidecar, params + rationale, raw_files.txt -- never
                            raw data
      output/               *.docx (the report of record), methods.md, AUDIT/SAMPLE_QUALITY, tables/,
                            small figures/, and search/ (provenance, stats, engine + SLURM logs,
                            FRAN receipt, and a LINK to report.parquet)
      scripts/              commands.log, reproduce.sh, REPRODUCE.md
    .index/                 out_<16 hex> / session_<16 hex> -> sessions/<name>: how a re-record
                            finds its folder. Folder names carry no hash; the identity (the search
                            out dir's real path) is in run_record.json, and a DIFFERENT search that
                            wants the same name gets <name>_2.

Where it goes -- the routing report_issue.sh uses:
  1. On HIVE, registry writable (a Proteomics Core account)       -> written directly.
  2. On HIVE, registry NOT writable                                -> not recorded (not_core_member).
     This is the gate, and it is a permission rather than a flag: a collaborator's HIVE account
     is outside proteomics-grp, their runs are theirs, and the filesystem keeps them out.
  3. Off HIVE with a HIVE login (hive.env -- HIVE_ENV_FILE -- or HIVE_USER + HIVE_KEY) -> over SSH
     via hive_exec.sh: what exists only on this machine is uploaded with this script, and HIVE
     writes the record.
  4. Anything else                                                 -> not recorded (not_on_hive).
RECORD_RUN=off (or SKILL_RUNS_DIR=off) switches it off before anything is read, written or sent.
SKILL_RUNS_DIR is the destination on every route.

Bounded: --timeout (45 s: the job-end hook allows 60), no directory walks, a per-file cap, a
per-record copy budget. Stdlib only. Nothing named like a credential and no text that looks like
one is ever copied; raw data and per-run .quant files never are -- the record says where they are.
"""
import argparse
import copy
import csv
import datetime
import getpass
import glob
import hashlib
import io
import json
import os
import random
import re
import shlex
import shutil
import signal
import socket
import statistics
import struct
import subprocess
import sys
import tempfile
import time
import zipfile

# Piped to a remote python (`python3 - list ...`) there is no __file__: the sibling modules are
# then simply unavailable, and every caller degrades to "not parsed here" instead of failing.
_FILE = globals().get("__file__")
HERE = os.path.dirname(os.path.abspath(_FILE)) if _FILE else os.getcwd()
if _FILE and HERE not in sys.path:
    sys.path.insert(0, HERE)

SCHEMA_VERSION = 2
GROUP_ROOT = "/quobyte/proteomics-grp"
MB, GB = 1 << 20, 1 << 30
OFF = ("off", "0", "no", "false")
ACTIVITY_COLUMNS = ["timestamp", "session", "action", "tool", "target", "status", "notes"]
MASTER_HEADER = (
    "# Proteomics Core run registry -- master log\n\n"
    "> One entry per search run by the ucdavis-proteomics-core-pipeline skill for the UC Davis "
    "Proteomics Core, written by `record_run.py`. Append-only: a later event (a failure fixed, the "
    "analysis finalized) is a new dated line, never an edit. Each run's full record is "
    "`sessions/<name>/SEARCH_LOG.md`; `activity_log.csv` is the machine-readable companion.\n\n"
    "---\n")
# A search that ran through run_search.py names its out dir one of these; the name says nothing
# about WHICH search, so the fallback name climbs to the folder above.
GENERIC_DIR_NAMES = {"search_out", "search", "out", "output", "outputs", "results",
                     "search_results"}
REPORT_CANDIDATES = ("report.parquet", "report.tsv", "delimp_report.parquet",
                     "dia-quant-output/report.tsv", "radiant_results/fulcrum-results",
                     "fulcrum-results", "results.sage.parquet", "results.sage.tsv")
# Per-run DIA-NN temporaries: single-shot `quant/`, the chain's `quant_step2/` (+ the `_orig`
# resume backup) and `quant_step4/`. Tens of MB per run, three copies in a chain -- never copied.
QUANT_DIRS = ("quant", "quant_step2", "quant_step2_orig", "quant_step4")
PARAMS_FALLBACK = ("params.resolved.cfg", "params_with_xic.cfg", "params.base.cfg", "params.cfg")
# SLURM log names the skill's job scripts use: <step>_<jobid>.log, <step>_<jobid>_<task>.log.
JOBLOG_RE = re.compile(r"_(\d{5,})(?:_(\d+))?\.(?:log|out|err)$")
BAD_STATES = {"FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED", "NODE_FAIL", "BOOT_FAIL",
              "DEADLINE", "PREEMPTED"}
LIVE_STATES = {"PENDING", "RUNNING", "REQUEUED", "RESIZING", "SUSPENDED", "CONFIGURING",
               "COMPLETING"}
# HIVE puts sacct on PATH only in a login shell (hive_exec.sh runs `bash -l -c` for that reason);
# a hook may run without one. Verified 2026-09-24: `which sacct` in a HIVE login shell.
SACCT_CANDIDATES = ("/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin/sacct",
                    "/usr/bin/sacct", "/usr/local/bin/sacct", "/opt/slurm/bin/sacct")
# The registry is readable by the whole Core group. Refuse by NAME anything that holds a
# credential, and by CONTENT anything that looks like one -- the patterns report_issue.sh refuses,
# plus Slack tokens and webhook URLs. A refused file is listed with the reason, never silent.
SECRET_NAME_RE = re.compile(
    r"token|webhook|secret|passw|credential|(^|\.)env$|^id_(rsa|dsa|ecdsa|ed25519)|"
    r"\.(pem|key|p12|pfx)$|^\.netrc$|^\.pgpass$", re.I)
SECRET_TEXT_RE = re.compile(
    rb"-----BEGIN [A-Z ]*PRIVATE KEY|ghp_[A-Za-z0-9]{20,}|github_pat_|hf_[A-Za-z0-9]{20,}|"
    rb"Authorization: *(Token|Bearer) +[A-Za-z0-9]|xox[abprs]-[A-Za-z0-9-]{10,}|"
    rb"hooks\.slack\.com/services/")
# CoreOmics identifiers (DataAnalysis CLAUDE.md, "Every session links to its CoreOmics
# submission"): PROT_#### plus a 12-hex id. The same forms core_submission.py accepts.
HEX_ID_RE = re.compile(r"^[0-9a-f]{12}$")
PROT_RE = re.compile(r"^#?(?:prot)?[\s_\-]*0*(\d{1,5})$", re.I)
UPLOAD_DIR = ".proteomics-pipeline/record_run_upload"
# The shared files' lock. NOT flock: measured 2026-09-24 on /quobyte with two SLURM jobs on
# different HIVE nodes writing 400 times each to one file, flock lost 578 of 800 updates -- the
# same as no lock at all (567) -- and readers saw torn JSON; an atomic `mkdir <file>.lock.d` with
# stale-lock breaking lost 0 of 800 (FRAN ingest/auto_ingest_state.py). Waiting is bounded: past
# it the write happens unlocked and the result says `lock_timeout` rather than losing the record.
STALE_LOCK_S = 60
LOCK_TIMEOUTS = []
# The registry's own README.md, written by ensure_readme() when it is missing or carries an older
# version marker than this. THE source: references/run-registry.md quotes it verbatim, and
# tests/test_record_run.py fails when the two differ. Raise README_VERSION with any edit.
README_VERSION = 1
README_TEXT = """# Skill run registry

Every search the `ucdavis-proteomics-core-pipeline` skill runs for the UC Davis Proteomics Core,
recorded by `scripts/record_run.py` when the search ends (or fails) and updated when the analysis
is finalized. Laid out like the Core's DataAnalysis sessions.

- `data_analysis.md` -- the master log: one entry per run; later events (a failure fixed, the
  analysis finalized) are appended as dated lines, never edited in.
- `activity_log.csv` -- `timestamp,session,action,tool,target,status,notes`, one row per action
  (search submitted / started / completed / failed, FRAN staged or skipped, skill issue recorded,
  analysis completed).
- `sessions/<YYYY-MM-DD_Short-Description>/` -- one folder per search, named after the skill's
  session (else `<date>_<search folder name>`; a different search wanting the same name gets
  `_2`). Start with `SEARCH_LOG.md`: the CoreOmics submission, Data Quality Notes, status, engine
  and the version that ran, key parameters and where each came from, results, and where every
  output is. Beside it: `run_record.json` (the same, machine-readable), the session `README.md`
  and zip, `input/` (conditions, FASTA sidecar, parameters, `raw_files.txt`), `output/` (the
  Word report of record, Methods, tables, `search/` logs and a link to `report.parquet`) and
  `scripts/` (commands, reproduce script).
- `.index/` -- how recording the same search again finds its folder (by the search folder's real
  path, not its name). Leave it alone.
- `*.lock.d` -- a writer's lock, held for a second or two; one older than 60 s is broken
  automatically. Do not use `flock` on these files: it does not lock across HIVE nodes.

Never copied: raw data, FASTAs, per-run `.quant` files, anything over 20 MB, or anything that
looks like a credential -- the record says where they are.

Only Proteomics Core (proteomics-grp) accounts can write here, so collaborators' runs are never
recorded. To keep a run out: `RECORD_RUN=off`. To list the registry:
`python3 <skill>/scripts/record_run.py list [--since YYYY-MM-DD] [--user U] [--status failed]`.
Skill problems go to the sibling folder `../skill_issues/`.

<!-- record_run.py README v%d: written by the skill; edit README_TEXT in record_run.py, not this
file -->
""" % README_VERSION


# ---------------------------------------------------------------------------- small helpers
def runs_dir():
    return os.environ.get("SKILL_RUNS_DIR") or f"{GROUP_ROOT}/skill_runs"


def disabled():
    return (os.environ.get("RECORD_RUN", "").strip().lower() in OFF
            or os.environ.get("SKILL_RUNS_DIR", "").strip().lower() == "off")


def issues_dir():
    """report_issue.sh's folder: the registry's sibling, unless set on its own."""
    return (os.environ.get("SKILL_ISSUES_DIR")
            or os.path.join(os.path.dirname(runs_dir().rstrip("/")), "skill_issues"))


def env_num(name, default):
    try:
        return float(os.environ.get(name) or default)
    except ValueError:
        return float(default)


def whoami():
    try:
        return getpass.getuser()
    except Exception:
        return os.environ.get("USER") or os.environ.get("USERNAME") or "unknown"


def now_iso():
    return datetime.datetime.now().astimezone().isoformat(timespec="seconds")


def ts_min(value=None):
    """ISO 8601 to the minute with the UTC offset -- 2026-09-24T12:45-07:00 -- the activity-log
    form. `value` is an ISO string (sacct's has no offset: it is local time) or None for now."""
    dt = None
    if value:
        try:
            dt = datetime.datetime.fromisoformat(str(value)[:25])
        except ValueError:
            dt = None
    dt = (dt or datetime.datetime.now()).astimezone()
    s = dt.strftime("%Y-%m-%dT%H:%M%z")
    return s[:-2] + ":" + s[-2:]


def load_json(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError, TypeError):
        return None


def fsize(path):
    try:
        return os.path.getsize(path)
    except (OSError, TypeError):
        return None


def mtime_iso(path):
    try:
        return datetime.datetime.fromtimestamp(os.path.getmtime(path)).astimezone() \
            .isoformat(timespec="seconds")
    except (OSError, TypeError, ValueError):
        return None


def listdir(d):
    try:
        return sorted(os.listdir(d))
    except (OSError, TypeError):
        return []


def owner(path):
    try:
        import pwd
        return pwd.getpwuid(os.stat(path).st_uid).pw_name
    except Exception:
        return None


def clean(s, n=40):
    """Folder-name safe, like report_issue.sh's clean(): the name lands in a shared folder."""
    s = re.sub(r"[^A-Za-z0-9._-]+", "_", str(s or "")).strip("._-")
    return re.sub(r"_+", "_", s)[:n].rstrip("._-")


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s or "").lower())


def key16(real):
    """16 hex of a RESOLVED path: the identity a re-record is matched by."""
    return hashlib.sha1(str(real).encode()).hexdigest()[:16]


def same_or_later(a, b):
    """a >= b for two local timestamps, whether ISO with an offset or sacct's bare form."""
    return (a or "")[:19].replace("T", " ") >= (b or "")[:19].replace("T", " ")


def fmt_n(x):
    if x is None:
        return "?"
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, float) and x.is_integer():
        x = int(x)
    return f"{x:,}" if isinstance(x, int) else f"{x:g}" if isinstance(x, float) else str(x)


def fmt_bytes(n):
    if n is None:
        return "?"
    for unit, f in (("GB", GB), ("MB", MB), ("KB", 1024)):
        if n >= f:
            return f"{n / f:.1f} {unit}"
    return f"{n} B"


def skill_version():
    pj = load_json(os.path.join(HERE, "..", ".claude-plugin", "plugin.json")) or {}
    return pj.get("version")


def say(msg):
    sys.stderr.write(f"[record_run] {msg}\n")


class Deadline:
    def __init__(self, seconds):
        self.end = time.monotonic() + max(1.0, float(seconds))

    def left(self):
        return self.end - time.monotonic()

    def check(self):
        if self.left() <= 0:
            raise TimeoutError("the --timeout was reached")


class Stop(Exception):
    """The SIGALRM backstop: a hung filesystem call must not hold a SLURM job open."""


def copy_n(src, dst, n, deadline):
    while n > 0:
        deadline.check()
        b = src.read(min(n, 8 * MB))
        if not b:
            raise IOError("unexpected end of file")
        dst.write(b)
        n -= len(b)


def copy_file(src, dst, deadline):
    """Chunked, deadline-checked, written beside the target and renamed into place: a copy cut
    short by the time limit leaves no half-file posing as the real one."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    part = f"{dst}.part{os.getpid()}"
    try:
        with open(src, "rb") as fi, open(part, "wb") as fo:
            copy_n(fi, fo, os.fstat(fi.fileno()).st_size, deadline)
        os.replace(part, dst)
    finally:
        if os.path.exists(part):
            os.remove(part)
    try:
        # the source's times, but NOT its mode: a copy takes the registry's group permissions
        st = os.stat(src)
        os.utime(dst, (st.st_atime, st.st_mtime))
    except OSError:
        pass


def write_text(path, text):
    """tmp + fsync + rename: a reader on another node sees the old file or the new one, never
    half of one."""
    tmp = f"{path}.tmp.{socket.gethostname()}.{os.getpid()}"
    try:
        with open(tmp, "w") as fh:
            fh.write(text)
            fh.flush()
            os.fsync(fh.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)


def read_json_retry(path, tries=4):
    """(object, "ok" | "missing" | "unreadable") for a JSON another writer may be replacing. A
    failed parse is retried -- it can be a moment of someone else's write -- and a file that stays
    unparseable is reported, never deleted or overwritten here."""
    last = None
    for i in range(tries):
        try:
            with open(path) as fh:
                return json.load(fh), "ok"
        except FileNotFoundError:
            return None, "missing"
        except (ValueError, UnicodeDecodeError, OSError) as e:
            last = e
            time.sleep(0.2 * (i + 1))
    say(f"{path} is unreadable after {tries} tries ({last})")
    return None, "unreadable"


def secret_reason(path):
    """Why this file must not be copied into a group-readable folder, or None."""
    if SECRET_NAME_RE.search(os.path.basename(path)):
        return "its name marks it as a credential (token / webhook / .env / key file)"
    try:
        with open(path, "rb") as fh:
            if SECRET_TEXT_RE.search(fh.read()):
                return "its contents look like a key, token or webhook URL"
    except OSError as e:
        return f"unreadable: {e}"
    return None


def json_tail(text, opener):
    """The JSON document at the end of `text` (a login shell may print before it)."""
    m = re.search(r"^" + re.escape(opener) + r".*\Z", text or "", re.S | re.M)
    try:
        return json.loads(m.group(0)) if m else None
    except ValueError:
        return None


# --------------------------------------------------------------- append-only shared files
class DirLock:
    """`mkdir <target>.lock.d` -- the lock that holds across HIVE nodes on /quobyte (see
    STALE_LOCK_S). Retries with jitter; waits at most RECORD_RUN_LOCK_WAIT seconds (10), then
    proceeds UNLOCKED and records the file in LOCK_TIMEOUTS -- a late record is better than none.

    Each acquisition writes a unique owner token into the lock directory (FRAN ad81863, the same
    lock in ingest/auto_ingest_state.py). Breaking a stale lock is judge -> rename -> VERIFY: the
    breaker notes whose lock it judged dead, renames the directory away, and checks the renamed
    directory still carries that owner. If not, the lock changed hands in between -- A died, C broke
    A's lock and took a fresh one, and B, who judged A's, renamed C's LIVE lock -- so B puts it back
    and waits like anyone else. Release removes the directory only if its owner is still us: a
    holder that outlived STALE_LOCK_S may have had its lock broken and re-taken by someone else."""

    def __init__(self, target, deadline, wait=None):
        self.target, self.dir, self.deadline = target, target + ".lock.d", deadline
        self.wait = env_num("RECORD_RUN_LOCK_WAIT", 10) if wait is None else wait
        self.held = False
        self.token = None
        self.break_hook = None               # tests: interleave another process in a break

    @staticmethod
    def owner_of(d):
        try:
            with open(os.path.join(d, "owner")) as fh:
                return fh.read().strip() or None
        except OSError:
            return None

    def __enter__(self):
        me = f"{socket.gethostname()}:{os.getpid()}"
        stop = time.monotonic() + max(0.0, min(self.wait, self.deadline.left() - 3))
        while True:
            token = f"{me}:{time.time_ns()}"
            try:
                os.mkdir(self.dir)
            except FileExistsError:
                pass
            except OSError as e:
                say(f"cannot lock {self.dir} ({e}); writing unlocked")
                LOCK_TIMEOUTS.append(os.path.basename(self.target))
                return self
            else:
                self.held = True
                try:
                    with open(os.path.join(self.dir, "owner"), "w") as fh:
                        fh.write(token + "\n")
                    self.token = token
                except OSError:
                    self.token = None            # held, but unmarked: see __exit__
                return self
            try:
                age = time.time() - os.stat(self.dir).st_mtime
            except OSError:
                continue                          # released between mkdir and stat: retry
            if age > STALE_LOCK_S:
                judged = self.owner_of(self.dir)
                if self.break_hook is not None:
                    self.break_hook(judged)
                grave = f"{self.dir}.stale.{me.replace(':', '.')}.{time.time_ns()}"
                try:
                    os.rename(self.dir, grave)
                except OSError:
                    continue                      # someone else broke or released it first
                if self.owner_of(grave) == judged:
                    shutil.rmtree(grave, ignore_errors=True)
                    say(f"broke a lock {age:.0f} s old held by {judged}: {self.dir}")
                else:
                    try:
                        os.rename(grave, self.dir)
                        say(f"the lock changed hands while being broken; restored it to its live "
                            f"holder: {self.dir}")
                    except OSError as e:
                        say(f"moved a live lock aside by mistake and could not restore it ({e}); "
                            f"it is at {grave}")
                continue
            if time.monotonic() >= stop:
                say(f"{self.dir} still held after {self.wait:.0f} s; writing unlocked")
                LOCK_TIMEOUTS.append(os.path.basename(self.target))
                return self
            time.sleep(0.05 + random.random() * 0.1)

    def __exit__(self, *exc):
        if self.held:
            owner = self.owner_of(self.dir)
            if owner == self.token:               # ours (or unmarked ours: both None)
                shutil.rmtree(self.dir, ignore_errors=True)
            else:
                say(f"lock {self.dir} now belongs to {owner}, not to this writer; left in place")
            self.held, self.token = False, None
        return False


def create_with_header(path, header):
    """Create `path` holding exactly `header`, or leave an existing file alone. Written to a
    temporary name and hard-linked into place, so even an unlocked writer can never append to a
    file that does not have its header yet."""
    if os.path.exists(path):
        return
    tmp = f"{path}.new.{socket.gethostname()}.{os.getpid()}"
    with open(tmp, "w") as fh:
        fh.write(header)
        fh.flush()
        os.fsync(fh.fileno())
    try:
        os.link(tmp, path)
    except FileExistsError:
        pass
    except OSError:                      # a filesystem without hard links
        try:
            fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o660)
            os.write(fd, header.encode())
            os.close(fd)
        except FileExistsError:
            pass
    finally:
        os.remove(tmp)


def append_locked(path, texts, deadline, header=None, marker=None):
    """Append each of `texts` with ONE write() apiece, under the mkdir lock, then fsync. The file
    is opened only once the lock is held, so on a network filesystem the append lands at the true
    end. With `marker`, nothing is written when the file already holds it (checked under the same
    lock), so one event is never logged twice. Returns how many texts were written."""
    with DirLock(path, deadline):
        if header:
            create_with_header(path, header)
        if marker:
            try:
                with open(path, "rb") as fh:
                    if marker.encode() in fh.read():
                        return 0
            except FileNotFoundError:
                pass
        fd = os.open(path, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o660)
        try:
            for t in texts:
                data = t.encode()
                n = os.write(fd, data)
                if n < len(data):            # a short write on a regular file: the disk is full
                    os.write(fd, data[n:])
            os.fsync(fd)
        finally:
            os.close(fd)
    return len(texts)


def ensure_readme(root, deadline):
    """Write the registry's README.md when it is missing or older than README_VERSION, so the
    directory never describes a layout it no longer has. Returns True when written."""
    path = os.path.join(root, "README.md")

    def current():
        try:
            with open(path) as fh:
                m = re.search(r"record_run\.py README v(\d+)", fh.read())
            return bool(m) and int(m.group(1)) >= README_VERSION
        except FileNotFoundError:
            return False
    if current():
        return False
    with DirLock(path, deadline):
        if current():
            return False
        write_text(path, README_TEXT)
    return True


def csv_row(values):
    """One activity-log row, quoted as csv would, on one line, short."""
    buf = io.StringIO()
    vals = [re.sub(r"\s+", " ", str(v if v is not None else "")).strip() for v in values]
    vals[-1] = vals[-1][:400]
    csv.writer(buf, lineterminator="\n").writerow(vals)
    return buf.getvalue()


# ---------------------------------------------------------------------------------- routing
def hive_exec_path():
    return os.environ.get("HIVE_EXEC") or os.path.join(HERE, "hive_exec.sh")


def hive_login():
    """The HIVE user when this machine can reach HIVE through hive_exec.sh, else None. The same
    two sources hive_exec.sh reads: the environment, then the saved hive.env at HIVE_ENV_FILE
    (read for these two names only; never copied anywhere)."""
    user, key = os.environ.get("HIVE_USER"), os.environ.get("HIVE_KEY")
    cfg = (os.environ.get("HIVE_ENV_FILE")
           or os.path.expanduser("~/.config/ucdavis-proteomics/hive.env"))
    if not (user and key) and os.path.isfile(cfg):
        try:
            with open(cfg) as fh:
                for ln in fh:
                    m = re.match(r"\s*(?:export\s+)?(HIVE_USER|HIVE_KEY)=(.*)$", ln)
                    if m and m.group(1) == "HIVE_USER":
                        user = user or m.group(2).strip().strip("'\"")
                    elif m:
                        key = key or m.group(2).strip().strip("'\"")
        except OSError:
            pass
    return user if (user and key and os.path.isfile(hive_exec_path())) else None


def route(remote_hop=False, read_only=False):
    """direct | ssh | not_core | none | off -- where a record from HERE can go."""
    if disabled():
        return "off"
    d = runs_dir().rstrip("/")
    parent = os.path.dirname(d)
    need = (os.R_OK | os.X_OK) if read_only else (os.W_OK | os.X_OK)
    if os.path.isdir(d) and os.access(d, need):
        return "direct"
    # Missing but creatable (the very first record) is a Core account too.
    if not read_only and not os.path.exists(d) and os.path.isdir(parent) \
            and os.access(parent, os.W_OK | os.X_OK):
        return "direct"
    if os.path.isdir(parent):
        return "not_core"          # on HIVE -- the group folder is visible -- but not writable
    if not remote_hop and hive_login():
        return "ssh"
    return "none"


REASONS = {"off": "disabled", "not_core": "not_core_member", "none": "not_on_hive"}


def where_text(r):
    d = runs_dir()
    return {
        "direct": f"direct: {d}",
        "ssh": f"ssh: {hive_login()}@hive:{d} (via hive_exec.sh)",
        "not_core": (f"not recorded: {d} is not writable by {whoami()} -- only Proteomics "
                     f"Core (proteomics-grp) runs are recorded"),
        "none": ("not recorded: not on HIVE, and no HIVE login is set up here "
                 "(hive.env, or HIVE_USER + HIVE_KEY)"),
        "off": "not recorded: RECORD_RUN=off",
    }[r]


def not_recorded(reason, detail, **extra):
    return dict({"recorded": False, "reason": reason, "detail": detail}, **extra)


# ------------------------------------------------------------------------ CoreOmics (PROT)
def parse_prot(values):
    """--prot values -> {"prot": "PROT_0807", "id": "<12 hex>"} (either may be None), or None.
    Accepts PROT_0807, prot-807, 0807, 807, #807 and a 12-hex CoreOmics id, several per value."""
    prot = hexid = None
    for v in values or []:
        for tok in re.split(r"[\s,;/:()]+", str(v)):
            t = tok.strip()
            if not t:
                continue
            if HEX_ID_RE.match(t.lower()):
                hexid = t.lower()
            else:
                m = PROT_RE.match(t)
                if m and int(m.group(1)) > 0:
                    prot = "PROT_%04d" % int(m.group(1))
    return {"prot": prot, "id": hexid} if (prot or hexid) else None


def prot_from_obj(obj, depth=0):
    """A PROT / CoreOmics id out of a metadata JSON -- ONLY from keys that name it (never from a
    path or a folder name: those are free-form). A `.core_submission.json` receipt carries
    internal_id (PROT_####) + id (12 hex)."""
    if not isinstance(obj, dict) or depth > 3:
        return None
    found = {}
    for k in ("prot", "prot_id", "prot_number", "internal_id", "coreomics_prot", "coreomics"):
        v = obj.get(k)
        if isinstance(v, (str, int)) and not isinstance(v, bool):
            p = parse_prot([v])
            if p and p["prot"] and (k != "internal_id" or str(v).upper().startswith("PROT")):
                found.setdefault("prot", p["prot"])
            if k == "coreomics" and p and p["id"]:
                found.setdefault("id", p["id"])
    for k in ("coreomics_id", "submission_id"):
        v = obj.get(k)
        if isinstance(v, str) and HEX_ID_RE.match(v.lower()):
            found.setdefault("id", v.lower())
    if str(obj.get("schema") or "").startswith("core_submission") \
            and isinstance(obj.get("id"), str) and HEX_ID_RE.match(obj["id"].lower()):
        found.setdefault("id", obj["id"].lower())
    for k in ("coreomics", "submission", "core_submission"):
        sub = prot_from_obj(obj.get(k), depth + 1) if isinstance(obj.get(k), dict) else None
        if isinstance(obj.get(k), dict):
            v = obj[k].get("id")
            if isinstance(v, str) and HEX_ID_RE.match(v.lower()):
                found.setdefault("id", v.lower())
        for kk, vv in (sub or {}).items():
            found.setdefault(kk, vv)
    return found or None


def find_prot(explicit, session, out):
    """{"prot", "id", "source"} or None. --prot first; then the metadata the session keeps."""
    p = parse_prot(explicit)
    if p:
        return dict(p, source="--prot")
    cands = []
    if session and os.path.isdir(session):
        cands += sorted(glob.glob(os.path.join(session, "input", "*.json")))
        cands += [os.path.join(session, f) for f in listdir(session) if f.endswith(".json")]
        d = session
        for _ in range(4):                # a core_submission receipt sits in the work dir above
            cands.append(os.path.join(d, ".core_submission.json"))
            d = os.path.dirname(d)
    if out and os.path.isdir(out):
        cands += [os.path.join(out, "fran_deposit.json"), os.path.join(out,
                                                                       "search_provenance.json")]
    for c in dict.fromkeys(cands):
        if (fsize(c) or 0) > 2 * MB:
            continue
        found = prot_from_obj(load_json(c))
        if found:
            return {"prot": found.get("prot"), "id": found.get("id"), "source": c}
    return None


def prot_label(p):
    if not p or not (p.get("prot") or p.get("id")):
        return "not recorded"
    return " / ".join(x for x in (p.get("prot"), f"`{p['id']}`" if p.get("id") else None) if x)


# --------------------------------------------------------------------- reading a search dir
def task_log(name):
    m = JOBLOG_RE.search(name)
    return bool(m and m.group(2))


def job_ids(out):
    """SLURM job ids of this search: jobs.txt (submit.sh writes one per line) and the ids in the
    log names -- a job submitted by hand has only the latter."""
    ids = []
    try:
        with open(os.path.join(out, "jobs.txt")) as fh:
            ids += re.findall(r"\d{5,}", fh.read())
    except OSError:
        pass
    for f in listdir(out):
        m = JOBLOG_RE.search(f)
        if m:
            ids.append(m.group(1))
    return list(dict.fromkeys(ids))


def sacct_jobs(ids, deadline):
    """([row], None) from ONE time-bounded sacct call -- a row per job and per array task -- or
    (None, why)."""
    if not ids:
        return None, "no SLURM job ids found (no jobs.txt, no <step>_<jobid>.log)"
    exe = os.environ.get("RECORD_RUN_SACCT")          # a stand-in, or "none" to skip it
    if exe == "none":
        return None, "sacct not consulted (RECORD_RUN_SACCT=none)"
    exe = exe or shutil.which("sacct") or next(
        (c for c in SACCT_CANDIDATES if os.access(c, os.X_OK)), None)
    if not exe:
        return None, "sacct is not available here"
    cols = ["id", "name", "state", "exit_code", "submit", "start", "end", "elapsed", "nodes",
            "partition", "account", "user"]
    try:
        p = subprocess.run([exe, "-j", ",".join(ids[:300]), "-X", "-n", "-P",
                            "--format=JobID,JobName,State,ExitCode,Submit,Start,End,Elapsed,"
                            "NodeList,Partition,Account,User"],
                           capture_output=True, text=True, stdin=subprocess.DEVNULL,
                           timeout=max(2, min(10, deadline.left() - 5)))
    except (OSError, subprocess.SubprocessError) as e:
        return None, f"sacct failed: {e}"
    rows = [dict(zip(cols, ln.split("|"))) for ln in p.stdout.splitlines()
            if ln.count("|") >= len(cols) - 1]
    for r in rows:
        r["state"] = (r.get("state") or "").split(" ")[0]     # "CANCELLED by 123" -> CANCELLED
    return (rows, None) if rows else (None, (p.stderr.strip() or "sacct returned no rows")[:200])


def _t(s):
    return s if s and s[:1].isdigit() else None


def summarize_jobs(rows):
    """One entry per job -- an array collapses to its task counts -- in submission order."""
    groups = {}
    for r in rows:
        base = r["id"].split("_")[0].split(".")[0]
        g = groups.setdefault(base, {"id": base, "name": r.get("name"), "tasks": 0, "states": {},
                                     "exit_codes": set(), "submit": None, "start": None,
                                     "end": None, "nodes": set(), "partition": r.get("partition"),
                                     "account": r.get("account"), "failed_tasks": []})
        g["tasks"] += 1
        g["states"][r["state"]] = g["states"].get(r["state"], 0) + 1
        g["exit_codes"].add(r.get("exit_code"))
        for k, pick in (("submit", min), ("start", min), ("end", max)):
            if _t(r.get(k)):
                g[k] = pick(x for x in (g[k], r[k]) if x)
        if r.get("nodes") and r["nodes"] != "None assigned":
            g["nodes"].add(r["nodes"])
        if r["state"] in BAD_STATES and len(g["failed_tasks"]) < 50:
            g["failed_tasks"].append({"id": r["id"], "state": r["state"],
                                      "exit_code": r.get("exit_code")})
    out = []
    for g in sorted(groups.values(), key=lambda g: (g["submit"] or "", g["id"])):
        g["exit_codes"] = sorted(x for x in g["exit_codes"] if x)
        nodes = sorted(g["nodes"])
        g["nodes"] = ", ".join(nodes[:3]) + (f" (+{len(nodes) - 3} more)" if len(nodes) > 3 else "")
        out.append(g)
    return out


def find_report(out, prov):
    """(path, bytes) of the quant of record. A report that is a DIRECTORY (Fulcrum) is not
    walked to size it: bytes is None."""
    res = prov.get("result") if isinstance(prov.get("result"), dict) else {}
    for c in [res.get("report")] + [os.path.join(out, c) for c in REPORT_CANDIDATES]:
        if c and os.path.exists(c):
            return c, (fsize(c) if os.path.isfile(c) else None)
    return None, 0


def input_files(out, prov):
    if isinstance(prov.get("files"), list) and prov["files"]:
        return [str(f) for f in prov["files"]]
    for name in ("job_input_files.txt", "parallel_input_files.txt", "file_list.txt"):
        try:
            with open(os.path.join(out, name)) as fh:
                files = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
            if files:
                return files
        except OSError:
            pass
    return []


def vendor_of(ext):
    return {".raw": "Thermo", ".d": "Bruker", ".wiff": "Sciex", ".wiff2": "Sciex",
            ".mzml": "converted mzML", ".dia": "DIA-NN .dia"}.get(ext.lower(), ext or "?")


def log_facts(path):
    """What an engine or job log says: the engine's version line, the mass accuracy it settled
    on, and the error lines. Read line by line, and bounded."""
    facts = {"version_line": None, "optimised_ms2_ppm": [], "recommended_ms1_ppm": [],
             "errors": []}
    try:
        with open(path, errors="replace") as fh:
            for i, ln in enumerate(fh):
                if i > 200000:
                    break
                s = ln.strip()
                if facts["version_line"] is None and re.match(r"(DIA-NN|Sage|FragPipe) \S+", s):
                    facts["version_line"] = s[:120]
                m = re.search(r"Optimised mass accuracy: ([\d.]+) ppm", s)
                if m:
                    facts["optimised_ms2_ppm"].append(float(m.group(1)))
                m = re.search(r"Recommended MS1 mass accuracy setting: ([\d.]+) ppm", s)
                if m:
                    facts["recommended_ms1_ppm"].append(float(m.group(1)))
                if re.search(r"\bERROR\b|error:|Killed|oom-kill|CANCELLED|DUE TO TIME LIMIT", s) \
                        and len(facts["errors"]) < 12:
                    facts["errors"].append(s[:300])
    except OSError:
        return None
    return facts


def stats_summary(report):
    """Median precursors / proteins per run from DIA-NN's stats file, read with
    check_report_runs.py -- the one reader of that file."""
    try:
        from check_report_runs import stats_path, stats_rows
    except Exception as e:
        return {"error": f"check_report_runs.py is not available here: {e}"}
    if not report:
        return {"error": "no report"}
    path = stats_path(report)
    if not os.path.isfile(path) and report.endswith(".tsv"):
        path = report[:-len(".tsv")] + ".stats.tsv"      # DIA-NN 1.x: report.tsv
    rows = stats_rows(path)
    if not rows:
        return {"stats_file": None, "error": "no report.stats.tsv beside the report (DIA-NN "
                                             "writes one; the other engines do not)"}

    def num(r, c):
        try:
            return float(r.get(c) or 0)
        except ValueError:
            return 0.0
    per = [{"run": k, "precursors": int(num(r, "Precursors.Identified")),
            "proteins": int(num(r, "Proteins.Identified"))} for k, r in rows.items()]
    prec, prot = [p["precursors"] for p in per], [p["proteins"] for p in per]
    return {"stats_file": path, "n_runs": len(per),
            "median_precursors": statistics.median(prec),
            "median_proteins": statistics.median(prot),
            "precursors_range": [min(prec), max(prec)], "proteins_range": [min(prot), max(prot)],
            "zero_id_runs": [p["run"] for p in per if p["precursors"] == 0],
            "per_run": per}


def parameters(out, prov, params_file, manifest_path, rationale):
    """The settings that decide a search, each with where it came from. The values are read by
    make_methods.search_record() -- THE reader of a search's parameters, which the Methods and
    the SDRF use as well, so the three cannot disagree (DE-LIMP rule 3). The rationale sidecar
    adds WHY each value is what it is."""
    p = {"record": None, "why": {}, "resolution": None, "mass_accuracy": {}, "scan_window": None}
    try:
        from make_methods import search_record
        prov_path = os.path.join(out, "search_provenance.json")
        p["record"] = search_record(params_file,
                                    prov_path if os.path.isfile(prov_path) else None,
                                    manifest_path)
    except Exception as e:
        p["record_error"] = f"not parsed here (make_methods.py: {type(e).__name__}: {e})"
    ra = (rationale or {}).get("rationale") or {}
    for flag in ("--qvalue", "--min-pr-mz", "--cut", "--missed-cleavages", "--mass-acc",
                 "--mass-acc-ms1", "--window", "--cont-quant-exclude", "enzyme"):
        if isinstance(ra.get(flag), dict) and ra[flag].get("source"):
            p["why"][flag] = ra[flag]["source"]
    man = load_json(manifest_path) or {}
    p["resolution"] = ((rationale or {}).get("resolution") or man.get("resolution")
                       or (man.get("search") or {}).get("resolution"))
    p["instrument_class"] = (rationale or {}).get("instrument_class") or man.get(
        "instrument_class")
    p["instrument_class_label"] = ((rationale or {}).get("class_label")
                                   or man.get("instrument_label"))
    ma = p["mass_accuracy"]
    ma["source"] = ((rationale or {}).get("mass_accuracy_source")
                    or (man.get("search") or {}).get("ppm_source"))
    ma["plan"] = (rationale or {}).get("mass_accuracy_plan")
    res = prov.get("result") if isinstance(prov.get("result"), dict) else {}
    if res.get("mass_acc"):
        ma["measured"] = res["mass_acc"]
    try:
        with open(os.path.join(out, "massacc.txt")) as fh:
            ma["massacc_txt"] = fh.read().strip()[:200]
    except OSError:
        pass
    p["scan_window"] = prov.get("scan_window")
    return p


def fasta_facts(meta_path):
    m = load_json(meta_path)
    if not isinstance(m, dict):
        return None
    sel = m.get("selected") if isinstance(m.get("selected"), dict) else m
    g = lambda k: sel.get(k) if sel.get(k) is not None else m.get(k)  # noqa: E731
    sf = g("staged_file") or {}
    return {"meta_file": meta_path, "fasta": g("fasta"), "organism": g("organism"),
            "taxid": g("taxid"), "organism_source": g("organism_source"),
            "proteome": g("proteome"), "proteome_type": g("proteome_type"),
            "source": g("source"), "n_proteome": g("n_proteome"),
            "n_sequences": g("n_sequences") or g("n_entries"), "md5": g("md5"),
            "content_requested": g("content_requested"), "content_used": g("content_used"),
            "uniprot_release": g("uniprot_release") or None,
            "staged_file_date": (sf.get("mtime_utc") or "")[:10] or None,
            "contaminant_set": g("contaminant_set"),
            "n_contaminants_appended": g("n_contaminants_appended"),
            "n_contaminants_already_present": g("n_contaminants_already_present"),
            "n_contaminants_dropped_as_target": g("n_contaminants_dropped_as_target"),
            "contaminants_dropped_note": g("contaminants_dropped_note"),
            "contaminant_target_rule": g("contaminant_target_rule"),
            "contaminant_source": g("contaminant_source"),
            "contaminant_citation": g("contaminant_citation"),
            "cont_quant_exclude": g("diann_cont_quant_exclude"),
            "digestion_enzymes_used": g("digestion_enzymes_used"),
            "warnings": g("warnings") or []}


def first_existing(paths):
    for p in paths:
        if p and os.path.isfile(p):
            return os.path.abspath(p)
    return None


def fran_facts(out):
    try:
        from fran_deposit import read_receipt
        r = read_receipt(out)
    except Exception:
        r = load_json(os.path.join(out, "fran_deposit.json"))
    if not r:
        return {"status": None, "detail": "no fran_deposit.json in the search folder -- not "
                                          "staged for FRAN (yet)"}
    return {"status": r.get("status"), "entry": r.get("entry"), "engine": r.get("engine"),
            "organism": r.get("organism"), "taxon": r.get("taxon"),
            "staged_by": r.get("staged_by"), "n_linked": len(r.get("linked") or []),
            "xic": (r.get("xic") or {}).get("n_files") if isinstance(r.get("xic"), dict) else None,
            "verified": r.get("verified"), "link_errors": r.get("link_errors"),
            "search_id": r.get("search_id")}


def detection_facts(paths):
    """The newest detect_acquisition.py output among `paths`: per-file warnings and the cohort
    flags (resolution unknown / mixed, ion-trap MS2, low confidence, mixed m/z range)."""
    best = None
    for p in paths:
        d = load_json(p) if (fsize(p) or 0) < 5 * MB else None
        if isinstance(d, dict) and isinstance(d.get("files"), list) and "overall" in d:
            if best is None or (os.path.getmtime(p) >= os.path.getmtime(best[0])):
                best = (p, d)
    if not best:
        return None
    p, d = best
    warn = {}
    for f in d["files"]:
        for w in (f.get("warnings") or []) if isinstance(f, dict) else []:
            warn.setdefault(str(w)[:300], []).append(os.path.basename(str(f.get("file", "?"))))
    return {"file": p, "overall": d.get("overall"), "instrument": d.get("instrument"),
            "needs_confirmation": d.get("needs_confirmation"),
            "low_confidence_files": d.get("low_confidence_files") or [],
            "precursor_mz_range_mixed": d.get("precursor_mz_range_mixed"),
            "orbitrap_resolution_unknown": d.get("orbitrap_resolution_unknown"),
            "resolution_mixed": d.get("resolution_mixed") or [],
            "ms2_ion_trap": d.get("ms2_ion_trap"),
            "file_warnings": [{"warning": w, "files": fs} for w, fs in warn.items()][:20]}


def detection_candidates(dirs):
    out = []
    for d in dirs:
        out += [os.path.join(d, f) for f in listdir(d)
                if f.endswith(".json") and re.search(r"detect|acquisition|acq", f, re.I)]
    return out


def label_for(out, session, explicit):
    """The short description part of the folder name."""
    if explicit:
        base = explicit
    elif session:
        base = os.path.basename(str(session).rstrip("/"))
    else:
        d = os.path.abspath(out).rstrip("/")
        base = os.path.basename(d)
        while base.lower() in GENERIC_DIR_NAMES and os.path.dirname(d) not in ("/", d):
            d = os.path.dirname(d)
            base = os.path.basename(d)
    base = re.sub(r"^\d{4}-\d{2}-\d{2}_", "", base)
    return clean(base, 70) or "search"


def slurm_context(a):
    """The job this call runs in: the job-end hook runs inside the job that ended, and for a
    failure that job IS the failing step (DataAnalysis: 'with the failing step in the log')."""
    jid = os.environ.get("SLURM_JOB_ID")
    if not jid and not getattr(a, "step", None):
        return None
    return {"job_id": jid, "job_name": os.environ.get("SLURM_JOB_NAME"),
            "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
            "array_task": os.environ.get("SLURM_ARRAY_TASK_ID"),
            "step": getattr(a, "step", None) or os.environ.get("SLURM_JOB_NAME"),
            "node": os.environ.get("SLURMD_NODENAME")}


# ------------------------------------------------------------------------------ the plan
def new_plan(event, a):
    return {"event": event, "out": None, "session": None,
            "identity": {"out": None, "session": None}, "user": None, "label": None,
            "date": None, "search": None, "analysis": None, "prot": None,
            "copies": [], "not_copied": [], "generated": {}, "zip": None, "links": {},
            "findings": [], "issue_tags": list(a.issues_tag or []),
            "skill_version": a.skill_version or skill_version(),
            "by": whoami(), "host": socket.gethostname(),
            "file_cap": int(a.file_cap_mb * MB), "zip_cap": int(a.zip_cap_gb * GB),
            "budget": int(env_num("RECORD_RUN_COPY_BUDGET_MB", 200) * MB)}


def add_copy(plan, src, rel, section, cap=None):
    """Queue `src` to be copied to `rel` in the record, or record why not. Deduplicated by real
    path; a second file wanting the same name is kept apart by its folder's name."""
    if not src or not os.path.isfile(src):
        return
    real = os.path.realpath(src)
    if any(c["real"] == real for c in plan["copies"] + plan["not_copied"]):
        return
    if any(c["rel"] == rel for c in plan["copies"]):
        d, n = os.path.split(rel)
        rel = f"{d}/{clean(os.path.basename(os.path.dirname(real)), 30)}__{n}"
    size = fsize(src)
    cap = min(cap or plan["file_cap"], plan["file_cap"])
    used = sum(c["bytes"] or 0 for c in plan["copies"])
    if size is None:
        why = "unreadable"
    elif size > cap:
        why = f"{fmt_bytes(size)} is over the {fmt_bytes(cap)} per-file cap"
    elif used + size > plan["budget"]:
        why = (f"over the per-record copy budget ({fmt_bytes(plan['budget'])}, "
               f"RECORD_RUN_COPY_BUDGET_MB)")
    else:
        why = secret_reason(src)
    entry = {"src": os.path.abspath(src), "real": real, "rel": rel, "bytes": size,
             "section": section}
    if why:
        entry["reason"] = why
        plan["not_copied"].append(entry)
    else:
        plan["copies"].append(entry)


def plan_search(plan, out, a, deadline):
    """Everything about the search that is readable from here, plus what to copy and link."""
    out = os.path.abspath(out)
    prov = load_json(os.path.join(out, "search_provenance.json")) or {}
    plan["out"] = out
    plan["identity"]["out"] = os.path.realpath(out)
    # Who ran it: the owner of the search folder (on HIVE that is the HIVE account, whoever is
    # recording it -- a maintainer's backfill must not re-attribute someone's search).
    plan["user"] = owner(out) or plan["user"] or whoami()
    s = {"out_dir": out, "out_dir_realpath": os.path.realpath(out),
         "gathered_at": now_iso(), "gathered_on": socket.gethostname(), "warnings": []}
    if not prov:
        s["warnings"].append("no search_provenance.json -- this search did not run through "
                             "run_search.py, or predates it; engine and version are inferred")

    # ---- engine + the version that ran
    engine, version, vsrc = (prov.get("engine") or "").lower() or None, prov.get("version"), None
    if engine:
        ev = prov.get("engine_version") if isinstance(prov.get("engine_version"), dict) else {}
        vsrc = ev.get("source") or "search_provenance.json"
    else:
        try:
            from fran_deposit import detect_engine
            engine, version, vsrc = detect_engine(out)
        except Exception as e:
            s["warnings"].append(f"engine not determined (fran_deposit.py: {e})")
    report, report_bytes = find_report(out, prov)
    logp = first_existing([os.path.splitext(report)[0] + ".log.txt" if report else None,
                           os.path.join(out, "report.log.txt")])
    lf = log_facts(logp) if logp else None
    if lf and lf["version_line"]:
        m = re.match(r"\S+ (\S+)", lf["version_line"])
        if m and version and m.group(1) != str(version):
            s["warnings"].append(f"search_provenance.json says {engine} {version} but the engine "
                                 f"log says '{lf['version_line']}'")
        if m and not version:
            version, vsrc = m.group(1), f"engine log ({os.path.basename(logp)})"
    try:
        from make_methods import ENGINE_LABEL
        label = ENGINE_LABEL.get(engine, engine)
    except Exception:
        label = engine
    s.update(engine=engine, engine_label=label, engine_version=version,
             engine_version_source=vsrc, engine_log=logp,
             engine_log_version_line=(lf or {}).get("version_line"),
             engine_log_mass_accuracy={k: (lf or {}).get(k) or [] for k in (
                 "optimised_ms2_ppm", "recommended_ms1_ppm")},
             engine_command=prov.get("resolved_command"), search_mode=prov.get("search_mode"),
             routing_reason=prov.get("parallel_routing_reason"))

    # ---- when, and did it work
    rows, why = sacct_jobs(job_ids(out), deadline)
    jobs = summarize_jobs(rows) if rows else []
    s["jobs"], s["jobs_note"] = jobs, why
    bad = [g for g in jobs if set(g["states"]) & BAD_STATES]
    live = [g for g in jobs if set(g["states"]) & LIVE_STATES]
    times = {"source": None}
    if jobs:
        times.update(submitted=min((g["submit"] for g in jobs if g["submit"]), default=None),
                     started=min((g["start"] for g in jobs if g["start"]), default=None),
                     finished=None if live else max((g["end"] for g in jobs if g["end"]),
                                                    default=None),
                     source=f"sacct ({len(jobs)} job(s))")
    if not times.get("submitted"):
        times["submitted"] = mtime_iso(os.path.join(out, "search_provenance.json"))
        times["source"] = ("file dates (submitted = search_provenance.json written, "
                           "finished = report written)")
    if not times.get("finished") and report and not live:
        times["finished"] = mtime_iso(report)
    has_report = bool(report and (report_bytes is None or report_bytes > 0))
    # A failed first attempt followed by a newer report is a search that completed.
    report_new = has_report and all(same_or_later(mtime_iso(report), g["end"]) for g in bad)
    reporter = slurm_context(a)
    if reporter:
        s["reported_by"] = reporter
    if a.status:
        status, ssrc = a.status, ("--status, from the job that ran it"
                                  + (f" ({reporter['job_name']} {reporter['job_id']})"
                                     if reporter and reporter.get("job_id") else ""))
        if live and not times.get("finished"):
            times["finished"] = now_iso()
    elif bad and not (has_report and report_new):
        status, ssrc = "failed", "sacct: " + ", ".join(
            f"{g['name'] or g['id']} {'/'.join(sorted(set(g['states']) & BAD_STATES))}"
            for g in bad)
    elif live:
        status, ssrc = "running", "sacct: " + ", ".join(
            f"{g['name'] or g['id']} {'/'.join(sorted(set(g['states']) & LIVE_STATES))}"
            for g in live)
    elif has_report:
        status = "completed"
        ssrc = ("sacct: every job COMPLETED, and the report is present" if jobs and not bad else
                "a non-empty report is present, newer than the failed attempt(s) above"
                if bad else "a non-empty report is present")
    elif jobs and all(set(g["states"]) == {"COMPLETED"} for g in jobs):
        status, ssrc = "failed", ("sacct says COMPLETED but there is no report -- DIA-NN exits 0 "
                                  "on fatal errors")
    else:
        status, ssrc = "unknown", why or "no report and no job state"
    if status == "completed" and not has_report:
        s["warnings"].append("recorded as completed, but no non-empty report was found")
    exit_code, esrc = a.exit_code, ("--exit-code" if a.exit_code is not None else None)
    if exit_code is None and jobs and jobs[-1]["exit_codes"]:
        codes = jobs[-1]["exit_codes"]
        exit_code = codes[0] if len(codes) == 1 else ",".join(codes)
        # sacct writes <exit>:<signal>; "0:0" is exit 0. A signal is kept as sacct wrote it.
        m = re.match(r"^(\d+):0$", str(exit_code))
        exit_code = int(m.group(1)) if m else exit_code
        esrc = f"sacct ExitCode of {jobs[-1]['name'] or jobs[-1]['id']}"
    failing = None
    if status == "failed":
        failing = (reporter or {}).get("step") or next(
            (g["name"] or g["id"] for g in bad), None)
    s.update(status=status, status_source=ssrc, exit_code=exit_code, exit_code_source=esrc,
             failing_step=failing, times=times)
    if status in ("failed", "unknown"):
        errs = list((lf or {}).get("errors") or [])
        logs = sorted((f for f in listdir(out) if JOBLOG_RE.search(f)),
                      key=lambda f: os.path.getmtime(os.path.join(out, f)))
        for f in logs[-5:]:
            errs += (log_facts(os.path.join(out, f)) or {}).get("errors") or []
        if errs:
            s["errors_excerpt"] = errs[-12:]

    # ---- the data
    files = input_files(out, prov)
    exts = {}
    for f in files:
        e = os.path.splitext(f.rstrip("/"))[1].lower() or "?"
        exts[e] = exts.get(e, 0) + 1
    dirs = sorted({os.path.dirname(f.rstrip("/")) for f in files})
    try:
        common = os.path.commonpath(dirs) if dirs else None
    except ValueError:
        common = None
    params_file = prov.get("params_file")
    if params_file and not os.path.isabs(params_file):
        # run_search.py records the path as given; relative means relative to where it ran,
        # which is the folder holding the out dir in every documented layout.
        params_file = first_existing([os.path.join(os.path.dirname(out), params_file),
                                      os.path.join(out, params_file)]) or params_file
    params_file = params_file or first_existing([os.path.join(out, f) for f in PARAMS_FALLBACK])
    rationale_path = first_existing([f"{params_file}.rationale.json" if params_file else None])
    rationale = load_json(rationale_path) or {}
    sess = plan.get("session") if plan.get("session") and os.path.isdir(plan["session"]) else None
    manifest_path = first_existing(
        [os.path.join(os.path.dirname(params_file), "workflow.manifest.json") if params_file
         else None,
         os.path.join(sess, "input", "wf", "workflow.manifest.json") if sess else None,
         os.path.join(os.path.dirname(out), "wf", "workflow.manifest.json"),
         os.path.join(os.path.dirname(os.path.dirname(out)), "input", "wf",
                      "workflow.manifest.json")])
    man = load_json(manifest_path) or {}
    instr = rationale.get("instrument") or ", ".join(man.get("instruments") or []) or None
    s["data"] = {"n_files": prov.get("n_files") or len(files), "extensions": exts,
                 "vendor": ", ".join(sorted({vendor_of(e) for e in exts})) or None,
                 "instrument": instr,
                 "instrument_source": (os.path.basename(rationale_path) if rationale.get(
                     "instrument") else "workflow.manifest.json" if instr else None),
                 "acquisition": rationale.get("acquisition") or man.get("acquisition"),
                 "raw_dirs": dirs[:50], "n_raw_dirs": len(dirs), "raw_common_dir": common,
                 "files": files}
    det_dirs = [os.path.dirname(out)] + ([sess, os.path.join(sess, "input"),
                                          os.path.join(sess, "output")] if sess else [])
    if getattr(a, "detect_json", None):
        s["detection"] = detection_facts([a.detect_json])
    else:
        s["detection"] = detection_facts(detection_candidates(det_dirs))

    # ---- parameters + the sequence database
    s["parameters"] = parameters(out, prov, params_file, manifest_path, rationale)
    fasta = prov.get("fasta")
    meta_path = first_existing(
        [f"{fasta}.meta.json" if fasta else None]
        + sorted(glob.glob(os.path.join(out, "*.fasta.meta.json")))
        + sorted(glob.glob(os.path.join(os.path.dirname(out), "*.fasta.meta.json")))
        + (sorted(glob.glob(os.path.join(sess, "input", "*.fasta.meta.json"))) if sess else [])
        + sorted(glob.glob(os.path.join(os.path.dirname(os.path.dirname(out)), "input",
                                        "*.fasta.meta.json"))))
    s["fasta"] = fasta_facts(meta_path) if meta_path else None
    if not s["fasta"]:
        s["warnings"].append("no FASTA sidecar (<fasta>.meta.json) was found -- the organism "
                             "and database are not recorded")

    # ---- results + where the big things are (listed and counted, never copied or walked)
    s["results"] = stats_summary(report)
    quant = []
    for q in QUANT_DIRS:
        d = os.path.join(out, q)
        fs = [f for f in listdir(d) if f.endswith(".quant")]
        if fs:
            quant.append({"dir": d, "n_files": len(fs), "bytes": sum(
                fsize(os.path.join(d, f)) or 0 for f in fs) if len(fs) <= 200 else None})
    xic = []
    rx = os.path.join(out, "report_xic")
    if os.path.isdir(rx):
        xic.append({"dir": rx, "n_files": sum(f.endswith(".xic.parquet") for f in listdir(rx))})
    xr = os.path.join(out, "xic")
    if os.path.isdir(xr):
        xic.append({"dir": xr, "n_task_dirs": sum(f.endswith("_xic") for f in listdir(xr))})
    libs = [{"file": os.path.join(out, f), "bytes": fsize(os.path.join(out, f))}
            for f in listdir(out) if f.endswith((".speclib", "-lib.parquet"))]
    s["outputs"] = {"report": report, "report_bytes": report_bytes, "quant": quant, "xic": xic,
                    "libraries": libs, "n_task_logs": sum(task_log(f) for f in listdir(out))}
    s["fran"] = fran_facts(out)
    plan["search"] = s

    # ---- copies: search products to output/search/, the search's inputs to input/
    for f in ("search_provenance.json", "report.manifest.txt", "fran_deposit.json", "jobs.txt",
              "submit.sh", "window.json", "window.txt", "massacc.txt", "job_input_files.txt",
              "parallel_input_files.txt", "file_list.txt", "results.json", "fragpipe.workflow"):
        add_copy(plan, os.path.join(out, f), f"output/search/{f}", "search")
    for f in listdir(out):
        if f.endswith((".log.txt", ".stats.tsv", ".cfg", ".rationale.json", ".fp-manifest",
                       ".toml")) or re.match(r"log_.*\.txt$", f):
            add_copy(plan, os.path.join(out, f), f"output/search/{f}", "search")
    logs = sorted((f for f in listdir(out) if JOBLOG_RE.search(f)),
                  key=lambda f: os.path.getmtime(os.path.join(out, f)))
    # One log per job always; per-file array-task logs only when the search failed (the newest
    # ten) -- a 399-file chain writes 800 of them, identical apart from the file name.
    for f in [f for f in logs if not task_log(f)] + (
            [f for f in logs if task_log(f)][-10:] if status == "failed" else []):
        add_copy(plan, os.path.join(out, f), f"output/search/{f}", "search")
    for f in (params_file, prov.get("resolved_params_file"), rationale_path, meta_path,
              manifest_path):
        if f:
            add_copy(plan, f, f"input/{os.path.basename(f)}", "search")
    if files:
        plan["generated"].setdefault("input/raw_files.txt", (
            "# Raw MS files this search read (never copied -- too large). From "
            "search_provenance.json.\n" + "\n".join(files) + "\n"))
    if report:
        plan["links"][f"output/search/{os.path.basename(report.rstrip('/'))}"] = report
    plan["links"]["output/search/search_out"] = out
    plan["issue_tags"] = list(dict.fromkeys(plan["issue_tags"] + [label_for(out, None, None)]))
    plan["date"] = plan["date"] or (times.get("submitted") or now_iso())[:10]


def zip_survey(path):
    with zipfile.ZipFile(path) as z:
        infos = z.infolist()
    quant = [i for i in infos if i.filename.lower().endswith(".quant")]
    secret = [i for i in infos if not i.is_dir()
              and SECRET_NAME_RE.search(os.path.basename(i.filename.rstrip("/")))]
    spec = [i for i in infos if i.filename.endswith(".speclib")]
    xic = [i for i in infos if "_xic/" in i.filename or i.filename.endswith(".xic.parquet")]
    dropped = {id(i) for i in quant + secret}
    grp = lambda xs: {"n": len(xs), "bytes_compressed": sum(i.compress_size for i in xs),  # noqa
                      "bytes": sum(i.file_size for i in xs)}
    return {"n_entries": len(infos), "bytes": fsize(path),
            "quant": dict(grp(quant), dirs=sorted({os.path.dirname(i.filename) for i in quant})),
            "secrets": [i.filename for i in secret], "speclib": grp(spec), "xic": grp(xic),
            "drop": [i.filename for i in quant + secret],
            "kept_bytes": sum(i.compress_size for i in infos if id(i) not in dropped)}


def session_search_out(session, explicit):
    """The search out dir behind a session: --out; else where the session's search_provenance.json
    says the report was written, when that is somewhere else (a hive_remote session holds a COPY
    pulled from HIVE); else output/search when the search ran there; else .recovery.json."""
    if explicit:
        return explicit
    ss = os.path.join(session, "output", "search")
    prov = load_json(os.path.join(ss, "search_provenance.json")) or {}
    rep = (prov.get("result") or {}).get("report") if isinstance(prov.get("result"), dict) \
        else None
    if rep and os.path.isabs(rep) and \
            os.path.realpath(os.path.dirname(rep)) != os.path.realpath(ss):
        return os.path.dirname(rep)
    if prov or any(os.path.exists(os.path.join(ss, c)) for c in REPORT_CANDIDATES):
        return ss
    rec = load_json(os.path.join(session, ".recovery.json")) or {}
    return os.path.dirname(rec["report"]) if rec.get("report") else None


def section_of(md_path, title):
    """The body of `## <title>` in a Markdown file (up to the next `## `), or None."""
    try:
        with open(md_path, errors="replace") as fh:
            text = fh.read()
    except OSError:
        return None
    m = re.search(r"^##\s+" + re.escape(title) + r"\s*$(.*?)(?=^##\s|\Z)", text, re.M | re.S)
    return m.group(1).strip() if m and m.group(1).strip() else None


def plan_analysis(plan, session, a, zip_cap):
    session = os.path.abspath(session)
    plan["session"] = session
    plan["identity"]["session"] = os.path.realpath(session)
    out_d = os.path.join(session, "output")
    de = load_json(os.path.join(out_d, "tables", "de_provenance.json")) or {}
    man_txt = os.path.join(session, "MANIFEST.txt")
    ok, skipped = 0, []
    try:
        with open(man_txt, errors="replace") as fh:
            for ln in fh:
                if ln.startswith("[OK]"):
                    ok += 1
                elif ln.startswith("[SKIPPED]"):
                    skipped.append(ln.strip()[:300])
    except OSError:
        man_txt = None
    dep = os.path.join(out_d, "DATA_SUBMISSION")
    docx = [os.path.join(out_d, f) for f in listdir(out_d) if f.lower().endswith(".docx")]
    audit = load_json(first_existing([os.path.join(out_d, "AUDIT.json"),
                                      os.path.join(session, "AUDIT.json")])) or {}
    sq = load_json(first_existing([os.path.join(out_d, "SAMPLE_QUALITY.json"),
                                   os.path.join(session, "SAMPLE_QUALITY.json")])) or {}
    expert = None
    for md in (os.path.join(session, "README.md"), os.path.join(out_d, "AI_Analysis_Report.md")):
        expert = expert or section_of(md, "Expert Review Notes")
    an = {"session": session, "finalized": mtime_iso(man_txt) if man_txt else None,
          "recorded_at": now_iso(),
          "de": {k: de.get(k) for k in ("method", "contrasts", "significant_per_contrast",
                                         "q_cutoff", "logfc", "adjp")} if de else None,
          "manifest": {"file": man_txt, "n_ok": ok, "skipped": skipped} if man_txt else None,
          "docx": [{"file": d, "kind": "Methods (Word)" if "method" in os.path.basename(d).lower()
                    else "Report (Word)"} for d in docx],
          "methods_md": first_existing([os.path.join(out_d, "methods.md"),
                                        os.path.join(out_d, "METHODS.md")]),
          "data_submission": dep if os.path.isdir(dep) else None,
          "how_to_submit": first_existing([os.path.join(dep, "HOW_TO_SUBMIT.md")]),
          "reproduce_md": first_existing([os.path.join(out_d, "reproducibility", "REPRODUCE.md")]),
          "audit": [{"check": f.get("check"), "status": f.get("status"),
                     "message": f.get("message")} for f in (audit.get("findings") or [])
                    if isinstance(f, dict) and f.get("status") in ("WARN", "FAIL")],
          "audit_overall": audit.get("overall"),
          "sample_quality_flags": [str(x) for x in (sq.get("flags") or [])][:30],
          "expert_review_notes": expert[:20000] if expert else None,
          "warnings": []}
    if not de:
        an["warnings"].append("no output/tables/de_provenance.json -- the DE step has not run, "
                              "or wrote somewhere else")
    if not man_txt:
        an["warnings"].append("no MANIFEST.txt at the session root -- session.py finalize has "
                              "not run on this session")

    # ---- copies, mirroring the session: top level, input/, output/, scripts/
    for rel in ("README.md", "MANIFEST.txt", "DIFFERENCES.md"):
        add_copy(plan, os.path.join(session, rel), rel, "analysis")
    inp = os.path.join(session, "input")
    for f in listdir(inp):
        if f in ("conditions.csv", "raw_files.txt") or f.endswith(
                (".meta.json", ".rationale.json", ".cfg")) or f.startswith("params."):
            add_copy(plan, os.path.join(inp, f), f"input/{f}", "analysis")
    for f in listdir(os.path.join(inp, "wf")):
        if f.endswith((".json", ".cfg")) or f.startswith("params."):
            add_copy(plan, os.path.join(inp, "wf", f), f"input/{f}", "analysis")
    for d in docx:                                     # the deliverable of record, first
        add_copy(plan, d, f"output/{os.path.basename(d)}", "analysis")
    for f in ("methods.md", "METHODS.md", "AI_Analysis_Report.md", "OUTPUT_FILES.md", "AUDIT.md",
              "AUDIT.json", "SAMPLE_QUALITY.md", "SAMPLE_QUALITY.json", "QC_Report.html"):
        add_copy(plan, os.path.join(out_d, f), f"output/{f}", "analysis")
    for f in listdir(os.path.join(out_d, "tables")):
        add_copy(plan, os.path.join(out_d, "tables", f), f"output/tables/{f}", "analysis")
    for f in listdir(os.path.join(out_d, "figures")):
        add_copy(plan, os.path.join(out_d, "figures", f), f"output/figures/{f}", "analysis",
                 cap=5 * MB)
    for f in ("HOW_TO_SUBMIT.md", "sdrf.tsv"):
        add_copy(plan, os.path.join(dep, f), f"output/DATA_SUBMISSION/{f}", "analysis")
    add_copy(plan, os.path.join(session, "logs", "commands.log"), "scripts/commands.log",
             "analysis")
    for f in ("reproduce.sh", "REPRODUCE.md", "run_manifest.json"):
        add_copy(plan, os.path.join(out_d, "reproducibility", f), f"scripts/{f}", "analysis")

    # ---- the session zip: "save the output" -- minus per-run .quant and anything secret-named
    zpath = os.path.abspath(a.zip) if a.zip else session.rstrip("/") + ".zip"
    z = {"path": zpath, "exists": os.path.isfile(zpath), "copy": False}
    if not z["exists"]:
        z["reason"] = "no session zip -- finalize ran without --zip, or the zip was moved"
    else:
        try:
            sv = zip_survey(zpath)
            z.update(sv)
            if sv["quant"]["n"]:
                plan["findings"].append({
                    "id": "session_zip_contains_quant", "zip": zpath,
                    "detail": (f"the session zip holds {sv['quant']['n']} per-run .quant files "
                               f"({fmt_bytes(sv['quant']['bytes_compressed'])} compressed, "
                               f"{fmt_bytes(sv['quant']['bytes'])} unpacked) under "
                               f"{', '.join(sv['quant']['dirs'][:4])}. session.py finalize --zip "
                               f"leaves out only output/raw_data and upload_staging, so a search "
                               f"run into output/search carries its .quant into the zip. They "
                               f"were left out of the registry's copy.")})
            if sv["kept_bytes"] > zip_cap:
                z["reason"] = (f"{fmt_bytes(sv['kept_bytes'])} (without .quant) is over the "
                               f"{fmt_bytes(zip_cap)} cap -- the record points at it instead")
            else:
                z.update(copy=True, dest_rel=os.path.basename(zpath))
                if sv["drop"]:
                    z["note"] = (f"copied without {sv['quant']['n']} .quant and "
                                 f"{len(sv['secrets'])} secret-named entr"
                                 f"{'y' if len(sv['secrets']) == 1 else 'ies'}")
        except (zipfile.BadZipFile, OSError, ValueError) as e:
            z["reason"] = f"could not read the zip: {e}"
    an["zip"] = {k: v for k, v in z.items() if k != "drop"}
    plan["zip"] = z
    plan["analysis"] = an
    plan["date"] = plan["date"] or (re.match(r"\d{4}-\d{2}-\d{2}",
                                             os.path.basename(session)) or [None])[0]


# ------------------------------------------------------------------------ zip without .quant
def _strip_zip64_extra(extra):
    out, i = b"", 0
    while i + 4 <= len(extra):
        hid, ln = struct.unpack("<HH", extra[i:i + 4])
        if hid != 1:
            out += extra[i:i + 4 + ln]
        i += 4 + ln
    return out


def copy_zip_without(src, dst, drop, deadline):
    """Copy a zip leaving out `drop`, moving each kept member's COMPRESSED bytes as they are.
    Re-compressing a multi-GB zip would be real CPU on a login node; this is only I/O. The result
    is re-opened and its member count checked before it replaces anything."""
    drop, kept = set(drop), 0
    part = f"{dst}.part{os.getpid()}"
    os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
    try:
        with zipfile.ZipFile(src) as zin, open(src, "rb") as raw, \
                zipfile.ZipFile(part, "w", allowZip64=True) as zout:
            fp = zout.fp
            for zi in zin.infolist():
                if zi.filename in drop:
                    continue
                raw.seek(zi.header_offset)
                head = raw.read(30)
                if len(head) < 30 or head[:4] != b"PK\x03\x04":
                    raise ValueError(f"no local header for {zi.filename}")
                n, m = struct.unpack("<HH", head[26:30])
                new = copy.copy(zi)
                new.header_offset = fp.tell()
                new.extra = _strip_zip64_extra(zi.extra)     # zipfile re-adds it where needed
                fp.write(head + raw.read(n + m))
                copy_n(raw, fp, zi.compress_size, deadline)
                if zi.flag_bits & 0x08:                      # a data descriptor follows
                    z64 = max(zi.compress_size, zi.file_size) >= 0xFFFFFFFF
                    sig = raw.read(4)
                    fp.write(sig + raw.read((16 if z64 else 8)
                                            + (4 if sig == b"PK\x07\x08" else 0)))
                zout.filelist.append(new)
                zout.NameToInfo[new.filename] = new
                kept += 1
            zout.start_dir = fp.tell()
        with zipfile.ZipFile(part) as chk:
            if len(chk.infolist()) != kept:
                raise ValueError("the filtered zip does not list every kept member")
        os.replace(part, dst)
    finally:
        if os.path.exists(part):
            os.remove(part)
    return kept


# ------------------------------------------------------------------------- skill issues
ISSUE_FILE_RE = re.compile(r"^(\d{4}-\d{2}-\d{2})_(.+)\.md$")


def find_issues(user, tags, days):
    """report_issue.sh files for this run: <date>_<user>_<tag>.md whose tag matches the session
    (compared on letters and digits only -- the agent may have spelled the name either way), from
    any day, since setup problems are often recorded the day before the search; plus the same
    user's untagged files from the days the search ran."""
    d = issues_dir()
    names = listdir(d)
    r = {"dir": d, "files": []}
    if not names and not os.path.isdir(d):
        r["note"] = "the issue folder is not readable from here"
        return r
    cu = clean(user)
    want = {norm(t) for t in tags if len(norm(t)) >= 3}
    for n in names:
        m = ISSUE_FILE_RE.match(n)
        if not m:
            continue
        day, rest = m.groups()
        if rest == cu and day in days:
            kind = "same user and day, no session tag"
        elif rest.startswith(cu + "_") and norm(rest[len(cu) + 1:]) in want:
            kind = "session"
        else:
            continue
        titles = []
        try:
            with open(os.path.join(d, n), errors="replace") as fh:
                titles = [ln[3:].strip() for ln in fh if ln.startswith("## ")]
        except OSError:
            pass
        r["files"].append({"file": os.path.join(d, n), "match": kind, "n_entries": len(titles),
                           "titles": titles[:30]})
    return r


# ------------------------------------------------------------------- Data Quality Notes
def dq(severity, what, where=None, why=None, cause=None, fix=None, source=None):
    return {"severity": severity, "what": what, "where": where, "why": why, "cause": cause,
            "fix": fix, "source": source}


def data_quality_notes(rec):
    """The anomalies the skill already found, in the DataAnalysis five-part form -- what, where,
    why it matters, likely cause, suggested fix -- where the source says enough to fill them."""
    s, an = rec.get("search") or {}, rec.get("analysis") or {}
    notes = []
    if s.get("status") == "failed":
        notes.append(dq("CRITICAL", f"the search FAILED (exit {s.get('exit_code')})",
                        f"step {s.get('failing_step') or '?'}; {s.get('out_dir')}",
                        "no report, or an incomplete one -- nothing downstream is valid",
                        "; ".join((s.get("errors_excerpt") or [])[-2:]) or None,
                        "read the step's log (output/search/ in this record), fix, and resubmit "
                        "from the failed step (references/watcher.md)", "search status"))
    for w in s.get("warnings") or []:
        notes.append(dq("WARNING", w, s.get("out_dir"), source="record_run"))
    res = s.get("results") or {}
    if res.get("zero_id_runs"):
        notes.append(dq("WARNING", f"{len(res['zero_id_runs'])} run(s) identified no precursors",
                        ", ".join(res["zero_id_runs"][:10]),
                        "an empty run contributes only missing values and lowers completeness",
                        "a blank, wash or failed injection -- or a file DIA-NN could not read",
                        "check the run's MS1/MS2 signal in report.stats.tsv; leave non-samples "
                        "out of the inputs", res.get("stats_file")))
    det = s.get("detection") or {}
    params = s.get("parameters") or {}
    oru = det.get("orbitrap_resolution_unknown")
    if oru or (params.get("instrument_class") == "orbitrap_generic"):
        files = (oru or {}).get("files") if isinstance(oru, dict) else None
        notes.append(dq(
            "WARNING", "Orbitrap resolution not recorded for this search",
            f"{len(files)} file(s): {', '.join(os.path.basename(f) for f in files[:5])}"
            if files else "the search parameters (instrument class orbitrap_generic)",
            "with no resolution there is no documented DIA-NN mass tolerance: DIA-NN "
            "auto-calibrates per run, results depend on file order, and the 5-step parallel "
            "chain is declined",
            "ThermoRawFileParser does not report the Orbitrap resolution",
            "ask for the MS1/MS2 resolution (instrument method) and re-run estimate_params.py "
            "with --ms1-resolution/--ms2-resolution", det.get("file") or "params rationale"))
    for key, label in (("resolution_mixed", "files acquired at different Orbitrap resolutions"),
                       ("ms2_ion_trap", "MS2 read in the ion trap")):
        v = det.get(key)
        if v:
            files = v.get("files") if isinstance(v, dict) else v
            notes.append(dq("WARNING", label,
                            ", ".join(os.path.basename(str(x.get("file", x)) if isinstance(
                                x, dict) else str(x)) for x in (files or [])[:8]) or None,
                            "one mass tolerance is being applied to data that needs two",
                            None, "split the cohort, or confirm the tolerances per group",
                            det.get("file")))
    if det.get("low_confidence_files"):
        notes.append(dq("WARNING", "acquisition detected with low confidence",
                        ", ".join(os.path.basename(str(f)) for f in det["low_confidence_files"][:8]),
                        "a DDA file searched as DIA (or the reverse) gives wrong results",
                        None, "confirm DIA/DDA with the user", det.get("file")))
    if det.get("precursor_mz_range_mixed"):
        notes.append(dq("WARNING", "precursor m/z range differs between files", None,
                        "one --min/--max-pr-mz for all files cuts some files' windows", None,
                        "confirm the acquisition methods match", det.get("file")))
    for fw in det.get("file_warnings") or []:
        notes.append(dq("NOTE", fw["warning"], ", ".join(fw["files"][:6]),
                        source=det.get("file")))
    fa = s.get("fasta") or {}
    for w in fa.get("warnings") or []:
        notes.append(dq("WARNING", str(w), fa.get("meta_file"), source="FASTA sidecar"))
    if fa.get("n_contaminants_dropped_as_target"):
        notes.append(dq("NOTE", f"{fa['n_contaminants_dropped_as_target']} contaminant entries "
                                f"were identical to a target protein and were dropped",
                        fa.get("meta_file"),
                        "left in, they take the target's peptides and --cont-quant-exclude "
                        "removes them from quant", "the universal contaminant set holds bovine/"
                        "human/mouse proteins identical to real ones",
                        None, "FASTA sidecar"))
    if fa and fa.get("n_contaminants_appended") and not fa.get("contaminant_target_rule"):
        notes.append(dq(
            "WARNING", "legacy database: built before contaminants identical to a target protein "
                       "were removed", fa.get("meta_file"),
            "a contaminant entry identical to a real protein (bovine ACTB = human ACTB, ...) takes "
            "its peptides; --cont-quant-exclude Cont_ then drops them from quant and "
            "normalisation, so real proteins go missing without an error",
            "the FASTA sidecar has no contaminant_target_rule: fetch_fasta.py predates the check",
            "rebuild the FASTA with the current fetch_fasta.py and re-search, or run "
            "check_contaminant_competition.py on the pg_matrix", "FASTA sidecar"))
    if not fa and s:
        notes.append(dq("WARNING", "no FASTA sidecar, so the organism and database are not "
                                   "recorded", s.get("out_dir"), source="record_run"))
    for f in an.get("audit") or []:
        notes.append(dq("CRITICAL" if f.get("status") == "FAIL" else "WARNING",
                        f.get("message"), f"AUDIT.md: {f.get('check')}", source="AUDIT.md"))
    for flag in an.get("sample_quality_flags") or []:
        notes.append(dq("CRITICAL" if "CONFOUNDED" in flag.upper() else "WARNING",
                        re.sub(r"\*\*", "", flag), "SAMPLE_QUALITY.md",
                        source="sample_quality.py"))
    for w in an.get("warnings") or []:
        notes.append(dq("NOTE", w, an.get("session"), source="record_run"))
    for f in rec.get("findings") or []:
        notes.append(dq("WARNING", f.get("detail"), f.get("zip"), source=f.get("id")))
    for f in (rec.get("issues") or {}).get("files") or []:
        notes.append(dq("NOTE", f"{f['n_entries']} skill issue(s) recorded during this run: "
                                + "; ".join(f.get("titles")[:6]), f["file"],
                        "a skill problem may have changed what ran",
                        fix="read the issue file", source="report_issue.sh"))
    for n in rec.get("record_notes") or []:
        notes.append(dq("WARNING", n, rec.get("folder"), source="record_run"))
    p = rec.get("prot") or {}
    if not p.get("prot") and not p.get("id"):
        notes.append(dq("NOTE", "CoreOmics submission: not recorded", "this record",
                        "the submission is the only identifier that ties the analysis to what "
                        "the client asked for, who they are and what they sent",
                        "no --prot given, and no coreomics/prot key in the session metadata",
                        "re-record with --prot PROT_#### (the PI-surname search in CoreOmics, "
                        "lab=PROTEOMICS)", "record_run"))
    return notes


# ------------------------------------------------------------------------ SEARCH_LOG.md
def _row(k, v, src=None):
    return f"| {k} | {v} | {src or ''} |"


def _range(v, plain=False):
    """lo–hi. Counts get thousands separators; settings (m/z 357–1105) must read as typed."""
    if not v:
        return None
    lo, hi = (list(v) + [None, None])[:2]
    f = lambda x: "?" if x is None else ("%g" % x if plain and isinstance(x, (int, float))  # noqa
                                         else fmt_n(x))
    return f"{f(lo)}–{f(hi)}"


def render_dq(notes):
    L = ["## Data Quality Notes"]
    if not notes:
        return L + ["", "Nothing anomalous observed in what the skill recorded."]
    for i, n in enumerate(notes, 1):
        parts = [f"{i}. **{n['severity']}** -- {n['what']}"]
        for k, label in (("where", "Where"), ("why", "Why it matters"), ("cause", "Likely cause"),
                         ("fix", "Suggested fix")):
            if n.get(k):
                parts.append(f"   *{label}:* {n[k]}")
        if n.get("source"):
            parts.append(f"   *Source:* {n['source']}")
        L.append("\n".join(parts))
    return L


def render_parameters(s):
    p = s.get("parameters") or {}
    rec, why = p.get("record") or {}, p.get("why") or {}
    ma = p.get("mass_accuracy") or {}
    L = ["| Setting | Value | Where it came from |", "|---|---|---|"]

    def src(item, flag=None):
        parts = [item.get("source")] if isinstance(item, dict) else []
        if flag and why.get(flag):
            parts.append(why[flag])
        return "; ".join(x for x in parts if x) or None
    if rec.get("pr_mz"):
        L.append(_row("Precursor m/z range", _range(rec["pr_mz"]["value"], True),
                      src(rec["pr_mz"], "--min-pr-mz")))
    for key, label, flag in (("ms2_tol", "MS2 mass accuracy", "--mass-acc"),
                             ("ms1_tol", "MS1 mass accuracy", "--mass-acc-ms1")):
        if rec.get(key):
            t = rec[key]
            L.append(_row(label, f"{fmt_n(t['value'])} {t.get('unit') or ''}".strip(),
                          "; ".join(x for x in (t.get("source"), ma.get("source")) if x)))
        elif rec:
            L.append(_row(label, "not fixed -- the engine optimised it per run",
                          "; ".join(x for x in (rec.get("tol_note"),
                                                why.get(flag) or ma.get("source")) if x)))
    if ma.get("massacc_txt"):
        L.append(_row("Mass accuracy the chain pinned", f"`{ma['massacc_txt']}`",
                      "massacc.txt (measured by step 1b)"))
    if ma.get("measured"):
        L.append(_row("Mass accuracy measured", f"`{json.dumps(ma['measured'])[:160]}`",
                      "search_provenance.json result.mass_acc"))
    lm = s.get("engine_log_mass_accuracy") or {}
    o2, r1 = lm.get("optimised_ms2_ppm") or [], lm.get("recommended_ms1_ppm") or []
    if o2 or r1:
        L.append(_row("Engine log: mass accuracy", "; ".join(x for x in (
            f"optimised MS2 {', '.join(fmt_n(v) for v in o2)} ppm" if o2 else None,
            f"recommended MS1 {_range([min(r1), max(r1)], True)} ppm ({len(r1)} calibrations "
            f"logged -- DIA-NN logs one per run per pass)" if r1 else None) if x),
            os.path.basename(s.get("engine_log") or "") or None))
    res = p.get("resolution")
    if isinstance(res, dict) and (res.get("ms1") or res.get("ms2") or res.get("ms2_analyzer")):
        v = ", ".join(x for x in (f"MS1 {fmt_n(res.get('ms1'))}" if res.get("ms1") else None,
                                  f"MS2 {fmt_n(res.get('ms2'))}" if res.get("ms2") else None,
                                  f"MS2 in {res['ms2_analyzer']}" if res.get("ms2_analyzer")
                                  else None) if x)
        L.append(_row("Orbitrap resolution", v, res.get("source_label") or res.get("source")))
    elif p.get("instrument_class_label"):
        L.append(_row("Orbitrap resolution", "not recorded", p["instrument_class_label"]))
    sw = p.get("scan_window")
    if isinstance(sw, dict):
        L.append(_row("Scan window", fmt_n(sw.get("value")) if sw.get("value") is not None else
                      "not pinned", (sw.get("source") or "")[:220]))
    if rec.get("precursor_fdr"):
        L.append(_row("FDR", f"{fmt_n(rec['precursor_fdr']['value'])} "
                             f"({rec['precursor_fdr'].get('level') or 'precursor'})",
                      src(rec["precursor_fdr"], "--qvalue")))
    if rec.get("cleavage"):
        c, mc = rec["cleavage"], rec.get("missed_cleavages") or {}
        L.append(_row("Enzyme (in silico)",
                      f"{c.get('name') or '?'} (`{c.get('rule')}`)"
                      + (f", {mc['value']} missed cleavage(s)" if mc else ""),
                      src(c, "--cut")))
    if rec.get("pep_len"):
        L.append(_row("Peptide length", _range(rec["pep_len"]["value"], True),
                      src(rec["pep_len"])))
    if rec.get("pr_charge"):
        L.append(_row("Precursor charge", _range(rec["pr_charge"]["value"], True),
                      src(rec["pr_charge"])))
    if rec.get("mods"):
        L.append(_row("Modifications", "; ".join(
            f"{m.get('name')} ({m.get('targets') or m.get('position')}, {m.get('type')})"
            for m in rec["mods"]), None))
    if rec.get("library"):
        L.append(_row("Library", rec["library"]["value"], src(rec["library"])))
    if rec.get("mbr"):
        L.append(_row("Match-between-runs", "on", src(rec["mbr"])))
    if len(L) == 2:
        L.append(_row("(none parsed)", p.get("record_error") or "no parameters file was found",
                      ""))
    for w in rec.get("warnings") or []:
        L.append(f"\n_Note: {w}_")
    return L


def render_search(s):
    d, o = s.get("data") or {}, s.get("outputs") or {}
    eng = " ".join(str(x) for x in (s.get("engine_label") or s.get("engine") or "engine ?",
                                    s.get("engine_version") or "(version not recorded)"))
    exts = ", ".join(f"{n} × {e}" for e, n in sorted((d.get("extensions") or {}).items()))
    L = ["## Data",
         f"- **Raw files:** {fmt_n(d.get('n_files'))}" + (f" ({exts}; {d.get('vendor')})"
                                                         if exts else ""),
         f"- **Instrument:** {d.get('instrument') or 'not recorded'}"
         + (f" (source: {d['instrument_source']})" if d.get("instrument_source") else ""),
         f"- **Acquisition:** {d.get('acquisition') or 'not recorded'}",
         f"- **Where the raw files live:** `{d.get('raw_common_dir') or '?'}`"
         + (f" ({d['n_raw_dirs']} folders)" if (d.get("n_raw_dirs") or 0) > 1 else "")
         + " -- never copied; the list is `input/raw_files.txt`"]
    files = d.get("files") or []
    if files:
        L.append("- **Files:** " + ", ".join(f"`{os.path.basename(f.rstrip('/'))}`"
                                             for f in files[:12])
                 + (f", ... ({len(files)} in all)" if len(files) > 12 else ""))
    L += ["", "## Engine",
          f"- **{eng}** -- version source: {s.get('engine_version_source') or 'not recorded'}"]
    if s.get("engine_log_version_line"):
        L.append(f"- Engine log says: `{s['engine_log_version_line']}`")
    if s.get("engine_command"):
        L.append(f"- Binary: `{s['engine_command']}`")
    if s.get("search_mode"):
        L.append(f"- Route: {s['search_mode']}"
                 + (f" -- {s['routing_reason']}" if s.get("routing_reason") else ""))
    L += ["", "## Key parameters"] + render_parameters(s) + ["", "## Sequence database (FASTA)"]
    fa = s.get("fasta")
    if fa:
        nc = (fa.get("n_contaminants_appended") or 0) + (
            fa.get("n_contaminants_already_present") or 0)
        L += [f"- **Organism:** {fa.get('organism') or 'not recorded'}"
              + (f" (taxid {fa['taxid']})" if fa.get("taxid") else "")
              + (f" -- {fa['organism_source']}" if fa.get("organism_source") else ""),
              f"- **Proteome:** {fa.get('proteome') or '?'}"
              + (f" ({fa['proteome_type']})" if fa.get("proteome_type") else "")
              + f", {fmt_n(fa.get('n_proteome'))} entries, content "
                f"{fa.get('content_requested') or '?'} (used: {fa.get('content_used') or '?'})",
              f"- **Source:** `{fa.get('source') or '?'}`; UniProt release: "
              + (fa.get("uniprot_release") or (
                  f"not recorded (pre-staged copy dated {fa['staged_file_date']})"
                  if fa.get("staged_file_date") else "not recorded")),
              f"- **Contaminants:** {fmt_n(nc)} ({fa.get('contaminant_set') or 'none'})"
              + (f", {fa['n_contaminants_dropped_as_target']} dropped as identical to a target "
                 f"protein" if fa.get("n_contaminants_dropped_as_target") else "")
              + (f"; excluded from quantification with `--cont-quant-exclude "
                 f"{fa['cont_quant_exclude']}`" if fa.get("cont_quant_exclude") else "")
              + (f"; digestion enzymes kept as contaminants (`fetch_fasta.py --enzyme`): "
                 f"{', '.join(fa['digestion_enzymes_used'])}"
                 if fa.get("digestion_enzymes_used") else ""),
              f"- **Searched file:** `{fa.get('fasta') or '?'}` -- "
              f"{fmt_n(fa.get('n_sequences'))} sequences, md5 {fa.get('md5') or '?'}"]
    else:
        L.append("- not recorded (no <fasta>.meta.json was found)")
    r = s.get("results") or {}
    L += ["", "## Results"]
    if r.get("n_runs"):
        L += [f"- **Median precursors per run:** {fmt_n(r['median_precursors'])} "
              f"(range {_range(r['precursors_range'])}; {r['n_runs']} runs)",
              f"- **Median proteins per run:** {fmt_n(r['median_proteins'])} "
              f"(range {_range(r['proteins_range'])})",
              f"- Source: `{r.get('stats_file')}`"]
        per = sorted(r.get("per_run") or [], key=lambda x: x["precursors"])
        L += ["", "| Run | Precursors | Proteins |", "|---|---|---|"]
        L += [f"| {p['run']} | {fmt_n(p['precursors'])} | {fmt_n(p['proteins'])} |"
              for p in (per if len(per) <= 40 else per[:10])]
        if len(per) > 40:
            L.append(f"| ... {len(per) - 10} more (the 10 lowest are shown) | | |")
    else:
        L.append(f"- no per-run numbers: {r.get('error') or 'not recorded'}")
    L += ["", "## Where the outputs are", f"- **Search folder:** `{s.get('out_dir')}`",
          f"- **Report:** `{o.get('report') or 'none'}`"
          + (f" ({fmt_bytes(o.get('report_bytes'))})" if o.get("report_bytes") else "")
          + (" -- linked at `output/search/`" if o.get("report") else "")]
    if s.get("engine_log"):
        L.append(f"- **Engine log:** `{s['engine_log']}`")
    for q in o.get("quant") or []:
        L.append(f"- **Per-run .quant (left in place, never copied):** `{q['dir']}` -- "
                 f"{q['n_files']} files"
                 + (f", {fmt_bytes(q['bytes'])}" if q.get("bytes") is not None else ""))
    for x in o.get("xic") or []:
        L.append(f"- **XIC chromatograms:** `{x['dir']}` -- "
                 + (f"{x['n_files']} files" if "n_files" in x else
                    f"{x.get('n_task_dirs')} per-task folders"))
    for lib in o.get("libraries") or []:
        L.append(f"- **Library:** `{lib['file']}` ({fmt_bytes(lib['bytes'])})")
    if o.get("n_task_logs"):
        L.append(f"- **Per-file task logs:** {o['n_task_logs']} in the search folder")
    L += ["", "## FRAN"]
    fr = s.get("fran") or {}
    if fr.get("status"):
        L.append(f"- **{fr['status']}**"
                 + (f" -> `{fr['entry']}`" if fr.get("entry") else "")
                 + (f"; organism {fr['organism']}" if fr.get("organism") else "")
                 + (f" (taxon {fr['taxon']})" if fr.get("taxon") else "")
                 + (f"; staged by {fr['staged_by']}" if fr.get("staged_by") else "")
                 + (f"; {fr['n_linked']} item(s) linked" if fr.get("n_linked") else "")
                 + (f"; verify: {fr['verified'].get('state')}"
                    if isinstance(fr.get("verified"), dict) else ""))
    else:
        L.append(f"- {fr.get('detail') or 'not recorded'}")
    return L


def render_analysis(rec, an):
    de = an.get("de") or {}
    L = ["", f"## Analysis -- finalized {an.get('finalized') or '?'}",
         f"- **Session:** `{an.get('session')}`"]
    for d in an.get("docx") or []:
        L.append(f"- **{d['kind']}:** `output/{os.path.basename(d['file'])}`"
                 f" (original `{d['file']}`)")
    if de:
        L.append(f"- **DE:** {de.get('method') or '?'}; q {fmt_n(de.get('q_cutoff'))}, "
                 f"|logFC| {fmt_n(de.get('logfc'))}, adj.P {fmt_n(de.get('adjp'))}")
        sig = de.get("significant_per_contrast") or {}
        if sig:
            L.append("- **Significant proteins per contrast:** "
                     + ", ".join(f"{k} = {fmt_n(v)}" for k, v in sig.items()))
        elif de.get("contrasts"):
            L.append(f"- **Contrasts:** {', '.join(de['contrasts'])}")
    if an.get("audit_overall"):
        L.append(f"- **Audit (AUDIT.md):** {an['audit_overall']}")
    if an.get("methods_md"):
        L.append(f"- **Methods:** `{an['methods_md']}`")
    if an.get("data_submission"):
        L.append(f"- **DATA_SUBMISSION:** `{an['data_submission']}`"
                 + (" (start with HOW_TO_SUBMIT.md)" if an.get("how_to_submit") else ""))
    if an.get("reproduce_md"):
        L.append(f"- **Reproduce:** `{an['reproduce_md']}` (copy: `scripts/REPRODUCE.md`)")
    m = an.get("manifest") or {}
    if m:
        L.append(f"- **MANIFEST.txt:** {m.get('n_ok', 0)} part(s) OK, "
                 f"{len(m.get('skipped') or [])} skipped")
        L += [f"  - {x}" for x in (m.get("skipped") or [])[:10]]
    z, zc = an.get("zip") or {}, rec.get("zip_copy") or {}
    L.append(f"- **Session zip:** `{z.get('path')}`"
             + (f" ({fmt_bytes(z.get('bytes'))}, {fmt_n(z.get('n_entries'))} entries)"
                if z.get("exists") else "")
             + (f" -- copied to `{zc['rel']}`" + (f" ({z['note']})" if z.get("note") else "")
                if zc.get("copied") else
                f" -- not copied: {zc.get('reason') or z.get('reason') or '?'}"))
    return L


def render_log(rec):
    s, an = rec.get("search") or {}, rec.get("analysis")
    t, d = s.get("times") or {}, s.get("data") or {}
    eng = " ".join(str(x) for x in (s.get("engine_label") or s.get("engine") or "engine ?",
                                    s.get("engine_version") or "(version not recorded)"))
    ec = s.get("exit_code")
    L = [f"# Search log -- {rec.get('name')}", "",
         f"**{(s.get('status') or 'not recorded').upper()}**"
         + (f" (exit {ec})" if ec not in (None, "") else "")
         + (f" at step `{s['failing_step']}`" if s.get("failing_step") else "")
         + f" · {eng} · {fmt_n(d.get('n_files'))} file(s)"
         + (f" · {d['instrument']}" if d.get("instrument") else "")
         + (f" {d['acquisition']}" if d.get("acquisition") else "")
         + f" · run by {rec.get('user')}", "",
         f"**CoreOmics submission:** {prot_label(rec.get('prot'))}"
         + (f" (source: {rec['prot']['source']})" if (rec.get("prot") or {}).get("source")
            else ""), ""]
    if not s:
        L += ["_The search itself is not recorded yet -- no search out dir was readable when "
              "this was written. The analysis is below._", ""]
    rb = s.get("reported_by") or {}
    L += ["| | |", "|---|---|",
          f"| Run by | {rec.get('user')} (HIVE account) |",
          f"| Submitted | {t.get('submitted') or '?'} |",
          f"| Started | {t.get('started') or '?'} |",
          f"| Finished | {t.get('finished') or '?'} |",
          f"| Status | {s.get('status') or '?'} -- {s.get('status_source') or ''} |",
          f"| Exit code | {ec if ec not in (None, '') else '?'}"
          + (f" ({s['exit_code_source']})" if s.get("exit_code_source") else "") + " |"]
    if rb:
        L.append(f"| Reported by | job {rb.get('job_id') or '?'} `{rb.get('job_name') or '?'}`"
                 + (f" task {rb['array_task']}" if rb.get("array_task") else "")
                 + (f" on {rb['node']}" if rb.get("node") else "") + " |")
    L += [f"| Times from | {t.get('source') or '?'} |", ""]
    L += render_dq(rec.get("data_quality_notes") or []) + [""]
    if s.get("errors_excerpt"):
        L += ["## Errors (from the engine and job logs)", "```"] + s["errors_excerpt"] + ["```",
                                                                                           ""]
    if s.get("jobs"):
        L += ["## SLURM jobs", "| Job | Name | Tasks | State | Exit | Submitted -> ended | "
              "Node(s) | Partition / account |", "|---|---|---|---|---|---|---|---|"]
        for g in s["jobs"]:
            L.append(f"| {g['id']} | {g.get('name') or ''} | {g['tasks']} | "
                     + ", ".join(f"{k} {v}" if g["tasks"] > 1 else k
                                 for k, v in sorted(g["states"].items()))
                     + f" | {', '.join(g.get('exit_codes') or [])} | {g.get('submit') or '?'} -> "
                       f"{g.get('end') or '?'} | {g.get('nodes') or ''} | "
                       f"{g.get('partition') or ''} / {g.get('account') or ''} |")
        L.append("")
    elif s.get("jobs_note"):
        L += [f"_SLURM: {s['jobs_note']}._", ""]
    if s:
        L += render_search(s)
    iss = rec.get("issues") or {}
    L += ["", "## Skill issues recorded for this session"]
    if iss.get("files"):
        L += [f"- `{f['file']}` ({f['match']}): {f['n_entries']} entr"
              f"{'y' if f['n_entries'] == 1 else 'ies'}"
              + (" -- " + "; ".join(f["titles"][:8]) if f.get("titles") else "")
              for f in iss["files"]]
    else:
        L.append(f"- none found in `{iss.get('dir') or issues_dir()}` for "
                 f"{', '.join(rec.get('issue_tags') or []) or 'this run'}"
                 + (f" ({iss['note']})" if iss.get("note") else ""))
    if an:
        L += render_analysis(rec, an)
        if an.get("expert_review_notes"):
            L += ["", "## Expert Review Notes", "", an["expert_review_notes"]]
    L += ["", "## Copies in this folder"]
    L += [f"- `{rel}` <- `{c.get('from')}` ({fmt_bytes(c.get('bytes'))})"
          for rel, c in sorted((rec.get("copies") or {}).items())] or ["- none"]
    nc = [x for sec in (rec.get("not_copied") or {}).values() for x in sec]
    if nc:
        L += ["", "## Not copied"] + [f"- `{x.get('src')}` -- {x.get('reason')}" for x in nc]
    if rec.get("links"):
        L += ["", "## Links (symlinks, not copies)"] + [
            f"- `{k}` -> `{v}`" for k, v in sorted(rec["links"].items())]
    L += ["", "---", "Record history:"] + [
        f"- {h.get('at')} {h.get('event')} by {h.get('by')} on {h.get('host')} "
        f"(route {h.get('route')}, skill {h.get('skill_version') or '?'})"
        for h in rec.get("history") or []]
    L += ["", f"_Written by record_run.py (schema {SCHEMA_VERSION}). The same facts, "
              f"machine-readable: `run_record.json`. Master log: `../../data_analysis.md`._"]
    return "\n".join(L) + "\n"


# --------------------------------------------------------------------- finding the folder
def index_links(identity):
    return [(k, f"{k}_{key16(v)}") for k, v in (("out", identity.get("out")),
                                               ("session", identity.get("session"))) if v]


def same_identity(rec, identity):
    ri = (rec or {}).get("identity") or {}
    if identity.get("out") and ri.get("out"):
        return identity["out"] == ri["out"]
    return bool(identity.get("session") and identity["session"] == ri.get("session"))


def find_record(root, identity):
    """The session folder already holding this search: through the index first (O(1)), then
    the folder whose run_record.json names the same identity."""
    for _, name in index_links(identity):
        p = os.path.join(root, ".index", name)
        if os.path.islink(p):
            tgt = os.path.normpath(os.path.join(os.path.dirname(p), os.readlink(p)))
            if os.path.isdir(tgt):
                rec, state = read_json_retry(os.path.join(tgt, "run_record.json"))
                if state != "ok" or same_identity(rec, identity):
                    return tgt                # the index is the claim; a bad record is not a miss
    return None


def claim_folder(root, base, identity, dry_run):
    """(folder, is_new): `sessions/<base>`, or `<base>_2`, `_3` ... when that name is a
    DIFFERENT search. mkdir is the claim, so two writers never take one name; a folder another
    writer has only just claimed is given a moment to show whose it is."""
    sd = os.path.join(root, "sessions")
    for i in range(1, 200):
        name = base if i == 1 else f"{base}_{i}"
        d = os.path.join(sd, name)
        if not os.path.exists(d):
            if dry_run:
                return d, True
            os.makedirs(sd, mode=0o2770, exist_ok=True)
            try:
                os.mkdir(d, 0o2770)
                return d, True
            except FileExistsError:
                pass
        for _ in range(1 if dry_run else 8):
            rec, state = read_json_retry(os.path.join(d, "run_record.json"))
            if state == "ok":
                break
            if find_record(root, identity) == d:
                return d, False
            time.sleep(0.25)
        if state == "ok" and same_identity(rec, identity):
            return d, False
    raise RuntimeError(f"no free folder name for {base}")


def link_index(root, identity, folder):
    idx = os.path.join(root, ".index")
    os.makedirs(idx, mode=0o2770, exist_ok=True)
    rel = os.path.relpath(folder, idx)
    for _, name in index_links(identity):
        p = os.path.join(idx, name)
        try:
            if os.path.islink(p) and os.readlink(p) == rel:
                continue
            tmp = f"{p}.tmp{os.getpid()}"
            os.symlink(rel, tmp)
            os.replace(tmp, p)
        except OSError as e:
            say(f"could not write index link {p}: {e}")


# ------------------------------------------------------------------------------ writing
def merge(old, plan, folder, route_name):
    rec = dict(old or {})
    rec["schema_version"] = SCHEMA_VERSION
    rec.setdefault("created", now_iso())
    rec["updated"] = now_iso()
    rec["folder"] = folder
    rec["name"] = os.path.basename(folder)
    ident = dict(rec.get("identity") or {})
    for k, v in (plan.get("identity") or {}).items():
        ident[k] = ident.get(k) or v
    rec["identity"] = ident
    for k in ("label", "user", "date"):
        rec[k] = rec.get(k) or plan.get(k)
    if plan.get("prot"):
        rec["prot"] = plan["prot"]
    rec.setdefault("prot", None)
    if plan.get("search"):
        s, prev = dict(plan["search"]), rec.get("search") or {}
        # An inference that cannot tell (no report, no job state) must not erase what the job
        # itself said when it finished.
        keep = ("status", "status_source", "exit_code", "exit_code_source", "failing_step",
                "reported_by")
        if s.get("status") in (None, "unknown", "running") and prev.get("status") in (
                "completed", "failed"):
            s["status_inferred_now"] = s.get("status")
            for k in keep:
                s[k] = prev.get(k)
            s["times"] = dict(s.get("times") or {}, finished=(
                (s.get("times") or {}).get("finished") or (prev.get("times") or {}).get("finished")))
        elif (s.get("status") == prev.get("status") and not str(s.get("status_source") or "")
              .startswith("--status") and str(prev.get("status_source") or "")
              .startswith("--status")):
            # The same outcome, re-inferred at finalize: the job's own word (its exit code, the
            # step that reported) is the better record of it.
            for k in keep:
                s[k] = prev.get(k)
        rec["search"] = s
    if plan.get("analysis"):
        rec["analysis"] = plan["analysis"]
    for k in ("out", "session"):
        if plan.get(k):
            rec["search_out" if k == "out" else "session"] = plan[k]
    rec["issue_tags"] = list(dict.fromkeys((rec.get("issue_tags") or [])
                                           + (plan.get("issue_tags") or [])))
    rec["copies"] = dict(rec.get("copies") or {})
    nc = dict(rec.get("not_copied") or {})
    for sec in {c["section"] for c in plan["copies"] + plan["not_copied"]}:
        nc[sec] = [{"src": x["src"], "reason": x["reason"]} for x in plan["not_copied"]
                   if x["section"] == sec]
    rec["not_copied"] = nc
    f = {x["id"]: x for x in rec.get("findings") or []}
    f.update({x["id"]: x for x in plan.get("findings") or []})
    rec["findings"] = list(f.values())
    rec["skill_version"] = plan.get("skill_version") or rec.get("skill_version")
    rec.setdefault("master_log", [])
    rec.setdefault("activity_logged", [])
    rec["history"] = ((rec.get("history") or []) + [{
        "at": now_iso(), "event": plan["event"], "by": plan.get("by"), "host": plan.get("host"),
        "route": route_name, "skill_version": plan.get("skill_version")}])[-100:]
    return rec


def master_entries(rec, event):
    """[(marker, text)] for data_analysis.md: the run's first entry is a `## date:` block; every
    later event is a dated `###` line under it. The marker keeps an event from being logged
    twice for the same status."""
    s, an = rec.get("search") or {}, rec.get("analysis") or {}
    name, key = rec["name"], key16(rec["identity"].get("out") or rec["identity"].get("session"))
    folder = f"sessions/{name}/"
    today = datetime.date.today().isoformat()
    t = s.get("times") or {}
    s_day = (t.get("finished") or t.get("submitted") or today)[:10]      # when it happened
    a_day = (an.get("finalized") or today)[:10]
    first = not rec.get("master_log")
    search_logged = any(" search " in m for m in rec.get("master_log") or [])
    d = s.get("data") or {}
    fa = s.get("fasta") or {}
    res = s.get("results") or {}
    out = []
    if s and (event == "search-done" or not search_logged):
        status = s.get("status") or "unknown"
        marker = f"<!-- record_run {key} search {status} -->"
        what = (f"{s.get('engine_label') or s.get('engine') or 'search'} "
                f"{s.get('engine_version') or ''} search, {fmt_n(d.get('n_files'))} × "
                f"{'/'.join(sorted(d.get('extensions') or {})) or 'files'}"
                + (f" ({d.get('instrument')}, {d.get('acquisition') or '?'})"
                   if d.get("instrument") else "")
                + (f", {fa.get('organism')}" if fa.get("organism") else "")
                + f" -- {status.upper()}"
                + (f" (exit {s.get('exit_code')}, step {s.get('failing_step')})"
                   if status == "failed" else ""))
        body = [f"- **Session:** [{folder}]({folder}) -- `SEARCH_LOG.md`",
                f"- **CoreOmics submission:** {prot_label(rec.get('prot'))}",
                f"- **Run by:** {rec.get('user')} · **Search folder:** `{s.get('out_dir')}`"]
        if res.get("n_runs"):
            body.append(f"- **Result:** median {fmt_n(res['median_precursors'])} precursors / "
                        f"{fmt_n(res['median_proteins'])} proteins per run")
        if first:
            text = f"\n## {s_day}: {what}\n\n" + "\n".join(body)
        else:                                  # a later event: one dated line, not a new entry
            text = (f"\n### {s_day} update -- {name}: search {status}"
                    + (f" (exit {s.get('exit_code')}, step {s.get('failing_step')})"
                       if status == "failed" else "") + "\n\n"
                    + "\n".join(b for b in body if b.startswith("- **Result")))
        out.append((marker, text.rstrip("\n") + f"\n{marker}\n"))
        first = False
    if event == "analysis-done" and an:
        marker = f"<!-- record_run {key} analysis {an.get('finalized') or 'complete'} -->"
        body = []
        for x in an.get("docx") or []:              # the report of record comes first
            body.append(f"- **{x['kind']}:** `{folder}output/{os.path.basename(x['file'])}`")
        sig = (an.get("de") or {}).get("significant_per_contrast") or {}
        if sig:
            body.append("- **Significant proteins:** " + ", ".join(
                f"{k} = {fmt_n(v)}" for k, v in sig.items()))
        zc, z = rec.get("zip_copy") or {}, an.get("zip") or {}
        if z.get("exists"):
            body.append(f"- **Reproducibility / zip:** "
                        + (f"`{folder}{zc['rel']}`" if zc.get("copied") else
                           f"`{z.get('path')}` (not copied: {zc.get('reason') or z.get('reason')})"))
        if an.get("reproduce_md"):
            body.append(f"- **Reproduce:** `{folder}scripts/REPRODUCE.md`")
        body.append(f"- **CoreOmics submission:** {prot_label(rec.get('prot'))}")
        head = (f"\n## {a_day}: analysis of {name} -- finalized\n\n- **Session:** "
                f"[{folder}]({folder}) -- `SEARCH_LOG.md`\n" if first else
                f"\n### {a_day} update -- {name}: analysis complete\n\n")
        out.append((marker, head + "\n".join(body) + f"\n{marker}\n"))
    return out


def activity_rows(rec, event, plan):
    """[(dedupe key or None, csv line)]. State events (submitted, started, FRAN, an issue) are
    logged once per run; the call's own event (completed / failed / analysis) on every call."""
    s, an = rec.get("search") or {}, rec.get("analysis") or {}
    name, t = rec["name"], (s.get("times") or {})
    prot = prot_label(rec.get("prot"))
    prot_note = f"; CoreOmics {prot}" if prot != "not recorded" else ""
    tool = " ".join(str(x) for x in (s.get("engine_label") or s.get("engine"),
                                     s.get("engine_version")) if x) or "search"
    out = []
    if s:
        jobs = ",".join(g["id"] for g in s.get("jobs") or [])
        if t.get("submitted") and (s.get("times") or {}).get("source", "").startswith("sacct"):
            out.append(("search_submitted", csv_row([ts_min(t["submitted"]), name,
                                                     "search_submitted", "SLURM", s["out_dir"],
                                                     "ok", f"jobs {jobs}{prot_note}"])))
        if t.get("started"):
            out.append(("search_started", csv_row([ts_min(t["started"]), name, "search_started",
                                                   "SLURM", s["out_dir"], "ok", f"jobs {jobs}"])))
        if event == "search-done" and s.get("status") in ("completed", "failed"):
            res = s.get("results") or {}
            again = f"search_{s['status']}" in rec.get("activity_logged", [])
            notes = "; ".join(x for x in (
                f"exit {s.get('exit_code')}" if s.get("exit_code") not in (None, "") else None,
                f"step {s.get('failing_step')}" if s.get("failing_step") else None,
                f"median {fmt_n(res.get('median_precursors'))} precursors / "
                f"{fmt_n(res.get('median_proteins'))} proteins per run" if res.get("n_runs")
                else None,
                "re-recorded" if again else None) if x) + prot_note
            out.append((None, csv_row([ts_min(t.get("finished")), name,
                                       f"search_{s['status']}", tool, s["out_dir"], s["status"],
                                       notes])))
            rec["activity_logged"].append(f"search_{s['status']}")
        fr = s.get("fran") or {}
        if fr.get("status") in ("staged", "ingested"):
            out.append((f"fran_{fr['status']}", csv_row([ts_min(), name, "fran_staged",
                                                          "fran_deposit.py", fr.get("entry"),
                                                          fr["status"],
                                                          fr.get("organism") or ""])))
        elif event == "analysis-done":
            out.append(("fran_skipped", csv_row([ts_min(), name, "fran_skipped",
                                                 "fran_deposit.py", s["out_dir"], "skipped",
                                                 fr.get("status") or fr.get("detail")
                                                 or "no receipt"])))
    for f in (rec.get("issues") or {}).get("files") or []:
        for title in f.get("titles") or []:
            out.append((f"issue:{os.path.basename(f['file'])}:{title}",
                        csv_row([ts_min(), name, "issue_recorded", "report_issue.sh", f["file"],
                                 "recorded", title])))
    if event == "analysis-done" and an:
        sig = (an.get("de") or {}).get("significant_per_contrast") or {}
        zc = rec.get("zip_copy") or {}
        notes = "; ".join(x for x in (
            ", ".join(f"{k}={v}" for k, v in sig.items()) or None,
            ("zip copied" if zc.get("copied") else f"zip not copied: {zc.get('reason')}")
            if (an.get("zip") or {}).get("exists") else "no session zip") if x) + prot_note
        out.append((None, csv_row([ts_min(), name, "analysis_completed",
                                   f"session.py finalize ({(an.get('de') or {}).get('method') or '?'})",
                                   an.get("session"), "completed", notes])))
    return out


def read_record(folder, notes):
    """This folder's run_record.json, or None. One that stays unparseable is MOVED ASIDE, never
    deleted or overwritten, and the move is noted in the new record."""
    path = os.path.join(folder, "run_record.json")
    # Longer than other reads: across HIVE nodes even a tmp + fsync + rename file can read torn
    # for a moment (measured 2026-09-24), and a record wrongly judged corrupt loses its history.
    rec, state = read_json_retry(path, tries=7)
    if state != "unreadable":
        return rec
    aside = f"{path}.unreadable-{int(time.time())}"
    try:
        os.replace(path, aside)
        notes.append(f"run_record.json could not be read and was kept as {os.path.basename(aside)};"
                     f" this record was rebuilt from what is readable now")
    except OSError as e:
        notes.append(f"run_record.json could not be read ({e}); it was left in place")
        raise RuntimeError(f"{path} is unreadable and could not be moved aside: {e}")
    return None


def execute(plan, a, deadline, route_name):
    """Write (or, with --dry-run, describe) the record. Runs where the registry is writable: on
    HIVE, or at the far end of the SSH route. A copy's `here` path, when set, is the uploaded
    copy of a file that exists only on the machine that sent it.

    Three phases, so the record's lock is never held through a long copy: (A) under the lock,
    read the record, merge this event in, write it; (B) copy, link, unlocked -- every file is
    written to a temporary name and renamed; (C) under the lock again, re-read the record --
    another writer may have updated it meanwhile -- apply this call's copy results to THAT, append
    to the shared logs, write."""
    root = runs_dir()
    ident = plan.get("identity") or {}
    if not (ident.get("out") or ident.get("session")):
        return not_recorded("nothing_to_record",
                            "no readable search out dir or session to key the record on")
    old_umask = os.umask(0o007)          # group read/write; the registry is setgid proteomics-grp
    try:
        return _execute(plan, a, deadline, route_name, root, ident)
    finally:
        os.umask(old_umask)


def _execute(plan, a, deadline, route_name, root, ident):
    folder = find_record(root, ident)
    new = False
    if folder is None:
        base = clean(plan.get("name") or f"{plan.get('date') or now_iso()[:10]}_"
                                          f"{plan.get('label') or 'search'}", 90)
        folder, new = claim_folder(root, base, ident, a.dry_run)
    z = plan.get("zip") or {}
    notes = []
    res = {"recorded": False, "path": folder, "new": new, "route": route_name,
           "event": plan["event"], "session": os.path.basename(folder)}

    def build(old):
        rec = merge(old, plan, folder, route_name)
        rec["record_notes"] = (rec.get("record_notes") or []) + notes
        t = (rec.get("search") or {}).get("times") or {}
        days = {x[:10] for x in (t.get("submitted"), t.get("finished")) if x}
        rec["issues"] = find_issues(rec.get("user") or "", rec.get("issue_tags") or [], days)
        rec["links"] = dict(rec.get("links") or {}, **links)
        if z.get("exists"):
            rec["zip_copy"] = ({"copied": True, "rel": z["dest_rel"], "planned": True}
                               if z.get("copy") else {"copied": False, "reason": z.get("reason")})
        for c in plan["copies"]:
            rec["copies"][c["rel"]] = {"from": c["src"], "bytes": c["bytes"], "planned": True}
        rec["data_quality_notes"] = data_quality_notes(rec)
        return rec

    links = {k: v for k, v in (plan.get("links") or {}).items() if v and os.path.exists(v)}
    if a.dry_run:
        rec = build(None if new else read_json_retry(os.path.join(folder, "run_record.json"))[0])
        say(f"DRY RUN ({route_name}) -- would {'create' if new else 'update'} {folder}")
        say("---------- SEARCH_LOG.md ----------\n" + render_log(rec))
        for c in plan["copies"]:
            say(f"would copy {c['rel']} ({fmt_bytes(c['bytes'])}) <- {c['src']}")
        for c in plan["not_copied"]:
            say(f"would NOT copy {c['src']} -- {c['reason']}")
        for k, v in links.items():
            say(f"would link {k} -> {v}")
        for marker, text in master_entries(rec, plan["event"]):
            say("data_analysis.md would get:" + text)
        return dict(res, reason="dry_run", detail=f"nothing written; would "
                                                  f"{'create' if new else 'update'} {folder}")

    os.makedirs(folder, mode=0o2770, exist_ok=True)
    link_index(root, ident, folder)
    try:
        ensure_readme(root, deadline)
    except OSError as e:
        say(f"README.md not refreshed: {e}")
    rec_path = os.path.join(folder, "run_record.json")
    # ---- (A) merge this event into the record, under its lock
    with DirLock(rec_path, deadline):
        rec = build(None if new else read_record(folder, notes))
        rec["pending"] = "copying files" if plan["copies"] or z.get("copy") else None
        write_text(rec_path, json.dumps(rec, indent=2, default=str))
        write_text(os.path.join(folder, "SEARCH_LOG.md"), render_log(rec))

    # ---- (B) copies, generated files and links: each written aside and renamed into place
    done, failed, zip_result, link_errors = {}, [], None, []
    for rel, text in (plan.get("generated") or {}).items():
        if rel in (rec.get("copies") or {}) and not rec["copies"][rel].get("generated"):
            continue                          # the session's own copy of this file wins
        os.makedirs(os.path.dirname(os.path.join(folder, rel)), exist_ok=True)
        write_text(os.path.join(folder, rel), text)
        done[rel] = {"from": "(written from search_provenance.json)",
                     "bytes": len(text.encode()), "generated": True}
    for c in plan["copies"]:
        try:
            copy_file(c.get("here") or c["src"], os.path.join(folder, c["rel"]), deadline)
            done[c["rel"]] = {"from": c["src"], "bytes": c["bytes"], "copied_at": now_iso()}
        except (OSError, TimeoutError, Stop) as e:
            failed.append((c, f"copy failed: {e}"))
    if z.get("copy"):
        dst = os.path.join(folder, z["dest_rel"])
        try:
            if z.get("here"):
                copy_file(z["here"], dst, deadline)
            elif z.get("drop"):
                copy_zip_without(z["path"], dst, z["drop"], deadline)
            else:
                copy_file(z["path"], dst, deadline)
            zip_result = {"copied": True, "rel": z["dest_rel"], "bytes": fsize(dst),
                          "copied_at": now_iso()}
        except (OSError, TimeoutError, ValueError, zipfile.BadZipFile, Stop) as e:
            zip_result = {"copied": False, "reason": f"copy failed: {e}"}
    for rel, target in links.items():
        p = os.path.join(folder, rel)
        try:
            os.makedirs(os.path.dirname(p), exist_ok=True)
            if os.path.islink(p) and os.readlink(p) == target:
                continue
            tmp = f"{p}.tmp.{os.getpid()}"
            os.symlink(target, tmp)
            os.replace(tmp, p)
        except OSError as e:
            link_errors.append(f"{rel}: {e}")

    # ---- (C) apply the results to the CURRENT record, append to the shared logs, write
    with DirLock(rec_path, deadline):
        cur, state = read_json_retry(rec_path)
        if state != "ok" or not same_identity(cur, ident):
            cur = rec
        for c in plan["copies"]:
            cur["copies"].pop(c["rel"], None)
        cur["copies"].update(done)
        for c, why in failed:
            cur["not_copied"].setdefault(c["section"], []).append({"src": c["src"],
                                                                    "reason": why})
        if zip_result:
            cur["zip_copy"] = zip_result
        if link_errors:
            cur["link_errors"] = link_errors
        cur.pop("pending", None)
        cur["data_quality_notes"] = data_quality_notes(cur)
        cur.setdefault("master_log", [])
        cur.setdefault("activity_logged", [])
        for marker, text in master_entries(cur, plan["event"]):
            if marker in cur["master_log"]:
                continue
            try:
                append_locked(os.path.join(root, "data_analysis.md"), [text], deadline,
                              header=MASTER_HEADER, marker=marker)
                cur["master_log"].append(marker)
            except OSError as e:
                say(f"data_analysis.md not updated: {e}")
        rows = [(k, r) for k, r in activity_rows(cur, plan["event"], plan)
                if not (k and k in cur["activity_logged"])]
        if rows:
            try:
                append_locked(os.path.join(root, "activity_log.csv"), [r for _, r in rows],
                              deadline, header=",".join(ACTIVITY_COLUMNS) + "\n")
                cur["activity_logged"] += [k for k, _ in rows if k]
            except OSError as e:
                say(f"activity_log.csv not updated: {e}")
        cur["activity_logged"] = list(dict.fromkeys(cur["activity_logged"]))[-500:]
        write_text(rec_path, json.dumps(cur, indent=2, default=str))
        write_text(os.path.join(folder, "SEARCH_LOG.md"), render_log(cur))
    zc = cur.get("zip_copy") or {}
    out = dict(res, recorded=True, search_log=os.path.join(folder, "SEARCH_LOG.md"),
               status=(cur.get("search") or {}).get("status"),
               prot=prot_label(cur.get("prot")), n_copied=len(cur["copies"]),
               zip_copied=zc.get("copied") if z.get("exists") else None,
               n_data_quality_notes=len(cur["data_quality_notes"]),
               findings=[f["id"] for f in cur.get("findings") or []])
    return out


# ------------------------------------------------------------------------ the SSH route
def ship(plan, a, deadline):
    """Upload what exists only here, plus this script and the modules it reads with, and have
    HIVE write the record. hive_exec.sh shares one multiplexed connection across the calls."""
    exe, user = hive_exec_path(), hive_login()
    stage = tempfile.mkdtemp(prefix="record_run_")
    try:
        sd = os.path.join(stage, "scripts")
        os.makedirs(sd)
        for f in glob.glob(os.path.join(HERE, "*.py")) + glob.glob(os.path.join(HERE, "*.sh")):
            shutil.copy2(f, sd)
        pj = os.path.join(HERE, "..", ".claude-plugin", "plugin.json")
        if os.path.isfile(pj):
            os.makedirs(os.path.join(stage, ".claude-plugin"))
            shutil.copy2(pj, os.path.join(stage, ".claude-plugin"))
        if not a.dry_run:
            for c in list(plan["copies"]):
                try:
                    copy_file(c["src"], os.path.join(stage, "files", c["rel"]), deadline)
                except (OSError, TimeoutError) as e:
                    plan["copies"].remove(c)
                    plan["not_copied"].append(dict(c, reason=f"upload failed: {e}"))
            z = plan.get("zip") or {}
            ssh_cap = int(env_num("RECORD_RUN_SSH_ZIP_CAP_GB", 1) * GB)
            if z.get("copy") and z.get("kept_bytes", 0) > ssh_cap:
                z.update(copy=False, reason=(
                    f"{fmt_bytes(z['kept_bytes'])} is over the {fmt_bytes(ssh_cap)} cap for "
                    f"uploading from {socket.gethostname()} (RECORD_RUN_SSH_ZIP_CAP_GB); the "
                    f"zip stays at {z['path']} on that machine"))
            elif z.get("copy"):
                dst = os.path.join(stage, "files", z["dest_rel"])
                if z.get("drop"):
                    copy_zip_without(z["path"], dst, z["drop"], deadline)
                else:
                    copy_file(z["path"], dst, deadline)
        with open(os.path.join(stage, "bundle.json"), "w") as fh:
            json.dump({"plan": plan, "uploaded_from": socket.gethostname(),
                       "args": {k: getattr(a, k, None) for k in (
                           "status", "exit_code", "name", "file_cap_mb", "zip_cap_gb", "step")}},
                      fh, indent=2, default=str)
        rd = f'"$HOME/{UPLOAD_DIR}/{os.path.basename(stage)}"'
        for cmd in (["bash", exe, f'mkdir -p "$HOME/{UPLOAD_DIR}"'],
                    ["bash", exe, "--put", stage, UPLOAD_DIR]):
            p = subprocess.run(cmd, capture_output=True, text=True, stdin=subprocess.DEVNULL,
                               timeout=max(5, deadline.left() - 2))
            if p.returncode != 0:
                return not_recorded("ssh_failed", f"upload to HIVE failed (exit {p.returncode}): "
                                                  f"{(p.stderr or '').strip()[-200:]}",
                                    route="ssh")
        env = " ".join(f"{k}={shlex.quote(os.environ[k])}" for k in (
            "SKILL_RUNS_DIR", "SKILL_ISSUES_DIR") if os.environ.get(k))
        cmd = (f"{env + ' ' if env else ''}python3 {rd}/scripts/record_run.py merge --bundle {rd} "
               f"--remote-hop --timeout {int(max(5, deadline.left() - 4))}"
               f"{' --dry-run' if a.dry_run else ''}; rc=$?; rm -rf {rd}; exit $rc")
        p = subprocess.run(["bash", exe, cmd], capture_output=True, text=True,
                           stdin=subprocess.DEVNULL, timeout=max(5, deadline.left()))
        if p.stderr.strip():
            sys.stderr.write(p.stderr)
        res = json_tail(p.stdout, "{")
        if not res:
            return not_recorded("ssh_failed", f"HIVE did not answer (exit {p.returncode}): "
                                              f"{((p.stderr or '') + (p.stdout or '')).strip()[-200:]}",
                                route="ssh")
        res["route"] = "ssh"
        if res.get("path"):
            res["path"] = f"{user}@hive:{res['path']}"
        return res
    finally:
        shutil.rmtree(stage, ignore_errors=True)


# ------------------------------------------------------------------------------ commands
def gather(event, a, deadline, r):
    """The plan, from whatever is readable HERE. What is not (a search out dir that lives on
    HIVE, seen from a laptop) stays a path for the HIVE side to read."""
    plan = new_plan(event, a)
    sess = os.path.abspath(a.session) if a.session and os.path.isdir(a.session) else a.session
    plan["session"] = sess
    if sess and os.path.isdir(sess):
        plan["identity"]["session"] = os.path.realpath(sess)
    out = a.out
    if event == "analysis-done" and sess and os.path.isdir(sess):
        out = session_search_out(sess, a.out)
    plan["out"] = os.path.abspath(out) if out and os.path.isdir(out) else out
    if sess:
        plan["issue_tags"] = list(dict.fromkeys(
            plan["issue_tags"] + [os.path.basename(str(sess).rstrip("/")),
                                  label_for("", sess, None)]))
    if out and os.path.isdir(out):
        plan_search(plan, out, a, deadline)
    if event == "analysis-done" and sess and os.path.isdir(sess):
        plan_analysis(plan, sess, a, plan["zip_cap"])
    plan["prot"] = find_prot(a.prot, sess if sess and os.path.isdir(str(sess)) else None,
                             out if out and os.path.isdir(str(out)) else None)
    # Over SSH the record is attributed to the HIVE account, never this machine's login name.
    if r == "ssh":
        plan["user"] = hive_login() or plan["user"]
    if not plan["user"]:
        plan["user"] = (owner(sess) if sess and os.path.isdir(str(sess)) else None) or whoami()
    plan["date"] = plan["date"] or now_iso()[:10]
    plan["label"] = label_for(plan["out"] or "", sess, a.name)
    plan["name"] = folder_name(plan, sess, a.name)
    return plan


def folder_name(plan, session, explicit):
    """YYYY-MM-DD_Short-Description: the skill's session name, which already has that form;
    else <date>_<out-dir name>."""
    if explicit:
        base = clean(explicit, 90)
    elif session:
        base = clean(os.path.basename(str(session).rstrip("/")), 90)
    else:
        base = clean(label_for(plan.get("out") or "", None, None), 70)
    if not re.match(r"^\d{4}-\d{2}-\d{2}_", base):
        base = f"{plan.get('date') or now_iso()[:10]}_{base}"
    return base


def do_event(event, a, deadline):
    r = route(remote_hop=a.remote_hop)
    if r in REASONS:
        return not_recorded(REASONS[r], where_text(r).split(": ", 1)[1], route=r)
    if event == "search-done" and not a.out:
        return not_recorded("bad_input", "--out is required")
    if event == "analysis-done" and not a.session:
        return not_recorded("bad_input", "--session is required")
    for opt, v in (("--out", a.out), ("--session", a.session)):
        if v and os.path.exists(v) and not os.path.isdir(v):
            return not_recorded("bad_input", f"{opt} {v} is not a directory")
    plan = gather(event, a, deadline, r)
    if r == "direct":
        if not plan.get("search") and not plan.get("analysis"):
            what = a.out if event == "search-done" else a.session
            return not_recorded("out_not_found" if event == "search-done" else
                                "session_not_found", f"{what} does not exist here", route=r)
        return execute(plan, a, deadline, r)
    return ship(plan, a, deadline)          # ssh: HIVE reads what is only on HIVE


def do_merge(a, deadline):
    """The HIVE end of the SSH route: point the uploaded copies at the upload, read what is only
    readable here, and write."""
    b = load_json(os.path.join(a.bundle, "bundle.json")) or {}
    plan = b.get("plan")
    if not plan:
        return not_recorded("bad_input", f"no bundle.json in {a.bundle}")
    # The gate again, on the side that writes: a HIVE account outside proteomics-grp that has a
    # working SSH login must still not get a record in.
    r = route(remote_hop=True)
    if r != "direct":
        return not_recorded(REASONS.get(r, r), where_text(r).split(": ", 1)[1], route=r)
    for k, v in (b.get("args") or {}).items():
        if getattr(a, k, None) in (None, False) and v is not None:
            setattr(a, k, v)
    for c in plan["copies"]:
        up = os.path.join(a.bundle, "files", c["rel"])
        if os.path.isfile(up):
            c["here"] = up
    z = plan.get("zip") or {}
    if z.get("copy") and os.path.isfile(os.path.join(a.bundle, "files", z.get("dest_rel", ""))):
        z["here"] = os.path.join(a.bundle, "files", z["dest_rel"])
    out, session = plan.get("out"), plan.get("session")
    if not plan.get("search"):
        if not out and session and os.path.isdir(session):
            out = session_search_out(session, None)
        if out and os.path.isdir(out):
            plan_search(plan, out, a, deadline)
    if not plan.get("analysis") and plan.get("event") == "analysis-done" \
            and session and os.path.isdir(session):
        plan_analysis(plan, session, a, plan.get("zip_cap") or int(a.zip_cap_gb * GB))
    if session and os.path.isdir(session):
        plan["identity"]["session"] = plan["identity"].get("session") or os.path.realpath(session)
    if not plan.get("prot"):
        plan["prot"] = find_prot(None, session if session and os.path.isdir(session) else None,
                                 plan.get("out") if plan.get("out") and os.path.isdir(
                                     plan["out"]) else None)
    plan["name"] = plan.get("name") or folder_name(plan, session, a.name)
    plan["by"] = f"{plan.get('by')} via SSH as {whoami()}"
    return execute(plan, a, deadline, "ssh")


def collect_rows(root):
    rows = []
    for p in sorted(glob.glob(os.path.join(root, "sessions", "*", "run_record.json"))):
        r, state = read_json_retry(p)
        if state != "ok":
            continue                      # reported on stderr; never touched here
        r = r or {}
        s, an = r.get("search") or {}, r.get("analysis") or {}
        d, res, fa = s.get("data") or {}, s.get("results") or {}, s.get("fasta") or {}
        sig = (an.get("de") or {}).get("significant_per_contrast") or {}
        rows.append({
            "date": r.get("date"), "session": r.get("name"), "user": r.get("user"),
            "prot": (r.get("prot") or {}).get("prot"), "status": s.get("status"),
            "exit_code": s.get("exit_code"), "engine": s.get("engine"),
            "version": s.get("engine_version"), "n_files": d.get("n_files"),
            "instrument": d.get("instrument"), "acquisition": d.get("acquisition"),
            "median_precursors": res.get("median_precursors"),
            "median_proteins": res.get("median_proteins"), "organism": fa.get("organism"),
            "finalized": bool(an),
            "significant": ";".join(f"{k}={v}" for k, v in sig.items()) or None,
            "fran": (s.get("fran") or {}).get("status"),
            "n_data_quality_notes": len(r.get("data_quality_notes") or []),
            "folder": os.path.dirname(p), "search_out": s.get("out_dir")})
    return rows


def do_list(a, deadline):
    r = route(remote_hop=a.remote_hop, read_only=True)
    if r == "direct":
        rows = collect_rows(runs_dir())
    elif r == "ssh" and _FILE:
        env = " ".join(f"{k}={shlex.quote(os.environ[k])}" for k in ("SKILL_RUNS_DIR",)
                       if os.environ.get(k))
        try:
            with open(_FILE) as fh:        # the list needs nothing but this file on HIVE
                p = subprocess.run(["bash", hive_exec_path(),
                                    f"{env + ' ' if env else ''}python3 - list --json --remote-hop"],
                                   stdin=fh, capture_output=True, text=True,
                                   timeout=max(10, deadline.left()))
        except (OSError, subprocess.SubprocessError) as e:
            print(f"record_run: could not list over SSH -- {e}")
            return
        rows = json_tail(p.stdout, "[")
        if rows is None:
            print(f"record_run: could not list over SSH -- {(p.stdout or p.stderr).strip()[-200:]}")
            return
    else:
        print(f"record_run: cannot list -- {where_text(r).split(': ', 1)[1]}")
        return
    if a.since:
        rows = [x for x in rows if (x.get("date") or "") >= a.since]
    if a.user:
        rows = [x for x in rows if x.get("user") in a.user]
    if a.status:
        rows = [x for x in rows if x.get("status") in a.status]
    if a.json:
        print(json.dumps(rows, indent=2, default=str))
        return
    cols = ["date", "session", "user", "prot", "status", "exit_code", "engine", "version",
            "n_files", "instrument", "acquisition", "median_precursors", "median_proteins",
            "organism", "finalized", "significant", "fran", "n_data_quality_notes", "folder",
            "search_out"]
    if a.tsv:
        print("\t".join(cols))
        for x in rows:
            print("\t".join("" if x.get(c) is None else str(x.get(c)) for c in cols))
        return
    show = ["date", "session", "user", "prot", "status", "engine", "version", "n_files",
            "median_precursors", "median_proteins", "organism"]
    table = [["" if x.get(c) is None else fmt_n(x[c]) for c in show] for x in rows]
    w = [max([len(c)] + [len(t[i]) for t in table]) for i, c in enumerate(show)]
    print("  ".join(c.ljust(w[i]) for i, c in enumerate(show)))
    for t in table:
        print("  ".join(v.ljust(w[i]) for i, v in enumerate(t)))
    print(f"({len(rows)} run(s) in {runs_dir() if r == 'direct' else 'the registry on HIVE'})")


def build_parser():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--where", action="store_true",
                    help="print where a record would go from here, and nothing else")
    sub = ap.add_subparsers(dest="cmd", metavar="{search-done,analysis-done,list}")

    def common(p):
        p.add_argument("--dry-run", action="store_true",
                       help="write nothing; the SEARCH_LOG.md it would write goes to stderr")
        p.add_argument("--prot", action="append", default=[],
                       help="the CoreOmics submission: PROT_0807 / 0807 / 807 and/or its 12-hex "
                            "id (repeatable)")
        p.add_argument("--name", default=None,
                       help="folder name for a NEW record (default: the session's name, else "
                            "<date>_<search dir name>)")
        p.add_argument("--issues-tag", action="append", default=[],
                       help="the --session tag given to report_issue.sh, when it was not the "
                            "session name (repeatable)")
        p.add_argument("--timeout", type=float, default=env_num("RECORD_RUN_TIMEOUT", 45),
                       help="give up after this many seconds (default 45; the job hook allows "
                            "60)")
        p.add_argument("--file-cap-mb", type=float, default=env_num("RECORD_RUN_FILE_CAP_MB", 20),
                       help="copy no single file larger than this (default 20 MB)")
        p.add_argument("--zip-cap-gb", type=float, default=env_num("RECORD_RUN_ZIP_CAP_GB", 5),
                       help="copy the session zip only below this size (default 5 GB)")
        p.add_argument("--remote-hop", action="store_true", help=argparse.SUPPRESS)
        p.add_argument("--skill-version", default=None, help=argparse.SUPPRESS)
    s = sub.add_parser("search-done", help="record a finished (or failed) search")
    s.add_argument("--out", help="the search output directory")
    s.add_argument("--status", choices=["completed", "failed"],
                   help="what the job that ran it saw (default: inferred from sacct + the report)")
    s.add_argument("--exit-code", default=None, help="the search's exit code")
    s.add_argument("--session", default=None,
                   help="the session dir, when the out dir is <session>/output/search")
    s.add_argument("--step", default=None,
                   help="the step that reported (default: $SLURM_JOB_NAME inside a job)")
    s.add_argument("--detect-json", default=None,
                   help="detect_acquisition.py output, when it is not beside the search")
    s.set_defaults(zip=None)
    common(s)
    f = sub.add_parser("analysis-done", help="add the finalized analysis to its search's record")
    f.add_argument("--session", help="the session dir session.py finalize ran on")
    f.add_argument("--out", default=None, help="the search output dir (default: found from the "
                                               "session)")
    f.add_argument("--zip", default=None, help="the session zip (default: <session>.zip)")
    f.add_argument("--detect-json", default=None, help=argparse.SUPPRESS)
    f.set_defaults(status=None, exit_code=None, step=None)
    common(f)
    m = sub.add_parser("merge", help=argparse.SUPPRESS)
    m.add_argument("--bundle", required=True)
    m.set_defaults(out=None, session=None, zip=None, status=None, exit_code=None, step=None,
                   detect_json=None)
    common(m)
    ls = sub.add_parser("list", help="list recorded runs (reads every run_record.json)")
    ls.add_argument("--since", help="YYYY-MM-DD: runs on or after this date")
    ls.add_argument("--user", action="append", help="only this HIVE user (repeatable)")
    ls.add_argument("--status", action="append",
                    help="completed | failed | running | unknown (repeatable)")
    ls.add_argument("--tsv", action="store_true", help="tab-separated, every column")
    ls.add_argument("--json", action="store_true")
    ls.add_argument("--timeout", type=float, default=120)
    ls.add_argument("--remote-hop", action="store_true", help=argparse.SUPPRESS)
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)          # a usage error exits 2 -- the one non-zero exit
    if a.where:
        print(where_text(route()))
        return 0
    if not a.cmd:
        build_parser().print_usage(sys.stderr)
        print(json.dumps(not_recorded("bad_input", "no command given")))
        return 0
    if a.cmd != "list" and disabled():
        # The kill switch: before any read, write or SSH attempt.
        print(json.dumps(not_recorded("disabled", "RECORD_RUN=off (or SKILL_RUNS_DIR=off)")))
        return 0
    deadline = Deadline(a.timeout)
    if hasattr(signal, "SIGALRM"):
        def _alarm(*_):
            raise Stop(f"stopped after {int(a.timeout) + 5} s")
        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(int(a.timeout) + 5)
    try:
        if a.cmd == "list":
            do_list(a, deadline)
            return 0
        if a.exit_code is not None:
            try:
                a.exit_code = int(a.exit_code)
            except ValueError:
                pass
        res = do_merge(a, deadline) if a.cmd == "merge" else do_event(a.cmd, a, deadline)
        if LOCK_TIMEOUTS:
            res["lock_timeout"] = sorted(set(res.get("lock_timeout") or []) | set(LOCK_TIMEOUTS))
        print(json.dumps(res, default=str))
    except KeyboardInterrupt:
        raise
    except BaseException as e:      # noqa: BLE001 -- never fatal, by contract
        reason = "timeout" if isinstance(e, (TimeoutError, Stop)) else "error"
        print(json.dumps(not_recorded(reason, f"{type(e).__name__}: {e}")))
    finally:
        if hasattr(signal, "SIGALRM"):
            signal.alarm(0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
