#!/usr/bin/env python3
"""
notify_slack.py  --  Tell the Core's Slack channel when a search ends and when an analysis is
finalized.

WHY THIS EXISTS
---------------
A search is hours long and the laptop that started it is usually closed by the time it ends,
so the watcher (watch_run.sh) that would have said "done" or "failed" is not running. The job
itself is. So the LAST job of every search route posts here from the compute node -- from
traps around the job's work (wrap_job_script), so a failure or a time limit posts too -- and
`session.py finalize` posts when the analysis is packaged. A job whose failure stops the rest of its chain (every `afterok` link) posts that
failure, because the chain then waits for ever and its last job never gets to say anything.

TWO RULES THAT MATTER MORE THAN THE FEATURE  (ported from STAN's stan/notify.py)
-----------------------------------------------------------------------------
1. **A notifier must never take down its caller.** It runs at the end of a multi-hour search
   and inside finalize. A dead webhook, a DNS blip or a Slack outage costs one log line, not a
   failed job: every public function returns and raises nothing, the POST is bounded (well
   under a minute even when DNS hangs), and the job trap never changes the job's exit status.
2. **The webhook URL is a bearer credential.** Anyone holding it can post into the channel. It
   is never printed, logged, put in an exception message, passed on a command line or written
   into a generated job script -- the job runs this script, and this script reads the URL
   itself. `_scrub` exists because urllib is perfectly willing to put the URL it was given
   into the text of an error.

WHERE THE WEBHOOK COMES FROM (first that is set)
------------------------------------------------
  $SKILL_SLACK_WEBHOOK
  ~/.config/ucdavis-proteomics/slack_webhook             (chmod 600)
  /quobyte/proteomics-grp/.config/skill_slack_webhook    (dir 2750, file 640, group
                                                          proteomics-grp: only Core members'
                                                          runs can read it, so only they post)
  none -> notifications are off: one info line says so and why, nothing else happens.
A value that is not a https://hooks.slack.com/ URL is treated as unconfigured, never posted to.
Off switch: SKILL_SLACK=0 (also false/no/off), or --no-notify on run_search.py / finalize.

When finalize runs on a laptop whose only webhook is the one on HIVE, the facts are handed to
this same script on HIVE through hive_exec.sh (`relay`), which posts from there.

    notify_slack.py search-done --out <search out dir> [--status ok|failed] [--exit-code N]
                                      # THE JOB-END HOOK: record_run.py -> fran_deposit.py
                                      # stage (success of the last job) -> Slack
    notify_slack.py analysis-done --session <session dir> [--zip <session>.zip]
    notify_slack.py --test            # one short test message; prints sent / not sent
    any of them + --dry-run           # print the JSON payload, send nothing
Setup, permissions and what gets posted: references/notifications.md.
"""
import argparse
import base64
import datetime
import getpass
import glob
import json
import os
import re
import shlex
import shutil
import socket
import statistics
import subprocess
import sys
import threading
import urllib.error
import urllib.parse
import urllib.request

# `python3 - relay` (the HIVE side of a relay) reads this file from stdin and has no __file__.
HERE = (os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else None)
if HERE and HERE not in sys.path:
    sys.path.insert(0, HERE)

WEBHOOK_ENV = "SKILL_SLACK_WEBHOOK"
OPT_OUT_ENV = "SKILL_SLACK"
USER_FILE = os.path.join("~", ".config", "ucdavis-proteomics", "slack_webhook")
# The Core's shared webhook. Env override for the tests only; Brett creates the real file.
GROUP_FILE = os.environ.get("SKILL_SLACK_GROUP_FILE",
                            "/quobyte/proteomics-grp/.config/skill_slack_webhook")
# Present only on HIVE. Its absence is how a laptop knows the group webhook is not local.
CORE_GROUP_DIR = os.environ.get("SKILL_CORE_GROUP_DIR", "/quobyte/proteomics-grp")

#: Slack's own host. Anything else is refused rather than posted to: a typo'd or planted
#: "webhook" would otherwise become an exfiltration channel for whatever the message holds.
_SLACK_HOST_PREFIX = "https://hooks.slack.com/"
#: Tests only: also accept a loopback http:// URL (a local mock server). Loopback cannot
#: carry anything off the machine, so this cannot become an exfiltration route.
_LOOPBACK_ENV = "SKILL_SLACK_TEST_LOOPBACK"
_LOOPBACK = re.compile(r"http://(127\.0\.0\.1|localhost)(:\d+)?/")

#: Long enough for a slow TLS handshake from HIVE, short enough that a job's exit does not
#: wait on an unresponsive Slack. urllib's timeout does not cover DNS, so the whole POST also
#: runs under a hard deadline (_DEADLINE_PAD more) on a daemon thread.
POST_TIMEOUT_SECONDS = 10
_DEADLINE_PAD = 5
#: ssh to HIVE (hive_exec.sh's ConnectTimeout is 20 s) plus the post on the far side.
RELAY_TIMEOUT_SECONDS = 90

_WEBHOOK_PATTERN = re.compile(r"(https://hooks\.slack\.com/|http://(127\.0\.0\.1|localhost)[:/])\S*")

ISSUE_LOOKBACK_DAYS = 14


def _issue_dirs():
    """Where report_issue.sh writes (its own defaults and overrides): the Core's folder on
    HIVE, and this machine's fallback."""
    return (os.environ.get("SKILL_ISSUES_DIR", "/quobyte/proteomics-grp/skill_issues"),
            os.path.expanduser(os.environ.get("SKILL_ISSUES_LOCAL_DIR",
                                              os.path.join("~", ".proteomics-pipeline",
                                                           "issues"))))


#: A TERM this close to the job's time limit is SLURM enforcing the limit. SLURM's timer
#: checks every 30-60 s, so the signal lands at or just after the limit.
_TIME_LIMIT_SLACK_S = 90


def say(msg):
    """The one line a call leaves in a log: secrets redacted, the webhook scrubbed, one line."""
    try:
        sys.stderr.write(f"[notify_slack] {clean_text(msg)}\n")
    except Exception:
        pass


#: What report_issue.sh refuses to write to a shared folder, plus a database DSN with a password
#: and the webhook itself. Redacted, not refused: an alert with "[redacted]" in it still tells the
#: channel something went wrong. The reviewer posted a webhook through send_alert().
_SECRET_PATTERNS = [
    re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----[\s\S]*?(?:-----END [A-Z ]*PRIVATE KEY-----|\Z)"),
    re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY[\s\S]*"),       # a header with no dashes after it
    re.compile(r"ghp_[A-Za-z0-9]{20,}"),
    re.compile(r"github_pat_[A-Za-z0-9_]+"),
    re.compile(r"hf_[A-Za-z0-9]{20,}"),
    re.compile(r"(?i)authorization:\s*(?:token|bearer)\s+\S+"),
    re.compile(r"(?i)password\s*[:=]\s*\S+"),
    re.compile(r"(?i)postgres(?:ql)?://[^\s:/@]+:[^\s@]+@"),
    re.compile(r"(?i)https?://hooks\.slack\.com/services/\S+"),
]


def redact(text):
    """Every secret-shaped substring of `text` replaced with [redacted]."""
    out = str(text)
    for pat in _SECRET_PATTERNS:
        out = pat.sub("[redacted]", out)
    return out


def one_line(text):
    """\r and \n flattened to spaces: a MANIFEST line or a status line stays one line."""
    return str(text).replace("\r", " ").replace("\n", " ")


def clean_text(text, webhook=None):
    """What may leave this module as a log/status/MANIFEST line: redacted, scrubbed, one line."""
    return one_line(_scrub(redact(text), webhook))


def _clean_payload(obj):
    """redact + _scrub every string in a Slack payload -- the one choke point for what is posted,
    whatever field (session name, stage error tail, alert body) it came in through."""
    if isinstance(obj, str):
        return _scrub(redact(obj))
    if isinstance(obj, list):
        return [_clean_payload(x) for x in obj]
    if isinstance(obj, dict):
        return {k: _clean_payload(v) for k, v in obj.items()}
    return obj


def _bounded(fn, seconds, on_timeout):
    """fn() on a daemon thread, abandoned after `seconds`: a stalled quobyte mount under an open()
    or a glob must not hold the caller. Never raises."""
    box = {}

    def work():
        try:
            box["r"] = fn()
        except Exception as e:  # noqa: BLE001
            box["r"] = (False, clean_text(f"not sent: {type(e).__name__}: {e}"))

    try:
        t = threading.Thread(target=work, daemon=True)
        t.start()
        t.join(seconds)
        if t.is_alive():
            return on_timeout
        return box.get("r", on_timeout)
    except Exception as e:  # noqa: BLE001
        return (False, clean_text(f"not sent: {type(e).__name__}: {e}"))


# ── webhook resolution ───────────────────────────────────────────


def opted_out():
    """The reason notifications are switched off by the user, or None."""
    v = (os.environ.get(OPT_OUT_ENV) or "").strip().lower()
    if v in ("0", "false", "no", "off"):
        return f"{OPT_OUT_ENV}={os.environ.get(OPT_OUT_ENV)}"
    return None


def _looks_like_webhook(url):
    if url.startswith(_SLACK_HOST_PREFIX):
        return True
    return os.environ.get(_LOOPBACK_ENV) == "1" and bool(_LOOPBACK.match(url))


def resolve_webhook():
    """(url, source) when a webhook is configured, else (None, why not). Never raises.

    `source` and `why` name WHERE a value came from or was looked for -- never the value."""
    env = (os.environ.get(WEBHOOK_ENV) or "").strip()
    if env:
        if _looks_like_webhook(env):
            return env, f"${WEBHOOK_ENV}"
        return None, f"${WEBHOOK_ENV} is set but is not a {_SLACK_HOST_PREFIX} URL"
    looked = [f"${WEBHOOK_ENV} unset"]
    for label, path in (("~/.config/ucdavis-proteomics/slack_webhook",
                         os.path.expanduser(USER_FILE)),
                        (GROUP_FILE, GROUP_FILE)):
        try:
            with open(path) as fh:
                url = fh.read().strip()
        except FileNotFoundError:
            looked.append(f"no {label}")
            continue
        except PermissionError:
            looked.append(f"{label} not readable by {_who()} (only proteomics-grp members post)")
            continue
        except OSError as e:
            looked.append(f"{label} unreadable ({type(e).__name__})")
            continue
        if not url:
            looked.append(f"{label} is empty")
            continue
        if not _looks_like_webhook(url):
            # Stop here, like STAN: a half-finished config leaves alerts off, it does not fall
            # through to a different channel.
            return None, f"{label} does not hold a {_SLACK_HOST_PREFIX} URL"
        return url, label
    return None, "no webhook configured (" + "; ".join(looked) + ")"


def _scrub(text, webhook=None):
    """Remove any webhook URL from text destined for a log or a caller -- the whole URL, and
    its path and token on their own, since an error body can echo back just the path."""
    out = str(text)
    if webhook:
        out = out.replace(webhook, "<webhook>")
        try:
            path = urllib.parse.urlsplit(webhook).path
        except Exception:
            path = ""
        if len(path) > 8:
            out = out.replace(path, "<webhook>")
        for seg in path.split("/"):
            if len(seg) >= 12:                  # Slack's token segment is 24 characters
                out = out.replace(seg, "<token>")
    return _WEBHOOK_PATTERN.sub("<webhook>", out)


# ── posting ──────────────────────────────────────────────────────


def _post_once(payload, hook, timeout):
    try:
        data = json.dumps(payload).encode("utf-8")
        req = urllib.request.Request(hook, data=data,
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            if 200 <= resp.status < 300:
                return True, "sent"
            return False, f"Slack returned HTTP {resp.status}"
    except urllib.error.HTTPError as e:
        # A revoked or mistyped webhook answers 403/404 with a one-word body ("no_service",
        # "invalid_token"). Worth having; the URL in e.url is not.
        body = ""
        try:
            body = e.read().decode("utf-8", "replace")[:120]
        except Exception:
            pass
        return False, _scrub(f"Slack rejected the message: HTTP {e.code} {body}".strip(), hook)
    except Exception as e:  # noqa: BLE001 -- a notifier never breaks its caller
        return False, _scrub(f"Slack post failed: {type(e).__name__}: {e}", hook)


def post(payload, hook, timeout=POST_TIMEOUT_SECONDS):
    """POST one message. (landed?, scrubbed detail); never raises, never blocks past the deadline."""
    box = {}

    def work():
        box["r"] = _post_once(payload, hook, timeout)

    try:
        t = threading.Thread(target=work, daemon=True)
        t.start()
        t.join(timeout + _DEADLINE_PAD)
        if t.is_alive():
            return False, f"no answer from Slack within {timeout + _DEADLINE_PAD} s"
        return box.get("r", (False, "Slack post failed"))
    except Exception as e:  # noqa: BLE001
        return False, _scrub(f"Slack post failed: {type(e).__name__}: {e}", hook)


# ── facts ────────────────────────────────────────────────────────


def _load(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError, TypeError):
        return None


def _who():
    """The HIVE user when there is one (the job's owner, or the saved HIVE login), else the
    local login."""
    try:
        return os.environ.get("SLURM_JOB_USER") or _hive_env_user() or getpass.getuser()
    except Exception:
        return "unknown"


def _skill_version():
    if not HERE:
        return None
    pj = _load(os.path.join(HERE, "..", ".claude-plugin", "plugin.json")) or {}
    return pj.get("version")


def _engine_label(engine):
    try:
        from make_methods import ENGINE_LABEL      # one definition of the display names
        return ENGINE_LABEL.get(engine, engine)
    except Exception:
        return engine


def _session_of(out):
    """The session directory a search output belongs to (<session>/output/search), or None."""
    parent = os.path.dirname(out)
    if os.path.basename(out) == "search" and os.path.basename(parent) == "output":
        return os.path.dirname(parent)
    return None


def _mtime_date(path):
    try:
        return datetime.date.fromtimestamp(os.path.getmtime(path)).isoformat()
    except OSError:
        return None


def _age_s(path):
    try:
        return max(0, int(datetime.datetime.now().timestamp() - os.path.getmtime(path)))
    except OSError:
        return None


def _run_stats(report):
    """Per-run headline numbers from DIA-NN's <report>.stats.tsv, or None."""
    try:
        import check_report_runs as crr            # one definition of where the stats file is
        rows = crr.stats_rows(crr.stats_path(report))
        if not rows:
            return None

        def num(r, col):
            try:
                return float(r.get(col) or 0)
            except ValueError:
                return 0.0

        prec = [num(r, "Precursors.Identified") for r in rows.values()]
        prot = [num(r, "Proteins.Identified") for r in rows.values()]
        return {"runs": len(rows), "median_precursors": int(statistics.median(prec)),
                "median_proteins": int(statistics.median(prot)),
                "runs_with_no_ids": sum(1 for p in prec if p <= 0)}
    except Exception:
        return None


def classify(exit_code, signal=None, job_seconds=None, time_limit_min=None):
    """(status, why) for how a job ended: status is ok, failed or stopped."""
    if signal:
        if time_limit_min and job_seconds is not None \
                and job_seconds >= time_limit_min * 60 - _TIME_LIMIT_SLACK_S:
            return "failed", f"hit its {_dur(time_limit_min * 60)} time limit"
        # scancel, or preemption on a requeue queue: both stop a job before its limit. A
        # cancel is deliberate, and a preempted job is requeued and reports when it finishes,
        # so neither is posted.
        return "stopped", (f"stopped by SIG{signal} before its time limit (cancelled, or "
                           "preempted and requeued)")
    if exit_code in (None, 0):
        return "ok", None
    if exit_code == 137:
        return "failed", "killed (exit 137, SIGKILL) -- most often out of memory"
    return "failed", f"exit {exit_code}"


def search_facts(out, exit_code=None, status=None, signal=None, started=None,
                 time_limit_min=None, stage=None, final=True):
    out = os.path.abspath(out)
    prov_path = os.path.join(out, "search_provenance.json")
    prov = _load(prov_path) or {}
    res = prov.get("result") if isinstance(prov.get("result"), dict) else {}
    session = _session_of(out)
    bundle = _load(prov.get("bundle") or "") if prov.get("bundle") else None
    if bundle is None and session:
        bundle = _load(os.path.join(session, "input", "wf", "workflow.manifest.json"))
    bundle = bundle or {}
    job_s = None
    if started:
        try:
            job_s = max(0, int(datetime.datetime.now().timestamp()) - int(started))
        except (TypeError, ValueError):
            job_s = None
    st, why = classify(exit_code, signal, job_s, time_limit_min)
    if status in ("ok", "failed") and not signal:
        st = status
        if status == "failed" and not why:
            why = "reported failed"
    engine = prov.get("engine") or (bundle.get("engine") or {}).get("name")
    mode = res.get("mode") or prov.get("search_mode")
    report = res.get("report") or os.path.join(out, "report.parquet")
    f = {
        "kind": "search", "status": st, "why": why, "final": bool(final), "stage": stage,
        "who": _who(), "host": socket.gethostname().split(".")[0],
        "session": os.path.basename(session) if session else None,
        "folder": None if session else "/".join(out.rstrip("/").split("/")[-2:]),
        "acquisition": bundle.get("acquisition"),
        "instrument": ((bundle.get("instruments") or [None])[0]
                       or bundle.get("instrument_label")),
        "runs": prov.get("n_files"),
        "engine": _engine_label(engine) if engine else None,
        "engine_version": prov.get("version"),
        "mode": mode,
        "job_seconds": job_s,
        "since_setup_seconds": _age_s(prov_path),
        "since": _mtime_date(prov_path) or datetime.date.today().isoformat(),
        "stats": _run_stats(report) if st == "ok" and final else None,
        "location": out,
        "log": _find_log(out),
        "job": {k: os.environ.get(v) for k, v in (("id", "SLURM_JOB_ID"),
                                                    ("name", "SLURM_JOB_NAME"),
                                                    ("array_id", "SLURM_ARRAY_JOB_ID"),
                                                    ("task", "SLURM_ARRAY_TASK_ID"))
                if os.environ.get(v)},
        "skill_version": _skill_version(),
    }
    return f


def _find_log(out):
    """This job's SLURM log, when it can be named without guessing."""
    jid, aid, task = (os.environ.get(k) for k in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID",
                                                   "SLURM_ARRAY_TASK_ID"))
    pats = ([f"*_{aid}_{task}.log"] if aid and task else []) + ([f"*_{jid}.log"] if jid else [])
    for pat in pats:
        hits = glob.glob(os.path.join(out, pat)) + glob.glob(os.path.join(out, "*", pat))
        if len(hits) == 1:
            return hits[0]
    return None


def analysis_facts(session_dir, zip_path=None):
    import session as sess                          # one definition of the session layout
    p = sess.paths_for(session_dir)
    sprov = _load(p["search_prov"]) or {}
    wf = _load(p["workflow_manifest"]) or {}
    de = _load(os.path.join(p["de_dir"], "de_provenance.json")) or {}
    raws = sess.read_raw_list(p["session_dir"])
    skipped = 0
    try:
        with open(p["manifest_txt"]) as fh:
            skipped = sum(1 for ln in fh if ln.startswith("[SKIPPED]"))
    except OSError:
        pass
    engine = sprov.get("engine") or (wf.get("engine") or {}).get("name")
    howto = os.path.join(p["deposit_dir"], "HOW_TO_SUBMIT.md")
    sig = de.get("significant_per_contrast")
    return {
        "kind": "analysis", "status": "ok", "final": True,
        "who": _who(), "host": socket.gethostname().split(".")[0],
        "on_hive": os.path.isdir(CORE_GROUP_DIR),
        "session": os.path.basename(p["session_dir"]),
        "acquisition": wf.get("acquisition"),
        "instrument": (wf.get("instruments") or [None])[0] or wf.get("instrument_label"),
        "runs": len(raws) or sprov.get("n_files"),
        "samples": de.get("n_samples"),
        "engine": _engine_label(engine) if engine else None,
        "engine_version": sprov.get("version") or (wf.get("engine") or {}).get("version"),
        "de_method": de.get("display_label") or de.get("method"),
        # run_de.R's own words for what "significant" meant (CLAUDE.md rule 1): counted there,
        # never re-derived here from the tables.
        "significant": sig if isinstance(sig, dict) else None,
        "significance_rule": de.get("significance_rule"),
        "adjp": de.get("adjp"),
        "since_setup_seconds": _age_s(p["raw_list"]),
        "since": _mtime_date(p["raw_list"]) or datetime.date.today().isoformat(),
        "location": p["session_dir"],
        "zip": zip_path if zip_path and os.path.isfile(zip_path) else None,
        "deposit": howto if os.path.isfile(howto) else None,
        "export_skipped": skipped,
        "skill_version": _skill_version(),
    }


def _clean(s):
    """report_issue.sh's clean(): the file-name form of a user or session name."""
    return re.sub(r"[^A-Za-z0-9._-]", "_", s)[:40]


def add_issues(facts):
    """Count the skill problems report_issue.sh recorded for this user since the run began.

    Its files are <date>_<user>[_<session>].md, one entry per `## ` heading. Counted where
    they are readable from here -- the Core's folder on HIVE, or this machine's fallback."""
    try:
        since = datetime.date.fromisoformat(facts.get("since") or "")
    except Exception:                               # a bad date must never cost the post
        since = datetime.date.today()
    files, n = [], 0
    try:
        today = datetime.date.today()
        since = max(since, today - datetime.timedelta(days=ISSUE_LOOKBACK_DAYS))
        user = _clean(facts.get("who") or "")
        for d in _issue_dirs():
            day = since
            while day <= today:
                for fpath in sorted(glob.glob(os.path.join(d, f"{day.isoformat()}_{user}*.md"))):
                    try:
                        with open(fpath) as fh:
                            k = sum(1 for ln in fh if ln.startswith("## "))
                    except OSError:
                        continue
                    if k:
                        n += k
                        files.append(fpath)
                day += datetime.timedelta(days=1)
    except Exception:
        return facts
    if n:
        facts["issues"] = {"entries": n, "files": files[:3], "more_files": max(0, len(files) - 3)}
    return facts


# ── rendering ────────────────────────────────────────────────────


def _dur(s):
    if s is None:
        return None
    s = int(s)
    d, rem = divmod(s, 86400)
    h, rem = divmod(rem, 3600)
    m = rem // 60
    if d:
        return f"{d}d {h}h"
    if h:
        return f"{h}h {m:02d}m"
    return f"{m}m" if m else f"{s}s"


def _esc(s):
    return str(s).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def _cut(s, n):
    s = str(s)
    return s if len(s) <= n else s[:n - 1] + "…"


_MODE_LABEL = {"parallel_5step": "5-step SLURM chain", "two_job_libfree": "library + search jobs",
               "single_shot": "one job", "radiant_parallel_3step": "3-step SLURM chain"}

_ICON = {"ok": ":white_check_mark:", "failed": ":x:", "stopped": ":double_vertical_bar:",
         "test": ":wave:"}


def _field(title, value):
    return {"type": "mrkdwn", "text": _cut(f"*{title}*\n{value}", 1900)}


def render(f):
    """facts -> Slack payload: `text` is what the phone notification shows, so it carries the
    whole headline; `blocks` are what the channel renders."""
    if f.get("kind") == "alert":
        head = f":warning: {_esc(f['title'])}: " if f.get("title") else ":warning: "
        text = _cut(head + _esc(f.get("body") or ""), 300)
        ctx = " · ".join(x for x in (
            f"ucdavis-proteomics-core-pipeline v{f['skill_version']}" if f.get("skill_version")
            else "ucdavis-proteomics-core-pipeline", f"{_esc(f.get('who'))} on {_esc(f.get('host'))}")
            if x)
        return {"text": text,
                "blocks": [{"type": "section", "text": {"type": "mrkdwn",
                                                        "text": _cut(head + _esc(f.get("body") or ""), 2900)}},
                           {"type": "context", "elements": [{"type": "mrkdwn", "text": ctx}]}]}
    if f.get("kind") == "test":
        text = _esc(f"{_ICON['test']} UC Davis proteomics pipeline: Slack notifications work "
                    f"(test sent by {f.get('who')} from {f.get('host')})")
        return {"text": text, "blocks": [{"type": "section",
                                          "text": {"type": "mrkdwn", "text": text}}]}
    name = f.get("session") or f.get("folder") or "(unnamed)"
    eng = " ".join(str(x) for x in (f.get("engine"), f.get("engine_version")) if x) or None
    runs = f.get("runs")
    if f["kind"] == "search":
        if f["status"] == "ok":
            head = "Search finished"
        else:
            head = "Search FAILED" + (f" at {f['stage']}" if f.get("stage") and not
                                      f.get("final") else "")
    else:
        head = "Analysis complete"
    icon = _ICON.get(f["status"], ":information_source:")
    bits = [x for x in (eng, f"{runs} runs" if runs else None) if x]
    st = f.get("stats")
    if st:
        bits.append(f"median {st['median_precursors']:,} precursors/run")
    if f["kind"] == "analysis" and f.get("significant"):
        bits.append(f"{len(f['significant'])} contrast(s)")
    if f.get("why") and f["status"] != "ok":
        bits.append(f["why"])
    # Escaped like the blocks: `<!channel>` or `<url|x>` in a session name must not ping or link.
    text = _cut(f"{icon} {head}: {_esc(name)}" + (f" ({_esc(', '.join(str(b) for b in bits))})"
                                                   if bits else "")
                + f" -- {_esc(f.get('who'))}", 300)

    fields = [_field("Who", _esc(f.get("who"))),
              _field("Instrument", _esc(" · ".join(x for x in (f.get("instrument"),
                                                                f.get("acquisition")) if x)
                                        or "not recorded"))]
    if f["kind"] == "analysis" and f.get("samples"):
        fields.append(_field("Runs / samples in DE", f"{runs or '?'} / {f['samples']}"))
    else:
        fields.append(_field("Runs", runs or "not recorded"))
    fields.append(_field("Engine", _esc(eng or "not recorded")
                         + (f" ({_esc(_MODE_LABEL.get(f['mode'], f['mode']))})"
                            if f.get("mode") else "")))
    if f["kind"] == "search":
        stage = f" -- {f['stage']}" if f.get("stage") else ""
        status = ("OK" if f["status"] == "ok" else f"FAILED: {f.get('why')}") + stage
        fields.append(_field("Status", _esc(status)))
        el = [x for x in ((f"{_dur(f['since_setup_seconds'])} since set-up"
                           if f.get("since_setup_seconds") is not None else None),
                          (f"this job {_dur(f['job_seconds'])}"
                           if f.get("job_seconds") is not None else None)) if x]
        fields.append(_field("Elapsed", " · ".join(el) or "not recorded"))
        if st:
            per = f"{st['median_precursors']:,} precursors · {st['median_proteins']:,} proteins"
            if st.get("runs_with_no_ids"):
                per += f"\n:warning: {st['runs_with_no_ids']} run(s) identified nothing"
            fields.append(_field("Median per run", per))
    else:
        fields.append(_field("DE method", _esc(f.get("de_method") or "not recorded")))
        if f.get("since_setup_seconds") is not None:
            fields.append(_field("Elapsed", f"{_dur(f['since_setup_seconds'])} since the "
                                            "session was created"))
    blocks = [{"type": "section",
               "text": {"type": "mrkdwn", "text": _cut(f"{icon} *{head}* — `{_esc(name)}`", 2900)}},
              {"type": "section", "fields": fields[:10]}]

    if f["kind"] == "analysis":
        sig = f.get("significant")
        if sig:
            rule = f.get("significance_rule") or (f"adj.P.Val < {f['adjp']}" if f.get("adjp")
                                                  else "as recorded by run_de.R")
            if f.get("adjp") is not None and "adjp" in rule:
                rule += f"; adjp = {f['adjp']}"      # run_de.R's words, with its value
            lines = [f"• {_esc(k)}: *{v}*" for k, v in list(sig.items())[:12]]
            if len(sig) > 12:
                lines.append(f"… and {len(sig) - 12} more contrasts")
            blocks.append({"type": "section", "text": {"type": "mrkdwn", "text": _cut(
                f"*Significant proteins per contrast* ({_esc(rule)})\n" + "\n".join(lines),
                2900)}})
        else:
            blocks.append({"type": "section", "text": {"type": "mrkdwn", "text":
                           "*Significant proteins per contrast:* not recorded "
                           "(no de_provenance.json)"}})

    here = "" if f["kind"] == "search" or f.get("on_hive") else f" (on {_esc(f.get('host'))})"
    where = [f"*Results:* `{_esc(f['location'])}`{here}"]
    if f.get("zip"):
        where.append(f"*Zip:* `{_esc(f['zip'])}`")
    if f.get("deposit"):
        where.append(f"*Deposit (PRIDE/MassIVE):* `{_esc(f['deposit'])}`")
    if f.get("export_skipped"):
        where.append(f":warning: {f['export_skipped']} export part(s) skipped -- see MANIFEST.txt")
    if f.get("log") and f["status"] != "ok":
        where.append(f"*Log:* `{_esc(f['log'])}`")
    steps = []
    if f.get("run_log") is not None:
        steps.append(f"*Run log:* {_esc(_run_log_text(f['run_log']))}")
    if f.get("fran") is not None:
        steps.append(f"*Staged for FRAN:* {_esc(_fran_text(f['fran']))}")
    if steps:
        where.append(" · ".join(steps))
    if f["kind"] == "search" and f["status"] == "failed" and not f.get("final"):
        where.append("The jobs after this one wait on it (`afterok`) and will not start "
                     "until it is fixed and resubmitted -- see references/watcher.md.")
    blocks.append({"type": "section", "text": {"type": "mrkdwn",
                                               "text": _cut("\n".join(where), 2900)}})

    ctx = []
    iss = f.get("issues")
    if iss:
        files = ", ".join(f"`{_esc(x)}`" for x in iss["files"])
        more = f" (+{iss['more_files']} more)" if iss.get("more_files") else ""
        ctx.append(f":memo: {iss['entries']} skill issue(s) recorded by {_esc(f.get('who'))} "
                   f"since {_esc(f.get('since'))}: {files}{more}")
    job = f.get("job") or {}
    tail = [f"ucdavis-proteomics-core-pipeline v{f['skill_version']}" if f.get("skill_version")
            else "ucdavis-proteomics-core-pipeline"]
    if job.get("id"):
        jid = (f"{job['array_id']}_{job['task']}" if job.get("array_id") and job.get("task")
               else job["id"])
        tail.append(f"job {jid}" + (f" ({_esc(job['name'])})" if job.get("name") else "")
                    + f" on {_esc(f.get('host'))}")
    ctx.append(" · ".join(tail))
    blocks.append({"type": "context", "elements": [{"type": "mrkdwn", "text": _cut(c, 2900)}
                                                   for c in ctx]})
    return {"text": text, "blocks": blocks}


# ── delivery ─────────────────────────────────────────────────────


def _hive_env_user():
    """HIVE_USER from the environment or hive_exec.sh's saved config, or None."""
    if os.environ.get("HIVE_USER"):
        return os.environ["HIVE_USER"]
    cfg = os.environ.get("HIVE_ENV_FILE",
                         os.path.expanduser("~/.config/ucdavis-proteomics/hive.env"))
    try:
        with open(cfg) as fh:
            for ln in fh:
                m = re.match(r"\s*(?:export\s+)?HIVE_USER=['\"]?([^'\"\s]+)", ln)
                if m:
                    return m.group(1)
    except OSError:
        pass
    return None


def _relay_route():
    """hive_exec.sh when this machine can hand a message to HIVE, else None."""
    if os.environ.get("SKILL_SLACK_NO_RELAY") == "1" or os.path.isdir(CORE_GROUP_DIR):
        return None                               # on HIVE already: nothing to relay to
    hx = os.environ.get("HIVE_EXEC") or (os.path.join(HERE, "hive_exec.sh") if HERE else None)
    if not hx or not os.path.isfile(hx) or not HERE or not shutil.which("bash"):
        return None
    return hx if _hive_env_user() else None


def _relay(facts, hx):
    """Post from HIVE: this same script, sent over hive_exec.sh's stdin, with the facts."""
    try:
        with open(os.path.join(HERE, "notify_slack.py")) as fh:
            src = fh.read()
        b64 = base64.b64encode(json.dumps(facts).encode("utf-8")).decode("ascii")
        r = subprocess.run(["bash", hx, f"python3 - relay --facts-b64 {b64}"], input=src,
                           capture_output=True, text=True, timeout=RELAY_TIMEOUT_SECONDS)
        last = (r.stdout.strip().splitlines() or [""])[-1].strip()
        if last == "sent":
            return True, f"sent through HIVE ({_hive_env_user()})"
        if last.startswith("not sent"):
            return False, clean_text(last + " (on HIVE)")
        err = (r.stderr.strip().splitlines() or [f"exit {r.returncode}"])[-1]
        return False, clean_text(f"not sent: could not reach HIVE to post: {err}")
    except subprocess.TimeoutExpired:
        return False, f"not sent: no answer from HIVE within {RELAY_TIMEOUT_SECONDS} s"
    except Exception as e:  # noqa: BLE001
        return False, clean_text(f"not sent: relay through HIVE failed: {type(e).__name__}: {e}")


#: The public words for "no webhook here": no path, no group name -- a collaborator's MANIFEST
#: and job log must not describe the Core's internals. The detail is what --dry-run is for.
PUBLIC_OFF = "Core notification not configured for this user"
#: deliver()'s whole budget, webhook lookup and issue count included (a stalled quobyte mount
#: under an open() used to sit outside the POST's own deadline).
DELIVER_DEADLINE_S = POST_TIMEOUT_SECONDS + _DEADLINE_PAD + 5


def payload(facts):
    """render() with every string redacted and scrubbed -- what is actually posted."""
    return _clean_payload(render(facts))


def _deliver(facts, relay_ok, dry_run):
    off = opted_out()
    hook, where = (None, None) if off else resolve_webhook()
    hx = None if (off or hook or not relay_ok) else _relay_route()
    if not hx and facts.get("kind") in ("search", "analysis"):
        add_issues(facts)                       # alerts and tests render no issue count
    if dry_run:
        print(json.dumps(payload(facts), indent=2))
        route = (f"off: {off}" if off else f"would send via {where}" if hook else
                 f"would relay through HIVE ({_hive_env_user()}) via hive_exec.sh" if hx
                 else f"off: {where}")
        return False, clean_text(f"dry run, nothing sent ({route})")
    if off:
        return False, f"off: {off}"
    if hook:
        ok, detail = post(payload(facts), hook)
        return ok, ("sent" if ok else clean_text(f"not sent: {detail}", hook))
    if hx:
        return _relay(facts, hx)
    return False, f"off: {PUBLIC_OFF}"


def deliver(facts, relay_ok=False, dry_run=False):
    """Send one message. (sent?, one-line status) -- never raises, never prints the URL, and
    never runs past its deadline, webhook lookup included.

    Local webhook first; with none, and `relay_ok`, hand the facts to HIVE through
    hive_exec.sh, where the Core's webhook lives. `dry_run` prints the payload and sends
    nothing (and names where it would have gone)."""
    limit = DELIVER_DEADLINE_S + (RELAY_TIMEOUT_SECONDS if relay_ok else 0)
    return _bounded(lambda: _deliver(facts, relay_ok, dry_run), limit,
                    (False, f"not sent: no answer within {limit} s"))


# ── the other job-end steps: the central run log and the FRAN hand-over ──

#: Both helpers are the other scripts' own CLIs, run as subprocesses, bounded, and never able to
#: change how the job ended. record_run.py appends one line to the Core's run log; fran_deposit.py
#: `stage` links a finished search into FRAN's drop directory (symlinks, no database, no
#: credential) after its own eligibility check.
RECORD_RUN_TIMEOUT_S = 60
FRAN_STAGE_TIMEOUT_S = 300


def _helper(name):
    """A sibling skill script, or None when this install does not have it (record_run.py is
    newer than some installs; a relayed copy of this file has no siblings at all)."""
    if not HERE:
        return None
    p = os.path.join(HERE, name)
    return p if os.path.isfile(p) else None


def _run_helper(argv, timeout):
    """(exit code, stdout, stderr) -- exit code None when it hung or could not start. Never raises."""
    try:
        r = subprocess.run(argv, capture_output=True, text=True, timeout=timeout,
                           stdin=subprocess.DEVNULL)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return None, "", f"no answer within {timeout} s"
    except Exception as e:  # noqa: BLE001
        return None, "", f"{type(e).__name__}: {e}"


def _json_out(text):
    """The JSON object a helper printed -- all of stdout, else its last line -- or None."""
    t = (text or "").strip()
    for cand in (t, (t.splitlines() or [""])[-1]):
        try:
            d = json.loads(cand)
        except ValueError:
            continue
        if isinstance(d, dict):
            return d
    return None


def _last_line(*texts):
    for t in texts:
        lines = [ln.strip() for ln in (t or "").strip().splitlines() if ln.strip()]
        if lines:
            return _cut(_scrub(lines[-1]), 160)
    return ""


def record_run(kind, *, out=None, session=None, status=None, exit_code=None):
    """Log this run in the Core's central run log with record_run.py. None when this install
    has no record_run.py; otherwise {"logged": True|False|None, "detail", "path"} (None = it ran
    but did not say). Never raises; record_run.py always exits 0 by contract, so a non-zero exit
    or a hang is reported as an error, never as logged."""
    rr = _helper("record_run.py")
    if not rr:
        return None
    if (os.environ.get("RECORD_RUN") or "").strip().lower() in ("off", "0", "no", "false"):
        return {"logged": False, "off": True, "detail": "RECORD_RUN=off"}
    argv = [sys.executable or "python3", rr, kind]
    if kind == "search-done":
        argv += ["--out", out, "--status", status, "--exit-code", str(exit_code)]
    if session:
        argv += ["--session", session]
    rc, so, se = _run_helper(argv, RECORD_RUN_TIMEOUT_S)
    if rc != 0:
        why = "did not finish" if rc is None else f"exit {rc}"
        return {"logged": False, "error": True,
                "detail": f"record_run.py {why}: {_last_line(se, so) or 'no output'}"}
    d = _json_out(so)
    if d is None:
        return {"logged": None, "detail": _last_line(so) or None}
    logged = next((d[k] for k in ("recorded", "logged", "ok") if isinstance(d.get(k), bool)), None)
    return {"logged": logged, "detail": d.get("reason") or d.get("detail"),
            "path": d.get("path") or d.get("file")}


def _stage_argv(fd, out, *, name=None, qc=None, fasta_meta=None):
    """`fran_deposit.py stage` argv, in exactly the shape of fran_deposit.stage_argv() (the FRAN
    side's own builder): stage --out O [--fasta-meta M] [--name N] [--qc | --not-qc]. Built here
    rather than imported, so the hook never depends on importing fran_deposit.py; the flags
    match (tests/test_slack_notify.py pins the shape). --qc / --not-qc / --skip exist only in
    the fran_deposit.py that ships with them -- an older stage rejects them, which the hook
    reports as an `error` line and nothing more."""
    argv = [sys.executable or "python3", fd, "stage", "--out", out]
    if fasta_meta:
        argv += ["--fasta-meta", fasta_meta]
    if name and str(name).strip():
        argv += ["--name", str(name).strip()]
    return argv + ({True: ["--qc"], False: ["--not-qc"]}.get(qc, []))


def fran_stage(out, session=None, timeout=FRAN_STAGE_TIMEOUT_S, name=None, qc=None,
               skip=False):
    """Hand a finished search to FRAN: `fran_deposit.py stage --out <out>` (+ --fasta-meta when
    the session has one). stage decides eligibility itself -- Core member (drop dir writable),
    a DIA engine, a complete report, not already staged, FRAN_DEPOSIT not off. None when this
    install has no fran_deposit.py; otherwise {"staged", "reason", "detail", "entry"}, with
    reason "error" when stage crashed or hung. Never raises.
    name / qc: the session's name and the QC decision made at generation (qc True -> --qc,
    False -> --not-qc). skip: --no-fran -- `stage --skip` stages nothing but records
    `opted_out` in <out>/fran_deposit.json, so a later `fran_deposit.py backfill` leaves it be."""
    fd = _helper("fran_deposit.py")
    if not fd:
        return None
    try:
        import session as sess                  # one definition of the session layout
        meta = sess.paths_for(session)["fasta_meta"] if session else None
    except Exception:
        meta = None
    argv = _stage_argv(fd, out, name=name, qc=qc,
                       fasta_meta=meta if meta and os.path.isfile(meta) else None)
    if skip:
        argv.append("--skip")
    rc, so, se = _run_helper(argv, timeout)
    d = _json_out(so)
    if rc != 0 or d is None:
        why = "did not finish" if rc is None else f"exit {rc}"
        return {"staged": False, "reason": "error",
                "detail": f"fran_deposit.py stage {why}: {_last_line(se, so) or 'no output'}"}
    # health_warning: stage's one line when FRAN's own ingest is unhealthy (stuck, not running,
    # stale code). The search is staged either way; the post says it will not be picked up yet.
    return {"staged": bool(d.get("staged")), "reason": d.get("reason"),
            "detail": _cut(d.get("detail") or "", 300) or None, "entry": d.get("entry"),
            "health_warning": _cut(d.get("health_warning") or "", 200) or None}


def _steps_summary(f):
    """'run log: yes; FRAN: staged' -- what the other job-end steps did, for the job log."""
    parts = []
    rl, fr = f.get("run_log"), f.get("fran")
    if rl is not None:
        parts.append(f"run log: {_run_log_text(rl)}")
    if fr is not None:
        parts.append(f"FRAN: {_fran_text(fr)}")
    return "; ".join(parts)


def _run_log_text(rl):
    if rl.get("logged") is True:
        return "yes"
    if rl.get("logged") is False:
        return ("error -- " if rl.get("error") else "no -- ") + str(rl.get("detail") or "no reason given")
    return "ran" + (f" -- {rl['detail']}" if rl.get("detail") else "")


def _fran_text(fr):
    if fr.get("staged"):
        return ("yes, but FRAN's ingest is unhealthy: " + fr["health_warning"]
                if fr.get("health_warning") else "yes")
    if fr.get("reason") == "error":
        return f"error: {fr.get('detail')}"
    words = {"opted_out": "off for this search",
             "left_to_agent": "left to the agent (no completeness guard on this route)",
             "near_time_limit": "left to the agent (too close to the job's time limit)",
             "qc_run": "instrument QC / standard run"}
    return f"skipped: {words.get(fr.get('reason'), fr.get('reason') or 'no reason given')}"


def slack_manifest(sent, detail):
    """(level, part, note) for finalize's MANIFEST.txt. [SKIPPED] only for a post that was
    ATTEMPTED and failed; not configured / opted out is [INFO] -- not a missing export part --
    and never names a path or a group."""
    d = str(detail or "")
    if sent:
        return "OK", "Slack notification (Core channel)", one_line(d or "sent")
    if d.startswith("--no-notify"):
        return "INFO", "Slack notification (Core channel)", "off (--no-notify)"
    m = re.search(r"off: (SKILL_SLACK=\S+)", d)
    if m:
        return "INFO", "Slack notification (Core channel)", f"off ({m.group(1)})"
    if d.startswith("off:") or d.startswith("not sent: off:"):
        return "INFO", "Core notification", "not configured for this user"
    return "SKIPPED", "Slack notification (Core channel)", clean_text(d)


def run_log_manifest(run_log):
    """(level, part, note) for the central run log's MANIFEST.txt line. Never a path: the zip is
    shared with collaborators."""
    part = "Core run log"
    if run_log is None:
        return "INFO", part, "not configured for this user"
    if run_log.get("off"):
        return "INFO", part, f"off ({run_log.get('detail')})"
    if run_log.get("logged") is True:
        return "OK", part, "logged"
    if run_log.get("error"):
        return "SKIPPED", part, "record_run.py failed -- the finalize output has the error"
    if run_log.get("logged") is False:
        return "INFO", part, "not configured for this user"
    return "INFO", part, "ran; record_run.py did not report whether it logged"


def _array_first_failure(out):
    """True for the first failing task of an array job (or any non-array job): one report per
    failed array, not one per task. The tasks run on DIFFERENT nodes, so the claim is an
    os.mkdir() in the (shared) search folder -- the primitive measured to be atomic across HIVE
    nodes on Quobyte (2026-09-24: flock lost 578/800 updates between 2 nodes, mkdir lost 0).
    FileExistsError = another task already reported."""
    aid = os.environ.get("SLURM_ARRAY_JOB_ID")
    if not aid:
        return True
    try:
        os.mkdir(os.path.join(out, f".slack_failed_{aid}"))
        return True
    except FileExistsError:
        return False
    except OSError:
        return True


#: Staging needs this much of the job's time left: stage takes seconds, but a job SIGKILLed
#: mid-stage leaves a half-linked drop entry. Step 7c's own check + stage picks it up instead.
NEAR_LIMIT_S = 120


def _seconds_left():
    """Seconds until this job's SLURM end time (`squeue -h -j $SLURM_JOB_ID -o %e`), or None when
    that cannot be read -- then the check is skipped, silently."""
    jid = os.environ.get("SLURM_JOB_ID")
    if not jid:
        return None
    rc, so, _ = _run_helper(["squeue", "-h", "-j", jid, "-o", "%e"], 10)
    if rc != 0:
        return None
    try:
        end = datetime.datetime.strptime(so.strip().splitlines()[0].strip(), "%Y-%m-%dT%H:%M:%S")
    except Exception:
        return None                              # N/A, Unknown, an array listing, anything else
    return int((end - datetime.datetime.now()).total_seconds())


def _fran_step(out, session, mode, fran_name=None, qc=None):
    """The FRAN part of the job-end hook, for a final job that succeeded.

    --qc and --no-fran never stage anything, on any route -- but they still CALL stage, which
    records the decision in <out>/fran_deposit.json, so a QC run or an opted-out search found
    later by `fran_deposit.py backfill` is left alone:
      --qc       stage --out O [--name N] --qc   (receipt: qc_run)
      --no-fran  stage --out O [--name N] --skip (receipt: opted_out)
    Every stage call -- recording or staging -- has the same bounds: not within NEAR_LIMIT_S of
    the job's time limit, and at most the time left minus 30 s."""
    record_only = qc is True or mode == "off"
    if not record_only and mode != "stage":
        say("fran: left_to_agent (no completeness guard on this route)")
        return {"staged": False, "reason": "left_to_agent",
                "detail": "no completeness guard on this route; step 7c stages it after its check"}
    left = _seconds_left()
    if left is not None and left < NEAR_LIMIT_S:
        if record_only:
            return {"staged": False, "reason": "qc_run" if qc is True else "opted_out",
                    "detail": f"{max(left, 0)} s left of the job's time limit: the decision was "
                              "not recorded in fran_deposit.json"}
        return {"staged": False, "reason": "near_time_limit",
                "detail": f"{max(left, 0)} s left of the job's time limit; step 7c stages it"}
    timeout = FRAN_STAGE_TIMEOUT_S if left is None else max(10, min(FRAN_STAGE_TIMEOUT_S, left - 30))
    return fran_stage(out, session, timeout=timeout, name=fran_name, qc=qc,
                      skip=record_only and qc is not True)


def search_done(out, exit_code=None, status=None, signal=None, started=None,
                time_limit_min=None, stage=None, final=True, from_job=False, dry_run=False,
                slack=True, fran="left_to_agent", fran_name=None, qc=None):
    """THE JOB-END HOOK: what a search job does when it ends, in this order --
      1. record_run.py search-done   (the run log; on success and on failure)
      2. fran_deposit.py stage        (the route's LAST job, on SUCCESS only)
      3. the Slack post               (unless `slack` is False / SKILL_SLACK=0)
    so the post can say what 1 and 2 did. Returns (posted?, one status line); never raises,
    and nothing here can change the job's exit status.

    Acted on: a final job's success; any job's failure (an `afterok` chain stops there, so its
    last job will never report); an array's first failing task only. Not acted on: a non-final
    job's success, and a job stopped before its time limit (cancelled, or preempted and
    requeued -- it runs again and reports then).

    `fran` is job_end_plan()'s decision, baked into the hook: "stage" (a guarded DIA-NN route,
    and not within NEAR_LIMIT_S of the time limit), "left_to_agent", or "off"."""
    try:
        if from_job and not os.environ.get("SLURM_JOB_ID"):
            return False, "nothing done: --from-job outside a SLURM job"
        f = search_facts(out, exit_code, status, signal, started, time_limit_min, stage, final)
        if f["status"] == "stopped":
            return False, f"nothing done: {f['why']}"
        if f["status"] == "ok" and not final:
            return False, "nothing to do: this job succeeded; the search's last job reports"
        if f["status"] == "failed" and not dry_run and not _array_first_failure(f["location"]):
            return False, "nothing done: another task of this array already reported the failure"
        if not dry_run:
            session = _session_of(f["location"])
            f["run_log"] = record_run("search-done", out=f["location"], session=session,
                                      status="completed" if f["status"] == "ok" else "failed",
                                      exit_code=143 if signal else (exit_code or 0))
            if f["status"] == "ok" and final:
                f["fran"] = _fran_step(f["location"], session, fran, fran_name, qc)
                if f["fran"] is not None and f["fran"].get("reason") != "left_to_agent":
                    say("fran_deposit stage: " + json.dumps(f["fran"], separators=(",", ":")))
        if slack:
            sent, msg = deliver(f, relay_ok=False, dry_run=dry_run)
        else:
            sent, msg = False, "Slack: off for this job (--no-notify when it was generated)"
        steps = _steps_summary(f)
        return sent, clean_text(f"{steps}; {msg}" if steps else msg)
    except Exception as e:  # noqa: BLE001
        return False, clean_text(f"job-end hook failed: {type(e).__name__}: {e}")


def analysis_done(session_dir, zip_path=None, dry_run=False, run_log=None):
    """Post that an analysis was finalized. (sent?, status line); never raises. `run_log` is
    what record_run.py said (session.py finalize runs it first), shown in the post."""
    try:
        f = analysis_facts(session_dir, zip_path)
        f["run_log"] = run_log
        return deliver(f, relay_ok=True, dry_run=dry_run)
    except Exception as e:  # noqa: BLE001
        return False, _scrub(f"not sent: {type(e).__name__}: {e}")


def send_alert(text, *, title=None, with_status=False):
    """Post a short free-text alert (one or two lines) to the Core's channel -- the API for
    other skill scripts, e.g. `fran_deposit.py health --alert`.

    Returns True when Slack accepted it, False otherwise: no webhook configured (off), opted
    out (SKILL_SLACK=0), an HTTP error, a timeout. With with_status=True it returns
    (sent, status line) instead -- the line never contains the webhook. Never raises. Prints
    NOTHING to stdout or stderr (a caller's stdout may be a JSON contract); log the status line
    yourself if you want one. Does not relay through hive_exec.sh: it runs where it is called.
    Bounded as a whole (~25 s), webhook lookup included. `text` and `title` are redacted
    ([redacted] for keys, tokens, passwords, DSNs, webhooks) before anything is posted."""
    def go():
        return deliver({"kind": "alert", "status": "alert",
                        "title": redact(title) if title else None, "body": redact(text),
                        "who": _who(), "host": socket.gethostname().split(".")[0],
                        "skill_version": _skill_version()}, relay_ok=False)

    # The facts (user, plugin.json) are read under the same deadline as the post: a stalled
    # mount anywhere must not hold `fran_deposit.py health --alert`.
    sent, status = _bounded(go, DELIVER_DEADLINE_S + 5,
                            (False, f"not sent: no answer within {DELIVER_DEADLINE_S + 5} s"))
    return (sent, status) if with_status else sent


def send_test(dry_run=False):
    return deliver({"kind": "test", "status": "test", "who": _who(),
                    "host": socket.gethostname().split(".")[0]},
                   relay_ok=True, dry_run=dry_run)


# ── the job-script wrapper ───────────────────────────────────────


def fran_opted_out():
    """FRAN_DEPOSIT=off (fran_deposit.py's own opt-out values) in THIS environment, or None."""
    v = (os.environ.get("FRAN_DEPOSIT") or "").strip().lower()
    return f"FRAN_DEPOSIT={os.environ.get('FRAN_DEPOSIT')}" if v in ("off", "0", "no", "false") \
        else None


def job_end_plan(*, slack=True, fran=True, fran_guarded=False, fran_name=None, qc=None):
    """What a generated job's end hook will do, decided NOW, at generation.

    In hive_remote every hive_exec.sh call is a fresh shell, so an opt-out in the generating
    environment (SKILL_SLACK=0, FRAN_DEPOSIT=off) never reaches the job -- it is baked into the
    hook as a flag, and recorded (run_search.py writes this into search_provenance.json).
      fran "stage"          the route has a completeness guard (DIA-NN: report_guard / step 5's
                            count), so its last job stages a finished search itself
      fran "left_to_agent"  no completeness guard (Radiant, FragPipe, Sage, AlphaDIA): the job
                            never stages; step 7c's check + stage does, after verifying
      fran "off"            --no-fran, FRAN_DEPOSIT=off when the job was generated, or qc=True
    fran_name: the session's descriptive name, handed to `stage --name` -- FRAN's corpus name,
      and one of the names its QC rule reads.
    qc: True (--qc: an instrument QC / standard run -- the hook NEVER stages it, whatever its
      name; the agent decides later), False (--not-qc, passed on to stage), None (stage decides
      from the names it has). Decided here, not later: a QC session whose name has no "QC" token
      would otherwise be staged by the hook before the agent's own `stage --name "... QC"`."""
    return {"run_log": "on",
            "slack": "on" if slack and not opted_out() else "off",
            "fran": ("off" if not fran or qc is True or fran_opted_out() else
                     "stage" if fran_guarded else "left_to_agent"),
            "fran_name": fran_name or None,
            "qc": qc}

_HEADER_LINE = re.compile(r"^(#|\s*$)")


def wrap_job_script(script, out, *, final, time_limit_h, stage, slack=True, fran=True,
                    fran_guarded=False, fran_name=None, qc=None, py=None, scripts_dir=None):
    """A SLURM job script that runs the job-end hook (search_done: run log -> FRAN -> Slack)
    however it ends.

    ONE definition, used by every generator (run_search.emit_sbatch, diann_parallel.py,
    radiant_parallel.py). Everything after the shebang/#SBATCH header runs, unchanged, in a
    background subshell that the batch shell `wait`s on:

    * The time limit. SLURM ends a job with SIGTERM to the BATCH SHELL ONLY, and bash runs a
      trap only after the foreground command returns -- so with the work in the foreground
      (DIA-NN, for hours) the trap never ran: measured on HIVE 2026-09-24 (job 23990731, a
      1-minute limit), TERM at 1:15, no trap, SIGKILL 130 s later (KillWait). `wait` IS
      interrupted by a trapped signal, so the TERM handler runs at once, reports, and re-raises
      TERM so the job still dies by the signal. (Without a TERM trap, bash's EXIT trap sees
      $? = 0 on SIGTERM -- measured, bash 5.1.16 -- and a timed-out search would read as done.)
    * Every other ending: the EXIT trap gets the subshell's status from `wait` and exits with
      it, so the job's exit code, `set -e` and the report guards are exactly what they were.
      The body keeps its own traps (a subshell does not inherit the batch shell's).
    * Only the notifier's PATH is written here; the notifier reads the webhook itself, so
      nothing secret is ever in a job script.
    * slack / fran / fran_guarded: job_end_plan() -- the opt-outs (--no-notify, --no-fran,
      SKILL_SLACK=0 / FRAN_DEPOSIT=off at generation) are baked in as --no-slack / --no-fran,
      and only a route with a completeness guard gets --fran-stage.
    * The header is every LEADING comment or blank line (not only #SBATCH), so a comment between
      two #SBATCH lines cannot push the second below the hook, where SLURM would ignore it.
    * scripts_dir: where notify_slack.py and its helpers live (default: beside this file).
    """
    lines = script.rstrip("\n").split("\n")
    n = 0
    while n < len(lines) and _HEADER_LINE.match(lines[n]):
        n += 1
    head, body = lines[:n], lines[n:]
    q = shlex.quote
    here = scripts_dir or HERE or os.path.dirname(os.path.abspath(sys.argv[0]))
    notifier = os.path.join(here, "notify_slack.py")
    py = py or sys.executable or "python3"
    mins = int(round(float(time_limit_h) * 60)) if time_limit_h else 0
    plan = job_end_plan(slack=slack, fran=fran, fran_guarded=fran_guarded,
                        fran_name=fran_name, qc=qc)
    flags = (("" if final else " --fail-only")
             + ("" if plan["slack"] == "on" else " --no-slack")
             + {"stage": " --fran-stage", "left_to_agent": "", "off": " --no-fran"}[plan["fran"]]
             + {True: " --qc", False: " --not-qc", None: ""}[plan["qc"]]
             + (f" --fran-name {shlex.quote(plan['fran_name'])}" if plan["fran_name"] else ""))
    hook = [
        "# Job end (notify_slack.py search-done): log the run (record_run.py), hand a finished Core",
        "# search to FRAN (fran_deposit.py stage), post to the Core's Slack channel. The work runs",
        "# in the ( ... ) below, which this shell waits on, so a time-limit SIGTERM is handled at",
        "# once. Nothing secret is in this file, and the job's exit status is unchanged.",
        "# references/notifications.md",
        "_JOB_T0=$(date +%s)",
        "_job_end_hook() {",
        f"  local py={q(py)}; [ -x \"$py\" ] || py=python3",
        f"  [ -f {q(notifier)} ] || return 0",
        '  local to=""; command -v timeout >/dev/null 2>&1 && to="timeout 480"',
        f"  $to \"$py\" {q(notifier)} search-done --from-job --out {q(out)} "
        f"--exit-code \"$1\" --signal \"$2\" --started \"$_JOB_T0\" "
        f"--time-limit-min {mins} --stage {q(stage)}{flags} </dev/null || true",
        "}",
        "_job_end_exit() {",
        "  local rc=$1",
        "  trap - EXIT TERM",
        # 143: the work itself died of SIGTERM -- a signal, not a plain failure
        '  if [ "$rc" -eq 143 ]; then _job_end_hook "$rc" TERM; else _job_end_hook "$rc" ""; fi',
        '  exit "$rc"',
        "}",
        "_job_end_term() {",
        "  trap - EXIT TERM",
        "  _job_end_hook 143 TERM",
        "  kill -TERM $$",
        "}",
        "trap '_job_end_exit $?' EXIT",
        "trap _job_end_term TERM",
        "(",
    ]
    tail = [")&", "_job_work=$!", 'wait "$_job_work"']
    return "\n".join(head + hook + body + tail) + "\n"


# ── CLI ──────────────────────────────────────────────────────────


def main(argv=None):
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--dry-run", action="store_true",
                        help="print the JSON payload and send nothing")
    ap = argparse.ArgumentParser(description=__doc__, parents=[common],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--test", action="store_true",
                    help="send one short test message; prints only sent / not sent")
    sub = ap.add_subparsers(dest="cmd")
    s = sub.add_parser("search-done", parents=[common],
                       help="the job-end hook: a search job ended (run by the job's own trap): "
                            "run log -> FRAN stage (on success) -> Slack")
    s.add_argument("--out", required=True, help="the search output directory")
    s.add_argument("--status", choices=["ok", "failed"],
                   help="default: from --exit-code (0 or absent = ok)")
    s.add_argument("--exit-code", type=int)
    s.add_argument("--signal", default="", help="TERM when the job was signalled")
    s.add_argument("--started", help="epoch seconds the job started")
    s.add_argument("--time-limit-min", type=int, default=0)
    s.add_argument("--stage", help="which job of the route this is")
    s.add_argument("--fail-only", action="store_true",
                   help="not the route's last job: post only a failure")
    s.add_argument("--from-job", action="store_true",
                   help="refuse outside a SLURM job (what the generated trap passes)")
    s.add_argument("--no-slack", action="store_true",
                   help="log the run and stage for FRAN, but post nothing to Slack")
    fr = s.add_mutually_exclusive_group()
    fr.add_argument("--fran-stage", action="store_true",
                    help="this route has a completeness guard: stage a finished search for FRAN")
    fr.add_argument("--no-fran", action="store_true",
                    help="never stage for FRAN (--no-fran / FRAN_DEPOSIT=off at generation)")
    s.add_argument("--fran-name", help="the session's descriptive name, for `stage --name`")
    qcg = s.add_mutually_exclusive_group()
    qcg.add_argument("--qc", action="store_true",
                     help="an instrument QC / standard run: never staged from the job")
    qcg.add_argument("--not-qc", action="store_true", help="not a QC run: passed on to stage")
    an = sub.add_parser("analysis-done", parents=[common],
                        help="an analysis was finalized: the Slack post only (session.py "
                             "finalize runs record_run.py itself, first)")
    an.add_argument("--session", required=True)
    an.add_argument("--zip")
    r = sub.add_parser("relay", help=argparse.SUPPRESS)   # the HIVE side of a relay
    r.add_argument("--facts-b64", required=True)
    a = ap.parse_args(argv)

    if a.test:
        ok, msg = send_test(dry_run=a.dry_run)
        print("sent" if ok else f"not sent: {msg}")
        return 0 if ok or a.dry_run else 1
    if a.cmd == "search-done":
        ok, msg = search_done(a.out, a.exit_code, a.status, a.signal or None, a.started,
                              a.time_limit_min or None, a.stage, final=not a.fail_only,
                              from_job=a.from_job, dry_run=a.dry_run, slack=not a.no_slack,
                              fran=("off" if a.no_fran else
                                    "stage" if a.fran_stage else "left_to_agent"),
                              fran_name=a.fran_name,
                              qc=(True if a.qc else False if a.not_qc else None))
    elif a.cmd == "analysis-done":
        ok, msg = analysis_done(a.session, a.zip, dry_run=a.dry_run)
    elif a.cmd == "relay":
        try:
            facts = json.loads(base64.b64decode(a.facts_b64).decode("utf-8"))
            os.environ["SKILL_SLACK_NO_RELAY"] = "1"
            ok, msg = deliver(facts, relay_ok=False)
        except Exception as e:  # noqa: BLE001
            ok, msg = False, _scrub(f"{type(e).__name__}: {e}")
        print("sent" if ok else (msg if msg.startswith("not sent") else f"not sent: {msg}"))
        return 0
    else:
        ap.print_help(sys.stderr)
        return 2
    say(msg)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception as e:  # noqa: BLE001 -- the last line of rule 1
        say(_scrub(f"not sent: {type(e).__name__}: {e}"))
        sys.exit(0)
