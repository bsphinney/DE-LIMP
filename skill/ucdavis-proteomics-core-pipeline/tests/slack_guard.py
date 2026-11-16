"""Keep every test away from the Core's real Slack channel (review A10, Ruling 109.1).

notify_slack.py reads its group webhook file (GROUP_FILE) and the Core's group folder
(CORE_GROUP_DIR) from the environment ONCE, at import. A test that sets SKILL_SLACK_GROUP_FILE
later -- mock.patch.dict(os.environ, ...) around an in-process call -- therefore did not move them,
and on HIVE an in-process "no webhook" call found the Core's real
/quobyte/proteomics-grp/.config/skill_slack_webhook and posted a test alert to the Core channel
on every suite run (test_slack_notify.AlertApi.test_sent_off_and_failed).

install() closes that for the whole test process, whichever module imports it first:
  * notify_slack.GROUP_FILE and CORE_GROUP_DIR point into a temp folder that holds nothing. A test
    that wants a group webhook patches them itself (mock.patch.object), with a file in its own
    temp dir, as test_slack_notify.Resolution does;
  * the guards, belt and braces, each RECORDING what it stopped:
      - resolve_webhook(): a webhook read from a FILE outside a temp dir is refused (None);
      - urllib.request.urlopen: a slack.com URL is refused unless it carries the tests' fake
        token mark (FAKE_MARK; the DNS tests post FAKE_SLACK to a patched getaddrinfo);
      - the relay through HIVE (notify_slack._relay): refused unless its hive_exec is a test's
        fake in a temp dir.
violations() lists what was stopped; check() raises with it. test_slack_notify's tearDownModule
calls check(), and test_slack_guard.py proves the guard with a webhook file at a non-temp path.

Imported by job_env.py (the shared isolation, imported by every module that runs a job or a
finalize), so a whole `unittest discover` run is covered from its discovery on. Stdlib only.
"""
import os
import re
import shutil
import sys
import tempfile
import urllib.request
from unittest import mock

#: what every fake webhook in the tests carries (test_slack_notify.FAKE_SLACK); a real one never
FAKE_MARK = "T0FAKE000"
_SLACK_URL = re.compile(r"https?://([a-z0-9.-]*\.)?slack\.com(/|:|$)", re.IGNORECASE)

_STATE = {"installed": False, "tmp": None, "violations": [], "patchers": []}


def _temp_roots():
    roots = {tempfile.gettempdir(), "/tmp", "/private/tmp", "/var/folders", "/private/var/folders"}
    if _STATE["tmp"]:
        roots.add(_STATE["tmp"])
    return {os.path.realpath(r) for r in roots}


def inside_temp(path):
    """True when `path` (resolved) lies in a temp dir -- where a test's own files live."""
    p = os.path.realpath(os.path.expanduser(str(path)))
    return any(p == r or p.startswith(r + os.sep) for r in _temp_roots())


def _refuse(what):
    _STATE["violations"].append(what)
    return what


def install():
    """Idempotent; -> the temp folder GROUP_FILE / CORE_GROUP_DIR now point into."""
    if _STATE["installed"]:
        return _STATE["tmp"]
    scripts = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts")
    if scripts not in sys.path:
        sys.path.insert(0, scripts)
    import notify_slack as ns
    tmp = tempfile.mkdtemp(prefix="slack_guard_")
    _STATE["tmp"] = tmp

    def start(p):
        p.start()
        _STATE["patchers"].append(p)

    start(mock.patch.object(ns, "GROUP_FILE", os.path.join(tmp, "no_group_webhook")))
    start(mock.patch.object(ns, "CORE_GROUP_DIR", os.path.join(tmp, "no_core_group")))

    real_resolve = ns.resolve_webhook

    def resolve_webhook():
        url, source = real_resolve()
        if url and not str(source).startswith("$"):
            # a file: the user's (~ = the test's HOME) or the group's
            path = os.path.expanduser(ns.USER_FILE) if str(source).startswith("~") else source
            if not inside_temp(path):
                return None, _refuse(f"slack_guard: refused a webhook read from {source}, "
                                     f"outside the test's temp dir")
        return url, source
    start(mock.patch.object(ns, "resolve_webhook", resolve_webhook))

    real_relay = ns._relay

    def relay(facts, hx):
        if not inside_temp(hx):
            return False, _refuse(f"slack_guard: refused a relay through HIVE ({hx})")
        return real_relay(facts, hx)
    start(mock.patch.object(ns, "_relay", relay))

    real_urlopen = urllib.request.urlopen

    def urlopen(req, *a, **k):
        url = str(getattr(req, "full_url", req))
        if _SLACK_URL.match(url) and FAKE_MARK not in url:
            raise OSError(_refuse("slack_guard: refused a request to Slack that is not a "
                                  "test's fake URL"))
        return real_urlopen(req, *a, **k)
    start(mock.patch.object(urllib.request, "urlopen", urlopen))
    _STATE["installed"] = True
    return tmp


def violations():
    return list(_STATE["violations"])


def clear():
    del _STATE["violations"][:]


def check():
    """Raise when anything was stopped since the last clear(): a test reached for real Slack."""
    v = violations()
    if v:
        raise AssertionError("a test reached for the Core's real Slack (slack_guard stopped "
                             "it):\n  " + "\n  ".join(v))


def uninstall():
    """Tests of the guard itself only."""
    for p in reversed(_STATE["patchers"]):
        p.stop()
    del _STATE["patchers"][:]
    if _STATE["tmp"]:
        shutil.rmtree(_STATE["tmp"], ignore_errors=True)
    _STATE.update(installed=False, tmp=None)
