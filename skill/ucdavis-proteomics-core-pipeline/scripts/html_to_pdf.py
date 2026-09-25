#!/usr/bin/env python3
"""
html_to_pdf.py -- print a self-contained HTML page to PDF with a headless Chromium browser.

Why: Brett (2026-09-25) wanted "a PDF, something that has the figures built in that google
llm can read and make a podcast about". NotebookLM reads a PDF's text AND its images, and a
PDF also prints. The HTML report is already self-contained with a print stylesheet
(report_style.py), so a headless Chrome/Chromium/Edge print of it is the PDF -- no LaTeX, no
Python PDF library to install.

Never fatal. HIVE compute nodes usually have no browser, and that is fine: in hive_remote,
session.py finalize runs on the laptop, which has one. When no browser is found, or it fails
or runs out of time, convert() returns (False, reason) and the caller records an [INFO] line
saying how to make the PDF by hand. The browser runs with a throw-away --user-data-dir, so it
never touches the user's own profile, and is killed after --timeout seconds (default 120).

  python3 html_to_pdf.py --html output/Analysis_Report.html [--pdf output/Analysis_Report.pdf]

Browser search order: $PDF_BROWSER; google-chrome / chromium / msedge ... on PATH; the macOS
and Windows install locations.
"""
import argparse
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import time

HOW_BY_HAND = ("to make one by hand, open the .html in any browser, choose Print, and "
               "Save as PDF")
PATH_NAMES = ("google-chrome", "google-chrome-stable", "chromium", "chromium-browser",
              "chrome", "microsoft-edge", "microsoft-edge-stable", "msedge")
MAC_APPS = ("Google Chrome.app/Contents/MacOS/Google Chrome",
            "Chromium.app/Contents/MacOS/Chromium",
            "Microsoft Edge.app/Contents/MacOS/Microsoft Edge",
            "Google Chrome Canary.app/Contents/MacOS/Google Chrome Canary")
WIN_EXES = (r"Google\Chrome\Application\chrome.exe", r"Microsoft\Edge\Application\msedge.exe",
            r"Chromium\Application\chrome.exe")


def candidates(env=None):
    """Install locations to try after PATH, most likely first."""
    env = os.environ if env is None else env
    out = []
    for root in ("/Applications", os.path.expanduser("~/Applications")):
        out += [os.path.join(root, a) for a in MAC_APPS]
    for var in ("PROGRAMFILES", "PROGRAMFILES(X86)", "LOCALAPPDATA"):
        if env.get(var):
            out += [os.path.join(env[var], e) for e in WIN_EXES]
    # Git Bash / MSYS spelling of the same Windows locations
    out += [f"/c/Program Files/{e.replace(chr(92), '/')}" for e in WIN_EXES]
    out += [f"/c/Program Files (x86)/{e.replace(chr(92), '/')}" for e in WIN_EXES]
    return out


def find_browser(env=None, which=shutil.which, isfile=os.path.isfile):
    env = os.environ if env is None else env
    if env.get("PDF_BROWSER") and isfile(env["PDF_BROWSER"]):
        return env["PDF_BROWSER"]
    for name in PATH_NAMES:
        p = which(name)
        if p:
            return p
    return next((p for p in candidates(env) if isfile(p)), None)


def page_count(data):
    return len(re.findall(rb"/Type\s*/Page(?![a-zA-Z])", data))


def _complete(path):
    """A finished PDF: starts %PDF, ends %%EOF."""
    try:
        with open(path, "rb") as fh:
            head = fh.read(5)
            fh.seek(max(0, os.path.getsize(path) - 1024))
            tail = fh.read()
        return head.startswith(b"%PDF") and b"%%EOF" in tail
    except OSError:
        return False


def _stop(proc):
    """Stop the browser and the helper processes it started (its own process group, which
    start_new_session gave it -- so nothing else is signalled)."""
    if proc.poll() is not None:
        return
    group = isinstance(proc, subprocess.Popen) and hasattr(os, "killpg")
    try:
        os.killpg(proc.pid, 15) if group else proc.terminate()
        proc.wait(timeout=5)
    except Exception:
        try:
            os.killpg(proc.pid, 9) if group else proc.kill()
        except Exception:
            pass


def convert(html_path, pdf_path, timeout=120, browser=None, popen=subprocess.Popen, env=None,
            poll=0.25):
    """-> (ok, note). Writes pdf_path only when the browser produced a complete PDF.

    The browser is watched rather than waited for: with a fresh --user-data-dir (so the
    user's own profile is never touched) headless Chrome on macOS writes the PDF within
    seconds and then does not exit (measured 2026-09-25: PDF done at ~3 s, process still up
    at 60 s). So a complete PDF whose size has stopped changing ends the run, and the
    browser is stopped; a browser that exits, or runs past `timeout`, ends it too."""
    b = browser or find_browser(env)
    if not b:
        return False, (f"no Chrome, Chromium or Edge found (PATH, the macOS and Windows install "
                       f"locations, $PDF_BROWSER); {HOW_BY_HAND}")
    html_abs = os.path.abspath(html_path)
    pdf_abs = os.path.abspath(pdf_path)
    profile = tempfile.mkdtemp(prefix="html2pdf-profile-")
    tmp_pdf = pdf_abs + ".part"
    cmd = [b, "--headless=new", "--disable-gpu", "--no-first-run", "--no-default-browser-check",
           "--disable-extensions", "--use-mock-keychain", "--password-store=basic",
           f"--user-data-dir={profile}",
           "--no-pdf-header-footer", "--print-to-pdf-no-header",   # new and old spellings
           f"--print-to-pdf={tmp_pdf}", pathlib.Path(html_abs).as_uri()]
    errlog = tempfile.TemporaryFile()
    proc, timed_out = None, False
    try:
        if os.path.exists(tmp_pdf):
            os.remove(tmp_pdf)
        proc = popen(cmd, stdout=subprocess.DEVNULL, stderr=errlog,
                     **({"start_new_session": True} if hasattr(os, "killpg") else {}))
        deadline, last = time.monotonic() + timeout, -1
        while True:
            done = proc.poll() is not None
            size = os.path.getsize(tmp_pdf) if os.path.exists(tmp_pdf) else -1
            if size > 0 and size == last and _complete(tmp_pdf):
                break                                   # finished writing; exited or not
            if done and size == last:
                break                                   # exited; whatever it wrote is final
            if time.monotonic() > deadline:
                timed_out = True
                break
            last = size
            time.sleep(poll)
    except OSError as e:
        return False, f"{b} could not be started ({e}); {HOW_BY_HAND}"
    finally:
        if proc is not None:
            _stop(proc)
        shutil.rmtree(profile, ignore_errors=True)
    ok = os.path.isfile(tmp_pdf) and _complete(tmp_pdf)
    if not ok:
        if os.path.exists(tmp_pdf):
            os.remove(tmp_pdf)
        if timed_out:
            return False, (f"{os.path.basename(b)} did not finish within {timeout}s and was "
                           f"stopped; {HOW_BY_HAND}")
        errlog.seek(0)
        err = errlog.read().decode("utf-8", "replace")
        tail = " ".join(ln for ln in err.splitlines() if "rror" in ln)[-200:]
        return False, (f"{os.path.basename(b)} exited {proc.returncode if proc else '?'} without "
                       f"a PDF{': ' + tail if tail else ''}; {HOW_BY_HAND}")
    with open(tmp_pdf, "rb") as fh:
        data = fh.read()
    os.replace(tmp_pdf, pdf_abs)
    return True, (f"{page_count(data)} pages, {len(data) / 1e6:.1f} MB, printed by "
                  f"{os.path.basename(b)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--html", required=True)
    ap.add_argument("--pdf", help="default: the .html path with .pdf")
    ap.add_argument("--timeout", type=int, default=120)
    a = ap.parse_args()
    pdf = a.pdf or os.path.splitext(a.html)[0] + ".pdf"
    ok, note = convert(a.html, pdf, timeout=a.timeout)
    print(json.dumps({"pdf": pdf if ok else None, "ok": ok, "note": note}, indent=2))
    sys.exit(0)            # never fatal: the note says what happened


if __name__ == "__main__":
    main()
