#!/usr/bin/env python3
"""html_to_pdf.py: the report PDF, printed by a headless browser -- never fatal.

Brett, 2026-09-25: "a PDF, something that has the figures built in that google llm can read
and make a podcast about". Guards (stdlib; the real-browser test skips without one):
  * a browser is found in the documented order ($PDF_BROWSER, PATH, install locations);
  * the found-browser path: a throw-away profile, headless print, the PDF moved into place --
    and a browser that writes the PDF but does not exit is stopped (headless Chrome with a
    fresh profile does exactly that on macOS);
  * no browser: no PDF, and a note saying how to make one by hand;
  * a browser that never finishes is stopped at the timeout;
  * make_analysis_html --no-pdf, and finalize's [INFO] line when there is no HTML to print.
"""
import os
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import html_to_pdf as h2p  # noqa: E402
import make_deposit        # noqa: E402

PDF = b"%PDF-1.4\n1 0 obj << /Type /Page >> endobj\n2 0 obj << /Type /Pages >> endobj\ntrailer\n%%EOF\n"


class FakeBrowser:
    """Stands in for subprocess.Popen. writes: put a PDF at --print-to-pdf; exits: whether
    poll() ever reports an exit (headless Chrome with a fresh profile does not)."""
    calls = []

    def __init__(self, writes=True, exits=False):
        self.writes, self.exits = writes, exits

    def __call__(self, cmd, **kw):
        FakeBrowser.calls.append(cmd)
        self.cmd, self.stopped, self.returncode = cmd, False, None
        out = next(c.split("=", 1)[1] for c in cmd if c.startswith("--print-to-pdf="))
        self.profile = next(c.split("=", 1)[1] for c in cmd if c.startswith("--user-data-dir="))
        if self.writes:
            with open(out, "wb") as fh:
                fh.write(PDF)
        return self

    def poll(self):
        if self.exits or self.stopped:
            self.returncode = 0 if not self.stopped else -15
            return self.returncode
        return None

    def terminate(self):
        self.stopped = True

    def kill(self):
        self.stopped = True

    def wait(self, timeout=None):
        return self.returncode


class Convert(unittest.TestCase):
    def setUp(self):
        self._td = tempfile.TemporaryDirectory()
        self.html = os.path.join(self._td.name, "Analysis_Report.html")
        with open(self.html, "w") as fh:
            fh.write("<!doctype html><title>t</title><p>x</p>")
        self.pdf = os.path.join(self._td.name, "Analysis_Report.pdf")

    def tearDown(self):
        self._td.cleanup()

    def test_browser_found_prints_and_is_stopped_once_the_pdf_is_complete(self):
        fake = FakeBrowser(writes=True, exits=False)
        ok, note = h2p.convert(self.html, self.pdf, timeout=10, browser="/x/chrome",
                               popen=fake, poll=0.01)
        self.assertTrue(ok, note)
        self.assertTrue(fake.stopped)                       # did not wait the 10 s out
        with open(self.pdf, "rb") as fh:
            self.assertTrue(fh.read().startswith(b"%PDF"))
        self.assertIn("1 pages", note)
        cmd = fake.cmd
        self.assertEqual(cmd[0], "/x/chrome")
        self.assertIn("--headless=new", cmd)
        self.assertIn("--no-pdf-header-footer", cmd)
        self.assertTrue(cmd[-1].startswith("file://"))
        # a throw-away profile, removed afterwards -- the user's own is never touched
        self.assertTrue(os.path.basename(fake.profile).startswith("html2pdf-profile-"))
        self.assertFalse(os.path.exists(fake.profile))
        self.assertFalse(os.path.exists(self.pdf + ".part"))

    def test_no_browser_writes_no_pdf_and_says_how(self):
        with mock.patch.object(h2p, "find_browser", return_value=None):
            ok, note = h2p.convert(self.html, self.pdf)
        self.assertFalse(ok)
        self.assertIn("no Chrome, Chromium or Edge found", note)
        self.assertIn("Save as PDF", note)
        self.assertFalse(os.path.exists(self.pdf))

    def test_no_browser_via_the_lookup(self):
        self.assertIsNone(h2p.find_browser({}, which=lambda n: None, isfile=lambda p: False))

    def test_timeout_stops_the_browser(self):
        fake = FakeBrowser(writes=False, exits=False)
        ok, note = h2p.convert(self.html, self.pdf, timeout=0.3, browser="/x/chrome",
                               popen=fake, poll=0.02)
        self.assertFalse(ok)
        self.assertIn("did not finish within 0.3s", note)
        self.assertTrue(fake.stopped)
        self.assertFalse(os.path.exists(self.pdf))

    def test_exit_without_a_pdf_is_reported(self):
        ok, note = h2p.convert(self.html, self.pdf, timeout=5, browser="/x/chrome",
                               popen=FakeBrowser(writes=False, exits=True), poll=0.01)
        self.assertFalse(ok)
        self.assertIn("without a PDF", note)

    def test_lookup_order(self):
        seen = {"/opt/my/chrome", "/usr/bin/chromium",
                "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"}
        isfile = seen.__contains__
        self.assertEqual(h2p.find_browser({"PDF_BROWSER": "/opt/my/chrome"},
                                          which=lambda n: None, isfile=isfile), "/opt/my/chrome")
        self.assertEqual(h2p.find_browser({}, which=lambda n: "/usr/bin/chromium"
                                          if n == "chromium" else None, isfile=isfile),
                         "/usr/bin/chromium")
        self.assertEqual(h2p.find_browser({}, which=lambda n: None, isfile=isfile),
                         "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome")
        win = h2p.candidates({"PROGRAMFILES": r"C:\Program Files"})
        self.assertIn(os.path.join(r"C:\Program Files", r"Google\Chrome\Application\chrome.exe"), win)

    @unittest.skipUnless(h2p.find_browser(), "no Chrome/Chromium/Edge on this machine")
    def test_real_browser(self):
        with open(self.html, "w") as fh:
            fh.write("<!doctype html><meta charset=utf-8><title>t</title><h1>Report</h1>"
                     "<img alt=x src='data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAf"
                     "FcSJAAAADUlEQVR42mP8z8BQDwAEhQGAhKmMIQAAAABJRU5ErkJggg=='>")
        ok, note = h2p.convert(self.html, self.pdf, timeout=120)
        self.assertTrue(ok, note)
        with open(self.pdf, "rb") as fh:
            self.assertGreaterEqual(h2p.page_count(fh.read()), 1)


class FinalizeKeepsThePdf(unittest.TestCase):
    """finalize's tidy step moved every loose .pdf in output/ into figures/ -- the report's PDF
    too -- so each finalize reprinted it and left a stray figures/Analysis_Report.pdf."""

    def test_an_up_to_date_pdf_stays_and_is_not_reprinted(self):
        import test_deposit_package as tdp
        with tempfile.TemporaryDirectory() as d:
            p = tdp.dia_session(d)
            html = os.path.join(p["output_dir"], "Analysis_Report.html")
            pdf = os.path.join(p["output_dir"], "Analysis_Report.pdf")
            with open(html, "w") as fh:
                fh.write("<html><body>r</body></html>")
            with open(pdf, "wb") as fh:
                fh.write(b"%PDF-1.4\nCANARY\n%%EOF\n")
            os.utime(html, (1000, 1000))                     # the PDF is newer than the HTML
            for _ in range(2):
                r = tdp.finalize(p["session_dir"])
                self.assertEqual(r.returncode, 0, r.stderr)
                with open(pdf, "rb") as fh:
                    self.assertTrue(b"CANARY" in fh.read(), "the PDF was reprinted")
                self.assertFalse(os.path.exists(os.path.join(p["figures_dir"],
                                                             "Analysis_Report.pdf")))
                with open(p["manifest_txt"]) as fh:
                    self.assertIn("made from the current HTML", fh.read())


class Integration(unittest.TestCase):
    def test_make_analysis_html_no_pdf(self):
        with tempfile.TemporaryDirectory() as tmp:
            rep = os.path.join(tmp, "AI_Analysis_Report.md")
            with open(rep, "w") as fh:
                fh.write("# T\n\n## Overview\n\nText.\n")
            out = os.path.join(tmp, "Analysis_Report.html")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "make_analysis_html.py"),
                                "--report", rep, "--out", out, "--no-pdf"],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn('"pdf": null', r.stdout)
            self.assertIn("INFO: no PDF -- --no-pdf was given", r.stderr)
            self.assertFalse(os.path.exists(os.path.join(tmp, "Analysis_Report.pdf")))

    def test_manifest_info_is_not_a_skip(self):
        m = make_deposit.Manifest()
        m.info("Report PDF (output/Analysis_Report.pdf)", "no browser; print it by hand")
        self.assertTrue(m.lines[0].startswith("[INFO]    Report PDF"))
        self.assertEqual(m.n_skipped, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
