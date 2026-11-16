#!/usr/bin/env python3
"""Nothing delivered for a default timsTOF DDA search names a licence, a jar or FRAGPIPE_TOOLS_FOLDER
(2.15 Ffix, review M1, Ruling 107.1).

From 2.15 FragPipe is the default for timsTOF DDA-PASEF, so every such delivery is a FragPipe one.
The review found detect_env.sh's FragPipe note ("... licence-gated -- set FRAGPIPE_TOOLS_FOLDER")
copied twice into workflow.manifest.json (its `notes` and its `platform`), which the reproducibility
bundle carries to the client; and the search record (search_provenance.json, delivered) named the
tools folder's source by the variable's name. Ruling 56 keeps the licence notes where they are --
detect_env.sh's wording and tools.json, staff files -- so the note now goes to the staff-only
record beside the manifest, and the delivered record says where the folder came from in its own
words.

The test runs step 4 (resolve_defaults.py, the default route), the search generation
(run_search.py --sbatch) and the bundle (provenance.py) as the skill runs them, with detect_env.sh's
real FragPipe note, then scans EVERY file of the bundle and every delivered search file
(core_submission.SEARCH_FILES) -- the folder, not one file. Stdlib only; no network.
"""
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import finish_fixture  # noqa: E402
import staff  # noqa: E402
from fake_fragpipe_install import fake_install, tools_json  # noqa: E402
from job_env import job_env  # noqa: E402
from test_fragpipe_dda_stage import FAKE_FRAGPIPE, _write  # noqa: E402

PY = sys.executable
# what no delivered file may carry (Ruling 107.1)
FORBIDDEN = re.compile(r"licen[cs]|\bjars?\b|FRAGPIPE_TOOLS_FOLDER", re.IGNORECASE)
# detect_env.sh's FragPipe note on Linux with Java, read from the script itself: its wording is
# unchanged by Ruling 107.1 (pinned below), only where it goes changed
NOTE_LINE = re.compile(r'^\s*FRAGPIPE_NOTE="(Linux build;[^"]*)"\s*$', re.MULTILINE)
DETECT_ENV_NOTE = ("Linux build; MSFragger/IonQuant/diaTracer are licence-gated -- set "
                   "FRAGPIPE_TOOLS_FOLDER")


def detect_env_note():
    with open(os.path.join(SCRIPTS, "detect_env.sh")) as fh:
        return NOTE_LINE.search(fh.read()).group(1)


def delivered_search_files():
    import core_submission
    return core_submission.SEARCH_FILES


class DefaultTimsTofDdaDelivery(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        self.env = job_env(self.d, PATH=os.pathsep.join([os.path.dirname(PY), "/usr/bin",
                                                         "/bin"]))
        # no Core install and no tools folder from the environment: tools.json is the one carrier
        self.env["FRAGPIPE_CORE_INSTALL"] = os.path.join(self.d, "no-core-install")
        self.env.pop("FRAGPIPE_TOOLS_FOLDER", None)
        self.fp_bin = _write(os.path.join(self.d, "bin", "fragpipe"), FAKE_FRAGPIPE, 0o755)
        fp = fake_install(self.d)
        # as acquire_tools.sh records a folder found through the variable (engine_choice's own
        # `why` wording), with the jars
        fp["fragpipe_tools"]["why"] = (f"MSFragger 4.4.1, IonQuant 1.11.20 in "
                                       f"{fp['tools_folder']} (FRAGPIPE_TOOLS_FOLDER)")
        self.tools = tools_json(self.d, self.fp_bin, fp)
        env_json = {"os": "linux", "arch": "x86_64", "has_java":
                    "openjdk version 21.0.4 2024-07-16",
                    "engines": {"diann": {"available": True, "how": "native",
                                          "note": "native Linux build"},
                                "sage": {"available": True, "how": "native",
                                         "note": "native build for this OS/arch (DDA only)"},
                                "fragpipe": {"available": True, "how": "native",
                                             "note": detect_env_note()},
                                "radiant": {"available": False, "how": None,
                                            "note": "Seer ships no native binary -- needs "
                                                    "Docker or Apptainer"}}}
        self.env_json = _write(os.path.join(self.d, "env.json"), json.dumps(env_json))
        self.raw = os.path.join(self.d, "raw")
        self.files = []
        for r in ("S1_dda_A_1", "S2_dda_B_1"):
            dd = os.path.join(self.raw, r + ".d")
            _write(os.path.join(dd, "analysis.tdf"), "")
            self.files.append(dd)

    def run_ok(self, argv, what):
        p = subprocess.run([PY, *argv], capture_output=True, text=True, env=self.env,
                           cwd=self.d, timeout=300)
        self.assertEqual(p.returncode, 0, f"{what}:\n{p.stdout[-2000:]}\n{p.stderr[-2000:]}")
        return p

    def deliver(self):
        wf = os.path.join(self.d, "wf")
        self.run_ok([os.path.join(SCRIPTS, "resolve_defaults.py"), "--acquisition", "DDA",
                     "--instrument", "timsTOF Pro 2", "--env", self.env_json, "--tools",
                     self.tools, "--dest", wf], "step 4")
        man_path = os.path.join(wf, "workflow.manifest.json")
        with open(man_path) as fh:
            man = json.load(fh)
        self.assertEqual((man["engine"]["name"], man["engine"]["selection"]["how"]),
                         ("fragpipe", "default"), "not the default FragPipe route")
        search = os.path.join(self.d, "search")
        fasta = _write(os.path.join(self.d, "db.fasta"),
                       ">sp|P1|ONE_HUMAN one\nPEPTIDER\n>rev_sp|P1|ONE_HUMAN one\nREDITPEP\n")
        self.run_ok([os.path.join(SCRIPTS, "run_search.py"), "--tools", self.tools,
                     "--bundle", man_path, "--params", man["search"]["params_file"],
                     "--fasta", fasta, "--out", search, "--files", *self.files,
                     "--engine", "fragpipe", "--threads", "4", "--no-notify", "--no-fran",
                     "--allow-damaged-tdf", "--sbatch", os.path.join(self.d, "job.sh")],
                    "the search generation")
        home = finish_fixture.fake_pipeline_home(self.d)
        bundle = os.path.join(self.d, "bundle")
        self.run_ok([os.path.join(SCRIPTS, "provenance.py"), "--outdir", bundle,
                     "--setup-json", os.path.join(home, "setup.json"),
                     "--workflow-manifest", man_path, "--params", man["search"]["params_file"],
                     "--tools-json", self.tools, "--fasta", fasta,
                     "--search-record", os.path.join(search, "search_provenance.json"),
                     "--engine", "fragpipe", "--acquisition", "DDA",
                     "--instrument", "timsTOF Pro 2", "--raw", *self.files], "the bundle")
        delivered = [os.path.join(dp, f) for dp, _, fs in os.walk(bundle) for f in fs]
        delivered += [os.path.join(search, f) for f in delivered_search_files()
                      if os.path.isfile(os.path.join(search, f))]
        return wf, search, delivered

    def test_no_delivered_file_names_a_licence_a_jar_or_the_tools_variable(self):
        wf, search, delivered = self.deliver()
        names = {os.path.basename(p) for p in delivered}
        for must in ("workflow.manifest.json", "search_provenance.json", "reproduce.sh",
                     "REPRODUCE.md", "MANIFEST.txt"):
            self.assertIn(must, names, "the scan would miss a delivered file")
        hits = []
        for p in delivered:
            with open(p, errors="replace") as fh:
                for m in FORBIDDEN.finditer(fh.read()):
                    hits.append(f"{os.path.relpath(p, self.d)}: {m.group(0)!r}")
        self.assertEqual(hits, [])

    def test_the_note_is_kept_for_the_staff(self):
        """Ruling 56: the licence note stays as detect_env.sh writes it, in the staff-only record
        beside the manifest (never delivered: staff.is_staff_only)."""
        wf, _, _ = self.deliver()
        self.assertEqual(detect_env_note(), DETECT_ENV_NOTE)
        rec = staff.read(os.path.join(wf, "workflow.manifest.json" + staff.STAFF_SUFFIX))
        self.assertIn(DETECT_ENV_NOTE, [e.get("platform_note") for e in rec])
        with open(os.path.join(wf, "workflow.manifest.json")) as fh:
            man = json.load(fh)
        self.assertEqual(man["platform"], {"available": True, "how": "native",
                                           "note_in": "workflow.manifest.json.staff.json"})
        self.assertEqual(man["notes"], [])


if __name__ == "__main__":
    unittest.main()
