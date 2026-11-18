#!/usr/bin/env python3
"""
Which search engine runs, and why -- engine_choice.py, the one route table (2.15, lane F2).

Rulings 55/56/67/86/90: FragPipe is the default for timsTOF DDA-PASEF (IDs and quant); Sage is
the automatic fallback when FragPipe cannot run on the machine (macOS, no Java, Java < 11) or its
licence-gated MSFragger/IonQuant jars are not installed for this account. Orbitrap and other DDA
stay on Sage, FragPipe on request. An IDs-only search, a non-specific digest and an ion-trap MS2
never change the engine (FragPipe has workflows and a fragment setting for each). DIA -> DIA-NN.

The manifest's selection record is client-safe (no version, path or licence word); the reason
and `why` go only into the staff-only record beside it (review M1, Ruling 78 F4).

One test per rule row, then the resolve_defaults wiring and the pins (the Lumos DDA and the DIA
manifests are byte-identical to 2f9b229c's apart from the added selection record).
Stdlib only; no network; resolve_defaults.py runs as the orchestrator runs it.
"""
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
PINS = os.path.join(HERE, "fixtures", "engine_choice_pins")
sys.path.insert(0, SCRIPTS)

try:
    import engine_choice as ec  # noqa: E402
except ImportError:      # a fail-first run on the base tree: each test that needs it fails alone
    ec = None
# the environment names engine_choice reads (spelled here so the base tree's tests still run)
CORE_ENV, TOOLS_ENV = "FRAGPIPE_CORE_INSTALL", "FRAGPIPE_TOOLS_FOLDER"
import resolve_defaults as rd  # noqa: E402
import estimate_params as ep  # noqa: E402
import staff  # noqa: E402

MAC_NOTE = ("FragPipe ships only a Windows installer and a Linux zip -- there is NO macOS build, "
            "Intel or Apple Silicon. This route (and diaTracer) is unavailable on any Mac; run "
            "it on HIVE.")


def env(os_="linux", fragpipe=True, java="openjdk version 25.0.4.1 2025-10-21", note=None):
    """A detect_env.sh JSON (the fields engine_choice reads)."""
    return {"os": os_, "arch": "x86_64" if os_ == "linux" else "arm64", "has_java": java,
            "engines": {"diann": {"available": True, "note": "native Linux build"},
                        "sage": {"available": True, "note": "native build"},
                        "fragpipe": {"available": fragpipe,
                                     "note": note or ("Linux build" if fragpipe else MAC_NOTE)},
                        "radiant": {"available": False, "note": "needs Docker or Apptainer"}}}


ENV_LINUX = env()
ENV_MAC = env("darwin", False, java=None)
ENV_JAVA8 = env(java="openjdk version 1.8.0_392")   # an env.json from before detect_env's Java floor
FP_LICENSED = {"licensed": True, "ready": True, "folder": "/x/tools", "source": "core_install",
               "java": {"version": "openjdk version 25.0.4.1", "major": 25, "ok": True},
               "why": "MSFragger 4.4.1, IonQuant 1.11.20 in /x/tools (core_install)"}
FP_ABSENT = {"licensed": False, "ready": False, "folder": None, "source": None,
             "java": {"version": "openjdk version 25.0.4.1", "major": 25, "ok": True},
             "why": "MSFragger and IonQuant (licence-gated) were not found for this account "
                    "(/home/u/fragpipe/tools: no MSFragger or IonQuant jar within 2 levels)"}

# licensed, but no runnable FragPipe (review M2's probe: FRAGPIPE_TOOLS_FOLDER set off HIVE and
# acquire_tools.sh's FragPipe download failed) -- engine_choice.fragpipe_tools' own `why` wording
FP_NOT_READY = dict(FP_LICENSED, ready=False,
                    why="MSFragger 4.4.1, IonQuant 1.11.20 in /x/tools (FRAGPIPE_TOOLS_FOLDER); "
                        "but no runnable FragPipe binary")

# what a client sentence must never carry (review M1): a version, a path, a licence word
CLIENT_FORBIDDEN = re.compile(r"\d|/|licen|MSFragger|IonQuant|jar|java", re.IGNORECASE)


class RuleRows(unittest.TestCase):
    """One test per row of the decision table."""

    def test_dia_is_dia_nn_on_every_family(self):
        for fam in ("timstof", "orbitrap", "other"):
            sel, staff_entry = ec.choose("DIA", fam, env=ENV_LINUX)
            self.assertEqual((sel["engine"], sel["how"], sel["sentence"]), ("diann", "default", None))
            self.assertIsNone(staff_entry)

    def test_timstof_dda_is_fragpipe_where_it_runs_and_is_licensed(self):
        sel, st = ec.choose("DDA", "timstof", env=ENV_LINUX, fp_tools=FP_LICENSED)
        self.assertEqual((sel["engine"], sel["default"], sel["how"]), ("fragpipe", "fragpipe", "default"))
        self.assertIsNone(sel["sentence"])
        self.assertEqual(ec.de_method(sel["engine"], "DDA"), "maxlfq")
        self.assertEqual((st["reason"], st["checks"]["tools"]["licensed"]), (None, True))

    def test_timstof_dda_on_a_mac_falls_back_to_sage_unavailable_here(self):
        sel, st = ec.choose("DDA", "timstof", env=ENV_MAC, fp_tools=FP_ABSENT)
        self.assertEqual((sel["engine"], sel["default"], sel["how"]), ("sage", "fragpipe", "fallback"))
        self.assertEqual(st["reason"], "unavailable_here")
        self.assertIn("darwin/arm64", st["why"])
        self.assertIn("NO macOS build", st["why"])
        self.assertEqual(sel["sentence"],
                         "FragPipe, the pipeline's default search engine for timsTOF DDA-PASEF "
                         "data, was not available where this search ran, so the data were "
                         "searched with Sage.")

    def test_timstof_dda_without_the_licensed_jars_falls_back_not_licensed_here(self):
        sel, st = ec.choose("DDA", "timstof", env=ENV_LINUX, fp_tools=FP_ABSENT)
        self.assertEqual((sel["engine"], sel["how"]), ("sage", "fallback"))
        self.assertEqual(st["reason"], "not_licensed_here")
        self.assertIn(FP_ABSENT["why"], st["why"])
        # the same neutral sentence for both reasons: the client never reads which
        self.assertEqual(sel["sentence"], ec.choose("DDA", "timstof", env=ENV_MAC)[0]["sentence"])

    def test_a_licence_nobody_could_check_is_not_licensed(self):
        """Ruling 103: tools.json written without python3 says licensed: null -- Sage fallback with
        'licence not checked on this machine' for the staff, the same client sentence; a requested
        FragPipe is refused on it."""
        unchecked = dict(FP_ABSENT, licensed=None, ready=None,
                         why="not checked: python3 is not on PATH where acquire_tools.sh ran")
        sel, st = ec.choose("DDA", "timstof", env=ENV_LINUX, fp_tools=unchecked)
        self.assertEqual((sel["engine"], sel["how"], st["reason"]), ("sage", "fallback", "licence_not_checked"))
        self.assertTrue(st["why"].startswith("licence not checked on this machine"), st["why"])
        self.assertEqual(sel["sentence"], ec.choose("DDA", "timstof", env=ENV_MAC)[0]["sentence"])
        with self.assertRaises(ec.Refused) as cm:
            ec.choose("DDA", "timstof", requested="fragpipe", env=ENV_LINUX, fp_tools=unchecked)
        self.assertIn("licence not checked on this machine", str(cm.exception))

    def test_licensed_but_not_runnable_falls_back_to_sage(self):
        """Review M2 / Ruling 107.2 (the reviewer's probe): the jars are found but there is no
        runnable FragPipe (tools.json licensed true, ready false). run_search.py refuses such a
        search (refuse_fragpipe_not_ready), so the default must fall back to Sage here --
        how=fallback, the reason staff-only -- never choose FragPipe and refuse at step 7."""
        sel, st = ec.choose("DDA", "timstof", env=ENV_LINUX, fp_tools=FP_NOT_READY)
        self.assertEqual((sel["engine"], sel["default"], sel["how"]), ("sage", "fragpipe", "fallback"))
        self.assertEqual(st["reason"], "not_ready_here")
        self.assertIn(FP_NOT_READY["why"], st["why"])
        self.assertIs(ec.fragpipe_ready({"fragpipe_tools": FP_NOT_READY})["ready"], False)
        # the client reads the same neutral sentence as for every other fallback reason
        self.assertEqual(sel["sentence"], ec.choose("DDA", "timstof", env=ENV_MAC)[0]["sentence"])
        blob = json.dumps(sel)
        for word in ("not_ready_here", "runnable", "/x/tools", "licen"):
            self.assertNotIn(word, blob)
        # a REQUEST for FragPipe is refused with the staff reason, never swapped for Sage
        with self.assertRaises(ec.Refused) as cm:
            ec.choose("DDA", "timstof", requested="fragpipe", env=ENV_LINUX, fp_tools=FP_NOT_READY)
        self.assertIn("no runnable FragPipe binary", str(cm.exception))

    def test_java_older_than_11_is_unavailable_here(self):
        sel, st = ec.choose("DDA", "timstof", env=ENV_JAVA8, fp_tools=FP_LICENSED)
        self.assertEqual((sel["engine"], sel["how"], st["reason"]), ("sage", "fallback", "unavailable_here"))
        self.assertIn("Java 8 found", st["why"])
        # the same floor read from tools.json's java when there is no env.json at all
        old = dict(FP_LICENSED, java={"version": "openjdk version 1.8.0_392", "major": 8, "ok": False})
        sel, st = ec.choose("DDA", "timstof", fp_tools=old)
        self.assertEqual((sel["engine"], st["reason"]), ("sage", "unavailable_here"))

    def test_no_de_timstof_dda_stays_fragpipe_ids_only(self):
        sel, _ = ec.choose("DDA", "timstof", env=ENV_LINUX, fp_tools=FP_LICENSED, ids_only=True)
        self.assertEqual((sel["engine"], sel["how"], sel["ids_only"]), ("fragpipe", "default", True))

    def test_ion_trap_ms2_with_fragpipe_is_not_refused(self):
        sel, st = ec.choose("DDA", "orbitrap", requested="fragpipe", env=ENV_LINUX,
                            fp_tools=FP_LICENSED, cls="orbitrap_iontrap")
        self.assertEqual((sel["engine"], sel["how"]), ("fragpipe", "requested"))
        self.assertEqual(st["instrument_class"], "orbitrap_iontrap")

    def test_non_specific_digest_stays_fragpipe(self):
        sel, st = ec.choose("DDA", "timstof", env=ENV_LINUX, fp_tools=FP_LICENSED, enzyme="none")
        self.assertEqual(sel["engine"], "fragpipe")
        self.assertEqual(st["enzyme"], "none")

    def test_orbitrap_dda_stays_sage_and_records_no_check(self):
        sel, st = ec.choose("DDA", "orbitrap", env=ENV_LINUX, fp_tools=FP_LICENSED)
        self.assertEqual((sel["engine"], sel["default"], sel["how"]), ("sage", "sage", "default"))
        self.assertIsNone(st)

    def test_non_specific_orbitrap_dda_is_fragpipe_by_default(self):
        """Ruling 120 (replaces Ruling 99's row): a non-specific Orbitrap DDA search (--enzyme
        none, or pepsin: no in-silico rule) goes to FragPipe's Nonspecific-peptidome workflow by
        default -- 16,203 peptides against Sage's 179 on the same data (pep22 acceptance)."""
        for enzyme in ("none", "pepsin", "Pepsin"):
            with self.subTest(enzyme):
                sel, st = ec.choose("DDA", "orbitrap", env=ENV_LINUX, fp_tools=FP_LICENSED,
                                    enzyme=enzyme)
                self.assertEqual((sel["engine"], sel["default"], sel["how"], sel["sentence"]),
                                 ("fragpipe", "fragpipe", "default", None))
                self.assertEqual((st["enzyme"], st["reason"]), (enzyme, None))

    def test_non_specific_orbitrap_dda_falls_back_to_sage_where_fragpipe_cannot_run(self):
        for e, ft, reason in ((ENV_MAC, FP_LICENSED, "unavailable_here"),
                              (ENV_LINUX, FP_ABSENT, "not_licensed_here"),
                              (ENV_LINUX, FP_NOT_READY, "not_ready_here")):
            with self.subTest(reason):
                sel, st = ec.choose("DDA", "orbitrap", env=e, fp_tools=ft, enzyme="none")
                self.assertEqual((sel["engine"], sel["default"], sel["how"]),
                                 ("sage", "fragpipe", "fallback"))
                self.assertEqual(st["reason"], reason)
                self.assertIn("non-specific", sel["sentence"])

    def test_sage_asked_for_on_non_specific_orbitrap_dda_is_honoured(self):
        """Sage stays reachable on request (Ruling 120): honoured, recorded as a request."""
        sel, st = ec.choose("DDA", "orbitrap", requested="sage", env=ENV_LINUX,
                            fp_tools=FP_LICENSED, enzyme="none")
        self.assertEqual((sel["engine"], sel["default"], sel["how"]), ("sage", "fragpipe", "requested"))
        self.assertEqual(sel["sentence"], "Sage was chosen for this analysis instead of FragPipe, "
                                          "the pipeline's default search engine for non-specific "
                                          "Orbitrap DDA data.")
        self.assertIsNone(st)

    def test_tryptic_orbitrap_dda_and_other_dda_are_unchanged(self):
        """Ruling 120 moves only the non-specific Orbitrap row: a tryptic (or unstated) digest on
        an Orbitrap stays Sage, and a non-Orbitrap DDA stays Sage whatever the digest."""
        for fam, enzyme in (("orbitrap", None), ("orbitrap", "trypsin"), ("orbitrap", "lysc"),
                            ("other", None), ("other", "none")):
            with self.subTest((fam, enzyme)):
                sel, st = ec.choose("DDA", fam, env=ENV_LINUX, fp_tools=FP_LICENSED, enzyme=enzyme)
                self.assertEqual((sel["engine"], sel["default"], sel["how"], st),
                                 ("sage", "sage", "default", None))
        sel, st = ec.choose("DDA", "orbitrap", requested="fragpipe", env=ENV_LINUX,
                            fp_tools=FP_LICENSED, enzyme="trypsin")
        self.assertEqual((sel["engine"], sel["how"]), ("fragpipe", "requested"))

    def test_timstof_routes_do_not_depend_on_the_digest(self):
        for acq in ("DDA", "DIA"):
            for enzyme in (None, "none", "trypsin"):
                with self.subTest((acq, enzyme)):
                    self.assertEqual(ec.routes_for(acq, "timstof", enzyme),
                                     ec.ROUTES[(acq, "timstof")])

    def test_an_unknown_enzyme_word_is_refused_by_the_route(self):
        with self.assertRaises(ec.Refused) as cm:
            ec.choose("DDA", "orbitrap", env=ENV_LINUX, fp_tools=FP_LICENSED, enzyme="papain")
        self.assertIn("papain", str(cm.exception))

    def test_the_route_rows_are_ruling_99s(self):
        self.assertEqual(ec.ROUTES, {
            ("DIA", "timstof"): ["diann", "fragpipe"],
            ("DIA", "orbitrap"): ["diann", "fragpipe", "radiant"],
            ("DIA", "other"): ["diann"],
            ("DDA", "timstof"): ["fragpipe", "sage"],
            ("DDA", "orbitrap"): ["sage"],
            ("DDA", "other"): ["sage"]})
        # Ruling 120: the one digest row
        self.assertEqual(ec.NONSPECIFIC_ROUTES, {("DDA", "orbitrap"): ["fragpipe", "sage"]})
        self.assertEqual(ec.routes_for("DDA", "orbitrap", "none"), ["fragpipe", "sage"])
        self.assertEqual(ec.routes_for("DDA", "orbitrap", None), ["sage"])

    def test_sage_asked_for_on_timstof_is_a_request_not_a_fallback(self):
        sel, st = ec.choose("DDA", "timstof", requested="sage", env=ENV_LINUX, fp_tools=FP_LICENSED)
        self.assertEqual((sel["engine"], sel["default"], sel["how"]), ("sage", "fragpipe", "requested"))
        self.assertEqual(sel["sentence"], "Sage was chosen for this analysis instead of FragPipe, "
                                          "the pipeline's default search engine for timsTOF "
                                          "DDA-PASEF data.")
        self.assertIsNone(st)

    def test_asking_for_the_default_is_the_default(self):
        """reproduce.sh replays step 4 with --engine <what ran>; the replayed manifest must equal
        the original (test_provenance_replay), so --engine diann on DIA is `default`, not `requested`."""
        self.assertEqual(ec.choose("DIA", "timstof", requested="diann", env=ENV_LINUX)[0],
                         ec.choose("DIA", "timstof", env=ENV_LINUX)[0])
        sel, st = ec.choose("DDA", "timstof", requested="fragpipe", env=ENV_LINUX, fp_tools=FP_LICENSED)
        self.assertEqual((sel["engine"], sel["how"], sel["sentence"]), ("fragpipe", "default", None))
        with self.assertRaises(ec.Refused):     # still a request: refused, never swapped for Sage
            ec.choose("DDA", "timstof", requested="fragpipe", env=ENV_MAC)

    def test_fragpipe_asked_for_where_it_cannot_run_is_refused_naming_sage(self):
        with self.assertRaises(ec.Refused) as cm:
            ec.choose("DDA", "timstof", requested="fragpipe", env=ENV_MAC)
        self.assertIn("fragpipe cannot run on this machine", str(cm.exception))
        self.assertIn("Available here for this data type: sage", str(cm.exception))

    def test_fragpipe_asked_for_without_the_licensed_jars_is_refused_with_why(self):
        with self.assertRaises(ec.Refused) as cm:
            ec.choose("DDA", "timstof", requested="fragpipe", env=ENV_LINUX, fp_tools=FP_ABSENT)
        self.assertIn(FP_ABSENT["why"], str(cm.exception))

    def test_tools_not_inspected_does_not_block_the_default(self):
        """No tools.json and no probe: FragPipe stays; run_search.py checks tools.json again."""
        sel, st = ec.choose("DDA", "timstof", env=ENV_LINUX, fp_tools=None)
        self.assertEqual((sel["engine"], sel["how"]), ("fragpipe", "default"))
        self.assertIn("not inspected", st["why"])

    def test_no_client_sentence_carries_a_version_path_or_licence_word(self):
        cases = [("DDA", "timstof", None, ENV_MAC, FP_ABSENT), ("DDA", "timstof", None, ENV_LINUX, FP_ABSENT),
                 ("DDA", "timstof", "sage", ENV_LINUX, FP_LICENSED),
                 ("DIA", "timstof", "fragpipe", ENV_LINUX, FP_LICENSED)]
        for acq, fam, req, e, ft in cases:
            sel, _ = ec.choose(acq, fam, requested=req, env=e, fp_tools=ft)
            blob = json.dumps(sel)
            self.assertIsNone(CLIENT_FORBIDDEN.search(sel["sentence"] or ""), sel["sentence"])
            for word in ("unavailable_here", "not_licensed_here", "why", "/x/tools"):
                self.assertNotIn(word, blob)


class OneTable(unittest.TestCase):
    def test_routes_are_one_table(self):
        self.assertIs(rd.ROUTES, ec.ROUTES)
        self.assertFalse(hasattr(rd, "DE_METHOD"), "DE_METHOD is engine_choice.de_method now")

    def test_de_method_by_engine_and_acquisition(self):
        self.assertEqual({(e, a): ec.de_method(e, a) for e in ("diann", "radiant", "sage", "fragpipe")
                          for a in ("DIA", "DDA")},
                         {("diann", "DIA"): "dpc", ("diann", "DDA"): "dpc",
                          ("radiant", "DIA"): "maxlfq", ("radiant", "DDA"): "maxlfq",
                          ("sage", "DIA"): "maxlfq", ("sage", "DDA"): "maxlfq",
                          ("fragpipe", "DIA"): "dpc", ("fragpipe", "DDA"): "maxlfq"})

    def test_every_route_engine_is_pinned(self):
        for engines in ec.ROUTES.values():
            for e in engines:
                self.assertIn(e, rd.ENGINE_VERSIONS)


class ResolveHarness(unittest.TestCase):
    """resolve_defaults.py as step 4 runs it (the helpers; no tests)."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, name, obj):
        p = os.path.join(self.d, name)
        with open(p, "w") as fh:
            json.dump(obj, fh)
        return p

    def resolve(self, *args, name="wf", ok=True):
        dest = os.path.join(self.d, name)
        envv = {k: v for k, v in os.environ.items() if k not in (TOOLS_ENV,)}
        envv[CORE_ENV] = os.path.join(self.d, "no-core-install")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            "--dest", dest, *args], capture_output=True, text=True, timeout=120,
                           env=envv)
        if not ok:
            return p, None
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(dest, "workflow.manifest.json")) as fh:
            return json.loads(p.stdout), json.load(fh)

class ResolveDefaults(ResolveHarness):
    """resolve_defaults.py as step 4 runs it."""

    def test_laptop_env_gives_sage_with_a_client_sentence_and_a_staff_reason(self):
        out, man = self.resolve("--acquisition", "DDA", "--instrument", "timsTOF Pro 2",
                                "--env", self.write("env.json", ENV_MAC))
        sel = man["engine"]["selection"]
        self.assertEqual((man["engine"]["name"], man["de"]["method"], man["alternatives"]),
                         ("sage", "maxlfq", ["fragpipe"]))
        self.assertEqual((sel["how"], sel["default"], sel["staff_record"]),
                         ("fallback", "fragpipe", "workflow.manifest.json.staff.json"))
        self.assertEqual(out["selection"], sel)
        rec = staff.read(os.path.join(self.d, "wf", "workflow.manifest.json" + staff.STAFF_SUFFIX))
        self.assertEqual([(e["event"], e["reason"]) for e in rec],
                         [("engine_selection", "unavailable_here")])
        # nothing of the reason in the manifest a client may receive
        blob = json.dumps(man)
        for word in ("unavailable_here", "macOS", "licen"):
            self.assertNotIn(word, blob)

    def test_tools_json_without_the_jars_gives_sage_not_licensed(self):
        tools = self.write("tools.json", {"fragpipe": "/x/bin/fragpipe", "fragpipe_tools": FP_ABSENT})
        _, man = self.resolve("--acquisition", "DDA", "--instrument", "timsTOF Pro 2",
                              "--env", self.write("env.json", ENV_LINUX), "--tools", tools)
        self.assertEqual((man["engine"]["name"], man["engine"]["selection"]["how"]), ("sage", "fallback"))
        rec = staff.read(os.path.join(self.d, "wf", "workflow.manifest.json" + staff.STAFF_SUFFIX))
        self.assertEqual(rec[-1]["reason"], "not_licensed_here")

    def test_tools_json_licensed_but_not_ready_gives_sage_not_ready(self):
        """Review M2: step 4 with a tools.json whose FragPipe is licensed but cannot run picks the
        Sage fallback (the manifest), the reason in the staff record only."""
        tools = self.write("tools.json", {"fragpipe": "/x/bin/fragpipe", "fragpipe_tools": FP_NOT_READY})
        _, man = self.resolve("--acquisition", "DDA", "--instrument", "timsTOF Pro 2",
                              "--env", self.write("env.json", ENV_LINUX), "--tools", tools)
        self.assertEqual((man["engine"]["name"], man["engine"]["selection"]["how"],
                          man["de"]["method"]), ("sage", "fallback", "maxlfq"))
        rec = staff.read(os.path.join(self.d, "wf", "workflow.manifest.json" + staff.STAFF_SUFFIX))
        self.assertEqual(rec[-1]["reason"], "not_ready_here")
        self.assertNotIn("not_ready_here", json.dumps(man))

    def test_ids_only_is_refused_off_dda_and_says_what_to_do(self):
        p, _ = self.resolve("--acquisition", "DIA", "--instrument", "timsTOF HT", "--ids-only",
                            ok=False)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("--ids-only applies to DDA searches only (Sage or FragPipe)", p.stderr)
        self.assertIn("the search's own report is the result", p.stderr)

    def test_a_replayed_fallback_keeps_the_original_choice(self):
        """Ruling 103: reproduce.sh replays step 4 with --engine sage. Without the original
        manifest a Sage FALLBACK came back `requested` -- a different story about why Sage ran. The
        bundle's step 2 passes --replay-of its copy of the original; the choice is kept, and the
        replay says it was replayed."""
        _, orig = self.resolve("--acquisition", "DDA", "--instrument", "timsTOF Pro 2",
                               "--env", self.write("env_laptop.json", ENV_MAC), name="orig")
        self.assertEqual(orig["engine"]["selection"]["how"], "fallback")
        bundle = os.path.join(self.d, "bundle")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                            "--outdir", bundle, "--workflow-manifest",
                            os.path.join(self.d, "orig", "workflow.manifest.json"),
                            "--engine", "sage", "--acquisition", "DDA",
                            "--instrument", "timsTOF Pro 2"], capture_output=True, text=True,
                           timeout=300)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(os.path.join(bundle, "reproduce.sh")) as fh:
            sh = fh.read()
        start = sh.index('python3 "$SKILL/scripts/resolve_defaults.py"')
        step2 = sh[start:sh.index("--dest ./wf", start) + len("--dest ./wf")]
        self.assertIn("--replay-of ./inputs/workflow.manifest.json", step2)
        # replay it where FragPipe WOULD run: the original choice must survive anyway
        self.write("bundle/env.json", ENV_LINUX)
        envv = dict(os.environ, SKILL=os.path.dirname(SCRIPTS))
        envv[CORE_ENV] = os.path.join(self.d, "no-core-install")
        r = subprocess.run(["bash", "-c", step2.replace("python3 ", f"{sys.executable} ", 1)],
                           cwd=bundle, capture_output=True, text=True, env=envv, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        with open(os.path.join(bundle, "wf", "workflow.manifest.json")) as fh:
            replayed = json.load(fh)["engine"]["selection"]
        # review L6: the original's staff record stayed with the original session -- the replay
        # keeps the choice but names no staff file that is not in its own folder
        self.assertEqual(orig["engine"]["selection"]["staff_record"],
                         "workflow.manifest.json.staff.json")
        self.assertEqual(replayed, dict(orig["engine"]["selection"],
                                        replayed="replayed with --engine sage",
                                        staff_record=None))

    def test_requested_fragpipe_that_cannot_run_exits_as_before(self):
        p, _ = self.resolve("--acquisition", "DDA", "--instrument", "timsTOF Pro 2", "--engine",
                            "fragpipe", "--env", self.write("env.json", ENV_MAC), ok=False)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("resolve_defaults: fragpipe cannot run on this machine", p.stderr)

    def test_requested_fragpipe_is_refused_on_tools_json_only_never_on_a_probe_here(self):
        """A workflow may be prepared here for a search elsewhere (CI, a laptop): with no
        tools.json, nothing installed here refuses a REQUESTED FragPipe (as in 2.14.2)."""
        _, man = self.resolve("--acquisition", "DIA", "--instrument", "timsTOF HT",
                              "--engine", "fragpipe")
        self.assertEqual((man["engine"]["name"], man["engine"]["selection"]["how"]),
                         ("fragpipe", "requested"))
        rec = staff.read(os.path.join(self.d, "wf", "workflow.manifest.json" + staff.STAFF_SUFFIX))
        self.assertIsNone(rec[-1]["checks"]["tools"])
        tools = self.write("tools.json", {"fragpipe_tools": FP_ABSENT})
        p, _ = self.resolve("--acquisition", "DIA", "--instrument", "timsTOF HT", "--engine",
                            "fragpipe", "--tools", tools, name="wf2", ok=False)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("MSFragger/IonQuant are not installed for this account", p.stderr)

    def test_timstof_dda_fragpipe_writes_its_preset_without_a_fasta(self):
        tools = self.write("tools.json", {"fragpipe_tools": FP_LICENSED})
        out, man = self.resolve("--acquisition", "DDA", "--instrument", "timsTOF Pro 2",
                                "--env", self.write("env.json", ENV_LINUX), "--tools", tools,
                                "--ids-only")
        self.assertEqual(man["engine"]["name"], "fragpipe")
        self.assertIsNone(man["de"], "an identifications-only search carries no DE block (X8 F5)")
        self.assertTrue(os.path.isfile(man["search"]["params_file"]))
        self.assertNotIn("preset_pending", json.dumps(man))
        self.assertIs(man["engine"]["selection"]["ids_only"], True)


class UnitResolutionHasNoAmounts(ResolveHarness):
    """Review M3 / Ruling 107.3: a unit-resolution instrument's MS1 cannot support label-free MS1
    quantities -- a fact about the DATA, so the refusal (Ruling 105.2) is ONE engine-independent
    rule, applied before the engine is chosen: the default Sage route of a Stellar DDA is refused
    exactly as a FragPipe one, with the same words, and nothing is written. --ids-only searches
    it; a DIA search is not covered (DIA-NN quantifies from fragments)."""
    WORDS = ("reads MS1 at unit resolution", "--ids-only", "PRM", "Skyline")

    def refused(self, *args, name):
        p, _ = self.resolve("--acquisition", "DDA", *args, name=name, ok=False)
        self.assertEqual(p.returncode, 1, p.stdout + p.stderr)
        for w in self.WORDS:
            self.assertIn(w, p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.d, name)), "something was written")
        return p.stderr.strip().splitlines()[-1]

    def test_the_default_sage_route_is_refused_too(self):
        for inst in ("Orbitrap Stellar", "Stellar"):
            with self.subTest(inst):
                msg = self.refused("--instrument", inst, "--env", self.write("env.json", ENV_LINUX),
                                   name="wf_" + inst.replace(" ", "_"))
                self.assertIn(inst, msg)

    def test_both_engines_give_the_same_words(self):
        tools = self.write("tools.json", {"fragpipe_tools": FP_LICENSED})
        envp = self.write("env.json", ENV_LINUX)
        sage = self.refused("--instrument", "Stellar", "--engine", "sage", "--env", envp,
                            name="wf_sage")
        frag = self.refused("--instrument", "Stellar", "--engine", "fragpipe", "--env", envp,
                            "--tools", tools, name="wf_fragpipe")
        self.assertEqual(sage, frag)
        # Ruling 105's message, unchanged
        self.assertEqual(sage, "resolve_defaults: " + ep.UNIT_RESOLUTION_NO_QUANT.format(
            instrument="Stellar") + ". Nothing was written.")

    def test_ids_only_is_searched_on_either_engine(self):
        _, man = self.resolve("--acquisition", "DDA", "--instrument", "Orbitrap Stellar",
                              "--ids-only", "--env", self.write("env.json", ENV_LINUX))
        self.assertEqual((man["engine"]["name"], man["engine"]["selection"]["ids_only"]),
                         ("sage", True))
        tools = self.write("tools.json", {"fragpipe_tools": FP_LICENSED})
        _, man = self.resolve("--acquisition", "DDA", "--instrument", "Stellar", "--ids-only",
                              "--engine", "fragpipe", "--tools", tools, name="wf2")
        self.assertEqual(man["engine"]["name"], "fragpipe")

    def test_the_sage_cfg_step_refuses_it_too(self):
        """step 6b (estimate_params.py --engine sage) reads the same rule (dda_args)"""
        out = os.path.join(self.d, "sage_config.json")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "sage", "--acquisition", "DDA", "--instrument", "Stellar",
                            "--out", out], capture_output=True, text=True, timeout=120)
        self.assertNotEqual(p.returncode, 0)
        for w in self.WORDS:
            self.assertIn(w, p.stdout + p.stderr)
        self.assertFalse(os.path.exists(out))
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "sage", "--acquisition", "DDA", "--instrument", "Stellar",
                            "--ids-only", "--out", out], capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out) as fh:
            self.assertIs(json.load(fh)["quant"]["lfq"], False)

    def test_dia_and_high_resolution_dda_are_not_touched(self):
        _, man = self.resolve("--acquisition", "DIA", "--instrument", "Stellar")
        self.assertEqual(man["engine"]["name"], "diann")
        for inst in ("Orbitrap Fusion Lumos", "timsTOF Pro 2", "Orbitrap Astral"):
            self.assertIsNone(ep.unit_resolution_amounts_refused(inst, "DDA", False), inst)
        self.assertIsNone(ep.unit_resolution_amounts_refused("Stellar", "DIA", False))
        self.assertIsNone(ep.unit_resolution_amounts_refused("Stellar", "DDA", True))

    def test_one_definition(self):
        import make_presets as mp
        self.assertFalse(hasattr(mp, "UNIT_RESOLUTION_NO_QUANT") and
                         mp.UNIT_RESOLUTION_NO_QUANT is not ep.UNIT_RESOLUTION_NO_QUANT,
                         "make_presets keeps a second copy of the refusal")


class Pins(unittest.TestCase):
    """The Lumos DDA and the DIA manifests: what 2f9b229c wrote, plus engine.selection only."""
    CASES = {
        "lumos_dda": ["--acquisition", "DDA", "--instrument", "Orbitrap Fusion Lumos"],
        "lumos_dda_itms": ["--acquisition", "DDA", "--instrument", "Orbitrap Fusion Lumos",
                           "--ms1-resolution", "120000", "--ms2-analyzer", "ITMS"],
        "timstofht_dia": ["--acquisition", "DIA", "--instrument", "timsTOF HT"],
        "astral_dia": ["--acquisition", "DIA", "--instrument", "Orbitrap Astral"],
        "exploris_dia": ["--acquisition", "DIA", "--instrument", "Orbitrap Exploris 480",
                         "--ms1-resolution", "120000", "--ms2-resolution", "30000",
                         "--resolution-source", "detected"],
    }

    def manifests(self):
        """{case: (manifest without engine.selection, the selection or None, staff record?)}"""
        out = {}
        for name, args in self.CASES.items():
            with tempfile.TemporaryDirectory() as d:
                dest = os.path.join(d, "out")
                p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                                    *args, "--dest", dest], capture_output=True, text=True,
                                   timeout=120)
                self.assertEqual(p.returncode, 0, p.stderr)
                with open(os.path.join(dest, "workflow.manifest.json")) as fh:
                    man = json.load(fh)
                sel = man["engine"].pop("selection", None)
                got = json.dumps(man, indent=2, sort_keys=True).replace(dest, "<DEST>") + "\n"
                out[name] = (got, sel, os.path.exists(os.path.join(
                    dest, "workflow.manifest.json" + staff.STAFF_SUFFIX)))
        return out

    def test_manifests_unchanged_but_for_the_selection_record(self):
        """(a pin: passes on 2f9b229c too)"""
        for name, (got, _sel, staff_file) in self.manifests().items():
            with self.subTest(name), open(os.path.join(PINS, f"{name}.json")) as fh:
                self.assertEqual(got, fh.read())
                self.assertFalse(staff_file)

    def test_lumos_dda_with_enzyme_none_routes_to_fragpipe(self):
        """Ruling 120 (replaces Ruling 99's pin): --enzyme none on a Lumos DDA search is FragPipe's
        Nonspecific-peptidome where FragPipe can run and is licensed; Sage, recorded as the
        fallback, where it is not. Without --enzyme the Lumos DDA manifest is the pin above."""
        fp_ok = {"fragpipe_tools": dict(FP_LICENSED, binary="/x/bin/fragpipe")}
        fp_no = {"fragpipe_tools": FP_ABSENT}
        for tools, engine, how in ((fp_ok, "fragpipe", "default"), (fp_no, "sage", "fallback")):
            with self.subTest(engine), tempfile.TemporaryDirectory() as d:
                tp = os.path.join(d, "tools.json")
                with open(tp, "w") as fh:
                    json.dump(tools, fh)
                dest = os.path.join(d, "out")
                p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                                    *self.CASES["lumos_dda"], "--enzyme", "none", "--tools", tp, "--dest", dest],
                                   capture_output=True, text=True, timeout=120)
                self.assertEqual(p.returncode, 0, p.stderr)
                with open(os.path.join(dest, "workflow.manifest.json")) as fh:
                    man = json.load(fh)
                sel = man["engine"]["selection"]
                self.assertEqual((man["engine"]["name"], sel["default"], sel["how"]),
                                 (engine, "fragpipe", how))
                self.assertEqual(man["alternatives"], ["sage"] if engine == "fragpipe"
                                 else ["fragpipe"])

    def test_each_carries_a_default_selection_with_no_sentence(self):
        for name, (got, sel, _) in self.manifests().items():
            with self.subTest(name):
                self.assertIsNotNone(sel)
                self.assertEqual((sel["engine"], sel["how"], sel["sentence"], sel["staff_record"]),
                                 (json.loads(got)["engine"]["name"], "default", None, None))

if __name__ == "__main__":
    unittest.main()
