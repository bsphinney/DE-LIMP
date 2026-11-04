#!/usr/bin/env python3
"""
--var-mods: an unknown word is refused, `deam` is supported, and Met-ox is on by default for a
Sage DDA search (I47, I47b; Core decision 2026-10-30, Q2).

A "like-for-like" re-search of a FragPipe study (2026-10-29) asked for N/Q deamidation and got
none: estimate_params.py acted on `ox` only and dropped every other word without a word. The
words now come from one table (estimate_params.VAR_MODS) and anything else stops the script
before a file is written.

Met-ox: a Sage DDA search with no --var-mods now searches oxidised methionine, as FragPipe and
MaxQuant do by default; `--var-mods none` turns it off. A DIA-NN search keeps none by default
(DIA-NN's README advises against variable mods for quantification). The Methods read the mods
from the parameters file, so they name Oxidation (M) for the one and none for the other.

Hermetic: estimate_params.py as a subprocess on a temp dir; make_methods.py in-process.
"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import estimate_params as ep  # noqa: E402


def estimate(d, engine, *extra, acquisition="DIA", ok=True):
    """Run estimate_params.py -> (cfg path, rationale sidecar); with ok=False, the process."""
    out = os.path.join(d, "params.cfg" if engine == "diann" else "sage.json")
    p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                        "--engine", engine, "--acquisition", acquisition,
                        "--instrument", "timsTOF HT",
                        "--precursor-mz-range", "299.5", "1200.5", "--out", out, *extra],
                       capture_output=True, text=True, timeout=60)
    if not ok:
        return p
    if p.returncode != 0:
        raise AssertionError(p.stderr)
    with open(out + ".rationale.json") as fh:
        return out, json.load(fh)


def cfg_lines(path):
    with open(path) as fh:
        return fh.read().splitlines()


def sage_vm(path):
    with open(path) as fh:
        return json.load(fh)["database"]["variable_mods"]


class Table(unittest.TestCase):
    def test_the_table(self):
        self.assertEqual(ep.VAR_MODS["ox"]["diann"], "UniMod:35,15.994915,M")
        self.assertEqual(ep.VAR_MODS["ox"]["sage"], {"M": [15.994915]})
        self.assertEqual(ep.VAR_MODS["deam"]["diann"], "UniMod:7,0.984016,NQ")
        self.assertEqual(ep.VAR_MODS["deam"]["sage"], {"N": [0.984016], "Q": [0.984016]})

    def test_parse(self):
        self.assertEqual(ep.parse_var_mods("OX, deam"), ["ox", "deam"])
        self.assertEqual(ep.parse_var_mods("ox,ox"), ["ox"])
        self.assertEqual(ep.parse_var_mods("none"), [])
        for bad in ("deamid", "ox,phospho", "ox,none"):
            with self.assertRaises(ValueError) as c:
                ep.parse_var_mods(bad)
            self.assertIn("ox, deam, none", str(c.exception))


class EmptyList(unittest.TestCase):
    """`--var-mods ","` has no word in it: it is NOT the explicit `none` (that would switch off
    Sage DDA's default Met-ox without anyone asking); it is blank, so the default applies."""

    def test_parse_gives_none_for_no_words(self):
        for blank in (",", " , ", ",,", ""):
            self.assertIsNone(ep.parse_var_mods(blank), repr(blank))
        self.assertEqual(ep.parse_var_mods("none"), [])           # only `none` is off

    def test_a_comma_only_value_keeps_the_sage_dda_default(self):
        for blank in (",", " , "):
            with tempfile.TemporaryDirectory() as d:
                out, _ = estimate(d, "sage", "--var-mods", blank, acquisition="DDA")
                self.assertEqual(sage_vm(out)["M"], [15.994915], repr(blank))


class ReplacesTheDefault(unittest.TestCase):
    """Giving --var-mods replaces the default (deam alone drops Met-ox): said in --help and in
    references/parameters.md."""
    SENTENCE = ("giving --var-mods replaces the default; to add deamidation and keep met-ox, "
                "give ox,deam")

    def test_help_and_reference_say_so(self):
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"), "--help"],
                           capture_output=True, text=True, timeout=60)
        flat = lambda t: " ".join(t.split()).replace("`", "").lower()   # noqa: E731
        self.assertIn(self.SENTENCE, flat(p.stdout))
        with open(os.path.join(os.path.dirname(HERE), "references", "parameters.md")) as fh:
            self.assertIn(self.SENTENCE, flat(fh.read()))


class Refused(unittest.TestCase):
    def test_unknown_word_is_refused_and_nothing_written(self):
        with tempfile.TemporaryDirectory() as d:
            cfg_path = os.path.join(d, "sage.json")
            # a cfg from an earlier run of the same command: a refusal leaves nothing to read
            for stale in (cfg_path, cfg_path + ".rationale.json"):
                with open(stale, "w") as fh:
                    fh.write("{}")
            p = estimate(d, "sage", "--var-mods", "deamid", acquisition="DDA", ok=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("ox, deam, none", p.stderr)
            self.assertIn("deamid", p.stderr)
            self.assertFalse(os.path.exists(cfg_path))
            self.assertFalse(os.path.exists(cfg_path + ".rationale.json"))

    def test_a_mistyped_word_never_deletes_the_earlier_cfg(self):
        """C2 errors_ux L1: the session's first search ran with the cfg at --out; a re-run of
        step 6 with a typo deleted it (and its rationale), saying only "No cfg was written". Now
        both are set aside (.stale-<time>), never deleted, and the message says where."""
        import glob
        with tempfile.TemporaryDirectory() as d:
            out, _side = estimate(d, "sage", acquisition="DDA")
            before = {}
            for path in (out, out + ".rationale.json"):
                with open(path) as fh:
                    before[path] = fh.read()
            p = estimate(d, "sage", "--var-mods", "ox,deamidation", acquisition="DDA", ok=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("No cfg was written", p.stderr)
            self.assertIn("was NOT deleted: set aside as", p.stderr)
            self.assertFalse(os.path.exists(out))                 # nothing at --out to misread
            for path, text in before.items():
                aside = glob.glob(path + ".stale-*")
                self.assertEqual(len(aside), 1, os.listdir(d))
                with open(aside[0]) as fh:
                    self.assertEqual(fh.read(), text)
                self.assertIn(aside[0], p.stderr)

    def test_unknown_word_is_refused_for_diann_too(self):
        with tempfile.TemporaryDirectory() as d:
            p = estimate(d, "diann", "--var-mods", "ox,phospho", ok=False)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("ox, deam, none", p.stderr)
            self.assertFalse(os.path.exists(os.path.join(d, "params.cfg")))


class DiaNN(unittest.TestCase):
    def test_deam_reaches_diann(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "diann", "--var-mods", "deam")
            lines = cfg_lines(cfg)
        self.assertIn("--var-mod UniMod:7,0.984016,NQ", lines)
        self.assertIn("--var-mods 1", lines)
        self.assertIn("--var-mods deam", side["rationale"]["variable_mods"]["source"])

    def test_ox_and_deam_both_reach_diann_and_count_two_kinds(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "diann", "--var-mods", "ox,deam")
            lines = cfg_lines(cfg)
        self.assertIn("--var-mod UniMod:35,15.994915,M", lines)
        self.assertIn("--var-mod UniMod:7,0.984016,NQ", lines)
        self.assertIn("--var-mods 2", lines)

    def test_diann_dia_default_unchanged(self):
        """(pin) A DIA-NN DIA search with no --var-mods searches none, as before."""
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "diann")
            lines = cfg_lines(cfg)
        self.assertFalse([ln for ln in lines if ln.startswith("--var-mod")], lines)
        vm = side["rationale"]["variable_mods"]
        self.assertEqual(vm["value"], "none")
        self.assertIn("DIA-NN README", vm["source"])


class Sage(unittest.TestCase):
    def test_deam_reaches_sage(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "sage", "--var-mods", "deam", acquisition="DDA")
            vm = sage_vm(cfg)
        self.assertEqual(vm["N"], [0.984016])
        self.assertEqual(vm["Q"], [0.984016])
        self.assertNotIn("M", vm, "deam alone was asked for, so Ox(M) is not added")

    def test_sage_dda_has_met_ox_by_default(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "sage", acquisition="DDA")
            vm = sage_vm(cfg)
        rat = side["rationale"]
        self.assertEqual(vm["M"], [15.994915])
        self.assertEqual(vm["["], [42.0106], "protein N-term acetyl as before")
        self.assertIn("on by default for DDA", rat["variable_mods"]["source"])
        self.assertIn("DEFAULT", rat["variable_mods"]["source"])

    def test_none_turns_it_off(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "sage", "--var-mods", "none", acquisition="DDA")
            vm = sage_vm(cfg)
        self.assertNotIn("M", vm)
        src = side["rationale"]["variable_mods"]["source"]
        self.assertIn("--var-mods none", src)
        self.assertNotIn("DEFAULT", src)

    def test_sage_dia_has_no_met_ox_by_default(self):
        """(pin) The Met-ox default is for DDA only."""
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "sage", acquisition="DIA")
            self.assertNotIn("M", sage_vm(cfg))


class Methods(unittest.TestCase):
    """The Methods read the mods from the parameters file (make_methods.search_record), so a
    default Sage DDA search names Oxidation (M) and a default DIA-NN DIA search names none."""

    def paragraph(self, d, engine, cfg):
        import make_methods
        prov = os.path.join(d, "search_provenance.json")
        with open(prov, "w") as fh:
            json.dump({"engine": engine, "version": "x", "params_file": cfg}, fh)
        return make_methods.search_paragraph(make_methods.search_record(search_prov=prov))

    def test_sage_dda_methods_name_oxidation(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "sage", acquisition="DDA")
            text = self.paragraph(d, "sage", cfg)
        self.assertIn("Variable modifications: Oxidation (M)", text)

    def test_diann_dia_methods_name_no_variable_mods(self):
        """(pin)"""
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "diann")
            text = self.paragraph(d, "diann", cfg)
        self.assertIn("No variable modifications were searched.", text)


PHOS = "UniMod:21,79.966331,STY"


class AnOverrideVarModAdds(unittest.TestCase):
    """Review F3: a --var-mod in --overrides replaced EVERY --var-mod line, so the documented way
    to search a modification not in the table (references/parameters.md: "a modification not in
    the table goes in --overrides") beside --var-mods ox,deam dropped Ox(M) and Deamidation
    without a word, and the rationale still named them. It now adds to the set (union,
    deduplicated), and --var-mods (one per kind) and the rationale follow the final lines."""

    def test_it_adds_to_the_table_mods(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "diann", "--var-mods", "ox,deam",
                                 "--overrides", json.dumps({"--var-mod": PHOS}))
            lines = cfg_lines(cfg)
        vm = [ln for ln in lines if ln.startswith("--var-mod ")]
        self.assertEqual(vm, ["--var-mod UniMod:35,15.994915,M", "--var-mod UniMod:7,0.984016,NQ",
                              f"--var-mod {PHOS}"])
        self.assertIn("--var-mods 3", lines)
        self.assertNotIn("--var-mods 2", lines)
        r = side["rationale"]
        self.assertEqual(r["--var-mod"]["value"], ["UniMod:35,15.994915,M",
                                                   "UniMod:7,0.984016,NQ", PHOS])
        self.assertIn(f"+ user override: --var-mod = {PHOS}", r["--var-mod"]["source"])
        self.assertEqual(r["--var-mods"]["value"], 3)
        self.assertIn(PHOS, r["--var-mods"]["source"])
        self.assertEqual(r["variable_mods"]["value"], ["ox", "deam", PHOS])

    def test_a_list_adds_several_and_a_repeat_is_kept_once(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "diann", "--var-mods", "ox", "--overrides", json.dumps(
                {"--var-mod": [PHOS, "UniMod:35,15.994915,M", PHOS]}))
            vm = [ln for ln in cfg_lines(cfg) if ln.startswith("--var-mod ")]
        self.assertEqual(vm, ["--var-mod UniMod:35,15.994915,M", f"--var-mod {PHOS}"])
        self.assertEqual(side["rationale"]["--var-mods"]["value"], 2)

    def test_with_no_table_mods_it_is_the_only_one(self):
        with tempfile.TemporaryDirectory() as d:
            cfg, side = estimate(d, "diann", "--overrides", json.dumps({"--var-mod": PHOS}))
            lines = cfg_lines(cfg)
        self.assertEqual([ln for ln in lines if ln.startswith("--var-mod ")], [f"--var-mod {PHOS}"])
        self.assertIn("--var-mods 1", lines)
        # review R2: the default's "left off" reason does not describe a mod the user added
        src = side["rationale"]["--var-mods"]["source"]
        self.assertNotIn("left off", src)
        self.assertIn(PHOS, src)
        self.assertIn("--var-mod added in --overrides", src)

    def test_an_explicit_var_mods_count_is_kept(self):
        """--var-mods given in --overrides is the user's maximum, and stays."""
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "diann", "--var-mods", "ox,deam", "--overrides",
                              json.dumps({"--var-mod": PHOS, "--var-mods": 1}))
            lines = cfg_lines(cfg)
        self.assertIn("--var-mods 1", lines)
        self.assertEqual(len([ln for ln in lines if ln.startswith("--var-mods ")]), 1)
        self.assertEqual(len([ln for ln in lines if ln.startswith("--var-mod ")]), 3)


if __name__ == "__main__":
    unittest.main()
