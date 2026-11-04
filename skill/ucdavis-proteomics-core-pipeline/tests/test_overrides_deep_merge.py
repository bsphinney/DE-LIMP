#!/usr/bin/env python3
"""
H2-4: estimate_params.py --overrides merges into a Sage cfg at every depth, and the rationale
describes the cfg as written.

The merge took top-level keys only, so `{"database": {"enzyme": {...}}}` REPLACED the whole
`database` block: the fixed carbamidomethyl-C, the protein N-term acetyl, the decoys and the
FASTA placeholder were gone, with no message. The rationale's enzyme entry still said
"trypsin/P, 2 missed cleavages" whatever the cfg held.

Now (estimate_params.deep_merge, the one merge): a dict that lands on a dict merges key by key;
anything else replaces -- a list such as isotope_errors is replaced, never appended, and a Sage
tolerance ({"ppm": ...} or {"da": ...}, a one-key enum) is replaced whole. Every rationale entry
that describes a cfg value reads it from the merged cfg, and an entry an override changed says so.

Hermetic: build_sage() in-process, and estimate_params.py as a subprocess on a temp dir.
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import estimate_params as ep  # noqa: E402
import staff  # noqa: E402

LYSC = {"database": {"enzyme": {"cleave_at": "K"}}}


def sage(overrides, acq="DDA", cls="orbitrap_generic"):
    text, rationale = ep.build_sage(acq, cls, None, overrides)
    return json.loads(text), rationale


class NoStaffList(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.tmp, True)
        patch = mock.patch.dict(os.environ, {
            staff.STAFF_FILE_ENV: os.path.join(self.tmp, "no_staff_list.txt")})
        patch.start()
        self.addCleanup(patch.stop)


class DeepMerge(NoStaffList):
    def test_an_enzyme_only_override_keeps_the_rest_of_the_database_block(self):
        base, _ = sage({})
        cfg, _ = sage(LYSC)
        db = cfg["database"]
        self.assertEqual(db["enzyme"], dict(ep.SAGE_ENZYME, cleave_at="K"))
        self.assertEqual(db["static_mods"], {"C": 57.0215}, "the fixed carbamidomethyl-C stays")
        self.assertEqual(db["variable_mods"], base["database"]["variable_mods"])
        self.assertEqual(db["variable_mods"]["["], [42.0106], "the protein N-term acetyl stays")
        self.assertIs(db["generate_decoys"], True)
        self.assertEqual(db["decoy_tag"], "rev_")
        self.assertEqual(db["fasta"], "REPLACED_AT_RUNTIME.fasta")
        self.assertEqual({k: v for k, v in db.items() if k != "enzyme"},
                         {k: v for k, v in base["database"].items() if k != "enzyme"})

    def test_a_mod_override_adds_to_the_mods_it_lands_on(self):
        base, _ = sage({})
        cfg, _ = sage({"database": {"variable_mods": {"S": [79.9663]}}})
        self.assertEqual(cfg["database"]["variable_mods"],
                         dict(base["database"]["variable_mods"], S=[79.9663]))

    def test_a_list_is_replaced_not_appended(self):
        """(pin) a list value replaces, as the top-level merge already did."""
        cfg, _ = sage({"isotope_errors": [-1, 3]})
        self.assertEqual(cfg["isotope_errors"], [-1, 3])

    def test_a_tolerance_is_replaced_whole(self):
        """(pin) {"da": ...} over a ppm window gives a Da window, never both units: a Sage
        Tolerance is a one-key enum."""
        cfg, _ = sage({"fragment_tol": {"da": [-0.5, 0.5]}})
        self.assertEqual(cfg["fragment_tol"], {"da": [-0.5, 0.5]})

    def test_deep_merge_leaves_its_inputs_alone_and_lists_what_it_set(self):
        base = {"a": {"b": 1, "c": [1, 2]}, "d": 2}
        over = {"a": {"c": [3]}, "e": {"f": 1}}
        merged, set_ = ep.deep_merge(base, over)
        self.assertEqual(merged, {"a": {"b": 1, "c": [3]}, "d": 2, "e": {"f": 1}})
        self.assertEqual(base, {"a": {"b": 1, "c": [1, 2]}, "d": 2})
        self.assertEqual(sorted(set_), [(("a", "c"), [3]), (("e",), {"f": 1})])


class TheRationaleReadsTheFinalCfg(NoStaffList):
    def test_the_rationale_names_the_overrides_enzyme(self):
        cfg, r = sage(LYSC)
        self.assertEqual(r["enzyme"]["value"], cfg["database"]["enzyme"])
        self.assertEqual(r["enzyme"]["value"]["cleave_at"], "K")
        self.assertIn("user override: database.enzyme.cleave_at = K", r["enzyme"]["source"])
        self.assertNotIn("trypsin/P", json.dumps(r["enzyme"]))
        self.assertIn("user override", r["database.enzyme.cleave_at"]["source"])

    def test_an_entry_no_override_touched_is_unchanged(self):
        cfg, r = sage(LYSC)
        self.assertEqual(r["static_mods"]["value"], {"C": 57.0215})
        self.assertNotIn("user override", r["static_mods"]["source"])
        self.assertNotIn("user override", r["variable_mods"]["source"])

    def test_the_default_enzyme_entry_is_the_cfg_block(self):
        cfg, r = sage({})
        self.assertEqual(r["enzyme"]["value"], cfg["database"]["enzyme"])
        self.assertNotIn("trypsin/P", json.dumps(r["enzyme"]),
                         "Sage's restrict 'P' is trypsin, not trypsin/P")

    def test_a_mod_override_shows_the_merged_mods(self):
        cfg, r = sage({"database": {"variable_mods": {"S": [79.9663]}}})
        self.assertEqual(r["variable_mods"]["value"], cfg["database"]["variable_mods"])
        self.assertIn("user override: database.variable_mods.S", r["variable_mods"]["source"])

    def test_a_quant_override_shows_the_final_lfq(self):
        cfg, r = sage({"quant": {"lfq": False}})
        self.assertIs(cfg["quant"]["lfq"], False)
        self.assertIs(r["lfq"]["value"], False)
        self.assertIn("user override: quant.lfq", r["lfq"]["source"])

    def test_a_top_level_override_still_replaces_its_derived_entry(self):
        """(pin) precursor_tol drops precursor_tol_ppm, as tests/test_overrides_are_the_users.py
        pins through the CLI."""
        tol = {"ppm": [-25.0, 25.0]}
        cfg, r = sage({"precursor_tol": tol})
        self.assertEqual(cfg["precursor_tol"], tol)
        self.assertNotIn("precursor_tol_ppm", r)
        self.assertEqual(r["precursor_tol"]["value"], tol)


class TheCommandLine(NoStaffList):
    def test_an_enzyme_override_through_the_cli_keeps_c57_and_decoys(self):
        out = os.path.join(self.tmp, "sage.json")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "sage", "--acquisition", "DDA",
                            "--instrument", "Orbitrap Exploris 480", "--out", out,
                            "--overrides", json.dumps({"database": {"enzyme": {
                                "missed_cleavages": 1}}}),
                            "--override-reason", "one missed cleavage for this study"],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out) as fh:
            db = json.load(fh)["database"]
        with open(out + ".rationale.json") as fh:
            side = json.load(fh)
        self.assertEqual(db["enzyme"], dict(ep.SAGE_ENZYME, missed_cleavages=1))
        self.assertEqual(db["static_mods"], {"C": 57.0215})
        self.assertIs(db["generate_decoys"], True)
        self.assertEqual(db["variable_mods"]["["], [42.0106])
        self.assertEqual(side["rationale"]["enzyme"]["value"]["missed_cleavages"], 1)
        self.assertIn("one missed cleavage for this study",
                      side["rationale"]["enzyme"]["source"])


class AnEmptyBlockIsRefused(NoStaffList):
    """An empty dict merged into a block changes nothing, so an override of {} -- meant, say, to
    clear the fixed modifications -- would be dropped without a word. It is refused, saying what
    to give instead (H2-4 follow-up)."""

    def test_deep_merge_refuses_an_empty_dict_on_a_dict(self):
        with self.assertRaises(ep.EmptyOverride) as cm:
            ep.deep_merge({"database": {"static_mods": {"C": 57.0215}}},
                          {"database": {"static_mods": {}}})
        self.assertIn("database.static_mods", str(cm.exception))

    def test_an_empty_dict_on_a_new_key_is_a_value(self):
        """(pin of the rule) nothing to merge into: it is set."""
        merged, set_ = ep.deep_merge({"a": 1}, {"b": {}})
        self.assertEqual(merged, {"a": 1, "b": {}})
        self.assertEqual(set_, [(("b",), {})])

    def test_the_cli_refuses_and_writes_no_cfg(self):
        out = os.path.join(self.tmp, "sage.json")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "sage", "--acquisition", "DDA",
                            "--instrument", "Orbitrap Exploris 480", "--out", out,
                            "--overrides", json.dumps({"database": {"static_mods": {}}})],
                           capture_output=True, text=True, timeout=60)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("database.static_mods", p.stderr)
        self.assertIn("--cys-alkylation none", p.stderr)
        self.assertIn("No cfg was written", p.stderr)
        self.assertFalse(os.path.exists(out))


class AKeyAtTheWrongLevelIsRefused(NoStaffList):
    """Review F2: the rationale names its entries flat (enzyme, static_mods, variable_mods, lfq),
    and an override written that way -- {"enzyme": {...}} at the top level -- landed where Sage
    0.14.7 ignores an unknown field (no `deny_unknown_fields` in crates/sage-cli/src/input.rs or
    crates/sage/src/database.rs), while the rationale, the sidecar and the run record said it was
    applied. A key that is not a Sage setting at its level is refused, naming where it belongs."""

    def refused(self, overrides):
        with self.assertRaises(ep.OverrideRefused) as cm:
            ep.build_sage("DDA", "orbitrap_generic", None, overrides)
        return str(cm.exception)

    def test_a_misplaced_key_names_its_nesting(self):
        for over, where in (({"enzyme": {"cleave_at": "K"}}, "database.enzyme"),
                            ({"lfq": False}, "quant.lfq"),
                            ({"static_mods": {"C": 0}}, "database.static_mods"),
                            ({"variable_mods": {"S": [79.9663]}}, "database.variable_mods"),
                            ({"database": {"cleave_at": "K"}}, "database.enzyme.cleave_at"),
                            ({"missed_cleavages": 1}, "database.enzyme.missed_cleavages")):
            msg = self.refused(over)
            self.assertIn(f"belongs at {where}", msg, over)

    def test_a_key_sage_does_not_have_is_refused(self):
        msg = self.refused({"databse": {"enzyme": {"cleave_at": "K"}}})
        self.assertIn("databse", msg)
        self.assertIn("not a Sage 0.14.7 setting", msg)
        self.assertIn("quant.lfq_bogus", self.refused({"quant": {"lfq_bogus": 1}}))

    def test_a_block_given_a_plain_value_is_refused(self):
        self.assertIn("database is a block", self.refused({"database": 5}))

    def test_sage_settings_the_default_cfg_does_not_write_are_kept(self):
        """(pin) real Sage 0.14.7 keys outside the generated cfg still merge."""
        over = {"report_psms": 2, "precursor_charge": [2, 4],
                "quant": {"lfq_settings": {"ppm_tolerance": 10.0}},
                "database": {"enzyme": {"semi_enzymatic": True}}}
        cfg, _ = sage(over)
        self.assertEqual((cfg["report_psms"], cfg["precursor_charge"]), (2, [2, 4]))
        self.assertEqual(cfg["quant"], {"lfq": True, "lfq_settings": {"ppm_tolerance": 10.0}})
        self.assertIs(cfg["database"]["enzyme"]["semi_enzymatic"], True)

    def test_the_cli_refuses_and_records_nothing(self):
        out = os.path.join(self.tmp, "sage.json")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "estimate_params.py"),
                            "--engine", "sage", "--acquisition", "DDA",
                            "--instrument", "Orbitrap Exploris 480", "--out", out,
                            "--overrides", json.dumps({"enzyme": {"cleave_at": "K"}})],
                           capture_output=True, text=True, timeout=60)
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("belongs at database.enzyme", p.stderr)
        self.assertIn("No cfg was written", p.stderr)
        for f in (out, out + ".rationale.json", out + ".staff.json"):
            self.assertFalse(os.path.exists(f), f)


class AToleranceMustBeOneWindow(NoStaffList):
    """Review F5: a Sage Tolerance is exactly one of {"ppm": [lo, hi]} or {"da": [lo, hi]}
    (crates/sage/src/mass.rs). {"precursor_tol": {}} wrote a cfg Sage cannot read, and the job
    failed only after it had queued. An empty or malformed window is refused."""

    def test_empty_or_malformed_windows_are_refused(self):
        for tol in ({}, {"ppm": [-10, 10], "da": [-0.5, 0.5]}, {"ppm": [10]}, {"mda": [-1, 1]},
                    {"ppm": "10"}, {"ppm": [10, -10]}, {"ppm": [True, 1]}, 25, [-10, 10]):
            for key in ("precursor_tol", "fragment_tol"):
                with self.assertRaises(ep.OverrideRefused, msg=f"{key}={tol!r}") as cm:
                    ep.build_sage("DDA", "orbitrap_generic", None, {key: tol})
                self.assertIn(key, str(cm.exception))
                self.assertIn('{"ppm": [lo, hi]}', str(cm.exception))

    def test_a_window_in_either_unit_is_kept(self):
        """(pin)"""
        for tol in ({"ppm": [-25.0, 25.0]}, {"da": [-0.5, 0.5]}, {"ppm": [-5, 20]}):
            cfg, _ = sage({"fragment_tol": tol})
            self.assertEqual(cfg["fragment_tol"], tol)


class TheUntouchedPartKeepsItsTag(NoStaffList):
    """Review F4 (rule 2): an override of part of an entry -- one residue of static_mods, the
    missed cleavages of the enzyme -- replaced the entry's source, so the rest of it (C+57,
    "DEFAULT, not user-confirmed"; trypsin, the universal default) lost its tag. The entry keeps
    its own source and adds the override's; an override of the whole value says only its own."""

    def test_part_of_static_mods(self):
        _, r = sage({"database": {"static_mods": {"K": 8.0142}}})
        src = r["static_mods"]["source"]
        self.assertEqual(r["static_mods"]["value"], {"C": 57.0215, "K": 8.0142})
        self.assertIn("DEFAULT, not user-confirmed", src)
        self.assertIn("+ user override: database.static_mods.K = 8.0142", src)
        self.assertLess(src.index("DEFAULT"), src.index("user override"))

    def test_part_of_the_enzyme(self):
        _, r = sage({"database": {"enzyme": {"missed_cleavages": 1}}})
        src = r["enzyme"]["source"]
        self.assertIn("universal trypsin/LFQ default", src)
        self.assertIn("+ user override: database.enzyme.missed_cleavages = 1", src)

    def test_a_user_set_c_is_not_tagged_default(self):
        """Review R2: {"C": 0} replaces every key the DEFAULT described (static_mods held only C),
        so nothing of it is a default any more: the source is the user's alone."""
        _, r = sage({"database": {"static_mods": {"C": 0}}})
        self.assertEqual(r["static_mods"]["value"], {"C": 0})
        self.assertNotIn("DEFAULT", r["static_mods"]["source"])
        self.assertTrue(r["static_mods"]["source"].startswith(
            "user override: database.static_mods.C = 0"))

    def test_the_whole_value_says_only_the_override(self):
        _, r = sage({"quant": {"lfq": False}})
        self.assertNotIn("label-free quantification", r["lfq"]["source"])
        self.assertTrue(r["lfq"]["source"].startswith("user override: quant.lfq = "))


if __name__ == "__main__":
    unittest.main()
