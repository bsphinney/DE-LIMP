#!/usr/bin/env python3
"""
2.15.1 follow-up pass: the cross-lane items the five audit lanes handed on.

  1. normalization_check.py decide accepts the analyst's OWN login for their own analysis (a
     session that is not Core client work), recorded as "the user"; Core client work stays
     staff-only, and so does a check with no session.

Hermetic: temp dirs, synthetic logins.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import experiment_type as et  # noqa: E402
import normalization_check as nc  # noqa: E402
import staff  # noqa: E402
import submission_report as sr  # noqa: E402
from test_staff_decisions import StaffEnv, check_record, REASON  # noqa: E402

ME = staff._uid_login()


class OwnAnalysisDecision(StaffEnv):
    """Item 1: a non-Core user on HIVE (a staff list is readable and they are not on it) records
    the normalisation decision for their own analysis with their own login."""

    def in_session(self, session_json=None):
        s = os.path.join(self.tmp, "S")
        os.makedirs(s, exist_ok=True)
        if session_json is not None:
            with open(os.path.join(s, "session.json"), "w") as fh:
                json.dump(session_json, fh)
        return s, check_record(os.path.join(s, "output", "norm_check", nc.RECORD), s)

    def test_own_login_decides_own_analysis_as_the_user(self):
        self.staff_list("jdoe")                         # readable list; this account not on it
        self.assertFalse(staff.check_by(ME)["ok"])
        s, path = self.in_session({"created": "2026-11-17"})
        d = nc.decide(path, et.RAW, ME, REASON)["decision"]
        self.assertEqual(d["by"], staff.USER_ROLE)
        self.assertTrue(d["statement"].startswith("the user chose non-normalised (raw)"), d)
        self.assertNotIn("Core staff", d["statement"])
        e = staff.read(os.path.join(s, "logs", nc.STAFF_RECORD))[-1]
        self.assertEqual((e["by"], e["role"]), (ME.lower(), staff.USER_ROLE))
        self.assertIn("own analysis", e["staff_check"])
        ok, _, _ = nc.gate(path, et.RAW)
        self.assertTrue(ok)

    def test_someone_elses_login_is_still_refused(self):
        self.staff_list("jdoe")
        _, path = self.in_session()
        with self.assertRaises(ValueError) as cm:
            nc.decide(path, et.RAW, "bjones" if ME != "bjones" else "cjones", REASON)
        self.assertIn("is not on the Core staff list", str(cm.exception))
        self.assertIn("your own login", str(cm.exception))
        self.assertIsNone(nc.load(path)["decision"])

    def test_core_client_work_stays_staff_only(self):
        self.staff_list("jdoe")
        for sj in ({"analysis_for": {"submissions": [], "internal": True}},
                   {"coreomics": {"internal_id": "PROT_9001"}},
                   {"core_client": {"marked": "2026-11-17T00:00:00Z"}}):
            with self.subTest(sj=sj):
                shutil.rmtree(os.path.join(self.tmp, "S"), True)
                _, path = self.in_session(sj)
                self.assertTrue(sr.core_work(os.path.join(self.tmp, "S")))
                with self.assertRaises(ValueError):
                    nc.decide(path, et.RAW, ME, REASON)
                self.assertIsNone(nc.load(path)["decision"])
                # staff on the list still decide it, as Core staff
                d = nc.decide(path, et.RAW, "jdoe", REASON)["decision"]
                self.assertEqual(d["by"], staff.ROLE)

    def test_a_damaged_session_json_fails_closed(self):
        self.staff_list("jdoe")
        s, path = self.in_session()
        with open(os.path.join(s, "session.json"), "w") as fh:
            fh.write("{not json")
        self.assertTrue(sr.core_work(s))
        with self.assertRaises(ValueError):
            nc.decide(path, et.RAW, ME, REASON)

    def test_a_missing_session_json_fails_closed(self):
        """2.15.1 review minor 8: --session <a folder with no session.json> opened the own
        signer; a missing session.json is now treated like a damaged one (staff-only)."""
        self.staff_list("jdoe")
        s, path = self.in_session()
        self.assertFalse(os.path.exists(os.path.join(s, "session.json")))
        self.assertTrue(sr.core_work(s))
        with self.assertRaises(ValueError):
            nc.decide(path, et.RAW, ME, REASON)
        self.assertIsNone(nc.load(path)["decision"])

    def test_no_session_stays_staff_only(self):
        self.staff_list("jdoe")
        path = check_record(os.path.join(self.tmp, "loose", nc.RECORD))
        with self.assertRaises(ValueError):
            nc.decide(path, et.RAW, ME, REASON)

    def test_an_agent_named_account_cannot_sign_its_own(self):
        self.staff_list("jdoe")
        _, path = self.in_session()
        with mock.patch.object(staff, "_uid_login", return_value="claude"):
            with self.assertRaises(ValueError):
                nc.decide(path, et.RAW, "claude", REASON)

    def test_other_staff_only_decisions_are_unchanged(self):
        """ack-legacy and qc_bracket ack go through staff.require: still refused for a login not
        on the list, own analysis or not."""
        self.staff_list("jdoe")
        with self.assertRaises(ValueError):
            staff.require(ME)


class PtmRefusalOnBothEngines(unittest.TestCase):
    """Item 11: SKILL.md "What this skill does not do" says a PTM is never added through
    --overrides. DIA-NN's --var-mod refused phospho and K-GG (lane dia); a Sage
    database.variable_mods override with the same masses is refused the same way."""

    def test_sage_phospho_and_k_gg_are_refused(self):
        from test_var_mods import estimate
        for over in ({"database": {"variable_mods": {"S": [79.966331]}}},
                     {"database": {"variable_mods": {"M": [15.994915], "K": [114.042927]}}}):
            with self.subTest(over=over), tempfile.TemporaryDirectory() as d:
                p = estimate(d, "sage", "--overrides", json.dumps(over), acquisition="DDA",
                             ok=False)
                self.assertNotEqual(p.returncode, 0)
                out = p.stdout + p.stderr
                self.assertIn("searches and reports proteins, not modification sites", out)
                self.assertEqual(out.count("No cfg was written"), 1, out)
                self.assertFalse(os.path.exists(os.path.join(d, "sage.json")))

    def test_sage_met_ox_override_still_works(self):
        from test_var_mods import estimate
        with tempfile.TemporaryDirectory() as d:
            cfg, _ = estimate(d, "sage", "--overrides",
                              json.dumps({"database": {"variable_mods": {"M": [15.994915]}}}),
                              acquisition="DDA")
            self.assertTrue(os.path.exists(cfg))

    def test_skill_md_names_the_refusal(self):
        with open(os.path.join(os.path.dirname(HERE), "SKILL.md")) as fh:
            skill = " ".join(fh.read().split())
        self.assertIn("never add a PTM through `--overrides`", skill)


class RadiantNormalisationUnknown(unittest.TestCase):
    """Item 3: whether Radiant Fulcrum normalised across runs is not documented where the skill
    can read it, so the declaration says "unknown" -- explicitly, because a peptide-level report
    that declares nothing was read as Sage's (not normalised), with no [confirm]."""

    def setUp(self):
        import run_search
        self.rs = run_search
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def test_both_levels_declare_unknown(self):
        for col in ("PG.Normalised", "Precursor.Normalised"):
            dec = self.rs.radiant_quantity_declaration(self.d, col, "Q.Value")
            self.assertEqual(dec["engine_normalised"], self.rs.ENGINE_NORMALISED_UNKNOWN)
            self.assertIn("NOT RECORDED", dec["engine_normalisation"])
            self.assertIn("[confirm]", dec["engine_normalisation"])

    R = r'''setwd("{scripts}"); source("build_maxlfq.R")
      d <- list(level = "{level}", label = "Radiant Fulcrum", engine_normalised = "{en}",
                engine_normalisation = "E")
      b <- maxlfq_between_run("auto", "normalised", d)
      cat(b$method, "|", is.na(b$engine_normalised$value), "|", b$engine_normalised$source, "\n")'''

    def r(self, level, en):
        import subprocess
        from test_run_de_contaminants import r_has
        if not r_has("jsonlite"):
            self.skipTest("needs Rscript")
        p = subprocess.run(["Rscript", "-e", self.R.format(scripts=SCRIPTS, level=level, en=en)],
                           capture_output=True, text=True, timeout=120)
        self.assertEqual(p.returncode, 0, p.stderr[-1500:])
        return [x.strip() for x in p.stdout.split("|")]

    def test_the_between_run_rule_reads_unknown_as_not_known(self):
        for level in ("peptide", "protein"):
            method, na, source = self.r(level, "unknown")
            self.assertEqual((method, na), ("median_ratio", "TRUE"), level)
            self.assertNotIn("Sage", source)
        # a peptide-level report that declares nothing is still Sage's (pin)
        method, na, source = self.r("peptide", "")
        self.assertEqual((method, na), ("median_ratio", "FALSE"))
        self.assertIn("Sage", source)


class MethodsReadTheNewRecords(unittest.TestCase):
    """Item 4: the Methods and the report print what the records now say -- the Radiant route's
    status (make_presets route_status), the values Sage's LFQ measured by matching between runs
    (run_de.R detection_matrix) and the between-run step actually applied (the DE's `normalisation`)."""

    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)

    def radiant_search(self, with_preset=True):
        import make_presets
        out = os.path.join(self.d, "search")
        os.makedirs(out, exist_ok=True)
        cfg = os.path.join(self.d, "run.radiantConfig")
        with open(cfg, "w") as fh:
            fh.write("[General]\nthreadCount=8\n")
        if with_preset:
            with open(cfg + ".provenance.json", "w") as fh:
                json.dump({"engine": "radiant",
                           "route_status": make_presets.RADIANT_ROUTE_STATUS}, fh)
        sp = os.path.join(out, "search_provenance.json")
        with open(sp, "w") as fh:
            json.dump({"engine": "radiant", "version": "2.3.3", "params_file": cfg}, fh)
        return sp

    def test_radiant_methods_and_report_say_experimental(self):
        import make_methods
        import make_analysis_html as mah
        import make_presets
        for with_preset in (True, False):
            with self.subTest(with_preset=with_preset):
                shutil.rmtree(self.d, True)
                os.makedirs(self.d)
                sp = self.radiant_search(with_preset)
                rec = make_methods.search_record(search_prov=sp)
                self.assertEqual(rec["route_status"]["status"], "experimental")
                para = make_methods.search_paragraph(rec)
                self.assertIn("This search route is experimental", para)
                self.assertIn(make_presets.RADIANT_ROUTE_STATUS["say"], para)
                report = os.path.join(os.path.dirname(sp), "report.parquet")
                note = mah.route_status_note({"input": report})
                self.assertEqual(note["kind"], "warning")
                self.assertIn("experimental", note["text"])

    def test_other_routes_say_nothing(self):
        import make_methods
        out = os.path.join(self.d, "s")
        os.makedirs(out)
        cfg = os.path.join(self.d, "p.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n")
        sp = os.path.join(out, "search_provenance.json")
        with open(sp, "w") as fh:
            json.dump({"engine": "diann", "version": "2.7.0", "params_file": cfg}, fh)
        rec = make_methods.search_record(search_prov=sp)
        self.assertIsNone(rec["route_status"])
        self.assertEqual(make_methods.route_status_sentence(rec), "")

    def test_matched_values_are_counted_in_the_methods(self):
        import make_methods
        import sage_lfq_check as slc
        prov = {"display_label": "Peptide median polish + limma",
                "normalisation": "median-ratio centring of the engine's quantities",
                "detection_matrix": {"matched_cells": 1234, "quantified_cells": 20000,
                                     "matched_definition": slc.MATCHED_BETWEEN_RUNS}}
        para = make_methods.de_paragraph(prov)
        self.assertIn("Of the 20,000 quantified protein × run values, 1,234 had no peptide of "
                      "the protein identified by MS2 in that run", para)
        self.assertIn(slc.MATCHED_BETWEEN_RUNS, para)          # quoted, not re-worded
        self.assertIn("runs with no MS1 signal were missing values", para)
        self.assertIn("(matched_between_runs_only_in)", para)
        self.assertNotIn("transferred", para)
        self.assertIn("Between-run normalisation: median-ratio centring", para)
        prov["detection_matrix"] = {"zero_means": "missing"}
        self.assertNotIn("matched_between_runs_only_in", make_methods.de_paragraph(prov))


class ContainerTagIsNotAVersion(unittest.TestCase):
    """Item 5: session.search_engine (the report header, the Slack post, the README) prints a
    Radiant container tag the way the Methods do, never as "Radiant 2.3.3"."""

    def test_one_rule(self):
        import resolve_defaults as rd
        self.assertEqual(rd.engine_version_text("diann", "2.7.0"), "2.7.0")
        self.assertIsNone(rd.engine_version_text("sage", None))
        self.assertEqual(rd.engine_version_text("radiant", "2.3.3", label="Radiant"),
                         "(radiant-fulcrum container 2.3.3; the Radiant version inside it was "
                         "not recorded)")
        self.assertEqual(rd.engine_version_text("radiant", "2.3.3", "1.4.2"), "1.4.2")
        self.assertEqual(rd.engine_version_text("radiant", "2.3.3", pinned=True),
                         "radiant-fulcrum container 2.3.3")

    def test_session_search_engine(self):
        import session
        import make_methods
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        p = session.paths_for(d)
        se = session.search_engine(p, record={"engine": "radiant", "version": "2.3.3"})
        self.assertIn("radiant-fulcrum container 2.3.3", se["version"])
        self.assertEqual(session.search_engine(p, record={"engine": "diann",
                                                          "version": "2.7.0"})["version"], "2.7.0")
        # the Methods say the same
        out = os.path.join(d, "s")
        os.makedirs(out)
        sp = os.path.join(out, "search_provenance.json")
        with open(sp, "w") as fh:
            json.dump({"engine": "radiant", "version": "2.3.3"}, fh)
        self.assertEqual(make_methods.search_record(search_prov=sp)["version"], se["version"])
        # a manifest pin, no record
        os.makedirs(os.path.dirname(p["workflow_manifest"]), exist_ok=True)
        with open(p["workflow_manifest"], "w") as fh:
            json.dump({"engine": {"name": "radiant", "version": "2.3.3"}}, fh)
        self.assertEqual(session.search_engine(p)["version"], "radiant-fulcrum container 2.3.3")


class InitRecordsTheAgent(unittest.TestCase):
    """Item 5: session.py init records `ai_agent` when the environment names the assistant, so a
    report written later in a job (no CLAUDECODE there) still names it; never a guess."""

    def init(self, env):
        import subprocess
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        e = {k: v for k, v in os.environ.items()
             if not k.startswith("CLAUDE") and k != "ANTHROPIC_MODEL"}
        e.update(env)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "session.py"), "init",
                            "--name", "t", "--base", d], capture_output=True, text=True,
                           timeout=120, env=e)
        self.assertEqual(r.returncode, 0, r.stderr)
        s = json.loads(r.stdout)["created"]
        try:
            with open(os.path.join(s, "session.json")) as fh:
                return s, json.load(fh)
        except FileNotFoundError:
            return s, {}

    def test_claude_code_is_recorded(self):
        import ai_agent
        s, sj = self.init({"CLAUDECODE": "1", "ANTHROPIC_MODEL": "m-9"})
        self.assertEqual(sj["ai_agent"],
                         "Claude Code (Anthropic's AI coding assistant, model m-9)")
        self.assertEqual(ai_agent.agent_label(s, env={}), sj["ai_agent"])

    def test_nothing_is_guessed(self):
        _, sj = self.init({})
        self.assertNotIn("ai_agent", sj)


class RadiantLogsAreAbsolute(unittest.TestCase):
    """Item 8: radiant_parallel.py's #SBATCH --output is absolute (the step scripts' folder), as
    diann_parallel.header's, so a step resubmitted from $HOME logs where the watcher looks."""

    def test_header(self):
        import radiant_parallel as rp
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        h = rp.header("r2_search", 4, 8, 2, "high", "acct", None, log_dir=d)
        out = [ln for ln in h.splitlines() if ln.startswith("#SBATCH --output")]
        self.assertEqual(out, [f'#SBATCH --output="{os.path.abspath(d)}/%x_%j.log"'])


class EngineChoiceAndUnknownInstrument(unittest.TestCase):
    """Item 11: Radiant's DE input is the adapted protein table (maxlfq), as the protein screen
    and SKILL.md say; an unrecognised instrument needs confirmation, with a plain note."""

    def test_radiant_de_method(self):
        import engine_choice as ec
        for acq in ("DIA", "DDA"):
            self.assertEqual(ec.de_method("radiant", acq), "maxlfq")
        self.assertEqual(ec.de_method("diann", "DIA"), "dpc")

    def test_unknown_instrument_needs_confirmation(self):
        import subprocess
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "resolve_defaults.py"),
                            "--acquisition", "DIA", "--instrument", "Mystery 9000", "--dest",
                            os.path.join(d, "wf")], capture_output=True, text=True, timeout=120)
        self.assertEqual(r.returncode, 0, r.stderr)
        out = json.loads(r.stdout)
        self.assertTrue(out["needs_confirmation"])
        self.assertIn("not recognised", out["ask_user"])
        self.assertIn("preliminary analyses only", out["ask_user"])
        self.assertTrue(any("not recognised" in n for n in out["notes"]))
        with open(out["manifest"]) as fh:
            self.assertTrue(json.load(fh)["needs_confirmation"])

    def test_a_known_instrument_is_unchanged(self):
        import resolve_defaults as rd
        self.assertIsNone(rd.unknown_instrument_note("timstof", "diann", "timsTOF HT"))


class DeliveryNamesTheTableOfRecord(unittest.TestCase):
    """Item 12: where a delivery README describes its files, it says the DE tables are the
    results of record and the engine's protein matrix is kept for reference."""

    def test_tables_description(self):
        import core_submission as cs
        t = cs.tables_description({"tables/DE_dpc_B-A.csv"})
        self.assertIn("results of record", t)
        self.assertNotIn("results of record", cs.tables_description({"tables/methods.txt"}))


class CitationsAsVerified(unittest.TestCase):
    """Items 6 and 7: the citations checked against their sources on 2026-11-17."""

    def read(self, *parts):
        with open(os.path.join(os.path.dirname(HERE), *parts), encoding="utf-8") as fh:
            return fh.read()

    def test_dpc_cites_limma(self):
        self.assertIn("limma: Ritchie ME et al. (2015) Nucleic Acids Res 43(7):e47",
                      self.read("scripts", "run_de.R"))

    def test_fry_is_slides_not_peer_reviewed(self):
        for parts in (("scripts", "run_sets.R"), ("scripts", "set_tests.R"),
                      ("references", "set-tests.md")):
            with self.subTest(parts=parts):
                t = self.read(*parts)
                self.assertIn("doi:10.7490/f1000research.1113351.1", t)
                self.assertIn("not peer", t)
                self.assertNotRegex(t.lower(), r"\bposter\b")     # "posterior.sd" is fine

    def test_hpa_citation(self):
        t = self.read("references", "protein-card-sources.md")
        self.assertIn("Science 347(6220):1260419", t)
        self.assertIn("doi:10.1126/science.1260419", t)


if __name__ == "__main__":
    unittest.main()
