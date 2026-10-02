#!/usr/bin/env python3
"""
Normalised or non-normalised quantities: the experiment type's default, and the data check
that catches it when the type is wrong (2.10, item 9).

Global normalisation assumes most proteins do not change. An IP's controls carry little
protein, so normalising them up hides the enrichment; a TurboID's controls carry plenty, so the
standard holds -- unless they came out low. experiment_type.py decides the type once (from the
submission, confirmed, or asked) with a default for the DE's quantities; normalization_check.py
checks the data for EVERY type -- normalisation factors against the groups, the shift before
normalisation, identifications per group, Normalisation.Instability, and each contrast's volcano
both ways -- and stops to ask when the data disagree. run_de.R --quantities reads DIA-NN's
non-normalised Precursor.Quantity (or a --no-norm report on maxlfq, with no quantile step), and
--normalization-check holds the final DE to the decision and records it.

Synthetic DIA-NN-shaped reports (normalised = Normalisation.Factor x non-normalised, each run
scaled to a common median, as a cross-run normalisation does):
  whole      two conditions, 20 proteins up and 20 down -- no trip
  ip         IP vs IgG, the IgG carrying 4x less protein and no interactors -- trips when the
             type says whole proteome; agrees when it says IP (default non-normalised)
  turbo      TurboID bait vs control with plenty of signal, bait + carboxylases -- no trip,
             carboxylase stability reported both ways
  turbo_low  the same with the controls 4x low -- trips although the type is proximity
Pure-Python parts run everywhere; the report readers need pyarrow; the end-to-end flow needs R.
"""
import csv
import json
import math
import os
import random
import shutil
import subprocess
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
from test_run_de_contaminants import r_has  # noqa: E402

CHECK = os.path.join(SCRIPTS, "normalization_check.py")
ETYPE = os.path.join(SCRIPTS, "experiment_type.py")
RUN_DE = os.path.join(SCRIPTS, "run_de.R")
NEEDS = ("limpa", "limma", "arrow", "dplyr", "tidyr", "jsonlite")

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False


def make_report(out, case, seed=7, n_bg=120, runs_per=3, no_norm=False):
    """A DIA-NN-shaped report.parquet + conditions.csv for `case` -> (report, conditions, contrast)."""
    rnd = random.Random(seed)
    n_int = 20 if case.startswith("turbo") else 80
    os.makedirs(out, exist_ok=True)
    g1, g2 = {"whole": ("B", "A"), "ip": ("IP", "IgG"), "turbo": ("Bait", "Ctrl"),
              "turbo_low": ("Bait", "Ctrl")}[case]
    runs = [(f"{g1}_{i}", g1) for i in range(1, runs_per + 1)] + \
           [(f"{g2}_{i}", g2) for i in range(1, runs_per + 1)]
    load = {r: rnd.gauss(0, 0.15) for r, _ in runs}
    prots = [(f"P{i:05d}", f"G{i}", rnd.uniform(15, 22), "bg") for i in range(n_bg)]
    prots += [(f"Q{i:05d}", f"I{i}", rnd.uniform(15, 22), "int") for i in range(n_int)]
    if case != "whole":
        prots.append(("B00001", "BAITX", 21.0, "bait"))
    if case.startswith("turbo"):
        prots += [(f"C{j:05d}", g, 20.0 + 0.3 * j, "carbox") for j, g in enumerate(et.CARBOXYLASES)]

    def level(kind, grp, base, i):
        if case == "whole":
            return base + (0 if kind != "int" or i >= 40 or grp != g1 else
                           1.5 if i % 2 == 0 else -1.5)
        if case == "ip":
            if grp == g2:          # IgG: 4x less protein, no interactors, no bait
                return None if kind in ("int", "bait") else base - 2.0
            return base + {"int": 3.0, "bait": 6.0}.get(kind, 0.0)
        if grp == g2:              # TurboID control: plenty of signal (or 4x low), no bait
            return None if kind == "bait" else base - (2.0 if case == "turbo_low" else 0.0)
        return base + {"int": 3.0, "bait": 6.0}.get(kind, 0.0)

    cols = ("Run", "Precursor.Id", "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
            "Proteotypic", "Precursor.Quantity", "Q.Value", "Lib.Q.Value", "Lib.PG.Q.Value",
            "PG.Q.Value", "Global.Q.Value", "Global.PG.Q.Value")
    rows = {k: [] for k in cols}
    for pi, (acc, gene, base, kind) in enumerate(prots):
        offs = [rnd.gauss(0, 0.8) for _ in range(3)]
        for r, grp in runs:
            lv = level(kind, grp, base, pi - n_bg)
            if lv is None:
                continue
            lv += rnd.gauss(0, 0.3)            # the protein's biological variation
            for k, o in enumerate(offs):
                v = lv + o + load[r] + rnd.gauss(0, 0.15)
                if rnd.random() < 1 / (1 + math.exp((v - 15.5) * 2)):
                    continue                     # low signal is missed more often
                for c, x in (("Run", r), ("Precursor.Id", f"{gene}PEP{k}K2"),
                             ("Protein.Group", acc), ("Protein.Ids", acc),
                             ("Protein.Names", gene + "_HUMAN"), ("Genes", gene),
                             ("Proteotypic", 1), ("Precursor.Quantity", 2 ** v)):
                    rows[c].append(x)
                for q in cols[8:]:
                    rows[q].append(0.001)
    by_run = {}
    for r, q in zip(rows["Run"], rows["Precursor.Quantity"]):
        by_run.setdefault(r, []).append(math.log2(q))
    meds = {r: sorted(v)[len(v) // 2] for r, v in by_run.items()}
    ref = sorted(meds.values())[len(meds) // 2]
    fac = {r: 1.0 if no_norm else 2 ** (ref - m) for r, m in meds.items()}
    rows["Normalisation.Factor"] = [fac[r] for r in rows["Run"]]
    rows["Precursor.Normalised"] = [q * fac[r] for r, q in zip(rows["Run"],
                                                              rows["Precursor.Quantity"])]
    pgq = {}
    for p, r, v in zip(rows["Protein.Group"], rows["Run"], rows["Precursor.Normalised"]):
        pgq[(p, r)] = pgq.get((p, r), 0) + v
    rows["PG.MaxLFQ"] = [pgq[(p, r)] for p, r in zip(rows["Protein.Group"], rows["Run"])]
    report = os.path.join(out, "report.parquet")
    pq.write_table(pa.table(rows), report)
    cond = os.path.join(out, "conditions.csv")
    with open(cond, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["File.Name", "Group"])
        w.writerows(runs)
    return report, cond, f"{g1}-{g2}"


def session(d, etype=None, added=(), **kw):
    s = os.path.join(d, "session")
    os.makedirs(os.path.join(s, "input"), exist_ok=True)
    if etype:
        et.record(s, etype, "user", **kw)
    if added:
        add_sequences(s, added)
    return s


def add_sequences(s, accessions):
    """The session's FASTA record as fetch_fasta.py --add-fasta writes it (added_sequences)."""
    entries = [{"accession": a, "name": f"IN|{a}|ADDED", "header": f"IN|{a}|ADDED test",
                "length": 239, "sequence_sha256": "0" * 64} for a in accessions]
    with open(os.path.join(s, et.FASTA_META), "w", encoding="utf-8") as fh:
        json.dump({"added_sequences": [{"file": "/x/bait.fasta", "sha256": "1" * 64,
                                        "n_entries": len(entries), "entries": entries}]}, fh)


# ------------------------------------------------------------------- the experiment type --
class ExperimentType(unittest.TestCase):
    def test_an_added_sequence_is_the_bait_candidate(self):
        """Integration 2.10: --add-fasta's sequences (a bait such as EGFP) are the bait for the
        data check when none was given or recorded -- one of them, never a guess between two."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s = session(os.path.join(d, "none"))
        self.assertEqual(et.bait_for_check(s), (None, "not given", []))
        s = session(os.path.join(d, "one"), added=("B00001",))
        bait, src, cands = et.bait_for_check(s)
        self.assertEqual((bait, len(cands)), ("B00001", 1))
        self.assertIn("--add-fasta", src)
        self.assertIn("not confirmed", src)
        self.assertEqual(et.bait_for_check(s, "BAITX")[:2], ("BAITX", "--bait"))
        self.assertEqual(et.bait_for_check(s, rec={"bait": "GFP"})[:2],
                         ("GFP", "the experiment-type record"))
        s = session(os.path.join(d, "two"), added=("B00001", "B00002"))
        bait, src, cands = et.bait_for_check(s)
        self.assertIsNone(bait)
        self.assertEqual([c["accession"] for c in cands], ["B00001", "B00002"])
        self.assertIn("several", src)
        with open(os.path.join(s, et.FASTA_META), "w") as fh:
            fh.write("{not json")
        self.assertEqual(et.bait_candidates(s), [])

    def test_propose_lists_the_added_sequences_as_bait_candidates(self):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s = session(d, added=("B00001",))
        p = subprocess.run([sys.executable, ETYPE, "propose", "--session", s],
                           capture_output=True, text=True, timeout=60)
        self.assertEqual(p.returncode, 0, p.stderr)
        out = json.loads(p.stdout)
        self.assertEqual([c["accession"] for c in out["bait_candidates"]], ["B00001"])
        self.assertIn("set-bait --bait <accession> (or --bait none)", out["ask_bait"])

    def test_the_bait_answer_is_recorded_and_never_left_unconfirmed(self):
        """2.10 safety review: a single --add-fasta sequence was used as the bait unconfirmed.
        set-bait records the user's answer to ask_bait -- an accession, or none."""
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        s = session(d, added=("B00001", "B00002"))

        def etype(*args):
            return subprocess.run([sys.executable, ETYPE, *args, "--session", s],
                                  capture_output=True, text=True, timeout=60)
        p = etype("set-bait", "--bait", "B00001")
        self.assertNotEqual(p.returncode, 0)
        self.assertIn("record it first", p.stderr)
        self.assertEqual(etype("set", "--type", "ip", "--source", "user").returncode, 0)
        p = etype("set-bait", "--bait", "B00009")
        self.assertIn("not one of the added sequences (B00001, B00002)", p.stderr)
        p = etype("set-bait", "--bait", "B00002")
        self.assertEqual(p.returncode, 0, p.stderr)
        r = et.load(s)
        self.assertEqual((r["bait"], r["bait_confirmed"]["answer"]), ("B00002", "B00002"))
        self.assertIn("asked ask_bait", r["bait_confirmed"]["by"])
        self.assertEqual(et.bait_for_check(s, rec=r)[:2], ("B00002", "the experiment-type record"))
        # re-recording the type keeps the answer
        self.assertEqual(etype("set", "--type", "ip", "--source", "user").returncode, 0)
        self.assertEqual(et.load(s)["bait"], "B00002")
        # "none": no added sequence is the bait -- the one-sequence fallback never applies
        s1 = session(os.path.join(d, "one"), added=("B00001",))
        et.record(s1, "ip", "user")
        et.record_bait(s1, "none")
        bait, src, _c = et.bait_for_check(s1, rec=et.load(s1))
        self.assertEqual((bait, src), (None, "the user said none of the added sequences is the bait"))

    def test_source_submission_needs_an_attached_submission(self):
        with tempfile.TemporaryDirectory() as d:
            s = session(d)
            p = subprocess.run([sys.executable, ETYPE, "set", "--session", s, "--type",
                                "whole_proteome", "--source", "submission"],
                               capture_output=True, text=True, timeout=60)
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("no submission is attached to this session", p.stderr)

    def test_the_submission_proposes_and_the_most_specific_wins(self):
        self.assertEqual(et.propose_from_submission(
            ["DIA (Quantitative)", "Affinity Purification (magnetic Beads)"])[0], "ip")
        self.assertEqual(et.propose_from_submission(
            ["Affinity Purification (Streptavidin)", "TurboID proximity labelling"])[0],
            "proximity")
        self.assertEqual(et.propose_from_submission(["Global proteomics"])[0], "whole_proteome")
        self.assertEqual(et.propose_from_submission(["Secretome"])[0], "secretome")
        self.assertEqual(et.propose_from_submission(["DIA (Quantitative)"]), (None, None))
        # fractions: compared (a separation) vs combined (depth); a bare "fractionation" asks
        self.assertEqual(et.propose_from_submission(["SEC fractions"])[0], "fractions_separation")
        self.assertEqual(et.propose_from_submission(["Sucrose gradient"])[0],
                         "fractions_separation")
        self.assertEqual(et.propose_from_submission(["BN-PAGE complexome"])[0],
                         "fractions_separation")
        self.assertEqual(et.propose_from_submission(["Organelle fractionation"])[0],
                         "fractions_separation")
        self.assertEqual(et.propose_from_submission(["High-pH fractionation"])[0],
                         "fractions_depth")
        self.assertEqual(et.propose_from_submission(["Fractionation"]), (None, None))
        self.assertEqual(et.propose_from_submission(["Global proteome profiling"])[0],
                         "whole_proteome")

    def test_each_type_has_a_default_and_unknown_is_tagged(self):
        for t in ("ip", "fractions_separation"):
            self.assertEqual(et.default_for(t)["quantities"], et.RAW, t)
        for t in ("whole_proteome", "proximity", "secretome", "fractions_depth", "other"):
            self.assertEqual(et.default_for(t)["quantities"], et.NORMALISED, t)
        self.assertIn("DEFAULT -- not user-confirmed", et.default_for(None)["why"])
        # the separation default cites DIA-NN's README
        self.assertIn("SEC, \"normalisation should not be used\"",
                      et.default_for("fractions_separation")["why"])
        self.assertNotIn("fractions", et.DEFAULTS, "a bare 'fractions' class says neither")

    def test_the_submitters_loading_is_reconciled(self):
        self.assertIn("by volume", et.reconcile("whole_proteome", "normalized by volume"))
        self.assertIn("protein amount", et.reconcile("ip", "equal protein mass (ug)"))
        self.assertIsNone(et.reconcile("proximity", "by protein amount"))
        self.assertIsNone(et.reconcile("ip", "by volume"))
        self.assertIn("not known", et.reconcile("ip", "lets discuss/no idea"))
        self.assertIn("each fraction", et.reconcile("fractions_separation", "10 ug per fraction"))

    def test_record_and_the_pulldown_design_come_from_one_place(self):
        with tempfile.TemporaryDirectory() as d:
            s = session(d, "proximity", bait="BAITX", controls="Ctrl")
            rec = et.load(s)
            self.assertEqual((rec["type"], rec["bait"], rec["controls"]),
                             ("proximity", "BAITX", ["Ctrl"]))
            pd = et.pulldown_design(rec, ["Bait", "Ctrl"], ["Bait-Ctrl"])
            self.assertEqual((pd["pulldown"], pd["vs_control"], pd["note"]),
                             (True, ["Bait-Ctrl"], None))
            # the type says whole proteome but the groups read as IP controls: raised
            et.record(s, "whole_proteome", "user")
            pd = et.pulldown_design(et.load(s), ["IP", "IgG"], ["IP-IgG"])
            self.assertFalse(pd["pulldown"])
            self.assertIn("confirm the type", pd["note"])
            # an IP type with no control group named: raised too
            et.record(s, "ip", "user")
            self.assertIn("record the control groups",
                          et.pulldown_design(et.load(s), ["A", "B"], ["A-B"])["note"])
        # no record: the old rule, said so
        pd = et.pulldown_design(None, ["IP", "IgG"], ["IP-IgG"])
        self.assertTrue(pd["pulldown"])
        self.assertIn("no experiment type recorded", pd["source"])

    def test_a_submission_whose_loading_disagrees_is_asked(self):
        """The CoreOmics record proposes the type; its Normalization answer (how the samples were
        loaded) contradicting the default is raised at `propose`, at `set`, and in the check."""
        import submission_report as sr
        from test_submission_report import fixture
        with tempfile.TemporaryDirectory() as d:
            s = session(d)
            rec = fixture()
            rec["submission_data"].update(proteomics_type=["Global proteomics"],
                                          volume_or_mass="normalized by volume")
            sr.attach(s, sr.sanitize(rec))
            p = subprocess.run([sys.executable, ETYPE, "propose", "--session", s],
                               capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            j = json.loads(p.stdout)
            self.assertEqual(j["proposed"], "whole_proteome")
            self.assertIn("Confirm with the user", j["ask"])
            self.assertIn("by volume", j["reconcile"])
            p = subprocess.run([sys.executable, ETYPE, "set", "--session", s, "--type",
                                "whole_proteome", "--source", "submission"],
                               capture_output=True, text=True, timeout=60)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertIn("[experiment_type] ASK:", p.stderr)
            r = et.load(s)
            self.assertEqual((r["source"], r["submission_normalisation"]),
                             ("the CoreOmics submission, confirmed by the user", "normalized by volume"))
            md = nc.summary_md({"experiment_type": {"label": r["label"], "source": r["source"],
                                                    "reconcile": r["reconcile"]},
                                "default": r["default"]})
            self.assertIn("Ask the user (the submission):", md)

    def test_the_cli_refuses_an_unknown_type(self):
        with tempfile.TemporaryDirectory() as d:
            s = session(d)
            p = subprocess.run([sys.executable, ETYPE, "set", "--session", s, "--type", "guess",
                                "--source", "user"], capture_output=True, text=True, timeout=60)
            self.assertNotEqual(p.returncode, 0)


# --------------------------------------------------------------------------- the volcano --
def de_rows(lfc_p):
    return [{"Protein.Group": f"P{i}", "Genes": g, "logFC": str(l), "P.Value": str(p),
             "adj.P.Val": str(min(1.0, p * 10))} for i, (g, l, p) in enumerate(lfc_p)]


class Volcano(unittest.TestCase):
    def test_a_balanced_volcano_has_no_anomaly(self):
        rnd = random.Random(1)
        rows = de_rows([(f"G{i}", rnd.gauss(0, 0.2), rnd.uniform(0.1, 1)) for i in range(300)]
                       + [(f"U{i}", 2.0, 1e-5) for i in range(15)]
                       + [(f"D{i}", -2.0, 1e-5) for i in range(15)])
        v = nc.volcano(rows, 0.05, False, None)
        self.assertEqual(v["anomalies"], [])
        self.assertEqual((v["up"], v["down"]), (15, 15))
        self.assertLess(abs(v["centre_nonsig"]), 0.1)

    def test_depleted_background_in_an_enrichment_contrast(self):
        rows = de_rows([(f"G{i}", -0.8, 1e-4) for i in range(60)]
                       + [(f"G{i}", -0.7, 0.5) for i in range(60, 200)]
                       + [("BAITX", 0.3, 0.2)] + [(f"I{i}", 3.0, 1e-6) for i in range(10)])
        v = nc.volcano(rows, 0.05, True, "BAITX")
        a = " | ".join(v["anomalies"])
        self.assertIn("significantly DOWN", a)
        self.assertIn("centre shift -0.70", a)
        self.assertIn("the bait BAITX is not among the 10 most enriched", a)
        self.assertEqual(v["bait"]["rank"], 11)

    def test_too_much_significant_lopsided_and_an_odd_p_histogram(self):
        rows = de_rows([(f"G{i}", 1.0, 1e-6) for i in range(150)]
                       + [(f"H{i}", 0.0, 0.99) for i in range(100)])
        v = nc.volcano(rows, 0.05, False, None)
        a = " | ".join(v["anomalies"])
        self.assertIn("60% of proteins significant", a)
        self.assertIn("lopsided: 150 up vs 0 down", a)
        self.assertIn("excess near 1", a)

    def test_the_raw_version_of_an_enrichment_contrast_expects_background_up(self):
        de = {"prov": {"adjp": 0.05}, "tables": {"IP-IgG": de_rows(
            [(f"G{i}", 2.0, 1e-6) for i in range(150)] + [(f"H{i}", 1.0, 0.3) for i in range(60)])}}
        design = {"vs_control": ["IP-IgG"]}
        self.assertEqual(nc.volcano_version(de, design, None, raw=True)["contrasts"]["IP-IgG"]
                         ["anomalies"], [])
        self.assertTrue(nc.volcano_version(de, design, None, raw=False)["contrasts"]["IP-IgG"]
                        ["anomalies"])


# ----------------------------------------------------------------- the verdict and gate --
def no_staff_list(test, d):
    """Hermetic `decide --by`: no Core staff list (staff.py reads CORE_STAFF_FILE), so a person's
    name passes, said to be unchecked -- tests/test_staff_decisions.py tests the list itself."""
    p = mock.patch.dict(os.environ, {"CORE_STAFF_FILE": os.path.join(d, "no_core_staff.txt")})
    p.start()
    test.addCleanup(p.stop)


def record(default, **checks):
    c = {"factor_confound": {}, "raw_shift": {}, "id_gap": {}, "instability": {}}
    c.update(checks)
    return {"default": {"quantities": default}, "checks": c, "volcano": {}}


class VerdictAndGate(unittest.TestCase):
    def test_what_trips_depends_on_the_default(self):
        big = {"large_and_confounded": True, "gap_fold": 4.0, "highest": "IgG", "lowest": "IP",
               "run_spread_log2": 2.2, "large_unexplained": False}
        self.assertTrue(nc.trips(record(et.NORMALISED, factor_confound=big)))
        self.assertEqual(nc.trips(record(et.RAW, factor_confound=big)), [],
                         "an IP whose controls are scaled up agrees with non-normalised")
        loose = dict(big, large_and_confounded=False, large_unexplained=True)
        self.assertTrue(nc.trips(record(et.RAW, factor_confound=loose)))
        self.assertTrue(nc.trips(record(et.NORMALISED, id_gap={"gap": True, "ratio": 2.0})))

    def test_no_silent_switch_either_way(self):
        with tempfile.TemporaryDirectory() as d:
            no_staff_list(self, d)
            path = os.path.join(d, nc.RECORD)
            rec = dict(record(et.NORMALISED), schema=nc.SCHEMA, schema_version=1, tripped=True,
                       trips=["B1 ..."], decision=None)
            with open(path, "w") as fh:
                json.dump(rec, fh)
            ok, msg, _ = nc.gate(path, et.NORMALISED)
            self.assertFalse(ok)
            self.assertIn("none is recorded", msg)
            with self.assertRaises(ValueError):
                nc.decide(path, et.RAW, "", "why")             # who chose is required
            nc.decide(path, et.RAW, "canalyst", "IgG controls carry little protein")
            self.assertFalse(nc.gate(path, et.NORMALISED)[0])
            ok, msg, rec = nc.gate(path, et.RAW)
            self.assertTrue(ok)
            # the record names the role; the person is in the staff-only record beside it
            self.assertEqual(rec["decision"]["by"], "Core staff")
            self.assertEqual(rec["decision"]["default_was"], et.NORMALISED)
            with open(os.path.join(d, nc.STAFF_RECORD)) as fh:
                self.assertEqual(json.load(fh)["entries"][0]["by"], "canalyst")


class MethodsSayTheTypeAndWhy(unittest.TestCase):
    def test_separation_fractions_are_non_normalised_with_the_readme_reason(self):
        import make_methods
        d = et.default_for("fractions_separation")
        prov = {"normalisation": "none: DIA-NN's non-normalised Precursor.Quantity",
                "normalization_check": {
                    "status": "decided", "tripped": False,
                    "experiment_type": {"type": "fractions_separation",
                                        "label": et.DEFAULTS["fractions_separation"][0],
                                        "source": "the user"},
                    "default": {"quantities": d["quantities"], "why": d["why"]},
                    "decision": {"quantities": d["quantities"], "by": "the skill"}}}
        s = make_methods.de_normalisation_sentence(prov)
        self.assertIn("Experiment type: separation / profiling fractions compared with each "
                      "other (SEC, density or sucrose gradient", s)
        self.assertIn("the default for it is non-normalised quantities (each fraction holds a "
                      "different part of the proteome", s)
        self.assertIn('"normalisation should not be used"', s)
        self.assertIn("agreed with it", s)


class ReproduceReplaysTheQuantities(unittest.TestCase):
    def test_reproduce_sh_names_the_quantities_the_de_read(self):
        with tempfile.TemporaryDirectory() as d:
            de = os.path.join(d, "tables")
            os.makedirs(de)
            with open(os.path.join(de, "de_provenance.json"), "w") as fh:
                json.dump({"method": "dpc", "quantities": "raw"}, fh)
            out = os.path.join(d, "repro")
            p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"), "--outdir",
                                out, "--de-dir", de, "--de-method", "dpc",
                                "--timestamp", "2026-10-01T00:00:00Z"],
                               capture_output=True, text=True, timeout=300)
            self.assertEqual(p.returncode, 0, p.stderr[-1500:])
            with open(os.path.join(out, "reproduce.sh")) as fh:
                self.assertIn("--quantities raw", fh.read())


# ---------------------------------------------------------- the report checks (pyarrow) --
@unittest.skipUnless(HAVE_ARROW, "pyarrow is needed to read the reports")
class ReportChecks(unittest.TestCase):
    def checks(self, case):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        report, cond, contrast = make_report(d, case)
        groups = nc.read_conditions(cond)
        data, why = nc.read_report(report, set(groups))
        self.assertIsNone(why)
        return (nc.factor_confound(data, groups), nc.raw_shift(data, groups, [contrast]),
                nc.id_gap(data, groups))

    def test_whole_proteome_factors_are_small(self):
        fc, rs, ig = self.checks("whole")
        self.assertFalse(fc["large_and_confounded"])
        self.assertLess(fc["gap_log2"], nc.FACTOR_GAP_LOG2)
        self.assertFalse(rs["per_contrast"]["B-A"]["lopsided"])
        self.assertFalse(ig["gap"])

    def test_an_ips_scaled_up_controls_are_the_signature(self):
        fc, rs, ig = self.checks("ip")
        self.assertTrue(fc["large_and_confounded"])
        self.assertEqual((fc["highest"], fc["lowest"]), ("IgG", "IP"))
        self.assertGreater(fc["gap_fold"], 3)
        self.assertTrue(rs["per_contrast"]["IP-IgG"]["lopsided"])
        self.assertTrue(ig["gap"])
        self.assertEqual(fc["permutation"], "exact")

    def test_an_adapted_report_is_not_assessed_and_says_why(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "report.parquet")
            pq.write_table(pa.table({"Run": ["a"], "Protein.Group": ["P1"], "PG.MaxLFQ": [1.0]}), p)
            data, why = nc.read_report(p, {"a"})
        self.assertIsNone(data)
        self.assertIn("not a DIA-NN precursor report", why)


# ------------------------------------------------------------------- end to end (R) --
@unittest.skipUnless(HAVE_ARROW and r_has(*NEEDS), "needs pyarrow and R with " + "/".join(NEEDS))
class EndToEnd(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        no_staff_list(self, self.d)

    def de(self, report, cond, contrast, quantities, out, *extra):
        return subprocess.run(["Rscript", RUN_DE, "--input", report, "--metadata", cond,
                               "--contrasts", contrast, "--quantities", quantities,
                               "--outdir", out, *extra], capture_output=True, text=True,
                              timeout=900, cwd=self.d)

    def check(self, case, etype=None, added=(), **kw):
        report, cond, contrast = make_report(os.path.join(self.d, case), case)
        s = session(os.path.join(self.d, case), etype, added=added, **kw)
        dirs = {}
        for q in (et.NORMALISED, et.RAW):
            dirs[q] = os.path.join(self.d, case, "norm_check", q)
            p = self.de(report, cond, contrast, q, dirs[q])
            self.assertEqual(p.returncode, 0, p.stderr[-2000:])
        p = subprocess.run([sys.executable, CHECK, "run", "--report", report, "--conditions",
                            cond, "--de-normalised", dirs[et.NORMALISED], "--de-raw",
                            dirs[et.RAW], "--session", s], capture_output=True, text=True,
                           timeout=300)
        out = os.path.join(self.d, case, "norm_check")
        with open(os.path.join(out, nc.RECORD)) as fh:
            rec = json.load(fh)
        with open(os.path.join(out, nc.SUMMARY)) as fh:
            md = fh.read()
        return p, rec, md, (report, cond, contrast), os.path.join(out, nc.RECORD)

    def test_whole_proteome_passes_and_the_final_de_records_it(self):
        p, rec, md, (report, cond, contrast), path = self.check("whole", "whole_proteome")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertEqual(rec["trips"], [])
        self.assertEqual(rec["decision"]["quantities"], et.NORMALISED)
        v = rec["volcano"]["normalised"]["contrasts"]["B-A"]
        for k in ("centre_nonsig", "up", "down", "frac_sig", "pvalue"):
            self.assertIn(k, v)                            # reported even when nothing trips
        self.assertIn("## Side by side: normalised vs non-normalised", md)
        final = os.path.join(self.d, "final")
        q = self.de(report, cond, contrast, et.NORMALISED, final, "--normalization-check", path)
        self.assertEqual(q.returncode, 0, q.stderr[-2000:])
        with open(os.path.join(final, "de_provenance.json")) as fh:
            prov = json.load(fh)
        n = prov["normalization_check"]
        self.assertEqual((n["schema"], n["schema_version"], n["status"], n["quantities_applied"]),
                         (nc.SCHEMA, 1, "decided", et.NORMALISED))
        for k in ("experiment_type", "default", "checks", "volcano", "trips", "decision"):
            self.assertIn(k, n)
        import make_methods
        sentence = make_methods.de_normalisation_sentence(prov)
        self.assertIn("Experiment type: whole-proteome (expression) comparison (per the user)", sentence)
        self.assertIn("agreed with it", sentence)

    def test_user_says_whole_proteome_but_the_data_say_enrichment(self):
        p, rec, md, (report, cond, contrast), path = self.check("ip", "whole_proteome")
        self.assertEqual(p.returncode, nc.TRIPPED_EXIT, p.stdout + p.stderr)
        self.assertIn("TRIPPED", p.stderr)
        self.assertTrue(any(t.startswith("B1 ") for t in rec["trips"]), rec["trips"])
        self.assertIsNone(rec["decision"])
        # the type's default stays the recommendation (2.10 review MED 6), and the trip is the
        # question to ask: is the type wrong?
        self.assertTrue(rec["recommendation"].startswith("normalised quantities"),
                        rec["recommendation"])
        self.assertIn("the experiment type is wrong and raw quantities fit", rec["recommendation"])
        self.assertIn("**Decision: not made.**", md)
        self.assertIn("confirm the type", rec["design_note"])
        # the final DE refuses until someone chooses -- in either direction
        final = os.path.join(self.d, "final")
        for q in (et.NORMALISED, et.RAW):
            r = self.de(report, cond, contrast, q, final, "--normalization-check", path)
            self.assertNotEqual(r.returncode, 0)
            self.assertIn("none is recorded", r.stderr)
        subprocess.run([sys.executable, CHECK, "decide", "--check", path, "--quantities", "raw",
                        "--by", "canalyst", "--reason", "the IgG lanes carry far less protein"],
                       check=True, capture_output=True, text=True, timeout=60)
        r = self.de(report, cond, contrast, et.RAW, final, "--normalization-check", path)
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        with open(os.path.join(final, "de_provenance.json")) as fh:
            n = json.load(fh)["normalization_check"]
        # the block keeps its schema: a one-element array stays an array
        self.assertEqual(n["design"]["vs_control"], ["IP-IgG"])
        self.assertIsInstance(n["trips"], list)
        self.assertEqual((n["decision"]["by"], n["decision"]["data_check_tripped"]),
                         ("Core staff", True))
        with open(os.path.join(final, "methods.txt")) as fh:
            methods = fh.read()
        self.assertIn("Quantities    : non-normalised (--quantities raw)", methods)
        self.assertIn("Core staff chose non-normalised (raw) quantities because the IgG lanes carry "
                      "far less protein.", " ".join(methods.split()))
        # who chose, by name, is staff-only: never in the DE's (delivered) record or Methods
        with open(os.path.join(final, "de_provenance.json")) as fh:
            self.assertNotIn("canalyst", fh.read() + methods)
        with open(os.path.join(final, "reproducibility_log.R")) as fh:
            self.assertIn("setNames(list('Precursor.Quantity'), int_arg)", fh.read())

    def test_an_ip_typed_as_ip_agrees_with_non_normalised(self):
        p, rec, md, (report, cond, contrast), path = self.check("ip", "ip", bait="BAITX")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual(rec["decision"]["quantities"], et.RAW)
        self.assertTrue(rec["checks"]["factor_confound"]["large_and_confounded"])
        self.assertEqual(rec["volcano"]["raw"]["contrasts"]["IP-IgG"]["bait"]["rank"], 1)
        r = self.de(report, cond, contrast, et.NORMALISED, os.path.join(self.d, "final"),
                    "--normalization-check", path)
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("never a silent switch", r.stderr)

    def test_the_one_added_sequence_is_the_bait_when_none_is_recorded(self):
        """Integration 2.10: the bait came only from set --bait; an IP searched with its bait
        added (--add-fasta) now has its bait checked with nothing more said."""
        p, rec, md, *_ = self.check("ip", "ip", added=("B00001",))
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual((rec["bait"], rec["bait_candidates"]), ("B00001", ["B00001"]))
        self.assertIn("--add-fasta", rec["bait_source"])
        self.assertEqual(rec["volcano"]["raw"]["contrasts"]["IP-IgG"]["bait"]["rank"], 1)
        self.assertIn("- **Bait:** B00001 (the one sequence the user added", md)

    def test_turboid_with_good_controls_passes_and_reports_the_carboxylases(self):
        p, rec, md, *_ = self.check("turbo", "proximity", bait="BAITX", controls="Ctrl")
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        self.assertEqual(rec["decision"]["quantities"], et.NORMALISED)
        cx = rec["carboxylases"]
        for k in (et.NORMALISED, et.RAW):
            self.assertEqual(cx[k]["n_found"], 5, cx)
            self.assertIsNotNone(cx[k]["median_sd_log2"])
        self.assertIn("never a normaliser", md)
        self.assertEqual(rec["volcano"]["normalised"]["contrasts"]["Bait-Ctrl"]["bait"]["rank"], 1)

    def test_turboid_whose_controls_came_out_low_still_trips(self):
        p, rec, *_ = self.check("turbo_low", "proximity", bait="BAITX", controls="Ctrl")
        self.assertEqual(p.returncode, nc.TRIPPED_EXIT, p.stdout + p.stderr)
        self.assertTrue(any(t.startswith("B1 ") for t in rec["trips"]), rec["trips"])
        self.assertIn("biotinylated carboxylases are steadier with normalised", rec["recommendation"])

    def test_maxlfq_raw_needs_a_no_norm_report_and_then_skips_the_quantile_step(self):
        report, cond, contrast = make_report(os.path.join(self.d, "n"), "whole")
        r = self.de(report, cond, contrast, et.RAW, os.path.join(self.d, "m1"), "--method", "maxlfq")
        self.assertNotEqual(r.returncode, 0)
        self.assertIn("NORMALISED protein quantity", r.stderr)
        # ... and says what to do next: dpc, or a --no-norm report given to the check
        flat = " ".join(r.stderr.split())
        self.assertIn("(1) --method dpc", flat)
        self.assertIn("--report-raw <no_norm_report.parquet>", flat)
        report, cond, contrast = make_report(os.path.join(self.d, "nn"), "whole", no_norm=True)
        out = os.path.join(self.d, "m2")
        r = self.de(report, cond, contrast, et.RAW, out, "--method", "maxlfq")
        self.assertEqual(r.returncode, 0, r.stderr[-2000:])
        with open(os.path.join(out, "de_provenance.json")) as fh:
            prov = json.load(fh)
        self.assertIn("no quantile step", prov["normalisation"])
        self.assertEqual(prov["normalization_check"]["status"], "not_run")
        with open(os.path.join(out, "reproducibility_log.R")) as fh:
            self.assertNotIn("normalizeBetweenArrays", fh.read())


if __name__ == "__main__":
    unittest.main()
