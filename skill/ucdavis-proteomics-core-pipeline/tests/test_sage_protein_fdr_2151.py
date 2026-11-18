#!/usr/bin/env python3
"""
2.15.1 audit (2026-11-17), the Sage quantified route (L2-4, L2-6).

L2-4  adapt_sage kept every peptide whatever its protein's FDR, and a peptide shared by several
      proteins became its own "protein group" (P1;P2) beside P1's unique peptides, so one protein
      could be tested twice. Now: protein groups at Sage's picked-protein q <= 0.01
      (results.sage.parquet protein_q, the identifications-only route's rule) and one stated
      shared-peptide rule (run_search.SAGE_SHARED_RULE; from the 2.15.1 review, protein-group
      parsimony with razor peptides, run_search.sage_protein_groups), both recorded in sage_adapt.json and the
      Methods sentence.
L2-6  Sage's LFQ integrates each identified peptide's MS1 signal in every run at its aligned
      retention time (lfq.rs, v0.14.7). Now each row says whether the peptide was identified by
      MS2 in THAT run (Identified.In.Run, the identifications-only route's own per-run rule).
      Wording corrected in 2.15.1: a value in a run without an identification is MEASURED by
      matching between runs (sage_lfq_check.MATCHED_BETWEEN_RUNS), never "transferred"; only its
      identity comes from other runs, with no per-run FDR. Rows are counted three ways
      (classify_rows): identified in the run, matched between runs (signal > 0), no MS1 signal
      (0 / missing -- never matched, never present).
Older results.sage.parquet without protein_q / spectrum_q: nothing filtered, said so.

Synthetic Sage 0.14.7-shaped outputs only. Needs pyarrow.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

try:
    import pyarrow as pa
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False

import sage_lfq_check as slc  # noqa: E402

RUNS = ["r1", "r2", "r3", "r4"]
# peptide -> Sage's proteins string. P1 has its own peptides; PEPS1 is shared P1;P2 (P2 has a
# peptide of its own too); PEPI1/2 belong to two indistinguishable isoforms I1;I2 (neither has
# its own peptide); PEPLOW1 belongs to a protein above 1% protein FDR.
PEPTIDES = {
    "PEPAAK": "sp|P1|ONE_HUMAN", "PEPABK": "sp|P1|ONE_HUMAN",
    "PEPBAK": "sp|P2|TWO_HUMAN",
    "PEPSHK": "sp|P1|ONE_HUMAN;sp|P2|TWO_HUMAN",
    "PEPIAK": "sp|I1|ISO1_HUMAN;sp|I2|ISO2_HUMAN", "PEPIBK": "sp|I1|ISO1_HUMAN;sp|I2|ISO2_HUMAN",
    "PEPLOK": "sp|L1|LOW_HUMAN",
}
PROTEIN_Q = {"sp|L1|LOW_HUMAN": 0.2}            # every other group 0.001
# (run, peptide) pairs with a confident rank-1 PSM in that run; the rest were not identified there
IDENTIFIED = {(r, p) for r in RUNS for p in PEPTIDES} - {("r4", "PEPAAK"), ("r3", "PEPBAK")}
# target (run, peptide) intensities Sage left at 0.0 / null: no MS1 signal in that run. One not
# identified there (r3 PEPBAK), one identified with 0 (r2 PEPABK), one identified with null.
NO_SIGNAL = {("r3", "PEPBAK"): 0.0, ("r2", "PEPABK"): 0.0, ("r1", "PEPIBK"): None}


def write_search(out, protein_q=True, spectrum_q=True):
    os.makedirs(out, exist_ok=True)
    res = {k: [] for k in ("filename", "is_decoy", "rank", "peptide", "proteins", "peptide_q",
                           "spectrum_q", "protein_q", "precursor_ppm", "expmass", "calcmass",
                           "isotope_error")}
    for r in RUNS:
        for i, (pep, prot) in enumerate(PEPTIDES.items()):
            for decoy in (False, True):
                if not decoy and (r, pep) not in IDENTIFIED:
                    # a poor PSM in this run: spectrum_q fails, so not this run's identification
                    sq = 0.3
                else:
                    sq = 0.001
                res["filename"].append(f"{r}.mzML")
                res["is_decoy"].append(decoy)
                res["rank"].append(1)
                res["peptide"].append(pep)
                res["proteins"].append(("rev_" if decoy else "") + prot)
                res["peptide_q"].append(0.001 if not decoy else 0.9)
                res["spectrum_q"].append(sq if not decoy else 0.9)
                res["protein_q"].append(PROTEIN_Q.get(prot, 0.001) if not decoy else 0.9)
                res["precursor_ppm"].append(1.0)
                cm = 1000.0 + 10 * i
                res["expmass"].append(cm)
                res["calcmass"].append(cm)
                res["isotope_error"].append(0.0)
    if not protein_q:
        del res["protein_q"]
    if not spectrum_q:
        del res["spectrum_q"]
    pq.write_table(pa.table(res), os.path.join(out, "results.sage.parquet"))
    lf = {k: [] for k in ("peptide", "stripped_peptide", "charge", "proteins", "is_decoy",
                          "q_value", "filename", "intensity")}
    for pep, prot in PEPTIDES.items():
        for decoy in (False, True):
            for j, r in enumerate(RUNS):
                lf["peptide"].append(pep)
                lf["stripped_peptide"].append(pep)
                lf["charge"].append(None)
                lf["proteins"].append(prot)
                lf["is_decoy"].append(decoy)
                lf["q_value"].append(0.001)
                lf["filename"].append(f"{r}.mzML")
                lf["intensity"].append(NO_SIGNAL[(r, pep)] if not decoy and (r, pep) in NO_SIGNAL
                                       else 1e6 * (j + 1) * (0.3 if decoy else 1.0))
    pq.write_table(pa.table({**{k: v for k, v in lf.items()
                                if k not in ("charge", "q_value", "intensity")},
                             "charge": pa.array(lf["charge"], pa.int32()),
                             "q_value": pa.array(lf["q_value"], pa.float32()),
                             "intensity": pa.array(lf["intensity"], pa.float32())}),
                   os.path.join(out, "lfq.parquet"))
    with open(os.path.join(out, "results.json"), "w") as fh:
        json.dump({"version": "0.14.6", "quant": {"lfq": True}}, fh)


class SharedPeptideRule(unittest.TestCase):
    """2.15.1 review finding 2: shared peptides are grouped by parsimony, never dropped."""

    def groups(self, g):
        import run_search as rs
        return rs.sage_protein_groups(g)

    def test_the_rule(self):
        got, rec = self.groups({
            "a": "P1", "b": "P1;P2", "c": "P2",          # razor: b goes to P1 (more peptides)
            "a2": "P1",
            "d": "I1;I2", "e": "I1;I2",                  # indistinguishable: one group
            "f": "X1;X2", "g": "X2;X3",                  # GSTM-like: X1, X3 fold into X2
            "h": "P1;Cont_P02769",                       # the contaminant rule's, not this one's
        })
        self.assertEqual(got, {"a": "P1", "a2": "P1", "b": "P1", "c": "P2",
                               "d": "I1;I2", "e": "I1;I2",
                               "f": "X2;X1;X3", "g": "X2;X1;X3",
                               "h": "P1;Cont_P02769"})
        self.assertEqual(rec["razor_peptides"], 1)
        self.assertEqual(rec["indistinguishable_groups"], 1)
        self.assertEqual(rec["subset_proteins_folded"], 2)
        self.assertEqual(rec["contaminant_shared_peptides"], 1)
        self.assertEqual(rec["groups_not_quantified"], 0)
        self.assertEqual(rec["protein_groups_quantified"], 4)

    def test_a_histone_like_family_is_quantified_not_dropped(self):
        # no member has an entry-unique peptide except H2AX; H2A1's peptides are all H2AX's
        got, rec = self.groups({
            "a": "H2A1;H2A2;H2A3;H2AX", "b": "H2A1;H2A2;H2AX", "c": "H2A1;H2A3;H2AX",
            "d": "H2A2;H2A3", "e": "H2AX",
        })
        self.assertEqual(got["a"], "H2AX;H2A1")
        self.assertEqual(got["b"], "H2AX;H2A1")
        self.assertEqual(got["c"], "H2AX;H2A1")
        self.assertEqual(got["e"], "H2AX;H2A1")
        # d: H2A2 and H2A3 tie on peptides (3) and own peptides (0) -> the first accession
        self.assertEqual(got["d"], "H2A2")
        # H2A3 keeps no peptide: not quantified, and named with where its peptides went
        self.assertEqual(rec["groups_not_quantified"], 1)
        self.assertEqual(rec["groups_not_quantified_list"],
                         [{"group": "H2A3", "peptides": 3,
                           "peptides_quantified_in": ["H2A2", "H2AX;H2A1"]}])

    def test_overlapping_lists_without_own_peptides_keep_the_family(self):
        # the review's GH1/GH2/CSH case: three lists, every one overlapping another
        got, rec = self.groups({"p": "GH1;GH2", "q": "GH1;CSH1", "r": "GH2;CSH1;CSH2",
                                "s": "CSH1;CSH2"})
        self.assertEqual(set(got), {"p", "q", "r", "s"})       # every peptide quantified
        quantified = {a for g in got.values() for a in g.split(";")}
        listed = {a for x in rec["groups_not_quantified_list"] for a in x["group"].split(";")}
        self.assertEqual(quantified | listed, {"GH1", "GH2", "CSH1", "CSH2"})

    def test_input_order_does_not_change_the_groups(self):
        import random
        g = {"a": "H2A1;H2A2;H2A3;H2AX", "b": "H2A1;H2A2;H2AX", "c": "H2A1;H2A3;H2AX",
             "d": "H2A2;H2A3", "e": "H2AX", "f": "X1;X2", "g": "X2;X3", "k": "P9"}
        ref = self.groups(g)
        items = list(g.items())
        for seed in range(5):
            random.Random(seed).shuffle(items)
            shuffled = {k: ";".join(reversed(v.split(";"))) for k, v in items}
            self.assertEqual(self.groups(shuffled), ref)

    def test_nothing_vanishes_on_a_simulated_digest_with_protein_families(self):
        """A synthetic proteome with histone-like families (members share most peptides, few or
        none their own), random detection: every protein with a peptide is in a quantified
        group, or named in groups_not_quantified_list, or in the contaminant list."""
        import random
        rnd = random.Random(7)
        proteins = {}
        pool = [f"PEP{i}K" for i in range(4000)]
        n = 0
        for fam in range(120):                       # families of 2-8 near-identical members
            core = rnd.sample(pool, 6)
            for m in range(rnd.randint(2, 8)):
                own = [f"F{fam}M{m}U{j}K" for j in range(rnd.choice((0, 0, 1, 2)))]
                proteins[f"FAM{fam}_{m}"] = set(rnd.sample(core, rnd.randint(2, 6))) | set(own)
                n += 1
        for i in range(600):                          # ordinary proteins, some sharing
            proteins[f"P{i}"] = {f"U{i}_{j}K" for j in range(rnd.randint(1, 8))} | \
                set(rnd.sample(pool, rnd.choice((0, 0, 1))))
        proteins["Cont_P02769"] = {"PEP1K", "CONTK"}
        pep2 = {}
        for a, ps in proteins.items():
            for p in ps:
                if rnd.random() < 0.6:
                    pep2.setdefault(p, set()).add(a)
        g = {p: ";".join(sorted(a)) for p, a in pep2.items()}
        got, rec = self.groups(g)
        self.assertEqual(set(got), set(g))           # every peptide quantified somewhere
        with_evidence = set().union(*pep2.values())
        in_rows = {a for grp in got.values() for a in grp.split(";")}
        listed = {a for x in rec["groups_not_quantified_list"] for a in x["group"].split(";")}
        self.assertEqual(with_evidence - in_rows - listed, set())
        self.assertEqual(rec["groups_not_quantified"], len(rec["groups_not_quantified_list"]))
        # one peptide, one group: no protein is quantified as two rows
        member = {}
        for grp in set(got.values()):
            if "Cont_" in grp:
                continue
            for a in grp.split(";"):
                self.assertNotIn(a, member, f"{a} in {member.get(a)} and {grp}")
                member[a] = grp
        self.assertGreater(rec["subset_proteins_folded"], 0)
        self.assertGreater(rec["razor_peptides"], 0)


@unittest.skipUnless(HAVE_ARROW, "needs pyarrow")
class AdaptSage(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        import run_search as rs
        self.rs = rs
        self.addCleanup(setattr, rs, "ENGINE_VERSION", rs.ENGINE_VERSION)
        rs.ENGINE_VERSION = None

    def adapt(self, **kw):
        write_search(self.d, **kw)
        report = self.rs.adapt_sage(self.d)
        with open(os.path.join(self.d, self.rs.SAGE_ADAPT)) as fh:
            rec = json.load(fh)
        return pq.read_table(report), rec

    def test_protein_fdr_shared_peptides_and_matches(self):
        t, rec = self.adapt()
        rows = t.to_pylist()
        groups = {r["Protein.Group"] for r in rows}
        # L1 (protein q 0.2) out; the P1;P2 peptide is P1's razor peptide; I1;I2 one group
        self.assertEqual(groups, {"P1", "P2", "I1;I2"})
        self.assertEqual({r["Protein.Group"] for r in rows if r["Peptide.Id"] == "PEPSHK"},
                         {"P1"})
        self.assertEqual({r["Protein.Names"] for r in rows if r["Protein.Group"] == "I1;I2"},
                         {"ISO1_HUMAN;ISO2_HUMAN"})
        self.assertEqual(rec["protein_q_max"], 0.01)
        self.assertEqual(rec["totals"]["protein_q_dropped"], len(RUNS))      # PEPLOK x 4 runs
        self.assertEqual(len(rows), 6 * len(RUNS))       # no shared peptide is left out
        pg = rec["protein_grouping"]
        self.assertEqual((pg["razor_peptides"], pg["indistinguishable_groups"],
                          pg["groups_not_quantified"]), (1, 1, 0))
        self.assertEqual(pg["rule"], rec["shared_peptide_rule"])
        # each value says whether it was identified in its run
        flag = {(r["Run"], r["Peptide.Id"]): r["Identified.In.Run"] for r in rows}
        self.assertEqual(flag[("r4", "PEPAAK")], 0)
        self.assertEqual(flag[("r3", "PEPBAK")], 0)
        self.assertEqual(flag[("r1", "PEPAAK")], 1)
        self.assertNotIn("transferred_rows", rec)
        meta = t.schema.metadata
        self.assertEqual(meta[b"delimp.quantity.engine_normalised"], b"false")
        ident = meta[b"delimp.quantity.identification"].decode()
        self.assertIn("picked-protein q-value <= 0.01", ident)
        self.assertIn(slc.MATCHED_BETWEEN_RUNS, ident)
        self.assertIn(slc.NO_MS1_SIGNAL, ident)
        # the one definition reaches run_de.R through the declaration
        self.assertEqual(meta[b"delimp.quantity.matched_between_runs"].decode(),
                         slc.MATCHED_BETWEEN_RUNS)
        s = rec["methods_sentence"]
        self.assertIn("protein groups were kept at Sage's picked-protein q-value <= 0.01", s)
        self.assertIn("Shared peptides:", s)
        self.assertIn(slc.MATCHED_BETWEEN_RUNS, s)
        self.assertIn("count as present whether identified or matched", s)
        self.assertNotIn("not counted as detections", s)
        self.assertNotIn("no protein-level FDR filter was applied", s)
        for text in (ident, s, rec["identified_in_run"]):
            for word in ("transferred", "filled in", "copied"):
                self.assertNotIn(word, text)

    def test_three_way_count_never_calls_a_zero_matched(self):
        """identified in the run / matched between runs (signal > 0) / no MS1 signal (0 or
        null), counted separately; a no-signal row is never matched, identified or not."""
        t, rec = self.adapt()
        ev = rec["row_evidence"]
        self.assertEqual(ev["rows"], len(t))
        self.assertEqual(ev["identified_in_run"] + ev["matched_between_runs"]
                         + ev["no_ms1_signal"], ev["rows"])
        # r4 PEPAAK: not identified there, signal -> the only matched row; r3 PEPBAK: not
        # identified there and 0 -> no signal, NOT matched
        self.assertEqual(ev["matched_between_runs"], 1)
        self.assertEqual(ev["no_ms1_signal"], 3)
        self.assertEqual(ev["no_ms1_signal_identified_in_run"], 2)
        self.assertEqual(ev["identified_in_run"], 6 * len(RUNS) - 1 - 3)
        self.assertEqual(ev["matched_between_runs_rule"], slc.MATCHED_BETWEEN_RUNS)
        self.assertEqual(ev["no_ms1_signal_rule"], slc.NO_MS1_SIGNAL)
        self.assertIn(slc.row_evidence_sentence(ev), rec["methods_sentence"])
        self.assertIn(f"{ev['matched_between_runs']} were measured by matching between runs",
                      rec["methods_sentence"])

    def test_an_older_results_file_filters_nothing_and_says_so(self):
        t, rec = self.adapt(protein_q=False, spectrum_q=False)
        self.assertNotIn("Identified.In.Run", t.column_names)
        self.assertIn("P1", set(t.column("Protein.Group").to_pylist()))
        self.assertIn("L1", set(t.column("Protein.Group").to_pylist()))
        self.assertIsNone(rec["protein_q_max"])
        self.assertTrue(rec["protein_rule"].startswith("NOT APPLIED"))
        self.assertTrue(rec["identified_in_run"].startswith("NOT RECORDED"))
        self.assertIn("no protein-level FDR filter was applied", rec["methods_sentence"])
        # what a value in a run without an identification is, said even when the report cannot
        # mark which values those are; the no-signal rows are still counted
        self.assertIn(slc.MATCHED_BETWEEN_RUNS, rec["methods_sentence"])
        self.assertIsNone(rec["row_evidence"]["matched_between_runs"])
        self.assertEqual(rec["row_evidence"]["no_ms1_signal"], 3)
        self.assertTrue(rec["row_evidence"]["not_recorded"])


class MethodsSentence(unittest.TestCase):
    def test_default_call_says_what_a_matched_value_is(self):
        s = slc.methods_sentence()
        self.assertIn("no protein-level FDR filter was applied.", s)
        self.assertIn(slc.MATCHED_BETWEEN_RUNS, s)
        self.assertIn("Sage gives no per-run FDR for that identity", s)
        self.assertIn(slc.NO_MS1_SIGNAL, s)
        self.assertNotIn("transferred", s)


class ClassifyRows(unittest.TestCase):
    """sage_lfq_check.classify_rows: the one three-way count."""

    def test_zero_null_nan_and_negative_are_no_signal(self):
        n = slc.classify_rows([1, 0, 0, 1, 0, 0, 1],
                              [5.0, 3.0, 0.0, 0.0, None, float("nan"), -1.0])
        self.assertEqual((n["identified_in_run"], n["matched_between_runs"], n["no_ms1_signal"],
                          n["no_ms1_signal_identified_in_run"], n["rows"]), (1, 1, 5, 2, 7))

    def test_sentence_needs_a_full_count(self):
        self.assertEqual(slc.row_evidence_sentence(None), "")
        self.assertEqual(slc.row_evidence_sentence({"rows": 3, "no_ms1_signal": 1,
                                                    "identified_in_run": None,
                                                    "matched_between_runs": None}), "")
        s = slc.row_evidence_sentence(slc.classify_rows([1, 0, 0], [2.0, 2.0, 0.0]))
        self.assertEqual(s, "Of the 3 peptide × run values kept, 1 were identified by MS2 in "
                            "their run, 1 were measured by matching between runs, and 1 had no "
                            "MS1 signal (missing values).")


if __name__ == "__main__":
    unittest.main()
