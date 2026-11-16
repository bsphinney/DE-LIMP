#!/usr/bin/env python3
"""The FragPipe DDA adapter's quantity (2.15 F4; Rulings 64, 90; review H4, H6).

run_search.adapt_fragpipe_dda() in 2.14.2 took EVERY column ending in "MaxLFQ Intensity", else
every column ending in "Intensity" -- so a Unique / Total variant became a sample of its own --
never checked the runs against the search's run map, and with 2 experiments shipped IonQuant's
protein Intensity, which a deeper run inflates (+1.67 log2 against directLFQ's +0.14 on a
two-run case). Now (fragpipe_quant.plan, the one quantity object):
  * exactly "<exp> MaxLFQ Intensity" for every experiment, else exactly "<exp> Intensity";
  * a run of the map missing from combined_protein.tsv, or with no quantity at all, stops the
    adapter, naming it, and an earlier report.parquet is set aside;
  * 2 experiments: PG.MaxLFQ = directLFQ (its guarded record), PG.IonQuant.Intensity beside it;
    no directLFQ venv or a refused guard = no report;
  * the declaration says what the quantity is, n_experiments, why it is not MaxLFQ (IonQuant's
    own log line), and a caveat when there is no run map;
  * fragpipe_quantity.json + search_provenance.json record it.
"""
import csv
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

import fragpipe_quant as fq                     # noqa: E402

try:
    import pyarrow.parquet as pq
    HAVE_ARROW = True
except ImportError:
    HAVE_ARROW = False


PROTEINS = [("sp|P02768|ALBU_HUMAN", "P02768", "ALBU_HUMAN", "ALB"),
            ("sp|P60709|ACTB_HUMAN", "P60709", "ACTB_HUMAN", "ACTB"),
            ("sp|Q99999|ONLYFP_HUMAN", "Q99999", "ONLYFP_HUMAN", "ONLYFP")]
LOG_LINE = "[WARNING] - There are only 2 experiments. Will not calculate MaxLFQ intensity."


def write_output(d, runs, maxlfq=True, variants=True, table_runs=None, zero_run=None):
    """A FragPipe output: combined_protein.tsv (+ Unique/Total variants with absurd values),
    the skill's fragpipe.fp-manifest (run names with dashes) and FragPipe's
    fragpipe-files.fp-manifest ('-' -> '_'), fragpipe.workflow. -> {run: FragPipe name}."""
    exps = {r: r.replace("-", "_") for r in runs}
    in_table = [exps[r] for r in (table_runs or runs)]
    hdr = ["Protein", "Protein ID", "Entry Name", "Gene"]
    for e in in_table:
        hdr += [f"{e} Spectral Count", f"{e} Intensity"]
        if variants:
            hdr += [f"{e} Unique Intensity", f"{e} Total Intensity"]
        if maxlfq:
            hdr += [f"{e} MaxLFQ Intensity"]
            if variants:
                hdr += [f"{e} Unique MaxLFQ Intensity"]
    hdr += ["Indistinguishable Proteins"]
    with open(os.path.join(d, "combined_protein.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(hdr)
        for i, (prot, pid, entry, gene) in enumerate(PROTEINS):
            row = [prot, pid, entry, gene]
            for j, e in enumerate(in_table):
                inten = 0 if zero_run and e == exps[zero_run] else 1000.0 * (i + 1) * (j + 1)
                row += ["3", str(inten)]
                if variants:
                    row += ["999999999", "888888888"]
                if maxlfq:
                    row += [str(inten * 2)]
                    if variants:
                        row += ["777777777"]
            w.writerow(row + [""])
    st = os.path.join(d, "raw_stage")
    with open(os.path.join(d, "fragpipe.fp-manifest"), "w") as fh:
        fh.write("".join(f"{st}/{r}.d\t{r}\t\tDDA\n" for r in runs))
    with open(os.path.join(d, "fragpipe-files.fp-manifest"), "w") as fh:
        fh.write("".join(f"{st}/{r}.d\t{exps[r]}\t\tDDA\n" for r in runs))
    with open(os.path.join(d, "fragpipe.workflow"), "w") as fh:
        fh.write("# FragPipe version 24.0\n# IonQuant version 1.11.20\nionquant.mbr=1\n"
                 "ionquant.normalization=1\nphi-report.filter=--sequential --prot 0.01 --picked\n")
    with open(os.path.join(d, "search_provenance.json"), "w") as fh:
        json.dump({"engine": "fragpipe", "version": "24.0"}, fh)
    return exps


def write_directlfq(d, exps, values, status="ok", passed=True):
    """combined_ion.tsv and a directLFQ record for it, as fragpipe_quant.py directlfq writes."""
    with open(os.path.join(d, fq.ION_TABLE), "w") as fh:
        fh.write("Protein\tModified Sequence\t" + "\t".join(f"{e} Intensity" for e in exps)
                 + "\n")
    os.makedirs(os.path.join(d, fq.DIRECTLFQ_DIR), exist_ok=True)
    rel = f"{fq.DIRECTLFQ_DIR}/combined_ion.tsv.protein_intensities.tsv"
    with open(os.path.join(d, rel), "w") as fh:
        fh.write("protein\t" + "\t".join(exps) + "\n")
        for prot, per in values.items():
            fh.write(prot + "\t" + "\t".join(str(per.get(e, 0)) for e in exps) + "\n")
    rec = {"status": status, "pins": list(fq.DIRECTLFQ_PINS), "samples": list(exps),
           "input": fq.ION_TABLE, "input_sha256": fq.sha256_file(os.path.join(d, fq.ION_TABLE)),
           "outputs": {"normalised": rel}, "proteins_quantified": len(values),
           "versions": {"directlfq": "0.3.3", "pandas": "2.2.3"},
           "guard": {"passed": passed,
                     "refusal": None if passed else "REFUSED: directLFQ's result is not used "
                                                    "-- pandas 3.1.1 (test)",
                     "checks": {"shift_applied_equals_found": {
                         "found": {e: 0.0 for e in exps}, "max_abs_difference": 0.0}}}}
    with open(os.path.join(d, fq.DIRECTLFQ_DIR, fq.DIRECTLFQ_RECORD), "w") as fh:
        json.dump(rec, fh)


@unittest.skipUnless(HAVE_ARROW, "pyarrow")
class _Case(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, self.d, True)
        env = mock.patch.dict(os.environ, {fq.VENV_ENV: os.path.join(self.d, "no_venv")})
        env.start()
        self.addCleanup(env.stop)

    def adapt(self):
        import run_search
        report = run_search.adapt_fragpipe_dda(self.d)
        t = pq.read_table(report)
        meta = {k.decode(): v.decode() for k, v in (t.schema.metadata or {}).items()}
        return t.to_pydict(), {k[len("delimp.quantity."):]: v for k, v in meta.items()
                                if k.startswith("delimp.quantity.")}

    def refused(self):
        import run_search
        with self.assertRaises(SystemExit) as cm:
            run_search.adapt_fragpipe_dda(self.d)
        return str(cm.exception.code)


class ExactColumnsAndTheRunCheck(_Case):
    def test_maxlfq_columns_exactly_and_the_runs_keep_their_dashes(self):
        write_output(self.d, ["S-1", "S-2", "S-3"])
        t, decl = self.adapt()
        self.assertEqual(sorted(set(t["Run"])), ["S-1", "S-2", "S-3"])
        vals = {(r, g): v for r, g, v in zip(t["Run"], t["Protein.Group"], t["PG.MaxLFQ"])}
        self.assertEqual(vals[("S-2", "P02768")], 4000.0)          # MaxLFQ = 2 x Intensity
        self.assertNotIn(777777777.0, t["PG.MaxLFQ"])                # never Unique MaxLFQ
        self.assertEqual(decl["quantity"], fq.QUANTITY_MAXLFQ)
        self.assertEqual(decl["n_experiments"], "3")
        self.assertNotIn("why_not_maxlfq", decl)
        self.assertNotIn("caveat", decl)

    def test_unique_and_total_intensity_are_never_read(self):
        write_output(self.d, ["S-1", "S-2", "S-3"], maxlfq=False)
        write_directlfq(self.d, ["S_1", "S_2", "S_3"],
                        {"sp|P02768|ALBU_HUMAN": {"S_1": 5, "S_2": 6, "S_3": 7}})
        t, decl = self.adapt()
        self.assertEqual(sorted(set(t["Run"])), ["S-1", "S-2", "S-3"])
        for bad in (999999999.0, 888888888.0):
            self.assertNotIn(bad, t["PG.MaxLFQ"])
            self.assertNotIn(bad, t[fq.COMPANION_COLUMN])
        self.assertEqual(sorted(v for v in t[fq.COMPANION_COLUMN] if v == v),
                         sorted(1000.0 * (i + 1) * (j + 1) for i in range(3) for j in range(3)))

    def test_three_runs_without_maxlfq_are_directlfq(self):
        """Ruling 102 (3): IonQuant's Intensity is inflated by depth whatever the run count."""
        write_output(self.d, ["S-1", "S-2", "S-3"], maxlfq=False)
        write_directlfq(self.d, ["S_1", "S_2", "S_3"],
                        {"sp|P02768|ALBU_HUMAN": {"S_1": 5, "S_2": 6, "S_3": 7}})
        t, decl = self.adapt()
        self.assertEqual(decl["quantity"], fq.QUANTITY_DIRECTLFQ)
        self.assertEqual(decl["n_experiments"], "3")
        self.assertEqual(decl["companion_column"], fq.COMPANION_COLUMN)
        self.assertIn("Ammar et al. 2023", decl["citation"])
        self.assertEqual(sorted(v for v in t["PG.MaxLFQ"] if v == v), [5.0, 6.0, 7.0])

    def test_a_run_missing_from_the_table_stops_the_adapter(self):
        write_output(self.d, ["S-1", "S-2", "S-3"], table_runs=["S-1", "S-2"])
        old = os.path.join(self.d, "report.parquet")
        with open(old, "w") as fh:
            fh.write("an earlier search's report")
        msg = self.refused()
        self.assertIn("S_3", msg)
        self.assertFalse(os.path.exists(old))                       # set aside, never used
        self.assertTrue(any(f.startswith("report.parquet.stale-") for f in os.listdir(self.d)))
        with open(os.path.join(self.d, fq.RECORD)) as fh:
            self.assertEqual(json.load(fh)["status"], "refused")

    def test_a_run_with_no_quantity_stops_the_adapter(self):
        write_output(self.d, ["S-1", "S-2", "S-3"], zero_run="S-3")
        self.assertIn("S-3", self.refused())

    def test_without_a_map_the_declaration_carries_the_caveat(self):
        write_output(self.d, ["S1", "S2", "S3"])
        for m in ("fragpipe.fp-manifest", "fragpipe-files.fp-manifest"):
            os.remove(os.path.join(self.d, m))
        t, decl = self.adapt()
        self.assertEqual(sorted(set(t["Run"])), ["S1", "S2", "S3"])
        self.assertEqual(decl["caveat"], fq.NO_MAP_CAVEAT)
        self.assertEqual(decl["run_map"], "False")


class TwoExperimentsAreDirectLfq(_Case):
    def setUp(self):
        super().setUp()
        self.exps = write_output(self.d, ["S-1", "S-2"], maxlfq=False)
        with open(os.path.join(self.d, "log_2026-11-01_20-59-15.txt"), "w") as fh:
            fh.write(LOG_LINE + "\n")

    def test_directlfq_is_the_quantity_and_ionquant_intensity_is_beside_it(self):
        write_directlfq(self.d, ["S_1", "S_2"], {"sp|P02768|ALBU_HUMAN": {"S_1": 50, "S_2": 55},
                                                  "sp|P60709|ACTB_HUMAN": {"S_1": 7, "S_2": 8}})
        t, decl = self.adapt()
        rows = {(r, g): (v, c) for r, g, v, c in zip(t["Run"], t["Protein.Group"], t["PG.MaxLFQ"],
                                                    t[fq.COMPANION_COLUMN])}
        self.assertEqual(rows[("S-2", "P02768")], (55.0, 2000.0))
        self.assertEqual(rows[("S-1", "P60709")], (7.0, 2000.0))
        v, c = rows[("S-1", "Q99999")]                     # no shared ions: no directLFQ value
        self.assertNotEqual(v, v)
        self.assertEqual(c, 3000.0)
        self.assertEqual(decl["quantity"], fq.QUANTITY_DIRECTLFQ)
        self.assertEqual(decl["n_experiments"], "2")
        self.assertIn(LOG_LINE.split(" - ", 1)[1], decl["why_not_maxlfq"])
        self.assertEqual(decl["companion_column"], fq.COMPANION_COLUMN)
        self.assertIn("directLFQ 0.3.3", decl["value"])
        self.assertIn("Ammar et al. 2023", decl["citation"])
        self.assertNotIn("Cox et al. 2014", decl["citation"])
        self.assertIn("passed", decl["directlfq_guard"])
        with open(os.path.join(self.d, "search_provenance.json")) as fh:
            prov = json.load(fh)[fq.PROVENANCE_KEY]
        self.assertEqual(prov["quantity"], fq.QUANTITY_DIRECTLFQ)
        self.assertTrue(prov["directlfq"]["guard"]["passed"])

    def test_without_the_venv_there_is_no_report(self):
        with open(os.path.join(self.d, fq.ION_TABLE), "w") as fh:
            fh.write("Protein\tModified Sequence\tS_1 Intensity\tS_2 Intensity\n")
        msg = self.refused()
        self.assertIn(fq.BUILDER, msg)
        self.assertFalse(os.path.exists(os.path.join(self.d, "report.parquet")))

    def test_a_refused_guard_means_no_report(self):
        write_directlfq(self.d, ["S_1", "S_2"], {"sp|P02768|ALBU_HUMAN": {"S_1": 5, "S_2": 6}},
                        status="refused", passed=False)
        msg = self.refused()
        self.assertIn("REFUSED", msg)
        self.assertFalse(os.path.exists(os.path.join(self.d, "report.parquet")))

    def test_ionquant_intensity_is_never_the_two_run_quantity(self):
        write_directlfq(self.d, ["S_1", "S_2"], {"sp|P02768|ALBU_HUMAN": {"S_1": 50, "S_2": 55}})
        t, _ = self.adapt()
        quantified = [v for v in t["PG.MaxLFQ"] if v == v]
        self.assertEqual(sorted(quantified), [50.0, 55.0])


@unittest.skipUnless(HAVE_ARROW, "pyarrow")
class TheMapIsOne(unittest.TestCase):
    def test_the_renamed_names_come_from_the_full_map(self):
        import run_search
        with tempfile.TemporaryDirectory() as d:
            write_output(d, ["S-1", "S2"])
            m = run_search.fragpipe_experiment_map(d)
            self.assertEqual([(x["experiment"], x["run"]) for x in m],
                             [("S_1", "S-1"), ("S2", "S2")])
            self.assertEqual(run_search.fragpipe_experiment_names(d), {"S_1": "S-1"})
            os.remove(os.path.join(d, "fragpipe-files.fp-manifest"))
            self.assertIsNone(run_search.fragpipe_experiment_map(d))


if __name__ == "__main__":
    unittest.main()
