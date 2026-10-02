#!/usr/bin/env python3
"""
Independent 2.10 review (wrong-results risks), the part fixed on fix/2.10-review-a: which run is
which sample on an HT plate, and re-injections kept as technical replicates reaching the DE. Each
test pins a path that finished with exit 0 while a sample carried another sample's run, or two
injections of one sample counted as two samples. Synthetic names only. (The review's other two
classes, on collect_conditions' identifier columns and numbered conditions, live with their fix.)
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from test_core_submission import (child_env, raw, record, run, tims, tsv_rows,  # noqa: E402
                                  write_summary)


class TestHtPlateWellIsNotASampleId(unittest.TestCase):
    """--files-from on an HT plate: a unique_id is matched against the WHOLE run name before the
    sample_name fallback is tried, so a submitter who labelled tubes by their own wells (A1..A4)
    gets the run in that well of the CORE's plate (`_S5-A1_`) -- every sample matched, exit 0, no
    gate. The HT field (ht_sample_field) is where the sample is; the plate position never is."""

    S = [("A1", "Liver", "ctrl"), ("A2", "Heart", "ctrl"), ("A3", "Brain", "treat"),
         ("A4", "Lung", "treat")]
    CORE_PLATE = ["Heart", "Brain", "Lung", "Liver"]          # the Core's wells A1..A4

    def test_each_sample_gets_the_run_that_carries_its_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = child_env(tmp)
            rec = record(samples=[(u, c) for u, _n, c in self.S])
            for smp, (_u, n, _c) in zip(rec["submission_data"]["samples"], self.S):
                smp["sample_name"] = n
            summary, _ = write_summary(tmp, rec)
            files = {n: raw(tmp, "tTOF_HT", "mar25",
                            f"20250315_PROT_0807_100spd_{n}_S5-A{i + 1}_1_{24100 + i}.d")
                     for i, n in enumerate(self.CORE_PLATE)}
            lst = os.path.join(tmp, "plate.txt")
            with open(lst, "w") as fh:
                fh.write("\n".join(files.values()) + "\n")
            out = os.path.join(tmp, "loc")
            rc, js, p = run(["locate", "--summary", summary, "--out", out, "--files-from", lst], env)
            got = {r["unique_id"]: r["file"] for r in tsv_rows(os.path.join(out, "sample_files.tsv"))
                   if r["status"] == "matched"}
            wrong = {u: os.path.basename(got[u]) for u, n, _c in self.S
                     if u in got and got[u] != files[n]}
            self.assertFalse(wrong and rc == 0,
                             f"exit {rc} with samples on another sample's run: {wrong}")


class TestReinjectionsAllReachTheDE(unittest.TestCase):
    """--reinjections all keeps every injection; `conditions` asks "average or block?" -- but the
    CSV it proposes has only File.Name,Group: nothing names the sample each run came from, so
    run_de.R --block has no column, run_de.R's within-group repeat note cannot fire, there is no
    averaging route, and no flag records the answer (re-running asks again). Run as proposed,
    the injections count as independent samples (n doubled)."""

    def test_the_proposed_csv_names_each_runs_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = child_env(tmp)
            ids = ["RJ1", "RJ2", "RJ3", "RJ4"]
            summary, _ = write_summary(tmp, record(samples=[(u, "ctrl" if i < 2 else "treat")
                                                            for i, u in enumerate(ids)]))
            for i, u in enumerate(ids):
                tims(tmp, "03122025", u, 10 + i)
                tims(tmp, "03122025", u, 50 + i)
            out = os.path.join(tmp, "loc")
            rc, _js, p = run(["locate", "--summary", summary, "--out", out,
                              "--reinjections", "all"], env)
            self.assertEqual(rc, 0, p.stderr)
            csv_path = os.path.join(tmp, "conditions.csv")
            run(["conditions", "--summary", summary, "--sample-files",
                 os.path.join(out, "sample_files.tsv"), "--out", csv_path], env)
            with open(csv_path, newline="") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual(len(rows), 8)
            per_sample = {u: [r for r in rows if f"DIA-{u}_" in r["File.Name"]] for u in ids}
            extra = [c for c in rows[0] if c not in ("File.Name", "Group")]
            blockable = [c for c in extra
                         if all(len({r[c] for r in v}) == 1 for v in per_sample.values())
                         and len({v[0][c] for v in per_sample.values()}) == len(ids)]
            self.assertTrue(blockable, f"conditions.csv columns {list(rows[0])}: none names the "
                                       f"sample, so --block cannot keep 2 injections of one sample "
                                       f"from counting as 2 samples")


if __name__ == "__main__":
    unittest.main()
