"""Rank-matched comparison: precision of each arm's top-N predictions by score.

Usage: rank_score.py <report.parquet> <run_name> <label=file.ssl> [...]
Score scales differ between models (a fine-tune can shift every score below 0.9), so a fixed
threshold is not a fair comparison. This takes each arm's N highest-scoring predictions instead.
"""
import re
import sys

import pandas as pd
import pyarrow.parquet as pq

report, run = sys.argv[1:3]
arms = [a.split("=", 1) for a in sys.argv[3:]]
cols = [c for c in ("Run", "Stripped.Sequence", "Q.Value", "Decoy") if c in pq.read_schema(report).names]
dn = pq.read_table(report, columns=cols).to_pandas()
dn = dn[(dn["Run"] == run) & (dn["Q.Value"] <= 0.01)]
if "Decoy" in dn:
    dn = dn[dn["Decoy"] == 0]
truth = set(dn["Stripped.Sequence"].str.replace("I", "L"))


def strip(seq):
    return re.sub(r"[^A-Z]", "", re.sub(r"\[.*?\]", "", str(seq))).replace("I", "L")


NS = (100, 300, 1000, 3000, 10000)
print(f"{run}  (DIA-NN peps {len(truth)})")
print(f"{'arm':12s} " + " ".join(f"{'top' + str(n):>16s}" for n in NS) + "   score@top1000")
for label, path in arms:
    s = pd.read_csv(path, sep="\t").sort_values("score", ascending=False)
    s["hit"] = s["sequence"].map(strip).isin(truth)
    cells = []
    for n in NS:
        t = s.head(n)
        cells.append(f"{t['hit'].mean():6.1%} ({t.loc[t['hit'], 'sequence'].map(strip).nunique():5d})")
    print(f"{label:12s} " + " ".join(f"{c:>16s}" for c in cells) + f"   {s['score'].iloc[min(999, len(s) - 1)]:.3f}")
print()
