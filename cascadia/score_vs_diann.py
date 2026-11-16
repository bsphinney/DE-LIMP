"""Score Cascadia .ssl predictions against DIA-NN identifications from the same run.

Usage: score_vs_diann.py <report.parquet> <run_name> <label=file.ssl> [label=file.ssl ...]
A prediction "matches" when its stripped sequence (I->L, mods removed) is a DIA-NN q<=0.01 peptide of
that run; "RT-consistent" additionally needs |RT_pred - RT_DIA-NN| <= 1 min. Reported per score
threshold: predictions, distinct predicted peptides, distinct DIA-NN peptides recovered, precision.
"""
import re
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

report, run = sys.argv[1:3]
arms = [a.split("=", 1) for a in sys.argv[3:]]
cols = [c for c in ("Run", "Stripped.Sequence", "RT", "Q.Value", "Decoy") if c in pq.read_schema(report).names]
dn = pq.read_table(report, columns=cols).to_pandas()
dn = dn[(dn["Run"] == run) & (dn["Q.Value"] <= 0.01)]
if "Decoy" in dn:
    dn = dn[dn["Decoy"] == 0]
dn["pep"] = dn["Stripped.Sequence"].str.replace("I", "L")
truth = set(dn["pep"])
rt_by_pep = dn.groupby("pep")["RT"].apply(lambda s: np.array(sorted(s))).to_dict()
print(f"DIA-NN q<=0.01 in {run}: {len(truth)} distinct peptides\n")


def strip(seq):
    return re.sub(r"[^A-Z]", "", re.sub(r"\[.*?\]", "", str(seq))).replace("I", "L")


print(f"{'arm':12s} {'score>=':>7s} {'preds':>8s} {'pred peps':>9s} {'DIA-NN peps hit':>15s} "
      f"{'precision':>9s} {'RT-consistent':>13s}")
for label, path in arms:
    s = pd.read_csv(path, sep="\t")
    s["pep"] = s["sequence"].map(strip)
    rt = s["retention-time"].astype(float)
    s["rt_min"] = rt / 60.0 if rt.max() > 200 else rt  # write_results emits minutes only when max RT > 500 s
    s["hit"] = s["pep"].isin(truth)
    s["rt_ok"] = [h and np.min(np.abs(rt_by_pep[p] - r)) <= 1.0 if h else False
                  for p, r, h in zip(s["pep"], s["rt_min"], s["hit"])]
    for thr in (0.0, 0.5, 0.8, 0.9, 0.95, 0.99):
        t = s[s["score"] >= thr]
        if t.empty:
            continue
        print(f"{label:12s} {thr:7.2f} {len(t):8d} {t['pep'].nunique():9d} "
              f"{t.loc[t['hit'], 'pep'].nunique():15d} {t['hit'].mean():9.1%} {t['rt_ok'].mean():13.1%}")
    print()
