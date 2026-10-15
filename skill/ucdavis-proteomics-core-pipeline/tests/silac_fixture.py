#!/usr/bin/env python3
"""A tiny synthetic DIA-NN plexDIA report (triple SILAC), for the SILAC tests.

Deterministic (seeded). 4 runs (replicates), channels L/M/H holding WT/KO/DKO. 60 protein groups
with 2-5 precursors each; per run a global offset (the light channel is ~3x the others, as in an
equal-volume mix of unequal cultures), protein effects KO/WT = +2 for P01-P10 and -2 for P11-P20,
DKO/KO = 0 throughout, precursor noise sd 0.15. P51-P55 are in WT only (channel-specific), P56 has
a ratio in 2 runs only (not tested), a Cont_ protein and unlabelled precursors (Channel "") are
present, and one channel row per run fails its channel q-value. Precursor ids carry the label once
-- or twice with doubled=True, as a search that declared it twice wrote them.
"""
import math
import random

RUNS = ["run_A1", "run_A2", "run_B1", "run_B2"]
REPS = ["Rep1", "Rep2", "Rep3", "Rep4"]
CHANNELS = {"L": "WT", "M": "KO", "H": "DKO"}
OFFSET = {"run_A1": -1.6, "run_A2": -1.5, "run_B1": -1.7, "run_B2": -1.6}   # log2 M/L and H/L


def design_rows(swap=False):
    rows = []
    for run, rep in zip(RUNS, REPS):
        chans = dict(CHANNELS)
        if swap and run in ("run_B1", "run_B2"):
            chans = {"L": "KO", "M": "WT", "H": "DKO"}
        for ch, g in chans.items():
            rows.append({"File.Name": run, "Replicate": rep, "Channel": ch, "Group": g})
    return rows


def write_design(path, swap=False, extra=()):
    import csv
    rows = design_rows(swap) + list(extra)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["File.Name", "Replicate", "Channel", "Group"])
        w.writeheader()
        w.writerows(rows)
    return path


def make_report(path, doubled=False, swap=False, label_check_run=None, seed=7):
    import pyarrow as pa
    import pyarrow.parquet as pq
    rnd = random.Random(seed)
    tag = "(SILAC)(SILAC)" if doubled else "(SILAC)"
    cols = {k: [] for k in ("Run", "Channel", "Protein.Group", "Protein.Ids", "Protein.Names",
                            "Genes", "Precursor.Id", "Precursor.Quantity", "Q.Value",
                            "PG.Q.Value", "Channel.Q.Value", "Precursor.Normalised")}

    def add(run, ch, pg, pid, qty, cq=0.001, q=0.001, pgq=0.001, name=None, gene=None):
        cols["Run"].append(run); cols["Channel"].append(ch)
        cols["Protein.Group"].append(pg); cols["Protein.Ids"].append(pg)
        cols["Protein.Names"].append(name or f"{pg}_TEST"); cols["Genes"].append(gene or f"G{pg[1:]}")
        cols["Precursor.Id"].append(pid); cols["Precursor.Quantity"].append(qty)
        cols["Q.Value"].append(q); cols["PG.Q.Value"].append(pgq); cols["Channel.Q.Value"].append(cq)
        cols["Precursor.Normalised"].append(qty)

    proteins = [f"P{i:02d}" for i in range(1, 61)]
    for pi, pg in enumerate(proteins, start=1):
        n_prec = 2 + (pi % 4)
        base = 10 + (pi % 7)
        ko = 2.0 if pi <= 10 else (-2.0 if pi <= 20 else 0.0)
        for k in range(n_prec):
            pid = f"PEPTIDE{pi:02d}{'ACDEF'[k]}K{tag}{2 + k % 2}"
            pbase = base + rnd.gauss(0, 0.5)
            for run, rep in zip(RUNS, REPS):
                if pg == "P56" and run in ("run_B1", "run_B2"):
                    continue                                  # ratio in 2 runs only
                swapped = swap and run in ("run_B1", "run_B2")
                lv = {"WT": pbase, "KO": pbase + OFFSET[run] + ko, "DKO": pbase + OFFSET[run] + ko}
                chmap = ({"L": "KO", "M": "WT", "H": "DKO"} if swapped else CHANNELS)
                for ch, g in chmap.items():
                    if pi > 50 and pi <= 55 and g != "WT":
                        add(run, ch, pg, pid, 2 ** (pbase - 6), cq=0.6)   # absent: fails channel q
                        continue
                    qty = 2 ** (lv[g] + rnd.gauss(0, 0.15))
                    cq = 0.001
                    if pi == 30 and k == 0 and ch == "H":
                        cq = 0.2                                     # one channel row fails
                    add(run, ch, pg, pid, qty, cq=cq)
    # a contaminant and an unlabelled precursor
    for run in RUNS:
        for ch in "LMH":
            add(run, ch, "Cont_P00761", f"CONTPEPK{tag}2", 2 ** 14, name="TRYP_PIG", gene="")
        add(run, "", "P01", "PEPTIDEC2", 2 ** 12)
    # a label-check run: only H loaded; light is ~2% of the heavy signal
    if label_check_run:
        for pi, pg in enumerate(proteins[:20], start=1):
            pid = f"PEPTIDE{pi:02d}AK{tag}2"
            h = 2 ** (12 + rnd.gauss(0, 0.3))
            add(label_check_run, "H", pg, pid, h)
            add(label_check_run, "L", pg, pid, h * 0.02, cq=0.4)
        add(label_check_run, "H", "P21", f"PEPTIDE21AK{tag}2", 2 ** 12)   # no light row at all
    pq.write_table(pa.table(cols), path)
    return path
