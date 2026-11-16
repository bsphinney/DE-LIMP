#!/usr/bin/env python3
"""Writes the hand-made mini screen_data folder beside this file (manifest.json, index.tsv, parts/).

Same shape as make_screen_data.R's output (schema screen_data/1), small enough to read:
3 groups (A, B, C) x 2 samples; contrasts B-A and C-A; four proteins in one part:
  P10001  GENEA  plain; shares, detection and weights present; peptides shared with P10003 and P10009
  P10002  GENEX  } two protein groups with the same gene name (an ambiguous query)
  P10003  GENEX  }
  P10004  GENED  not tested in B-A (measured in neither A nor B); shares withheld; no weights
The sha256 values are computed from the files written here. Run `python3 make.py` after any edit.
"""
import gzip
import hashlib
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLES = ["a1", "a2", "b1", "b2", "c1", "c2"]
GROUPS = ["A", "A", "B", "B", "C", "C"]
# The delivered DE table is not part of the fixture; tests write this text where the manifest says it is.
DE_TABLE_TEXT = "Protein.Group,logFC\nP10001,0.5\n"
DE_TABLE_SHA = hashlib.sha256(DE_TABLE_TEXT.encode()).hexdigest()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        h.update(fh.read())
    return h.hexdigest()


def de_entry(tested, lfc, det, p=0.01, adj=0.04):
    if not tested:
        return {"tested": False, "logFC": None, "CI.L": None, "CI.R": None, "t": None, "P.Value": None,
                "adj.P.Val": None, "df_total": None, "evidence": "not measured in either group",
                "detected": det}
    return {"tested": True, "logFC": lfc, "CI.L": lfc - 0.3, "CI.R": lfc + 0.3, "t": lfc / 0.15,
            "P.Value": p, "adj.P.Val": adj, "df_total": 9.5, "evidence": "measured in both",
            "detected": det}


def rows_for(pids, protein_ids, used=True):
    out = []
    for pid, ids in zip(pids, protein_ids):
        for k, s in enumerate(SAMPLES):
            out.append({"Run": s, "sample": s, "Precursor.Id": pid, "Protein.Ids": ids, "Proteotypic": 1,
                        "Precursor.Quantity": 1000000 + k, "Normalisation.Factor": 1.0,
                        "Precursor.Normalised": 1000000 + k, "Quantity.Quality": 0.9,
                        "Empirical.Quality": 0.8, "Ms1.Area": 2000000, "RT": 12.5 + k / 10,
                        "Q.Value": 0.001, "Lib.Q.Value": 0.001, "PG.MaxLFQ": 3000000,
                        "PG.MaxLFQ.Quality": 0.9, "failed": [], "contaminant": False, "used": used})
    return out


def record(pg, genes, names, pids, protein_ids, de, shares, detection, weights, base=20.0):
    n = len(SAMPLES)
    return {
        "pg": pg, "genes": genes, "names": names, "n_precursors": len(pids), "sigma": 0.4, "de": de,
        "values": {"protein": [base + 0.1 * i for i in range(n)], "se": [0.1] * n, "n_obs": [2] * n},
        "precursors": [{"id": p, "log2": [None if (j == 1 and i == 5) else round(base + 0.2 * i, 4)
                                          for i in range(n)]} for j, p in enumerate(pids)],
        "rows": rows_for(pids, protein_ids), "shares": shares, "detection": detection, "weights": weights}


def build():
    n = len(SAMPLES)
    det_ab = {"A": "2/2", "B": "2/2", "C": "0/2"}
    shares_ok = {"own": [[0.4] * n, [0.5 if i != 5 else None for i in range(n)]], "other": [0.05] * n,
                 "total": [0.95] * n, "refit_max_diff": 1.0e-9}
    det_ok = {"pf_fit": [[0.9] * n, [0.8] * n], "pf_obs": [[0.9] * n, [0.8 if i != 5 else None for i in range(n)]]}
    # 2 samples per group with equal weights: each carries half of its group's weighted mean
    weights_ok = {"w": [100.0] * n, "share_in_group": [0.5] * n, "df_residual": 3, "df_prior": 8.5}
    p1 = ["AEEATEAQEVVEATPEGAC(UniMod:4)TEPR3", "LGSEQGGGK2"]
    recs = {
        "P10001": record("P10001", "GENEA", "AAA_HUMAN", p1, ["P10001;P10003", "P10001;P10009"],
                         {"B-A": de_entry(True, 0.5, {"A": "2/2", "B": "2/2"}, 0.001, 0.03),
                          "C-A": de_entry(True, -0.2, {"A": "2/2", "C": "2/2"}, 0.02, 0.04)},
                         shares_ok, det_ok, weights_ok),
        "P10002": record("P10002", "GENEX", "XXA_HUMAN", ["PEPTIDEK2"], ["P10002"],
                         {"B-A": de_entry(True, 0.1, {"A": "2/2", "B": "2/2"}, 0.01, 0.05),
                          "C-A": de_entry(True, 0.0, {"A": "2/2", "C": "1/2"}, 0.5, 0.9)},
                         dict(shares_ok, own=[[0.4] * n], other=[0.0] * n, total=[0.4] * n),
                         {"pf_fit": [[0.9] * n], "pf_obs": [[0.9] * n]}, weights_ok),
        "P10003": record("P10003", "GENEX", "XXB_HUMAN", ["OTHERPEPK2"], ["P10003"],
                         {"B-A": de_entry(True, -0.1, {"A": "2/2", "B": "2/2"}, 0.004, 0.04),
                          "C-A": de_entry(True, 0.0, {"A": "2/2", "C": "1/2"}, 0.5, 0.9)},
                         dict(shares_ok, own=[[0.4] * n], other=[0.0] * n, total=[0.4] * n),
                         {"pf_fit": [[0.9] * n], "pf_obs": [[0.9] * n]}, weights_ok),
        "P10004": record("P10004", "GENED", "DDD_HUMAN", ["ONLYCK2"], ["P10004"],
                         {"B-A": de_entry(False, None, {"A": "0/2", "B": "0/2"}),
                          "C-A": de_entry(True, 1.2, {"A": "0/2", "C": "2/2"}, 0.0005, 0.01)},
                         {"withheld": "refit differs from the delivered value by 5.15e-05"}, None, None,
                         base=18.0),
    }
    recs["P10004"]["weights"] = None
    return recs


def main():
    recs = build()
    os.makedirs(os.path.join(HERE, "parts"), exist_ok=True)
    part = os.path.join(HERE, "parts", "part-0001.json.gz")
    with open(part, "wb") as raw:
        with gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
            gz.write(json.dumps(recs, separators=(",", ":"), allow_nan=False).encode())
    with open(os.path.join(HERE, "index.tsv"), "w", newline="") as fh:
        cons = ["B-A", "C-A"]
        cols = ["%s__%s" % (k, c) for c in cons for k in ("P.Value", "adj.P.Val", "logFC")]
        fh.write("\t".join(["Protein.Group", "Genes", "Protein.Names", "part", "n_precursors"] + cols) + "\n")
        for pg, r in recs.items():
            vals = []
            for c in cons:
                d = r["de"][c]
                vals += [("" if not d["tested"] else repr(d[k])) for k in ("P.Value", "adj.P.Val", "logFC")]
            fh.write("\t".join([pg, r["genes"], r["names"], "part-0001.json.gz", str(len(r["precursors"]))]
                                + vals) + "\n")
    manifest = {
        "schema": "screen_data/1", "written_by": "make_screen_data.R (hand-made fixture)",
        "written_at": "2026-10-28T00:00:00Z", "supported": True, "unsupported_reason": None, "partial": False,
        "analysis": {"pipeline_id": "dpc", "display_label": "DPC-Quant + limma (limpa)",
                     "normalisation": "DIA-NN's cross-run normalisation (Precursor.Normalised); limpa applies none of its own",
                     "quantities": "normalised", "q_columns": ["Q.Value", "Lib.Q.Value"], "q_cutoffs": [0.01, 0.01],
                     "filters_applied": ["Q.Value <= 0.01", "Lib.Q.Value <= 0.01"],
                     "rollup_method": "DPC-Quant (Detection Probability Curve quantification, dpcCN) of Precursor.Normalised",
                     "de_engine": "limpa::dpcDE (voomaLmFitWithImputation) -> contrasts.fit -> eBayes",
                     "significance_rule": "adj.P.Val < adjp (BH); no fold-change filter"},
        "samples": SAMPLES, "groups": GROUPS,
        "contrasts": [
            {"name": "B-A", "numerator": "B", "denominator": "A", "table": "DE_dpc_B.A.csv",
             "n_tested": 3, "n_sig": 2, "adjp": 0.05},
            {"name": "C-A", "numerator": "C", "denominator": "A", "table": "DE_dpc_C.A.csv",
             "n_tested": 4, "n_sig": 2, "adjp": 0.05}],
        "dpc": [-12.5, 0.7], "prior": {"prior.mean": 20.0, "prior.sd": 1.8, "prior.logFC": 1},
        "n_precursors": 5, "run_stats": None,
        "sources": {"tables/DE_dpc_B.A.csv": DE_TABLE_SHA},
        "parts": {"parts/part-0001.json.gz": sha256(part),
                  "index.tsv": sha256(os.path.join(HERE, "index.tsv"))}, "index": "index.tsv",
        "skipped": [{"block": "run_stats", "reason": "search/report.stats.tsv does not exist"}]
                   + [{"block": b, "reason": "added in a later part (design section 12)"}
                      for b in ("xic", "library_uniqueness", "check_routes", "fran", "report_text")],
        "shares": {"method": "derivative", "tolerance": 1e-6, "settle_tolerance": 1e-5,
                   "withheld": [{"pg": "P10004", "reason": "refit differs from the delivered value by 5.15e-05"}],
                   "nudge_check": {"proteins": ["P10001"], "cells_checked": 11, "cells_found": 11, "step": 0.001,
                                   "max_diff": 4e-10, "tolerance": 1e-4, "passed": True, "elapsed_s": 0.1},
                   "limpa_version": "1.4.0"},
        "detection": {"withheld": [{"pg": "P10004", "reason": "refit differs from the delivered value by 5.15e-05"}]},
        "elapsed_s": {"total": 1.0, "read": 0.1, "shares": 0.1, "records": 0.1, "write": 0.1},
        "compute": {"de": "zen4", "here": "zen4"}}
    with open(os.path.join(HERE, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=2)
        fh.write("\n")


if __name__ == "__main__":
    main()
