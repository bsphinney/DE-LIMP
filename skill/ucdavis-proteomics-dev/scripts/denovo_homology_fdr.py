#!/usr/bin/env python3
"""
denovo_homology_fdr.py -- FDR for de novo -> homology species ID.

Controls the two error modes that matter when calling a species from de novo
peptides BLASTed against nr: chance homology, and de novo sequencing error.
Method: NovoBoard decoy spectra (Tran NH et al., Mol Cell Proteomics 2024,
doi:10.1016/j.mcpro.2024.100849) + a linear mokapot rescoring.

WHAT THIS SHIPS, AND WHY -- every choice below was measured, not assumed. The
measurements are on the UC Davis ocelot (Leopardus pardalis) hair cohort, 9x
Exploris-480 DDA, Casanovo v5 -> DIAMOND, plus a HeLa entrapment where a Sage
database search provides ground truth.

1. SCORE WITH A LINEAR MODEL, NEVER GRADIENT BOOSTING.
   Reproduced at 1% FDR on the clean FRAC-0.8 decoy, EXTRA nr-Eukaryota:
       raw bitscore                        10,404
       E-value                             13,074
       linear  [bits+len]                  14,925
       linear  [bits+E+len+%id]            19,757
       GBT     [bits+len]                   9,554   <- worse than raw bitscore
       GBT     [bits+E+len+%id]            11,325   <- worse than E-value alone
   The non-linear model overfits the decoy in the semi-supervised loop. This is
   why Percolator/mokapot default to a linear model.

2. CALIBRATE ON E-VALUE, NEVER RAW BITSCORE. Bitscore is length-confounded, so a
   global cut starves short peptides (11 peptides at <=15 aa vs 1,661 for the
   linear model). E-value is length-normalised.

3. FIVE FEATURES, NOT THIRTEEN. On HeLa ground truth (Sage-confirmed), adding
   granular confidence/traceback features raises YIELD without raising
   CORRECTNESS:
       accept all hits        5,068 accepted   42.4% precision   2,151 correct
       5-feature  (default)   2,400            59.8%             1,436
       8-feature              2,506            55.7%             1,397
       12-feature             2,532            55.5%             1,405
   More features accept more peptides, less precisely, and recover FEWER correct
   ones. The 5-feature model is the best of the models tested.

4. CONFIDENCE IS LOAD-BEARING; HOMOLOGY ALONE FAILS. The homology-only feature
   set cannot separate at all on a standard search ("No PSMs found below the
   eval_fdr") -- too few decoy peptides get hits to train on. Every model that
   works includes mean Casanovo confidence.

5. DO NOT GATE THE SPECIES CALL ON THIS MODEL. Filtering lifts peptide precision
   42.4% -> 59.8% but moves primate attribution only 29.5% -> 32.9%, leaves
   mammal flat at 68.5%, and discards 715 genuinely-correct peptides. Species
   mis-assignment comes from CONSERVED peptides matching many organisms, which is
   an LCA problem, not a scoring problem. Use this score as a reported confidence
   tier and run LCA over the unfiltered hits.

TWO FDR DEFINITIONS, AND THEY ARE NOT INTERCHANGEABLE:
  * single raw score (bitscore / E-value): RATE-NORMALISED target-decoy, because
    the target and decoy peptide universes differ in size (312,820 vs 424,552 on
    ocelot). FDR = (decoys>=t / N_decoy) / (targets>=t / N_target).
  * mokapot score: mokapot's OWN q-values. Applying the rate-normalised formula
    to a mokapot score returns ~everything (298,889 of 299,099) because the score
    already embodies target-decoy competition and normalising again double-counts.

TIES MUST MOVE TOGETHER. Accepting hits one at a time down the ranked list gives
11,558 for bitscore instead of 10,404 (+11%); bitscore is coarse and heavily
tied. Threshold on the score VALUE so tied hits are accepted as a group.
"""
import argparse
import csv
import json
import math
import os
import sys

# Feature set. Order is fixed so a rebuilt model is comparable to a recorded one.
DEFAULT_FEATURES = ["bitscore", "log_evalue", "pep_len", "pident", "mean_conf"]

# Extended sets, kept ONLY so a user can reproduce the comparison above. They are
# not better -- see point 3. Never make one of these the default.
FEATURE_SETS = {
    "homology": ["bitscore", "log_evalue", "pep_len", "pident"],
    "default": DEFAULT_FEATURES,
    "8feat": DEFAULT_FEATURES + ["n_hi", "longest_run", "pred_count"],
    "12feat": DEFAULT_FEATURES + ["n_hi", "longest_run", "pred_count",
                                  "cterm", "nterm", "min_conf", "pep_score"],
}


def rate_normalised_accept(scores, is_target, n_target, n_decoy, fdr=0.01):
    """Tie-aware, rate-normalised target-decoy threshold.

    Returns (accepted_target_count, threshold). Higher score = better.
    Validated: reproduces 10,404 (bitscore) and 13,074 (E-value) exactly on the
    ocelot EXTRA nr-Eukaryota search against the clean FRAC-0.8 decoy.
    """
    order = sorted(range(len(scores)), key=lambda i: -scores[i])
    nt = nd = 0
    best_n = 0
    best_thr = None
    i = 0
    while i < len(order):
        s = scores[order[i]]
        # consume the whole tie group before evaluating
        while i < len(order) and scores[order[i]] == s:
            if is_target[order[i]]:
                nt += 1
            else:
                nd += 1
            i += 1
        if nt == 0:
            continue
        f = (nd / n_decoy) / (nt / n_target)
        if f <= fdr and nt > best_n:
            best_n, best_thr = nt, s
    return best_n, best_thr


def read_hits(path, layout="auto"):
    """DIAMOND tabular hits -> list of dicts, best hit per peptide.

    Handles both layouts seen in practice:
      7-col  qseqid sseqid pident length evalue bitscore staxids
      16-col qseqid sseqid pident len mm gap qs qe ss se evalue bitscore
             staxids qseq sseq btop
    """
    rows = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 7:
                continue
            if layout == "auto":
                wide = len(f) >= 13
            else:
                wide = layout == "16"
            try:
                pident = float(f[2])
                ev = float(f[10] if wide else f[4])
                bit = float(f[11] if wide else f[5])
                tax = f[12] if wide else f[6]
            except (ValueError, IndexError):
                continue
            pep = "".join(c for c in f[0].upper() if c.isalpha())
            prev = rows.get(pep)
            if prev is None or bit > prev["bitscore"]:
                rows[pep] = {"pep": pep, "subject": f[1], "pident": pident,
                             "evalue": ev, "bitscore": bit, "staxids": tax}
    for r in rows.values():
        r["pep_len"] = len(r["pep"])
        r["log_evalue"] = -math.log10(max(r["evalue"], 1e-300))
    return list(rows.values())


def read_confidence(path):
    """Casanovo per-peptide confidence features.

    Accepts either an mzTab (PSM rows with opt_ms_run[1]_aa_scores) or a TSV with
    columns peptide / peptide_score / aa_scores.
    """
    feats = {}
    counts = {}
    is_mztab = path.lower().endswith(".mztab")
    with open(path) as fh:
        hdr = None
        for line in fh:
            if is_mztab:
                if line.startswith("PSH"):
                    hdr = line.rstrip("\n").split("\t")
                    continue
                if not line.startswith("PSM") or hdr is None:
                    continue
                d = dict(zip(hdr, line.rstrip("\n").split("\t")))
                seq = d.get("sequence", "")
                aa = d.get("opt_ms_run[1]_aa_scores", "")
                score = d.get("search_engine_score[1]", "0")
            else:
                if hdr is None:
                    hdr = line.rstrip("\n").split("\t")
                    continue
                d = dict(zip(hdr, line.rstrip("\n").split("\t")))
                seq = d.get("peptide", "")
                aa = d.get("aa_scores", "")
                score = d.get("peptide_score", "0")
            pep = "".join(c for c in seq.upper() if c.isalpha())
            if not pep or not aa:
                continue
            try:
                v = [float(x) for x in aa.split(",") if x]
                sc = float(score)
            except ValueError:
                continue
            if not v:
                continue
            hi = [x >= 0.9 for x in v]
            run = best = 0
            for h in hi:
                run = run + 1 if h else 0
                best = max(best, run)
            cur = {"mean_conf": sum(v) / len(v), "min_conf": min(v),
                   "longest_run": best, "n_hi": sum(hi),
                   "cterm": v[-1], "nterm": v[0], "pep_score": sc}
            counts[pep] = counts.get(pep, 0) + 1
            prev = feats.get(pep)
            if prev is None or cur["mean_conf"] > prev["mean_conf"]:
                feats[pep] = cur
    for pep, c in counts.items():
        feats[pep]["pred_count"] = math.log1p(c)
    return feats


def attach(hits, conf):
    zero = {k: 0.0 for k in ("mean_conf", "min_conf", "longest_run", "n_hi",
                             "cterm", "nterm", "pep_score", "pred_count")}
    for h in hits:
        h.update(conf.get(h["pep"], zero))
    return hits


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target-hits", required=True, help="DIAMOND hits, real spectra")
    ap.add_argument("--decoy-hits", required=True, help="DIAMOND hits, decoy spectra")
    ap.add_argument("--target-conf", help="Casanovo mzTab/TSV for the real spectra")
    ap.add_argument("--decoy-conf", help="Casanovo mzTab/TSV for the decoy spectra")
    ap.add_argument("--n-target", type=int, required=True,
                    help="TOTAL target peptides searched (not just those with a hit)")
    ap.add_argument("--n-decoy", type=int, required=True,
                    help="TOTAL decoy peptides searched")
    ap.add_argument("--features", default="default", choices=sorted(FEATURE_SETS),
                    help="'default' is the 5-feature linear model; the others exist "
                         "only to reproduce the comparison and are NOT better")
    ap.add_argument("--fdr", type=float, default=0.01)
    ap.add_argument("--out", required=True, help="accepted peptides TSV")
    ap.add_argument("--no-model", action="store_true",
                    help="skip mokapot; threshold on E-value alone")
    a = ap.parse_args()

    t = read_hits(a.target_hits)
    d = read_hits(a.decoy_hits)
    if a.target_conf:
        t = attach(t, read_confidence(a.target_conf))
    if a.decoy_conf:
        d = attach(d, read_confidence(a.decoy_conf))
    for r in t:
        r["target"] = True
    for r in d:
        r["target"] = False
    allr = t + d
    notes = []

    feats = FEATURE_SETS[a.features]
    have_conf = any(r.get("mean_conf") for r in allr)
    use_model = not a.no_model and have_conf
    if not a.no_model and not have_conf:
        notes.append("no Casanovo confidence supplied -- falling back to E-value. "
                     "The homology-only model does not separate on its own "
                     "(measured: 'No PSMs found below the eval_fdr').")

    scored = None
    if use_model:
        try:
            import numpy as np
            for _o, _n in (("float_", "float64"), ("unicode_", "str_"), ("int_", "int64")):
                if not hasattr(np, _o) and hasattr(np, _n):
                    setattr(np, _o, getattr(np, _n))
            import pandas as pd
            import mokapot
            from sklearn.svm import LinearSVC
            df = pd.DataFrame(allr)
            df["psm_id"] = range(len(df))
            ds = mokapot.LinearPsmDataset(psms=df, target_column="target",
                                          spectrum_columns=["psm_id"],
                                          peptide_column="pep", feature_columns=feats)
            # LinearSVC on purpose. Gradient boosting overfits the decoy here and
            # scores WORSE than raw bitscore -- see the header.
            res, _ = mokapot.brew(ds, mokapot.Model(LinearSVC(dual="auto")),
                                  test_fdr=a.fdr)
            pep = res.peptides if hasattr(res, "peptides") else \
                res.confidence_estimates["peptides"]
            qcol = [c for c in pep.columns if "q-value" in c.lower()][0]
            keep = set(pep.loc[pep[qcol] <= a.fdr, "pep"])
            scored = [r for r in t if r["pep"] in keep]
            notes.append(f"mokapot linear, features={feats}; mokapot q-values "
                         f"(NOT rate-normalised -- see header)")
        except ImportError as e:
            notes.append(f"mokapot/sklearn unavailable ({e}); falling back to E-value")
            use_model = False
        except Exception as e:
            notes.append(f"mokapot failed ({type(e).__name__}: {e}); "
                         f"falling back to E-value")
            use_model = False

    if scored is None:
        s = [r["log_evalue"] for r in allr]
        tg = [r["target"] for r in allr]
        n, thr = rate_normalised_accept(s, tg, a.n_target, a.n_decoy, a.fdr)
        scored = [r for r in t if thr is not None and r["log_evalue"] >= thr]
        notes.append("scored on E-value with the rate-normalised, tie-aware "
                     "target-decoy threshold")

    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        cols = ["pep", "subject", "staxids", "pident", "evalue", "bitscore",
                "pep_len", "mean_conf"]
        w.writerow(cols)
        for r in sorted(scored, key=lambda r: -r["bitscore"]):
            w.writerow([r.get(c, "") for c in cols])

    print(json.dumps({
        "accepted": len(scored),
        "target_peptides_with_hit": len(t),
        "decoy_peptides_with_hit": len(d),
        "n_target": a.n_target, "n_decoy": a.n_decoy,
        "features": feats if use_model else None,
        "model": "mokapot linear (LinearSVC)" if use_model else "E-value threshold",
        "fdr": a.fdr,
        "out": os.path.abspath(a.out),
        "notes": notes,
        "species_call_guidance":
            "Do NOT gate the species call on this set. Run LCA over the UNFILTERED "
            "hits; filtering lifts peptide precision 42.4%->59.8% but moves primate "
            "attribution only 29.5%->32.9% and discards 715 correct peptides. "
            "Species error is conserved peptides -- an LCA problem, not a scoring one.",
    }, indent=2))


if __name__ == "__main__":
    main()
