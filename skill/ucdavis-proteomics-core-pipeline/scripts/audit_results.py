#!/usr/bin/env python3
"""
audit_results.py  --  Sanity-check the analysis for common proteomics mistakes a
new user can make, BEFORE they over-interpret the results. Deterministic checks
grounded in the real data (never fabricated). Emits AUDIT.md + audit.json with a
PASS / WARN / FAIL per check.

The orchestrator runs this after DE and **surfaces every WARN/FAIL to the user**;
FAILs should stop interpretation until resolved (e.g. a group with no replicate).
The findings also become the report's "Audit & caveats" section.

Checks:
  replication       groups with <2 (FAIL) or <3 (WARN) replicates
  group_balance     very unequal group sizes (WARN)
  confounding       a covariate (Batch) perfectly confounded with Group (FAIL)
  acquisition_mix   DIA and DDA mixed in one analysis (FAIL)
  instrument_mix    >1 instrument model in one analysis — batch effect (WARN)
  id_depth          suspiciously few proteins quantified (WARN)
  missingness       very high fraction of missing values (WARN)
  contamination     keratin/trypsin contaminants among the proteins (WARN)
  target_contaminants
                    quantified proteins whose sequence is ALSO a common-contaminant
                    entry (fetch_fasta.py removed those entries so the protein is
                    quantified under its own accession): possible contamination, KEPT
                    in quantification -- reported, never excluded (WARN)
  contaminant_overlap
                    real proteins present only as identical Cont_ entries in a
                    database used as-is -- excluded from quant, missing here (WARN)
  de_signal         0 significant (WARN: underpowered) or >50% significant
                    (WARN: likely batch/normalization/confounding artefact)

Usage:
  python3 audit_results.py --out AUDIT.md \
      --conditions input/conditions.csv --de-dir output/tables \
      [--acquisition-json acq.json] [--adjp 0.05] [--logfc 1] \
      [--min-proteins 500] [--max-missing 0.5] \
      [--fasta-meta search.fasta.meta.json]   # default: ./search.fasta.meta.json if present
"""
import sys, os, csv, json, glob, argparse
from collections import Counter, defaultdict

# The list's ONE definition is the FASTA sidecar; the wording for lost proteins lives there too.
from fetch_fasta import target_contaminants, seen_only_as_cont, lost_to_contaminants_message

CONTAMINANT_PATTERNS = ("KRT", "KRTAP",            # keratins (skin/hair)
                        "TRYP", "PRSS1", "TRY1",   # trypsin (digestion)
                        "CASEIN", "CSN1", "CSN2",  # casein (milk)
                        )
CONTAMINANT_WORDS = ("keratin", "trypsin", "casein", "contaminant", "cRAP")


def add(findings, cid, status, message, detail=None):
    findings.append({"check": cid, "status": status, "message": message, "detail": detail or {}})


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def audit_conditions(findings, rows, adjp):
    groups = Counter(r.get("Group", "").strip() for r in rows if r.get("Group", "").strip())
    if not groups:
        add(findings, "replication", "FAIL", "No groups assigned in conditions.csv.")
        return groups
    singletons = [g for g, n in groups.items() if n < 2]
    small = [g for g, n in groups.items() if n == 2]
    if singletons:
        add(findings, "replication", "FAIL",
            f"Group(s) with <2 replicates have no within-group variance — differential statistics are not valid: {singletons}.",
            {"group_sizes": dict(groups)})
    elif small:
        add(findings, "replication", "WARN",
            f"Group(s) with only 2 replicates: {small}. Usable but low power; 3+ is recommended.",
            {"group_sizes": dict(groups)})
    else:
        add(findings, "replication", "PASS", f"All groups have ≥3 replicates.", {"group_sizes": dict(groups)})
    # balance
    if len(groups) >= 2:
        hi, lo = max(groups.values()), min(groups.values())
        if lo and hi / lo >= 3:
            add(findings, "group_balance", "WARN",
                f"Group sizes are very unequal ({dict(groups)}). Large imbalance reduces power and can bias results.")
        else:
            add(findings, "group_balance", "PASS", "Group sizes are reasonably balanced.")
    # confounding: Batch (or covariate) perfectly nested within Group
    for cov in ("Batch", "Covariate1", "Covariate2"):
        if rows and cov in rows[0]:
            pairs = defaultdict(set)
            for r in rows:
                g, c = r.get("Group", "").strip(), r.get(cov, "").strip()
                if g and c:
                    pairs[g].add(c)
            # confounded if each group maps to exactly one distinct cov value AND covs differ across groups
            covsets = [next(iter(v)) for v in pairs.values() if len(v) == 1]
            if pairs and all(len(v) == 1 for v in pairs.values()) and len(set(covsets)) == len(pairs) and len(pairs) > 1:
                add(findings, "confounding", "FAIL",
                    f"'{cov}' is perfectly confounded with Group — the batch/covariate effect cannot be "
                    f"separated from the biological effect. {dict((g, sorted(v)) for g, v in pairs.items())}")
            else:
                add(findings, "confounding", "PASS", f"'{cov}' is not confounded with Group.")
    return groups


def audit_acquisition(findings, acq_json):
    data = None
    try:
        data = json.load(open(acq_json))
    except Exception:
        return
    overall = data.get("overall")
    if overall == "mixed":
        add(findings, "acquisition_mix", "FAIL",
            "DIA and DDA files are mixed in one analysis. They must be searched separately "
            "(DIA→DIA-NN, DDA→Sage); mixing them produces invalid quantification.",
            {"files": [(f.get("file"), f.get("acquisition")) for f in data.get("files", [])]})
    elif overall in ("DIA", "DDA"):
        add(findings, "acquisition_mix", "PASS", f"All files are {overall}.")
    instr = data.get("instruments_seen") or []
    if len(instr) > 1:
        add(findings, "instrument_mix", "WARN",
            f"More than one instrument model in this analysis ({instr}). This introduces a batch "
            "effect — account for it (add a Batch covariate) or analyze per instrument.")
    elif instr:
        add(findings, "instrument_mix", "PASS", f"Single instrument: {instr[0]}.")


def audit_matrix(findings, em_path, min_proteins, max_missing, keratin_sample=False):
    rows = read_csv(em_path)
    if not rows:
        add(findings, "id_depth", "WARN", "Expression matrix is empty.")
        return
    idcols = [c for c in ("Protein.Group", "Genes", "Protein.Names") if c in rows[0]]
    sample_cols = [c for c in rows[0].keys() if c not in idcols]
    n_prot = len(rows)
    if n_prot < min_proteins:
        add(findings, "id_depth", "WARN",
            f"Only {n_prot} proteins quantified — lower than typical (<{min_proteins}). "
            "Check the search parameters and that the FASTA matches the sample's organism.")
    else:
        add(findings, "id_depth", "PASS", f"{n_prot} proteins quantified.")
    # missingness
    total = miss = 0
    for r in rows:
        for c in sample_cols:
            total += 1
            v = (r.get(c) or "").strip()
            if v in ("", "NA", "NaN", "nan"):
                miss += 1
    frac = miss / total if total else 0
    if frac > max_missing:
        add(findings, "missingness", "WARN",
            f"{frac*100:.0f}% of values are missing (> {max_missing*100:.0f}%). High missingness "
            "weakens quantification; check sample quality and whether missingness differs by group.")
    else:
        add(findings, "missingness", "PASS", f"{frac*100:.0f}% missing values.")
    # contamination. Keratin is a contaminant for MOST matrices, but is the ANALYTE
    # for keratin samples (nail / hair / wool / skin / feather) — never flag it there.
    # (see references/anomaly-checks.md "Keratin-matrix samples").
    pats = tuple(p for p in CONTAMINANT_PATTERNS if not (keratin_sample and p in ("KRT", "KRTAP")))
    words = tuple(w for w in CONTAMINANT_WORDS if not (keratin_sample and w == "keratin"))
    glab = "Genes" if "Genes" in rows[0] else ("Protein.Names" if "Protein.Names" in rows[0] else None)
    if glab:
        names = [(r.get(glab) or "").upper() for r in rows]
        hits = [n for n in names if any(n.startswith(p) or p in n for p in pats)
                or any(w.upper() in n for w in words)]
        label = "trypsin/casein" if keratin_sample else "keratins/trypsin/casein"
        if hits:
            frac_c = len(hits) / n_prot
            status = "WARN" if frac_c >= 0.02 or len(hits) >= 10 else "PASS"
            add(findings, "contamination", status,
                f"{len(hits)} likely contaminant protein(s) detected ({label}). "
                + ("A high contaminant load can distort normalization and quant — consider filtering them."
                   if status == "WARN" else "Low level; usually fine.")
                + (" [keratin-matrix sample: keratins treated as the analyte, not flagged]" if keratin_sample else ""),
                {"examples": sorted(set(hits))[:10]})
        elif keratin_sample:
            add(findings, "contamination", "PASS",
                "Keratin-matrix sample — keratins are the analyte (not flagged); no trypsin/casein issue.")


def _tokens(cell):
    return {t.strip().upper() for t in (cell or "").split(";") if t.strip()}


def audit_target_contaminants(findings, meta_path, em_path, de_dir, adjp, keratin_sample=False):
    """Proteins that are both a <organism> protein and a common contaminant.

    Why (HIVE, 2026-09-24): contaminant entries identical to target proteins made DIA-NN
    report ACTB, EEF1A1 and KRT8 only as Cont_ groups, excluded from quant. fetch_fasta.py
    now removes those entries, so the proteins are quantified under their own accession --
    which also means skin/hair keratins (usually handling contamination) count toward
    normalisation. That trade-off is surfaced here, never hidden and never "fixed" by
    excluding them: the list comes from the FASTA sidecar, not from a copy kept here.
    """
    if not meta_path:
        add(findings, "target_contaminants", "INFO",
            "Not assessed: no <fasta>.meta.json from fetch_fasta.py (pass --fasta-meta), so "
            "proteins that are also common-contaminant sequences cannot be identified.")
        return
    try:
        with open(meta_path) as fh:
            # Keratin-matrix sample: keratin is the analyte, so never flag it here.
            tc = target_contaminants(json.load(fh), keratin_sample)
    except (OSError, ValueError) as e:
        add(findings, "target_contaminants", "INFO",
            f"Not assessed: could not read {meta_path} ({e}).")
        return
    org = tc["organism"] or "target-organism"

    em_rows = read_csv(em_path) if em_path and os.path.exists(em_path) else []

    # Real proteins that exist only as Cont_ entries: a database used as-is, one built with
    # --keep-target-contaminants, or one built BEFORE the overlap check (legacy sidecar --
    # the databases that lost ACTB/EEF1A1/KRT8). The matrix makes it concrete where it can.
    kept = tc["kept_as_contaminant"]
    seen = seen_only_as_cont(kept, [(_tokens(r.get("Protein.Group")), r.get("Genes") or "?")
                                    for r in em_rows])
    msg = lost_to_contaminants_message(tc, seen)
    if msg:
        add(findings, "contaminant_overlap", "WARN", msg,
            {"fasta_meta": meta_path, "legacy_database": bool(tc.get("legacy_note")),
             "genes": sorted({r.get("gene") or r.get("target_acc") or "?" for r in kept}),
             "seen_only_as_cont": seen})

    if not tc["dropped"]:
        return
    genes, accs = tc["genes"], tc["accessions"]

    def match(r):
        # Gene OR accession: an NCBI database has no GN= in its headers, so its records
        # carry no gene and only the accession can match.
        return bool(_tokens(r.get("Genes")) & genes) or bool(_tokens(r.get("Protein.Group")) & accs)

    def label(r):
        return (r.get("Genes") or r.get("Protein.Group") or "?").split(";")[0]

    hits = sorted({label(r) for r in em_rows if match(r)})
    sig = {}
    for f in sorted(glob.glob(os.path.join(de_dir or "", "DE_*.csv"))):
        s = set()
        for r in read_csv(f):
            try:
                if float(r.get("adj.P.Val", "nan")) < adjp and match(r):
                    s.add(label(r))
            except ValueError:
                continue
        if s:
            sig[os.path.basename(f)] = sorted(s)
    if not hits and not sig:
        add(findings, "target_contaminants", "PASS",
            f"None of the {len(tc['dropped'])} {org} proteins that are also common-contaminant "
            f"sequences were quantified.", {"fasta_meta": meta_path})
        return
    msg = (f"{len(hits)} quantified {org} protein group(s) are also common-contaminant sequences "
           f"({', '.join(hits[:12])}{', ...' if len(hits) > 12 else ''}): possible contamination, "
           f"KEPT in quantification and normalisation. Skin/hair keratins (KRT1/2/9/10, KRTAPs, "
           f"FLG, HRNR) are usually handling contamination; ACTB/EEF1A1/tubulins -- and the "
           f"simple-epithelial keratins KRT8/18/19/7 in epithelial cells such as HeLa -- are "
           f"usually endogenous.")
    if sig:
        msg += (" Significant in DE -- check these are not contamination before interpreting "
                "them: " + "; ".join(f"{k}: {', '.join(v[:10])}" for k, v in sig.items()) + ".")
    add(findings, "target_contaminants", "WARN", msg,
        {"fasta_meta": meta_path, "quantified": hits, "significant": sig})


def audit_de(findings, de_dir, adjp, logfc):
    for f in sorted(glob.glob(os.path.join(de_dir, "DE_*.csv"))):
        rows = read_csv(f)
        ct = os.path.basename(f)
        n = sig = beyond = 0
        for r in rows:
            try:
                a = float(r.get("adj.P.Val", "nan")); l = float(r.get("logFC", "nan"))
            except ValueError:
                continue
            if a != a:  # nan
                continue
            n += 1
            # Significance is the adjusted p-value alone -- the same rule run_de.R and the
            # figures use. |logFC| is only counted descriptively, for the report.
            if a < adjp:
                sig += 1
                if abs(l) >= logfc:
                    beyond += 1
        if n == 0:
            continue
        frac = sig / n
        if sig == 0:
            add(findings, "de_signal", "WARN",
                f"{ct}: 0 proteins significant (adj.P<{adjp}). The experiment may be "
                "underpowered, the effect small, or the groups mislabeled.")
        elif frac > 0.5:
            add(findings, "de_signal", "WARN",
                f"{ct}: {sig}/{n} ({frac*100:.0f}%) proteins significant — implausibly high. This usually "
                "means a batch effect, a normalization problem, or confounded/mislabeled groups, not real biology.")
        else:
            add(findings, "de_signal", "PASS",
                f"{ct}: {sig}/{n} proteins significant at adj.P<{adjp} ({frac*100:.0f}%); "
                f"{beyond} of those are also ≥{2**logfc:.3g}-fold.")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="AUDIT.md")
    ap.add_argument("--conditions")
    ap.add_argument("--de-dir")
    ap.add_argument("--acquisition-json")
    ap.add_argument("--adjp", type=float, default=0.05)
    ap.add_argument("--logfc", type=float, default=1.0)
    ap.add_argument("--min-proteins", type=int, default=500)
    ap.add_argument("--max-missing", type=float, default=0.5)
    ap.add_argument("--keratin-sample", action="store_true",
                    help="sample IS keratin (nail/hair/wool/skin/feather) — keratin is the analyte, "
                         "not a contaminant; do not flag KRT/KRTAP/keratin")
    ap.add_argument("--fasta-meta",
                    help="fetch_fasta.py's <fasta>.meta.json -- the list of proteins that are also "
                         "common-contaminant sequences. Default: ./search.fasta.meta.json if present")
    a = ap.parse_args()
    fasta_meta = a.fasta_meta or ("search.fasta.meta.json"
                                  if os.path.exists("search.fasta.meta.json") else None)

    findings = []
    if a.conditions and os.path.exists(a.conditions):
        audit_conditions(findings, read_csv(a.conditions), a.adjp)
    if a.acquisition_json and os.path.exists(a.acquisition_json):
        audit_acquisition(findings, a.acquisition_json)
    em = os.path.join(a.de_dir, "Expression_Matrix.csv") if a.de_dir else None
    if em and os.path.exists(em):
        audit_matrix(findings, em, a.min_proteins, a.max_missing, a.keratin_sample)
    if a.de_dir and os.path.isdir(a.de_dir):
        audit_target_contaminants(findings, fasta_meta, em, a.de_dir, a.adjp, a.keratin_sample)
        audit_de(findings, a.de_dir, a.adjp, a.logfc)

    n_fail = sum(1 for f in findings if f["status"] == "FAIL")
    n_warn = sum(1 for f in findings if f["status"] == "WARN")
    overall = "FAIL" if n_fail else ("WARN" if n_warn else "PASS")

    # INFO = a check that could not run; it never changes the overall status.
    icon = {"PASS": "✅", "WARN": "⚠️", "FAIL": "⛔", "INFO": "ℹ️"}
    lines = ["# Results audit — common proteomics pitfalls", "",
             f"**Overall: {icon[overall]} {overall}** — {n_fail} blocking, {n_warn} warning(s).", ""]
    if n_fail:
        lines += ["> ⛔ Resolve the blocking issues before trusting the differential-expression results.", ""]
    for f in findings:
        lines.append(f"- {icon[f['status']]} **{f['check']}** — {f['message']}")
    lines.append("")
    with open(a.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    with open(os.path.splitext(a.out)[0] + ".json"
              if a.out.endswith(".md") else a.out + ".json", "w") as fh:
        json.dump({"overall": overall, "n_fail": n_fail, "n_warn": n_warn, "findings": findings}, fh, indent=2)

    print(json.dumps({"overall": overall, "n_fail": n_fail, "n_warn": n_warn,
                      "report": os.path.abspath(a.out),
                      "blocking": [f["message"] for f in findings if f["status"] == "FAIL"],
                      "warnings": [f["message"] for f in findings if f["status"] == "WARN"]}, indent=2))


if __name__ == "__main__":
    main()
