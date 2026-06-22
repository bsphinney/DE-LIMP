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
  de_signal         0 significant (WARN: underpowered) or >50% significant
                    (WARN: likely batch/normalization/confounding artefact)

Usage:
  python3 audit_results.py --out AUDIT.md \
      --conditions input/conditions.csv --de-dir output/tables \
      [--acquisition-json acq.json] [--adjp 0.05] [--logfc 1] \
      [--min-proteins 500] [--max-missing 0.5]
"""
import sys, os, csv, json, glob, argparse
from collections import Counter, defaultdict

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


def audit_matrix(findings, em_path, min_proteins, max_missing):
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
    # contamination
    glab = "Genes" if "Genes" in rows[0] else ("Protein.Names" if "Protein.Names" in rows[0] else None)
    if glab:
        names = [(r.get(glab) or "").upper() for r in rows]
        hits = [n for n in names if any(n.startswith(p) or p in n for p in CONTAMINANT_PATTERNS)
                or any(w.upper() in n for w in CONTAMINANT_WORDS)]
        if hits:
            frac_c = len(hits) / n_prot
            status = "WARN" if frac_c >= 0.02 or len(hits) >= 10 else "PASS"
            add(findings, "contamination", status,
                f"{len(hits)} likely contaminant protein(s) detected (keratins/trypsin/casein). "
                + ("A high contaminant load can distort normalization and quant — consider filtering them."
                   if status == "WARN" else "Low level; usually fine."),
                {"examples": sorted(set(hits))[:10]})


def audit_de(findings, de_dir, adjp, logfc):
    for f in sorted(glob.glob(os.path.join(de_dir, "DE_*.csv"))):
        rows = read_csv(f)
        ct = os.path.basename(f)
        n = sig = 0
        for r in rows:
            try:
                a = float(r.get("adj.P.Val", "nan")); l = float(r.get("logFC", "nan"))
            except ValueError:
                continue
            if a != a:  # nan
                continue
            n += 1
            if a < adjp and abs(l) >= logfc:
                sig += 1
        if n == 0:
            continue
        frac = sig / n
        if sig == 0:
            add(findings, "de_signal", "WARN",
                f"{ct}: 0 proteins significant (adj.P<{adjp}, |logFC|≥{logfc}). The experiment may be "
                "underpowered, the effect small, or the groups mislabeled.")
        elif frac > 0.5:
            add(findings, "de_signal", "WARN",
                f"{ct}: {sig}/{n} ({frac*100:.0f}%) proteins significant — implausibly high. This usually "
                "means a batch effect, a normalization problem, or confounded/mislabeled groups, not real biology.")
        else:
            add(findings, "de_signal", "PASS", f"{ct}: {sig}/{n} proteins significant ({frac*100:.0f}%).")


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
    a = ap.parse_args()

    findings = []
    if a.conditions and os.path.exists(a.conditions):
        audit_conditions(findings, read_csv(a.conditions), a.adjp)
    if a.acquisition_json and os.path.exists(a.acquisition_json):
        audit_acquisition(findings, a.acquisition_json)
    if a.de_dir and os.path.exists(os.path.join(a.de_dir, "Expression_Matrix.csv")):
        audit_matrix(findings, os.path.join(a.de_dir, "Expression_Matrix.csv"), a.min_proteins, a.max_missing)
    if a.de_dir and os.path.isdir(a.de_dir):
        audit_de(findings, a.de_dir, a.adjp, a.logfc)

    n_fail = sum(1 for f in findings if f["status"] == "FAIL")
    n_warn = sum(1 for f in findings if f["status"] == "WARN")
    overall = "FAIL" if n_fail else ("WARN" if n_warn else "PASS")

    icon = {"PASS": "✅", "WARN": "⚠️", "FAIL": "⛔"}
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
