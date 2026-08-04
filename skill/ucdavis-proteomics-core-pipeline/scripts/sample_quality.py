#!/usr/bin/env python3
"""
sample_quality.py -- biological sample-quality / contamination diagnostics.

Goes beyond audit_results.py's lab-contaminant check (keratin/trypsin) to the
sample-level biological contaminations that repeatedly produce WRONG-but-plausible
DE in real core-facility work (see references/sample-quality.md, distilled from the
UC Davis Proteomics Core analysis log):

  * HEMOLYSIS      red-cell lysis in plasma/serum (hemoglobin, carbonic anhydrase,
                   catalase, peroxiredoxin, spectrin ...). Drove artefactual
                   plasma DE (Taha dog P1-vs-P2; Pcal GZ1).
  * SKELETAL_MUSCLE muscle debris in a non-muscle tissue biopsy (myosins, troponins,
                   CK-M, myoglobin ...). In the Vining cow-liver study, uniform
                   muscle contamination of one breed's biopsies produced 3,045 "DE"
                   proteins that were contamination, NOT biology.
  * EPIDERMIS      skin/hair squames (epidermal keratins + FLG/LOR/IVL).

For each panel it computes a per-sample abundance score (z across samples) and --
crucially -- checks whether that score is CONFOUNDED WITH GROUP. When a
contamination panel separates the groups with little/no overlap, DE between those
groups may just be the contamination gradient; protein-level marker removal does
NOT fix a confounded contrast (proven in the Vining study -- dropping muscle
markers INCREASED the DE count). The fix is at the sample/design level.

Also flags the limpa/DPC-Quant "complete matrix" depth trap: a DPC-Quant expression
matrix has ~no NAs (the detection model fills every cell), so counting non-empty
cells is a CONSTANT, not per-sample depth -- real depth needs detected-precursor
counts. Pass --report report.parquet for true per-run detected protein-group depth.

Usage:
  python3 sample_quality.py --matrix Expression_Matrix.csv [--conditions conditions.csv]
      [--report report.parquet] [--out SAMPLE_QUALITY.md] [--z 1.5]

Reads/writes plain files; stdlib only (pyarrow optional, just for --report).
"""
import sys, os, csv, json, math, argparse, re

# Curated gene-symbol panels (case-insensitive; matched against the matrix's gene /
# protein-name column). Deliberately specific markers -- avoid ubiquitous glycolytic
# enzymes. Extend per tissue as needed; document additions in references/sample-quality.md.
PANELS = {
    "HEMOLYSIS": ["HBA1", "HBA2", "HBA", "HBB", "HBD", "CA1", "CA2", "CAT", "PRDX2",
                  "BLVRB", "SPTA1", "SPTB", "ANK1", "SLC4A1", "EPB42", "PKLR", "BPGM",
                  "ALAS2", "AHSP", "CATALASE"],
    "SKELETAL_MUSCLE": ["MYH1", "MYH2", "MYH7", "MYBPC1", "MYBPC2", "ACTN2", "ACTN3",
                        "CKM", "MB", "TNNI1", "TNNI2", "TNNT3", "TNNC2", "PYGM", "ENO3",
                        "TPM1", "TPM2", "MYOM1", "MYOM2", "MYOM3", "DES", "MYL1", "MYL2",
                        "ATP2A1", "NEB", "TTN", "CASQ1", "SLN"],
    "EPIDERMIS": ["KRT1", "KRT2", "KRT9", "KRT10", "KRT5", "KRT14", "KRT16", "KRT6A",
                  "FLG", "LOR", "IVL", "DSP", "JUP", "SBSN"],
}


def read_matrix(path):
    """Return (samples, rows) where each row = {'genes': set, 'vals': {sample: float|None}}."""
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        # gene/id column: prefer an explicit gene column, else first text column
        gi = next((i for i, h in enumerate(header)
                   if h.strip().lower() in ("genes", "gene", "protein.names", "protein_names",
                                            "protein.group", "protein.ids", "protein")), 0)
        # sample columns = everything that parses as numeric on the first data row
        rows, sample_idx = [], None
        for rec in rd:
            if sample_idx is None:
                sample_idx = [i for i in range(len(rec))
                              if i != gi and _num(rec[i]) is not None]
                if not sample_idx:                    # header-only numeric detection fallback
                    sample_idx = [i for i in range(len(rec)) if i != gi]
            genes = _genes(rec[gi]) if gi < len(rec) else set()
            vals = {header[i]: _num(rec[i]) for i in sample_idx if i < len(rec)}
            rows.append({"genes": genes, "vals": vals})
        samples = [header[i] for i in sample_idx]
    return samples, rows


def _num(x):
    try:
        v = float(x)
        return v if not math.isnan(v) else None
    except (ValueError, TypeError):
        return None


def _genes(cell):
    return {g.upper() for g in re.split(r"[;,/\s]+", (cell or "").strip()) if g and g not in (".", "NA")}


def _to_log2(rows, samples):
    """DE matrices may be raw intensity or log2. If values look like raw intensity
    (95th pct > 100), log2-transform so panel means are comparable."""
    allv = [v for r in rows for v in r["vals"].values() if v is not None and v > 0]
    if not allv:
        return
    allv.sort()
    if allv[int(0.95 * (len(allv) - 1))] > 100:        # raw intensities -> log2
        for r in rows:
            for s in list(r["vals"]):
                v = r["vals"][s]
                r["vals"][s] = math.log2(v) if (v is not None and v > 0) else None


def panel_scores(rows, samples, panel_genes):
    pg = set(g.upper() for g in panel_genes)
    hits = [r for r in rows if r["genes"] & pg]
    per = {}
    for s in samples:
        vals = [r["vals"].get(s) for r in hits]
        vals = [v for v in vals if v is not None]
        per[s] = (sum(vals) / len(vals)) if vals else None
    matched = sorted({g for r in hits for g in (r["genes"] & pg)})
    return per, len(hits), matched


def zscore(per, samples):
    xs = [per[s] for s in samples if per[s] is not None]
    if len(xs) < 2:
        return {s: None for s in samples}
    m = sum(xs) / len(xs)
    sd = (sum((x - m) ** 2 for x in xs) / (len(xs) - 1)) ** 0.5 or 1e-9
    return {s: ((per[s] - m) / sd if per[s] is not None else None) for s in samples}


def load_conditions(path, samples):
    """Map matrix sample columns -> group via conditions.csv (fuzzy basename-stem)."""
    if not path or not os.path.exists(path):
        return {}
    pairs = []
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh)
        fcol = next((c for c in rd.fieldnames if "file" in c.lower() or "run" in c.lower() or "sample" in c.lower()), rd.fieldnames[0])
        gcol = next((c for c in rd.fieldnames if "group" in c.lower() or "condition" in c.lower()), rd.fieldnames[-1])
        for r in rd:
            pairs.append((_stem(r.get(fcol, "")), r.get(gcol, "").strip()))
    gmap = {}
    for s in samples:
        ss = _stem(s)
        hit = next((g for stem, g in pairs if stem and (stem == ss or stem in ss or ss in stem)), None)
        if hit:
            gmap[s] = hit
    return gmap


def _stem(x):
    b = os.path.basename(str(x).strip().rstrip("/"))
    for ext in (".mzml", ".raw", ".d", ".dia", ".parquet"):
        if b.lower().endswith(ext):
            b = b[: -len(ext)]
    return b.lower()


def confound_check(z, gmap, samples, thr):
    """Return (confounded: bool, detail) -- is the panel score separated by group with
    little overlap? That is the danger signal: DE may be contamination, not biology."""
    if not gmap:
        return False, "no conditions.csv -- group-confounding not assessed"
    groups = {}
    for s in samples:
        if s in gmap and z.get(s) is not None:
            groups.setdefault(gmap[s], []).append(z[s])
    if len(groups) < 2:
        return False, "fewer than 2 groups with panel data"
    means = {g: sum(v) / len(v) for g, v in groups.items()}
    hi = max(means, key=means.get); lo = min(means, key=means.get)
    gap = means[hi] - means[lo]
    overlap = max(min(groups[hi]), min(groups[lo])) <= min(max(groups[hi]), max(groups[lo]))  # crude
    # strong signal: group means differ by > ~1.5 SD of z (i.e. > thr) AND ranges barely overlap
    separated = (gap >= thr) and (min(groups[hi]) > max(groups[lo]) - 0.25 * gap)
    detail = "; ".join(f"{g}: mean z {means[g]:+.2f} (n={len(groups[g])})" for g in sorted(means))
    return bool(separated), f"{detail}  [gap {gap:+.2f}]"


def detected_depth(report):
    """Per-run distinct protein groups at Q<=0.01 from a DIA-NN report.parquet -- the
    REAL per-sample depth (a DPC-Quant expression matrix is complete, so counting its
    non-empty cells gives a constant, not depth)."""
    try:
        import pyarrow.parquet as pq
    except ImportError:
        return None, "pyarrow not available"
    t = pq.read_table(report)
    cols = {c.lower(): c for c in t.column_names}
    run = cols.get("run"); pg = cols.get("protein.group") or cols.get("protein.ids")
    q = cols.get("q.value") or cols.get("global.q.value")
    if not (run and pg):
        return None, "report.parquet missing Run/Protein.Group"
    runs = t.column(run).to_pylist(); pgs = t.column(pg).to_pylist()
    qs = t.column(q).to_pylist() if q else [0.0] * len(runs)
    seen = {}
    for r, p, qv in zip(runs, pgs, qs):
        if qv is None or qv <= 0.01:
            seen.setdefault(str(r), set()).add(str(p))
    return {r: len(v) for r, v in seen.items()}, None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--matrix", required=True, help="expression matrix CSV (proteins x samples; a gene/id column + numeric sample columns)")
    ap.add_argument("--conditions", help="conditions.csv (File.Name,Group) to test group-confounding")
    ap.add_argument("--report", help="DIA-NN report.parquet for TRUE per-sample detected depth")
    ap.add_argument("--out", default="SAMPLE_QUALITY.md")
    ap.add_argument("--z", type=float, default=1.5, help="|z| threshold to flag an elevated sample (default 1.5)")
    ap.add_argument("--keratin-sample", action="store_true",
                    help="sample IS keratin (nail/hair/wool/skin/feather) — keratin is the ANALYTE, "
                         "so the EPIDERMIS panel is reported for QC but never flagged as contamination")
    a = ap.parse_args()

    samples, rows = read_matrix(a.matrix)
    _to_log2(rows, samples)
    gmap = load_conditions(a.conditions, samples)

    na = sum(1 for r in rows for s in samples if r["vals"].get(s) is None)
    complete = na / max(1, len(rows) * len(samples)) < 0.005

    results, flags = {}, []
    for name, genes in PANELS.items():
        per, nhit, matched = panel_scores(rows, samples, genes)
        z = zscore(per, samples)
        elevated = sorted(s for s in samples if z.get(s) is not None and z[s] >= a.z)
        confounded, detail = confound_check(z, gmap, samples, a.z)
        expected = a.keratin_sample and name == "EPIDERMIS"
        results[name] = {"n_panel_proteins": nhit, "matched_genes": matched,
                         "z": {s: (round(z[s], 2) if z[s] is not None else None) for s in samples},
                         "elevated_samples": elevated, "group_confounded": confounded,
                         "group_detail": detail, "expected_analyte": expected}
        if expected:
            pass   # keratin IS the analyte for a keratin-matrix sample: report for QC, never flag
        elif confounded:
            flags.append(f"**{name} is CONFOUNDED WITH GROUP** ({detail}). DE between these "
                         "groups may be contamination, not biology; protein-level marker removal "
                         "will NOT fix a confounded contrast -- resolve at the sample/design level.")
        elif elevated:
            flags.append(f"{name}: elevated in {', '.join(elevated)} (|z|>={a.z}) -- possible "
                         "per-sample contamination; check before interpreting these samples.")

    depth, depth_note = (None, None)
    if a.report and os.path.exists(a.report):
        depth, depth_note = detected_depth(a.report)

    # ---- write report ----
    with open(a.out, "w") as fh:
        fh.write("# Sample-quality & contamination diagnostics\n\n")
        if not flags:
            fh.write("No contamination panel is group-confounded or per-sample elevated at the "
                     f"current threshold (|z|>={a.z}).\n\n")
        for f in flags:
            fh.write(f"- {f}\n")
        fh.write("\n")
        for name, r in results.items():
            fh.write(f"## {name}  ({r['n_panel_proteins']} panel proteins detected)\n\n")
            if r.get("expected_analyte"):
                fh.write("_This is the **analyte** for a keratin-matrix sample (nail/hair/wool/skin/"
                         "feather) — reported for QC, **not** treated as contamination._\n\n")
            if not r["matched_genes"]:
                fh.write("_None of this panel's markers were detected — not assessable._\n\n")
                continue
            fh.write(f"markers: {', '.join(r['matched_genes'])}\n\n")
            fh.write("| sample | group | z (panel abundance) |\n|---|---|---|\n")
            for s in samples:
                zz = r["z"][s]
                fh.write(f"| {s} | {gmap.get(s,'?')} | {zz if zz is not None else 'NA'} |\n")
            fh.write(f"\n_group-confounding:_ {r['group_detail']}\n\n")
        if complete:
            fh.write("## ⚠ Per-sample depth\n\nThe expression matrix is **complete "
                     "(~no missing values)** — consistent with DPC-Quant/limpa, whose detection "
                     "model fills every cell. **Do not** use non-empty cell counts as per-sample "
                     "depth (it is a constant). ")
            if depth:
                fh.write("True detected protein groups per run (Q≤0.01, from report.parquet):\n\n")
                fh.write("| run | detected protein groups |\n|---|---|\n")
                for r in sorted(depth):
                    fh.write(f"| {r} | {depth[r]} |\n")
            else:
                fh.write("Pass `--report report.parquet` for true detected-precursor depth"
                         + (f" ({depth_note})" if depth_note else "") + ".\n")
            fh.write("\n")
    with open(os.path.splitext(a.out)[0].lower().replace("sample_quality", "sample_quality") + ".json"
             if False else a.out.replace(".md", ".json"), "w") as jf:
        json.dump({"samples": samples, "groups": gmap, "matrix_complete": complete,
                   "panels": results, "detected_depth": depth, "flags": flags}, jf, indent=2)

    print(json.dumps({"out": a.out, "flags": flags,
                      "group_confounded": [n for n, r in results.items() if r["group_confounded"]],
                      "matrix_complete": complete}, indent=2))


if __name__ == "__main__":
    main()
