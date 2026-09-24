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
  * CONTAMINANT_IDENTICAL
                   proteins of the searched organism whose sequence is ALSO a common-
                   contaminant entry (human keratins, ACTB = bovine ACTB ...). fetch_fasta.py
                   removed those entries so the proteins are quantified under their own
                   accession -- they are KEPT in quantification, so a per-sample or
                   group-confounded excess is flagged as possible contamination. Built from
                   the FASTA sidecar (--fasta-meta), the list's one source.

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
      [--fasta-meta search.fasta.meta.json]   # default: ./search.fasta.meta.json if present

Reads/writes plain files; stdlib only (pyarrow optional, just for --report).
"""
import sys, os, csv, json, math, argparse, re

# The list's ONE definition is the FASTA sidecar; the wording for lost proteins lives there too.
from fetch_fasta import target_contaminants, seen_only_as_cont, lost_to_contaminants_message

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
# Not a curated panel: built per run from the FASTA sidecar (--fasta-meta).
TARGET_PANEL = "CONTAMINANT_IDENTICAL"


ID_COLUMNS = ("genes", "gene", "protein.names", "protein_names", "protein.group", "protein.ids",
              "protein")
PROTEIN_ID_COLUMNS = ("protein.group", "protein.ids", "protein")


def read_matrix(path):
    """Return (samples, rows); each row = {'genes': set, 'ids': set, 'vals': {sample: float|None}}.

    `genes` pools the tokens of EVERY id column. It used to read only the first one found,
    and run_de.R writes Protein.Group, Genes, Protein.Names (R's merge puts the key first) --
    so the gene-symbol panels were compared against UniProt accessions and never matched
    (review 2026-09-24: HBB/HBA1 +5 log2 in one sample -> HEMOLYSIS 0 proteins, no flag).
    `ids` is the protein-accession column(s) alone."""
    with open(path, newline="") as fh:
        rd = csv.reader(fh)
        header = next(rd)
        id_idx = [i for i, h in enumerate(header) if h.strip().lower() in ID_COLUMNS] or [0]
        pid_idx = [i for i in id_idx if header[i].strip().lower() in PROTEIN_ID_COLUMNS]
        # sample columns = everything that parses as numeric on the first data row
        rows, sample_idx = [], None
        for rec in rd:
            if sample_idx is None:
                sample_idx = [i for i in range(len(rec))
                              if i not in id_idx and _num(rec[i]) is not None]
                if not sample_idx:                    # header-only numeric detection fallback
                    sample_idx = [i for i in range(len(rec)) if i not in id_idx]
            genes = set().union(*(_genes(rec[i]) for i in id_idx if i < len(rec)))
            ids = set().union(*(_genes(rec[i]) for i in pid_idx if i < len(rec)))
            vals = {header[i]: _num(rec[i]) for i in sample_idx if i < len(rec)}
            rows.append({"genes": genes, "ids": ids, "vals": vals})
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
    ap.add_argument("--fasta-meta",
                    help="fetch_fasta.py's <fasta>.meta.json -- supplies the proteins that are also "
                         "common-contaminant sequences. Default: ./search.fasta.meta.json if present")
    a = ap.parse_args()
    fasta_meta = a.fasta_meta or ("search.fasta.meta.json"
                                  if os.path.exists("search.fasta.meta.json") else None)

    samples, rows = read_matrix(a.matrix)
    _to_log2(rows, samples)
    gmap = load_conditions(a.conditions, samples)

    na = sum(1 for r in rows for s in samples if r["vals"].get(s) is None)
    complete = na / max(1, len(rows) * len(samples)) < 0.005

    # Proteins of the searched organism that are also common-contaminant sequences --
    # read from the sidecar, never kept as a list here.
    tc, tc_note = None, None
    if fasta_meta:
        try:
            with open(fasta_meta) as fh:
                tc = target_contaminants(json.load(fh), a.keratin_sample)
        except (OSError, ValueError) as e:
            tc_note = f"not assessed: could not read {fasta_meta} ({e})"
    else:
        tc_note = "not assessed: no <fasta>.meta.json from fetch_fasta.py (pass --fasta-meta)"
    org = (tc or {}).get("organism") or "target-organism"
    panels = dict(PANELS)
    if tc and (tc["genes"] or tc["accessions"]):
        # Genes AND accessions: the matrix's id column may be either, and an NCBI database
        # carries no gene names at all.
        panels[TARGET_PANEL] = sorted(tc["genes"] | tc["accessions"])

    results, flags = {}, []
    for name, genes in panels.items():
        per, nhit, matched = panel_scores(rows, samples, genes)
        z = zscore(per, samples)
        elevated = sorted(s for s in samples if z.get(s) is not None and z[s] >= a.z)
        confounded, detail = confound_check(z, gmap, samples, a.z)
        expected = a.keratin_sample and name == "EPIDERMIS"
        kept = name == TARGET_PANEL
        results[name] = {"n_panel_proteins": nhit, "matched_genes": matched,
                         "z": {s: (round(z[s], 2) if z[s] is not None else None) for s in samples},
                         "elevated_samples": elevated, "group_confounded": confounded,
                         "group_detail": detail, "expected_analyte": expected,
                         "kept_in_quantification": kept}
        why = (f" These are {org} proteins that are also common-contaminant sequences: possible "
               f"contamination, KEPT in quantification and normalisation." if kept else "")
        if expected:
            pass   # keratin IS the analyte for a keratin-matrix sample: report for QC, never flag
        elif confounded:
            flags.append(f"**{name} is CONFOUNDED WITH GROUP** ({detail}). DE between these "
                         "groups may be contamination, not biology; protein-level marker removal "
                         "will NOT fix a confounded contrast -- resolve at the sample/design level."
                         + why)
        elif elevated:
            flags.append(f"{name}: elevated in {', '.join(elevated)} (|z|>={a.z}) -- possible "
                         "per-sample contamination; check before interpreting these samples." + why)
    # Real proteins that exist only as Cont_ entries (database used as-is, built with
    # --keep-target-contaminants, or built before the overlap check) -- same wording as
    # audit_results.py, from fetch_fasta.py.
    lost_msg = None
    if tc:
        seen = seen_only_as_cont(tc["kept_as_contaminant"],
                                 [(r["ids"], ";".join(sorted(r["genes"] - r["ids"])) or "?")
                                  for r in rows])
        lost_msg = lost_to_contaminants_message(tc, seen)
        if lost_msg:
            flags.append(lost_msg)

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
            if r.get("kept_in_quantification"):
                fh.write(f"_{org} proteins whose sequence is also a common-contaminant entry "
                         f"(from `{fasta_meta}`). fetch_fasta.py removed those contaminant entries "
                         f"so the proteins are quantified under their own accessions: they are "
                         f"**kept in quantification** — possible contamination (skin/hair keratins "
                         f"usually are), reported here, not excluded._\n\n")
            if not r["matched_genes"]:
                fh.write("_None of this panel's markers were detected — not assessable._\n\n")
                continue
            fh.write(f"markers: {', '.join(r['matched_genes'])}\n\n")
            fh.write("| sample | group | z (panel abundance) |\n|---|---|---|\n")
            for s in samples:
                zz = r["z"][s]
                fh.write(f"| {s} | {gmap.get(s,'?')} | {zz if zz is not None else 'NA'} |\n")
            fh.write(f"\n_group-confounding:_ {r['group_detail']}\n\n")
        if tc_note:
            fh.write(f"## {TARGET_PANEL}\n\n_Proteins that are also common-contaminant "
                     f"sequences: {tc_note}._\n\n")
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
                   "panels": results, "detected_depth": depth, "flags": flags,
                   "fasta_meta": fasta_meta, "target_contaminants_note": tc_note,
                   "lost_to_contaminants": lost_msg}, jf, indent=2)

    print(json.dumps({"out": a.out, "flags": flags,
                      "group_confounded": [n for n, r in results.items() if r["group_confounded"]],
                      "matrix_complete": complete}, indent=2))


if __name__ == "__main__":
    main()
