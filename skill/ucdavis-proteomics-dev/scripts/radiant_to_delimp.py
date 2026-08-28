#!/usr/bin/env python3
"""
radiant_to_delimp.py -- turn Radiant/Fulcrum output into a report.parquet the DE-LIMP
app (Hugging Face / local Shiny) will load.

WHY THIS IS NOT `adapt_radiant`
-------------------------------
The two consumers want different SHAPES, and conflating them breaks one of them:

  adapt_radiant   -> one row per protein x run, with PG.MaxLFQ. That is the skill's
                     own DE contract (references/de-analysis.md).
  DE-LIMP         -> `fileInput("report_file", "DIA-NN Report (.parquet)")` feeds
                     `limpa::readDIANN()`, which is PRECURSOR-level and needs:
                       run.column        = Run
                       feature.column    = Precursor.Id
                       intensity.column  = Precursor.Normalised
                       annotation.columns= Protein.Group, Protein.Names, Genes, Proteotypic
                       q.columns         = Q.Value, Lib.Q.Value, Lib.PG.Q.Value

Fulcrum's `combined` output is already precursor-level with DIA-NN-style names, so most
columns pass straight through. Verified against a real 3-run HeLa result (2026-08-03),
these four had to be derived:

  Precursor.Id    <- Modified.Sequence + Precursor.Charge  (how DIA-NN builds it)
  Lib.Q.Value     <- Global.Precursor.Q.Value
  Lib.PG.Q.Value  <- Global.PG.Q.Value
  Protein.Names / Genes  <- joined from the spectral library, which carries them;
                            Fulcrum drops them. Without the join DE-LIMP loses all
                            gene-level annotation (GSEA, gene labels on the volcano),
                            so we fill them rather than emitting blanks.

DE-LIMP also compares Precursor.Quantity vs Precursor.Normalised to detect whether
normalisation was applied (server_data.R) — both are carried through.

Usage:
  python3 radiant_to_delimp.py --results <out>/radiant_results --out delimp_report.parquet
  python3 radiant_to_delimp.py --results ... --library radiant_library.tsv --out ...
"""
import sys, os, json, argparse, csv

# Fulcrum column -> DIA-NN/limpa column, where a straight rename suffices.
PASSTHROUGH = {
    "Protein.Group":        "Protein.Group",
    "Proteotypic":          "Proteotypic",
    "Q.Value":              "Q.Value",
    "Precursor.Quantity":   "Precursor.Quantity",
    "Precursor.Normalised": "Precursor.Normalised",
    "Precursor.Charge":     "Precursor.Charge",
    "Modified.Sequence":    "Modified.Sequence",
    "RT":                   "RT",
    "PG.Quantity":          "PG.Quantity",
    "PG.Normalised":        "PG.Normalised",
}


def run_name(raw):
    """Fulcrum reports Run as file:///.../Sample.mzML.radiantDIA — reduce to `Sample`."""
    s = str(raw)
    if "://" in s:
        s = s.split("://", 1)[1]
    s = os.path.basename(s.rstrip("/"))
    for _ in range(3):
        stem, ext = os.path.splitext(s)
        if ext.lower() in (".radiantdia", ".mzml", ".raw", ".d", ".gz", ".parquet"):
            s = stem
        else:
            break
    return s


def ms1_from_per_file(results_root):
    """(run, peptide, charge) -> MS1 intensity, read from Radiant's own per-file output.

    Fulcrum's `combined` backend DROPS every MS1 column, so a report built from it alone
    has no MS1 signal and DE-LIMP's MS1_Signal QC metric (sum of `Ms1.Apex.Area`) comes
    out blank. Radiant itself does compute MS1 -- its `.radiantDIA` per-file result is a
    PARQUET table carrying 27 MS1 columns. `Ms1IntensityFound100` is the monoisotopic
    intensity and the closest analogue of DIA-NN's Ms1.Apex.Area, so we join it back.

    Verified: the (PeptideStringWithMods, Charge) key matches Fulcrum's
    (Modified.Sequence, Precursor.Charge) for 100% of rows.
    """
    import pyarrow.parquet as pq
    per_file = os.path.join(results_root, "radiant-results")
    if not os.path.isdir(per_file):
        return {}
    out = {}
    for fn in sorted(os.listdir(per_file)):
        if not fn.endswith(".radiantDIA"):
            continue
        run = run_name(fn)
        try:
            pf = pq.ParquetFile(os.path.join(per_file, fn))
            names = pf.schema_arrow.names
            cols = [c for c in ("PeptideStringWithMods", "Charge", "Ms1IntensityFound100")
                    if c in names]
            if len(cols) < 3:
                continue
            t = pf.read(columns=cols)
            for p, c, m in zip(t.column("PeptideStringWithMods").to_pylist(),
                               t.column("Charge").to_pylist(),
                               t.column("Ms1IntensityFound100").to_pylist()):
                if m:                       # keep the best MS1 seen for the precursor
                    k = (run, str(p), c)
                    if m > out.get(k, 0):
                        out[k] = float(m)
        except Exception as e:
            sys.stderr.write(f"[radiant_to_delimp] could not read MS1 from {fn}: {e}\n")
    return out


def library_annotation(lib_tsv):
    """ModifiedPeptide -> (ProteinName, Genes) from the Radiant TSV library.

    Streamed and de-duplicated on the peptide, so a 50M-row / 10 GB library costs
    memory proportional to the number of PEPTIDES, not rows.
    """
    ann = {}
    with open(lib_tsv, newline="") as fh:
        rd = csv.reader(fh, delimiter="\t")
        hdr = next(rd)
        idx = {c: i for i, c in enumerate(hdr)}
        pep = idx.get("ModifiedPeptide", idx.get("ModifiedPeptideSequence"))
        pn, gn = idx.get("ProteinName"), idx.get("Genes")
        if pep is None or (pn is None and gn is None):
            return ann
        for row in rd:
            k = row[pep]
            if k not in ann:
                ann[k] = (row[pn] if pn is not None and pn < len(row) else "",
                          row[gn] if gn is not None and gn < len(row) else "")
    return ann


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", required=True,
                    help="the run's radiant_results dir (containing fulcrum-results/)")
    ap.add_argument("--library", help="Radiant TSV library, to recover Protein.Names/Genes")
    ap.add_argument("--out", required=True, help="report.parquet to upload to DE-LIMP")
    ap.add_argument("--no-ms1", action="store_true",
                    help="skip the MS1 join from the per-file .radiantDIA results")
    a = ap.parse_args()

    try:
        import pyarrow as pa, pyarrow.parquet as pq, pyarrow.dataset as ds
    except ImportError:
        sys.exit("pyarrow required. pip install pyarrow")

    root = a.results
    if os.path.basename(root) != "fulcrum-results":
        cand = os.path.join(root, "fulcrum-results")
        if os.path.isdir(cand):
            root = cand
    if not os.path.isdir(root):
        sys.exit(f"no fulcrum-results under {a.results}")

    t = ds.dataset(root, format="parquet").to_table()
    cols = set(t.column_names)
    if "Run" not in cols:
        sys.exit(f"Fulcrum output has no Run column; saw {sorted(cols)}")

    out = {}
    out["Run"] = [run_name(r) for r in t.column("Run").to_pylist()]
    for src, dst in PASSTHROUGH.items():
        if src in cols:
            out[dst] = t.column(src).to_pylist()

    # Precursor.Id: DIA-NN's feature key is modified sequence + charge.
    mods = out.get("Modified.Sequence") or [""] * t.num_rows
    chgs = out.get("Precursor.Charge") or [0] * t.num_rows
    out["Precursor.Id"] = [f"{m}{c}" for m, c in zip(mods, chgs)]

    # limpa filters on all three q-columns; map Fulcrum's global estimates onto the
    # library-level names DIA-NN uses, rather than leaving them absent.
    q = out.get("Q.Value") or [0.0] * t.num_rows
    out["Lib.Q.Value"] = (t.column("Global.Precursor.Q.Value").to_pylist()
                          if "Global.Precursor.Q.Value" in cols else list(q))
    out["Lib.PG.Q.Value"] = (t.column("Global.PG.Q.Value").to_pylist()
                             if "Global.PG.Q.Value" in cols else list(q))
    if "Global.PG.Q.Value" in cols:
        out["PG.Q.Value"] = t.column("Global.PG.Q.Value").to_pylist()

    # Protein.Names / Genes — Fulcrum drops them; recover from the library if given.
    ann = library_annotation(a.library) if a.library and os.path.exists(a.library) else {}
    if ann:
        pairs = [ann.get(m, ("", "")) for m in mods]
        out["Protein.Names"] = [p[0] for p in pairs]
        out["Genes"] = [p[1] for p in pairs]
        n_ann = sum(1 for p in pairs if p[0] or p[1])
    else:
        # Never invent annotation. Empty strings are honest; a fake gene name would
        # silently corrupt GSEA and every gene label downstream.
        out["Protein.Names"] = [""] * t.num_rows
        out["Genes"] = [""] * t.num_rows
        n_ann = 0

    # MS1: recovered from Radiant's per-file output, since Fulcrum drops it.
    ms1 = {} if a.no_ms1 else ms1_from_per_file(a.results)
    if ms1:
        out["Ms1.Apex.Area"] = [
            ms1.get((r, m, c), 0.0)
            for r, m, c in zip(out["Run"], mods, chgs)]
        n_ms1 = sum(1 for v in out["Ms1.Apex.Area"] if v > 0)
    else:
        n_ms1 = 0

    pq.write_table(pa.table(out), a.out)
    runs = sorted(set(out["Run"]))
    print(json.dumps({
        "report": os.path.abspath(a.out),
        "rows": t.num_rows,
        "runs": len(runs),
        "run_names": runs,
        "protein_groups": len(set(out.get("Protein.Group", []))),
        "precursors": len(set(out["Precursor.Id"])),
        "annotated_rows": n_ann,
        "ms1_rows": n_ms1,
        "columns": list(out.keys()),
        "upload": "DE-LIMP -> 'DIA-NN Report (.parquet)'",
        "ms1_note": ("Ms1.Apex.Area recovered from Radiant's per-file output — Fulcrum's "
                     "combined backend drops all MS1 columns, which is why DE-LIMP's "
                     "MS1_Signal is blank without this join.") if n_ms1 else
                    ("NO MS1 — DE-LIMP's MS1_Signal QC metric will be blank. Radiant does "
                     "compute MS1; check that the .radiantDIA per-file results are still "
                     "beside the fulcrum-results directory."),
        "warning": None if ann else
        "Protein.Names/Genes are EMPTY — pass --library to recover them, or gene-level "
        "features in DE-LIMP (GSEA, gene labels) will have nothing to work with.",
    }, indent=2))


if __name__ == "__main__":
    main()
