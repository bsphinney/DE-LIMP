#!/usr/bin/env python3
"""
make_report.py  --  Catalog every file the run produced and explain what each is.

After a run, the user gets a pile of files (search output, DE tables, the
reproducibility bundle, the AI report). This walks the output directories and
writes OUTPUT_FILES.md: one row per file with its size and a plain-language
description of what it is and how to use it. Files it doesn't recognize are still
listed (honest — never silently omit), tagged "unrecognized output".

Usage:
  python3 make_report.py --out OUTPUT_FILES.md \
      --search-out ./search_out --de-dir ./de_results \
      --repro ./reproducibility --extra ./conditions.csv ./search.fasta ./wf
"""
import sys, os, json, glob, argparse, re

# (regex on basename, category, description). First match wins.
CATALOG = [
    (r"^report\.parquet$", "Search output",
     "Normalized search result in the DIA-NN report format (protein × run, with PG.MaxLFQ and Q-values). This is the exact input to the DE step."),
    (r"^report\.tsv$", "Search output", "DIA-NN precursor report (tab-separated)."),
    (r"^report\.stats\.tsv$", "Search output", "DIA-NN per-run summary stats (IDs, proteins, precursors)."),
    (r"^report\.log\.txt$|.*\.log$", "Search output", "Search-engine run log (parameters, timing, warnings)."),
    (r".*\.speclib$|.*lib\.parquet$|.*\.predicted\.speclib$", "Search output",
     "Spectral library generated/used during the library-free search."),
    (r"^lfq\.parquet$|^results\.sage\.parquet$", "Search output", "Sage output (LFQ intensities / PSMs)."),
    (r"^combined_protein\.tsv$", "Search output", "FragPipe/IonQuant protein-level MaxLFQ table."),
    (r"^search_provenance\.json$", "Search output", "Exact search engine, version, and command used (reproducibility)."),

    (r".*\.png$", "Figures", "Publication-quality figure (volcano / PCA / heatmap / QC) embedded in the analysis report."),
    (r"^figures\.json$", "Figures", "Figure manifest: each figure's file, type, and caption."),

    (r"^DE_.*\.csv$", "Differential expression",
     "DE results for one comparison: Protein.Group, logFC, AveExpr, t, P.Value, adj.P.Val (BH), B, plus gene annotation. Sorted by adjusted p-value."),
    (r"^methods\.txt$", "Differential expression",
     "Self-describing methods paragraph (pipeline, quantification, DE engine, normalization, thresholds, citation). Paste verbatim into a paper's Methods."),
    (r"^sessionInfo\.txt$", "Differential expression",
     "Exact R + package versions (limpa/limma/arrow/...) that produced the DE — provenance."),
    (r"^de_provenance\.json$", "Differential expression",
     "Machine-readable DE record: method, design, contrasts, thresholds, per-contrast significant counts, package versions."),
    (r"^Expression_Matrix\.csv$", "Differential expression",
     "Log2 protein expression per sample (proteins × runs), the matrix DE was run on."),

    (r"^run_manifest\.json$", "Reproducibility bundle",
     "The master record — registry commit, engine + versions, all parameters, environment, input/output checksums."),
    (r"^REPRODUCE\.md$", "Reproducibility bundle", "Human-readable methods + step-by-step how to re-run."),
    (r"^reproduce\.sh$", "Reproducibility bundle",
     "Runnable script that rebuilds the env, re-fetches the pinned workflow, and re-runs search + DE."),
    (r"^MANIFEST\.txt$", "Reproducibility bundle",
     "Capture log: [OK]/[SKIPPED] for each artifact, so you can trust what the bundle contains."),
    (r"^conda-explicit\.txt$", "Reproducibility bundle", "Fully pinned conda environment lock (URL + md5 per package)."),
    (r"^pip-freeze\.txt$", "Reproducibility bundle", "Installed Python packages + versions."),
    (r"^r-sessionInfo\.txt$", "Reproducibility bundle", "R + package versions captured for the bundle."),
    (r"^versions\.txt$", "Reproducibility bundle", "Search-engine versions + resolved commands."),
    (r"^skill\.txt$", "Reproducibility bundle", "Which skill produced this analysis (name + version) and how it was installed."),
    (r"^checksums\.json$", "Reproducibility bundle", "sha256 / structural fingerprints of raw inputs, FASTA, report, and DE outputs."),
    (r".*\.rationale\.json$", "Reproducibility bundle",
     "Per-setting provenance for the estimated search parameters (which value came from the data type vs a default)."),

    (r"^conditions\.csv$", "Inputs",
     "Experimental design: File.Name → Group (+ optional Batch/Covariates) used in the DE model."),
    (r".*\.fasta$|.*\.fa$", "Inputs", "Protein sequence database used for the search (proteome + contaminants)."),
    (r"^workflow\.manifest\.json$", "Inputs", "The validated workflow that drove the run (engine, version, FASTA spec, DE method, pinned registry commit)."),
    (r".*\.cfg$|^sage_config.*\.json$|.*\.workflow$|^params\..*$", "Inputs", "Engine search parameters actually used (estimated from the data type, or a validated SOP config)."),
    (r"^commands\.log$", "Inputs", "Verbatim log of every command the run executed (audit trail)."),

    (r"^AI_Analysis_Report\.md$", "Analysis report", "The biological + QC interpretation of the results (the AI analysis)."),
    (r"^AI_Analysis_Report\.docx$", "Analysis report", "The analysis report as a Word document (same content as the .md)."),
    (r"^methods\.md$|^methods\.docx$", "Analysis report", "Publication-ready LC-MS/MS Methods section (from raw metadata) + instrument grant acknowledgment."),
    (r"^methods_params\.json$", "Analysis report", "Acquisition parameters extracted from the raw data for the Methods section."),
    (r"^ANALYSIS_PROMPT\.md$", "Analysis report", "The analysis brief the agent followed to write the report."),
    (r"^OUTPUT_FILES\.md$", "Analysis report", "This file — the catalog of all outputs."),
]

CATEGORY_ORDER = ["Analysis report", "Differential expression", "Search output",
                  "Reproducibility bundle", "Inputs", "Other"]


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def describe(basename):
    for pat, cat, desc in CATALOG:
        if re.match(pat, basename):
            return cat, desc
    return "Other", "unrecognized output (no description available)"


def collect(paths):
    seen = {}
    for p in paths:
        if not p or not os.path.exists(p):
            continue
        if os.path.isfile(p):
            seen[os.path.abspath(p)] = p
        else:
            for dp, _, fns in os.walk(p):
                for fn in fns:
                    fp = os.path.join(dp, fn)
                    seen[os.path.abspath(fp)] = fp
    return sorted(seen.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="OUTPUT_FILES.md")
    ap.add_argument("--search-out")
    ap.add_argument("--de-dir")
    ap.add_argument("--repro")
    ap.add_argument("--extra", nargs="*", default=[])
    ap.add_argument("--root", default=".", help="base dir to show paths relative to")
    a = ap.parse_args()

    files = collect([a.search_out, a.de_dir, a.repro, *a.extra])
    root = os.path.abspath(a.root)
    rows, by_cat = [], {}
    for f in files:
        if os.path.basename(f) == os.path.basename(a.out):
            continue
        cat, desc = describe(os.path.basename(f))
        try:
            size = human(os.path.getsize(f))
        except OSError:
            size = "?"
        rel = os.path.relpath(f, root)
        by_cat.setdefault(cat, []).append((rel, size, desc))
        rows.append({"file": rel, "category": cat, "size": size, "description": desc})

    lines = ["# Output files — what each one is", "",
             f"This run produced {len(rows)} file(s). Each is described below, grouped by purpose.", ""]
    n_unknown = 0
    for cat in CATEGORY_ORDER:
        items = by_cat.get(cat)
        if not items:
            continue
        lines.append(f"## {cat}")
        lines.append("")
        lines.append("| File | Size | What it is |")
        lines.append("|---|---|---|")
        for rel, size, desc in sorted(items):
            if "unrecognized" in desc:
                n_unknown += 1
            lines.append(f"| `{rel}` | {size} | {desc} |")
        lines.append("")
    if n_unknown:
        lines.append(f"> {n_unknown} file(s) were not recognized and are listed under **Other** "
                     "without a specific description.")
        lines.append("")
    lines.append("## Where to start")
    lines.append("")
    lines.append("- **`AI_Analysis_Report.md`** — read this first: the biological interpretation.")
    lines.append("- **`de_results/DE_*.csv`** — the differentially expressed proteins per comparison.")
    lines.append("- **`de_results/methods.txt`** — the Methods paragraph for your paper.")
    lines.append("- **`reproducibility/REPRODUCE.md`** — how to reproduce everything.")
    lines.append("")

    with open(a.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")

    print(json.dumps({"report": os.path.abspath(a.out), "n_files": len(rows),
                      "categories": {c: len(by_cat.get(c, [])) for c in CATEGORY_ORDER if by_cat.get(c)},
                      "unrecognized": n_unknown}, indent=2))


if __name__ == "__main__":
    main()
