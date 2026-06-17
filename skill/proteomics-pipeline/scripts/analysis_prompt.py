#!/usr/bin/env python3
"""
analysis_prompt.py  --  Build the analysis brief the agent follows to interpret
the results, ported from DE-LIMP's "Export Prompt for Claude" (R/server_ai.R).

In DE-LIMP the prompt is shipped to an external Claude/Gemini. Here the agent IS
Claude (running in Claude Code / cowork), so this writes ANALYSIS_PROMPT.md and
the agent then reads it + the data CSVs and writes AI_Analysis_Report.md inline —
no API key, no external call.

The prompt structure is faithful to DE-LIMP's. The pipeline self-description and
the educational/stats sections are chosen from the actual engine + DE method
(read from de_provenance.json) rather than hardcoded, per DE-LIMP rule #1.

Usage:
  python3 analysis_prompt.py --out ANALYSIS_PROMPT.md \
      --de-dir ./de_results --report ./search_out/report.parquet \
      --conditions ./conditions.csv [--qc ./QC_Metrics.csv] \
      --engine diann --acquisition DIA --instrument "Orbitrap Astral" \
      --workflow-manifest ./wf/workflow.manifest.json
"""
import sys, os, json, glob, argparse


def load(path):
    try:
        return json.load(open(path))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="ANALYSIS_PROMPT.md")
    ap.add_argument("--de-dir", required=True)
    ap.add_argument("--report")
    ap.add_argument("--conditions")
    ap.add_argument("--qc")
    ap.add_argument("--gsea")
    ap.add_argument("--engine", default="")
    ap.add_argument("--acquisition", default="")
    ap.add_argument("--instrument", default="")
    ap.add_argument("--workflow-manifest")
    ap.add_argument("--report-out", default="AI_Analysis_Report.md")
    a = ap.parse_args()

    prov = load(os.path.join(a.de_dir, "de_provenance.json")) or {}
    wfman = load(a.workflow_manifest) if a.workflow_manifest else {}
    method = prov.get("method") or (wfman.get("de", {}) or {}).get("method", "")
    engine = a.engine or (wfman.get("engine", {}) or {}).get("name", "")
    acq = (a.acquisition or wfman.get("acquisition", "")).upper()
    de_files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(a.de_dir, "DE_*.csv")))
    contrasts = prov.get("contrasts") or []

    # pipeline self-description (authoritative — from the DE step, not hardcoded)
    desc_line = prov.get("display_label") or "the configured quantification + limma pipeline"
    rollup = prov.get("rollup_method", "")
    de_engine = prov.get("de_engine", "")
    missing_policy = prov.get("missing_policy", "")
    citation = prov.get("citation", "")

    has_qc = bool(a.qc and os.path.exists(a.qc))
    has_gsea = bool(a.gsea and os.path.exists(a.gsea))
    is_dia = acq == "DIA" or engine == "diann"

    L = []
    w = L.append

    w("# Analysis brief — proteomics differential expression")
    w("")
    w("You are a senior proteomics and systems biology consultant. Analyze the "
      f"differential expression results from this **{desc_line}** run (engine: "
      f"`{engine or 'unknown'}`, acquisition: `{acq or 'unknown'}`). Read the data files "
      "below, then write the report described under OUTPUT.")
    w("")
    w("## Pipeline (authoritative — describe it exactly this way)")
    w(f"- Quantification: {rollup or '(see methods.txt)'}")
    w(f"- DE engine: {de_engine or '(see methods.txt)'}")
    if missing_policy:
        w(f"- Missing values: {missing_policy}")
    if citation:
        w(f"- Citation: {citation}")
    w("- The exact, self-describing methods text is in `de_results/methods.txt` — do not "
      "contradict it or invent a different pipeline.")
    w("")

    w("## Attached data files — read these")
    for f in de_files:
        w(f"- `de_results/{f}` — DE results for one comparison: Protein.Group, logFC, "
          "AveExpr, t, P.Value, adj.P.Val (BH), B, plus gene annotation.")
    w("- `de_results/methods.txt` — the methods paragraph.")
    if a.report:
        w(f"- `{a.report}` — the normalized search report (protein × run).")
    if a.conditions:
        w(f"- `{a.conditions}` — experimental design (File.Name → Group [+ Batch/Covariates]).")
    if has_qc:
        w(f"- `{a.qc}` — per-sample QC metrics.")
    if has_gsea:
        w(f"- `{a.gsea}` — Gene Set Enrichment results.")
    w("")
    w("Compute what you need directly from these files: significant proteins per "
      "contrast (apply the thresholds in `de_provenance.json`: adj.P.Val and |logFC|), "
      "up/down splits, proteins significant in ≥2 contrasts, and — using the expression "
      "matrix/report + conditions — the most stable proteins (lowest CV across replicates). "
      "Use gene names where available; cite specific proteins to support every claim. "
      "Never fabricate a value you can't find in the data.")
    w("")

    w("## OUTPUT — write `" + a.report_out + "` with these markdown sections")
    w("")
    w("### Overview")
    w("Number of comparisons analyzed, total significant proteins per comparison (up/down "
      "split). Overall assessment of the experiment's quality and scope.")
    w("")
    if has_qc:
        w("### QC Assessment")
        w("Evaluate technical quality from the QC metrics. Comment on consistency of "
          "precursor/protein identifications across replicates and groups, and flag any "
          "outlier samples or systematic biases.")
        w("")
    w("### Key Findings Per Comparison")
    w("For each comparison: highlight the top upregulated and downregulated proteins by "
      "fold-change (use gene names). Note any comparison with unusually few or many hits.")
    w("")
    w("### Cross-Comparison Biomarkers")
    w("Proteins significant in multiple comparisons are highest-confidence candidates. "
      "Discuss consistency of direction (always up, always down, or mixed).")
    w("")
    w("### High-Confidence Biomarker Insights")
    w("For the most stable proteins (lowest CV): discuss known biological functions, "
      "pathway involvement, and disease associations where you recognize the gene. Assess "
      "biomarker potential from the combination of low CV, significant p-value, and "
      "meaningful fold-change.")
    w("")
    if has_gsea:
        w("### Pathway & Gene Set Enrichment Analysis")
        w("Summarize the top enriched pathways by ontology (highest NES). Connect them to "
          "the DE findings above.")
        w("")
    w("### Biological Interpretation")
    w("Suggest what biological processes or pathways may be affected based on the protein "
      "lists. Note well-known protein families, complexes, or signaling cascades. If the "
      "data suggests a clear biological narrative, describe it.")
    w("")
    w("### How This Analysis Works")
    w("Write an educational background a biologist with no MS/bioinformatics background can "
      "follow, using analogies. Cover, in plain language:")
    w("- **LC-MS/MS**: proteins digested to peptides, separated by liquid chromatography, "
      "ionized and measured by mass; MS1 measures intact peptide mass, MS2 fragments them "
      "to read the sequence.")
    if is_dia:
        w("- **Data-Independent Acquisition (DIA)**: vs DDA (picks the loudest signals one "
          "at a time), DIA systematically scans all peptides in m/z windows — more complete, "
          "reproducible quantification; needs specialized software to deconvolve.")
        if engine == "diann":
            w("- **DIA-NN**: turns the raw DIA data into protein quantities per sample; uses "
              "neural networks to score IDs and a library-free predicted library.")
    else:
        w("- **Data-Dependent Acquisition (DDA)**: the instrument picks the most abundant "
          "precursors to fragment; explain why that under-samples low-abundance peptides.")
        if engine == "sage":
            w("- **Sage**: a fast database search engine that matches MS2 spectra to peptides "
              "from the FASTA, with FDR control and label-free quantification.")
    if method == "dpc":
        w("- **LIMPA / limma (DPC-Quant)**: limma borrows information across all proteins for "
          "better variance estimates (empirical Bayes), powerful with few replicates; limpa "
          "adds a detection-probability model so missing precursors are modelled, not imputed "
          "or dropped.")
    elif method == "maxlfq":
        w("- **MaxLFQ + limma**: MaxLFQ builds protein quantities from peptide signals; limma "
          "then tests for differences using empirical-Bayes-moderated statistics. Missing "
          "values are left in place and handled per protein.")
    w("- **Key terms** (define with examples from this dataset): log2 fold change, p-value, "
      "adjusted p-value/FDR (Benjamini-Hochberg — explain the multiple-testing problem with "
      "an intuitive coin-flip example), volcano plot, coefficient of variation, normalization.")
    w("Keep the tone approachable and encouraging; define jargon when unavoidable.")
    w("")
    w("### Methods & Reproducibility")
    w("Write a concise, publication-ready Methods paragraph from `methods.txt` and "
      "`de_provenance.json` (engine + version, FASTA, FDR, normalization, statistical test, "
      "multiple-testing correction, software versions). Note that a full reproducibility "
      "bundle (`reproducibility/`) with a pinned environment and `reproduce.sh` accompanies "
      "this analysis. Include the citation above.")
    w("")
    w("---")
    w("Use markdown headers. Be scientific but accessible. Reference specific proteins from "
      "the CSV files to support your analysis. Do not invent data, pathways, or citations "
      "you cannot ground in the attached files or established biology.")
    w("")

    with open(a.out, "w") as fh:
        fh.write("\n".join(L) + "\n")

    print(json.dumps({
        "prompt": os.path.abspath(a.out),
        "report_to_write": a.report_out,
        "engine": engine, "acquisition": acq, "de_method": method,
        "de_files": de_files, "contrasts": contrasts,
        "has_qc": has_qc, "has_gsea": has_gsea,
        "next": f"Read {a.out} and the listed data files, then write {a.report_out} following the OUTPUT sections.",
    }, indent=2))


if __name__ == "__main__":
    main()
