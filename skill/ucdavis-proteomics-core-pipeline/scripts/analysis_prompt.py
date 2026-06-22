#!/usr/bin/env python3
"""
analysis_prompt.py  --  Build the analysis brief the agent follows to interpret
the results. Faithful, complete port of DE-LIMP's "Export Prompt for Claude"
(R/server_ai.R), so the report is as thorough as the DE-LIMP app's AI export.

In DE-LIMP the prompt is shipped to an external Claude/Gemini. Here the agent IS
Claude (Claude Code / Desktop), so this writes ANALYSIS_PROMPT.md and the agent
reads it + the data files + the figures and writes a complete, figure-rich,
expert AI_Analysis_Report.md (then to_docx.py also saves it as Word).

The pipeline self-description and the educational/stats sections are chosen from
the actual engine + DE method (from de_provenance.json), never hardcoded (rule #1).

Usage:
  python3 analysis_prompt.py --out ANALYSIS_PROMPT.md \
      --de-dir output/tables --report output/search/report.parquet \
      --conditions input/conditions.csv --figures-dir output/figures \
      [--qc QC_Metrics.csv] [--gsea GSEA_Results.csv] \
      --engine diann --acquisition DIA --instrument "Orbitrap Astral" \
      --workflow-manifest input/wf/workflow.manifest.json
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
    ap.add_argument("--figures-dir")
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
    eng_ver = (wfman.get("engine", {}) or {}).get("version", "")
    acq = (a.acquisition or wfman.get("acquisition", "")).upper()
    de_files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(a.de_dir, "DE_*.csv")))
    contrasts = prov.get("contrasts") or []
    q_cut, lfc, adjp = prov.get("q_cutoff", 0.01), prov.get("logfc", 1.0), prov.get("adjp", 0.05)

    desc_line = prov.get("display_label") or "the configured quantification + limma pipeline"
    rollup = prov.get("rollup_method", "")
    de_engine = prov.get("de_engine", "")
    missing_policy = prov.get("missing_policy", "")
    citation = prov.get("citation", "")

    # figures
    figs = []
    fdir_rel = ""
    if a.figures_dir and os.path.isdir(a.figures_dir):
        fj = load(os.path.join(a.figures_dir, "figures.json"))
        if isinstance(fj, list):
            figs = fj
        fdir_rel = os.path.basename(a.figures_dir.rstrip("/")) or "figures"

    has_qc = bool(a.qc and os.path.exists(a.qc))
    has_gsea = bool(a.gsea and os.path.exists(a.gsea))
    is_dia = acq == "DIA" or engine == "diann"

    L = []
    w = L.append

    w("# Analysis brief — proteomics differential expression")
    w("")
    w("**You are a senior proteomics and systems biology consultant.** Write a "
      "comprehensive, publication-grade analysis of the differential expression "
      f"results from this **{desc_line}** run (engine: `{engine or 'unknown'}"
      f"{(' ' + eng_ver) if eng_ver else ''}`, acquisition: `{acq or 'unknown'}`). "
      "Be rigorous and specific — this report goes to the scientist who generated "
      "the data. Read every data file and figure below, then write the full report "
      "described under OUTPUT. Do not produce a short summary; produce the complete "
      "multi-section report.")
    w("")

    w("## Pipeline (authoritative — describe it exactly this way)")
    w(f"- Quantification: {rollup or '(see methods.txt)'}")
    w(f"- DE engine: {de_engine or '(see methods.txt)'}")
    if missing_policy:
        w(f"- Missing values: {missing_policy}")
    if citation:
        w(f"- Citation: {citation}")
    w(f"- Significance thresholds used: adj.P.Val < {adjp} and |logFC| ≥ {lfc} (ID FDR q ≤ {q_cut}).")
    w("- The exact, self-describing methods text is in `tables/methods.txt` — do not "
      "contradict it or invent a different pipeline.")
    w("")

    w("## Attached data files — read these")
    for f in de_files:
        w(f"- `tables/{f}` — DE results for one comparison: Protein.Group, logFC, "
          "AveExpr, t, P.Value, adj.P.Val (BH), B, gene annotation.")
    w("- `tables/Expression_Matrix.csv` — log2 protein abundance per sample.")
    w("- `tables/methods.txt`, `tables/de_provenance.json` — methods + exact versions.")
    if a.conditions:
        w(f"- `{a.conditions}` — experimental design (File.Name → Group [+ Batch/Covariates]).")
    if has_qc:
        w(f"- `{a.qc}` — per-sample QC metrics.")
    if has_gsea:
        w(f"- `{a.gsea}` — Gene Set Enrichment results.")
    w("")
    w("Compute everything you cite directly from these files: significant proteins per "
      f"contrast (apply adj.P.Val < {adjp} and |logFC| ≥ {lfc}), up/down splits, the "
      "top up/down proteins by fold change, proteins significant in ≥2 contrasts, and — "
      "from the expression matrix + conditions — the most stable proteins (lowest CV "
      "across replicates). Use gene names where available. **Never fabricate a value, "
      "protein, pathway, or citation you cannot ground in the data or established biology.**")
    w("")

    # ---- figures ----
    if figs:
        w("## Figures — embed and interpret each one")
        w(f"These publication-quality figures were generated for you in `{fdir_rel}/`. "
          "**Embed every figure** in the relevant section using markdown image syntax "
          f"(e.g. `![caption]({fdir_rel}/<file>)`) so it renders in the Word document, "
          "and **write an expert interpretation of what each shows for THIS dataset** — "
          "not a generic caption. Available figures:")
        for fig in figs:
            w(f"- `{fdir_rel}/{fig.get('file')}` ({fig.get('type')}) — {fig.get('caption')}")
        w("")
        w("Placement: volcano + p-value figures in **Key Findings Per Comparison**; PCA + "
          "per-sample counts in **QC Assessment**; the heatmap in **Cross-Comparison "
          "Biomarkers** or **Biological Interpretation**.")
        w("")

    w("## OUTPUT — write `" + a.report_out + "` with ALL of these sections (markdown)")
    w("Use markdown headers, tables for numbers, and embed the figures. Be scientific "
      "but accessible. Reference specific proteins/genes throughout.")
    w("")
    w("### Overview")
    w("Number of comparisons analyzed, total significant proteins per comparison "
      "(up/down split, as a table). Overall assessment of the experiment's quality, "
      "depth (proteins quantified), and scope.")
    w("")
    w("### QC Assessment")
    w("Evaluate technical quality. Use the PCA and per-sample protein-count figures, and "
      "the QC metrics if present. Comment on consistency of identifications across "
      "replicates and groups, whether replicates cluster, and flag any outlier samples or "
      "systematic biases (e.g. a group quantifying far fewer proteins). State clearly "
      "whether the data are fit for differential analysis.")
    w("")
    w("### Key Findings Per Comparison")
    w("For each comparison: embed its volcano plot, then highlight the top up- and "
      "down-regulated proteins by fold change (use gene names, give logFC and adj.P). "
      "Note any comparison with unusually few or many hits. Use the p-value distribution "
      "to comment on whether the statistics are well-calibrated.")
    w("")
    w("### Cross-Comparison Biomarkers")
    w("Proteins significant in multiple comparisons are the highest-confidence candidates. "
      "Discuss consistency of direction (always up, always down, or mixed). Embed the "
      "heatmap and use it to show which proteins drive the separation. (If there is only "
      "one comparison, say so and focus on the strongest, most reproducible hits.)")
    w("")
    w("### High-Confidence Biomarker Insights")
    w("For the most stable significant proteins (lowest CV): discuss their known "
      "biological functions, pathway involvement, and disease associations where you "
      "recognize the gene. Assess each one's potential as a reliable biomarker from the "
      "combination of low CV, significant adj.P, and meaningful fold change.")
    w("")
    if has_gsea:
        w("### Pathway & Gene Set Enrichment Analysis")
        w("Summarize the top enriched pathways by ontology (highest |NES|). Connect "
          "enriched pathways to the DE protein findings above.")
        w("")
    w("### Biological Interpretation")
    w("Synthesize what biological processes or pathways are affected based on the protein "
      "lists. Name well-known protein families, complexes, or signaling cascades that are "
      "represented. If the data tell a coherent biological story, describe it — and state "
      "your confidence and the caveats (sample size, missingness, single comparison, etc.).")
    w("")
    w("### How This Analysis Works")
    w("Write an educational background a biologist with no MS/bioinformatics background "
      "can follow, using analogies. Cover, in plain language:")
    w("- **LC-MS/MS**: proteins are digested into peptides, separated by liquid "
      "chromatography (sorting by 'stickiness'), then ionized and weighed; MS1 measures "
      "intact peptide mass, MS2 fragments them to read the sequence.")
    if is_dia:
        w("- **Data-Independent Acquisition (DIA)**: vs DDA (picks the loudest signals one "
          "at a time), DIA systematically scans all peptides in m/z windows every cycle — "
          "more complete, reproducible quantification; it needs specialized software to "
          "deconvolve the complex spectra.")
        if engine == "diann":
            w("- **DIA-NN**: turns the raw DIA data into protein quantities per sample, using "
              "neural networks to score identifications and a library-free predicted spectral "
              "library (it predicts what peptides should look like rather than needing a "
              "pre-built library).")
    else:
        w("- **Data-Dependent Acquisition (DDA)**: the instrument repeatedly picks the most "
          "abundant precursors to fragment, which under-samples low-abundance peptides and "
          "makes quantification less reproducible run-to-run.")
        if engine == "sage":
            w("- **Sage**: a very fast database search engine that matches each MS2 spectrum to "
              "peptides from the FASTA, with target-decoy FDR control and label-free quant.")
    if method == "dpc":
        w("- **LIMPA / limma (DPC-Quant)**: limma borrows information across all proteins for "
          "better variance estimates (empirical Bayes moderation), which is powerful with the "
          "few replicates typical in proteomics; limpa adds a detection-probability model so "
          "missing precursors are modelled rather than imputed or dropped.")
    elif method == "maxlfq":
        w("- **MaxLFQ + limma**: MaxLFQ reconstructs protein quantities from peptide signals; "
          "limma then tests for differences with empirical-Bayes-moderated statistics. Missing "
          "values are left in place and handled per protein.")
    w("- **Key statistical concepts** — define each in plain language with an example from "
      "THIS dataset:")
    w("  - **log2 fold change (logFC)**: how much a protein goes up/down between groups "
      "(logFC 1 = doubled, −1 = halved).")
    w("  - **p-value**: the probability of seeing this difference by chance alone.")
    w("  - **adjusted p-value (FDR, Benjamini–Hochberg)**: p-values corrected for testing "
      "thousands of proteins at once — explain the multiple-testing problem with an intuitive "
      "coin-flip example.")
    w("  - **volcano plot**: why it's volcano-shaped and how to read it (x = effect size, "
      "y = significance).")
    w("  - **coefficient of variation (CV)**: a measure of measurement reproducibility — "
      "lower is more reliable.")
    w("  - **normalization**: why raw intensities need correcting (sample loading differences) "
      "and, briefly, how this pipeline normalizes.")
    w("Keep the tone approachable and encouraging; define jargon when unavoidable.")
    w("")
    w("### Methods (publication-ready) & Reproducibility")
    w("Write a concise Methods paragraph in third-person past tense suitable for a journal, "
      "drawn from `methods.txt` and `de_provenance.json`: the search engine + version, FASTA "
      "database, enzyme/missed cleavages and modifications if available, mass-accuracy "
      "settings, FDR threshold, quantification method, normalization, the statistical test, "
      "and multiple-testing correction. Cite the pipeline appropriately (the citation above).")
    if is_dia and engine == "diann":
        w("Cite DIA-NN (Demichev et al., Nature Methods, 2020).")
    w("Note that a full reproducibility bundle accompanies this analysis "
      "(`output/reproducibility/`): which skill produced it and how it was installed, a "
      "pinned conda environment lock, R sessionInfo, exact tool versions, input/output "
      "checksums, and a runnable `reproduce.sh`. State that the analysis was produced by the "
      "**UC Davis Proteomics Core pipeline** Claude skill and can be reproduced from that "
      "bundle (see `reproducibility/REPRODUCE.md`).")
    if a.instrument:
        w("")
        w("### Instrument & Acquisition")
        w(f"Instrument detected: **{a.instrument}**. Write a short Sample Preparation & Data "
          "Acquisition note (instrument model, acquisition type) for the Methods, in "
          "third-person past tense.")
    w("")
    w("---")
    w("Reference specific proteins from the CSVs to support every claim. Embed every figure. "
      "Do not invent data, pathways, or citations you cannot ground in the attached files or "
      "established biology. After writing the report, it will also be saved as a Word document.")
    w("")

    with open(a.out, "w") as fh:
        fh.write("\n".join(L) + "\n")

    print(json.dumps({
        "prompt": os.path.abspath(a.out),
        "report_to_write": a.report_out,
        "engine": engine, "engine_version": eng_ver, "acquisition": acq, "de_method": method,
        "de_files": de_files, "contrasts": contrasts,
        "figures": [f.get("file") for f in figs], "n_figures": len(figs),
        "has_qc": has_qc, "has_gsea": has_gsea,
        "next": f"Read {a.out} + the data files + figures, then write {a.report_out} "
                f"(ALL sections, embed all {len(figs)} figures, expert interpretation), "
                "then convert it to .docx with to_docx.py.",
    }, indent=2))


if __name__ == "__main__":
    main()
