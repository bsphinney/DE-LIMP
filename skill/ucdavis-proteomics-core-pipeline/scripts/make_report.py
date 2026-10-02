#!/usr/bin/env python3
"""
make_report.py  --  Catalog every file the run produced and explain what each is.

After a run, the user gets a pile of files (search output, DE tables, the
reproducibility bundle, the AI report). This walks the output directories and
writes OUTPUT_FILES.md: one row per file with its size and a plain-language
description of what it is and how to use it. Files it doesn't recognize are still
listed (honest — never silently omit), tagged "unrecognized output". The search's
per-run working files (the 5-step chain's .quant and xic/ folders) are one line per
kind (INTERNAL_KINDS), and a podcast/ folder beside --out is listed when present.

Usage:
  python3 make_report.py --out OUTPUT_FILES.md \
      --search-out ./search_out --de-dir ./de_results \
      --repro ./reproducibility --extra ./conditions.csv ./search.fasta ./wf
"""
import sys, os, json, glob, argparse, re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from report_files import SHARE_NAME  # noqa: E402  defined once; no make_podcast / notify_slack

# (regex on basename, category, description). First match wins.
CATALOG = [
    (r"^report\.parquet$", "Search output",
     "The search result the DE step reads, in the DE-input column layout (Run, Protein.Group, "
     "PG.MaxLFQ, q-values). From DIA-NN it is DIA-NN's own report; from another engine it is "
     "that engine's quantities under those column names -- de_provenance.json / methods.txt say "
     "what they are."),
    (r"^report\.tsv$", "Search output", "DIA-NN precursor report (tab-separated)."),
    (r".*\.stats\.tsv$", "Search output", "DIA-NN per-run summary stats (IDs, proteins, precursors)."),
    (r".*\.log\.txt$|.*\.log$", "Search output", "Search-engine run log (parameters, timing, warnings)."),
    (r"^report\..*_matrix\.tsv$", "Search output",
     "DIA-NN quantity matrix, runs as columns (pg = protein groups, pr = precursors, gg and "
     "unique_genes = genes)."),
    (r"^report\.protein_description\.tsv$", "Search output", "DIA-NN protein descriptions."),
    (r"^report\.manifest\.txt$", "Search output", "DIA-NN report manifest."),
    (r"^empirical\.parquet$", "Search output",
     "The empirical spectral library assembled from these runs (DIA-NN 5-step chain, step 3); "
     "the final pass searched against it."),
    (r".*\.predicted\.speclib$", "Search output",
     "DIA-NN's in-silico predicted spectral library (from the FASTA). Left out of the session "
     "zip: the FASTA, the pinned engine and the parameters rebuild it exactly."),
    (r".*\.speclib$|.*lib\.parquet$", "Search output",
     "Spectral library generated/used during the search."),
    (r"^lfq\.parquet$|^results\.sage\.parquet$", "Search output", "Sage output (LFQ intensities / PSMs)."),
    (r"^combined_protein\.tsv$", "Search output", "FragPipe/IonQuant protein-level MaxLFQ table."),
    (r"^search_provenance\.json$", "Search output", "Exact search engine, version, and command used (reproducibility)."),
    (r"^sage_adapt\.json$", "Search output",
     "Sage only: which lfq.parquet rows became report.parquet -- targets at q_value <= 0.05, "
     "decoys and failing rows dropped -- with the counts per file."),
    (r"^sage_lfq_check\.json$", "Search output",
     "Sage only: whether Sage's LFQ window (quant.lfq_settings.ppm_tolerance) fit each run's "
     "median precursor mass error, and how many target MS1 peaks Sage kept at 5% FDR "
     "(sage_lfq_check.py). A warning here means the identifications stand but the quantities "
     "are unreliable."),

    (r"^qc_pvalue_panel\.png$", "Figures",
     "QC panel: the raw p-value distribution of every contrast as small multiples "
     "(make_figures.R) -- a calibration check for the appendix, where a flat background with a "
     "peak near 0 means real signal on a well-behaved model, and a skew toward 1, a mid-range "
     "hump or a U shape flags a model/QC problem in that contrast."),
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
     "Log2 protein expression per sample (proteins × runs), the matrix DE was run on. With DPC "
     "every cell has a value, measured or inferred -- Detection_Matrix.csv says which."),
    (r"^Detection_Matrix\.csv$", "Differential expression",
     "Which values of Expression_Matrix.csv were measured and which inferred, per protein and "
     "sample (same rows and columns). What a 0 means is recorded in de_provenance.json "
     "`detection_matrix` (DPC: precursors observed, 0 = inferred by the model)."),
    (r"^QC_detected_vs_inferred\.csv$", "Differential expression",
     "Per sample: proteins detected (at least one precursor observed) vs inferred by the DPC "
     "model, as counts and percentages -- the real depth of each run."),
    (r"^QC_contaminant_share\.csv$", "Differential expression",
     "Per run: the contaminants' share of the signal -- Contaminant.Pct, the "
     "contaminant and total intensity, and the intensity column used."),
    (r"^contaminants_removed\.csv$", "Differential expression",
     "The contaminant protein groups removed before quantification. A real protein that sat "
     "in the database only as a contaminant entry is listed here, not in the DE tables."),
    (r"^Sets_baitnorm_.*\.csv$", "Protein-set tests",
     "Pulldown only: protein-set tests between conditions for one bait, relative to its complex "
     "(each IP offset by its bait's interactome median; the bait protein as a cross-check). "
     "Columns: sets_provenance.json `columns`."),
    (r"^Sets_.*\.csv$", "Protein-set tests",
     "Protein-set tests for one comparison: camera (competitive) and fry (self-contained) side "
     "by side, the same tests with run depth in the model, presence-call fractions and a "
     "measured-only version. Columns: sets_provenance.json `columns`."),
    (r"^sets_provenance\.json$", "Protein-set tests",
     "The protein-set tests' record: set sources and versions, mapping rate, settings (defaults "
     "tagged), per-contrast counts, the pulldown record, how to read them, methods."),
    (r"^sets_methods\.txt$", "Protein-set tests", "The protein-set tests' Methods paragraph."),
    (r"^sets_members\.csv$", "Protein-set tests",
     "Every tested set's proteins (gene symbols), the same for every comparison."),
    (r"^set_test_inputs\.rds$", "Differential expression",
     "The DE model's inputs for the protein-set tests (run_sets.R reads it): expression, "
     "weights, design, contrasts, block, statistics."),
    (r"^DE-LIMP_session\.rds$", "Differential expression",
     "The analysis as a DE-LIMP session: load it in the DE-LIMP app "
     "(https://delimp.stan-proteomics.org/) to explore the results interactively."),
    (r"^sample_labels\.csv$", "Figures",
     "The short sample names the figures use, and the run file each one stands for "
     "(File.Name, Label, Source)."),
    (r"^reproducibility_log\.R$", "Differential expression",
     "THE ANALYSIS AS PLAIN R — the whole DE written out top to bottom with every value "
     "literal (report path, FDR cutoff, sample→group map, design, contrasts). Read it to "
     "see exactly what was done; run it with `Rscript reproducibility_log.R` to redo it "
     "using only R and limpa/limma — no conda, no skill install."),

    (r"^run_manifest\.json$", "Reproducibility bundle",
     "The master record — skill version and search-defaults version, engine + versions, all "
     "parameters, environment, input/output checksums."),
    (r"^REPRODUCE\.md$", "Reproducibility bundle", "Human-readable methods + step-by-step how to re-run."),
    (r"^reproduce\.sh$", "Reproducibility bundle",
     "Runnable script that rebuilds the env, re-derives the search defaults from the data "
     "type (they ship with the skill -- nothing is fetched), re-resolves the pinned engine, "
     "rebuilds the FASTA and re-runs search + DE."),
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
    (r"^workflow\.manifest\.json$", "Inputs", "The search defaults that drove the run, derived from the data type (engine, version, FASTA spec, DE method, the skill's defaults version)."),
    (r".*\.cfg$|^sage_config.*\.json$|.*\.workflow$|^params\..*$", "Inputs", "Engine search parameters actually used (estimated from the data type, or a validated SOP config)."),
    (r"^commands\.log$", "Inputs", "Verbatim log of every command the run executed (audit trail)."),

    (r"^AI_Analysis_Report\.md$", "Analysis report", "The biological + QC interpretation of the results (the AI analysis)."),
    (r"^Analysis_Report\.html$", "Analysis report",
     "THE REPORT: QC panels, figures and the interpretation in one self-contained page "
     "(double-click it)."),
    (r"^AI_Analysis_Report\.docx$", "Analysis report",
     "An older Word copy of the analysis report, from before the skill stopped making one "
     "(Word mangled the figures); Analysis_Report.html is the report."),
    (r"^Analysis_Report\.pdf$", "Analysis report", "The full report with its figures as a PDF, for NotebookLM or printing (the HTML through its print stylesheet)."),
    (r"^Analysis_Report\.md$", "Analysis report", "The full report as plain text, for NotebookLM or other AI notebooks: same sections as the HTML, each figure's caption and the numbers it shows written out, top proteins per contrast."),
    (r"^methods\.md$|^methods\.docx$", "Analysis report", "Publication-ready LC-MS/MS Methods section (from raw metadata) + instrument grant acknowledgment."),
    (r"^methods_params\.json$", "Analysis report", "Acquisition parameters extracted from the raw data for the Methods section."),
    (r"^ANALYSIS_PROMPT\.md$", "Analysis report", "The analysis brief the agent followed to write the report."),
    (r"^OUTPUT_FILES\.md$", "Analysis report", "This file — the catalog of all outputs."),
    (r"^README\.html$", "Analysis report",
     "START HERE: the session README as a web page (double-click it) — summary, links to the "
     "report, Methods and tables, and where everything lives on HIVE."),
    (r"^README\.md$", "Analysis report", "The same README as plain text (Markdown)."),
    (r"^AGENTS\.md$", "Analysis report",
     "For an AI agent given this folder: the study, which file is authoritative for what, the "
     "table columns, the traps, and what it must not do."),
    ("^" + re.escape(SHARE_NAME) + "$", "Analysis report",
     "The report with the podcast's audio and transcript built in: one file to send on by itself "
     "(make_podcast.py share). The same report as Analysis_Report.html, plus the audio."),
    (r"^podcast\.(m4a|wav)$", "Analysis report",
     "OPTIONAL: an AI-generated audio discussion of the report (synthetic voices; the report is "
     "the record). make_podcast.py."),
    (r"^conversation\.md$", "Analysis report",
     "The analysis conversation, readable (save_transcript.py): the user's messages, the "
     "assistant's replies and each command. CORE-INTERNAL: never delivered, not in the zip."),
    (r"^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}\.jsonl$", "Analysis report",
     "A Claude Code transcript of the analysis, redacted (save_transcript.py). CORE-INTERNAL."),
    (r"^decisions\.md$", "Analysis report",
     "The decisions log (log_decision.py): what was decided at each decision point, and why."),
    (r"^podcast_script\.md$", "Analysis report",
     "The podcast's script, with its 'Claims beyond the report' ledger (what it says that the "
     "report does not)."),
    (r"^transcript\.html$", "Analysis report", "The podcast transcript as a web page."),
    (r"^podcast\.json$", "Analysis report",
     "How the podcast was made: TTS service and model, voices, consent, the script's checksum."),
    (r"^check\.txt$", "Analysis report",
     "The podcast script's token check against the report (numbers, symbols, names)."),
    (r"^verify\.txt$", "Analysis report",
     "The podcast audio heard back by an ASR: word match, dropped spans, numbers misread."),
    (r"^verify_transcript\.txt$", "Analysis report",
     "What the ASR heard in the podcast audio (for verify.txt)."),
    (r"^raw_files\.txt$", "Inputs",
     "Where the raw data files are (full paths); raw data is never copied into the session."),
    (r"^submission\.json$", "Inputs",
     "The CoreOmics submission this analysis answers (submission_report.py): the sample sheet, "
     "who prepared the samples, the description as written. Contacts and billing are left out."),
    (r"^samples\.tsv$", "Inputs",
     "The submission's sample sheet as a table, for people to read (submission.json is the record)."),
    (r"^session\.json$", "Inputs",
     "Which CoreOmics submission this session answers (its id and link)."),
    (r"^acquisition\.json$", "Inputs",
     "Step 2's detection (detect_acquisition.py): each raw file's acquisition, instrument and "
     "acquisition time, as read from the file."),

    (r"^AUDIT\.(md|json)$", "Analysis report",
     "The pitfall audit (audit_results.py): PASS / WARN / FAIL per check -- replication, "
     "confounding, contamination, DE signal. Read every WARN and FAIL before the results."),
    (r"^SAMPLE_QUALITY\.(md|json)$", "Analysis report",
     "Sample quality and contamination flags per sample (sample_quality.py): hemolysis, muscle, "
     "skin, contaminant-identical proteins."),
    (r"^DIFFERENCES\.md$", "Analysis report",
     "What this re-analysis changed from the original: settings, database, significant counts."),

    (r"^HOW_TO_SUBMIT\.(md|html)$", "Data deposit (PRIDE / MassIVE)",
     "START HERE to deposit the data: the steps for PRIDE (or MassIVE), what to upload and "
     "what to fill in first."),
    (r"^sdrf\.tsv$", "Data deposit (PRIDE / MassIVE)",
     "SDRF-Proteomics sample metadata for the deposit. Fill every TO-FILL cell first."),
    (r"^protocols\.txt$", "Data deposit (PRIDE / MassIVE)",
     "The sample- and data-processing protocols to paste into the submission form."),
    (r"^files_to_upload\.tsv$", "Data deposit (PRIDE / MassIVE)",
     "Every file to upload: its PRIDE file type, where it is, and what to do to it first."),
    (r"^prepare_upload\.sbatch$", "Data deposit (PRIDE / MassIVE)",
     "A SLURM job (not run by the skill) that archives each .d run and checksums the upload."),
    (r"^methods_complete_draft\.md$", "Analysis report",
     "A complete Methods draft, written beside a hand-edited methods.md that lacks a section."),

    # the search's own working files (run_search.py / diann_parallel.py's 5-step chain)
    (r"^step3_assembly\.parquet$", "Search internals",
     "DIA-NN's report of the library-assembly pass (step 3); the results are report.parquet."),
    (r".*\.sbatch$|^submit\.sh$", "Search internals",
     "A SLURM job script of the search, as submitted (submit.sh submitted the chain in order)."),
    (r"^window\.(json|txt)$", "Search internals",
     "The DIA-NN scan window measured on these runs (step 1b): window.json records every probe, "
     "window.txt is the value passed as --window -- or `auto` when the measurement failed and "
     "DIA-NN chose the window per run (probe_fallback.json)."),
    (r".*\.stale-\d{8}T\d{6}(\.\d+)?$", "Search internals",
     "An earlier search's file, set aside (renamed, never deleted) when this search was generated "
     "or re-run into the same folder. It does not describe this search."),
    (r"^probe_fallback\.json$", "Search internals",
     "The pre-search measurement failed: what the search ran with instead (not measured)."),
    (r"^(file_list|parallel_input_files)\.txt$", "Search internals",
     "The raw files the search read, one path per line."),
    (r"^jobs\.txt$", "Search internals", "The search's SLURM job ids (watch_run.sh --all reads it)."),
    (r"^fran_deposit\.json$", "Search internals",
     "This search's FRAN hand-off receipt (fran_deposit.py)."),
]

# Kinds that come by the hundred -- the 5-step chain's per-run files. Each is ONE line of
# OUTPUT_FILES.md (how many files, how much space), never hundreds of rows: (regex on the path
# relative to --root, the line's name for them from that path, what they are).
INTERNAL_KINDS = [
    (r"(^|/)xic/", lambda rel: rel[:rel.index("xic/") + 4],
     "DIA-NN extracted ion chromatograms and mobilograms (step 4, --xic): per array task a "
     "t<N>.parquet report and a t<N>_xic/ folder with each run's .xic / mobilogram parquet."),
    (r"\.quant$", lambda rel: (os.path.dirname(rel) or ".") + "/*.quant",
     "DIA-NN per-run .quant intermediates (quant_step2 = first pass, quant_step4 = final pass, "
     "quant_step2_orig = step 3's copy of quant_step2). Kept out of the session zip."),
]

CATEGORY_ORDER = ["Analysis report", "Differential expression", "Protein-set tests", "Figures",
                  "Search output", "Search internals",
                  "Reproducibility bundle", "Inputs", "Data deposit (PRIDE / MassIVE)", "Other"]


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


FINALIZE_BEGIN = "<!-- written at finalize: BEGIN -->"
FINALIZE_END = "<!-- written at finalize: END -->"


def add_finalize_files(output_files_md, files, root):
    """Add the files session.py finalize writes after this catalog was made (README.html,
    AGENTS.md, ...) to OUTPUT_FILES.md, described by the same CATALOG. Re-running replaces the
    block. Returns the number of files listed."""
    with open(output_files_md, encoding="utf-8") as fh:
        text = fh.read()
    if FINALIZE_BEGIN in text:
        head, _, rest = text.partition(FINALIZE_BEGIN)
        text = head.rstrip("\n") + "\n" + rest.partition(FINALIZE_END)[2].lstrip("\n")
    rows = []
    for f in files:
        if os.path.isfile(f):
            _, desc = describe(os.path.basename(f))
            rows.append(f"| `{os.path.relpath(f, root)}` | {human(os.path.getsize(f))} | {desc} |")
    block = [FINALIZE_BEGIN, "## Written when the session was finalized", "",
             "| File | Size | What it is |", "|---|---|---|", *rows, "", FINALIZE_END]
    with open(output_files_md, "w", encoding="utf-8") as fh:
        fh.write(text.rstrip("\n") + "\n\n" + "\n".join(block) + "\n")
    return len(rows)


def collect(paths):
    """Every file under `paths`, less scratch (scratch_files.py: .cache folders, *.part,
    podcast.wav beside podcast.m4a)."""
    import scratch_files
    seen = {}
    for p in paths:
        if not p or not os.path.exists(p):
            continue
        if os.path.isfile(p):
            if not scratch_files.is_scratch_file(os.path.basename(p), os.listdir(
                    os.path.dirname(os.path.abspath(p)))):
                seen[os.path.abspath(p)] = p
        else:
            for dp, dns, fns in os.walk(p):
                fns, _ = scratch_files.prune(dp, dns, fns)
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

    # the podcast (make_podcast.py) sits in podcast/ beside OUTPUT_FILES.md: listed when present
    podcast = os.path.join(os.path.dirname(os.path.abspath(a.out)), "podcast")
    files = collect([a.search_out, a.de_dir, a.repro, *a.extra,
                     podcast if os.path.isdir(podcast) else None])
    root = os.path.abspath(a.root)
    rows, by_cat, kinds = [], {}, {}
    for f in files:
        if os.path.basename(f) == os.path.basename(a.out):
            continue
        rel = os.path.relpath(f, root)
        try:
            nbytes = os.path.getsize(f)
        except OSError:
            nbytes = None
        kind = next(((name(rel.replace(os.sep, "/")), desc) for pat, name, desc in INTERNAL_KINDS
                     if re.search(pat, rel.replace(os.sep, "/"))), None)
        if kind:                                     # one line per kind, below
            k = kinds.setdefault(kind, [0, 0])
            k[0] += 1
            k[1] += nbytes or 0
            continue
        cat, desc = describe(os.path.basename(f))
        size = human(nbytes) if nbytes is not None else "?"
        by_cat.setdefault(cat, []).append((rel, size, desc))
        rows.append({"file": rel, "category": cat, "size": size, "description": desc})
    for (name, desc), (n, nbytes) in sorted(kinds.items()):
        size = f"{human(nbytes)} in {n} file{'s' if n != 1 else ''}"
        by_cat.setdefault("Search internals", []).append((name, size, desc))
        rows.append({"file": name, "category": "Search internals", "size": size,
                     "description": desc, "n_files": n})

    n_files = sum(r.get("n_files", 1) for r in rows)
    lines = ["# Output files — what each one is", "",
             f"This run produced {n_files} file(s). Each is described below, grouped by purpose"
             + (" (the search's per-run working files one line per kind)" if kinds else "")
             + ".", ""]
    n_unknown = 0
    # every category that has files, the known ones in order: a category missing from the order
    # (Figures was) must never drop its files from the catalog
    for cat in CATEGORY_ORDER + sorted(set(by_cat) - set(CATEGORY_ORDER)):
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
    lines.append("- **`Analysis_Report.html`** — read this first: the report (QC, figures and the "
                 "interpretation in one page; the text alone is `AI_Analysis_Report.md`).")
    lines.append("- **`de_results/DE_*.csv`** — the differentially expressed proteins per comparison.")
    lines.append("- **`de_results/methods.txt`** — the Methods paragraph for your paper.")
    lines.append("- **`de_results/reproducibility_log.R`** — the analysis as plain R. "
                 "Read it to see what was done; `Rscript` it to redo it with only R + limpa/limma.")
    lines.append("- **`reproducibility/REPRODUCE.md`** — the full bundle, for reproducing the "
                 "*search* too (pinned engine, environment lock, checksums).")
    lines.append("")

    with open(a.out, "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print(json.dumps({"report": os.path.abspath(a.out), "n_files": n_files,
                      "categories": {c: len(by_cat.get(c, [])) for c in CATEGORY_ORDER if by_cat.get(c)},
                      "unrecognized": n_unknown}, indent=2))


if __name__ == "__main__":
    main()
