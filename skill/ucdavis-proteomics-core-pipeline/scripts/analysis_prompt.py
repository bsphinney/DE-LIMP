#!/usr/bin/env python3
"""
analysis_prompt.py  --  Build the analysis brief the agent follows to interpret
the results. Faithful, complete port of DE-LIMP's "Export Prompt for Claude"
(R/server_ai.R), so the report is as thorough as the DE-LIMP app's AI export.

In DE-LIMP the prompt is shipped to an external Claude/Gemini. Here the agent IS
Claude (Claude Code / Desktop), so this writes ANALYSIS_PROMPT.md and the agent
reads it + the data files + the figures and writes a complete, figure-rich,
expert AI_Analysis_Report.md (then make_analysis_html.py renders it into
Analysis_Report.html, the report of record; no Word copy of the report is made).

The pipeline self-description and the educational/stats sections are chosen from
the actual engine + DE method (from de_provenance.json), never hardcoded (rule #1).

Usage:
  python3 analysis_prompt.py --out ANALYSIS_PROMPT.md \
      --de-dir output/tables --report output/search/report.parquet \
      --conditions input/conditions.csv --figures-dir output/figures \
      [--qc QC_Metrics.csv] [--gsea GSEA_Results.csv] \
      --engine diann --acquisition DIA --instrument "Orbitrap Astral" \
      --workflow-manifest input/wf/workflow.manifest.json --session <session dir>
"""
import sys, os, json, glob, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# ONE test for "the matrix is complete by construction", ONE contrast display form, ONE
# direction sentence and ONE default tag, shared with the HTML report.
from make_analysis_html import (matrix_complete, SUPPRESS_WHEN_COMPLETE,  # noqa: E402
                                contrast_label, LOGFC_DIRECTION, DEFAULT_TAG, background_flag,
                                _make_names, read_figures_json, APPENDIX, APPENDIX_FIGURES,
                                load_record, read_text, csv_text, _UNREADABLE, gene_label)
# what the DE columns mean -- session_docs' own (which groups are pull-down controls is
# experiment_type.pulldown_design's).
from session_docs import COLUMNS as DE_COLUMNS  # noqa: E402
import csv  # noqa: E402

# Abundant ER proteins that ride along in almost any membrane pull-down: background, not
# interactors (HGNC symbols; the mouse genes are the same names in title case).
ER_BACKGROUND = ("HSPA5", "HSP90B1", "CANX", "CALR", "P4HB", "PDIA3", "PDIA4", "PDIA6")
# The names of run_de.R's `Evidence` categories (its EVIDENCE_LABELS: every run of every group
# measured / never measured in a group / otherwise, "partly" + the record's zero_means word).
# They name the tiers; what they mean is quoted from de_provenance.json when the run recorded it.
EVIDENCE = ("measured in both", "presence call", "partly")


def submission_brief(w, source):
    """The CoreOmics submission, quoted for the writer: the report must describe the samples in
    the submitter's words. A report once turned the record's "cross-linked" into "chemically
    cross-linked"; every added word like that is a claim about the lab's samples nobody made."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import submission_report as sr
    rec, session = sr.resolve(source)
    if rec is None:
        return
    who, why = sr.prepared_by(rec)
    w(f"## The submission ({sr.label(rec)}) — quote it, add nothing")
    w("The HTML report opens with this record (make_analysis_html.py renders it from the "
      "session). Wherever you describe the samples, how they were prepared or what the study "
      "is, use the submission's own words and add no detail it does not state (if it says "
      "\"cross-linked\", do not write \"chemically cross-linked\"). Never copy a contact or "
      "billing detail into the report.")
    w({"lab": "- **The submitting lab prepared the samples"
              + (" and sent peptides" if sr.sent_as_peptides(rec) else "") + ".** Any Sample "
              "Preparation note must say so; do not describe extraction, reduction, alkylation "
              "or digestion as work the Core did.",
       "core": "- **The Core prepared the samples.** Their protocol belongs in the Methods "
               "(methods.md); do not invent one here."}.get(
        who, f"- **Who prepared the samples is not clear** ({why}). Say so; do not guess."))
    notes = sr.quality_notes(rec, session)
    if notes:
        w("- The HTML report's Submission section already lists these notes, so do not copy "
          "them. In **Data Quality Notes**, say what each one means for THESE results (e.g. "
          "whether paired samples were analysed as independent):")
        for n in notes:
            w(f"  - {n['text']}")
    w("")
    w(sr.render_markdown(rec, (), heading="### The record, as submitted"))


def sample_type_brief(w, session):
    """What the samples are, as the user stated it (sample_type.py, the session's one record) --
    or NOT STATED, and never to be guessed. A draft report once named a tissue from an
    abbreviation in the file names, and the data did not support it (a staff report, 2.10)."""
    import sample_type as st
    try:
        rec = st.load(session) if session else None
    except (OSError, ValueError) as e:
        print(f"[analysis_prompt] the sample-type record could not be read: {e}", file=sys.stderr)
        rec, unreadable = None, str(e)
    else:
        unreadable = None
    w("## What the samples are — as the user stated it")
    w(f"- **Sample type: {st.describe(rec)}**"
      + (f" (its record could not be read: {unreadable})" if unreadable else ""))
    if rec and rec.get("stated"):
        w("- Describe the samples in these words, adding nothing they do not state (a tissue "
          "region, a cell type, a disease state).")
        w("- Check it against the data: if the most abundant proteins or the sample-quality "
          "panels (SAMPLE_QUALITY.md) point to a different tissue or biofluid (e.g. muscle, "
          "plasma or skin proteins dominating a sample stated as brain), say so in **Data "
          "Quality Notes** as a question for the user; do not resolve it either way yourself.")
    else:
        w("- Do not name a tissue, cell type or biofluid anywhere in the report or Methods; "
          f"{st.NEVER_INFER}. Say in **Data Quality Notes** that the sample type was not stated, "
          "and describe the samples by their groups only. If the data suggest one (e.g. muscle "
          "proteins dominating), you may say what the proteins suggest, as an observation, "
          "never as the sample's identity.")
    w("")


def normalisation_brief(w, prov, etype, etype_err, design):
    """The experiment type, the quantities the DE read, and the normalisation check -- quoted
    from run_de.R's record (de_provenance.json normalization_check), never re-judged here."""
    nc = prov.get("normalization_check") if isinstance(prov.get("normalization_check"), dict) else {}
    w("## Experiment type and normalisation — as recorded")
    w(f"- **Experiment type:** " + (f"{etype['label']} ({etype['source']})" if etype else
                                    "NOT RECORDED" + (f" (its record could not be read: {etype_err})"
                                                      if etype_err else " -- do not assume one")))
    w(f"- **Pull-down design:** {'yes' if design['pulldown'] else 'no'} ({design['source']})"
      + (f". **Raise with the user:** {design['note']}" if design.get("note") else ""))
    w(f"- **Quantities the DE read:** {prov.get('normalisation') or 'NOT RECORDED'}")
    if nc.get("status") != "decided":
        w("- **Normalisation check: NOT RUN.** Say in Data Quality Notes that whether "
          "normalisation suited this experiment was not checked.")
    else:
        d = nc.get("decision") or {}
        w(f"- **Normalisation check:** default {nc.get('default', {}).get('quantities')} for this "
          "type; the data check " + ("agreed." if not nc.get("tripped") else
                                    "did NOT agree (" + "; ".join(nc.get("trips") or []) +
                                    f"); {d.get('by')} chose {d.get('quantities')}: {d.get('reason')}."
                                    " Say so in Data Quality Notes."))
        vol = (nc.get("volcano") or {}).get(nc.get("quantities_applied") or "") or {}
        for c, v in (vol.get("contrasts") or {}).items():
            if "n" not in v:
                continue
            ctr = "n/a" if v.get("centre_nonsig") is None else f"{v['centre_nonsig']:+.2f}"
            w(f"  - Volcano {c}: centre (non-significant) {ctr} log2, "
              f"{v['up']} up : {v['down']} down"
              + (f", bait rank {v['bait']['rank']}" if (v.get("bait") or {}).get("rank") else "")
              + ". Caption this contrast's volcano with the centre shift and the up:down ratio.")
    w("")


def _detection(de_dir, conditions):
    """-> (per-sample detection {protein: {run: bool}}, groups {group: [runs]}) from
    Detection_Matrix.csv + conditions.csv, or ({}, {}) when either is missing."""
    path = os.path.join(de_dir, "Detection_Matrix.csv")
    det, groups = {}, {}
    if not (os.path.exists(path) and conditions and os.path.exists(conditions)):
        return det, groups
    rd = csv.reader(csv_text(path))
    head = next(rd)
    for rec in rd:
        det[rec[0]] = {head[i]: _on(rec[i]) for i in range(1, len(rec))}
    for r in csv.DictReader(csv_text(conditions)):
        groups.setdefault((r.get("Group") or "").strip(), []).append(
            (r.get("File.Name") or "").strip())
    return det, groups


def _on(v):
    try:
        return float(v) > 0
    except (TypeError, ValueError):
        return False


def _k_of_n(row, group, det, groups):
    """(k measured, n) for one protein in one group: the DE table's own Detected_<group>
    column when present, else Detection_Matrix.csv. None when neither can say."""
    v = row.get(f"Detected_{group}")
    if v and "/" in v:
        k, n = v.split("/", 1)
        try:
            return int(k), int(n)
        except ValueError:
            pass
    d = det.get(row.get("Protein.Group"))
    if d is not None and group in groups:
        return sum(d.get(s, False) for s in groups[group]), len(groups[group])
    return None


def never_in_control(de_dir, de_file, contrast, control, adjp, det, groups, k=15):
    """Significant, enriched proteins of a bait-vs-control contrast that were never measured
    in the control -> (count, [gene names by adj.P]), or None when detection is unknown."""
    rows, known = [], False
    with csv_text(os.path.join(de_dir, de_file)) as fh:
        for r in csv.DictReader(fh):
            try:
                p, lfc = float(r.get("adj.P.Val")), float(r.get("logFC"))
            except (TypeError, ValueError):
                continue
            kn = _k_of_n(r, control, det, groups)
            known = known or kn is not None
            if p < adjp and lfc > 0 and kn is not None and kn[0] == 0:
                if not background_flag(r.get("Genes"), r.get("Protein.Group"),
                                       r.get("Protein.Names")):
                    rows.append((p, gene_label(r)))
    if not known:
        return None
    rows.sort()
    return len(rows), [g for _, g in rows[:k]]


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
    ap.add_argument("--submission", help="the session dir: its CoreOmics submission "
                                         "(submission_report.py) is quoted in the brief")
    ap.add_argument("--session", help="the session dir: its sample-type record "
                                      "(sample_type.py) is quoted in the brief. Default: "
                                      "--submission when that is a session dir")
    a = ap.parse_args()

    # the page's loaders: UTF-8 whatever the locale, and a record that exists but cannot be
    # read is said (stderr + a line in this brief), never a silent {}
    prov = load_record(os.path.join(a.de_dir, "de_provenance.json"))
    wfman = load_record(a.workflow_manifest) if a.workflow_manifest else {}
    method = prov.get("method") or (wfman.get("de", {}) or {}).get("method", "")
    engine = a.engine or (wfman.get("engine", {}) or {}).get("name", "")
    eng_ver = (wfman.get("engine", {}) or {}).get("version", "")
    acq = (a.acquisition or wfman.get("acquisition", "")).upper()
    de_files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(a.de_dir, "DE_*.csv")))
    contrasts = prov.get("contrasts") or []

    def recorded(key):
        v = prov.get(key)
        return v if isinstance(v, (int, float)) and not isinstance(v, bool) else None
    # As recorded, or said to be unrecorded (rule 2): never a silent 0.01 / 1.0 / 0.05.
    q_cut, lfc, adjp_rec = recorded("q_cutoff"), recorded("logfc"), recorded("adjp")
    adjp = adjp_rec if adjp_rec is not None else 0.05
    adjp_txt = f"{adjp:g}" if adjp_rec is not None else f"{adjp:g} {DEFAULT_TAG}"
    groups_all = list((prov.get("groups") or {}).keys()) if isinstance(prov.get("groups"), dict) else []
    # THE pull-down decision (experiment_type.pulldown_design): from the session's experiment
    # type when one is recorded -- the same record that set the DE's quantities -- else from the
    # group names, said so. The brief and the DE cannot disagree about the design.
    session = a.session or (a.submission if a.submission and os.path.isdir(
        os.path.join(a.submission, "input")) else None)
    import experiment_type as etm
    try:
        etype, etype_err = etm.load(session), None
    except (OSError, ValueError) as e:
        etype, etype_err = None, str(e)
        print(f"[analysis_prompt] the experiment-type record could not be read: {e}",
              file=sys.stderr)
    design = etm.pulldown_design(etype, groups_all, contrasts)
    controls, vs_control, pulldown = design["controls"], design["vs_control"], design["pulldown"]
    det, det_groups = _detection(a.de_dir, a.conditions)
    has_detmat = os.path.exists(os.path.join(a.de_dir, "Detection_Matrix.csv"))
    zero_word = ("missing" if (prov.get("detection_matrix") or {}).get("zero_means") == "missing"
                 else "inferred")

    desc_line = prov.get("display_label") or "the configured quantification + limma pipeline"
    rollup = prov.get("rollup_method", "")
    de_engine = prov.get("de_engine", "")
    missing_policy = prov.get("missing_policy", "")
    citation = prov.get("citation", "")

    # figures: make_figures.R's figures.json, through the reader the HTML report uses (both
    # the 2.8.0 object {"figures", "failed"} and an older session's bare list)
    figs, failed, fig_error = [], [], None
    fdir_rel = ""
    if a.figures_dir and os.path.isdir(a.figures_dir):
        fjd = read_figures_json(a.figures_dir)
        figs, failed = fjd["figures"], fjd["failed"]
        fig_error = fjd["error"] or (None if fjd["found"] else "not found")
        if fig_error:
            print(f"[analysis_prompt] WARNING: {a.figures_dir}/figures.json {fig_error}: the "
                  f"brief lists no figures", file=sys.stderr)
        fdir_rel = os.path.basename(a.figures_dir.rstrip("/")) or "figures"

    # A complete-by-construction matrix (DPC/limpa) makes "proteins quantified per sample" a
    # row of identical bars: never ask the writer to embed or describe it -- make_analysis_html
    # drops it, and a paragraph written about it, anyway (the same test, one definition).
    complete, complete_why = matrix_complete(a.de_dir, prov)
    if complete:
        figs = [f for f in figs if not str(f.get("file", "")).startswith(SUPPRESS_WHEN_COMPLETE)]
    # The p-value panel is a calibration check: the HTML report adds it as its appendix, so the
    # narrative does not embed it (or interleave per-contrast p-value plots with the findings).
    appendix = [f for f in figs if f["file"].startswith(APPENDIX_FIGURES)]
    embed = [f for f in figs if f not in appendix]
    legacy_pv = [f for f in embed if f["file"].startswith("pvalue_")]

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
    # run_de.R's contaminant record -- stated as recorded, never assumed.
    cont = prov.get("contaminants") if isinstance(prov.get("contaminants"), dict) else {}
    if cont.get("policy") == "removed":
        from make_methods import _count_of          # the record's unit, singular for one
        n_items, were = _count_of(cont.get("n_precursors"), cont)
        w(f"- Contaminants: {n_items} mapping to a {cont.get('tag')} "
          f"entry {were} removed before quantification ({cont.get('n_protein_groups')} contaminant "
          f"protein groups); they are NOT in the DE tables or the expression matrix.")
    elif cont.get("policy") == "kept":
        w(f"- Contaminants: KEPT (--keep-contaminants) — {cont.get('n_protein_groups')} "
          f"{cont.get('tag')} protein groups are in the DE tables; call any that are "
          f"significant contamination, not biology.")
    elif prov and not cont:
        w("- Contaminants: not recorded by this DE run (older run_de.R, which did NOT remove "
          "them) — any `Cont_` protein in the DE tables is contamination, not biology.")
    if cont.get("database_risk") is True:
        w(f"- **Contaminant-filter caveat (say this in Data Quality Notes, naming the proteins):** "
          f"{cont.get('database_note')}")
    # A keratin sample (hair, wool, feather, skin, nail): make_methods' sentence is its one
    # description. Keratins are the analyte there -- never call them contamination.
    from make_methods import de_keratin_sentence
    ker = de_keratin_sentence(cont)
    if ker and (cont.get("keratin_sample") or {}).get("value") is not True:
        w(f"- Keratin: {ker} Say so in Data Quality Notes if keratins matter to the question.")
    elif ker:
        kn = (cont.get("keratin_sample") or {}).get("note")
        w(f"- Keratin sample: {ker} Keratins in the results are the ANALYTE -- do not call them "
          f"contamination." + (f" Caveat (say it in Data Quality Notes): {kn}"
                               if (cont.get("keratin_sample") or {}).get("caution") is True
                               and kn else ""))
    # run_de.R's block record (--block): make_methods' sentence is its one description.
    blk = prov.get("block") if isinstance(prov.get("block"), dict) else {}
    if blk.get("applied") is True:
        from make_methods import de_block_sentence
        w(f"- Blocking: {de_block_sentence(prov)}")
        for warn in blk.get("warnings") or []:
            w(f"- **Blocking caveat (say this in Data Quality Notes):** {warn}")
    elif blk.get("column") and blk.get("note"):
        w(f"- Blocking: none used — {blk['note']}")
    if adjp_rec is not None:
        w(f"- Significance rule: adj.P.Val < {adjp:g} (Benjamini-Hochberg within each contrast) "
          f"— the adjusted p-value ALONE (identification FDR: "
          + (f"q ≤ {q_cut:g}" if q_cut is not None else "not recorded") + "). "
          "**No fold-change filter is applied.**"
          + (f" |log2FC| = {lfc:g} ({2 ** float(lfc):.3g}-fold) is drawn on the volcano as a "
             f"reference line only." if lfc is not None else ""))
        w("- Do not re-impose a fold-change cutoff when you count or describe significant "
          "proteins, and never call a result significant \"because it is ≥2-fold\". The test "
          "was of H0: fold change = 0, so that is the claim the error rate covers. You may say "
          "how many significant proteins also exceed the reference line — as a description of "
          "effect size, not as a second significance criterion.")
    else:
        w(f"- Significance rule: **not recorded** (de_provenance.json has no adjp). Count at "
          f"adj.P.Val < {adjp_txt} and say in the report that this cutoff is a default, not the "
          f"rule the DE applied. Do not claim that no fold-change filter was used unless "
          f"`methods.txt` says so.")
    w(f"- Contrasts — name each exactly this way everywhere (text, tables, figure notes): "
      + "; ".join(f"`{c}` = **{contrast_label(c)}**" for c in contrasts)
      + f". {LOGFC_DIRECTION}" if contrasts else f"- {LOGFC_DIRECTION}")
    w("- The exact, self-describing methods text is in `tables/methods.txt` — do not "
      "contradict it or invent a different pipeline.")
    # qc_bracket.py's verdict is staff-only (Brett, 2026-10-01): the writer of a client report
    # must not carry it over from logs/ or from what staff were told
    w("- The instrument's QC check around this project (`logs/qc_bracket.*`) is STAFF-ONLY: "
      "do not mention its verdict, the QC runs, IPS or that check in this report.")
    w("")

    # the sample type, as the user stated it -- or NOT STATED, which forbids a guess
    sample_type_brief(w, session)
    normalisation_brief(w, prov, etype, etype_err, design)
    if a.submission:
        submission_brief(w, a.submission)

    unread_at = len(L)                       # records that could not be read: listed here, last
    w("## Attached data files — read these")
    heads = {}                               # each DE table's own header: its own groups
    for f in de_files:
        try:
            heads[f] = next(csv.reader(csv_text(os.path.join(a.de_dir, f))))
        except (OSError, StopIteration):
            heads[f] = []
    de_cols = [c for f in de_files for c in heads[f]]
    det_cols = [c for c in de_cols if c.startswith("Detected_")]
    has_evidence = "Evidence" in de_cols
    # run_de.R's own descriptions of Detected_<group> / Evidence (de_provenance.json): read,
    # never restated -- session_docs' AGENTS.md quotes the same record.
    dm_rec = prov.get("detection_matrix") if isinstance(prov.get("detection_matrix"), dict) else {}
    de_col_desc = dm_rec.get("de_columns") if isinstance(dm_rec.get("de_columns"), dict) else {}
    for f in de_files:
        own = [c for c in heads[f] if c.startswith("Detected_")]
        w(f"- `tables/{f}` — DE results for one comparison: Protein.Group, logFC, "
          "AveExpr, t, P.Value, adj.P.Val (BH within the contrast), B, gene annotation"
          + (", PropObs/NPeptides (limpa)" if "PropObs" in heads[f] else "")
          + (f", per-group detection {', '.join(f'`{c}`' for c in own)} (k/n samples measured)"
             if own else "")
          + (", `Evidence`" if "Evidence" in heads[f] else "") + ".")
    for col, desc in (("Detected_<group>", de_col_desc.get("Detected_<group>")),
                      ("Evidence", de_col_desc.get("Evidence"))):
        if desc and (det_cols if col.startswith("Detected_") else has_evidence):
            w(f"  - `{col}`: {desc} (as recorded in `de_provenance.json`).")
    w("- `tables/Expression_Matrix.csv` — log2 protein abundance per sample"
      + (" (complete: every cell has a value, measured OR inferred)." if complete else "."))
    if has_detmat:
        w(f"- `tables/Detection_Matrix.csv` — per protein x sample: precursors observed; "
          f"**0 = {zero_word}** (not measured in that sample). Same rows and columns as the "
          f"expression matrix. This, not the expression matrix, says what was measured.")
    if os.path.exists(os.path.join(a.de_dir, "QC_detected_vs_inferred.csv")):
        w("- `tables/QC_detected_vs_inferred.csv` — per sample: proteins measured vs inferred.")
    w("- `tables/methods.txt`, `tables/de_provenance.json` — methods + exact versions.")
    if cont.get("share_table"):
        w(f"- `tables/{cont['share_table']}` — per-run contaminant share of "
          f"{cont.get('intensity_column')} (QC): report a high or group-confounded share in "
          f"Data Quality Notes.")
    if cont.get("removed_table"):
        w(f"- `tables/{cont['removed_table']}` — the protein groups the contaminant filter "
          f"removed (or trimmed) before quantification.")
    if a.conditions:
        w(f"- `{a.conditions}` — experimental design (File.Name → Group [+ Batch/Covariates]).")
    if has_qc:
        w(f"- `{a.qc}` — per-sample QC metrics.")
    if has_gsea:
        w(f"- `{a.gsea}` — Gene Set Enrichment results (gene-permutation GSEA: exploratory, see "
          "below).")
    # run_sets.R: its record, its tables, its own reading rules -- quoted, not restated
    sets_rec = load_record(os.path.join(a.de_dir, "sets_provenance.json"))
    if sets_rec:
        w("- `tables/Sets_<contrast>.csv` — protein-set tests per comparison (camera + fry, the "
          "run-depth check, presence-call fractions); `tables/sets_provenance.json` — their "
          "record and methods" + ("; `tables/Sets_baitnorm_<contrast>.csv` — the pulldown's "
                                  "tests relative to each bait's complex"
                                  if (sets_rec.get("pulldown") or {}).get("baits") else "") + ".")
        cols = sets_rec.get("columns") if isinstance(sets_rec.get("columns"), dict) else {}
        for col, desc in cols.items():
            w(f"  - `{col}`: {desc}")
    elif (prov.get("set_tests") or {}).get("suggested_for"):
        w(f"- Protein-set tests were suggested by run_de for "
          f"{', '.join(prov['set_tests']['suggested_for'])} (few significant proteins) but have "
          "not been run: say so if the report discusses pathways, and do not stand in for them "
          "by eye.")
    w("")
    w("Compute everything you cite directly from these files:")
    w(f"- significant proteins per contrast (adj.P.Val < {adjp_txt} only) and the up/down split;")
    w("- the top up/down proteins per contrast **ranked by adj.P.Val (ties by P.Value) — never "
      "by |logFC|**: where a group was never measured, the largest fold changes are detection "
      "events set by the missing-value model, not the strongest effects;")
    w("- proteins significant in ≥2 contrasts;")
    w("- reproducibility (CV) per group — on the **linear** scale (2^value), from **measured "
      "values only** (Detection_Matrix.csv > 0 for that protein and sample); skip a group with "
      "fewer than 2 measured values, and call a protein with no measured value in a group "
      "*all-inferred* there. Never compute a CV from inferred values: a group never measured "
      "gets near-identical model values, which reads as CV ≈ 0 — \"most stable\" — when nothing "
      "was measured"
      + ("." if has_detmat else ". **Detection_Matrix.csv is not attached, so you cannot tell "
                                "measured from inferred: do not report CVs as reproducibility.**"))
    w("Use gene names where available; a protein with no gene name goes by its UniProt entry "
      "name (`Protein.Names`, e.g. HVM51_MOUSE), not its accession, as the HTML report's "
      "tables do. **Never fabricate a value, protein, pathway, or citation you cannot ground "
      "in the data or established biology.**")
    w("")
    # measured vs inferred: how to tier, and what PropObs is NOT
    w("## Measured vs inferred — how to tier hits")
    if "PropObs" in de_cols:
        w(f"- `PropObs`: {DE_COLUMNS['PropObs']}. It is limpa's rowMeans(n.observations) / "
          "NPeptides over **all runs of the study** — study-wide, not per comparison. **Never "
          "tier hits by PropObs and never describe it as detection in either group of a "
          "comparison.**")
    ev_all, ev_none, ev_part = EVIDENCE[0], EVIDENCE[1], f"{EVIDENCE[2]} {zero_word}"
    source = ("the DE tables' `Evidence` column" if has_evidence else
              "the DE tables' `Detected_<group>` columns (k/n samples measured)" if det_cols else
              "Detection_Matrix.csv + conditions.csv (count each group's runs with the protein "
              "measured, > 0)" if has_detmat else None)
    recorded_ev = has_evidence and bool(de_col_desc.get("Evidence"))
    if source is None:
        w("- Per-group detection is not available in this run — say so rather than guessing, "
          "and do not sort hits into measured and detection-event tiers.")
    else:
        w(f"- Tier each hit by **per-group detection**, from {source}"
          + ("" if has_evidence else
             ", sorted into the three categories run_de.R writes as `Evidence`") + ":")
        w(f"  - **{ev_all}**"
          + ("" if recorded_ev else " (measured in every run of every group in the contrast)")
          + ": the fold change is a measured ratio;")
        w(f"  - **{ev_none}**"
          + ("" if recorded_ev else " (never measured in at least one group)")
          + ": a detection event — present vs absent. Report presence, not the size of the fold "
            f"change (the absent side is {zero_word});")
        w(f"  - **{ev_part}**" + ("" if recorded_ev else " (every other case)")
          + ": a measured ratio with gaps — give each group's k/n.")
    w("")
    if pulldown:
        w("## Pull-down design — enrichment over a control")
        w(f"Control group(s), by name: {', '.join(controls)}. Contrasts against them "
          f"({', '.join(contrast_label(c) for c in vs_control)}) measure **enrichment over the "
          "control**, not abundance change.")
        w("- Do NOT write a blanket sentence such as \"the IgG value is a model-inferred floor\". "
          "Name the proteins that were never measured in the control, per contrast:")
        for c in vs_control:
            ctrl = c.split("-", 1)[1]
            f = next((x for x in de_files if x.endswith(f"_{_make_names(c)}.csv")), None)
            res = never_in_control(a.de_dir, f, c, ctrl, adjp, det, det_groups) if f else None
            w(f"  - {contrast_label(c)}: "
              + ("detection in the control not recorded — say so." if res is None else
                 f"{res[0]} significant enriched protein{'' if res[0] == 1 else 's'} never "
                 f"measured in {contrast_label(ctrl)}"
                 + (f" (top by adj.P: {', '.join(res[1])})" if res[1] else "") + "."))
        w("- Background, not interactors: the antibody's own Ig heavy/light chains, "
          "contaminant (`Cont_`) entries, and abundant ER proteins (" + ", ".join(ER_BACKGROUND)
          + "). The HTML report's top-protein tables flag Ig chains and contaminants.")
        w("")

    # ---- figures ----
    if figs or failed or fig_error:
        w("## Figures — embed and interpret each one")
    if fig_error:
        w(f"**No figure list:** `{fdir_rel}/figures.json` {fig_error}, so no figures are listed "
          "here. Do not embed image files by guessing their names; say in QC Assessment that the "
          "figures are not available.")
        w("")
    if embed:
        w(f"These publication-quality figures were generated for you in `{fdir_rel}/`. "
          "**Embed every figure** in the relevant section using markdown image syntax "
          f"(e.g. `![caption]({fdir_rel}/<file>)`) so it renders in the HTML report, "
          "and **write an expert interpretation of what each shows for THIS dataset** — "
          "not a generic caption. Available figures:")
        for fig in embed:
            w(f"- `{fdir_rel}/{fig['file']}` ({fig['type']}) — {fig['caption']}")
        w("")
        w("Placement: each volcano in **Key Findings Per Comparison**; PCA + "
          + ("detected-vs-inferred" if complete else "per-sample counts")
          + " in **QC Assessment**; the heatmap in "
          + ("**Specificity**" if pulldown else "**Cross-Comparison Biomarkers**")
          + " or **Biological Interpretation**.")
        if any(f["type"] == "violin" for f in embed):
            w(f"Place each `{fdir_rel}/violin_top_<contrast>.png` in that comparison's section, "
              "right after its volcano. Call it the **top-protein plot**: a group with few runs "
              "is drawn as its points and a mean bar, not a violin. Add one line on how to read "
              "it, in its caption's words (filled points were measured in that run; the caption "
              "says how unmeasured values are shown). Name every protein its subtitle flags as "
              "never measured in one group: that fold change is a detection event (seen in one "
              "group, not the other), not a measured magnitude, so do not quote its size as an "
              "effect. If the figure says detection status was not recorded, say so; never "
              "describe those points as measured.")
        if legacy_pv:
            w(f"The per-contrast p-value histograms (`{fdir_rel}/pvalue_<contrast>.png`) are a "
              "calibration check, not findings: embed them together in a final "
              f"`## {APPENDIX}` section, never in the per-comparison sections.")
        if complete:
            w(f"Do NOT embed, reference or describe a proteins-quantified-per-sample plot "
              f"(`qc_protein_counts.png`): {complete_why}, so every sample shows the same "
              f"count and the plot says nothing. Per-sample depth is the *detected* part of "
              f"`qc_detected_vs_inferred.png`.")
        w("")
    for fig in appendix:
        w(f"**Appendix figure — do not embed:** `{fdir_rel}/{fig['file']}` — {fig['caption']} "
          f"The HTML report adds it as its last section (*{APPENDIX}*). In **QC Assessment**, "
          "name any contrast whose histogram is not flat with a peak near 0 (skewed toward 1, "
          "a hump mid-range or a U shape) and point the reader to the appendix; if all look "
          "well calibrated, say so in one sentence.")
        w("")
    if failed:
        w("**Not drawn in this run** (make_figures.R's `failed` list) — do not embed or refer "
          "to these as images. Where the report would have used one, say in one line that it "
          "is not available and why (the HTML report also lists them at the top):")
        for f in failed:
            w(f"- `{fdir_rel}/{f['file']}`" + (f" ({f['type']})" if f["type"] else "")
              + (f" — {f['reason']}" if f["reason"] else " — no reason recorded"))
        w("")

    w("## OUTPUT — write `" + a.report_out + "` with ALL of these sections (markdown)")
    w("Use markdown headers, tables for numbers, and embed the figures. Be scientific "
      "but accessible. Reference specific proteins/genes throughout.")
    w("")
    w("### Overview")
    w("Number of comparisons analyzed and the overall picture. The HTML report's **Results at "
      "a glance** already shows significant / up / down per contrast (counted from the DE "
      "tables) — **do not add a second per-contrast table**; refer to it. Overall assessment of "
      "the experiment's quality, depth (proteins measured, not merely quantified), and scope.")
    w("")
    w("### QC Assessment")
    w("Evaluate technical quality. Use the PCA and "
      + ("detected-vs-inferred (per-sample depth)" if complete else "per-sample protein-count")
      + " figures, and "
      "the QC metrics if present. Comment on consistency of identifications across "
      "replicates and groups, whether replicates cluster, and flag any outlier samples or "
      "systematic biases (e.g. a group quantifying far fewer proteins). State clearly "
      "whether the data are fit for differential analysis.")
    w("")
    w("### Key Findings Per Comparison")
    w("For each comparison: embed its volcano plot, then highlight the top up- and "
      "down-regulated proteins **ranked by adj.P** (use gene names; give logFC, adj.P and each "
      "group's detection, k/n measured). Say which are detection events. Note any comparison "
      "with unusually few or many hits"
      + (f", and name it if its p-value histogram (*{APPENDIX}*) looks miscalibrated — no "
         "p-value figure in this section." if appendix or legacy_pv else "."))
    w("")
    if pulldown:
        w("### Specificity: bait-specific, shared, background")
        w("Sort the enriched proteins into three groups, naming them: **bait-specific** "
          "(enriched over the control for one bait only), **shared across baits** (enriched for "
          "two or more — complex partners, or common background), and **background** (Ig "
          "chains, contaminants, abundant ER proteins). For each bait-specific or shared "
          "protein, give its detection in the bait and control groups (k/n measured). Embed the "
          "heatmap here if there is one.")
        w("")
    else:
        w("### Cross-Comparison Biomarkers")
        w("Proteins significant in several comparisons are worth listing, with their direction "
          "in each (always up, always down, or mixed) — but being significant in many "
          "comparisons is **not** by itself higher confidence: comparisons that share a group "
          "(e.g. several treatments against one control) are not independent, so one noisy "
          "shared group can make a protein \"significant everywhere\". Embed the heatmap and use "
          "it to show which proteins drive the separation. (If there is only one comparison, say "
          "so.)")
        w("")
        w("### High-Confidence Findings")
        w("Pick the strongest findings by the tier rule above — significant AND measured in "
          "both groups first — not by fold change. Where you quote reproducibility, use the CV "
          "from measured values only. Discuss known functions and pathway involvement of the "
          "genes you recognize, and state for each whether it is a measured change or a "
          "detection event.")
        w("")
    if has_gsea:
        w("### Pathway & Gene Set Enrichment Analysis")
        # gseGO / fgsea permute GENES: they treat a set's proteins as independent -- the flaw
        # that keeps limma's geneSetTest out of run_sets.R (Dog CSF, 2026-09-21)
        w("The GSEA file comes from gene-permutation GSEA (gseGO / fgsea), which treats a set's "
          "proteins as independent: co-regulated proteins are not, so its p-values are too small "
          "for exactly the sets that matter.")
        if sets_rec:
            w("Do not summarise it as findings, and do not put it beside the protein-set tests as "
              "if it were a second test: the protein-set tests (camera + fry, above) are the "
              "result. Mention the GSEA file at most as exploratory, with that caveat.")
        else:
            w("Summarise it only as exploratory, with that caveat stated: the top pathways by "
              "|NES|, each connected to the DE proteins it rests on. Never as a finding on its own.")
        w("")
    if sets_rec:
        fdr_rec = ((sets_rec.get("settings") or {}).get("fdr") or {}).get("value")
        w("### Protein-set tests")
        w("The HTML report adds its own 'Protein-set tests' section from `sets_provenance.json` "
          "(counts, the sets significant in both tests, the pulldown table, the methods). In your "
          "text, interpret; do not re-tabulate. Per comparison, from the `Sets_*.csv` tables "
          f"(significant = FDR < {fdr_rec if fdr_rec is not None else 'not recorded'}):")
        w("- start from the record's own reading of each comparison (below): use its words, and "
          "do not make it stronger — 'not detectable in this design' never becomes 'unchanged', a "
          "broad shift that cannot be separated from run depth never becomes a named signature, "
          "and where camera has little power its count of zero is not a finding;")
        w("- lead with sets significant in BOTH camera and fry that hold with run depth in the "
          "model, whose fry call holds on every reference basis and that are not mostly presence "
          "calls; name their proteins from the DE table;")
        w("- a one-test result is reported as such (which test, and what that means below), never "
          "as a pathway finding.")
        w("The record's reading of each comparison:")
        for cn, c in (sets_rec.get("contrasts") or {}).items():
            if c.get("reading"):
                w(f"- {contrast_label(cn)}: {c['reading']}")
        for b, r in ((sets_rec.get("pulldown") or {}).get("baits") or {}).items():
            for cn, c in (r.get("contrasts") or {}).items():
                if c.get("reading"):
                    w(f"- {b}, relative to its complex, {contrast_label(cn)}: {c['reading']}")
        if (sets_rec.get("pulldown") or {}).get("direction_bias"):
            w("Relative to a bait's complex, what is not found is not evidence of no loss: "
              + sets_rec["pulldown"]["direction_bias"] + " Never write that the complex lost nothing.")
        w("Read them by run_sets.R's own rules:")
        for rule in sets_rec.get("reading") or []:
            w(f"- {rule}")
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
    if prov.get("plain_language"):
        # the pipeline's own words for itself (build_maxlfq.R descriptor), never a guess from
        # the --method flag: a Sage DE runs --method maxlfq and is not MaxLFQ
        w(f"- {prov['plain_language']}")
    elif method == "dpc":
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
    w("  - **p-value**: if the protein did not truly change, the probability of a difference "
      "at least this large from measurement noise alone.")
    w("  - **adjusted p-value (FDR, Benjamini–Hochberg)**: p-values corrected for testing "
      "thousands of proteins at once — BH within each contrast; explain the multiple-testing "
      "problem with an intuitive coin-flip example.")
    w("  - **volcano plot**: why it's volcano-shaped and how to read it (x = effect size, "
      "y = significance).")
    w("  - **coefficient of variation (CV)**: SD / mean on the linear scale, from measured "
      "values only — lower is more reproducible; an inferred value is not a measurement.")
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
    w("")
    w("### Next steps")
    w("Close with what to do next, using this explicit tier rule (the detection categories "
      "above), and name the proteins in each tier (top few by adj.P):")
    w(f"1. **Follow up first** — significant and *{ev_all}*: a measured ratio.")
    w(f"2. **Check the counts** — significant and *{ev_part}*: give each group's k/n; those "
      "measured in most runs of both groups come first.")
    w(f"3. **Confirm presence** — significant with a *{ev_none}*: a detection event; confirm "
      "with an orthogonal method (western, targeted MS) before quoting a fold change"
      + (" — for a pull-down, these are the interactor candidates, those seen in most bait "
         "runs first." if pulldown else "."))
    w("Whatever its tier, a hit measured in fewer than half the runs of every group is a lead "
      "only.")
    w("Then point the reader to `README.html` at the top of the results folder for where every "
      "file is and how to reuse them.")
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
      "established biology. After writing the report, it is rendered into "
      "Analysis_Report.html (the report of record).")
    w("")

    if _UNREADABLE:
        L[unread_at:unread_at] = (
            ["## Records that could not be read — say so, do not fill in"]
            + [f"- `{p}` {why}" for p, why in _UNREADABLE.items()]
            + ["What these records hold (a cutoff, a contaminant caveat, a study fact) may be "
               "missing from this brief. Say in Data Quality Notes which could not be read; never "
               "describe what they would have said. The HTML report says so too.", ""])
    with open(a.out, "w", encoding="utf-8") as fh:        # "≥", "—": not cp1252-safe
        fh.write("\n".join(L) + "\n")

    print(json.dumps({
        "prompt": os.path.abspath(a.out),
        "report_to_write": a.report_out,
        "engine": engine, "engine_version": eng_ver, "acquisition": acq, "de_method": method,
        "de_files": de_files, "contrasts": contrasts,
        "figures": [f["file"] for f in figs], "n_figures": len(embed),
        "figures_appendix": [f["file"] for f in appendix],
        "figures_failed": [f["file"] for f in failed], "figures_error": fig_error,
        "has_qc": has_qc, "has_gsea": has_gsea, "records_unreadable": dict(_UNREADABLE),
        "next": f"Read {a.out} + the data files + figures, then write {a.report_out} "
                f"(ALL sections, embed all {len(embed)} figures, expert interpretation), "
                "then render it with make_analysis_html.py into Analysis_Report.html "
                "(no Word copy of the report).",
    }, indent=2))


if __name__ == "__main__":
    main()
