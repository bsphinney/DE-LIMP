"""
report_files.py -- the names of report files that more than one script needs, defined once.
No imports: make_report.py (OUTPUT_FILES.md) reads it without pulling in make_podcast.py and
its dependencies, so a partial copy of scripts/ still catalogues.

  * SHARE_NAME  -- the report with the podcast's audio and transcript built in, one file to
    send (make_podcast.py share writes it; core_submission.py deliver ships it);
  * SHARE_STALE -- the name an out-of-date copy is renamed to, so it is never sent
    (scratch_files.py leaves every *.stale.html out);
  * is_de_output -- the differential expression's output, by name (below).
"""
SHARE_NAME = "Analysis_Report_with_audio.html"
SHARE_STALE = "Analysis_Report_with_audio.stale.html"

# The differential expression's output, by name: what run_de.R (or, for SILAC, run_silac.R)
# writes into its --outdir (the
# session's output/tables/: the DE tables, the matrices it ran on, its record and methods, and
# contaminants.R's tables), what run_sets.R's protein-set tests on its model write beside them,
# and what make_figures.R and run_sets.R draw from them into output/figures/. A delivery without
# differential expression (core_submission.py deliver --mode no-de) refuses a session holding
# any of these: a DE goes out only in an analysis delivery, behind the normalisation check.
# tests/test_de_output_names.py reads the R scripts against this list.
DE_RECORD = "de_provenance.json"                  # run_de.R's record: a session with one has a DE
DE_OUTPUT_NAMES = frozenset((
    DE_RECORD, "Expression_Matrix.csv", "Detection_Matrix.csv", "methods.txt", "sessionInfo.txt",
    "reproducibility_log.R", "set_test_inputs.rds", "DE-LIMP_session.rds",
    "QC_detected_vs_inferred.csv", "QC_contaminant_share.csv", "contaminants_removed.csv",
    "sets_members.csv", "sets_methods.txt", "sets_provenance.json",
    "figures.json", "sample_labels.csv", "pca.png", "heatmap_top.png", "qc_pvalue_panel.png",
    "qc_protein_counts.png", "qc_detected_vs_inferred.png",
    # run_silac.R (a SILAC session's ratio analysis, in the same folders)
    "SILAC_channel_specific_detections.csv", "SILAC_channel_summary.csv", "SILAC_ratio_offsets.csv",
    "SILAC_label_incorporation.csv", "SILAC_sample_matrix.csv", "SILAC_samples.csv",
    "silac_channel_composition.png", "silac_channel_detections.png"))
# One file per contrast, (prefix, suffix): DE_<method>_<contrast>.csv, Sets_[baitnorm_]<contrast>.csv,
# volcano_ / violin_top_ / sets_[baitnorm_]<contrast>.png, and pvalue_<contrast>.png (before 2.8.0).
DE_OUTPUT_PATTERNS = (("DE_", ".csv"), ("Sets_", ".csv"), ("volcano_", ".png"),
                      ("violin_top_", ".png"), ("sets_", ".png"), ("pvalue_", ".png"),
                      ("silac_scatter_", ".png"))
DE_OUTPUT_DIRS = ("tables", "figures")            # output/<dir>/, where they are written


def is_de_output(rel):
    """`rel`, a path relative to the session's output/ (as deliver plans it: "tables/x.csv"),
    is the differential expression's output."""
    parts = rel.replace("\\", "/").strip("/").split("/")
    if len(parts) < 2 or parts[0] not in DE_OUTPUT_DIRS:
        return False
    name = parts[-1]
    return name in DE_OUTPUT_NAMES or any(name.startswith(p) and name.endswith(s)
                                          for p, s in DE_OUTPUT_PATTERNS)
