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
  numeric_covariate a covariate fitted as a number (run_de.R --numeric-covariates) that nearly
                    matches the groups of a comparison (de_provenance.json `covariates`: its
                    `confounded` flag, the comparisons in `confounded_contrasts`, and run_de.R's
                    own caution naming them, quoted) -- those comparisons are adjusted for it
                    (WARN); PASS for one that is not; nothing without a numeric covariate
  acquisition_mix   DIA and DDA mixed in one analysis (FAIL)
  instrument_mix    >1 instrument model in one analysis — batch effect (WARN)
  id_depth          suspiciously few proteins quantified (WARN)
  missingness       very high fraction of missing values (WARN)
  contamination     keratin/trypsin contaminants among the proteins (WARN)
  target_contaminants
                    quantified proteins whose sequence is ALSO a common-contaminant
                    entry (fetch_fasta.py removed those entries so the protein is
                    quantified under its own accession): possible contamination, KEPT
                    in quantification -- reported, never excluded (WARN)
  contaminant_overlap
                    real proteins present only as identical (or peptide-
                    indistinguishable) Cont_ entries in a
                    database used as-is -- excluded from quant, missing here (WARN)
  contaminant_evidence
                    (a DIA-NN report) contaminant protein groups whose every peptide in the
                    DE's input report (after its q-value filter) is also one of a target
                    protein's digest peptides (the build's rule) -- shared with it, or the same
                    but for I/L, which no spectrum tells apart (horse myoglobin on human heart):
                    not attributable to the contaminant's species -- the sample's own protein or
                    handling contamination (WARN); INFO when it could not be checked
  de_signal         0 significant (WARN: underpowered) or >50% significant
                    (WARN: likely batch/normalization/confounding artefact)
  first_vs_final_pass
                    5-step DIA-NN chain only: runs whose final pass kept less than half
                    their first-pass precursors (pass_comparison.py, recorded by step 5 in
                    search_provenance.json; counts before protein-group FDR) -- WARN,
                    offering the first-pass report only for runs where it has more rows
                    after the q-value filter; INFO when the DE already used it, or when no
                    comparison was recorded
  search_measurement
                    (DIA-NN) the scan window and/or mass accuracy the search was to measure on
                    these runs were NOT measured: the measurement failed and the search fell back
                    (search_provenance.json scan_window.mode fallback_auto / mass_acc.mode
                    fallback_default) -- a CAUTION, worded by probe_fallback.caution() (WARN);
                    and a mass accuracy pinned at the SOP floor over probed runs that disagreed
                    (the probe's warnings, probe_fallback.measurement_note) (WARN)
  pipeline checks   the DE pipeline's own checks, as de_provenance.json `checks` records them
                    (run_silac.R: silac_channel_balance, silac_replicates,
                    silac_label_incorporation, silac_label_declaration, ...), quoted
  quantification    the DE pipeline's own caveat about its protein quantities
                    (de_provenance.json `caveat`, e.g. a Sage DE's highest-peptide rollup),
                    quoted as the pipeline states it (WARN)
  sage_lfq          (--search-out, Sage searches) the runs' MS1 mass error does not fit
                    Sage's LFQ window, or Sage kept (almost) no MS1 peaks: quantities
                    unreliable (WARN) -- sage_lfq_check.py's record, quoted as written
  normalization_check
                    the DE's quantities (normalised / non-normalised), the experiment type's
                    default and whether normalization_check.py's data check agreed -- quoted
                    from de_provenance.json: WARN when it did not and someone chose, when no
                    type was recorded, and when no check was run (the delivery is then held)
  normalisation_stability
                    (a DIA-NN report) DIA-NN's per-run Normalisation.Instability, from the
                    <report>.stats.tsv beside the DE input, for the runs in conditions.csv:
                    WARN naming every run above check_report_runs.NORM_INSTABILITY_WARN, with
                    the alternatives DIA-NN documents (--global-norm / --no-norm, and the
                    report's non-normalised Precursor.Quantity); every value in AUDIT.json

Usage:
  python3 audit_results.py --out AUDIT.md \
      --conditions input/conditions.csv --de-dir output/tables \
      [--acquisition-json acq.json] [--adjp 0.05] [--logfc 1] \
      [--min-proteins 500] [--max-missing 0.5] \
      [--fasta-meta search.fasta.meta.json]   # default: ./search.fasta.meta.json if present
      [--search-prov output/search/search_provenance.json]
                                              # default: beside de_provenance.json's input
      [--search-out search_out]               # the search's output dir (Sage LFQ check)
"""
import sys, os, re, csv, json, glob, argparse
from collections import Counter, defaultdict

# The list's ONE definition is the FASTA sidecar; the wording for lost proteins lives there too.
from fetch_fasta import (target_contaminants, seen_only_as_cont, lost_to_contaminants_message,
                         keratin_sample_recorded, CONT_TAG, CONTAMINANT_TAGS,
                         is_contaminant_accession, contaminant_groups_on_target_peptides)
# What a subject / animal column is called: one definition, shared with --map.
from collect_conditions import subject_header, subject_assessment

CONTAMINANT_PATTERNS = ("KRT", "KRTAP",            # keratins (skin/hair)
                        "TRYP", "PRSS1", "TRY1",   # trypsin (digestion)
                        "CASEIN", "CSN1", "CSN2",  # casein (milk)
                        )
CONTAMINANT_WORDS = ("keratin", "trypsin", "casein", "contaminant", "cRAP")


def add(findings, cid, status, message, detail=None, plain=None, label=None):
    """One finding. `message` is the technical record (AUDIT.md, staff); `plain` is what it means
    for the results in words a client reads, when it differs from its check's usual words
    (PLAIN); `label` names a check PLAIN does not know (a pipeline's own)."""
    f = {"check": cid, "status": status, "message": message, "detail": detail or {}}
    if plain:
        f["plain"] = plain
    if label:
        f["label"] = label
    findings.append(f)


# ---- the client's view: one plain row per check ----------------------------------------------
# The report's "Audit & caveats" is a table a client can read: what was checked, the result, what
# it means for the results (client_section_md). A staff report: the raw AUDIT.md -- check codes,
# the normalisation check's B1-B4 strings, Normalisation.Instability values, "implausibly high
# ... not real biology" -- went to clients, and neither staff nor clients could read it. The words
# are written HERE, beside each check, never by the model writing the report; a finding that
# means something else than its check's usual words carries its own (`plain`), and a pipeline's
# own checks (de_provenance.json `checks`, e.g. run_silac.R's) carry their own label and words.
# The technical AUDIT.md stays the full record.
RESULT_WORDS = {"PASS": "✅ Passed", "WARN": "⚠️ Caution", "FAIL": "⛔ Problem", "INFO": "ℹ️ Note"}
_RANK = {"INFO": 0, "PASS": 1, "WARN": 2, "FAIL": 3}
# A check with no words of its own (a pipeline's check from before 2.13.2): never its message.
GENERIC_MEANING = {
    "PASS": "No problem found.",
    "WARN": "A caution about these results; the Core's technical record (AUDIT.md) says what it "
            "is, and the Core can explain it.",
    "FAIL": "A problem to resolve before relying on these results; the Core's technical record "
            "(AUDIT.md) says what it is.",
    "INFO": "For information; the Core's technical record (AUDIT.md) has the details.",
}
NOT_CHECKED = "Not checked: the information this check needs was not available here."
PLAIN_DE_SIGNAL_NONE = ("No protein changed significantly in at least one comparison: the study "
                        "may be too small to detect the changes, or the changes are small.")
PLAIN_DE_SIGNAL_HIGH = ("More than half of the proteins changed significantly in at least one "
                        "comparison. That is unusual: it can be real when the groups are very "
                        "different materials, or it can come from a technical difference between "
                        "the groups (loading, batch or normalisation), so read the results with "
                        "that in mind.")
PLAIN_DE_SIGNAL_GENERIC = ("The number of proteins called significant is unusual in at least one "
                           "comparison (none, or more than half of them), so read the results "
                           "with that in mind.")
# A check that emits one status for different reasons: each variant's words, used where it is
# emitted (plain=) and, for a record from before 2.13.2 that has none, picked by OLD_RULES below.
P_NO_GROUPS = "No groups were assigned to the samples, so nothing could be compared."
P_FEW_REPLICATES = ("At least one group has fewer than two replicates, so its variability cannot be "
                    "measured and its comparisons are not valid.")
P_EMPTY_MATRIX = "No protein quantities were found."
P_FEW_PROTEINS = ("Fewer proteins were quantified than is typical, which can point to a sample, "
                  "search or database problem.")
P_CONT_FEW = ("Few common lab contaminants (skin keratins, the digestion enzyme, milk proteins) were "
              "found.")
P_CONT_MANY = ("Common lab contaminants (skin keratins, the digestion enzyme or milk proteins) are "
               "present at a noticeable level, which can affect the quantification.")
P_KERATIN_FEW = ("The sample itself is keratin, so keratins are not counted as contaminants; few "
                 "other common lab contaminants (the digestion enzyme, milk proteins) were found.")
P_KERATIN_MANY = ("The sample itself is keratin, so keratins are not counted as contaminants; other "
                  "common lab contaminants (the digestion enzyme or milk proteins) are present at a "
                  "noticeable level, which can affect the quantification.")
P_KERATIN_NONE = ("The sample itself is keratin, so keratins are not counted as contaminants; no "
                  "other common contaminant was found.")
P_PASSES_LOST = ("Some samples kept less than half of their first-pass identifications in the "
                 "search's final pass, so their results rest on fewer peptides.")
P_PASSES_FIRST_MORE = ("For some samples the search's first pass kept more identifications than its "
                       "final pass, which these results use; the first pass could be used instead.")
P_PASSES_USED_FIRST = ("The analysis used the search's first-pass results, so this does not affect "
                       "them.")
P_PASSES_UNRECORDED = "Not checked: the search did not record this comparison."
P_MEASURE_FALLBACK = ("Some search settings that are normally measured from these runs could not be "
                      "measured, so standard values were used instead; the results may be somewhat "
                      "less complete.")
P_MEASURE_SPREAD = ("The runs disagreed about the instrument's mass accuracy, and the search used the "
                    "standard setting for it.")
P_NORM_NA = "Not applicable to this kind of analysis."
P_NORM_AGREED = ("The usual normalisation for this kind of experiment was used, and a check of the "
                 "data agreed with it.")
P_NORM_UNDECIDED = "The choice of normalisation was not yet checked against the data."
P_NORM_NO_TYPE = ("No experiment type was recorded, so there was no usual normalisation to compare "
                  "with.")
P_NORM_TRIPPED = ("The usual normalisation for this kind of experiment did not fit these samples, so "
                  "the Core chose which quantities to use.")
P_STAB_RAW = "Not relevant here: the analysis did not use DIA-NN's normalised amounts."
P_STAB_UNCHECKED = "Not checked: the search's statistics this check reads were not available."

PLAIN = {
    "replication": ("Replicates per group", {
        "PASS": "Every group has at least three replicates, enough for the statistics.",
        "WARN": "At least one group has only two replicates: its comparisons still work, but can "
                "miss real changes.",
        "FAIL": "At least one group has fewer than two replicates, so its variability cannot be "
                "measured and its comparisons are not valid."}),
    "group_balance": ("Group sizes", {
        "PASS": "The groups are of similar size.",
        "WARN": "The groups differ a lot in size, which lowers the power of the comparisons."}),
    "confounding": ("Batch and groups", {
        "PASS": "No batch or other recorded factor lines up exactly with the groups.",
        "FAIL": "A batch or other recorded factor lines up exactly with the groups, so its effect "
                "cannot be told apart from a real difference between the groups."}),
    "numeric_covariate": ("Measured factors in the statistics", {
        "PASS": "A measured factor was included in the statistics; the groups explain less of it "
                "than the level at which the Core cautions, so it is not expected to distort the "
                "comparisons.",
        "WARN": "A measured factor included in the statistics largely follows the groups, so "
                "adjusting for it can remove part of the real difference; the comparisons should "
                "also be checked without it."}),
    "acquisition_mix": ("Acquisition method", {
        "PASS": "All runs used the same acquisition method.",
        "FAIL": "The runs mix two acquisition methods (DIA and DDA), which cannot be quantified "
                "together, so the comparison is not valid as run."}),
    "instrument_mix": ("Instruments", {
        "PASS": "All runs come from one instrument model.",
        "WARN": "The runs come from more than one instrument model, which can add a technical "
                "difference between samples."}),
    "id_depth": ("Proteins quantified", {
        "PASS": "The number of proteins quantified is not unusually low.",
        "WARN": "Fewer proteins were quantified than is typical, which can point to a sample, "
                "search or database problem."}),
    "missingness": ("Missing values", {
        "PASS": "Missing values are within the usual range.",
        "WARN": "A large share of values is missing, which weakens the quantification."}),
    "contamination": ("Common lab contaminants", {
        "PASS": "Few common lab contaminants (skin keratins, the digestion enzyme, milk proteins) "
                "were found.",
        "WARN": "Common lab contaminants (skin keratins, the digestion enzyme or milk proteins) are "
                "present at a noticeable level, which can affect the quantification."}),
    "target_contaminants": ("Sample proteins identical to contaminants", {
        "PASS": "None of the sample's own proteins that share a sequence with a common contaminant "
                "were quantified.",
        "WARN": "Some of the sample's own proteins have the same sequence as a common lab "
                "contaminant. They are reported as the sample's proteins, but part of their "
                "signal could come from contamination.",
        "INFO": NOT_CHECKED}),
    "contaminant_overlap": ("Proteins lost to the contaminant list", {
        "WARN": "Some real proteins of the sample are missing from the results because the search "
                "database listed them only as contaminants; a corrected database and a new search "
                "would recover them."}),
    "contaminant_evidence": ("Contaminants explained by the sample's proteins", {
        "PASS": "Each contaminant found has evidence of its own.",
        "WARN": "Some proteins listed as contaminants rest only on peptides the sample's own "
                "proteins also have: they are the sample's own proteins or handling contamination "
                "(such as skin keratin), not evidence of the contaminant species.",
        "INFO": NOT_CHECKED}),
    "de_signal": ("Number of significant proteins", {
        "PASS": "The number of proteins called significant is in a plausible range in every "
                "comparison.",
        "WARN": PLAIN_DE_SIGNAL_GENERIC}),
    "first_vs_final_pass": ("Search consistency between passes", {
        "PASS": "No sample lost a large share of its identifications between the search's two "
                "passes.",
        "WARN": "Some samples kept less than half of their first-pass identifications in the "
                "search's final pass, so their results rest on fewer peptides.",
        "INFO": "Not checked: the search did not record this comparison."}),
    "search_measurement": ("Search settings measured from the data", {
        "WARN": "Some search settings that are normally measured from these runs could not be "
                "measured, so standard values were used instead; the results may be somewhat less "
                "complete."}),
    "quantification": ("How protein amounts were calculated", {
        "WARN": "The method used to turn peptide signals into protein amounts has a known "
                "limitation for this kind of data, so treat small differences with care."}),
    "sage_lfq": ("Quantification settings (Sage)", {
        "PASS": "The quantification settings fit these runs.",
        "WARN": "The quantification settings did not fit these runs' measured mass error, so the "
                "protein amounts are less reliable; the identifications are not affected.",
        "INFO": "Not checked: the search's own check could not run."}),
    "normalization_check": ("Normalisation choice", {
        "PASS": "The usual normalisation for this kind of experiment was used, and a check of the "
                "data agreed with it.",
        "WARN": "The choice of normalisation was not confirmed against the data."}),
    "normalisation_stability": ("DIA-NN normalisation stability", {
        "PASS": "DIA-NN's normalisation between runs was stable in every run.",
        "WARN": "DIA-NN's normalisation between runs was unstable for some runs, and the analysis "
                "used its normalised amounts, so those runs' values may be distorted.",
        "INFO": "Not checked: the search's statistics this check reads were not available."}),
    "search_checks": ("The search's own checks", {"INFO": NOT_CHECKED}),
}


# (check, status) emitted for more than one reason: PLAIN's words are not used for them; every
# site gives its own (tests/test_audit_plain.py holds the sites to this).
AMBIGUOUS = {("replication", "FAIL"), ("id_depth", "WARN"), ("contamination", "PASS"),
             ("contamination", "WARN"), ("de_signal", "WARN"), ("first_vs_final_pass", "INFO"),
             ("first_vs_final_pass", "WARN"), ("search_measurement", "WARN"),
             ("normalization_check", "PASS"), ("normalization_check", "WARN"),
             ("normalisation_stability", "INFO"), ("contaminant_evidence", "INFO"),
             ("target_contaminants", "INFO")}
_STATEMENT = re.compile(r"(Core staff chose [^.]*\.)")


def _with_statement(words, msg):
    m = _STATEMENT.search(msg)
    return words + (" " + m.group(1) if m else "")


# A record from before 2.13.2 has no `plain`: its variant is read from the message, as that
# version wrote it. No rule matches -> the generic caution, never words that may be false.
OLD_RULES = (
    ("replication", "FAIL", r"^No groups assigned", P_NO_GROUPS),
    ("replication", "FAIL", r"<2 ", P_FEW_REPLICATES),
    ("id_depth", "WARN", r"^Expression matrix is empty", P_EMPTY_MATRIX),
    ("id_depth", "WARN", r"^Only \d+ proteins", P_FEW_PROTEINS),
    ("contamination", "PASS", r"^Keratin-matrix sample", P_KERATIN_NONE),
    ("contamination", "PASS", r"keratin-matrix sample", P_KERATIN_FEW),
    ("contamination", "WARN", r"keratin-matrix sample", P_KERATIN_MANY),
    ("contamination", "PASS", r"likely contaminant protein", P_CONT_FEW),
    ("contamination", "WARN", r"likely contaminant protein", P_CONT_MANY),
    ("de_signal", "WARN", r": 0 proteins significant", PLAIN_DE_SIGNAL_NONE),
    ("de_signal", "WARN", r"proteins significant —", PLAIN_DE_SIGNAL_HIGH),
    ("first_vs_final_pass", "INFO", r"used the first-pass report", P_PASSES_USED_FIRST),
    ("first_vs_final_pass", "INFO", r"comparison is not in", P_PASSES_UNRECORDED),
    ("first_vs_final_pass", "WARN", r"^\d+ of \d+ runs kept less than half", P_PASSES_LOST),
    ("first_vs_final_pass", "WARN", r"^no run lost more than half", P_PASSES_FIRST_MORE),
    ("search_measurement", "WARN", r"^NOTE: the mass accuracy was pinned|SOP floor", P_MEASURE_SPREAD),
    ("search_measurement", "WARN", r"^CAUTION", P_MEASURE_FALLBACK),
    ("normalization_check", "PASS", r"^Not applicable", P_NORM_NA),
    ("normalization_check", "PASS", r"the data check agreed", P_NORM_AGREED),
    ("normalization_check", "WARN", r"^This DE is one of|^The normalisation check was not run",
     P_NORM_UNDECIDED),
    ("normalization_check", "WARN", r"^No experiment type was recorded", P_NORM_NO_TYPE),
    ("normalization_check", "WARN", r"did not agree with its default", P_NORM_TRIPPED),
    ("normalisation_stability", "INFO", r"^Not relevant", P_STAB_RAW),
    ("normalisation_stability", "INFO", r"^Not assessed", P_STAB_UNCHECKED),
    ("contaminant_evidence", "INFO", r"^Not assessed", NOT_CHECKED),
    ("target_contaminants", "INFO", r"^Not assessed", NOT_CHECKED),
)


def _old_words(check, status, msg, detail):
    if check == "normalization_check" and status == "PASS" and \
            isinstance(detail, dict) and detail.get("status") == "not_applicable":
        return P_NORM_NA
    for c, st, rx, words in OLD_RULES:
        if c == check and st == status and re.search(rx, msg):
            return (_with_statement(words, msg) if words in (P_NORM_NO_TYPE, P_NORM_TRIPPED)
                    else words)
    return None


def client_row(f):
    """(plain label, what it means for the results) for one finding: its own words (`plain`,
    written where it was found); for a record without them, the variant its message names
    (OLD_RULES), else the check's words when the check has one meaning per status (PLAIN), else
    the generic caution -- never words that may not fit."""
    check = str(f.get("check") or "")
    label, words = PLAIN.get(check, (check.replace("_", " ").strip().capitalize(), {}))
    status = f.get("status")
    meaning = f.get("plain") or _old_words(check, status, str(f.get("message") or ""),
                                           f.get("detail"))
    if not meaning:
        meaning = (None if (check, status) in AMBIGUOUS else words.get(status)) \
            or GENERIC_MEANING.get(status, "")
    return str(f.get("label") or label), str(meaning)


def client_section_md(findings):
    """The report's "Audit & caveats": one row per check, the worst result of its findings and
    their meanings at that result. Markdown; the technical record stays AUDIT.md."""
    order, by = [], {}
    for f in findings or []:
        if isinstance(f, dict) and f.get("check") and f.get("status") in RESULT_WORDS:
            if f["check"] not in by:
                order.append(f["check"])
            by.setdefault(f["check"], []).append(f)
    if not order:
        return ""
    cell = lambda t: " ".join(str(t).replace("|", "/").split())      # noqa: E731
    L = ["Every analysis goes through a set of automatic checks. Each row says what was checked, "
         "the result, and what it means for these results; the full technical record is "
         "AUDIT.md.", "", "| Check | Result | What it means for the results |", "|---|---|---|"]
    for check in order:
        fs = by[check]
        status = max((f["status"] for f in fs), key=_RANK.get)
        rows = [client_row(f) for f in fs if f["status"] == status]
        meanings = list(dict.fromkeys(m for _, m in rows if m))
        L.append(f"| {cell(rows[0][0])} | {RESULT_WORDS[status]} | {cell(' '.join(meanings))} |")
    return "\n".join(L) + "\n"


_MD_FINDING = re.compile(r"^- (\S+)\s+\*\*([A-Za-z0-9_.-]+)\*\*\s+—\s+(.*)$")
_MD_ICON = {"✅": "PASS", "⚠": "WARN", "⛔": "FAIL", "❌": "FAIL", "ℹ": "INFO"}


def findings_from_md(text):
    """An AUDIT.md's findings (a session from before AUDIT.json, or without it): the lines main()
    writes, "- <icon> **<check>** — <message>". [] when there are none."""
    out = []
    for ln in (text or "").splitlines():
        m = _MD_FINDING.match(ln.strip())
        if m and m.group(1)[:1] in _MD_ICON:
            out.append({"check": m.group(2), "status": _MD_ICON[m.group(1)[:1]],
                        "message": m.group(3)})
    return out


def read_csv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def replication_unit(rows, block_col=None):
    """The column whose DISTINCT values are the replicates: the block run_de.R recorded
    (de_provenance block_column), else a conditions.csv column named for the animal /
    subject (collect_conditions.subject_header). None = every run is its own replicate.
    Three injections of one mouse are one mouse: counting runs called them three."""
    if not rows:
        return None
    if block_col and block_col in rows[0]:
        return block_col
    groups = {r["File.Name"]: (r.get("Group") or "").strip() for r in rows if r.get("File.Name")}
    for c in rows[0]:
        if c in ("File.Name", "Group") or not subject_header(c):
            continue
        vals = {r["File.Name"]: (r.get(c) or "").strip() for r in rows if r.get("File.Name")}
        # the header is not enough (a 'Subject' of M/F): the values must look like subjects
        if all(vals.values()) and not subject_assessment(vals, groups):
            return c
    return None


def audit_conditions(findings, rows, adjp, block_col=None):
    runs = Counter(r.get("Group", "").strip() for r in rows if r.get("Group", "").strip())
    if not runs:
        add(findings, "replication", "FAIL", "No groups assigned in conditions.csv.",
            plain=P_NO_GROUPS)
        return runs
    unit = replication_unit(rows, block_col)
    if unit:
        members = defaultdict(set)
        for r in rows:
            g = r.get("Group", "").strip()
            if g:
                members[g].add(r[unit].strip())
        groups = Counter({g: len(v) for g, v in members.items()})
        what = f"distinct {unit} (n = {unit}, not runs)"
    else:
        groups, what = runs, "replicates"
    detail = {"group_sizes": dict(groups), "unit": unit or "run", "runs_per_group": dict(runs)}
    singletons = [g for g, n in groups.items() if n < 2]
    small = [g for g, n in groups.items() if n == 2]
    if singletons:
        add(findings, "replication", "FAIL",
            f"Group(s) with <2 {what} have no within-group variance — differential statistics are not valid: {singletons}.",
            detail, plain=P_FEW_REPLICATES)
    elif small:
        add(findings, "replication", "WARN",
            f"Group(s) with only 2 {what}: {small}. Usable but low power; 3+ is recommended.",
            detail)
    else:
        add(findings, "replication", "PASS", f"All groups have ≥3 {what}.", detail)
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


def audit_numeric_covariates(findings, prov):
    """A covariate fitted as a number (run_de.R --numeric-covariates) that the groups nearly
    determine: every comparison is adjusted for it, so they rest on its fitted slope, and if the
    groups themselves change it, adjusting removes part of the effect being tested. Read from
    de_provenance.json `covariates` -- `confounded` per term and run_de.R's caution in `warnings`,
    quoted -- never re-judged here. Nothing to say without a numeric covariate (a factor covariate
    is the `confounding` check's)."""
    rec = prov.get("covariates") if isinstance(prov.get("covariates"), dict) else {}
    terms = [t for t in (rec.get("terms") or []) if isinstance(t, dict) and t.get("type") == "numeric"]
    if not terms:
        return
    thr = rec.get("confound_r2_threshold")
    cautions = [str(w) for w in (rec.get("warnings") or []) if w]

    def se_of(t):
        se = t.get("se_inflation")
        return {k: v for k, v in se.items() if isinstance(v, (int, float))} if isinstance(se, dict) else {}

    def r2(t):
        v = t.get("group_r2")
        return f"{v * 100:.0f}%" if isinstance(v, (int, float)) else "an unrecorded share"
    for t in terms:
        col, se = t.get("column") or "?", se_of(t)
        detail = {"column": col, "group_r2": t.get("group_r2"), "threshold": thr,
                  "se_inflation": t.get("se_inflation")}
        if t.get("confounded") is True:
            # run_de.R's own caution, as it wrote it (naming each comparison whose groups the
            # covariate nearly matches); composed from the term only when the record has none
            pc = t.get("per_contrast") if isinstance(t.get("per_contrast"), dict) else {}
            conf = [c for c in (t.get("confounded_contrasts") or []) if isinstance(pc.get(c), dict)]
            detail["confounded_contrasts"] = conf
            caution = next((w for w in cautions if w.startswith(f"{col} ")), None) or (
                (f"{col} nearly matches the groups compared in "
                 + ", and in ".join(f"{c} ({'; '.join(str(w) for w in pc[c].get('why') or [])})"
                                    for c in conf)
                 + f". {'That comparison is' if len(conf) == 1 else 'Those comparisons are'} "
                   f"adjusted for it, so {'it rests' if len(conf) == 1 else 'they rest'} largely on "
                   f"its fitted slope" if conf else
                 f"{col} nearly matches the groups: Group explains {r2(t)} of its variance. The "
                 f"group comparisons are adjusted for it, so they rest largely on its fitted slope"
                 + (f", and the design alone widens their standard errors up to "
                    f"{max(se.values()):.1f}x ({', '.join(f'{k} {v:.2f}x' for k, v in se.items())})"
                    if se else ""))
                + f". If {col} differs between the groups because of the grouping itself (a "
                  f"consequence of the treatment), adjusting for it removes part of the effect "
                  f"being tested: ask whether it belongs in the model.")
            add(findings, "numeric_covariate", "WARN",
                f"{caution} Run the DE without {col} as well, and compare, before trusting the "
                f"adjusted comparisons.", detail)
        else:
            add(findings, "numeric_covariate", "PASS",
                f"{col} was fitted as a number; Group explains {r2(t)} of its variance"
                + (f" (cautioned from {thr * 100:.0f}%)" if isinstance(thr, (int, float)) else "")
                + (f", and the design alone widens the comparisons' standard errors by at most "
                   f"{max(se.values()):.2f}x" if se else "") + ".", detail)


def audit_quantification(findings, prov):
    """The pipeline's own caveat about its quantities (build_maxlfq.R maxlfq_descriptor ->
    de_provenance.json `caveat`), quoted -- this check knows no pipeline by name."""
    if prov.get("caveat"):
        add(findings, "quantification", "WARN",
            f"{prov.get('display_label') or 'This pipeline'}: {prov['caveat']}",
            {"rollup_method": prov.get("rollup_method")})


def audit_search(findings, search_out):
    """What the search itself recorded about its own output. Sage: sage_lfq_check.py's record
    (the ONE judgement of whether the LFQ window fit the runs' MS1 mass error), quoted, not
    re-derived."""
    import sage_lfq_check
    hit = sage_lfq_check.audit_finding(sage_lfq_check.load(search_out))
    if hit:
        add(findings, "sage_lfq", *hit)


def audit_matrix(findings, em_path, min_proteins, max_missing, keratin_sample=False):
    rows = read_csv(em_path)
    if not rows:
        add(findings, "id_depth", "WARN", "Expression matrix is empty.", plain=P_EMPTY_MATRIX)
        return
    idcols = [c for c in ("Protein.Group", "Genes", "Protein.Names") if c in rows[0]]
    sample_cols = [c for c in rows[0].keys() if c not in idcols]
    n_prot = len(rows)
    if n_prot < min_proteins:
        add(findings, "id_depth", "WARN",
            f"Only {n_prot} proteins quantified — lower than typical (<{min_proteins}). "
            "Check the search parameters and that the FASTA matches the sample's organism.",
            plain=P_FEW_PROTEINS)
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
    # contamination. Keratin is a contaminant for MOST matrices, but is the ANALYTE
    # for keratin samples (nail / hair / wool / skin / feather) — never flag it there.
    # (see references/anomaly-checks.md "Keratin-matrix samples").
    pats = tuple(p for p in CONTAMINANT_PATTERNS if not (keratin_sample and p in ("KRT", "KRTAP")))
    words = tuple(w for w in CONTAMINANT_WORDS if not (keratin_sample and w == "keratin"))
    glab = "Genes" if "Genes" in rows[0] else ("Protein.Names" if "Protein.Names" in rows[0] else None)
    if glab:
        names = [(r.get(glab) or "").upper() for r in rows]
        hits = [n for n in names if any(n.startswith(p) or p in n for p in pats)
                or any(w.upper() in n for w in words)]
        label = "trypsin/casein" if keratin_sample else "keratins/trypsin/casein"
        if hits:
            frac_c = len(hits) / n_prot
            status = "WARN" if frac_c >= 0.02 or len(hits) >= 10 else "PASS"
            add(findings, "contamination", status,
                f"{len(hits)} likely contaminant protein(s) detected ({label}). "
                + ("A high contaminant load can distort normalization and quant — consider filtering them."
                   if status == "WARN" else "Low level; usually fine.")
                + (" [keratin-matrix sample: keratins treated as the analyte, not flagged]" if keratin_sample else ""),
                {"examples": sorted(set(hits))[:10]},
                plain={(True, "WARN"): P_KERATIN_MANY, (True, "PASS"): P_KERATIN_FEW,
                       (False, "WARN"): P_CONT_MANY, (False, "PASS"): P_CONT_FEW}[
                           (bool(keratin_sample), status)])
        elif keratin_sample:
            add(findings, "contamination", "PASS",
                "Keratin-matrix sample — keratins are the analyte (not flagged); no trypsin/casein issue.",
                plain=P_KERATIN_NONE)


def _tokens(cell):
    return {t.strip().upper() for t in (cell or "").split(";") if t.strip()}


def removed_contaminant_groups(de_dir):
    """-> (record, rows): run_de.R's contaminant record and the protein groups its filter
    removed. The filter takes Cont_ groups out of Expression_Matrix.csv, so a real protein
    that sat in the database only as a Cont_ entry is visible only in this table now."""
    try:
        with open(os.path.join(de_dir or "", "de_provenance.json")) as fh:
            rec = json.load(fh).get("contaminants") or {}
    except (OSError, ValueError):
        return {}, []
    table = rec.get("removed_table") if rec.get("removed") else None
    path = os.path.join(de_dir, table) if table else None
    return rec, (read_csv(path) if path and os.path.exists(path) else [])


def audit_target_contaminants(findings, meta_path, em_path, de_dir, adjp, keratin_sample=False):
    """Proteins that are both a <organism> protein and a common contaminant.

    Why (HIVE, 2026-09-24): contaminant entries identical to target proteins made DIA-NN
    report ACTB, EEF1A1 and KRT8 only as Cont_ groups, excluded from quant. fetch_fasta.py
    now removes those entries, so the proteins are quantified under their own accession --
    which also means skin/hair keratins (usually handling contamination) count toward
    normalisation. That trade-off is surfaced here, never hidden and never "fixed" by
    excluding them: the list comes from the FASTA sidecar, not from a copy kept here.
    """
    if not meta_path:
        add(findings, "target_contaminants", "INFO",
            "Not assessed: no <fasta>.meta.json from fetch_fasta.py (pass --fasta-meta), so "
            "proteins that are also common-contaminant sequences cannot be identified.", plain=NOT_CHECKED)
        return
    try:
        with open(meta_path) as fh:
            # Keratin-matrix sample: keratin is the analyte, so never flag it here.
            tc = target_contaminants(json.load(fh), keratin_sample)
    except (OSError, ValueError) as e:
        add(findings, "target_contaminants", "INFO",
            f"Not assessed: could not read {meta_path} ({e}).", plain=NOT_CHECKED)
        return
    org = tc["organism"] or "target-organism"

    em_rows = read_csv(em_path) if em_path and os.path.exists(em_path) else []

    # Real proteins that exist only as Cont_ entries: a database used as-is, one built with
    # --keep-target-contaminants, one built BEFORE the overlap check (legacy sidecar -- the
    # databases that lost ACTB/EEF1A1/KRT8), or one built by the identity rule alone (near-
    # identical entries: bovine EEF1A1/YWHAZ vs mouse). The matrix makes it concrete where it can.
    kept = tc["kept_as_contaminant"]
    seen = seen_only_as_cont(kept, [(_tokens(r.get("Protein.Group")), r.get("Genes") or "?")
                                    for r in em_rows])
    msg = lost_to_contaminants_message(tc, seen)
    if msg:
        rec, removed = removed_contaminant_groups(de_dir)
        gone = seen_only_as_cont(kept, [(_tokens(r.get("Protein.Group")), r.get("Genes") or "?")
                                        for r in removed])
        if rec.get("removed"):
            msg += (f" run_de.R's contaminant filter then removed every {rec.get('tag') or CONT_TAG}"
                    f" group from the DE ({rec.get('n_protein_groups')} groups, listed in "
                    f"{rec.get('removed_table')})"
                    + (f", among them {', '.join(gone[:12])}" if gone else "")
                    + "; --keep-contaminants keeps them, at the price of testing the true "
                      "contaminants too.")
        add(findings, "contaminant_overlap", "WARN", msg,
            {"fasta_meta": meta_path, "legacy_database": bool(tc.get("legacy_note")),
             "database_state": tc.get("state"),
             "genes": sorted({r.get("gene") or r.get("target_acc") or "?" for r in kept}),
             "seen_only_as_cont": seen, "removed_by_de_filter": gone})

    if not tc["dropped"]:
        return
    genes, accs = tc["genes"], tc["accessions"]

    def match(r):
        # Gene OR accession: an NCBI database has no GN= in its headers, so its records
        # carry no gene and only the accession can match.
        return bool(_tokens(r.get("Genes")) & genes) or bool(_tokens(r.get("Protein.Group")) & accs)

    def label(r):
        return (r.get("Genes") or r.get("Protein.Group") or "?").split(";")[0]

    hits = sorted({label(r) for r in em_rows if match(r)})
    sig = {}
    for f in sorted(glob.glob(os.path.join(de_dir or "", "DE_*.csv"))):
        s = set()
        for r in read_csv(f):
            try:
                if float(r.get("adj.P.Val", "nan")) < adjp and match(r):
                    s.add(label(r))
            except ValueError:
                continue
        if s:
            sig[os.path.basename(f)] = sorted(s)
    if not hits and not sig:
        add(findings, "target_contaminants", "PASS",
            f"None of the {len(tc['dropped'])} {org} proteins that are also common-contaminant "
            f"sequences were quantified.", {"fasta_meta": meta_path})
        return
    msg = (f"{len(hits)} quantified {org} protein group(s) are also common-contaminant sequences "
           f"({', '.join(hits[:12])}{', ...' if len(hits) > 12 else ''}): possible contamination, "
           f"KEPT in quantification and normalisation. Skin/hair keratins (KRT1/2/9/10, KRTAPs, "
           f"FLG, HRNR) are usually handling contamination; ACTB/EEF1A1/tubulins -- and the "
           f"simple-epithelial keratins KRT8/18/19/7 in epithelial cells such as HeLa -- are "
           f"usually endogenous.")
    if sig:
        msg += (" Significant in DE -- check these are not contamination before interpreting "
                "them: " + "; ".join(f"{k}: {', '.join(v[:10])}" for k, v in sig.items()) + ".")
    add(findings, "target_contaminants", "WARN", msg,
        {"fasta_meta": meta_path, "quantified": hits, "significant": sig})


REPORT_PEPTIDE_COLUMNS = ("Run", "Protein.Group", "Protein.Ids", "Stripped.Sequence")


def _split_ids(cell):
    return {t.strip() for t in str(cell or "").split(";") if t.strip()}


def report_contaminant_peptides(report, q_cutoff, q_columns=None):
    """({contaminant protein group: {peptide: accessions (Protein.Ids)}}, {group: runs}, q-columns
    applied, None) from a DIA-NN report (parquet or tsv), its rows kept as run_de.R keeps them: the
    q-columns it recorded (else diann_q_columns.fdr_columns of the report), each at cutoff_for(),
    a row dropped only when a value exceeds its cutoff -- or (None, None, [], why not). Only rows
    naming a contaminant accession are read (a pyarrow filter), so a 400-run report is not loaded
    whole. A report that cannot be read is a why-not, never a traceback that leaves no AUDIT.md
    (error-handling review, 2.12)."""
    if not report or not os.path.isfile(report):
        return None, None, [], f"the DE's input report ({report}) is not readable here"
    if report.endswith(".parquet"):
        try:
            import pyarrow  # noqa: F401
        except ImportError:
            return None, None, [], "pyarrow is not installed in this Python, so the report could not be read"
    try:
        return _contaminant_peptides(report, q_cutoff, q_columns)
    except (OSError, ValueError, csv.Error) as e:          # pyarrow's ArrowInvalid is a ValueError
        name = os.path.basename(report)
        why = " ".join(str(e).replace(os.path.abspath(report), name).replace(report, name).split())
        return None, None, [], (f"{name} could not be read ({type(e).__name__}: "
                                f"{why[:200]}{'...' if len(why) > 200 else ''})")


def _contaminant_peptides(report, q_cutoff, q_columns):
    from diann_q_columns import fdr_columns, cutoff_for
    name = os.path.basename(report)
    groups, runs = defaultdict(dict), defaultdict(set)

    def take(run, pg, ids, pep):
        if pep and pg and any(is_contaminant_accession(t) for t in _split_ids(pg)):
            groups[pg].setdefault(pep, set()).update(_split_ids(ids))
            runs[pg].add(run)

    if report.endswith(".parquet"):
        import pyarrow as pa
        import pyarrow.compute as pc
        import pyarrow.dataset as pads
        d = pads.dataset(report, format="parquet")
        names = d.schema.names
    else:
        with open(report, newline="") as fh:
            names = next(csv.reader(fh, delimiter="\t"), [])
    missing = [c for c in REPORT_PEPTIDE_COLUMNS if c not in names]
    if missing:
        return None, None, [], (f"{name} has no {', '.join(missing)} column (not a DIA-NN "
                                f"precursor report)")
    qcols = [c for c in (q_columns or fdr_columns(names)) if c in names]
    cuts = {c: cutoff_for(c, q_cutoff) for c in qcols}
    if report.endswith(".parquet"):
        tags = "|".join(re.escape(t) for t in CONTAMINANT_TAGS)
        # cast: a dictionary-encoded column has no regex kernel (DIA-NN writes plain strings)
        keep = pc.match_substring_regex(pads.field("Protein.Ids").cast(pa.string()),
                                        f"(^|;)({tags})")
        for c, cut in cuts.items():
            keep = keep & ((pads.field(c) <= cut) | pads.field(c).is_null())
        for b in d.to_batches(columns=list(REPORT_PEPTIDE_COLUMNS), filter=keep):
            for row in zip(*(b.column(c).to_pylist() for c in REPORT_PEPTIDE_COLUMNS)):
                take(*row)
    else:
        with open(report, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if not any(is_contaminant_accession(t) for t in _split_ids(r.get("Protein.Ids"))):
                    continue
                try:
                    if any((r.get(c) or "").strip() not in ("", "NA", "NaN")
                           and float(r[c]) > cut for c, cut in cuts.items()):
                        continue
                except ValueError:
                    continue
                take(*(r.get(c) for c in REPORT_PEPTIDE_COLUMNS))
    return dict(groups), dict(runs), qcols, None


def searched_fasta(meta, meta_path):
    """(path, None) of the FASTA the sidecar describes, checked by its sha256 when the sidecar has
    one -- or (None, why not)."""
    want = (meta or {}).get("sha256")
    cands = list(dict.fromkeys(p for p in ((meta or {}).get("fasta"),
                                           meta_path[:-len(".meta.json")] if meta_path and
                                           meta_path.endswith(".meta.json") else None) if p))
    changed = []
    for p in cands:
        if os.path.isfile(p):
            if want:
                import hashlib
                h = hashlib.sha256()
                with open(p, "rb") as fh:
                    for chunk in iter(lambda: fh.read(1 << 20), b""):
                        h.update(chunk)
                if h.hexdigest() != want:
                    changed.append(p)
                    continue
            return p, None
    names = ", ".join(dict.fromkeys(os.path.basename(p) for p in cands))
    if changed:
        return None, (f"the FASTA ({', '.join(dict.fromkeys(os.path.basename(p) for p in changed))}) "
                      f"has changed since the search (its sha256 is not the sidecar's "
                      f"{want[:12]}...)")
    return None, (f"the searched FASTA ({names or 'no path in the sidecar'}"
                  f"{', sha256 ' + want[:12] + '...' if want else ''}) is not readable here")


def audit_contaminant_evidence(findings, de_prov, meta_path):
    """Contaminant protein groups that rest only on peptides of the sample's own proteins
    (fetch_fasta.contaminant_groups_on_target_peptides: shared with one, or the same with I = L)
    in the report the DE read. The build keeps a contaminant entry with peptides of its own; which
    of them a sample shows only the report can say (horse myoglobin in human heart, 2026-10-14)."""
    cont = (de_prov or {}).get("contaminants")
    inp = (de_prov or {}).get("input")
    if not inp or not isinstance(cont, dict) or not cont.get("n_protein_groups"):
        return                      # no DE record, or no contaminant group in it: nothing to check
    meta_path = meta_path or cont.get("fasta_meta")
    try:
        with open(meta_path) as fh:
            meta = json.load(fh)
    except (OSError, TypeError, ValueError):
        meta = None
    if not isinstance(meta, dict):
        add(findings, "contaminant_evidence", "INFO",
            "Not assessed: no readable <fasta>.meta.json (pass --fasta-meta), so whether a "
            "contaminant group rests only on the sample's own peptides cannot be checked.", plain=NOT_CHECKED)
        return
    fasta, why = searched_fasta(meta, meta_path)
    groups = runs = None
    if fasta:
        q = de_prov.get("q_cutoff")
        from diann_q_columns import DEFAULT_Q_CUTOFF
        groups, runs, qcols, why = report_contaminant_peptides(
            inp, float(q) if isinstance(q, (int, float)) else DEFAULT_Q_CUTOFF,
            de_prov.get("q_columns") or None)
    if groups is None:
        add(findings, "contaminant_evidence", "INFO",
            f"Not assessed: {why}, so whether a contaminant group rests only on the sample's own "
            "peptides cannot be checked.", plain=NOT_CHECKED)
        return
    if not groups:
        # the DE counted contaminant groups and this reading of its report finds none: something
        # differs (the report, its filter) -- never a PASS over nothing (error-handling review)
        add(findings, "contaminant_evidence", "INFO",
            f"Not assessed: the DE record counts {cont.get('n_protein_groups')} contaminant protein "
            f"group(s), but none passes the DE's q-value filter in {os.path.basename(inp)} as read "
            f"here (q-columns {', '.join(qcols) or 'none'}), so their peptides could not be checked.",
            {"report": inp, "q_columns": qcols}, plain=NOT_CHECKED)
        return
    # A target's peptides are those the build's digest makes from it (the sidecar's
    # contaminant_digest: the rule that kept these entries), never any string inside it
    from estimate_params import DIANN_DIGEST
    dig = meta.get("contaminant_digest")
    digest, digest_src = ((dig, "the sidecar's contaminant_digest") if isinstance(dig, dict) and
                          {"cut", "missed_cleavages", "min_pep_len", "max_pep_len"} <= set(dig)
                          else (DIANN_DIGEST, "estimate_params.DIANN_DIGEST (the sidecar records "
                                              "no contaminant_digest)"))
    with open(fasta, encoding="utf-8", errors="replace") as fh:
        flagged = contaminant_groups_on_target_peptides(groups, fh.read(), digest)
    org = meta.get("organism") or "target-organism"
    detail = {"report": inp, "fasta": fasta, "q_columns": qcols, "n_groups_checked": len(groups),
              "digest": dict(digest), "digest_source": digest_src,
              "groups": [dict(g, n_runs=len(runs.get(g["group"], ()))) for g in flagged]}
    if not flagged:
        add(findings, "contaminant_evidence", "PASS",
            f"Each of the {len(groups)} contaminant protein group(s) in "
            f"{os.path.basename(inp)} has a peptide that is not a {org} protein's (I = L).",
            detail)
        return

    def pep(p):
        who = f"{p['target_acc']}{' ' + p['gene'] if p.get('gene') else ''}"
        return (f"{p['peptide']} differs from {who}'s {p['target_peptide']} only by I/L (the "
                f"same mass and fragments)"
                if p["how"] == "il_twin" else f"{p['peptide']} is also {who}'s")

    def group(g):
        ps = g["peptides"]
        return (f"{g['group']} ({', '.join(e for e in g['cont_entries'] if e) or 'contaminant'}; "
                f"{len(ps)} peptide{'s' if len(ps) != 1 else ''} in "
                f"{len(runs.get(g['group'], ()))} run(s)): "
                + ", ".join(pep(p) for p in ps[:5]) + (f", +{len(ps) - 5} more" if len(ps) > 5 else ""))
    genes = sorted({p.get("gene") or p["target_acc"] for g in flagged for p in g["peptides"]})
    # the contaminant accessions of the flagged groups: what fetch_fasta.py --drop-contaminant takes
    drop = list(dict.fromkeys(t.strip() for g in flagged for t in g["group"].split(";")
                              if is_contaminant_accession(t)))
    # whose filter it was, in the DE record's own words (rule 1: never one pipeline's name for
    # all): its rule_origin; never a "[... -- confirm]" tag cut out of the rule text
    owner = str(cont.get("rule_origin") or "").strip() or "the DE's contaminant filter"
    fate = (f"Their precursors were removed with the contaminants ({owner}; DIA-NN's "
            f"--cont-quant-exclude too, when set), so {', '.join(genes[:12])} may be "
            f"under-quantified and the contaminant share overstated."
            if cont.get("removed") else
            f"They were quantified and tested as contaminant groups (--keep-contaminants), beside "
            f"{', '.join(genes[:12])}.")
    # What the data can say: these peptides are the sample organism's, so the group is not
    # evidence of the contaminant's species. It may be the sample's own protein, or handling
    # contamination of that organism (human skin keratin behind a sheep wool keratin group) -- the
    # data cannot tell which (2.12 statistics review), so neither is asserted.
    add(findings, "contaminant_evidence", "WARN",
        f"{len(flagged)} contaminant protein group(s) rest only on peptides of {org} proteins, so "
        f"the data cannot attribute them to the contaminant's species: "
        + "; ".join(group(g) for g in flagged[:12]) + (f"; +{len(flagged) - 12} more" if len(flagged) > 12 else "")
        + f". {fate} Report each as the {org} protein -- the sample's own, or handling "
          f"contamination (skin and hair keratins) -- never as a signal of the contaminant species "
          f"or an added standard. To quantify it as the {org} protein, rebuild the FASTA with the "
          f"same fetch_fasta.py command plus "
        + " ".join(f"--drop-contaminant {a}" for a in drop)
        + " and search again.", dict(detail, drop_contaminant=drop))


def audit_pipeline_checks(findings, prov):
    """The DE pipeline's own checks, as it recorded them (de_provenance.json `checks`: {check,
    status, message, detail}) -- quoted; this audit knows no pipeline by name (run_silac.R's channel
    balance, label incorporation and declaration are its own, for example)."""
    for c in prov.get("checks") or []:
        if isinstance(c, dict) and c.get("check") and c.get("status") in ("PASS", "WARN", "FAIL",
                                                                         "INFO") and c.get("message"):
            # its own plain words and label (run_silac.R chk(plain=, label=)): this audit
            # writes none for a pipeline it does not know
            add(findings, str(c["check"]), c["status"], str(c["message"]),
                c.get("detail") if isinstance(c.get("detail"), dict) else {},
                plain=str(c.get("plain") or "") or None, label=str(c.get("label") or "") or None)


def audit_de(findings, de_dir, adjp, logfc, prov=None):
    for f in sorted(glob.glob(os.path.join(de_dir, "DE_*.csv"))):
        rows = read_csv(f)
        ct = os.path.basename(f)
        n = sig = beyond = 0
        for r in rows:
            try:
                a = float(r.get("adj.P.Val", "nan")); l = float(r.get("logFC", "nan"))
            except ValueError:
                continue
            if a != a:  # nan
                continue
            n += 1
            # Significance is the adjusted p-value alone -- the same rule run_de.R and the
            # figures use. |logFC| is only counted descriptively, for the report.
            if a < adjp:
                sig += 1
                if abs(l) >= logfc:
                    beyond += 1
        if n == 0:
            continue
        frac = sig / n
        if sig == 0:
            add(findings, "de_signal", "WARN",
                f"{ct}: 0 proteins significant (adj.P<{adjp}). The experiment may be "
                "underpowered, the effect small, or the groups mislabeled.",
                plain=PLAIN_DE_SIGNAL_NONE)
        elif frac > 0.5:
            note = (prov or {}).get("high_significance_note")
            add(findings, "de_signal", "WARN",
                f"{ct}: {sig}/{n} ({frac*100:.0f}%) proteins significant — " + (
                    f"a high share. {note}" if note else
                    "implausibly high. This usually means a batch effect, a normalization problem, or "
                    "confounded/mislabeled groups, not real biology."), plain=PLAIN_DE_SIGNAL_HIGH)
        else:
            add(findings, "de_signal", "PASS",
                f"{ct}: {sig}/{n} proteins significant at adj.P<{adjp} ({frac*100:.0f}%); "
                f"{beyond} of those are also ≥{2**logfc:.3g}-fold.")


def load_search_prov(de_prov, search_prov=None):
    """(path, search_provenance.json) -- --search-prov, else the one beside de_provenance.json's
    `input` -- or (path, None) when there is none to read."""
    inp = (de_prov or {}).get("input")
    path = search_prov or (os.path.join(os.path.dirname(inp), "search_provenance.json")
                           if inp else None)
    try:
        with open(path) as fh:
            return path, json.load(fh)
    except (OSError, TypeError, ValueError):
        return path, None


def audit_search_measurement(findings, de_prov, search_prov=None):
    """A search that was to measure the scan window / mass accuracy on these runs and fell back
    instead (probe_fallback.py): the CAUTION, in probe_fallback.caution()'s one wording."""
    path, sprov = load_search_prov(de_prov, search_prov)
    if not isinstance(sprov, dict):
        return                      # no search record here: not a question this audit can ask
    from probe_fallback import fallback_modes, measurement_note, measurement_warnings, \
        provenance_caution
    msg = provenance_caution(sprov)
    if msg:
        window, mass_acc = fallback_modes(sprov)
        fb = sprov.get("probe_fallback") if isinstance(sprov.get("probe_fallback"), dict) else {}
        add(findings, "search_measurement", "WARN", msg,
            {"scan_window_mode": window, "mass_acc_mode": mass_acc, "reason": fb.get("reason"),
             "search_provenance": path}, plain=P_MEASURE_FALLBACK)
    # ...and a measured mass accuracy pinned at the SOP floor over a spread the probe warned about
    note = measurement_note(sprov)
    if note:
        add(findings, "search_measurement", "WARN", note,
            {"warnings": measurement_warnings(sprov), "search_provenance": path},
            plain="The runs disagreed about the instrument's mass accuracy, and the search used "
                  "the standard setting for it.")


def audit_normalisation(findings, de_prov, rows=None, search_prov=None):
    """DIA-NN's Normalisation.Instability for the runs the DE used (check_report_runs.
    normalisation_instability: the one reader, threshold and its basis). The DE quantities are
    DIA-NN's normalised ones (Precursor.Normalised, PG.MaxLFQ), so a normalisation DIA-NN itself
    calls unstable reaches every result -- and nothing said so (a staff report: 0.67-1.00 in
    every run of a cohort, passed silently into the DE)."""
    inp = (de_prov or {}).get("input")
    if not inp:
        return                      # no DE record of its input: nothing to read
    if (de_prov or {}).get("quantities") == "raw":
        add(findings, "normalisation_stability", "INFO",
            "Not relevant: the DE read DIA-NN's non-normalised quantities (de_provenance.json "
            "quantities: raw), so DIA-NN's cross-run normalisation does not enter it.",
            plain="Not relevant here: the analysis did not use DIA-NN's normalised amounts.")
        return
    from check_report_runs import (normalisation_instability, stats_path,
                                   NORM_INSTABILITY_COLUMN)
    runs = [r["File.Name"] for r in rows or [] if (r.get("File.Name") or "").strip()]
    ni = normalisation_instability(inp, runs or None)
    if ni is None:
        _, sprov = load_search_prov(de_prov, search_prov)
        if isinstance(sprov, dict) and (sprov.get("engine") or "").lower() == "diann":
            add(findings, "normalisation_stability", "INFO",
                f"Not assessed: there is no {os.path.basename(stats_path(str(inp)))} beside "
                f"{os.path.basename(str(inp))} (a cfg with --no-stats?), so DIA-NN's per-run "
                f"{NORM_INSTABILITY_COLUMN} could not be read.", plain=P_STAB_UNCHECKED)
        return                      # another engine's report: DIA-NN's question does not apply
    where = os.path.basename(ni["stats_file"])
    if not ni["column"]:
        add(findings, "normalisation_stability", "INFO",
            f"Not assessed: {where} has no {NORM_INSTABILITY_COLUMN} column (a DIA-NN version "
            "that does not write it).", {"stats_file": ni["stats_file"]}, plain=P_STAB_UNCHECKED)
        return
    if not ni["per_run"]:
        add(findings, "normalisation_stability", "INFO",
            f"Not assessed: none of the analysed runs has identified precursors in {where}.",
            {"stats_file": ni["stats_file"]}, plain=P_STAB_UNCHECKED)
        return
    detail = {k: ni[k] for k in ("stats_file", "per_run", "median", "max", "flagged",
                                 "threshold", "basis")}
    n, thr = len(ni["per_run"]), ni["threshold"]
    if not ni["flagged"]:
        add(findings, "normalisation_stability", "PASS",
            f"DIA-NN's {NORM_INSTABILITY_COLUMN} is at most {thr:g} in all {n} runs (median "
            f"{ni['median']:.2f}, max {ni['max']:.2f}; {where}).", detail)
        return
    runs_txt = ", ".join(f"{k} {v:.2f}" for k, v in sorted(ni["flagged"].items(),
                                                            key=lambda kv: -kv[1]))
    add(findings, "normalisation_stability", "WARN",
        f"DIA-NN's {NORM_INSTABILITY_COLUMN} is above {thr:g} in {len(ni['flagged'])} of {n} runs "
        f"({runs_txt}; median of all {n}: {ni['median']:.2f}; {where}). DIA-NN's own "
        "normalisation was unstable for these runs, and the DE used its normalised quantities "
        "(Precursor.Normalised / PG.MaxLFQ), so a distortion of that normalisation is in every "
        f"result. ({ni['basis']}.) The alternatives DIA-NN documents: re-quantify with "
        "--global-norm (one global factor per run instead of RT-dependent normalisation) or "
        "--no-norm (no cross-run normalisation: Precursor.Normalised and PG.MaxLFQ are then "
        "not normalised). Either can reuse this search's .quant files: DIA-NN's README lists "
        "cross-run normalisation among the settings that may differ when .quant files are "
        "reused. The "
        "report's Precursor.Quantity is DIA-NN's non-normalised quantity. Run the DE both "
        "ways and compare before trusting either.", detail)


def audit_normalization_check(findings, de_prov):
    """The normalisation decision, quoted from run_de.R's record (normalization_check.py's block
    in de_provenance.json): which quantities, the experiment type's default, whether the data
    check agreed, and who chose otherwise -- never re-judged here."""
    nc = (de_prov or {}).get("normalization_check")
    if not isinstance(nc, dict):
        return                      # a DE record from before 2.10: nothing to quote
    from normalization_check import not_applicable_reason
    na_reason = not_applicable_reason(nc)
    if na_reason:
        # the pipeline's own record of why there is no between-run normalisation to check
        add(findings, "normalization_check", "PASS",
            f"Not applicable to {(de_prov or {}).get('display_label') or 'this pipeline'}: "
            f"{na_reason}", {"status": "not_applicable"},
            plain="Not applicable to this kind of analysis.")
        return
    if nc.get("status") != "decided":
        # a WARN, and core_submission.py deliver holds the delivery (normalization_check.
        # delivery_gate): an IP read with normalised quantities must not pass as a footnote
        add(findings, "normalization_check", "WARN",
            ("This DE is one of the normalisation check's inputs, not the final DE" if
             nc.get("status") == "check_input" else "The normalisation check was not run")
            + f": {nc.get('note') or 'no normalisation check'} (quantities: "
            f"{nc.get('quantities_applied')}). Run step 8's check and the final DE with "
            "--normalization-check before delivering.", {"status": nc.get("status")},
            plain="The choice of normalisation was not yet checked against the data.")
        return
    et_ = nc.get("experiment_type") or {}
    d = nc.get("decision") or {}
    # AUDIT.md is delivered: a person's decision by role (normalization_check.statement), never
    # by name -- the name is in the staff-only record
    chose = d.get("statement") or (f"{d.get('by')} chose {d.get('quantities')} quantities: "
                                   f"{d.get('reason')}.")
    detail = {"experiment_type": et_.get("type"), "default": (nc.get("default") or {}).get("quantities"),
              "applied": nc.get("quantities_applied"), "trips": nc.get("trips"),
              "decision_by": d.get("by")}
    if not et_.get("type"):
        # no experiment type recorded: nothing to agree with, and never a PASS
        add(findings, "normalization_check", "WARN",
            "No experiment type was recorded, so no default applied"
            + (f" (the data check: {'; '.join(nc.get('trips') or [])})" if nc.get("trips") else "")
            + f"; {chose}", detail,
            plain=f"No experiment type was recorded, so there was no usual normalisation to "
                  f"compare with. {chose}")
    elif nc.get("tripped"):
        add(findings, "normalization_check", "WARN",
            f"Experiment type {et_.get('label')} ({et_.get('source')}): the data check did not "
            f"agree with its default ({'; '.join(nc.get('trips') or [])}); {chose}", detail,
            plain=f"The usual normalisation for this kind of experiment did not fit these "
                  f"samples, so the Core chose which quantities to use. {chose}")
    elif d.get("source") == "user" and nc.get("quantities_applied") != detail["default"]:
        # the data agreed with the default, and a person chose otherwise: their decision, said
        add(findings, "normalization_check", "WARN",
            f"Experiment type {et_.get('label')} ({et_.get('source')}): the data check agreed "
            f"with its default ({detail['default']} quantities), and "
            f"{nc.get('quantities_applied')} quantities were used; {chose}", detail,
            plain=f"The usual normalisation for this kind of experiment fitted these samples, and "
                  f"the Core chose other quantities. {chose}")
    else:
        add(findings, "normalization_check", "PASS",
            f"Experiment type {et_.get('label')} ({et_.get('source')}): "
            f"{nc.get('quantities_applied')} quantities, its default; the data check agreed"
            + (f"; {chose}" if d.get("source") == "user" else "."),
            detail, plain=P_NORM_AGREED)


def audit_passes(findings, de_prov, search_prov=None):
    """The chain's first pass vs its final pass, as step 5 recorded it (pass_comparison.py)."""
    inp = (de_prov or {}).get("input")
    path, sprov = load_search_prov(de_prov, search_prov)
    if sprov is None:
        return                      # no search record here: not a question this audit can ask
    if sprov.get("search_mode") != "parallel_5step":
        return                      # single-shot: DIA-NN's own MBR, no separate first pass
    from diann_parallel import FIRST_PASS_REPORT
    used_first = os.path.basename(str(inp or "")) == FIRST_PASS_REPORT
    pc = sprov.get("pass_comparison")
    if not isinstance(pc, dict) or not isinstance(pc.get("runs"), list):
        add(findings, "first_vs_final_pass", "INFO",
            "the 5-step chain's first-pass / final-pass comparison is not in "
            f"{os.path.basename(path)} (a chain generated before skill 2.9, or step 5 did not "
            f"get that far): compare {FIRST_PASS_REPORT[:-len('.parquet')]}.stats.tsv with "
            "report.stats.tsv by hand before trusting the final pass", plain=P_PASSES_UNRECORDED)
        return
    flagged = [r for r in pc["runs"] if r.get("flagged")]
    better = pc.get("switch_recommended_for")          # None: post-filter rows not counted
    if not flagged and not better:
        add(findings, "first_vs_final_pass", "PASS",
            f"no run lost more than half its first-pass precursors in the final pass "
            f"({pc.get('n_runs')} runs compared)", pc)
        return
    caveat = ("counts before protein-group FDR; compare the rows that pass the q-value filter "
              "before switching")
    lost = ", ".join(
        f"{r['run']} {r['first']['precursors']} -> {r['final']['precursors']}"
        + ("" if r["first"].get("rows_after_q") is None else
           f" ({r['first']['rows_after_q']} vs {r['final']['rows_after_q']} rows after the q "
           "filter)") for r in flagged)
    what = (f"{len(flagged)} of {pc.get('n_runs')} runs kept less than half their first-pass "
            f"precursors in the final pass ({lost}; {caveat})" if flagged else
            f"no run lost more than half its first-pass precursors ({caveat})")
    if better is None:
        verdict = (" The rows after the q-value filter were not counted (no pyarrow where step 5 "
                   "ran): count them with pass_comparison.py before offering the first-pass "
                   f"report ({FIRST_PASS_REPORT}).")
    elif better:
        verdict = (f" After the q-value filter the first pass has more rows for "
                   f"{', '.join(better)}: offer the user the first-pass report "
                   f"({FIRST_PASS_REPORT}: each run searched once, no match-between-runs) as the "
                   "DE input -- DIA-NN's README: use the first pass when it performs better.")
    else:
        verdict = (" After the q-value filter the first pass has no more rows than the final pass "
                   "for any run, so switching would not recover them.")
    if used_first:
        add(findings, "first_vs_final_pass", "INFO",
            f"{what}. The DE used the first-pass report ({FIRST_PASS_REPORT}), so the results "
            "do not depend on the second pass.", pc,
            plain=P_PASSES_USED_FIRST)
    else:
        add(findings, "first_vs_final_pass", "WARN", what + "." + verdict, pc,
            plain=P_PASSES_LOST if flagged else P_PASSES_FIRST_MORE)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default="AUDIT.md")
    ap.add_argument("--conditions")
    ap.add_argument("--de-dir")
    ap.add_argument("--acquisition-json")
    # Default: the cutoffs run_de.R recorded in de_provenance.json -- the one definition the
    # tables, figures and methods use. 0.05 / 1 only when there is no record, and said so.
    ap.add_argument("--adjp", type=float, default=None)
    ap.add_argument("--logfc", type=float, default=None)
    ap.add_argument("--min-proteins", type=int, default=500)
    ap.add_argument("--max-missing", type=float, default=0.5)
    ap.add_argument("--keratin-sample", action="store_true",
                    help="sample IS keratin (nail/hair/wool/skin/feather) — keratin is the analyte, "
                         "not a contaminant; do not flag KRT/KRTAP/keratin")
    ap.add_argument("--search-prov",
                    help="run_search.py's search_provenance.json, for the 5-step chain's "
                         "first-pass / final-pass comparison. Default: beside "
                         "de_provenance.json's `input`")
    ap.add_argument("--fasta-meta",
                    help="fetch_fasta.py's <fasta>.meta.json -- the list of proteins that are also "
                         "common-contaminant sequences. Default: ./search.fasta.meta.json if present")
    ap.add_argument("--search-out",
                    help="the search's output directory: a Sage search's LFQ mass-window check "
                         "(sage_lfq_check.json) becomes the sage_lfq finding")
    a = ap.parse_args()
    fasta_meta = a.fasta_meta or ("search.fasta.meta.json"
                                  if os.path.exists("search.fasta.meta.json") else None)

    prov = {}
    if a.de_dir:
        try:
            with open(os.path.join(a.de_dir, "de_provenance.json")) as fh:
                prov = json.load(fh)
        except (OSError, ValueError):
            prov = {}
    cutoff_source = {}
    for key, default in (("adjp", 0.05), ("logfc", 1.0)):
        if getattr(a, key) is not None:
            cutoff_source[key] = "command line"
        elif isinstance(prov.get(key), (int, float)):
            setattr(a, key, float(prov[key])); cutoff_source[key] = "de_provenance.json"
        else:
            setattr(a, key, default)
            cutoff_source[key] = "DEFAULT -- not recorded in de_provenance.json"

    findings = []
    if a.conditions and os.path.exists(a.conditions):
        audit_conditions(findings, read_csv(a.conditions), a.adjp, block_col=prov.get("block_column"))
    if a.acquisition_json and os.path.exists(a.acquisition_json):
        audit_acquisition(findings, a.acquisition_json)
    # A keratin sample: --keratin-sample, or a database built with fetch --keratin-sample (its
    # sidecar says so), so the step-8c flag is no longer the only way the audit learns it.
    keratin = a.keratin_sample
    if fasta_meta and not keratin:
        try:
            with open(fasta_meta) as fh:
                keratin = keratin_sample_recorded(json.load(fh)) is True
        except (OSError, ValueError):
            pass            # audit_target_contaminants reports the unreadable sidecar
    audit_quantification(findings, prov)
    audit_pipeline_checks(findings, prov)
    audit_numeric_covariates(findings, prov)
    if a.search_out and os.path.isdir(a.search_out):
        audit_search(findings, a.search_out)
    elif a.search_out:
        add(findings, "search_checks", "INFO", f"--search-out {a.search_out} is not a directory, "
            "so the search's own checks (Sage LFQ window) were not read.")
    em = os.path.join(a.de_dir, "Expression_Matrix.csv") if a.de_dir else None
    if em and os.path.exists(em):
        audit_matrix(findings, em, a.min_proteins, a.max_missing, keratin)
    if a.de_dir and os.path.isdir(a.de_dir):
        audit_target_contaminants(findings, fasta_meta, em, a.de_dir, a.adjp, keratin)
        audit_contaminant_evidence(findings, prov, fasta_meta)
        audit_de(findings, a.de_dir, a.adjp, a.logfc, prov)
    audit_search_measurement(findings, prov, a.search_prov)
    audit_passes(findings, prov, a.search_prov)
    audit_normalization_check(findings, prov)
    audit_normalisation(findings, prov,
                        read_csv(a.conditions) if a.conditions and os.path.exists(a.conditions)
                        else None, a.search_prov)

    n_fail = sum(1 for f in findings if f["status"] == "FAIL")
    n_warn = sum(1 for f in findings if f["status"] == "WARN")
    overall = "FAIL" if n_fail else ("WARN" if n_warn else "PASS")

    # INFO = a check that could not run; it never changes the overall status.
    icon = {"PASS": "✅", "WARN": "⚠️", "FAIL": "⛔", "INFO": "ℹ️"}
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
        json.dump({"overall": overall, "n_fail": n_fail, "n_warn": n_warn, "findings": findings,
                   "cutoffs": {"adjp": a.adjp, "logfc": a.logfc, "source": cutoff_source}}, fh, indent=2)

    print(json.dumps({"overall": overall, "n_fail": n_fail, "n_warn": n_warn,
                      "report": os.path.abspath(a.out),
                      "blocking": [f["message"] for f in findings if f["status"] == "FAIL"],
                      "warnings": [f["message"] for f in findings if f["status"] == "WARN"]}, indent=2))


if __name__ == "__main__":
    main()
