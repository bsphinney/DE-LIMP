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

Each curated panel is ONE list of human gene symbols (PANELS), matched case-insensitively,
plus the ortholog names that differ by more than case for the searched organism
(PANEL_ORTHOLOGS: mouse/rat adult globins, verified against MGI/RGD and the reference
proteomes; the organism comes from --fasta-meta's taxid or --taxid). A panel that matches
nothing is reported as "check could not run" when NO panel matched anything -- the names
did not match -- rather than as an absence. Cont_-tagged proteins (common-contaminant
entries, e.g. bovine serum haemoglobin from an antibody prep) are listed as "contaminant,
not sample" and never scored; so are the ones run_de.R already removed (read from its
de_provenance.json record beside the matrix).

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
      [--taxid 10090]                          # organism, when there is no sidecar

Reads/writes plain files; stdlib only (pyarrow optional, just for --report).
"""
import sys, os, csv, json, math, argparse, random, re

# The list's ONE definition is the FASTA sidecar; the wording for lost proteins lives there too.
from fetch_fasta import (target_contaminants, seen_only_as_cont, lost_to_contaminants_message,
                         CONT_TAG)

# Curated marker panels, as HUMAN gene symbols (matched case-insensitively against every
# id column of the matrix). Deliberately specific markers -- avoid ubiquitous glycolytic
# enzymes. Extend per tissue as needed; document additions in references/anomaly-checks.md.
# LORICRIN: HGNC's current symbol for loricrin, and the GN= UniProt uses for human, mouse
# and rat (UP000005640 / UP000000589 / UP000002494, 2026_03) -- "LOR" alone never matched.
PANELS = {
    "HEMOLYSIS": ["HBA1", "HBA2", "HBA", "HBB", "HBD", "CA1", "CA2", "CAT", "PRDX2",
                  "BLVRB", "SPTA1", "SPTB", "ANK1", "SLC4A1", "EPB42", "PKLR", "BPGM",
                  "ALAS2", "AHSP", "CATALASE"],
    "SKELETAL_MUSCLE": ["MYH1", "MYH2", "MYH7", "MYBPC1", "MYBPC2", "ACTN2", "ACTN3",
                        "CKM", "MB", "TNNI1", "TNNI2", "TNNT3", "TNNC2", "PYGM", "ENO3",
                        "TPM1", "TPM2", "MYOM1", "MYOM2", "MYOM3", "DES", "MYL1", "MYL2",
                        "ATP2A1", "NEB", "TTN", "CASQ1", "SLN"],
    "EPIDERMIS": ["KRT1", "KRT2", "KRT9", "KRT10", "KRT5", "KRT14", "KRT16", "KRT6A",
                  "FLG", "LOR", "LORICRIN", "IVL", "DSP", "JUP", "SBSN"],
}

# Ortholog names that differ from the human symbol by MORE than letter case, per NCBI taxid.
# The Genes column DIA-NN writes is the GN= of the searched FASTA, so these are the gene
# names UniProt gives the orthologs in each reference proteome, plus the MGI/RGD symbols
# of the same genes (a full proteome or a later release may carry either).
# Why a table and not a prefix rule like HBB-*/HBA-*: a prefix also pulls in the embryonic
# chains (Hbb-y, Hbb-bh0/1/2, Hba-x -- embryonic erythropoiesis, not lysis of adult blood),
# still misses rat Hbbl1 and HBB2_RAT (which has no GN= at all), and cannot express
# anything that is not a prefix. A table is exact, and a test re-reads it.
# Verified 2026-09-24 against the MGI HOM_MouseHumanSequence.rpt and RGD_ORTHOLOGS.txt
# ortholog tables and the UP000000589 (mouse) / UP000002494 (rat) one-per-gene FASTAs:
# every other panel gene's ortholog IS the human symbol in title case (Ca1, Cat, Prdx2,
# Myh1, Ckm, Krt10, ...), so case-insensitive matching already covers it. Only the adult
# globins differ. (msalemi 2026-09-24, mouse brain: Hbb-bs at ~15 log2 matched nothing.)
PANEL_ORTHOLOGS = {
    10090: {   # Mus musculus. UP000000589 GN=: Hba (P01942), Hbb-b1, Hbb-b2, Hbb-bs (A8DUK4)
        "HBA1": ("Hba", "Hba-a1", "Hba-a2"),
        "HBA2": ("Hba", "Hba-a1", "Hba-a2"),
        "HBB": ("Hbb-b1", "Hbb-b2", "Hbb-bs", "Hbb-bt"),
        "HBD": ("Hbb-b1", "Hbb-b2", "Hbb-bs", "Hbb-bt"),
    },
    10116: {   # Rattus norvegicus. UP000002494 GN=: Hba1 (P01946), Hbb (P02091); P11517
               # HBB2_RAT has no GN=, so only its accession can match. Full proteome adds
               # Hba-a1/-a2/-a3, Hbb-b1, Hbb-bs, Hbbl1.
        "HBA1": ("Hba-a1", "Hba-a2", "Hba-a3"),
        "HBA2": ("Hba1", "Hba-a1", "Hba-a2", "Hba-a3"),
        "HBB": ("Hbb-b1", "Hbb-b2", "Hbb-bs", "Hbb-bt", "Hbbl1", "P11517"),
        "HBD": ("Hbb", "Hbb-b1", "Hbb-b2", "Hbb-bs", "Hbb-bt", "Hbbl1", "P11517"),
    },
}
# Organisms whose names the panels are verified for: human symbols, and PANEL_ORTHOLOGS.
PANEL_ORGANISMS = {9606: "Homo sapiens", 10090: "Mus musculus", 10116: "Rattus norvegicus"}


def panel_genes(name, taxid=None):
    """The ONE place a panel's match set is built: its human symbols plus the organism's
    ortholog names. taxid None (organism unknown) -> every organism's names; they cannot
    collide, since no mouse/rat globin name is a human symbol."""
    human = PANELS[name]
    tables = ([PANEL_ORTHOLOGS[taxid]] if taxid in PANEL_ORTHOLOGS else
              [] if taxid is not None else list(PANEL_ORTHOLOGS.values()))
    names = set(human)
    for t in tables:
        for g in human:
            names.update(t.get(g, ()))
    return {n.upper() for n in names}


# How a Cont_-tagged marker is reported: it measures a common-contaminant entry, not the sample.
CONTAMINANT_NOT_SAMPLE = ("contaminant, not sample (a common-contaminant entry: reagent, serum or "
                          "handling -- e.g. bovine serum haemoglobin from an antibody prep); not "
                          "counted in the panel score")
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
            rows.append({"genes": genes, "ids": ids, "vals": vals,
                         "cont": _cont_accessions(ids)})
        samples = [header[i] for i in sample_idx]
    return samples, rows


def _cont_accessions(ids):
    """The Cont_-tagged accessions of a row (upper-cased tokens). A protein group naming one
    is a common-contaminant entry -- bovine serum haemoglobin from an antibody prep, say --
    not the sample's own protein, so it never counts toward a sample panel score."""
    return sorted(CONT_TAG + t[len(CONT_TAG):] for t in ids if t.startswith(CONT_TAG.upper()))


def removed_contaminant_rows(matrix_path):
    """Contaminant groups run_de.R removed before building the matrix (its record in
    de_provenance.json names the table), as panel rows -- so a Cont_ haemoglobin is still
    reported as "contaminant, not sample" after the filter took it out of the matrix."""
    d = os.path.dirname(os.path.abspath(matrix_path))
    try:
        with open(os.path.join(d, "de_provenance.json")) as fh:
            rec = json.load(fh).get("contaminants") or {}
        if not rec.get("removed_table"):
            return []
        with open(os.path.join(d, rec["removed_table"]), newline="") as fh:
            out = []
            for r in csv.DictReader(fh):
                if (r.get("Contaminant.Group") or "").upper() != "TRUE":
                    continue
                ids = _genes(r.get("Protein.Group"))
                out.append({"genes": _genes(r.get("Genes")) | ids, "ids": ids, "vals": {},
                            "cont": _cont_accessions(ids), "removed": True})
            return out
    except (OSError, ValueError):
        return []


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


def panel_scores(rows, samples, panel_genes, removed=()):
    """-> (per-sample mean, n sample proteins, matched names, contaminant matches). A
    Cont_-tagged row is reported under contaminant matches and never scored: it measures
    the reagent, not the sample. `removed` = rows run_de.R already filtered out."""
    pg = set(g.upper() for g in panel_genes)
    hits = [r for r in rows if r["genes"] & pg and not r.get("cont")]
    per = {}
    for s in samples:
        vals = [r["vals"].get(s) for r in hits]
        vals = [v for v in vals if v is not None]
        per[s] = (sum(vals) / len(vals)) if vals else None
    matched = sorted({g for r in hits for g in (r["genes"] & pg)})
    cont = sorted({f"{'/'.join(sorted(r['genes'] & pg))} ({';'.join(r['cont'])})"
                   + (" -- removed before DE" if r.get("removed") else "")
                   for r in list(rows) + list(removed) if r.get("cont") and r["genes"] & pg})
    return per, len(hits), matched, cont


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


# Is a panel score confounded with group? A permutation test of the one-way between-group
# F statistic: how often does a random relabelling of the samples (same group sizes) put at
# least as much of the panel's variation between groups as the real labels do?
# Why not the old rule (flag when the highest and lowest GROUP MEANS differ by >= --z SD and
# barely overlap): with k groups it compares the extremes of k noisy means, and their spread
# grows with k even when nothing differs -- the expected range of 10 means is ~3.1 standard
# errors against ~1.1 for 2. Silva08172026 (10 groups x 3) flagged all three panels that way
# (gaps 1.55-2.74). The test's null distribution comes from the design itself, so k, group
# sizes and non-normal scores are all accounted for. 1% because three panels are tested.
CONFOUND_P = 0.01
N_PERM = 9999


def _between(values, labels):
    """sum over groups of (group sum)^2 / group size -- monotone in the one-way F for fixed
    data, so permutations can be compared on it directly."""
    sums, sizes = {}, {}
    for v, g in zip(values, labels):
        sums[g] = sums.get(g, 0.0) + v
        sizes[g] = sizes.get(g, 0) + 1
    return sum(sums[g] ** 2 / sizes[g] for g in sums)


def _min_p(sizes):
    """Smallest p-value any relabelling can give: 1 / number of distinct partitions."""
    labellings = math.factorial(sum(sizes))
    for n in sizes:
        labellings //= math.factorial(n)
    sym = 1
    for n in set(sizes):
        sym *= math.factorial(sizes.count(n))
    return sym / labellings


def confound_check(z, gmap, samples, thr, n_perm=N_PERM, seed=20260925):
    """Return (confounded: bool, detail, p) -- does the panel score differ between groups
    more than relabelling the samples at random would? That is the danger signal: DE may be
    contamination, not biology. p is None when the design is too small for any relabelling
    to reach CONFOUND_P (3 vs 3: 1 in 10); then complete separation of the extreme groups
    by >= thr is reported instead, labelled as such."""
    if not gmap:
        return False, "no conditions.csv -- group-confounding not assessed", None
    vals, labs = [], []
    for s in samples:
        if s in gmap and z.get(s) is not None:
            vals.append(z[s])
            labs.append(gmap[s])
    groups = {}
    for v, g in zip(vals, labs):
        groups.setdefault(g, []).append(v)
    if len(groups) < 2:
        return False, "fewer than 2 groups with panel data", None
    means = {g: sum(v) / len(v) for g, v in groups.items()}
    detail = "; ".join(f"{g}: mean z {means[g]:+.2f} (n={len(groups[g])})" for g in sorted(means))
    sizes = [len(v) for v in groups.values()]
    if _min_p(sizes) > CONFOUND_P:
        hi = max(means, key=means.get); lo = min(means, key=means.get)
        gap = means[hi] - means[lo]
        separated = gap >= thr and min(groups[hi]) > max(groups[lo])
        return bool(separated), (f"{detail}  [too few samples for a test (best possible p "
                                 f"{_min_p(sizes):.2g}); extreme groups "
                                 f"{'completely separated' if separated else 'overlap'}, "
                                 f"gap {gap:+.2f}]"), None
    obs = _between(vals, labs)
    rng = random.Random(seed)
    perm, ge = list(labs), 0
    for _ in range(n_perm):
        rng.shuffle(perm)
        if _between(vals, perm) >= obs - 1e-9:
            ge += 1
    p = (ge + 1) / (n_perm + 1)
    return p < CONFOUND_P, (f"{detail}  [permutation F-test across {len(groups)} groups: "
                            f"p = {p:.2g}, {n_perm:,} relabellings; flagged at p < {CONFOUND_P}]"), p


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
                         "common-contaminant sequences, and the organism (taxid) the panels' "
                         "ortholog names are chosen for. Default: ./search.fasta.meta.json if present")
    ap.add_argument("--taxid", type=int,
                    help="NCBI taxid of the searched organism, when there is no --fasta-meta "
                         "(9606 human, 10090 mouse, 10116 rat have verified panel names)")
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
    tc, tc_note, meta = None, None, {}
    if fasta_meta:
        try:
            with open(fasta_meta) as fh:
                meta = json.load(fh)
            tc = target_contaminants(meta, a.keratin_sample)
        except (OSError, ValueError) as e:
            tc_note = f"not assessed: could not read {fasta_meta} ({e})"
    else:
        tc_note = "not assessed: no <fasta>.meta.json from fetch_fasta.py (pass --fasta-meta)"
    org = (tc or {}).get("organism") or "target-organism"
    taxid = a.taxid
    if taxid is None:
        try:
            taxid = int(meta.get("taxid")) if meta.get("taxid") else None
        except (TypeError, ValueError):
            taxid = None
    org_name = (tc or {}).get("organism") or PANEL_ORGANISMS.get(taxid) or (
        f"taxid {taxid}" if taxid else "an unknown organism")
    panels = {name: panel_genes(name, taxid) for name in PANELS}
    if tc and (tc["genes"] or tc["accessions"]):
        # Genes AND accessions: the matrix's id column may be either, and an NCBI database
        # carries no gene names at all.
        panels[TARGET_PANEL] = sorted(tc["genes"] | tc["accessions"])
    removed = removed_contaminant_rows(a.matrix)

    results, flags = {}, []
    for name, genes in panels.items():
        per, nhit, matched, cont_hits = panel_scores(rows, samples, genes, removed)
        z = zscore(per, samples)
        elevated = sorted(s for s in samples if z.get(s) is not None and z[s] >= a.z)
        confounded, detail, confound_p = confound_check(z, gmap, samples, a.z)
        expected = a.keratin_sample and name == "EPIDERMIS"
        kept = name == TARGET_PANEL
        if kept or expected:
            # A Cont_ row here is a real protein lost to a Cont_ entry (the lost-protein note
            # below says so), or the analyte of a keratin sample -- never "not sample".
            cont_hits = []
        results[name] = {"n_panel_proteins": nhit, "matched_genes": matched,
                         "contaminant_matches": cont_hits,
                         "z": {s: (round(z[s], 2) if z[s] is not None else None) for s in samples},
                         "elevated_samples": elevated, "group_confounded": confounded,
                         "confound_p": confound_p,
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
    # Did the check RUN? A panel that matched nothing is only evidence of absence when the
    # matrix's gene names are ones the panels know. If NO curated panel matched a single
    # sample protein -- each lists markers (CAT, PRDX2, TPM1, DSP ...) found in nearly any
    # cell or tissue proteome -- the names did not match, and "0 detected" means nothing.
    # (msalemi 2026-09-24: all three read "0 panel proteins -- not assessable" on a mouse
    # brain matrix with haemoglobin at ~15 log2.)
    any_match = any(results[n]["n_panel_proteins"] for n in PANELS)
    covered = taxid in PANEL_ORGANISMS
    for name in PANELS:
        r = results[name]
        if r["n_panel_proteins"]:
            r["status"], r["status_note"] = "assessed", None
        elif not any_match:
            r["status"] = "not_run"
            r["status_note"] = (
                f"check could not run (no panel gene matched this organism's symbols): none "
                f"of the {len(PANELS)} panels matched a single protein of {org_name}"
                + ("" if covered else
                   f" -- its gene names are not in the verified set ({', '.join(PANEL_ORGANISMS.values())})")
                + ". Check that the matrix has a Genes column, then pass --taxid / "
                  "--fasta-meta or add the organism's names to PANEL_ORTHOLOGS.")
        elif covered:
            r["status"] = "none_detected"
            r["status_note"] = (f"none of this panel's markers was detected; the check ran "
                                f"({org_name} names are covered, and other panels matched)")
        else:
            r["status"] = "none_detected_unverified"
            r["status_note"] = (f"none of this panel's markers matched; {org_name} gene names "
                                f"are not in the verified set, so this may be a naming mismatch "
                                f"rather than an absence")
    not_run = [n for n in PANELS if results[n]["status"] == "not_run"]
    if not_run:
        flags.append(f"Contamination panels {', '.join(not_run)}: "
                     + results[not_run[0]]["status_note"])
    # A database built before fetch_fasta.py removed target-identical contaminant entries can
    # hold the sample's OWN protein as a Cont_ entry; say so rather than call it a reagent.
    lost_accs = {(r.get("cont_acc") or "").upper() for r in (tc or {}).get("kept_as_contaminant", [])}
    for name in PANELS:
        cm = results[name]["contaminant_matches"]
        if not cm:
            continue
        msg = (f"{name}: {len(cm)} `{CONT_TAG}`-tagged marker(s) -- {CONTAMINANT_NOT_SAMPLE}: "
               + "; ".join(cm[:8]))
        own = [m for m in cm if any(acc and acc in m.upper() for acc in lost_accs)]
        if own:
            msg += (f". Except: {', '.join(own)} -- identical to {org} proteins in this database, "
                    f"so possibly the sample's own (see the database note).")
        elif (tc or {}).get("legacy_note"):
            msg += (". Caveat: this database predates the target-identical contaminant check, so a "
                    f"{CONT_TAG} marker may be the sample's own protein (see the database note).")
        flags.append(msg)

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
            if r.get("contaminant_matches"):
                fh.write(f"_`{CONT_TAG}`-tagged matches -- {CONTAMINANT_NOT_SAMPLE}:_ "
                         f"{'; '.join(r['contaminant_matches'])}\n\n")
            if not r["matched_genes"]:
                note = (r.get("status_note") if name != TARGET_PANEL else None) or \
                    "none of these proteins was quantified"
                fh.write(f"_{note[0].upper()}{note[1:]}._\n\n")
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
                   "taxid": taxid, "panel_organisms_verified": taxid in PANEL_ORGANISMS,
                   "panels": results, "detected_depth": depth, "flags": flags,
                   "fasta_meta": fasta_meta, "target_contaminants_note": tc_note,
                   "lost_to_contaminants": lost_msg}, jf, indent=2)

    print(json.dumps({"out": a.out, "flags": flags,
                      "group_confounded": [n for n, r in results.items() if r["group_confounded"]],
                      "not_run": not_run,
                      "matrix_complete": complete}, indent=2))


if __name__ == "__main__":
    main()
