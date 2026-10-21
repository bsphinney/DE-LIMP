#!/usr/bin/env python3
"""
contaminant_competition.py  --  Which of the sample's proteins did the contaminant list take away?

The Core searches every sample with its universal contaminant list (human and bovine entries,
Cont_-tagged). Several entries share identical or near-identical tryptic peptides with abundant
human and mouse proteins. DIA-NN can then report those peptides under the Cont_ entry, and
`--cont-quant-exclude Cont_` (and run_de.R's contaminant filter, the same rule) removes every
peptide that maps to a Cont_ entry from quantification and normalisation. Nothing errors and
nothing warns; the protein count just looks a little low. Measured incidents:
  * a mouse tissue search: mouse Actb and Eef1a1 had no protein group at all, and 28% of the ion
    current sat on Cont_ entries, most of it the mouse proteins' own;
  * a human IP: 68 human genes were lost the same way.
fetch_fasta.py already drops contaminant entries identical to a target, and those the search
cannot tell apart from one (its peptide rule). Which peptides a sample shows only the search can
say, so this check reads the report after every DIA-NN search. It reports; it never changes the
contaminant list or the search, and never fails a job (the Core director's decision, 2026-10-20:
keep the universal list for every search, and name what it takes).

WHAT IT READS. DIA-NN's report (report.parquet, or a .tsv report), rows kept as run_de.R keeps
them (diann_q_columns: the q-columns present, cutoff_for() each, a row dropped only when a value
exceeds its cutoff), and the searched FASTA (gene names; I/L twins), with its fetch_fasta.py
sidecar when there is one (organism, the digest the build used).

THE RULES. A precursor is a contaminant's when ANY accession in its Protein.Ids carries the tag:
DIA-NN's --cont-quant-exclude rule, which run_de.R applies too.
  * A contaminant entry SHOWS ITSELF when at least fetch_fasta.MIN_UNIQUE_PEPTIDES non-overlapping
    peptides of its own were found: peptides that name no sample protein in Protein.Ids and are
    not one of a sample protein's digest peptides, as written or with I = L
    (fetch_fasta.target_peptide_hits). The build keeps an entry by the same count in its digest.
  * A precursor shared between sample protein(s) and contaminant entries that do NOT show
    themselves is the sample's: the contaminant rule removed the sample's own signal.
  * Genes are grouped as protein inference groups them (parsimony: the fewest genes that explain
    the precursors), so a peptide shared by paralogs (actins, keratins) names one gene, not each.
  * COMPETED OUT: a sample gene with such precursors and no precursor of its own left: it has no
    quantity, so it is missing from the results.
  * SIPHONED: a sample gene that kept precursors of its own, but whose removed precursors carry
    more signal than those it kept: it is quantified from the minority of its signal.
  * A precursor shared with a contaminant entry that does show itself (bovine albumin in a
    cell-culture sample) cannot be attributed to either: such a gene is listed for staff as
    `masked`, never as lost.
  * The contaminant share of the precursor signal, per run, as contaminants.R's
    contaminant_share() computes it (Precursor.Quantity, else Precursor.Normalised).

Usage:
  contaminant_competition.py --report <out>/report.parquet --fasta <searched FASTA>
                             [--fasta-meta <fasta>.meta.json] [--q-cutoff 0.01]
                             [--out <out>/contaminant_competition.json]
                             [--provenance <out>/search_provenance.json]
  contaminant_competition.py --search-out <out>     # an existing search: the report and FASTA
                                                    # from its search_provenance.json; writes
                                                    # <out>/contaminant_competition.json there
Exit 0 whenever a record was made, findings or not (a check that could not run is a record that
says why); 2 for a usage error.
"""
import argparse
import csv
import json
import math
import os
import re
import shlex
import statistics
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# The tags, the evidence threshold, the digest and the target-peptide test are fetch_fasta.py's,
# never a copy here (DE-LIMP rule 3).
from fetch_fasta import (CONTAMINANT_TAGS, MIN_UNIQUE_PEPTIDES, REBUILD_ADVICE,  # noqa: E402
                         _contaminant_ids, _digestion_enzyme, _fasta_records, _header_ids,
                         _independent_count, fasta_accessions_and_targets,
                         is_contaminant_accession, sidecar_state, target_peptide_hits)
from diann_q_columns import DEFAULT_Q_CUTOFF, cutoff_for, fdr_columns  # noqa: E402

RECORD = "contaminant_competition.json"
PROVENANCE = "search_provenance.json"
TAG = "[contaminant_competition]"
# contaminants.R's CONTAMINANT_SHARE_COLUMNS, in its order (tests/test_contaminant_competition.py
# keeps the two equal): the measured signal first; DIA-NN's RT-dependent normalisation does not
# keep a run's signal fractions.
SHARE_COLUMNS = ("Precursor.Quantity", "Precursor.Normalised")
REPORT_COLUMNS = ("Run", "Precursor.Id", "Protein.Ids", "Stripped.Sequence")
# run_search.QUANTITY_META: an adapter's declaration (Sage, FragPipe, Radiant ...). Such a report
# is not DIA-NN's protein inference, so DIA-NN's question does not apply to it.
ADAPTED_META = b"delimp.quantity."
LIST_MAX = 12
RULE = ("a precursor is a contaminant's when any accession in its Protein.Ids carries the tag "
        "(DIA-NN's --cont-quant-exclude rule, applied by run_de.R too); a contaminant entry shows "
        f"itself with at least {MIN_UNIQUE_PEPTIDES} non-overlapping peptides that name no sample "
        "protein and are not a sample protein's digest peptide (I = L); a precursor shared only "
        "with entries that do not show themselves is the sample's")


def _num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _ids(cell):
    return [t.strip() for t in str(cell or "").split(";") if t.strip()]


def _is_cont(cell):
    return any(is_contaminant_accession(a) for a in _ids(cell))


def _new():
    return {"precursors": set(), "signal": 0.0}


# ---- reading the report ------------------------------------------------------------------------
def _tsv_rows(report, cuts):
    with open(report, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            try:
                if any((r.get(c) or "").strip() not in ("", "NA", "NaN")
                       and float(r[c]) > cut for c, cut in cuts.items()):
                    continue
            except ValueError:
                continue
            yield r


def _scan_tsv(report, cuts, qcol, ev, wanted):
    for r in _tsv_rows(report, cuts):
        if _is_cont(r.get("Protein.Ids")):
            _take_cont(ev, r.get("Precursor.Id"), r.get("Protein.Ids"), r.get("Stripped.Sequence"),
                       r.get(qcol))
    want = wanted(ev)
    for r in _tsv_rows(report, cuts):
        ids, q = r.get("Protein.Ids"), r.get(qcol)
        cont = _is_cont(ids)
        _take_run(ev, r.get("Run"), q, cont)
        if not cont and want:
            _take_kept(ev, want, ids, r.get("Precursor.Id"), q)


def _scan_parquet(report, cuts, qcol, ev, wanted):
    """Two scans, filtered by the q-columns: the contaminant rows (few) row by row; then every
    row's signal per run, and the kept rows of the wanted sample accessions -- aggregated by
    arrow, so a 400-run report is never held whole."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.dataset as pads
    d = pads.dataset(report, format="parquet")
    keep = None
    for c, cut in cuts.items():
        f = (pads.field(c) <= cut) | pads.field(c).is_null()
        keep = f if keep is None else keep & f
    tag_rx = "(^|;)(" + "|".join(re.escape(t) for t in CONTAMINANT_TAGS) + ")"
    # cast: a dictionary-encoded column has no regex kernel (DIA-NN writes plain strings)
    cont = pc.match_substring_regex(pads.field("Protein.Ids").cast(pa.string()), tag_rx)
    cols = ["Precursor.Id", "Protein.Ids", "Stripped.Sequence", qcol]
    for b in d.to_batches(columns=cols, filter=cont if keep is None else keep & cont):
        for prec, ids, pep, q in zip(*(b.column(c).to_pylist() for c in cols)):
            _take_cont(ev, prec, ids, pep, q)
    want = wanted(ev)
    want_rx = ("(^|;)(" + "|".join(re.escape(a) for a in sorted(want)) + ")(;|$)") if want else None
    for b in d.to_batches(columns=["Run", "Protein.Ids", "Precursor.Id", qcol], filter=keep):
        if not b.num_rows:
            continue
        t = pa.Table.from_batches([b])
        ids = pc.cast(t["Protein.Ids"], pa.string())
        is_cont = pc.fill_null(pc.match_substring_regex(ids, tag_rx), False)
        for sub, flag in ((t, False), (t.filter(is_cont), True)):
            if sub.num_rows:
                g = sub.group_by("Run").aggregate([(qcol, "sum")])
                for run, s in zip(g["Run"].to_pylist(), g[f"{qcol}_sum"].to_pylist()):
                    _take_run(ev, run, s, flag, summed=not flag)
        if want_rx:
            hit = pc.and_(pc.invert(is_cont), pc.fill_null(
                pc.match_substring_regex(ids, want_rx, ignore_case=True), False))
            sub = t.filter(hit)
            for ids_, prec, q in zip(*(sub.column(c).to_pylist()
                                       for c in ("Protein.Ids", "Precursor.Id", qcol))):
                _take_kept(ev, want, ids_, prec, q)


def _take_cont(ev, prec, ids, pep, q):
    rec = ev["cont"].setdefault((ids or "", (pep or "").upper()), _new())
    rec["precursors"].add(prec)
    rec["signal"] += _num(q) or 0.0


def _take_run(ev, run, q, is_cont, summed=True):
    """`summed`: a total over every row (the TSV path adds each row to the total itself)."""
    r = ev["runs"].setdefault(str(run), [0.0, 0.0])
    if summed or not is_cont:
        r[0] += _num(q) or 0.0
    if is_cont:
        r[1] += _num(q) or 0.0


def _take_kept(ev, want, ids, prec, q):
    if {x.upper() for x in _ids(ids)} & want:
        rec = ev["kept"].setdefault(ids, _new())
        rec["precursors"].add(prec)
        rec["signal"] += _num(q) or 0.0


def read_report(report, wanted, q_cutoff=DEFAULT_Q_CUTOFF, q_columns=None):
    """(evidence, None), ({"not_applicable": True}, None) or (None, why not). evidence: the
    q-columns and cutoffs applied, the quantity column, per-run [total, contaminant] signal, the
    contaminant precursors by (Protein.Ids, peptide), and by Protein.Ids the kept
    (contaminant-free) precursors of rows naming a sample accession wanted(evidence) returns
    (upper case) once the contaminant rows are read. A report that cannot be read is a why-not,
    never a traceback."""
    if not report or not os.path.isfile(report):
        return None, f"the report ({report or 'none given'}) is not readable here"
    name = os.path.basename(report)
    parquet = report.endswith(".parquet")
    try:
        if parquet:
            try:
                import pyarrow.parquet as pq
            except ImportError:
                return None, ("pyarrow is not installed in this Python, so the report could not "
                              "be read")
            schema = pq.read_schema(report)
            if any(k.startswith(ADAPTED_META) for k in (schema.metadata or {})):
                return {"not_applicable": True}, None
            names = schema.names
        else:
            with open(report, newline="") as fh:
                names = next(csv.reader(fh, delimiter="\t"), [])
        missing = [c for c in REPORT_COLUMNS if c not in names]
        if missing:
            return None, (f"{name} has no {', '.join(missing)} column (not a DIA-NN precursor "
                          "report)")
        qcol = next((c for c in SHARE_COLUMNS if c in names), None)
        if not qcol:
            return None, f"{name} has neither {' nor '.join(SHARE_COLUMNS)} (no precursor signal)"
        qcols = [c for c in (q_columns or fdr_columns(names)) if c in names]
        cuts = {c: cutoff_for(c, q_cutoff) for c in qcols}
        ev = {"report": report, "quantity_column": qcol, "q_cutoff": q_cutoff, "q_columns": qcols,
              "cutoffs": cuts, "runs": {}, "cont": {}, "kept": {}}
        (_scan_parquet if parquet else _scan_tsv)(report, cuts, qcol, ev, wanted)
    except (OSError, ValueError, csv.Error, AttributeError) as e:
        # ValueError: pyarrow's ArrowInvalid; AttributeError: a pyarrow older than group_by (7.0)
        why = " ".join(str(e).replace(report, name).split())
        return None, (f"{name} could not be read ({type(e).__name__}: {why[:200]}"
                      f"{'...' if len(why) > 200 else ''})")
    return ev, None


# ---- the searched FASTA and the sidecar --------------------------------------------------------
def read_fasta(fasta):
    """{"accs": {accession (upper): (entry, gene)}, "targets": [target records], "by_gene":
    {gene: {its target accessions (upper)}}, "cont_seq": {contaminant accession (upper):
    sequence}} of the searched FASTA."""
    with open(fasta, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    accs, targets = fasta_accessions_and_targets(text)
    by_gene = {}
    for header, _lines in targets:
        acc, _entry, gene = _header_ids(header)
        if gene:
            by_gene.setdefault(gene, set()).add(acc.upper())
    cont_seq = {}
    for header, lines in _fasta_records(text):
        c = _contaminant_ids(header)
        if c:
            cont_seq[c[1].upper()] = "".join(ln.strip() for ln in lines).upper()
    return {"accs": accs, "targets": targets, "by_gene": by_gene, "cont_seq": cont_seq}


def resolve_inputs(search_out=None, report=None, fasta=None, fasta_meta=None):
    """(report, fasta, sidecar) -- what was given, else what the search recorded in its
    search_provenance.json (`result.report`, `fasta`), and <fasta>.meta.json."""
    if search_out:
        try:
            with open(os.path.join(search_out, PROVENANCE)) as fh:
                prov = json.load(fh)
        except (OSError, ValueError):
            prov = {}
        res = prov.get("result") if isinstance(prov.get("result"), dict) else {}
        cands = [report, res.get("report"), os.path.join(search_out, "report.parquet")]
        report = next((p for p in cands if p and os.path.isfile(p)), report or cands[-1])
        fasta = fasta or prov.get("fasta")
    meta = None
    meta_path = fasta_meta or (fasta + ".meta.json" if fasta else None)
    if meta_path:
        try:
            with open(meta_path) as fh:
                meta = json.load(fh)
        except (OSError, ValueError):
            meta = None
    meta = meta if isinstance(meta, dict) else None
    return report, fasta or (meta or {}).get("fasta"), meta


# ---- the check ---------------------------------------------------------------------------------
def check(report, fasta=None, meta=None, q_cutoff=DEFAULT_Q_CUTOFF, q_columns=None,
          no_fasta_reason=None):
    """The record (a dict, always): status ok | not_run (with `reason`) | not_applicable.
    `no_fasta_reason`: why there is no FASTA, when the caller knows (its sha256 changed ...)."""
    rec = {"check": "contaminant_competition", "report": report, "fasta": fasta,
           "organism": (meta or {}).get("organism") or None, "rule": RULE,
           # which fetch_fasta.py rule built the database (legacy / identity_only / current)
           "database_state": sidecar_state(meta) if meta else None}
    from estimate_params import DIANN_DIGEST
    dig = (meta or {}).get("contaminant_digest")
    digest = dig if isinstance(dig, dict) and {"cut", "missed_cleavages", "min_pep_len",
                                               "max_pep_len"} <= set(dig) else DIANN_DIGEST
    fa, rec["fasta_note"] = None, None
    if fasta:
        try:
            fa = read_fasta(fasta)
        except OSError as e:
            rec["fasta_note"] = f"the searched FASTA could not be read ({e})"
    else:
        rec["fasta_note"] = no_fasta_reason or "no searched FASTA was given"
    if fa is None:
        rec["fasta_note"] += (": genes are named by accession, and a contaminant's peptide is "
                              "checked against Protein.Ids only (no I/L twins)")
    found = {}

    def wanted(ev):
        # the peptides that name no sample protein: is any of them a sample digest peptide?
        lone = {p for (ids, p) in ev["cont"] if not any(not is_contaminant_accession(a)
                                                         for a in _ids(ids))}
        found["on_target"] = {}
        if fa:
            try:
                found["on_target"] = target_peptide_hits(lone, fa["targets"], digest)
            except ValueError as e:              # a --cut this digest cannot read
                rec["fasta_note"] = (f"the search's digest could not be applied ({e}): a "
                                     "contaminant's peptide is checked against Protein.Ids only "
                                     "(no I/L twins)")
        accs = {a.upper() for ids, _p in ev["cont"] for a in _ids(ids)
                if not is_contaminant_accession(a)}
        accs |= {h["target_acc"].upper() for h in found["on_target"].values()}
        if fa:        # every accession of those genes: an isoform may keep the gene's peptides
            for g in {_gene(a, fa["accs"]) for a in accs}:
                accs |= fa["by_gene"].get(g, set())
        return accs

    ev, why = read_report(report, wanted, q_cutoff, q_columns)
    if ev is None:
        rec.update(status="not_run", reason=why)
        return rec
    if ev.get("not_applicable"):
        rec.update(status="not_applicable",
                   reason="the report is an adapted engine's (it declares its quantities), not "
                          "DIA-NN's protein inference")
        return rec
    rec.update({k: ev[k] for k in ("quantity_column", "q_cutoff", "q_columns", "cutoffs")})
    if not ev["runs"]:
        rec.update(status="not_run",
                   reason=f"no row of {os.path.basename(report)} passes the q-value filter "
                          f"(q-columns {', '.join(ev['q_columns']) or 'none'})")
        return rec
    rec.update(classify(ev, (fa or {}).get("accs") or {}, (fa or {}).get("cont_seq") or {},
                        found.get("on_target") or {}))
    rec["status"] = "ok"
    return rec


def _gene(acc, accs):
    return (accs.get(acc.upper()) or ("", ""))[1] or acc


def _lead_key(gene):
    # the lead of genes the precursors cannot tell apart: the shortest name, then A-Z
    # (EEF1A1 before its pseudogene EEF1A1P5, ACTB before ACTBL2)
    return (len(gene), gene)


def parsimony(feats, kind):
    """Which genes the precursors evidence, as protein inference groups them (the fewest genes
    that explain every precursor): repeatedly the gene that names the most precursors not yet
    explained takes all of them -- with every gene naming exactly the same ones, which the data
    cannot tell apart from it. A precursor shared by paralogs (actins, keratins) then counts for
    one gene, never for each paralog it also matches. -> [{"genes": [lead, indistinguishable ...],
    "feats": [...]}] for each group that holds a feature of `kind`."""
    left, groups = list(feats), []
    while any(f["kind"] == kind for f in left):
        by_gene = {}
        for i, f in enumerate(left):
            for g in f["genes"]:
                by_gene.setdefault(g, set()).add(i)
        n = {g: len(set().union(*(left[i]["precursors"] for i in idx)))
             for g, idx in by_gene.items()}
        best = max(n.values())
        lead = min((g for g in n if n[g] == best), key=_lead_key)
        same = sorted((g for g in by_gene if by_gene[g] == by_gene[lead] and g != lead),
                      key=_lead_key)
        mine = [left[i] for i in sorted(by_gene[lead])]
        left = [f for i, f in enumerate(left) if i not in by_gene[lead]]
        if any(f["kind"] == kind for f in mine):
            groups.append({"genes": [lead] + same, "feats": mine})
    return groups


def classify(ev, accs, cont_seq, on_target):
    """The findings from the evidence: the share, which contaminant entries show themselves, and
    the sample genes COMPETED OUT / SIPHONED / masked (and how many only lost a minority)."""
    # 1. each contaminant precursor's entries and sample proteins (an I/L twin of a sample digest
    # peptide names that protein), and the entries' own peptides
    own, rows = {}, []
    for (ids, pep), ag in ev["cont"].items():
        conts = [a for a in _ids(ids) if is_contaminant_accession(a)]
        targets = [a for a in _ids(ids) if not is_contaminant_accession(a)]
        twin = on_target.get(pep) if not targets else None
        if twin:
            targets = [twin["target_acc"]]
        if not targets:
            for c in conts:
                own.setdefault(c, set()).add(pep)
        rows.append((conts, targets, pep, ag, twin))

    def n_own(c):
        """Non-overlapping own peptides (fetch_fasta._independent_count, as the build counts them):
        missed-cleavage variants of one stretch count once. Distinct peptides without a FASTA."""
        seq, peps = cont_seq.get(c.upper()), sorted(own.get(c, ()))
        if not seq:
            return len(peps)
        spans = []
        for i, p in enumerate(peps):
            at = seq.find(p)
            # not in the entry as written (should not happen): its own stretch, past the end
            spans.append((at, at + len(p)) if at >= 0 else (len(seq) + 2 * i, len(seq) + 2 * i + 1))
        return _independent_count(spans)
    seen = {c for conts, *_ in rows for c in conts}
    shows = {c for c in seen if n_own(c) >= MIN_UNIQUE_PEPTIDES}

    # 2. the features: the sample's precursors the rule removed ("lost"), those shared with an
    # entry that shows itself ("masked"), and the kept (contaminant-free) ones -- each with the
    # genes it names
    feats = []
    for conts, targets, pep, ag, twin in rows:
        if targets:
            feats.append({"kind": "masked" if any(c in shows for c in conts) else "lost",
                          "genes": {_gene(t, accs) for t in targets}, "accs": set(targets),
                          "conts": set(conts), "twin": pep if twin else None, **ag})
    for ids, ag in ev["kept"].items():
        targets = [a for a in _ids(ids) if not is_contaminant_accession(a)]
        feats.append({"kind": "kept", "genes": {_gene(t, accs) for t in targets},
                      "accs": set(targets), "conts": set(), "twin": None, **ag})

    def entry(grp):
        lost = [f for f in grp["feats"] if f["kind"] == "lost"]
        kept = [f for f in grp["feats"] if f["kind"] == "kept"]
        ls, ks = sum(f["signal"] for f in lost), sum(f["signal"] for f in kept)
        mine = set(grp["genes"])
        return {"gene": grp["genes"][0], "indistinguishable": grp["genes"][1:],
                "accessions": sorted({a for f in grp["feats"] for a in f["accs"]
                                      if _gene(a, accs) in mine}),
                "contaminants": sorted({c for f in lost for c in f["conts"]}),
                "n_removed_precursors": len(set().union(*(f["precursors"] for f in lost))),
                "n_kept_precursors": len(set().union(*(f["precursors"] for f in kept))),
                "removed_signal_pct": round(100 * ls / (ls + ks), 1) if ls + ks else None,
                "removed_signal": ls, "kept_signal": ks,
                **({"il_twin_peptides": sorted({f["twin"] for f in lost if f["twin"]})}
                   if any(f["twin"] for f in lost) else {})}
    competed, siphoned, minority = [], [], 0
    for grp in parsimony([f for f in feats if f["kind"] in ("lost", "kept")], "lost"):
        e = entry(grp)
        if e["n_removed_precursors"] and not e["n_kept_precursors"]:
            competed.append(e)
        elif e["removed_signal"] > e["kept_signal"]:
            siphoned.append(e)
        else:
            minority += 1
    # genes seen ONLY through precursors shared with an entry that shows itself
    named = {g for f in feats if f["kind"] != "masked" for g in f["genes"]}
    masked = [{"gene": grp["genes"][0], "indistinguishable": grp["genes"][1:],
               "contaminants": sorted({c for f in grp["feats"] for c in f["conts"]}),
               "n_shared_precursors": len(set().union(*(f["precursors"] for f in grp["feats"])))}
              for grp in parsimony([f for f in feats if f["kind"] == "masked"
                                    and not f["genes"] & named], "masked")]
    competed.sort(key=lambda x: (-x["removed_signal"], x["gene"]))
    siphoned.sort(key=lambda x: (-x["removed_signal"], x["gene"]))
    masked.sort(key=lambda x: (-x["n_shared_precursors"], x["gene"]))

    # what --drop-contaminant takes for a re-search: the entries behind the lost genes -- never a
    # digestion enzyme in use, which fetch_fasta.py keeps on purpose
    behind = sorted({c for x in competed + siphoned for c in x["contaminants"]})
    share = {run: round(100 * c / t, 2) for run, (t, c) in sorted(ev["runs"].items()) if t > 0}
    pcts = list(share.values())
    return {"n_runs": len(ev["runs"]),
            "contaminant_share_pct": {"per_run": share,
                                      "median": round(statistics.median(pcts), 2) if pcts else None,
                                      "max": max(pcts) if pcts else None,
                                      "column": ev["quantity_column"]},
            "n_contaminant_entries_seen": len(seen),
            "contaminants_showing_themselves": sorted(shows),
            "contaminants_not_showing_themselves": sorted(seen - shows),
            "competed_out": competed, "siphoned": siphoned, "masked": masked,
            "n_genes_lost_minority": minority,
            "drop_contaminant": [c for c in behind if not _digestion_enzyme(c)],
            "kept_digestion_enzymes": [c for c in behind if _digestion_enzyme(c)]}


# ---- what people read --------------------------------------------------------------------------
def _names(lst, n=LIST_MAX):
    names = [x["gene"] for x in lst]
    return ", ".join(names[:n]) + (f" and {len(names) - n} more" if len(names) > n else "")


def _n(n, one, many):
    return f"{n} {one if n == 1 else many}"


def _same(x):
    same = x.get("indistinguishable") or []
    return f"= {', '.join(same[:3])}{' ...' if len(same) > 3 else ''}; " if same else ""


def _via(x):
    return ", ".join(x["contaminants"][:4]) + (" ..." if len(x["contaminants"]) > 4 else "")


def competed_text(rec):
    co, org = rec.get("competed_out") or [], rec.get("organism") or "sample"
    if not co:
        return ""
    return (f"COMPETED OUT: {_n(len(co), org + ' gene has', org + ' genes have')} no quantity "
            "left: "
            "every precursor of theirs that passes the q-value filter is shared with a contaminant "
            "entry that shows no peptides of its own here, so DIA-NN's --cont-quant-exclude and "
            "run_de.R's contaminant filter removed them all, and they are missing from the "
            "results: "
            + "; ".join(f"{x['gene']} ({_same(x)}"
                        f"{_n(x['n_removed_precursors'], 'precursor', 'precursors')} via {_via(x)})"
                        for x in co[:LIST_MAX])
            + (f"; +{len(co) - LIST_MAX} more" if len(co) > LIST_MAX else "") + ".")


def siphoned_text(rec):
    si, org = rec.get("siphoned") or [], rec.get("organism") or "sample"
    if not si:
        return ""
    return (f"SIPHONED: {_n(len(si), org + ' gene', org + ' genes')} kept precursors of their own, "
            "but most of their signal was on precursors removed as contaminant (shared with an "
            "entry that shows no peptides of its own here): "
            + "; ".join(f"{x['gene']} ({_same(x)}{x['removed_signal_pct']:g}% of its signal "
                        f"removed; {x['n_removed_precursors']} removed, "
                        f"{x['n_kept_precursors']} kept; via {_via(x)})" for x in si[:LIST_MAX])
            + (f"; +{len(si) - LIST_MAX} more" if len(si) > LIST_MAX else "") + ".")


def fix_text(rec):
    drop = rec.get("drop_contaminant") or []
    if not drop:
        return ""
    old = (f"The database was built by an older fetch_fasta.py rule ({rec['database_state']}), "
           "which kept contaminant entries that are sample proteins: " + REBUILD_ADVICE +
           " That alone removes most of them; check again on the new search. "
           if rec.get("database_state") in ("legacy", "identity_only") else "")
    return (old + "To recover them: rebuild the FASTA with the same fetch_fasta.py command plus "
            + " ".join(f"--drop-contaminant {c}" for c in drop)
            + ", and search again (every other contaminant entry stays)."
            + (f" Not dropped: {', '.join(rec['kept_digestion_enzymes'])} (the digestion "
               "enzyme)." if rec.get("kept_digestion_enzymes") else ""))


def notes_text(rec):
    ma, out = rec.get("masked") or [], []
    if rec.get("n_genes_lost_minority"):
        out.append(f"{_n(rec['n_genes_lost_minority'], 'gene', 'genes')} lost a minority of "
                   "their signal to shared precursors (in the record).")
    if ma:
        out.append(f"Not counted as lost: {_n(len(ma), 'gene is', 'genes are')} seen only "
                   "through precursors shared with a contaminant entry that does show peptides "
                   "of its own, so the data cannot say whether "
                   f"{'it is' if len(ma) == 1 else 'they are'} in the sample: {_names(ma)}.")
    if rec.get("fasta_note"):
        out.append(f"({rec['fasta_note']}.)")
    s = rec.get("contaminant_share_pct") or {}
    if s.get("median") is not None:
        out.append(f"Contaminant entries carry a median {s['median']:g}% (max {s['max']:g}%) of "
                   f"a run's precursor signal ({s.get('column')}, {rec.get('n_runs')} runs).")
    return " ".join(out)


def staff_message(rec):
    """The technical record: AUDIT.md, the job log, the watcher. Check codes and accessions
    allowed; never the client's words."""
    st = rec.get("status")
    if st == "not_run":
        return f"Not run: {rec.get('reason')}."
    if st == "not_applicable":
        return f"Not applicable: {rec.get('reason')}."
    head = " ".join(x for x in (competed_text(rec), siphoned_text(rec)) if x)
    if not head:
        org = rec.get("organism") or "sample"
        quiet = len(rec.get("contaminants_not_showing_themselves") or [])
        head = (f"Checked: no {org} gene lost most of its signal, or all of it, to a contaminant "
                f"entry ({_n(quiet, 'contaminant entry', 'contaminant entries')} seen without "
                "peptides of their own).")
    return " ".join(x for x in (head, fix_text(rec), notes_text(rec)) if x)


def plain_competed(rec):
    co = rec.get("competed_out") or []
    one = len(co) == 1
    it, its = ("it", "its") if one else ("them", "their")
    return (f"{_names(co)} {'is' if one else 'are'} missing from the results because the "
            f"contaminant list absorbed {its} peptides: every peptide of {it} that the search "
            "found is shared with an entry of the contaminant list, and those peptides are left "
            "out of the measurements. The Core can search again without those contaminant "
            f"entries to recover {it}.")


def plain_siphoned(rec):
    si = rec.get("siphoned") or []
    one = len(si) == 1
    return (f"Most of the signal of {_names(si)} was on peptides "
            f"{'it shares' if one else 'they share'} with an entry of the contaminant list, which "
            f"are left out of the measurements, so {'it is' if one else 'they are'} measured from "
            f"{'its' if one else 'their'} remaining peptides only. The Core can search again "
            "without those contaminant entries.")


def audit_findings(rec):
    """[(status, message, detail, plain)] for audit_results.py; [] for a record that does not
    apply (another engine's report). COMPETED OUT and SIPHONED are each a WARN, never a FAIL: the
    check reports and never stops an analysis. A check that could not run is an INFO saying why.
    plain None = the check's usual words (audit_results.PLAIN)."""
    if not isinstance(rec, dict) or rec.get("status") == "not_applicable":
        return []
    if rec.get("status") != "ok":
        return [("INFO", staff_message(rec), {"reason": rec.get("reason"),
                                              "report": rec.get("report")}, None)]
    keys = ("report", "fasta", "organism", "contaminant_share_pct", "drop_contaminant",
            "contaminants_not_showing_themselves", "q_columns", "cutoffs", "record_file")
    detail = {k: rec.get(k) for k in keys if rec.get(k) is not None}
    detail.update({k: [x["gene"] for x in rec.get(k) or []]
                   for k in ("competed_out", "siphoned", "masked")})
    tail = " ".join(x for x in (fix_text(rec), notes_text(rec)) if x)
    out = []
    if rec.get("competed_out"):
        out.append(("WARN", f"{competed_text(rec)} {tail}", detail, plain_competed(rec)))
    if rec.get("siphoned"):
        out.append(("WARN", f"{siphoned_text(rec)} {fix_text(rec) if out else tail}".strip(),
                    detail, plain_siphoned(rec)))
    return out or [("PASS", staff_message(rec), detail, None)]


def say(rec, stream=None):
    """The job log's line: WARNING when a gene was lost, else INFO (NOTE when it could not run)."""
    level = ("WARNING" if rec.get("competed_out") or rec.get("siphoned") else
             "NOTE" if rec.get("status") == "not_run" else "INFO")
    print(f"{TAG} {level}: {staff_message(rec)}", file=stream or sys.stderr, flush=True)


# ---- the record --------------------------------------------------------------------------------
def write(rec, out=None, provenance=None):
    """<out> (JSON) and search_provenance.json `contaminant_competition` (merged in place: a .tmp
    then a rename, so a reader never sees half a file). -> the files written."""
    written = []
    rec["message"] = staff_message(rec)            # what the watcher and a person read
    if out:
        rec["record_file"] = os.path.abspath(out)
        with open(out + ".tmp", "w") as fh:
            json.dump(rec, fh, indent=2)
        os.replace(out + ".tmp", out)
        written.append(out)
    if provenance and os.path.isfile(provenance):
        with open(provenance) as fh:
            prov = json.load(fh)
        prov["contaminant_competition"] = rec
        with open(provenance + ".tmp", "w") as fh:
            json.dump(prov, fh, indent=2)
        os.replace(provenance + ".tmp", provenance)
        written.append(provenance)
    return written


def job_lines(report, fasta, out_dir):
    """The bash a DIA-NN search job ends with (diann_parallel.py step 5, run_search.py's
    single-shot job): this check on its report, under the generator's own python (the pipeline
    env's, which has pyarrow), else python3. It never fails the job; a check that cannot start is
    said in the log, and AUDIT.md then checks the report itself."""
    report, fasta, out_dir = (os.path.abspath(p) for p in (report, fasta, out_dir))
    return [f'PY={shlex.quote(sys.executable)}; [ -x "$PY" ] || PY=python3; '
            f'"$PY" {shlex.quote(os.path.abspath(__file__))} --report {shlex.quote(report)} '
            f'--fasta {shlex.quote(fasta)} --out {shlex.quote(os.path.join(out_dir, RECORD))} '
            f'--provenance {shlex.quote(os.path.join(out_dir, PROVENANCE))} '
            '|| echo "WARNING: the contaminant competition check could not run -- the messages '
            'above say why" >&2']


def load(search_out):
    """The record a search wrote in its output folder, or None."""
    try:
        with open(os.path.join(search_out or "", RECORD)) as fh:
            rec = json.load(fh)
    except (OSError, ValueError):
        return None
    return rec if isinstance(rec, dict) else None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", help="DIA-NN's report (report.parquet, or a .tsv report)")
    ap.add_argument("--fasta", help="the FASTA the search used")
    ap.add_argument("--fasta-meta", help="its fetch_fasta.py sidecar (default <fasta>.meta.json)")
    ap.add_argument("--search-out", help="a finished search's output folder: the report and FASTA "
                                         "from its search_provenance.json, the record written "
                                         "there")
    ap.add_argument("--q-cutoff", type=float, default=DEFAULT_Q_CUTOFF,
                    help="run_de.R's --q-cutoff (default %(default)s)")
    ap.add_argument("--out", help="write the record here (JSON)")
    ap.add_argument("--provenance", help="search_provenance.json to record it in")
    a = ap.parse_args(argv)
    if not (a.report or a.search_out):
        ap.error("give --report or --search-out")
    report, fasta, meta = resolve_inputs(a.search_out, a.report, a.fasta, a.fasta_meta)
    rec = check(report, fasta, meta, a.q_cutoff)
    out = a.out or (os.path.join(a.search_out, RECORD) if a.search_out else None)
    prov = a.provenance or (os.path.join(a.search_out, PROVENANCE) if a.search_out else None)
    try:
        written = write(rec, out, prov)
    except (OSError, ValueError) as e:
        written = []
        print(f"{TAG} NOTE: the record could not be written ({e}); it is in this log only.",
              file=sys.stderr)
    say(rec)
    if written:
        print(f"{TAG} recorded in {', '.join(written)}", file=sys.stderr)
    summary = {k: rec.get(k) for k in ("status", "reason", "drop_contaminant")}
    summary.update({k: [x["gene"] for x in rec.get(k) or []]
                    for k in ("competed_out", "siphoned", "masked")})
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
