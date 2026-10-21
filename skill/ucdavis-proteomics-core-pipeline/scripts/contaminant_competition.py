#!/usr/bin/env python3
"""
contaminant_competition.py  --  Which of the sample's proteins did the contaminant list take away?

The Core searches every sample with its universal contaminant list (human and bovine entries,
Cont_-tagged). Several entries share identical or near-identical tryptic peptides with abundant
human and mouse proteins. DIA-NN can then report those peptides under the Cont_ entry, and the
contaminant rule -- DIA-NN's `--cont-quant-exclude Cont_`, and a DE that removes contaminant
precursors -- takes every peptide that maps to a Cont_ entry out of quantification and
normalisation. Nothing errors and nothing warns; the protein count just looks a little low.
Measured incidents:
  * a mouse tissue search: mouse Actb and Eef1a1 had no protein group at all, and 28% of the ion
    current sat on Cont_ entries, most of it the mouse proteins' own;
  * a human IP: 68 human genes were lost the same way.
fetch_fasta.py already drops contaminant entries identical to a target, and those the search
cannot tell apart from one (its peptide rule). Which peptides a sample shows only the search can
say, so this check reads the report after every DIA-NN search. It reports; it never changes the
contaminant list or the search, and never fails a job or an audit (the Core director's decision,
2026-10-20: keep the universal list for every search, and name what it takes).

WHAT IT READS. DIA-NN's report (report.parquet, or a .tsv report), rows kept as run_de.R keeps
them (diann_q_columns: the q-columns present, cutoff_for() each, keep_expression / row_passes),
and the searched FASTA (gene names; I/L twins), with its fetch_fasta.py sidecar when there is one
(organism, the digest the build used, which rule built the database). What took the shared
precursors is said from the records, never assumed (DE-LIMP rule 1): the DE's own contaminant
record when an audit has one, else DIA-NN's --cont-quant-exclude as make_methods.py reads it.

THE RULES. A precursor is a contaminant's when ANY accession in its Protein.Ids carries the tag:
DIA-NN's --cont-quant-exclude rule, which run_de.R applies too.
  * A contaminant entry SHOWS ITSELF when at least fetch_fasta.MIN_UNIQUE_PEPTIDES non-overlapping
    peptides of its own were found: peptides that name no sample protein in Protein.Ids and are
    not one of a sample protein's digest peptides, as written or with I = L
    (fetch_fasta.target_peptide_hits). The build keeps an entry by the same count in its digest.
  * A precursor shared between sample protein(s) and contaminant entries that do NOT show
    themselves is the sample's: the contaminant rule takes the sample's own signal.
  * The kept (contaminant-free) precursors are read for every protein connected to those through
    shared precursors, and the genes are grouped as protein inference groups them (parsimony: the
    fewest genes that explain the precursors). A peptide shared by paralogs (actins, keratins)
    then names one gene, not each, and a paralog's own peptides take the ones it shares.
  * COMPETED OUT: a sample gene with such precursors and no precursor of its own left (one it
    shares with another sample protein counts for that protein): under the contaminant rule it
    has no quantity.
  * SIPHONED: a sample gene that kept precursors of its own, but whose taken precursors carry
    more signal than those it kept: it is quantified from the minority of its signal.
  * A precursor shared with a contaminant entry that does show itself (bovine albumin in a
    cell-culture sample) cannot be attributed to either: a gene seen only, or mostly, through
    such precursors is listed for staff as `masked`, never as lost.
  * A keratin sample's keratin-family entries keep their precursors, as the DE keeps them
    (contaminants.R is_contaminant(exempt=), mirrored by exempted()).
  * Not applied: run_de.R's optional QuantUMS cutoffs and a run subset.
  * The contaminant share of the precursor signal, per run, as contaminants.R's
    contaminant_share() computes it (Precursor.Quantity, else Precursor.Normalised).
The record carries the report's size and modification time (`report_stat`): a record is read
only for the report it describes (describes()), never for a later search in the same folder.

Usage:
  contaminant_competition.py --report <out>/report.parquet --fasta <searched FASTA>
                             [--fasta-meta <fasta>.meta.json] [--q-cutoff 0.01]
                             [--out <out>/contaminant_competition.json]
                             [--provenance <out>/search_provenance.json]
  contaminant_competition.py --search-out <out>     # an existing search (in a job on a cluster:
                                                    # it reads report.parquet): the report and
                                                    # FASTA from its search_provenance.json; writes
                                                    # <out>/contaminant_competition.json there
Exit 0 whenever a record was made, findings or not (a check that could not run is a record that
says why); 2 for a usage error.
"""
import argparse
import csv
import heapq
import json
import math
import os
import re
import shlex
import signal
import statistics
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

# The tags, the evidence threshold, the digest and the target-peptide test are fetch_fasta.py's,
# the q-value filter diann_q_columns.py's -- never a copy here (DE-LIMP rule 3).
from fetch_fasta import (CONTAMINANT_TAGS, MIN_UNIQUE_PEPTIDES, REBUILD_ADVICE,  # noqa: E402
                         _contaminant_ids, _digestion_enzyme, _fasta_records, _header_ids,
                         _independent_count, fasta_accessions_and_targets,
                         is_contaminant_accession, keratin_database_check,
                         keratin_sample_recorded, sidecar_state, target_peptide_hits)
from diann_q_columns import (DEFAULT_Q_CUTOFF, cutoff_for, fdr_columns,  # noqa: E402
                             keep_expression, row_passes)

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
# The most a search job spends on this check (seconds): 710k rows take 4 s on HIVE. Inside a
# SLURM job it also stops two minutes before the job's own time limit (job_budget): a search job
# that ran out of time over this report would read as a failed search, and its fix -- more time,
# resubmit -- would redo a search that was complete (review of 2.14). Stopped, it records why,
# and the audit checks the report itself.
JOB_TIME_LIMIT = 900
# What the job still needs after the check: the job-end hook's own limit, and a minute
JOB_MARGIN_EXTRA_S = 60
MIN_BUDGET_S = 30


class OutOfTime(BaseException):
    """The job budget ran out. A BaseException, so no `except OSError` / `except Exception` inside
    the check takes it for an ordinary error and carries on with a degraded "ok" record (review
    of 2.14: TimeoutError is an OSError)."""
RULE = ("a precursor is a contaminant's when any accession in its Protein.Ids carries the tag "
        "(DIA-NN's --cont-quant-exclude rule, applied by run_de.R too); a contaminant entry shows "
        f"itself with at least {MIN_UNIQUE_PEPTIDES} non-overlapping peptides that name no sample "
        "protein and are not a sample protein's digest peptide (I = L); a precursor shared only "
        "with entries that do not show themselves is the sample's; genes grouped by parsimony")


def _num(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _ids(cell):
    return [t.strip() for t in str(cell or "").split(";") if t.strip()]


def _targets(cell):
    return [a for a in _ids(cell) if not is_contaminant_accession(a)]


def _is_cont(cell):
    return any(is_contaminant_accession(a) for a in _ids(cell))


def exempted(cell, exempt):
    """True when a keratin sample's exemption keeps this precursor -- contaminants.R's
    is_contaminant(exempt=), mirrored: every tagged accession is one of the keratin-family entries
    (`shared`), and it also names a sample protein or every tagged one is the searched species'
    own (`alone`)."""
    if not exempt or not exempt.get("shared"):
        return False
    ids = _ids(cell)
    tagged = [a for a in ids if is_contaminant_accession(a)]
    return bool(tagged) and all(a in exempt["shared"] for a in tagged) and (
        len(tagged) < len(ids) or all(a in exempt.get("alone", ()) for a in tagged))


def _new():
    return {"precursors": set(), "signal": 0.0}


def _add(store, key, prec, q, run=None):
    rec = store.setdefault(key, _new())
    rec["precursors"].add(prec)
    rec["signal"] += _num(q) or 0.0
    if run is not None:                      # a contaminant row's signal per run (the share)
        by = rec.setdefault("by_run", {})
        by[run] = by.get(run, 0.0) + (_num(q) or 0.0)


# ---- reading the report ------------------------------------------------------------------------
def _tsv_rows(report, cuts):
    with open(report, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if row_passes(r, cuts):
                yield r


def _scan_tsv(report, cuts, qcol, ev):
    for r in _tsv_rows(report, cuts):
        ids, q = r.get("Protein.Ids") or "", r.get(qcol)
        run = ev["runs"].setdefault(str(r.get("Run")), [0.0, 0.0])
        run[0] += _num(q) or 0.0
        if _is_cont(ids):
            run[1] += _num(q) or 0.0
            _add(ev["cont"], (ids, (r.get("Stripped.Sequence") or "").upper()),
                 r.get("Precursor.Id"), q, str(r.get("Run")))
        elif ids:
            _add(ev["kept"], ids, r.get("Precursor.Id"), q)


def _scan_parquet(report, cuts, qcol, ev):
    """Two scans of the rows run_de.R keeps: the contaminant rows (few) row by row; then every
    row's signal per run, and the kept rows' distinct precursors and signal per Protein.Ids --
    aggregated by arrow batch by batch, so a 400-run report is never held whole."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.dataset as pads
    d = pads.dataset(report, format="parquet")
    keep = keep_expression(cuts)
    tag_rx = "(^|;)(" + "|".join(re.escape(t) for t in CONTAMINANT_TAGS) + ")"
    # cast: a dictionary-encoded column has no regex kernel (DIA-NN writes plain strings)
    cont = pc.match_substring_regex(pads.field("Protein.Ids").cast(pa.string()), tag_rx)
    cols = ["Precursor.Id", "Protein.Ids", "Stripped.Sequence", qcol, "Run"]
    for b in d.to_batches(columns=cols, filter=cont if keep is None else keep & cont):
        for prec, ids, pep, q, run in zip(*(b.column(c).to_pylist() for c in cols)):
            _add(ev["cont"], (ids or "", (pep or "").upper()), prec, q, str(run))
    for b in d.to_batches(columns=["Run", "Protein.Ids", "Precursor.Id", qcol], filter=keep):
        if not b.num_rows:
            continue
        t = pa.Table.from_batches([b])
        # a NaN or infinite quantity counts as nothing, as in the text reader: arrow's sum would
        # carry a NaN into the run's whole total (review of 2.14)
        q = pc.cast(t[qcol], pa.float64())
        t = t.set_column(t.schema.get_field_index(qcol), qcol,
                         pc.if_else(pc.is_finite(q), q, pa.scalar(None, pa.float64())))
        is_cont = pc.fill_null(pc.match_substring_regex(
            pc.cast(t["Protein.Ids"], pa.string()), tag_rx), False)
        for sub, slot in ((t, 0), (t.filter(is_cont), 1)):
            if sub.num_rows:
                g = sub.group_by("Run").aggregate([(qcol, "sum")])
                for run, s in zip(g["Run"].to_pylist(), g[f"{qcol}_sum"].to_pylist()):
                    ev["runs"].setdefault(str(run), [0.0, 0.0])[slot] += _num(s) or 0.0
        sub = t.filter(pc.invert(is_cont))
        if sub.num_rows:
            # one entry per protein set and batch, not per row (a 1,000-run report)
            g = sub.group_by("Protein.Ids").aggregate([("Precursor.Id", "distinct"),
                                                       (qcol, "sum")])
            for ids, precs, s in zip(g["Protein.Ids"].to_pylist(),
                                     g["Precursor.Id_distinct"].to_pylist(),
                                     g[f"{qcol}_sum"].to_pylist()):
                if ids:
                    rec = ev["kept"].setdefault(ids, _new())
                    rec["precursors"].update(precs)
                    rec["signal"] += _num(s) or 0.0


def read_report(report, q_cutoff=DEFAULT_Q_CUTOFF, q_columns=None):
    """(evidence, None), ({"not_applicable": True}, None) or (None, why not). evidence: the
    q-columns and cutoffs applied, the quantity column, per-run [total, contaminant] signal, the
    contaminant precursors by (Protein.Ids, peptide), and the kept (contaminant-free) precursors
    by Protein.Ids. A report that cannot be read is a why-not."""
    if not report or not os.path.isfile(report):
        return None, f"the report ({report or 'none given'}) is not readable here"
    name = os.path.basename(report)
    parquet = report.endswith(".parquet")
    if parquet:
        try:
            import pyarrow.parquet as pq
        except ImportError:
            return None, "pyarrow is not installed in this Python, so the report could not be read"
        schema = pq.read_schema(report)
        if any(k.startswith(ADAPTED_META) for k in (schema.metadata or {})):
            return {"not_applicable": True}, None
        names = schema.names
    else:
        with open(report, newline="") as fh:
            names = next(csv.reader(fh, delimiter="\t"), [])
    missing = [c for c in REPORT_COLUMNS if c not in names]
    if missing:
        return None, f"{name} has no {', '.join(missing)} column (not a DIA-NN precursor report)"
    qcol = next((c for c in SHARE_COLUMNS if c in names), None)
    if not qcol:
        return None, f"{name} has neither {' nor '.join(SHARE_COLUMNS)} (no precursor signal)"
    qcols = [c for c in (q_columns or fdr_columns(names)) if c in names]
    cuts = {c: cutoff_for(c, q_cutoff) for c in qcols}
    ev = {"report": report, "quantity_column": qcol, "q_cutoff": q_cutoff, "q_columns": qcols,
          "cutoffs": cuts, "runs": {}, "cont": {}, "kept": {}}
    (_scan_parquet if parquet else _scan_tsv)(report, cuts, qcol, ev)
    return ev, None


def report_stat(report):
    """{size, mtime} of a report as it is now, or None: what describes() compares."""
    try:
        st = os.stat(report)
    except (OSError, TypeError):
        return None
    return {"size": st.st_size, "mtime": int(st.st_mtime)}


def describes(rec, report):
    """True when `rec` is this check's record of the report at `report` as it is now: the same
    file name, size and modification time. A re-search into the same folder writes a new report,
    so an earlier search's record is never read for it -- nor a record without `report_stat`."""
    if not isinstance(rec, dict) or not report:
        return False
    st = rec.get("report_stat")
    return (isinstance(st, dict) and st == report_stat(report)
            and os.path.basename(str(rec.get("report") or "")) == os.path.basename(str(report)))


def agrees_with_de(rec, de_prov):
    """True when the record's q-value filter is the DE's (de_provenance.json q_cutoff /
    q_columns, where recorded) and so is its keratin exemption: a search-time record uses
    run_de.R's defaults. (Not compared: run_de.R's optional QuantUMS cutoffs and a run subset,
    which this check does not apply.)"""
    q, cols = (de_prov or {}).get("q_cutoff"), (de_prov or {}).get("q_columns")
    if isinstance(q, (int, float)) and rec.get("q_cutoff") != float(q):
        return False
    if isinstance(cols, list) and cols and set(cols) != set(rec.get("q_columns") or ()):
        return False
    ex = de_exemption(de_prov)
    mine = rec.get("keratin_exempt") or {}
    return ex is None or (ex["shared"] == set(mine.get("shared") or ())
                          and ex["alone"] == set(mine.get("alone") or ()))


def de_exemption(de_prov):
    """The keratin exemption the DE applied (de_provenance.json contaminants.keratin_sample:
    exempt_accessions / exempt_alone_accessions) -- {"shared", "alone"}, empty when it applied
    none -- or None when the DE record does not say."""
    cont = (de_prov or {}).get("contaminants")
    k = cont.get("keratin_sample") if isinstance(cont, dict) else None
    if not isinstance(k, dict):
        return None
    if k.get("under_quantified") is True:
        # the DE kept those precursors but its quantities are DIA-NN's own, which left them out
        # (contaminants.R diann_kept_quant): nothing is exempt from what this check reports
        return {"shared": set(), "alone": set()}
    return {"shared": set(k.get("exempt_accessions") or ()),
            "alone": set(k.get("exempt_alone_accessions") or ())}


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
    """(report, fasta, sidecar, why there is no FASTA) -- what was given, else what the search
    recorded in its search_provenance.json (`result.report`, `fasta`), and <fasta>.meta.json."""
    why = None
    if search_out:
        try:
            with open(os.path.join(search_out, PROVENANCE)) as fh:
                prov = json.load(fh)
        except (OSError, ValueError) as e:
            prov = {}
            why = f"{PROVENANCE} in {search_out} could not be read ({type(e).__name__})"
        prov = prov if isinstance(prov, dict) else {}
        res = prov.get("result") if isinstance(prov.get("result"), dict) else {}
        cands = [report, res.get("report"), os.path.join(search_out, "report.parquet")]
        report = next((p for p in cands if p and os.path.isfile(p)), report or cands[-1])
        fasta = fasta or prov.get("fasta")
        if not fasta and not why:
            why = f"{PROVENANCE} in {search_out} names no FASTA"
    meta = None
    meta_path = fasta_meta or (fasta + ".meta.json" if fasta else None)
    if meta_path:
        try:
            with open(meta_path) as fh:
                meta = json.load(fh)
        except (OSError, ValueError):
            meta = None
    meta = meta if isinstance(meta, dict) else None
    fasta = fasta or (meta or {}).get("fasta")
    return report, fasta, meta, (None if fasta else why)


def keratin_exemption(meta, meta_path, report, keratin_sample=None):
    """(exemption or None, where it comes from) for a search not yet analysed: a keratin sample
    (search_provenance.json `keratin_sample`) whose database was not built with fetch_fasta.py
    --keratin-sample keeps the precursors of its keratin-family contaminant entries, which
    contaminants.R reads from fetch_fasta.keratin_database_check -- read here from the same place.
    A DE says its own (de_provenance.json contaminants.keratin_sample), and the audit uses that.
    `keratin_sample`: the search's own answer when the caller has it, else its record."""
    if keratin_sample is None:
        try:
            with open(os.path.join(os.path.dirname(os.path.abspath(report)), PROVENANCE)) as fh:
                sp = json.load(fh)
        except (OSError, ValueError, TypeError):
            sp = {}
        keratin_sample = sp.get("keratin_sample") if isinstance(sp, dict) else None
    if keratin_sample is not True:
        return None, None
    if keratin_sample_recorded(meta) is True:
        return None, "a keratin sample's database: its keratin entries were removed at build"
    q = keratin_database_check(meta_path, os.path.dirname(os.path.abspath(report)))
    if not q.get("checked"):
        return None, (f"a keratin sample, but its keratin-family entries are unknown "
                      f"({q.get('why')})")
    return {"shared": set(q["accessions"]), "alone": set(q["own_species"])}, q.get("why")


def cont_quant_exclude(report):
    """DIA-NN's --cont-quant-exclude for this report, as make_methods.py's one reader reads it
    (the DIA-NN log's command line, else the search's parameters file): {"value": tag or None,
    "source"} or {"not_recorded": why}."""
    try:
        from make_methods import diann_cont_quant_exclude_why
        got, why = diann_cont_quant_exclude_why(report)
    except Exception as e:           # its reasons are said, never a stop
        return {"not_recorded": f"the reader failed ({type(e).__name__}: {e})"}
    return got if isinstance(got, dict) else {"not_recorded": why or "nothing recorded it"}


# ---- the check ---------------------------------------------------------------------------------
def check(report, fasta=None, meta=None, q_cutoff=DEFAULT_Q_CUTOFF, q_columns=None,
          no_fasta_reason=None, exempt=None, meta_path=None, keratin_sample=None):
    """The record (a dict, always): status ok | not_run (with `reason`) | not_applicable. Never
    raises: this check reports, so whatever goes wrong is a not_run record that says what (a
    parquet column of an unexpected type made pyarrow raise ArrowNotImplementedError, review of
    2.14). `no_fasta_reason`: why there is no FASTA, when the caller knows. `exempt`: a keratin
    sample's exemption as the DE recorded it ({"shared", "alone"}; empty: none) -- None reads it
    from the search's records (keratin_exemption). `meta_path`: the sidecar's path.
    `keratin_sample`: whether this search is a keratin sample, when the caller knows it (an inline
    search, whose search_provenance.json is written after this check)."""
    rec = {"check": "contaminant_competition", "report": report, "fasta": fasta,
           "report_stat": report_stat(report)}
    try:
        return _check(rec, report, fasta, meta, q_cutoff, q_columns, no_fasta_reason, exempt,
                      meta_path, keratin_sample)
    except OutOfTime as e:
        rec.update(status="not_run", reason=str(e))
        return rec
    except Exception as e:
        name = os.path.basename(str(report or ""))
        why = " ".join(str(e).replace(str(report), name).split())
        rec.update(status="not_run",
                   reason=(f"the check failed on {name or 'the report'} ({type(e).__name__}: "
                           f"{why[:200]}{'...' if len(why) > 200 else ''})"))
        return rec


def _check(rec, report, fasta, meta, q_cutoff, q_columns, no_fasta_reason, exempt, meta_path,
           keratin_sample):
    rec.update({"organism": (meta or {}).get("organism") or None, "rule": RULE,
                # which fetch_fasta.py rule built the database (legacy / identity_only / current)
                "database_state": sidecar_state(meta) if meta else None})
    from estimate_params import DIANN_DIGEST
    dig = (meta or {}).get("contaminant_digest")
    digest = dig if isinstance(dig, dict) and {"cut", "missed_cleavages", "min_pep_len",
                                               "max_pep_len"} <= set(dig) else DIANN_DIGEST
    ev, why = read_report(report, q_cutoff, q_columns)
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
    # a keratin sample's keratin entries keep their precursors in the DE: not taken (review of 2.14)
    why = "as the DE recorded it"
    if exempt is None:
        meta_path = meta_path or (fasta + ".meta.json" if fasta else None)
        exempt, why = keratin_exemption(meta, meta_path if meta_path and os.path.isfile(meta_path)
                                        else None, report, keratin_sample)
    if exempt and exempt.get("shared"):
        for key in [k for k in ev["cont"] if exempted(k[0], exempt)]:
            ag = ev["cont"].pop(key)
            kept = ev["kept"].setdefault(key[0], _new())
            kept["precursors"] |= ag["precursors"]
            kept["signal"] += ag["signal"]
            for run, q in (ag.get("by_run") or {}).items():     # not a contaminant's signal
                if run in ev["runs"]:
                    ev["runs"][run][1] -= q
    rec["keratin_exempt"] = ({"shared": sorted(exempt["shared"]),
                              "alone": sorted(exempt.get("alone", ())), "source": why}
                             if exempt and exempt.get("shared") else
                             ({"shared": [], "source": why} if why else None))
    fa, rec["fasta_note"] = None, None
    if fasta:
        try:
            fa = read_fasta(fasta)
        except OSError as e:
            rec["fasta_note"] = f"the searched FASTA could not be read ({e})"
    else:
        rec["fasta_note"] = no_fasta_reason or "no searched FASTA was given"
    on_target = {}
    if fa:
        # the peptides that name no sample protein: is any of them a sample digest peptide?
        lone = {p for (ids, p) in ev["cont"] if not _targets(ids)}
        try:
            on_target = target_peptide_hits(lone, fa["targets"], digest)
        except ValueError as e:              # a --cut this digest cannot read
            rec["fasta_note"] = f"the search's digest could not be applied ({e})"
    if rec["fasta_note"]:
        rec["fasta_note"] += (": " + ("genes are named by accession, and " if fa is None else "")
                              + "a contaminant's peptide is checked against Protein.Ids only "
                                "(no I/L twins)")
    accs = (fa or {}).get("accs") or {}
    seeds = {a.upper() for ids, _p in ev["cont"] for a in _targets(ids)}
    seeds |= {h["target_acc"].upper() for h in on_target.values()}
    kept = connected(ev["kept"], seeds, accs, (fa or {}).get("by_gene") or {})
    rec.update(classify(ev, kept, accs, (fa or {}).get("cont_seq") or {}, on_target))
    rec["cont_quant_exclude"] = cont_quant_exclude(report)
    rec["status"] = "ok"
    return rec


def _gene(acc, accs):
    return (accs.get(acc.upper()) or ("", ""))[1] or acc


def connected(kept, seeds, accs, by_gene):
    """The kept rows (Protein.Ids keys of `kept`) connected to the seed accessions: a row naming
    one, then every row naming an accession those rows name, and a gene's other accessions, to
    closure. A paralog that shares a precursor with a taken gene brings its own peptides too, so
    the grouping sees all the evidence (review of 2.14: reading only the rows that named a seed
    let a paralog's shared peptides count for the taken gene)."""
    index = {}
    for ids in kept:
        for a in _targets(ids):
            index.setdefault(a.upper(), []).append(ids)
    seen, rows, todo = set(), set(), list(seeds)
    while todo:
        a = todo.pop()
        if a in seen:
            continue
        seen.add(a)
        todo.extend(by_gene.get((accs.get(a) or ("", ""))[1], ()))
        for ids in index.get(a, ()):
            if ids not in rows:
                rows.add(ids)
                todo.extend(x.upper() for x in _targets(ids))
    return rows


def _lead_key(gene):
    # the lead of genes the precursors cannot tell apart: the shortest name, then A-Z
    # (EEF1A1 before its pseudogene EEF1A1P5, ACTB before ACTBL2)
    return (len(gene), gene)


def parsimony(feats):
    """Which genes the precursors evidence, as protein inference groups them (the fewest genes
    that explain every precursor): repeatedly the gene that names the most precursors not yet
    explained takes all of them -- with every gene naming exactly the same ones, which the data
    cannot tell apart from it. A precursor shared by paralogs then counts for one gene, never for
    each paralog it also matches. Features hold disjoint precursors (a precursor has one
    Protein.Ids), so a gene's count is the sum of its features'; a heap with lazy updates keeps
    this fast for a large connected family. -> [{"genes": [lead, indistinguishable ...],
    "feats": [...]}], every feature in exactly one group."""
    size = [len(f["precursors"]) for f in feats]
    of_gene = {}
    for i, f in enumerate(feats):
        for g in f["genes"]:
            of_gene.setdefault(g, set()).add(i)
    count = {g: sum(size[i] for i in idx) for g, idx in of_gene.items()}
    heap = [(-n, _lead_key(g), g) for g, n in count.items()]
    heapq.heapify(heap)
    left = set(range(len(feats)))
    groups = []
    while left and heap:
        n, _key, g = heapq.heappop(heap)
        if -n != count[g]:                     # stale: counts only fall, so re-queue it
            heapq.heappush(heap, (-count[g], _lead_key(g), g))
            continue
        mine = of_gene[g] & left
        if not mine:
            continue
        cand = {h for i in mine for h in feats[i]["genes"]} - {g}
        same = sorted((h for h in cand if count[h] == count[g] and of_gene[h] & left == mine),
                      key=_lead_key)
        left -= mine
        for i in mine:
            for h in feats[i]["genes"]:
                count[h] -= size[i]
        groups.append({"genes": [g] + same, "feats": [feats[i] for i in sorted(mine)]})
    return groups


def classify(ev, kept_rows, accs, cont_seq, on_target):
    """The findings from the evidence: the share, which contaminant entries show themselves, and
    the sample genes COMPETED OUT / SIPHONED / masked (and how many only lost a minority)."""
    # 1. each contaminant precursor's entries and sample proteins (an I/L twin of a sample digest
    # peptide names that protein), and the entries' own peptides
    own, rows = {}, []
    for (ids, pep), ag in ev["cont"].items():
        conts = [a for a in _ids(ids) if is_contaminant_accession(a)]
        targets = _targets(ids)
        twin = on_target.get(pep) if not targets else None
        if twin:
            targets = [twin["target_acc"]]
            twin = twin if twin["how"] == "il_twin" else None
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

    # 2. the features: the sample's precursors the rule takes ("lost"), those shared with an entry
    # that shows itself ("masked"), and the connected kept (contaminant-free) ones -- each with
    # the genes it names
    feats = []
    for conts, targets, pep, ag, twin in rows:
        if targets:
            feats.append({"kind": "masked" if any(c in shows for c in conts) else "lost",
                          "genes": {_gene(t, accs) for t in targets}, "accs": set(targets),
                          "conts": set(conts), "twin": pep if twin else None, **ag})
    for ids in sorted(kept_rows):
        targets = _targets(ids)
        feats.append({"kind": "kept", "genes": {_gene(t, accs) for t in targets},
                      "accs": set(targets), "conts": set(), "twin": None, **ev["kept"][ids]})

    # 3. the genes, grouped as protein inference groups them; each feature belongs to one group
    groups = parsimony(feats)
    owner = {id(f): gi for gi, grp in enumerate(groups) for f in grp["feats"]}
    of_gene = {}
    for i, f in enumerate(feats):
        for g in f["genes"]:
            of_gene.setdefault(g, []).append(i)

    def entry(gi, grp):
        lost = [f for f in grp["feats"] if f["kind"] == "lost"]
        kept = [f for f in grp["feats"] if f["kind"] == "kept"]
        masked = [f for f in grp["feats"] if f["kind"] == "masked"]
        ls, ks = sum(f["signal"] for f in lost), sum(f["signal"] for f in kept)
        mine = set(grp["genes"])
        # precursors naming this gene that another gene explains (a paralog's): not this gene's,
        # so "every peptide of it" would be false (review of 2.14)
        elsewhere = [feats[i] for i in {i for g in mine for i in of_gene.get(g, ())}
                     if owner[id(feats[i])] != gi]
        return {"gene": grp["genes"][0], "indistinguishable": grp["genes"][1:],
                "accessions": sorted({a for f in grp["feats"] for a in f["accs"]
                                      if _gene(a, accs) in mine}),
                "contaminants": sorted({c for f in lost for c in f["conts"]}),
                "n_removed_precursors": sum(len(f["precursors"]) for f in lost),
                "n_kept_precursors": sum(len(f["precursors"]) for f in kept),
                "n_masked_precursors": sum(len(f["precursors"]) for f in masked),
                "masked_signal": sum(f["signal"] for f in masked),
                "masked_by": sorted({c for f in masked for c in f["conts"]}),
                "n_precursors_explained_by_other_genes":
                    sum(len(f["precursors"]) for f in elsewhere),
                "other_genes": sorted({groups[owner[id(f)]]["genes"][0] for f in elsewhere}),
                "removed_signal_pct": round(100 * ls / (ls + ks), 1) if ls + ks else None,
                "removed_signal": ls, "kept_signal": ks,
                **({"il_twin_peptides": sorted({f["twin"] for f in lost if f["twin"]})}
                   if any(f["twin"] for f in lost) else {})}
    competed, siphoned, masked, minority = [], [], [], 0
    for gi, grp in enumerate(groups):
        if all(f["kind"] == "kept" for f in grp["feats"]):
            continue                     # a connected protein the contaminant rule left alone
        e = entry(gi, grp)
        if e["n_removed_precursors"]:
            if not e["n_kept_precursors"]:
                competed.append(e)
            elif e["removed_signal"] > e["kept_signal"]:
                siphoned.append(e)
            else:
                minority += 1
        elif e["n_masked_precursors"] and e["masked_signal"] >= e["kept_signal"]:
            # most of it, or all, rests on precursors shared with an entry that shows itself:
            # cannot be attributed, so never "lost" -- but said (review of 2.14: a gene that kept a
            # few peptides was reported nowhere)
            masked.append({k: e[k] for k in ("gene", "indistinguishable", "n_kept_precursors")}
                          | {"contaminants": e["masked_by"],
                             "n_shared_precursors": e["n_masked_precursors"]})
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
            "n_kept_protein_sets_read": len(kept_rows),
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


def fate(rec, de=None):
    """(staff sentence, taken) -- what was done with the shared precursors, as the records say
    (DE-LIMP rule 1: never one pipeline's words for every case). `de`: the DE's own contaminant
    record (de_provenance.json `contaminants`) when there is one; else DIA-NN's
    --cont-quant-exclude as the search recorded it. taken: True when they were left out of the
    quantities, else False (the client then reads the neutral words)."""
    if isinstance(de, dict) and de.get("policy") in ("removed", "kept"):
        if de.get("removed"):
            owner = (str(de.get("rule_origin") or "").strip()
                     or "[the filter's origin not recorded -- confirm]")
            return (f"The DE removed every contaminant precursor before quantification ({owner}), "
                    "so they are not in its results.", True)
        return ("The DE kept contaminant precursors (--keep-contaminants), so it did not remove "
                "these; which protein group DIA-NN put each one in is in the report's "
                "Protein.Group.", False)
    cq = rec.get("cont_quant_exclude") if isinstance(rec.get("cont_quant_exclude"), dict) else {}
    if cq.get("value"):
        return (f"DIA-NN's --cont-quant-exclude {cq['value']} ({cq.get('source')}) left them out "
                "of its own quantities and normalisation, and a DE that applies the same rule "
                "leaves them out too.", True)
    if "value" in cq:
        return (f"DIA-NN ran without --cont-quant-exclude ({cq.get('source')}), so its own "
                "quantities kept them; a DE that removes contaminant precursors leaves them out.",
                False)
    return ("Whether DIA-NN left them out of its own quantities: [not recorded -- confirm]"
            + (f" ({cq['not_recorded']})" if cq.get("not_recorded") else "") + ".", False)


def _others(x):
    k, who = x.get("n_precursors_explained_by_other_genes") or 0, x.get("other_genes") or []
    if not k:
        return ""
    return (f"; {_n(k, 'more precursor', 'more precursors')} shared with "
            f"{', '.join(who[:3])}{' ...' if len(who) > 3 else ''}, counted there")


def competed_text(rec):
    """COMPETED OUT, for staff. What was then done with the precursors is fate()'s to say."""
    co, org = rec.get("competed_out") or [], rec.get("organism") or "sample"
    if not co:
        return ""
    one = len(co) == 1
    return (f"COMPETED OUT: {_n(len(co), org + ' gene has', org + ' genes have')} no precursor "
            f"of {'its' if one else 'their'} own left: every precursor that identifies "
            f"{'it' if one else 'them'} and passes the q-value filter (one shared with another "
            "sample protein counts for that protein) is shared with a contaminant entry that "
            "shows no peptides of its own here, so where the contaminant rule (any contaminant "
            f"accession in Protein.Ids) is applied {'it has' if one else 'they have'} no "
            "quantity: "
            + "; ".join(f"{x['gene']} ({_same(x)}"
                        f"{_n(x['n_removed_precursors'], 'precursor', 'precursors')} via {_via(x)}"
                        f"{_others(x)})" for x in co[:LIST_MAX])
            + (f"; +{len(co) - LIST_MAX} more" if len(co) > LIST_MAX else "") + ".")


def siphoned_text(rec):
    """SIPHONED, for staff. What was then done with the precursors is fate()'s to say."""
    si, org = rec.get("siphoned") or [], rec.get("organism") or "sample"
    if not si:
        return ""
    one = len(si) == 1
    return (f"SIPHONED: {_n(len(si), org + ' gene', org + ' genes')} kept precursors of "
            f"{'its' if one else 'their'} own, but most of {'its' if one else 'their'} signal "
            "is on precursors shared with an entry that shows no peptides of its own here, which "
            f"the contaminant rule, where it is applied, takes from {'it' if one else 'them'}: "
            + "; ".join(f"{x['gene']} ({_same(x)}{x['removed_signal_pct']:g}% of its signal "
                        f"taken; {x['n_removed_precursors']} taken, "
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
        out.append(f"Not counted as lost: {_n(len(ma), 'gene is', 'genes are')} seen only, or "
                   "mostly, through precursors shared with a contaminant entry that does show "
                   "peptides of its own, so the data cannot say how much of that signal is the "
                   f"sample's: {_names(ma)}.")
    if (rec.get("keratin_exempt") or {}).get("shared"):
        out.append(f"A keratin sample: precursors of its {len(rec['keratin_exempt']['shared'])} "
                   "keratin-family contaminant entries are kept, as the DE keeps them "
                   f"({rec['keratin_exempt'].get('source')}).")
    if rec.get("fasta_note"):
        out.append(f"({rec['fasta_note']}.)")
    s = rec.get("contaminant_share_pct") or {}
    if s.get("median") is not None:
        out.append(f"Contaminant entries carry a median {s['median']:g}% (max {s['max']:g}%) of "
                   f"a run's precursor signal ({s.get('column')}, {rec.get('n_runs')} runs).")
    return " ".join(out)


def staff_message(rec, de=None):
    """The technical record: AUDIT.md, the job log, the watcher. Check codes and accessions
    allowed; never the client's words."""
    st = rec.get("status")
    if st == "not_run":
        return f"Not run: {rec.get('reason')}."
    if st == "not_applicable":
        return f"Not applicable: {rec.get('reason')}."
    head = " ".join(x for x in (competed_text(rec), siphoned_text(rec)) if x)
    if head:
        head += " " + fate(rec, de)[0]             # once, for both
    else:
        org = rec.get("organism") or "sample"
        quiet = len(rec.get("contaminants_not_showing_themselves") or [])
        head = (f"Checked: no {org} gene lost all or most of its signal to a contaminant entry "
                f"({_n(quiet, 'contaminant entry', 'contaminant entries')} seen without "
                "peptides of their own).")
    return " ".join(x for x in (head, fix_text(rec), notes_text(rec)) if x)


def plain_competed(rec, de=None):
    co = rec.get("competed_out") or []
    one = len(co) == 1
    it, its = ("it", "its") if one else ("them", "their")
    # "the peptides that tell it apart": true also when a paralog explains its other peptides
    if fate(rec, de)[1]:
        return (f"{_names(co)} {'is' if one else 'are'} missing from the results because the "
                f"contaminant list absorbed {its} peptides: the peptides that tell {it} apart from "
                "other proteins are all shared with an entry of the contaminant list, and those "
                "peptides are left out of the measurements. The Core can search again without "
                f"those contaminant entries to recover {it}.")
    return (f"The peptides that tell {_names(co)} apart from other proteins are all shared with an "
            f"entry of the contaminant list, so the data cannot tell {its} signal from that "
            f"entry's. The Core can search again without those contaminant entries to measure "
            f"{it} on {its} own.")


def plain_siphoned(rec, de=None):
    si = rec.get("siphoned") or []
    one = len(si) == 1
    shares, its = ("it shares", "its") if one else ("they share", "their")
    if fate(rec, de)[1]:
        return (f"Most of the signal of {_names(si)} was on peptides {shares} with an entry of the "
                "contaminant list, which are left out of the measurements, so "
                f"{'it is' if one else 'they are'} measured from {its} remaining peptides only. "
                "The Core can search again without those contaminant entries.")
    return (f"Most of the signal of {_names(si)} was on peptides {shares} with an entry of the "
            f"contaminant list, so that part of {its} signal cannot be told apart from the "
            "contaminant. The Core can search again without those contaminant entries.")


def audit_findings(rec, de=None):
    """[(status, message, detail, plain)] for audit_results.py; [] for a record that does not
    apply (another engine's report). `de`: the DE's contaminant record, which says what was done
    with the shared precursors. COMPETED OUT and SIPHONED are each a WARN, never a FAIL: the
    check reports and never stops an analysis. A check that could not run is an INFO saying why.
    plain None = the check's usual words (audit_results.PLAIN)."""
    if not isinstance(rec, dict) or rec.get("status") == "not_applicable":
        return []
    if rec.get("status") != "ok":
        return [("INFO", staff_message(rec), {"reason": rec.get("reason"),
                                              "report": rec.get("report")}, None)]
    keys = ("report", "fasta", "organism", "contaminant_share_pct", "drop_contaminant",
            "contaminants_not_showing_themselves", "q_columns", "cutoffs", "record_file",
            "checked_by")
    detail = {k: rec.get(k) for k in keys if rec.get(k) is not None}
    detail.update({k: [x["gene"] for x in rec.get(k) or []]
                   for k in ("competed_out", "siphoned", "masked")})
    tail = " ".join(x for x in (fix_text(rec), notes_text(rec)) if x)
    said = fate(rec, de)[0]
    out = []
    if rec.get("competed_out"):
        out.append(("WARN", f"{competed_text(rec)} {said} {tail}", detail, plain_competed(rec, de)))
    if rec.get("siphoned"):
        # the fix and the notes once, in the first of the two
        out.append(("WARN", f"{siphoned_text(rec)} {said}" + ("" if out else f" {tail}"), detail,
                    plain_siphoned(rec, de)))
    if out:
        return out
    # a gene whose signal rests mostly on an entry that shows itself (bovine serum proteins in
    # cell culture): nothing attributable was taken, but the client is told what that means here
    return [("PASS", staff_message(rec, de), detail,
             plain_masked(rec, de) if rec.get("masked") else None)]


def plain_masked(rec, de=None):
    """The client's words for genes resting mostly on a contaminant that is present (masked), as
    fate() says what was done with the shared peptides -- never "no protein lost its signal" when
    a removed one left a protein with nothing (review of 2.14)."""
    head = ("Some proteins share peptides with a contaminant that is present in the samples "
            "(such as serum proteins from cell culture), so their signal cannot be told apart "
            "from it")
    if fate(rec, de)[1]:
        return (head + ": those shared peptides are left out of the measurements, so these "
                "proteins are measured from their remaining peptides only, and one with none "
                "left is not in the results. No other protein lost its peptides to the "
                "contaminant list.")
    return head + ". No other protein lost its peptides to the contaminant list."


def say(rec, stream=None):
    """The job log's line: WARNING when a gene was lost, else INFO (NOTE when it could not run)."""
    level = ("WARNING" if rec.get("competed_out") or rec.get("siphoned") else
             "NOTE" if rec.get("status") == "not_run" else "INFO")
    print(f"{TAG} {level}: {staff_message(rec)}", file=stream or sys.stderr, flush=True)


# ---- the record --------------------------------------------------------------------------------
def write(rec, out=None, provenance=None):
    """<out> (JSON) and search_provenance.json `contaminant_competition` (merged in place: a .tmp
    then a rename, so a reader never sees half a file). -> ([files written], [what could not be
    written, and why]); never raises."""
    rec["message"] = staff_message(rec)            # what the watcher and a person read
    written, problems = [], []
    if out:
        rec["record_file"] = os.path.abspath(out)
        try:
            with open(out + ".tmp", "w") as fh:
                json.dump(rec, fh, indent=2)
            os.replace(out + ".tmp", out)
            written.append(out)
        except (OSError, TypeError, ValueError) as e:
            problems.append(f"{out} ({e})")
    if provenance and os.path.isfile(provenance):
        try:
            with open(provenance) as fh:
                prov = json.load(fh)
            prov["contaminant_competition"] = rec
            with open(provenance + ".tmp", "w") as fh:
                json.dump(prov, fh, indent=2)
            os.replace(provenance + ".tmp", provenance)
            written.append(provenance)
        except (OSError, TypeError, ValueError, AttributeError) as e:
            problems.append(f"{provenance} ({e})")
    return written, problems


def _slurm_seconds(text):
    """squeue's %L ([days-]hours:minutes:seconds, minutes:seconds, or seconds) -> seconds, or
    None (UNLIMITED, NOT_SET, INVALID, nothing)."""
    try:
        days, _, rest = text.strip().rpartition("-")
        parts = [int(x) for x in rest.split(":")]
        if not 1 <= len(parts) <= 3:
            return None
        parts = [0] * (3 - len(parts)) + parts
        return int(days or 0) * 86400 + parts[0] * 3600 + parts[1] * 60 + parts[2]
    except ValueError:
        return None


def job_margin():
    """Seconds the job still needs after the check: the job-end hook's own limit
    (notify_slack.JOB_END_HOOK_S: FRAN staging, the run log, Slack) and a minute."""
    from notify_slack import JOB_END_HOOK_S
    return JOB_END_HOOK_S + JOB_MARGIN_EXTRA_S


def job_budget(limit):
    """(seconds this check may take, what bounded it): at most `limit`, and inside a SLURM job
    never into the time the job still needs after it (job_margin; squeue -o %L for
    SLURM_JOB_ID). Never raises: what could not be read is said."""
    jid = os.environ.get("SLURM_JOB_ID")
    if not jid:
        return limit, f"the check's own limit of {limit} s"
    try:
        import shutil
        from on_compute import _slurm_bin
        squeue = shutil.which("squeue") or os.path.join(_slurm_bin(os.environ), "squeue")
        got = subprocess.run([squeue, "-h", "-j", jid, "-o", "%L"], capture_output=True,
                             text=True, timeout=30).stdout.strip()
        margin = job_margin()
    except Exception as e:           # noqa: BLE001 -- the job's time is then not known: said
        return limit, f"the check's own limit of {limit} s (the job's time left: {e})"
    left = _slurm_seconds(got.splitlines()[0] if got else "")
    if left is None or left - margin >= limit:
        return limit, f"the check's own limit of {limit} s"
    return left - margin, f"{left} s were left in job {jid}, {margin} s of them for its end"


def set_aside_record(out):
    """A search generated into `out`: an earlier search's record there does not describe it --
    set aside (diann_parallel.set_aside: renamed .stale-<time>, never deleted), and said. The
    watcher would otherwise show the earlier search's genes while this one queues and runs (review
    of 2.14). -> the new path, or None."""
    from diann_parallel import set_aside
    path = os.path.join(out, RECORD)
    try:
        new = set_aside(path)
    except (ValueError, OSError) as e:   # never a reason to stop a generated search: it is said
        sys.stderr.write(f"{TAG} WARNING: {e}; left in place\n")
        return None
    if new:
        sys.stderr.write(f"{TAG} an earlier search's {RECORD} in {out} does not describe this "
                         f"search: set aside as {os.path.basename(new)}\n")
    return new


def job_lines(report, fasta, out_dir):
    """The bash a DIA-NN search job ends with (diann_parallel.py step 5, run_search.py's
    single-shot job): this check on its report, under the generator's own python (the pipeline
    env's, which has pyarrow), else python3, for at most JOB_TIME_LIMIT seconds where `timeout`
    exists. It never fails the job; a check that does not finish is said in the log, and the
    audit then checks the report itself."""
    report, fasta, out_dir = (os.path.abspath(p) for p in (report, fasta, out_dir))
    # --job-budget stops it in time and records why; `timeout` is only the backstop
    hard = JOB_TIME_LIMIT + JOB_MARGIN_EXTRA_S
    return [f'PY={shlex.quote(sys.executable)}; [ -x "$PY" ] || PY=python3; '
            f'T=""; command -v timeout >/dev/null 2>&1 && T="timeout {hard}"; '
            f'$T "$PY" {shlex.quote(os.path.abspath(__file__))} --report {shlex.quote(report)} '
            f'--fasta {shlex.quote(fasta)} --out {shlex.quote(os.path.join(out_dir, RECORD))} '
            f'--provenance {shlex.quote(os.path.join(out_dir, PROVENANCE))} '
            f'--job-budget {JOB_TIME_LIMIT} '
            '|| echo "WARNING: the contaminant competition check did not finish (exit $?; at most '
            f'{JOB_TIME_LIMIT // 60} min in the search job) -- the messages above say why; the '
            'audit checks the report itself" >&2']


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
    ap.add_argument("--job-budget", type=int,
                    help="inside a search job: stop after this many seconds, or two minutes "
                         "before the job's own time limit, whichever comes first; the record "
                         "then says the check did not run, and why")
    a = ap.parse_args(argv)
    if not (a.report or a.search_out):
        ap.error("give --report or --search-out")
    report, fasta, meta, why = resolve_inputs(a.search_out, a.report, a.fasta, a.fasta_meta)
    budget, bound = job_budget(a.job_budget) if a.job_budget else (None, None)
    if budget is not None and budget < MIN_BUDGET_S:
        rec = {"check": "contaminant_competition", "report": report,
               "report_stat": report_stat(report), "status": "not_run",
               "reason": f"too little of the search job's time was left for it ({bound})"}
    else:
        if budget:
            def stop(_signum, _frame):
                raise OutOfTime(f"stopped after {budget} s, before the search job's time ran out "
                                f"({bound}); the audit checks the report itself")
            signal.signal(signal.SIGALRM, stop)
            signal.alarm(budget)
        try:
            rec = check(report, fasta, meta, a.q_cutoff, no_fasta_reason=why,
                        meta_path=a.fasta_meta)
        finally:
            if budget:
                signal.alarm(0)
                signal.signal(signal.SIGALRM, signal.SIG_DFL)
    out = a.out or (os.path.join(a.search_out, RECORD) if a.search_out else None)
    prov = a.provenance or (os.path.join(a.search_out, PROVENANCE) if a.search_out else None)
    written, problems = write(rec, out, prov)
    say(rec)
    for p in problems:
        print(f"{TAG} NOTE: could not write {p}.", file=sys.stderr)
    if written:
        print(f"{TAG} recorded in {', '.join(written)}", file=sys.stderr)
    summary = {k: rec.get(k) for k in ("status", "reason", "drop_contaminant")}
    summary.update({k: [x["gene"] for x in rec.get(k) or []]
                    for k in ("competed_out", "siphoned", "masked")})
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
