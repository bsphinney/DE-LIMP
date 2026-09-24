"""
fran_deposit.py  --  Hand a finished HIVE search to FRAN, automatically.

Every search this skill runs on HIVE **for the UC Davis Proteomics Core** should end up in
FRAN (the Core's cross-search corpus, https://fran.stan-proteomics.org). Doing that by hand
means it happens for the searches somebody remembered and for none of the others, which is
how a corpus develops holes that look like biology. So the orchestrator calls this at the
end of every HIVE search and it decides for itself whether the run is eligible.

**This script does not ingest anything.** FRAN runs its own cron on HIVE that scans for new
searches and ingests them (`ingest/find_uningested.py` -> `auto_ingest.py`). All the skill has
to do is put the search where that scan will find it. So this **stages** the search into
FRAN's drop directory and stops. No database, no credential, no SLURM job — which is why it
works for every Core member and not just whoever owns the corpus token.

  # 1. is this run eligible, and what would be handed over?  (login-node safe)
  python3 fran_deposit.py check --out <search out dir>
  # 2. stage it for the cron (symlinks -- nothing is copied)
  python3 fran_deposit.py stage --out <search out dir> [--name "..."] \
      [--organism "Homo sapiens" --taxon 9606] [--fasta-meta <search.fasta.meta.json>]
  # 3. later: did the cron actually take it?  (DB-free: answers from the cron's own logs)
  python3 fran_deposit.py verify --out <search out dir>
  # can a staged search reach FRAN right now?  (logs + a few small files; login-node safe)
  python3 fran_deposit.py health [--alert]
  # past searches that were never staged  (dry run; --sbatch writes the job, --apply stages)
  python3 fran_deposit.py backfill --sbatch

All of them run **on HIVE** (in hive_remote mode through `hive_exec.sh`, like every other
HIVE-side step). Each prints one JSON object on stdout; warnings go to stderr.

## Staged is not ingested

Staging only puts the search where FRAN's cron looks. Whether the cron is actually getting to
it is a separate fact, and it has been false for a week at a time while every staged entry
looked fine: the cron kept running, ingested nothing, and never reached the skill's entries.
`health` reads what the cron itself wrote (its `auto_ingest_<jobid>.out` logs and the submit
log) and compares the ingest code HIVE runs with FRAN's GitHub `main`, because that code
changes often and a stale copy writes rows that look fine. `stage` prints the one-line result
as a warning when the answer is bad; it never fails because of it.

## Symlinks, never copies

A search directory is tens of GB. The drop entry is a **real directory** holding **symlinks**
to the search's outputs — `report.parquet`, `report_xic/` (DIA-NN `--xic` chromatograms),
the logs and `search_provenance.json` — plus a small real `fran_manifest.json`.

It has to be a real directory of links rather than one symlink to the search dir, because
FRAN's scanner walks with `os.walk(..., followlinks=False)`: a symlinked *directory* is never
descended into, so a bare symlink would be silently invisible to the cron. A real directory
is walked normally and `os.path.exists()` follows the links inside it, so the scanner detects
the engine exactly as it would in the original tree.

The entry name is **deterministic** — `<search dir name>__<8 hex of its real path>` — so
re-staging the same search reuses the same path instead of creating a second candidate.

## Who gets staged — and who must not

Only **Proteomics Core** searches. A collaborator with their own HIVE account is running their
own data through this skill; their results are theirs and never enter the Core corpus. The
gate is membership of the `proteomics-grp` group, tested by whether the drop directory (which
lives inside `/quobyte/proteomics-grp`) is writable. A non-member cannot stage: the group
permission on the directory is the enforcement, not a flag in this file.

`check` is refused (with a reason code, never an exception) when:
  not_on_hive          the search did not run on HIVE -- nothing to stage from here
  search_incomplete    no report, or a zero-byte one: a failed search must never be ingested
  engine_unsupported   FRAN's corpus schema covers the DIA engines (DIA-NN / FragPipe /
                       Radiant / Spectronaut). Sage + AlphaDIA runs are skipped, not faked.
  not_core_facility    the drop directory is not writable -> not a Core account
  no_drop_dir          FRAN's drop directory does not exist and cannot be created
  already_staged       this search is already staged (--force to redo)
  opted_out            FRAN_DEPOSIT=off, or --skip
  qc_run               a QC run (is_qc_run: --qc/--not-qc, session `qc`, FRAN's QC trees, or
                       FRAN's QC name rule). Never staged; one staged earlier is marked qc: true
`backfill` adds three of its own, for a directory it found but will not hand over:
  not_a_skill_search   a search this skill did not run (e.g. a DE-LIMP app search)
  already_ingested     the cron's logs show FRAN already ingested it (by any route)
  needs_agent_check    a FragPipe/Radiant search with no completion marker to read
Not decisions, so never recorded -- permission problems, reported and retried next time:
  drop_dir_not_writable  a proteomics-grp member cannot write the drop dir
  entry_not_writable     the entry was staged by another account without group write

Set `FRAN_DEPOSIT=off` (or pass `--skip`) to keep a run out of FRAN. `stage` writes that
decision into the receipt (`status: opted_out`), so a later `backfill` honours it instead of
handing the search over anyway.

## Pass the organism through — it is the one thing only the skill knows

A DIA-NN `report.parquet` carries **no organism column**, so a search ingested without one
sits in the corpus with `organism` NULL and is invisible on FRAN's species page. The skill
already had the user *confirm* the organism, and `fetch_fasta.py` wrote it to
`<fasta>.meta.json` — read automatically here, overridable with `--organism`/`--taxon`, and
written into `fran_manifest.json` for the cron to pass to `corpus_ingest.py`. Never invented:
an unresolved organism stays absent from the manifest rather than being guessed.

Overrides: FRAN_DROP_DIR, FRAN_DEPOSIT=off, DELIMP_PG_TOKEN_FILE (verify only),
FRAN_INGEST_DIR, FRAN_INGEST_LOG_DIR, FRAN_HEALTH=off (stage does not even read the status
file), FRAN_HEALTH_FILE (where `health` writes it and `stage` reads it; default
/quobyte/proteomics-grp/fran/ingest_health.json), FRAN_HEALTH_CACHE.
"""
import argparse
import collections
import datetime
import getpass
import glob
import hashlib
import json
import math
import os
import re
import shlex
import subprocess
import sys
import threading
import time
import urllib.error
import urllib.request

GROUP_ROOT = "/quobyte/proteomics-grp"
FRAN_URL = "https://fran.stan-proteomics.org"
RECEIPT = "fran_deposit.json"
MANIFEST = "fran_manifest.json"
# FRAN's drop directory: what its ingest cron scans. Inside the group root on purpose -- the group
# write permission IS the Core-member gate, so a collaborator cannot stage even by accident.
DROP_DIR = os.environ.get("FRAN_DROP_DIR", f"{GROUP_ROOT}/fran/incoming")
# Receipt statuses that mean "FRAN already has this search, or is about to".
BLOCKING_STATUS = {"staged", "ingested"}

# Used ONLY by `verify`, to ask the corpus directly whether the cron has ingested a staged search.
# Staging needs none of this -- which is the point of the drop-directory design: a Core member with
# no corpus credential can still hand a search over. When these are unavailable, verify reports
# what it could see (the drop entry) and says so, rather than calling a pending search a failure.
INGEST_DIRS = [os.environ.get("FRAN_INGEST_DIR"), f"{GROUP_ROOT}/brett/glendon/fran_ingest"]
PY_CANDIDATES = [os.environ.get("FRAN_INGEST_PYTHON"),
                 f"{GROUP_ROOT}/brett/envs/alphadia2/bin/python"]
TOKEN_CANDIDATES = [os.environ.get("DELIMP_PG_TOKEN_FILE"),
                    f"{GROUP_ROOT}/fran/.pgfarm_token",
                    os.path.expanduser("~/.pgfarm_token")]

# What gets linked into the drop entry. A curated list, not the whole directory: the search dir
# also holds tens of GB of per-file .quant temporaries and a predicted library, and linking those
# would make the handover unreadable to a person without helping either ingester. Everything the
# corpus ingest, the XIC lane, or a later audit actually reads is here.
LINK_ITEMS = [
    "report.parquet", "report.tsv", "delimp_report.parquet",   # the quant of record
    "report.log.txt", "report.stats.tsv",                      # what the engine says it did
    "dia-quant-output",                                        # FragPipe's DIA-NN output
    "radiant_results", "fulcrum-results",                      # Radiant / Fulcrum
]
# `report_xic` is NOT in this list: stage() builds it from find_xic_dirs(), because the parallel
# chain leaves chromatograms in `xic/t<N>_xic/` rather than `report_xic/` and a plain link would
# miss every one of them.
#
# NOT linked, deliberately: `search_provenance.json`. FRAN's scanner lists it as a **Radiant**
# marker and tests Radiant before DIA-NN, so an entry containing it is detected as Radiant whatever
# engine actually ran -- and run_search.py writes one into EVERY search directory. Linking it would
# have relabelled every staged DIA-NN and FragPipe search as Radiant. Its contents are embedded in
# fran_manifest.json instead, where nothing can mistake them for a marker.
#
# What each engine is then detected by, with this link set (verified against FRAN's ENGINE_MARKERS):
#   DIA-NN    report.parquet / report.tsv
#   FragPipe  dia-quant-output/report.tsv   (tested first, so it wins over the DIA-NN report inside it)
#   Radiant   radiant_results/fulcrum-results, or fulcrum-results/_SUCCESS

# search out dir -> (engine, report path relative to it). Mirrors corpus_ingest.ingest()'s own
# candidate list; kept in the same order so `check` reports the file the ingester will pick.
ENGINE_REPORTS = {
    "diann":    ("report.parquet", "report.tsv"),
    "radiant":  ("radiant_results/fulcrum-results", "fulcrum-results", "delimp_report.parquet"),
    "fragpipe": ("dia-quant-output/report.tsv", "report.tsv"),
}
# Engines FRAN's corpus schema does not model. Skipping them is the honest outcome; inventing a
# DIA row for a DDA search would put a wrong acquisition label on real data.
UNSUPPORTED = {"sage": "Sage is DDA; the FRAN corpus is a DIA corpus",
               "alphadia": "AlphaDIA output is not one of corpus_ingest.py's supported schemas"}

# Markers that IDENTIFY an engine, in order — distinct from ENGINE_REPORTS, which says where that
# engine's report lives. The two cannot be the same list: FragPipe's report is also called
# `report.tsv`, so testing "fragpipe first" against the report candidates would relabel every
# DIA-NN 1.9 search (report.tsv, no parquet) as FragPipe. Order and specificity both matter —
# FragPipe's tree CONTAINS a DIA-NN report, so it must be tested first and only on markers unique
# to it. TWIN: FRAN ingest/find_uningested.py ENGINE_MARKERS (c4838fe), so both ends agree on what
# a directory is -- with ONE deliberate difference: FRAN also lists search_provenance.json under
# Radiant, but run_search.py writes that file for EVERY engine, and detect_engine() reads its
# `engine` field before sniffing, so it is never a marker here. delimp_report.parquet is a Radiant
# REPORT location (ENGINE_REPORTS, as corpus_ingest's candidates) but not a marker on either side.
DETECT_MARKERS = [
    ("fragpipe", ("dia-quant-output/report.tsv", "fragpipe.fp-manifest")),
    ("radiant",  ("radiant_results/fulcrum-results", "fulcrum-results/_SUCCESS")),
    ("diann",    ("report.parquet", "report.tsv")),
]


def jout(d, code=0):
    print(json.dumps(d, indent=2))
    sys.exit(code)


def first_readable(paths, isdir=False, executable=False):
    """First of `paths` this account can actually use. `executable` matters for the interpreter:
    a python that is readable but not executable passes an isfile() test and then fails inside the
    job, an hour later and in a log nobody is watching."""
    need = os.R_OK | (os.X_OK if executable else 0)
    for p in paths:
        if not p:
            continue
        p = os.path.expanduser(p)
        if (os.path.isdir(p) if isdir else os.path.isfile(p)) and os.access(p, need):
            return p
    return None


def dir_size(path):
    """Bytes under a path -- a Spark/Fulcrum 'report' is a DIRECTORY of parquet parts, so a
    plain getsize() on it returns the inode size and a broken result looks non-empty."""
    if os.path.isfile(path):
        return os.path.getsize(path)
    tot = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                tot += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return tot


def detect_engine(out):
    """Engine that produced this search dir. search_provenance.json is written by run_search.py
    for every run and is authoritative; the file sniff is the fallback for a dir we did not
    produce (an older run, or one a person made by hand)."""
    prov = os.path.join(out, "search_provenance.json")
    if os.path.isfile(prov):
        try:
            with open(prov) as fh:
                p = json.load(fh)
            eng = (p.get("engine") or "").lower()
            if eng:
                return eng, p.get("version"), "search_provenance.json"
        except (OSError, ValueError):
            pass                       # a truncated provenance file is not a reason to stop
    for eng, markers in DETECT_MARKERS:
        for mk in markers:
            if os.path.exists(os.path.join(out, mk)):
                return eng, None, f"found {mk}"
    return None, None, "no engine could be determined"


def completion_marker(out, engine, report):
    """(done, evidence) for engines whose report can exist before the search FINISHED.
      done True   the engine's own completion marker is there
      done False  its marker is missing, or its log says it was cancelled: a partial search
      done None   no reliable marker for this layout -- an agent must confirm it finished
    DIA-NN needs none here: every skill route writes the report last, behind must_exist().

    Found by LOOKING at every Radiant and FragPipe output under the Core trees on HIVE (read-only
    find, job 23992239, 2026-09-24):
      FragPipe  16 workdirs; dia-quant-output/report.tsv exists in ALL of them, including 2 whose
                log says "Cancelling N remaining tasks" and 3 that never finished -- so a report
                is no evidence at all. The newest log_<date>.txt in the workdir ends "ALL JOBS
                DONE IN <n> MINUTES" in the 11 that finished, and in none of the other 5.
      Fulcrum   20 fulcrum-results/ dirs, all carrying _SUCCESS: the Spark output committer's
                marker, written only once the parquet write commits (FRAN's ENGINE_MARKERS use it
                too). A Radiant search with only delimp_report.parquet has no such marker."""
    if engine == "radiant":
        if report and os.path.isdir(report) and os.path.basename(report.rstrip("/")) == "fulcrum-results":
            ok = os.path.isfile(os.path.join(report, "_SUCCESS"))
            return ok, (f"{report}/_SUCCESS present" if ok else
                        f"{report} has no _SUCCESS: the Fulcrum write never committed")
        return None, (f"Radiant output without fulcrum-results/ ({os.path.basename(report or '')}): "
                      f"no completion marker to check")
    if engine == "fragpipe":
        logs = sorted(glob.glob(os.path.join(out, "log_*.txt")))    # log_YYYY-MM-DD_HH-MM-SS.txt
        if not logs:
            return None, f"no FragPipe log_*.txt in {out}: no completion marker to check"
        try:
            with open(logs[-1], "rb") as fh:
                fh.seek(0, os.SEEK_END)
                fh.seek(max(0, fh.tell() - 8192))
                tail = fh.read().decode(errors="replace")
        except OSError as e:
            return None, f"cannot read {logs[-1]} ({e.strerror or e})"
        if "ALL JOBS DONE" in tail:
            return True, f"{logs[-1]} ends ALL JOBS DONE"
        return False, (f"{logs[-1]} has no ALL JOBS DONE"
                       + (" (it was cancelled)" if "remaining tasks" in tail else "")
                       + ": FragPipe did not finish, whatever is in dia-quant-output/")
    return True, None


def find_report(out, engine):
    for c in ENGINE_REPORTS.get(engine, ()):
        p = os.path.join(out, c)
        if os.path.exists(p):
            return p
    return None


_FASTA_ARG = re.compile(r"--fasta\s+(\"[^\"]+\"|'[^']+'|\S+)")


def search_fastas(out):
    """The FASTA path(s) this search actually used: `fasta` in search_provenance.json, then every
    `--fasta` on the engine command line DIA-NN echoes at the top of its log. Read from the head
    of the log only (2 MB): a 399-file command line puts `--fasta` ~60 KB in, and the rest of a
    log is the search itself. Empty when neither says -- which is the caller's cue not to guess."""
    found = []
    try:
        with open(os.path.join(out, "search_provenance.json")) as fh:
            p = json.load(fh)
        if isinstance(p, dict) and p.get("fasta"):
            found.append(str(p["fasta"]))
    except (OSError, ValueError):
        pass
    for log in ("report.log.txt", "dia-quant-output/report.log.txt"):
        try:
            with open(os.path.join(out, log), errors="replace") as fh:
                head = fh.read(2 << 20)
        except OSError:
            continue
        found += [m.strip("\"'") for m in _FASTA_ARG.findall(head)]
    return tuple(dict.fromkeys(found))


def explicit_meta_mismatch(out, meta):
    """Why an explicit --fasta-meta does NOT describe this search's database, or None (it matches,
    or the search's FASTA is unknown so there is nothing to check against). The same tie as for a
    sidecar found nearby: an explicit one got no check at all, so a wrong --fasta-meta put another
    search's database (and organism) into the manifest. A mismatched one is dropped and the
    sidecar tied to the search's own FASTA is used instead, if there is one."""
    used = search_fastas(os.path.abspath(out))
    if not meta or not used:
        return None
    names = {os.path.basename(f) for f in used}
    reals = {os.path.realpath(f) for f in used}
    sib = meta[:-len(".meta.json")] if meta.endswith(".meta.json") else None
    if sib and (os.path.realpath(sib) in reals or os.path.basename(sib) in names):
        return None
    mf = _meta_fasta(meta)
    if mf and (mf in reals or os.path.basename(mf) in names):
        return None
    return (f"--fasta-meta {meta} describes {mf or sib or 'an unknown FASTA'}, but the search read "
            f"{', '.join(used)}")


def _meta_candidates(out, explicit_meta=None):
    """The <fasta>.meta.json files that describe THIS search's database, best first.

    Globbing the neighbourhood and taking the first hit is how a mouse search was handed to FRAN
    with the HUMAN database: PROT_0793/ holds human_, mouse_ and mouse_mousecont .fasta.meta.json
    side by side, `sorted()` puts human first, and search_mouse_mousecont's manifest recorded
    human_UP000005640.fasta (seen in the drop dir on 2026-09-24). So a neighbouring meta is only
    used when it can be tied to a FASTA the search really used (search_fastas), and when the
    search's FASTA is unknown, only a LONE meta is used -- several are a guess, and a guessed
    database or organism is a claim about the data (architectural rule #2)."""
    if explicit_meta and explicit_meta_mismatch(out, explicit_meta):
        explicit_meta = None       # the wrong database: dropped; the search's own sidecar decides
    cands = [explicit_meta] if explicit_meta else []
    used = search_fastas(os.path.abspath(out))
    cands += [f + ".meta.json" for f in used if os.path.isfile(f + ".meta.json")]
    parent = os.path.dirname(os.path.abspath(out))
    near = glob.glob(os.path.join(out, "*.meta.json")) \
        + glob.glob(os.path.join(parent, "*.fasta.meta.json")) \
        + glob.glob(os.path.join(parent, "input", "*.fasta.meta.json"))
    if os.path.basename(parent) == "output":     # a session: <session>/output/search
        near += glob.glob(os.path.join(os.path.dirname(parent), "input", "*.fasta.meta.json"))
    near = sorted(dict.fromkeys(os.path.realpath(p) for p in near))
    if used:
        names = {os.path.basename(f) for f in used}
        reals = {os.path.realpath(f) for f in used}
        for m in near:
            stem = os.path.basename(m)[:-len(".meta.json")]
            if stem in names or _meta_fasta(m) in reals:
                cands.append(m)
    elif len(near) == 1:
        cands += near
    return list(dict.fromkeys(cands))


def _meta_fasta(meta_path):
    """Real path of the FASTA a meta.json describes, or None."""
    try:
        with open(meta_path) as fh:
            m = json.load(fh)
    except (OSError, ValueError):
        return None
    sel = m.get("selected") if isinstance(m.get("selected"), dict) else m
    f = sel.get("fasta") or m.get("fasta")
    return os.path.realpath(f) if f else None


def organism_from_meta(out, explicit_meta=None):
    """(organism, taxid, source). DIA-NN's report carries NO organism column, so unless we pass
    it the corpus row is left NULL -- the species page then simply cannot see this search. The
    skill already had the user CONFIRM the organism at step 3, and fetch_fasta.py wrote it to
    <fasta>.meta.json, so pass that through rather than leaving a hole. Never guessed."""
    for c in _meta_candidates(out, explicit_meta):
        if not c or not os.path.isfile(c):
            continue
        try:
            with open(c) as fh:
                m = json.load(fh)
        except (OSError, ValueError):
            continue
        sel = m.get("selected") if isinstance(m.get("selected"), dict) else m
        org = sel.get("organism") or m.get("organism")
        tax = sel.get("taxid") or m.get("taxid")
        if org:
            return org, (int(tax) if tax else None), os.path.abspath(c)
    return None, None, None


def fasta_from_meta(out, explicit_meta=None):
    """Which DATABASE this search used, from the <fasta>.meta.json fetch_fasta.py wrote.

    The cron cannot derive this. FRAN's own detector (ingest/engine_fasta.py) falls back to
    parsing --fasta out of report.log.txt, which works for DIA-NN and is best-effort for
    Spectronaut -- but a search that came through this skill KNOWS its database exactly, so
    hand it over rather than making the corpus guess.

    Why the corpus wants it: fasta_n_proteins divided by the distinct gene count gives
    entries-per-gene, and that is what separates a real depth difference from database
    redundancy. A one-protein-per-gene proteome sits near 1.00; a full proteome with
    unreviewed isoforms can exceed 2, which alone can move a cross-engine protein-group gap
    from a few percent to tens of percent.

    Returns (path, md5, n_entries) with any unknown field None. Older meta files predate the
    md5/n_entries fields, so those are recomputed here when the FASTA is still readable --
    and left None rather than invented when it is not.
    """
    for c in _meta_candidates(out, explicit_meta):
        if not c or not os.path.isfile(c):
            continue
        try:
            with open(c) as fh:
                m = json.load(fh)
        except (OSError, ValueError):
            continue
        sel = m.get("selected") if isinstance(m.get("selected"), dict) else m
        path = sel.get("fasta") or m.get("fasta")
        if not path:
            continue
        # The PATH is the file the search READ, when the search says so; the sidecar describes the
        # content (md5, count) but records wherever fetch_fasta.py first wrote it. Two ways that
        # differs, both seen in the drop dir on 2026-09-24: hive_remote writes the sidecar on the
        # laptop (Dupanloup: a /Users/... path), and a session's sidecar can name the user's
        # staging copy (~/proteomics-pipeline/staging/search.fasta, overwritten by the next
        # session) while the search read <session>/input/search.fasta. Same md5 in all three.
        sib = os.path.realpath(c[:-len(".meta.json")]) if c.endswith(".meta.json") else None
        here = [f for f in search_fastas(os.path.abspath(out)) if os.path.isfile(f)
                and (os.path.realpath(f) == sib
                     or os.path.basename(f) == os.path.basename(path))]
        if here:
            path = here[0]
        md5 = sel.get("md5") or m.get("md5")
        # `n_sequences` is what fetch_fasta.py has ALWAYS written (proteome + appended
        # contaminants), and it is the same number: verified on a real sidecar, n_sequences 34306
        # against 34306 headers counted in the file. Reading it means no meta written before the
        # md5/n_entries fields has to be re-scanned at all -- which matters because this runs
        # inside `check`, on a login node, and the file is hundreds of MB (measured 0.6-1.3 s).
        n = (sel.get("n_entries") or m.get("n_entries")
             or sel.get("n_sequences") or m.get("n_sequences"))
        # Only a sidecar carrying NEITHER count falls through to counting the file, and only when
        # the FASTA is still there. Left None rather than invented when it is not: a guessed
        # database size is a claim about comparability.
        if (md5 is None or n is None) and os.path.isfile(path):
            helpers = _fasta_helpers()
            if helpers:
                _m, _c = helpers
                md5 = md5 if md5 is not None else _m(path)
                n = n if n is not None else _c(path)
        return path, md5, n
    return None, None, None


def _fasta_helpers():
    """(md5, entry-count) from fetch_fasta.py, or None if it cannot be imported.

    ONE definition, imported rather than copied -- the same rule radiant_parallel.py follows for
    slurm_queue(). The entry counter walks chunk boundaries (a per-chunk count(b"\n>") silently
    drops every header landing on one), and a second copy of logic like that is a bug fixed in one
    place and left in the other. Degrades to None rather than raising: the fields are then absent
    from the manifest, which is honest, and FRAN's own log-parsing detector is the fallback."""
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from fetch_fasta import _count_entries, _md5
        return _md5, _count_entries
    except Exception as e:                                          # noqa: BLE001
        sys.stderr.write(f"[fran_deposit] could not import the FASTA helpers from "
                         f"fetch_fasta.py ({type(e).__name__}: {e}); recording the database "
                         f"path without md5/entry count\n")
        return None


def entry_name(out):
    """Deterministic drop-entry name for a search dir: `<dir name>__<8 hex of its real path>`.

    Deterministic on purpose. If re-staging invented a new name each time, the cron would see a
    second candidate directory for a search it already ingested, and — since it keys idempotency on
    the path it ingested from — would write a duplicate search rather than replacing one. The hash
    is of the RESOLVED path so two directories with the same basename never collide."""
    real = os.path.realpath(out)
    base = re.sub(r"[^A-Za-z0-9._-]+", "_", os.path.basename(real.rstrip("/"))) or "search"
    return f"{base}__{hashlib.sha1(real.encode()).hexdigest()[:8]}"


def find_xic_dirs(out):
    """Every DIA-NN `--xic` chromatogram directory this search produced, whichever route ran it.

    There is no single location, and the difference is invisible unless you go looking:

      <out>/report_xic/                 single-shot DIA-NN (--out <out>/report.parquet)
      <out>/xic/t<N>_xic/               the 5-STEP PARALLEL CHAIN -- one dir per array task,
                                        because step 4 runs per file with --out <out>/xic/t<N>.parquet
      <out>/dia-quant-output/report_xic FragPipe's bundled DIA-NN

    DIA-NN names the directory `<--out stem>_xic`, so the parallel chain scatters chromatograms
    across as many directories as there were files -- 399 of them on a real 399-run cohort, with
    NOTHING at <out>/report_xic. A finder that only knew the single-shot layout reported "no XICs"
    for every parallel run, which is the default route above 5 files. Verified against real HIVE
    output: 399 tasks -> 399 *.xic.parquet, one per run, no basename collisions.
    """
    hits, seen = [], set()

    def add(d, rel):
        if d in seen or not os.path.isdir(d):
            return
        n = len([f for f in os.listdir(d) if f.endswith(".xic.parquet")])
        if n:
            seen.add(d)
            hits.append({"rel": rel, "dir": d, "n_files": n})

    for rel in ("report_xic", "dia-quant-output/report_xic"):
        add(os.path.join(out, rel), rel)
    xic_root = os.path.join(out, "xic")
    if os.path.isdir(xic_root):
        for name in sorted(os.listdir(xic_root)):
            if name.endswith("_xic"):
                add(os.path.join(xic_root, name), f"xic/{name}")
    # Anything else DIA-NN named `<stem>_xic` at the top level (a search run with a custom --out).
    for name in sorted(os.listdir(out)) if os.path.isdir(out) else []:
        if name.endswith("_xic"):
            add(os.path.join(out, name), name)
    return hits


def xic_files(out):
    """Every chromatogram/mobilogram file across all of this search's XIC directories, as
    (basename, full path). Run-named by DIA-NN and unique across array tasks (verified), so they
    can be presented as one flat directory."""
    files = []
    for d in find_xic_dirs(out):
        for f in sorted(os.listdir(d["dir"])):
            if f.endswith((".xic.parquet", "mobilogram.parquet")):
                files.append((f, os.path.join(d["dir"], f)))
    return files


def read_receipt(out):
    p = os.path.join(out, RECEIPT)
    if os.path.isfile(p):
        try:
            with open(p) as fh:
                return json.load(fh)
        except (OSError, ValueError):
            return {"malformed": True, "path": p}
    return None


def write_receipt(out, data):
    p = os.path.join(out, RECEIPT)
    try:
        with open(p, "w") as fh:
            json.dump(data, fh, indent=2)
        try:
            os.chmod(p, 0o664)        # the next Core member's stage/verify rewrites it
        except OSError:
            pass
        return p
    except OSError as e:
        # A receipt we could not write is a resume hazard, not a failure of the deposit -- report
        # it instead of raising, so the caller learns the deposit itself still stands.
        data["receipt_error"] = str(e)
        return None


# ==================================================================================== QC ==
# QC runs are not handed to FRAN (Brett, 2026-09-24): a HeLa QC series is instrument monitoring,
# not a customer search, and FRAN keeps QC out of the corpus for the same reason it excludes STAN.
# FRAN's scanner can only see QC by PATH (find_uningested.DEFAULT_EXCLUDES), and a drop entry in
# incoming/ has no QC path, so the decision has to be made here, before anything is staged.
#
# ONE definition: is_qc_run(). backfill and stage both call it; nothing else decides QC.
# TWIN RULE in FRAN: ingest/find_uningested.py `policy_exclusion` / `qc_reason` / `QC_NAME_RE`
# (branch fix/auto-ingest-starvation). The two must agree byte for byte -- change one, change the
# other. Pinned by the same vectors on both sides (tests/test_fran_health_backfill.py QcRuleTests):
#   excluded: "chkLUppm_HeLa50_2026 Lumos QC", "QC_run_01", "hela_qc_2", "Exploris QC2"
#   kept:     "HeLa_digest_timecourse", "aqc_buffer_study", "QCM_study", "Plasma_liver2"
# "HeLa" alone is deliberately NOT a QC signal: HeLa digests are real experiments too.
QC_NAME_RE = re.compile(r"(?i)(?<![a-z0-9])qc(?![a-z])")
# FRAN's find_uningested.DEFAULT_EXCLUDES, verbatim: trees whose engine output is never a corpus
# search (STAN QC, the QC watcher, smoke tests, FRAN's scratch dir, the ToF QC series). Substring
# match, as FRAN does it. A search under one is excluded by POLICY -- even with --not-qc, because
# FRAN's scanner would refuse it anyway. backfill's walk prunes the same list (BACKFILL_EXCLUDES).
FRAN_DEFAULT_EXCLUDES = ("/quobyte/proteomics-grp/STAN/", "/quobyte/proteomics-grp/hela_qcs/",
                         "/quobyte/proteomics-grp/brett/v1_smoke",
                         "/quobyte/proteomics-grp/brett/glendon/", "/Data/lab/ToFEvoQC/")
QC_OVERRIDE = "user override"          # the qc_rule a --not-qc stage writes beside "qc": false
# FRAN's _excludes_path (find_uningested.py, c4838fe): the laptop's SMB spelling of the group share
# is mapped to HIVE's before DEFAULT_EXCLUDES is tested, and a trailing "/" lets a root itself match.
_SMB_PREFIX, _HIVE_PREFIX = "/Volumes/proteomics-grp", "/quobyte/proteomics-grp"


def _excludes_path(path):
    s = str(path or "").replace("\\", "/").rstrip("/")
    if s == _SMB_PREFIX or s.startswith(_SMB_PREFIX + "/"):
        s = _HIVE_PREFIX + s[len(_SMB_PREFIX):]
    return s + "/"
# Session metadata that can carry an explicit `"qc": true|false` for a whole session.
SESSION_QC_FILES = ("session.json", os.path.join("input", "session.json"),
                    os.path.join("input", "wf", "workflow.manifest.json"))


def session_for(out):
    """The session directory a search out dir belongs to, or None: `<session>/output/search`
    (session.py's layout), or a job dir holding the search next to its `input/` (hive_remote)."""
    out = os.path.abspath(out)
    parent = os.path.dirname(out)
    cands = [os.path.dirname(parent)] if os.path.basename(parent) == "output" else []
    cands.append(parent)
    for c in cands:
        if os.path.isfile(os.path.join(c, "README.md")) or os.path.isdir(os.path.join(c, "input")):
            return c
    return None


def _session_title(session):
    """The name the session was created with: README.md's `# <name>` (session.py init)."""
    try:
        with open(os.path.join(session, "README.md"), errors="replace") as fh:
            first = fh.readline().strip()
        return first[2:].strip() if first.startswith("# ") else None
    except OSError:
        return None


def _session_qc_marker(session):
    """(True|False, file) when session metadata says explicitly, else (None, None)."""
    for rel in SESSION_QC_FILES:
        p = os.path.join(session, rel)
        try:
            with open(p) as fh:
                v = json.load(fh).get("qc")
        except (OSError, ValueError, AttributeError):
            continue
        if isinstance(v, bool):
            return v, p
    return None, None


def is_qc_run(out, session=None, *, names=(), override=None):
    """(is_qc, why) for a search out dir. FRAN's precedence (policy_exclusion), first match wins:
      1. an explicit QC marker -- `override=True` (stage --qc) or session metadata "qc": true
      2. the out dir is under one of FRAN's DEFAULT_EXCLUDES trees -- excluded EVEN with --not-qc:
         FRAN refuses anything there, so staging it would only put a refusal in its queue
      3. an explicit NOT-QC marker -- `override=False` (--not-qc) or session metadata "qc": false
      4. QC_NAME_RE on the search name(s) in `names` and the session name, then on the last three
         components of the out dir path
    `why` mirrors FRAN's reason text: "QC run: excluded by policy (<what matched>)". It goes into
    the receipt, and for a staged search into the manifest's `qc_rule`."""
    out = os.path.abspath(out)
    # FRAN judges the manifest's output_dir, which is realpath(out): the path components the name
    # rule reads must be the REAL ones, or a symlinked search dir escapes the rule here and is
    # excluded there. DEFAULT_EXCLUDES is tested on both spellings.
    real = os.path.realpath(out)
    session = session or session_for(out)
    marker, src = _session_qc_marker(session) if session else (None, None)
    if override is True:
        return True, "QC run: excluded by policy (qc: true, user override)"
    if marker is True:
        return True, f"QC run: excluded by policy (session metadata qc: true, {src})"
    for p in dict.fromkeys((out, real)):
        norm = _excludes_path(p)
        for root in FRAN_DEFAULT_EXCLUDES:
            if root in norm:
                return True, f"QC run: excluded by policy (output_dir is under {root} (DEFAULT_EXCLUDES))"
    if override is False:
        return False, QC_OVERRIDE
    if marker is False:
        return False, f"not QC: session metadata qc: false ({src})"
    labelled = [("search_name", n) for n in ([names] if isinstance(names, str) else names) if n]
    if session:
        labelled += [("session_name", n) for n in (_session_title(session),
                                                   os.path.basename(session)) if n]
    labelled += [("output_dir", c) for c in [x for x in real.split("/") if x][-3:]]
    for field, text in labelled:
        if QC_NAME_RE.search(text):
            return True, f"QC run: excluded by policy ({field} {text!r} matches QC_NAME_RE)"
    return False, "not QC: no qc marker, not under a FRAN excluded tree, no QC token in the names or path"


def _recorded_qc(receipt, manifest):
    """The QC decision an earlier stage left behind: (True, why) | (False, why) | (None, None).

    A withdrawal must STAY a withdrawal. A later plain `stage --out X` (no --name, no flag) used to
    re-stage a withdrawn QC run with qc: false -- which FRAN honours -- so the run got ingested
    after all. A qc_run receipt, or a staged manifest saying qc/exclude: true, is therefore QC
    until someone says --not-qc explicitly. A recorded --not-qc ("user override") is not-QC."""
    if receipt.get("status") == "qc_run":
        return True, (f"QC run: excluded by policy (recorded by {receipt.get('decided_by') or '?'} "
                      f"at {receipt.get('at') or '?'}: {receipt.get('qc_rule') or 'qc_run'})")
    if manifest.get("qc") is True or manifest.get("exclude") is True:
        return True, (f"QC run: excluded by policy (its staged manifest says qc: true: "
                      f"{manifest.get('qc_rule') or 'no reason recorded'})")
    if QC_OVERRIDE in (receipt.get("qc_rule"), manifest.get("qc_rule")):
        return False, QC_OVERRIDE
    return None, None


def decide_qc(out, *, names=(), override=None, receipt=None, manifest=None):
    """THE QC decision for one search, for stage and backfill alike: an explicit flag, else what
    an earlier stage recorded, else is_qc_run's rule. (is_qc_run still puts FRAN's excluded trees
    above any not-QC decision.)"""
    if override is None:
        rec, why = _recorded_qc(receipt or {}, manifest or {})
        if rec is True:
            return True, why
        override = rec
    return is_qc_run(out, names=names, override=override)


def _write_json_group(path, data):
    """Write JSON atomically and leave it group-writable (0664), whatever the umask. Anything the
    skill writes under incoming/ must be rewritable by the NEXT Core member: with a 022 umask a
    manifest came out 0644, and another member's withdrawal then failed."""
    tmp = f"{path}.{os.getpid()}.tmp"
    with open(tmp, "w") as fh:
        json.dump(data, fh, indent=2)
    os.chmod(tmp, 0o664)
    os.replace(tmp, path)


def _withdraw_qc_entry(entry, why, user):
    """A QC run that was staged BEFORE it was known to be QC -- the job-end hook stages without
    the analysis name, and "Lumos QC" may only arrive with the agent's --name later -- is marked
    rather than deleted: the manifest gets `"qc": true, "exclude": true`, which FRAN's ingester
    honours. Nothing in the shared drop dir is removed. Returns (entry, None) when the entry is
    marked, or (None, why-not) -- a failure is reported as one, never as QC kept out."""
    if os.path.islink(entry):
        return None, (f"{entry} is a symlink, not a drop entry with a manifest: nothing is written "
                      f"through it (re-stage it to get a real entry)")
    mp = os.path.join(entry, MANIFEST)
    try:
        with open(mp) as fh:
            man = json.load(fh)
    except (OSError, ValueError) as e:
        return None, f"cannot read {mp} ({type(e).__name__}: {getattr(e, 'strerror', None) or e})"
    if man.get("qc") is True and man.get("exclude") is True:
        return entry, None
    man.update(qc=True, exclude=True, qc_rule=why, withdrawn_by=user,
               withdrawn_at=_utc_now())
    try:
        _write_json_group(mp, man)
    except OSError as e:
        return None, f"cannot rewrite {mp} ({e.strerror or e})"
    return entry, None


def _utc_now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def check(a):
    out = os.path.abspath(a.out)
    r = {"eligible": False, "reason": None, "search_dir": out, "fran": FRAN_URL,
         "user": getpass.getuser()}

    if a.skip or os.environ.get("FRAN_DEPOSIT", "").lower() in ("off", "0", "no", "false"):
        r.update(reason="opted_out", detail="FRAN_DEPOSIT=off (or --skip)")
        return r
    if not os.path.isdir(out):
        r.update(reason="not_on_hive",
                 detail=f"{out} does not exist here. Run this ON HIVE (hive_exec.sh) with the "
                        f"HIVE-side search directory.")
        return r
    # QC before anything else that could make it look eligible: a QC run is never staged, and one
    # staged before it was known to be QC is withdrawn by stage() (see _withdraw_qc_entry).
    override = (True if getattr(a, "qc", False) else
                False if getattr(a, "not_qc", False) else None)
    # A decision an EARLIER stage recorded outlives it (decide_qc): the job-end hook stages with
    # the --qc / --not-qc baked in at generation time, and neither the agent's later flagless stage
    # nor a plain `stage --out X` may quietly overturn it -- least of all a withdrawal.
    prior_man = {}
    try:
        with open(os.path.join(os.environ.get("FRAN_DROP_DIR", DROP_DIR), entry_name(out),
                               MANIFEST)) as fh:
            prior_man = json.load(fh)
    except (OSError, ValueError):
        pass
    name = (a.name or "").strip() or None
    qc, why = decide_qc(out, names=[name] if name else (), override=override,
                        receipt=read_receipt(out) or {}, manifest=prior_man)
    r["qc"], r["qc_rule"] = qc, why
    if qc:
        r.update(reason="qc_run", name=name,
                 entry=os.path.join(os.environ.get("FRAN_DROP_DIR", DROP_DIR), entry_name(out)),
                 detail=f"QC run ({why}): QC runs are not handed to FRAN. If this is NOT a QC "
                        f"run, stage it with --not-qc.")
        return r
    # Coursework is never Core work: ONE rule (_teaching_reason) for stage and backfill alike, so
    # the job-end hook of a teaching account (proteomics-class-NN, members of proteomics-grp and so
    # able to write the drop dir) never stages a class exercise.
    teach = _teaching_reason(out, r["user"])
    if teach:
        r.update(reason="not_core_facility", core_member=False, detail=teach)
        return r
    # The Core gate. Not a flag: the write credential lives under this directory, so a HIVE user
    # outside proteomics-grp physically cannot deposit. Collaborators' searches stay theirs.
    r["proteomics_grp_access"] = os.path.isdir(GROUP_ROOT) and os.access(GROUP_ROOT, os.R_OK)

    engine, engine_version, engine_src = detect_engine(out)
    r.update(engine=engine, engine_version=engine_version, engine_source=engine_src)
    # Carried into the manifest rather than linked into the entry -- see LINK_ITEMS on why linking
    # this file would relabel the search as Radiant.
    try:
        with open(os.path.join(out, "search_provenance.json")) as fh:
            r["search_provenance"] = json.load(fh)
    except (OSError, ValueError):
        r["search_provenance"] = None
    if engine in UNSUPPORTED:
        r.update(reason="engine_unsupported", detail=UNSUPPORTED[engine])
        return r
    if engine not in ENGINE_REPORTS:
        r.update(reason="engine_unsupported",
                 detail=f"engine '{engine}' has no FRAN corpus adapter ({engine_src})")
        return r

    report = find_report(out, engine)
    # Deliberately a stat, not a parse: `check` must stay login-node safe, and reading a
    # multi-GB report to count rows is precisely the compute that belongs in the job below.
    r["report"] = report
    r["report_bytes"] = dir_size(report) if report else 0
    if not report or r["report_bytes"] == 0:
        r.update(reason="search_incomplete",
                 detail=f"no non-empty report in {out}. A failed or zero-ID search must never be "
                        f"deposited -- fix the search first (references/watcher.md).")
        return r
    # A non-empty report is not a finished search for FragPipe or Radiant (completion_marker).
    done, evidence = completion_marker(out, engine, report)
    r["completion"] = {"done": done, "evidence": evidence}
    if done is False:
        r.update(reason="search_incomplete",
                 detail=f"{evidence}. A partial search must never be deposited.")
        return r
    # No marker to read. The agent stages only after watching the job to COMPLETED (step 7c), so
    # its stage goes ahead; an unattended backfill has no such witness and does not guess.
    if done is None and getattr(a, "require_completion_marker", False):
        r.update(reason="needs_agent_check",
                 detail=f"{evidence}. backfill cannot tell a finished search from a partial one "
                        f"here: an agent confirms it finished, then runs stage --out {out}")
        return r

    # Checked BEFORE the corpus environment, and deliberately: "this search is already in" is
    # both cheap and decisive, and a resumed session must learn it even on a box where the token or
    # the ingester is momentarily unreadable. (Found by tests/test_fran_deposit_gate.py, which got
    # `ingest_code_missing` for a search whose receipt was sitting right there.)
    prior = read_receipt(out)
    # An ALLOWLIST, not a denylist: only a receipt claiming the search is staged or already
    # ingested may stand in the way. Every other status means FRAN does NOT have it, and blocking
    # on those would drop the search from the corpus forever -- the one receipt proving it is
    # missing would be the thing preventing the hand-over that fixes it.
    if prior and prior.get("status") not in BLOCKING_STATUS:
        prior = None
    if prior and not a.force:
        r.update(reason="already_staged", prior=prior,
                 detail=f"{os.path.join(out, RECEIPT)} already records this search as "
                        f"'{prior.get('status')}'. Pass --force to stage it again (the drop entry "
                        f"name is deterministic, so re-staging replaces it rather than adding a "
                        f"second candidate).")
        return r
    r["prior"] = prior

    drop = os.environ.get("FRAN_DROP_DIR", DROP_DIR)
    r["drop_dir"] = drop
    r["entry"] = os.path.join(drop, entry_name(out))
    # Create it on demand: the first Core search to finish should not fail because nobody had made
    # the directory yet. 2775 = group-writable + setgid, so every later Core member can stage into
    # it and the group ownership sticks.
    gate = drop
    if not os.path.isdir(drop) and getattr(a, "dry_run", False):
        # A dry run (`stage --dry-run`, `backfill`) creates nothing in a shared tree. Whoever could
        # create the directory is whoever could write it, so the gate moves to its parent.
        r["drop_dir_missing"] = True
        gate = os.path.dirname(drop.rstrip("/")) or drop
    elif not os.path.isdir(drop):
        try:
            os.makedirs(drop, mode=0o2775, exist_ok=True)
            os.chmod(drop, 0o2775)        # makedirs' mode is filtered by the umask
        except OSError as e:
            r.update(reason="no_drop_dir",
                     detail=f"FRAN's drop directory {drop} does not exist and could not be created "
                            f"({e}). Create it once: mkdir -p {drop} && chmod 2775 {drop}")
            return r
    # THE CORE GATE, and it is a real permission rather than a flag: the drop directory lives
    # inside /quobyte/proteomics-grp, so a HIVE account outside proteomics-grp cannot write here.
    # A collaborator's search physically cannot be staged.
    if not os.access(gate, os.W_OK | os.X_OK):
        # Not writable means "not a Core member" ONLY when the account really is outside the group.
        # A member who cannot write it has hit a permission problem -- reported, never recorded as
        # a decision, or the search would be locked out of FRAN by a transient fault.
        member = _in_core_group()
        r["core_member"] = member
        if member:
            r.update(reason="drop_dir_not_writable",
                     detail=f"{r['user']} is in {GROUP_NAME} but cannot write {gate}: a permission "
                            f"problem, not a decision. Fix it (chmod 2775 {gate}) and stage again.")
        else:
            r.update(reason="not_core_facility",
                     detail=f"{drop} is not writable by {r['user']}, so this is not a Proteomics "
                            f"Core run. Collaborator searches are never handed to the Core corpus.")
        return r

    bad_meta = explicit_meta_mismatch(out, a.fasta_meta)
    org, tax, org_src = organism_from_meta(out, a.fasta_meta)
    if bad_meta:
        r["fasta_meta_ignored"] = bad_meta
        sys.stderr.write(f"[fran_deposit] WARNING: --fasta-meta ignored ({bad_meta}); "
                         + (f"using the sidecar of the FASTA the search read ({org_src})"
                            if org_src else "no sidecar of the FASTA the search read, so the "
                                            "organism and database are left blank")
                         + "\n")
    # None, never "": FRAN's read_manifest rejects an empty search_name or organism outright
    r["organism"] = (a.organism or "").strip() or org or None
    fp, fmd5, fn = fasta_from_meta(out, a.fasta_meta)
    r["fasta_path"], r["fasta_md5"], r["fasta_n_proteins"] = fp, fmd5, fn
    r["taxon"] = a.taxon or tax
    r["organism_source"] = "--organism (given)" if a.organism else org_src
    r["name"] = (a.name or "").strip() or None
    # Where a name came from when nobody typed it: `backfill` derives one from the folder, and
    # the manifest says so rather than passing a folder name off as the analysis name.
    if r["name"] and getattr(a, "name_source", None):
        r["name_source"] = a.name_source

    # XIC chromatograms, if the search asked DIA-NN for them (`--xic` in the cfg; diann_parallel.py
    # puts it on step 4 only). Reported either way -- their absence is a fact about the search, not
    # an error, and FRAN's XIC lane simply has nothing to ingest for this one.
    xics = find_xic_dirs(out)
    r["xic"] = {"present": bool(xics), "n_files": sum(x["n_files"] for x in xics), "dirs": xics}
    r["links"] = [i for i in LINK_ITEMS if os.path.exists(os.path.join(out, i))]

    r.update(eligible=True, reason="ok")
    r["next_command"] = f"python3 {os.path.abspath(__file__)} stage --out {shlex.quote(out)}"
    return r


def stage_argv(out, *, name=None, qc=None, fasta_meta=None, python=None):
    """THE argv for `fran_deposit.py stage` -- for the job-end hook, which bakes it into the job at
    GENERATION time. That is what closes the QC race: a hook that stages without the analysis name
    or the user's QC decision hands a QC run to FRAN, and the cron (every 4 h) can ingest it before
    the agent's later `stage --name "... Lumos QC"` withdraws it.
      name  the session's descriptive name (the corpus name, and what the QC name rule reads)
      qc    True -> --qc, False -> --not-qc, None -> let the rule decide
    Shell-quote each element when baking it into a script (names have spaces and em dashes)."""
    argv = [python or sys.executable or "python3", os.path.abspath(__file__), "stage", "--out", out]
    if fasta_meta:
        argv += ["--fasta-meta", fasta_meta]
    if name and str(name).strip():
        argv += ["--name", str(name).strip()]
    if qc is True:
        argv.append("--qc")
    elif qc is False:
        argv.append("--not-qc")
    return argv


def stage(a):
    """Hand the search to FRAN by linking it into the drop directory. Nothing is copied, nothing is
    ingested here — FRAN's cron does the ingest when it next scans.

    Prints ONE JSON object (the contract the orchestrator and the job-end hook parse). When the
    last `health` run found FRAN's ingest unhealthy it also carries `health_warning` and prints that
    line to stderr -- the search IS handed over either way, and a health problem never changes the
    exit status. stage runs inside every search job, so it only READS health's status file: no
    network, no database, no log scan (attach_health_status)."""
    res = do_stage(a)
    if res.get("staged") or res.get("reason") == "already_staged":
        attach_health_status(res)
    jout(res, 0)


# Receipt statuses stage() records for a DECISION not to hand a search over, so a later `backfill`
# (run by someone else, months on) honours it instead of staging the search anyway. None of them
# blocks an explicit `stage`: BLOCKING_STATUS is an allowlist and these are not on it.
DECISION_STATUS = ("opted_out", "not_core_facility", "qc_run")


def _record_decision(a, c):
    out = c.get("search_dir")
    if (c.get("reason") not in DECISION_STATUS or getattr(a, "dry_run", False)
            or not out or not os.path.isdir(out)):
        return
    # "Not a Core member" is recorded only when the account is KNOWN to be outside the group (or is
    # a teaching account). When membership cannot be read, the refusal stands for this call only.
    if c["reason"] == "not_core_facility" and c.get("core_member") is not False:
        return
    prior = read_receipt(out) or {}
    rec = {"status": c["reason"], "search_dir": out, "decided_by": c.get("user"),
           "at": datetime.datetime.now().isoformat(timespec="seconds"), "detail": c.get("detail")}
    if c["reason"] == "qc_run":
        rec["qc_rule"] = c.get("qc_rule")
        rec["search_name"] = c.get("name") or prior.get("search_name")
        # Staged before anyone knew it was QC (the job-end hook stages without the analysis name):
        # withdraw it by marking the manifest, which FRAN's ingester honours.
        if c.get("entry") and os.path.isdir(c["entry"]):
            done, err = _withdraw_qc_entry(c["entry"], c.get("qc_rule"), c.get("user"))
            if done:
                c["withdrawn"] = rec["withdrawn_entry"] = c["entry"]
                c["detail"] += (f" It had already been staged; {c['entry']}/{MANIFEST} is now "
                                f"marked qc: true, so FRAN's ingester skips it.")
            else:
                c["withdraw_failed"] = rec["withdraw_failed"] = err
                c["detail"] = (f"withdraw FAILED: {err}. This QC run ({c.get('qc_rule')}) is "
                               f"STILL STAGED at {c['entry']} and FRAN may ingest it. A Core member "
                               f"who can write {c['entry']}/{MANIFEST} must run: python3 "
                               f"{os.path.abspath(__file__)} stage --out {shlex.quote(out)} --qc")
        if prior.get("status") == "ingested":
            c["detail"] += " It is ALREADY in FRAN's corpus; taking it out is a FRAN-side step."
            return                   # never overwrite the record that it IS in FRAN
    if (prior.get("status") in BLOCKING_STATUS and not c.get("withdrawn")
            and not (c["reason"] == "qc_run" and c.get("withdraw_failed"))):
        return                       # never overwrite the record that it IS in FRAN
    write_receipt(out, rec)


def do_stage(a):
    """stage() without the printing: returns the result dict. `backfill --apply` calls this per
    search, so one search failing cannot take the rest of the batch down with a SystemExit."""
    c = check(a)
    if not c["eligible"]:
        _record_decision(a, c)
        return {**c, "staged": False}
    out, entry = c["search_dir"], c["entry"]
    if a.dry_run:
        return {**c, "staged": False, "dry_run": True,
                "would_link": c["links"], "would_write": os.path.join(entry, MANIFEST)}

    # WHEN it was first handed over, which FRAN's runner uses to take drop entries oldest-first. A
    # re-stage keeps the ORIGINAL time -- otherwise re-staging would push a search to the back of
    # the queue it has been waiting in. An entry staged before this field existed gets the time its
    # manifest was written, which is the same fact.
    staged_at = staged_by = None
    mp = os.path.join(entry, MANIFEST)
    if not os.path.islink(entry) and os.path.isfile(mp):
        try:
            with open(mp) as fh:
                prior_man = json.load(fh)
            staged_at, staged_by = prior_man.get("staged_at"), prior_man.get("staged_by")
        except (OSError, ValueError):
            pass
        if not staged_at:
            try:
                staged_at = datetime.datetime.fromtimestamp(
                    os.path.getmtime(mp), datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
            except OSError:
                pass
    restaged = bool(staged_at)
    staged_at = staged_at or _utc_now()
    staged_by = staged_by or c["user"]

    # Re-staging must converge, not accumulate: an entry from an earlier attempt may hold links to
    # files the search has since replaced (a resumed chain rewrites report.parquet). Relink from
    # scratch rather than leaving a stale mixture of both runs.
    # Everything under incoming/ is made group-writable EXPLICITLY (2775 dirs, 0664 files):
    # makedirs' mode and open() are filtered by the umask, and with a 022 umask the next Core
    # member could neither withdraw nor re-stage this entry -- the relink below died in os.unlink.
    try:
        # An entry that is itself a SYMLINK (a legacy bare-link entry, or one pointing at the search
        # dir) must never be followed: clearing it "relinked" report.parquet and report.log.txt
        # INSIDE THE REAL SEARCH DIR -- unlinked, replaced by self-referential links -- and the
        # chmod below then set the search dir to 2775. Remove the link itself; build a real entry.
        if os.path.islink(entry):
            os.unlink(entry)
        if os.path.isdir(entry):
            for f in os.listdir(entry):
                fp = os.path.join(entry, f)
                if os.path.islink(fp) or os.path.isfile(fp):
                    os.unlink(fp)
                elif f == "report_xic" and os.path.isdir(fp):
                    # A real directory of links we built. Clear it too: a re-run that dropped a
                    # pathological file would otherwise leave that file's trace behind forever.
                    for g in os.listdir(fp):
                        os.unlink(os.path.join(fp, g))
                    os.rmdir(fp)
        else:
            os.makedirs(entry, mode=0o2775, exist_ok=True)
        try:
            os.chmod(entry, 0o2775)
        except OSError:
            pass                      # another member's entry: its owner already set it
    except OSError as e:
        return {**c, "eligible": False, "staged": False, "reason": "entry_not_writable",
                "detail": f"cannot rewrite {entry} ({e.strerror or e}): it was staged by another "
                          f"account without group write. A permission problem, not a decision -- "
                          f"its owner can run: chmod -R g+w {entry}"}

    linked = []
    # Normalise the chromatograms into ONE `report_xic/` directory of links, whatever layout the
    # search used. FRAN's diann_xic_to_lance.py defaults to `<dir>/report_xic`, so the parallel
    # chain's 399 scattered `xic/t<N>_xic/` directories would need special handling at the other
    # end; flattening them here means the ingest side sees the same shape from every route. Safe
    # because DIA-NN names each file after its RUN, so there are no collisions across tasks.
    xf = [x for x in xic_files(out) if not os.path.exists(os.path.join(entry, "report_xic", x[0]))] \
        if c["xic"]["present"] else []
    if xf:
        xdir = os.path.join(entry, "report_xic")
        os.makedirs(xdir, exist_ok=True)
        os.chmod(xdir, 0o2775)
        for name, src in xf:
            dst = os.path.join(xdir, name)
            if not os.path.lexists(dst):
                try:
                    os.symlink(src, dst)
                except OSError as e:
                    c.setdefault("link_errors", []).append(f"report_xic/{name}: {e}")
        linked.append(f"report_xic/ ({len(xf)} files)")

    for item in c["links"]:
        src, dst = os.path.join(out, item), os.path.join(entry, item)
        try:
            os.symlink(src, dst)
            linked.append(item)
        except OSError as e:
            # One link failing must not abandon the rest: a staged entry missing its log is still
            # ingestable, and the manifest records exactly what did and did not make it.
            c.setdefault("link_errors", []).append(f"{item}: {e}")

    manifest = {
        "fran_manifest_version": 1,
        # THE AUTHORITATIVE FACTS. Read these rather than inferring from the entry: `output_dir` is
        # the real search directory (use it as corpus_ingest --output-dir so the corpus records
        # where the search actually lives and stays idempotent across re-stages), and `engine` is
        # what genuinely ran — marker sniffing cannot tell FragPipe's bundled DIA-NN from DIA-NN.
        "output_dir": os.path.realpath(out),
        "engine": c["engine"],
        "engine_version": c.get("engine_version"),
        "report": c["report"],
        "search_name": c.get("name"),
        **({"search_name_source": c["name_source"]} if c.get("name_source") else {}),
        # Absent, not guessed, when unknown -- a NULL organism is honest, an invented one is a
        # species claim. DIA-NN reports carry no organism column, so this is the only source.
        "organism": c.get("organism"),
        "taxon": c.get("taxon"),
        "organism_source": c.get("organism_source"),
        # The search DATABASE, for delimp_searches.fasta_path / fasta_md5 / fasta_n_proteins.
        # Absent, not guessed, when the meta.json is missing -- FRAN's own log-parsing detector
        # is the fallback for searches that did not come through this skill.
        "fasta_path": c.get("fasta_path"),
        "fasta_md5": c.get("fasta_md5"),
        "fasta_n_proteins": c.get("fasta_n_proteins"),
        "xic": c["xic"],
        "linked": linked,
        "staged_by": staged_by,
        # ISO-8601 UTC, the ORIGINAL hand-over time (kept across re-stages); FRAN orders its drop
        # box by this instead of the entry directory's mtime.
        "staged_at": staged_at,
        **({"restaged_at": _utc_now(), "restaged_by": c["user"]} if restaged else {}),
        # The QC decision, made here because FRAN cannot make it for a drop entry (is_qc_run).
        # Written on every staged search so the ingester can see a decision WAS made; `qc: false`
        # is authoritative there, and --not-qc records qc_rule "user override".
        "qc": False,
        "qc_rule": c.get("qc_rule"),
        "staged_by_pipeline": "ucdavis-proteomics-core-pipeline/fran_deposit.py",
        "search_provenance": c.get("search_provenance"),
        # TWO forms of the same thing, on purpose. `suggested_ingest` is an argv ARRAY, correct for
        # subprocess. `_shell` is the same list already shell-quoted, because the natural-looking
        # " ".join(args) breaks on every two-word species name -- "Homo sapiens" arrives as two
        # arguments and corpus_ingest exits "unrecognized arguments: sapiens". Every organism name
        # has a space in it, so that bug would fire on the first real ingest, not on an edge case.
        "suggested_ingest": suggested_ingest(c),
        "suggested_ingest_shell": " ".join(shlex.quote(x) for x in suggested_ingest(c)),
    }
    _write_json_group(os.path.join(entry, MANIFEST), manifest)

    receipt = {"status": "staged", "search_dir": out, "entry": entry, "engine": c["engine"],
               "organism": c.get("organism"), "taxon": c.get("taxon"),
               "search_name": c.get("name"), "staged_at": staged_at, "qc_rule": c.get("qc_rule"),
               "xic": c["xic"], "linked": linked, "staged_by": c["user"], "fran": FRAN_URL}
    if c.get("link_errors"):
        receipt["link_errors"] = c["link_errors"]
    write_receipt(out, receipt)
    return {**c, "staged": True, "entry": entry, "linked": linked,
            "manifest": os.path.join(entry, MANIFEST),
            "receipt": os.path.join(out, RECEIPT),
            "note": "FRAN's ingest cron picks this up on its next scan; nothing was copied.",
            "then": f"python3 {os.path.abspath(__file__)} verify --out {shlex.quote(out)}"}


def suggested_ingest(c):
    """The exact corpus_ingest.py arguments this search should be ingested with — so the cron does
    not have to re-derive the organism (it cannot: DIA-NN reports have no organism column) or guess
    the engine from markers. Advisory: the cron owns ingestion, this only removes the guesswork.

    Returned as an argv LIST. Pass it to subprocess as a list, or use the manifest's
    `suggested_ingest_shell` if it has to go through a shell — never " ".join() this."""
    args = ["--engine", c["engine"], "--output-dir", os.path.realpath(c["search_dir"])]
    if c.get("name"):
        args += ["--name", c["name"]]
    if c.get("organism"):
        args += ["--organism-name", c["organism"]]
    if c.get("taxon"):
        args += ["--taxon", str(c["taxon"])]
    return args


# Asked of the corpus itself, over the ingester's own connection helper, so the answer comes from
# the database rather than from the job having exited 0. A SLURM COMPLETED does not mean rows landed.
# Asked of the corpus itself, over the ingester's own connection helper, so the answer comes from
# the database rather than from a marker file. Only usable by an account that can read a PG Farm
# token -- staging deliberately needs none, so this half of verify is best-effort.
VERIFY_PY = r'''
import json, os, sys
sys.path.insert(0, os.environ["ING"])
from corpus_ingest import _conn
out = os.environ["OUTDIR"]
c = _conn(); cur = c.cursor()
cur.execute("""SELECT id, search_name, search_engine, search_engine_version, n_raw_files,
                      n_precursors_total, n_proteins_total, status, submitted_at
               FROM delimp_searches WHERE output_dir IN (%s, %s)
               ORDER BY submitted_at DESC""", (out, os.environ.get("ENTRY") or out))
rows = cur.fetchall()
res = {"in_corpus": bool(rows), "n_rows": len(rows)}
if rows:
    r = rows[0]
    res["search"] = {"search_id": r[0], "name": r[1], "engine": r[2], "engine_version": r[3],
                     "n_raw_files": r[4], "n_precursors": r[5], "n_proteins": r[6],
                     "status": r[7], "submitted_at": str(r[8])}
    cur.execute("SELECT COUNT(*) FROM delimp_precursors WHERE search_id=%s", (r[0],))
    res["precursor_rows"] = cur.fetchone()[0]
    cur.execute("SELECT COUNT(*) FROM delimp_proteins WHERE search_id=%s", (r[0],))
    res["protein_rows"] = cur.fetchone()[0]
c.close()
print(json.dumps(res))
'''

# What the cron drops into an entry when it has ingested it. The skill does not write this and does
# not require it -- an absent marker means "not ingested yet", which is the normal state for the
# first minutes after staging and must never be reported as a failure.
INGESTED_MARKER = "fran_ingested.json"


def corpus_query(out, entry=None):
    """Ask the corpus whether this search is in it. Returns None when this account has no way to
    ask -- which is the COMMON case now that staging needs no credential, and is not an error.

    Asked under BOTH names the search can have in the corpus: its real directory, and the drop
    entry. FRAN's auto_ingest passes realpath(<the dir it scanned>) as --output-dir, and a drop
    entry is a real directory, so a staged search is recorded under incoming/<entry> -- a lookup
    by the real search dir alone would never find one."""
    ing = first_readable(INGEST_DIRS, isdir=True)
    py = first_readable(PY_CANDIDATES, executable=True)
    tok = first_readable(TOKEN_CANDIDATES)
    if not (ing and py and tok):
        return None
    env = {**os.environ, "ING": ing, "OUTDIR": os.path.realpath(out), "DELIMP_PG_TOKEN_FILE": tok,
           "ENTRY": entry or ""}
    try:
        p = subprocess.run([py, "-c", VERIFY_PY], capture_output=True, text=True,
                           env=env, timeout=600)
        if p.returncode != 0:
            return {"error": (p.stderr or p.stdout).strip()[-400:]}
        return json.loads(p.stdout.strip().splitlines()[-1])
    except (OSError, subprocess.SubprocessError, ValueError, IndexError) as e:
        return {"error": str(e)[:200]}


def verify(a):
    """Did the handover actually land? Four questions, in order of how much they prove:
      1. is the entry still in the drop directory, with a live link to a real report?
      2. has the cron marked it ingested?
      3. (only if this account can read a corpus token) does the corpus actually hold the rows?
      4. otherwise, what do the cron's own logs say about this entry -- ingested, failed, or never
         reached? That needs no credential, so it is the answer most accounts get.
    States: ingested | ingest_failed | staged_pending_cron | not_staged. A staged-but-not-yet-
    ingested search is the NORMAL state right after a run — report it as pending, never as a
    failure; when the cron itself is stuck, the detail says so."""
    out = os.path.abspath(a.out)
    entry = os.path.join(os.environ.get("FRAN_DROP_DIR", DROP_DIR), entry_name(out))
    r = {"search_dir": out, "entry": entry, "staged": os.path.isdir(entry), "fran": FRAN_URL}

    if r["staged"]:
        # A link whose target has been deleted or moved still LOOKS like a staged file: os.listdir
        # shows it, and only following it reveals the entry is hollow. The cron would find the
        # directory, detect no engine, and silently skip it.
        broken = [f for f in os.listdir(entry)
                  if os.path.islink(os.path.join(entry, f))
                  and not os.path.exists(os.path.join(entry, f))]
        r["broken_links"] = broken
        r["links"] = sorted(f for f in os.listdir(entry) if f != MANIFEST)
        try:
            with open(os.path.join(entry, MANIFEST)) as fh:
                r["qc_excluded"] = json.load(fh).get("qc") is True
        except (OSError, ValueError):
            pass
        mk = os.path.join(entry, INGESTED_MARKER)
        if os.path.isfile(mk):
            try:
                with open(mk) as fh:
                    r["cron_marker"] = json.load(fh)
            except (OSError, ValueError):
                r["cron_marker"] = {"malformed": True}

    corpus = corpus_query(out, entry)
    r["corpus"] = corpus
    # The DB-free half, for the common case of no corpus token: what the cron's own logs say about
    # this entry. A dict with "error" is no answer at all, not a "no".
    db = corpus.get("in_corpus") if isinstance(corpus, dict) and "in_corpus" in corpus else None
    runs, _info = load_ingest_logs()
    log = entry_log_state(runs, entry, os.path.realpath(out)) if r["staged"] else None
    r["ingest_log"] = log
    log_state = (log or {}).get("state")
    r["ingested"] = bool(db or r.get("cron_marker") or (db is None and log_state == "ingested"))
    if db is False and log_state == "ingested":
        r["note"] = ("the cron's log records this entry as ingested, but the corpus has no row for "
                     "it -- the database answer is the one reported")
    r["state"] = ("ingested" if r["ingested"] else
                  "qc_excluded" if r.get("qc_excluded") else
                  "ingest_failed" if r["staged"] and log_state == "failed" else
                  "staged_pending_cron" if r["staged"] else
                  "not_staged")
    cron = None
    if r["state"] == "staged_pending_cron":
        cron = progress_health(runs, submit=submit_log_info(), info=_info)
        r["cron"] = {"verdict": cron["verdict"], "detail": cron.get("detail")}
    via = ("" if db is not None else
           " (No corpus token here, so the database was not asked; the cron's own logs were read "
           "instead.)")
    r["detail"] = {
        "qc_excluded": "A QC run: its entry is marked qc: true, so FRAN's ingester skips it. "
                       "Expected, not a failure.",
        "ingested": "FRAN has this search"
                    + (f" (the cron's log, {log['when']}: {log['outcome']})"
                       if db is None and log_state == "ingested" else "") + ".",
        "ingest_failed": (f"FRAN's cron tried this entry and it failed ({(log or {}).get('when')}: "
                          f"{(log or {}).get('detail')}). That is on FRAN's side -- the search is "
                          f"handed over and does not need re-staging; report the reason."),
        "staged_pending_cron": "Handed over; FRAN's ingest cron takes it on its next scan. This is "
                               "the expected state immediately after a run — not a failure."
                               + via,
        "not_staged": f"No entry at {entry}. Run `stage` first, or check why `check` refused.",
    }[r["state"]]
    if cron and cron["verdict"] in ("stuck", "not_running"):
        r["detail"] += (f" BUT FRAN's cron is {cron['verdict'].replace('_', ' ')}: "
                        f"{cron.get('detail')}. The search is safely staged; it waits on FRAN.")
    if r["staged"] and r.get("broken_links"):
        r["detail"] += (f" WARNING: {len(r['broken_links'])} link(s) point at files that no longer "
                        f"exist — the cron will skip this entry. Re-run `stage --force`.")

    # verify only RECORDS what it found. It never invents a "staged" receipt: written for a search
    # that was never staged, that receipt made stage answer already_staged and backfill skip the
    # search for ever. No receipt, nothing staged, nothing ingested -> nothing is written.
    prior = read_receipt(out)
    if r["ingested"]:
        status = "ingested"
    elif prior and prior.get("status"):
        status = prior["status"]
    elif r["staged"]:
        status = "staged"                  # the entry exists; only its receipt was missing
    else:
        r["receipt"] = None
        jout(r)
    receipt = prior or {"search_dir": out}
    receipt["status"] = status
    receipt["verified"] = {k: r[k] for k in ("state", "entry", "ingested") if k in r}
    if log and log.get("outcome"):
        receipt["verified"]["log"] = {k: log.get(k) for k in ("outcome", "when", "log")}
    if corpus and corpus.get("search"):
        receipt["search_id"] = corpus["search"]["search_id"]
    r["receipt"] = write_receipt(out, receipt)
    jout(r)


# ================================================================================ health ==
# "Is it staged?" and "will it be ingested?" are different questions, and for a week the second
# answer was no while the first looked fine: after 2026-09-17 13:55 the cron ingested nothing in 41
# consecutive runs, the recent ones all "0 ingested, 3 duplicate-skipped, 2 failed, ~181 still
# queued", because the same five FRAN_reports exports sort first, fail or duplicate, are never
# marked, and are picked again next run. None of the skill's drop entries in incoming/ appeared in
# a single log. Nothing on the skill's side could see that. This section reads what the cron itself
# wrote and says so.
#
# All of it is READ-ONLY: it never writes FRAN's code, its logs, or its database, and it needs no
# credential. It does not walk anything: one listdir of the log dir, the logs themselves (~14 KB
# each), one listdir of the drop dir and a stat per link in each entry.

INGEST_LOG_DIR = "/quobyte/proteomics-grp/de-limp/fran_refresh/logs"
SUBMIT_LOG = "auto_ingest_submit.log"
# The cron's own cadence and limits (brettsp crontab `23 */4 * * *`; fran_auto_ingest.sbatch
# `--time=08:00:00`). A run is "overdue" after three missed 4-hourly ticks.
NOT_RUNNING_AFTER_H = 12
JOB_TIME_LIMIT_H = 8
# Consecutive runs that ingested nothing while work was still queued before it counts as stuck.
# One or two such runs happen (a batch of duplicates); three in a row at 4 h each is half a day.
STUCK_AFTER_RUNS = 3
# A staged entry the cron has not reached in two days (~12 runs) is starved, even when the cron
# is ingesting OTHER searches -- which is the head-of-line case, where only the front moves.
STARVED_AFTER_H = 48
MAX_LOGS = 400
UNHEALTHY = ("stuck", "not_running", "stale_code")

_TS = r"(\d{4}-\d\d-\d\d \d\d:\d\d:\d\d)"
_RUN_START = re.compile(r"^===== (?:fran auto-ingest|auto_ingest) " + _TS + r" on (\S+)")
_RUN_DONE = re.compile(r"^===== done: (\d+) ingested, (\d+) duplicate-skipped, (\d+) failed, "
                       r"(\d+) still queued\W+" + _TS)
_ITEM = re.compile(r"^\[(\d+)/(\d+)\] (\S+) (.*)$")
_SKIP = re.compile(r"^  SKIP (.+?)  \((.*)\)\s*$")
_QCLAIM = re.compile(r"^  Q(\d+)\s+(\S+)\s+(/.*?)\s*$")
_CAND = re.compile(r"^  ([a-z]+)\s+(/.*?)\s*$")
_SLURM_STATE = re.compile(r"^State\s+:\s+(\S+)")
_LOG_NAME = re.compile(r"^auto_ingest_(\d+)\.out$")
_ENTRY_NAME = re.compile(r"__[0-9a-f]{8}$")          # entry_name()'s suffix
# Lines in a failed ingest's tail that are noise, not the reason (PG Farm prints the collation
# warning on every connection).
_NOISE = ("collation", "HINT:", "DETAIL:  The database", "--- stderr ---", "--- last output ---")


def _log_dir():
    return os.environ.get("FRAN_INGEST_LOG_DIR", INGEST_LOG_DIR)


def _epoch(ts):
    try:
        return time.mktime(time.strptime(ts, "%Y-%m-%d %H:%M:%S"))
    except (ValueError, OverflowError):
        return None


def _when(epoch):
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(epoch)) if epoch else None


def parse_ingest_log(text):
    """One auto_ingest_<jobid>.out -> what that run did.

    Parses the lines FRAN's auto_ingest.py prints (verified against its source on GitHub main and
    against real logs, 2026-09-24):
      ===== auto_ingest <ts> on <host> =====                     run start
        <engine>  <dir>                                          a scan candidate (first 40 only,
                                                                 and the scan output is cut to its
                                                                 last 4000 chars -- so absence
                                                                 from this list proves nothing)
        SKIP <search>  (<why>)                                   selected out, never attempted
        Q<id> <engine> <dir>                                     a queue row claimed this run
      [i/N] <engine> <search>                                    an attempt, followed by
            <dir>  /  -> <identity>  /  report: ...              ...where it came from
            OK in Ns | SKIPPED-DUPLICATE ... | FAILED ... | TIMEOUT ...
      ===== done: A ingested, B duplicate-skipped, C failed, D still queued — <ts> =====
    A log with no `done` line is a run still going, or one that died (PG Farm unreachable, the
    scan failed, SLURM killed it)."""
    run = {"started": None, "host": None, "finished": None, "complete": False,
           "ingested": None, "duplicate": None, "failed": None, "queued": None,
           "aborted": None, "slurm_state": None,
           "items": [], "skips": [], "candidates": [], "queue_claims": []}
    item = None
    for line in text.splitlines():
        m = _RUN_START.match(line)
        if m:
            if run["started"] is None:
                run["started"], run["host"] = _epoch(m.group(1)), m.group(2)
            item = None
            continue
        m = _RUN_DONE.match(line)
        if m:
            run.update(complete=True, ingested=int(m.group(1)), duplicate=int(m.group(2)),
                       failed=int(m.group(3)), queued=int(m.group(4)),
                       finished=_epoch(m.group(5)))
            item = None
            continue
        m = _ITEM.match(line)
        if m:
            item = {"engine": m.group(3), "search": m.group(4).split("  (")[0].strip(),
                    "dir": None, "identity": None, "outcome": None, "detail": None}
            run["items"].append(item)
            continue
        if item is not None and line.startswith("      "):
            s = line.strip()
            if item["dir"] is None and line[6:].startswith("/"):
                item["dir"] = line[6:].rstrip().rstrip("/")
            elif s.startswith("-> "):
                item["identity"] = s[3:].rstrip("/")
            elif s.startswith("OK in"):
                item["outcome"], item["detail"] = "ok", s
            elif s.startswith("SKIPPED-DUPLICATE"):
                item["outcome"], item["detail"] = "duplicate", s
            elif s.startswith(("FAILED", "TIMEOUT")):
                item["outcome"], item["detail"] = "failed", s
            elif s.startswith("DRY RUN"):
                item["outcome"] = "dry_run"
            elif (s.startswith("| ") and item["outcome"] == "failed"
                  and not any(n in s for n in _NOISE) and s[2:].strip()):
                item["error"] = s[2:].strip()[:300]      # the LAST real line wins: the reason
            continue
        if not line.startswith("      "):
            item = None
        m = _SKIP.match(line)
        if m:
            run["skips"].append((m.group(1).strip(), m.group(2)))
            continue
        m = _QCLAIM.match(line)
        if m:
            run["queue_claims"].append(m.group(3).rstrip("/"))
            continue
        m = _CAND.match(line)
        if m:
            run["candidates"].append(m.group(2).rstrip("/"))
            continue
        if line.startswith(("ABORT", "SCAN FAILED")) or "ABORT:" in line:
            run["aborted"] = line.strip()[:200]
            continue
        m = _SLURM_STATE.match(line)
        if m:
            run["slurm_state"] = m.group(1)
    return run


def load_ingest_logs(log_dir=None, max_logs=MAX_LOGS, max_bytes=4 << 20):
    """(runs newest first, info) -- or (None, info) when the log dir cannot be read.

    Newest by SLURM job id from the filename, not by mtime: ids only increase on one cluster, and
    sorting by name needs no stat of every file on a network filesystem."""
    d = log_dir or _log_dir()
    info = {"log_dir": d}
    try:
        names = os.listdir(d)
    except OSError as e:
        info["error"] = f"cannot read {d}: {e.strerror or e}"
        return None, info
    logs = sorted(((int(m.group(1)), n) for n in names for m in [_LOG_NAME.match(n)] if m),
                  reverse=True)
    runs = []
    for jobid, n in logs[:max_logs]:
        p = os.path.join(d, n)
        try:
            st = os.stat(p)
            with open(p, errors="replace") as fh:
                text = fh.read(max_bytes)
        except OSError:
            continue
        r = parse_ingest_log(text)
        r.update(jobid=jobid, log=p, mtime=st.st_mtime)
        runs.append(r)
    info.update(n_logs=len(logs), n_read=len(runs))
    if runs:
        info["newest"] = _when(runs[0]["started"] or runs[0]["mtime"])
        info["oldest"] = _when(runs[-1]["started"] or runs[-1]["mtime"])
    return runs, info


def _tail_line(path, nbytes=4096):
    try:
        with open(path, "rb") as fh:
            fh.seek(0, os.SEEK_END)
            fh.seek(max(0, fh.tell() - nbytes))
            lines = [x for x in fh.read().decode(errors="replace").splitlines() if x.strip()]
        return lines[-1] if lines else None
    except OSError:
        return None


def submit_log_info(log_dir=None, now=None):
    """Is the cron still SUBMITTING? cron_auto_ingest.sh appends one line per tick -- `Submitted
    batch job N`, or `skip: 1 ... already pending/running` when the last one has not started."""
    now = now or time.time()
    p = os.path.join(log_dir or _log_dir(), SUBMIT_LOG)
    try:
        st = os.stat(p)
    except OSError:
        return {"path": p, "exists": False}
    return {"path": p, "exists": True, "age_h": round((now - st.st_mtime) / 3600, 1),
            "last_line": _tail_line(p)}


def progress_health(runs, now=None, submit=None, info=None):
    """Is the cron making progress? healthy | stuck | not_running | unknown, with the evidence."""
    now = now or time.time()
    res = {"verdict": "unknown", "logs": info or {}, "cron": submit or {}}
    if runs is None:
        res["detail"] = (info or {}).get("error") or "the cron's logs could not be read"
        return res
    if not runs:
        res.update(verdict="not_running", detail=f"no auto_ingest logs in {(info or {}).get('log_dir')}")
        return res

    def t(r):
        return r["started"] or r["mtime"]

    newest = runs[0]
    running = (not newest["complete"] and not newest["aborted"] and not newest["slurm_state"]
               and now - t(newest) < JOB_TIME_LIMIT_H * 3600)
    finished = runs[1:] if running else runs
    streak = aborted = 0
    for r in finished:
        if r["complete"] and r["ingested"]:
            break
        streak += 1
        aborted += not r["complete"]
    last_ok = next((r for r in finished if r["complete"] and r["ingested"]), None)
    last_done = next((r for r in finished if r["complete"]), None)
    queued = last_done["queued"] if last_done else None
    age_h = (now - t(newest)) / 3600

    def brief(r):
        return None if r is None else {
            "jobid": r["jobid"], "started": _when(t(r)), "finished": _when(r["finished"]),
            "state": ("running" if r is newest and running else
                      "complete" if r["complete"] else "aborted"),
            "ingested": r["ingested"], "duplicate": r["duplicate"], "failed": r["failed"],
            "queued": r["queued"], "aborted": r["aborted"] or (
                None if r["complete"] or (r is newest and running) else
                f"no summary line (SLURM state {r['slurm_state'] or 'unknown'})"),
            "log": r["log"]}

    res.update(last_run=brief(newest), last_run_age_h=round(age_h, 1),
               last_ingest=brief(last_ok), queued=queued,
               consecutive_runs_without_ingest=streak, of_which_aborted=aborted)
    if last_ok:
        res["last_ingest_age_days"] = round((now - (last_ok["finished"] or t(last_ok))) / 86400, 1)
    sub_age = (submit or {}).get("age_h")
    res["cron"]["submitting"] = bool(sub_age is not None and sub_age <= NOT_RUNNING_AFTER_H)

    if age_h > NOT_RUNNING_AFTER_H or (sub_age is not None and sub_age > NOT_RUNNING_AFTER_H):
        last_word = (submit or {}).get("last_line")
        res.update(verdict="not_running",
                   detail=f"last cron run started {age_h:.0f} h ago"
                          + (f"; cron last wrote its submit log {sub_age:.0f} h ago"
                             if sub_age is not None else "")
                          # `skip: 1 ... already pending/running` = a job stuck in the queue
                          + (f" ({last_word[20:] if last_word[:4].isdigit() else last_word})"
                             if last_word else ""))
    elif streak >= STUCK_AFTER_RUNS and (queued is None or queued > 0):
        since = (f"last ingest {_when(last_ok['finished'] or t(last_ok))}" if last_ok
                 else f"no ingest in the {len(finished)} logs kept")
        res.update(verdict="stuck",
                   detail=f"{streak} consecutive cron runs ingested nothing"
                          + (f" ({aborted} of them died before finishing)" if aborted else "")
                          + (f"; {queued} searches still queued" if queued is not None else "")
                          + f"; {since}")
    else:
        res.update(verdict="healthy",
                   detail=("a run is in progress; " if running else "")
                          + (f"last ingest {_when(last_ok['finished'] or t(last_ok))}" if last_ok
                             else "nothing ingested in the logs kept, and nothing waiting"))
    return res


def entry_log_state(runs, entry, output_dir=None):
    """What the cron's logs say about ONE staged entry: ingested | failed | never_reached.

    Matched on the entry's path (what find_uningested reports for a drop entry -- it is a real
    directory, so realpath does not change it), its name (what auto_ingest prints for SKIP), and
    the search's real directory (the identity a queue row or a future manifest-aware cron prints).
    A duplicate counts as ingested: the guard refused it because the corpus already holds it."""
    entry = entry.rstrip("/")
    # The NAME is only a key when it is an entry name (`<dir>__<8 hex>`, unique by construction).
    # A bare search folder name is not: every parallel-chain search is called `search_out`.
    name = os.path.basename(entry) if _ENTRY_NAME.search(os.path.basename(entry)) else None
    keys = {entry}
    if output_dir:
        keys.add(output_dir.rstrip("/"))
    res = {"state": "never_reached", "outcome": None, "when": None, "log": None, "detail": None,
           "attempts": 0, "listed_as_candidate": 0, "logs_searched": len(runs or [])}
    if not runs:
        res["state"] = "unknown" if runs is None else "never_reached"
        return res
    hits = []                                     # (run, outcome, detail), newest first
    for r in runs:
        when = r["finished"] or r["started"] or r["mtime"]
        for it in r["items"]:
            if (it["dir"] in keys or it["identity"] in keys
                    or (name and it["dir"] and os.path.basename(it["dir"]) == name)):
                hits.append((r, when, it["outcome"] or "attempted",
                             it.get("error") or it["detail"]))
        for sname, why in r["skips"]:
            if name and sname == name:
                hits.append((r, when, "skipped", why))
        res["listed_as_candidate"] += sum(1 for c in r["candidates"] if c in keys)
        res["listed_as_candidate"] += sum(1 for c in r["queue_claims"] if c in keys)
    res["attempts"] = sum(1 for h in hits if h[2] not in ("skipped", "dry_run"))
    done = next((h for h in hits if h[2] in ("ok", "duplicate")), None)
    last = next((h for h in hits if h[2] != "dry_run"), None)
    pick = done or last
    if pick:
        r, when, outcome, detail = pick
        res.update(state="ingested" if done else "failed", outcome=outcome, when=_when(when),
                   log=r["log"], detail=detail)
    return res


def incoming_health(runs, now=None, drop=None):
    """Every entry in the drop dir: age, who staged it, broken links, and what the logs say."""
    now = now or time.time()
    drop = drop or os.environ.get("FRAN_DROP_DIR", DROP_DIR)
    res = {"drop_dir": drop, "verdict": "unknown", "entries": []}
    try:
        ents = sorted((e for e in os.scandir(drop) if e.is_dir(follow_symlinks=False)),
                      key=lambda e: e.name)
    except OSError as e:
        res["detail"] = f"cannot read {drop}: {e.strerror or e}"
        return res
    for e in ents:
        try:
            age_d = (now - e.stat(follow_symlinks=False).st_mtime) / 86400
        except OSError:
            age_d = None
        man = {}
        try:
            with open(os.path.join(e.path, MANIFEST)) as fh:
                man = json.load(fh)
        except (OSError, ValueError):
            pass
        try:
            broken = [f for f in os.listdir(e.path)
                      if os.path.islink(os.path.join(e.path, f))
                      and not os.path.exists(os.path.join(e.path, f))]
        except OSError:
            broken = []
        st = entry_log_state(runs, e.path, man.get("output_dir"))
        # Entries staged before the QC rule existed carry no `qc` key: apply the rule now, so a QC
        # run already sitting in the drop box is visible (stage --qc marks it for the ingester).
        if isinstance(man.get("qc"), bool):
            qc, qrule = man["qc"], man.get("qc_rule")
        else:
            qc, qrule = is_qc_run(man.get("output_dir") or e.path, names=[man.get("search_name")])
        res["entries"].append({
            "entry": e.name, "age_days": None if age_d is None else round(age_d, 1),
            "staged_by": man.get("staged_by"), "staged_at": man.get("staged_at"),
            "engine": man.get("engine"), "output_dir": man.get("output_dir"),
            "search_name": man.get("search_name"), "state": st["state"],
            "outcome": st["outcome"], "when": st["when"], "detail": st["detail"],
            "attempts": st["attempts"], "broken_links": broken,
            "qc": qc, "qc_rule": qrule, "qc_marked": man.get("qc") is True})
    counts = collections.Counter(x["state"] for x in res["entries"])
    res["n_entries"] = len(res["entries"])
    res["by_state"] = dict(counts)
    waiting = [x for x in res["entries"] if x["state"] == "never_reached"]
    oldest = max((x["age_days"] or 0 for x in waiting), default=0)
    res["oldest_never_reached_days"] = round(oldest, 1) if waiting else None
    if runs is None:
        res["detail"] = "the cron's logs could not be read, so no entry's state is known"
    elif waiting and oldest * 24 > STARVED_AFTER_H:
        res.update(verdict="starved",
                   detail=f"{len(waiting)} of {len(res['entries'])} staged entries never reached "
                          f"by the cron (oldest {oldest:.0f} d)")
    else:
        res.update(verdict="ok", detail=f"{len(res['entries'])} staged entries; "
                   + ", ".join(f"{n} {k}" for k, n in sorted(counts.items())))
    if counts.get("failed"):
        res["failed"] = [{"entry": x["entry"], "detail": x["detail"]}
                         for x in res["entries"] if x["state"] == "failed"]
    unmarked = [x for x in res["entries"] if x["qc"] and not x["qc_marked"]]
    if unmarked:
        res["qc_unmarked"] = [{"entry": x["entry"], "output_dir": x["output_dir"],
                               "qc_rule": x["qc_rule"]} for x in unmarked]
    return res


# ------------------------------------------------------------------------- ingest code --
# FRAN's ingest changes often, and the HIVE copy is a folder of files, not a checkout, so nothing
# tells it when it falls behind. FRAN's own guard (publish_manifest.py) keeps content md5s in PG
# Farm -- reading those needs a database credential, which this file never touches. GitHub `main` is
# public, so compare against that: raw.githubusercontent.com for content (no API quota), and the
# commits API (60 requests/h unauthenticated, shared by everyone behind a login node's address)
# only for a file that DIFFERS, to say whether HIVE runs an older commit or an edit nobody pushed.
FRAN_REPO = "bsphinney/FRAN"
RAW_URL = "https://raw.githubusercontent.com/{repo}/{ref}/ingest/{name}"
COMMITS_URL = "https://api.github.com/repos/{repo}/commits?path=ingest/{name}&sha=main&per_page={n}"
# publish_manifest.REFUSE_FILES: "scripts whose staleness silently corrupts the corpus".
REFUSE_FILES = ("corpus_ingest.py", "spectronaut_to_corpus.py", "diann_to_corpus.py", "versions.py")
INGEST_CODE_FILES = REFUSE_FILES + (
    # what runs the cron and decides which searches it reaches
    "auto_ingest.py", "find_uningested.py", "fran_queue.py",
    "cron_auto_ingest.sh", "fran_auto_ingest.sbatch",
    # what corpus_ingest.py imports to read the engines this skill stages, and the XIC lane
    "radiant_to_corpus.py", "raw_metadata.py", "engine_fasta.py", "engine_version.py",
    "diann_xic_to_lance.py",
)
HISTORY_CAP = 15           # commits of one file searched for the HIVE copy
MAIN_TTL_S = 600           # main moves; a cached answer older than this is re-fetched
COMMITS_TTL_S = 3600
USER_AGENT = "ucdavis-proteomics-core-pipeline/fran_deposit.py"


def _default_fetch(url, timeout):
    """(status | None, body, headers). Never raises: no network is an answer, not a crash."""
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:          # noqa: S310
            return resp.status, resp.read(), {k.lower(): v for k, v in resp.headers.items()}
    except urllib.error.HTTPError as e:
        try:
            body = e.read()
        except Exception:                                                   # noqa: BLE001
            body = b""
        return e.code, body, {k.lower(): v for k, v in (e.headers or {}).items()}
    except Exception as e:                                                  # noqa: BLE001
        return None, f"{type(e).__name__}: {e}"[:200], {}


def _cache_path():
    return os.environ.get("FRAN_HEALTH_CACHE") or os.path.expanduser(
        "~/.cache/ucdavis-proteomics-core-pipeline/fran_health.json")


def _load_cache(p):
    try:
        with open(p) as fh:
            c = json.load(fh)
        return c if isinstance(c, dict) else {}
    except (OSError, ValueError):
        return {}


def _save_cache(p, c):
    try:
        text = json.dumps(c)        # RuntimeError if a late fetch thread is still adding to it
        os.makedirs(os.path.dirname(p), exist_ok=True)
        tmp = f"{p}.{os.getpid()}.tmp"
        with open(tmp, "w") as fh:
            fh.write(text)
        os.replace(tmp, p)
    except (OSError, RuntimeError, TypeError, ValueError):
        pass                        # a cache we cannot keep only costs a re-fetch


def _md5_file(p):
    try:
        h = hashlib.md5()                                                   # noqa: S324
        with open(p, "rb") as fh:
            for b in iter(lambda: fh.read(1 << 20), b""):
                h.update(b)
        return h.hexdigest()
    except OSError:
        return None


def ingest_code_health(ingest_dir=None, fetch=None, cache_path=None, now=None, budget_s=10,
                       files=INGEST_CODE_FILES, history_cap=HISTORY_CAP):
    """Each ingest file on HIVE vs FRAN GitHub main:
         current | stale (matches an older commit) | local_modification (matches none of the last
         `history_cap` commits) | differs (history not checked) | missing (on main, not on HIVE) |
         not_on_main | unknown (no network)
    Verdict: stale if any file is stale or missing, or a REFUSE_FILES file differs at all;
    modified if only other files carry local edits; current; or unknown."""
    now = now or time.time()
    ing = ingest_dir or first_readable(INGEST_DIRS, isdir=True)
    res = {"verdict": "unknown", "ingest_dir": ing, "compared_with": f"github.com/{FRAN_REPO} main",
           "files": []}
    if not ing:
        res["detail"] = ("FRAN's ingest directory is not readable here (not on HIVE, or not a "
                         "Core account), so there is nothing to compare")
        return res
    fetch = fetch or _default_fetch
    cp = cache_path or _cache_path()
    cache = _load_cache(cp)
    for k in ("main", "commits", "blob"):
        cache.setdefault(k, {})
    deadline = time.monotonic() + budget_s

    def left():
        return deadline - time.monotonic()

    def md5_at(ref, name):
        """(status, md5) of ingest/<name> at a ref; status None = could not ask."""
        if ref != "main" and f"{ref}:{name}" in cache["blob"]:
            v = cache["blob"][f"{ref}:{name}"]
            return (200, v) if v else (404, None)
        if ref == "main":
            c = cache["main"].get(name)
            if c and now - c.get("at", 0) < MAIN_TTL_S:
                return c["status"], c.get("md5")
        if left() < 0.5:
            return None, None
        st, body, _ = fetch(RAW_URL.format(repo=FRAN_REPO, ref=ref, name=name), min(8.0, left()))
        if st == 200 and isinstance(body, (bytes, bytearray)):
            md5 = hashlib.md5(body).hexdigest()                            # noqa: S324
        elif st == 404:
            md5 = None
        else:
            return None, None
        if ref == "main":
            cache["main"][name] = {"status": st, "md5": md5, "at": now}
        else:
            cache["blob"][f"{ref}:{name}"] = md5          # a commit's content never changes
        return st, md5

    def history(name):
        """([{sha, date}] newest first, None) or (None, why-not)."""
        c = cache["commits"].get(name)
        if c and now - c.get("at", 0) < COMMITS_TTL_S:
            return c["list"], None
        rl = cache.get("rate_limited_until") or 0
        if rl > now:
            return None, f"GitHub API rate limit until {_when(rl)}"
        if left() < 0.5:
            return None, "time budget spent"
        st, body, hdr = fetch(COMMITS_URL.format(repo=FRAN_REPO, name=name, n=history_cap),
                              min(8.0, left()))
        if st == 200:
            try:
                lst = [{"sha": x["sha"], "date": (x.get("commit") or {}).get("committer", {})
                        .get("date")} for x in json.loads(body)][:history_cap]
            except (ValueError, TypeError, KeyError):
                return None, "unreadable commit list"
            cache["commits"][name] = {"list": lst, "at": now}
            return lst, None
        if st in (403, 429) and (hdr.get("x-ratelimit-remaining") == "0" or st == 429):
            try:
                cache["rate_limited_until"] = float(hdr.get("x-ratelimit-reset") or now + 3600)
            except ValueError:
                cache["rate_limited_until"] = now + 3600
            return None, f"GitHub API rate limit until {_when(cache['rate_limited_until'])}"
        return None, (f"GitHub API answered {st}" if st else "no network")

    # Everything that touches the network runs in DAEMON threads under one deadline. urlopen's
    # timeout covers each socket read but not a hung DNS lookup, and a ThreadPoolExecutor's threads
    # are joined at interpreter exit -- either could hold `health` well past its budget. This
    # returns at the deadline (<= budget_s + 0.5 s) whatever GitHub is doing, and a file with no
    # answer by then is `unknown`, never an error.
    mains, done = {}, {}

    def fetch_main(n):
        mains[n] = md5_at("main", n)

    def compare():
        workers = [threading.Thread(target=fetch_main, args=(n,), daemon=True) for n in files]
        for w in workers:
            w.start()
        for w in workers:
            w.join(max(0.0, left()))
        for name in files:
            if name not in mains:
                continue
            hive = _md5_file(os.path.join(ing, name))
            st, main_md5 = mains[name]
            f = {"file": name, "hive_md5": hive, "main_md5": main_md5, "refuse_gate": name in REFUSE_FILES}
            if st is None:
                f.update(status="unknown", detail="could not reach GitHub")
            elif st == 404:
                f.update(status="not_on_main",
                         detail=("on HIVE but not on main (deleted or renamed upstream?)" if hive
                                 else "on neither -- nothing to compare"))
            elif hive is None:
                f.update(status="missing", detail="on main but not in the HIVE ingest folder")
            elif hive == main_md5:
                f.update(status="current")
            else:
                commits, why = history(name)
                if commits is None:
                    f.update(status="differs", detail=f"differs from main; history not checked ({why})")
                else:
                    match = None
                    for i, c in enumerate(commits):
                        cst, cmd5 = md5_at(c["sha"], name)
                        if cst is None:
                            why = "could not fetch an older commit"
                            break
                        if cmd5 == hive:
                            match = (i, c)
                            break
                    if match:
                        i, c = match
                        f.update(status="stale", matches={"sha": c["sha"][:10], "date": c["date"]},
                                 detail=f"stale: matches {c['sha'][:7]} from {(c['date'] or '')[:10]}"
                                        f" ({i} newer commit(s) of this file on main)")
                    elif why:
                        f.update(status="differs", detail=f"differs from main; {why}")
                    else:
                        f.update(status="local_modification",
                                 detail=f"local modification: matches none of the last "
                                        f"{len(commits)} commit(s) of this file on main")
            done[name] = f

    t = threading.Thread(target=compare, daemon=True)
    t.start()
    t.join(budget_s + 0.5)
    finished = not t.is_alive()
    for name in files:
        f = done.get(name)
        if f is None:
            f = {"file": name, "hive_md5": None, "main_md5": None, "refuse_gate": name in REFUSE_FILES,
                 "status": "unknown", "detail": f"no answer from GitHub within {budget_s:g} s"}
        res["files"].append(f)
    if finished:
        _save_cache(cp, cache)      # a thread still running may still be writing to it

    by = collections.defaultdict(list)
    for f in res["files"]:
        by[f["status"]].append(f["file"])
    res["by_status"] = {k: sorted(v) for k, v in by.items()}
    risky = [f["file"] for f in res["files"] if f["refuse_gate"]
             and f["status"] in ("local_modification", "differs")]
    if by.get("stale") or by.get("missing") or risky:
        res["verdict"] = "stale"
    elif by.get("local_modification") or by.get("differs"):
        res["verdict"] = "modified"
    elif by.get("unknown"):
        res["verdict"] = "unknown"
    else:
        res["verdict"] = "current"
    bits = []
    for f in res["files"]:
        if f["status"] in ("stale", "local_modification", "differs", "missing"):
            bits.append(f"{f['file']} {f['detail']}")
    res["detail"] = ("; ".join(bits) if bits else
                     "every ingest file matches main" if res["verdict"] == "current" else
                     "could not reach GitHub")
    return res


def health_report(now=None, fetch=None, check_code=True, code_budget_s=10, cache_path=None):
    """The one answer to "can a staged search reach FRAN right now?"."""
    now = now or time.time()
    runs, info = load_ingest_logs()
    prog = progress_health(runs, now, submit_log_info(now=now), info)
    inc = incoming_health(runs, now)
    code = (ingest_code_health(fetch=fetch, now=now, budget_s=code_budget_s, cache_path=cache_path)
            if check_code else {"verdict": "skipped", "detail": "not checked"})
    if prog["verdict"] == "not_running":
        verdict = "not_running"
    elif prog["verdict"] == "stuck" or inc["verdict"] == "starved":
        verdict = "stuck"
    elif code["verdict"] == "stale":
        verdict = "stale_code"
    elif prog["verdict"] == "unknown":
        verdict = "unknown"
    else:
        verdict = "healthy"

    parts = []
    if prog["verdict"] == "not_running":
        parts.append(f"FRAN ingest NOT RUNNING: {prog['detail']}")
    elif prog["verdict"] == "stuck":
        parts.append(f"FRAN ingest STUCK: {prog['detail']}")
    elif prog["verdict"] == "healthy":
        parts.append(f"FRAN ingest cron running ({prog['detail']})")
    else:
        parts.append(f"FRAN ingest progress unknown ({prog.get('detail')})")
    if inc["verdict"] == "starved" or (inc.get("by_state") or {}).get("never_reached"):
        parts.append(inc["detail"])
    if inc.get("failed"):
        parts.append(f"{len(inc['failed'])} staged entr{'y' if len(inc['failed']) == 1 else 'ies'}"
                     f" FAILED ingest: " + "; ".join(f"{x['entry']}: {x['detail']}"
                                                     for x in inc["failed"][:2]))
    if inc.get("qc_unmarked"):
        n = len(inc["qc_unmarked"])
        parts.append(f"{n} staged entr{'y is a QC run' if n == 1 else 'ies are QC runs'} not yet "
                     f"marked ({', '.join(x['entry'] for x in inc['qc_unmarked'][:3])}; "
                     f"stage --out <its out dir> --qc marks it)")
    parts.append({"current": "ingest code current with FRAN main",
                  "stale": f"ingest code STALE: {code.get('detail')}",
                  "modified": f"ingest code modified on HIVE: {code.get('detail')}",
                  "skipped": "ingest code not checked"}.get(
                      code["verdict"], f"ingest code not checked ({code.get('detail')})"))
    return {"verdict": verdict, "healthy": verdict == "healthy", "summary": "; ".join(parts),
            "checked_at": _when(now), "progress": prog, "incoming": inc, "ingest_code": code}


# ------------------------------------------------------------------- health status file --
# stage() runs inside EVERY search job, on a compute node, inside the job's time limit. So it must
# not reach GitHub, the database, or scan 170 logs on a network mount. The full check (`health`)
# writes its verdict to one small file and stage only reads that. The file lives in the drop dir's
# PARENT (/quobyte/proteomics-grp/fran/, group-writable), never inside incoming/: nothing in the drop
# dir may be anything but a drop entry.
HEALTH_FILE = "ingest_health.json"
HEALTH_STALE_H = 12          # three missed cron ticks: older than this is "unknown", said so
HEALTH_READ_TIMEOUT_S = 2.0  # a hung mount must not hold the search job


def health_file(drop=None):
    """/quobyte/proteomics-grp/fran/ingest_health.json: beside the drop dir, not in it.
    FRAN_HEALTH_FILE moves it."""
    if not drop and os.environ.get("FRAN_HEALTH_FILE"):
        return os.environ["FRAN_HEALTH_FILE"]
    drop = drop or os.environ.get("FRAN_DROP_DIR", DROP_DIR)
    return os.path.join(os.path.dirname(drop.rstrip("/")) or "/", HEALTH_FILE)


def write_health_status(h, drop=None):
    """Record a `health` verdict for stage() to read. Group-writable, replaced atomically, so any
    Core member's next `health` run can refresh it. Never raises."""
    p = health_file(drop)
    rec = {"checked_at": _utc_now(), "verdict": h["verdict"], "healthy": h["healthy"],
           "summary": h["summary"], "ingest_code": (h.get("ingest_code") or {}).get("verdict"),
           "checked_by": getpass.getuser(), "host": os.uname().nodename,
           "written_by": "fran_deposit.py health"}
    try:
        tmp = f"{p}.{os.getpid()}.tmp"
        with open(tmp, "w") as fh:
            json.dump(rec, fh, indent=2)
        os.chmod(tmp, 0o664)
        os.replace(tmp, p)
        return {"written": p}
    except OSError as e:
        return {"written": None, "why": f"{type(e).__name__}: {e.strerror or e}"}


def read_health_status(now=None, drop=None, timeout_s=HEALTH_READ_TIMEOUT_S):
    """(fran_health, warning) from the status file, or (None, None) when it is missing,
    unreadable, or the read does not finish within `timeout_s` -- then stage says nothing.
    Older than HEALTH_STALE_H: verdict unknown and NO warning (the summary says when it was last
    checked). Only a fresh unhealthy verdict warns."""
    box = {}

    def work():
        try:
            with open(health_file(drop)) as fh:
                box["rec"] = json.load(fh)
        except (OSError, ValueError):
            pass

    t = threading.Thread(target=work, daemon=True)
    t.start()
    t.join(timeout_s)
    rec = box.get("rec")
    if not isinstance(rec, dict) or not rec.get("verdict"):
        return None, None
    now = now or time.time()
    try:
        checked = datetime.datetime.strptime(rec["checked_at"], "%Y-%m-%dT%H:%M:%SZ").replace(
            tzinfo=datetime.timezone.utc).timestamp()
    except (KeyError, TypeError, ValueError):
        checked = None
    when = _when(checked) if checked else "at an unknown time"
    if checked is None or now - checked > HEALTH_STALE_H * 3600:
        # Old news is not bad news: a stale verdict is UNKNOWN, with no warning -- exactly like a
        # missing file -- so the job's Slack post does not call a healthy FRAN unhealthy.
        return {"verdict": "unknown", "checked_at": rec.get("checked_at"),
                "summary": f"FRAN ingest health unknown (last checked {when})"}, None
    fh = {"verdict": rec["verdict"], "summary": rec.get("summary"), "checked_at": rec["checked_at"]}
    return fh, (rec.get("summary") if rec["verdict"] in UNHEALTHY else None)


def attach_health_status(res, now=None):
    """Put the last `health` verdict on a stage() result, and say it on stderr when it is bad or
    stale. A health problem is FRAN-side and never fails the stage -- the search IS handed over."""
    if os.environ.get("FRAN_HEALTH", "").lower() in ("off", "0", "no", "false"):
        return res
    try:
        fh, warn = read_health_status(now=now)
    except Exception:                                                       # noqa: BLE001
        return res
    if fh:
        res["fran_health"] = fh
    if warn:
        res["health_warning"] = warn
        sys.stderr.write(f"[fran_deposit] WARNING: {warn}. The search IS handed over; this is "
                         f"on FRAN's side. Detail: python3 {os.path.abspath(__file__)} health\n")
    return res


ALERT_REPEAT_S = 24 * 3600


def _notify(mod, text):
    """(sent, status) through notify_slack.send_alert -- its API for other skill scripts. It
    resolves the webhook itself, never raises, prints nothing, and its status line never carries
    the URL."""
    sent, status = mod.send_alert(text, title="FRAN ingest code", with_status=True)
    return bool(sent), status


def alert_reason(h):
    """What `health --alert` pages about, or None. ONLY what FRAN's own runner cannot see about
    itself: the ingest CODE on HIVE is stale or modified against GitHub main, or lacks a file the
    runner imports. A stuck cron, unreached entries and FASTA mismatches are FRAN's runner's to
    alert on (FRAN branch fix/auto-ingest-starvation, 78374dc + 2163070) -- paging on them here too
    would be a duplicate page. They are still in the verdict file and the printed report."""
    code = h.get("ingest_code") or {}
    if code.get("verdict") not in ("stale", "modified"):
        return None
    return f"FRAN's ingest code on HIVE is {code['verdict']} vs GitHub main: {code.get('detail')}"


def send_alert(h, cache_path=None, now=None):
    """Post a short Slack alert when alert_reason() finds something, through scripts/
    notify_slack.py when it exists. A no-op otherwise. The webhook is that module's business: this
    never reads or prints it. The same alert is not re-posted within 24 h, because `health` runs on
    a schedule and a stale file stays stale until someone syncs it."""
    now = now or time.time()
    reason = alert_reason(h)
    if not reason:
        return {"sent": False, "why": "nothing to page about: the ingest code matches GitHub main "
                                      "(cron progress is FRAN's runner's to alert on)"}
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import notify_slack                                                 # noqa: PLC0415
    except Exception:                                                       # noqa: BLE001
        return {"sent": False, "why": "notify_slack.py is not installed"}
    cp = cache_path or _cache_path()
    cache = _load_cache(cp)
    key = hashlib.sha1(reason.encode()).hexdigest()                          # noqa: S324
    last = cache.get("last_alert") or {}
    if last.get("key") == key and now - last.get("at", 0) < ALERT_REPEAT_S:
        return {"sent": False, "why": f"same alert already sent at {_when(last['at'])}"}
    try:
        ok, status = _notify(notify_slack, reason)
    except Exception as e:                                                  # noqa: BLE001
        return {"sent": False, "why": f"notify_slack failed ({type(e).__name__})"}
    if ok:
        cache["last_alert"] = {"key": key, "at": now}
        _save_cache(cp, cache)
    return {"sent": ok, "why": None if ok else status}


def health(a):
    """The full check. Also the ONLY writer of the status file stage() reads."""
    h = health_report(check_code=not a.no_code)
    h["status_file"] = write_health_status(h)
    h["alert_reason"] = alert_reason(h)
    if a.alert:
        h["alert"] = send_alert(h)
    jout(h)


# ============================================================================== backfill ==
# Staging became automatic partway through the skill's life, so searches it ran before that (or
# whose hand-over was skipped because a session ended early) sit on HIVE and FRAN never hears of
# them: its cron only scans incoming/ and FRAN_reports/, never the Core's service trees. `backfill`
# finds those searches and stages the eligible ones through the same stage() path, with the same
# reason codes. Dry run unless --apply.
#
# Finding them means walking NFS, and a walk of the Core service trees is exactly the kind of
# metadata storm that gets an account flagged on a login node. So a walk runs inside SLURM
# (--sbatch writes the job), depth-limited, time-boxed, never following a symlink.

# Where Core searches live. The Quobyte SERVICE tree and the Flinders service tree hold the facility's
# customer work; ~/proteomics-pipeline is the skill's install on HIVE (verified 2026-09-24 to hold
# only presets/, references/, scripts/ for brettsp -- kept so a member's sessions there are found).
BACKFILL_ROOTS = ["/quobyte/proteomics-grp/SERVICE",
                  "/nfs/lssc0/flinders/proteomics/Data/lab/service"]
# A search whose real path is inside one of these is a Core search whatever account ran it.
CORE_PREFIXES = (GROUP_ROOT + "/", "/nfs/lssc0/flinders/proteomics/")
GROUP_NAME = "proteomics-grp"
# The group also holds ~30 teaching accounts (proteomics-class-NN). Their work is coursework, not
# facility searches, and FRAN deliberately keeps teaching data out of the corpus.
TEACHING_ACCOUNT = re.compile(r"-class-\d+$")
# Trees that produce engine output which is NOT a corpus search. Mirrors FRAN's own
# find_uningested.DEFAULT_EXCLUDES (STAN QC, QC watchers, smoke tests, scratch), plus the drop dir.
BACKFILL_EXCLUDES = FRAN_DEFAULT_EXCLUDES + ("/fran/incoming/",)
# A directory holding any of these is a search output: the walk stops there (its children are its
# own outputs) and the directory becomes a candidate.
SEARCH_HINTS = {"search_provenance.json", RECEIPT, "report.parquet", "report.tsv",
                "step5_report.sbatch", "step3_fulcrum.sbatch", "dia-quant-output",
                "radiant_results", "fulcrum-results"}
# What proves THIS skill ran a search, from file names alone. search_provenance.json is not enough
# by itself: in hive_remote mode run_search.py writes it on the laptop, and the Dupanloup search on
# HIVE has none. The DE-LIMP app writes the same step*_*.sbatch names, but its jobs are called
# diann_<name>_<step>, while diann_parallel.py's are s1_libpred / s1b_window / s5_report -- so the
# SLURM logs tell the two apart without opening a file.
SKILL_FILES = ("search_provenance.json", RECEIPT, "step3_fulcrum.sbatch",
               "parallel_input_files.txt", "radiant_input_files.txt")
SKILL_LOG = re.compile(r"^s(?:1|1b|5)_(?:libpred|window|report)_\d+(?:_\d+)?\.log$")
PRUNE_SUFFIX = (".d", ".raw", ".wiff", ".wiff2", ".mzml", ".mzxml", ".lance", ".sne", ".zip")
PRUNE_NAMES = {".snapshot", ".git", "__pycache__", ".Trash", "lost+found", ".ipynb_checkpoints",
               "node_modules"}
REC_JSON = ".recovery.json"          # checkpoint.py's session record
LOGIN_LIST_MAX = 10                  # --list entries a login node may check without SLURM
WALK_BUDGET_S = 1200
WALK_MAX_DEPTH = 9                   # Flinders: service/on_campus/<lab>/<yr>/<proj>/<set>/<session>/output/search = 8


def skill_marker(names):
    """The file name that proves this skill ran the search in a directory, or None."""
    for f in SKILL_FILES:
        if f in names:
            return f
    return next((n for n in sorted(names) if SKILL_LOG.match(n)), None)


def discover(roots, max_depth=WALK_MAX_DEPTH, budget_s=WALK_BUDGET_S, excludes=BACKFILL_EXCLUDES,
             clock=time.monotonic):
    """Breadth-first walk of `roots` for search output directories. Never follows a symlink,
    prunes raw-data containers (.d/.raw/...), stops descending at a search directory, and stops
    altogether when the time budget is spent (`truncated`).

    The budget is shared in rounds: each unfinished root gets an equal part of what is left, and
    time a small root does not use goes back into the pot for the next round. A single equal
    split starved the one tree that mattered: on HIVE (2026-09-24) the Flinders service tree got
    1/7 of 18 min because six near-empty ~/proteomics-pipeline roots each held a share, and the
    walk stopped with 138,779 directories still queued."""
    t0 = clock()
    found, sessions = [], []
    stats = {"dirs_walked": 0, "unreadable": 0, "truncated": False, "max_depth": max_depth,
             "budget_s": budget_s, "per_root": {}}
    pending = []
    for root in roots:
        rstat = {"dirs_walked": 0, "found": 0, "truncated": False}
        stats["per_root"][root] = rstat
        if not os.path.isdir(root):
            rstat["missing"] = True
            continue
        pending.append((collections.deque([(root.rstrip("/") or "/", 0)]), rstat))

    def walk(queue, rstat, stop):
        while queue:
            if clock() > stop:
                return
            d, depth = queue.popleft()
            try:
                with os.scandir(d) as it:
                    entries = list(it)
            except OSError:
                stats["unreadable"] += 1
                continue
            stats["dirs_walked"] += 1
            rstat["dirs_walked"] += 1
            names = {e.name for e in entries}
            if names & SEARCH_HINTS or any(SKILL_LOG.match(n) for n in names):
                found.append(d)
                rstat["found"] += 1
                continue
            if REC_JSON in names:
                sessions.append(d)
            if depth >= max_depth:
                continue
            for e in entries:
                if e.name in PRUNE_NAMES or e.name.lower().endswith(PRUNE_SUFFIX):
                    continue
                try:
                    if not e.is_dir(follow_symlinks=False):
                        continue
                except OSError:
                    continue
                if any(x in e.path + "/" for x in excludes):
                    continue
                queue.append((e.path, depth + 1))

    while pending and clock() - t0 < budget_s:
        rnd = list(pending)
        for i, (queue, rstat) in enumerate(rnd):
            left = budget_s - (clock() - t0)
            if left <= 0:
                break
            walk(queue, rstat, clock() + left / (len(rnd) - i))
            if not queue:
                pending.remove((queue, rstat))
    for queue, rstat in pending:
        rstat["truncated"] = stats["truncated"] = True
        rstat["left_in_queue"] = len(queue)
    stats["elapsed_s"] = round(clock() - t0, 1)
    return found, sessions, stats


def session_search_dirs(session):
    """Search dirs a checkpoint.py session points at: its `report` (the expected search output)
    lives in the search's out dir, which in hive_remote mode is not under the session at all."""
    try:
        with open(os.path.join(session, REC_JSON)) as fh:
            rec = json.load(fh)
    except (OSError, ValueError):
        return []
    rep = rec.get("report") if isinstance(rec, dict) else None
    if not rep:
        return []
    d = rep if os.path.isdir(rep) else os.path.dirname(rep)
    return [d] if os.path.isdir(d) else []


def core_members():
    """Non-teaching members of proteomics-grp, or None when the group cannot be read."""
    try:
        import grp                                                          # noqa: PLC0415
        return {u for u in grp.getgrnam(GROUP_NAME).gr_mem if not TEACHING_ACCOUNT.search(u)}
    except (ImportError, KeyError, OSError):
        return None


def member_home_roots(members):
    """~/proteomics-pipeline of the current user and of every non-teaching group member."""
    roots = [os.path.expanduser("~/proteomics-pipeline")]
    try:
        import pwd                                                          # noqa: PLC0415
        for u in sorted(members or ()):
            try:
                roots.append(os.path.join(pwd.getpwnam(u).pw_dir, "proteomics-pipeline"))
            except KeyError:
                continue
    except ImportError:
        pass
    return list(dict.fromkeys(roots))


def _owner(path):
    try:
        import pwd                                                          # noqa: PLC0415
        return pwd.getpwuid(os.stat(path).st_uid).pw_name
    except (ImportError, KeyError, OSError):
        return None


def _teaching_reason(out, user=None):
    """Why a search is coursework, or None: the account running the skill, or the account that
    owns the search, is a teaching account (TEACHING_ACCOUNT). ONE rule for stage and backfill."""
    for who, acct in (("run by", user), ("owned by", _owner(os.path.realpath(out)))):
        if acct and TEACHING_ACCOUNT.search(acct):
            return (f"{who} teaching account {acct}: coursework is never handed to the Core "
                    f"corpus (FRAN keeps teaching data out)")
    return None


def _in_core_group():
    """True / False when this account's proteomics-grp membership can be read, else None."""
    try:
        import grp                                                          # noqa: PLC0415
        g = grp.getgrnam(GROUP_NAME)
    except (ImportError, KeyError, OSError):
        return None
    return (g.gr_gid in os.getgroups() or os.getgid() == g.gr_gid
            or getpass.getuser() in g.gr_mem)


def core_search(out, members, prefixes=CORE_PREFIXES):
    """(is_core, why). Coursework never is (_teaching_reason). Otherwise a search is the Core's if
    it lives in a Core tree, or if a non-teaching member of proteomics-grp owns it. Anything else
    is a collaborator's and is never handed over."""
    real = os.path.realpath(out)
    teach = _teaching_reason(out)
    if teach:
        return False, teach
    if any(real.startswith(p) for p in prefixes):
        return True, "in a Core tree"
    owner = _owner(real)
    if owner and members and owner in members:
        return True, f"owned by Core member {owner}"
    return False, (f"outside the Core trees and owned by {owner or 'an unknown account'}, who is "
                   f"not a Core member of {GROUP_NAME}")


GENERIC_DIRS = {"search", "search_out", "output", "out", "results", "search_results", "diann"}


def derive_name(out):
    """A readable search name when nobody typed one: the folder, or for a generic folder
    (`<session>/output/search`, `<job>/search_out`) the nearest ancestor that says what it was."""
    parts = [p for p in os.path.realpath(out).split("/") if p]
    for p in reversed(parts):
        if p.lower() not in GENERIC_DIRS:
            return p
    return parts[-1] if parts else None


class _StageArgs:
    """argparse-shaped arguments for check()/do_stage() on one backfill candidate."""
    def __init__(self, out, dry_run, name=None):
        self.out, self.dry_run, self.name = out, dry_run, name
        self.name_source = "derived from the folder name by `fran_deposit.py backfill`"
        self.skip = self.force = self.qc = self.not_qc = False
        self.organism = self.taxon = self.fasta_meta = None
        self.require_completion_marker = True     # unattended: no marker, no stage


def plan_backfill(dirs, apply=False, members=None, prefixes=CORE_PREFIXES, runs=None):
    """Classify each candidate directory; with apply=True, stage the eligible ones.

    Decision order, cheapest and most decisive first: a search this skill did not run, a search
    that is not the Core's, a recorded decision (opted out / not Core), an entry already in the
    drop dir, a search the cron's logs show FRAN already ingested (`runs`, from load_ingest_logs),
    and then check() -- the same gate `stage` uses, with its reason codes.

    The log check matters: searches have reached FRAN by other routes than the drop dir. On
    2026-09-08 FRAN's own queue ingested PROT_0793/search_mouse and search_hela by their real paths,
    and neither has a receipt -- staging them again would only put two duplicates in front of the
    cron's duplicate guard."""
    drop = os.environ.get("FRAN_DROP_DIR", DROP_DIR)
    seen, rows = set(), []
    for d in dirs:
        real = os.path.realpath(d)
        if real in seen:
            continue
        seen.add(real)
        row = {"out": d}
        try:
            names = set(os.listdir(d))
        except OSError as e:
            rows.append({**row, "decision": "skip", "reason": "not_on_hive",
                         "detail": f"cannot list {d}: {e.strerror or e}"})
            continue
        marker = skill_marker(names)
        if not marker:
            rows.append({**row, "decision": "skip", "reason": "not_a_skill_search",
                         "detail": "no file this skill writes (search_provenance.json, "
                                   "fran_deposit.json, s5_report_<job>.log, ...)"})
            continue
        row["marker"] = marker
        # QC first: the same rule as stage (is_qc_run), with every name this search has been given --
        # the folder, and the analysis name a previous stage recorded. An explicit --qc / --not-qc
        # recorded by an earlier stage wins over the rule, both ways.
        entry = os.path.join(drop, entry_name(d))
        prior = read_receipt(d) or {}
        man = {}
        if os.path.isdir(entry):
            try:
                with open(os.path.join(entry, MANIFEST)) as fh:
                    man = json.load(fh)
            except (OSError, ValueError):
                pass
        names = [n for n in (derive_name(d), man.get("search_name"), prior.get("search_name")) if n]
        qc, qwhy = decide_qc(d, names=names, receipt=prior, manifest=man)
        if qc:
            qrow = {**row, "decision": "excluded", "reason": "qc_run", "detail": qwhy}
            if os.path.isdir(entry) and not (man.get("qc") is True and man.get("exclude") is True):
                qrow["staged_entry"] = entry
                c = {"search_dir": d, "reason": "qc_run", "qc_rule": qwhy, "entry": entry,
                     "user": getpass.getuser(), "detail": f"QC run ({qwhy})"}
                if apply:
                    _record_decision(_StageArgs(d, dry_run=False), c)
                    qrow["withdrawn"] = c.get("withdrawn")
                else:
                    qrow["would_withdraw"] = entry
            rows.append(qrow)
            continue
        core, why = core_search(d, members, prefixes)
        if not core:
            rows.append({**row, "decision": "skip", "reason": "not_core_facility", "detail": why})
            continue
        if prior.get("status") in DECISION_STATUS:
            rows.append({**row, "decision": "excluded" if prior["status"] == "qc_run" else "skip",
                         "reason": prior["status"],
                         "detail": f"recorded by {prior.get('decided_by') or '?'} at "
                                   f"{prior.get('at') or '?'}: {prior.get('detail') or ''}".strip()})
            continue
        if os.path.isdir(entry):
            rows.append({**row, "decision": "skip", "reason": "already_staged",
                         "detail": f"{entry} already exists", "entry": entry})
            continue
        if runs:
            st = entry_log_state(runs, entry, real)
            if st["state"] == "ingested":
                rows.append({**row, "decision": "skip", "reason": "already_ingested",
                             "detail": f"the cron's log records it {st['outcome']} at "
                                       f"{st['when']} ({st['log']})"})
                continue
        a = _StageArgs(d, dry_run=not apply, name=derive_name(d))
        try:
            res = do_stage(a)
        except Exception as e:                                              # noqa: BLE001
            rows.append({**row, "decision": "error", "reason": "stage_error",
                         "detail": f"{type(e).__name__}: {e}"})
            continue
        keep = {k: res.get(k) for k in ("engine", "engine_version", "organism", "taxon",
                                        "organism_source", "fasta_path", "name", "entry",
                                        "report_bytes", "link_errors") if res.get(k) is not None}
        keep["xic_files"] = (res.get("xic") or {}).get("n_files")
        try:
            keep["report_mtime"] = _when(os.path.getmtime(res["report"])) if res.get("report") else None
        except OSError:
            pass
        if res.get("staged"):
            rows.append({**row, **keep, "decision": "staged", "reason": "ok"})
        elif res.get("reason") == "needs_agent_check":
            rows.append({**row, **keep, "decision": "skip", "reason": "needs_agent_check",
                         "detail": res.get("detail")})
        elif res.get("eligible"):
            rows.append({**row, **keep, "decision": "would_stage", "reason": "ok"})
        else:
            rows.append({**row, **keep, "decision": "skip", "reason": res.get("reason"),
                         "detail": res.get("detail")})
    return rows


def login_guard(a, n_list, env=None):
    """A refusal dict when this would walk NFS on a login node, else None."""
    env = os.environ if env is None else env
    if env.get("SLURM_JOB_ID") or a.allow_login_node:
        return None
    walking = bool(a.roots) or not a.list
    if not walking and n_list <= LOGIN_LIST_MAX:
        return None
    what = (f"walk {', '.join(a.roots or BACKFILL_ROOTS)}" if walking
            else f"check {n_list} directories")
    return {"action": "backfill", "refused": "login_node",
            "detail": f"This would {what} over NFS, and that is not allowed on a HIVE login node. "
                      f"Write the job and submit it: python3 {os.path.abspath(__file__)} backfill "
                      f"--sbatch (then `sbatch <the file it prints>`). A --list of up to "
                      f"{LOGIN_LIST_MAX} directories is fine here."}


def _queue_lines(partition=None, account=None, qos=None):
    """#SBATCH queue lines from the skill's ONE queue rule, run_search.slurm_queue()."""
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import run_search                                                   # noqa: PLC0415
        part, acct, q = run_search.slurm_queue(partition, account, qos, preemptible_ok=True)
        try:
            import diann_parallel                                           # noqa: PLC0415
            req = diann_parallel.needs_requeue(part, q)
        except Exception:                                                   # noqa: BLE001
            req = False
        src = "run_search.slurm_queue()"
    except (Exception, SystemExit) as e:                                    # noqa: BLE001
        part, acct, q, req = partition, account, qos, False
        src = f"queue not detected ({e}); the values given are used as-is"
    lines = [f"#SBATCH --{k}={v}" for k, v in (("partition", part), ("account", acct),
                                                ("qos", q)) if v]
    if req:
        lines.append("#SBATCH --requeue")
    return lines, src


def write_backfill_sbatch(a):
    d = os.path.abspath(os.path.expanduser(a.sbatch_dir))
    os.makedirs(d, exist_ok=True)
    stamp = time.strftime("%Y%m%d_%H%M%S")
    qlines, qsrc = _queue_lines(a.partition, a.account, a.qos)
    minutes = a.sbatch_minutes or (math.ceil(a.time_budget / 60) + 10)
    args = ["backfill", "--time-budget", str(a.time_budget), "--max-depth", str(a.max_depth)]
    for r in a.roots or []:
        args += ["--roots", os.path.abspath(r)]
    if a.list:
        args += ["--list", os.path.abspath(a.list)]
    if a.no_homes:
        args.append("--no-homes")
    if a.apply:
        args.append("--apply")
    cmd = " ".join(shlex.quote(x) for x in ["python3", os.path.abspath(__file__)] + args)
    report = os.path.join(d, "fran_backfill_${SLURM_JOB_ID}.json")
    path = os.path.join(d, f"fran_backfill_{stamp}.sbatch")
    script = "\n".join([
        "#!/bin/bash -l",
        "#SBATCH --job-name=fran_backfill",
        f"#SBATCH --output={d}/fran_backfill_%j.log",
        "#SBATCH --cpus-per-task=1",
        "#SBATCH --mem=2G",
        f"#SBATCH --time={minutes}",
        *qlines,
        f"# fran_backfill -- written by fran_deposit.py backfill --sbatch on {stamp}.",
        f"# {'STAGES eligible searches (--apply).' if a.apply else 'DRY RUN: stages nothing.'}",
        f"# Queue: {qsrc}. Override on the command line: sbatch --partition=... --account=...",
        "set -uo pipefail",
        f'{cmd} --report {shlex.quote(d)}/"fran_backfill_${{SLURM_JOB_ID}}.json"',
        ""])
    with open(path, "w") as fh:
        fh.write(script)
    os.chmod(path, 0o755)
    return {"action": "backfill", "sbatch": path, "mode": "apply" if a.apply else "dry_run",
            "report_will_be": report.replace("${SLURM_JOB_ID}", "<jobid>"),
            "queue": qlines, "queue_source": qsrc,
            "submit": f"sbatch {shlex.quote(path)}"}


def backfill(a):
    if a.sbatch:
        jout(write_backfill_sbatch(a))
    listed = []
    if a.list:
        try:
            with open(a.list) as fh:
                listed = [x.strip() for x in fh if x.strip() and not x.startswith("#")]
        except OSError as e:
            jout({"action": "backfill", "error": f"cannot read --list {a.list}: {e}"}, 2)
    refusal = login_guard(a, len(listed))
    if refusal:
        jout(refusal, 2)

    members = core_members()
    walk_roots = []
    if a.roots or not a.list:
        walk_roots = list(a.roots or BACKFILL_ROOTS)
        if not a.no_homes and not a.roots:
            walk_roots += member_home_roots(members)
    found, sessions, stats = discover(walk_roots, a.max_depth, a.time_budget) if walk_roots \
        else ([], [], {"dirs_walked": 0, "truncated": False})
    extra = [x for s in sessions for x in session_search_dirs(s)]
    runs, info = load_ingest_logs()
    rows = plan_backfill(listed + found + extra, apply=a.apply, members=members, runs=runs)

    by = collections.Counter(r["reason"] for r in rows if r["decision"] == "skip")
    res = {"action": "backfill", "mode": "apply" if a.apply else "dry_run",
           "host": os.uname().nodename, "slurm_job": os.environ.get("SLURM_JOB_ID"),
           "user": getpass.getuser(), "roots": walk_roots, "listed": len(listed),
           "walk": stats, "sessions_seen": len(sessions),
           "n_search_dirs": len(rows),
           "n_would_stage": sum(r["decision"] == "would_stage" for r in rows),
           "n_staged": sum(r["decision"] == "staged" for r in rows),
           "n_errors": sum(r["decision"] == "error" for r in rows),
           "skipped_by_reason": dict(sorted(by.items())),
           "would_stage": [r for r in rows if r["decision"] == "would_stage"],
           "staged": [r for r in rows if r["decision"] == "staged"],
           "errors": [r for r in rows if r["decision"] == "error"],
           # QC runs are listed on their own: Brett's rule, not a defect of the search
           "n_excluded_qc_run": sum(r["decision"] == "excluded" for r in rows),
           "excluded_qc_run": [r for r in rows if r["decision"] == "excluded"],
           "skipped": [r for r in rows if r["decision"] == "skip"]}
    if stats.get("truncated"):
        res["note_truncated"] = ("the walk hit its time budget before finishing; raise "
                                 "--time-budget (and the job's --time) or pass --roots to cover "
                                 "the rest")
    # The cron's state, without the network half: staging into a stuck cron adds to its queue.
    try:
        p = progress_health(runs, submit=submit_log_info(), info=info)
        res["fran_cron"] = {"verdict": p["verdict"], "detail": p.get("detail")}
    except Exception as e:                                                  # noqa: BLE001
        res["fran_cron"] = {"verdict": "unknown", "detail": f"{type(e).__name__}: {e}"}
    if not a.apply and res["n_would_stage"]:
        res["next"] = ("re-run with --apply to stage these (or `backfill --sbatch --apply`); "
                       "each becomes one entry in the drop dir, named deterministically")
    if a.report:
        try:
            with open(a.report, "w") as fh:
                json.dump(res, fh, indent=2)
            res["report"] = a.report
        except OSError as e:
            res["report_error"] = str(e)
    jout(res)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["check", "stage", "verify", "health", "backfill"])
    ap.add_argument("--out", help="the search output directory (on HIVE); check/stage/verify")
    ap.add_argument("--name", default=None, help="corpus search name (default: from the raw files)")
    ap.add_argument("--organism", default=None,
                    help="e.g. 'Homo sapiens' -- the organism the user CONFIRMED")
    ap.add_argument("--taxon", type=int, default=None)
    ap.add_argument("--fasta-meta", default=None, help="<fasta>.meta.json from fetch_fasta.py")
    ap.add_argument("--force", action="store_true", help="re-stage even if a receipt exists")
    ap.add_argument("--skip", action="store_true", help="opt this run out of FRAN")
    ap.add_argument("--dry-run", action="store_true", help="show what would be linked, link nothing")
    q = ap.add_mutually_exclusive_group()
    q.add_argument("--qc", action="store_true",
                   help="this IS a QC run: never hand it to FRAN (withdraws it if already staged)")
    q.add_argument("--not-qc", action="store_true",
                   help="this is NOT a QC run, whatever its name says (a false positive of the rule)")
    h = ap.add_argument_group("health")
    h.add_argument("--alert", action="store_true",
                   help="post a short Slack alert when not healthy (needs notify_slack.py)")
    h.add_argument("--no-code", action="store_true",
                   help="skip the GitHub comparison of FRAN's ingest code (no network)")
    b = ap.add_argument_group("backfill")
    b.add_argument("--roots", action="append", default=None,
                   help=f"walk these instead of the Core trees (default: {', '.join(BACKFILL_ROOTS)} "
                        f"and members' ~/proteomics-pipeline); repeatable")
    b.add_argument("--list", default=None, help="file of search out dirs, one per line")
    b.add_argument("--apply", action="store_true", help="stage what is eligible (default: dry run)")
    b.add_argument("--max-depth", type=int, default=WALK_MAX_DEPTH)
    b.add_argument("--time-budget", type=int, default=WALK_BUDGET_S, help="seconds for the walk")
    b.add_argument("--no-homes", action="store_true", help="do not walk members' ~/proteomics-pipeline")
    b.add_argument("--report", default=None, help="also write the JSON result here")
    b.add_argument("--allow-login-node", action="store_true",
                   help="walk even outside SLURM (don't, on a HIVE login node)")
    b.add_argument("--sbatch", action="store_true", help="write the SLURM job for this backfill")
    b.add_argument("--sbatch-dir", default="~/fran_backfill")
    b.add_argument("--sbatch-minutes", type=int, default=None)
    b.add_argument("--partition", default=None)
    b.add_argument("--account", default=None)
    b.add_argument("--qos", default=None)
    a = ap.parse_args()
    if a.action in ("check", "stage", "verify") and not a.out:
        ap.error(f"{a.action} needs --out <search out dir>")
    {"check": lambda: jout(check(a)), "stage": lambda: stage(a), "verify": lambda: verify(a),
     "health": lambda: health(a), "backfill": lambda: backfill(a)}[a.action]()


if __name__ == "__main__":
    main()
