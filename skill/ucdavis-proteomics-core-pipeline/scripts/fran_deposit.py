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
  # 3. later: did the cron actually take it?
  python3 fran_deposit.py verify --out <search out dir>

All three run **on HIVE** (in hive_remote mode through `hive_exec.sh`, like every other
HIVE-side step). Each prints one JSON object.

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

Set `FRAN_DEPOSIT=off` (or pass `--skip`) to keep a run out of FRAN.

## Pass the organism through — it is the one thing only the skill knows

A DIA-NN `report.parquet` carries **no organism column**, so a search ingested without one
sits in the corpus with `organism` NULL and is invisible on FRAN's species page. The skill
already had the user *confirm* the organism, and `fetch_fasta.py` wrote it to
`<fasta>.meta.json` — read automatically here, overridable with `--organism`/`--taxon`, and
written into `fran_manifest.json` for the cron to pass to `corpus_ingest.py`. Never invented:
an unresolved organism stays absent from the manifest rather than being guessed.

Overrides: FRAN_DROP_DIR, FRAN_DEPOSIT=off, DELIMP_PG_TOKEN_FILE (verify only).
"""
import argparse
import getpass
import glob
import hashlib
import json
import os
import re
import shlex
import subprocess
import sys

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
# to it. Mirrors FRAN's own ENGINE_MARKERS so both ends agree on what a directory is.
DETECT_MARKERS = [
    ("fragpipe", ("dia-quant-output/report.tsv", "fragpipe.fp-manifest")),
    ("radiant",  ("radiant_results/fulcrum-results", "fulcrum-results/_SUCCESS",
                  "delimp_report.parquet")),
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


def find_report(out, engine):
    for c in ENGINE_REPORTS.get(engine, ()):
        p = os.path.join(out, c)
        if os.path.exists(p):
            return p
    return None


def organism_from_meta(out, explicit_meta=None):
    """(organism, taxid, source). DIA-NN's report carries NO organism column, so unless we pass
    it the corpus row is left NULL -- the species page then simply cannot see this search. The
    skill already had the user CONFIRM the organism at step 3, and fetch_fasta.py wrote it to
    <fasta>.meta.json, so pass that through rather than leaving a hole. Never guessed."""
    cands = [explicit_meta] if explicit_meta else []
    cands += sorted(glob.glob(os.path.join(out, "*.meta.json"))) \
        + sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(out)), "*.fasta.meta.json"))) \
        + sorted(glob.glob(os.path.join(out, "..", "input", "*.fasta.meta.json")))
    for c in cands:
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
        return p
    except OSError as e:
        # A receipt we could not write is a resume hazard, not a failure of the deposit -- report
        # it instead of raising, so the caller learns the deposit itself still stands.
        data["receipt_error"] = str(e)
        return None


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
    if not os.path.isdir(drop):
        try:
            os.makedirs(drop, mode=0o2775, exist_ok=True)
        except OSError as e:
            r.update(reason="no_drop_dir",
                     detail=f"FRAN's drop directory {drop} does not exist and could not be created "
                            f"({e}). Create it once: mkdir -p {drop} && chmod 2775 {drop}")
            return r
    # THE CORE GATE, and it is a real permission rather than a flag: the drop directory lives
    # inside /quobyte/proteomics-grp, so a HIVE account outside proteomics-grp cannot write here.
    # A collaborator's search physically cannot be staged.
    if not os.access(drop, os.W_OK | os.X_OK):
        r.update(reason="not_core_facility",
                 detail=f"{drop} is not writable by {r['user']}, so this is not a Proteomics Core "
                        f"run. Collaborator searches are never handed to the Core corpus.")
        return r

    org, tax, org_src = organism_from_meta(out, a.fasta_meta)
    r["organism"] = a.organism or org
    r["taxon"] = a.taxon or tax
    r["organism_source"] = "--organism (given)" if a.organism else org_src
    r["name"] = a.name

    # XIC chromatograms, if the search asked DIA-NN for them (`--xic` in the cfg; diann_parallel.py
    # puts it on step 4 only). Reported either way -- their absence is a fact about the search, not
    # an error, and FRAN's XIC lane simply has nothing to ingest for this one.
    xics = find_xic_dirs(out)
    r["xic"] = {"present": bool(xics), "n_files": sum(x["n_files"] for x in xics), "dirs": xics}
    r["links"] = [i for i in LINK_ITEMS if os.path.exists(os.path.join(out, i))]

    r.update(eligible=True, reason="ok")
    r["next_command"] = f"python3 {os.path.abspath(__file__)} stage --out {shlex.quote(out)}"
    return r


def stage(a):
    """Hand the search to FRAN by linking it into the drop directory. Nothing is copied, nothing is
    ingested here — FRAN's cron does the ingest when it next scans."""
    c = check(a)
    if not c["eligible"]:
        jout({**c, "staged": False}, 0)
    out, entry = c["search_dir"], c["entry"]
    if a.dry_run:
        jout({**c, "staged": False, "dry_run": True,
              "would_link": c["links"], "would_write": os.path.join(entry, MANIFEST)})

    # Re-staging must converge, not accumulate: an entry from an earlier attempt may hold links to
    # files the search has since replaced (a resumed chain rewrites report.parquet). Relink from
    # scratch rather than leaving a stale mixture of both runs.
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
        # Absent, not guessed, when unknown -- a NULL organism is honest, an invented one is a
        # species claim. DIA-NN reports carry no organism column, so this is the only source.
        "organism": c.get("organism"),
        "taxon": c.get("taxon"),
        "organism_source": c.get("organism_source"),
        "xic": c["xic"],
        "linked": linked,
        "staged_by": c["user"],
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
    with open(os.path.join(entry, MANIFEST), "w") as fh:
        json.dump(manifest, fh, indent=2)

    receipt = {"status": "staged", "search_dir": out, "entry": entry, "engine": c["engine"],
               "organism": c.get("organism"), "taxon": c.get("taxon"),
               "xic": c["xic"], "linked": linked, "staged_by": c["user"], "fran": FRAN_URL}
    if c.get("link_errors"):
        receipt["link_errors"] = c["link_errors"]
    write_receipt(out, receipt)
    jout({**c, "staged": True, "entry": entry, "linked": linked,
          "manifest": os.path.join(entry, MANIFEST),
          "receipt": os.path.join(out, RECEIPT),
          "note": "FRAN's ingest cron picks this up on its next scan; nothing was copied.",
          "then": f"python3 {os.path.abspath(__file__)} verify --out {shlex.quote(out)}"})


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
               FROM delimp_searches WHERE output_dir=%s ORDER BY submitted_at DESC""", (out,))
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


def corpus_query(out):
    """Ask the corpus whether this search is in it. Returns None when this account has no way to
    ask -- which is the COMMON case now that staging needs no credential, and is not an error."""
    ing = first_readable(INGEST_DIRS, isdir=True)
    py = first_readable(PY_CANDIDATES, executable=True)
    tok = first_readable(TOKEN_CANDIDATES)
    if not (ing and py and tok):
        return None
    env = {**os.environ, "ING": ing, "OUTDIR": os.path.realpath(out), "DELIMP_PG_TOKEN_FILE": tok}
    try:
        p = subprocess.run([py, "-c", VERIFY_PY], capture_output=True, text=True,
                           env=env, timeout=600)
        if p.returncode != 0:
            return {"error": (p.stderr or p.stdout).strip()[-400:]}
        return json.loads(p.stdout.strip().splitlines()[-1])
    except (OSError, subprocess.SubprocessError, ValueError, IndexError) as e:
        return {"error": str(e)[:200]}


def verify(a):
    """Did the handover actually land? Three questions, in order of how much they prove:
      1. is the entry still in the drop directory, with a live link to a real report?
      2. has the cron marked it ingested?
      3. (only if this account can read a corpus token) does the corpus actually hold the rows?
    A staged-but-not-yet-ingested search is the NORMAL state right after a run — report it as
    pending, never as a failure."""
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
        mk = os.path.join(entry, INGESTED_MARKER)
        if os.path.isfile(mk):
            try:
                with open(mk) as fh:
                    r["cron_marker"] = json.load(fh)
            except (OSError, ValueError):
                r["cron_marker"] = {"malformed": True}

    corpus = corpus_query(out)
    r["corpus"] = corpus
    in_corpus = bool(corpus and corpus.get("in_corpus"))
    r["ingested"] = bool(in_corpus or r.get("cron_marker"))
    r["state"] = ("ingested" if r["ingested"] else
                  "staged_pending_cron" if r["staged"] else
                  "not_staged")
    r["detail"] = {
        "ingested": "FRAN has this search.",
        "staged_pending_cron": "Handed over; FRAN's ingest cron takes it on its next scan. This is "
                               "the expected state immediately after a run — not a failure."
                               + ("" if corpus is not None else
                                  " (This account cannot read a corpus token, so only the drop "
                                  "entry was checked, not the database.)"),
        "not_staged": f"No entry at {entry}. Run `stage` first, or check why `check` refused.",
    }[r["state"]]
    if r["staged"] and r.get("broken_links"):
        r["detail"] += (f" WARNING: {len(r['broken_links'])} link(s) point at files that no longer "
                        f"exist — the cron will skip this entry. Re-run `stage --force`.")

    receipt = read_receipt(out) or {"search_dir": out}
    receipt["status"] = "ingested" if r["ingested"] else receipt.get("status", "staged")
    receipt["verified"] = {k: r[k] for k in ("state", "entry", "ingested") if k in r}
    if corpus and corpus.get("search"):
        receipt["search_id"] = corpus["search"]["search_id"]
    r["receipt"] = write_receipt(out, receipt)
    jout(r)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("action", choices=["check", "stage", "verify"])
    ap.add_argument("--out", required=True, help="the search output directory (on HIVE)")
    ap.add_argument("--name", default=None, help="corpus search name (default: from the raw files)")
    ap.add_argument("--organism", default=None,
                    help="e.g. 'Homo sapiens' -- the organism the user CONFIRMED")
    ap.add_argument("--taxon", type=int, default=None)
    ap.add_argument("--fasta-meta", default=None, help="<fasta>.meta.json from fetch_fasta.py")
    ap.add_argument("--force", action="store_true", help="re-stage even if a receipt exists")
    ap.add_argument("--skip", action="store_true", help="opt this run out of FRAN")
    ap.add_argument("--dry-run", action="store_true", help="show what would be linked, link nothing")
    a = ap.parse_args()
    {"check": lambda: jout(check(a)), "stage": lambda: stage(a), "verify": lambda: verify(a)}[a.action]()


if __name__ == "__main__":
    main()
