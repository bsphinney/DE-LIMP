#!/usr/bin/env python3
"""
fragpipe_db.py -- THE database a FragPipe search reads, and its decoys (2.15, review H3).

MSFragger refuses a FASTA without decoys ("No decoys found in the FASTA file"), and the skill's
FASTA (fetch_fasta.py) is targets only -- so is the one a bundle's reproduce.sh rebuilds in its
step 4. Before 2.15 the user had to add decoys by hand with Philosopher and give the .fas, and a
replay of a FragPipe search could not run at all. Now:

  * --fasta with a COMPLETE decoy set (every target `<id>` has a `<tag><id>` entry and every decoy
    its target): searched as it is;
  * --fasta with NO decoys: the search's job builds them first, with the Philosopher of the
    FragPipe install it runs (`database --custom <target> --prefix <tag>`, never --contam: the
    skill's contaminants are already in the FASTA, tagged Cont_), into
    <out>/fragpipe_db/<YYYY-MM-DD>-decoys-<target name>.fas -- Philosopher's own name, the date
    the search was generated -- and the workflow FragPipe runs names that file. Philosopher
    5.1.3's build is deterministic: two builds of a Core submission's 2,403-entry FASTA, three
    days apart, gave the same sha256 (2199868e..., HIVE, 2026-11-04);
  * --fasta with SOME decoys but not a complete set (decoys without their target, targets without
    their decoy): REFUSED before anything is written -- FragPipe would search it, and the
    target-decoy FDR would be estimated on half a decoy set.

The decoy tag is the workflow's database.decoy-tag, else FragPipe's own default (rev_), and the
record says which (a default, said so).

The sha256 of what was searched is not computed on the login node: it is the fetch_fasta.py
sidecar's when that sidecar describes --fasta itself, and in every case the job measures the
database MSFragger reads (sha256, entries, decoys; a half set fails the job) into
<out>/fragpipe_job.json, which run_search.py --adapt-only merges into search_provenance.json
`fragpipe_database`. Stdlib only.

  python3 fragpipe_db.py measure --db <searched db> --tag rev_ --record <out>/fragpipe_job.json \\
      [--target <target FASTA>] [--java-text "<java -version>"]
"""
import argparse
import datetime
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scratch_files  # noqa: E402
from safe_write import sha256_file as _sha256  # noqa: E402  the ONE chunked sha256
# FragPipe's Java floor and the reading of `java -version`: defined once, in engine_choice.py
# (detect_env.sh, acquire_tools.sh, resolve_defaults and this job record all read it)
from engine_choice import FRAGPIPE_JAVA_MIN as MIN_JAVA, java_major  # noqa: E402,F401

DEFAULT_DECOY_TAG = "rev_"                 # FragPipe's own default (its workflows set rev_)
DECOY_TAG_KEY = "database.decoy-tag"
DB_PATH_KEY = "database.db-path"
DB_DIR = scratch_files.FRAGPIPE_DB_DIR     # <out>/fragpipe_db: scratch, never zipped or delivered
JOB_RECORD = "fragpipe_job.json"           # what the job measured (database, java), in <out>
# the sha256 of a database the job will build, until the job has measured it (review L1)
NOT_YET_MEASURED = ("not yet measured: the job measures the database it builds (fragpipe_job.json) "
                    "before FragPipe starts")
TAG_FROM_WORKFLOW = f"workflow ({DECOY_TAG_KEY})"
TAG_DEFAULT = (f"default: FragPipe's own {DEFAULT_DECOY_TAG} (the workflow sets no "
               f"{DECOY_TAG_KEY})")
REFUSED = "[run_search] REFUSED -- nothing was written or submitted."


def workflow_value(params, key):
    """A FragPipe workflow's value for `key` as FragPipe reads it (Java .properties: a backslash
    escapes the next character), or None when it is absent or empty -- THE reading of one key
    (run_search._workflow_value is this)."""
    try:
        with open(params, encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                if ln.startswith(key + "="):
                    v = re.sub(r"\\(.)", r"\1", ln.split("=", 1)[1].strip())
                    return v or None
    except (OSError, TypeError):
        pass
    return None


def decoy_tag(params):
    """(the decoy tag the search's workflow uses, where it came from): the workflow's
    database.decoy-tag, else FragPipe's own rev_ -- tagged as that default (L11 note)."""
    tag = workflow_value(params, DECOY_TAG_KEY) if params else None
    return (tag, TAG_FROM_WORKFLOW) if tag else (DEFAULT_DECOY_TAG, TAG_DEFAULT)


def _entry_id(header, tag=""):
    """A FASTA header line's identifier (its first word, without '>' and without `tag`)."""
    word = header[1:].strip().split(None, 1)
    word = word[0] if word else ""
    return word[len(tag):] if tag and word.startswith(tag) else word


def inspect(path, tag):
    """A FASTA's decoy state, from its headers only (a proteome's sequences are never held):
    {"n_entries", "n_targets", "n_decoys", "state", "targets_without_decoy",
    "decoys_without_target" (the first few of each)}. state: "target_only" (no entry starts with
    `tag`), "target_decoy" (every target `<id>` has `<tag><id>` and every decoy its target),
    "partial" (anything else: decoys, but not one for each target). OSError propagates."""
    targets, decoys = [], []
    import gzip
    opn = gzip.open if str(path).endswith(".gz") else open
    with opn(path, "rt", errors="replace") as fh:
        for ln in fh:
            if ln.startswith(">"):
                (decoys if ln[1:].lstrip().startswith(tag) else targets).append(
                    _entry_id(ln, tag if ln[1:].lstrip().startswith(tag) else ""))
    t, d = set(targets), set(decoys)
    state = ("target_only" if not decoys else
             "target_decoy" if t == d and len(targets) == len(decoys) else "partial")
    return {"n_entries": len(targets) + len(decoys), "n_targets": len(targets),
            "n_decoys": len(decoys), "state": state,
            "targets_without_decoy": sorted(t - d)[:5], "decoys_without_target": sorted(d - t)[:5]}


def planned_database(out, fasta, today=None):
    """Where the job builds a target-only FASTA's target+decoy database: Philosopher's own name,
    <YYYY-MM-DD>-decoys-<target file name>.fas (fetch_fasta.PHILOSOPHER_DECOY_FAS), dated the day
    the search was generated, in <out>/fragpipe_db/."""
    day = (today or datetime.date.today()).isoformat()
    return os.path.join(os.path.abspath(out), DB_DIR, f"{day}-decoys-{os.path.basename(fasta)}.fas")


def sidecar_sha256(fasta, sidecar):
    """(sha256, where it came from) of --fasta WITHOUT reading it: the fetch_fasta.py sidecar's
    `sha256` when that sidecar describes this very file (its `fasta` field, or it is
    <fasta>.meta.json), else (None, why not). The job measures the file either way."""
    if not sidecar:
        return None, "no fetch_fasta.py sidecar describes --fasta: the job measures it"
    try:
        with open(sidecar) as fh:
            meta = json.load(fh)
    except (OSError, ValueError) as e:
        return None, f"{sidecar} could not be read ({e}): the job measures --fasta"
    own = os.path.abspath(sidecar) == os.path.abspath(fasta) + ".meta.json"
    named = isinstance(meta, dict) and meta.get("fasta") and \
        os.path.realpath(str(meta["fasta"])) == os.path.realpath(fasta)
    sha = meta.get("sha256") if isinstance(meta, dict) else None
    if sha and (own or named):
        return sha, (f"the fetch_fasta.py sidecar {sidecar} (not re-hashed when the search was "
                     f"generated; the job measures the database it searches)")
    return None, (f"{sidecar} records no sha256 for {fasta}: the job measures it"
                  if (own or named) else
                  f"{sidecar} describes another file (the target FASTA of this .fas): the job "
                  f"measures --fasta")


def database_record(params, fasta, out, philosopher=None, sidecar=None, today=None):
    """The record of the database a FragPipe search reads (search_provenance.json
    `fragpipe_database`), decided before anything is written; sys.exit (REFUSED) when --fasta
    cannot be read, holds half a decoy set, names another database than the workflow's, or needs
    decoys and no Philosopher is known. `philosopher`: tools.json fragpipe_tools.philosopher
    ({"path", "version"}) -- the one the FragPipe install runs."""
    fas = os.path.abspath(fasta)
    given = workflow_value(params, DB_PATH_KEY)
    if given and os.path.realpath(os.path.abspath(given)) != os.path.realpath(fas):
        sys.exit(f"{REFUSED} The workflow {params} names {DB_PATH_KEY}={given}, but --fasta is "
                 f"{fas}: FragPipe would search the workflow's database, not the one recorded for "
                 f"this search. Give --fasta {given}, or a workflow without {DB_PATH_KEY} "
                 f"(run_search.py sets it from --fasta).")
    tag, tag_source = decoy_tag(params)
    try:
        st = inspect(fas, tag)
    except OSError as e:
        sys.exit(f"{REFUSED} --fasta {fas} could not be read ({e}).")
    if st["state"] == "partial":
        miss = ", ".join(st["targets_without_decoy"]) or "-"
        orph = ", ".join(st["decoys_without_target"]) or "-"
        sys.exit(f"{REFUSED} --fasta {fas} holds {st['n_decoys']:,} decoys (entries starting "
                 f"{tag}, the workflow's decoy tag; {tag_source}) for {st['n_targets']:,} "
                 f"targets, but not one for each: targets without their decoy include {miss}; "
                 f"decoys without their target include {orph}. MSFragger would search it, and "
                 f"the FDR would be estimated on half a decoy set. Give the target FASTA alone "
                 f"(the search's job adds the decoys with FragPipe's Philosopher), or a complete "
                 f"target+decoy .fas built from it (philosopher database --custom <target> "
                 f"--prefix {tag}, without --contam).")
    sha, sha_source = sidecar_sha256(fas, sidecar)
    rec = {"path": fas, "target_fasta": fas, "decoys_built_in_job": False, "decoy_build": None,
           "sha256": sha, "sha256_source": sha_source,
           "n_entries": st["n_entries"], "n_targets": st["n_targets"], "n_decoys": st["n_decoys"],
           "decoy_tag": tag, "decoy_tag_source": tag_source, "workflow_db_path": given,
           "source": "--fasta, written as the workflow's database.db-path by run_search.py"}
    if st["state"] == "target_decoy":
        return rec
    exe = (philosopher or {}).get("path") if isinstance(philosopher, dict) else philosopher
    if not exe:
        sys.exit(f"{REFUSED} --fasta {fas} has no decoys (no entry starts with {tag}; "
                 f"{tag_source}), and MSFragger refuses a target-only FASTA (\"No decoys found in "
                 f"the FASTA file\"). The search's job adds them with the FragPipe install's "
                 f"Philosopher, but tools.json records none (fragpipe_tools.philosopher): re-run "
                 f"acquire_tools.sh, or give a complete target+decoy .fas.")
    built = planned_database(out, fas, today)
    # not measured yet: the job measures what it builds, and merge_job_record() then says
    # "measured in the job"; a job that never got that far leaves this, true as it stands (L1)
    rec.update(path=built, sha256=None,
               sha256_source=NOT_YET_MEASURED,
               n_entries=st["n_targets"] * 2, n_decoys=st["n_targets"], decoys_built_in_job=True,
               decoy_build={"tool": "Philosopher", "path": exe,
                            "version": (philosopher or {}).get("version")
                            if isinstance(philosopher, dict) else None,
                            "command": f"philosopher database --custom <target FASTA> --prefix "
                                       f"{tag}",
                            "contaminants_added": False,
                            "why": ("--fasta is targets only (fetch_fasta.py's database, and the "
                                    "one a replay rebuilds): MSFragger needs decoys, so the job "
                                    "adds them first, never with --contam (the skill's "
                                    "contaminants are already in the FASTA, tagged Cont_)")},
               source=("the target+decoy database the job builds from --fasta with FragPipe's "
                       "Philosopher, written as the workflow's database.db-path by "
                       "run_search.py"))
    return rec


def job_lines(rec, out, python_line, script=None):
    """The job's lines before FragPipe: build the decoys (a target-only --fasta), then measure
    the database it will search into <out>/fragpipe_job.json -- which fails the job on half a
    decoy set. `python_line`: run_search._job_python()'s line (defines $_py). The java text
    ($_fp_java_text, set by the job's Java check before these lines) goes into the record."""
    import shlex
    q = shlex.quote
    script = script or os.path.abspath(__file__)
    lines = []
    if rec.get("decoys_built_in_job"):
        db = rec["path"]
        dbdir = os.path.dirname(db)
        build = os.path.join(dbdir, ".build")
        phi = q(rec["decoy_build"]["path"])
        lines += [
            "# [fragpipe_db] --fasta is targets only: the decoys MSFragger needs, built with the "
            "install's Philosopher (never --contam)",
            f"if [ -L {q(dbdir)} ]; then echo \"[fragpipe_db] REFUSED: \"{q(dbdir)}\" is a link, "
            f"not a folder of its own -- nothing was built or searched\" >&2; exit 1; fi",
            f"rm -rf -- {q(build)} && mkdir -p {q(build)}",
            f"( cd {q(build)} && {phi} workspace --init --nocheck && {phi} database --custom "
            f"{q(rec['target_fasta'])} --prefix {q(rec['decoy_tag'])} && {phi} workspace --clean "
            f"--nocheck )",
            f"_fp_built=$(find {q(build)} -maxdepth 1 -name '*-decoys-*.fas' -type f -print -quit)",
            f'[ -n "$_fp_built" ] || {{ echo "[fragpipe_db] Philosopher wrote no .fas in" '
            f'{q(build)} >&2; exit 1; }}',
            # the same bytes as a database already there: keep it (and MSFragger's index of it)
            f'if [ -f {q(db)} ] && cmp -s "$_fp_built" {q(db)}; then rm -f -- "$_fp_built"; else '
            f'rm -f -- {q(db)}.*.pepindex; mv -f -- "$_fp_built" {q(db)}; fi',
            f"rm -rf -- {q(build)}",
            f'echo "[fragpipe_db] target+decoy database:" {q(db)}']
    lines += [python_line,
              f'"$_py" {q(script)} measure --db {q(rec["path"])} --tag {q(rec["decoy_tag"])} '
              f'--target {q(rec["target_fasta"])} --record {q(os.path.join(out, JOB_RECORD))} '
              f'--java-text "${{_fp_java_text:-}}"']
    return lines


def measure(db, tag, target=None, java_text=None):
    """What the job measured: the searched database's sha256 and decoy state (and its target
    FASTA's sha256 when it was built from one), the java FragPipe ran on."""
    st = inspect(db, tag)
    rec = {"database": {"path": os.path.abspath(db), "sha256": _sha256(db),
                        "n_entries": st["n_entries"], "n_targets": st["n_targets"],
                        "n_decoys": st["n_decoys"], "state": st["state"], "decoy_tag": tag,
                        "measured": "in the job, before FragPipe started"}}
    if target and os.path.abspath(target) != os.path.abspath(db):
        rec["database"]["target_fasta"] = {"path": os.path.abspath(target),
                                           "sha256": _sha256(target)}
    if java_text is not None:
        text = java_text.strip()
        rec["java"] = {"text": text, "major": java_major(text), "minimum": MIN_JAVA}
    return rec, st


def merge_job_record(out, rec):
    """The generation record `rec` (fragpipe_database) with what the job measured
    (<out>/fragpipe_job.json) merged in -> (merged record, java record or None, note or None).
    The measured sha256 wins; a sidecar sha256 that disagrees with it is said, never kept."""
    try:
        with open(os.path.join(out, JOB_RECORD)) as fh:
            job = json.load(fh)
    except (OSError, ValueError):
        return rec, None, None
    meas = job.get("database") if isinstance(job, dict) else None
    if not isinstance(meas, dict):
        return rec, job.get("java") if isinstance(job, dict) else None, None
    merged, note = dict(rec or {}), None
    if merged.get("path") and os.path.abspath(merged["path"]) != meas.get("path"):
        return rec, job.get("java"), (f"{JOB_RECORD} measured {meas.get('path')}, not the "
                                      f"database recorded ({merged.get('path')}): not merged")
    if merged.get("sha256") and merged["sha256"] != meas.get("sha256"):
        note = (f"the sidecar's sha256 ({merged['sha256'][:12]}...) is not the database's the "
                f"job measured ({str(meas.get('sha256'))[:12]}...): the file changed since its "
                f"sidecar was written; the measured one is recorded")
    merged.update({k: meas[k] for k in ("sha256", "n_entries", "n_targets", "n_decoys")
                   if k in meas})
    merged["sha256_source"] = "measured in the job, before FragPipe started"
    if meas.get("target_fasta"):
        merged["target_fasta_sha256"] = meas["target_fasta"].get("sha256")
    if note:
        merged["sha256_note"] = note
    return merged, job.get("java"), note


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("measure", help="the job's record of the database it searches")
    m.add_argument("--db", required=True)
    m.add_argument("--tag", default=DEFAULT_DECOY_TAG)
    m.add_argument("--target")
    m.add_argument("--record", required=True)
    m.add_argument("--java-text")
    a = ap.parse_args(argv)
    try:
        rec, st = measure(a.db, a.tag, a.target, a.java_text)
    except OSError as e:
        sys.exit(f"[fragpipe_db] the database {a.db} could not be read ({e}): nothing was "
                 f"searched")
    tmp = a.record + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(rec, fh, indent=2)
    os.replace(tmp, a.record)
    d = rec["database"]
    print(f"[fragpipe_db] {d['path']}: {d['n_entries']:,} entries, {d['n_decoys']:,} decoys "
          f"tagged {a.tag}, sha256 {d['sha256']}")
    if st["state"] != "target_decoy":
        sys.exit(f"[fragpipe_db] REFUSED: {d['path']} is not a complete target+decoy database "
                 f"({st['state']}: {st['n_decoys']:,} decoys tagged {a.tag} for "
                 f"{st['n_targets']:,} targets) -- nothing was searched")


if __name__ == "__main__":
    main()
