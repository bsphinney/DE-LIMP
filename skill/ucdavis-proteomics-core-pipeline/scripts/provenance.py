#!/usr/bin/env python3
"""
provenance.py  --  Assemble a COMPLETE, self-contained reproducibility bundle.

A result that can't be reproduced isn't a result. After a run finishes, the
orchestrator MUST call this to capture everything needed to reproduce the
analysis byte-for-byte: exact tool + package versions, the pinned registry
commit, every parameter, input/output checksums, and a runnable `reproduce.sh`.

It never throws on a missing piece — like DE-LIMP's safe_section(), it records
`[SKIPPED] <what> -- <why>` in MANIFEST.txt and keeps going, so the bundle is
honest about what it could and couldn't capture.

Usage (the orchestrator fills these from earlier steps):
  python3 provenance.py \
    --outdir ./reproducibility \
    --workflow-manifest ./wf/workflow.manifest.json \
    --setup-json   ~/.proteomics-pipeline/setup.json \
    --tools-json   ~/.proteomics-pipeline/tools/tools.json \
    --params       ./wf/diann.cfg \
    --conditions   ./conditions.csv \
    --fasta        ./search.fasta \
    --fasta-info-file ./search.fasta.meta.json \
    --raw          /data/*.d \
    --report       ./search_out/report.parquet \
    --search-record ./search_out/search_provenance.json \
    --de-dir       ./de_results \
    --engine diann --de-method dpc --contrasts "B-A,C-A" \
    --q-cutoff 0.01 --logfc 1.0 --adjp 0.05 \
    --organism-taxid 9606 --instrument "Orbitrap Astral" --acquisition DIA \
    --commands ./commands.log         # optional: a log of the exact commands run
    --qc-bracket <session>/logs/qc_bracket.json   # optional: the instrument's QC record
                                      # (qc_bracket.py) -- referenced by path + sha256 only
    --peak-check <session>/output/skyline/skyline_provenance.json
                                      # optional (default: the one beside --de-dir's folder):
                                      # the peak check's record (skyline_peaks.py), copied in

Outputs under --outdir:
  run_manifest.json          everything, machine-readable
  REPRODUCE.md               human-readable methods + how-to-rerun
  reproduce.sh               re-creates env, rebuilds the FASTA, re-runs the search with the
                             parameters it ran with (inputs/), and the DE when one ran
  MANIFEST.txt               [OK]/[SKIPPED] log of what was captured
  environment/               conda-explicit.txt, pip-freeze.txt, r-sessionInfo.txt, versions.txt
  inputs/                    copies of params, conditions.csv, the workflow manifest (and a
                             FragPipe search's workflow as it ran, fragpipe_input.workflow)
  checksums/                 sha256 of inputs, report, and DE outputs
"""
import sys, os, json, glob, shutil, hashlib, argparse, subprocess, platform, shlex

# one definition each, where it is written
from fetch_fasta import (KEEP_TARGET_CONTAMINANTS_RULE, MIN_UNIQUE_PEPTIDES, sidecar_state,
                         keratin_sample_recorded, enzyme_arg)
from skill_version import skill_version, plugin_meta, label as skill_label

MANIFEST_LINES = []
def ok(msg):      MANIFEST_LINES.append(f"[OK]      {msg}")
def skip(w, why): MANIFEST_LINES.append(f"[SKIPPED] {w} -- {why}")
def note(msg):    MANIFEST_LINES.append(f"[NOTE]    {msg}")

MAX_HASH_BYTES = 5 * 1024**3   # don't sha256 files bigger than 5 GB (record size instead)


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fingerprint(path):
    """sha256 for a normal file; for big files / directories (.d), a structural
    fingerprint (sorted name+size list) so re-runs can detect input drift."""
    p = path.rstrip("/")
    if os.path.isdir(p):
        entries = []
        for dp, _, fns in os.walk(p):
            for fn in sorted(fns):
                fp = os.path.join(dp, fn)
                try:
                    entries.append(f"{os.path.relpath(fp, p)}\t{os.path.getsize(fp)}")
                except OSError:
                    pass
        blob = "\n".join(sorted(entries)).encode()
        return {"path": p, "type": "dir", "n_files": len(entries),
                "size_bytes": sum(os.path.getsize(os.path.join(dp, fn))
                                  for dp, _, fns in os.walk(p) for fn in fns),
                "structure_sha256": hashlib.sha256(blob).hexdigest()}
    try:
        size = os.path.getsize(p)
    except OSError as e:
        return {"path": p, "type": "missing", "error": str(e)}
    if size > MAX_HASH_BYTES:
        return {"path": p, "type": "file", "size_bytes": size,
                "sha256": None, "note": "too large to hash; size recorded"}
    return {"path": p, "type": "file", "size_bytes": size, "sha256": sha256_file(p)}


def run_capture(cmd, timeout=120):
    try:
        r = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=timeout)
        return (r.stdout or "") + (r.stderr or "")
    except Exception as e:
        return f"(could not run `{cmd}`: {e})"


def expand(patterns):
    out = []
    for p in patterns or []:
        out.extend(sorted(glob.glob(p)) or [p])
    return out


def _copy_into(src, in_dir, said):
    """Copy src into in_dir and return the name it has there. A file of that name already there
    is never replaced: the same file or the same bytes are left as they are; a different file is
    kept, and src gets a distinct name (<stem>.de_design<ext>, then .2, .3 ...), which is said."""
    import filecmp
    name = os.path.basename(src)
    dest = os.path.join(in_dir, name)
    if os.path.exists(dest):
        if os.path.samefile(src, dest) or filecmp.cmp(src, dest, shallow=False):
            said(f"the DE's design ({name}, de_provenance.json metadata) is already in inputs/")
            return name
        stem, ext = os.path.splitext(name)
        n, cand = 1, f"{stem}.de_design{ext}"
        while os.path.exists(os.path.join(in_dir, cand)) and not filecmp.cmp(
                src, os.path.join(in_dir, cand), shallow=False):
            n += 1
            cand = f"{stem}.de_design.{n}{ext}"
        shutil.copy2(src, os.path.join(in_dir, cand))
        note(f"inputs/{name} is another file than the DE's design ({src}): both kept, the "
             f"design as inputs/{cand}, which reproduce.sh step 6 reads")
        return cand
    shutil.copy2(src, dest)
    said(f"copied the DE's design ({name}, de_provenance.json metadata)")
    return name


def _utc_now():
    """ISO-8601 UTC. Separate helper so the default is obvious and testable."""
    import datetime
    return datetime.datetime.now(datetime.timezone.utc).replace(
        microsecond=0).isoformat().replace("+00:00", "Z")


def load_json(path):
    try:
        return json.load(open(path))
    except Exception:
        return None


def read_fasta_info(text=None, path=None):
    """fetch_fasta.py's sidecar (<fasta>.meta.json) for the bundle -> (dict, None), or (None, the
    problem). From --fasta-info-file <path> or --fasta-info @<path>, or --fasta-info '<JSON>' (the
    older form). A file, because the JSON on a command line can exceed the system's argument limit:
    a database with 8,685 added sequences has a 2 MB sidecar, and `--fasta-info "$(cat ...)"` failed
    with "Argument list too long" (exit 126; staff report 2026-10-09)."""
    if text and text.startswith("@"):
        path, text = text[1:], None
    if path:
        try:
            with open(path, encoding="utf-8") as fh:
                text = fh.read()
        except OSError as e:
            return None, f"cannot read {path}: {e}"
    try:
        info = json.loads(text)
    except (TypeError, ValueError) as e:
        return None, f"not valid JSON ({e})"
    if not isinstance(info, dict):
        return None, f"a JSON {type(info).__name__}, not fetch_fasta.py's sidecar (a JSON object)"
    added = info.get("added_sequences")
    if added is not None and not (isinstance(added, list) and
                                  all(isinstance(f, dict) for f in added)):
        return None, ("added_sequences must be a list of objects, one per --add-fasta file "
                      "({file, sha256, n_entries, entries}), as fetch_fasta.py writes it; got "
                      f"{type(added).__name__}. Pass the sidecar fetch_fasta.py wrote, unedited")
    return info, None


def fasta_info_for_bundle(info):
    """The sidecar as the bundle records it: each --add-fasta file by path, sha256 and entry count.
    The per-entry records stay in the sidecar (thousands of them for an added proteome); the file's
    sha256 already pins every entry."""
    if not isinstance(info, dict) or not info.get("added_sequences"):
        return info
    return dict(info, added_sequences=[
        {k: f.get(k) for k in ("file", "sha256", "n_entries")} for f in info["added_sequences"]])


# PROT_0756 v2 (2026-09-28): identical inputs, software and flags on a zen4 and a zen2 HIVE node.
CPU_FAMILY_EVIDENCE = ("|ΔlogFC| ≤ 0.0022 and |Δt| ≤ 0.0064 across 12 contrasts, no significance "
                       "call changed; each family reproduced its own tables exactly")


def de_runtime(de_dir):
    """(de_provenance.json, its `runtime`) -- the R, libraries, conda env and container the DE ran
    in, as run_de.R recorded them in that process. runtime None: no such record (no --de-dir, or a
    run_de.R before skill 2.10)."""
    rec = load_json(os.path.join(de_dir, "de_provenance.json")) if de_dir else None
    rec = rec if isinstance(rec, dict) else {}
    rt = rec.get("runtime")
    return rec, (rt if isinstance(rt, dict) else None)


def _relabel_wording(text, flag, role):
    """`text` with each override_source wording of `flag` saying `role` where it said staff.ROLE.
    The wording's two halves around the value come from override_source itself (the value as
    shown differs: resolve_defaults' ppm_source shows 25, its record 25.0), so nothing here
    restates it."""
    import staff
    from estimate_params import OVERRIDE_REASON_MISSING, override_source
    mark = "\x00"
    head, _, old = override_source(flag, mark, None, staff.ROLE).partition(mark)
    new = override_source(flag, mark, None, role).partition(mark)[2]
    cut = len(OVERRIDE_REASON_MISSING)
    old, new = old[:-cut], new[:-cut]                    # ", set by <role>; reason: "
    i = text.find(head)
    while i >= 0:
        j = text.find(old, i + len(head))
        if j < 0:
            break
        text = text[:j] + new + text[j + len(old):]
        i = text.find(head, j + len(new))
    return text


def relabel_old_overrides(doc, overrides, staff_path, name):
    """A record copied into the bundle -- a params sidecar, a workflow manifest -- whose overrides
    say "set by Core staff" (staff.ROLE): before 2.14.2 every override did, whoever set it (S4).
    2.14.2 writes the role it gave beside the setter in the staff-only record `staff_path`
    (estimate_params.record_override_setter `role`); an entry that records one leaves `doc` as it
    is. Otherwise the role is the entry's setter's now (staff.role_for: Core staff only for a login
    on a usable Core staff list), or the user when there is no entry to check, and every label is
    rebuilt with it (set_by, override_source's wording). -> (the record to bundle, a MANIFEST NOTE
    naming no one -- or None when nothing changed). The session's own file is never changed."""
    import staff
    flags = [k for k, o in (overrides or {}).items()
             if isinstance(o, dict) and o.get("set_by") == staff.ROLE]
    if not flags:
        return doc, None
    try:
        entries = [e for e in staff.read(staff_path) if e.get("event") == "override"]
        why = None if entries else ("has no override entry" if os.path.exists(staff_path)
                                    else "is not here")
    except (OSError, ValueError) as e:
        entries, why = [], f"could not be read ({e})"
    entry = entries[-1] if entries else None
    if entry is not None and "role" in entry:
        return doc, None
    role = staff.role_for(entry.get("by")) if entry else staff.USER_ROLE
    if role == staff.ROLE:
        return doc, None
    for k in flags:
        overrides[k]["set_by"] = role

    def walk(x):
        if isinstance(x, dict):
            return {key: walk(v) for key, v in x.items()}
        if isinstance(x, list):
            return [walk(v) for v in x]
        if isinstance(x, str):
            for k in flags:
                x = _relabel_wording(x, k, role)
        return x
    doc = walk(doc)
    return doc, (
        f"{name}: its override{'s' if len(flags) > 1 else ''} {', '.join(flags)} said \"set by "
        f"{staff.ROLE}\", as every override did before 2.14.2 whoever set it; re-labelled "
        f"\"{role}\" here -- "
        + ("the setter its staff-only record names is not on a usable Core staff list "
           "(staff.role_for)" if entry else
           f"its staff-only record ({os.path.basename(staff_path)}) {why}, so who set it cannot "
           f"be checked")
        + ". The session's own file is unchanged.")


def copy_overrides_record(src, dst, overrides_of, staff_path):
    """Copy a record that may carry overrides (`overrides_of(doc)`: its `overrides` mapping) to
    `dst`, re-labelled by relabel_old_overrides when it is a pre-2.14.2 one. `staff_path`: its
    staff-only record when its overrides name none (their `staff_record`, beside `src`). -> the
    NOTE, or None. A file that is not JSON is copied as it is, said."""
    import staff
    doc = load_json(src)
    ov = overrides_of(doc) if isinstance(doc, dict) else None
    if not isinstance(ov, dict) or not ov:
        shutil.copy2(src, dst)
        return None if isinstance(doc, dict) else (
            f"{os.path.basename(dst)} is not readable JSON: copied as it is, its override labels "
            f"not checked")
    rec = next((o.get("staff_record") for o in ov.values()
                if isinstance(o, dict) and o.get("staff_record")), None)
    if rec:
        staff_path = os.path.join(os.path.dirname(os.path.abspath(src)), rec)
    doc, msg = relabel_old_overrides(doc, ov, staff_path, os.path.basename(dst))
    if msg is None:
        shutil.copy2(src, dst)
    else:
        with open(dst, "w", encoding="utf-8") as fh:
            json.dump(doc, fh, indent=2)
    return msg


def _contrast_list(v):
    """Contrasts as a list, trimmed: run_de.R's record (a list; a bare string for one contrast in
    records written before it was kept a list) or a --contrasts value (split at commas, as run_de.R
    splits it)."""
    items = v if isinstance(v, list) else str(v or "").split(",")
    return [str(x).strip() for x in items if str(x).strip()]


def replayed_contrasts(recorded, given):
    """(the contrasts reproduce.sh replays, as run_de.R --contrasts takes them, "" for none; a
    MANIFEST.txt NOTE, or None). The DE's own record (de_provenance.json `contrasts`) is what ran,
    and is replayed whatever --contrasts says: SKILL.md step 10 shows --contrasts as a placeholder,
    and since 2.14.2 the default contrasts are "<group>-<control>" (experiment_type
    default_contrasts) while run_de.R's own default is not, so a blank or stale value replayed
    the contrasts reversed -- every fold change's sign flipped, silently (C1 review of 2.14.2, F2).
    A --contrasts that differs from the record is named beside it. A record with none: --contrasts,
    said; neither: run_de.R's own default, said."""
    rec, giv = _contrast_list(recorded), _contrast_list(given)
    if rec:
        return ",".join(rec), (
            f"--contrasts \"{given}\" differs from the contrasts the DE ran ({', '.join(rec)}, its "
            f"de_provenance.json): reproduce.sh replays the DE's, run_manifest.json and "
            f"REPRODUCE.md name them" if giv and giv != rec else None)
    if giv:
        return ",".join(giv), (
            f"the DE record names no contrasts (no de_provenance.json `contrasts`): reproduce.sh "
            f"replays --contrasts \"{given}\" as given, not checked against what ran [confirm]")
    return "", ("neither the DE record nor --contrasts names the contrasts: reproduce.sh's DE step "
                "runs run_de.R's own default contrasts, which may not be the ones that ran "
                "[confirm]")


def conda_prefix_of_r_home(r_home):
    """A conda env's R lives at <prefix>/lib/R: that prefix when it is a conda env, else None."""
    if not r_home:
        return None
    cand = os.path.dirname(os.path.dirname(os.path.normpath(r_home)))
    return cand if os.path.isdir(os.path.join(cand, "conda-meta")) else None


def env_of_rscript(path):
    """The conda env an Rscript belongs to (<prefix>/bin/Rscript or <prefix>/lib/R/bin/Rscript),
    or None."""
    if not path:
        return None
    d = os.path.dirname(os.path.realpath(path))                    # .../bin
    for cand in (os.path.dirname(d), os.path.dirname(os.path.dirname(os.path.dirname(d)))):
        if os.path.isdir(os.path.join(cand, "conda-meta")):
            return cand
    return None


def runtime_section(rec, rt, setup_env):
    """REPRODUCE.md's 'environment the DE ran in' section, from run_de.R's record."""
    head = "## The environment the DE ran in\n\n"
    if not rt:
        return (head + "NOT RECORDED -- this DE record has no `runtime` (a run_de.R before skill "
                "2.10), so `environment/` describes setup.json's environment, which may not be "
                "the one the DE ran in. The DE's own `sessionInfo.txt` (beside its tables) is the "
                "record of its R packages.\n")
    pk = rec.get("packages") or {}
    pkgs = ", ".join(f"{k} {v}" for k, v in pk.items() if v) or "not recorded"
    where = (f"inside the container `{rt['container']}`" if rt.get("container") else
             "in a Docker container (image not recorded)" if rt.get("docker") else
             "outside any container")
    env = rt.get("conda_prefix") or conda_prefix_of_r_home(rt.get("r_home"))
    lines = [f"- {rt.get('r_version') or 'R version not recorded'} (`{rt.get('r_home')}`), {where}",
             f"- packages: {pkgs}",
             f"- conda env: `{env}`" if env else "- conda env: none recorded",
             f"- library paths: {', '.join(f'`{x}`' for x in rt.get('lib_paths') or [])}"]
    if setup_env:
        lines.append(f"- setup.json names a different environment ({setup_env}); the bundle "
                     "describes the one above, which is the one the DE ran in")
    return (head + "As run_de.R recorded it in the process that ran the DE "
            "(`de_provenance.json` `runtime`); `environment/` was captured from it:\n\n"
            + "\n".join(lines) + "\n")


def compute_section(compute):
    """REPRODUCE.md's 'same numbers need the same kind of CPU' section, from run_de.R's record."""
    if not compute:
        return ("## Exact numbers need the same kind of CPU\n\n"
                "This DE record does not name the machine it ran on (a run_de.R from before the "
                "`compute` record). DPC-Quant's numbers are exact only on the same CPU family "
                "(`references/reproducibility.md`).\n")
    fam = compute.get("cpu_family")
    where = (f"`{compute.get('cpu_model')}`" +
             (f" on HIVE node `{compute.get('slurm_node')}` (CPU family `{fam}`)" if fam else ""))
    how = (f"On HIVE, submit the re-run to the same family: `sbatch --constraint={fam} ...` "
           "(`sinfo -N -o '%N %f'` lists each node's features)." if fam else
           "On HIVE, find the family of the node that ran it with `sinfo -N -o '%N %f'` (its "
           "`zen2` / `zen3` / `zen4` ... feature) and submit the re-run with "
           "`sbatch --constraint=<family>`; elsewhere, re-run on the same kind of CPU.")
    return f"""## Exact numbers need the same kind of CPU

The DE ran on {where}, BLAS `{compute.get('blas')}`, LAPACK `{compute.get('lapack')}`
(de_provenance.json `compute`). DPC-Quant fits each protein with an optimiser that stops at a
tolerance, and the BLAS library's CPU-specific kernels round the last bit differently, so the
same inputs and software give slightly different numbers on another CPU family. Measured on
PROT_0756 v2 (a zen4 and a zen2 HIVE node): {CPU_FAMILY_EVIDENCE}. {how} Do not set
`OPENBLAS_CORETYPE` to force one kind of kernel: it crashed limpa's dpcQuant on an EPYC 9734.
"""


def de_ran_of(de_rec, de_method, de_dir):
    """Whether a differential expression ran in this analysis -> (bool, why): its record
    (de_provenance.json in --de-dir), a --de-method, or DE tables in --de-dir. With none of them
    the bundle replays the search alone -- no step 6 with a literal `--method <method>`, no
    "the environment the DE ran in" (2.14.2 C2, E2E F2)."""
    if de_rec:
        return True, "de_provenance.json"
    if de_method:
        return True, f"--de-method {de_method}"
    if de_dir and glob.glob(os.path.join(de_dir, "DE_*.csv")):
        return True, f"DE tables in {de_dir}"
    return False, "no de_provenance.json, no --de-method and no DE tables"


def search_record_of(given, report):
    """(path, record) of the search's own record: --search-record, else the search_provenance.json
    beside --report; (path or None, None) when there is none or it is not a JSON object."""
    path = given or (os.path.join(os.path.dirname(os.path.abspath(report)),
                                  "search_provenance.json") if report else None)
    rec = load_json(path) if path and os.path.isfile(path) else None
    return path, (rec if isinstance(rec, dict) else None)


def _record_key(sp, key):
    """A search record's entry `key`, top level first, else in its `result` (run_search.py lifts
    some result entries to the top: result_extras)."""
    sp = sp if isinstance(sp, dict) else {}
    v = sp.get(key)
    if v is None and isinstance(sp.get("result"), dict):
        v = sp["result"].get(key)
    return v


# The line of a FragPipe workflow that names its database. run_search.py sets it from --fasta
# (fragpipe_workflow_for_run) and refuses a workflow that names another, so the replay's copy of
# the workflow that ran leaves it out: the replay's --fasta is the database rebuilt in step 4.
FRAGPIPE_DB_LINE = "database.db-path="


def replay_params(params, sp, engine, in_dir):
    """Step 5's --params: the parameters the search RAN with, copied into inputs/ -- never
    re-derived from the data type, which loses the user's overrides, an identifications-only
    cfg and a non-trypsin enzyme, and for Sage or DIA-NN (no preset in step 2) left `--params
    ./wf/` (2.14.2 C2, E2E F2). -> (the --params word, the note on step 5's heading, what
    MANIFEST.txt says):
      * FragPipe: the workflow it ran (search_provenance.json fragpipe_workflow.run, i.e.
        fragpipe_input.workflow) without its database line (FRAGPIPE_DB_LINE);
      * else the file the search was given (--params; finalize passes the record's
        params_file), already in inputs/;
      * none recorded: a preset engine (resolve_defaults.PRESET_ENGINES) searches with the preset
        step 2 derives, said to be the skill's default; any other engine has nothing to replay,
        and the word says so."""
    import run_search
    from resolve_defaults import PRESET_ENGINES
    fp = _record_key(sp, "fragpipe_workflow")
    run_wf = fp.get("run") if isinstance(fp, dict) else None
    if isinstance(run_wf, str) and os.path.isfile(run_wf):
        name = run_search.FRAGPIPE_RUN_WORKFLOW
        if params and os.path.basename(params) == name:
            name = name.replace(".workflow", ".as_run.workflow")   # never over the given one
        with open(run_wf, encoding="utf-8", errors="replace", newline="") as fh:
            lines = fh.read().splitlines(keepends=True)
        with open(os.path.join(in_dir, name), "w", encoding="utf-8", newline="") as fh:
            fh.write("".join(ln for ln in lines if not ln.startswith(FRAGPIPE_DB_LINE)))
        return (f"./inputs/{name}",
                f"\n#    The workflow FragPipe ran ({os.path.basename(run_wf)}, in inputs/ as "
                f"{name}) is searched with, its\n#    database line left out: run_search.py sets "
                f"it from --fasta, the database step 4 rebuilds.",
                ("OK", f"copied the workflow FragPipe ran ({os.path.basename(run_wf)}) as "
                       f"inputs/{name}, without its {FRAGPIPE_DB_LINE.rstrip('=')}: step 5 "
                       f"searches with it"))
    if params and os.path.isfile(params):
        name = os.path.basename(params)
        return (f"./inputs/{name}",
                f"\n#    The parameters the search was given (inputs/{name}) are searched with, as "
                f"they ran -- its\n#    overrides, enzyme and quantification setting; step 2 "
                f"only rebuilds the workflow manifest.",
                ("OK", f"step 5 searches with the parameters the search was given "
                       f"(inputs/{name}), not re-derived ones"))
    if engine in PRESET_ENGINES:
        return ('./wf/$(basename "$(ls wf | grep -vi manifest | head -n1)")',
                "\n#    No parameters file was recorded: this searches with the preset step 2 "
                "derives for the\n#    data type -- the skill's default, which may not be what "
                "ran [confirm].",
                ("SKIPPED", "no --params and no record of a FragPipe workflow: step 5 searches "
                            "with step 2's preset, the skill's default, which may not be what "
                            "ran"))
    return ("'<the search parameters file: not recorded>'",
            "\n#    No parameters file was recorded, so this search cannot be replayed as it ran: "
            "give --params\n#    the file the search was given.",
            ("SKIPPED", f"no --params, and step 2 writes none for {engine or 'this engine'}: "
                        f"step 5 names none (fill it in)"))


def identified_only(sp_path, sp):
    """True when the search this bundle replays was for identifications only: NOT
    sage_lfq_check.search_quantified -- THE predicate, which the DE gate, the README, the email,
    the report and the brief share -- asked of the record's folder with the record read here
    (--search-record). False with no record (nothing says it did not quantify). One definition:
    step 5's notes and REPRODUCE.md read this, never a key of the record themselves (C2, rule 3)."""
    if not isinstance(sp, dict) or not sp_path:
        return False
    import sage_lfq_check
    return not sage_lfq_check.search_quantified(os.path.dirname(os.path.abspath(sp_path)),
                                                record=sp)


def analysis_section(sp, engine, de_ran, ids_only=False):
    """REPRODUCE.md's "What this analysis was", for an analysis with no DE: the search, how its
    inputs went in, and -- for one that did not quantify (`ids_only`, identified_only) -- where
    its identifications are and how they were kept, each from the search record's own words. ""
    for an analysis with a DE (its sections say what it was)."""
    if de_ran:
        return ""
    from make_methods import ENGINE_LABEL
    sp = sp if isinstance(sp, dict) else {}
    eng = sp.get("engine") or engine
    label = ENGINE_LABEL.get(str(eng or "").lower(), eng or "the search engine (not recorded)")
    ver = sp.get("version")
    lines = [f"- Search: {label} {ver or '(version not recorded)'}"
             + (" (its record: search_provenance.json)" if sp else
                " (no search record was given: --search-record)")]
    route = _record_key(sp, "sage_input")
    if isinstance(route, dict) and route.get("methods_sentence"):
        lines.append(f"- How the raw files went in: {route['methods_sentence']}")
    ids = next((x for x in (sp.get("sage_identifications"), sp.get("fragpipe_identifications"))
                if isinstance(x, dict)), None)
    if ids_only:
        lines.append("- Identifications only: the search did not quantify, so there is no "
                     "report.parquet and no differential expression. Its results are "
                     "`identifications.tsv` and the record of how its rows were kept, beside the "
                     "search record.")
        if ids:
            lines += [f"- {k}: {ids[key]}" for k, key in (("Rows", "row_unit"),
                                                           ("Kept", "kept_rule"))
                      if ids.get(key)]
    else:
        lines.append("- No differential expression ran: the search's report is the result.")
    return "## What this analysis was\n\n" + "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="./reproducibility")
    ap.add_argument("--workflow-manifest")
    ap.add_argument("--setup-json", default=os.path.expanduser("~/.proteomics-pipeline/setup.json"),
                    help="defaults to ~/.proteomics-pipeline/setup.json (where setup.sh writes it)")
    ap.add_argument("--tools-json")
    ap.add_argument("--params")
    ap.add_argument("--conditions")
    ap.add_argument("--fasta")
    ap.add_argument("--fasta-info-file", metavar="PATH",
                    help="fetch_fasta.py's <fasta>.meta.json -- the database as it was built")
    ap.add_argument("--fasta-info", help="the same, as @<path> or (older form) its JSON text; "
                                         "use the file form: a big sidecar exceeds the "
                                         "command-line limit")
    ap.add_argument("--raw", nargs="*")
    ap.add_argument("--report")
    ap.add_argument("--search-record", metavar="PATH",
                    help="the search's own record (<search out>/search_provenance.json): step 5 "
                         "replays the parameters it ran with, and REPRODUCE.md says what it was. "
                         "Default: the one beside --report")
    ap.add_argument("--search-output", action="append", default=[], metavar="PATH",
                    help="a file the search wrote, fingerprinted in checksums.json "
                         "`search_outputs` (repeat it): identifications.tsv and its record for "
                         "a search for identifications only, the engine's own tables "
                         "(results.sage.parquet). finalize passes session.search_result_files")
    ap.add_argument("--de-dir")
    ap.add_argument("--engine")
    ap.add_argument("--de-method")
    ap.add_argument("--contrasts", default="")
    ap.add_argument("--q-cutoff"); ap.add_argument("--logfc"); ap.add_argument("--adjp")
    ap.add_argument("--organism-taxid"); ap.add_argument("--instrument", default="")
    ap.add_argument("--acquisition")
    ap.add_argument("--commands", help="optional log file of the exact commands run")
    ap.add_argument("--qc-bracket", help="qc_bracket.py's record (logs/qc_bracket.json): "
                                         "referenced in run_manifest.json by path and sha256 "
                                         "only -- the bundle is delivered, the verdict is "
                                         "staff-only")
    ap.add_argument("--peak-check", help="skyline_peaks.py's record (output/skyline/"
                                         "skyline_provenance.json); default: the one beside the "
                                         "--de-dir folder, when there is one")
    ap.add_argument("--timestamp", default="", help="ISO timestamp (script can't call clock)")
    a = ap.parse_args()
    # sha256 of every input raw file (up to 5 GB each) and the report: I/O for a compute node.
    # On a cluster login node this runs itself again there (on_compute.py; HPC@UCD's ucd-hpc rules).
    try:
        import on_compute
        # restartable: a run stopped half way is simply run again (preemptible_ok, pick_queue)
        on_compute.ensure_compute_node("provenance.py hashes every input file", minutes=60, mem="4G",
                                       preemptible_ok=True, script=__file__)
    except ImportError:
        pass                              # a partial scripts/ without on_compute.py: as before
    # The database record is checked before anything is written: a bundle whose database replay
    # silently fell back to the workflow defaults would reproduce a different database.
    fasta_info = None
    # `is not None`: an empty value (`--fasta-info-file "$EMPTY_VAR"`) is a mistake to say, never
    # "not given" -- that quietly replayed the workflow's default database (error-handling review)
    if a.fasta_info_file is not None and a.fasta_info is not None:
        ap.error("give the FASTA sidecar once: --fasta-info-file <path> (or --fasta-info @<path>)")
    if a.fasta_info_file is not None or a.fasta_info is not None:
        flag = "--fasta-info-file" if a.fasta_info_file is not None else "--fasta-info"
        if not (a.fasta_info_file if a.fasta_info_file is not None else a.fasta_info).strip():
            ap.error(f"{flag} is empty: give fetch_fasta.py's <fasta>.meta.json (or leave the flag "
                     f"out, and the bundle says the database was not recorded)")
        fasta_info, problem = read_fasta_info(a.fasta_info, a.fasta_info_file)
        if problem:
            ap.error(f"{flag}: {problem}")
        fasta_info = fasta_info_for_bundle(fasta_info)

    out = os.path.abspath(a.outdir)
    for sub in ("environment", "inputs", "checksums"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    setup = load_json(a.setup_json) if a.setup_json else None
    tools = load_json(a.tools_json) if a.tools_json else None
    wfman = load_json(a.workflow_manifest) if a.workflow_manifest else None

    # ---- environment capture -------------------------------------------------
    # The environment the DE RAN in, as run_de.R recorded it (de_provenance.json `runtime`):
    # its Rscript and conda env. setup.json says only which env setup.sh built, and a DE run in
    # another (a separate R 4.6 env for limpa 1.4, staff report 2026-09-28) made the bundle name
    # the wrong R and limpa. setup.json / PATH only when the DE recorded no runtime -- said so.
    env_dir = os.path.join(out, "environment")
    de_rec, rt = de_runtime(a.de_dir)
    # Whether a DE ran at all (de_ran): a session with no DE gets no DE step and no DE text
    de_ran, de_ran_why = de_ran_of(de_rec, a.de_method, a.de_dir)
    # The search's own record (search_provenance.json): what step 5 replays, and what REPRODUCE.md
    # says the analysis was. --search-record (finalize passes it), else the one beside --report.
    sp_path, sp = search_record_of(a.search_record, a.report)
    if a.search_record and sp is None:
        skip("search record", f"not readable JSON or not found: {a.search_record}")
    # ...and whether it was for identifications only: the one predicate (identified_only)
    ids_only = identified_only(sp_path, sp)
    py = (setup or {}).get("python") or shutil.which("python3") or sys.executable
    conda = (setup or {}).get("conda") or shutil.which("micromamba") or shutil.which("mamba") or shutil.which("conda")
    setup_rscript, setup_prefix = (setup or {}).get("rscript"), (setup or {}).get("env_prefix")
    setup_differs = None
    if rt:
        rscript = rt.get("rscript")
        env_prefix = rt.get("conda_prefix") or conda_prefix_of_r_home(rt.get("r_home"))
        env_source = ("de_provenance.json runtime: recorded by run_de.R in the process that ran "
                      "the DE")
        # A value the DE's record does not carry is NOT RECORDED, never a difference: a SILAC
        # record without `rscript` (run_silac.R before 2.13.1) had the user told the DE ran in
        # another environment than setup.json's when it had not (2026-10).
        pairs = (("rscript", setup_rscript, rscript), ("env_prefix", setup_prefix, env_prefix))
        unrecorded = [k for k, sv, dv in pairs if sv and not dv]
        ok(f"environment: the one the DE ran in ({rt.get('r_version')}, {rt.get('r_home')})"
           + (f"; the DE's record has no {' / '.join(unrecorded)}, so setup.json's is not "
              f"compared" if unrecorded else ""))

        def _same(k, x, y):
            if os.path.realpath(x) == os.path.realpath(y):
                return True
            # one conda env's two Rscripts: setup.json names <env>/bin/Rscript, R.home("bin")
            # is <env>/lib/R/bin
            return k == "rscript" and env_of_rscript(x) is not None and \
                env_of_rscript(x) == env_of_rscript(y)
        diff = {k: (sv, dv) for k, sv, dv in pairs if sv and dv and not _same(k, sv, dv)}
        if diff:
            setup_differs = {k: {"setup_json": sv, "de": dv} for k, (sv, dv) in diff.items()}
            note("setup.json names a different environment than the one the DE ran in ("
                 + "; ".join(f"{k}: setup.json {sv}, DE {dv}" for k, (sv, dv) in diff.items())
                 + "): the bundle describes the DE's")
    elif not de_ran:
        rscript = setup_rscript or shutil.which("Rscript")
        env_prefix = setup_prefix
        env_source = "setup.json / PATH -- the skill's environment (no DE ran)"
        note(f"no differential expression ran ({de_ran_why}): reproduce.sh has no DE step, and "
             f"REPRODUCE.md describes the search alone")
    else:
        rscript = setup_rscript or shutil.which("Rscript")
        env_prefix = setup_prefix
        env_source = ("setup.json / PATH -- NOT necessarily the environment the DE ran in: "
                      "de_provenance.json records no runtime (run_de.R before skill 2.10)")
        skip("DE runtime", "de_provenance.json records none (run_de.R before skill 2.10, or no "
             "--de-dir): environment/ describes setup.json's environment, which may not be the "
             "one the DE ran in")
    if not rt and (not env_prefix or not os.path.isdir(env_prefix)) and py:
        # infer the env prefix from the interpreter location (…/<env>/bin/python)
        cand = os.path.dirname(os.path.dirname(os.path.realpath(py)))
        if os.path.isdir(os.path.join(cand, "conda-meta")):
            env_prefix = cand

    # conda explicit lock (fully pinned, URL+hash per package) — the gold standard
    if conda and env_prefix and os.path.isdir(env_prefix):
        txt = run_capture(f"{conda} list -p {env_prefix} --explicit --md5")
        if txt and "://" in txt:
            open(os.path.join(env_dir, "conda-explicit.txt"), "w").write(txt)
            ok(f"conda explicit lock of {env_prefix}" + (" (the env the DE ran in)" if rt else ""))
        else:
            skip("conda-explicit.txt", "conda list returned no URLs")
    else:
        skip("conda-explicit.txt",
             f"no conda / mamba / micromamba here to list {env_prefix}"
             if env_prefix and os.path.isdir(env_prefix) else
             f"the DE's R ({rt.get('r_home')}) is in no conda env it recorded" if rt else
             "no conda env found (setup.json + PATH)")

    # pip freeze (the env's python)
    if py:
        txt = run_capture(f"{py} -m pip freeze")
        open(os.path.join(env_dir, "pip-freeze.txt"), "w").write(txt); ok("pip freeze")

    # R sessionInfo with every package version (the DE step's exact stack): the one run_de.R
    # wrote as the DE ran, when it is there -- it IS that stack, whatever R this runs under
    si_de = os.path.join(a.de_dir or "", "sessionInfo.txt")
    if a.de_dir and os.path.isfile(si_de):
        shutil.copy2(si_de, os.path.join(env_dir, "r-sessionInfo.txt"))
        ok("R sessionInfo + package versions, as the DE recorded them (sessionInfo.txt)")
    elif rscript and (os.path.exists(rscript) or not rt):
        txt = run_capture(f"{rscript} -e 'sink(stdout()); "
                          f"cat(R.version.string,\"\\n\"); "
                          f"for(p in c(\"limpa\",\"limma\",\"arrow\",\"dplyr\",\"tidyr\")) "
                          f"try(cat(p, as.character(packageVersion(p)), \"\\n\")); "
                          f"cat(\"\\n\"); print(sessionInfo())'")
        open(os.path.join(env_dir, "r-sessionInfo.txt"), "w").write(txt)
        ok(f"R sessionInfo + package versions, from {'the DE' if rt else 'setup.json / PATH'}'s "
           f"Rscript ({rscript})" + ("" if rt or not de_ran else
                                     " -- may not be the R the DE ran in"))
    else:
        skip("r-sessionInfo.txt", "no sessionInfo.txt in --de-dir and no Rscript found ("
             + ("the DE's recorded Rscript is not here" if rt else "setup.json + PATH") + ")")

    # engine versions
    versions = {"os": platform.platform(), "python": platform.python_version()}
    if rt:
        # the DE's R stack, as run_de.R recorded it -- not whatever R is on PATH here
        versions["de_r"] = {"r_version": rt.get("r_version"), "packages": de_rec.get("packages"),
                            "container": rt.get("container"), "conda_prefix": env_prefix,
                            "source": "de_provenance.json runtime"}
    if tools:
        versions["tools_versions"] = tools.get("versions")
        for eng in ("diann", "sage"):
            cmd = tools.get(eng)
            if cmd:
                versions[f"{eng}_cmd"] = cmd
    if a.engine == "sage" and (setup or {}).get("sage"):
        # NOT "sage_version": `sage --version` is the binary's own claim about itself, and it
        # can be wrong about the release it came from -- the v0.14.7 release binary prints
        # "sage 0.14.6" (measured on HIVE 2026-09-16). Under the plain name, this file and
        # search_provenance.json contradicted each other for the same run, with nothing to say
        # which was which. Kept, under a name that says what it is: it is still the only thing
        # that reports the binary actually on disk, and a disagreement with the release is
        # itself worth seeing.
        versions["sage_self_reported_version"] = {
            "value": run_capture(f"{setup['sage']} --version").strip(),
            "note": "`sage --version` as the binary prints it. This is NOT necessarily the "
                    "release it came from: the v0.14.7 release binary prints 0.14.6. The "
                    "release of record is `tools_versions.sage` above (read from the release "
                    "tarball) and search_provenance.json `engine_version`.",
        }
    open(os.path.join(env_dir, "versions.txt"), "w").write(json.dumps(versions, indent=2)); ok("tool versions")

    # ---- which skill produced this + how it was installed --------------------
    skill_meta = plugin_meta()
    skill_info = {
        "name": skill_meta.get("name", "ucdavis-proteomics-core-pipeline"),
        "version": skill_version(),                 # skill_version.py: the one reader
        "title": "UC Davis Proteomics Core pipeline",
        "repository": skill_meta.get("repository", "https://github.com/bsphinney/DE-LIMP"),
        "marketplace": "ucdavis-proteomics-core",
        "install": [
            # by its HTTPS URL: for the owner/repo shorthand Claude Code tries SSH first, and a
            # laptop with no GitHub SSH key then cannot refresh it (references/install.md)
            "claude plugin marketplace add https://github.com/bsphinney/DE-LIMP.git",
            "claude plugin install ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core",
        ],
    }
    open(os.path.join(env_dir, "skill.txt"), "w", encoding="utf-8").write(
        "Produced by the {title} Claude skill.\n\n"
        "Skill:       {name} {shown}\n"
        "Repository:  {repository}\n"
        "Marketplace: {marketplace}\n\n"
        "Installed with:\n  {i0}\n  {i1}\n".format(
            i0=skill_info["install"][0], i1=skill_info["install"][1],
            shown=skill_label(skill_info["version"]), **skill_info)); ok("skill identity + install")

    # ---- copy inputs ---------------------------------------------------------
    import staff                  # who decided or set a value is staff-only: the bundle says roles
    in_dir = os.path.join(out, "inputs")
    for label, src in (("params", a.params), ("conditions.csv", a.conditions),
                       ("workflow.manifest.json", a.workflow_manifest),
                       ("commands.log", a.commands)):
        if src and os.path.exists(src) and label == "commands.log":
            # every command verbatim -- but who decided is staff-only: names become the role
            with open(src, encoding="utf-8", errors="replace") as fh:
                logged = fh.read()
            with open(os.path.join(in_dir, os.path.basename(src)), "w", encoding="utf-8") as fh:
                fh.write(staff.redact(logged))
            ok(f"copied {label} (names after --by replaced by \"<{staff.ROLE}>\", after --override-by by \"{staff.REDACTED_SETTER}\")")
        elif src and os.path.exists(src):
            if label == "workflow.manifest.json":
                # its --ms1-ppm/--ms2-ppm overrides (resolve_defaults.py), labelled by role:
                # a pre-2.14.2 manifest's are checked again (relabel_old_overrides)
                relabelled = copy_overrides_record(
                    src, os.path.join(in_dir, os.path.basename(src)),
                    lambda d: (d.get("search") or {}).get("overrides"),
                    src + staff.STAFF_SUFFIX)
                ok(f"copied {label}")
                if relabelled:
                    note(relabelled)
            else:
                shutil.copy2(src, os.path.join(in_dir, os.path.basename(src))); ok(f"copied {label}")
            # params estimated by estimate_params.py carry a sibling rationale — capture it
            if label == "params" and os.path.exists(src + ".rationale.json"):
                relabelled = copy_overrides_record(
                    src + ".rationale.json",
                    os.path.join(in_dir, os.path.basename(src) + ".rationale.json"),
                    lambda d: d.get("overrides"), src + staff.STAFF_SUFFIX)
                ok("copied params rationale (per-setting provenance)")
                if relabelled:
                    note(relabelled)
            # the user's sample-identity answers behind the design (collect_conditions.py
            # --confirm-multi), when there were any
            if label == "conditions.csv" and os.path.exists(src + ".decisions.json"):
                shutil.copy2(src + ".decisions.json",
                             os.path.join(in_dir, os.path.basename(src) + ".decisions.json"))
                ok("copied conditions decisions (labels the user said name several runs)")
        elif src:
            skip(label, f"not found: {src}")
    # The design the DE pipeline names in its own record (de_provenance.json `metadata`) when it
    # is not the conditions.csv given -- run_silac.R's channel design, for one: reproduce.sh's
    # step 6 reads it from inputs/<de_design_name>. The same file reached by another path (a
    # symlink, a relative path) is not copied again; a DIFFERENT file whose name is already taken
    # in inputs/ never replaces it -- it is kept beside it under a distinct name.
    de_design = de_rec.get("metadata")
    de_design_name = None
    if isinstance(de_design, str) and os.path.isfile(de_design):
        if a.conditions and os.path.isfile(a.conditions) and os.path.samefile(de_design, a.conditions):
            de_design_name = os.path.basename(a.conditions)
        else:
            de_design_name = _copy_into(de_design, in_dir, ok)

    # ---- checksums -----------------------------------------------------------
    checks = {}
    raw_files = expand(a.raw)
    checks["raw_inputs"] = [fingerprint(f) for f in raw_files]
    ok(f"fingerprinted {len(raw_files)} raw input(s)")
    if a.fasta and os.path.exists(a.fasta):
        checks["fasta"] = fingerprint(a.fasta); checks["fasta_info"] = fasta_info; ok("fingerprinted FASTA")
    elif a.fasta:
        skip("fasta", f"not found: {a.fasta}")
    if a.report and os.path.exists(a.report):
        checks["report"] = fingerprint(a.report); ok("fingerprinted search report")
    elif a.report:
        skip("report", f"not found: {a.report}")
    # What else the search wrote -- for a search for identifications only, its whole result
    # (identifications.tsv and its record) -- and the search's own record: their sha256, so a
    # replay can be checked against them (2.14.2 C2, E2E F3)
    outs = [x for x in dict.fromkeys(a.search_output) if x]
    found = [x for x in outs if os.path.isfile(x)]
    for x in outs:
        if x not in found:
            skip("search output", f"not found: {x}")
    if found:
        checks["search_outputs"] = [fingerprint(x) for x in found]
        ok(f"fingerprinted {len(found)} search output(s) "
           f"({', '.join(os.path.basename(x) for x in found)})")
    if sp is not None:
        checks["search_record"] = fingerprint(sp_path)
        ok(f"fingerprinted the search record ({os.path.basename(sp_path)})")
    de_outputs = []
    if a.de_dir and os.path.isdir(a.de_dir):
        for f in sorted(glob.glob(os.path.join(a.de_dir, "*"))):
            if os.path.isfile(f):
                de_outputs.append(fingerprint(f))
        checks["de_outputs"] = de_outputs; ok(f"fingerprinted {len(de_outputs)} DE output(s)")
    elif a.de_dir:
        skip("de_outputs", f"not a dir: {a.de_dir}")
    open(os.path.join(out, "checksums", "checksums.json"), "w").write(json.dumps(checks, indent=2))

    # ---- the instrument's QC around the project (qc_bracket.py), when it was checked --------
    # STAFF-ONLY: this bundle is delivered to the client, so it carries where the record is and
    # its checksum -- never the verdict, its summary or the record itself (logs/ in the session
    # and the run registry hold those).
    instrument_qc = None
    if a.qc_bracket and os.path.exists(a.qc_bracket):
        q = load_json(a.qc_bracket)
        if isinstance(q, dict) and str(q.get("schema", "")).startswith("qc_bracket/"):
            instrument_qc = {"record": os.path.abspath(a.qc_bracket),
                             "sha256": sha256_file(a.qc_bracket),
                             "schema_version": q.get("schema_version"),
                             "checked_at": q.get("checked_at"),
                             "note": "the instrument's QC around this project was checked; the "
                                     "record is Core-internal and not in this bundle"}
            ok("instrument QC record referenced (Core-internal; not copied)")
        else:
            skip("instrument QC (qc_bracket.json)", "not a qc_bracket.py record, or unreadable")
    elif a.qc_bracket:
        skip("instrument QC (qc_bracket.json)", f"not found: {a.qc_bracket}")

    # ---- the peak check (skyline_peaks.py, step 8f), when there is a record ---------------------
    # The tool, its container and checksum, the tolerance and its source, the proteins' rule and
    # N, and the outcome -- the record itself is copied in (it names no person).
    peak_check = None
    pk = a.peak_check or (os.path.join(os.path.dirname(os.path.abspath(a.de_dir)), "skyline",
                                       "skyline_provenance.json") if a.de_dir else None)
    if pk and os.path.exists(pk):
        rec = load_json(pk)
        if isinstance(rec, dict) and str(rec.get("schema", "")).startswith("skyline_peak_check/"):
            tool = rec.get("tool") or {}
            peak_check = {"record": os.path.abspath(pk), "sha256": sha256_file(pk),
                          "status": rec.get("status"), "reason": rec.get("reason"),
                          "skyline_version": tool.get("version"),
                          "container": tool.get("container"),
                          "container_sha256": tool.get("container_sha256"),
                          "container_verified": tool.get("verified"),
                          "docker_digest": tool.get("docker_digest"),
                          "mass_accuracy": rec.get("mass_accuracy"),
                          "proteins_n": (rec.get("top") or {}).get("value"),
                          "proteins_n_source": (rec.get("top") or {}).get("source"),
                          "copied_to": "inputs/skyline_provenance.json"}
            shutil.copy2(pk, os.path.join(out, "inputs", "skyline_provenance.json"))
            ok(f"peak check record copied (status {rec.get('status')})")
        else:
            skip("peak check (skyline_provenance.json)", "not a skyline_peaks.py record, or unreadable")
    elif a.peak_check:
        skip("peak check (skyline_provenance.json)", f"not found: {a.peak_check}")

    # ---- the contrasts the DE ran: what reproduce.sh, run_manifest.json and REPRODUCE.md name --
    de_contrasts, contrasts_note = replayed_contrasts(de_rec.get("contrasts"), a.contrasts)
    if contrasts_note and de_ran:                # no DE ran: there are no contrasts to replay
        note(contrasts_note)

    # ---- the master run manifest --------------------------------------------
    reg = (wfman or {}).get("registry")
    manifest = {
        # Default to now (UTC) rather than null. A reproducibility bundle whose
        # whole job is to make a run auditable must be able to say WHEN it ran;
        # recording null because the orchestrator forgot the flag is a gap, not
        # an honest absence.
        "timestamp": a.timestamp or _utc_now(),
        "skill": skill_info,
        "registry": reg,
        "workflow": {k: (wfman or {}).get(k) for k in
                     ("id", "name", "path", "engine", "fasta", "de", "validated")} if wfman else None,
        "query": {"acquisition": a.acquisition, "organism_taxid": a.organism_taxid,
                  "instrument": a.instrument},
        "engine": a.engine,
        # `ran`: whether a DE ran at all (de_ran_of), and from what that was read
        "de": {"method": a.de_method, "contrasts": de_contrasts,
               "q_cutoff": a.q_cutoff, "logfc": a.logfc, "adjp": a.adjp,
               "ran": de_ran, "ran_source": de_ran_why},
        # the search's own record, which step 5 replays and REPRODUCE.md describes
        "search_record": sp_path if sp is not None else None,
        "environment": {"os": platform.platform(), "python": py, "rscript": rscript,
                        "conda": conda, "env_prefix": env_prefix,
                        # where rscript / env_prefix come from: the DE's own record, or setup.json
                        "source": env_source,
                        "de_runtime": rt, "setup_json_differs": setup_differs,
                        "tool_versions": versions},
        "inputs": {"raw": [f.rstrip("/") for f in raw_files],
                   "fasta": a.fasta, "fasta_info": fasta_info,
                   "conditions": a.conditions, "params": a.params},
        # qc_bracket.py's record: where and which (sha256), never its verdict; null = not given
        "instrument_qc": instrument_qc,
        # skyline_peaks.py's record (step 8f): the tool, container sha256, tolerance + source, N
        "peak_check": peak_check,
        "checksums_file": "checksums/checksums.json",
        "files_in_bundle": "see MANIFEST.txt",
    }
    open(os.path.join(out, "run_manifest.json"), "w").write(json.dumps(manifest, indent=2)); ok("run_manifest.json")

    # ---- reproduce.sh --------------------------------------------------------
    commit = (reg or {}).get("commit") or "main"
    wf_id = (wfman or {}).get("id", "<workflow-id>")
    # Parameters now ship with the skill; a run record names the defaults table
    # version rather than a commit in a repo that no longer holds them.
    defaults_version = (reg or {}).get("defaults_version") or "<defaults_version>"
    acq_repro = (wfman or {}).get("acquisition") or "<DIA|DDA>"
    _instr = ((wfman or {}).get("instruments") or [None])[0]
    instr_repro = shlex.quote(_instr) if _instr else "'<instrument>'"
    engine_repro = ((wfman or {}).get("engine") or {}).get("name") or "<engine>"
    # Replay the ORIGINAL invocation, not a subset of it. Omitting these dropped
    # the organism and the platform block from the regenerated manifest: the
    # engine params still came back byte-identical (they do not depend on either),
    # but the run RECORD lost the species -- the single most important contextual
    # fact about a proteomics run -- and could not say what machine it ran on.
    _tax = (wfman or {}).get("organism_taxid") or a.organism_taxid
    tax_repro = f" \\\n  --organism-taxid {int(_tax)}" if _tax else ""
    # The Orbitrap resolution decides the instrument class and so the mass accuracy: a replay
    # without it reclassifies a 60k/15k Lumos as orbitrap_generic and derives different
    # tolerances. Replayed with the SOURCE it had (detected / user / cfg), so the rebuilt
    # manifest labels the numbers the same way.
    _res = (wfman or {}).get("resolution") or {}
    _res_flags = [f"--{lvl}-resolution {int(float(_res[lvl]))}" for lvl in ("ms1", "ms2")
                  if _res.get(lvl)]
    if _res.get("ms2_analyzer") in ("ITMS", "mixed"):
        # an ion-trap (or partly ion-trap) MS2 has no single MS2 resolution; the analyzer is
        # what selects its tolerances
        _res_flags.append(f"--ms2-analyzer {_res['ms2_analyzer']}")
    if _res_flags and _res.get("source"):
        _res_flags.append(f"--resolution-source {shlex.quote(str(_res['source']))}")
    res_repro = (" \\\n  " + " ".join(_res_flags)) if _res_flags else ""
    raw_arg = " ".join(f"'{f}'" for f in raw_files) or "/path/to/raw/*"
    # A run whose file name holds '--' was searched under a new name (run_search.py
    # --rename-dashed-runs; search_provenance.json `input_links`): the replay must rename it the
    # same way, or it is refused before it starts -- and the design names the new run.
    # A Sage .d converted in a container other than the pin's (run_search.py --pwiz-sif,
    # `pwiz_sif_arg`): the replay converts with the same one.
    # Both read from the search's own record (search_record_of: --search-record, else the one
    # beside --report); none -> neither flag.
    rename_repro = ""
    _sp = sp or {}
    _links = _sp.get("input_links") if isinstance(_sp.get("input_links"), list) else []
    if any(isinstance(x, dict) and x.get("renamed") for x in _links):
        rename_repro = " --rename-dashed-runs"
    if isinstance(_sp.get("pwiz_sif_arg"), list) and _sp["pwiz_sif_arg"]:
        rename_repro += " " + " ".join(shlex.quote(str(x)) for x in _sp["pwiz_sif_arg"])

    # Rebuild the FASTA from what actually ran (fetch_fasta.py's output), falling back
    # to the bundle only if --fasta-info-file was not supplied. The bundle's proteome is the
    # workflow DEFAULT — the user may have confirmed a different organism in step 3, and
    # regenerating from the bundle would silently reproduce a different database.
    # `or {}` on every lookup, not just the outer one: a manifest with an explicit
    # "fasta": null used to raise AttributeError here — after run_manifest.json was
    # written but before MANIFEST.txt/reproduce.sh, leaving a silently partial bundle.
    fi = fasta_info if isinstance(fasta_info, dict) else {}
    bundle_fasta = (wfman or {}).get("fasta") or {}
    fasta_repro_proteome = fi.get("proteome") or bundle_fasta.get("uniprot_proteome") or "<PROTEOME>"
    fasta_repro_content = fi.get("content_used") or fi.get("content_requested") or "one_per_gene"
    # Derive the contaminant set from what the FASTA actually contains. `contaminant_set`
    # is absent in older fasta_info payloads, and falling through to the bundle flag then
    # emitted `--contaminants none` for a database that demonstrably had 381 of them.
    if fi.get("contaminant_set"):
        fasta_repro_contam = fi["contaminant_set"]
    elif fi.get("n_contaminants_appended"):
        fasta_repro_contam = "universal"
    elif fi:
        fasta_repro_contam = "none"
    else:
        fasta_repro_contam = "universal" if bundle_fasta.get("add_contaminants", True) else "none"
    if fasta_repro_contam == "already_in_supplied_database":
        # Contaminants came in with the supplied database, not from a set we can name.
        fasta_repro_contam = "none"
    if fi:
        fasta_repro_note = (
            f"{fi.get('organism') or 'organism not recorded'} "
            f"(taxid {fi.get('taxid') or '?'}), {fi.get('n_proteome', '?')} proteome + "
            f"{(fi.get('n_contaminants_appended') or 0) + (fi.get('n_contaminants_already_present') or 0)}"
            f" contaminant sequences, "
            f"UniProt release {fi.get('uniprot_release') or 'unrecorded'}. "
            f"Re-running now uses the CURRENT release, so counts may differ slightly; "
            f"the searched FASTA's sha256 is in checksums/.")
        # A HIVE pre-staged copy has no knowable release (fetch_fasta.py leaves it empty on
        # purpose), so name the copy itself: its date and hash are what identify it.
        sf = fi.get("staged_file") or {}
        if sf:
            fasta_repro_note += (f" Built from a pre-staged copy ({sf.get('path', '?')}, "
                                 f"dated {(sf.get('mtime_utc') or '?')[:10]}, sha256 "
                                 f"{(sf.get('sha256') or '?')[:12]}...), release unknown.")
    else:
        fasta_repro_note = ("NOT RECORDED — --fasta-info-file was not passed to provenance.py, "
                            "so these values come from the workflow bundle's defaults and "
                            "may not be what was searched. Verify against checksums/.")
    # The digestion enzyme(s) decide which protease contaminant entries stay Cont_ when they
    # match a target protein, so a non-trypsin digest rebuilt with the default would produce a
    # different database. Sidecars from before --enzyme existed lack the field; their build was
    # the default, which is what the fallback names. An empty list is --enzyme none (L11 N3).
    fasta_repro_enzyme = shlex.quote(enzyme_arg(fi))
    # fetch_fasta.py now removes contaminant entries identical to a target protein (bovine
    # ACTB = human ACTB ...). A sidecar with no contaminant_target_rule was written BEFORE that
    # check, so its database still holds them; replayed as-is, today's fetch_fasta.py would
    # build a DIFFERENT database (153 fewer human Cont_ entries). Replay it faithfully -- and
    # a sidecar built with the check disabled (a replay of a replay) likewise.
    _rule = fi.get("contaminant_target_rule")
    fasta_repro_keep = ""
    # Which rule built the database: fetch_fasta.sidecar_state(), the one definition (legacy /
    # identity_only / current). A database this replay rebuilds without contaminants has none.
    _state = sidecar_state(fi) if fi and fasta_repro_contam != "none" else "current"
    if _state == "legacy" or (fi and fasta_repro_contam != "none"
                              and _rule == KEEP_TARGET_CONTAMINANTS_RULE):
        fasta_repro_keep = " --keep-target-contaminants"
        fasta_repro_note += (
            " The original database was built before target-identical contaminants were "
            "removed; this replays that faithfully -- drop --keep-target-contaminants to get "
            "the corrected database." if not _rule else
            " The original database was built with --keep-target-contaminants; this replays "
            "that faithfully -- drop the flag to get the corrected database.")
    # Built by the identity rule alone (sidecar_state "identity_only": the rule, no
    # min_unique_peptides -- before the peptide rule, so near-identical entries such as bovine
    # EEF1A1 vs mouse stayed). Today's fetch_fasta.py would drop more; --min-unique-peptides 0
    # rebuilds that database. A recorded threshold other than the default (0 included) is
    # replayed as recorded.
    _k = fi.get("min_unique_peptides")
    if _state == "identity_only" and _k is None:
        fasta_repro_keep += " --min-unique-peptides 0"
        fasta_repro_note += (
            " The original database was built before near-identical contaminants were "
            "removed; this replays that faithfully -- drop --min-unique-peptides 0 to get "
            "the corrected database.")
    elif (_state != "legacy" and _k is not None and _k != MIN_UNIQUE_PEPTIDES and
          _rule != KEEP_TARGET_CONTAMINANTS_RULE and fasta_repro_contam != "none"):
        fasta_repro_keep += f" --min-unique-peptides {int(_k)}"
    # A keratin sample's database (fetch --keratin-sample) has no keratin-family contaminant
    # entries; rebuilt without the flag it would have them, and run_search.py --keratin-sample
    # would refuse it. The one reader of the field is fetch_fasta.keratin_sample_recorded().
    fasta_repro_keratin = keratin_sample_recorded(fi) is True
    if fasta_repro_keratin:
        fasta_repro_keep += " --keratin-sample"
        fasta_repro_note += (" A keratin sample: its keratin-family contaminant entries were "
                             "removed (--keratin-sample).")
    elif keratin_sample_recorded(fi) is False and fi.get("keratin_sample_source") == "user":
        # the user's "not keratin" is replayed as an answer, never left to the default
        fasta_repro_keep += " --no-keratin-sample"
    # Contaminant entries the user left out by name (fetch --drop-contaminant), as given.
    _named = [str(x) for x in (fi or {}).get("drop_contaminant_requested") or [] if x]
    for acc in _named:
        fasta_repro_keep += f" --drop-contaminant {shlex.quote(acc)}"
    if _named:
        fasta_repro_note += (f" Contaminant entr{'y' if len(_named) == 1 else 'ies'} "
                             f"{', '.join(_named)} {'was' if len(_named) == 1 else 'were'} left "
                             f"out by name (--drop-contaminant).")
    # The user's own target sequences (fetch --add-fasta), from the files the sidecar records.
    for f in (fi or {}).get("added_sequences") or []:
        fasta_repro_keep += f" --add-fasta {shlex.quote(f.get('file') or '<added sequences file>')}"
        fasta_repro_note += (f" Target sequences were added from {f.get('file') or '?'} "
                             f"({f.get('n_entries', '?')} entries, sha256 "
                             f"{(f.get('sha256') or '?')[:12]}...): check that file is unchanged.")
    if fasta_repro_content in ("unknown", "as_staged"):
        # A --path override or a HIVE-staged file: not reconstructible from a proteome ID.
        # fetch_fasta.py's entry-count check (content_inferred) is the best guess at what a
        # staged file was; it only ever names a real --content value, so it is safe here.
        inferred = (fi or {}).get("content_inferred")
        fasta_repro_content = inferred if inferred in ("one_per_gene", "full",
                                                       "full_isoforms") else "one_per_gene"
        fasta_repro_note += (" NOTE: the original FASTA was supplied directly (override or "
                             "pre-staged), not downloaded — this command approximates it.")

    # the quantities the DE read (normalised / non-normalised, normalization_check.py's decision):
    # a replay without them would quietly read the other ones
    de_quantities = (f" \\\n  --quantities {de_rec['quantities']}"
                     if de_rec.get("quantities") in ("normalised", "raw") else "")
    # the covariates the DE fitted as numbers (run_de.R --numeric-covariates, its `covariates`
    # record): replayed without them each would be a factor -- another model, or no fit at all
    _cov = de_rec.get("covariates") if isinstance(de_rec.get("covariates"), dict) else {}
    de_numeric = [c for c in (_cov.get("numeric_requested") or []) if c in ("Covariate1", "Covariate2")]
    de_numeric_arg = f" \\\n  --numeric-covariates {','.join(de_numeric)}" if de_numeric else ""
    # the MaxLFQ route's between-run step THAT RAN (run_de.R --between-run, build_maxlfq.R's
    # between_run.method), always, for normalised quantities: a replay leaning on the replaying
    # skill's auto would silently run another step whenever auto changes (2.14.2 changed it once).
    # A maxlfq DE from before 2.14.2 records none: normalization_check.between_run_of's legacy rule.
    import normalization_check
    _br_flag = (normalization_check.between_run_of(de_rec, "method")
                if de_rec.get("quantities") != "raw" else None)
    de_between_arg = f" \\\n  --between-run {_br_flag}" if _br_flag else ""
    # A route with a between-run step of its own (its record's `between_run`, read by
    # between_run_of -- never by the route's name, rule 1: MaxLFQ) on raw quantities read a
    # --no-norm report's PG.MaxLFQ. run_search.py has no --no-norm, so the replayed search below
    # writes a NORMALISED report, on which run_de.R refuses --quantities raw: said here plainly,
    # with what to do, rather than left to fail at step 6.
    _own_step = normalization_check.between_run_of(de_rec, "method") is not None
    _route = de_rec.get("method") or a.de_method or "<method>"
    de_nonorm_note = ""
    if _own_step and de_rec.get("quantities") == "raw":
        de_nonorm_note = (
            f"\n#    NOTE: the DE read a --no-norm report's PG.MaxLFQ (--method {_route} --quantities "
            "raw: no\n#    between-run normalisation). This search replay passes no --no-norm, so "
            "its report is\n#    normalised and step 6 would be refused. Re-quantify with "
            "--no-norm instead -- diann_parallel.py\n#    --no-norm (it writes "
            "no_norm_report.parquet), or --no-norm in a single-shot cfg -- and point\n#    "
            "step 6's --input at that report.\n"
            "echo 'reproduce.sh: the DE read a --no-norm report: re-quantify with --no-norm "
            "before step 6 (see the NOTE above)' >&2")
    elif _own_step and isinstance(de_rec.get("dia_nn_normalisation"), dict) \
            and de_rec["dia_nn_normalisation"].get("no_norm") is True:
        # normalised quantities of a --no-norm report: the DE's own --between-run step was its only
        # between-run normalisation; on the replayed (DIA-NN-normalised) report it would sit on
        # top of DIA-NN's -- the double normalisation 2.14.2 removed (I20). Not refused: said.
        _had = ("no between-run normalisation at all" if _br_flag == "none" else
                f"--between-run {_br_flag} as its only between-run normalisation")
        _gets = ("DIA-NN's normalisation instead" if _br_flag == "none" else
                 "that step on top of DIA-NN's")
        de_nonorm_note = (
            f"\n#    NOTE: the DE read a --no-norm report's PG.MaxLFQ (--method {_route} --quantities "
            f"normalised, with\n#    {_had}). This search replay passes no --no-norm,\n#    so "
            f"its report is DIA-NN-normalised and step 6 would carry {_gets} --\n#    another "
            "analysis. Re-quantify with --no-norm instead -- diann_parallel.py --no-norm (it\n#    "
            "writes no_norm_report.parquet), or --no-norm in a single-shot cfg -- and point step "
            "6's\n#    --input at that report.\n"
            "echo 'reproduce.sh: the DE read a --no-norm report: re-quantify with --no-norm "
            "before step 6 (see the NOTE above)' >&2")
    # Step 5's parameters: the ones the search RAN with, from inputs/ (replay_params) -- step 2's
    # defaults for the data type know none of the user's overrides, an identifications-only cfg,
    # a non-trypsin enzyme or SILAC labels (estimate_params.py --labels). Its notes are comments
    # on step 5's heading line: placed after a line's trailing backslash, a comment ends the
    # continuation, and bash then runs the next line (--bundle ...) as a command of its own.
    search_params, search_params_note, (_lvl, _said) = replay_params(
        a.params, sp, a.engine or engine_repro, os.path.join(out, "inputs"))
    if _lvl == "OK":
        ok(_said)
    else:
        skip("search parameters", _said)
    side = load_json(a.params + ".rationale.json") if a.params else None
    if isinstance(side, dict) and isinstance(side.get("labels"), dict):
        search_params_note += ("\n#    They carry the SILAC labels (estimate_params.py --labels), which "
                               "step 2's defaults for the data\n#    type do not.")
    # Where it runs: run_search.py refuses an inline search on a cluster login node, and a search
    # for identifications only writes identifications.tsv (after the job: --adapt-only)
    search_params_note += ("\n#    On a SLURM cluster's login node run_search.py refuses to search "
                           "inline: add --sbatch ./search_job.sh\n#    and submit what it prints"
                           + ("; when the job ends, run this line again with --adapt-only."
                              if ids_only else "."))
    if ids_only:
        search_params_note += ("\n#    Identifications only (its record says it did not quantify): "
                               "the result is ./search_out/identifications.tsv,\n#    not a "
                               "report.parquet, and no DE follows.")
    # Step 6 only when a DE ran (de_ran_of); the last line names what to compare
    done_line = ("Compare ./de_results against checksums/checksums.json to confirm reproduction."
                 if de_ran else "Compare ./search_out against checksums/checksums.json to "
                                "confirm reproduction.")
    # Step 6: a DE pipeline that records its own command (de_provenance.json `replay`, e.g.
    # run_silac.R) is replayed with it -- its paths pointed at this bundle's -- else run_de.R.
    rp = de_rec.get("replay") if isinstance(de_rec.get("replay"), dict) else None
    if rp and rp.get("script") and isinstance(rp.get("args"), list):
        here = {"--input": "./search_out/report.parquet", "--outdir": "./de_results",
                "--figures-dir": "./de_results/figures", "--fasta": "./search.fasta"}
        words, args, i = [], [str(x) for x in rp["args"]], 0
        while i < len(args):
            flag = args[i]
            val = args[i + 1] if i + 1 < len(args) and not args[i + 1].startswith("--") else None
            if val is None:
                words.append(flag)
                i += 1
                continue
            if flag in here:
                val = here[flag]
            elif (flag in ("--design", "--metadata") and de_design_name and
                  (val == de_design or (os.path.isfile(val) and os.path.samefile(val, de_design)))):
                val = f"inputs/{de_design_name}"
            elif flag in ("--design", "--metadata", "--params"):
                val = f"inputs/{os.path.basename(val)}"
            words += [flag, val]
            i += 2
        de_step = (f'Rscript "$SKILL/scripts/{os.path.basename(rp["script"])}" '
                   + " ".join(shlex.quote(w) for w in words))
    else:
        de_step = None
    # where the DE ran (run_de.R's runtime record), so the re-run can be put in the same place
    de_where = ""
    if rt and rt.get("container"):
        de_where = (f"\n#    The DE ran inside the container {rt['container']} "
                    f"({rt.get('r_version')}): run this step inside it (apptainer exec ...).")
    elif rt:
        de_where = (f"\n#    The DE ran in {rt.get('r_version')} at {rt.get('r_home')}"
                    + (f" (conda env {env_prefix})" if env_prefix else "") + ".")
    # the run_de.R replay, exactly as reproduce.sh has always written it
    contrasts_arg = ('--contrasts "' + de_contrasts + '"') if de_contrasts else ""
    run_de_step = ('Rscript "$SKILL/scripts/run_de.R" --input ./search_out/report.parquet \\\n'
                   "  --metadata inputs/"
                   + (os.path.basename(a.conditions) if a.conditions else "conditions.csv") + " \\\n"
                   "  --method " + (a.de_method or "<method>") + " --outdir ./de_results \\\n"
                   "  " + contrasts_arg + " \\\n"
                   "  --q-cutoff " + (a.q_cutoff or "0.01") + " --logfc " + (a.logfc or "1.0")
                   + " --adjp " + (a.adjp or "0.05") + de_between_arg + de_quantities
                   + de_numeric_arg)
    de_section = (f"\n# 6. Re-run differential expression with identical settings.{de_where}\n"
                  f"{de_step or run_de_step}\n" if de_ran else
                  "\n# No differential expression ran in this analysis, so there is no step 6.\n")
    repro = f"""#!/usr/bin/env bash
# Auto-generated by provenance.py — re-creates this exact analysis.
# Requires: this skill's scripts/ on $SKILL, and internet for the registry + UniProt.
set -euo pipefail
SKILL="${{SKILL:?set SKILL to the ucdavis-proteomics-core-pipeline skill dir}}"

# 1. Recreate the analysis environment from the exact lock (byte-identical packages).
if [ -f environment/conda-explicit.txt ] && command -v micromamba >/dev/null 2>&1; then
  micromamba create -y -n proteomics-pipeline-repro --file environment/conda-explicit.txt
  micromamba activate proteomics-pipeline-repro
else
  echo "Run \\$SKILL/scripts/setup.sh to build the environment, then re-run."; bash "$SKILL/scripts/setup.sh"
  source ~/.proteomics-pipeline/activate.sh
fi

# 2. Re-derive the search defaults from the data type. These ship with the skill, so
#    the same skill version reproduces them exactly -- nothing is fetched.
#    Original defaults_version: {defaults_version}
bash "$SKILL/scripts/detect_env.sh" > ./env.json
python3 "$SKILL/scripts/resolve_defaults.py" --acquisition {acq_repro} \\
  --instrument {instr_repro} --engine {engine_repro}{tax_repro}{res_repro} \\
  --env ./env.json --dest ./wf

# 3. Resolve the same engine + version.
PIN_ENGINE={a.engine or '<engine>'} PIN_VERSION={(wfman or {}).get('engine',{}).get('version','')} \\
  bash "$SKILL/scripts/acquire_tools.sh" "$(bash "$SKILL/scripts/detect_env.sh" | python3 -c 'import sys,json;print(json.load(sys.stdin)["platform_class"])')"

# 4. Rebuild the FASTA — same proteome, database type, and contaminant set that were
#    ACTUALLY searched (from fetch_fasta.py's output), not the workflow bundle's default.
#    Original: {fasta_repro_note}
python3 "$SKILL/scripts/fetch_fasta.py" fetch --proteome {fasta_repro_proteome} \\
  --content {fasta_repro_content} --contaminants {fasta_repro_contam} --enzyme {fasta_repro_enzyme}{fasta_repro_keep} --out ./search.fasta

# 5. Re-run the search (inputs from inputs/checksums; verify against checksums/checksums.json).{search_params_note}{de_nonorm_note}
python3 "$SKILL/scripts/run_search.py" --tools ~/.proteomics-pipeline/tools/tools.json \\
  --bundle ./wf/workflow.manifest.json --params {search_params} \\
  --fasta ./search.fasta --out ./search_out --files {raw_arg}{' --keratin-sample' if fasta_repro_keratin else ''}{rename_repro}
{de_section}
echo "Done. {done_line}"
"""
    rp = os.path.join(out, "reproduce.sh")
    open(rp, "w").write(repro); os.chmod(rp, 0o755); ok("reproduce.sh")

    # ---- human-readable REPRODUCE.md ----------------------------------------
    methods = ""
    mt = os.path.join(a.de_dir or "", "methods.txt")
    if a.de_dir and os.path.exists(mt):
        methods = open(mt).read()
        ok("included DE methods.txt")
    # The plain-R transcript of the DE, written by run_de.R. It is the thing most
    # readers actually want, so REPRODUCE.md leads with it rather than with the
    # environment lock. Only advertise it if it is really there.
    repro_r = os.path.join(a.de_dir or "", "reproducibility_log.R")
    if a.de_dir and os.path.exists(repro_r):
        ok("found the DE's plain-R transcript (reproducibility_log.R)")
        r_section = f"""## Just show me the code

The whole differential-expression analysis, as plain R with every value written
out literally, is **`{os.path.relpath(repro_r, out)}`**. Read it to see exactly what
was done. To re-run it you need only R and limpa/limma — no conda environment, no
skill install, no shell scripts:

```
Rscript reproducibility_log.R
```

That covers the DE. Everything below is for reproducing the **search** as well —
the pinned engine build, the exact software environment, and the checksums that
make input drift visible.
"""
    else:
        if de_ran:                              # no DE ran: nothing missing (the NOTE says so)
            skip("reproducibility_log.R", "not found in --de-dir (run run_de.R to generate it)")
        r_section = ""

    # The machine the DE ran on (run_de.R's de_provenance.json `compute`): DPC-Quant's numbers
    # are exact only on the same CPU family -- say which, and how to get it again on HIVE.
    compute = {}
    dp = os.path.join(a.de_dir or "", "de_provenance.json")
    if a.de_dir and os.path.exists(dp):
        try:
            with open(dp, encoding="utf-8") as fh:
                compute = json.load(fh).get("compute") or {}
        except (OSError, ValueError) as e:
            skip("de_provenance.json compute", f"unreadable: {e}")
    # An analysis with no DE (de_ran_of) is described as what it was -- the search, how its inputs
    # went in, its identifications -- with none of the DE's sections (2.14.2 C2, E2E F2)
    cpu_section = compute_section(compute) if de_ran else ""
    rt_section = runtime_section(de_rec, rt, "; ".join(
        f"{k} {v['setup_json']}" for k, v in (setup_differs or {}).items())) if de_ran else ""
    what_section = analysis_section(sp, a.engine or engine_repro, de_ran, ids_only)
    de_line = (f"- DE: method=`{a.de_method}`, contrasts=`{de_contrasts}`, q≤{a.q_cutoff}, "
               f"adj.P<{a.adjp}\n  (significance is the BH adjusted p-value alone; "
               f"|logFC|={a.logfc} is a volcano reference line, not a filter)" if de_ran else
               f"- DE: none ran ({de_ran_why})")
    how_line = ("rebuilds the FASTA, and re-runs search + DE. Compare outputs to\n"
                "`checksums/checksums.json`. This is the heavyweight path — it re-runs a "
                "multi-hour\nsearch. If you only want the statistics, use the R script above."
                if de_ran else
                "rebuilds the FASTA, and re-runs the search with the parameters it ran with\n"
                "(`inputs/`). Compare its results to `checksums/checksums.json`. No DE ran, so "
                "there is\nno DE step.")
    methods_text = methods or ("(methods.txt not found — run run_de.R to generate it)" if de_ran
                               else "No differential expression ran, so there is no DE "
                                    "methods.txt. The publication Methods are the session's "
                                    "`output/methods.md` (make_methods.py).")

    head = f"{r_section}\n{rt_section}\n{cpu_section}" if de_ran else what_section
    # what checksums.json holds -- said from what was fingerprinted, never a fixed list
    _held = [w for k, w in (("raw_inputs", "the raw inputs"), ("fasta", "the FASTA"),
                            ("report", "the search report"),
                            ("search_outputs", "the search's outputs ("
                             + ", ".join(os.path.basename(x["path"])
                                         for x in checks.get("search_outputs") or []) + ")"),
                            ("search_record", "the search record"),
                            ("de_outputs", "the DE outputs")) if checks.get(k)]
    checksums_line = ("sha256 of " + (", ".join(_held[:-1]) + " and " + _held[-1]
                                      if len(_held) > 1 else _held[0])
                      if _held else "nothing was fingerprinted (see MANIFEST.txt)")
    md = f"""# Reproducibility — {(wfman or {}).get('name', 'proteomics analysis')}

{head}
## Skill that produced this
- **{skill_info['title']}** — `{skill_info['name']}` {skill_label(skill_info['version'])}
- Repository: {skill_info['repository']}
- This analysis was run by the above Claude skill (in Claude Code / Claude Desktop).
- Installed with:
  ```
  {skill_info['install'][0]}
  {skill_info['install'][1]}
  ```
- A copy of the exact skill scripts that ran is in the session's `scripts/` folder.

## Validated workflow
- id: `{wf_id}`
- registry: {REGISTRY_LINE(reg)}
- engine: {(wfman or {}).get('engine')}
{de_line}
- query: acquisition=`{a.acquisition}`, organism_taxid=`{a.organism_taxid}`, instrument=`{a.instrument}`

## How to reproduce the whole thing, search included
```
SKILL=/path/to/ucdavis-proteomics-core-pipeline bash reproduce.sh
```
`reproduce.sh` rebuilds the conda env from `environment/conda-explicit.txt`,
re-derives the search defaults from the data type (they ship with the skill version above;
defaults table `{defaults_version}`), re-resolves the engine,
{how_line}

## What's captured
- `run_manifest.json` — full machine-readable record
- `environment/` — conda lock, pip freeze, R sessionInfo (all package versions), tool versions
- `inputs/` — the exact params file, conditions.csv, workflow manifest{', commands.log' if a.commands else ''}
- `checksums/` — {checksums_line}

## Methods
{methods_text}

## Capture log
See `MANIFEST.txt` for exactly what was and wasn't captured.
"""
    open(os.path.join(out, "REPRODUCE.md"), "w", encoding="utf-8").write(md); ok("REPRODUCE.md")

    open(os.path.join(out, "MANIFEST.txt"), "w").write(
        "Reproducibility bundle — capture log\n" + "=" * 40 + "\n" + "\n".join(MANIFEST_LINES) + "\n")

    n_skip = sum(1 for l in MANIFEST_LINES if l.startswith("[SKIPPED]"))
    print(json.dumps({"bundle": out, "captured": sum(1 for l in MANIFEST_LINES
                                                     if l.startswith("[OK]")),
                      "skipped": n_skip, "manifest": os.path.join(out, "MANIFEST.txt"),
                      "reproduce": rp}, indent=2))


def REGISTRY_LINE(reg):
    if not reg:
        return "(not recorded)"
    if not reg.get("repo"):
        # the defaults ship with the skill (resolve_defaults.py: no remote registry), so there is
        # no repo or commit -- "None @ commit `None` (None)" said nothing (2.14.2 C2, E2E)
        return (f"{reg.get('source') or 'no remote registry'} (defaults table "
                f"`{reg.get('defaults_version') or 'not recorded'}`)")
    return f"{reg.get('repo')} @ commit `{reg.get('commit')}` ({reg.get('tree_url')})"


if __name__ == "__main__":
    main()
