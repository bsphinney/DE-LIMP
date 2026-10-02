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
    --fasta-info   "$(cat ./search.fasta.meta.json)" \
    --raw          /data/*.d \
    --report       ./search_out/report.parquet \
    --de-dir       ./de_results \
    --engine diann --de-method dpc --contrasts "B-A,C-A" \
    --q-cutoff 0.01 --logfc 1.0 --adjp 0.05 \
    --organism-taxid 9606 --instrument "Orbitrap Astral" --acquisition DIA \
    --commands ./commands.log         # optional: a log of the exact commands run
    --qc-bracket <session>/logs/qc_bracket.json   # optional: the instrument's QC record
                                      # (qc_bracket.py) -- referenced by path + sha256 only

Outputs under --outdir:
  run_manifest.json          everything, machine-readable
  REPRODUCE.md               human-readable methods + how-to-rerun
  reproduce.sh               re-creates env, re-derives the shipped defaults, re-runs
  MANIFEST.txt               [OK]/[SKIPPED] log of what was captured
  environment/               conda-explicit.txt, pip-freeze.txt, r-sessionInfo.txt, versions.txt
  inputs/                    copies of params, conditions.csv, the workflow manifest
  checksums/                 sha256 of inputs, report, and DE outputs
"""
import sys, os, json, glob, shutil, hashlib, argparse, subprocess, platform, shlex

# one definition each, where it is written
from fetch_fasta import (KEEP_TARGET_CONTAMINANTS_RULE, MIN_UNIQUE_PEPTIDES, sidecar_state,
                         keratin_sample_recorded)
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


def conda_prefix_of_r_home(r_home):
    """A conda env's R lives at <prefix>/lib/R: that prefix when it is a conda env, else None."""
    if not r_home:
        return None
    cand = os.path.dirname(os.path.dirname(os.path.normpath(r_home)))
    return cand if os.path.isdir(os.path.join(cand, "conda-meta")) else None


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
    ap.add_argument("--fasta-info", help="JSON from fetch_fasta.py")
    ap.add_argument("--raw", nargs="*")
    ap.add_argument("--report")
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
    ap.add_argument("--timestamp", default="", help="ISO timestamp (script can't call clock)")
    a = ap.parse_args()

    out = os.path.abspath(a.outdir)
    for sub in ("environment", "inputs", "checksums"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)

    setup = load_json(a.setup_json) if a.setup_json else None
    tools = load_json(a.tools_json) if a.tools_json else None
    wfman = load_json(a.workflow_manifest) if a.workflow_manifest else None
    fasta_info = None
    if a.fasta_info:
        try: fasta_info = json.loads(a.fasta_info)
        except json.JSONDecodeError: skip("fasta-info", "not valid JSON")

    # ---- environment capture -------------------------------------------------
    # The environment the DE RAN in, as run_de.R recorded it (de_provenance.json `runtime`):
    # its Rscript and conda env. setup.json says only which env setup.sh built, and a DE run in
    # another (a separate R 4.6 env for limpa 1.4, staff report 2026-09-28) made the bundle name
    # the wrong R and limpa. setup.json / PATH only when the DE recorded no runtime -- said so.
    env_dir = os.path.join(out, "environment")
    de_rec, rt = de_runtime(a.de_dir)
    py = (setup or {}).get("python") or shutil.which("python3") or sys.executable
    conda = (setup or {}).get("conda") or shutil.which("micromamba") or shutil.which("mamba") or shutil.which("conda")
    setup_rscript, setup_prefix = (setup or {}).get("rscript"), (setup or {}).get("env_prefix")
    setup_differs = None
    if rt:
        rscript = rt.get("rscript")
        env_prefix = rt.get("conda_prefix") or conda_prefix_of_r_home(rt.get("r_home"))
        env_source = ("de_provenance.json runtime: recorded by run_de.R in the process that ran "
                      "the DE")
        ok(f"environment: the one the DE ran in ({rt.get('r_version')}, {rt.get('r_home')})")
        def _same(x, y):
            return bool(x and y) and os.path.realpath(x) == os.path.realpath(y)
        diff = {k: (sv, dv) for k, sv, dv in (("rscript", setup_rscript, rscript),
                                              ("env_prefix", setup_prefix, env_prefix))
                if sv and not _same(sv, dv)}
        if diff:
            setup_differs = {k: {"setup_json": sv, "de": dv} for k, (sv, dv) in diff.items()}
            note("setup.json names a different environment than the one the DE ran in ("
                 + "; ".join(f"{k}: setup.json {sv}, DE {dv}" for k, (sv, dv) in diff.items())
                 + "): the bundle describes the DE's")
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
           f"Rscript ({rscript})" + ("" if rt else " -- may not be the R the DE ran in"))
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
            "claude plugin marketplace add bsphinney/DE-LIMP",
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
    in_dir = os.path.join(out, "inputs")
    for label, src in (("params", a.params), ("conditions.csv", a.conditions),
                       ("workflow.manifest.json", a.workflow_manifest),
                       ("commands.log", a.commands)):
        if src and os.path.exists(src) and label == "commands.log":
            # every command verbatim -- but who decided is staff-only: names become the role
            import staff
            with open(src, encoding="utf-8", errors="replace") as fh:
                logged = fh.read()
            with open(os.path.join(in_dir, os.path.basename(src)), "w", encoding="utf-8") as fh:
                fh.write(staff.redact(logged))
            ok(f"copied {label} (names after --by/--override-by replaced by \"{staff.ROLE}\")")
        elif src and os.path.exists(src):
            shutil.copy2(src, os.path.join(in_dir, os.path.basename(src))); ok(f"copied {label}")
            # params estimated by estimate_params.py carry a sibling rationale — capture it
            if label == "params" and os.path.exists(src + ".rationale.json"):
                shutil.copy2(src + ".rationale.json",
                             os.path.join(in_dir, os.path.basename(src) + ".rationale.json"))
                ok("copied params rationale (per-setting provenance)")
            # the user's sample-identity answers behind the design (collect_conditions.py
            # --confirm-multi), when there were any
            if label == "conditions.csv" and os.path.exists(src + ".decisions.json"):
                shutil.copy2(src + ".decisions.json",
                             os.path.join(in_dir, os.path.basename(src) + ".decisions.json"))
                ok("copied conditions decisions (labels the user said name several runs)")
        elif src:
            skip(label, f"not found: {src}")

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
        "de": {"method": a.de_method, "contrasts": a.contrasts,
               "q_cutoff": a.q_cutoff, "logfc": a.logfc, "adjp": a.adjp},
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

    # Rebuild the FASTA from what actually ran (fetch_fasta.py's output), falling back
    # to the bundle only if --fasta-info was not supplied. The bundle's proteome is the
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
        fasta_repro_note = ("NOT RECORDED — --fasta-info was not passed to provenance.py, "
                            "so these values come from the workflow bundle's defaults and "
                            "may not be what was searched. Verify against checksums/.")
    # The digestion enzyme(s) decide which protease contaminant entries stay Cont_ when they
    # match a target protein, so a non-trypsin digest rebuilt with the default would produce a
    # different database. Sidecars from before --enzyme existed lack the field; their build was
    # the default, which is what the fallback names.
    fasta_repro_enzyme = shlex.quote(",".join((fi or {}).get("digestion_enzymes_used")
                                              or ["trypsin", "lysc"]))
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
    # maxlfq + raw read a --no-norm report's PG.MaxLFQ. run_search.py has no --no-norm, so the
    # replayed search below writes a NORMALISED report, on which run_de.R refuses --quantities raw:
    # said here plainly, with what to do, rather than left to fail at step 6.
    de_nonorm_note = ""
    if de_rec.get("method") == "maxlfq" and de_rec.get("quantities") == "raw":
        de_nonorm_note = (
            "\n#    NOTE: the DE read a --no-norm report's PG.MaxLFQ (--method maxlfq --quantities "
            "raw: no\n#    between-run normalisation). This search replay passes no --no-norm, so "
            "its report is\n#    normalised and step 6 would be refused. Re-quantify with "
            "--no-norm instead -- diann_parallel.py\n#    --no-norm (it writes "
            "no_norm_report.parquet), or --no-norm in a single-shot cfg -- and point\n#    "
            "step 6's --input at that report.\n"
            "echo 'reproduce.sh: the DE read a --no-norm report: re-quantify with --no-norm "
            "before step 6 (see the NOTE above)' >&2")
    # where the DE ran (run_de.R's runtime record), so the re-run can be put in the same place
    de_where = ""
    if rt and rt.get("container"):
        de_where = (f"\n#    The DE ran inside the container {rt['container']} "
                    f"({rt.get('r_version')}): run this step inside it (apptainer exec ...).")
    elif rt:
        de_where = (f"\n#    The DE ran in {rt.get('r_version')} at {rt.get('r_home')}"
                    + (f" (conda env {env_prefix})" if env_prefix else "") + ".")
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

# 5. Re-run the search (inputs from inputs/checksums; verify against checksums/checksums.json).{de_nonorm_note}
python3 "$SKILL/scripts/run_search.py" --tools ~/.proteomics-pipeline/tools/tools.json \\
  --bundle ./wf/workflow.manifest.json --params ./wf/$(basename "$(ls wf | grep -vi manifest | head -n1)") \\
  --fasta ./search.fasta --out ./search_out --files {raw_arg}{' --keratin-sample' if fasta_repro_keratin else ''}

# 6. Re-run differential expression with identical settings.{de_where}
Rscript "$SKILL/scripts/run_de.R" --input ./search_out/report.parquet \\
  --metadata inputs/{os.path.basename(a.conditions) if a.conditions else 'conditions.csv'} \\
  --method {a.de_method or '<method>'} --outdir ./de_results \\
  {('--contrasts "'+a.contrasts+'"') if a.contrasts else ''} \\
  --q-cutoff {a.q_cutoff or '0.01'} --logfc {a.logfc or '1.0'} --adjp {a.adjp or '0.05'}{de_quantities}

echo "Done. Compare ./de_results against checksums/checksums.json to confirm reproduction."
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
    cpu_section = compute_section(compute)
    rt_section = runtime_section(de_rec, rt, "; ".join(
        f"{k} {v['setup_json']}" for k, v in (setup_differs or {}).items()))

    md = f"""# Reproducibility — {(wfman or {}).get('name', 'proteomics analysis')}

{r_section}
{rt_section}
{cpu_section}
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
- DE: method=`{a.de_method}`, contrasts=`{a.contrasts}`, q≤{a.q_cutoff}, adj.P<{a.adjp}
  (significance is the BH adjusted p-value alone; |logFC|={a.logfc} is a volcano reference line, not a filter)
- query: acquisition=`{a.acquisition}`, organism_taxid=`{a.organism_taxid}`, instrument=`{a.instrument}`

## How to reproduce the whole thing, search included
```
SKILL=/path/to/ucdavis-proteomics-core-pipeline bash reproduce.sh
```
`reproduce.sh` rebuilds the conda env from `environment/conda-explicit.txt`,
re-derives the search defaults from the data type (they ship with the skill version above;
defaults table `{defaults_version}`), re-resolves the engine,
rebuilds the FASTA, and re-runs search + DE. Compare outputs to
`checksums/checksums.json`. This is the heavyweight path — it re-runs a multi-hour
search. If you only want the statistics, use the R script above.

## What's captured
- `run_manifest.json` — full machine-readable record
- `environment/` — conda lock, pip freeze, R sessionInfo (all package versions), tool versions
- `inputs/` — the exact params file, conditions.csv, workflow manifest{', commands.log' if a.commands else ''}
- `checksums/` — sha256 of raw inputs, FASTA, search report, and DE outputs

## Methods
{methods or '(methods.txt not found — run run_de.R to generate it)'}

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
    return f"{reg.get('repo')} @ commit `{reg.get('commit')}` ({reg.get('tree_url')})"


if __name__ == "__main__":
    main()
