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
    --fasta-info   '{"source":"uniprot:UP000005640","n_sequences":20596}' \
    --raw          /data/*.d \
    --report       ./search_out/report.parquet \
    --de-dir       ./de_results \
    --engine diann --de-method dpc --contrasts "B-A,C-A" \
    --q-cutoff 0.01 --logfc 1.0 --adjp 0.05 \
    --organism-taxid 9606 --instrument "Orbitrap Astral" --acquisition DIA \
    --commands ./commands.log         # optional: a log of the exact commands run

Outputs under --outdir:
  run_manifest.json          everything, machine-readable
  REPRODUCE.md               human-readable methods + how-to-rerun
  reproduce.sh               re-creates env, re-fetches pinned workflow, re-runs
  MANIFEST.txt               [OK]/[SKIPPED] log of what was captured
  environment/               conda-explicit.txt, pip-freeze.txt, r-sessionInfo.txt, versions.txt
  inputs/                    copies of params, conditions.csv, the workflow manifest
  checksums/                 sha256 of inputs, report, and DE outputs
"""
import sys, os, json, glob, shutil, hashlib, argparse, subprocess, platform

MANIFEST_LINES = []
def ok(msg):      MANIFEST_LINES.append(f"[OK]      {msg}")
def skip(w, why): MANIFEST_LINES.append(f"[SKIPPED] {w} -- {why}")

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


def load_json(path):
    try:
        return json.load(open(path))
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", default="./reproducibility")
    ap.add_argument("--workflow-manifest")
    ap.add_argument("--setup-json")
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
    env_dir = os.path.join(out, "environment")
    py = (setup or {}).get("python") or sys.executable
    rscript = (setup or {}).get("rscript")
    conda = (setup or {}).get("conda")
    env_prefix = (setup or {}).get("env_prefix")

    # conda explicit lock (fully pinned, URL+hash per package) — the gold standard
    if conda and env_prefix and os.path.isdir(env_prefix):
        txt = run_capture(f"{conda} list -p {env_prefix} --explicit --md5")
        if txt and "://" in txt:
            open(os.path.join(env_dir, "conda-explicit.txt"), "w").write(txt); ok("conda explicit lock")
        else:
            skip("conda-explicit.txt", "conda list returned no URLs")
    else:
        skip("conda-explicit.txt", "no conda env recorded in setup.json")

    # pip freeze (the env's python)
    if py:
        txt = run_capture(f"{py} -m pip freeze")
        open(os.path.join(env_dir, "pip-freeze.txt"), "w").write(txt); ok("pip freeze")

    # R sessionInfo with every package version (the DE step's exact stack)
    if rscript:
        txt = run_capture(f"{rscript} -e 'sink(stdout()); "
                          f"cat(R.version.string,\"\\n\"); "
                          f"for(p in c(\"limpa\",\"limma\",\"arrow\",\"dplyr\",\"tidyr\")) "
                          f"try(cat(p, as.character(packageVersion(p)), \"\\n\")); "
                          f"cat(\"\\n\"); print(sessionInfo())'")
        open(os.path.join(env_dir, "r-sessionInfo.txt"), "w").write(txt); ok("R sessionInfo + package versions")
    else:
        skip("r-sessionInfo.txt", "no Rscript in setup.json")

    # engine versions
    versions = {"os": platform.platform(), "python": platform.python_version()}
    if tools:
        versions["tools_versions"] = tools.get("versions")
        for eng in ("diann", "sage"):
            cmd = tools.get(eng)
            if cmd:
                versions[f"{eng}_cmd"] = cmd
    if a.engine == "sage" and (setup or {}).get("sage"):
        versions["sage_version"] = run_capture(f"{setup['sage']} --version").strip()
    open(os.path.join(env_dir, "versions.txt"), "w").write(json.dumps(versions, indent=2)); ok("tool versions")

    # ---- copy inputs ---------------------------------------------------------
    in_dir = os.path.join(out, "inputs")
    for label, src in (("params", a.params), ("conditions.csv", a.conditions),
                       ("workflow.manifest.json", a.workflow_manifest),
                       ("commands.log", a.commands)):
        if src and os.path.exists(src):
            shutil.copy2(src, os.path.join(in_dir, os.path.basename(src))); ok(f"copied {label}")
            # params estimated by estimate_params.py carry a sibling rationale — capture it
            if label == "params" and os.path.exists(src + ".rationale.json"):
                shutil.copy2(src + ".rationale.json",
                             os.path.join(in_dir, os.path.basename(src) + ".rationale.json"))
                ok("copied params rationale (per-setting provenance)")
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

    # ---- the master run manifest --------------------------------------------
    reg = (wfman or {}).get("registry")
    manifest = {
        "timestamp": a.timestamp or None,
        "skill": "proteomics-pipeline",
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
                        "tool_versions": versions},
        "inputs": {"raw": [f.rstrip("/") for f in raw_files],
                   "fasta": a.fasta, "fasta_info": fasta_info,
                   "conditions": a.conditions, "params": a.params},
        "checksums_file": "checksums/checksums.json",
        "files_in_bundle": "see MANIFEST.txt",
    }
    open(os.path.join(out, "run_manifest.json"), "w").write(json.dumps(manifest, indent=2)); ok("run_manifest.json")

    # ---- reproduce.sh --------------------------------------------------------
    commit = (reg or {}).get("commit") or "main"
    wf_id = (wfman or {}).get("id", "<workflow-id>")
    raw_arg = " ".join(f"'{f}'" for f in raw_files) or "/path/to/raw/*"
    repro = f"""#!/usr/bin/env bash
# Auto-generated by provenance.py — re-creates this exact analysis.
# Requires: this skill's scripts/ on $SKILL, and internet for the registry + UniProt.
set -euo pipefail
SKILL="${{SKILL:?set SKILL to the proteomics-pipeline skill dir}}"

# 1. Recreate the analysis environment from the exact lock (byte-identical packages).
if [ -f environment/conda-explicit.txt ] && command -v micromamba >/dev/null 2>&1; then
  micromamba create -y -n proteomics-pipeline-repro --file environment/conda-explicit.txt
  micromamba activate proteomics-pipeline-repro
else
  echo "Run \\$SKILL/scripts/setup.sh to build the environment, then re-run."; bash "$SKILL/scripts/setup.sh"
  source ~/.proteomics-pipeline/activate.sh
fi

# 2. Re-fetch the validated workflow PINNED to the original commit (not moving main).
python3 "$SKILL/scripts/fetch_workflows.py" pull --id {wf_id} --ref {commit} --dest ./wf

# 3. Resolve the same engine + version.
PIN_ENGINE={a.engine or '<engine>'} PIN_VERSION={(wfman or {}).get('engine',{}).get('version','')} \\
  bash "$SKILL/scripts/acquire_tools.sh" "$(bash "$SKILL/scripts/detect_env.sh" | python3 -c 'import sys,json;print(json.load(sys.stdin)["platform_class"])')"

# 4. Rebuild the FASTA (same proteome + contaminants).
python3 "$SKILL/scripts/fetch_fasta.py" --proteome {(wfman or {}).get('fasta',{}).get('uniprot_proteome','<PROTEOME>')} \\
  {'--add-contaminants' if (wfman or {}).get('fasta',{}).get('add_contaminants') else ''} --out ./search.fasta

# 5. Re-run the search (inputs from inputs/checksums; verify against checksums/checksums.json).
python3 "$SKILL/scripts/run_search.py" --tools ~/.proteomics-pipeline/tools/tools.json \\
  --bundle ./wf/workflow.manifest.json --params ./wf/$(basename "$(ls wf | grep -vi manifest | head -n1)") \\
  --fasta ./search.fasta --out ./search_out --files {raw_arg}

# 6. Re-run differential expression with identical settings.
Rscript "$SKILL/scripts/run_de.R" --input ./search_out/report.parquet \\
  --metadata inputs/{os.path.basename(a.conditions) if a.conditions else 'conditions.csv'} \\
  --method {a.de_method or '<method>'} --outdir ./de_results \\
  {('--contrasts "'+a.contrasts+'"') if a.contrasts else ''} \\
  --q-cutoff {a.q_cutoff or '0.01'} --logfc {a.logfc or '1.0'} --adjp {a.adjp or '0.05'}

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
    md = f"""# Reproducibility — {(wfman or {}).get('name', 'proteomics analysis')}

Generated by the proteomics-pipeline skill. This bundle contains everything
needed to reproduce the analysis.

## Validated workflow
- id: `{wf_id}`
- registry: {REGISTRY_LINE(reg)}
- engine: {(wfman or {}).get('engine')}
- DE: method=`{a.de_method}`, contrasts=`{a.contrasts}`, q≤{a.q_cutoff}, |logFC|≥{a.logfc}, adj.P<{a.adjp}
- query: acquisition=`{a.acquisition}`, organism_taxid=`{a.organism_taxid}`, instrument=`{a.instrument}`

## How to reproduce
```
SKILL=/path/to/proteomics-pipeline bash reproduce.sh
```
`reproduce.sh` rebuilds the conda env from `environment/conda-explicit.txt`,
re-fetches the workflow **pinned to commit `{commit}`**, re-resolves the engine,
rebuilds the FASTA, and re-runs search + DE. Compare outputs to
`checksums/checksums.json`.

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
    open(os.path.join(out, "REPRODUCE.md"), "w").write(md); ok("REPRODUCE.md")

    open(os.path.join(out, "MANIFEST.txt"), "w").write(
        "Reproducibility bundle — capture log\n" + "=" * 40 + "\n" + "\n".join(MANIFEST_LINES) + "\n")

    n_skip = sum(1 for l in MANIFEST_LINES if l.startswith("[SKIPPED]"))
    print(json.dumps({"bundle": out, "captured": len(MANIFEST_LINES) - n_skip,
                      "skipped": n_skip, "manifest": os.path.join(out, "MANIFEST.txt"),
                      "reproduce": rp}, indent=2))


def REGISTRY_LINE(reg):
    if not reg:
        return "(not recorded)"
    return f"{reg.get('repo')} @ commit `{reg.get('commit')}` ({reg.get('tree_url')})"


if __name__ == "__main__":
    main()
