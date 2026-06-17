#!/usr/bin/env python3
"""
fetch_workflows.py  --  The fetch/match layer (the keystone of the skill).

Pulls the compiled workflow index from the public DE-LIMP repo, hard-filters by
acquisition + organism (the two things that are never negotiable), scores the
remaining candidates by instrument, and returns the best match for the
orchestrator to CONFIRM with the user. No auth, stdlib only.

Reproducibility: every fetch resolves and records the registry's exact git
**commit SHA**. `pull` fetches params at that pinned commit (not a moving
branch), so re-running `pull --ref <sha>` later returns byte-identical
parameters. Always record the returned `registry.commit` in the run manifest.

Two modes:

  # 1. match: pick the validated workflow for this run
  python3 fetch_workflows.py match \
      --acquisition DIA --organism-taxid 9606 --instrument "Orbitrap Astral"
  #   -> {selected, candidates, needs_menu, registry, all} on stdout

  # 2. pull: download a chosen workflow's params + extra files, pinned to a commit
  python3 fetch_workflows.py pull --id diann_astral_dia_human --dest ./wf
  python3 fetch_workflows.py pull --id ... --ref <commit_sha> --dest ./wf  # exact re-run

Match key (locked, see PLAN.md §2): acquisition + instrument + organism_taxid.
Organism is a hard filter (it defines the FASTA); instrument is a tiebreaker.
"""
import sys, os, json, argparse, tempfile, urllib.request, urllib.error

REPO = "bsphinney/DE-LIMP"
RAW_BASE = f"https://raw.githubusercontent.com/{REPO}"
API_BASE = f"https://api.github.com/repos/{REPO}"


def raw_url(ref, relpath):
    return f"{RAW_BASE}/{ref}/{relpath}"


def _get(url, binary=False):
    req = urllib.request.Request(url, headers={"User-Agent": "proteomics-pipeline-skill"})
    with urllib.request.urlopen(req, timeout=30) as r:
        data = r.read()
    return data if binary else data.decode("utf-8")


def resolve_commit(ref="main"):
    """Resolve a ref (branch/tag/sha) to a concrete commit SHA for pinning."""
    try:
        data = json.loads(_get(f"{API_BASE}/commits/{ref}"))
        return data.get("sha")
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, ValueError):
        return None


def load_index(ref="main", cache_dir=None):
    """Fetch index.json at a ref (cache to temp so repeated calls are cheap)."""
    cache = os.path.join(cache_dir or tempfile.gettempdir(),
                         f"delimp_workflow_index_{ref}.json")
    url = raw_url(ref, "workflows/index.json")
    try:
        txt = _get(url)
        try:
            with open(cache, "w") as fh:
                fh.write(txt)
        except OSError:
            pass
        return json.loads(txt)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError) as e:
        if os.path.exists(cache):
            sys.stderr.write(f"[fetch_workflows] network failed ({e}); using cached index\n")
            return json.loads(open(cache).read())
        sys.exit(f"Could not fetch workflow index from {url} and no cache exists: {e}")


def registry_info(ref):
    commit = resolve_commit(ref)
    return {"repo": REPO, "ref": ref, "commit": commit,
            "tree_url": f"https://github.com/{REPO}/tree/{commit or ref}"}


def score_instrument(match, instrument):
    """exact model = 2; alias substring = 1; no instrument info = 0."""
    if not instrument:
        return 0
    instr = instrument.strip().lower()
    models = [m.lower() for m in match.get("instruments", [])]
    if instr in models or any(instr == m for m in models):
        return 2
    if any(instr == m or instr in m or m in instr for m in models):
        return 2
    for alias in match.get("instrument_aliases", []):
        if alias.lower() in instr:
            return 1
    return 0


def do_match(args):
    reg = registry_info(args.ref)
    idx = load_index(args.ref)
    acq = args.acquisition.upper()
    taxid = int(args.organism_taxid)
    instrument = args.instrument or ""

    pool = [w for w in idx["workflows"]
            if w["match"]["acquisition"].upper() == acq
            and int(w["match"]["organism_taxid"]) == taxid]

    scored = []
    for w in pool:
        s = score_instrument(w["match"], instrument)
        scored.append({"score": s, **_summary(w, reg)})
    scored.sort(key=lambda x: x["score"], reverse=True)

    selected = scored[0] if scored else None
    top = scored[0]["score"] if scored else None
    tie = len(scored) > 1 and scored[0]["score"] == scored[1]["score"]
    needs_menu = (len(scored) == 0) or tie or (top == 0)

    out = {
        "query": {"acquisition": acq, "organism_taxid": taxid, "instrument": instrument},
        "registry": reg,
        "selected": selected,
        "candidates": scored,
        "needs_menu": needs_menu,
        "reason": (
            "no validated workflow matches this acquisition+organism — add one to the registry"
            if not scored else
            "tie at the top score — ask the user which instrument/SOP applies" if tie else
            "no instrument info — confirm the auto-pick with the user" if top == 0 else
            "confident auto-pick — still confirm with the user before a multi-hour search"
        ),
    }
    print(json.dumps(out, indent=2))


def _summary(w, registry=None):
    s = {
        "id": w["id"], "name": w["name"], "path": w["path"],
        "acquisition": w["match"]["acquisition"],
        "instruments": w["match"].get("instruments", []),
        "organism": w["match"].get("organism"),
        "organism_taxid": w["match"]["organism_taxid"],
        "engine": w["engine"], "fasta": w.get("fasta", {}),
        "de": w["de"], "search": w["search"], "validated": w.get("validated", {}),
    }
    if registry:
        s["registry"] = registry
    return s


def do_pull(args):
    reg = registry_info(args.ref)
    # Pin to the resolved commit so files are byte-stable across re-runs.
    pin = reg["commit"] or args.ref
    idx = load_index(pin)
    wf = next((w for w in idx["workflows"] if w["id"] == args.id), None)
    if wf is None:
        sys.exit(f"No workflow with id '{args.id}' in the index. "
                 f"Available: {[w['id'] for w in idx['workflows']]}")
    dest = os.path.abspath(args.dest)
    os.makedirs(dest, exist_ok=True)

    # params_file is optional: when absent, params are estimated from the data
    # type (estimate_params.py) rather than fetched.
    want = [*([wf["search"]["params_file"]] if wf["search"].get("params_file") else []),
            *wf["search"].get("extra_files", [])]
    written = []
    for rel in want:
        url = raw_url(pin, f"{wf['path']}/{rel}")
        local = os.path.join(dest, os.path.basename(rel))
        try:
            with open(local, "wb") as fh:
                fh.write(_get(url, binary=True))
            written.append(local)
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            sys.exit(f"Failed to fetch {url}: {e}")

    summary = _summary(wf, reg)
    manifest = os.path.join(dest, "workflow.manifest.json")
    with open(manifest, "w") as fh:
        json.dump(summary, fh, indent=2)

    pf = wf["search"].get("params_file")
    print(json.dumps({
        "id": wf["id"], "dir": dest, "registry": reg,
        "params_file": os.path.join(dest, os.path.basename(pf)) if pf else None,
        "estimate_params": wf["search"].get("estimate_params", not pf),
        "var_mods": wf["search"].get("var_mods", ""),
        "param_overrides": wf["search"].get("param_overrides", {}),
        "files": written, "manifest": manifest,
        "engine": wf["engine"], "fasta": wf.get("fasta", {}), "de": wf["de"],
    }, indent=2))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    m = sub.add_parser("match", help="pick the validated workflow for a run")
    m.add_argument("--acquisition", required=True, choices=["DIA", "DDA", "dia", "dda"])
    m.add_argument("--organism-taxid", required=True)
    m.add_argument("--instrument", default="")
    m.add_argument("--ref", default="main", help="registry ref to read (branch/tag/sha)")
    m.set_defaults(func=do_match)

    p = sub.add_parser("pull", help="download a chosen workflow's files, pinned to a commit")
    p.add_argument("--id", required=True)
    p.add_argument("--dest", default="./workflow_bundle")
    p.add_argument("--ref", default="main",
                   help="registry ref; pass the recorded commit SHA to reproduce an old run exactly")
    p.set_defaults(func=do_pull)

    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
