#!/usr/bin/env python3
"""
resolve_defaults.py -- pick the search defaults for a run FROM THE DATA TYPE.

Replaces the old fetch_workflows.py registry. Same output contract
(workflow.manifest.json + a JSON summary on stdout), so run_search.py and
provenance.py consume it unchanged -- but nothing is fetched over the network
and there is no menu.

WHY THE REGISTRY WENT AWAY
--------------------------
It fetched per-organism "validated workflow" bundles from GitHub at run time.
Three things were wrong with that:

  * Search parameters do not depend on species. Only the FASTA does, and that
    always came from the user's own organism answer -- never from the bundle.
    So the organism dimension bought nothing and forced a near-duplicate bundle
    per species.
  * A remote, mutable registry meant one skill version produced different
    parameters on different days. In practice two bundles pinned two different
    DIA-NN versions (2.3.0 and 2.6.0) and a run silently got whichever one its
    instrument matched.
  * Every bundle was unvalidated or placeholder-validated, so the "is this
    validated?" gate fired on essentially every run and turned the first
    interaction into a menu of warnings.

Defaults now ship WITH the skill, so the skill version fully determines the
search. Site-specific overrides are explicit flags, not a remote lookup.

  python3 resolve_defaults.py --acquisition DIA --instrument "timsTOF HT" \
      --engine diann --dest ./wf [--env env.json] [--fasta db.fasta] [--threads 32]
"""
import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from estimate_params import classify_instrument  # noqa: E402  (one ppm table)

HERE = os.path.dirname(os.path.abspath(__file__))

# Version the DEFAULTS TABLE itself, so a run record can name what produced its
# parameters without depending on a git commit in another repo.
DEFAULTS_VERSION = "2026-08-14"

# Pinned engine versions. One place. Bump here, not in six files.
ENGINE_VERSIONS = {
    "diann": "2.6.1",
    "sage": "0.14.7",
    "fragpipe": "24.0",
    "radiant": "2.3.3",   # container release, NOT the Radiant DIA version inside
}

# Which engines are sensible for a data type, best first. Instrument class comes
# from estimate_params.classify_instrument, so the ppm table and this table can
# never disagree about what an instrument is.
ROUTES = {
    ("DIA", "timstof"): ["diann", "fragpipe"],
    ("DIA", "orbitrap"): ["diann", "fragpipe", "radiant"],
    ("DIA", "other"): ["diann"],
    ("DDA", "timstof"): ["sage"],
    ("DDA", "orbitrap"): ["sage"],
    ("DDA", "other"): ["sage"],
}

DE_METHOD = {"diann": "dpc", "fragpipe": "dpc", "radiant": "dpc", "sage": "maxlfq"}

# Engines configured by a whole config file rather than by flags. These are
# generated per run by make_presets.py.
PRESET_ENGINES = {"fragpipe": ".workflow", "radiant": ".radiantConfig"}


def family(cls):
    if cls == "timstof":
        return "timstof"
    if cls.startswith("orbitrap"):
        return "orbitrap"
    return "other"


def load_env(path):
    if not path:
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError) as e:
        sys.stderr.write(f"[resolve_defaults] could not read --env {path}: {e}\n")
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--acquisition", required=True, choices=["DIA", "DDA", "dia", "dda"])
    ap.add_argument("--instrument", default="")
    ap.add_argument("--engine", default="", help="force an engine; else the data-type default")
    ap.add_argument("--organism-taxid", default=None,
                    help="recorded for provenance only -- it does NOT affect any search parameter")
    ap.add_argument("--ms1-res", type=float, default=None)
    ap.add_argument("--ms2-res", type=float, default=None)
    ap.add_argument("--ms1-ppm", type=float, default=None, help="site SOP override")
    ap.add_argument("--ms2-ppm", type=float, default=None, help="site SOP override")
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--env", default=None, help="detect_env.sh JSON, to flag engines that "
                                                "cannot run on this machine")
    ap.add_argument("--dest", default="./wf")
    args = ap.parse_args()

    acq = args.acquisition.upper()
    cls, ms1, ms2, label, src = classify_instrument(args.instrument, args.ms1_res, args.ms2_res)
    fam = family(cls)
    routes = ROUTES[(acq, fam)]

    engine = (args.engine or "").strip().lower()
    if engine and engine not in ENGINE_VERSIONS:
        sys.exit(f"resolve_defaults: unknown engine '{engine}'. "
                 f"Known: {', '.join(sorted(ENGINE_VERSIONS))}")
    notes = []
    if engine and engine not in routes:
        # Honour it, but never silently: the user asked for something this data
        # type does not normally use.
        notes.append(f"{engine} is not a default route for {acq} on {label} "
                     f"(defaults: {', '.join(routes)}) -- running it because it was requested")
    if not engine:
        engine = routes[0]

    # Platform gate. Verified facts live in detect_env.sh; we only read them.
    env = load_env(args.env)
    availability = None
    if env and isinstance(env.get("engines"), dict):
        info = env["engines"].get(engine) or {}
        availability = info
        if info and not info.get("available", True):
            alt = [e for e in routes
                   if (env["engines"].get(e) or {}).get("available")]
            sys.exit(f"resolve_defaults: {engine} cannot run on this machine -- "
                     f"{info.get('note', 'unavailable')}"
                     + (f"\n  Available here for this data type: {', '.join(alt)}"
                        if alt else "\n  No engine for this data type runs here; use HIVE."))
        if info.get("note"):
            notes.append(f"{engine}: {info['note']}")

    dest = os.path.abspath(args.dest)
    os.makedirs(dest, exist_ok=True)

    # Engines driven by a config file get one generated now, before the search.
    params_file = None
    preset_prov = None
    if engine in PRESET_ENGINES:
        params_file = os.path.join(dest, f"{engine}{PRESET_ENGINES[engine]}")
        cmd = [sys.executable, os.path.join(HERE, "make_presets.py"),
               "--engine", engine, "--instrument", args.instrument,
               "--acquisition", acq, "--out", params_file]
        for flag, val in (("--ms1-res", args.ms1_res), ("--ms2-res", args.ms2_res),
                          ("--ms1-ppm", args.ms1_ppm), ("--ms2-ppm", args.ms2_ppm),
                          ("--fasta", args.fasta), ("--threads", args.threads)):
            if val is not None:
                cmd += [flag, str(val)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"resolve_defaults: preset generation failed:\n{r.stderr.strip()}")
        preset_prov = json.loads(r.stdout)

    # Same manifest shape the old registry emitted, so downstream steps are
    # unchanged. Fields that only made sense for a remote bundle are explicit
    # about being local now rather than being filled with plausible-looking junk.
    manifest = {
        "id": f"{engine}_{fam}_{acq.lower()}",
        "name": f"{engine} — {label} ({acq})",
        "path": None,
        "acquisition": acq,
        "instruments": [args.instrument] if args.instrument else [],
        "instrument_class": cls,
        "instrument_label": label,
        "organism": None,
        "organism_taxid": int(args.organism_taxid) if args.organism_taxid else None,
        "engine": {"name": engine, "version": ENGINE_VERSIONS[engine]},
        "fasta": {"path": os.path.abspath(args.fasta) if args.fasta else None,
                  "add_contaminants": True},
        "de": {"method": DE_METHOD[engine], "q_cutoff": 0.01, "logfc": 1.0, "adjp": 0.05},
        "search": {
            "estimate_params": engine not in PRESET_ENGINES,
            "params_file": params_file,
            "ms1_ppm": args.ms1_ppm if args.ms1_ppm is not None else ms1,
            "ms2_ppm": args.ms2_ppm if args.ms2_ppm is not None else ms2,
            "ppm_source": ("site SOP override" if args.ms1_ppm is not None else src),
            "preset_provenance": preset_prov,
        },
        "validated": {
            "date": None, "by": None, "dataset": None,
            "notes": ("Parameters are this skill's data-type defaults "
                      f"(resolve_defaults.py {DEFAULTS_VERSION}), not a site validation run. "
                      "They are reproducible from the skill version alone."),
        },
        "registry": {"source": "local skill defaults",
                     "defaults_version": DEFAULTS_VERSION,
                     "note": "no remote registry; parameters ship with the skill"},
        "platform": availability,
        "alternatives": [e for e in routes if e != engine],
        "notes": notes,
    }
    mpath = os.path.join(dest, "workflow.manifest.json")
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=2)

    print(json.dumps({
        "id": manifest["id"], "dir": dest, "manifest": mpath,
        "engine": manifest["engine"], "de": manifest["de"],
        "params_file": params_file,
        "estimate_params": manifest["search"]["estimate_params"],
        "ms1_ppm": manifest["search"]["ms1_ppm"],
        "ms2_ppm": manifest["search"]["ms2_ppm"],
        "ppm_source": manifest["search"]["ppm_source"],
        "instrument_label": label,
        "alternatives": manifest["alternatives"],
        "notes": notes,
    }, indent=2))


if __name__ == "__main__":
    main()
