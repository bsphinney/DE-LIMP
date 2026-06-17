#!/usr/bin/env python3
"""
build_index.py  --  Compile every workflows/<id>/workflow.yaml into index.json,
the single JSON file the skill fetches (so the skill needs no YAML parser and
one HTTP request). Run after editing a workflow, or let the GitHub Action do it.

    python3 build_index.py            # writes ./index.json
    python3 build_index.py --check    # verify index.json is up to date (CI)
"""
import os, sys, json, glob, datetime

try:
    import yaml
except ImportError:
    sys.exit("PyYAML required: pip install pyyaml")

HERE = os.path.dirname(os.path.abspath(__file__))

def load(wf_dir):
    with open(os.path.join(wf_dir, "workflow.yaml")) as fh:
        m = yaml.safe_load(fh)
    folder = os.path.basename(wf_dir)
    if m.get("id") != folder:
        sys.exit(f"id '{m.get('id')}' != folder name '{folder}' in {wf_dir}")
    # denormalize exactly the fields the skill needs to match + run
    return {
        "id": m["id"],
        "name": m["name"],
        "path": f"workflows/{folder}",
        "match": {
            "acquisition": m["match"]["acquisition"],
            "instruments": m["match"].get("instruments", []),
            "instrument_aliases": m["match"].get("instrument_aliases", []),
            "organism": m["match"].get("organism"),
            "organism_taxid": m["match"]["organism_taxid"],
        },
        "engine": {"name": m["engine"]["name"], "version": m["engine"]["version"]},
        "fasta": m.get("fasta", {}),
        "search": {"estimate_params": m["search"].get("estimate_params", "params_file" not in m["search"]),
                   "var_mods": m["search"].get("var_mods", ""),
                   "param_overrides": m["search"].get("param_overrides", {}),
                   "params_file": m["search"].get("params_file"),
                   "extra_files": m["search"].get("extra_files", [])},
        "de": m["de"],
        "validated": m.get("validated", {}),
    }

def build():
    dirs = sorted(d for d in glob.glob(os.path.join(HERE, "*"))
                  if os.path.isfile(os.path.join(d, "workflow.yaml")))
    return {
        "schema_version": 1,
        "repo": "bsphinney/DE-LIMP",
        "generated": datetime.date.today().isoformat(),
        "workflows": [load(d) for d in dirs],
    }

def main():
    idx = build()
    out = os.path.join(HERE, "index.json")
    new = json.dumps(idx, indent=2, ensure_ascii=False) + "\n"
    if "--check" in sys.argv:
        cur = open(out).read() if os.path.exists(out) else ""
        # ignore the generated date when comparing
        import re
        norm = lambda s: re.sub(r'"generated": "[^"]*"', '', s)
        if norm(cur) != norm(new):
            sys.exit("index.json is stale — run: python3 build_index.py")
        print("index.json is up to date.")
        return
    with open(out, "w") as fh:
        fh.write(new)
    print(f"Wrote {out} with {len(idx['workflows'])} workflow(s).")

if __name__ == "__main__":
    main()
