#!/usr/bin/env python3
"""
diatracer_stage.py  --  Stage dia-PASEF .d files for FragPipe/diaTracer so the run
never writes into the shared raw-data directory, and so an existing conversion is
reused instead of redone.

WHY THIS EXISTS
---------------
diaTracer writes its pseudo-MS/MS mzML NEXT TO THE INPUT FILE, not into --workdir.
FragPipe builds that output path like this (CmdDiaTracer.java:138):

    f.getPath().toAbsolutePath().normalize().toString().replace(".d", "_diatracer.mzML")

Two consequences drive this script:

  1. Two people searching the SAME shared .d files -- a class working from one
     dataset -- would both write <name>_diatracer.mzML into the same directory and
     race each other. On a read-only share the run just fails.

  2. `normalize()` is NOT `toRealPath()`: it collapses "." and ".." but does NOT
     resolve symlinks. So if the manifest points at a SYMLINK to the .d, diaTracer's
     output lands next to the SYMLINK -- inside the user's own staging directory.
     That is what this script sets up. Each user gets their own pseudo-spectra, the
     shared raw directory is never touched, and it may be read-only.

Reuse: diaTracer only needs to run once per file. The FragPipe preset says so itself
-- if the mzML already exists, load it as 'DDA' and the original .d as 'DIA-Quant' to
skip reconversion. This script detects existing mzML (in the staging dir, or beside
the real .d from an earlier run) and writes that two-row form automatically.

  fresh file   -> <stage>/<name>.d              DIA         (diaTracer will convert)
  already done -> <mzML path>                   DDA         (reuse; skips conversion)
                  <stage>/<name>.d              DIA-Quant   (DIA-NN quantifies vs .d)

Usage:
  python3 diatracer_stage.py --raw /data/*.d --stage <session>/input/diatracer_stage \\
      --manifest <session>/input/fragpipe.fp-manifest [--conditions conditions.csv]
"""
import argparse, csv, glob, json, os, sys

# diaTracer/FragPipe have used more than one spelling across versions; match them all.
MZML_SUFFIXES = ("_diatracer.mzml", ".diatracer.mzml", "_diaTracer.mzML".lower())


def expand(patterns):
    out = []
    for p in patterns:
        hits = sorted(glob.glob(p))
        out.extend(hits or [p])
    return [os.path.abspath(x.rstrip("/")) for x in out]


def find_existing_mzml(*dirs_and_base):
    """Look for an already-converted pseudo-MS/MS mzML for this .d, case-insensitively.
    Returns its path, or None. Checking beside the real .d as well means a conversion
    someone did earlier (or a facility pre-conversion) is picked up rather than redone."""
    *dirs, base = dirs_and_base
    for d in dirs:
        if not d or not os.path.isdir(d):
            continue
        try:
            entries = os.listdir(d)
        except OSError:
            continue
        for e in entries:
            el = e.lower()
            if not el.endswith(".mzml"):
                continue
            for suf in MZML_SUFFIXES:
                if el == (base + suf).lower():
                    return os.path.join(d, e)
    return None


def check_stage_path(stage):
    """FragPipe rewrites the output name with Java's String.replace, which replaces EVERY
    occurrence of '.d'. A '.d' anywhere in a parent directory name corrupts the path
    (/data/proj.d/run.d -> /data/proj_diatracer.mzML/run_diatracer.mzML), so refuse."""
    bad = [p for p in os.path.abspath(stage).split(os.sep) if ".d" in p and not p.endswith(".d")]
    # a component that IS "<name>.d" is fine (that's a raw file); one that merely CONTAINS
    # ".d" mid-name is what breaks the rewrite.
    bad = [p for p in bad if not p.endswith(".d")]
    if bad:
        sys.exit(f"Staging path contains '.d' inside a directory name {bad} — FragPipe's "
                 "path rewrite would mangle the output. Choose a path without '.d' in it.")


def load_groups(path):
    """Optional conditions.csv -> {run_basename: group} for the manifest's experiment column."""
    if not path or not os.path.exists(path):
        return {}
    out = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            name = (row.get("File.Name") or row.get("Run") or "").strip()
            grp = (row.get("Group") or "").strip()
            if name:
                out[os.path.splitext(os.path.basename(name.rstrip("/")))[0]] = grp
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--raw", nargs="+", required=True, help="dia-PASEF .d paths or globs")
    ap.add_argument("--stage", required=True, help="per-user staging dir for the symlinks")
    ap.add_argument("--manifest", required=True, help="where to write the .fp-manifest")
    ap.add_argument("--conditions", help="optional conditions.csv to fill the experiment column")
    a = ap.parse_args()

    check_stage_path(a.stage)
    os.makedirs(a.stage, exist_ok=True)
    os.makedirs(os.path.dirname(os.path.abspath(a.manifest)) or ".", exist_ok=True)

    raws = expand(a.raw)
    if not raws:
        sys.exit("No .d files matched --raw.")
    groups = load_groups(a.conditions)

    rows, fresh, reused, notes = [], [], [], []
    for real in raws:
        if not os.path.exists(real):
            sys.exit(f"No such raw file: {real}")
        name = os.path.basename(real)
        if not name.endswith(".d"):
            sys.exit(f"{name} is not a .d directory — diaTracer only runs on Bruker .d "
                     "(FragPipe hard-errors if a DIA input is not .d).")
        base = name[:-2]                                  # strip the trailing '.d'
        link = os.path.join(a.stage, name)                # symlink KEEPS the .d suffix,
        # otherwise FragPipe skips it (CmdDiaTracer requires endsWith(".d")).
        if os.path.islink(link):
            if os.path.realpath(link) != os.path.realpath(real):
                os.unlink(link)
                os.symlink(real, link)
        elif os.path.exists(link):
            sys.exit(f"{link} exists and is not a symlink — refusing to overwrite it.")
        else:
            os.symlink(real, link)

        grp = groups.get(base, "")
        mz = find_existing_mzml(a.stage, os.path.dirname(real), base)
        if mz:
            # Reuse: the converted spectra go in as DDA, and the .d as DIA-Quant so
            # DIA-NN still quantifies against the real chromatograms.
            rows.append((mz, grp, "", "DDA"))
            rows.append((link, grp, "", "DIA-Quant"))
            reused.append({"raw": real, "mzml": mz})
        else:
            rows.append((link, grp, "", "DIA"))
            fresh.append({"raw": real, "symlink": link,
                          "will_write": os.path.join(a.stage, base + "_diatracer.mzML")})

    with open(a.manifest, "w", newline="") as fh:
        for path, exp, rep, dtype in rows:
            fh.write(f"{path}\t{exp}\t{rep}\t{dtype}\n")

    raw_dirs = {os.path.dirname(r) for r in raws}
    unwritable = [d for d in raw_dirs if not os.access(d, os.W_OK)]
    if unwritable:
        notes.append("Raw directory is not writable — which is exactly why the symlink "
                     "staging is used; diaTracer will write into the staging dir instead.")
    if reused:
        notes.append(f"{len(reused)} file(s) already had pseudo-MS/MS spectra and will NOT "
                     "be reconverted (~20 min each saved).")
    if not fresh:
        notes.append("Nothing to convert — every file already has its pseudo-MS/MS spectra.")

    print(json.dumps({
        "stage_dir": os.path.abspath(a.stage),
        "manifest": os.path.abspath(a.manifest),
        "n_files": len(raws),
        "n_to_convert": len(fresh),
        "n_reused": len(reused),
        "to_convert": fresh,
        "reused": reused,
        "notes": notes,
    }, indent=2))


if __name__ == "__main__":
    main()
