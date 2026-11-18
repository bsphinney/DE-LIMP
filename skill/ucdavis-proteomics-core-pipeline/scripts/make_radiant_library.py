#!/usr/bin/env python3
"""
make_radiant_library.py -- build a spectral library Radiant DIA can actually read.

WHY THIS EXISTS (all three facts verified on HIVE, 2026-08-03)
-------------------------------------------------------------
Radiant needs a library, and DIA-NN is the natural way to predict one. But DIA-NN 2.x
and this Radiant build share NO library format directly:

1. `--out-lib foo.tsv --predictor` does NOT write foo.tsv. DIA-NN ignores the
   extension for a PREDICTED library and always writes `foo.predicted.speclib`
   (its compact binary format).
2. Radiant DOES accept `.speclib` -- but only format version >= -3
   (FileReadersLib/src/SpecLibSrc/Library.h: `if (version < -3) ... not supported`).
   DIA-NN 2.6 writes v-11 and 2.x generally writes v-10, so every DIA-NN 2.x speclib
   is rejected with `ERROR: version is not supported11`. Every speclib on the UCD
   share is v-10 or v-11.
3. DIA-NN cannot emit TSV at all. Asked to convert a library to `foo.tsv`, it
   silently writes `foo.parquet` instead.

So the only working path is: predict -> convert to parquet -> rename columns to the
schema Radiant's FragLibTsvReader keys on -> TSV. DIA-NN's parquet uses DOT-separated
names (`Precursor.Mz`, `Relative.Intensity`); Radiant wants the OpenSWATH-style
concatenated names (`PrecursorMz`, `LibraryIntensity`). A straight dump is unreadable.

Usage:
  python3 make_radiant_library.py --diann '<diann cmd>' --fasta db.fasta --out-dir ./lib
  python3 make_radiant_library.py --from-speclib existing.predicted.speclib \
      --diann '<diann cmd>' --out-dir ./lib
  python3 make_radiant_library.py --from-parquet lib.parquet --out-dir ./lib   # no DIA-NN
"""
import sys, os, json, shlex, argparse, subprocess, csv, hashlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diann_parallel import diann_output_path  # noqa: E402  (how DIA-NN names its outputs)
# what each alkylating agent leaves on cysteine, and DIA-NN's flags for it: the one map
from estimate_params import CYS_ALKYLATION_CHOICES, cys_diann_flags  # noqa: E402
# the digest, charge and fragment ranges the DIA-NN route searches: the one definition
from estimate_params import diann_library_flags  # noqa: E402
from safe_write import write_file  # noqa: E402

# DIA-NN parquet column -> the name Radiant's FragLibTsvReader looks for.
# The reader accepts either spelling of the aliased ones (ProductMz|FragmentMz,
# LibraryIntensity|RelativeIntensity, ModifiedPeptide|ModifiedPeptideSequence,
# FragmentSeriesNumber|FragmentNumber, Tr_recalibrated|NormalizedRetentionTime);
# we emit the DIA-NN-style spelling, which matches its own test fixture.
COLMAP = {
    "Precursor.Mz":           "PrecursorMz",
    "Product.Mz":             "ProductMz",
    "Relative.Intensity":     "LibraryIntensity",
    "Modified.Sequence":      "ModifiedPeptide",
    "Stripped.Sequence":      "PeptideSequence",
    "Precursor.Charge":       "PrecursorCharge",
    "Fragment.Type":          "FragmentType",
    "Fragment.Charge":        "FragmentCharge",
    "Fragment.Series.Number": "FragmentSeriesNumber",
    "Fragment.Loss.Type":     "FragmentLossType",
    "RT":                     "Tr_recalibrated",
    "IM":                     "IonMobility",
    "Decoy":                  "decoy",
    "Protein.Group":          "ProteinGroup",
    "Protein.Names":          "ProteinName",
    "Genes":                  "Genes",
    "Precursor.Id":           "transition_group_id",
    "Q.Value":                "QValue",
    "PG.Q.Value":             "PGQValue",
    "Proteotypic":            "Proteotypic",
    "Exclude.From.Quant":     "ExcludeFromAssay",
}
# The library's name in --out-dir: the ONE definition (run_search.radiant_library_lines checks a
# job wrote it, and the parallel route's step 2 reads it).
LIBRARY_TSV = "radiant_library.tsv"
# Without these Radiant cannot build a fragment library at all.
REQUIRED = ["PrecursorMz", "ProductMz", "LibraryIntensity", "ModifiedPeptide",
            "PrecursorCharge", "FragmentType", "FragmentCharge", "FragmentSeriesNumber"]


# What a library predicted here was built from, recorded beside it (<tsv>.meta.json) and compared
# before one is reused. Up to 2.15.0 only the cysteine choice was, so a re-search into the same
# --out after the FASTA changed searched the OLD database (audit 2026-11-17, lane 1 #1).
META_SCHEMA = 2
PRECURSOR_MZ_NOTE = ("DIA-NN's own default precursor m/z range: the Radiant route does not read "
                     "the acquired isolation windows, so the library is not narrowed to them")


def sha256_file(path, block=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def library_recipe(diann, fasta, cys):
    """Everything a predicted library depends on, as recorded and compared: the FASTA's bytes,
    the DIA-NN command (its path carries the pinned version), and the prediction flags."""
    return {"schema": META_SCHEMA, "cys_alkylation": cys,
            "fasta": os.path.abspath(fasta), "fasta_sha256": sha256_file(fasta),
            "diann": diann,
            "library_flags": diann_library_flags() + cys_diann_flags(cys),
            "precursor_mz": PRECURSOR_MZ_NOTE}


def reuse_mismatch(meta, recipe):
    """None when the recorded library matches `recipe`; else what differs (said, then the library
    is rebuilt). A record from before 2.15.1 has no FASTA hash and is always rebuilt."""
    if not isinstance(meta, dict):
        return "no record of what it was built from"
    diff = [k for k in ("schema", "cys_alkylation", "fasta_sha256", "diann", "library_flags")
            if meta.get(k) != recipe.get(k)]
    return ("built with a different " + ", ".join(diff)) if diff else None


def sh(cmd):
    print(f"  $ {cmd}", flush=True)
    subprocess.run(cmd, shell=True, check=True)


def find_written(d, stem):
    """DIA-NN renames the output; find what it really wrote."""
    for ext in (".predicted.speclib", ".speclib", ".parquet", ".tsv"):
        p = os.path.join(d, stem + ext)
        if os.path.exists(p) and os.path.getsize(p) > 1000:
            return p
    return None


def parquet_to_tsv(src, dst, batch=200_000):
    try:
        import pyarrow.parquet as pq
    except ImportError:
        sys.exit("pyarrow required. pip install pyarrow")
    pf = pq.ParquetFile(src)
    cols = [c for c in pf.schema_arrow.names if c in COLMAP]
    out_names = [COLMAP[c] for c in cols]
    missing = [r for r in REQUIRED if r not in out_names]
    if missing:
        sys.exit(f"{src} lacks columns Radiant requires: {missing}\n"
                 f"  parquet had: {pf.schema_arrow.names}")
    n = 0
    with open(dst, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(out_names)
        for b in pf.iter_batches(batch_size=batch, columns=cols):
            d = b.to_pydict()
            w.writerows(zip(*(d[c] for c in cols)))
            n += b.num_rows
    return n, out_names


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--diann", help="DIA-NN command (needed unless --from-parquet)")
    ap.add_argument("--fasta", help="FASTA to predict from")
    ap.add_argument("--from-speclib", help="skip prediction; convert this .speclib")
    ap.add_argument("--from-parquet", help="skip DIA-NN entirely; convert this .parquet")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--cys-alkylation", choices=CYS_ALKYLATION_CHOICES, default="cam",
                    help="the fixed cysteine modification the predicted library carries "
                         "(estimate_params.py's choices; default cam = carbamidomethyl). Radiant "
                         "takes its modifications from the library")
    a = ap.parse_args()

    os.makedirs(a.out_dir, exist_ok=True)
    tsv = os.path.join(a.out_dir, LIBRARY_TSV)
    # What the library was predicted from, beside it: a library is reused only when the FASTA
    # (sha256), the DIA-NN command, the prediction flags and the cysteine choice all match
    # (library_recipe). Anything else -- an older record included -- is rebuilt, and said.
    meta = tsv + ".meta.json"
    predict = not (a.from_speclib or a.from_parquet)
    recipe = None
    if predict:
        if not a.fasta:
            sys.exit("--fasta is required to predict a library.")
        if not a.diann:
            sys.exit("--diann is required unless --from-parquet is given.")
        if not os.path.isfile(a.fasta):
            sys.exit(f"--fasta {a.fasta} is not a readable file.")
        recipe = library_recipe(a.diann, a.fasta, a.cys_alkylation)
    if predict and os.path.exists(tsv) and os.path.getsize(tsv) > 1_000_000:
        try:
            with open(meta) as fh:
                built = json.load(fh)
        except (OSError, ValueError):
            built = None
        why = reuse_mismatch(built, recipe)
        if why is None:
            print(json.dumps({"library": tsv, "reused": True, **recipe}, indent=2))
            return
        sys.stderr.write(f"[make_radiant_library] the library already in {a.out_dir} is not "
                         f"reused ({why}); predicting a new one\n")

    parquet = a.from_parquet
    if not parquet:
        if not a.diann:
            sys.exit("--diann is required unless --from-parquet is given.")
        speclib = a.from_speclib
        if not speclib:
            if not a.fasta:
                sys.exit("--fasta is required to predict a library.")
            # Both --out-lib names carry an extension: DIA-NN cuts the name it writes at the
            # last '.' of the whole path, so a bare stem in an --out-dir with a '.' in its name
            # put the library one folder up (diann_parallel.diann_output_path).
            out_lib = diann_output_path(os.path.join(a.out_dir, "diann_predicted_lib.speclib"),
                                        "--out-lib")
            print("[1/3] DIA-NN: predicting the library from the FASTA")
            # the same digest, charges and fragment range as the DIA-NN route's cfg, and the
            # cysteine choice (recipe["library_flags"]); DIA-NN's defaults otherwise
            flags = " ".join(shlex.quote(x) for x in recipe["library_flags"])
            sh(f"{a.diann} --fasta {shlex.quote(a.fasta)} --fasta-search --predictor "
               f"--gen-spec-lib --out-lib {shlex.quote(out_lib)} --threads {a.threads}"
               + (f" {flags}" if flags else ""))
            speclib = find_written(a.out_dir, "diann_predicted_lib")
            if not speclib:
                sys.exit(f"DIA-NN wrote no library into {a.out_dir}: "
                         f"{sorted(os.listdir(a.out_dir))}")
        # Radiant rejects DIA-NN 2.x speclib (v-10/-11 < the v>=-3 it supports), so we
        # cannot hand the binary straight over -- convert it to a readable table.
        print("[2/3] DIA-NN: speclib -> parquet (it cannot write TSV; it writes parquet)")
        out_lib2 = diann_output_path(os.path.join(a.out_dir, "radiant_library_src.parquet"),
                                     "--out-lib")
        sh(f"{a.diann} --lib {shlex.quote(speclib)} --out-lib {shlex.quote(out_lib2)} "
           f"--gen-spec-lib --threads {a.threads}")
        parquet = find_written(a.out_dir, "radiant_library_src")
        if not parquet or not parquet.endswith(".parquet"):
            sys.exit(f"expected a .parquet from the conversion; got {parquet}")

    print("[3/3] parquet -> Radiant TSV (renaming to the schema its reader keys on)")
    n, names = parquet_to_tsv(parquet, tsv)
    # predicted here: what it was predicted from; converted from a given library: unknown, so
    # never reused (no fasta_sha256)
    write_file(meta, json.dumps(recipe if predict else {"schema": META_SCHEMA,
                                                         "cys_alkylation": None},
                                indent=2) + "\n")
    print(json.dumps({
        "library": tsv,
        "rows": n,
        "size_bytes": os.path.getsize(tsv),
        "columns": names,
        "source_parquet": parquet,
        "cys_alkylation": a.cys_alkylation if predict else None,
        **({k: recipe[k] for k in ("fasta_sha256", "library_flags", "precursor_mz")}
           if predict else {}),
        "note": "Radiant reads .fragLibFF/.tsv/.csv/.speclib, but rejects DIA-NN 2.x "
                ".speclib (format v-10/-11; it supports v>=-3) — hence the TSV.",
    }, indent=2))


if __name__ == "__main__":
    main()
