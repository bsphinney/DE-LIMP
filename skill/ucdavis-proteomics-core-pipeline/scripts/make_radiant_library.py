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
import sys, os, json, shlex, argparse, subprocess, csv

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from diann_parallel import diann_output_path  # noqa: E402  (how DIA-NN names its outputs)
# what each alkylating agent leaves on cysteine, and DIA-NN's flags for it: the one map
from estimate_params import CYS_ALKYLATION_CHOICES, cys_diann_flags  # noqa: E402
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
# Without these Radiant cannot build a fragment library at all.
REQUIRED = ["PrecursorMz", "ProductMz", "LibraryIntensity", "ModifiedPeptide",
            "PrecursorCharge", "FragmentType", "FragmentCharge", "FragmentSeriesNumber"]


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
    tsv = os.path.join(a.out_dir, "radiant_library.tsv")
    # What the library was predicted with, beside it: a library is reused only for the same
    # cysteine choice. One from before 2.12 has no record -- and fixed nothing on C, whatever the
    # prep -- so it is rebuilt.
    meta = tsv + ".meta.json"
    try:
        with open(meta) as fh:
            built = json.load(fh).get("cys_alkylation")
    except (OSError, ValueError, AttributeError):
        built = None
    if (os.path.exists(tsv) and os.path.getsize(tsv) > 1_000_000 and built == a.cys_alkylation
            and not (a.from_speclib or a.from_parquet)):
        print(json.dumps({"library": tsv, "reused": True, "cys_alkylation": built}, indent=2))
        return

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
            cys = " ".join(shlex.quote(x) for x in cys_diann_flags(a.cys_alkylation))
            sh(f"{a.diann} --fasta {shlex.quote(a.fasta)} --fasta-search --predictor "
               f"--gen-spec-lib --out-lib {shlex.quote(out_lib)} --threads {a.threads}"
               + (f" {cys}" if cys else ""))
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
    # predicted here: the choice it was predicted with; converted from a given library: unknown
    write_file(meta, json.dumps({"cys_alkylation": (None if (a.from_speclib or a.from_parquet)
                                                    else a.cys_alkylation)}) + "\n")
    print(json.dumps({
        "library": tsv,
        "rows": n,
        "size_bytes": os.path.getsize(tsv),
        "columns": names,
        "source_parquet": parquet,
        "cys_alkylation": None if (a.from_speclib or a.from_parquet) else a.cys_alkylation,
        "note": "Radiant reads .fragLibFF/.tsv/.csv/.speclib, but rejects DIA-NN 2.x "
                ".speclib (format v-10/-11; it supports v>=-3) — hence the TSV.",
    }, indent=2))


if __name__ == "__main__":
    main()
