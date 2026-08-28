#!/usr/bin/env python3
"""
denovo_search_cmd.py -- emit the DIAMOND command for de novo -> homology species ID.

Two settings here are not defaults and both were measured. Getting either wrong
silently degrades the result rather than failing.

1. THE SHORT-PEPTIDE ("EXTRA") CONFIGURATION.
   De novo peptides are 7-30 aa. DIAMOND's defaults find essentially nothing in
   that range -- its heuristics are tuned for full-length proteins. The EXTRA
   settings below disable the filters that discard short alignments. Measured
   effect at 1% FDR on nr-Eukaryota: standard finds 9,801 peptides with confirmed
   homology, EXTRA finds 19,671, and the gain is concentrated exactly where it
   should be -- 168 -> 1,661 peptides at <=15 aa.
   Cost: ~5x standard (~217 core-hours vs ~44 for a 9-file cohort).

2. --max-target-seqs 25, NOT 1.
   The species call is made by LCA over ALL hits for a peptide. With one hit per
   peptide there is nothing to take an ancestor of, so LCA degrades silently to
   best-hit -- the 36%-mis-assignment error mode it exists to prevent. The ocelot
   production search used --max-target-seqs 1 (0.9 hits/peptide) and could not
   run the LCA step its own method section recommends. The HeLa benchmark used 25
   (22.5 hits/peptide).
   Cost: ~22x the hit rows. This is a real trade against (1) on a large cohort,
   which is why it is stated rather than hidden.

The binary matters too: use the Riffle build (diamond_mriffle_2.1.10.sif on
HIVE). Standard DIAMOND releases do not accept --short-query-ungapped-bitscore.
"""
import argparse
import json
import shlex

EXTRA = ["-b8.0", "--matrix", "BLOSUM62", "--ultra-sensitive", "-s2", "--id2", "1",
         "--short-query-ungapped-bitscore", "1", "--algo", "0", "--masking", "0",
         "--gapped-filter-evalue", "0", "--min-score", "1",
         "--gapopen", "6", "--gapextend", "2"]
STANDARD = ["--evalue", "1"]
OUTFMT = ["--outfmt", "6", "qseqid", "sseqid", "pident", "length", "mismatch",
          "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
          "staxids"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query", required=True)
    ap.add_argument("--db", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--mode", choices=["extra", "standard"], default="extra",
                    help="extra = short-peptide sensitive (default, ~5x cost)")
    ap.add_argument("--max-target-seqs", type=int, default=25,
                    help="25 keeps LCA possible; 1 silently reduces it to best-hit")
    ap.add_argument("--taxonlist", default="2759",
                    help="2759 = Eukaryota; empty string for all of nr")
    ap.add_argument("--threads", type=int, default=32)
    ap.add_argument("--sif", default="/quobyte/proteomics-grp/apptainers/diamond_mriffle_2.1.10.sif",
                    help="Riffle build; standard DIAMOND lacks --short-query-ungapped-bitscore")
    a = ap.parse_args()

    warnings = []
    if a.max_target_seqs < 2:
        warnings.append(
            f"--max-target-seqs {a.max_target_seqs} makes LCA impossible: with one hit "
            f"per peptide the species step degrades to best-hit, which mis-assigns "
            f"~36% of peptides. Use 25 unless you are deliberately skipping LCA.")
    if a.mode == "standard":
        warnings.append(
            "standard mode finds about half the homology of EXTRA on short peptides "
            "(9,801 vs 19,671 at 1% FDR; 168 vs 1,661 peptides at <=15 aa). Use it "
            "only when compute-bound.")

    cmd = ["diamond", "blastp"]
    cmd += EXTRA if a.mode == "extra" else STANDARD
    if a.taxonlist:
        cmd += ["--taxonlist", a.taxonlist]
    cmd += ["--ignore-warnings"] + OUTFMT
    cmd += ["--max-target-seqs", str(a.max_target_seqs),
            "--threads", str(a.threads),
            "--query", a.query, "--db", a.db, "--out", a.out]
    full = ["apptainer", "exec", "--bind", "/quobyte:/quobyte", a.sif] + cmd if a.sif else cmd

    print(json.dumps({
        "mode": a.mode,
        "max_target_seqs": a.max_target_seqs,
        "lca_possible": a.max_target_seqs >= 2,
        "command": " ".join(shlex.quote(c) for c in full),
        "warnings": warnings,
    }, indent=2))


if __name__ == "__main__":
    main()
