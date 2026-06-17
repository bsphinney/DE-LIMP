# v4.0.1 — de novo→homology FDR fix (findings #1 + #2)

Branch: `fix/denovo-homology-fdr-v4.0.1`. **Gated on HIVE (down until Fri 2026-06-19).**
Nothing merges to `main` until the ocelot validation confirms the numbers.

Authoritative method: `~/Downloads/denovo_decoy_method.html` (NovoBoard decoy-spectra,
FRAC-0.8, ocelot). Decision: **EXTRA DIAMOND flags** (user, 2026-06-17).

## Root cause (confirmed by reading the code)
The Target-Decoy FDR compares two non-comparable arms:
- **Real arm** = `blast_results.tsv` → `values$dda_casanovo_blast`, produced by an
  app job vs **SwissProt** `--sensitive --id 50 --max-target-seqs 5`
  (`server_dda.R:1722/1745` and `2655/2676`).
- **Decoy arm** = `blast_results_decoy_spectra.tsv` → `values$dda_decoy_blast`, vs
  **nr** `--evalue 1 --max-target-seqs 25` (`helpers_dda.R:1951`,
  `generate_denovo_decoy_fdr_sbatch`).

There is **no app-generated real nr search** — the nr/LCA table is produced
externally on HIVE. So the in-app FDR is SwissProt-real vs nr-decoy = finding #2.

Finding #1: call site `server_dda.R:4359` passes `decoy_n = real_n`, forcing the
population scale (`real_n/decoy_n`, `helpers_dda.R:897-899`) to 1. The doc's ocelot
numbers show the decoy peptide population is NOT 1:1 with real, so scale≠1; the
direction is unknowable without the true decoy count → must use the real count.

New since v3.11.66: doc says **calibrate on E-value, never raw bitscore** (bitscore
length-confounded: 10,404 vs 13,074 peptides @1% FDR). App uses raw bitscore.

## Fix design (execute Friday on HIVE)
1. **One job, both arms, identical settings.** Extend
   `generate_denovo_decoy_fdr_sbatch` so it also BLASTs the **real** de novo peptides
   (extracted from the main Casanovo `*_sequence.mztab`, ≥7aa, deduped, I/L via
   `build_dda_canonical_peptide`) against the **same nr DB** with the **same EXTRA
   flags**, writing `blast_results_real_nr.tsv`. Guarantees apples-to-apples.
   - **TODO Friday:** find where the main Casanovo job writes real `*_sequence.mztab`
     (need that path as a new param `real_casanovo_dir`).
   - EXTRA §11 flags (define once in the script):
     `-b8.0 --matrix BLOSUM62 --ultra-sensitive -s2 --id2 1
      --short-query-ungapped-bitscore 1 --algo 0 --masking 0 --gapped-filter-evalue 0
      --min-score 1 --gapopen 6 --gapextend 2 --ignore-warnings --taxonlist 2759
      --outfmt 6 qseqid sseqid pident length evalue bitscore staxids
      --max-target-seqs 1`  (nr-Eukaryota = `--taxonlist 2759`).
   - Print `real peptides >=7aa: N` and `decoy peptides >=7aa: N`.
2. **Capture counts.** Loader (`server_dda.R:~1664`) reads `blast_results_real_nr.tsv`
   → `values$dda_real_nr_blast`, and the two `>=7aa` counts → `values$dda_real_nr_n`,
   `values$dda_decoy_n`.
3. **Feed FDR from the matched arms.** FDR reactive (`server_dda.R:~4357`):
   `build_denovo_*_fdr(real = values$dda_real_nr_blast, decoy = values$dda_decoy_blast,
    real_n = values$dda_real_nr_n, decoy_n = values$dda_decoy_n)`. Never `decoy_n=real_n`.
4. **E-value axis.** Add E-value calibration to the helper (rank ascending by E-value;
   threshold direction flips vs bitscore). Make it the primary; keep bitscore as a
   secondary view. Single definition (`helpers_dda.R`).
5. **Retire the legacy Casanovo-score FDR pane** (`server_dda.R:~4258-4334`,
   `build_denovo_score_calibration`) — doc says score-based gating is wrong. Relabel
   as diagnostic-only or remove.
6. **Fix the stale comment** `helpers_dda.R:1870` ("SAME pipeline as the real arm" —
   only true after this fix).
7. **Validate on ocelot (HIVE):** decoy-spectra null, EXTRA, nr-Eukaryota; confirm the
   recovery@1%FDR is in the doc's ballpark; confirm decoy-spectra ≈ reversed-DB null
   on the clean f08 decoy. Then bump VERSION 4.0.1 + CHANGELOG, merge, push.

## Deferred (separate future feature, not v4.0.1)
Full linear **mokapot** control (E-value+length+bitscore+Casanovo confidence) — the
doc's fully-validated method. Bigger build (Python mokapot dep). Doc shows it raises
*yield* at fixed FDR, not ground-truth precision, so the E-value + same-DB + true
decoy_n fix is the correct, defensible control to ship first.
