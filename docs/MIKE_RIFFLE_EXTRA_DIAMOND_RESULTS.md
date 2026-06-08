# EXTRA-DIAMOND vs standard DIAMOND on de novo peptides + a randomized-decoy null

**Test set:** ocelot (*Leopardus pardalis*) hair, 9× Thermo Exploris-480 DDA. Casanovo v5.0.0 de novo →
DIAMOND blastp → NCBI nr. Question: effect of the `-DEXTRA=ON` sensitivity settings on a target/decoy
test where the decoy is **randomized spectra** (NovoBoard), not a decoy database.

## 1. Materials

| Item | Value |
|---|---|
| DIAMOND (EXTRA) | `mriffle/diamond:2.1.10` → v2.1.10.164, run via Apptainer |
| DIAMOND (standard) | diamond/2.1.7 |
| Database | NCBI nr, DIAMOND db v3, **707,028,945 seqs / 272,881,947,790 letters** |
| nr subset (EXTRA run) | `--taxonlist 2759` (Eukaryota); taxonomy filter built in 46 s |
| Real query | 312,820 Casanovo peptides (mods stripped, ≥7 aa, deduped) |
| Decoy query | 416,673 peptides from NovoBoard decoy spectra (FRAC 0.5: 50% peaks replaced from global pool, precursor m/z+charge preserved) |
| Hardware | 32-core / 256 GB nodes (chunked array for the EXTRA run) |

## 2. Commands (verbatim)

**Standard** (full nr):
```
diamond blastp -d nr -q all_casanovo_peptides.fasta \
  --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids qseq sseq btop \
  --evalue 1 --max-target-seqs 25 --threads 32 --ignore-warnings
```

**EXTRA** (Mike's flags verbatim; only resource flags + db scope changed — see §5):
```
diamond blastp -b8.0 --matrix BLOSUM62 --ultra-sensitive -s2 --id2 1 \
  --short-query-ungapped-bitscore 1 --algo 0 --masking 0 --gapped-filter-evalue 0 \
  --min-score 1 --gapopen 6 --gapextend 2 --ignore-warnings --taxonlist 2759 \
  --outfmt 6 qseqid sseqid pident length evalue bitscore staxids \
  --max-target-seqs 1 --threads 32 -d nr -q <chunk>.fasta
```

## 3. Results — summary

| condition | real hit | decoy hit | enrichment (real/decoy) | peptides @1% FDR |
|---|---|---|---|---|
| **Standard** (E ≤ 1, full nr) | 9,852 / 312,820 = **3.15%** | 356 / 416,673 = **0.085%** | **37×** | **463** |
| **EXTRA** (min-score 1, ultra-sens, Eukaryota) | 299,099 / 312,820 = **95.6%** | 413,409 / 416,673 = **99.2%** | **0.96×** | **0** |

FDR by target/decoy competition on Casanovo score; "@1% FDR" = max real peptides accepted with
cumulative FDR ≤ 0.01.

## 4. Results — hit-rate by peptide length (%)

| length (aa) | n real | standard real | EXTRA real | n decoy | standard decoy | EXTRA decoy |
|---|---|---|---|---|---|---|
| ≤8 | 31,797 | 0.0 | 95.8 | 61,982 | 0.00 | 97.94 |
| 8–10 | 44,932 | 0.0 | 97.8 | 75,677 | 0.00 | 99.02 |
| 10–12 | 43,291 | 0.0 | 98.7 | 62,358 | 0.00 | 99.37 |
| 12–15 | 53,560 | 0.3 | 99.5 | 70,977 | 0.02 | 99.66 |
| 15–20 | 64,162 | 3.0 | 99.6 | 75,910 | 0.17 | 99.72 |
| 20+ | 65,534 | 11.9 | 98.8 | 69,769 | 0.31 | 99.42 |

Standard: hits are length-gated (≤12 aa ≈ 0%, because a short exact peptide cannot reach E ≤ 1 against
7×10⁸ sequences). EXTRA: both populations saturate at every length; **decoy ≥ real in every bin**.

## 5. Adaptations to Mike's command (and why)

| Change | Reason | Affects results? |
|---|---|---|
| Docker → Apptainer (`apptainer pull docker://mriffle/diamond:2.1.10`) | HPC is Apptainer-only | no |
| `+ --ignore-warnings` | de novo poly-A/T peptides trip the DNA-letters guard (aborts at query load) | no |
| `- -c1` | single index chunk OOMs on the 273-Gletter nr at 250 GB; auto-chunk instead | no (RAM/speed only) |
| `-b50 → -b8`, `--threads 90 → 32` | fit nodes | no (RAM/speed only) |
| `+ --taxonlist 2759` (Eukaryota) | drop bacterial/viral bulk; ~2× faster | yes (db scope) |
| all sensitivity flags | **unchanged, verbatim** | — |

## 6. Interpretation

- **Mechanism of the FDR collapse:** `--min-score 1` with no E-value filter reports any alignment with
  bitscore ≥ 1. For peptide-length amino-acid strings against ~2×10⁸ proteins, a chance match scoring ≥ 1
  exists for essentially every query — real, randomized-decoy, or random. Both populations therefore
  saturate (95–99%) and the score carries no target/decoy discrimination. The E-value threshold was the
  sole source of separation in this scheme (37× → ~1×). Confirmed independent of db scope: the same
  saturation appears on a 5,000-peptide full-nr EXTRA run (96.6% real).
- **Why this is expected given the decoy design:** these decoys are *realistic*, not random — ~20–50% real
  peaks retained, noise peaks drawn from the real fragment-ion pool, real precursor mass. The decoy
  peptides are valid-mass, plausible-composition sequences, so they behave like real peptides under
  permissive search. This is a conservative (hard) null by construction.
- **Compatibility with Mike's pipeline:** Mike's FDR is a **decoy database**
  (`uniprot_bacteria…plusdecoys.dmnd` + `--max-target-seqs 1`): single best hit, FDR = fraction of
  best-hits landing on a decoy sequence. That competition is robust to permissive matching (a real
  peptide's best hit is its true target above any decoy; noise is ~50/50 target/decoy). So the EXTRA
  settings + decoy-database are internally consistent. The two approaches null different error modes:
  decoy-database = "real protein vs chance database match"; decoy-spectra = "real sequence vs de novo
  sequencing artifact."

## 7. Proposed follow-up

Run the EXTRA settings against a **target+decoy nr** with `--max-target-seqs 1` and best-hit target/decoy
competition. Expectation: retains the recall gains (§4) while restoring FDR control via the decoy
database, decoupling FDR from the E-value filter.

## 8. Feasibility notes

- Standard full-nr search of 312,820 peptides: ~32 min / 32 threads (DIAMOND 2.1.7).
- EXTRA ultra-sensitive (Eukaryota) is db-streaming-bound; parallelized by query-chunking (per-chunk
  ~40–70 min, run in parallel). `--taxonlist` filters alignments but still streams the reference (~2×, not
  proportional to subset size). DIAMOND ≈ 500× faster than NCBI blastp at comparable sensitivity
  (Buchfink et al., Nat Methods 2015), which is what makes nr-scale target/decoy feasible at all.
