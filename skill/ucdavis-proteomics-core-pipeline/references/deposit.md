# Depositing the data in a public repository (`make_deposit.py`)

`session.py finalize` (step 12) ends every session with two things a paper needs:

1. **Publication Methods** — `output/methods.md` (+ `methods.docx`): LC-MS acquisition, the
   database search (engine, the version that ran, parameters, FDR), the sequence database and
   contaminants, the differential-expression analysis, and the UC Davis instrument-grant
   acknowledgment. If the file is missing, finalize writes it with `make_methods.py` from what the
   session holds; an existing `methods.md` is never overwritten (it may be hand-polished) — if it
   lacks a section, a complete draft is written beside it as `methods_complete_draft.md`.
2. **`output/DATA_SUBMISSION/`** — everything to deposit the data in PRIDE (default) or MassIVE:

| File | What it is |
|---|---|
| `HOW_TO_SUBMIT.md` / `.html` | Step by step, written for this session (counts, paths, file table). The `.html` opens by double-click. |
| `sdrf.tsv` | SDRF-Proteomics v1.1.0 sample sheet, one row per raw file. Filled from the session; `TO-FILL` where only the user knows. |
| `protocols.txt` | "Sample processing protocol" and "Data processing protocol" texts for the submission form, from `methods.md`, with their character counts. |
| `files_to_upload.tsv` | Every file to upload: PRIDE file type, MassIVE category, requirement, source path + where it lives (HIVE /quobyte, Flinders, or local), path inside the session zip, size, md5/sha1 (session files ≤ 256 MB), what to do before upload. |
| `prepare_upload.sbatch` | **Written, never run by the skill.** Packs each `.d` run into one `<run>.d.tar.gz`, links/copies everything into one staging folder, writes `checksum.txt` (SHA-1, PRIDE format) + `md5sums.txt`. Refuses to run outside SLURM unless `RUN_HERE=1` (laptop with local raws). Resumable. |

Every part is recorded in **`MANIFEST.txt`** at the session root (the top of the zip) as
`[OK]` or `[SKIPPED] <name> -- <reason>` — the Python mirror of DE-LIMP's `safe_section()`
(`Manifest` in `make_deposit.py`). A session without DE (QC-only), without conditions, from
Sage/FragPipe, or finalized away from the raw files still finalizes: the missing parts are
`TO-FILL` cells or `[SKIPPED]` lines with the reason, never a crash or a silent gap.
`--no-deposit` skips the package (the Methods are still ensured). The zip excludes
`output/DATA_SUBMISSION/upload_staging/` (raw archives), like `output/raw_data/`.

Run it by hand: `python3 scripts/make_deposit.py --session <session_dir>`.

## What is filled, and from where

**Never guessed.** A value the session did not record is `TO-FILL` — sex, age, disease,
organism part, cell type, cell line, sample preparation. The SDRF reserved words
`not available` / `not applicable` are *claims* ("unknown", "does not apply"), so the skill
leaves them for the user to choose; `HOW_TO_SUBMIT.md` lists each `TO-FILL` column with what to
enter and why it is empty. `[... — confirm]` tags from `make_methods.py` are carried into
`protocols.txt` and counted.

| SDRF column | Source |
|---|---|
| organism | `input/search.fasta.meta.json` organism (the organism the user confirmed; strain text in parentheses removed, flagged) |
| biological replicate | numbered within each `conditions.csv` Group (the DE treated each run as a sample) |
| acquisition method | workflow manifest (step 2 detection); diaPASEF noted but written as the parent DIA term (see below) |
| label | `label free sample` only when the search parameters show no `--channels`/label mods (DIA-NN) or no TMT (Sage); otherwise `TO-FILL` |
| instrument | raw metadata / manifest, mapped to a PSI-MS term from a table verified in OLS4 (2026-09-24); unmapped → `TO-FILL` with the recorded name |
| cleavage agent | the search's in-silico rule: DIA-NN `--cut K*,R*` → Trypsin/P (MS:1001313), `K*,R*,!*P` → Trypsin (MS:1001251); Sage `cleave_at`/`restrict` likewise. Flagged "confirm it matches the wet-lab enzyme" |
| modification parameters | the search parameters; Unimod names/accessions only for ids verified in `unimod.xml` (1, 4, 7, 21, 28, 35); others `TO-FILL` |
| precursor/fragment tolerance | only when fixed in the resolved parameters; DIA-NN auto-calibration → column omitted, reason noted |
| ms1 scan range | Bruker `analysis.tdf` MzAcqRange, only when read from the raw files |
| fraction / technical replicate | `1` (one raw file per sample in the DE design) — flagged to change if fractionated / re-injected |
| factor value | `conditions.csv` Group values under the header `factor value[TO-FILL]`: the user names the variable |

The search parameters come from `make_methods.search_record()` — the one reader, so the
Methods paragraph and the SDRF cannot disagree.

## Sources (read 2026-09-24)

| Source | What it settled |
|---|---|
| PRIDE "How to submit data" — https://www.ebi.ac.uk/pride/markdownpage/submitdatapage (markdown at `/pride/markdown/submitdatapage/content.md`) | Routes: "The PRIDE Submission Tool, a desktop wizard … uploads it over Aspera or FTP. This is the standard route" and Globus "for very large datasets or when Aspera/FTP traffic is blocked". File types table (RAW required; SEARCH "Required when no RESULT files are provided"; EXPERIMENTAL DESIGN = SDRF "Strongly recommended"; FASTA recommended; SPECTRUM_LIBRARY "Recommended for DIA/library searches"). Relations: "Every RAW file must be related to at least one RESULT or SEARCH file … Every SEARCH file must be related to at least one RAW file." Registration sends no confirmation e-mail (24 h note). Submission reference `1-XXXXXXXX-X`; "Processing can take up to five working days". Private by default + reviewer account; **Publish** button after acceptance. File names "letters, digits, underscore and dot". |
| PRIDE Submission Tool page — …/markdownpage/pridesubmissiontool | The tool is still current ("also called the PX Submission Tool"), Java 21, `start.sh`; SHA-1 checksum step "strongly recommended … you may skip it"; Aspera default, FTP fallback; resumes after a crash; hyphen allowed in names. |
| PRIDE Submission Tool **2.11.6** (release 2026-09-15; `SubmissionValidator` strings in the jar) | Protocol limits: "Sample processing protocol must be both more than 50 and less than 5000 characters" (same for data processing protocol and description); "Project title must be less than 500 and more than 30 characters"; "Filenames must contains only -_.A-Za-z0-9", "should start with alpha numeric"; species/tissue/instrument/modifications required; "Submitter name must have a space". File-type enum: RESULT PEAK SEARCH RAW QUANT GEL FASTA SPECTRUM_LIBRARY MS_IMAGE_DATA OPTICAL_IMAGE OTHER EXPERIMENTAL_DESIGN. No command-line mode found. Built-in SDRF validation against the PRIDE SDRF validator API. |
| PRIDE Data submission guidelines v2.2.0 — …/markdownpage/datasubmissionguidelines | PARTIAL vs COMPLETE: PRIDE "no longer emphasize[s]" it; COMPLETE only with mzIdentML/mzTab, otherwise registered as PARTIAL "behind the scenes". DIA-NN: `report.parquet` mandatory (ANALYSIS), log mandatory (OTHER), site report mandatory if produced, library mandatory if used, matrices recommended. FragPipe: psm.tsv, protein.tsv, fragpipe.workflow, fp-manifest mandatory. SDRF: "Recommended". `.d` folders "must be compressed … one folder per compressed file"; ZIP/GZ/TAR.GZ only. FASTA: public UniProt → name + release suffice; custom databases must be provided. |
| PRIDE "Software-specific recommendations" — …/specificsoftwareformats | "DIA results are currently submitted with the software's native output as SEARCH files; there is no widely adopted standard result export yet for DIA." |
| PRIDE Globus / checksum / submission.px pages | Globus needs `submission.px` (tool: Export summary) + `checksum.txt`; checksum format "file name (without path) and the checksum separated by a tab", SHA-1. |
| PRIDE SDRF page — …/markdownpage/sdrf | SDRF uploaded as the EXPERIMENTAL DESIGN file; editor/validator URLs. |
| PRIDE citation page — …/markdownpage/citationpage | The data-availability sentence and the 2025 NAR reference, used verbatim. |
| PRIDE data policy — …/markdownpage/datapolicy (v1.0, May 2022) | Private up to 2 years, then a release date is required. |
| ProteomeXchange submission page — https://www.proteomexchange.org/submission | Still points PRIDE submitters to the PX Submission Tool download. |
| ProteomeXchange guidelines v3.0.1 (13 Oct 2019) — https://www.proteomexchange.org/docs/guidelines_px.pdf | Complete vs partial definitions; Table 3: DIA MS/MS — PRIDE "Partial only", MassIVE "Partial and complete"; PXD is the identifier to cite. |
| MassIVE documentation — https://ccms-ucsd.github.io/MassIVEDocumentation/ (repo CCMS-UCSD/MassIVEDocumentation, 2026-06-01) | Account required; FTP to `massive-ftp.ucsd.edu` "with TLS/SSL explicit encryption … port 21"; workflow "MassIVE Dataset Submission"; required metadata (species, instrument, PTMs, keywords, PI) and a non-empty dataset password; file categories; reviewer `{MSV}_reviewer`; PX announcement on **Make Public**. No SDRF category. Does not say whether `.d` folders should be zipped, nor when a PXD is shown for a private dataset. |
| SDRF-Proteomics v1.1.0 (2026-01) — https://github.com/bigbio/proteomics-sample-metadata (`sdrf-proteomics/README.adoc`, SAMPLE-GUIDELINES.adoc; commit of 2026-08-10) | Required columns; reserved words must be lower case; `.d` folders named directly in `comment[data file]`; templates listed as leaves (`human` + `dia-acquisition`); modification key order NT, AC, …; factor values after comments. |
| SDRF templates — https://github.com/bigbio/sdrf-templates (2026-08-27) | Which columns are required and whether `not available`/`not applicable` is allowed (e.g. human `age`: not applicable disallowed; `cell line`: neither allowed). |
| `sdrf-pipelines` 0.1.6 (PyPI, the released validator; run locally) | The generated SDRF, with its TO-FILL cells filled with example values, passes `parse_sdrf validate-sdrf -t human -t dia-acquisition` with ontology checks ("Everything seems to be fine"). Its bundled dia-acquisition 1.1.0 accepts **only** `Data-independent acquisition` (rejects `NT=…;AC=PRIDE:0000450` and diaPASEF, which the spec and the GitHub template allow) — so the plain label is written. |
| OLS4 (https://www.ebi.ac.uk/ols4) and unimod.xml | Instrument, enzyme, acquisition and label accessions; Unimod names, accessions and monoisotopic masses. |
| DIA-NN README (vdemichev/DiaNN) / Sage DOCS.md v0.14.7 | `--cut` syntax ("K*,R*,!*P" = canonical tryptic), `--var-mod` sites (`*n` protein N-term); Sage `restrict`, terminus keys `^ $ [ ]`. |
| HPC@UCD docs — https://docs.hpc.ucdavis.edu/data-transfer/ and …/software/ondemand/ | Hive has Globus v5 ("UC Davis Hive home"; PI shares exported on request) and Open OnDemand with a Hive Desktop (https://ondemand.hive.hpc.ucdavis.edu), which runs as a Slurm job. |
| Globus FAQ — https://docs.globus.org/faq/transfer-sharing/ | Symbolic links inside a transferred folder are skipped (hence `COPY_RAW=1` for Globus). |
| Measured on HIVE, 2026-09-24 | Compute node (publicgrp/low): OpenJDK 21.0.12.1; TCP to `hx-fasp-1.ebi.ac.uk:33001`, `ftp-pride-private.ebi.ac.uk:21`, `massive-ftp.ucsd.edu:21` open (login node too). `lftp`, `zip`, `sha1sum`, `md5sum` present; `module pigz/2.8`, `aspera-cli/3.7.7`. `low` allows 7 days, `high` 30. |

## Not verified

- Aspera's UDP data channel from HIVE (only the TCP control port was tested); FTP passive data
  ports.
- Whether the PRIDE Submission Tool's file chooser follows the symlinks the script stages for
  single-file raws (it reads files through the OS, so it should; `COPY_RAW=1` avoids the question).
- Whether MassIVE reserves a PXD for a private dataset; whether it prefers `.d` folders zipped.
- The PRIDE experiment-type drop-down values (the tool only says they are PRIDE CV terms).

## When things change

Re-read the pages above, then update `PX_TOOL_VERSION`, the protocol/title limits and
`ACQ_TERMS` in `make_deposit.py` (and this table's dates). If a newer `sdrf-pipelines` accepts
`NT=diaPASEF;AC=PRIDE:0000650` in a DIA SDRF, write that instead of the parent term.
