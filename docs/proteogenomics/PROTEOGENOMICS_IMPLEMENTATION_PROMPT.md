# Implementation Prompt: Proteogenomics Database Builder — DE-LIMP Shiny Integration

> **For**: Claude Code session implementing the validated proteogenomics pipeline into DE-LIMP
> **Author**: Brett Phinney / UC Davis Proteomics Core
> **Status**: Ready to execute. Architecture validated end-to-end May 20, 2026.
>             Reconciled against actual DE-LIMP codebase May 21, 2026.
> **Pairs with**: `PROTEOGENOMICS_DB_BUILDER_SPEC_v1.1.md`, `CLAUDE_EXPORT_PROTEOGENOMICS_ADDENDUM_v1.1.md`

## Codebase reconciliation note (May 21, 2026)

The original draft of this spec assumed file names from a Galaxy-P-style template
(`server_diann_search.R`, `helpers_slurm.R`, `ui_search.R`, `inst/extdata/`).
DE-LIMP's actual layout differs:

- DIA-NN HPC orchestration lives in **`R/server_search.R`** (~8,500 lines).
- SLURM/SSH helpers live in **`R/helpers_search.R`** (~5,000 lines, not a separate
  `helpers_slurm.R`). Do not refactor — add new helpers to dedicated files
  (`helpers_rnaseq.R`, `helpers_proteog_qc.R`) and reuse the existing helpers
  in `R/helpers_search.R` directly.
- UI is monolithic in **`R/ui.R`** (~2,600 lines, single `build_ui()`); there is
  no `ui_search.R`.
- DE-LIMP is run via `shiny::runApp()` and is **not an installed R package**.
  `inst/extdata` and `system.file("extdata", ..., package = "delimp")` conventions
  do not apply. The glossary text and helper scripts live under `scripts/` at the
  repo root, resolved relative to `app_dir` using the same pattern as
  `get_contaminant_fasta(library_name, app_dir = NULL)` in `R/helpers_search.R`.
- Navigation decision (May 21, 2026): the existing top-level **"New Search"
  `nav_panel` is being converted into a `nav_menu("New Search", …)` dropdown**
  containing the existing Run Search workflow plus the new Build Database 🧬 panel.
  This mirrors the Comparator's home-in-a-dropdown pattern.

---

## Context for the implementing agent

You are implementing a major feature in the DE-LIMP Shiny app — a "Build Database"
tab that lets users construct proteogenomics search databases from RNA-seq data on
the UC Davis Hive HPC cluster.

The architecture is **fully validated**. A complete end-to-end test run on May 20, 2026
produced a working `predicted_orfs.fasta` (67,386 entries, 100% parse-clean headers,
1,340 NOVEL_GENE discoveries). Every tool choice, every parameter, every file path,
and every quality gate has been verified against real data on Hive.

Your job is to take that validated architecture and turn it into Shiny code that
orchestrates the pipeline from a clickable UI. **You are not researching, you are
implementing.** When in doubt, refer to the spec — it captures the answers to every
ambiguous question that came up during validation.

### What you should NOT do

- Do NOT re-design the pipeline shape. fastp → bowtie2 → STAR → stringtie → gffcompare → gffread → transdecoder → header rewrite is fixed.
- Do NOT introduce new tool choices. The spec specifies module versions verified on Hive.
- Do NOT build an Apptainer container. Validation confirmed Hive's central modules cover everything.
- Do NOT write your own header rewriter from scratch. The validated script at `/quobyte/proteomics-grp/de-limp/pipeline_test/proteog_v3_relaxed_rrnafilt/` produced 100% parse-clean output — port that, don't reinvent it.
- Do NOT skip quality gates. The spec specifies four mandatory gates; all must be enforced.

### What you SHOULD do

- Read the spec in full before writing any code
- Use the existing DIA-NN integration as your architectural template — same SLURM submission pattern, same status polling, same registry idiom
- Mirror the existing addendum patterns (Spectronaut, DDA/DIA, MOFA2) for consistency
- Test incrementally, smallest unit first
- Halt and ask if you encounter anything the spec doesn't cover

---

## Reading order before starting

Read these files in this order:

1. **`PROTEOGENOMICS_DB_BUILDER_SPEC_v1.1.md`** — the architecture and parameters
2. **`CLAUDE_EXPORT_PROTEOGENOMICS_ADDENDUM_v1.1.md`** — the downstream report integration
3. **`DIANN_HPC_INTEGRATION_SPEC.md`** — the architectural pattern you're mirroring
4. **`DIANN_SEARCH_INTEGRATION_SPEC.md`** — the New Search tab you're integrating with
5. **`DE-LIMP_CORE_FACILITY_SPEC.md`** — the multi-user infrastructure you're plugging into
6. **`NAV_RESTRUCTURE_SPEC.md`** — where in the nav this tab lives

Then check the existing DE-LIMP codebase for:

- `R/server_search.R` — existing pattern for SLURM job submission (DE-LIMP's DIA-NN search orchestrator; ~8,500 lines)
- `R/helpers_search.R` — existing SLURM/SSH helpers (`ssh_exec()`, `generate_sbatch_script()`, `parse_sbatch_output()`, `check_slurm_status()`, `scan_prestaged_databases()`); ~5,000 lines
- `R/ui.R` — monolithic UI file containing every `nav_panel`; the Build Database panel is added to `build_ui()` here, not in a separate file

**Do not refactor SLURM helpers into a separate file.** Earlier versions of this spec
proposed extracting them into `R/helpers_slurm.R`; the maintainer has confirmed that
churn is out of scope. Add new helpers to dedicated files (`helpers_rnaseq.R`,
`helpers_proteog_qc.R`, etc.) and *reuse* the existing SSH/SLURM helpers in
`R/helpers_search.R` directly.

---

## Project conventions discovered during validation (lessons #14–#17)

The May 20–21, 2026 validation cycle established four project conventions
that future Phase D/E work must follow. Full text in
`NOTES_spec_lessons.md` on Hive (`/quobyte/proteomics-grp/de-limp/pipeline_test/NOTES_spec_lessons.md`).

### Lesson #14 — Functions that write output must integrity-check inputs

Any helper that writes output FASTAs and takes input FASTAs MUST snapshot
its inputs at function entry (size + mtime + first-1KB md5) and verify
them before return. Cheap defense (~5 ms/input) against the class of
"function corrupted my source data" bugs invisible at code-review time.

**Concrete pattern** (already implemented in `assemble_proteogenomics_fasta()`):

```r
input_snapshots <- list(
  uniprot  = .snapshot_input(uniprot_fasta),
  predict  = .snapshot_input(predicted_orfs_fasta),
  ...
)
on.exit({
  for (s in input_snapshots) .verify_input_unchanged(s, strict = FALSE)
}, add = TRUE)
# ... do work ...
for (s in input_snapshots) .verify_input_unchanged(s, strict = TRUE)
```

`.snapshot_input`, `.verify_input_unchanged`, `.head_md5`, `.log_disk_write`
all live in `R/helpers_proteog_assembly.R`. Reuse them in any new helper
that writes output. Also reuse `.validate_fasta_input(path, label)` to
catch malformed (e.g., 2-byte "0\n") inputs AT ENTRY before downstream
tools see them.

### Lesson #15 — Document the upstream source of every static asset

Contaminant FASTAs, reference proteomes, and any other curated file shipped
with DE-LIMP must carry provenance metadata: a `provenance.json` sibling in
the same directory, OR at minimum a `# SOURCE: <url>` comment at the top
of the file. When asked "where did this file come from?" the answer should
be one filesystem lookup, not "I don't remember."

**Concrete pattern**:
- `contaminants/provenance.json` documents all 6 HaoGroup contaminant FASTAs
- `/quobyte/.../fasta/UP000000589_mus_musculus_opg_2026_05.provenance.json`
  documents the recovered UniProt OPG with download URL, md5, and notes
- `load_asset_provenance(asset_path)` in `helpers_proteog_assembly.R` is
  the canonical reader

### Lesson #16 — Prefer R native over `system2()` shell-outs

When R has a builtin that does the work (`grepl`, `readLines`, `file.size`,
`nchar`, `Biostrings::fasta.index`), use it. Reach for `system2()` only
when the work genuinely requires an external binary (`sbatch`, `samtools`,
`seqkit`, etc.).

When `system2()` IS necessary, every string argument that isn't a known-safe
flag or constant must go through `shQuote()`. No exceptions.

**Why this matters**: validation surfaced a serious data-corruption
incident where `system2("grep", c("-c", "^>", f))` was parsed by bash as
`grep -c ^> /path/to/f` — interpreting `^>` as redirect operator. Three
production FASTAs were silently truncated to "0\n" before the calling
function even ran. R-native equivalents (`sum(grepl("^>", readLines(f)))`)
make this class of bug impossible.

### Lesson #17 — Maintain a canonical demo dataset

The "demo_mouse_reproducibility" build is the project's known-good
reference. Any change to assemble / orchestrator / rewriter / qc gates
should reproduce its post-recovery composition (~75,465 entries total,
0 UNPARSED) before merge.

| Field | Value |
|---|---|
| Raw FASTQ source | `/quobyte/proteomics-grp/de-limp/pipeline_test/sra_data/` |
| Samples | `SRR1303776`, `SRR1303777` (mouse, ENCODE 2014) |
| Reference | `mm39_GRCm39` (GENCODE vM38) |
| Contaminants | `Mouse_Tissue_Contaminants.fasta` (HaoGroup) |
| STAR tier | `significantly_relaxed` (92 bp reads) |
| Expected final FASTA | ~75,465 entries, 0 UNPARSED, ~58 NOVEL_GENE, ~1,360 NOVEL_ISOFORM |
| Expected wall | ~40 min on `high` partition |

UNPARSED must always be 0; class counts may drift ±100 between runs but
not orders of magnitude.

### Header rewriter — gffcompare class code map (final, locked)

The rewriter (`scripts/rewrite_transdecoder_headers.py`) maps gffcompare
class codes to `source=` classes. As of Phase B.5 the map covers all 16
gffcompare class codes observed in real mouse data — including the four
added during validation:

```python
CLASS_CODE_MAP = {
    "=": "REF",            "c": "REF",
    "o": "REF",            "x": "REF",  "s": "REF",
    "p": "REF",            "r": "REF",
    "j": "NOVEL_ISOFORM",  "e": "NOVEL_ISOFORM",
    "k": "NOVEL_ISOFORM",  # query CONTAINS reference (UTR extension)
    "m": "NOVEL_ISOFORM",  # retained intron, full chain
    "n": "NOVEL_ISOFORM",  # retained intron, partial chain
    "y": "NOVEL_ISOFORM",  # contains reference within intron (rare)
    "u": "NOVEL_GENE",     "i": "NOVEL_GENE",  # i = intronic (uORF/sORF candidates)
    ".": "UNPARSED",       # gffcompare couldn't classify
}
```

Any unmapped class code → UNPARSED → exit non-zero per Rule 4. If a future
gffcompare run surfaces a new code, the rewriter's pre-flight check WARNS
before classification; extend the map with documented reasoning (lesson
#13: locked constants are conservative defaults, not immutable laws).

---

## Phased implementation plan

This is a substantial feature. Implement in phases and **verify each phase works
before starting the next**. Each phase should end with a working, testable state.

### Phase A — Infrastructure setup (Day 1, ~2-3 hours)

Goal: directories, conda env, reference data registration. No UI yet, no Shiny code.

#### A.1 — Verify Hive prerequisites

SSH to Hive (you have access via the user's account). Verify:

```bash
# Check that the canonical directory tree exists
ls -la /quobyte/proteomics-grp/de-limp/

# Should see at minimum:
#   containers/de-limp.sif
#   references/   (may be empty)
#   envs/         (may be empty)

# Check modules are available
module avail fastp bowtie2 star stringtie gffcompare gffread transdecoder diamond seqkit 2>&1
```

If any of these are missing, **stop and ask** — that would mean the validation
environment has changed and the spec needs updating before continuing.

#### A.2 — Set up the proteog_helpers conda env

If `/quobyte/proteomics-grp/de-limp/envs/proteog_helpers/` already exists from
validation, skip this. Otherwise:

```bash
module load conda

# IMPORTANT: -c conda-forge -c bioconda --strict-channel-priority is required.
# Reverse order picks python-2-only biopython and fails.
conda create -p /quobyte/proteomics-grp/de-limp/envs/proteog_helpers \
  -c conda-forge -c bioconda --strict-channel-priority -y \
  python=3.11 biopython gffutils

# Verify
source activate /quobyte/proteomics-grp/de-limp/envs/proteog_helpers
python -c "import gffutils, Bio; print('OK:', gffutils.__version__, Bio.__version__)"
```

#### A.3 — Set up reference data tree

Create the directory structure documented in spec §13:

```bash
mkdir -p /quobyte/proteomics-grp/de-limp/references/{genomes,rrna_index,star_index,gtf,scripts}
mkdir -p /quobyte/proteomics-grp/de-limp/databases/proteogenomics
mkdir -p /quobyte/proteomics-grp/de-limp/rnaseq

# Initialize registry
echo '{}' > /quobyte/proteomics-grp/de-limp/references/registry.json
echo '{}' > /quobyte/proteomics-grp/de-limp/databases/proteogenomics/registry.json
```

For mouse (already complete in bioinfocore-grp), create symlinks:

```bash
ln -s /quobyte/bioinfocore-grp/genomes/mouse/GRCm39/STAR_GRCm39_vM38 \
  /quobyte/proteomics-grp/de-limp/references/star_index/mm39

ln -s /quobyte/bioinfocore-grp/genomes/mouse/GRCm39/gencode.vM38.basic.annotation.gtf \
  /quobyte/proteomics-grp/de-limp/references/gtf/mm39.gtf

ln -s /quobyte/bioinfocore-grp/genomes/mouse/GRCm39/GRCm39.primary_assembly.genome.fa \
  /quobyte/proteomics-grp/de-limp/references/genomes/mm39_GRCm39_genome.fa
```

For human, download the genome FASTA (bioinfocore-grp doesn't stage it):

```bash
cd /quobyte/proteomics-grp/de-limp/references/genomes
mkdir -p hg38_GRCh38.p14
cd hg38_GRCh38.p14
wget -O genome.fna.gz \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip genome.fna.gz
mv genome.fna genome.fa
module load samtools/1.19.2
samtools faidx genome.fa

# Symlinks for human
cd /quobyte/proteomics-grp/de-limp/references
ln -s /quobyte/bioinfocore-grp/genomes/human/GRCh38.p14/STAR_2.7.11b_index star_index/hg38
ln -s /quobyte/bioinfocore-grp/genomes/human/GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf gtf/hg38.gtf
```

#### A.4 — Build the rRNA bowtie2 indices

```bash
module load bowtie2/2.5.2
mkdir -p /quobyte/proteomics-grp/de-limp/references/rrna_index/{hg38,mm39}

bowtie2-build \
  /quobyte/bioinfocore-grp/genomes/human/GRCh38.p14/rRNA_human_03-12-2026.fasta \
  /quobyte/proteomics-grp/de-limp/references/rrna_index/hg38/rrna

bowtie2-build \
  /quobyte/bioinfocore-grp/genomes/mouse/GRCm39/rRNA_mouse_03-12-2026.fasta \
  /quobyte/proteomics-grp/de-limp/references/rrna_index/mm39/rrna
```

#### A.5 — Populate the references registry

Write to `/quobyte/proteomics-grp/de-limp/references/registry.json`:

```json
{
  "hg38_GRCh38.p14": {
    "organism": "Homo sapiens",
    "build": "GRCh38.p14",
    "annotation_source": "RefSeq",
    "annotation_release": "GCF_000001405.40",
    "genome_fasta": "/quobyte/proteomics-grp/de-limp/references/genomes/hg38_GRCh38.p14/genome.fa",
    "star_index": "/quobyte/proteomics-grp/de-limp/references/star_index/hg38",
    "gtf": "/quobyte/proteomics-grp/de-limp/references/gtf/hg38.gtf",
    "rrna_index": "/quobyte/proteomics-grp/de-limp/references/rrna_index/hg38/rrna",
    "completeness": "complete",
    "registered": "2026-05-21",
    "registered_by": "<your user>"
  },
  "mm39_GRCm39": {
    "organism": "Mus musculus",
    "build": "GRCm39",
    "annotation_source": "GENCODE",
    "annotation_release": "vM38",
    "genome_fasta": "/quobyte/proteomics-grp/de-limp/references/genomes/mm39_GRCm39_genome.fa",
    "star_index": "/quobyte/proteomics-grp/de-limp/references/star_index/mm39",
    "gtf": "/quobyte/proteomics-grp/de-limp/references/gtf/mm39.gtf",
    "rrna_index": "/quobyte/proteomics-grp/de-limp/references/rrna_index/mm39/rrna",
    "completeness": "complete",
    "registered": "2026-05-21",
    "registered_by": "<your user>"
  }
}
```

#### A.6 — Stage the header rewriter script

Copy the validated header rewriter from the test run to its production location:

```bash
# Source: validation output
SRC=/quobyte/proteomics-grp/de-limp/pipeline_test/proteog_v3_relaxed_rrnafilt/
DEST=/quobyte/proteomics-grp/de-limp/references/scripts/

# Find the actual rewriter script in the validation output
find $SRC -name "rewrite_*.py" -o -name "*_headers*.py"

# Copy it
cp <found_script> $DEST/rewrite_transdecoder_headers.py
chmod 755 $DEST/rewrite_transdecoder_headers.py

# Verify it runs (should print usage)
source activate /quobyte/proteomics-grp/de-limp/envs/proteog_helpers
python $DEST/rewrite_transdecoder_headers.py --help
```

If the validation script needs documentation/cleanup before being production-grade,
do that now. Required interface:

```
python rewrite_transdecoder_headers.py \
  --transdecoder <input.pep> \
  --merged-gtf <stringtie_merged.gtf> \
  --gffcompare-tmap <gffcmp.merged.gtf.tmap> \
  --ref-gtf <reference.gtf> \
  --project-tag <TAG> \
  --output <output.fasta>
```

**Verification at end of Phase A**: print a summary of what was created:

```
✓ Conda env: /quobyte/proteomics-grp/de-limp/envs/proteog_helpers (biopython 1.87, gffutils 0.14)
✓ Reference registry: 2 entries (hg38, mm39)
✓ rRNA indices: 2 (hg38, mm39)
✓ Header rewriter: /quobyte/proteomics-grp/de-limp/references/scripts/rewrite_transdecoder_headers.py
✓ Directory tree created for rnaseq inputs, databases output
```

If you cannot get to this clean checkpoint, halt and report what went wrong.

---

### Phase B — R helpers and sbatch generators (Day 1-2, ~4-6 hours)

Goal: pure-function R code that generates sbatch scripts, parses outputs, manages
registries. No Shiny reactivity yet, no UI.

#### B.1 — Create `R/helpers_proteogenomics.R`

This file replaces/supersedes whatever the existing addendum stub is. Contains:

```r
# Classification of proteins by source
classify_proteins <- function(diann_report) { ... }

# Build helpers for Claude export prompts
build_proteog_note         <- function(values) { ... }
build_proteog_file_note    <- function(values) { ... }
build_proteog_section      <- function(values, template_type) { ... }  
build_biosynth_proteog_note <- function(values) { ... }
build_proteog_inline       <- function(values, input) { ... }
```

Follow the addendum v1.1 spec exactly for the prompt-block contents. Each helper
returns `""` when `!isTRUE(values$is_proteogenomics)`.

Write unit tests for `classify_proteins()` covering:
- Pure UniProt input → all UNIPROT class
- Pipeline output FASTA → REF/NOVEL_GENE/NOVEL_ISOFORM classes via source= parsing
- Mixed input with VARIANT entries → VARIANT class detected by prefix
- Malformed descriptions → graceful handling, default to UNKNOWN or UNIPROT

#### B.2 — Create `R/helpers_slims.R`

```r
# SLIMS URL scanning (synchronous, fast HTTP)
scan_slims_url <- function(slims_url) { ... }

# ENA accession metadata verification
verify_sra_accession <- function(accession) { ... }

# Reference registry I/O
load_reference_registry <- function() { ... }
load_proteog_registry   <- function() { ... }

# Background download launcher (login-node only)
launch_slims_download   <- function(slims_url, project_name) { ... }
launch_ena_download     <- function(accessions, project_name, subsample_reads = NULL) { ... }
```

The download launchers run on the login node via `nohup` since compute nodes
may not have outbound HTTP. They write a status file the UI can poll.

Test: with a fake SLIMS URL pattern, verify `scan_slims_url()` returns the
expected error structure. With a real ENA accession (any small one), verify
`verify_sra_accession()` returns parseable metadata.

#### B.3 — Create `R/helpers_rnaseq.R`

The big one. Generates sbatch scripts for each pipeline stage:

```r
generate_fastp_sbatch <- function(project_dir, sample_names, slurm_account, slurm_partition) { ... }
generate_rrna_sbatch  <- function(project_dir, sample_names, rrna_index_path, slurm_account, slurm_partition) { ... }
generate_star_sbatch  <- function(project_dir, sample_names, star_index, tier_params, slurm_account, slurm_partition) { ... }
generate_stringtie_sbatch <- function(project_dir, sample_names, ref_gtf, strand_flag, slurm_account, slurm_partition) { ... }
generate_merge_sbatch <- function(project_dir, ref_gtf, slurm_account, slurm_partition) { ... }
generate_gffcompare_sbatch <- function(project_dir, ref_gtf, slurm_account, slurm_partition) { ... }
generate_transdecoder_sbatch <- function(project_dir, genome_fasta, diamond_db = NULL, slurm_account, slurm_partition) { ... }
generate_rewrite_sbatch <- function(project_dir, project_tag, rewriter_path, conda_env_path, slurm_account, slurm_partition) { ... }
```

Each takes the project directory and produces a string (sbatch script content)
that can be written to disk and submitted. **Use the exact module versions from
the spec.** Use `glue::glue()` for templating, NOT `paste0()` with `sprintf` — the
spec's sbatch examples assume `glue` semantics.

Critical: implement `select_star_params()` exactly as specified in spec §5.3.
The tier thresholds are validation-derived and must not drift.

Tests:
- Generate each sbatch type for a fake project with 3 samples
- Verify the output sbatch is syntactically valid (e.g., `sbatch --test-only`
  via SSH to Hive if possible)
- Verify all `#SBATCH` directives are present
- Verify module load lines are correct

#### B.4 — Create `R/helpers_proteog_qc.R`

```r
# Parse STAR Log.final.out for uniquely-mapped rate
parse_star_log <- function(log_path) { ... }

# Parse fastp JSON for read length, total reads, %passed
parse_fastp_json <- function(json_path) { ... }

# Parse bowtie2 rRNA filter log
parse_rrna_log <- function(log_path) { ... }

# QC gate enforcement
check_alignment_quality <- function(star_log_path, tier_params) { ... }
check_pipeline_gates <- function(project_dir, gate_results) { ... }
```

`check_pipeline_gates()` aggregates all gate results for a run and returns a
combined pass/fail with per-gate diagnostics. The UI consumes this for the
status panel.

#### B.5 — Create `R/helpers_proteog_assembly.R`

```r
assemble_proteogenomics_fasta <- function(...) { ... }  # spec §7
count_proteog_classes         <- function(fasta_path) { ... }
register_proteogenomics_fasta <- function(...) { ... }
load_proteog_registry         <- function() { ... }  # may be duplicate of helpers_slims; pick one location
```

The composition counter parses `source=` tags from headers. Test against the
validation output FASTA to confirm it produces the expected 66046 REF / 1340
NOVEL_GENE / 0 NOVEL_ISOFORM / 0 UNPARSED counts.

#### B.6 — Verification

Before moving to Phase C, run all helpers against the validation output:

```r
# In an R session on Hive (inside the de-limp.sif container)
source("R/helpers_proteogenomics.R")
source("R/helpers_proteog_assembly.R")

fasta <- "/quobyte/proteomics-grp/de-limp/pipeline_test/proteog_v3_relaxed_rrnafilt/mm39_test_92bp_relaxed_RRNAfilt.fasta"
counts <- count_proteog_classes(fasta)
stopifnot(counts$total == 67386)
stopifnot(counts$REF == 66046)
stopifnot(counts$NOVEL_GENE == 1340)
stopifnot(counts$UNPARSED == 0)
```

If these pass, Phase B is done.

---

### Phase C — Server logic and SLURM orchestration (Day 2-3, ~6-8 hours)

Goal: a callable `submit_proteogenomics_build()` function that takes a project
config and submits the full SLURM dependency chain. Still no UI.

#### C.1 — Create `R/server_proteog_builder.R`

```r
# Main entry point — submits the full pipeline as a chained SLURM job array
submit_proteogenomics_build <- function(
  project_name,
  rnaseq_dir,
  reference_key,          # e.g., "mm39_GRCm39" — looked up in registry
  sample_names,
  library_type,            # "polyA" | "ribo_depleted" | "stranded"
  strand_flag,             # "--rf" | "--fr" | "" (unstranded)
  uniprot_fasta,           # path to UniProt reference to merge with predicted ORFs
  slurm_account,
  slurm_partition,
  values                   # reactiveValues for status updates
) {
  # Returns a list of job IDs (one per stage) that the UI can poll
}

# Status polling — called periodically by reactivePoll in app.R
poll_proteog_build_status <- function(job_chain_id) {
  # Returns list with current stage, % complete, any errors
}

# Cancellation — for the "Cancel build" button
cancel_proteog_build <- function(job_chain_id) {
  # scancel all jobs in the chain
}
```

The orchestration:

1. Validate inputs (project name sanitized, all sample files exist, reference exists in registry)
2. Create project subdirectory in `/quobyte/proteomics-grp/de-limp/rnaseq/<project>/`
3. Generate all sbatch scripts using helpers from Phase B
4. Submit them with SLURM dependency flags:
   ```
   jid_fastp = sbatch fastp.sbatch
   jid_rrna  = sbatch --dependency=afterok:$jid_fastp rrna.sbatch
   jid_star  = sbatch --dependency=afterok:$jid_rrna star.sbatch
   ... etc
   ```
5. Insert a QC gate job between STAR and stringtie that fails if uniquely-mapped is low
6. Record the job chain in a SQLite table (or in-memory list for single-user mode)
7. Return the job chain ID

#### C.2 — QC gate as a SLURM job

The cleanest way to enforce QC gates within SLURM is a small Python or bash
script that runs as its own sbatch job, checks the upstream logs, and exits
non-zero if gates fail. Downstream jobs with `--dependency=afterok` will not run
if it fails.

```bash
# qc_gate_alignment.sbatch
#!/bin/bash -l
#SBATCH --job-name=qc_gate_${PROJECT}
#SBATCH --time=10:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Parse STAR logs, check unique-mapped %
python /quobyte/proteomics-grp/de-limp/references/scripts/check_alignment_qc.py \
  --star-logs star_out/*_Log.final.out \
  --threshold ${QC_GATE_THRESHOLD} \
  --output qc_gate_result.json

# Exit code from the Python script propagates to SLURM
```

The Python helper writes a JSON result that DE-LIMP can read to show diagnostics
to the user when a gate fails. The downstream jobs are still queued but never
run; DE-LIMP detects the cancellation via `squeue` and surfaces the failure.

#### C.3 — Status polling integration

Use the existing DIA-NN integration's status polling pattern. The status file
written to disk is what the Shiny `reactivePoll` watches:

```
/quobyte/proteomics-grp/de-limp/rnaseq/<project>/status.json
```

Updated by a SLURM finalization job that runs after each stage:

```json
{
  "project_name": "my_experiment_2026_05",
  "current_stage": "stringtie",
  "stages": [
    {"stage": "fastp",       "status": "complete", "started": "...", "finished": "..."},
    {"stage": "rrna_filter", "status": "complete", "started": "...", "finished": "..."},
    {"stage": "star",        "status": "complete", "started": "...", "finished": "..."},
    {"stage": "qc_gate",     "status": "complete", "passed": true},
    {"stage": "stringtie",   "status": "running",  "started": "..."},
    {"stage": "merge",       "status": "pending"},
    {"stage": "gffcompare",  "status": "pending"},
    {"stage": "transdecoder","status": "pending"},
    {"stage": "rewrite",     "status": "pending"},
    {"stage": "assemble",    "status": "pending"}
  ],
  "qc_metrics": {
    "rrna_pct_mean": 4.2,
    "uniquely_mapped_pct_mean": 78.3,
    "read_length_median": 148,
    "tier": "default"
  }
}
```

#### C.4 — Verification

End of Phase C, you should be able to call `submit_proteogenomics_build()`
from an R console on Hive and have the full pipeline run on a small test case
(2 samples, ~5M reads each). Use the validation data still on disk if it's
convenient (mouse, 5M-pair subsample).

The expected result: a final FASTA in
`/quobyte/proteomics-grp/de-limp/databases/proteogenomics/<project_name>_proteogenomics_<YYYY_MM>.fasta`
with composition counts matching the test data.

**This is the proof-of-life checkpoint.** If you can run the pipeline end-to-end
from R without any UI, the hard work is done.

---

### Phase D — UI integration (Day 3-4, ~4-6 hours)

Goal: a working Build Database tab in the DE-LIMP nav.

#### D.1 — Create `R/ui_proteog_builder.R`

UI structure follows spec §10. Use `bslib`/`shinydashboard` matching the existing
DE-LIMP nav style. Sections:

1. **Step 1: Source selection** — SLIMS URL OR SRA/ENA accessions
2. **Step 2: Sample scan / metadata verification** — shown after Step 1 submit
3. **Step 3: Reference selection** — dropdown of organisms from the registry
4. **Step 4: Pipeline parameters** — library type, strand, project tag, project name
5. **Step 5: Submit** — disabled until all prior steps validated
6. **Active builds panel** — table of in-flight and recent builds with status

Active builds panel uses `reactivePoll` against the per-project `status.json`
files in the rnaseq dir.

#### D.2 — Modify `R/ui.R` (monolithic `build_ui()`)

The existing navbar has **"New Search" as a top-level `nav_panel`**, not a dropdown.
Convert it to a `nav_menu("New Search", ...)` dropdown containing the existing
search workflow plus the new Build Database panel, mirroring the Analysis/Comparator
pattern. Visibility of Build Database is gated by HPC mode detection (the
existing `is_hpc` / `sbatch_available()` check used elsewhere in `R/ui.R`):

```r
nav_menu("New Search", icon = icon("rocket"),
  nav_panel("Run Search",     value = "search_tab", icon = icon("magnifying-glass"),
            uiOutput("run_search_content")),  # the existing "New Search" body
  if (is_hpc && !is_hf_space) {
    nav_panel("Build Database", value = "build_database_tab", icon = icon("dna"),
              uiOutput("build_database_content"))
  }
)
```

Add both `value` strings to the protected list in `CLAUDE.md` (the "Tab values
that MUST NOT change" block):
- `"search_tab"` (renamed/preserved value for the existing search UI body)
- `"build_database_tab"` (new)

Any existing code calling `nav_select("main_tabs", "New Search")` must be
updated to `nav_select("main_tabs", "search_tab")`. Grep for these references
before merging.

#### D.3 — Modify `app.R`

Add reactiveValues per spec §11:

```r
values <- reactiveValues(
  # ... existing ...
  is_proteogenomics      = FALSE,
  protein_classification = NULL,
  proteog_build_jobs     = list(),
  proteog_active_fasta   = NULL
)
```

Add the proteog status poller (one per active build):

```r
observe({
  invalidateLater(15000)  # poll every 15s
  values$proteog_build_jobs <- lapply(values$proteog_build_jobs, function(job) {
    if (job$status %in% c("complete", "failed")) return(job)
    poll_proteog_build_status(job$chain_id)
  })
})
```

#### D.4 — Wire detection at DIA-NN result load

In whatever handler currently runs after `readDIANN()`, append:

```r
values$protein_classification <- classify_proteins(values$diann_report)
values$is_proteogenomics <- any(
  values$protein_classification$source %in% c("REF", "NOVEL_GENE", "NOVEL_ISOFORM", "VARIANT")
)
```

#### D.5 — Wire the Claude export integration

In `R/server_ai.R`, wherever `build_claude_prompt()` is defined, add the
proteog conditionals per addendum v1.1. Build the seven helper functions in
`R/helpers_proteogenomics.R` if not already done.

In the ZIP assembly block, add:

```r
if (isTRUE(values$is_proteogenomics)) {
  write_proteogenomics_glossary(zip_dir)
}
```

The glossary text from addendum v1.1 goes in `scripts/proteogenomics_glossary.txt` at
the repo root, resolved relative to `app_dir` at runtime (DE-LIMP is not an installed
R package — `inst/extdata` package conventions do not apply). Use the same pattern
as `get_contaminant_fasta(library_name, app_dir = NULL)` in `R/helpers_search.R`.

#### D.6 — Update `R/helpers_search.R` (`scan_prestaged_databases()`)

The current `scan_prestaged_databases(fasta_dir)` at `R/helpers_search.R:581` returns a
plain `character()` vector of FASTA file paths. Extend it to read the proteogenomics
registry and emit labeled choices with the 🧬 tag. The actual call site in
`R/server_search.R` consumes the labeled vector directly — no other changes needed
there beyond passing the labels through to `updateSelectInput()`.

```r
scan_prestaged_databases <- function(fasta_dir) {
  files <- list.files(fasta_dir, pattern = "\\.fasta$",
                      full.names = TRUE, recursive = TRUE)
  
  registry_path <- "/quobyte/proteomics-grp/de-limp/databases/proteogenomics/registry.json"
  registry <- if (file.exists(registry_path)) {
    jsonlite::read_json(registry_path)
  } else list()
  
  choices <- setNames(
    files,
    sapply(files, function(f) {
      proteog <- Filter(function(r) r$path == f, registry)
      if (length(proteog) > 0) {
        r <- proteog[[1]]
        sprintf("🧬 %s — Proteogenomics (composition: %s UniProt + %s REF + %s NOVEL_GENE + %s NOVEL_ISOFORM)",
                r$project_name,
                format(r$composition$UNIPROT, big.mark = ","),
                format(r$composition$REF, big.mark = ","),
                format(r$composition$NOVEL_GENE, big.mark = ","),
                format(r$composition$NOVEL_ISOFORM, big.mark = ","))
      } else {
        sprintf("◇ %s", basename(f))
      }
    })
  )
  choices
}
```

Add the auto-warning observer when proteogenomics FASTA is selected (spec §8).

#### D.7 — Verification

End of Phase D, a user should be able to:

1. Open DE-LIMP, navigate to Search → Build Database
2. Paste a SLIMS URL (or SRA accessions)
3. See sample scan + metadata verification
4. Choose mouse reference, default parameters
5. Click Build
6. Watch the status panel update through 8 stages
7. After completion, navigate to Search → New Search and see the new proteogenomics
   FASTA in the pre-staged dropdown with the 🧬 tag
8. Run a DIA-NN search against it
9. Load the results, see `is_proteogenomics = TRUE`
10. Generate a Claude export and confirm `Proteogenomics_Glossary.txt` is in the ZIP

If all 10 steps work end-to-end, the feature is complete.

---

### Phase E — Polish (Day 4-5, ~2-4 hours)

Goal: edge cases, error handling, user-facing improvements.

- Session save/load round-trip for `is_proteogenomics` and `protein_classification`
- Resume-from-stage logic for interrupted builds (use existing SLURM checkpoints)
- Disk space pre-check (warn if /quobyte/proteomics-grp has <250 GB free)
- MultiQC summary generation as a final stage (optional, off by default)
- "View RNA-seq QC report" button in the active builds panel
- Help text / tooltips on every UI element
- Test on Docker backend (everything should gracefully hide; "HPC mode required" message)

---

## Things you will probably get wrong on first attempt

Based on patterns from validation, watch for these:

1. **`set +o pipefail` in stream-subsample sbatch.** Without it, `curl | zcat | head | gzip` 
   fails with exit 141 when head closes the pipe. This was a real validation bug.

2. **`source activate` order in mixed module+conda sbatch scripts.** Module loads
   must precede conda env activation. If you `source activate` first, the module
   load can clobber the PATH.

3. **bgzf crash in stringtie/2.2.1 module.** Use `conda/stringtie/3.0.3`. Do not
   try to work around the 2.2.1 module.

4. **Conda channel order.** `-c conda-forge -c bioconda --strict-channel-priority`. 
   Reverse order picks python-2-only biopython.

5. **`#SBATCH --account` flag.** Must be set explicitly. Default is `publicgrp` 
   which doesn't have `high` partition access.

6. **Compute nodes lack outbound HTTP.** SLIMS download, ENA download, NCBI genome 
   download must run on login nodes (background process via `nohup`), not as sbatch jobs.

7. **`stringtie --merge` collapses novel isoforms into REF.** This is by design. 
   The gffcompare step is what recovers them. Don't skip gffcompare.

8. **Header rewriter must produce ZERO UNPARSED entries.** If validation runs of 
   the rewriter produce any UNPARSED, fix the rewriter before proceeding. Do not 
   accept "we'll fix it later."

9. **Don't trust SRA metadata in project descriptions.** Always verify via ENA 
   XML API. Validation surfaced this when SRRs claimed to be K562 human turned 
   out to be Mus musculus.

10. **PE150 is default for DNA Tech Core, but validation used 92bp data.** 
    Make sure the read-length adaptive tier selection works on both — test on 
    real DNA Tech Core data (PE150) if at all possible during validation, not 
    only on the historical 92bp test data.

---

## Status reporting back

At the end of each phase, produce a brief report:

```
PHASE A COMPLETE — Infrastructure setup
  Conda env: ✓
  References: ✓ (hg38, mm39)
  rRNA indices: ✓
  Header rewriter staged: ✓
  Time elapsed: 2h 15min

PHASE B COMPLETE — R helpers and sbatch generators
  Files created: 5
  Unit tests passing: 23/23
  Composition counter validates against test FASTA: ✓
  Time elapsed: 5h 30min

[etc.]
```

If any phase produces unexpected results, **halt and report** rather than improvise.
The validation work surfaced 11 spec-quality bugs by halting at the first sign of
trouble. Continue that practice.

---

## Final deliverable

When all phases complete, the feature should be:

- ✅ A working "Build Database" tab in DE-LIMP's Search nav
- ✅ End-to-end SLIMS-URL → predicted_orfs.fasta workflow on Hive
- ✅ Adaptive STAR thresholds with QC gating
- ✅ rRNA pre-filter as a mandatory step
- ✅ Species verification for SRA inputs
- ✅ Header rewriter producing 100% parse-clean output
- ✅ Proteogenomics FASTA auto-detected as 🧬 in pre-staged dropdown
- ✅ Claude export ships `Proteogenomics_Glossary.txt` when proteogenomic session loaded
- ✅ Three-template (Brief / Full / Manuscript) prompt blocks active
- ✅ Inline summaries of non-canonical hits in the export context

The user (Brett) should be able to demo the feature to a bench scientist:
"Submit your RNA-seq to DNA Tech Core, paste the SLIMS URL into DE-LIMP,
click build, wait a few hours, then search your proteomics data against the
expanded database — no command line, no bioinformatics expertise required."

That's the user value. Everything else is plumbing.

---

*Implementation prompt v1.0 — Brett Phinney / UC Davis Proteomics Core — May 2026*
*Built on validated pipeline architecture from May 20, 2026 end-to-end test run*
