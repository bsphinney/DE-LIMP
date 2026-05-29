# Proteogenomics Workflow (v3.11.0)

Architecture reference for the RNA-seq → custom-FASTA proteogenomics pipeline. Referenced from CLAUDE.md. Detailed implementation patterns (SSH-aware launchers, per-user catalog, status-JSON gotchas, deferred submit) live in `docs/PATTERNS.md` under "Proteogenomics Patterns". Subsystem gotchas are in `docs/GOTCHAS.md`.

## Module files
- `R/server_proteog_builder.R` — orchestrator + UI
- `R/helpers_proteog_assembly.R` — FASTA concat + dedupe
- `R/helpers_rnaseq.R` — stage sbatch generators
- `R/helpers_slims.R` — SLIMS + ENA download launchers

## Public functions
`submit_proteogenomics_build()`, `submit_assemble_only()`, `poll_proteog_build_status()`, `cancel_proteog_build()`, `assemble_proteogenomics_fasta()`. All accept `ssh_config = NULL`.

## 11-stage SLURM chain
`fastp → rrna_filter → star → qc_gate → stringtie → merge → gffcompare → gffread → transdecoder → rewrite → assemble`. Each stage is its own sbatch chained via `--dependency=afterok`. `PROTEOG_STAGE_ORDER` constant is the contract.

## Per-build output
`/quobyte/proteomics-grp/de-limp/rnaseq/<project_name>/` with `status.json` (self-describing schema: `pipeline_id`, `project_name`, `sample_names`, `reference_key`, `read_length_tier`, `stages[]`, `build_metadata`), `sbatch/` scripts, `logs/`, and per-stage output dirs.

## Final FASTA output
`/quobyte/proteomics-grp/de-limp/databases/proteogenomics/<project>_proteogenomics_<YYYY_MM>.fasta`. Predicted ORFs + optional UniProt FASTA concatenated by the `assemble` stage on Hive.

## Reference registry
`/quobyte/proteomics-grp/de-limp/references/registry.json` — one entry per species with `organism`, `build`, `genome_fasta`, `gtf`, `star_index`, `rrna_index` paths. Add references via `references/scripts/build_reference_genome.sh` (downloads from Ensembl, builds STAR + bowtie2 rRNA indexes) → writes to `references/registry_pending/` → merge via `references/scripts/merge_registry_pending.sh`.

## FASTA library catalog auto-registration
`poll_proteog_build_status()` detects `current_stage == "complete"` and writes a new entry to local `~/.delimp_fasta_library/catalog.rds` with `content_type = "proteogenomics"` + 9 extension fields (`proteog_pipeline_id`, `proteog_project_dir`, `proteog_methods_paragraph`, `proteog_sample_names`, `proteog_reference_key`, `proteog_uniprot_fasta`, `proteog_read_length_tier`, `proteog_status_path`). The `library_entry_id` is stored back in `status.json` so polling doesn't re-register.

## Shared FASTA + per-user catalog
FASTAs live on shared Hive storage (anyone can read). Catalog entries live in each user's local `~/.delimp_fasta_library/catalog.rds` (no shared-write conflicts). "Discover from Hive" on the search page's Proteogenomics DBs modal scans `PROTEOG_RNASEQ_ROOT/*/status.json` and populates the local catalog with any completed build it lacks. "Restore from Hive" on the Build Database page restores in-progress builds.

## Active builds persistence
`~/.delimp_proteog_builds.rds` mirrors `values$proteog_build_jobs`. The save observer is gated against overwriting a non-empty file with an empty list (so a fresh Shiny session doesn't erase prior state before Discover runs). Restore is wrapped in `isolate({})` because reactiveValues access outside a reactive consumer throws.

## SSH-aware download launchers
`launch_slims_download()`, `launch_ena_download()`, `poll_download_status()` accept `ssh_config = NULL`. When set: `ssh_exec mkdir`, `scp_upload` the status JSON + script, then `nohup bash ... </dev/null &` over SSH. The `</dev/null` is mandatory — without it the SSH connection waits on the background process's stdin.

## Deferred build submit (sra/slims modes)
The submit handler does NOT call `submit_proteogenomics_build()` immediately; it queues in the `pending_build_submits` reactiveVal and adds a placeholder row with `download_pending=TRUE` (blue "downloading" badge) to Active Builds. A 15s `reactivePoll` observer watches `download_status.json` and fires the build when `state="complete"`. Skip-if-present: before launching any download, `find`s for `*.fastq.gz` under the target project_dir; if present, treats it as local mode.

## Per-row Assemble button
Shown in the Active Builds Action column when stage 11 (`assemble`) has `status: "unknown"` and no `job_id`. Opens a modal with the UniProt source dropdown (None / Download from UniProt / Download from NCBI / Enter path on Hive). When the user picks "Download from UniProt" and the upload completes, the assemble job auto-submits without a second click (gated on `proteog_assemble_target` being non-NULL).

## UniProt/NCBI download integration
Reuses `search_uniprot_proteomes()`, `download_uniprot_fasta()`, `ncbi_search_assemblies()`, `ncbi_download_proteome()` via proteog-prefixed observers (`proteog_uniprot_state`, `proteog_ncbi_state`). Downloads land in `/quobyte/proteomics-grp/de-limp/databases/uniprot/` on Hive (created via SCP). NCBI also uploads the side-car `_gene_map.tsv` since RefSeq headers don't embed gene symbols.
