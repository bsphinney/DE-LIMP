# FRAN handover — every Core HIVE search reaches the corpus, automatically

[FRAN](https://fran.stan-proteomics.org) is the UC Davis Proteomics Core's cross-search
corpus: every precursor the facility has ever identified, queryable across searches,
species and instruments. Its value is entirely a function of coverage — a corpus with
holes answers "has anyone ever seen this peptide?" with a *no* that is really a *don't
know*. Handing searches over by hand guarantees holes, because it happens for the ones
somebody remembered.

So: **the orchestrator hands over every eligible HIVE search, without being asked.** The
user does not have to request it and is not prompted for it; they are *told* it happened.

## The skill does not ingest anything

FRAN runs its own ingest cron on HIVE that scans for new searches. All the skill has to do
is put the search where that scan will find it:

> **The drop directory is `/quobyte/proteomics-grp/fran/incoming/`** (`FRAN_DROP_DIR`).

So `fran_deposit.py` **stages** — symlinks, never copies — and stops. No database, no
credential, no SLURM job. That is why it works for **every** Core member rather than only
whoever owns the corpus token.

## Who is handed over — and who must never be

**Proteomics Core searches only.** A collaborator with their own HIVE account is running
their own data through this skill; their results are theirs and never enter the Core corpus.

The gate is **write permission on the drop directory**, which lives inside
`/quobyte/proteomics-grp`. It is not a courtesy flag: a HIVE account outside `proteomics-grp`
physically cannot create an entry there. The filesystem enforces the policy; `fran_deposit.py`
only reports it. (Same signal `check_access.sh` reports as `proteomics_grp_access`.)

## The three commands

All run **on HIVE** — in `hive_remote` mode through `hive_exec.sh`, like every other
HIVE-side step. Each prints one JSON object.

```bash
# 1. eligible? (stat only — login-node safe, parses nothing)
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/fran_deposit.py check --out <hive search out dir>'

# 2. hand it over (symlinks; instant, regardless of search size)
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/fran_deposit.py stage \
    --out <hive search out dir> --organism "<the organism the USER confirmed>" --taxon <taxid> \
    --name "<analysis name>"'

# 3. later — did the cron take it?
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/fran_deposit.py verify --out <hive search out dir>'
```

`verify` returns one of three states. **`staged_pending_cron` is a success**, not a failure —
it is the expected state for the minutes or hours between the run finishing and the cron's
next scan. Report it as "handed to FRAN; it ingests on the next pass". The only bad state is
a `broken_links` warning: an entry whose targets have moved looks staged but the cron will
skip it in silence. Re-run `stage --force`.

## What gets linked

| | |
|---|---|
| DIA-NN | `report.parquet` / `report.tsv`, `report.log.txt`, `report.stats.tsv` |
| FragPipe | `dia-quant-output/` (holds its DIA-NN `report.tsv`) |
| Radiant | `radiant_results/` (holds `fulcrum-results`), `delimp_report.parquet` |
| any DIA-NN | `report_xic/` — the chromatograms, **flattened** (see below) |
| plus | `fran_manifest.json`, the only real file in the entry |

**`fran_manifest.json` carries what the cron cannot derive** — the real `output_dir`, the engine
that genuinely ran, the confirmed **organism**, and the search **database**:

| field | why it cannot be inferred |
|---|---|
| `fasta_path` | FRAN parses `--fasta` out of an engine log as a fallback, but that is best-effort and, for Spectronaut, depends on an `ExperimentSetupOverview` the export may not have kept |
| `fasta_md5` | fingerprints *which build* of a proteome, which the filename does not |
| `fasta_n_proteins` | ÷ distinct genes = **entries-per-gene**, and that is what separates a real depth difference from database redundancy — near 1.00 for one-protein-per-gene, above 2 for a full proteome with unreviewed isoforms |

All three come from the `<fasta>.meta.json` `fetch_fasta.py` writes, so the skill hands over what
it already knows rather than making the corpus guess. Absent, never invented, when that file is
missing: `delimp_searches.fasta_*` were populated for 157 / 0 / 0 of 2,014 searches before this,
and a NULL is honest where a guessed database is a claim about comparability.

The entry count is read from the sidecar's `n_entries`, falling back to **`n_sequences`** — the
count `fetch_fasta.py` has always written (proteome + appended contaminants). They are the same
number, verified against a real sidecar: `n_sequences` 34,306 against 34,306 headers counted in
the file. That matters for `check`, which is meant to stat rather than parse: reading the count
that is already there means no sidecar written before these fields existed has to be re-scanned.
Only a sidecar carrying neither count causes the FASTA to be read (measured on a HIVE login node:
0.6 s for 40 MB, 1.3 s for 118 MB), and then only if the file is still where the sidecar says.

**XICs ride along, but they are not in one place.** Every DIA-NN search this skill runs
extracts chromatograms (`--xic` is forced into the cfg — SKILL.md step 6b), and where they land
depends on the route:

| route | where DIA-NN writes them |
|---|---|
| single-shot | `<out>/report_xic/` |
| **5-step parallel chain** (the default above 5 files) | `<out>/xic/t<N>_xic/` — **one directory per array task** |
| FragPipe | `<out>/dia-quant-output/report_xic/` |

DIA-NN names the directory after `--out`, and step 4 runs per file, so a 399-run cohort leaves
chromatograms in **399 separate directories with nothing at `report_xic/`**. `stage` flattens
them into a single `report_xic/` of symlinks in the drop entry, so FRAN's
`diann_xic_to_lance.py --dir <entry>` works with its default path from every route. Safe because
DIA-NN names each file after its run — verified on a real 399-run cohort: 399 files, no basename
collisions, 27 GB of chromatograms handed over as 201 KB of links.

If a search genuinely has no XICs (a cfg from before this was enforced), `check` reports
`xic.present: false` — a fact about that search, not an error.

**`search_provenance.json` is deliberately not linked.** FRAN's scanner lists it as a
*Radiant* marker and tests Radiant before DIA-NN, while `run_search.py` writes one into every
search directory — so linking it would relabel every staged DIA-NN and FragPipe search as
Radiant. Its full contents go into the manifest instead. Verified against FRAN's own
`detect_engine()` on real staged entries: DIA-NN → `diann`, FragPipe → `fragpipe`, Radiant →
`radiant`.

## Pass the organism — it is the one thing only the skill knows

A DIA-NN `report.parquet` carries **no organism column**. Staged without one, the corpus row
is `NULL` and the search is invisible on FRAN's species page. Spectronaut reports self-resolve
from `PEP.AllOccurringOrganisms`; DIA-NN cannot.

The user already confirmed the organism at step 3, and `fetch_fasta.py` wrote it to
`<fasta>.meta.json` — read automatically, overridable with `--organism`/`--taxon`. Never
invent one: an unknown organism is **absent** from the manifest, not guessed (architectural
rule #2).

## Why `check` refuses (stable codes, never an exception)

| reason | meaning |
|---|---|
| `not_core_facility` | the drop directory is not writable → a collaborator's run. **Correct behaviour, not an error.** Do not work around it. |
| `not_on_hive` | the search directory does not exist here — you are not on HIVE, or the path is the local one |
| `search_incomplete` | no report, or a zero-byte one. A failed or zero-ID search must never be ingested |
| `engine_unsupported` | Sage (DDA) and AlphaDIA have no FRAN corpus adapter. The corpus is DIA |
| `no_drop_dir` | the drop directory does not exist and could not be created |
| `already_staged` | a receipt exists; `--force` to re-stage |
| `opted_out` | `FRAN_DEPOSIT=off` or `--skip` |

An ineligible run is **not** a failure of the analysis. Note it in one line and carry on with
DE — never block, retry, or ask the user to fix it.

## Idempotency

The entry name is deterministic — `<search dir name>__<8 hex of its real path>` — so
re-staging reuses the same path rather than presenting the cron with a second candidate that
would ingest as a duplicate search. Re-staging relinks from scratch, so an entry never ends up
holding a mixture of two runs. A receipt at `<out>/fran_deposit.json` stops a resumed session
staging twice.

The manifest's `output_dir` is the **real** search directory, so `corpus_ingest.py --output-dir`
keys idempotency and provenance on where the search actually lives, not on the handover path.

## Follow-ups
- **`--xic` is not on by default.** Searches only carry chromatograms if the cfg asked for
  them. Turning it on for Core runs is a storage/time decision for the facility, not something
  this skill should switch on silently.
- **Local (non-HIVE) searches** have no drop directory to write to and are not handed over.
