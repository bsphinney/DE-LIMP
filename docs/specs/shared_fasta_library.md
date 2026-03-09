# Shared FASTA Database Library — Spec

## Problem
Every time someone runs a search, they download the same FASTA from UniProt or manually point to a file. There's no central catalog of what databases exist, what's in them, or which search settings they're compatible with. Lab members waste time re-downloading and can't easily share curated databases that include contaminants or custom proteins.

## Overview
A shared, browsable library of pre-built FASTA databases stored on the proteomics network volume. Users can select a database from a catalog instead of downloading from UniProt, with full metadata about what's in each database.

## User Flow

### New dropdown option
The existing FASTA source selector gets a third option:

```
FASTA Source: [UniProt Download] [Upload File] [Database Library]  <-- new
```

### Database Library modal
Selecting "Database Library" opens a modal/panel showing a searchable table:

| Name | Organism | Proteins | Contaminants | Custom | Created | Created By |
|------|----------|----------|--------------|--------|---------|------------|
| Human OPG + Contam | Homo sapiens | 20,598 | Universal | — | 2026-03-01 | Brett |
| Human Full + Contam | Homo sapiens | 83,526 | Universal | — | 2026-02-15 | Brett |
| Mouse OPG + Contam | Mus musculus | 17,142 | Universal | — | 2026-02-20 | Taha |
| Human OPG + Contam + iRT | Homo sapiens | 20,609 | Universal | 11 iRT peptides | 2026-03-02 | Brett |

### Detail panel (on row select)
Selecting a row expands details:

- **Organism**: Homo sapiens (Human)
- **UniProt proteome**: UP000005640
- **Content type**: One protein per gene
- **Protein count**: 20,598 sequences
- **File size**: 12.3 MB
- **Contaminant library**: Universal (246 proteins)
- **Custom sequences**: None
- **FASTA files**:
  - `UP000005640_homo_sapiens_opg_2026_03.fasta`
  - `Universal_Contaminants.fasta`
- **Compatible search settings**:
  - Enzyme: Trypsin (K*, R*)
  - Missed cleavages: 1
  - Variable mods: Met oxidation (UniMod:35)
  - Fixed mods: Carbamidomethylation (UniMod:4)
- **Predicted spectral library**: Yes — `step1.predicted.speclib` available (skips Step 1)
- **Created**: 2026-03-01 by Brett
- **Notes**: "Standard human database for routine DIA searches"

### Actions
- **[Use This Database]** — Sets the FASTA path(s) for the current search. If a predicted speclib exists, offers to use it too (skip Step 1).
- **[Add to Library]** — Available after downloading from UniProt or uploading. Registers the current FASTA + metadata in the catalog.
- **[Delete]** — Removes entry from catalog (optionally deletes files too). Confirm dialog.
- **[Edit Notes]** — Update description/notes for an entry.

## Storage

### Shared volume location
```
/Volumes/proteomics-grp/dia-nn/fasta_library/
├── catalog.rds                          # R list of database entries
├── human_opg_2026_03/
│   ├── UP000005640_homo_sapiens_opg_2026_03.fasta
│   ├── Universal_Contaminants.fasta
│   └── metadata.json                    # Redundant copy for non-R tools
├── human_full_2026_02/
│   ├── UP000005640_homo_sapiens_full_2026_02.fasta
│   ├── Universal_Contaminants.fasta
│   └── metadata.json
└── mouse_opg_2026_02/
    ├── ...
```

On the HPC this maps to:
```
/quobyte/proteomics-grp/dia-nn/fasta_library/
```

### Catalog entry schema
```r
list(
  id = "uuid-string",
  name = "Human OPG + Contam",
  organism = "Homo sapiens",
  organism_common = "Human",
  proteome_id = "UP000005640",            # UniProt proteome ID (if applicable)
  content_type = "one_per_gene",          # one_per_gene | reviewed | full | full_isoforms | custom
  protein_count = 20598L,
  file_size_bytes = 12300000L,
  contaminant_library = "Universal",       # NULL if none
  contaminant_count = 246L,
  custom_sequences = NULL,                 # Character: pasted FASTA text, or NULL
  custom_sequence_count = 0L,
  fasta_files = c(                         # Paths relative to this entry's directory
    "UP000005640_homo_sapiens_opg_2026_03.fasta",
    "Universal_Contaminants.fasta"
  ),
  fasta_dir = "human_opg_2026_03",        # Subdirectory name

  # Search compatibility info
  search_settings = list(
    enzyme = "K*,R*",
    missed_cleavages = 1L,
    var_mods = "UniMod:35 (Met oxidation)",
    fixed_mods = "UniMod:4 (Carbamidomethylation)",
    min_pep_len = 7L,
    max_pep_len = 30L,
    min_pr_mz = 300,
    max_pr_mz = 1200
  ),

  # Predicted library link (ties into speclib cache)
  speclib_path = "/quobyte/.../step1.predicted.speclib",  # NULL if not available
  speclib_search_mode = "libfree",

  # Metadata
  created_at = "2026-03-01 10:30:00 PST",
  created_by = "Brett",
  notes = "Standard human database for routine DIA searches",

  # Remote paths (HPC)
  remote_dir = "/quobyte/proteomics-grp/dia-nn/fasta_library/human_opg_2026_03"
)
```

## Integration Points

### With UniProt download
After downloading from UniProt, a toast/button offers: "Add to shared library?" This pre-fills all metadata from the download.

### With speclib cache
When a search completes and the speclib is cached, the library entry gets updated with the `speclib_path`. Future users selecting this database see "Predicted library available — Step 1 will be skipped."

### With contaminant libraries
The bundled contaminant FASTA files (Universal, Cell Culture, etc.) are automatically included in the catalog entry when the user adds a database that included contaminants.

### With custom FASTA sequences
If the user pasted custom sequences (e.g., iRT standards), those are stored in the entry and appended to the FASTA at search time.

## UI Considerations

- Modal should be searchable/filterable (organism, created by)
- Show a badge or icon when a predicted speclib is available
- "Last used" column to surface frequently-used databases
- Graceful fallback if the shared volume isn't mounted — show message, fall back to UniProt/upload
- Read-only for non-creators? Or fully shared (anyone can delete)?

## Expiration Policy (6-Month Freshness)

Databases older than 6 months cannot be used for new searches. UniProt updates regularly (monthly releases), and stale databases miss new annotations, gene names, and sequence corrections.

### Behavior
- **< 6 months old**: Normal — green indicator, fully usable
- **5-6 months old**: Warning — yellow/amber indicator, "This database expires in X days. Consider refreshing."
- **> 6 months old**: Expired — red indicator, **[Use This Database] button disabled**. Shows: "This database is over 6 months old. Please download a fresh version from UniProt."
- **Expired + has predicted speclib**: The speclib is also invalidated. A fresh FASTA means a fresh prediction.

### Age calculation
Based on `created_at` in the catalog entry (when the FASTA was downloaded from UniProt), not when it was added to the library.

### Refresh workflow
When viewing an expired entry, offer a **[Refresh Database]** button that:
1. Re-downloads the same proteome/content_type from UniProt
2. Creates a new catalog entry with the updated FASTA
3. Archives the old entry (mark as `expired = TRUE`, keep files for reproducibility)
4. Clears the linked speclib (new FASTA = new prediction needed)

### Table display
| Name | Organism | Proteins | Age | Status |
|------|----------|----------|-----|--------|
| Human OPG + Contam | Homo sapiens | 20,598 | 2 months | Fresh |
| Mouse OPG + Contam | Mus musculus | 17,142 | 5 months | Expiring soon |
| Human Full + Contam | Homo sapiens | 83,526 | 8 months | Expired |

### Configuration
The 6-month default is hardcoded but could be made configurable via `DELIMP_FASTA_MAX_AGE_DAYS` env var for labs with different policies.

## Fallback / Edge Cases

- **Volume not mounted**: Show "Shared library unavailable — volume not mounted" and disable the option
- **Catalog file locked/corrupted**: Fall back to empty catalog, log warning
- **FASTA files deleted from disk but catalog entry remains**: Verify file exists on select, mark as "missing" with option to remove entry
- **Concurrent writes**: Use file locking (`flock` or R `lockfile`) when writing catalog.rds
- **Speclib path stale**: Verify speclib exists before offering to skip Step 1

## Future Extensions

- **Version tracking**: When UniProt updates, detect that the same proteome has a newer version and flag it
- **Auto-refresh**: Periodically scan the library directory for new entries added by other tools
- **Search history link**: Show which searches used each database
- **Export catalog**: Generate a summary table for lab documentation
