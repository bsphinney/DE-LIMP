# Validated proteomics workflows

This folder is the source of truth for the **proteomics-pipeline** skill. When a
user says "analyze my proteomics data," the skill fetches `index.json` from here
(raw GitHub, no auth), matches their data to one of these validated workflows,
downloads the pinned software, and runs the search + differential expression
exactly as validated.

You maintain the workflows; the skill never has to be re-shipped. Add or revise
a workflow by committing to this folder.

```
workflows/
├── index.json                  # compiled list the skill reads first (DO NOT hand-edit)
├── build_index.py              # regenerates index.json from the workflow.yaml files
├── schema/workflow.schema.json # JSON Schema documenting/validating a manifest
└── <workflow_id>/
    ├── workflow.yaml           # human-authored manifest (the source of truth)
    ├── VALIDATION.md           # provenance: what data, when, by whom, results
    └── <engine params file>    # OPTIONAL: a validated config to override estimation
```

**Search parameters are estimated from the data type by default** — you do NOT
hand-maintain a `.cfg`/`.json` per workflow. The skill's `estimate_params.py`
derives mass tolerances from the detected instrument (DIA-NN's known-good table)
and the DIA/DDA window mode from acquisition. A workflow only needs a params file
if you want to lock a fully validated SOP config (see "Adding a workflow" below).

## How a workflow is chosen (matching)

The skill keys on **acquisition + instrument + organism**:

- **Acquisition (DDA/DIA)** and **instrument** are detected from the raw files.
- **Organism** cannot be detected — it defines the FASTA — so the skill asks the
  user. That question is folded into the same step where it collects the
  experimental design.

Selection logic (implemented in the skill's `fetch_workflows.py`):

1. **Hard filters** — a workflow is a candidate only if its `acquisition`
   matches the detected acquisition AND its `organism_taxid` matches the
   user's organism. Running a human workflow on yeast data is never right.
2. **Score** the candidates by instrument: exact model (alias-aware) beats a
   family/partial match beats no instrument info.
3. **Auto-pick** the top score and **show it to the user to confirm**, with the
   other candidates offered as fallback. (This is the "auto-match, then confirm"
   behavior.) Ties or a weak top score → present the menu instead.

The match-relevant fields are denormalized into `index.json` so the skill needs
only one fetch and **no YAML parser at runtime**. The engine param files are in
each engine's native format, so they are fetched and handed straight to the tool.

## Adding a validated workflow

1. `cp -r diann_astral_dia_human workflows/<your_new_id>`
2. Edit `workflow.yaml` (see the schema and the worked examples). Pin the engine
   version you validated against.
3. **Parameters:** leave `search.estimate_params: true` to derive them from the
   data type (the default). To nudge specific values, set `search.param_overrides`
   (e.g. `{"--mass-acc": 8}`) — merged on top of the estimate and tagged as a
   user-override. To lock a complete validated config, drop the file in the folder
   and set `search.params_file:` to it (estimation is then skipped).
4. Write `VALIDATION.md` — the dataset you tested on, date, expected result.
5. Regenerate the index: `python3 build_index.py` (or just let the GitHub
   Action do it on push — see `.github/workflows/build-index.yml`).
6. Commit. The skill sees it on the next run.

## Manifest fields (workflow.yaml)

| field | meaning |
|---|---|
| `id` | unique slug; must equal the folder name |
| `name` | human label shown in the confirm prompt |
| `match.acquisition` | `DDA` or `DIA` (hard filter) |
| `match.instruments` | validated instrument model(s) |
| `match.instrument_aliases` | lowercase substrings to fuzzy-match detected names |
| `match.organism` / `match.organism_taxid` | hard filter; taxid is authoritative |
| `engine.name` / `engine.version` | `diann`\|`sage`\|`fragpipe`, pinned version the acquire step installs |
| `fasta.uniprot_proteome` | UniProt proteome ID the skill downloads (e.g. `UP000005640`) |
| `fasta.add_contaminants` | append the universal contaminants FASTA |
| `search.estimate_params` | derive params from data type (default true) |
| `search.var_mods` | optional extra variable mods for estimation, e.g. `ox` |
| `search.param_overrides` | validated values forced on top of the estimate |
| `search.params_file` | OPTIONAL validated config; if set, used verbatim (skips estimation) |
| `de.method` | `dpc` (limpa) or `maxlfq` |
| `de.q_cutoff` / `de.logfc` / `de.adjp` | DE thresholds |
| `validated.*` | provenance, for trust and reproducibility |

## Design note

Versions are pinned per workflow, not in the skill. That is deliberate: it keeps
results reproducible and under your control, and lets you validate a new engine
version in one folder without disturbing the others.
