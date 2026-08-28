# This is a sandbox, not a second product

Forked from `ucdavis-proteomics-core-pipeline` **v2.4.3**. It exists so new capability can
be built without touching the skill that runs real searches on real data — a skill that
now submits SLURM chains, stages results into FRAN, and is used for facility work.

## The rules that keep the fork from becoming a liability

1. **Never run a real analysis with this skill.** Its description tells Claude not to
   select it automatically, and it is not a fallback if the stable skill is unavailable.
   That is deliberate: a sandbox that quietly does production work is worse than no
   sandbox, because nobody knows which one produced a result.

2. **Edit only what is genuinely new here.** Everything else is a copy, and a copy edited
   in two places is the divergence this repo has been fixing all week (24 hardcoded
   contaminant checks, two copies of a chunk-boundary counter). If you find yourself
   fixing a *stable* bug in this tree, fix it in the stable skill and re-fork.

3. **Port files back, never merge trees.** When something is ready it moves to the stable
   skill as the specific new files plus their tests and doc sections. Merging this whole
   directory back would drag along whatever drifted.

4. **Re-fork rather than rebase.** When the stable skill moves, delete the unchanged files
   here and copy them again. Track what is genuinely new in the list below so that is cheap.

## New here, not in stable
_(nothing yet — this is the fork point)_

## Work in progress
Porting DE-LIMP's phosphosite capability into the skill:

| from DE-LIMP | what it does | notes |
|---|---|---|
| `extract_phosphosites()` | report.parquet → site-level matrix at a localization threshold | the threshold decides which sites exist; it must be a stated parameter, never a silent default |
| site-level DE | limma on the site matrix, contrast-wise | shares the DE contract with the protein path |
| KSEA | kinase activity from site fold-changes | needs the PSP + NetworKIN table; PSP is licence-encumbered, so the skill must not redistribute it — see DE-LIMP `R/helpers.R` `load_ksea_database()` |
| residue / motif analysis | S/T/Y distribution, sequence motifs around the site | needs the FASTA to fetch flanking sequence |
| the visualisations | volcano, heatmap, completeness, expression grid | DE-LIMP renders these in Shiny; the skill needs them as static figures for the HTML report |

**The hard part is not the plots.** DE-LIMP's versions are Shiny reactives that read
`values$...`; the skill needs standalone scripts that take files in and write figures out,
with every parameter stated. Expect to rewrite rather than copy.
