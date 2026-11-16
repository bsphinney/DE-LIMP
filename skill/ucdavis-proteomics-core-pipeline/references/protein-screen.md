# The protein screen ("show me the data for <protein>"), part 1

Design: `docs/specs/protein-data-screen-design.md` (sections 1, 3, 8, 10-12). Step 8g builds the data
behind the page; "Show me the data for <protein>" below builds the page. Where each protein and gene
fact comes from: `protein-card-sources.md`.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 19: 8g. Screen data for "show me the data" (limpa analyses)
- 67: Which analyses part 1 covers
- 98: Show me the data for <protein>
- 143: What the folder holds
- 155: The protein record
- 168: How the page is built
- 180: What is not in part 1, and why
- 190: The CPU-family rule
<!-- /contents -->

## 8g. Screen data for "show me the data" (limpa analyses)

After the final DE (and 8f), build the data behind the per-protein page, on a compute node and on the
DE's CPU family: limpa's refit of each protein must reproduce the delivered value to 1e-6, and on a
different family limpa's hyperparameters differ, so no refit describes the analysis and every protein's
shares are withheld ("The CPU-family rule", below). Read `compute.cpu_family` in
`output/tables/de_provenance.json`. A login shell on HIVE does not put the skill's R on PATH (the setup
step sources `activate.sh`), so the job sources it and names the script by its full path under
`~/proteomics-pipeline`.

When the family is recorded (the normal case), pin it with sbatch; the blocking job further down cannot.
Write the job file as step 8 does (`<S>` = the session; `session.py init` already made `<S>/logs`):
```
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/on_compute.py --header <S>/logs --name screen_data --time 120 --mem 16G --cpus 8 > <S>/logs/screen_data.sbatch'
bash scripts/hive_exec.sh 'cat >> <S>/logs/screen_data.sbatch' <<'EOF'
source ~/.proteomics-pipeline/activate.sh
Rscript ~/proteomics-pipeline/scripts/make_screen_data.R --session <S> --cpus "$SLURM_CPUS_PER_TASK"
EOF
```
Lint it, then submit with the family, both on HIVE (`U` is the `ucd-hpc` scripts folder, as in step 8).
Command-line options add to the header, whose log path already has `%j` (an `#SBATCH` line appended
after `set -euo pipefail` would be ignored):
```
bash scripts/hive_exec.sh 'bash -s -- <S>/logs/screen_data.sbatch' < "$U"/lint-jobscript.sh
bash scripts/hive_exec.sh 'sbatch --constraint=<that family> <S>/logs/screen_data.sbatch'
```
Only when no family is recorded, run it as one blocking job; `manifest.json` `compute.de` and
`compute.here` then say whether the two matched:
```
bash scripts/hive_exec.sh --job --time 120 --mem 16G --cpus 8 'source ~/.proteomics-pipeline/activate.sh && Rscript ~/proteomics-pipeline/scripts/make_screen_data.R --session <S> --cpus 8'
```
Log each command to `commands.log` like every other. It writes `output/screen_data/` (`manifest.json`,
`index.tsv`, gzipped `parts/`). On PXD063416 it took 2-3 minutes for 8,500-10,300 proteins and the
folder was 36-52 MB. A protein whose refit missed is withheld with the reason shown on its page; relay
the count (`manifest.json` `shares.withheld`). Not part of the report, and not yet delivered with it;
step 12 names the folder (outputs.md, step 12).

The run records its state in `output/screen_data_status.json` (`running`, `built`, `refused` with the
route and reason, or `failed` with one line), and `session.py status` shows it as an optional note
("protein screen data: not built (optional)" when 8g was not run). A `running` record that no job is
still working on means the job stopped part-way (killed, out of memory, timed out): say so and rerun.
It is never a finish item.

**Exit 2 = refused.** The reason is on stderr and in `manifest.json` (`supported: false`,
`unsupported_reason`, and `analysis`: the route as the analysis describes itself). Say so in one line,
with its reason; do not work round it. It also refuses a contrast that is not one group against one other
(`(B+C)/2-A`, an interaction), because the page names a numerator and a denominator group. What each route gets is in "Which analyses part 1 covers".

## Which analyses part 1 covers

One decision, `screen_support()` in `make_screen_data.R`, read from what the analysis says it is (repo
rule 1): `de_provenance.json` `pipeline_id`, named by its own `display_label` and `rollup_method`, and
whether the report is DIA-NN's (an adapter's report declares its own quantity in the parquet metadata
`delimp.quantity.*`, read by `build_maxlfq.R` `declared_quantity()`; DIA-NN's declares nothing). Every
block of part 1 reads limpa's own objects (`raw_data`, `y_protein`, `dpc_fit`, the vooma weights) and
DIA-NN's precursor rows.

| route (`pipeline_id`) | engine and quantity | part 1 |
|---|---|---|
| `dpc` | DIA-NN report (native, or FragPipe's DIA route), limpa DPC-Quant + limma, plain group design | the whole page: DE, protein values, precursors and report rows, shares, detection chances, weights, run statistics |
| `dpc` over a report that declares its own quantity | a Sage / FragPipe DDA / AlphaDIA / Radiant adapter's report, then limpa | refused: the precursor blocks would describe DIA-NN rows that are not there |
| `dpc` with `--block` or a numeric covariate | DIA-NN, limpa, blocked or covariate design | refused (the reduced screen, part 5) |
| `dpc`, every protein one precursor | limpa's `dpcQuantByRow` | refused (no protein model to show) |
| `maxlfq` | DIA-NN `PG.MaxLFQ` + limma | refused (part 5) |
| `engine_protein` | an engine's own protein quantity + limma: FragPipe IonQuant MaxLFQ, the skill's directLFQ, AlphaDIA, Radiant | refused (part 5) |
| `peptide_max` | Sage: highest peptide per protein and run + limma | refused (part 5) |
| SILAC (`run_silac.R`) | ratios | refused: its record is not `run_de.R`'s limpa record (part 5) |

A refusal never leaves a half page: `make_protein_screen.py` reads `supported: false` and stops with the
reason and the route ("It used MaxLFQ + limma (pipeline_id maxlfq)"). Repeat that; never describe the
protein's data another way.

The report read is always the one the DE read: `de_provenance.json` `input` when it lies under
`output/` and still has the DE's `input_sha256`; for a session moved or copied since, the file at the same
relative path under this session's `output/` (the tail after `/output/`), if it has that checksum. There is
no search through the folders. Its
`report.stats.tsv` is read from the same folder. The manifest's `files` names both, relative to
`output/`, and the page shows those names. A record without the path or hash, a report that changed, or none at either place is a refusal.

## Show me the data for <protein>

When the user asks to see the data behind a protein ("show me BASP1", "why is this one up?", "show
me the top 10 changed proteins"), build the page. It is the answer.
```
python3 scripts/make_protein_screen.py --screen-data <session>/output/screen_data --protein BASP1 --out BASP1.html
python3 scripts/make_protein_screen.py --screen-data <session>/output/screen_data --top 10 --out-dir top10/
python3 scripts/make_protein_screen.py --screen-data <session>/output/screen_data --top 10 --direction up --out-dir top10up/
```
It runs on the user's own computer (Python only, no R); then open the page for them.
**For a session on HIVE, build the page on a compute node, never on the login node:** it reads and
hashes the search report (gigabytes), which golden rule 3 and the `ucd-hpc` rules keep off it. Run it
as one blocking job, as step 8g's own fallback does, write the page into the session, then fetch it:
```
bash scripts/hive_exec.sh --job --time 60 --mem 4G 'python3 ~/proteomics-pipeline/scripts/make_protein_screen.py --screen-data <S>/output/screen_data --protein BASP1 --out <S>/output/protein_screens/BASP1.html'
bash scripts/hive_exec.sh --get '<S>/output/protein_screens/BASP1.html' ./
```
`--top N` takes `--out-dir <S>/output/protein_screens/top10` instead; `--get` that folder. In
`hive_local` the same command goes through `python3 scripts/on_compute.py --time 60 --mem 4G -- ...`
(the Flow section's split). A source the compute node cannot reach is marked "not reached" on the
page; say so. Like `screen_data/`, the pages are not delivered with the results. If there is no
`screen_data/` folder, run step 8g for a session you can reach; otherwise say it has not been built
and who can build it (the Core). If 8g refused, say why in one line (its `unsupported_reason`). A gene
that matches several proteins is listed in full: ask which.

`--top N` takes the proteins the comparison calls changed (`adj.P.Val` alone, no fold-change cut),
best `P.Value` first, at most 100. `--direction up|down|both` (default both) keeps the up or down
ones: "top 10 up" is `--top 10 --direction up`. If fewer are called than asked, the page says how
many; repeat that, never pad. `--contrast NAME` picks the comparison (the first by default).
**Never pick, rank or list proteins yourself, from memory or by reading the tables.** A "top N" answer
is the builder's list, quoted as it stands, including its "only K" message.

The page checks every file against the checksums in `manifest.json` and refuses an analysis that has
changed since. Protein and gene facts come from UniProt, NCBI, HGNC, Ensembl and the Human Protein
Atlas, looked up when the page is made, each shown as the source wrote it with its release, date and
link; `--offline` looks nothing up and the page says so. The page makes no network request.

**Never describe a protein's biology, and never write a number about it, from memory.** Do not
add "this protein is known for ...", a function, a pathway, a disease link or a fold change. If the
user wants a summary, quote only what is on the page, attributed to its source on the page. A fact
the page lacks is a fact you do not have: say the source had nothing and offer to look it up with
`scripts/protein_facts.py`. Wording on the page: a peptide or precursor is *found* or *not found* in a
sample; only a protein value is *inferred*. XICs, check routes, FRAN facts, tracing the report's text
and delivery with the results are not built yet: say so if asked.

## What the folder holds

`output/screen_data/`, written by `scripts/make_screen_data.R --session <dir> [--out <dir>] [--cpus N] [--proteins PG1,PG2]`:

- `manifest.json`: `schema: "screen_data/1"`, `written_by`, `written_at`, `supported`, `unsupported_reason`, `partial` (true with `--proteins`), `analysis` (`pipeline_id`, `display_label`, `normalisation`, `quantities`, `q_columns`, `q_cutoffs`, `filters_applied`, `rollup_method`, `de_engine`, `significance_rule`, copied from `de_provenance.json`; a refusal carries it too), `files` (`session_rds`, `report`, `stats`: where the page's files are, relative to `output/`), `samples`, `groups` (one per sample), `contrasts` (`{name, numerator, denominator, table, n_tested, n_sig, adjp}`), `dpc`, `prior`, `n_precursors`, `run_stats`, `sources` (`{path under output/: sha256}` of every file read), `parts` (`{path: sha256}`), `index`, `skipped` (`{block, reason}`), `shares` (`method: "derivative"`, `nudge_check`, `limpa_version`, `withheld`), `detection`, `elapsed_s`, `compute` (`de`: the DE's CPU family, `here`: this node's).
- `index.tsv`: `Protein.Group`, `Genes`, `Protein.Names`, `part`, `n_precursors`, then `P.Value__<contrast>`, `adj.P.Val__<contrast>`, `logFC__<contrast>` per contrast (empty = not tested).
- `parts/part-NNNN.json.gz`: `{Protein.Group: record}`, at most 250 proteins per part.

`--proteins` writes a partial folder into its own `--out`; it refuses `<session>/output/screen_data`, which holds the full screen data. A refusal for wrong arguments leaves the folder's `manifest.json` as it was. The builder's `--top` refuses a partial folder, and a protein it cannot find there is reported as possibly in the full folder. Screen data written before the manifest had `files` is read from `search/` (the builder's defaults).

All sample-indexed arrays follow `manifest.samples`. Numbers: log2 values to 4 decimals, quantities to 6 significant digits, q-values to 3. NA is null.

## The protein record

`{pg, genes, names, n_precursors, sigma, de, values, precursors, rows, shares, detection, weights}`

- `n_precursors`: limpa's count of the protein's precursors, `y_protein$genes$NPrec` (limpa 1.4 and later; `NPeptides` before, the same count). Null, with the reason in `manifest.skipped`, when the analysis has neither.
- `de`: `{<contrast>: {tested, logFC, CI.L, CI.R, t, P.Value, adj.P.Val, df_total, evidence, detected: {<group>: Detected_<group> as written}}}`.
- `values`: `{protein: [...], se: [...], n_obs: [...]}`, the delivered protein values.
- `precursors`: `[{id, log2: [...]}]`, null for a gap.
- `rows`: one per DIA-NN report row: `Run, sample, Precursor.Id, Protein.Ids, Proteotypic, Precursor.Quantity, Normalisation.Factor, Precursor.Normalised, Quantity.Quality, Empirical.Quality, Ms1.Area, RT`, each q-value column, `PG.MaxLFQ, PG.MaxLFQ.Quality`, `failed` (the q columns over their cutoff), `contaminant`, `used`.
- `shares`: `{own: [[precursor][sample]], other: [sample], total: [sample], refit_max_diff}`, or `{withheld: <reason>}`. A share is the derivative of a protein value with respect to one found precursor value, taken at limpa's optimum from limpa's own Hessian (design 9b).
- `detection`: `{pf_fit, pf_obs}`, the chance limpa's fit gives each precursor of being found in each sample; withheld with the shares, for the same reason.
- `weights`: `{w, share_in_group, df_residual, df_prior}`, limma's precision weights of the delivered fit, never refitted; null, with the reason in `manifest.skipped`, when `set_test_inputs.rds` is missing.

## How the page is built

`scripts/make_protein_screen.py` (Python standard library only; runs on a laptop). For a session on HIVE it runs on a compute node, never on the login node, because it hashes the search report ("Show me the data for <protein>", above).

1. Reads `manifest.json` and checks the sha256 of the parts, the index and every source file in `output/`; any change since the screen data was made stops it ("rebuild it with make_screen_data.R"). `supported: false` stops it with the reason and the route.
2. Finds the protein by gene, accession or Protein.Group. An ambiguous gene lists every match and picks none.
3. `--top N [--direction up|down|both]`: the proteins the comparison calls changed (`adj.P.Val` under the analysis's cutoff, no fold-change cut), ranked by `P.Value`, at most 100. Fewer called than asked: the page says how many.
4. Fetches protein and gene facts through `scripts/protein_facts.py` from UniProt, NCBI, HGNC, Ensembl and the Human Protein Atlas; each fact shown verbatim with source, release, fetch date and link. Sources, endpoints, licences and limits: `protein-card-sources.md`. `--offline` asks nothing. NCBI email and key come only from `NCBI_EMAIL` and `NCBI_API_KEY`. A group of several accessions is looked up by its first, and the card says so. An isoform accession (`P80723-2`) has no UniProt entry of its own: the canonical entry is shown, with a note saying so, and sequence coverage and peptide positions are left out, with the reason (they would be counted on the canonical sequence). Isoform sequences are not fetched.
5. Writes one HTML file that makes no network request. `--dump-data` also writes the page data as JSON.

Nothing on the page about a protein comes from the AI's memory (Brett, 2026-10-28).

## What is not in part 1, and why

Design section 12 builds the screen in parts, each working alone.

- XICs (peak traces) and library uniqueness: part 2. They need the raw-data excerpts, which `make_screen_data.R` does not read yet.
- Delivery with the results (`core_submission.py` shipping `screen_data/`): part 2. Until then the folder stays in the session.
- The check routes (MaxLFQ and Moschem) and the flags: part 3. The referee route's normalisation must be fixed first (design 10b).
- FRAN facts, the short answer, the picture, and tracing numbers in the report text: part 4.
- The reduced screen for other routes and designs: part 5. Until then `make_screen_data.R` refuses every route but limpa DPC-Quant over a DIA-NN report with a plain group design ("Which analyses part 1 covers").

## The CPU-family rule

limpa's refit has to reproduce every delivered protein value to 1e-6. The value differs between CPU families (reproducibility.md, "Exact DE numbers need the same kind of CPU (DPC-Quant)": DPC-Quant's numbers shift in about the third decimal between zen2 and zen4 nodes, |ΔlogFC| ≤ 0.0022 across 12 contrasts of one cohort, with `dpcQuant` the stage that differed). On another family limpa's hyperparameters differ, so no refit describes the analysis and `make_screen_data.R` withholds every protein's shares with that reason.

So when `compute.cpu_family` is recorded in `de_provenance.json`, the job is a job file from `on_compute.py --header`, submitted with `sbatch --constraint=<that family>` (step 8g, above); `hive_exec.sh --job` cannot pin a family, and is used only when no family is recorded, when the manifest's `compute.de` and `compute.here` say what happened. A protein whose refit misses is withheld with the reason, never shown. Heavy work stays on a compute node (golden rule 3); under `ucd-hpc` the job names its account, partition, time and memory.

Measured on PXD063416 (HIVE, zen2, 8 CPUs, 16G): 8,500-10,300 proteins took 106-155 s per analysis; screen data was 36-52 MB.
