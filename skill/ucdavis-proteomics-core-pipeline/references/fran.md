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
only reports it. (Same group `check_access.sh` reports as `core_member`.)

### QC runs are never handed over

A QC run (a HeLa series watching an instrument) is not a customer search, and FRAN keeps QC
out of the corpus, as it does STAN. FRAN's scanner can only recognise QC by path, and a drop
entry has no QC path, so the skill decides. There is **one** definition,
`fran_deposit.is_qc_run(out, session)`, used by `check`, `stage` and `backfill`. Its twin in
FRAN is `ingest/find_uningested.py` `policy_exclusion` / `qc_reason` / `QC_NAME_RE`: change one,
change the other. The precedence is FRAN's; the first match wins:

1. **Explicit QC:** `stage --qc`, or `"qc": true` in the session's `session.json` (or
   `input/session.json`, `input/wf/workflow.manifest.json`).
2. **FRAN's DEFAULT_EXCLUDES trees:** `/quobyte/proteomics-grp/STAN/`, `…/hela_qcs/`,
   `…/brett/v1_smoke`, `…/brett/glendon/` and `/Data/lab/ToFEvoQC/`. These win **even over
   `--not-qc`**: FRAN refuses anything there, so staging it would only queue a refusal.
3. **Explicit not-QC:** `--not-qc`, or `"qc": false` in session metadata.
4. **FRAN's name rule** `(?i)(?<![a-z0-9])qc(?![a-z])`, applied to the search name
   (`--name`), the session's name (its README title and folder) and the last three
   components of the out dir's **real** path. FRAN judges `output_dir`, which is the realpath,
   so a symlinked search dir must not escape the rule here. DEFAULT_EXCLUDES is tested on both
   spellings, with `/Volumes/proteomics-grp` mapped to `/quobyte/proteomics-grp` as FRAN does. It catches `chkLUppm_HeLa50_2026 Lumos QC`, `QC_run_01`,
   `hela_qc_2` and `Exploris QC2`. It keeps `HeLa_digest_timecourse`, `aqc_buffer_study`,
   `QCM_study` and `Plasma_liver2`, the same pinned vectors as FRAN. "HeLa" alone is **not**
   QC.

A QC run gets reason `qc_run`, with a `why` in FRAN's wording, e.g. `QC run: excluded by policy
(search_name 'chkLUppm_HeLa50_2026 Lumos QC' matches QC_NAME_RE)`. The receipt records it. `--not-qc` corrects a
false positive, and the refusal says so. Every staged manifest carries `"qc": false` plus
`"qc_rule": "<why>"` (`"user override"` for `--not-qc`), so FRAN's ingester sees that a
decision was made.

**Decide at generation, not afterwards.** On 2026-09-23 the QC session was called
`2026-09-23_chkLUppm_HeLa50_2026`, with no QC token and no `conditions.csv`. "Lumos QC" only
arrived later, with the agent's `--name`. But the job-end hook stages the moment the search
ends, and FRAN's cron (every 4 h) can ingest before a later `stage` says otherwise. So the name
and the QC decision are baked into the hook when the job is written. `run_search.py` /
`diann_parallel.py` / `radiant_parallel.py` take `--fran-name "<descriptive name>"` and
`--qc` / `--not-qc`, and the hook's argv comes from the one helper, `fran_deposit.stage_argv()`.
A QC run named at generation never reaches the drop dir. An explicit `--qc`/`--not-qc` is
recorded in the receipt, and a later stage without a flag honours it rather than overturning it.

**Withdrawal is the backstop.** When a later `stage` still finds that a staged search is QC, it
**withdraws** it: the entry's manifest gets `"qc": true, "exclude": true`, which FRAN's ingester
skips at ingest time. Nothing in the shared drop dir is deleted, and `verify` then reports
`qc_excluded`.

- **A withdrawal sticks.** A `qc_run` receipt, or a staged manifest that says `qc: true`, means
  QC for every later `stage` until someone passes `--not-qc` explicitly. A plain `stage --out X`
  used to re-stage a withdrawn run with `qc: false`, which FRAN honours.
- **A withdrawal that fails says so.** For example, the entry was staged by another account
  without group write. The result then starts `withdraw FAILED: <why>`, names the entry that is
  still staged, and never claims the run was kept out.
- **`--no-fran` should run `stage --skip`,** so the opt-out is recorded and `backfill` never
  picks the search up later.

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

# is FRAN's cron actually taking what is staged?  (reads logs; login-node safe, ~5 s)
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/fran_deposit.py health'
```

`verify` returns one of five states:

| state | meaning | say |
|---|---|---|
| `ingested` | the corpus has it (database, cron marker, or the cron's log: `OK` or `SKIPPED-DUPLICATE`) | "in FRAN" |
| `staged_pending_cron` | handed over, not reached yet. **A success**, the normal state for hours after a run | "handed to FRAN; it ingests on the next pass" — and if `cron.verdict` is `stuck`/`not_running`, add that FRAN's cron is stuck, which is FRAN-side |
| `ingest_failed` | the cron tried this entry and it failed; `detail` has the reason from its log | the reason, in one line. FRAN-side: do **not** re-stage |
| `qc_excluded` | a QC run whose entry is marked `qc: true` | "a QC run, kept out of FRAN" |
| `not_staged` | no entry | why `check` refused. `verify` writes **no** receipt then: an invented `staged` receipt used to make `stage` say `already_staged` and `backfill` skip the search for ever |

Without a corpus token (the usual case), `verify` answers from the cron's own logs rather than
the database. `broken_links` is the one entry problem to act on: an entry whose targets have
moved looks staged but the cron will skip it in silence. Re-run `stage --force`.

## Staged is not ingested — `health`

Staging puts the search where FRAN's cron looks. Whether the cron gets to it is a separate
fact, and for a week it was false while everything looked fine. After 2026-09-17 13:55 the
cron ran every 4 h, each run a COMPLETED SLURM job, and ingested nothing in 41 runs in a row.
The recent ones all ended `0 ingested, 3 duplicate-skipped, 2 failed, ~181 still queued`. The
same five `FRAN_reports` exports sort first, fail or duplicate, are never marked done, and are
picked again next run, so nothing behind them moves. None of the skill's drop entries appeared
in a single log. (Two of those searches were in FRAN anyway: FRAN's queue ingested them by
their real paths on 2026-09-08.) `health` is how this becomes visible from the skill's side.

It is read-only and needs no credential. It reads what the cron itself writes
(`/quobyte/proteomics-grp/de-limp/fran_refresh/logs/auto_ingest_<jobid>.out` and
`auto_ingest_submit.log`), lists the drop dir, and compares FRAN's ingest code on HIVE with
GitHub `main`:

| part | answers | verdicts |
|---|---|---|
| `progress` | last run, last run that ingested anything, consecutive runs without an ingest, queue size, whether the cron is still submitting | `healthy` · `stuck` (≥ 3 runs in a row ingested nothing while work was queued; a run that died counts) · `not_running` (no run, or no submission, for > 12 h) · `unknown` |
| `incoming` | each drop entry: age, who staged it and when, broken links, what the logs say (`ingested` / `failed` / `never_reached`), and whether it is a QC run not yet marked (`qc_unmarked`) | `starved` when an entry has gone unreached for > 48 h, even if the cron is ingesting *other* searches |
| `ingest_code` | each ingest file's md5 on HIVE vs `main`: `current`, `stale` (matches an older commit: *which one, from when*), `local_modification` (matches none of that file's recent commits), `missing`, `not_on_main`, `unknown` | `stale` if any file is stale or missing, or a file on FRAN's own refuse list (`corpus_ingest.py`, `spectronaut_to_corpus.py`, `diann_to_corpus.py`, `versions.py`) differs at all · `modified` · `current` |

Overall `verdict`: `healthy` · `stuck` · `not_running` · `stale_code` · `unknown`, plus a one-line
`summary`. The ingest-code check fetches file content from raw.githubusercontent.com (no API
quota) and uses GitHub's commits API only for a file that differs. That API allows 60 requests
an hour per address, shared by everyone on a login node, so answers are cached in
`~/.cache/ucdavis-proteomics-core-pipeline/fran_health.json` and a rate limit is remembered
until it resets. No network gives `unknown`, never an error. FRAN's own staleness guard
(`publish_manifest.py`, content md5s in PG Farm) needs a database credential, so the skill does
not read it.

**`stage` never runs this check.** It runs inside every search job, on a compute node, within
the job's time limit, so it must not reach GitHub or the database, or scan 170 logs on a network
mount. Instead `health` writes its verdict to **`/quobyte/proteomics-grp/fran/ingest_health.json`**
(`checked_at` in UTC, `verdict`, a one-line `summary`; group-writable, replaced atomically). It
sits in the drop dir's parent, never inside `incoming/`, where nothing but drop entries belongs.
`FRAN_HEALTH_FILE` moves it.
`stage` only reads that file, and gives up on the read after 2 s. What it adds to its JSON:

| status file | `stage` adds |
|---|---|
| verdict `stuck` / `not_running` / `stale_code`, checked < 12 h ago | `fran_health` + `health_warning` (the summary); the same line on stderr |
| `healthy` or `unknown`, checked < 12 h ago | `fran_health` only |
| older than 12 h | `fran_health: {verdict: "unknown", summary: "… (last checked <when>)"}`, and **no** warning: old news is not bad news |
| missing, unreadable, or a read that does not finish | nothing |

The search **is** staged either way. Health never changes stage's exit status, and stdout
stays one JSON object. `FRAN_HEALTH=off` skips even the read. The file is only as fresh as
the last `health` run: the agent runs it at step 7c, and the cron line below keeps it current.

The GitHub check sends no credential. It tolerates the rate limit, and it is bounded: all of
it runs in daemon threads under a 10 s deadline. That also covers a DNS lookup that hangs,
which a socket timeout does not.

**`health --alert` pages about one thing only:** what FRAN's own runner cannot see about
itself. That is FRAN's ingest code on HIVE being `stale` or `modified` against GitHub `main`,
or missing a file the runner imports. A stuck cron, entries that are never reached, and
manifest-vs-search FASTA mismatches are alerted by FRAN's runner itself (FRAN branch
`fix/auto-ingest-starvation`, 78374dc + 2163070). Paging on them here as well would be a
duplicate page, so they go only into the verdict file and the printed report (`alert_reason`
says what would page). The alert goes through `scripts/notify_slack.py` when that module is
installed (a no-op otherwise). The same alert is not re-posted within 24 h. The webhook is
that module's business and is never read or printed here.

**Suggested schedule** (in the deploy bundle, for Brett's OK; not installed). Add to brettsp's
crontab, half an hour after FRAN's own `23 */4` submit, on the login node. It reads logs and
a few small files, plus at most 10 s of GitHub:

```cron
# FRAN ingest health for the skill: refreshes /quobyte/proteomics-grp/fran/ingest_health.json
# (read by every search job's stage) and pages only if HIVE's ingest code drifts from GitHub main.
53 */4 * * * flock -n /tmp/fran_skill_health.lock bash -lc "python3 $HOME/proteomics-pipeline/scripts/fran_deposit.py health --alert > /dev/null" >> /quobyte/proteomics-grp/de-limp/fran_refresh/logs/skill_health_cron.log 2>&1
```

What to do with an unhealthy answer: **nothing on the search**. Tell the user in one line that
the search is safely staged and that FRAN's ingest is stuck/stale on FRAN's side (quote the
summary), and carry on. Do not re-stage, and do not touch FRAN's code, its HIVE copy, or the
database.

## Backfill — searches that were never staged

Staging became automatic partway through the skill's life, and a session that ended early
never reached step 7c, so some Core searches on HIVE were never handed over. FRAN's cron scans
only `incoming/` and `FRAN_reports/`, never the service trees, so it will not find them by
itself.

```bash
# on HIVE: write the job (dry run), then submit it -- never walk NFS on a login node
python3 ~/proteomics-pipeline/scripts/fran_deposit.py backfill --sbatch
sbatch ~/fran_backfill/fran_backfill_<stamp>.sbatch
# read ~/fran_backfill/fran_backfill_<jobid>.json; then, if the list is right:
python3 ~/proteomics-pipeline/scripts/fran_deposit.py backfill --sbatch --apply
```

- **Where it looks:** `/quobyte/proteomics-grp/SERVICE`, `/nfs/lssc0/flinders/proteomics/Data/lab/service`,
  and `~/proteomics-pipeline` of each non-teaching `proteomics-grp` member (`--no-homes` to skip).
  `--roots` replaces the trees; `--list <file>` checks named out dirs instead (up to 10 are fine
  on a login node). Checkpoint sessions (`.recovery.json`) found on the way add the search out dir
  they point at.
- **How it walks:** breadth-first, `--max-depth 9` (a Flinders session's `output/search` sits at
  depth 8), a shared `--time-budget` (default 1200 s, split across roots), no symlink ever
  followed, raw-data containers (`.d`, `.raw`, `.wiff`, ...) pruned, and the trees FRAN's own
  scanner excludes (STAN QC, QC watchers, smoke tests, scratch) skipped. A walk outside SLURM is
  refused unless `--allow-login-node`. A walk that runs out of time says `truncated`.
- **What counts as this skill's search**, from file names alone: `search_provenance.json`,
  `fran_deposit.json`, the 5-step chain's SLURM logs (`s1_libpred_<job>.log`,
  `s5_report_<job>.log`), `step3_fulcrum.sbatch`. `search_provenance.json` alone is not enough:
  in hive_remote mode `run_search.py` writes it on the laptop, and a real parallel search on HIVE
  had none. The DE-LIMP app writes the same `step*_*.sbatch` names but its jobs are
  `diann_<name>_<step>`, so its searches are `not_a_skill_search`.
- **Who is handed over:** a search inside `/quobyte/proteomics-grp/` or
  `/nfs/lssc0/flinders/proteomics/`, or owned by a non-teaching Core member. Anything else is
  `not_core_facility`. **Coursework never is.** A search run by, or owned by, a teaching account
  (`proteomics-class-NN`) is refused by `backfill` and by `stage` alike, through one rule
  (`_teaching_reason`). That matters because those accounts are in `proteomics-grp` and can
  write the drop dir, so their job-end hook could otherwise stage a class exercise.
- **A finished search, not just a report.** For FragPipe and Radiant a report exists before
  the search has finished. On HIVE, 16 of 16 FragPipe workdirs had `dia-quant-output/report.tsv`,
  including two cancelled runs. So `check` requires the engine's own marker. For FragPipe that
  is the newest workdir `log_*.txt` ending `ALL JOBS DONE IN <n> MINUTES` (all 11 finished runs
  had it; none of the 5 unfinished did). For Radiant it is `_SUCCESS` in `fulcrum-results/` (the
  Spark commit marker, on 20 of 20). A missing or cancelled marker is `search_incomplete`.
  Where there is no marker to read (no FragPipe log; a Radiant search with only
  `delimp_report.parquet`), the agent's own `stage` goes ahead, because it watched the job
  finish. `backfill` reports `needs_agent_check` instead, and `--apply` never stages it.
- **QC runs are listed on their own** under `excluded_qc_run`, never staged. The rule is the
  same `is_qc_run`, applied to the folder name and to any analysis name an earlier `stage`
  recorded. An `--qc`/`--not-qc` recorded by an earlier stage wins. A QC run that is already in
  the drop dir shows `would_withdraw`, and `--apply` marks its manifest `qc: true`.
- **Skipped, with the same reason codes as `check`,** plus two of backfill's own:
  `not_a_skill_search`, and `already_ingested` when the cron's logs show FRAN ingested the
  search by *any* route. FRAN's queue ingested some skill searches by their real paths with no
  receipt, and staging those again would only feed the duplicate guard. A recorded `opted_out`
  is honoured: `stage --skip` / `FRAN_DEPOSIT=off` now writes it into the receipt.
- **With `--apply`,** each eligible search goes through the same `stage` code. The search name
  is derived from the folder (for `<session>/output/search`, the session's name), and the
  manifest says so in `search_name_source`. The organism and database come only from a
  `<fasta>.meta.json` that can be tied to the search (below); with no sidecar they stay absent.

The report lists `would_stage` (engine, organism, name, report size and date), `skipped` with
reasons, the walk statistics, and the cron's current verdict. Staging into a stuck cron only
adds to its queue, so it is worth a look before `--apply`.

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
`radiant`. The skill's own file-sniff markers (`DETECT_MARKERS`, used when there is no
provenance) are FRAN's `ENGINE_MARKERS` (c4838fe) minus `search_provenance.json`.
`delimp_report.parquet` is a Radiant *report location* on both sides, but a *marker* on
neither.

## Pass the organism — it is the one thing only the skill knows

A DIA-NN `report.parquet` carries **no organism column**. Staged without one, the corpus row
is `NULL` and the search is invisible on FRAN's species page. Spectronaut reports self-resolve
from `PEP.AllOccurringOrganisms`; DIA-NN cannot.

The user already confirmed the organism at step 3, and `fetch_fasta.py` wrote it to
`<fasta>.meta.json` — read automatically, overridable with `--organism`/`--taxon`. Never
invent one: an unknown organism is **absent** from the manifest, not guessed (architectural
rule #2).

**Which sidecar is read.** It has to be tied to the FASTA the search actually used: `fasta` in
`search_provenance.json`, or `--fasta` on the command line DIA-NN echoes at the top of
`report.log.txt`. The old rule took the first `*.fasta.meta.json` near the search in sort order,
and `PROT_0793/` holds human, mouse and mouse+contaminant sidecars side by side. The mouse
search's manifest recorded the **human** database (seen in the drop dir, 2026-09-24). Now:

- the search's own FASTA + `.meta.json`, then a nearby sidecar whose file name or recorded
  `fasta` matches it. Nearby means the search dir, its parent, `<parent>/input/`, and for a
  session's `output/search` the session's `input/`.
- if the search's FASTA cannot be found, a **lone** nearby sidecar is used. Two or more is a
  guess, so organism and database stay absent.
- `fasta_path` is **the file the search read**, whenever the search names it. The md5 and
  entry count come from the sidecar. A sidecar records wherever `fetch_fasta.py` first wrote the
  FASTA, and in the drop dir on 2026-09-24 that differed three ways, all with identical md5s:
  a laptop `/Users/...` path (hive_remote), and the user's staging copy
  `~/proteomics-pipeline/staging/search.fasta`, which the next session overwrites (twice), where
  the search read `<session>/input/search.fasta`.

A sidecar that itself records `organism: ""` (PROT_0793's human one does) gives no organism.
That is honest, not a bug.

## Why `check` refuses (stable codes, never an exception)

| reason | meaning |
|---|---|
| `not_core_facility` | the drop directory is not writable → a collaborator's run. **Correct behaviour, not an error.** Do not work around it. |
| `not_on_hive` | the search directory does not exist here — you are not on HIVE, or the path is the local one |
| `search_incomplete` | no report, a zero-byte one, or FragPipe/Radiant without its completion marker (a cancelled FragPipe run keeps its report). A failed or partial search must never be ingested |
| `engine_unsupported` | Sage (DDA) and AlphaDIA have no FRAN corpus adapter. The corpus is DIA |
| `no_drop_dir` | the drop directory does not exist and could not be created |
| `already_staged` | a receipt exists; `--force` to re-stage |
| `opted_out` | `FRAN_DEPOSIT=off` or `--skip` — recorded in the receipt so `backfill` honours it |
| `qc_run` | a QC run (see *QC runs are never handed over*); `--not-qc` if the rule is wrong |
| `not_a_skill_search` | `backfill` only: a search this skill did not run |
| `already_ingested` | `backfill` only: the cron's logs show FRAN already has it |
| `needs_agent_check` | `backfill` only: FragPipe/Radiant with no completion marker to read; never staged unattended |
| `drop_dir_not_writable` | a `proteomics-grp` member cannot write the drop dir: a permission problem, reported and never recorded |
| `entry_not_writable` | the entry was staged by another account without group write: reported, never a crash |

`not_core_facility` is recorded in the receipt only when the account is **known** to be
outside `proteomics-grp`, or is a teaching account. A member who cannot write the drop dir gets
`drop_dir_not_writable`. When membership cannot be read, the refusal stands for that call only
and is not recorded.

**Permissions.** Everything the skill creates under `incoming/` is set explicitly to `2775`
(dirs) and `0664` (manifests), and receipts to `0664`, whatever the umask. With a `022` umask
the next Core member could neither withdraw nor re-stage an entry.

**An explicit `--fasta-meta` must describe the search's database too.** If it does not match
the search's own FASTA, it is ignored with a warning (`fasta_meta_ignored`) and organism and
database stay blank. They are never taken from another search.

An ineligible run is **not** a failure of the analysis. Note it in one line and carry on with
DE — never block, retry, or ask the user to fix it.

`stage` records the three *decisions* (`opted_out`, `not_core_facility`, `qc_run`) in
`<out>/fran_deposit.json` so that a `backfill` months later, run by someone else, does not
hand the search over anyway. Neither blocks an explicit `stage`, and neither overwrites a
receipt that says the search is staged or ingested.

## Idempotency

Every manifest carries **`staged_at`**, the ISO-8601 UTC time the search was first handed over,
plus `staged_by`. FRAN's runner orders the drop box oldest-first by it (its fallback is the entry
directory's mtime, which every re-stage resets). A re-stage keeps the original `staged_at` and
`staged_by` and adds `restaged_at` / `restaged_by`. An entry staged before the field existed
gets the time its manifest was written. All other manifest keys are unchanged. Every manifest `stage` writes passes a copy of FRAN's own
`read_manifest` validation (tests/test_fran_health_backfill.py `FranManifestContractTests`):
`fran_manifest_version` 1, `qc`/`exclude` JSON booleans, `staged_at` ISO 8601, and
`search_name`/`organism` either a non-empty string or absent. An empty `--name ""` is written as
absent, because FRAN rejects `""` as malformed.

The entry name is deterministic — `<search dir name>__<8 hex of its real path>` — so
re-staging reuses the same path rather than presenting the cron with a second candidate that
would ingest as a duplicate search. Re-staging relinks from scratch, so an entry never ends up
holding a mixture of two runs. A receipt at `<out>/fran_deposit.json` stops a resumed session
staging twice.

The manifest's `output_dir` is the **real** search directory, so `corpus_ingest.py --output-dir`
keys idempotency and provenance on where the search actually lives, not on the handover path.
FRAN's `auto_ingest.py` does not read the manifest yet (checked on `main`, 2026-09-24). It
passes `realpath(<the dir it scanned>)` as `--output-dir`, and a drop entry is a real
directory, so a staged search lands in the corpus under `incoming/<entry>`. `verify`'s
database lookup therefore asks under both names.

## Follow-ups
- **XIC storage.** Every DIA-NN search extracts chromatograms — `--xic` is forced by
  `run_search.ensure_xic()` since v2.3.0 (see above), so it is not something a cfg can leave
  out. The cost is disk: 27 GB for the 399-run cohort above. How long they are kept is the
  facility's retention decision, not the skill's.
- **Local (non-HIVE) searches** have no drop directory to write to and are not handed over.
