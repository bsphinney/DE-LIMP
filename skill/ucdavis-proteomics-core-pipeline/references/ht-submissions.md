# High-throughput submissions — searching a whole Core plate run

**When this applies:** a UC Davis Proteomics Core member gives you a **submission number**
("search 0793", "run the HT plate for 0793") instead of a folder of raw files. Everything
here replaces step 1 only. From step 2 onward the run is completely ordinary — same
acquisition detection, same defaults, same confirmation, same FRAN handover.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 21: Why a submission is not a folder
- 69: 1. Get the file list
- 93: 2. Honour the exit code — this is the whole point of the gates
- 194: 3. Organism and FASTA — still asked, never inferred
- 208: 4. Search — the DIA-NN parallel chain, automatically
- 235: 5. Deposit to FRAN, then tell STAN
- 273: Deployment status
- 286: Two ways in — and which one a Core member can actually use
- 361: 1a. The HT steps as the pipeline runs them (step 1a)
<!-- /contents -->

## Why a submission is not a folder

Globbing a directory gets an HT submission wrong in four ways, none of which announce
themselves:

* **A submission can span two trays**, and the second tray's filenames may never mention
  the submission number. Its extent is inferred from the acquisition counter, not stored.
* **Not every run is a customer sample.** Well blanks and HeLa standards sit on the same
  plate. Searching them as samples pollutes the protein matrix and the corpus.
* **Some samples are already known to be bad.** STAN flags them as needing another
  injection.
* **Re-injections live elsewhere in the name.** A re-run is
  `20260901_0793_rerun_SI-48_S1-A1_1_24200.d`; STAN maps it back to the same submission.

STAN knows which files a submission is. This skill knows how to search them. Neither
learns the other's job — so ask STAN for the list, then proceed exactly as normal.

### Measured, on submission 0793 (2026-08-31)

Do not talk yourself into `ls | grep 793`. Against STAN's answer of **120 files**, a naive
substring grep over the acquisition directory returns **97** — and is wrong twice over:

| | files | |
|---|---|---|
| STAN (truth) | **120** | 88 on tray S6, 32 on tray S5 |
| naive `grep 793` | 97 | |
| **missed** | **32** | the S5 continuation — `20260828_100spd_COH-12_S5-D2_1_24157.d` and friends. Their names never contain 793. |
| **wrongly included** | **9** | `20260827_793_100spd_Hel50_S6-A12_...` — **HeLa standards**, not customer samples |

So the shortcut drops a quarter of the submission *and* searches nine standards as if they
were samples. Neither shows up as an error; the search just succeeds on the wrong files.

**Three traps behind that:**

* **A submission spans trays that never name it.** 0793 ran S6 on 2026-08-27
  (`..._1_24026.d` … `_1_24125.d`) and continued onto S5 the next day. The link is the
  **acquisition counter** — S5 begins at 24126, one after S6 ends — not the filename.
* **Substring matching hits other customers.** `793` also matches `..._1_23793.d` and
  `..._1_22793.d`. STAN matches the submission as a **delimited token with leading zeros
  ignored**, which is why `0793` and `793` both work and neither over-matches.
* **Tray labels repeat.** `S6` is a tray *position*, reused every month since December. It
  is not a submission identifier.

Since 2.12.0 a tray named without the number (0793's S5 above) is **not** searched from STAN's
list: `named_for_submission` fails and lists its runs, and the operator confirms them through a
curated list (`core_submission.py locate --files-from`) — because a tray STAN joins by position
can just as well be another project's (2026-10-05, below).

## 1. Get the file list

`ht_manifest.py` runs **on HIVE** (like `fran_deposit.py`), because that is where STAN,
its credential, and the raw files are. Drive it through `hive_exec.sh`:

```
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py fetch 0793 \
    --summary ~/core/PROT_0793/submission_summary.json --out ~/ht0793'
```

`--summary` is the submission's CoreOmics record, written by `core_submission.py fetch` and put
on HIVE (`references/core-submissions.md`). It is what lets `ht_manifest.py` see a sample STAN
has no run for (gate `sample_list`, below); without it that check is not made, and says so.

`--include` takes `samples` (**default** — blanks and standards excluded), `rerun` (only
what STAN flagged for re-injection), `standards`, or `all`. Both `0793` and `793` work.

It writes two files into `--out`:

| file | what it is |
|---|---|
| `files.txt` | one absolute, resolved `/nfs/...` path per line — **only the runs named for the submission** (`named_for_submission`, below); fed to `core_submission.py locate --files-from`, or to `--raw` / `--files` |
| `ht_manifest.json` | the full STAN payload (plates, wells, per-run class, verdicts) **plus** the gate results |

## 2. Honour the exit code — this is the whole point of the gates

| code | meaning | what to do |
|---|---|---|
| **0** | safe | proceed to step 2 of the normal flow |
| **2** | **HARD GATE FAILED** | **do not search.** Show the operator the failing gate and stop |
| **3** | STAN unreachable, or too old to have `ht-manifest` | fix STAN; nothing was searched |
| **4** | (`record` only) write-back failed | the **search is fine**; only STAN's HT tab is stale |

The hard gates exist because each one otherwise produces a search that **succeeds while
covering the wrong files**:

* **`missing_paths`** — runs STAN knows about but has no raw path for. They are *excluded*
  from `files`, so a non-empty list means searching a subset and reporting success.
* **`n_files`** — nothing matched. Usually a mistyped submission number.
* **`paths_exist`** — a path STAN resolved that the filesystem does not have. Checked here
  because otherwise a 120-file SLURM array dies partway through, hours in.
* **`repeated_paths` FAIL** — STAN lists one file more than once with *different* details
  (well, run name, class, verdict): which sample it is cannot be read from the manifest.
* **`named_for_submission`** — a run whose name does not carry the submission number right after
  its date (`<date>_<n>_…`, `<date>_PROT_<n>_…`; the same rule `core_submission.py locate` uses to
  recognise a plate). Every such file is listed, with STAN's class, plate, well and injection, and
  none of them is written to `files.txt`.
* **`one_method`** — the file names declare more than one acquisition method: another throughput
  (`100spd` beside `60SPD`) or another mode (`DIA` beside `DDA`), each group listed with its count
  and how many of it are named for the submission. In a run named for the submission the method
  is read **only right after the plate token** (`<date>_[PROT_]<n>_<speed>_…`), never in the
  sample name after it: a re-injection `<date>_<n>_rerun_Liver_24m_3_…` declares no method (WARN
  below), not `24m` (before 2.12 it stopped the plate as two methods). DIA vs DDA is read again
  from the data at step 2 (`detect_acquisition.py`).
* **`sample_list`** (with `--summary`) — a sample on the CoreOmics record with no run named for
  the submission, or such a run that names no sample on the record. Each is named. Runs are
  matched to samples by the rule `core_submission.py locate --files-from` applies (`sample_name`
  inside the run's sample field), so the two never disagree. A run that names **several** samples
  (two samples share a `sample_name`, or a label looks like a plate well) is not a run of no
  sample: it is `sample_list` **INFO** (`runs_of_several_samples`), and staff say whose it is when
  the list is staged (`locate --files-from … --choose <unique_id>=<file>`). Not checked for
  `--include rerun|standards`. A `--summary` that cannot be read is a FAIL that says to put
  `submission_summary.json` on HIVE and pass its path — nothing to ask the operator.

### Why STAN's list is checked, not trusted (2026-10-05)

A 96-sample plate on tray S3 came back from STAN as **136 files with every gate PASS**: 95 of
the plate's 96 runs, one real run **missing**, and **41 runs of other projects** that had used tray
S3 on earlier days — named `<MMDDYYYY>__60SPD_DIA-<code>…` (5 of them `__30SPD_DDA…`), other
dates, other methods. Searched as given, other clients' data would have gone into the search, the
report and FRAN, and a sample would have been left out silently. It was caught by hand.

How STAN produced it (STAN `expand_submission_runs`, read 2026-10-05): it takes the runs whose
names carry the number, then **every run on the same tray** within ±300 injections, from the
**tray's earliest run in that window** — not from the submission's first run. Tray labels are
reused constantly, so a tray another project used a few days earlier is swept in whole. (Its rule
that a tray naming another submission never joins applies only to the further trays it adds, not
to the submission's own tray, and recognises only names shaped `<8 digits>_<number>_`.) The
missing run was not in STAN's records at all, so no STAN-based gate could see it; only the
CoreOmics sample list can.

So: **when `named_for_submission` fails, treat STAN's tray assignment as suspect** and ask the
operator whose those runs are. Runs that ran *before* the first run named for the submission
(by the acquisition counter) are called out: they cannot be a later tray of it.

**Getting past a failed gate is staff choosing the files.** There is no flag that accepts STAN's
list as it is. Write the runs the operator confirms — `files.txt` plus, say, a continuation tray
named without the number, or a sample's run that STAN never recorded — to a list, and stage it with
`core_submission.py locate --files-from <list>` (`references/core-submissions.md`), which matches
each run to its sample and records the staff choices in `locate.json`. When no listed run is
named for the submission at all, `n_files` fails and points there too.

And the warnings, which need a human rather than a stop:

**Repeats are okay, but they are flagged** (Brett, 2026-10-01) — never dropped silently, never a
stop on their own:

* **`repeated_paths` WARN** — STAN listed the same file more than once (a plate of 96 came back
  as 100 lines, one run 4 times, all gates PASS). No engine can take a file twice, so
  `files.txt` holds it **once**; `ht_manifest.json` records `n_files_listed`, `n_files` and
  `repeated_paths` (each path with its count), and stderr flags each. Say so — STAN's `counts`
  still include the repeats.
* **`repeated_names` WARN** — different files share a run name (a rerun in another folder, one
  name on two plates). **All are kept** in `files.txt` and listed under `repeated_names` (each
  name with its paths). But a search cannot take them as given — DIA-NN names a run by its file
  name without the folder, so they would become one run — so `run_search.py` stops before
  searching, naming them. Ask the operator: give one another name (a link under another name
  works), or search them separately.
* Every search route follows the same rule (`check_report_runs`: `distinct_inputs`,
  `repeated_names`, one wording): `run_search.py` for every engine, `diann_parallel.py`,
  `radiant_parallel.py`. A repeated path is searched once and recorded in
  `search_provenance.json` `repeated_paths`, which the analysis report shows.
* **`plates`** — more than two trays. Confirm the extent; it is inferred, not recorded.
* **`one_method` WARN** — a file name declares no method (a re-injection is named
  `<date>_<n>_rerun_…`), so it could not be compared with the rest. Confirm with the operator.
* **`sample_list` WARN** — not checked, because no `--summary` was given.
* **`counts`** — fewer than 12 customer samples. A mistyped submission looks exactly like
  this. A genuinely small submission is legal, so this warns rather than blocks.
* **`needs_rerun`** — samples STAN flagged for another injection. They **are** included in
  the default `samples` set. Say so; the operator may want them excluded, or may want
  `--include rerun` to search only those.

Surface every non-PASS gate to the user before committing compute. This is the same rule
as golden rule #1, just with STAN's evidence attached.

## 3. Organism and FASTA — still asked, never inferred

**STAN does not know the organism and does not guess one.** The skill's rule against
fabricating parameters applies unchanged: ask the operator, resolve with
`fetch_fasta.py resolve`, confirm the proteome.

Pre-staged FASTAs live in `/quobyte/proteomics-grp/de-limp/fasta/` (20 files as of
2026-08-31 — human, mouse, chicken, pig, cow, dog and others, most as `_opg_` one-per-gene
sets with a `.provenance.json` beside them), plus `/quobyte/proteomics-grp/MRS/` for human
± contaminants. `fetch_fasta.py --hive` reuses only a proteome in `MRS/` (a file named
`<UPID>*.fasta` there); it does not look in `de-limp/fasta/`. For any other organism let
`fetch` download it (seconds) rather than passing a `de-limp/fasta/` file with `--path`,
which records no organism or content type in the sidecar.

## 4. Search — the DIA-NN parallel chain, automatically

Nothing HT-specific. Pass `files.txt` to the normal flow and `run_search.py` routes itself:
a plate is far more than five files on a machine with SLURM, so it takes the **5-step
parallel chain** (per-file passes as a job array) rather than one long single-node job.

The chain keeps the two settings DIA-NN would otherwise optimise per file the same across
the plate (steps 3/5 reuse the per-file `.quant` files, so per-file values would be stitched
together inconsistently). Nobody measures anything by hand first:

* **Mass accuracy** is pinned in the cfg by `estimate_params.py` (step 6b). When it planned
  `measure_with_diann` instead (an Orbitrap level outside DIA-NN's table), the chain's step 1b
  measures it on the same representative runs as the window and pins max(measured, SOP).
* **The scan-window radius** is measured by the chain's step 1b: `probe_window.py` is handed the
  whole plate, measures three representative runs (the median and quartile runs; never a blank,
  wash, failed injection or damaged `.d`) and pins the median radius into steps 2–5. It never
  probes just one file: left to DIA-NN per file, an 18-file run got radius 7 for seventeen files
  and 8 for one, and with 120 files a split is more likely, not less. An omitted `--window` does
  not decline the chain.

Both are described in SKILL.md step 7 and `references/diann_parallel.md` ("Scan window MUST be
pinned", "Which runs step 1b measures").

Watch it with `watch_run.sh --all` (step 7b). A plate-sized array will usually lose a file
or two to a stalled node; the documented playbook is to retry once, then drop that file and
record it in **Data Quality Notes** rather than restarting the cohort.

## 5. Deposit to FRAN, then tell STAN

FRAN handover is unchanged — `fran_deposit.py stage` symlinks the output into
`/quobyte/proteomics-grp/fran/incoming/` and FRAN's cron ingests it (`references/fran.md`).

> ⚠ **Never write HT search output under `/quobyte/proteomics-grp/STAN/`.** That subtree is
> on FRAN's `DEFAULT_EXCLUDES` list deliberately: STAN writes a DIA-NN `report.parquet` for
> every QC run, and ingesting those as customer searches would corrupt every corpus count.
> Write somewhere else and let `fran_deposit.py stage` link it in.

**There is no write-back, and that is deliberate.** STAN v1.0.42 briefly added an
`ht_searches` table and an `ht-record-search` command; **v1.0.43 reverted both**. The
reasoning is worth keeping because it decides how this step works:

> A table in STAN would be a second copy of what FRAN already knows, and the copy is the
> one that goes stale — it depends on the skill remembering to call it, and says nothing
> when that step is skipped. It would also mean re-implementing FRAN's authorization, or
> proxying its confidential data, when FRAN already decides who may see a search.

So **STAN links and FRAN answers.** The submission number is in the raw filenames, so
FRAN's `/api/internal/submission/{id}` already returns the submission plus every search
under it — search id, engine, organism, precursor and protein counts. Nothing to record.

Do confirm the loop closed, so a plate that silently failed to ingest is visible:

```
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py link 0793'
```

**Exit 4 is not a search failure.** It means "no search visible under this submission
yet", which is the *expected* answer while `fran_deposit.py` reports
`staged_pending_cron` — FRAN ingests on its next scan. Re-run it later.

⚠ **A 404 from that endpoint does not mean "not ingested".** It is internal-only: FRAN
returns 404 to a caller who is not signed in, and to a lab user asking about another lab's
submission, precisely so it never confirms another lab's data exists. Sign in (Entra) or
pass `--cookie` before treating a 404 as evidence of anything.

## Deployment status

**STAN v1.0.43 is deployed on HIVE** at `/quobyte/proteomics-grp/brett/stan_venv/bin/stan`,
with `ht-manifest` available — verified 2026-08-31 end to end against submission 0793:
120 files, plates S5/S6, `counts {sample:120, standard:8, blank:7}`, 6 flagged for
re-injection, every hard gate PASS. `--include rerun` returns 6; `--include all` returns
135 (120+8+7); a mistyped `9999` hard-fails with exit 2.

There is no database migration to apply — the `ht_searches` table was reverted in v1.0.43
(see above). If you meet an older STAN, `ht-manifest` returns *"No such command"* and
`ht_manifest.py` exits **3** naming the version and the fix, rather than looking like an
empty submission.

## Two ways in — and which one a Core member can actually use

**The auth is Microsoft Entra, not CAS.** STAN and FRAN are both Azure App Service Easy
Auth (`/.auth/login/aad`), gated on an Entra **group object id**. There is no CAS anywhere
in either app, so "log in with your UC Davis account" means Entra.

| | how it authenticates | works headless on HIVE? |
|---|---|---|
| **CLI** (default) — `stan ht-manifest` | Postgres credential | **only for the token's owner** |
| **HTTP** — `--http https://ucd.stan-proteomics.org` + `--share-token-file` | per-submission HMAC share link | **yes** |
| **HTTP** + `--cookie` | Entra session cookie from a signed-in browser | awkward — cookie expires |
| **HTTP** with neither | Entra sign-in (browser redirect) | no |

### The CLI credential is the service account's secret: never copy or share it

`/quobyte/proteomics-grp/brett/.pgfarm_token` is **not** a short-lived token, whatever its
name says. It holds the **long-lived 512-character secret of the
`genome-proteomics-service-account`** (STAN's `CLAUDE.md`, "PG Farm auth"), which STAN exchanges
for a fresh JWT on every use. It is mode `0600` on purpose. Anyone who can read it can act as
the service account, which holds `SELECT, INSERT, UPDATE, DELETE` on STAN's tables, with no
per-person audit trail, for as long as the secret is valid.

**Never copy it, and never make it group-readable**: not into `/quobyte/proteomics-grp/etc/`,
not into a home directory, not into a session folder. Skill 2.9 and earlier called it a
"7-day token" and suggested publishing a group-readable copy at
`/quobyte/proteomics-grp/etc/pgfarm_token`. That was wrong, and `ht_manifest.py` no longer
looks there.

`ht_manifest.py` looks for a credential, in order:

1. `$PGPASSWORD`
2. `$STAN_PG_TOKEN` (a file)
3. `~/.pgfarm_token`
4. `/quobyte/proteomics-grp/brett/.pgfarm_token` (readable by its owner only)

Every other Core member gets *permission denied* at (4). That is expected, and it does not mean
"that submission does not exist".

**So everyone who is not the owner uses the HTTP path with a share token:**

```
printf '%s\n' '<tok>' | bash scripts/hive_exec.sh 'umask 077; mkdir -p ~/.stan && cat > ~/.stan/share_0793'
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py fetch 0793 \
    --summary ~/core/PROT_0793/submission_summary.json --http https://ucd.stan-proteomics.org \
    --share-token-file ~/.stan/share_0793 --out ~/ht0793'
```

**The token goes in a file, never on a command line.** A token in argv is in the process list,
and in `commands.log`, which the run registry copies into a folder the whole Core group can
read. The first line above passes it through stdin (`printf` is a shell builtin, so it never
shows in a process list either) into a mode-600 file. Do not write that line into
`commands.log`. `--share-token <tok>` still works, with a warning, and so does `STAN_HT_SHARE_TOKEN`.
`--cookie-file` does the same for an Entra session cookie. `ht_manifest.py` masks the token,
the cookie and the PG credential in every line it prints and in `ht_manifest.json`, whatever
the server echoes back, in an error or in a normal reply's plates, counts and example paths.

Get `<tok>` from the submission's HT tab in the dashboard — the token is an HMAC of the
submission number, so a link for 0793 opens 0793 and nothing else, and rotating
`STAN_HT_SHARE_SECRET` invalidates every outstanding link at once. There is currently **no
CLI to mint one**, so it comes from the web UI.

Without either credential the endpoint answers with a clean 403 naming the login URL —
verified live against `ucd.stan-proteomics.org`:

> `{"detail": "High-throughput submission data requires an authorized sign-in.",
> "login_url": "/.auth/login/aad?post_login_redirect_uri=/"}`

The durable fix for the CLI path is a group-readable token under
`/quobyte/proteomics-grp`; until then the HTTP + share-token route is what makes this
workflow usable by the Core rather than by one person.

Less-used `ht_manifest.py` flags: `--stan <path>` (the stan binary), `--token <file>` (the STAN
Postgres credential; or `STAN_PG_TOKEN` / `PGPASSWORD`), `--fran <url>` (the FRAN base URL for
`link`; or `FRAN_URL`).

## 1a. The HT steps as the pipeline runs them (step 1a)

If the user gives a **submission number** instead of a folder ("search 0793", "run the HT
plate"), the file list comes from STAN, not from a directory listing. **Do not glob.** A
submission can span two trays whose second tray's filenames never mention the number; it
contains well blanks and HeLa standards that must not be searched as samples; and STAN
already knows which samples want re-injecting. Globbing gets all of that wrong, silently.

```
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py fetch 0793 \
    --summary ~/core/PROT_0793/submission_summary.json --out ~/ht0793'
```

`--summary` is the CoreOmics record (`core_submission.py fetch`, put on HIVE as in step 1c):
with it, every sample on the record must have a run and every run a sample.
**If you are not the owner of STAN's Postgres token** (mode 0600 — i.e. almost everyone),
that fails with *permission denied*, which reads exactly like "no such submission". Use the
hosted dashboard instead, with a per-submission share token from the HT tab — the one auth
path that works headless. Sign-in is **Microsoft Entra, not CAS**:

```
printf '%s\n' '<tok>' | bash scripts/hive_exec.sh 'umask 077; mkdir -p ~/.stan && cat > ~/.stan/share_0793'
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py fetch 0793 \
    --summary ~/core/PROT_0793/submission_summary.json --http https://ucd.stan-proteomics.org \
    --share-token-file ~/.stan/share_0793 --out ~/ht0793'
```

The token goes into a mode-600 file on HIVE through stdin, and the search reads it from
there. **Never put it on a command line or in `commands.log`**, and don't log the `printf` line
above. The run registry copies `commands.log` into a folder the whole Core can read, and a
token on a command line is in the process list too. `ht_manifest.py` masks the token in
everything it prints.
Writes `files.txt` (absolute resolved paths → feed to `--raw`/`--files`; **only the runs named
for the submission**) and `ht_manifest.json` (STAN's payload + gate results). `--include` = `samples` (**default**,
excludes blanks/standards), `rerun`, `standards`, `all`. Both `0793` and `793` work.

**Honour the exit code.** `0` proceed · **`2` HARD GATE FAILED — do not search** · `3` STAN
unreachable or too old. A hard gate means the search would *succeed while covering the wrong
files*: `missing_paths` (runs STAN knows but has no path for — silently excluded from the
list), `n_files` (nothing matched — usually a mistyped number), `paths_exist` (a resolved
path the filesystem lacks — otherwise a 120-file array dies hours in). **STAN's list is
checked, not trusted:** STAN assigns runs to a submission by *tray*, whoever's they are, and in
2.10.0 a 96-sample plate came back as 136 files — 41 runs of other projects that used the tray
earlier, at other methods (5 DDA), and one real run missing — with every gate PASS. So
`named_for_submission` FAILS on any run whose name does not carry the number (`<date>_[PROT_]<n>_`)
and lists every one; they are never in `files.txt`. `one_method` FAILS when the names declare
different methods (`100spd` beside `60SPD DIA`, DIA beside DDA). `sample_list` (with `--summary`)
FAILS naming each sample with no run and each run of no sample. When no listed run is named for
the submission, `n_files` FAILS and says so. **The way past any of these is staff choosing the
files**: write the confirmed list and stage it with `core_submission.py locate --files-from`
(their choices are recorded) — never search STAN's list as given. **Surface every
non-PASS gate before committing compute**, warnings included: `plates` >2, `counts` <12
samples, `sample_list` not checked (no `--summary`), `one_method` for a name that declares no
method (a `_rerun_` run), `needs_rerun` (those samples ARE in the default set — say so), `repeated_paths` (STAN
listed a file more than once; `files.txt` holds it once and it is searched once) and
`repeated_names` (different files share a run name: all kept and flagged, but `run_search.py`
stops before searching them as given — ask the operator to rename one or split the search).
**A CoreOmics HT submission is then staged from that list**: `core_submission.py locate
--files-from ~/ht0793/files.txt` (step 1c). HT run names carry the submitter's `sample_name`,
not the `unique_id`, so `--files-from` matches this submission's plate runs on `sample_name`
(gate `ht_plate_runs`; pass `--ht-manifest ~/ht0793/ht_manifest.json` so `--reinjections
latest` knows STAN's `needs_rerun`); samples that share a name get `choose_run` — ask staff, then
`--choose <unique_id>=<file>`. Never `--allow-partial` for samples that were run.

Then continue with step 2 exactly as normal — nothing after this point is HT-specific.
**Organism is still asked, never inferred**: STAN does not know it. Reuse a pre-staged
proteome from `/quobyte/proteomics-grp/de-limp/fasta/` (step 6, `--hive`). The plate will
route itself to the **5-step DIA-NN parallel chain** at step 7, whose step 1b measures the scan
window (and a planned mass accuracy) on representative runs itself; nothing is measured by hand
first (section 4 above).
