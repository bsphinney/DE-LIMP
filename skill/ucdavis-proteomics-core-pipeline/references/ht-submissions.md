# High-throughput submissions — searching a whole Core plate run

**When this applies:** a UC Davis Proteomics Core member gives you a **submission number**
("search 0793", "run the HT plate for 0793") instead of a folder of raw files. Everything
here replaces step 1 only. From step 2 onward the run is completely ordinary — same
acquisition detection, same defaults, same confirmation, same FRAN handover.

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

## 1. Get the file list

`ht_manifest.py` runs **on HIVE** (like `fran_deposit.py`), because that is where STAN,
its credential, and the raw files are. Drive it through `hive_exec.sh`:

```
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py fetch 0793 --out ~/ht0793'
```

`--include` takes `samples` (**default** — blanks and standards excluded), `rerun` (only
what STAN flagged for re-injection), `standards`, or `all`. Both `0793` and `793` work.

It writes two files into `--out`:

| file | what it is |
|---|---|
| `files.txt` | one absolute, resolved `/nfs/...` path per line — feed straight to `--raw` / `--files` |
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

And the warnings, which need a human rather than a stop:

* **`plates`** — more than two trays. Confirm the extent; it is inferred, not recorded.
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

Two things the chain requires, both of which matter more at plate scale — read
`references/diann_parallel.md` before submitting:

* **Mass accuracy must be pinned**, not auto. Steps 3/5 reuse the per-file `.quant` files,
  so anything DIA-NN auto-optimises per file gets stitched together inconsistently.
* **Scan-window radius must be measured**, not left to auto — `probe_window.py`, once, on
  one file. On an 18-file run DIA-NN inferred radius 7 for seventeen files and 8 for one,
  and the chain combined them. With 120 files the chance of a split is higher, not lower.

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
    --http https://ucd.stan-proteomics.org --share-token-file ~/.stan/share_0793 --out ~/ht0793'
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
