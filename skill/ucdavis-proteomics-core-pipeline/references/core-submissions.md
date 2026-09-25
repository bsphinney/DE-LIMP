# CoreOmics submissions — from "search submission 807" to a shared Bioshare folder

**When this applies:** a UC Davis Proteomics Core staff member names a **CoreOmics
submission** — "search the data from submission 807", "analyze PROT_0807", "deliver the
results to the collaborator", "put it in Bioshare". SKILL.md steps 1c and 12b are the short
form; this is the runbook and the evidence behind it. Everything marked *measured* or
*verified* was checked live against CoreOmics, HIVE, the Flinders share or upstream source on
2026-09-16.

A **high-throughput plate** is a different route: STAN knows its file list
(`references/ht-submissions.md`). `core_submission.py locate` notices a plate and exits 4.

## Why this needs its own script

A service run attaches a PI, a sample sheet and a Bioshare folder to a set of raw files.
When any link in that chain is wrong, nothing errors:

- **The raw files have to be found.** A submission records each sample's `unique_id`, not a
  path. Ids are short, reused across submissions, and collide with plate wells — a glob
  searches another lab's runs and the search succeeds.
- **The service directory is organized by people, not a formula** (`McDonald karen`,
  `UCSF/Feeley_lab`), so "where does PROT_0807 live" is a lookup that can be wrong.
- **A share built from links depends on server settings nobody watches.** Links into
  `/quobyte` served nothing over https, and SMB hides absolute links (the PROT_0793 lesson,
  below).

`scripts/core_submission.py` does the bookkeeping deterministically and hands every judgment
call back as an **exit code plus the exact question**. It never runs a search engine; steps
2–9 of the normal flow sit in the middle.

## Where each step runs

| step | runs | why there |
|---|---|---|
| `identify` | **local** | reads names offline; a sample-id lookup needs the CoreOmics token |
| `fetch` | **local** | the CoreOmics token is on the staff member's computer; `~/.coreomics_token` does **not** exist on HIVE (checked) |
| `locate` | **HIVE** | reads the Flinders `raw_data` tree |
| `stage` | **HIVE** | writes links into the service directory |
| `conditions` | either | wraps `collect_conditions.py`; run it locally on the pulled `sample_files.tsv` |
| search, DE, figures, audit, `make_methods.py` | **HIVE (SLURM)** | in the session under the printed `work_dir` — the session of record |
| writing the report, `make_analysis_html.py`, `to_docx.py` | **local** | on a pulled copy; the finished files are pushed back |
| `deliver` | **HIVE** | copies into the Flinders share dir (big copies → `sbatch`) |
| `bioshare`, `email-draft` | **local** | CoreOmics API again; drafts never send |

Every subcommand prints one JSON object to stdout and notes to stderr. `stage`, `deliver` and
`bioshare ensure|send` are dry runs until `--apply`.

## What staff need first

- **A CoreOmics API token with staff access**, saved on **their own computer** as
  `~/.coreomics_token` (`chmod 600`), or exported as `COREOMICS_TOKEN`. Not on HIVE. How to
  obtain one is not documented here yet: ask the Core's CoreOmics administrator. Without it,
  `identify` still reads ids in names and the agent asks the rest (section 0).
- **A HIVE account in `proteomics-grp`.** The service directory is
  `gc-prot-core-user:proteomics-grp`, mode `2775`.
- **The whole scripts directory on HIVE at `~/proteomics-pipeline/scripts/`**
  (`references/access.md`) — not `~/pipeline/`, which does not exist. `core_submission.py`
  imports `session.py`, `collect_conditions.py` and `run_search.py`; a partial copy exits 3
  and says to sync the directory.
- **For HT plates only:** a STAN per-submission share token (`references/ht-submissions.md`).
- `hive_exec.sh` reuses one SSH connection for 10 minutes. HIVE's sshd throttles rapid new
  connections (MaxStartups): eight quick separate `ssh` calls were refused with
  `kex_exchange_identification: read: Operation timed out` for 20+ minutes on 2026-09-16,
  while the HIVE status page said operational. Native Windows OpenSSH has no ControlMaster,
  so reuse is off there; `HIVE_SSH_MUX=0` turns it off anywhere. `--get` copies with `-rlt`,
  not `-a`: from HIVE's setgid (2775) folders, `-a` made macOS exit 23 with `fchmodat ...
  Operation not permitted` even though every file arrived (measured).

Environment overrides (the tests use them; staff never need to): `COREOMICS_BASE_URL`,
`COREOMICS_TOKEN`, `CORE_FLINDERS_ROOT` (default `/nfs/lssc0/flinders/proteomics` — where
*this machine* does file work), `CORE_WORK_ROOT` (default `/quobyte/proteomics-grp/SERVICE`).

**Server-side paths never come from `CORE_FLINDERS_ROOT`.** Bioshare knows a share only by its
HIVE path, so `share_dir` in the summary and Bioshare's `link_to_path` are always built with
forward slashes from `/nfs/lssc0/flinders/proteomics`. A staff member on Windows, or with the
root pointed at an SMB mount, still registers the right path.

## 0. `identify` — which submission is this data? (local)

Every report carries its submission, so the number is settled at SKILL.md step 1 — for any
Core data, not only a staff "search submission 807" run. Never guessed:

```
python3 scripts/core_submission.py identify <raw files and/or their folder> --text "<the user's message>"
```

1. **Named ids.** A `PROT_####` token (`PROT_0756`, `prot-756`) or a 12-character CoreOmics id
   (with a letter and a digit — a 12-digit timestamp is not one) anywhere in the paths, or
   "submission 756" in the message (3-4 digits; `prot_10ug` is an amount, not PROT_0010). A bare number in a FILE name never counts: Exploris
   runs carry counters (`Ex08312026_380_JE21`), and a 12-hex tail of a UUID/GUID is not an id.
   With a token every named id is looked up: several that are one submission are merged; one
   CoreOmics does not know is `not_found`; one whose sample IDs are in none of the file names
   is `named_unconfirmed`. Without a token a single named id is `named` (unverified).
2. **Sample ids (token).** Otherwise it lists the submissions made up to 240 days before the
   runs (`--max-days`) and matches their sample ids in the file names with `locate`'s rules run
   in reverse: delimited tokens, the timsTOF sample field only, the longest id owns a run, weak
   ids (`A3`, `001`) are no evidence, and a run counts only for a submission made on or before
   the date in its name. Two submissions that could own one run → `ambiguous` (`locate`'s
   `ambiguous_label`). A name with no date matches every submission using the id. `matched`
   needs the one candidate to own at least half the files and two sheet IDs, and the whole
   window to have been listed; otherwise `weak`. A folder is listed one level deep for names.

| exit | status | what to do |
|---|---|---|
| 0 | `named`, `matched` | confirm in one line — the JSON's `ask` (it says how many file names carry the sample IDs) |
| 2 | `none`, `ambiguous`, `weak`, `named_unconfirmed`, `not_found` | ask the user for the number — the JSON's `ask` |
| 3 | `needs_token`, `lookup_failed` (CoreOmics refused or unreachable) | ask for the number, relay `token_help`; ask the key facts and `attach --given` (below) |

**Never search `/quobyte/proteomics-grp/coreomics/.submissions_db`.** It is a stale snapshot
(March 2026), and a free-text search for `0756` there matched an unrelated 2019 record.

**The record goes into the session once** (SKILL.md step 3b):
```
python3 scripts/submission_report.py attach --session <S> --record ~/core/PROT_0756
python3 scripts/submission_report.py attach --session <S> --given \
    '{"internal_id": "PROT_0756", "organism": "mouse", "prot_or_pep": "peptides", "sample_prep": "lab"}'
```
`input/submission.json` holds an **allowlisted** copy (identity, PI name/department/institution,
submitter name, submitted date, organism, UniProt, description, experiment type, proteins or
peptides, sample prep, buffer, beads, normalisation, analysis requested, the sample sheet) —
never an email, phone, payment, PPMS, contact or internal field; free text is scrubbed of
anything shaped like an email or phone number. `session.json` says `coreomics: {internal_id,
id, url, source}`. From there the report's **Submission** section (`make_analysis_html.py`),
the Methods' **Sample preparation** (`make_methods.py --submission`), the analysis brief
(`analysis_prompt.py --submission`), the session README and `record_run.py`'s PROT lookup all
read that one record. `--given` records are labelled "given by the user" everywhere.

`submission_report.py notes --session <S>` lists the **Data Quality Notes** it finds: organism
blank or disagreeing with the FASTA searched; UniProt blank (the Core chose the database); sheet
ids with no raw file, and raw files with no sheet id; conditions analysed that pool or split the
sheet's; beads blank on an affinity experiment; a form that contradicts itself about who
prepared the samples; and **pairing** — PROT_0756's names are `Old - JPH3 - Mouse 1`: every mouse
gave all five IPs (JPH3, JPH4, Kv2.1, RyR, IgG) and mice 1–3 are Old, 4–6 Young, so samples from
one mouse are not independent, and the note says whether the design analysed carries the mouse.

## 1. `fetch` — the submission, its neighbours, its shares (local)

```
python3 scripts/core_submission.py fetch 807 --out ~/core/PROT_0807
bash scripts/hive_exec.sh 'mkdir -p ~/core/PROT_0807'
bash scripts/hive_exec.sh --put ~/core/PROT_0807/hive/submission_summary.json '~/core/PROT_0807/'
bash scripts/hive_exec.sh --put ~/core/PROT_0807/hive/submission.json '~/core/PROT_0807/'
```

**The API (verified).** Base `https://ucdavis.coreomics.com/server/api`, header
`Authorization: Token <tok>`.

- **Look up by number with the exact filter:**
  `GET submissions/?lab=PROTEOMICS&internal_id=PROT_0807` returns exactly one record.
  **Never `search=`** — `search=0807` returns two. Numbers are `PROT_%04d` (PROT_0763 …
  PROT_0807 seen). The script accepts `PROT_0807`, `prot-807`, `0807`, `807`, `#807`, or the
  12-hex CoreOmics id (`GET submissions/<id>/`), and re-checks `internal_id` client-side.
- **List:** `GET submissions/?lab=PROTEOMICS&page=N&page_size=100&ordering=-submitted`
  (DRF `count/next/results`). List results **already include** `submission_data.samples`,
  which is what makes the neighbour scan cheap.
- **Fields used:** `id` (12-hex), `internal_id` (**can be null** — folders then use the hex
  id), `url`, `status`, `submitted` (e.g. `2026-09-10T15:02:53.863095-07:00`),
  `first_name/last_name/email` (submitter), `pi_first_name/pi_last_name/pi_email`,
  `pi.institution.name`, `pi.department`, `institute`, `contacts`, and
  `submission_data.{organism, description, sample_prep, data_analysis, proteomics_type,
  mass_spec_wanted, samples[{unique_id, sample_name, condition_name}]}`.
- `data_analysis` is either `"I want the proteomics core to do the data analysis"` or
  `"I only require raw data and will do my own data analysis"` → `raw_data_only`.
- `mass_spec_wanted` is `"timsTOF HT"`, `"Exploris 480"`, `"Fusion Lumos"`,
  `"No idea what this means, I'll let you decide what's best"`, or null.
- **Campus:** on campus only when `pi.institution.name` is exactly `"UC Davis"` (153 of the 300
  newest records); anything else is off campus. With no PI, `institute` is used; with
  neither, the summary says `campus_basis: missing` and warns.

**Neighbours.** `fetch` pages the list newest-first and keeps every submission within ±240
days (`--neighbor-days`), with its `unique_id`s, stopping once a page is older than the window.
`locate` needs them to tell *our* BN1 from another lab's, and refuses a `--max-days` wider
than the window recorded here.

**Outputs:** `submission.json` (raw record), `samples.tsv`, `submission_summary.json`
(identity, PI, contacts, campus, conditions present/groups/missing, `canonical_project_dir`,
`share_dir`, neighbours, existing Bioshare shares). A failure listing shares is recorded as
`existing_shares_error` and does not block. **Only `hive/` leaves this computer:**
`hive/submission_summary.json` is the summary with every email and the contacts removed and
free text scrubbed (`hive_summary()`), and `hive/submission.json` is the allowlisted record
(`submission_report.py`). The raw record holds emails, phones and PPMS/payment fields and the
full summary holds emails; both stay local, where `bioshare` and `email-draft` need them.
`stage`'s staff-facing `SUBMISSION.md` names people but shows no email.

**The organism is not resolved here.** `organism_as_submitted` is the submitter's free text,
labelled `confirmed: false`. Put it to the staff member as the proposed answer and confirm
it (golden rule #4), then `fetch_fasta.py resolve`.

## 2. `locate` — find each sample's raw file (HIVE)

```
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/core_submission.py locate \
    --summary ~/core/PROT_0807/submission_summary.json --out ~/core/PROT_0807'
bash scripts/hive_exec.sh --get '~/core/PROT_0807/sample_files.tsv' ~/core/PROT_0807/
bash scripts/hive_exec.sh --get '~/core/PROT_0807/locate.json' ~/core/PROT_0807/
```

**Where raw data lives.** Flinders root on HIVE is `/nfs/lssc0/flinders/proteomics` (NFS; on
a Mac the SMB mount is `/Volumes/proteomics/...`, on Windows `R:\...`). Raw files are in
`Data/raw_data/{tTOF_HT,Exploris480,Lumos1}/<month folder>/`. **Month folder names are
free-form** — `sep26`, `Aug26`, `JUL26`, `June26`, `jan25AndPM`, `Std_He_...` — so every
folder is scanned. About 20k `.raw`/`.d` entries in total; a listing is fine on the login
node. The requested instrument's folder is scanned first but never exclusively: staff
sometimes run a sample on another instrument (gate `instrument` warns). Compute nodes can
read this tree (the PROT_0793 search did).

**Filename conventions (measured on the full listing):**

| instrument | pattern | acquisition date |
|---|---|---|
| timsTOF | `MMDDYYYY__60SPD_DIA-<unique_id>_S3-<well>_1_<acq#>.d` | leading `MMDDYYYY` |
| timsTOF HT plate | `YYYYMMDD_<num>_...`, `MMDDYYYY_<num>_rerun_...` | leading 8 digits |
| Exploris | `Ex<MMDDYYYY>_<n>_<unique_id>.raw` (`Ex08312026_380_JE21.raw`) | `Ex` + 8 = MMDDYYYY |
| Exploris | `Ex<DDMMYY>_...` (`Ex040826_HeL50…` in `aug26`, `Ex100926_…` in `sep26`) | `Ex` + 6 = DDMMYY |
| Lumos | `FL<DDMMYY>_...` (`FL280826_HeL50…` in `aug26`) | `FL` + 6 = DDMMYY |
| Lumos, undated | `FLsep26_wa_20260909005054.raw`, `FL9sep_waID_3.raw` | modification time |

Only a real calendar date in [2015-01-01, today+1] counts. Typos exist — `070622026__`,
`062420266__` (nine digits) — and fall back to modification time, flagged `date_source:
mtime` and warned by the `mtime_dates` gate.

**Matching.** Each `unique_id` must appear as a **delimited, case-insensitive token**, with
`-`, `_` and space interchangeable on both sides: `EB_001` finds `EB-001`; `KG1` does **not**
find `KG13`. Four further rules, each from a failure found by review:

- **timsTOF names are matched only in their sample field** (`DIA-<unique_id>_S<n>-<well>`).
  Anywhere else, an id like `A1-1` or `S3-A1` matches the `_S3-A1_1_` plate position that
  every timsTOF run carries — and succeeds with another lab's files.
- **The longest id owns a run.** Every id in play — this submission's and every neighbour's
  — is tried against each candidate file, and the longest match wins: `DH1` does not take
  `DH1-1`'s run (gate `shadowed_by_longer_id` shows what was left to the longer id).
- **An id that names more than 10 raw files across the whole listing**, any instrument and
  any year, is weak: too common to pick out one submission's runs.
- **Only files acquired on or after the submission date and within `--max-days`** (default
  240, never more than fetch's neighbour window) count.

**Measured on the 45 newest submissions:** token matching found 1:1 files for many timsTOF
submissions (PROT_0772 4/4, 0771 12/12, 0768 10/10, 0764 6/6). It failed in four ways, and
each has a gate:

| failure mode seen | example | handling |
|---|---|---|
| short / well-like / numeric / too-common ids match hundreds–thousands of files | `A3`, `H10`, `001` | **weak id** (alnum < 3, all digits, `^[A-H](1-12)$`, or > 10 files overall): never auto-assigned → `weak_ids` FAIL |
| the same ids reused by another submission | `BN1–6` in 0794 **and** 0776; `SG001`; `GV1` | **ambiguous label**: a run is ambiguous when ANY other submission — older, same day or newer — uses the label and was submitted on or before the run. Only ambiguous runs → `ambiguous_label` FAIL; an unambiguous run also exists → it is chosen, `ambiguous_files_excluded` WARN |
| older runs of the same label | 2025 files for a 2026 submission | date window; counted in `out_of_window` INFO |
| several files per sample | re-injections | most recent chosen, the rest listed → `alternates` WARN |

**Labels alone cannot decide an ambiguous run, so staff do.** The `ambiguous_label` detail
names each other submission and whether it already has **earlier unambiguous runs** of that
label (runs acquired after it was submitted and before anyone else used the label). Measured
example: PROT_0776 and PROT_0794 both use BN1–6. For 0776, the Aug-24 BN runs are unambiguous
(0794 had not been submitted) and are chosen. For 0794, the Sep-8 BN runs are ambiguous —
0776 was submitted first — so `locate` fails, and the detail says 0776 already has its own
Aug-24 runs. That fact is what lets staff decide: re-run with `--accept-ambiguous` (recorded in
`locate.json` under `accepted`) if these runs are this submission's, or `--files-from` their
own list.

Also hard gates: `no_files`, and `duplicate_assignment` (one file claimed by two samples).
`unmatched_samples` is a hard gate unless staff confirm those samples were never run
(`--allow-partial`).

**HT plates.** If a file in the window is named `<8-digit date>_<num>_…` for this submission
(`_807_` or `_0807_`), `locate` exits **4**: use `ht_manifest.py`. The number is matched only
right after the leading date, because anywhere else it collides with run counters —
Exploris `Ex08312026_380_JE21.raw` would otherwise impersonate submission 0380.

**Outputs:** `files.txt` (chosen absolute paths), `sample_files.tsv` (`unique_id,
sample_name, condition_name, file, acquired, date_source, status, alternates, note`; status
`matched|weak|unmatched|ambiguous_label`, plus `unassigned` in `--files-from` mode), and
`locate.json` (every gate with detail, the file list, and the `accepted` flags). Proposal
files are written even on exit 2.

**Resolving a hard gate:** show the staff member the proposal in one message. When they pick
files by hand, write them to a list and re-run with `--files-from FILE`: matching is skipped,
paths are still checked for existence (`paths_exist`), and each file is mapped back to the
sample whose id is the longest one it carries — which makes even a weak id usable inside a
curated list.

## 3. `stage` — link the raw files into the service directory (HIVE)

```
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/core_submission.py stage \
    --summary ~/core/PROT_0807/submission_summary.json --files ~/core/PROT_0807/files.txt'   # dry run
# ... ONE staff confirmation of the file list AND this folder, then the same with --apply
```

**The service directory (verified)** is `<flinders>/Data/lab/service/{on_campus,off_campus}/`
— the active one, updated daily (newest change 2026-09-15). Layout:
`on_campus/<PI folder>/<project>` and `off_campus/<Institution folder>/[<PI/lab folder>/]<project>`.
Folder names are human and inconsistent: `Isseroff`, `McDonald karen`, `Serapio-Palacios-lab`,
`UCSF/Feeley_lab`, `Stanford/Dixon-lab`, `Meneses-Erica_National-Autonm-Uni-MX`, `Reckitt`.
Existing project folders hold **real raw copies** today (0 symlinks); this flow links instead.
`/quobyte/proteomics-grp/SERVICE/` is a small HIVE-side compute mirror, **not** the service
directory.

**Group folder resolution** (case-insensitive, accent-folded alphanumeric tokens). Each rule
exists because the looser version picked a wrong folder with exit 0:

- **The PI's surname:** all of its tokens must be in the folder name. Surname particles
  (`de`, `la`, `van`, `von`, `da`, `di`, `le` …) never match on their own: a surname that has
  them must appear whole (`delacruz`) — "de la Cruz" once matched `UC_Santa_Cruz`. Two-letter
  surnames (Li, Wu) match whole.
- **A folder naming somebody else is not a match.** Other personal-name tokens in the folder
  that are not the PI's first name make it ambiguous: PI Ying Wang is not `Wang Wei`.
  (`McDonald karen` for Karen McDonald is fine.)
- **on campus:** top-level folders. **off campus:** first- and second-level folders; a
  second-level match must sit under a folder that looks like the PI's institution.
- **off campus with no PI folder:** a top-level folder whose *distinctive* institution words
  are **exactly** the institution's (after dropping words like university/of/the/college/lab/
  uc/medical/center, and mapping UC campuses: "UC San Francisco" → `ucsf`, "UC Santa Barbara"
  → `ucsb`, "UC Berkeley" → `berkeley`, "UC Santa Cruz" → `ucsc`, San Diego → `ucsd`, Los
  Angeles → `ucla`, Irvine → `uci`, Riverside → `ucr`, Merced → `ucm`), then `<that>/<PI last>`;
  otherwise a new `<Institution>/<PI last>`. "University of Washington" is not
  `Washington_State_Univ`, and "UT Southwestern" is not `Texas_AM`.
- **One clean candidate** is reused; **none** proposes a new folder; **anything else** exits 2
  with the candidates and reasons — re-run with `--service-dir <folder>`, which must resolve
  under the service root.

**The chosen folder is part of the single staff confirmation**, together with the file list —
a new-folder proposal is often a spelling of a folder that already exists.

**The project** is `<group>/<internal_id>`. A `.core_submission.json` there naming a different
submission → exit 2; the same submission → idempotent re-run. `--apply` refuses a `files.txt`
whose `locate.json` (beside it) hard-failed, unless it came from `--files-from`. It creates:

- `raw/<basename>` — **relative** symlinks to the raw files (absolute only when a file is
  outside the Flinders root, with a warning that it will not be visible). A real file of the
  same name is never overwritten; a link pointing elsewhere is replaced and reported. Every
  link is checked to resolve.
- `SUBMISSION.md` — staff-facing: CoreOmics link, PI and institution, submitter, dates,
  organism *as submitted (unconfirmed)*, experiment type, analysis requested, the sample table
  with raw file and acquisition date, the work dir and the server-side share dir. Never
  delivered. It starts with a generated-by marker; a `SUBMISSION.md` without it was written by
  hand and is left alone — the record goes to `SUBMISSION.core.md` instead.
- `.core_submission.json` — the machine record (ids, files, service/work/share dirs, skill
  version when resolvable), also copied into the work dir so `deliver` can find the project.
- the compute **work dir** `CORE_WORK_ROOT/<campus>/<group path>/<internal_id>/`, where the
  session of record lives.

Folders it creates are `g+rwxs,o+rx` (the tree is shared; umask 077 would lock colleagues out).

## Symlink visibility — measured, and the rules that follow

| link | on HIVE | over SMB (Mac/Windows) | Bioshare |
|---|---|---|---|
| **relative**, source and target both under the Flinders root (CoreOmics `views/`, e.g. `../../../../../../projects/2019/02/<id>`) | works | **visible and followable** — verified for a relative link to a `.raw` file and to a `.d` folder: both visible and readable | untested for `Data/raw_data` (see step 6) |
| **absolute** `/nfs/lssc0/...` | works | **invisible** — PROT_0793 `share/raw`: 210 links on HIVE, 0 entries on the Mac mount; an absolute link placed beside the relative ones above was invisible too | — |
| **into `/quobyte/...`** | works | **invisible** | **works only if Bioshare's file-streaming allows `/quobyte`** — PROT_0793 `share/search`: 59 links → `/quobyte/proteomics-grp/brett/PROT_0793` served nothing because Bioshare streams files through an Apache module whose allowed-folder list had `/quobyte` for http (port 80) but not https (443) — fixed 2026-09-16 (Adam Schaal) |

Bioshare (`bioshareX/sendfile.py`) also serves a file only if its **realpath** is under the
server's `DIRECTORY_WHITELIST`, and `check_symlinks_dfs` raises `IllegalPathException` for a
symlink whose target is outside it — on the older code path that **locks the share**. Which
paths are allowed is server configuration, and it can differ between http and https: Bioshare
hands downloads to an Apache module with its own allowed-folder list, which had `/quobyte` for
port 80 but not 443 until 2026-09-16. Brett reports PROT_0793's absolute `Data/raw_data` links
download.

So: **links between two Flinders paths are relative. Deliverables in a share are real files,
so a delivery doesn't depend on a server setting. Never create a link from a Flinders directory
into /quobyte — and never "repair" one that is already there; `deliver` leaves it alone.**

## 4. `conditions` — the design from CoreOmics

```
python3 scripts/core_submission.py conditions --summary ~/core/PROT_0807/submission_summary.json \
    --sample-files ~/core/PROT_0807/sample_files.tsv --out ~/core/PROT_0807/conditions.csv
bash scripts/hive_exec.sh --put ~/core/PROT_0807/conditions.csv "$S/input/"
```

Builds `{raw run name: condition_name}` for located samples and calls `collect_conditions.py
--map --runs --mapping-json`, so `File.Name` has one definition (architectural rule #3). It
sets `needs_user_input` (exit 2) when a located sample has a blank condition, conditions
differ only in case or punctuation (`Control` / `control`), every sample shares one
condition, every condition is unique (sample names typed into the condition column), any
group has one sample, or `collect_conditions.py` reports ambiguities. The output's
`questions` are exactly what to ask — nothing more. Otherwise the conditions go into the one
compute confirmation as they are.

## 5. The search and DE — session on HIVE, report written locally

**"I only require raw data" submissions get no search and no DE** unless staff explicitly ask:
go straight to step 6 with `--mode raw-only`.

Otherwise run SKILL.md steps 2–9 with the **session of record on HIVE**, under the work dir:
```
bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/session.py init \
    --name PROT_0807 --base <work_dir> --raw $(cat ~/core/PROT_0807/files.txt)'
S=<the session dir it printed>
```
Search, DE, figures, audit, `sample_quality.py`, `make_methods.py` (it reads raw metadata) and
provenance all run there — heavy steps as SLURM jobs. Writing the report needs the data
locally, so pull a working copy:
```
mkdir -p ~/core/PROT_0807/session/output/search
for d in tables figures; do bash scripts/hive_exec.sh --get "$S/output/$d" ~/core/PROT_0807/session/output/; done
bash scripts/hive_exec.sh --get "$S/output/search/report.parquet" ~/core/PROT_0807/session/output/search/
for f in AUDIT.md SAMPLE_QUALITY.md methods.md; do bash scripts/hive_exec.sh --get "$S/output/$f" ~/core/PROT_0807/session/output/; done
```
Attach the submission to `$S` right after `init` (`submission_report.py attach --session "$S"
--record ~/core/PROT_0807`, on HIVE: the allowlisted `hive/submission.json` put there in step 1), and pull
`$S/session.json` and `$S/input/{submission.json,samples.tsv,raw_files.txt,conditions.csv,
search.fasta.meta.json}` with the rest. Then locally: write `AI_Analysis_Report.md`;
`make_analysis_html.py --session ~/core/PROT_0807/session --out
~/core/PROT_0807/session/output/Analysis_Report.html` (its Submission section comes from the
attached record, never an email address); `to_docx.py` for the report and methods. **Push the finished files back before delivering** — `deliver` copies from `$S`:
```
for f in AI_Analysis_Report.md AI_Analysis_Report.docx Analysis_Report.html methods.md methods.docx; do
  bash scripts/hive_exec.sh --put ~/core/PROT_0807/session/output/$f "$S/output/"; done
```
**Step 8d's expert review always runs**: the result goes to a collaborator.

## 6. `deliver` — real files into the Bioshare folder (HIVE)

```
python3 ~/proteomics-pipeline/scripts/core_submission.py deliver \
    --summary ~/core/PROT_0807/submission_summary.json --session "$S"            # dry run first
# ... then the same with --apply
python3 ~/proteomics-pipeline/scripts/core_submission.py deliver \
    --summary ~/core/PROT_0807/submission_summary.json --mode raw-only          # "raw data only"
```

**The CoreOmics project tree (verified)** is `<flinders>/coreomics/projects/<YYYY>/<MM>/<id>/share/`,
with YYYY and MM the **literal** first two `-`-separated fields of `submitted` (coreomics_fs
`build_views.ensure_canonical`: `y, m = submitted.split()[0].split("-")[:2]`) — no timezone
conversion, so `2026-08-31T23:30:00-07:00` is `2026/08`. **Verified live:** all 6 CoreOmics
submissions made on the last evening of a month (already the next month in UTC) are filed under
their LOCAL, literal month. Directories are `amschaal:proteomics-grp 2775`. The builder had
**not** created `2026/09/` by 2026-09-16, so `deliver --apply` creates `<id>/share` itself —
and never writes `<id>/.submission/`, which the builder owns. `/quobyte/proteomics-grp/coreomics/`
is a **stale** copy (stops at 2026/03); never use it. Existing shares (4 of the 200 newest
submissions) have
`link_to_path = /nfs/lssc0/flinders/proteomics/coreomics/projects/2026/08/<id>/share`.

**Before writing anything, `deliver` refuses (exit 2)** when:

- **the session is not this submission's.** It must sit under the staged `work_dir`, or its
  `input/raw_files.txt` must list exactly the staged files. A stage record for a DIFFERENT
  submission above the session is refused outright. With no stage record at all, only
  `--force`, after staff confirm (recorded as `session_ownership: NOT CHECKED`). A PROT_0806
  session once delivered into PROT_0807's share with `verified: true`.
- **the delivery folder already has files.** Two deliveries mixed in one folder leave stale
  results beside new ones — use a new `--label`. (A size-guard job resuming its own interrupted
  folder is recognised by the `.core_delivery.json` it left.)
- **the path passes through a symlink or leaves the Flinders root** — e.g. `<id>/share` itself
  linked into /quobyte would have put the delivery there, verified.

**Analysis mode** (`--mode analysis`, the default unless CoreOmics says raw data only) fills
`<share>/<internal_id>_analysis_<YYYY-MM-DD>/` (`--label` replaces the suffix), copied
**dereferenced**:

- `Analysis_Report.html` — **required**; exit 2 without it
- `AI_Analysis_Report.docx/.md`, `methods.md/.docx`, `OUTPUT_FILES.md`, `AUDIT.md`,
  `SAMPLE_QUALITY.md`
- `tables/`, `figures/`, `reproducibility/`
- from `search/`: `report.parquet`, `report.pg_matrix.tsv`, `report.pr_matrix.tsv`,
  `report.gg_matrix.tsv`, `report.unique_genes_matrix.tsv`, `report.stats.tsv`,
  `report.log.txt`, `search_provenance.json`
- never: XIC folders, `.quant`, `.speclib`, temp folders, `output/raw_data`

**Raw-only mode** (`--mode raw-only`, automatic for "I only require raw data") needs no session
and no report: `<share>/<internal_id>_raw_data_<date>/` holds a README that says raw files
only, the MANIFEST, checksums, and `methods.md` when a `--session` with one is given; the raw
files are relative links in `<share>/raw/`.

Every delivery also gets:

- **`MANIFEST.txt`** — `[OK] <name>` for every item copied and `[SKIPPED] <name> -- <reason>` for
  every item missing, unreadable (including an unreadable subfolder), excluded, or not reached
  because of an error (architectural rule #4).
- **`README.md`** for the collaborator — every claim traced to a delivered file: "searched" only
  with `search/` output, "compared between sample groups" only with `DE_*.csv` tables. No
  internal paths, no email addresses.
- **`checksums.sha256`**.
- permissions `g+rw,o+r` on files and `g+rwxs,o+rx` on folders — Bioshare reads as another
  user, and a copied 0600 file would be unservable.

Writes never follow a symlink on the destination side (a planted `tables -> /elsewhere` link
once let a copy overwrite a file outside the share).

**Raw data** (`raw-only`, or `--include-raw yes`, or `auto` with "raw data only") goes in as
relative links in `<share>/raw/` to the staged files. The service project gets a relative
`results_<date>` link to the delivery (a failure there is `[SKIPPED]`, not fatal).

> ⚠ **Check once: relative `raw/` links in Bioshare.** Brett reports PROT_0793's absolute raw
> links download from Bioshare, and a relative link resolves to the same files, but nobody has
> opened a relative one in Bioshare yet. For the first raw-link delivery, **open the Bioshare
> link yourself and download one raw file before `bioshare send`**. `deliver` prints this
> warning and records `raw_whitelist_unverified: true`.

**Size guard:** above `--max-gb` (default 5) it writes `deliver_job.sh` and exits **5** —
`sbatch` it; the login node is for small copies only (golden rule #3).

**Verification on `--apply` walks the WHOLE share:** a symlink inside this delivery, or a
`raw/` link this run created that is absolute, broken or outside the Flinders root, fails.
Links already in the share that this run did not create (PROT_0793's `search/ -> /quobyte`)
are listed as warnings and left untouched. Inside the delivery folder every file must be one this
run delivered and readable by group and other. An error part-way still writes the MANIFEST,
checksums and `delivery.json` with `verified: false`. **Exit 2 means do not share.**
`delivery.json` (in the session dir; in the work dir for raw-only) records the folder, mode,
file count, bytes, skipped items, raw links, the server-side `share_dir`, and `verified`.

## 7. `bioshare` — register and share (local)

Through CoreOmics' Bioshare plugin (source: `amschaal/bioshare_coreomics_plugin`,
`amschaal/coreomics_fs`). Every call uses the summary's server-side `share_dir`.

```
bash scripts/hive_exec.sh --get '<delivery_json printed by deliver>' ~/core/PROT_0807/delivery.json
python3 scripts/core_submission.py bioshare status --summary ~/core/PROT_0807/submission_summary.json
python3 scripts/core_submission.py bioshare ensure --summary ~/core/PROT_0807/submission_summary.json [--apply]
python3 scripts/core_submission.py bioshare send   --summary ~/core/PROT_0807/submission_summary.json \
    --delivery ~/core/PROT_0807/delivery.json [--apply] [--email]
```

- **`status`** — `GET plugins/bioshare/submissions/<id>/submission_shares/` (verified with a
  staff token; paginated or a bare list, both handled). Keys: `id, submission, bioshare_id,
  name, notes, sub_folder, link_to_path, url`. The share whose `link_to_path` equals this
  submission's share dir is marked linked (verified live on PROT_0793's real share). (`GET
  plugins/bioshare/shares/?lab_id=PROTEOMICS` returns **403** for lab members — not used.)
- **`ensure`** — if nothing is linked, `POST` the same endpoint with
  `{"submission": <id>, "name": ..., "notes": ..., "link_to_path": <share_dir>}`, exactly what
  coreomics_fs `api.py:create_submission_share` sends. The name is the plugin model's default,
  `"<PI last>, <PI first>: <internal_id>"`, reduced to the serializer's
  `^[\w\d\s'".!?\-:,]+$`. Run `deliver --apply` first: the directory must exist. The response
  must be a share record with an `id`. *The create call has not been exercised live.* Errors come
  back as DRF JSON, e.g. `{"link_to_path": ["Path not allowed."]}`, and are printed with exit 3.
- **`send`** — `POST .../submission_shares/<pk>/share/` with `{"email": true|false}`, where `<pk>`
  is the share record's `id`. It grants **view + download** to the submitter, the PI and the
  submission's contacts; the dry run lists all three. **Outward-facing: only after the staff
  member says yes to "Share with <submitter>, <PI> and the contacts now?"** — a second yes,
  separate from the compute confirmation. `--apply` requires `--delivery delivery.json` with
  `verified: true`, this `internal_id` and this `share_dir`. It sends `email: false` unless
  `--email` is passed. *Its success response shape is unverified*; a list or a listing-shaped
  reply exits 3.
- **A redirect on any write is an error (exit 3).** urllib silently turns a redirected POST into
  a GET — `send --apply` once reported `applied: true` with the share LISTING as its response.

## 8. `email-draft`

```
python3 scripts/core_submission.py email-draft --summary ~/core/PROT_0807/submission_summary.json \
    --delivery ~/core/PROT_0807/delivery.json --share-url <url> --out ~/core/PROT_0807/EMAIL_DRAFT.md
```

A plain, friendly draft to the submitter (cc the PI): the Bioshare link, "start with
Analysis_Report.html" (or, for raw-only, where the raw files are), what is in the folder, the
acknowledgment request. Numbers appear only when they are in its inputs. It never sends; exit 2
means a placeholder is left (no link, or no submitter address).

## Troubleshooting by exit code

| exit | from | usual cause | fix |
|---|---|---|---|
| 2 | `fetch` | no such number; garbage id | check the number; pass the 12-hex id |
| 2 | `locate` | `weak_ids`, `ambiguous_label`, `duplicate_assignment`, `unmatched_samples`, `no_files`; `--max-days` wider than the neighbour window | show the proposal; `--files-from` a corrected list, `--accept-ambiguous` or `--allow-partial` as staff decide; re-fetch with `--neighbor-days` |
| 2 | `stage` | several candidate folders, or one naming someone else; project owned by another submission; `--apply` on a hard-failed locate | `--service-dir <folder>`; resolve the locate gates |
| 2 | `conditions` | blank / case-variant / single / all-unique / singleton conditions | ask exactly the `questions` |
| 2 | `deliver` | no `Analysis_Report.html`; session not this submission's; folder not empty; symlink in the path or the share; raw requested but no staged project; any verification failure | finish step 9 and push it; use the right session (or `--force` with staff); `--label`; remove the offending link; `stage --apply` |
| 2 | `bioshare send` | no share linked; no or unverified `delivery.json`, or one for another share | `bioshare ensure --apply`; deliver again and `--get` the new `delivery.json` |
| 3 | `fetch`, `bioshare` | no/invalid token; CoreOmics down; DRF validation error; a redirect or unexpected reply to a write | fix `~/.coreomics_token`; read the detail |
| 3 | `locate`, `stage`, `deliver`, `conditions` | not on a machine with the Flinders tree; incomplete scripts directory | run through `hive_exec.sh`; sync the whole `scripts/` |
| 4 | `locate` | the submission is an HT plate | `ht_manifest.py` (step 1a) |
| 5 | `deliver` | more than `--max-gb` to copy | `sbatch deliver_job.sh` |

## The PROT_0793 lesson

PROT_0793 was searched and "shared" by linking: `share/raw` held 210 absolute `/nfs/...` links
and `share/search` held 59 links into `/quobyte/proteomics-grp/brett/PROT_0793`. On HIVE the
share looked complete. Over SMB it showed **0** raw entries. The `search/` links served the
collaborator nothing, because Bioshare streams files through an Apache module whose allowed-folder list had `/quobyte` for http (port 80) but not https (443) — fixed 2026-09-16 (Adam Schaal). Nothing about the share was wrong on
HIVE, and nothing reported an error.
That is why `deliver` copies real files, makes only relative Flinders-to-Flinders links, and
refuses to call a delivery done until it has walked the whole share and found no other symlink.
