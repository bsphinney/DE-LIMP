# Access & where the skill runs (HIVE runbook)

**Model: Claude Code runs LOCALLY on the user's computer.** When a step should run on
HIVE, the local Claude Code drives HIVE **over SSH using the user's private key**
(`ssh -i <key> <user>@hive.hpc.ucdavis.edu '<command>'`, wrapped by `hive_exec.sh`).
The skill is never installed *on* HIVE for the user; it submits work there.

**Why local (not Claude Code on HIVE):** running Claude Code on a HIVE interactive
node ties up that node and the session **times out** mid-analysis. Keeping Claude Code
on the laptop means it stays alive while the actual compute runs as detached **SLURM
jobs** on HIVE — nothing depends on an interactive session staying open.

## Step 0a asks two questions → pick a mode
1. **"Do you have access to UC Davis HIVE (account + SSH private key)?"**
   If yes, **ask where the private key is** (e.g. `~/.ssh/id_ed25519`).
2. **"Are you a member of the UC Davis Proteomics Core?"**

Verify before trusting the answers:
```
bash scripts/check_access.sh <hive_user> <private_key_path>
```
Reads `recommended_mode` + `core_member` (= `facility_software_available`). If
`hive_ssh` is `"failed"`, `hive_ssh_error.kind` says why — relay it, don't guess:

| `kind` | Meaning → what to tell the user |
|---|---|
| `host_key` | First contact with HIVE, or its key changed. Compare `hive_host_key_fingerprints` with "HIVE host key" below. |
| `permission_denied` | Wrong username or key. On Windows, never the `DOMAIN+user` login name (see "Windows"). |
| `timeout` | VPN off / off-campus, or HIVE's connection throttle after many quick connections — wait ~20 min, then retry once. |
| `other` | Read `hive_ssh_error.detail` (e.g. a DNS failure = no VPN). |

Then:

| HIVE access | Core member | Mode | What happens |
|---|---|---|---|
| no | no | **local** | `setup.sh` installs the toolchain on the user's machine; public engines (DIA-NN Academia, Sage). |
| **yes** | **yes** | **hive_remote** | Drive HIVE over SSH. **Reuse the Core software already installed** in `/quobyte/proteomics-grp` (DIA-NN builds, pre-staged FASTAs). Search runs as a SLURM job. |
| **yes** | no | **hive_remote** | Drive HIVE over SSH, but you **rebuild the toolchain in your own HIVE home** (you can't read the Core group dir). See "Rebuild on HIVE" below. |
| no | yes | **local** | The Core software is on HIVE; without HIVE access you can't reach it → run locally with public engines and tell the user to request a HIVE account. |

## hive_remote runbook (Claude Code local → HIVE over SSH)
Set the connection once, then everything HIVE-side goes through `hive_exec.sh`:
```
export HIVE_USER=<hive_user>
export HIVE_KEY=<private_key_path>        # the path the user gave you
bash scripts/hive_exec.sh 'hostname; sbatch --version | head -1'   # confirm
```
1. **Put the skill on HIVE, and put it again after every skill update** (the scripts run
   there, and every record they write names the skill version from `.claude-plugin/
   plugin.json`, which is copied with them):
   ```
   bash scripts/hive_exec.sh --put-skill     # scripts/ + .claude-plugin/ -> ~/proteomics-pipeline/
   ```
   That is the same as `--put ./scripts '~/proteomics-pipeline/'` plus `--put ./.claude-plugin
   '~/proteomics-pipeline/'`. A `scripts/` folder put up without `.claude-plugin/` records the
   version as unknown (the run log, provenance, the deposit package), and one not re-put after
   an update runs the old scripts. Every session checks this first (SKILL.md step 0; "Keeping
   the skill current" below).
2. **Toolchain on HIVE:**
   - **Core member:** `acquire_tools.sh` (run on HIVE) finds the group's DIA-NN
     builds (`/quobyte/proteomics-grp/dia-nn/build_*/diann-<version>/`, plus an older
     `.sif`); `fetch_fasta.py --hive`
     reuses `/quobyte/proteomics-grp/MRS/` FASTAs. Build the R/Python/DE env once with
     `setup.sh` (it's the same micromamba env):
     ```
     bash scripts/hive_exec.sh 'bash ~/proteomics-pipeline/scripts/setup.sh'
     ```
   - **Non-Core HIVE user:** same `setup.sh`, plus you must acquire DIA-NN/Sage
     yourself (no access to the group's builds). See "Rebuild on HIVE".
3. **Stage the raw data** — usually there is nothing to stage. See "Data already on a
   network drive" below: run `hive_path.sh` on the raw folder FIRST, and use the HIVE
   path it verifies in place. Upload only when the data is on the laptop's own disk:
   ```
   bash scripts/hive_path.sh '<the path the user gave>'     # exit 3 = local disk
   bash scripts/hive_exec.sh --put /path/to/raw '~/proteomics-pipeline/data/'
   ```
4. **Run the search as SLURM jobs** (never the login node). Run `run_search.py --sbatch
   job.sh` ON HIVE, then submit exactly what it prints under "submit with:":
   ```
   bash scripts/hive_exec.sh 'cd ~/proteomics-pipeline && python3 scripts/run_search.py ... --sbatch job.sh'
   bash scripts/hive_exec.sh 'bash <out>/submit.sh'     # DIA-NN from a FASTA, or the 5-step chain
   bash scripts/hive_exec.sh 'cd ~/proteomics-pipeline && sbatch job.sh'   # only if job.sh was written
   bash scripts/watch_run.sh --all <out> --hive          # submit.sh wrote <out>/jobs.txt
   ```
   A DIA-NN search that predicts its library from the FASTA is two jobs (library, then
   search) chained by `<out>/submit.sh`, and writes no `job.sh`; a search routed to the
   5-step chain (>5 files) also writes none, and exits 3.
5. **DE + figures + report:** run on HIVE (`run_de.R`, `make_figures.R`, …) or pull
   `report.parquet` back and run them locally (DE/figures are light).
6. **Retrieve results** into the session folder on the user's machine:
   ```
   bash scripts/hive_exec.sh --get '~/proteomics-pipeline/out' ./<session>/output/
   ```
   The search's per-run `.quant` files sit in `<out>/quant` (single-shot) or
   `<out>/quant_step2`/`quant_step4` (5-step chain) — kept there, not beside the raw data on
   the instrument share, and they count toward your HIVE home quota. They are only needed to
   re-run the search; for DE and the report, fetch the report files instead of the whole
   folder (e.g. `--get '<out>/report.parquet' ./<session>/output/`).

## Keeping the skill current (hive_remote)
Three copies of the skill must agree: the plugin on the laptop (it runs SKILL.md), its copy
in `~/proteomics-pipeline/` on HIVE (it runs the scripts), and the Core's current release. On
2026-09-29 two Core staff were on 2.6.0 in both places while 2.8.0 was out, so none of the
2.7/2.8 fixes had run for them; one HIVE copy had no `.claude-plugin/`, so its `sdrf.tsv` said
"v0.0.0". Nothing compared the copies.

`bash scripts/skill_version.sh --check-hive --mode hive_remote` (SKILL.md step 0) compares all
three in one SSH call and prints one JSON object: `status`, `local`, `hive`, `release`,
`release_state`, `behind_release`, `jobs`, `next`, `say`.

| Exit | `status` | What to do |
|---|---|---|
| 0 | `in_step` | Nothing, unless `say` is set: this computer is behind the release. Relay it. (`release_state` `untrusted`: the release file is not an admin's, and is ignored.) |
| 0 | `skipped` | Nothing, and say nothing: no HIVE login, and no `--mode hive_remote` (local mode, no HIVE). |
| 3 | `hive_missing`, `hive_behind`, `hive_differs` | Run `next` (`hive_exec.sh --put-skill`), then the check again. |
| 4 | `hive_ahead` | Stop. Relay `say`: update the plugin here, or, for an unreleased build on HIVE, ask the user which to keep. |
| 5 | `no_login`, `unreachable`, `local_unknown` | One line saying the check could not run; carry on. `no_login`: run `check_access.sh`, which saves the login. |
| 6 | as 3 | `jobs` of the user's `jobs_total` SLURM jobs use the copy, or squeue answered with an error: do not put now. Ask the user; put after they finish. |

- **The login.** The check needs the HIVE login in the environment or in
  `~/.config/ucdavis-proteomics/hive.env`, and the environment does not survive from one tool
  call to the next. `check_access.sh` writes the file once SSH works (user and key path, mode
  600; `hive_login_saved` in its JSON). With `--mode hive_remote` a missing login is exit 5
  (`no_login`), never a silent skip.
- **The version** is `.claude-plugin/plugin.json`'s, read the same way as everywhere else:
  `skill_version.sh` is the bash mirror of `skill_version.py` (tests keep them equal), and the
  HIVE copy is read by that same function, sent along with the command. A HIVE copy too old to
  have `skill_version.sh` is still read.
- **Which is newer** is decided on the release number first, without a `v` or a `-dev`/`+build`
  suffix: HIVE's `2.10.0` is newer than this computer's `2.9.0-dev` (exit 4, never a
  downgrade). Only when the release numbers are equal does a different string, or different
  files, count as `hive_differs`, and this computer's copy wins, because its SKILL.md is the one
  running. **Different files** means a `cksum` digest over `plugin.json` and every file in this
  computer's `scripts/`, taken the same way on both sides. Files an older put left on HIVE, which
  this computer no longer has, are not counted (`--put-skill` never deletes).
- **Running jobs.** Before an automatic put, the check lists this user's SLURM jobs
  (`squeue -u <user> -o %o`, the batch script of each) and counts the ones whose script names the
  copy's `scripts/`: `$HOME/proteomics-pipeline/scripts/` expanded, its `pwd -P` form, or
  literally `~/`, `$HOME/` or `${HOME}/proteomics-pipeline/scripts/`. A job's end-of-job steps
  (run log, Slack, FRAN) run those scripts hours after it started, so a copy is never swapped
  under them (exit 6). The bare word `proteomics-pipeline` is not enough: it also matches
  `~/.proteomics-pipeline` (the toolchain) and `~/proteomics-pipeline-280-*`. On 2026-09-29, 1 of
  Brett's 225 jobs used the copy (a PROT_0002 step, through its job-end hook), so counting every
  job would have stopped every session with a false message.
  - A script that cannot be read counts: fail closed. A job with no script path (an
    interactive `srun`) does not count.
  - squeue present but answering with an error ("Socket timed out") is exit 6 too, with
    `jobs: null` and the error in `say`. Only a HIVE with no squeue at all puts without
    counting.
  - `jobs_total` is every job of the user, for the message.
- **Updating the plugin** (checked on Claude Code 2.1.285): `/plugin` → **Installed** →
  **ucdavis-proteomics-core-pipeline** → **Update now**, or in a terminal
  `claude plugin update ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core`, then
  `/reload-plugins`. `/plugin marketplace update ucdavis-proteomics-core` only refreshes the
  catalogue: the installed plugin stays where it was. Automatic updates are **off** for this
  marketplace until the user turns them on: `/plugin` → **Marketplaces** →
  **ucdavis-proteomics-core** → **Enable auto-update**. That default is the likely reason staff
  sat on 2.6.0. Claude Desktop: **+** → **Plugins** → **Manage plugins**.

Step 0 also reads the Core's notes, in `/quobyte/proteomics-grp/skill_notes/` (`references/notes.md`).

### The Core's current release: `CURRENT_VERSION`
`/quobyte/proteomics-grp/skill_release/CURRENT_VERSION`, beside `skill_runs/` and
`skill_issues/` (`SKILL_RELEASE_DIR` overrides the folder). Line 1 is the version; the lines
after it are `#` comments saying who wrote it and when. Only the Core group can read it (2770),
so for anyone else `release_state` is `no_access` and nothing is said.

**Who writes it: a maintainer, as the last step of a release** — after the new version is
merged to GitHub `main`, which is where `claude plugin update` gets it:
```
bash scripts/skill_version.sh --publish-release     # from the released copy, with a HIVE login
```
It writes this copy's version (numbers and dots only, a leading `v` dropped, never a test build
or an unknown one), through `mktemp` in that folder and a rename. From a git checkout it also
refuses a version that `origin/main`, as last fetched, does not ship (no network is used:
`git fetch origin` first). It prints `main_checked: false` when this copy is not a checkout,
such as the installed plugin. It refuses to move the release backwards unless given
`--allow-older`, which is for withdrawing a bad release.

The reader takes the first line that is neither blank nor a `#` comment, and reads `v2.9.0` as
`2.9.0`; anything else is `release_state: unreadable`, and nothing is said.

**Only an admin's file is trusted.** `/quobyte/proteomics-grp` is group-writable and not sticky,
so any group member can rename or replace a folder at its top, `skill_release/` included. The
check trusts `CURRENT_VERSION` only while the file AND its folder are owned by a Core admin.
Anything else reads `release_state: untrusted`: no release, no warning, as if nothing were
published.
- **The admins** are listed in `scripts/core_admins.txt`, one HIVE username per line (`#`
  comments and blank lines allowed), currently `brettsp`. It is the one list, with two
  readers kept equal by `tests/test_skill_version.py`: `skill_version.sh`'s
  `skill_core_admins` (bash) and `skill_version.py`'s `core_admins()` (Python: `notes.py`,
  `slack_collab.py`).
- **Ownership** is read on HIVE with GNU `stat -c %U`, falling back to BSD `stat -f %Su`. It
  reads the link itself, never its target, so a member's symlink to an admin's file stays the
  member's.
- **`--publish-release` refuses** unless the HIVE user running it is an admin, or when
  `skill_release/` exists but belongs to someone else. Move that folder aside first; the
  publish then makes it anew, owned by the admin.

`hive_exec.sh --put-skill` never writes it. Everyone runs `--put-skill`, maintainers included,
and a maintainer's HIVE copy often holds an unreleased build. Publishing from there would tell
staff to update to a version the marketplace does not have yet.

## Data already on a network drive (check BEFORE any upload)
Core data lives on shares that HIVE mounts. A path the user gives from a mapped drive
(`T:\Data\lab\service\...`) or a Mac SMB mount (`/Volumes/proteomics/...`) is usually
already on HIVE under a different name. Real case (2026-09-23): `T:` was
`\\128.120.208.24\proteomics`, HIVE's `/nfs/lssc0/flinders/proteomics`; nothing mapped
the two, HIVE's `/quobyte` was searched instead, and 7.1 GB that was already on HIVE was
uploaded (~20 min).
```
bash scripts/hive_path.sh 'T:\Data\lab\service\PROT_0807'
```
It resolves the drive locally (`net use`, else PowerShell `Get-PSDrive` on Windows;
`mount` on macOS/Linux), maps the share, and checks the candidates on HIVE in one SSH
call:

| Share (server) | HIVE path |
|---|---|
| `proteomics` (`128.120.208.24`, the Flinders share) | `/nfs/lssc0/flinders/proteomics` |
| `proteomics-grp` (the Core group dir) | `/quobyte/proteomics-grp` |
| any other `<share>` | tries `/nfs/lssc0/flinders/<share>`, then `/quobyte/<share>` |

- **`verified: true` (exit 0)** → use `hive_path` in place in the search. Upload nothing.
  "Verified" means HIVE has the same file size, or every top-level name of the folder
  (a case difference like `t:\data` vs `Data` is resolved).
- **`verified: false` (exit 1)** → it IS on a share, but no candidate matched. Do not
  search HIVE for it (the Flinders mount is NFS, and walking it takes forever). Show
  `candidates` + `how` and **ask the user** where that drive lives on HIVE, then re-run
  `hive_path.sh` or check the path they give with `hive_exec.sh 'ls <path>'`.
- **exit 3** → the laptop's own disk; uploading with `--put` is the only way.

`hive_exec.sh --put` runs the same check and **refuses** a source it verifies on HIVE,
printing the path to use (override: `HIVE_PUT_FORCE=1`). For an ordinary local path the
check makes no SSH call.

## Windows (Git Bash)
`hive_remote` from Windows needs only Git Bash + the private key. Nothing runs locally:
- **Every script runs on HIVE** through `hive_exec.sh` — including the Python ones. On
  Windows `python3` is often the Microsoft Store stub (`...\WindowsApps\python3`), which
  opens the Store instead of running; `check_access.sh` reports it as
  `local_python3.usable: false`. Never run a skill script locally with it. The one
  exception is Core staff's CoreOmics lookups (`core_submission.py check` / `identify` /
  `fetch` / `bioshare`): the key stays on this computer, so they need a real local Python 3
  (section "Core staff" below).
- **No rsync** in Git Bash: `--put/--get` fall back to `scp` on their own. A trailing
  slash on the source (`dir/`, rsync's "copy the contents") is refused — drop the slash.
- **Username:** ssh's default user on Windows is the AD login, e.g. `AD3+gabrig`, and
  HIVE answers "Permission denied". HIVE usernames are the plain UC Davis id
  (`gabrig`). Always set `HIVE_USER`; `hive_exec.sh` and `check_access.sh` refuse a name
  containing `+`, `\` or `@`.
- **Mapped drives** (`T:`, `R:`) are Flinders or quobyte shares — run `hive_path.sh`
  before any upload (section above).
- ControlMaster connection reuse is off on Windows (unsupported), so each call is a new
  SSH connection — batch commands into one `hive_exec.sh` call where possible.

## HIVE host key (first contact)
On a machine that has never connected, the first non-interactive SSH call fails with
"Host key verification failed" (nothing can answer the yes/no prompt). The scripts use
`StrictHostKeyChecking=accept-new`: an unknown key is recorded, a CHANGED key is still
refused. `check_access.sh` reports `hive_host_key_known`, the fingerprints, and
`hive_host_key_matches_published`. HIVE's keys as of 2026-09-23 (`ssh-keyscan`):
```
ED25519  SHA256:b5nv86Ciaqg1yrUVai6bZ0Hk4IpzAFLWtIPDBdacbQM
ECDSA    SHA256:AmJ+z2miIMXlSAcm7k8YwKlIWk5+VqxyT0R4q3fvpcA
```
If an unknown key does NOT match these, `check_access.sh` does not auto-accept it
(`hive_ssh_error.kind: host_key`). Stop and have the user confirm the key with UC Davis
HPC support before connecting — it may have been rotated, or someone is in the way.

## Core staff — submission workflow prerequisites
Searching and delivering a **CoreOmics submission** ("search the data from submission 807",
SKILL.md steps 1c + 12b) splits across two machines, so it needs three things in place:
- **Your own CoreOmics API key, saved on your own computer.** `check`, `fetch` and `bioshare`
  run locally; HIVE has no CoreOmics key. Steps:
  1. In CoreOmics (https://ucdavis.coreomics.com): **Profile → API Key**, click to show it,
     **Create** if there is none. Copy it.
  2. In your own terminal, never in the chat, run this ONE line. It then waits, showing
     nothing, while you paste the key (in Git Bash: Shift+Insert); then press Enter:
     ```
     read -rs TOK && [ -n "$TOK" ] && F="$(cygpath "$USERPROFILE")/.coreomics_token" && (umask 077; printf '%s\n' "$TOK" > "$F") && chmod 600 "$F"; unset TOK
     ```
     That is Git Bash on Windows: it saves where Python looks (`USERPROFILE`), which Git
     Bash's own `~` need not be on a domain account. There `chmod` may not restrict access
     (Git Bash does not map it onto NTFS permissions); the profile folder is private to you by
     default. macOS/Linux: `F="$HOME/.coreomics_token"`, and `chmod 600` makes it yours alone.
  3. `python3 scripts/core_submission.py check` → `ok`. Anything else: its `fix`, and the table
     in `references/core-submissions.md` "Troubleshooting `check`". `no_lab_access` means the
     key works but the CoreOmics account is not in the Proteomics lab: ask Brett (brettsp).

  Never paste the key into the chat. These run in Python on your computer, so they need a
  working local Python 3 — not the Windows Store stub (`local_python3.usable: false`; see
  "Python on this computer" there).
- **HIVE membership in `proteomics-grp`** — `locate`, `stage` and `deliver` read the Flinders
  raw data and write the service directory and Bioshare share (`/nfs/lssc0/flinders/proteomics`).
- **The whole skill at `~/proteomics-pipeline/`** on HIVE (`hive_exec.sh --put-skill`, step 1
  above) — `core_submission.py` imports its neighbours, and a partial copy exits 3 asking for a re-sync.
- For an **HT plate** instead: a STAN share token (`references/ht-submissions.md`).

`hive_exec.sh` reuses one SSH connection for 10 minutes, because HIVE throttles rapid new
connections (`HIVE_SSH_MUX=0` to disable; off automatically on native Windows). Its `--get`
copies with `rsync -rlt`: `-a` fails with exit 23 on macOS when pulling from HIVE's setgid
folders, even though the files arrive.
→ runbook: `references/core-submissions.md`.

## Rebuild on HIVE (non-Core users) — exact steps, no guessing
You have HIVE compute but not the Core's `/quobyte/proteomics-grp` software, so build
your own copy in your HIVE home. Run all of this **on HIVE** (via `hive_exec.sh`):
1. **Skill scripts + base env:**
   ```
   bash scripts/hive_exec.sh 'bash ~/proteomics-pipeline/scripts/setup.sh'
   ```
   This installs micromamba + R + limpa + limma + arrow + Sage + Python/pyarrow into
   `~/.proteomics-pipeline/` (no admin). DE, Sage (DDA), and figures now work.
2. **DIA-NN (for DIA data):** you can't read the group's DIA-NN builds, so let the skill
   fetch the free academic Linux build into your home. Pin the version the workflow
   manifest pins (`engine.version`; `resolve_defaults.py` pins 2.7.0) — any other version
   makes `run_search.py` print an engine version mismatch WARNING:
   ```
   bash scripts/hive_exec.sh 'PIN_ENGINE=diann PIN_VERSION=2.7.0 bash ~/proteomics-pipeline/scripts/acquire_tools.sh hpc'
   ```
   `acquire_tools.sh` downloads the DIA-NN Academia Linux zip to
   `~/.proteomics-pipeline/tools/diann/<version>/` (needs glibc ≥ Mint 21.2 / .NET 8;
   if the native binary won't run on the node, build the Apptainer image from the
   zip's Dockerfile). Check its `tools.json` `notes`.
3. **FASTA:** without the Core's pre-staged proteomes, download from UniProt.
   Resolve the user's organism first — never assume human:
   ```
   python3 scripts/fetch_fasta.py resolve --organism "<what the user said>"
   # no usable local python3 (check_access.sh local_python3.usable=false, e.g. Windows)?
   #   run the same resolve on HIVE:  bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/fetch_fasta.py resolve --organism "<...>"'
   bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/fetch_fasta.py fetch --proteome <confirmed UPID> --content one_per_gene --contaminants universal --out ~/proteomics-pipeline/search.fasta'
   ```
4. Then run the search via SLURM as above. Everything else (DE/figures/report/audit/
   reproducibility) is identical to a local run.

## Notes
- **Never run computationally intensive work on the login/head node** — the search
  and any heavy DE/figure/conversion step go through `sbatch`. The login node is only
  for submitting jobs, `squeue`/`sacct` polling, and small file moves. Holds on HIVE
  and any other cluster/scheduler.
- "Core member" = read access to `/quobyte/proteomics-grp`, reported by `check_access.sh`
  as `core_member` (readable locally — `local_proteomics_grp_access` — or over SSH —
  `hive_ssh_has_proteomics_grp`). `local_proteomics_grp_access: false` on a laptop is
  normal and says nothing about membership. If `hive_ssh` is `ok` and
  `hive_ssh_has_proteomics_grp` is false, the HIVE account isn't in the group yet —
  request it from the Core.
- Licensed software (e.g. Spectronaut) only runs where licensed; Sage + DIA-NN Academia
  run anywhere.
- Laptop→HIVE staging of large raw data is real bandwidth; prefer pointing at data
  already on HIVE/the proteomics share — `hive_path.sh` finds the HIVE path of a folder
  on a mapped drive or SMB mount.
