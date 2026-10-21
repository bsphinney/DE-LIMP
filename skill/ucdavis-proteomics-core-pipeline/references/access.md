# Access & where the skill runs (HIVE runbook)

**Model: Claude Code runs LOCALLY on the user's computer.** When a step should run on
HIVE, the local Claude Code drives HIVE **over SSH using the user's private key**
(`ssh -i <key> <user>@hive.hpc.ucdavis.edu '<command>'`, wrapped by `hive_exec.sh`).
The skill is never installed *on* HIVE for the user; it submits work there.

**Why local (not Claude Code on HIVE):** running Claude Code on a HIVE interactive
node ties up that node and the session **times out** mid-analysis. Keeping Claude Code
on the laptop means it stays alive while the actual compute runs as detached **SLURM
jobs** on HIVE — nothing depends on an interactive session staying open.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 32: Step 0a asks two questions (Core staff, three) → pick a mode
- 65: HPC@UCD's `ucd-hpc` skill (anyone on a UC Davis cluster)
- 90: hive_remote runbook (Claude Code local → HIVE over SSH)
- 159: Keeping the skill current (hive_remote)
- 277: Data already on a network drive (check BEFORE any upload)
- 310: Windows (Git Bash)
- 330: HIVE host key (first contact)
- 344: Core staff — submission workflow prerequisites
- 385: Rebuild on HIVE (non-Core users) — exact steps, no guessing
- 419: Notes
- 438: On HIVE: the login node or a compute node (golden rule 3, step 0a)
- 486: Step 0a: the questions, the check and the modes
- 558: Step 0: is this skill current? (every conversation)
- 680: Step 1: locating the raw files
- 685: Step 1: a path the user gave on a network drive
<!-- /contents -->

## Step 0a asks two questions (Core staff, three) → pick a mode
1. **"Do you have access to UC Davis HIVE (account + SSH private key)?"**
   If yes, **ask where the private key is** (e.g. `~/.ssh/id_ed25519`).
2. **"Do you work for the UC Davis Proteomics Core (staff)?"**
3. **Core staff only** (yes to 2, and `check_access.sh` says `core_member: true`), right after
   that check: **"Which CoreOmics submission is this analysis for?"** — a number, several, or
   "none" for internal work (SKILL.md step 0a; `references/board.md` section 2). Anyone else is
   never asked, and never hears of CoreOmics or the Core Project Board.

Verify before trusting the answers:
```
bash scripts/check_access.sh <hive_user> <private_key_path>
```
Reads `recommended_mode` + `core_member` (= `facility_software_available`). With
`core_member: false`, `coreomics_runs_on` is `none` and there is no `coreomics_key_advice`. If
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

## HPC@UCD's `ucd-hpc` skill (anyone on a UC Davis cluster)
HPC@UCD (HPCCF) made its own agent skill mandatory on Farm, Franklin and Hive on 2026-10-16:
"agents that do not follow it put the user's account at risk of suspension". On cluster policy
(what may run on a login node, `--account`/`--partition`/`--time`/memory, software installs,
home vs group storage, one writer per file on Quobyte) it supersedes this skill: where they
disagree, follow `ucd-hpc`, tell the user which instruction you set aside, and record it
(`report_issue.sh`).

Claude reads a skill only from the computer Claude Code runs on, so it belongs on the laptop in
`hive_remote` and on HIVE in `hive_local`. `check_access.sh` reports `ucd_hpc_skill`
(`installed`, `required` = this user works on a UC Davis cluster, `install`, `note`). The
install, from the ucd-hpc README (ask the user first; it goes in their own `~/.claude`):
```
git clone https://github.com/ucdavis/ucdhpc-cluster-skill ~/ucdhpc-cluster-skill
mkdir -p ~/.claude/skills
ln -s ~/ucdhpc-cluster-skill/ucd-hpc ~/.claude/skills/ucd-hpc
```
- **Windows (Git Bash):** the same lines, under the user's home. Where `ln -s` cannot make a
  link, copy instead: `cp -r ~/ucdhpc-cluster-skill/ucd-hpc ~/.claude/skills/ucd-hpc`, and copy
  again after each update.
- Then `/reload-plugins` (or restart Claude Code), and load it (Skill `ucd-hpc`) before the
  first command on the cluster.
- **Updates:** HPCCF updates it on GitHub. `git -C ~/ucdhpc-cluster-skill pull --ff-only` every
  few days.

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
     `setup.sh` (it's the same micromamba env). On the login node it runs itself on a compute
     node (one blocking `srun`, up to 2 h; HIVE's compute nodes have internet), so run the call
     in the background and wait:
     ```
     bash scripts/hive_exec.sh 'bash ~/proteomics-pipeline/scripts/setup.sh'
     ```
   - **Non-Core HIVE user:** same `setup.sh`, plus you must acquire DIA-NN/Sage
     yourself (no access to the group's builds). See "Rebuild on HIVE".
   **Where data and results go: group storage, never home.** HIVE's home is 20 GB per user
   (HPC@UCD's `ucd-hpc`: "conda environments and every package cache belong on group
   storage"; data and results too). Raw uploads, the search output and its per-run `.quant`
   files go under the user's group share — `/quobyte/<pi>grp/$USER/proteomics-pipeline/`
   (`id -Gn` lists the `…grp` groups; `qinfo quota /quobyte/<grp>` its space); a Core session
   is already on `/quobyte/proteomics-grp`. Below, `<work>` is that folder. Only a user with no
   group share (free tier only) uses `~/proteomics-pipeline`, and then watches `df -h ~`.
3. **Stage the raw data** — usually there is nothing to stage. See "Data already on a
   network drive" below: run `hive_path.sh` on the raw folder FIRST, and use the HIVE
   path it verifies in place. Upload only when the data is on the laptop's own disk:
   ```
   bash scripts/hive_path.sh '<the path the user gave>'     # exit 3 = local disk
   bash scripts/hive_exec.sh --put /path/to/raw '<work>/data/'
   ```
4. **Run the search as SLURM jobs** (never the login node). Run `run_search.py --sbatch
   job.sh` ON HIVE, then submit exactly what it prints under "submit with:":
   ```
   bash scripts/hive_exec.sh 'cd <work> && python3 ~/proteomics-pipeline/scripts/run_search.py ... --sbatch job.sh'
   bash scripts/hive_exec.sh 'bash <out>/submit.sh'     # DIA-NN from a FASTA, or the 5-step chain
   bash scripts/hive_exec.sh 'cd <work> && sbatch job.sh'   # only if job.sh was written
   bash scripts/watch_run.sh --all <out> --hive          # submit.sh wrote <out>/jobs.txt
   ```
   A DIA-NN search that predicts its library from the FASTA is two jobs (library, then
   search) chained by `<out>/submit.sh`, and writes no `job.sh`; a search routed to the
   5-step chain (>5 files) also writes none, and exits 3.
5. **DE + figures + report:** run on HIVE as a job (`run_de.R`, `make_figures.R`, …; the job's
   top from `on_compute.py --header`, SKILL.md step 8), never on the login node, or pull
   `report.parquet` back and run them locally.
6. **Retrieve results** into the session folder on the user's machine:
   ```
   bash scripts/hive_exec.sh --get '<work>/out' ./<session>/output/
   ```
   The search's per-run `.quant` files sit in `<out>/quant` (single-shot) or
   `<out>/quant_step2`/`quant_step4` (5-step chain) — kept there, not beside the raw data on
   the instrument share; in home they would count toward its 20 GB. They are only needed to
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
- **Updating the plugin**: `bash scripts/skill_version.sh --update`. Step 0 runs it once the user
  agrees, when `--check-update` (this computer vs GitHub main, every mode) or `--check-hive`
  (vs the Core's release) says this computer is behind. It refreshes the plugin catalogue (the
  marketplace clone) over HTTPS when it is behind main:
  `git -C ~/.claude/plugins/marketplaces/ucdavis-proteomics-core pull --ff-only`. A catalogue
  with an SSH remote is pulled from `https://github.com/bsphinney/DE-LIMP.git` instead. Then it
  runs `claude plugin update ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core` and
  compares the version that landed with main. Exit 1 means the update did not happen, and it
  never says "latest" while behind. Then `/reload-plugins`.
  - Why: on 2026-10-07 (msalemi, Windows, Git Bash, no GitHub SSH key) `claude plugin update`
    alone warned "marketplace not refreshed: SSH host key is not in your known_hosts", then
    "SSH authentication failed", and said "already at the latest version (2.10.0)" while main
    shipped 2.11.2. The catalogue's remote was HTTPS. The same `git pull` and update took it
    to 2.11.2.
  - Checked on Claude Code 2.1.285: `/plugin marketplace update ucdavis-proteomics-core` only
    refreshes the catalogue, and the installed plugin stays where it was.
  - Automatic updates are **off** for this marketplace until the user turns them on: `/plugin`
    → **Marketplaces** → **ucdavis-proteomics-core** → **Enable auto-update**. That default is
    the likely reason staff sat on 2.6.0.
  - A new install adds the marketplace by its HTTPS URL (`references/install.md`); Claude Code
    probes SSH only for the `bsphinney/DE-LIMP` shorthand. Claude Desktop: **+** → **Plugins**
    → **Manage plugins**.

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
- **Your own CoreOmics API key, where the CoreOmics steps run.** `check`, `identify`, `fetch`,
  `bioshare` and `email-draft` run on your computer when it has a working Python, and **on HIVE**
  (through `hive_exec.sh`, reading `~/.coreomics_token` there, mode 600) when it has none —
  `check_access.sh` says which (`coreomics_runs_on`) and what to do (`coreomics_key_advice`).
  On HIVE the agent sets the key up once: it sends a key saved on your computer with
  `bash scripts/coreomics_key_to_hive.sh` (over ssh standard input, never shown), or you paste
  one at that script's `--paste` prompt in your own Git Bash
  (`references/core-submissions.md` "The key on HIVE"). Saving it on your computer:
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

  Never paste the key into the chat. With no working local Python 3 (the Windows Store stub:
  `local_python3.usable: false`) they run on HIVE instead, with the key there — nothing to
  install (see "Python on this computer" there).
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
   `~/.proteomics-pipeline/` (no admin; about 3 GB of the 20 GB home — HPC@UCD's `ucd-hpc`
   would have environments on group storage, which the skill does not do yet). DE, Sage (DDA),
   and figures now work.
2. **DIA-NN (for DIA data):** you can't read the group's DIA-NN builds, so let the skill
   fetch the free academic Linux build into your home. Pin the version the workflow
   manifest pins (`engine.version`; `resolve_defaults.py` pins 2.7.0) — any other version
   makes `run_search.py` print an engine version mismatch WARNING:
   ```
   bash scripts/hive_exec.sh 'PIN_ENGINE=diann PIN_VERSION=2.7.0 bash ~/proteomics-pipeline/scripts/acquire_tools.sh hpc'
   ```
   `acquire_tools.sh` downloads the DIA-NN Academia Linux zip to
   `~/.proteomics-pipeline/tools/diann/<version>/` — on the login node the download runs itself
   on a compute node (one blocking `srun`), so run the call in the background — (needs glibc ≥ Mint 21.2 / .NET 8;
   if the native binary won't run on the node, build the Apptainer image from the
   zip's Dockerfile). Check its `tools.json` `notes`.
3. **FASTA:** without the Core's pre-staged proteomes, download from UniProt.
   Resolve the user's organism first — never assume human:
   ```
   python3 scripts/fetch_fasta.py resolve --organism "<what the user said>"
   # no usable local python3 (check_access.sh local_python3.usable=false, e.g. Windows)?
   #   run the same resolve on HIVE:  bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/fetch_fasta.py resolve --organism "<...>"'
   bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/fetch_fasta.py fetch --proteome <confirmed UPID> --content one_per_gene --contaminants universal --out <work>/search.fasta'   # <work>: group storage (runbook above); it runs itself on a compute node
   ```
4. Then run the search via SLURM as above. Everything else (DE/figures/report/audit/
   reproducibility) is identical to a local run.

## Notes
- **Never run computationally intensive work on the login/head node** — the search
  and any heavy DE/figure/conversion step go through `sbatch`. The login node is only
  for submitting jobs, `squeue`/`sacct` polling, small file moves and the bookkeeping
  scripts SKILL.md step 0a lists; anything else runs on a compute node through
  `hive_exec.sh --job` (one blocking `srun`) or a job. Holds on HIVE and any other
  cluster/scheduler; on HIVE, Farm and Franklin HPC@UCD's `ucd-hpc` skill is the rule.
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

## On HIVE: the login node or a compute node (golden rule 3, step 0a)

3. **Never run computationally intensive work on a cluster login/head node — EVER.**
   On any cluster (HIVE/SLURM, or any other scheduler), **every** heavy step — the
   search, and any large DE / figure / conversion (e.g. msconvert) — must run as a
   **scheduled job** (`run_search.py --sbatch` → `sbatch`; a DIA-NN search of >5 files on
   SLURM instead routes to the 5-step chain, exits 3 without writing `--sbatch`, and is
   submitted with `bash <out>/submit.sh`), not inline on the login node. Login nodes are shared; running compute there gets the user flagged/killed.
   The **only** things allowed on the login node are tiny orchestration commands —
   submitting jobs, `squeue`/`sacct` polling (`watch_run.sh`), small file moves and the
   bookkeeping scripts step 0a lists. Everything else runs on a compute node:
   `hive_exec.sh --job '<cmd>'` (one blocking `srun`) or an sbatch job (step 0a, "On HIVE").
   If you're unsure whether a step is heavy, submit it as a job. (No SLURM but a big
   dataset locally? Warn the user it may be slow rather than hammering a shared host.)
   **On a UC Davis cluster (HIVE, Farm, Franklin), HPC@UCD's own `ucd-hpc` skill is
   mandatory, and on cluster policy it wins over this skill**: what may run on a login node,
   resource requests (`--account`/`--partition`/`--time`/memory), software installs and
   storage. Load it (Skill `ucd-hpc`) before the first command on the cluster. Where the two
   disagree, follow `ucd-hpc`, tell the user which instruction of this skill you set aside, and
   record it (`report_issue.sh`). Step 0a checks that it is installed.


- **On HIVE: the login node or a compute node** (golden rule #3, HPC@UCD's `ucd-hpc` rules).
  A plain `hive_exec.sh '<cmd>'` runs on HIVE's shared login node, where only seconds of light
  work belong: `sbatch`/`squeue`/`sacct`/`scancel`, `watch_run.sh`, `ls`/`cat` of one folder or
  a small file, `mkdir`, `--put`/`--get`, writing a job (`run_search.py --sbatch`,
  `diann_parallel.py`, `make_podcast.py check`/`job`), and the bookkeeping scripts (`board.sh`,
  `notes.py`, `log_decision.py`, `report_issue.sh`, `skill_version.sh`, `session.py init`,
  `checkpoint.py status`, `submission_report.py attach`, `core_submission.py check`/`identify`/
  `fetch`/`bioshare`/`email-draft`, `node_fault.py retry`, `cpu_hours.py`). **Everything else
  goes to a compute node** — whatever reads `report.parquet` or raw files, runs R, walks a
  share, installs, renders or zips: `bash scripts/hive_exec.sh --job [--time MIN] [--mem SIZE]
  [--cpus N] [--preemptible] '<cmd>'` (one blocking `srun` on a queue that does not preempt,
  spelled out; defaults 30 min, 4G, 1 CPU; its output and exit status are the command's; when
  your own searches fill your 64-CPU cap on `high` it says it waits behind them, and
  `--preemptible` — only for a command that just reads and can be run again, as step 2 — runs it
  on `low` instead) for minutes of work, or an sbatch job for
  longer work and anything that must survive a closed laptop (step 8: `on_compute.py --header`).
  These move themselves when started on the login node: `setup.sh`, `ensure_dotnet8.sh`,
  `fetch_fasta.py fetch`/`check-db`, `core_submission.py locate`/`stage`/`deliver`,
  `provenance.py`, `session.py finalize`, `checkpoint.py find`, and `acquire_tools.sh` when it
  downloads. **Each of these blocks until its job ends (queue wait + the work): run the call
  with `run_in_background` and wait for it, or with a Bash timeout of at least its `--time` plus
  10 minutes — never the 2-minute default. A killed call kills the ssh, which cancels the srun
  half way (in the middle of an `--apply` copy or an install).** Unsure → `--job`. In
  `hive_local` (Claude Code on HIVE) the same split, with `python3 scripts/on_compute.py [--time
  MIN] [--mem SIZE] -- <cmd>` for `--job`.

## Step 0a: the questions, the check and the modes

**Claude Code runs locally on the user's machine.** Heavy steps either run locally or
are driven on **HIVE over SSH with the user's private key** — they do NOT run Claude
Code on HIVE. Ask two questions, verify, and ask Core staff a third:

3. **Core staff only — yes to 2 AND `check_access.sh` below says `core_member: true` — right
   after that check: "Which CoreOmics submission is this analysis for?"** A number (807 or
   PROT_0807), several numbers when the analysis combines submissions, or "none" for internal
   work (a method test, QC, training, a tool comparison). Skip it when they already said (the
   request names a submission, or says it is a method test) and when the user is coming back to
   an earlier analysis (step 0c checks it). Every analysis by Core staff goes on the Core
   Project Board before its long steps: a number → step 1c, internal work → step 1. Pass the
   answer to `session.py init` at step 3b: `--submission PROT_0807` (once per number, the main
   one first) or `--internal '<its title on the board>'`. → `references/board.md` section 2.
   **Anyone else (no to 2, or `core_member: false`): never ask it, skip every board step, and
   never mention the board, CoreOmics or a submission number** (`check_access.sh` then says
   `coreomics_runs_on: none`).

Verify (don't just trust the answers):

```
bash scripts/check_access.sh <hive_user> <private_key_path>
```

Read `recommended_mode` + `core_member`, then pick the mode below. **`ucd_hpc_skill`**: when
`required` is true and `installed` is false, tell the user in one line that HPC@UCD requires its
`ucd-hpc` skill for any agent on its clusters, ask, then run its `install` line here (it goes in
their own `~/.claude/skills`; relay `note` for Windows) and ask them to type `/reload-plugins`.
Until it loads, follow golden rule #3 to the letter. Once the SSH login works,
`check_access.sh` saves the user and the key's path (never the key) to
`~/.config/ucdavis-proteomics/hive.env` (`hive_login_saved`). Later calls read the login there,
because the environment does not carry over from one command to the next.

- **HIVE = yes → `hive_remote`:** drive HIVE over SSH from the local Claude Code
  (`export HIVE_USER=… HIVE_KEY=…`; use `bash scripts/hive_exec.sh '<cmd>'`). **Every
  script runs ON HIVE through `hive_exec.sh` — the laptop needs no Python, R or engines.**
  **Free text in a command** (a title, a `--text`, `--reason` or `--why`) is read by the shell
  here and, through `hive_exec.sh`, again on HIVE, so a backtick or `$(...)` in it would run:
  on this computer put it in single quotes with no `'` in it; inside a `hive_exec.sh` command
  use none of `$` `` ` `` `\` `"` `'` in it. Longer text goes in a file, on standard input.
  On Windows, Git Bash + the key is enough: `python3` there is often the Microsoft Store
  stub (`check_access.sh` → `local_python3.usable: false` → never run a skill script
  locally; the CoreOmics steps then run on HIVE too, `coreomics_runs_on: hive`, step 1c.1), rsync is absent (`--put/--get` fall back to scp by themselves), and
  `HIVE_USER` is the plain UC Davis id, never the Windows `DOMAIN+user` login. If
  `hive_ssh` is `failed`, relay `hive_ssh_error.kind` (host_key / permission_denied /
  timeout) as the `kind` table at the top of this file describes. The search runs as **SLURM jobs**, never
  the login node: `run_search.py --sbatch job.sh`, then submit **what it prints under
  "submit with:"** — `hive_exec.sh 'bash <out>/submit.sh'` for a DIA-NN search that
  predicts its library from the FASTA (two jobs: library, then search — the usual case) or
  that routed to the 5-step chain (>5 files; exit 3), `hive_exec.sh 'sbatch job.sh'` only
  when a single `job.sh` was written. **The skill must be there, at this computer's
  version:** `bash scripts/hive_exec.sh --put-skill` (`scripts/` and `.claude-plugin/` into
  `~/proteomics-pipeline/`; the latter carries the version every record states). Step 0's
  `skill_version.sh --check-hive` says when that is needed. HIVE gives **compute**; the Core
  software is separate (next).


- **Core member = yes (with HIVE) → reuse the installed software** in
  `/quobyte/proteomics-grp/`: `acquire_tools.sh` finds the Core's DIA-NN builds,
  `fetch_fasta.py --hive` reuses pre-staged FASTAs. No rebuilding.
- **HIVE = yes but NOT Core → rebuild the toolchain in your own HIVE home** (you can't
  read the group dir). Follow "Rebuild on HIVE" (this file) (run
  `setup.sh` + `acquire_tools.sh` + `fetch_fasta.py` on HIVE via `hive_exec.sh`). **List
  these steps for the user — do not guess.**
- **Core = yes but NO HIVE →** the Core software is on HIVE, unreachable without HIVE
  access → fall back to **local**; suggest requesting a HIVE account.
- **Both = no → `local`:** `setup.sh` installs the toolchain on the user's machine
  (no admin); public engines (DIA-NN Academia, Sage). Self-contained, works fully.

Tell the user which mode you're using and why.

## Step 0: is this skill current? (every conversation)

**One call does both checks** (SKILL.md step 0 runs it at the start of every conversation):
```
bash scripts/skill_version.sh --session-start [--mode hive_remote]
```
It runs `--check-update` and then, only with `--mode hive_remote`, `--check-hive` (not while a
reload is pending), and prints one JSON object. Act on its fields in this order:

| Field | What to do |
|---|---|
| `up_to_date` | `true` only when `--check-update` exited 0 and HIVE's release check is not behind: the one case in which the skill may be called up to date. |
| `relay_first` | Relay it before any other work (older than `main` or the Core's release; or the `/reload-plugins` line of `reload_needed`). |
| `ask_update` | Ask whether to update now. Yes → `bash scripts/skill_version.sh --update`, then `/reload-plugins`, then start again. No → record that they declined, once the session exists (below). |
| `ask` | A question for the user: HIVE has a newer build the update would not bring here (never put the older copy over it unasked), or HIVE's copy is in use by the user's jobs (put it now or after they finish). |
| `next` | Run it at once, without asking (`hive_exec.sh --put-skill`: the user's own copy), after any question, then `--session-start` again. |
| `say` | One line to relay: a check that could not run, HIVE's other notes. |
| `update`, `hive` | The two checks' own objects, each with its `exit` (`hive` is null when it did not run). |

Exit 2: the user must answer or act first; 3: `next` is set; 0: carry on. Every exit code of the
two checks, and why each exists, follows; `--check-update` alone is what step 0's "Check again"
runs before steps 7 and 12.

**The version check comes first, in every conversation: also when this computer is already set up
(the setup, install.md "Step 0: the one-time setup", then only verifies), also on a resume (0c), and again before step 7 and step 12
(SKILL.md step 0, "Check again").** Staff reports of 2026-10-16 still came from 2.11.2, the day after 2.12.0 and
2.13.0 were released: conversations begun before a release ran on into the next day without
checking again, and the check sat in the middle of a step titled "One-time setup", easy to skip
on a computer already set up.

**Every mode, first: is this computer's skill behind GitHub `main`?** One HTTPS read of
`main`'s `plugin.json` (bash and curl, which Git Bash has; no Python, no GitHub SSH key):

```
bash scripts/skill_version.sh --check-update
```

- `0` → this copy is what `main` ships, or newer (a test build). Say nothing.
- `3` → `behind_main` true: the update question below, with `main` as the newer version.
  `status: reload_needed` means the new version is installed already and only this session runs
  the old one: ask the user to type `/reload-plugins` (or restart Claude Code), then start again.
- `5` → `main` could not be read (offline, no curl) → relay `say` in one line and carry on.

**Never tell the user the skill is up to date** unless `--check-update` exited 0, and never take
it from `claude plugin update`. On 2026-10-07 a Windows laptop with no GitHub SSH key ran
`claude plugin update ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core`: it warned
"marketplace not refreshed: SSH host key is not in your known_hosts", then "SSH authentication
failed", and printed "already at the latest version (2.10.0)" while `main` shipped 2.11.2. It
had compared against its own stale catalogue (the marketplace clone).

**In `hive_remote`, then check that HIVE runs this same skill** (one SSH call). Core staff once
ran 2.6.0 on their laptops and on HIVE two releases after it, so none of the fixes in between
ran, and a HIVE copy without `.claude-plugin/` wrote "v0.0.0" into sdrf.tsv:

```
bash scripts/skill_version.sh --check-hive --mode hive_remote
```

**`behind_main` true (above) or `behind_release` true** means this computer's skill is older
than the Core's current release. `behind_release` can come with exit 0, 3, 4 or 6. **Relay it
FIRST, before any other work**, in plain words: "Your copy of the skill is older than the Core's
current release (you have `local`, the Core is on `main` / `release`). Updating takes a minute."
Then **ask whether to update now**:

- **Yes** → run the update yourself; the user types nothing:
  ```
  bash scripts/skill_version.sh --update
  ```

  When the catalogue is behind `main` it refreshes it over HTTPS:
  `git -C ~/.claude/plugins/marketplaces/ucdavis-proteomics-core pull --ff-only`, the route that
  worked on that laptop. A catalogue with an SSH remote is pulled from
  `https://github.com/bsphinney/DE-LIMP.git` instead. Then it runs
  `claude plugin update ucdavis-proteomics-core-pipeline@ucdavis-proteomics-core` with
  `CLAUDE_CODE_PLUGIN_PREFER_HTTPS=1` (Claude Code's documented switch for machines without a
  GitHub SSH key). Last, it reads the version that landed and compares it with `main`:
  - exit 0, `updated` → ask the user to type `/reload-plugins` (or restart Claude Code), then
    start again with `--check-update`. `current` or `reload_needed`: nothing, or only the
    reload, was needed.
  - exit 1, `not_updated` → **The update did NOT happen**. Say it in those words and relay `say`
    (why, plus the by-hand line in `manual`), then record it (`report_issue.sh`). Never relay
    claude's "already at the latest version" as the answer; `claude_said` keeps its words for
    the record.
  - exit 5, `unverified` → it ran, but the result could not be checked against `main`. Relay
    `say`: it is not a confirmation.
- **No** → carry on. Once the session exists (step 3b), record that they declined:

```
python3 scripts/log_decision.py --session <S> --step '0. skill version' --what-file - --why-file - <<'EOF'
Skill update declined: ran <local>, the Core's release is <release>
---
<the user's words>
EOF
```

Either way, **always recommend turning on auto-update**: `/plugin` → Marketplaces →
ucdavis-proteomics-core → Enable auto-update. It is off by default for this marketplace, which
is how staff fell two releases behind. A new install should add the marketplace by its HTTPS
URL, which Claude Code clones over HTTPS without probing SSH (`references/install.md`).

Then, by exit code:

- **exit 0** → carry on (after the update question above, when `behind_release` is true).
- **exit 3** → `~/proteomics-pipeline` on HIVE is missing, has no `.claude-plugin/`, is older, or
  is another build (its files differ) → run `bash scripts/hive_exec.sh --put-skill` without
  asking (it is the user's own copy), then the check again, which must exit 0; if it does not,
  record it (`report_issue.sh`) and tell the user before going on. Relay `say` in one line.
- **exit 4** → HIVE has a NEWER skill than this computer → stop and relay `say`. With
  `behind_release` true, that is the update question above. Otherwise HIVE holds an unreleased
  build: ask the user which to keep. `--put-skill` the older one over it only if the user says
  so.
- **exit 5** → the check could not run (no saved HIVE login, HIVE unreachable) → say so in one
  line and carry on. With no login, `check_access.sh` (step 0a) saves one; then check again.
- **exit 6** → as exit 3, but `jobs` of the user's SLURM jobs have batch scripts that use
  `~/proteomics-pipeline`, and their end-of-job steps run its scripts hours later. Or squeue
  answered with an error, so they could not be counted. → Do NOT put it now. Relay `say`, and
  ask the user whether to put it now or after those jobs finish (`watch_run.sh` shows when);
  until then HIVE runs the older scripts.

Only `hive_remote` runs `--check-hive`. Without `--mode hive_remote`, a computer with no HIVE
login skips it silently; the release is compared only for Core accounts.

## Step 1: locating the raw files

Ask the user for a directory or file list if not given. Recognized: `.d` (Bruker),
`.raw` (Thermo), `.mzML[.gz]`, `.wiff` (convert first). Glob to a concrete list.

## Step 1: a path the user gave on a network drive

**`hive_remote`: before any upload, run `bash scripts/hive_path.sh '<path the user gave>'`.**
A mapped drive (`T:\…`) or SMB mount (`/Volumes/proteomics/…`) is usually a Flinders or
quobyte share that HIVE already mounts. `verified: true` → use `hive_path` in place and
glob on HIVE (`hive_exec.sh 'ls <hive_path>'`) — upload nothing. `verified: false` → show
`candidates` and **ask the user** where that share lives on HIVE; do not search HIVE for
it (a `find` over the Flinders NFS mount does not finish). Exit 3 (the laptop's own disk)
→ `hive_exec.sh --put`, which itself refuses a source it verifies on HIVE.

- **Any** drive letter (`R:`, `T:`, …), UNC path (`\\server\share\…`) or SMB mount gets this
  check first — never search `/quobyte` (or anywhere on HIVE) instead: not finding the files
  there does not mean HIVE lacks them. `R:` was `/nfs/lssc0/flinders/proteomics`, and 70 GB
  already on HIVE was uploaded (msalemi, 2026-09-24).
