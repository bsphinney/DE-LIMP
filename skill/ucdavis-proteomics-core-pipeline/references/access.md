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
Reads `recommended_mode` + `facility_software_available`. Then:

| HIVE access | Core member | Mode | What happens |
|---|---|---|---|
| no | no | **local** | `setup.sh` installs the toolchain on the user's machine; public engines (DIA-NN Academia, Sage). |
| **yes** | **yes** | **hive_remote** | Drive HIVE over SSH. **Reuse the Core software already installed** in `/quobyte/proteomics-grp` (DIA-NN `.sif`, pre-staged FASTAs). Search runs as a SLURM job. |
| **yes** | no | **hive_remote** | Drive HIVE over SSH, but you **rebuild the toolchain in your own HIVE home** (you can't read the Core group dir). See "Rebuild on HIVE" below. |
| no | yes | **local** | The Core software is on HIVE; without HIVE access you can't reach it → run locally with public engines and tell the user to request a HIVE account. |

## hive_remote runbook (Claude Code local → HIVE over SSH)
Set the connection once, then everything HIVE-side goes through `hive_exec.sh`:
```
export HIVE_USER=<hive_user>
export HIVE_KEY=<private_key_path>        # the path the user gave you
bash scripts/hive_exec.sh 'hostname; sbatch --version | head -1'   # confirm
```
1. **Put the skill's scripts on HIVE once** (they run there):
   ```
   bash scripts/hive_exec.sh 'mkdir -p ~/proteomics-pipeline'
   bash scripts/hive_exec.sh --put ./scripts '~/proteomics-pipeline/'
   ```
2. **Toolchain on HIVE:**
   - **Core member:** `acquire_tools.sh` (run on HIVE) finds the group DIA-NN
     container (`/quobyte/proteomics-grp/dia-nn/*.sif`); `fetch_fasta.py --hive`
     reuses `/quobyte/proteomics-grp/MRS/` FASTAs. Build the R/Python/DE env once with
     `setup.sh` (it's the same micromamba env):
     ```
     bash scripts/hive_exec.sh 'bash ~/proteomics-pipeline/scripts/setup.sh'
     ```
   - **Non-Core HIVE user:** same `setup.sh`, plus you must acquire DIA-NN/Sage
     yourself (no group `.sif`). See "Rebuild on HIVE".
3. **Stage the raw data** (skip if it's already on HIVE — Core data usually is):
   ```
   bash scripts/hive_exec.sh --put /path/to/raw '~/proteomics-pipeline/data/'
   ```
4. **Run the search as a SLURM job** (never the login node). Generate the sbatch
   script with `run_search.py --sbatch`, put it on HIVE, submit, and poll:
   ```
   # build job.sh locally (or on HIVE), then:
   bash scripts/hive_exec.sh --put ./job.sh '~/proteomics-pipeline/job.sh'
   bash scripts/hive_exec.sh 'cd ~/proteomics-pipeline && sbatch job.sh'
   bash scripts/hive_exec.sh 'squeue -u $USER'     # poll until done
   ```
5. **DE + figures + report:** run on HIVE (`run_de.R`, `make_figures.R`, …) or pull
   `report.parquet` back and run them locally (DE/figures are light).
6. **Retrieve results** into the session folder on the user's machine:
   ```
   bash scripts/hive_exec.sh --get '~/proteomics-pipeline/out' ./<session>/output/
   ```

## Rebuild on HIVE (non-Core users) — exact steps, no guessing
You have HIVE compute but not the Core's `/quobyte/proteomics-grp` software, so build
your own copy in your HIVE home. Run all of this **on HIVE** (via `hive_exec.sh`):
1. **Skill scripts + base env:**
   ```
   bash scripts/hive_exec.sh 'bash ~/proteomics-pipeline/scripts/setup.sh'
   ```
   This installs micromamba + R + limpa + limma + arrow + Sage + Python/pyarrow into
   `~/.proteomics-pipeline/` (no admin). DE, Sage (DDA), and figures now work.
2. **DIA-NN (for DIA data):** you can't use the group `.sif`, so let the skill fetch
   the free academic Linux build into your home:
   ```
   bash scripts/hive_exec.sh 'PIN_ENGINE=diann PIN_VERSION=2.6.0 bash ~/proteomics-pipeline/scripts/acquire_tools.sh hpc'
   ```
   `acquire_tools.sh` downloads the DIA-NN Academia Linux zip to
   `~/.proteomics-pipeline/tools/diann/<version>/` (needs glibc ≥ Mint 21.2 / .NET 8;
   if the native binary won't run on the node, build the Apptainer image from the
   zip's Dockerfile). Check its `tools.json` `notes`.
3. **FASTA:** without the Core's pre-staged proteomes, download from UniProt:
   ```
   bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/fetch_fasta.py --proteome UP000005640 --add-contaminants --out ~/proteomics-pipeline/search.fasta'
   ```
4. Then run the search via SLURM as above. Everything else (DE/figures/report/audit/
   reproducibility) is identical to a local run.

## Notes
- **Never run computationally intensive work on the login/head node** — the search
  and any heavy DE/figure/conversion step go through `sbatch`. The login node is only
  for submitting jobs, `squeue`/`sacct` polling, and small file moves. Holds on HIVE
  and any other cluster/scheduler.
- "Core member" = read access to `/quobyte/proteomics-grp`. If `check_access.sh` shows
  `proteomics_grp_access: false` over SSH, the account isn't in the group yet — request
  it from the Core.
- Licensed software (e.g. Spectronaut) only runs where licensed; Sage + DIA-NN Academia
  run anywhere.
- Laptop→HIVE staging of large raw data is real bandwidth; prefer pointing at data
  already on HIVE/the proteomics share when possible.
