# Peak check (Skyline) and the targeted follow-up list

`scripts/skyline_peaks.py` (SKILL.md step 8f). Skyline's command line, `SkylineCmd`, runs headless
under Wine in ProteoWizard's vendor-reader container, as SLURM jobs. For the top proteins of a
finished DIA-NN analysis it extracts each peptide's chromatograms from the raw files again, draws
them, and lists the peaks where Skyline and DIA-NN disagree.

## What it is, and what it is not (the wording is fixed)

- **The pictures are a visual check of DIA-NN's peaks.** Skyline reads DIA-NN's spectral library
  for the search, and its library setting "use explicit peak bounds" is on by default (SkylineCmd
  has no switch for it). So Skyline integrates inside DIA-NN's peak edges, and **a peak that
  matches is not independent confirmation** of it.
- **The real warning signs are the peaks where Skyline picked a different peak**:
  `disagreements.tsv` lists every precursor-run where
  - Skyline's apex is more than 6 s from DIA-NN's (`APEX_DIFF_S`),
  - or Skyline's apex lies outside DIA-NN's peak edges (`RT.Start`-`RT.Stop`),
  - or Skyline found no peak where DIA-NN reported one.
- **Never write** "confirmed by Skyline", "validated", "verified", "corroborated" or "independent
  evidence" about it. `client_text.peak_check_overclaims()` finds those (a negation in the same
  clause clears a sentence: "not independent confirmation" is fine). `make_analysis_html.py` warns
  on them in `AI_Analysis_Report.md`. The report section (`report_section()`) is tested to have
  none.
- What Skyline does measure on its own: the apex, the areas, the library dot product (dotp) and
  the isotope dot product (idotp).

## When it runs

- **By default after DE in a Core analysis**, and for anyone else whose HIVE (or other SLURM
  cluster) can reach the pinned container. Run `plan` right after step 8's final DE; the jobs take
  minutes, so they are done long before the report (step 9).
- **The user can turn it off:** `skyline_peaks.py skip --session <S> --reason '<their words>'`.
  The report then says it was turned off and why. Jobs still running are cancelled, and an earlier
  run's pictures and document are removed, so none is delivered as current.
- **When it cannot run**, `plan` exits 3 and writes the record with the reason. The report's
  "Peak check (Skyline)" section starts "**Not run.**" with that reason, and so does finalize's
  `MANIFEST.txt` line. Reasons include:
  - no container here;
  - not a DIA-NN 2.x search, or no Skyline library beside the report;
  - DDA data (DIA-NN `--dda`): the document extracts DIA;
  - DIA-NN's mass accuracy is not recorded (see below);
  - the search's command does not state a setting the document must copy;
  - no carbamidomethyl cysteine (`--unimod4`) in the search;
  - not a Bruker timsTOF `.d` or Thermo `.raw` search;
  - ion-trap MS2;
  - no SLURM;
  - DIA-NN reports no precursor of the chosen proteins at the q-value filter, or none of
    Skyline's precursors matched DIA-NN's by peptide, modifications, charge and run name, so
    nothing could be compared (never an "ok" with 0 of 0).

**Matching the two programs' precursors** (`precursor_key`): the stripped sequence, the charge,
and each modified residue with its mass shift, summed per residue and rounded to a whole dalton.
Skyline writes `PEPM[+16]K` and an N-terminal acetyl with an oxidised M as one shift (`M[+58]`);
DIA-NN writes `PEPM(UniMod:35)K` and `(UniMod:1)M(UniMod:35)`. Their masses come from the
search's own `--var-mod` / `--fixed-mod` flags (`UniMod:<id>,<mass>,<sites>`), then from the skill's
UniMod tables (`make_methods.UNIMOD`, `silac.LABELS`). A residue with one of two modifications
never matches the residue with both, and a modification neither names is not compared rather than
compared wrongly. The report counts the Skyline precursor-runs that had no DIA-NN identification
to compare with, and the labels in the tables and pictures carry the modifications
(`PEPM[+16]K/2+`).
- **A session with a DE but no record** gets `[SKIPPED] Peak check (Skyline, output/skyline) --
  not run` in `MANIFEST.txt`: run `plan`, or `skip` with the user's reason.

```
python3 scripts/skyline_peaks.py plan --session <S> [--top N]     # writes the jobs; prints "submit with"
bash <S>/output/skyline/skyline_work/submit.sh                     # through hive_exec.sh in hive_remote
python3 scripts/skyline_peaks.py status --session <S>             # the record, summarised
python3 scripts/skyline_peaks.py skip --session <S> --reason '<the user's words>'
```

`plan` is light (it reads small records and writes the jobs): it belongs with the job-writing
commands on a login node. Everything heavy runs in the three jobs.

### Options (`plan`)

- `--top N`: the number of proteins (default 10, tagged `(DEFAULT — not user-confirmed)` in the
  report and the Methods when it is the default).
- `--name`: the document's name (default: the session folder's name).
- `--python`: the Python the jobs run. Default: setup.json's, the pipeline env, which has pyarrow,
  and matplotlib from 2.14 (`bash setup.sh` adds it to an older env). The build job checks both
  packages and writes the reason when one is missing.
- `--sif`: another copy of the container (see "The container").
- `--ppm-ms1 / --ppm-ms2 / --ppm-source`: only when the search records no mass accuracy, from the
  search's own record, with where it came from.
- `--max-parallel K`: at most K import tasks at once (default 8; 8 x 4 CPUs = 32, half of a user's
  64-CPU cap on HIVE's `high`).
- `--partition / --account / --qos`: the queue (default: the skill's queue rule,
  `run_search.slurm_queue`).
- `--replan`: plan again although the last plan has not written its result. Its jobs still
  running are cancelled first: only the ids `skyline_work/jobs.txt` names, and only those SLURM
  (`sacct`) says are yours, named as this step's jobs (`sky_*`) and not finished. When that cannot
  be checked, or a recorded job is not this step's, nothing is cancelled and `plan` refuses.
  Without `--replan`, `plan` refuses. `skip` cancels the same way. Each plan has an
  id, and a job writes the record only while it is its own plan's, so a job of an earlier plan
  never overwrites a newer record or a `skip`.
- **A DE run again** after the check makes the record **out of date**: the record keeps the
  sha256 of the `de_provenance.json` its proteins were chosen from. The report then says "Out of
  date", and `MANIFEST.txt` says to plan again. A record that cannot be read is said too ("Could
  not be read").
- **Every plan starts clean.** It removes the last run's pictures, tables, `.sky.zip` and its whole
  `skyline_work/`, so a stale per-file cache or status file is never taken for this run's. Until
  the new plan is written, the record says the planning did not finish; a plan stopped by an
  error records the error.
- **Jobs belong to their plan.** Each job is written with its plan's id (`--plan-id`) and does
  nothing when `plan.json` is another plan's; the record and the status files carry the id too.
  `submit.sh` writes each job to `jobs.txt` the moment it is submitted, so a submission refused
  half way still names the jobs it made, and it refuses to submit a plan a second time.
- Without `--session`: `--search <DIA-NN output folder> --out <folder>`, plus `--de-dir`,
  `--manifest` and `--conditions` as available.

## Which proteins: the ONE rule (`select_targets`)

1. **With DE hits:** the up-to-N significant proteins (adjusted p below the DE's own cutoff, from
   `de_provenance.json`), ranked by their smallest adjusted p in any comparison, ties broken by
   the raw p-value. This is the same order as the report's top-protein tables
   (`make_analysis_html.Tables`). Contaminant entries (`Cont_`, `contam_`) are left out.
2. **No significant protein:** the N most abundant proteins measured in every run. Abundance is the
   mean of the DE's log2 `Expression_Matrix.csv`, and "measured" is `Detection_Matrix.csv` > 0. The
   report says these were not chosen for a change.
3. **No DE at all** (a search-only session): the N most abundant protein groups of one accession,
   measured in every run with at least 3 precursors (DIA-NN's `PG.MaxLFQ`).

DIA-NN's rows for them pass the DE's q-value filter (`diann_q_columns`). A chosen protein that is
not in the searched FASTA by accession is dropped and listed in the record (`not_in_fasta`).

## Mass tolerance: DIA-NN's own, never a fixed number

In the feasibility trial, a guessed 10 ppm gave gappy, broken peaks on a Lumos run that DIA-NN had
measured at 24 ppm; 24-25 ppm fixed it. `mass_accuracy()` takes, in this order:

1. the value the search **fixed**: report.log.txt's "Mass accuracy will be fixed to ... (MS2) and
   ... (MS1)", or `--mass-acc` / `--mass-acc-ms1` in its command;
2. DIA-NN's **measured** "Optimised mass accuracy: N ppm" (MS2). DIA-NN logs only a
   *recommendation* for MS1, so MS1 uses the same value, as the trial did;
3. `search_provenance.json` `mass_acc`;
4. otherwise the step is **not run**, with that reason (or the user gives `--ppm-ms1 / --ppm-ms2`
   from the search's record, with `--ppm-source`).

Both programs read the value as +- ppm around the m/z: DIA-NN's `Scan::level` (src/diann.cpp,
tag 1.8) and Skyline's centroided window (`TransitionFullScan.GetFilterWindow`).

- **Thermo** (Orbitrap, Astral): centroided analysers at that ppm.
- **timsTOF**: TOF analysers. Skyline's TOF window is +- m/z / R, so +-P ppm is resolving power
  R = 1e6 / P (15 ppm -> 66,667). The library's ion mobility (1/K0) is used at
  `--ims-library-res=50`.
- The analyser is **always set with its resolution**. Setting the acquisition method alone left
  Skyline's default product analyser, "qit" at 0.7 (the trial).

## What the document copies from the search

From the DIA-NN command in `report.log.txt`, never guessed (a missing flag means the step is not
run):

- the cleavage rule, `--cut`, mapped to Skyline's predefined enzyme that cuts the same way
  (`K*,R*` is Trypsin/P, `K*,R*,!*P` is Trypsin; table `SKYLINE_ENZYMES`, from Skyline's
  `Settings.cs`);
- missed cleavages, peptide length, precursor charges, the fragment m/z range;
- carbamidomethyl C (Skyline's default; the step is not run without `--unimod4`);
- variable modifications, which Skyline adds from the library (`--import-search-add-mods`).

**Not carried over** (the report says so):
- DIA-NN's N-terminal methionine excision (`--met-excision`): no SkylineCmd option.
- DIA-NN's precursor m/z range (`--min-pr-mz`/`--max-pr-mz`): SkylineCmd has no precursor m/z
  filter of its own. The checked precursors are DIA-NN's library precursors, all inside it.

**Raw files** are the ones the search read (`search_provenance.json` `files`). A run whose path
holds `--` was searched through a link (`input_links`), which may rename it; Skyline imports the
same link, so its run name matches DIA-NN's.

**Fixed settings**, the ones the trial's agreement numbers were measured with:
- DIA, with the isolation scheme read from the first raw file;
- 3 precursor isotopes;
- product ions y, b and p at charge 1-2, the 6 best library ions (at least 3), matched at
  0.05 m/z;
- chromatograms within +-1 min of DIA-NN's identification times, or +-2 min when the runs go
  past 20 min (`rt_window`).

## The container

Pinned in `resolve_defaults.SKYLINE_PIN`:

- **Image:** Skyline **26.1.0.266** (578b175), the release. It was skyline.ms's current download on
  2026-10-20 (GitHub tag `Skyline-26.1.0.266`).
- **Docker:** `proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses@sha256:448a833eb92eac2ca731f2461c046963c0c2afafb3550b4dfb6fc3914cec7b0e`
  (tag `skyline_26.1.0.266-578b175`).
- **On HIVE:** `/quobyte/proteomics-grp/apptainers/pwiz-skyline-i-agree-to-the-vendor-licenses_skyline_26.1.0.266-578b175.sif`,
  sha256 `b5f135587e9b589dd2d4686ff520573ea237b84f53b1343d5b76213c16f302bf`.
- **Checked before every use:** the build job hashes the Core's `.sif` and refuses one that does
  not match. It also asks `SkylineCmd --version`, and refuses any other Skyline version. The record
  keeps the sha256, the version line and how it was checked.
- **Anyone else** (no access to the Core's folder): pull the same digest once, in a job, never on a
  login node. The pull took 7 min and 16.8 GB, so ask for 32 GB:
  ```
  srun --account=<yours> --partition=<yours> --time=60 --mem=32G --cpus-per-task=4 \
    apptainer pull <your folder>/skyline_26.1.0.266.sif \
    docker://proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses@sha256:448a833eb92eac2ca731f2461c046963c0c2afafb3550b4dfb6fc3914cec7b0e
  ```
  Then `plan ... --sif <that file>`. A `.sif` built elsewhere has another sha256, so it is checked
  by its `SkylineCmd --version` instead, and the record says so.
- **Licence:** the image includes the instrument vendors' reader DLLs. Using it means agreeing to
  their licences (proteowizard.sourceforge.io/licenses.html). This is the same image family the
  Core uses for `wine msconvert`. Skyline and ProteoWizard themselves are Apache 2.0.
- **A laptop** cannot run it: there is no SLURM, and the import runs one raw file per job. The
  record says so.

## The jobs

The jobs are written into `<S>/output/skyline/skyline_work/`. Each is a thin bash wrapper:
modules, `module load apptainer`, then one subcommand of `skyline_peaks.py`. Every job has an
explicit account, partition, QOS, time and memory, and one log per job (`%j`, `%A_%a`). They pass
ucd-hpc's `lint-jobscript.sh`.

| Job | What | Size |
|---|---|---|
| `sky_1_build.sbatch` | checks the container, chooses the proteins, cuts the FASTA to them, reads DIA-NN's rows for them, then **one** SkylineCmd call builds the document | 4 CPUs, 16G, 60 min |
| `sky_2_import.sbatch` | a SLURM array, **one raw file per task** (`--import-no-join`) | 4 CPUs, 16G, 60 min each, at most `--max-parallel` at once |
| `sky_3_join.sbatch` | **one** `--batch-commands` call joins the per-file caches and exports the precursor and transition reports, the chromatograms of the targets and the `.sky.zip` (`--share-type=minimal`); then it draws the pictures and writes the record | 4 CPUs, 16G, 90 min |

`submit.sh` chains the jobs:
- the imports wait for the build (`afterok`, `--kill-on-invalid-dep=yes`: a failed build cancels
  them);
- the join runs after the imports whatever happened (`afterany`), and writes the record of what
  did.

A run that one import task could not read is left out of the pictures and tables, and named in the
report.

**Rules from the trial, encoded:**
- **Never two readers on one timsTOF `.d`.** It holds a SQLite database, and a second reader gets
  "database is locked". So there is one task per file, and nothing else reads the raw data while
  the imports run.
- **Each import task works in its own copy of the document** (`skyline_work/tasks/<i>/`).
  Skyline keeps an import's spill files in `<document folder>/xic/` and deletes that folder when
  the import ends. When the tasks shared the document's folder, the first to finish deleted the
  others' spill files mid-import ("Unable to find the specified file"; the skill's first test on
  HIVE, 2026-10-20). The per-file cache is named after the raw file and its folder only, so the
  task moves it beside the document, where the join finds it.
- **Never `--import-process-count`**: it hung under Wine. `run_skyline()` refuses it.
- **Each SkylineCmd call costs 8-10 s of Wine start-up**, so the build is one call and the join's
  five steps are one `--batch-commands` file.
- **Paths go in as Windows paths**, `Z:\...` (Wine maps `Z:` to `/`).
- **Binds:**
  - the output folder is bound read-write;
  - the search folder, the FASTA's folder and every raw-data folder are bound read-only, each as
    given and as its real path (a raw file reached through a link is visible both ways);
  - parents are bound first. The session usually sits inside the raw-data folder, and a later
    read-only bind of that folder would cover the read-write output folder below it;
  - DIA-NN's library and report are seen through symlinks in `skyline_work/search/`. Skyline finds
    the report beside the library by name, and nothing is written into the search folder.
- **Protein accessions come from the FASTA names** (`sp|P02768|ALBU_HUMAN`): offline, Skyline
  leaves its accession and gene columns empty.
- **Harmless noise in the logs:** `INFO: Terminating squashfuse_ll after timeout`,
  `libgcc_s.so.1 must be installed`, and a Wine assertion at exit. None of them changed an exit
  code or an output.

**Measured times (trial, Skyline 26.1).**

| Data | Build | Import, per raw file | Join + exports |
|---|---|---|---|
| Thermo `.raw` | 23-27 s | about 40 s (0.7 GB file), 2.3-2.5 GB | about 1 min |
| timsTOF `.d` | about 47 s | 2.5-4 min, about 2 GB per task | about 2.3 min |

Plus queue waits. A 12-run timsTOF project should take about 10-15 min of wall time.

## What it delivers (`<S>/output/skyline/`)

- `peaks_<gene>.png`, one per protein. Its 3 strongest precursors (rows) in each run (columns, at
  most 12, shared across the groups). Each coloured line is a fragment ion. The grey band is the
  peak Skyline integrated; the dashed lines are DIA-NN's peak edges. An empty panel says why
  ("DIA-NN: not reported in this run", "Skyline: no peak").
- `peaks_disagreements.png`: the precursors with the largest apex difference, when there are any,
  in the runs where Skyline disagreed first (at most 12).
- `areas_<gene>.png`: fragment peak areas per run, stacked by ion, with Skyline's dotp above each
  bar.
- `skyline_precursors.tsv`: Skyline's precursor report (apex, peak edges, areas, dotp, idotp,
  mass error, ion mobility).
- `disagreements.tsv`: the table above, one row per flagged precursor-run. The report's median
  apex difference is over the precursor-runs where Skyline picked the same peak, so a flagged 50 s
  never moves it.
- `<name>_top_hits.sky.zip`: the Skyline document with every run's results. It opens in Skyline
  26.1 or later and can be uploaded to Panorama. It is shared "minimal": the library and cache are
  cut to the checked proteins. On the HeLa test it was 3.2 MB against 12.5 MB for "complete", and
  reopened with the same 348 precursor rows.
- `skyline_provenance.json` (schema `skyline_peak_check/1`), the record. It holds:
  - the status (`ok`, `skipped` with the reason and stage, or `planned`), and the plan's id: a job
    writes the record only while it is its own plan's;
  - the sha256 of the `de_provenance.json` the proteins were chosen from;
  - the tool, its container, sha256 and how it was checked;
  - the tolerance and its source;
  - N and its source, and the selection mode;
  - every setting and where it came from;
  - the runs imported and failed;
  - the results, the figures with their captions, and the Methods paragraph.

**What the report says for each state** (`record_state`): "Not run." with the reason (skipped,
or planned but never submitted); "Not finished when this report was made." (submitted, no result
yet); "Out of date." (the DE was run again); "Could not be read."; otherwise the section above.
`MANIFEST.txt` says the same, with what to run.

**Where it goes:**
- The report adds a "Peak check (Skyline)" section from the record, after the top-protein tables,
  whatever the outcome. The report's own section of that title wins, as for the set tests.
- `make_methods.py` adds the "Peak check (Skyline)" paragraph to the Methods when it ran.
- `provenance.py` copies the record into `reproducibility/inputs/` and summarises it in
  `run_manifest.json` (`peak_check`).
- finalize writes its `MANIFEST.txt` line.
- `deliver` ships `skyline/`.
- `skyline_work/` (caches, per-file imports, job scripts, logs) is scratch (`scratch_files.py`): it
  is kept on disk but never zipped or delivered.

**How to read it** (for the report and the user): a section with few disagreements says DIA-NN's
peaks look as expected in the raw data. A protein whose peptides have several disagreements, or a
"Skyline: no peak" in one group only, is one to look at in the pictures (or in Skyline itself)
before relying on its number.

## Targeted follow-up list (PRM / SRM): Core staff, on request

Only when a Core staff member asks ("make a follow-up assay for the top hits"). It needs the peak
check's record to be `ok`; `prm` refuses anyone not on the Core staff list (`staff.py`; fail
closed).

```
python3 scripts/skyline_peaks.py prm --session <S> --instrument bruker_timstof,thermo_fusion \
    [--same-gradient] [--run-length <minutes of the NEW method>]
sbatch <S>/output/skyline/skyline_work/sky_prm.sbatch          # what it prints under "submit with"
```

**What the job does:**
1. **Refine** the peak-check document from the measured areas: the 3 best peptides, 6 best
   fragments and one charge per protein (`--refine-max-peptide-peak-rank=3
   --refine-max-transition-peak-rank=6 --refine-max-precursor-only`, to a new document; the
   original is untouched).
2. **Switch DIA -> PRM** (`--full-scan-acquisition-method=PRM --full-scan-precursor-isotopes=None`).
   Exported from the DIA document, an isolation list holds the DIA windows, not the peptides (the
   trial).
3. **Export one list per instrument:**

| `--instrument` | List | SkylineCmd |
|---|---|---|
| `bruker_timstof` | prm-PASEF isolation list, with RT and 1/K0 windows from the search's runs; **collision energy left at 0 for the method to set** | `--exp-isolationlist-instrument="Bruker timsTOF"` |
| `thermo_fusion` | PRM isolation list (Fusion, Lumos, Eclipse) | `"Thermo Fusion"` |
| `thermo_qe` | PRM isolation list (Q Exactive family) | `"Thermo Q Exactive"` |
| `thermo_stellar` | PRM isolation list | `"Thermo Stellar"` |
| `thermo_srm` | SRM transition list (TSQ Quantiva / Altis): the "p" (precursor-as-product) ion type dropped **before** the transitions are chosen, collision energy predicted ("Thermo TSQ Quantiva") | `--exp-translist-instrument="Thermo Quantiva"` |

- **Exploris and Astral:** SkylineCmd 26.1's `--help` offers no isolation list for them, only a
  native method export, which needs a method template from the instrument PC. `prm` refuses those
  keys and says so.

**Retention-time scheduling:**
- A list is scheduled only when the target instrument is the **same family** as the search
  (timsTOF / Orbitrap) **and** staff confirm the **same gradient** (`--same-gradient`). The window
  is Skyline's default around the measured time; there is no option to change it.
- Otherwise it is unscheduled, and its README says why. An unscheduled Thermo list needs
  `--run-length`: the new method's gradient length, never the search's.

**Outputs** go to `<S>/output/skyline_followup/`:
- the lists;
- a `README.md` (each list, its rows, scheduled or not and why, the collision-energy notes);
- `followup_assay.json`.

The folder is not in a collaborator delivery by default; staff hand the lists over. Each list
still needs a test injection on the target instrument.

## Where the numbers come from

The feasibility trial (2026-10-20, on HIVE) ran the steps above by hand on the two standards
searches this step was then tested on: a timsTOF HT NIST plasma search and an Orbitrap Fusion Lumos
HeLa search. Its findings:

- **Peak edges:** 90.5% of DIA-NN's peaks were matched within 1.5 s at both edges, and 98.4%
  within 5 s.
- **Different peaks:** about 2% of timsTOF precursor-runs had an apex more than 6 s from
  DIA-NN's: Skyline picked a different peak.
- **Areas:**
  - MS1 area r = 0.94-0.99;
  - fragment area against DIA-NN's quantity r = 0.89-0.91;
  - run-to-run changes on HeLa r = 0.97.
- **Failures:** two, both from how Skyline was run (two readers on one `.d`; the process-count
  option).
