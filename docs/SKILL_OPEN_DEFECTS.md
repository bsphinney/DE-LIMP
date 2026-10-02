# ucdavis-proteomics-core-pipeline — open defects

**Audited:** 2026-08-17 against `c29f269`, re-checked 2026-08-20 against `a735fbb`.
**Scope:** `skill/ucdavis-proteomics-core-pipeline/`. Every claim below was checked against the
tree, not recalled — file and line references are from that commit, and the audited files are
byte-identical between `main` and the `skill-v2.0.0` branch, so both carry what follows.

This file is for **open** defects. `docs/GOTCHAS.md` is the quick reference for problems that
already have solutions; when something here is fixed, delete its entry and, if the lesson is
reusable, add a row there instead. A shrinking file is the point.

---

## Status at a glance

**Open (still, after 2.10):**
- a Sage LFQ window that does not fit the MS1 mass error: gated in 2.10 with the corrected
  re-run command; the automatic re-run is still backlog (below);
- the protein rollup of a Sage report is max(peptide). Since 2.9 it is described as exactly
  that; a real rollup is still missing;
- a multi-file `.raw` Radiant chain converts nothing;
- Sage `.d`: no route gives LFQ on HIVE;
- upstream Sage stores one LFQ q per target/decoy pair (issue text written, not filed);
- AlphaDIA's and Radiant's protein ids are read as bare accessions without a real output to
  confirm it (below).
Details below.

Also for 2.10, from sage-review's check of b1b94f5 (2026-09-30):
- **N2: an adapted report made before 2.9 must be re-adapted before its DE is re-run.** It
  declares no quantity level and names no peptide. A DE re-run on it still says "MaxLFQ + limma"
  / "DIA-NN PG.MaxLFQ", as for a DIA-NN report. Since this fix its contaminant counts say
  "report rows", and which sample proteins lost rows is "not computed". Run the search's own
  `run_search.py` command again with `--adapt-only` first (it takes the search's usual required
  arguments and only reads the existing output), then the DE. A guard that refuses an
  undeclared report whose folder holds no DIA-NN log would close it.
- **N6: `Protein.Names` is not position-aligned with `Protein.Group`.**
  `protein_ids.group_entry_names()` (b1b94f5) drops repeated names, and members whose header
  has no `db|ACC|NAME` entry name. `group_accessions()` drops repeats separately. So name *i* is
  not always the name of accession *i*. `sp|P02769|ALBU_BOVIN;sp|Cont_P02769|ALBU_BOVIN` gives
  `Protein.Group` `P02769;Cont_P02769` but `Protein.Names` `ALBU_BOVIN`, and a proteogenomics
  header with no entry name drops out of the names list. No harm today: `make_analysis_html`,
  `analysis_prompt`, `audit_results` and build_maxlfq's annotation use only the whole string or
  its first member. **Fix (2.10):** write one name per accession, "" where there is none, with no
  dedup, as DIA-NN's own `Protein.Names` does. Otherwise, document that the two lists are not
  aligned.

From keratin-fixer's real-data test of 979dfd5 (a Q Exactive Plus hair DDA set, 2026-09-30):
- **A sample protein seen only through peptides it shares with a non-keratin contaminant is
  removed entirely.** When every detected peptide of a human protein is shared with a `Cont_`
  entry that has unique peptides of its own (ALDOA and rabbit ALDOA, GSN, CYCS), all of them go
  with the contaminant, and the protein leaves the DE. The count is disclosed in the methods.
  **Consider (2.10):** naming those proteins in the report, not only counting them.
- **The Q Exactive Plus is not in `make_methods`' acknowledgment registry** (`ACKS`: the
  Fusion Lumos, the Exploris 480 and the timsTOFs), so its Methods carry a placeholder
  acknowledgment. Pending Brett's answer on which grant, if any, covers it.

From notes-reviewer's final verdict on `feat/skill-2.9-notes` 8e4da34 (SHIP, 2026-09-30), for
2.9.x:
- **A stray name in the inbox can still carry a few words into Brett's Claude.** `notes.py`
  `safe()` quotes a stray entry's name as `ascii(name)[:40]`, so a member who plants a file in
  `skill_notes/` can put about 38 characters of their own words into what `replies` and `list`
  show. **Fix:** show strays as counts by kind only, never their names.
- **Two name-rule edge cases (the review's N2/N3)** need the right to rename a note inside its
  folder, which the folder rule (sticky, owned by an admin or the sender) now blocks; they stay
  documented, not fixed.
- **Pre-planted receipt names can block an ack.** `ack` writes `<id>.read.<user>.<ts>` (the
  second it runs) with `os.replace`, which in a sticky folder cannot replace a file someone else
  owns. A member who pre-plants those names (about 86,000 per note per day) makes the
  recipient's `ack` fail. **Fix:** create receipts with `link()` / `O_EXCL`, and on a taken name
  retry with the next timestamp.

From the 2.9 integration's check of `fix/2.9-coreomics-key` 9fb14df (2026-09-30), for 2.10:
- **The CoreOmics key-save line is quoted in three docs, and no test keeps the quotes equal to
  the code.** `core_submission.SAVE_KEY_LINE` is the one definition; SKILL.md,
  `references/access.md` and `references/core-submissions.md` quote it word for word (all four
  quotes matched at 9fb14df). An edit to one of them alone would leave staff saving the key with
  an old line. **Fix:** a test in `test_core_submission.py` that finds every `read -rs TOK`
  line in the docs and asserts it equals one of `SAVE_KEY_LINE`'s values.

The 2.10 backlog from the DIA-NN DDA review (2026-09-29) is listed directly below: these
are measured leads to follow up, not bugs in shipped behaviour.

### 2.10 backlog — DIA-NN on DDA (proteomics review of `fix/2.9-dda-chain`, 2026-09-29)

Evidence: the SET28 test chain on HIVE (6 Exploris 480 DDA hair runs, DIA-NN 2.7.0, jobs
24197069–24197073, `/quobyte/proteomics-grp/claude/dda_chain_test_20260929/`) and msalemi's
28-run SET28 chain. Do not implement these without the A/B each one names.

1. **`--fwhm` into steps 4 and 5, behind an A/B.** Step 3 and step 5 both log
   "INFO: consider re-running DDA quantification with --fwhm 0.175602 for optimal results"; the
   DIA-NN README says `--fwhm` "may affect quantification, in particular in DDA mode". Proposal:
   step 3 writes the mean FWHM it logs to `<out>/fwhm.txt`, steps 4 and 5 pass
   `--fwhm $(cat fwhm.txt)`. Ship it only if an A/B on the same runs shows it helps. On a
   Q Exactive Plus hair DDA set (2026-09-30) DIA-NN suggested `--fwhm 0.023`: the value depends
   on the instrument and gradient, so it has to come from the run, never a constant.
2. **A `--fix-scoring` A/B on the same 6 SET28 runs.** The README: "in DDA mode, during second MBR
   pass or when analysing with empirical library, DIA-NN may auto-switch to Proteoforms
   regardless of the initial Scoring setting". The final pass lost more than half of the
   first-pass precursors in 4 of 6 runs (`pass_comparison.py` now reports this); test whether
   `--fix-scoring` changes that.
3. **F6 — the timsTOF DDA range from the precursors actually isolated.** A ddaPASEF run's
   precursor range is now GlobalMetadata `MzAcqRangeLower/Upper` (99.99–1700 on one of the SET1-28 runs), much
   wider than what was picked (PasefFrameMsMsInfo IsolationMz 273.75–1699.2). A tighter,
   data-derived bound from `PasefFrameMsMsInfo` would shrink the predicted library; measure the
   library size and IDs both ways first.
4. **A registry flag for test runs.** A test chain's job-end hook still logs the run in
   `/quobyte/proteomics-grp/skill_runs` (`--no-notify`/`--no-fran` do not stop the run log): the
   SET28 test landed as `skill_runs/sessions/2026-09-29_session`. Add a `--test-run` (or
   `RECORD_RUN=off` baked in at generation) so validation runs stay out of the Core's registry.
5. **Chain steps 3 and 5 ask for more than the per-user CPU cap leaves** (the staff-issue
   backlog's #8: the chain ignores the per-user 64-CPU cap). Each requests 64 CPUs, 128 GB and
   12 h, so while other jobs of the same user hold part of the cap they wait in
   `QOSMaxCpuPerUserLimit`. On the Q Exactive Plus hair set they finished in 22 s at 32 CPUs.
   Size them from the cohort, with headroom under the cap.
6. **The probe fallback, from dda-review's pass on probe-estale** (none blocking; its N1, a
   stale `probe_fallback.json`, is fixed in 2.9.0):
   - **N2: an out-of-band value from a failed run is ignored** when two other runs pin the mass
     accuracy in band.
   - **N3: one answered run plus a transient `io_error` fails without a retry.**
   - **N4: a DIA-NN killed for memory (OOM) reads as "not measured".** DIA-NN's return code
     is not recorded, so an OOM-killed run is treated as a data failure.
   - **N5: a repeatable crash in the probe would make every search fall back.** If the CAUTION
     starts appearing on every run, treat it as a probe bug, not bad luck. **Consider:** a
     registry-level alarm that fires when fallbacks exceed a rate.
   - **`set_aside`'s refusal says "not a job script"** even when the path is a probe output
     that is not a regular file (`diann_parallel.set_aside`). Cosmetic: name what was refused.
   - **`watch_run.sh --all` and `fallback_modes` read a fallback differently.** `watch_run.sh
     --all` ORs the `probe_fallback` record with the mode, while `fallback_modes` lets the mode
     win; they disagree only on an inconsistent provenance no writer produces. **Fix (2.10):**
     unify by matching only the `scan_window` / `mass_acc` mode vocabulary, not any `"mode":`
     (a first attempt, 3a37ff0, matched every `"mode": "parallel_5step"` and was not shipped).
7. **The chain's recovery notes can land outside the session.** `submit.sh` writes
   `RECOVERY.md` and `.recovery.json` two folders above the output folder
   (`diann_parallel.py`: `sess = os.path.dirname(os.path.dirname(out))`), which is the session
   root only for the usual `<session>/output/search` layout. For any other output folder they
   land in an unrelated directory: dda-fixer's HIVE chain test (2026-09-30) wrote them into
   Brett's `Data/lab/4brett/`. **Fix (2.10):** write them to the session root when there is one
   (a `session.json` above the output folder), else to the output folder itself.

Every entry from the 2026-08-17 audit has been fixed; those are listed below and in
`docs/GOTCHAS.md`. The run-restriction defect was found afterwards, by an end-to-end
confirmation run of the #48 fix, and is fixed as of skill v2.1.1.

---

## Different files that share a run name cannot be searched together (2.11)

**Decided for 2.10 (option A):** the list stages (`ht_manifest.py`, `core_submission.py locate`)
keep every file and flag the shared name (`repeated_names`, with each path), but the engine
step stops before writing anything (`check_report_runs.names_stop`, in `run_search.py` for every
engine, `diann_parallel.py` and `radiant_parallel.py`). Passed as given, DIA-NN names a run by
its file name without the folder, so `/plate1/s1.d` and `/plate2/s1.d` become ONE run (two
samples summed into one column; in the 5-step chain, two array tasks write one `.quant`), and
Sage converts both to one mzML name. Ways out today: keep one (`locate --reinjections latest`
for re-injections of one sample), rename one, or search them separately.

**2.11 (option B):** search both under unique run names -- `run_search.py` links each under a
name that cannot collide (`plate1__s1.d`), searches the links, and records the original path ->
run name map in `search_provenance.json` and the analysis report; conditions and the report
check then use the new names. Mind the container routes, whose bind mounts are derived from the
input paths (the link and its target must both be bound). **`locate --reinjections all` with
same-named re-injections needs this**: until it lands, such a list hits the clear stop above,
never a merge.

## Sage LFQ with an MS1 offset near its ±5 ppm window: gated in 2.10, not re-run automatically (backlog)

**What.** Sage integrates MS1 for LFQ only within ±`quant.lfq_settings.ppm_tolerance` of the
theoretical mass (default 5.0; Sage 0.14.7 `crates/sage/src/lfq.rs` `build_feature_map`).
`estimate_params.py` never sets it, and the search's `precursor_tol` is ±10 ppm. So a cohort
whose MS1 sits a few ppm off identifies normally and quantifies badly. gabrig 2026-09-29
(HeLa, Fusion Lumos, `skill_issues/2026-09-29_gabrig_HeL50_UnvPe_HCDOT_28sep26.md`):
- at +7.1/+7.7 ppm: 29,846 PSMs, but "discovered 0 target MS1 peaks at 5% FDR";
- at +1.8 to +5.8 ppm: 17,288 peaks, protein CV 50% against MSFragger's 14%.

**Done in 2.9.** `sage_lfq_check.py` runs after every Sage search (in the job, and again at
`run_search.py --adapt-only`). It reads each run's median `precursor_ppm` and Sage's
MS1-peak count, and warns when a median + 2 ppm exceeds the window, or when the count is 0 or
under 10% of the target peptides. The warning names a window that would fit and changes
nothing. It reaches `search_provenance.json`, `watch_run.sh --out`, `checkpoint.py status` and
`audit_results.py --search-out` (AUDIT.md, and from there the report).

**Done in 2.10: a gate with the way through.** When the check warns because of the mass error,
it writes the corrected window into a copy of the config (`<out>/sage_config.lfq_ppm<N>.json`)
and records the exact `run_search.py` command that repeats the search with it into
`<out>_lfq_ppm<N>`. `--adapt-only`, an inline search and `run_de.R` then refuse the quantities
until that re-run is used or the user accepts them (`--accept-lfq-window "<who, why>"`,
recorded). Tests: `tests/test_sage_lfq_gate.py`.

**Backlog: the automatic re-run.** When the check warns *because of the mass error* (not a
low peak count with a fitting window, whose cause is unknown):
1. write a copy of the Sage config with `quant.lfq_settings.ppm_tolerance` =
   `suggested_ppm_tolerance`, or with a per-run recalibration if a single window cannot cover
   the spread of per-run offsets (gabrig's cohort spans +1.8 to +7.7 ppm);
2. re-run LFQ as a new job (Sage has no LFQ-only mode, so this is a full re-search);
3. re-check, and adapt from whichever run passes.

It needs a decision first. A wider window also admits more co-eluting interference, so the
target/decoy separation of the re-run must be compared, not assumed better. And a changed
search setting is a new search to confirm (golden rule 1), not a silent retry. Validate it on
gabrig's HeL50 cohort against the MSFragger CVs before shipping.

---


**What.** `--method dpc` reads the entire report, then subsets:

```r
# scripts/run_de.R:248-267
dat <- limpa::readDIANN(dpc_input, format = format, q.cutoffs = ..., q.columns = q_use)
...
.keep <- .rep_runs %in% meta$File.Name
if (any(!.keep)) dat <- dat[, .keep]        # <-- columns dropped AFTER the matrix is built
```

`readDIANN()` builds the precursor matrix over **all** runs in the report, so the precursor set
— and therefore the protein set that `dpcQuant()` rolls up — is decided using runs the metadata
excludes. Subsetting afterwards drops those columns but keeps the rows they justified.

`build_maxlfq.R` does the opposite, and is right: `keep_runs` is applied to the arrow query at
line 88, **before** `collect()` and the pivot, so the protein set is decided on the analysed
runs only. **This is the same class of defect as #48 — the two `--method` values disagreeing
about the same report — this time about which runs define the protein set rather than which
q-columns are applied.**

**Measured.** Found by the end-to-end confirmation run of the #48 fix, which tested 6,629
protein groups where an otherwise identical arm tested 6,531. With the q-columns now equal, that
gap is entirely this: the 6,531 set is an exact subset — **98 extra, 0 lost**.

The report has 399 runs; the analysis used 373 (22 UE, 2 DDA and 2 others excluded by design).

| protein groups | median peptides | median fraction of analysed runs observed | called significant (BH < 0.05) |
|---|---|---|---|
| in both sets (6,531) | 12 | 25.6% | 71.4% |
| **admitted only by the 399-run read (98)** | **5** | **6.8%** | **40.8% — 40 proteins** |

Eleven of the 98 are observed in ≤1% of the analysed samples. **Forty of them reach the results
table**, tested on data that barely contains them.

**Why it matters.** Excluding runs from an analysis is a deliberate act — here a fourth
preparation present in only one cohort, plus two DDA acquisitions inside a DIA search. The
exclusion should mean those runs cannot influence the result. At present they still decide which
proteins are eligible to be tested, which is a quieter kind of influence than contributing data
and harder to notice.

**Proposed fix.** Restrict before the protein set is decided, matching the maxlfq path: pre-filter
the report to `meta$File.Name` and hand `readDIANN()` the filtered input, as the session wrapper
that produced 6,531 did. If pre-filtering the input is awkward for the tsv route, the fallback is
to drop precursor rows with no observations in the retained runs immediately after `dat[, .keep]`
— narrower, and it would not catch precursors that are merely sparse rather than absent.

**Not urgent for existing results.** Any workflow whose report contains only the runs being
analysed is unaffected, which is the common case. It bites exactly when a report deliberately
contains more runs than the design uses.

---

---

## ~~1. The precursor *m/z* range is hardcoded~~ — FIXED 2026-08-17

Fixed in `fix/skill-precursor-mz-range`. `detect_acquisition.py` now returns the acquired
bounds it was already reading (`precursor_mz_range`), and `estimate_params.py` takes
`--precursor-mz-range LO HI` and rounds outward. When no range is available the old 380–980 is
still used but tagged **`FALLBACK — acquired range unknown for this input; NOT measured`**
rather than `universal trypsin/LFQ default`, so a reader can tell a measurement from a guess.

Verified against a real `.d` on HIVE: acquired 299.5–1200.5 → emits `--min-pr-mz 299`
/ `--max-pr-mz 1201`. Covered by `skill/ucdavis-proteomics-core-pipeline/tests/`, now run in CI.

Kept here as a stub only until the next edit of this file; the lesson is in `docs/GOTCHAS.md`.

---

## ~~2. The q-value column set is duplicated~~ — FIXED 2026-08-17

Fixed in `fix/skill-q-column-definition`. The set now lives in
`scripts/diann_q_columns.py`, mirrored by `scripts/diann_q_columns.R` (the two R scripts are
standalone Rscripts that treat `jsonlite` as optional, so reading a shared JSON file would put
a hard dependency in front of the identification filter). `tests/test_q_columns.py` parses the R
source and asserts it equals the Python values, so drift is a CI failure rather than two subtly
different FDR filters.

The audit undercounted: `run_search.py` alone held **six** hand-written copies, not one — four
of them in output adapters (AlphaDIA, Sage, FragPipe DIA, FragPipe combined_protein,
Radiant/Fulcrum) that emit the DE contract. A test asserts no call site restates the set inline,
which is what found them.

It also conflated **two distinct concepts**, and merging them would have been a bug:

* `FDR_REQUIRED` + `FDR_OPTIONAL` — the **filter set**: every column ANDed at the q-cutoff
  (`run_de.R`, `build_maxlfq.R`).
* `PROTEIN_Q_PREFERENCE` — a **preference chain**: the first available protein-level q-column
  wins and the rest are ignored (`compare_searches.py`, `run_search.py`).

Applying the preference order as a filter would over-filter; filtering on only the first
available column would under-filter.

---

## ~~3. A single `--q-cutoff` for all six columns~~ — FIXED 2026-08-17

`COLUMN_CUTOFFS` in `scripts/diann_q_columns.py` (mirrored in the `.R`) gives `PG.Q.Value`
DIA-NN's own recommended **0.05**; the other five keep `--q-cutoff`. `limpa::readDIANN()`
recycles `q.cutoffs` against `q.columns` element-wise, verified in
`EListFromLongFormatFile`, so the dpc path passes a vector; `build_maxlfq.R` applies the same map
and labels any column whose cutoff differs as `PG.Q.Value@0.050` in `filters_applied`.

**This is a behaviour change that LOOSENS identification FDR.** Measured on a 2-run HeLa report:
32,559 → 32,741 rows and 2,695 → 2,715 protein groups. `cutoff_for(..., uniform=True)` restores
one cutoff for all six.

Deliberately *not* clever about a tightened `--q-cutoff`: the pipeline default (0.01) is already
stricter than 0.05, so a "never loosen what the caller asked for" rule would stop the
recommendation ever applying. `--q-cutoff` governs the other five.

---

## ~~4. `--min-pr-charge 2` narrows DIA-NN's default without saying why~~ — FIXED 2026-08-17

Tag only; no behaviour change. It now reads: *"z=1 excluded: rarely informative for tryptic
bottom-up. DIA-NN's own default is 1-4; this drops ~19% of the predicted library (measured:
10,899 → 8,805 precursors on a 60-protein FASTA)"*. The narrowing is still applied — it is the
right call — but a reader of the emitted parameters can now see that it was a decision.

---

## ~~The dpc path decides the protein set before restricting to analysed runs~~ — FIXED 2026-08-20

Fixed in `fix/dpc-restrict-runs-before-rollup`. `run_de.R` now applies
`Run %in% meta$File.Name` to the arrow query **before** `readDIANN()`, matching what
`build_maxlfq.R` has always done, instead of subsetting columns afterwards.

**The audit understated it.** It measured the protein-set effect (98 extra, 0 lost). It did not
measure what happens to the proteins that were *already* there. Reproduced on a 12-run report
analysed at 7 runs:

| | before fix | after fix |
|---|---|---|
| all-NA precursor rows surviving the subset | 158 | 0 |
| extra protein groups vs a pre-filtered report | 3 (exact superset) | 0 |
| **shared proteins with a different `logFC`** | **4,645 of 4,645 (100%)** | 0 structural |
| **significance calls flipped** | **54** | **0** |

So the excluded runs were not only deciding which proteins were *eligible* — those all-NA rows
reach `dpcCN()`/`dpcQuant()`, so they shifted the detection model **every retained protein was
quantified against**. On this data the extra proteins were not themselves significant; the harm
was 54 flipped calls among proteins that would have been analysed either way.

A residual numeric difference remains vs a natively pre-filtered report (median |Δ logFC|
1.1e-07, max 3.0e-03, **0** significance flips), consistent with row ordering changing an
iterative fit's path rather than anything structural. Before the fix the same comparison was
median 0.0021 / max 0.11 / 54 flips.

Guarded by `tests/test_reproducibility_contract.py::TestRunRestrictionHappensBeforeRollup`,
verified to fail when the restriction is removed. The post-hoc `dat[, .keep]` is retained as a
backstop for the tsv route, which cannot pre-filter.

---

## Upstream Sage: the stored LFQ `q_value` is shared by a target and its decoy (issue to write up, NOT filed)

**What (Sage v0.14.7 source):**
- `crates/sage/src/fdr.rs` `picked_precursor()` computes a q per (PrecursorId, decoy) entry
  and counts targets at `q <= 0.05` (the "discovered N target MS1 peaks at 5% FDR" line).
- It then builds `scores: FnvHashMap<PrecursorId, f32>` keyed by the precursor ONLY, from the
  score-sorted list, so the later (lower-scoring) member of each target/decoy pair overwrites
  the other.
- It then writes `peak.q_value = scores[ix]` into BOTH entries.

So `lfq.parquet` / `lfq.tsv` report, for every target, the q of whichever of it and its decoy
scored lower. That is always ≥ its own q. The file cannot reproduce Sage's own logged count.

**Measured:** gabrig's HeL50 UnvPe `lfq.parquet` (2026-09-29): identical q in 1,489 of 1,489
target/decoy pairs. Filtering targets at stored q ≤ 0.05 gives 16,649 against Sage's logged
17,288: 639 (3.5%) that Sage counted are lost. It errs conservative, so there is no FDR harm.

**Issue text, for Brett to file at github.com/lazear/sage (not filed):**
- Title: "picked_precursor writes one q_value per PrecursorId into both target and decoy LFQ
  peaks".
- Body: the numbers above, the three lines of `fdr.rs`, and the fix. Key the map by
  `(PrecursorId, bool)`, the same key `peaks` uses, so each entry keeps its own q.

**Skill side (2.9):** the adapter keeps the conservative subset and records Sage's logged count
beside `target_precursors_kept`. `sage_lfq_check.ms1_peaks_from_lfq()` is labelled a lower
bound.

## The protein rollup of an adapted Sage report is max(peptide): LABEL fixed in 2.9, rollup backlog

**What.** `adapt_sage` emits one row per **peptide** × file (Sage quantifies peptides).
`build_maxlfq.R` takes `max(PG.MaxLFQ)` per (Protein.Group, Run), so a protein's value is its
single most intense peptide in each run, not MaxLFQ. Until 2.9 the DE nonetheless called it
"MaxLFQ + limma" / "DIA-NN PG.MaxLFQ" / "Quantification: DIA-NN MaxLFQ (Demichev 2020)". That
wording reached methods.txt, de_provenance.json, the Methods, the AI brief and
reproducibility_log.R (architectural rule 1).

**Fixed in 2.9: the description.**
- The adapted report declares its quantity (parquet metadata `delimp.quantity.*`), and
  `build_maxlfq.R` `maxlfq_descriptor()` describes it: "Highest-peptide intensity + limma",
  pipeline `peptide_max`, a rollup naming the highest-peptide rule, and Sage's own FDR.
- Its `caveat` becomes a `quantification` WARN in AUDIT.md (and the report).
- Every consumer reads the descriptor, and none branches on "sage".
- `tests/test_sage_de_describes_itself.py` fails if any of them says "DIA-NN" or "MaxLFQ" for
  a Sage DE.

**Backlog: a real peptide → protein rollup.** Two options:
- MaxLFQ over the Sage peptide rows (e.g. `iq::maxLFQ` or `limpa`'s own);
- feed Sage peptide intensities to limpa's DPC path (`dpcQuant` on a peptide matrix), which
  would also model the missingness.

The adapter would then declare a level the descriptor turns into that method's name. Validate
against MSFragger/IonQuant on gabrig's HeL50 UnvPe. She reported protein CV 50% for Sage against
MSFragger's 14% there; her rollup for that is in `sessions/comparisons/harmonised_compare.py`.

## Contaminant ids from the other adapters: what 2.9 checked and what it did not (open)

- **AlphaDIA and Radiant ids.** AlphaDIA's `pg` comes from alphabase's FASTA parsing (UniProt
  accessions), so `adapt_alphadia` is unchanged. `adapt_radiant` now passes Fulcrum's groups
  through `protein_ids.group_accessions`, which leaves a bare accession as it is. Neither has
  been checked on a real output whose FASTA carries `Cont_` entries.

## Multi-file `.raw` Radiant chain: no conversion at all (backlog)

`run_radiant_parallel` / `radiant_parallel.py` (Radiant, >1 file, on SLURM) hand the input list
straight to Radiant, which reads only mzML or Parquet. A `.raw` chain therefore fails in its
step-2 array. It fails loudly, which is acceptable for now. Only the single-job `--sbatch` route
converts (in its job, since 2.9). **Fix:** plan the conversions with `run_search.mzml_plan()`
and run `conversion_lines()` in the chain, either as a step before the array or per array task
for its own file.

## Sage + Bruker `.d`: no route gives LFQ (backlog)

`mzml_plan` sends `.d` to msconvert, but bioconda's Linux msconvert has no vendor readers (built
from `pwiz-src-without-v`), so it cannot open a `.d` on HIVE. Sage reads `.d` natively since
v0.14.4, but "no MS1/LFQ support yet" (Sage CHANGELOG v0.14.4). The v0.14.7 `tdf.rs` emits only
`ms_level: 2` spectra, so there is no MS1 to quantify. Skipping msconvert for `.d` would search
fine and **lose all quantification**. **Options:** a vendor-enabled msconvert (ProteoWizard's
wine container, as DE-LIMP's `$MSCONVERT_SIF`) for LFQ; native `.d` for identification-only
searches; or DIA-NN `--dda` for timsTOF DDA.

## Fixed since the last audit — recorded so they are not re-litigated

- **`adapt_sage` passed Sage's decoy MS1 peaks, and targets failing Sage's own 5% line, into
  the DE** — fixed in 2.9. `lfq.parquet` (sage-cloudpath `serialize_lfq`) keeps each target's
  decoy, with the target's `proteins` string, at any `q_value`. Sage's own `lfq.tsv` writer drops
  the decoys (`write_lfq`), and Sage counts a peak only at its own `q_value` ≤ 0.05 (`fdr.rs
  picked_precursor`, compared in f32). The adapter now keeps non-decoy rows whose STORED q ≤ 0.05
  and records the counts per file (`sage_adapt.json`). The stored q is shared with the decoy
  (below), so this is a conservative subset of Sage's count, and both are recorded. Measured read-only on gabrig's HeL50 Sage outputs
  (2026-09-29, Flinders `deLimpClawd_28sep26/sessions`; no DE had been run on them):
  - **UnvPe**, 4 runs, well calibrated: 1,561 decoy + 1,620 failing rows per file. With them
    dropped, 19–20% of proteins per file change their summed intensity by more than 10%.
    158–236 protein×run cells per file took their value from a decoy or failing row, and 177 of
    3,765 protein groups existed only through them.
  - **HCDOT**, +7 ppm: 2,507 protein groups, 100% from decoy or failing rows.

- **No Sage DE ever had a contaminant removed, and a keratin sample's kept keratin was
  described as DIA-NN's for every engine** — fixed in 2.9 (post-merge).
  - The filter's rule (contaminants.R: an accession *starts* with `Cont_`/`contam_`) tests bare
    accessions. `adapt_sage` wrote Sage's `proteins` as full FASTA ids
    (`sp|Cont_P02769|ALBU_BOVIN`), so nothing matched. On gabrig's HeL50 UnvPe the Sage DE
    limma-tested 171 `Cont_` groups (BSA, gelsolin, mouse keratins) while methods.txt said
    "Contaminants : none", and a keratin sample's exemption could never apply.
  - `adapt_sage` (and `adapt_radiant`) now write bare accessions through
    `protein_ids.group_accessions`, with the entry names in `Protein.Names`. FragPipe DDA was
    already right after drift's `Protein ID` change.
  - The keratin lines (methods.txt, make_methods, the AI brief) now quote each pipeline's
    `kept_contaminant_quant`. DIA-NN's `--cont-quant-exclude` is named only for DIA-NN
    reports; a Sage, FragPipe, AlphaDIA or Radiant DE says what its own engine did.
  - An older FragPipe DIA output searched on a Philosopher `--contam` database is WARNED about
    when adapted: its DIA route drops the tag, so those contaminants would be tested.
  - **Counts and units (sage-review B1).** b1b94f5 counted an adapted report's item × run rows as
    "precursors": its Sage DE printed 2,744 for 686 peptides × 4 runs, its FragPipe DE 288 for 72
    protein groups × 4 runs. It also printed "0 lost some" for Sage's mixed groups. Now:
    - the census counts distinct items in one unit (`contaminant_unit`, from the declared
      level): precursors, peptides (`Peptide.Id`, which `adapt_sage` now writes) or protein
      groups;
    - mixed `Cont_` groups are counted, and for a peptide-level report so are the sample proteins
      they named;
    - a protein-level report says "not computed" for what it cannot show, never 0.
  - **The DIA-NN sentence is read, not assumed (N1).** The dpc and undeclared-maxlfq descriptors
    name DIA-NN's `--cont-quant-exclude` only when the run recorded it. The source is the command
    line in the DIA-NN log beside the report, else the parameters file the search ran with
    (`make_methods.diann_cont_quant_exclude`, which `search_record` uses too). Otherwise they say
    `[not recorded -- confirm]`. FragPipe's DIA route is not assumed to have had the flag.
    - **Fixed after sage-review's re-check of ca9aff2:** DIA-NN logs `--cfg <file>`
      unexpanded. This was seen in DIA-NN 1.8.1 (FragPipe), 2.6.1 (taha_dog) and 2.7.0
      (PROT_0002). So every `run_search.py --one-step` search, whose flags sit in the `--cfg`
      file, read "not set" in the DE, the Methods and the run record.
    - The reader now opens each `--cfg` file where it stands on the command line, as DIA-NN
      does. A relative path resolves in DIA-NN's working directory: the log's folder minus the
      folder part of a relative `--out` (FragPipe), otherwise the log's folder (`cd <out>`), then
      its parent.
    - A `--cfg` that can't be read falls back to the parameters file, then NOT RECORDED, never
      "not set".
    - Checked read-only on real HIVE logs: taha_dog's `diann_lib.log.txt` (`--cfg`) gives
      `Cont_` from `params.cfg`; FragPipe's `dia-quant-output/report.log.txt` gives "no flag"
      from its `filelist_diann.txt`.
  - Also from that review:
    - N3: make_methods quotes the record's rule text instead of its own copy;
    - N4: FragPipe/AlphaDIA/Radiant keratin is "under-quantified: not known" (Philosopher's razor
      assignment can move shared peptides), not FALSE;
    - N5: with no keratin kept, the engine sentence is left out.
  - Tests: `tests/test_contaminants_every_engine.py` (8 in b1b94f5, of which 6 fail on 979dfd5,
    one of them on a TypeError from `keratin_sample_check`'s new engine argument, not on an
    assertion). Sixteen after the review: 9 of them fail on b1b94f5, plus
    `test_keratin_sample.py`'s DIA-NN log test.
  - Confirmed on real data (HIVE jobs 24210171 and 24210605, the latter on the final scripts; read-only copies of gabrig's HeL50 UnvPe
    outputs, `/quobyte/proteomics-grp/claude/contam_check_20260930/`):
    - Sage: 686 peptides removed, all 171 `Cont_` groups, 101 of them also naming a human
      protein. Those named 120 human proteins: 57 keep other peptides (36 as their own group,
      21 only inside a group shared with other proteins), and 63 have none left. An independent
      count from the adapted report agrees.
    - FragPipe: 72 protein groups removed (72 `Cont_` groups, none mixed). Whether human
      proteins lost peptides inside IonQuant's quantities is "not computed".
    - Neither contaminant block says "precursor". No tagged group reaches either DE matrix
      (3,414 and 3,164 groups).
  - **The run record read the flag from the FASTA's recommendation.** record_run.py's SEARCH_LOG
    said "excluded from quantification with `--cont-quant-exclude Cont_`" whenever the FASTA
    sidecar recommended it (`diann_cont_quant_exclude`). It did so for any engine, and whether or
    not the flag ran, and the transcript and AI review read that record. It now asks
    `make_methods.diann_cont_quant_exclude` (`search.cont_quant_exclude` in run_record.json):
    - a DIA-NN search says what its log or parameters recorded, or "not recorded";
    - another engine says nothing about the flag;
    - the sidecar's value is kept as `fasta.cont_quant_exclude_recommended`.
    The reader's parameters-file fallback now applies only to a DIA-NN search. run_record.json's
    `schema_version` goes from `2` to `3` (the history is beside `record_run.SCHEMA_VERSION`), so
    records before and after the rename can be told apart.
    Tests: `tests/test_record_run.py` ContQuantExcludeIsReadFromTheRun (5).
- **FragPipe-DDA, AlphaDIA and Radiant DEs were described as "DIA-NN PG.MaxLFQ" / Demichev** —
  fixed in 2.9. Their numbers were right (one engine-made value per protein and run), but the
  description was not. Seen in Michelle's PI_Example SET1-28 `2026-09-25_SET28_DDA_FragPipe23`, whose
  de_provenance.json says "DIA-NN PG.MaxLFQ". Each adapter now declares level `protein` with the
  engine's own quantity, version, FDR and citation, and `maxlfq_descriptor()` passes it through
  (`tests/test_adapters_describe_themselves.py`). To re-render the methods of such a session,
  re-adapt it first (`run_search.py --adapt-only --engine fragpipe`), because its old
  `report.parquet` declares nothing, and then re-run the DE. The quantities do not change.

  Every Sage DE made before 2.9 should be re-run. A search the same day found no other adapted
  Sage report (an `lfq.parquet` with a `report.parquet` beside it). It covered every
  `lfq.parquet` under Flinders `Data/lab/service`, `/quobyte/proteomics-grp` and gabrig's and
  msalemi's homes (srun, depth 11). The other eight are `brett/zach_hair_gvp` identification
  runs that were never adapted. Not searched: other users' homes, and runs made off HIVE.

All verified present in `c29f269`:

- **Experiment-wide FDR columns on the dpc path** — PR #48. `--method dpc` treated
  `Global.Q.Value` / `Global.PG.Q.Value` / `PG.Q.Value` as fallbacks and so never applied them on
  any modern DIA-NN report. Worth 227 protein groups (median one precursor, present in 9.4% of
  runs) on a 373-run report.
- **Executable bits.** Every directly-invoked `scripts/*.sh` is `100755` in the index.
  `diann_release.sh` is `100644` and correctly so — it is `. sourced`, never executed.
- **`--requeue`** is emitted for generated SLURM steps, so a preempted array task on a
  preemptible partition no longer breaks an `afterok` chain hours in.
- **The step-4 shared-library write race.** Each array task now gets a private library copy;
  previously all tasks re-saved the same file and the losers blocked indefinitely.
- **The glob hazard** in generated commands — `shlex.split` / `shlex.quote` now handle
  `--cut K*,R*` and `--var-mod ...,*n`, which would otherwise be expanded by the shell against
  the working directory.
- **The deprecated one-stage DIA-NN invocation.** `run_search.py` detects the library-free flag
  combination, strips it, and emits two dependent SLURM jobs — the two-stage workflow DIA-NN 2.3+
  requires.
