# SILAC with DIA-NN plexDIA (skill 2.13)

SILAC samples are searched with DIA-NN's plexDIA and analysed as **ratios inside each run**: each
run holds two or three samples, one per isotope label, so a protein is compared between labels of
the same injection, not between runs. This page is the detail behind SKILL.md steps 3, 6b and 8.

**Supported:** DIA data, DIA-NN (2.7.0 verified), Lys/Arg SILAC with 2 (L/H) or 3 (L/M/H)
channels, independent samples per channel, label swaps, label-check injections.
**Not supported (refused, with the reason):** Sage or DIA-NN's DDA mode; a seed or label-free
spectral library (`--seed-lib`); chemical labels (mTRAQ, dimethyl); a turnover / pulsed-SILAC
*analysis* (the search runs, `run_silac.R` refuses); protein-set tests on the ratios.

## 1. What to ask (step 3)

Ask; never guess. Each answer goes into `logs/decisions.md` (`log_decision.py`).

1. **The labels.** Which labelled lysine and arginine each labelled channel was grown in. The
   presets (`python3 scripts/silac.py presets`), masses from UniMod (read 2026-10-15):

   | Preset | Lys | Arg |
   |---|---|---|
   | K0R0 (light, always channel L) | 0 | 0 |
   | K4R6 | +4.025107 (2H(4), UniMod:481) | +6.020129 (13C(6), UniMod:188) |
   | K6R6 | +6.020129 (13C(6), UniMod:188) | +6.020129 (13C(6), UniMod:188) |
   | K8R10 | +8.014199 (13C(6)15N(2), UniMod:259) | +10.008269 (13C(6)15N(4), UniMod:267) |

   `silac3` takes two presets, lightest first (`K4R6,K8R10`: M then H); `silac2` takes one
   (`K8R10`). A label not in the table is not guessed: say so, and do not search.
2. **Which sample is in which channel of which run.** Written as `input/silac_design.csv`, one
   row per run and channel that holds a sample:

   ```
   File.Name,Replicate,Channel,Group
   run_01,Rep1,L,WT
   run_01,Rep1,M,KO
   run_01,Rep1,H,DKO
   run_02,Rep2,L,KO          <- a label swap: just another assignment
   run_02,Rep2,H,WT
   run_check,Check1,H,KO     <- a label-check run: the labelled culture injected alone
   ```

   `File.Name` is the run as the report names it (its `Run` column, as in `conditions.csv`).
   `Replicate` defaults to `File.Name`; one replicate is one run. A channel that holds no sample
   has no row. Group names: letters, digits, `_` or `.`, starting with a letter.
   `python3 scripts/silac.py design --validate input/silac_design.csv --against
   ./search_out/report.parquet` checks it (every report run needs rows and every row a run, no
   group twice in a run, at least one run with two groups). A searched run the user decides to
   leave out (a failed injection) is recorded, never just dropped: add `--exclude <run> --reason
   "<why>"` to that command. It goes into `<design>.decisions.json` (`excluded_runs`, the record
   `collect_conditions.py --exclude` writes), which every later check and `run_silac.R` read; the
   Methods and the report then list the run with its reason.
3. **Independent samples or turnover?** The channel-normalisation question DIA-NN cannot
   answer from the data ("This is necessary as the raw data alone does not contain sufficient
   information", README):
   - `--channel-norm spec` (**the default, recorded as one**): independent samples in each
     channel, e.g. different strains or treatments, each grown in its own label, then mixed.
   - `--channel-norm run`: one sample per run whose channels are old and new protein (pulsed
     SILAC, protein turnover). The ratio analysis below compares independent samples, so it
     refuses such a search: a turnover analysis asks another question.

   The ratio analysis reads raw `Precursor.Quantity`, so this setting does not change its
   numbers; it changes DIA-NN's normalised columns (`Precursor.Normalised`, `PG.MaxLFQ`), which
   a client may use.
4. **The comparisons** (step 8): which group over which. The default, every pair with the group
   in the heavier channel over the lighter, is recorded as a DEFAULT.

## 2. The search (steps 6b and 7)

```
python3 scripts/estimate_params.py --engine diann --acquisition DIA --instrument "<detected>" \
    --precursor-mz-range <LO> <HI> --fasta-meta ./search.fasta.meta.json \
    --labels silac3 --label-preset K4R6,K8R10 --channel-norm spec --out ./wf/params.cfg
```

writes, beside the usual library-free cfg (each line tagged in the rationale):

```
--fixed-mod SILAC,0.0,KR,label
--channels SILAC,L,KR,0:0;SILAC,M,KR,4.025107:6.020129;SILAC,H,KR,8.014199:10.008269
--original-mods
--channel-spec-norm
```

and the sidecar's `labels` record (channels, masses, UniMod accessions, the normalisation answer
with `given: true/false`). The labels are the study design, never `--overrides` (refused: as an
override they were recorded as one person's setting). `silac.py flags --labels ...` previews the
lines without writing a cfg.

**Why there is no `--lib-fixed-mod`.** DIA-NN's README example for SILAC (plexDIA, Scenario 1)
is `--fixed-mod SILAC,0.0,KR,label --lib-fixed-mod SILAC`, for a *label-free library* that the
label must be pasted into. Every route of this skill predicts its library from the FASTA, and
`--fixed-mod` already applies the label during that digest ("applying fixed modifications during
FASTA digest", README `--fixed-mod`). Measured on DIA-NN 2.7.0 (HIVE srun 24808559,
2026-10-15, a 6-protein FASTA): `--fixed-mod` alone gives `K(SILAC)`; with `--lib-fixed-mod`,
`K(SILAC)(SILAC)`. A hand-set triple-SILAC search (2026-10-14) carried both: its library read
`K(SILAC)(SILAC)`, the search step applied the label again (`K(SILAC)(SILAC)(SILAC)`), and DIA-NN
logged "contains more than one channel family label, this is not supported". Every later step
reads a library that already carries the label (README Scenario 4: "--lib-fixed-mod must not be
supplied"). The cfg as written here, run through the two-job route on that dataset, logged no
such error, registered the three channels, and its library reads `K(SILAC)` once.

**Which step gets which flag.** The whole declaration is in the cfg, so it reaches every step:

| Route | Library prediction | Search steps |
|---|---|---|
| two jobs (≤ 5 files) | job 1 reads the cfg (`--cfg`); `--channels` is ignored without raw data ("no raw data to analyse") | job 2: the cfg's flags minus the library-only ones |
| one step (`--one-step`) | one DIA-NN command reads the cfg | the same command |
| 5-step chain (> 5 files) | step 1 | steps 1b (probes), 2, 3, 4, 5 |

`run_search.py` (both single-shot shapes) and `diann_parallel.py` refuse a cfg that would
mis-search a multiplexed sample: `--lib-fixed-mod` on a FASTA-predicted library; `--channels`
without a normalisation flag (DIA-NN 2.7.0 stops with "channel normalisation strategy must be
specified" and falls back to `--channel-run-norm`) or without `--original-mods`; a label group no
`--fixed-mod` declares; masses not non-decreasing in channel order; a seed library.
`silac.py check-cfg <cfg>` says the same for any cfg, including one made before 2.13.

**Output.** Two DIA-NN 2.7.0 plexDIA searches wrote no pg/pr/gg matrices although their logs said
they would (`references/search-engines.md`). Everything is read from `report.parquet`, which
has one row per run, precursor and channel: `Channel` (L/M/H; empty for a precursor with no K or
R), `Channel.Q.Value`, a channel-specific `PG.Q.Value`, and `Precursor.Id`, which is the same for
every channel of a precursor (README: "channel-oblivious").

## 3. The ratio analysis (step 8)

```
Rscript scripts/run_silac.R --input ./search_out/report.parquet --design input/silac_design.csv \
    --outdir <S>/output/tables --figures-dir <S>/output/figures --params ./wf/params.cfg \
    --fasta ./search.fasta [--comparisons "KO-WT,DKO-WT"] [--min-replicates 3] \
    [--q-cutoff 0.01] [--pg-q-cutoff 0.05] [--channel-q-cutoff 0.01] [--scatter "KO-WT:DKO-KO"]
```

**The method** (`scripts/silac_ratios.R`; every setting recorded with its source):

1. Precursors at `Q.Value ≤ 0.01` and channel-specific `PG.Q.Value ≤ 0.05`; precursors mapping
   to a contaminant entry (`Cont_`, `contam_`; `contaminants.R`) removed.
2. A precursor ratio log2(A/B) needs both channels at `Channel.Q.Value ≤ 0.01` with a non-zero
   raw `Precursor.Quantity`. DIA-NN's normalised quantities are not used.
3. Per protein group and replicate: the **median** of its precursor log2 ratios.
4. Per replicate: **median-centring** (the uncentred ratios are kept).
5. Protein groups with a ratio in **≥ 3 replicates** (`--min-replicates`): a moderated
   one-sample t-test against 0 (limma `lmFit` with an intercept-only design, `eBayes(robust =
   TRUE)`), Benjamini–Hochberg within each comparison. Significant: adj.P.Val < 0.05, no
   fold-change filter.
6. Protein groups detected (`Channel.Q.Value ≤ 0.01`) in ≥ 3 replicates of some channels and in
   no replicate of another: **channel-specific detections**, listed, not tested.

A label swap needs nothing special: comparisons are between groups, and each replicate's ratio is
taken between the channels its groups are in. Validated on a real triple-SILAC experiment (4
replicates, DIA-NN 2.7.0): the tables of the hand-made analysis reproduced to floating-point
precision (1,224 / 1,153 / 1,093 proteins tested; 816 / 799 / 0 significant; 404/412 and 386/413
up/down; every tested protein's ratio and adjusted p within 1e-14). Its 287 channel-specific
proteins were all found, plus 2 in a pattern the hand-made list did not include.

**Outputs** (in `--outdir`; figures and `figures.json` in `--figures-dir`):

| File | What it holds |
|---|---|
| `DE_silac_<A.B>.csv` | per protein group: `logFC` = mean **median-centred** log2 ratio A/B, `t`, `P.Value`, `adj.P.Val`, `n_replicates`, `Detected_<group>` (k/n replicates), `log2_ratio_uncentred_mean`, `Evidence` (tested, or why not), then per replicate `centred_*`, `uncentred_*`, `n_precursors_*`, and `Description` from the FASTA (UniProt and NCBI headers) |
| `SILAC_channel_specific_detections.csv` | present in some channels, absent from another: the pattern and the replicates per group |
| `SILAC_channel_summary.csv` | per replicate and channel: precursors, protein groups, the share of the signal |
| `SILAC_ratio_offsets.csv` | per comparison and replicate: the median uncentred log2 ratio (the mixing offset) |
| `SILAC_label_incorporation.csv` | label-check runs only |
| `SILAC_sample_matrix.csv`, `SILAC_samples.csv` | per replicate × channel, log2 summed quantity median-centred within the sample: for `sample_quality.py`, not what was tested |
| `de_provenance.json`, `methods.txt`, `sessionInfo.txt`, `reproducibility_log.R` | the record, the Methods paragraph, the R stack, and the statistics as they ran (`silac_ratios.R` verbatim, every value literal; replays the tables exactly) |
| figures | `volcano_<A.B>.png` (one x-axis for all comparisons), `silac_channel_composition.png`, `silac_channel_detections.png`, `silac_scatter_<X>_vs_<Y>.png` (3 comparisons: the first against the last), `qc_pvalue_panel.png` (appendix) |

## 4. Reading the results (say these in the report)

- **Centred, not "amount per volume".** `logFC` is the ratio relative to the typical protein of
  that replicate. Up and down counts come out close to even by construction. The uncentred
  ratios are amount per volume and carry the mixing offset (`SILAC_ratio_offsets.csv`).
- **Unbalanced channels** (`silac_channel_balance`: a channel outside 0.5–1.5 × an equal share):
  equal volume is not equal protein. On the validation dataset the light channel held 68–73% of
  the signal (median uncentred log2 M/L and H/L −1.6 to −1.8, i.e. the light sample put 3–3.5
  times more protein into an equal-volume mix). Centring removes the offset; incomplete labelling
  could not explain it there (incorporation 93–99%, lower bounds, from a DDA search of the same
  samples).
- **Presence / absence.** A ratio needs both channels. A protein lost in the smaller channel has
  no ratio and appears only among the channel-specific detections, so the strongest decreases can
  be missing from the tables. Read that list for whole pathways before calling it a loading or
  detection-limit effect: on the validation dataset a whole regulon, including one of the most
  abundant proteins of the light sample, was there.
- **Shared channels couple comparisons.** B/A and C/A share A's measurement; B/A and C/B share
  B's with opposite sign. Their agreement is not independent evidence (the scatter's caption says
  which channel is shared).
- **Not significant is not unchanged.** The brief carries, per comparison, the centred |log2
  ratio| a typical protein needs to reach p < 0.05 before the adjustment, and to be significant
  after it when it is the only change (from the median posterior variance; a rough reading).
- **Label incorporation** is measured only from label-check runs (a labelled culture injected
  alone: the labelled share H/(H+L) of each precursor; a lower bound, since only precursors with a
  light row are estimated from, and light signal can be noise). In a mixed sample it cannot be
  measured from plexDIA, which gives each peptide one channel; a DDA search with the labels as
  variable modifications can (peptides with two labelable sites: p = 2·full / (2·full + half)).
  Below 95% (`silac_label_incorporation` WARN), unlabelled protein in the labelled cells pulls
  ratios towards the light channel.

## 5. The rest of the pipeline

Nothing else tests for SILAC by name: each step reads what `run_silac.R`'s `de_provenance.json`
says about itself (DE-LIMP rule 1).

| Step | What it reads |
|---|---|
| report (`make_analysis_html.py`) | the DE tables and figures as usual; `glance_tiles` (replicates, channels compared, channel-specific proteins); `report_sections` (channel mixing, channel-specific proteins); Submission section, audit and the Core survey as for any Core run |
| Methods (`make_methods.py`) | the search paragraph's plexDIA sentence and channel rows from the sidecar's `labels` (or the cfg's `--channels`); the `methods_heading` / `methods_paragraph` it wrote |
| audit (`audit_results.py`) | its `checks` (channel balance, replicates, incorporation, a label declared twice, a defaulted normalisation, channel-specific count); `normalization_check` not applicable, quoted; normalisation stability not relevant (`quantities: raw`); a high share significant read with its `high_significance_note` |
| brief (`analysis_prompt.py`) | `brief_notes` and the not-applicable normalisation check |
| delivery (`core_submission.py deliver`) | `normalization_check.status: not_applicable` with a reason passes the gate; one without a reason does not |
| reproducibility (`provenance.py`) | step 5 searches with the recorded cfg (it holds the labels); step 6 is its `replay` (`run_silac.R` with the bundle's paths); the design file is copied into `inputs/` |
| figures, set tests | `make_figures.R` leaves its figures alone (`figures.written_by`); `run_sets.R` stops with `set_tests.reason` |

`run_silac.R` refuses to write into a folder that holds another pipeline's `de_provenance.json`.

## 6. Open questions (for the statistician review)

- **Minimum replicates.** 3 of 4 was the validated choice. With 3 replicates every protein needs
  all three; a protein with 2 is listed untested. Should the minimum scale with n?
- **Precursor aggregation.** The median of precursor log2 ratios weights every precursor equally
  regardless of intensity; an intensity-weighted or MS1-assisted ratio may be more precise.
- **Centring.** The median of all tested proteins' ratios; when most proteins change (a global
  regulator), "typical" may itself move. Centring on a reference set (housekeeping proteins, or
  proteins with the least variance) is an alternative.
- **Presence / absence** has no test or error rate; a detection-limit model (DIA-NN's
  `Channel.Evidence`, or the light channel's intensity distribution) could give one.
- **Power reading** uses the median posterior variance and BH's first step; a simulation from
  the fitted prior would be better.
- **Run-level or global q-values.** The method filters each run's rows on `Q.Value` and the
  channel-specific `PG.Q.Value`; DIA-NN's README also recommends `Global.Q.Value` and
  `Global.PG.Q.Value` ≤ 0.01 across runs. Which suits a within-run ratio is open.
