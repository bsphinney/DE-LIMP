# Study questions (SKILL.md step 3)

What to ask the user before the database is built and before the one compute confirmation,
and what each answer changes. Every one of these is asked, never inferred from file or folder
names, a sample ID or an earlier session. Each answer is recorded once the session exists
(step 3b: `experiment_type.py set`, `sample_type.py set`, `log_decision.py`).

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 21: Golden rule 4: the organism and the sample type are asked
- 42: Organism
- 63: Database type
- 85: Contaminants
- 98: SILAC (plexDIA)
- 112: Experiment type
- 141: Reduced and alkylated? (enrichment experiments only)
- 155: Sample type
- 177: Keratin samples
<!-- /contents -->

## Golden rule 4: the organism and the sample type are asked

4. **Organism decides the FASTA and nothing else**, and is **always asked, never
   assumed** — not from a folder name, not "probably human". The **sample type** (tissue,
   cell line, biofluid, IP …) is likewise **asked, never inferred** — not from a file-name
   token or abbreviation, a sample ID or a folder — though Core staff (any analysis: a
   submission or internal work) are asked it only when a step depends on it (step 3,
   `sample_type.py`). Resolve the organism with
   `fetch_fasta.py resolve` and have the user confirm the proteome. It does **not**
   affect a single search parameter: mass accuracy keys on instrument, and search
   behaviour on acquisition. Both are auto-detected and confirmed. So any organism
   is supported out of the box — there is no such thing as an unsupported species.


**Every organism works.** There is no "unsupported species" case any more —
parameters key on instrument and acquisition, and the FASTA comes from step 3. Never
tell a user their organism isn't covered.

⚠ **The organism always comes from step 3**, never from anything cached or inferred.
Searching rat data against a human proteome silently voids the whole run.

## Organism

- **Organism** cannot be detected — you MUST ask the user, every run. **Never** take
  it from a previous session, the folder name, or an assumption that it's human.
  Nothing in the skill supplies a default organism, by design.
  Getting this wrong silently searches the wrong species and every
  downstream number is void. Resolve their answer — don't guess the proteome ID:

```
python3 scripts/fetch_fasta.py resolve --organism "<what they said>"   # or --taxid <n>
```

  Accepts a common name, a scientific name, an NCBI taxid, or a `UP…` accession.
  (`scripts/resolve_organism.py --organism …` is a thin alias for the same resolver,
  kept for older call sites — there is only one implementation.)
  Show the user `selected` (organism, proteome ID, protein count) and **confirm it**.
  If `needs_menu` is true, present `candidates` as a menu instead of auto-picking —
  it fires when the top hit isn't a clear reference proteome or several compete
  (e.g. "mouse" also matches *Myotis myotis* and mouse-ear cress). Carry the
  confirmed `taxid` into step 4 and the `proteome_id` into step 6.

## Database type

- **Database type** — ask what kind of protein database they want. Don't just
  take a one-word answer if they sound unsure: **most users don't know, and that's
  expected.** Explain the trade-off in their terms and recommend — every extra
  sequence enlarges the search space, which costs sensitivity through FDR
  correction, so bigger is not better:
  - `one_per_gene` — **the default; recommend this.** One canonical protein per
    gene (human ≈ 20.6k). Right for essentially all standard expression /
    differential-abundance work, which is what almost everyone is doing.
  - `reviewed` — Swiss-Prot only (manually curated). Use for organisms where the
    reference proteome is padded with unreviewed predictions.
  - `full` — every entry incl. unreviewed TrEMBL (human ≈ 147.5k, **7× bigger**).
    Only for poorly-annotated organisms where curation is thin.
  - `*_isoforms` variants — add splice isoforms. Only if the biological question
    is *about* isoforms. Isoforms share most peptides, so they mostly add
    ambiguous protein groups rather than answers.

  If they don't know, ask what they're trying to find out (standard "which
  proteins changed between my groups?" → `one_per_gene`), then say which you're
  using and why. Don't make them learn UniProt vocabulary to answer.

## Contaminants

- **Contaminants** — ask which set to append (they are all `Cont_`-tagged, so
  DIA-NN's `--cont-quant-exclude Cont_` keeps them out of its normalisation, and
  `run_de.R` removes them before the DE either way):
  - `universal` — **default**, and what the UC Davis Core stages on HIVE. Use it
    unless the user has a reason not to.
  - sample-type matched: `cell_culture`, `mouse_tissue`, `rat_tissue`,
    `neuron_culture`, `stem_cell_culture` — offer the one matching the sample type the user
    gave (below), never one matched to a file name.
  - `none` — only if the user explicitly declines. Record it; the contaminant
    anomaly check in step 10 is meaningless without contaminants in the database.

## SILAC (plexDIA)

- **SILAC (plexDIA) — only when the submission or the user says SILAC / labelled.** Ask four
  things, never guess them: (1) the Lys/Arg labels of each labelled channel (`silac.py presets`:
  K4R6, K6R6, K8R10); (2) which sample is in which channel of which run (label swaps are fine),
  written as `input/silac_design.csv` (`File.Name,Replicate,Channel,Group`; a labelled-only
  label-check injection is a run with that channel alone); (3) are the channels independent
  samples (strains, treatments: `--channel-norm spec`, the default) or old and new protein of
  one culture (pulsed SILAC / turnover: `run`, which the ratio analysis does not cover); (4) is
  each replicate a separate culture and mix, or are some runs re-injections of one mix (those
  share a `Replicate` and are averaged, never counted twice; confirmed separate →
  `--replicates biological` in step 8). Record the answers (`log_decision.py`). Detail:
  `references/silac.md`.

## Experiment type

- **Experiment type — decided once, confirmed with the user, never inferred from file
  names.** It decides the DE's quantities (step 8) and the pull-down analysis (8a, the brief).
  With a CoreOmics submission attached (step 3b), `python3 scripts/experiment_type.py propose
  --session <S>` reads its experiment type and proposes one; **confirm it with the user**.
  Without one, ask: whole proteome / IP · AP-MS · pull-down with IgG or bead controls /
  proximity labelling (TurboID · BioID · APEX) / secretome / fractionated for depth (offline
  high-pH, fractions combined per sample) / separation or profiling fractions compared with each
  other (SEC, density or sucrose gradient, BN-PAGE complexome, organelle or subcellular) / other.
  A bare "fractions" is not enough: ask which. Each has a default for the DE's quantities
  (`experiment_type.DEFAULTS`):
  - **IP / AP-MS / pull-down:** **non-normalised** — the controls carry little protein and are
    not loaded by amount, so normalisation would scale them up and hide the enrichment;
  - **proximity labelling:** normalised — the controls carry plenty of signal and the samples
    are usually loaded by equal protein;
  - **separation / profiling fractions:** **non-normalised** — each fraction holds a different
    part of the proteome; DIA-NN's README: for fractions from a separation technique such as SEC,
    "normalisation should not be used";
  - **whole proteome, fractionated for depth:** normalised; **secretome:** normalised, leaning on
    the data check.
  The default is applied only after step 8's data check, which runs for EVERY type and stops
  to ask when the data disagree. Record the answer once the session exists (step 3b), with the
  bait (the IP target or the TurboID fusion; asked ONCE — when its sequence will be added with
  `--add-fasta`, ask it after step 6 instead: SKILL.md step 6) and, for controls not named IgG/beads,
  `--controls`.
  When the submission's "Normalization" answer (how the samples were loaded) contradicts the
  default, `set` prints `ASK:` — ask before going on.

## Reduced and alkylated? (enrichment experiments only)

- **Reduced and alkylated? Ask ONLY for an enrichment experiment** — IP / AP-MS / pull-down or
  proximity labelling (`experiment_type.ENRICHMENT`), or a submission that says on-bead
  digestion. On-bead digests are often not reduced and alkylated, and the default fixed
  carbamidomethyl-C then leaves no cysteine peptide searchable (0.03–0.5% of identified peptides
  held a C on three such Core searches, against ~8%). **If yes, ask which agent.** Pass the
  answer to `estimate_params.py` (step 6b; `make_presets.py` for FragPipe / Radiant, step 6c) and
  log it (`log_decision.py`) as `--cys-alkylation`: iodoacetamide or chloroacetamide → `cam`;
  acrylamide → `propionamide`; N-ethylmaleimide → `nem`; iodoacetic acid → `carboxymethyl`;
  MMTS → `mmts`; 4-vinylpyridine → `pyridylethyl`; **not alkylated → `none`** — each recorded as
  the user's answer. For every other experiment type do not ask and **omit the flag**: fixed
  carbamidomethyl stands, recorded in the sidecar as a DEFAULT.

## Sample type

- **Sample type — never infer it; ask outside users, and Core staff only when a choice
  depends on it.** **Never** take it from file or folder names, sample IDs or abbreviations: a
  file-name abbreviation was once written up as a tissue the data did not support.
  - **Outside users: ask**, with the other step-3 questions, what the samples are (cell line,
    tissue — which one —, plasma, IP, hair …), in their words; "don't know" is an answer
    (`--not-stated`, and the report then names no tissue). No submission record describes their
    samples, and they are making the choices below without the Core's defaults.
  - **Core staff, any analysis (a CoreOmics submission or internal work): do not ask it.**
    Record `--not-asked` (step 3b) — or `--stated` with what they already said ("HeLa digest
    standard") — and go on; a submission record, when there is one, describes the samples for
    the report (step 9, `--submission`). Ask only at a step whose choice depends on the answer:
    the keratin question below (one yes/no, before the FASTA is built), and a flagged
    sample-quality panel (8c). The contaminant set is the Core's `universal`, whatever the
    tissue.
  The answer is used at every step that depends on it: the contaminant set below, the keratin
  decision below, the sample-quality panels (8c: a muscle panel is the analyte in a muscle
  sample; hemolysis applies to plasma/serum), the analysis brief and report (9), and the SDRF's
  organism part (`make_deposit.py`, still filled by the user). Record it once the session
  exists (step 3b, `sample_type.py set`).

## Keratin samples

- **Is the sample keratin?** Decide it from the sample type the user gave (above). Core staff,
  not asked the sample type: ask only "are the samples keratin tissue (hair, wool, skin, nail,
  feather …)?" when the submission's description or anything they said names one of those
  tissues; otherwise pass neither `--keratin-sample` nor `--no-keratin-sample` (step 6): the
  sidecar then records the default, never an answer nobody gave. When the
  samples are **hair, wool, fur, feather, skin (epidermis /
  stratum corneum), nail, claw, hoof, horn, beak, baleen, quill or scale**, this is a
  **keratin sample**: keratin is the *analyte*, not handling contamination. Decide it here,
  once, and carry it through with `--keratin-sample`: `fetch_fasta.py fetch` (step 6) then
  removes every keratin-family `Cont_` entry and records `keratin_sample: true` in the
  sidecar; `run_search.py` (step 7) refuses a database that still holds them; `run_de.R`
  (step 8) and the auditors (step 8c) read it from the sidecar / `search_provenance.json`.
  **Deciding it after the search is too late**: the contaminant set holds ~190 keratins and
  KRTAPs (human, mouse hair, sheep wool), and every peptide the sample's keratins share with
  one is taken out by DIA-NN's `--cont-quant-exclude Cont_` and by `run_de.R`'s contaminant
  filter (the Core's hair benchmark SET28 lost them this way). Tell the user the trade-off:
  keratin from handling can then no longer be told apart from the sample's own. Three limits of
  this release, worth saying to the user when they apply:
  - **Non-human hair, wool, feather:** human skin keratins (KRT1/2/9/10) are genuine handling
    contamination there, and dropping every keratin entry loses that readout — watch KRT1/2/9/10
    in the results rather than assuming they are the sample's.
  - **Skin:** filaggrin, filaggrin-2, hornerin and dermokine (FLG/FLG2/HRNR/DMKN) are analytes
    too, but only keratins are treated so; for human they are already quantified (the identity
    rule drops those entries), for other species they stay contaminants.
  - **Poorly annotated species:** when the target proteome lacks KRTAPs, the contaminant set's
    KRTAPs may be the only KRTAP sequences in the search, and `--keratin-sample` removes them —
    those peptides then go unidentified. Say so, or add a related species' KRTAPs to the database.
  The list of
  tissues is `fetch_fasta.KERATIN_SAMPLE_TISSUES`; for anything else (a keratinocyte culture,
  a tissue with some skin) ask, and treat it as a keratin sample only if keratin is what they
  want to measure.
