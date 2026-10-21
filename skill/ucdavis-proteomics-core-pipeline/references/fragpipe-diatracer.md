# Running FragPipe + diaTracer — what actually blocks you

Everything here was found by running the route end-to-end on nine mouse dia-PASEF files
(timsTOF HT, 100 SPD) as a fresh cluster account. Each item cost a failed job, so none
of it is theoretical. **Nothing below is HIVE-specific except where a path is given as
an example** — a laptop or a generic Linux box hits the same four gates.

<!-- contents: tests/test_skill_structure.py --toc -->
**Contents** (section: line)
- 19: The four gates, in the order you hit them
- 124: Parallelising the conversion
- 147: Output contract
- 160: When comparing against the DIA-NN route, say this out loud
- 167: The verified working invocation
- 189: 6c. FragPipe / Radiant: the config is generated, not stored (step 6c)
- 227: 7. Running the FragPipe + diaTracer route (step 7)
<!-- /contents -->

## The four gates, in the order you hit them

### 1. The three engines are licence-gated and cannot be scripted

MSFragger, IonQuant and diaTracer come from a web form plus a licence key (one shared
academic licence; commercial through the tools' commercial licensor). `acquire_tools.sh` can fetch FragPipe
itself but never these. Point FragPipe at wherever you unpacked them.

**The tools folder must be the one that also contains `speclib/`.** A directory holding
only the three jars is not enough — FragPipe looks for `speclib/gen_con_spec_lib.py`
alongside them. On a stock install that is `<fragpipe>/tools`, which already contains all
three jars plus `speclib/`, `Philosopher/` and `diann/`. Handing it a curated
"licensed jars only" folder fails.

> The DIA-NN that FragPipe bundles is explicitly *"available for both academic and
> commercial use within FragPipe"* — DIA-NN is **not** the restriction. Never tell a
> commercial user this route is closed to them; they must licence the three gated tools.

### 2. MSFragger needs a target+decoy FASTA — DIA-NN does not

DIA-NN generates decoys internally. FragPipe does not, and fails immediately:

```
ERROR - No decoys found in the FASTA file.
```

Build one with the Philosopher that ships with FragPipe, using the decoy tag the
workflow expects (`database.decoy-tag=rev_` in the diaPASEF preset):

```bash
PH=<fragpipe>/tools/Philosopher/philosopher-v5.1.3-RC9
cd <somewhere writable>
cp search.fasta target.fasta
"$PH" workspace --init
"$PH" database --custom ./target.fasta --prefix rev_
# -> <date>-decoys-target.fasta.fas   (targets + decoys)
```

Point `database.db-path` at the `.fas`, not the original FASTA. A 55,245-sequence
proteome becomes 110,490 entries.

### 3. Headless FragPipe needs a per-user config cache — THE BIG ONE

On any account that has never opened the GUI:

```
ERROR - Spectral Library Generation module was not configured correctly.
        Please make sure that Python and FragPipe-SpecLib have been installed.
```

**This cannot be fixed with command-line flags.** It survives a correct
`--config-tools-folder`, a correct `--config-diann`, a python with `easypqp` first on
`PATH`, and `--config-python` pointed at either the venv or its binary. Five separate
invocations, five identical failures.

The state FragPipe actually reads is:

```
~/.config/FragPipe/fragpipe/<version>/fragpipe-ui.cache
    fragpipe-config.bin-python=/path/to/python     <-- the key that matters
    fragpipe-config.bin-diann=/path/to/diann
```

Write it without a GUI:

```bash
python3 scripts/fragpipe_bootstrap.py \
    --python <interpreter that can import easypqp> \
    --diann  <fragpipe>/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8
python3 scripts/fragpipe_bootstrap.py --check     # verify, writes nothing
```

`fragpipe_bootstrap.py` verifies the interpreter can actually `import easypqp` before
writing, because a python that can't fails later and less legibly. If you need to build
that interpreter: `python3 -m venv pqp && pqp/bin/pip install easypqp`.

**Every user on a fresh account hits this** — for a class, it has to be part of setup,
not troubleshooting.

`--config-python` is also mis-documented: its help says *"the Python directory"*, but
FragPipe execs the path directly, so a directory gives
`Cannot run program ...: error 13 (Permission denied)`. Pass a binary, or better, set
the cache and omit the flag — the verified working invocation does not use it at all.

### 4. diaTracer standalone does not apply its own defaults

Running the jar directly (which is what parallelising requires) with only `-d`, `-w`
and `-t` dies instantly:

```
Exception in thread "main" java.lang.NullPointerException:
    Cannot invoke "String.trim()" because "in" is null
        at java.base/java.lang.Double.parseDouble
        at diatracer.diaTracerMainClass.main
```

`main()` parses **every** numeric option unconditionally and ignores the defaults its
own `--help` advertises. All seven must be passed:

```
-dI 0.01   -dR 3   -mC 0.3   -mF 1   -mO 0.1   -r 0   -rM 500
```

`diatracer_parallel.py` encodes these; don't "simplify" them away.

## Parallelising the conversion

FragPipe runs diaTracer serially inside its own job — ~20 min per file, so ~3 hours for
nine. The conversion is per-file and independent, so on a cluster it should be an array:

```bash
python3 scripts/diatracer_stage.py --raw *.d --stage ./stage \
    --manifest ./fragpipe.fp-manifest --conditions conditions.csv
python3 scripts/diatracer_parallel.py --stage ./stage --out ./dt \
    --partition low --account <acct> --threads 16 --mem 96
bash ./dt/submit_diatracer.sh
# then RE-STAGE: the manifest becomes the reuse form
python3 scripts/diatracer_stage.py --raw *.d --stage ./stage \
    --manifest ./fragpipe.fp-manifest --conditions conditions.csv
```

**Measured: nine files in ~25 min instead of ~3 h.** The mzML are reusable, so
re-staging emits the two-row reuse form (mzML as `DDA` + the original `.d` as
`DIA-Quant`) and FragPipe skips conversion entirely, going straight to MSFragger.

No cluster? Run the conversion serially and warn the user about the wall-clock; the
mzML still only need building once.

## Output contract

FragPipe 24 bundles **DIA-NN 1.8.2 beta 8**, which has **no parquet writer**. The output
is `<workdir>/dia-quant-output/report.tsv` — there is no `report.parquet` unless someone
configured DIA-NN 2.x. `adapt_fragpipe` handles either, and `compare_searches.py` reads
both. `run_de.R --format tsv` consumes the TSV directly, and the `Run` column matches
the original `.d` basenames, so a `conditions.csv` written for the DIA-NN route works
unchanged.

`combined_protein.tsv` also appears (IonQuant over the pseudo-spectra). **Do not use it
as the quant of record** — it is not equivalent to the DIA-NN numbers and the paper
never uses it.

## When comparing against the DIA-NN route, say this out loud

FragPipe's bundled DIA-NN 1.8.2 beta 8 is **two major versions behind** the 2.x the
`diann_*` workflows pin. Some of any difference is version lag, not the spectrum-centric
approach. To control for it, pass an external DIA-NN 2.x via `--config-diann`. Reporting
a head-to-head without noting this attributes the version gap to the method.

## The verified working invocation

```bash
export PATH=<venv>/bin:$PATH          # python that can import easypqp
"$FP" --headless \
  --workflow  workflow_with_db.workflow \    # database.db-path -> the .fas decoy DB
  --manifest  fragpipe.fp-manifest \         # reuse form after conversion
  --workdir   out \
  --ram 220 --threads 32 \
  --config-tools-folder <fragpipe>/tools \   # must contain speclib/
  --config-diann <fragpipe>/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8
```

Chain: diaTracer → MSFragger → MSBooster → Percolator → ProteinProphet → EasyPQP →
DIA-NN. On nine files with conversion already done: **63.8 minutes**, exit 0.

Never trust SLURM `State` alone — verify `dia-quant-output/report.tsv` exists. A
segfault masked by a trailing `echo` has reported `COMPLETED` with no output.

**Cite:** K. Li et al., *diaTracer enables spectrum-centric analysis of diaPASEF
proteomics data*, Nat Commun 16, 95 (2025).

## 6c. FragPipe / Radiant: the config is generated, not stored (step 6c)

These two are driven by a whole config file rather than flags, so
`resolve_defaults.py` (step 4) already generated one into `./wf/` by calling:

```
python3 scripts/make_presets.py --engine <fragpipe|radiant> \
    --instrument "<detected instrument>" --acquisition DIA \
    [--fasta ./search.fasta] [--threads N] [--cys-alkylation <agent>|none] \
    --out ./wf/<engine>.<workflow|radiantConfig>
```

Re-run it directly if the FASTA or thread count changed after step 4, or with
`--cys-alkylation` when step 3 asked about alkylation: FragPipe's presets fix carbamidomethyl-C
(MSFragger's C row is set to the agent's mass, or 0 for `none`), and Radiant searches whatever
fixed C modification its DIA-NN library carries (`run_search.py` builds that library with it, from
`<out>.provenance.json`).

It starts from the **vendor's own preset for that data type** and patches only what
is run-specific. That matters: FragPipe's two DIA presets differ in **20 keys**, not
one — MBR, missed cleavages, precursor charge and tolerance, topN peaks, and the
var-mod table — so flipping `diatracer.run-diatracer` on the wrong template yields a
config no vendor ever validated. Never hand-edit a preset to switch routes; pick the
right template by data type.

- **timsTOF** → `DIA_SpecLib_Quant_diaPASEF.workflow`, diaTracer **on**.
- **Thermo** → `DIA_SpecLib_Quant.workflow`, diaTracer **off** (MSFragger-DIA
  searches the spectra directly; there are no pseudo-spectra to trace).
- **Radiant on Bruker `.d` is refused**, not adapted — the container reads only
  mzML/Parquet. It also narrows Seer's flat 20/20 ppm default using the instrument,
  tagged `[DERIVED]` because that inference is ours, not Seer's. A level with no documented
  DIA-NN value (e.g. a 15k MS2) keeps the vendor 20 ppm, and the provenance says which.
- Mass tolerances otherwise keep the **vendor's** values: FragPipe already tunes them
  per data type, and DIA-NN's ppm table describes DIA-NN's matcher. Override only
  with a value the user gives (`--ms1-ppm/--ms2-ppm`, with `--override-by`/`--override-reason`).

`<out>.provenance.json` records every key changed and why — it feeds the methods text.

## 7. Running the FragPipe + diaTracer route (step 7)

- **FragPipe + diaTracer (timsTOF dia-PASEF, spectrum-centric DIA).** A second DIA route,
  used when the bundle names it or the user asks for it. diaTracer traces 3-D features and
  writes pseudo-MS/MS spectra, MSFragger searches those to build a library from the data
  itself, then **DIA-NN quantifies against it**. Worth reaching for as a cross-check on the
  DIA-NN route, or when the search needs semi-tryptic / nonspecific / open-modification
  behaviour that library-free DIA-NN handles poorly.
  - **⚠ Licensing — get this right, it is easy to state wrongly.** The gate is
    **MSFragger / IonQuant / diaTracer** (one shared academic licence; commercial users
    license through the tools' commercial licensor). The DIA-NN that FragPipe bundles is explicitly *"available
    for both academic and commercial use within FragPipe"*, so DIA-NN is **not** the
    restriction here. A commercial user is not locked out — they must license the three
    gated tools. **Never tell them the route is closed to them.**
  - **⚠ The clause that matters for a core facility:** the academic licence forbids *"any
    work product, report or deliverable by LICENSEE, for a commercial or for-profit
    organization"* — which reads on fee-for-service work for an industry client. If asked
    whether their work qualifies, send them to the tools' commercial licensor and their licensing
    office. Do **not** rule on it yourself. → `references/search-engines.md`.
  - Output is DIA-NN's own, in `<workdir>/dia-quant-output/`. With FragPipe's bundled
    DIA-NN 1.8.2 Beta 8 that is **`report.tsv` — there is no `report.parquet`** unless
    someone configured DIA-NN 2.x. `adapt_fragpipe` handles either.
  - **Never point FragPipe at shared raw files directly — `run_search.py` stages symlinks
    for you.** diaTracer writes its mzML *next to the input*, so two people searching the
    same dataset (a class working from one folder) would race to write the same file, and a
    read-only share fails outright. `diatracer_stage.py` builds a per-run directory of
    symlinks; FragPipe normalizes paths without resolving them, so the output lands in
    **your** directory and the shared raw folder is never written to. This happens
    automatically for DIA + `.d` inputs; nothing to pass.
  - **Existing pseudo-spectra are reused, never regenerated.** The stager looks for the
    mzML in the staging dir *and* beside the real `.d`, and when it finds one emits the
    reuse form — the mzML as `DDA` plus the `.d` as `DIA-Quant` — so DIA-NN still
    quantifies against the real chromatograms. That saves ~20 min per file, so a re-search
    with different parameters is cheap. Report how many were reused.
  - A manifest data type that isn't one of `DIA`/`GPF-DIA`/`DIA-Quant`/`DIA-Lib`/`DDA+`/
    `DDA` is **silently treated as DDA**, quietly turning a DIA run into a DDA one.
  - FragPipe 24.0 bundles **DIA-NN 1.8.2_beta_8**, older than the 2.x the `diann_*`
    workflows pin. Expect different numbers from that alone; say so rather than
    attributing every difference to the spectrum-centric approach.
  - **Four things block this route on a fresh account — read "The four gates" above
    BEFORE running it.** (1) the three engines are licence-gated and the tools folder must be the
    one that also contains `speclib/`; (2) MSFragger needs a **target+decoy** FASTA, which DIA-NN
    does not; (3) headless FragPipe rejects the spectral-library module until
    `fragpipe-config.bin-python` exists in `~/.config/FragPipe/.../fragpipe-ui.cache` — no CLI flag
    substitutes, run `scripts/fragpipe_bootstrap.py`; (4) diaTracer's standalone `main()` ignores its
    own documented defaults and needs all seven numeric options passed.
  - **On a cluster, parallelise the conversion automatically** with `scripts/diatracer_parallel.py`
    — one SLURM task per file instead of FragPipe's serial loop. Measured 9 files in ~25 min vs ~3 h.
    Re-stage afterwards so the manifest takes the reuse form and FragPipe skips conversion.
  - Budget roughly **20 min per file for the diaTracer conversion alone** (paper: 34 files,
    32 cores, ~19.4 min each), before MSFragger and DIA-NN. The mzML are **reusable**, so a
    re-search with different parameters skips it. Needs Java 11+, MSFragger 4.4+,
    IonQuant 1.11.18+, diaTracer 2.2.1+. The three gated jars are **not** in the FragPipe
    zip and **cannot be scripted** — a human accepts the licence once and mints a key, then
    point `FRAGPIPE_TOOLS_FOLDER` at them. → `references/search-engines.md`.
