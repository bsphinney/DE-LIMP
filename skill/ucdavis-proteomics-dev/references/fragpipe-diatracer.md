# Running FragPipe + diaTracer — what actually blocks you

Everything here was found by running the route end-to-end on nine mouse dia-PASEF files
(timsTOF HT, 100 SPD) as a fresh cluster account. Each item cost a failed job, so none
of it is theoretical. **Nothing below is HIVE-specific except where a path is given as
an example** — a laptop or a generic Linux box hits the same four gates.

## The four gates, in the order you hit them

### 1. The three engines are licence-gated and cannot be scripted

MSFragger, IonQuant and diaTracer come from a web form plus a licence key (one shared
academic licence; commercial via fragmatics.com). `acquire_tools.sh` can fetch FragPipe
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
