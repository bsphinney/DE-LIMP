#!/usr/bin/env python3
"""
diann_parallel.py  --  Generate DIA-NN's canonical 5-step PARALLEL search as a SLURM
job chain. This is the high-throughput DIA-NN workflow that is poorly documented
upstream; ported faithfully from DE-LIMP's generate_parallel_scripts() (R/helpers_
search.R) and the facility's usage.

The 5 steps (each chained `afterok` on the previous):
  1  library prediction  single job, no raw — predict a spectral library from the FASTA
  2  first pass          SLURM array (1 file/task) — search vs the predicted lib -> .quant
  3  empirical assembly  single job — `--use-quant` over step-2 .quant -> empirical lib
  4  final pass          SLURM array — search vs the empirical lib -> .quant
  5  cross-run report    single job — `--use-quant --matrices` -> report.parquet

Why it's faster: the per-file passes (steps 2 & 4) run as a SLURM **array** across many
nodes at once instead of one long single-node job; MBR is replaced by the empirical-
library round-trip. **Mass accuracy is FIXED, not auto** — steps 3/5 reuse the .quant
files and auto-calibration would be inconsistent (per DIA-NN dev guidance). It is fixed in
the cfg, or -- for an Orbitrap with no documented DIA-NN value -- measured once in step 1b
and pinned for steps 2-5 (see parallel_safe).

Writes into <out>: `file_list.txt`, `step{1..5}_*.sbatch`, and `submit.sh` (submits the
chain with dependencies). Run `submit.sh` on the cluster (or via `hive_exec.sh`). All
heavy work runs on compute nodes through the array — never the login node.

Usage:
  python3 diann_parallel.py --diann '<diann binary | apptainer exec ... diann-linux>' \
      --raw /data/*.d --fasta /path/search.fasta --out ./diann_parallel \
      --cfg params.cfg [--threads-max 16 | --threads-per-file N] [--mem-per-file 32] \
      [--time-per-file 2] \
      [--assembly-cpus 64] [--assembly-mem 128] [--assembly-time 12] \
      [--partition <auto>] [--account <auto>] [--max-simultaneous 20] [--no-norm]
"""
import os, re, sys, glob, argparse, json, shlex, subprocess, math, stat, time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# ONE definition of what a plausible Orbitrap mass accuracy is: probe_window.py measures it and
# refuses to pin a value outside the band, and needs_measured() below refuses to let one that
# reached massacc.txt some other way onto a DIA-NN command line.
from probe_window import MASS_ACC_BAND, SOP_MASS_ACC, band_text        # noqa: E402
# what a probe's exit status allows: the ONE rule (probe_window: fall back only when the probe's
# own machinery failed -- never on a refusal, the environment, the arguments or the data)
from probe_window import RETRY_ON, FALLBACK_ON, EXIT_MEANING            # noqa: E402
# ONE definition of "this cfg searches DDA" (the flag estimate_params.py writes for DDA), and of
# what an unset --window means under it.
from estimate_params import DIANN_DDA_FLAG, DDA_WINDOW_NOTE, is_dda     # noqa: E402

# The same SOP floor, keyed by the DIA-NN flag rather than by level, for the provenance records.
SOP_MASS_ACC_FLAGS = {"--mass-acc": SOP_MASS_ACC["ms2_ppm"],
                      "--mass-acc-ms1": SOP_MASS_ACC["ms1_ppm"]}


def floor_note(evidence_file, value_file):
    """Why `measured: true` does not mean "the search ran at the measured value".

    A measured level is floored at the SOP, so the pinned tolerance is max(measured, SOP). The
    two can differ and the record must not let a reader assume they do not."""
    return ("`measured: true` means the tolerance IS measured with DIA-NN, not that the search "
            "runs at the measured number: a MEASURED level is floored at the SOP "
            + ", ".join(f"{f} {v:g}" for f, v in sorted(SOP_MASS_ACC_FLAGS.items()))
            + f", so the pinned value is max(measured, SOP) per level. {evidence_file} records "
              "the raw median as mass_acc.measured_ms2_ppm / .measured_ms1_ppm, what was pinned "
              "as mass_acc.pinned_ms2_ppm / .pinned_ms1_ppm, and whether the floor applied as "
              f"mass_acc.floored. {value_file} holds the pinned pair, as passed to DIA-NN. A "
              "level given from DIA-NN's resolution table is pinned as given and never floored.")

# flags that are step-specific or auto-determined — never carry them into every step.
# NOTE: --dda is intentionally NOT stripped: estimate_params.py writes it into a DDA cfg and it
# flows into every step, library prediction included (DIA-NN 2.6 searches DDA per file exactly
# as it does DIA). run_search.py refuses a cfg whose --dda disagrees with the bundle's
# acquisition (dda_mismatch) before it generates anything.
STRIP = ("--fasta-search", "--predictor", "--gen-spec-lib", "--matrices", "--reanalyse",
         "--rt-profiling", "--no-norm", "--xic", "--mobilograms", "--out-lib", "--lib", "--out", "--f",
         # NOTE: --xic is stripped here on purpose and re-added to step 4 ONLY (see
         # xic_flag() below) -- step 2 IDs are not final, and step 5 runs --use-quant,
         # which never re-reads the raw spectra so --xic is silently a no-op there.
         "--fasta", "--threads", "--temp")

# Step 3's cross-run report: the FIRST pass (each run searched once against the predicted
# library), which step 5 compares its own report with (pass_comparison.py). The one name of it:
# the comparison, the Methods (make_methods.py, when it is the deliverable) and SKILL.md use it.
FIRST_PASS_REPORT = "step3_assembly.parquet"
# Step 5's report under --no-norm: DIA-NN's quantities with no cross-run normalisation. Its name
# is the record of that (normalization_check.searched_no_norm reads it, so run_de.R's MaxLFQ
# says what ran without re-reading the report).
NO_NORM_REPORT = "no_norm_report.parquet"

# The per-task --out of the array steps. Without one, DIA-NN writes report.parquet,
# report.stats.tsv, report-lib.parquet and report.log.txt into the WORKING directory -- <out>,
# from every task, concurrently, on several nodes: a one-run report at the very path run_de.R is
# pointed at if steps 3-5 never finish, and a stray report-lib.parquet that record_run.py lists
# as a library (SET28, HIVE 2026-09-29). One folder per step, one file set per task.
TASK_OUT_DIRS = {"step2": "firstpass", "step4": "finalpass"}

# Step 1b measures the radius on this many REPRESENTATIVE runs (the median and quartile runs
# of the cohort, never a blank, wash, failed injection or a .d with a damaged index -- see
# probe_window.py) and pins the median. DIA-NN's README: auto-optimised values "depend on which
# run is first in the list"; one timsTOF cohort gave 10, 11 or 14 depending on the run probed.
PROBE_CANDIDATES = 3
# A run that logs no radius is replaced by the next run nearest the median, so one bad run does
# not fail the cohort. After this many such runs step 1b gives up (the radius is a property of
# the method, so repeated failures mean something is wrong beyond one run).
PROBE_MAX_FAILURES = 3
PROBE_TIMEOUT_S = 3600      # per probe: probe_window.py's own default, and step 1b's effective
                            # limit before it tried more than one file. Nothing shorter has been
                            # measured on a large Astral .raw, so it is not shortened.
PROBE_WALL_HOURS = -(-PROBE_CANDIDATES * PROBE_TIMEOUT_S // 3600) + 1   # 3 full probes + 1 h
# ONE budget for all probes, replacements included (they can outnumber what the wall clock covers
# at the full per-probe timeout): the wall clock less 10 minutes, so the probe stops itself and
# writes window.json before SLURM kills the job with no evidence.
PROBE_BUDGET_S = PROBE_WALL_HOURS * 3600 - 600
# A probe that measured nothing is retried ONCE in the same job, from what is left of that one
# budget -- and only when at least this much is left: a retry that cannot finish one probe only
# delays the fallback (probe_attempts()).
PROBE_RETRY_MIN_S = 900
# ZERO is `invalid` for both flags, never `auto`: `auto` means the flag was not set, and DIA-NN
# does not read a literal 0 on the command line as "optimise automatically" for either (its
# README's "set to 0 ... optimise them automatically" describes the GUI fields, which omit the
# flag at 0 -- https://github.com/vdemichev/DiaNN, "Changing default settings"):
#   --window 0     DIA-NN logs "scan window radius should be a positive integer" and then
#                  chooses a radius per file (the 18-file poplar run: 7 for seventeen, 8 for
#                  one; _window_value)
#   --mass-acc 0   a literal 0 ppm tolerance: "Mass accuracy will be fixed to 0 (MS2) and 0
#                  (MS1)", 0 IDs at 1% FDR on a 28-run Lumos search (_mass_acc_value, PR #38)
# search_provenance.json `scan_window.mode` and `mass_acc.mode` (also on result.scan_window /
# result.mass_acc): STABLE, machine-readable values -- FRAN ingests them, from fran_manifest.json's
# copy of the provenance, as a variable of its DIA-NN vs Spectronaut comparison. Never rename or
# reuse one; a new case gets a new value, added here, in references/environment.md
# ("search_provenance.json: scan_window.mode / mass_acc.mode") and in
# tests/test_probe_estale.py's ProvenanceModeTests, which pin them.
SCAN_WINDOW_MODES = {
    "measured": "measured on these runs by step 1b and pinned for every step (set at generation, "
                "as the plan; replaced by fallback_auto if the measurement fails)",
    "fallback_auto": "step 1b's measurement failed; DIA-NN chose the radius itself, per run",
    "pinned": "given in the cfg and passed to every step",
    "auto": "not set, by design (single-shot search, DDA): DIA-NN chose the radius itself",
    "invalid": "the cfg passes a --window that is not one positive integer, 0 included (see "
               "below); what DIA-NN does with anything but 0 is unverified",
    "unknown": "the cfg could not be read"}
MASS_ACC_MODES = {
    "measured": "measured on these runs before the search: a measured level is floored at the "
                "SOP, a documented level passed as documented (set at generation, as the plan; "
                "replaced by fallback_default if the measurement fails)",
    "fallback_default": "the measurement failed; the documented level as given, the other at the "
                        "facility SOP -- DEFAULT, not measured",
    "pinned": "given in the cfg",
    "pinned_default": "pinned by estimate_params.py at the facility SOP for a level that cannot "
                      "be measured (DDA) -- DEFAULT (`default` names the levels)",
    "auto": "neither level set: DIA-NN optimised it itself",
    "partial": "one level given in the cfg: DIA-NN 2.7.0 then fixes both, the other at 20 ppm "
               "(estimate_params.LONE_FLAG_NOTE)",
    "invalid": "the cfg passes a value DIA-NN does not read as a tolerance: 0 (see below), "
               "negative, not a number, or set twice with different values",
    "unknown": "the cfg could not be read"}

# probe_fallback.py's record, beside window.txt / massacc.txt: steps 2-5 accept window.txt `auto`
# only with it (needs_measured)
FALLBACK_RECORD = "probe_fallback.json"
# What a pre-search probe leaves in <out> -- step 1b's, or the single-shot search's
PROBE_OUTPUTS = ("window.txt", "massacc.txt", "window.json", "window.json.attempt1",
                 "mass_acc.json", "mass_acc.json.attempt1", FALLBACK_RECORD)


def _not_a_regular_file(path):
    """None if `path` is absent or a regular file; otherwise what it is (for the refusal)."""
    try:
        mode = os.lstat(path).st_mode
    except FileNotFoundError:
        return None
    if stat.S_ISREG(mode):
        return None
    return ("a directory" if stat.S_ISDIR(mode) else "a symlink" if stat.S_ISLNK(mode)
            else "not a regular file")


def set_aside(path):
    """Rename an existing REGULAR FILE out of the way -- never delete it. Returns the new path,
    or None when there was nothing there. The name records that it is stale and when it was
    moved. Anything else (a directory, a symlink, a device) raises ValueError: `--sbatch proj`
    once renamed a whole project folder, including the cfg the search was about to read. The one
    definition: run_search.py's re-run and --sbatch paths use it too."""
    kind = _not_a_regular_file(path)
    if kind:
        raise ValueError(f"{path} is {kind}, not a job script -- refusing to move it")
    if not os.path.lexists(path):
        return None
    base = f"{path}.stale-{time.strftime('%Y%m%dT%H%M%S')}"
    new, k = base, 1
    while os.path.lexists(new):
        new, k = f"{base}.{k}", k + 1
    os.rename(path, new)
    return new


def set_aside_probe_outputs(out):
    """Set an earlier search's probe outputs in `out` aside (set_aside: renamed `.stale-<time>`,
    never deleted) when the next search is GENERATED there, and say so on stderr. Each job
    removes them before it probes, but a search that does not probe (a pinned --window, a pinned
    mass accuracy) never did: dda-review N1 (2026-09-30) generated a pinned-window chain into a
    folder holding an earlier chain's probe_fallback.json, and checkpoint.py status, the Slack
    post and watch_run.sh --all then reported THIS search as fallen back. Callers run this only
    once every refusal has passed (dda-review R1: when it ran first, a REFUSED generation deleted
    a completed search's window.txt / massacc.txt / window.json and left its provenance pointing
    at files that were gone). Returns [(name, new path)]."""
    moved = []
    for name in PROBE_OUTPUTS:
        path = os.path.join(out, name)
        try:
            new = set_aside(path)
        except ValueError as e:                  # never a reason to stop a generated search:
            sys.stderr.write(f"[probe outputs] WARNING: {e}; left in place\n")   # it is said
            continue
        if new:
            moved.append((name, new))
    if moved:
        sys.stderr.write("[probe outputs] an earlier search's " + ", ".join(n for n, _ in moved)
                         + f" in {out} do not describe this search: set aside as "
                         + ", ".join(os.path.basename(p) for _, p in moved) + "\n")
    return moved


# What the search runs with when it did not measure (probe_fallback.py writes it), as the plan
# records it at generation.
FALLBACK_PLAN = {
    "window": "retried once; if it still measures nothing, window.txt says `auto`, steps 2-5 pass "
              "no --window and DIA-NN chooses the radius itself, per run -- recorded as a "
              "fallback (search_provenance.json probe_fallback), never as measured",
    "mass-acc": "retried once; if it still measures nothing, the documented level as given and "
                "the other at the facility SOP, tagged DEFAULT -- recorded as a fallback "
                "(search_provenance.json probe_fallback), never as measured"}


def dotnet_prefix(raws):
    """If any input is Thermo .raw, the DIA-NN 2.6 native binary needs a .NET 8 runtime
    (>= 8.0.17) on PATH to read it. Resolve/install via ensure_dotnet8.sh and return an
    'export DOTNET_ROOT=...; export PATH=...; ' prefix to put in front of every DIA-NN
    invocation (each array task/step runs it; the shared install is read on the node).
    Returns "" for mzML/.d-only inputs. Run the generator on a login node (internet)."""
    if not any(r.lower().endswith(".raw") for r in raws):
        return ""
    helper = os.path.join(os.path.dirname(os.path.abspath(__file__)), "ensure_dotnet8.sh")
    try:
        root = subprocess.check_output(["bash", helper], text=True).strip().splitlines()[-1]
    except Exception as e:
        sys.stderr.write(f"[diann_parallel] ensure_dotnet8.sh failed ({e}); DIA-NN may not "
                         "read .raw. Provide mzML or install .NET 8 >= 8.0.17.\n")
        return ""
    return f'export DOTNET_ROOT={root}; export PATH={root}:"$PATH"; '


class CfgError(ValueError):
    """The cfg cannot be read as flags. `code` is "cfg_missing" (no such regular file) or
    "cfg_unparseable" (e.g. an unclosed quote) -- parallel_safe() reports them separately, so a
    missing file never reads as "mass accuracy is not pinned"."""

    def __init__(self, msg, code="cfg_unparseable"):
        super().__init__(msg)
        self.code = code


def _split_cfg_text(text, where):
    """Split cfg text into words by BASH's quoting and comment rules, with no expansion.

    shlex (comments=True) was the rule before, and it is not bash: it ends a word at ANY `#`,
    so `/data/run#1` read as `/data/run` while bash -- which the chain splices these words into
    -- keeps the `#` and starts a comment only at the START of a word. The gate and the steps
    must read one file the same way, so this follows bash:
      * space, tab and newline separate words (so does CR, which bash would keep -- a cfg saved
        with Windows line endings must not turn `--window 7` into `7\\r`)
      * '...' is literal; "..." is literal except \\$ \\` \\" \\\\ and \\<newline>
      * outside quotes a backslash escapes the next character; backslash-newline joins lines
      * `#` starts a comment only where a word would start
    Every other character is literal -- `(`, `{`, `~`, `*`, `$` included. Expansion is decided
    on the way OUT, by _shield(), not here."""
    toks, cur, in_word, i, n = [], [], False, 0, len(text)
    while i < n:
        c = text[i]
        if c in " \t\r\n":
            if in_word:
                toks.append("".join(cur))
                cur, in_word = [], False
            i += 1
        elif c == "#" and not in_word:
            j = text.find("\n", i)
            i = n if j < 0 else j
        elif c == "'":
            j = text.find("'", i + 1)
            if j < 0:
                raise CfgError(f"{where} cannot be parsed as DIA-NN flags (unclosed ')")
            cur.append(text[i + 1:j])
            in_word, i = True, j + 1
        elif c == '"':
            in_word, i = True, i + 1
            while True:
                if i >= n:
                    raise CfgError(f'{where} cannot be parsed as DIA-NN flags (unclosed ")')
                c = text[i]
                if c == '"':
                    i += 1
                    break
                if c == "\\" and i + 1 < n and text[i + 1] in '$`"\\\n':
                    if text[i + 1] != "\n":
                        cur.append(text[i + 1])
                    i += 2
                else:
                    cur.append(c)
                    i += 1
        elif c == "\\":
            if i + 1 >= n:
                cur.append(c)
                in_word, i = True, i + 1
            elif text[i + 1] == "\n":
                i += 2
            else:
                cur.append(text[i + 1])
                in_word, i = True, i + 2
        else:
            cur.append(c)
            in_word, i = True, i + 1
    if in_word:
        toks.append("".join(cur))
    return toks


def cfg_tokens(cfg):
    """THE tokeniser for a DIA-NN cfg. Every reader goes through it: the parallel gate, the
    step flags, params.base.cfg, the single-shot command and ensure_xic in run_search.py.

    There used to be several rules for one file: read_cfg_flags() dropped flags per LINE
    (`line.startswith("--window ")`), mass_acc_status() found them per shlex TOKEN, the
    params.base.cfg writer used `split()[:1]`, and run_search used `txt.split()`. They
    disagreed on ordinary cfgs:

        --qvalue 0.01 --window 0    the token rule saw `--window 0` and sent the cfg to step
                                    1b; the line rule left it in the step flags, next to the
                                    measured `--window $(cat window.txt)`
        --qvalue 0.01  # 1% FDR     the token rule counted every later flag as set; the line
                                    rule spliced the `#` into the joined bash line, which
                                    comments out EVERY flag after it
        # --xic 10                  split() saw --xic, so ensure_xic added none, while the
                                    chain saw no --xic -- and step 4 extracted no XICs

    so the gate approved flags the generated steps did not carry (CLAUDE.md rule 3). The rule
    is bash's (_split_cfg_text), because bash is what finally reads these words. Returns []
    when no cfg was given; raises CfgError for a path that is not a regular file or that cannot
    be parsed -- never a fallback split, which is how a second answer creeps back in."""
    if not cfg:
        return []
    if not os.path.isfile(cfg):
        raise CfgError(f"cfg not found: {cfg}" if not os.path.exists(cfg)
                       else f"cfg is not a regular file: {cfg}", code="cfg_missing")
    with open(cfg) as fh:
        return _split_cfg_text(fh.read(), cfg)


def cfg_groups(tokens):
    """[(flag, [values])]: a `--flag` owns every following token up to the next `--flag`,
    wherever the line breaks fall. DIA-NN values never start with `--` (a negative number is
    `-3`, which stays a value)."""
    groups = []
    for t in tokens:
        if t.startswith("--") or not groups:
            groups.append((t, []))
        else:
            groups[-1][1].append(t)
    return groups


# `$NAME` / `${NAME}` in a cfg value still expand in the chain, deliberately and ONLY in that
# form: tests/test_hive_submission_guards.py pins `--lib-dir $HOME/libs`, and a hand-written cfg
# may rely on it. `$(...)`, backticks and every other `$` are emitted literally.
_SHELL_VAR = re.compile(r"\$(?:[A-Za-z_][A-Za-z0-9_]*|\{[A-Za-z_][A-Za-z0-9_]*\})")


def _shield(token):
    """One cfg token as ONE bash word whose value is exactly the token.

    A token made only of characters bash never reinterprets (shlex.quote's own safe set:
    letters, digits, `@%+=:,./-_`) is emitted bare, so a cfg with nothing at risk -- every cfg
    estimate_params.py writes, apart from the glob below -- yields the same command line as
    before. Anything else is quoted:

      * globs. DIA-NN's own recommended `--cut K*,R*` is the trypsin rule and an N-terminal mod
        is `UniMod:1,42.010565,*n`. Bare, one matching file rewrites them. Verified:
            no matching file      --cut K*,R*      -> --cut K*,R*
            file 'Kfoo,Rbar'      --cut K*,R*      -> --cut Kfoo,Rbar
      * parentheses. `--var-mod "Phospho(STY),79.966331,STY"` loses its quotes in the
        tokeniser; the previous rule re-emitted it bare, and EVERY step died on bash's
        "syntax error near unexpected token `('" -- after the gate had approved the cfg.
      * braces and tildes: `x{1,2}` would become two words and `~/libs` a home directory.
      * whitespace, `#`, `;&|<>`, quotes, backslash, `!`, `$`: re-split, comment, control
        operator, or expansion.

    A token carrying `$NAME`/`${NAME}` is double-quoted with everything else escaped, so the
    variable expands and nothing else does; any other token is single-quoted."""
    if token and shlex.quote(token) == token:
        return token
    if not _SHELL_VAR.search(token):
        return shlex.quote(token)
    out, pos = ['"'], 0
    esc = lambda s: re.sub(r'([\\"`$])', r"\\\1", s)
    for m in _SHELL_VAR.finditer(token):
        out += [esc(token[pos:m.start()]), m.group(0)]
        pos = m.end()
    out += [esc(token[pos:]), '"']
    return "".join(out)


def bash_flags(groups, drop=()):
    """Flag groups -> one bash-safe string. The one emitter for cfg flags spliced into a
    command line: the chain's steps (read_cfg_flags) and run_search's single-shot search."""
    names = set(drop)
    return " ".join(_shield(t) for flag, vals in groups if flag not in names
                    for t in (flag, *vals))


def read_cfg_flags(cfg, drop=()):
    """Read a diann.cfg into a flat, bash-safe flag string, dropping step-specific flags.

    `drop` removes further flags on top of STRIP. It exists for --window: when step 1b
    measures the radius, steps 2-5 get it PREFIXED as `--window $(cat window.txt)`, so a
    --window still in the cfg lands on the same command line twice. Which of the two DIA-NN
    honours is not verified -- and if it is the cfg's, step 1b is silently undone. A cfg
    `--window 0` is not even a radius: on the 18-file poplar run DIA-NN logged "scan window
    radius should be a positive integer" and chose a radius per file (7 for seventeen, 8 for
    one), which is exactly what the chain exists to prevent.

    Flags are matched by TOKEN (cfg_tokens) and removed with their values, wherever they sit
    on a line; `# comments` never reach bash."""
    return bash_flags(cfg_groups(cfg_tokens(cfg)), drop=tuple(STRIP) + tuple(drop))


def write_cfg(cfg, dest, drop=()):
    """Write `cfg` minus the `drop` flags to `dest`, one flag per line.

    Written from the same tokens every other reader sees, so what params.base.cfg leaves out
    is exactly what the gate and the step flags left out. Comments are not carried over. This
    file is read back by DIA-NN (`--cfg`), not by bash, so glob and parenthesis values stay
    bare -- quoting `K*,R*` here would hand DIA-NN the quote characters. Only a value the
    tokeniser itself would split or strip (whitespace, a quote, a backslash, a leading `#`, or
    empty) is quoted, so cfg_tokens reads the file back identically; how DIA-NN's own cfg
    reader treats quotes is not verified."""
    names = set(drop)
    special = re.compile(r"""[\s'"\\]|^#""")
    with open(dest, "w") as fh:
        for flag, vals in cfg_groups(cfg_tokens(cfg)):
            if flag not in names:
                fh.write(" ".join([flag] + [shlex.quote(v) if not v or special.search(v)
                                            else v for v in vals]) + "\n")


def xic_flag(cfg):
    """Return the '--xic N' flag (plus --mobilograms when requested), or '' if not.

    XIC extraction must happen on the FINAL per-file pass (step 4): step 2 works
    against the predicted library so its IDs are not final, and step 5 runs with
    --use-quant, which reuses .quant files without re-reading the raw spectra --
    DIA-NN accepts --xic there and logs that it will extract, but writes nothing.
    """
    groups = cfg_groups(cfg_tokens(cfg))
    out = ""
    for flag, vals in groups:
        if flag == "--xic":
            out = f"--xic {_shield(vals[0])}" if vals and not vals[0].startswith("-") else "--xic"
    # --mobilograms must ride along with --xic or the mobilogram parquets are written
    # full of zeros. Emitted only when XICs are requested; DIA-NN ignores it on
    # instruments without ion mobility.
    if out and any(flag == "--mobilograms" for flag, _ in groups):
        out += " --mobilograms"
    return out


MASS_ACC_FLAGS = ("--mass-acc", "--mass-acc-ms1")


def _passed(groups, flag):
    """Each occurrence of `flag` exactly as the cfg passes it, e.g. ["--window 7.0"]."""
    return [" ".join([flag, *vals]) for f, vals in groups if f == flag]


def _mass_acc_value(groups, flag):
    """(state, ppm, problem) for one mass-accuracy flag; state is 'unset' | 'ok' | 'invalid'.

    0 is INVALID, not auto: DIA-NN reads `--mass-acc 0` as a literal 0 ppm tolerance ("Mass
    accuracy will be fixed to 0 (MS2) and 0 (MS1)") and a 28-run Lumos search returned 0 IDs
    at 1% FDR (estimate_params.py, PR #38). Auto-calibration is the flag ABSENT."""
    seen = [vals for f, vals in groups if f == flag]
    if not seen:
        return "unset", None, None
    ppms = []
    for vals in seen:
        shown = f"{flag} {' '.join(vals)}".strip()
        if len(vals) != 1:
            return "invalid", None, f"`{shown}` needs exactly one value"
        try:
            v = float(vals[0])
        except ValueError:
            return "invalid", None, f"`{shown}` is not a number"
        if not math.isfinite(v):
            return "invalid", None, f"`{shown}` is not a finite number"
        if v == 0:
            return "invalid", None, (f"`{shown}` is a literal 0 ppm tolerance in DIA-NN "
                                     "(0 IDs), not auto")
        if v < 0:
            return "invalid", None, f"`{shown}` is negative"
        ppms.append(v)
    if len(set(ppms)) > 1:
        return "invalid", None, (f"{flag} is set {len(ppms)} times with different values "
                                 f"({', '.join(map(str, ppms))}); which one DIA-NN uses is "
                                 "not verified")
    return "ok", ppms[0], None


def _window_value(groups):
    """(state, radius, problem) for --window; state is 'unset' | 'zero' | 'ok' | 'invalid'.

    DIA-NN requires a POSITIVE INTEGER ("scan window radius should be a positive integer").
    So only digits count: `0.5` used to read as a truthy float and route straight to DIA-NN,
    and `nan` crashed int() with a traceback. 0 is its own state because it behaves like the
    flag being absent -- on the poplar run DIA-NN logged that warning and chose a radius per
    file -- which is recoverable by measuring, not a typo."""
    seen = [vals for f, vals in groups if f == "--window"]
    if not seen:
        return "unset", None, None
    radii = []
    for vals in seen:
        if len(vals) != 1 or not re.fullmatch(r"[0-9]+", vals[0]):
            shown = f"--window {' '.join(vals)}".strip()
            return "invalid", None, f"`{shown}` is not a positive integer"
        radii.append(int(vals[0]))
    if len(set(radii)) > 1:
        return "invalid", None, (f"--window is set {len(radii)} times with different values "
                                 f"({', '.join(map(str, radii))})")
    if radii[0] == 0:
        return "zero", 0, ("`--window 0` is not a positive integer (on the poplar run DIA-NN "
                           "warned \"scan window radius should be a positive integer\" and "
                           "chose a radius per file)")
    return "ok", radii[0], None


def mass_acc_status(cfg):
    """Read the per-file-optimised settings out of a cfg, validated. Raises CfgError.

    Steps 3/5 reuse the .quant files from steps 2/4, so anything DIA-NN auto-optimises per
    file is applied inconsistently between passes and then stitched together. DIA-NN says so
    itself:

        WARNING: combining reuse of .quant files with automatic optimisation of mass
        accuracies OR SCAN WINDOW will lead to results that are different from those
        of the original analysis that produced the .quant files and is strongly not
        recommended

    So this reads BOTH mass accuracy and --window. Checking only mass accuracy was a real
    defect: on an 18-file poplar run with mass accuracy correctly pinned, DIA-NN still
    inferred a scan-window radius of 7 for seventeen files and 8 for one, emitted the warning
    above, and the chain combined them anyway.

    This only READS; parallel_safe() decides. Mass accuracy and the window are reported
    SEPARATELY (mass_acc_fixed/mass_acc_reason vs window_state/window_passed/window_reason):
    one combined `fixed: false, reason: "not set: --window"` for a cfg whose mass accuracy WAS
    pinned is the wording that caused the original routing bug, and it had reappeared in the
    provenance. `state` maps each flag to unset/ok/invalid (and zero, for --window)."""
    groups = cfg_groups(cfg_tokens(cfg))
    ms2 = _mass_acc_value(groups, "--mass-acc")
    ms1 = _mass_acc_value(groups, "--mass-acc-ms1")
    win = _window_value(groups)
    per = {"--mass-acc": ms2, "--mass-acc-ms1": ms1, "--window": win}
    state = {k: v[0] for k, v in per.items()}
    unset = [k for k, s in state.items() if s == "unset"]
    bad = [k for k, s in state.items() if s == "invalid"]

    ma_unset = [k for k in MASS_ACC_FLAGS if state[k] == "unset"]
    ma_problems = [per[k][2] for k in MASS_ACC_FLAGS if per[k][2]]
    if ma_unset:
        ma_problems.insert(0, "not in the cfg: " + ", ".join(ma_unset)
                           + " (DIA-NN calibrates it itself)")
    ma_fixed = not ma_problems
    ma_reason = ("; ".join(ma_problems) if ma_problems else
                 f"MS1 {ms1[1]} ppm / MS2 {ms2[1]} ppm, pinned in the cfg")
    win_reason = win[2] or ("--window not in the cfg" if win[0] == "unset" else None)

    reason = "; ".join(x for x in ([] if ma_fixed else [ma_reason]) + [win_reason] if x) or \
        f"fixed (MS1 {ms1[1]} ppm / MS2 {ms2[1]} ppm / window {win[1]})"
    return {"ms2": ms2[1], "ms1": ms1[1], "window": win[1], "state": state,
            "bad": bad, "unset": unset,
            "mass_acc_fixed": ma_fixed, "mass_acc_reason": ma_reason,
            "window_state": win[0], "window_passed": _passed(groups, "--window"),
            "window_reason": win_reason, "reason": reason,
            # the cfg searches DDA: nothing in it can be measured by step 1b (DDA_WINDOW_NOTE)
            "dda": any(f == DIANN_DDA_FLAG for f, _ in groups)}


def mass_acc_defaults(cfg, ma):
    """{flag: ppm} of the cfg's mass-accuracy levels that estimate_params.py pinned at the SOP as
    a DEFAULT (a DDA level with no DIA-NN table value: estimate_params.dda_sop_levels), from the
    cfg's rationale sidecar -- and only while the cfg still passes that very value, so a cfg edited
    since cannot inherit a label that no longer describes it."""
    try:
        with open(cfg + ".rationale.json") as fh:
            side = json.load(fh)
    except (OSError, ValueError, TypeError):
        return {}
    got = side.get("mass_accuracy_default") if isinstance(side, dict) else None
    now = {"--mass-acc": ma.get("ms2"), "--mass-acc-ms1": ma.get("ms1")}
    return {f: v for f, v in (got if isinstance(got, dict) else {}).items()
            if f in now and isinstance(v, (int, float)) and not isinstance(v, bool)
            and now[f] is not None and float(now[f]) == float(v)}


def mass_acc_record(ma, cfg=None):
    """Mass accuracy only, for the generator's output and search_provenance.json. A level the
    cfg's sidecar calls a DEFAULT (mass_acc_defaults) is named as one: rule 2 of CLAUDE.md."""
    levels = [ma["state"][f] for f in MASS_ACC_FLAGS]
    default = mass_acc_defaults(cfg, ma) if cfg else {}
    rec = {"mode": ("invalid" if "invalid" in levels else
                    "auto" if levels.count("unset") == 2 else "partial" if "unset" in levels else
                    "pinned_default" if default else "pinned"),
           "fixed": ma["mass_acc_fixed"], "ms1": ma["ms1"], "ms2": ma["ms2"],
           "reason": ma["mass_acc_reason"]}
    if default:
        rec["default"] = default
        rec["default_note"] = (
            "DEFAULT, not user-confirmed: " + ", ".join(f"{f} {v:g}" for f, v in
                                                        sorted(default.items()))
            + " pinned at the facility SOP by estimate_params.py for a DDA search, which cannot "
              "measure it -- see mass_accuracy_default in the cfg's .rationale.json")
    return rec


def window_record(ma):
    """What the cfg hands DIA-NN for --window, for provenance -- exactly what was passed, and
    "unverified" wherever DIA-NN's behaviour has not been measured. Used for the chain when it
    does not probe, and by run_search.py for the single-shot search."""
    st, passed = ma["window_state"], ma["window_passed"]
    if st == "ok":
        return {"mode": "pinned", "source": f"pinned in the cfg ({'; '.join(passed)})",
                "value": ma["window"], "passed": passed}
    if st == "unset" and ma.get("dda"):
        return {"mode": "auto", "source": DDA_WINDOW_NOTE, "value": None, "passed": []}
    if st == "unset":
        return {"mode": "auto",
                "source": "not in the cfg -- DIA-NN chooses the radius itself (on the 18-file "
                          "poplar chain it chose per file: 7 for seventeen, 8 for one; how it "
                          "chooses within one multi-file search is unverified)",
                "value": None, "passed": []}
    if st == "zero":
        # not `auto`: the flag IS set, to a value DIA-NN itself rejects (SCAN_WINDOW_MODES)
        return {"mode": "invalid",
                "source": "passed as `--window 0`, which is not a positive integer -- on the "
                          "poplar run DIA-NN warned and chose a radius per file; unverified "
                          "beyond that run",
                "value": None, "passed": passed}
    return {"mode": "invalid",
            "source": f"passed as given ({'; '.join(passed)}) -- not one positive integer, so "
                      "what DIA-NN does with it is unverified",
            "value": None, "passed": passed}


def cfg_acquisition(cfg):
    """The acquisition estimate_params.py wrote `cfg` for (`acquisition` in its
    <cfg>.rationale.json), or None when there is no readable sidecar."""
    try:
        with open(cfg + ".rationale.json") as fh:
            side = json.load(fh)
    except (OSError, ValueError, TypeError):
        return None
    acq = side.get("acquisition") if isinstance(side, dict) else None
    return acq if isinstance(acq, str) else None


def dda_mismatch(cfg, acquisition, source="the bundle's acquisition"):
    """Why this cfg must not search data of this acquisition, or None. THE check, for both of
    run_search.py's DIA-NN routes (the 5-step chain and the single-shot search).

    DIA-NN's README: "--dda process data as DDA -- must be used with DDA data, must not be used
    with DIA data". The flag lives in the cfg (estimate_params.py writes it for DDA) and nowhere
    else: run_search.py used to append it to the single-shot command only, so a DDA cohort of more
    than 5 files went to the chain without it and was searched as DIA, silently (SET28, 28 Exploris
    DDA .raw, caught by a manual grep). A cfg that disagrees with the bundle is refused, not
    patched: patching the command would leave the recorded cfg describing a different search.
    An acquisition that is not DIA or DDA (empty, "unknown", "mixed") is not checked.
    Raises CfgError for a cfg that cannot be read."""
    acq = (acquisition or "").strip().upper()
    if acq not in ("DIA", "DDA"):
        return None
    has = any(f == DIANN_DDA_FLAG for f, _ in cfg_groups(cfg_tokens(cfg)))
    if is_dda(acq) and not has:
        return (f"{source} is DDA, but {cfg} has no {DIANN_DDA_FLAG}, so DIA-NN "
                "would search the DDA spectra as DIA -- with no error. Fix: re-run "
                "estimate_params.py --engine diann --acquisition DDA (it writes "
                f"{DIANN_DDA_FLAG}, the MS1 survey range and a DDA mass accuracy), or add "
                f"{DIANN_DDA_FLAG} to the cfg.")
    if acq == "DIA" and has:
        return (f"{source} is DIA, but {cfg} has {DIANN_DDA_FLAG}, which DIA-NN "
                "says must not be used with DIA data. Fix: re-run estimate_params.py --engine "
                f"diann --acquisition DIA, or remove {DIANN_DDA_FLAG} from the cfg -- or, if the "
                "data really is DDA, re-run estimate_params.py with --acquisition DDA and correct "
                "the bundle's acquisition (step 4).")
    return None


def mass_acc_measure_plan(cfg):
    """Is this cfg's mass accuracy to be MEASURED with DIA-NN before the search? Returns
    {"documented": {flag: ppm}} -- the levels that have a documented DIA-NN value, pinned as
    documented -- or None.

    The plan travels with the cfg as its rationale sidecar, `<cfg>.rationale.json` -- the file
    estimate_params.py already writes next to every cfg and provenance.py already copies into the
    reproducibility bundle. It says `measure_with_diann` only for an Orbitrap with a level outside
    DIA-NN's table (see estimate_params.mass_acc_plan). None, and so nothing measured, when:
      * there is no sidecar, or it says anything else -- a hand-written cfg that merely forgot
        mass accuracy is a mistake to report, not something to measure over (on a timsTOF the
        documented 15/15 is right, and a measured value would quietly replace it);
      * the cfg sets either flag after all, validly or not (mass_acc_status "unset" for BOTH is
        required) -- the plan no longer describes it, and a lone flag fixes BOTH levels in DIA-NN
        2.7.0, so there would be nothing left to measure;
      * the cfg cannot be read (CfgError -- parallel_safe reports that itself);
      * a documented value is not one positive number per flag, or documents both levels.
    Both callers -- parallel_safe() for the chain, run_search.run_diann() for a single-shot
    search -- read the plan here, so they cannot disagree about it."""
    if not cfg:
        return None
    try:
        with open(cfg + ".rationale.json") as fh:
            side = json.load(fh)
    except (OSError, ValueError):
        return None
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    from estimate_params import MEASURE_WITH_DIANN     # one name for the plan, one file
    if not isinstance(side, dict) or side.get("mass_accuracy_plan") != MEASURE_WITH_DIANN:
        return None
    try:
        st = mass_acc_status(cfg)["state"]
    except CfgError:
        return None
    if any(st[f] != "unset" for f in MASS_ACC_FLAGS):
        return None
    doc = side.get("mass_accuracy_documented") or {}
    if not isinstance(doc, dict) or set(doc) - set(MASS_ACC_FLAGS) or len(doc) == 2:
        return None
    for v in doc.values():
        if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) \
                or not v > 0:
            return None
    return {"documented": dict(doc)}


MASS_ACC_DDA_REASON = (
    "mass accuracy is to be measured with DIA-NN before the search (planned by "
    f"estimate_params.py), but the cfg searches DDA ({DIANN_DDA_FLAG}) and the probe cannot "
    "measure anything under it: no scan-window radius is logged in DDA mode, and mass-accuracy "
    "optimisation did not finish within 3600 s on the SET28 runs")


def mass_acc_dda_refusal(cfg):
    """{code, reason, remedy} when `cfg` plans a DIA-NN mass-accuracy measurement AND searches
    DDA, else None. The one rule for both routes: parallel_safe() declines the chain with it, and
    run_search.py refuses the single-shot search with it at GENERATION -- the single-shot job
    used to carry the probe into the job, where probe_window.py refused --dda at run time.
    The plan comes from a sidecar written before 2.9 (estimate_params.py pins a DDA cfg now) or
    a cfg given --dda by hand."""
    if not mass_acc_measure_plan(cfg):
        return None
    try:
        if not mass_acc_status(cfg)["dda"]:
            return None
    except CfgError:
        return None
    return {"code": "mass_acc_dda", "reason": MASS_ACC_DDA_REASON,
            "remedy": _remedy("mass_acc_dda")}


def window_flag(path):
    """Bash for the --window steps 2-5 pass: `--window N` from `path`, or nothing when it says
    `auto` (the probe's fallback: DIA-NN chooses the radius per run). needs_measured() has
    already refused anything else, a missing file included."""
    return f"$(sed -n 's/^\\([1-9][0-9]*\\)$/--window \\1/p' \"{path}\") "


def keep_attempt(workdir):
    """probe_attempts() `reset` lines that keep attempt 1's probe logs, as `<workdir>.attempt1`,
    beside its evidence (`<evidence>.attempt1`) -- they are what says why it measured nothing."""
    q = shlex.quote
    return [f"rm -rf {q(workdir + '.attempt1')}",
            f"mv -f {q(workdir)} {q(workdir + '.attempt1')} 2>/dev/null || true"]


def _rc_case(codes):
    """A bash `case` pattern for these exit statuses."""
    return "|".join(str(c) for c in codes)


def probe_failure_lines(indent="  "):
    """Bash: why the probe stopped the job, from $PROBE_RC (probe_window.EXIT_MEANING)."""
    q = shlex.quote
    return ([f'{indent}case "$PROBE_RC" in']
            + [f"{indent}  {code}) echo {q('  ' + text)} >&2 ;;"
               for code, text in sorted(EXIT_MEANING.items())]
            + [f'{indent}  {_rc_case(FALLBACK_ON)}) echo "  the fallback could not be written '
               '(above)" >&2 ;;',
               f'{indent}  *) if [ "$PROBE_RC" -ge 128 ]; then echo "  the probe was stopped by '
               'a signal (exit $PROBE_RC)" >&2; else echo "  the probe exited $PROBE_RC" >&2; '
               'fi ;;',
               f"{indent}esac"])


def probe_attempts(head, tail, evidence, measure, documented, reset, window_file=None,
                   massacc_file=None, write_cfg=None, provenance=None, fallback_out=None):
    """Bash that runs the probe, runs it ONCE more when its own machinery failed, and then falls
    back (probe_fallback.py) instead of failing the search -- only then. The one definition, for
    the chain's step 1b and the single-shot search's probe.

    A failed probe used to fail its job, and in the chain afterok then left steps 2-5
    DependencyNeverSatisfied -- 62 of fran-5b's 261 step-1b jobs overnight 2026-09-29/30 died
    that way, over an ESTALE in the probe's log tail while DIA-NN itself succeeded. The probe's
    exit status now says why it pinned nothing (probe_window.EXIT_*), and:
      0                   measured -- the caller takes the values from `evidence`
      RETRY_ON            a crash, or its own log unreadable (io_error): retried once when at
                          least PROBE_RETRY_MIN_S of the budget is left, then probe_fallback.py
      FALLBACK_ON, else   a time limit: probe_fallback.py, no retry (a second attempt from what
                          is left of the same budget would hit it again)
      anything else       the job FAILS as before -- a refused measurement, the environment
                          (no .NET, DIA-NN cannot start), the probe's arguments, runs DIA-NN
                          finished without logging it, a signal: no fallback can fix those, and
                          one would hide them (dda-review, 2026-09-30)
    `head` is the probe's command up to (not including) `--budget`; `tail` the DIA-NN flags
    after `--`. `reset` clears an attempt's leftovers before the retry. Afterwards PROBE_RC is
    the probe's last exit status and PROBE_FALLBACK is 1 when the fallback was written. Safe
    under `set -e` (the single-shot inline route runs these lines with it)."""
    q = shlex.quote
    fb = os.path.join(os.path.dirname(os.path.abspath(__file__)), "probe_fallback.py")
    doc = f" {probe_mass_acc_args(documented)}" if documented else ""
    opt = "".join(f" {flag} {q(v)}" for flag, v in (
        ("--window-file", window_file), ("--massacc-file", massacc_file),
        ("--write-cfg", write_cfg), ("--provenance", provenance)) if v)
    return [
        f"PROBE_DEADLINE=$(( $(date +%s) + {PROBE_BUDGET_S} )); PROBE_RC=1; PROBE_FALLBACK=0; "
        "PROBE_ATTEMPTS=0",
        "for PROBE_ATTEMPT in 1 2; do",
        "  PROBE_BUDGET=$(( PROBE_DEADLINE - $(date +%s) ))",
        '  if [ "$PROBE_ATTEMPT" -gt 1 ]; then',
        f'    if [ "$PROBE_BUDGET" -lt {PROBE_RETRY_MIN_S} ]; then echo "[probe] not retried: '
        'only $PROBE_BUDGET s of the budget left" >&2; break; fi',
        '    echo "[probe] WARNING: attempt 1 failed in the probe\'s own machinery (exit '
        f'$PROBE_RC); retrying once. Its evidence: {evidence}.attempt1" >&2',
        f"    mv -f {q(evidence)} {q(evidence + '.attempt1')} 2>/dev/null || true",
        *("    " + r for r in reset),
        "  fi",
        "  PROBE_ATTEMPTS=$PROBE_ATTEMPT",
        f"  {head} --budget $PROBE_BUDGET -- {tail} > {q(evidence)} && PROBE_RC=0 || PROBE_RC=$?",
        f'  case "$PROBE_RC" in {_rc_case(RETRY_ON)}) ;; *) break ;; esac',
        "done",
        f'case "$PROBE_RC" in {_rc_case(FALLBACK_ON)})',
        f"  if python3 {q(fb)} --measure {' '.join(measure)}{doc} --exit-code $PROBE_RC "
        f"--attempts $PROBE_ATTEMPTS --evidence {q(evidence)}{opt}"
        + (f" --out {q(fallback_out)}" if fallback_out else "") + "; then PROBE_FALLBACK=1; fi ;;",
        "esac",
    ]


def probe_mass_acc_args(documented):
    """probe_window.py flags that pin the documented levels: {"--mass-acc-ms1": 7} -> "--ms1-ppm 7"."""
    names = {"--mass-acc-ms1": "--ms1-ppm", "--mass-acc": "--ms2-ppm"}
    return " ".join(f"{names[f]} {v:g}" for f, v in sorted(documented.items(), reverse=True))


def documented_levels_text(documented):
    """The levels a plan pins as given, for a verdict or a provenance record. "Documented" here
    means a value from DIA-NN's Orbitrap resolution table -- an exact tier, OR interpolated between
    tiers (MS1 at 90k gives 7.5 ppm) -- so the text must not call it a README value: it said
    "--mass-acc-ms1 7.5 as documented". The cfg's rationale sidecar names which each level is."""
    return (", ".join(f"{f} {v:g}" for f, v in sorted(documented.items()))
            + " from DIA-NN's Orbitrap resolution table (an exact tier, or interpolated between "
              "tiers -- `mass_accuracy_documented` in the cfg's .rationale.json says which)")


# What step 1b writes, as steps 2-5 read it. A pattern for `grep -Eqx` and Python's re alike.
# It checks the SHAPE of the line -- the two flags, in order, each with a number. A shape check
# cannot tell 14 ppm from 999999, so massacc.txt gets massacc_band_check() as well.
MEASURED_FILE_RE = {
    # a radius. `auto` (probe_fallback.py: DIA-NN chooses per run, and steps 2-5 pass no
    # --window -- window_flag()) is accepted by needs_measured() only beside the fallback's own
    # record, never on its own
    "window.txt": "[1-9][0-9]*",
    "massacc.txt": "--mass-acc [0-9]*[.]?[0-9]+ --mass-acc-ms1 [0-9]*[.]?[0-9]+",
}


def massacc_band_check(path, producer):
    """Bash that fails the job unless both numbers in `path` are inside MASS_ACC_BAND.

    The probe cannot write an out-of-band value (probe_window.pin_mass_acc refuses to pin one),
    so this is the guard for every OTHER way the file can hold one: a hand-edited massacc.txt, a
    file left by an older skill version, a copy from another cohort. Steps 2-5 splice
    `$(cat massacc.txt)` straight onto a DIA-NN command line, so the number the shell expands is
    the number that gets searched at, and nothing downstream looks at it again."""
    (ms2_lo, ms2_hi), (ms1_lo, ms1_hi) = MASS_ACC_BAND["ms2_ppm"], MASS_ACC_BAND["ms1_ppm"]
    band = (f"MS2 {band_text('ms2_ppm')}, MS1 {band_text('ms1_ppm')}")
    awk = ("awk 'NR==1{ok = ($2 >= %g && $2 <= %g && $4 >= %g && $4 <= %g)} END{exit !ok}'"
           % (ms2_lo, ms2_hi, ms1_lo, ms1_hi))
    return (f'if ! {awk} "{path}" 2>/dev/null; then '
            f'echo "FAILED: $(cat {shlex.quote(path)}) in {path} is outside the plausible band '
            f'for an Orbitrap ({band}). {producer} cannot produce that, so the file was edited, '
            f'copied from another cohort or written by an older version. Delete it and re-run '
            f'{producer}; DIA-NN would search every precursor at that tolerance." >&2; exit 1; fi')


def needs_measured(path, what, producer="step 1b (step1b_window.sbatch)"):
    """Bash that fails the job unless step 1b's `path` holds a measurement.

    Steps 2-5 splice `$(cat massacc.txt)` and `--window N` from window.txt (window_flag();
    nothing for the fallback's `auto`) into DIA-NN's command line. A missing file expands to NOTHING: DIA-NN then optimises mass accuracy per file -- the
    very thing the chain exists to prevent -- and still writes its .quant, so must_exist passes
    and steps 3/5 stitch auto-calibrated passes together. That is not hypothetical: step 1b
    deletes both files before it probes, and references/watcher.md tells the orchestrator to
    resubmit the downstream steps of a stalled chain (review reproduction: `sbatch
    step2_firstpass.sbatch` after a failed step 1b ran DIA-NN with no mass-accuracy flag). The
    content is checked too, so a "None" printed into the file cannot reach DIA-NN."""
    base = os.path.basename(path)
    pat = MEASURED_FILE_RE[base]
    ok = f'grep -Eqx -- {shlex.quote(pat)} "{path}" 2>/dev/null'
    if base == "window.txt":
        # `auto` only with probe_fallback.py's record of a window fallback beside it: a
        # hand-written or stale `auto` must not make steps 2-5 run without a window unrecorded
        fb = os.path.join(os.path.dirname(path), FALLBACK_RECORD)
        marker = shlex.quote('"mode": "fallback_auto"')
        ok = (f'{{ {ok} || {{ grep -qx auto "{path}" 2>/dev/null && '
              f'grep -q {marker} "{fb}" 2>/dev/null; }}; }}')
    shape = (f'if ! {ok}; then '
             f'echo "FAILED: {path} does not hold {what} -- {producer} measures '
             f'it and has not completed successfully. Re-run it first; running DIA-NN '
             f'without it would let DIA-NN optimise it itself." >&2; exit 1; fi')
    # The shape is not the value: `--mass-acc 999999 --mass-acc-ms1 0.001` matches the pattern.
    if base == "massacc.txt":
        return shape + "\n" + massacc_band_check(path, producer)
    return shape


# The refusals that mean "mass accuracy is OMITTED and nothing will measure it" -- the only ones
# --allow-auto-mass-acc may override (upstream #70: never an invalid value, a bad --window or an
# unreadable cfg). A cfg estimate_params.py planned to measure is omitted mass accuracy too, when
# the chain has no step 1b to measure it in.
MASS_ACC_OMITTED_CODES = ("mass_acc_unset", "mass_acc_seeded", "mass_acc_no_probe",
                          "mass_acc_dda")


def _remedy(code):
    """How to fix a refusal -- derived from WHY it was refused, in one place, so the router's
    decline and the generator's exit say the same thing. Both used to hardcode "re-run
    estimate_params.py with the instrument table" whatever the cause, and the generator also
    told users to hand-run probe_window.py for a window step 1b measures by itself."""
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    from estimate_params import instrument_ppm_summary         # the ONE ppm table
    table = instrument_ppm_summary()
    return {
        "cfg_missing":
            "check the cfg path (diann_parallel --cfg / run_search --params) -- nothing "
            "could be read from it",
        "cfg_unparseable":
            "fix the quoting in the cfg (every quote must be closed), then re-run",
        "mass_acc_unset":
            "pin mass accuracy: re-run estimate_params.py with the real instrument "
            f"({table}); for an Orbitrap pass --ms1-resolution/--ms2-resolution (an ion-trap "
            "MS2 has no MS2 resolution: there, pin --mass-acc/--mass-acc-ms1 from a validated "
            "SOP instead). A plan to "
            "measure it is read from the <cfg>.rationale.json estimate_params.py writes beside "
            "the cfg -- keep the two together, and never add just one of the two flags. "
            "Left on auto it can only run as the single-shot search",
        "mass_acc_seeded":
            "pin --mass-acc and --mass-acc-ms1 in the cfg (a validated SOP value): a seeded "
            "chain has no step 1 for step 1b to follow, so nothing can measure them -- or "
            "drop --seed-lib and let step 1b measure them",
        "mass_acc_no_probe":
            "drop --no-probe-window (step 1b then measures mass accuracy), or pin --mass-acc "
            "and --mass-acc-ms1 in the cfg",
        "mass_acc_dda":
            "re-run estimate_params.py --engine diann --acquisition DDA: for DDA it pins the "
            "level it cannot measure at the facility SOP (tagged DEFAULT) instead of planning a "
            "measurement -- or pin --mass-acc and --mass-acc-ms1 in the cfg yourself",
        "mass_acc_invalid":
            "correct --mass-acc/--mass-acc-ms1 in the cfg to ONE positive ppm value each "
            f"({table}), or re-run estimate_params.py with the real instrument. For "
            "auto-calibration delete the flags -- never 0 -- and run single-shot",
        "window_invalid":
            "set --window to one positive integer, or delete it and the chain measures it "
            "itself (step 1b) -- for a DDA cfg (--dda) just delete it: step 1b does not run for "
            "DDA. Do not guess a value: it depends on the acquisition scheme",
        "window_seeded":
            "pin --window in the cfg: a seeded chain has no step 1 for step 1b to follow, so "
            "measure it once with probe_window.py against the seed library, handing it the "
            "cohort's runs (it picks representative ones)",
        "window_no_probe":
            "drop --no-probe-window (step 1b then measures it), or pin --window in the cfg "
            "with a value measured by probe_window.py",
    }.get(code)


def parallel_safe(cfg, probe_window=True, seed_lib=None):
    """THE definition of "may this cfg run as the 5-step chain?" -- for BOTH callers.

    diann_parallel.main() gates generation on `ok`; run_search.py auto-routes on the same
    `ok`. They used to answer separately and drifted: run_search read
    `mass_acc_status()["fixed"]` literally, so it declined every cfg estimate_params.py
    produces (which omits --window by design) and silently demoted the cohort to ONE
    single-shot search. --threads parallelises WITHIN a run, not across runs, so at
    ~30 min/file a 310-file cohort is ~155 h against a few hours for the chain -- and
    nothing errors, SLURM reports success and the user just waits a week. Then the generator
    gated on `ma["fixed"]` while the router used `ok`, which let a duplicated, unparseable
    --mass-acc through one and not the other. Keep the rule here, never in a caller: a second
    copy is what caused both (CLAUDE.md rule 3).

    What is recoverable, and what is not:
      * mass accuracy unset -- recoverable ONLY when step 1b will measure it: both flags
        absent, the cfg's estimate_params.py sidecar plans `measure_with_diann` (an Orbitrap
        with a level outside DIA-NN's table; see mass_acc_measure_plan), and the chain has a
        step 1b (probing on, no seed library -- else mass_acc_seeded / mass_acc_no_probe).
        Step 1b then runs DIA-NN in automatic mode on representative runs and pins the median
        of each measured level, and the documented value of a level that has one, for steps
        2-5 -- so every pass uses the same tolerance, which is what reusing .quant files needs.
        That is a TOLERANCE, not the per-run mass CALIBRATION: DIA-NN recalibrates every run
        whether or not the tolerance is fixed (its log prints "Calibrating with mass accuracies
        25 (MS1), 25 (MS2)" under a pinned 20/7 too). Any other unset mass accuracy -- no
        sidecar, instrument not identified, only one of the two flags -- is NOT recoverable:
        DIA-NN would calibrate it per run, so there is no single value to carry into steps 3/5.
      * mass accuracy 0 / negative / non-numeric / non-finite / set twice differently --
        NOT recoverable, plan or no plan, and not "auto": 0 is a literal 0 ppm tolerance.
      * --window unset or 0 -- recoverable when probing. Step 1b runs probe_window.py and
        pins one radius into steps 2-5. estimate_params.py cannot supply it: the radius is a
        property of the acquisition scheme and has to be MEASURED on a real file. DIA-NN
        does not accept 0 ("scan window radius should be a positive integer") and optimises
        per file instead -- the very inconsistency step 1b removes.
      * --window anything but a non-negative integer (0.5, nan, -1, 7.0, wide), or set twice
        differently -- NOT recoverable. A typo is a mistake to report, not something to
        quietly measure over.
      * an unparseable cfg (unbalanced quote) -- NOT recoverable.
      * a DDA cfg (--dda): nothing is measured, because step 1b's probe cannot measure anything
        under --dda (DDA_WINDOW_NOTE; probe_window.py refuses it). Mass accuracy must be pinned
        -- estimate_params.py pins a DDA cfg -- and a plan to measure it is refused
        (mass_acc_dda). An unset --window is accepted as it is (dda_window_unset): there is no
        radius to measure and nothing to derive one from, and window_record() says so; a
        `--window 0` is refused (window_invalid) -- nothing would replace it.

    Returns {ok, probe, code, ma, reason, remedy, measure, mass_acc_documented}: `probe` says
    step 1b is needed, `measure` lists what it measures ("window", "mass-acc"),
    `mass_acc_documented` the mass-accuracy levels it pins as documented instead ({flag: ppm}),
    `code` names the outcome (probe | pinned | dda_window_unset | cfg_missing | cfg_unparseable |
    mass_acc_unset | mass_acc_seeded | mass_acc_no_probe | mass_acc_dda | mass_acc_invalid |
    window_invalid | window_seeded | window_no_probe), `remedy` is how to fix a refusal. A cfg
    path that is not a file is `cfg_missing`, never "mass accuracy is not pinned": `--sbatch
    proj` once renamed the folder holding the cfg, and the refusal that followed blamed mass
    accuracy.
    `ma` is mass_acc_status() untouched -- on the probe path what step 1b measures is still
    unset, because it IS unset until step 1b runs.
    """
    measure, documented = [], {}

    def verdict(ok, probe, code, ma, reason):
        return {"ok": ok, "probe": probe, "code": code, "ma": ma, "reason": reason,
                "remedy": None if ok else _remedy(code),
                "measure": list(measure) if ok else [],
                "mass_acc_documented": dict(documented) if ok else {}}

    try:
        ma = mass_acc_status(cfg)
    except CfgError as e:
        return verdict(False, False, e.code, None, str(e))
    st, dda = ma["state"], ma["dda"]
    # Invalid values before unset ones: the MASS_ACC_OMITTED_CODES are the ones
    # --allow-auto-mass-acc may override, so they must never hide a junk value behind them --
    # and a plan to measure never rescues a junk value either.
    if any(st[f] == "invalid" for f in MASS_ACC_FLAGS):
        return verdict(False, False, "mass_acc_invalid", ma,
                       f"mass accuracy is set but not usable ({ma['reason']})")
    if st["--window"] == "invalid":
        return verdict(False, False, "window_invalid", ma,
                       f"--window is set but is not a usable radius ({ma['reason']})")
    if dda and st["--window"] == "zero":
        return verdict(False, False, "window_invalid", ma,
                       f"{ma['window_reason']}, and in a DDA cfg ({DIANN_DDA_FLAG}) step 1b does "
                       "not run, so nothing would replace it")
    if any(st[f] == "unset" for f in MASS_ACC_FLAGS):
        plan = mass_acc_measure_plan(cfg)          # None unless BOTH are unset and planned
        if not plan:
            both = all(st[f] == "unset" for f in MASS_ACC_FLAGS)
            return verdict(False, False, "mass_acc_unset", ma,
                           f"mass accuracy is not pinned ({ma['reason']})"
                           + (f"; {cfg}.rationale.json does not plan to measure it with DIA-NN"
                              if both else ""))
        refusal = mass_acc_dda_refusal(cfg)       # SET28's step 1b burned ~3 h of probes on it
        if refusal:
            return verdict(False, False, refusal["code"], ma, refusal["reason"])
        if seed_lib:
            return verdict(False, False, "mass_acc_seeded", ma,
                           "mass accuracy is to be measured in step 1b (planned by "
                           "estimate_params.py), but the first pass is seeded from an existing "
                           "library, so there is no step 1 for step 1b to follow")
        if not probe_window:
            return verdict(False, False, "mass_acc_no_probe", ma,
                           "mass accuracy is to be measured in step 1b (planned by "
                           "estimate_params.py), and --no-probe-window was given")
        measure.append("mass-acc")
        documented = plan["documented"]
    if st["--window"] != "ok" and dda:
        # "unset" (zero was refused above): nothing measures it and nothing is pinned for it --
        # window_record() describes exactly that in the provenance.
        return verdict(True, False, "dda_window_unset", ma,
                       f"MS1 {ma['ms1']} ppm / MS2 {ma['ms2']} ppm, pinned in the cfg; --window "
                       f"not in the cfg and not measured: the cfg searches DDA ({DIANN_DDA_FLAG}), "
                       "where DIA-NN logs no scan-window radius for step 1b to measure")
    if st["--window"] != "ok":
        if seed_lib:
            return verdict(False, False, "window_seeded", ma,
                           "--window is unpinned and the first pass is seeded from an existing "
                           "library, so there is no step 1 for step 1b to follow")
        if not probe_window:
            return verdict(False, False, "window_no_probe", ma,
                           "--window is unpinned and --no-probe-window was given")
        measure.insert(0, "window")
    if not measure:
        return verdict(True, False, "pinned", ma, ma["reason"])
    doc = f"; {documented_levels_text(documented)}" if documented else ""
    got = (f"MS1 {ma['ms1']} ppm / MS2 {ma['ms2']} ppm" if "mass-acc" not in measure else
           "mass accuracy is unpinned but recoverable -- DIA-NN measures it on representative "
           "runs in step 1b (planned by estimate_params.py: no documented value for this "
           f"Orbitrap level{doc})")
    win = (f"window {ma['window']}" if "window" not in measure else
           "--window is unpinned but recoverable -- step 1b measures it")
    return verdict(True, True, "probe", ma, f"{got}; {win}; pinned for steps 2-5")


# DIA-NN 2.6 EXITS 0 ON A FATAL ERROR -- verified: a run against a nonexistent .mzML and
# a nonexistent library prints "ERROR: ..." and returns exit code 0. Since DIA-NN is the
# last command in every step, SLURM records COMPLETED, watch_run.sh reports success, and
# the afterok dependency releases the next step. A step-4 task that dies this way simply
# leaves no .quant, and step 5 then builds the cross-run report from the survivors and
# calls it done -- a silently DROPPED SAMPLE. Exit status is therefore not trustworthy;
# every step must assert that the artefact it was supposed to produce actually exists.
def must_exist(path, what):
    """Bash that fails the job unless `path` exists and is non-empty."""
    return (f'if [ ! -s "{path}" ]; then '
            f'echo "FAILED: DIA-NN exited 0 but did not write {what}: {path}" >&2; '
            f'echo "(DIA-NN 2.6 returns 0 on fatal errors -- check the log above for ERROR:)" >&2; '
            f'exit 1; fi')


def clear_stale(*paths):
    """Bash that deletes `paths` BEFORE DIA-NN runs, so must_exist() can only pass on a file
    this run wrote.

    must_exist() checks that a file is there, not that this run made it. A search re-run into
    the same --out still has the previous run's artefacts, and a DIA-NN that exits 0 having
    written nothing leaves them untouched. Measured on HIVE, DIA-NN 2.7.0 (review srun
    23512013): a re-run with no DOTNET_ROOT logged "ERROR: cannot read .raw files", exited 0,
    the old report.parquet and report.stats.tsv were unchanged by md5, and every guard passed
    on them -- a COMPLETED job reporting results from different parameters. Deleting beats a
    marker file checked with `-newer`: after a failed run the marker approach still leaves a
    plausible report.parquet in --out for the DE step (or anyone listing the folder) to pick
    up, and `-nt` compares whole seconds in some shells (macOS /bin/bash 3.2), so a fast
    failure can look new. Paths go in double quotes, like must_exist(), so a `$VAR` in an
    array-task path still expands -- which is also why a path handed IN has to be checked
    first: see refuse_unsafe_path()."""
    return "rm -f -- " + " ".join(f'"{p}"' for p in paths)


# The double quotes clear_stale() and must_exist() put around a path are deliberate (an array
# task's path carries $QUANT / ${SLURM_ARRAY_TASK_ID}), and bash re-reads `$`, a backtick and
# `\` inside them. Those expansions are built by THIS generator; a path handed in on the
# command line is not entitled to any of them. Reproduced: `--out '/tmp/x$(touch pwned)'` put
# the substitution inside `rm -f -- "..."`, so it EXECUTED when the job ran -- and moved the rm
# target with it. A `"` closes the quote outright. This is refused at the entry point rather
# than quoted at the emitter, because the emitter cannot tell a path it built from one it was
# handed, and `cd "<out>"` in submit.sh has the same hole.
_UNSAFE_IN_PATH = re.compile(r'[$`"\\\n]')


def refuse_unsafe_path(path, flag="--out", prog="diann_parallel"):
    """Exit unless `path` is safe to splice into those double-quoted guards. Returns it."""
    bad = sorted({c for c in _UNSAFE_IN_PATH.findall(path or "")})
    if bad:
        shown = ", ".join(repr(c) for c in bad)
        sys.exit(f"[{prog}] REFUSING {flag} {path!r}: it contains {shown}. Every generated "
                 f"guard puts DOUBLE quotes around a path so that an array task's $QUANT still "
                 f"expands, so a `$(...)` in {flag} would EXECUTE when the job runs and a `\"` "
                 f"would break the command line. Use a path without them.")
    return path


# Wall clock per array task (steps 2 and 4), hours, for a task of TIME_REFERENCE_CPUS CPUs -- the
# default every chain ran at before 2.10 -- and SCALED UP for fewer (array_task_hours). Real Core
# runs (HIVE sacct, every brettsp step-2 task since 2026-07-01 at 16 CPUs, 6,182 completed):
# median 20 min, p95 65 min, p99 160 min, p99.5 248 min; 414 over 60 min, 84 over 120 min, 54
# over 180 min, 32 over 240 min (blank and failed injections are the slowest). The base is 4 h:
# all but those 32 -- the 2 h of 2.9.1 killed ~1% of real files even at 16 CPUs -- and staff pass
# --time-per-file for the rest. 2.10 sizes big arrays to 8 CPUs, measured 1.96x slower per file
# than 16 (a 1.9 GB HT HeLa QC run, 2026-10-01), so the limit scales with them. One timed-out task
# loses the whole chain (steps 3-5 DependencyNeverSatisfied), and a time limit costs no
# throughput (it reserves nothing, at most some backfill), so it is never the thing to save on.
TIME_PER_FILE_HOURS = 4
TIME_REFERENCE_CPUS = 16


def array_task_hours(cpus, base=TIME_PER_FILE_HOURS):
    """Hours per array task of `cpus` CPUs: `base` at TIME_REFERENCE_CPUS, scaled by
    TIME_REFERENCE_CPUS / cpus when the task has fewer (8 h at 8, 16 h at 4), never below `base`
    (16 -> 32 CPUs is only 1.62x faster, so more CPUs earn no shorter limit)."""
    return max(base, math.ceil(base * TIME_REFERENCE_CPUS / max(1, int(cpus))))
# Memory per task, GB: steps 1b and 2 (the predicted library), and step 4 (the smaller empirical
# one). From real Core runs, not the one QC file 2.10 first measured (7.4-9.7 GB RSS): HIVE's sacct
# MaxRSS counts page cache, so the lower bound of a task's own memory is MaxRSS - MaxDiskRead -
# MaxDiskWrite, and for every brettsp chain task since 2026-07-01 that bound is over 32 GB for 65
# step-2 tasks and 9 step-1b jobs, up to 61 GB (a blank-like 1.2 GB .d against the standard human
# library hit its 64 GB limit, 2026-10-01); step 4's stays under 32. 64 GB costs no concurrency on
# genome-center-grp/high: 8 tasks x 64 GB = 512 GB, under the 1 TB per-user cap. --mem-per-file
# sets all three.
MEM_PER_FILE_GB = 64
MEM_FINAL_PASS_GB = 48


def task_memory(a):
    """(GB for steps 1b and 2, GB for step 4): --mem-per-file for all three when given, else
    MEM_PER_FILE_GB and MEM_FINAL_PASS_GB."""
    return ((a.mem_per_file, a.mem_per_file) if a.mem_per_file
            else (MEM_PER_FILE_GB, MEM_FINAL_PASS_GB))


def size_cpus(a, n, array_queue, single_queue):
    """Size the chain's jobs to their queues: a.threads_per_file becomes the array tasks' CPUs
    (run_search.array_task_cpus, the one rule), and a.libpred_cpus / a.assembly_cpus are lowered
    to what one job there can ever get (run_search.fit_to_queue) -- on a queue whose per-job or
    per-user cap is below them they would wait for ever. Says what it chose on stderr. Returns
    (the record for search_provenance.json `cpu_sizing`, step 1b's CPUs)."""
    mem_first, mem_final = task_memory(a)
    try:
        from run_search import user_limits, array_task_cpus, fit_to_queue
        lim = user_limits(*array_queue)
        rec = array_task_cpus(n, a.threads_max, limits=lim, mem_per_task_gb=mem_first,
                              max_simultaneous=a.max_simultaneous, pinned=a.threads_per_file)
        single = user_limits(*single_queue) if single_queue != array_queue else lim
        probe = fit_to_queue(a.threads_per_file or a.threads_max, single)
        lowered = {}
        for flag, attr in (("--libpred-cpus", "libpred_cpus"), ("--assembly-cpus", "assembly_cpus")):
            v = fit_to_queue(getattr(a, attr), single)
            if v < getattr(a, attr):
                lowered[flag] = {"asked": getattr(a, attr), "used": v}
                setattr(a, attr, v)
        if lowered:
            rec["single_jobs_lowered"] = lowered
    except Exception as e:                     # said, never silent: the old fixed sizing stands
        t = a.threads_per_file or a.threads_max
        rec = {"cpus": t, "concurrent": None, "tasks": n, "requested": t,
               "pinned": bool(a.threads_per_file),
               "reason": f"not sized to the queue ({type(e).__name__}: {e}): {t} CPUs per task",
               "summary": f"{t} CPUs per file (not sized to the queue: {e})"}
        probe = t
    a.threads_per_file = rec["cpus"]
    rec["step1b_cpus"] = probe
    if a.time_per_file is None:
        a.time_per_file = array_task_hours(rec["cpus"])
        rule = (f"{TIME_PER_FILE_HOURS} h at {TIME_REFERENCE_CPUS} CPUs per task, scaled to "
                f"{rec['cpus']} (array_task_hours)")
    else:
        rule = "--time-per-file, as given"
    rec["time_per_file_hours"] = {"hours": a.time_per_file, "rule": rule}
    rec["mem_gb"] = {"step1b": mem_first, "step2": mem_first, "step4": mem_final,
                     "rule": ("--mem-per-file, as given" if a.mem_per_file else
                              "MEM_PER_FILE_GB / MEM_FINAL_PASS_GB (sacct, real Core runs)")}
    sys.stderr.write(f"[diann_parallel] array steps 2 and 4: {rec['summary']}, "
                     f"{a.time_per_file} h each. {rec['reason']}."
                     + "".join(f" {f} lowered {v['asked']} -> {v['used']} to fit the queue."
                               for f, v in (rec.get("single_jobs_lowered") or {}).items()) + "\n")
    return rec, probe


def needs_requeue(partition, qos):
    """Should a job on this queue carry `#SBATCH --requeue`?

    Shared by this header() and run_search.emit_sbatch(), which used to disagree (emit_sbatch
    keyed on `partition == "low"` alone and ignored a public QOS). radiant_parallel.header()
    (public QOS only) and diatracer_parallel.py (`low` only) still carry their own rules.

    What it protects, measured on HIVE 2026-09-16 (`scontrol show partition/config`):
    preemption is by partition (PreemptType=preempt/partition_prio; `low` PreemptMode=REQUEUE,
    `high` OFF) and JobRequeue=1, so HIVE requeues a preempted batch job even without this
    line. It is written anyway because JobRequeue=0 on another cluster would make a preempted
    job simply lost. The public-QOS clause also marks publicgrp jobs on `high`, which are not
    preempted there; on such a job the line only allows a requeue after a node failure."""
    return (qos or "").startswith("public") or partition == "low"


def header(name, cpus, mem_gb, hours, partition, account, qos=None, array=None):
    h = ["#!/bin/bash -l",
         f"#SBATCH --job-name={name}",
         f"#SBATCH --cpus-per-task={cpus}",
         f"#SBATCH --mem={mem_gb}G",
         f"#SBATCH --time={hours}:00:00",
         f"#SBATCH --partition={partition}",
         f"#SBATCH --account={account}"]
    if qos:
        h.append(f"#SBATCH --qos={qos}")
    if needs_requeue(partition, qos):
        h.append("#SBATCH --requeue")
    h += [f"#SBATCH -o {name}_%j.log", f"#SBATCH -e {name}_%j.log"]
    if array:
        h.insert(2, f"#SBATCH --array={array}")
        h = [x.replace("_%j.log", "_%A_%a.log") for x in h]
    return "\n".join(h)


def step8_next(report, sess):
    """The command a resumed session runs once the chain is COMPLETED: step 8's check
    (normalization_check.step8_commands), never the final DE without it."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import normalization_check
    return normalization_check.step8_commands(report, f"{sess}/input/conditions.csv", sess,
                                              "dpc")["next"]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diann", required=True, help="DIA-NN command (native binary path, or 'apptainer exec --bind … <sif> /diann-*/diann-linux')")
    ap.add_argument("--raw", nargs="+", default=[], help="raw paths/globs (or use --raw-list)")
    ap.add_argument("--raw-list", help="file with one raw path per line — handles spaces in paths")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cfg", help="diann.cfg with the search params (estimate_params.py output)")
    # CPUs per array task (steps 2 and 4): sized by run_search.array_task_cpus() to the queue's
    # per-user cap unless pinned. 16 was the fixed default; it is now the ceiling.
    ap.add_argument("--threads-per-file", type=int, default=None,
                    help="PIN the CPUs of each array task (steps 2/4) and of step 1b. Default: "
                         "sized to the queue's per-user CPU cap, at most --threads-max")
    ap.add_argument("--threads-max", type=int, default=16,
                    help="the most CPUs an array task is sized to (default 16); step 1b's CPUs")
    ap.add_argument("--mem-per-file", type=int, default=None,
                    help="GB per array task (steps 2/4) and step 1b. Default: MEM_PER_FILE_GB "
                         "(64) for steps 1b and 2, MEM_FINAL_PASS_GB (48) for step 4 -- from "
                         "real Core runs' sacct memory")
    ap.add_argument("--time-per-file", type=int, default=None,
                    help="hours per array task (steps 2/4). Default: TIME_PER_FILE_HOURS at "
                         "TIME_REFERENCE_CPUS CPUs per task, scaled up for fewer "
                         "(array_task_hours: 8 h at 8 CPUs)")
    ap.add_argument("--assembly-cpus", type=int, default=64)
    ap.add_argument("--assembly-mem", type=int, default=128)
    ap.add_argument("--assembly-time", type=int, default=12)
    ap.add_argument("--libpred-cpus", type=int, default=16)
    ap.add_argument("--libpred-mem", type=int, default=64)
    ap.add_argument("--libpred-time", type=int, default=4)
    ap.add_argument("--seed-lib", help="Skip step-1 prediction; use this existing library "
                    "(e.g. an InfinDIA empirical .parquet/.speclib) as the first-pass seed. "
                    "This is how you PARALLELIZE a semi-tryptic / non-specific / InfinDIA search: "
                    "build the small empirical library once (InfinDIA --pre-search), then fan the "
                    "per-file passes out across the cluster against it.")
    ap.add_argument("--seed-dep", help="SLURM job id the first pass should wait for (afterok) — "
                    "e.g. the InfinDIA lib-build job that produces --seed-lib.")
    ap.add_argument("--partition")
    ap.add_argument("--account")
    ap.add_argument("--qos", default=None, help="SLURM QOS. Needed for publicgrp/low "
                    "(publicgrp-low-qos); high/genome-center-grp uses its default.")
    ap.add_argument("--max-simultaneous", type=int, default=20)
    ap.add_argument("--no-norm", action="store_true")
    ap.add_argument("--probe-window", action=argparse.BooleanOptionalAction, default=True,
                    help="measure the scan-window radius after step 1 and pin it for all "
                         "steps (default: on). --no-probe-window requires --window in the cfg.")
    ap.add_argument("--allow-auto-mass-acc", action="store_true",
                    help="proceed even though the cfg OMITS mass accuracy (auto). Results "
                         "across steps become inconsistent -- only for deliberate testing. "
                         "It does not override an INVALID value (0, negative, non-numeric, "
                         "set twice), a bad --window, or an unparseable cfg.")
    ap.add_argument("--no-notify", action="store_true",
                    help="no Slack post from the chain's jobs (same as SKILL_SLACK=0); the "
                         "run log and the FRAN hand-over still happen. "
                         "See references/notifications.md")
    ap.add_argument("--no-fran", action="store_true",
                    help="step 5 never hands the search to FRAN (same as FRAN_DEPOSIT=off when "
                         "generating); step 7c can still stage it")
    ap.add_argument("--fran-name", help="the session's descriptive name, for stage --name")
    qcx = ap.add_mutually_exclusive_group()
    qcx.add_argument("--qc", action="store_true",
                     help="an instrument QC / standard run: never staged from a job")
    qcx.add_argument("--not-qc", action="store_true", help="not a QC run: passed on to stage")
    a = ap.parse_args()
    # The cfg path is recorded (search_provenance.json, params.resolved.cfg's origin) and read by
    # later jobs and readers in other working directories: absolute from here on, and it must
    # exist now -- a relative one is otherwise resolved against whatever directory reads it.
    if a.cfg:
        given, a.cfg = a.cfg, os.path.abspath(a.cfg)
        if not os.path.isfile(a.cfg):
            sys.exit(f"[diann_parallel] cfg not found: {a.cfg}"
                     + (f" (given as {given})" if given != a.cfg else "")
                     + ". Nothing was generated.")

    raws = []
    if a.raw_list:
        with open(a.raw_list) as fh:
            raws.extend(line.strip() for line in fh if line.strip())
    for p in a.raw:
        raws.extend(sorted(glob.glob(p)) or [p])
    raws = [os.path.abspath(r.rstrip("/")) for r in raws]
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from check_report_runs import distinct_inputs, names_stop, repeated_names, repeats_note
    # a file listed twice is searched once and flagged; different files sharing a run name are
    # flagged and stop the generator (check_report_runs: the one rule, as in run_search.py)
    raws, repeats = distinct_inputs(raws)
    shared = repeated_names(raws)
    if repeats or shared:
        sys.stderr.write(repeats_note(repeats, "diann_parallel", shared))

    # DIA-NN names a run -- and this generator names its .quant -- by the file name alone, so
    # /plate1/s1.mzML and /plate2/s1.mzML are ONE Run in the report and TWO array tasks writing
    # the same .quant. run_search.main() refuses that before it routes anywhere, but
    # diann_parallel.py is run directly too (this module's own docstring names that as the way
    # to size the chain by hand), and that route had no check at all: the chain would be
    # generated, burn its SLURM hours and merge two samples into one column. Knowable from the
    # input list, so it stops here.
    if shared:
        sys.exit(names_stop(shared, "diann_parallel"))
    refuse_unsafe_path(a.out)
    # run_search.py checks --dda against the bundle; the documented direct call -- and a
    # --seed-lib phase 2 -- never pass through it. So the generator checks it too, against the
    # acquisition estimate_params.py wrote the cfg for. An unreadable cfg is parallel_safe's to
    # report, below.
    if a.cfg:
        try:
            why = dda_mismatch(a.cfg, cfg_acquisition(a.cfg),
                               source=f"{a.cfg}.rationale.json's acquisition")
        except CfgError:
            why = None
        if why:
            sys.exit(f"[diann_parallel] REFUSING: {why} Nothing was generated.")

    # Detect the queue from the submitting user's SLURM associations rather than
    # assuming facility membership. genome-center-grp/high for members; publicgrp/low
    # for everyone else (incl. class accounts) — where `high` caps at 8 CPUs/job, so a
    # 32-CPU request would never start.
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from run_search import slurm_queue
        # a.qos goes in so that a partial queue (e.g. --qos alone) is completed from the
        # association it belongs to; the detected QOS itself is still not used (see below).
        a.partition, a.account, _q = slurm_queue(a.partition, a.account, a.qos)
        if not a.qos and a.partition == "low" and a.account == "publicgrp":
            a.qos = "publicgrp-low-qos"
    except Exception as e:
        # NEVER silently. This fallback puts the run on the PREEMPTIBLE queue, and for a facility
        # member that is a real demotion -- a multi-hour search exposed to preemption because an
        # import or sacctmgr call failed. Silent was how it read as "the skill just uses low".
        a.partition = a.partition or "low"
        a.account = a.account or "publicgrp"
        if not a.qos and a.partition == "low" and a.account == "publicgrp":
            a.qos = "publicgrp-low-qos"
        sys.stderr.write(f"[diann_parallel] WARNING: queue detection failed ({type(e).__name__}: "
                         f"{e}); falling back to {a.account}/{a.partition}, which is PREEMPTIBLE. "
                         f"If you are in genome-center-grp, pass --partition high --account "
                         f"genome-center-grp to avoid preemption.\n")
    n = len(raws)
    if n < 2:
        sys.exit("Parallel search needs >= 2 raw files (pass --raw or --raw-list); "
                 "use the single-shot run_search.py for 1.")
    out = os.path.abspath(a.out); os.makedirs(out, exist_ok=True)
    # Every path below is spliced into the job scripts in DOUBLE quotes (so an array task's $QUANT
    # still expands, and a folder with a space -- 2,557 of them in the Core's service tree --
    # stays one word): refused here if it holds what double quotes cannot carry.
    fasta = refuse_unsafe_path(os.path.abspath(a.fasta), "--fasta")
    dnet = dotnet_prefix(raws)             # .NET 8 export prefix when inputs are Thermo .raw
    DN = dnet + a.diann
    # The chain is only valid with pinned mass accuracy AND a scan window that is the same
    # in every step (steps 3/5 reuse .quant files; see mass_acc_status for DIA-NN's own
    # warning). parallel_safe()["ok"] is the shared rule -- run_search.py routes on this same
    # value, so the router and the generator cannot disagree. Gate on "ok", never on a
    # field of `ma`: gating on ma["fixed"] here while the router used "ok" is how a cfg that
    # one approved could be refused (or waved through) by the other.
    safe = parallel_safe(a.cfg, probe_window=a.probe_window, seed_lib=a.seed_lib)
    if not safe["ok"]:
        # The override is for the one deliberate-testing case it is named after: mass accuracy
        # left to DIA-NN. An invalid value is a mistake (0 is a literal 0 ppm tolerance), and a
        # bad --window or an unparseable cfg would be spliced into every step as-is.
        if safe["code"] in ("cfg_missing", "cfg_unparseable"):
            sys.exit(f"{safe['reason']}.\nFix: {safe['remedy']}.")
        if not (a.allow_auto_mass_acc and safe["code"] in MASS_ACC_OMITTED_CODES):
            sys.exit(
                f"Not parallel-safe: {a.cfg or '(no --cfg given)'} -- {safe['reason']}.\n"
                "The 5-step chain reuses .quant files across steps, so anything DIA-NN\n"
                "auto-optimises PER FILE (mass accuracy AND scan window) is applied\n"
                "inconsistently between passes and then stitched together.\n"
                f"Fix: {safe['remedy']}.\n"
                "Or run the single-shot search instead (run_search.py --no-parallel)."
                + ("\nTo override deliberately: --allow-auto-mass-acc."
                   if safe["code"] in MASS_ACC_OMITTED_CODES else ""))
        sys.stderr.write(f"[diann_parallel] WARNING: proceeding with auto mass accuracy "
                         f"({safe['reason']}) -- steps will not be mutually consistent.\n")
    win_probe, ma = safe["probe"], safe["ma"]
    # What step 1b measures ("window", "mass-acc"; empty when it does not run), and the
    # mass-accuracy levels it pins as documented rather than measured ({flag: ppm}).
    measure, documented = safe["measure"], safe["mass_acc_documented"]

    # What step 1b measures is PREFIXED onto steps 2-5 at run time, so the same flags still in
    # the cfg have to come out or two values land on the same command line. (A planned mass
    # accuracy has neither flag in the cfg -- mass_acc_measure_plan requires it -- so dropping
    # them only makes that explicit.)
    measured_flags = ((("--window",) if "window" in measure else ())
                      + (MASS_ACC_FLAGS if "mass-acc" in measure else ()))
    flags = read_cfg_flags(a.cfg, drop=measured_flags)
    D = out  # all DIA-NN intermediate/output lives here (real paths; native binary reads them directly)
    from check_report_runs import stats_path   # one definition of <report>.stats.tsv
    report = NO_NORM_REPORT if a.no_norm else "report.parquet"
    norm = "--no-norm" if a.no_norm else ""
    xic = xic_flag(a.cfg)          # step 4 only -- see xic_flag() docstring

    # Every refusal is behind us: an earlier search's probe outputs must not describe this one
    # (set aside, not deleted -- they are that search's evidence)
    set_aside_probe_outputs(out)
    # file list (1 raw path per line) — array tasks index into it
    open(os.path.join(out, "file_list.txt"), "w").write("\n".join(raws) + "\n")
    all_f = " ".join(f"--f {shlex.quote(r)}" for r in raws)   # quote — data paths may contain spaces
    array = f"0-{n-1}%{a.max_simultaneous}"
    seed = refuse_unsafe_path(os.path.abspath(a.seed_lib), "--seed-lib") if a.seed_lib else None
    predicted = seed if seed else f"{D}/step1.predicted.speclib"
    empirical = f"{D}/empirical.parquet"
    first_pass = f"{D}/{FIRST_PASS_REPORT}"
    pass_cmp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "pass_comparison.py")

    # QUEUE CHOICE, ported from DE-LIMP's select_best_partition()
    # (R/helpers_search.R) + docs/QUEUE_SWITCHING.md:
    #   steps 2 & 4 are ARRAYS -- embarrassingly parallel, one file per task, so a
    #     preemption costs one task. They are also exactly what the per-user CPU cap on
    #     the priority queue throttles, so prefer publicgrp/low when it has idle CPUs
    #     (it routinely has thousands).
    #   steps 1, 3 & 5 are SINGLE jobs that cannot restart mid-way (3 consumes every
    #     .quant file, 5 needs all of step 4), so a preemption throws the stage away --
    #     keep them on the priority queue.
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from run_search import slurm_queue as _sq
        _pa, _aa, _qa = _sq(a.partition, a.account, a.qos,
                            peak_cpus=a.threads_per_file or a.threads_max, preemptible_ok=True)
        _ps, _as_, _qs = _sq(a.partition, a.account, a.qos, peak_cpus=a.assembly_cpus)
    except Exception:
        _pa, _aa, _qa = a.partition, a.account, a.qos
        _ps, _as_, _qs = a.partition, a.account, a.qos
    cpu_sizing, probe_cpus = size_cpus(a, n, (_pa, _aa, _qa), (_ps, _as_, _qs))

    # The job-end hook (notify_slack.wrap_job_script, the one definition: run log -> FRAN ->
    # Slack). Step 5 ends the search: it reports success and failure, and stages a finished
    # search for FRAN. Every earlier step reports only a failure: the rest of the chain waits on
    # it with afterok and never starts, so step 5 would never get to say anything. An array
    # reports its first failing task only.
    import notify_slack
    # Every job first checks that its node reaches the engine, the FASTA, the raw data and this
    # folder (node_fault.py): a node whose /quobyte mount has dropped fails as a NODE fault (exit
    # 75), which watch_run.sh and node_fault.py retry act on, never as "DIA-NN exited 0 but ...".
    from node_fault import preflight_lines, chain_checks
    node_check = preflight_lines(chain_checks(a.diann, out, fasta, raws)
                                 + ([("-r", seed, "the seed library")]
                                    if seed and os.path.exists(seed) else []))

    def write(name, body, stage=None, hours=None, final=False):
        if stage:
            # final = step 5, whose quant count is the chain's completeness guard: the one step
            # that may stage for FRAN
            body = notify_slack.wrap_job_script(body, out, final=final, time_limit_h=hours,
                                                stage=stage, slack=not a.no_notify,
                                                fran=not a.no_fran, fran_guarded=final,
                                                fran_name=a.fran_name,
                                                qc=(True if a.qc else
                                                    False if a.not_qc else None),
                                                preflight=node_check).rstrip("\n")
        p = os.path.join(out, name)
        open(p, "w").write(body + "\n"); os.chmod(p, 0o755); return name

    # array preamble: pick this task's raw file
    pick = ('FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ' + f'"{D}/file_list.txt")\n'
            'if [ -z "$FILE" ]; then echo "no file for task $SLURM_ARRAY_TASK_ID"; exit 1; fi\n'
            'echo "Processing: $FILE"\n')

    # Step 1 — library prediction (single job) — SKIPPED when --seed-lib is given
    s1 = None
    if not seed:
        s1 = write("step1_libpred.sbatch", "\n".join([
            header("s1_libpred", a.libpred_cpus, a.libpred_mem, a.libpred_time, _ps, _as_, qos=_qs), "",
            f'echo "Step 1/5 library prediction"; date',
            clear_stale(predicted),       # see clear_stale(): a re-run must not pass on the old one
            f'{DN} --fasta "{fasta}" --fasta-search --predictor --gen-spec-lib \\',
            f'  --out-lib "{D}/step1.speclib" --out "{D}/step1_lib.parquet" \\',
            f'  --threads {a.libpred_cpus} {flags}',
            must_exist(predicted, "the predicted spectral library")]),
            stage="step 1/5 library prediction", hours=a.libpred_time)

    # Step 1b — measure the scan-window radius (and, when planned, mass accuracy) ONCE, so
    # steps 2-5 share it.
    # DIA-NN optimises the radius per file when --window is absent, and steps 3/5 then
    # combine .quant files produced under different windows -- which DIA-NN's own
    # warning calls "strongly not recommended". Measuring beats guessing: the radius
    # depends on the acquisition scheme (cycle time vs peak width), not the instrument.
    #
    # MASS ACCURACY ("mass-acc" in measure): only for a cfg estimate_params.py planned that way
    # -- an Orbitrap with a level outside DIA-NN's table (see parallel_safe). The same DIA-NN run
    # per probe logs it: measured on HIVE with DIA-NN 2.7.0 (srun job 23522741, Exploris 480
    # 120k/15k, full mouse library) the radius came at 1:34, "Recommended MS1 mass accuracy
    # setting: 4.1 ppm" at 1:35 and "Optimised mass accuracy: 14 ppm" at 2:27, of a 5:13 search.
    # The probes run with BOTH flags omitted -- DIA-NN 2.7.0 fixes both levels when either is
    # given -- and massacc.txt gets the two flags: the median of a measured level, the
    # documented value of a level that has one (--ms1-ppm 7 at 120k; a measured 4.2 there drew
    # DIA-NN's "deviates significantly from the value recommended (7 ppm)" warning on every
    # pass). Both also go into params.resolved.cfg.
    s1b = None
    resolved_cfg = a.cfg
    what = " and ".join({"window": "scan-window radius", "mass-acc": "mass accuracy"}[m]
                        for m in measure)
    wtxt, mtxt = f"{D}/window.txt", f"{D}/massacc.txt"
    if win_probe:
        probe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "probe_window.py")
        # The measured radius has to end up in a PARAMETER FILE, not just window.txt, or the
        # run is not reproducible from what we recorded: search_provenance.json would name a
        # cfg with no --window, and replaying it would re-optimise per file and not reproduce
        # the numbers (SKILL.md golden rule 5). params.base.cfg is the cfg minus any --window,
        # by the same token rule as the step flags. params.resolved.cfg is NOT created here:
        # step 1b builds it in a .tmp and moves it into place only once a radius is measured,
        # so a "resolved" cfg with no --window can never exist to be replayed -- unless the probe
        # measured nothing and fell back (probe_attempts()): then it has no --window BECAUSE
        # steps 2-5 ran without one, and replaying it lets DIA-NN choose per run, as they did.
        base_cfg = os.path.join(out, "params.base.cfg")
        resolved_cfg = os.path.join(out, "params.resolved.cfg")
        tmp_cfg = resolved_cfg + ".tmp"
        probe_dir = os.path.join(out, "window_probe")
        write_cfg(a.cfg, base_cfg, drop=measured_flags)
        q = shlex.quote
        mess = "mass-acc" in measure
        # The window-only script is exactly upstream's; mass accuracy adds its own lines.
        doc_args = f" {probe_mass_acc_args(documented)}" if documented else ""
        failed_if = " || ".join(f'[ -z "${v}" ]' for m, v in (("window", "W"), ("mass-acc", "M"))
                                if m in measure)
        s1b = write("step1b_window.sbatch", "\n".join([
            header("s1b_window", probe_cpus, task_memory(a)[0], PROBE_WALL_HOURS,
                   _ps, _as_, qos=_qs), "",
            f'echo "Step 1b/5 measuring {what} on representative runs"; date',
            # Every other step reaches DIA-NN through DN, which carries the .NET 8 exports a
            # Thermo .raw needs. probe_window.py runs DIA-NN as its own subprocess (no shell),
            # so the prefix cannot ride on --diann -- it has to be in the ENVIRONMENT the probe
            # inherits. Without it DIA-NN cannot read .raw, no radius is logged, window.txt is
            # never written, and steps 2-5 sit on afterok for ever.
            *([f"{dnet.strip()}   # .NET 8 for Thermo .raw, inherited by probe_window.py's DIA-NN"]
              if dnet else []),
            # A resubmitted step 1b must never find the previous run's answer and carry on.
            f'rm -f "{wtxt}" "{D}/window.json" "{D}/window.json.attempt1" "{D}/{FALLBACK_RECORD}" '
            f"{q(resolved_cfg)} {q(tmp_cfg)}" + (f' "{mtxt}"' if mess else ""),
            f"rm -rf {q(probe_dir)} {q(probe_dir + '.attempt1')}",
            f"cp {q(base_cfg)} {q(tmp_cfg)}",
            # WHICH runs: not the first files of the listing. The probe gets the whole cohort
            # (file_list.txt) and, at run time on this node, keeps out blanks, washes, failed
            # injections and any .d whose analysis.tdf index is damaged (WAL mode, a stale
            # -wal/-journal, or an index that stops short of tdf_bin -- 342 such .d on HIVE, and
            # a bytes rule picked one in the pilot). It measures the median and quartile runs,
            # replaces a run that does not log everything asked with the next one nearest the
            # median, and pins the MEDIAN. window.json records every probe, and the damaged runs.
            # --timeout bounds one hung run; --budget bounds them all inside this job's limit.
            *(['W=""'] if "window" in measure else []),
            *(['M=""'] if mess else []),
            # the probe, once more if it measured nothing, then probe_fallback.py -- see
            # probe_attempts(): a failed probe no longer blocks steps 2-5
            *probe_attempts(
                "python3 " + " \\\n    ".join([
                    f"{q(probe)} --diann {q(a.diann)}",
                    f"--raw-list {q(os.path.join(out, 'file_list.txt'))}",
                    f"--fasta {q(fasta)} --lib {q(predicted)} --threads {probe_cpus}",
                    f"--max-probes {PROBE_CANDIDATES} --max-failures {PROBE_MAX_FAILURES}",
                    f"--timeout {PROBE_TIMEOUT_S}"
                    + (f" --measure {' '.join(measure)}{doc_args}" if mess else ""),
                    f"--workdir {q(probe_dir)} --write-cfg {q(tmp_cfg)}"]),
                # the flags as bash words, after `--`: the probe's DIA-NN gets the same argv as
                # steps 2-5, not a second parse of them through shlex
                flags, f"{D}/window.json", measure, documented,
                reset=keep_attempt(probe_dir) + [f"cp {q(base_cfg)} {q(tmp_cfg)}"],
                window_file=wtxt if "window" in measure else None,
                massacc_file=mtxt if mess else None, write_cfg=tmp_cfg,
                provenance=f"{D}/search_provenance.json",
                fallback_out=f"{D}/{FALLBACK_RECORD}"),
            'if [ "$PROBE_RC" -eq 0 ]; then',
            *(['  W=$(python3 -c "import json,sys; w=json.load(open(sys.argv[1]))[\'window_radius\']; '
               f'assert isinstance(w, int) and w > 0; print(w)" "{D}/window.json")']
              if "window" in measure else []),
            # the two flags, in exactly the shape steps 2-5's guard (needs_measured) accepts
            *(['  M=$(python3 -c "import json,re,sys; m=json.load(open(sys.argv[1]))[\'mass_acc\'][\'pin_as\']; '
               f'assert re.fullmatch(sys.argv[2], m); print(m)" "{D}/window.json" '
               f'{q(MEASURED_FILE_RE["massacc.txt"])})'] if mess else []),
            'elif [ "$PROBE_FALLBACK" -eq 1 ]; then',
            # probe_fallback.py wrote window.txt / massacc.txt and the cfg's mass-acc flags
            *([f'  W=$(cat "{wtxt}")'] if "window" in measure else []),
            *([f'  M=$(cat "{mtxt}")'] if mess else []),
            "fi",
            f"if {failed_if}; then",
            f'  echo "FAILED: step 1b measured no {what} -- the reason and each probe\'s DIA-NN '
            'log tail are above -- and did not fall back (only a failure of the probe\'s own '
            'machinery may):" >&2',
            *probe_failure_lines(),
            f'  echo "Every probe is recorded in {D}/window.json." >&2',
            # Resubmitting step 1b ALONE does not restart the chain: steps 2-5 were submitted
            # afterok on THIS job id, so they sit PENDING (DependencyNeverSatisfied) for ever.
            '  echo "Steps 2-5 were submitted afterok on THIS job, so they are now PENDING with '
            'DependencyNeverSatisfied and will never start -- even if step 1b is resubmitted '
            'and succeeds." >&2',
            '  echo "Recover: fix the cause (do NOT guess a '
            + ("--window or a mass accuracy" if mess and "window" in measure else
               "mass accuracy" if mess else "--window")
            + f'), scancel steps 2-5 (ids in '
            f'{D}/jobs.txt), then resubmit step1b_window.sbatch and steps 2-5 chained afterok on '
            'the new ids, reusing step1.predicted.speclib -- or re-run submit.sh, which also '
            'repeats step 1. See references/watcher.md (dependency_failed)." >&2',
            # window.json stays: it is the evidence. window.txt, massacc.txt and the resolved cfg
            # never exist.
            f'  rm -f {q(tmp_cfg)} "{wtxt}"' + (f' "{mtxt}"' if mess else ""),
            "  exit 1",
            "fi",
            # These are written by this script, not by DIA-NN, so must_exist()'s "DIA-NN exited
            # 0 but did not write" would name the wrong culprit. Say what actually failed.
            *([f'if ! echo "$W" > "{wtxt}" || [ ! -f "{wtxt}" ] || [ ! -s "{wtxt}" ]; then',
               f'  echo "FAILED: radius $W was measured but could not be written to {wtxt} '
               '(disk full? permissions?)" >&2',
               "  exit 1",
               "fi"] if "window" in measure else []),
            *([f'if ! echo "$M" > "{mtxt}" || [ ! -f "{mtxt}" ] || [ ! -s "{mtxt}" ]; then',
               f'  echo "FAILED: mass accuracy $M was measured but could not be written to {mtxt} '
               '(disk full? permissions?)" >&2',
               "  exit 1",
               "fi"] if mess else []),
            # -f as well as -s: `mv` INTO a directory of that name succeeds, and a directory
            # is non-empty.
            f"if ! mv -f {q(tmp_cfg)} {q(resolved_cfg)} || [ ! -f {q(resolved_cfg)} ] "
            f"|| [ ! -s {q(resolved_cfg)} ]; then",
            ("  echo \"FAILED: radius $W was measured" if measure == ["window"] else
             f"  echo \"FAILED: {what} {'were' if len(measure) > 1 else 'was'} measured")
            + f' but {resolved_cfg} could not be moved into '
            f'place from {tmp_cfg} (disk full? permissions?)" >&2',
            "  exit 1",
            "fi",
            'if [ "$PROBE_FALLBACK" -eq 1 ]; then',
            *(['  echo "scan window = auto -- FALLBACK: the probe measured nothing, so DIA-NN '
               'chooses the radius itself, per run (search_provenance.json probe_fallback)"']
              if "window" in measure else []),
            *(['  echo "mass accuracy = $M -- FALLBACK: the documented level as given, the other '
               'at the facility SOP (DEFAULT, not measured)"'] if mess else []),
            'else',
            *(['  echo "scan window radius = $W (median of the runs listed above; pinned for '
               'steps 2-5)"'] if "window" in measure else []),
            *(['  echo "mass accuracy = $M (per level: the median of the runs listed above, or the '
               'documented value; pinned for steps 2-5)"'] if mess else []),
            'fi',
            f'echo "fully-resolved parameters -> {resolved_cfg}"']),
            stage=f"step 1b/5 measuring {what}", hours=PROBE_WALL_HOURS)
    # steps 2-5 read the measured values at RUNTIME so every pass uses the identical ones --
    # after checking they are there (needs_measured: a missing file expands to nothing)
    wflag = window_flag(wtxt) if "window" in measure else ''
    mflag = f'$(cat "{mtxt}") ' if "mass-acc" in measure else ''
    measured_guard = (([needs_measured(wtxt, "a scan-window radius")]
                       if "window" in measure else [])
                      + ([needs_measured(mtxt, "the two mass-accuracy flags")]
                         if "mass-acc" in measure else []))

    # DIA-NN aborts with "cannot find the temp folder" if --temp does not exist -- it will NOT
    # create it -- and it does so BEFORE doing any work, so the whole submission cycle is lost to
    # a missing directory. submit.sh makes them, but the watcher playbook (references/watcher.md)
    # tells the orchestrator to resubmit individual steps after a failure, and `sbatch
    # step4_finalpass.sbatch` never goes through submit.sh. So each step makes its own: mkdir -p
    # is idempotent and free, and it means no step can be submitted into that error.
    def tmpguard(d):
        return f'mkdir -p "{D}/{d}"   # DIA-NN will NOT create --temp and aborts without it'

    # Step 2 — first pass (array): predicted lib -> per-file .quant
    s2 = write("step2_firstpass.sbatch", "\n".join([
        header("s2_firstpass", a.threads_per_file, task_memory(a)[0], a.time_per_file, _pa, _aa, qos=_qa, array=array), "",
        f'echo "Step 2/5 first pass, task ${{SLURM_ARRAY_TASK_ID}} of {n}"; date',
        *measured_guard, pick, tmpguard("quant_step2"),
        f'mkdir -p "{D}/{TASK_OUT_DIRS["step2"]}"   # this task\'s report files (TASK_OUT_DIRS)',
        'QOUT="${FILE##*/}"; QOUT="${QOUT%.*}.quant"',
        clear_stale(f'{D}/quant_step2/$QOUT'),
        f'{DN} --f "$FILE" --fasta "{fasta}" --lib "{predicted}" \\',
        f'  --temp "{D}/quant_step2" --rt-profiling --gen-spec-lib --quant-ori-names \\',
        f'  --out "{D}/{TASK_OUT_DIRS["step2"]}/t${{SLURM_ARRAY_TASK_ID}}.parquet" \\',
        f'  --threads {a.threads_per_file} {wflag}{mflag}{flags}',
        must_exist(f'{D}/quant_step2/$QOUT', "this file's .quant")]),
        stage="step 2/5 first pass (array)", hours=a.time_per_file)

    # Step 3 — empirical library assembly (single job, --use-quant)
    s3 = write("step3_assembly.sbatch", "\n".join([
        header("s3_assembly", a.assembly_cpus, a.assembly_mem, a.assembly_time, _ps, _as_, qos=_qs), "",
        f'echo "Step 3/5 empirical library assembly"; date', *measured_guard,
        tmpguard("quant_step2"),
        f'cp -r "{D}/quant_step2" "{D}/quant_step2_orig" 2>/dev/null || true   # backup for resume',
        # the first-pass report and its stats too: step 5 compares them with its own
        # (pass_comparison.py), and a previous run's must not stand in for this one's
        clear_stale(empirical, first_pass, stats_path(first_pass)),
        f'{DN} {all_f} --fasta "{fasta}" --lib "{predicted}" --use-quant --quant-ori-names \\',
        f'  --rt-profiling --gen-spec-lib --out-lib "{empirical}" \\',
        f'  --temp "{D}/quant_step2" --out "{first_pass}" \\',
        f'  --threads {a.assembly_cpus} {wflag}{mflag}{flags}',
        must_exist(empirical, "the empirical spectral library")]),
        stage="step 3/5 empirical library assembly", hours=a.assembly_time)

    # Step 4 — final pass (array): empirical lib -> per-file .quant
    # --xic alone is not enough. DIA-NN names the XIC folder after the --out report
    # basename; with no --out it resolves that against the filesystem ROOT and dies with
    # "cannot create directory: Permission denied [/report_xic]" -- AFTER writing the
    # .quant, and STILL EXITING 0. Verified on DIA-NN 2.6.0 (single file, --xic, no
    # --out): zero .xic.parquet produced, SLURM records COMPLETED, afterok advances.
    # So step 4 needs its own --out whenever XICs are requested. DIA-NN then writes
    # <out>/xic/t<TASKID>_xic/<run>.xic.parquet, one folder per array task.
    # Both carry their own leading space so the command line has no double space (and no
    # trailing space before the line continuation) when XICs are off.
    xic_arg = f' {xic}' if xic else ''
    # With XICs off the task still gets its own --out (TASK_OUT_DIRS), for the same reason as
    # step 2: otherwise its report files land in <out> itself.
    out4_dir = "xic" if xic else TASK_OUT_DIRS["step4"]
    task_out4 = f' --out "{D}/{out4_dir}/t${{SLURM_ARRAY_TASK_ID}}.parquet"'

    s4 = write("step4_finalpass.sbatch", "\n".join([
        header("s4_finalpass", a.threads_per_file, task_memory(a)[1], a.time_per_file, _pa, _aa, qos=_qa, array=array), "",
        f'echo "Step 4/5 final pass, task ${{SLURM_ARRAY_TASK_ID}} of {n}"; date',
        *measured_guard, pick, tmpguard("quant_step4"),
        'QUANT="${FILE##*/}"; QUANT="${QUANT%.*}.quant"',
        # cleared BEFORE the skip: a skipped task must leave no .quant for step 5 to count
        clear_stale(f'{D}/quant_step4/$QUANT'),
        f'if [ ! -f "{D}/quant_step2/$QUANT" ]; then echo "SKIP: no step-2 quant for $QUANT"; exit 0; fi',
        f'mkdir -p "{D}/{out4_dir}"',
        # DIA-NN re-saves the library it is handed as "<lib>.skyline.speclib", written
        # NEXT TO --lib. With one shared path every concurrent array task writes the same
        # file; most win the race in seconds, the losers block until the wall clock kills
        # them. Give each task its own copy so there is nothing to contend on.
        f'LIBPRIV="{D}/libpriv/t${{SLURM_ARRAY_TASK_ID}}"',
        'mkdir -p "$LIBPRIV"',
        f'cp -f "{empirical}" "$LIBPRIV/lib.parquet"',
        'trap \'rm -rf "$LIBPRIV"\' EXIT',
        f'{DN} --f "$FILE" --fasta "{fasta}" --lib "$LIBPRIV/lib.parquet" \\',
        f'  --temp "{D}/quant_step4" --quant-ori-names{xic_arg}{task_out4} \\',
        f'  --threads {a.threads_per_file} {wflag}{mflag}{flags}',
        must_exist(f'{D}/quant_step4/$QUANT', "this file's final-pass .quant")]),
        stage="step 4/5 final pass (array)", hours=a.time_per_file)

    # Step 5 — cross-run report (single job, --use-quant --matrices)
    s5 = write("step5_report.sbatch", "\n".join([
        header("s5_report", a.assembly_cpus, a.assembly_mem, a.assembly_time, _ps, _as_, qos=_qs), "",
        f'echo "Step 5/5 cross-run report"; date', *measured_guard, tmpguard("quant_step4"),
        # the report AND its stats file: check_report_runs.stats_path() names the latter
        clear_stale(f'{D}/{report}', stats_path(f'{D}/{report}')),
        # A step-4 task that failed silently leaves no .quant, and step 5 happily
        # reports on whatever survived. Count them: fewer quants than inputs means a
        # sample was dropped, which must never pass as success. BEFORE DIA-NN runs: counted
        # after, a report from N-1 runs was already on disk when the job failed (review of
        # 2.10: a node_fault retry that left out an OOM task), at the path run_de.R reads.
        # step 5 is afterok on step 4, so this fires only on a chain resubmitted by hand or
        # by a retry -- or a step-4 task that skipped its file.
        #
        # Count FILES, not input lines. This chain's names come from file_list.txt, not from
        # `ls quant_step4/*.quant` -- a previous search's .quant for another run would
        # otherwise make up the number for a missing one -- but a name derived per INPUT LINE
        # re-introduces the very bug: two inputs whose basenames collide (/plate1/s1.mzML,
        # /plate2/s1.mzML) map to ONE s1.quant, and incrementing once per line counts that
        # single surviving file twice. `NQ=2, n=2 -> PASS` for a report in which two samples
        # were merged into one Run. `sort -u` collapses them to the files that can actually
        # exist, so the count falls short and the job fails -- which is what `ls | wc -l` gave
        # before. main() refuses colliding names outright; this stays the backstop for a chain
        # generated before that check, or whose file_list.txt was edited by hand.
        f'QUANTS=$(while IFS= read -r f; do [ -n "$f" ] || continue; b="${{f##*/}}"; '
        f'printf "%s\\n" "${{b%.*}}.quant"; done < "{D}/file_list.txt" | sort -u)',
        'NQ=0; NNAMES=0',
        'while IFS= read -r q; do [ -n "$q" ] || continue; NNAMES=$((NNAMES + 1)); '
        f'if [ -s "{D}/quant_step4/$q" ]; then NQ=$((NQ + 1)); '
        f'else echo "MISSING: $q -- no step-4 task wrote it" >&2; fi; done <<< "$QUANTS"',
        f'if [ "$NNAMES" -ne {n} ]; then '
        f'echo "FAILED: {n} inputs map to only $NNAMES distinct run names -- DIA-NN names a run '
        f'by its file name alone, so inputs that share one are merged into a single Run. '
        f'Rename them, or search them separately." >&2; exit 1; fi',
        f'if [ "$NQ" -ne {n} ]; then '
        f'echo "FAILED: only $NQ of {n} runs have a final-pass .quant -- a step-4 task did not '
        f'succeed, so no report was written. Resubmit the step-4 tasks of the runs listed '
        f'MISSING above, then this step." >&2; '
        f'exit 1; fi',
        f'{DN} {all_f} --fasta "{fasta}" --lib "{empirical}" --use-quant --quant-ori-names \\',
        f'  --temp "{D}/quant_step4" --matrices --out "{D}/{report}" \\',
        f'  --threads {a.assembly_cpus} {norm} {wflag}{mflag}{flags}',
        must_exist(f'{D}/{report}', "the cross-run report"),
        f'echo "OK: report built from all {n} runs"',
        # Did the final pass keep what the first pass found? A table, a WARNING per run that
        # lost more than half its precursors, and the record in search_provenance.json -- see
        # pass_comparison.py. It never fails the job: the report is not wrong, the numbers are a
        # decision for the user. A comparison that cannot be made is said, never skipped quietly.
        # Under the generator's own python (the pipeline env's, which has pyarrow for the rows
        # after run_de.R's q filter), else python3 -- which still gives the stats table.
        f'PY={shlex.quote(sys.executable)}; [ -x "$PY" ] || PY=python3; '
        f'"$PY" {shlex.quote(pass_cmp)} --first-report {shlex.quote(first_pass)} '
        f'--final-report {shlex.quote(D + "/" + report)} '
        f'--out {shlex.quote(D + "/pass_comparison.json")} '
        f'--provenance {shlex.quote(D + "/search_provenance.json")} '
        '|| echo "WARNING: the first-pass / final-pass comparison could not be made -- the '
        'messages above say why" >&2']),
        stage="step 5/5 cross-run report", hours=a.assembly_time, final=True)

    # submit.sh — chain the steps with afterok dependencies.
    #
    # It runs under `set -euo pipefail`, and the orchestrator often reads only the first line of
    # what it prints (`bash submit.sh | head -1`). The reader then closes the pipe, and the next
    # echo killed bash with SIGPIPE -- after every sbatch had succeeded, but before the checkpoint
    # and jobs.txt were written (msalemi, SET28 2026-09-25: a whole chain queued with no jobs.txt
    # and no RECOVERY.md). So SIGPIPE is ignored, every message goes through say(), which cannot
    # fail, and jobs.txt and the checkpoint are written straight after the last sbatch, before
    # anything else is printed -- the order run_search.py's two-job submit.sh already uses.
    sub_lines = ["#!/bin/bash", "set -euo pipefail",
                 "trap '' PIPE   # a reader that stops early must not kill this script",
                 # messages only: it never fails, so a closed stdout cannot stop the script
                 "say() { printf '%s\\n' \"$*\" 2>/dev/null || true; }",
                 f'cd "{out}"',
                 f'mkdir -p "{D}/quant_step2" "{D}/quant_step4"   # DIA-NN --temp dirs MUST pre-exist']
    if seed:
        # Step 1 skipped — the InfinDIA/empirical seed library IS the first-pass lib.
        # First pass optionally waits (afterok) on the lib-build job that produces it.
        dep2 = f'--dependency=afterok:{a.seed_dep} ' if a.seed_dep else ''
        sub_lines += [
            f'say "Step 1/5 SKIPPED — seeding first pass with {predicted}"',
            'jid2=$(sbatch --parsable %s%s)' % (dep2, s2)]
    else:
        sub_lines += ['jid1=$(sbatch --parsable %s)' % s1]
        if s1b:
            sub_lines += [
                'jid1b=$(sbatch --parsable --dependency=afterok:$jid1 %s)' % s1b,
                'jid2=$(sbatch --parsable --dependency=afterok:$jid1b %s)' % s2]
        else:
            sub_lines += ['jid2=$(sbatch --parsable --dependency=afterok:$jid1 %s)' % s2]
    sub_lines += [
        'jid3=$(sbatch --parsable --dependency=afterok:$jid2 %s)' % s3,
        'jid4=$(sbatch --parsable --dependency=afterok:$jid3 %s)' % s4,
        'jid5=$(sbatch --parsable --dependency=afterok:$jid4 %s)' % s5]

    # Record the submission so a disconnected session can pick the run back up.
    # SLURM keeps the jobs alive after the user closes their terminal; without this
    # the next session has no idea what was submitted or what is still outstanding.
    # <out> is normally <session>/output/search, so the session is two levels up.
    sess = os.path.dirname(os.path.dirname(out))
    ck = os.path.join(os.path.dirname(os.path.abspath(__file__)), "checkpoint.py")
    # $jid1b belongs here too. It is an afterok link like any other, so if it dies the whole
    # chain stalls forever on a dependency that will never fire -- and left out of jobs.txt
    # that shows up to `watch_run.sh --all` as "nothing running, nothing failed", with a
    # resuming session having no record the job was ever submitted (golden rule 7b).
    jobs = ('$jid2,$jid3,$jid4,$jid5' if seed else
            '$jid1,$jid1b,$jid2,$jid3,$jid4,$jid5' if s1b else
            '$jid1,$jid2,$jid3,$jid4,$jid5')
    sub_lines += [
        # the record first: nothing may be printed between the last sbatch and these two
        f'printf "%s\\n" {jobs.replace(",", " ")} > "{D}/jobs.txt"',
        f'python3 "{ck}" record --session "{sess}" --stage search \\',
        f'  --jobs "{jobs}" --desc "DIA-NN 5-step parallel chain ({n} files)" \\',
        f'  --report "{D}/{report}" --watch-job "$jid5" --watch-log "{D}/s5_report_${{jid5}}.log" \\',
        # what follows the chain is step 8's normalisation check, then the final DE held to it --
        # never the final DE alone (normalization_check.step8_commands: the one wording)
        f'  --next {shlex.quote(step8_next(f"{D}/{report}", sess))} \\',
        '  >/dev/null 2>&1 || true',
        'say "submitted: firstpass=$jid2 assembly=$jid3 finalpass=$jid4 report=$jid5"',
        f'say "final report will be {D}/{report}; watch with: watch_run.sh --slurm $jid5 --log {D}/s5_report_${{jid5}}.log"',
        f'say "all chain job ids -> {D}/jobs.txt  (watch the WHOLE chain: watch_run.sh --all {D})"',
        f'say "recovery notes written to {sess}/RECOVERY.md — you can safely close your terminal"']
    write("submit.sh", "\n".join(sub_lines))

    # Describe what WILL run, not what the cfg says (CLAUDE.md rule 1). On the probe path the
    # radius and the resolved cfg do not exist yet -- they are produced at run time by step
    # 1b -- so they are recorded as such, not as if they were already resolved.
    if win_probe:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from probe_window import SELECTION_RULE
    if "mass-acc" in measure:
        # Upstream's mass_acc_record() describes the CFG, which omits both flags -- "not in the
        # cfg (DIA-NN calibrates it itself)". That is not what runs: step 1b measures it and every
        # step gets the same pinned value. The documented level is known now; the measured one
        # only at run time, in massacc.txt.
        mass_acc = {"mode": "measured", "fixed": True, "measured": True,
                    "documented": documented,
                    "ms1": documented.get("--mass-acc-ms1"), "ms2": documented.get("--mass-acc"),
                    "source": "measured at run time by step 1b (probe_window.py --measure "
                              "mass-acc) and pinned for steps 2-5; planned by estimate_params.py "
                              "(measure_with_diann) for an Orbitrap level with no documented "
                              "DIA-NN value" + (f"; {documented_levels_text(documented)}"
                                                if documented else ""),
                    "value_file": mtxt, "evidence_file": f"{D}/window.json",
                    "probe_rule": SELECTION_RULE,
                    "sop_floor": dict(SOP_MASS_ACC_FLAGS),
                    "floor_note": floor_note(f"{D}/window.json", mtxt),
                    "reason": "not in the cfg; measured with DIA-NN on representative runs in "
                              "step 1b, the same value for every step",
                    # a probe that measures nothing twice falls back instead of failing the chain;
                    # step 1b's probe_fallback.py then REPLACES this record (result.mass_acc)
                    "if_probe_fails": FALLBACK_PLAN["mass-acc"]}
    else:
        mass_acc = mass_acc_record(ma, a.cfg)
    if "window" in measure:
        scan_window = {"mode": "measured",
                       "source": "measured at run time by step 1b (probe_window.py) and pinned "
                                 "for steps 2-5",
                       "value": None, "value_file": wtxt,
                       # which runs were measured is decided on the compute node, against the
                       # files as they are then; window.json is the record of it
                       "evidence_file": f"{D}/window.json",
                       "probe_rule": SELECTION_RULE,
                       # ...and replaced by probe_fallback.py's when the probe measures nothing
                       "if_probe_fails": FALLBACK_PLAN["window"]}
    else:
        # --window pinned in the cfg (also when step 1b measures only mass accuracy), or, under
        # --allow-auto-mass-acc, whatever the cfg hands every step
        scan_window = dict(window_record(ma), value_file=None)
    if win_probe:
        resolved = {"file": resolved_cfg, "produced": "runtime", "by": s1b,
                    "note": f"written by step 1b only after the {what} "
                            f"{'are' if len(measure) > 1 else 'is'} measured -- or set by the "
                            "fallback when the probe measured nothing (probe_fallback.json: "
                            "no --window in it, DIA-NN chooses per run); absent until then, so "
                            "a missing file after step 1b means the probe was refused or "
                            "stopped"}
    elif safe["code"] == "dda_window_unset":
        resolved = {"file": resolved_cfg, "produced": "generation", "by": None,
                    "note": "the cfg as given pins mass accuracy; --window is not in it, and for "
                            "a DDA search nothing measures it (see scan_window)"}
    elif safe["ok"]:
        resolved = {"file": resolved_cfg, "produced": "generation", "by": None,
                    "note": "the cfg as given already pins mass accuracy and --window"}
    else:
        # --allow-auto-mass-acc. Describe what the steps are actually handed -- a `--window 7`
        # in the cfg IS passed to every step -- rather than assuming the override unpinned it.
        resolved = {"file": resolved_cfg, "produced": "generation", "by": None,
                    "note": "NOT fully resolved: mass accuracy is not in the cfg "
                            "(--allow-auto-mass-acc), so DIA-NN chooses it at run time and "
                            "this cfg does not record the value used"}

    import json
    print(json.dumps({
        "out": out, "n_files": n, "report": f"{D}/{report}",
        "parallel_safe": {k: safe[k] for k in ("ok", "probe", "code", "reason")},
        # what step 1b measures at run time ("window", "mass-acc"); [] when nothing is measured
        "step1b_measures": measure,
        "mass_acc": mass_acc,
        "scan_window": scan_window,
        "resolved_params": resolved,
        # made at run time by step 5 (pass_comparison.py), which also copies it into
        # search_provenance.json as `pass_comparison`
        "pass_comparison": {"file": f"{D}/pass_comparison.json", "produced": "runtime",
                            "by": s5, "first_pass_report": first_pass,
                            "final_report": f"{D}/{report}"},
        "seeded": bool(seed), "seed_lib": predicted if seed else None,
        # CPUs per array task, how many run at once, and why (run_search.array_task_cpus)
        "cpu_sizing": cpu_sizing,
        "scripts": [x for x in [s1, s1b, s2, s3, s4, s5, "submit.sh"] if x],
        "submit": f"bash {out}/submit.sh   (or: hive_exec.sh 'bash {out}/submit.sh')",
        "report_jobid_var": "jid5",
        "note": "5-step DIA-NN parallel chain. Submit with submit.sh, watch EVERY step with "
                "watch_run.sh --all, then point run_de.R at the report.",
    }, indent=2))


if __name__ == "__main__":
    main()
