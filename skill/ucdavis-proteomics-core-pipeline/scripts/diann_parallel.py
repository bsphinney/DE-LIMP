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
library round-trip. **Mass accuracy is FIXED (manual), not auto** — steps 3/5 reuse the
.quant files and auto-calibration would be inconsistent (per DIA-NN dev guidance).

Writes into <out>: `file_list.txt`, `step{1..5}_*.sbatch`, and `submit.sh` (submits the
chain with dependencies). Run `submit.sh` on the cluster (or via `hive_exec.sh`). All
heavy work runs on compute nodes through the array — never the login node.

Usage:
  python3 diann_parallel.py --diann '<diann binary | apptainer exec ... diann-linux>' \
      --raw /data/*.d --fasta /path/search.fasta --out ./diann_parallel \
      --cfg params.cfg [--threads-per-file 16] [--mem-per-file 64] [--time-per-file 2] \
      [--assembly-cpus 64] [--assembly-mem 128] [--assembly-time 12] \
      [--partition <auto>] [--account <auto>] [--max-simultaneous 20] [--no-norm]
"""
import os, re, sys, glob, argparse, shlex, subprocess, math

# flags that are step-specific or auto-determined — never carry them into every step.
# NOTE: --dda is intentionally NOT stripped — for DDA data put --dda in the --cfg and it
# flows into every step (DIA-NN 2.6 searches DDA per file exactly as it does DIA).
STRIP = ("--fasta-search", "--predictor", "--gen-spec-lib", "--matrices", "--reanalyse",
         "--rt-profiling", "--no-norm", "--xic", "--mobilograms", "--out-lib", "--lib", "--out", "--f",
         # NOTE: --xic is stripped here on purpose and re-added to step 4 ONLY (see
         # xic_flag() below) -- step 2 IDs are not final, and step 5 runs --use-quant,
         # which never re-reads the raw spectra so --xic is silently a no-op there.
         "--fasta", "--threads", "--temp")

# Step 1b tries this many files before giving up: the radius is a property of the method, so
# one blank or wash first in the list must not fail the cohort. Each attempt normally costs
# minutes (DIA-NN logs the radius during calibration); the per-attempt timeout bounds a file
# that never gets there, and the wall clock covers every attempt.
PROBE_CANDIDATES = 3
PROBE_TIMEOUT_S = 3600      # probe_window.py's own default, and step 1b's effective limit before
                            # it tried more than one file. Nothing shorter has been measured on a
                            # large Astral .raw, so it is not shortened to fit three attempts.
PROBE_WALL_HOURS = -(-PROBE_CANDIDATES * PROBE_TIMEOUT_S // 3600) + 1   # every attempt + 1 h


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
            "window_reason": win_reason, "reason": reason}


def mass_acc_record(ma):
    """Mass accuracy only, for the generator's output and search_provenance.json."""
    return {"fixed": ma["mass_acc_fixed"], "ms1": ma["ms1"], "ms2": ma["ms2"],
            "reason": ma["mass_acc_reason"]}


def window_record(ma):
    """What the cfg hands DIA-NN for --window, for provenance -- exactly what was passed, and
    "unverified" wherever DIA-NN's behaviour has not been measured. Used for the chain when it
    does not probe, and by run_search.py for the single-shot search."""
    st, passed = ma["window_state"], ma["window_passed"]
    if st == "ok":
        return {"source": f"pinned in the cfg ({'; '.join(passed)})", "value": ma["window"],
                "passed": passed}
    if st == "unset":
        return {"source": "not in the cfg -- DIA-NN chooses the radius itself (on the 18-file "
                          "poplar chain it chose per file: 7 for seventeen, 8 for one; how it "
                          "chooses within one multi-file search is unverified)",
                "value": None, "passed": []}
    if st == "zero":
        return {"source": "passed as `--window 0`, which is not a positive integer -- on the "
                          "poplar run DIA-NN warned and chose a radius per file; unverified "
                          "beyond that run",
                "value": None, "passed": passed}
    return {"source": f"passed as given ({'; '.join(passed)}) -- not one positive integer, so "
                      "what DIA-NN does with it is unverified",
            "value": None, "passed": passed}


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
            f"({table}); for an Orbitrap pass --ms1-resolution/--ms2-resolution. "
            "Left on auto it can only run as the single-shot search",
        "mass_acc_invalid":
            "correct --mass-acc/--mass-acc-ms1 in the cfg to ONE positive ppm value each "
            f"({table}), or re-run estimate_params.py with the real instrument. For "
            "auto-calibration delete the flags -- never 0 -- and run single-shot",
        "window_invalid":
            "set --window to one positive integer, or delete it and the chain measures it "
            "itself (step 1b). Do not guess a value: it depends on the acquisition scheme",
        "window_seeded":
            "pin --window in the cfg: a seeded chain has no step 1 for step 1b to follow, so "
            "measure it once with probe_window.py against the seed library on one file",
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
      * mass accuracy unset -- NOT recoverable. DIA-NN calibrates it per run against the
        library, so there is no single value to carry into steps 3/5, which reuse .quant.
      * mass accuracy 0 / negative / non-numeric / non-finite / set twice differently --
        NOT recoverable, and not "auto": 0 is a literal 0 ppm tolerance (0 IDs).
      * --window unset or 0 -- recoverable when probing. Step 1b runs probe_window.py and
        pins one radius into steps 2-5. estimate_params.py cannot supply it: the radius is a
        property of the acquisition scheme and has to be MEASURED on a real file. DIA-NN
        does not accept 0 ("scan window radius should be a positive integer") and optimises
        per file instead -- the very inconsistency step 1b removes.
      * --window anything but a non-negative integer (0.5, nan, -1, 7.0, wide), or set twice
        differently -- NOT recoverable. A typo is a mistake to report, not something to
        quietly measure over.
      * an unparseable cfg (unbalanced quote) -- NOT recoverable.

    Returns {ok, probe, code, ma, reason, remedy}: `probe` says step 1b is needed, `code`
    names the outcome (probe | pinned | cfg_missing | cfg_unparseable | mass_acc_unset |
    mass_acc_invalid | window_invalid | window_seeded | window_no_probe), `remedy` is how to fix
    a refusal. A cfg path that is not a file is `cfg_missing`, never "mass accuracy is not
    pinned": `--sbatch proj` once renamed the folder holding the cfg, and the refusal that
    followed blamed mass accuracy.
    `ma` is mass_acc_status() untouched -- on the probe path its window is still unset,
    because it IS unset until step 1b runs.
    """
    def verdict(ok, probe, code, ma, reason):
        return {"ok": ok, "probe": probe, "code": code, "ma": ma, "reason": reason,
                "remedy": None if ok else _remedy(code)}

    try:
        ma = mass_acc_status(cfg)
    except CfgError as e:
        return verdict(False, False, e.code, None, str(e))
    st = ma["state"]
    # Invalid values before unset ones: `mass_acc_unset` is the one code --allow-auto-mass-acc
    # may override, so it must never be what hides a junk --window behind it.
    if any(st[f] == "invalid" for f in MASS_ACC_FLAGS):
        return verdict(False, False, "mass_acc_invalid", ma,
                       f"mass accuracy is set but not usable ({ma['reason']})")
    if st["--window"] == "invalid":
        return verdict(False, False, "window_invalid", ma,
                       f"--window is set but is not a usable radius ({ma['reason']})")
    if any(st[f] == "unset" for f in MASS_ACC_FLAGS):
        return verdict(False, False, "mass_acc_unset", ma,
                       f"mass accuracy is not pinned ({ma['reason']})")
    if st["--window"] == "ok":
        return verdict(True, False, "pinned", ma, ma["reason"])
    if seed_lib:
        return verdict(False, False, "window_seeded", ma,
                       "--window is unpinned and the first pass is seeded from an existing "
                       "library, so there is no step 1 for step 1b to follow")
    if not probe_window:
        return verdict(False, False, "window_no_probe", ma,
                       "--window is unpinned and --no-probe-window was given")
    return verdict(True, True, "probe", ma,
                   f"MS1 {ma['ms1']} ppm / MS2 {ma['ms2']} ppm; --window is unpinned but "
                   "recoverable -- step 1b measures it and pins it for steps 2-5")


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
    # publicgrp/low is PREEMPTIBLE: without --requeue a preempted task is simply lost.
    if (qos or "").startswith("public") or partition == "low":
        h.append("#SBATCH --requeue")
    h += [f"#SBATCH -o {name}_%j.log", f"#SBATCH -e {name}_%j.log"]
    if array:
        h.insert(2, f"#SBATCH --array={array}")
        h = [x.replace("_%j.log", "_%A_%a.log") for x in h]
    return "\n".join(h)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--diann", required=True, help="DIA-NN command (native binary path, or 'apptainer exec --bind … <sif> /diann-*/diann-linux')")
    ap.add_argument("--raw", nargs="+", default=[], help="raw paths/globs (or use --raw-list)")
    ap.add_argument("--raw-list", help="file with one raw path per line — handles spaces in paths")
    ap.add_argument("--fasta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cfg", help="diann.cfg with the search params (estimate_params.py output)")
    ap.add_argument("--threads-per-file", type=int, default=16)
    ap.add_argument("--mem-per-file", type=int, default=64)
    ap.add_argument("--time-per-file", type=int, default=2)
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
    a = ap.parse_args()

    raws = []
    if a.raw_list:
        with open(a.raw_list) as fh:
            raws.extend(line.strip() for line in fh if line.strip())
    for p in a.raw:
        raws.extend(sorted(glob.glob(p)) or [p])
    raws = [os.path.abspath(r.rstrip("/")) for r in raws]

    # Detect the queue from the submitting user's SLURM associations rather than
    # assuming facility membership. genome-center-grp/high for members; publicgrp/low
    # for everyone else (incl. class accounts) — where `high` caps at 8 CPUs/job, so a
    # 32-CPU request would never start.
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from run_search import slurm_queue
        a.partition, a.account, _q = slurm_queue(a.partition, a.account, None)
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
    fasta = os.path.abspath(a.fasta)
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
        if not (a.allow_auto_mass_acc and safe["code"] == "mass_acc_unset"):
            sys.exit(
                f"Not parallel-safe: {a.cfg or '(no --cfg given)'} -- {safe['reason']}.\n"
                "The 5-step chain reuses .quant files across steps, so anything DIA-NN\n"
                "auto-optimises PER FILE (mass accuracy AND scan window) is applied\n"
                "inconsistently between passes and then stitched together.\n"
                f"Fix: {safe['remedy']}.\n"
                "Or run the single-shot search instead (run_search.py --no-parallel)."
                + ("\nTo override deliberately: --allow-auto-mass-acc."
                   if safe["code"] == "mass_acc_unset" else ""))
        sys.stderr.write(f"[diann_parallel] WARNING: proceeding with auto mass accuracy "
                         f"({safe['reason']}) -- steps will not be mutually consistent.\n")
    win_probe, ma = safe["probe"], safe["ma"]

    # When step 1b measures the radius it is PREFIXED onto steps 2-5, so any --window
    # still in the cfg has to come out or both land on the same command line.
    flags = read_cfg_flags(a.cfg, drop=("--window",) if win_probe else ())
    D = out  # all DIA-NN intermediate/output lives here (real paths; native binary reads them directly)
    report = "no_norm_report.parquet" if a.no_norm else "report.parquet"
    norm = "--no-norm" if a.no_norm else ""
    xic = xic_flag(a.cfg)          # step 4 only -- see xic_flag() docstring

    # file list (1 raw path per line) — array tasks index into it
    open(os.path.join(out, "file_list.txt"), "w").write("\n".join(raws) + "\n")
    all_f = " ".join(f"--f {shlex.quote(r)}" for r in raws)   # quote — data paths may contain spaces
    array = f"0-{n-1}%{a.max_simultaneous}"
    seed = os.path.abspath(a.seed_lib) if a.seed_lib else None
    predicted = seed if seed else f"{D}/step1.predicted.speclib"
    empirical = f"{D}/empirical.parquet"

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
                            peak_cpus=a.threads_per_file, preemptible_ok=True)
        _ps, _as_, _qs = _sq(a.partition, a.account, a.qos, peak_cpus=a.assembly_cpus)
    except Exception:
        _pa, _aa, _qa = a.partition, a.account, a.qos
        _ps, _as_, _qs = a.partition, a.account, a.qos

    def write(name, body):
        p = os.path.join(out, name)
        open(p, "w").write(body + "\n"); os.chmod(p, 0o755); return name

    # array preamble: pick this task's raw file
    pick = ('FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" ' + f'{D}/file_list.txt)\n'
            'if [ -z "$FILE" ]; then echo "no file for task $SLURM_ARRAY_TASK_ID"; exit 1; fi\n'
            'echo "Processing: $FILE"\n')

    # Step 1 — library prediction (single job) — SKIPPED when --seed-lib is given
    s1 = None
    if not seed:
        s1 = write("step1_libpred.sbatch", "\n".join([
            header("s1_libpred", a.libpred_cpus, a.libpred_mem, a.libpred_time, _ps, _as_, qos=_qs), "",
            f'echo "Step 1/5 library prediction"; date',
            f'{DN} --fasta {fasta} --fasta-search --predictor --gen-spec-lib \\',
            f'  --out-lib {D}/step1.speclib --out {D}/step1_lib.parquet \\',
            f'  --threads {a.libpred_cpus} {flags}',
            must_exist(predicted, "the predicted spectral library")]))

    # Step 1b — measure the scan-window radius ONCE, so steps 2-5 share it.
    # DIA-NN optimises the radius per file when --window is absent, and steps 3/5 then
    # combine .quant files produced under different windows -- which DIA-NN's own
    # warning calls "strongly not recommended". Measuring beats guessing: the radius
    # depends on the acquisition scheme (cycle time vs peak width), not the instrument.
    s1b = None
    resolved_cfg = a.cfg
    probe_cands = raws[:PROBE_CANDIDATES]
    if win_probe:
        probe = os.path.join(os.path.dirname(os.path.abspath(__file__)), "probe_window.py")
        # The measured radius has to end up in a PARAMETER FILE, not just window.txt, or the
        # run is not reproducible from what we recorded: search_provenance.json would name a
        # cfg with no --window, and replaying it would re-optimise per file and not reproduce
        # the numbers (SKILL.md golden rule 5). params.base.cfg is the cfg minus any --window,
        # by the same token rule as the step flags. params.resolved.cfg is NOT created here:
        # step 1b builds it in a .tmp and moves it into place only once a radius is measured,
        # so a "resolved" cfg with no --window can never exist to be replayed.
        base_cfg = os.path.join(out, "params.base.cfg")
        resolved_cfg = os.path.join(out, "params.resolved.cfg")
        tmp_cfg = resolved_cfg + ".tmp"
        write_cfg(a.cfg, base_cfg, drop=("--window",))
        q = shlex.quote
        s1b = write("step1b_window.sbatch", "\n".join([
            header("s1b_window", a.threads_per_file, a.mem_per_file, PROBE_WALL_HOURS,
                   _ps, _as_, qos=_qs), "",
            'echo "Step 1b/5 measuring scan-window radius"; date',
            # Every other step reaches DIA-NN through DN, which carries the .NET 8 exports a
            # Thermo .raw needs. probe_window.py runs DIA-NN as its own subprocess (no shell),
            # so the prefix cannot ride on --diann -- it has to be in the ENVIRONMENT the probe
            # inherits. Without it DIA-NN cannot read .raw, no radius is logged, window.txt is
            # never written, and steps 2-5 sit on afterok for ever.
            *([f"{dnet.strip()}   # .NET 8 for Thermo .raw, inherited by probe_window.py's DIA-NN"]
              if dnet else []),
            # A resubmitted step 1b must never find the previous run's answer and carry on.
            f"rm -f {D}/window.txt {D}/window.json {q(resolved_cfg)} {q(tmp_cfg)}",
            # One blank, wash or failed injection first in the list used to fail the whole
            # cohort. The radius is a property of the method, so any good file answers it.
            'W=""',
            f'for RAW in {" ".join(q(r) for r in probe_cands)}; do',
            f"  cp {q(base_cfg)} {q(tmp_cfg)}",
            f'  if python3 {q(probe)} --diann {q(a.diann)} --raw "$RAW" \\',
            f"      --fasta {fasta} --lib {predicted} --threads {a.threads_per_file} \\",
            f"      --timeout {PROBE_TIMEOUT_S} --write-cfg {q(tmp_cfg)} \\",
            # the flags as bash words, after `--`: the probe's DIA-NN gets the same argv as
            # steps 2-5, not a second parse of them through shlex
            f"      -- {flags} > {D}/window.json; then",
            '    W=$(python3 -c "import json,sys; w=json.load(open(sys.argv[1]))[\'window_radius\']; '
            f'assert isinstance(w, int) and w > 0; print(w)" {D}/window.json) && break',
            "  fi",
            '  echo "step 1b: no scan-window radius from $RAW -- trying the next file" >&2',
            '  W=""',
            "done",
            'if [ -z "$W" ]; then',
            f'  echo "FAILED: no scan-window radius from any of the first {len(probe_cands)} '
            'file(s) (DIA-NN log tails above)." >&2',
            # Resubmitting step 1b ALONE does not restart the chain: steps 2-5 were submitted
            # afterok on THIS job id, so they sit PENDING (DependencyNeverSatisfied) for ever.
            '  echo "Steps 2-5 were submitted afterok on THIS job, so they are now PENDING with '
            'DependencyNeverSatisfied and will never start -- even if step 1b is resubmitted '
            'and succeeds." >&2',
            '  echo "Recover: fix the cause (do NOT guess a --window), scancel steps 2-5 (ids in '
            f'{D}/jobs.txt), then resubmit step1b_window.sbatch and steps 2-5 chained afterok on '
            'the new ids, reusing step1.predicted.speclib -- or re-run submit.sh, which also '
            'repeats step 1. See references/watcher.md (dependency_failed)." >&2',
            f"  rm -f {q(tmp_cfg)} {D}/window.json",
            "  exit 1",
            "fi",
            # These are written by this script, not by DIA-NN, so must_exist()'s "DIA-NN exited
            # 0 but did not write" would name the wrong culprit. Say what actually failed.
            f'if ! echo "$W" > {D}/window.txt || [ ! -f {D}/window.txt ] || [ ! -s {D}/window.txt ]; then',
            f'  echo "FAILED: radius $W was measured but could not be written to {D}/window.txt '
            '(disk full? permissions?)" >&2',
            "  exit 1",
            "fi",
            # -f as well as -s: `mv` INTO a directory of that name succeeds, and a directory
            # is non-empty.
            f"if ! mv -f {q(tmp_cfg)} {q(resolved_cfg)} || [ ! -f {q(resolved_cfg)} ] "
            f"|| [ ! -s {q(resolved_cfg)} ]; then",
            f'  echo "FAILED: radius $W was measured but {resolved_cfg} could not be moved into '
            f'place from {tmp_cfg} (disk full? permissions?)" >&2',
            "  exit 1",
            "fi",
            'echo "scan window radius = $W (pinned for steps 2-5)"',
            f'echo "fully-resolved parameters -> {resolved_cfg}"']))
    # steps 2-5 read the measured radius at RUNTIME so every pass uses the identical value
    wflag = f'--window $(cat {D}/window.txt) ' if win_probe else ''

    # DIA-NN aborts with "cannot find the temp folder" if --temp does not exist -- it will NOT
    # create it -- and it does so BEFORE doing any work, so the whole submission cycle is lost to
    # a missing directory. submit.sh makes them, but the watcher playbook (references/watcher.md)
    # tells the orchestrator to resubmit individual steps after a failure, and `sbatch
    # step4_finalpass.sbatch` never goes through submit.sh. So each step makes its own: mkdir -p
    # is idempotent and free, and it means no step can be submitted into that error.
    def tmpguard(d):
        return f'mkdir -p {D}/{d}   # DIA-NN will NOT create --temp and aborts without it'

    # Step 2 — first pass (array): predicted lib -> per-file .quant
    s2 = write("step2_firstpass.sbatch", "\n".join([
        header("s2_firstpass", a.threads_per_file, a.mem_per_file, a.time_per_file, _pa, _aa, qos=_qa, array=array), "",
        f'echo "Step 2/5 first pass, task ${{SLURM_ARRAY_TASK_ID}} of {n}"; date', pick,
        tmpguard("quant_step2"),
        f'{DN} --f "$FILE" --fasta {fasta} --lib {predicted} \\',
        f'  --temp {D}/quant_step2 --rt-profiling --gen-spec-lib --quant-ori-names \\',
        f'  --threads {a.threads_per_file} {wflag}{flags}',
        'QOUT="${FILE##*/}"; QOUT="${QOUT%.*}.quant"',
        must_exist(f'{D}/quant_step2/$QOUT', "this file's .quant")]))

    # Step 3 — empirical library assembly (single job, --use-quant)
    s3 = write("step3_assembly.sbatch", "\n".join([
        header("s3_assembly", a.assembly_cpus, a.assembly_mem, a.assembly_time, _ps, _as_, qos=_qs), "",
        f'echo "Step 3/5 empirical library assembly"; date', tmpguard("quant_step2"),
        f'cp -r {D}/quant_step2 {D}/quant_step2_orig 2>/dev/null || true   # backup for resume',
        f'{DN} {all_f} --fasta {fasta} --lib {predicted} --use-quant --quant-ori-names \\',
        f'  --rt-profiling --gen-spec-lib --out-lib {empirical} \\',
        f'  --temp {D}/quant_step2 --out {D}/step3_assembly.parquet \\',
        f'  --threads {a.assembly_cpus} {wflag}{flags}',
        must_exist(empirical, "the empirical spectral library")]))

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
    xic_out = f' --out {D}/xic/t${{SLURM_ARRAY_TASK_ID}}.parquet' if xic else ''

    s4 = write("step4_finalpass.sbatch", "\n".join([
        header("s4_finalpass", a.threads_per_file, a.mem_per_file, a.time_per_file, _pa, _aa, qos=_qa, array=array), "",
        f'echo "Step 4/5 final pass, task ${{SLURM_ARRAY_TASK_ID}} of {n}"; date', pick,
        tmpguard("quant_step4"),
        'QUANT="${FILE##*/}"; QUANT="${QUANT%.*}.quant"',
        f'if [ ! -f "{D}/quant_step2/$QUANT" ]; then echo "SKIP: no step-2 quant for $QUANT"; exit 0; fi',
        # Splat an empty list, not an empty string: a conditional string leaves a stray
        # blank line in the generated sbatch when XICs are off.
        *([f'mkdir -p {D}/xic'] if xic else []),
        # DIA-NN re-saves the library it is handed as "<lib>.skyline.speclib", written
        # NEXT TO --lib. With one shared path every concurrent array task writes the same
        # file; most win the race in seconds, the losers block until the wall clock kills
        # them. Give each task its own copy so there is nothing to contend on.
        f'LIBPRIV={D}/libpriv/t${{SLURM_ARRAY_TASK_ID}}',
        'mkdir -p "$LIBPRIV"',
        f'cp -f {empirical} "$LIBPRIV/lib.parquet"',
        'trap \'rm -rf "$LIBPRIV"\' EXIT',
        f'{DN} --f "$FILE" --fasta {fasta} --lib "$LIBPRIV/lib.parquet" \\',
        f'  --temp {D}/quant_step4 --quant-ori-names{xic_arg}{xic_out} \\',
        f'  --threads {a.threads_per_file} {wflag}{flags}',
        must_exist(f'{D}/quant_step4/$QUANT', "this file's final-pass .quant")]))

    # Step 5 — cross-run report (single job, --use-quant --matrices)
    s5 = write("step5_report.sbatch", "\n".join([
        header("s5_report", a.assembly_cpus, a.assembly_mem, a.assembly_time, _ps, _as_, qos=_qs), "",
        f'echo "Step 5/5 cross-run report"; date', tmpguard("quant_step4"),
        f'{DN} {all_f} --fasta {fasta} --lib {empirical} --use-quant --quant-ori-names \\',
        f'  --temp {D}/quant_step4 --matrices --out {D}/{report} \\',
        f'  --threads {a.assembly_cpus} {norm} {wflag}{flags}',
        must_exist(f'{D}/{report}', "the cross-run report"),
        # A step-4 task that failed silently leaves no .quant, and step 5 happily
        # reports on whatever survived. Count them: fewer quants than inputs means a
        # sample was dropped, which must never pass as success.
        f'NQ=$(ls -1 {D}/quant_step4/*.quant 2>/dev/null | wc -l | tr -d " ")',
        f'if [ "$NQ" -ne {n} ]; then '
        f'echo "FAILED: report built from $NQ of {n} runs -- a step-4 task produced no .quant." >&2; '
        f'echo "Find it: for f in \\$(cat {D}/file_list.txt); do b=\\$(basename \\"\\$f\\"); '
        f'[ -f {D}/quant_step4/\\${{b%.*}}.quant ] || echo MISSING \\$b; done" >&2; '
        f'exit 1; fi',
        f'echo "OK: report built from all {n} runs"']))

    # submit.sh — chain the steps with afterok dependencies
    sub_lines = ["#!/bin/bash", "set -euo pipefail", f'cd "{out}"',
                 f'mkdir -p "{D}/quant_step2" "{D}/quant_step4"   # DIA-NN --temp dirs MUST pre-exist']
    if seed:
        # Step 1 skipped — the InfinDIA/empirical seed library IS the first-pass lib.
        # First pass optionally waits (afterok) on the lib-build job that produces it.
        dep2 = f'--dependency=afterok:{a.seed_dep} ' if a.seed_dep else ''
        sub_lines += [
            f'echo "Step 1/5 SKIPPED — seeding first pass with {predicted}"',
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
        'jid5=$(sbatch --parsable --dependency=afterok:$jid4 %s)' % s5,
        'echo "submitted: firstpass=$jid2 assembly=$jid3 finalpass=$jid4 report=$jid5"',
        f'echo "final report will be {D}/{report}; watch with: watch_run.sh --slurm $jid5 --log {D}/s5_report_${{jid5}}.log"']

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
        f'python3 "{ck}" record --session "{sess}" --stage search \\',
        f'  --jobs "{jobs}" --desc "DIA-NN 5-step parallel chain ({n} files)" \\',
        f'  --report "{D}/{report}" --watch-job "$jid5" --watch-log "{D}/s5_report_${{jid5}}.log" \\',
        f'  --next "Rscript run_de.R --input {D}/{report} --metadata {sess}/input/conditions.csv --method dpc --outdir {sess}/output/tables" \\',
        '  >/dev/null 2>&1 || true',
        f'printf "%s\\n" {jobs.replace(",", " ")} > "{D}/jobs.txt"',
        f'echo "all chain job ids -> {D}/jobs.txt  (watch the WHOLE chain: watch_run.sh --all {D})"',
        f'echo "recovery notes written to {sess}/RECOVERY.md — you can safely close your terminal"']
    write("submit.sh", "\n".join(sub_lines))

    # Describe what WILL run, not what the cfg says (CLAUDE.md rule 1). On the probe path the
    # radius and the resolved cfg do not exist yet -- they are produced at run time by step
    # 1b -- so they are recorded as such, not as if they were already resolved.
    if win_probe:
        scan_window = {"source": "measured at run time by step 1b (probe_window.py) and pinned "
                                 "for steps 2-5",
                       "value": None, "value_file": f"{D}/window.txt",
                       "probe_candidates": probe_cands}
        resolved = {"file": resolved_cfg, "produced": "runtime", "by": s1b,
                    "note": "written by step 1b only after a radius is measured; absent until "
                            "then, so a missing file after step 1b means step 1b failed"}
    elif safe["ok"]:
        scan_window = dict(window_record(ma), value_file=None)
        resolved = {"file": resolved_cfg, "produced": "generation", "by": None,
                    "note": "the cfg as given already pins mass accuracy and --window"}
    else:
        # --allow-auto-mass-acc. Describe what the steps are actually handed -- a `--window 7`
        # in the cfg IS passed to every step -- rather than assuming the override unpinned it.
        scan_window = dict(window_record(ma), value_file=None)
        resolved = {"file": resolved_cfg, "produced": "generation", "by": None,
                    "note": "NOT fully resolved: mass accuracy is not in the cfg "
                            "(--allow-auto-mass-acc), so DIA-NN chooses it at run time and "
                            "this cfg does not record the value used"}

    import json
    print(json.dumps({
        "out": out, "n_files": n, "report": f"{D}/{report}",
        "parallel_safe": {k: safe[k] for k in ("ok", "probe", "code", "reason")},
        "mass_acc": mass_acc_record(ma),
        "scan_window": scan_window,
        "resolved_params": resolved,
        "seeded": bool(seed), "seed_lib": predicted if seed else None,
        "scripts": [x for x in [s1, s1b, s2, s3, s4, s5, "submit.sh"] if x],
        "submit": f"bash {out}/submit.sh   (or: hive_exec.sh 'bash {out}/submit.sh')",
        "report_jobid_var": "jid5",
        "note": "5-step DIA-NN parallel chain. Submit with submit.sh, watch EVERY step with "
                "watch_run.sh --all, then point run_de.R at the report.",
    }, indent=2))


if __name__ == "__main__":
    main()
