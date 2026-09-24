#!/usr/bin/env python3
"""
estimate_params.py  --  Derive search-engine parameters from the DATA TYPE,
using community known-good settings, instead of hand-maintaining a static cfg.

The values that genuinely vary by data type are the MASS TOLERANCES (set by the
mass analyzer / resolution) and the DIA-vs-DDA window mode. Everything else is a
stable trypsin/LFQ default. Every emitted value is TAGGED with its provenance
(data-type-default / auto-calibration / universal-default) so a reviewer never
mistakes a derived default for a user-confirmed setting (DE-LIMP rule #2).

Known-good mass tolerances are taken from DIA-NN's own recommended-settings
table (verified against the DIA-NN README, June 2026):
    Orbitrap Astral (240k MS1):  MS1 4 ppm,  MS2 10 ppm
    Orbitrap by MS2 resolution:  240k->4, 120k->7, 60k->10, 30k->15 ppm
    Bruker timsTOF (dia-PASEF):  MS1 15 ppm, MS2 15 ppm
    SCIEX TripleTOF / ZenoTOF:   MS1 20 ppm, MS2 20 ppm
    Orbitrap level outside 30k-240k (the other level has a tier):
                                 MEASURED with DIA-NN on representative runs
                                 before the search; the level that has a tier keeps it
    Orbitrap of unknown resolution (no level has a tier):
                                 automatic calibration -- NOT measured, because both
                                 levels would be and DIA-NN warns on a measured MS1;
                                 pass --ms1-resolution/--ms2-resolution to get a tier
    Orbitrap MS1 / ion-trap MS2 (--ms2-analyzer ITMS, e.g. Fusion Lumos OT/IT):
                                 no MS2 value (not an Orbitrap level); automatic
                                 calibration, since either flag fixes both levels
    Unidentified instrument:     automatic calibration (DIA-NN's own default)
Sage's docs give NO instrument-specific tolerances, so Sage ppm windows are
DERIVED from the same per-instrument logic and tagged as such -- except an ion-trap
MS2, which takes Sage's own documented low-res MS/MS fragment window (+/-0.4 Da).

Usage:
  python3 estimate_params.py --engine diann --acquisition DIA \
      --instrument "Orbitrap Astral" --out diann.cfg
  python3 estimate_params.py --engine sage --acquisition DDA \
      --instrument "timsTOF Pro" --out sage_config.json [--var-mods ox]
  python3 estimate_params.py --engine diann --acquisition DIA \
      --instrument "Orbitrap Fusion Lumos" --ms1-resolution 120000 --ms2-resolution 15000 \
      --resolution-source detected --out diann.cfg
      (--ms1-res/--ms2-res are the same flags; --resolution-source is detected when the values
       came from detect_acquisition.py, and defaults to user)
  python3 estimate_params.py --engine sage --acquisition DDA \
      --instrument "Orbitrap Fusion Lumos" --ms1-resolution 120000 --ms2-analyzer ITMS \
      --resolution-source detected --out sage_config.json      # an OT/IT method: MS2 in Da

Writes the engine params file to --out and prints a rationale JSON to stdout
(one entry per setting: value + source). Pass --overrides '<json>' to force
specific fields (e.g. a validated SOP value) — those are tagged "user-override".
"""
import sys, os, json, math, argparse

# --- known-good mass accuracy by instrument class (DIA-NN README) ------------
ORBITRAP_RES_PPM = {240000: 4, 120000: 7, 60000: 10, 30000: 15}
SRC_DIANN = "DIA-NN README recommended settings (verified 2026-06)"


# DIA-NN's own Orbitrap resolution -> mass-accuracy table, verbatim from its README
# ("Changing default settings" §4). Resolution is the thing that actually sets the
# achievable tolerance, which is why a model name alone is not enough.
DIANN_ORBITRAP_PPM = {240000: 4, 120000: 7, 60000: 10, 30000: 15}
SRC_TABLE = "DIA-NN README, Orbitrap resolution->accuracy table"

# The per-instrument (MS1 ppm, MS2 ppm) rows of the same README table. classify_instrument()
# reads them from here, and so does every message that tells a user what to pin -- the
# parallel-chain refusal in diann_parallel.py and the routing decline in run_search.py used
# to spell "timsTOF 15/15, Astral 4/10" out by hand, i.e. a second copy of this table.
DIANN_INSTRUMENT_PPM = {"orbitrap_astral": (4, 10), "timstof": (15, 15), "sciex_tof": (20, 20)}

# The facility's validated SOP mass accuracy: the ONE definition of it in this skill.
# probe_window.pin_mass_acc imports it and uses it as a FLOOR under a measured tolerance.
#
# There was no stored SOP before this. Everything the code called an SOP was a value the CALLER
# supplies at run time -- `resolve_defaults.py --ms1-ppm/--ms2-ppm`, `estimate_params.py
# --overrides` -- and a supplied one pins the cfg (PLAN_PINNED, below), so nothing is measured
# and no floor is needed. make_presets.py's 20/20 is Radiant's and FragPipe's VENDOR default,
# which that file itself calls "visibly wrong ... too wide for narrow-window data": it is not
# this, and grafting it here would have been a second source of truth rather than the one the
# floor needs. These numbers are the FRAN pilot's hand-set override for the Orbitraps it
# re-searches (Exploris 480 and Fusion Lumos, 120k MS1 / 15k MS2), benchmarked in
# references/diann_parallel.md: 19,592 precursors and 2,618 protein groups at 1% FDR on the
# validation cohort -- the widest margin of any candidate tried -- and one of only three
# candidates DIA-NN 2.7.0 logged no deviation warning for.
#
# Change them HERE and the floor moves with them; nothing else hard-codes a tolerance.
SOP_MASS_ACC = {"ms1_ppm": 7.0, "ms2_ppm": 20.0}
SOP_MASS_ACC_SOURCE = ("the facility's validated SOP tolerance for the Orbitraps in the FRAN "
                       "re-search pilot (estimate_params.SOP_MASS_ACC), benchmarked in "
                       "references/diann_parallel.md")


def instrument_ppm_summary():
    """One line of the table, MS1/MS2 ppm, for remediation text."""
    t = DIANN_INSTRUMENT_PPM
    orb = ", ".join(f"{r // 1000}k->{p}" for r, p in sorted(DIANN_ORBITRAP_PPM.items(), reverse=True))
    return (f"MS1/MS2 ppm: timsTOF {t['timstof'][0]}/{t['timstof'][1]}, "
            f"Astral {t['orbitrap_astral'][0]}/{t['orbitrap_astral'][1]}, "
            f"SCIEX {t['sciex_tof'][0]}/{t['sciex_tof'][1]}, "
            f"Orbitrap by resolution ({orb}; a level outside 30k-240k is not extrapolated: "
            f"estimate_params.py plans a DIA-NN measurement of it. An Orbitrap of UNKNOWN "
            f"resolution is not measured either -- both levels would be, and a measured MS1 is "
            f"what DIA-NN warns about -- so pass --ms1-resolution/--ms2-resolution, or DIA-NN "
            f"calibrates per run)")


# How the mass accuracy of a DIA-NN cfg is settled. Written to the rationale sidecar
# (<cfg>.rationale.json) as `mass_accuracy_plan`, which diann_parallel.mass_acc_measure_plan()
# reads for both the 5-step chain and the single-shot search:
#   pinned              -- --mass-acc/--mass-acc-ms1 are in the cfg (a documented value or an SOP)
#   measure_with_diann  -- an Orbitrap with at least one level outside DIA-NN's table; the cfg
#                          omits BOTH flags, and before the search DIA-NN is run on representative
#                          runs (chain step 1b, or run_search.py's single-shot probe). The level
#                          that has a documented tier is recorded under
#                          `mass_accuracy_documented` and pinned as given. An SOP override of
#                          one flag makes the plan pinned (the other level from the table) or
#                          is refused (no table value for it: LoneMassAccOverride).
#   auto                -- instrument not identified; flags omitted, DIA-NN optimises on its own
PLAN_PINNED, MEASURE_WITH_DIANN, PLAN_AUTO = "pinned", "measure_with_diann", "auto"

# Orbitrap classes with a level that has no documented value, measured instead.
#
# `orbitrap_generic` -- an Orbitrap whose RESOLUTION is unknown -- is deliberately NOT here. Since
# 2.5.1 detect_acquisition.py reads a .raw's resolution from its scan trailer and the orchestrator
# passes it in, so this is now the fallback for a .raw whose trailer could not be read, or a Thermo
# mzML with no MS:1000800 term for read_mzml_resolution() to find. With no resolution
# NEITHER level has a tier, so a measure_with_diann plan would measure BOTH -- and measuring MS1
# is the one thing this branch's own evidence says not to do. The first cut measured it and
# pinned 4.2 ppm at 120k; every DIA-NN pass then logged "WARNING: the MS1 mass accuracy setting
# (4.2 ppm) deviates significantly from the value recommended (7 ppm) for the Orbitrap resolution
# of this run (120000)". DIA-NN reads the run's resolution itself and this script does not, so a
# measured MS1 would be pinned against DIA-NN's own advice, silently, on the commonest input the
# facility has.
#
# So an Orbitrap of unknown resolution falls to PLAN_AUTO: both flags omitted, DIA-NN calibrates
# per run, and the 5-step chain declines it (mass_acc_unset) until someone pins a value or passes
# --ms1-resolution/--ms2-resolution -- which reclassifies it as orbitrap_measured (both tiers) or
# orbitrap_untabled (one tier, the other measured) and gets the measurement with a documented MS1.
MEASURE_CLASSES = ("orbitrap_untabled",)

# Where a measured level comes from -- for the DIA-NN cfg only (build_diann). NOT for
# classify_instrument()'s source: resolve_defaults.py and make_presets.py write that into the
# manifest for Radiant and FragPipe too, which measure nothing.
#
# The README's item 6 of "Changing default settings", in full: "One can also optimise all
# parameters to achieve the best possible performance from the data. For this, run DIA-NN on
# several representative runs (best to use any suitable empirical library, as this is the
# quickest) with Unrelated runs option checked and review the 'Averaged recommended settings for
# this experiment' values reported at the end of the log." What is measured here is NOT that
# procedure: it uses the predicted library, one DIA-NN per run rather than one "Unrelated runs"
# search, and pins the median of what each run printed rather than DIA-NN's averaged line. On the
# HIVE validation cohort the two disagreed (MS2 14 here, 16 averaged) -- references/
# diann_parallel.md has the numbers.
SRC_MEASURE = ("measured with DIA-NN before the search: one DIA-NN per representative run (the "
               "median run by size, plus the quartile runs) in automatic mode, stopped once it "
               "logs 'Optimised mass accuracy' (MS2) and 'Recommended MS1 mass accuracy "
               "setting' (MS1); the median of each measured level (a level with a documented "
               "value keeps it) is pinned for every pass -- step 1b of the 5-step chain, or "
               "run_search.py's probe before a single-shot search. Adapted "
               "from DIA-NN README 'Changing default settings' item 6, which averages over an "
               "'Unrelated runs' search instead")

# A lone --mass-acc or --mass-acc-ms1 FIXES BOTH levels in DIA-NN 2.7.0: "WARNING: note the mass
# accuracy settings used by DIA-NN, automatic optimisation will not be performed as at least one
# of MS1/MS2 mass accuracies is user-provided", and the other level falls to 20 ppm -- `--window 7
# --mass-acc-ms1 7` logged "Mass accuracy will be fixed to 2e-05 (MS2) and 7e-06 (MS1)" (HIVE srun
# job 23528991, 2026-09-16). So a documented level is never written into the cfg on its own.
LONE_FLAG_NOTE = ("not written into the cfg on its own: DIA-NN 2.7.0 fixes BOTH levels when either "
                  "flag is given, the other at 20 ppm; it is pinned together with the measured "
                  "level")

MASS_ACC_FLAGS = ("--mass-acc", "--mass-acc-ms1")        # the same pair as diann_parallel's
MASS_ACC_LEVEL = {"--mass-acc": "MS2", "--mass-acc-ms1": "MS1"}


class LoneMassAccOverride(ValueError):
    """An --overrides mass-accuracy flag whose other level has no value: it could only be written
    as a lone flag (LONE_FLAG_NOTE). references/parameters.md's own example, {"--mass-acc": 8},
    once went into the cfg alone for a 120k/15k Orbitrap while the sidecar still planned to measure
    MS2 -- so MS1 ran at 20 ppm instead of its documented 7, and nothing said so."""


def lone_override_message(given, missing, instr_class, label):
    """Why a one-flag override is refused, and what to give instead."""
    lvl, gone = MASS_ACC_LEVEL[given], MASS_ACC_LEVEL[missing]
    why = {"orbitrap_untabled": f"its resolution is outside DIA-NN's 30k-240k Orbitrap table "
                                f"({label})",
           "orbitrap_generic": "the Orbitrap resolution is unknown"}.get(instr_class, label)
    fixes = [f"give both {given} and {missing} in --overrides"]
    if instr_class == "orbitrap_generic":
        fixes.append("pass --ms1-resolution/--ms2-resolution, so DIA-NN's table can supply "
                     f"{gone} if it has a tier")
    if instr_class in MEASURE_CLASSES:
        fixes.append(f"or override neither, and DIA-NN measures the level with no value on "
                     "representative runs before the search")
    elif instr_class == "orbitrap_generic":
        # Not measured: with no resolution BOTH levels would be, and a measured MS1 is what
        # DIA-NN warns about (see MEASURE_CLASSES). Nothing is left but its own calibration.
        fixes.append("or override neither, and DIA-NN calibrates both levels itself, per run "
                     "(the 5-step chain then declines the cfg until --allow-auto-mass-acc)")
    return (f"--overrides sets {given} ({lvl}) but not {missing} ({gone}), and {gone} has no "
            f"value here: {why}. Written alone, {given} makes DIA-NN fix {gone} at 20 ppm too "
            "(DIA-NN 2.7.0: 'automatic optimisation will not be performed as at least one of "
            "MS1/MS2 mass accuracies is user-provided'), so no cfg was written. To fix: "
            + "; ".join(fixes) + ".")


def ppm_for_resolution(res):
    """Map an Orbitrap resolving power to mass accuracy (ppm), per DIA-NN's table.

    Returns (ppm, source), or (None, why) when the table has nothing to say.

    Exact tiers use the documented value. Between 30k and 240k the value is interpolated on a
    log-log fit of the table and TAGGED as interpolated. Outside that range there is NO value:
    the old extrapolation turned an Exploris 480's 15,000 MS2 into 23.3 ppm and pinned it for a
    whole cohort, with nothing but a curve fit behind it -- and 15k is not a corner case, it is
    what both Orbitraps in the FRAN re-search pilot (2026-09-16) acquire MS2 at. A number with no
    evidence is a fabricated default; the caller measures it instead (MEASURE_WITH_DIANN).
    """
    if not res or res <= 0:
        return None, None
    res = float(res)
    if int(res) in DIANN_ORBITRAP_PPM:
        return DIANN_ORBITRAP_PPM[int(res)], SRC_TABLE
    pts = sorted(DIANN_ORBITRAP_PPM.items())
    if not pts[0][0] <= res <= pts[-1][0]:
        return None, (f"{int(res):,} is outside DIA-NN's documented 30k-240k Orbitrap table, "
                      "not extrapolated")
    # log(ppm) is close to linear in log(resolution) across the documented tiers
    (r1, p1), (r2, p2) = pts[0], pts[-1]
    slope = (math.log(p2) - math.log(p1)) / (math.log(r2) - math.log(r1))
    ppm = round(math.exp(math.log(p1) + slope * (math.log(res) - math.log(r1))), 1)
    return ppm, f"{SRC_TABLE}, interpolated for {int(res):,} resolution"


# Where an Orbitrap's MS1/MS2 resolution came from. The phrase goes into the class label, and
# from there into workflow.manifest.json, <cfg>.rationale.json and every mass-accuracy rationale
# line. gabrig 2026-09-23 (Fusion Lumos, skill 2.5.0): resolutions the user TYPED IN came out as
# "MS1 60,000 / MS2 15,000 resolution read from the data" -- a value presented as something it
# is not (DE-LIMP rule #2). The caller says where the numbers came from; a number given on the
# command line is "user" unless the caller says otherwise.
RESOLUTION_SOURCES = {
    "detected": "read from the raw file (scan trailer)",       # detect_acquisition.py's values
    "user": "supplied by the user",
    "mzml": "read from the mzML (MS:1000800 resolving power)",  # estimate_params.py --from-mzml
    "cfg": "taken from a saved configuration, not read from the data",
}
RESOLUTION_SOURCE_UNRECORDED = "(source not recorded)"

# The MS2 of an Orbitrap-MS1 / ion-trap-MS2 method (a Fusion Lumos OT/IT DDA, say): the scan filter
# reads "ITMS" and there is no Orbitrap resolution for MS2 at all, so DIA-NN's resolution table has
# nothing for it. detect_acquisition.py reports it per file as `ms2_analyzer` and lists the files
# under `ms2_ion_trap`. Before --ms2-analyzer existed, such a run got a 10 ppm Sage FRAGMENT window
# (an ion trap needs ~0.4 Da) and, given MS1 alone, a DIA-NN label saying the resolution was unknown.
# "mixed": a file with MS2 scans from BOTH analyzers (e.g. Orbitrap HCD + ion-trap CID in one
# method). One cfg has one fragment window, and it has to fit the ion-trap spectra, so mixed is
# treated as ITMS for every tolerance and labelled as mixed.
MS2_ANALYZERS = ("FTMS", "ITMS", "mixed")
ION_TRAP_ANALYZERS = ("ITMS", "mixed")
ION_TRAP_MS2 = "read in the ion trap (ITMS) -- no Orbitrap-table MS2 value"
MIXED_MS2 = ("read partly in the ion trap and partly in the Orbitrap (mixed ITMS/FTMS) -- no "
             "one Orbitrap-table MS2 value covers it")


def ms2_analyzer_arg(value):
    """argparse type: FTMS/ITMS in capitals, "mixed" as detect_acquisition.py writes it."""
    v = (value or "").strip()
    return "mixed" if v.lower() == "mixed" else v.upper()


def ms2_in_ion_trap(analyzer):
    """True when some or all MS2 is read in the ion trap -- the tolerances must fit it."""
    return analyzer in ION_TRAP_ANALYZERS


def add_resolution_args(ap):
    """The resolution flags, defined once for estimate_params.py, resolve_defaults.py and
    make_presets.py. gabrig 2026-09-23: resolve_defaults.py took --ms1-res/--ms2-res and this
    script --ms1-resolution/--ms2-resolution, so one of the two commands failed whichever
    spelling was used. Both spellings are accepted everywhere; --ms1-resolution is canonical."""
    for lvl in ("ms1", "ms2"):
        ap.add_argument(f"--{lvl}-resolution", f"--{lvl}-res", dest=f"{lvl}_resolution",
                        type=float, default=None,
                        help=f"Orbitrap {lvl.upper()} resolving power (e.g. 120000); maps to ppm "
                             f"via DIA-NN's table. --{lvl}-res is the same flag")
    ap.add_argument("--ms2-analyzer", type=ms2_analyzer_arg, choices=MS2_ANALYZERS, default=None,
                    help="where MS2 was read: FTMS = the Orbitrap, ITMS = the ion trap (step 2's "
                         "ms2_ion_trap files), mixed = both in one file (treated as ITMS for "
                         "tolerances). ITMS has no MS2 resolution: do not pass --ms2-resolution "
                         "with it")
    ap.add_argument("--resolution-source", choices=sorted(RESOLUTION_SOURCES), default=None,
                    help="where the resolutions and --ms2-analyzer came from: detected = "
                         "detect_acquisition.py read them from the .raw; user = the user said so "
                         "(the default for a value given on the command line); cfg = a saved "
                         "configuration")


def resolution_args_error(ms2_res, ms2_analyzer):
    """Why this combination cannot be right, or None. An ion-trap MS2 has no Orbitrap resolution,
    so both at once means one of them is wrong -- and nothing here can tell which."""
    if ms2_analyzer == "ITMS" and ms2_res:
        return ("--ms2-analyzer ITMS and --ms2-resolution contradict each other: an MS2 read in "
                "the ion trap has no Orbitrap resolution. Pass the one step 2 reported for these "
                "files (ms2_ion_trap -> --ms2-analyzer ITMS; ms2_resolution -> --ms2-resolution).")
    return None


def resolution_source(given, ms1_res, ms2_res, ms2_analyzer=None):
    """The source to record: None with nothing given, else the caller's word, else "user".
    The analyzer shares it: detect_acquisition.py reads both from the same scan filters."""
    if not (ms1_res or ms2_res or ms2_analyzer):
        return None
    return given or "user"


def resolution_record(ms1_res, ms2_res, source, ms2_analyzer=None):
    """The `resolution` block of the manifest and the rationale sidecar, or None."""
    if not (ms1_res or ms2_res or ms2_analyzer):
        return None
    as_int = lambda v: int(v) if v and float(v).is_integer() else v  # noqa: E731
    return {"ms1": as_int(ms1_res), "ms2": as_int(ms2_res), "ms2_analyzer": ms2_analyzer,
            "source": source,
            "source_label": RESOLUTION_SOURCES.get(source, RESOLUTION_SOURCE_UNRECORDED)}


# Engines whose settings the Orbitrap resolution changes: DIA-NN (its documented mass-accuracy
# table) and Radiant (extraction widths derived from that table by make_presets.py). NOT
# FragPipe -- make_presets.py keeps its vendor preset tolerances either way -- and NOT Sage,
# whose window sage_ppm() picks by instrument class. Asking for a number the route ignores
# would be a question with no effect on the search.
RESOLUTION_ENGINES = ("diann", "radiant")


def resolution_question(instr_class, engine, ms1_res=None, ms2_res=None):
    """What to ask the user when an Orbitrap's resolution is unknown and the engine uses it,
    else None. gabrig 2026-09-23: a Fusion Lumos with no resolution got a manifest saying
    "resolution unknown", exit 0, and no prompt -- so nobody asked. The class test is
    classify_instrument()'s own: orbitrap_generic IS "an Orbitrap with no usable resolution"
    (the Astral is its own class and assumes 240k).

    An ion-trap MS2 (orbitrap_iontrap) is never asked for: it HAS no resolution. Its MS1 is
    asked for only by Radiant, which narrows the MS1 extraction width from it; DIA-NN pins
    neither level for it (either flag fixes both, and there is no MS2 value), so an MS1 answer
    would change nothing there."""
    if instr_class == "orbitrap_iontrap":
        missing = "MS1" if (engine == "radiant" and not ms1_res) else ""
    elif instr_class == "orbitrap_generic" and engine in RESOLUTION_ENGINES:
        missing = " and ".join(lvl for lvl, v in (("MS1", ms1_res), ("MS2", ms2_res)) if not v)
    else:
        missing = ""
    if not missing:
        return None
    levels = [lvl for lvl in ("MS1", "MS2") if lvl in missing]
    flags = "/".join(f"--{lvl.lower()}-resolution" for lvl in levels)
    example = " / ".join({"MS1": "120,000 MS1", "MS2": "15,000 MS2"}[lvl] for lvl in levels)
    keep = ", keeping --ms2-analyzer ITMS" if instr_class == "orbitrap_iontrap" else ""
    return (f"Ask the user for the Orbitrap {missing} resolution (it is in the instrument "
            f"method, e.g. {example}), or run detect_acquisition.py, which reads it from the "
            f".raw scan trailer. Then re-run with {flags}{keep} (add --resolution-source "
            "detected when the values came from detect_acquisition.py). Without it DIA-NN's "
            "documented Orbitrap mass-accuracy table cannot be applied.")


def classify_instrument(name, ms1_res=None, ms2_res=None, res_source=None, ms2_analyzer=None):
    """Return (class, ms1_ppm, ms2_ppm, label, source). None ppm => not pinned.

    `res_source` (a RESOLUTION_SOURCES key) says where ms1_res/ms2_res came from; the label
    names it, and says "(source not recorded)" rather than guess when it is None.

    `ms2_analyzer` "ITMS" -- MS2 read in the ion trap -- is class orbitrap_iontrap: the MS1 level
    from the table as usual, and no MS2 value whatever ms2_res says (there is no Orbitrap MS2).
    It used to be dropped: the MS1-only resolution fell through to orbitrap_generic, labelled as
    though nothing were known.

    DIA-NN asks for these to be FIXED rather than auto-optimised, and not only for
    speed: "This optimisation is inherently noisy: even replicate injections may not
    produce identical results, and therefore the analysis results will depend on which
    run is first in the list." Auto therefore makes a result depend on FILE ORDER,
    which also blocks the 5-step parallel chain (steps 3/5 reuse .quant files).

    `source` names each level separately when both resolutions are read. It used to be ONE
    level's source stamped on both flags -- MS2's when MS2 was extrapolated, else MS1's -- so a
    documented 120k MS1 read "EXTRAPOLATED", and an interpolated MS2 read as documented.

    An Orbitrap with a level outside the table (class orbitrap_untabled) returns the documented
    value for the level that has one and None for the other; orbitrap_generic (resolution
    unknown) returns None for both. A None level is measured with DIA-NN (see mass_acc_plan).
    The documented level is NOT measured away: an earlier cut of this branch measured both, and
    the pinned MS1 of 4.2 ppm at 120k made every DIA-NN pass log "WARNING: the MS1 mass accuracy
    setting (4.2 ppm) deviates significantly from the value recommended (7 ppm) for the Orbitrap
    resolution of this run (120000)" (HIVE, compare/probe_median).

    `source` says only where each value comes from, never what a search route does about a
    missing one: resolve_defaults.py writes it as the manifest's ppm_source for every engine.
    """
    n = (name or "").strip().lower()
    origin = RESOLUTION_SOURCES.get(res_source, RESOLUTION_SOURCE_UNRECORDED)
    if ms2_in_ion_trap(ms2_analyzer):
        p1, s1 = ppm_for_resolution(ms1_res)
        if ms2_analyzer == "mixed":
            ms2_label, ms2_src = "MS2 read partly in the ion trap (mixed ITMS/FTMS)", MIXED_MS2
        else:
            ms2_label, ms2_src = "MS2 read in the ion trap (ITMS)", ION_TRAP_MS2
        if not ms1_res:
            label = f"Orbitrap, {ms2_label}, {origin}; MS1 resolution unknown"
            ms1_src = "MS1: resolution unknown"
        else:
            label = f"Orbitrap, MS1 {int(ms1_res):,} resolution / {ms2_label}, {origin}"
            ms1_src = (f"MS1 {int(ms1_res):,}: {s1}" if p1 is not None else
                       f"MS1 {int(ms1_res):,}: no documented DIA-NN value ({s1})")
        return ("orbitrap_iontrap", p1, None, label, f"{ms1_src}; MS2: {ms2_src}")
    # A given resolution beats any model-name guess.
    if ms1_res or ms2_res:
        p1, s1 = ppm_for_resolution(ms1_res)
        p2, s2 = ppm_for_resolution(ms2_res)
        if p1 and p2:
            src = s1 if s1 == s2 else f"MS1: {s1}; MS2: {s2}"
            return ("orbitrap_measured", p1, p2,
                    f"Orbitrap, MS1 {int(ms1_res):,} / MS2 {int(ms2_res):,} resolution "
                    f"{origin}", src)
        if ms1_res and ms2_res:
            src = "; ".join(
                f"{lvl} {int(res):,}: {s}" if p is not None else
                f"{lvl} {int(res):,}: no documented DIA-NN value ({s})"
                for lvl, res, p, s in (("MS1", ms1_res, p1, s1), ("MS2", ms2_res, p2, s2)))
            return ("orbitrap_untabled", p1, p2,
                    f"Orbitrap, MS1 {int(ms1_res):,} / MS2 {int(ms2_res):,} resolution "
                    f"{origin}", src)
    if not n:
        return ("unknown", None, None, "instrument not detected", "auto-calibration fallback")
    if "astral" in n:
        return ("orbitrap_astral", *DIANN_INSTRUMENT_PPM["orbitrap_astral"], "Orbitrap Astral (assumes 240k MS1)", SRC_DIANN)
    if "tims" in n:
        return ("timstof", *DIANN_INSTRUMENT_PPM["timstof"], "Bruker timsTOF (dia-PASEF / ddaPASEF)", SRC_DIANN)
    if "tripletof" in n or "zenotof" in n or "sciex" in n:
        return ("sciex_tof", *DIANN_INSTRUMENT_PPM["sciex_tof"], "SCIEX TripleTOF / ZenoTOF", SRC_DIANN)
    # Orbitrap family with no resolution to work from -> measured with DIA-NN (DIA-NN route).
    # Passing --ms1-resolution/--ms2-resolution (read_mzml_resolution() gets them straight out
    # of the mzML) lets a documented tier pin it instead.
    if any(k in n for k in ("orbitrap", "exploris", "exactive", "fusion", "lumos",
                            "eclipse", "velos", "hf-x", "hf", "qe", "astral")):
        # One level given: say which is known rather than "resolution unknown". Still no value
        # for either -- a lone DIA-NN flag fixes both levels (LONE_FLAG_NOTE).
        for have, res, miss in (("MS1", ms1_res, "MS2"), ("MS2", ms2_res, "MS1")):
            if res:
                return ("orbitrap_generic", None, None,
                        f"Orbitrap, {have} {int(res):,} resolution {origin}; {miss} resolution "
                        f"unknown -- pass --{miss.lower()}-resolution to use DIA-NN's "
                        "documented table",
                        f"no documented DIA-NN value ({miss} resolution unknown)")
        return ("orbitrap_generic", None, None,
                "Orbitrap (resolution unknown -- pass --ms1-resolution/--ms2-resolution to "
                "use DIA-NN's documented table)",
                "no documented DIA-NN value (resolution unknown)")
    return ("unknown", None, None, f"unrecognized instrument '{name}'", "auto-calibration fallback")


def mass_acc_plan(instr_class, ms1, ms2):
    """pinned | measure_with_diann | auto -- see PLAN_PINNED above. The one place that decides."""
    if ms1 and ms2:
        return PLAN_PINNED
    return MEASURE_WITH_DIANN if instr_class in MEASURE_CLASSES else PLAN_AUTO


def read_mzml_resolution(path, max_bytes=4_000_000):
    """(ms1_res, ms2_res) from an mzML's `mass resolving power` (MS:1000800) terms.

    msconvert preserves the per-scan resolving power, so the value DIA-NN needs is
    already in the file -- no vendor reader required. First value seen is MS1, the
    most common subsequent one is MS2.
    """
    import re, collections
    try:
        with open(path, "rb") as fh:
            head = fh.read(max_bytes).decode("utf-8", "replace")
    except OSError:
        return None, None
    vals = [float(v) for v in
            re.findall(r'MS:1000800"[^/]*?value="([0-9.]+)"', head)]
    if not vals:
        return None, None
    ms1 = max(vals)                                  # MS1 is acquired at higher res
    others = [v for v in vals if v != ms1]
    ms2 = collections.Counter(others).most_common(1)[0][0] if others else ms1
    return ms1, ms2


def tagged(value, source):
    return {"value": value, "source": source}


# --- DIA-NN cfg --------------------------------------------------------------
def build_diann(acq, instr_class, ms1, ms2, label, src, var_mods, overrides,
                cont_tag=None, mz_range=None, level_src=None):
    UNIV = "universal trypsin/LFQ default"
    r = {}  # rationale
    lines = []

    def add(flag, value, source, render=None):
        r[flag] = tagged(value, source)
        if render is False:
            return
        if value is True:
            lines.append(flag)
        else:
            lines.append(f"{flag} {value}")

    add("--qvalue", 0.01, "standard 1% precursor FDR")
    add("--matrices", True, UNIV)
    # Extracted ion chromatograms: needed to visually validate an identification
    # (DIA-NN XIC Viewer / Skyline) rather than trusting a q-value alone. 10s is
    # DIA-NN's own default window and is cheap; the documented disk-space caveat
    # applies to --xic-theoretical-fr and large windows, which are NOT enabled here.
    add("--xic", 10, "extract XICs for identification validation (DIA-NN default 10s window)")
    # --xic ALONE writes mobilogram files containing nothing but zeros. DIA-NN
    # allocates ms1_mobilogram/ms2_mobilogram parquets but never populates them
    # unless --mobilograms is ALSO passed -- silently, at plausible file size.
    # On timsTOF the mobilograms are the ion-mobility axis: a real peptide should
    # form a coherent peak in BOTH retention time and mobility, which a
    # chromatogram alone cannot show. DIA-NN ignores the flag on instruments
    # without ion mobility, so it is safe to emit unconditionally.
    add("--mobilograms", True, "populate the mobilogram files -- --xic alone leaves them all zeros")
    add("--fasta-search", True, "library-free (predicted spectral library)")
    add("--gen-spec-lib", True, UNIV)
    add("--predictor", True, "deep-learning predictor for library-free DIA")
    add("--reanalyse", True, "MBR across the run set")
    add("--rt-profiling", True, UNIV)
    add("--cut", "K*,R*", "trypsin/P")
    add("--missed-cleavages", 1, "DIA-NN default")
    add("--min-pep-len", 7, UNIV)
    add("--max-pep-len", 30, UNIV)
    # Precursor m/z search range.
    #
    # This MUST follow the acquisition. Searching narrower than was acquired
    # discards real data, does not error, and is invisible in the output --
    # measured on a UC Davis two-platform pilot where the old hardcoded 380-980
    # was emitted for both instruments and matched neither:
    #     timsTOF HT dia-PASEF   acquired 299.5-1200.5  -> lost 300-380, 980-1200
    #     Orbitrap Lumos DIA     acquired  350 -1200    -> lost 350-380, 980-1200
    # detect_acquisition.py already reads the isolation windows to classify
    # DIA vs DDA; it now returns their bounds too, and run_search.py passes them
    # in. Round outward to whole m/z so we never clip the edge window.
    if mz_range:
        lo, hi = float(mz_range[0]), float(mz_range[1])
        add("--min-pr-mz", int(math.floor(lo)),
            f"measured from the acquired isolation windows ({lo:.1f}-{hi:.1f} m/z)")
        add("--max-pr-mz", int(math.ceil(hi)),
            f"measured from the acquired isolation windows ({lo:.1f}-{hi:.1f} m/z)")
    else:
        # Tagged FALLBACK, not UNIV: a reader must be able to tell a measured
        # range from a guess, and this guess can silently cost identifications.
        add("--min-pr-mz", 380,
            "FALLBACK -- acquired range unknown for this input; NOT measured. "
            "If the method acquired outside 380-980, widen it")
        add("--max-pr-mz", 980,
            "FALLBACK -- acquired range unknown for this input; NOT measured. "
            "If the method acquired outside 380-980, widen it")
    # NOT a universal default: DIA-NN's own default is 1-4. Measured on the 2.6.0
    # binary with a 60-protein FASTA -- no charge flags and --min-pr-charge 1 both
    # give 10,899 precursors, --min-pr-charge 2 gives 8,805. So this discards ~19%
    # of the predicted library. That is the right call for tryptic bottom-up work
    # (z=1 precursors are rarely informative), but a reader must be able to see
    # that a deliberate narrowing was applied on their behalf.
    add("--min-pr-charge", 2,
        "z=1 excluded: rarely informative for tryptic bottom-up. DIA-NN's own "
        "default is 1-4; this drops ~19% of the predicted library (measured: "
        "10,899 -> 8,805 precursors on a 60-protein FASTA)")
    add("--max-pr-charge", 4, UNIV)
    add("--min-fr-mz", 200, UNIV)
    add("--max-fr-mz", 1800, UNIV)
    add("--met-excision", True, UNIV)
    add("--unimod4", True, "fixed carbamidomethyl (C); DIA-NN recommended fixed mod")

    # variable mods: DIA-NN README says var mods don't help pure quant; default OFF.
    if var_mods and "ox" in var_mods:
        add("--var-mods", 1, "user requested Ox(M)")
        add("--var-mod", "UniMod:35,15.994915,M", "oxidation of methionine (user requested)")
    else:
        r["variable_mods"] = tagged("none",
            "DIA-NN README: variable mods do not improve depth for relative quant; left off")

    # mass accuracy — the data-type-dependent part.
    #
    # CRITICAL: `--mass-acc 0` is NOT "automatic" in DIA-NN. It fixes the tolerance
    # at a literal 0 ppm, so no fragment can ever match and the run yields 0 IDs
    # ("Mass accuracy will be fixed to 0 (MS2) and 0 (MS1)" in the log). Automatic
    # calibration is what you get by OMITTING the flags entirely. So when we have no
    # instrument-derived value, record the rationale but emit nothing (render=False).
    # diann_parallel.parallel_safe() declines the chain for both an absent flag (auto) and
    # a 0 (literal 0 ppm) -- for different reasons, which its message now names -- EXCEPT an
    # absent pair this function planned to measure (MEASURE_WITH_DIANN, below).
    plan = mass_acc_plan(instr_class, ms1, ms2)
    documented, basis = {}, {}
    # An SOP override of a mass-accuracy level (applied below, with the other overrides) is a
    # value for that level. Both levels known -- overridden, or from the table -- means both are
    # written and nothing is measured. A one-flag override whose other level has NO value would
    # be a lone flag (LONE_FLAG_NOTE: the other level silently fixed at 20 ppm), so it is refused.
    given = [f for f in MASS_ACC_FLAGS if f in (overrides or {})]
    if given:
        table = {"--mass-acc-ms1": ms1, "--mass-acc": ms2}
        missing = [f for f in MASS_ACC_FLAGS if f not in given and table[f] is None]
        if missing:
            raise LoneMassAccOverride(lone_override_message(given[0], missing[0], instr_class,
                                                            label))
        plan = PLAN_PINNED            # an SOP value wins; it must not be re-measured
    if plan == MEASURE_WITH_DIANN:
        # No curve fit for the level outside the table (see ppm_for_resolution): it is measured
        # with DIA-NN before the search -- chain step 1b, or run_search.py's single-shot probe --
        # and written to params.resolved.cfg. diann_parallel.mass_acc_measure_plan() accepts the
        # omitted flags ONLY because of this plan. A level WITH a tier keeps it, recorded under
        # `mass_accuracy_documented`, but neither flag is rendered (LONE_FLAG_NOTE). Without a
        # library before the search (--one-step) nothing can be measured; run_search.py says so.
        s1, s2 = level_src or (src, src)
        for flag, level, value, lsrc in (("--mass-acc-ms1", "MS1", ms1, s1),
                                         ("--mass-acc", "MS2", ms2, s2)):
            if value is not None:
                documented[flag] = value
                basis[flag] = lsrc
                add(flag, value, f"{label}: {level} {value} ppm [{lsrc}] -- {LONE_FLAG_NOTE}",
                    render=False)
            else:
                add(flag, "measured with DIA-NN before the search (flag omitted)",
                    f"{label}: {level}: {lsrc}; {SRC_MEASURE}", render=False)
    elif plan == PLAN_AUTO:
        # An ion-trap MS2 lands here with MS1 known: DIA-NN fixes BOTH levels when either flag is
        # given (LONE_FLAG_NOTE) and there is no MS2 value, so neither is written.
        why = ("MS1 not pinned alone: DIA-NN 2.7.0 fixes BOTH levels when either flag is given, "
               "the other at 20 ppm, and there is no MS2 value to pin with it; "
               if instr_class == "orbitrap_iontrap" and ms1 is not None else "")
        auto_src = (f"{src}; {why}flags omitted so DIA-NN auto-calibrates per run "
                    "(--mass-acc 0 would pin the tolerance at 0 ppm -> 0 IDs)")
        add("--mass-acc", "auto (flag omitted)", auto_src, render=False)
        add("--mass-acc-ms1", "auto (flag omitted)", auto_src, render=False)
    else:
        # level_src = (MS1 source, MS2 source) when both came from resolutions, so each flag
        # names ITS level's evidence (a documented tier vs an interpolation), not a merged one.
        # An overridden level is written by the override loop below, with its own tag.
        s1, s2 = level_src or (src, src)
        if "--mass-acc" not in given:
            add("--mass-acc", ms2, f"{label}: MS2 {ms2} ppm [{s2}]")
        if "--mass-acc-ms1" not in given:
            add("--mass-acc-ms1", ms1, f"{label}: MS1 {ms1} ppm [{s1}]")
    # --window 0 is likewise rejected ("scan window radius should be a positive
    # integer"); omitting it lets DIA-NN set the radius from the observed peak width.
    #
    # What happens to the omitted flag depends on the ROUTE, which is decided later, by
    # run_search.py -- so this names both instead of asserting one. A single-shot search
    # leaves DIA-NN to choose the radius (how, within one multi-file search, is unverified);
    # the 5-step parallel chain measures it ONCE in step 1b and pins it for steps 2-5. Saying "auto-optimised per run"
    # here was wrong for the chain, which is the route most real cohorts take. The route
    # that actually ran is recorded in search_provenance.json (`scan_window`).
    add("--window", "not set here (flag omitted)",
        "radius depends on the acquisition scheme and must be measured, not guessed. "
        "Single-shot search: DIA-NN chooses it itself (how, within one multi-file search, "
        "is unverified). 5-step parallel chain: step 1b "
        "measures it once (probe_window.py) and pins it for steps 2-5, writing "
        "params.resolved.cfg. search_provenance.json `scan_window` records which ran",
        render=False)

    # Contaminants: identify them, but keep them out of quant + normalisation.
    # The tag comes from fetch_fasta.py's sidecar so the cfg self-describes the
    # database it was built for -- rather than us assuming contaminants are present.
    if cont_tag:
        add("--cont-quant-exclude", cont_tag,
            f"contaminants ('{cont_tag}'-tagged) excluded from quant + normalisation "
            f"(DIA-NN README --cont-quant-exclude)")

    # apply overrides (e.g. a validated SOP value) — re-render the file
    for k, v in (overrides or {}).items():
        r[k] = tagged(v, "user-override (validated SOP)")
        lines = [ln for ln in lines if not (ln == k or ln.startswith(k + " "))]
        lines.append(k if v is True else f"{k} {v}")
    r["mass_accuracy_plan"] = tagged(plan, {
        PLAN_PINNED: "--mass-acc/--mass-acc-ms1 are in the cfg",
        MEASURE_WITH_DIANN: SRC_MEASURE,
        # It said "instrument not identified" for every auto cfg -- also for a named Orbitrap of
        # unknown resolution, and for an ion-trap MS2, neither of which is unidentified.
        PLAN_AUTO: {"orbitrap_iontrap": "MS2 read in the ion trap (all of it, or part: see the "
                                        "class label) -- no Orbitrap-table MS2 value; mass "
                                        "accuracy left to DIA-NN calibration (not parallel-safe)",
                    "orbitrap_generic": "Orbitrap resolution unknown; DIA-NN optimises it itself "
                                        "(not parallel-safe)",
                    }.get(instr_class, "instrument not identified; DIA-NN optimises it itself "
                                       "(not parallel-safe)"),
    }[plan])
    # Each level's basis, not one claim for all: a value between tiers (MS1 at 90k -> 7.5) is
    # interpolated from the table, and the gate and provenance once called it "as documented".
    r["mass_accuracy_documented"] = tagged(
        documented,
        ("levels pinned as given alongside the measured level, each a value from DIA-NN's Orbitrap "
         "resolution table: " + "; ".join(f"{f} {v} [{basis[f]}]" for f, v in documented.items()))
        if documented else
        "none: only a measure_with_diann plan pins a level from the table alongside a measured one")

    return "\n".join(lines) + "\n", r


# --- Sage config -------------------------------------------------------------
def sage_ppm(instr_class):
    """Sage docs give no instrument tolerances; derive from DIA-NN per-instrument."""
    src = "derived from DIA-NN per-instrument recommendation (Sage docs give none)"
    if instr_class == "orbitrap_astral":   return (10, 10, src)   # high-res Orbitrap
    if instr_class == "timstof":           return (20, 20, src)
    if instr_class == "sciex_tof":         return (40, 40, src)
    # Every other Orbitrap class, not only orbitrap_generic: GIVING the resolution reclassifies an
    # Orbitrap as orbitrap_measured / orbitrap_untabled, and those used to fall through to 20/20
    # "instrument not identified" -- a wider window, and a label that was false, for the same
    # instrument that got 10/10 when nobody knew its resolution.
    if instr_class.startswith("orbitrap"): return (10, 10, "high-res Orbitrap default (derived)")
    return (20, 20, "safe high-res default (instrument not identified)")


# Ion-trap MS2 (--ms2-analyzer ITMS) in Sage. Every value is Sage's own documented low-res MS/MS
# setting, read 2026-09-24 -- not derived, and not DIA-NN's (whose table is Orbitrap-only):
#   https://sage-docs.vercel.app/docs/configuration/tolerance -- "For high-res MS/MS:
#     { "fragment_tol": { "ppm": [-10, 10] } } Or for low-res MS/MS: { "fragment_tol":
#     { "da": [-0.4, 0.4] } }"
#   https://sage-docs.vercel.app/docs/configuration/spectra -- "Recommended settings for low-res
#     MS/MS { "deisotope": false, "min_peaks": 15, "max_peaks": 150, "min_matched_peaks": 4,
#     "max_fragment_charge": 2 }" (the last four are what every Sage cfg here already uses)
#   https://sage-docs.vercel.app/docs/configuration -- bucket_size: "Use a lower value (8192) for
#     high-res MS/MS, and higher values for low-res MS/MS", 32768 in the same example
# The Da form is valid in the pinned 0.14.7: crates/sage/src/mass.rs has
# `#[serde(rename_all = "lowercase")] pub enum Tolerance { Ppm(f32, f32), Da(f32, f32) }`, and
# crates/sage-cli/src/input.rs reads `fragment_tol: Tolerance`. With a 10 ppm fragment window
# an ion-trap spectrum (unit resolution) matches almost no fragments.
SAGE_ITMS_FRAGMENT_DA = 0.4
SAGE_ITMS_BUCKET_SIZE = 32768
SRC_SAGE_LOWRES = "Sage docs' low-res MS/MS setting (sage-docs.vercel.app/docs/configuration)"


def build_sage(acq, instr_class, var_mods, overrides, ms2_analyzer=None):
    prec_ppm, frag_ppm, ppm_src = sage_ppm(instr_class)
    ion_trap = ms2_in_ion_trap(ms2_analyzer)
    UNIV = "universal trypsin/LFQ default"
    variable = {"M": [15.9949]} if (var_mods and "ox" in var_mods) else {}
    variable["["] = [42.0106]  # protein N-term acetyl is a common, cheap variable mod
    cfg = {
        "database": {
            "bucket_size": SAGE_ITMS_BUCKET_SIZE if ion_trap else 8192,
            "enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 30,
                       "cleave_at": "KR", "restrict": "P"},
            "fragment_min_mz": 200.0, "fragment_max_mz": 1800.0,
            "peptide_min_mass": 500.0, "peptide_max_mass": 5000.0,
            "ion_kinds": ["b", "y"], "min_ion_index": 2,
            "static_mods": {"C": 57.0215},
            "variable_mods": variable,
            "max_variable_mods": 2,
            "decoy_tag": "rev_", "generate_decoys": True,
            "fasta": "REPLACED_AT_RUNTIME.fasta",
        },
        "precursor_tol": {"ppm": [-float(prec_ppm), float(prec_ppm)]},
        "fragment_tol": ({"da": [-SAGE_ITMS_FRAGMENT_DA, SAGE_ITMS_FRAGMENT_DA]} if ion_trap
                         else {"ppm": [-float(frag_ppm), float(frag_ppm)]}),
        "isotope_errors": [0, 1],
        "deisotope": not ion_trap,
        "chimera": acq.upper() == "DDA",
        "wide_window": acq.upper() == "DIA",
        "predict_rt": True,
        "min_peaks": 15, "max_peaks": 150, "min_matched_peaks": 4,
        "max_fragment_charge": 2,
        "quant": {"lfq": True},
        "output_directory": "sage_out",
    }
    # overrides: shallow-merge top-level keys
    for k, v in (overrides or {}).items():
        cfg[k] = v

    rationale = {"precursor_tol_ppm": tagged(prec_ppm, ppm_src)}
    if ion_trap:
        it = ("MS2 read partly in the ion trap (mixed ITMS/FTMS)" if ms2_analyzer == "mixed"
              else "MS2 read in the ion trap (ITMS)")
        rationale.update({
            "fragment_tol_da": tagged(SAGE_ITMS_FRAGMENT_DA, f"{it}: {SRC_SAGE_LOWRES}"),
            "deisotope": tagged(False, f"{it}: {SRC_SAGE_LOWRES}"),
            "bucket_size": tagged(SAGE_ITMS_BUCKET_SIZE, f"{it}: {SRC_SAGE_LOWRES} ('higher "
                                                          "values for low-res MS/MS')"),
        })
    else:
        rationale["fragment_tol_ppm"] = tagged(frag_ppm, ppm_src)
    rationale.update({
        "wide_window": tagged(cfg["wide_window"], f"{acq.upper()} acquisition"),
        "chimera": tagged(cfg["chimera"], f"{acq.upper()} acquisition"),
        "static_mods": tagged({"C": 57.0215}, "fixed carbamidomethyl (standard)"),
        "variable_mods": tagged(variable,
            "Ox(M) " + ("on (user requested)" if var_mods and "ox" in var_mods else "off (quant default)")
            + " + protein N-term acetyl"),
        "enzyme": tagged("trypsin/P, 2 missed cleavages", UNIV),
        "lfq": tagged(True, "label-free quantification"),
    })
    return json.dumps(cfg, indent=2) + "\n", rationale


def sage_fragment_mismatch(ms2_analyzer, sage_cfg):
    """("refuse" | "warn", why) when a Sage config's fragment window does not fit the MS2 analyzer
    the manifest records, else None. run_search.py calls it before anything is generated.

    SKILL.md passes --ms2-analyzer in step 4 (resolve_defaults.py -> the manifest) AND step 6b
    (this script -> the Sage cfg); only the second changes the search. Given in step 4 and
    forgotten in 6b, the manifest said ITMS while sage_config.json kept a +/-10 ppm fragment
    window -- which on an ion-trap spectrum matches nothing (a +0.25 Da synthetic spectrum: 0 PSMs
    at 10 ppm, the target at 0.4 Da; Sage 0.14.7, HIVE job 23990274)."""
    tol = (sage_cfg or {}).get("fragment_tol") or {}
    unit = "da" if "da" in tol else "ppm" if "ppm" in tol else None
    rerun = ("re-run estimate_params.py --engine sage with the same --ms2-analyzer that "
             "resolve_defaults.py was given, then run_search.py again")
    if ms2_in_ion_trap(ms2_analyzer) and unit == "ppm":
        return ("refuse", f"workflow.manifest.json records ms2_analyzer {ms2_analyzer} (MS2 read "
                f"in the ion trap{' for part of the run' if ms2_analyzer == 'mixed' else ''}), "
                f"but the Sage config's fragment_tol is {json.dumps(tol)} -- a ppm window matches "
                f"almost no ion-trap fragments. Fix: {rerun} (--ms2-analyzer {ms2_analyzer} "
                f"writes Sage's low-res window, fragment_tol {{\"da\": [-{SAGE_ITMS_FRAGMENT_DA}, "
                f"{SAGE_ITMS_FRAGMENT_DA}]}}).")
    # Refused too, not warned: estimate_params.py only writes a Da window for an ion-trap MS2, so a
    # Da cfg beside an FTMS manifest was built for other data. On Orbitrap spectra it is ~40x
    # Sage's documented high-res +/-10 ppm -- far more candidate matches for the scorer to reject.
    if ms2_analyzer == "FTMS" and unit == "da":
        return ("refuse", "workflow.manifest.json records ms2_analyzer FTMS (MS2 in the Orbitrap), "
                f"but the Sage config's fragment_tol is {json.dumps(tol)} -- the ion-trap window. "
                f"Fix: {rerun}; if the MS2 really is ion trap, re-run resolve_defaults.py with "
                "--ms2-analyzer ITMS instead.")
    if ms2_analyzer is None and unit == "da":
        return ("warn", f"the Sage config's fragment_tol is {json.dumps(tol)} (an ion-trap window), "
                "but workflow.manifest.json records no ms2_analyzer, so it cannot be checked "
                "against the data. If the MS2 is ion trap, re-run resolve_defaults.py with "
                "--ms2-analyzer ITMS so the record says so.")
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--engine", required=True, choices=["diann", "sage"])
    ap.add_argument("--acquisition", required=True, choices=["DIA", "DDA", "dia", "dda"])
    ap.add_argument("--instrument", default="")
    ap.add_argument("--var-mods", default="", help="comma list, e.g. 'ox' to add Ox(M)")
    ap.add_argument("--overrides", default="", help="JSON of fields to force (validated SOP)")
    add_resolution_args(ap)
    ap.add_argument("--from-mzml", default="",
                    help="read both resolutions straight out of this mzML")
    ap.add_argument("--precursor-mz-range", nargs=2, type=float, metavar=("LO", "HI"),
                    default=None,
                    help="ACQUIRED precursor m/z bounds, from detect_acquisition.py's "
                         "precursor_mz_range. Without it the range falls back to 380-980 "
                         "and is tagged FALLBACK -- which silently discards anything the "
                         "method acquired outside that window.")
    ap.add_argument("--fasta-meta", default="",
                    help="fetch_fasta.py's <fasta>.meta.json — supplies the contaminant "
                         "tag so DIA-NN excludes contaminants from quant")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    cont_tag = None
    if a.fasta_meta:
        try:
            with open(a.fasta_meta) as fh:
                cont_tag = json.load(fh).get("diann_cont_quant_exclude")
        except (OSError, json.JSONDecodeError) as e:
            sys.exit(f"--fasta-meta could not be read: {e}")

    overrides = {}
    if a.overrides:
        try:
            overrides = json.loads(a.overrides)
        except json.JSONDecodeError as e:
            sys.exit(f"--overrides is not valid JSON: {e}")

    r1, r2, analyzer = a.ms1_resolution, a.ms2_resolution, a.ms2_analyzer
    bad = resolution_args_error(r2, analyzer)
    if bad:
        sys.exit(f"[estimate_params] {bad}")
    res_src = resolution_source(a.resolution_source, r1, r2, analyzer)
    if a.from_mzml and not (r1 and r2):
        m1, m2 = read_mzml_resolution(a.from_mzml)
        if m1:      # an mzML with no MS:1000800 terms leaves the given values alone
            # an ion-trap MS2 has no resolving power of its own; read_mzml_resolution() would
            # hand back the MS1 value for it
            r1, r2, res_src = m1, (None if analyzer == "ITMS" else m2), "mzml"
            print(f"[estimate_params] read resolution from {os.path.basename(a.from_mzml)}: "
                  f"MS1 {int(r1):,}" + (f" / MS2 {int(r2):,}" if r2 else ""), file=sys.stderr)
    cls, ms1, ms2, label, src = classify_instrument(a.instrument, r1, r2, res_src, analyzer)
    ask = resolution_question(cls, a.engine, r1, r2)
    if ask:
        print(f"[estimate_params] NEEDS CONFIRMATION: {ask}", file=sys.stderr)
    var_mods = [v.strip().lower() for v in a.var_mods.split(",") if v.strip()]
    level_src = ((ppm_for_resolution(r1)[1], ppm_for_resolution(r2)[1])
                 if cls in ("orbitrap_measured", "orbitrap_untabled") else None)

    if a.engine == "diann":
        try:
            text, rationale = build_diann(a.acquisition, cls, ms1, ms2, label, src, var_mods,
                                          overrides, cont_tag, a.precursor_mz_range,
                                          level_src=level_src)
        except LoneMassAccOverride as e:
            # The message says "no cfg was written", so leave nothing at --out that contradicts
            # it. build_diann() raises BEFORE the write below, so THIS run wrote neither file --
            # but an earlier run of the same command may have left both, and the next step reads
            # them by path, not by mtime. A stale cfg with a stale mass_accuracy_plan beside a
            # refusal on stdout is the worst of both.
            for stale in (a.out, a.out + ".rationale.json"):
                try:
                    os.unlink(stale)
                except OSError:
                    pass
            sys.exit(f"[estimate_params] {e}")
    else:
        text, rationale = build_sage(a.acquisition, cls, var_mods, overrides, analyzer)

    with open(a.out, "w") as fh:
        fh.write(text)

    out_payload = {
        "engine": a.engine, "acquisition": a.acquisition.upper(),
        "instrument": a.instrument, "instrument_class": cls, "class_label": label,
        "resolution": resolution_record(r1, r2, res_src, analyzer),
        # an Orbitrap with no resolution, on an engine that uses it: ask before searching
        "needs_confirmation": bool(ask), "ask_user": ask,
        "mass_accuracy_source": src,
        # read by diann_parallel.mass_acc_measure_plan(): measure_with_diann is the ONLY way an
        # unpinned mass accuracy is accepted by the 5-step chain, and what makes the single-shot
        # search measure it; the documented levels are pinned as documented
        "mass_accuracy_plan": (rationale.get("mass_accuracy_plan") or {}).get("value"),
        "mass_accuracy_documented": (rationale.get("mass_accuracy_documented") or {}).get("value"),
        "params_file": os.path.abspath(a.out),
        "rationale": rationale,
        "note": "Every value is tagged with its provenance. Mass tolerances are the "
                "data-type-dependent settings; the rest are standard trypsin/LFQ defaults.",
    }
    # persist the rationale next to the params file so it travels into the
    # reproducibility bundle automatically (provenance.py picks up the sibling).
    with open(a.out + ".rationale.json", "w") as fh:
        json.dump(out_payload, fh, indent=2)
    print(json.dumps(out_payload, indent=2))


if __name__ == "__main__":
    main()
