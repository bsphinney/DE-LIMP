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
    Orbitrap level with no documented tier (resolution unknown, or outside 30k-240k):
                                 MEASURED with DIA-NN on representative runs
                                 before the search; a level that has a tier keeps it
    Unidentified instrument:     automatic calibration (DIA-NN's own default)
Sage's docs give NO instrument-specific tolerances, so Sage ppm windows are
DERIVED from the same per-instrument logic and tagged as such.

Usage:
  python3 estimate_params.py --engine diann --acquisition DIA \
      --instrument "Orbitrap Astral" --out diann.cfg
  python3 estimate_params.py --engine sage --acquisition DDA \
      --instrument "timsTOF Pro" --out sage_config.json [--var-mods ox]

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


def instrument_ppm_summary():
    """One line of the table, MS1/MS2 ppm, for remediation text."""
    t = DIANN_INSTRUMENT_PPM
    orb = ", ".join(f"{r // 1000}k->{p}" for r, p in sorted(DIANN_ORBITRAP_PPM.items(), reverse=True))
    return (f"MS1/MS2 ppm: timsTOF {t['timstof'][0]}/{t['timstof'][1]}, "
            f"Astral {t['orbitrap_astral'][0]}/{t['orbitrap_astral'][1]}, "
            f"SCIEX {t['sciex_tof'][0]}/{t['sciex_tof'][1]}, "
            f"Orbitrap by resolution ({orb}; a level outside 30k-240k, or of unknown resolution, "
            f"is not extrapolated: estimate_params.py plans a DIA-NN measurement of it)")


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
MEASURE_CLASSES = ("orbitrap_generic", "orbitrap_untabled")

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


def classify_instrument(name, ms1_res=None, ms2_res=None):
    """Return (class, ms1_ppm, ms2_ppm, label, source). None ppm => not pinned.

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
    # Measured resolution beats any model-name guess.
    if ms1_res or ms2_res:
        p1, s1 = ppm_for_resolution(ms1_res)
        p2, s2 = ppm_for_resolution(ms2_res)
        if p1 and p2:
            src = s1 if s1 == s2 else f"MS1: {s1}; MS2: {s2}"
            return ("orbitrap_measured", p1, p2,
                    f"Orbitrap, MS1 {int(ms1_res):,} / MS2 {int(ms2_res):,} resolution "
                    f"read from the data", src)
        if ms1_res and ms2_res:
            src = "; ".join(
                f"{lvl} {int(res):,}: {s}" if p is not None else
                f"{lvl} {int(res):,}: no documented DIA-NN value ({s})"
                for lvl, res, p, s in (("MS1", ms1_res, p1, s1), ("MS2", ms2_res, p2, s2)))
            return ("orbitrap_untabled", p1, p2,
                    f"Orbitrap, MS1 {int(ms1_res):,} / MS2 {int(ms2_res):,} resolution "
                    f"read from the data", src)
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
        auto_src = (f"{src}; flags omitted so DIA-NN auto-calibrates per run "
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
        PLAN_AUTO: "instrument not identified; DIA-NN optimises it itself (not parallel-safe)",
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
    if instr_class == "orbitrap_generic":  return (10, 10, "high-res Orbitrap default (derived)")
    return (20, 20, "safe high-res default (instrument not identified)")


def build_sage(acq, instr_class, var_mods, overrides):
    prec_ppm, frag_ppm, ppm_src = sage_ppm(instr_class)
    UNIV = "universal trypsin/LFQ default"
    variable = {"M": [15.9949]} if (var_mods and "ox" in var_mods) else {}
    variable["["] = [42.0106]  # protein N-term acetyl is a common, cheap variable mod
    cfg = {
        "database": {
            "bucket_size": 8192,
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
        "fragment_tol":  {"ppm": [-float(frag_ppm), float(frag_ppm)]},
        "isotope_errors": [0, 1],
        "deisotope": True,
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

    rationale = {
        "precursor_tol_ppm": tagged(prec_ppm, ppm_src),
        "fragment_tol_ppm": tagged(frag_ppm, ppm_src),
        "wide_window": tagged(cfg["wide_window"], f"{acq.upper()} acquisition"),
        "chimera": tagged(cfg["chimera"], f"{acq.upper()} acquisition"),
        "static_mods": tagged({"C": 57.0215}, "fixed carbamidomethyl (standard)"),
        "variable_mods": tagged(variable,
            "Ox(M) " + ("on (user requested)" if var_mods and "ox" in var_mods else "off (quant default)")
            + " + protein N-term acetyl"),
        "enzyme": tagged("trypsin/P, 2 missed cleavages", UNIV),
        "lfq": tagged(True, "label-free quantification"),
    }
    return json.dumps(cfg, indent=2) + "\n", rationale


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--engine", required=True, choices=["diann", "sage"])
    ap.add_argument("--acquisition", required=True, choices=["DIA", "DDA", "dia", "dda"])
    ap.add_argument("--instrument", default="")
    ap.add_argument("--var-mods", default="", help="comma list, e.g. 'ox' to add Ox(M)")
    ap.add_argument("--overrides", default="", help="JSON of fields to force (validated SOP)")
    ap.add_argument("--ms1-resolution", type=float, default=0,
                    help="Orbitrap MS1 resolving power; maps to ppm via DIA-NN's table")
    ap.add_argument("--ms2-resolution", type=float, default=0,
                    help="Orbitrap MS2 resolving power; maps to ppm via DIA-NN's table")
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

    r1, r2 = a.ms1_resolution, a.ms2_resolution
    if a.from_mzml and not (r1 and r2):
        r1, r2 = read_mzml_resolution(a.from_mzml)
        if r1:
            print(f"[estimate_params] read resolution from {os.path.basename(a.from_mzml)}: "
                  f"MS1 {int(r1):,} / MS2 {int(r2):,}", file=sys.stderr)
    cls, ms1, ms2, label, src = classify_instrument(a.instrument, r1, r2)
    var_mods = [v.strip().lower() for v in a.var_mods.split(",") if v.strip()]
    level_src = ((ppm_for_resolution(r1)[1], ppm_for_resolution(r2)[1])
                 if cls in ("orbitrap_measured", "orbitrap_untabled") else None)

    if a.engine == "diann":
        try:
            text, rationale = build_diann(a.acquisition, cls, ms1, ms2, label, src, var_mods,
                                          overrides, cont_tag, a.precursor_mz_range,
                                          level_src=level_src)
        except LoneMassAccOverride as e:
            sys.exit(f"[estimate_params] {e}")
    else:
        text, rationale = build_sage(a.acquisition, cls, var_mods, overrides)

    with open(a.out, "w") as fh:
        fh.write(text)

    out_payload = {
        "engine": a.engine, "acquisition": a.acquisition.upper(),
        "instrument": a.instrument, "instrument_class": cls, "class_label": label,
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
