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
import sys, os, json, argparse

# --- known-good mass accuracy by instrument class (DIA-NN README) ------------
ORBITRAP_RES_PPM = {240000: 4, 120000: 7, 60000: 10, 30000: 15}
SRC_DIANN = "DIA-NN README recommended settings (verified 2026-06)"


# DIA-NN's own Orbitrap resolution -> mass-accuracy table, verbatim from its README
# ("Changing default settings" §4). Resolution is the thing that actually sets the
# achievable tolerance, which is why a model name alone is not enough.
DIANN_ORBITRAP_PPM = {240000: 4, 120000: 7, 60000: 10, 30000: 15}
SRC_TABLE = "DIA-NN README, Orbitrap resolution->accuracy table"


def ppm_for_resolution(res):
    """Map an Orbitrap resolving power to mass accuracy (ppm), per DIA-NN's table.

    Exact tiers use the documented value. Anything between/outside is interpolated on
    a log-log fit of the table and TAGGED as extrapolated -- DIA-NN documents 30k-240k
    only, and real acquisitions sit outside it (an Exploris 480 DIA method commonly
    runs MS2 at 15k). Silently emitting a number as if it were documented is exactly
    the fabricated-default failure mode; the caller renders the tag.
    """
    if not res or res <= 0:
        return None, None
    res = float(res)
    if int(res) in DIANN_ORBITRAP_PPM:
        return DIANN_ORBITRAP_PPM[int(res)], SRC_TABLE
    import math
    pts = sorted(DIANN_ORBITRAP_PPM.items())
    # log(ppm) is close to linear in log(resolution) across the documented tiers
    (r1, p1), (r2, p2) = pts[0], pts[-1]
    slope = (math.log(p2) - math.log(p1)) / (math.log(r2) - math.log(r1))
    ppm = math.exp(math.log(p1) + slope * (math.log(res) - math.log(r1)))
    ppm = round(ppm, 1)
    inside = pts[0][0] <= res <= pts[-1][0]
    tag = (f"{SRC_TABLE}, interpolated for {int(res):,} resolution" if inside else
           f"{SRC_TABLE}, EXTRAPOLATED — {int(res):,} is outside the documented "
           f"30k-240k range; verify against a real run")
    return ppm, tag


def classify_instrument(name, ms1_res=None, ms2_res=None):
    """Return (class, ms1_ppm, ms2_ppm, label, source). None ppm => auto-calibrate.

    DIA-NN asks for these to be FIXED rather than auto-optimised, and not only for
    speed: "This optimisation is inherently noisy: even replicate injections may not
    produce identical results, and therefore the analysis results will depend on which
    run is first in the list." Auto therefore makes a result depend on FILE ORDER,
    which also blocks the 5-step parallel chain (steps 3/5 reuse .quant files).
    """
    n = (name or "").strip().lower()
    # Measured resolution beats any model-name guess.
    if ms1_res or ms2_res:
        p1, s1 = ppm_for_resolution(ms1_res)
        p2, s2 = ppm_for_resolution(ms2_res)
        if p1 and p2:
            return ("orbitrap_measured", p1, p2,
                    f"Orbitrap, MS1 {int(ms1_res):,} / MS2 {int(ms2_res):,} resolution "
                    f"read from the data", s2 if "EXTRAPOL" in (s2 or "") else s1)
    if not n:
        return ("unknown", None, None, "instrument not detected", "auto-calibration fallback")
    if "astral" in n:
        return ("orbitrap_astral", 4, 10, "Orbitrap Astral (assumes 240k MS1)", SRC_DIANN)
    if "tims" in n:
        return ("timstof", 15, 15, "Bruker timsTOF (dia-PASEF / ddaPASEF)", SRC_DIANN)
    if "tripletof" in n or "zenotof" in n or "sciex" in n:
        return ("sciex_tof", 20, 20, "SCIEX TripleTOF / ZenoTOF", SRC_DIANN)
    # Orbitrap family with no resolution to work from -> DIA-NN auto-calibration.
    # Prefer passing --ms1-resolution/--ms2-resolution (read_mzml_resolution() gets
    # them straight out of the mzML) so the run is reproducible and parallel-capable.
    if any(k in n for k in ("orbitrap", "exploris", "exactive", "fusion", "lumos",
                            "eclipse", "velos", "hf-x", "hf", "qe", "astral")):
        return ("orbitrap_generic", None, None,
                "Orbitrap (resolution unknown — using DIA-NN automatic calibration; "
                "pass --ms1-resolution/--ms2-resolution to pin it)",
                "DIA-NN automatic calibration (mass-accuracy flags omitted)")
    return ("unknown", None, None, f"unrecognized instrument '{name}'", "auto-calibration fallback")


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
                cont_tag=None):
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
    add("--fasta-search", True, "library-free (predicted spectral library)")
    add("--gen-spec-lib", True, UNIV)
    add("--predictor", True, "deep-learning predictor for library-free DIA")
    add("--reanalyse", True, "MBR across the run set")
    add("--rt-profiling", True, UNIV)
    add("--cut", "K*,R*", "trypsin/P")
    add("--missed-cleavages", 1, "DIA-NN default")
    add("--min-pep-len", 7, UNIV)
    add("--max-pep-len", 30, UNIV)
    add("--min-pr-mz", 380, UNIV)
    add("--max-pr-mz", 980, UNIV)
    add("--min-pr-charge", 2, UNIV)
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
    # diann_parallel.parallel_safe() treats an absent flag exactly like 0 (both are
    # "not fixed"), so the parallel-chain gate is unaffected.
    if ms2 is None:
        auto_src = (f"{src}; flags omitted so DIA-NN auto-calibrates per run "
                    "(--mass-acc 0 would pin the tolerance at 0 ppm -> 0 IDs)")
        add("--mass-acc", "auto (flag omitted)", auto_src, render=False)
        add("--mass-acc-ms1", "auto (flag omitted)", auto_src, render=False)
    else:
        add("--mass-acc", ms2, f"{label}: MS2 {ms2} ppm [{src}]")
        add("--mass-acc-ms1", ms1, f"{label}: MS1 {ms1} ppm [{src}]")
    # --window 0 is likewise rejected ("scan window radius should be a positive
    # integer"); omitting it lets DIA-NN set the radius from the observed peak width.
    add("--window", "auto (flag omitted)",
        "scan window radius auto-optimised per run from observed peak width",
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

    if a.engine == "diann":
        text, rationale = build_diann(a.acquisition, cls, ms1, ms2, label, src, var_mods,
                                      overrides, cont_tag)
    else:
        text, rationale = build_sage(a.acquisition, cls, var_mods, overrides)

    with open(a.out, "w") as fh:
        fh.write(text)

    out_payload = {
        "engine": a.engine, "acquisition": a.acquisition.upper(),
        "instrument": a.instrument, "instrument_class": cls, "class_label": label,
        "mass_accuracy_source": src,
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
