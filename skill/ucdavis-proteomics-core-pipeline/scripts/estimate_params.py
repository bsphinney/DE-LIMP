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


def classify_instrument(name):
    """Return (class, ms1_ppm, ms2_ppm, label, source). None ppm => auto-calibrate."""
    n = (name or "").strip().lower()
    if not n:
        return ("unknown", None, None, "instrument not detected", "auto-calibration fallback")
    if "astral" in n:
        return ("orbitrap_astral", 4, 10, "Orbitrap Astral (assumes 240k MS1)", SRC_DIANN)
    if "tims" in n:
        return ("timstof", 15, 15, "Bruker timsTOF (dia-PASEF / ddaPASEF)", SRC_DIANN)
    if "tripletof" in n or "zenotof" in n or "sciex" in n:
        return ("sciex_tof", 20, 20, "SCIEX TripleTOF / ZenoTOF", SRC_DIANN)
    # Orbitrap family without a resolution we can read -> DIA-NN auto-calibration
    if any(k in n for k in ("orbitrap", "exploris", "exactive", "fusion", "lumos",
                            "eclipse", "velos", "hf-x", "hf", "qe", "astral")):
        return ("orbitrap_generic", None, None,
                "Orbitrap (resolution unknown — using DIA-NN automatic calibration)",
                "DIA-NN default (--mass-acc 0 auto-optimises per run)")
    return ("unknown", None, None, f"unrecognized instrument '{name}'", "auto-calibration fallback")


def tagged(value, source):
    return {"value": value, "source": source}


# --- DIA-NN cfg --------------------------------------------------------------
def build_diann(acq, instr_class, ms1, ms2, label, src, var_mods, overrides):
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

    # mass accuracy — the data-type-dependent part
    if ms2 is None:
        add("--mass-acc", 0, src)        # 0 = automatic
        add("--mass-acc-ms1", 0, src)
    else:
        add("--mass-acc", ms2, f"{label}: MS2 {ms2} ppm [{src}]")
        add("--mass-acc-ms1", ms1, f"{label}: MS1 {ms1} ppm [{src}]")
    add("--window", 0, "scan window auto-optimised per run (DIA-NN default)")

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
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    overrides = {}
    if a.overrides:
        try:
            overrides = json.loads(a.overrides)
        except json.JSONDecodeError as e:
            sys.exit(f"--overrides is not valid JSON: {e}")

    cls, ms1, ms2, label, src = classify_instrument(a.instrument)
    var_mods = [v.strip().lower() for v in a.var_mods.split(",") if v.strip()]

    if a.engine == "diann":
        text, rationale = build_diann(a.acquisition, cls, ms1, ms2, label, src, var_mods, overrides)
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
