#!/usr/bin/env python3
"""
make_presets.py -- generate the engine preset file for THIS run, from the data type.

FragPipe and Radiant are not configured by command-line flags the way DIA-NN and
Sage are: each is driven by a whole config file. So "estimate parameters from the
data type" means, for them, "write the right config file before the search".

WHY TEMPLATE-BASED, NOT GENERATED FROM SCRATCH
----------------------------------------------
FragPipe's two DIA presets differ in 20 keys, not one. Diffed on 2026-08-14:

  diatracer.run-diatracer  true / false        msfragger.use_topN_peaks   300 / 1000
  ionquant.mbr             1 / 0               percolator.min-prob        0.5 / 0.7
  ionquant.minisotopes     2 / 1               msfragger.precursor_mass   +/-10 / +/-20
  allowed_missed_cleavage  2 / 1               precursor-charge-lo        1 / 2
  ...plus the var-mod table (SILAC rows only in the diaPASEF one), ptmprophet,
  protein-prophet, phi-report and four ptmshepherd keys.

Flipping diatracer.run-diatracer on the wrong template does NOT convert one route
into the other -- it produces a config no vendor ever validated. So we start from
the vendor's own preset for the detected data type and patch only what is
genuinely run-specific (FASTA, threads) or explicitly overridden.

Editing is LINE-BASED on purpose: these are Java .properties files whose values
carry escapes (\\=, \\:, embedded HTML). A parse/re-serialise round-trip would
silently rewrite that escaping across ~400 lines. We touch only the lines we mean
to change and copy every other byte through unaltered.

The two run-specific FragPipe keys we append are absent from the shipped presets
by design, and both were verified against FragPipe source on 2026-08-14 rather
than assumed:
  database.db-path   PropsFile.java get/setProperty; present in jobs/example_job
                     .workflow; PropertiesUtils.java: "Move database.db-path to
                     the top. Other run specific parameters are provided through
                     command line."
  workflow.threads   Fragpipe.java sets it from nThreadsHeadlessOnly, and
                     defaults it when the key is null.

Usage
  python3 make_presets.py --engine fragpipe --instrument "timsTOF HT" \
      --acquisition DIA --fasta /path/db.fasta --threads 32 --out run/fp.workflow

  python3 make_presets.py --engine radiant --instrument "Orbitrap Astral" \
      --acquisition DIA --out run/run.radiantConfig

Prints a JSON provenance record naming every key changed from the vendor default
and why. Writes the same record next to the output as <out>.provenance.json.
"""
import argparse
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# ONE instrument -> ppm table for the whole skill. Do not restate it here.
from estimate_params import classify_instrument  # noqa: E402

PRESETS = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "presets")

FRAGPIPE_TEMPLATES = {
    "timstof": ("fragpipe/DIA_SpecLib_Quant_diaPASEF.workflow",
                "FragPipe 24.0 diaPASEF preset (diaTracer ON) -- Bruker dia-PASEF"),
    "orbitrap": ("fragpipe/DIA_SpecLib_Quant.workflow",
                 "FragPipe 24.0 DIA_SpecLib_Quant preset (diaTracer OFF) -- MSFragger-DIA direct"),
}


def _esc(v):
    """Escape a value for a Java .properties file. Only backslashes need it in
    a value position; ':' and '=' are literal after the first separator."""
    return str(v).replace("\\", "\\\\")


def patch_properties(text, changes):
    """Set each key in `changes`, in place, preserving every other line verbatim.
    Keys absent from the template are appended (FragPipe accepts extra keys)."""
    lines = text.splitlines(keepends=True)
    seen = set()
    for i, line in enumerate(lines):
        s = line.lstrip()
        if not s or s.startswith("#") or "=" not in line:
            continue
        key = line.split("=", 1)[0].strip()
        if key in changes:
            eol = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
            lines[i] = f"{key}={_esc(changes[key])}{eol}"
            seen.add(key)
    tail = "".join(f"{k}={_esc(v)}\n" for k, v in changes.items() if k not in seen)
    out = "".join(lines)
    if tail:
        if out and not out.endswith("\n"):
            out += "\n"
        out += tail
    return out, sorted(seen), sorted(k for k in changes if k not in seen)


def patch_ini(text, changes):
    """Set `section.key = value` in an INI file, preserving layout and comments."""
    lines = text.splitlines(keepends=True)
    section = None
    applied = []
    for i, line in enumerate(lines):
        s = line.strip()
        m = re.match(r"^\[(.+)\]$", s)
        if m:
            section = m.group(1)
            continue
        if not s or s.startswith(("#", ";")) or "=" not in line:
            continue
        key = line.split("=", 1)[0].strip()
        want = changes.get((section, key))
        if want is not None:
            eol = "\r\n" if line.endswith("\r\n") else "\n" if line.endswith("\n") else ""
            lines[i] = f"{key} = {want}{eol}"
            applied.append(f"[{section}] {key}")
    return "".join(lines), applied


def build_fragpipe(args, cls, ms1, ms2, label, src, prov):
    if cls == "timstof":
        family = "timstof"
    elif cls.startswith("orbitrap"):
        family = "orbitrap"
    else:
        sys.exit(f"make_presets: FragPipe has no validated DIA preset for instrument class "
                 f"'{cls}' ({label}). FragPipe ships presets for Bruker dia-PASEF and Thermo "
                 f"Orbitrap only -- pick one explicitly with --template if you know better.")
    if args.template:
        family = args.template
    rel, desc = FRAGPIPE_TEMPLATES[family]
    path = os.path.join(PRESETS, rel)
    text = open(path).read()

    changes = {}
    why = {}
    # Guard rather than trust: the whole point of picking a template by data type
    # is that diaTracer is ON for Bruker and OFF for Thermo.
    want_dt = "true" if family == "timstof" else "false"
    changes["diatracer.run-diatracer"] = want_dt
    why["diatracer.run-diatracer"] = (
        f"{want_dt} -- diaTracer traces Bruker dia-PASEF .d into pseudo-MS/MS; it has no role "
        f"on Thermo DIA, where MSFragger-DIA searches the spectra directly")
    if args.fasta:
        changes["database.db-path"] = os.path.abspath(args.fasta)
        why["database.db-path"] = "this run's FASTA"
    if args.threads:
        changes["workflow.threads"] = args.threads
        why["workflow.threads"] = "this run's CPU allocation"

    # Mass tolerances: the vendor already tuned these PER DATA TYPE (that is what
    # the template selection buys us). Only override when the caller supplies a
    # measured/SOP number -- do NOT silently graft DIA-NN's ppm table onto
    # MSFragger, which is a different matcher with different semantics.
    if args.ms1_ppm is not None:
        changes["msfragger.precursor_true_tolerance"] = args.ms1_ppm
        changes["msfragger.precursor_true_units"] = 1          # 1 = ppm
        why["msfragger.precursor_true_tolerance"] = f"{args.ms1_ppm} ppm (caller override)"
    if args.ms2_ppm is not None:
        changes["msfragger.fragment_mass_tolerance"] = args.ms2_ppm
        changes["msfragger.fragment_mass_units"] = 1
        why["msfragger.fragment_mass_tolerance"] = f"{args.ms2_ppm} ppm (caller override)"

    out, patched, added = patch_properties(text, changes)
    prov.update({
        "engine": "fragpipe",
        "template": rel,
        "template_note": desc,
        "template_is_vendor_default": True,
        "keys_changed": {k: why[k] for k in why},
        "keys_patched_in_place": patched,
        "keys_appended": added,
        "tolerances": (
            "vendor preset values kept -- FragPipe tunes these per data type and no "
            "override was supplied" if args.ms1_ppm is None and args.ms2_ppm is None
            else "caller-supplied override applied"),
    })
    return out


def build_radiant(args, cls, ms1, ms2, label, src, prov):
    if cls == "timstof":
        sys.exit("make_presets: Radiant/Fulcrum cannot read Bruker .d at all -- its search "
                 "backend takes mzML or Parquet only. Use the DIA-NN or FragPipe+diaTracer "
                 "route for timsTOF data.")
    path = os.path.join(PRESETS, "radiant/default.radiantConfig")
    text = open(path).read()

    changes = {}
    why = {}
    m1 = args.ms1_ppm if args.ms1_ppm is not None else ms1
    m2 = args.ms2_ppm if args.ms2_ppm is not None else ms2
    # The shipped default is 20/20 ppm for every Orbitrap, which is visibly wrong
    # for an Astral's narrow windows. Narrow it when we know the instrument.
    # Tagged DERIVED: the numbers come from DIA-NN's documented instrument table,
    # and applying them to Radiant's extraction width is our inference, not Seer's.
    if m1 is not None:
        changes[("MS1Params", "ms1ExtractionWidthPPM")] = float(m1)
        why["[MS1Params] ms1ExtractionWidthPPM"] = (
            f"{m1} ppm for {label} [DERIVED from DIA-NN's instrument table; Seer ships a flat "
            f"20.0 for every Orbitrap, which is too wide for narrow-window data]"
            if args.ms1_ppm is None else f"{m1} ppm (caller override)")
    if m2 is not None:
        changes[("MS2Params", "ms2ExtractionWidthPPM")] = float(m2)
        why["[MS2Params] ms2ExtractionWidthPPM"] = (
            f"{m2} ppm for {label} [DERIVED as above]"
            if args.ms2_ppm is None else f"{m2} ppm (caller override)")
    if args.threads:
        changes[("General", "threadCount")] = int(args.threads)
        why["[General] threadCount"] = "this run's CPU allocation (0 = all cores)"

    out, applied = patch_ini(text, changes)
    prov.update({
        "engine": "radiant",
        "template": "radiant/default.radiantConfig",
        "template_note": "radiant_fulcrum_workflow 0.3.3 shipped default, copied verbatim",
        "template_is_vendor_default": True,
        "keys_changed": why,
        "keys_patched_in_place": applied,
        "tolerances": ("vendor 20/20 ppm kept -- instrument not recognised"
                       if m1 is None else "narrowed from the instrument"),
        "licence": "Apache-2.0 + Commons Clause + grant-back; restricts fee-for-service use",
    })
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--engine", required=True, choices=["fragpipe", "radiant"],
                    help="DIA-NN and Sage take flags/JSON from estimate_params.py instead")
    ap.add_argument("--instrument", default="", help="detected instrument name")
    ap.add_argument("--acquisition", default="DIA", choices=["DIA", "DDA"])
    ap.add_argument("--ms1-res", type=float, default=None)
    ap.add_argument("--ms2-res", type=float, default=None)
    ap.add_argument("--ms1-ppm", type=float, default=None, help="override; else vendor/derived")
    ap.add_argument("--ms2-ppm", type=float, default=None, help="override; else vendor/derived")
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--template", choices=["timstof", "orbitrap"],
                    help="force a FragPipe template instead of choosing by instrument")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    if args.acquisition == "DDA":
        sys.exit("make_presets: these presets are the DIA routes. DDA goes through Sage "
                 "(estimate_params.py --engine sage).")

    cls, ms1, ms2, label, src = classify_instrument(args.instrument, args.ms1_res, args.ms2_res)
    prov = {"instrument": args.instrument or None, "instrument_class": cls,
            "instrument_label": label, "acquisition": args.acquisition,
            "ppm_source": src, "output": os.path.abspath(args.out)}

    text = (build_fragpipe if args.engine == "fragpipe" else build_radiant)(
        args, cls, ms1, ms2, label, src, prov)

    d = os.path.dirname(os.path.abspath(args.out))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(text)
    with open(args.out + ".provenance.json", "w") as fh:
        json.dump(prov, fh, indent=2)
    print(json.dumps(prov, indent=2))


if __name__ == "__main__":
    main()
