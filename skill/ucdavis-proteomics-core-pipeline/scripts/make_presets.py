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

FRAGPIPE DDA (2.15)
------------------
--acquisition DDA starts from one of FragPipe 24.0's own DDA workflows, chosen by the digest and
by whether amounts are wanted (choose_dda_template): LFQ-MBR (an enzyme rule; IonQuant label-free
quantification with match-between-runs), Basic-Search (--ids-only), Nonspecific-peptidome or
Nonspecific-HLA (--enzyme none or pepsin; label-free quantification switched on only when
amounts are wanted). presets/fragpipe/SOURCES.json records each file's origin and sha256. The
DDA flags are estimate_params.py's (add_dda_args: --var-mods, --ids-only, --cys-alkylation,
--enzyme, --overrides ...), with the same words and refusals. An enzymatic template searches the
skill's DDA digest, shared with Sage (estimate_params.DDA_DIGEST, Ruling 57.1: trypsin not before
P, 2 missed cleavages, 7-30), or the --enzyme rule. NO database.db-path and NO workflow.threads
are written: run_search.py sets them at search time. A rule the skill
states (--enzyme, or the rule keys in --overrides) is written with the name `nonspecific` and NO
enzyme dropdown key: FragPipe applies the dropdown over the rule when it loads the workflow
(cleavage_rules.FRAGPIPE_DROPDOWN_KEY; C0, 2026-11-04). An ion-trap MS2 takes FragPipe's own
ion-trap setting (ION_TRAP_KEYS, from its TMT10-MS3 workflow); a unit-resolution instrument
(Thermo Stellar) its Stellar-DDA setting (UNIT_RESOLUTION_KEYS).

Usage
  python3 make_presets.py --engine fragpipe --instrument "timsTOF HT" \
      --acquisition DIA --fasta /path/db.fasta --threads 32 --out run/fp.workflow

  python3 make_presets.py --engine fragpipe --instrument "timsTOF HT" --acquisition DDA \
      [--ids-only | --directlfq] [--enzyme trypsin|none|...] [--var-mods ox] \
      [--cys-alkylation cam] \
      --out run/fragpipe.workflow

  python3 make_presets.py --engine radiant --instrument "Orbitrap Astral" \
      --acquisition DIA --out run/run.radiantConfig

Prints a JSON provenance record naming every key changed from the vendor default
and why. Writes the same record next to the output as <out>.provenance.json.
"""
import argparse
import hashlib
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# ONE instrument -> ppm table for the whole skill. Do not restate it here.
from estimate_params import (classify_instrument, add_resolution_args,  # noqa: E402
                             resolution_args_error, resolution_source)
# ONE map of what each alkylating agent leaves on cysteine (the DIA-NN and Sage routes use it too)
from estimate_params import (CYS_REAGENTS, cys_alkylation_record, cys_diann_flags)  # noqa: E402
# ONE set of DDA flags and their shared settings, with estimate_params.py (Sage)
from estimate_params import (add_dda_args, dda_args, DdaArgsError, refuse_out,  # noqa: E402
                             DDA_DIGEST, DDA_MAX_VAR_MODS, PROTEIN_NTERM_ACETYL, VAR_MODS,
                             var_mods_choice,
                             ms2_in_ion_trap, unit_resolution, unit_resolution_amounts_refused,
                             enzyme_source, override_source,
                             override_record, record_override_setter, ENZYME_FLAG, SEMI_FLAG, IDS_ONLY_FLAG,
                             DIRECTLFQ_FLAG,
                             FRAGPIPE_C_ROW, cys_searched, cys_record_as_written)
import cleavage_rules  # noqa: E402  the --enzyme rule, the dropdown rule, the digest record
import fragpipe_params as fpp  # noqa: E402  THE reader of a FragPipe workflow

# MSFragger's fixed-modification table names cysteine by this row; the shipped DIA presets fix it
# at 57.02146 (carbamidomethyl).
FRAGPIPE_FIX_C = re.compile(r"(?<![0-9.\-])-?[0-9.]+(,C \(cysteine\),)")

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
    # Cysteine: the vendor presets fix carbamidomethyl (57.02146) on C. Another prep's
    # modification (or none) replaces that one row of MSFragger's fixed-mod table; the rest of the
    # table is copied through. The library FragPipe builds -- and so its DIA-NN step -- carries it.
    cys = prov["cys_alkylation"]
    if cys["value"] != "cam":
        key = "msfragger.table.fix-mods"
        line = next((ln for ln in text.splitlines() if ln.split("=", 1)[0].strip() == key), None)
        if line is None or not FRAGPIPE_FIX_C.search(line):
            sys.exit(f"make_presets: {rel} has no cysteine row in {key}, so --cys-alkylation "
                     f"{cys['value']} cannot be applied -- the template changed; check it.")
        mass = 0.0 if cys["mass"] is None else cys["mass"]
        changes[key] = FRAGPIPE_FIX_C.sub(lambda m: f"{mass}{m.group(1)}",
                                          line.split("=", 1)[1].strip(), count=1)
        why[key] = (f"C (cysteine) fixed at {mass} -- {cys['searched']}: {cys['source']} "
                    "(the vendor preset fixes carbamidomethyl, 57.02146)")

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


# ------------------------------------------------------------------ FragPipe DDA (2.15) -------
# FragPipe 24.0's own DDA workflows (presets/fragpipe/SOURCES.json: byte for byte the FragPipe
# git tag 24.0 files). The template is the vendor's validated starting point for the digest and
# the scope; the skill patches only the lines its record names.
DDA_TEMPLATES = {
    "lfq": ("fragpipe/LFQ-MBR.workflow", "LFQ-MBR",
            "FragPipe 24.0 LFQ-MBR workflow: closed search, IonQuant label-free quantification "
            "with match-between-runs"),
    "ids": ("fragpipe/Basic-Search.workflow", "Basic-Search",
            "FragPipe 24.0 Basic-Search workflow: closed search, identifications only"),
    "peptidome": ("fragpipe/Nonspecific-peptidome.workflow", "Nonspecific-peptidome",
                  "FragPipe 24.0 Nonspecific-peptidome workflow: non-specific search, no "
                  "quantification as shipped"),
    "hla": ("fragpipe/Nonspecific-HLA.workflow", "Nonspecific-HLA",
            "FragPipe 24.0 Nonspecific-HLA workflow: non-specific search for HLA peptides, no "
            "quantification as shipped"),
}
NONSPECIFIC_WORKFLOWS = ("peptidome", "hla")
NONSPECIFIC_DEFAULT = "peptidome"

# MS2 read in the ion trap: FragPipe 24.0's own ion-trap setting, from its TMT10-MS3 and TMT16-MS3
# workflows ("identification from low mass accuracy MS2 (ion trap)"; tag 24.0, read 2026-11-04 --
# C0 fact 6). Against TMT10, TMT10-MS3 changes exactly the fragment tolerance (20 ppm -> 0.6 Da)
# and the isotope errors; the precursor window stays +/-20 ppm (an Orbitrap MS1). FragPipe ships
# no label-free workflow with a Da fragment tolerance, so this is ADAPTED from those workflows,
# and the record says so (Ruling 90). A test reads each value against TMT10-MS3.workflow.
ION_TRAP_KEYS = {
    "msfragger.fragment_mass_tolerance": "0.6", "msfragger.fragment_mass_units": "0",
    "msfragger.isotope_error": "0/1/2/3",
    "msfragger.precursor_mass_lower": "-20", "msfragger.precursor_mass_upper": "20",
    "msfragger.precursor_mass_units": "1",
}
ION_TRAP_SOURCE = ("FragPipe 24.0's own ion-trap MS2 setting, adapted from its TMT10-MS3 / "
                   "TMT16-MS3 workflows ('identification from low mass accuracy MS2 (ion trap)'): "
                   "fragment 0.6 Da, precursor +/-20 ppm, isotope errors 0/1/2/3. FragPipe ships "
                   "no label-free workflow with an ion-trap fragment tolerance")
# A unit-resolution instrument (MS1 and MS2; estimate_params.unit_resolution: a Thermo Stellar):
# FragPipe 24.0's Stellar-DDA workflow ("Stellar data (unit resolution in MS1 and MS2) acquired in
# DDA mode"; tag 24.0). Every msfragger key it changes from LFQ-MBR, except its digest (missed
# cleavages 1) and its variable modifications, which the skill's flags set. A test reads each
# value against Stellar-DDA.workflow.
UNIT_RESOLUTION_KEYS = {
    "msfragger.calibrate_mass": "0", "msfragger.deisotope": "0", "msfragger.deneutralloss": "0",
    "msfragger.fragment_ion_series": "y", "msfragger.fragment_mass_tolerance": "0.3",
    "msfragger.fragment_mass_units": "0", "msfragger.isotope_error": "0",
    "msfragger.max_fragment_charge": "1", "msfragger.minimum_ratio": "0.00",
    "msfragger.precursor_mass_lower": "-1.5", "msfragger.precursor_mass_mode": "isolated",
    "msfragger.precursor_mass_units": "0", "msfragger.precursor_mass_upper": "1.5",
    "msfragger.precursor_true_tolerance": "1", "msfragger.precursor_true_units": "0",
}
UNIT_RESOLUTION_SOURCE = ("FragPipe 24.0's own unit-resolution setting, adapted from its "
                          "Stellar-DDA workflow (unit resolution in MS1 and MS2): fragment 0.3 Da, "
                          "precursor +/-1.5 Da, no mass calibration or deisotoping, y ions")

# The short name of each adapted tolerance setting, by the record's tolerances.mode: the preset's
# change notes and the Methods (make_methods.py) both say it from here
TOLERANCE_METHODS = {
    "unit_resolution": "FragPipe 24.0's unit-resolution setting, adapted from its Stellar-DDA "
                       "workflow",
    "ion_trap": "FragPipe 24.0's ion-trap MS2 setting, adapted from its TMT10-MS3 workflow",
}

# Overrides the preset refuses: each is set elsewhere, or would undo what the record says
DIRECTLFQ_KEY = "ionquant.maxlfq"
FRAGPIPE_OVERRIDE_REFUSED = {
    "database.db-path": "run_search.py writes the database into the run's own copy of the "
                        "workflow at search time (the FASTA searched, with its decoys)",
    cleavage_rules.FRAGPIPE_DROPDOWN_KEY: "FragPipe applies the enzyme dropdown over the rule "
                                          "when it loads the workflow (C0, 2026-11-04): give the "
                                          "rule's own keys (msfragger.search_enzyme_cut_1, ...) "
                                          "or --enzyme",
    "workflow.threads": "run_search.py passes the job's CPUs (--threads) at search time",
}
DATABASE_NOTE = ("not written: run_search.py sets database.db-path in the run's own copy of this "
                 "workflow at search time (fragpipe_workflow_for_run), from its --fasta")
THREADS_NOTE = ("not written: run_search.py passes the job's CPUs to FragPipe (--threads) at "
                "search time")
# The default digest's tag (no --enzyme): the skill's, not FragPipe's vendor stricttrypsin
DDA_DIGEST_DEFAULT = ("DEFAULT, not user-confirmed: the skill's DDA digest, trypsin (cut after K/R, "
                      "not before P), shared with Sage (estimate_params.DDA_DIGEST; Ruling 57.1); "
                      f"pass {ENZYME_FLAG} with the enzyme(s) added to the samples to set it")
# A MSFragger variable-modification row: mass, site, on, max per peptide
VAR_ROW = re.compile(r"^\s*(-?[0-9.]+),([^,]*),(true|false),(-?\d+)\s*$")
# the row --var-mods deam adds: FragPipe 24.0's own N/Q deamidation row (citrullination.workflow
# @ tag 24.0, 0.984016,NQ,true,2), VAR_MODS' mass
DEAM_ROW = f"{VAR_MODS['deam']['mass']},NQ,true,{DDA_MAX_VAR_MODS}"


def _same(a, b, tol=0.001):
    try:
        return abs(float(a) - float(b)) < tol
    except (TypeError, ValueError):
        return False


def delete_properties(text, keys):
    """`text` without the lines that set `keys` -> (text, [keys deleted])."""
    out, gone = [], []
    for ln in text.splitlines(keepends=True):
        s = ln.lstrip()
        key = ln.split("=", 1)[0].strip() if "=" in ln and not s.startswith("#") else None
        if key in keys:
            gone.append(key)
            continue
        out.append(ln)
    return "".join(out), sorted(gone)


def choose_dda_template(enzyme, ids_only, nonspecific=None):
    """(DDA_TEMPLATES key, why) for a FragPipe DDA search. An enzyme rule (or --enzyme not given:
    the vendor's trypsin) -> LFQ-MBR when amounts are wanted, Basic-Search for --ids-only (also
    the no-DE timsTOF DDA-PASEF route, Ruling 86). No rule (--enzyme none, pepsin) ->
    Nonspecific-peptidome, or -HLA on request; label-free quantification is switched on in it only
    when amounts are wanted. Raises ValueError for --nonspecific-workflow beside an enzyme rule."""
    none = bool(enzyme) and enzyme["specificity"] == "none"
    if not none:
        if nonspecific:
            raise ValueError(f"--nonspecific-workflow {nonspecific} picks the workflow of a search "
                             f"with no enzyme rule (--enzyme none, or pepsin); this one has "
                             f"{'the rule ' + enzyme['rule'] if enzyme else 'the default trypsin'}")
        if ids_only:
            return "ids", (f"identifications only ({IDS_ONLY_FLAG}): FragPipe's Basic-Search "
                           f"workflow, with no label-free quantification")
        return "lfq", ("amounts wanted (no --ids-only): FragPipe's LFQ-MBR workflow, IonQuant "
                       "label-free quantification with match-between-runs")
    key = nonspecific or NONSPECIFIC_DEFAULT
    which = (f"--nonspecific-workflow {key}" if nonspecific else
             "the default non-specific workflow; --nonspecific-workflow hla for HLA peptides")
    scope = (f"identifications only ({IDS_ONLY_FLAG}): the workflow as shipped, with no "
             f"quantification" if ids_only else
             "amounts wanted: label-free quantification switched on "
             f"({fpp.LFQ_KEY}=true, the one change H7 validated)")
    return key, f"{enzyme['why']} -> FragPipe's {DDA_TEMPLATES[key][1]} workflow ({which}); {scope}"


def edit_var_mods(table, words, fixed_c):
    """MSFragger's msfragger.table.var-mods with the skill's variable modifications: the Met-ox
    row on or off as `words` (VAR_MODS) say, the N/Q deamidation row likewise (DEAM_ROW into the
    first free row when the template has none), protein N-term acetyl (PROTEIN_NTERM_ACETYL)
    on; with a fixed modification on C, a variable one on C alone is switched off (a blocked
    cysteine carries neither cysteinylation nor anything else). The template's other rows are
    kept as shipped. -> (table, [what changed], [the template's own rows still on]). Raises
    ValueError when a row the skill relies on is missing (the template changed)."""
    rows = [r.strip() for r in table.split(";")]
    parsed = [VAR_ROW.match(r) for r in rows]
    changed, kept_on = [], []
    seen = {"ox": False, "deam": False, "acetyl": False}

    def put(i, on, why):
        m = parsed[i]
        new = f"{m.group(1)},{m.group(2)},{'true' if on else 'false'},{m.group(4)}"
        if new != rows[i]:
            changed.append(f"{rows[i]} -> {new} ({why})")
            rows[i] = new
            parsed[i] = VAR_ROW.match(new)

    for i, m in enumerate(parsed):
        if not m:
            continue
        mass, site = m.group(1), m.group(2)
        if site == "M" and _same(mass, VAR_MODS["ox"]["mass"]):
            seen["ox"] = True
            put(i, "ox" in words, "Ox(M)" + (" on" if "ox" in words else " off"))
        elif site == "NQ" and _same(mass, VAR_MODS["deam"]["mass"]):
            seen["deam"] = True
            put(i, "deam" in words, "Deamidation (NQ)" + (" on" if "deam" in words else " off"))
        elif site == "[^" and _same(mass, PROTEIN_NTERM_ACETYL):
            seen["acetyl"] = True
            put(i, True, "protein N-term acetyl, always (estimate_params.PROTEIN_NTERM_ACETYL)")
        elif site == "C" and fixed_c and m.group(3) == "true":
            put(i, False, "a variable modification on C, off: cysteine carries a fixed one")
        elif m.group(3) == "true":
            kept_on.append(rows[i])
    if "ox" in words and not seen["ox"]:
        raise ValueError("the template has no Met-ox (15.9949,M) row to switch on")
    if not seen["acetyl"]:
        raise ValueError(f"the template has no protein N-term acetyl ({PROTEIN_NTERM_ACETYL},[^) row")
    if "deam" in words and not seen["deam"]:
        free = next((i for i, r in enumerate(rows) if re.match(r"^0\.0,site_\d+,false,\d+$", r)),
                    None)
        if free is None:
            raise ValueError("the template has no free variable-modification row for deamidation")
        changed.append(f"{rows[free]} -> {DEAM_ROW} (Deamidation (NQ) on: FragPipe 24.0's own "
                       f"N/Q row, citrullination.workflow)")
        rows[free] = DEAM_ROW
    return "; ".join(rows), changed, kept_on


def fixed_c(props, cys):
    """The fixed-C row of a template for `cys` (cys_alkylation_record) -> (new fix-mods value or
    None when unchanged, the record of what is searched, the fixed mass). A --cys-alkylation that
    was GIVEN sets C's row on at its agent's mass (0.0 for none, as the DIA presets write it). Not
    given, the template's own row stands: carbamidomethyl in LFQ-MBR, Basic-Search and
    Nonspecific-peptidome (the skill's default, recorded as one), and NOTHING fixed in
    Nonspecific-HLA, whose row is off -- recorded as that workflow's default, not as cam."""
    key = "msfragger.table.fix-mods"
    line = props.get(key)
    m = FRAGPIPE_C_ROW.search(line or "")
    if not m:
        raise ValueError(f"no cysteine row in {key}")
    have = float(m.group(1)) if m.group(2) == "true" else 0.0
    if not cys["given"]:
        if _same(have, CYS_REAGENTS["cam"]["mass"]):
            return None, cys, have
        rec = cys_alkylation_record("none")
        rec.update(given=False, source=(
            "DEFAULT of the vendor workflow, not user-confirmed: its fixed C row is off (HLA "
            "peptides are usually not reduced and alkylated), so cysteine is searched unmodified"
            " -- the workflow searches cysteinylation as a variable modification. Pass "
            "--cys-alkylation cam (or the agent used) if the samples were alkylated"))
        return None, rec, have
    want = 0.0 if cys["mass"] is None else cys["mass"]
    if _same(have, want):
        return None, cys, have
    new = FRAGPIPE_C_ROW.sub(lambda _m: f"{want},C (cysteine),true,", line, count=1)
    return new, cys, want


def build_fragpipe_dda(args, label, prov, dda):
    """The FragPipe DDA workflow for this run (see the module docstring); fills `prov`.
    Raises ValueError with what cannot be honoured (main refuses it, nothing written)."""
    enzyme, ids_only = dda["enzyme"], dda["ids_only"]
    key, choice_why = choose_dda_template(enzyme, ids_only, args.nonspecific_workflow)
    rel, vendor_name, desc = DDA_TEMPLATES[key]
    path = os.path.join(PRESETS, rel)
    with open(path, "rb") as fh:
        raw = fh.read()
    text = raw.decode("utf-8")
    props, _ = fpp.read_workflow(path)
    changes, why, warnings = {}, {}, []

    # tolerances: the vendor's, unless the data need FragPipe's ion-trap or unit-resolution setting
    analyzer = args.ms2_analyzer
    if unit_resolution(args.instrument):
        mode, src, tol_keys = "unit_resolution", UNIT_RESOLUTION_SOURCE, UNIT_RESOLUTION_KEYS
        short = TOLERANCE_METHODS[mode]
        if args.ms1_ppm is not None or args.ms2_ppm is not None:
            raise ValueError(f"--ms1-ppm/--ms2-ppm on {args.instrument}: a unit-resolution "
                             "instrument is searched with Da windows (FragPipe's Stellar-DDA); "
                             "give a tolerance with --overrides (msfragger.* keys) instead")
        # Ruling 105.2: identifications only -- THE engine-independent rule (estimate_params),
        # which dda_args() has already applied; again here for a direct caller
        no_amounts = unit_resolution_amounts_refused(args.instrument, "DDA", ids_only)
        if no_amounts:
            raise ValueError(no_amounts)
    elif ms2_in_ion_trap(analyzer):
        mode, tol_keys = "ion_trap", ION_TRAP_KEYS
        short = TOLERANCE_METHODS[mode]
        src = ION_TRAP_SOURCE + (
            "; mixed analyzers (MS2 partly in the ion trap, partly in the Orbitrap): one fragment "
            "window per search, so every MS2 is searched at the ion-trap 0.6 Da, the Orbitrap "
            "ones too (wider than they need, nothing dropped). MSFragger's analyzer filter "
            "(msfragger.analyzer_types) stays 'all': filtering to one analyzer would discard the "
            "other's spectra, and it works on mzML/.raw only" if analyzer == "mixed" else
            "; MS2 read in the ion trap (ITMS)")
        if args.ms2_ppm is not None:
            raise ValueError(f"--ms2-ppm {args.ms2_ppm:g} with --ms2-analyzer {analyzer}: MS2 read "
                             f"in the ion trap needs a Da fragment window (FragPipe's 0.6 Da), "
                             f"not ppm; give another with --overrides "
                             f"'{{\"msfragger.fragment_mass_tolerance\": ...}}'")
    else:
        mode, tol_keys = "vendor", {}
        src = (f"the vendor {vendor_name} workflow's values kept: precursor "
               f"{props.get('msfragger.precursor_mass_lower')}/+{props.get('msfragger.precursor_mass_upper')} "
               f"{fpp.UNITS.get(props.get('msfragger.precursor_mass_units', ''), '?')}, fragment "
               f"{props.get('msfragger.fragment_mass_tolerance')} "
               f"{fpp.UNITS.get(props.get('msfragger.fragment_mass_units', ''), '?')} -- the "
               f"workflow's own values, used on every instrument (timsTOF too); "
               + ("starting values for MSFragger's mass calibration (calibrate_mass="
                  f"{props.get('msfragger.calibrate_mass')})"
                  if (props.get('msfragger.calibrate_mass') or '0') not in ('0', '')
                  else "searched as given (calibrate_mass=0: no mass calibration)")
               + "; DIA-NN's instrument table is not applied to MSFragger")
    for k, v in tol_keys.items():
        if props.get(k) != v:
            changes[k] = v
            why[k] = f"{v} -- {short} (the record's tolerances.source says why)"
    if args.ms1_ppm is not None:
        v = f"{args.ms1_ppm:g}"
        for k, val in (("msfragger.precursor_mass_lower", f"-{v}"),
                       ("msfragger.precursor_mass_upper", v),
                       ("msfragger.precursor_mass_units", "1"),
                       ("msfragger.precursor_true_tolerance", v),
                       ("msfragger.precursor_true_units", "1")):
            changes[k], why[k] = val, f"{val} -- precursor +/-{v} ppm (caller override, --ms1-ppm)"
        mode, src = "caller", f"--ms1-ppm {v}" + (f"; {src}" if tol_keys else "")
    if args.ms2_ppm is not None:
        v = f"{args.ms2_ppm:g}"
        for k, val in (("msfragger.fragment_mass_tolerance", v),
                       ("msfragger.fragment_mass_units", "1")):
            changes[k], why[k] = val, f"{val} -- fragment +/-{v} ppm (caller override, --ms2-ppm)"
        mode, src = "caller", (f"{src}; " if mode == "caller" else "") + f"--ms2-ppm {v}"

    # the digest of an enzymatic template: the skill's DDA digest, shared with Sage (Ruling 57.1:
    # estimate_params.DDA_DIGEST -- trypsin, not before P; its missed cleavages and lengths), or
    # the --enzyme rule with the same missed cleavages and lengths; written as FragPipe's own
    # custom-rule form. A non-specific template keeps its own digest and lengths.
    if key not in NONSPECIFIC_WORKFLOWS:
        if enzyme and enzyme["specificity"] in ("full", "semi"):
            rule, rule_why = enzyme, enzyme_source(enzyme)
        else:
            rule = {"cut": DDA_DIGEST["cleave_at"], "nocut": DDA_DIGEST.get("restrict") or "",
                    "sense": "C" if DDA_DIGEST.get("c_terminal", True) else "N"}
            rule_why = DDA_DIGEST_DEFAULT
        for k, v in cleavage_rules.fragger_enzyme_keys(rule).items():
            if props.get(k) != v:
                changes[k], why[k] = v, rule_why
        for k, v in (("msfragger.allowed_missed_cleavage_1", DDA_DIGEST["missed_cleavages"]),
                     ("msfragger.digest_min_length", DDA_DIGEST["min_len"]),
                     ("msfragger.digest_max_length", DDA_DIGEST["max_len"])):
            if props.get(k) != str(v):
                changes[k] = str(v)
                why[k] = (f"{v} -- the skill's DDA digest, shared with Sage (estimate_params."
                          f"DDA_DIGEST; Ruling 57.1); the vendor workflow has {props.get(k)}")
    # the scope: label-free quantification on in a non-specific workflow when amounts are wanted
    if key in NONSPECIFIC_WORKFLOWS and not ids_only:
        changes[fpp.LFQ_KEY] = "true"
        why[fpp.LFQ_KEY] = ("true -- amounts wanted (no --ids-only): label-free quantification "
                            "with IonQuant, the one change to Nonspecific-peptidome H7 validated "
                            "(combined_protein.tsv written)")
    # variable modifications (estimate_params.VAR_MODS; the DDA default Met-ox) and the cysteine
    vm_words, vm_src = var_mods_choice("fragpipe", "DDA", dda["var_mods"])
    c_line, cys, c_mass = fixed_c(props, dda["cys"])
    if c_line is not None:
        changes["msfragger.table.fix-mods"] = c_line
        why["msfragger.table.fix-mods"] = (f"C (cysteine) fixed at {c_mass} -- {cys['searched']}: "
                                           f"{cys['source']} (the vendor workflow's row differs)")
    table, vm_changed, vm_kept = edit_var_mods(props["msfragger.table.var-mods"], vm_words,
                                               bool(c_mass))
    if vm_changed:
        changes["msfragger.table.var-mods"] = table
        why["msfragger.table.var-mods"] = ("; ".join(vm_changed)
                                           + f" -- variable modifications: {vm_src}")
    # at most DDA_MAX_VAR_MODS per peptide in an enzymatic preset, as Sage; a non-specific
    # workflow keeps FragPipe's own number (Ruling 105: what H7 validated, what FragPipe ships
    # for peptidomics), recorded as FragPipe's default
    vendor_mv = props.get("msfragger.max_variable_mods_per_peptide")
    if key in NONSPECIFIC_WORKFLOWS:
        max_vm = {"value": int(vendor_mv) if str(vendor_mv).isdigit() else vendor_mv,
                  "source": f"DEFAULT, not user-confirmed: FragPipe's own {vendor_name} workflow "
                            f"(the skill's limit of {DDA_MAX_VAR_MODS} applies to enzymatic "
                            f"searches only; Ruling 105)"}
    else:
        mv = str(DDA_MAX_VAR_MODS)
        max_vm = {"value": DDA_MAX_VAR_MODS,
                  "source": "the skill's DDA limit, as Sage's max_variable_mods "
                            "(estimate_params.DDA_MAX_VAR_MODS)"}
        if vendor_mv != mv:
            changes["msfragger.max_variable_mods_per_peptide"] = mv
            why["msfragger.max_variable_mods_per_peptide"] = (
                f"{mv} -- {max_vm['source']}; the vendor workflow has {vendor_mv}")

    # --directlfq: IonQuant's MaxLFQ off, so the report's quantity is directLFQ (fragpipe_quant)
    if dda.get("directlfq"):
        changes[DIRECTLFQ_KEY] = "0"
        why[DIRECTLFQ_KEY] = (f"0 -- {DIRECTLFQ_FLAG}: directLFQ requested for the protein "
                              f"quantities instead of IonQuant's MaxLFQ (the vendor workflow has "
                              f"{props.get(DIRECTLFQ_KEY)})")
    # --overrides: workflow keys of this template, as the user chose them; applied last
    over, om = dda["overrides"], dda["override_meta"]
    if not isinstance(over, dict):
        raise ValueError("--overrides for a FragPipe preset is a JSON object of workflow keys, "
                         "e.g. {\"msfragger.allowed_missed_cleavage_1\": 1}")
    for k, v in over.items():
        if k in FRAGPIPE_OVERRIDE_REFUSED:
            raise ValueError(f"--overrides sets {k}: {FRAGPIPE_OVERRIDE_REFUSED[k]}")
        if k not in props:
            raise ValueError(f"--overrides sets {k}, which is not a key of {rel} (FragPipe 24.0's "
                             f"{vendor_name} workflow): check its spelling there")
        if isinstance(v, (dict, list)) or v is None:
            raise ValueError(f"--overrides sets {k} to {json.dumps(v)}: a workflow value is one "
                             "string, number or true/false")
        if dda.get("directlfq") and k == DIRECTLFQ_KEY and str(v).lower() not in ("0", "false"):
            raise ValueError(f"{DIRECTLFQ_FLAG} sets {k}=0, and --overrides sets it to {v!r}: "
                             "give one or the other")
        if ids_only and k == fpp.LFQ_KEY:
            raise ValueError(f"{IDS_ONLY_FLAG} searches without label-free quantification, and "
                             f"--overrides sets {k}: give one or the other")
        val = ("true" if v else "false") if isinstance(v, bool) else v
        was = why.get(k)
        changes[k] = val
        why[k] = override_source(k, val, om.get("reason"), om.get("role")) + (
            f" (instead of: {was})" if was else "")
    # a rule the workflow states itself is named `nonspecific` and has NO dropdown key beside it
    deleted = []
    rule_keys = [k for k in cleavage_rules.FRAGPIPE_RULE_KEYS if k in changes]
    if rule_keys and key not in NONSPECIFIC_WORKFLOWS:
        nk = "msfragger.search_enzyme_name_1"
        # the rule as written: a row of FragPipe's enzyme list is written under its own name AND
        # its dropdown (FragPipe reloads the same rule from it), so the workflow opened in
        # FragPipe shows the enzyme searched, not `nonspecific` (X8 F6); any other rule is the
        # custom form, no dropdown key
        final = {k: changes.get(k, props.get(k)) for k in (
            "msfragger.search_enzyme_cut_1", "msfragger.search_enzyme_nocut_1",
            "msfragger.search_enzyme_sense_1", "msfragger.num_enzyme_termini")}
        entry = (cleavage_rules.fragpipe_list_entry(
            final["msfragger.search_enzyme_cut_1"], final["msfragger.search_enzyme_nocut_1"],
            final["msfragger.search_enzyme_sense_1"]) if final["msfragger.num_enzyme_termini"] == "2"
            and nk not in over else None)
        if entry:
            changes[nk] = entry
            changes[cleavage_rules.FRAGPIPE_DROPDOWN_KEY] = entry
            why[nk] = why[cleavage_rules.FRAGPIPE_DROPDOWN_KEY] = (
                f"{entry} -- the rule written is that row of FragPipe 24.0's enzyme list; its name "
                f"and dropdown say so, and FragPipe reloads the same cut, no-cut and sense from it "
                f"(TabMsfragger.java 24.0)")
        else:
            if nk not in changes:
                changes[nk] = cleavage_rules.FRAGPIPE_CUSTOM_NAME
                why[nk] = (f"{cleavage_rules.FRAGPIPE_CUSTOM_NAME} -- FragPipe 24.0's own form for a "
                           f"rule the workflow states (TabMsfragger), for "
                           f"{', '.join(rule_keys)}")
            deleted.append(cleavage_rules.FRAGPIPE_DROPDOWN_KEY)
    out, patched, added = patch_properties(text, changes)
    out, gone = delete_properties(out, deleted)

    tag = (f"{cleavage_rules.FRAGPIPE_DROPDOWN_KEY} deleted: loading it, FragPipe reloads the "
           f"enzyme from its list over the rule written here -- `nonspecific` would make the "
           f"search non-specific, `custom` ran as argc (C0, 2026-11-04; the rule reaches "
           f"fragger.params only without it)")
    prov.update({
        "engine": "fragpipe",
        "template": rel,
        "template_note": desc,
        "template_is_vendor_default": True,
        "template_sha256": hashlib.sha256(raw).hexdigest(),
        "template_choice": {"workflow": vendor_name, "why": choice_why},
        "ids_only": ids_only,
        "quantity_requested": ("directlfq" if dda.get("directlfq") else None),
        "var_mods": {"words": vm_words, "source": vm_src + "; + protein N-term acetyl (always)",
                     "kept_from_template": vm_kept},
        "max_var_mods": max_vm,
        "enzyme": enzyme,
        "cys_alkylation": cys,
        "tolerances": {"mode": mode, "source": src},
        "keys_changed": {k: why[k] for k in changes},
        "keys_deleted": gone,
        "keys_deleted_why": {k: tag for k in gone},
        "keys_patched_in_place": patched,
        "keys_appended": added,
        "database": DATABASE_NOTE + (f" (--fasta {os.path.abspath(args.fasta)} was given here "
                                     f"and is recorded only)" if args.fasta else ""),
        "threads": THREADS_NOTE + (f" (--threads {args.threads} was given here and is recorded "
                                   f"only)" if args.threads else ""),
        "warnings": warnings,
    })
    return out


def dda_record_as_written(prov, out, dda):
    """What the WRITTEN workflow searches, read back with fragpipe_params (the one reader), into
    `prov`: its cleavage (cleavage_rules.fragger_cleavage, tagged with where the rule came from),
    whether it quantifies, and any enzyme dropdown that would overwrite its rule."""
    props, _ = fpp.read_workflow(out)
    kv = {k[len("msfragger."):]: v for k, v in props.items() if k.startswith("msfragger.")}
    cl = cleavage_rules.fragger_cleavage(fpp.digest(kv), os.path.basename(out))
    enzyme = dda["enzyme"]
    stated = [k for k in cleavage_rules.FRAGPIPE_RULE_KEYS if k in dda["overrides"]]
    if stated:
        src = "user override (--overrides " + ", ".join(stated) + ")"
    elif enzyme:
        src = enzyme_source(enzyme)
    else:
        src = DDA_DIGEST_DEFAULT
    cl["source"] = src
    q = fpp.workflow_quantification(props)
    prov.update({"cleavage": cl, "quantified": q["quantified"], "quantification_why": q["why"],
                 # whether the written filter sets a protein-level FDR, and why not when it does
                 # not (the Nonspecific workflows' --prot 1 is kept: Ruling 120 condition 1)
                 "protein_fdr": fpp.protein_fdr_record(props.get(fpp.PHILOSOPHER_FILTER_KEY))})
    # what the written workflow fixes on C (estimate_params.cys_searched, the one reader): an
    # --overrides of the fixed-mod table replaces what --cys-alkylation said (review R1's rule)
    if "msfragger.table.fix-mods" in dda["overrides"]:
        prov["cys_alkylation"] = cys_record_as_written(
            prov["cys_alkylation"], cys_searched(out), dda["override_meta"],
            how="an msfragger.table.fix-mods override sets C's row")
    # the skill never writes a rule under an enzyme dropdown that would change it (Ruling 104):
    # its own preset is refused on it (main sets the written file aside); a user's workflow never is
    conflict = cleavage_rules.fragpipe_dropdown_conflict(props)
    if conflict:
        raise ValueError(f"the preset written would not search the rule it states -- {conflict}")


def radiant_tolerances(cls, m1, m2, label):
    """One line for the Radiant provenance: which extraction widths were narrowed, and why any
    level kept the vendor 20 ppm."""
    if m1 is not None and m2 is not None:
        return "narrowed from the instrument"
    why = ("instrument not recognised" if cls == "unknown" else
           f"no documented DIA-NN tolerance for {label}")
    if m1 is None and m2 is None:
        return f"vendor 20/20 ppm kept -- {why}"
    kept = "MS2" if m2 is None else "MS1"
    why_kept = ("MS2 is read in the ion trap (all of it, or part: see the label), which has no "
                "Orbitrap-table value"
                if cls == "orbitrap_iontrap" and kept == "MS2" else
                f"no documented DIA-NN tolerance for this {kept} resolution")
    return (f"{'MS1' if kept == 'MS2' else 'MS2'} narrowed from the instrument; {kept} vendor "
            f"20 ppm kept -- {why_kept} ({label})")


# The Radiant route is EXPERIMENTAL, and its records say so (route_status, read by the Methods
# and the report: rule 1, the route describes itself). Audit 2026-11-17, lane 1 #3: its
# extraction widths are DIA-NN's documented mass-accuracy guidance applied to Radiant's
# ms1/ms2ExtractionWidthPPM -- our inference, not Seer's, and DIA-NN's README calls its own values
# "guidance", not extraction windows. Whether Radiant's parameter is a full or a half width is
# unknown, and no benchmark on Core data backs the narrowed values.
RADIANT_ROUTE_STATUS = {
    "status": "experimental",
    "reasons": [
        "extraction widths narrowed from DIA-NN's documented mass-accuracy table, which DIA-NN "
        "calls guidance rather than an extraction window; applying it to Radiant's "
        "ms1/ms2ExtractionWidthPPM is the skill's inference, not the vendor's, and has not been "
        "benchmarked on Core data (whether Radiant's value is a full or a half width is "
        "unknown)",
        "the spectral library is DIA-NN's prediction converted to a TSV, not a vendor library",
    ],
    "say": "Radiant results are experimental: use the DIA-NN route for results of record."}


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
            f"{m1} ppm for {label} [DERIVED from DIA-NN's instrument table, UNVALIDATED for "
            f"Radiant (route_status: experimental); Seer ships a flat 20.0 for every Orbitrap]"
            if args.ms1_ppm is None else f"{m1} ppm (caller override)")
    if m2 is not None:
        changes[("MS2Params", "ms2ExtractionWidthPPM")] = float(m2)
        why["[MS2Params] ms2ExtractionWidthPPM"] = (
            f"{m2} ppm for {label} [DERIVED as above, UNVALIDATED for Radiant]"
            if args.ms2_ppm is None else f"{m2} ppm (caller override)")
    if args.threads:
        changes[("General", "threadCount")] = int(args.threads)
        why["[General] threadCount"] = "this run's CPU allocation (0 = all cores)"

    out, applied = patch_ini(text, changes)
    # Radiant takes its modifications from the DIA-NN library it searches (its config has none):
    # run_search.py / radiant_parallel.py build that library with these flags, read from here.
    why["cysteine (DIA-NN library prediction)"] = (
        f"{' '.join(cys_diann_flags(prov['cys_alkylation']['value'])) or 'no fixed modification'}"
        f" -- {prov['cys_alkylation']['searched']}: {prov['cys_alkylation']['source']}")
    prov.update({
        "engine": "radiant",
        "template": "radiant/default.radiantConfig",
        "template_note": "radiant_fulcrum_workflow 0.3.3 shipped default, copied verbatim",
        "template_is_vendor_default": True,
        "keys_changed": why,
        "keys_patched_in_place": applied,
        # An Orbitrap level with no documented DIA-NN tier (e.g. a 15k MS2) keeps the vendor
        # width. Its DIA-NN value is measured with DIA-NN before a DIA-NN search, which Radiant
        # does not run, and the table is no longer extrapolated -- so say which level stayed at
        # the vendor default and why, rather than "narrowed" or "not recognised".
        "tolerances": radiant_tolerances(cls, m1, m2, label),
        "route_status": RADIANT_ROUTE_STATUS,
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
    add_resolution_args(ap)       # --ms1-res/--ms1-resolution, and --resolution-source
    ap.add_argument("--ms1-ppm", type=float, default=None, help="override; else vendor/derived")
    ap.add_argument("--ms2-ppm", type=float, default=None, help="override; else vendor/derived")
    ap.add_argument("--fasta", default=None,
                    help="DIA: written as database.db-path. FragPipe DDA: recorded only -- "
                         "run_search.py sets the database at search time")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--template", choices=["timstof", "orbitrap"],
                    help="force a FragPipe DIA template instead of choosing by instrument")
    # --var-mods --ids-only --cys-alkylation --enzyme --overrides --override-by --override-reason:
    # estimate_params.py's, one definition. --cys-alkylation: FragPipe's fixed C row (DIA and DDA)
    # and Radiant's DIA-NN library; the others are the FragPipe DDA preset's.
    add_dda_args(ap)
    ap.add_argument("--nonspecific-workflow", choices=NONSPECIFIC_WORKFLOWS, default=None,
                    help="FragPipe DDA with no enzyme rule (--enzyme none, or pepsin): "
                         "peptidome = FragPipe's Nonspecific-peptidome workflow (the default, "
                         "recorded as one), hla = its Nonspecific-HLA workflow (HLA peptides)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    def refuse(msg, code=1):
        refuse_out(args.out, msg, code=code, sidecars=(".provenance.json",), who="make_presets")

    fp_dda = args.engine == "fragpipe" and args.acquisition == "DDA"
    if args.engine == "radiant" and args.acquisition == "DDA":
        sys.exit("make_presets: Radiant is a DIA engine. DDA goes through FragPipe "
                 "(make_presets.py --engine fragpipe --acquisition DDA) or Sage "
                 "(estimate_params.py --engine sage).")
    dda_only = [flag for flag, given in (
        ("--var-mods", args.var_mods is not None), (IDS_ONLY_FLAG, args.ids_only),
        (DIRECTLFQ_FLAG, args.directlfq),
        (ENZYME_FLAG, args.enzyme is not None), (SEMI_FLAG, args.semi_specific),
        ("--overrides", bool(args.overrides)),
        ("--override-by", args.override_by is not None),
        ("--override-reason", args.override_reason is not None),
        ("--nonspecific-workflow", args.nonspecific_workflow is not None)) if given]
    if dda_only and not fp_dda:
        refuse(f"{', '.join(dda_only)}: for the FragPipe DDA preset (--engine fragpipe "
               f"--acquisition DDA). The {args.engine} {args.acquisition} preset is the vendor's, "
               f"set by --cys-alkylation and the tolerance flags only. Nothing was written.", 2)
    if fp_dda and args.template:
        refuse("--template picks a FragPipe DIA preset; a DDA preset is chosen by --ids-only and "
               "--enzyme. Nothing was written.", 2)

    # resolve_defaults.py refuses this pair before calling here; a direct call must not get
    # past it either
    bad = resolution_args_error(args.ms2_resolution, args.ms2_analyzer)
    if bad:
        sys.exit(f"make_presets: {bad}")
    # the source goes into the label, and so into every "keys_changed" reason below
    cls, ms1, ms2, label, src = classify_instrument(
        args.instrument, args.ms1_resolution, args.ms2_resolution,
        resolution_source(args.resolution_source, args.ms1_resolution, args.ms2_resolution,
                          args.ms2_analyzer), args.ms2_analyzer)
    prov = {"instrument": args.instrument or None, "instrument_class": cls,
            "instrument_label": label, "acquisition": args.acquisition,
            "ppm_source": src, "output": os.path.abspath(args.out),
            # read back by run_search.py (Radiant's library) and make_methods.py
            "cys_alkylation": cys_alkylation_record(args.cys_alkylation)}

    dda = None
    if fp_dda:
        # the DDA flags, read and refused as estimate_params.py reads them (dda_args)
        try:
            dda = dda_args(args, "fragpipe", who="make_presets")
        except DdaArgsError as e:
            if e.refuse:
                refuse(e.msg.replace("No cfg was written", "Nothing was written"), e.code)
            sys.exit(e.msg)
        try:
            text = build_fragpipe_dda(args, label, prov, dda)
        except ValueError as e:
            refuse(f"{str(e).rstrip('.')}. Nothing was written.")
    else:
        text = (build_fragpipe if args.engine == "fragpipe" else build_radiant)(
            args, cls, ms1, ms2, label, src, prov)

    d = os.path.dirname(os.path.abspath(args.out))
    if d:
        os.makedirs(d, exist_ok=True)
    with open(args.out, "w") as fh:
        fh.write(text)
    if dda is not None:
        try:
            dda_record_as_written(prov, args.out, dda)
        except ValueError as e:
            refuse(f"{str(e).rstrip('.')}. Nothing usable was written.")
        over = dda["overrides"]
        # who set the overrides, by name: staff-only, beside the preset (never in the record)
        staff_rec = (record_override_setter(args.out, over, dda["by"], dda["by_how"],
                                            args.override_reason, dda["role"]) if over else None)
        prov["overrides"] = override_record(over, args.override_reason, staff_rec,
                                            role=dda["role"])
    with open(args.out + ".provenance.json", "w") as fh:
        json.dump(prov, fh, indent=2)
    print(json.dumps(prov, indent=2))


if __name__ == "__main__":
    main()
