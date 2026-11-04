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
                                 before the search; the level that has a tier keeps it.
                                 DDA (--dda) cannot be measured that way, so there the
                                 level is pinned at the facility SOP, tagged DEFAULT
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
      --instrument "timsTOF Pro" --out sage_config.json [--var-mods ox] [--ids-only]
      (--ids-only: identifications only, no differential expression -- the cfg's quant.lfq is
       false, and run_search.py then hands a timsTOF .d to Sage directly)
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
specific fields the user chose, with --override-by "<who chose them>" and
--override-reason "<why, in their words>": each is tagged as a USER OVERRIDE with its
value, who set it by role, and why (override_source) -- never as a validated SOP, which the
skill cannot know it is. The role is "Core staff" only when the setter is a login on the Core
staff list (staff.role_for), else "the user". Who, by name, goes only to the staff-only
<out>.staff.json (staff.py): the sidecar travels to the client.
"""
import sys, os, json, math, argparse

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import staff  # noqa: E402  the role client documents name, and the staff-only record
import silac  # noqa: E402  SILAC labels (--labels): the one definition of the plexDIA declaration

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

# The in-silico digest of the DIA-NN search: --cut (DIA-NN's syntax), --missed-cleavages,
# --min/--max-pep-len and --met-excision, written into every cfg by build_diann(). The ONE
# definition: fetch_fasta.py digests contaminant and target entries with it to decide whether
# a contaminant can be told apart from a target protein by the peptides this search looks for.
DIANN_DIGEST = {"cut": "K*,R*", "missed_cleavages": 1, "min_pep_len": 7, "max_pep_len": 30,
                "met_excision": True}

# DIA-NN's DDA switch. The ONE definition of "this cfg searches DDA": build_diann() writes it for
# --acquisition DDA, and diann_parallel.py reads it back (the chain's routing, the refusal of a
# cfg that does not match the bundle's acquisition, probe_window.py's refusal to measure).
# It used to live nowhere: run_search.py appended it to the single-shot command only, so every
# DDA cohort of more than 5 files went to the 5-step chain WITHOUT it and DIA-NN searched DDA
# spectra as DIA, with no error (msalemi, SET28 / PROT_0000, 28 Exploris .raw, 2026-09-25).
DIANN_DDA_FLAG = "--dda"
SRC_DDA = ("DDA acquisition: DIA-NN README, '--dda process data as DDA -- must be used with DDA "
           "data, must not be used with DIA data' (beta DDA support since DIA-NN 2.3)")


def is_dda(acquisition):
    return (acquisition or "").strip().upper() == "DDA"


# What the search assumes the sample prep did to cysteine. THE one definition: build_diann() and
# build_sage() write it, the sidecar records it (`cys_alkylation`), cys_searched() reads it back
# from a params file (sample_quality.py's cysteine check), make_presets.py applies it to FragPipe
# and to Radiant's DIA-NN library, make_methods.py names the modification, and SKILL.md step 3
# asks for it. Each alkylating agent leaves its own FIXED modification on C; masses and names are
# UniMod's (unimod.org/obo/unimod.obo, read 2026-10-15):
#   cam            iodoacetamide or chloroacetamide -> Carbamidomethyl (the standard prep)
#   none           not reduced and alkylated: cysteine is searched unmodified
#   nem            N-ethylmaleimide                 -> Nethylmaleimide
#   carboxymethyl  iodoacetic acid                  -> Carboxymethyl
#   propionamide   acrylamide                       -> Propionamide (the Core has seen it)
#   mmts           S-methyl methanethiosulfonate    -> Methylthio
#   pyridylethyl   4-vinylpyridine                  -> Pyridylethyl
# On-bead, IP and pull-down digests are often not reduced and alkylated, and a fixed
# carbamidomethyl then leaves no cysteine peptide searchable. Staff reports (2026-10-09): 0.03-0.5%
# of the identified peptides held a C on three such searches, against 7.8% on an S-Trap prep with
# DTT/IAA; the only way out was to delete --unimod4 from params.cfg by hand.
CYS_REAGENTS = {
    "cam": {"reagent": "iodoacetamide or chloroacetamide", "mod": "Carbamidomethyl",
            "unimod": 4, "mass": 57.021464},
    "none": None,
    "nem": {"reagent": "N-ethylmaleimide (NEM)", "mod": "Nethylmaleimide",
            "unimod": 108, "mass": 125.047679},
    "carboxymethyl": {"reagent": "iodoacetic acid", "mod": "Carboxymethyl",
                      "unimod": 6, "mass": 58.005479},
    "propionamide": {"reagent": "acrylamide", "mod": "Propionamide",
                     "unimod": 24, "mass": 71.037114},
    "mmts": {"reagent": "S-methyl methanethiosulfonate (MMTS)", "mod": "Methylthio",
             "unimod": 39, "mass": 45.987721},
    "pyridylethyl": {"reagent": "4-vinylpyridine", "mod": "Pyridylethyl",
                     "unimod": 31, "mass": 105.057849},
}
CYS_ALKYLATION_CHOICES = tuple(CYS_REAGENTS)
CYS_ALKYLATION_DEFAULT = "cam"
# Sage's own carbamidomethyl value, as the config always carried it
SAGE_CAM = {"C": 57.0215}


def cys_diann_flags(choice):
    """The DIA-NN flags for cysteine: --unimod4 (DIA-NN's own carbamidomethyl-C) for cam, nothing
    for none, else `--fixed-mod UniMod:<n>,<mass>,C` (DIA-NN README: "--fixed-mod
    [name],[mass],[sites]", applied during the FASTA digest)."""
    r = CYS_REAGENTS[choice]
    if r is None:
        return []
    if choice == "cam":
        return ["--unimod4"]
    return ["--fixed-mod", f"UniMod:{r['unimod']},{r['mass']},C"]


def cys_sage_static(choice):
    """Sage's database.static_mods for cysteine."""
    r = CYS_REAGENTS[choice]
    return {} if r is None else dict(SAGE_CAM) if choice == "cam" else {"C": r["mass"]}


def cys_alkylation_record(choice=None):
    """The sidecar's `cys_alkylation`: what was searched for cysteine, and on whose word. `choice`
    is --cys-alkylation as given; None is the default, tagged as one (DE-LIMP rule 2)."""
    value = choice or CYS_ALKYLATION_DEFAULT
    r = CYS_REAGENTS[value]
    if choice is None:
        source = ("DEFAULT, not user-confirmed: the samples are assumed reduced and alkylated "
                  "(iodoacetamide or chloroacetamide), the standard prep. Pass --cys-alkylation "
                  "none for samples that were not -- on-bead, IP and pull-down digests often are "
                  "not -- or the agent the lab used")
    elif r is None:
        source = ("given (--cys-alkylation none): the samples were not reduced and alkylated, so "
                  "cysteine is searched unmodified")
    elif value == "cam":
        source = "given (--cys-alkylation cam): the samples were reduced and alkylated"
    else:
        source = (f"given (--cys-alkylation {value}): the samples were alkylated with "
                  f"{r['reagent']}, searched as fixed {r['mod']} (UniMod:{r['unimod']}, "
                  f"+{r['mass']}) on C")
    return {"value": value, "given": choice is not None,
            "reagent": r and r["reagent"], "mod": r and r["mod"],
            "unimod": r and r["unimod"], "mass": r and r["mass"],
            "searched": (f"{r['mod'].lower()} (C), fixed" if r else
                         "cysteine unmodified (no fixed modification)"),
            "source": source}


def cys_record_as_written(cys, found, override_meta=None):
    """`cys` (cys_alkylation_record) when the cfg as written (`found`: cys_searched()) fixes what it
    says on C, else the record of what the cfg does fix, with `overridden` naming the choice an
    --overrides value replaced (review R1). The reagent is then not known -- an override names a
    modification, not what the lab used -- so it is None, and so is `value` when nothing is fixed
    (`none` would claim the samples were not alkylated)."""
    if not found:
        return cys
    same = (bool(found["fixed"]) == bool(cys.get("mod")) and
            (not found["fixed"] or abs(float(found["mass"]) - float(cys["mass"])) < 0.01))
    if same:
        return cys
    om = override_meta or {}
    reason = om.get("reason")
    choice = found.get("choice") if found["fixed"] else None
    r = CYS_REAGENTS.get(choice) if choice else None
    searched = (f"{found['mod'].lower()} (C), fixed" if found["fixed"] else
                "cysteine unmodified (no fixed modification)")
    return {"value": choice, "given": True, "reagent": None,
            "mod": found.get("mod") if found["fixed"] else None,
            "unimod": r["unimod"] if r else None,
            "mass": found.get("mass") if found["fixed"] else None,
            "searched": searched,
            "overridden": {"was": cys.get("value"), "was_searched": cys.get("searched")},
            "source": (f"user override: --overrides replaced what --cys-alkylation "
                       f"{cys.get('value')} fixes on C ({cys.get('searched')}; a --fixed-mod "
                       f"override replaces every --fixed-mod line), so the cfg as written "
                       f"searches {searched}; set by {om.get('role', staff.USER_ROLE)}; reason: "
                       f"{reason.strip() if reason and reason.strip() else OVERRIDE_REASON_MISSING}")}


def _cys_choice(mass):
    """The CYS_REAGENTS key whose mass this is (within 0.01 Da), or None."""
    return next((k for k, r in CYS_REAGENTS.items() if r and abs(r["mass"] - mass) < 0.01), None)


def cys_searched(params):
    """What the params file fixes on C: {"fixed": False} (unmodified), {"fixed": True, "choice":
    a CYS_REAGENTS key or None for an unlisted one, "mod", "mass"}, or None when the file cannot
    be read. A DIA-NN cfg: --unimod4, or a --fixed-mod on C; a Sage config: database.static_mods C."""
    from diann_parallel import cfg_groups, cfg_tokens, CfgError   # here: it imports this module

    def fixed(mass, name=None):
        choice = _cys_choice(mass)
        r = CYS_REAGENTS.get(choice) if choice else None
        return {"fixed": True, "choice": choice, "mass": mass,
                "mod": r["mod"] if r else (name or f"+{mass:g} on C")}
    try:
        if str(params).lower().endswith(".json"):
            with open(params) as fh:
                mods = (json.load(fh).get("database") or {}).get("static_mods") or {}
            for k, v in mods.items():
                if "C" in str(k) and float(v):
                    return fixed(float(v))
            return {"fixed": False}
        groups = cfg_groups(cfg_tokens(params))
    except (OSError, ValueError, TypeError, AttributeError, CfgError):
        return None
    for flag, vals in groups:
        if flag == "--unimod4":
            return fixed(CYS_REAGENTS["cam"]["mass"])
        if flag == "--fixed-mod" and vals:
            parts = [x.strip() for x in ",".join(vals).split(",")]
            try:
                mass = float(parts[1]) if len(parts) > 1 else None
            except ValueError:
                mass = None
            if mass and len(parts) > 2 and "C" in parts[2] and len(parts) < 4:
                return fixed(mass, parts[0])
    return {"fixed": False}


def cam_searched(params):
    """True when the params file fixes carbamidomethyl on C, False when it does not (unmodified, or
    another fixed C modification), None when it cannot be read (cys_searched)."""
    r = cys_searched(params)
    return None if r is None else (r["fixed"] and r["choice"] == "cam")


# The variable modifications --var-mods can ask for: THE one table. word -> its label, its UniMod
# record (id, title, monoisotopic delta mass -- unimod.org: Oxidation 35, Deamidated 7) and
# sites, and from those DIA-NN's `--var-mod UniMod:<id>,<mass>,<sites>` value (README: "--var-mod
# [name],[mass],[sites]") and Sage's database.variable_mods entries ({residue: [mass]});
# make_methods.UNIMOD names them from the same record, as it does CYS_REAGENTS. Any other word is
# refused (parse_var_mods): it used to be dropped without a word, so a "like-for-like" re-search
# of a FragPipe study that asked for N/Q deamidation ran without it (I47, 2026-10-29).
VAR_MODS = {
    "ox": {"label": "Ox(M)", "unimod": 35, "mod": "Oxidation", "mass": 15.994915, "sites": "M"},
    "deam": {"label": "Deamidation (NQ)", "unimod": 7, "mod": "Deamidated", "mass": 0.984016,
             "sites": "NQ"},
}
for _vm in VAR_MODS.values():
    _vm["diann"] = f"UniMod:{_vm['unimod']},{_vm['mass']},{_vm['sites']}"
    _vm["sage"] = {site: [_vm["mass"]] for site in _vm["sites"]}
#: --var-mods' explicit off switch: none of the table's modifications, whatever the default.
VAR_MODS_NONE = "none"
VAR_MODS_WORDS = ", ".join(list(VAR_MODS) + [VAR_MODS_NONE])     # "ox, deam, none"
# What is searched when --var-mods is not given, by (engine, acquisition); anything not listed
# searches none. Met-ox for Sage DDA is the FragPipe/MaxQuant default, and without it those
# peptides are lost (Brett, 2026-10-30, Q2). DIA-NN keeps none: its README advises against
# variable mods for quantification.
VAR_MODS_DEFAULT = {("sage", "DDA"): ["ox"]}
SRC_VAR_MODS_SAGE_DDA = ("Ox(M) on by default for DDA (the FragPipe/MaxQuant default; Core "
                         "decision 2026-10-30)")
SRC_VAR_MODS_DIANN = "DIA-NN README: variable mods do not improve depth for relative quant; left off"


def parse_var_mods(text):
    """--var-mods as given ("ox", "ox,deam", "none") -> the VAR_MODS words, in order, once each;
    `none` -> []. No word at all ("," or " , ") -> None, the same as not given: the default
    applies (it is not `none`, which would silently switch off Sage DDA's Met-ox). Raises
    ValueError naming the known words for any other word, and for `none` beside another word
    (they contradict)."""
    words = []
    for w in (x.strip().lower() for x in str(text or "").split(",")):
        if w and w not in words:
            words.append(w)
    if not words:
        return None
    bad = [w for w in words if w not in VAR_MODS and w != VAR_MODS_NONE]
    if bad:
        raise ValueError(f"--var-mods {', '.join(bad)}: not a modification this skill can "
                         f"search. Known words: {VAR_MODS_WORDS} (ox = Ox(M), deam = "
                         f"Deamidation (NQ), none = no variable modifications). A modification "
                         f"not in the list goes in --overrides.")
    if VAR_MODS_NONE in words and len(words) > 1:
        raise ValueError(f"--var-mods {','.join(words)}: `none` turns variable modifications "
                         f"off and cannot be combined with another word. Known words: "
                         f"{VAR_MODS_WORDS}")
    return [] if words == [VAR_MODS_NONE] else words


def var_mods_choice(engine, acquisition, var_mods=None):
    """(words, source) -- the variable modifications a search gets, and on whose word. `var_mods`
    is parse_var_mods()'s list (given), a --var-mods string, or None (not given: the default for
    this engine and acquisition, tagged DEFAULT)."""
    if isinstance(var_mods, str):
        var_mods = parse_var_mods(var_mods)
    if var_mods is not None:
        words = list(var_mods)
        if not words:
            return [], ("none: --var-mods none (variable modifications turned off on request)")
        return words, (f"user requested (--var-mods {','.join(words)}): "
                       + ", ".join(f"{VAR_MODS[w]['label']} ({VAR_MODS[w]['diann'].split(',')[0]})"
                                   for w in words))
    acq = (acquisition or "").strip().upper()
    words = list(VAR_MODS_DEFAULT.get((engine, acq), []))
    if words:
        return words, (f"DEFAULT, not user-confirmed: {SRC_VAR_MODS_SAGE_DDA}; "
                       f"--var-mods none turns it off")
    if engine == "diann":
        return [], SRC_VAR_MODS_DIANN
    return [], (f"DEFAULT: Ox(M) off for {acq or 'this acquisition'} (quant default; it is on "
                f"by default only for DDA); --var-mods ox adds it")


# What a precursor m/z range is, by acquisition. detect_acquisition.py reads DIA's from the
# isolation windows and DDA's from the MS1 survey scans, and says which per file
# (`precursor_mz_range_source`); build_diann() tags the cfg's --min/--max-pr-mz with the same words.
PR_MZ_RANGE_BASIS = {"DIA": "isolation_windows", "DDA": "ms1_survey_scan"}
PR_MZ_RANGE_TEXT = {
    "isolation_windows": "measured from the acquired isolation windows",
    "ms1_survey_scan": ("measured from the acquired MS1 survey scan range -- a DDA run picks its "
                        "precursors from its survey scans, so none lies outside it"),
}

# Scan window under --dda. DIA-NN logs NO scan-window radius in DDA mode: not on the 67-run nail
# DDA chain (DIA-NN 2.6.0, HIVE, /quobyte/proteomics-grp/brett/diann26_dda_parallel_nail: "RT window
# set to 3.55907" and "Recommended MS1 mass accuracy setting: 9 ppm", and no radius line in the
# step 2, 3, 4 and 5 logs read), and not in SET28's step 1b probes (DIA-NN 2.7.0), two of which
# waited out their 3600 s for it before the third was cancelled. So nothing can measure it, and
# nothing known before the search can derive it: the flag is left out. The one description, for
# the cfg's rationale and diann_parallel.window_record().
DDA_WINDOW_NOTE = (
    "not set (flag omitted), not probed: DDA (--dda), so DIA-NN chooses. DIA-NN logs no "
    "scan-window radius in DDA mode (nail DDA chain, DIA-NN 2.6.0; SET28 probes, 2.7.0), so the "
    "5-step chain's step 1b cannot measure one and does not run for DDA, and nothing known before "
    "the search can derive one. Whether DIA-NN applies a scan window to DDA data at all is "
    "unverified. EXPECTED for DDA: the chain's steps 3 and 5 print DIA-NN's 'combining reuse of "
    ".quant files with automatic optimisation of mass accuracies or scan window' warning, "
    "because the flag is absent (seen on the nail chain with mass accuracy pinned)")
# How published text (the Methods paragraph and its parameter table, the run record) names a
# mass-accuracy level pinned at the SOP as a DEFAULT: the one wording, read by make_methods.py
# wherever the run's provenance marks a level as one (diann_parallel.mass_acc_defaults).
SOP_DEFAULT_PUBLISHED = "Core SOP default for this instrument, not measured on these data"
# Why a DDA search cannot have its mass accuracy measured, for the tags that say so.
DDA_SOP_WHY = ("no scan-window radius is logged in DDA mode, and mass-accuracy optimisation did "
               "not finish within 3600 s on the SET28 runs (DIA-NN 2.7.0), so step 1b's probe "
               "cannot measure it")


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
#   pinned              -- --mass-acc/--mass-acc-ms1 are in the cfg (a documented value, the
#                          facility SOP's for a DDA level, or a user override)
#   measure_with_diann  -- an Orbitrap with at least one level outside DIA-NN's table; the cfg
#                          omits BOTH flags, and before the search DIA-NN is run on representative
#                          runs (chain step 1b, or run_search.py's single-shot probe). The level
#                          that has a documented tier is recorded under
#                          `mass_accuracy_documented` and pinned as given. A user override of
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


# What a value forced with --overrides is: the USER's choice for this dataset, which may be a
# validated method or a one-off workaround -- the skill cannot tell. Every override was tagged
# "user-override (validated SOP)", so a 25 ppm chosen to cover a timsTOF ~20 ppm off calibration
# went into the provenance, and on to the Methods, as a validated SOP (staff report, 2026-09-25).
# The ONE wording, with the value, who set it and why; resolve_defaults.py --ms1-ppm/--ms2-ppm
# uses it too. Who is said by role: the sidecar and the manifest reach the client
# (reproducibility/inputs/, the session zip). The role is staff.ROLE ("Core staff") only when the
# setter -- --override-by when given, else the login that ran the command -- is on a usable Core
# staff list, else staff.USER_ROLE ("the user"): staff.role_for(), which both main()s call. Every
# override used to say "set by Core staff", whoever set it (S4). The person is in the staff-only
# record (record_override_setter).
OVERRIDE_REASON_MISSING = "reason not given [ask the user -- confirm]"


def override_setter(given=None):
    """(who, how it is known): --override-by as given, else the login that ran this command --
    which is who TYPED it, not necessarily who chose the value, so it says so."""
    if given and str(given).strip():
        return str(given).strip(), "given"
    import getpass
    try:
        return getpass.getuser(), "login"
    except (KeyError, OSError):          # no login name to read (no USER/LOGNAME, no pwd entry)
        return None, "unknown"


def override_source(flag, value, reason=None, role=staff.USER_ROLE):
    """The tag of one overridden setting: a user override, its value, by role, and why. `role` is
    staff.role_for(<the setter>): main() and resolve_defaults.py pass it; a caller that forgets
    gets the less privileged label, staff.USER_ROLE."""
    shown = value if isinstance(value, (int, float, str)) and not isinstance(value, bool) \
        else json.dumps(value)
    return (f"user override: {flag} = {shown}, set by {role}; reason: "
            f"{reason.strip() if reason and reason.strip() else OVERRIDE_REASON_MISSING}")


def override_record(overrides, reason=None, staff_record=None, role=staff.USER_ROLE):
    """{setting: {value, set_by, reason, source, staff_record}} -- the sidecar's `overrides`, by
    role (staff.role_for); `staff_record` names the staff-only file that holds who, by name."""
    return {k: {"value": v, "set_by": role,
                "reason": reason.strip() if reason and reason.strip() else None,
                "source": override_source(k, v, reason, role), "staff_record": staff_record}
            for k, v in (overrides or {}).items()}


def staff_record_name(path):
    """The staff-only record beside `path` (a cfg, a manifest) -- its file name."""
    return os.path.basename(path) + staff.STAFF_SUFFIX


def record_override_setter(path, overrides, by, by_how, reason, role=None):
    """Who set `overrides`, by name and how that is known ("given": --override-by; "login": the
    login that ran the command, which is who TYPED it), into the staff-only record beside `path`,
    with the role the client documents name for them (staff.role_for). -> its file name, for
    override_record()."""
    staff.append(path + staff.STAFF_SUFFIX, {
        "event": "override", "overrides": overrides, "by": by, "by_how": by_how,
        "account": staff.login(), "at": staff.now(), "role": role,
        "reason": reason.strip() if reason and reason.strip() else None})
    return staff_record_name(path)


class LoneMassAccOverride(ValueError):
    """An --overrides mass-accuracy flag whose other level has no value: it could only be written
    as a lone flag (LONE_FLAG_NOTE). references/parameters.md's own example, {"--mass-acc": 8},
    once went into the cfg alone for a 120k/15k Orbitrap while the sidecar still planned to measure
    MS2 -- so MS1 ran at 20 ppm instead of its documented 7, and nothing said so."""


class NoPrecursorRange(ValueError):
    """A DIA-NN DDA cfg with no --precursor-mz-range. It used to fall back to 380-980, and
    SKILL.md told the agent to omit the range for DDA: SET28's Exploris DDA surveyed 350-1500 and
    its FragPipe PSMs spanned 360-1315 m/z, so about a third of the identifiable precursors would
    have been left out of the predicted library, silently."""


NO_DDA_RANGE_MESSAGE = (
    "DIA-NN in DDA mode needs the acquired precursor m/z range and none was given "
    "(--precursor-mz-range). It is not guessed: the old 380-980 fallback would have cut about a "
    "third of SET28's precursors (survey scan 350-1500, FragPipe PSMs 360-1315 m/z) without an "
    "error. Pass --precursor-mz-range <LO> <HI> from step 2's detect_acquisition.py "
    "`precursor_mz_range` -- for a DDA run it is the MS1 survey scan range, read from the Thermo "
    "scan filter ('Full ms [lo-hi]'), the mzML MS1 scan window or a timsTOF .d's MzAcqRange -- or, "
    "when step 2 could not read it, the survey scan range from the instrument method (ask the "
    "user). No cfg was written.")


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


def dda_sop_levels(instr_class, ms1, ms2, acquisition):
    """{flag: ppm} -- the levels a DIA-NN DDA cfg pins at the facility SOP, else {}.

    For DIA these levels (an Orbitrap level outside DIA-NN's table: SET28's 15k MS2) are MEASURED
    with DIA-NN before the search. Under --dda that cannot work: step 1b's probe waits for a
    scan-window radius DIA-NN never logs in DDA mode (DDA_WINDOW_NOTE); SET28's probes each ran to
    their 3600 s timeout -- about 3 h of step 1b, then a failed step 1b and steps 2-5 stuck on
    DependencyNeverSatisfied. So the level is pinned at the SOP -- the floor a measured level would
    have been raised to anyway (SOP_MASS_ACC) -- and tagged as a DEFAULT, not user-confirmed. The
    level that has a tier keeps it, as in DIA. Classes with no usable tier at all (resolution
    unknown, ion-trap MS2) are not covered: the SOP is an Orbitrap-MS2 number, and those stay on
    DIA-NN's own calibration as before."""
    if not is_dda(acquisition) or instr_class not in MEASURE_CLASSES:
        return {}
    whole = lambda v: int(v) if float(v).is_integer() else v   # noqa: E731 -- "20", not "20.0"
    return {flag: whole(SOP_MASS_ACC[key])
            for flag, key, value in (("--mass-acc-ms1", "ms1_ppm", ms1),
                                     ("--mass-acc", "ms2_ppm", ms2))
            if value is None}


def mass_acc_plan(instr_class, ms1, ms2, acquisition="DIA"):
    """pinned | measure_with_diann | auto -- see PLAN_PINNED above. The one place that decides.

    DDA never measures (dda_sop_levels): a level DIA would measure is pinned at the SOP instead."""
    if ms1 and ms2:
        return PLAN_PINNED
    if instr_class in MEASURE_CLASSES:
        return PLAN_PINNED if is_dda(acquisition) else MEASURE_WITH_DIANN
    return PLAN_AUTO


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
# DIA-NN flags the cfg carries once per value (one --var-mod line per modification searched,
# from VAR_MODS): an --overrides value ADDS to those lines (union, deduplicated; a list adds
# several), never replaces them -- {"--var-mod": "UniMod:21,79.966331,STY"} beside --var-mods
# ox,deam used to drop Ox(M) and Deamidation without a word (review F3).
DIANN_ADDED_FLAGS = ("--var-mod",)


def build_diann(acq, instr_class, ms1, ms2, label, src, var_mods, overrides,
                cont_tag=None, mz_range=None, level_src=None, override_meta=None, cys=None,
                labels=None):
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

    dda = is_dda(acq)
    add("--qvalue", 0.01, "standard 1% precursor FDR")
    # In the cfg, not on a command line: every route reads the cfg (the chain's steps 1-5 and
    # step-1 library prediction, the single-shot library and search jobs), and the cfg is what
    # the provenance records and a replay re-runs. Library prediction with --dda is proven: the
    # nail DDA chain's step 1 (DIA-NN 2.6.0) logged "All runs will be analysed as DDA runs" and
    # saved its predicted library.
    if dda:
        add(DIANN_DDA_FLAG, True, SRC_DDA)
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
    add("--cut", DIANN_DIGEST["cut"], "trypsin/P")
    add("--missed-cleavages", DIANN_DIGEST["missed_cleavages"], "DIA-NN default")
    add("--min-pep-len", DIANN_DIGEST["min_pep_len"], UNIV)
    add("--max-pep-len", DIANN_DIGEST["max_pep_len"], UNIV)
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
    #
    # DDA has no isolation-window scheme: its range is the MS1 SURVEY scan range, which
    # detect_acquisition.py reads for a DDA run instead (PR_MZ_RANGE_BASIS). With none, a DDA
    # cfg is refused, not given the fallback below -- see NoPrecursorRange.
    if mz_range:
        lo, hi = float(mz_range[0]), float(mz_range[1])
        what = PR_MZ_RANGE_TEXT[PR_MZ_RANGE_BASIS["DDA" if dda else "DIA"]]
        add("--min-pr-mz", int(math.floor(lo)), f"{what} ({lo:.1f}-{hi:.1f} m/z)")
        add("--max-pr-mz", int(math.ceil(hi)), f"{what} ({lo:.1f}-{hi:.1f} m/z)")
    elif dda:
        raise NoPrecursorRange(NO_DDA_RANGE_MESSAGE)
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
    if DIANN_DIGEST["met_excision"]:
        add("--met-excision", True, UNIV)
    # cysteine: what the sample prep did to it (CYS_REAGENTS), tagged with who said so
    cys = cys or cys_alkylation_record()
    cys_flags = cys_diann_flags(cys["value"])
    if cys["value"] == "cam":
        add("--unimod4", True, "fixed carbamidomethyl (C); DIA-NN recommended fixed mod -- "
            + cys["source"])
    elif not cys_flags:
        add("--unimod4", "not set (flag omitted): cysteine searched unmodified", cys["source"],
            render=False)
    else:
        add("--fixed-mod", cys_flags[1], f"fixed {cys['mod']} (C) -- {cys['source']}")

    # variable mods (VAR_MODS): DIA-NN README says var mods don't help pure quant; default OFF.
    # --var-mods is DIA-NN's maximum per peptide: one for each kind asked for.
    vm_words, vm_src = var_mods_choice("diann", acq, var_mods)
    if vm_words:
        names = ", ".join(VAR_MODS[w]["label"] for w in vm_words)
        add("--var-mods", len(vm_words), f"one per variable-modification kind searched "
                                         f"({names}) -- {vm_src}")
        r["--var-mod"] = tagged([VAR_MODS[w]["diann"] for w in vm_words], f"{names} -- {vm_src}")
        lines.extend(f"--var-mod {VAR_MODS[w]['diann']}" for w in vm_words)
        r["variable_mods"] = tagged(vm_words, vm_src)
    else:
        r["variable_mods"] = tagged("none", vm_src)

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
    plan = mass_acc_plan(instr_class, ms1, ms2, acq)
    # DDA: the levels DIA would measure, pinned at the SOP instead ({} for DIA) -- dda_sop_levels
    sop_default = dda_sop_levels(instr_class, ms1, ms2, acq)
    documented, basis = {}, {}
    # A user override of a mass-accuracy level (applied below, with the other overrides) is a
    # value for that level. Both levels known -- overridden, or from the table -- means both are
    # written and nothing is measured. A one-flag override whose other level has NO value would
    # be a lone flag (LONE_FLAG_NOTE: the other level silently fixed at 20 ppm), so it is refused.
    given = [f for f in MASS_ACC_FLAGS if f in (overrides or {})]
    if given:
        table = {"--mass-acc-ms1": ms1, "--mass-acc": ms2}
        table.update(sop_default)
        missing = [f for f in MASS_ACC_FLAGS if f not in given and table[f] is None]
        if missing:
            raise LoneMassAccOverride(lone_override_message(given[0], missing[0], instr_class,
                                                            label))
        plan = PLAN_PINNED            # the user's value wins; it must not be re-measured
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
        for flag, level, value, lsrc in (("--mass-acc", "MS2", ms2, s2),
                                         ("--mass-acc-ms1", "MS1", ms1, s1)):
            if flag in given:
                continue
            if flag in sop_default:
                # rule 2: a value nobody confirmed, in a cfg every export describes
                add(flag, sop_default[flag],
                    f"{label}: {level} {sop_default[flag]:g} ppm -- DEFAULT, not user-confirmed: "
                    f"{lsrc}, and a DDA search (--dda) cannot measure it ({DDA_SOP_WHY}), so it "
                    f"is pinned at {SOP_MASS_ACC_SOURCE}. Set it with --overrides if the "
                    "instrument needs another value")
            else:
                add(flag, value, f"{label}: {level} {value} ppm [{lsrc}]")
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
        DDA_WINDOW_NOTE if dda else
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

    # SILAC (--labels): the plexDIA declaration, from silac.py -- the one definition of it
    if labels:
        for flag, value, source, render in silac.diann_cfg_lines(labels):
            if flag in r:
                # a flag the cfg already carries -- --fixed-mod for a cysteine agent other than
                # iodoacetamide (CYS_REAGENTS): both lines are written, each with its own rationale
                r[f"{flag} (SILAC label)"] = tagged(value, source)
                if render:
                    lines.append(f"{flag} {value}")
            else:
                add(flag, value, source, render=None if render else False)

    # apply the user's overrides -- re-render the file; each tagged with who set it and why. A
    # flag replaces its line, except a repeatable one (DIANN_ADDED_FLAGS), which adds.
    om = override_meta or {}
    added, added_src = [], []
    for k, v in (overrides or {}).items():
        src_o = override_source(k, v, om.get("reason"), om.get("role", staff.USER_ROLE))
        if k in DIANN_ADDED_FLAGS:
            have = [ln[len(k) + 1:] for ln in lines if ln.startswith(k + " ")]
            new = []
            for x in (v if isinstance(v, (list, tuple)) else [v]):
                if str(x) not in have + new:
                    new.append(str(x))
            lines.extend(f"{k} {x}" for x in new)
            prev = r.get(k)
            r[k] = tagged(have + new, (f"{prev['source']}; + " if prev else "") + src_o)
            added += new
            added_src.append(src_o)
            continue
        r[k] = tagged(v, src_o)
        lines = [ln for ln in lines if not (ln == k or ln.startswith(k + " "))]
        lines.append(k if v is True else f"{k} {v}")
    if added:
        # the modifications searched are now the table's AND the added ones: the rationale and
        # --var-mods (one per kind, unless the user set it) follow the final lines
        final = r["--var-mod"]["value"]
        names = ", ".join([VAR_MODS[w]["label"] for w in vm_words] + added)
        r["variable_mods"] = tagged(list(vm_words) + added,
                                    (f"{vm_src}; + " if vm_words else "") + "; ".join(added_src))
        if "--var-mods" not in (overrides or {}):
            lines = [ln for ln in lines if not ln.startswith("--var-mods ")]
            lines.append(f"--var-mods {len(final)}")
            # with no table mods, vm_src is why there were none ("left off"): not this one's
            r["--var-mods"] = tagged(len(final), f"one per variable-modification kind searched "
                                                 f"({names})"
                                                 + (f" -- {vm_src}; +" if vm_words else ":")
                                                 + " --var-mod added in --overrides")
    # The DDA levels pinned at the SOP, by flag -- what an export must call a default. A level
    # the caller overrode is theirs, not a default, so it is not listed.
    defaulted = {f: v for f, v in sop_default.items() if f not in given}
    r["mass_accuracy_plan"] = tagged(plan, {
        PLAN_PINNED: "--mass-acc/--mass-acc-ms1 are in the cfg" + (
            "; " + ", ".join(f"{f} {v:g}" for f, v in sorted(defaulted.items()))
            + " pinned at the facility SOP -- DEFAULT, not user-confirmed (a DDA search cannot "
              "measure it: " + DDA_SOP_WHY + ")" if defaulted else ""),
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
    r["mass_accuracy_default"] = tagged(
        defaulted,
        (f"DEFAULT, not user-confirmed: DDA levels with no DIA-NN resolution-table value, pinned "
         f"at {SOP_MASS_ACC_SOURCE}, because {DDA_SOP_WHY}") if defaulted else
        "none: only a DDA cfg for an Orbitrap level outside DIA-NN's table pins one")

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


# Sage's in-silico digest (its database.enzyme block): trypsin with the proline rule, 2 missed
# cleavages. The ONE definition: build_sage writes it, and run_search.py recognises it as Sage's
# KNOWN difference from DIANN_DIGEST (which fetch_fasta.py judges contaminants on).
SAGE_ENZYME = {"missed_cleavages": 2, "min_len": 7, "max_len": 30,
               "cleave_at": "KR", "restrict": "P"}


# --ids-only: a Sage search for identifications only, with no LFQ (no differential expression).
# The cfg's quant.lfq is false, which is the ONE carrier of the choice: run_search.py reads it
# (sage_lfq_enabled) to send a timsTOF .d to Sage directly and to skip the LFQ check. A DIA-NN
# search always quantifies, so the flag is refused for it (a Core submission, 2026-10-29).
IDS_ONLY_FLAG = "--ids-only"
SRC_IDS_ONLY = ("identifications only (--ids-only: no differential expression); Sage 0.14.7 reads "
                "a timsTOF .d directly but has no MS1 quantification on it")
IDS_ONLY_DIANN = "--ids-only is for Sage; a DIA-NN search always quantifies"


# --overrides into a Sage cfg: the ONE merge (H2-4). It took top-level keys only, so an
# enzyme-only override, {"database": {"enzyme": {...}}}, REPLACED the whole database block and the
# fixed carbamidomethyl-C, the protein N-term acetyl, the decoys and the FASTA placeholder were
# gone, with no message. A Sage tolerance is a one-key enum ({"ppm": [lo, hi]} or {"da": [lo, hi]}:
# crates/sage/src/mass.rs `enum Tolerance { Ppm(f32, f32), Da(f32, f32) }`), so precursor_tol and
# fragment_tol are replaced whole: merged, {"da": ...} over a ppm window would hold both units.
SAGE_REPLACED_WHOLE = (("precursor_tol",), ("fragment_tol",))

# Sage 0.14.7's settings, by block: crates/sage-cli/src/input.rs (`Input`, `QuantOptions`,
# `LfqOptions`, `TmtOptions`) and crates/sage/src/database.rs (`Builder`, `EnzymeBuilder`), read at
# tag v0.14.7 (2026-11-03). None of them denies unknown fields, so Sage IGNORES a key it does not
# know, without a word: an override at the wrong level (the rationale names its entries flat --
# enzyme, static_mods, lfq -- and {"enzyme": ...} at the top level is the easy mistake) was
# recorded as applied and never searched (review F2). static_mods / variable_mods are maps of
# residues, so their keys are not checked here.
SAGE_CFG_KEYS = {
    (): {"database", "precursor_tol", "fragment_tol", "report_psms", "chimera", "wide_window",
         "min_peaks", "max_peaks", "max_fragment_charge", "min_matched_peaks", "precursor_charge",
         "isotope_errors", "deisotope", "quant", "predict_rt", "output_directory", "mzml_paths",
         "annotate_matches", "write_pin"},
    ("database",): {"bucket_size", "enzyme", "fragment_min_mz", "fragment_max_mz",
                    "peptide_min_mass", "peptide_max_mass", "ion_kinds", "min_ion_index",
                    "static_mods", "variable_mods", "max_variable_mods", "decoy_tag",
                    "generate_decoys", "fasta"},
    ("database", "enzyme"): {"missed_cleavages", "min_len", "max_len", "cleave_at", "restrict",
                             "c_terminal", "semi_enzymatic"},
    ("quant",): {"tmt", "tmt_settings", "lfq", "lfq_settings"},
    ("quant", "lfq_settings"): {"peak_scoring", "integration", "spectral_angle", "ppm_tolerance",
                                "combine_charge_states"},
    ("quant", "tmt_settings"): {"level", "sn"},
}


class OverrideRefused(ValueError):
    """An --overrides value the cfg cannot carry as asked: refused (refuse_out: no cfg, no
    sidecar), so it is never recorded as applied."""


class EmptyOverride(OverrideRefused):
    """An override of {} on a block: merged key by key it changes nothing, so it is refused
    rather than dropped without a word."""


def _dotted(path):
    return ".".join(str(p) for p in path)


def check_sage_overrides(over, _at=()):
    """Refuse (OverrideRefused) a Sage --overrides key that is not a Sage 0.14.7 setting at its
    level (SAGE_CFG_KEYS) -- naming where a key of that name belongs -- and a plain value given
    for a block. Nothing else is judged here: the values go to Sage as given."""
    known = SAGE_CFG_KEYS.get(_at)
    if known is None or not isinstance(over, dict):
        return
    for k, v in over.items():
        at = _at + (k,)
        if at in SAGE_REPLACED_WHOLE and not _one_window(v):
            # review F5: {} (or two units, or one number) wrote a cfg Sage cannot read, and the
            # job failed only once it had queued
            raise OverrideRefused(
                f"--overrides sets {_dotted(at)} to {json.dumps(v)}: a Sage tolerance is exactly "
                f'one window, {{"ppm": [lo, hi]}} or {{"da": [lo, hi]}}, with lo <= hi '
                f"(crates/sage/src/mass.rs `enum Tolerance`).")
        if k not in known:
            where = sorted(_dotted(lvl + (k,)) for lvl, keys in SAGE_CFG_KEYS.items()
                           if k in keys and lvl != _at)
            raise OverrideRefused(
                f"--overrides sets {_dotted(at)}, which is not a Sage 0.14.7 setting at that level "
                + (f"-- it belongs at {' or '.join(where)} (write it nested: "
                   f"{json.dumps(_nested(where[0].split('.'), v))})" if where else
                   "(crates/sage-cli/src/input.rs, crates/sage/src/database.rs)")
                + ". Sage would ignore it without a word, and the search would run without it.")
        if at in SAGE_CFG_KEYS:
            if not isinstance(v, dict):
                raise OverrideRefused(
                    f"--overrides sets {_dotted(at)} to {json.dumps(v)}, but {_dotted(at)} is a "
                    f"block of settings ({', '.join(sorted(SAGE_CFG_KEYS[at]))}): give the keys "
                    f"to change, e.g. {{\"{k}\": {{\"<key>\": <value>}}}}.")
            check_sage_overrides(v, at)


def _one_window(tol):
    """True for a Sage Tolerance: one key, "ppm" or "da", holding two numbers lo <= hi."""
    if not (isinstance(tol, dict) and len(tol) == 1 and next(iter(tol)) in ("ppm", "da")):
        return False
    w = next(iter(tol.values()))
    return (isinstance(w, (list, tuple)) and len(w) == 2
            and all(isinstance(x, (int, float)) and not isinstance(x, bool) for x in w)
            and w[0] <= w[1])


def _nested(path, value):
    for k in reversed(path):
        value = {k: value}
    return value


def deep_merge(base, over, whole=(), _at=()):
    """`over` merged into a copy of `base` (neither is changed): a dict that lands on a dict merges
    into it key by key, at any depth; anything else replaces -- a number, a string, a list (a
    list is replaced, never appended), a dict that lands on a non-dict or on a path in `whole`.
    An EMPTY dict that lands on a dict would change nothing: EmptyOverride.
    -> (merged, [(path, value), ...]): each setting `over` set, at the depth it was merged to."""
    import copy
    merged, set_ = copy.deepcopy(base), []
    for k, v in (over or {}).items():
        at = _at + (k,)
        if v == {} and isinstance(merged.get(k), dict) and merged[k] and at not in whole:
            raise EmptyOverride(
                f"--overrides sets {'.'.join(str(p) for p in at)} to {{}}: a block merges key by "
                f"key, so {{}} would change nothing (it holds {json.dumps(merged[k])}). Give the "
                f"keys to change; to search cysteine unmodified use --cys-alkylation none, and "
                f"--var-mods for the variable modifications.")
        if isinstance(v, dict) and isinstance(merged.get(k), dict) and at not in whole:
            merged[k], below = deep_merge(merged[k], v, whole, at)
            set_ += below
        else:
            merged[k] = copy.deepcopy(v)
            set_.append((at, copy.deepcopy(v)))
    return merged, set_


# The cfg value each of build_sage's rationale entries describes, so an override that lands on it,
# or under it, is named in that entry. The tolerance entries hold the NUMBER a window was derived
# from (10 for {"ppm": [-10, 10]}): an override of the window drops them (SAGE_DERIVED_NUMBERS).
# The others hold the cfg's value itself, and keep it as merged.
SAGE_RATIONALE_PATHS = {
    "precursor_tol_ppm": ("precursor_tol",), "fragment_tol_ppm": ("fragment_tol",),
    "fragment_tol_da": ("fragment_tol",), "deisotope": ("deisotope",),
    "bucket_size": ("database", "bucket_size"), "wide_window": ("wide_window",),
    "chimera": ("chimera",), "static_mods": ("database", "static_mods"),
    "variable_mods": ("database", "variable_mods"), "enzyme": ("database", "enzyme"),
    "lfq": ("quant", "lfq"),
}
SAGE_DERIVED_NUMBERS = ("precursor_tol_ppm", "fragment_tol_ppm", "fragment_tol_da")


def _cfg_at(cfg, path):
    """The value at `path` in a nested cfg, else None."""
    for k in path:
        cfg = cfg.get(k) if isinstance(cfg, dict) else None
    return cfg


def build_sage(acq, instr_class, var_mods, overrides, ms2_analyzer=None, override_meta=None,
               cys=None, ids_only=False):
    prec_ppm, frag_ppm, ppm_src = sage_ppm(instr_class)
    cys = cys or cys_alkylation_record()
    static = cys_sage_static(cys["value"])
    ion_trap = ms2_in_ion_trap(ms2_analyzer)
    UNIV = "universal trypsin/LFQ default"
    # VAR_MODS: what was asked for, else this acquisition's default (Ox(M) on for DDA)
    vm_words, vm_src = var_mods_choice("sage", acq, var_mods)
    variable = {}
    for w in vm_words:
        for site, masses in VAR_MODS[w]["sage"].items():
            have = variable.setdefault(site, [])
            have.extend(m for m in masses if m not in have)
    variable["["] = [42.0106]  # protein N-term acetyl is a common, cheap variable mod
    cfg = {
        "database": {
            "bucket_size": SAGE_ITMS_BUCKET_SIZE if ion_trap else 8192,
            "enzyme": dict(SAGE_ENZYME),
            "fragment_min_mz": 200.0, "fragment_max_mz": 1800.0,
            "peptide_min_mass": 500.0, "peptide_max_mass": 5000.0,
            "ion_kinds": ["b", "y"], "min_ion_index": 2,
            "static_mods": static,
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
        "quant": {"lfq": not ids_only},
        "output_directory": "sage_out",
    }
    # the user's overrides: each a Sage setting at its level (check_sage_overrides), merged at
    # every depth (deep_merge); `set_` is what each one set
    check_sage_overrides(overrides)
    built = cfg                        # before the overrides: what each derived source describes
    cfg, set_ = deep_merge(cfg, overrides, SAGE_REPLACED_WHOLE)
    db = cfg.get("database") if isinstance(cfg.get("database"), dict) else {}

    # Every entry that describes a cfg value reads it from the MERGED cfg, never the pre-override
    # one (H2-4). The enzyme entry is the cfg's own database.enzyme block: it was the text
    # "trypsin/P, 2 missed cleavages" whatever the cfg held -- and SAGE_ENZYME's restrict "P" is
    # trypsin, not trypsin/P.
    rationale = {"precursor_tol_ppm": tagged(prec_ppm, ppm_src)}
    if ion_trap:
        it = ("MS2 read partly in the ion trap (mixed ITMS/FTMS)" if ms2_analyzer == "mixed"
              else "MS2 read in the ion trap (ITMS)")
        rationale.update({
            "fragment_tol_da": tagged(SAGE_ITMS_FRAGMENT_DA, f"{it}: {SRC_SAGE_LOWRES}"),
            "deisotope": tagged(cfg.get("deisotope"), f"{it}: {SRC_SAGE_LOWRES}"),
            "bucket_size": tagged(db.get("bucket_size"), f"{it}: {SRC_SAGE_LOWRES} ('higher "
                                                         "values for low-res MS/MS')"),
        })
    else:
        rationale["fragment_tol_ppm"] = tagged(frag_ppm, ppm_src)
    rationale.update({
        "wide_window": tagged(cfg.get("wide_window"), f"{acq.upper()} acquisition"),
        "chimera": tagged(cfg.get("chimera"), f"{acq.upper()} acquisition"),
        "static_mods": tagged(db.get("static_mods"),
                              (f"fixed {cys['mod'].lower()} (C) -- " if static else
                               "none: cysteine searched unmodified -- ") + cys["source"]),
        "variable_mods": tagged(db.get("variable_mods"),
                                vm_src + "; + protein N-term acetyl (always)"),
        "enzyme": tagged(db.get("enzyme"), UNIV),
        "lfq": tagged(_cfg_at(cfg, ("quant", "lfq")),
                      SRC_IDS_ONLY if ids_only else "label-free quantification"),
    })
    # The user's overrides, each tagged with who set it and why, under the path deep_merge set it
    # at ("precursor_tol", "database.enzyme.cleave_at"). An entry it lands on or under keeps the
    # merged value with the override's source instead of its derived one; a derived number whose
    # window it replaced (precursor_tol_ppm under an overridden precursor_tol) goes. They used to
    # be merged into the config with no rationale entry at all.
    # An override of PART of an entry (one residue of static_mods, the enzyme's missed cleavages)
    # keeps the entry's own source, DEFAULT tag included, and adds its own (review F4, rule 2):
    # the rest of the value is still what that source says. One that replaces the whole value
    # says only its own.
    om = override_meta or {}
    changed = {}                    # entry -> {"orig": its source, "srcs": [...], "whole": bool}
    for path, v in set_:
        k = ".".join(str(p) for p in path)
        src = override_source(k, v, om.get("reason"), om.get("role", staff.USER_ROLE))
        for x in [x for x in rationale if x in SAGE_RATIONALE_PATHS]:
            xp = SAGE_RATIONALE_PATHS[x]
            n = min(len(xp), len(path))
            if tuple(xp[:n]) != tuple(path[:n]):
                continue
            if x in SAGE_DERIVED_NUMBERS:
                del rationale[x]
            else:
                c = changed.setdefault(x, {"orig": rationale[x]["source"], "srcs": [],
                                           "whole": False, "keys": set()})
                c["srcs"].append(src)
                c["whole"] = c["whole"] or len(path) <= len(xp)
                if len(path) > len(xp):
                    c["keys"].add(path[len(xp)])
        rationale[k] = tagged(v, src)
    for x, c in changed.items():
        if x in rationale:
            # every key the derived source described was replaced ({"C": 0} on a static_mods
            # that held only C): nothing of it is left to that source (review R2)
            was = _cfg_at(built, SAGE_RATIONALE_PATHS[x])
            if isinstance(was, dict) and was and set(was) <= c["keys"]:
                c["whole"] = True
            said = "; ".join(c["srcs"])
            rationale[x] = tagged(_cfg_at(cfg, SAGE_RATIONALE_PATHS[x]),
                                  said if c["whole"] else f"{c['orig']}; + {said}")
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


# The flags a SILAC search is declared with: set by --labels, never by --overrides. As overrides
# (a triple-SILAC search, 2026-10-14) they were recorded as one person's override rather than
# the study design, and --lib-fixed-mod went in with them and labelled the library twice.
# --channels and the channel normalisations exist only for plexDIA: refused as overrides always.
# --original-mods and --lib-fixed-mod also have label-free uses: refused only beside SILAC (with
# --labels, --channels or a label --fixed-mod in the same request).
SILAC_OVERRIDE_FLAGS = ("--channels", *silac.CHANNEL_NORM.values())
SILAC_SHARED_FLAGS = ("--lib-fixed-mod", "--original-mods")


def refuse_out(out, msg, code=1):
    """Exit with `msg`, leaving nothing at --out. The message says no cfg was written, and THIS run
    wrote none -- but an earlier run of the same command may have left the cfg and its sidecar,
    and the next step reads them by path, not by mtime. A stale cfg beside a refusal on stdout is
    the worst of both. Every refusal that says "no cfg was written" goes through here.

    The earlier cfg is never DELETED (C2 errors_ux L1: a mistyped --var-mods word deleted the cfg
    the session's first search ran with, which the Methods and the bundle read): it and its
    sidecar are set aside (diann_parallel.set_aside, the one definition: renamed .stale-<time>),
    the message names where, and says so when one could not be moved."""
    moved, kept = [], []
    for stale in (out, out + ".rationale.json"):
        if not os.path.lexists(stale):
            continue
        try:
            import diann_parallel
            new = diann_parallel.set_aside(stale)
            if new:
                moved.append(new)
        except (ImportError, ValueError, OSError) as e:
            kept.append(f"{stale} ({e})")
    if moved:
        msg += (" The cfg an earlier run left at --out was NOT deleted: set aside as "
                + " and ".join(moved) + " (rename it back to use it again).")
    if kept:
        msg += (" WARNING: could not set aside " + "; ".join(kept) + " -- it is still at --out "
                "and does NOT have what this command asked for.")
    print(f"[estimate_params] {msg}", file=sys.stderr)
    sys.exit(code)


def silac_labels(a, overrides):
    """--labels/--label-preset/--channel-norm -> silac.labels_record(), or None. Exits with the
    reason for a request this skill cannot search: not DIA-NN DIA, a preset without --labels, or
    SILAC flags passed as --overrides. Like the other refusals, it leaves no cfg at --out from an
    earlier run to contradict "no cfg was written"."""
    def refuse(msg):
        refuse_out(a.out, f"{msg} No cfg was written.")
    label_mod = lambda k: k == "--fixed-mod" and ",label" in str(overrides[k])  # noqa: E731
    silac_request = bool(a.labels) or any(k in SILAC_OVERRIDE_FLAGS or label_mod(k)
                                          for k in overrides or {})
    pass_as = [k for k in overrides or {} if k in SILAC_OVERRIDE_FLAGS or label_mod(k)
               or (k in SILAC_SHARED_FLAGS and silac_request)]
    if pass_as:
        refuse(f"{', '.join(pass_as)} in --overrides: SILAC is declared with --labels "
               "silac2|silac3 --label-preset ... [--channel-norm spec|run], which records it as "
               "the study design and declares the label once (references/silac.md).")
    if a.labels and "--fixed-mod" in (overrides or {}):
        refuse("--fixed-mod in --overrides with --labels: it would replace the SILAC label's "
               "--fixed-mod line. A cysteine agent goes in --cys-alkylation.")
    if not a.labels:
        if a.label_preset or a.channel_norm:
            refuse("--label-preset/--channel-norm describe --labels, and it was not given.")
        return None
    if a.engine != "diann" or is_dda(a.acquisition):
        refuse("--labels is DIA-NN plexDIA: DIA data searched with --engine diann. SILAC with "
               "Sage or DIA-NN's DDA mode is not supported (references/silac.md).")
    try:
        return silac.labels_record(a.labels, a.label_preset, a.channel_norm)
    except silac.SilacError as e:
        refuse(str(e).rstrip(".") + ".")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--engine", required=True, choices=["diann", "sage"])
    ap.add_argument("--acquisition", required=True, choices=["DIA", "DDA", "dia", "dda"])
    ap.add_argument("--instrument", default="")
    ap.add_argument("--var-mods", default=None,
                    help=f"variable modifications, a comma list of: {VAR_MODS_WORDS}. ox = Ox(M), "
                         "deam = Deamidation (NQ), none = off. Any other word is refused. Not "
                         "given: Ox(M) for a Sage DDA search, nothing otherwise (DIA-NN advises "
                         "against variable mods for quantification) -- recorded as a DEFAULT. "
                         "Sage always adds protein N-term acetyl. Giving --var-mods replaces the "
                         "default; to add deamidation and keep Met-ox, give ox,deam")
    ap.add_argument(IDS_ONLY_FLAG, action="store_true",
                    help="Sage only: identifications only, no differential expression. The cfg's "
                         "quant.lfq is false, so run_search.py searches a timsTOF .d directly "
                         "(Sage 0.14.7 reads it, but has no MS1 quantification on it). Refused "
                         "with --engine diann (exit 2): a DIA-NN search always quantifies")
    ap.add_argument("--cys-alkylation", choices=CYS_ALKYLATION_CHOICES, default=None,
                    help="what the sample prep did to cysteine. cam (the default): reduced and "
                         "alkylated with iodoacetamide or chloroacetamide -- fixed carbamidomethyl "
                         "C (DIA-NN --unimod4, Sage static_mods C). none: not alkylated (often so "
                         "for on-bead, IP and pull-down digests) -- cysteine searched unmodified. "
                         "nem / carboxymethyl / propionamide / mmts / pyridylethyl: alkylated with "
                         "N-ethylmaleimide / iodoacetic acid / acrylamide / MMTS / 4-vinylpyridine "
                         "-- that fixed modification on C (UniMod). Not given, cam is used and "
                         "recorded as a DEFAULT")
    ap.add_argument("--labels", choices=sorted(silac.SCHEMES), default=None,
                    help="SILAC multiplexing (DIA-NN plexDIA, DIA only): silac2 = L/H, silac3 = "
                         "L/M/H on Lys and Arg. Needs --label-preset; references/silac.md")
    ap.add_argument("--label-preset", default=None,
                    help="the labelled channels' Lys/Arg labels, lightest first, e.g. K4R6,K8R10 "
                         "(silac3) or K8R10 (silac2); `silac.py presets` lists them")
    ap.add_argument("--channel-norm", choices=sorted(silac.CHANNEL_NORM), default=None,
                    help="with --labels: spec = independent samples in each channel (the "
                         "default, recorded as one); run = old and new protein of one culture "
                         "(pulsed SILAC / turnover)")
    ap.add_argument("--overrides", default="",
                    help="JSON of fields to force, as the user chose them; each is recorded as a "
                         "user override with its value, who set it and why -- never as an SOP")
    ap.add_argument("--override-by", default=None,
                    help="who chose the --overrides values (default: the login running this, "
                         "recorded as such) -- in the staff-only <out>.staff.json only; the "
                         "sidecar says Core staff when that login is on the Core staff list "
                         "(staff.py), else the user")
    ap.add_argument("--override-reason", default=None,
                    help="why, in the user's words (e.g. 'instrument ~20 ppm off calibration; "
                         "user chose 25 ppm'); without it the record says the reason was not "
                         "given")
    add_resolution_args(ap)
    ap.add_argument("--from-mzml", default="",
                    help="read both resolutions straight out of this mzML")
    ap.add_argument("--precursor-mz-range", nargs=2, type=float, metavar=("LO", "HI"),
                    default=None,
                    help="ACQUIRED precursor m/z bounds, from detect_acquisition.py's "
                         "precursor_mz_range (DIA: the isolation windows; DDA: the MS1 survey "
                         "scan range). Without it a DIA range falls back to 380-980 and is "
                         "tagged FALLBACK -- which silently discards anything the method "
                         "acquired outside that window -- and a DIA-NN DDA cfg is refused.")
    ap.add_argument("--fasta-meta", default="",
                    help="fetch_fasta.py's <fasta>.meta.json — supplies the contaminant "
                         "tag so DIA-NN excludes contaminants from quant")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    # --var-mods: an unknown word stops here, before anything is written (VAR_MODS). Blank is not
    # given (so is a list with no word in it, "," ): the default for this engine and acquisition applies.
    try:
        var_mods = parse_var_mods(a.var_mods)
    except ValueError as e:
        refuse_out(a.out, f"{e} No cfg was written.")
    if a.ids_only and a.engine != "sage":
        refuse_out(a.out, f"{IDS_ONLY_DIANN}. No cfg was written.", code=2)

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
    if a.ids_only and "quant" in overrides:
        refuse_out(a.out, f"{IDS_ONLY_FLAG} sets the Sage cfg's quant block (lfq false), and "
                          f"--overrides sets quant too: give one or the other. No cfg was written.")
    if (a.override_by or a.override_reason) and not overrides:
        sys.exit("[estimate_params] --override-by/--override-reason describe --overrides, and "
                 "none were given")
    by, by_how = override_setter(a.override_by) if overrides else (None, "unknown")
    # "Core staff" only for a login on the Core staff list, else "the user" (S4)
    role = staff.role_for(by) if overrides else staff.USER_ROLE
    override_meta = {"reason": a.override_reason, "role": role}
    labels = silac_labels(a, overrides)

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
    level_src = ((ppm_for_resolution(r1)[1], ppm_for_resolution(r2)[1])
                 if cls in ("orbitrap_measured", "orbitrap_untabled") else None)

    cys = cys_alkylation_record(a.cys_alkylation)
    if a.engine == "diann":
        try:
            text, rationale = build_diann(a.acquisition, cls, ms1, ms2, label, src, var_mods,
                                          overrides, cont_tag, a.precursor_mz_range,
                                          level_src=level_src, override_meta=override_meta,
                                          cys=cys, labels=labels)
        except (LoneMassAccOverride, NoPrecursorRange) as e:
            # build_diann() raises BEFORE the write below; a stale cfg with a stale
            # mass_accuracy_plan must not outlive the refusal (refuse_out)
            refuse_out(a.out, str(e))
    else:
        try:
            text, rationale = build_sage(a.acquisition, cls, var_mods, overrides, analyzer,
                                         override_meta=override_meta, cys=cys,
                                         ids_only=a.ids_only)
        except OverrideRefused as e:
            refuse_out(a.out, f"{e} No cfg was written.")

    # a new session's wf/ may not exist yet: make --out's folder rather than crash on it
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write(text)
    # The cysteine record says what the WRITTEN cfg searches (review R1): a DIA-NN --fixed-mod
    # override replaces every --fixed-mod line, the --cys-alkylation one included, and a Sage
    # static_mods C override replaces its C -- read back by cys_searched(), the one reader
    if overrides:
        cys = cys_record_as_written(cys, cys_searched(a.out), override_meta)
    # who set the overrides, by name: staff-only, beside the cfg (never in the sidecar)
    staff_rec = (record_override_setter(a.out, overrides, by, by_how, a.override_reason, role)
                 if overrides else None)

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
        # DDA levels pinned at the facility SOP because they cannot be measured under --dda:
        # DEFAULT, not user-confirmed ({flag: ppm}; {} when nothing was defaulted)
        "mass_accuracy_default": (rationale.get("mass_accuracy_default") or {}).get("value"),
        "params_file": os.path.abspath(a.out),
        # each --overrides value: who set it and why ({} when nothing was overridden)
        "overrides": override_record(overrides, a.override_reason, staff_rec, role=role),
        # what was searched for cysteine (fixed carbamidomethyl or unmodified), and whether the
        # user said so or it is the DEFAULT -- read by make_methods.py
        "cys_alkylation": cys,
        # SILAC channels, their masses and the channel-normalisation answer (silac.labels_record);
        # absent for a label-free search -- read by make_methods.py, run_silac.R and the audit
        **({"labels": labels} if labels else {}),
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
