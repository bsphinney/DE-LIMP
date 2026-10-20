#!/usr/bin/env python3
"""
skyline_peaks.py -- the peak check after DE (SKILL.md step 8f; references/skyline.md): Skyline
draws DIA-NN's peaks for the top proteins. And, for Core staff on request, a targeted follow-up
list (PRM/SRM) from the same Skyline document.

WHAT IT IS, AND WHAT IT IS NOT
  Skyline (SkylineCmd, headless under Wine in ProteoWizard's vendor-reader container) extracts each
  target peptide's chromatograms from the raw files again, from DIA-NN's own spectral library and
  the search's digest settings. The pictures are a visual check of DIA-NN's peaks. Skyline reuses
  DIA-NN's peak boundaries from the library by default (its "use explicit peak bounds" library
  setting is on, and SkylineCmd has no switch for it), so a peak that matches is NOT independent
  confirmation. The real warning signs are the precursor-runs where Skyline picked a different
  peak: its apex more than APEX_DIFF_S seconds from DIA-NN's, or outside DIA-NN's peak edges, or no
  peak at all (disagreements.tsv). report_section() says exactly that, and client_text.py's
  peak_check_overclaims() flags "confirmed by Skyline", "validated by Skyline" and the like.

USAGE
  python3 skyline_peaks.py plan --session <S> [--top N] [--name NAME] [--sif PATH] [--python PY]
          [--ppm-ms1 X --ppm-ms2 Y --ppm-source TEXT] [--max-parallel K]
          [--partition P --account A --qos Q]
          [--search DIR] [--de-dir DIR] [--manifest FILE] [--conditions CSV] [--out DIR]
      Checks what the step needs and writes three SLURM jobs and submit.sh into
      <out>/skyline_work/ (out: <S>/output/skyline). Light: run it where jobs are written (a HIVE
      login node is fine). Exit 0: jobs written, "submit with:" says how. Exit 3: the step cannot
      run here -- the record (<out>/skyline_provenance.json) says why, and the report, the
      Methods and MANIFEST.txt say so. Exit 2: usage.
  python3 skyline_peaks.py skip --session <S> --reason '<why, in the user's words>'
      The user turned the step off: the record says so, and so does the report.
  python3 skyline_peaks.py status --session <S>          the record, summarised (JSON)
  python3 skyline_peaks.py prm --session <S> --instrument <key>[,<key>...] [--same-gradient]
          [--run-length MIN]
      Core staff only, on request: a follow-up PRM/SRM assay from the peak-check document
      (PRM_INSTRUMENTS lists the keys). Writes one job; its lists go to output/skyline_followup/.
  Inside the jobs only (never by hand): build, import, join, prm-run.

THE JOBS (each a thin bash wrapper around one subcommand here)
  1 build   verify the container (the Core's .sif by its pinned sha256; any other .sif by
            SkylineCmd --version), choose the proteins (select_targets: the ONE rule), cut the
            FASTA to them, read DIA-NN's rows for them, then ONE SkylineCmd call builds the
            document (build_args).
  2 import  a SLURM array, ONE raw file per task (--import-no-join; import_args): never two
            readers on one timsTOF .d ("database is locked"), and never --import-process-count,
            which hangs under Wine.
  3 join    one --batch-commands call joins the per-file caches and exports the precursor and
            transition reports, the targets' chromatograms and the shareable .sky.zip (each
            SkylineCmd call costs 8-10 s of Wine start-up, so they are batched); then the pictures
            are drawn (matplotlib) and the record is written.
  A step that fails writes the record with the reason (status "skipped"), so the report says it;
  nothing is dropped silently (architectural rule 4).

Stdlib only at import; pyarrow (reading report.parquet) and matplotlib (drawing) are imported
inside the jobs that need them.
"""
import argparse
import csv
import datetime
import hashlib
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

from resolve_defaults import SKYLINE_PIN, family as _instrument_family  # noqa: E402  the pins
from protein_ids import header_accession, group_accessions              # noqa: E402
from fetch_fasta import CONTAMINANT_TAGS                                 # noqa: E402
from probe_window import MS2_ACC_RE, MS1_ACC_RE                         # noqa: E402  DIA-NN's log

SCHEMA = "skyline_peak_check/1"
RECORD = "skyline_provenance.json"
WORK = "skyline_work"               # scratch_files.SCRATCH_DIRS: never zipped or delivered
SECTION_TITLE = "Peak check (Skyline)"
FOLLOWUP_DIR = "skyline_followup"   # output/<this>: the staff-only PRM/SRM lists

# ---- the rules, each in one place ----------------------------------------------------------------
TOP_N_DEFAULT = 10                  # proteins checked (--top)
APEX_DIFF_S = 6.0                   # Skyline's apex this far from DIA-NN's: a different peak
PRECURSORS_DRAWN = 3                # per protein, the strongest by Skyline's fragment area
MAX_RUNS_DRAWN = 12                 # columns in a peak picture (every run is in the tables)
FRAGMENTS_DRAWN = 8
# Chromatogram extraction around DIA-NN's identification times (--full-scan-rt-filter=ms2_ids).
# The feasibility trial's two settings, both checked against DIA-NN: +-1 min on an 11-min timsTOF
# gradient and +-2 min on a 35-min Lumos one. Runs longer than LONG_GRADIENT_MIN get the wider one.
LONG_GRADIENT_MIN = 20.0
RT_WINDOW_MIN = {"short": 1.0, "long": 2.0}
IMS_LIBRARY_RES = 50                # ion-mobility filter from the library's 1/K0 (timsTOF)
# Product ions: the trial's settings, the ones its agreement numbers were measured with.
PRODUCT_ION_ARGS = ["--tran-product-ion-charges=1,2", "--tran-product-ion-types=y, b, p",
                    "--library-pick-product-ions=all", "--library-product-ions=6",
                    "--library-min-product-ions=3", "--library-match-tolerance=0.05mz"]
# Job sizes, from the trial's sacct: build 23-47 s / <= 2.2 GB, one import 0.6-4 min / <= 3.9 GB
# (a whole proteome), join + exports ~1-2.5 min / <= 3.8 GB. Generous on time, modest on CPUs:
# MAX_PARALLEL_DEFAULT tasks x JOB_CPUS = 32 of a user's 64 on HIVE's high.
JOB_CPUS = 4
JOB_MEM = "16G"
JOB_MINUTES = {"build": 60, "import": 60, "join": 90, "prm": 30}
MAX_PARALLEL_DEFAULT = 8
# The follow-up assay: the 3 best peptides, 6 best fragments and one charge per protein, judged
# from the measured areas (Skyline's own refinement).
PRM_REFINE_ARGS = ["--refine-max-peptide-peak-rank=3", "--refine-max-transition-peak-rank=6",
                   "--refine-max-precursor-only"]
# key -> (list kind, Skyline's instrument name, instrument family). Only what SkylineCmd 26.1's
# --help offers: isolation lists for Bruker timsTOF, Thermo Fusion (Fusion, Lumos, Eclipse),
# Thermo Q Exactive and Thermo Stellar; transition lists for the TSQ Quantiva/Altis family.
# It offers NO isolation list for the Exploris or the Astral -- only a native method export, which
# needs a method template from the instrument PC.
PRM_INSTRUMENTS = {
    "bruker_timstof": ("isolation", "Bruker timsTOF", "timstof"),
    "thermo_fusion": ("isolation", "Thermo Fusion", "orbitrap"),
    "thermo_qe": ("isolation", "Thermo Q Exactive", "orbitrap"),
    "thermo_stellar": ("isolation", "Thermo Stellar", "ion_trap"),
    "thermo_srm": ("transition", "Thermo Quantiva", "triple_quad"),
}
PRM_NOT_OFFERED = {"thermo_exploris": "Exploris", "thermo_astral": "Astral"}
SRM_CE_PREDICTOR = "Thermo TSQ Quantiva"

# Skyline 26.1's predefined enzymes, name -> (cleaves, not before, N-terminal), from its source
# (pwiz_tools/Skyline/Properties/Settings.cs at tag Skyline-26.1.0.266), in the order a DIA-NN
# --cut rule is matched against them (the first that cuts the same way wins: TrypsinK = LysC).
SKYLINE_ENZYMES = [
    ("Trypsin", "KR", "P", False), ("Trypsin/P", "KR", "", False), ("LysC", "K", "P", False),
    ("LysC/P", "K", "", False), ("ArgC", "R", "P", False), ("Clostripain", "R", "", False),
    ("Chymotrypsin", "FWYL", "P", False), ("GluC", "DE", "P", False),
    ("GluC bicarb", "E", "P", False), ("Trypsin-GluC", "DEKR", "P", False),
    ("Trypsin-CNBr", "KRM", "P", False), ("CNBr", "M", "P", False),
    ("Elastase", "GVLIA", "P", False), ("Formic Acid", "D", "P", False),
    ("Iodosobenzoate", "W", "", False), ("PepsinA", "FL", "", False),
    ("Protein endopeptidase", "P", "", False), ("Staph protease", "E", "", False),
    ("AspN", "D", "", True), ("LysN", "K", "", True), ("LysN promisc", "KASR", "", True),
]
# The DIA-NN flags the document must copy; without one the step does not guess its value.
SEARCH_FLAGS = ("cut", "missed-cleavages", "min-pep-len", "max-pep-len", "min-pr-charge",
                "max-pr-charge", "min-fr-mz", "max-fr-mz")

# Report definitions (Skyline .skyr), written into the work folder: the trial's columns.
REPORT_PRECURSORS = "Skill Peak Check Precursors"
REPORT_TRANSITIONS = "Skill Peak Check Transitions"
_PREC_COLS = ["Peptide.Protein.Name", "Peptide.Protein.Accession", "Peptide.Protein.Gene",
              "Peptide.Sequence", "Peptide.ModifiedSequence", "Charge", "Mz", "IsDecoy",
              "LibraryName", "LibraryScore1"] + ["Results!*.Value." + c for c in (
                  "PeptideResult.ResultFile.Replicate.Name", "PeptideResult.ResultFile.FileName",
                  "BestRetentionTime", "MinStartTime", "MaxEndTime", "MaxFwhm", "TotalArea",
                  "TotalAreaMs1", "TotalAreaFragment", "TotalBackground", "MaxHeight",
                  "LibraryDotProduct", "IsotopeDotProduct", "AverageMassErrorPPM",
                  "CountTruncated", "PrecursorPeakFoundRatio", "Identified", "IonMobilityMS1",
                  "IonMobilityFragment", "IonMobilityWindow", "CollisionalCrossSection",
                  "DetectionQValue", "DetectionZScore")]
_TRAN_COLS = ["Precursor.Peptide.Protein.Name", "Precursor.Peptide.ModifiedSequence",
              "Precursor.Charge", "Precursor.Mz", "FragmentIon", "ProductCharge", "ProductMz",
              "IsotopeDistIndex", "LibraryRank", "LibraryIntensity", "Quantitative"] + [
                  "Results!*.Value." + c for c in (
                      "PrecursorResult.PeptideResult.ResultFile.FileName", "RetentionTime",
                      "StartTime", "EndTime", "Area", "Background", "Height", "Fwhm",
                      "MassErrorPPM", "PeakRank", "Truncated", "Coeluting")]


def _skyr(name, rowsource, cols):
    lines = ['<?xml version="1.0"?>', "<views>",
             f'  <view name="{name}" rowsource="pwiz.Skyline.Model.Databinding.Entities.'
             f'{rowsource}" sublist="Results!*">']
    lines += [f'    <column name="{c}" />' for c in cols]
    return "\n".join(lines + ["  </view>", "</views>", ""])


SKYR = {"precursor.skyr": _skyr(REPORT_PRECURSORS, "Precursor", _PREC_COLS),
        "transition.skyr": _skyr(REPORT_TRANSITIONS, "Transition", _TRAN_COLS)}


class Skip(Exception):
    """The step cannot run here. The message is the reason the report, the Methods and
    MANIFEST.txt give, so it is written for a reader, not a log."""


# ---- small helpers ---------------------------------------------------------------------------------
def now():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def win(path):
    """The Windows path SkylineCmd sees under Wine, which maps Z: to /."""
    return "Z:" + os.path.abspath(path).replace("/", "\\")


def load_json(path):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def write_json(path, obj):
    """Whole or not at all: <path>.part, then rename (scratch_files: *.part is scratch)."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, path)


def write_text(path, text):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    tmp = path + ".part"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.replace(tmp, path)


def sha256_file(path, chunk=1 << 22):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def fmt_num(x):
    return f"{x:g}" if isinstance(x, (int, float)) else str(x)


def skill_version():
    try:
        import skill_version as sv
        return sv.skill_version()
    except Exception:                    # recorded as unknown, never guessed
        return None


def _float(v):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return f if math.isfinite(f) else None


def default_tag():
    from make_analysis_html import DEFAULT_TAG     # the one wording (architectural rule 2)
    return DEFAULT_TAG


# ---- the DIA-NN search, read -------------------------------------------------------------------------
def diann_command(log_text):
    """The DIA-NN command line its report.log.txt starts with (the first line that names the
    program and has flags), or None."""
    for line in (log_text or "").splitlines()[:20]:
        if re.search(r"diann", line, re.I) and re.search(r"(?:^|\s)--[a-z]", line):
            return line.strip()
    return None


def diann_version(log_text):
    m = re.match(r"\s*DIA-NN\s+(\d+(?:\.\d+)+)", log_text or "")
    return m.group(1) if m else None


def command_flags(cmd):
    """{flag: [value, ...]} of a DIA-NN command line; a switch has [""]. DIA-NN's log prints the
    command space-joined, so a value is everything up to the next --flag (a path with a space
    stays whole)."""
    toks, out, i = (cmd or "").split(), {}, 0
    while i < len(toks):
        if toks[i].startswith("--"):
            j = i + 1
            while j < len(toks) and not toks[j].startswith("--"):
                j += 1
            out.setdefault(toks[i][2:], []).append(" ".join(toks[i + 1:j]))
            i = j
        else:
            i += 1
    return out


_FIXED_RE = re.compile(r"Mass\s+accuracy\s+will\s+be\s+fixed\s+to\s+([0-9.eE+-]+)\s*\(MS2\)\s*"
                       r"and\s+([0-9.eE+-]+)\s*\(MS1\)", re.I)


def mass_accuracy(log_text, search_prov=None, ms1=None, ms2=None, source=None):
    """The extraction tolerance (+- ppm) for Skyline: DIA-NN's own mass accuracy for this search,
    never a fixed number (the trial: 10 ppm gave gappy, broken peaks on a run DIA-NN had measured
    at 24 ppm; 24-25 fixed it). Both programs use it as +- ppm around the m/z: DIA-NN's Scan::level
    (src/diann.cpp at tag 1.8) and Skyline's centroided window (TransitionFullScan).
    In order: the user's value; the value the search FIXED (its log: "Mass accuracy will be fixed
    to ... (MS2) and ... (MS1)", or --mass-acc/--mass-acc-ms1 in its command); DIA-NN's measured
    "Optimised mass accuracy" (MS2; DIA-NN logs only a recommendation for MS1, so MS1 uses the
    same value, as the trial did); search_provenance.json `mass_acc`. -> {ms1_ppm, ms2_ppm,
    source}; Skip when none of them gives one."""
    if ms2 is not None or ms1 is not None:
        if ms2 is None or ms1 is None:
            raise Skip("only one of --ppm-ms1 / --ppm-ms2 was given: give both")
        out = {"ms1_ppm": float(ms1), "ms2_ppm": float(ms2),
               "source": source or "given by the user (--ppm-ms1 / --ppm-ms2)"}
    else:
        out = None
        m = _FIXED_RE.search(log_text or "")
        if m:
            out = {"ms1_ppm": round(float(m.group(2)) * 1e6, 4),
                   "ms2_ppm": round(float(m.group(1)) * 1e6, 4),
                   "source": "fixed in the DIA-NN search (report.log.txt: \"Mass accuracy will be "
                             "fixed to ...\", its --mass-acc / --mass-acc-ms1)"}
        if out is None:
            fl = command_flags(diann_command(log_text) or "")
            a2, a1 = _float((fl.get("mass-acc") or [None])[-1]), \
                _float((fl.get("mass-acc-ms1") or [None])[-1])
            if a2:
                out = {"ms1_ppm": a1 or a2, "ms2_ppm": a2,
                       "source": "the DIA-NN command's --mass-acc" + (" / --mass-acc-ms1" if a1
                                                                       else " (no --mass-acc-ms1: "
                                                                       "MS1 uses the same value)")}
        if out is None:
            opt = [float(x) for x in MS2_ACC_RE.findall(log_text or "")]
            if opt:
                rec = sorted({float(x) for x in MS1_ACC_RE.findall(log_text or "")})
                out = {"ms1_ppm": max(opt), "ms2_ppm": max(opt),
                       "source": f"measured by DIA-NN for this search (report.log.txt "
                                 f"\"Optimised mass accuracy: {fmt_num(max(opt))} ppm\"); DIA-NN "
                                 f"logs only a recommendation for MS1"
                                 + (f" ({', '.join(fmt_num(x) for x in rec)} ppm)" if rec else "")
                                 + ", so MS1 uses the same value"}
        if out is None and isinstance((search_prov or {}).get("mass_acc"), dict):
            ma = search_prov["mass_acc"]
            a2, a1 = _float(ma.get("ms2")), _float(ma.get("ms1"))
            if a2:
                mode = ma.get("mode") or "mode not recorded"
                why = ma.get("reason") or "no reason recorded"
                out = {"ms1_ppm": a1 or a2, "ms2_ppm": a2,
                       "source": f"search_provenance.json mass_acc ({mode}: {why})"}
    if out is None:
        raise Skip("DIA-NN's mass accuracy for this search is not in its report.log.txt or "
                   "search_provenance.json, and Skyline's tolerance must match it (a guessed "
                   "10 ppm gave broken peaks in the feasibility trial). Give the search's value "
                   "with --ppm-ms1 / --ppm-ms2 and where it came from with --ppm-source")
    for k in ("ms1_ppm", "ms2_ppm"):
        if not (1.0 <= out[k] <= 1000.0):
            raise Skip(f"mass accuracy {fmt_num(out[k])} ppm ({k[:3].upper()}) is outside "
                       f"Skyline's 1-1000 ppm range")
    return out


def enzyme_for(cut):
    """Skyline's predefined enzyme that cuts as DIA-NN's --cut rule does ("K*,R*" -> Trypsin/P,
    "K*,R*,!*P" -> Trypsin), or Skip: a rule with no Skyline twin is never approximated."""
    cterm, nterm, nobefore = set(), set(), set()
    for tok in (t.strip() for t in (cut or "").split(",") if t.strip()):
        m = re.fullmatch(r"([A-Z])\*", tok) or None
        if m:
            cterm.add(m.group(1))
            continue
        m = re.fullmatch(r"!\*([A-Z])", tok)
        if m:
            nobefore.add(m.group(1))
            continue
        m = re.fullmatch(r"\*([A-Z])", tok)
        if m:
            nterm.add(m.group(1))
            continue
        raise Skip(f"the search's cleavage rule --cut {cut} has a part ({tok}) the skill cannot "
                   f"match to a Skyline enzyme")
    if not (cterm or nterm) or (cterm and nterm) or (nterm and nobefore):
        raise Skip(f"the search's cleavage rule --cut {cut} has no Skyline enzyme to match")
    for name, cleave, restrict, n_terminal in SKYLINE_ENZYMES:
        if n_terminal == bool(nterm) and set(cleave) == (nterm or cterm) \
                and set(restrict) == nobefore:
            return name
    raise Skip(f"the search's cleavage rule --cut {cut} has no Skyline enzyme to match")


def search_settings(flags):
    """What the Skyline document copies from the DIA-NN command (the trial's way): the enzyme,
    missed cleavages, peptide lengths, precursor charges and fragment m/z range. A flag the
    command does not state is not guessed: Skip."""
    missing = [f"--{f}" for f in SEARCH_FLAGS if not (flags.get(f) or [""])[-1].strip()]
    if missing:
        raise Skip(f"the DIA-NN command in report.log.txt does not state {', '.join(missing)}, "
                   f"which the Skyline document must copy")
    if "unimod4" not in flags:
        raise Skip("the search did not fix carbamidomethyl cysteine (--unimod4); the skill's "
                   "Skyline document assumes it (Skyline's default), so the peak check is not "
                   "run for this search")

    def num(f, kind=int):
        v = flags[f][-1].strip()
        try:
            return kind(float(v)) if kind is int else kind(v)
        except ValueError:
            raise Skip(f"--{f} {v} in the DIA-NN command is not a number")
    zmin, zmax = num("min-pr-charge"), num("max-pr-charge")
    if not (1 <= zmin <= zmax <= 10):
        raise Skip(f"precursor charges {zmin}-{zmax} in the DIA-NN command are not usable")
    out = {"cut": flags["cut"][-1].strip(), "enzyme": enzyme_for(flags["cut"][-1].strip()),
           "missed_cleavages": num("missed-cleavages"), "min_length": num("min-pep-len"),
           "max_length": num("max-pep-len"), "charges": list(range(zmin, zmax + 1)),
           "fragment_mz": [num("min-fr-mz", float), num("max-fr-mz", float)],
           "var_mods": [v for v in flags.get("var-mod", []) if v],
           "max_var_mods": num("var-mods") if flags.get("var-mods") else None,
           "fixed_mods": ["Carbamidomethyl (C)"],
           "not_carried_over": []}
    if "met-excision" in flags:
        out["not_carried_over"].append(
            "DIA-NN's N-terminal methionine excision (--met-excision): SkylineCmd has no option "
            "for it, so peptides that exist only after the excision are not targets")
    return out


def instrument_family(manifest, raw_files):
    """-> (family, source): "timstof" (Bruker .d, TOF analysers with ion mobility) or "orbitrap"
    (Thermo .raw: Orbitrap and Astral, centroided). Skip for anything else -- the skill has run
    Skyline on those two only."""
    man = manifest or {}
    res = man.get("resolution") if isinstance(man.get("resolution"), dict) else {}
    if res.get("ms2_analyzer") in ("ITMS", "mixed"):
        raise Skip("the MS2 spectra were recorded in the ion trap, which the skill's Skyline "
                   "settings (high-resolution, centroided) do not cover")
    exts = {os.path.splitext(f.rstrip("/"))[1].lower() for f in raw_files}
    cls = man.get("instrument_class")
    if cls:
        fam, src = _instrument_family(cls), f"workflow.manifest.json instrument_class {cls}"
    elif exts == {".d"}:
        fam, src = "timstof", "the raw files are Bruker .d folders"
    elif exts == {".raw"}:
        fam, src = "orbitrap", "the raw files are Thermo .raw files"
    else:
        raise Skip(f"the instrument is not recorded and the raw files ({', '.join(sorted(exts))}) "
                   f"do not say it; the skill runs Skyline on Bruker timsTOF .d and Thermo .raw")
    want = {"timstof": ".d", "orbitrap": ".raw"}.get(fam)
    if want is None:
        raise Skip(f"the instrument ({cls}) is neither a Bruker timsTOF nor a Thermo Orbitrap / "
                   f"Astral; the skill has run Skyline on those two only")
    if exts != {want}:
        raise Skip(f"the raw files ({', '.join(sorted(exts))}) are not {want} files, which a "
                   f"{fam} search reads; the skill imports vendor files only")
    return fam, src


def tof_resolving_power(ppm):
    """Skyline's TOF extraction window is m/z / (R / 2) wide (TransitionFullScan.GetFilterWindow
    with ResPerFilter 2), i.e. +- m/z / R: +- ppm needs R = 1e6 / ppm."""
    return int(round(1e6 / ppm))


def analyzer_args(fam, ms1_ppm, ms2_ppm):
    """The mass analysers, ALWAYS with their resolution: setting the acquisition method alone left
    Skyline's default product analyser, "qit" at 0.7 (the trial)."""
    if fam == "timstof":
        return ["--full-scan-precursor-analyzer=tof",
                f"--full-scan-precursor-res={tof_resolving_power(ms1_ppm)}",
                "--full-scan-product-analyzer=tof",
                f"--full-scan-product-res={tof_resolving_power(ms2_ppm)}",
                f"--ims-library-res={IMS_LIBRARY_RES}"]
    return ["--full-scan-precursor-analyzer=centroided", f"--full-scan-precursor-res={fmt_num(ms1_ppm)}",
            "--full-scan-product-analyzer=centroided", f"--full-scan-product-res={fmt_num(ms2_ppm)}"]


def rt_window(gradient_min):
    """+- minutes of chromatogram around DIA-NN's identification times -> (minutes, source)."""
    key = "long" if (gradient_min or 0) > LONG_GRADIENT_MIN else "short"
    return RT_WINDOW_MIN[key], (
        f"the skill's setting for runs {'longer' if key == 'long' else 'no longer'} than "
        f"{fmt_num(LONG_GRADIENT_MIN)} min (the last DIA-NN retention time here: "
        f"{fmt_num(round(gradient_min or 0, 1))} min); the feasibility trial's two checked settings")


# ---- the Skyline command lines -----------------------------------------------------------------------
def build_args(plan, rt_minutes):
    """ONE SkylineCmd call: new document, settings, DIA-NN's library, the targets' FASTA, save."""
    s, ma = plan["settings"], plan["mass_accuracy"]
    args = [f"--new={win(plan['doc'])}", "--overwrite",
            "--full-scan-acquisition-method=DIA",
            f"--full-scan-isolation-scheme={win(plan['raw_files'][0])}",
            "--full-scan-precursor-isotopes=Count", "--full-scan-precursor-threshold=3"]
    args += analyzer_args(plan["family"], ma["ms1_ppm"], ma["ms2_ppm"])
    args += ["--full-scan-rt-filter=ms2_ids", f"--full-scan-rt-filter-tolerance={fmt_num(rt_minutes)}",
             f"--pep-digest-enzyme={s['enzyme']}", f"--pep-max-missed-cleavages={s['missed_cleavages']}",
             f"--pep-min-length={s['min_length']}", f"--pep-max-length={s['max_length']}",
             f"--tran-precursor-ion-charges={','.join(str(z) for z in s['charges'])}"]
    args += PRODUCT_ION_ARGS
    args += [f"--instrument-min-mz={fmt_num(s['fragment_mz'][0])}",
             f"--instrument-max-mz={fmt_num(s['fragment_mz'][1])}"]
    if s.get("max_var_mods"):
        args.append(f"--pep-max-variable-mods={s['max_var_mods']}")
    args += [f"--import-search-file={win(plan['speclib_link'])}", "--import-search-add-mods",
             f"--import-fasta={win(plan['targets_fasta'])}", "--save"]
    return args


def import_args(plan, i):
    """Array task i (1-based): ONE raw file, its own cache, no join and no save."""
    return [f"--in={win(plan['doc'])}", f"--import-file={win(plan['raw_files'][i - 1])}",
            "--import-no-join"]


def batch_line(args):
    """One --batch-commands line: a path, or any value with a space, in double quotes (the
    trial's form); a plain word as it is."""
    out = []
    for a in args:
        if '"' in a or "\n" in a:
            raise Skip(f"a path holds a double quote or a line break, which a SkylineCmd batch "
                       f"file cannot carry: {a}")
        k, eq, v = a.partition("=")
        out.append((f'{k}="{v}"' if re.search(r"[\s\\:,]", v) else f"{k}={v}") if eq else k)
    return " ".join(out)


def join_batch(plan, imported):
    """The join job's ONE SkylineCmd call (--batch-commands): join the per-file caches of the
    files that imported, then the two reports, the targets' chromatograms and the .sky.zip."""
    doc, w, out = plan["doc"], plan["work"], plan["out"]
    lines = [[f"--in={win(doc)}"] + [f"--import-file={win(f)}" for f in imported] + ["--save"]]
    for skyr, name, dest in (("precursor.skyr", REPORT_PRECURSORS,
                              os.path.join(out, "skyline_precursors.tsv")),
                             ("transition.skyr", REPORT_TRANSITIONS,
                              os.path.join(w, "skyline_transitions.tsv"))):
        lines.append([f"--in={win(doc)}", f"--report-add={win(os.path.join(w, skyr))}",
                      "--report-conflict-resolution=overwrite", f"--report-name={name}",
                      f"--report-file={win(dest)}", "--report-format=tsv", "--report-invariant"])
    lines.append([f"--in={win(doc)}",
                  f"--chromatogram-file={win(os.path.join(w, 'skyline_chromatograms.tsv'))}",
                  "--chromatogram-precursors", "--chromatogram-products"])
    lines.append([f"--in={win(doc)}", f"--share-zip={win(plan['share_zip'])}",
                  "--share-type=complete"])
    return "\n".join(batch_line(x) for x in lines) + "\n"


def apptainer_argv(plan, rw=(), tmp=None):
    """apptainer exec with the binds the step needs and no more: the output folder (and `rw`)
    read-write, the search folder, the FASTA's and every raw-data folder read-only (nested binds
    go on top, so a raw folder under the output folder is still read-only). `tmp`: this call's own
    temporary folder, as Windows' TMP and TEMP -- the image's TMP is the node's shared /tmp, and
    two imports started together on one node lost a spill file to each other there (the skill's
    test on HIVE, 2026-10-20: "Unable to find the specified file" in ChromGroups)."""
    argv = ["apptainer", "exec", "--cleanenv", "--env", "WINEDEBUG=-all"]
    if tmp:
        argv += ["--env", f"TMP={tmp}", "--env", f"TEMP={tmp}"]
    for d in [plan["out"]] + [x for x in rw if x != plan["out"]] + ([tmp] if tmp else []):
        argv += ["--bind", f"{d}:{d}"]
    for d in plan["ro_binds"]:
        argv += ["--bind", f"{d}:{d}:ro"]
    return argv + [plan["sif"]]


def ro_binds(search_dir, raw_files, fasta=None):
    """The real folders to bind read-only: the search's, the raw files', the FASTA's -- a folder
    inside another is covered by it."""
    dirs = {os.path.realpath(search_dir)}
    for f in raw_files:
        dirs.add(os.path.realpath(os.path.dirname(os.path.abspath(f.rstrip("/")))))
    if fasta:
        dirs.add(os.path.realpath(os.path.dirname(os.path.abspath(fasta))))
    keep = []
    for d in sorted(dirs, key=lambda x: (len(x), x)):
        if not any(d == k or d.startswith(k.rstrip("/") + "/") for k in keep):
            keep.append(d)
    keep.sort()
    for d in keep:
        if ":" in d or "," in d:
            raise Skip(f"the folder {d} holds ':' or ',', which an apptainer --bind cannot carry")
    return keep


# ---- the jobs -------------------------------------------------------------------------------------------
def sbatch_header(name, log_dir, minutes, queue, cpus=JOB_CPUS, mem=JOB_MEM, array=None):
    """ucd-hpc's canonical job: --account/--partition/--qos from the skill's queue rule, --time,
    --mem with its unit, --ntasks=1, --cpus-per-task, one log per job (%j; %A_%a for the array),
    modules.sh before `module`. No --exclusive."""
    part, acct, qos = queue
    log = os.path.join(log_dir, f"{name}_%A_%a.log" if array else f"{name}_%j.log")
    if re.search(r'["\n]', log):
        raise Skip(f"the output folder holds a quote or a line break, which #SBATCH --output "
                   f"cannot carry: {log_dir}")
    h = ["#!/bin/bash -l", f"#SBATCH --job-name={name}"]
    h += [f"#SBATCH --account={acct}"] if acct else []
    h += [f"#SBATCH --partition={part}"] if part else []
    h += [f"#SBATCH --qos={qos}"] if qos else []
    h += [f"#SBATCH --array={array}"] if array else []
    h += [f"#SBATCH --time={minutes}", "#SBATCH --ntasks=1", f"#SBATCH --cpus-per-task={cpus}",
          f"#SBATCH --mem={mem}", f'#SBATCH --output="{log}"', f'#SBATCH --error="{log}"']
    h += ["#SBATCH --requeue"] if part == "low" else []
    return h


def job_script(plan, step, queue, array=None):
    """A thin wrapper: modules, apptainer, then this script's `step` subcommand."""
    name = {"build": "sky_build", "import": "sky_import", "join": "sky_join", "prm": "sky_prm"}[step]
    h = sbatch_header(name, os.path.join(plan["work"], "logs"), JOB_MINUTES[step], queue,
                      array=array)
    q = shlex.quote
    sub = f"{q(plan['python'])} {q(plan['script'])} {step} --work {q(plan['work'])}"
    if step == "import":
        sub += ' --task "$SLURM_ARRAY_TASK_ID"'
    body = ["# written by skyline_peaks.py plan: the peak check (references/skyline.md)",
            "set -uo pipefail",
            "if [ -r /etc/profile.d/modules.sh ]; then",
            "  source /etc/profile.d/modules.sh",
            "  module load apptainer",
            "fi",
            f"export MPLCONFIGDIR={q(os.path.join(plan['work'], '.mpl'))}",
            f"exec {sub}"]
    return "\n".join(h + body) + "\n"


def submit_script(plan):
    q = shlex.quote
    w = plan["work"]
    return "\n".join([
        "#!/bin/bash",
        "# written by skyline_peaks.py plan: build -> import (one raw file per array task) -> join",
        "set -euo pipefail",
        f"cd {q(w)}",
        f"j1=$(sbatch --parsable {q(os.path.join(w, 'sky_1_build.sbatch'))} | cut -d';' -f1)",
        # a failed build cancels the imports (--kill-on-invalid-dep); the join runs after the
        # imports whatever happened, and writes the record of what did
        f"j2=$(sbatch --parsable --dependency=afterok:$j1 --kill-on-invalid-dep=yes "
        f"{q(os.path.join(w, 'sky_2_import.sbatch'))} | cut -d';' -f1)",
        f"j3=$(sbatch --parsable --dependency=afterany:$j2 "
        f"{q(os.path.join(w, 'sky_3_join.sbatch'))} | cut -d';' -f1)",
        f"printf 'build\\t%s\\nimport\\t%s\\njoin\\t%s\\n' \"$j1\" \"$j2\" \"$j3\" > {q(os.path.join(w, 'jobs.txt'))}",
        'echo "submitted: build $j1, import $j2 (array), join $j3"',
        ""])


# ---- choosing the proteins: the ONE rule ---------------------------------------------------------------
SELECTION_RULE = (
    "the {n} proteins with the strongest evidence of change: the significant proteins (adjusted p "
    "below the DE's own cutoff) ranked by their smallest adjusted p in any comparison, ties broken "
    "by the raw p-value; contaminant entries left out. When no protein is significant, the {n} "
    "most abundant proteins measured in every run instead (the mean of the DE's log2 matrix; "
    "without a DE, the mean of DIA-NN's PG.MaxLFQ over the runs, protein groups with at least 3 "
    "precursors).")


def accessions(group):
    """The members of a Protein.Group (protein_ids.group_accessions), as a list."""
    return [x for x in group_accessions(group).split(";") if x]


def is_contaminant(group):
    return any(acc.startswith(CONTAMINANT_TAGS) for acc in accessions(group))


def first_gene(genes, fallback):
    g = (genes or "").split(";")[0].strip()
    return g or fallback


def select_de_hits(de_dir, n):
    """-> (targets, adjp, adjp_source): the significant proteins, best first (make_analysis_html's
    Tables: the same rows, cutoff and order as the report's top-protein tables)."""
    import make_analysis_html as mah
    prov = mah.load_record(os.path.join(de_dir, "de_provenance.json"))
    adjp, src = mah.significance_rule(de_dir, 0.05)
    tables = mah.Tables(de_dir, prov, adjp, src)
    best = {}
    for c in tables.contrasts:
        for r in tables.rows(c):
            pg = r.get("Protein.Group") or ""
            if r["_p"] is None or r["_p"] >= adjp or not pg or is_contaminant(pg):
                continue
            key = (r["_p"], r["_raw_p"] if r["_raw_p"] is not None else float("inf"))
            if pg not in best or key < best[pg][0]:
                best[pg] = (key, {"protein_group": pg, "gene": first_gene(r.get("Genes"), pg),
                                  "contrast": tables.display(c), "adj_p": r["_p"],
                                  "p": r["_raw_p"], "logFC": r["_lfc"],
                                  "evidence": r.get("Evidence") or ""})
    ranked = sorted(best.values(), key=lambda kv: (kv[0], kv[1]["protein_group"]))
    return [t for _, t in ranked[:n]], adjp, src


def select_abundant_de(de_dir, n):
    """The most abundant proteins measured in every run, from the DE's own matrices."""
    em = os.path.join(de_dir, "Expression_Matrix.csv")
    dm = os.path.join(de_dir, "Detection_Matrix.csv")
    if not os.path.isfile(em):
        raise Skip("no protein is significant and the DE wrote no Expression_Matrix.csv to "
                   "choose the most abundant proteins from")
    measured = None
    if os.path.isfile(dm):
        measured = {}
        with open(dm, newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                vals = [_float(v) for k, v in r.items() if k != "Protein.Group"]
                measured[r.get("Protein.Group")] = bool(vals) and all(v and v > 0 for v in vals)
    rows = []
    with open(em, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        samples = [c for c in rd.fieldnames or [] if c not in ("Protein.Group", "Genes",
                                                                  "Protein.Names")]
        for r in rd:
            pg = r.get("Protein.Group") or ""
            if not pg or is_contaminant(pg) or (measured is not None and not measured.get(pg)):
                continue
            vals = [_float(r.get(s)) for s in samples]
            vals = [v for v in vals if v is not None]
            if vals:
                rows.append((-sum(vals) / len(vals), pg, first_gene(r.get("Genes"), pg)))
    rows.sort()
    return [{"protein_group": pg, "gene": g, "mean_log2": round(-m, 3)} for m, pg, g in rows[:n]]


def _q_filtered(report, cols, q_cutoff):
    """report.parquet's `cols` (those it has), only the rows that pass the DE's q-value filter
    (diann_q_columns: every q column the report has, each at its cutoff) -- filtered in Arrow, so
    a report of millions of rows is never turned into Python objects. -> (table, columns read)."""
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    import diann_q_columns as dq
    have = pq.read_schema(report).names
    qcols = dq.fdr_columns(have)
    cols = [c for c in cols if c in have]
    t = pq.read_table(report, columns=list(dict.fromkeys(cols + qcols)))
    mask = None
    for c in qcols:
        m = pc.less_equal(t[c], dq.cutoff_for(c, q_cutoff))
        mask = m if mask is None else pc.and_(mask, m)
    return (t.filter(mask) if mask is not None else t), cols


def select_abundant_report(report, n, q_cutoff):
    """Without a DE: DIA-NN's PG.MaxLFQ mean over the runs, protein groups of one accession
    measured in every run with at least 3 precursors (the trial's choice). Needs pyarrow."""
    import pyarrow.compute as pc
    t, cols = _q_filtered(report, ["Run", "Protein.Group", "Genes", "Precursor.Id", "PG.MaxLFQ"],
                          q_cutoff)
    if len(cols) < 5:
        raise Skip("report.parquet lacks Run / Protein.Group / Genes / Precursor.Id / PG.MaxLFQ")
    n_runs = len(pc.unique(t["Run"]))
    per_run = t.group_by(["Protein.Group", "Run"]).aggregate([("PG.MaxLFQ", "max")])
    per_pg = per_run.group_by("Protein.Group").aggregate([("PG.MaxLFQ_max", "mean"),
                                                          ("Run", "count")])
    prec = {r["Protein.Group"]: r["Precursor.Id_count_distinct"] for r in
            t.group_by("Protein.Group").aggregate([("Precursor.Id", "count_distinct")]).to_pylist()}
    genes = {r["Protein.Group"]: r["Genes_min"] for r in
             t.group_by("Protein.Group").aggregate([("Genes", "min")]).to_pylist()}
    rows = []
    for r in per_pg.to_pylist():
        pg = r["Protein.Group"]
        if (not pg or ";" in pg or is_contaminant(pg) or r["Run_count"] != n_runs
                or prec.get(pg, 0) < 3 or r["PG.MaxLFQ_max_mean"] is None):
            continue
        rows.append((-r["PG.MaxLFQ_max_mean"], pg, first_gene(genes.get(pg), pg)))
    rows.sort()
    return [{"protein_group": pg, "gene": g, "mean_maxlfq": -m} for m, pg, g in rows[:n]]


def select_targets(plan):
    """THE rule (SELECTION_RULE) -> {"mode", "targets", "rule", ...}."""
    n = plan["top"]["value"]
    de_dir = plan.get("de_dir")
    if de_dir and os.path.isfile(os.path.join(de_dir, "de_provenance.json")):
        hits, adjp, src = select_de_hits(de_dir, n)
        if hits:
            return {"mode": "de_hits", "targets": hits, "adjp": adjp, "adjp_source": src}
        return {"mode": "abundance_no_hits", "targets": select_abundant_de(de_dir, n),
                "adjp": adjp, "adjp_source": src}
    return {"mode": "abundance_no_de",
            "targets": select_abundant_report(plan["report"], n, plan["q_cutoff"])}


def fasta_subset(fasta, accessions, dest):
    """The targets' entries of the searched FASTA -> {accession: header} written."""
    want, found, keep, out = set(accessions), {}, False, []
    with open(fasta, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                acc = header_accession(line[1:].split()[0]) if line[1:].strip() else ""
                keep = acc in want and acc not in found
                if keep:
                    found[acc] = line[1:].strip()
            if keep:
                out.append(line if line.endswith("\n") else line + "\n")
    write_text(dest, "".join(out))
    return found


DIANN_COLS = ("Run", "Protein.Group", "Genes", "Precursor.Id", "Modified.Sequence",
              "Stripped.Sequence", "Precursor.Charge", "RT", "RT.Start", "RT.Stop", "IM",
              "Precursor.Quantity", "Ms1.Area", "Q.Value")


def diann_rows(report, groups, q_cutoff, dest):
    """DIA-NN's rows for the target protein groups (the DE's q-value filter: diann_q_columns),
    written as a TSV -> (rows written, the last retention time of the whole report, minutes)."""
    import pyarrow as pa
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    last = pc.max(pq.read_table(report, columns=["RT"])["RT"]).as_py() \
        if "RT" in pq.read_schema(report).names else None
    t, cols = _q_filtered(report, list(DIANN_COLS), q_cutoff)
    t = t.filter(pc.is_in(t["Protein.Group"], value_set=pa.array(sorted(set(groups)))))
    n = 0
    with open(dest + ".part", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(cols)
        for r in t.select(cols).to_pylist():
            w.writerow([r[c] for c in cols])
            n += 1
    os.replace(dest + ".part", dest)
    return n, last


# ---- reading Skyline's reports ------------------------------------------------------------------------
_MOD = re.compile(r"(\[[^\]]*\]|\([^)]*\))")


def precursor_key(modified_seq, charge):
    """(stripped sequence, charge, modified positions) -- the same for Skyline's PEPC[+57]K and
    DIA-NN's PEPC(UniMod:4)K. A mod before the first residue counts on residue 1, as Skyline
    writes an N-terminal one."""
    seq, pos, n = [], [], 0
    for part in _MOD.split(modified_seq or ""):
        if not part:
            continue
        if _MOD.fullmatch(part):
            pos.append(max(n, 1))
        else:
            letters = re.sub(r"[^A-Z]", "", part)
            seq.append(letters)
            n += len(letters)
    try:
        z = int(float(charge))
    except (TypeError, ValueError):
        z = None
    return ("".join(seq), z, tuple(sorted(pos)))


def key_label(key):
    seq, z, _ = key
    return f"{seq}/{z}+"


def run_of(name):
    """A run's name as DIA-NN writes it: the file's base name without .d/.raw/.mzML."""
    base = os.path.basename(str(name).replace("\\", "/").rstrip("/"))
    for ext in (".d", ".raw", ".mzML", ".mzml"):
        if base.lower().endswith(ext.lower()):
            return base[: -len(ext)]
    return base


def read_tsv(path):
    with open(path, newline="", encoding="utf-8", errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _col(row, *names):
    for nm in names:
        if nm in row:
            return row[nm]
    return None


def skyline_precursors(rows):
    """Skyline's precursor report -> [{key, run, accession, protein, apex, start, end, ...}]."""
    out = []
    for r in rows:
        if str(_col(r, "IsDecoy") or "").lower() == "true":
            continue
        prot = _col(r, "ProteinName") or ""
        key = precursor_key(_col(r, "PeptideModifiedSequence", "ModifiedSequence")
                            or _col(r, "PeptideSequence"), _col(r, "PrecursorCharge", "Charge"))
        out.append({"key": key, "run": run_of(_col(r, "FileName") or ""),
                    "accession": header_accession(prot) if prot else "", "protein": prot,
                    "modified": _col(r, "PeptideModifiedSequence", "ModifiedSequence") or "",
                    "apex": _float(_col(r, "BestRetentionTime")),
                    "start": _float(_col(r, "MinStartTime")), "end": _float(_col(r, "MaxEndTime")),
                    "area_fragment": _float(_col(r, "TotalAreaFragment")),
                    "area": _float(_col(r, "TotalArea")),
                    "dotp": _float(_col(r, "LibraryDotProduct")),
                    "idotp": _float(_col(r, "IsotopeDotProduct"))})
    return out


def diann_index(rows):
    """DIA-NN's rows -> {(key, run): row}; the first row of a duplicate wins."""
    out = {}
    for r in rows:
        key = precursor_key(r.get("Modified.Sequence") or r.get("Stripped.Sequence"),
                            r.get("Precursor.Charge"))
        d = {"group": r.get("Protein.Group"), "rt": _float(r.get("RT")),
             "start": _float(r.get("RT.Start")), "stop": _float(r.get("RT.Stop")),
             "quantity": _float(r.get("Precursor.Quantity"))}
        out.setdefault((key, r.get("Run")), d)
    return out


def disagreements(sky, dia, apex_diff_s=APEX_DIFF_S):
    """Every precursor-run both programs have, compared once (a peptide Skyline lists under two
    proteins counts once). Flagged when Skyline picked a different peak: no peak, its apex more
    than apex_diff_s seconds from DIA-NN's, or its apex outside DIA-NN's [RT.Start, RT.Stop].
    -> (rows compared, flagged rows), each row a dict ready for disagreements.tsv."""
    seen, compared, flagged = set(), [], []
    for s in sky:
        k = (s["key"], s["run"])
        if k in seen or k not in dia:
            continue
        seen.add(k)
        d = dia[k]
        row = {"protein": s["accession"], "precursor": key_label(s["key"]),
               "modified_sequence": s["modified"], "run": s["run"],
               "diann_apex_min": d["rt"], "diann_start_min": d["start"], "diann_end_min": d["stop"],
               "skyline_apex_min": s["apex"], "skyline_start_min": s["start"],
               "skyline_end_min": s["end"], "apex_difference_s": None, "skyline_dotp": s["dotp"],
               "reason": ""}
        why = []
        if s["apex"] is None:
            why.append("Skyline found no peak")
        elif d["rt"] is not None:
            diff = (s["apex"] - d["rt"]) * 60.0
            row["apex_difference_s"] = round(diff, 2)
            if abs(diff) > apex_diff_s:
                why.append(f"apex {abs(diff):.1f} s from DIA-NN's")
            if d["start"] is not None and d["stop"] is not None \
                    and not (d["start"] <= s["apex"] <= d["stop"]):
                why.append("apex outside DIA-NN's peak edges")
        row["reason"] = "; ".join(why)
        compared.append(row)
        if why:
            flagged.append(row)
    flagged.sort(key=lambda r: (r["apex_difference_s"] is not None,
                                -abs(r["apex_difference_s"] or 0), r["precursor"], r["run"]))
    return compared, flagged


DISAGREEMENT_COLUMNS = ["gene", "protein", "precursor", "modified_sequence", "run", "reason",
                        "apex_difference_s", "diann_apex_min", "skyline_apex_min",
                        "diann_start_min", "diann_end_min", "skyline_start_min",
                        "skyline_end_min", "skyline_dotp"]


def write_disagreements(path, rows, gene_of):
    with open(path + ".part", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(DISAGREEMENT_COLUMNS)
        for r in rows:
            r = dict(r, gene=gene_of.get(r["protein"], ""))
            w.writerow(["" if r.get(c) is None else r.get(c) for c in DISAGREEMENT_COLUMNS])
    os.replace(path + ".part", path)


def agreement_summary(compared):
    diffs = sorted(abs(r["apex_difference_s"]) for r in compared
                   if r["apex_difference_s"] is not None)
    if not diffs:
        return {"compared": len(compared), "apex_median_abs_s": None, "within_s_pct": None}
    mid = len(diffs) // 2
    med = diffs[mid] if len(diffs) % 2 else (diffs[mid - 1] + diffs[mid]) / 2
    return {"compared": len(compared), "apex_median_abs_s": round(med, 2),
            "within_s_pct": round(100.0 * sum(d <= APEX_DIFF_S for d in diffs) / len(diffs), 1)}


# ---- the pictures (matplotlib, inside the join job) --------------------------------------------------
# A fixed categorical order, never cycled; text in neutral ink (the trial's palette).
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
INK, INK2, GRID, SHADE, SURF = "#0b0b0b", "#52514e", "#e4e3df", "#d9d8d3", "#fcfcfb"


def safe_name(text):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text).strip("_") or "protein"


def runs_to_draw(runs, group_of, limit=MAX_RUNS_DRAWN):
    """Every run, or `limit` of them shared across the groups (round robin) when there are more."""
    runs = sorted(runs)
    if len(runs) <= limit:
        return runs
    by = {}
    for r in runs:
        by.setdefault(group_of.get(r, ""), []).append(r)
    out, i = [], 0
    while len(out) < limit:
        added = False
        for g in sorted(by):
            if i < len(by[g]) and len(out) < limit:
                out.append(by[g][i])
                added = True
        if not added:
            break
        i += 1
    return sorted(out, key=lambda r: (group_of.get(r, ""), r))


def chromatograms(rows):
    """Skyline's chromatogram export -> {(key, run): [(label, times, intensities)]}, fragments
    only (the MS1 isotopes and the precursor-as-product are left out of the pictures)."""
    out = {}
    for r in rows:
        ion = str(r.get("FragmentIon") or "")
        if ion.startswith("precursor"):
            continue
        key = precursor_key(r.get("PeptideModifiedSequence"), r.get("PrecursorCharge"))
        try:
            t = [float(x) for x in str(r.get("Times") or "").split(",") if x]
            y = [float(x) for x in str(r.get("Intensities") or "").split(",") if x]
        except ValueError:
            continue
        z = int(_float(r.get("ProductCharge")) or 1)
        out.setdefault((key, run_of(r.get("FileName") or "")), []).append(
            (f"{ion}{'+' * z}", t, y))
    return out


def draw_peaks(path, title, keys, runs, chrom, sky_by, dia, caption_runs):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 7, "axes.edgecolor": INK2, "axes.labelcolor": INK2,
                         "xtick.color": INK2, "ytick.color": INK2, "axes.titlesize": 7,
                         "figure.facecolor": SURF, "axes.facecolor": SURF})
    fig, axes = plt.subplots(len(keys), len(runs), squeeze=False,
                             figsize=(max(4.0, 2.3 * len(runs)), 1.9 * len(keys) + 0.8))
    for i, key in enumerate(keys):
        for j, run in enumerate(runs):
            a = axes[i][j]
            s = sky_by.get((key, run))
            d = dia.get((key, run))
            lo = hi = None
            if s and s["start"] is not None and s["end"] is not None:
                a.axvspan(s["start"], s["end"], color=SHADE, lw=0, zorder=0)
                lo, hi = s["start"], s["end"]
            for n, (lab, t, y) in enumerate((chrom.get((key, run)) or [])[:FRAGMENTS_DRAWN]):
                if t and len(t) == len(y):
                    a.plot(t, y, lw=1.0, color=SERIES[n], label=lab)
            if d and d["start"] is not None and d["stop"] is not None:
                for x in (d["start"], d["stop"]):
                    a.axvline(x, color=INK, lw=0.8, ls="--")
                lo = d["start"] if lo is None else min(lo, d["start"])
                hi = d["stop"] if hi is None else max(hi, d["stop"])
            if lo is not None:
                w = max(hi - lo, 0.1)
                a.set_xlim(lo - 1.5 * w, hi + 1.5 * w)
            note = [x for x, miss in (("DIA-NN: not reported in this run", d is None),
                                      ("Skyline: no peak", not s or s["apex"] is None)) if miss]
            if note:
                a.text(0.02, 0.97, "\n".join(note), transform=a.transAxes, fontsize=5.5,
                       color=INK2, va="top")
            a.grid(color=GRID, lw=0.5)
            a.set_axisbelow(True)
            for sp in ("top", "right"):
                a.spines[sp].set_visible(False)
            a.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
            if i == 0:
                a.set_title(run[-24:], color=INK2)
            if j == 0:
                a.set_ylabel(key_label(key), color=INK)
            if i == len(keys) - 1:
                a.set_xlabel("retention time (min)")
            if j == len(runs) - 1 and a.get_legend_handles_labels()[0]:
                a.legend(fontsize=5.5, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
    fig.suptitle(f"{title}\nFragment chromatograms extracted by Skyline. Grey band = the peak "
                 f"Skyline integrated; dashed lines = DIA-NN's peak edges.{caption_runs}",
                 color=INK, fontsize=7.5)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(path, dpi=110)
    plt.close(fig)


def draw_areas(path, title, keys, runs, tran_rows, sky_by):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(keys), squeeze=False,
                             figsize=(max(3.4, 0.35 * len(runs) + 1.6) * len(keys), 3.3))
    for a, key in zip(axes[0], keys):
        frag = {}
        for r in tran_rows:
            if str(r.get("FragmentIon") or "").startswith("precursor"):
                continue
            if precursor_key(r.get("PeptideModifiedSequence"), r.get("PrecursorCharge")) != key:
                continue
            ion = f"{r.get('FragmentIon')}{'+' * int(_float(r.get('ProductCharge')) or 1)}"
            run = run_of(r.get("FileName") or "")
            frag.setdefault(ion, {})[run] = frag.get(ion, {}).get(run, 0.0) + (_float(r.get("Area")) or 0.0)
        bottom = [0.0] * len(runs)
        for n, ion in enumerate(list(frag)[:FRAGMENTS_DRAWN]):
            vals = [frag[ion].get(r, 0.0) for r in runs]
            a.bar(range(len(runs)), vals, bottom=bottom, width=0.6, color=SERIES[n], label=ion,
                  edgecolor=SURF, linewidth=0.8)
            bottom = [b + v for b, v in zip(bottom, vals)]
        for x, r in enumerate(runs):
            s = sky_by.get((key, r))
            if s and s["dotp"] is not None:
                a.text(x, bottom[x], f"{s['dotp']:.2f}", ha="center", va="bottom", fontsize=5,
                       color=INK2)
        a.set_xticks(range(len(runs)))
        a.set_xticklabels([r[-14:] for r in runs], rotation=60, ha="right", fontsize=5.5)
        a.set_title(key_label(key), color=INK, fontsize=7)
        a.set_ylabel("fragment peak area", fontsize=7)
        a.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
        a.grid(axis="y", color=GRID, lw=0.6)
        a.set_axisbelow(True)
        for sp in ("top", "right"):
            a.spines[sp].set_visible(False)
        if frag:
            a.legend(fontsize=5.5, frameon=False, loc="upper left", bbox_to_anchor=(1.0, 1.0))
        a.set_ylim(0, (max(bottom) * 1.15) if max(bottom, default=0) > 0 else 1)
    fig.suptitle(f"{title}\nPeak areas by run from Skyline (stacked by fragment ion; Skyline's "
                 f"library dot product above each bar)", color=INK, fontsize=7.5)
    fig.tight_layout(rect=(0, 0, 1, 0.9))
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ---- the record and its words -------------------------------------------------------------------------
def record_path(out_dir):
    return os.path.join(out_dir, RECORD)


def load_record(out_dir):
    rec = load_json(record_path(out_dir))
    return rec if isinstance(rec, dict) and str(rec.get("schema", "")).startswith("skyline_peak_check/") \
        else None


def skipped_record(reason, stage, base=None):
    rec = dict(base or {})
    rec.update({"schema": SCHEMA, "status": "skipped", "stage": stage, "reason": reason,
                "written": now()})
    rec.setdefault("skill_version", skill_version())
    rec.setdefault("tool", {"name": "Skyline", "version": SKYLINE_PIN["version"],
                            "docker_digest": SKYLINE_PIN["docker_digest"]})
    return rec


def selection_sentence(rec):
    top = rec.get("top") or {}
    n = top.get("value", TOP_N_DEFAULT)
    tag = "" if top.get("source") == "user" else f" {default_tag()}"
    sel = rec.get("selection") or {}
    mode = sel.get("mode")
    if mode == "de_hits":
        cut = sel.get("adjp")
        cut_s = f"{fmt_num(cut)}" + ("" if sel.get("adjp_source") == "de_provenance.json"
                                     else f" {default_tag()}")
        return (f"The proteins checked are the up to {n}{tag} significant proteins with the "
                f"smallest adjusted p in any comparison (adjusted p < {cut_s}, the DE's own "
                f"cutoff; ties broken by the raw p-value; contaminant entries left out).")
    if mode == "abundance_no_hits":
        return (f"No protein was significant, so the proteins checked are the {n}{tag} most "
                f"abundant proteins measured in every run (the mean of the DE's log2 matrix; "
                f"contaminant entries left out). They were not chosen for a change.")
    if mode == "abundance_no_de":
        return (f"This search has no differential expression, so the proteins checked are the "
                f"{n}{tag} most abundant protein groups measured in every run (DIA-NN's PG.MaxLFQ, "
                f"at least 3 precursors; contaminant entries left out).")
    return f"The proteins checked: {SELECTION_RULE.format(n=n)}"


CAVEAT = ("These pictures are a visual check of DIA-NN's peaks. Skyline reuses DIA-NN's peak "
          "edges by default, so a peak that matches is not independent confirmation of it. The "
          "real warning signs are the peaks where Skyline picked a different peak.")


def methods_paragraph(rec):
    s, ma, t = rec["settings"], rec["mass_accuracy"], rec["tool"]
    res = rec.get("results") or {}
    fam = "TOF analysers, with the library's ion mobility" if rec.get("family") == "timstof" \
        else "centroided spectra"
    return (
        f"Peak check. For the proteins described below, chromatograms were extracted again from "
        f"the raw files with Skyline {t['version']} (MacLean et al., Bioinformatics 2010;26:966-968; "
        f"SkylineCmd in the ProteoWizard container {SKYLINE_PIN['image']}, image digest "
        f"{SKYLINE_PIN['docker_digest']}, run under Wine), from DIA-NN's spectral library for this "
        f"search and the search's digest settings ({s['enzyme']}, {s['missed_cleavages']} missed "
        f"cleavage(s), {s['min_length']}-{s['max_length']} residues, precursor charges "
        f"{s['charges'][0]}-{s['charges'][-1]}, fragments {fmt_num(s['fragment_mz'][0])}-"
        f"{fmt_num(s['fragment_mz'][1])} m/z), {fam}, with a tolerance of +-{fmt_num(ma['ms1_ppm'])} "
        f"ppm (MS1) and +-{fmt_num(ma['ms2_ppm'])} ppm (MS2): {ma['source']}. "
        f"{selection_sentence(rec)} Skyline integrates within DIA-NN's peak boundaries from the "
        f"library by default, so this is a visual check of DIA-NN's peaks, not an independent "
        f"identification; precursor-runs where Skyline's apex was more than "
        f"{fmt_num(APEX_DIFF_S)} s from DIA-NN's, or outside DIA-NN's peak boundaries, or where "
        f"Skyline found no peak, are listed"
        + (f" ({res.get('flagged', 0)} of {res.get('compared', 0)})." if res else "."))


def report_section(rec, rel="skyline"):
    """The report's "Peak check (Skyline)" section, from the record -> {"md", "figures",
    "captions"}. The wording is fixed here (CAVEAT), for every outcome; client_text's
    peak_check_overclaims() finds nothing in it (tests/test_skyline_peaks.py)."""
    status = rec.get("status")
    if status == "skipped":
        md = (f"**Not run.** {rec.get('reason') or 'No reason was recorded.'}\n\n"
              f"Without it, the peaks behind the results were not drawn; the quantities in this "
              f"report are DIA-NN's, as for any analysis.\n")
        return {"md": md, "figures": [], "captions": {}}
    if status != "ok":
        md = ("**Not finished when this report was made.** The Skyline jobs had been submitted "
              "but had not written their results; make the report again once they finish.\n")
        return {"md": md, "figures": [], "captions": {}}
    t, ma, res = rec["tool"], rec["mass_accuracy"], rec.get("results") or {}
    L = [f"Skyline {t['version']}, a second program, extracted the signal of the checked "
         f"proteins' peptides from the raw files again, using DIA-NN's spectral library and the "
         f"search's settings. **{CAVEAT}** Skyline picked a different peak when its apex is more "
         f"than {fmt_num(APEX_DIFF_S)} s from DIA-NN's, falls outside DIA-NN's peak edges, or "
         f"there is no peak at all.", "", selection_sentence(rec), ""]
    flagged, compared = res.get("flagged", 0), res.get("compared", 0)
    pct = f" ({100.0 * flagged / compared:.1f}%)" if compared else ""
    L += [f"**Result.** Of {compared} precursor-runs compared, Skyline picked a different peak in "
          f"{flagged}{pct}"
          + (f"; for the others the median apex difference was "
             f"{fmt_num(res['apex_median_abs_s'])} s." if res.get("apex_median_abs_s") is not None
             else ".")
          + (" Look at these before relying on the numbers they feed." if flagged else ""), ""]
    L += ["| Protein | Why checked | Peptides drawn | Precursor-runs compared | Skyline picked a "
          "different peak |", "|---|---|---|---|---|"]
    for p in rec.get("proteins") or []:
        why = (f"{p['contrast']}, adj. p {p['adj_p']:.2g}" if p.get("adj_p") is not None
               else "abundance")
        L.append(f"| {p['gene']} ({p['protein_group']}) | {why} | {p.get('drawn', 0)} | "
                 f"{p.get('compared', 0)} | {p.get('flagged', 0)} |")
    L.append("")
    worst = res.get("worst") or []
    if worst:
        L += [f"Where Skyline picked a different peak (the first {len(worst)} of {flagged}; all "
              f"in `{rel}/disagreements.tsv`):", "",
              "| Protein | Precursor | Run | Apex difference (s) | Why flagged |",
              "|---|---|---|---|---|"]
        for w in worst:
            d = w.get("apex_difference_s")
            L.append(f"| {w.get('gene') or w.get('protein')} | {w['precursor']} | {w['run']} | "
                     f"{'' if d is None else fmt_num(d)} | {w['reason']} |")
        L.append("")
    L += ["How to read the pictures: one row per peptide precursor and one column per run; each "
          "coloured line is a fragment ion; the grey band is the peak Skyline integrated and the "
          "dashed lines are DIA-NN's peak edges. A grey band away from the dashed lines is a peak "
          "to look at before trusting its number.", ""]
    figs, caps = [], {}
    for f in rec.get("figures") or []:
        if f.get("kind") != "peaks":
            continue
        figs.append(f["file"])
        caps[f["file"]] = f.get("caption") or ""
        L += [f"![{f.get('title') or f['file']}]({rel}/{f['file']})", ""]
    files = rec.get("files") or {}
    L += [f"In `{rel}/`: the Skyline document `{files.get('share_zip', '')}` (open it in "
          f"Skyline {t['version']} or later), Skyline's precursor report "
          f"`{files.get('precursor_report', '')}`, `{files.get('disagreements', '')}`, and a "
          f"peak-area chart per protein (`areas_<gene>.png`). Tolerance: +-{fmt_num(ma['ms1_ppm'])} "
          f"ppm (MS1) and +-{fmt_num(ma['ms2_ppm'])} ppm (MS2), {ma['source']}."]
    runs = rec.get("runs") or {}
    if runs.get("failed"):
        L += ["", "Runs Skyline could not import, left out of the pictures and tables: "
              + "; ".join(f"{os.path.basename(x['file'])} ({x['reason']})" for x in runs["failed"])
              + "."]
    if (rec.get("settings") or {}).get("not_carried_over"):
        L += ["", "Not carried over from the search: "
              + "; ".join(rec["settings"]["not_carried_over"]) + "."]
    return {"md": "\n".join(L) + "\n", "figures": figs, "captions": caps}


def report_data(output_dir, base=None):
    """For make_analysis_html: the section of output/skyline/'s record, its pictures referenced
    relative to `base` (the report's folder, default output/), or None without a record."""
    sky = os.path.join(output_dir, "skyline")
    rec = load_record(sky)
    if not rec:
        return None
    rel = os.path.relpath(sky, base or output_dir).replace(os.sep, "/")
    return report_section(rec, rel)


def manifest_line(output_dir, has_de):
    """finalize's MANIFEST.txt line -> (level, name, note), or None when the step does not apply
    (no DE and no record)."""
    name = "Peak check (Skyline, output/skyline)"
    rec = load_record(os.path.join(output_dir, "skyline"))
    if rec is None:
        if not has_de:
            return None
        return ("SKIPPED", name, "not run: no output/skyline/skyline_provenance.json (SKILL.md "
                                 "step 8f: skyline_peaks.py plan, or skip with the user's reason)")
    if rec.get("status") == "ok":
        res = rec.get("results") or {}
        return ("OK", name, f"{len(rec.get('proteins') or [])} protein(s), "
                            f"{len((rec.get('runs') or {}).get('imported') or [])} run(s); Skyline "
                            f"picked a different peak in {res.get('flagged', 0)} of "
                            f"{res.get('compared', 0)} precursor-runs")
    if rec.get("status") == "skipped":
        return ("SKIPPED", name, rec.get("reason") or "no reason recorded")
    return ("SKIPPED", name, "the Skyline jobs were submitted but have not written their results "
                             "(skyline_work/jobs.txt names them)")


# ---- plan ----------------------------------------------------------------------------------------------
def _session_paths(session):
    import session as sess
    return sess.paths_for(session)


def _first_existing(*paths):
    return next((p for p in paths if p and os.path.exists(p)), None)


def resolve_inputs(a):
    """-> dict of the inputs, from --session and the explicit flags (which win)."""
    p = _session_paths(a.session) if a.session else {}
    if not a.session and not (a.search and a.out):
        raise SystemExit("skyline_peaks.py plan: give --session, or --search and --out")
    out = os.path.abspath(a.out or p["skyline_dir"])
    search = os.path.abspath(a.search or p["search_out"])
    de_dir = a.de_dir or (p.get("de_dir") if p else None)
    manifest = a.manifest or _first_existing(p.get("workflow_manifest"),
                                             os.path.join(p.get("input_dir", ""),
                                                          "workflow.manifest.json"))
    conditions = a.conditions or (p.get("conditions") if p else None)
    name = a.name or (os.path.basename(os.path.abspath(a.session)) if a.session
                      else os.path.basename(os.path.dirname(search)))
    return {"out": out, "search": search, "de_dir": os.path.abspath(de_dir) if de_dir else None,
            "manifest": manifest, "conditions": conditions, "session": a.session,
            "name": safe_name(name)}


def find_speclib(search, report):
    stem = os.path.splitext(os.path.basename(report))[0]
    cand = os.path.join(search, f"{stem}-lib.parquet.skyline.speclib")
    if os.path.isfile(cand):
        return cand
    found = sorted(f for f in os.listdir(search) if f.endswith(".skyline.speclib"))
    if len(found) == 1:
        return os.path.join(search, found[0])
    raise Skip("the search folder has no DIA-NN library for Skyline (*-lib.parquet.skyline."
               "speclib, which DIA-NN 2.x writes beside the report)" if not found else
               f"the search folder has {len(found)} Skyline libraries and none named for "
               f"{os.path.basename(report)}")


def raw_files_of(search_prov, session_paths):
    files = [f for f in (search_prov or {}).get("files") or [] if f]
    if not files and session_paths and os.path.isfile(session_paths.get("raw_list", "")):
        with open(session_paths["raw_list"], encoding="utf-8") as fh:
            files = [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    return [f.rstrip("/") for f in files]


def check_python(py):
    """The jobs' Python: it exists here. Whether it has pyarrow and matplotlib is asked by the
    build job, on a compute node (importing them from group storage took over 2 minutes on a HIVE
    login node) -> None, or what is wrong."""
    path = shutil.which(py) if not os.path.isabs(py) else py
    if not path or not os.access(path, os.X_OK):
        return f"{py} is not an executable Python here"
    return None


def needs_for_pictures():
    """In the build job: Skip unless this Python can read the report and draw."""
    missing = []
    for mod in ("pyarrow", "matplotlib"):
        try:
            __import__(mod)
        except ImportError:
            missing.append(mod)
    if missing:
        raise Skip(f"the Python the jobs use ({sys.executable}) lacks {' and '.join(missing)}, "
                   f"which the peak check needs; `bash setup.sh` installs them into the pipeline "
                   f"env (or plan again with --python)")


def job_python(arg):
    if arg:
        return arg
    setup = load_json(os.path.expanduser("~/.proteomics-pipeline/setup.json")) or {}
    return setup.get("python") or sys.executable


def make_plan(a, inp):
    """Every check that can be made before a job runs; Skip with the reason when one fails."""
    search = inp["search"]
    if not os.path.isdir(search):
        raise Skip(f"there is no search folder at {search}")
    sprov = load_json(os.path.join(search, "search_provenance.json")) or {}
    if sprov.get("engine") and sprov.get("engine") != "diann":
        raise Skip(f"the search was run with {sprov.get('engine')}; the peak check reads DIA-NN's "
                   f"spectral library and report")
    report = os.path.join(search, "report.parquet")
    log = os.path.join(search, "report.log.txt")
    for f, what in ((report, "DIA-NN's report.parquet"), (log, "DIA-NN's report.log.txt")):
        if not os.path.isfile(f):
            raise Skip(f"the search folder has no {what}")
    with open(log, encoding="utf-8", errors="replace") as fh:
        log_text = fh.read()
    ver = diann_version(log_text)
    if not ver or int(ver.split(".")[0]) < 2:
        raise Skip(f"the search was not run with DIA-NN 2.x (its log says {ver or 'no version'}); "
                   f"Skyline {SKYLINE_PIN['version']} reads DIA-NN 2.x's parquet report and library")
    speclib = find_speclib(search, report)
    flags = command_flags(diann_command(log_text) or "")
    settings = search_settings(flags)
    ma = mass_accuracy(log_text, sprov, a.ppm_ms1, a.ppm_ms2, a.ppm_source)
    sp = _session_paths(inp["session"]) if inp["session"] else None
    raw = raw_files_of(sprov, sp)
    if not raw:
        raise Skip("the raw files the search read are not recorded (search_provenance.json "
                   "files, or input/raw_files.txt)")
    missing = [f for f in raw if not os.path.exists(f)]
    if missing:
        raise Skip(f"{len(missing)} of the {len(raw)} raw files the search read are not reachable "
                   f"here (first: {missing[0]})")
    manifest = load_json(inp["manifest"]) if inp["manifest"] else None
    fam, fam_src = instrument_family(manifest, raw)
    fasta = sprov.get("fasta") or (flags.get("fasta") or [None])[-1]
    if not fasta or not os.path.isfile(fasta):
        raise Skip(f"the FASTA the search read is not reachable here ({fasta or 'not recorded'})")
    sif = os.path.abspath(a.sif or SKYLINE_PIN["sif"])
    if not os.path.isfile(sif):
        raise Skip(f"the pinned Skyline container is not reachable here ({sif}). It is the "
                   f"ProteoWizard image {SKYLINE_PIN['image']}@{SKYLINE_PIN['docker_digest']}; "
                   f"references/skyline.md says how to pull it and pass --sif")
    import on_compute
    if not (shutil.which("sbatch") or on_compute.on_login_node()):
        raise Skip("there is no SLURM cluster here: the peak check imports one raw file per "
                   "SLURM job, in the ProteoWizard container")
    py = job_python(a.python)
    bad = check_python(py)
    if bad:
        raise Skip(f"the Python the jobs would use is not usable ({bad}); pass --python")
    de_dir = inp["de_dir"] if inp["de_dir"] and os.path.isfile(
        os.path.join(inp["de_dir"], "de_provenance.json")) else None
    de_prov = load_json(os.path.join(de_dir, "de_provenance.json")) if de_dir else {}
    q = _float((de_prov or {}).get("q_cutoff")) or 0.01
    out, work = inp["out"], os.path.join(inp["out"], WORK)
    doc = os.path.join(work, f"{inp['name']}_top_hits.sky")
    group_of = {}
    if inp["conditions"] and os.path.isfile(inp["conditions"]):
        with open(inp["conditions"], newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                if r.get("File.Name"):
                    group_of[run_of(r["File.Name"])] = r.get("Group") or ""
    plan = {
        "schema": SCHEMA, "planned": now(), "skill_version": skill_version(),
        "session": os.path.abspath(inp["session"]) if inp["session"] else None,
        "name": inp["name"], "out": out, "work": work, "doc": doc,
        "share_zip": os.path.join(out, f"{inp['name']}_top_hits.sky.zip"),
        "search": {"dir": search, "report": report, "speclib": speclib, "log": log,
                   "engine": "diann", "version": ver},
        "report": report, "speclib_link": os.path.join(work, "search", os.path.basename(speclib)),
        "fasta": fasta, "targets_fasta": os.path.join(work, "targets.fasta"),
        "de_dir": de_dir, "q_cutoff": q, "raw_files": raw, "group_of": group_of,
        "family": fam, "family_source": fam_src, "settings": settings, "mass_accuracy": ma,
        "top": {"value": a.top if a.top else TOP_N_DEFAULT,
                "source": "user" if a.top else "skill default"},
        "sif": sif, "ro_binds": ro_binds(search, raw, fasta),
        "python": py, "script": os.path.abspath(__file__),
        "max_parallel": a.max_parallel or MAX_PARALLEL_DEFAULT,
    }
    return plan


def write_jobs(plan, queue):
    w = plan["work"]
    os.makedirs(os.path.join(w, "logs"), exist_ok=True)
    os.makedirs(os.path.join(w, "status"), exist_ok=True)
    os.makedirs(os.path.join(w, "search"), exist_ok=True)
    # DIA-NN's library and its report, seen from the document's folder (Skyline finds the
    # report beside the library by name); the search folder itself stays read-only
    for src in (plan["search"]["speclib"], plan["report"],
                os.path.join(plan["search"]["dir"], "report-lib.parquet")):
        if os.path.isfile(src):
            link = os.path.join(w, "search", os.path.basename(src))
            if os.path.lexists(link):
                os.remove(link)
            os.symlink(os.path.realpath(src), link)
    for fn, text in SKYR.items():
        write_text(os.path.join(w, fn), text)
    n = len(plan["raw_files"])
    write_text(os.path.join(w, "raw_files.txt"), "\n".join(plan["raw_files"]) + "\n")
    scripts = {"sky_1_build.sbatch": job_script(plan, "build", queue),
               "sky_2_import.sbatch": job_script(plan, "import", queue,
                                                 array=f"1-{n}%{min(n, plan['max_parallel'])}"),
               "sky_3_join.sbatch": job_script(plan, "join", queue),
               "submit.sh": submit_script(plan)}
    for fn, text in scripts.items():
        write_text(os.path.join(w, fn), text)
        os.chmod(os.path.join(w, fn), 0o755)
    write_json(os.path.join(w, "plan.json"), plan)


def planned_record(plan):
    return {"schema": SCHEMA, "status": "planned", "written": now(),
            "skill_version": plan["skill_version"], "name": plan["name"],
            "tool": {"name": "Skyline", "version": SKYLINE_PIN["version"],
                     "docker_digest": SKYLINE_PIN["docker_digest"], "container": plan["sif"]},
            "top": plan["top"], "family": plan["family"], "family_source": plan["family_source"],
            "settings": plan["settings"], "mass_accuracy": plan["mass_accuracy"],
            "search": plan["search"], "runs": {"n": len(plan["raw_files"])},
            "work": plan["work"]}


def cmd_plan(a):
    inp = resolve_inputs(a)
    out = inp["out"]
    try:
        plan = make_plan(a, inp)
        import run_search
        queue = run_search.slurm_queue(a.partition, a.account, a.qos, peak_cpus=JOB_CPUS,
                                       preemptible_ok=True)
        write_jobs(plan, queue)
    except Skip as e:
        write_json(record_path(out), skipped_record(str(e), "plan"))
        print(json.dumps({"status": "skipped", "reason": str(e), "record": record_path(out)},
                         indent=2))
        return 3
    write_json(record_path(out), planned_record(plan))
    submit = os.path.join(plan["work"], "submit.sh")
    print(json.dumps({"status": "planned", "record": record_path(out), "work": plan["work"],
                      "runs": len(plan["raw_files"]), "family": plan["family"],
                      "mass_accuracy": plan["mass_accuracy"], "settings": plan["settings"],
                      "queue": {"partition": queue[0], "account": queue[1], "qos": queue[2]},
                      "submit with": f"bash {shlex.quote(submit)}"}, indent=2))
    return 0


def cmd_skip(a):
    out = os.path.abspath(a.out or _session_paths(a.session)["skyline_dir"])
    reason = " ".join((a.reason or "").split())
    if not reason:
        raise SystemExit("skyline_peaks.py skip: --reason is required (why, in the user's words)")
    write_json(record_path(out), skipped_record(f"Turned off for this analysis: {reason}", "user"))
    print(json.dumps({"status": "skipped", "record": record_path(out)}, indent=2))
    return 0


def cmd_status(a):
    out = os.path.abspath(a.out or _session_paths(a.session)["skyline_dir"])
    rec = load_record(out)
    jobs = os.path.join(out, WORK, "jobs.txt")
    print(json.dumps({"record": record_path(out), "status": (rec or {}).get("status"),
                      "reason": (rec or {}).get("reason"),
                      "results": (rec or {}).get("results"),
                      "jobs": open(jobs).read().split() if os.path.isfile(jobs) else None},
                     indent=2))
    return 0 if rec else 1


# ---- inside the jobs -----------------------------------------------------------------------------------
def _plan_of(work):
    plan = load_json(os.path.join(work, "plan.json"))
    if not isinstance(plan, dict):
        raise SystemExit(f"skyline_peaks.py: no plan.json in {work}")
    return plan


def run_skyline(plan, args, log_name, minutes, rw=()):
    """SkylineCmd in the container, with a temporary folder of its own on the node (removed
    after) -> (exit code, log path). Never --import-process-count."""
    assert not any(x.startswith("--import-process-count") for x in args)
    import tempfile
    log = os.path.join(plan["work"], "logs", log_name)
    base = os.environ.get("TMPDIR")
    tmp = tempfile.mkdtemp(prefix="skyline_", dir=base if base and os.path.isdir(base) else None)
    argv = apptainer_argv(plan, rw, tmp) + ["wine", "SkylineCmd"] + list(args)
    with open(os.path.join(plan["work"], "logs", log_name.replace(".txt", ".cmd.txt")), "w") as fh:
        fh.write(" ".join(shlex.quote(x) for x in argv) + "\n")
    try:
        with open(log, "w") as fh:
            rc = subprocess.run(argv, stdout=fh, stderr=subprocess.STDOUT,
                                timeout=max(60, (minutes - 5) * 60)).returncode
    except subprocess.TimeoutExpired:
        rc = 124
    except OSError as e:
        with open(log, "a") as fh:
            fh.write(f"cannot start apptainer: {e}\n")
        rc = 127
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return rc, log


def log_tail(path, n=3):
    """What went wrong, from a SkylineCmd log: its "Error:" line and the exception's own message
    when there are any (a .NET stack trace is not a reason), else the last few lines."""
    try:
        with open(path, errors="replace") as fh:
            lines = [ln.strip() for ln in fh if ln.strip() and not re.match(
                r"^(INFO:|wineserver|libgcc|wine client|Settings|\[\d+\] \d+%|\d+%)", ln.strip())]
    except OSError:
        return "no log"
    said = [ln for ln in lines if ln.startswith("Error:")]
    said += [m.group(1).strip() for ln in lines
             for m in [re.search(r"---> [\w.]+Exception: (.*)", ln)] if m]
    if said:
        return " | ".join(dict.fromkeys(said))[:400]
    return " | ".join(lines[-n:])[:400] or "the log is empty"


def verify_container(plan):
    """The Core's .sif by its pinned sha256; any .sif by the version SkylineCmd reports."""
    sif = plan["sif"]
    got = sha256_file(sif)
    pinned = os.path.basename(sif) == os.path.basename(SKYLINE_PIN["sif"])
    if pinned and got != SKYLINE_PIN["sif_sha256"]:
        raise Skip(f"the Skyline container {sif} does not match its pinned sha256 (it is "
                   f"{got[:12]}..., the pin is {SKYLINE_PIN['sif_sha256'][:12]}...): it is not "
                   f"the image the skill was tested with")
    rc, log = run_skyline(plan, ["--version"], "skyline_version.txt", 10)
    with open(log, errors="replace") as fh:
        text = fh.read()
    m = re.search(r"Skyline[^\n]*?(\d+\.\d+\.\d+\.\d+)\s*\(([0-9a-f]+)\)", text)
    if not m:
        raise Skip(f"SkylineCmd --version did not answer in the container ({log_tail(log)})")
    if m.group(1) != SKYLINE_PIN["version"]:
        raise Skip(f"the container holds Skyline {m.group(1)}, not the pinned "
                   f"{SKYLINE_PIN['version']}")
    return {"name": "Skyline", "version": m.group(1), "commit": m.group(2),
            "version_reported": m.group(0).strip(), "container": sif, "container_sha256": got,
            "verified": ("sha256 matches the pin" if pinned else
                         "a .sif other than the Core's: its sha256 is recorded, and its version "
                         "is checked from SkylineCmd --version"),
            "docker_digest": SKYLINE_PIN["docker_digest"], "image": SKYLINE_PIN["image"]}


def _fail(plan, stage, reason, extra=None):
    base = {"skill_version": plan.get("skill_version"), "name": plan.get("name"),
            "top": plan.get("top"), "family": plan.get("family"),
            "settings": plan.get("settings"), "mass_accuracy": plan.get("mass_accuracy"),
            "search": plan.get("search"), "work": plan.get("work")}
    base.update(extra or {})
    write_json(record_path(plan["out"]), skipped_record(reason, stage, base))


def cmd_build(a):
    plan = _plan_of(a.work)
    w = plan["work"]
    st = os.path.join(w, "status", "build.json")
    try:
        needs_for_pictures()
        tool = verify_container(plan)
        sel = select_targets(plan)
        if not sel["targets"]:
            raise Skip("no protein qualified for the check (" + {
                "de_hits": "no significant protein",
                "abundance_no_hits": "no protein measured in every run",
                "abundance_no_de": "no protein group measured in every run with 3 precursors",
            }.get(sel["mode"], sel["mode"]) + ")")
        accs = sorted({x for t in sel["targets"] for x in accessions(t["protein_group"])})
        found = fasta_subset(plan["fasta"], accs, plan["targets_fasta"])
        lost = [t["protein_group"] for t in sel["targets"]
                if not any(x in found for x in accessions(t["protein_group"]))]
        sel["targets"] = [t for t in sel["targets"] if t["protein_group"] not in lost]
        if not sel["targets"]:
            raise Skip("none of the chosen proteins is in the searched FASTA by accession")
        sel["not_in_fasta"] = lost
        nrows, last_rt = diann_rows(plan["report"], [t["protein_group"] for t in sel["targets"]],
                                    plan["q_cutoff"], os.path.join(w, "diann_targets.tsv"))
        rt_min, rt_src = rt_window(last_rt)
        args = build_args(plan, rt_min)
        write_text(os.path.join(w, "build_args.txt"), "\n".join(args) + "\n")
        rc, log = run_skyline(plan, args, "build.txt", JOB_MINUTES["build"])
        if rc != 0 or not os.path.isfile(plan["doc"]):
            raise Skip(f"Skyline could not build the document (exit {rc}: {log_tail(log)})")
        with open(plan["doc"], errors="replace") as fh:
            n_pep = len(re.findall(r"<peptide\b", fh.read()))
        if not n_pep:
            raise Skip("the Skyline document has no peptides: no DIA-NN library peptide matched "
                       "the chosen proteins' digest")
    except Skip as e:
        _fail(plan, "build", str(e))
        write_json(st, {"ok": False, "reason": str(e)})
        print(f"[skyline_peaks] build: not run -- {e}", file=sys.stderr)
        return 1
    write_json(st, {"ok": True, "tool": tool, "selection": sel, "diann_rows": nrows,
                    "rt_window": {"minutes": rt_min, "source": rt_src}, "peptides": n_pep,
                    "finished": now()})
    print(f"[skyline_peaks] build: {len(sel['targets'])} protein(s), {n_pep} peptide(s)",
          file=sys.stderr)
    return 0


def cmd_import(a):
    plan = _plan_of(a.work)
    i = int(a.task)
    raw = plan["raw_files"][i - 1]
    rc, log = run_skyline(plan, import_args(plan, i), f"import_{i}.txt", JOB_MINUTES["import"])
    write_json(os.path.join(plan["work"], "status", f"import_{i}.json"),
               {"file": raw, "rc": rc, "reason": None if rc == 0 else log_tail(log),
                "finished": now()})
    return 0 if rc == 0 else 1


def cmd_join(a):
    plan = _plan_of(a.work)
    w = plan["work"]
    build = load_json(os.path.join(w, "status", "build.json"))
    if not build or not build.get("ok"):
        if not build:                  # the build job died before it could say why
            _fail(plan, "build", "the build job ended without a result (see "
                                 f"{os.path.join(w, 'logs')})")
        return 0                       # the build wrote the record of why
    imported, failed = [], []
    for i, f in enumerate(plan["raw_files"], 1):
        s = load_json(os.path.join(w, "status", f"import_{i}.json"))
        if s and s.get("rc") == 0:
            imported.append(f)
        else:
            failed.append({"file": f, "reason": (s or {}).get("reason")
                           or "its import job ended without a result"})
    if not imported:
        _fail(plan, "import", f"Skyline could not import any of the {len(failed)} raw files "
                              f"(first: {os.path.basename(failed[0]['file'])}: {failed[0]['reason']})")
        return 1
    try:
        batch = join_batch(plan, imported)
        write_text(os.path.join(w, "join_batch.txt"), batch)
        rc, log = run_skyline(plan, [f"--batch-commands={win(os.path.join(w, 'join_batch.txt'))}"],
                              "join.txt", JOB_MINUTES["join"] - 15)
        need = {"precursor report": os.path.join(plan["out"], "skyline_precursors.tsv"),
                "transition report": os.path.join(w, "skyline_transitions.tsv"),
                "chromatograms": os.path.join(w, "skyline_chromatograms.tsv"),
                "Skyline document (.sky.zip)": plan["share_zip"]}
        gone = [k for k, p in need.items() if not (os.path.isfile(p) and os.path.getsize(p))]
        if gone:
            raise Skip(f"Skyline did not write the {', '.join(gone)} (exit {rc}: {log_tail(log)})")
        rec = finish(plan, build, imported, failed)
    except Skip as e:
        _fail(plan, "join", str(e), {"runs": {"imported": imported, "failed": failed}})
        return 1
    write_json(record_path(plan["out"]), rec)
    print(f"[skyline_peaks] done: {rec['results']['flagged']} of {rec['results']['compared']} "
          f"precursor-runs where Skyline picked a different peak", file=sys.stderr)
    return 0


def finish(plan, build, imported, failed):
    """Disagreements, pictures and the record, from Skyline's exports and DIA-NN's rows."""
    w, out = plan["work"], plan["out"]
    sky = skyline_precursors(read_tsv(os.path.join(out, "skyline_precursors.tsv")))
    dia = diann_index(read_tsv(os.path.join(w, "diann_targets.tsv")))
    targets = build["selection"]["targets"]
    acc_to_group = {}
    for t in targets:
        for x in accessions(t["protein_group"]):
            acc_to_group.setdefault(x, t["protein_group"])
    gene_of = {t["protein_group"]: t["gene"] for t in targets}
    gene_by_acc = {x: gene_of[g] for x, g in acc_to_group.items()}
    compared, flagged = disagreements(sky, dia)
    write_disagreements(os.path.join(out, "disagreements.tsv"), flagged, gene_by_acc)
    sky_by = {}
    for s in sky:
        sky_by.setdefault((s["key"], s["run"]), s)
    chrom = chromatograms(read_tsv(os.path.join(w, "skyline_chromatograms.tsv")))
    tran = read_tsv(os.path.join(w, "skyline_transitions.tsv"))
    runs_all = sorted({run_of(f) for f in imported})
    runs = runs_to_draw(runs_all, plan.get("group_of") or {})
    shown = "" if len(runs) == len(runs_all) else \
        f" {len(runs)} of {len(runs_all)} runs shown, chosen across the groups."
    flag_keys = {(r["precursor"], r["run"]) for r in flagged}
    figures, proteins, used = [], [], set()
    for t in targets:
        accs = set(accessions(t["protein_group"]))
        mine = [s for s in sky if s["accession"] in accs]
        area = {}
        for s in mine:
            area[s["key"]] = area.get(s["key"], 0.0) + (s["area_fragment"] or 0.0)
        in_dia = {k for (k, _r), d in dia.items() if d["group"] == t["protein_group"]}
        keys = sorted(area, key=lambda k: (k not in in_dia, -area[k], key_label(k)))[:PRECURSORS_DRAWN]
        n_cmp = sum(1 for r in compared if r["protein"] in accs)
        n_flag = sum(1 for r in flagged if r["protein"] in accs)
        p = dict(t, drawn=len(keys), compared=n_cmp, flagged=n_flag)
        proteins.append(p)
        if not keys:
            continue
        base = safe_name(t["gene"])
        if base in used:
            base = safe_name(f"{t['gene']}_{t['protein_group'].split(';')[0]}")
        used.add(base)
        title = f"{t['gene']} ({t['protein_group']})"
        draw_peaks(os.path.join(out, f"peaks_{base}.png"), title, keys, runs, chrom, sky_by, dia,
                   shown)
        draw_areas(os.path.join(out, f"areas_{base}.png"), title, keys, runs_all, tran, sky_by)
        why = (f"{t['contrast']}, adjusted p {t['adj_p']:.2g}" if t.get("adj_p") is not None
               else "one of the most abundant proteins")
        nf = sum(1 for k in keys for r in runs if (key_label(k), r) in flag_keys)
        figures.append({"kind": "peaks", "file": f"peaks_{base}.png", "protein": t["protein_group"],
                        "title": f"Peaks of {title}",
                        "caption": (f"Chosen as {why}. Its {len(keys)} strongest peptide "
                                    f"precursors in each run, as Skyline extracted them. Grey band "
                                    f"= the peak Skyline integrated; dashed lines = DIA-NN's peak "
                                    f"edges. Skyline picked a different peak in {nf} of these "
                                    f"panels.{shown}")})
        figures.append({"kind": "areas", "file": f"areas_{base}.png", "protein": t["protein_group"],
                        "title": f"Peak areas of {title}"})
    if flagged:
        worst_keys, seen = [], set()
        by_label = {key_label(s["key"]): s["key"] for s in sky}
        for r in flagged:
            k = by_label.get(r["precursor"])
            if k and k not in seen:
                seen.add(k)
                worst_keys.append(k)
            if len(worst_keys) == PRECURSORS_DRAWN:
                break
        if worst_keys:
            draw_peaks(os.path.join(out, "peaks_disagreements.png"),
                       "Where Skyline picked a different peak (the largest differences)",
                       worst_keys, runs, chrom, sky_by, dia, shown)
            figures.append({"kind": "peaks", "file": "peaks_disagreements.png",
                            "title": "Where Skyline picked a different peak",
                            "caption": "The precursors with the largest difference between "
                                       "Skyline's apex and DIA-NN's: grey band = Skyline's peak, "
                                       "dashed lines = DIA-NN's peak edges." + shown})
    summ = agreement_summary(compared)
    worst = [dict(r, gene=gene_by_acc.get(r["protein"], "")) for r in flagged[:15]]
    rec = {"schema": SCHEMA, "status": "ok", "written": now(),
           "skill_version": plan.get("skill_version"), "name": plan["name"],
           "tool": build["tool"], "top": plan["top"],
           "selection": {k: v for k, v in build["selection"].items() if k != "targets"},
           "proteins": proteins, "family": plan["family"], "family_source": plan["family_source"],
           "settings": dict(plan["settings"], rt_window_min=build["rt_window"]["minutes"],
                            rt_window_source=build["rt_window"]["source"],
                            analyzer_args=analyzer_args(plan["family"],
                                                        plan["mass_accuracy"]["ms1_ppm"],
                                                        plan["mass_accuracy"]["ms2_ppm"]),
                            product_ions=PRODUCT_ION_ARGS),
           "mass_accuracy": plan["mass_accuracy"], "search": plan["search"],
           "runs": {"n": len(plan["raw_files"]), "imported": imported, "failed": failed,
                    "drawn": runs},
           "results": dict(summ, flagged=len(flagged), apex_diff_s=APEX_DIFF_S, worst=worst),
           "figures": figures,
           "files": {"precursor_report": "skyline_precursors.tsv",
                     "disagreements": "disagreements.tsv",
                     "share_zip": os.path.basename(plan["share_zip"])},
           "work": w}
    rec["methods_paragraph"] = methods_paragraph(rec)
    return rec


# ---- the follow-up assay (PRM/SRM): Core staff, on request --------------------------------------------
def core_staff():
    """(login, None) when the account running this is on the Core staff list, else (None, why).
    Fails closed: no usable list is no (make_podcast.core_staff's rule)."""
    import pwd
    import staff
    login = pwd.getpwuid(os.getuid()).pw_name
    logins, why = staff.load_staff()
    if logins is None:
        return None, f"the Core staff list cannot be used ({why})"
    if login not in logins:
        return None, f"{login} is not on the Core staff list ({staff.staff_file()})"
    return login, None


def prm_plan(rec, keys, same_gradient, run_length):
    """-> [{key, kind, instrument, scheduled, why, file, extra args}] or raises ValueError."""
    out = []
    for k in keys:
        if k in PRM_NOT_OFFERED:
            raise ValueError(f"Skyline {SKYLINE_PIN['version']}'s command line has no isolation "
                             f"list for the {PRM_NOT_OFFERED[k]} (its --help offers Bruker "
                             f"timsTOF, Thermo Fusion, Thermo Q Exactive and Thermo Stellar); "
                             f"only a native method export, which needs a method template from "
                             f"the instrument PC")
        if k not in PRM_INSTRUMENTS:
            raise ValueError(f"unknown --instrument {k}: one of {', '.join(PRM_INSTRUMENTS)}")
        kind, name, fam = PRM_INSTRUMENTS[k]
        same_family = fam == rec.get("family")
        scheduled = same_family and same_gradient
        if scheduled:
            why = (f"scheduled: the same instrument family as the search ({fam}) and, as staff "
                   f"stated, the same gradient; windows from the measured retention times")
        elif not same_family:
            why = (f"unscheduled: the search ran on a {rec.get('family')} instrument and this list "
                   f"is for a {fam} one, so its retention times do not transfer")
        else:
            why = ("unscheduled: the same instrument family, but the same gradient was not "
                   "confirmed (--same-gradient)")
        extra = ["--exp-method-type=Scheduled"] if scheduled else ["--exp-method-type=Standard"]
        if not scheduled and name.startswith("Thermo"):
            if not run_length:
                raise ValueError(f"an unscheduled {name} list needs the new method's run length: "
                                 f"--run-length <minutes> (the length of ITS gradient, not the "
                                 f"search's)")
            extra.append(f"--exp-run-length={fmt_num(run_length)}")
        label = ("SRM" if kind == "transition" else "PRM") + "_" + safe_name(name.replace(" ", ""))
        out.append({"key": k, "kind": kind, "instrument": name, "family": fam,
                    "scheduled": scheduled, "why": why, "extra": extra,
                    "file": f"{label}_{'scheduled' if scheduled else 'unscheduled'}.csv"})
    return out


def prm_batch(plan, doc, targets, dest_dir):
    """The follow-up assay as one --batch-commands file: refine to the top 3 peptides, 6
    fragments and one charge per protein from the measured areas, switch DIA -> PRM, export each
    isolation list; for SRM, drop the "p" ion type and predict CE before the same refinement."""
    w = os.path.join(plan["work"], "prm")
    refined, prm, srm = (os.path.join(w, f) for f in ("refined.sky", "prm.sky", "srm.sky"))
    lines = []
    iso = [t for t in targets if t["kind"] == "isolation"]
    if iso:
        lines.append([f"--in={win(doc)}"] + PRM_REFINE_ARGS + [f"--out={win(refined)}"])
        lines.append([f"--in={win(refined)}", "--full-scan-acquisition-method=PRM",
                      "--full-scan-precursor-isotopes=None", f"--out={win(prm)}"])
        for t in iso:
            lines.append([f"--in={win(prm)}", f"--exp-isolationlist-instrument={t['instrument']}"]
                         + t["extra"] + ["--exp-strategy=Single",
                                         f"--exp-file={win(os.path.join(dest_dir, t['file']))}"])
    for t in (t for t in targets if t["kind"] == "transition"):
        lines.append([f"--in={win(doc)}", "--tran-product-ion-types=y, b",
                      f"--tran-predict-ce={SRM_CE_PREDICTOR}"] + PRM_REFINE_ARGS
                     + [f"--out={win(srm)}"])
        lines.append([f"--in={win(srm)}", f"--exp-translist-instrument={t['instrument']}"]
                     + t["extra"] + ["--exp-strategy=Single",
                                     f"--exp-file={win(os.path.join(dest_dir, t['file']))}"])
    return "\n".join(batch_line(x) for x in lines) + "\n"


def cmd_prm(a):
    out = os.path.abspath(a.out or _session_paths(a.session)["skyline_dir"])
    rec = load_record(out)
    if not rec or rec.get("status") != "ok":
        raise SystemExit("skyline_peaks.py prm: the peak check has not finished for this session "
                         "(its record is not 'ok'); the follow-up assay is made from its document")
    if not a.no_staff_check:
        login, why = core_staff()
        if not login:
            raise SystemExit(f"skyline_peaks.py prm: Core staff only -- {why}")
    keys = [k.strip() for k in (a.instrument or "").split(",") if k.strip()]
    if not keys:
        raise SystemExit("skyline_peaks.py prm: --instrument <key>[,<key>...]: "
                         + ", ".join(PRM_INSTRUMENTS))
    try:
        targets = prm_plan(rec, keys, a.same_gradient, a.run_length)
    except ValueError as e:
        raise SystemExit(f"skyline_peaks.py prm: {e}")
    plan = _plan_of(rec.get("work") or os.path.join(out, WORK))
    dest = os.path.abspath(a.dest or os.path.join(os.path.dirname(out), FOLLOWUP_DIR))
    os.makedirs(os.path.join(plan["work"], "prm"), exist_ok=True)
    os.makedirs(dest, exist_ok=True)
    spec = {"targets": targets, "dest": dest, "requested": now(),
            "batch": prm_batch(plan, plan["doc"], targets, dest)}
    write_json(os.path.join(plan["work"], "prm", "prm.json"), spec)
    write_text(os.path.join(plan["work"], "prm", "prm_batch.txt"), spec["batch"])
    import run_search
    queue = run_search.slurm_queue(a.partition, a.account, a.qos, peak_cpus=JOB_CPUS,
                                   preemptible_ok=True)
    job = job_script(plan, "prm", queue)
    path = os.path.join(plan["work"], "sky_prm.sbatch")
    write_text(path, job)
    os.chmod(path, 0o755)
    print(json.dumps({"targets": [{k: t[k] for k in ("instrument", "scheduled", "why", "file")}
                                  for t in targets], "dest": dest,
                      "submit with": f"sbatch {shlex.quote(path)}"}, indent=2))
    return 0


def followup_readme(spec, rec, counts):
    L = ["# Targeted follow-up assay (PRM / SRM) from the peak check", "",
         f"Made from the peak-check Skyline document ({rec.get('name')}): for each protein the 3 "
         f"peptides with the largest measured areas, their 6 largest fragments and one charge "
         f"state (Skyline's refinement on the DIA runs), with Skyline {SKYLINE_PIN['version']}.",
         "", "| File | Instrument | Targets (rows) | Retention-time scheduling |", "|---|---|---|---|"]
    for t in spec["targets"]:
        L.append(f"| {t['file']} | {t['instrument']} | {counts.get(t['file'], '?')} | {t['why']} |")
    L += ["", "Notes:",
          "- A scheduled list uses Skyline's default window around the measured retention time; "
          "SkylineCmd has no option to change it.",
          "- The Bruker timsTOF list carries the 1/K0 window from the search's runs and leaves the "
          "collision energy at 0: the prm-PASEF method sets it.",
          "- The Thermo isolation lists give Skyline's default normalised collision energy; the "
          "SRM list's collision energies are predicted per precursor "
          f"(\"{SRM_CE_PREDICTOR}\"), and the precursor-as-product (\"p\") ion type was left out "
          "before the transitions were chosen.",
          "- The proteins and peptides come from this study's DIA data. They are a starting point "
          "for an assay; each list still needs a test injection on the target instrument.", ""]
    return "\n".join(L)


def cmd_prm_run(a):
    plan = _plan_of(a.work)
    spec = load_json(os.path.join(plan["work"], "prm", "prm.json"))
    rec = load_record(plan["out"]) or {}
    rc, log = run_skyline(plan, [f"--batch-commands={win(os.path.join(plan['work'], 'prm', 'prm_batch.txt'))}"],
                          "prm.txt", JOB_MINUTES["prm"], rw=[spec["dest"]])
    counts, srm_self = {}, 0
    for t in spec["targets"]:
        p = os.path.join(spec["dest"], t["file"])
        if os.path.isfile(p):
            with open(p, newline="", encoding="utf-8", errors="replace") as fh:
                rows = list(csv.reader(fh))
            counts[t["file"]] = max(len(rows) - 1, 0)
            if t["kind"] == "transition" and rows:
                h = rows[0]
                try:
                    ip, iq = h.index("Precursor (m/z)"), h.index("Product (m/z)")
                    srm_self = sum(1 for r in rows[1:] if abs(float(r[ip]) - float(r[iq])) < 0.01)
                except (ValueError, IndexError):
                    pass
        else:
            counts[t["file"]] = "not written"
    write_text(os.path.join(spec["dest"], "README.md"), followup_readme(spec, rec, counts))
    write_json(os.path.join(spec["dest"], "followup_assay.json"),
               {"schema": "skyline_followup/1", "written": now(), "skyline": SKYLINE_PIN["version"],
                "exit": rc, "log": log, "lists": spec["targets"], "rows": counts,
                "srm_precursor_as_product_rows": srm_self})
    return 0 if rc == 0 and all(isinstance(v, int) and v > 0 for v in counts.values()) else 1


# ---- command line ----------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    pl = sub.add_parser("plan")
    for s in (pl,):
        s.add_argument("--session")
        s.add_argument("--search")
        s.add_argument("--de-dir")
        s.add_argument("--manifest")
        s.add_argument("--conditions")
        s.add_argument("--out")
        s.add_argument("--name")
        s.add_argument("--top", type=int)
        s.add_argument("--sif")
        s.add_argument("--python")
        s.add_argument("--ppm-ms1", type=float)
        s.add_argument("--ppm-ms2", type=float)
        s.add_argument("--ppm-source")
        s.add_argument("--max-parallel", type=int)
        s.add_argument("--partition")
        s.add_argument("--account")
        s.add_argument("--qos")
    sk = sub.add_parser("skip")
    sk.add_argument("--session")
    sk.add_argument("--out")
    sk.add_argument("--reason", required=True)
    stt = sub.add_parser("status")
    stt.add_argument("--session")
    stt.add_argument("--out")
    pr = sub.add_parser("prm")
    pr.add_argument("--session")
    pr.add_argument("--out")
    pr.add_argument("--dest")
    pr.add_argument("--instrument", required=True)
    pr.add_argument("--same-gradient", action="store_true")
    pr.add_argument("--run-length", type=float)
    pr.add_argument("--partition")
    pr.add_argument("--account")
    pr.add_argument("--qos")
    pr.add_argument("--no-staff-check", action="store_true", help=argparse.SUPPRESS)
    for name in ("build", "join", "prm-run"):
        x = sub.add_parser(name)
        x.add_argument("--work", required=True)
    im = sub.add_parser("import")
    im.add_argument("--work", required=True)
    im.add_argument("--task", required=True)
    a = ap.parse_args(argv)
    if a.cmd in ("plan", "skip", "status", "prm") and not (a.session or a.out or
                                                             getattr(a, "search", None)):
        ap.error("give --session (or --out)")
    if a.cmd == "plan" and a.top is not None and a.top < 1:
        ap.error("--top must be at least 1")
    return {"plan": cmd_plan, "skip": cmd_skip, "status": cmd_status, "prm": cmd_prm,
            "build": cmd_build, "import": cmd_import, "join": cmd_join,
            "prm-run": cmd_prm_run}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
