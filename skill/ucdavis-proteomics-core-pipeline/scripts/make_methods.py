#!/usr/bin/env python3
"""
make_methods.py  --  Generate a publication-ready LC-MS/MS Methods section from
facility raw data, plus the correct UC Davis Proteomics Core instrument-grant
acknowledgment.

It reads what it can directly from the raw metadata and fills the rest from facility
defaults that are CLEARLY TAGGED `[facility default — confirm]` so nothing is
silently fabricated (DE-LIMP rule #2). From a Bruker .d (bruker_method.py) it reads the LC
system and method, the ion source, TIMS settings, the dia-PASEF window scheme, cycle time and
collision-energy ramp, and writes them in the order published timsTOF Methods use; Thermo .raw
is identified by facility filename prefix. The analytical column comes, in order, from
--lc-column, from HyStar's ColumnInfo when an operator entered one, from --column-log (an
export of STAN's column-change log, matched to the acquisition dates), or else the facility's
standard column, tagged. It writes:

  methods.md          drop-in Methods prose (LC, MS, database search, sequence
                      database, differential expression) + a parameter table
                      (value + where each value came from) + an
                      instrument-specific Acknowledgments section
  methods_params.json the extracted parameters, machine-readable

Then to_docx.py can render methods.md to Word. The agent should verify the draft
against the extracted params and polish the prose (keep the acknowledgment exact).

Acknowledgments are from https://proteomics.ucdavis.edu/instrument-grant-acknowledgments
(verified 2026-06). Confirm exact wording there before publishing.

Usage:
  python3 make_methods.py --raw '/data/*.d' --out methods.md \
      [--lc-column "PepSep MAX C18, 10 cm × 150 µm, 1.5 µm"] \
      [--column-log stan_column_changes.csv [--column-log-instrument NAME]] \
      [--params wf/params.cfg --search-prov search/search_provenance.json \
       --workflow-manifest wf/workflow.manifest.json]   # adds the Database-search section
      [--de-dir output/tables]      # optional: adds a Differential-expression paragraph
      [--instrument "timsTOF HT" --acquisition DIA]   # used only when the raw files
                                                      # cannot be read from here
      [--submission <session dir>]  # its CoreOmics submission (submission_report.py):
                                    # adds Sample preparation -- who prepared the samples
      [--qc-bracket logs/qc_bracket.json]   # qc_bracket.py's record: adds its one fixed
                                    # Instrument performance sentence (never the verdict)
"""
import sys, os, re, csv, json, glob, argparse, statistics
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# analysis.tdf (and every other sqlite file in a .d) is opened read-only AND immutable -- see
# bruker_tdf.py for how a read-write open truncates a tdf (the state of 342 on HIVE).
import bruker_method  # noqa: E402
# the contaminant tag and the sidecar's keratin_sample reader: one definition each (rule 3)
from fetch_fasta import CONT_TAG, keratin_sample_recorded, FEW_ADDED_SEQUENCES  # noqa: E402

ACK_SOURCE = "https://proteomics.ucdavis.edu/instrument-grant-acknowledgments"
# (instrument-name substrings, facility filename prefixes, label, acknowledgment).
# Verified against the UC Davis Proteomics Core grant-acknowledgment page (2026-06). The page
# names only the timsTOF Pro 2 for the HHMI acknowledgment; Brett confirmed on 2026-09-25 that it
# covers the timsTOF HT too, so every timsTOF gets it.
# The filename prefixes are the facility's Thermo naming (FL*.raw, Ex*.raw) and are only a
# fallback for a .raw whose instrument no record names -- never for a .d (pick_ack).
ACKS = [
    (("fusion lumos", "lumos"), ("FL",), "Thermo Orbitrap Fusion Lumos",
     "Mass spectrometry was performed at the UC Davis Proteomics Core on an "
     "Orbitrap Fusion Lumos mass spectrometer acquired through NIH S10 grant "
     "S10OD021801."),
    (("exploris",), ("Ex",), "Thermo Orbitrap Exploris 480",
     "Mass spectrometry was performed at the UC Davis Proteomics Core on an "
     "Orbitrap Exploris 480 mass spectrometer acquired through NIH S10 grant "
     "S10OD026918-01A1."),
    (("timstof",), (), "Bruker timsTOF",
     "Mass spectrometry was performed at the UC Davis Proteomics Core on a Bruker "
     "timsTOF mass spectrometer. We thank Dr. Neil Hunter and the Howard Hughes "
     "Medical Institute for the timsTOF instrument."),
]
# THE vendor rule (2.14.2, S7/I09): the maker of an instrument from the words of its name. A
# name matching neither list, or both, has no vendor, and the Methods print VENDOR_TAG instead
# of guessing: a Bruker impact II's mzML was once called "(Thermo Fisher Scientific)" because
# every run that was not a .d took the Thermo sentence.
VENDOR_WORDS = {
    "Bruker": ("bruker", "timstof", "tims", "impact", "maxis", "compact"),
    "Thermo": ("thermo", "orbitrap", "exploris", "exactive", "lumos", "fusion", "eclipse",
               "astral", "ascend", "ltq", "velos"),
}
VENDOR_NAME = {"Bruker": "Bruker Daltonics", "Thermo": "Thermo Fisher Scientific"}
VENDOR_TAG = "[vendor — confirm]"
# Outside data (an mzML, an instrument the Core does not have) whose instrument name is in the
# Core's registry: the acknowledgment would say the Core ran it, which nothing here shows.
ACK_OUTSIDE_TAG = "[acquired at the UC Davis Proteomics Core? — confirm]"
# The instrument when it comes from the session's record rather than the raw file.
REC_SRC = "session record (workflow manifest) — not read from the raw file"


def instrument_vendor(name):
    """"Bruker", "Thermo" or None for an instrument name: case-insensitive, on words (a word may
    run straight into a model number, as in "Exploris480"). None when neither list matches, or
    both do -- the caller tags it, never guesses."""
    low = (name or "").lower()
    hits = {vendor for vendor, words in VENDOR_WORDS.items()
            if any(re.search(r"(?<![a-z])" + w + r"(?![a-z])", low) for w in words)}
    return hits.pop() if len(hits) == 1 else None


# The facility's standard column and emitter, as STAN's column catalogue names them
# (STAN config/columns.yml `default_column_id` and the default emitter, 2026-09-03). Printed
# only with the DEF tag: they say what is usually fitted, not what was fitted for a given run,
# and only for a Core timsTOF .d -- they are the timsTOF HT's (column_record's `core`).
LC_COLUMN_DEFAULT = ("PepSep MAX C18 column (10 cm × 150 µm i.d., 1.5 µm particles; "
                     "Bruker PepSep, part no. 1893483)")
EMITTER_DEFAULT = "a 20 µm-bore CaptiveSpray emitter (Bruker)"
# STAN config/columns.yml `defaults.oven_c` (oven_c_source: operator-reported, read off the
# timsTOF HT's Bruker Column Toaster on 2026-09-02). Not recorded per run, so always DEF-tagged.
COLUMN_TEMP_DEFAULT = "50 °C"
DEF = "[facility default — confirm]"
# A value this script could not find in any record of the run. Printed in place of the value,
# never replaced by a plausible default (DE-LIMP rule #2).
NR_TAG = "[not recorded — confirm]"
NOT_RECORDED = f"____ {NR_TAG}"
# A run left out of the DE whose reason exists but could not be read: run_de.R's runs_left_out
# record has `reasons_unreadable` (its decisions record was damaged), so "not recorded" would be
# wrong -- the user did record one (2.12 review).
REASON_NOT_READ = "reason not read: the decisions record could not be read [confirm]"
# Evosep loads every sample from an Evotip; the tip type and the amount on it are the Core's
# bench record, not the instrument's.
EVOTIP_TAG = "[Evotip type, loading protocol and peptide amount — confirm]"
# The ddaPASEF precursor-selection settings are in the .d's method, but this script does not
# read them: say so, rather than "not recorded".
DDA_TAG = "[not extracted by this script — take from the MS method; confirm]"
# A .d that could not be read: its values are unknown here, not unrecorded.
UNREADABLE_TAG = "[raw file not readable here — confirm]"

ENGINE_LABEL = {"diann": "DIA-NN", "sage": "Sage", "fragpipe": "FragPipe",
                "radiant": "Radiant", "alphadia": "AlphaDIA"}
# In-silico cleavage rule -> PSI-MS cleavage agent (name, accession). Accessions from OLS4,
# 2026-09-24. DIA-NN: `--cut K*,R*,!*P` is "canonical tryptic specificity" (DIA-NN README), so
# K*,R* without the !*P exception also cleaves before proline: Trypsin/P.
DIANN_CUT = {"K*,R*": ("Trypsin/P", "MS:1001313"), "K*,R*,!*P": ("Trypsin", "MS:1001251")}
# Sage's cleave_at + restrict and MSFragger's cut + nocut -> PSI-MS term, and how each engine's
# enzyme block is read: cleavage_rules.py, the one definition (run_search.py sizes a job from it)
from cleavage_rules import (SAGE_CUT, NONSPECIFIC, NO_CLEAVAGE, sage_cleavage,  # noqa: E402,F401
                            fragger_cleavage)
# Unimod record -> (title, monoisotopic delta mass), from unimod.xml (2026-09-24). Only these are
# named from a mass; anything else is reported as the engine gave it, flagged for the user.
UNIMOD = {1: ("Acetyl", 42.010565), 4: ("Carbamidomethyl", 57.021464),
          21: ("Phospho", 79.966331), 28: ("Gln->pyro-Glu", -17.026549)}
# ...and every variable modification and cysteine alkylation estimate_params.py can search, from
# its one table each (UniMod too)
from estimate_params import CYS_REAGENTS, VAR_MODS  # noqa: E402
UNIMOD.update({r["unimod"]: (r["mod"], r["mass"])
               for r in list(VAR_MODS.values()) + list(CYS_REAGENTS.values()) if r})


def _num(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def bruker_meta(d):
    """Acquisition parameters of a Bruker .d, each with the file and field it came from
    (bruker_method.read_run: analysis.tdf plus the HyStar/timsControl side files)."""
    if not os.path.exists(os.path.join(d, "analysis.tdf")):
        return None
    r = bruker_method.read_run(d)
    m = dict(r["values"], vendor="Bruker", file=os.path.basename(d.rstrip("/")),
             sources=r["sources"])
    m["software"] = " ".join(str(x) for x in (m.get("acquisition_software"),
                                              m.get("acquisition_software_version")) if x) or None
    s = bruker_method.window_scheme(m)
    if s:
        m["scheme"] = s
        m["n_windows"], m["n_window_groups"] = s["n_windows"], s.get("n_ramps")
        widths = [w["IsolationWidth"] for w in m["windows"] if w.get("IsolationWidth") is not None]
        if widths:
            m["isolation_width"] = round(statistics.median(widths), 1)
        if s.get("ce_range"):
            m["ce_low"], m["ce_high"] = (round(x, 1) for x in s["ce_range"])
    m["ce_check"] = bruker_method.check_ce_ramp(m)
    # the scheme and the check summarise the per-window lists; the JSON need not carry them
    m.pop("windows", None)
    m.pop("window_im", None)
    return m


def _wall(ts, like=None):
    """A timestamp as wall-clock time. AcquisitionDateTime carries the instrument PC's UTC
    offset; a STAN event_date carries none and is the lab's local time. An offset-aware event
    is moved to the acquisition's offset first, so the two are compared on one clock."""
    try:
        t = datetime.fromisoformat(str(ts).strip().replace("Z", "+00:00"))
    except ValueError:
        return None
    if t.tzinfo is not None and like is not None and like.tzinfo is not None:
        t = t.astimezone(like.tzinfo)
    return t.replace(tzinfo=None)


def read_column_log(path, instrument=None):
    """The column changes in a column log: an export (CSV, or JSON list) of STAN's
    maintenance_events rows -- event_type, event_date, column_vendor, column_model and
    optionally instrument, column_serial, first_run. STAN keeps these in PG Farm, which needs
    credentials, so this script never connects to it; the log has to be exported and passed.
    Returns (rows, problem): rows are the column_change events of one instrument."""
    try:
        with open(path, encoding="utf-8-sig", newline="") as fh:
            rows = (json.load(fh) if path.lower().endswith(".json") else
                    list(csv.DictReader(fh)))
    except (OSError, ValueError) as e:
        return [], f"the column log {os.path.basename(path)} could not be read ({e})"
    if isinstance(rows, dict):
        rows = rows.get("events") or rows.get("rows") or []
    rows = [r for r in rows if isinstance(r, dict)
            and str(r.get("event_type") or "").strip() == "column_change"]
    names = sorted({str(r.get("instrument")).strip() for r in rows if r.get("instrument")})
    if instrument:
        rows = [r for r in rows if str(r.get("instrument") or "").strip().lower()
                == instrument.strip().lower()]
        if not rows:
            return [], (f"the column log has no column change for instrument {instrument!r} "
                        f"(it names: {', '.join(names) or 'none'})")
    elif len(names) > 1:
        return [], (f"the column log covers {len(names)} instruments ({', '.join(names)}); "
                    f"pass --column-log-instrument to choose one")
    return rows, None


def column_record(user_column=None, runs=(), log_path=None, log_instrument=None,
                  first=None, last=None, core=True):
    """THE analytical-column statement for the Methods, with where it came from (DE-LIMP
    rule 3: the prose and the parameter table both read it). Most authoritative first:
      1. --lc-column, given by the user for this run;
      2. HyStar ColumnInfo, when an operator entered the column in the LC method;
      3. the column log's latest column_change at or before the first acquisition;
      4. the facility's standard column, tagged DEF -- what is usually fitted, not a record --
         only when `core` (a Core timsTOF .d: the default is the timsTOF HT's column); for any
         other run, a blank tagged NR_TAG (a Thermo cohort once got the Bruker PepSep column).
    `first`/`last` are the series' first and last AcquisitionDateTime."""
    warn = []
    if user_column:
        return {"text": user_column, "source": "--lc-column (user-given)", "tag": None,
                "warnings": warn}
    infos = sorted({m.get("column_info") for m in runs if m.get("column_info")})
    if len(infos) == 1:
        return {"text": infos[0], "source": "hystar.method ColumnInfo (entered in HyStar)",
                "tag": None, "warnings": warn}
    if len(infos) > 1:
        warn.append(f"the runs name {len(infos)} different columns in HyStar ColumnInfo "
                    f"({'; '.join(infos)})")
    if log_path:
        rows, problem = read_column_log(log_path, log_instrument)
        t0 = _wall(first) if first else None
        t1 = _wall(last) if last else None
        like = datetime.fromisoformat(first) if first else None
        if problem:
            warn.append(problem)
        elif t0 is None:
            warn.append("no acquisition time could be read from the raw files, so the column "
                        "log could not be matched to them")
        else:
            dated = [(w, r) for r in rows if (w := _wall(r.get("event_date"), like)) is not None]
            before = [x for x in dated if x[0] <= t0]
            during = sorted(x[0] for x in dated if t1 is not None and t0 < x[0] <= t1)
            if during:
                warn.append("the column log records a column change during this series ("
                            + ", ".join(w.isoformat(sep=" ") for w in during) + ")")
            if not before:
                warn.append(f"the column log has no column change before the first run "
                            f"({t0.isoformat(sep=' ')})")
            else:
                when, row = max(before, key=lambda x: x[0])
                vendor = str(row.get("column_vendor") or "").strip()
                model = str(row.get("column_model") or "").strip()
                serial = str(row.get("column_serial") or "").strip()
                if not model:
                    warn.append(f"the column change of {when.date()} in the column log does "
                                f"not name the column")
                elif not during:
                    paren = [x for x in (vendor if vendor and not model.lower().startswith(
                        vendor.lower()) else None, f"serial/LOT {serial}" if serial else None)
                        if x]
                    anchor = str(row.get("first_run") or "").strip()
                    return {"text": model + (f" ({'; '.join(paren)})" if paren else ""),
                            "source": (f"column log {os.path.basename(log_path)}: column_change "
                                       f"of {when.isoformat(sep=' ')}"
                                       + (f", first run on it {anchor}" if anchor else "")),
                            "tag": None, "warnings": warn, "installed": when.isoformat()}
    if not core:
        return {"text": "____", "source": "not recorded: the facility's standard column is a "
                "timsTOF one and is not given for this instrument — confirm", "tag": NR_TAG,
                "warnings": warn}
    return {"text": LC_COLUMN_DEFAULT, "source": "facility default — confirm", "tag": DEF,
            "warnings": warn}


def _r(x, nd=2):
    """A recorded number for prose: 85.05 -> '85', 0.7 -> '0.70' (nd decimals, or an integer
    when it is one to within 0.1). The parameter table keeps the value as recorded."""
    if x is None:
        return None
    return f"{x:.0f}" if abs(x - round(x)) < 0.1 else f"{x:.{nd}f}"


def _unit(u):
    return {"l/min": "L/min"}.get(u, u)


def _with_tag(text, tag):
    return f"{text} {tag}" if tag else text


def _source_name(name):
    """The ion source as Methods write it: timsTOF files name code 11 "Captive Spray", Bruker's
    product is "CaptiveSpray" -- one spelling in the prose (the table keeps the file's)."""
    return re.sub(r"(?i)\bcaptive\s+spray\b", "CaptiveSpray", name or "")


def lc_paragraph(rep, col, is_bruker, lc_known, facility=True):
    """The Liquid chromatography paragraph: LC system and method as the .d recorded them,
    the column from column_record(), and tagged placeholders for what no file records.
    `is_bruker`: a Core timsTOF .d, whose CaptiveSpray emitter and column temperature are the
    facility's (tagged) defaults. `facility`: a run on a Core instrument (.d or .raw); for data
    from elsewhere (an mzML, an instrument the Core does not have) no Core LC default is
    printed, only tagged blanks (2.14.2, I09)."""
    column = (f"an analytical column {NOT_RECORDED}" if col.get("tag") == NR_TAG else
              _with_tag(f"a {col['text']}", col.get("tag")))
    phases = ("Mobile phase A was 0.1% (v/v) formic acid in water and mobile phase B 0.1% "
              f"(v/v) formic acid in acetonitrile {DEF}." if facility else
              f"The mobile phases were ____ {NR_TAG}.")
    if not lc_known and not facility:
        return (f"Peptides were separated by liquid chromatography on {column}. The mobile "
                f"phases were ____ {NR_TAG}, and the column was interfaced to the mass "
                f"spectrometer by ____ {NR_TAG}. The LC system and gradient were "
                f"[LC system / gradient — confirm] {NR_TAG}.")
    if not lc_known:
        return (f"Peptides were separated by reversed-phase nano-LC on {column}, using water "
                "containing 0.1% (v/v) formic acid as mobile phase A and acetonitrile containing "
                f"0.1% (v/v) formic acid as mobile phase B {DEF}. "
                + ("The column was interfaced to the mass spectrometer through a Bruker "
                   f"CaptiveSpray source with {EMITTER_DEFAULT} {DEF}. " if is_bruker else
                   f"The column was interfaced to the mass spectrometer by a nanospray source "
                   f"{DEF}. ")
                + f"The LC system and gradient were [LC system / gradient — confirm] {NR_TAG}.")
    system = rep["lc_system"] + (f" LC system ({rep['lc_vendor']})" if rep.get("lc_vendor")
                                 else " LC system")
    meth = rep.get("lc_method") or rep.get("lc_method_name")
    run = f" (run time {_r(rep['lc_run_min'], 1)} min)" if rep.get("lc_run_min") else ""
    evosep = "evosep" in rep["lc_system"].lower()
    # The published order (literature survey, 2026-09-24): LC coupled to the timsTOF via its
    # source; then the column and its temperature; then the LC method; then mobile phases.
    s = (f"Peptides were loaded onto Evotips {EVOTIP_TAG} and analysed"
         if evosep and "evotip" in (rep.get("tray_type") or "").lower() else
         "Peptides were analysed")
    s += f" on {'an' if system[0] in 'AEIOU' else 'a'} {system} coupled online to a " \
         f"{rep.get('instrument') or NOT_RECORDED} mass spectrometer (Bruker Daltonics)"
    s += (f" via a {_source_name(rep['source_type'])} ion source."
          if rep.get("source_type") and is_bruker else ".")
    temp = f"{COLUMN_TEMP_DEFAULT} {DEF}" if is_bruker else f"____ °C {NR_TAG}"
    s += f" Peptides were separated on {column}, at a column temperature of {temp},"
    spd = re.search(r"(\d+)\s*samples?\s*per\s*day", meth or "", re.I)
    s += (f" with the {meth} ({spd.group(1)} SPD) method{run}." if spd else
          f" with the '{meth}' method{run}." if meth else
          f" with the LC method ____ {NR_TAG}.")
    if not evosep:
        s += (f" The gradient (time, %B, flow) was ____ {NR_TAG}: this LC method's gradient "
              "table is not read from the raw file.")
    return s + " " + phases


def ms_paragraph(rep, v, coupled=False):
    """The Mass spectrometry paragraph for a timsTOF run, in the order published dia-PASEF
    Methods report it (literature survey, 2026-09-24): acquisition mode and scan range, TIMS,
    window scheme, cycle time, collision energy, then source settings and software. `coupled`:
    the LC paragraph has already named the instrument and its source. One unit form throughout,
    "1/K₀ 0.70–1.30 V·s/cm²"; a cycle time is never called a duty cycle. A sentence whose values
    the files lack is left out, or carries the blank tag `v()` gives it; nothing is filled from
    memory."""
    out = []
    pol = f"{rep['polarity']}-ion " if rep.get("polarity") else ""
    meth = f" (method {rep['ms_method']})" if rep.get("ms_method") else ""
    out.append(f"The mass spectrometer was operated in {pol}{v(rep.get('mode'))} mode{meth}."
               if coupled else
               f"Mass spectra were acquired on a {v(rep.get('instrument'))} mass spectrometer "
               f"(Bruker Daltonics) operated in {pol}{v(rep.get('mode'))} mode{meth}.")
    lo, hi = rep.get("mz_low"), rep.get("mz_high")
    scan = "MS1 and MS2 spectra were" if rep.get("mode") in ("dia-PASEF", "ddaPASEF") else \
        "Spectra were"
    out.append(f"{scan} recorded over m/z {v(_r(lo, 0) if lo is not None else None)}–"
               f"{v(_r(hi, 0) if hi is not None else None)}.")
    ramp, acc = rep.get("ramp_ms"), rep.get("accumulation_ms")
    tims = (f"The trapped ion mobility (TIMS) ramp spanned 1/K₀ "
            f"{v(_r(rep.get('im_low')))}–{v(_r(rep.get('im_high')))} V·s/cm²")
    if ramp and acc:
        duty = 100.0 * acc / ramp
        tims += (f", with ramp and accumulation times of {_r(ramp)} ms each"
                 if abs(ramp - acc) < 0.01 else
                 f", with a ramp time of {_r(ramp)} ms and an accumulation time of {_r(acc)} ms")
        tims += f" ({duty:.0f}% duty cycle)"
    out.append(tims + ".")
    s = rep.get("scheme") or {}
    if rep.get("mode") == "dia-PASEF" and s:
        n_r, fpc = s.get("n_ramps"), rep.get("frames_per_cycle")
        if n_r and fpc == n_r + 1:
            cyc = f"Each acquisition cycle comprised one MS1 frame and {n_r} dia-PASEF frames"
        elif fpc:
            cyc = f"Each acquisition cycle comprised {fpc} frames"
        else:
            cyc = "The dia-PASEF method used"
        w = s.get("width")
        wtxt = (f" of {_r(w, 1)} Th" if isinstance(w, float) else
                f" of {_r(w[0], 1)}–{_r(w[1], 1)} Th" if w else "")
        detail = []
        if s.get("spacing") is not None and s.get("overlap") is not None:
            ov = s["overlap"]
            detail.append(f"{_r(s['spacing'], 1)} Th spacing, " +
                          (f"{_r(ov, 1)} Th overlap" if ov > 0 else
                           f"{_r(-ov, 1)} Th gap" if ov < 0 else "no overlap"))
        if s.get("per_ramp"):
            a, b = s["per_ramp"]
            detail.append(f"{a if a == b else f'{a}–{b}'} per TIMS ramp")
        place = (" placing" if cyc.startswith("Each") else "") + \
            f" {s['n_windows']} isolation windows{wtxt}" + \
            (f" ({'; '.join(detail)})" if detail else "")
        area = []
        if s.get("mz_lo") is not None:
            area.append(f"m/z {_r(s['mz_lo'], 1)}–{_r(s['mz_hi'], 1)}")
        if s.get("im_lo") is not None:
            area.append(f"1/K₀ {_r(s['im_lo'])}–{_r(s['im_hi'])} V·s/cm²")
        out.append(cyc + ("," if cyc.startswith("Each") else "") + place
                   + (f" across {' and '.join(area)}" if area else "") + ".")
    if rep.get("mode") == "ddaPASEF":
        out.append("Precursor selection (PASEF ramps per cycle, target intensity, charge and "
                   f"mobility filters, dynamic exclusion) was ____ {DDA_TAG}.")
    if rep.get("cycle_s"):
        out.append(f"The cycle time was {rep['cycle_s']:.2f} s.")
    ramp_ce, chk = rep.get("ce_ramp"), rep.get("ce_check") or {}
    if ramp_ce and chk.get("status") != "mismatch":
        pts = ramp_ce["points"]
        if len(pts) == 2:
            (x0, y0), (x1, y1) = pts
            out.append(f"The collision energy was ramped linearly with ion mobility from "
                       f"{_r(y0, 1)} eV at 1/K₀ {_r(x0)} V·s/cm² to {_r(y1, 1)} eV at "
                       f"1/K₀ {_r(x1)} V·s/cm².")
        else:
            out.append("The collision energy followed a mobility-dependent ramp through "
                       + ", ".join(f"{_r(y, 1)} eV at 1/K₀ {_r(x)} V·s/cm²" for x, y in pts)
                       + ".")
    elif rep.get("ce_low") is not None:
        out.append(f"Collision energies of {rep['ce_low']}–{rep['ce_high']} eV were applied "
                   "across the isolation windows"
                   + (" [they do not match the method's recorded ramp — confirm]"
                      if chk.get("status") == "mismatch" else "") + ".")
    if rep.get("source_type"):
        bits = [f"{label} {_r(rep[key]['value'], 1)} {_unit(rep[key]['unit'])}".strip()
                for key, label in (("capillary_v", "a capillary voltage of"),
                                   ("dry_gas", "a dry gas flow of"),
                                   ("dry_temp", "a dry temperature of")) if rep.get(key)]
        joined = (", ".join(bits[:-1]) + " and " + bits[-1]) if len(bits) > 1 else \
            (bits[0] if bits else "")
        out.append(f"The {_source_name(rep['source_type'])} source was fitted with "
                   f"{EMITTER_DEFAULT} {DEF}"
                   + (f" and operated at {joined}" if joined else "") + ".")
    ctrl = rep.get("ms_control") or rep.get("control_software")
    if ctrl or rep.get("acquisition_software_version"):
        sw = f"Data were acquired with {ctrl or rep.get('acquisition_software')}"
        if rep.get("acquisition_software_version"):
            sw += f" (acquisition software version {rep['acquisition_software_version']})"
        if rep.get("hystar_version"):
            sw += f" and HyStar {rep['hystar_version']}"
        out.append(sw + ".")
    elif rep.get("software"):
        out.append(f"Data were acquired with {rep['software']}.")
    return " ".join(out)


def thermo_meta(f):
    """Thermo .raw: the model is not readable here without a vendor reader. The facility
    filename prefix (FL*, Ex*) is kept as a GUESS, apart from `instrument`: a session record
    (--instrument, read from the file by detect_acquisition.py) outranks it, and a .raw from
    elsewhere named FLAG_... must not become a Fusion Lumos run."""
    base = os.path.basename(f)
    return {"vendor": "Thermo", "file": base, "mode": None,
            "prefix_instrument": prefix_instrument([f])}


def prefix_instrument(files):
    """The instrument the facility's Thermo filename prefix implies, when EVERY file is a .raw
    carrying the same entry's prefix; otherwise None. Never for a .d: a timsTOF run renamed
    FLAG_IP_1.d or Exp3_HeLa.d is not a Fusion Lumos or an Exploris run."""
    raws = [os.path.basename(f.rstrip("/")) for f in files]
    if not raws or not all(b.lower().endswith(".raw") for b in raws):
        return None
    for subs, prefixes, label, _ in ACKS:
        if prefixes and all(any(b.startswith(p) for p in prefixes) for b in raws):
            return label
    return None


def detect(files, instrument=None):
    """One metadata dict per readable run. A .d is read (bruker_meta), a .raw is Thermo; any other
    file (an mzML from elsewhere) records nothing of its maker, so its vendor is
    instrument_vendor() of the recorded `instrument` -- None, and tagged, when that names
    neither."""
    metas = []
    for f in files:
        low = f.lower().rstrip("/")
        if low.endswith(".d"):
            mm = bruker_meta(f)
        elif low.endswith(".raw"):
            mm = thermo_meta(f)
        else:
            mm = {"vendor": instrument_vendor(instrument), "file": os.path.basename(f),
                  "instrument": instrument or None}
            if instrument:
                mm["instrument_source"] = REC_SRC
        if mm: metas.append(mm)
    return metas


def _session_acquisition_json(a):
    """input/acquisition.json (step 2's detect_acquisition.py output) of the session this run
    belongs to: the first folder that is a session (session.is_session_dir) and holds one,
    walking up from --out, --workflow-manifest, --de-dir, --params, --search-prov, then each
    --raw. None when there is none."""
    from session import is_session_dir, paths_for
    starts = [getattr(a, k, None) for k in ("out", "workflow_manifest", "de_dir", "params",
                                             "search_prov")] + list(getattr(a, "raw", None) or [])
    seen = set()
    for start in starts:
        if not start:
            continue
        d = os.path.dirname(os.path.abspath(start))
        for _ in range(4):
            if d in seen:
                break
            seen.add(d)
            if is_session_dir(d):
                cand = os.path.join(paths_for(d)["input_dir"], "acquisition.json")
                if os.path.isfile(cand):
                    return cand
            up = os.path.dirname(d)
            if up == d:
                break
            d = up
    return None


def acquisition_from_records(a):
    """THE acquisition mode (DIA/DDA) for the Methods, with where it came from, in this order:
      1. --acquisition (as finalize passes it from the session record);
      2. the workflow manifest's `acquisition` (resolve_defaults.py records step 2's answer);
      3. the session's input/acquisition.json (`overall`, step 2's detection);
      4. else (None, None): the MS sentence then carries "[DDA/DIA — confirm]".
    A record that says neither DIA nor DDA ("mixed", "unknown") is passed over. The step 9d run
    once wrote "[DDA/DIA — confirm]" with step 2's answer in the manifest beside it (S7)."""
    flag = (getattr(a, "acquisition", None) or "").strip().upper()
    if flag:
        return flag, "--acquisition (the session record)"
    wf = getattr(a, "workflow_manifest", None)
    rec = _load_json(wf) if wf else None
    val = str(rec.get("acquisition") or "").strip().upper() if isinstance(rec, dict) else ""
    if val in ("DIA", "DDA"):
        return val, f"workflow manifest ({os.path.basename(wf)}): step 2 detection"
    path = _session_acquisition_json(a)
    rec = _load_json(path) if path else None
    if isinstance(rec, dict):
        val = str(rec.get("overall") or rec.get("acquisition") or "").strip().upper()
        if val in ("DIA", "DDA"):
            return val, "input/acquisition.json: step 2 detection"
    return None, None


def pick_ack(instrument, files):
    """The acknowledgment for the instrument NAME, matched across every registry entry first. The
    Thermo filename prefix is only a fallback for .raw files whose instrument nothing names -- a
    real timsTOF HT run renamed FLAG_IP_1.d once got the Fusion Lumos S10 grant."""
    missing = (None, f"[Instrument not in the UC Davis acknowledgment registry — check "
                     f"{ACK_SOURCE} and insert the correct instrument-grant acknowledgment.]")
    instr = (instrument or "").lower()
    if instr:
        for subs, _prefixes, label, text in ACKS:
            if any(s in instr for s in subs):
                return label, text
        return missing               # a named instrument the registry lacks: no filename guess
    guess = prefix_instrument(files)
    for _subs, _prefixes, label, text in ACKS:
        if label == guess:
            return label, text
    return missing


def _load_json(path):
    if not path or not os.path.isfile(path):
        return None
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def _unimod_for_mass(mass):
    for rid, (_, m) in UNIMOD.items():
        if mass is not None and abs(m - mass) <= 0.002:
            return rid
    return None


def _mod(unimod, name, mtype, position, targets, mass, source):
    if unimod in UNIMOD:
        name = UNIMOD[unimod][0]
    return {"name": name, "unimod": unimod, "type": mtype, "position": position,
            "targets": targets, "mass": mass, "source": source}


def _diann_mod(vals, mtype, where):
    """`--var-mod/--fixed-mod name,mass,sites[,label]` (DIA-NN README): sites are residues, `n`
    for the peptide N-terminus and `*n` for the protein N-terminus."""
    parts = ",".join(vals).split(",")
    name = parts[0].strip() if parts else ""
    try:
        mass = float(parts[1])
    except (IndexError, ValueError):
        mass = None
    sites = parts[2].strip() if len(parts) > 2 else ""
    rid = None
    if name.lower().startswith("unimod:") and name[7:].isdigit():
        rid = int(name[7:])
    position = ("Protein N-term" if "*n" in sites else
                "Any N-term" if "n" in sites else "Anywhere")
    targets = "".join(c for c in sites if c.isalpha() and c.isupper())
    m = _mod(rid, name, mtype, position, targets, mass, where)
    m["label"] = len(parts) > 3 and parts[3].strip() == "label"
    return m


def _sage_mods(mods, mtype, where):
    """Sage static_mods/variable_mods: key = residue, or a terminus symbol (Sage DOCS.md v0.14.7:
    `^` peptide N-term, `$` peptide C-term, `[` protein N-term, `]` protein C-term), optionally
    followed by a residue ("^E")."""
    out = []
    for key, v in (mods or {}).items():
        masses = v if isinstance(v, list) else [v]
        pos = {"^": "Any N-term", "$": "Any C-term", "[": "Protein N-term",
               "]": "Protein C-term"}.get(key[:1], "Anywhere")
        targets = key[1:] if pos != "Anywhere" else key
        for mass in masses:
            try:
                mass = float(mass)
            except (TypeError, ValueError):
                mass = None
            rid = _unimod_for_mass(mass)
            out.append(_mod(rid, f"{mass:+.4f} Da" if mass is not None else "?", mtype, pos,
                            targets, mass, where))
    return out


def _sage_results_json(search_prov):
    """Sage's own record of the search as it ran -- <out>/results.json, beside
    search_provenance.json (sage-cli input.rs `Search`: "Actual search parameters") -- or None.
    It wins over the config given (N5, review F16)."""
    if not search_prov:
        return None
    ran = _load_json(os.path.join(os.path.dirname(os.path.abspath(search_prov)), "results.json"))
    return ran if isinstance(ran, dict) else None


def _sage_lfq_fdr(sa):
    """The FDR sentence of a quantifying Sage search, from its sage_adapt record: the sentence it
    recorded, else -- a record from before 2.14.2 -- sage_lfq_check.methods_sentence() of the
    record's own q_value_max. None for a refused search (no report) or no record."""
    if not isinstance(sa, dict) or not sa.get("report"):
        return None
    if sa.get("methods_sentence"):
        return sa["methods_sentence"]
    q = sa.get("q_value_max")
    if isinstance(q, (int, float)) and not isinstance(q, bool):
        import sage_lfq_check
        return sage_lfq_check.methods_sentence(q)
    return None


def _fragpipe_record(rec, prov, pf, search_prov):
    """search_record()'s FragPipe part: the MSFragger parameters FragPipe wrote for this search
    (fragger.params) and the Philosopher filter of the workflow it saved (fragpipe.workflow),
    both in the search folder beside search_provenance.json (fragpipe_params.py, the one
    reader). What cannot be read is said so, never filled (N4)."""
    import fragpipe_params as fpp
    out = os.path.dirname(os.path.abspath(search_prov)) if search_prov else None
    given = pf if os.path.basename(pf) == fpp.FRAGGER_PARAMS else None
    params_path = (fpp.find(out, fpp.FRAGGER_PARAMS) or given
                   or fpp.find(os.path.dirname(os.path.abspath(pf)), fpp.FRAGGER_PARAMS))
    params = None
    if params_path:
        try:
            params = fpp.read_fragger_params(params_path)
        except OSError as e:
            rec["warnings"].append(f"{params_path} could not be read ({e.strerror or e})")
    if params is None:
        rec["warnings"].append(f"{fpp.FRAGGER_PARAMS}, the MSFragger parameters FragPipe writes "
                               f"into the search folder, was not found here, so the enzyme, "
                               f"modifications and mass tolerances are not stated")
    else:
        w = os.path.basename(params_path)
        dg = fpp.digest(params)
        rec["cleavage"] = fragger_cleavage(dg, w)
        if dg["specificity"] in ("semi", "full") and dg["enzymes"] \
                and dg["enzymes"][0]["missed"] is not None:
            rec["missed_cleavages"] = {"value": int(dg["enzymes"][0]["missed"]), "source": w}
        lo, hi = (_num(params.get(k)) for k in ("digest_min_length", "digest_max_length"))
        if lo is not None or hi is not None:
            rec["pep_len"] = {"value": (int(lo) if lo is not None else None,
                                        int(hi) if hi is not None else None), "source": w}
        for mtype, mass, pos, targets, _per in fpp.mods(params):
            rec["mods"].append(_mod(_unimod_for_mass(mass), f"{mass:+.4f} Da", mtype, pos,
                                    targets, mass, w))
        rec["mods_read"] = True
        mv = _num(params.get("max_variable_mods_per_peptide"))
        if mv is not None:
            rec["max_var_mods"] = {"value": int(mv), "source": w}
        rec["met_excision"] = params.get("clip_nTerm_M") == "1"
        # the database MSFragger read (fragger.params database_name), with the record's sha256
        # when it is the one run_search.py recorded (H7 bug 1)
        dbn = (params.get("database_name") or "").strip()
        pdb = prov.get("fragpipe_database") if isinstance(prov.get("fragpipe_database"),
                                                          dict) else {}
        if dbn:
            same = pdb.get("path") and os.path.basename(pdb["path"]) == os.path.basename(dbn)
            rec["database"] = {"name": os.path.basename(dbn), "source": f"{w} database_name",
                               "sha256": pdb.get("sha256") if same else None,
                               "n_entries": pdb.get("n_entries") if same else None,
                               "n_decoys": pdb.get("n_decoys") if same else None,
                               "decoy_tag": pdb.get("decoy_tag") if same else None}
        for key, which in (("ms1_tol", "precursor"), ("ms2_tol", "fragment")):
            t = fpp.tolerance(params, which)
            if t:
                rec[key] = {"value": t["value"], "unit": t["unit"], "symmetric": t["symmetric"],
                            "source": w}
        # mass offsets, labile mode ... -- never described as a plain closed search (review F3)
        rec["fragger_extras"] = fpp.search_extras(params)
        cal = _num(params.get("calibrate_mass"))
        if cal:
            rec["warnings"].append(f"MSFragger mass calibration was on (calibrate_mass="
                                   f"{int(cal)} in {w}), so the tolerances it searched with may "
                                   f"differ from these starting values; its log has the ones used")
    saved = fpp.find(out, fpp.SAVED_WORKFLOW)
    wf_path = saved or (pf if pf.endswith(".workflow") and os.path.isfile(pf) else None)
    if wf_path:
        try:
            props, _ = fpp.read_workflow(wf_path)
        except OSError as e:
            props = {}
            rec["warnings"].append(f"{wf_path} could not be read ({e.strerror or e})")
        filt = props.get(fpp.PHILOSOPHER_FILTER_KEY)
        # a workflow that ran no quantification: the search was for identifications only
        # (fragpipe_params.workflow_quantification, the one reading; Ruling 73a)
        quant = fpp.workflow_quantification(props)
        if quant["quantified"] is False and not rec.get("route_from_record"):
            rec["route_sentence"] = ("The search was for identifications only, with no "
                                     f"quantification: {fpp.methods_why(quant['why'])}.")
        if filt is not None and fpp.report_runs(props):       # the one reading (FU2)
            # the filter AS RUN, from FragPipe's own log beside the record, else the workflow's
            # line, said so; and over which runs its rates hold (review F9) --
            # fragpipe_params.fdr_statement, the one statement, the DE's too (C2 proteomics M1)
            st = fpp.fdr_statement(out, filt)
            if st:
                rec["fdr_sentence"] = st["sentence"]
                rec["fdr_source"] = ("FragPipe's log (the filter as run)" if st["from_log"] else
                                     f"{os.path.basename(wf_path)} "
                                     + ("(the workflow FragPipe saved with the search)" if saved
                                        else "(the workflow given to the search)"))


# DIA-NN echoes its argv, joined with single spaces and unquoted, near the top of its log; an
# option's value runs to the next "--" (fran_deposit._FASTA_ARG has the why: FragPipe's DIA-NN
# 1.8.1 logs even `--cfg <workdir>/filelist_diann.txt-- `). The command line is the first line of
# the head that carries an option. A `--cfg <file>` on it is logged UNEXPANDED (DIA-NN 1.8.1,
# 2.6.1, 2.7.0 on HIVE; run_search.py --one-step passes the whole cfg that way), so its flags are
# read from the file, spliced in where --cfg stands, as DIA-NN reads them.
_DIANN_CMD_LINE = re.compile(r"^\S.*?\s--[a-z]", re.M)


def _logged_options(line):
    """[(flag, value)] of a logged DIA-NN command line, in order; the binary is dropped."""
    out = []
    for seg in line.split("--")[1:]:
        flag, _, val = seg.strip().partition(" ")
        if flag:
            out.append(("--" + flag, val.strip().strip("\"'")))
    return out


def _diann_workdirs(logp, opts):
    """Where a relative path on the logged command line resolves: DIA-NN's working directory,
    which the log does not name. A relative --out tells it exactly (the log is written beside
    --out, so the log's folder minus --out's own folder part; FragPipe's DIA route runs in its
    workdir with --out dia-quant-output/report.tsv); otherwise the log's folder (run_search.py
    runs every DIA-NN job with `cd <out>`), then its parent -- CANDIDATES, of which
    _expand_cfgs() reads one only when the others hold no different file of that name. (The
    skill's own jobs pass --cfg as an absolute path, run_search.run_diann.)"""
    d = os.path.dirname(os.path.abspath(logp))
    outs = [v for f, v in opts if f == "--out" and v]
    if outs and not os.path.isabs(outs[-1]):
        rel = os.path.dirname(os.path.normpath(outs[-1]))
        if not rel:
            return [d]
        if d.endswith(os.sep + rel):
            return [d[:-len(rel) - 1]]
    return [d, os.path.dirname(d)]


def _same_file_contents(paths):
    """True when every path holds the same bytes (one file seen from two folders counts)."""
    blobs = set()
    for p in paths:
        with open(p, "rb") as fh:
            blobs.add(fh.read())
    return len(blobs) <= 1


def _expand_cfgs(opts, workdirs, depth=0):
    """opts with each --cfg replaced by the flags of the file it names (cfg_groups, the one
    tokeniser) -> (expanded [(flag, [values])], [cfg files read], [why each cfg NOT read]).

    A relative --cfg is looked for in every candidate working directory (_diann_workdirs). Found
    in more than one, with different contents, it is NOT read: the log does not say which one
    DIA-NN opened, and reading the first would describe a run from a file it may never have seen
    (sage-review, 6d49ec3: the log's folder won even when DIA-NN ran in its parent)."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from diann_parallel import cfg_tokens, cfg_groups, CfgError
    out, read, unread = [], [], []
    for flag, val in opts:
        if flag != "--cfg":
            out.append((flag, val.split() if isinstance(val, str) else list(val)))
            continue
        path = val if isinstance(val, str) else " ".join(val)
        cands = [path] if os.path.isabs(path) else [os.path.join(w, path) for w in workdirs]
        hits = [c for c in cands if os.path.isfile(c)]
        try:
            if len(hits) > 1 and not _same_file_contents(hits):
                raise CfgError(f"the log names --cfg {path}, a relative path, and both "
                               f"{' and '.join(hits)} exist and differ: which one DIA-NN read "
                               "is not recorded")
            if not hits:
                raise CfgError(f"the log names --cfg {path}, which could not be read (not found"
                               + ("" if os.path.isabs(path) else f" in {', '.join(workdirs)}")
                               + ")")
            if depth > 3:
                raise CfgError(f"the log names --cfg {path}, nested more than 3 deep")
            sub, r, u = _expand_cfgs(cfg_groups(cfg_tokens(hits[0])), workdirs, depth + 1)
        except CfgError as e:
            unread.append(str(e))
            continue
        out += sub
        read += [hits[0]] + r
        unread += u
    return out, read, unread


def diann_cont_quant_exclude(report=None, prov=None, params_file=None):
    """Did DIA-NN run with --cont-quant-exclude for this report? -> the record, or None when
    nothing answered (diann_cont_quant_exclude_why() says why)."""
    return diann_cont_quant_exclude_why(report, prov, params_file)[0]


def diann_cont_quant_exclude_why(report=None, prov=None, params_file=None):
    """Did DIA-NN run with --cont-quant-exclude for this report, and where does that come from?

    THE reader of that flag (DE-LIMP rule 3): search_record() takes it from here for the
    database-search paragraph, record_run.py for the run record, and contaminants.R asks it
    (`make_methods.py cont-quant-exclude <report>`) before a DE descriptor says what DIA-NN left
    out of its quantities. Why (sage-review N1, 2026-09-30): the dpc and maxlfq descriptors said
    "DIA-NN --cont-quant-exclude" left keratin peptides out for every DIA-NN report -- FragPipe's
    DIA-NN step and a hand-run DIA-NN included, which nothing had checked.

    Most authoritative first: the command line in the DIA-NN log that wrote the report
    (<report>.log.txt: what ran), with every `--cfg` file it names read in place; then the
    parameters file the search ran with (`params_file`, else search_provenance.json
    resolved_params_file, then params_file -- only for a DIA-NN search: a Sage or FragPipe
    parameters file says nothing about DIA-NN). A command line whose --cfg cannot be read is not
    taken as "absent" (sage-review, ca9aff2: every --one-step search then read "not set").
    -> ({"value": tag or None, "source": where}, None) -- value None when that source was read
    and has no such flag -- or (None, why) when no source could be read, which the caller prints
    as NOT RECORDED with that reason (rule 2): "no DIA-NN log" was said even when the log named a
    --cfg that could not be read (sage-review, 6d49ec3)."""
    why = []
    if prov is None and report:
        # the report's search folder: its own, or its parent for FragPipe's DIA route
        from session import search_folder_of, SEARCH_RECORD     # the one rule (T44 F1)
        prov = _load_json(os.path.join(search_folder_of(report), SEARCH_RECORD))
    prov = prov if isinstance(prov, dict) else {}
    if report:
        logp = os.path.splitext(report)[0] + ".log.txt"
        try:
            with open(logp, errors="replace") as fh:
                head = fh.read(2 << 20).replace("\r", "")
        except OSError:
            head = ""
        cmd = _DIANN_CMD_LINE.search(head)
        if cmd:
            opts = _logged_options(head[cmd.start():].split("\n", 1)[0])
            groups, read, unread = _expand_cfgs(opts, _diann_workdirs(logp, opts))
            vals = [v for f, v in groups if f == "--cont-quant-exclude"]
            src = f"the DIA-NN command line in {os.path.basename(logp)}" + (
                f" and the --cfg file{'s' if len(read) > 1 else ''} it names "
                f"({', '.join(os.path.basename(r) for r in read)})" if read else "")
            if vals:
                return {"value": (vals[-1][0].strip("\"'") if vals[-1] else None) or None,
                        "source": src}, None
            if not unread:
                return {"value": None, "source": src}, None
            # the flag may be in the --cfg that could not be read: not "absent" -- ask the
            # parameters file the search ran with, else NOT RECORDED
            why += [f"{os.path.basename(logp)}: {u}" for u in unread]
        else:
            why.append(f"no DIA-NN log beside the report ({os.path.basename(logp)})" if not head
                       else f"no DIA-NN command line in {os.path.basename(logp)}")
    else:
        why.append("no report was named")
    cands = [params_file] if params_file else (
        [prov.get("resolved_params_file"), prov.get("params_file")]
        if (prov.get("engine") or "diann").lower() == "diann" else [])
    for pf in cands:
        if not (pf and os.path.isfile(pf)):
            continue
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from diann_parallel import cfg_tokens, cfg_groups, CfgError
        try:
            groups = cfg_groups(cfg_tokens(pf))
        except CfgError:
            continue
        vals = [v for f, v in groups if f == "--cont-quant-exclude"]
        return {"value": (vals[-1][0] if vals[-1] else None) if vals else None,
                "source": os.path.basename(pf)}, None
    why.append("no DIA-NN parameters file was found" if (prov.get("engine") or "diann").lower()
               == "diann" else f"the search's parameters file is {prov.get('engine')}'s, not "
               "DIA-NN's")
    return None, "; ".join(why)


def _cys_alkylation(paths):
    """estimate_params.py's `cys_alkylation` record (what was searched for cysteine, and on whose
    word) from the first of these params files whose .rationale.json sidecar holds one, else None.
    The resolved cfg step 1b writes has no sidecar, so the given params file is asked too."""
    for path in paths:
        side = _load_json(f"{path}.rationale.json") if path else None
        if isinstance(side, dict) and isinstance(side.get("cys_alkylation"), dict):
            return side["cys_alkylation"]
    return None


def _labels(paths):
    """estimate_params.py's `labels` record (SILAC channels, masses and the channel-normalisation
    answer with its given/DEFAULT tag) from the first of these params files whose sidecar holds
    one, else None."""
    for path in paths:
        side = _load_json(f"{path}.rationale.json") if path else None
        if isinstance(side, dict) and isinstance(side.get("labels"), dict):
            return side["labels"]
    return None


def search_record(params=None, search_prov=None, manifest=None):
    """What the database search ran with, each value with where it came from.

    THE reader of a search's parameters for publication text: make_methods.py writes the
    Database-search paragraph from it and make_deposit.py writes the SDRF search columns from
    it, so the Methods and the repository metadata cannot disagree (DE-LIMP rule 3).

    Sources, most authoritative first: the parameters the search actually resolved to
    (search_provenance.json `resolved_params_file`), the params file it was given, then the
    workflow manifest (a pin, not a record of the run). A value in none of them is None --
    never a plausible default."""
    prov = _load_json(search_prov) or {}
    wf = _load_json(manifest) or {}
    engine = (prov.get("engine") or (wf.get("engine") or {}).get("name") or "").lower() or None
    rec = {"engine": engine, "engine_label": ENGINE_LABEL.get(engine, engine),
           "version": None, "version_source": None, "params_file": None,
           "params_source": None, "cleavage": None, "missed_cleavages": None,
           "pep_len": None, "pr_charge": None, "pr_mz": None, "mods": [],
           # True once the parameters' modifications were read: only then may an empty list be
           # said as "No fixed modifications were set" (L11 N4)
           "mods_read": False, "fdr_sentence": None,
           "max_var_mods": None, "met_excision": False, "ms1_tol": None, "ms2_tol": None,
           "tol_note": None, "precursor_fdr": None, "library": None, "mbr": None, "dda": None,
           "search_mode": prov.get("search_mode"),
           "labelled": None, "labels": None, "cont_quant_exclude": None, "cys_alkylation": None,
           "protein_inference": None,
           "probe_fallback": probe_fallback_record(prov),
           # a --dda search whose second pass crashed: its first-pass report was used
           "dda_first_pass": dda_first_pass_caution(prov),
           # ...and what the DE paragraph needs to say about it (dda_first_pass_de_sentence)
           "first_pass_record": {k: prov.get(k) for k in ("report_used", "dda_first_pass")},
           "warnings": []}
    # how a FragPipe search ran, from its record (rule 1; review F8): identifications only when
    # its quantification record says so -- the workflow file is only the fallback
    q = prov.get("quantification")
    if engine == "fragpipe" and isinstance(q, dict) and q.get("status") == "not_applicable":
        import fragpipe_params as fpp
        rec["route_sentence"] = ("The search was for identifications only, with no "
                                 f"quantification: "
                                 f"{fpp.methods_why(q.get('why')) or 'no quantification ran'}.")
    rec["route_from_record"] = isinstance(q, dict) and q.get("status") in ("not_applicable",
                                                                           "quantified")
    if prov.get("version"):
        rec["version"] = str(prov["version"])
        rec["version_source"] = "search_provenance.json (the version that ran)"
    elif prov.get("engine"):
        # a record whose version is null (run_search.engine_version_record: every conda-env
        # Sage): not recorded -- the manifest's pin is a request, never what ran (T44 F4;
        # session.search_engine, the same rule for the report, Slack and the README)
        rec["version_source"] = ("search_provenance.json names no version (the engine did not "
                                 "report one)")
    elif (wf.get("engine") or {}).get("version"):
        rec["version"] = f"{wf['engine']['version']} [pinned version — confirm it is what ran]"
        rec["version_source"] = ("workflow.manifest.json pin; the run's own record "
                                 "(search_provenance.json) was not found")

    # the parameters file: resolved (what ran) > given > manifest's
    cands = []
    if prov.get("resolved_params_file"):
        cands.append((prov["resolved_params_file"], "the resolved parameters the search ran "
                      "with (search_provenance.json)"))
    if search_prov and prov:
        # the record's own path is often HIVE's, unreadable where finalize runs: the resolved cfg
        # beside a copied search_provenance.json is the same file
        cands.append((os.path.join(os.path.dirname(os.path.abspath(search_prov)),
                                   "params.resolved.cfg"),
                      "the resolved parameters the search ran with (params.resolved.cfg beside "
                      "search_provenance.json)"))
    if prov.get("params_file"):
        cands.append((prov["params_file"], "the parameters file given to the search "
                      "(search_provenance.json)"))
    if params:
        cands.append((params, "the session's parameters file"))
    if (wf.get("search") or {}).get("params_file"):
        cands.append((wf["search"]["params_file"], "the workflow manifest's parameters file"))
    for path, src in cands:
        if path and os.path.isfile(path):
            rec["params_file"], rec["params_source"] = path, src
            break
    pf = rec["params_file"]
    if not pf:
        rec["warnings"].append("no search parameters file could be found from here")
        return rec
    where = os.path.basename(pf)
    # cysteine searched unmodified: the Methods say why (paragraph), from the params' sidecar
    rec["cys_alkylation"] = _cys_alkylation([pf] + [c for c, _ in cands])

    if engine == "diann" or pf.endswith(".cfg"):
        # the cfg tokeniser every other DIA-NN reader uses (diann_parallel.cfg_tokens)
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from diann_parallel import cfg_tokens, cfg_groups, CfgError, mass_acc_defaults, \
            mass_acc_status, protein_inference
        from estimate_params import SOP_DEFAULT_PUBLISHED, DIANN_DDA_FLAG
        try:
            groups = cfg_groups(cfg_tokens(pf))
        except CfgError as e:
            rec["warnings"].append(str(e))
            return rec
        flags = {}
        for flag, vals in groups:
            flags.setdefault(flag, []).append(vals)

        def one(f):                       # the last occurrence's first value, as DIA-NN reads it
            v = flags.get(f)
            return v[-1][0] if v and v[-1] else None

        def num(f):
            return _num(one(f))
        cut = one("--cut")
        if cut is not None:
            nm = DIANN_CUT.get(cut)
            rec["cleavage"] = {"rule": f"--cut {cut}", "name": nm[0] if nm else None,
                               "ac": nm[1] if nm else None, "source": where}
        if num("--missed-cleavages") is not None:
            rec["missed_cleavages"] = {"value": int(num("--missed-cleavages")), "source": where}
        for key, lo, hi in (("pep_len", "--min-pep-len", "--max-pep-len"),
                            ("pr_charge", "--min-pr-charge", "--max-pr-charge"),
                            ("pr_mz", "--min-pr-mz", "--max-pr-mz")):
            if num(lo) is not None or num(hi) is not None:
                rec[key] = {"value": (num(lo), num(hi)), "source": where}
        if "--unimod4" in flags:
            rec["mods"].append(_mod(4, "Carbamidomethyl", "fixed", "Anywhere", "C", 57.021464,
                                    f"--unimod4 in {where}"))
        for vals in flags.get("--fixed-mod", []):
            rec["mods"].append(_diann_mod(vals, "fixed", where))
        for vals in flags.get("--var-mod", []):
            rec["mods"].append(_diann_mod(vals, "variable", where))
        rec["mods_read"] = True
        if num("--var-mods") is not None:
            rec["max_var_mods"] = {"value": int(num("--var-mods")), "source": where}
        rec["met_excision"] = "--met-excision" in flags
        # --pg-level / --no-prot-inf as the cfg set them, or DIA-NN's measured default
        rec["protein_inference"] = protein_inference(groups)
        # DIA-NN searched the spectra as DDA: the paragraph must say so (DE-LIMP rule 1)
        if DIANN_DDA_FLAG in flags:
            rec["dda"] = {"value": True, "source": f"{DIANN_DDA_FLAG} in {where}"}
        ms2, ms1 = num("--mass-acc"), num("--mass-acc-ms1")
        if ms2 is None and ms1 is None:
            ma = ((prov.get("result") or {}).get("mass_acc") or {})
            # how the run set it: result.mass_acc, else the top-level record (run_search.py)
            how = ma if isinstance(ma, dict) and ma else (prov.get("mass_acc") or {})
            how = how if isinstance(how, dict) else {}
            ran_this_cfg = str(rec["params_source"] or "").startswith("the resolved parameters")
            if ma.get("fixed") and ma.get("ms2") is not None:
                ms2, ms1 = _num(ma.get("ms2")), _num(ma.get("ms1"))
                where_ma = "search_provenance.json result.mass_acc"
            elif prov and (how.get("mode") == "auto" or how.get("fixed") is False
                           or (not how and ran_this_cfg)):
                # the record says DIA-NN chose it -- or the cfg that RAN does not fix it
                where_ma = None
                rec["tol_note"] = ("mass accuracy was not fixed in the parameters; DIA-NN "
                                   "optimised it automatically per run")
            elif prov and (how.get("fixed") or how.get("mode") in ("measured", "pinned",
                                                                    "pinned_default",
                                                                    "fallback_default")):
                # 2.14.2 review: the chain and the single-shot probe record {mode: measured,
                # fixed: True, ms1/ms2: None}; the value is in a resolved cfg not readable here
                where_ma = None
                rec["tol_note"] = ("mass accuracy was measured before the search and pinned; "
                                   f"the value is not readable here {NR_TAG}")
            elif prov:
                where_ma = None
                rec["tol_note"] = ("the mass accuracy the search used is not recorded: the "
                                   "parameters read here do not fix it, and the run's record "
                                   f"does not say how it was set {NR_TAG}")
            else:
                # no record of the run: the parameters read here may be the ones the search was
                # GIVEN, not what it resolved to -- a search pinned at 10/20 ppm was once
                # described as optimised automatically (2.14.2, I04)
                where_ma = None
                rec["tol_note"] = ("the mass accuracy the search used is not recorded: the "
                                   "parameters read here do not fix it, and no record of the "
                                   f"run (search_provenance.json) was found {NR_TAG}")
        else:
            where_ma = where
        if ms2 is not None:
            rec["ms2_tol"] = {"value": ms2, "unit": "ppm", "source": where_ma}
        if ms1 is not None:
            rec["ms1_tol"] = {"value": ms1, "unit": "ppm", "source": where_ma}
        # DE-LIMP rule 2: a level estimate_params.py pinned at the SOP -- a DDA level it cannot
        # measure -- was nobody's choice for these data, and the Methods must not read as if it
        # were. Which levels those are comes from the run's provenance (search_provenance.json
        # result.mass_acc.default); with no run record, from the params file's own rationale
        # sidecar, through the same check (diann_parallel.mass_acc_defaults). A level is tagged
        # only while its value is the one the provenance names.
        ma_rec = (prov.get("result") or {}).get("mass_acc") or {}
        if isinstance(ma_rec.get("default"), dict):
            defaults = ma_rec["default"]
            default_src = "search_provenance.json result.mass_acc.default"
        else:
            try:
                defaults = mass_acc_defaults(pf, mass_acc_status(pf))
            except CfgError:
                defaults = {}
            default_src = f"{where}.rationale.json mass_accuracy_default"
        for key, flag in (("ms2_tol", "--mass-acc"), ("ms1_tol", "--mass-acc-ms1")):
            t, v = rec.get(key), _num(defaults.get(flag))
            if t and v is not None and float(t["value"]) == v:
                t["default"] = SOP_DEFAULT_PUBLISHED
                t["source"] = f"{t['source']}; DEFAULT per {default_src}"
        if num("--qvalue") is not None:
            rec["precursor_fdr"] = {"value": num("--qvalue"), "level": "precursor",
                                    "source": where}
        cmd_words = str(prov.get("resolved_command") or "").split()
        if one("--lib"):
            rec["library"] = {"value": f"spectral library {os.path.basename(one('--lib'))}",
                              "source": where}
        elif "--fasta-search" in flags or "--predictor" in flags:
            rec["library"] = {"value": "library-free: an in silico spectral library was "
                                       "predicted from the sequence database with DIA-NN's "
                                       "deep-learning predictor", "source": where}
        # DE-LIMP rule 1: match-between-runs is described from how the search RAN
        # (search_provenance.json search_mode), not from --reanalyse in the cfg. The 5-step
        # chain strips --reanalyse and does its own two-pass round trip as separate jobs; a
        # cfg that carries the flag said "with match-between-runs" for a DIA-NN MBR that never
        # ran as such.
        if rec["search_mode"] == "parallel_5step":
            seeded = bool((prov.get("result") or {}).get("seeded"))
            rec["mbr"] = {"value": "two_pass_chain", "label": "two-pass, as separate jobs "
                          "(empirical-library round trip)", "seeded": seeded,
                          "source": "search_provenance.json search_mode (parallel_5step)"}
        elif "--reanalyse" in flags or "--reanalyse" in cmd_words:
            rec["mbr"] = {"value": True, "label": "on", "source": where if "--reanalyse" in flags
                          else "search_provenance.json resolved_command"}
        # DIA-NN's own contaminant handling, as the run's log (else these parameters) set it --
        # {"value": None} when that was read and the flag is absent. The sidecar's
        # diann_cont_quant_exclude is only a recommendation, so it is never taken as proof that
        # the flag ran.
        rec["cont_quant_exclude"], rec["cont_quant_exclude_why"] = diann_cont_quant_exclude_why(
            (prov.get("result") or {}).get("report"), prov, pf)
        labelled = "--channels" in flags or any(m.get("label") for m in rec["mods"])
        rec["labelled"] = {"value": labelled, "source": where + (
            " (--channels / label mods present)" if labelled else
            " (no --channels or label modifications)")}
        # SILAC (plexDIA): the channels and the normalisation answer -- the params' sidecar (with
        # its given/DEFAULT tags), else as the cfg itself declares them (a search before 2.13)
        import silac
        rec["labels"] = _labels([pf] + [c for c, _ in cands]) or silac.labels_from_cfg(groups)
    elif engine == "sage" or pf.endswith(".json"):
        cfg = _load_json(pf) or {}
        # the search as Sage RAN it: its database block (enzyme, mods) and tolerances from the
        # results.json beside the record when there is one (review F16), else the config. quant
        # stays the config's: results.json always writes tmt_settings, which would read as TMT.
        ran = _sage_results_json(search_prov)
        run_where = "results.json (Sage's record of the search as run)"
        if isinstance((ran or {}).get("database"), dict):
            db, dwhere = ran["database"], run_where
        else:
            db, dwhere = cfg.get("database") or {}, where
        tcfg, twhere = ((ran, run_where) if ran and (ran.get("precursor_tol") or
                                                     ran.get("fragment_tol")) else (cfg, where))
        # the enzyme read as Sage reads it; an absent restrict is no proline rule, never 'P'
        enz = db.get("enzyme") if isinstance(db.get("enzyme"), dict) and db.get("enzyme") else None
        enz_where = dwhere
        if enz:
            rec["cleavage"] = sage_cleavage(enz, enz_where)
            if enz.get("missed_cleavages") is not None and \
                    rec["cleavage"]["specificity"] in ("semi", "full"):
                rec["missed_cleavages"] = {"value": int(enz["missed_cleavages"]),
                                           "source": enz_where}
            if enz.get("min_len") is not None or enz.get("max_len") is not None:
                rec["pep_len"] = {"value": (enz.get("min_len"), enz.get("max_len")),
                                  "source": enz_where}
        rec["mods"] += _sage_mods(db.get("static_mods"), "fixed", dwhere)
        rec["mods"] += _sage_mods(db.get("variable_mods"), "variable", dwhere)
        rec["mods_read"] = True
        if db.get("max_variable_mods") is not None:
            rec["max_var_mods"] = {"value": int(db["max_variable_mods"]), "source": dwhere}
        for key, name in (("ms1_tol", "precursor_tol"), ("ms2_tol", "fragment_tol")):
            tol = tcfg.get(name) or {}
            for unit in ("ppm", "da"):
                if isinstance(tol.get(unit), list) and len(tol[unit]) == 2:
                    lo, hi = (_num(x) for x in tol[unit])
                    rec[key] = {"value": hi if lo is not None and abs(lo) == hi else
                                f"{lo} to +{hi}", "unit": "ppm" if unit == "ppm" else "Da",
                                "symmetric": lo is not None and abs(lo) == hi,
                                "source": twhere}
        q = cfg.get("quant") or {}
        labelled = bool(q.get("tmt")) or bool(q.get("tmt_settings"))
        rec["labelled"] = {"value": labelled, "source": where + (
            " (quant.tmt set)" if labelled else
            " (quant.lfq false: identifications only)" if q.get("lfq") is False else
            " (quant.lfq, no TMT)")}
        # the search's own sentences, as run_search.py wrote them: quoted, never rebuilt here
        rec["sage_route"] = (prov.get("sage_input") or {}).get("methods_sentence")
        rec["sage_ids_fdr"] = (prov.get("sage_identifications") or {}).get("methods_sentence")
        # a quantifying search: the rule adapt_sage kept its rows with (sage_adapt), never a blank
        sa = prov.get("sage_adapt")
        if not isinstance(sa, dict) and search_prov:
            sa = _load_json(os.path.join(os.path.dirname(os.path.abspath(search_prov)),
                                         "sage_adapt.json"))
        rec["fdr_sentence"] = _sage_lfq_fdr(sa)
    elif engine == "fragpipe" or pf.endswith((".workflow", ".params")):
        _fragpipe_record(rec, prov, pf, search_prov)
    else:
        rec["warnings"].append(f"the {rec['engine_label'] or 'search'} parameter file "
                               f"{os.path.basename(pf)} is not parsed here; its settings are "
                               "in that file")
    return rec


def _g(x):
    """15.0 -> '15', 0.5 -> '0.5'; text passes through."""
    return ("%g" % x) if isinstance(x, (int, float)) else str(x)


def _fmt_range(v, unit=""):
    lo, hi = v
    f = lambda x: ("%g" % x) if isinstance(x, (int, float)) else str(x)
    if lo is not None and hi is not None:
        return f"{f(lo)}–{f(hi)}{unit}"
    return f"{'≥ ' + f(lo) if lo is not None else '≤ ' + f(hi)}{unit}"


def cleavage_phrase(cl):
    """The Methods wording of a search record's cleavage: non-specific, no digestion, a named
    enzyme (semi-specific when the search was), or the rule itself, tagged, when it has no name."""
    spec = cl.get("specificity")
    if spec == "none":
        return "no enzyme specificity (non-specific search)"
    if spec == "no_digestion":
        return "no in silico digestion (each database entry searched as one peptide)"
    if not cl.get("name"):
        return (f"the cleavage rule {cl['rule']} [enzyme name not mapped — confirm]"
                + (" (semi-specific)" if spec == "semi" else ""))
    return f"semi-specific {cl['name']} cleavage" if spec == "semi" else f"{cl['name']} specificity"


def mod_phrase(m):
    """'Oxidation (M)', 'Acetyl (Protein N-term)' -- the Methods wording of one modification."""
    where = m["position"] if m["position"] != "Anywhere" else ""
    tgt = m["targets"] or ""
    site = " ".join(x for x in (where, tgt) if x)
    name = m["name"] if m["unimod"] in UNIMOD else f"{m['name']} [name not verified — confirm]"
    return f"{name} ({site})" if site else name


def _first_pass_library(rec):
    """What the chain's first pass searched against: its seed library, or the predicted one."""
    return ("the seed spectral library" if (rec.get("mbr") or {}).get("seeded")
            else "the predicted library")


# The 5-step chain's match-between-runs, in the words of the proteomics review of 2026-09-29.
def chain_mbr_sentence(rec):
    first = _first_pass_library(rec)
    return ("Match-between-runs used DIA-NN's two-pass procedure, run as separate jobs: each run "
            f"was first searched against {first}; precursors identified across the experiment "
            "were assembled into an experiment-specific empirical spectral library; every run was "
            "then searched again against it. This differs from the retention-time-alignment-based "
            "match-between-runs of DDA software such as MaxQuant.")


def first_pass_used(de_prov):
    """Did the DE step read the 5-step chain's FIRST-pass report (diann_parallel.
    FIRST_PASS_REPORT), per de_provenance.json `input`? SKILL.md offers it when the final pass
    lost identifications (pass_comparison.py); the Methods must then describe that report."""
    inp = (de_prov or {}).get("input")
    if not inp:
        return False
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from diann_parallel import FIRST_PASS_REPORT
    return os.path.basename(str(inp)) == FIRST_PASS_REPORT


def first_pass_sentence(rec):
    """What was quantified, as a fact -- no "because", which would imply the first pass recovered
    what the second lost. It may not have: pass_comparison's flag is a pre-FDR count, and on SET28
    the flagged runs had no rows at 1% on every q-column in either report (what the first pass
    kept for them under run_de.R's filter were protein groups with a q between 1% and 5%)."""
    return ("Quantification used the experiment-wide first-pass report: each run was searched "
            f"once against {_first_pass_library(rec)}, with no match-between-runs; the second "
            "pass against an experiment-specific empirical library was not used.")


def probe_fallback_record(prov):
    """What a pre-search measurement that FAILED left the search with (probe_fallback.py: the
    stable scan_window.mode / mass_acc.mode, read by probe_fallback.fallback_modes()), or None
    when nothing fell back. DE-LIMP rule 2: the Methods must say the value was not measured,
    never present it as a measurement."""
    # the stable modes, read by the one reader the CAUTION uses too (DE-LIMP rule 3)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from probe_fallback import fallback_modes, WINDOW_FALLBACK, MASS_ACC_FALLBACK
    window, mass_acc = fallback_modes(prov)
    if WINDOW_FALLBACK != window and MASS_ACC_FALLBACK != mass_acc:
        return None
    fb = (prov or {}).get("probe_fallback")
    fb = fb if isinstance(fb, dict) else {}
    return {"window": window == WINDOW_FALLBACK, "mass_acc": mass_acc == MASS_ACC_FALLBACK,
            "reason": fb.get("reason"),
            "source": "search_provenance.json scan_window.mode / mass_acc.mode"}


def quantified_report_row(srec, de_prov):
    """The parameter table's row for the report that was quantified when it is a first pass --
    the chain's, chosen by the DE (de_provenance.json input), or a --dda search's whose second
    pass failed (search_provenance.json report_used) -- else the match-between-runs row; None
    when neither is recorded."""
    if first_pass_used(de_prov):
        return ("Report quantified", "the first pass (no match-between-runs)",
                "de_provenance.json input")
    if srec.get("dda_first_pass"):
        return ("Report quantified", "DIA-NN's first-pass report: the --dda second pass "
                "(match-between-runs) failed", "search_provenance.json report_used")
    if srec.get("mbr"):
        return ("Match-between-runs", srec["mbr"].get("label") or "on", srec["mbr"]["source"])
    return None


def dda_first_pass_caution(prov):
    """dda_first_pass.caution(): the one wording, for the search the Methods describe (this
    script's folder is on sys.path from its imports at the top) -- its Methods form, which states
    only that the report has no engine protein quantities; whether a DE estimated them is the DE
    paragraph's (dda_first_pass_de_sentence; C2 biologist B2 = statistician H2)."""
    from dda_first_pass import caution
    return caution(prov, methods=True)


def dda_first_pass_de_sentence(srec, de_prov):
    """dda_first_pass.de_sentence() for the search the Methods describe (the part of its record
    that says so, srec["first_pass_record"]) and the DE that ran: the DPC-Quant sentence, or
    None."""
    from dda_first_pass import de_sentence
    return de_sentence((srec or {}).get("first_pass_record"), de_prov)


def probe_fallback_sentence(rec):
    """The Methods sentence for a failed pre-search measurement, or "" when none failed."""
    fb = rec.get("probe_fallback")
    if not fb or not (fb["window"] or fb["mass_acc"]):
        return ""
    what = " and ".join(x for x in ("the scan window" if fb["window"] else None,
                                    "the mass accuracy" if fb["mass_acc"] else None) if x)
    then = []
    if fb["window"]:
        then.append("DIA-NN set the scan window automatically, for each run")
    if fb["mass_acc"]:
        then.append("the mass tolerances given here were not measured on these data")
    return (f"{what[0].upper() + what[1:]} {'were' if ' and ' in what else 'was'} to be measured "
            "on representative runs before the search, but that measurement failed; "
            + ", and ".join(then) + ".")


def ambiguous_search_record(candidates, why=None, how=None):
    """The search record when the session holds searches and none was chosen
    (session.search_choice; `why`, its reason): nothing is read -- not a params file, not the
    manifest's pin, either of which may describe any of them -- so the Methods and the SDRF carry
    tagged blanks rather than one search's values (2.14.2, I04). `how`: search_choice's
    ("unrecorded_input": the search the DE read has no record)."""
    rec = search_record(None, None, None)
    rec["ambiguous"], rec["warnings"] = len(candidates), []
    rec["unchosen_why"], rec["unchosen_how"] = why, how
    return rec


def search_fdr_sentence(rec, de_prov=None):
    """THE search's identification-FDR sentence, from its record (search_record): the search
    paragraph's, and the no-DE analysis brief's (C2 verify: the brief referred to an FDR rule it
    never stated). None when a DE that states its own makes it a repeat; "Identification FDR:
    NOT_RECORDED." when the record holds none."""
    fdr = rec.get("precursor_fdr")
    if fdr:
        return (f"Precursor identifications were filtered at {fdr['value'] * 100:g}% FDR "
                f"(q ≤ {fdr['value']:g}).")
    if rec.get("sage_ids_fdr"):
        return rec["sage_ids_fdr"]
    if not (de_prov or {}).get("q_columns"):
        # the search's own FDR rule when it is recorded (a quantifying Sage search's sage_adapt,
        # FragPipe's Philosopher filter); a DE that states its own is not repeated here
        return rec.get("fdr_sentence") or f"Identification FDR: {NOT_RECORDED}."
    if rec.get("engine") == "fragpipe" and rec.get("fdr_sentence") and \
            (de_prov or {}).get("q_columns_role") == "placeholder":
        # a FragPipe DDA DE's q-columns are 0.0 placeholders that filter nothing: Philosopher's
        # filter as run, and over which runs its rates hold, is the FDR statement (C2 prot M1).
        # Only then: a FragPipe DIA DE reads DIA-NN's report, whose q-values are real (its record
        # names no role), and there the Philosopher filter only built the library -- beside
        # DIA-NN's q-values it would read as the final filter (C2 verify, proteomics nit)
        return rec["fdr_sentence"]
    return None


def search_paragraph(rec, de_prov=None):
    """The Database-search paragraph. Every number comes from `rec`; anything missing prints
    as NOT_RECORDED rather than as a default."""
    if rec and rec.get("unchosen_how") == "unrecorded_input":
        return (f"Raw data were searched with {NOT_RECORDED} (the search the differential "
                f"expression read has no search record); the search parameters were "
                f"{NOT_RECORDED}.")
    if rec and rec.get("ambiguous"):
        return (f"Raw data were searched with {NOT_RECORDED} ({rec['ambiguous']} search records "
                f"were found for this analysis and none was chosen); the search parameters "
                f"were {NOT_RECORDED}.")
    if not rec or not rec.get("engine"):
        return (f"Raw data were searched with {NOT_RECORDED} (no search record — "
                "search_provenance.json / workflow manifest — was found).")
    eng = rec["engine_label"] or rec["engine"]
    ver = rec["version"] or NOT_RECORDED
    s = [f"Raw data were processed with {eng} {ver}"]
    if rec.get("library"):
        s[0] += f" ({rec['library']['value']})"
    mbr = (rec.get("mbr") or {}).get("value")
    first_pass = first_pass_used(de_prov)
    # its second pass failed, so no match-between-runs reached the report the caveat describes
    dda_fp = rec.get("dda_first_pass")
    s[0] += ", with match-between-runs." if mbr is True and not (first_pass or dda_fp) else "."
    for key in ("sage_route", "route_sentence"):       # how the search ran, in its own words
        if rec.get(key):
            s.append(rec[key])
    if rec.get("dda"):
        s.append("DIA-NN searched the spectra in its DDA mode (--dda), which it describes as "
                 "beta-stage support.")
    if dda_fp:
        s.append(dda_fp)
    if first_pass:
        s.append(first_pass_sentence(rec))
    elif mbr == "two_pass_chain":
        s.append(chain_mbr_sentence(rec))
    if not rec.get("params_file"):
        s.append(f"The search parameters could not be read from here: {NOT_RECORDED}.")
        return " ".join(s)
    digest = []
    cl = rec.get("cleavage")
    if cl:
        digest.append(cleavage_phrase(cl))
    if rec.get("missed_cleavages"):
        digest.append(f"up to {rec['missed_cleavages']['value']} missed cleavage"
                      + ("s" if rec["missed_cleavages"]["value"] != 1 else ""))
    if rec.get("pep_len"):
        digest.append(f"peptide length {_fmt_range(rec['pep_len']['value'])} residues")
    if rec.get("pr_charge"):
        digest.append(f"precursor charge {_fmt_range(rec['pr_charge']['value'])}")
    if rec.get("pr_mz"):
        digest.append(f"precursor m/z {_fmt_range(rec['pr_mz']['value'])}")
    if digest:
        s.append("In silico digestion used " + ", ".join(digest) + ".")
    if rec.get("met_excision"):
        s.append("N-terminal methionine excision was enabled.")
    fixed = [mod_phrase(m) for m in rec["mods"] if m["type"] == "fixed" and not m.get("label")]
    var = [mod_phrase(m) for m in rec["mods"] if m["type"] == "variable" and not m.get("label")]
    if not rec.get("mods_read"):
        # the parameters' modifications were not read: an empty list is not "none" (L11 N4)
        s.append(f"Fixed and variable modifications: {NOT_RECORDED}.")
    else:
        s.append(("Fixed modifications: " + "; ".join(fixed) + ". ") if fixed else
                 "No fixed modifications were set. ")
        # The fixed list above comes from the parameters, so it never names carbamidomethylation
        # that was not searched; when the samples were not alkylated (estimate_params
        # --cys-alkylation none), say why cysteine carries none.
        cys = rec.get("cys_alkylation") or {}
        cys_fixed = [m for m in rec["mods"]
                     if m["type"] == "fixed" and "C" in str(m.get("targets") or "")]
        if cys.get("value") == "none" and not cys_fixed:
            s[-1] += ("Cysteine was searched unmodified (no carbamidomethylation), as the samples "
                      "were not alkylated. ")     # what the flag says, no more (C2 proteomics nit)
        elif cys.get("reagent") and cys.get("value") != "cam" and any(
                m.get("unimod") == cys.get("unimod") for m in cys_fixed):
            s[-1] += f"Cysteines were alkylated with {cys['reagent']}. "
        s[-1] += (("Variable modifications: " + "; ".join(var)
                   + (f" (at most {rec['max_var_mods']['value']} per peptide)"
                      if rec.get("max_var_mods") else "") + ".") if var else
                  "No variable modifications were searched.")
    if rec.get("labels"):
        import silac
        s.append(silac.methods_search_sentence(rec["labels"]))
    tol = []
    for key, lvl in (("ms1_tol", "precursor (MS1)"), ("ms2_tol", "fragment (MS2)")):
        t = rec.get(key)
        if t:
            sym = ("±" if t.get("symmetric", True) and rec["engine"] in ("sage", "fragpipe")
                   else "")
            tol.append(f"{lvl} {sym}{_g(t['value'])} {t['unit']}"
                       + (f" ({t['default']})" if t.get("default") else ""))
    if tol:
        s.append("Mass tolerances were " + " and ".join(tol) + ".")
    db = rec.get("database")
    if db:
        s.append(f"MSFragger searched {db['name']}"
                 + (f" ({db['n_entries']:,} entries, {db['n_decoys']:,} of them decoys tagged "
                    f"{db['decoy_tag']}; sha256 {db['sha256'][:12]}…)" if db.get("sha256") else
                    " [its checksum is not recorded for this search — confirm]") + ".")
    if rec.get("fragger_extras"):
        s.append("MSFragger also searched " + "; ".join(rec["fragger_extras"])
                 + " [settings as fragger.params records them — confirm].")
    elif rec.get("tol_note"):
        s.append(rec["tol_note"][0].upper() + rec["tol_note"][1:] + ".")
    if probe_fallback_sentence(rec):
        s.append(probe_fallback_sentence(rec))
    if protein_inference_sentence(rec):
        s.append(protein_inference_sentence(rec))
    fdr = search_fdr_sentence(rec, de_prov)
    if fdr:
        s.append(fdr)
    for w in rec.get("warnings") or []:
        s.append(f"[{w} — confirm]")
    return " ".join(s)


def protein_inference_sentence(rec):
    """DIA-NN's protein inference, as the search's parameters set it (diann_parallel.
    protein_inference), in DIA-NN's own terms; "" for another engine. The default is named as
    the default, and as genes only for the DIA-NN version it was measured on."""
    pi = (rec or {}).get("protein_inference")
    if not isinstance(pi, dict) or rec.get("engine") not in (None, "diann"):
        return ""
    if pi.get("inference") is False:
        return "Protein inference was switched off (--no-prot-inf)."
    what = ("which decides which peptides count as proteotypic for the protein-level FDR; protein "
            "groups were assembled from the protein accessions each precursor maps to")
    if pi.get("passed"):
        basis = pi.get("basis") or f"{NOT_RECORDED} (an undocumented level)"
        return f"Protein inference used DIA-NN's implicit protein grouping by {basis} ({pi['passed']}), {what}."
    if str(rec.get("version") or "").startswith("2.7.0"):
        return (f"Protein inference used DIA-NN's implicit protein grouping by "
                f"{pi.get('basis')} (DIA-NN's default; no --pg-level was set), {what}.")
    return ("Protein inference used DIA-NN's default implicit protein grouping (no --pg-level was "
            f"set; measured as {pi.get('basis')} for DIA-NN 2.7.0), {what}.")


def diann_contaminant_sentence(srec):
    """What DIA-NN did with the contaminant entries, from the parameters the search ran with.
    It used to be read off the FASTA sidecar's diann_cont_quant_exclude -- a recommendation,
    not a record -- and worded as "excluded from quantification and normalisation" although
    run_de.R then re-quantified them (msalemi, 2026-09-24). Empty for a non-DIA-NN search."""
    srec = srec or {}
    if srec.get("engine") not in (None, "diann"):
        return ""
    cq = srec.get("cont_quant_exclude")
    if cq is None:
        return (f" Whether DIA-NN's --cont-quant-exclude was set: {NOT_RECORDED} "
                f"({srec.get('cont_quant_exclude_why') or 'no DIA-NN parameters file was read'}).")
    if not cq.get("value"):
        return (" DIA-NN's --cont-quant-exclude was not set, so contaminant peptides took part "
                "in DIA-NN's own normalisation.")
    tag = cq["value"]
    return (f" In DIA-NN (--cont-quant-exclude {tag}), peptides of {tag}-tagged entries were "
            f"excluded from normalisation and from the quantification of protein groups "
            f"containing no {tag} entry.")


# At most this many added sequences (and contaminant entries removed or kept beside them) are
# named one by one in the Methods: the bait / tag case (fetch_fasta.FEW_ADDED_SEQUENCES). Above
# it, the count and the file with its SHA-256 identify them -- an added second proteome of 4,506
# entries once put a 145 KB line of accessions into a client's Methods (staff report 2026-10-10).
ADDED_SEQUENCES_NAMED_MAX = FEW_ADDED_SEQUENCES


def added_sequences_sentence(fmeta):
    """The user's own target sequences (fetch_fasta.py --add-fasta: a bait, a tag), named with
    the file and SHA-256 they came from, the contaminant entries removed because most of each is
    an added protein, and those kept whose shared peptides are ambiguous
    (fetch_fasta.ADDED_SEQUENCE_RULE). Each list is named entry by entry only up to
    ADDED_SEQUENCES_NAMED_MAX, else counted. "" when none were added."""
    files = fmeta.get("added_sequences") or []
    names = [e.get("name") or e.get("accession") or "____"
             for f in files for e in (f.get("entries") or [])]
    # the entry count each file records (an entries list trimmed from a copy of the sidecar
    # still leaves n_entries)
    n = max(len(names), sum(f.get("n_entries") or 0 for f in files))
    if not n:
        return ""
    many = n != 1
    listed = len(names) == n and n <= ADDED_SEQUENCES_NAMED_MAX
    src = "; ".join(f"{os.path.basename(f.get('file') or '') or '____'}, SHA-256 "
                    f"{f.get('sha256') or '____'}" for f in files)
    sent = (f" {n:,} user-supplied sequence{'s' if many else ''}"
            + (f" ({', '.join(names)})" if listed else "")
            + f" {'were' if many else 'was'} added to the database as target "
              f"entr{'ies' if many else 'y'} ({src}).")
    gone = fmeta.get("contaminants_dropped_for_added_sequences") or []
    if gone:
        one = len(gone) == 1
        sent += (f" {len(gone):,} contaminant entr{'y' if one else 'ies'} sharing most of "
                 f"{'its' if one else 'their'} peptides with {'it' if not many else 'them'}"
                 + (" (" + ", ".join(f"{r.get('cont_acc')} ({r.get('cont_entry') or r.get('cont_gene') or '____'})"
                                     for r in gone) + ")"
                    if len(gone) <= ADDED_SEQUENCES_NAMED_MAX else "")
                 + f" {'was' if one else 'were'} removed from the library, so the added "
                   f"sequence{'s keep' if many else ' keeps'} those peptides.")
    kept = fmeta.get("contaminants_sharing_peptides_with_added_sequences") or []
    if kept:
        one = len(kept) == 1
        sent += (f" {len(kept):,} contaminant entr{'y' if one else 'ies'} sharing fewer of "
                 f"{'its' if one else 'their'} peptides with {'it' if not many else 'them'} "
                 f"{'was' if one else 'were'} kept, and the peptides shared are ambiguous between "
                 f"the contaminant and the added sequence")
        if len(kept) <= ADDED_SEQUENCES_NAMED_MAX:
            sent += (": " + "; ".join(
                f"{r.get('cont_acc')} ({r.get('cont_entry') or r.get('cont_gene') or '____'}) "
                f"and {r.get('target_acc')}: {', '.join(r.get('shared_peptides') or [])}"
                for r in kept) + ".")
        else:
            n_pep = len({p for r in kept for p in (r.get("shared_peptides") or [])})
            sent += f" ({n_pep:,} peptide{'' if n_pep == 1 else 's'})."
    return sent


def named_contaminants_sentence(fmeta):
    """The contaminant entries the analyst left out by name (fetch_fasta.py --drop-contaminant),
    named one by one up to ADDED_SEQUENCES_NAMED_MAX, else counted. "" when there were none -- the
    paragraph is then what it was before the flag existed."""
    gone = [r for r in (fmeta.get("contaminants_dropped_by_name") or []) if isinstance(r, dict)]
    if not gone:
        return ""
    one = len(gone) == 1
    return (f" {len(gone):,} contaminant entr{'y' if one else 'ies'}"
            + (" (" + ", ".join(f"{r.get('cont_acc') or '____'} "
                                f"({r.get('cont_entry') or r.get('cont_gene') or '____'})"
                                for r in gone) + ")"
               if len(gone) <= ADDED_SEQUENCES_NAMED_MAX else "")
            + f" {'was' if one else 'were'} left out of the database by name.")


def keratin_database_sentence(fmeta):
    """The FASTA sidecar's keratin-sample build (fetch_fasta.py --keratin-sample), or "" for any
    other database -- whose paragraph is then exactly what it was before keratin samples were
    handled."""
    if keratin_sample_recorded(fmeta) is not True:
        return ""
    n = fmeta.get("n_contaminants_dropped_keratin_sample")
    n = f"{n:,}" if isinstance(n, int) else "____"
    return (f" The samples are keratinous (hair, wool, feather, skin or nail: keratin is the "
            f"analyte), so the {n} keratin-family entries (keratins and keratin-associated "
            f"proteins) were also removed from the contaminant sequences, and keratin peptides "
            f"were identified and quantified as sample proteins; keratin introduced during "
            f"sample handling cannot be distinguished from the sample's own.")


def _count_of(n, c, missing="____"):
    """('<n> <unit>', verb) as contaminants.R count_of() words a count: the record's own unit
    (contaminant_unit), its singular for exactly one -- "1 protein groups" was printed
    (sage-review) -- and the verb that agrees. A record from before 2.9 counted precursors."""
    unit, unit1 = c.get("unit") or "precursors", c.get("unit_singular") or "precursor"
    if not isinstance(n, int):
        return f"{missing} {unit}".strip(), "were"
    return f"{n:,} {unit1 if n == 1 else unit}", "was" if n == 1 else "were"


def de_keratin_sentence(c):
    """A keratin sample's step of the DE, from the contaminant record's `keratin_sample`
    (contaminants.R keratin_sample_status); for a sample whose "not keratin" nobody confirmed,
    one tagged sentence when keratin precursors were removed on that default (rule 2); "" for a
    sample the user said is not keratin."""
    k = c.get("keratin_sample") if isinstance(c, dict) else None
    if not isinstance(k, dict):
        return ""
    # the unit the record counted in (contaminants.R contaminant_unit); older records: precursors
    unit = c.get("unit") or "precursors"
    if k.get("value") is not True:
        n = k.get("n_keratin_precursors_removed")
        if k.get("default_removed") is not True:
            return ""
        t, were = _count_of(n, c, missing="")
        return (f"{t} mapping to keratin-family contaminant entries {were} removed as "
                f"contaminants; whether samples are keratinous was not recorded (DEFAULT — not "
                f"user-confirmed).")
    head = ("The samples are keratinous (a keratin sample), so keratin is the analyte, not a "
            "contaminant: ")
    db = k.get("database")
    if db == "keratins_removed_at_build":
        n = k.get("n_removed_at_build")
        out = (head + f"the search database was built without its "
               f"{f'{n:,} ' if isinstance(n, int) else ''}keratin-family contaminant entries, so "
               f"keratin {unit} were quantified and tested as sample proteins.")
    elif db == "no_keratin_contaminants":
        out = (head + f"the search database held no keratin-family contaminant entry, so keratin "
               f"{unit} were quantified and tested as sample proteins.")
    elif db == "keratins_in_database":
        n_kept = k.get("n_precursors_kept")
        kept, were = _count_of(n_kept, c)
        n_ent = len(k.get("exempt_accessions") or [])
        # the pipeline's own descriptor says whether a kept item reaches a quantity
        # (requantified), and what the SEARCH ENGINE's own quantities did with those peptides
        # (kept_quant, per engine -- never assumed to be DIA-NN's; contaminants.R)
        if k.get("requantified") is True:
            out = (head + f"{kept} mapping only to keratin-family contaminant entries "
                   f"({n_ent:,} in the search database) {were} kept and quantified as sample "
                   f"proteins rather than removed.")
        else:
            out = (head + f"{kept} mapping only to keratin-family contaminant entries "
                   f"({n_ent:,} in the search database) {were} kept, so their protein groups stay "
                   f"in the analysis.")
        # none kept: what the engine did with kept ones is moot (sage-review N5)
        if n_kept != 0:
            eng = k.get("kept_quant")
            if isinstance(eng, str) and eng:
                out += f" {eng}"
            elif k.get("requantified") is not True:
                out += (f" What the search engine's own protein quantities did with those "
                        f"peptides: {NOT_RECORDED} (this DE record predates it).")
            under = k.get("under_quantified")
            if under is True:
                out += " [keratin under-quantified — resolve before publication]"
            elif under is None and k.get("requantified") is not True and eng:
                # the engine's own quantities, and whether its protein inference moved the
                # shared peptides elsewhere is not recorded (sage-review N4)
                out += (" [whether keratin is under-quantified is not known — check before "
                        "publication]")
    else:
        out = (head + f"the keratin-family entries of the search database could not be "
               f"identified, so {unit} mapping to them were removed with the other contaminants "
               f"[keratin under-counted — resolve before publication].")
    tags = [t for t in (k.get("unchecked_tags") or []) if isinstance(t, str)]
    if tags:
        t, were = _count_of(k.get("n_precursors_unchecked"), c, missing="")
        out += (f" {t} of "
                f"{'/'.join(tags)}-tagged contaminant entries -- contaminants the skill did not "
                f"add, never checked for keratins -- {were} removed with the other contaminants "
                f"[keratin may be under-counted — resolve before publication].")
    return out


def de_contaminant_sentence(prov):
    """The contaminant step of the DE, from run_de.R's `contaminants` record -- never assumed.
    A record older than the filter says so, tagged, instead of implying either answer. A keratin
    sample's own sentence (de_keratin_sentence) follows it."""
    c = prov.get("contaminants")
    ker = de_keratin_sentence(c)
    out = _de_contaminant_sentence(prov)
    return f"{out} {ker}" if ker else out


def _de_contaminant_sentence(prov):
    c = prov.get("contaminants")
    if not isinstance(c, dict):
        return (f"Contaminant handling in the differential-expression step: {NOT_RECORDED} "
                f"(this DE record predates it; run_de.R versions that did not record it did not "
                f"remove contaminants).")
    tag = c.get("tag") or CONT_TAG
    fmt = lambda k: f"{c.get(k):,}" if isinstance(c.get(k), int) else "____"  # noqa: E731
    policy = c.get("policy")
    # what one counted item is (contaminants.R contaminant_unit): DIA-NN precursors, Sage
    # peptides, a protein-level adapter's protein groups -- a record without it predates 2.9
    # and counted DIA-NN precursors (_count_of)
    unit1 = c.get("unit_singular") or "precursor"
    if policy == "removed":
        # THE rule text, as run_de.R recorded it (contaminants.R) -- not restated here
        rule = c.get("rule") or (f"any accession in {c.get('id_column') or '____'} starts with "
                                 f"{tag}")
        mixed = c.get("n_protein_groups_mixed")
        t, were = _count_of(c.get("n_precursors"), c)
        out = (f"Before protein quantification, {t} mapping to a "
               f"{tag}-tagged contaminant entry {were} removed ({rule}), taking out "
               f"{fmt('n_protein_groups')} contaminant protein groups"
               + (f" ({mixed:,} of them also naming a sample protein)"
                  if isinstance(mixed, int) and mixed else "")
               + ", so contaminants entered neither normalisation, the linear model nor the "
                 "multiple-testing correction.")
        if c.get("sample_loss"):
            out += f" {c['sample_loss']}"
        elif c.get("n_sample_groups_sharing") or c.get("n_sample_groups_all_shared"):
            out += (f" Sample protein groups sharing precursors with a contaminant entry lost "
                    f"those precursors ({fmt('n_sample_groups_sharing')} lost some, "
                    f"{fmt('n_sample_groups_all_shared')} lost all).")
        return out
    if policy == "kept":
        return (f"Contaminant protein groups ({fmt('n_protein_groups')} {tag}-tagged groups) were "
                f"kept in the differential-expression analysis (--keep-contaminants): they were "
                f"quantified, normalised and tested together with the sample proteins.")
    if policy == "none_present":
        return f"No identified {unit1} mapped to a {tag}-tagged contaminant entry."
    return (f"Contaminant filtering in the differential-expression step: {NOT_RECORDED} "
            f"({c.get('note') or 'not checked'}).")


def de_block_sentence(prov):
    """The random blocking factor (run_de.R --block), from its `block` record. None when the
    run fitted samples as independent -- or predates --block, which could only do that."""
    b = prov.get("block")
    if not isinstance(b, dict) or not b.get("applied"):
        return None
    col = b.get("column") or NOT_RECORDED
    rho = b.get("consensus_correlation")
    rho_s = f"{rho:.3f}" if isinstance(rho, (int, float)) else NOT_RECORDED
    n_est, n_all = b.get("n_proteins_estimated"), b.get("n_proteins")
    over = (f", estimated from {n_est:,} of {n_all:,} proteins"
            if isinstance(n_est, int) and isinstance(n_all, int) else "")
    levels = f" ({b['n_blocks']} levels)" if isinstance(b.get("n_blocks"), int) else ""
    if b.get("effect") == "fixed":
        absorbed = [str(x) for x in (b.get("absorbed_covariates") or [])]
        return (f"Samples sharing a {col} were paired: {col}{levels} was included in the linear "
                f"model as a fixed effect ({prov.get('design') or NOT_RECORDED}), since it is "
                f"crossed with the groups and every contrast compares samples within one {col} "
                f"-- the exact paired analysis."
                + (f" {', '.join(absorbed)} was left out of the design: every {col} sits in one "
                   f"{' / '.join(absorbed)}, so the {col} effect absorbs it." if absorbed else ""))
    out = (f"Samples sharing a {col} were modelled as correlated rather than independent: {col}"
           f"{levels} was fitted as a random blocking factor, with a consensus within-{col} "
           f"correlation of {rho_s} (limma duplicateCorrelation{over}) used in the linear-model "
           f"fit ({b.get('fit') or NOT_RECORDED}).")
    tr = b.get("technical_replicates") if isinstance(b.get("technical_replicates"), dict) else None
    if tr:
        # run_de.R's Sample column: injections of one biological sample (--reinjections all)
        out = (f"{tr.get('n_samples', NOT_RECORDED)} sample(s) were injected more than once "
               f"({tr.get('n_runs', NOT_RECORDED)} runs); these technical replicates were not "
               f"counted as independent samples. " + out)
    # Which fit reported each contrast (--block-scope): stated from the record, per contrast.
    model = b.get("contrast_model") if isinstance(b.get("contrast_model"), dict) else None
    if model is None:
        return out + (f" Which contrasts were reported from the blocked fit: {NOT_RECORDED} "
                      f"(this DE record predates --block-scope; that version reported all of them "
                      f"from it).")
    ind = [c for c, m in model.items() if m == "independent"]
    blk = [c for c, m in model.items() if m == "blocked"]
    if not ind:
        return out + " All contrasts were reported from this fit."
    return out + ((f" This fit reported {', '.join(blk)};" if blk else "")
                  + f" {', '.join(ind)} -- contrasts between different {col} levels using at most "
                  f"one sample per {col} -- were reported from the same data fitted with samples "
                  f"as independent, since there is no pairing to model and a single consensus "
                  f"correlation can understate their variance for proteins with strong "
                  f"{col}-to-{col} variation.")


def de_covariate_sentence(prov):
    """How the covariates entered the linear model (run_de.R's `covariates` record), when one was
    numeric (--numeric-covariates): its range, that one slope is shared by all groups, and -- for
    each comparison whose groups it nearly matches (run_de.R's per_contrast, confounded_contrasts)
    -- that comparison and why. None when every covariate
    was a factor (the design formula says it all), and for a record from before 2.12, which could
    only fit factors. The one description: the Methods, the analysis brief and the report's
    callout all read it."""
    rec = prov.get("covariates")
    terms = [t for t in (rec.get("terms") or []) if isinstance(t, dict)] if isinstance(rec, dict) else []
    numeric = [t for t in terms if t.get("type") == "numeric"]
    if not numeric:
        return None

    def rng(t):
        lo, hi = t.get("min"), t.get("max")
        return (f"values {lo:.4g} to {hi:.4g}" if isinstance(lo, (int, float))
                and isinstance(hi, (int, float)) else f"range {NOT_RECORDED}")
    names = [str(t.get("column") or NOT_RECORDED) for t in numeric]
    one = len(names) == 1
    out = (f"{' and '.join(names)} {'was' if one else 'were'} fitted as "
           f"{'a continuous covariate' if one else 'continuous covariates'} ("
           + (rng(numeric[0]) if one else
              "; ".join(f"{n}: {rng(t)}" for n, t in zip(names, numeric)))
           + ", as given in the conditions file), "
           + ("with one slope shared by all groups, so each group comparison is adjusted for it "
              "and compares the groups at the same covariate value." if one else
              "each with one slope shared by all groups, so each group comparison is adjusted "
              "for them and compares the groups at the same covariate values."))
    factors = [t for t in terms if t.get("type") == "factor"]
    if factors:
        out += " " + "; ".join(
            f"{t.get('column') or NOT_RECORDED} was a categorical factor"
            + (f" ({t['n_levels']} levels)" if isinstance(t.get("n_levels"), int) else "")
            for t in factors) + "."
    for t in numeric:
        if t.get("confounded") is not True:
            continue
        pc = t.get("per_contrast") if isinstance(t.get("per_contrast"), dict) else {}
        conf = [c for c in (t.get("confounded_contrasts") or []) if isinstance(pc.get(c), dict)]
        if conf:
            # per contrast (2.12): the comparisons whose groups the covariate separates, and why
            out += " " + " ".join(
                f"In {c}, {t.get('column')} nearly matched the groups compared ("
                + "; ".join(str(w) for w in (pc[c].get("why") or [])) + "), so that comparison "
                "rests largely on the fitted slope." for c in conf)
            continue
        # a record without per-contrast detail: the overall share
        r2 = t.get("group_r2")
        se = t.get("se_inflation") if isinstance(t.get("se_inflation"), dict) else {}
        se_max = max((v for v in se.values() if isinstance(v, (int, float))), default=None)
        out += (f" Group explained {r2 * 100:.0f}% of the variance of {t.get('column')}"
                if isinstance(r2, (int, float)) else
                f" {t.get('column')} was strongly confounded with Group")
        out += (", so the adjusted group comparisons rest largely on its fitted slope"
                + (f" (the design alone widens their standard errors up to {se_max:.1f}x)"
                   if se_max is not None else "") + ".")
    return out


def de_runs_left_out_sentence(prov):
    """The searched runs the DE left out, each with the reason the user gave, from run_de.R's
    `runs_left_out` record (the reasons come from collect_conditions.py --validate --exclude). A
    run with no recorded reason says NOT RECORDED -- never a reason supplied here -- or, when
    the record says its decisions file could not be read (`reasons_unreadable`),
    REASON_NOT_READ. A run of the design that had no quantities (`not_quantified`, 2.12) is said
    on its own, with run_de.R's reason. None when every searched run was analysed, or the record
    predates it. The one wording: the Methods and the report's callout
    (make_analysis_html.runs_left_out_note) both use it."""
    r = prov.get("runs_left_out")
    if not isinstance(r, dict):
        return None
    # design runs with no quantities (2.12): their own sentence, with run_de.R's reason
    nq = [x for x in (r.get("not_quantified") or []) if isinstance(x, dict) and x.get("run")]
    nq_text = (f"{len(nq)} run{' was' if len(nq) == 1 else 's were'} in the design but had no "
               f"quantities, so {'it was' if len(nq) == 1 else 'they were'} not analysed: "
               + "; ".join(f"{x['run']} ({str(x.get('reason') or 'reason ' + NOT_RECORDED).strip().rstrip('.')})"
                           for x in nq) + ".") if nq else None
    if not r.get("determined"):
        return nq_text
    runs = [x for x in (r.get("runs") or []) if isinstance(x, dict) and x.get("run")]
    if not runs:
        return nq_text

    unread = bool(r.get("reasons_unreadable"))

    def one(x):
        why = x.get("reason")
        why = why.strip().rstrip(".") if isinstance(why, str) and why.strip() else None
        return f"{x['run']} ({why or (REASON_NOT_READ if unread else 'reason ' + NOT_RECORDED)})"

    n_rep, n_an = r.get("n_report_runs"), r.get("n_analysed")
    n = len(runs)
    tail = (f"{'was' if n == 1 else 'were'} left out of the differential-expression analysis: "
            f"{'; '.join(one(x) for x in runs)}.")
    if isinstance(n_rep, int) and isinstance(n_an, int):
        out = f"Of the {n_rep} runs searched, {n_an} were analysed and {n} {tail}"
    else:
        out = f"{n} searched run{'' if n == 1 else 's'} {tail}"
    return out + (f" {nq_text}" if nq_text else "")


def de_normalisation_sentence(prov):
    """Which quantities the DE read and why, from run_de.R's record (`normalisation`, and the
    normalization_check block normalization_check.py wrote): the experiment type, the default
    applied for it, and whether the data check agreed -- or who chose otherwise, and why. A
    record without the check says NOT RECORDED; one from before 2.10 says nothing was recorded."""
    norm = prov.get("normalisation")
    nc = prov.get("normalization_check") if isinstance(prov.get("normalization_check"), dict) else None
    if not norm and not nc:
        return f"Between-run normalisation of the quantities: {NOT_RECORDED} (a DE record from before skill 2.10)."
    out = f"Between-run normalisation: {norm or NOT_RECORDED}."
    if not nc or nc.get("status") != "decided":
        return out + (f" Whether this suited the experiment was not checked {NOT_RECORDED} (no "
                      "normalisation check was run).")
    et_ = nc.get("experiment_type") or {}
    d = nc.get("decision") or {}
    # a person's decision, by role (normalization_check.statement); their name is staff-only
    chose = d.get("statement") or (
        f"{d.get('by') or NOT_RECORDED} chose "
        f"{(d.get('quantities') or NOT_RECORDED).replace('raw', 'non-normalised')} "
        f"quantities: {d.get('reason') or NOT_RECORDED}.")
    if not et_.get("type"):
        return out + (f" No experiment type was recorded, so no default applied; {chose}"
                      + (" A check of the data found: " + "; ".join(nc.get("trips")) + "."
                         if nc.get("trips") else ""))
    out += (f" Experiment type: {et_.get('label') or NOT_RECORDED} (per "
            f"{et_.get('source') or NOT_RECORDED}); the default for it is "
            f"{((nc.get('default') or {}).get('quantities') or NOT_RECORDED).replace('raw', 'non-normalised')} "
            "quantities" + (f" ({(nc.get('default') or {})['why']})"
                            if (nc.get("default") or {}).get("why") else "") + ".")
    if nc.get("tripped"):
        out += (" A check of the data (normalisation factors against the groups, the fold-change "
                "balance before normalisation, identifications per group, DIA-NN's normalisation "
                "stability and the volcano's shape) did not agree with that default ("
                + "; ".join(nc.get("trips") or []) + f"), and {chose}")
    else:
        out += (" A check of the data (normalisation factors against the groups, the fold-change "
                "balance before normalisation, identifications per group, DIA-NN's normalisation "
                "stability and the volcano's shape) agreed with it.")
    return out


def de_not_measured_sentence(prov):
    """run_de.R's `not_measured` rule (2.12.0), with each contrast's count from `de_tables`: a
    protein measured in no run of the groups a contrast compares was not tested in it. None for
    a record without the rule -- that run tested every protein, which the BH sentence covers."""
    rec = prov.get("not_measured") if isinstance(prov.get("not_measured"), dict) else None
    if not rec:
        return None
    tabs = prov.get("de_tables") if isinstance(prov.get("de_tables"), dict) else {}

    def count(t):
        n, of = t.get("n_not_measured"), t.get("n_proteins")
        if not isinstance(n, int):
            return NOT_RECORDED
        return f"{n:,} of {of:,}" if isinstance(of, int) else f"{n:,}"
    per = "; ".join(f"{c}: {count(t)}" for c, t in tabs.items() if isinstance(t, dict))
    return (f"In each contrast, a protein measured in no run of either compared group (DE "
            f"table Evidence '{rec.get('label') or NOT_RECORDED}') was not tested in that "
            f"contrast: it is listed without a fold change or p-value and was left out of that "
            f"contrast's Benjamini–Hochberg adjustment"
            + (f" ({per} proteins)." if per else f" (counts per contrast: {NOT_RECORDED})."))


DE_HEADING = "Differential expression"


def de_heading(de_prov):
    """The Methods' heading for the DE section: the pipeline's own (de_provenance.json
    `methods_heading`; run_silac.R's is "SILAC ratio analysis"), else "Differential expression".
    make_deposit.py checks a kept methods.md for this same heading -- one definition."""
    return str((de_prov or {}).get("methods_heading") or "").strip() or DE_HEADING


def identification_filter(prov):
    """THE reading of how a DE record's identifications were filtered -> (the pipeline's own
    FDR statement or None, the q-value filter it applied as text or None). A pipeline that
    states its own FDR (build_maxlfq.R descriptor, from the report's declaration) is quoted;
    its q-columns are named with their cutoffs only when they are real (q_columns_role) --
    placeholder columns (a Sage or FragPipe report adapted to DIA-NN's shape) filter nothing,
    and naming them would describe a filter that did not happen. The Methods and the analysis
    brief both read it (one definition)."""
    cols, cuts = prov.get("q_columns") or [], prov.get("q_cutoffs") or []
    ident = prov.get("identification_fdr") or None
    num = lambda x: f"{x:g}" if isinstance(x, (int, float)) else str(x)       # noqa: E731
    if ident and prov.get("q_columns_role") != "real":
        filt = None
    elif cols and len(cols) == len(cuts):
        filt = ", ".join(f"{c} ≤ {num(x)}" for c, x in zip(cols, cuts))
    elif prov.get("q_cutoff") is not None:
        filt = f"q ≤ {num(prov['q_cutoff'])}"
    else:
        filt = None
    return ident, filt


def de_paragraph(prov, srec=None):
    """The Differential-expression paragraph, from run_de.R's de_provenance.json. Significance
    is described exactly as run_de.R applied it: an adjusted-p cutoff, with |log2FC| only a
    reference line on the volcano (logfc_role) -- never as a second filter it did not apply.
    `srec`: the search record (search_record), for what the DE did about its report (a
    first-pass report with no engine protein quantities: dda_first_pass_de_sentence)."""
    def fmt_pkgs(p):
        pk = p.get("packages") or {}
        bits = [f"{n} {pk[n]}" for n in ("limpa", "limma") if pk.get(n)]
        if p.get("R_version"):
            bits.append(f"R {p['R_version']}")
        return f" ({', '.join(bits)})" if bits else ""
    s = [f"Differential expression was analysed with "
         f"{prov.get('display_label') or NOT_RECORDED}{fmt_pkgs(prov)}."]
    left_out = de_runs_left_out_sentence(prov)
    if left_out:
        s.append(left_out)
    if prov.get("rollup_method"):
        s.append(f"Protein quantities: {prov['rollup_method']}.")
    fp = dda_first_pass_de_sentence(srec, prov)
    if fp:
        s.append(fp)
    s.append(de_normalisation_sentence(prov))
    if prov.get("missing_policy"):
        s.append(prov["missing_policy"].rstrip(".") + ".")
    ident, filt = identification_filter(prov)
    if ident:
        s.append(f"Identifications entering quantification: {ident}.")
    if ident and prov.get("q_columns_role") != "real":
        pass
    elif filt:
        s.append(f"Identifications entering quantification were filtered at {filt}.")
    else:
        s.append(f"Identification q-value filter: {NOT_RECORDED}.")
    s.append(de_contaminant_sentence(prov))
    if prov.get("design"):
        # A record written before run_de.R listed contrasts holds a lone contrast as a bare
        # string; joining that would spell it out character by character.
        cons = prov.get("contrasts")
        cons = [cons] if isinstance(cons, str) else (cons or [])
        s.append(f"The linear model was {prov['design']}"
                 + (f", with contrasts {', '.join(cons)}" if cons else "") + ".")
    cov = de_covariate_sentence(prov)
    if cov:
        s.append(cov)
    blk = de_block_sentence(prov)
    if blk:
        s.append(blk)
    eng = prov.get("de_engine")
    adjp = prov.get("adjp")
    sig = (f"adj.P.Val < {adjp:g}" if isinstance(adjp, (int, float)) else
           f"adj.P.Val < {NOT_RECORDED}")
    s.append((f"Moderated t-statistics ({eng}) were computed and p-values adjusted by the "
              f"Benjamini–Hochberg method; proteins with {sig} were called significant."
              if eng else f"Proteins with {sig} (Benjamini–Hochberg) were called significant."))
    nm = de_not_measured_sentence(prov)
    if nm:
        s.append(nm)
    if prov.get("logfc_role") == "reference_line_only":
        lf = prov.get("logfc")
        s.append("No fold-change filter was applied"
                 + (f"; |log2FC| = {lf:g} is drawn on volcano plots for reference only."
                    if isinstance(lf, (int, float)) else "."))
    else:
        s.append(f"Fold-change filter: {NOT_RECORDED} (the DE record does not say whether "
                 f"one was applied).")
    if prov.get("citation"):
        s.append(f"Citation: {prov['citation']}.")
    return " ".join(s)


def acquisition_rows(rep, col, ser, n_files, record_source=None, acquisition=(None, None)):
    """(parameter, value, source) rows for the acquisition table -- every value with the file
    and field it was read from (bruker_method's `sources`), or the tag it carries.
    `acquisition`: acquisition_from_records()'s (mode, source), for a run whose own file does
    not record its mode (a .raw, an mzML, a run known only from the session record)."""
    src = rep.get("sources") or {}
    mode, mode_src = rep.get("mode"), record_source or src.get("mode", "Frames MsMsType")
    if (mode is None or record_source) and acquisition[0]:
        mode, mode_src = acquisition
    is_d = str(rep.get("file") or "").lower().rstrip("/").endswith(".d")
    s = rep.get("scheme") or {}

    def rng(lo, hi):
        return f"{lo}–{hi}" if lo is not None and hi is not None else None

    def unit(key):
        x = rep.get(key)
        return f"{_g(x['value'])} {_unit(x['unit'])}".strip() if x else None
    ramp, acc = rep.get("ramp_ms"), rep.get("accumulation_ms")
    chk = rep.get("ce_check") or {}
    ce = None
    if rep.get("ce_ramp"):
        ce = "; ".join(f"{_g(y)} eV at 1/K₀ {_g(x)}" for x, y in rep["ce_ramp"]["points"])
        ce += {"ok": f" (matches all {chk.get('n_windows')} window energies to within "
                     f"{chk.get('max_diff_ev')} eV)",
               "mismatch": f" [window energies differ by up to {chk.get('max_diff_ev')} eV — "
                           f"confirm]",
               "unchecked": " (not cross-checked: no window 1/K₀ bounds in the method)"
               }.get(chk.get("status"), "")
    per = s.get("per_ramp")
    rows = [
        ("Instrument", rep.get("instrument"), rep.get("instrument_source")),
        ("Instrument serial number", rep.get("instrument_serial"), src.get("instrument_serial")),
        ("Acquisition software", rep.get("software"), "analysis.tdf GlobalMetadata "
         "AcquisitionSoftware/-Version"),
        ("Control software (method file)", " ".join(x for x in (rep.get("control_software"),
         rep.get("control_software_version")) if x) or None, src.get("control_software")),
        ("HyStar version", rep.get("hystar_version"), src.get("hystar_version")),
        ("MS method", rep.get("ms_method"), src.get("ms_method")),
        ("Acquisition mode", mode, mode_src),
        ("Polarity", rep.get("polarity"), src.get("polarity")),
        ("m/z range", rng(rep.get("mz_low"), rep.get("mz_high")), src.get("mz_low")),
        ("1/K₀ range (V·s/cm²)", rng(rep.get("im_low"), rep.get("im_high")), src.get("im_low")),
        ("TIMS ramp / accumulation (ms)", f"{ramp} / {acc}" if ramp else None, src.get("ramp_ms")),
        ("Duty cycle", f"{100.0 * acc / ramp:.0f}%" if ramp and acc else None,
         "accumulation / ramp time"),
        ("Frames per cycle", rep.get("frames_per_cycle"), src.get("frames_per_cycle")),
        ("Cycle time (s)", rep.get("cycle_s"), src.get("cycle_s")),
        ("Isolation windows", rep.get("n_windows"), src.get("windows", "DiaFrameMsMsWindows")),
        ("TIMS ramps with windows (windows per ramp)",
         f"{s['n_ramps']} ({per[0] if per[0] == per[1] else f'{per[0]}–{per[1]}'})"
         if s.get("n_ramps") and per else None, src.get("windows")),
        ("Isolation width (Th)", rep.get("isolation_width"), "DiaFrameMsMsWindows IsolationWidth"),
        ("Window spacing / overlap (Th)", f"{_g(s['spacing'])} / {_g(s['overlap'])}"
         if s.get("spacing") is not None else None, "DiaFrameMsMsWindows IsolationMz"),
        ("Windows cover m/z", rng(s.get("mz_lo"), s.get("mz_hi")), src.get("windows")),
        ("Windows cover 1/K₀ (V·s/cm²)", rng(s.get("im_lo"), s.get("im_hi")),
         src.get("window_im")),
        ("Collision-energy ramp", ce, src.get("ce_ramp")),
        ("Collision energy per window (eV)", rng(rep.get("ce_low"), rep.get("ce_high")),
         "DiaFrameMsMsWindows CollisionEnergy"),
        ("Ion source", rep.get("source_type"), src.get("source_type")),
        ("Capillary voltage", unit("capillary_v"), src.get("capillary_v")),
        ("Dry gas", unit("dry_gas"), src.get("dry_gas")),
        ("Dry temperature", unit("dry_temp"), src.get("dry_temp")),
        ("LC system", " ".join(x for x in (rep.get("lc_system"), f"({rep['lc_vendor']})"
         if rep.get("lc_vendor") else None, f"S/N {rep['lc_serial']}" if rep.get("lc_serial")
         else None) if x) or None, src.get("lc_system")),
        ("LC method", rep.get("lc_method") or rep.get("lc_method_name"), src.get("lc_method")
         or src.get("lc_method_name")),
        ("LC run time (min)", rep.get("lc_run_min"), src.get("lc_run_min")),
        ("Gradient / flow", "the Evosep method's own (fixed, vendor-defined); not recorded in "
         "the .d" if "evosep" in (rep.get("lc_system") or "").lower() else None, "—"),
        ("Autosampler tray", rep.get("tray_type"), src.get("tray_type")),
        ("LC procedure log (on the LC PC)", rep.get("lc_log"), src.get("lc_log")),
        ("Analytical column", _with_tag(col["text"], col.get("tag")), col["source"]),
        ("Column temperature / emitter / mobile phases",
         f"not recorded in the {'.d' if is_d else 'raw file'} — confirm", "—"),
        ("Acquired", f"{ser['acquired_first']} to {ser['acquired_last']}"
         if ser.get("acquired_first") else None,
         "analysis.tdf GlobalMetadata AcquisitionDateTime (first – last run)"),
        ("Files in series", n_files, "this run"),
        ("Series consistency", ("all runs share the acquisition method" if not
         ser.get("differences") else "runs differ in: " + ", ".join(ser["differences"]))
         if ser.get("acquired_first") else None,
         "bruker_method.series (" + ", ".join(bruker_method.SERIES_KEYS[:4]) + ", …)"),
    ]
    return [(n, val, sc or "—") for n, val, sc in rows]


def qc_lines(rec, path):
    """The Instrument performance section: ONE fixed sentence about the Core's practice
    (qc_bracket.METHODS_SENTENCE -- the record's own copy when it carries one), printed whenever
    the project has a QC record, whatever it says. The verdict, and whether the check could run
    at all, are staff-only (Brett, 2026-10-01): a sentence that came and went with the check's
    outcome -- missing when STAN was unreachable or the record unreadable -- would tell the client
    what the verdict does not (2.10 safety review)."""
    readable = isinstance(rec, dict) and str(rec.get("schema", "")).startswith("qc_bracket/")
    if not readable and path:
        print(f"make_methods: {os.path.basename(path)} could not be read; the fixed Instrument "
              f"performance sentence is printed anyway", file=sys.stderr)
    sentence = rec.get("methods_sentence") if readable else None
    if not sentence:
        from qc_bracket import METHODS_SENTENCE     # the one wording
        sentence = METHODS_SENTENCE
    return ["## Instrument performance", "", sentence, ""]


def sample_prep_lines(rec, sr):
    """The Sample preparation section, from the CoreOmics submission (`sr` is
    submission_report). Who prepared the samples is sr.prepared_by()'s reading of the form.
    When the submitting lab sent peptides, every step before LC-MS/MS was theirs: the prose
    says so and carries no placeholder the Core could never fill. Nothing the form does not
    state is added -- its own words are quoted, in a note for the author, not in the prose."""
    who, why = sr.prepared_by(rec)
    src = (f"CoreOmics submission {sr.label(rec)}" if rec["source"] == sr.SOURCE_COREOMICS
           else f"submission {sr.label(rec)}, details given by the user")
    if who == "lab":
        # "peptides" only when the form's proteins/peptides answer says so.
        lines = ["Samples were prepared by the submitting laboratory and provided to the UC Davis "
                 "Proteomics Core" + (" as peptides ready for LC-MS/MS" if sr.sent_as_peptides(rec)
                                       else "") + f" ({src})."]
        # One line, no stray "*": the note must stay ONE italic line, which make_deposit drops
        # from the PRIDE protocol -- a multi-line quote would leak into it.
        said = [f"{k} “{' '.join(rec[f].split()).replace('*', '')}”"
                for k, f in (("buffer", "buffer"), ("beads", "beads")) if rec.get(f)]
        lines += ["", "*Describe the preparation from the submitting laboratory's own protocol."
                  + (f" As submitted: {'; '.join(said)}." if said else "") + "*"]
        return lines
    if who == "core":
        return [f"Samples were prepared by the UC Davis Proteomics Core ({src}); protocol: "
                f"{NOT_RECORDED}."]
    return [f"Who prepared the samples is not recorded: {why} ({src}). {NOT_RECORDED}"]


def database_build_warnings(fmeta):
    """fetch_fasta.py's build warnings that the Sequence-database paragraph does not already
    describe (its drop notes are records, written out in the sentences). They are for the agent
    and the user to resolve -- printed on stderr and in the JSON -- never in methods.md: a client's
    Methods once carried "Database build warnings (resolve before publication)" with a server path
    in it (staff report 2026-10-10)."""
    described = (fmeta.get("contaminants_dropped_note"),
                 fmeta.get("contaminants_dropped_for_added_sequences_note"),
                 fmeta.get("contaminants_sharing_peptides_with_added_sequences_note"))
    return [str(x) for x in (fmeta.get("warnings") or []) if x and x not in described]


# An absolute path: /quobyte/..., /nfs/..., /home/<user>/..., ~/..., C:\Users\... -- a server's or
# a laptop's file system, which means nothing to the Methods' reader and names internal storage.
# It starts the text, or follows a space or an opening quote/bracket -- never a word, a ")" or a
# "%": "S-Trap (ProtiFi)/Trypsin/Lys-C", "0.1%/80%/20%" and "(Biognosys)/DIA-NN/FragPipe" are
# prose (error-handling review, 2.12), as are m/z, 1/K0, Trypsin/P and URLs. And it starts with a
# root a file system has (FS_ROOTS), or ~/, or a drive letter -- or else has at least three
# segments made of file-name characters (/a/b/c.d).
FS_ROOTS = ("quobyte", "nfs", "home", "Users", "private", "tmp", "Volumes", "var", "mnt",
            "scratch", "group", "share", "data", "opt", "srv", "media", "lustre", "gpfs", "work")
_SEG = r"[^\s/\\'\"`<>|()\[\]{},;]+"
_NAME = r"[\w.@+-]+"
_ABS_PATH = re.compile(
    r"(?<![^\s(\[{\"'`])(?:"
    rf"/(?:{'|'.join(FS_ROOTS)})(?:/{_SEG})+/?"
    rf"|~(?:/{_SEG})+/?"
    rf"|/{_NAME}(?:/{_NAME}){{2,}}/?"
    rf"|[A-Za-z]:\\{_SEG}(?:\\{_SEG})*\\?)")


def client_text(text):
    """-> (text with every absolute path cut to its file name, [(path, what it became)]).
    Applied to the whole of methods.md as it is written, so no path from any record reaches the
    client; the caller echoes each cut on stderr, so a wrong one is seen."""
    cuts = []

    def base(m):
        after = re.split(r"[/\\]", m.group(0).rstrip("/\\"))[-1]
        cuts.append((m.group(0), after))
        return after
    return _ABS_PATH.sub(base, text), cuts


def _write_text(path, text):
    """UTF-8 whatever the platform's locale ("1/K₀" and "—" fail under Windows cp1252), and
    atomic: written to <path>.part and renamed, so a failure never leaves a 0-byte methods.md."""
    part = path + ".part"
    try:
        with open(part, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        os.replace(part, path)
    except BaseException:
        if os.path.exists(part):
            os.remove(part)
        raise


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw", nargs="+", required=True, help="raw file paths/globs (.d or .raw)")
    ap.add_argument("--out", default="methods.md")
    ap.add_argument("--lc-column", help="the analytical column, as it should read in the "
                                        "Methods (overrides every other source)")
    ap.add_argument("--column-log", help="a column-change log (CSV or JSON export of STAN's "
                                         "maintenance_events): the column in place at the first "
                                         "acquisition is used, with its install date as source")
    ap.add_argument("--column-log-instrument", help="the instrument to take from a column log "
                                                    "that covers several")
    ap.add_argument("--de-dir", help="optional: de_provenance.json for a Differential-expression paragraph")
    ap.add_argument("--fasta-meta", help="fetch_fasta.py's <fasta>.meta.json — writes the "
                                         "sequence-database sentence journals require")
    ap.add_argument("--params", help="the search parameters file (DIA-NN .cfg / Sage .json)")
    ap.add_argument("--search-prov", action="append",
                    help="run_search.py's search_provenance.json (engine, the version that ran, "
                         "the resolved parameters). Given more than once, none is chosen: the "
                         "Database-search paragraph is written with tagged blanks")
    ap.add_argument("--workflow-manifest", help="resolve_defaults.py's workflow.manifest.json")
    ap.add_argument("--instrument", help="instrument as the session recorded it; used ONLY when "
                                         "the raw files cannot be read from here")
    ap.add_argument("--acquisition", help="DIA/DDA as the session recorded it (step 2 detection)")
    ap.add_argument("--submission", help="the session dir (or a record) holding its CoreOmics "
                                         "submission: writes the Sample preparation section")
    ap.add_argument("--session", help="the analysis session: with no --search-prov, its search "
                                      "record is the one session.search_choice gives (wherever "
                                      "the search is; the one the DE read). Records and none "
                                      "chosen: tagged blanks, and why")
    ap.add_argument("--qc-bracket", help="qc_bracket.py's record (the session's "
                                         "logs/qc_bracket.json): its one fixed Instrument "
                                         "performance sentence -- the verdict stays staff-only")
    a = ap.parse_args()
    # several search records and none chosen (finalize passes every one it found): nothing of any
    # of them is described, and neither --params nor the manifest's pin stands in for them
    search_provs = a.search_prov or []
    unchosen = None                   # session.search_choice when it leaves the choice open
    if a.session and not search_provs:
        import session as sess
        choice = sess.search_choice(sess.paths_for(a.session))
        if choice["path"]:
            search_provs = [choice["path"]]
        elif sess.search_unchosen(choice):
            unchosen = choice
    # a record named by hand that is not there is refused, never read as "no record": the
    # paragraph would describe the manifest's pin and the session's params instead (T44 F8)
    gone = [x for x in search_provs if not os.path.isfile(x)]
    if gone:
        sys.exit(f"make_methods.py: --search-prov {gone[0]}: no such file. Name the session's "
                 "search record (output/<folder>/search_provenance.json -- finalize's "
                 "MANIFEST.txt names it), or give --session <session> to take the session's own")
    a.search_prov = search_provs[0] if len(search_provs) == 1 else None

    fmeta, build_warnings = None, []
    if a.fasta_meta:
        try:
            with open(a.fasta_meta, encoding="utf-8") as fh:
                fmeta = json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            sys.exit(f"--fasta-meta could not be read: {e}")

    files = []
    for p in a.raw:
        files.extend(sorted(glob.glob(p)) or [p])
    rec_instr = (a.instrument or "").strip() or None
    rec_src = REC_SRC
    metas = detect(files, rec_instr)
    acq, acq_src = acquisition_from_records(a)
    for m in metas:
        if m.get("vendor") == "Bruker" and "sources" in m:          # read from a .d
            m.setdefault("instrument_source", "GlobalMetadata InstrumentName")
        elif m.get("vendor") == "Thermo" and not m.get("instrument"):
            # a recorded instrument outranks the facility filename prefix, which is a guess
            if rec_instr:
                m["instrument"], m["instrument_source"] = rec_instr, rec_src
            elif m.get("prefix_instrument"):
                m["instrument"] = m["prefix_instrument"]
                m["instrument_source"] = "facility filename prefix (a guess — confirm)"
    # A .d that could not be read says so -- on stderr, in the note under Mass spectrometry,
    # and in the tag on its blank values -- rather than passing for an unrecorded value.
    problems = []
    for m in metas:
        for msg in ([f"could not be read ({m['error']})"] if m.get("error") else []) + \
                [f"read with warnings ({w})" for w in m.get("warnings") or []]:
            problems.append(f"{m.get('file')} {msg}")
            print(f"make_methods: {m.get('file')} {msg}", file=sys.stderr)
    unreadable = [m for m in metas if m.get("error")]
    from_record = False
    if not metas and rec_instr:
        # The raw files cannot be read from here (e.g. a session finalized away from the data).
        # Say so, and take ONLY instrument and acquisition mode from the session record: every
        # acquisition value the raw metadata would have given stays blank and tagged.
        metas = [{"vendor": instrument_vendor(rec_instr), "file": os.path.basename(f.rstrip("/")),
                  "instrument": rec_instr, "instrument_source": rec_src, "mode": acq}
                 for f in files]
        from_record = True
    if not metas:
        sys.exit("No raw files found (they may not be reachable from here); pass --instrument "
                 "(and --acquisition) from the session record to write the Methods anyway.")

    # representative metadata (facility usually acquires a series identically)
    bru = [m for m in metas if m.get("vendor") == "Bruker" and "sources" in m
           and m.get("instrument") and not m.get("error")]
    rep = bru[0] if bru else metas[0]
    if not rep.get("instrument") and rec_instr:
        rep = dict(rep, instrument=rec_instr, instrument_source=rec_src)
    instrument = rep.get("instrument") or next((m.get("instrument") for m in metas if m.get("instrument")), None)
    ack_label, ack_text = pick_ack(instrument, files)

    srec = (ambiguous_search_record(unchosen["candidates"], unchosen["why"], unchosen["how"])
            if unchosen else
            ambiguous_search_record(search_provs) if len(search_provs) > 1 else
            search_record(a.params, a.search_prov, a.workflow_manifest))
    de_prov = _load_json(os.path.join(a.de_dir, "de_provenance.json")) if a.de_dir else None

    # One paragraph describes the whole series only if every run was acquired the same way.
    ser = (bruker_method.series([{"values": m} for m in metas]) if not from_record
           else {"differences": {}, "acquired_first": None, "acquired_last": None})
    # Which paragraphs and defaults apply (2.14.2, S7/I09). The timsTOF MS paragraph is for a
    # run read from a .d, or one the record names a timsTOF; the Core's LC defaults (column,
    # emitter, column temperature) for a timsTOF .d only; the facility's mobile phases and
    # source for a run on a Core instrument (.d or .raw) -- never for an mzML or data from an
    # instrument the Core does not have, which get tagged blanks.
    vendor = rep.get("vendor")
    rep_file = str(rep.get("file") or "").lower().rstrip("/")
    is_tims = vendor == "Bruker" and ("sources" in rep or "tims" in (instrument or "").lower())
    core_tims = is_tims and rep_file.endswith(".d")
    facility = core_tims or (ack_label is not None and rep_file.endswith((".d", ".raw")))
    col = column_record(a.lc_column, metas if not from_record else (), a.column_log,
                        a.column_log_instrument, ser["acquired_first"], ser["acquired_last"],
                        core=core_tims)

    _write_text(os.path.splitext(a.out)[0] + "_params.json", json.dumps(
        {"files": [m.get("file") for m in metas], "representative": rep,
         "instrument": instrument, "vendor": vendor, "acknowledgment_for": ack_label,
         "all": metas, "from_session_record": from_record, "acquisition": acq,
         "acquisition_source": acq_src, "series": ser, "column": col, "search": srec,
         "read_problems": problems}, indent=2))

    # A blank acquisition value is one no record holds -- never a facility default (rule 2) --
    # or, when the representative raw file could not be read, simply unknown from here.
    blank_tag = UNREADABLE_TAG if (from_record or rep.get("error")) else NR_TAG

    def v(x, unit="", default=None):
        if x is None:
            return f"{default} {blank_tag}" if default is not None else f"____ {blank_tag}"
        return f"{x}{unit}"

    L, w = [], lambda s="": L.append(s)

    w("# Materials and Methods — LC-MS/MS")
    w("")
    if from_record:
        w(f"*Generated by the UC Davis Proteomics Core pipeline skill from the session record: "
          f"the {len(metas)} raw file(s) could not be read from where this was run, so the "
          f"acquisition values below are blank and tagged {UNREADABLE_TAG}, and the instrument "
          f"and acquisition mode come from the session record. Values marked {DEF} are "
          f"facility defaults to confirm; values marked {NR_TAG} were not in any record. Re-run "
          f"make_methods.py where the raw files are readable to fill them in.*")
    else:
        n_bad = len(unreadable)
        w(f"*Generated by the UC Davis Proteomics Core pipeline skill from the raw data "
          f"({len(metas)} file(s)"
          + (f"; {n_bad} could not be read — see the note under Mass spectrometry" if n_bad
             else "") + f"). Values marked {DEF} are facility defaults to confirm and values "
          f"marked {NR_TAG} were not in any record; every other value is listed with its "
          "source in the parameter tables below.*")
    w("")
    if a.submission:
        import submission_report
        rec, _session = submission_report.resolve(a.submission)
        if rec is None:
            sys.exit(f"--submission: no CoreOmics submission is attached to {a.submission}")
        w("## Sample preparation")
        w("")
        for line in sample_prep_lines(rec, submission_report):
            w(line)
        w("")
    w("## Liquid chromatography")
    w("")
    w(lc_paragraph(rep, col, core_tims, lc_known=bool(rep.get("lc_system")),
                   facility=facility))
    w("")
    w("## Mass spectrometry")
    w("")
    if is_tims:
        # a timsTOF not read from its .d (an mzML, the session record) states the recorded mode,
        # as the table does -- not "____ [not recorded]" beside a table row that says DIA
        w(ms_paragraph(rep if "sources" in rep else dict(rep, mode=rep.get("mode") or acq), v,
                       coupled=bool(rep.get("lc_system"))))
    else:
        mode = (f"{acq} mode (as detected from the data in step 2)" if acq in ("DIA", "DDA")
                else f"[DDA/DIA — confirm] mode {NR_TAG}")
        # "Thermo Orbitrap ... (Thermo Fisher Scientific)" names the vendor twice
        name = v(re.sub(r"(?i)^(Thermo|Bruker)\s+", "", rep.get("instrument") or "") or None,
                 default="[instrument]")
        art = "an" if name[:1].upper() in "AEIOU" else "a"
        maker = f"({VENDOR_NAME[vendor]})" if vendor in VENDOR_NAME else VENDOR_TAG
        what = ("resolution, AGC, isolation width, NCE, gradient" if vendor == "Thermo" else
                "scan range, isolation width, collision energy, gradient")
        w(f"Mass spectra were acquired on {art} {name} mass spectrometer {maker} operated in "
          f"{mode}. Full acquisition parameters ({what}) should be taken from the instrument "
          f"method file {NR_TAG}.")
    notes = [f"the runs differ in {k} ({', '.join(vals)})"
             for k, vals in ser["differences"].items()] + col["warnings"] + problems
    if notes:
        w("")
        w("> One paragraph cannot describe every run as it stands — resolve before publication: "
          + "; ".join(notes) + ".")
    w("")
    if a.qc_bracket:
        for line in qc_lines(_load_json(a.qc_bracket), a.qc_bracket):
            w(line)

    # Sequence database — journals require source, release, entry count, and how
    # contaminants were handled. Never invent these: if the sidecar wasn't passed,
    # emit a blank tagged line rather than a plausible-looking default.
    w("## Sequence database")
    w("")
    if fmeta:
        content_phrase = {
            "one_per_gene": "one canonical protein sequence per gene",
            "reviewed": "reviewed (Swiss-Prot) entries only",
            "reviewed_isoforms": "reviewed (Swiss-Prot) entries including splice isoforms",
            "full": "all entries including unreviewed (TrEMBL)",
            "full_isoforms": "all entries including unreviewed (TrEMBL) and splice isoforms",
        }.get(fmeta.get("content_used"))
        rel = f"release {rel}" if (rel := fmeta.get("uniprot_release")) else f"release ____ {NR_TAG}"
        n_p = fmeta.get("n_proteome")
        n_p = f"{n_p:,}" if isinstance(n_p, int) else "____"
        # Only call it a *reference* proteome when UniProt says it is one: a strain
        # assembly ("Non Reference proteome") or a user-supplied file is not, and
        # asserting otherwise puts a false claim in a published Methods section.
        kind = ("reference proteome"
                if (fmeta.get("proteome_type") or "").strip().lower() == "reference proteome"
                else "proteome")
        staged = fmeta.get("staged_file") if isinstance(fmeta.get("staged_file"), dict) else None
        if content_phrase is None and staged:
            # 'as_staged' (--hive), from a sidecar that describes the copy (gabrig,
            # 2026-09-23). Proteome and organism are known; the release it was cut from
            # is not -- name the copy's date rather than invent a release, and label a
            # composition inferred from entry counts as inferred.
            tax = f", taxid {fmeta['taxid']}" if fmeta.get("taxid") else ""
            sent = (f"Spectra were searched against a pre-staged copy of the UniProt "
                    f"{fmeta.get('organism') or '____'} {kind} "
                    f"({fmeta.get('proteome') or '____'}{tax}; copy dated "
                    f"{(staged.get('mtime_utc') or '')[:10] or '____'}, "
                    f"release ____ {NR_TAG}), comprising {n_p} sequences")
            g = (fmeta.get("content_check") or {}).get("uniprot_gene_count")
            if fmeta.get("content_inferred") == "one_per_gene" and isinstance(g, int):
                sent += (f"; the entry count is consistent with one canonical protein "
                         f"sequence per gene (inferred, not verified: UniProt lists "
                         f"{g:,} genes).")
            else:
                sent += f". Database composition: ____ {NR_TAG}."
        elif content_phrase is None:
            # 'unknown' (--path) / 'as_staged' (--hive): we did not build this database,
            # so we cannot describe its composition. Leave it tagged for the user. Its
            # organism is the user's answer (fetch --organism/--taxid; organism_source says so).
            org = fmeta.get("organism") or ""
            none = (fmeta.get("organism_source") or "").startswith("user: no single organism")
            tax = f"taxid {fmeta['taxid']}; " if fmeta.get("taxid") else ""
            sent = (f"Spectra were searched against a supplied "
                    f"{org + ' ' if org else ''}sequence database"
                    f"{' with no single source organism' if none else ''} "
                    f"({os.path.basename(fmeta.get('fasta', '') ) or '____'}; {tax}"
                    f"{n_p} sequences). Database composition and version: ____ {NR_TAG}.")
        else:
            sent = (f"Spectra were searched against the UniProt "
                    f"{fmeta.get('organism') or '____'} {kind} "
                    f"({fmeta.get('proteome') or '____'}, {rel}), comprising "
                    f"{content_phrase} ({n_p} sequences).")
        n_c = fmeta.get("n_contaminants_appended") or 0
        n_already = fmeta.get("n_contaminants_already_present") or 0
        if not n_c and n_already:
            sent += (f" The database already included {n_already} common-contaminant "
                     f"sequences.")
            sent += diann_contaminant_sentence(srec)
        elif n_c:
            sent += (f" A common-contaminant library ({n_c} sequences; "
                     f"{fmeta.get('contaminant_set')} set of Frankenfield et al., "
                     f"J Proteome Res 2022, 21:2104-2113) was appended.")
            sent += diann_contaminant_sentence(srec)
            # fetch_fasta.py removes contaminant entries whose sequence IS a target protein
            # (bovine ACTB = human ACTB, human keratins); a reader must know those proteins
            # were quantified, not excluded as contaminants.
            n_drop = fmeta.get("n_contaminants_dropped_as_target") or 0
            # min_unique_peptides > 0: built with the peptide rule too (near-identical
            # entries such as bovine EEF1A1 vs mouse). Absent: the identity rule alone.
            k = fmeta.get("min_unique_peptides") or 0
            if n_drop and k:
                sent += (f" {n_drop} contaminant entries that the search could not tell apart "
                         f"from {fmeta.get('organism') or '____'} proteins (identical or "
                         f"contained sequence, or fewer than {k} peptides of their own) were "
                         f"removed from the library first, so those proteins are quantified "
                         f"under their own accessions.")
            elif n_drop:
                sent += (f" {n_drop} contaminant entries identical to (or contained in) "
                         f"{fmeta.get('organism') or '____'} proteins were removed from the "
                         f"library first, so those proteins are quantified under their own "
                         f"accessions.")
        else:
            sent += " No contaminant database was appended."
        sent += named_contaminants_sentence(fmeta)
        sent += added_sequences_sentence(fmeta)
        sent += keratin_database_sentence(fmeta)
        w(sent)
        build_warnings = database_build_warnings(fmeta)
    else:
        w(f"Spectra were searched against {NOT_RECORDED} "
          f"(run `fetch_fasta.py` and pass `--fasta-meta <fasta>.meta.json` to fill "
          f"this in automatically).")
    w("")

    # the database search: engine, the version that ran, and the parameters it ran with
    if a.params or a.search_prov or a.workflow_manifest or srec.get("ambiguous"):
        w("## Database search")
        w("")
        w(search_paragraph(srec, de_prov))
        w("")

    # The DE paragraph from the skill's own run. It used to call the DE pipeline label the
    # search engine ("Raw files were searched and quantified with DPC-Quant + limma") and to
    # state "|log2FC| >= 1" as a significance threshold that run_de.R never applies.
    if de_prov:
        # A pipeline that writes its own Methods paragraph (run_silac.R) is quoted as written:
        # it describes itself (DE-LIMP rule 1); run_de.R's record is described by de_paragraph()
        w(f"## {de_heading(de_prov)}")
        w("")
        own = de_prov.get("methods_paragraph")
        left_out = de_runs_left_out_sentence(de_prov) if own else None
        w(" ".join(x for x in (own, left_out) if x) if own else de_paragraph(de_prov, srec))
        w("")
        cont = de_prov.get("contaminants") if isinstance(de_prov.get("contaminants"), dict) else {}
        if cont.get("database_risk") is True and cont.get("database_note"):
            w(f"> Contaminant filter caveat (resolve before publication): "
              f"{cont['database_note']}")
            w("")
        ker = cont.get("keratin_sample") if isinstance(cont.get("keratin_sample"), dict) else {}
        if ker.get("value") is True and ker.get("caution") is True and ker.get("note"):
            w(f"> Keratin-sample caveat (resolve before publication): {ker['note']}")
            w("")

    # The protein-set tests' paragraph, as run_sets.R wrote it (sets_provenance.json): its
    # defaults carry their "(DEFAULT — not user-confirmed)" tags from there.
    sets_prov = _load_json(os.path.join(a.de_dir, "sets_provenance.json")) if a.de_dir else None
    if sets_prov and sets_prov.get("methods_paragraph"):
        w("## Protein-set tests")
        w("")
        w(sets_prov["methods_paragraph"])
        w("")

    # The peak check's paragraph, as skyline_peaks.py wrote it (output/skyline/'s record, beside
    # the DE tables' folder): only when it ran -- a check that did not run has no method.
    # Its state is the one the delivery reads (skyline_peaks.record_state, rule 3): a STALE check
    # (its DE run again since) is held back from deliveries, so its Methods paragraph is too --
    # a delivered methods.md once described a check the delivery did not hold (2.14.2 review).
    if a.de_dir:
        import skyline_peaks
        pk_state, pk = skyline_peaks.record_state(os.path.dirname(os.path.abspath(a.de_dir)))
        if pk_state == "ok" and pk and pk.get("methods_paragraph"):
            w("## " + skyline_peaks.SECTION_TITLE)
            w("")
            w(pk["methods_paragraph"])
            w("")
        elif pk_state in skyline_peaks.STATE_WORDS:
            why = skyline_peaks.STATE_WORDS[pk_state]
            print(f"make_methods: {skyline_peaks.SECTION_TITLE} not described: {why}",
                  file=sys.stderr)
            w(f"> {skyline_peaks.SECTION_TITLE} not described: {why} [run it again and re-write "
              f"the Methods, or leave it out — confirm].")
            w("")

    # parameter table (value + source)
    w("## Acquisition parameters (extracted from the raw data)" if not from_record else
      "## Acquisition parameters (from the session record — raw files not readable here)")
    w("")
    w("| Parameter | Value | Source |")
    w("|---|---|---|")
    rows = acquisition_rows(rep, col, ser, len(metas), rec_src if from_record else None,
                            (acq, acq_src))
    for name, val, src in rows:
        if val is None: continue
        w(f"| {name} | {val} | {src} |")
    w("")

    if srec.get("engine"):
        w("## Search parameters (from the search record)")
        w("")
        w("| Parameter | Value | Source |")
        w("|---|---|---|")
        srows = [("Search engine", f"{srec['engine_label']} {srec['version'] or NOT_RECORDED}",
                  srec.get("version_source") or "not recorded"),
                 ("Parameters file", srec.get("params_file") and
                  os.path.basename(srec["params_file"]), srec.get("params_source"))]
        if srec.get("dda"):
            srows.append(("Spectra searched as", "DDA (--dda)", srec["dda"]["source"]))
        row = quantified_report_row(srec, de_prov)
        if row:
            srows.append(row)
        cl = srec.get("cleavage")
        if cl:
            srows.append(("Cleavage", f"{cl.get('name') or '[not mapped — confirm]'} "
                                      f"({cl['rule']})", cl["source"]))
        for key, label in (("missed_cleavages", "Missed cleavages"),
                           ("max_var_mods", "Max variable modifications")):
            if srec.get(key):
                srows.append((label, srec[key]["value"], srec[key]["source"]))
        for key, label in (("pep_len", "Peptide length"), ("pr_charge", "Precursor charge"),
                           ("pr_mz", "Precursor m/z")):
            if srec.get(key):
                srows.append((label, _fmt_range(srec[key]["value"]), srec[key]["source"]))
        for m in srec["mods"]:
            if m.get("label") and (srec.get("labels") or {}).get("channels"):
                continue                      # the SILAC label: the channel rows below say it
            srows.append((f"{m['type'].capitalize()} modification", mod_phrase(m), m["source"]))
        lab = srec.get("labels")
        if lab and lab.get("channels"):
            srows.append(("SILAC channels (plexDIA)", "; ".join(
                f"{c['name']} {c.get('preset') or 'K/R +' + '/'.join(str(c.get(x, {}).get('delta')) for x in 'KR')}"
                for c in lab["channels"]), lab.get("source") or "the search parameters"))
            cn = lab.get("channel_norm") or {}
            srows.append(("Channel normalisation", cn.get("flag") or NOT_RECORDED,
                          cn.get("source") or "not in the search parameters"))
        for key, label in (("ms1_tol", "Precursor (MS1) tolerance"),
                           ("ms2_tol", "Fragment (MS2) tolerance")):
            t = srec.get(key)
            srows.append((label, (f"{_g(t['value'])} {t['unit']}"
                                  + (f" ({t['default']})" if t.get("default") else "")) if t else
                          (srec.get("tol_note") or NOT_RECORDED),
                          t["source"] if t else "search record"))
        fb = srec.get("probe_fallback")
        if fb:
            if fb["window"]:
                srows.append(("Scan window", "set automatically by DIA-NN, per run (NOT "
                              "measured: the pre-search measurement failed)",
                              f"{fb['source']}: {fb.get('reason') or 'no reason recorded'}"))
            if fb["mass_acc"]:
                srows.append(("Mass accuracy measurement", "failed -- the tolerances above were "
                              "not measured on these data",
                              f"{fb['source']}: {fb.get('reason') or 'no reason recorded'}"))
        if srec.get("precursor_fdr"):
            f = srec["precursor_fdr"]
            srows.append(("Precursor FDR", f"q ≤ {f['value']:g}", f["source"]))
        for name, val, src in srows:
            if val is None or val == "":
                continue
            w(f"| {name} | {val} | {src} |")
        w("")

    w("## Acknowledgments")
    w("")
    # I09's rule: outside data gets no untagged Core statement -- the name matching a Core
    # instrument does not make the Core the place it was acquired (2.14.2 review)
    w(ack_text if facility or ack_label is None else f"{ack_text} {ACK_OUTSIDE_TAG}")
    w("")
    w(f"*Acknowledgment source: {ACK_SOURCE} (confirm the exact current wording before publishing).*")
    w("")

    text, cuts = client_text("\n".join(L) + "\n")
    for before, after in cuts:
        print(f"make_methods: a file path cut to its name in {os.path.basename(a.out)} (a server "
              f"path is not for the client): {before} -> {after}", file=sys.stderr)
    _write_text(a.out, text)
    for x in build_warnings:
        print(f"make_methods: database build warning (not in {os.path.basename(a.out)}): {x}",
              file=sys.stderr)
    print(json.dumps({"methods": os.path.abspath(a.out), "instrument": instrument,
                      "acknowledgment_for": ack_label, "n_files": len(metas),
                      "params_json": os.path.splitext(a.out)[0] + "_params.json",
                      # for the agent, never the client's Methods: tell the user, resolve them
                      "database_build_warnings": build_warnings,
                      "next": ("Tell the user the database build warnings above and resolve "
                               "them (they are not in methods.md). " if build_warnings else "")
                              + "Verify the draft against the params table, polish the prose, "
                                "then convert to .docx with to_docx.py."}, indent=2))


if __name__ == "__main__":
    if sys.argv[1:2] == ["cont-quant-exclude"]:
        # contaminants.R's question: make_methods.py cont-quant-exclude <report> -> JSON
        rec, why = diann_cont_quant_exclude_why(sys.argv[2] if len(sys.argv) > 2 else None)
        print(json.dumps(dict(rec or {"why": why}, recorded=rec is not None)))
        sys.exit(0)
    main()
