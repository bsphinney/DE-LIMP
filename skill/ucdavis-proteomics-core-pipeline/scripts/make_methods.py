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
from fetch_fasta import CONT_TAG, keratin_sample_recorded  # noqa: E402

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
# The facility's standard column and emitter, as STAN's column catalogue names them
# (STAN config/columns.yml `default_column_id` and the default emitter, 2026-09-03). Printed
# only with the DEF tag: they say what is usually fitted, not what was fitted for a given run.
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
# Sage: cleave_at + restrict (Sage DOCS.md v0.14.7: restrict = "do not cleave if this AA follows").
SAGE_CUT = {("KR", "P"): ("Trypsin", "MS:1001251"), ("KR", ""): ("Trypsin/P", "MS:1001313")}
# Unimod record -> (title, monoisotopic delta mass), from unimod.xml (2026-09-24). Only these are
# named from a mass; anything else is reported as the engine gave it, flagged for the user.
UNIMOD = {1: ("Acetyl", 42.010565), 4: ("Carbamidomethyl", 57.021464),
          7: ("Deamidated", 0.984016), 21: ("Phospho", 79.966331),
          28: ("Gln->pyro-Glu", -17.026549), 35: ("Oxidation", 15.994915)}


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
                  first=None, last=None):
    """THE analytical-column statement for the Methods, with where it came from (DE-LIMP
    rule 3: the prose and the parameter table both read it). Most authoritative first:
      1. --lc-column, given by the user for this run;
      2. HyStar ColumnInfo, when an operator entered the column in the LC method;
      3. the column log's latest column_change at or before the first acquisition;
      4. the facility's standard column, tagged DEF -- what is usually fitted, not a record.
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


def lc_paragraph(rep, col, is_bruker, lc_known):
    """The Liquid chromatography paragraph: LC system and method as the .d recorded them,
    the column from column_record(), and tagged placeholders for what no file records."""
    column = _with_tag(f"a {col['text']}", col.get("tag"))
    phases = ("Mobile phase A was 0.1% (v/v) formic acid in water and mobile phase B 0.1% "
              f"(v/v) formic acid in acetonitrile {DEF}.")
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


def detect(files):
    metas = []
    for f in files:
        low = f.lower().rstrip("/")
        if low.endswith(".d"):
            mm = bruker_meta(f)
        elif low.endswith(".raw"):
            mm = thermo_meta(f)
        else:
            mm = {"vendor": "?", "file": os.path.basename(f), "instrument": None}
        if mm: metas.append(mm)
    return metas


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
        d = os.path.dirname(os.path.abspath(report))
        for cand in (d, os.path.dirname(d)):
            prov = _load_json(os.path.join(cand, "search_provenance.json"))
            if prov is not None:
                break
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
           "max_var_mods": None, "met_excision": False, "ms1_tol": None, "ms2_tol": None,
           "tol_note": None, "precursor_fdr": None, "library": None, "mbr": None, "dda": None,
           "search_mode": prov.get("search_mode"),
           "labelled": None, "cont_quant_exclude": None,
           "probe_fallback": probe_fallback_record(prov), "warnings": []}
    if prov.get("version"):
        rec["version"] = str(prov["version"])
        rec["version_source"] = "search_provenance.json (the version that ran)"
    elif (wf.get("engine") or {}).get("version"):
        rec["version"] = f"{wf['engine']['version']} [pinned version — confirm it is what ran]"
        rec["version_source"] = ("workflow.manifest.json pin; the run's own record "
                                 "(search_provenance.json) was not found")

    # the parameters file: resolved (what ran) > given > manifest's
    cands = []
    if prov.get("resolved_params_file"):
        cands.append((prov["resolved_params_file"], "the resolved parameters the search ran "
                      "with (search_provenance.json)"))
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

    if engine == "diann" or pf.endswith(".cfg"):
        # the cfg tokeniser every other DIA-NN reader uses (diann_parallel.cfg_tokens)
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from diann_parallel import cfg_tokens, cfg_groups, CfgError, mass_acc_defaults, \
            mass_acc_status
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
        if num("--var-mods") is not None:
            rec["max_var_mods"] = {"value": int(num("--var-mods")), "source": where}
        rec["met_excision"] = "--met-excision" in flags
        # DIA-NN searched the spectra as DDA: the paragraph must say so (DE-LIMP rule 1)
        if DIANN_DDA_FLAG in flags:
            rec["dda"] = {"value": True, "source": f"{DIANN_DDA_FLAG} in {where}"}
        ms2, ms1 = num("--mass-acc"), num("--mass-acc-ms1")
        if ms2 is None and ms1 is None:
            ma = ((prov.get("result") or {}).get("mass_acc") or {})
            if ma.get("fixed") and ma.get("ms2") is not None:
                ms2, ms1 = _num(ma.get("ms2")), _num(ma.get("ms1"))
                where_ma = "search_provenance.json result.mass_acc"
            else:
                where_ma = None
                rec["tol_note"] = ("mass accuracy was not fixed in the parameters; DIA-NN "
                                   "optimised it automatically per run")
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
    elif engine == "sage" or pf.endswith(".json"):
        cfg = _load_json(pf) or {}
        db = cfg.get("database") or {}
        enz = db.get("enzyme") or {}
        if enz:
            # an absent key is Sage's documented default (DOCS.md: cleave_at 'KR', restrict 'P')
            key = (enz.get("cleave_at", "KR"), enz.get("restrict", "P") or "")
            nm = SAGE_CUT.get(key)
            rec["cleavage"] = {"rule": f"cleave_at={key[0]!r}, restrict={key[1]!r}",
                               "name": nm[0] if nm else None, "ac": nm[1] if nm else None,
                               "source": where}
            if enz.get("missed_cleavages") is not None:
                rec["missed_cleavages"] = {"value": int(enz["missed_cleavages"]), "source": where}
            if enz.get("min_len") is not None or enz.get("max_len") is not None:
                rec["pep_len"] = {"value": (enz.get("min_len"), enz.get("max_len")),
                                  "source": where}
        rec["mods"] += _sage_mods(db.get("static_mods"), "fixed", where)
        rec["mods"] += _sage_mods(db.get("variable_mods"), "variable", where)
        if db.get("max_variable_mods") is not None:
            rec["max_var_mods"] = {"value": int(db["max_variable_mods"]), "source": where}
        for key, name in (("ms1_tol", "precursor_tol"), ("ms2_tol", "fragment_tol")):
            tol = cfg.get(name) or {}
            for unit in ("ppm", "da"):
                if isinstance(tol.get(unit), list) and len(tol[unit]) == 2:
                    lo, hi = (_num(x) for x in tol[unit])
                    rec[key] = {"value": hi if lo is not None and abs(lo) == hi else
                                f"{lo} to +{hi}", "unit": "ppm" if unit == "ppm" else "Da",
                                "symmetric": lo is not None and abs(lo) == hi,
                                "source": where}
        q = cfg.get("quant") or {}
        labelled = bool(q.get("tmt")) or bool(q.get("tmt_settings"))
        rec["labelled"] = {"value": labelled, "source": where + (
            " (quant.tmt set)" if labelled else " (quant.lfq, no TMT)")}
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


def search_paragraph(rec, de_prov=None):
    """The Database-search paragraph. Every number comes from `rec`; anything missing prints
    as NOT_RECORDED rather than as a default."""
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
    s[0] += ", with match-between-runs." if mbr is True and not first_pass else "."
    if rec.get("dda"):
        s.append("DIA-NN searched the spectra in its DDA mode (--dda), which it describes as "
                 "beta-stage support.")
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
        digest.append(f"{cl['name']} specificity" if cl.get("name") else
                      f"the cleavage rule {cl['rule']} [enzyme name not mapped — confirm]")
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
    s.append(("Fixed modifications: " + "; ".join(fixed) + ". ") if fixed else
             "No fixed modifications were set. ")
    s[-1] += (("Variable modifications: " + "; ".join(var)
               + (f" (at most {rec['max_var_mods']['value']} per peptide)"
                  if rec.get("max_var_mods") else "") + ".") if var else
              "No variable modifications were searched.")
    tol = []
    for key, lvl in (("ms1_tol", "precursor (MS1)"), ("ms2_tol", "fragment (MS2)")):
        t = rec.get(key)
        if t:
            sym = "±" if t.get("symmetric", True) and rec["engine"] == "sage" else ""
            tol.append(f"{lvl} {sym}{_g(t['value'])} {t['unit']}"
                       + (f" ({t['default']})" if t.get("default") else ""))
    if tol:
        s.append("Mass tolerances were " + " and ".join(tol) + ".")
    elif rec.get("tol_note"):
        s.append(rec["tol_note"][0].upper() + rec["tol_note"][1:] + ".")
    if probe_fallback_sentence(rec):
        s.append(probe_fallback_sentence(rec))
    fdr = rec.get("precursor_fdr")
    if fdr:
        s.append(f"Precursor identifications were filtered at {fdr['value'] * 100:g}% FDR "
                 f"(q ≤ {fdr['value']:g}).")
    elif not (de_prov or {}).get("q_columns"):
        s.append(f"Identification FDR: {NOT_RECORDED}.")
    for w in rec.get("warnings") or []:
        s.append(f"[{w} — confirm]")
    return " ".join(s)


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


def added_sequences_sentence(fmeta):
    """The user's own target sequences (fetch_fasta.py --add-fasta: a bait, a tag), named with
    the file and SHA-256 they came from, the contaminant entries removed because most of each is
    an added protein, and those kept whose shared peptides are ambiguous
    (fetch_fasta.ADDED_SEQUENCE_RULE). "" when none were added."""
    files = fmeta.get("added_sequences") or []
    names = [e.get("name") or e.get("accession") or "____"
             for f in files for e in (f.get("entries") or [])]
    if not names:
        return ""
    many = len(names) != 1
    src = "; ".join(f"{os.path.basename(f.get('file') or '') or '____'}, SHA-256 "
                    f"{f.get('sha256') or '____'}" for f in files)
    sent = (f" {len(names)} user-supplied sequence{'s' if many else ''} ({', '.join(names)}) "
            f"{'were' if many else 'was'} added to the database as target "
            f"entr{'ies' if many else 'y'} ({src}).")
    gone = fmeta.get("contaminants_dropped_for_added_sequences") or []
    if gone:
        one = len(gone) == 1
        sent += (f" {len(gone)} contaminant entr{'y' if one else 'ies'} sharing most of "
                 f"{'its' if one else 'their'} peptides with {'it' if not many else 'them'} ("
                 + ", ".join(f"{r.get('cont_acc')} ({r.get('cont_entry') or r.get('cont_gene') or '____'})"
                             for r in gone)
                 + f") {'was' if one else 'were'} removed from the library, so the added "
                   f"sequence{'s keep' if many else ' keeps'} those peptides.")
    kept = fmeta.get("contaminants_sharing_peptides_with_added_sequences") or []
    if kept:
        one = len(kept) == 1
        sent += (f" {len(kept)} contaminant entr{'y' if one else 'ies'} sharing fewer of "
                 f"{'its' if one else 'their'} peptides with {'it' if not many else 'them'} "
                 f"{'was' if one else 'were'} kept, and the peptides shared are ambiguous between "
                 f"the contaminant and the added sequence: "
                 + "; ".join(f"{r.get('cont_acc')} ({r.get('cont_entry') or r.get('cont_gene') or '____'}) "
                             f"and {r.get('target_acc')}: {', '.join(r.get('shared_peptides') or [])}"
                             for r in kept) + ".")
    return sent


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


def de_runs_left_out_sentence(prov):
    """The searched runs the DE left out, each with the reason the user gave, from run_de.R's
    `runs_left_out` record (the reasons come from collect_conditions.py --validate --exclude). A
    run with no recorded reason says NOT RECORDED -- never a reason supplied here. None when every
    searched run was analysed, or the record predates it. The one wording: the Methods and the
    report's callout (make_analysis_html.runs_left_out_note) both use it."""
    r = prov.get("runs_left_out")
    if not isinstance(r, dict) or not r.get("determined"):
        return None
    runs = [x for x in (r.get("runs") or []) if isinstance(x, dict) and x.get("run")]
    if not runs:
        return None

    def one(x):
        why = x.get("reason")
        why = why.strip().rstrip(".") if isinstance(why, str) and why.strip() else None
        return f"{x['run']} ({why if why else 'reason ' + NOT_RECORDED})"

    n_rep, n_an = r.get("n_report_runs"), r.get("n_analysed")
    n = len(runs)
    tail = (f"{'was' if n == 1 else 'were'} left out of the differential-expression analysis: "
            f"{'; '.join(one(x) for x in runs)}.")
    if isinstance(n_rep, int) and isinstance(n_an, int):
        return f"Of the {n_rep} runs searched, {n_an} were analysed and {n} {tail}"
    return f"{n} searched run{'' if n == 1 else 's'} {tail}"


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


def de_paragraph(prov):
    """The Differential-expression paragraph, from run_de.R's de_provenance.json. Significance
    is described exactly as run_de.R applied it: an adjusted-p cutoff, with |log2FC| only a
    reference line on the volcano (logfc_role) -- never as a second filter it did not apply."""
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
    s.append(de_normalisation_sentence(prov))
    if prov.get("missing_policy"):
        s.append(prov["missing_policy"].rstrip(".") + ".")
    cols, cuts = prov.get("q_columns") or [], prov.get("q_cutoffs") or []
    ident = prov.get("identification_fdr")
    if ident:
        # a pipeline that states its own FDR (build_maxlfq.R descriptor, from the report's
        # declaration) is quoted; its q-columns are named with their cutoffs only when they are
        # real (q_columns_role) -- placeholder columns filter nothing, and naming them would
        # describe a filter that did not happen
        s.append(f"Identifications entering quantification: {ident}.")
    if ident and prov.get("q_columns_role") != "real":
        pass
    elif cols and len(cols) == len(cuts):
        s.append("Identifications entering quantification were filtered at "
                 + ", ".join(f"{c} ≤ {x:g}" for c, x in zip(cols, cuts)) + ".")
    elif prov.get("q_cutoff") is not None:
        s.append(f"Identifications entering quantification were filtered at q ≤ "
                 f"{prov['q_cutoff']:g}.")
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


def acquisition_rows(rep, col, ser, n_files, record_source=None):
    """(parameter, value, source) rows for the acquisition table -- every value with the file
    and field it was read from (bruker_method's `sources`), or the tag it carries."""
    src = rep.get("sources") or {}
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
        ("Acquisition mode", rep.get("mode"), record_source or src.get("mode", "Frames MsMsType")),
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
        ("Column temperature / emitter / mobile phases", "not recorded in the .d — confirm",
         "—"),
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
    ap.add_argument("--search-prov", help="run_search.py's search_provenance.json (engine, the "
                                          "version that ran, the resolved parameters)")
    ap.add_argument("--workflow-manifest", help="resolve_defaults.py's workflow.manifest.json")
    ap.add_argument("--instrument", help="instrument as the session recorded it; used ONLY when "
                                         "the raw files cannot be read from here")
    ap.add_argument("--acquisition", help="DIA/DDA as the session recorded it (step 2 detection)")
    ap.add_argument("--submission", help="the session dir (or a record) holding its CoreOmics "
                                         "submission: writes the Sample preparation section")
    ap.add_argument("--qc-bracket", help="qc_bracket.py's record (the session's "
                                         "logs/qc_bracket.json): its one fixed Instrument "
                                         "performance sentence -- the verdict stays staff-only")
    a = ap.parse_args()

    fmeta = None
    if a.fasta_meta:
        try:
            with open(a.fasta_meta, encoding="utf-8") as fh:
                fmeta = json.load(fh)
        except (OSError, json.JSONDecodeError) as e:
            sys.exit(f"--fasta-meta could not be read: {e}")

    files = []
    for p in a.raw:
        files.extend(sorted(glob.glob(p)) or [p])
    metas = detect(files)
    rec_instr = (a.instrument or "").strip() or None
    rec_src = "session record (workflow manifest) — not read from the raw file"
    for m in metas:
        if m.get("vendor") == "Bruker":
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
        low = rec_instr.lower()
        vendor = ("Bruker" if "tims" in low else
                  "Thermo" if any(k in low for k in ("orbitrap", "exploris", "exactive",
                                                     "lumos", "fusion", "eclipse", "astral",
                                                     "ascend")) else None)
        metas = [{"vendor": vendor, "file": os.path.basename(f.rstrip("/")),
                  "instrument": rec_instr, "instrument_source": rec_src,
                  "mode": (a.acquisition or "").upper() or None} for f in files]
        from_record = True
    if not metas:
        sys.exit("No raw files found (they may not be reachable from here); pass --instrument "
                 "(and --acquisition) from the session record to write the Methods anyway.")

    # representative metadata (facility usually acquires a series identically)
    bru = [m for m in metas if m.get("vendor") == "Bruker" and m.get("instrument")
           and not m.get("error")]
    rep = bru[0] if bru else metas[0]
    if not rep.get("instrument") and rec_instr:
        rep = dict(rep, instrument=rec_instr, instrument_source=rec_src)
    instrument = rep.get("instrument") or next((m.get("instrument") for m in metas if m.get("instrument")), None)
    ack_label, ack_text = pick_ack(instrument, files)

    srec = search_record(a.params, a.search_prov, a.workflow_manifest)
    de_prov = _load_json(os.path.join(a.de_dir, "de_provenance.json")) if a.de_dir else None

    # One paragraph describes the whole series only if every run was acquired the same way.
    ser = (bruker_method.series([{"values": m} for m in metas]) if not from_record
           else {"differences": {}, "acquired_first": None, "acquired_last": None})
    col = column_record(a.lc_column, metas if not from_record else (), a.column_log,
                        a.column_log_instrument, ser["acquired_first"], ser["acquired_last"])

    _write_text(os.path.splitext(a.out)[0] + "_params.json", json.dumps(
        {"files": [m.get("file") for m in metas], "representative": rep,
         "instrument": instrument, "acknowledgment_for": ack_label, "all": metas,
         "from_session_record": from_record, "acquisition": a.acquisition,
         "series": ser, "column": col, "search": srec, "read_problems": problems}, indent=2))

    # A blank acquisition value is one no record holds -- never a facility default (rule 2) --
    # or, when the representative raw file could not be read, simply unknown from here.
    blank_tag = UNREADABLE_TAG if (from_record or rep.get("error")) else NR_TAG

    def v(x, unit="", default=None):
        if x is None:
            return f"{default} {blank_tag}" if default is not None else f"____ {blank_tag}"
        return f"{x}{unit}"

    is_bruker = rep.get("vendor") == "Bruker"
    L, w = [], lambda s="": L.append(s)

    w("# Materials and Methods — LC-MS/MS")
    w("")
    if from_record:
        w(f"*Generated by the UC Davis Proteomics Core pipeline skill from the session record: "
          f"the {len(metas)} raw file(s) could not be read from where this was run, so the "
          f"acquisition values below are blank and tagged {UNREADABLE_TAG}, and the instrument "
          f"and acquisition mode come from the workflow manifest. Values marked {DEF} are "
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
    w(lc_paragraph(rep, col, is_bruker, lc_known=bool(rep.get("lc_system"))))
    w("")
    w("## Mass spectrometry")
    w("")
    if is_bruker:
        w(ms_paragraph(rep, v, coupled=bool(rep.get("lc_system"))))
    else:
        acq = (a.acquisition or "").upper()
        mode = (f"{acq} mode (as detected from the data in step 2)" if acq in ("DIA", "DDA")
                else f"[DDA/DIA — confirm] mode {NR_TAG}")
        # "Thermo Orbitrap ... (Thermo Fisher Scientific)" names the vendor twice
        thermo_name = v(re.sub(r"^Thermo\s+", "", rep.get("instrument") or "") or None,
                        default="[instrument]")
        art = "an" if thermo_name[:1] in "AEIOU" else "a"
        w(f"Mass spectra were acquired on {art} {thermo_name} mass "
          f"spectrometer (Thermo Fisher Scientific) operated in {mode}. "
          "Full acquisition parameters (resolution, AGC, isolation width, NCE, gradient) should be "
          f"taken from the instrument method file {NR_TAG}.")
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
        sent += added_sequences_sentence(fmeta)
        sent += keratin_database_sentence(fmeta)
        w(sent)
        # The drop notes are described in the sentences above -- they are records, not
        # something to resolve before publication.
        build_warnings = [x for x in (fmeta.get("warnings") or [])
                          if x not in (fmeta.get("contaminants_dropped_note"),
                                       fmeta.get("contaminants_dropped_for_added_sequences_note"),
                                       fmeta.get("contaminants_sharing_peptides_with_added_sequences_note"))]
        if build_warnings:
            w("")
            w(f"> Database build warnings (resolve before publication): "
              f"{'; '.join(build_warnings)}")
    else:
        w(f"Spectra were searched against {NOT_RECORDED} "
          f"(run `fetch_fasta.py` and pass `--fasta-meta <fasta>.meta.json` to fill "
          f"this in automatically).")
    w("")

    # the database search: engine, the version that ran, and the parameters it ran with
    if a.params or a.search_prov or a.workflow_manifest:
        w("## Database search")
        w("")
        w(search_paragraph(srec, de_prov))
        w("")

    # The DE paragraph from the skill's own run. It used to call the DE pipeline label the
    # search engine ("Raw files were searched and quantified with DPC-Quant + limma") and to
    # state "|log2FC| >= 1" as a significance threshold that run_de.R never applies.
    if de_prov:
        w("## Differential expression")
        w("")
        w(de_paragraph(de_prov))
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

    # parameter table (value + source)
    w("## Acquisition parameters (extracted from the raw data)" if not from_record else
      "## Acquisition parameters (from the session record — raw files not readable here)")
    w("")
    w("| Parameter | Value | Source |")
    w("|---|---|---|")
    rows = acquisition_rows(rep, col, ser, len(metas), rec_src if from_record else None)
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
        if first_pass_used(de_prov):
            srows.append(("Report quantified", "the first pass (no match-between-runs)",
                          "de_provenance.json input"))
        elif srec.get("mbr"):
            srows.append(("Match-between-runs", srec["mbr"].get("label") or "on",
                          srec["mbr"]["source"]))
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
            srows.append((f"{m['type'].capitalize()} modification", mod_phrase(m), m["source"]))
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
    w(ack_text)
    w("")
    w(f"*Acknowledgment source: {ACK_SOURCE} (confirm the exact current wording before publishing).*")
    w("")

    _write_text(a.out, "\n".join(L) + "\n")
    print(json.dumps({"methods": os.path.abspath(a.out), "instrument": instrument,
                      "acknowledgment_for": ack_label, "n_files": len(metas),
                      "params_json": os.path.splitext(a.out)[0] + "_params.json",
                      "next": "Verify the draft against the params table, polish the prose, then "
                              "convert to .docx with to_docx.py."}, indent=2))


if __name__ == "__main__":
    if sys.argv[1:2] == ["cont-quant-exclude"]:
        # contaminants.R's question: make_methods.py cont-quant-exclude <report> -> JSON
        rec, why = diann_cont_quant_exclude_why(sys.argv[2] if len(sys.argv) > 2 else None)
        print(json.dumps(dict(rec or {"why": why}, recorded=rec is not None)))
        sys.exit(0)
    main()
