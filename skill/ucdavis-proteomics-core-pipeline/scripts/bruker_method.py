#!/usr/bin/env python3
"""
bruker_method.py  --  How a Bruker timsTOF run was acquired, read from the run itself.

The Methods paragraph (make_methods.py) needs the LC system and method, the ion source, the
TIMS settings, the dia-PASEF window scheme and the collision-energy ramp. A .d records all of
them, in small files beside the spectra:

  analysis.tdf                  GlobalMetadata; Frames (timing, cycle, polarity);
                                DiaFrameMsMsWindows; the method's settings, as the Properties
                                of the first MS/MS frame (ion source, collision-energy ramp)
  <run>.m/diaSettings.diasqlite the dia-PASEF windows with their 1/K0 bounds
  <run>.m/microTOFQImpacTemAcquisition.method   the control software that wrote the method
  <run>.m/hystar.method         the LC method, and ColumnInfo when an operator entered a column
  HyStarMetadata.xml            the LC modules (device, vendor) and the LC method they ran
  SampleInfo.xml (UTF-16)       HyStar version, LC/MS method names, autosampler tray type

Verified against a timsTOF HT + Evosep One series (timsControl 6.0, HyStar 6.3.1.8; PROT_0756,
2026-08-13). Nothing here opens analysis.tdf_bin, and every sqlite file is opened read-only and
immutable (bruker_tdf.connect_tdf). A value the files do not record is None -- never a plausible
default (DE-LIMP rule 2) -- and each value carries the file and field it came from, in
`sources`. A newer or older .d that lacks a table or side file loses only those values.

What the files do NOT record, and so this never returns: the analytical column (unless an
operator typed it into HyStar's ColumnInfo), column temperature, emitter, mobile phases, and
the %B gradient of an Evosep method (it is Evosep's fixed, named method).

Stdlib only; imported by make_methods.py.
"""
import glob
import os
import re
import sqlite3
import statistics
import struct
import xml.etree.ElementTree as ET
from datetime import datetime

from bruker_tdf import connect_tdf

# Settings read from the first MS/MS frame (PropertyDefinitions.PermanentName).
PROPS = ("Source_Type", "Source_CapillarySetValue", "Source_DryGasSetValue",
         "Source_DryHeaterSetValue", "Mode_IonPolarity",
         "Energy_Ramping_Collision_Energy_Active", "Energy_Ramping_Advanced_Settings_Active",
         "Energy_Ramping_Mobility_StartEnd", "Energy_Ramping_Collision_Energy_StartEnd",
         "Energy_Ramping_Advanced_ListMobilityValues",
         "Energy_Ramping_Advanced_ListCollisionEnergyValues")
# A per-window collision energy may differ from the method's ramp at the window's 1/K0 midpoint
# by this much and the ramp is still what was applied (PROT_0756: at most 0.3 eV; a 20->59 eV
# ramp instead of the 20->65 eV one that ran misses the same windows by up to 3.7 eV).
CE_CHECK_EV = 1.0
MS_FRAME_TYPES = {0: "MS1", 8: "ddaPASEF", 9: "dia-PASEF"}
# The m/z range timsControl acquired the run over (analysis.tdf GlobalMetadata), e.g. 100-1700.
# One reader: the Methods text reports it (read_tdf), and detect_acquisition.py takes it as a
# ddaPASEF run's precursor m/z range -- every precursor it picked came from an MS1 frame
# acquired over it (a dia-PASEF run's range comes from its isolation windows instead).
MS1_RANGE_FIELDS = ("MzAcqRangeLower", "MzAcqRangeUpper")


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def doubles(blob):
    """A Properties array value: little-endian IEEE doubles packed into a blob."""
    if isinstance(blob, (bytes, bytearray)) and blob and len(blob) % 8 == 0:
        return list(struct.unpack("<%dd" % (len(blob) // 8), bytes(blob)))
    return None


def _enum(display_value_text):
    """PropertyDefinitions.DisplayValueText '0:No Source;1:ESI;...;11:Captive Spray' -> dict.
    The .d carries its own code table, so no code is interpreted from memory here."""
    out = {}
    for part in (display_value_text or "").split(";"):
        k, _, v = part.partition(":")
        if k.strip().lstrip("-").isdigit() and v:
            out[int(k)] = v.strip()
    return out


def ms1_acq_range(gm):
    """(lo, hi) from a GlobalMetadata {Key: Value} dict, or None when either end is missing,
    not a number, or not above the other."""
    lo, hi = (_num((gm or {}).get(f)) for f in MS1_RANGE_FIELDS)
    return (lo, hi) if lo is not None and hi is not None and hi > lo else None


def _tables(cur):
    return {r[0] for r in cur.execute("SELECT name FROM sqlite_master WHERE type IN "
                                      "('table','view')")}


def method_dir(d):
    """The .m method folder HyStar copied into the run (top level only -- the backup-*.m
    folders nested inside it are timsControl's own copies)."""
    ms = sorted(p for p in glob.glob(os.path.join(d, "*.m")) if os.path.isdir(p))
    return ms[0] if ms else None


def read_tdf(d, out, src):
    """analysis.tdf: identity, timing, windows and the method's settings."""
    tdf = os.path.join(d, "analysis.tdf")
    if not os.path.exists(tdf):
        return
    con = connect_tdf(tdf)
    try:
        cur = con.cursor()
        tabs = _tables(cur)
        if "GlobalMetadata" in tabs:
            gm = dict(cur.execute("SELECT Key, Value FROM GlobalMetadata"))
            for key, field in (("instrument", "InstrumentName"),
                               ("instrument_serial", "InstrumentSerialNumber"),
                               ("acquisition_software", "AcquisitionSoftware"),
                               ("acquisition_software_version", "AcquisitionSoftwareVersion"),
                               ("acquired_at", "AcquisitionDateTime"),
                               ("ms_method", "MethodName")):
                if gm.get(field):
                    out[key] = str(gm[field])
                    src[key] = f"analysis.tdf GlobalMetadata {field}"
            for key, field in (("mz_low", MS1_RANGE_FIELDS[0]), ("mz_high", MS1_RANGE_FIELDS[1]),
                               ("im_low", "OneOverK0AcqRangeLower"),
                               ("im_high", "OneOverK0AcqRangeUpper")):
                if _num(gm.get(field)) is not None:
                    out[key] = _num(gm[field])
                    src[key] = f"analysis.tdf GlobalMetadata {field}"
        if "Frames" in tabs:
            _frames(cur, out, src)
        if "DiaFrameMsMsWindows" in tabs:
            _windows(cur, out, src)
        if {"Properties", "PropertyDefinitions", "Frames"} <= tabs:
            _settings(cur, out, src)
    finally:
        con.close()


def _frames(cur, out, src):
    types = dict(cur.execute("SELECT MsMsType, COUNT(*) FROM Frames GROUP BY MsMsType"))
    out["mode"] = ("dia-PASEF" if types.get(9) else "ddaPASEF" if types.get(8) else
                   "MS" if types else None)
    if out["mode"]:
        src["mode"] = "analysis.tdf Frames.MsMsType"
    row = cur.execute("SELECT AccumulationTime, RampTime FROM Frames WHERE MsMsType IN (8,9) "
                      "LIMIT 1").fetchone()
    if row:
        out["accumulation_ms"], out["ramp_ms"] = _num(row[0]), _num(row[1])
        src["accumulation_ms"] = src["ramp_ms"] = "analysis.tdf Frames AccumulationTime/RampTime"
    try:
        pol = [r[0] for r in cur.execute("SELECT DISTINCT Polarity FROM Frames")]
        if len(pol) == 1 and pol[0] in ("+", "-"):
            out["polarity"] = "positive" if pol[0] == "+" else "negative"
            src["polarity"] = "analysis.tdf Frames.Polarity"
    except sqlite3.Error:
        pass
    try:
        t = [r[0] for r in cur.execute("SELECT Time FROM Frames WHERE MsMsType = 0 ORDER BY Id")]
    except sqlite3.Error:                     # a Frames table without Time (minimal fixtures)
        t = []
    n_ms1, n_all = types.get(0, 0), sum(types.values())
    if len(t) > 2 and n_ms1:
        out["cycle_s"] = round(statistics.median(b - a for a, b in zip(t, t[1:])), 2)
        src["cycle_s"] = "analysis.tdf Frames.Time (median MS1-to-MS1 interval)"
        out["frames_per_cycle"] = round(n_all / n_ms1)
        src["frames_per_cycle"] = "analysis.tdf Frames (all frames / MS1 frames)"
        out["run_min"] = round(max(t) / 60.0, 1)
        src["run_min"] = "analysis.tdf Frames.Time (last MS1 frame)"


def _windows(cur, out, src):
    cols = [c[1] for c in cur.execute("PRAGMA table_info(DiaFrameMsMsWindows)")]
    want = [c for c in ("WindowGroup", "IsolationMz", "IsolationWidth", "CollisionEnergy")
            if c in cols]
    rows = [dict(zip(want, r)) for r in
            cur.execute(f"SELECT {', '.join(want)} FROM DiaFrameMsMsWindows")]
    if not rows:
        return
    out["windows"] = rows
    src["windows"] = "analysis.tdf DiaFrameMsMsWindows"


def _settings(cur, out, src):
    """The method's settings as they applied to the first MS/MS frame (Properties view:
    GroupProperties overlaid by FrameProperties)."""
    first = cur.execute("SELECT MIN(Id) FROM Frames WHERE MsMsType IN (8,9)").fetchone()
    if not first or first[0] is None:
        first = cur.execute("SELECT MIN(Id) FROM Frames").fetchone()
    if not first or first[0] is None:
        return
    q = ("SELECT pd.PermanentName, p.Value, pd.DisplayValueText, pd.DisplayDimension "
         "FROM Properties p JOIN PropertyDefinitions pd ON pd.Id = p.Property "
         f"WHERE p.Frame = ? AND pd.PermanentName IN ({','.join('?' * len(PROPS))})")
    got = {name: (val, enum, dim) for name, val, enum, dim in cur.execute(q, (first[0],) + PROPS)}
    where = "analysis.tdf Properties (first MS/MS frame)"

    def num(name):
        return _num(got[name][0]) if name in got else None

    if num("Source_Type") is not None:
        code = int(num("Source_Type"))
        name = _enum(got["Source_Type"][1]).get(code)
        out["source_type"] = name or f"source type code {code}"
        src["source_type"] = f"{where} Source_Type = {code}, named by its DisplayValueText"
    for key, prop in (("capillary_v", "Source_CapillarySetValue"),
                      ("dry_gas", "Source_DryGasSetValue"),
                      ("dry_temp", "Source_DryHeaterSetValue")):
        if num(prop) is not None:
            out[key] = {"value": num(prop), "unit": got[prop][2] or ""}
            src[key] = f"{where} {prop}"
    if "polarity" not in out and num("Mode_IonPolarity") is not None:
        pol = _enum(got["Mode_IonPolarity"][1]).get(int(num("Mode_IonPolarity")))
        if pol:
            out["polarity"] = pol.lower()
            src["polarity"] = f"{where} Mode_IonPolarity"
    if num("Energy_Ramping_Collision_Energy_Active") == 0:
        out["ce_ramp"] = None
        src["ce_ramp"] = f"{where} Energy_Ramping_Collision_Energy_Active = 0 (no ramp)"
        return
    advanced = num("Energy_Ramping_Advanced_Settings_Active") == 1
    im_key, ce_key = (("Energy_Ramping_Advanced_ListMobilityValues",
                       "Energy_Ramping_Advanced_ListCollisionEnergyValues") if advanced else
                      ("Energy_Ramping_Mobility_StartEnd", "Energy_Ramping_Collision_Energy_StartEnd"))
    ims = doubles(got.get(im_key, (None,))[0])
    ces = doubles(got.get(ce_key, (None,))[0])
    if ims and ces and len(ims) == len(ces) >= 2:
        out["ce_ramp"] = {"points": sorted(zip(ims, ces)), "advanced": advanced}
        src["ce_ramp"] = f"{where} {im_key} / {ce_key}"


def read_dia_settings(mdir, out, src):
    """The dia-PASEF windows as the method defined them, with their 1/K0 bounds."""
    path = os.path.join(mdir, "diaSettings.diasqlite") if mdir else None
    if not path or not os.path.isfile(path):
        return
    con = connect_tdf(path)          # any sqlite file in a run: read-only and immutable
    try:
        rows = con.execute("SELECT CycleId, IsolationMz, OneOverK0Start, OneOverK0End "
                           "FROM DiaWindowsSpecification WHERE Type = 1").fetchall()
    except sqlite3.Error:
        rows = []
    finally:
        con.close()
    rows = [r for r in rows if None not in r]
    if rows:
        out["window_im"] = [{"group": int(g), "mz": float(mz), "im_lo": float(a),
                             "im_hi": float(b)} for g, mz, a, b in rows]
        src["window_im"] = f"{os.path.basename(mdir)}/diaSettings.diasqlite DiaWindowsSpecification"


def read_method_header(mdir, out, src):
    """The control software named in the MS method file's header (e.g. 'Bruker timsControl')."""
    path = os.path.join(mdir, "microTOFQImpacTemAcquisition.method") if mdir else None
    if not path or not os.path.isfile(path):
        return
    with open(path, "rb") as fh:
        head = fh.read(4096).decode("latin-1")
    tag = re.search(r"<fileinfo\b[^>]*>", head)
    if not tag:
        return
    for key, attr in (("control_software", "appname"), ("control_software_version", "appversion")):
        m = re.search(attr + r'="([^"]*)"', tag.group(0))
        if m and m.group(1).strip():
            out[key] = m.group(1).strip()
            src[key] = f"{os.path.basename(mdir)}/microTOFQImpacTemAcquisition.method fileinfo @{attr}"


def read_hystar_method(mdir, out, src):
    """ColumnInfo: the column an operator entered in HyStar, when one was entered."""
    path = os.path.join(mdir, "hystar.method") if mdir else None
    if not path or not os.path.isfile(path):
        return
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError):
        return
    node = root.find(".//ColumnInfo")
    if node is not None:
        text = " ".join("".join(node.itertext()).split())
        out["column_info"] = text or None
        src["column_info"] = (f"{os.path.basename(mdir)}/hystar.method ColumnInfo"
                              + ("" if text else " (empty -- no column entered)"))


def _local(tag):
    return tag.rsplit("}", 1)[-1]


def read_hystar_metadata(d, out, src):
    """HyStarMetadata.xml: the LC device, its vendor, and the method it ran."""
    path = os.path.join(d, "HyStarMetadata.xml")
    if not os.path.isfile(path):
        return
    try:
        root = ET.parse(path).getroot()
    except (ET.ParseError, OSError):
        return
    params = {}
    for p in root.iter():
        if _local(p.tag) != "Parameter" or not p.get("ID"):
            continue
        val = next((c.text for c in p if _local(c.tag) in ("Value", "Int", "Double", "DateTime")
                    and c.text), None)
        params.setdefault(p.get("ID"), val)
    where = "HyStarMetadata.xml Parameter"
    for key, pid in (("lc_system", "ConfigurationObject_DeviceName"),
                     ("lc_method", "MethodDataObject_MethodName"),
                     ("lc_run_min", "MethodDataObject_RunTimeMinutes"),
                     ("lc_serial", "ConfigurationObject_SerialNumber"),
                     ("lc_log", "LOG_FILE_PATH"),
                     ("ms_control", "MS Control")):
        if params.get(pid):
            out[key] = _num(params[pid]) if key == "lc_run_min" else params[pid].strip()
            src[key] = f"{where} {pid}"
    # the module table names each LC module's vendor
    for table in root.iter():
        if _local(table.tag) != "Table":
            continue
        names = {c.get("ColIndex"): c.get("Name") for c in table if _local(c.tag) == "Col"}
        if "Vendor" not in names.values():
            continue
        for row in (r for r in table if _local(r.tag) == "Row"):
            cells = {names.get(c.get("ColIndex")): "".join(c.itertext()).strip()
                     for c in row if _local(c.tag) == "Cell"}
            if out.get("lc_system") and cells.get("Name") == out["lc_system"] and cells.get("Vendor"):
                out["lc_vendor"] = cells["Vendor"]
                src["lc_vendor"] = "HyStarMetadata.xml module table (Vendor)"


def read_sample_info(d, out, src):
    """SampleInfo.xml (UTF-16): HyStar version, method names, autosampler tray type."""
    path = os.path.join(d, "SampleInfo.xml")
    if not os.path.isfile(path):
        return
    try:
        with open(path, "rb") as fh:
            raw = fh.read()
        text = raw.decode("utf-16") if raw[:2] in (b"\xff\xfe", b"\xfe\xff") else raw.decode("utf-8")
        root = ET.fromstring(re.sub(r"^\s*<\?xml[^>]*\?>", "", text))
    except (ET.ParseError, OSError, UnicodeDecodeError):
        return
    head = root.find("SampleTableHeader")
    if head is not None and head.get("HyStarVersion"):
        out["hystar_version"] = head.get("HyStarVersion")
        src["hystar_version"] = "SampleInfo.xml SampleTableHeader@HyStarVersion"
    props = {p.get("Name"): p.get("Value") for p in root.iter("Property")}
    for key, name in (("tray_type", "TrayType"), ("lc_method_name", "HyStar_LC_Method_Name")):
        if props.get(name):
            out[key] = props[name]
            src[key] = f"SampleInfo.xml Property {name}"


def read_run(d):
    """Everything this module reads from one .d: {'values': {...}, 'sources': {...}}."""
    out, src = {}, {}
    try:
        read_tdf(d, out, src)
    except sqlite3.Error as e:
        out["error"] = str(e)
    mdir = method_dir(d)
    for reader in (read_dia_settings, read_method_header, read_hystar_method):
        try:
            reader(mdir, out, src)
        except (OSError, sqlite3.Error) as e:
            out.setdefault("warnings", []).append(f"{reader.__name__}: {e}")
    for reader in (read_hystar_metadata, read_sample_info):
        try:
            reader(d, out, src)
        except OSError as e:
            out.setdefault("warnings", []).append(f"{reader.__name__}: {e}")
    return {"values": out, "sources": src}


def window_scheme(v):
    """The dia-PASEF scheme a reader of the Methods needs: window count, TIMS ramps, width,
    spacing/overlap, and the m/z and 1/K0 area covered. Only what the windows say."""
    wins = v.get("windows") or []
    if not wins:
        return None
    widths = sorted({round(float(w["IsolationWidth"]), 3) for w in wins
                     if w.get("IsolationWidth") is not None})
    s = {"n_windows": len(wins)}
    groups = {}
    for w in wins:
        groups.setdefault(w.get("WindowGroup"), []).append(w)
    if None not in groups:
        s["n_ramps"] = len(groups)
        per = [len(g) for g in groups.values()]
        s["per_ramp"] = (min(per), max(per))
    if widths:
        s["width"] = widths[0] if len(widths) == 1 else (widths[0], widths[-1])
    if all(w.get("IsolationMz") is not None and w.get("IsolationWidth") is not None for w in wins):
        s["mz_lo"] = min(w["IsolationMz"] - w["IsolationWidth"] / 2 for w in wins)
        s["mz_hi"] = max(w["IsolationMz"] + w["IsolationWidth"] / 2 for w in wins)
        centres = sorted({round(float(w["IsolationMz"]), 4) for w in wins})
        steps = {round(b - a, 3) for a, b in zip(centres, centres[1:])}
        if len(steps) == 1 and len(widths) == 1:
            s["spacing"] = steps.pop()
            s["overlap"] = round(widths[0] - s["spacing"], 3)
    ces = [float(w["CollisionEnergy"]) for w in wins if w.get("CollisionEnergy") is not None]
    if ces:
        s["ce_range"] = (min(ces), max(ces))
    wim = v.get("window_im") or []
    if wim:
        s["im_lo"] = min(w["im_lo"] for w in wim)
        s["im_hi"] = max(w["im_hi"] for w in wim)
    return s


def ce_at(points, im):
    """The collision energy a piecewise-linear ramp gives at 1/K0 `im` (held flat outside)."""
    pts = sorted(points)
    if im <= pts[0][0]:
        return pts[0][1]
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        if im <= x1:
            return y0 + (y1 - y0) * (im - x0) / (x1 - x0) if x1 != x0 else y1
    return pts[-1][1]


def check_ce_ramp(v):
    """Does each window's recorded collision energy match the ramp at the window's 1/K0
    midpoint? 'ok' / 'mismatch' (with the worst difference) / 'unchecked' (no 1/K0 bounds)."""
    ramp, wins, wim = v.get("ce_ramp"), v.get("windows") or [], v.get("window_im") or []
    if not ramp or not wins or not wim:
        return {"status": "unchecked"}
    bounds = {(w["group"], round(w["mz"], 3)): (w["im_lo"] + w["im_hi"]) / 2 for w in wim}
    diffs = []
    for w in wins:
        key = (w.get("WindowGroup"), round(float(w.get("IsolationMz") or 0), 3))
        if key in bounds and w.get("CollisionEnergy") is not None:
            diffs.append(abs(float(w["CollisionEnergy"]) - ce_at(ramp["points"], bounds[key])))
    if not diffs:
        return {"status": "unchecked"}
    worst = round(max(diffs), 2)
    return {"status": "ok" if worst <= CE_CHECK_EV else "mismatch", "max_diff_ev": worst,
            "n_windows": len(diffs)}


# What must match across a series for one Methods paragraph to describe every run.
SERIES_KEYS = ("instrument", "instrument_serial", "ms_method", "lc_system", "lc_method",
               "lc_run_min", "mz_low", "mz_high", "im_low", "im_high", "ramp_ms",
               "accumulation_ms", "column_info", "source_type")


def _when(s):
    try:
        return datetime.fromisoformat(s)
    except (TypeError, ValueError):
        return None


def series(runs):
    """Summarise a list of read_run() results: the differences that one paragraph would hide,
    and the acquisition window (first/last AcquisitionDateTime)."""
    diffs = {}
    for key in SERIES_KEYS:
        vals = {str(r["values"].get(key)) for r in runs if r["values"].get(key) is not None}
        if len(vals) > 1:
            diffs[key] = sorted(vals)
    schemes = {repr(window_scheme(r["values"])) for r in runs if r["values"].get("windows")}
    if len(schemes) > 1:
        diffs["dia-PASEF window scheme"] = [f"{len(schemes)} different schemes"]
    stamps = sorted((w, r["values"]["acquired_at"]) for r in runs
                    if (w := _when(r["values"].get("acquired_at"))) is not None)
    return {"differences": diffs,
            "acquired_first": stamps[0][1] if stamps else None,
            "acquired_last": stamps[-1][1] if stamps else None}
