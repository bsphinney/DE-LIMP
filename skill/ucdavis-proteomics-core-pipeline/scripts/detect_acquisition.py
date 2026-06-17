#!/usr/bin/env python3
"""
detect_acquisition.py  --  Decide DDA vs DIA for each input file, so the
pipeline can auto-route to the right engine (DIA -> DIA-NN, DDA -> Sage by
default; FragPipe selectable for either).

Detection per format (best-effort, with a confidence score):
  Bruker .d   read analysis.tdf (SQLite). Presence of a DiaFrameMsMsInfo /
              DiaFrameMsMsWindowGroups table, or Frames.MsMsType == 9,
              => dia-PASEF. PasefFrameMsMsInfo / MsMsType == 8 => ddaPASEF.
              (Mirrors DE-LIMP's helpers_instrument.R logic.)
  mzML(.gz)   stream isolation windows of MS2 scans. Wide windows (median
              >= 3 Da) over a small repeating set of centers => DIA; narrow
              windows (<= 2 Da) with many distinct centers => DDA.
  Thermo .raw not natively readable here. If ThermoRawFileParser is on PATH we
              convert a header sample; otherwise we return 'unknown' and ask.

Output: JSON to stdout. The orchestrator MUST confirm with the user whenever
confidence != "high" before launching a multi-hour search.

Usage: python3 detect_acquisition.py FILE [FILE ...]
"""
import sys, os, json, glob, gzip, sqlite3, statistics, shutil, subprocess

# ---------------------------------------------------------------------------
# Instrument extraction (PLAN.md §7d) — mirrors DE-LIMP R/helpers_instrument.R.
# The matcher in fetch_workflows.py scores on instrument; null is allowed and
# falls back to score 0 (user confirms). Best-effort, never fatal.
# ---------------------------------------------------------------------------
def instrument_bruker_d(path):
    tdf = os.path.join(path, "analysis.tdf")
    if not os.path.exists(tdf):
        return None
    try:
        con = sqlite3.connect(f"file:{tdf}?mode=ro", uri=True)
        cur = con.cursor()
        # GlobalMetadata is a key/value table; InstrumentName holds e.g. "timsTOF Pro"
        rows = dict(cur.execute("SELECT Key, Value FROM GlobalMetadata"))
        con.close()
        for k in ("InstrumentName", "InstrumentSourceType", "Instrument"):
            if rows.get(k):
                return rows[k]
    except sqlite3.Error:
        return None
    return None

def instrument_thermo_raw(path):
    # Thermo .raw model string requires a reader. Try ThermoRawFileParser metadata.
    parser = shutil.which("ThermoRawFileParser") or shutil.which("thermorawfileparser")
    if not parser:
        return None
    try:
        out = subprocess.run([parser, "metadata", "-i", path],
                             capture_output=True, text=True, timeout=120).stdout
        for line in out.splitlines():
            low = line.lower()
            if "instrument model" in low or "instrument name" in low:
                return line.split(":", 1)[-1].strip() or None
    except Exception:
        return None
    return None

def detect_instrument(path):
    low = path.lower().rstrip("/")
    if low.endswith(".d"):
        return instrument_bruker_d(path)
    if low.endswith(".raw"):
        return instrument_thermo_raw(path)
    return None

def detect_bruker_d(path):
    tdf = os.path.join(path, "analysis.tdf")
    if not os.path.exists(tdf):
        return ("unknown", "low", "no analysis.tdf in .d folder")
    tmp = tdf  # read-only open
    try:
        con = sqlite3.connect(f"file:{tmp}?mode=ro", uri=True)
        cur = con.cursor()
        tables = {r[0] for r in cur.execute(
            "SELECT name FROM sqlite_master WHERE type='table'")}
        if {"DiaFrameMsMsInfo", "DiaFrameMsMsWindowGroups"} & tables:
            return ("DIA", "high", "dia-PASEF window tables present")
        if "PasefFrameMsMsInfo" in tables:
            return ("DDA", "high", "PasefFrameMsMsInfo present (ddaPASEF)")
        # fall back to MsMsType histogram
        try:
            rows = dict(cur.execute(
                "SELECT MsMsType, COUNT(*) FROM Frames GROUP BY MsMsType"))
            if rows.get(9, 0) > 0:  return ("DIA", "high", "Frames.MsMsType==9")
            if rows.get(8, 0) > 0:  return ("DDA", "medium", "Frames.MsMsType==8")
        except sqlite3.Error:
            pass
        return ("unknown", "low", "no DIA/DDA markers in tdf")
    except sqlite3.Error as e:
        return ("unknown", "low", f"tdf read error: {e}")
    finally:
        try: con.close()
        except Exception: pass

def _iter_mzml(path, cap=3000):
    """Yield (ms_level, iso_target, iso_width) for spectra, streaming."""
    import xml.etree.ElementTree as ET
    opn = gzip.open if path.endswith(".gz") else open
    ms_level = iso_target = lower = upper = None
    n = 0
    with opn(path, "rb") as fh:
        for ev, el in ET.iterparse(fh, events=("end",)):
            tag = el.tag.rsplit("}", 1)[-1]
            if tag == "cvParam":
                acc = el.get("accession"); val = el.get("value")
                if acc == "MS:1000511": ms_level = int(float(val))
                elif acc == "MS:1000827" and val: iso_target = float(val)
                elif acc == "MS:1000828" and val: lower = float(val)
                elif acc == "MS:1000829" and val: upper = float(val)
            elif tag == "spectrum":
                if ms_level == 2:
                    w = (lower or 0) + (upper or 0)
                    yield (2, iso_target, w if w else None)
                    n += 1
                ms_level = iso_target = lower = upper = None
                el.clear()
                if n >= cap: break

def detect_mzml(path):
    widths, targets = [], []
    try:
        for lvl, tgt, w in _iter_mzml(path):
            if w: widths.append(w)
            if tgt is not None: targets.append(round(tgt * 2) / 2)
    except Exception as e:
        return ("unknown", "low", f"mzML parse error: {e}")
    if not widths:
        return ("unknown", "low", "no MS2 isolation windows found")
    med = statistics.median(widths)
    distinct = len(set(targets))
    n = len(widths)
    if med >= 3.0 and distinct <= max(80, n // 50):
        return ("DIA", "high" if med >= 4 else "medium",
                f"median isolation width {med:.1f} Da over {distinct} distinct centers")
    if med <= 2.0:
        return ("DDA", "high" if distinct > n // 10 else "medium",
                f"median isolation width {med:.1f} Da, {distinct} distinct precursors")
    return ("unknown", "low", f"ambiguous: median width {med:.1f} Da, {distinct} centers")

def detect_thermo_raw(path):
    parser = shutil.which("ThermoRawFileParser") or shutil.which("thermorawfileparser")
    if not parser:
        return ("unknown", "low",
                "Thermo .raw needs ThermoRawFileParser (-> mzML) or DIA-NN's own reader to classify; confirm manually")
    try:
        out = subprocess.run([parser, "query", "-i", path], capture_output=True,
                             text=True, timeout=120).stdout.lower()
        if "dia" in out: return ("DIA", "medium", "ThermoRawFileParser metadata mentions DIA")
        if "dda" in out: return ("DDA", "medium", "ThermoRawFileParser metadata mentions DDA")
    except Exception as e:
        return ("unknown", "low", f"ThermoRawFileParser error: {e}")
    return ("unknown", "low", "could not classify from header")

def classify(path):
    p = path.rstrip("/")
    low = p.lower()
    if low.endswith(".d"):
        kind, conf, why = detect_bruker_d(p); vendor = "Bruker"
    elif low.endswith((".mzml", ".mzml.gz")):
        kind, conf, why = detect_mzml(p); vendor = "mzML"
    elif low.endswith(".raw"):
        kind, conf, why = detect_thermo_raw(p); vendor = "Thermo"
    elif low.endswith((".wiff", ".wiff2")):
        kind, conf, why, vendor = "unknown", "low", "SCIEX .wiff: convert to mzML to classify", "SCIEX"
    else:
        kind, conf, why, vendor = "unknown", "low", "unrecognized extension", "?"
    instrument = detect_instrument(p)
    return {"file": p, "vendor": vendor, "acquisition": kind,
            "confidence": conf, "reason": why, "instrument": instrument}

def main(argv):
    files = []
    for a in argv:
        files.extend(sorted(glob.glob(a)) or [a])
    results = [classify(f) for f in files]
    kinds = {r["acquisition"] for r in results}
    overall = (next(iter(kinds)) if len(kinds) == 1 else "mixed")
    low_conf = [r["file"] for r in results if r["confidence"] != "high"]
    instruments = sorted({r["instrument"] for r in results if r.get("instrument")})
    # one instrument string for the matcher: unambiguous only if all files agree
    instrument = instruments[0] if len(instruments) == 1 else ""
    print(json.dumps({
        "overall": overall,
        "instrument": instrument,
        "instruments_seen": instruments,
        "needs_confirmation": bool(low_conf) or overall in ("mixed", "unknown") or len(instruments) > 1,
        "low_confidence_files": low_conf,
        "files": results,
    }, indent=2))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: detect_acquisition.py FILE [FILE ...]", file=sys.stderr); sys.exit(2)
    main(sys.argv[1:])
