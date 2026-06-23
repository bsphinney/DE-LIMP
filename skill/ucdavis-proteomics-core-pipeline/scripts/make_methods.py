#!/usr/bin/env python3
"""
make_methods.py  --  Generate a publication-ready LC-MS/MS Methods section from
facility raw data, plus the correct UC Davis Proteomics Core instrument-grant
acknowledgment.

It reads what it can directly from the raw metadata (Bruker .d analysis.tdf;
Thermo .raw by facility filename prefix / reader) and fills the rest from facility
defaults that are CLEARLY TAGGED `[facility default — confirm]` so nothing is
silently fabricated (DE-LIMP rule #2). The default LC column is a PepSep C18
10 cm × 150 µm, 1.5 µm column (override with --lc-column). It writes:

  methods.md          drop-in Methods prose (LC, MS, Data processing) + a
                      parameter table (value + where each value came from) +
                      an instrument-specific Acknowledgments section
  methods_params.json the extracted parameters, machine-readable

Then to_docx.py can render methods.md to Word. The agent should verify the draft
against the extracted params and polish the prose (keep the acknowledgment exact).

Acknowledgments are from https://proteomics.ucdavis.edu/instrument-grant-acknowledgments
(verified 2026-06). Confirm exact wording there before publishing.

Usage:
  python3 make_methods.py --raw '/data/*.d' --out methods.md \
      [--lc-column "PepSep C18, 10 cm × 150 µm, 1.5 µm"] \
      [--de-dir output/tables]      # optional: adds a Data-processing paragraph
"""
import sys, os, json, glob, sqlite3, argparse, statistics

ACK_SOURCE = "https://proteomics.ucdavis.edu/instrument-grant-acknowledgments"
# (instrument-name substrings, facility filename prefixes, label, acknowledgment).
# Verified against the UC Davis Proteomics Core grant-acknowledgment page (2026-06).
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
LC_COLUMN_DEFAULT = "PepSep C18, 10 cm × 150 µm i.d., 1.5 µm reversed-phase particles (Bruker/Dr. Maisch)"
DEF = "[facility default — confirm]"


def _num(x):
    try: return float(x)
    except (TypeError, ValueError): return None


def bruker_meta(d):
    """Extract acquisition parameters from a Bruker .d analysis.tdf (best-effort)."""
    tdf = os.path.join(d, "analysis.tdf")
    if not os.path.exists(tdf):
        return None
    m = {"vendor": "Bruker", "file": os.path.basename(d.rstrip("/"))}
    try:
        con = sqlite3.connect(f"file:{tdf}?mode=ro", uri=True)
        cur = con.cursor()
        gm = dict(cur.execute("SELECT Key, Value FROM GlobalMetadata"))
        m["instrument"] = gm.get("InstrumentName")
        sw = gm.get("AcquisitionSoftware", "")
        ver = gm.get("AcquisitionSoftwareVersion", "")
        m["software"] = (sw + (" " + ver if ver else "")).strip() or None
        m["mz_low"], m["mz_high"] = _num(gm.get("MzAcqRangeLower")), _num(gm.get("MzAcqRangeUpper"))
        m["im_low"], m["im_high"] = _num(gm.get("OneOverK0AcqRangeLower")), _num(gm.get("OneOverK0AcqRangeUpper"))
        types = dict(cur.execute("SELECT MsMsType, COUNT(*) FROM Frames GROUP BY MsMsType"))
        m["mode"] = "dia-PASEF" if types.get(9) else ("ddaPASEF" if types.get(8) else "MS")
        row = cur.execute("SELECT AccumulationTime, RampTime FROM Frames WHERE MsMsType IN (8,9) LIMIT 1").fetchone()
        if row:
            m["accumulation_ms"], m["ramp_ms"] = _num(row[0]), _num(row[1])
        tbls = {r[0] for r in cur.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        if "DiaFrameMsMsWindows" in tbls:
            widths = [r[0] for r in cur.execute("SELECT IsolationWidth FROM DiaFrameMsMsWindows") if r[0] is not None]
            ces = [r[0] for r in cur.execute("SELECT CollisionEnergy FROM DiaFrameMsMsWindows") if r[0] is not None]
            n = cur.execute("SELECT COUNT(*) FROM DiaFrameMsMsWindows").fetchone()[0]
            grps = cur.execute("SELECT COUNT(DISTINCT WindowGroup) FROM DiaFrameMsMsWindows").fetchone()[0] \
                if "WindowGroup" in [c[1] for c in cur.execute("PRAGMA table_info(DiaFrameMsMsWindows)")] else None
            m["n_windows"] = n; m["n_window_groups"] = grps
            if widths: m["isolation_width"] = round(statistics.median(widths), 1)
            if ces: m["ce_low"], m["ce_high"] = round(min(ces), 1), round(max(ces), 1)
        con.close()
    except sqlite3.Error as e:
        m["error"] = str(e)
    return m


def thermo_meta(f):
    """Thermo .raw: identify by facility filename prefix (FL*, Ex*) — the model is
    not reliably readable without a vendor reader."""
    base = os.path.basename(f)
    m = {"vendor": "Thermo", "file": base, "mode": None}
    for subs, prefixes, label, _ in ACKS:
        if any(base.startswith(p) for p in prefixes):
            m["instrument"] = label
            break
    return m


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
    instr = (instrument or "").lower()
    bn = [os.path.basename(f) for f in files]
    for subs, prefixes, label, text in ACKS:
        if any(s in instr for s in subs) or any(b.startswith(p) for b in bn for p in prefixes):
            return label, text
    return None, (f"[Instrument not in the UC Davis acknowledgment registry — "
                  f"check {ACK_SOURCE} and insert the correct instrument-grant acknowledgment.]")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw", nargs="+", required=True, help="raw file paths/globs (.d or .raw)")
    ap.add_argument("--out", default="methods.md")
    ap.add_argument("--lc-column", default=LC_COLUMN_DEFAULT)
    ap.add_argument("--de-dir", help="optional: de_provenance.json/methods.txt for a Data-processing paragraph")
    a = ap.parse_args()

    files = []
    for p in a.raw:
        files.extend(sorted(glob.glob(p)) or [p])
    metas = detect(files)
    if not metas:
        sys.exit("No raw files found.")

    # representative metadata (facility usually acquires a series identically)
    bru = [m for m in metas if m.get("vendor") == "Bruker" and m.get("instrument")]
    rep = bru[0] if bru else metas[0]
    instrument = rep.get("instrument") or next((m.get("instrument") for m in metas if m.get("instrument")), None)
    ack_label, ack_text = pick_ack(instrument, files)

    json.dump({"files": [m.get("file") for m in metas], "representative": rep,
               "instrument": instrument, "acknowledgment_for": ack_label, "all": metas},
              open(os.path.splitext(a.out)[0] + "_params.json", "w"), indent=2)

    def v(x, unit="", default=None):
        if x is None:
            return f"{default} {DEF}" if default is not None else f"____ {DEF}"
        return f"{x}{unit}"

    is_bruker = rep.get("vendor") == "Bruker"
    L, w = [], lambda s="": L.append(s)

    w("# Materials and Methods — LC-MS/MS")
    w("")
    w(f"*Generated by the UC Davis Proteomics Core pipeline skill from the raw data "
      f"({len(metas)} file(s)). Values marked {DEF} are facility defaults to confirm; "
      "all other values were extracted from the raw acquisition metadata.*")
    w("")
    w("## Liquid chromatography")
    w("")
    w(f"Peptides were separated by reversed-phase nano-LC on a {a.lc_column} "
      f"{DEF if a.lc_column == LC_COLUMN_DEFAULT else ''}, using water containing 0.1% "
      "(v/v) formic acid as mobile phase A and acetonitrile containing 0.1% (v/v) "
      f"formic acid as mobile phase B {DEF}. "
      + ("The column was interfaced to the mass spectrometer through a Bruker "
         f"CaptiveSpray source with a 20 µm i.d. PepSep emitter {DEF}. "
         if is_bruker else
         f"The column was interfaced to the mass spectrometer by a nanospray source {DEF}. ")
      + f"The LC system and gradient were [LC system / gradient — confirm] {DEF}.")
    w("")
    w("## Mass spectrometry")
    w("")
    if is_bruker:
        w(f"Mass spectra were acquired on a {v(rep.get('instrument'))} mass spectrometer "
          f"(Bruker Daltonics)" + (f", operated with {rep['software']}" if rep.get("software") else "")
          + f" in positive-ion {v(rep.get('mode'))} mode. "
          f"Spectra were recorded over m/z {v(rep.get('mz_low'))}–{v(rep.get('mz_high'))}, "
          f"and the trapped-ion-mobility analyzer was scanned over 1/K₀ = "
          f"{v(rep.get('im_low'))}–{v(rep.get('im_high'))} V·s/cm²"
          + (f", with a TIMS ramp/accumulation time of {v(rep.get('ramp_ms'))}/{v(rep.get('accumulation_ms'))} ms"
             if rep.get("ramp_ms") else "") + ".")
        if rep.get("n_windows"):
            w("")
            w(f"The {v(rep.get('mode'))} method used {v(rep.get('n_windows'))} isolation windows"
              + (f" across {rep['n_window_groups']} window groups" if rep.get("n_window_groups") else "")
              + (f" (≈{rep['isolation_width']} Th wide)" if rep.get("isolation_width") else "")
              + (f", with collision energy ramped from ≈{rep['ce_low']} to ≈{rep['ce_high']} eV with ion mobility"
                 if rep.get("ce_low") is not None else "") + ".")
    else:
        w(f"Mass spectra were acquired on a {v(rep.get('instrument'), default='[instrument]')} mass "
          f"spectrometer (Thermo Fisher Scientific) operated in [DDA/DIA — confirm] mode {DEF}. "
          "Full acquisition parameters (resolution, AGC, isolation width, NCE, gradient) should be "
          f"taken from the instrument method file {DEF}.")
    w("")

    # optional data-processing paragraph from the skill's own run
    if a.de_dir and os.path.exists(os.path.join(a.de_dir, "de_provenance.json")):
        prov = json.load(open(os.path.join(a.de_dir, "de_provenance.json")))
        w("## Data processing")
        w("")
        w(f"Raw files were searched and quantified with {prov.get('display_label', 'the configured pipeline')} "
          f"({prov.get('rollup_method','')}; {prov.get('de_engine','')}). Identifications were filtered to "
          f"{prov.get('q_cutoff', 0.01)*100:.0f}% FDR; differential expression used the thresholds "
          f"adj.P.Val < {prov.get('adjp', 0.05)} and |log2FC| ≥ {prov.get('logfc', 1)}. "
          f"{prov.get('missing_policy','')} {prov.get('citation','')}")
        w("")

    # parameter table (value + source)
    w("## Acquisition parameters (extracted from the raw data)")
    w("")
    w("| Parameter | Value | Source |")
    w("|---|---|---|")
    rows = [("Instrument", rep.get("instrument"), "GlobalMetadata InstrumentName" if is_bruker else "filename prefix"),
            ("Acquisition software", rep.get("software"), "GlobalMetadata"),
            ("Acquisition mode", rep.get("mode"), "Frames MsMsType"),
            ("m/z range", f"{rep.get('mz_low')}–{rep.get('mz_high')}" if rep.get("mz_low") else None, "GlobalMetadata MzAcqRange*"),
            ("1/K₀ range (V·s/cm²)", f"{rep.get('im_low')}–{rep.get('im_high')}" if rep.get("im_low") else None, "GlobalMetadata OneOverK0AcqRange*"),
            ("TIMS ramp / accumulation (ms)", f"{rep.get('ramp_ms')} / {rep.get('accumulation_ms')}" if rep.get("ramp_ms") else None, "Frames RampTime/AccumulationTime"),
            ("Isolation windows", rep.get("n_windows"), "DiaFrameMsMsWindows"),
            ("Isolation width (Th)", rep.get("isolation_width"), "DiaFrameMsMsWindows IsolationWidth"),
            ("Collision energy (eV)", f"{rep.get('ce_low')}–{rep.get('ce_high')}" if rep.get("ce_low") is not None else None, "DiaFrameMsMsWindows CollisionEnergy"),
            ("Analytical column", a.lc_column, "facility default — confirm"),
            ("Files in series", len(metas), "this run")]
    for name, val, src in rows:
        if val is None: continue
        w(f"| {name} | {val} | {src} |")
    w("")

    w("## Acknowledgments")
    w("")
    w(ack_text)
    w("")
    w(f"*Acknowledgment source: {ACK_SOURCE} (confirm the exact current wording before publishing).*")
    w("")

    open(a.out, "w").write("\n".join(L) + "\n")
    print(json.dumps({"methods": os.path.abspath(a.out), "instrument": instrument,
                      "acknowledgment_for": ack_label, "n_files": len(metas),
                      "params_json": os.path.splitext(a.out)[0] + "_params.json",
                      "next": "Verify the draft against the params table, polish the prose, then "
                              "convert to .docx with to_docx.py."}, indent=2))


if __name__ == "__main__":
    main()
