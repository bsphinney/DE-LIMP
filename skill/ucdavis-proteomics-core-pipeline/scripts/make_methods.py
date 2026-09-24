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
      [--lc-column "PepSep C18, 10 cm × 150 µm, 1.5 µm"] \
      [--params wf/params.cfg --search-prov search/search_provenance.json \
       --workflow-manifest wf/workflow.manifest.json]   # adds the Database-search section
      [--de-dir output/tables]      # optional: adds a Differential-expression paragraph
      [--instrument "timsTOF HT" --acquisition DIA]   # used only when the raw files
                                                      # cannot be read from here
"""
import sys, os, json, glob, sqlite3, argparse, statistics

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
# analysis.tdf is opened read-only AND immutable -- see bruker_tdf.py for how a read-write
# open truncates a tdf (the state of 342 on HIVE), and why mode=ro alone is not enough.
from bruker_tdf import connect_tdf  # noqa: E402

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
# A value this script could not find in any record of the run. Printed in place of the value,
# never replaced by a plausible default (DE-LIMP rule #2).
NOT_RECORDED = "____ [not recorded — confirm]"

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
    """Extract acquisition parameters from a Bruker .d analysis.tdf (best-effort)."""
    tdf = os.path.join(d, "analysis.tdf")
    if not os.path.exists(tdf):
        return None
    m = {"vendor": "Bruker", "file": os.path.basename(d.rstrip("/"))}
    try:
        con = connect_tdf(tdf)
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


def _load_json(path):
    if not path or not os.path.isfile(path):
        return None
    try:
        with open(path) as fh:
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
           "tol_note": None, "precursor_fdr": None, "library": None, "mbr": None,
           "labelled": None, "warnings": []}
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
        from diann_parallel import cfg_tokens, cfg_groups, CfgError
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
        if "--reanalyse" in flags or "--reanalyse" in cmd_words:
            rec["mbr"] = {"value": True, "source": where if "--reanalyse" in flags else
                          "search_provenance.json resolved_command"}
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
    if rec.get("mbr"):
        s[0] += ", with match-between-runs"
    s[0] += "."
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
            tol.append(f"{lvl} {sym}{_g(t['value'])} {t['unit']}")
    if tol:
        s.append("Mass tolerances were " + " and ".join(tol) + ".")
    elif rec.get("tol_note"):
        s.append(rec["tol_note"][0].upper() + rec["tol_note"][1:] + ".")
    fdr = rec.get("precursor_fdr")
    if fdr:
        s.append(f"Precursor identifications were filtered at {fdr['value'] * 100:g}% FDR "
                 f"(q ≤ {fdr['value']:g}).")
    elif not (de_prov or {}).get("q_columns"):
        s.append(f"Identification FDR: {NOT_RECORDED}.")
    for w in rec.get("warnings") or []:
        s.append(f"[{w} — confirm]")
    return " ".join(s)


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
    if prov.get("rollup_method"):
        s.append(f"Protein quantities: {prov['rollup_method']}.")
    if prov.get("missing_policy"):
        s.append(prov["missing_policy"].rstrip(".") + ".")
    cols, cuts = prov.get("q_columns") or [], prov.get("q_cutoffs") or []
    if cols and len(cols) == len(cuts):
        s.append("Identifications entering quantification were filtered at "
                 + ", ".join(f"{c} ≤ {x:g}" for c, x in zip(cols, cuts)) + ".")
    elif prov.get("q_cutoff") is not None:
        s.append(f"Identifications entering quantification were filtered at q ≤ "
                 f"{prov['q_cutoff']:g}.")
    else:
        s.append(f"Identification q-value filter: {NOT_RECORDED}.")
    if prov.get("design"):
        s.append(f"The linear model was {prov['design']}"
                 + (f", with contrasts {', '.join(prov['contrasts'])}"
                    if prov.get("contrasts") else "") + ".")
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


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw", nargs="+", required=True, help="raw file paths/globs (.d or .raw)")
    ap.add_argument("--out", default="methods.md")
    ap.add_argument("--lc-column", default=LC_COLUMN_DEFAULT)
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
    a = ap.parse_args()

    fmeta = None
    if a.fasta_meta:
        try:
            fmeta = json.load(open(a.fasta_meta))
        except (OSError, json.JSONDecodeError) as e:
            sys.exit(f"--fasta-meta could not be read: {e}")

    files = []
    for p in a.raw:
        files.extend(sorted(glob.glob(p)) or [p])
    metas = detect(files)
    for m in metas:
        m.setdefault("instrument_source", "GlobalMetadata InstrumentName" if m.get("vendor") ==
                     "Bruker" else "facility filename prefix")
    rec_instr = (a.instrument or "").strip() or None
    rec_src = "session record (workflow manifest) — not read from the raw file"
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
    bru = [m for m in metas if m.get("vendor") == "Bruker" and m.get("instrument")]
    rep = bru[0] if bru else metas[0]
    if not rep.get("instrument") and rec_instr:
        rep = dict(rep, instrument=rec_instr, instrument_source=rec_src)
    instrument = rep.get("instrument") or next((m.get("instrument") for m in metas if m.get("instrument")), None)
    ack_label, ack_text = pick_ack(instrument, files)

    srec = search_record(a.params, a.search_prov, a.workflow_manifest)
    de_prov = _load_json(os.path.join(a.de_dir, "de_provenance.json")) if a.de_dir else None

    json.dump({"files": [m.get("file") for m in metas], "representative": rep,
               "instrument": instrument, "acknowledgment_for": ack_label, "all": metas,
               "from_session_record": from_record, "acquisition": a.acquisition,
               "search": srec},
              open(os.path.splitext(a.out)[0] + "_params.json", "w"), indent=2)

    # A blank acquisition value is a facility default to confirm when the raw file was read, but
    # simply unknown when it could not be -- it must not then be labelled a facility default.
    blank_tag = "[raw file not readable here — confirm]" if from_record else DEF

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
          f"acquisition values below are blank and tagged, and the instrument and acquisition "
          f"mode come from the workflow manifest. Values marked {DEF} are facility defaults to "
          f"confirm; values marked {NOT_RECORDED} were not in any record. Re-run "
          f"make_methods.py where the raw files are readable to fill them in.*")
    else:
        w(f"*Generated by the UC Davis Proteomics Core pipeline skill from the raw data "
          f"({len(metas)} file(s)). Values marked {DEF} are facility defaults to confirm; "
          "all other values were extracted from the raw acquisition metadata"
          + (" and the search and analysis records" if (srec.get("engine") or de_prov)
             else "") + ".*")
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
        acq = (a.acquisition or "").upper()
        mode = (f"{acq} mode (as detected from the data in step 2)" if acq in ("DIA", "DDA")
                else f"[DDA/DIA — confirm] mode {DEF}")
        w(f"Mass spectra were acquired on a {v(rep.get('instrument'), default='[instrument]')} mass "
          f"spectrometer (Thermo Fisher Scientific) operated in {mode}. "
          "Full acquisition parameters (resolution, AGC, isolation width, NCE, gradient) should be "
          f"taken from the instrument method file {DEF}.")
    w("")

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
        rel = f"release {rel}" if (rel := fmeta.get("uniprot_release")) else f"release ____ {DEF}"
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
                    f"release ____ {DEF}), comprising {n_p} sequences")
            g = (fmeta.get("content_check") or {}).get("uniprot_gene_count")
            if fmeta.get("content_inferred") == "one_per_gene" and isinstance(g, int):
                sent += (f"; the entry count is consistent with one canonical protein "
                         f"sequence per gene (inferred, not verified: UniProt lists "
                         f"{g:,} genes).")
            else:
                sent += f". Database composition: ____ {DEF}."
        elif content_phrase is None:
            # 'unknown' (--path) / 'as_staged' (--hive): we did not build this database,
            # so we cannot describe its composition. Leave it tagged for the user.
            sent = (f"Spectra were searched against a supplied sequence database "
                    f"({os.path.basename(fmeta.get('fasta', '') ) or '____'}; "
                    f"{n_p} sequences). Database composition and version: ____ {DEF}.")
        else:
            sent = (f"Spectra were searched against the UniProt "
                    f"{fmeta.get('organism') or '____'} {kind} "
                    f"({fmeta.get('proteome') or '____'}, {rel}), comprising "
                    f"{content_phrase} ({n_p} sequences).")
        n_c = fmeta.get("n_contaminants_appended") or 0
        n_already = fmeta.get("n_contaminants_already_present") or 0
        if not n_c and n_already:
            sent += (f" The database already included {n_already} common-contaminant "
                     f"sequences")
            sent += (" and these entries were excluded from quantification and "
                     "normalisation." if fmeta.get("diann_cont_quant_exclude") else ".")
        elif n_c:
            sent += (f" A common-contaminant library ({n_c} sequences; "
                     f"{fmeta.get('contaminant_set')} set of Frankenfield et al., "
                     f"J Proteome Res 2022, 21:2104-2113) was appended")
            sent += (" and these entries were excluded from quantification and "
                     "normalisation."
                     if fmeta.get("diann_cont_quant_exclude") else ".")
            # fetch_fasta.py removes contaminant entries whose sequence IS a target protein
            # (bovine ACTB = human ACTB, human keratins); a reader must know those proteins
            # were quantified, not excluded as contaminants.
            n_drop = fmeta.get("n_contaminants_dropped_as_target") or 0
            if n_drop:
                sent += (f" {n_drop} contaminant entries identical to (or contained in) "
                         f"{fmeta.get('organism') or '____'} proteins were removed from the "
                         f"library first, so those proteins are quantified under their own "
                         f"accessions.")
        else:
            sent += " No contaminant database was appended."
        w(sent)
        # The drop note is described in the sentence above -- it is a record, not
        # something to resolve before publication.
        build_warnings = [x for x in (fmeta.get("warnings") or [])
                          if x != fmeta.get("contaminants_dropped_note")]
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

    # parameter table (value + source)
    w("## Acquisition parameters (extracted from the raw data)" if not from_record else
      "## Acquisition parameters (from the session record — raw files not readable here)")
    w("")
    w("| Parameter | Value | Source |")
    w("|---|---|---|")
    rows = [("Instrument", rep.get("instrument"), rep.get("instrument_source")),
            ("Acquisition software", rep.get("software"), "GlobalMetadata"),
            ("Acquisition mode", rep.get("mode"), rec_src if from_record else "Frames MsMsType"),
            ("m/z range", f"{rep.get('mz_low')}–{rep.get('mz_high')}" if rep.get("mz_low") else None, "GlobalMetadata MzAcqRange*"),
            ("1/K₀ range (V·s/cm²)", f"{rep.get('im_low')}–{rep.get('im_high')}" if rep.get("im_low") else None, "GlobalMetadata OneOverK0AcqRange*"),
            ("TIMS ramp / accumulation (ms)", f"{rep.get('ramp_ms')} / {rep.get('accumulation_ms')}" if rep.get("ramp_ms") else None, "Frames RampTime/AccumulationTime"),
            ("Isolation windows", rep.get("n_windows"), "DiaFrameMsMsWindows"),
            ("Isolation width (Th)", rep.get("isolation_width"), "DiaFrameMsMsWindows IsolationWidth"),
            ("Collision energy (eV)", f"{rep.get('ce_low')}–{rep.get('ce_high')}" if rep.get("ce_low") is not None else None, "DiaFrameMsMsWindows CollisionEnergy"),
            ("Analytical column", a.lc_column, "facility default — confirm"
             if a.lc_column == LC_COLUMN_DEFAULT else "--lc-column (user-given)"),
            ("Files in series", len(metas), "this run")]
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
            srows.append((label, f"{_g(t['value'])} {t['unit']}" if t else
                          (srec.get("tol_note") or NOT_RECORDED),
                          t["source"] if t else "search record"))
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

    open(a.out, "w").write("\n".join(L) + "\n")
    print(json.dumps({"methods": os.path.abspath(a.out), "instrument": instrument,
                      "acknowledgment_for": ack_label, "n_files": len(metas),
                      "params_json": os.path.splitext(a.out)[0] + "_params.json",
                      "next": "Verify the draft against the params table, polish the prose, then "
                              "convert to .docx with to_docx.py."}, indent=2))


if __name__ == "__main__":
    main()
