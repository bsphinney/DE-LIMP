#!/usr/bin/env python3
"""
session_docs.py  --  The session's README (README.md and README.html, rendered from ONE source)
and AGENTS.md, written by `session.py finalize` (or `session.py docs`) from the session's own
records -- and input/raw_files.txt when a session does not have one.

  README.html  what collaborators open: a self-contained page (inline CSS, no external assets)
               that double-clicks open on Windows and Macs, with links to the report, the Word
               files, the deposit guide and the tables
  README.md    the same content as text; the DataAnalysis rules and record_run.py copy it
  AGENTS.md    for an AI agent handed the whole folder: the study, which file is authoritative
               for what, the table columns, the traps, and what it must not do

Everything comes from files in the session: search_provenance.json, de_provenance.json,
methods.txt, the FASTA sidecar, conditions.csv, figures.json, AUDIT / SAMPLE_QUALITY,
MANIFEST.txt, fran_deposit.json, the attached CoreOmics submission (submission_report.py). The pipeline is described in de_provenance.json's own words --
name, roll-up, missing-value policy, significance rule -- never from memory here (CLAUDE.md rule
1), and a location that cannot be determined is "not recorded", never a guess (rule 2). Where a
path is on HIVE, and how a collaborator browses to it from Windows or a Mac, is share_map.py
(hive_shares.tsv, the table hive_path.sh reads).
"""
import ast
import csv
import datetime
import glob
import json
import os
import re
import sys
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
from session import paths_for, read_raw_list      # noqa: E402  the session layout, one place
import share_map                                    # noqa: E402

NOT_RECORDED = "not recorded"
RAW_HEADER = "# Raw MS files used in this analysis (not copied — too large)."

# DE_*.csv / Expression_Matrix.csv columns. Sources: limma's topTable() (logFC ... B), limpa's
# dpcQuant() documentation (NPeptides, PropObs; limpa 1.5.0 man/dpcQuant.Rd) and dpcQuant.R (an
# annotation column survives to the protein row only when all its precursors agree), DIA-NN's
# README "Main output reference" (Protein.Group, Protein.Names, Genes, Proteotypic).
COLUMNS = {
    "Protein.Group": "the protein group DIA-NN inferred (UniProt accession(s)); the row key in "
                     "every table",
    "Protein.Names": "UniProt names of the proteins in the group (DIA-NN)",
    "Genes": "gene names of the proteins in the group (DIA-NN)",
    "Proteotypic": "DIA-NN's proteotypic flag (0/1). limpa keeps it on the protein row only when "
                   "all of that protein's precursors share the value",
    "NPeptides": "number of precursors quantifying the protein (limpa)",
    "PropObs": "fraction of the protein's precursor intensities that were OBSERVED rather than "
               "missing, over all its precursors and samples (limpa) -- low = mostly inferred",
    "logFC": "log2 fold change for the contrast: first group minus second (contrast `A-B`; the "
             "file name writes it `A.B`). Positive = higher in the first group",
    "AveExpr": "average log2 expression of the protein across all samples (limma)",
    "t": "moderated t-statistic (limma eBayes)",
    "P.Value": "raw p-value of the moderated t-test",
    "adj.P.Val": "Benjamini-Hochberg adjusted p-value -- THE significance column",
    "B": "log-odds that the protein is differentially expressed (limma)",
}
QC_COLUMNS = {
    "Sample": "the run (conditions.csv File.Name)", "Group": "its group",
    "Detected": "proteins with at least one precursor actually observed in this run",
    "Inferred": "proteins whose value in this run the detection-probability model supplied",
    "Total": "proteins in the matrix", "PctDetected": "Detected / Total, %",
    "PctInferred": "Inferred / Total, %",
}
# Group names that are, by name, the negative control of a pull-down / IP (IgG, beads-only).
# Only ever used to POINT at contrasts against them; the user's design says what they are.
IP_CONTROL_NAME = re.compile(r"(?i)(^|[_\-. ])(igg|beads?)($|[_\-. ])")


# ---------------------------------------------------------------------------- reading
def _load(path):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except (OSError, ValueError, TypeError):
        return None


def _text(path):
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            return fh.read()
    except (OSError, TypeError):
        return ""


def _csv_header(path):
    try:
        with open(path, newline="", encoding="utf-8", errors="replace") as fh:
            return next(csv.reader(fh), [])
    except OSError:
        return []


def _sig_counts(prov):
    """significant_per_contrast as a dict. Older jsonlite-less de_provenance.json files carry it
    as the Python repr of a dict; anything unreadable is {} (no counts shown, none invented)."""
    v = (prov or {}).get("significant_per_contrast")
    if isinstance(v, dict):
        return v
    if isinstance(v, str):
        try:
            d = ast.literal_eval(v)
            return d if isinstance(d, dict) else {}
        except (ValueError, SyntaxError):
            return {}
    return {}


def raw_record(p):
    """(raw paths, where the list came from): input/raw_files.txt, else the search's own record
    of what it read -- search_provenance.json `files`, output/search/file_list.txt -- else the
    reproducibility bundle's run_manifest.json."""
    raws = read_raw_list(p["session_dir"])
    if raws:
        return raws, "input/raw_files.txt"
    sp = _load(p["search_prov"]) or {}
    if isinstance(sp.get("files"), list) and sp["files"]:
        return [str(x) for x in sp["files"]], "output/search/search_provenance.json"
    fl = os.path.join(p["search_out"], "file_list.txt")
    lines = [ln.strip() for ln in _text(fl).splitlines() if ln.strip() and not ln.startswith("#")]
    if lines:
        return lines, "output/search/file_list.txt"
    rm = _load(os.path.join(p["repro_dir"], "run_manifest.json")) or {}
    raw = (rm.get("inputs") or {}).get("raw")
    if isinstance(raw, list) and raw:
        return [str(x) for x in raw], "output/reproducibility/run_manifest.json"
    return [], None


def ensure_raw_list(p):
    """Write input/raw_files.txt from the search's record when the session has none (a
    hive_remote session is initialised without --raw). Returns (level, note) for MANIFEST.txt."""
    if read_raw_list(p["session_dir"]):
        return "OK", f"present ({len(read_raw_list(p['session_dir']))} files)"
    raws, src = raw_record(p)
    if not raws:
        return "SKIPPED", ("no record of the raw files: no search_provenance.json `files`, "
                           "output/search/file_list.txt or run_manifest.json")
    os.makedirs(os.path.dirname(p["raw_list"]), exist_ok=True)
    with open(p["raw_list"], "w", encoding="utf-8") as fh:
        fh.write(f"{RAW_HEADER}\n# Written by session.py from {src}: the paths are where the "
                 f"search read them.\n" + "".join(r.rstrip("/") + "\n" for r in raws))
    return "OK", f"written from {src} ({len(raws)} files)"


def _session_root_of(path, levels):
    """The session folder `levels` above a recorded file (…/input/conditions.csv -> 2)."""
    if not path or not os.path.isabs(str(path).replace("\\", "/")):
        return None
    d = str(path)
    for _ in range(levels):
        d = os.path.dirname(d)
    return d or None


def _audit_notes(p):
    """(overall, [notes]) from AUDIT.json/AUDIT.md: every WARN/FAIL finding, verbatim."""
    out = p["output_dir"]
    j = _load(os.path.join(out, "AUDIT.json")) or _load(os.path.join(p["session_dir"],
                                                                     "AUDIT.json"))
    if isinstance(j, dict):
        notes = [f"{f.get('status')}: {f.get('check')} -- {f.get('message')}"
                 for f in (j.get("findings") or []) if isinstance(f, dict)
                 and f.get("status") in ("WARN", "FAIL")]
        return j.get("overall"), notes
    md = _text(os.path.join(out, "AUDIT.md"))
    if not md:
        return None, []
    m = re.search(r"\*\*Overall:\s*([^*]+)\*\*", md)
    notes = [re.sub(r"\*\*", "", ln.lstrip("- ").strip()) for ln in md.splitlines()
             if ln.lstrip().startswith("- ") and ("⚠" in ln or "❌" in ln)]
    return (m.group(1).strip() if m else None), notes


def _quality_notes(p):
    """SAMPLE_QUALITY flags (json) or its warning sections (md): heading + first paragraph."""
    out = p["output_dir"]
    j = _load(os.path.join(out, "SAMPLE_QUALITY.json"))
    if isinstance(j, dict) and j.get("flags"):
        return [str(x) for x in j["flags"]][:20]
    md = _text(os.path.join(out, "SAMPLE_QUALITY.md"))
    notes, lines = [], md.splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("#") and "⚠" in ln:
            para = []
            for nx in lines[i + 1:]:
                if nx.startswith("#") or nx.startswith("|"):
                    break
                if nx.strip():
                    para.append(nx.strip())
                elif para:
                    break
            text = " ".join(para)
            if text.endswith(":"):
                text += " (the table is in SAMPLE_QUALITY.md)"
            notes.append(ln.lstrip("# ").strip() + (": " + text if text else ""))
    return notes


def _manifest(p):
    lines = _text(p["manifest_txt"]).splitlines()
    return [ln for ln in lines if ln.startswith("[")]


def significance_text(de):
    """The significance rule as de_provenance.json records it. An older record without
    `significance_rule` is described only from the fields it has (adjp, logfc_role)."""
    rule, adjp, role = de.get("significance_rule"), de.get("adjp"), de.get("logfc_role")
    if rule:
        return rule + (f" (adjp = {adjp})" if adjp is not None else "")
    if adjp is not None and role == "reference_line_only":
        return f"adj.P.Val < {adjp} (BH); |log2FC| is a volcano reference line only"
    return NOT_RECORDED + (f" (adjp = {adjp}; whether a fold-change filter was applied is not "
                           "recorded)" if adjp is not None else "")


def design_terms(design):
    """['groups', 'Batch'] from '~ 0 + groups + Batch'."""
    return [t.strip() for t in str(design or "").lstrip("~").split("+")
            if t.strip() and t.strip() not in ("0", "1")]


def submission_line(p):
    """The CoreOmics submission line (submission_report.py attach), or None when the session
    has none. A record that is there but unreadable is said so, not left out."""
    if not (os.path.isfile(p["submission_record"]) or os.path.isfile(p["session_json"])):
        return None
    try:
        import submission_report
        rec = submission_report.load(p["session_dir"])
    except Exception as e:
        return f"CoreOmics submission: could not be read ({type(e).__name__}: {e})"
    return submission_report.one_line(rec) if rec else None


def gather(session_dir, registry=None, registry_note=None, pending=(), located_at=None):
    """Every fact the three documents state, read from the session. `registry`: the Core
    run-registry folder of this session, when known (record_run.locate() or its result).
    `pending`: keys of files the caller writes right after (finalize: "manifest"), linked as if
    present. `located_at`: where the session really is, when these documents are written in a
    copy of it (session.py docs --as)."""
    p = paths_for(session_dir)
    sd = p["session_dir"]
    de = _load(os.path.join(p["de_dir"], "de_provenance.json")) or {}
    sp = _load(p["search_prov"]) or {}
    wf = _load(p["workflow_manifest"]) or {}
    rm = _load(os.path.join(p["repro_dir"], "run_manifest.json")) or {}
    fm = _load(p["fasta_meta"]) or {}
    q = rm.get("query") or {}
    f = {"p": p, "name": os.path.basename(sd), "de": de, "sp": sp, "wf": wf, "fm": fm}
    rel = lambda path: os.path.relpath(path, sd).replace(os.sep, "/")
    has = lambda path: os.path.exists(path)

    # --- the study
    f["submission"] = submission_line(p)
    f["organism"] = fm.get("organism") or None
    f["taxid"] = fm.get("taxid") or wf.get("organism_taxid") or q.get("organism_taxid")
    f["instrument"] = q.get("instrument") or next(iter(wf.get("instruments") or []), None)
    f["acquisition"] = (wf.get("acquisition") or q.get("acquisition") or None)
    eng = sp.get("engine") or (wf.get("engine") or {}).get("name")
    try:
        from make_methods import ENGINE_LABEL          # the one engine-name table
    except Exception:
        ENGINE_LABEL = {}
    f["engine"] = ENGINE_LABEL.get(str(eng).lower(), eng) if eng else None
    f["engine_version"] = sp.get("version")
    f["engine_version_src"] = "search_provenance.json" if sp.get("version") else None
    if not f["engine_version"] and (wf.get("engine") or {}).get("version"):
        f["engine_version"] = wf["engine"]["version"]
        f["engine_version_src"] = "workflow manifest pin -- the search's own record is not here"
    raws, raw_src = raw_record(p)
    f["raws"], f["raw_src"] = raws, raw_src
    conds = []
    try:
        with open(p["conditions"], newline="", encoding="utf-8") as fh:
            conds = list(csv.DictReader(fh))
    except OSError:
        pass
    groups = de.get("groups") if isinstance(de.get("groups"), dict) else {}
    if not groups and conds and "Group" in conds[0]:
        for r in conds:
            groups[r["Group"]] = groups.get(r["Group"], 0) + 1
    f["groups"] = groups
    f["n_samples"] = de.get("n_samples") or (len(conds) if conds else None)
    f["contrasts"] = de.get("contrasts") or []
    f["sig"] = _sig_counts(de)
    f["controls"] = sorted(g for g in groups if IP_CONTROL_NAME.search(g))

    # --- files that exist (links are only ever to these)
    out = p["output_dir"]
    f["files"] = {k: v for k, v in {
        "report_html": os.path.join(out, "Analysis_Report.html"),
        "report_docx": os.path.join(out, "AI_Analysis_Report.docx"),
        "report_md": p["analysis_report"],
        "methods_docx": p["methods_docx"], "methods_md": p["methods_md"],
        "howto_html": os.path.join(p["deposit_dir"], "HOW_TO_SUBMIT.html"),
        "howto_md": os.path.join(p["deposit_dir"], "HOW_TO_SUBMIT.md"),
        "tables": p["de_dir"] if glob.glob(os.path.join(p["de_dir"], "*")) else None,
        "expr": os.path.join(p["de_dir"], "Expression_Matrix.csv"),
        "qc_di": os.path.join(p["de_dir"], "QC_detected_vs_inferred.csv"),
        "methods_txt": os.path.join(p["de_dir"], "methods.txt"),
        "de_prov": os.path.join(p["de_dir"], "de_provenance.json"),
        "repro_R": os.path.join(p["de_dir"], "reproducibility_log.R"),
        "reproduce_md": os.path.join(p["repro_dir"], "REPRODUCE.md"),
        "output_files": p["output_files_md"], "manifest": p["manifest_txt"],
        "conditions": p["conditions"], "raw_list": p["raw_list"], "fasta_meta": p["fasta_meta"],
        "search_prov": p["search_prov"], "report": os.path.join(p["search_out"], "report.parquet"),
        "params": next((x for x in (os.path.join(p["search_out"], "params.resolved.cfg"),
                                    os.path.join(p["workflow_dir"], "params.cfg"),
                                    os.path.join(p["workflow_dir"], "params.json"))
                        if has(x)), None),
        "figures_json": os.path.join(p["figures_dir"], "figures.json"),
        "audit": os.path.join(out, "AUDIT.md"),
        "quality": os.path.join(out, "SAMPLE_QUALITY.md"),
        "agents": os.path.join(sd, "AGENTS.md"),
    }.items() if v and (has(v) or k in pending)}
    f["de_files"] = sorted(glob.glob(os.path.join(p["de_dir"], "DE_*.csv")))
    f["rel"] = rel
    f["audit_overall"], f["audit_notes"] = _audit_notes(p)
    f["quality_notes"] = _quality_notes(p)
    f["manifest_lines"] = _manifest(p)
    fran = _load(os.path.join(p["search_out"], "fran_deposit.json")) or {}
    f["fran"] = fran

    # --- where it all lives on HIVE
    rows = share_map.load_table()
    loc = lambda path: share_map.locate(path, rows)
    L = [dict(loc(located_at or sd), what="This session folder")]
    ran = (_session_root_of(de.get("metadata"), 2) if str(de.get("metadata") or "").endswith(
        "conditions.csv") else None) or \
        (_session_root_of(sp.get("params_file"), 3) if "/input/wf/" in str(
            sp.get("params_file") or "") else None)
    if ran:
        r = loc(ran)
        if r["hive"] and r["hive"] != L[0]["hive"]:
            L.append(dict(r, what="The session where the analysis ran (per the DE and search "
                                  "records; this folder is a later copy or move of it)"))
    by_dir = {}
    for r in raws:
        by_dir.setdefault(os.path.dirname(r.rstrip("/").replace("\\", "/")), []).append(r)
    for i, (d, members) in enumerate(sorted(by_dir.items())):
        if i == 3:
            L.append({"what": f"Raw data: {len(by_dir) - 3} more folder(s)", "path": None,
                      "hive": None, "how": "see input/raw_files.txt", "windows": None,
                      "mac": None})
            break
        L.append(dict(loc(d), what=f"Raw data ({len(members)} of {len(raws)} files; full list: "
                                   f"input/raw_files.txt)"))
    if not raws:
        L.append({"what": "Raw data", "path": None, "hive": None, "windows": None, "mac": None,
                  "how": f"{NOT_RECORDED}: no raw-file list in this session"})
    res = sp.get("result") if isinstance(sp.get("result"), dict) else {}
    s_out = res.get("out") or (os.path.dirname(res["report"]) if res.get("report") else None)
    if not s_out and sp:                        # the search ran in this session's output/search
        s_out = p["search_out"]
    L.append(dict(loc(s_out), what="Search output (report.parquet, logs)") if s_out else
             {"what": "Search output", "path": None, "hive": None, "windows": None, "mac": None,
              "how": f"{NOT_RECORDED}: search_provenance.json names no output folder"})
    fasta = sp.get("fasta") or fm.get("fasta")
    L.append(dict(loc(fasta), what="The FASTA the search read") if fasta else
             {"what": "The FASTA the search read", "path": None, "hive": None, "windows": None,
              "mac": None, "how": f"{NOT_RECORDED}"})
    L.append({"what": "Core run-registry record (Proteomics Core members)", "path": registry,
              "hive": registry if registry and registry.startswith("/") else None,
              "how": ("record_run.py" if registry else
                      f"{NOT_RECORDED}{': ' + registry_note if registry_note else ''}"),
              **share_map.views_of(registry, rows)})
    if fran.get("entry"):
        status = fran.get("status") or "status not recorded"
        L.append(dict(loc(fran["entry"]), what=f"FRAN hand-off entry ({status})"))
    f["locations"] = L
    return f


# ---------------------------------------------------------------------------- rendering
def _link(f, key, text):
    path = f["files"].get(key)
    if not path:
        return None
    rel = f["rel"](path) + ("/" if os.path.isdir(path) else "")
    return f"[{text}]({urllib.parse.quote(rel, safe='/')})"


def _esc(s):
    return str(s).replace("|", "\\|").replace("\n", " ")


def locations_table(f):
    rows = []
    for r in f["locations"]:
        hive = f"`{r['hive']}`" if r.get("hive") else NOT_RECORDED
        win = f"`{r['windows']}`" if r.get("windows") else ("—" if not r.get("hive") else
                                                             NOT_RECORDED)
        mac = f"`{r['mac']}`" if r.get("mac") else ("—" if not r.get("hive") else NOT_RECORDED)
        how = r.get("how") or ""
        rows.append(f"| {_esc(r['what'])} | {_esc(hive)} | {_esc(win)} | {_esc(mac)} | "
                    f"{_esc(how)} |")
    smb = [f"`smb://{t['server']}/{t['share']}`" for t in share_map.load_table()
           if t["server"] and t["server"] != "*" and t["mac"]]
    return "\n".join(
        ["| What | On HIVE | Windows | Mac | How it was found |", "|---|---|---|---|---|", *rows,
         "", "> [!NOTE]", "> To browse there: on **Windows**, paste the Windows path into File Explorer's "
             "address bar (a PC that maps the share to a drive letter, e.g. `R:`, has the same "
             "folders below that letter). On a **Mac**, Finder → Go → Connect to Server"
         + (f" ({', '.join(smb)})" if smb else "") + ", then open the Mac path. On **HIVE**, "
         "`cd` to the HIVE path. \"not recorded\" means the session's records do not say — it is "
         "never filled with a guess. Raw data is not copied into this folder."])


def summary_lines(f):
    de = f["de"]
    org = (f"{f['organism']} (taxid {f['taxid']})" if f["organism"] else
           f"taxid {f['taxid']}" if f["taxid"] else NOT_RECORDED)
    eng = f"{f['engine']} {f['engine_version'] or ''}".strip() if f["engine"] else NOT_RECORDED
    L = ([f"- {f['submission']}"] if f.get("submission") else []) + [
         f"- Organism: {org}",
         f"- Instrument / acquisition: {f['instrument'] or NOT_RECORDED} / "
         f"{f['acquisition'] or NOT_RECORDED}",
         f"- Search engine: {eng}" + (f" ({f['engine_version_src']})"
                                      if f["engine_version_src"] and "pin" in
                                      f["engine_version_src"] else ""),
         f"- Raw files: {len(f['raws']) or NOT_RECORDED}"]
    if de:
        L.append(f"- Differential expression: {de.get('display_label') or NOT_RECORDED}")
        L.append(f"- Significance: {significance_text(de)}")
    else:
        L.append("- Differential expression: not run in this session (no "
                 "output/tables/de_provenance.json)")
    if f["groups"]:
        L.append("- Groups: " + ", ".join(f"{g} ({n})" for g, n in sorted(f["groups"].items())))
    return L


def readme_md(f, for_html=False):
    """THE README source. README.md is this text; README.html is this text rendered (for_html
    drops only the line pointing the Markdown reader at the .html)."""
    de, files = f["de"], f["files"]
    L = [f"# {f['name']}", "",
         "Proteomics results from the UC Davis Proteomics Core pipeline "
         "(the `ucdavis-proteomics-core-pipeline` skill)."]
    if not for_html:
        L += ["", "**Open `README.html`** (double-click it) — the same page, with working links."]
    start = []
    for key, text, what in (
            ("report_html", "Analysis report", " — figures, QC and the interpretation in one page"),
            ("methods_docx", "Methods for the paper (Word)", ""),
            ("methods_md", "Methods for the paper (text)", ""),
            ("tables", "Results tables", " — one `DE_*.csv` per comparison, "
                                         "`Expression_Matrix.csv`, `methods.txt`"),
            ("howto_html", "How to deposit the data in PRIDE / MassIVE", ""),
            ("agents", "AGENTS.md", " — for an AI assistant: give it this file with the folder"),
            ("manifest", "MANIFEST.txt", " — what this export contains, and anything that could "
                                         "not be made (with the reason)")):
        link = _link(f, key, text)
        if link:
            start.append(f"- {link}{what}")
    if not files.get("agents"):
        start.append("- `AGENTS.md` — for an AI assistant: give it this file with the folder")
    L += ["", "## Start here", "", *start]
    L += ["", "## Summary", "", *summary_lines(f)]
    if f["contrasts"]:
        L += ["", "| Comparison | Significant proteins |", "|---|---|"]
        for c in f["contrasts"]:
            n = f["sig"].get(c)
            L.append(f"| {_esc(c)} | {n if n is not None else NOT_RECORDED} |")
    L += ["", "## Where this lives on HIVE", "", locations_table(f)]
    L += ["", "## Where everything is (in this folder)", "", "| Path | What it is |", "|---|---|"]
    for path, what in (
            ("input", "conditions.csv (the design), the FASTA sidecar, the search parameters"
                      + (", raw_files.txt (where the raw data are)" if files.get("raw_list")
                         else "")),
            ("output/search", "the search output (report, engine logs, search_provenance.json)"),
            ("output/tables", "DE results (DE_*.csv), Expression_Matrix.csv, methods.txt, "
                              "de_provenance.json, reproducibility_log.R (the analysis as R)"),
            ("output/figures", "plots; captions in figures.json"),
            ("output/reproducibility", "the pinned bundle for re-running the search too"),
            ("output/DATA_SUBMISSION", "everything to deposit the data in PRIDE / MassIVE"),
            ("output/OUTPUT_FILES.md", "a catalog of every file"),
            ("scripts", "the skill scripts that ran"), ("logs", "commands.log + engine logs")):
        if os.path.exists(os.path.join(f["p"]["session_dir"], path)):
            L.append(f"| `{path}` | {what} |")
    L += ["", "## Reproduce", ""]
    L.append("**The analysis, as code:** `output/tables/reproducibility_log.R` — the whole "
             "differential-expression analysis in plain R with every value written out."
             if files.get("repro_R") else
             "**The analysis, as code:** not in this folder (no reproducibility_log.R).")
    L += ["", ("**The whole run, pinned** (search included, takes hours, needs the raw data on "
               "HIVE): `output/reproducibility/REPRODUCE.md`." if files.get("reproduce_md") else
               "**The whole run, pinned:** not in this folder (no REPRODUCE.md).")]
    L += ["", "## Methods", ""]
    L.append((f"**For the paper:** `{f['rel'](files['methods_md'])}` (Word: the `.docx` beside it)"
              " — LC-MS acquisition, the search, the database, the DE and the UC Davis "
              "instrument-grant acknowledgment. Resolve every `[... — confirm]` tag before "
              "publishing.") if files.get("methods_md") else
             "**For the paper:** no publication Methods in this folder — `MANIFEST.txt` says why.")
    if files.get("methods_txt"):
        L += ["", "The DE step's own record is `output/tables/methods.txt`."]
    L += ["", "## Deposit the data (PRIDE / MassIVE)", ""]
    L.append("Everything to deposit the data is in `output/DATA_SUBMISSION/` — start with "
             "**`HOW_TO_SUBMIT.html`** (the same as `HOW_TO_SUBMIT.md`)." if files.get("howto_md")
             else "The deposit package was not written — `MANIFEST.txt` says why."
             if files.get("manifest") else
             "No deposit package in this folder (the session was finalized before the skill "
             "wrote one; `session.py finalize` writes it).")
    L += ["", "## What is in this export", "",
          ("`MANIFEST.txt` lists every part of the Methods and deposit package as [OK], or "
           "[SKIPPED] with the reason." if files.get("manifest") else
           "This session has no `MANIFEST.txt` (it was finalized before the skill wrote one).")
          + " `AGENTS.md` is a guide to this folder for an AI agent.", "",
          f"*Written {datetime.date.today().isoformat()} by `session.py` from this session's "
          "records.*", ""]
    return "\n".join(L)


def agents_md(f):
    """AGENTS.md: what an AI agent handed this folder needs, from the records only."""
    de, files, rel = f["de"], f["files"], f["rel"]
    L = ["# AGENTS.md — guide to this folder for an AI agent", "",
         "Generated by `session.py` from this session's own records. Every fact here is read "
         "from a file in this folder; where a file below is named authoritative, it wins over "
         "this summary and over `AI_Analysis_Report.md` (an interpretation, not a record).", ""]
    L += ["## The study", "", *summary_lines(f)]
    if f["contrasts"]:
        adjp = de.get("adjp")
        L.append("- Contrasts (significant = adj.P.Val < " + (str(adjp) if adjp is not None
                                                              else NOT_RECORDED) + "): "
                 + "; ".join(f"{c} = {f['sig'].get(c, NOT_RECORDED)}" for c in f["contrasts"]))
    if de.get("design"):
        terms = design_terms(de["design"])
        L.append(f"- Design: `{de['design']}`" + (f" — covariates: {', '.join(terms[1:])}"
                                                  if len(terms) > 1 else ""))
    if de.get("input"):
        L.append(f"- The DE read `{os.path.basename(str(de['input']))}` "
                 "(de_provenance.json `input`).")

    L += ["", "## Which file answers what (authoritative first)", "",
          "| Question | File |", "|---|---|"]
    for q, key in (("What the DE did: pipeline, filters, design, thresholds, versions",
                    "de_prov"),
                   ("The Methods text of the DE (verbatim, do not paraphrase)", "methods_txt"),
                   ("DE results, one table per contrast", None),
                   ("Protein abundance per sample (log2)", "expr"),
                   ("Which values were measured vs inferred, per sample", "qc_di"),
                   ("Samples and their groups", "conditions"),
                   ("Search engine, the version that ran, exact command", "search_prov"),
                   ("Search parameters as run", "params"),
                   ("Precursor-level search results", "report"),
                   ("The database: organism, UniProt release, contaminants", "fasta_meta"),
                   ("Where the raw data are", "raw_list"),
                   ("Figure captions", "figures_json"),
                   ("Pitfall audit", "audit"), ("Sample quality / contamination", "quality"),
                   ("Publication Methods (LC-MS + search + DE + acknowledgment)", "methods_md"),
                   ("The report of record (QC, figures, interpretation)", "report_html"),
                   ("What this export contains, and what is missing", "manifest"),
                   ("The DE as runnable R", "repro_R"),
                   ("Re-running everything", "reproduce_md")):
        if key is None:
            n = len(f["de_files"])
            v = (f"`output/tables/DE_<pipeline>_<contrast>.csv` ({n} file{'s' if n != 1 else ''})"
                 if n else "not in this folder")
        else:
            v = f"`{rel(files[key])}`" if files.get(key) else "not in this folder"
        L.append(f"| {q} | {v} |")

    L += ["", "## Table columns", ""]
    if f["de_files"]:
        hdr = _csv_header(f["de_files"][0])
        L.append(f"`DE_*.csv` (header of `{os.path.basename(f['de_files'][0])}`):")
        L += [f"- `{c}` — {COLUMNS.get(c, 'no description recorded')}" for c in hdr]
    if files.get("expr"):
        hdr = _csv_header(files["expr"])
        ann = [c for c in hdr if c in COLUMNS]
        ann_s = ", ".join(f"`{c}`" for c in ann) or "no"
        L += ["", f"`Expression_Matrix.csv`: {ann_s} annotation columns, then "
                  f"{len(hdr) - len(ann)} sample columns (named as "
                  "`File.Name` in conditions.csv). Values: log2 protein quantities from "
                  f"{de.get('rollup_method') or 'the quantification (not recorded)'}."]
    if files.get("qc_di"):
        hdr = _csv_header(files["qc_di"])
        L += ["", "`QC_detected_vs_inferred.csv`: "
              + "; ".join(f"`{c}` = {QC_COLUMNS.get(c, 'no description recorded')}" for c in hdr)]

    L += ["", "## Traps — read before you compute anything", ""]
    role = de.get("logfc_role")
    if de:
        L.append(f"- **Significance:** {significance_text(de)}."
                 + (f" |log2FC| = {de.get('logfc')} is a reference line on the volcano only: do "
                    "not add a fold-change cutoff when you count, rank or describe significant "
                    "proteins, and never call a protein significant because it is ≥2-fold."
                    if role == "reference_line_only" else
                    " Whether a fold-change filter was applied is not recorded — do not assume "
                    "one."))
        if de.get("missing_policy"):
            L.append(f"- **Missing values:** {de['missing_policy']}")
    if files.get("qc_di"):
        pct = []
        try:
            with open(files["qc_di"], newline="") as fh:
                pct = [float(r["PctInferred"]) for r in csv.DictReader(fh)
                       if r.get("PctInferred") not in (None, "", "NA")]
        except (OSError, ValueError, KeyError):
            pct = []
        span = f" ({min(pct):g}–{max(pct):g}% of values per sample)" if pct else ""
        L.append(f"- **Inferred is not measured.** `Expression_Matrix.csv` has a value in "
                 f"cells where no precursor was observed{span}; `QC_detected_vs_inferred.csv` "
                 "counts them per sample and `PropObs` per protein. Never report an inferred "
                 "value as a detection, and never use non-empty cell counts as depth — a "
                 "\"0% missing\" in the audit reflects this filled matrix, not detection.")
    fm = f["fm"]
    if fm.get("n_contaminants_appended") or fm.get("n_contaminants_already_present"):
        n = (fm.get("n_contaminants_appended") or 0) + (fm.get("n_contaminants_already_present")
                                                        or 0)
        tag = fm.get("diann_cont_quant_exclude")
        L.append(f"- **Contaminants:** the database holds {n} contaminant sequences"
                 + (f" ({fm.get('contaminant_set')} set)" if fm.get("contaminant_set") else "")
                 + (f"; their protein IDs start with `{tag}` and the search kept them out of "
                    "quantification and normalisation" if tag else "")
                 + ". A protein you expect but cannot find may be listed under a contaminant "
                   "entry.")
    if f["controls"] and f["contrasts"]:
        vs = [c for c in f["contrasts"] if any(re.search(rf"(^|-){re.escape(g)}$", c)
                                               for g in f["controls"])]
        if vs:
            L.append(f"- **Pull-down controls, by group name: {', '.join(f['controls'])}.** "
                     f"{len(vs)} contrast(s) compare against them ({', '.join(vs[:4])}"
                     f"{', …' if len(vs) > 4 else ''}): a positive logFC there is enrichment "
                     "over that control, not a change in abundance. A protein never observed in "
                     "the control gets its control values from the missing-value policy above, "
                     "so its logFC rests on that policy, not on a measured baseline. (Check "
                     "conditions.csv: the name is the only evidence these are controls.)")
    skipped = [ln for ln in f["manifest_lines"] if ln.startswith("[SKIPPED]")]
    if skipped:
        L.append("- **Missing parts** (MANIFEST.txt): a [SKIPPED] part is absent, not empty — "
                 + "; ".join(re.sub(r"\s+", " ", s[len("[SKIPPED]"):]).strip() for s in skipped))

    L += ["", "## Data quality notes to respect", ""]
    if f["audit_overall"] or f["audit_notes"]:
        L.append(f"- Audit (`output/AUDIT.md`): {f['audit_overall'] or NOT_RECORDED}")
        L += [f"  - {n}" for n in f["audit_notes"]]
    if f["quality_notes"]:
        L.append("- Sample quality (`output/SAMPLE_QUALITY.md`):")
        L += [f"  - {n}" for n in f["quality_notes"]]
    if not (f["audit_overall"] or f["audit_notes"] or f["quality_notes"]):
        L.append("- No AUDIT or SAMPLE_QUALITY record in this folder.")

    L += ["", "## Reproduce", ""]
    L.append(f"- The DE: `{rel(files['repro_R'])}` (plain R, every value literal; it reads "
             f"`{os.path.basename(str(de.get('input') or 'the search report'))}`)."
             if files.get("repro_R") else "- The DE as R: not in this folder.")
    L.append(f"- Search + DE: `{rel(files['reproduce_md'])}` — needs the raw data, which are on "
             "HIVE (below), not here." if files.get("reproduce_md") else
             "- The full run: no REPRODUCE.md in this folder.")

    L += ["", "## Where this lives on HIVE", "", locations_table(f)]

    L += ["", "## Do not", "",
          "- Invent a value, protein, count or threshold that is not in these files.",
          "- Treat an inferred value as a measured one, or a missing/[SKIPPED] part as empty.",
          "- Claim a threshold that was not applied (see Significance above) or describe the "
          "pipeline differently from `de_provenance.json` / `methods.txt`.",
          "- Re-run the search from this folder: the raw data are not in it (see the HIVE "
          "table), and a re-run needs HIVE and the pinned engine (`REPRODUCE.md`).",
          "- Edit the records (`*_provenance.json`, `methods.txt`, `MANIFEST.txt`); write new "
          "files instead.", ""]
    return "\n".join(L)


def _render_html(md, title):
    import make_deposit                               # the one Markdown -> HTML renderer
    return make_deposit.md_to_html(md, title, semantic=True)


def write_docs(session_dir, man=None, registry=None, registry_note=None, pending=(),
               located_at=None):
    """README.md, README.html and AGENTS.md at the session root. Each is its own [OK]/[SKIPPED]
    line when `man` (make_deposit.Manifest) is given -- one that fails never stops the others.
    Returns {"readme", "readme_html", "agents"} with the paths written (None if not)."""
    written = {"readme": None, "readme_html": None, "agents": None}

    def section(name, fn):
        if man is not None:
            return man.section(name, fn)
        fn()
        return True

    f = gather(session_dir, registry, registry_note, pending, located_at)
    sd = f["p"]["session_dir"]

    def _md():
        with open(f["p"]["readme"], "w", encoding="utf-8") as fh:
            fh.write(readme_md(f))
        written["readme"] = f["p"]["readme"]

    def _html():
        path = os.path.join(sd, "README.html")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(_render_html(readme_md(f, for_html=True), f["name"]))
        written["readme_html"] = path

    def _agents():
        path = os.path.join(sd, "AGENTS.md")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write(agents_md(f))
        written["agents"] = path

    # AGENTS.md first, so the README's link to it is to a file that exists
    section("AGENTS.md (guide for an AI agent)", _agents)
    if written["agents"]:
        f["files"]["agents"] = written["agents"]
    section("README.md", _md)
    section("README.html (for collaborators: open this)", _html)
    return written
