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
from session import (paths_for, read_raw_list, raw_list_encoding_note,  # noqa: E402  one place
                     params_file, NOT_RECORDED)
import share_map                                    # noqa: E402
import make_podcast                                 # noqa: E402  the optional audio discussion
from fetch_fasta import CONT_TAG                    # noqa: E402  the contaminant tag, one place

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
    have = read_raw_list(p["session_dir"])
    if have:
        note = raw_list_encoding_note(p["session_dir"])
        return "OK", f"present ({len(have)} files)" + (f"; {note}" if note else "")
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


# How the rule run_de.R records ("adj.P.Val < adjp (BH); no fold-change filter") reads for people.
# topTable() adjusts each contrast's p-values on its own (run_de.R, adjust.method = "BH"), hence
# "within each comparison" -- the wording the report uses too.
_RULE_WORDS = (("adj.P.Val < adjp", "adjusted p < {adjp}"),
               ("(BH)", "(Benjamini–Hochberg, within each comparison)"))


def significance_text(de):
    """The significance rule as de_provenance.json records it, written for people: its
    `significance_rule` with the recorded adjp in place of the variable ("adjusted p < 0.05
    (Benjamini–Hochberg, within each comparison); no fold-change filter"). A threshold the record
    lacks is tagged, never assumed. An older record without `significance_rule` is described only
    from the fields it has (adjp, logfc_role)."""
    rule, adjp, role = de.get("significance_rule"), de.get("adjp"), de.get("logfc_role")
    thr = (f"{adjp:g}" if isinstance(adjp, (int, float)) and not isinstance(adjp, bool) else
           str(adjp) if adjp is not None else f"[threshold {NOT_RECORDED}]")
    if rule:
        text = str(rule)
        for recorded, said in _RULE_WORDS:
            text = text.replace(recorded, said.format(adjp=thr))
        # a rule worded some other way keeps its own words, and the value it names
        return text + (f" (adjp = {thr})" if "adjp" in text else "")
    if adjp is not None and role == "reference_line_only":
        return (f"adjusted p < {thr} (Benjamini–Hochberg, within each comparison); |log2FC| is a "
                "volcano reference line only")
    return NOT_RECORDED + (f" (adjp = {thr}; whether a fold-change filter was applied is not "
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


def feedback(p):
    """For a UC Davis Core run (submission_report.core_run: the test that shows the report's
    Submission section) -> {"prot"}: the README and AGENTS.md then carry the Core's feedback
    survey (core_submission.feedback_line / feedback_url). None otherwise -- and None when that
    cannot be told, since a user outside the Core must never be sent the Core's survey."""
    try:
        import submission_report
        return submission_report.core_run(p["session_dir"])
    except Exception:
        return None


def _submission_facts(p):
    """Who prepared the samples and the record's data-quality notes, from the attached CoreOmics
    record -- submission_report's readings (prepared_by, quality_notes), never re-derived here.
    None when there is no record; {"error"} when it cannot be read."""
    if not os.path.isfile(p["submission_record"]):
        return None
    try:
        import submission_report as sr
        rec = sr.load(p["session_dir"])
        if not rec:
            return None
        who, why = sr.prepared_by(rec)
        return {"who": who, "why": why, "peptides": sr.sent_as_peptides(rec),
                "notes": [n["text"] for n in sr.quality_notes(rec, p["session_dir"])]}
    except Exception as e:                      # said in AGENTS.md, not dropped
        return {"error": f"{type(e).__name__}: {e}"}


def _de_contaminants(de):
    """The DE's contaminant step in make_methods' words (de_contaminant_sentence reads run_de.R's
    `contaminants` record -- the one description of it), or None when make_methods is absent."""
    try:
        from make_methods import de_contaminant_sentence
    except Exception:
        return None
    return de_contaminant_sentence(de)


def _de_block(de):
    """The DE's blocking factor in make_methods' words (de_block_sentence reads run_de.R's
    `block` record), or None when the run fitted samples as independent."""
    try:
        from make_methods import de_block_sentence
    except Exception:
        return None
    return de_block_sentence(de)


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
    # a 0-byte file (a write that died) is not linked as if it held the thing it is named for
    has = lambda path: os.path.isdir(path) or (os.path.isfile(path) and
                                                os.path.getsize(path) > 0)

    # --- the study
    f["submission"] = submission_line(p)
    f["feedback"] = feedback(p)
    f["submission_facts"] = _submission_facts(p)
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
        # collect_conditions.py writes it in the computer's own encoding: a Windows sample name
        # must not cost the README (the groups shown are de_provenance.json's when it has them)
        with open(p["conditions"], newline="", encoding="utf-8", errors="replace") as fh:
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
    cont = de.get("contaminants") if isinstance(de.get("contaminants"), dict) else {}
    f["files"] = {k: v for k, v in {
        "report_html": os.path.join(out, "Analysis_Report.html"),
        "report_pdf": os.path.join(out, "Analysis_Report.pdf"),
        "report_twin": os.path.join(out, "Analysis_Report.md"),
        "report_docx": os.path.join(out, "AI_Analysis_Report.docx"),
        "report_md": p["analysis_report"],
        "methods_docx": p["methods_docx"], "methods_md": p["methods_md"],
        "howto_html": os.path.join(p["deposit_dir"], "HOW_TO_SUBMIT.html"),
        "howto_md": os.path.join(p["deposit_dir"], "HOW_TO_SUBMIT.md"),
        "tables": p["de_dir"] if glob.glob(os.path.join(p["de_dir"], "*")) else None,
        "expr": os.path.join(p["de_dir"], "Expression_Matrix.csv"),
        "qc_di": os.path.join(p["de_dir"], "QC_detected_vs_inferred.csv"),
        "det_matrix": os.path.join(p["de_dir"], (de.get("detection_matrix") or {}).get(
            "file") or "Detection_Matrix.csv"),
        "cont_removed": os.path.join(p["de_dir"], cont["removed_table"])
        if cont.get("removed_table") else None,
        "cont_share": os.path.join(p["de_dir"], cont["share_table"])
        if cont.get("share_table") else None,
        "submission": p["submission_record"],
        "methods_txt": os.path.join(p["de_dir"], "methods.txt"),
        "sets_prov": os.path.join(p["de_dir"], "sets_provenance.json"),
        "de_prov": os.path.join(p["de_dir"], "de_provenance.json"),
        "repro_R": os.path.join(p["de_dir"], "reproducibility_log.R"),
        "reproduce_md": os.path.join(p["repro_dir"], "REPRODUCE.md"),
        "output_files": p["output_files_md"], "manifest": p["manifest_txt"],
        "conditions": p["conditions"], "raw_list": p["raw_list"], "fasta_meta": p["fasta_meta"],
        "search_prov": p["search_prov"], "report": os.path.join(p["search_out"], "report.parquet"),
        "params": params_file(p),
        "figures_json": os.path.join(p["figures_dir"], "figures.json"),
        "audit": os.path.join(out, "AUDIT.md"),
        "quality": os.path.join(out, "SAMPLE_QUALITY.md"),
        "agents": os.path.join(sd, "AGENTS.md"),
        "podcast": os.path.join(out, "podcast", "podcast.json"),
        "conversation": os.path.join(sd, "logs", "conversation", "conversation.md"),
        "decisions": os.path.join(sd, "logs", "decisions.md"),
        "commands_log": p["commands_log"],
    }.items() if v and (has(v) or k in pending)}
    f["de_files"] = sorted(glob.glob(os.path.join(p["de_dir"], "DE_*.csv")))
    f["sets_files"] = sorted(glob.glob(os.path.join(p["de_dir"], "Sets_*.csv")))
    f["sets"] = _load(os.path.join(p["de_dir"], "sets_provenance.json")) or {}
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
        acc = r.get("access")               # hive_shares.tsv: who can open that share
        hive = f"`{r['hive']}`" if r.get("hive") else NOT_RECORDED
        win = (f"`{r['windows']}`" if r.get("windows") else acc if acc else
               "—" if not r.get("hive") else NOT_RECORDED)
        mac = (f"`{r['mac']}`" + (" (Core members)" if acc else "") if r.get("mac") else
               "—" if not r.get("hive") else NOT_RECORDED)
        how = r.get("how") or ""
        rows.append(f"| {_esc(r['what'])} | {_esc(hive)} | {_esc(win)} | {_esc(mac)} | "
                    f"{_esc(how)} |")
    table = share_map.load_table()
    # Connect to Server only for a share with a known server, and said for its own Mac paths:
    # smb://128.120.208.24/proteomics is Flinders, and does not reach /Volumes/proteomics-grp.
    smb = [f"`smb://{t['server']}/{t['share']}` for paths under `{t['mac']}`" for t in table
           if t["server"] and t["server"] != "*" and t["mac"]]
    used = {r.get("access") for r in f["locations"] if r.get("access")}
    limited = [f"Paths under `{t['hive']}`" + (f" (`{t['mac']}`)" if t["mac"] else "")
               + f": {t['access']}." for t in table if t["access"] in used]
    return "\n".join(
        ["| What | On HIVE | Windows | Mac | How it was found |", "|---|---|---|---|---|", *rows,
         "", "> [!NOTE]", "> To browse there: on **Windows**, paste the Windows path into File Explorer's "
             "address bar (a PC that maps the share to a drive letter, e.g. `R:`, has the same "
             "folders below that letter). On a **Mac**, Finder → Go → Connect to Server"
         + (f" ({'; '.join(smb)})" if smb else "") + ", then open the Mac path. On **HIVE**, "
         "`cd` to the HIVE path. " + "".join(x + " " for x in limited)
         + "\"not recorded\" means the session's records do not say — it is never filled with a "
           "guess. Raw data is not copied into this folder."])


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


def _named_tables(f):
    """(files key, what it is) for the tables and records the README names one by one -- each
    in its own record's words where it has them (de_provenance.json)."""
    dm = f["de"].get("detection_matrix") if isinstance(f["de"].get("detection_matrix"),
                                                         dict) else {}
    both = lambda k: " (and `.json` beside it)" if os.path.isfile(
        os.path.splitext(f["files"].get(k) or "")[0] + ".json") else ""
    return (("det_matrix", "which values were measured and which inferred, per protein and "
                           "sample" + (f": {dm['values']}" if dm.get("values") else "")),
            ("qc_di", "per sample: proteins detected (at least one precursor observed) vs "
                      "inferred by the model"),
            ("cont_removed", "the contaminant protein groups removed before the DE"),
            ("sets_prov", "the protein-set tests' record: sets, versions, methods, how to read "
                          "the `Sets_*.csv` tables"),
            ("audit", "the pitfall audit: PASS / WARN / FAIL per check" + both("audit")),
            ("quality", "sample quality and contamination flags" + both("quality")))


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
            ("report_pdf", "Analysis report (PDF)", " — the full report with figures, for NotebookLM "
                                                    "or printing (the HTML through its print "
                                                    "stylesheet)"),
            ("report_twin", "Analysis report (plain text)",
             " — the full report as plain text, for NotebookLM or other AI notebooks: the same "
             "sections as the HTML, with each figure's caption and the numbers it shows written "
             "out, and the top proteins per contrast"),
            ("methods_docx", "Methods for the paper (Word)", ""),
            ("methods_md", "Methods for the paper (text)", ""),
            ("tables", "Results tables", " — one `DE_*.csv` per comparison, "
                                         "`Expression_Matrix.csv`, `methods.txt`"
                                         + (", protein-set tests `Sets_*.csv`"
                                            if f.get("sets_files") else "")),
            ("howto_html", "How to deposit the data in PRIDE / MassIVE", ""),
            ("agents", "AGENTS.md", " — for an AI assistant: give it this file with the folder"),
            ("manifest", "MANIFEST.txt", " — what this export contains, and anything that could "
                                         "not be made (with the reason)")):
        link = _link(f, key, text)
        if link:
            start.append(f"- {link}{what}")
    if files.get("report_html") and not files.get("report_pdf"):
        # why is MANIFEST.txt's "Report PDF" line (finalize, html_to_pdf.py) -- not assumed here
        start.append("- `output/Analysis_Report.pdf` — not made"
                     + (" (the \"Report PDF\" line of `MANIFEST.txt` says why)"
                        if files.get("manifest") else "")
                     + ": open `Analysis_Report.html`, Print, Save as PDF")
    if not files.get("agents"):
        start.append("- `AGENTS.md` — for an AI assistant: give it this file with the folder")
    items = ([make_podcast.readme_item_md(f["p"]["output_dir"], f["p"]["session_dir"]),
              make_podcast.share_item_md(f["p"]["output_dir"], f["p"]["session_dir"])]
             if files.get("podcast") else [])
    at = 1 if files.get("report_html") else 0    # after the report it discusses
    start[at:at] = [x for x in items if x]
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
        if path == "output/tables":                  # the files a reader needs by name
            L += [f"| `{f['rel'](files[k])}` | {_esc(what)} |" for k, what in _named_tables(f)
                  if files.get(k)]
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
          + " `AGENTS.md` is a guide to this folder for an AI agent.", ""]
    if f.get("feedback"):                        # a Core run: the Core's survey
        import core_submission
        L += [core_submission.feedback_line(
            "readme", f["feedback"]["prot"], make_podcast.has_podcast(f["p"]["output_dir"])), ""]
    L += [f"*Written {datetime.date.today().isoformat()} by `session.py` from this session's "
          "records.*", ""]
    return "\n".join(L)


def agents_md(f, for_delivery=False):
    """AGENTS.md: what an AI agent handed this folder needs, from the records only.
    `for_delivery`: the collaborator's copy (core_submission deliver) -- without "Reviewing this
    analysis", whose records (the conversation, the decisions log) are Core-internal and stay in
    the session."""
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
    blk = _de_block(de) if de else None
    if blk:
        L.append(f"- Blocking: {blk}")
        L += [f"  - **Caveat:** {w}" for w in (de["block"].get("warnings") or [])]
    if de.get("input"):
        L.append(f"- The DE read `{os.path.basename(str(de['input']))}` "
                 "(de_provenance.json `input`).")

    L += ["", "## Which file answers what (authoritative first)", "",
          "| Question | File |", "|---|---|"]
    for q, key in (("What the DE did: pipeline, filters, design, thresholds, versions",
                    "de_prov"),
                   ("The Methods text of the DE (verbatim, do not paraphrase)", "methods_txt"),
                   ("DE results, one table per contrast", None),
                   ("Protein-set tests (camera + fry) per contrast, and how they were run",
                    "sets"),
                   ("Protein abundance per sample (log2)", "expr"),
                   ("Which values were measured vs inferred, per protein and sample",
                    "det_matrix"),
                   ("How many values were measured vs inferred, per sample", "qc_di"),
                   ("The contaminant protein groups removed before the DE", "cont_removed"),
                   ("Each run's contaminant share of the signal (QC)", "cont_share"),
                   ("Samples and their groups", "conditions"),
                   ("The CoreOmics submission: sample sheet, who prepared the samples, the "
                    "description as written (contacts and billing left out)", "submission"),
                   ("Search engine, the version that ran, exact command", "search_prov"),
                   ("Search parameters as run", "params"),
                   ("Precursor-level search results", "report"),
                   ("The database: organism, UniProt release, contaminants", "fasta_meta"),
                   ("Where the raw data are", "raw_list"),
                   ("Figure captions", "figures_json"),
                   ("Pitfall audit", "audit"), ("Sample quality / contamination", "quality"),
                   ("Publication Methods (LC-MS + search + DE + acknowledgment)", "methods_md"),
                   ("The report of record (QC, figures, interpretation)", "report_html"),
                   ("The same report as plain text: sections, figure captions and the numbers "
                    "each figure shows, top proteins per contrast", "report_twin"),
                   ("The same report as a PDF with its figures", "report_pdf"),
                   ("What this export contains, and what is missing", "manifest"),
                   ("The DE as runnable R", "repro_R"),
                   ("Re-running everything", "reproduce_md")):
        if key is None:
            n = len(f["de_files"])
            v = (f"`output/tables/DE_<pipeline>_<contrast>.csv` ({n} file{'s' if n != 1 else ''})"
                 if n else "not in this folder")
        elif key == "sets":
            n = len(f["sets_files"])
            v = (f"`output/tables/Sets_<contrast>.csv` ({n} file{'s' if n != 1 else ''}) + "
                 f"`{rel(files['sets_prov'])}`" if n and files.get("sets_prov") else
                 "not run for this analysis")
        else:
            v = f"`{rel(files[key])}`" if files.get(key) else "not in this folder"
        L.append(f"| {q} | {v} |")

    L += ["", "## Table columns", ""]
    if f["de_files"]:
        hdr = _csv_header(f["de_files"][0])
        L.append(f"`DE_*.csv` (header of `{os.path.basename(f['de_files'][0])}`):")
        # Detected_<group> / Evidence: described by run_de.R in de_provenance.json
        # (detection_matrix.de_columns) -- read, never restated here.
        dcols = ((de.get("detection_matrix") or {}).get("de_columns") or {}) if isinstance(
            de.get("detection_matrix"), dict) else {}

        def _desc(c):
            if c in COLUMNS:
                return COLUMNS[c]
            if c.startswith("Detected_") and dcols.get("Detected_<group>"):
                return f"group `{c[len('Detected_'):]}`: {dcols['Detected_<group>']}"
            return dcols.get(c, "no description recorded")
        L += [f"- `{c}` — {_desc(c)}" for c in hdr]
    if files.get("expr"):
        hdr = _csv_header(files["expr"])
        ann = [c for c in hdr if c in COLUMNS]
        ann_s = ", ".join(f"`{c}`" for c in ann) or "no"
        L += ["", f"`Expression_Matrix.csv`: {ann_s} annotation columns, then "
                  f"{len(hdr) - len(ann)} sample columns (named as "
                  "`File.Name` in conditions.csv). Values: log2 protein quantities from "
                  f"{de.get('rollup_method') or 'the quantification (not recorded)'}."]
    sets_cols = (f.get("sets") or {}).get("columns")
    if f.get("sets_files") and isinstance(sets_cols, dict):
        # run_sets.R's own column descriptions (sets_provenance.json), read, never restated
        L += ["", "`Sets_*.csv` (protein-set tests; `sets_provenance.json` `columns`):"]
        L += [f"- `{c}` — {d}" for c, d in sets_cols.items()]
    if files.get("qc_di"):
        hdr = _csv_header(files["qc_di"])
        L += ["", "`QC_detected_vs_inferred.csv`: "
              + "; ".join(f"`{c}` = {QC_COLUMNS.get(c, 'no description recorded')}" for c in hdr)]

    dm = de.get("detection_matrix") if isinstance(de.get("detection_matrix"), dict) else {}
    if files.get("det_matrix"):
        L += ["", f"`{os.path.basename(files['det_matrix'])}`: `Protein.Group`, then one column "
                  "per sample (the rows and columns of `Expression_Matrix.csv`). Values: "
                  f"{dm.get('values') or NOT_RECORDED}; 0 = {dm.get('zero_means') or NOT_RECORDED}"
                  + (f"; empty = {dm['na_means']}" if dm.get("na_means") else "")
                  + " (de_provenance.json `detection_matrix`)."]

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
            with open(files["qc_di"], newline="", encoding="utf-8") as fh:
                pct = [float(r["PctInferred"]) for r in csv.DictReader(fh)
                       if r.get("PctInferred") not in (None, "", "NA")]
        except (OSError, ValueError, KeyError):
            pct = []
        span = f" ({min(pct):g}–{max(pct):g}% of values per sample)" if pct else ""
        L.append(f"- **Inferred is not measured.** `Expression_Matrix.csv` has a value in "
                 f"cells where no precursor was observed{span}; "
                 + (f"`{os.path.basename(files['det_matrix'])}` marks every cell, "
                    if files.get("det_matrix") else "")
                 + "`QC_detected_vs_inferred.csv` "
                 "counts them per sample and `PropObs` per protein. Never report an inferred "
                 "value as a detection, and never use non-empty cell counts as depth — a "
                 "\"0% missing\" in the audit reflects this filled matrix, not detection.")
    fm = f["fm"]
    cont = de.get("contaminants") if isinstance(de.get("contaminants"), dict) else {}
    if fm.get("n_contaminants_appended") or fm.get("n_contaminants_already_present") or cont:
        n = (fm.get("n_contaminants_appended") or 0) + (fm.get("n_contaminants_already_present")
                                                        or 0)
        tag = fm.get("diann_cont_quant_exclude") or cont.get("tag")
        db = (f"the database holds {n} contaminant sequences"
              + (f" ({fm.get('contaminant_set')} set)" if fm.get("contaminant_set") else "")
              + (f", protein IDs starting `{tag}`" if tag else "") + ". " if n else "")
        de_step = _de_contaminants(de) if de else None
        L.append(f"- **Contaminants:** {db}"
                 + (f"{de_step} " if de_step else "")
                 + (f"A `{tag or CONT_TAG}` protein in the DE tables is contamination, not "
                    "biology. " if cont.get("policy") == "kept" or (de and not cont) else "")
                 + "A protein you expect but cannot find may be listed under a contaminant "
                   "entry.")
        if cont.get("database_risk") is True and cont.get("database_note"):
            L.append(f"  - **Caveat:** {cont['database_note']}")
    sf = f.get("submission_facts")
    if sf:
        if sf.get("error"):
            L.append(f"- **The CoreOmics submission** (`input/submission.json`) could not be read "
                     f"({sf['error']}); do not describe the samples from memory.")
        else:
            who = {"lab": "the submitting lab prepared the samples"
                          + (" and sent peptides" if sf["peptides"] else "")
                          + ": do not describe extraction, reduction, alkylation or digestion as "
                            "work the Core did",
                   "core": "the Core prepared the samples; the protocol belongs to the Methods "
                           "(`methods.md`)"}.get(sf["who"], f"who prepared the samples is not "
                                                            f"clear ({sf['why']}): say so, do "
                                                            f"not guess")
            L.append("- **The samples, in the submitter's words:** `input/submission.json` is the "
                     "submission as recorded (its `source`: the CoreOmics form, or facts the user "
                     "gave). Describe the samples only as it does and add no detail it does not "
                     f"state; {who}.")
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
    if f.get("sets_files") and (f.get("sets") or {}).get("reading"):
        L.append("- **Protein-set tests** (`Sets_*.csv`): read them by run_sets.R's rules "
                 "(`sets_provenance.json` `reading`):")
        L += [f"  - {r}" for r in f["sets"]["reading"]]
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
    sub_notes = (f.get("submission_facts") or {}).get("notes") or []
    if sub_notes:
        L.append("- CoreOmics submission (`input/submission.json`, submission_report.py):")
        L += [f"  - {n}" for n in sub_notes]
    if not (f["audit_overall"] or f["audit_notes"] or f["quality_notes"] or sub_notes):
        L.append("- No AUDIT or SAMPLE_QUALITY record in this folder.")

    L += ["", "## Reproduce", ""]
    L.append(f"- The DE: `{rel(files['repro_R'])}` (plain R, every value literal; it reads "
             f"`{os.path.basename(str(de.get('input') or 'the search report'))}`)."
             if files.get("repro_R") else "- The DE as R: not in this folder.")
    L.append(f"- Search + DE: `{rel(files['reproduce_md'])}` — needs the raw data, which are on "
             "HIVE (below), not here." if files.get("reproduce_md") else
             "- The full run: no REPRODUCE.md in this folder.")

    L += ["", "## Where this lives on HIVE", "", locations_table(f)]
    if files.get("podcast"):
        L += ["", *make_podcast.agents_md_lines(f["p"]["output_dir"], f["p"]["session_dir"])]
    if f.get("feedback"):                        # a Core run: the Core's survey
        import core_submission
        L += ["", "**Feedback:** the UC Davis Proteomics Core's survey for this project is "
                  f"{core_submission.feedback_url('readme', f['feedback']['prot'])} -- when the "
                  "person you are helping has comments on the report or the podcast, point them "
                  "to it; never fill it in for them."]
    if not for_delivery:
        L += ["", *reviewing_lines(f)]

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


def reviewing_lines(f):
    """AGENTS.md's checklist for an AI reviewing the analysis (reproducibility, not the biology),
    naming the records this session has: the decisions log, the saved conversation
    (save_transcript.py -- Core-internal, never delivered), commands.log."""
    files, rel = f["files"], f["rel"]
    have = [f"`{rel(files['decisions'])}` (the decisions, and why)" if files.get("decisions")
            else None,
            f"`{rel(files['conversation'])}` (the conversation, readable; the redacted "
            "transcripts are beside it -- Core-internal, never delivered)"
            if files.get("conversation") else None,
            f"`{rel(files['commands_log'])}` (every command run)" if files.get("commands_log")
            else None]
    have = [x for x in have if x]
    return ["## Reviewing this analysis", "",
            "For an AI asked to review how this analysis was done (a reproducibility check, not "
            "the biology). Read " + ("; ".join(have) if have else
                                     "`MANIFEST.txt` and the provenance files (no decisions log, "
                                     "conversation or commands.log is in this folder)")
            + ", then check:", "",
            "- **Groups and contrasts:** were the groups and contrasts analysed the ones the user "
            "confirmed? Compare `input/conditions.csv` and the contrasts in "
            "`output/tables/de_provenance.json` with what the user agreed to.",
            "- **Every number traces to a command:** does each count, threshold and result in "
            "the report trace to a logged command (`logs/commands.log`) and a table in "
            "`output/tables/`?",
            "- **Nothing unrecorded:** was any step skipped, changed or re-run without being "
            "recorded -- a command in the conversation that is not in `commands.log`, a "
            "parameter that differs from `search_provenance.json` / `de_provenance.json`, a "
            "re-run with no note?",
            "- **No warning ignored:** were any warnings ignored -- `[WARN]`/`WARNING` lines in "
            "the conversation and logs, `[SKIPPED]` lines in `MANIFEST.txt`, AUDIT.md findings "
            "the report does not mention?",
            "",
            "Say what you checked and what you could not (a record that is missing is a "
            "finding, not a pass)."]


def _render_html(md, title):
    import make_deposit                               # the one Markdown -> HTML renderer
    return make_deposit.md_to_html(md, title, semantic=True)


def _write_whole(path, text):
    """`text` into `path` whole or not at all: <path>.part, then os.replace (scratch_files.py
    keeps a stray .part out of the zip and the catalog)."""
    part = path + ".part"
    try:
        with open(part, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(part, path)
    finally:
        if os.path.exists(part):
            os.remove(part)


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

    # Each document is rendered to a string first, then written to <name>.part and renamed over
    # the old one: a render that fails leaves the previous file whole (or none) -- never the
    # 0-byte README.html that an open-then-render left, and zipped, before.
    def _md():
        _write_whole(f["p"]["readme"], readme_md(f))
        written["readme"] = f["p"]["readme"]

    def _html():
        path = os.path.join(sd, "README.html")
        _write_whole(path, _render_html(readme_md(f, for_html=True), f["name"]))
        written["readme_html"] = path

    def _agents():
        path = os.path.join(sd, "AGENTS.md")
        _write_whole(path, agents_md(f))
        written["agents"] = path

    # AGENTS.md first, so the README's link to it is to a file that exists
    section("AGENTS.md (guide for an AI agent)", _agents)
    if written["agents"]:
        f["files"]["agents"] = written["agents"]
    section("README.md", _md)
    section("README.html (for collaborators: open this)", _html)
    return written
