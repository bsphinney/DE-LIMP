#!/usr/bin/env python3
"""
session.py  --  Package a run's inputs and outputs into a tidy, browsable session
directory. By DEFAULT the session is created IN THE FOLDER WITH THE RAW DATA being
analyzed; pass --base to put it in a central location instead (e.g. the user's
Documents). The orchestrator asks the user which they want (see SKILL.md).

    <YYYY-MM-DD>_<DescriptiveName>/    # next to the raw files, or under <base>/sessions/
      README.md                 # what this analysis was + where everything is
      input/                    # conditions, FASTA, params, workflow manifest, raw-file list
      output/
        search/                 # the normalized search report (+ engine logs)
        tables/                 # DE_*.csv, methods.txt, sessionInfo.txt, de_provenance.json, QC
        figures/                # plots (reserved)
        reproducibility/        # the full reproducibility bundle
        AI_Analysis_Report.md   # the interpretation
        OUTPUT_FILES.md         # catalog of every file
      scripts/                  # copy of the skill scripts actually used (self-contained)
      logs/                     # commands.log + engine logs

Two subcommands:

  # at the start — make the folders, get paths.
  #   default (results live with the raw data):
  python3 session.py init --name "HeLa QC DIA" --raw /data/HeLaQC/*.d
  #   central location instead (user chose Documents / a custom folder):
  python3 session.py init --name "HeLa QC DIA" --raw /data/HeLaQC/*.d --base ~/Documents/DataAnalysis
  #   -> prints JSON with every canonical path + "placement"; route later steps into them

  # at the end — write README, catalog outputs, optionally zip
  python3 session.py finalize --dir <session_dir> [--zip]

Raw MS files are NOT copied (they're huge and live elsewhere) — their paths are
recorded in input/raw_files.txt instead.
"""
import sys, os, json, re, glob, shutil, argparse, datetime

SUBDIRS = ["input", "output", "output/search", "output/tables", "output/figures",
           "output/reproducibility", "scripts", "logs"]


def slugify(name):
    s = re.sub(r"[^A-Za-z0-9]+", "_", name.strip()).strip("_")
    return s or "proteomics_run"


def paths_for(session_dir):
    d = os.path.abspath(session_dir)
    return {
        "session_dir": d,
        "readme": os.path.join(d, "README.md"),
        "input_dir": os.path.join(d, "input"),
        "conditions": os.path.join(d, "input", "conditions.csv"),
        "fasta": os.path.join(d, "input", "search.fasta"),
        "workflow_dir": os.path.join(d, "input", "wf"),
        "raw_list": os.path.join(d, "input", "raw_files.txt"),
        "output_dir": os.path.join(d, "output"),
        "search_out": os.path.join(d, "output", "search"),
        "de_dir": os.path.join(d, "output", "tables"),
        "figures_dir": os.path.join(d, "output", "figures"),
        "repro_dir": os.path.join(d, "output", "reproducibility"),
        "analysis_report": os.path.join(d, "output", "AI_Analysis_Report.md"),
        "analysis_prompt": os.path.join(d, "output", "ANALYSIS_PROMPT.md"),
        "output_files_md": os.path.join(d, "output", "OUTPUT_FILES.md"),
        "scripts_dir": os.path.join(d, "scripts"),
        "logs_dir": os.path.join(d, "logs"),
        "commands_log": os.path.join(d, "logs", "commands.log"),
    }


def raw_set(session_dir):
    """The set of raw file paths recorded for a session (for same-dataset detection)."""
    rl = paths_for(session_dir)["raw_list"]
    if not os.path.exists(rl):
        return set()
    return {ln.strip() for ln in open(rl) if ln.strip() and not ln.startswith("#")}


def do_find_prior(a):
    """Scan existing sessions for ones covering the same raw files (same dataset)."""
    mine, raw_dirs = set(), set()
    for pat in (a.raw or []):
        hits = [os.path.abspath(p.rstrip("/")) for p in glob.glob(pat)]
        if hits:
            mine.update(hits)
            raw_dirs.update(os.path.dirname(h) for h in hits)
        else:
            mine.add(pat)
    # look where sessions can live: alongside the raw data (default), and in a
    # central --base/sessions if the user used one. Include reanalysis subfolders.
    roots = list(raw_dirs)
    if a.base:
        roots.append(os.path.join(os.path.abspath(os.path.expanduser(a.base)), "sessions"))
    candidates = []
    for root in roots:
        candidates += glob.glob(os.path.join(root, "*"))
        candidates += glob.glob(os.path.join(root, "*", "reanalysis", "*"))
    hits = []
    for c in candidates:
        if not os.path.isdir(c):
            continue
        rs = raw_set(c)
        if not rs or not mine:
            continue
        inter = mine & rs
        if inter:
            hits.append({"session": c, "overlap": len(inter),
                         "of_mine": len(mine), "of_theirs": len(rs),
                         "same_dataset": inter == mine == rs})
    hits.sort(key=lambda h: h["overlap"], reverse=True)
    print(json.dumps({"query_raw_count": len(mine), "matches": hits,
                      "suggestion": ("re-analysis of " + hits[0]["session"]) if hits else
                                    "no prior session covers these raw files — this is a fresh analysis"},
                     indent=2))


def _resolve_raws(patterns):
    raws = []
    for pat in (patterns or []):
        raws.extend(sorted(glob.glob(pat)) or [pat])
    return raws


def _raw_dir(raws):
    """The directory that contains the raw files (their common parent)."""
    if not raws:
        return None
    dirs = [os.path.dirname(os.path.abspath(r.rstrip("/"))) for r in raws]
    try:
        return os.path.commonpath(dirs)
    except ValueError:
        return dirs[0]


def do_init(a):
    date = a.date or datetime.date.today().isoformat()
    slug = slugify(a.name)
    raws = _resolve_raws(a.raw)
    raw_dir = _raw_dir(raws)

    # WHERE the results go (the orchestrator asks the user; see SKILL.md):
    #   --reanalysis-of <prior>  -> nested under the original
    #   --base <path>            -> a central location the user chose (e.g. ~/Documents/DataAnalysis)
    #   (neither, with --raw)    -> DEFAULT: in the folder with the raw data being analyzed
    if a.reanalysis_of:
        prior = os.path.abspath(os.path.expanduser(a.reanalysis_of))
        if not os.path.isdir(prior):
            sys.exit(f"--reanalysis-of: prior session not found: {prior}")
        session_dir = os.path.join(prior, "reanalysis", f"{date}_{slug}")
        placement = "reanalysis"
    elif a.base:
        base = os.path.abspath(os.path.expanduser(a.base))
        session_dir = os.path.join(base, "sessions", f"{date}_{slug}")
        placement = "central"
    elif raw_dir:
        session_dir = os.path.join(raw_dir, f"{date}_{slug}")
        placement = "with-raw-data"
    else:
        session_dir = os.path.join(os.path.abspath("."), f"{date}_{slug}")
        placement = "cwd"
    for sd in SUBDIRS:
        os.makedirs(os.path.join(session_dir, sd), exist_ok=True)
    p = paths_for(session_dir)

    # self-contained: copy the skill scripts that ran this analysis
    skill_scripts = os.path.dirname(os.path.abspath(__file__))
    try:
        for f in glob.glob(os.path.join(skill_scripts, "*")):
            if os.path.isfile(f):
                shutil.copy2(f, os.path.join(p["scripts_dir"], os.path.basename(f)))
    except Exception as e:
        sys.stderr.write(f"[session] could not copy skill scripts: {e}\n")

    # record raw file locations (not the files themselves)
    if raws:
        with open(p["raw_list"], "w") as fh:
            fh.write("# Raw MS files used in this analysis (not copied — too large).\n")
            for r in raws:
                fh.write(os.path.abspath(r.rstrip("/")) + "\n")

    # record the parent when this is a re-analysis
    parent = os.path.abspath(os.path.expanduser(a.reanalysis_of)) if a.reanalysis_of else None
    if parent:
        with open(os.path.join(session_dir, ".reanalysis_of"), "w") as fh:
            fh.write(parent + "\n")

    # starter README (finalize fills in results)
    with open(p["readme"], "w") as fh:
        fh.write(f"# {a.name}\n\n- Date: {date}\n- Status: in progress\n")
        if parent:
            fh.write(f"- **Re-analysis of:** `{parent}` — see `DIFFERENCES.md` (written at finalize) "
                     "for exactly what changed.\n")
        fh.write("\nLayout: `input/` (conditions, FASTA, params), `output/` "
                 "(search, tables, figures, reproducibility, report), `scripts/`, `logs/`.\n")

    print(json.dumps({"created": session_dir, "date": date, "name": a.name,
                      "placement": placement, "reanalysis_of": parent, "paths": p}, indent=2))


def _load(path):
    try: return json.load(open(path))
    except Exception: return None


def do_finalize(a):
    p = paths_for(a.dir)
    if not os.path.isdir(p["session_dir"]):
        sys.exit(f"session dir not found: {p['session_dir']}")

    # tidy: move any loose tables/figures left in output/ root into their subdirs
    for f in glob.glob(os.path.join(p["output_dir"], "*")):
        if not os.path.isfile(f):
            continue
        ext = f.lower().rsplit(".", 1)[-1]
        base = os.path.basename(f)
        if base in ("AI_Analysis_Report.md", "ANALYSIS_PROMPT.md", "OUTPUT_FILES.md"):
            continue
        if ext in ("csv", "tsv"):
            shutil.move(f, os.path.join(p["de_dir"], base))
        elif ext in ("png", "svg", "pdf", "jpg", "jpeg"):
            shutil.move(f, os.path.join(p["figures_dir"], base))

    # gather run facts for the README
    manifest = _load(os.path.join(p["repro_dir"], "run_manifest.json")) or {}
    prov = _load(os.path.join(p["de_dir"], "de_provenance.json")) or {}
    wfman = _load(os.path.join(p["workflow_dir"], "workflow.manifest.json")) or {}
    engine = (manifest.get("engine") or (wfman.get("engine", {}) or {}).get("name") or "?")
    eng_ver = (wfman.get("engine", {}) or {}).get("version", "")
    reg = manifest.get("registry") or wfman.get("registry") or {}
    method = prov.get("method") or (wfman.get("de", {}) or {}).get("method", "?")
    sig = prov.get("significant_per_contrast") or {}
    contrasts = prov.get("contrasts") or []
    q = manifest.get("query") or {}

    de_files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(p["de_dir"], "DE_*.csv")))
    lines = [
        f"# {os.path.basename(p['session_dir'])}", "",
        "Proteomics search + differential expression, run by the ucdavis-proteomics-core-pipeline skill.", "",
        "## Summary",
        f"- Organism (taxid): {q.get('organism_taxid', '?')}",
        f"- Acquisition / instrument: {q.get('acquisition', '?')} / {q.get('instrument') or '?'}",
        f"- Search engine: {engine} {eng_ver}".rstrip(),
        f"- DE method: {method}",
        f"- Contrasts: {', '.join(contrasts) if contrasts else '?'}",
    ]
    if sig:
        lines.append("- Significant proteins per contrast: "
                     + ", ".join(f"{k}={v}" for k, v in sig.items()))
    if reg.get("commit"):
        lines.append(f"- Validated workflow: {reg.get('repo','')} @ `{reg['commit']}`")
    lines += [
        "", "## Where everything is",
        "```",
        "input/                 conditions.csv, search.fasta, params, workflow manifest, raw_files.txt",
        "output/search/         normalized search report (DE input)",
        "output/tables/         DE results (DE_*.csv), methods.txt, sessionInfo.txt, de_provenance.json",
        "output/figures/        plots",
        "output/reproducibility/ full reproducibility bundle (reproduce.sh, env lock, checksums)",
        "output/AI_Analysis_Report.md   the biological interpretation (read this first)",
        "output/OUTPUT_FILES.md         catalog of every file",
        "scripts/               copy of the skill scripts used",
        "logs/                  commands.log + engine logs",
        "```",
        "", "## Reproduce",
        "See `output/reproducibility/REPRODUCE.md`. The DE results are "
        f"{', '.join(de_files) if de_files else '(none found)'}.",
        "", "## Methods",
        "See `output/tables/methods.txt` (self-describing) and `output/AI_Analysis_Report.md`.", "",
    ]
    with open(p["readme"], "w") as fh:
        fh.write("\n".join(lines) + "\n")

    result = {"session_dir": p["session_dir"], "readme": p["readme"], "de_files": de_files}

    # re-analysis: write DIFFERENCES.md vs the parent
    parent = a.reanalysis_of
    marker = os.path.join(p["session_dir"], ".reanalysis_of")
    if not parent and os.path.exists(marker):
        parent = open(marker).read().strip()
    if parent:
        diff_path = write_differences(parent, p["session_dir"])
        result["differences"] = diff_path
        result["reanalysis_of"] = os.path.abspath(os.path.expanduser(parent))

    if a.zip:
        archive = shutil.make_archive(p["session_dir"], "zip",
                                      root_dir=os.path.dirname(p["session_dir"]),
                                      base_dir=os.path.basename(p["session_dir"]))
        result["zip"] = archive

    print(json.dumps(result, indent=2))


def _facts(session_dir):
    """Pull the comparable facts of a run from its session files."""
    p = paths_for(session_dir)
    man = _load(os.path.join(p["repro_dir"], "run_manifest.json")) or {}
    prov = _load(os.path.join(p["de_dir"], "de_provenance.json")) or {}
    wf = _load(os.path.join(p["workflow_dir"], "workflow.manifest.json")) or {}
    params = sorted(glob.glob(os.path.join(p["input_dir"], "params.*"))) + \
             sorted(glob.glob(os.path.join(p["input_dir"], "*.cfg"))) + \
             sorted(glob.glob(os.path.join(p["input_dir"], "sage_config*.json")))
    ptext = ""
    for pf in params:
        if not pf.endswith(".rationale.json"):
            try: ptext = open(pf).read(); break
            except OSError: pass
    fi = (man.get("inputs") or {}).get("fasta_info") or {}
    return {
        "engine": man.get("engine") or (wf.get("engine", {}) or {}).get("name"),
        "engine_version": (wf.get("engine", {}) or {}).get("version"),
        "de_method": prov.get("method") or (wf.get("de", {}) or {}).get("method"),
        "q_cutoff": prov.get("q_cutoff"), "logfc": prov.get("logfc"), "adjp": prov.get("adjp"),
        "contrasts": prov.get("contrasts") or [],
        "fasta": f"{fi.get('source','?')} (n={fi.get('n_sequences','?')})",
        "registry_commit": (man.get("registry") or wf.get("registry") or {}).get("commit"),
        "significant_per_contrast": prov.get("significant_per_contrast") or {},
        "params_text": ptext,
        "conditions": p["conditions"],
    }


def write_differences(prior_dir, new_dir):
    prior_dir = os.path.abspath(os.path.expanduser(prior_dir))
    old, new = _facts(prior_dir), _facts(new_dir)
    L = [f"# What changed in this re-analysis", "",
         f"Re-analysis of `{prior_dir}`.", "",
         "Same raw data; the table below is exactly what differs. Unchanged settings are omitted.",
         "", "| Aspect | Original | This re-analysis |", "|---|---|---|"]
    fields = [("Search engine", "engine"), ("Engine version", "engine_version"),
              ("DE method", "de_method"), ("ID FDR (q)", "q_cutoff"),
              ("logFC threshold", "logfc"), ("adj.P threshold", "adjp"),
              ("Contrasts", "contrasts"), ("FASTA", "fasta"),
              ("Validated workflow commit", "registry_commit")]
    n_changes = 0
    for label, key in fields:
        a_v, b_v = old.get(key), new.get(key)
        if a_v != b_v:
            n_changes += 1
            L.append(f"| {label} | {a_v} | {b_v} |")
    if n_changes == 0:
        L.append("| (settings) | — | identical settings; difference is data/environment only |")

    # search-parameter text diff (mass tolerances etc.)
    if old["params_text"] and new["params_text"] and old["params_text"] != new["params_text"]:
        import difflib
        d = list(difflib.unified_diff(old["params_text"].splitlines(),
                                      new["params_text"].splitlines(),
                                      "original/params", "reanalysis/params", lineterm=""))
        L += ["", "## Search-parameter diff", "```diff", *d[:200], "```"]

    # results delta
    so, sn = old["significant_per_contrast"], new["significant_per_contrast"]
    if so or sn:
        L += ["", "## Significant-protein counts", "",
              "| Contrast | Original | This re-analysis |", "|---|---|---|"]
        for ct in sorted(set(so) | set(sn)):
            L.append(f"| {ct} | {so.get(ct, '—')} | {sn.get(ct, '—')} |")
        L += ["", "_For a protein-level comparison (overlap, concordance, logFC correlation), "
              "run `compare_analyses.R` across the two sessions' `output/tables` dirs._"]

    out = os.path.join(new_dir, "DIFFERENCES.md")
    with open(out, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    i = sub.add_parser("init", help="scaffold a session directory and print canonical paths")
    i.add_argument("--name", required=True, help="short descriptive study name")
    i.add_argument("--base", default="", help="central location for the session (e.g. ~/Documents/DataAnalysis). OMIT to put results in the folder with the raw data (default).")
    i.add_argument("--date", default="", help="YYYY-MM-DD (default: today)")
    i.add_argument("--raw", nargs="*", help="raw file paths/globs; their folder is where results go by default")
    i.add_argument("--reanalysis-of", default="", help="prior session dir; nests this run under <prior>/reanalysis/")
    i.set_defaults(func=do_init)
    fp = sub.add_parser("find-prior", help="find existing sessions covering the same raw files")
    fp.add_argument("--base", default="", help="also scan this central location's sessions/ (optional)")
    fp.add_argument("--raw", nargs="*", required=True)
    fp.set_defaults(func=do_find_prior)
    f = sub.add_parser("finalize", help="write README, tidy output/, optionally zip")
    f.add_argument("--dir", required=True, help="the session directory")
    f.add_argument("--reanalysis-of", default="", help="prior session dir to diff against (auto-detected if omitted)")
    f.add_argument("--zip", action="store_true", help="also produce <session>.zip")
    f.set_defaults(func=do_finalize)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
