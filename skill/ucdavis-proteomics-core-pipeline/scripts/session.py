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
        tables/                 # DE_*.csv, methods.txt, sessionInfo.txt, de_provenance.json,
                                #   reproducibility_log.R (the analysis as plain R), QC
        figures/                # plots (reserved)
        reproducibility/        # the full reproducibility bundle
        AI_Analysis_Report.md   # the interpretation
        OUTPUT_FILES.md         # catalog of every file
        methods.md / .docx      # publication Methods (ensured at finalize)
        DATA_SUBMISSION/        # PRIDE/MassIVE deposit package (written at finalize)
      MANIFEST.txt              # [OK]/[SKIPPED] log of what finalize produced, and why not
      scripts/                  # copy of the skill scripts actually used (self-contained)
      logs/                     # commands.log + engine logs

Two subcommands:

  # at the start — make the folders, get paths.
  #   default (results live with the raw data):
  python3 session.py init --name "HeLa QC DIA" --raw /data/HeLaQC/*.d
  #   central location instead (user chose Documents / a custom folder):
  python3 session.py init --name "HeLa QC DIA" --raw /data/HeLaQC/*.d --base ~/Documents/DataAnalysis
  #   -> prints JSON with every canonical path + "placement"; route later steps into them

  # at the end — ensure the publication Methods, write the deposit package
  # (output/DATA_SUBMISSION, see make_deposit.py), README + MANIFEST.txt, optionally zip,
  # then log the run (record_run.py) and post "analysis complete" to the Core's Slack channel
  # (notify_slack.py)
  python3 session.py finalize --dir <session_dir> [--zip] [--no-deposit] [--no-notify]

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
        "fasta_meta": os.path.join(d, "input", "search.fasta.meta.json"),
        "workflow_dir": os.path.join(d, "input", "wf"),
        "workflow_manifest": os.path.join(d, "input", "wf", "workflow.manifest.json"),
        "raw_list": os.path.join(d, "input", "raw_files.txt"),
        "output_dir": os.path.join(d, "output"),
        "search_out": os.path.join(d, "output", "search"),
        "search_prov": os.path.join(d, "output", "search", "search_provenance.json"),
        "de_dir": os.path.join(d, "output", "tables"),
        "figures_dir": os.path.join(d, "output", "figures"),
        "repro_dir": os.path.join(d, "output", "reproducibility"),
        "analysis_report": os.path.join(d, "output", "AI_Analysis_Report.md"),
        "analysis_prompt": os.path.join(d, "output", "ANALYSIS_PROMPT.md"),
        "output_files_md": os.path.join(d, "output", "OUTPUT_FILES.md"),
        "methods_md": os.path.join(d, "output", "methods.md"),
        "methods_docx": os.path.join(d, "output", "methods.docx"),
        "deposit_dir": os.path.join(d, "output", "DATA_SUBMISSION"),
        "manifest_txt": os.path.join(d, "MANIFEST.txt"),
        "scripts_dir": os.path.join(d, "scripts"),
        "logs_dir": os.path.join(d, "logs"),
        "commands_log": os.path.join(d, "logs", "commands.log"),
    }


def read_raw_list(session_dir):
    """The raw file paths recorded for a session, in the order init wrote them ([] if none)."""
    rl = paths_for(session_dir)["raw_list"]
    if not os.path.exists(rl):
        return []
    with open(rl) as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


def raw_set(session_dir):
    """The set of raw file paths recorded for a session (for same-dataset detection)."""
    return set(read_raw_list(session_dir))


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

        # Symlink the raw files next to the results. A path in a text file is easy to
        # lose track of; a directory you can open is not. These are links, not copies —
        # the directory costs a few KB. finalize --zip skips it (see do_finalize), or a
        # .d cohort would turn a 1.6 GB archive into 16 GB.
        link_dir = os.path.join(p["output_dir"], "raw_data")
        os.makedirs(link_dir, exist_ok=True)
        for r in raws:
            src = os.path.abspath(r.rstrip("/"))
            dst = os.path.join(link_dir, os.path.basename(src))
            try:
                if os.path.islink(dst) or os.path.exists(dst):
                    os.remove(dst) if os.path.islink(dst) else None
                os.symlink(src, dst)
            except OSError:
                pass   # e.g. Windows without developer mode — the text list still has it

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


def _append_manifest(path, level, name, note):
    """One more line in MANIFEST.txt, in make_deposit.Manifest's layout: [OK] made or sent,
    [SKIPPED] attempted and failed, [INFO] a notice that is not an export part (the orchestrator
    does not relay it as missing). The note is flattened to one line -- a \r or \n in an error
    would otherwise start a line that is not a MANIFEST entry."""
    note = " ".join(str(note).replace("\r", " ").replace("\n", " ").split())
    if len(note) > 200:
        note = note[:197] + "..."
    with open(path, "a") as fh:
        fh.write(f"{'[' + level + ']':<10}{name:<50} -- {note}\n")


def _finish_hooks(a, session_dir, zip_path):
    """After the zip: log the run in the Core's run log (record_run.py analysis-done, when this
    install has it), then post "analysis complete" to the Core's Slack channel -- in that order,
    so the post can say whether the run was logged. Never fatal (notify_slack.py rule 1).
    Returns {"run_log", "slack": {"sent", "level", "detail"}, "manifest": [(level, part, note)]};
    the words are notify_slack's (run_log_manifest / slack_manifest), which never name a path or
    a group -- the zip goes to collaborators."""
    try:
        import notify_slack
    except Exception as e:                      # recorded, never swallowed
        why = f"notify_slack.py could not be loaded: {type(e).__name__}"
        sys.stderr.write(f"[session] {why}\n")
        return {"run_log": None, "slack": {"sent": False, "level": "SKIPPED", "detail": why},
                "manifest": [("SKIPPED", "Core run log + Slack notification", why)]}
    run_log = notify_slack.record_run("analysis-done", session=session_dir)
    rl = notify_slack.run_log_manifest(run_log)
    if a.no_notify:
        sent, detail = False, "--no-notify was given"
    else:
        sent, detail = notify_slack.analysis_done(session_dir, zip_path, run_log=run_log)
    sl = notify_slack.slack_manifest(sent, detail)
    # The terminal (the Core member's own) gets the full text; MANIFEST.txt, which travels in the
    # zip to collaborators, gets only the reason and the kind of failure.
    full = (run_log or {}).get("detail") if (run_log or {}).get("error") else None
    sys.stderr.write(f"[session] run log: {full or rl[2]}\n[session] slack: {sl[2]}\n")
    return {"run_log": run_log, "slack": {"sent": bool(sent), "level": sl[0], "detail": sl[2]},
            "manifest": [rl, sl]}


def _zip_manifest(zip_path, arcname, text, pre_hook):
    """Put MANIFEST.txt into the zip LAST. If that fails, retry once with the manifest as it was
    before the run-log/Slack lines: a zip with no MANIFEST.txt at all would be a regression
    (CLAUDE.md rule 4 -- a missing part must be visible). Returns what happened, for the result."""
    import zipfile
    try:
        with zipfile.ZipFile(zip_path, "a", zipfile.ZIP_DEFLATED) as z:
            z.writestr(arcname, text)
        return "added"
    except Exception as e:
        first = f"{type(e).__name__}: {e}"
    if pre_hook is None:
        return f"MISSING: {first}; no earlier copy of MANIFEST.txt to fall back on"
    try:
        with zipfile.ZipFile(zip_path, "a", zipfile.ZIP_DEFLATED) as z:
            z.writestr(arcname, pre_hook)
        return f"added without the run-log/Slack lines ({first})"
    except Exception as e:
        return f"MISSING: {first}; retry: {type(e).__name__}: {e}"


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

    # Publication Methods + the repository-deposit package (make_deposit.py). Every part is
    # recorded [OK] or [SKIPPED] -- with the reason -- in MANIFEST.txt at the session root, the
    # top of the zip: a part that could not be made is visible, never silently absent
    # (CLAUDE.md rule 4). Nothing here may stop finalize from writing the README and the zip.
    methods_md, deposit = None, None
    try:
        import make_deposit
        man = make_deposit.Manifest()
    except Exception as e:                      # recorded below, not swallowed
        make_deposit, man = None, None
        import_error = f"make_deposit.py could not be loaded: {type(e).__name__}: {e}"
    if man is not None:
        try:
            methods_md = make_deposit.ensure_methods(p["session_dir"], man)
        except Exception as e:
            man.skip("Publication methods (output/methods.md)", f"{type(e).__name__}: {e}")
        if a.no_deposit:
            man.skip("Deposit package (output/DATA_SUBMISSION)", "--no-deposit was given")
        else:
            try:
                deposit = make_deposit.build(p["session_dir"], man, methods_md)
            except Exception as e:
                man.skip("Deposit package (output/DATA_SUBMISSION)", f"{type(e).__name__}: {e}")
        man.write(p["manifest_txt"], "Session export manifest")
    else:
        with open(p["manifest_txt"], "w") as fh:
            fh.write("Session export manifest\n=======================\n"
                     f"[SKIPPED] {'Publication methods + deposit package':<50} -- "
                     f"{import_error}\n")
    has_deposit = os.path.isfile(os.path.join(p["deposit_dir"], "HOW_TO_SUBMIT.md"))
    methods_rel = (os.path.relpath(methods_md, p["session_dir"]) if methods_md else None)

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
        "output/tables/         DE results (DE_*.csv), methods.txt, sessionInfo.txt, de_provenance.json,",
        "                       reproducibility_log.R  <-- the analysis as plain R",
        "output/figures/        plots",
        "output/reproducibility/ pinned bundle for re-running the search too (reproduce.sh, env lock, checksums)",
        "output/AI_Analysis_Report.md   the biological interpretation (read this first)",
        "output/OUTPUT_FILES.md         catalog of every file",
        "output/methods.md (+ .docx)    publication Methods: LC-MS, search, database, DE, grant acknowledgment",
        "output/DATA_SUBMISSION/        deposit the data in PRIDE / MassIVE -- start with HOW_TO_SUBMIT.md (.html)",
        "MANIFEST.txt                   what this export contains, and anything skipped (with the reason)",
        "scripts/               copy of the skill scripts used",
        "logs/                  commands.log + engine logs",
        "```",
        "", "## Reproduce",
        "**The analysis, as code:** `output/tables/reproducibility_log.R` — the whole "
        "differential-expression analysis in plain R with every value written out. "
        "Read it to see what was done, or `Rscript` it to redo it with nothing but R "
        "and limpa/limma.",
        "",
        "**The whole run, pinned** (search included, takes hours): "
        "`output/reproducibility/REPRODUCE.md`.",
        "",
        f"The DE results are {', '.join(de_files) if de_files else '(none found)'}.",
        "", "## Methods",
        (f"**For the paper:** `{methods_rel}` (Word: the `.docx` beside it) — LC-MS "
         "acquisition, the database search (engine, pinned version, parameters, FDR), the "
         "sequence database and contaminants, the differential-expression analysis, and the "
         "UC Davis instrument-grant acknowledgment. Resolve every `[... — confirm]` tag before "
         "publishing." if methods_rel else
         "**For the paper:** no publication Methods could be written — `MANIFEST.txt` says "
         "why; run `scripts/make_methods.py` where the raw files are readable."),
        "",
        "The DE step's own record is `output/tables/methods.txt`; the interpretation is "
        "`output/AI_Analysis_Report.md`.", "",
        "## Deposit the data (PRIDE / MassIVE)",
        ("Journals ask for the raw data in a public repository. Everything to do that is in "
         "`output/DATA_SUBMISSION/` — start with **`HOW_TO_SUBMIT.md`** (or double-click "
         "`HOW_TO_SUBMIT.html`): a pre-filled SDRF sample sheet, the protocol texts, the list "
         "of files to upload, and a SLURM script that packs and checksums the raw files on "
         "HIVE." if has_deposit else
         "The deposit package was not written — `MANIFEST.txt` says why."), "",
        "## What is in this export",
        "`MANIFEST.txt` lists every part of the Methods and deposit package as [OK], or "
        "[SKIPPED] with the reason.", "",
    ]
    with open(p["readme"], "w") as fh:
        fh.write("\n".join(lines) + "\n")

    result = {"session_dir": p["session_dir"], "readme": p["readme"], "de_files": de_files,
              "methods": methods_md, "manifest": p["manifest_txt"],
              "deposit": deposit, "skipped": (man.n_skipped if man is not None else 1)}

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
        # shutil.make_archive FOLLOWS symlinks, so output/raw_data (links to the raw
        # cohort) would be dereferenced and inlined — turning a 1.6 GB archive into
        # 16 GB for a 9-file .d study. Zip by hand and skip that one directory; the
        # raw paths are still recorded in input/raw_files.txt. The deposit package's
        # upload_staging/ (the .d archives prepare_upload.sbatch builds) is raw data too.
        import zipfile
        sdir = p["session_dir"]
        base = os.path.basename(sdir)
        skips = {os.path.abspath(os.path.join(sdir, "output", "raw_data")):
                 "output/raw_data (symlinks to raw files)",
                 os.path.abspath(os.path.join(p["deposit_dir"], "upload_staging")):
                 "output/DATA_SUBMISSION/upload_staging (raw archives staged for upload)"}
        archive = sdir + ".zip"
        excluded = {label: 0 for label in skips.values()}
        manifest_txt = os.path.abspath(p["manifest_txt"])
        # DIA-NN's per-run .quant intermediates: ~30 MB each (measured on HIVE, 28-34 MB), so
        # ~3 GB for a 100-file cohort, and nothing a reader of the zip can use. The single-shot
        # search writes them to <search out>/quant (its --temp, inside the session since #79);
        # the 5-step chain to quant_step2/ and quant_step4/. Covered: every *.quant file anywhere
        # in the session, plus the whole `quant` directory directly under any search out dir (the
        # session's output/search, or a folder holding a search_provenance.json). They stay on
        # disk where they are. The key is short on purpose: the finalize JSON carries it and the
        # agent relays it; which files it covers is said here.
        quant_label = "DIA-NN .quant intermediates (kept on disk where they are)"
        n_quant = 0
        # The in-silico predicted library (step1.predicted.speclib, <lib>.predicted.speclib):
        # 687 MB for one 15-file Lumos session on HIVE. It is regenerated exactly from the FASTA,
        # the pinned engine and the params the zip already holds, so it carries nothing a reader
        # needs. Empirical libraries (.parquet/.speclib without "predicted") stay in.
        # "where they are": a predicted library may sit anywhere in the session, not only in
        # the search out dir (a seeded chain, a Radiant library dir).
        speclib_label = ("predicted spectral libraries (*.predicted.speclib anywhere in the "
                         "session; rebuilt from the FASTA + params, kept on disk where they are)")
        n_speclib = 0

        def is_search_out(d):
            return (os.path.abspath(d) == os.path.abspath(p["search_out"])
                    or os.path.isfile(os.path.join(d, "search_provenance.json")))

        with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as z:
            for root, dirs, files in os.walk(sdir):
                if os.path.abspath(root) in skips:
                    excluded[skips[os.path.abspath(root)]] = len(files) + len(dirs)
                    dirs[:] = []
                    continue
                if is_search_out(root) and "quant" in dirs:
                    dirs.remove("quant")
                    n_quant += sum(len(fs) for _, _, fs in os.walk(os.path.join(root, "quant")))
                for fn in files:
                    full = os.path.join(root, fn)
                    if fn.endswith(".quant"):
                        n_quant += 1
                        continue
                    if fn.endswith(".predicted.speclib"):
                        n_speclib += 1
                        continue
                    if os.path.islink(full) or os.path.abspath(full) == manifest_txt:
                        continue                 # MANIFEST.txt goes in last, below
                    z.write(full, os.path.join(base, os.path.relpath(full, sdir)))
        excluded[quant_label] = n_quant
        excluded[speclib_label] = n_speclib
        result["zip"] = archive
        result["zip_excluded"] = excluded

    # The run log and Slack go once everything they point at exists -- after the zip. Their
    # outcomes are the last two lines of MANIFEST.txt, so the zip's copy (added now) says them.
    try:
        with open(p["manifest_txt"]) as fh:
            pre_hook = fh.read()
    except OSError:
        pre_hook = None
    hooks = _finish_hooks(a, p["session_dir"], result.get("zip"))
    result["run_log"], result["slack"] = hooks["run_log"], hooks["slack"]
    for level, part, note in hooks["manifest"]:
        try:
            _append_manifest(p["manifest_txt"], level, part, note)
        except Exception as e:
            sys.stderr.write(f"[session] could not record '{part}' in MANIFEST.txt: {e}\n")
    if a.zip:
        try:
            with open(p["manifest_txt"]) as fh:
                text = fh.read()
        except OSError:
            text = pre_hook
        result["zip_manifest"] = _zip_manifest(
            result["zip"], os.path.join(os.path.basename(p["session_dir"]), "MANIFEST.txt"),
            text if text is not None else "", pre_hook)
        if result["zip_manifest"] != "added":
            sys.stderr.write(f"[session] MANIFEST.txt in the zip: {result['zip_manifest']}\n")

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
        # Spell the database out: a re-analysis that switched organism, database type,
        # or contaminant set must show up as a difference, not hide behind one count.
        "fasta": (f"{fi.get('organism') or '?'} "
                  f"({fi.get('proteome') or fi.get('source', '?')}), "
                  f"{fi.get('content_used') or '?'}, "
                  f"{fi.get('n_proteome', '?')} proteome + "
                  f"{(fi.get('n_contaminants_appended') or 0) + (fi.get('n_contaminants_already_present') or 0)}"
                  f" contaminants "
                  f"[{fi.get('contaminant_set') or ('universal' if fi.get('n_contaminants_appended') else 'none')}]"
                  + (f", UniProt {fi['uniprot_release']}" if fi.get("uniprot_release") else "")
                  ) if fi else "?",
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
    f.add_argument("--no-deposit", action="store_true",
                   help="skip the repository-deposit package (output/DATA_SUBMISSION); the "
                        "publication Methods are still ensured")
    f.add_argument("--no-notify", action="store_true",
                   help="do not post 'analysis complete' to the Core's Slack channel (same as "
                        "SKILL_SLACK=0; see references/notifications.md)")
    f.set_defaults(func=do_finalize)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
