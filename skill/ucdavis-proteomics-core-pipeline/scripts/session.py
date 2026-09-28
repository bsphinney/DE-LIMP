#!/usr/bin/env python3
"""
session.py  --  Package a run's inputs and outputs into a tidy, browsable session
directory. By DEFAULT the session is created IN THE FOLDER WITH THE RAW DATA being
analyzed; pass --base to put it in a central location instead (e.g. the user's
Documents). The orchestrator asks the user which they want (see SKILL.md).

    <YYYY-MM-DD>_<DescriptiveName>/    # next to the raw files, or under <base>/sessions/
      README.html               # OPEN THIS: what the analysis was, links, where it lives on HIVE
      README.md                 # the same, as text (both rendered from one source: session_docs.py)
      AGENTS.md                 # a guide to the folder for an AI agent, from the session's records
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

Subcommands:

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

  # README.md + README.html + AGENTS.md alone (e.g. for a session finalized before they existed)
  python3 session.py docs --dir <session_dir>

Raw MS files are NOT copied (they're huge and live elsewhere) — their paths are
recorded in input/raw_files.txt instead (finalize writes it from the search's record when a
session was initialised without --raw).
"""
import sys, os, json, re, glob, shutil, argparse, datetime

NOT_RECORDED = "not recorded"      # a fact the session's records do not give -- never a guess

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
        # the CoreOmics submission this session answers (submission_report.py attach)
        "session_json": os.path.join(d, "session.json"),
        "submission_record": os.path.join(d, "input", "submission.json"),
        "submission_samples": os.path.join(d, "input", "samples.tsv"),
    }


def read_raw_list(session_dir):
    """The raw file paths recorded for a session, in the order init wrote them ([] if none).

    Skill 2.7 and older wrote the file in the computer's own encoding (cp1252 on Windows, where
    even the header's dash is not UTF-8): a byte that is not UTF-8 reads as U+FFFD, never a
    UnicodeDecodeError that cost finalize the README, AGENTS.md and the Methods.
    raw_list_encoding_note() says whether any path was affected."""
    rl = paths_for(session_dir)["raw_list"]
    if not os.path.exists(rl):
        return []
    with open(rl, encoding="utf-8", errors="replace") as fh:
        return [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]


NOT_UTF8 = ("not UTF-8 (written in the computer's own encoding -- skill 2.7 and older did, and "
            "cp1252 on Windows is not UTF-8): each such byte reads as U+FFFD")


def encoding_note(path, noun="line", consequence="it is not as written"):
    """None when `path` is UTF-8 (or cannot be read); else a MANIFEST note saying so and which
    entries -- non-comment lines, as `noun` -- read with U+FFFD, and what that means."""
    try:
        with open(path, "rb") as fh:
            data = fh.read()
    except OSError:
        return None
    try:
        data.decode("utf-8")
        return None
    except UnicodeDecodeError:
        pass
    bad = [ln.strip() for ln in data.decode("utf-8", "replace").splitlines()
           if ln.strip() and not ln.startswith("#") and "\ufffd" in ln]
    return (f"{NOT_UTF8}; {len(bad)} {noun}(s) have one, so {consequence} (first: {bad[0]})"
            if bad else f"{NOT_UTF8}, only in its comment lines -- every {noun} reads as written")


def raw_list_encoding_note(session_dir):
    """encoding_note() for input/raw_files.txt: a path shown with U+FFFD is not the path on
    disk."""
    return encoding_note(paths_for(session_dir)["raw_list"], "path",
                         "they are not the path on disk -- check them against the raw data")


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
        with open(p["raw_list"], "w", encoding="utf-8") as fh:
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
        with open(os.path.join(session_dir, ".reanalysis_of"), "w", encoding="utf-8") as fh:
            fh.write(parent + "\n")

    # starter README (finalize fills in results)
    with open(p["readme"], "w", encoding="utf-8") as fh:
        fh.write(f"# {a.name}\n\n- Date: {date}\n- Status: in progress\n")
        if parent:
            fh.write(f"- **Re-analysis of:** `{parent}` — see `DIFFERENCES.md` (written at finalize) "
                     "for exactly what changed.\n")
        fh.write("\nLayout: `input/` (conditions, FASTA, params), `output/` "
                 "(search, tables, figures, reproducibility, report), `scripts/`, `logs/`.\n")

    print(json.dumps({"created": session_dir, "date": date, "name": a.name,
                      "placement": placement, "reanalysis_of": parent, "paths": p}, indent=2))


def _load(path):
    try:
        with open(path, encoding="utf-8") as fh:
            return json.load(fh)
    except Exception:
        return None


def _append_manifest(path, level, name, note):
    """One more line in MANIFEST.txt, in make_deposit.Manifest's layout: [OK] made or sent,
    [SKIPPED] attempted and failed, [INFO] a notice that is not an export part (the orchestrator
    does not relay it as missing). The note is flattened to one line -- a \r or \n in an error
    would otherwise start a line that is not a MANIFEST entry."""
    note = " ".join(str(note).replace("\r", " ").replace("\n", " ").split())
    if len(note) > 200:
        note = note[:197] + "..."
    with open(path, "a", encoding="utf-8") as fh:
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


LATE_DOCS = ("README.md", "README.html", "AGENTS.md")

# DIA-NN's in-silico predicted library (step1.predicted.speclib, <lib>.predicted.speclib): 687 MB
# for one 15-file Lumos session on HIVE. It is rebuilt exactly from the FASTA, the pinned engine
# and the params the zip already holds, so the session zip leaves it out -- wherever it sits in the
# session (a seeded chain, a Radiant library dir). ONE rule: finalize's zip and make_deposit.py
# (files_to_upload.tsv says it is not in the zip) both read it here. Empirical libraries
# (.parquet, a .speclib without "predicted") are not covered and stay in.
PREDICTED_SPECLIB = ".predicted.speclib"


def is_predicted_speclib(name):
    return os.path.basename(name).lower().endswith(PREDICTED_SPECLIB)


def _registry_lookup(session_dir, in_finalize=True):
    """(the Core run-registry folder already holding this session, or None; why not). Read-only
    and local: record_run.locate() -- never an ssh call from here."""
    # A subprocess, like every other call into record_run.py: nothing it does -- or a broken
    # copy of it -- can reach this process.
    import subprocess
    rr = os.path.join(os.path.dirname(os.path.abspath(__file__)), "record_run.py")
    if not os.path.isfile(rr):
        return None, "record_run.py is not installed here"
    try:
        r = subprocess.run([sys.executable, rr, "locate", "--session", session_dir],
                           capture_output=True, text=True, timeout=30,
                           stdin=subprocess.DEVNULL)
        got = json.loads((r.stdout.strip().splitlines() or ["{}"])[-1])
    except Exception as e:                      # reported in the README line, not swallowed
        return None, f"lookup failed ({type(e).__name__})"
    if not isinstance(got, dict) or "located" not in got:
        return None, "lookup failed (record_run.py gave no answer)"
    path = got.get("located")
    if path:
        return path, None
    if in_finalize:
        return None, ("no record found from where finalize ran -- the 'Core run log' line of "
                      "MANIFEST.txt says whether this run was recorded")
    return None, "no record found from this machine (the registry is on HIVE)"


def _write_docs(session_docs, p, man, registry, registry_note, ok_lines=True, pending=(),
                located_at=None):
    """README.md + README.html + AGENTS.md (session_docs.py). Never fatal. With a Manifest each
    file is its own [OK]/[SKIPPED] line; without one the lines come back in "manifest" (the
    [OK] ones only when ok_lines -- a re-render after the hook reports failures only)."""
    out = {"written": {}, "manifest": []}
    if session_docs is None:
        why = "session_docs.py could not be loaded"
        (man.skip("README / AGENTS.md", why) if man is not None
         else out["manifest"].append(("SKIPPED", "README / AGENTS.md", why)))
        return out
    rec = man
    if rec is None:
        class _Lines:                           # the same [OK]/[SKIPPED] outcome, as tuples
            def section(self, name, fn):
                try:
                    fn()
                    if ok_lines:
                        out["manifest"].append(("OK", name, "written"))
                    return True
                except Exception as e:
                    out["manifest"].append(("SKIPPED", name, f"{type(e).__name__}: {e}"))
                    return False
        rec = _Lines()
    try:
        out["written"] = session_docs.write_docs(p["session_dir"], rec, registry, registry_note,
                                                 pending=pending, located_at=located_at)
    except Exception as e:                      # gather() itself failed: say so, keep going
        why = f"{type(e).__name__}: {e}"
        (man.skip("README / AGENTS.md", why) if man is not None
         else out["manifest"].append(("SKIPPED", "README / AGENTS.md", why)))
        return out
    try:
        import make_report
        listed = [os.path.join(p["session_dir"], n) for n in LATE_DOCS + ("MANIFEST.txt",)]
        listed.append(p["raw_list"])
        if os.path.isfile(p["output_files_md"]) and man is not None:
            n = make_report.add_finalize_files(p["output_files_md"], listed, p["session_dir"])
            man.ok("OUTPUT_FILES.md: files written at finalize", f"{n} listed")
        elif man is not None:
            man.info("OUTPUT_FILES.md: files written at finalize",
                     "no output/OUTPUT_FILES.md (step 11 was not run)")
    except Exception as e:
        if man is not None:
            man.skip("OUTPUT_FILES.md: files written at finalize", f"{type(e).__name__}: {e}")
    return out


def _zip_docs(zip_path, session_dir, written):
    """Add README.md, README.html and AGENTS.md to the zip (after the run-log hook) -- only the
    ones this finalize wrote (`written`: the paths session_docs.write_docs returned). One it could
    not write stays out even when an older copy is on disk: that copy describes an earlier state.
    Returns "added", or what went wrong -- recorded in MANIFEST.txt, which goes in after."""
    import zipfile
    base = os.path.basename(session_dir)
    done = {os.path.abspath(w) for w in written if w}
    missing = []
    try:
        with zipfile.ZipFile(zip_path, "a", zipfile.ZIP_DEFLATED) as z:
            for n in LATE_DOCS:
                full = os.path.join(session_dir, n)
                if os.path.abspath(full) in done and os.path.isfile(full):
                    z.write(full, os.path.join(base, n))
                else:
                    missing.append(n)
    except Exception as e:
        return f"not added: {type(e).__name__}: {e}"
    return "added" if not missing else f"added; not written, so not in the zip: {', '.join(missing)}"


def report_pdf_step(p, man, print_report=None):
    """The report PDF (html_to_pdf.py): made here when step 9 ran where no browser was -- in
    hive_remote, on HIVE -- or the HTML changed since. Never fatal. A *.stale.pdf (an earlier
    report that print_report could not re-print) never reaches a collaborator: it is scratch
    (scratch_files.py), so the zip and the catalog leave it out, and once a current PDF exists it
    is deleted. `print_report` is html_to_pdf.print_report (a stand-in in tests)."""
    name = "Report PDF (output/Analysis_Report.pdf)"
    html_rep = os.path.join(p["output_dir"], "Analysis_Report.html")
    pdf_rep = os.path.join(p["output_dir"], "Analysis_Report.pdf")
    stale = os.path.splitext(pdf_rep)[0] + ".stale.pdf"
    if not os.path.isfile(html_rep):
        man.info(name, "no output/Analysis_Report.html to print -- step 9 (make_analysis_html.py) "
                       "has not run")
        return
    if os.path.isfile(pdf_rep) and os.path.getmtime(pdf_rep) >= os.path.getmtime(html_rep):
        status, note = "OK", "made from the current HTML"
    else:
        try:
            if print_report is None:
                import html_to_pdf
                print_report = html_to_pdf.print_report
            # an older PDF that could not be re-printed is renamed *.stale.pdf -> [SKIPPED]
            status, note = print_report(html_rep, pdf_rep)
        except Exception as e:                      # recorded, never swallowed
            status, note = "INFO", f"{type(e).__name__}: {e}"
    removed = None
    if os.path.isfile(stale):
        if status == "OK":
            try:
                os.remove(stale)
                removed = ("OK", "removed the stale copy superseded by the current PDF")
            except OSError as e:
                removed = ("SKIPPED", f"the stale copy could not be removed "
                                      f"({e.strerror or e}); it is left out of the zip")
        else:
            # First, so the 200-character cap on a [SKIPPED] reason never cuts it off.
            status = "SKIPPED"
            note = ("no current PDF; Analysis_Report.stale.pdf (an earlier report) is kept on disk "
                    "but left out of the zip and the file catalog, never send it -- " + note)
    {"OK": man.ok, "SKIPPED": man.skip}.get(status, man.info)(name, note)
    if removed:
        (man.ok if removed[0] == "OK" else man.skip)(
            "Report PDF (output/Analysis_Report.stale.pdf)", removed[1])


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
        # The report's own files stay beside it: Analysis_Report.pdf (and the .stale.pdf that
        # html_to_pdf.print_report leaves when it cannot re-print) are .pdf files, and moving
        # them into figures/ left stray copies there and made every finalize reprint.
        if (base in ("AI_Analysis_Report.md", "ANALYSIS_PROMPT.md", "OUTPUT_FILES.md")
                or base.startswith("Analysis_Report.")):
            continue
        if ext in ("csv", "tsv"):
            shutil.move(f, os.path.join(p["de_dir"], base))
        elif ext in ("png", "svg", "pdf", "jpg", "jpeg"):
            shutil.move(f, os.path.join(p["figures_dir"], base))

    # Publication Methods + the repository-deposit package (make_deposit.py). Every part is
    # recorded [OK] or [SKIPPED] -- with the reason -- in MANIFEST.txt at the session root, the
    # top of the zip: a part that could not be made is visible, never silently absent
    # (CLAUDE.md rule 4). Nothing here may stop finalize from writing the README and the zip.
    # input/raw_files.txt first: a hive_remote session is initialised without --raw, so it has
    # none, and the deposit package, README and AGENTS.md all read it. Written from the search's
    # own record of what it read (session_docs.raw_record).
    # The import and the call fail apart: a raw list that cannot be written must not also cost
    # the README and AGENTS.md (nor blame session_docs.py for it in MANIFEST.txt).
    try:
        import session_docs
    except Exception as e:                      # recorded below, not swallowed
        session_docs = None
        raw_list = ("SKIPPED", f"session_docs.py could not be loaded: {type(e).__name__}: {e}")
    else:
        try:
            raw_list = session_docs.ensure_raw_list(p)
        except Exception as e:                  # recorded below, not swallowed
            raw_list = ("SKIPPED", f"{type(e).__name__}: {e}")
    methods_md, deposit = None, None
    try:
        import make_deposit
        man = make_deposit.Manifest()
    except Exception as e:                      # recorded below, not swallowed
        make_deposit, man = None, None
        import_error = f"make_deposit.py could not be loaded: {type(e).__name__}: {e}"
    if man is not None:
        (man.ok if raw_list[0] == "OK" else man.skip)("input/raw_files.txt (where the raw data are)",
                                                       raw_list[1])
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
        report_pdf_step(p, man)
    registry, registry_note = _registry_lookup(p["session_dir"])
    docs = _write_docs(session_docs, p, man, registry, registry_note,
                       pending=("manifest",))        # MANIFEST.txt is written right below
    if man is not None:
        man.write(p["manifest_txt"], "Session export manifest")
    else:
        with open(p["manifest_txt"], "w", encoding="utf-8") as fh:
            fh.write("Session export manifest\n=======================\n"
                     f"[SKIPPED] {'Publication methods + deposit package':<50} -- "
                     f"{import_error}\n")
        _append_manifest(p["manifest_txt"], raw_list[0], "input/raw_files.txt (where the raw data are)",
                         raw_list[1])
        for level, part, note in docs["manifest"]:
            _append_manifest(p["manifest_txt"], level, part, note)

    de_files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(p["de_dir"], "DE_*.csv")))
    result = {"session_dir": p["session_dir"], "readme": p["readme"],
              "readme_html": docs["written"].get("readme_html"),
              "agents": docs["written"].get("agents"), "de_files": de_files,
              "methods": methods_md, "manifest": p["manifest_txt"],
              "deposit": deposit, "skipped": (man.n_skipped if man is not None else 1)}

    # re-analysis: write DIFFERENCES.md vs the parent
    parent = a.reanalysis_of
    marker = os.path.join(p["session_dir"], ".reanalysis_of")
    if not parent and os.path.exists(marker):
        with open(marker, encoding="utf-8", errors="replace") as fh:   # 2.7 and older: any encoding
            parent = fh.read().strip()
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
        # README.md, README.html and AGENTS.md go in after the run-log hook too: they name the
        # Core registry record, which record_run.py may only create then.
        late = {os.path.abspath(os.path.join(sdir, n)) for n in LATE_DOCS}
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
        # The predicted library: is_predicted_speclib() above.
        speclib_label = ("predicted spectral libraries (*.predicted.speclib anywhere in the "
                         "session; rebuilt from the FASTA + params, kept on disk where they are)")
        n_speclib = 0
        # Scratch (scratch_files.py, the one rule): every .cache folder -- make_podcast.py's TTS
        # chunks, ~60 MB per episode -- every *.part, and podcast.wav beside podcast.m4a. The
        # podcast itself (podcast.m4a, transcript, script, check.txt, podcast.json) goes in.
        import scratch_files
        n_scratch = 0

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
                files, k = scratch_files.prune(root, dirs, files)
                n_scratch += k
                for fn in files:
                    full = os.path.join(root, fn)
                    if fn.endswith(".quant"):
                        n_quant += 1
                        continue
                    if is_predicted_speclib(fn):
                        n_speclib += 1
                        continue
                    if os.path.islink(full) or os.path.abspath(full) == manifest_txt \
                            or os.path.abspath(full) in late:
                        continue                 # MANIFEST.txt and the docs go in last, below
                    z.write(full, os.path.join(base, os.path.relpath(full, sdir)))
        excluded[quant_label] = n_quant
        excluded[speclib_label] = n_speclib
        if n_scratch:
            excluded[scratch_files.LABEL] = n_scratch
        result["zip"] = archive
        result["zip_excluded"] = excluded

    # The run log and Slack go once everything they point at exists -- after the zip. Their
    # outcomes are the last two lines of MANIFEST.txt, so the zip's copy (added now) says them.
    try:
        with open(p["manifest_txt"], encoding="utf-8") as fh:
            pre_hook = fh.read()
    except OSError:
        pre_hook = None
    hooks = _finish_hooks(a, p["session_dir"], result.get("zip"))
    result["run_log"], result["slack"] = hooks["run_log"], hooks["slack"]
    rl = hooks["run_log"] or {}
    written = list(docs["written"].values())
    if rl.get("logged") is True and rl.get("path") and rl["path"] != registry:
        # the registry record exists now: name it in README / AGENTS (the registry's own copy of
        # the README, taken during the hook, is the one written above). A document this re-render
        # cannot write keeps the whole one written above (session_docs writes via a .part).
        again = _write_docs(session_docs, p, None, rl["path"], None, ok_lines=False)
        hooks["manifest"] = again["manifest"] + hooks["manifest"]
        written += list(again["written"].values())
    if a.zip:
        how = _zip_docs(result["zip"], p["session_dir"], written)
        result["zip_docs"] = how
        if how != "added":
            hooks["manifest"].insert(0, ("SKIPPED", "README / AGENTS.md in the zip", how))
    for level, part, note in hooks["manifest"]:
        try:
            _append_manifest(p["manifest_txt"], level, part, note)
        except Exception as e:
            sys.stderr.write(f"[session] could not record '{part}' in MANIFEST.txt: {e}\n")
    if a.zip:
        try:
            with open(p["manifest_txt"], encoding="utf-8") as fh:
                text = fh.read()
        except OSError:
            text = pre_hook
        result["zip_manifest"] = _zip_manifest(
            result["zip"], os.path.join(os.path.basename(p["session_dir"]), "MANIFEST.txt"),
            text if text is not None else "", pre_hook)
        if result["zip_manifest"] != "added":
            sys.stderr.write(f"[session] MANIFEST.txt in the zip: {result['zip_manifest']}\n")

    print(json.dumps(result, indent=2))


def do_docs(a):
    """README.md + README.html + AGENTS.md for a session, e.g. one finalized before they existed.
    Touches nothing else: no MANIFEST.txt line, no zip, no run log."""
    p = paths_for(a.dir)
    if not os.path.isdir(p["session_dir"]):
        sys.exit(f"session dir not found: {p['session_dir']}")
    import session_docs
    try:
        raw = session_docs.ensure_raw_list(p)
    except Exception as e:                      # reported below; the documents still get written
        raw = ("SKIPPED", f"{type(e).__name__}: {e}")
    registry, note = _registry_lookup(p["session_dir"], in_finalize=False)
    docs = _write_docs(session_docs, p, None, registry, note,
                       located_at=os.path.abspath(a.as_path) if a.as_path else None)
    print(json.dumps({"session_dir": p["session_dir"], "raw_files": {"level": raw[0],
                      "note": raw[1]}, "registry": registry, **docs["written"],
                      "parts": [{"level": l, "part": n, "note": t} for l, n, t in
                                docs["manifest"]]}, indent=2))


def params_file(p, resolved=True):
    """The search parameters file a session holds (paths_for dict -> path or None): with
    `resolved`, the cfg the search wrote as it ran (<search out>/params.resolved.cfg) first;
    then the one the workflow step staged (DIA-NN cfg / Sage json / FragPipe .workflow / Radiant
    config), newest convention first; then a params.* / *.cfg / sage_config*.json directly in
    input/ (older sessions). resolved=False is the file the search was GIVEN (make_deposit's
    Methods read)."""
    wf, inp = p["workflow_dir"], p["input_dir"]
    cands = [os.path.join(p["search_out"], "params.resolved.cfg")] if resolved else []
    cands += [os.path.join(wf, n) for n in ("params.cfg", "params.json")]
    for pat in ("*.cfg", "params*.json", "sage*.json", "*.workflow", "*.radiantConfig"):
        cands += sorted(glob.glob(os.path.join(wf, pat)))
    for pat in ("params.*", "*.cfg", "sage_config*.json"):
        cands += sorted(glob.glob(os.path.join(inp, pat)))
    return next((c for c in cands if os.path.isfile(c)
                 and not c.endswith((".rationale.json", "manifest.json"))), None)


def _block(prov):
    """de_provenance.json `block` (run_de.R / blocking.R) as column, effect, scope -- or None
    when the record predates it."""
    b = prov.get("block")
    if not isinstance(b, dict):
        return None
    if b.get("applied"):
        return (f"{b.get('column')} ({b.get('effect') or 'random'} effect, "
                f"scope {b.get('scope') or NOT_RECORDED})")
    if b.get("column"):
        return f"{b['column']} given, not applied" + (f": {b['note']}" if b.get("note") else "")
    return "none (samples modelled as independent)"


def _facts(session_dir):
    """Pull the comparable facts of a run from its session files. None = not recorded there."""
    p = paths_for(session_dir)
    man = _load(os.path.join(p["repro_dir"], "run_manifest.json")) or {}
    prov = _load(os.path.join(p["de_dir"], "de_provenance.json")) or {}
    wf = _load(os.path.join(p["workflow_dir"], "workflow.manifest.json")) or {}
    ptext = ""
    pf = params_file(p)
    if pf:
        try:
            with open(pf, encoding="utf-8", errors="replace") as fh:
                ptext = fh.read()
        except OSError:
            pass
    # the session's own sidecar first; the reproducibility bundle's copy of it for older ones
    fi = _load(p["fasta_meta"]) or (man.get("inputs") or {}).get("fasta_info") or {}
    try:
        from fetch_fasta import sidecar_state          # the one definition of the states
        fasta_state = sidecar_state(fi) if fi else None
    except Exception as e:                             # said in the table, not dropped
        fasta_state = f"could not be read ({type(e).__name__})"
    try:
        from session_docs import _sig_counts           # also reads an old repr() record
        sig = _sig_counts(prov)
    except Exception:
        sig = prov.get("significant_per_contrast") if isinstance(
            prov.get("significant_per_contrast"), dict) else {}
    cont = prov.get("contaminants") if isinstance(prov.get("contaminants"), dict) else None
    return {
        "engine": man.get("engine") or (wf.get("engine", {}) or {}).get("name"),
        "engine_version": (wf.get("engine", {}) or {}).get("version"),
        "de_method": prov.get("method") or (wf.get("de", {}) or {}).get("method"),
        "q_cutoff": prov.get("q_cutoff"), "logfc": prov.get("logfc"),
        "logfc_role": prov.get("logfc_role"), "adjp": prov.get("adjp"),
        "contrasts": prov.get("contrasts") or None,
        "design": prov.get("design"),
        "block": _block(prov),
        "contaminants": (cont.get("policy") + (f" ({cont['tag']} entries)" if cont.get("tag")
                                               else "")) if cont and cont.get("policy") else None,
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
                  ) if fi else None,
        "fasta_state": fasta_state,
        "skill_version": (man.get("skill") or {}).get("version"),
        "limpa": (prov.get("packages") or {}).get("limpa"),
        "raw": read_raw_list(session_dir),
        "significant_per_contrast": sig,
        "params_text": ptext,
        "conditions": p["conditions"],
    }


# (label, _facts key) -- every setting DIFFERENCES.md compares, in the order it lists them
DIFF_FIELDS = [("Search engine", "engine"), ("Engine version", "engine_version"),
               ("DE method", "de_method"), ("ID FDR (q)", "q_cutoff"),
               ("logFC value", "logfc"), ("logFC role (de_provenance.json)", "logfc_role"),
               ("adj.P threshold", "adjp"), ("Contrasts", "contrasts"),
               ("Design (with covariates)", "design"),
               ("Blocking (column, effect, scope)", "block"),
               ("Contaminant policy (DE)", "contaminants"), ("FASTA", "fasta"),
               ("FASTA sidecar state (fetch_fasta.sidecar_state)", "fasta_state"),
               ("Skill version", "skill_version"), ("limpa version", "limpa")]


def write_differences(prior_dir, new_dir):
    prior_dir = os.path.abspath(os.path.expanduser(prior_dir))
    old, new = _facts(prior_dir), _facts(new_dir)
    shown = lambda v: NOT_RECORDED if v is None else (
        ", ".join(map(str, v)) if isinstance(v, list) else str(v).replace("|", "\\|"))
    same_raw = bool(old["raw"]) and sorted(old["raw"]) == sorted(new["raw"])
    L = [f"# What changed in this re-analysis", "",
         f"Re-analysis of `{prior_dir}`.", "",
         (f"Same raw files ({len(new['raw'])}, input/raw_files.txt). " if same_raw else "")
         + "The table below is exactly what differs among the settings both sessions record. "
           "Unchanged settings are omitted.",
         "", "| Aspect | Original | This re-analysis |", "|---|---|---|"]
    if not same_raw:
        common = len(set(old["raw"]) & set(new["raw"]))
        L.append(f"| Raw files (input/raw_files.txt) | {len(old['raw']) or NOT_RECORDED} | "
                 f"{len(new['raw']) or NOT_RECORDED}"
                 + (f" ({common} in common)" if old["raw"] and new["raw"] else "") + " |")
    n_changes, neither = 0, []
    for label, key in DIFF_FIELDS:
        a_v, b_v = old.get(key), new.get(key)
        if a_v is None and b_v is None:
            neither.append(label)
        elif a_v != b_v:
            n_changes += 1
            L.append(f"| {label} | {shown(a_v)} | {shown(b_v)} |")
    if n_changes == 0:
        L.append("| (settings) | — | none of the settings compared here differ"
                 + ("" if same_raw else "; check the raw-file row") + " |")
    if neither:
        L += ["", f"_Not recorded in either session, so not compared: {', '.join(neither)}._"]

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
    with open(out, "w", encoding="utf-8") as fh:
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
    d = sub.add_parser("docs", help="(re)write README.md, README.html and AGENTS.md -- and "
                                    "input/raw_files.txt if missing -- without finalizing")
    d.add_argument("--dir", required=True, help="the session directory")
    d.add_argument("--as", dest="as_path", default="",
                   help="where the session really is, when --dir is a copy of it (the README's "
                        "'Where this lives on HIVE' table is then about that location)")
    d.set_defaults(func=do_docs)
    a = ap.parse_args()
    a.func(a)


if __name__ == "__main__":
    main()
