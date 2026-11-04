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
  #   a re-analysis BESIDE the earlier session (which is not written to): <date>_<name>_v2
  python3 session.py init --name "HeLa QC DIA" --raw ... --reanalysis-of <prior> --beside
  #   Core staff: what the analysis is for (SKILL.md step 0a), kept in session.json's
  #   "analysis_for" -- the CoreOmics submission(s), the first the main one, or internal work
  #   with its title on the Core Project Board (references/board.md section 2)
  python3 session.py init --name PROT_0807 --raw ... --submission PROT_0807 [--submission 9002]
  python3 session.py init --name "Astral HeLa gradient test" --raw ... --internal "Astral HeLa gradient test"
  #   -> prints JSON with every canonical path + "placement"; route later steps into them

  # at the end — ensure the publication Methods, write the deposit package
  # (output/DATA_SUBMISSION, see make_deposit.py), README + MANIFEST.txt, optionally zip,
  # then log the run (record_run.py) and post "analysis complete" to the Core's Slack channel
  # (notify_slack.py)
  python3 session.py finalize --dir <session_dir> [--zip] [--no-deposit] [--no-notify]
  #   it refuses (exit 2, nothing written) without the skill's own report, or for a Claude Code
  #   analysis with no saved conversation; it BUILDS a missing or stale reproducibility bundle
  #   (refusing only when a record the bundle needs is absent); --override-finish "<reason>" /
  #   --no-transcript "<reason>" go ahead and are recorded. Its last act is the finalize stamp,
  #   output/FINALIZED.json (exit 2 without AGENTS.md, README.md or the Methods)

  # before telling the user an analysis is finished: the finish checklist (exit 0 complete, 3 not)
  python3 session.py status --session <session_dir> [--json]     # --dir works too, as finalize

  # README.md + README.html + AGENTS.md alone (e.g. for a session finalized before they existed)
  python3 session.py docs --dir <session_dir>

Raw MS files are NOT copied (they're huge and live elsewhere) — their paths are
recorded in input/raw_files.txt instead (finalize writes it from the search's record when a
session was initialised without --raw).
"""
import sys, os, json, re, glob, shutil, argparse, datetime

NOT_RECORDED = "not recorded"      # a fact the session's records do not give -- never a guess
# session.json: what a Core staff analysis is for (init --submission / --internal)
ANALYSIS_FOR = "analysis_for"

SUBDIRS = ["input", "output", "output/search", "output/tables", "output/figures",
           "output/reproducibility", "scripts", "logs"]


def is_session_dir(path):
    """THE test of "this folder is an analysis session", for a writer that must put session files
    (RECOVERY.md, .recovery.json) only in one: it has session.json (init --submission/--internal,
    submission_report.py attach), or the layout init always makes -- input/ and output/ (SUBDIRS).
    A chain whose --out sat outside a session once wrote its recovery notes into a home folder."""
    p = paths_for(path)
    return (os.path.isfile(p["session_json"])
            or (os.path.isdir(p["input_dir"]) and os.path.isdir(os.path.join(p["session_dir"],
                                                                            "output"))))


def recovery_dir(out):
    """THE folder a search's recovery notes (checkpoint.py: RECOVERY.md, .recovery.json) live in,
    for every writer and reader (diann_parallel.py's submit.sh, node_fault.py's retry): the
    session two levels above the search's --out (<session>/output/search) when that folder IS a
    session (is_session_dir), else beside --out -- never a home folder that happens to sit there."""
    out = os.path.abspath(out)
    sess = os.path.dirname(os.path.dirname(out))
    return sess if is_session_dir(sess) else out


def session_of_search(out):
    """THE session a search folder belongs to (the Slack posts, the run log), the other side of
    find_search_provenance: <S> for a search in <S>/output/<folder> or <S>/output/<folder>/<sub>
    (the two depths find_search_provenance reads) when <S> is a session (is_session_dir) --
    output/search_dda/ and output/search_sage/ as much as output/search/. <S>/output/search/,
    the default folder, names <S> as it always has. None for a folder anywhere else."""
    out = os.path.abspath(out)
    up = os.path.dirname(out)
    if os.path.basename(out) == "search" and os.path.basename(up) == "output":
        return os.path.dirname(up)
    for o in (up, os.path.dirname(up)):
        if os.path.basename(o) == "output":
            sess = os.path.dirname(o)
            return sess if is_session_dir(sess) else None
    return None


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
        # the DEFAULT search folder (init makes it; a new search goes here). A search may be
        # elsewhere (output/search_dda/, output/search_sage/): THE search record is
        # search_record(p) / find_search_provenance, never these two paths read alone (T44)
        "search_out": os.path.join(d, "output", "search"),
        "search_prov": os.path.join(d, "output", "search", "search_provenance.json"),
        "de_dir": os.path.join(d, "output", "tables"),
        "figures_dir": os.path.join(d, "output", "figures"),
        # the peak check (skyline_peaks.py, step 8f): pictures, tables, the Skyline document
        "skyline_dir": os.path.join(d, "output", "skyline"),
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
        # the staff decisions behind the file list (core_submission.py locate: run choices and
        # the re-injection policy), carried in by submission_report.py attach
        "locate_decisions": os.path.join(d, "input", "locate_decisions.json"),
    }


def analysis_for_block(submissions, internal):
    """session.json's "analysis_for" from init's --submission / --internal: {"submissions":
    ["PROT_0807", "PROT_9002"]} (as staff named them, the first the main one) or {"internal":
    "<its title on the Core Project Board>"}. None when neither was given: only Core staff are
    asked (SKILL.md step 0a). ValueError for anything that is not one or the other."""
    if submissions and internal is not None:
        raise ValueError("--submission and --internal: an analysis is for CoreOmics submissions "
                         "or internal work, not both")
    if internal is not None:
        title = " ".join(internal.split())
        if not title:
            raise ValueError("--internal: give the title the analysis has on the Core Project Board")
        return {"internal": title}
    if not submissions:
        return None
    try:
        from core_submission import normalize_submission  # the one reading of a submission number
    except ImportError as e:
        raise ValueError(f"--submission needs core_submission.py beside session.py ({e}): "
                         "copy the whole scripts folder")
    out = []
    for s in submissions:
        kind, val = normalize_submission(s)               # ValueError says what it accepts
        if kind != "internal_id":
            raise ValueError(f"--submission {s}: give the submission's PROT number (PROT_0807), "
                             "not its CoreOmics id")
        if val not in out:
            out.append(val)
    return {"submissions": out}


def read_analysis_for(session_dir):
    """(the session's "analysis_for" or None, why it could not be read or None). The second is
    set when session.json EXISTS but cannot be read -- damaged, not a JSON object, or its
    "analysis_for" not an object -- so a resume (checkpoint.py status) says so instead of
    silently losing which submission or board thread the analysis is for (2.12 review). No
    session.json, or one without the key, is (None, None)."""
    p = paths_for(session_dir)["session_json"]
    try:
        with open(p, encoding="utf-8") as fh:
            data = json.load(fh)
    except FileNotFoundError:
        return None, None
    except (OSError, ValueError) as e:              # UnicodeDecodeError is a ValueError
        return None, f"{p} could not be read ({type(e).__name__}: {e})"
    if not isinstance(data, dict):
        return None, f"{p} is not a JSON object"
    v = data.get(ANALYSIS_FOR)
    if v is None or isinstance(v, dict):
        return v, None
    return None, f"{p} has an \"{ANALYSIS_FOR}\" that is not an object"


def analysis_for(session_dir):
    """The session's "analysis_for" (see analysis_for_block), or None when it has none or
    session.json cannot be read (read_analysis_for says why)."""
    return read_analysis_for(session_dir)[0]


def _record_analysis_for(p, block):
    """Add "analysis_for" to session.json, keeping what is there (submission_report.py attach
    writes "coreomics" into the same file). Refuses a session.json that is not a JSON object."""
    sj = {}
    if os.path.isfile(p["session_json"]):
        sj = _load(p["session_json"])
        if not isinstance(sj, dict):
            sys.exit(f"{p['session_json']} is not readable JSON: look at it before init writes "
                     "to it again")
    sj[ANALYSIS_FOR] = block
    tmp = p["session_json"] + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        json.dump(sj, fh, indent=2)
        fh.write("\n")
    os.replace(tmp, p["session_json"])


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
    """Scan existing sessions for ones covering the same raw files (same dataset).

    A match is a fact to mention, not a decision: a new session is FRESH by default. It
    reads nothing from the earlier one and does not nest under it unless the user asks for a
    re-analysis (init --reanalysis-of). It used to print "re-analysis of <session>" as its
    verdict, and a fresh run was steered into reusing old results (Core staff, skill 2.9)."""
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
    if hits:
        suggestion = ("a fresh analysis (the default). These raw files were analysed before, in "
                      + hits[0]["session"] + ": mention that in one line and read nothing from "
                      "it. Re-analyse it (init --reanalysis-of, plus --beside to put the new run "
                      "beside it as _v<N> rather than inside it) only if the user asks to.")
    else:
        suggestion = "a fresh analysis: no earlier session covers these raw files"
    print(json.dumps({"query_raw_count": len(mine), "matches": hits, "default": "fresh",
                      "suggestion": suggestion}, indent=2))


# A re-analysis's link to the analysis it re-does: the marker (one line, the prior session's
# path; what 2.9 and older wrote) and the record beside it (the version and the prior's own
# provenance files). reanalysis_of() is the one reader.
REANALYSIS_MARKER = ".reanalysis_of"
REANALYSIS_JSON = "reanalysis.json"
_VERSION_RE = re.compile(r"_v(\d+)$")


def reanalysis_of(session_dir):
    """{"session": prior dir, "version": N or None, "placement": ..., ...} for a re-analysis
    session, else None. From reanalysis.json, else the marker alone (sessions made before it)."""
    rec = _load(os.path.join(session_dir, REANALYSIS_JSON))
    if isinstance(rec, dict) and rec.get("reanalysis_of"):
        return dict(rec, session=rec["reanalysis_of"])
    marker = os.path.join(session_dir, REANALYSIS_MARKER)
    if os.path.exists(marker):
        with open(marker, encoding="utf-8", errors="replace") as fh:   # 2.7 and older: any encoding
            prior = fh.read().strip()
        if prior:
            return {"session": prior, "reanalysis_of": prior, "version": None,
                    "placement": None, "source": REANALYSIS_MARKER}
    return None


def version_of(session_dir):
    """A session's version: its re-analysis record's, else a trailing _v<N> in its name, else 1."""
    rec = reanalysis_of(session_dir) or {}
    if isinstance(rec.get("version"), int):
        return rec["version"]
    m = _VERSION_RE.search(os.path.basename(os.path.normpath(session_dir)))
    return int(m.group(1)) if m else 1


def beside_dir(prior, date, slug):
    """The folder for a re-analysis BESIDE `prior`: <prior's parent>/<date>_<slug>_v<N>, N one
    above the prior's version and above every <any date>_<slug>_v<M> already there -- never an
    existing folder (two re-analyses of one session do not share a v2). Returns (dir, N)."""
    parent = os.path.dirname(os.path.normpath(prior))
    stem = _VERSION_RE.sub("", slug)
    taken = [version_of(prior)]
    pat = re.compile(r"^\d{4}-\d{2}-\d{2}_" + re.escape(stem) + r"_v(\d+)$")
    try:
        names = os.listdir(parent)
    except OSError:
        names = []
    taken += [int(m.group(1)) for m in map(pat.match, names) if m]
    n = max(taken) + 1
    while os.path.exists(os.path.join(parent, f"{date}_{stem}_v{n}")):
        n += 1
    return os.path.join(parent, f"{date}_{stem}_v{n}"), n


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
    #   --reanalysis-of <prior> --beside -> beside the original, <date>_<name>_v<N>; nothing is
    #                                       written into the original (a delivered folder)
    #   --reanalysis-of <prior>  -> nested under the original
    #   --base <path>            -> a central location the user chose (e.g. ~/Documents/DataAnalysis)
    #   (neither, with --raw)    -> DEFAULT: in the folder with the raw data being analyzed
    version = None
    if getattr(a, "beside", False) and not a.reanalysis_of:
        sys.exit("--beside places a re-analysis beside the session it re-does: give that session "
                 "with --reanalysis-of <prior session dir>")
    try:
        wanted = analysis_for_block(getattr(a, "submission", None), getattr(a, "internal", None))
    except ValueError as e:
        sys.exit(str(e))
    if a.reanalysis_of:
        prior = os.path.abspath(os.path.expanduser(a.reanalysis_of))
        if not os.path.isdir(prior):
            sys.exit(f"--reanalysis-of: prior session not found: {prior}")
        if getattr(a, "beside", False):
            # A Core re-analysis (2026-09-25): the prior was a delivered, read-only-by-policy folder,
            # and nesting wrote into it; the agent had to patch _raw_dir from a wrapper.
            session_dir, version = beside_dir(prior, date, slug)
            placement = "beside-prior"
        else:
            session_dir = os.path.join(prior, "reanalysis", f"{date}_{slug}")
            placement = "reanalysis"
            version = max([version_of(prior)] + [
                version_of(d) for d in glob.glob(os.path.join(prior, "reanalysis", "*"))
                if os.path.isdir(d) and os.path.abspath(d) != os.path.abspath(session_dir)]) + 1
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
    if wanted:
        _record_analysis_for(p, wanted)

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

    # record the parent when this is a re-analysis: the marker, and the record that links this
    # session's provenance to the prior's (DIFFERENCES.md, the README, AGENTS.md and the run log
    # read it through reanalysis_of())
    parent = os.path.abspath(os.path.expanduser(a.reanalysis_of)) if a.reanalysis_of else None
    if parent:
        with open(os.path.join(session_dir, REANALYSIS_MARKER), "w", encoding="utf-8") as fh:
            fh.write(parent + "\n")
        pp = paths_for(parent)
        # the prior's search record wherever its search is (search_record, T44)
        prior_search = search_folder(pp)
        pp = dict(pp, search_prov=os.path.join(prior_search, SEARCH_RECORD) if prior_search
                  else "")
        link = {"reanalysis_of": parent, "version": version, "prior_version": version_of(parent),
                "placement": placement, "created": datetime.datetime.now().isoformat(
                    timespec="seconds"),
                # the prior's own records, as they are now: what DIFFERENCES.md compares with
                "prior_records": {k: pp[k] for k in ("search_prov", "readme", "raw_list")
                                  if os.path.isfile(pp[k])}}
        de_prov = os.path.join(pp["de_dir"], "de_provenance.json")
        if os.path.isfile(de_prov):
            link["prior_records"]["de_prov"] = de_prov
        with open(os.path.join(session_dir, REANALYSIS_JSON), "w", encoding="utf-8") as fh:
            json.dump(link, fh, indent=2)

    # starter README (finalize fills in results)
    with open(p["readme"], "w", encoding="utf-8") as fh:
        fh.write(f"# {a.name}\n\n- Date: {date}\n- Status: in progress\n")
        if parent:
            fh.write(f"- **Re-analysis of:** `{parent}` (this is version {version}) — see "
                     "`DIFFERENCES.md` (written at finalize) for exactly what changed.\n")
        import report_files                      # the skill's own: finalize writes over it
        fh.write(f"\n{report_files.STARTER_README_LAYOUT}\n")
        fh.write("\n" + report_files.SKILL_DOC_MARKER + "\n")
    # ... while it is as written here: one edited by hand is kept (report_files.is_skill_doc)
    report_files.record_skill_doc(p["readme"], "session.py init")

    # this Claude Code conversation runs this analysis: the transcript hook saves it here
    # (save_transcript.py; a no-op outside Claude Code, and never a reason to fail init). Whether
    # it was recorded is in the JSON, not only on stderr: a refusal (subagents at work and the
    # call not traced to one conversation) must reach the agent, which tells the user.
    try:
        import save_transcript
        transcript = save_transcript.remember_result(*save_transcript.remember(
            session_dir, from_now=getattr(a, "transcript_from_now", False),
            call_id=getattr(a, "call_id", None)))
    except Exception as e:
        sys.stderr.write(f"[session] conversation not recorded for the transcript hook: {e}\n")
        transcript = {"recorded": False, "reason": f"{type(e).__name__}: {e}"}

    print(json.dumps({"created": session_dir, "date": date, "name": a.name,
                      "placement": placement, "reanalysis_of": parent, "version": version,
                      ANALYSIS_FOR: wanted, "paths": p, "transcript": transcript}, indent=2))


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


def _finish_hooks(a, session_dir, zip_path, not_finalized=None, finish_override=None):
    """After the zip: log the run in the Core's run log (record_run.py analysis-done, when this
    install has it), then post "analysis complete" to the Core's Slack channel -- in that order,
    so the post can say whether the run was logged. Never fatal (notify_slack.py rule 1).
    Returns {"run_log", "slack": {"sent", "level", "detail"}, "manifest": [(level, part, note)]};
    the words are notify_slack's (run_log_manifest / slack_manifest), which never name a path or
    a group -- the zip goes to collaborators. `finish_override`: the --override-finish reason
    when this finalize went ahead without a finish item; the post says NOT FINALIZED with it."""
    try:
        import notify_slack
    except Exception as e:                      # recorded, never swallowed
        why = f"notify_slack.py could not be loaded: {type(e).__name__}"
        sys.stderr.write(f"[session] {why}\n")
        return {"run_log": None, "slack": {"sent": False, "level": "SKIPPED", "detail": why},
                "manifest": [("SKIPPED", "Core run log + Slack notification", why)]}
    if not_finalized:
        # a finalize that refused to stamp is not logged as a finished analysis: record_run's
        # analysis-done would call it finalized and complete (C2 errors_ux M2), as Slack is not
        # told either (below). A finalize that stamps, later, logs it.
        run_log = None
        rl = ("SKIPPED", "Core run log", "not logged as a finished analysis: the session is not "
                                         "finalized -- finalize again once it is")
    else:
        run_log = notify_slack.record_run("analysis-done", session=session_dir)
        rl = notify_slack.run_log_manifest(run_log)
    if a.no_notify:
        sent, detail = False, "--no-notify was given"
    elif not_finalized:
        # "analysis complete" is for a finalized session only (report_files.FINISH_NEEDS
        # "announce"): this finalize refused to stamp, so the channel is not told
        sent, detail = False, f"not posted: the session is not finalized ({not_finalized})"
    elif (run_log or {}).get("analysis_logged") == "unchanged":
        # finalize run again on the same analysis (to add a podcast, say): the run log already
        # has it, unchanged, so the Core channel is not told a second time
        sent, detail = False, ("not posted again: this analysis was already logged and posted, "
                               "and nothing in it changed")
    else:
        sent, detail = notify_slack.analysis_done(session_dir, zip_path, run_log=run_log,
                                                  finish_override=finish_override)
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


def skill_doc_paths(session_dir, kept=False):
    """Where the skill's README.md, README.html and AGENTS.md are written (report_files.
    skill_doc_targets: README.skill.md ... beside one written by hand) -- and, with `kept`, the
    hand-written ones too."""
    import report_files
    t = report_files.skill_doc_targets(session_dir)
    return [t[n] for n in LATE_DOCS] + ([os.path.join(session_dir, n) for n in t["kept"]]
                                        if kept else [])

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

            def info(self, name, why):          # a document kept as written by hand
                if ok_lines:
                    out["manifest"].append(("INFO", name, why))
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
        import make_podcast
        listed = skill_doc_paths(p["session_dir"], kept=True) + [p["manifest_txt"]]
        listed.append(p["raw_list"])
        listed.append(os.path.join(p["output_dir"], make_podcast.SHARE_NAME))   # when there is one
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
            for full in skill_doc_paths(session_dir):
                n = os.path.basename(full)
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


def has_de(p):
    """The session has a differential expression: run_de.R's record (report_files.DE_RECORD)."""
    import report_files
    return os.path.isfile(os.path.join(p["de_dir"], report_files.DE_RECORD))


def peak_check_step(p, man):
    """The peak check (skyline_peaks.py, step 8f): its outcome as one MANIFEST.txt line -- [OK]
    with its counts, [SKIPPED] with the record's reason (the user turned it off, the container was
    not reachable, a job failed) or because no record exists for a session with a DE. A session
    with neither has no DE to check: [SKIPPED] with NO_DE_REASON. Never fatal."""
    name = PEAK_CHECK_PART
    try:
        import skyline_peaks
        got = skyline_peaks.manifest_line(p["output_dir"], has_de(p))
    except Exception as e:                      # recorded, never swallowed
        man.skip(name, f"{type(e).__name__}: {e}")
        return
    if got:
        {"OK": man.ok, "SKIPPED": man.skip}.get(got[0], man.info)(got[1], got[2])
    else:
        man.skip(name, NO_DE_REASON)


# ------------------------------------------------------------------ a session of any type --
# finalize is every analysis's last step, not the DE's (2.14.2, Amendment 2 to T13: the beer
# re-analysis and three earlier HT submissions took it for the DE path's, and never ran it). A
# session with no differential expression -- a search-only analysis, an identification list, a
# re-analysis or a benchmark with its own question -- gets AGENTS.md, the README, MANIFEST.txt
# and its (search-only) Methods like any other, and MANIFEST.txt names each part only a DE makes
# as [SKIPPED], with NO_DE_REASON: never a silent gap (DE-LIMP rule 4). The bundle is not one of
# them: finalize builds it without --de-dir when step 10 did not (make_deposit.ensure_bundle).
NO_DE_REASON = ("no differential expression in this session (no output/tables/de_provenance.json"
                "): a search-only analysis. If one was wanted, run step 8, then finalize again")
DE_ONLY_PARTS = ("Differential expression tables (output/tables)",
                 "The DE as R (output/tables/reproducibility_log.R)",
                 "Publication methods: the differential-expression section")
PEAK_CHECK_PART = "Peak check (Skyline, output/skyline)"      # skyline_peaks.manifest_line's name
SEARCH_RECORD_PART = "Search record (output/<search>/search_provenance.json)"


def de_only_step(p, man):
    """A session without a DE: one [SKIPPED] line per DE_ONLY_PARTS part. Nothing with one."""
    if not has_de(p):
        for part in DE_ONLY_PARTS:
            man.skip(part, NO_DE_REASON)


def search_record_step(p, man):
    """The search record the README and AGENTS.md describe (search_choice), when it is not
    plain: none at all, or records and none chosen -- one [SKIPPED] line. Nothing for one."""
    c = search_choice(p)
    if c["how"] == "none":
        man.skip(SEARCH_RECORD_PART, "no search record in this session (no output/<folder>/"
                                     "search_provenance.json): the README, AGENTS.md and the "
                                     "Methods describe the search from the workflow manifest, "
                                     "tagged")
    elif search_unchosen(c):
        # the key phrase first: a MANIFEST reason is cut at 200 characters
        man.skip(SEARCH_RECORD_PART, "the README, AGENTS.md and the Methods describe none of "
                                     f"them -- {c['why']}")


def custom_step(p, gate):
    """Work done outside the session (finalize --report/--methods/--code-dir/--extra-dir, planned
    and checked by finalize_gate): copied into output/custom/, its checksums into the bundle's
    run_manifest.json; with none given, an earlier finalize's items are checked again and put back
    into a bundle provenance.py wrote anew. -> [(level, part, note)] for MANIFEST.txt. Never
    fatal: a copy that fails is a [SKIPPED] line, it leaves no partial copy and the earlier
    copies and record as they were (custom_work.apply), and those earlier items are still listed
    (C2 errors_ux M4)."""
    import custom_work
    failed = []
    try:
        items = custom_work.apply(p["session_dir"], p["output_dir"], gate.get("custom") or [],
                                  gate.get("override_reason"))
    except Exception as e:                      # recorded, never swallowed
        failed = [("SKIPPED", custom_work.MANIFEST_PART,
                   f"not folded in, no partial copy, earlier copies kept: {type(e).__name__}: {e}")]
        items = custom_work.load(p["output_dir"])
    if not items:
        return failed
    lines = failed + custom_work.manifest_lines(p["session_dir"], items)
    try:
        level, note = custom_work.bundle_merge(p["output_dir"], items)
    except Exception as e:                      # recorded, never swallowed
        level, note = "SKIPPED", f"{type(e).__name__}: {e}"
    return lines + [(level, f"{custom_work.MANIFEST_PART} in the bundle", note)]


def podcast_share_step(p, man):
    """With a podcast (step 9e): output/Analysis_Report_with_audio.html, the report with the audio
    and transcript built in -- the one file a collaborator can send on (make_podcast.ensure_share:
    left alone when current, else built; a copy that cannot be brought up to date is renamed
    *.stale.html, which scratch_files keeps out of the zip). Never fatal; nothing without a
    podcast."""
    name = "Shareable report with audio"
    try:
        import make_podcast
        name = f"{name} (output/{make_podcast.SHARE_NAME})"
        got = make_podcast.ensure_share(p["output_dir"])
    except Exception as e:                      # recorded, never swallowed
        man.skip(name, f"{type(e).__name__}: {e}")
        return
    if got:
        {"OK": man.ok, "SKIPPED": man.skip}.get(got[0], man.info)(name, got[1])


def conversation_step(p, man, call_id=None, no_transcript=None):
    """The analysis conversation (save_transcript.py): this Claude Code conversation saved,
    redacted, into logs/conversation/ beside any saved earlier (a resume, a new chat, a HIVE
    session saved from the laptop). Core-internal: deliver never ships it and the zip leaves it
    out. [OK] when a conversation is saved; [SKIPPED] when this IS Claude Code and nothing was
    saved (a real failure, whatever else is there); [INFO] when it is not Claude Code but the
    decisions log (logs/decisions.md) is there -- another agent, or finalize on HIVE in
    hive_remote, where save_transcript.py --hive saves it from the laptop; [SKIPPED] when there
    is neither: a run with no record of how it was done. With `no_transcript` (finalize
    --no-transcript REASON, which waives finalize_gate's transcript check) and nothing saved:
    [SKIPPED] with that reason. Never fatal."""
    name = "Analysis conversation (logs/conversation)"
    try:
        import save_transcript
        got = save_transcript.save(p["session_dir"], call_id=call_id)
        n = save_transcript.saved_count(p["session_dir"])
    except Exception as e:                      # recorded, never swallowed
        man.skip(name, f"{type(e).__name__}: {e}")
        return
    decisions = os.path.isfile(os.path.join(p["logs_dir"], "decisions.md"))
    if n:
        man.ok(name, f"{n} conversation(s), redacted -- Core-internal: never delivered, not in "
                     "the zip" + ("" if got.get("saved") else f" (this one: {got.get('skipped')})"))
    elif no_transcript:                              # waived, with the reason (I35)
        man.skip(name, f"--no-transcript: {no_transcript}")
    elif os.environ.get("CLAUDE_CODE_SESSION_ID"):   # Claude Code, and nothing saved: a failure
        man.skip(name, got.get("skipped") or "nothing saved")
    elif decisions:
        man.info(name, "no Claude Code transcript (" + (got.get("skipped") or "none") + "); the "
                       "decisions log (logs/decisions.md) is the record")
    else:
        man.skip(name, "no conversation and no decisions log (log_decision.py): nothing records "
                       "how this analysis was done -- " + (got.get("skipped") or "nothing saved"))


# ------------------------------------------------------------------------- the finish gates --
# 2.14.2 (I35; Brett 2026-10-31: "an agents.md file and all the reproducibility and
# methods are always done"). What each tool needs is report_files.FINISH_NEEDS (one definition);
# finalize also owes the conversation of a Claude Code analysis, and writes the finalize stamp
# (output/FINALIZED.json) as its last act. `session.py status` reads all of it back as a list.
FINISH_GATE_ENV = "SKILL_FINISH_GATE_DONE"      # finalize's gate already spoke (main, before srun)
FINALIZE_TOOL = "session.py finalize"
# The stamp's outputs and which it requires: report_files.STAMP_OUTPUTS (one definition).
_DECISIONS_CONVERSATION = re.compile(r"(?m)^- \*\*Conversation:\*\* `([^`]+)`")


def _reason(text):
    """A reason given on the command line, on one line; None when none (or only blanks)."""
    t = " ".join(str(text or "").split())
    return t or None


def claude_code_record(session_dir):
    """(the Claude Code conversation this analysis names, where that is recorded), or
    (None, None). Only records that exist already are read: CLAUDE_CODE_SESSION_ID here; the
    transcript hook's record that `session.py init` (or checkpoint --resume) wrote for this
    folder (save_transcript's map, on the machine that ran it); a decisions-log entry made
    inside Claude Code (log_decision.py's "- **Conversation:**" line). A hive_remote session
    finalized on HIVE has none of them unless the laptop's conversation logged a decision."""
    sid = os.environ.get("CLAUDE_CODE_SESSION_ID")
    if sid:
        return sid, "CLAUDE_CODE_SESSION_ID is set here"
    me = os.path.abspath(session_dir)
    m = None
    try:
        import save_transcript
    except ImportError as e:                     # said: the hook's record goes unchecked
        sys.stderr.write(f"[session] WARNING: save_transcript.py could not be loaded ({e}): the "
                         "transcript hook's record of this folder is not checked\n")
        save_transcript = None
    if save_transcript is not None:
        path = save_transcript.map_path()
        try:
            with open(path, encoding="utf-8") as fh:
                m = json.load(fh)
        except FileNotFoundError:
            m = None                             # no map here: init recorded nothing on this machine
        except (OSError, ValueError) as e:       # said: an unreadable map is not "no record"
            sys.stderr.write(f"[session] WARNING: the transcript hook's record {path} cannot be "
                             f"read ({type(e).__name__}): whether init recorded a Claude Code "
                             "conversation for this folder is not checked\n")
    if isinstance(m, dict):
        rec_of = getattr(save_transcript, "_record", None)
        for k in sorted(m):
            sessions = ((rec_of(m, k) or {}).get("sessions") or {}).values() if rec_of else \
                save_transcript.lookup(k)
            if any(isinstance(x, dict) and not x.get("hive") and x.get("session_dir")
                   and os.path.abspath(x["session_dir"]) == me for x in sessions):
                return k, ("the transcript hook's record of this folder (session.py init / "
                           "checkpoint.py --resume)")
    try:
        with open(os.path.join(session_dir, "logs", "decisions.md"), encoding="utf-8",
                  errors="replace") as fh:
            found = _DECISIONS_CONVERSATION.findall(fh.read())
    except OSError:
        found = []                               # no decisions log: nothing recorded there
    if found:
        return found[-1], "logs/decisions.md names a Claude Code conversation"
    return None, None


def transcript_gap(p, call_id=None, try_save=True):
    """The conversation finalize owes (I35), as one sentence with the commands that close it, or
    None: nothing owed (not a Claude Code analysis, or a conversation is saved). With try_save,
    an unsaved conversation is saved first (save_transcript.save, as conversation_step does)."""
    sid, where = claude_code_record(p["session_dir"])
    if not sid:
        return None
    import save_transcript
    n = save_transcript.saved_count(p["session_dir"])
    said = None
    if not n and try_save:
        said = (save_transcript.save(p["session_dir"], call_id=call_id) or {}).get("skipped")
        n = save_transcript.saved_count(p["session_dir"])
    if n:
        return None
    import shlex
    s = shlex.quote(p["session_dir"])
    return (f"the analysis conversation (logs/conversation): this is a Claude Code analysis "
            f"({where}) and no conversation is saved"
            + (f" (the save here said: {said})" if said else "")
            + f" -- save it with: python3 scripts/save_transcript.py {s} [--call-id <token>] "
              f"(from the laptop, for a HIVE session: python3 scripts/save_transcript.py --hive "
              f"{s}), or pass --no-transcript \"<why>\"")


def finalize_gate(a, p, speak=True):
    """finalize's finish gate, before anything is written (and, on a login node, before the
    compute-node re-exec, so a refusal costs no job and is said once): the skill's report
    (report_files.FINISH_NEEDS["finalize"]); the reproducibility bundle only when finalize cannot
    build it (a record it needs is absent: _bundle_plan) -- a missing or stale one it builds; for
    a Claude Code analysis, a saved conversation; a --report the skill did not write is a gap
    too, and work that cannot be folded in (custom_work.plan) is refused outright. Refuses (exit
    2, nothing written -- the earlier stamp is removed first -- no zip) listing every gap with its
    command, unless --override-finish REASON. A mistake in the command (an empty reason, a --code-dir
    that is not there) is refused first, and keeps the earlier stamp.
    -> {"gaps", "warnings", "override_reason", "override", "no_transcript", "custom",
    "build_bundle"}."""
    import report_files as rf
    # Mistakes in the command itself first (an empty reason, a --code-dir that is not there):
    # refused with nothing changed, the earlier stamp included -- it is not a finalize that ran
    # (L5b review L4).
    reason = _reason(getattr(a, "override_finish", None))
    no_tr = _reason(getattr(a, "no_transcript", None))
    for flag, given, val in (("--override-finish", getattr(a, "override_finish", None), reason),
                             ("--no-transcript", getattr(a, "no_transcript", None), no_tr)):
        if given is not None and not val:
            sys.stderr.write(f"[session] {flag} needs a reason: {flag} \"<why>\" (nothing was "
                             "changed)\n")
            sys.exit(2)
    # work done outside the session (--report/--methods/--code-dir/--extra-dir): what cannot be
    # folded in is refused here, with nothing written; a report the skill did not write is a gap
    import custom_work
    try:
        custom = custom_work.plan(p["session_dir"], p["output_dir"], getattr(a, "report", None),
                                  getattr(a, "methods", None), getattr(a, "code_dir", None),
                                  getattr(a, "extra_dir", None), getattr(a, "custom_max_gb", None))
    except custom_work.Refused as e:
        sys.stderr.write(f"[{FINALIZE_TOOL}] refusing, nothing written (the earlier finalize "
                         f"stamp is kept): {e}\n")
        sys.exit(2)
    # The stamp says what the LAST finalize made: a finalize that runs ends it first, so one that
    # stops or refuses for a finish gap (below) never leaves an earlier stamp to vouch for what
    # is there now.
    old_stamp = os.path.join(p["output_dir"], rf.FINALIZED_RECORD)
    if os.path.lexists(old_stamp):
        try:
            os.remove(old_stamp)
        except OSError as e:
            sys.stderr.write(f"[session] cannot remove the earlier finalize stamp {old_stamp} "
                             f"({e.strerror or e}); nothing was changed\n")
            sys.exit(2)
    need = rf.FINISH_NEEDS["finalize"]
    gaps = rf.finish_gaps(p["output_dir"], need)
    warnings = rf.finish_warnings(p["output_dir"], need)
    # logs/commands.log missing or behind the analysis's steps: a warning, never a gap
    cmd_log = commands_log_warning(p)
    if cmd_log:
        warnings.append(cmd_log)
    # The reproducibility bundle: finalize builds it when it is missing or stale (L5b review H1,
    # Ruling 59) -- a gap only when a record it needs is absent, with the exact command.
    bundle = _bundle_plan(p)
    if bundle["missing"]:
        gaps.append(rf.gap_line(rf._ITEM_LABEL["reproducibility"],
                                f"{bundle['why']}; finalize cannot build it without "
                                + " and ".join(bundle["missing"]), bundle["command"]
                                + " (step 10), or give finalize the record and run it again"))
    elif bundle["build"] and speak:
        sys.stderr.write(f"[session] the reproducibility bundle: {bundle['why']} -- finalize "
                         "builds it (provenance.py, with step 10's flags read from the session's "
                         "records)\n")
    tgap = None if no_tr else transcript_gap(p, getattr(a, "call_id", None), try_save=speak)
    if tgap:
        gaps.append(tgap)
    gaps += custom_work.report_gaps(custom)
    if speak:
        for w in warnings:
            sys.stderr.write(f"[session] {w}\n")
    if gaps and not reason:
        sys.stderr.write(rf.refusal(FINALIZE_TOOL, gaps) + "\n")
        sys.exit(2)
    if gaps and speak:
        sys.stderr.write(f"[session] going ahead under --override-finish ({reason}) without: "
                         + "; ".join(gaps) + "\n")
    return {"gaps": gaps, "warnings": warnings, "override_reason": reason,
            "override": rf.override_record(FINALIZE_TOOL, gaps, reason) if gaps else None,
            "no_transcript": no_tr, "custom": custom, "build_bundle": bool(bundle["build"]),
            "commands_log": cmd_log}


def _bundle_plan(p):
    """make_deposit.bundle_plan, or -- when make_deposit cannot be loaded (a partial scripts/) --
    the bundle's state with no way to build it: then a missing bundle is a gap, as before."""
    try:
        import make_deposit
        return make_deposit.bundle_plan(p["session_dir"])
    except Exception as e:                       # said in the gap, never swallowed
        import report_files as rf
        state, why = rf.bundle_state(p["output_dir"])
        return {"state": state, "why": why, "build": False,
                "missing": [] if state == "present" else
                [f"make_deposit.py, which could not be loaded ({type(e).__name__}: {e})"],
                "command": rf._command("reproducibility", p["session_dir"])}


def _override_line(name, note):
    """MANIFEST.txt's line for a step taken under --override-finish: [OVERRIDE], in
    make_deposit.Manifest's layout."""
    note = " ".join(str(note).split())
    return f"[OVERRIDE] {name:<50} -- {note}"


def _manifest_reason(manifest_txt, part):
    """The reason MANIFEST.txt gives for `part` -- a [SKIPPED] line first, else an [INFO] one
    (a README.md kept as written by hand has an [INFO] line too) -- or None."""
    try:
        with open(manifest_txt, encoding="utf-8", errors="replace") as fh:
            lines = [ln for ln in fh if part in ln and " -- " in ln]
    except OSError:
        return None                              # no MANIFEST.txt: the caller says "not written"
    for level in ("[SKIPPED]", "[INFO]"):
        for ln in lines:
            if ln.startswith(level):
                return ln.split(" -- ", 1)[1].strip()
    return None


def stamp_outputs(p, written):
    """The stamp's outputs: {session-relative path: sha256 or "[SKIPPED] <reason>"}, and the
    required ones that are missing: [(path, reason)]. AGENTS.md and the READMEs count only when
    THIS finalize wrote them (`written`, session_docs' paths): an older copy describes an earlier
    state (init's starter README.md, for one). One written by hand (report_files.
    skill_doc_targets) is stamped by its checksum -- it is in the zip, and kept as it is -- and
    the skill's copy beside it (README.skill.md, ...) is the one this finalize must have written."""
    import report_files as rf
    done = {os.path.abspath(w) for w in written if w}
    kept = rf.skill_doc_targets(p["session_dir"])["kept"]
    docs = set(LATE_DOCS) | set(rf.SKILL_DOC_TWIN.values())
    outputs, missing = {}, []
    for rel, part, required in rf.STAMP_OUTPUTS:
        if rel in kept:
            try:
                outputs[rel] = rf._sha256(os.path.join(p["session_dir"], rel))
            except OSError as e:
                outputs[rel] = f"[SKIPPED] written by hand, and cannot be read ({e.strerror or e})"
            rel = rf.SKILL_DOC_TWIN[rel]
        path = os.path.join(p["session_dir"], rel)
        ok = os.path.isfile(path) and (rel not in docs or os.path.abspath(path) in done)
        if ok and rel == "output/reproducibility/run_manifest.json":
            ok = rf.bundle_state(p["output_dir"])[0] == "present"     # current, not just there
        if ok:
            outputs[rel] = rf._sha256(path)
            continue
        why = (_manifest_reason(p["manifest_txt"], part) if part else None) or (
            "not written by this finalize" if rel in docs else "not in the session")
        outputs[rel] = f"[SKIPPED] {why}"
        if required:
            missing.append((rel, why))
    import custom_work                           # work folded in (finalize --code-dir ...)
    record = custom_work.record_path(p["output_dir"])
    if os.path.isfile(record):
        outputs[os.path.relpath(record, p["session_dir"]).replace(os.sep, "/")] = rf._sha256(record)
    return outputs, missing


def write_stamp(p, gate, written, zip_path, missing):
    """output/FINALIZED.json, finalize's last act (after MANIFEST.txt is complete, whose checksum
    it holds). -> its path."""
    import report_files as rf
    outputs, _ = stamp_outputs(p, written)
    try:
        import skill_version as sv
        v = sv.skill_version()
        v = None if v == sv.UNKNOWN else v
    except ImportError:
        v = None                                 # a partial scripts/: the stamp says null
    override = gate.get("override")
    if missing and gate.get("override_reason"):
        owed = [f"{rel}: {why}" for rel, why in missing]
        override = rf.override_record(FINALIZE_TOOL, (override or {}).get("missing", []) + owed,
                                      gate["override_reason"])
    stamp = {"schema": rf.FINALIZED_SCHEMA, "finalized_at": rf._now(), "skill_version": v,
             "session": p["session_dir"], "outputs": outputs, "zip": zip_path,
             "finish_override": override, "no_transcript": gate.get("no_transcript"),
             "warnings": gate.get("warnings") or []}
    path = os.path.join(p["output_dir"], rf.FINALIZED_RECORD)
    rf._write_json(path, stamp)
    return path


COMMANDS_LOG_PART = "Command log (logs/commands.log)"


def step_records(p):
    """The records the analysis's own steps write, each as its step runs (never finalize's own
    outputs): the workflow manifest, the FASTA sidecar, every search record, the DE record and
    the report's record. Existing files only, session-relative path -> mtime."""
    import report_files as rf
    cands = ([p["workflow_manifest"], p["fasta_meta"]] + list(search_choice(p)["candidates"])
             + [os.path.join(p["de_dir"], rf.DE_RECORD),
                os.path.join(p["output_dir"], rf.SKILL_REPORT_RECORD)])
    return {os.path.relpath(c, p["session_dir"]): os.path.getmtime(c)
            for c in dict.fromkeys(cands) if c and os.path.isfile(c)}


def commands_log_warning(p):
    """logs/commands.log against the analysis's steps (SKILL.md golden rule 5: every command
    appended as it is run): a WARNING when it is missing or empty, or older than the newest step
    record (step_records) -- commands run since are then missing from it; None when it is at
    least as new. A warning only, never a finish gap: finalize says it (stderr, its JSON's
    finish_warnings, MANIFEST.txt [WARN], the stamp) and goes on (2.14.2 C2, E2E F11)."""
    import datetime
    log = p["commands_log"]
    # (a bundle already made is kept by a later finalize, so it is not promised a later copy)
    fix = ("append each command to it as it is run (SKILL.md golden rule 5); the reproducibility "
           "bundle copies the log as it is when the bundle is built")
    if not os.path.isfile(log) or os.path.getsize(log) == 0:
        return (f"WARNING: logs/commands.log is {'empty' if os.path.isfile(log) else 'missing'}"
                f", so the commands of this analysis are in neither its record nor its bundle: "
                f"{fix}")
    recs = step_records(p)
    if not recs:
        return None
    newest = max(recs, key=recs.get)
    if os.path.getmtime(log) >= recs[newest]:
        return None

    def when(t):
        return datetime.datetime.fromtimestamp(t).isoformat(sep=" ", timespec="minutes")
    return (f"WARNING: logs/commands.log (last written {when(os.path.getmtime(log))}) is older "
            f"than {newest} ({when(recs[newest])}), the newest step record: the commands run "
            f"since may be missing from it -- {fix}")


def do_finalize(a):
    p = paths_for(a.dir)
    if not os.path.isdir(p["session_dir"]):
        sys.exit(f"session dir not found: {p['session_dir']}")
    # The finish gate first (2.14.2): nothing is written, and nothing zipped, for a session whose
    # report the skill did not write, with no reproducibility bundle, or a Claude Code analysis
    # with no saved conversation -- unless overridden with a reason. main() already ran it once
    # before a compute-node re-exec; here it is read again quietly for what it recorded.
    # (the gate also removed an earlier stamp: the stamp is the last act of a finalize that made
    # everything, so one that stops, or refuses to stamp, never leaves a stamp that says otherwise)
    gate = finalize_gate(a, p, speak=os.environ.get(FINISH_GATE_ENV) != p["session_dir"])
    import report_files

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
        if gate["override"]:                    # first: what finalize went ahead without, and why
            man.lines.append(_override_line(f"Finish check ({FINALIZE_TOOL})",
                                            f"{gate['override_reason']} -- went ahead without: "
                                            + "; ".join(g.split(":", 1)[0]
                                                        for g in gate["gaps"])))
        (man.ok if raw_list[0] == "OK" else man.skip)("input/raw_files.txt (where the raw data are)",
                                                       raw_list[1])
        # The reproducibility bundle (step 10): kept when current, else built now -- after the
        # raw list it reads, before the work folded in below, whose checksums go into it
        # (L5b review H1: finalize MAKES it, for a session with a DE or without)
        try:
            make_deposit.ensure_bundle(p["session_dir"], man)
        except Exception as e:                  # recorded, never swallowed
            man.skip(make_deposit.BUNDLE_PART, f"{type(e).__name__}: {e}")
    # Work done outside the session (--report/--methods/--code-dir/--extra-dir): copied into
    # output/custom/ now, before the documents that list it (custom_work.py)
    custom_lines = custom_step(p, gate)
    if man is not None:
        search_record_step(p, man)
        if gate.get("commands_log"):          # the gate's warning, in MANIFEST.txt too
            man.warn(COMMANDS_LOG_PART, gate["commands_log"].replace("WARNING: ", "", 1))
        for level, part, note in custom_lines:
            {"OK": man.ok, "SKIPPED": man.skip, "NOTE": man.note}.get(level, man.info)(part, note)
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
        de_only_step(p, man)
        peak_check_step(p, man)
        podcast_share_step(p, man)             # before the README, which lists it
        conversation_step(p, man, getattr(a, "call_id", None), gate["no_transcript"])
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
        _append_manifest(p["manifest_txt"], "SKIPPED", "Reproducibility bundle (output/"
                         "reproducibility)", "not built: " + import_error)
        for level, part, note in custom_lines + docs["manifest"]:
            _append_manifest(p["manifest_txt"], level, part, note)
        if gate["override"]:
            with open(p["manifest_txt"], "a", encoding="utf-8") as fh:
                fh.write(_override_line(f"Finish check ({FINALIZE_TOOL})",
                                        f"{gate['override_reason']} -- went ahead without: "
                                        + "; ".join(g.split(":", 1)[0] for g in gate["gaps"]))
                         + "\n")
    # The finalize stamp (Brett, 2026-10-31): AGENTS.md, README.md and the Methods must have been
    # written by this finalize, or it refuses to stamp (unless --override-finish) -- and then
    # the Core channel is not told the analysis is complete either. Its MANIFEST line goes in
    # now, before the zip; the stamp itself is written last.
    stamp_name = f"Finalize stamp (output/{report_files.FINALIZED_RECORD})"
    _, stamp_missing = stamp_outputs(p, list(docs["written"].values()))
    not_finalized = None
    if stamp_missing and not gate["override_reason"]:
        not_finalized = "; ".join(f"{rel} was not written ({why})" for rel, why in stamp_missing)
        _append_manifest(p["manifest_txt"], "SKIPPED", stamp_name, f"not written: {not_finalized}")
    elif stamp_missing:
        with open(p["manifest_txt"], "a", encoding="utf-8") as fh:
            fh.write(_override_line(stamp_name, f"written under --override-finish "
                                                f"({gate['override_reason']}) without: "
                                    + ", ".join(rel for rel, _ in stamp_missing)) + "\n")
    else:
        _append_manifest(p["manifest_txt"], "OK", stamp_name,
                         "written as finalize's last act, with this file's checksum")

    de_files = sorted(os.path.basename(f) for f in glob.glob(os.path.join(p["de_dir"], "DE_*.csv")))
    result = {"session_dir": p["session_dir"], "readme": docs["written"].get("readme") or p["readme"],
              "readme_html": docs["written"].get("readme_html"),
              "agents": docs["written"].get("agents"), "de_files": de_files,
              "methods": methods_md, "manifest": p["manifest_txt"],
              "deposit": deposit, "skipped": (man.n_skipped if man is not None else 1)}
    # for the agent, never the zip: make_methods' database build warnings, when finalize built the
    # Methods (make_deposit.methods_build_warnings) -- tell the user
    result.update(getattr(man, "for_agent", None) or {})

    # re-analysis: write DIFFERENCES.md vs the parent (read, never written to)
    parent = a.reanalysis_of or (reanalysis_of(p["session_dir"]) or {}).get("session")
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
                 # the analysis conversation (save_transcript.py): Core-internal -- it holds
                 # internal paths and remarks, and a zip travels
                 os.path.abspath(os.path.join(sdir, "logs", "conversation")):
                 "logs/conversation (the analysis conversation: Core-internal, kept on disk)",
                 # ... and the decisions log, internal for the same reason (log_decision.py)
                 os.path.abspath(os.path.join(sdir, "logs", "decisions.md")):
                 "logs/decisions.md (the decisions log: Core-internal, kept on disk)",
                 os.path.abspath(os.path.join(p["deposit_dir"], "upload_staging")):
                 "output/DATA_SUBMISSION/upload_staging (raw archives staged for upload)"}
        # The rest of logs/ is Core-internal too (I11): expert reviews, job logs, an agent-named
        # copy of qc_bracket's output with the verdict and no marker. Only logs/commands.log goes
        # in, redacted (below). The two labels above keep their own counts; this one counts the
        # files the walk would otherwise have written (a staff-only or scratch file keeps its own).
        logs_dir = os.path.abspath(p["logs_dir"])
        logs_label = ("logs/ (Core-internal: reviews, job logs, decisions -- kept on disk; "
                      "commands.log is in, redacted)")
        n_logs = 0
        archive = sdir + ".zip"
        excluded = {label: 0 for label in skips.values()}
        manifest_txt = os.path.abspath(p["manifest_txt"])
        # README.md, README.html and AGENTS.md go in after the run-log hook too: they name the
        # Core registry record, which record_run.py may only create then.
        # (the skill's: one written by hand is zipped where it is, like any other file)
        late = {os.path.abspath(x) for x in skill_doc_paths(sdir)}
        # A skill copy (README.skill.md ...) left from an earlier finalize, when its original was
        # hand-written and is not now: it describes an earlier state, so it stays out (L5b N4)
        import report_files as _rf
        superseded = {os.path.abspath(os.path.join(sdir, twin))
                      for twin in _rf.SKILL_DOC_TWIN.values()} - late
        superseded_label = ("earlier skill copies of the README / AGENTS.md (README.skill.md "
                            "...; the original is the skill's again -- kept on disk)")
        n_superseded = 0
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
        # Staff-only records (staff.py, *.staff.json): who decided, by name and login -- wherever
        # they are, never in a zip that goes to collaborators.
        import staff
        staff_label = "staff-only records (*.staff.json: who decided, by name -- kept on disk)"
        n_staff = 0
        # The QC verdict around the project (qc_bracket.py) is staff-only: ANY file carrying its
        # marker, wherever it is and whatever it is called (--out can name it anything), and
        # deliver's delivery.json (its record of a delivery, staff-side) stay out. Without
        # qc_bracket.py beside this (a partial copy), its file names are the rule.
        qc_label = ("staff-only QC records (qc_bracket.py's marker, any name, anywhere) and "
                    "delivery.json (kept on disk)")
        n_qc = 0
        # Core-internal by name wherever they are (staff.is_core_internal, the rule deliver uses
        # too; C2 security M1): the scheduler's job logs and the Core's own records
        internal_label = ("Core-internal files outside logs/ (scheduler job logs <job>_<id>.log, "
                          "fran_deposit.json, node_faults.json -- kept on disk)")
        n_internal = 0
        try:
            from qc_bracket import is_staff_only_file
        except ImportError:
            def is_staff_only_file(path):
                return bool(re.match(r"^(?:qc_bracket(?:_ack)?\.(?:json|md)|qc_gate\.json)$",
                                     os.path.basename(path)))

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
                    if staff.is_staff_only(fn):
                        n_staff += 1
                        continue
                    if os.path.abspath(full) in skips:           # a file kept out (decisions.md)
                        excluded[skips[os.path.abspath(full)]] = 1
                        continue
                    if os.path.abspath(full) in superseded:
                        n_superseded += 1
                        continue
                    if fn == "delivery.json" or (not os.path.islink(full)
                                                 and is_staff_only_file(full)):
                        n_qc += 1
                        continue
                    if os.path.islink(full) or os.path.abspath(full) == manifest_txt \
                            or os.path.abspath(full) in late:
                        continue                 # MANIFEST.txt and the docs go in last, below
                    if (os.path.abspath(full).startswith(logs_dir + os.sep)
                            and os.path.abspath(full) != os.path.abspath(p["commands_log"])):
                        n_logs += 1              # Core-internal: kept on disk, never zipped
                        continue
                    if staff.is_core_internal(fn):
                        n_internal += 1
                        continue
                    redacted = staff.client_copy(full)
                    if redacted is not None:     # who accepted, as the role (staff.client_copy)
                        z.writestr(os.path.join(base, os.path.relpath(full, sdir)), redacted)
                        continue
                    if fn == os.path.basename(p["commands_log"]):
                        # every command verbatim (logs/, or a copy anywhere), but who decided is
                        # staff-only (staff.redact)
                        with open(full, encoding="utf-8", errors="replace") as fh:
                            z.writestr(os.path.join(base, os.path.relpath(full, sdir)),
                                       staff.redact(fh.read()))
                        continue
                    z.write(full, os.path.join(base, os.path.relpath(full, sdir)))
        excluded[quant_label] = n_quant
        excluded[speclib_label] = n_speclib
        excluded[qc_label] = n_qc
        excluded[logs_label] = n_logs
        if n_internal:
            excluded[internal_label] = n_internal
        if n_scratch:
            excluded[scratch_files.LABEL] = n_scratch
        if n_staff:
            excluded[staff_label] = n_staff
        if n_superseded:
            excluded[superseded_label] = n_superseded
        result["zip"] = archive
        result["zip_excluded"] = excluded

    # The run log and Slack go once everything they point at exists -- after the zip. Their
    # outcomes are the last two lines of MANIFEST.txt, so the zip's copy (added now) says them.
    try:
        with open(p["manifest_txt"], encoding="utf-8") as fh:
            pre_hook = fh.read()
    except OSError:
        pre_hook = None
    # a finalize that went ahead without something (--override-finish) is announced as such
    hooks = _finish_hooks(a, p["session_dir"], result.get("zip"), not_finalized=not_finalized,
                          finish_override=(gate["override_reason"]
                                           if gate["override"] or stamp_missing else None))
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

    # Last: the stamp. `session.py status`, deliver, and the outward steps after it read it.
    result["finish_override"] = gate["override"]
    result["finish_warnings"] = gate["warnings"]
    if not_finalized:
        result["finalized"] = None
        result["not_finalized"] = [f"{rel}: {why}" for rel, why in stamp_missing]
    else:
        try:
            result["finalized"] = write_stamp(p, gate, written, result.get("zip"), stamp_missing)
        except OSError as e:
            result["finalized"] = None
            result["not_finalized"] = [f"output/{report_files.FINALIZED_RECORD} could not be "
                                       f"written: {e.strerror or e}"]
    print(json.dumps(result, indent=2))
    if not result.get("finalized"):
        hints = _not_finalized_hints(p, a, [rel for rel, _ in stamp_missing])
        sys.stderr.write("[session] NOT FINALIZED -- no output/FINALIZED.json: "
                         + "; ".join(result["not_finalized"]) + ". Make what is missing and run "
                         "finalize again (or --override-finish \"<reason>\")."
                         + "".join(f"\n  - {h}" for h in hints) + "\n")
        sys.exit(2)


def _not_finalized_hints(p, a, missing):
    """How to make each document the stamp still needs (L5b review L5: a --methods fold-in with
    no raw-file list was refused with no way forward). -> [sentence]."""
    import shlex
    import report_files as rf
    q = shlex.quote
    hints = []
    if "output/methods.md" in missing:
        given = getattr(a, "methods", None)
        hints.append("output/methods.md: "
                     + (f"use the Methods you gave with --methods as the session's own: cp "
                        f"{q(os.path.abspath(os.path.expanduser(given)))} {q(p['methods_md'])} "
                        "(finalize keeps a methods.md that is there); or " if given else "")
                     + "write input/raw_files.txt (one raw file path per line) so make_methods.py "
                       f"can write it ({rf._command('methods', p['session_dir'])})")
    if "output/reproducibility/run_manifest.json" in missing:
        hints.append("output/reproducibility/run_manifest.json: "
                     + rf._command("reproducibility", p["session_dir"]))
    for rel in missing:
        if rel not in ("output/methods.md", "output/reproducibility/run_manifest.json"):
            hints.append(f"{rel}: MANIFEST.txt's [SKIPPED] line for it says why it could not be "
                         "written")
    return hints


STATUS_COMPLETE, STATUS_INCOMPLETE = 0, 3


def finish_status(session_dir):
    """Every finish item of a session, in order: [{"item", "label", "done", "state", "detail",
    "command"}]. "done" is True for an item that is done or not needed here. The items: the
    skill's report, Methods, the bundle, AGENTS.md, the README, the saved conversation (a Claude
    Code analysis only), the finalize stamp (current: nothing it stamped has changed since), and
    the delivery (a CoreOmics submission only; delivery.json `owed`, an interim delivery's list of
    what is still to come, keeps it open)."""
    import shlex
    import report_files as rf
    p = paths_for(session_dir)
    s, out = p["session_dir"], p["output_dir"]
    q = shlex.quote(s)
    finalize_cmd = f"python3 scripts/session.py finalize --dir {q} --zip (step 12)"
    items = []
    # What the last finalize waived (--no-transcript) or went ahead without (--override-finish),
    # with the reason, from its stamp: such an item is done, said as such (L5a review F2). The
    # stamp's `finish_override.missing` names each gap as finalize did: a finish label
    # (report_files), the conversation's sentence, or "<path>: <why>" for a document not written.
    rec = rf.finalized_record(out)[0] or {}
    fo = rec.get("finish_override") if isinstance(rec.get("finish_override"), dict) else {}
    waived = {"conversation": rec.get("no_transcript")}
    named = {"report": (rf._ITEM_LABEL["report"],),
             "methods": (rf._ITEM_LABEL["methods"], "output/methods.md:"),
             "reproducibility": (rf._ITEM_LABEL["reproducibility"],
                                 "output/reproducibility/run_manifest.json:"),
             "agents": ("AGENTS.md:", "AGENTS.skill.md:"),
             "readme": ("README.md:", "README.skill.md:", "README.html:", "README.skill.html:"),
             "conversation": ("the analysis conversation (logs/conversation)",)}

    def add(item, label, done, state, detail, command=None):
        reason = None
        if not done and waived.get(item):
            done, state, reason = True, "waived", waived[item]
            detail = f"waived ({reason})"
        elif not done and any(str(g).lower().startswith(n.lower()) for g in fo.get("missing") or []
                              for n in named.get(item, ())):
            done, state, reason = True, "overridden", fo.get("reason")
            detail = f"overridden ({reason})"
        items.append({"item": item, "label": label, "done": bool(done), "state": state,
                      "detail": detail, "command": None if done else command,
                      **({"reason": reason} if reason is not None else {})})

    for e in rf.finish_check(out, ("report", "methods", "reproducibility")):
        add(e["item"], e["label"][0].upper() + e["label"][1:], e["ok"], e["state"],
            (e["warning"] or e["why"]) if e["ok"] else e["why"], e["command"])
    # the skill's own documents: beside one written by hand, its .skill copy (skill_doc_targets)
    t = rf.skill_doc_targets(s)
    hand = lambda *names: "; ".join(f"{n} is written by hand, so the skill's is "   # noqa: E731
                                    f"{os.path.basename(t[n])}" for n in names if n in t["kept"])
    # required as the stamp requires them (STAMP_OUTPUTS, one definition): README.html is not
    required = {rel: req for rel, _part, req in rf.STAMP_OUTPUTS}

    def doc(item, label, names, why_missing):
        absent = [n for n in names if not os.path.isfile(t[n])
                  or (n == "README.md" and rf.is_starter_readme(t[n]))]
        done = not [n for n in absent if required.get(n)]
        optional = "; ".join(f"{os.path.basename(t[n])} not written (optional)" for n in absent)
        add(item, label, done, "present" if done else "missing",
            ("present" + (f"; {optional}" if optional else "")) if done else
            ("not written yet " + (f"({hand(*names)})" if hand(*names) else why_missing)).rstrip(),
            finalize_cmd)

    doc("agents", "AGENTS.md (a guide for an AI agent)", ("AGENTS.md",), "")
    doc("readme", "README.html + README.md (for collaborators)", ("README.md", "README.html"),
        "(init's starter README.md is not it)")
    sid, where = claude_code_record(s)
    if sid:
        try:
            import save_transcript
            n = save_transcript.saved_count(s)
        except ImportError as e:
            n, why = 0, f"save_transcript.py could not be loaded ({e})"
        else:
            why = f"a Claude Code analysis ({where}) with no saved conversation"
        add("conversation", "The analysis conversation (logs/conversation)", n > 0,
            "saved" if n else "missing", f"{n} saved" if n else why,
            f"python3 scripts/save_transcript.py {q} (from the laptop, for a HIVE session: "
            f"python3 scripts/save_transcript.py --hive {q}), or finalize --no-transcript "
            "\"<why>\"")
    else:
        dec = os.path.isfile(os.path.join(p["logs_dir"], "decisions.md"))
        add("conversation", "The analysis conversation (logs/conversation)", True, "not_needed",
            "not a Claude Code analysis" + ("; the decisions log (logs/decisions.md) is the "
                                            "record" if dec else ""))
    # the stamp, current or not: report_files.finish_check, the one rule deliver uses too
    st = rf.finish_check(out, ("finalized",))[0]
    add("finalized", "The finalize stamp (output/FINALIZED.json)", st["ok"], st["state"],
        st["why"] + (" under --override-finish" if st["ok"] and fo else ""), finalize_cmd)
    core = None
    try:
        import submission_report
        core = submission_report.core_run(s)
    except Exception as e:                        # said in the item, never swallowed
        add("delivery", "The delivery (core_submission.py deliver)", False, "unknown",
            f"cannot tell whether this is a CoreOmics submission: {type(e).__name__}: {e}",
            "python3 scripts/submission_report.py ... (see SKILL.md step 1c)")
        return items
    if not core:
        add("delivery", "The delivery (core_submission.py deliver)", True, "not_needed",
            "not a CoreOmics submission: the session itself is the result")
        return items
    dj = os.path.join(s, "delivery.json")
    d = _load(dj) if os.path.isfile(dj) else None
    deliver_cmd = (f"python3 scripts/core_submission.py deliver --summary <submission_summary.json> "
                   f"--session {q} --mode <analysis|no-de> --apply (step 12b)")
    if not isinstance(d, dict):
        add("delivery", f"The delivery of {core.get('prot') or 'the submission'}", False,
            "missing", "not delivered yet" if not os.path.isfile(dj) else
            "delivery.json cannot be read", deliver_cmd)
    elif d.get("verified") is not True:
        # before `owed`: an interim that did not verify is not "delivered" (C2 errors_ux L2)
        add("delivery", f"The delivery of {core.get('prot') or 'the submission'}", False,
            "unverified", "the last delivery did not verify (delivery.json `verified` is not "
                          "true): do not share it", deliver_cmd)
    elif d.get("owed"):
        add("delivery", f"The delivery of {core.get('prot') or 'the submission'}", False,
            "owed", f"{d.get('mode') or 'interim'} delivered; owed: "
            + ", ".join(str(x) for x in d["owed"]), deliver_cmd)
    elif rf.delivery_finish_override(d):
        # delivered under --override-finish (deliver's, or the stamp's): done, said as such --
        # never "every item is done" (C2 errors_ux M1), as for finalize's own override
        dfo = rf.delivery_finish_override(d)
        items.append({"item": "delivery",
                      "label": f"The delivery of {core.get('prot') or 'the submission'}",
                      "done": True, "state": "overridden",
                      "detail": (f"{d.get('mode')} delivery at {d.get('delivered_at')}, NOT "
                                 f"FINALIZED: under --override-finish "
                                 f"({dfo.get('reason') or 'no reason recorded'})"),
                      "command": None, "reason": dfo.get("reason")})
    elif (d.get("mode") not in rf.UNGATED_DELIVERY_MODES and isinstance(d.get("finalized"), dict)
          and d["finalized"].get("sha256")
          and d["finalized"]["sha256"] != _sha256_or_none(os.path.join(out,
                                                                        rf.FINALIZED_RECORD))):
        # finalized again since: the delivered package is not the one finalized now (M1)
        add("delivery", f"The delivery of {core.get('prot') or 'the submission'}", False,
            "stale", f"{d.get('mode')} delivery at {d.get('delivered_at')} is of the finalize "
                     f"at {d['finalized'].get('finalized_at')}; the session was finalized again "
                     "since -- deliver again", deliver_cmd)
    else:
        add("delivery", f"The delivery of {core.get('prot') or 'the submission'}", True,
            "delivered", f"{d.get('mode')} delivery at {d.get('delivered_at')}")
    return items


def _sha256_or_none(path):
    import report_files as rf
    try:
        return rf._sha256(path)
    except OSError:
        return None


def do_status(a):
    """`session.py status --session <S>`: the finish checklist (golden rule 9: run it before
    telling the user an analysis is finished, and show its result). Exit 0 when every item is
    done, 3 when one is not."""
    p = paths_for(a.session)
    if not os.path.isdir(p["session_dir"]):
        sys.stderr.write(f"session dir not found: {p['session_dir']}\n")
        sys.exit(2)
    items = finish_status(p["session_dir"])
    missing = [i for i in items if not i["done"]]
    # done only because finalize was told to go ahead without it (--override-finish) or to waive
    # it (--no-transcript): never "every item is done" while one is (L5b review H1)
    excused = [i for i in items if i["done"] and i["state"] in ("overridden", "waived")]
    if a.json:
        print(json.dumps({"session_dir": p["session_dir"], "complete": not missing,
                          "overridden_or_waived": excused, "items": items}, indent=2))
    else:
        print(f"Finish checklist for {p['session_dir']}:")
        for i in items:
            if i["done"]:
                print(f"  done     {i['label']} -- {i['detail']}")
            else:
                print(f"  missing  {i['label']} -- {i['detail']}: {i['command']}")
        if missing:
            print(f"NOT complete: {len(missing)} item(s) missing. Finish each one before telling "
                  "the user the analysis is finished (SKILL.md golden rule 9).")
        elif excused:
            print(f"Complete, with {len(excused)} overridden/waived: " + "; ".join(
                f"{i['label']} ({i['state']}: {i.get('reason')})" for i in excused))
        else:
            print("Complete: every item is done.")
    sys.exit(STATUS_COMPLETE if not missing else STATUS_INCOMPLETE)


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
    `resolved`, the cfg the search wrote as it ran (<search out>/params.resolved.cfg, the search
    out being the folder of the session's search record -- search_record: output/search_sage/ as
    much as output/search/) first; then the one the workflow step staged (DIA-NN cfg / Sage json
    / FragPipe .workflow / Radiant config), newest convention first; then a params.* / *.cfg /
    sage_config*.json directly in input/ (older sessions). resolved=False is the file the search
    was GIVEN (make_deposit's Methods read)."""
    wf, inp = p["workflow_dir"], p["input_dir"]
    cands = []
    search_out = search_folder(p) if resolved else None
    if search_out:
        cands.append(os.path.join(search_out, "params.resolved.cfg"))
    cands += [os.path.join(wf, n) for n in ("params.cfg", "params.json")]
    for pat in ("*.cfg", "params*.json", "sage*.json", "*.workflow", "*.radiantConfig"):
        cands += sorted(glob.glob(os.path.join(wf, pat)))
    for pat in ("params.*", "*.cfg", "sage_config*.json"):
        cands += sorted(glob.glob(os.path.join(inp, pat)))
    import staff                                 # a cfg's staff-only record is not a cfg
    return next((c for c in cands if os.path.isfile(c) and not staff.is_staff_only(c)
                 and not c.endswith((".rationale.json", "manifest.json"))), None)


SEARCH_RECORD = "search_provenance.json"
# search_choice's `how` when there are records and none is chosen: every reader treats the two
# alike (no search described, and the reason said) -- search_unchosen()
SEARCH_UNCHOSEN = ("ambiguous", "unrecorded_input")
# How the user chooses then (T44 follow-up N1): the refusals and SKILL.md step 12b say it
SEARCH_HOW_TO_CHOOSE = ("to choose one, rename each other search folder to "
                        "output/superseded_<name>/ (a folder whose name holds \"superseded\" is "
                        "not a candidate), or deliver every search folder with "
                        "core_submission.py deliver --mode interim")


def search_folder_of(report):
    """THE search folder of a report (a DE's input): the report's own folder, or its parent when
    only the parent holds search_provenance.json -- FragPipe's DIA route writes its report to
    <workdir>/dia-quant-output/ -- else the report's folder. contaminants.R search_folder() is
    its R twin (run_de.R reads the search's record through it)."""
    d = os.path.dirname(os.path.abspath(str(report)))
    for cand in (d, os.path.dirname(d)):
        if os.path.isfile(os.path.join(cand, SEARCH_RECORD)):
            return cand
    return d


def search_records(output_dir):
    """Every search record of a session: output/*/search_provenance.json and
    output/*/*/search_provenance.json, one per real path, leaving out any folder whose name says
    "superseded"."""
    cands, seen = [], set()
    for depth in (("*",), ("*", "*")):
        for path in sorted(glob.glob(os.path.join(output_dir, *depth, SEARCH_RECORD))):
            rel = os.path.relpath(os.path.dirname(path), output_dir)
            if any("superseded" in part.lower() for part in rel.split(os.sep)):
                continue
            real = os.path.realpath(path)
            if real not in seen:
                seen.add(real)
                cands.append(path)
    return cands


def choose_search(output_dir, de_input=None):
    """THE answer to "which search" for a session's output/ (DE-LIMP rule 3) -- every reader of a
    search record gives this one, through search_choice(p) or find_search_provenance -> dict:
      path        the chosen search_provenance.json, or None
      candidates  every record (search_records)
      how         "de_input": the record of the folder the DE read (search_folder_of its input:
                  the report's folder, or its parent for FragPipe's DIA route);
                  "only": the one record, and no DE input here says otherwise;
                  "none": no record at all (readers fall back to the default output/search/ and
                  the workflow manifest's pin, tagged);
                  "unrecorded_input": the DE read a report that is here and that no record
                  describes (a recordless output/search/ beside a recorded search) -- never the
                  lone record (T44 review F2);
                  "ambiguous": several records and the DE's input names none of them (no DE, an
                  older DE record without `input`, or an input that is not here -- a laptop copy
                  of a HIVE session) -- never output/search/'s by default (F3)
      why         for "none", "unrecorded_input" and "ambiguous", one sentence naming the records;
                  None when a record is chosen
      short       the same in a few words, with no paths (a client-facing page says this)
    A search in output/search_dda/ was once missed, and finalize described the session's params
    file instead (2.14.2, I04)."""
    cands = search_records(output_dir)
    sd = os.path.dirname(os.path.abspath(output_dir))

    def rel(x):
        r = os.path.relpath(os.path.abspath(str(x)), sd)
        return str(x) if r.startswith("..") else r.replace(os.sep, "/")

    names = ", ".join(rel(c) for c in cands)
    res = {"path": None, "candidates": cands, "how": None, "why": None, "short": None}
    if de_input:
        folder = os.path.realpath(search_folder_of(de_input))
        hit = [c for c in cands if os.path.realpath(os.path.dirname(c)) == folder]
        if hit:
            return dict(res, path=hit[0], how="de_input")
        if cands and os.path.exists(str(de_input)):
            return dict(res, how="unrecorded_input", why=(
                f"the DE read {rel(de_input)}, which no search record describes (the session's: "
                f"{names})"), short="the search the DE read has no search record")
    if len(cands) == 1:
        return dict(res, path=cands[0], how="only")
    if not cands:
        return dict(res, how="none", why="no search record in this session (no output/<folder>/"
                                         "search_provenance.json)", short="no search record")
    return dict(res, how="ambiguous", short=f"{len(cands)} search records, none chosen", why=(
        f"{len(cands)} search records ({names}), none chosen: " +
        (f"the report the DE read ({rel(de_input)}) is in none of their folders" if de_input
         else "the DE's record names no report it read" if de_input is not None
         else "there is no DE to choose by")))


def find_search_provenance(output_dir, de_input=None):
    """(path or None, candidates): choose_search's answer for a session's output/, kept for the
    readers that hold only the output folder. Here, beside paths_for and params_file, so a reader
    that imports session need not import make_deposit (which re-exports it)."""
    c = choose_search(output_dir, de_input)
    return c["path"], c["candidates"]


def de_input_of(p):
    """The report the session's DE read (output/tables/de_provenance.json `input`): the
    tie-break between a session's search records (choose_search). None when there is no DE
    record; "" when the record names no input (an older run_de.R)."""
    de = _load(os.path.join(p["de_dir"], "de_provenance.json"))
    if not isinstance(de, dict):
        return None
    v = de.get("input")
    return v if isinstance(v, str) and v.strip() else ""


def search_choice(p, de_input=None):
    """choose_search for a session (paths_for dict), tie-broken by `de_input` -- by default the
    report the session's own DE read (de_input_of). THE answer every reader of a session's search
    gives: the README and AGENTS.md (session_docs), the report, the Slack posts, the delivery,
    the deposit list, the Methods, the bundle, the run log, the peak check, the resolved
    parameters."""
    return choose_search(p["output_dir"], de_input_of(p) if de_input is None else de_input)


def search_unchosen(choice):
    """True when the session has search records and none is chosen (search_choice `how` in
    SEARCH_UNCHOSEN): no reader describes one of them, and each says so with choice["why"]."""
    return choice.get("how") in SEARCH_UNCHOSEN


def search_record(p, de_input=None):
    """(path or None, candidates): search_choice's record and every candidate."""
    c = search_choice(p, de_input)
    return c["path"], c["candidates"]


ENGINE_PIN_SOURCE = "workflow manifest pin -- the search's own record is not here"


def search_engine(p, de_input=None, record=None):
    """THE engine and version of a session's search, for every reader that names them (the
    report's header and "How this report was made", the Slack post, the README and AGENTS.md, the
    bundle's --engine) -> {"engine", "version", "source", "note", "choice"}:
      a chosen record (search_choice; or `record`, a record the caller read beside the DE's
        input): its engine and its version. A version the record leaves null stays None
        (run_search.engine_version_record: every conda-env Sage) -- never the manifest's pin,
        which is a request, not what ran; source "search_provenance.json".
      records and none chosen: engine and version None, `note` = choice["why"] -- never the
        manifest's engine (T44 review F5).
      no record at all: the workflow manifest's engine and version, source ENGINE_PIN_SOURCE.
    make_methods.search_record applies the same version rule to the Methods."""
    c = search_choice(p, de_input)
    if record is None and c["path"]:
        record = _load(c["path"])
    wf = _load(p["workflow_manifest"]) or {}
    pin = wf.get("engine") if isinstance(wf, dict) and isinstance(wf.get("engine"), dict) else {}
    res = {"engine": None, "version": None, "source": None, "note": None, "choice": c}
    if isinstance(record, dict) and record.get("engine"):
        v = record.get("version")
        return dict(res, engine=str(record["engine"]).lower(), version=str(v) if v else None,
                    source="search_provenance.json")
    if record is None and search_unchosen(c):
        return dict(res, note=c["why"])
    if pin.get("name"):
        return dict(res, engine=str(pin["name"]).lower(),
                    version=str(pin["version"]) if pin.get("version") else None,
                    source=ENGINE_PIN_SOURCE)
    return res


def search_result_files(p, choice=None, engine_tables=False):
    """THE result files of the session's chosen search (search_choice) -> absolute paths, existing
    files only, each once:
      * the report the DE read when it lies in the search folder, else <folder>/report.parquet:
        what a search that quantified writes;
      * identifications.tsv and the record of how its rows were kept (sage_identifications.json,
        fragpipe_identifications.json): what a search for identifications only writes instead
        (run_search.py adapt_sage_ids / adapt_fragpipe_ids -- a search of the other kind into the
        same folder sets them aside, so a folder holds one kind);
      * with `engine_tables`, the engine's own tables those were made from, as the record's
        adapter entries name them (`sage_identifications` / `sage_adapt`: `source`,
        `sage_results` -- Sage's results.sage.parquet, lfq.parquet), found in the search folder.
    With no search record at all, the default folder's; with records and none chosen, none.
    The working copy's pull (search_pull_paths) and the bundle's fingerprints
    (make_deposit.bundle_command) read it (2.14.2 C2, E2E F3/F4)."""
    c = choice or search_choice(p)
    folder = search_folder(p, c)
    if not folder:
        return []
    import run_search                                  # the names its adapters write
    real = os.path.realpath(folder)

    def _inside(f):
        return bool(f) and os.path.isfile(f) and os.path.commonpath(
            [os.path.realpath(f), real]) == real
    inp = de_input_of(p)
    files = [inp if _inside(inp) else os.path.join(folder, "report.parquet")]
    files += [os.path.join(folder, n) for n in (run_search.SAGE_IDS, run_search.SAGE_IDS_RECORD,
                                                run_search.FRAGPIPE_IDS_RECORD)]
    if engine_tables:
        rec = _load(c["path"]) if c["path"] else None
        for key in ("sage_identifications", "sage_adapt"):
            r = rec.get(key) if isinstance(rec, dict) else None
            for k in ("source", "sage_results"):
                v = r.get(k) if isinstance(r, dict) else None
                if isinstance(v, str) and v:
                    # a copied session's record names the original path: the same file here
                    files.append(v if _inside(v) else os.path.join(folder, os.path.basename(v)))
    out = []
    for f in files:
        if f and os.path.isfile(f) and os.path.abspath(f) not in out:
            out.append(os.path.abspath(f))
    return out


def search_pull_paths(p, choice=None):
    """The files a working copy of the session needs to describe its search -- the hive_remote
    route writes the report on the laptop from a copy (SKILL.md "Pull what the report needs") --
    session-relative and '/'-separated, only those that exist: the chosen search's
    search_provenance.json, its results (search_result_files: the report the DE read, or
    <folder>/report.parquet; for a search for identifications only, identifications.tsv and its
    record), and the workflow manifest. With no record at all, the default folder's results;
    with records and none chosen, no search file (search_unchosen)."""
    c = choice or search_choice(p)
    files = ([c["path"]] if c["path"] else []) + search_result_files(p, c)
    files.append(p["workflow_manifest"])
    return [os.path.relpath(f, p["session_dir"]).replace(os.sep, "/") for f in files
            if f and os.path.isfile(f)]


def do_which_search(a):
    """`which-search`: the session's search (search_choice) -- JSON: how it was chosen or why
    none was, the folder, the record, every candidate and the files a working copy pulls; with
    --paths, only those files, one per line. Exit 0 when one is chosen (or the session has no
    record at all), 3 when it has records and none is chosen (the reason on stderr)."""
    p = paths_for(a.dir)
    if not os.path.isdir(p["session_dir"]):
        sys.exit(f"session dir not found: {p['session_dir']}")
    c = search_choice(p)
    rel = lambda x: os.path.relpath(x, p["session_dir"]).replace(os.sep, "/") if x else None  # noqa: E731
    unchosen = search_unchosen(c)
    if a.paths:
        if unchosen:
            sys.stderr.write(f"none chosen, nothing to pull: {c['why']}; "
                             f"{SEARCH_HOW_TO_CHOOSE}\n")
        else:
            print("\n".join(search_pull_paths(p, c)))
    else:
        folder = search_folder(p, c)
        print(json.dumps({"session_dir": p["session_dir"], "how": c["how"],
                          "folder": rel(folder) if folder and os.path.isdir(folder) else None,
                          "record": rel(c["path"]), "why": c["why"], "short": c["short"],
                          "candidates": [rel(x) for x in c["candidates"]],
                          "pull": [] if unchosen else search_pull_paths(p, c)}, indent=2))
    sys.exit(STATUS_INCOMPLETE if unchosen else STATUS_COMPLETE)


def search_folder(p, choice=None):
    """THE search folder of a session: the chosen record's folder (search_choice); the default
    output/search/ when the session has no search record at all (as before 2.14.2); None when it
    has records and none is chosen (search_unchosen -- the reader says choice["why"])."""
    c = choice or search_choice(p)
    if c["path"]:
        return os.path.dirname(c["path"])
    return None if search_unchosen(c) else p["search_out"]


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
    i.add_argument("--call-id", metavar="TOKEN",
                   help="a fresh token (8+ letters and digits) written into this command line, "
                        "for save_transcript.py: it tells which conversation -- the main one or "
                        "a subagent's -- is running it. Needed while subagents are at work")
    i.add_argument("--transcript-from-now", action="store_true",
                   help="this analysis's part of the conversation starts now, not at the "
                        "skill's first load (a conversation about several clients)")
    i.add_argument("--name", required=True, help="short descriptive study name")
    i.add_argument("--base", default="", help="central location for the session (e.g. ~/Documents/DataAnalysis). OMIT to put results in the folder with the raw data (default).")
    i.add_argument("--date", default="", help="YYYY-MM-DD (default: today)")
    i.add_argument("--raw", nargs="*", help="raw file paths/globs; their folder is where results go by default")
    i.add_argument("--reanalysis-of", default="", help="prior session dir; nests this run under <prior>/reanalysis/")
    i.add_argument("--beside", action="store_true",
                   help="with --reanalysis-of: put this run BESIDE the prior session, as "
                        "<date>_<name>_v<N> in the prior's parent folder, instead of inside it -- "
                        "nothing is written into the prior (a delivered folder)")
    i.add_argument("--submission", action="append", metavar="PROT_NUMBER",
                   help="Core staff: the CoreOmics submission this analysis is for (807 or "
                        "PROT_0807); repeat it for an analysis that combines submissions, the "
                        "main one first. Kept in session.json's analysis_for")
    i.add_argument("--internal", metavar="TITLE",
                   help="Core staff: internal work (a method test, QC, training, a tool "
                        "comparison), with the title it has on the Core Project Board")
    i.set_defaults(func=do_init)
    fp = sub.add_parser("find-prior", help="find existing sessions covering the same raw files "
                                           "(to mention; a new session is fresh by default)")
    fp.add_argument("--base", default="", help="also scan this central location's sessions/ (optional)")
    fp.add_argument("--raw", nargs="*", required=True)
    fp.set_defaults(func=do_find_prior)
    f = sub.add_parser("finalize", help="write README, tidy output/, optionally zip")
    f.add_argument("--dir", required=True, help="the session directory")
    f.add_argument("--reanalysis-of", default="", help="prior session dir to diff against (auto-detected if omitted)")
    f.add_argument("--zip", action="store_true", help="also produce <session>.zip")
    f.add_argument("--call-id", metavar="TOKEN",
                   help="a fresh token (8+ letters and digits) written into this command line, "
                        "for save_transcript.py: it tells which conversation -- the main one or "
                        "a subagent's -- is running it. Needed while subagents are at work")
    f.add_argument("--no-deposit", action="store_true",
                   help="skip the repository-deposit package (output/DATA_SUBMISSION); the "
                        "publication Methods are still ensured")
    f.add_argument("--no-notify", action="store_true",
                   help="do not post 'analysis complete' to the Core's Slack channel (same as "
                        "SKILL_SLACK=0; see references/notifications.md)")
    f.add_argument("--override-finish", metavar="REASON", default=None,
                   help="go ahead although the skill's own report, the reproducibility bundle or "
                        "(a Claude Code analysis) the saved conversation is missing -- or "
                        "AGENTS.md, README.md or the Methods could not be written; the reason is "
                        "recorded in MANIFEST.txt and output/FINALIZED.json")
    f.add_argument("--no-transcript", metavar="REASON", default=None,
                   help="a Claude Code analysis whose conversation cannot be saved: finalize "
                        "without it; MANIFEST.txt records [SKIPPED] Analysis conversation -- "
                        "--no-transcript: REASON")
    # work done outside the session, folded in (custom_work.py): output/custom/<kind>/<name>
    f.add_argument("--report", metavar="HTML", default=None,
                   help="a report made outside the skill's steps, copied to output/custom/report/"
                        " (one the skill did not write goes in only with --override-finish)")
    f.add_argument("--methods", metavar="MD", default=None,
                   help="Methods written outside the skill's steps, copied to "
                        "output/custom/methods/ (the skill's output/methods.md is still made)")
    f.add_argument("--code-dir", metavar="DIR", action="append", default=[],
                   help="analysis code written outside the session, copied to "
                        "output/custom/code/<name>/ (repeatable)")
    f.add_argument("--extra-dir", metavar="DIR", action="append", default=[],
                   help="other work done outside the session, copied to "
                        "output/custom/extra/<name>/ (repeatable)")
    f.add_argument("--custom-max-gb", metavar="GB", type=float, default=None,
                   help="the most --report/--methods/--code-dir/--extra-dir may hold in all "
                        "(default 10)")
    f.set_defaults(func=do_finalize)
    st = sub.add_parser("status", help="the finish checklist: report, Methods, bundle, AGENTS.md, "
                                       "README, conversation, finalize stamp, delivery -- each "
                                       "done or missing with the command that makes it (exit 0 "
                                       "complete, 3 not)")
    # --dir too, as finalize, docs and which-search take it (C2 docs M2)
    st.add_argument("--session", "--dir", dest="session", required=True,
                    help="the session directory (--dir works too)")
    st.add_argument("--json", action="store_true", help="the checklist as JSON")
    st.set_defaults(func=do_status)
    d = sub.add_parser("docs", help="(re)write README.md, README.html and AGENTS.md -- and "
                                    "input/raw_files.txt if missing -- without finalizing")
    d.add_argument("--dir", required=True, help="the session directory")
    d.add_argument("--as", dest="as_path", default="",
                   help="where the session really is, when --dir is a copy of it (the README's "
                        "'Where this lives on HIVE' table is then about that location)")
    d.set_defaults(func=do_docs)
    w = sub.add_parser("which-search", help="the session's search (the one the DE read, wherever "
                                            "it is): how it was chosen, or why none was; exit 3 "
                                            "when none is chosen")
    w.add_argument("--dir", required=True, help="the session directory")
    w.add_argument("--paths", action="store_true",
                   help="only the files a working copy pulls (the search record, the report the "
                        "DE read -- or, for a search for identifications only, "
                        "identifications.tsv and its record -- the workflow manifest), one per "
                        "line, session-relative")
    w.set_defaults(func=do_which_search)
    a = ap.parse_args()
    if a.cmd == "finalize":
        # The finish gate before anything else -- on a login node, before the re-exec below, so a
        # refusal costs no compute job. Its warnings and its try at saving the conversation are
        # said once: here on the first run; the compute-node re-exec inherits FINISH_GATE_ENV
        # and runs it quietly (a refusal there is still said), and so does do_finalize.
        p = paths_for(a.dir)
        gate = {}
        if os.path.isdir(p["session_dir"]):
            gate = finalize_gate(a, p, speak=os.environ.get(FINISH_GATE_ENV) != p["session_dir"])
            os.environ[FINISH_GATE_ENV] = p["session_dir"]
        # Methods, the Word copy, the PDF and a zip of the whole session: a compute node's work.
        # On a cluster login node it runs itself again there (on_compute.py; HPC@UCD's ucd-hpc rules).
        # Building the bundle too (provenance.py hashes every input file: its own 60 minutes)
        # doubles the time asked for.
        try:
            import on_compute
            on_compute.ensure_compute_node("session.py finalize builds the documents and zips the "
                                           "session", minutes=120 if gate.get("build_bundle")
                                           else 60, mem="8G", script=__file__)
        except ImportError:
            pass                          # a partial scripts/ without on_compute.py: as before
    a.func(a)


if __name__ == "__main__":
    main()
