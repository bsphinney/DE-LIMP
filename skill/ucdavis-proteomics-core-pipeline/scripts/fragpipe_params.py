#!/usr/bin/env python3
"""
fragpipe_params.py -- THE reader of what a FragPipe search ran with (2.14.2, L11 N4).

A FragPipe search folder (--workdir) holds two records of its settings, both written by FragPipe:

  * fragger.params   -- the MSFragger parameters it ran with: `key = value  # comment`, one per
                        line; a line starting with '#' is a comment or a switched-off entry
                        (FragPipe 24.0 writes a disabled variable mod as "# variable_mod_06 = ...").
  * fragpipe.workflow -- the workflow it saved there: Java properties (`key=value`, '#' comment
                        lines, `\\:` and `\\=` escaped), headed "# <Tool> version <v>". Its
                        `phi-report.filter` is the Philosopher filter that set the FDR.

make_methods.py (the Methods and the SDRF, through search_record) and run_search.py (the quantity
declaration) both read them here. Stdlib only.
"""
import glob
import os
import re
import shlex

FRAGGER_PARAMS = "fragger.params"            # FragPipe writes it into the search folder
SAVED_WORKFLOW = "fragpipe.workflow"         # ...and the workflow it ran
PHILOSOPHER_FILTER_KEY = "phi-report.filter"

# Philosopher's filter defaults (Nesvilab/philosopher cmd/filter.go: --psm, --pep, --ion and
# --prot all 0.01; --sequential, --picked, --razor off), read 2026-11-03. A threshold the filter
# line does not set is one of these, and is said to be.
PHILOSOPHER_DEFAULTS = {"psm": 0.01, "pep": 0.01, "ion": 0.01, "prot": 0.01}
PHILOSOPHER_SWITCHES = ("sequential", "picked", "razor", "2d")

# MSFragger's mass units (fragger.params comments: "0 for Da, 1 for ppm")
UNITS = {"0": "Da", "1": "ppm"}


def read_fragger_params(path):
    """fragger.params -> {key: value} (values stripped, inline '# ...' comments removed); a
    commented-out line is not a setting. Raises OSError when it cannot be read."""
    out = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or "=" not in line:
                continue
            key, _, val = line.partition("=")
            out[key.strip()] = val.split("#", 1)[0].strip()
    return out


def _unescape(v):
    return re.sub(r"\\(.)", r"\1", v)


def read_workflow(path):
    """A .workflow -> ({key: value}, {tool: version} from its '# <Tool> version <v>' header).
    Raises OSError when it cannot be read."""
    props, versions = {}, {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if line.startswith("#"):
                m = re.match(r"#\s*(\S+(?:\s\S+)*?) version (\S.*)$", line)
                if m:
                    versions[m.group(1)] = m.group(2).strip()
                continue
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            props[_unescape(key.strip())] = _unescape(val.strip())
    return props, versions


def _words(text):
    """A workflow's filter line as words: shell-style when it parses, else on whitespace (an
    unpaired quote must not stop the Methods; follow-up 2)."""
    try:
        return shlex.split(text or "")
    except ValueError:
        return (text or "").split()


def philosopher_filter(text):
    """A Philosopher filter line (`--sequential --prot 1`) -> {"raw", "psm", "pep", "ion", "prot",
    "set": the thresholds the line sets, the switches as booleans, "other": anything else}."""
    words = _words(text)
    rec = {"raw": (text or "").strip(), "set": set(), "other": []}
    rec.update(PHILOSOPHER_DEFAULTS)
    rec.update({s: False for s in PHILOSOPHER_SWITCHES})
    i = 0
    while i < len(words):
        w = words[i].lstrip("-")
        if w in PHILOSOPHER_DEFAULTS and i + 1 < len(words):
            try:
                rec[w] = float(words[i + 1])
                rec["set"].add(w)
                i += 2
                continue
            except ValueError:
                pass
        if w in PHILOSOPHER_SWITCHES:
            rec[w] = True
        else:
            rec["other"].append(words[i])
        i += 1
    return rec


def _pct(x):
    return f"{x * 100:g}%"


# Ruling 120 condition 1: a filter with no protein-level FDR (--prot 1, as FragPipe 24.0's
# Nonspecific-peptidome and Nonspecific-HLA workflows ship) is KEPT and disclosed, never "fixed".
# Philosopher's protein FDR is not a filter on the protein list alone: with a protein inference
# both its --sequential and its 2D filter keep only the PSMs of proteins that pass it
# (Nesvilab/philosopher lib/fil/fil.go: Run, extractPSMfromPepXML -- read 2026-11-11), so
# --prot 0.01 would drop every endogenous peptide of a protein seen by one or two peptides, which
# is most of a peptidome. The vendor's workflows leave it off for that reason; so does the skill.
PROTEIN_FDR_OFF_WHY = (
    "Philosopher's protein-level FDR also removes the PSMs of every protein that fails it (with a "
    "protein inference both its --sequential and its 2D filter keep only the PSMs of passing "
    "proteins: philosopher lib/fil/fil.go, extractPSMfromPepXML), so in a non-specific "
    "(peptidome) search it would remove the endogenous peptides of the many proteins seen by one "
    "or two peptides. FragPipe's Nonspecific workflows ship --prot 1 for that reason, and the "
    "skill keeps it; the peptide-level FDR is the result")
PROTEIN_GROUPS_UNCONTROLLED = ("the proteins the peptides map to are therefore not FDR-controlled, "
                               "and no protein count is reported as a result")


def protein_fdr_record(text):
    """Whether a Philosopher filter line applies a protein-level FDR -- THE reading, for the
    preset's record (make_presets: the manifest), fragpipe_identifications.json and the analysis
    brief. -> {"applied": True | False | None (no filter recorded), "threshold", "filter", "why",
    "protein_filter": the sentence a reader is given}."""
    if not (text or "").strip():
        return {"applied": None, "threshold": None, "filter": None,
                "why": "the Philosopher filter this search ran is not recorded",
                "protein_filter": "not recorded (the Philosopher filter is not known): give no "
                                  "protein-level FDR and no protein count as a result"}
    f = philosopher_filter(text)
    if f["prot"] >= 1:
        return {"applied": False, "threshold": None, "filter": f["raw"],
                "why": PROTEIN_FDR_OFF_WHY,
                "protein_filter": (f"none (the Philosopher filter {f['raw']} sets --prot "
                                   f"{f['prot']:g}): {PROTEIN_GROUPS_UNCONTROLLED}")}
    how = "picked protein FDR" if f["picked"] else "protein FDR"
    dflt = "" if "prot" in f["set"] else " (Philosopher's default; not set in the filter)"
    return {"applied": True, "threshold": f["prot"], "filter": f["raw"],
            "why": f"the filter sets a protein-level FDR ({how}){dflt}",
            "protein_filter": f"protein {_pct(f['prot'])} FDR ({how}, Philosopher filter "
                              f"{f['raw']}){dflt}"}


def no_protein_fdr(pfdr):
    """THE predicate on a protein_fdr_record(): True when the record says the protein groups rest
    on no protein-level FDR (applied False, or None: the filter is not recorded); False when it
    says one was applied, or there is no record at all (an older search: nothing is claimed
    either way). The adapters, the analysis brief and the report's callouts all ask this, never
    a test of their own (C2 proteomics H1)."""
    return isinstance(pfdr, dict) and "applied" in pfdr and pfdr["applied"] is not True


def recorded_protein_fdr(*records):
    """The protein_fdr_record() a search's own record carries: the first of `records`
    (fragpipe_identifications.json for an identifications-only search, search_provenance.json's
    `fragpipe_quantity` for a quantified one) whose `protein_fdr` is a dict, else None."""
    for rec in records:
        p = rec.get("protein_fdr") if isinstance(rec, dict) else None
        if isinstance(p, dict):
            return p
    return None


def philosopher_methods_sentence(text):
    """The Methods sentence for a Philosopher filter line, every threshold as the line set it or
    as Philosopher's default (said so)."""
    f = philosopher_filter(text)
    levels = (("psm", "PSM"), ("ion", "peptide ion"), ("pep", "peptide"))
    given = [f"{label} {_pct(f[k])}" for k, label in levels if k in f["set"]]
    default = [f"{label} {_pct(f[k])}" for k, label in levels if k not in f["set"]]

    def join(xs):
        return xs[0] if len(xs) == 1 else ", ".join(xs[:-1]) + " and " + xs[-1]
    parts = []
    if given:
        parts.append(join(given) + " FDR")
    if default:
        parts.append(join(default) + " FDR (Philosopher's defaults; not set in the filter)")
    how = []
    if f["picked"]:
        how.append("picked protein FDR")
    if f["razor"]:
        how.append("razor peptides for protein scoring")
    if f["prot"] >= 1:
        parts.append(f"no protein-level FDR filter (--prot {f['prot']:g}): "
                     f"{PROTEIN_GROUPS_UNCONTROLLED}")
    else:
        parts.append(f"protein {_pct(f['prot'])} FDR" + (f" ({'; '.join(how)})" if how else "")
                     + ("" if "prot" in f["set"] else " (Philosopher's default)"))
    if f["sequential"]:
        # Philosopher's own description of the flag (cmd/filter.go)
        parts.append("with Philosopher's sequential algorithm (--sequential), which estimates "
                     "FDR using both the filtered PSM and protein lists")
    extra = f"; other filter options: {' '.join(f['other'])}" if f["other"] else ""
    return (f"Identifications were filtered with Philosopher (filter {f['raw']}): "
            + "; ".join(parts) + extra + ".")


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def digest(params):
    """fragger.params' digest -> {"specificity": "none" | "semi" | "full" | None, "termini",
    "enzymes": [{"name", "cut", "nocut", "sense", "missed"}]}. MSFragger's num_enzyme_termini:
    0 non-enzymatic, 1 semi-enzymatic, 2 fully enzymatic (the file's own comment); a cleavage
    named `nonspecific` with 0 termini is FragPipe's non-specific search."""
    enzymes = []
    for n in (1, 2):
        name = params.get(f"search_enzyme_name_{n}")
        if name is None or name.lower() in ("", "null"):
            continue
        enzymes.append({"name": name, "cut": params.get(f"search_enzyme_cut_{n}", ""),
                        "nocut": params.get(f"search_enzyme_nocut_{n}", ""),
                        "sense": params.get(f"search_enzyme_sense_{n}", ""),
                        "missed": _num(params.get(f"allowed_missed_cleavage_{n}"))})
    t = params.get("num_enzyme_termini")
    termini = int(_num(t)) if _num(t) is not None else None
    spec = {0: "none", 1: "semi", 2: "full"}.get(termini)
    return {"specificity": spec, "termini": termini, "enzymes": enzymes}


# MSFragger site notation (fragger.params: "* for any amino acid, [ and ] specifies protein
# termini, n and c specifies peptide termini"; '^' after a terminus = any residue there)
SITE_POSITION = {"n": "Any N-term", "c": "Any C-term", "[": "Protein N-term",
                 "]": "Protein C-term"}
# the add_* keys of a fixed modification at a terminus
FIXED_TERMINI = {"Nterm_peptide": "Any N-term", "Cterm_peptide": "Any C-term",
                 "Nterm_protein": "Protein N-term", "Cterm_protein": "Protein C-term"}


def mods(params):
    """fragger.params' modifications -> [(type, mass, position, targets, max per peptide)]:
    fixed = every add_* that is not 0, variable = every variable_mod_NN that is switched on."""
    out = []
    for key, val in params.items():
        # a residue (add_C_cysteine) or a terminus (add_Nterm_peptide): MSFragger has other add_*
        # keys that are not modifications (add_topN_complementary; review F4)
        m = re.match(r"add_([A-Z]_\w+|[NC]term_(?:peptide|protein))$", key)
        mass = _num(val)
        if not m or not mass:
            continue
        sub = m.group(1)
        if sub in FIXED_TERMINI:
            out.append(("fixed", mass, FIXED_TERMINI[sub], "", None))
        else:
            out.append(("fixed", mass, "Anywhere", sub.split("_", 1)[0], None))
    for key in sorted(k for k in params if re.match(r"variable_mod_\d+$", k)):
        bits = params[key].split()
        mass = _num(bits[0]) if bits else None
        if mass is None:
            continue
        site = bits[1] if len(bits) > 1 else ""
        per = int(_num(bits[2])) if len(bits) > 2 and _num(bits[2]) is not None else None
        for pos, targets in site_parts(site):
            out.append(("variable", mass, pos, targets, per))
    return out


def site_parts(site):
    """An MSFragger site string -> [(position, residues)], one per position it names, in order:
    each terminus marker (n, c, [, ]) starts a part that runs to the next marker, residues before
    any marker are "Anywhere"; '^' (any residue at that terminus) and '*' (any amino acid) add no
    residue. "nQnC" is the N-terminal Q or C ("Any N-term", "QC"), never one target "QnC"
    (2.14.2 C2 proteomics L3); "STY" is ("Anywhere", "STY")."""
    parts, pos, res = [], "Anywhere", ""
    for ch in site:
        if ch in SITE_POSITION:
            if pos != "Anywhere" or res:
                parts.append((pos, res))
            pos, res = SITE_POSITION[ch], ""
        elif ch not in "^*":
            res += ch
    if pos != "Anywhere" or res or not parts:
        parts.append((pos, res))
    merged = {}                                  # one part per position, residues in order
    for p_, r in parts:
        merged[p_] = merged.get(p_, "") + "".join(c for c in r if c not in merged.get(p_, ""))
    return list(merged.items())


def search_extras(params):
    """What makes an MSFragger search more than a closed search, as fragger.params sets it: mass
    offsets, a labile search mode, mass differences turned into variable mods, a delta-mass
    residue restriction (review F3) -> [phrase], [] for a plain search."""
    out = []
    offs = [x for x in re.split(r"[/\s]+", params.get("mass_offsets", "") or "") if x]
    if any(_num(x) for x in offs):
        vals = [x for x in offs]
        out.append("mass offsets of " + (", ".join(vals[:-1]) + " and " + vals[-1]
                                          if len(vals) > 1 else vals[0]) + " Da (mass_offsets)")
    if params.get("use_detailed_offsets") == "1":
        out.append("the detailed mass-offset list (mass_offsets_detailed)")
    mode = (params.get("labile_search_mode") or "off").strip()
    if mode.lower() != "off":
        out.append(f"labile search mode {mode} (labile_search_mode)")
    md = (params.get("mass_diff_to_variable_mod") or "0").strip()
    if md not in ("0", ""):
        out.append(f"unassigned mass differences as variable modifications "
                   f"(mass_diff_to_variable_mod={md})")
    rd = (params.get("restrict_deltamass_to") or "all").strip()
    if rd.lower() != "all" and out:
        out.append(f"mass offsets restricted to {rd} (restrict_deltamass_to)")
    return out


# FragPipe's console logs: log_<date>.txt, written into the search's workdir. The ONE list of
# them every reader walks (the Philosopher filter, IonQuant's MaxLFQ line, MSBooster's bandwidth):
# newest first -- the date in the name sorts as text.
LOG_GLOB = "log_*.txt"


def logs(out):
    """FragPipe's console logs in the search folder `out`, NEWEST FIRST; [] when there are none."""
    return sorted(glob.glob(os.path.join(out, LOG_GLOB)), reverse=True) if out else []


# FragPipe's console log (log_<date>.txt in the workdir): each Philosopher filter it ran, as
# "PhilosopherFilter [Work dir: <dir>]" and then the command line (FragPipe 24.0). The filter AS
# RUN can differ from the workflow's phi-report.filter: FragPipe added --razor and --tag in both
# runs read (a 38-run hydrolysate search; H5's LFQ-MBR).
_FILTER_DIR = re.compile(r"^PhilosopherFilter \[Work dir: (.*)\]\s*$")
_FILTER_CMD = re.compile(r"philosopher\S*\s+filter\s+(.*)$")


def read_filter_runs(log):
    """FragPipe's log -> [{"args", "workdir"}], one per Philosopher filter run. Raises OSError."""
    runs, wd = [], None
    with open(log, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = _FILTER_DIR.match(line.strip())
            if m:
                wd = m.group(1).strip()
                continue
            m = _FILTER_CMD.search(line)
            if m and wd is not None:
                runs.append({"args": m.group(1).strip(), "workdir": wd})
                wd = None
    return runs


def filter_as_run(args):
    """A filter command line's options without its inputs and decoy tag (--pepxml, --protxml,
    --dbbin, --tag and their values): the filter in a Methods section's terms. FragPipe logs its
    argv joined with spaces and UNQUOTED, so a path may hold spaces or quotes (a work dir like
    .../it's here/ once crashed shlex, or was silently mis-paired): the line is split on
    whitespace, and a dropped option's value runs to the next option (follow-up 2)."""
    words, out, i = (args or "").split(), [], 0
    while i < len(words):
        if words[i] in ("--pepxml", "--protxml", "--dbbin", "--tag"):
            i += 1
            while i < len(words) and not words[i].startswith("-"):
                i += 1
            continue
        out.append(words[i])
        i += 1
    return " ".join(out)


def filter_scope_sentence(runs, n_runs=None):
    """Whether Philosopher's rates hold per FragPipe experiment or over all runs together, from
    the filter runs FragPipe's log shows (review F9; read, never assumed: with LFQ off FragPipe
    filtered all 38 runs of a search together although experiments were set). None = unknown."""
    dirs = sorted({r["workdir"] for r in runs})
    if not dirs:
        return None
    if len(dirs) == 1:
        return "Philosopher filtered all runs together, so these rates hold over all runs together."
    per_run = n_runs is not None and n_runs == len(dirs)
    per = "run" if per_run else "experiment"
    head = (f"Philosopher filtered each FragPipe experiment separately ({len(dirs)} experiments"
            + (", one per run" if per_run else "") + ")")
    # every experiment's filter reads the ONE combined protein inference (--protxml
    # combined.prot.xml): a protein FDR it applies is over all runs together (H7 bug 2, its log:
    # the same protein counts in each experiment's filter), only PSM/ion/peptide are per run
    args = runs[0].get("args") or ""
    prot = philosopher_filter(filter_as_run(args))["prot"]
    if "--protxml" in args.split() and prot < 1:
        return (head + f": the PSM, peptide-ion and peptide rates hold within each {per}, and "
                f"the protein-level {_pct(prot)} FDR was estimated over all runs together (every "
                f"experiment's filter read the one combined protein inference).")
    return head + f", so these rates hold within each {per}, not over the runs combined."


def find_filter_runs(out):
    """The Philosopher filter runs FragPipe's newest log in `out` (log_*.txt) shows
    (read_filter_runs), [] when none is read."""
    for log in logs(out):
        try:
            runs = read_filter_runs(log)
        except OSError:
            continue
        if runs:
            return runs
    return []


def filter_as_recorded(out, workflow_filter=None):
    """The ONE reading of the Philosopher filter a FragPipe search in `out` ran, for the Methods'
    FDR sentence and the identifications record alike (follow-up 3, H7: the record quoted the
    workflow's line, --sequential --prot 1, while the Methods quoted the log's, ... --razor).
    -> {"filter", "from_log", "runs", "log_error"}: the filter AS RUN (filter_as_run of the first
    filter run in FragPipe's newest log, find_filter_runs) when its line is read and parses; else
    `workflow_filter`, the workflow's phi-report.filter, which FragPipe may add options to when it
    runs it (from_log False; log_error says why a log line found was not used). "filter" is None
    when neither is there."""
    runs, err = find_filter_runs(out), None
    if runs:
        try:
            text = filter_as_run(runs[0]["args"])
            philosopher_filter(text)                  # it must parse, or the workflow's line
            return {"filter": text, "from_log": True, "runs": runs, "log_error": None}
        except (ValueError, IndexError, TypeError, KeyError) as e:
            err = f"{type(e).__name__}: {e}"
    return {"filter": workflow_filter, "from_log": False, "runs": [], "log_error": err}


# MSBooster (FragPipe 24.0 runs v1.4.14) fits each run's retention times to the predicted ones
# with a LOESS whose bandwidth it CHOOSES for that run by a grid search; FragPipe's log says so per
# run: "Processing pin <...>/<experiment>/<run>.pin", "RT regression using 126 PSMs", then "Best
# average bandwidth for mass  from grid search of 0.01,0.05,0.1,0.2 after 5 iterations is 0.1"
# (the same words follow each IM regression after it, which is not the RT one). It is no setting
# of the skill or the workflow, and two searches of the same input can choose differently: the
# two-run timsTOF DDA reference gave 0.1 once and 0.11 another time for one run, and with it
# 6,797 PSMs / 187 protein groups or 6,762 / 186 (everything before MSBooster byte-identical;
# Ruling 106). So the value is recorded as MSBooster's choice, from the log, or as not stated.
MSBOOSTER_CHOSEN = "chosen by MSBooster"
# its key in search_provenance.json, fragpipe_quantity.json and fragpipe_identifications.json
MSBOOSTER_KEY = "msbooster_rt_bandwidth"
NOT_STATED = "not stated in the log"
_MSB_PIN = re.compile(r"\bProcessing pin (.+\.pin)\s*$")
_MSB_RT = re.compile(r"\bRT regression using (\d+) PSMs")
_MSB_IM = re.compile(r"\bIM regression\b")
_MSB_BW = re.compile(r"\bBest average bandwidth\b.*?(?:grid search of (\S+).*?)?\bis\s+"
                     r"([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)\s*$")


def msbooster_rt_bandwidths(out):
    """MSBooster's retention-time LOESS bandwidth for each run, as FragPipe's newest log in `out`
    states it (logs()[0]: the search that left these files). -> {"what", "source", "log" (its file
    name), "status": "stated" | NOT_STATED, "runs": [{"run" (the .pin's name), "experiment" (its
    folder), "bandwidth": the number MSBooster chose | NOT_STATED, "source": MSBOOSTER_CHOSEN |
    NOT_STATED, "rt_regression_psms", "grid"}]}. A run's bandwidth is the line after ITS "RT
    regression", never an IM regression's. File names only, never a path."""
    found = logs(out)
    log = found[0] if found else None
    rec = {"what": "MSBooster's retention-time LOESS bandwidth, per run",
           "source": (f"{MSBOOSTER_CHOSEN}: it picks the bandwidth for each run by a grid search "
                      f"-- not a setting of the skill or the workflow; read from FragPipe's log"
                      + (f" ({os.path.basename(log)})" if log else "")),
           "log": os.path.basename(log) if log else None, "status": NOT_STATED, "runs": []}
    if not log:
        rec["source"] = f"{NOT_STATED}: no FragPipe log (log_*.txt) in the search folder"
        return rec
    runs, cur, state = [], None, None
    try:
        with open(log, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\r\n")
                m = _MSB_PIN.search(line)
                if m:
                    pin = m.group(1).replace("\\", "/")
                    cur = {"run": os.path.basename(pin)[:-len(".pin")],
                           "experiment": os.path.basename(os.path.dirname(pin)) or None,
                           "bandwidth": NOT_STATED, "source": NOT_STATED,
                           "rt_regression_psms": None, "grid": None}
                    runs.append(cur)
                    state = None
                    continue
                if cur is None:
                    continue
                m = _MSB_RT.search(line)
                if m:
                    cur["rt_regression_psms"], state = int(m.group(1)), "rt"
                    continue
                if _MSB_IM.search(line):
                    state = "im"
                    continue
                m = _MSB_BW.search(line)
                if m and state == "rt" and cur["source"] == NOT_STATED:
                    cur.update(bandwidth=float(m.group(2)), source=MSBOOSTER_CHOSEN,
                               grid=m.group(1))
                    state = None
    except OSError as e:
        rec["source"] = f"{NOT_STATED}: FragPipe's log could not be read ({type(e).__name__})"
        return rec
    rec["runs"] = runs
    if any(r["source"] == MSBOOSTER_CHOSEN for r in runs):
        rec["status"] = "stated"
    elif not runs:
        rec["source"] = (f"{NOT_STATED}: FragPipe's log ({os.path.basename(log)}) shows no "
                         f"MSBooster run")
    return rec


def msbooster_rerun_sentence(bw):
    """THE sentence (Ruling 106) for the Methods (make_methods.py) and REPRODUCE.md
    (provenance.py): a re-run can differ slightly because MSBooster chose its RT bandwidth per run,
    with the values this search ran with, from the record msbooster_rt_bandwidths() wrote. None
    when the record shows no MSBooster run. It names no run file (Ruling 110.2: file names can
    carry sample or lab names, and Methods travel into papers): the range across the runs, how
    many runs the log does not state, and where each run's value is (the search record)."""
    if not isinstance(bw, dict) or not bw.get("runs"):
        return None
    runs = bw["runs"]
    vals = [r["bandwidth"] for r in runs if r.get("source") == MSBOOSTER_CHOSEN
            and isinstance(r.get("bandwidth"), (int, float))]
    n_missing = len(runs) - len(vals)
    n = len(runs)
    head = ("A re-run of this FragPipe search can differ slightly from it, because MSBooster "
            "chooses its retention-time smoothing (the LOESS bandwidth) for each run itself; ")
    where = f"each run's value is in the search record, {MSBOOSTER_KEY} in search_provenance.json"
    if not vals:
        return (head + "FragPipe's log does not state the value for "
                + ("its run" if n == 1 else f"any of its {n} runs") + f" ({where}).")
    lo, hi = min(vals), max(vals)
    if len(vals) == 1:
        given = f"{lo:g} for one run"
    elif lo == hi:
        given = f"{lo:g} for {'both runs' if len(vals) == 2 else f'each of the {len(vals)} runs'}"
    else:
        given = f"{lo:g}–{hi:g} across the {len(vals)} runs"
    if n_missing:
        given += f"; the log does not state it for {n_missing} of the {n} runs"
    return head + f"this search ran with {given} ({MSBOOSTER_CHOSEN}, as FragPipe's log states; {where})."


# The steps FragPipe ran, as its console log closes each one: "Process '<step>' finished, exit
# code: <n>" (FragPipe 24.0; the command list at the top of the log is only the plan). The Methods
# name and cite a tool only when it ran (make_methods.py), with the version the log's header
# states ("<Tool> version <v>").
_PROCESS_DONE = re.compile(r"Process '([^']+)' finished, exit code: (-?\d+)")
_LOG_VERSION = re.compile(r"^([A-Za-z][\w .-]*?) version (\S.*?)\s*$")


def tools_run(out):
    """What FragPipe's newest log in `out` says ran: {"log" (its file name), "ran" {step: times it
    finished with exit code 0}, "failed" {step: times it finished otherwise}, "versions" {tool:
    version as the log's header states it; "N/A" left out}}; None with no log or an unreadable
    one. File names only, never a path."""
    found = logs(out)
    if not found:
        return None
    rec = {"log": os.path.basename(found[0]), "ran": {}, "failed": {}, "versions": {}}
    try:
        with open(found[0], encoding="utf-8", errors="replace") as fh:
            head = True
            for line in fh:
                line = line.rstrip("\r\n")
                if head:
                    m = _LOG_VERSION.match(line.strip())
                    if m and m.group(2) != "N/A":
                        rec["versions"].setdefault(m.group(1).strip(), m.group(2))
                    if line.startswith("LCMS files"):
                        head = False
                m = _PROCESS_DONE.search(line)
                if m:
                    key = "ran" if m.group(2) == "0" else "failed"
                    rec[key][m.group(1)] = rec[key].get(m.group(1), 0) + 1
    except OSError:
        return None
    return rec


def manifest_rows(out):
    """The rows of the search's own fragpipe.fp-manifest (one per run), [] when unreadable."""
    try:
        with open(os.path.join(out or "", "fragpipe.fp-manifest"), errors="replace") as fh:
            return [ln for ln in fh if ln.strip()]
    except OSError:
        return []


WORKFLOW_LINE_TAG = ("[This is the workflow's filter line{why}; FragPipe may add options when it "
                     "runs it, and whether its rates hold per run or over all runs together is in "
                     "its log (log_*.txt), not read here — confirm.]")


def fdr_statement(out, workflow_filter):
    """THE FDR statement of a FragPipe search in `out`, for the Methods' search paragraph and the
    DE's identification_fdr alike (2.14.2 C2 proteomics M1: the DE quoted the workflow's line
    while the Methods quoted the log's, and dropped the scope): the filter AS RUN
    (filter_as_recorded) in words (philosopher_methods_sentence) with over which runs its rates
    hold (filter_scope_sentence), else the workflow's line, tagged WORKFLOW_LINE_TAG with why.
    -> {"filter", "from_log", "sentence", "scope"} or None with no filter at all."""
    fr = filter_as_recorded(out, workflow_filter)
    log_err = fr["log_error"]
    if fr["filter"] and fr["from_log"]:
        try:
            scope = filter_scope_sentence(fr["runs"], len(manifest_rows(out)) or None)
            return {"filter": fr["filter"], "from_log": True, "scope": scope,
                    "sentence": philosopher_methods_sentence(fr["filter"])
                    + (" " + scope if scope else "")}
        except (ValueError, IndexError, TypeError) as e:
            # never a crash of finalize: the workflow's line, tagged, with why (FU2)
            log_err = f"{type(e).__name__}: {e}"
    if not workflow_filter:
        return None
    return {"filter": workflow_filter, "from_log": False, "scope": None,
            "sentence": philosopher_methods_sentence(workflow_filter) + " " + WORKFLOW_LINE_TAG.format(
                why=(f": FragPipe's log line for it could not be read ({log_err})"
                     if log_err else ""))}


def tolerance(params, which):
    """("precursor" | "fragment") -> {"value", "unit", "symmetric", "lo", "hi"} or None."""
    if which == "precursor":
        lo, hi = _num(params.get("precursor_mass_lower")), _num(params.get("precursor_mass_upper"))
        unit = UNITS.get(params.get("precursor_mass_units", ""))
        if lo is None or hi is None or unit is None:
            return None
        sym = abs(lo) == hi
        return {"value": hi if sym else f"{lo:g} to +{hi:g}", "unit": unit, "symmetric": sym,
                "lo": lo, "hi": hi}
    tol = _num(params.get("fragment_mass_tolerance"))
    unit = UNITS.get(params.get("fragment_mass_units", ""))
    if tol is None or unit is None:
        return None
    return {"value": tol, "unit": unit, "symmetric": True, "lo": -tol, "hi": tol}


# MSFragger's calibrate_mass, in fragger.params' own comment: "0 for OFF, 1 for ON, 2 for ON and
# find optimal parameters, 4 for ON and find the optimal fragment mass tolerance". FragPipe 24.0's
# LFQ-MBR (the timsTOF DDA default) ships 2, its Nonspecific workflows 0.
CALIBRATE_MASS = {
    1: "recalibrated the precursor and fragment masses",
    2: "recalibrated the precursor and fragment masses and chose its search parameters from the "
       "data",
    4: "recalibrated the precursor and fragment masses and chose the fragment mass tolerance from "
       "the data",
}


def calibration_sentence(value, source):
    """The Methods sentence for MSFragger's calibrate_mass as the search's workflow set it
    (`source`: the file it was read from) -- a statement of what ran, never a warning (C2
    proteomics M5). "" when off (0) or not recorded; an unknown value is named as such."""
    v = _num(value)
    if not v:
        return ""
    what = CALIBRATE_MASS.get(int(v))
    head = (f"MSFragger {what} (calibrate_mass={int(v)}, as the search's workflow sets it; "
            if what else
            f"MSFragger ran with calibrate_mass={int(v)}, a value not described here ("
            f"the search's workflow sets it; ")
    if int(v) == 1:
        return head + f"{source})."
    return (head + f"{source}), so the tolerances given here are its starting values; the values "
            f"it searched with are in FragPipe's log.")


# FragPipe 24.0's quantification switches, as its shipped workflows set them (read 2026-11-03):
# LFQ-MBR runs label-free quantification; TMT10 TMT-Integrator; DIA_SpecLib_Quant DIA-NN;
# Nonspecific-peptidome and Nonspecific-HLA none of them. ionquant.run-ionquant is true in ALL of
# them: it picks IonQuant over FreeQuant inside label-free quantification, it does not switch
# quantification on. THE reading of whether a FragPipe search quantifies (L11, Ruling 73a):
# run_search.py's job check, its identifications-only route and record, and make_methods.py.
LFQ_KEY = "quantitation.run-label-free-quant"
QUANT_SWITCHES = ((LFQ_KEY, "label-free quantification (IonQuant)"),
                  ("tmtintegrator.run-tmtintegrator", "isobaric quantification (TMT-Integrator)"),
                  ("diann.run-dia-nn", "DIA-NN quantification"),
                  ("diann.run-dia-plex", "plexDIA quantification"))
# What a finished FragPipe search has written: each experiment's identification tables, and
# IonQuant's protein table when label-free quantification ran (references/search-engines.md)
ID_TABLES = ("psm.tsv", "peptide.tsv", "protein.tsv")
LFQ_TABLE = "combined_protein.tsv"


def _tri(v):
    return None if v is None else str(v).strip().lower() == "true"


REPORT_KEY = "phi-report.run-report"


def report_runs(props):
    """Whether Philosopher's report (psm.tsv, peptide.tsv, protein.tsv and the FDR filter it
    applies) runs in this workflow -- THE one reading, for the job's table check
    (expected_outputs) and the Methods' FDR sentence (follow-up 2, verify R3). An absent key
    counts as true: all 78 of FragPipe 24.0's shipped workflows set it true (read 2026-11-03)."""
    v = props.get(REPORT_KEY)
    return True if v is None else str(v).strip().lower() == "true"


# The Methods' words for a workflow that ran no quantification: publication text, so no config key
# (the parameter table keeps the key; 2.14.2 C2 biologist L4)
NO_QUANT_METHODS = "FragPipe's workflow ran no label-free (IonQuant), isobaric or DIA-NN quantification"
_NO_QUANT_WHY = "the workflow ran no label-free quantification ("


def methods_why(why):
    """A workflow_quantification `why` (or a record's copy of it) as the Methods say it:
    NO_QUANT_METHODS for "ran no quantification", else as written."""
    return NO_QUANT_METHODS if str(why or "").startswith(_NO_QUANT_WHY) else why


def workflow_quantification(props):
    """A workflow's properties -> {"quantified": True | False | None (not said), "label_free":
    bool, "runs": [what it quantifies with], "switches": {key: value}, "reports": {"report",
    "protein_prophet": True | False | None}, "why"}. Identifications only (False) needs the
    label-free switch itself to say false: a workflow without it is None, decided from what the
    search wrote (run_search.fragpipe_run_quantification; review F7)."""
    seen = {k: props[k] for k, _ in QUANT_SWITCHES if k in props}
    on = [label for k, label in QUANT_SWITCHES if str(props.get(k, "")).strip().lower() == "true"]
    reports = {"report": report_runs(props),
               "protein_prophet": _tri(props.get("protein-prophet.run-protein-prophet"))}
    if not on and LFQ_KEY not in props:
        return {"quantified": None, "label_free": False, "runs": [], "switches": seen,
                "reports": reports,
                "why": (f"the workflow does not say whether it ran label-free quantification "
                        f"(no {LFQ_KEY})" if seen else
                        "the workflow names none of FragPipe's quantification switches")}
    why = (f"the workflow ran {'; '.join(on)}" if on else
           f"the workflow ran no label-free quantification ({LFQ_KEY}="
           f"{props.get(LFQ_KEY, 'not set')}; no isobaric or DIA-NN quantification either)")
    return {"quantified": bool(on), "label_free": LFQ_KEY in seen and
            str(props[LFQ_KEY]).strip().lower() == "true", "runs": on, "switches": seen,
            "reports": reports, "why": why}


def workflow_file_quantification(path):
    """workflow_quantification() of a .workflow file; quantified None when it cannot be read."""
    try:
        props, _ = read_workflow(path)
    except (OSError, TypeError) as e:
        return {"quantified": None, "label_free": False, "runs": [], "switches": {},
                "why": f"{path} could not be read ({e})"}
    return workflow_quantification(props)


def expected_outputs(quant):
    """The tables a finished search of this workflow must have written (depth <= 2 under its
    workdir), as its switches say (review F6): psm.tsv and peptide.tsv when Philosopher's report
    ran (report_runs: an absent key counts as true), protein.tsv when ProteinProphet ran too (an
    absent protein-prophet key requires nothing), LFQ_TABLE when label-free quantification ran."""
    r = (quant or {}).get("reports") or {}
    out = []
    if r.get("report"):
        out += ["psm.tsv", "peptide.tsv"] + (["protein.tsv"] if r.get("protein_prophet") else [])
    if (quant or {}).get("label_free"):
        out.append(LFQ_TABLE)
    return tuple(out)


def find(folder, name):
    """<folder>/<name> when it is a file, else None."""
    p = os.path.join(folder, name) if folder else None
    return p if p and os.path.isfile(p) else None
