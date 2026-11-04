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
        parts.append(f"no protein-level FDR filter (--prot {f['prot']:g})")
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
    logs = sorted(glob.glob(os.path.join(out, "log_*.txt"))) if out else []
    for log in reversed(logs):
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
