#!/usr/bin/env python3
"""
metamorpheus_params.py -- the one reading of a MetaMorpheus search-task TOML, for
resolve_defaults.py (which ships the validated one) and run_search.py --engine metamorpheus
(which checks it, writes the copy the job runs, and records it).

No TOML library: the skill runs on Python 3.9 on HIVE (no tomllib) and a MetaMorpheus task file
is flat `Key = value` lines under [Section] headers, which is all this reads. Values are kept as
written (a quoted string keeps its quotes off; anything else is the raw text), keyed
"Section.Key" -- e.g. "CommonParameters.PrecursorMassTolerance" -> "±10.0000 PPM".

What MetaMorpheus 1.1.12 reads (verified on HIVE 2026-11-16, `dotnet CMD.dll --help` in the Core's
image): -t task TOML(s), -d database(s) (.fasta/.xml[.gz]), -s spectra (.raw/.mzML/.mgf/
ms2.msalign, or folders), -o output folder, -v none|minimal|normal, --acceptThermoLicence (agrees
to the Thermo RawFileReader licence, needed to read .raw, and records it in the app's
settings.toml). It does not read a Bruker .d.
"""
import hashlib
import os
import re

# The spectra MetaMorpheus 1.1.12 reads (its --help), lower-case suffixes.
SPECTRA_SUFFIXES = (".raw", ".mzml", ".mgf", ".msalign")
TASK_SUFFIX = ".toml"
# The key the run's copy changes: MetaMorpheus' threads per file. The validated TOML says 16; the
# job gets --threads CPUs, and HPC@UCD's rules say a job uses what it asked for, never a
# hard-coded count.
THREADS_KEY = "CommonParameters.MaxThreadsToUsePerFile"
_SECTION = re.compile(r"^\s*\[([^\]]+)\]\s*$")
_KEYVAL = re.compile(r"^\s*('[^']*'|[A-Za-z0-9_.\-]+)\s*=\s*(.*?)\s*$")


def _unquote(v):
    if len(v) >= 2 and v[0] == v[-1] and v[0] in "\"'":
        return v[1:-1]
    return v


def read_task(path):
    """{"Section.Key": value as written, unquoted} of a MetaMorpheus task TOML; top-level keys
    (TaskType) have no section prefix. Raises OSError for an unreadable file."""
    out, section = {}, ""
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.lstrip().startswith("#"):
                continue
            m = _SECTION.match(line)
            if m:
                section = m.group(1).strip()
                continue
            m = _KEYVAL.match(line)
            if m:
                key = _unquote(m.group(1))
                out[f"{section}.{key}" if section else key] = _unquote(m.group(2))
    return out


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tolerance_ppm(text):
    """'±10.0000 PPM' -> 10.0; an Absolute (Da) tolerance or anything else -> None."""
    m = re.search(r"([0-9]+(?:\.[0-9]+)?)\s*PPM", text or "", re.I)
    return float(m.group(1)) if m else None


def mods(value):
    """MetaMorpheus' ListOfModsVariable / ListOfModsFixed value -> [(category, name)]. Entries are
    separated by two tabs, category and name by one ("Common Variable\\tOxidation on M\\t\\t...")."""
    found = []
    for entry in (value or "").split("\t\t"):
        if "\t" in entry:
            cat, name = entry.split("\t", 1)
            if cat.strip() and name.strip():
                found.append((cat.strip(), name.strip()))
    return found


def custom_mods_defined(path):
    """{"<ID> on <target>"} defined in a MetaMorpheus modification file (ID / TG lines), the way a
    task names a mod: 'Gln to PyroGlu peptide N-term on Q'."""
    names, ident = set(), None
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ID "):
                ident = line[3:].strip()
            elif line.startswith("TG ") and ident:
                for tgt in line[3:].strip().rstrip(".").split(" or "):
                    names.add(f"{ident} on {tgt.strip()}")
            elif line.startswith("//"):
                ident = None
    return names


def summary(task):
    """The settings a record and a Methods section quote, from read_task()."""
    def g(k):
        return task.get(k)
    var = mods(g("CommonParameters.ListOfModsVariable"))
    fixed = mods(g("CommonParameters.ListOfModsFixed"))
    return {
        "task_type": g("TaskType"),
        "search_type": g("SearchParameters.SearchType"),
        "protease": g("CommonParameters.DigestionParams.Protease"),
        "min_peptide_length": g("CommonParameters.DigestionParams.MinPeptideLength"),
        "max_peptide_length": g("CommonParameters.DigestionParams.MaxPeptideLength"),
        "max_mods_per_peptide": g("CommonParameters.DigestionParams.MaxModsForPeptide"),
        "precursor_tolerance": g("CommonParameters.PrecursorMassTolerance"),
        "product_tolerance": g("CommonParameters.ProductMassTolerance"),
        "dissociation": g("CommonParameters.DissociationType"),
        "variable_mods": [f"{c}: {n}" for c, n in var],
        "fixed_mods": [f"{c}: {n}" for c, n in fixed],
        "decoy_type": g("SearchParameters.DecoyType"),
        "q_value_threshold": g("CommonParameters.QValueThreshold"),
        "label_free_quantification": g("SearchParameters.DoLabelFreeQuantification"),
        "match_between_runs": g("SearchParameters.MatchBetweenRuns"),
        "parsimony": g("SearchParameters.DoParsimony"),
        "write_decoys": g("SearchParameters.WriteDecoys"),
        "threads_per_file": g(THREADS_KEY),
    }


def problems(task, mod_files=()):
    """What makes a task TOML one run_search.py will not search with -> [message]. `mod_files`:
    the custom modification files the run installs (each one's path)."""
    out = []
    if task.get("TaskType") != "Search":
        out.append(f"TaskType is {task.get('TaskType')!r}, not \"Search\": run_search.py runs "
                   f"MetaMorpheus' search task only (no calibration, G-PTM-D or crosslink task)")
    defined = set()
    for f in mod_files:
        defined |= custom_mods_defined(f)
    for key in ("CommonParameters.ListOfModsVariable", "CommonParameters.ListOfModsFixed"):
        for cat, name in mods(task.get(key)):
            if cat == "Custom" and name not in defined:
                out.append(f"{key.split('.')[-1]} names the custom modification {name!r}, which "
                           f"no installed modification file defines "
                           f"({', '.join(os.path.basename(f) for f in mod_files) or 'none'}); "
                           f"MetaMorpheus' own Mods/ has no Custom entries to fall back on")
    if THREADS_KEY not in task:
        out.append(f"there is no {THREADS_KEY.split('.')[-1]} line to set the job's CPUs in")
    return out


def run_copy_text(text, threads):
    """The task TOML's text with MaxThreadsToUsePerFile set to `threads` -> (text, old value)."""
    key = THREADS_KEY.split(".")[-1]
    pat = re.compile(rf"^(\s*{key}\s*=\s*)(\S+)(.*)$", re.M)
    m = pat.search(text)
    if not m:
        raise ValueError(f"no {key} line")
    return pat.sub(lambda mm: f"{mm.group(1)}{int(threads)}{mm.group(3)}", text, count=1), \
        m.group(2)
