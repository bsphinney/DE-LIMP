#!/usr/bin/env python3
"""
make_deposit.py  --  Build the repository-deposit package for a finished session, so the
user can put the data in a public ProteomeXchange repository (PRIDE by default, MassIVE as
the alternative) without working out the requirements themselves.

    output/DATA_SUBMISSION/
      HOW_TO_SUBMIT.md / .html   step by step, written for THIS session (the .html opens by
                                 double-click)
      sdrf.tsv                   SDRF-Proteomics v1.1.0, pre-filled from the session; every
                                 value it cannot know is TO-FILL (or a reserved word the spec
                                 allows), never a guess
      protocols.txt              the "Sample processing protocol" and "Data processing
                                 protocol" texts the submission form asks for, from methods.md
      files_to_upload.tsv        every file to upload: PRIDE file type, MassIVE category, path,
                                 size, md5/sha1 where cheap
      prepare_upload.sbatch      NOT run here: packs each Bruker .d run into one archive (PRIDE
                                 requires it) and checksums everything (checksum.txt, PRIDE's
                                 SHA-1 format) as a SLURM job

Nothing here reads raw data: raw sizes come from stat(), and raw checksums are left to
prepare_upload.sbatch -- never computed inline, where finalize may be running on a login node.
Session files up to HASH_MAX_BYTES are hashed inline.

session.py finalize calls ensure_methods() and build() by default. Every requirement quoted
here is cited, with the date it was read, in references/deposit.md.

Usage:
  python3 make_deposit.py --session <session_dir>     # (re)build the package + MANIFEST.txt
"""
import argparse
import csv
import datetime
import glob
import hashlib
import html
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from session import paths_for, read_raw_list      # noqa: E402  the session layout, one place
import make_methods as mm                           # noqa: E402  search_record(): one reader

TO_FILL = "TO-FILL"
HASH_MAX_BYTES = 256 * 1024 ** 2        # hash session files inline only up to this size
SDRF_VERSION = "v1.1.0"                 # SDRF-Proteomics spec version (README.adoc, 2026-01)
TEMPLATE_VERSION = "v1.1.0"             # human / ms-proteomics / dia-acquisition templates
PX_TOOL_VERSION = "2.11.6"              # PRIDE Submission Tool release read 2026-09-24
# PRIDE Submission Tool 2.11.6, SubmissionValidator: "Sample processing protocol must be both
# more than 50 and less than 5000 characters" (same for the data processing protocol and the
# description); "Project title must be less than 500 and more than 30 characters".
PROTOCOL_MIN, PROTOCOL_MAX = 50, 5000
TITLE_MIN, TITLE_MAX = 30, 500
# The same validator: "Filenames must contains only -_.A-Za-z0-9" and "should start with alpha
# numeric characters". (PRIDE's How-to page is stricter and omits the hyphen.)
PRIDE_NAME_OK = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")

# Mass spectrometers -> PSI-MS term, looked up in OLS4 on 2026-09-24. A name not here goes into
# the SDRF as TO-FILL with the recorded name in the notes, never as a nearest match.
INSTRUMENTS = {
    "orbitrap fusion lumos": ("Orbitrap Fusion Lumos", "MS:1002732"),
    "orbitrap exploris 480": ("Orbitrap Exploris 480", "MS:1003028"),
    "orbitrap eclipse": ("Orbitrap Eclipse", "MS:1003029"),
    "orbitrap ascend": ("Orbitrap Ascend", "MS:1003356"),
    "orbitrap astral": ("Orbitrap Astral", "MS:1003378"),
    "orbitrap astral zoom": ("Orbitrap Astral Zoom", "MS:1003442"),
    "q exactive hf": ("Q Exactive HF", "MS:1002523"),
    "q exactive hf-x": ("Q Exactive HF-X", "MS:1002877"),
    "timstof pro": ("timsTOF Pro", "MS:1003005"),
    "timstof pro 2": ("timsTOF Pro 2", "MS:1003230"),
    "timstof scp": ("timsTOF SCP", "MS:1003231"),
    "timstof ht": ("timsTOF HT", "MS:1003404"),
    "timstof ultra": ("timsTOF Ultra", "MS:1003383"),
    "timstof ultra 2": ("timsTOF Ultra 2", "MS:1003412"),
    "timstof flex": ("timsTOF fleX", "MS:1003124"),
}
# comment[proteomics data acquisition method]: the plain OLS labels of PRIDE:0000450 / PRIDE:0000627.
# The spec allows the label or NT=...;AC=..., but the released validator (sdrf-pipelines 0.1.6,
# its bundled dia-acquisition 1.1.0) accepts ONLY the plain "Data-independent acquisition" in a
# DIA SDRF -- measured 2026-09-24, it rejects both NT= forms and diaPASEF. The label passes both.
ACQ_TERMS = {"DIA": "Data-independent acquisition", "DDA": "Data-dependent acquisition"}
# What to type in a TO-FILL column, from the SDRF spec, its SAMPLE-GUIDELINES and the templates'
# validators (bigbio/proteomics-sample-metadata + sdrf-templates, read 2026-09-24).
FILL_HELP = {
    "characteristics[organism]": "NCBI Taxonomy name, lower case (e.g. homo sapiens)",
    "characteristics[organism part]": "UBERON/BTO anatomy term, lower case (e.g. liver). For a "
                                      "cell line: its tissue of origin (e.g. uterine cervix "
                                      "for HeLa) or not applicable",
    "characteristics[disease]": "MONDO/EFO/DOID term; healthy/control samples are normal (not "
                                "'control'); not available if unknown",
    "characteristics[cell type]": "Cell Ontology/BTO/CLO term (e.g. epithelial cell), or not "
                                  "available",
    "characteristics[sex]": "male, female or intersex (a cell line: the donor's, e.g. female "
                            "for HeLa), or not available",
    "characteristics[age]": "e.g. 45Y, 6M, 30Y6M or 40Y-50Y; not available if unknown (the human "
                            "template does not allow not applicable here)",
    "characteristics[biological replicate]": "an integer from 1 within each group, or pooled",
    "comment[proteomics data acquisition method]": "Data-independent acquisition or "
                                                   "Data-dependent acquisition",
    "comment[label]": "label free sample, or one label per row for labelled samples "
                      "(e.g. TMT126)",
    "comment[instrument]": "the PSI-MS term, NT=<name>;AC=MS:<number> (look it up at "
                           "https://www.ebi.ac.uk/ols4/ontologies/ms)",
    "comment[cleavage agent details]": "the PSI-MS enzyme term, e.g. NT=Trypsin;AC=MS:1001251",
    "comment[modification parameters]": "NT=<Unimod name>;AC=UNIMOD:<n>;TA=<residues>;"
                                        "MT=fixed|variable;PP=Anywhere|Protein N-term|...",
    "comment[data file]": "the raw file name exactly as uploaded (a .d folder by its folder "
                          "name)",
    f"factor value[{TO_FILL}]": "rename the header to the variable your groups are (e.g. "
                                "factor value[treatment] or factor value[disease]) and add the "
                                "matching characteristics[...] column with the same values",
}
MASSIVE_CATEGORY = {"RAW": "Raw Spectrum Files", "SEARCH": "Search Engine Files",
                    "FASTA": "Sequence Databases", "SPECTRUM_LIBRARY": "Spectral Libraries",
                    "EXPERIMENTAL_DESIGN": "Supplementary Files",
                    "OTHER": "Supplementary Files"}


# --------------------------------------------------------------------------- manifest --
class Skip(Exception):
    """Raise inside a section to record `[SKIPPED] <name> -- <reason>` with a clean reason."""


class Manifest:
    """The [OK]/[SKIPPED] capture log -- DE-LIMP's safe_section() (R/helpers.R) in Python.

    Silent catch is banned in export paths (CLAUDE.md rule 4): every step of the export runs
    through section(), which records `[OK]      <name> (<s>)` on success and `[SKIPPED] <name>
    -- <reason>` on any failure, then carries on. The log is written as MANIFEST.txt at the
    session root, which is the top level of the session zip."""

    def __init__(self):
        self.lines = []

    def ok(self, name, note="", elapsed=None):
        t = f" ({elapsed:.1f}s)" if elapsed is not None else ""
        self.lines.append(f"[OK]      {name:<50}{t}" + (f" -- {note}" if note else ""))

    def skip(self, name, why):
        why = " ".join(str(why).split())
        if len(why) > 200:
            why = why[:197] + "..."
        self.lines.append(f"[SKIPPED] {name:<50} -- {why}")

    def section(self, name, fn, *args, **kwargs):
        t0 = time.time()
        try:
            note = fn(*args, **kwargs)
        except Skip as e:
            self.skip(name, str(e))
            return False
        except Exception as e:                       # recorded, never swallowed
            self.skip(name, f"{type(e).__name__}: {e}")
            return False
        self.ok(name, note or "", time.time() - t0)
        return True

    @property
    def n_skipped(self):
        return sum(1 for ln in self.lines if ln.startswith("[SKIPPED]"))

    def write(self, path, title):
        with open(path, "w") as fh:
            fh.write(f"{title}\n{'=' * len(title)}\n"
                     f"Written {datetime.datetime.now().isoformat(timespec='seconds')} by "
                     f"make_deposit.py. [OK] = produced; [SKIPPED] = not produced, and why.\n\n"
                     + "\n".join(self.lines) + "\n")
        return path


# ----------------------------------------------------------------------- session facts --
def _load(path):
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError, TypeError):
        return None


def _first(paths):
    for p in paths:
        if p and os.path.isfile(p):
            return p
    return None


def find_params(p):
    """The search parameters file the session holds (DIA-NN cfg / Sage json / FragPipe
    .workflow / Radiant config), newest convention first."""
    wf, inp = p["workflow_dir"], p["input_dir"]
    cands = [os.path.join(wf, n) for n in ("params.cfg", "params.json")]
    for pat in ("*.cfg", "params*.json", "sage*.json", "*.workflow", "*.radiantConfig"):
        cands += sorted(glob.glob(os.path.join(wf, pat)))
    for pat in ("params.*", "*.cfg", "sage_config*.json"):
        cands += sorted(glob.glob(os.path.join(inp, pat)))
    return _first(c for c in cands if not c.endswith((".rationale.json", "manifest.json")))


def read_conditions(path):
    """conditions.csv (collect_conditions.py schema: File.Name,Group[,Batch,...]) -> rows."""
    if not path or not os.path.isfile(path):
        return None
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows or "File.Name" not in rows[0]:
        return None
    return rows


def gather(session_dir):
    """Everything the package and the Methods need, read from the session's own files."""
    p = paths_for(session_dir)
    f = {"p": p, "session": os.path.basename(p["session_dir"])}
    f["raws"] = read_raw_list(session_dir)
    f["wf"] = _load(p["workflow_manifest"]) or {}
    f["wf_path"] = p["workflow_manifest"] if os.path.isfile(p["workflow_manifest"]) else None
    f["run_manifest"] = _load(os.path.join(p["repro_dir"], "run_manifest.json")) or {}
    f["search_prov_path"] = _first([p["search_prov"]] + sorted(
        glob.glob(os.path.join(p["search_out"], "*", "search_provenance.json"))))
    f["params"] = find_params(p)
    f["fasta_meta_path"] = _first([p["fasta_meta"]] + sorted(
        glob.glob(os.path.join(p["input_dir"], "*.meta.json"))))
    f["fasta_meta"] = _load(f["fasta_meta_path"]) if f["fasta_meta_path"] else None
    f["de_prov_path"] = _first([os.path.join(p["de_dir"], "de_provenance.json")])
    f["de_prov"] = _load(f["de_prov_path"]) if f["de_prov_path"] else None
    f["conditions"] = read_conditions(p["conditions"])
    f["methods_params"] = _load(os.path.join(p["output_dir"], "methods_params.json")) or {}
    f["srec"] = mm.search_record(f["params"], f["search_prov_path"], f["wf_path"])
    q = f["run_manifest"].get("query") or {}
    f["acquisition"] = ((f["wf"].get("acquisition") or q.get("acquisition") or "").upper()
                        or None)
    rep = f["methods_params"].get("representative") or {}
    f["mode"] = rep.get("mode") if not f["methods_params"].get("from_session_record") else None
    f["instrument"] = (f["methods_params"].get("instrument")
                       or next(iter(f["wf"].get("instruments") or []), None)
                       or q.get("instrument") or None)
    f["instrument_source"] = ("methods_params.json (read from the raw files)"
                              if f["methods_params"].get("instrument")
                              and not f["methods_params"].get("from_session_record")
                              else "the session record (workflow manifest)")
    fm = f["fasta_meta"] or {}
    f["taxid"] = fm.get("taxid") or f["wf"].get("organism_taxid") or q.get("organism_taxid")
    f["organism"] = fm.get("organism") or None
    return f


def _size(path):
    """Bytes of a file, or of everything under a directory (stat only; nothing is read)."""
    if os.path.isdir(path):
        total = 0
        for dp, _, fns in os.walk(path):
            for fn in fns:
                try:
                    total += os.lstat(os.path.join(dp, fn)).st_size
                except OSError:
                    pass
        return total
    return os.path.getsize(path)


def _hash2(path):
    s, m = hashlib.sha1(), hashlib.md5()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 22), b""):
            s.update(b)
            m.update(b)
    return s.hexdigest(), m.hexdigest()


def _human(n):
    if n is None:
        return "?"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0


def location_of(path):
    """Where a recorded absolute path lives, in words the user can act on."""
    if path.startswith("/quobyte/"):
        return "HIVE (/quobyte)"
    if path.startswith("/nfs/lssc0/flinders/"):
        return "HIVE (Flinders share)"
    # finalize running ON HIVE (its Quobyte storage is mounted here): a path it can see is a
    # HIVE path, home directories included
    if os.path.isdir("/quobyte") and os.path.exists(path):
        return "HIVE"
    if path.startswith("/Volumes/") or re.match(r"^[A-Za-z]:", path):
        return ("the machine that ran the analysis (if this is a share HIVE mounts, "
                "`bash scripts/hive_path.sh <path>` gives the HIVE path)")
    return "the machine that ran the analysis"


def safe_name(name):
    """A name PRIDE's Submission Tool accepts ([A-Za-z0-9][-_.A-Za-z0-9]*)."""
    s = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
    return s if re.match(r"^[A-Za-z0-9]", s) else "f" + s


# ----------------------------------------------------------------------- upload plan --
# Search outputs by engine: (regex on the file name, PRIDE file type, requirement, note).
# PRIDE Submission Tool 2.11.6 file types: RAW RESULT PEAK SEARCH QUANT GEL FASTA
# SPECTRUM_LIBRARY MS_IMAGE_DATA OPTICAL_IMAGE OTHER EXPERIMENTAL_DESIGN. PRIDE's guidelines
# v2.2.0 call SEARCH "ANALYSIS" and RESULT "STANDARD"; the tool still shows the old names.
SEARCH_FILES = [
    # DIA-NN (PRIDE guidelines v2.2.0, "DIA-NN")
    (r"^(no_norm_)?report\.(parquet|tsv)$", "SEARCH", "required",
     "DIA-NN main report (PRIDE: mandatory)"),
    (r"\.site_report\.parquet$", "OTHER", "required", "DIA-NN site report (mandatory if produced)"),
    (r"\.log\.txt$", "OTHER", "required", "DIA-NN log (PRIDE: mandatory for DIA-NN)"),
    (r"_matrix\.tsv$", "SEARCH", "recommended", "DIA-NN quantity matrix"),
    (r"\.stats\.tsv$", "OTHER", "optional", "DIA-NN run statistics"),
    (r"\.protein_description\.tsv$", "OTHER", "optional", "DIA-NN protein descriptions"),
    (r"\.manifest\.txt$", "OTHER", "optional", "DIA-NN report manifest"),
    (r"^step3_assembly\.parquet$", "SPECTRUM_LIBRARY", "recommended",
     "empirical library assembled from these runs (PRIDE: recommended)"),
    (r"(lib|library)[^/]*\.parquet$", "SPECTRUM_LIBRARY", "recommended",
     "spectral library (PRIDE: mandatory if a library search was performed)"),
    (r"\.speclib$", "SPECTRUM_LIBRARY", "recommended",
     "predicted spectral library; can also be regenerated from the FASTA and pinned DIA-NN "
     "version"),
    (r"^params\.(resolved|base)\.cfg$", "OTHER", "recommended",
     "the DIA-NN parameters the search ran with"),
    (r"^search_provenance\.json$", "OTHER", "recommended",
     "engine, version and exact command of this search"),
    # Sage
    (r"^results\.sage\.(tsv|parquet)$", "SEARCH", "required", "Sage PSM results"),
    (r"^lfq\.(tsv|parquet)$", "SEARCH", "recommended", "Sage label-free quantification"),
    (r"^results\.json$", "OTHER", "recommended", "Sage run record (parameters + version)"),
    # FragPipe (PRIDE guidelines v2.2.0, "FragPipe/MSFragger")
    (r"^psm\.tsv$", "SEARCH", "required", "FragPipe PSMs (PRIDE: mandatory)"),
    (r"^protein\.tsv$", "SEARCH", "required", "FragPipe proteins (PRIDE: mandatory)"),
    (r"^combined_.*\.tsv$", "SEARCH", "recommended", "FragPipe combined results"),
    (r"^fragpipe\.workflow$", "OTHER", "required", "FragPipe parameters (PRIDE: mandatory)"),
    (r"\.fp-manifest$", "OTHER", "required", "FragPipe manifest (PRIDE: mandatory)"),
    (r"^log_.*\.txt$", "OTHER", "optional", "FragPipe log"),
]
# never listed: intermediates, and anything the search recreates
SKIP_DIRS = {"quant", "quant_step2", "quant_step4", "quant_step2_orig", "xic", "libpriv",
             "window_probe", "logs", "tmp", "temp"}


def _classify(name):
    for pat, ftype, req, note in SEARCH_FILES:
        if re.search(pat, name, re.I):
            return ftype, req, note
    return None


def plan_uploads(f):
    """Every file to deposit, with its type, where it is, and what must happen first."""
    p, rows = f["p"], []

    def add(kind, src, upload, ftype, req, note, rel_zip, size=None, action=""):
        rows.append({"kind": kind, "source_path": src, "upload_name": upload,
                     "pride_file_type": ftype,
                     "massive_category": MASSIVE_CATEGORY.get(ftype, "Supplementary Files"),
                     "requirement": req, "in_session_zip": rel_zip, "size_bytes": size,
                     "md5": "", "sha1": "", "before_upload": action, "notes": note,
                     "location": location_of(src) if kind in ("dir", "file") else
                     "this session folder"})

    for path in f["raws"]:
        base = os.path.basename(path.rstrip("/"))
        reach = os.path.exists(path)
        is_dir = os.path.isdir(path) if reach else base.lower().endswith(".d")
        sname = safe_name(base)
        size = _size(path) if reach else None
        act, note = [], "raw instrument data (PRIDE: mandatory)"
        if sname != base:
            act.append(f"rename: '{base}' has characters PRIDE rejects; renamed to {sname} in "
                       f"the SDRF -- rename the {'folder' if is_dir else 'file'} to match")
        if is_dir:
            act.append("compress to ONE archive per run (prepare_upload.sbatch)")
            note += "; directory format: PRIDE requires each .d folder compressed on its own"
        if not reach:
            note += "; not reachable from where finalize ran -- size not measured"
        add("dir" if is_dir else "file", path, sname + (".tar.gz" if is_dir else ""), "RAW",
            "required", note, "no (raw data is not copied into the session)", size,
            "; ".join(act))

    search = p["search_out"]
    seen = {}
    if os.path.isdir(search):
        for root, dirs, fns in os.walk(search):
            depth = os.path.relpath(root, search).count(os.sep) + (root != search)
            dirs[:] = [d for d in sorted(dirs) if d.lower() not in SKIP_DIRS
                       and not d.lower().endswith(".d") and depth < 2]
            for fn in sorted(fns):
                c = _classify(fn)
                if not c:
                    continue
                full = os.path.join(root, fn)
                up = fn if root == search else safe_name(
                    os.path.relpath(root, search).replace(os.sep, "_") + "_" + fn)
                add("copy", full, safe_name(up), *c,
                    os.path.relpath(full, p["session_dir"]), _size(full))
    if f["params"]:
        pb = os.path.basename(f["params"])
        add("copy", f["params"], safe_name(pb),
            *(_classify(pb) or ("OTHER", "recommended", "search parameters file given to the "
                                                        "engine")),
            os.path.relpath(f["params"], p["session_dir"]), _size(f["params"]))
    fasta = _first([p["fasta"], (f["fasta_meta"] or {}).get("fasta")])
    if fasta:
        fm = f["fasta_meta"] or {}
        custom = (fm.get("n_contaminants_appended") or fm.get("staged_file")
                  or fm.get("content_used") in (None, "unknown", "as_staged"))
        add("copy", fasta, safe_name(os.path.basename(fasta)), "FASTA",
            "required" if custom else "optional",
            ("searched database. Upload it: it has contaminants appended or is a staged/"
             "supplied copy, so nobody can download the exact file elsewhere (PRIDE: \"if a "
             "third party cannot obtain the exact same database independently, you must "
             "provide it\")") if custom else
            "searched database (a plain public UniProt proteome: PRIDE accepts name + release "
            "in the data processing protocol instead)",
            os.path.relpath(fasta, p["session_dir"]) if fasta.startswith(p["session_dir"])
            else "no", _size(fasta) if os.path.isfile(fasta) else None)
    add("copy", os.path.join(p["deposit_dir"], "sdrf.tsv"), "sdrf.tsv", "EXPERIMENTAL_DESIGN",
        "recommended", "SDRF-Proteomics sample metadata (PRIDE: strongly recommended) -- fill "
        "every TO-FILL first", "output/DATA_SUBMISSION/sdrf.tsv", None)
    for pat in ("DE_*.csv", "reproducibility_log.R", "de_provenance.json"):
        for full in sorted(glob.glob(os.path.join(p["de_dir"], pat))):
            add("copy", full, safe_name(os.path.basename(full)), "OTHER", "optional",
                "differential-expression result / record",
                os.path.relpath(full, p["session_dir"]), _size(full))

    # unique upload names: PRIDE stores a dataset's files in one folder
    for r in rows:
        seen[r["upload_name"]] = seen.get(r["upload_name"], 0) + 1
    for r in rows:
        if seen[r["upload_name"]] > 1:
            r["before_upload"] = "; ".join(x for x in (
                r["before_upload"], "duplicate upload name -- rename one copy") if x)
        if not PRIDE_NAME_OK.match(r["upload_name"]):
            r["before_upload"] = "; ".join(x for x in (
                r["before_upload"], "name has characters PRIDE rejects -- rename") if x)
    return rows


def hash_small(rows):
    """md5 + sha1 of session files small enough to hash here. Raw data: never (see module doc).
    sdrf.tsv: not yet -- it is edited before upload; prepare_upload.sbatch hashes it."""
    n = 0
    for r in rows:
        src = r["source_path"]
        if r["kind"] != "copy" or r["upload_name"] == "sdrf.tsv":
            r["md5"] = r["sha1"] = "computed by prepare_upload.sbatch"
            continue
        if not os.path.isfile(src):
            continue
        if (r["size_bytes"] or 0) > HASH_MAX_BYTES:
            r["md5"] = r["sha1"] = (f"computed by prepare_upload.sbatch (> "
                                    f"{HASH_MAX_BYTES // 1024 ** 2} MB)")
            continue
        r["sha1"], r["md5"] = _hash2(src)
        n += 1
    return n


UPLOAD_COLS = ["upload_name", "pride_file_type", "massive_category", "requirement",
               "source_path", "location", "in_session_zip", "size_bytes", "md5", "sha1",
               "before_upload", "notes"]


def write_upload_list(out, rows):
    path = os.path.join(out, "files_to_upload.tsv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(UPLOAD_COLS)
        for r in rows:
            w.writerow(["" if r.get(c) is None else r.get(c) for c in UPLOAD_COLS])
    return path


# ----------------------------------------------------------------------------- SDRF --
def _instrument_term(name):
    if not name:
        return None
    key = re.sub(r"^(thermo|bruker)\s+", "", name.strip().lower())
    hit = INSTRUMENTS.get(key)
    return f"NT={hit[0]};AC={hit[1]}" if hit else None


def _mod_term(m):
    """SDRF comment[modification parameters]: NT first, then AC (README 'Table Cell values')."""
    if m["unimod"] not in mm.UNIMOD:
        return None
    bits = [f"NT={m['name']}", f"AC=UNIMOD:{m['unimod']}"]
    if m["targets"]:
        bits.append("TA=" + ",".join(m["targets"]))
    bits += [f"MT={m['type']}", f"PP={m['position']}"]
    return ";".join(bits)


def _run_name(path):
    base = os.path.basename(path.rstrip("/"))
    for ext in (".d", ".raw", ".mzml", ".wiff", ".mzxml", ".mgf"):
        if base.lower().endswith(ext):
            return base[: -len(ext)]
    return base


def build_sdrf(f):
    """(header, rows, sources, to_fill). `sources`: one (column, value, where-from) per column;
    `to_fill`: what the user must supply, column by column."""
    raws, conds = f["raws"], f["conditions"] or []
    notes, to_fill = [], []
    cond_by = {}
    for r in conds:
        fn = (r.get("File.Name") or "").strip()
        for k in {fn, fn.lower(), _run_name(fn), _run_name(fn).lower()}:
            cond_by.setdefault(k, r)
    if raws:
        items = [(p, os.path.basename(p.rstrip("/"))) for p in raws]
    elif conds:
        items = [(None, None) for _ in conds]
        notes.append(("comment[data file]", TO_FILL, "input/raw_files.txt is missing, so rows "
                      "come from conditions.csv and the raw file names are unknown"))
    else:
        raise Skip("no raw file list (input/raw_files.txt) and no conditions.csv: "
                   "nothing to describe")

    human = str(f["taxid"] or "") == "9606"
    acq = f["acquisition"]
    srec = f["srec"] or {}

    # organism: the organism the user confirmed at step 3 (it chose the FASTA)
    org = f["organism"]
    org_val = (re.sub(r"\s*\(.*\)\s*$", "", org).strip().lower() if org else None)
    org_src = (f"search database organism '{org}' (taxid {f['taxid']}), the organism you "
               f"confirmed" if org else None)
    if org and "(" in org:
        org_src += " -- strain text removed; check it is the NCBI Taxonomy name"
    if not org and human:
        # the one taxid named without a lookup (NCBITaxon:9606 = Homo sapiens, OLS4)
        org_val, org_src = "homo sapiens", "organism_taxid 9606 recorded for this session"
    elif not org and f["taxid"]:
        org_src = (f"taxid {f['taxid']} is recorded but not its name -- enter the NCBI "
                   f"Taxonomy name for it")
    instr = _instrument_term(f["instrument"])
    mode = f["mode"] or ""
    if acq == "DIA" and "dia-pasef" in mode.lower():
        # The spec (README, 2026-08) and sdrf-templates main allow NT=diaPASEF;AC=PRIDE:0000650
        # here; the released validator does not (see ACQ_TERMS). Say diaPASEF in the notes.
        acq_val = ACQ_TERMS["DIA"]
        acq_src = ("dia-PASEF read from the raw files; written as the parent DIA term, the only "
                   "value the released SDRF validator (sdrf-pipelines 0.1.6) accepts in a "
                   "dia-acquisition SDRF (the spec also allows NT=diaPASEF;AC=PRIDE:0000650)")
    elif acq in ("DIA", "DDA"):
        acq_val, acq_src = ACQ_TERMS[acq], "acquisition detected from the data (step 2)"
    else:
        acq_val, acq_src = None, None
    cl = srec.get("cleavage") or {}
    enz = f"NT={cl['name']};AC={cl['ac']}" if cl.get("name") else None
    mods = [m for m in srec.get("mods") or [] if not m.get("label")]
    mod_vals = [(_mod_term(m), m) for m in mods]
    lab = srec.get("labelled")
    label_val = ("label free sample" if lab and lab["value"] is False else None)

    templates = (["human"] if human else []) + (["dia-acquisition"] if acq == "DIA"
                                                else ["ms-proteomics"])
    cols = ["source name", "characteristics[organism]", "characteristics[organism part]",
            "characteristics[disease]", "characteristics[cell type]"]
    if human:
        cols += ["characteristics[sex]", "characteristics[age]"]
    cols += ["characteristics[biological replicate]", "assay name", "technology type",
             "comment[proteomics data acquisition method]", "comment[label]",
             "comment[instrument]", "comment[cleavage agent details]"]
    cols += ["comment[modification parameters]"] * len(mod_vals)
    tol_cols = [(c, srec.get(k)) for c, k in (("comment[precursor mass tolerance]", "ms1_tol"),
                                              ("comment[fragment mass tolerance]", "ms2_tol"))
                if srec.get(k) and srec[k].get("symmetric", True)]
    cols += [c for c, _ in tol_cols]
    mp_rep = f["methods_params"].get("representative") or {}
    ms1_range = None
    if (not f["methods_params"].get("from_session_record") and mp_rep.get("mz_low") is not None
            and mp_rep.get("mz_high") is not None):
        ms1_range = f"{mm._g(mp_rep['mz_low'])}m/z-{mm._g(mp_rep['mz_high'])}m/z"
        cols.append("comment[ms1 scan range]")
    cols += ["comment[fraction identifier]", "comment[technical replicate]",
             "comment[data file]", "comment[sdrf version]"]
    cols += ["comment[sdrf template]"] * len(templates)
    cols += ["comment[sdrf annotation tool]", f"factor value[{TO_FILL}]"]

    skill_ver = (_load(os.path.join(HERE, "..", ".claude-plugin", "plugin.json")) or {}).get(
        "version") or "0.0.0"
    rows, rep_n = [], {}
    for i, (path, base) in enumerate(items):
        if path:
            run = _run_name(path)
            c = cond_by.get(base) or cond_by.get(run) or cond_by.get(run.lower())
        else:
            c = conds[i]
            run = _run_name(c.get("File.Name", f"run_{i + 1}"))
        group = (c or {}).get("Group", "").strip() or None
        if group:
            rep_n[group] = rep_n.get(group, 0) + 1
        row = [run, org_val or TO_FILL, TO_FILL, TO_FILL, TO_FILL]
        if human:
            row += [TO_FILL, TO_FILL]
        row += [str(rep_n[group]) if group else TO_FILL, run,
                "proteomic profiling by mass spectrometry",
                acq_val or TO_FILL, label_val or TO_FILL, instr or TO_FILL, enz or TO_FILL]
        row += [v or TO_FILL for v, _ in mod_vals]
        row += [f"{mm._g(t['value'])} {t['unit']}" for _, t in tol_cols]
        if ms1_range:
            row.append(ms1_range)
        row += ["1", "1", safe_name(base) if base else TO_FILL, SDRF_VERSION]
        row += [f"{t} {TEMPLATE_VERSION}" for t in templates]
        row += [f"ucdavis-proteomics-core-pipeline v{skill_ver}", group or TO_FILL]
        rows.append(row)
        if path and not c and conds:
            notes.append(("factor value", TO_FILL, f"{base} has no row in conditions.csv"))

    # where every column came from, and what the user must fill
    src = [("source name", "run name", "the raw file name without extension -- replace with "
            "your own sample IDs if you have them"),
           ("characteristics[organism]", org_val or TO_FILL, org_src or
            "no search-database record (input/search.fasta.meta.json)"),
           ("characteristics[organism part]", TO_FILL, "not in any session record"),
           ("characteristics[disease]", TO_FILL, "not in any session record"),
           ("characteristics[cell type]", TO_FILL, "not in any session record")]
    if human:
        src += [("characteristics[sex]", TO_FILL, "not in any session record"),
                ("characteristics[age]", TO_FILL, "not in any session record")]
    src += [("characteristics[biological replicate]", "1..n per group" if conds else TO_FILL,
             "numbered within each conditions.csv Group, one per raw file -- the DE design "
             "treated each run as a separate biological sample; correct it if some runs are "
             "re-injections" if conds else "no conditions.csv"),
            ("assay name", "run name", "the raw file name without extension"),
            ("technology type", "proteomic profiling by mass spectrometry",
             "fixed value required by the spec"),
            ("comment[proteomics data acquisition method]", acq_val or TO_FILL,
             acq_src or "acquisition not recorded in the session"),
            ("comment[label]", label_val or TO_FILL,
             (lab or {}).get("source", "search parameters not readable") if label_val else
             ("the search used labels -- enter each channel's label term and one row per "
              "channel" if lab and lab["value"] else
              "search parameters could not be read, so labelling is unknown")),
            ("comment[instrument]", instr or TO_FILL,
             f"{f['instrument']} ({f['instrument_source']})" if instr else
             (f"recorded as '{f['instrument']}', which has no PSI-MS term in this skill's "
              f"verified table -- enter its NT=...;AC=MS:... term" if f["instrument"] else
              "instrument not recorded")),
            ("comment[cleavage agent details]", enz or TO_FILL,
             f"the search's in-silico cleavage rule ({cl['rule']}, {cl['source']}) -- confirm "
             f"it matches the wet-lab enzyme (e.g. a Trypsin/Lys-C mix)" if enz else
             (f"the search's cleavage rule {cl['rule']} has no mapped enzyme term" if cl else
              "search parameters could not be read"))]
    for v, m in mod_vals:
        src.append(("comment[modification parameters]", v or TO_FILL,
                    m["source"] if v else f"'{m['name']}' ({m['source']}) has no verified "
                    f"Unimod term here -- enter NT=...;AC=UNIMOD:...;TA=...;MT=...;PP=..."))
    if not mod_vals:
        notes.append(("comment[modification parameters]", "(column omitted)",
                      "no modifications could be read from the search parameters"))
    for c, t in tol_cols:
        src.append((c, f"{mm._g(t['value'])} {t['unit']}", t["source"]))
    if len(tol_cols) < 2:
        notes.append(("comment[precursor/fragment mass tolerance]", "(column omitted)",
                      srec.get("tol_note") or "no single fixed tolerance in the search record"))
    if ms1_range:
        src.append(("comment[ms1 scan range]", ms1_range, "read from the raw files "
                    "(GlobalMetadata MzAcqRange)"))
    src += [("comment[fraction identifier]", "1", "assumed unfractionated: one raw file per "
             "sample in the DE design -- change it if you fractionated"),
            ("comment[technical replicate]", "1", "no technical replicates recorded -- change "
             "it if some runs are repeat injections"),
            ("comment[data file]", "raw file name", "input/raw_files.txt (a .d folder is "
             "named directly, per the SDRF spec)"),
            ("comment[sdrf version]", SDRF_VERSION, "SDRF-Proteomics spec version"),
            ("comment[sdrf template]", ", ".join(f"{t} {TEMPLATE_VERSION}" for t in templates),
             ("human" if human else "non-human: add the matching sample template (vertebrates, "
              "invertebrates or plants) and its required columns") + "; "
             + ("dia-acquisition implies ms-proteomics" if acq == "DIA" else "ms-proteomics")),
            ("comment[sdrf annotation tool]", f"ucdavis-proteomics-core-pipeline v{skill_ver}",
             "this skill"),
            (f"factor value[{TO_FILL}]", "conditions.csv Group" if conds else TO_FILL,
             "conditions.csv gives each run's group label (filled in) but not what variable "
             "the groups are" if conds else "no conditions.csv")]
    seen = set()
    for col, val, why in src + notes:
        if (val == TO_FILL or TO_FILL in str(col)) and col not in seen:
            seen.add(col)
            to_fill.append((col, FILL_HELP.get(col, why), why))
    return cols, rows, src + notes, to_fill


def write_sdrf(f, out, info):
    cols, rows, sources, to_fill = build_sdrf(f)
    path = os.path.join(out, "sdrf.tsv")
    with open(path, "w", newline="") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    n_cells = sum(r.count(TO_FILL) for r in rows)
    info.update(path=path, n_rows=len(rows), n_to_fill=n_cells, sources=sources,
                to_fill=to_fill, columns=cols)
    return f"{len(rows)} row(s), {n_cells} TO-FILL cell(s) in {len(to_fill)} column(s)"


# ------------------------------------------------------------------------ protocols --
def md_sections(text):
    """{heading: body} for the '## ' sections of a Markdown file, in order."""
    out, cur = {}, None
    for ln in text.splitlines():
        if ln.startswith("## "):
            cur = ln[3:].strip()
            out[cur] = []
        elif cur is not None:
            out[cur].append(ln)
    return {k: "\n".join(v).strip() for k, v in out.items()}


def _plain(md):
    """Markdown body -> the plain paragraphs a web form takes (tables and notes dropped)."""
    keep = []
    for ln in md.splitlines():
        s = ln.strip()
        if s.startswith("|") or (s.startswith("*") and s.endswith("*") and not s.startswith("**")):
            continue
        s = re.sub(r"^>\s*", "", s)
        s = s.replace("**", "").replace("`", "")
        keep.append(s)
    txt = "\n".join(keep)
    return re.sub(r"\n{3,}", "\n\n", txt).strip()


SAMPLE_PREP_TO_FILL = (f"[{TO_FILL}: sample preparation -- how proteins were extracted, "
                       "reduced and alkylated (reagents), digested (enzyme, enzyme:protein ratio, "
                       "time), cleaned up, and how much peptide was injected.]")


def build_protocols(methods_md):
    if not methods_md or not os.path.isfile(methods_md):
        raise Skip("no methods file (output/methods.md) to draw the protocols from -- see "
                   "the methods line above")
    with open(methods_md, encoding="utf-8") as fh:
        sec = md_sections(fh.read())
    sample = [SAMPLE_PREP_TO_FILL] + [_plain(sec[h]) for h in ("Liquid chromatography",
                                                              "Mass spectrometry") if h in sec]
    data = [_plain(sec[h]) for h in ("Sequence database", "Database search",
                                     "Differential expression") if h in sec]
    if not data:
        data = [f"[{TO_FILL}: how the data were searched and analysed -- methods.md has no "
                "Sequence database / Database search sections]"]
    return "\n\n".join(x for x in sample if x), "\n\n".join(x for x in data if x)


def write_protocols(out, methods_md, info):
    sample, data = build_protocols(methods_md)
    tag = re.compile(r"\[[^\]]*(confirm|TO-FILL)[^\]]*\]")
    blocks, warn = [], []
    for title, body in (("SAMPLE PROCESSING PROTOCOL", sample),
                        ("DATA PROCESSING PROTOCOL", data)):
        n = len(body)
        if not PROTOCOL_MIN < n < PROTOCOL_MAX:
            warn.append(f"{title.lower()} is {n} characters; PRIDE accepts more than "
                        f"{PROTOCOL_MIN} and fewer than {PROTOCOL_MAX}")
        blocks.append(f"===== {title} ({n} characters as generated) =====\n\n{body}\n")
    n_tags = len(tag.findall(sample + data))
    head = (f"Protocol texts for the repository submission form\n"
            f"Generated {datetime.date.today().isoformat()} by make_deposit.py from "
            f"{os.path.relpath(methods_md, os.path.dirname(os.path.dirname(out)))}.\n\n"
            f"Paste each block into the matching field of the PRIDE Submission Tool. PRIDE's "
            f"tool (v{PX_TOOL_VERSION}) requires each protocol to be more than {PROTOCOL_MIN} "
            f"and fewer than {PROTOCOL_MAX} characters.\n"
            f"BEFORE PASTING: replace every {TO_FILL} and resolve every '[... — confirm]' tag "
            f"({n_tags} in this file). Tagged values are defaults or blanks, not facts about "
            f"your samples.\n"
            + "".join(f"WARNING: {w}\n" for w in warn) + "\n")
    path = os.path.join(out, "protocols.txt")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(head + "\n".join(blocks))
    info.update(path=path, n_tags=n_tags, lengths=(len(sample), len(data)), warnings=warn)
    return f"{n_tags} tag(s) to resolve" + (f"; {'; '.join(warn)}" if warn else "")


# ------------------------------------------------------------------- prepare script --
def queue_lines():
    """The #SBATCH queue lines from the skill's one queue rule (run_search.slurm_queue)."""
    try:
        import run_search
        import diann_parallel
        part, acct, qos = run_search.slurm_queue(preemptible_ok=True)
        req = diann_parallel.needs_requeue(part, qos)
        src = ("detected from your SLURM associations (run_search.slurm_queue)"
               if shutil.which("sacctmgr") else
               "run_search.slurm_queue() default for a machine without SLURM (publicgrp/low, "
               "open to every HIVE allocation) -- check it")
    except (Exception, SystemExit) as e:
        part = acct = qos = None
        req, src = False, f"not detected ({e}); pass --partition/--account/--qos to sbatch"
    lines = [f"#SBATCH --{k}={v}" for k, v in (("partition", part), ("account", acct),
                                                ("qos", qos)) if v]
    if req:
        lines.append("#SBATCH --requeue")
    return lines, src


def write_prep_script(f, out, rows):
    items = [(r["kind"], r["source_path"], r["upload_name"]) for r in rows]
    bad = [i for i in items if "\t" in i[1] or "\n" in i[1]]
    if bad:
        raise Skip(f"a path contains a tab/newline: {bad[0][1]!r}")
    qlines, qsrc = queue_lines()
    skill_ver = (_load(os.path.join(HERE, "..", ".claude-plugin", "plugin.json")) or {}).get(
        "version") or "?"
    listing = "\n".join("\t".join(i) for i in items)
    script = f"""#!/bin/bash -l
#SBATCH --job-name=deposit_prep
#SBATCH --output={out}/prepare_upload_%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=24:00:00
{chr(10).join(qlines)}
# prepare_upload.sbatch -- written by make_deposit.py (ucdavis-proteomics-core-pipeline
# v{skill_ver}) on {datetime.date.today().isoformat()} for session {f['session']}.
# The skill did NOT run it. Submit it yourself when the SDRF is filled in:
#
#   sbatch {out}/prepare_upload.sbatch
#   STAGE=/another/folder sbatch {out}/prepare_upload.sbatch   # stage elsewhere
#   COPY_RAW=1 sbatch {out}/prepare_upload.sbatch    # copy single-file raws instead of
#                                                    # linking them (Globus skips links)
#
# Per raw run: a Bruker/Agilent .d folder becomes ONE archive <run>.d.tar.gz (PRIDE: .d folders
# "must be compressed", one run per archive); a single-file raw (.raw) is linked into STAGE
# unchanged. Session files (search results, FASTA, sdrf.tsv, ...) are copied in. Then SHA-1 and
# MD5 of every file: STAGE/checksum.txt is PRIDE's format (name<TAB>sha1, needed for a Globus
# submission) and STAGE/md5sums.txt is the same list as md5.
# Safe to re-run (e.g. after editing sdrf.tsv, or if SLURM stopped it): finished archives are
# kept (STAGE/.done); session files are re-copied and re-hashed every time.
# Queue: {qsrc}.
# Override on the command line: sbatch --partition=... --account=... --qos=...
set -euo pipefail
if [ -z "${{SLURM_JOB_ID:-}}" ] && [ "${{RUN_HERE:-0}}" != "1" ]; then
  echo "prepare_upload.sbatch reads every raw file: submit it with sbatch, never run it on a" >&2
  echo "login node. (On a laptop whose raw files are local: RUN_HERE=1 bash <this file>)" >&2
  exit 2
fi
PKG={shlex.quote(out)}
STAGE="${{STAGE:-$PKG/upload_staging}}"
COPY_RAW="${{COPY_RAW:-0}}"
THREADS="${{SLURM_CPUS_PER_TASK:-4}}"
mkdir -p "$STAGE/.done"
if command -v module >/dev/null 2>&1; then module load pigz >/dev/null 2>&1 || true; fi
if command -v pigz >/dev/null 2>&1; then GZ="pigz -p $THREADS"; else GZ="gzip"; fi

hash2() {{  # "<sha1> <md5>" of a file, read once when python3 is available
  if command -v python3 >/dev/null 2>&1; then
    python3 - "$1" <<'PY'
import hashlib, sys
s, m = hashlib.sha1(), hashlib.md5()
with open(sys.argv[1], "rb") as fh:
    for b in iter(lambda: fh.read(1 << 22), b""):
        s.update(b)
        m.update(b)
print(s.hexdigest(), m.hexdigest())
PY
  elif command -v sha1sum >/dev/null 2>&1; then
    echo "$(sha1sum "$1" | cut -d' ' -f1) $(md5sum "$1" | cut -d' ' -f1)"
  else
    echo "$(shasum -a 1 "$1" | cut -d' ' -f1) $(md5 -q "$1")"
  fi
}}

: > "$STAGE/checksum.txt.part"
: > "$STAGE/md5sums.txt.part"
n=0; missing=0
while IFS=$'\\t' read -r kind src name; do
  [ -n "$kind" ] || continue
  n=$((n + 1))
  dest="$STAGE/$name"; done_f="$STAGE/.done/$name"
  if [ ! -e "$src" ]; then
    echo "MISSING: $src" >&2; missing=$((missing + 1)); continue
  fi
  case "$kind" in
    dir)
      if [ ! -s "$done_f" ] || [ ! -s "$dest" ]; then
        echo "[$n] archiving $src -> $name"
        tar -C "$(dirname "$src")" -cf - "$(basename "$src")" | $GZ > "$dest.part"
        $GZ -t "$dest.part"        # PRIDE checks every archive extracts; check first
        mv -f "$dest.part" "$dest"
        hash2 "$dest" > "$done_f"
      fi ;;
    file)
      if [ "$COPY_RAW" = "1" ]; then
        if [ -L "$dest" ] || [ ! -f "$dest" ]; then rm -f "$dest"; cp -p "$src" "$dest"; fi
      else
        ln -sfn "$src" "$dest"
      fi
      [ -s "$done_f" ] || hash2 "$src" > "$done_f" ;;
    copy)
      cp -p "$src" "$dest"
      hash2 "$dest" > "$done_f" ;;
  esac
  read -r sha md5 < "$done_f"
  printf '%s\\t%s\\n' "$name" "$sha" >> "$STAGE/checksum.txt.part"
  printf '%s  %s\\n' "$md5" "$name" >> "$STAGE/md5sums.txt.part"
done <<'LIST'
{listing}
LIST
mv -f "$STAGE/checksum.txt.part" "$STAGE/checksum.txt"
mv -f "$STAGE/md5sums.txt.part" "$STAGE/md5sums.txt"
if grep -q '{TO_FILL}' "$PKG/sdrf.tsv" 2>/dev/null; then
  echo "WARNING: $PKG/sdrf.tsv still has {TO_FILL} cells. Fill them, then re-run this script" >&2
  echo "so the staged copy and its checksum match." >&2
fi
echo "Staged $((n - missing)) of $n file(s) in $STAGE ($missing missing)."
echo "Upload everything in $STAGE except md5sums.txt and .done/; checksum.txt goes with a"
echo "Globus submission (the PRIDE Submission Tool writes its own)."
[ "$missing" -eq 0 ]
"""
    path = os.path.join(out, "prepare_upload.sbatch")
    with open(path, "w") as fh:
        fh.write(script)
    os.chmod(path, 0o755)
    n_dir = sum(1 for i in items if i[0] == "dir")
    return f"{len(items)} file(s), {n_dir} .d run(s) to archive; queue: {qsrc}"


# ---------------------------------------------------------------------- HOW_TO_SUBMIT --
def _md_table(header, rows):
    esc = lambda s: str(s).replace("|", "\\|").replace("\n", " ")
    return "\n".join(["| " + " | ".join(header) + " |", "|" + "---|" * len(header)]
                     + ["| " + " | ".join(esc(c) for c in r) + " |" for r in rows])


def build_howto(f, rows, sdrf, prot, out):
    srec = f["srec"] or {}
    raw_rows = [r for r in rows if r["pride_file_type"] == "RAW"]
    raw_bytes = sum(r["size_bytes"] or 0 for r in raw_rows)
    dir_bytes = sum(r["size_bytes"] or 0 for r in raw_rows if r["kind"] == "dir")
    unsized = sum(1 for r in raw_rows if r["size_bytes"] is None)
    n_dirs = sum(1 for r in raw_rows if r["kind"] == "dir")
    eng = f"{srec.get('engine_label') or '?'} {srec.get('version') or ''}".strip()
    rel = lambda p: os.path.relpath(p, f["p"]["session_dir"])
    on_hive = any(r["location"].startswith("HIVE") for r in raw_rows)
    fill = sdrf.get("to_fill") or []
    today = datetime.date.today().isoformat()
    L = []
    w = L.append
    w(f"# How to deposit this dataset in a public repository")
    w("")
    w(f"*Session `{f['session']}` — package written {today} by the UC Davis Proteomics Core "
      f"pipeline skill. Everything referred to below is in this folder "
      f"(`output/DATA_SUBMISSION/`) unless a path says otherwise. Requirements are quoted from "
      f"the repositories' own documentation as read on 2026-09-24; check the linked pages if "
      f"you submit much later.*")
    w("")
    w("## What you are submitting")
    w("")
    w(_md_table(["Item", "This session"], [
        ["Raw files", f"{len(raw_rows)} ({_human(raw_bytes)}"
         + (f"; {unsized} not measurable from here" if unsized else "") + ")"
         + (f", {n_dirs} Bruker/Agilent `.d` folder(s) to compress" if n_dirs else "")],
        ["Where the raw files are", ", ".join(sorted({r['location'] for r in raw_rows})) or "?"],
        ["Instrument / acquisition", f"{f['instrument'] or '?'} / {f['acquisition'] or '?'}"],
        ["Search engine", eng],
        ["Organism", f"{f['organism'] or '?'} (taxid {f['taxid'] or '?'})"],
        ["Files to upload", f"{len(rows)} — listed in `files_to_upload.tsv`"],
    ]))
    w("")
    w("**Repository.** Use **PRIDE** (EMBL-EBI) unless you have a reason not to: it is the "
      "largest ProteomeXchange repository and PRIDE's guidelines name DIA-NN, FragPipe and "
      "Sage-style native outputs explicitly. **MassIVE** (UCSD) is the alternative — see the "
      "last section for when to pick it. Both issue a **PXD** accession, the identifier "
      "journals ask for.")
    w("")
    w("**Order of work:** (1) fill in what only you know → (2) register → (3) stage the files "
      "on HIVE → (4) upload with the PRIDE Submission Tool → (5) give the reviewer "
      "credentials to the journal → (6) make it public on acceptance → cite the PXD.")
    w("")
    w("## Step 1 — Fill in what only you know")
    w("")
    w(f"The skill filled in everything the analysis recorded. The rest is marked `{TO_FILL}`. "
      f"It is never guessed: an SDRF that says `female` or `liver` because a program guessed "
      f"would be worse than an empty cell.")
    w("")
    w(f"**`sdrf.tsv`** (sample metadata, SDRF-Proteomics {SDRF_VERSION}): "
      f"{sdrf.get('n_rows', 0)} row(s), {sdrf.get('n_to_fill', 0)} `{TO_FILL}` cell(s). Open it "
      "in Excel or any spreadsheet program, fill these columns, and save it as tab-separated "
      "text with the same name:")
    w("")
    if fill:
        w(_md_table(["Column", "What to enter", "Why it is not filled"], fill))
    else:
        w("Nothing — every column was filled from the session record (still check them).")
    w("")
    w("Rules from the SDRF-Proteomics specification: a value you do not know is "
      "`not available`; one that does not apply (e.g. organism part of a whole-organism "
      "culture) is `not applicable` -- where the column allows it (the table above says when "
      "it does not); both must be lower-case. Healthy/control samples have disease `normal`, "
      "not `control`. For a cell line (e.g. HeLa), use the `cell-lines` template columns "
      "(`characteristics[cell line]`, `characteristics[cellosaurus accession]`). The easiest "
      "check is the **PRIDE SDRF editor/validator**: "
      "https://www.ebi.ac.uk/pride/services/sdrf-editor/ (open `sdrf.tsv`, it flags every "
      "problem).")
    w("")
    w("How each filled column was derived:")
    w("")
    w(_md_table(["Column", "Value", "Where it came from"],
                [(c, v, why) for c, v, why in (sdrf.get("sources") or [])]))
    w("")
    w(f"**`protocols.txt`**: the two protocol texts. The *sample processing* block starts with "
      f"a `{TO_FILL}` for your sample preparation (lysis, reduction/alkylation, digestion, "
      f"clean-up) — the skill never saw that. Both blocks carry the Methods' "
      f"`[... — confirm]` tags ({prot.get('n_tags', '?')} in total): resolve each before "
      f"pasting. PRIDE's Submission Tool rejects a protocol of 50 characters or fewer, or of "
      f"5000 or more.")
    w("")
    w("You will also type in: a **title** (more than 30, fewer than 500 characters), a "
      "**description** like your abstract (more than 50, fewer than 5000), **keywords**, and "
      "the **lab head** (full name with a space, e-mail, affiliation).")
    w("")
    w("## Step 2 — Create a PRIDE account")
    w("")
    w("Register at https://www.ebi.ac.uk/pride/register (the submitter and the lab head can "
      "be different people). PRIDE notes that registration sends no confirmation e-mail; if "
      "you cannot log in 24 h after registering, write to pride-support@ebi.ac.uk.")
    w("")
    w("## Step 3 — Stage the files on HIVE (no download to your laptop)")
    w("")
    w("`prepare_upload.sbatch` was written for this session but **not run**. It packs each "
      "`.d` run into one `.tar.gz` (PRIDE: Bruker `.d` folders must be compressed, one run per "
      "archive, the whole folder unchanged), links or copies every other file into one staging "
      "folder, and writes `checksum.txt` (SHA-1, PRIDE's format) and `md5sums.txt`. It runs as "
      "a SLURM job because it reads every raw file — never on a login node.")
    w("")
    w("```")
    w(f"sbatch {os.path.join(out, 'prepare_upload.sbatch')}")
    w("```")
    w("")
    w(f"The staging folder is `{os.path.join(out, 'upload_staging')}` unless you set "
      f"`STAGE=/somewhere/else` in front of `sbatch`. "
      + (f"The `.d` archives need up to {_human(dir_bytes)} there (the size of the `.d` "
         f"folders){' plus the runs not measurable from here' if unsized else ''}; "
         if n_dirs else "")
      + "single-file raws are only linked unless you set `COPY_RAW=1`. The session zip never "
      "includes the staging folder. Re-run the job after editing `sdrf.tsv`: finished "
      "archives are kept, small files are re-copied and re-hashed.")
    if not on_hive and raw_rows:
        w("")
        w("> The raw paths were recorded on the machine that ran the analysis, not on HIVE. If "
          "that is a network share HIVE mounts, `bash scripts/hive_path.sh <path>` prints the "
          "HIVE path; edit the paths at the end of `prepare_upload.sbatch` to match. If the raw "
          "files are only on your computer, run it there instead: "
          "`RUN_HERE=1 bash prepare_upload.sbatch`.")
    w("")
    w("## Step 4 — Upload with the PRIDE Submission Tool, running on HIVE")
    w("")
    w("PRIDE's standard route is its desktop **PRIDE Submission Tool** (also called the PX "
      f"Submission Tool; current release {PX_TOOL_VERSION}). It is a Java program with a window, "
      "so run it on a HIVE **Open OnDemand desktop** — that runs as a SLURM job on a compute "
      "node next to the data, and keeps running if you close the browser:")
    w("")
    w("1. Go to https://ondemand.hive.hpc.ucdavis.edu → **Interactive Apps → Hive Desktop**. "
      "Ask for enough hours for the upload (PRIDE says Aspera \"can transfer terabytes within "
      "a day\"; if the session ends first, restart the tool — the upload resumes).")
    w("2. In the desktop, open a terminal and run:")
    w("")
    w("   ```")
    w("   cd ~ && curl -L -o px-tool.zip https://github.com/PRIDE-Archive/px-submission-tool/"
      "releases/latest/download/px-submission-tool-latest.zip")
    w("   unzip -q px-tool.zip && cd px-submission-tool-*/ && ./start.sh")
    w("   ```")
    w("")
    w("   (≈217 MB. HIVE compute nodes have Java 21, which the tool needs; checked 2026-09-24.)")
    w("3. **Log in** with your PRIDE account.")
    w("4. **Dataset description**: title, keywords, description, then paste the two blocks of "
      "`protocols.txt`. Pick the experiment type from the drop-down.")
    w("5. **Add files**: open the staging folder and add everything in it except "
      "`md5sums.txt`, `checksum.txt` and `.done/`. Check each file's type against the table "
      "below — the tool guesses from the extension; spectral libraries must be set by hand. "
      "Relate every RAW file to the SEARCH file(s): PRIDE requires every RAW file to be "
      "related to at least one SEARCH or RESULT file, and every SEARCH file to at least one "
      "RAW file.")
    w("6. **Checksums**: the tool offers to compute SHA-1s itself (PRIDE strongly recommends "
      "it; for very large datasets you may choose *No*).")
    w("7. **Metadata**: species, tissue (`not applicable` where tissue does not apply), "
      "instrument and modifications — pick the controlled-vocabulary terms that match "
      "`sdrf.tsv`; then the lab head.")
    w("8. **Upload**: Aspera is the default and the tool falls back to FTP by itself. From a "
      "HIVE compute node both PRIDE upload servers answered on 2026-09-24 "
      "(hx-fasp-1.ebi.ac.uk:33001, ftp-pride-private.ebi.ac.uk:21; Aspera's UDP data channel "
      "was not tested). If neither works, use Globus (below).")
    w("")
    w("**File types** (the tool's names; PRIDE's newer guideline text calls SEARCH "
      "\"ANALYSIS\" and RESULT \"STANDARD\"):")
    w("")
    w(_md_table(["File", "PRIDE type", "Needed", "Note"],
                [(r["upload_name"], r["pride_file_type"], r["requirement"],
                  r["before_upload"] or r["notes"]) for r in rows]))
    w("")
    w("This is a DIA/native-output submission: there is no mzIdentML/mzTab (RESULT) file, "
      "which PRIDE's guidelines say is normal for DIA — native DIA-NN/FragPipe/Sage output is "
      "the mandatory part. ProteomeXchange will list it as a \"partial\" submission; PRIDE "
      "states that this label says only that no standard-format file was included, not that "
      "anything is missing." if srec.get("engine") in ("diann", "fragpipe", "radiant", "alphadia")
      or f["acquisition"] == "DIA" else
      "No mzIdentML/mzTab (RESULT) file is produced by this pipeline; the native search output "
      "(SEARCH) is what PRIDE requires. ProteomeXchange will list the dataset as \"partial\", "
      "which PRIDE says reflects only the absence of a standard-format file.")
    w("")
    w("**Very large dataset, or Aspera/FTP blocked? Use Globus.** In your PRIDE profile, "
      "**New Submission** → request a Globus folder (you get its name by e-mail). You still "
      "need `submission.px` (from the Submission Tool: go through the wizard and click "
      "**Export summary** on the last panel, then close it without uploading) and "
      "`checksum.txt` (written by `prepare_upload.sbatch`). HIVE runs Globus: search for the "
      "collection *UC Davis Hive home* (your home directory); a lab (PI) share is visible only "
      "after HPC support exports it (ask them, with your PI in CC). Stage into a folder Globus "
      "can see (`STAGE=...`) and use `COPY_RAW=1`: Globus skips symbolic links inside a folder "
      "it transfers, which is how the script stages single-file raws by default. Transfer the "
      "staging folder plus `submission.px` to the *PRIDE Submissions* collection, then tell "
      "PRIDE the upload is done under **New Submission** in your profile.")
    w("")
    w("## Step 5 — After the upload: accession and reviewer access")
    w("")
    w("You get a submission reference `1-XXXXXXXX-X` at once (not the accession). A PRIDE "
      "curator validates the files; PRIDE says this can take up to five working days. Then an "
      "e-mail brings the **PXD accession** and a **reviewer account**. The dataset is private "
      "until you publish it: put the reviewer credentials in your cover letter or manuscript "
      "so editors and reviewers can log in at https://www.ebi.ac.uk/pride/login. You can edit "
      "metadata, and add or replace files (Submission Tool resubmission mode), while it is "
      "private; a leaked reviewer token can be reset from the dataset page.")
    w("")
    w("## Step 6 — Make it public, and cite it")
    w("")
    w("Once the paper is accepted, **you** release the data: log in, open the private dataset, "
      "click **Publish** and give the PubMed ID or DOI. (PRIDE keeps data private for up to two "
      "years before asking for a release date, and releases datasets it finds cited in "
      "published papers.)")
    w("")
    w("In the paper's data-availability statement, PRIDE asks for:")
    w("")
    w("> The mass spectrometry proteomics data have been deposited to the ProteomeXchange "
      "Consortium via the PRIDE [1] partner repository with the dataset identifier PXDxxxxxx.")
    w("")
    w("[1] Perez-Riverol Y, Bandla C, Kundu DJ, et al. The PRIDE database at 20 years: 2025 "
      "update. Nucleic Acids Res. 2025;53(D1):D543-D553. doi:10.1093/nar/gkae1011")
    w("")
    w("If the dataset page shows a DOI, you can add \"and DOI 10.6019/PXDxxxxxx\". Methods for "
      "the paper itself are in `output/methods.md` (and `methods.docx`), including the UC "
      "Davis instrument-grant acknowledgment.")
    w("")
    w("## Alternative: MassIVE")
    w("")
    w("Choose MassIVE (https://massive.ucsd.edu) if your lab or journal already uses it, or "
      "if you would rather not run a desktop program: MassIVE takes plain **FTP** uploads to "
      "your own account and a web form, so everything can be done from a HIVE terminal. In the "
      "ProteomeXchange guidelines' table, MassIVE is also the repository listing both partial "
      "and complete DIA submissions (PRIDE: partial only).")
    w("")
    w("1. Register an account on the MassIVE site (the link is under the login box).")
    w("2. Stage the files with `prepare_upload.sbatch` as above.")
    w("3. Upload from HIVE with `lftp` (installed there), from an Open OnDemand desktop "
      "terminal or an interactive SLURM session — not the login node — using **FTP with "
      "explicit TLS** on port 21, your MassIVE user name and password:")
    w("")
    w("   ```")
    w(f"   lftp -u <massive-user> -e 'set ftp:ssl-force true; set ftp:ssl-protect-data true; "
      f"mirror -R --parallel=4 {os.path.join(out, 'upload_staging')} {f['session']}; bye' "
      f"massive-ftp.ucsd.edu")
    w("   ```")
    w("")
    w("   (Add `-L` to `mirror` if you staged single-file raws as links.) MassIVE's "
      "documentation does not say whether it wants Bruker `.d` folders compressed; the "
      "archives are what the script makes — ask ccms-web@cs.ucsd.edu if unsure.")
    w("4. On the site: **Submit Data** → workflow **MassIVE Dataset Submission**. Required: "
      "species, instrument, post-translational modifications (controlled terms — take them "
      "from `sdrf.tsv`), at least one keyword, the principal investigator, and a non-empty "
      "**dataset password**.")
    w("5. Assign the uploaded files to categories — the `massive_category` column of "
      "`files_to_upload.tsv` (MassIVE has no SDRF category; `sdrf.tsv` goes under "
      "Supplementary Files, and `methods.docx` can go under Methods and Protocols).")
    w("6. Tick **Submit to ProteomeXchange**, then **Submit**. You get an **MSV** accession.")
    w("7. Reviewers log in as `MSV…_reviewer` with the dataset password (FTP: the MSV "
      "accession as user name). MassIVE's documentation says the ProteomeXchange announcement "
      "is sent when you click **Make Public**; whether the PXD is shown before that is not "
      "stated there — check the dataset page, and cite the PXD (and MSV) once you have it.")
    w("")
    w("## Where these instructions come from")
    w("")
    w("PRIDE: https://www.ebi.ac.uk/pride/markdownpage/submitdatapage, "
      "…/pridesubmissiontool, …/datasubmissionguidelines (v2.2.0), …/globus, …/checksum, "
      "…/sdrf, …/citationpage · PRIDE Submission Tool 2.11.6 (its built-in validator) · "
      "ProteomeXchange guidelines v3.0.1 https://www.proteomexchange.org/docs/guidelines_px.pdf "
      "· MassIVE https://ccms-ucsd.github.io/MassIVEDocumentation/ · SDRF-Proteomics v1.1.0 "
      "https://github.com/bigbio/proteomics-sample-metadata and templates "
      "https://github.com/bigbio/sdrf-templates · HIVE https://docs.hpc.ucdavis.edu/ (data "
      "transfer, Open OnDemand). The skill's `references/deposit.md` quotes each requirement "
      "with the date it was read.")
    w("")
    return "\n".join(L) + "\n"


def write_howto(f, out, rows, sdrf, prot):
    path = os.path.join(out, "HOW_TO_SUBMIT.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(build_howto(f, rows, sdrf, prot, out))
    return None


# ------------------------------------------------------------- minimal Markdown -> HTML --
def _inline(s):
    s = html.escape(s, quote=False)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<![\w*])\*([^*\s][^*]*)\*(?![\w*])", r"<em>\1</em>", s)
    s = re.sub(r"\[([^\]]+)\]\((https?://[^)\s]+)\)", r'<a href="\2">\1</a>', s)
    s = re.sub(r"(?<![\"'>=])(https?://[^\s<)]+[^\s<).,;])", r'<a href="\1">\1</a>', s)
    return s


def md_to_html(md, title):
    """Enough Markdown for HOW_TO_SUBMIT.md: headings, paragraphs, lists, tables, code blocks,
    block quotes. Stdlib only, so the .html opens anywhere by double-click."""
    out, para, lst, table, code = [], [], None, [], None
    lines = md.splitlines()

    def flush():
        nonlocal para, lst, table
        if para:
            out.append("<p>" + _inline(" ".join(para)) + "</p>")
            para = []
        if lst:
            tag, items = lst
            out.append(f"<{tag}>" + "".join(f"<li>{i}</li>" for i in items)
                       + f"</{tag.split()[0]}>")
            lst = None
        if table:
            rows = [r for r in table if not re.match(r"^\|(\s*-+\s*\|)+$", r.strip())]
            cells = [[c.strip().replace("\\|", "|") for c in
                      re.split(r"(?<!\\)\|", r.strip().strip("|"))] for r in rows]
            if cells:
                head = "".join(f"<th>{_inline(c)}</th>" for c in cells[0])
                body = "".join("<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in r) + "</tr>"
                               for r in cells[1:])
                out.append(f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>")
            table = []

    for ln in lines:
        s = ln.strip()
        if code is not None:
            if s.startswith("```"):
                out.append("<pre><code>" + html.escape("\n".join(code)) + "</code></pre>")
                code = None
            else:
                code.append(ln.strip() if ln.startswith("   ") else ln)
            continue
        if s.startswith("```"):
            flush()
            code = []
            continue
        if not s:
            flush()
            continue
        m = re.match(r"^(#{1,4})\s+(.*)$", s)
        if m:
            flush()
            n = len(m.group(1))
            out.append(f"<h{n}>{_inline(m.group(2))}</h{n}>")
            continue
        if s.startswith("|"):
            if para or lst:
                flush()
            table.append(s)
            continue
        if s.startswith(">"):
            flush()
            out.append("<blockquote>" + _inline(s.lstrip("> ")) + "</blockquote>")
            continue
        m = re.match(r"^(\d+)\.\s+(.*)$", s) or re.match(r"^[-*]\s+(.*)$", s)
        if m and not ln.startswith("   "):
            tag = "ol" if s[0].isdigit() else "ul"
            if table or para:
                flush()
            if not lst or lst[0].split()[0] != tag:
                flush()
                # a numbered list interrupted by a code block continues at its own number
                start = int(m.group(1)) if tag == "ol" and m.group(1) != "1" else None
                lst = (f'ol start="{start}"' if start else tag, [])
            lst[1].append(_inline(m.groups()[-1]))
            continue
        if lst and ln.startswith("   "):
            lst[1][-1] += " " + _inline(s)
            continue
        if table:
            flush()
        para.append(s)
    flush()
    css = ("body{font:15px/1.55 -apple-system,Segoe UI,Helvetica,Arial,sans-serif;max-width:"
           "980px;margin:2em auto;padding:0 16px;color:#1d2330;background:#fff}"
           "h1,h2{color:#1b3a5c}h2{border-bottom:1px solid #d8dee8;padding-bottom:.2em;"
           "margin-top:1.8em}code{background:#f1f3f7;padding:1px 4px;border-radius:3px}"
           "pre{background:#f1f3f7;padding:10px;overflow-x:auto}pre code{background:none}"
           "table{border-collapse:collapse;margin:.6em 0;font-size:13px;display:block;"
           "overflow-x:auto}th,td{border:1px solid #d8dee8;padding:4px 8px;text-align:left;"
           "vertical-align:top}th{background:#eef2f7}blockquote{border-left:4px solid #9fb4cc;"
           "margin:.8em 0;padding:.2em 1em;background:#f7f9fc}")
    return (f"<!DOCTYPE html>\n<html lang=\"en\"><head><meta charset=\"utf-8\">"
            f"<meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">"
            f"<title>{html.escape(title)}</title><style>{css}</style></head><body>\n"
            + "\n".join(out) + "\n</body></html>\n")


def write_html(out):
    src = os.path.join(out, "HOW_TO_SUBMIT.md")
    if not os.path.isfile(src):
        raise Skip("HOW_TO_SUBMIT.md was not written (see above)")
    with open(src, encoding="utf-8") as fh:
        page = md_to_html(fh.read(), "How to deposit this dataset")
    with open(os.path.join(out, "HOW_TO_SUBMIT.html"), "w", encoding="utf-8") as fh:
        fh.write(page)
    return None


# ------------------------------------------------------------------------- methods --
def required_sections(f):
    need = ["Liquid chromatography", "Mass spectrometry", "Sequence database"]
    if f["srec"].get("engine") or f["params"] or f["search_prov_path"]:
        need.append("Database search")
    if f["de_prov"]:
        need.append("Differential expression")
    return need + ["Acknowledgments"]


def methods_command(f, out):
    """make_methods.py's command line for this session, or Skip with the reason it can't run."""
    if not f["raws"]:
        raise Skip("no raw file list (input/raw_files.txt): make_methods.py needs the raw "
                   "files, or at least their names")
    cmd = [sys.executable, os.path.join(HERE, "make_methods.py"), "--raw", *f["raws"],
           "--out", out]
    for flag, val in (("--fasta-meta", f["fasta_meta_path"]),
                      ("--de-dir", f["p"]["de_dir"] if f["de_prov"] else None),
                      ("--params", f["params"]), ("--search-prov", f["search_prov_path"]),
                      ("--workflow-manifest", f["wf_path"]),
                      ("--instrument", next(iter(f["wf"].get("instruments") or []), None)
                       or (f["run_manifest"].get("query") or {}).get("instrument")),
                      ("--acquisition", f["acquisition"])):
        if val:
            cmd += [flag, str(val)]
    return cmd


def _run(cmd, what):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    except (OSError, subprocess.TimeoutExpired) as e:
        raise Skip(f"{what} could not run: {e}")
    if r.returncode != 0:
        msg = (r.stderr or r.stdout or "").strip().splitlines()
        raise Skip(f"{what} failed: {msg[-1] if msg else f'exit {r.returncode}'}")
    return r


def ensure_methods(session_dir, man):
    """Make sure the session carries publication Methods (output/methods.md + .docx). A
    methods.md already there is never overwritten -- it may be hand-polished; if it lacks a
    required section, a complete draft is written beside it. Returns the file the protocols
    should be drawn from (or None)."""
    f = gather(session_dir)
    p = f["p"]
    chosen = {}

    def _md():
        need = required_sections(f)
        md = p["methods_md"]
        if os.path.isfile(md):
            with open(md, encoding="utf-8") as fh:
                have = md_sections(fh.read())
            missing = [s for s in need if s not in have]
            chosen["path"] = md
            if not missing:
                return "present -- kept as written"
            draft = os.path.join(p["output_dir"], "methods_complete_draft.md")
            try:
                _run(methods_command(f, draft), "make_methods.py")
            except Skip as e:
                return (f"present, kept as written, but it lacks: {', '.join(missing)}; a "
                        f"complete draft could not be written ({e})")
            chosen["path"] = draft
            return (f"present, kept as written, but it lacks: {', '.join(missing)} -- a "
                    f"complete draft is in output/methods_complete_draft.md")
        _run(methods_command(f, md), "make_methods.py")
        chosen["path"] = md
        with open(md, encoding="utf-8") as fh:
            have = md_sections(fh.read())
        missing = [s for s in need if s not in have]
        return "generated by make_methods.py" + (
            f"; no record to write: {', '.join(missing)}" if missing else "")

    man.section("Publication methods (output/methods.md)", _md)

    def _docx():
        src = chosen.get("path")
        if not src:
            raise Skip("no methods .md to convert (see the line above)")
        out = os.path.splitext(src)[0] + ".docx"
        if os.path.isfile(out) and os.path.getmtime(out) >= os.path.getmtime(src):
            return f"{os.path.basename(out)} up to date -- kept"
        _run([sys.executable, os.path.join(HERE, "to_docx.py"), "--in", src, "--out", out],
             "to_docx.py")
        return f"wrote {os.path.relpath(out, p['session_dir'])}"

    man.section("Publication methods, Word (.docx)", _docx)
    return chosen.get("path")


# --------------------------------------------------------------------------- build --
def build(session_dir, man, methods_md=None):
    """Write output/DATA_SUBMISSION/. Each part is its own manifest section: one that fails is
    recorded as [SKIPPED] with its reason, and the rest are still written."""
    f = gather(session_dir)
    out = f["p"]["deposit_dir"]
    os.makedirs(out, exist_ok=True)
    rows, sdrf, prot = [], {}, {}

    def _plan():
        rows.extend(plan_uploads(f))
        n = hash_small(rows)
        return f"{len(rows)} file(s); {n} small file(s) hashed here"

    def _raws():
        raws = [r for r in rows if r["pride_file_type"] == "RAW"]
        if not raws:
            raise Skip("no raw files recorded (input/raw_files.txt): every submission needs "
                       "them -- add them to files_to_upload.tsv and the SDRF by hand")
        far = sum(1 for r in raws if r["size_bytes"] is None)
        return f"{len(raws)} raw file(s)" + (f"; {far} not reachable from here" if far else "")

    def _need_plan():
        if not rows:
            raise Skip("no upload plan (see the upload-plan line above)")

    def _list():
        _need_plan()
        write_upload_list(out, rows)

    def _prep():
        _need_plan()
        return write_prep_script(f, out, rows)

    planned = man.section("Deposit: upload plan", _plan)
    man.section("Deposit: raw files in the upload plan", _raws)
    man.section("Deposit: sdrf.tsv (SDRF-Proteomics v1.1.0)", write_sdrf, f, out, sdrf)
    man.section("Deposit: protocols.txt", write_protocols, out, methods_md
                or (f["p"]["methods_md"] if os.path.isfile(f["p"]["methods_md"]) else None),
                prot)
    man.section("Deposit: files_to_upload.tsv", _list)
    man.section("Deposit: prepare_upload.sbatch (written, not run)", _prep)
    man.section("Deposit: HOW_TO_SUBMIT.md", write_howto, f, out, rows, sdrf, prot)
    man.section("Deposit: HOW_TO_SUBMIT.html", write_html, out)
    return {"dir": out, "planned": planned, "n_files": len(rows),
            "sdrf_to_fill": sdrf.get("n_to_fill"), "protocol_tags": prot.get("n_tags")}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--session", required=True, help="the session directory")
    ap.add_argument("--skip-methods", action="store_true",
                    help="do not create output/methods.md when it is missing")
    a = ap.parse_args()
    if not os.path.isdir(a.session):
        sys.exit(f"session dir not found: {a.session}")
    man = Manifest()
    md = None if a.skip_methods else ensure_methods(a.session, man)
    res = build(a.session, man, md)
    mpath = man.write(paths_for(a.session)["manifest_txt"], "Session export manifest")
    res.update(manifest=mpath, skipped=man.n_skipped)
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
