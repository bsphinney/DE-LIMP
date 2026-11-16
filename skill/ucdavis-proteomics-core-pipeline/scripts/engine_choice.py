#!/usr/bin/env python3
"""
engine_choice.py -- which search engine runs, why, and whether FragPipe can run here. The ONE
route table and the ONE place that knows where the Core's licensed FragPipe lives.

  * ROUTES / de_method() / choose(): the engine for a data type (resolve_defaults.py step 4
    imports them; nothing else words the choice). timsTOF DDA-PASEF -> FragPipe when it can run
    here and its licence-gated tools are installed for this account, else Sage, with the reason
    (Rulings 55/56). Orbitrap DDA with a non-specific digest (--enzyme none, or pepsin: no
    in-silico rule) -> FragPipe's Nonspecific-peptidome workflow, else Sage (Ruling 120;
    routes_for). Other Orbitrap and other DDA -> Sage; FragPipe on request (Ruling 67). DIA ->
    DIA-NN.
  * fragpipe_tools(): what tools.json `fragpipe_tools` says -- the tools folder FragPipe is given
    (--config-tools-folder), the MSFragger/IonQuant/diaTracer jars in it, Philosopher and the
    DIA-NN FragPipe bundles (both from the RUNNING binary's own install: FragPipe takes
    Philosopher from there whatever the tools folder, C0 2026-11-04), the Java, and `licensed`.
    acquire_tools.sh writes it; run_search.py reads the folder from tools.json only.
  * CORE_FRAGPIPE_INSTALL: the Core's FragPipe 24.0 on HIVE. diatracer_parallel.py,
    fragpipe_bootstrap.py and acquire_tools.sh (through this script's CLI) read it from here.

The selection record (resolve_defaults.py puts it in workflow.manifest.json `engine.selection`)
is client-safe: one sentence with no version, path or licence word. The reason code, `why` and
the checks go only into the staff-only record beside the manifest (staff.STAFF_SUFFIX).

  python3 engine_choice.py fragpipe-tools [--binary <bin/fragpipe>] [--pin 24.0]   # tools.json object
  python3 engine_choice.py core-binary [--pin 24.0]      # the Core's bin/fragpipe, or nothing
  python3 engine_choice.py java [--line "<java -version line>"] [--tsv]
  python3 engine_choice.py choose --acquisition DDA --instrument "timsTOF Pro 2" \\
      [--env env.json] [--tools tools.json] [--engine sage] [--ids-only]
Stdlib only (acquire_tools.sh and detect_env.sh run it with the system python3).
"""
import argparse
import json
import os
import re
import subprocess
import sys

# ---------------------------------------------------------------- the Core install -------------
# The Core's licensed FragPipe 24.0 (MSFragger 4.4.1, IonQuant 1.11.20, diaTracer 2.2.1 at depth 1
# of tools/; Philosopher v5.1.3-RC9; bundled DIA-NN 1.8.2 beta 8). Read on HIVE by C0, 2026-11-04.
# The only place this path is written: everything else calls core_install()/core_tools().
CORE_FRAGPIPE_INSTALL = "/quobyte/proteomics-grp/fragpipe24/fragpipe-24.0"
# Another site (or a test) points the same logic at its own install.
CORE_ENV = "FRAGPIPE_CORE_INSTALL"
# What a user sets to point FragPipe at their own licensed jars. Read here only (at acquire time);
# the search reads tools.json `fragpipe_tools.folder` (review M4).
TOOLS_ENV = "FRAGPIPE_TOOLS_FOLDER"
# Where the tools folder came from, as a DELIVERED record says it (search_provenance.json is in
# the client's delivery): never the variable's name -- client records carry no licence, jar or
# FRAGPIPE_TOOLS_FOLDER wording (review M1, Ruling 107.1). tools.json, a staff file, keeps
# TOOLS_ENV (Ruling 56).
DELIVERED_SOURCE = {TOOLS_ENV: "account_tools_folder"}


def delivered_source(source):
    """fragpipe_tools `source` as a delivered record carries it (DELIVERED_SOURCE)."""
    return DELIVERED_SOURCE.get(source, source)

# The Python FragPipe hands EasyPQP (--config-python) for a spectral-library workflow: a venv kept
# beside the Core install (the shared settings cache's bin-python on HIVE, 2026-11-04). A user's own
# is FRAGPIPE_PYTHON, or the bin-python of their FragPipe settings cache (fragpipe_bootstrap.py).
CORE_EASYPQP_VENV = "easypqp_venv2"
PYTHON_ENV = "FRAGPIPE_PYTHON"

# FragPipe 24.0 is built for Java 11 (its build.gradle); HIVE's nodes run OpenJDK 25.0.4.1 (H7, C0).
FRAGPIPE_JAVA_MIN = 11

# FragPipe finds its jars by these names, case-insensitively, keeping the highest version
# (TabConfig.java, tag 24.0). It walks the whole folder; the skill looks two levels down, which
# covers the Core's layouts (tools/MSFragger-4.4.1.jar, tools_licensed/diatracer-2.2.1/
# diaTracer-2.2.1.jar) without a deep walk of a shared folder on a login node.
JARS = {"msfragger": "MSFragger", "ionquant": "IonQuant", "diatracer": "diaTracer"}
JAR_DEPTH = 2
LICENSED_NEEDS = ("msfragger", "ionquant")


def core_install():
    return os.environ.get(CORE_ENV) or CORE_FRAGPIPE_INSTALL


def core_tools():
    return os.path.join(core_install(), "tools")


def core_binary():
    return os.path.join(core_install(), "bin", "fragpipe")


def _vkey(v):
    """'4.4.1' -> (4, 4, 1): numbers compared as numbers, so 4.10 > 4.9."""
    return tuple(int(x) for x in re.findall(r"\d+", v or ""))


def _walk(folder, depth):
    """(path, name) of every file at depth 1..`depth` under `folder`; unreadable subfolders skipped."""
    try:
        entries = list(os.scandir(folder))
    except OSError:
        return
    for e in entries:
        try:
            if e.is_file():
                yield e.path, e.name
            elif depth > 1 and e.is_dir():
                yield from _walk(e.path, depth - 1)
        except OSError:
            continue


def find_jar(folder, tool, depth=JAR_DEPTH):
    """The `tool` jar ("msfragger" | "ionquant" | "diatracer") under `folder`, two levels down,
    case-insensitive, highest version (a shallower copy wins a tie) -> {"path", "version"} | None."""
    if not folder or tool not in JARS:
        return None
    pat = re.compile(rf"^{tool}-(.+)\.jar$", re.IGNORECASE)
    best = None
    for path, name in _walk(folder, depth):
        m = pat.match(name)
        if not m:
            continue
        key = (_vkey(m.group(1)), -path.count(os.sep))
        if best is None or key > best[0]:
            best = (key, {"path": path, "version": m.group(1)})
    return best[1] if best else None


def unusable(folder):
    """Why `folder` cannot be searched (absent, not a folder, not readable by this account), or None."""
    if not os.path.exists(folder):
        # a folder inside one this account cannot open looks absent; say which it is
        up = os.path.dirname(os.path.abspath(folder))
        while not os.path.exists(up) and os.path.dirname(up) != up:
            up = os.path.dirname(up)
        if os.path.isdir(up) and not os.access(up, os.R_OK | os.X_OK):
            return f"not readable by this account ({up} is closed to it)"
        return "does not exist"
    if not os.path.isdir(folder):
        return "is not a folder"
    try:
        os.listdir(folder)
    except OSError as e:
        return f"not readable by this account ({e.strerror or e})"
    return None


def _release_of(binary):
    """<...>/fragpipe-24.0/bin/fragpipe -> "24.0" (the install folder's name), else None."""
    m = re.match(r"^fragpipe-v?(\d+(?:\.\d+)+)$", os.path.basename(install_of(binary) or ""),
                 re.IGNORECASE)
    return m.group(1) if m else None


def install_of(binary):
    """<install>/bin/fragpipe -> <install>."""
    if not binary:
        return None
    return os.path.dirname(os.path.dirname(os.path.abspath(binary)))


def _same(a, b):
    return bool(a and b) and os.path.normpath(os.path.abspath(a)) == os.path.normpath(os.path.abspath(b))


def philosopher_of(install):
    """The Philosopher an install runs: <install>/tools/Philosopher/philosopher-v5.1.3-RC9 ->
    {"path", "version": "5.1.3-RC9"} | None. FragPipe always runs its OWN install's Philosopher."""
    if not install:
        return None
    pat = re.compile(r"^philosopher-v?(\d[^/]*?)(\.exe)?$", re.IGNORECASE)
    best = None
    for path, name in _walk(os.path.join(install, "tools"), 2):
        m = pat.match(name)
        if m and (best is None or _vkey(m.group(1)) > best[0]):
            best = (_vkey(m.group(1)), {"path": path, "version": m.group(1)})
    return best[1] if best else None


def bundled_diann_of(install):
    """The DIA-NN an install bundles (FragPipe's default --config-diann):
    <install>/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8 -> {"path", "version": "1.8.2 beta 8"}.
    The version is the bundle FOLDER's: the file is named 1.8.1.8 but prints "DIA-NN 1.8.2 beta 8"."""
    if not install:
        return None
    base = os.path.join(install, "tools", "diann")
    try:
        vers = sorted(os.listdir(base), key=_vkey)
    except OSError:
        return None
    for v in reversed(vers):
        d = os.path.join(base, v, "linux")
        try:
            names = sorted(n for n in os.listdir(d) if n.lower().startswith("diann")
                           and not n.lower().endswith((".txt", ".so")))
        except OSError:
            continue
        for n in names:
            p = os.path.join(d, n)
            if os.path.isfile(p):
                return {"path": p, "version": v.replace("_", " ")}
    return None


def core_easypqp_python():
    return os.path.join(os.path.dirname(core_install()), CORE_EASYPQP_VENV, "bin", "python")


def _cache_python(release):
    """bin-python from this account's FragPipe settings cache ($XDG_CONFIG_HOME, else ~/.config:
    FragPipe/fragpipe/<release>/fragpipe-ui.cache -- what fragpipe_bootstrap.py writes)."""
    if not release:
        return None
    base = os.environ.get("XDG_CONFIG_HOME") or os.path.join(os.path.expanduser("~"), ".config")
    try:
        with open(os.path.join(base, "FragPipe", "fragpipe", release, "fragpipe-ui.cache")) as fh:
            for line in fh:
                k, _, v = line.partition("=")
                if k.strip() == "fragpipe-config.bin-python" and v.strip():
                    return v.strip()
    except OSError:
        return None
    return None


def easypqp_python(release=None):
    """The Python FragPipe hands EasyPQP for a spectral-library workflow -> {"path", "version",
    "easypqp", "source", "resolves_to", "runnable", "why"} | None. Read from files only -- the venv's
    pyvenv.cfg and site-packages/easypqp-*.dist-info; nothing is imported (a login node). Order:
    $FRAGPIPE_PYTHON, this account's settings-cache bin-python, the Core's venv. `runnable` is for
    THIS account: the Core's venv links its python into one person's home (2026-11-04)."""
    cands = [(os.environ.get(PYTHON_ENV), PYTHON_ENV), (_cache_python(release), "fragpipe settings cache"),
             (core_easypqp_python(), "core_install")]
    for path, source in cands:
        if not path or not os.path.lexists(path):
            continue
        root = os.path.dirname(os.path.dirname(os.path.abspath(path)))
        version = None
        try:
            with open(os.path.join(root, "pyvenv.cfg")) as fh:
                for line in fh:
                    k, _, v = line.partition("=")
                    if k.strip() in ("version", "version_info") and v.strip():
                        version = v.strip()
                        break
        except OSError:
            pass
        easypqp = None
        try:
            for py in sorted(os.listdir(os.path.join(root, "lib"))):
                sp = os.path.join(root, "lib", py, "site-packages")
                for n in sorted(os.listdir(sp)) if os.path.isdir(sp) else []:
                    m = re.match(r"^easypqp-(.+)\.dist-info$", n, re.IGNORECASE)
                    if m:
                        easypqp = m.group(1)
                if version is None and re.match(r"^python\d", py):
                    version = py[len("python"):]
        except OSError:
            pass
        real = os.path.realpath(path)
        runnable = os.path.isfile(real) and os.access(real, os.X_OK)
        why = None
        if not runnable:
            why = f"{path} resolves to {real}, which this account cannot run"
        elif easypqp is None:
            why = f"no easypqp installed beside {path}"
        return {"path": path, "version": version, "easypqp": easypqp, "source": source,
                "resolves_to": real, "runnable": runnable, "why": why}
    return None


def core_fragpipe_binary(pin=None):
    """(the Core's bin/fragpipe, None) when this account can run it and it is the pinned release;
    else (None, why). The pin is matched against the install folder's name (fragpipe-24.0)."""
    b = core_binary()
    inst = install_of(b)
    why = unusable(inst)
    if why:
        return None, f"the Core install {inst} {why}"
    if not (os.path.isfile(b) and os.access(b, os.X_OK)):
        return None, f"the Core install has no runnable {b}"
    if pin and os.path.basename(inst).lower() != f"fragpipe-{pin}".lower():
        return None, (f"the Core install {inst} is not FragPipe {pin} (the pinned release) -- "
                      f"not substituted")
    return b, None


# ---------------------------------------------------------------- Java -------------------------
def java_major(line):
    """The Java major version in a `java -version` line, quoted or not (detect_env.sh strips the
    quotes): 'openjdk version "1.8.0_392"' -> 8, 'openjdk version 25.0.4.1 2025-10-21' -> 25.
    None when the line names no version (macOS's 'Unable to locate a Java Runtime')."""
    m = re.search(r'\b(?:version|openjdk|java)\s+"?(\d+)(?:\.(\d+))?', line or "", re.IGNORECASE)
    if not m:
        return None
    major = int(m.group(1))
    return int(m.group(2)) if major == 1 and m.group(2) else major


def java_record(line):
    """{"version", "major", "ok", "why"} for a `java -version` first line (None = no java)."""
    line = (line or "").strip() or None
    major = java_major(line)
    if line is None:
        return {"version": None, "major": None, "ok": False,
                "why": f"no java on PATH (FragPipe needs Java {FRAGPIPE_JAVA_MIN} or newer)"}
    if major is None:
        return {"version": line, "major": None, "ok": None,
                "why": f"could not read a Java version from {line!r}"}
    ok = major >= FRAGPIPE_JAVA_MIN
    return {"version": line, "major": major, "ok": ok,
            "why": None if ok else (f"Java {major} found ({line}); FragPipe needs Java "
                                    f"{FRAGPIPE_JAVA_MIN} or newer")}


def java_status(java_cmd="java"):
    """Run `java -version` (it prints on stderr) -> java_record of its first line."""
    try:
        r = subprocess.run([java_cmd, "-version"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return java_record(None)
    out = (r.stderr or "") + (r.stdout or "")
    first = next((ln for ln in out.splitlines() if ln.strip()), "")
    return java_record(first.replace('"', "") if first else None)


# ---------------------------------------------------------------- tools.json fragpipe_tools ----
def fragpipe_tools(binary=None, env_folder=None, pin=None, java=None):
    """The tools.json `fragpipe_tools` object (see INTERFACE / the module docstring).

    binary      the bin/fragpipe tools.json runs; None = the Core's when this account can run it
    env_folder  $FRAGPIPE_TOOLS_FOLDER (the caller reads the environment; None = not set)
    pin         the pinned FragPipe release, for the Core binary
    java        a java_record (default: run `java -version` here)"""
    skipped = []
    if not binary:
        binary, why = core_fragpipe_binary(pin)
        if why and os.path.exists(core_install()):
            skipped.append({"folder": core_install(), "source": "core_install", "why": why})
    inst = install_of(binary)
    cands = []
    if env_folder:
        cands.append((env_folder, TOOLS_ENV))
    if inst:
        cands.append((os.path.join(inst, "tools"),
                      "core_install" if _same(inst, core_install()) else "install"))
    if not any(_same(f, core_tools()) for f, _ in cands):
        cands.append((core_tools(), "core_install"))
    chosen = None
    for folder, source in cands:
        why = unusable(folder)
        if why:
            # an absent Core path is simply not this machine's; anything else is worth saying
            if not (source == "core_install" and why == "does not exist"):
                skipped.append({"folder": folder, "source": source, "why": why})
            continue
        jars = {t: find_jar(folder, t) for t in JARS}
        missing = [JARS[t] for t in LICENSED_NEEDS if not jars[t]]
        if not missing:
            chosen = (folder, source, jars)
            break
        skipped.append({"folder": folder, "source": source,
                        "why": f"no {' or '.join(missing)} jar within {JAR_DEPTH} levels"})
    # Philosopher and the bundled DIA-NN come from the RUNNING binary's install; with no binary,
    # from the Core install when its tools folder is the one chosen
    runs = inst or (core_install() if chosen and chosen[1] == "core_install" else None)
    java = java if java is not None else java_status()
    licensed = chosen is not None
    binary_ok = bool(binary) and os.path.isfile(binary) and os.access(binary, os.X_OK)
    ready = bool(licensed and binary_ok and java.get("ok"))
    if licensed:
        folder, source, jars = chosen
        found = ", ".join(f"{JARS[t]} {jars[t]['version']}" for t in JARS if jars[t])
        why = f"{found} in {folder} ({source})"
        if not binary_ok:
            why += "; but no runnable FragPipe binary"
        if java.get("ok") is False:
            why += f"; but {java['why']}"
    else:
        tried = "; ".join(f"{s['folder']}: {s['why']}" for s in skipped) or "no tools folder to look in"
        why = (f"MSFragger and IonQuant (licence-gated) were not found for this account ({tried}). "
               f"Set {TOOLS_ENV} to the folder holding your licensed jars, or accept the licence "
               f"in the FragPipe GUI once, then re-run acquire_tools.sh")
    return {
        "binary": binary or None,
        "folder": chosen[0] if chosen else None,
        "source": chosen[1] if chosen else None,
        "jars": chosen[2] if chosen else {t: None for t in JARS},
        "philosopher": philosopher_of(runs),
        "diann": bundled_diann_of(runs),
        # only a spectral-library workflow needs it (LFQ-MBR / Basic-Search do not, C0)
        "python": easypqp_python(pin or _release_of(binary)),
        "java": {k: java.get(k) for k in ("version", "major", "ok")},
        "licensed": licensed,
        "ready": ready,
        "why": why,
        "skipped": skipped,
    }


def fragpipe_ready(tools_json):
    """{"ready": True|False|None, "why"} from a tools.json dict's `fragpipe_tools` (None: no such
    key -- a tools.json from before 2.15 -- or no tools.json at all)."""
    ft = (tools_json or {}).get("fragpipe_tools") if isinstance(tools_json, dict) else None
    if not isinstance(ft, dict):
        return {"ready": None, "why": "tools.json has no fragpipe_tools record (written by "
                                      "acquire_tools.sh from 2.15); re-run acquire_tools.sh"}
    return {"ready": ft.get("ready"), "why": ft.get("why")}


# ---------------------------------------------------------------- the route table -------------
# Which engines are sensible for a data type, best first: [0] is the default. Instrument family
# from resolve_defaults.family(estimate_params.classify_instrument(...)).
ROUTES = {
    ("DIA", "timstof"): ["diann", "fragpipe"],
    ("DIA", "orbitrap"): ["diann", "fragpipe", "radiant"],
    ("DIA", "other"): ["diann"],
    ("DDA", "timstof"): ["fragpipe", "sage"],     # Ruling 55: FragPipe, else Sage (the fallback)
    ("DDA", "orbitrap"): ["sage"],                # Rulings 67/99: 2.14.2's route; FragPipe on request
    ("DDA", "other"): ["sage"],
}

# Ruling 120: a NON-SPECIFIC digest (cleavage_rules.enzyme_rule's specificity "none": --enzyme
# none, or pepsin, which has no in-silico rule) on Orbitrap DDA goes to FragPipe's
# Nonspecific-peptidome workflow, Sage the fallback (and reachable on request). Evidence: the same
# 30-run Q Exactive peptidome gave FragPipe 16,203 peptides and Sage 179 (paper: 14,270); Sage's
# +/-10 ppm window cut off an uncalibrated +8.5 ppm offset (sage_lfq_check.precursor_window_check
# now measures it). FragPipe's workflow does not calibrate either (calibrate_mass=0): its +/-20 ppm
# window and MSBooster/Percolator rescoring are the difference. Sage with a centred window gave
# 4,034 PSMs at 1% against 882 (8 runs), so the window is a large part; the rest is not measured. Only this row moves: a tryptic Orbitrap search and other DDA stay on Sage.
NONSPECIFIC_ROUTES = {("DDA", "orbitrap"): ["fragpipe", "sage"]}

# What the data type is called in the client sentence.
DATA_TYPE = {
    ("DIA", "timstof"): "timsTOF dia-PASEF data", ("DIA", "orbitrap"): "Orbitrap DIA data",
    ("DIA", "other"): "DIA data", ("DDA", "timstof"): "timsTOF DDA-PASEF data",
    ("DDA", "orbitrap"): "Orbitrap DDA data", ("DDA", "other"): "DDA data",
}

# ...and a non-specific row's (NONSPECIFIC_ROUTES)
NONSPECIFIC_DATA_TYPE = {("DDA", "orbitrap"): "non-specific Orbitrap DDA data"}

NAMES = {"diann": "DIA-NN", "sage": "Sage", "fragpipe": "FragPipe", "radiant": "Radiant"}


def non_specific(enzyme):
    """True when --enzyme `enzyme` (fetch_fasta's words) is searched with no in-silico rule:
    cleavage_rules.enzyme_rule's specificity "none" (the ONE reading, which make_presets' workflow
    choice uses too). None/empty -> False (the skill's default digest, trypsin). Raises Refused
    for a word no rule can be written for."""
    if not str(enzyme or "").strip():
        return False
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import cleavage_rules
    try:
        return cleavage_rules.enzyme_rule(enzyme)["specificity"] == "none"
    except cleavage_rules.EnzymeError as e:
        raise Refused(str(e)) from None


def routes_for(acq, fam, enzyme=None):
    """The engines for a data type, best first ([0] the default): ROUTES' row, or -- for a
    non-specific digest on a row NONSPECIFIC_ROUTES names -- that row (Ruling 120)."""
    key = (str(acq).upper(), fam)
    if key in NONSPECIFIC_ROUTES and non_specific(enzyme):
        return list(NONSPECIFIC_ROUTES[key])
    return ROUTES[key]


def data_type(acq, fam, enzyme=None):
    """What the data type is called in the client sentence (DATA_TYPE, or a non-specific row's)."""
    key = (str(acq).upper(), fam)
    if key in NONSPECIFIC_DATA_TYPE and non_specific(enzyme):
        return NONSPECIFIC_DATA_TYPE[key]
    return DATA_TYPE[key]


class Refused(Exception):
    """A REQUESTED engine cannot run here (or the digest names no rule); the message says why
    and what can."""


def de_method(engine, acq):
    """The DE input for an engine's output: "dpc" (limpa from precursors: DIA-NN, Radiant, and
    FragPipe's DIA route, whose quantities are DIA-NN's) or "maxlfq" (a protein table: Sage and
    FragPipe DDA, through run_search.py's adapter)."""
    if engine == "sage":
        return "maxlfq"
    if engine == "fragpipe":
        return "maxlfq" if str(acq).upper() == "DDA" else "dpc"
    return "dpc"



def _platform(env, engine):
    """The environment's verdict for `engine` (detect_env.sh), with FragPipe's Java floor applied
    to an env.json older than that check. None when no env was given."""
    if not env or not isinstance(env.get("engines"), dict):
        return None
    info = env["engines"].get(engine) or {}
    ok = info.get("available", True) if info else True
    note = info.get("note")
    if engine == "fragpipe" and ok and env.get("has_java"):
        j = java_record(env.get("has_java"))
        if j["ok"] is False:
            ok, note = False, j["why"]
    return {"available": bool(ok), "note": note, "os": env.get("os"), "arch": env.get("arch")}


def _sentence(how, engine, default, dtype):
    if how == "fallback":
        return (f"{NAMES[default]}, the pipeline's default search engine for {dtype}, was not "
                f"available where this search ran, so the data were searched with {NAMES[engine]}.")
    if how == "requested" and engine != default:
        return (f"{NAMES[engine]} was chosen for this analysis instead of {NAMES[default]}, the "
                f"pipeline's default search engine for {dtype}.")
    return None


def choose(acq, fam, requested=None, env=None, fp_tools=None, ids_only=False, enzyme=None,
           cls=None):
    """The engine for (acq, fam) -> (selection, staff_entry | None).

    requested  --engine (honoured, never overridden; refused only when it cannot run here)
    env        detect_env.sh's JSON (None: no platform check, as before 2.15)
    fp_tools   tools.json `fragpipe_tools` (or a live fragpipe_tools() probe); None = not inspected
               (does not block the default; run_search.py checks tools.json again). A record whose
               `licensed` is null (nobody could look) counts as NOT licensed (Ruling 103); one
               that is licensed but not `ready` (FragPipe cannot run) is a fallback too, never
               a choice run_search.py then refuses (Ruling 107.2)
    ids_only, cls   recorded; they never change the engine: an IDs-only search and an ion-trap
               MS2 stay on FragPipe (its Basic-Search workflow and its own ion-trap setting,
               Ruling 90)
    enzyme     recorded, and the digest picks the row (routes_for): a non-specific digest on
               Orbitrap DDA defaults to FragPipe (Ruling 120); on timsTOF DDA FragPipe was the
               default already (its Nonspecific workflows)
    Raises Refused for a requested engine that cannot run here, or an --enzyme word with no
    rule."""
    acq = str(acq).upper()
    routes = routes_for(acq, fam, enzyme)
    default, dtype = routes[0], data_type(acq, fam, enzyme)
    req = (requested or "").strip().lower() or None
    engine = req or default
    # asking for the data type's default is the default: a bundle's reproduce.sh replays step 4
    # with --engine <what ran>, and the replayed record must equal the original (provenance replay)
    how = "requested" if req and req != default else "default"
    reason = why = None
    checks = None

    def gate(e):
        p = _platform(env, e)
        if p and not p["available"]:
            alt = [x for x in routes if x != e and (_platform(env, x) or {}).get("available")]
            return p, (f"{e} cannot run on this machine -- {p['note'] or 'unavailable'}"
                       + (f"\n  Available here for this data type: {', '.join(alt)}" if alt
                          else "\n  No engine for this data type runs here; use HIVE."))
        return p, None

    plat, refusal = gate(engine)
    if engine == "fragpipe":
        ft = fp_tools if isinstance(fp_tools, dict) else None
        lic = ft.get("licensed") if ft else None
        java = (ft or {}).get("java") or {}
        checks = {"platform": plat,
                  "tools": None if ft is None else
                  {k: ft.get(k) for k in ("licensed", "ready", "folder", "source", "java", "why")}}
        if refusal is None and java.get("ok") is False:
            # the Java where tools.json was written (an env.json from before the Java check)
            plat = dict(plat or {"os": None, "arch": None}, available=False,
                        note=java_record(java.get("version"))["why"])
            refusal = f"fragpipe cannot run on this machine -- {plat['note']}"
        if refusal is not None:
            reason = "unavailable_here"
            why = (f"FragPipe cannot run on the machine that runs the search "
                   f"({plat.get('os')}/{plat.get('arch')}: {plat.get('note')})")
        elif lic is False:
            reason = "not_licensed_here"
            why = f"MSFragger/IonQuant are not installed for this account ({ft.get('why')})"
            refusal = f"fragpipe: {why}"
        elif ft is not None and lic is None:
            # tools.json says nobody could look (acquire_tools.sh without python3): not licensed
            # until shown otherwise (Ruling 103)
            reason = "licence_not_checked"
            why = f"licence not checked on this machine ({ft.get('why')})"
            refusal = f"fragpipe: {why}"
        elif ft is not None and ft.get("ready") is False:
            # licensed, but FragPipe cannot run (no runnable binary -- a failed download off
            # HIVE --, or a Java it could not read): run_search.py refuses a search whose
            # `ready` is false (refuse_fragpipe_not_ready), so the default falls back here
            # instead of choosing FragPipe and refusing at step 7 (review M2, Ruling 107.2)
            reason = "not_ready_here"
            why = f"the licensed tools are here but FragPipe cannot run ({ft.get('why')})"
            refusal = f"fragpipe: {why}"
        else:
            why = ("FragPipe can run here" if lic else
                   "FragPipe chosen; its tools were not inspected here (run_search.py checks "
                   "tools.json again before it writes anything)")
        if refusal and req:
            raise Refused(refusal)
        if refusal:
            # the default could not be used: the next route, recorded as a fallback
            engine, how = routes[1], "fallback"
            _, refusal = gate(engine)
    if refusal:
        raise Refused(refusal)

    staff_entry = None
    if checks is not None:
        staff_entry = {"event": "engine_selection", "engine": engine, "default": default,
                       "how": how, "reason": reason, "why": why, "checks": checks,
                       "ids_only": bool(ids_only), "enzyme": enzyme, "instrument_class": cls}
    selection = {"engine": engine, "default": default, "how": how, "data_type": dtype,
                 "sentence": _sentence(how, engine, default, dtype),
                 "ids_only": bool(ids_only), "staff_record": None,
                 "decided_by": "engine_choice.py"}
    return selection, staff_entry


# ---------------------------------------------------------------- CLI --------------------------
def _load(path):
    if not path:
        return None
    with open(path) as fh:
        return json.load(fh)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("fragpipe-tools", help="the tools.json fragpipe_tools object, as JSON")
    p.add_argument("--binary", default=None)
    p.add_argument("--pin", default=None)
    p = sub.add_parser("core-binary", help="print the Core's bin/fragpipe if usable (else nothing)")
    p.add_argument("--pin", default=None)
    p = sub.add_parser("java", help="the Java check, as JSON (or ok<TAB>why with --tsv)")
    p.add_argument("--line", default=None, help="a `java -version` first line; else java is run")
    p.add_argument("--tsv", action="store_true")
    p = sub.add_parser("choose", help="the selection record + staff entry, as JSON")
    p.add_argument("--acquisition", required=True)
    p.add_argument("--instrument", default="")
    p.add_argument("--engine", default=None)
    p.add_argument("--env", default=None)
    p.add_argument("--tools", default=None)
    p.add_argument("--ids-only", action="store_true")
    p.add_argument("--enzyme", default=None)
    a = ap.parse_args(argv)

    if a.cmd == "fragpipe-tools":
        print(json.dumps(fragpipe_tools(a.binary, os.environ.get(TOOLS_ENV) or None, a.pin),
                         indent=2))
    elif a.cmd == "core-binary":
        b, _ = core_fragpipe_binary(a.pin)
        if b:
            print(b)
    elif a.cmd == "java":
        j = java_record(a.line) if a.line is not None else java_status()
        if a.tsv:
            ok = {True: "true", False: "false"}.get(j["ok"], "unknown")
            print(f"{ok}\t{j['why'] or ''}")
        else:
            print(json.dumps(j, indent=2))
    elif a.cmd == "choose":
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from estimate_params import classify_instrument
        from resolve_defaults import family, ENGINE_VERSIONS
        cls = classify_instrument(a.instrument)[0]
        tools = _load(a.tools)
        ft = (tools or {}).get("fragpipe_tools") or fragpipe_tools(
            None, os.environ.get(TOOLS_ENV) or None, ENGINE_VERSIONS["fragpipe"])
        try:
            sel, staff = choose(a.acquisition, family(cls), a.engine, _load(a.env), ft,
                                a.ids_only, a.enzyme, cls)
        except Refused as e:
            sys.exit(f"engine_choice: {e}")
        print(json.dumps({"selection": sel, "staff": staff}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
