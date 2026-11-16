#!/usr/bin/env python3
"""
resolve_defaults.py -- pick the search defaults for a run FROM THE DATA TYPE.

Replaces the old fetch_workflows.py registry. Same output contract
(workflow.manifest.json + a JSON summary on stdout), so run_search.py and
provenance.py consume it unchanged -- but nothing is fetched over the network
and there is no menu.

WHY THE REGISTRY WENT AWAY
--------------------------
It fetched per-organism "validated workflow" bundles from GitHub at run time.
Three things were wrong with that:

  * Search parameters do not depend on species. Only the FASTA does, and that
    always came from the user's own organism answer -- never from the bundle.
    So the organism dimension bought nothing and forced a near-duplicate bundle
    per species.
  * A remote, mutable registry meant one skill version produced different
    parameters on different days. In practice two bundles pinned two different
    DIA-NN versions (2.3.0 and 2.6.0) and a run silently got whichever one its
    instrument matched.
  * Every bundle was unvalidated or placeholder-validated, so the "is this
    validated?" gate fired on essentially every run and turned the first
    interaction into a menu of warnings.

Defaults now ship WITH the skill, so the skill version fully determines the
search. Site-specific overrides are explicit flags, not a remote lookup.

  python3 resolve_defaults.py --acquisition DIA --instrument "timsTOF HT" \
      --engine diann --dest ./wf [--env env.json] [--fasta db.fasta] [--threads 32]
  Orbitrap: add --ms1-resolution 120000 --ms2-resolution 15000 [--resolution-source detected]
  (--ms1-res/--ms2-res are the same flags). Without them a DIA-NN or Radiant search of a
  non-Astral Orbitrap sets needs_confirmation + ask_user (exit 0). MS2 read in the ion trap
  (step 2's ms2_ion_trap): --ms1-resolution 120000 --ms2-analyzer ITMS, no --ms2-resolution.
"""
import argparse
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from estimate_params import (classify_instrument, MEASURE_CLASSES,  # noqa: E402  (one ppm table)
                             add_resolution_args, resolution_args_error, resolution_source,
                             resolution_record, resolution_question, SAGE_ITMS_FRAGMENT_DA,
                             override_setter, override_source, override_record,
                             staff_record_name, record_override_setter)
import staff  # noqa: E402  who may be called Core staff (role_for)

HERE = os.path.dirname(os.path.abspath(__file__))
MANIFEST = "workflow.manifest.json"

# Version the DEFAULTS TABLE itself, so a run record can name what produced its
# parameters without depending on a git commit in another repo.
DEFAULTS_VERSION = "2026-09-23"

# Pinned engine versions. One place. Bump here, not in six files.
# diann 2.6.1 -> 2.7.0 on 2026-09-23 (Brett). The Core's copy is
# /quobyte/proteomics-grp/dia-nn/build_270/diann-2.7.0/, unpacked from the Academia Linux zip
# whose sha256 matched the GitHub release asset (246ca8ba...cb223); off HIVE,
# acquire_tools.sh downloads that same zip.
ENGINE_VERSIONS = {
    "diann": "2.7.0",
    "sage": "0.14.7",
    "fragpipe": "24.0",
    "radiant": "2.3.3",   # container release, NOT the Radiant DIA version inside
    # The Docker Hub tag of the Core's image (METAMORPHEUS_PIN). The binary inside it does not
    # call itself 1.1.12: `dotnet CMD.dll --version` prints "CMD 1.0.0+2591d745..." and its
    # results.txt says "MetaMorpheus: version Not a release version" (HIVE, 2026-11-16).
    "metamorpheus": "1.1.12",
}

# The Core's MetaMorpheus container (run_search.py --engine metamorpheus). Built on HIVE
# 2026-11-04 from docker://smithchemwisc/metamorpheus:1.1.12 (MIT); its labels and sha256 were read
# on HIVE 2026-11-16 (job 25758472). The app is /MetaMorpheus/CMD.dll on the image's .NET 10.0.12.
# MetaMorpheus WRITES INTO ITS OWN FOLDER on every start (Data/CustomCrosslinkers.tsv,
# settings.toml where --acceptThermoLicence records the agreement), so run read-only from the
# image it aborts ("Read-only file system: /MetaMorpheus/Data/CustomCrosslinkers.tsv.tmp", job
# 25758490): the job copies /MetaMorpheus (not its duplicate publish/ folder) to node-local
# $TMPDIR and runs that copy, inside the same image. The sif's sha256 is checked by the job first.
METAMORPHEUS_PIN = {
    "version": "1.1.12",
    "image": "smithchemwisc/metamorpheus",
    "tag": "1.1.12",
    "docker_base_digest": "sha256:33ac82b46374d6b5616fcd218f3fd4f0c5f1a574e28a04533d51a5431e5ab92d",
    "sif": "/quobyte/proteomics-grp/apptainers/metamorpheus_1.1.12.sif",
    "sif_sha256": "221d0d78237860c2422ff099317c476c48cfe7086bd9cd868a2cd64730a8327c",
    "app_dir": "/MetaMorpheus",
    "app_skip": "publish",
    # what `dotnet CMD.dll --version` prints; results.txt says only "Not a release version".
    # The commit IS the GitHub tag 1.1.12 (git/ref/tags/1.1.12 -> 2591d745..., checked
    # 2026-11-16), so the version rests on the source, not only on the Docker tag.
    "build": "CMD 1.0.0+2591d74550a78618dbb7198b46591e45fc86e8cc",
    "source_commit": "2591d74550a78618dbb7198b46591e45fc86e8cc",
    "source_tag": "https://github.com/smith-chem-wisc/MetaMorpheus/releases/tag/1.1.12",
    "runtime": ".NET 10.0.12 (in the image)",
    "license": "MIT",
}

# MetaMorpheus is configured by a whole search-task TOML, shipped in presets/metamorpheus/ (not
# generated: the one there is the file a comparison validated). Each preset names the custom
# modification files its TOML needs: MetaMorpheus 1.1.12's own Mods/ has no peptide-N-terminal
# pyroGlu, so the validated search defined it (Unimod 28 / 27) in its own copy of the app; the job
# installs these into the run's copy of the app (run_search.py metamorpheus_job_lines).
METAMORPHEUS_PRESET_DIR = os.path.join(os.path.dirname(HERE), "presets", "metamorpheus")
METAMORPHEUS_PRESETS = {
    "nonspecific_peptidomics": {
        "file": "NonSpecific_Peptidomics_HCD.toml",
        "md5": "7f5d2d874b0cb577afe5c8d4f901b284",
        "mods": ["Mods/CustomPyroGluPeptideNterm.txt"],
        "description": ("non-specific digest (peptidomics), 7-25 residues, precursor +/-10 ppm, "
                        "product +/-20 ppm (HCD, Orbitrap MS2), variable Oxidation (M), "
                        "N-terminal acetylation (any residue), peptide-N-terminal pyroGlu from Q "
                        "and E; reverse decoys; PSM q-value 0.01; no LFQ, no MBR, no parsimony"),
        "validated": {
            "date": "2026-11-04",
            "by": "UC Davis Proteomics Core",
            "dataset": ("a Core submission's 4-engine comparison: Orbitrap Exploris 480 HCD "
                        "DDA, non-specific milk peptidomics, harmonised with the FragPipe and "
                        "Sage settings of the same comparison"),
            "notes": ("Run as MM_SearchTask_common.toml in that 2026-11-04 four-engine session "
                      "(md5 7f5d2d874b0cb577afe5c8d4f901b284, the same bytes as this file). "
                      "Validated for that data type only: Orbitrap HCD MS2 read in the Orbitrap, "
                      "non-specific digest."),
        },
    },
}
METAMORPHEUS_DEFAULT_PRESET = "nonspecific_peptidomics"

# Pinned Skyline, for the peak check after DE (skyline_peaks.py; references/skyline.md). Not a
# search engine: SkylineCmd runs headless under Wine in ProteoWizard's vendor-reader container.
# The release (not a daily): Skyline 26.1.0.266, skyline.ms's current download on 2026-10-20 and
# GitHub tag Skyline-26.1.0.266. Pulled by digest into the Core's container folder on HIVE; the
# .sif's sha256 is checked before every use (skyline_peaks.verify_container). Anyone else pulls
# the same digest themselves (references/skyline.md); a .sif built elsewhere has another sha256,
# so its version is checked from SkylineCmd --version instead.
SKYLINE_PIN = {
    "version": "26.1.0.266",
    "commit": "578b175",
    "image": "proteowizard/pwiz-skyline-i-agree-to-the-vendor-licenses",
    "tag": "skyline_26.1.0.266-578b175",
    "docker_digest": "sha256:448a833eb92eac2ca731f2461c046963c0c2afafb3550b4dfb6fc3914cec7b0e",
    "sif": "/quobyte/proteomics-grp/apptainers/"
           "pwiz-skyline-i-agree-to-the-vendor-licenses_skyline_26.1.0.266-578b175.sif",
    "sif_sha256": "b5f135587e9b589dd2d4686ff520573ea237b84f53b1343d5b76213c16f302bf",
    # ProteoWizard msconvert's own version inside this image (run_search.py converts a timsTOF .d
    # with it): read on HIVE 2026-11-01 (job 25387830) from the sif above -- `wine msconvert
    # --help` prints "ProteoWizard release: 3.0.26266 (578b175)" near its end (lines 613-614, so
    # `head` does not find it), and the mzML it writes records pwiz_Reader_Bruker 3.0.26266. Not
    # Skyline's 26.1.0.266, and not 3.0.26058 (the Windows build whose settings it runs with).
    "msconvert_version": "3.0.26266",
}

# Which engines are sensible for a data type, best first. Instrument class comes
# from estimate_params.classify_instrument, so the ppm table and this table can
# never disagree about what an instrument is.
ROUTES = {
    ("DIA", "timstof"): ["diann", "fragpipe"],
    ("DIA", "orbitrap"): ["diann", "fragpipe", "radiant"],
    ("DIA", "other"): ["diann"],
    ("DDA", "timstof"): ["sage"],
    ("DDA", "orbitrap"): ["sage"],
    ("DDA", "other"): ["sage"],
}

DE_METHOD = {"diann": "dpc", "fragpipe": "dpc", "radiant": "dpc", "sage": "maxlfq",
             # identifications only: no quantities, so no DE method (run_de.R is refused by the
             # search's not_applicable record, sage_lfq_check.not_applicable)
             "metamorpheus": None}

# What each engine does with an Orbitrap mass-accuracy level that has no documented DIA-NN value
# (ms1_ppm/ms2_ppm None: resolution unknown, or outside the README's 30k-240k table). Appended to
# ppm_source, which SKILL.md tells the agent to quote when it asks the user to confirm -- so it
# has to be true for THIS engine. It once said "measured in step 1b of the 5-step chain" for
# Radiant and FragPipe as well, which have no step 1b and run at their vendor tolerances.
UNDOCUMENTED_LEVEL = {
    "diann": "DIA-NN: the missing level is measured with DIA-NN on representative runs before "
             "the search (estimate_params.py plan measure_with_diann)",
    "radiant": "Radiant: nothing to derive the missing level from, so its vendor 20 ppm "
               "extraction width is kept for it",
    "fragpipe": "FragPipe: its vendor preset tolerances are used (the DIA-NN table is not "
                "applied to MSFragger)",
    "sage": "Sage: its own derived tolerances are used (estimate_params.py --engine sage)",
    "metamorpheus": "MetaMorpheus: the tolerances of its validated task TOML are used",
}

# The same, for an MS2 read in the ion trap (class orbitrap_iontrap): there is no Orbitrap MS2
# level to look up, and "measured with DIA-NN" would be false too -- an ion trap's fragment error
# is Da-scale, far outside the 3-30 ppm a measurement is accepted in.
ION_TRAP_MS2_ROUTE = {
    "diann": "DIA-NN: neither level is pinned (either flag fixes both, and there is no MS2 "
             "value), so mass accuracy is left to DIA-NN calibration",
    "radiant": "Radiant: nothing to derive an ion-trap MS2 width from, so its vendor 20 ppm MS2 "
               "extraction width (an Orbitrap-scale value) is kept",
    "fragpipe": "FragPipe: its vendor Thermo preset tolerances are used",
    "sage": f"Sage: fragment tolerance +/-{SAGE_ITMS_FRAGMENT_DA} Da, Sage's documented low-res "
            "MS/MS setting (estimate_params.py --engine sage --ms2-analyzer ITMS)",
    "metamorpheus": "MetaMorpheus: refused -- its validated task TOML is for Orbitrap MS2",
}

# Engines configured by a whole config file rather than by flags. These are
# generated per run by make_presets.py.
PRESET_ENGINES = {"fragpipe": ".workflow", "radiant": ".radiantConfig"}


def family(cls):
    if cls == "timstof":
        return "timstof"
    if cls.startswith("orbitrap"):
        return "orbitrap"
    return "other"


def metamorpheus_refusal(engine, acq, cls, label, args):
    """Why a MetaMorpheus search of this data type is not run with the shipped TOML, or None.
    MetaMorpheus is DDA here; the validated TOML reads an Orbitrap MS2 (+/-20 ppm product); its
    tolerances are the validated file's, so --ms1-ppm/--ms2-ppm are not applied to it."""
    if acq != "DDA":
        return ("MetaMorpheus is wired for DDA only (identifications); this is "
                f"{acq}. Search DIA with diann (or fragpipe / radiant).")
    if cls == "orbitrap_iontrap":
        return ("the validated MetaMorpheus task TOML matches fragments at +/-20 ppm, an "
                "Orbitrap MS2; this data's MS2 was read in the ion trap. Search it with sage "
                "(estimate_params.py --engine sage --ms2-analyzer ITMS), or give run_search.py a "
                "task TOML of your own with ProductMassTolerance in Da (--params).")
    if args.ms1_ppm is not None or args.ms2_ppm is not None:
        return ("--ms1-ppm/--ms2-ppm are not applied to a MetaMorpheus task TOML: its tolerances "
                "are the validated file's. Copy "
                f"{os.path.join(METAMORPHEUS_PRESET_DIR, METAMORPHEUS_PRESETS[args.metamorpheus_preset]['file'])}"
                ", edit PrecursorMassTolerance / ProductMassTolerance, and give run_search.py "
                "that copy as --params (it is then recorded as not the validated file).")
    return None


def install_metamorpheus_preset(name, dest):
    """Copy the shipped, validated task TOML into `dest` as metamorpheus.toml -- refused if its
    bytes are not the validated md5 -- and return (its path, the manifest record)."""
    import shutil
    import metamorpheus_params as mmp
    pre = METAMORPHEUS_PRESETS[name]
    src = os.path.join(METAMORPHEUS_PRESET_DIR, pre["file"])
    got = mmp.md5(src)
    if got != pre["md5"]:
        sys.exit(f"resolve_defaults: REFUSED -- {src} is not the validated MetaMorpheus task "
                 f"TOML (md5 {got}, validated {pre['md5']}). The skill's copy was changed; "
                 f"restore it from the skill's git history.")
    path = os.path.join(dest, "metamorpheus.toml")
    shutil.copyfile(src, path)
    task = mmp.read_task(path)
    st = mmp.summary(task)
    return path, {
        "preset": name, "shipped_file": src, "md5": got, "description": pre["description"],
        "settings": st,
        "ms1_ppm": mmp.tolerance_ppm(st["precursor_tolerance"]),
        "ms2_ppm": mmp.tolerance_ppm(st["product_tolerance"]),
        "custom_mods": [os.path.join(METAMORPHEUS_PRESET_DIR, m) for m in pre["mods"]],
        "validated": pre["validated"],
        "container": METAMORPHEUS_PIN,
    }


def load_env(path):
    if not path:
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError) as e:
        sys.stderr.write(f"[resolve_defaults] could not read --env {path}: {e}\n")
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--acquisition", required=True, choices=["DIA", "DDA", "dia", "dda"])
    ap.add_argument("--instrument", default="")
    ap.add_argument("--engine", default="", help="force an engine; else the data-type default")
    ap.add_argument("--organism-taxid", default=None,
                    help="recorded for provenance only -- it does NOT affect any search parameter")
    add_resolution_args(ap)
    ap.add_argument("--ms1-ppm", type=float, default=None,
                    help="the user's MS1 tolerance (ppm), recorded as a user override with who "
                         "set it and why -- never as an SOP")
    ap.add_argument("--ms2-ppm", type=float, default=None,
                    help="the user's MS2 tolerance (ppm); recorded the same way")
    ap.add_argument("--override-by", default=None,
                    help="who chose --ms1-ppm/--ms2-ppm (default: the login running this, "
                         "recorded as such) -- in the staff-only workflow.manifest.json.staff.json "
                         "only; the manifest says Core staff when that login is on the Core staff "
                         "list (staff.py), else the user")
    ap.add_argument("--override-reason", default=None,
                    help="why, in the user's words; without it the record says it was not given")
    ap.add_argument("--fasta", default=None)
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--env", default=None, help="detect_env.sh JSON, to flag engines that "
                                                "cannot run on this machine")
    ap.add_argument("--dest", default="./wf")
    ap.add_argument("--metamorpheus-preset", default=METAMORPHEUS_DEFAULT_PRESET,
                    choices=sorted(METAMORPHEUS_PRESETS),
                    help="--engine metamorpheus: which shipped, validated task TOML to search with "
                         f"(default {METAMORPHEUS_DEFAULT_PRESET}: "
                         f"{METAMORPHEUS_PRESETS[METAMORPHEUS_DEFAULT_PRESET]['description']})")
    args = ap.parse_args()

    acq = args.acquisition.upper()
    r1, r2, analyzer = args.ms1_resolution, args.ms2_resolution, args.ms2_analyzer
    bad = resolution_args_error(r2, analyzer)
    if bad:
        sys.exit(f"resolve_defaults: {bad}")
    res_src = resolution_source(args.resolution_source, r1, r2, analyzer)
    cls, ms1, ms2, label, src = classify_instrument(args.instrument, r1, r2, res_src, analyzer)
    fam = family(cls)
    routes = ROUTES[(acq, fam)]

    engine = (args.engine or "").strip().lower()
    if engine and engine not in ENGINE_VERSIONS:
        sys.exit(f"resolve_defaults: unknown engine '{engine}'. "
                 f"Known: {', '.join(sorted(ENGINE_VERSIONS))}")
    notes = []
    if engine and engine not in routes:
        # Honour it, but never silently: the user asked for something this data
        # type does not normally use.
        notes.append(f"{engine} is not a default route for {acq} on {label} "
                     f"(defaults: {', '.join(routes)}) -- running it because it was requested")
    if not engine:
        engine = routes[0]

    # Platform gate. Verified facts live in detect_env.sh; we only read them.
    env = load_env(args.env)
    availability = None
    if env and isinstance(env.get("engines"), dict):
        info = env["engines"].get(engine) or {}
        availability = info
        if info and not info.get("available", True):
            alt = [e for e in routes
                   if (env["engines"].get(e) or {}).get("available")]
            sys.exit(f"resolve_defaults: {engine} cannot run on this machine -- "
                     f"{info.get('note', 'unavailable')}"
                     + (f"\n  Available here for this data type: {', '.join(alt)}"
                        if alt else "\n  No engine for this data type runs here; use HIVE."))
        if info.get("note"):
            notes.append(f"{engine}: {info['note']}")

    # MetaMorpheus: refused here, before anything is written, for what its validated TOML is not
    mm = metamorpheus_refusal(engine, acq, cls, label, args) if engine == "metamorpheus" else None
    if mm:
        sys.exit(f"resolve_defaults: REFUSED -- nothing was written. {mm}")

    dest = os.path.abspath(args.dest)
    os.makedirs(dest, exist_ok=True)

    # Engines driven by a config file get one generated now, before the search.
    params_file = None
    preset_prov = None
    if engine in PRESET_ENGINES:
        params_file = os.path.join(dest, f"{engine}{PRESET_ENGINES[engine]}")
        cmd = [sys.executable, os.path.join(HERE, "make_presets.py"),
               "--engine", engine, "--instrument", args.instrument,
               "--acquisition", acq, "--out", params_file]
        for flag, val in (("--ms1-res", r1), ("--ms2-res", r2), ("--ms2-analyzer", analyzer),
                          ("--resolution-source", res_src),
                          ("--ms1-ppm", args.ms1_ppm), ("--ms2-ppm", args.ms2_ppm),
                          ("--fasta", args.fasta), ("--threads", args.threads)):
            if val is not None:
                cmd += [flag, str(val)]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"resolve_defaults: preset generation failed:\n{r.stderr.strip()}")
        preset_prov = json.loads(r.stdout)

    mm_preset = None
    if engine == "metamorpheus":
        params_file, mm_preset = install_metamorpheus_preset(args.metamorpheus_preset, dest)
        ms1 = mm_preset["ms1_ppm"]
        ms2 = mm_preset["ms2_ppm"]
        src = (f"the validated MetaMorpheus task TOML ({mm_preset['preset']}): precursor "
               f"{mm_preset['settings']['precursor_tolerance']}, product "
               f"{mm_preset['settings']['product_tolerance']}")
        notes.append(f"metamorpheus: identifications only (no quantification, no DE); the task "
                     f"TOML {mm_preset['preset']} was validated on "
                     f"{mm_preset['validated']['dataset']}")

    # A user's --ms1-ppm/--ms2-ppm REPLACES that level in the manifest below, so ppm_source has
    # to account for BOTH flags. Gated on --ms1-ppm alone, an override that sets only --ms2-ppm
    # left ppm_source saying MS2 "is measured with DIA-NN before the search" beside the user's
    # own number -- and SKILL.md tells the agent to read this line out to the user for
    # confirmation. Each given level is a USER override (estimate_params.override_source: the
    # value, who set it and why); it was called a "site SOP override", which the skill cannot
    # know it is (staff report, 2026-09-25: a calibration workaround recorded as an SOP).
    given = {f"--{lvl.lower()}-ppm": v for lvl, v in (("MS1", args.ms1_ppm), ("MS2", args.ms2_ppm))
             if v is not None}
    if (args.override_by or args.override_reason) and not given:
        sys.exit("resolve_defaults: --override-by/--override-reason describe --ms1-ppm/--ms2-ppm, "
                 "and neither was given")
    by, by_how = override_setter(args.override_by) if given else (None, "unknown")
    # "Core staff" only for a login on the Core staff list, else "the user" (S4)
    role = staff.role_for(by) if given else staff.USER_ROLE
    sop = [f"{f[2:5].upper()}: {override_source(f, f'{v:g}', args.override_reason, role)}"
           for f, v in given.items()]
    # What the manifest will actually carry: the user's number where they gave one, the table's
    # otherwise. UNDOCUMENTED_LEVEL applies to a level that is STILL without a value -- an
    # override that supplies the untabled level means nothing is left to measure, so saying it
    # "is measured with DIA-NN before the search" beside the user's own number is simply false.
    eff_ms1 = args.ms1_ppm if args.ms1_ppm is not None else ms1
    eff_ms2 = args.ms2_ppm if args.ms2_ppm is not None else ms2
    derived = (f"{src}; {ION_TRAP_MS2_ROUTE[engine]}" if cls == "orbitrap_iontrap" else
               f"{src}; {UNDOCUMENTED_LEVEL[engine]}"
               if cls in MEASURE_CLASSES and (eff_ms1 is None or eff_ms2 is None) else src)
    ppm_source = ("; ".join(sop) if len(sop) == 2 else derived if not sop else
                  f"{sop[0]}; the instrument's own values: {derived}")

    # gabrig 2026-09-23: a Fusion Lumos with no resolution exited 0 with "resolution unknown" in
    # the manifest and nothing asked the user. Still exit 0 -- callers rely on it -- but say so
    # on stderr and at the top of the manifest and stdout, where step 4's confirmation reads.
    ask = resolution_question(cls, engine, r1, r2)
    if ask:
        sys.stderr.write(f"[resolve_defaults] NEEDS CONFIRMATION: {ask}\n")

    # Same manifest shape the old registry emitted, so downstream steps are
    # unchanged. Fields that only made sense for a remote bundle are explicit
    # about being local now rather than being filled with plausible-looking junk.
    manifest = {
        "id": f"{engine}_{fam}_{acq.lower()}",
        "name": f"{engine} — {label} ({acq})",
        "path": None,
        "acquisition": acq,
        "instruments": [args.instrument] if args.instrument else [],
        "instrument_class": cls,
        "instrument_label": label,
        "resolution": resolution_record(r1, r2, res_src, analyzer),
        "needs_confirmation": bool(ask),
        "ask_user": ask,
        "organism": None,
        "organism_taxid": int(args.organism_taxid) if args.organism_taxid else None,
        "engine": {"name": engine, "version": ENGINE_VERSIONS[engine]},
        "fasta": {"path": os.path.abspath(args.fasta) if args.fasta else None,
                  "add_contaminants": True},
        "de": {"method": DE_METHOD[engine], "q_cutoff": 0.01, "logfc": 1.0, "adjp": 0.05},
        "search": {
            "estimate_params": engine not in PRESET_ENGINES and engine != "metamorpheus",
            "params_file": params_file,
            "ms1_ppm": eff_ms1,
            "ms2_ppm": eff_ms2,
            "ppm_source": ppm_source,
            # each --ms1-ppm/--ms2-ppm the user gave: by role, and why ({} when none); who, by
            # name, is in the staff-only record beside the manifest
            "overrides": override_record(given, args.override_reason,
                                         staff_record_name(os.path.join(dest, MANIFEST))
                                         if given else None, role=role),
            "preset_provenance": preset_prov,
            # --engine metamorpheus: the shipped task TOML copied, its md5, settings and the
            # custom modification files the job installs (install_metamorpheus_preset)
            "metamorpheus_preset": mm_preset,
        },
        "validated": (dict(mm_preset["validated"]) if mm_preset else {
            "date": None, "by": None, "dataset": None,
            "notes": ("Parameters are this skill's data-type defaults "
                      f"(resolve_defaults.py {DEFAULTS_VERSION}), not a site validation run. "
                      "They are reproducible from the skill version alone."),
        }),
        "registry": {"source": "local skill defaults",
                     "defaults_version": DEFAULTS_VERSION,
                     "note": "no remote registry; parameters ship with the skill"},
        "platform": availability,
        "alternatives": [e for e in routes if e != engine],
        "notes": notes,
    }
    mpath = os.path.join(dest, MANIFEST)
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=2)
    if given:
        record_override_setter(mpath, given, by, by_how, args.override_reason, role)

    print(json.dumps({
        "id": manifest["id"], "dir": dest, "manifest": mpath,
        "engine": manifest["engine"], "de": manifest["de"],
        "params_file": params_file,
        "estimate_params": manifest["search"]["estimate_params"],
        "ms1_ppm": manifest["search"]["ms1_ppm"],
        "ms2_ppm": manifest["search"]["ms2_ppm"],
        "ppm_source": manifest["search"]["ppm_source"],
        "instrument_label": label,
        "resolution": manifest["resolution"],
        "needs_confirmation": manifest["needs_confirmation"],
        "ask_user": ask,
        "alternatives": manifest["alternatives"],
        "notes": notes,
    }, indent=2))


if __name__ == "__main__":
    main()
