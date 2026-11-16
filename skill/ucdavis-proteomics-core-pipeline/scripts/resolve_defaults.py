#!/usr/bin/env python3
"""
resolve_defaults.py -- pick the search defaults for a run FROM THE DATA TYPE.

Replaces the old fetch_workflows.py registry (removed in 2.15). Same output contract
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
      [--tools tools.json] [--ids-only] [--enzyme none]
  The engine is engine_choice.choose()'s (the one route table): timsTOF DDA-PASEF, and Orbitrap
  DDA with --enzyme none (or pepsin; Ruling 120), -> FragPipe where it runs and is licensed, else
  Sage, recorded in engine.selection (client-safe) and the staff-only
  workflow.manifest.json.staff.json (the reason).
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
                             staff_record_name, record_override_setter,
                             unit_resolution_amounts_refused)
import staff  # noqa: E402  who may be called Core staff (role_for)
# the ONE route table, the DE method by engine + acquisition, and the engine choice with its
# FragPipe checks (2.15): resolve_defaults only reads them
from engine_choice import (ROUTES, routes_for, de_method, choose, Refused, fragpipe_tools,  # noqa: E402
                           java_record, TOOLS_ENV)
# FragPipe's ion-trap MS2 setting is defined once, in make_presets (which writes the values and
# records them); make_presets does not import this module, so the import is not circular
from make_presets import ION_TRAP_SOURCE  # noqa: E402

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
}

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
    # with it): read on HIVE 2026-11-01 (a test job) from the sif above -- `wine msconvert
    # --help` prints "ProteoWizard release: 3.0.26266 (578b175)" near its end (lines 613-614, so
    # `head` does not find it), and the mzML it writes records pwiz_Reader_Bruker 3.0.26266. Not
    # Skyline's 26.1.0.266, and not 3.0.26058 (the Windows build whose settings it runs with).
    "msconvert_version": "3.0.26266",
}

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
}

# The same, for an MS2 read in the ion trap (class orbitrap_iontrap): there is no Orbitrap MS2
# level to look up, and "measured with DIA-NN" would be false too -- an ion trap's fragment error
# is Da-scale, far outside the 3-30 ppm a measurement is accepted in.
ION_TRAP_MS2_ROUTE = {
    "diann": "DIA-NN: neither level is pinned (either flag fixes both, and there is no MS2 "
             "value), so mass accuracy is left to DIA-NN calibration",
    "radiant": "Radiant: nothing to derive an ion-trap MS2 width from, so its vendor 20 ppm MS2 "
               "extraction width (an Orbitrap-scale value) is kept",
    # Ruling 90: FragPipe's own ion-trap MS2 setting (its TMT10-MS3/TMT16-MS3 workflows), not the
    # 20 ppm of its label-free presets; make_presets.py writes the values and records them.
    "fragpipe": f"FragPipe: {ION_TRAP_SOURCE}",
    "sage": f"Sage: fragment tolerance +/-{SAGE_ITMS_FRAGMENT_DA} Da, Sage's documented low-res "
            "MS/MS setting (estimate_params.py --engine sage --ms2-analyzer ITMS)",
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


def load_env(path):
    if not path:
        return None
    try:
        with open(path) as fh:
            return json.load(fh)
    except (OSError, ValueError) as e:
        sys.stderr.write(f"[resolve_defaults] could not read --env {path}: {e}\n")
        return None


def load_fragpipe_tools(tools_path, env, probe=True):
    """tools.json's fragpipe_tools, or -- with no tools.json, or one from before 2.15 -- the same
    check made here (engine_choice.fragpipe_tools), with the env's Java when it has one; None
    without a tools.json record when `probe` is False."""
    if tools_path:
        try:
            with open(tools_path) as fh:
                ft = json.load(fh).get("fragpipe_tools")
            if isinstance(ft, dict):
                return ft
        except (OSError, ValueError, AttributeError) as e:
            sys.stderr.write(f"[resolve_defaults] could not read --tools {tools_path}: {e}\n")
    if not probe:
        return None
    java = java_record(env.get("has_java")) if env and "has_java" in env else None
    return fragpipe_tools(None, os.environ.get(TOOLS_ENV) or None,
                          ENGINE_VERSIONS["fragpipe"], java=java)


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
    ap.add_argument("--tools", default=None,
                    help="acquire_tools.sh's tools.json: its fragpipe_tools says whether FragPipe's "
                         "licence-gated tools are here. Without it (or with a tools.json from "
                         "before 2.15) the same check is made here (engine_choice.fragpipe_tools)")
    ap.add_argument("--replay-of", default=None,
                    help="a bundle's copy of the ORIGINAL workflow.manifest.json (reproduce.sh "
                         "step 2): its engine.selection is kept, with `replayed` added")
    ap.add_argument("--ids-only", action="store_true",
                    help="identifications only, no differential expression: the engine is the "
                         "same (FragPipe's Basic-Search workflow, or Sage without LFQ)")
    ap.add_argument("--enzyme", default=None,
                    help="the digest, in fetch_fasta's vocabulary (none = non-specific: "
                         "FragPipe's Nonspecific workflow); passed to the FragPipe DDA preset")
    ap.add_argument("--semi-specific", action="store_true",
                    help="search the digest semi-specifically (one end of each peptide follows the "
                         "--enzyme rule, else trypsin); passed to the Sage / FragPipe DDA preset. "
                         "See make_presets.py --help")
    ap.add_argument("--dest", default="./wf")
    args = ap.parse_args()

    acq = args.acquisition.upper()
    r1, r2, analyzer = args.ms1_resolution, args.ms2_resolution, args.ms2_analyzer
    bad = resolution_args_error(r2, analyzer)
    if bad:
        sys.exit(f"resolve_defaults: {bad}")
    res_src = resolution_source(args.resolution_source, r1, r2, analyzer)
    cls, ms1, ms2, label, src = classify_instrument(args.instrument, r1, r2, res_src, analyzer)
    fam = family(cls)
    # the data type's row, the digest included: a non-specific Orbitrap DDA search defaults to
    # FragPipe (engine_choice.routes_for, Ruling 120)
    try:
        routes = routes_for(acq, fam, args.enzyme)
    except Refused as e:
        sys.exit(f"resolve_defaults: --enzyme {args.enzyme}: {e}")

    engine = (args.engine or "").strip().lower()
    if engine and engine not in ENGINE_VERSIONS:
        sys.exit(f"resolve_defaults: unknown engine '{engine}'. "
                 f"Known: {', '.join(sorted(ENGINE_VERSIONS))}")
    # A unit-resolution DDA search reports no amounts on ANY engine: a fact about the data, so it
    # is refused before the engine is chosen, nothing written (Ruling 107.3; the one rule is
    # estimate_params.unit_resolution_amounts_refused)
    no_amounts = unit_resolution_amounts_refused(args.instrument, acq, args.ids_only)
    if no_amounts:
        sys.exit(f"resolve_defaults: {no_amounts}. Nothing was written.")
    notes = []
    if engine and engine not in routes:
        # Honour it, but never silently: the user asked for something this data
        # type does not normally use.
        notes.append(f"{engine} is not a default route for {acq} on {label} "
                     f"(defaults: {', '.join(routes)}) -- running it because it was requested")

    # The choice (engine_choice.choose): a requested engine is honoured, and refused only when it
    # cannot run here; the default timsTOF DDA route falls back from FragPipe to Sage, recorded.
    # Verified platform facts live in detect_env.sh; we only read them.
    env = load_env(args.env)
    fp_tools = None
    if (engine or routes[0]) == "fragpipe":
        # a REQUESTED FragPipe is refused only on tools.json's word (a workflow may be prepared
        # here for a search elsewhere); the DEFAULT route also looks here when there is none
        fp_tools = load_fragpipe_tools(args.tools, env, probe=not engine)
    try:
        selection, staff_entry = choose(acq, fam, engine or None, env, fp_tools,
                                        ids_only=args.ids_only, enzyme=args.enzyme, cls=cls)
    except Refused as e:
        sys.exit(f"resolve_defaults: {e}")
    engine = selection["engine"]
    availability = platform_note = None
    if env and isinstance(env.get("engines"), dict):
        info = env["engines"].get(engine) or {}
        # what the machine offers goes in the manifest; detect_env.sh's note about it (its
        # wording names licences, jars, the tools-folder variable) goes to the staff-only record
        # beside it: the manifest is in the bundle a client receives (review M1, Ruling 107.1)
        availability = {k: v for k, v in info.items() if k != "note"}
        platform_note = info.get("note") or None
    if args.ids_only and (acq != "DDA" or engine not in ("fragpipe", "sage")):
        sys.exit(f"resolve_defaults: --ids-only applies to DDA searches only (Sage or FragPipe); "
                 f"a {acq} search with {engine} always quantifies. For a DIA analysis with no "
                 f"differential expression, leave --ids-only out: the search's own report is the "
                 f"result (references/no-de-route.md)")
    # A bundle's reproduce.sh replays step 4 with --engine <what ran> and --replay-of its copy of
    # the original manifest: the original choice is kept (a Sage FALLBACK is never relabelled
    # "requested"), and the replay says so (Ruling 103)
    replay_of = load_env(args.replay_of) if args.replay_of else None
    orig_sel = ((replay_of or {}).get("engine") or {}).get("selection")
    if isinstance(orig_sel, dict) and orig_sel.get("engine") == engine:
        # its staff record stayed with the original session: this one names one only when this
        # replay writes its own (below), never a file that is not here (review L6)
        selection = dict(orig_sel, replayed=f"replayed with --engine {engine}",
                         staff_record=None)

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
        # the DDA digest flags (make_presets: Basic-Search for --ids-only, Nonspecific for
        # --enzyme none)
        if acq == "DDA":
            if args.ids_only:
                cmd.append("--ids-only")
            if args.enzyme:
                cmd += ["--enzyme", args.enzyme]
            if args.semi_specific:
                cmd.append("--semi-specific")
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit(f"resolve_defaults: preset generation failed:\n{r.stderr.strip()}")
        preset_prov = json.loads(r.stdout)

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
    # FragPipe searches what its WORKFLOW says, not DIA-NN's instrument table: the manifest records
    # the tolerances written there (read back, the one reader), in ppm where the workflow's unit is
    # ppm, else None with the unit in ppm_source -- never the README's value beside a search that
    # ran another (X8 F5; the preset's own record of why is tolerances.source)
    if engine == "fragpipe" and params_file and preset_prov:
        import fragpipe_params as fpp
        wf_props, _ = fpp.read_workflow(params_file)
        kv = {k[len("msfragger."):]: v for k, v in wf_props.items() if k.startswith("msfragger.")}
        said = []
        for which, label_ in (("precursor", "MS1"), ("fragment", "MS2")):
            t = fpp.tolerance(kv, which)
            val = t["hi"] if t and t["symmetric"] and t["unit"] == "ppm" else None
            if which == "precursor":
                eff_ms1 = val
            else:
                eff_ms2 = val
            said.append(f"{label_} " + (f"+/-{t['hi']:g} {t['unit']}" if t and t["symmetric"] else
                                        "as the workflow states" if t else "not read"))
        tol_rec = preset_prov.get("tolerances")
        tol_why = tol_rec.get("source") if isinstance(tol_rec, dict) else tol_rec
        ppm_source = (f"as written in the FragPipe workflow ({', '.join(said)}); "
                      f"{tol_why or 'FragPipe preset'}")

    # A Core member, 2026-09-23: a Fusion Lumos with no resolution exited 0 with "resolution unknown" in
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
        "engine": {"name": engine, "version": ENGINE_VERSIONS[engine], "selection": selection},
        "fasta": {"path": os.path.abspath(args.fasta) if args.fasta else None,
                  "add_contaminants": True},
        # an identifications-only search has no differential expression: no DE block to carry
        # into the delivered run_manifest.json (X8 F5)
        "de": (None if args.ids_only else
               {"method": de_method(engine, acq), "q_cutoff": 0.01, "logfc": 1.0, "adjp": 0.05}),
        "search": {
            "estimate_params": engine not in PRESET_ENGINES,
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
        },
        "validated": {
            "date": None, "by": None, "dataset": None,
            "notes": ("Parameters are this skill's data-type defaults "
                      f"(resolve_defaults.py {DEFAULTS_VERSION}), not a site validation run. "
                      "They are reproducible from the skill version alone."),
        },
        "registry": {"source": "local skill defaults",
                     "defaults_version": DEFAULTS_VERSION,
                     "note": "no remote registry; parameters ship with the skill"},
        "platform": (dict(availability, note_in=staff_record_name(
            os.path.join(dest, MANIFEST))) if platform_note else availability),
        "alternatives": [e for e in routes if e != engine],
        "notes": notes,
    }
    mpath = os.path.join(dest, MANIFEST)
    # detect_env.sh's note on the engine that runs: in the selection's entry when there is one,
    # else an entry of its own (one entry per run of this step)
    note = {"platform_note": platform_note} if platform_note else {}
    if staff_entry:
        # the reason, `why` and the checks name paths, versions and the licence: staff-only
        staff.append(mpath + staff.STAFF_SUFFIX,
                     dict(staff_entry, **note, account=staff.login(), at=staff.now()))
        selection["staff_record"] = staff_record_name(mpath)
    elif platform_note:
        staff.append(mpath + staff.STAFF_SUFFIX,
                     {"event": "platform_note", "engine": engine, **note,
                      "source": "detect_env.sh (--env)", "account": staff.login(),
                      "at": staff.now()})
    with open(mpath, "w") as fh:
        json.dump(manifest, fh, indent=2)
    if given:
        record_override_setter(mpath, given, by, by_how, args.override_reason, role)

    print(json.dumps({
        "id": manifest["id"], "dir": dest, "manifest": mpath,
        "engine": manifest["engine"], "de": manifest["de"], "selection": selection,
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
