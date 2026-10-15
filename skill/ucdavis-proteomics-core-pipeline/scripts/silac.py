#!/usr/bin/env python3
"""
silac.py -- SILAC (Lys/Arg) multiplexing with DIA-NN plexDIA: THE one definition in this skill.

Why it exists (PROT_0807-type triple SILAC, 2026-10-14): the skill had no SILAC option, so the
labels went in as estimate_params.py --overrides and were recorded as a user override, the
declaration that DIA-NN's README gives for a label-free library (`--lib-fixed-mod SILAC`) was
applied on top of a FASTA digest that had already labelled the library -- the precursor ids read
K(SILAC)(SILAC), and K(SILAC)(SILAC)(SILAC) once the search step repeated it -- the search stopped
on "channel normalisation strategy must be specified", the ratio analysis was written by hand and
the client report was made outside the skill. Everything SILAC now comes from here:

  LABELS, PRESETS, SCHEMES  the label masses (UniMod) and the presets staff choose from
  labels_record()           estimate_params.py --labels: the params sidecar's `labels` record
  diann_cfg_lines()         the cfg lines build_diann() writes for it (never --lib-fixed-mod)
  cfg_problems()            run_search.py / diann_parallel.py: a cfg that would mis-search, refused
  labels_from_cfg()         the channels as any cfg declares them (a search made before 2.13 too)
  read_design()             the channel design (silac_design.csv) run_silac.R analyses
  methods_search_sentence() the Methods' plexDIA sentence (make_methods.py)

THE DECLARATION, verified on DIA-NN 2.7.0 (HIVE srun 24808559, 2026-10-15, a 6-protein FASTA): a
library predicted from a FASTA with `--fixed-mod SILAC,0.0,KR,label` alone reads K(SILAC) -- the
fixed mod is applied during the digest (README, --fixed-mod: "applying fixed modifications during
FASTA digest") -- and adding `--lib-fixed-mod SILAC` makes it K(SILAC)(SILAC). Every DIA-NN route
the skill runs predicts its library from the FASTA (two jobs, one step, the 5-step chain), and every
later step reads a library that already carries the label (README plexDIA Scenario 4: "--lib-fixed-mod
must not be supplied"), so --lib-fixed-mod is never written. --channels, --original-mods and a
channel-normalisation flag ("In all scenarios above, an extra option specifying the normalisation
strategy must be included") go to every step through the cfg.

Usage:
  python3 silac.py presets
  python3 silac.py flags --labels silac3 --label-preset K4R6,K8R10 [--channel-norm spec|run]
  python3 silac.py check-cfg <params.cfg>
  python3 silac.py design --validate silac_design.csv [--against report.parquet]
      [--channels L,M,H] [--comparisons "B-A,C-A"] [--json]
"""
import argparse
import csv
import itertools
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# --- the labels ------------------------------------------------------------------------------
# Monoisotopic mass shifts, each read from its UniMod entry on 2026-10-15
# (https://www.unimod.org/modifications_view.php?editid1=<id>):
#   481 Label:2H(4)          4.025107  "D4 Lysine"
#   188 Label:13C(6)         6.020129  "C(6) Silac label" (K and R)
#   259 Label:13C(6)15N(2)   8.014199  "C(6) 15N(2) Silac label", heavy lysine
#   267 Label:13C(6)15N(4)  10.008269  "C(6) 15N(4) Silac label", heavy arginine
# tests/test_silac.py holds these numbers too, so an edit here that drifts from UniMod fails.
LABELS = {
    "Lys0": {"residue": "K", "delta": 0.0, "unimod": None, "composition": None},
    "Arg0": {"residue": "R", "delta": 0.0, "unimod": None, "composition": None},
    "Lys4": {"residue": "K", "delta": 4.025107, "unimod": 481, "composition": "2H(4)"},
    "Lys6": {"residue": "K", "delta": 6.020129, "unimod": 188, "composition": "13C(6)"},
    "Lys8": {"residue": "K", "delta": 8.014199, "unimod": 259, "composition": "13C(6)15N(2)"},
    "Arg6": {"residue": "R", "delta": 6.020129, "unimod": 188, "composition": "13C(6)"},
    "Arg10": {"residue": "R", "delta": 10.008269, "unimod": 267, "composition": "13C(6)15N(4)"},
}
UNIMOD_SOURCE = "UniMod entries 481, 188, 259 and 267, read 2026-10-15"
# A channel's (lysine, arginine) labels. K0R0 is the light channel of every scheme.
PRESETS = {"K0R0": ("Lys0", "Arg0"), "K4R6": ("Lys4", "Arg6"), "K6R6": ("Lys6", "Arg6"),
           "K8R10": ("Lys8", "Arg10")}
LIGHT = "K0R0"
# --labels: the channel names, lightest first (DIA-NN wants each site's masses non-decreasing in
# the order the channels are listed, README --channels), and how many presets --label-preset takes.
SCHEMES = {"silac2": ("L", "H"), "silac3": ("L", "M", "H")}
LABEL_GROUP = "SILAC"
SITES = "KR"
FIXED_MOD = f"{LABEL_GROUP},0.0,{SITES},label"

# The channel-normalisation question (README: "--channel-run-norm (pulsed SILAC, protein
# turnover) or --channel-spec-norm (multiplexing of independent samples)"). Without one DIA-NN
# 2.7.0 logs an ERROR and falls back to --channel-run-norm (PROT_0807-type search, 2026-10-14).
CHANNEL_NORM = {"spec": "--channel-spec-norm", "run": "--channel-run-norm"}
CHANNEL_NORM_DESIGN = {
    "spec": "independent samples in each channel (different strains, treatments or cell lines, "
            "each grown in its own label, then mixed)",
    "run": "one sample per run whose channels are old and new protein (pulsed SILAC, protein "
           "turnover)",
}
CHANNEL_NORM_DEFAULT = "spec"
# What the ratio analysis (run_silac.R) supports. A turnover experiment asks a different question
# (the fraction of new protein over time), not "does this protein differ between the channels".
RATIO_ANALYSIS_NORMS = ("spec",)

LIB_FIXED_MOD_NOTE = (
    "not used: every route predicts its library from the FASTA, where --fixed-mod already "
    "labels K and R (DIA-NN 2.7.0, HIVE 2026-10-15: --fixed-mod alone gives K(SILAC), with "
    "--lib-fixed-mod K(SILAC)(SILAC)); later steps read a library that carries the label, and "
    "DIA-NN's README says --lib-fixed-mod must not be supplied then (plexDIA Scenario 4)")

# Citations, said once (the Methods sentence here; run_silac.R's record names the statistics')
PLEXDIA_CITATION = ("Derks J, et al. Increasing the throughput of sensitive proteomics by "
                    "plexDIA. Nat Biotechnol 2023;41:50-59. doi:10.1038/s41587-022-01389-w")


class SilacError(ValueError):
    """A --labels request, a cfg or a design that cannot be searched or analysed as given."""


def _mass(x):
    """4.025107 -> '4.025107', 0.0 -> '0': as DIA-NN's README writes them, never rounded."""
    s = f"{x:.6f}".rstrip("0").rstrip(".")
    return s or "0"


def parse_presets(scheme, presets):
    """--labels/--label-preset -> [(channel name, preset)], lightest first. Raises SilacError."""
    if scheme not in SCHEMES:
        raise SilacError(f"--labels {scheme!r}: use one of {', '.join(SCHEMES)}")
    names = SCHEMES[scheme]
    given = [p.strip().upper() for p in (presets or "").split(",") if p.strip()]
    if not given:
        raise SilacError(
            f"--labels {scheme} needs --label-preset with {len(names) - 1} label(s), one per "
            f"labelled channel ({', '.join(names[1:])}), e.g. "
            + ("K4R6,K8R10" if scheme == "silac3" else "K8R10")
            + f". Known: {', '.join(p for p in PRESETS if p != LIGHT)}. The masses are not "
              "guessed: ask which labelled lysine and arginine the lab used.")
    if len(given) != len(names) - 1:
        raise SilacError(f"--labels {scheme} has {len(names) - 1} labelled channel(s) "
                         f"({', '.join(names[1:])}); --label-preset gave {len(given)}: "
                         f"{','.join(given)}")
    unknown = [p for p in given if p not in PRESETS or p == LIGHT]
    if unknown:
        raise SilacError(f"--label-preset {','.join(unknown)} unknown. Known: "
                         f"{', '.join(p for p in PRESETS if p != LIGHT)}")
    chans = [(names[0], LIGHT)] + list(zip(names[1:], given))
    # DIA-NN: "channels must be listed in such an order that the channel mass shifts are
    # non-decreasing for each modification site" -- and two channels with the same masses are
    # one channel.
    for (n1, p1), (n2, p2) in zip(chans, chans[1:]):
        k1, r1 = (LABELS[x]["delta"] for x in PRESETS[p1])
        k2, r2 = (LABELS[x]["delta"] for x in PRESETS[p2])
        if k2 < k1 or r2 < r1 or (k1, r1) == (k2, r2):
            raise SilacError(f"--label-preset: channel {n2} ({p2}) must be heavier than "
                             f"{n1} ({p1}) on both K and R; list the lighter label first")
    return chans


def channel_entry(name, preset):
    lys, arg = PRESETS[preset]
    return {"name": name, "preset": preset,
            "K": dict(LABELS[lys], label=lys), "R": dict(LABELS[arg], label=arg)}


def channels_flag(channels):
    """The --channels value: 'SILAC,L,KR,0:0;SILAC,M,KR,4.025107:6.020129;...'."""
    return ";".join(f"{LABEL_GROUP},{c['name']},{SITES},{_mass(c['K']['delta'])}:"
                    f"{_mass(c['R']['delta'])}" for c in channels)


def channel_norm_record(choice=None):
    """The sidecar's `channel_norm`: the design question's answer and on whose word -- the
    default is tagged as one (DE-LIMP rule 2)."""
    value = choice or CHANNEL_NORM_DEFAULT
    if value not in CHANNEL_NORM:
        raise SilacError(f"--channel-norm {value!r}: use one of {', '.join(CHANNEL_NORM)}")
    if choice is None:
        source = (f"DEFAULT, not user-confirmed: {CHANNEL_NORM_DESIGN[value]} is assumed. Ask "
                  "whether the channels are independent samples (--channel-norm spec) or old and "
                  "new protein of one culture (--channel-norm run, pulsed SILAC / turnover)")
    else:
        source = f"given (--channel-norm {value}): {CHANNEL_NORM_DESIGN[value]}"
    return {"value": value, "flag": CHANNEL_NORM[value], "given": choice is not None,
            "design": CHANNEL_NORM_DESIGN[value], "source": source,
            "ratio_analysis": value in RATIO_ANALYSIS_NORMS}


def labels_record(scheme, presets, channel_norm=None):
    """--labels/--label-preset/--channel-norm -> the params sidecar's `labels` record: the
    channels with their masses, the DIA-NN declaration and the normalisation choice."""
    chans = [channel_entry(n, p) for n, p in parse_presets(scheme, presets)]
    norm = channel_norm_record(channel_norm)
    given = f"--labels {scheme} --label-preset " + ",".join(c["preset"] for c in chans[1:])
    return {
        "scheme": scheme, "plex": len(chans), "label_group": LABEL_GROUP,
        "channels": chans, "light_channel": chans[0]["name"],
        "source": f"given ({given})", "unimod": UNIMOD_SOURCE,
        "channel_norm": norm,
        "diann": {"--fixed-mod": FIXED_MOD, "--channels": channels_flag(chans),
                  "--original-mods": True, norm["flag"]: True,
                  "--lib-fixed-mod": LIB_FIXED_MOD_NOTE},
        "citation": PLEXDIA_CITATION,
    }


def diann_cfg_lines(rec):
    """[(flag, value, source, render)] for build_diann()'s add(): the declaration, every flag
    tagged with why, --lib-fixed-mod recorded but not written."""
    d = rec["diann"]
    norm = rec["channel_norm"]
    plex = ", ".join(f"{c['name']} {c['preset']}" for c in rec["channels"])
    return [
        ("--fixed-mod", d["--fixed-mod"],
         f"SILAC label declared on K and R as a zero-mass label (DIA-NN README, plexDIA); applied "
         f"during the FASTA digest, so the predicted library carries it once -- {rec['source']}",
         True),
        ("--channels", d["--channels"],
         f"SILAC channels {plex}: masses from {rec['unimod']} -- {rec['source']}", True),
        ("--original-mods", True,
         "keeps DIA-NN from converting the declared label to a UniMod name (README, plexDIA)",
         True),
        (norm["flag"], True, f"channel normalisation, required with --channels (README): "
                             f"{norm['source']}", True),
        ("--lib-fixed-mod", d["--lib-fixed-mod"], LIB_FIXED_MOD_NOTE, False),
    ]


# --- reading a cfg -----------------------------------------------------------------------------
_CHANNEL_RE = re.compile(r"^\s*([^,]+),([^,]+),([^,]+),([0-9.:eE+-]+)\s*$")


def parse_channels(value):
    """'SILAC,L,KR,0:0;SILAC,H,KR,8.014199:10.008269' -> [{group, name, sites, masses}]."""
    out = []
    for part in str(value or "").split(";"):
        if not part.strip():
            continue
        m = _CHANNEL_RE.match(part)
        if not m:
            raise SilacError(f"--channels entry {part.strip()!r} is not group,name,sites,masses")
        group, name, sites, masses = (x.strip() for x in m.groups())
        try:
            ms = [float(x) for x in masses.split(":")]
        except ValueError:
            raise SilacError(f"--channels entry {part.strip()!r}: masses are not numbers")
        letters = [c for c in sites if c.isalpha()]
        if len(ms) != len(letters):
            raise SilacError(f"--channels entry {part.strip()!r}: {len(letters)} site(s) "
                             f"({sites}) but {len(ms)} mass(es)")
        out.append({"group": group, "name": name, "sites": sites, "masses": ms})
    return out


def _flags(groups):
    out = {}
    for flag, vals in groups:
        out.setdefault(flag, []).append(vals)
    return out


def labels_from_cfg(groups):
    """The channels as a cfg declares them (cfg_groups() output) -> a record shaped like
    labels_record()'s, read from the flags alone, or None when the cfg has no --channels. A
    search made before 2.13 (labels passed as --overrides) is described the same way."""
    f = _flags(groups)
    if "--channels" not in f:
        return None
    chans = parse_channels(" ".join(str(v) for v in (f["--channels"][-1] or [])))
    out = []
    for c in chans:
        e = {"name": c["name"], "preset": None}
        for site, mass in zip([x for x in c["sites"] if x.isalpha()], c["masses"]):
            lab = next((k for k, v in LABELS.items()
                        if v["residue"] == site and abs(v["delta"] - mass) < 1e-4), None)
            e[site] = dict(LABELS[lab], label=lab) if lab else {
                "residue": site, "delta": mass, "unimod": None, "composition": None,
                "label": None}
        e["preset"] = next((p for p, (k, r) in PRESETS.items()
                            if e.get("K", {}).get("label") == k
                            and e.get("R", {}).get("label") == r), None)
        out.append(e)
    norm = next((k for k, fl in CHANNEL_NORM.items() if fl in f), None)
    return {"scheme": {2: "silac2", 3: "silac3"}.get(len(out)), "plex": len(out),
            "label_group": chans[0]["group"] if chans else None, "channels": out,
            "light_channel": out[0]["name"] if out else None,
            "source": "read from the search's --channels",
            "channel_norm": ({"value": norm, "flag": CHANNEL_NORM[norm], "given": None,
                              "design": CHANNEL_NORM_DESIGN[norm],
                              "source": f"{CHANNEL_NORM[norm]} in the search's parameters "
                                        "(whether it was the user's answer is not recorded)",
                              "ratio_analysis": norm in RATIO_ANALYSIS_NORMS}
                             if norm else None),
            "lib_fixed_mod": "--lib-fixed-mod" in f}


def cfg_problems(groups, seed_lib=False):
    """Why a DIA-NN cfg (cfg_groups() output) would mis-search a multiplexed sample: [] when it
    is fine or carries no --channels. run_search.py and diann_parallel.py refuse on any."""
    f = _flags(groups)
    probs = []
    norms = [fl for fl in CHANNEL_NORM.values() if fl in f]
    if "--channels" not in f:
        if norms:
            probs.append(f"{' and '.join(norms)} without --channels: no channels to normalise")
        if "--lib-fixed-mod" in f and "--fasta-search" in f:
            probs.append("--lib-fixed-mod with --fasta-search applies the fixed mod twice: "
                         + LIB_FIXED_MOD_NOTE)
        return probs
    try:
        chans = parse_channels(" ".join(str(v) for v in (f["--channels"][-1] or [])))
    except SilacError as e:
        return [str(e)]
    if not chans:
        return ["--channels is empty"]
    if not norms:
        probs.append("--channels without --channel-spec-norm or --channel-run-norm: DIA-NN "
                     "2.7.0 stops with 'channel normalisation strategy must be specified' and "
                     "falls back to --channel-run-norm. Pass estimate_params.py --channel-norm")
    elif len(norms) > 1:
        probs.append("both --channel-spec-norm and --channel-run-norm: choose one")
    if "--original-mods" not in f:
        probs.append("--channels without --original-mods (DIA-NN README, plexDIA: it keeps the "
                     "declared label from being converted to a UniMod name)")
    groups_named = {c["group"] for c in chans}
    declared = {str((v or [""])[0]).split(",")[0].strip()
                for k in ("--fixed-mod", "--mod") for v in f.get(k, [])}
    for g in sorted(groups_named - declared):
        probs.append(f"--channels names the label group {g!r}, which no --fixed-mod declares "
                     f"(e.g. --fixed-mod {g},0.0,{SITES},label)")
    if len(groups_named) > 1:
        probs.append("--channels uses more than one label group; DIA-NN requires one")
    for a, b in zip(chans, chans[1:]):
        if a["sites"] == b["sites"] and any(y < x for x, y in zip(a["masses"], b["masses"])):
            probs.append(f"--channels {b['name']} is lighter than {a['name']} on a site; DIA-NN "
                         "wants each site's masses non-decreasing in the order listed")
    if "--lib-fixed-mod" in f:
        probs.append("--lib-fixed-mod with --channels: " + LIB_FIXED_MOD_NOTE
                     + ". Remove it from the parameters")
    if seed_lib:
        probs.append("a seed library with --channels: v2.13 searches SILAC only from a library "
                     "predicted from the FASTA (the label is then applied once); a label-free "
                     "seed library would need --lib-fixed-mod at exactly the steps that read it")
    return probs


def cfg_labels(cfg_path):
    """(labels record or None, problems) of a cfg FILE, through the one cfg reader."""
    from diann_parallel import cfg_groups, cfg_tokens
    groups = cfg_groups(cfg_tokens(cfg_path))
    return labels_from_cfg(groups), cfg_problems(groups)


def refusal(cfg_path, seed_lib=False, prog="run_search"):
    """The message a search generator exits with for a cfg cfg_problems() objects to, or None.
    A cfg that cannot be read is left to the generator's own reader to report."""
    from diann_parallel import cfg_groups, cfg_tokens, CfgError
    try:
        groups = cfg_groups(cfg_tokens(cfg_path))
    except (CfgError, OSError):
        return None
    probs = cfg_problems(groups, seed_lib=seed_lib)
    if not probs:
        return None
    return (f"[{prog}] {cfg_path} would not search a SILAC (plexDIA) sample correctly: "
            + "; ".join(probs) + ". Nothing was written or submitted. Build the parameters with "
            "estimate_params.py --labels (references/silac.md).")


# --- the Methods ---------------------------------------------------------------------------------
def _label_phrase(site):
    if not site or not site.get("delta"):
        return None
    aa = {"K": "Lys", "R": "Arg"}[site["residue"]]
    comp = f"{site['composition']}, " if site.get("composition") else ""
    um = f"UniMod:{site['unimod']}" if site.get("unimod") else "UniMod accession not matched"
    return f"{aa} +{_mass(site['delta'])} Da ({comp}{um})"


def methods_search_sentence(rec):
    """The Database-search paragraph's plexDIA sentence from a labels record (labels_record() or
    labels_from_cfg()); "" for None."""
    if not rec or not rec.get("channels"):
        return ""
    parts = []
    for c in rec["channels"]:
        lab = [x for x in (_label_phrase(c.get("K")), _label_phrase(c.get("R"))) if x]
        parts.append(f"{c['name']} ({'; '.join(lab) if lab else 'unlabelled Lys and Arg'})")
    norm = rec.get("channel_norm") or {}
    if norm.get("value"):
        tag = {True: "", False: "; DEFAULT, not user-confirmed"}.get(
            norm.get("given"), "; whether this was the user's answer is not recorded")
        norm_s = f" Channels were normalised with {norm['flag']} ({norm['design']}{tag})."
    else:
        norm_s = " The channel-normalisation setting is not recorded."
    return (f"Samples were multiplexed by SILAC and analysed as plexDIA (Derks et al., Nat "
            f"Biotechnol 2023) with {len(parts)} channels: " + ", ".join(parts) + ". The label "
            f"was declared as a zero-mass fixed label on K and R ({rec.get('label_group') or LABEL_GROUP}) "
            "and the channel masses with DIA-NN's --channels; channel-specific confidence was "
            "estimated from a decoy channel (Channel.Q.Value)." + norm_s)


# --- the channel design ----------------------------------------------------------------------------
DESIGN_COLUMNS = ("File.Name", "Replicate", "Channel", "Group")
_GROUP_OK = re.compile(r"^[A-Za-z][A-Za-z0-9_.]*$")


def _read_rows(path):
    with open(path, newline="", encoding="utf-8-sig") as fh:
        text = fh.read()
    dialect = "excel-tab" if text.count("\t") > text.count(",") else "excel"
    rows = list(csv.DictReader(text.splitlines(), dialect=dialect))
    return [{(k or "").strip(): (v or "").strip() for k, v in r.items()} for r in rows]


def parse_comparisons(text):
    out = []
    for c in [x.strip() for x in (text or "").split(",") if x.strip()]:
        a, sep, b = c.partition("-")
        if not sep or not a.strip() or not b.strip() or "-" in b:
            raise SilacError(f"comparison {c!r}: write it as <numerator>-<denominator>, "
                             "e.g. KO-WT for log2(KO/WT)")
        out.append((a.strip(), b.strip()))
    return out


def read_design(path, channels=None, comparisons=None, report_runs=None):
    """silac_design.csv -> the normalised design run_silac.R analyses, or SilacError naming
    every problem. One row per (run, channel) that holds a sample:

        File.Name,Replicate,Channel,Group
        run_01,Rep1,L,WT
        run_01,Rep1,H,KO          (a label swap is just another assignment: run_02,Rep2,L,KO ...)

    File.Name is the run as the report names it (its Run column; conditions.csv's convention).
    Replicate defaults to File.Name. A channel that holds no sample has no row. A run with one
    labelled channel only is a label-check run (incorporation is measured on it, never compared).
    `channels`: the search's channel names, lightest first (from the params; else the design's
    own, in L, M, H order). `comparisons`: 'A-B,...' or None for every pair, the group in the
    heavier channel over the lighter (DEFAULT, said so). `report_runs`: the report's runs."""
    rows = _read_rows(path)
    probs = []
    if not rows:
        raise SilacError(f"{path}: no rows")
    cols = set(rows[0])
    missing = [c for c in ("File.Name", "Channel", "Group") if c not in cols]
    if missing:
        raise SilacError(f"{path}: missing column(s) {', '.join(missing)}; the columns are "
                         f"{','.join(DESIGN_COLUMNS)} (Replicate optional)")
    order = list(channels) if channels else None
    if order is None:
        seen = []
        for r in rows:
            if r["Channel"] and r["Channel"] not in seen:
                seen.append(r["Channel"])
        known = ["L", "M", "H"]
        order = sorted(seen, key=lambda c: (known.index(c) if c in known else 99, c))
    runs, rep_of, groups = {}, {}, []
    for i, r in enumerate(rows, start=2):
        run, ch, g = r.get("File.Name", ""), r.get("Channel", ""), r.get("Group", "")
        rep = r.get("Replicate") or run
        if not run or not ch or not g:
            probs.append(f"line {i}: File.Name, Channel and Group are all needed")
            continue
        if ch not in order:
            probs.append(f"line {i}: channel {ch!r} is not one of the search's channels "
                         f"({', '.join(order)})")
            continue
        if not _GROUP_OK.match(g):
            probs.append(f"line {i}: group {g!r} -- use letters, digits, '_' or '.', starting "
                         "with a letter (it names columns and comparisons)")
        e = runs.setdefault(run, {"replicate": rep, "channels": {}})
        if e["replicate"] != rep:
            probs.append(f"line {i}: run {run} is replicate {e['replicate']!r} on another line")
        if ch in e["channels"]:
            probs.append(f"line {i}: run {run} channel {ch} is assigned twice")
        if g in e["channels"].values():
            probs.append(f"line {i}: group {g} is in two channels of run {run}")
        e["channels"][ch] = g
        if rep in rep_of and rep_of[rep] != run:
            probs.append(f"line {i}: replicate {rep!r} names two runs ({rep_of[rep]}, {run}); "
                         "one replicate is one run")
        rep_of[rep] = run
        if g not in groups:
            groups.append(g)
    if report_runs is not None:
        rr = set(report_runs)
        for run in runs:
            if run not in rr:
                probs.append(f"run {run!r} is not in the report (its runs: "
                             f"{', '.join(sorted(rr)[:6])}{' ...' if len(rr) > 6 else ''})")
        left = sorted(rr - set(runs))
        if left:
            probs.append(f"report run(s) with no row: {', '.join(left)} -- add them, or "
                         "remove them from the search")
    mixed = {r: e for r, e in runs.items() if len(e["channels"]) >= 2}
    checks = {r: e for r, e in runs.items() if len(e["channels"]) == 1
              and next(iter(e["channels"])) != order[0]}
    unused = sorted(r for r, e in runs.items() if len(e["channels"]) == 1 and r not in checks)
    if len(groups) < 2 or not mixed:
        probs.append("no run holds two groups: a SILAC ratio needs two channels of one run")
    # comparisons
    given = comparisons is not None and bool(str(comparisons).strip())
    try:
        pairs = parse_comparisons(comparisons) if given else []
    except SilacError as e:
        probs.append(str(e))
        pairs = []
    if not given:
        rank = {}
        for e in mixed.values():
            for ch, g in e["channels"].items():
                rank.setdefault(g, order.index(ch))
        for lo, hi in itertools.combinations(sorted(rank, key=lambda g: (rank[g], groups.index(g))), 2):
            pairs.append((hi, lo))
        pairs.sort(key=lambda p: (rank[p[1]], rank[p[0]]))
    out_pairs = []
    for a, b in pairs:
        if a == b or a not in groups or b not in groups:
            probs.append(f"comparison {a}-{b}: both groups must be in the design ({', '.join(groups)})"
                         if a != b else f"comparison {a}-{b} compares a group with itself")
            continue
        n = sum(1 for e in mixed.values() if a in e["channels"].values()
                and b in e["channels"].values())
        if not n:
            probs.append(f"comparison {a}-{b}: no run holds both")
            continue
        out_pairs.append({"contrast": f"{a}-{b}", "numerator": a, "denominator": b,
                          "n_runs": n})
    if probs:
        raise SilacError(f"{path}: " + "; ".join(probs))
    swaps = sorted({g for g in groups
                    if len({ch for e in mixed.values() for ch, gg in e["channels"].items()
                            if gg == g}) > 1})
    return {
        "design_file": os.path.abspath(path), "channels": order, "groups": groups,
        "runs": [{"run": r, "replicate": e["replicate"], "channels": e["channels"],
                  "role": "ratio" if r in mixed else "label_check" if r in checks else "unused"}
                 for r, e in runs.items()],
        "label_swapped_groups": swaps,
        "comparisons": out_pairs,
        "comparisons_source": ("given (--comparisons)" if given else
                               "DEFAULT, not user-confirmed: every pair of groups, the group in "
                               "the heavier channel over the lighter"),
        "label_check_runs": sorted(checks), "unused_runs": unused,
    }


def report_runs(report):
    """The runs of a DIA-NN report (its Run column), read with pyarrow."""
    import pyarrow.parquet as pq
    t = pq.read_table(report, columns=["Run"])
    return sorted(set(t.column("Run").to_pylist()))


# --- CLI ---------------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("presets", help="list the label presets and their masses")
    f = sub.add_parser("flags", help="the DIA-NN flags --labels writes (a preview; "
                                     "estimate_params.py writes the cfg)")
    f.add_argument("--labels", required=True, choices=sorted(SCHEMES))
    f.add_argument("--label-preset", required=True)
    f.add_argument("--channel-norm", choices=sorted(CHANNEL_NORM))
    c = sub.add_parser("check-cfg", help="what would make a cfg mis-search a SILAC sample")
    c.add_argument("cfg")
    d = sub.add_parser("design", help="validate silac_design.csv (one row per run and channel)")
    d.add_argument("--validate", required=True, metavar="CSV")
    d.add_argument("--against", help="the report: every run in it needs rows, and every row a run")
    d.add_argument("--channels", help="the search's channels, lightest first, e.g. L,M,H")
    d.add_argument("--comparisons", help="'A-B,C-B' (A over B); default every pair")
    d.add_argument("--json", action="store_true", help="print the normalised design as JSON")
    a = ap.parse_args(argv)

    if a.cmd == "presets":
        for p, (k, r) in PRESETS.items():
            print(f"{p:6s} K {_mass(LABELS[k]['delta']):>9s}  R {_mass(LABELS[r]['delta']):>9s}"
                  + ("   (light)" if p == LIGHT else ""))
        print(f"masses: {UNIMOD_SOURCE}")
        return 0
    if a.cmd == "flags":
        try:
            rec = labels_record(a.labels, a.label_preset, a.channel_norm)
        except SilacError as e:
            print(f"[silac] {e}", file=sys.stderr)
            return 2
        for flag, value, source, render in diann_cfg_lines(rec):
            print(("# not written: " if not render else "")
                  + (flag if value is True else f"{flag} {value}"))
        return 0
    if a.cmd == "check-cfg":
        try:
            rec, probs = cfg_labels(a.cfg)
        except Exception as e:                       # noqa: BLE001 -- said, with its kind
            print(f"[silac] {a.cfg} could not be read ({type(e).__name__}: {e})", file=sys.stderr)
            return 2
        print(json.dumps({"labels": rec, "problems": probs}, indent=2))
        return 1 if probs else 0
    # design
    try:
        runs = report_runs(a.against) if a.against else None
        chans = [x.strip() for x in a.channels.split(",")] if a.channels else None
        des = read_design(a.validate, chans, a.comparisons, runs)
    except SilacError as e:
        if a.json:
            print(json.dumps({"ok": False, "error": str(e)}))
        print(f"[silac] {e}", file=sys.stderr)
        return 1
    except (OSError, ImportError) as e:
        if a.json:
            print(json.dumps({"ok": False, "error": f"{type(e).__name__}: {e}"}))
        print(f"[silac] could not read the design or the report ({type(e).__name__}: {e})",
              file=sys.stderr)
        return 1
    if a.json:
        print(json.dumps(dict(des, ok=True), indent=2))
    else:
        print(f"{len(des['runs'])} runs, groups {', '.join(des['groups'])}; comparisons "
              f"{', '.join(c['contrast'] for c in des['comparisons'])} ({des['comparisons_source']})")
        if des["label_check_runs"]:
            print(f"label-check runs (incorporation measured, not compared): "
                  f"{', '.join(des['label_check_runs'])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
