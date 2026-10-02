#!/usr/bin/env python3
"""
experiment_type.py -- what kind of experiment this is, decided in ONE place, and what that means
for the DE's quantities: DIA-NN's normalised ones, or its non-normalised ones.

Why (2.10, Brett): DIA-NN's cross-run normalisation -- like any global normalisation -- assumes
most proteins do not change between samples. An IP / pull-down violates that: its IgG or bead
controls carry little protein, so normalisation scales them up several-fold and hides the
enrichment. The pieces existed in 2.9 and nothing connected them: the CoreOmics record's
experiment type and "Normalization" answer were only displayed, the pull-down flag came from
group names alone (analysis_prompt: contrasts against a group called IgG/beads), DIA-NN's
--no-norm was never offered, and the DE always used normalised input (the maxlfq path added
quantile normalisation on top). Now:

  * the type is PROPOSED from the attached submission (its proteomics_type), never from file
    names, and CONFIRMED with the user -- or asked when there is no record -- and recorded in
    <session>/input/experiment_type.json, with where it came from;
  * each type has a default for the DE's quantities (DEFAULTS, below) -- a default, applied
    only after normalization_check.py's data check, which runs for EVERY type and decides
    whether to stop and ask;
  * the pull-down flag the analysis uses (pulldown_design()) comes from the same record, so the
    brief and the DE cannot disagree about the design;
  * the submitter's "Normalization" answer (how the samples were loaded) is set beside the
    default, and a contradiction is raised as a question (reconcile()).

  python3 experiment_type.py propose --session <S>      # from the submission; ask the user to confirm
  python3 experiment_type.py set --session <S> --type ip --source submission|user \\
          [--stated "<the user's words>"] [--bait <gene>] [--controls IgG_A,IgG_B]
  python3 experiment_type.py show --session <S>

The bait: --bait, else a sequence the user added to the database (fetch_fasta.py --add-fasta,
read from the session's FASTA record) -- bait_candidates() / bait_for_check(), the one rule the
data check uses.

Stdlib only.
"""
import argparse
import datetime
import json
import os
import re
import sys

FILE = "experiment_type.json"
SCHEMA = "experiment_type/1"
# The session's copy of fetch_fasta.py's <fasta>.meta.json (session.py's "fasta_meta").
FASTA_META = os.path.join("input", "search.fasta.meta.json")

# The group names that read as a pull-down's controls. THE one definition: session_docs.py,
# analysis_prompt.py and normalization_check.py import it from here.
IP_CONTROL_NAME = re.compile(r"(?i)(^|[_\-. ])(igg|beads?)($|[_\-. ])")

RAW, NORMALISED = "raw", "normalised"
# type -> (label, default quantities, why). The defaults are Brett's (Core director, 2026-10-01):
# a classic IP's controls carry little protein and are not loaded by amount, so normalisation
# OFF; proximity labelling's controls carry plenty of signal and the samples are usually loaded
# by equal protein ("for TurboID it seems not to matter so much"), so ON; secretome ON but
# leaning on the data check. Fractions are two classes, by what is compared:
#   fractions_depth       offline high-pH (or similar) fractions COMBINED per sample, to see more
#                         of one proteome: every sample is still the whole proteome -- ON;
#   fractions_separation  fractions COMPARED with each other (SEC, density / sucrose gradient,
#                         BN-PAGE complexome, organelle / subcellular): each fraction holds a
#                         different part of the proteome -- OFF. DIA-NN's README (2.7.0, FAQ "What is
#                         normalisation and how does it work?"): global normalisation assumes most
#                         peptides are not differentially abundant, and "If this condition is not
#                         satisfied (e.g. when analysing fractions obtained with some separation
#                         technique, like SEC, for instance), then normalisation should not be
#                         used"; and ("Disabling normalisation") "any kind of protein
#                         fractionation" among the scenarios to disable it for.
# The data check runs for both.
DEFAULTS = {
    "whole_proteome": ("whole-proteome (expression) comparison", NORMALISED,
                       "most proteins are expected not to change between the groups, the "
                       "assumption cross-run normalisation rests on"),
    "ip": ("immunoprecipitation / affinity purification / pull-down with IgG or bead controls",
           RAW,
           "the controls carry little protein and are not loaded by amount, so normalisation "
           "would scale them up several-fold and hide the enrichment"),
    "proximity": ("proximity labelling (TurboID / BioID / APEX)", NORMALISED,
                  "the controls carry plenty of signal and the samples are usually loaded by "
                  "equal protein, so standard normalisation holds -- unless the data check finds "
                  "the controls came out low"),
    "secretome": ("secretome / conditioned medium", NORMALISED,
                  "standard, but secreted output can differ between conditions, so the data "
                  "check decides"),
    "fractions_depth": ("fractionated for depth (offline high-pH or similar, fractions combined "
                        "per sample)", NORMALISED,
                        "every sample is still the whole proteome, only seen in more depth, so "
                        "standard normalisation holds"),
    "fractions_separation": ("separation / profiling fractions compared with each other (SEC, "
                             "density or sucrose gradient, BN-PAGE complexome, organelle or "
                             "subcellular fractionation)", RAW,
                             "each fraction holds a different part of the proteome, so the "
                             "assumption normalisation rests on fails -- DIA-NN's README: for "
                             "fractions from a separation technique such as SEC, \"normalisation "
                             "should not be used\""),
    "other": ("other", NORMALISED, "standard; the data check decides"),
}
ENRICHMENT = ("ip", "proximity")
# Proximity labelling only: the endogenously biotinylated carboxylases, which streptavidin
# captures in every sample, so they should be roughly constant across samples -- a reference
# check of a normalisation, never a normaliser. Verified 2026-10-01: UniProt reviewed entries
# with keyword KW-0092 (Biotin) are exactly these carboxylases (plus the transporter SLC5A6) in
# human (P11498, P05165, Q96RQ3, Q13085, O00763) and in mouse (Pc, Pcca, Mccc1, Acaca, Acacb).
CARBOXYLASES = ("PC", "PCCA", "MCCC1", "ACACA", "ACACB")

# The CoreOmics form's proteomics_type answers -> a type to PROPOSE (first match wins: a record
# naming both "Affinity Purification" and "BioID" is proximity labelling). Unmatched -> ask.
SUBMISSION_PATTERNS = (
    ("proximity", re.compile(r"turbo|bio-?id\b|\bapex|proximity", re.I)),
    ("ip", re.compile(r"affinity|immuno-?precip|\bco-?ip\b|\bip\b|pull[- ]?down|ap-?ms|\bbeads?\b",
                      re.I)),
    ("secretome", re.compile(r"secretom|conditioned med|supernatant", re.I)),
    # a separation compares its fractions; high-pH fractions are combined per sample. A bare
    # "fractionation" says neither, so it proposes nothing: ask
    ("fractions_separation", re.compile(r"\bsec\b|size.exclusion|gradient|sucrose|complexome|"
                                        r"bn-?page|organell|subcellular|co-?fractionation", re.I)),
    ("fractions_depth", re.compile(r"high.?ph|offline|fractionat\w* for depth|deep proteome",
                                   re.I)),
    ("whole_proteome", re.compile(r"global|whole|total proteome|expression|discovery", re.I)),
)
# The form's "Normalization" answer is how the SUBMITTER normalised the samples they sent.
LOADED_BY = (("protein amount", re.compile(r"mass|protein|amount|\bu?g\b|µg|microgram|bca", re.I)),
             ("volume", re.compile(r"volume|\bu?l\b|µl|cell (count|number)", re.I)))


def path_for(session):
    return os.path.join(os.path.abspath(os.path.expanduser(session)), "input", FILE)


def propose_from_submission(types):
    """(type or None, the answer it came from) for the record's proteomics_type answers."""
    for t, pat in SUBMISSION_PATTERNS:
        for answer in types or []:
            if pat.search(str(answer)):
                return t, answer
    return None, None


def loaded_by(answer):
    """'protein amount' | 'volume' | None (not said, or "no idea") for the form's answer."""
    for what, pat in LOADED_BY:
        if answer and pat.search(str(answer)):
            return what
    return None


def reconcile(etype, submission_answer):
    """A question for the user when the submitter's loading contradicts the default, else None."""
    if etype not in DEFAULTS or submission_answer is None:
        return None
    how = loaded_by(submission_answer)
    q = DEFAULTS[etype][1]
    if how is None:
        return (f"The submission's Normalization answer is \"{submission_answer}\": how the "
                f"samples were loaded is not known. Confirm it -- the default for this type "
                f"({'non-normalised' if q == RAW else 'normalised'} quantities) assumes "
                + ("the controls were NOT loaded to equal protein." if q == RAW else
                   "the samples were loaded to equal protein."))
    if q == NORMALISED and how == "volume":
        return (f"The submission says the samples were normalised by volume (\"{submission_answer}\"), "
                "but normalised quantities assume equal protein was loaded. Ask whether the "
                "groups carried similar amounts of protein; the data check will also show it.")
    if q == RAW and how == "protein amount":
        return (f"The submission says the samples were loaded by protein amount (\"{submission_answer}\"). "
                "Loading every sample (an IP's controls, each fraction) to equal protein already "
                "equalises their total signal, so non-normalised quantities would carry that "
                "loading normalisation. Ask how the samples were loaded before choosing.")
    return None


def default_for(etype):
    """{"quantities", "why", "label"} for a type; the standard default, tagged, when unknown."""
    if etype in DEFAULTS:
        label, q, why = DEFAULTS[etype]
        return {"quantities": q, "why": why, "label": label}
    return {"quantities": NORMALISED, "label": "NOT RECORDED",
            "why": "DEFAULT -- not user-confirmed: no experiment type was recorded, so the "
                   "standard (normalised) quantities; the data check decides whether to ask"}


def _kept_bait_answer(session, bait):
    """{"bait", "bait_confirmed"} from an earlier record_bait answer, when this `set` names no
    --bait; {} otherwise."""
    if bait or not os.path.isfile(path_for(session)):
        return {}
    prior = load(session) or {}
    ans = prior.get("bait_confirmed")
    return {"bait": prior.get("bait"), "bait_confirmed": ans} if ans else {}


def record(session, etype, source, stated=None, bait=None, controls=None, submission=None,
           when=None):
    """Write the session's record -> its path. `source`: "submission" (proposed from the record
    and confirmed by the user) or "user"."""
    if etype not in DEFAULTS:
        raise ValueError(f"--type must be one of {', '.join(DEFAULTS)}")
    if source not in ("submission", "user"):
        raise ValueError("--source must be submission (confirmed by the user) or user")
    # The submission's ASK (its normalisation vs the type's default) needs the record: without
    # one, "the CoreOmics submission" would be claimed as the source of nothing (2.10 safety review).
    if source == "submission" and not submission:
        raise ValueError("--source submission, but no submission is attached to this session: "
                         "attach it first (submission_report.py attach), or record --source user")
    p = path_for(session)
    if not os.path.isdir(os.path.dirname(p)):
        raise ValueError("not a session directory (no input/): "
                         f"{os.path.dirname(os.path.dirname(p))}")
    sub = submission or {}
    rec = {"schema": SCHEMA, "type": etype, "label": DEFAULTS[etype][0],
           "source": ("the CoreOmics submission, confirmed by the user" if source == "submission"
                      else "the user"),
           "stated": " ".join(stated.split()) if stated else None,
           "submission_types": sub.get("experiment_types") or [],
           "submission_normalisation": sub.get("normalisation"),
           "default": default_for(etype),
           "reconcile": reconcile(etype, sub.get("normalisation")) if sub else None,
           "bait": bait or None,
           "controls": [c.strip() for c in (controls or "").split(",") if c.strip()] or None,
           # the user's answer to ask_bait (record_bait) survives a re-recorded type
           **_kept_bait_answer(session, bait),
           "recorded_at": when or datetime.datetime.now(datetime.timezone.utc).strftime(
               "%Y-%m-%dT%H:%M:%SZ")}
    with open(p + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=2)
        fh.write("\n")
    os.replace(p + ".tmp", p)
    return p


def load(session):
    """The record, or None (never asked, or a session from before 2.10). Unreadable raises."""
    if not session:
        return None
    p = path_for(session)
    if not os.path.isfile(p):
        return None
    with open(p, encoding="utf-8") as fh:
        rec = json.load(fh)
    if not isinstance(rec, dict) or rec.get("schema") != SCHEMA or rec.get("type") not in DEFAULTS:
        raise ValueError(f"{p} is not a {SCHEMA} record")
    return rec


def session_of(*paths, depth=4):
    """The analysis session a DE's files belong to: the nearest folder, at most `depth` levels up
    from any of `paths` (its --metadata <S>/input/conditions.csv, its --input
    <S>/output/search/report.parquet), that holds an experiment-type record. None when there is
    none -- a DE outside a session, or a session whose type was never recorded."""
    for p in paths:
        if not p:
            continue
        d = os.path.dirname(os.path.abspath(p))
        for _ in range(depth):
            if os.path.isfile(os.path.join(d, "input", FILE)):
                return d
            parent = os.path.dirname(d)
            if parent == d:
                break
            d = parent
    return None


def bait_candidates(session):
    """The sequences the user added to the search database (fetch_fasta.py --add-fasta: a bait
    such as EGFP, a tag), as the candidates for the bait: [{"accession", "name", "file"}], the
    accession being the one the report gives the protein (fetch_fasta records it). [] when there
    is no session, no FASTA record, or no added sequence; an unreadable record is [] too -- the
    bait is then asked, as before."""
    if not session:
        return []
    try:
        with open(os.path.join(session, FASTA_META), encoding="utf-8") as fh:
            meta = json.load(fh)
    except (OSError, ValueError):
        return []
    out = []
    for f in (meta.get("added_sequences") or []) if isinstance(meta, dict) else []:
        for e in (f.get("entries") or []) if isinstance(f, dict) else []:
            if isinstance(e, dict) and e.get("accession"):
                out.append({"accession": e["accession"], "name": e.get("name") or e["accession"],
                            "file": f.get("file")})
    return out


# The answer to ask_bait that NONE of the added sequences is the bait (a tag, a spike-in).
NO_BAIT = "none"


def record_bait(session, bait, when=None):
    """The user's answer to ask_bait (asked after step 6, once the database is built): the bait's
    accession, or NO_BAIT. Written into the session's record with who answered and when, so the
    one added sequence is never used as an unconfirmed bait (2.10 safety review). -> path."""
    rec = load(session)
    if rec is None:
        raise ValueError("no experiment type recorded for this session: record it first "
                         "(experiment_type.py set), then the bait")
    answer = " ".join(str(bait or "").split())
    if not answer:
        raise ValueError("--bait needs the bait's accession, or 'none'")
    cands = [c["accession"] for c in bait_candidates(session)]
    if answer.casefold() != NO_BAIT and cands and answer not in cands:
        raise ValueError(f"--bait {answer}: not one of the added sequences ({', '.join(cands)}); "
                         f"give one of them, or 'none'")
    none = answer.casefold() == NO_BAIT
    rec.update(bait=None if none else answer,
               bait_confirmed={"answer": NO_BAIT if none else answer,
                               "by": "the user, asked ask_bait after the database was built",
                               "at": when or datetime.datetime.now(datetime.timezone.utc).strftime(
                                   "%Y-%m-%dT%H:%M:%SZ")})
    p = path_for(session)
    with open(p + ".tmp", "w", encoding="utf-8") as fh:
        json.dump(rec, fh, indent=2)
        fh.write("\n")
    os.replace(p + ".tmp", p)
    return p


def bait_for_check(session, bait=None, rec=None):
    """THE bait the data check uses, and where it came from: (bait or None, source, candidates).
    Given (--bait) wins, then the recorded one (experiment_type.py set --bait / set-bait), then
    the user's recorded answer that none of the added sequences is the bait, then -- only when
    the user added exactly ONE sequence and was never asked -- that sequence, said to be
    unconfirmed. Several added sequences are listed, never chosen between."""
    cands = bait_candidates(session)
    if bait:
        return bait, "--bait", cands
    if (rec or {}).get("bait"):
        return rec["bait"], "the experiment-type record", cands
    if ((rec or {}).get("bait_confirmed") or {}).get("answer") == NO_BAIT:
        return None, "the user said none of the added sequences is the bait", cands
    if len(cands) == 1:
        return (cands[0]["accession"], "the one sequence the user added to the database "
                "(--add-fasta), not confirmed as the bait", cands)
    return None, ("not given, and the user added several sequences (--add-fasta): record "
                  "which is the bait with experiment_type.py set --bait" if cands else
                  "not given"), cands


def controls_in(groups, rec=None):
    """The control groups: the recorded --controls, else the groups named like IgG / beads."""
    named = (rec or {}).get("controls")
    return [g for g in groups if g in named] if named else \
        [g for g in groups if IP_CONTROL_NAME.search(g)]


def pulldown_design(rec, groups, contrasts):
    """THE pull-down decision the analysis uses: {"pulldown", "source", "controls", "vs_control",
    "note"}. From the experiment type when recorded -- an enrichment type is a pull-down design --
    else from the group names, said so; a disagreement between the two is a note to raise."""
    controls = controls_in(groups, rec)
    vs_control = [c for c in contrasts if "-" in c and c.split("-", 1)[1].strip() in controls]
    if not rec:
        return {"pulldown": bool(vs_control), "controls": controls, "vs_control": vs_control,
                "source": "group names only (no experiment type recorded)", "note": None}
    enr = rec["type"] in ENRICHMENT
    note = None
    if enr and not vs_control:
        note = (f"the experiment type is {rec['label']}, but no contrast is against a control "
                "group (none is named like IgG/beads, and no --controls were recorded): "
                "record the control groups (experiment_type.py set --controls ...)")
    elif not enr and vs_control:
        note = (f"the groups {', '.join(controls)} read as a pull-down's controls, but the "
                f"experiment type is {rec['label']}: confirm the type with the user")
    return {"pulldown": enr, "controls": controls, "vs_control": vs_control,
            "source": f"experiment type ({rec['source']})", "note": note}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("propose", help="the type the attached submission suggests -- confirm it")
    p.add_argument("--session", required=True)
    s = sub.add_parser("set", help="record the confirmed type")
    s.add_argument("--session", required=True)
    s.add_argument("--type", required=True, choices=sorted(DEFAULTS))
    s.add_argument("--source", required=True, choices=("submission", "user"))
    s.add_argument("--stated", help="the user's own words for the experiment")
    s.add_argument("--bait", help="the bait / fusion protein (gene or accession), when known")
    s.add_argument("--controls", help="comma-separated control groups, when not named IgG/beads")
    sb = sub.add_parser("set-bait", help="record the user's answer to ask_bait (after step 6)")
    sb.add_argument("--session", required=True)
    sb.add_argument("--bait", required=True,
                    help="the bait's accession (one of the added sequences), or 'none'")
    sh = sub.add_parser("show")
    sh.add_argument("--session", required=True)
    a = ap.parse_args()
    try:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import submission_report as sr
        try:
            submission = sr.load(a.session)
        except sr.RecordError as e:     # an unreadable submission is said, never skipped
            submission = None
            print(f"[experiment_type] the submission record could not be read: {e}",
                  file=sys.stderr)
        if a.cmd == "propose":
            t, answer = propose_from_submission((submission or {}).get("experiment_types"))
            cands = bait_candidates(a.session)
            print(json.dumps({
                "proposed": t, "from": answer,
                "submission_types": (submission or {}).get("experiment_types"),
                "submission_normalisation": (submission or {}).get("normalisation"),
                "default": default_for(t) if t else None,
                "reconcile": reconcile(t, (submission or {}).get("normalisation")) if t else None,
                "ask": ("Confirm with the user: is this " + DEFAULTS[t][0] + "?") if t else
                       "No submission type to go on: ask the user what kind of experiment this "
                       "is (" + "; ".join(v[0] for v in DEFAULTS.values()) + ").",
                # the sequences the user added (--add-fasta): a bait, if this is an
                # enrichment; confirm which one with set --bait
                "bait_candidates": cands,
                "ask_bait": (("The search database has sequence(s) the user added: "
                              + ", ".join(c["accession"] for c in cands) + ". If this is an "
                              "enrichment, is one of them the bait? Record the answer with "
                              "set-bait --bait <accession> (or --bait none).")
                             if cands else None),
                "types": {k: v[0] for k, v in DEFAULTS.items()}}, indent=2))
        elif a.cmd == "set-bait":
            path = record_bait(a.session, a.bait)
            print(json.dumps({"written": path, "record": load(a.session)}, indent=2))
        elif a.cmd == "set":
            path = record(a.session, a.type, a.source, a.stated, a.bait, a.controls, submission)
            rec = load(a.session)
            print(json.dumps({"written": path, "record": rec}, indent=2))
            if rec.get("reconcile"):
                print(f"[experiment_type] ASK: {rec['reconcile']}", file=sys.stderr)
        else:
            rec = load(a.session)
            print(json.dumps({"record": rec, "default": default_for((rec or {}).get("type"))},
                             indent=2))
    except (ValueError, OSError) as e:
        sys.exit(f"[experiment_type] {e}")


if __name__ == "__main__":
    main()
