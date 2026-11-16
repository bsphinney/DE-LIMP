#!/usr/bin/env python3
"""
cleavage_rules.py -- THE reading of a search's in-silico cleavage, from the engine's own record
(2.14.2, L11 N1/N4/N5). make_methods.py writes the Methods and the SDRF cleavage agent from it;
run_search.py sizes a non-specific search's job from it. Stdlib only.

  * sage_cleavage(enzyme block)        -- Sage's database.enzyme, as Sage 0.14.7 reads it
  * fragger_cleavage(fragpipe_params.digest(...)) -- MSFragger's fragger.params / workflow keys
  * cleavage_agent(cut, no_cut, c_terminal) -- the PSI-MS term of one rule, or None
  * enzyme_rule(--enzyme)              -- the rule a search uses for the enzyme(s) a lab added, in
                                          fetch_fasta.py's words (2.15): estimate_params.py writes
                                          it into a Sage cfg, make_presets.py into a FragPipe preset
  * fragpipe_digest(path)              -- a FragPipe workflow's / fragger.params' digest in
                                          estimate_params.DIANN_DIGEST's shape (contaminant check)
  * fragpipe_dropdown_conflict(props)  -- what a FragPipe enzyme dropdown would change in a rule:
                                          make_presets.py refuses its OWN preset on it; a user's or
                                          a vendor's workflow is never refused on it (Ruling 104)

Each of the first two returns {"rule", "name", "ac", "specificity", "source"}: `specificity` is
"none" (a non-specific search), "no_digestion" (Sage "$"), "semi" or "full"; `name`/`ac` the
PSI-MS cleavage term when one fits, else None (the caller tags the rule).
"""
import os
import re

# In-silico cleavage rule -> PSI-MS cleavage agent (name, accession). One rule is (the residues it
# cuts at, as a sorted set; the residue after which it does NOT cut, "" for none; True when it cuts
# C-terminal of them). Sage spells it cleave_at + restrict + c_terminal (Sage DOCS.md v0.14.7:
# restrict = "do not cleave if this AA follows the cleavage site"; cleave_at "FWYL" for
# chymotrypsin; enzyme.rs reads cleave_at as a character class, so its order does not matter);
# MSFragger spells it search_enzyme_cut/nocut/sense (sense C = C-terminal). Names and accessions:
# psi-ms.obo 4.2.14 (HUPO-PSI/psi-ms-CV, read 2026-11-03), children of MS:1001045 "cleavage agent
# name", each term's has_regexp. Residue sets checked against FragPipe 24.0's enzyme list
# (FragPipe-GUI .../tools/fragger/enzyme-defs-supported.txt @ tag 24.0) and Skyline 26.1's
# predefined enzymes (skyline_peaks.SKYLINE_ENZYMES); each agreeing source is named per row.
CLEAVAGE_AGENTS = {
    # (?<=[KR])(?!P) -- FragPipe trypsin, Skyline Trypsin
    ("KR", "P", True): ("Trypsin", "MS:1001251"),
    # (?<=[KR]) -- FragPipe stricttrypsin, Skyline Trypsin/P
    ("KR", "", True): ("Trypsin/P", "MS:1001313"),
    # (?<=[FYWL])(?!P) -- FragPipe chymotrypsin (FLWY), Skyline Chymotrypsin (FWYL), Sage's docs
    ("FLWY", "P", True): ("Chymotrypsin", "MS:1001306"),
    # ...and its low-specificity form, which also cuts after Met (FWYLM, not before P): the same
    # agent. PSI-MS has no separate term; the rule itself stays in the record ("rule") and in the
    # Methods' parameter table, so what was searched is still said.
    ("FLMWY", "P", True): ("Chymotrypsin", "MS:1001306"),
    # (?<=K)(?!P) -- FragPipe lysc, Skyline LysC
    ("K", "P", True): ("Lys-C", "MS:1001309"),
    # (?<=K) -- FragPipe lysc-p, Skyline LysC/P
    ("K", "", True): ("Lys-C/P", "MS:1001310"),
    # (?<=R)(?!P) -- FragPipe argc, Skyline ArgC
    ("R", "P", True): ("Arg-C", "MS:1001303"),
    # (?=[BD]): N-terminal of Asp (B = Asx) -- FragPipe aspn (sense N), Skyline AspN
    ("D", "", False): ("Asp-N", "MS:1001304"),
    # Lys-N cuts N-terminal of Lys -- FragPipe lysn (sense N), Skyline LysN. (psi-ms.obo 4.2.14
    # gives this term the regexp (?<=K)(?!P), Lys-C's; the term is matched by its name here.)
    ("K", "", False): ("Lys-N", "MS:1003093"),
    # Glu-C. (?<=[EZ])(?!P): in bicarbonate, after Glu -- FragPipe gluc_bicarb, Skyline GluC bicarb
    ("E", "P", True): ("V8-E", "MS:1001315"),
    # (?<=[BDEZ])(?!P): in phosphate, after Asp and Glu -- FragPipe gluc, Skyline GluC
    ("DE", "P", True): ("V8-DE", "MS:1001314"),
    # (?<=[FYWLKR])(?!P) -- FragPipe trypsin/chymotrypsin
    ("FKLRWY", "P", True): ("TrypChymo", "MS:1001312"),
}
# The table under the name make_methods.py imports it by (before 2.14.2 it held trypsin only)
SAGE_CUT = CLEAVAGE_AGENTS
# A digest with no enzyme rule (psi-ms.obo, read 2026-11-03): every peptide in the length range
# (Sage cleave_at "", MSFragger num_enzyme_termini 0), or each database entry as it is (Sage "$").
NONSPECIFIC = ("unspecific cleavage", "MS:1001956")
NO_CLEAVAGE = ("no cleavage", "MS:1001955")
# The specificities a search is sized as broad for (run_search.single_job_sizing)
BROAD = {"none": "non-specific", "semi": "semi-specific"}
# What a reader of a semi-specific search's numbers must know (C2 bio M3): said once here, and
# quoted by the Methods (make_methods.search_paragraph), the report's callout
# (make_analysis_html.semi_specific_note) and the analysis brief.
SEMI_SPECIFIC_CAVEAT = ("Peptides with one non-specific end were allowed (a semi-specific "
                        "search), so the peptide counts are not comparable with a fully specific "
                        "search, and some semi-specific peptides may be in-source fragments "
                        "rather than products of the digestion.")


def semi_specific_caveat(cleavage):
    """SEMI_SPECIFIC_CAVEAT when a search record's cleavage (make_methods.search_record
    `cleavage`) is semi-specific, else None."""
    return SEMI_SPECIFIC_CAVEAT if isinstance(cleavage, dict) and \
        cleavage.get("specificity") == "semi" else None


def cleavage_agent(cut, no_cut, c_terminal=True):
    """The PSI-MS cleavage agent (name, accession) of one rule -- the residues `cut` at, in any
    order and case; the residue(s) `no_cut` before; cutting C-terminal of them or not -- or None
    when CLEAVAGE_AGENTS has no row for exactly that rule (the caller tags it)."""
    key = ("".join(sorted(set(str(cut or "").upper()))), "".join(sorted(set(
        str(no_cut or "").upper()))), bool(c_terminal))
    return CLEAVAGE_AGENTS.get(key)


def sage_cleavage(enz, source):
    """A present Sage enzyme block, read as Sage 0.14.7 reads it: database.rs From<EnzymeBuilder>
    takes an absent/null cleave_at as "KR", c_terminal as true, semi_enzymatic as false, and
    passes `restrict` through AS IS -- absent/null is no proline rule (DOCS.md's "default 'P'"
    is EnzymeBuilder::default, used only when the whole block is absent); enzyme.rs Enzyme::new
    makes cleave_at "" a non-specific digest (every peptide of min_len..max_len, missed
    cleavages forced to 0) and "$" no digestion (each FASTA entry as one peptide)."""
    cut = enz.get("cleave_at")
    cut = "KR" if cut is None else str(cut)
    if cut == "":
        return {"rule": "cleave_at=''", "name": NONSPECIFIC[0], "ac": NONSPECIFIC[1],
                "specificity": "none", "source": source}
    if cut == "$":
        return {"rule": "cleave_at='$'", "name": NO_CLEAVAGE[0], "ac": NO_CLEAVAGE[1],
                "specificity": "no_digestion", "source": source}
    restrict = enz.get("restrict") or ""
    c_term = True if enz.get("c_terminal") is None else bool(enz["c_terminal"])
    semi = bool(enz.get("semi_enzymatic"))
    nm = cleavage_agent(cut, restrict, c_term)
    return {"rule": f"cleave_at={cut!r}" + (f", restrict={restrict!r}" if restrict else "")
                    + ("" if c_term else ", c_terminal=false")
                    + (", semi_enzymatic=true" if semi else ""),
            "name": nm[0] if nm else None, "ac": nm[1] if nm else None,
            "specificity": "semi" if semi else "full", "source": source}


def fragger_cleavage(dg, source):
    """fragpipe_params.digest() -> the record's cleavage: non-specific when MSFragger searched
    non-enzymatic (num_enzyme_termini 0); otherwise the rule, named when it is one
    CLEAVAGE_AGENTS knows (one enzyme; its cut, no-cut and sense). The enzyme's NAME is neither
    read nor quoted: FragPipe's own way to search a rule its list lacks is to enter the cuts with
    the name `nonspecific` and the cleavage ENZYMATIC (TabMsfragger.java 24.0), and FragPipe
    renames the rule after its enzyme dropdown when the workflow has one -- H5's chymotrypsin-like
    rule ran as `argc` (its dropdown said `custom`, which FragPipe's list lacks: C0, 2026-11-04).
    So the cuts are the rule, and the record says only them (2.15, Ruling 90)."""
    enz, spec = dg["enzymes"], dg["specificity"]
    if spec == "none":
        return {"rule": "num_enzyme_termini 0", "name": NONSPECIFIC[0], "ac": NONSPECIFIC[1],
                "specificity": "none", "source": source}
    rule = "; ".join(f"cut {e['cut'] or '-'}, no cut {e['nocut'] or '-'}, sense "
                     f"{e['sense'] or '-'}" for e in enz) + f"; num_enzyme_termini {dg['termini']}"
    nm = None
    if len(enz) == 1 and str(enz[0]["sense"] or "").upper() in ("C", "N"):
        nm = cleavage_agent(enz[0]["cut"], enz[0]["nocut"], str(enz[0]["sense"]).upper() == "C")
    return {"rule": rule, "name": nm[0] if nm else None, "ac": nm[1] if nm else None,
            "specificity": spec, "source": source}


# ------------------------------------------------------------------- --enzyme (2.15) ----------
# FragPipe 24.0's enzyme list: FragPipe-GUI/src/main/java/org/nesvilab/fragpipe/tools/fragger/
# enzyme-defs-supported.txt @ tag 24.0 (sha256 c1c58a79..., read 2026-11-04; a copy is in
# tests/fixtures/fragpipe24_dda_preset/, and a test reads it against this table). name -> (the
# residues it cuts at, the residues it does not cut before ("" none), the side it cuts: C or N).
FRAGPIPE_ENZYMES = {
    "null": ("", "", "C"),
    "trypsin": ("KR", "P", "C"), "stricttrypsin": ("KR", "", "C"), "argc": ("R", "P", "C"),
    "aspn": ("D", "", "N"), "chymotrypsin": ("FLWY", "P", "C"),
    "cnbr": ("M", "", "C"), "elastase": ("AGILV", "P", "C"), "formicacid": ("D", "P", "C"),
    "gluc": ("DE", "P", "C"), "gluc_bicarb": ("E", "P", "C"),
    "lysc": ("K", "P", "C"), "lysc-p": ("K", "", "C"), "lysn": ("K", "", "N"),
    "lysn_promisc": ("KR", "", "N"),
    "nonspecific": ("@", "", "C"), "nocleavage": ("@", "@", "C"),
    "trypsin/chymotrypsin": ("FYWLKR", "P", "C"),
    "trypsin/cnbr": ("KR", "P", "C"), "trypsin_gluc": ("DEKR", "P", "C"),
    "trypsin_k": ("K", "P", "C"), "trypsin_r": ("R", "P", "C"),
    "thermolysin": ("ALIVFM", "DE", "N"),
}

# The --enzyme vocabulary -- fetch_fasta.py's words for the digestion enzyme(s) a lab added
# (fetch_fasta.ENZYME_FAMILIES: their autolysis entries stay contaminants) -- -> the row of
# FragPipe's list that is the word's in-silico rule. One word, one rule, for Sage and FragPipe
# alike; estimate_params.py and make_presets.py take the same words. This closes 2.14.2's known
# limit "estimate_params.py has no --enzyme, so a non-trypsin Sage search needs an --overrides
# block".
#   trypsin  FragPipe's `trypsin`: after K/R, not before P (PSI-MS Trypsin). Trypsin/P (no
#            proline rule, FragPipe's stricttrypsin) is spelt with the suffix: --enzyme trypsin/p.
#   gluc     FragPipe's `gluc`: after D and E (V8-DE, the phosphate-buffer specificity). A digest
#            in ammonium bicarbonate cuts after E only (FragPipe gluc_bicarb, V8-E): give that
#            rule with --overrides.
#   pepsin   None: FragPipe 24.0's list has no pepsin, and pepsin cleaves too broadly for one rule
#            (PSI-MS PepsinA's regexp (?<=[FL]) is its narrowest form); a pepsin digest is
#            searched non-specifically, and the record says why.
ENZYME_RULES = {"argc": "argc", "aspn": "aspn", "chymotrypsin": "chymotrypsin", "gluc": "gluc",
                "lysc": "lysc", "lysn": "lysn", "pepsin": None, "trypsin": "trypsin"}

# The name a FragPipe workflow gives a rule it states itself: FragPipe 24.0's own way (the
# TabMsfragger tooltip) -- the cuts entered under `nonspecific`, the cleavage ENZYMATIC
# (num_enzyme_termini 2). C0 ran it (2026-11-04, R1): fragger.params name nonspecific, the cuts as
# written, 2 termini, every PSM with 2 enzymatic termini.
FRAGPIPE_CUSTOM_NAME = "nonspecific"
# ...and the key that must NOT be beside it. Loading a workflow, FragPipe fires its enzyme
# dropdown's handler, which reloads the enzyme's name, no-cut and sense from its list entry; for
# `nonspecific` it makes the search non-specific (num_enzyme_termini 0), and a value the list lacks
# (`custom`) falls to its first sorted entry, argc (C0 fact 1, TabMsfragger.java 24.0 ~l.699-720).
# A rule the workflow states survives only without it.
FRAGPIPE_DROPDOWN_KEY = "msfragger.misc.fragger.enzyme-dropdown-1"
# The workflow keys that state enzyme 1's rule
FRAGPIPE_RULE_KEYS = ("msfragger.search_enzyme_name_1", "msfragger.search_enzyme_cut_1",
                      "msfragger.search_enzyme_nocut_1", "msfragger.search_enzyme_sense_1",
                      "msfragger.num_enzyme_termini")


class EnzymeError(ValueError):
    """An --enzyme value no rule can be written for: an unknown word, or enzymes that cannot be
    searched as one rule."""


def _vocabulary():
    """fetch_fasta, the one holder of the --enzyme words (imported here, not at the top: it
    imports estimate_params, which imports this module)."""
    import fetch_fasta
    return fetch_fasta


def _residues(text):
    return "".join(sorted(set(str(text or "").upper())))


def enzyme_rule(text):
    """--enzyme as given (fetch_fasta.py's words, comma-separated; `none` alone = no enzyme was
    added) -> the rule the search uses:
      {"given", "families" (sorted), "specificity": "full" | "none", "cut", "nocut", "sense",
       "c_terminal", "rule", "name", "ac", "fragpipe" (the FragPipe list entries), "why"}.
    Enzymes that cut on the same side with the same no-cut residue are one rule, the union of
    their cut residues (trypsin + Lys-C = trypsin's KR, not before P). A trailing /P (Trypsin/P,
    Lys-C/P), the engines' spelling of "no proline rule", drops the no-cut residue of the
    enzyme(s) it is written with (fetch_fasta reads the same word as the same enzyme). An enzyme
    with no rule (pepsin), or `none`, makes the search non-specific. Raises EnzymeError for an
    unknown word, `none` beside an enzyme, and enzymes that cut on different sides or with
    different no-cut rules (one rule cannot say both)."""
    ff = _vocabulary()
    words = f"{', '.join(ff.ENZYME_FAMILIES)} (comma-separated), or {ff.NO_ENZYME}"
    raw = str(text or "").strip()
    if not raw:
        raise EnzymeError(f"--enzyme was given empty: name the enzyme(s) added to the samples -- "
                          f"{words}")
    chunks = [c.strip() for c in raw.split(",") if c.strip()]
    if any(c.lower() == ff.NO_ENZYME for c in chunks):
        if len(chunks) > 1:
            raise EnzymeError(f"--enzyme {raw}: `{ff.NO_ENZYME}` says no enzyme was added, and "
                              f"cannot stand beside an enzyme. Give one or the other: {words}")
        return {"given": raw, "families": [], "specificity": "none", "cut": "", "nocut": "",
                "sense": None, "c_terminal": None, "rule": "no enzyme rule (non-specific)",
                "name": NONSPECIFIC[0], "ac": NONSPECIFIC[1], "fragpipe": [],
                "why": "--enzyme none: no digestion enzyme was added (a peptidome, a "
                       "hydrolysate or HLA peptides), so every peptide in the length range is "
                       "searched (non-specific)"}
    rows = {}                               # family -> (cut, nocut, sense, FragPipe entry)
    for chunk in chunks:
        no_p = bool(re.search(r"/\s*p\b", chunk, flags=re.I))
        try:
            fams = ff.parse_enzymes(chunk)
        except Exception as e:      # argparse.ArgumentTypeError: an unknown word, said with the list
            raise EnzymeError(f"--enzyme {raw}: {e}; or {ff.NO_ENZYME}") from None
        for fam in fams:
            entry = ENZYME_RULES.get(fam)
            if entry is None:
                rows[fam] = None
                continue
            cut, nocut, sense = FRAGPIPE_ENZYMES[entry]
            if no_p and nocut:
                entry = {"trypsin": "stricttrypsin", "lysc": "lysc-p"}.get(
                    entry, f"{entry} without its no-cut {nocut}")
                nocut = ""
            rows[fam] = (cut, nocut, sense, entry)
    families = sorted(rows)
    if any(r is None for r in rows.values()):
        ruleless = [f for f in families if rows[f] is None]
        return {"given": raw, "families": families, "specificity": "none", "cut": "",
                "nocut": "", "sense": None, "c_terminal": None,
                "rule": "no enzyme rule (non-specific)", "name": NONSPECIFIC[0],
                "ac": NONSPECIFIC[1], "fragpipe": [],
                "why": (f"--enzyme {raw}: {', '.join(ruleless)} has no in-silico rule (FragPipe "
                        f"24.0's enzyme list has none; pepsin cleaves too broadly for one), so "
                        f"the search is non-specific: every peptide in the length range")}
    groups = sorted({(r[2], r[1]) for r in rows.values()})
    if len(groups) > 1:
        said = "; ".join(f"{f}: cut {rows[f][0]}, no cut {rows[f][1] or '-'}, "
                         f"{'C' if rows[f][2] == 'C' else 'N'}-terminal" for f in families)
        raise EnzymeError(f"--enzyme {raw}: these enzymes cannot be searched as one rule ({said}): "
                          f"they cut on different sides or with different no-cut rules. Search "
                          f"the rule you want with --overrides, or name the enzymes that share "
                          f"one")
    sense, nocut = groups[0]
    cut = _residues("".join(r[0] for r in rows.values()))
    nm = cleavage_agent(cut, nocut, sense == "C")
    entries = sorted({r[3] for r in rows.values()})
    rule = f"cut {cut}, no cut {nocut or '-'}, sense {sense}"
    return {"given": raw, "families": families, "specificity": "full", "cut": cut,
            "nocut": nocut, "sense": sense, "c_terminal": sense == "C", "rule": rule,
            "name": nm[0] if nm else None, "ac": nm[1] if nm else None, "fragpipe": entries,
            "why": (f"--enzyme {raw}: {rule}"
                    + (f" -- the union of {', '.join(families)}'s rules" if len(families) > 1
                       else "")
                    + f" (FragPipe 24.0's enzyme list: {', '.join(entries)})")}


def semi_specific(rule):
    """`rule` (an enzymatic enzyme_rule) searched SEMI-specifically: one end of a peptide must
    follow the rule, the other may be anywhere (MSFragger num_enzyme_termini 1, Sage
    semi_enzymatic true). specificity "semi", which sage_cleavage / fpp.digest read back and the
    Methods name. A non-specific rule is already broader than that: EnzymeError."""
    if rule["specificity"] != "full":
        raise EnzymeError("--semi-specific makes an enzymatic digest semi-specific, and this "
                          f"search has no enzyme rule ({rule['why']}): it is already "
                          "non-specific, so give one or the other")
    return dict(rule, specificity="semi",
                rule=rule["rule"] + "; semi-specific",
                why=rule["why"] + "; semi-specific (--semi-specific): one end of each peptide "
                    "must follow the rule, the other may be anywhere")


def sage_enzyme_block(rule, base):
    """Sage's database.enzyme for `rule` (enzyme_rule), on `base` (estimate_params.DDA_DIGEST:
    the missed cleavages and lengths it keeps). A non-specific rule only empties cleave_at, the
    form H7 ran (Sage enzyme.rs: cleave_at "" is a non-specific digest); a rule sets cleave_at,
    restrict (Sage takes ONE residue) and c_terminal, leaving out what Sage's default already is
    -- so --enzyme trypsin writes exactly DDA_DIGEST."""
    blk = dict(base)
    if rule["specificity"] == "none":
        blk["cleave_at"] = ""
        return blk
    if len(rule["nocut"]) > 1:
        raise EnzymeError(f"--enzyme {rule['given']}: Sage's restrict is one residue, and this rule "
                          f"has {rule['nocut']}")
    blk["cleave_at"] = rule["cut"]
    if rule["specificity"] == "semi":
        blk["semi_enzymatic"] = True
    if rule["nocut"]:
        blk["restrict"] = rule["nocut"]
    else:
        blk.pop("restrict", None)
    if rule["c_terminal"]:
        blk.pop("c_terminal", None)
    else:
        blk["c_terminal"] = False
    return blk


def fragger_enzyme_keys(rule):
    """The FragPipe workflow keys that state `rule` (an enzymatic enzyme_rule) as enzyme 1:
    FRAGPIPE_CUSTOM_NAME + the cut, no-cut and side + fully enzymatic (one enzymatic terminus
    for a semi-specific rule). The caller deletes
    FRAGPIPE_DROPDOWN_KEY beside them (C0 fact 1)."""
    return {"msfragger.search_enzyme_name_1": FRAGPIPE_CUSTOM_NAME,
            "msfragger.search_enzyme_cut_1": rule["cut"],
            "msfragger.search_enzyme_nocut_1": rule["nocut"],
            "msfragger.search_enzyme_sense_1": rule["sense"],
            "msfragger.num_enzyme_termini": "1" if rule.get("specificity") == "semi" else "2"}


FRAGPIPE_NOT_AN_ENZYME = ("null", "nonspecific", "nocleavage")


def fragpipe_list_entry(cut, nocut, sense):
    """The name of the row of FragPipe 24.0's enzyme list (FRAGPIPE_ENZYMES) that is exactly the
    rule `cut`, `nocut`, `sense`, else None. The cuts and no-cuts compare as sets of residues, as
    FragPipe's own focus listener does (TabMsfragger.java 24.0 l.768: sorted chars), and the
    first row of the list wins, as it does there (matching.get(0)). A rule that is a list entry
    is written under that name and its dropdown, which FragPipe reloads to the same rule (l.707-
    714: name, cuts, no-cuts, sense from the entry, nothing else); any other rule is entered
    under `nonspecific` with no dropdown key (FRAGPIPE_CUSTOM_NAME)."""
    key = ("".join(sorted(set(str(cut or "").upper()))),
           "".join(sorted(set(str(nocut or "").upper()))), str(sense or "").upper())
    for name, (c, n, sd) in FRAGPIPE_ENZYMES.items():
        if name in FRAGPIPE_NOT_AN_ENZYME:
            continue
        if ("".join(sorted(set(c))), "".join(sorted(set(n))), sd) == key:
            return name
    return None


def fragpipe_dropdown_conflict(props):
    """A FragPipe workflow's properties (`msfragger.`-prefixed keys) -> what its enzyme dropdown
    would change in the rule its keys state when FragPipe loads it (C0 fact 1), else None.
    FragPipe applies the dropdown's list entry: `nonspecific` makes the search non-specific
    (num_enzyme_termini 0); a value its list lacks falls to its first sorted entry, argc (H5's
    `custom` ran as argc); any other entry renames the rule and reloads its no-cut residue and
    side (C0 saw the cut key survive). A dropdown naming the entry the keys already state (every
    shipped workflow) is no conflict; neither is a workflow without one.

    Ruling 104: this is for the presets the SKILL writes -- make_presets.py refuses its own preset
    when it finds one there (a custom rule must never sit under a dropdown name). A user-supplied
    or vendor workflow is NEVER refused on it: FragPipe runs it as it runs it, and the record and
    the Methods describe what ran from that run's fragger.params (fragger_cleavage). Elsewhere the
    message is information only. FragPipe's own SAVED workflow writes a dropdown back beside the
    rule it ran (`trypsin` beside nonspecific/FLMWY in HIVE dry runs), so replay the run's input
    workflow, not that file."""
    dd = str(props.get(FRAGPIPE_DROPDOWN_KEY) or "").strip()
    if not dd or dd.lower() == "null":
        return None
    name = str(props.get("msfragger.search_enzyme_name_1") or "").strip()
    cut = str(props.get("msfragger.search_enzyme_cut_1") or "").strip()
    nocut = str(props.get("msfragger.search_enzyme_nocut_1") or "").strip()
    sense = str(props.get("msfragger.search_enzyme_sense_1") or "").strip().upper()
    termini = str(props.get("msfragger.num_enzyme_termini") or "").strip()
    lead = (f"{FRAGPIPE_DROPDOWN_KEY}={dd} beside the rule {name or '-'} (cut {cut or '-'}, no "
            f"cut {nocut or '-'}, sense {sense or '-'}, num_enzyme_termini {termini or '-'}): ")
    tail = ("; FragPipe applies the dropdown when it loads the workflow (C0, 2026-11-04), and "
            "fragger.params records what it then searched.")
    defn, entry, why = FRAGPIPE_ENZYMES.get(dd), dd, ""
    if defn is None:
        defn, entry = FRAGPIPE_ENZYMES["argc"], "argc"
        why = f"{dd} is not in FragPipe 24.0's enzyme list, so FragPipe takes its first entry, argc; "
    if dd == "nonspecific":
        if termini == "0":
            return None
        return lead + "the nonspecific dropdown makes the search non-specific (num_enzyme_termini 0)" + tail
    changes = []
    if _residues(nocut) != _residues(defn[1]):
        changes.append(f"no cut {nocut or '-'} -> {defn[1] or '-'}")
    if sense != defn[2]:
        changes.append(f"sense {sense or '-'} -> {defn[2]}")
    if name != entry:
        changes.append(f"name {name or '-'} -> {entry}")
    if not changes:
        return None
    return lead + why + f"FragPipe reloads {entry} from its list: {'; '.join(changes)}" + tail


def diann_cut(cut, nocut="", c_terminal=True):
    """One rule in DIA-NN's --cut syntax, the form estimate_params.DIANN_DIGEST and
    run_search.sage_digest compare digests in: cutting after K/R, not before P = "K*,R*,!*P";
    cutting before D = "*D". Residues in the order given."""
    if c_terminal:
        return ",".join([f"{c}*" for c in cut] + [f"!*{n}" for n in nocut])
    return ",".join([f"*{c}" for c in cut] + [f"!{n}*" for n in nocut])


NONSPECIFIC_CUT = "(none: a non-specific digest)"     # run_search.sage_digest's words


def fragpipe_digest(path):
    """A FragPipe search's in-silico digest in estimate_params.DIANN_DIGEST's shape -- `cut`
    (diann_cut; NONSPECIFIC_CUT for num_enzyme_termini 0; "(semi-specific)" added for 1),
    `missed_cleavages` (enzyme 1's), `min_pep_len`, `max_pep_len`, `met_excision` (clip_nTerm_M)
    -- from a .workflow (its msfragger.* keys) or a fragger.params, the keys it states only; None
    when it cannot be read. For run_search.warn_contaminant_digest (review L7), which judges it
    against the digest fetch_fasta.py judged the contaminants on. It reads the keys as written: a
    workflow whose dropdown would overwrite them (fragpipe_dropdown_conflict) is searched
    otherwise, and fragger.params, written by the search, is then the record."""
    import fragpipe_params as fpp
    try:
        if str(path).endswith(fpp.FRAGGER_PARAMS):           # fragger.params, or a copy's name
            kv = fpp.read_fragger_params(path)
        else:
            props, _ = fpp.read_workflow(path)
            kv = {k[len("msfragger."):]: v for k, v in props.items()
                  if k.startswith("msfragger.")}
    except (OSError, TypeError, ValueError):
        return None
    dg = fpp.digest(kv)
    got = {}
    if dg["specificity"] == "none":
        got["cut"] = NONSPECIFIC_CUT
    elif dg["enzymes"]:
        cuts = [diann_cut("" if e["cut"] == "@" else e["cut"], e["nocut"],
                          str(e["sense"] or "C").upper() != "N") for e in dg["enzymes"]]
        got["cut"] = ",".join(c for c in cuts if c) + (
            " (semi-specific)" if dg["specificity"] == "semi" else "")
        if dg["enzymes"][0]["missed"] is not None:
            got["missed_cleavages"] = int(dg["enzymes"][0]["missed"])
    for key, name in (("digest_min_length", "min_pep_len"), ("digest_max_length", "max_pep_len")):
        try:
            got[name] = int(float(kv[key]))
        except (KeyError, TypeError, ValueError):
            pass
    if "clip_nTerm_M" in kv:
        got["met_excision"] = str(kv["clip_nTerm_M"]).strip().lower() in ("1", "true")
    return got
