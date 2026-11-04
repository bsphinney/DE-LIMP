#!/usr/bin/env python3
"""
cleavage_rules.py -- THE reading of a search's in-silico cleavage, from the engine's own record
(2.14.2, L11 N1/N4/N5). make_methods.py writes the Methods and the SDRF cleavage agent from it;
run_search.py sizes a non-specific search's job from it. Stdlib only.

  * sage_cleavage(enzyme block)        -- Sage's database.enzyme, as Sage 0.14.7 reads it
  * fragger_cleavage(fragpipe_params.digest(...)) -- MSFragger's fragger.params / workflow keys
  * cleavage_agent(cut, no_cut, c_terminal) -- the PSI-MS term of one rule, or None

Each of the first two returns {"rule", "name", "ac", "specificity", "source"}: `specificity` is
"none" (a non-specific search), "no_digestion" (Sage "$"), "semi" or "full"; `name`/`ac` the
PSI-MS cleavage term when one fits, else None (the caller tags the rule).
"""

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
    CLEAVAGE_AGENTS knows (one enzyme; its cut, no-cut and sense). The enzyme's NAME is not read:
    FragPipe's own way to search a rule its list lacks is to enter the cuts with the name
    `nonspecific` and the cleavage ENZYMATIC (TabMsfragger.java 24.0), so the cuts are the rule."""
    enz, spec = dg["enzymes"], dg["specificity"]
    if spec == "none":
        return {"rule": f"num_enzyme_termini 0 (cleavage {enz[0]['name'] if enz else 'none'})",
                "name": NONSPECIFIC[0], "ac": NONSPECIFIC[1], "specificity": "none",
                "source": source}
    rule = "; ".join(f"{e['name']}: cut {e['cut'] or '-'}, no cut {e['nocut'] or '-'}, sense "
                     f"{e['sense'] or '-'}" for e in enz) + f"; num_enzyme_termini {dg['termini']}"
    nm = None
    if len(enz) == 1 and str(enz[0]["sense"] or "").upper() in ("C", "N"):
        nm = cleavage_agent(enz[0]["cut"], enz[0]["nocut"], str(enz[0]["sense"]).upper() == "C")
    return {"rule": rule, "name": nm[0] if nm else None, "ac": nm[1] if nm else None,
            "specificity": spec, "source": source}
