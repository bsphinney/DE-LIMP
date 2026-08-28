#!/usr/bin/env python3
"""
resolve_organism.py -- turn what the user SAYS into the proteome you actually search.

Organism cannot be detected from a raw file, so the skill asks for it. But asking is
only half the job: the search needs a UniProt *proteome accession* (UP000000589), not
the word "mouse". Searching mouse data against a human FASTA does not error; it just
quietly loses most of the proteome, which is the worst possible failure mode.

THIS IS A THIN WRAPPER. The resolver lives in `fetch_fasta.py` (`resolve`), which is
also what builds the FASTA -- one concept, one implementation. Prefer calling it:

    python3 fetch_fasta.py resolve --organism "mouse"
    python3 fetch_fasta.py fetch --proteome <confirmed UPID> --out search.fasta

This wrapper exists so the older `resolve_organism.py --organism ...` call keeps working
and keeps emitting the same flat JSON shape (`uniprot_proteome`, `taxid`, `organism`).

The previous standalone implementation ranked proteomes with `"reference" in
proteomeType`, which is TRUE for "Non Reference proteome" -- so `--organism
"Saccharomyces cerevisiae"` returned UP000243666, the *S. cerevisiae killer virus M1*
proteome. Delegating removes that class of bug rather than fixing it twice.

Accepts a common name, a Latin name, a taxid, or a UP accession.

Usage
  python3 resolve_organism.py --organism "mouse"
  python3 resolve_organism.py --organism 10090
  python3 resolve_organism.py --organism UP000000589
  python3 resolve_organism.py --list
"""
import argparse, json, sys, os, io, contextlib

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fetch_fasta as ff


def resolve(text):
    """-> the flat dict this script has always returned."""
    ns = argparse.Namespace(organism=str(text), taxid=None, size=25, list=False)
    # cmd_resolve prints the rich payload; capture it rather than duplicating the logic.
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            rc = ff.cmd_resolve(ns)
    except SystemExit as e:
        return {"error": f"could not resolve {text!r}: {e}"}
    if rc != 0:
        return {"error": f"no UniProt proteome match for {text!r}"}
    d = json.loads(buf.getvalue())
    sel = d.get("selected")
    if not sel:
        return {"error": f"no UniProt proteome match for {text!r}"}
    out = {
        "taxid": sel.get("taxid"),
        "organism": sel.get("organism"),
        "uniprot_proteome": sel.get("proteome_id"),
        "proteome_type": sel.get("proteome_type"),
        "protein_count": sel.get("protein_count"),
        "source": "uniprot",
        "n_candidates": d.get("n_candidates"),
        "alternatives": [c["proteome_id"] for c in d.get("candidates", [])[1:4]],
        "needs_menu": d.get("needs_menu"),
        "notes": d.get("notes") or [],
    }
    if d.get("needs_menu"):
        out["ambiguous"] = ("More than one plausible proteome matched — show the user "
                            "`candidates` from `fetch_fasta.py resolve` and let them pick "
                            "instead of accepting this one.")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--organism", help="common name, Latin name, taxid, or UPxxxxxxxxx")
    ap.add_argument("--list", action="store_true", help="show the curated organism table")
    ap.add_argument("--offline", action="store_true",
                    help="(accepted for back-compat; the resolver falls back to the "
                         "curated table automatically when UniProt is unreachable)")
    a = ap.parse_args()

    if a.list:
        print(json.dumps([{"taxid": tx, "organism": n, "offline_fallback_proteome": up,
                           "aliases": al}
                          for tx, (n, up, al) in sorted(ff.ORGANISM_TAXIDS.items())], indent=2))
        return
    if not a.organism:
        sys.exit("Need --organism (or --list). The organism cannot be detected from raw "
                 "files — ASK the user, then resolve it here.")

    out = resolve(a.organism)
    out["query"] = a.organism
    if "error" not in out:
        out["confirm_with_user"] = (
            f"Searching {out.get('organism') or 'this organism'} "
            f"(taxid {out.get('taxid')}) against UniProt {out['uniprot_proteome']}. "
            "Confirm before the search starts — a wrong proteome does not error, it "
            "silently loses most of the identifications.")
    print(json.dumps(out, indent=2))
    if "error" in out:
        sys.exit(2)


if __name__ == "__main__":
    main()
