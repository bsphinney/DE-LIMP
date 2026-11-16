"""A small REAL-looking reference database for the FRAN tests: UniProt-style headers, one real
selenoprotein (U), and the fetch_fasta.py sidecar -- so fran_deposit's database guard passes it.
Made in a temp dir only; no client data."""
import json
import os

REAL_FASTA = (">sp|P07203|GPX1_HUMAN Glutathione peroxidase 1 GN=GPX1 PE=1 SV=2\n"
              "MCAAQRSAAALAAAAPRTVYAFSARPLAGGEPVSLGSLRGKVLLIENVASLUGTTVRDYTQMNELQRRLGPRGLVVLGFPCNQFGHQEN\n"
              ">sp|P04637|P53_HUMAN Cellular tumor antigen p53 GN=TP53 PE=1 SV=4\n"
              "MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP\n")


def write_real_db(directory, name="real_human.fasta"):
    """-> the FASTA path; `<path>.meta.json` is fetch_fasta.py-shaped."""
    p = os.path.join(directory, name)
    with open(p, "w") as fh:
        fh.write(REAL_FASTA)
    with open(p + ".meta.json", "w") as fh:
        # what fetch_fasta.py writes when it was not told the organism: the fixture must not
        # introduce one (the FRAN manifest tests assert an unknown organism stays absent)
        json.dump({"organism": "", "taxid": 0, "organism_source": "none"}, fh)
    return p
