#!/usr/bin/env python3
"""denovo_lca.py -- per-peptide LCA species attribution from a DIAMOND nr search.

Ported from DE-LIMP scripts/lca_attribute.py. This is the step that actually
decides the species; the FDR model does not.

WHY LCA AND NOT BEST-HIT. A de novo peptide can be perfectly sequenced and
perfectly homologous and still have its single best hit land on the wrong
organism, because conserved peptides match many taxa. Measured on the HeLa
entrapment (true answer: human), by best hit: only 29.5% primate, 68.5% mammal,
13.1% not even metazoan. Rescoring does not fix this -- the best 5-feature FDR
model moved primate attribution to just 32.9% and left mammal flat at 68.5%,
while discarding 715 genuinely-correct peptides. Species mis-assignment is an
LCA problem, so run this over the UNFILTERED hits and use the FDR score as a
reported confidence tier.

REQUIRES MULTIPLE HITS PER PEPTIDE. "LCA over all hits" is meaningless if the
search kept one hit per peptide. The ocelot production search used
--max-target-seqs 1 (0.9 hits/peptide), which silently degrades this step to
best-hit -- exactly the failure it exists to prevent. Use --max-target-seqs 25
(the HeLa benchmark setting, 22.5 hits/peptide); denovo_search_cmd.py emits it.
This script warns when the input averages < 2 hits per peptide.

Reads NCBI taxonomy from the PRISTINE nodes.dmp.preDmnd (real ranks -- the
diamond DB itself was built with ranks neutralised to "no rank", which affects
diamond's labels but not the tree).

Usage: denovo_lca.py <blast_results.tsv> <out_prefix>
  hits tsv: qseqid sseqid pident len mm gap qs qe ss se evalue bitscore staxids
"""
import sys, collections

NODES = "/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/ncbi_nr/nodes.dmp.preDmnd"
NAMES = "/quobyte/proteomics-grp/bioinformatics_programs/blast_dbs/ncbi_nr/names.dmp"
HITS, OUTP = sys.argv[1], sys.argv[2]

# anchor taxa for kingdom/domain classification
METAZOA, BACTERIA, ARCHAEA, VIRUSES, FUNGI, VIRIDIPLANTAE = 33208, 2, 2157, 10239, 4751, 33090
DIAG_RANKS = {"species", "subspecies", "species group", "species subgroup",
              "genus", "subgenus", "strain", "isolate", "varietas", "forma"}

def _warn_if_single_hit(path):
    """One hit per peptide means this script cannot do what it says it does."""
    q, n = set(), 0
    with open(path) as fh:
        for i, line in enumerate(fh):
            if i >= 200000: break
            q.add(line.split("\t", 1)[0]); n += 1
    if q and n / len(q) < 2.0:
        sys.stderr.write(
            f"[lca] WARNING: {n/len(q):.1f} hits per peptide. LCA over all hits "
            f"needs several; with one hit each this degrades to best-hit, which "
            f"is the 36%-mis-assignment error mode LCA exists to fix. Re-run "
            f"DIAMOND with --max-target-seqs 25.\n")

_warn_if_single_hit(HITS)
sys.stderr.write("[lca] loading taxonomy...\n")
parent, rank = {}, {}
for line in open(NODES):
    p = line.split("\t|\t")
    t = int(p[0]); parent[t] = int(p[1]); rank[t] = p[2]
name = {}
for line in open(NAMES):
    p = line.split("\t|\t")
    if len(p) >= 4 and p[3].rstrip("\t|\n") == "scientific name":
        name[int(p[0])] = p[1]
sys.stderr.write(f"[lca] {len(parent)} nodes, {len(name)} names\n")

def lineage(t):
    out, seen = [], set()
    while t and t not in seen:
        seen.add(t); out.append(t)
        pt = parent.get(t)
        if pt is None or pt == t:
            break
        t = pt
    return out  # taxon -> ... -> root

def lca(taxids):
    taxids = [t for t in set(taxids) if t in parent]
    if not taxids:
        return None
    lins = [lineage(t) for t in taxids]
    common = set.intersection(*[set(l) for l in lins])
    for t in lins[0]:          # deepest (lowest) shared node
        if t in common:
            return t
    return 1

# group hits by peptide
hits = collections.defaultdict(list)   # pep -> [(bitscore, [taxids], pident)]
for line in open(HITS):
    f = line.rstrip("\n").split("\t")
    if len(f) < 13:
        continue
    tids = [int(x) for x in f[12].replace(";", " ").split() if x.isdigit() and int(x) > 0]
    hits[f[0]].append((float(f[11]), tids, float(f[2])))

cat = collections.Counter(); host_sp = collections.Counter(); micro = collections.Counter()
n_diag = 0
with open(OUTP + "_peptide_lca.tsv", "w") as o:
    o.write("peptide\tn_hits\ttop_pident\tlca_taxid\tlca_name\tlca_rank\tcategory\tdiagnostic\n")
    for pep, hl in hits.items():
        topbs = max(h[0] for h in hl)
        keep = [h for h in hl if h[0] >= 0.9 * topbs]   # MEGAN-style top-10% window
        l = lca([t for h in keep for t in h[1]])
        if l is None:
            continue
        lin = set(lineage(l))
        c = ("host" if METAZOA in lin else
             "microbiome" if (BACTERIA in lin or ARCHAEA in lin or VIRUSES in lin) else
             "plant/fungal" if (FUNGI in lin or VIRIDIPLANTAE in lin) else "other/conserved")
        rk = rank.get(l, "no rank"); nm = name.get(l, str(l))
        diag = rk in DIAG_RANKS
        n_diag += diag
        o.write(f"{pep}\t{len(hl)}\t{max(h[2] for h in keep):.1f}\t{l}\t{nm}\t{rk}\t{c}\t{int(diag)}\n")
        cat[c] += 1
        if c == "host" and diag:
            host_sp[f"{nm} ({rk})"] += 1
        if c == "microbiome":
            micro[nm] += 1

print(f"=== peptides with LCA: {sum(cat.values())}  (diagnostic species/genus: {n_diag}) ===")
print("--- category ---")
for c, n in cat.most_common():
    print(f"  {n:>6}  {c}")
print("--- top HOST species (diagnostic, species/genus LCA) ---")
for s, n in host_sp.most_common(20):
    print(f"  {n:>6}  {s}")
print("--- top microbiome taxa ---")
for s, n in micro.most_common(12):
    print(f"  {n:>6}  {s}")
