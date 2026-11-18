# Where the protein card's facts come from

The protein page ("show me the data for <protein>") looks every protein and gene fact up at build time in the sources below, through `scripts/protein_facts.py`. It never writes a fact from memory. Each value on the card is the source's own, and `tests/test_protein_facts.py` checks that every value appears in a response recorded from the source it is shown under. A source that cannot be reached, has no record, or does not cover the species says so on the card in one line.

Verified by live calls on 2026-10-28. Responses recorded for the tests: `tests/fixtures/protein_facts/` (re-record with `record.py` there; it keeps UniProt entries without their literature `references`, which the card never reads and whose author lists the repository's identifier guard would have to read).

Identifiers come from UniProt's cross-references (GeneID, HGNC, Ensembl, HPA), never from a search on a gene name. The Human Protein Atlas search matches loosely ("ALB" also returns CALB1 and LALBA).

Every request: `User-Agent: ucdavis-proteomics-skill (+https://github.com/bsphinney/UCD-Proteomics-AI-Skill)`, 20 s timeout, https only, hosts limited to the five below (a redirect elsewhere is refused). Offline (`--offline`) asks nothing; every source shows "Not reached: --offline".

| Source | Endpoint | Release or version | Licence | Rate limit | Species |
|---|---|---|---|---|---|
| UniProt | `GET https://rest.uniprot.org/uniprotkb/<ACC>.json` | response header `x-uniprot-release` (not in the body) | CC BY 4.0 | none stated; none hit | all |
| NCBI Gene | `esummary.fcgi?db=gene&id=<GeneID>&retmode=json` | none; the summary text carries "[provided by RefSeq, <date>]"; the card gives the fetch date | no restrictions (keep NCBI's disclaimer, acknowledge NLM) | 3 requests/s per IP without a key (10/s with one); calls here are 0.34 s apart | all |
| NCBI Protein | `esummary.fcgi?db=protein&id=<ACC>`, then `efetch.fcgi?db=protein&id=<ACC>&rettype=gp&retmode=text` for the gene (`/gene=`, `/db_xref="GeneID:"`) | none | as NCBI Gene | as NCBI Gene | all |
| HGNC | `GET https://rest.genenames.org/fetch/hgnc_id/<HGNC:n>` with `Accept: application/json` (default is XML) | `GET /info` -> `lastModified` (HGNC updates continuously) | CC0 | 10 requests/s; overload can bring an IP block | human only |
| Ensembl | `GET https://rest.ensembl.org/lookup/id/<id>?content-type=application/json` (unversioned id) | `GET /info/data` -> `releases[0]` | no statement found on the REST pages | 55,000 requests/hour; 429 with `Retry-After` | vertebrates and some model organisms, including S. cerevisiae (the old `rest.ensemblgenomes.org` has a broken TLS certificate and is not used) |
| Human Protein Atlas | `GET https://www.proteinatlas.org/<ENSG>.json` | none: HPA states no version or date in its answer, so the card says "version not stated by the source" and never types one in | CC BY 4.0 (cite Uhlén M et al. (2015) Tissue-based map of the human proteome. Science 347(6220):1260419, doi:10.1126/science.1260419 (PubMed 25613900, checked 2026-11-17); give the item URL) | none documented | human genes only |

## What each block shows

- **UniProt.** The card's UniProt columns (name, function, location, coverage) are from the entry. The release is the header. An accession UniProt has merged is answered with a redirect to the surviving entry; the card uses that entry and says so. A demerged accession (for example P04277 -> P02769 and P02768) lists every target and picks none. The UniParc id is given when UniProt gives one.
- **NCBI Gene.** Symbol, description, aliases, organism, chromosome, map location and the RefSeq summary, verbatim.
- **HGNC.** Symbol, approved name, HGNC id, status, locus group and type, location, aliases and previous symbols and names, gene groups. Optional fields that HGNC omits are omitted.
- **Ensembl.** Gene id, name, description, biotype, species, assembly, region, strand, source.
- **Human Protein Atlas.** `Blood concentration - Conc. blood IM [pg/L]` and `... MS [pg/L]`, `Secretome location`, `Evidence`, `RNA blood cell specificity`. Units are pg/L as HPA states them. Thousands separators are added for display only. Types follow HPA (a number in one record, a string in another). A null is "no HPA value", not zero and not "absent from blood".
- **An accession UniProt does not have** (a RefSeq id, for example NP_006308.3): NCBI Protein (title, organism, length, update date, and the gene found in the GenBank record), then NCBI Gene for that GeneID.

## Not used

- **neXtProt is retired.** `api.nextprot.org` and `www.nextprot.org` redirect to `expasy.org/archives/nextprot`, which says neXtProt "has reached its end of life and no longer provides data or services"; its last release was September 2023. UniProt entries no longer list a neXtProt cross-reference. It is not fetched, and the card makes no statement about it. (HPA still prints a "NeXtProt evidence" field derived from it; the card does not show that field.)
- **Human Plasma PeptideAtlas** has no documented public API (its pages ask robots to use bulk downloads), so it is not queried. HPA is the plasma source.
