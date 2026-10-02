#!/usr/bin/env python3
"""
protein_ids.py  --  The ONE reading of a protein identifier as the accession two searches are
compared on. run_search.py's FragPipe adapter writes it; compare_searches.py compares on it.

FragPipe's combined_protein.tsv names each protein two ways: `Protein` (the FASTA header
field, `sp|P12345|ALBU_HUMAN`) and `Protein ID` (`P12345`). The DDA adapter took `Protein`, so
its Protein.Group never equalled DIA-NN's (`P12345`), and compare_searches.py reported 0
proteins shared between FragPipe and DIA-NN on the same raw files.

    DIA-NN Protein.Group            P12345;Q67890                  -> P12345
    FragPipe Protein / Protein ID   sp|P12345|ALBU_HUMAN / P12345  -> P12345
    FragPipe's own contaminant      contam_sp|O43790|KRT86_HUMAN   -> contam_O43790
                                    (its Protein ID, O43790, drops the tag)
    the skill's contaminant FASTA   sp|Cont_P00761|TRYP_PIG        -> Cont_P00761
    a FragPipe group (fragpipe_group: Protein ID + Indistinguishable Proteins)
                                    Q15582 + "sp|Cont_P55906|BGH3_BOVIN"  -> Q15582;Cont_P55906
    Sage proteins                   sp|P12345|ALBU_HUMAN;sp|...    -> P12345
    an isoform                      P12345-2                       -> P12345

compare_analyses.R's normalize_protein_id() (a port of DE-LIMP's Run Comparator) does the same
for finished DE tables; tests/test_fragpipe_protein_ids.py checks the two agree on DIA-NN's and
FragPipe's identifiers. Where they part, the R one is the narrower: out of a FASTA header it
reads only 6-character accessions, from a group of headers the LAST one, and it strips an
isoform suffix only when it ends the whole group ("P04637-2;Q1" stays P04637-2). This one reads
any accession, the first member, and its isoform suffix.

A whole group, written the way DIA-NN writes Protein.Group, is group_accessions(): every member
read by header_accession(), isoform suffixes KEPT (the accession as the database has it). Sage's
`proteins` is a ';' list of full FASTA IDs, and the contaminant filter (contaminants.R
contaminant_regex: an accession STARTS with Cont_/contam_) matched none of them -- a Sage DE
tested 171 Cont_ groups on gabrig's HeL50 UnvPe (2026-09-29) while its methods said "none":

    Sage proteins   sp|P12345|ALBU_HUMAN;sp|Cont_P02769|ALBU_BOVIN  -> P12345;Cont_P02769
    its entry names                                                  -> ALBU_HUMAN;ALBU_BOVIN
"""
import re

_DB = re.compile(r"^(?P<tag>.*?)(?:sp|tr)$")      # "sp", "tr", "contam_sp", "rev_tr", ...
_ISOFORM = re.compile(r"-\d+$")


def header_accession(token):
    """The accession of one identifier: the middle field of a `db|ACCESSION|NAME` FASTA header,
    keeping a tag written in front of the database field (`contam_sp|...` -> `contam_...`);
    anything else as given, trimmed."""
    t = str(token or "").strip()
    parts = t.split("|")
    if len(parts) >= 3 and parts[1].strip():
        m = _DB.match(parts[0].strip())
        return (m.group("tag") if m else "") + parts[1].strip()
    return t


def fragpipe_protein_id(protein, protein_id):
    """FragPipe's protein as the accession DIA-NN would report: `Protein ID`, plus the tag its
    `Protein` carries in front of the database field, which `Protein ID` drops -- FragPipe's
    `contam_sp|P00167|CYB5_HUMAN` would otherwise read as the sample's own P00167. Without a
    `Protein ID`, the accession read from `Protein`."""
    pid = str(protein_id or "").strip()
    if not pid:
        return header_accession(protein)
    m = _DB.match(str(protein or "").strip().split("|", 1)[0])
    tag = m.group("tag") if (m and "|" in str(protein or "")) else ""
    return pid if not tag or pid.startswith(tag) else tag + pid


# How FragPipe's combined_protein.tsv separates the members of `Indistinguishable Proteins`:
# ", " (FragPipe 24.0 on HIVE, 2026-09-30: "sp|Q9HAP6|LIN7B_HUMAN, sp|Q9NUP9|LIN7C_HUMAN").
FRAGPIPE_INDISTINGUISHABLE_SEP = ","


def fragpipe_group(protein, protein_id, indistinguishable=""):
    """A FragPipe protein group written as DIA-NN writes Protein.Group: the leading protein
    (fragpipe_protein_id), then every protein FragPipe lists under `Indistinguishable Proteins`
    -- proteins with the same peptides, which DIA-NN would report in the same group -- each read
    by header_accession() (keeping a `contam_` tag), each once. The DDA adapter wrote the leading
    protein alone, so a Cont_/contam_ entry listed only there was never seen by the contaminant
    rule (any accession of the group carries the tag), which the DIA-NN path applies to every
    member. '' when there is no leading protein."""
    lead = fragpipe_protein_id(protein, protein_id)
    if not lead:
        return ""
    rest = [header_accession(t) for t in
            str(indistinguishable or "").split(FRAGPIPE_INDISTINGUISHABLE_SEP)]
    return ";".join(dict.fromkeys(a for a in [lead] + rest if a))


def normalize_protein_id(group):
    """The accession a protein group is compared on: its first member, read by
    header_accession(), without an isoform suffix. '' for a blank group."""
    first = str(group or "").split(";", 1)[0]
    return _ISOFORM.sub("", header_accession(first)).strip()


def group_accessions(group):
    """A ';' list of identifiers as DIA-NN's Protein.Group writes one: each member read by
    header_accession() (a bare accession stays as it is), in order, each once."""
    out = [header_accession(t) for t in str(group or "").split(";")]
    return ";".join(dict.fromkeys(a for a in out if a))


def header_entry_name(token):
    """The entry name of one identifier: the third field of a `db|ACCESSION|NAME` FASTA header
    (`sp|P12345|ALBU_HUMAN` -> ALBU_HUMAN), '' for anything else."""
    parts = str(token or "").strip().split("|")
    return parts[2].strip() if len(parts) >= 3 else ""


def group_entry_names(group):
    """The entry names of a ';' list, as DIA-NN's Protein.Names writes them; '' when no member
    carries one."""
    names = [header_entry_name(t) for t in str(group or "").split(";")]
    return ";".join(dict.fromkeys(n for n in names if n))
