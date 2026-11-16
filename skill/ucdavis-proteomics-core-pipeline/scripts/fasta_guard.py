"""
fasta_guard.py -- the ONE definition of "this database holds sequences that do not exist", and of
the `<fasta>.nofran.json` marker every script that WRITES such a database leaves beside it.

Why it exists: the Core's peptide corpus (FRAN) must never learn a shuffled, recoded or variant
sequence. fran_deposit.py `check` calls scan_fasta() on every FASTA a search used; the writers of
test databases call mark_nofran(); FRAN's ingest mirrors PATTERNS below (a test in
tests/test_fasta_guard.py ties the two copies together).

Three verdicts (classify()):
  test_database      HARD refusal, recorded in the receipt. Any of:
                       * a `<fasta>.nofran.json` sidecar with corpus_eligible false;
                       * a header matching ENTRAP_ (any case), or an entry id containing decoy /
                         shuffle (a description may say "Decoy": Arabidopsis has such a protein),
                         or a rev_ entry whose target is not in the file (a complete target +
                         rev_ pair set is FragPipe's ordinary database and passes);
                       * a header carrying a variant tag (_VAR_ / _CON_ / _ENC_ / _ENT_<n>,
                         upper case only);
                       * J anywhere in a sequence (Xle: a real database never writes it).
  needs_agent_check  NOT a refusal, nothing recorded until a person decides
                     (`stage --database-ok "<why>"`, which is recorded). Any of:
                       * U in an entry whose header does not name a selenoprotein family;
                       * O in an entry whose header does not name pyrrolysine / a pyrrolysine
                         enzyme;
                       * B or Z in more than BZ_TOLERANCE (0.2%) of entries;
                       * sequences added with fetch_fasta --add-fasta (named in the sidecar), which
                         are scanned with the same rules (a hit there is a hard refusal);
                       * no fetch_fasta.py sidecar (`<fasta>.meta.json`): a hand-built database.
  ok
U (selenocysteine), O (pyrrolysine) and B/Z/X (ambiguity codes) are REAL in UniProt, which is why
only J is a hard refusal. Stdlib only.
"""
import datetime
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import safe_write  # noqa: E402

SIDECAR_SUFFIX = ".nofran.json"
KINDS = ("entrapment", "recoded", "synthetic_variant", "denovo_candidate", "decoy", "other")

# --- the header patterns (case-insensitive). FRAN's corpus_ingest.py mirrors these two strings.
# Measured on the 17 non-empty reference FASTAs of the Core's cache (2026-11-09): "decoy" in a
# DESCRIPTION is real biology (Arabidopsis has a gene named DECOY and a protein named "Decoy"), so the
# two English words are tested on the entry's IDENTIFIER (the first token) only; the two
# distinctive tags are tested on the whole header.
HEADER_TEST = r"ENTRAP_"
ID_TEST = r"decoy|shuffle"
# FragPipe's target+decoy database (Philosopher) is a REAL database: every target `<id>` has a
# `rev_<id>` entry. So rev_ is test material only when it is not part of such a pair set: a decoy
# whose target is absent (a decoys-only database, or rev_rev_). Per-file, in scan_fasta().
DECOY_TAG = "rev_"
# CASE-SENSITIVE: the test databases write upper-case `P12345_CON_7QE` / `X_ENT_3`, and a real custom
# proteome (a kelp genome annotation, measured in the Lance corpus) names contigs
# `...Genewise1Plus.C_con_440104`; `scaffold_ent_12` is the same hazard for _ENT_<n>.
HEADER_VARIANT = r"_VAR_|_CON_|_ENC_|_ENT_\d"
HEADER_TEST_RE = re.compile(HEADER_TEST, re.I)
ID_TEST_RE = re.compile(ID_TEST, re.I)
HEADER_VARIANT_RE = re.compile(HEADER_VARIANT)      # case-sensitive: see HEADER_VARIANT
# Accession-level form for FRAN's protein-group column (no leading '>').
# Two halves, because the variant tags are case-sensitive (see HEADER_VARIANT) and SQL regex has no
# scoped flags. FRAN's corpus_ingest.py pins both literals (SYNTHETIC_ACCESSION_CI / _CS).
ACCESSION_TEST_CI = r"ENTRAP_|shuffle"
ACCESSION_TEST_CS = r"_VAR_|_CON_|_ENC_|_ENT_\d"
ACCESSION_TEST = ACCESSION_TEST_CI + "|" + ACCESSION_TEST_CS

# --- residues. J is never real. U/O are real only in named families.
HARD_RESIDUES = "J"
AMBIGUOUS_RESIDUES = "BZ"
# B/Z (and X) are ambiguity codes old sequencing left in a few real UniProt entries. Measured on the 17
# non-empty reference FASTAs: at most 6 entries in 54,825 (mouse, 0.011%). A database where more than
# this FRACTION of entries carry B or Z is recoded or machine-made, so it asks; one stray entry in a
# real proteome does not.
BZ_TOLERANCE = 0.002
# The eukaryotic selenoprotein families by name (UniProt/RefSeq names are not uniform: SELENOI is
# "Ethanolaminephosphotransferase 1", SELENOO "protein adenylyltransferase SelO", SEPHS2 "selenide,
# water dikinase 2"), plus the bacterial Sec enzymes. Extended from the measurement on real
# proteomes; a miss is a needs_agent_check, never a refusal.
SELENO_HEADER = re.compile(
    r"seleno|selenium|selenide|\bSel[A-Z]\b|glutathione peroxidase|thioredoxin reductase|"
    r"iodothyronine|deiodinase|methionine-R-sulfoxide reductase B1|ethanolaminephosphotransferase|"
    r"formate dehydrogenase|glycine reductase|hydrogenase|"
    r"GN=(?:SEL[A-Z0-9]+|SEP[A-Z0-9]+|GPX\d|GPX[A-Z]*\d|TXNRD\d|DIO\d|MSRB1|EPT1)\b", re.I)
PYRROLYSINE_HEADER = re.compile(
    r"pyrrolysine|methylamine methyltransferase|methyltransferase corrinoid|"
    r"(?:mono|di|tri)methylamine", re.I)


def sidecar_path(fasta):
    return fasta + SIDECAR_SUFFIX


def mark_nofran(fasta, kind, made_by, why):
    """Write `<fasta>.nofran.json` (corpus_eligible false). EVERY script that writes a non-reference
    FASTA calls this right after writing the FASTA. Never `<fasta>.meta.json`: that name is
    fetch_fasta.py's, and fran_deposit.py reads it for organism and database identity. Raises on a
    bad `kind` or an unwritable sidecar -- a test database that cannot be marked must say so."""
    if kind not in KINDS:
        raise ValueError(f"kind must be one of {KINDS}, not {kind!r}")
    rec = {"corpus_eligible": False, "kind": kind, "made_by": made_by,
           "made_at": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
           "why": why}
    safe_write.write_file(sidecar_path(fasta), json.dumps(rec, indent=2) + "\n")
    return sidecar_path(fasta)


def read_nofran(fasta):
    """The sidecar dict when it says corpus_eligible false, else None. A sidecar that exists but
    cannot be parsed counts as 'says false' (returned with unreadable=True): an unreadable marker
    must not read as permission."""
    p = sidecar_path(fasta)
    if not os.path.isfile(p):
        return None
    d = safe_write.load_json(p)
    if not isinstance(d, dict):
        return {"corpus_eligible": False, "unreadable": True}
    return d if d.get("corpus_eligible") is False else None


def is_fetch_fasta_sidecar(fasta):
    """True when a fetch_fasta.py sidecar describes this database (see _fetch_meta)."""
    return _fetch_meta(fasta) is not None


def _examples(lst, n=3):
    return lst[:n]


def scan_fasta(path):
    """Read the FASTA once. Returns counts and a few example headers per rule; never raises on
    content (raises OSError if the file cannot be read -- the caller reports that)."""
    r = {"n_entries": 0, "test_header": [], "variant_header": [], "j_entries": [],
         "u_unexplained": [], "u_explained": 0, "o_unexplained": [], "o_explained": 0,
         "bz_entries": [], "decoy_orphans": []}
    head, seq = None, []
    targets, decoys = set(), {}

    def flush():
        if head is None:
            return
        r["n_entries"] += 1
        s = "".join(seq).upper()
        ident = head[1:].split(None, 1)[0] if head[1:].strip() else ""
        if ident.lower().startswith(DECOY_TAG):
            decoys[ident[len(DECOY_TAG):]] = head
        else:
            targets.add(ident)
        if HEADER_TEST_RE.search(head) or ID_TEST_RE.search(ident):
            r["test_header"].append(head)
        if HEADER_VARIANT_RE.search(head):
            r["variant_header"].append(head)
        if any(c in s for c in HARD_RESIDUES):
            r["j_entries"].append(head)
        if "U" in s:
            if SELENO_HEADER.search(head):
                r["u_explained"] += 1
            else:
                r["u_unexplained"].append(head)
        if "O" in s:
            if PYRROLYSINE_HEADER.search(head):
                r["o_explained"] += 1
            else:
                r["o_unexplained"].append(head)
        if any(c in s for c in AMBIGUOUS_RESIDUES):
            r["bz_entries"].append(head)

    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                flush()
                head, seq = line, []
            elif head is not None:
                seq.append(line.strip())
    flush()
    r["decoy_orphans"] = [h for i, h in decoys.items() if i not in targets]
    r["n_decoy_pairs"] = len(decoys) - len(r["decoy_orphans"])
    return r


def _apply_rules(s, label, hard, soft, detail):
    """The content rules over one scan_fasta() result (the database itself, or a file added to it:
    `label` prefixes the reason and the detail key)."""
    for key, why in (("test_header", "header matches ENTRAP_, or its id says decoy/shuffle"),
                     ("variant_header", "header carries a variant tag (_VAR_/_CON_/_ENC_/_ENT_<n>)"),
                     ("decoy_orphans", f"a {DECOY_TAG} entry with no target in the file"),
                     ("j_entries", "J (not a residue) in a sequence")):
        if s[key]:
            hard.append(f"{label}{len(s[key])} entries: {why}")
            detail[label + key] = {"n": len(s[key]), "examples": _examples(s[key])}
    bz = s["bz_entries"]
    if bz and len(bz) <= BZ_TOLERANCE * s["n_entries"]:
        detail[label + "bz_tolerated"] = len(bz)
        bz = []
    for key, rows, why in (("u_unexplained", s["u_unexplained"], "U in an entry that names no selenoprotein"),
                           ("o_unexplained", s["o_unexplained"], "O in an entry that names no pyrrolysine protein"),
                           ("bz_entries", bz, f"B or Z in more than {BZ_TOLERANCE:.1%} of entries")):
        if rows:
            soft.append(f"{label}{len(rows)} entries: {why}")
            detail[label + key] = {"n": len(rows), "examples": _examples(rows)}


def _fetch_meta(fasta, meta_path=None):
    """The fetch_fasta.py sidecar dict that describes this database, or None: `meta_path` (the one
    the search's provenance recorded), `<fasta>.meta.json`, or for Philosopher's .fas the target's
    (fetch_fasta.sidecar_for, the one definition of that rule)."""
    for p in (meta_path, fasta + ".meta.json"):
        d = safe_write.load_json(p) if p else None
        if isinstance(d, dict) and ("organism_source" in d or "organism" in d):
            return d
    try:
        import fetch_fasta
        sc = fetch_fasta.sidecar_for(fasta, look_in=(os.path.dirname(os.path.dirname(fasta)),
                                                     os.path.join(os.path.dirname(os.path.dirname(
                                                         fasta)), "input")))
    except Exception:                                   # noqa: BLE001 -- unresolved = ask a person
        return None
    d = safe_write.load_json(sc["path"]) if sc.get("path") else None
    return d if isinstance(d, dict) and ("organism_source" in d or "organism" in d) else None


def classify(path, check_sidecar=True, sidecar_known=None, meta_path=None):
    """(verdict, reasons, detail) for one FASTA. verdict in ok | test_database | needs_agent_check.
    `reasons` are short strings; `detail` carries counts and up to 3 example headers per rule.
    `sidecar_known`: why a fetch_fasta.py sidecar may be taken as existing although it is not beside
    the FASTA now (the search's own provenance recorded one: run_search.py `fasta_meta.path`, whose
    file is `meta_path`: read for the sequences added with --add-fasta)."""
    reasons, hard, soft, detail = [], [], [], {"fasta": path}
    mark = read_nofran(path) if check_sidecar else None
    if mark:
        hard.append(f"{os.path.basename(path)}{SIDECAR_SUFFIX} says corpus_eligible false"
                    + (f" ({mark.get('kind')}: {mark.get('why')})" if not mark.get("unreadable")
                       else " (unreadable marker)"))
        detail["nofran_sidecar"] = mark
    elif check_sidecar and not os.path.lexists(sidecar_path(path)):
        # The FASTA's own fetch_fasta.py record names a no-FRAN marker that is not beside it: a
        # --test-database FASTA copied or moved with its .meta.json but not its .nofran.json. Its
        # kinds (recoded, denovo_candidate, other) have no content signature, so the record is the
        # only evidence left -- hard, as the marker itself is (C2 review ux M7). A marker that IS
        # there and says corpus_eligible true is a person's decision and is not overruled here.
        named = (_fetch_meta(path, meta_path) or {}).get("nofran_sidecar")
        if named:
            hard.append(f"its fetch_fasta.py record ({os.path.basename(path)}.meta.json) says it "
                        f"was made as a test database with a no-FRAN marker ({named}), and that "
                        f"marker is not beside it (copied without it?)")
            detail["meta_nofran_sidecar"] = named
    try:
        s = scan_fasta(path)
    except OSError as e:
        # An unreadable database cannot be proven real: ask, do not refuse.
        soft.append(f"cannot read {path} to check it ({e.strerror or e})")
        return ("test_database" if hard else "needs_agent_check"), hard + soft, detail
    detail["n_entries"] = s["n_entries"]
    _apply_rules(s, "", hard, soft, detail)
    if s["n_entries"] == 0:
        soft.append("the FASTA has no entries")
    if sidecar_known and _fetch_meta(path) is None:
        detail["sidecar_from_provenance"] = sidecar_known
    if not mark and not sidecar_known and _fetch_meta(path) is None:
        soft.append("no fetch_fasta.py sidecar (a hand-built database)")
        detail["no_fetch_fasta_sidecar"] = True
    # Sequences the user ADDED (fetch_fasta --add-fasta: a bait, a tag) carry the sidecar's blessing
    # but nobody has looked at them: scan them with the same rules, and ask once either way.
    meta = _fetch_meta(path, meta_path)
    for af in [x.get("file") for x in (meta or {}).get("added_sequences") or [] if isinstance(x, dict)]:
        label = f"added file {os.path.basename(str(af))}: "
        soft.append(f"{label}sequences were added with fetch_fasta --add-fasta ({af}); a person "
                    f"confirms they are real (a bait or a tag), not shuffled or recoded")
        detail.setdefault("added_files", []).append(str(af))
        try:
            _apply_rules(scan_fasta(str(af)), label, hard, soft, detail)
        except OSError as e:
            soft.append(f"{label}cannot read it to check it ({e.strerror or e})")
    if hard:
        return "test_database", hard, detail
    if soft:
        return "needs_agent_check", soft, detail
    return "ok", [], detail



# ============================================================================== E: detection ==
# What an agent runs to find out whether anything synthetic is ALREADY in the corpus. Nothing here
# touches PG Farm: audit_sql() only PRINTS statements for a person with a database login (agents never
# mint or use PG Farm tokens); lance_scan() reads the Lance spectrum corpus on HIVE, read-only.
SELENO_GENE_PG = r"^(SEL[A-Z0-9]+|SEP[A-Z0-9]+|GPX[0-9]+|TXNRD[0-9]|DIO[0-9]|MSRB1|EPT1)"
SELENO_GENE_RE = re.compile(SELENO_GENE_PG, re.I)
ACCESSION_ALTERNATIVES = {**{alt: re.compile(alt, re.I) for alt in ACCESSION_TEST_CI.split("|")},
                          **{alt: re.compile(alt) for alt in ACCESSION_TEST_CS.split("|")}}


def is_synthetic_accession(group):
    """The ONE accession rule: a protein group naming an entrapment / shuffled / variant entry."""
    return any(rx.search(group) for rx in ACCESSION_ALTERNATIVES.values())
LANCE_COLUMNS = ("search_id", "stripped_seq", "protein_group", "genes")


def _pg(pattern):
    return pattern.replace(r"\d", "[0-9]")


def audit_sql():
    """[{name, expect, sql}] -- the audit, as SQL for a person to run. Every statement is a SELECT
    (a read-only audit) and every `expect` is what a clean corpus returns."""
    acc_ci, acc_cs = _pg(ACCESSION_TEST_CI), _pg(ACCESSION_TEST_CS)
    return [
        {"name": "peptides_with_J", "expect": "0 rows",
         "sql": "SELECT search_id, count(*) AS n FROM delimp_precursors "
                "WHERE stripped_seq ~ '[J]' GROUP BY search_id ORDER BY n DESC;"},
        {"name": "peptides_with_O_B_Z", "expect": "a handful at most (real in UniProt)",
         "sql": "SELECT search_id, count(DISTINCT stripped_seq) AS n FROM delimp_precursors "
                "WHERE stripped_seq ~ '[OBZ]' GROUP BY search_id ORDER BY n DESC;"},
        {"name": "U_peptides_outside_selenoprotein_genes", "expect": "0 rows, or rows a person "
         "recognises as selenoproteins (the gene regex is the human/mouse family list)",
         "sql": "SELECT p.search_id, p.stripped_seq, pr.gene FROM delimp_precursors p "
                "LEFT JOIN delimp_proteins pr ON pr.search_id = p.search_id "
                "AND pr.protein_group = p.protein_group "
                f"WHERE p.stripped_seq ~ 'U' AND coalesce(pr.gene, '') !~* '{_pg(SELENO_GENE_PG)}' "
                "LIMIT 200;"},
        {"name": "synthetic_protein_groups", "expect": "0 rows",
         "sql": "SELECT search_id, protein_group FROM delimp_proteins "
                f"WHERE protein_group ~* '{acc_ci}' OR protein_group ~ '{acc_cs}' LIMIT 200;"},
        {"name": "searches_by_database_md5", "expect": "every md5 a person recognises as a real "
         "database; a bad md5 is purged by search_id list",
         "sql": "SELECT fasta_md5, fasta_path, count(*) AS n_searches FROM delimp_searches "
                "GROUP BY fasta_md5, fasta_path ORDER BY n_searches DESC;"},
    ]


def lance_scan(root, max_datasets=None, budget_s=None):
    """Read-only scan of the Lance spectrum corpus under `root` (one `<search>.lance` per search).
    Reads only LANCE_COLUMNS. -> {"n_datasets", "n_scanned", "rows", "J": n, "OBZ": n, "U_outside":
    n, "synthetic_accessions": n, "hits": [per-dataset dicts, hits only], "errors": [...],
    "complete": bool}. Needs the `lance` and `pyarrow` modules (HIVE: /quobyte/proteomics-grp/brett/envs/alphadia2 (pylance 8; casanovo5 has 0.15, which cannot read the corpus))."""
    import time
    try:
        import lance
        import pyarrow.compute as pc
    except ImportError as e:
        raise SystemExit(f"lance_scan needs the lance and pyarrow modules ({e}); on HIVE use the "
                         f"alphadia2 env, inside a compute job")
    t0 = time.time()
    names = sorted(n for n in os.listdir(root) if n.endswith(".lance"))
    out = {"root": root, "n_datasets": len(names), "n_scanned": 0, "rows": 0, "J": 0, "OBZ": 0,
           "U_outside": 0, "synthetic_accessions": 0, "hits": [], "errors": [], "complete": True,
           # per alternative of ACCESSION_TEST: datasets hit, distinct groups, examples
           "by_alternative": {alt: {"datasets": 0, "groups": 0, "examples": []}
                              for alt in ACCESSION_ALTERNATIVES}}
    for n in names[:max_datasets] if max_datasets else names:
        if budget_s and time.time() - t0 > budget_s:
            out["complete"] = False
            break
        try:
            ds = lance.dataset(os.path.join(root, n))
            cols = [c for c in LANCE_COLUMNS if c in ds.schema.names]
            t = ds.to_table(columns=cols)
            out["rows"] += t.num_rows
            seq = t["stripped_seq"]
            row = {"dataset": n, "rows": t.num_rows}
            j = pc.sum(pc.match_substring_regex(seq, "J").cast("int64")).as_py() or 0
            obz = pc.sum(pc.match_substring_regex(seq, "[OBZ]").cast("int64")).as_py() or 0
            u_mask = pc.match_substring_regex(seq, "U")
            u_out = 0
            u_ex = []
            if "genes" in t.schema.names:
                u_rows = t.filter(u_mask)
                for seqv, gv in zip(u_rows["stripped_seq"].to_pylist(), u_rows["genes"].to_pylist()):
                    if not SELENO_GENE_RE.search(gv or ""):
                        u_out += 1
                        if len(u_ex) < 3:
                            u_ex.append([seqv, gv])
            else:
                u_out = pc.sum(u_mask.cast("int64")).as_py() or 0
            syn = 0
            syn_ex = []
            if "protein_group" in t.schema.names:
                groups = pc.unique(t["protein_group"]).to_pylist()
                bad = [g for g in groups if g and is_synthetic_accession(g)]
                syn, syn_ex = len(bad), bad[:3]
                for alt, rx in ACCESSION_ALTERNATIVES.items():
                    hit = [g for g in bad if rx.search(g)]
                    if hit:
                        b = out["by_alternative"][alt]
                        b["datasets"] += 1
                        b["groups"] += len(hit)
                        if len(b["examples"]) < 5:
                            b["examples"] += [[n, g] for g in hit[:2]][:5 - len(b["examples"])]
            out["J"] += j
            out["OBZ"] += obz
            out["U_outside"] += u_out
            out["synthetic_accessions"] += syn
            out["n_scanned"] += 1
            if j or obz or u_out or syn:
                row.update(J=j, OBZ=obz, U_outside=u_out, U_examples=u_ex,
                           synthetic_accessions=syn, accession_examples=syn_ex)
                out["hits"].append(row)
        except Exception as e:                                       # noqa: BLE001
            out["errors"].append({"dataset": n, "error": f"{type(e).__name__}: {str(e)[:160]}"})
    return out


def main(argv=None):
    """  fasta_guard.py check <fasta>...                      verdict per FASTA (read-only)
      fasta_guard.py mark <fasta> --kind K --why "..."     write <fasta>.nofran.json (for a script
                                                          or a person who made a test database)"""
    import argparse
    ap = argparse.ArgumentParser(description=main.__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("check")
    c.add_argument("fasta", nargs="+")
    m = sub.add_parser("mark")
    m.add_argument("fasta")
    m.add_argument("--kind", required=True, choices=KINDS)
    m.add_argument("--why", required=True)
    m.add_argument("--made-by", default="fasta_guard.py mark")
    a = ap.parse_args(argv)
    if a.cmd == "mark":
        print(mark_nofran(a.fasta, a.kind, a.made_by, a.why))
        return 0
    out = []
    for f in a.fasta:
        v, reasons, d = classify(f)
        out.append({"fasta": f, "verdict": v, "reasons": reasons, "n_entries": d.get("n_entries")})
    print(json.dumps(out, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
