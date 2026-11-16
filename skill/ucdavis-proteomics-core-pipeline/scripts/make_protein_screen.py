#!/usr/bin/env python3
"""Find one protein in an analysis's screen_data/ folder, check the data, assemble the page data.

    make_protein_screen.py --screen-data <dir> --protein <gene|accession|Protein.Group>
                           [--contrast NAME] --out <file.html> [--offline] [--dump-data FILE]
    make_protein_screen.py --screen-data <dir> --top N [--direction up|down|both] [--contrast NAME]
                           --out-dir <dir> [--offline]

--top writes <out-dir>/index.html (the list: the N most significant proteins the comparison calls changed,
ranked by p-value) and one <out-dir>/<label>.html per protein, the same page --protein writes.

Reads the folder make_screen_data.R wrote (manifest.json, index.tsv, parts/part-NNNN.json.gz), refuses
when a checksum does not match, looks the protein up in UniProt, NCBI, HGNC, Ensembl and the Human Protein
Atlas (protein_facts.py is the only code that talks to them), and builds `D`, the object the page template renders. Exit 0 = written; exit 2 = refused, the reason
on stderr (unknown or ambiguous protein, a changed file, supported: false, no manifest.json).

`peptide_of` is the one place a DIA-NN precursor id becomes a peptide sequence; import it, don't copy it.
Standard library only; must run on Python 3.9.
"""
import argparse
import csv
import gzip
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from protein_facts import lookup, offline_facts  # noqa: E402  (the one place a protein fact is fetched)
from protein_ids import first_accession, strip_isoform  # noqa: E402  (the one reading of an accession and its isoform suffix)
from safe_write import sha256_file  # noqa: E402  (THE chunked sha256, not a copy of it here: repo rule 3)

TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "protein_screen.tmpl.html")
DATA_TOKEN = "__SCREEN_DATA__"
# Where the files the page names are, relative to output/. The manifest says (`files`: the report is the
# one the DE read, wherever it is); these are the places screen data written before it said so read from.
FILES = {"session_rds": "tables/DE-LIMP_session.rds", "report": "search/report.parquet",
         "stats": "search/report.stats.tsv"}
CHANGED = ("the analysis changed after its screen data was made; rebuild it with make_screen_data.R")


class Refused(Exception):
    """The request cannot be answered safely; main() prints the message and exits 2."""


# ---- precursor -> peptide ----------------------------------------------------------------

def peptide_of(precursor_id):
    """DIA-NN Precursor.Id -> peptide sequence: drop `(UniMod:n)` and the trailing charge digits."""
    return re.sub(r"\d+$", "", re.sub(r"\(UniMod:\d+\)", "", precursor_id))


# ---- the folder ----------------------------------------------------------------------------

def _no_manifest_reason(screen_dir):
    """Why a screen_data folder has no manifest.json, from make_screen_data.R's run record
    (output/screen_data_status.json, beside the folder): what happened, never a guess. The record is
    absent only for screen data that was never built here, or built before the record existed."""
    head = "%s has no manifest.json" % screen_dir
    rec = None
    try:
        with open(os.path.join(os.path.dirname(os.path.abspath(screen_dir)), "screen_data_status.json")) as fh:
            rec = json.load(fh)
    except (OSError, ValueError):
        pass
    state = rec.get("state") if isinstance(rec, dict) else None
    why = " ".join(str(rec.get("reason") or "").split()) if isinstance(rec, dict) else ""
    route = rec.get("route") if isinstance(rec, dict) else None
    again = "Run make_screen_data.R on it again (step 8g)."
    if state == "failed":
        return "%s: the last build of the screen data failed (%s). %s" % (head, why or "no reason recorded", again)
    if state == "running":
        return ("%s: the last build of the screen data started and did not finish (still running, or it was "
                "killed, ran out of memory or timed out; look at its job log). %s" % (head, again))
    if state == "refused":
        return "%s: the last build of the screen data was refused%s: %s" % (
            head, " (%s)" % route if route else "", why or "no reason recorded")
    if state == "built":
        return ("%s although its run record says it was built: the folder was changed or removed since. %s"
                % (head, again))
    return ("%s and no run record: the screen data was never built for this analysis (step 8g not run), or it "
            "was built by an older version of the skill that kept no record. Run make_screen_data.R on it "
            "(step 8g)." % head)


def load_manifest(screen_dir):
    path = os.path.join(screen_dir, "manifest.json")
    if not os.path.isfile(path):
        raise Refused(_no_manifest_reason(screen_dir))
    damaged = "%s is %s: the screen data is damaged. Rebuild it with make_screen_data.R."
    try:
        with open(path) as fh:
            m = json.load(fh)
    except ValueError as e:
        raise Refused(damaged % (path, "not valid JSON (%s)" % e))
    if not isinstance(m, dict):
        raise Refused(damaged % (path, "not a JSON object"))
    if m.get("supported") is False:
        a = m.get("analysis") or {}
        route = (" It used %s (pipeline_id %s)." % (a.get("display_label") or "an unnamed quantification",
                                                   a.get("pipeline_id") or "not recorded")
                 if a.get("pipeline_id") or a.get("display_label") else "")
        raise Refused("the protein screen does not cover this analysis: %s.%s No page was made; say so, and "
                      "do not describe the protein's data another way." %
                      ((m.get("unsupported_reason") or "no reason recorded").rstrip("."), route))
    return m


def parse_index(path):
    """index.tsv -> list of dicts (Protein.Group, Genes, Protein.Names, part, n_precursors, ...).

    The one place the index is parsed, so columns added later arrive here as extra keys. make_screen_data.R
    writes it with quote = FALSE, so it is read without quoting: a `"` in a name is plain text.
    """
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t", quoting=csv.QUOTE_NONE))


def check_sources(manifest, output_dir):
    """(verified, not_checked, mismatched) relpaths of manifest.sources under `output_dir`."""
    verified, not_checked, bad = [], [], []
    for rel, sha in (manifest.get("sources") or {}).items():
        p = os.path.join(output_dir, rel)
        if not os.path.isfile(p):
            not_checked.append(rel)
        elif sha256_file(p) == sha:
            verified.append(rel)
        else:
            bad.append(rel)
    return verified, not_checked, bad


def verify(screen_dir, manifest, output_dir=None, sources=None):
    """The problems found: index or part files missing or altered, source files changed. [] = clean.

    Source files are read from `output_dir` (default: the folder holding screen_dir, where
    make_screen_data.R puts screen_data/); those that are not there are skipped here and listed by
    check_sources() as not checked. `sources`, check_sources()'s result when the caller has it already,
    keeps a source from being hashed twice (report.parquet can be gigabytes).
    """
    problems = []
    idx = manifest.get("index") or "index.tsv"
    if not os.path.isfile(os.path.join(screen_dir, idx)):
        problems.append("%s is missing" % idx)
    for rel, sha in (manifest.get("parts") or {}).items():
        p = os.path.join(screen_dir, rel)
        if not os.path.isfile(p):
            problems.append("%s is missing" % rel)
        elif sha256_file(p) != sha:
            problems.append("%s does not match its checksum in manifest.json" % rel)
    if sources is None:
        sources = check_sources(manifest, output_dir or os.path.dirname(os.path.abspath(screen_dir)))
    for rel in sources[2]:
        problems.append("%s: %s" % (rel, CHANGED))
    return problems


def resolve_protein(index_rows, query):
    """Index rows whose gene (any of `;`-separated), accession (any in the group) or whole group is `query`."""
    q = query.strip().lower()
    hits = []
    for r in index_rows:
        pg = r["Protein.Group"].lower()
        names = [g.lower() for g in (r.get("Genes") or "").split(";") if g] + pg.split(";")
        if q == pg or q in names:
            hits.append(r)
    return hits


def read_record(screen_dir, index_row, manifest):
    """The protein's record from its part. The part must be one manifest.json lists (and verify() hashed), and
    the record filed under the protein's key must name that protein (its `pg`)."""
    rel = "parts/" + index_row["part"]
    if rel not in (manifest.get("parts") or {}):
        raise Refused("index.tsv puts %s in %s, which manifest.json does not list; its checksum was never "
                      "recorded. %s" % (index_row["Protein.Group"], rel, CHANGED))
    with gzip.open(os.path.join(screen_dir, rel), "rt") as fh:
        part = json.load(fh)
    pg = index_row["Protein.Group"]
    if pg not in part:
        raise Refused("%s is not in %s although index.tsv says it is; the screen data is inconsistent. %s"
                      % (pg, rel, CHANGED))
    got = part[pg].get("pg") if isinstance(part[pg], dict) else None
    if got != pg:
        raise Refused("%s in %s holds the record of %s, not of %s: the screen data is inconsistent (a worker of "
                      "make_screen_data.R may have died). %s" % (pg, rel, got, pg, CHANGED))
    return part[pg]


# ---- page data -------------------------------------------------------------------------------

def _shared(pg, record, index_rows):
    """{peptide: {proteotypic, others: [{acc, gene}]}} from the report rows' Protein.Ids; genes from index.tsv."""
    gene_of = {}
    for r in index_rows or []:
        accs, genes = r["Protein.Group"].split(";"), (r.get("Genes") or "").split(";")
        for n, a in enumerate(accs):
            gene_of[a] = genes[n] if len(genes) == len(accs) else (r.get("Genes") or None)
    own = set(pg.split(";"))
    shared = {}
    for r in record.get("rows", []):
        s = shared.setdefault(peptide_of(r["Precursor.Id"]),
                              {"proteotypic": r.get("Proteotypic"), "others": []})
        for a in (r.get("Protein.Ids") or "").split(";"):
            if a and a not in own and all(o["acc"] != a for o in s["others"]):
                s["others"].append({"acc": a, "gene": gene_of.get(a)})
    return shared


def lookup_accession(pg):
    """The accession of protein group `pg` whose facts are looked up: its first member (protein_ids.first_accession;
    the one place this choice is made), isoform suffix kept."""
    return first_accession(pg)


def _isoform_reason(pg):
    """Why the peptides are not placed on UniProt's sequence, when the looked-up accession is an isoform; else None."""
    acc = lookup_accession(pg)
    if strip_isoform(acc) == acc:
        return None
    return ("this protein group's accession is the isoform %s, and UniProt's entry holds the canonical sequence of %s, "
            "not the isoform's" % (acc, strip_isoform(acc)))


def _info(pg, record, uniprot, index_rows):
    """The card's UniProt block. Coverage and positions are None when the group's accession is an isoform."""
    cov, pos = None, None
    if _isoform_reason(pg) is None:
        pep_order = list(dict.fromkeys([peptide_of(p["id"]) for p in record["precursors"]]
                                       + [peptide_of(r["Precursor.Id"]) for r in record.get("rows", [])]))
        seq = uniprot["sequence"]
        used = {peptide_of(p["id"]) for p in record["precursors"]}
        cov, pos = [0] * len(seq), {}
        for pep in pep_order:
            i = seq.find(pep)
            pos[pep] = [i + 1, i + len(pep)] if i >= 0 else None
            if i >= 0 and pep in used:
                for k in range(i, i + len(pep)):
                    cov[k] = 1
    shared = _shared(pg, record, index_rows)
    info = {k: uniprot.get(k) for k in ("acc", "id", "reviewed", "name", "gene", "organism", "length", "function",
                                         "location", "uniprot_version", "uniprot_updated", "release")}
    info.update({"fetched_at": uniprot.get("fetched_at"), "cov": cov, "pos": pos, "shared": shared})
    return info


def _share_method(manifest):
    """How the shares were made, as manifest.shares records it (its per-protein withheld list left out).

    The recalculation check is passed on as `check`: the page describes it in its own words.
    """
    sh = manifest.get("shares") or {}
    out = {k: sh.get(k) for k in ("method", "limpa_version", "tolerance") if k in sh}
    if "nudge_check" in sh:
        out["check"] = sh["nudge_check"]
    return out


def _withheld(block):
    if isinstance(block, dict) and "withheld" in block:
        return str(block["withheld"])
    return None


def chosen_contrast(manifest, name):
    """(the comparison named `name`, True when nobody chose it): with no name, the analysis's first comparison,
    a default the page tags (repo rule 2) when there were 2 or more to choose from."""
    cons = manifest["contrasts"]
    c = next((x for x in cons if x["name"] == (name or cons[0]["name"])), None)
    if c is None:
        raise Refused("no comparison %r in this analysis; choose one of: %s" % (name, ", ".join(x["name"] for x in cons)))
    return c, name is None and len(cons) > 1


def page_data(manifest, record, contrast, uniprot, index_rows=None, verified=None, not_checked=None,
              uniprot_reason=None, facts=None):
    """The template's `D` for one protein record, in one contrast. `uniprot` None = no UniProt block.

    `facts` (protein_facts.Fact list) go to info[pg].facts; with no UniProt block they go to
    facts_no_uniprot[pg], which the card reads when it has no info.
    """
    pg = record["pg"]
    c, by_default = chosen_contrast(manifest, contrast)
    cname = c["name"]
    de = (record.get("de") or {}).get(cname) or {}
    tested = bool(de.get("tested"))
    vals = record.get("values") or {}
    ids = [p["id"] for p in record["precursors"]]
    table = c.get("table") or ""
    proteins = {pg: {
        "gene": record.get("genes"), "name": record.get("names"), "tested": tested,
        "logFC": de.get("logFC"), "P": de.get("P.Value"), "adjP": de.get("adj.P.Val"),
        "CI": [de.get("CI.L"), de.get("CI.R")], "t": de.get("t"), "df_total": de.get("df_total"),
        "evidence": de.get("evidence"), "det": de.get("detected") or {}, "nprec": record.get("n_precursors"),
        "prot": vals.get("protein"), "se": vals.get("se"), "nobs": vals.get("n_obs"),
        "prec_ids": ids, "prec": [p["log2"] for p in record["precursors"]], "sigma": record.get("sigma")}}

    missing = {m["block"]: m["reason"] for m in manifest.get("skipped") or []}
    contrib, detprob, routes = {}, {}, {}
    sh = record.get("shares")
    if isinstance(sh, dict) and "own" in sh:
        contrib[pg] = {"ids": ids, "own": sh["own"], "other": sh["other"], "total": sh["total"],
                       "refit_max_diff": sh.get("refit_max_diff")}
    else:
        missing["shares"] = _withheld(sh) or "this protein's shares are not in the screen data"
    det = record.get("detection")
    if isinstance(det, dict) and "pf_fit" in det:
        detprob[pg] = {"ids": ids, "pf_fit": det["pf_fit"], "pf_obs": det["pf_obs"]}
    else:
        missing["detection"] = (_withheld(det) or (missing.get("shares") if "shares" in missing else None)
                                or "this protein's detection chances are not in the screen data")
    w = record.get("weights")
    if isinstance(w, dict) and "w" in w:
        routes[pg] = {"w": w["w"], "wshare": w.get("share_in_group"), "df_res": w.get("df_residual"),
                      "df_prior": w.get("df_prior")}
    elif "weights" not in missing:
        missing["weights"] = "this protein's weights are not in the screen data"
    info, facts_no_uniprot = {}, {}
    if uniprot is not None:
        info[pg] = _info(pg, record, uniprot, index_rows)
        info[pg]["facts"] = facts or []
        if _isoform_reason(pg):
            missing["coverage"] = _isoform_reason(pg)
    else:
        if facts:
            facts_no_uniprot[pg] = facts
        missing["uniprot"] = uniprot_reason or "UniProt was not consulted (offline) or did not answer"

    rows = []
    for r in record.get("rows", []):
        r = dict(r)
        r["Protein.Group"] = pg
        rows.append(r)
    files = dict(FILES)
    files.update({k: v for k, v in (manifest.get("files") or {}).items() if isinstance(v, str) and v})
    files["de_table"] = table if table.startswith("tables/") else "tables/" + table
    return {
        "samples": manifest["samples"], "groups": manifest["groups"], "dpc": manifest.get("dpc"),
        "prior": manifest.get("prior"), "n_prec": manifest.get("n_precursors"),
        "n_tested": c.get("n_tested"), "n_sig": c.get("n_sig"),
        "contrast": {"name": cname, "num": c.get("numerator"), "den": c.get("denominator"), "adjp": c.get("adjp"),
                     "default": by_default},
        "proteins": proteins, "rows": rows, "contrib": contrib, "detprob": detprob, "routes": routes,
        "stats": manifest.get("run_stats") or {}, "info": info, "facts_no_uniprot": facts_no_uniprot, "xic": {}, "byrow": {},
        "files": files, "analysis": manifest.get("analysis"), "missing": missing,
        "shared": {pg: _shared(pg, record, index_rows)}, "share_method": _share_method(manifest),
        "lookup": {"accession": lookup_accession(pg), "accessions": pg.split(";")},
        "peptides": {pid: peptide_of(pid) for pid in ids + [r["Precursor.Id"] for r in record.get("rows", [])]},
        "provenance": {"written_by": manifest.get("written_by"), "written_at": manifest.get("written_at"),
                       "verified": list(verified or []), "not_checked": list(not_checked or [])}}


def render(template_text, D):
    """Put D into the template at DATA_TOKEN as inert JSON (no `</script>` can end the script early)."""
    # every `<` becomes the JSON escape \u003c: covers `</script>` and `<!--`, and JSON.parse reads it back
    blob = json.dumps(D, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c")
    return template_text.replace(DATA_TOKEN, blob)


# ---- the top N list ------------------------------------------------------------------------------

TOP_MAX = 100
LIST_TEMPLATE = os.path.join(HERE, "protein_screen_list.tmpl.html")
INDEX_COLS = ("P.Value", "adj.P.Val", "logFC")


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if v == v and v not in (float("inf"), float("-inf")) else None


def top_proteins(index_rows, contrast, n, direction, adjp):
    """-> (rows, n_called): the proteins `contrast` calls changed, best first, at most `n` of them.

    Called = adj.P.Val < adjp (strict, as run_de.R's significance_rule); there is no fold-change cut.
    Ranked by adjusted p ascending, then P.Value (BH gives many proteins one adjusted p; the raw p breaks
    that tie), then Protein.Group, so the order is the same every time. `direction` "up"/"down" keeps logFC > 0 / < 0 first; n_called counts
    the called proteins that pass it. An untested protein (empty cells) is never called, and an uncalled
    protein never fills the list. Raises Refused when the index lacks the contrast's columns or adjp is None.
    """
    if adjp is None:
        raise Refused("this analysis's manifest.json records no adjusted-p cutoff (adjp) for %s, so no protein "
                      "can be called changed. Rebuild the screen data with make_screen_data.R." % contrast)
    cols = ["%s__%s" % (k, contrast) for k in INDEX_COLS]
    if index_rows and any(c not in index_rows[0] for c in cols):
        raise Refused("index.tsv has no %s columns for %s: this screen data was made before the top-N list. "
                      "Rebuild it with make_screen_data.R." % (" / ".join(INDEX_COLS), contrast))
    if direction not in ("both", "up", "down"):
        raise Refused("direction must be up, down or both")
    called = []
    for r in index_rows:
        adj, p, lfc = (_num(r.get(c)) for c in (cols[1], cols[0], cols[2]))
        if adj is None or p is None or lfc is None or not adj < adjp:
            continue
        if (direction == "up" and not lfc > 0) or (direction == "down" and not lfc < 0):
            continue
        called.append((adj, p, r["Protein.Group"], r))
    called.sort(key=lambda t: (t[0], t[1], t[2]))
    return [t[3] for t in called[:n]], len(called)


def file_labels(rows):
    """{Protein.Group: file stem}: the gene, else the leading accession, safe for a file name.

    A stem two proteins share (a gene name twice, or one that differs only by case, which a Mac or
    Windows disk treats as the same file) takes the leading accession as well; `index` is the list page's.
    """
    def accession(r):
        return r["Protein.Group"].split(";")[0]

    def safe(t):
        return re.sub(r"[^A-Za-z0-9._-]+", "_", t).strip("._-")[:60]
    base = {}
    for r in rows:
        genes = [g for g in (r.get("Genes") or "").split(";") if g.strip()]
        base[r["Protein.Group"]] = safe(genes[0]) if genes and safe(genes[0]) else safe(accession(r))
    count = {}
    for b in base.values():
        count[b.lower()] = count.get(b.lower(), 0) + 1
    out, used = {}, {"index"}
    for r in rows:
        pg = r["Protein.Group"]
        b = base[pg]
        if count[b.lower()] > 1 or b.lower() == "index":
            b = (b + "_" + safe(accession(r))).strip("_")
        stem, k = b, 2
        while stem.lower() in used:
            stem = "%s_%d" % (b, k)
            k += 1
        used.add(stem.lower())
        out[pg] = stem
    return out


def _recommended_name_of(facts):
    """(name, why not): UniProt's recommended name, verbatim from the UniProt Fact."""
    for f in facts or []:
        if f.get("source") == "UniProt":
            if f.get("status") == "ok":
                for label, value in f.get("fields") or []:
                    if label == "Protein name" and value:
                        return value, None
                return None, "UniProt gave no recommended name"
            return None, {"not_found": "not in UniProt", "not_reached": "not reached",
                          "not_applicable": "not applicable"}.get(f.get("status"), f.get("status"))
    return None, "not reached"


def _fact_note(pg):
    """Which UniProt entry names protein group `pg` on the list, when that is not plain: None for one canonical accession."""
    acc = lookup_accession(pg)
    first = ", the group's first accession" if len(pg.split(";")) > 1 else ""
    if strip_isoform(acc) != acc:
        return ("UniProt's canonical entry %s; the group's accession is the isoform %s%s"
                % (strip_isoform(acc), acc, first))
    return "UniProt entry of %s%s" % (acc, first) if first else None


def list_data(manifest, contrast, adjp, direction, n, shown, n_called, record_of, index_rows, facts_of, labels,
              provenance, by_default=False):
    """The list template's `D`: one entry per shown protein, from its record, its facts and its page's file name.
    `by_default`: nobody chose the comparison (chosen_contrast), so the list tags it."""
    cname = contrast["name"]
    num, den = contrast.get("numerator"), contrast.get("denominator")
    groups = manifest["groups"]
    out = []
    for rank, r in enumerate(shown, 1):
        pg = r["Protein.Group"]
        rec = record_of[pg]
        de = (rec.get("de") or {}).get(cname) or {}
        nobs = (rec.get("values") or {}).get("n_obs") or []
        inferred = any(g in (num, den) and k < len(nobs) and nobs[k] == 0 for k, g in enumerate(groups))
        shared = any(v["others"] for v in _shared(pg, rec, index_rows).values())
        name, why = _recommended_name_of(facts_of[pg])
        uf = next((f for f in facts_of[pg] or [] if f.get("source") == "UniProt"), None) if name else None
        lfc = _num(r["logFC__" + cname])
        out.append({"rank": rank, "pg": pg, "gene": r.get("Genes") or None, "names": r.get("Protein.Names") or None,
                    "logFC": lfc, "fold": 2.0 ** lfc, "P": _num(r["P.Value__" + cname]),
                    "adjP": _num(r["adj.P.Val__" + cname]),
                    "det": [[g, (de.get("detected") or {}).get(g)] for g in (den, num)],
                    "inferred": inferred, "shared": shared, "fact": name, "fact_why": why, "fact_note": _fact_note(pg),
                    "fact_src": ({"source": "UniProt", "release": uf.get("release"), "fetched_at": uf.get("fetched_at"),
                                  "url": uf.get("url"), "acc": next((v for l, v in uf.get("fields") or []
                                                                      if l == "Accession"), None)} if uf else None),
                    "link": labels[pg] + ".html"})
    return {"contrast": {"name": cname, "num": num, "den": den, "adjp": adjp, "default": by_default},
            "direction": direction,
            "asked": n, "called": n_called, "proteins": out, "analysis": manifest.get("analysis"),
            "n_tested": contrast.get("n_tested"), "provenance": provenance}


# ---- CLI -----------------------------------------------------------------------------------------

def _open_folder(args):
    """(manifest, index_rows, verified, not_checked): the folder, checked; each source file hashed once."""
    screen_dir = args.screen_data
    manifest = load_manifest(screen_dir)
    sources = check_sources(manifest, args.output_dir or os.path.dirname(os.path.abspath(screen_dir)))
    problems = verify(screen_dir, manifest, args.output_dir, sources=sources)
    if problems:
        raise Refused("the screen data cannot be trusted:\n  " + "\n  ".join(problems))
    index_rows = parse_index(os.path.join(screen_dir, manifest.get("index") or "index.tsv"))
    names = [x["name"] for x in manifest["contrasts"]]
    if args.contrast and args.contrast not in names:
        raise Refused("no comparison %r in this analysis; choose one of: %s" % (args.contrast, ", ".join(names)))
    verified, not_checked, _ = sources
    return manifest, index_rows, verified, not_checked


def facts_of_accession(acc, offline, cache):
    """(uniprot or None, facts, reason): looked up once per accession for the whole run, kept in `cache`."""
    if acc not in cache:
        if offline:
            cache[acc] = (None, offline_facts("--offline"), "offline: UniProt was not consulted")
        else:
            errs = []
            uniprot, facts = lookup(acc, errors=errs)
            reason = None
            if uniprot is None:
                reason = (errs[0] if errs else "UniProt did not answer") + " (page built without it)"
                print("warning: " + reason, file=sys.stderr)
            cache[acc] = (uniprot, facts, reason)
    return cache[acc]


def _page_text(D, template_path=TEMPLATE):
    if not os.path.isfile(template_path):
        raise Refused("the page template %s is not installed" % os.path.basename(template_path))
    with open(template_path, encoding="utf-8") as fh:
        return render(fh.read(), D)


def _protein_D(args, manifest, index_rows, verified, not_checked, row, record, cache):
    uniprot, facts, reason = facts_of_accession(lookup_accession(row["Protein.Group"]), args.offline, cache)
    return page_data(manifest, record, args.contrast, uniprot, index_rows=index_rows, verified=verified,
                     not_checked=not_checked, uniprot_reason=reason, facts=facts)


def _partial(manifest, screen_dir, index_rows):
    """The partial folder (make_screen_data.R --proteins) named for a message, or None for a full one."""
    if not manifest.get("partial"):
        return None
    return ("%s is a partial screen data folder (made by make_screen_data.R --proteins) that holds only %d proteins"
            % (screen_dir, len(index_rows)))


def _run(args):
    manifest, index_rows, verified, not_checked = _open_folder(args)
    hits = resolve_protein(index_rows, args.protein)
    if not hits:
        part = _partial(manifest, args.screen_data, index_rows)
        raise Refused("no protein matches %r in this analysis (gene, accession or Protein.Group)%s"
                      % (args.protein, "; but " + part + ": it may be in the analysis's full screen_data folder"
                         if part else ""))
    if len(hits) > 1:
        raise Refused("%r matches %d protein groups; give the accession or Protein.Group of one:\n  %s"
                      % (args.protein, len(hits), "\n  ".join("%s  (%s)  %s" % (
                          h["Protein.Group"], h.get("Genes") or "no gene", h.get("Protein.Names") or "")
                          for h in hits)))
    record = read_record(args.screen_data, hits[0], manifest)
    D = _protein_D(args, manifest, index_rows, verified, not_checked, hits[0], record, {})
    if args.dump_data:
        with open(args.dump_data, "w") as fh:
            json.dump(D, fh, allow_nan=False)
    if os.path.isfile(TEMPLATE):
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(_page_text(D))
    elif not args.dump_data:
        raise Refused("the page template %s is not installed" % os.path.basename(TEMPLATE))
    return 0


def _run_top(args):
    if args.top < 1:
        raise Refused("--top must be 1 or more")
    if args.top > TOP_MAX:
        raise Refused("ask for %d or fewer; the list page links each one (--top %d is too many pages to build)"
                      % (TOP_MAX, args.top))
    if args.out:
        raise Refused("--top writes into --out-dir (index.html and one page per protein); drop --out")
    if not args.out_dir:
        raise Refused("--top needs --out-dir")
    manifest, index_rows, verified, not_checked = _open_folder(args)
    part = _partial(manifest, args.screen_data, index_rows)
    if part:
        raise Refused(part + ", so it cannot list the proteins a comparison calls changed; use the analysis's full "
                      "screen_data folder")
    # args.contrast stays as given: each protein page picks (and tags) the same comparison through chosen_contrast
    contrast, by_default = chosen_contrast(manifest, args.contrast)
    adjp = contrast.get("adjp")
    shown, n_called = top_proteins(index_rows, contrast["name"], args.top, args.direction, adjp)
    labels = file_labels(shown)
    os.makedirs(args.out_dir, exist_ok=True)
    cache, record_of, facts_of = {}, {}, {}
    for k, row in enumerate(shown, 1):
        pg = row["Protein.Group"]
        print("protein %d of %d: %s" % (k, len(shown), pg), file=sys.stderr, flush=True)
        record_of[pg] = read_record(args.screen_data, row, manifest)
        D = _protein_D(args, manifest, index_rows, verified, not_checked, row, record_of[pg], cache)
        facts_of[pg] = cache[lookup_accession(pg)][1]
        with open(os.path.join(args.out_dir, labels[pg] + ".html"), "w", encoding="utf-8") as fh:
            fh.write(_page_text(D))
    prov = {"written_by": manifest.get("written_by"), "written_at": manifest.get("written_at"),
            "verified": list(verified), "not_checked": list(not_checked)}
    LD = list_data(manifest, contrast, adjp, args.direction, args.top, shown, n_called, record_of, index_rows,
                   facts_of, labels, prov, by_default=by_default)
    with open(os.path.join(args.out_dir, "index.html"), "w", encoding="utf-8") as fh:
        fh.write(_page_text(LD, LIST_TEMPLATE))
    print("%d of the %d proteins called changed in %s -> %s" % (len(shown), n_called, contrast["name"],
                                                                 os.path.join(args.out_dir, "index.html")),
          file=sys.stderr)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--screen-data", required=True, help="the analysis's output/screen_data folder")
    which = ap.add_mutually_exclusive_group(required=True)
    which.add_argument("--protein", help="gene, accession or Protein.Group")
    which.add_argument("--top", type=int, help="the N most significant proteins the comparison calls changed "
                                               "(at most %d)" % TOP_MAX)
    ap.add_argument("--direction", choices=("up", "down", "both"), default="both", help="with --top (default both)")
    ap.add_argument("--contrast", help="comparison to show (default: the first)")
    ap.add_argument("--out", help="with --protein: the HTML page to write")
    ap.add_argument("--out-dir", help="with --top: the folder for index.html and the protein pages")
    ap.add_argument("--offline", action="store_true", help="look nothing up (UniProt, NCBI, HGNC, Ensembl, HPA)")
    ap.add_argument("--output-dir", help="the analysis's output/ folder, where manifest sources are checked "
                                         "(default: the folder holding --screen-data)")
    ap.add_argument("--dump-data", help="also write the page data D as JSON to this file")
    args = ap.parse_args(argv)
    try:
        if args.top is not None:
            return _run_top(args)
        if not args.out:
            raise Refused("--protein needs --out")
        return _run(args)
    except Refused as e:
        print("make_protein_screen: " + str(e), file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
