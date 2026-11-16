#!/usr/bin/env python3
"""Protein and gene facts looked up from established sources: the one place a protein fact is fetched.

    facts_for("P80723")   -> [Fact, ...]   UniProt, NCBI Gene, HGNC, Ensembl, Human Protein Atlas
    lookup("P80723")      -> (UniProt fields for the page, [Fact, ...])
    fetch_uniprot("P80723")                 the UniProt fields the page's card shows, or None
    offline_facts()                         every source "not reached: --offline"

Brett's rule: look the protein up, never rely on memory. Every value in a Fact is the source's own
value, copied from the response, never rephrased and never filled in from what we think we know.
Where a source says nothing (HPA states no version; a null concentration) the Fact says so.

A Fact is {source, status, url, release, fetched_at, licence, fields: [[label, value]], note}.
`status` is "ok", "not_found", "not_reached" or "not_applicable"; only an "ok" Fact has fields.
The identifiers sent to a source come from UniProt's cross-references (GeneID, HGNC, Ensembl, HPA),
never from a search on a gene name: the Human Protein Atlas search matches loosely ("ALB" also
returns CALB1 and LALBA).

Hosts, headers and rate limits: references/protein-card-sources.md. Standard library only; Python 3.9.
"""
import datetime
import http.client
import json
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

from protein_ids import strip_isoform  # the one reading of an isoform suffix

ALLOWED_HOSTS = {"rest.uniprot.org", "eutils.ncbi.nlm.nih.gov", "rest.genenames.org", "rest.ensembl.org",
                 "www.proteinatlas.org"}
USER_AGENT = "ucdavis-proteomics-skill (+https://github.com/bsphinney/UCD-Proteomics-AI-Skill)"
TIMEOUT = 20
NCBI_HOST = "eutils.ncbi.nlm.nih.gov"
NCBI_TOOL = "ucdavis_proteomics_skill"
MIN_GAP = {NCBI_HOST: 0.34}          # seconds between calls to a host (NCBI: 3 requests/s without a key)
_LAST_CALL = {}

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb/"
# UniProt's own pattern for an accession, from https://rest.uniprot.org/help/accession_numbers (the help
# page "Accession", last modified 2026-06-12): "The three patterns can be combined into the following
# regular expression: [OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
UNIPROT_ACCESSION = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")

SOURCES = ["UniProt", "NCBI Gene", "HGNC", "Ensembl", "Human Protein Atlas"]
LICENCE = {
    "UniProt": "CC BY 4.0 (rest.uniprot.org/help/license)",
    "NCBI": "NCBI places no restrictions on use or distribution of its data; keep its disclaimer and "
            "acknowledge NLM (ncbi.nlm.nih.gov/home/about/policies/)",
    "HGNC": "CC0 (genenames.org/about/license/)",
    "Ensembl": "no licence statement found on the REST pages",
    "Human Protein Atlas": "CC BY 4.0 (proteinatlas.org/about/licence)",
}
HPA_FIELDS = ["Blood concentration - Conc. blood IM [pg/L]", "Blood concentration - Conc. blood MS [pg/L]",
              "Secretome location", "Evidence", "RNA blood cell specificity"]
HUMAN = 9606


# ---- the one network call ----------------------------------------------------------------------

def _check_url(url):
    p = urllib.parse.urlparse(url)
    if p.scheme != "https" or p.hostname not in ALLOWED_HOSTS or p.username or p.password:
        raise ValueError("%s is not an https URL on an allowed host (%s)" % (url, ", ".join(sorted(ALLOWED_HOSTS))))
    return p.hostname


class _AllowedRedirect(urllib.request.HTTPRedirectHandler):
    """Follow redirects only to https URLs on ALLOWED_HOSTS (a merged UniProt accession answers 303 to another path)."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        try:
            _check_url(newurl)
        except ValueError as e:
            raise urllib.error.URLError("redirect refused: %s" % e)
        return super().redirect_request(req, fp, code, msg, headers, newurl)


_OPENER = urllib.request.build_opener(_AllowedRedirect)


def _get(url, headers=None, timeout=TIMEOUT):
    """GET `url` -> (status, headers, body). Raises ValueError for a host not in ALLOWED_HOSTS.

    Header names come back lower-case, plus `x-final-url`: where the request ended after redirects.
    An HTTP error status is returned, not raised; a network failure raises URLError / OSError.
    Calls to NCBI are spaced at least 0.34 s apart.
    """
    host = _check_url(url)
    gap = MIN_GAP.get(host)
    if gap:
        last = _LAST_CALL.get(host)
        now = time.monotonic()
        if last is not None and now - last < gap:
            time.sleep(gap - (now - last))
        _LAST_CALL[host] = time.monotonic()
    h = {"User-Agent": USER_AGENT}
    h.update(headers or {})
    req = urllib.request.Request(url, headers=h)
    try:
        with _OPENER.open(req, timeout=timeout) as resp:
            final = resp.geturl()
            try:
                _check_url(final)
            except ValueError:
                raise urllib.error.URLError("redirected to %s; only %s are allowed" % (final, ", ".join(sorted(ALLOWED_HOSTS))))
            status = getattr(resp, "status", None) or resp.getcode()
            out = {k.lower(): v for k, v in resp.headers.items()}
            out["x-final-url"] = final
            return status, out, resp.read()
    except urllib.error.HTTPError as e:
        try:
            body = e.read()
        except (OSError, http.client.HTTPException):
            body = b""
        out = {k.lower(): v for k, v in (e.headers or {}).items()}
        out["x-final-url"] = url
        return e.code, out, body


def _ncbi_url(endpoint, **params):
    """An E-utilities URL. tool is always sent; email and api_key only from NCBI_EMAIL / NCBI_API_KEY."""
    import os
    q = list(params.items()) + [("tool", NCBI_TOOL)]
    if os.environ.get("NCBI_EMAIL"):
        q.append(("email", os.environ["NCBI_EMAIL"]))
    if os.environ.get("NCBI_API_KEY"):
        q.append(("api_key", os.environ["NCBI_API_KEY"]))
    return "https://%s/entrez/eutils/%s?%s" % (NCBI_HOST, endpoint, urllib.parse.urlencode(q))


# ---- Facts ----------------------------------------------------------------------------------------

def _today():
    return datetime.date.today().isoformat()


def _display_url(url):
    """`url` without the `email` and `api_key` parameters: a Fact's url ends up on the page and in dumps."""
    if not url or "?" not in url:
        return url
    base, query = url.split("?", 1)
    keep = [kv for kv in query.split("&") if kv.split("=", 1)[0] not in ("email", "api_key")]      # others stay byte for byte
    return base + ("?" + "&".join(keep) if keep else "")


def _fact(source, status, url=None, release=None, licence=None, fields=None, note=None):
    return {"source": source, "status": status, "url": _display_url(url), "release": release, "fetched_at": _today(),
            "licence": licence, "fields": fields or [], "note": note}


def offline_facts(reason="--offline"):
    return [_fact(s, "not_reached", note=reason) for s in SOURCES]


def _host_of(url):
    return urllib.parse.urlparse(url).hostname


SHAPE_ERRORS = (AttributeError, TypeError, KeyError, IndexError)


def _guarded(source, tuple_result=False):
    """A fetcher whose source answered valid JSON of an unexpected shape reports "not_reached", not a crash."""
    def deco(fn):
        def run(*a, **k):
            try:
                return fn(*a, **k)
            except SHAPE_ERRORS:
                f = _fact(source, "not_reached", note="the source answered with an unexpected response shape")
                return (f, None) if tuple_result else f
        run.__name__ = fn.__name__
        return run
    return deco


def _release(source, url, headers, pick):
    """The source's release from a second call -> the value, or "not retrieved: <why>" (never silently None)."""
    data, _, fail = _call_json(source, url, headers)
    if fail:
        return "not retrieved: " + (fail["note"] or fail["status"])
    try:
        rel = pick(data)
    except SHAPE_ERRORS:
        rel = None
    return str(rel) if rel not in (None, "") else "not retrieved: the answer held no release (unexpected response shape)"


def _call_json(source, url, headers=None):
    """-> (data, response headers, None) or (None, None, a failure Fact)."""
    try:
        status, hdrs, body = _get(url, headers)
    except (urllib.error.URLError, http.client.HTTPException, OSError) as e:
        return None, None, _fact(source, "not_reached", url,
                                 note="%s could not be reached (%s)" % (_host_of(url), e))
    if status in (404, 400):
        return None, None, _fact(source, "not_found", url, note="%s answered HTTP %s" % (_host_of(url), status))
    if status != 200:
        return None, None, _fact(source, "not_reached", url, note="%s answered HTTP %s" % (_host_of(url), status))
    try:
        return json.loads(body.decode("utf-8")), hdrs, None
    except ValueError:
        return None, None, _fact(source, "not_reached", url,
                                 note="%s answered with something that is not JSON" % _host_of(url))


def _rows(pairs, keep_empty=()):
    """[[label, value]] without the pairs the source left empty (labels in `keep_empty` stay, even as None)."""
    return [[k, v] for k, v in pairs if v not in (None, "", [], {}) or k in keep_empty]


# ---- UniProt --------------------------------------------------------------------------------------

def entry_accession(acc):
    """The accession whose UniProt entry describes `acc`: an isoform id (P80723-2) has no entry of its own,
    so its canonical entry's (P80723); the entry's sequence is the canonical one, not the isoform's. The
    suffix is read by protein_ids.strip_isoform, so a dashed id that is no isoform (ABC123-ALT) stays."""
    return strip_isoform(acc)


def _recommended_name(up):
    pd = up.get("proteinDescription") or {}
    return (pd.get("recommendedName") or (pd.get("submissionNames") or [{}])[0]).get("fullName", {}).get("value")


def _parse_entry(up, acc, release):
    """The UniProt fields the page's card shows, from a live entry's JSON."""
    seq = up["sequence"]["value"]
    genes = [g.get("geneName", {}).get("value") for g in up.get("genes", []) if g.get("geneName")]
    comments = up.get("comments", [])
    audit = up.get("entryAudit") or {}
    return {"acc": up.get("primaryAccession", acc), "id": up.get("uniProtkbId"), "reviewed": up.get("entryType"),
            "name": _recommended_name(up), "gene": genes[0] if genes else None,
            "organism": (up.get("organism") or {}).get("scientificName"), "length": len(seq), "sequence": seq,
            "function": next((t["value"] for c in comments if c.get("commentType") == "FUNCTION"
                              for t in c.get("texts", [])), None),
            "location": [l["location"]["value"] for c in comments if c.get("commentType") == "SUBCELLULAR LOCATION"
                         for l in c.get("subcellularLocations", []) if l.get("location")],
            "uniprot_version": audit.get("entryVersion"), "uniprot_updated": audit.get("lastAnnotationUpdateDate"),
            "release": release, "fetched_at": _today()}


def _uniprot_fact_fields(up):
    genes = [g.get("geneName", {}).get("value") for g in up.get("genes", []) if g.get("geneName")]
    audit = up.get("entryAudit") or {}
    return _rows([["Accession", up.get("primaryAccession")], ["Entry name", up.get("uniProtkbId")],
                  ["Status", up.get("entryType")], ["Protein name", _recommended_name(up)],
                  ["Gene", genes[0] if genes else None],
                  ["Organism", (up.get("organism") or {}).get("scientificName")],
                  ["Protein existence", up.get("proteinExistence")],
                  ["Entry version", audit.get("entryVersion")],
                  ["Annotation updated", audit.get("lastAnnotationUpdateDate")]])


def _uniprot_step(acc):
    """_uniprot_step_raw, but an answer of an unexpected shape is a not_reached Fact, not a crash."""
    try:
        fact, entry, why = _uniprot_step_raw(acc)
        if entry is not None:
            _parse_entry(entry, acc, None)
            _ids(entry)
        return fact, entry, why
    except SHAPE_ERRORS:
        why = "UniProt's answer for %s was of an unexpected response shape" % acc
        return _fact("UniProt", "not_reached", note=why), None, why


def _uniprot_step_raw(acc):
    """Ask UniProt about `acc` -> (Fact, live entry JSON or None, why there is no entry).

    A merged accession is answered by UniProt with a redirect to the surviving entry; that entry is
    used and the Fact's note says so. A demerged one has several targets: they are listed, none is picked.
    An isoform id is answered with its canonical entry, and the Fact's note says so.
    """
    asked, acc = acc, entry_accession(acc)

    def no(status, why, url=None, **kw):
        return _fact("UniProt", status, url, note=why, **kw), None, why
    if not UNIPROT_ACCESSION.fullmatch(acc):
        return no("not_applicable", "%r is not a UniProt accession, so UniProt was not asked" % acc)
    url = UNIPROT_BASE + acc + ".json"
    try:
        status, hdrs, body = _get(url, {"Accept": "application/json"})
    except (urllib.error.URLError, http.client.HTTPException, OSError) as e:
        return no("not_reached", "UniProt could not be reached for %s (%s)" % (acc, e), url)
    if status == 404:
        return no("not_found", "UniProt has no entry %s (HTTP 404)" % acc, url)
    if status != 200:
        return no("not_reached", "UniProt answered HTTP %s for %s" % (status, acc), url)
    release = hdrs.get("x-uniprot-release")
    try:
        up = json.loads(body.decode("utf-8"))
    except ValueError:
        return no("not_reached", "UniProt's answer for %s was not JSON" % acc, url, release=release)
    if not isinstance(up, dict):
        return no("not_reached", "UniProt's answer for %s was not an entry" % acc, url, release=release)
    if up.get("entryType") == "Inactive" or "sequence" not in up:
        ir = up.get("inactiveReason") or {}
        kind, targets = ir.get("inactiveReasonType", "no sequence"), ir.get("mergeDemergeTo") or []
        why = "UniProt entry %s is inactive (%s)" % (acc, kind)
        if kind == "DELETED":
            return no("not_found", why + (": %s" % ir["deletedReason"] if ir.get("deletedReason") else ""),
                      url, release=release)
        upi = (up.get("extraAttributes") or {}).get("uniParcId")
        fields = _rows([["Accession", up.get("primaryAccession")], ["Status", up.get("entryType")],
                        ["Reason", kind], ["Merged or demerged into", targets], ["UniParc", upi]])
        if kind == "DEMERGED" or len(targets) != 1:
            note = ("%s is inactive in UniProt: %s%s. It maps to more than one entry, and none is picked here; "
                    "ask for the accession you want" % (acc, kind, " into " + ", ".join(targets) if targets else ""))
        else:
            note = "%s is inactive in UniProt: %s into %s" % (acc, kind, targets[0])
        return _fact("UniProt", "ok", url, release, LICENCE["UniProt"], fields, note), None, why
    note = None
    if up.get("primaryAccession") and up["primaryAccession"] != acc:
        note = ("%s is not a current UniProt entry; UniProt answered with %s (merged), and everything here is "
                "%s's" % (acc, up["primaryAccession"], up["primaryAccession"]))
    if asked != acc:
        iso = ("%s is an isoform, which has no UniProt entry of its own; the canonical entry, %s, is shown"
               % (asked, up.get("primaryAccession") or acc))
        note = iso if note is None else note + ". " + iso
    return _fact("UniProt", "ok", url, release, LICENCE["UniProt"], _uniprot_fact_fields(up), note), up, None


def fetch_uniprot(acc, timeout=TIMEOUT, errors=None):
    """The fields the page's card shows from one UniProtKB entry, or None; the reason is appended to `errors`.

    `timeout` is kept for callers; _get always times out after TIMEOUT seconds.
    """
    return lookup(acc, errors=errors, _only_uniprot=True)[0]


def _xrefs(up, database):
    return [x for x in up.get("uniProtKBCrossReferences", []) if x.get("database") == database]


def _prop(xref, key):
    return next((p.get("value") for p in xref.get("properties", []) if p.get("key") == key), None)


def _ids(up):
    """The identifiers UniProt itself lists for this entry: GeneID, HGNC id, an Ensembl gene id, the HPA gene id."""
    def first(db):
        x = _xrefs(up, db)
        return x[0] if x else None
    gid, hgnc, hpa = first("GeneID"), first("HGNC"), first("HPA")
    ensembl = None
    for db in ("Ensembl", "EnsemblFungi"):
        for x in _xrefs(up, db):
            g = _prop(x, "GeneId")
            if g:
                ensembl = re.sub(r"\.\d+$", "", g)       # UniProt appends a version Ensembl's lookup may not resolve
                break
        if ensembl:
            break
    return {"taxid": (up.get("organism") or {}).get("taxonId"),
            "organism": (up.get("organism") or {}).get("scientificName"),
            "geneid": gid["id"] if gid else None, "hgnc": hgnc["id"] if hgnc else None,
            "ensembl": ensembl, "hpa": hpa["id"] if hpa else None}


# ---- NCBI -----------------------------------------------------------------------------------------

@_guarded("NCBI Gene")
def _ncbi_gene(geneid):
    src = "NCBI Gene"
    if not geneid or not re.fullmatch(r"\d{1,12}", str(geneid)):
        return _fact(src, "not_applicable", note="no NCBI GeneID is known for this protein, so NCBI Gene was not asked")
    url = _ncbi_url("esummary.fcgi", db="gene", id=geneid, retmode="json")
    data, _, fail = _call_json(src, url)
    if fail:
        return fail
    rec = ((data or {}).get("result") or {}).get(str(geneid))
    if not isinstance(rec, dict) or rec.get("error"):
        return _fact(src, "not_found", url, note="NCBI Gene has no record for GeneID %s" % geneid)
    org = rec.get("organism") or {}
    fields = _rows([["Symbol", rec.get("name")], ["Description", rec.get("description")],
                    ["Also known as", rec.get("otheraliases")], ["Organism", org.get("scientificname")],
                    ["Chromosome", rec.get("chromosome")], ["Map location", rec.get("maplocation")],
                    ["Summary", rec.get("summary")]])
    note = None
    if rec.get("currentid"):
        note = "GeneID %s was replaced by GeneID %s" % (geneid, rec["currentid"])
    return _fact(src, "ok", url, None, LICENCE["NCBI"], fields, note)


@_guarded("NCBI Protein", tuple_result=True)
def _ncbi_protein(acc):
    """-> (Fact, GeneID or None) for a protein accession NCBI knows (RefSeq, GenBank)."""
    src = "NCBI Protein"
    if not re.fullmatch(r"[A-Za-z0-9_.]{4,30}", acc):
        return _fact(src, "not_applicable", note="%r is not a protein accession NCBI could be asked about" % acc), None
    url = _ncbi_url("esummary.fcgi", db="protein", id=acc, retmode="json")
    data, _, fail = _call_json(src, url)
    if fail:
        return fail, None
    res = (data or {}).get("result") or {}
    uids = res.get("uids") or []
    rec = res.get(str(uids[0])) if uids else None
    if not isinstance(rec, dict) or rec.get("error"):
        return _fact(src, "not_found", url, note="NCBI Protein has no record for %s" % acc), None
    fields = _rows([["Title", rec.get("title")], ["Organism", rec.get("organism")],
                    ["Accession version", rec.get("accessionversion")], ["Length", rec.get("slen")],
                    ["Updated", rec.get("updatedate")], ["Database", rec.get("sourcedb")]])
    # the docsum names no gene; the GenBank flat record does (/gene= and /db_xref="GeneID:...")
    gp_url = _ncbi_url("efetch.fcgi", db="protein", id=acc, rettype="gp", retmode="text")
    geneid, note = None, None
    try:
        status, _, body = _get(gp_url)
        text = body.decode("utf-8", "replace") if status == 200 else ""
        if status != 200:
            note = "the gene is not named: NCBI answered HTTP %s for the GenBank record" % status
    except (urllib.error.URLError, http.client.HTTPException, OSError) as e:
        text, note = "", "the gene is not named: the GenBank record could not be fetched (%s)" % e
    g, i = re.search(r'/gene="([^"]+)"', text), re.search(r'/db_xref="GeneID:(\d+)"', text)
    if g:
        fields.append(["Gene", g.group(1)])
    if i:
        geneid = i.group(1)
        fields.append(["GeneID", geneid])
    return _fact(src, "ok", url, None, LICENCE["NCBI"], fields, note), geneid


# ---- HGNC, Ensembl, Human Protein Atlas -------------------------------------------------------------

def _not_human(src, ids, what):
    return _fact(src, "not_applicable",
                 note="%s covers human genes only; this protein is from %s" % (what, ids["organism"] or "another organism"))


@_guarded("HGNC")
def _hgnc(ids):
    src = "HGNC"
    if ids["taxid"] != HUMAN:
        return _not_human(src, ids, "HGNC")
    if not ids["hgnc"] or not re.fullmatch(r"HGNC:\d{1,8}", ids["hgnc"]):
        return _fact(src, "not_applicable", note="UniProt lists no HGNC id for this protein, so HGNC was not asked")
    acc = {"Accept": "application/json"}          # without it HGNC answers XML
    url = "https://rest.genenames.org/fetch/hgnc_id/" + ids["hgnc"]
    data, _, fail = _call_json(src, url, acc)
    if fail:
        return fail
    docs = ((data or {}).get("response") or {}).get("docs") or []
    if not docs:
        return _fact(src, "not_found", url, note="HGNC has no record for %s" % ids["hgnc"])
    d = docs[0]
    fields = _rows([["Symbol", d.get("symbol")], ["Name", d.get("name")], ["HGNC ID", d.get("hgnc_id")],
                    ["Status", d.get("status")], ["Locus group", d.get("locus_group")],
                    ["Locus type", d.get("locus_type")], ["Location", d.get("location")],
                    ["Alias symbols", d.get("alias_symbol")], ["Alias names", d.get("alias_name")],
                    ["Previous symbols", d.get("prev_symbol")], ["Previous names", d.get("prev_name")],
                    ["Gene groups", d.get("gene_group")]])
    release = _release(src, "https://rest.genenames.org/info", acc, lambda d: d["lastModified"])
    return _fact(src, "ok", url, release, LICENCE["HGNC"], fields)


@_guarded("Ensembl")
def _ensembl(ids):
    src = "Ensembl"
    eid = ids["ensembl"]
    if not eid or not re.fullmatch(r"[A-Za-z0-9_.-]{3,40}", eid):
        return _fact(src, "not_applicable", note="UniProt lists no Ensembl gene id for this protein, so Ensembl was not asked")
    url = "https://rest.ensembl.org/lookup/id/%s?content-type=application/json" % urllib.parse.quote(eid)
    data, _, fail = _call_json(src, url)
    if fail:
        return fail
    if not isinstance(data, dict) or not data.get("id"):
        return _fact(src, "not_found", url, note="Ensembl has no gene %s" % eid)
    fields = _rows([["Gene ID", data.get("id")], ["Name", data.get("display_name")],
                    ["Description", data.get("description")], ["Biotype", data.get("biotype")],
                    ["Species", data.get("species")], ["Assembly", data.get("assembly_name")],
                    ["Chromosome", data.get("seq_region_name")], ["Start", data.get("start")],
                    ["End", data.get("end")], ["Strand", data.get("strand")], ["Source", data.get("source")]])
    release = _release(src, "https://rest.ensembl.org/info/data?content-type=application/json", None,
                       lambda d: d["releases"][0])
    return _fact(src, "ok", url, release, LICENCE["Ensembl"], fields)


@_guarded("Human Protein Atlas")
def _hpa(ids):
    src = "Human Protein Atlas"
    if ids["taxid"] != HUMAN:
        return _not_human(src, ids, "The Human Protein Atlas")
    gene = ids["hpa"] or ids["ensembl"]
    if not gene or not re.fullmatch(r"ENSG\d{11}", gene):
        return _fact(src, "not_applicable", note="UniProt lists no Ensembl gene id (ENSG) for this protein, so the "
                                                 "Human Protein Atlas was not asked")
    url = "https://www.proteinatlas.org/%s.json" % gene
    data, _, fail = _call_json(src, url)
    if fail:
        return fail
    if not isinstance(data, dict) or not data:
        return _fact(src, "not_found", url, note="the Human Protein Atlas has no record for %s" % gene)
    # HPA states no version or date in this answer, so release stays None. A null here means "no HPA
    # value", not "absent from blood"; the types are the source's (a number in one record, a string in another).
    fields = [[k, data.get(k)] for k in HPA_FIELDS]
    return _fact(src, "ok", url, None, LICENCE["Human Protein Atlas"], fields)


# ---- the entry points ---------------------------------------------------------------------------------

def lookup(acc, errors=None, _only_uniprot=False):
    """-> (the UniProt fields for the page or None, [Fact, ...]); the reason UniProt gave none goes to `errors`."""
    ufact, entry, why = _uniprot_step(acc)
    facts = [ufact]
    if entry is None:
        if errors is not None and why:
            errors.append(why)
        if not _only_uniprot and ufact["status"] in ("not_found", "not_applicable"):
            # not in UniProt: ask NCBI about the accession itself, then about its gene
            prot, geneid = _ncbi_protein(acc)
            facts.append(prot)
            if geneid:
                facts.append(_ncbi_gene(geneid))
        return None, facts
    page = _parse_entry(entry, acc, ufact["release"])
    if _only_uniprot:
        return page, facts
    ids = _ids(entry)
    facts += [_ncbi_gene(ids["geneid"]), _hgnc(ids), _ensembl(ids), _hpa(ids)]
    return page, facts


def facts_for(accession):
    """The Facts for one accession, in order: UniProt, then (NCBI Protein for an accession UniProt lacks) NCBI Gene, HGNC, Ensembl, HPA."""
    return lookup(accession)[1]


if __name__ == "__main__":
    out = facts_for(sys.argv[1]) if len(sys.argv) > 1 else sys.exit("usage: protein_facts.py ACCESSION")
    json.dump(out, sys.stdout, indent=1)
    print()
