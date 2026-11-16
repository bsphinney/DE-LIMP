#!/usr/bin/env python3
"""Record the live responses tests/test_protein_facts.py serves. Run by hand, rarely:

    python3 record.py            # needs the network; writes index.json and one body file per call

It runs protein_facts.facts_for() for each accession below with protein_facts._get wrapped, so the
URLs recorded are exactly the ones the module asks for. Calls to one host are at least 1 s apart
(NCBI also keeps protein_facts' own 0.34 s). No e-mail address and no API key are sent: NCBI_EMAIL
and NCBI_API_KEY are cleared first. A URL asked twice is fetched once.

Accessions: P80723 (BASP1), P02768 (ALB), P15108 (yeast HSC82), O43596 (merged into P80723), P04277
(demerged into P02769 and P02768), and the RefSeq protein NP_006308.3 (BASP1's; not a UniProt id).
"""
import datetime
import hashlib
import json
import os
import re
import sys
import time
import urllib.parse

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(HERE))), "scripts"))
import protein_facts as pf  # noqa: E402

ACCESSIONS = ["P80723", "P02768", "P15108", "O43596", "P04277", "NP_006308.3"]
KEEP_HEADERS = ("content-type", "x-uniprot-release", "x-uniprot-release-date", "x-final-url")
os.environ.pop("NCBI_EMAIL", None)
os.environ.pop("NCBI_API_KEY", None)

real_get = pf._get
responses, last = {}, {}


def slug(url):
    p = urllib.parse.urlparse(url)
    s = re.sub(r"[^A-Za-z0-9]+", "_", (p.hostname.split(".")[-2] + p.path + "_" + p.query)).strip("_")[:70]
    return "%s_%s" % (s, hashlib.sha1(url.encode()).hexdigest()[:6])


def strip_unread(host, body):
    """A UniProt entry without its `references` (literature citations, with their author lists):
    protein_facts never reads them, and the repository's identifier guard reads every text file
    (a citation's author can share a name with a client). Everything else is kept as served."""
    if host != "rest.uniprot.org" or body[:1] != b"{":
        return body
    d = json.loads(body)
    d.pop("references", None)
    return json.dumps(d, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def recording_get(url, headers=None, timeout=pf.TIMEOUT):
    if url in responses:
        return responses[url]["_t"]
    host = urllib.parse.urlparse(url).hostname
    wait = 1.0 - (time.monotonic() - last.get(host, -10))
    if wait > 0:
        time.sleep(wait)
    last[host] = time.monotonic()
    status, hdrs, body = real_get(url, headers, timeout)
    body = strip_unread(host, body)
    ext = ".txt" if (body[:1] not in (b"{", b"[")) else ".json"
    name = slug(url) + ext
    with open(os.path.join(HERE, name), "wb") as fh:
        fh.write(body)
    responses[url] = {"file": name, "status": status, "headers": {k: v for k, v in hdrs.items() if k in KEEP_HEADERS},
                      "recorded_at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
                      "_t": (status, hdrs, body)}
    print(status, len(body), url)
    return status, hdrs, body


pf._get = recording_get
for acc in ACCESSIONS:
    print("==", acc)
    for f in pf.facts_for(acc):
        print("  ", f["source"], f["status"], f["note"] or "")
index = {"recorded_by": "record.py", "accessions": ACCESSIONS,
         "responses": {u: {k: v for k, v in r.items() if k != "_t"} for u, r in responses.items()}}
with open(os.path.join(HERE, "index.json"), "w") as fh:
    json.dump(index, fh, indent=1, sort_keys=True)
