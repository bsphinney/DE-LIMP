#!/usr/bin/env python3
"""
fetch_fasta.py  --  Resolve the search FASTA for a workflow's fasta spec.

Priority (cheapest / most-trusted first):
  1. fasta.path override  -> use it verbatim (e.g. a pre-staged HIVE proteome).
  2. UC Davis HIVE        -> reuse /quobyte/proteomics-grp/MRS/ proteomes +
                             contaminants instead of downloading (PLAN.md §7e).
  3. UniProt download     -> reference proteome FASTA for fasta.uniprot_proteome,
                             then append the universal contaminants if requested.

The contaminants set is the cRAP/universal contaminants FASTA. On HIVE we reuse
the copy the core already keeps; off HIVE we download the maxquant/cRAP set.

Usage:
  python3 fetch_fasta.py --proteome UP000005640 --add-contaminants \
      --out ./search.fasta [--hive] [--path /abs/override.fasta]

Emits JSON: {fasta, source, proteome, n_contaminants_appended}.
"""
import sys, os, json, argparse, glob, gzip, shutil, urllib.request, urllib.error

HIVE_MRS = "/quobyte/proteomics-grp/MRS"
# UniProt REST: stream a reference proteome as uncompressed FASTA.
UNIPROT_URL = ("https://rest.uniprot.org/uniprotkb/stream"
               "?format=fasta&compressed=false"
               "&query=%28proteome%3A{proteome}%29")
# cRAP universal contaminants (GPM). Small, stable, widely used.
CRAP_URL = "https://www.thegpm.org/crap/crap.fasta"


def _download(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": "proteomics-pipeline-skill"})
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "wb") as fh:
        shutil.copyfileobj(r, fh)


def _read_fasta_text(path):
    opn = gzip.open if path.endswith(".gz") else open
    with opn(path, "rt") as fh:
        return fh.read()


def _count_entries(text):
    return text.count(">")


def hive_proteome(proteome):
    """Look for a pre-staged proteome FASTA on HIVE by UniProt ID."""
    if not os.path.isdir(HIVE_MRS):
        return None
    for pat in (f"*{proteome}*.fasta", f"*{proteome}*.fa", f"*{proteome}*.fasta.gz"):
        hits = sorted(glob.glob(os.path.join(HIVE_MRS, "**", pat), recursive=True))
        if hits:
            return hits[0]
    return None


def hive_contaminants():
    if not os.path.isdir(HIVE_MRS):
        return None
    for pat in ("*ontaminant*.fasta", "*cRAP*.fasta", "*crap*.fasta", "*ontaminant*.fa"):
        hits = sorted(glob.glob(os.path.join(HIVE_MRS, "**", pat), recursive=True))
        if hits:
            return hits[0]
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--proteome", help="UniProt proteome ID, e.g. UP000005640")
    ap.add_argument("--path", help="explicit FASTA override; used verbatim if set")
    ap.add_argument("--add-contaminants", action="store_true")
    ap.add_argument("--out", required=True)
    ap.add_argument("--hive", action="store_true",
                    help="prefer pre-staged HIVE FASTAs (set when env is uc_davis_hive)")
    a = ap.parse_args()

    source = None
    base_text = None

    # 1. explicit override
    if a.path:
        if not os.path.exists(a.path):
            sys.exit(f"--path given but not found: {a.path}")
        base_text = _read_fasta_text(a.path)
        source = f"override:{a.path}"

    # 2. HIVE pre-staged
    if base_text is None and a.hive and a.proteome:
        staged = hive_proteome(a.proteome)
        if staged:
            base_text = _read_fasta_text(staged)
            source = f"hive:{staged}"

    # 3. UniProt download
    if base_text is None:
        if not a.proteome:
            sys.exit("Need --proteome (or --path). The workflow's fasta.uniprot_proteome.")
        tmp = a.out + ".proteome.tmp"
        url = UNIPROT_URL.format(proteome=a.proteome)
        try:
            _download(url, tmp)
        except (urllib.error.URLError, urllib.error.HTTPError) as e:
            sys.exit(f"UniProt download failed for {a.proteome}: {e}\n  URL: {url}")
        base_text = _read_fasta_text(tmp)
        os.remove(tmp)
        source = f"uniprot:{a.proteome}"

    n_base = _count_entries(base_text)
    if n_base == 0:
        sys.exit(f"Resolved FASTA has 0 sequences (source={source}). Refusing to proceed.")

    # contaminants
    n_contam = 0
    contam_text = ""
    if a.add_contaminants:
        cpath = hive_contaminants() if a.hive else None
        if cpath:
            contam_text = _read_fasta_text(cpath)
        else:
            ctmp = a.out + ".crap.tmp"
            try:
                _download(CRAP_URL, ctmp)
                contam_text = _read_fasta_text(ctmp)
                os.remove(ctmp)
            except (urllib.error.URLError, urllib.error.HTTPError) as e:
                sys.stderr.write(f"[fetch_fasta] WARNING: could not fetch contaminants ({e}); "
                                 "proceeding WITHOUT them\n")
        n_contam = _count_entries(contam_text)

    with open(a.out, "w") as fh:
        fh.write(base_text)
        if not base_text.endswith("\n"):
            fh.write("\n")
        if contam_text:
            fh.write(contam_text)

    print(json.dumps({
        "fasta": os.path.abspath(a.out),
        "source": source,
        "proteome": a.proteome,
        "n_sequences": n_base + n_contam,
        "n_proteome": n_base,
        "n_contaminants_appended": n_contam,
    }, indent=2))


if __name__ == "__main__":
    main()
