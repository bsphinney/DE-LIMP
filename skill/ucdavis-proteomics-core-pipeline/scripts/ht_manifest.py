"""
ht_manifest.py  --  Turn a Core high-throughput submission number into a search file list.

A high-throughput submission is not a folder. Submission 0793 can occupy two trays, and
the second tray's filenames may never mention 0793 at all -- the extent of a submission is
inferred from the acquisition counter, not stored as a directory. It also contains runs
that must NOT be searched as customer samples (well blanks, HeLa standards), and runs the
operator already knows came out badly. Globbing a directory gets all of this wrong,
silently, and you find out six hours later.

STAN knows which files a submission is. This skill knows how to search them. Neither should
learn the other's job, so this script does exactly one thing: it asks STAN, checks the
answer is safe to search, and writes a file list the rest of the flow consumes normally.
After that, nothing about the run is HT-specific -- step 2 onward is the ordinary flow.

**Runs ON HIVE**, like fran_deposit.py -- the orchestrator invokes it through
`hive_exec.sh`, because that is where STAN, its database credential and the raw files are.

    bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py fetch 0793 --out ~/ht0793'
    bash scripts/hive_exec.sh 'python3 ~/proteomics-pipeline/scripts/ht_manifest.py link 0793'

`fetch` writes `files.txt` (one absolute .d path per line, for --raw/--files) and
`ht_manifest.json` (the full STAN payload plus the gate results).

`link` only READS. There is deliberately no write-back: STAN v1.0.42 added an
`ht_searches` table and an `ht-record-search` command, and v1.0.43 reverted both, because
a table in STAN is a second copy of what FRAN already knows and the copy is the one that
goes stale. The submission number is in the raw filenames, so FRAN's
`/api/internal/submission/{id}` already answers "what was searched for 0793". `link` just
confirms that loop closed.

EXIT CODES -- the orchestrator must honour these:
    0  safe to proceed
    2  HARD GATE FAILED. Do not search. The run would silently cover a subset.
    3  STAN could not be reached, or is too old to have `ht-manifest`.
    4  (`link` only) not visible in FRAN yet. NOT a search failure.
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

# STAN's own venv on HIVE. Overridable because nothing should hardcode one person's path.
DEFAULT_STAN = "/quobyte/proteomics-grp/brett/stan_venv/bin/stan"

# Where to look for the PG Farm credential, in order. The point of the SHARED entry is that
# this is not a personal secret at all: it is a 7-day token minted from the
# `genome-proteomics-service-account` secret -- one facility identity that STAN and FRAN
# both already authenticate as. It only *looks* personal because it lives under one user's
# directory at mode 0600. A group-readable copy under proteomics-grp makes the CLI path
# work for the whole Core, gated by exactly the group membership fran_deposit.py already
# treats as "is this a Core search".
SHARED_TOKEN = "/quobyte/proteomics-grp/etc/pgfarm_token"
OWNER_TOKEN = "/quobyte/proteomics-grp/brett/.pgfarm_token"
TOKEN_CANDIDATES = (SHARED_TOKEN, "~/.pgfarm_token", OWNER_TOKEN)
DEFAULT_TOKEN = OWNER_TOKEN  # kept for the error message's "tried" list

# A tray is 96 wells. A submission reporting a handful of samples has probably been
# mistyped -- warn, don't block, because a genuine small submission is legal.
MIN_PLAUSIBLE_SAMPLES = 12
MAX_PLAUSIBLE_PLATES = 2


def _stan_bin(explicit: str | None) -> str:
    for cand in (explicit, os.environ.get("STAN_BIN"), DEFAULT_STAN, "stan"):
        if not cand:
            continue
        if cand == "stan" or os.path.exists(cand):
            return cand
    return DEFAULT_STAN


def _env(token_path: str | None) -> dict:
    """STAN's Postgres backend needs a credential. Fail loudly and specifically if the
    caller cannot read it -- the default token is mode 0600 and owned by one person, so
    'permission denied' is the EXPECTED failure for every other Core member, and it must
    not be mistaken for 'the submission does not exist'."""
    env = dict(os.environ)
    env["STAN_DB_BACKEND"] = env.get("STAN_DB_BACKEND", "pg")
    if env.get("PGPASSWORD"):
        return env

    explicit = token_path or env.get("STAN_PG_TOKEN")
    candidates = [explicit] if explicit else list(TOKEN_CANDIDATES)

    tried = []
    for cand in candidates:
        p = os.path.expanduser(cand)
        try:
            tok = open(p).read().strip()
        except PermissionError:
            tried.append(f"{p} — permission denied (owned by another user, mode 0600)")
            continue
        except (FileNotFoundError, NotADirectoryError):
            tried.append(f"{p} — not found")
            continue
        if tok:
            env["PGPASSWORD"] = tok
            return env
        tried.append(f"{p} — empty")

    detail = "\n".join(f"    {t}" for t in tried)
    sys.exit(
        f"[ht_manifest] no usable PG Farm credential. Tried:\n{detail}\n"
        f"\n"
        f"  This is NOT a personal secret — it is a 7-day token minted from the\n"
        f"  genome-proteomics-service-account, the identity STAN and FRAN both use. It\n"
        f"  reads as personal only because it lives under one user's directory at 0600.\n"
        f"\n"
        f"  Fixes, best first:\n"
        f"    1. Have the Core publish a group-readable copy at\n"
        f"         {SHARED_TOKEN}   (mode 0640, group proteomics-grp)\n"
        f"       and this script finds it with no configuration at all.\n"
        f"       NOTE: a plain chmod will NOT hold — pgfarm_refresh_token.py rewrites the\n"
        f"       file every 5 minutes via os.chmod(tmp, 0o600) + os.replace(), so the mode\n"
        f"       must be set BY the refresher, not after it.\n"
        f"    2. export STAN_PG_TOKEN=<a token file you can read>\n"
        f"    3. Skip the database entirely and use the hosted dashboard:\n"
        f"         --http https://ucd.stan-proteomics.org --share-token <tok>\n"
        f"  Do NOT copy the token into your home directory — it expires in 7 days and\n"
        f"  yours will not be refreshed.")


def _run_stan(argv: list, env: dict) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(argv, capture_output=True, text=True, env=env, timeout=180)
    except FileNotFoundError:
        sys.exit(f"[ht_manifest] STAN not found at {argv[0]}. Pass --stan or set STAN_BIN.")
    except subprocess.TimeoutExpired:
        sys.exit(f"[ht_manifest] STAN timed out after 180s running: {' '.join(argv[:3])}")


def _manifest_over_http(a) -> dict | int:
    """Ask the hosted dashboard instead of the local database.

    This is the path that makes the workflow usable by the whole Core rather than by
    whoever owns the Postgres token. `/api/ht/manifest` returns the same payload under the
    same access rules, and it accepts a per-submission **share token** — which is the only
    one of its three auth paths a headless HIVE job can actually use. Entra sign-in needs a
    browser; a share link does not.
    """
    import urllib.error
    import urllib.parse
    import urllib.request

    base = (a.http or os.environ.get("STAN_HT_URL") or "").rstrip("/")
    tok = a.share_token or os.environ.get("STAN_HT_SHARE_TOKEN")
    qs = {"q": a.submission, "include": a.include}
    if tok:
        qs["token"] = tok
    url = f"{base}/api/ht/manifest?" + urllib.parse.urlencode(qs)
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    if a.cookie:
        req.add_header("Cookie", a.cookie)
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        body = ""
        try:
            body = e.read().decode()[:400]
        except Exception:
            pass
        if e.code in (401, 403):
            print(f"[ht_manifest] {base} refused the request ({e.code}).\n"
                  f"  HT data is not public — it carries customer sample names and paths.\n"
                  f"  Use ONE of:\n"
                  f"    --share-token <tok>   a per-submission share link from the HT tab\n"
                  f"                          (the only option that works headless)\n"
                  f"    --cookie '<session>'  an Entra session cookie from a signed-in browser\n"
                  f"  Sign-in is Microsoft Entra ({base}/.auth/login/aad), not CAS.\n"
                  f"  server said: {body}", file=sys.stderr)
        else:
            print(f"[ht_manifest] {url} returned HTTP {e.code}: {body}", file=sys.stderr)
        return 3
    except urllib.error.URLError as e:
        print(f"[ht_manifest] cannot reach {base}: {e.reason}", file=sys.stderr)
        return 3
    except json.JSONDecodeError:
        print(f"[ht_manifest] {base} returned non-JSON", file=sys.stderr)
        return 3


def _manifest_over_cli(a) -> dict | int:
    env = _env(a.token)
    stan = _stan_bin(a.stan)
    p = _run_stan([stan, "ht-manifest", a.submission, "--include", a.include], env)
    if p.returncode != 0:
        err = (p.stderr or p.stdout or "").strip()
        if "No such command" in err:
            print(f"[ht_manifest] the STAN on this machine has no 'ht-manifest' command.\n"
                  f"  {stan}\n"
                  f"  HT support landed in STAN v1.0.43; this venv predates it. Update the\n"
                  f"  venv (pip install -U from the STAN repo) and re-run. Nothing was searched.",
                  file=sys.stderr)
            return 3
        print(f"[ht_manifest] stan ht-manifest failed (rc={p.returncode}):\n{err}", file=sys.stderr)
        return 3
    try:
        return json.loads(p.stdout)
    except json.JSONDecodeError:
        print(f"[ht_manifest] STAN returned non-JSON:\n{p.stdout[:500]}", file=sys.stderr)
        return 3


def fetch(a) -> int:
    use_http = bool(a.http or os.environ.get("STAN_HT_URL"))
    m = _manifest_over_http(a) if use_http else _manifest_over_cli(a)
    if isinstance(m, int):
        return m
    m.setdefault("source", "http" if use_http else "cli")

    files = m.get("files") or []
    gates, hard_fail = [], False

    # --- HARD: runs STAN knows about but has no raw path for. They are EXCLUDED from
    # `files`, so proceeding means searching a subset and then reporting success.
    missing = m.get("missing_paths") or []
    if missing:
        hard_fail = True
        gates.append({"gate": "missing_paths", "status": "FAIL", "n": len(missing),
                      "detail": f"{len(missing)} run(s) have no raw path and are excluded "
                                f"from the file list", "examples": missing[:5]})
    else:
        gates.append({"gate": "missing_paths", "status": "PASS", "n": 0})

    # --- HARD: nothing to search.
    if not files:
        hard_fail = True
        gates.append({"gate": "n_files", "status": "FAIL", "n": 0,
                      "detail": f"no files for submission {a.submission} with --include {a.include}"})
    else:
        gates.append({"gate": "n_files", "status": "PASS", "n": len(files)})

    # --- HARD: a path STAN resolved but the filesystem does not have. Cheap to check here;
    # otherwise a 120-file SLURM array dies partway through, hours in.
    unreadable = [f for f in files if not os.path.exists(f)]
    if unreadable:
        hard_fail = True
        gates.append({"gate": "paths_exist", "status": "FAIL", "n": len(unreadable),
                      "detail": f"{len(unreadable)} listed path(s) do not exist on this "
                                f"filesystem", "examples": unreadable[:5]})
    else:
        gates.append({"gate": "paths_exist", "status": "PASS", "n": len(files)})

    # --- WARN: shape checks. A submission fills one or two trays of 96.
    plates = m.get("plates") or []
    if len(plates) > MAX_PLAUSIBLE_PLATES:
        gates.append({"gate": "plates", "status": "WARN", "n": len(plates),
                      "detail": f"{len(plates)} plates ({', '.join(map(str, plates))}) — more "
                                f"than a submission usually spans. Confirm the extent with the "
                                f"operator; it is inferred from the acquisition counter."})
    else:
        gates.append({"gate": "plates", "status": "PASS", "n": len(plates), "plates": plates})

    counts = m.get("counts") or {}
    n_sample = counts.get("sample", len(files))
    if n_sample < MIN_PLAUSIBLE_SAMPLES:
        gates.append({"gate": "counts", "status": "WARN", "n": n_sample,
                      "detail": f"only {n_sample} customer samples — a mistyped submission "
                                f"number looks exactly like this. Confirm before searching."})
    else:
        gates.append({"gate": "counts", "status": "PASS", "n": n_sample, "counts": counts})

    n_rerun = m.get("n_needs_rerun") or 0
    if n_rerun:
        gates.append({"gate": "needs_rerun", "status": "INFO", "n": n_rerun,
                      "detail": f"{n_rerun} sample(s) flagged by STAN as wanting another "
                                f"injection. They ARE included here (--include samples). Use "
                                f"--include rerun to search only those."})

    m["gates"] = gates
    m["hard_fail"] = hard_fail

    outdir = a.out or "."
    os.makedirs(outdir, exist_ok=True)
    fl = os.path.join(outdir, "files.txt")
    mf = os.path.join(outdir, "ht_manifest.json")
    with open(fl, "w") as fh:
        fh.write("".join(f + "\n" for f in files))
    with open(mf, "w") as fh:
        json.dump(m, fh, indent=2)

    print(f"submission {m.get('submission', a.submission)}  include={m.get('include', a.include)}")
    print(f"  files      : {len(files)}")
    print(f"  plates     : {', '.join(map(str, plates)) or '(none reported)'}")
    print(f"  counts     : {counts}")
    print(f"  files.txt  : {fl}")
    print(f"  manifest   : {mf}")
    for g in gates:
        if g["status"] != "PASS":
            print(f"  [{g['status']}] {g['gate']}: {g.get('detail', '')}")
            for ex in g.get("examples", [])[:5]:
                print(f"        {ex}")
    if hard_fail:
        print("\nHARD GATE FAILED — do not search. Resolve the above with the operator first.",
              file=sys.stderr)
        return 2
    return 0


def link(a) -> int:
    """Confirm the finished search is visible under its submission in FRAN.

    **There is no write-back, by design.** STAN v1.0.42 briefly had an `ht_searches` table
    and an `ht-record-search` command; v1.0.43 reverted both. The reasoning is worth
    keeping, because it is the reason this function only reads: a table in STAN would be a
    second copy of what FRAN already knows, and the copy is the one that goes stale --- it
    depends on the skill remembering to call it, and says nothing when that step is
    skipped. It would also mean re-implementing FRAN's authorization.

    So the link is implicit: the submission number is in the raw filenames, so FRAN's
    `/api/internal/submission/{id}` already answers "what was searched for 0793" with the
    search id, engine, organism and counts. Nothing needs recording. This just *checks*
    the loop closed, so a plate that silently failed to ingest is visible.
    """
    import urllib.error
    import urllib.request

    base = (a.fran or os.environ.get("FRAN_URL") or "https://fran.stan-proteomics.org").rstrip("/")
    url = f"{base}/api/internal/submission/{urllib_quote(a.submission)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    if a.cookie:
        req.add_header("Cookie", a.cookie)
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            payload = json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        if e.code == 404:
            print(f"[ht_manifest] FRAN returned 404 for submission {a.submission}.\n"
                  f"  That endpoint is INTERNAL — 404 is also what it returns to a caller\n"
                  f"  who is not signed in, or to a lab user asking about another lab's\n"
                  f"  submission. So this is 'not visible to you', not proof of 'not ingested'.\n"
                  f"  Sign in at {base} (Microsoft Entra) or pass --cookie.", file=sys.stderr)
            return 4
        print(f"[ht_manifest] FRAN returned HTTP {e.code} for {url}", file=sys.stderr)
        return 4
    except urllib.error.URLError as e:
        print(f"[ht_manifest] cannot reach {base}: {e.reason}", file=sys.stderr)
        return 4

    body = payload.get("data", payload) if isinstance(payload, dict) else {}
    searches = (body or {}).get("searches") or []
    print(f"submission {a.submission} — FRAN knows {len(searches)} search(es)")
    for s in searches:
        print(f"  {str(s.get('search_name') or s.get('id'))[:44]:<46}"
              f"{str(s.get('search_engine') or '-'):<13}"
              f"prec={s.get('n_precursors_total')} pg={s.get('n_protein_groups_total')}")
    if not searches:
        print("\n  No search under this submission yet. If fran_deposit.py reported\n"
              "  staged_pending_cron, that is expected — FRAN ingests on its next scan.\n"
              "  Re-run this later; it is not a failure.", file=sys.stderr)
        return 4
    return 0


def urllib_quote(s: str) -> str:
    import urllib.parse
    return urllib.parse.quote(str(s), safe="")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--stan", default=None, help=f"path to the stan binary (default {DEFAULT_STAN})")
    ap.add_argument("--token", default=None,
                    help=f"file holding the STAN Postgres credential (default {DEFAULT_TOKEN}; "
                         f"or set STAN_PG_TOKEN / PGPASSWORD)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    f = sub.add_parser("fetch", help="get + validate the file list for a submission")
    f.add_argument("submission", help="submission number, e.g. 0793 (793 also works)")
    f.add_argument("--include", default="samples",
                   choices=["samples", "rerun", "standards", "all"],
                   help="samples (default) excludes blanks and HeLa standards")
    f.add_argument("--out", default=None, help="directory for files.txt + ht_manifest.json")
    f.add_argument("--http", default=None, metavar="BASE_URL",
                   help="ask the hosted dashboard instead of the local DB, e.g. "
                        "https://ucd.stan-proteomics.org (or set STAN_HT_URL). Use this when "
                        "you cannot read the Postgres token — which is every Core member "
                        "except its owner.")
    f.add_argument("--share-token", default=None,
                   help="per-submission share token from the HT tab (or STAN_HT_SHARE_TOKEN). "
                        "The only auth path that works headless.")
    f.add_argument("--cookie", default=None,
                   help="Entra session cookie from a signed-in browser, if you have one")
    f.set_defaults(func=fetch)

    r = sub.add_parser("link", help="check the finished search is visible under its "
                                    "submission in FRAN (read-only; there is no write-back)")
    r.add_argument("submission")
    r.add_argument("--fran", default=None, metavar="BASE_URL",
                   help="FRAN base URL (or FRAN_URL; default https://fran.stan-proteomics.org)")
    r.add_argument("--cookie", default=None,
                   help="Entra session cookie — the submission endpoint is internal-only")
    r.set_defaults(func=link)

    a = ap.parse_args()
    return a.func(a)


if __name__ == "__main__":
    sys.exit(main())
