#!/usr/bin/env python3
"""
share_map.py  --  Where a path lives on HIVE, and how a collaborator browses to a HIVE path from
Windows or a Mac. For the session README and AGENTS.md ("where is this on HIVE, so we can find it
again").

Which share is which HIVE path is hive_shares.tsv -- the one table, also read by hive_path.sh. A
path is resolved, first match wins:
  1. already a HIVE path (/quobyte/..., /nfs/lssc0/...)                      -> itself
  2. running ON HIVE (its Quobyte storage is mounted here) and the path exists -> its real path
  3. a UNC path (\\\\server\\share\\... or //server/share/...) of a share in the table
  4. under a Mac mount point in the table (/Volumes/proteomics/...)
  5. anything else -- a drive letter, a mount under another name: hive_path.sh --no-verify says
     which share it is (net use / mount, no ssh); the table then says where HIVE mounts it
A share that is not in the table is "not recorded" -- never hive_path.sh's unverified guess.

  python3 share_map.py <path>     # prints {"path", "hive", "how", "windows", "mac", "access"}
"""
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TABLE = os.path.join(HERE, "hive_shares.tsv")
HIVE_ROOTS = ("/quobyte/", "/nfs/lssc0/")
# The Core's Flinders share, by its name in the table, and the trees inside it the skill works in.
# ONE place: core_submission.py (raw files, service folders, Bioshare shares) and fran_deposit.py
# (which searches are the Core's, where the backfill walks) read them here.
FLINDERS_SHARE = "proteomics"
FLINDERS_RAW = ("Data", "raw_data")              # instrument raw files
FLINDERS_SERVICE = ("Data", "lab", "service")    # the Core's service folders


def load_table(path=TABLE):
    """[{server, share, hive, mac, windows, access}] -- `-` becomes None (not recorded; for
    `access`: nothing to add to the paths)."""
    keys = ("server", "share", "hive", "mac", "windows", "access")
    rows = []
    with open(path, encoding="utf-8") as fh:
        for ln in fh:
            if not ln.strip() or ln.startswith("#"):
                continue
            f = (ln.rstrip("\r\n").split("\t") + [""] * len(keys))[:len(keys)]
            rows.append({k: (v if v and v != "-" else None) for k, v in zip(keys, f)})
    return rows


def hive_root(share, rows=None):
    """The HIVE path of `share` in the table, or None when it records none."""
    rows = load_table() if rows is None else rows
    return next((r["hive"] for r in rows
                 if (r["share"] or "").lower() == share.lower() and r["hive"]), None)


def _match(server, share, rows):
    for r in rows:
        if (r["share"] or "").lower() == (share or "").lower() and \
                (r["server"] == "*" or (r["server"] or "").lower() == (server or "").lower()):
            return r
    return None


def _join(base, rest, sep="/"):
    rest = (rest or "").replace("\\", "/").strip("/")
    return base.rstrip("/\\") + (sep + rest.replace("/", sep) if rest else "")


def _unc(path):
    """(server, share, rest) of \\\\server\\share\\rest or //[user@]server/share/rest, else None."""
    p = path.replace("\\", "/")
    m = re.match(r"^//([^/]+)/([^/]+)/?(.*)$", p)
    if not m:
        return None
    return m.group(1).split("@")[-1], m.group(2), m.group(3)


def _ask_hive_path_sh(path, timeout=20):
    """(server, share, rest) as hive_path.sh --no-verify resolves it (drive letters via net use,
    other mounts via mount) -- or None. No ssh is made."""
    bash = shutil.which("bash")
    script = os.path.join(HERE, "hive_path.sh")
    if not bash or not os.path.isfile(script):
        return None
    try:
        r = subprocess.run([bash, script, "--no-verify", path], capture_output=True, text=True,
                           timeout=timeout)
        j = json.loads(r.stdout)
    except (OSError, ValueError, subprocess.TimeoutExpired):
        return None
    if not j.get("share"):
        return None
    return j.get("server") or "", j["share"], j.get("rest") or ""


def on_hive():
    return os.path.isdir("/quobyte")


def hive_of(path, rows=None, resolver=_ask_hive_path_sh, here_is_hive=None):
    """{"hive": <HIVE path or None>, "how": <how it was found, or why not>}."""
    if not path:
        return {"hive": None, "how": "not recorded"}
    rows = load_table() if rows is None else rows
    here_is_hive = on_hive() if here_is_hive is None else here_is_hive
    p = str(path)
    if p.startswith(HIVE_ROOTS) or p.rstrip("/") in ("/quobyte", "/nfs/lssc0"):
        return {"hive": p.rstrip("/") or "/", "how": "recorded as a HIVE path"}
    if here_is_hive and os.path.isabs(p) and os.path.exists(p):
        return {"hive": os.path.realpath(p), "how": "on HIVE (written here)"}
    u = _unc(p)
    if u:
        r = _match(u[0], u[1], rows)
        return ({"hive": _join(r["hive"], u[2]), "how": f"the \\\\{u[0]}\\{u[1]} share"} if r else
                {"hive": None, "how": f"not recorded: \\\\{u[0]}\\{u[1]} is not a share "
                                      f"hive_shares.tsv knows"})
    for r in rows:
        m = r["mac"]
        if m and (p == m or p.startswith(m.rstrip("/") + "/")):
            return {"hive": _join(r["hive"], p[len(m):]), "how": f"the share mounted at {m}"}
    got = resolver(p) if resolver else None
    if got:
        r = _match(got[0], got[1], rows)
        if r:
            return {"hive": _join(r["hive"], got[2]),
                    "how": f"the \\\\{got[0]}\\{got[1]} share (hive_path.sh)"}
        return {"hive": None, "how": f"not recorded: \\\\{got[0]}\\{got[1]} is not a share "
                                     f"hive_shares.tsv knows"}
    return {"hive": None, "how": "not recorded: not on a network share HIVE mounts"}


def views_of(hive, rows=None):
    """How a collaborator browses to a HIVE path: {"windows": UNC path or None, "mac": path or
    None, "access": the table's note on who can open that share, or None}. None where the table
    does not record that share's Windows/Mac name."""
    out = {"windows": None, "mac": None, "access": None}
    if not hive:
        return out
    rows = load_table() if rows is None else rows
    for r in rows:
        base = (r["hive"] or "").rstrip("/")
        if base and (hive == base or hive.startswith(base + "/")):
            rest = hive[len(base):]
            if r["windows"]:
                out["windows"] = _join(r["windows"], rest, sep="\\")
            if r["mac"]:
                out["mac"] = _join(r["mac"], rest)
            out["access"] = r.get("access")
            break
    return out


def locate(path, rows=None, resolver=_ask_hive_path_sh, here_is_hive=None):
    """{"path", "hive", "how", "windows", "mac", "access"} for one path."""
    rows = load_table() if rows is None else rows
    h = hive_of(path, rows, resolver, here_is_hive)
    return dict({"path": path}, **h, **views_of(h["hive"], rows))


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit("usage: share_map.py <path>")
    print(json.dumps(locate(sys.argv[1]), indent=2))
