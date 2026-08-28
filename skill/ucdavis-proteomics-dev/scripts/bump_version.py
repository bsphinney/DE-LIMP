#!/usr/bin/env python3
"""
bump_version.py  --  Move the skill's version in every place that has to agree.

The version lives in THREE places, and CI fails the whole run if any one of them
disagrees with the others:

    skill/ucdavis-proteomics-core-pipeline/.claude-plugin/plugin.json   "version"
    .claude-plugin/marketplace.json                                    "version"  (top level)
    .claude-plugin/marketplace.json                                    plugins[].version

Two consecutive PRs bumped `plugin.json` alone and were caught by CI with
`version mismatch: plugin.json X vs marketplace ['Y']`. That is the check doing its job, but
it costs a push, a failed run, and a second commit every time. Nothing about the third
location is discoverable from the file you are editing, so this exists to make "bump the
version" a single action rather than something to remember.

    python3 scripts/bump_version.py 2.5.0        # set an exact version
    python3 scripts/bump_version.py patch        # 2.4.1 -> 2.4.2
    python3 scripts/bump_version.py minor        # 2.4.1 -> 2.5.0
    python3 scripts/bump_version.py major        # 2.4.1 -> 3.0.0
    python3 scripts/bump_version.py --check      # report agreement, change nothing

Edits are surgical -- the version LINE is rewritten by regex and nothing else in the file is
touched. Rewriting these files with json.dump() would reflow them and convert every escaped
`\\u2014` to a literal em-dash, producing a diff full of changes nobody made.
"""
import argparse
import json
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)                       # skill/ucdavis-proteomics-core-pipeline
REPO = os.path.dirname(os.path.dirname(SKILL))      # repo root
PLUGIN = os.path.join(SKILL, ".claude-plugin", "plugin.json")
MARKET = os.path.join(REPO, ".claude-plugin", "marketplace.json")

SEMVER = re.compile(r"^\d+\.\d+\.\d+$")


def current():
    """(plugin version, [marketplace entry versions]) -- read, never assumed to agree.

    Only the `plugins[]` entries for THIS plugin count. The marketplace's own top-level
    `version` is not this plugin's version: once the marketplace lists more than one plugin
    it cannot equal all of them, and CI does not require it to (it compares plugins[] only).
    Including it made `--check` fail for the dev sandbox purely for existing."""
    p = json.load(open(PLUGIN))
    m = json.load(open(MARKET))
    mv = [x.get("version") for x in m.get("plugins", []) if x.get("name") == p["name"]]
    return p["version"], mv


def bumped(v, part):
    major, minor, patch = (int(x) for x in v.split("."))
    if part == "major":
        return f"{major + 1}.0.0"
    if part == "minor":
        return f"{major}.{minor + 1}.0"
    return f"{major}.{minor}.{patch + 1}"


def rewrite(path, old, new):
    """Replace `"version": "<old>"` with the new one, leaving every other byte alone."""
    src = open(path).read()
    pat = re.compile(r'("version"\s*:\s*")' + re.escape(old) + r'(")')
    out, n = pat.subn(r"\g<1>" + new + r"\g<2>", src)
    if n:
        open(path, "w").write(out)
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("version", nargs="?",
                    help="an exact version (2.5.0) or one of: major, minor, patch")
    ap.add_argument("--check", action="store_true",
                    help="report whether the three locations agree; change nothing")
    a = ap.parse_args()

    pv, mv = current()
    if a.check or not a.version:
        agree = all(v == pv for v in mv)
        print(f"plugin.json: {pv} | marketplace: {mv} | agree: {agree}")
        # Exit non-zero on disagreement so this can be used as a pre-push check.
        sys.exit(0 if agree else 1)

    new = bumped(pv, a.version) if a.version in ("major", "minor", "patch") else a.version
    if not SEMVER.match(new):
        sys.exit(f"not a version: {new!r} (want X.Y.Z, or major/minor/patch)")

    total = rewrite(PLUGIN, pv, new)
    name = json.load(open(PLUGIN))["name"]
    src = open(MARKET).read()

    # Find THIS plugin's entry by its name, then the first "version" after it. A regex for
    # the whole `{...}` block does not work: the entries contain a nested `author` object,
    # so a [^{}]* block match stops at the inner brace and silently matches nothing --
    # which is exactly how the first version of this failed (it moved plugin.json and left
    # the marketplace behind, the very state this script exists to prevent).
    i = src.find(f'"name": "{name}"')
    if i != -1:
        m = re.compile(r'("version"\s*:\s*")([^"]+)(")').search(src, i)
        if m:
            src = src[:m.start()] + m.group(1) + new + m.group(3) + src[m.end():]
            total += 1
    # The marketplace's own top-level version tracks the PRIMARY plugin (the first entry);
    # with more than one plugin listed it cannot equal all of them.
    mp = json.load(open(MARKET))
    if mp.get("plugins") and mp["plugins"][0].get("name") == name:
        m0 = re.compile(r'("version"\s*:\s*")([^"]+)(")').search(src)
        if m0 and m0.start() < src.find('"plugins"'):
            src = src[:m0.start()] + m0.group(1) + new + m0.group(3) + src[m0.end():]
            total += 1
    open(MARKET, "w").write(src)

    pv2, mv2 = current()
    ok = pv2 == new and all(v == new for v in mv2)
    print(f"{pv} -> {new}: {total} line(s) updated in 2 files")
    print(f"plugin.json: {pv2} | marketplace: {mv2} | agree: {ok}")
    if not ok:
        sys.exit("version did not land in every location -- check both files by hand")


if __name__ == "__main__":
    main()
