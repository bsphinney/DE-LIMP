#!/usr/bin/env python3
"""
skill_version.py  --  The ONE reader of this skill's version: .claude-plugin/plugin.json beside
scripts/. In the installed skill it is always there. On HIVE it is there when both folders were
put up together (`hive_exec.sh --put-skill`, again after every skill update).

Every record that names the skill version reads it here: record_run.py, provenance.py,
make_deposit.py (sdrf.tsv, prepare_upload.sbatch), core_submission.py. Three cannot import it and
carry a mirror, kept equal by tests/test_skill_version.py: notify_slack.py (on HIVE it also runs
from stdin, with no sibling to import), report_issue.sh (bash; on Windows `python3` is often the
Microsoft Store stub) and skill_version.R (R: run_de.R's, for the DE-LIMP session).

When plugin.json is not there -- a scripts/ folder copied up on its own -- the version is UNKNOWN,
a tagged value, never a guessed one (CLAUDE.md rule 2): it used to read "0.0.0" in sdrf.tsv.

It is also the Python reader of core_admins.txt beside it, the Proteomics Core admins:
core_admins(), line for line what skill_version.sh's skill_core_admins prints (bash). notes.py
imports it; notes.py's copy sent over stdin to HIVE carries a copy, kept equal by
tests/test_notes.py.

    python3 skill_version.py        # prints the version, or UNKNOWN
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
UNKNOWN = "(unknown — plugin.json not found)"
# The same tag in ASCII, for a file whose validators may take nothing else (sdrf.tsv)
UNKNOWN_ASCII = UNKNOWN.strip("()").replace(" — ", " - ")


def skill_version(here=HERE):
    """The version in <here>/../.claude-plugin/plugin.json, or UNKNOWN when it cannot be read."""
    try:
        with open(os.path.join(here, "..", ".claude-plugin", "plugin.json"),
                  encoding="utf-8") as fh:
            v = json.load(fh).get("version")
    except (OSError, ValueError, AttributeError, TypeError):
        v = None
    return v.strip() if isinstance(v, str) and v.strip() else UNKNOWN


def plugin_meta(here=HERE):
    """The whole of plugin.json ({} when it cannot be read), for its other fields (name,
    repository). The version is skill_version()'s: it says when it is unknown."""
    try:
        with open(os.path.join(here, "..", ".claude-plugin", "plugin.json"),
                  encoding="utf-8") as fh:
            m = json.load(fh)
    except (OSError, ValueError, TypeError):
        return {}
    return m if isinstance(m, dict) else {}


#: The characters sed's [[:space:]] matches (skill_version.sh strips them from both ends).
_SPACE = " \t\n\r\x0b\x0c"


def admin_lines(text):
    """core_admins.txt's text -> its admins, in file order: `#` to the end of a line dropped,
    blanks at both ends and every CR removed, empty lines skipped -- skill_version.sh's
    `sed 's/#.*//; s/[[:space:]]*$//; s/^[[:space:]]*//' | tr -d '\r' | grep -v '^$'`."""
    out = []
    for ln in text.split("\n"):
        ln = ln.split("#", 1)[0].rstrip(_SPACE).lstrip(_SPACE).replace("\r", "")
        if ln:
            out.append(ln)
    return out


def core_admins(path=None):
    """The Proteomics Core admins in core_admins.txt (`path`, default beside this file). [] when
    it is missing or unreadable: no admins, and whatever they would vouch for is not trusted."""
    try:
        with open(path or os.path.join(HERE, "core_admins.txt"), "rb") as fh:
            data = fh.read()
    except OSError:
        return []
    return admin_lines(data.decode("utf-8", "replace"))


def label(version, ascii=False):
    """How a version reads after the skill's name: "v2.8.0", or the UNKNOWN tag as it is --
    UNKNOWN_ASCII with ascii=True."""
    if version == UNKNOWN:
        return UNKNOWN_ASCII if ascii else UNKNOWN
    return f"v{version}"


if __name__ == "__main__":
    print(skill_version())
