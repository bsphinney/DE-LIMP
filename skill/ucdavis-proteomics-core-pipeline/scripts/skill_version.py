#!/usr/bin/env python3
"""
skill_version.py  --  The ONE reader of this skill's version: .claude-plugin/plugin.json beside
scripts/. In the installed skill it is always there. On HIVE it is there when both folders were
put up together (`hive_exec.sh --put-skill`, again after every skill update).

Every record that names the skill version reads it here: record_run.py, provenance.py,
make_deposit.py (sdrf.tsv, prepare_upload.sbatch), core_submission.py. Two cannot import it and
carry a mirror, kept equal by tests/test_skill_version.py: notify_slack.py (on HIVE it also runs
from stdin, with no sibling to import) and report_issue.sh (bash; on Windows `python3` is often the
Microsoft Store stub).

When plugin.json is not there -- a scripts/ folder copied up on its own -- the version is UNKNOWN,
a tagged value, never a guessed one (CLAUDE.md rule 2): it used to read "0.0.0" in sdrf.tsv.

    python3 skill_version.py        # prints the version, or UNKNOWN
"""
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
UNKNOWN = "(unknown — plugin.json not found)"


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


def label(version):
    """How a version reads after the skill's name: "v2.8.0", or the UNKNOWN tag as it is."""
    return version if version == UNKNOWN else f"v{version}"


if __name__ == "__main__":
    print(skill_version())
