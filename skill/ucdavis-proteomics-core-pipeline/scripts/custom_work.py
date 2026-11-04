#!/usr/bin/env python3
"""
custom_work.py  --  work done outside the skill's steps, folded into a session by
`session.py finalize` (2.14.2, Amendment 2 to Task 13).

The beer re-analysis wrote its species-composition code and its report in a SECOND folder on the
laptop, outside the HIVE session, and delivered by hand: the session's record never saw them.
finalize now takes them in:

    session.py finalize --dir <S> --report <html> --methods <md> \
        --code-dir <dir> [--code-dir <dir> ...] --extra-dir <dir> [--extra-dir <dir> ...]

Each is copied to output/custom/<kind>/<name> (<name>: the file's or folder's own name), and
output/custom/custom.json records where it came from and every file's sha256. finalize then lists
each item in MANIFEST.txt, AGENTS.md and the README, and puts the checksums in the bundle's
run_manifest.json (`custom`). The rules:
  * a report the skill did not write (report_files.skill_report_state) goes in only under
    --override-finish "<reason>" (finalize's gate), and the README calls it custom;
  * refused before anything is written: a missing source, a file given as a folder (or the
    reverse), a source inside output/custom/ itself, a source that holds the session, a folder
    that is a home folder or holds one (C2 security H1), a folder or file inside the logs/, input/
    or output/ of this session or of any other (session.is_session_dir; C2 security M2: logs/ holds
    the Core-internal conversation and decisions the zip leaves out), a file that cannot be read
    (C2 errors_ux M4), a --report/--methods file in which a secret is found, two sources of one
    kind with the same name, and more than --custom-max-gb (default MAX_GB) in all;
  * a symlink inside a folder is neither followed nor copied: MANIFEST.txt lists it;
  * secrets, keys, credentials, environments, version control and caches stay out of a folder's
    copy (left_out): the folders in LEFT_OUT_DIRS (.git, .venv, node_modules, renv/library, and
    .ssh, .aws, .config, .gnupg, .kube, .docker, ...), the files LEFT_OUT_FILE_RULES name (.env*,
    *token* -- the skill's own .coreomics_token and .pgfarm_token --, .netrc, .pgpass, *.key,
    *.pem, id_*, ...), matched without regard to case, and any file up to SECRET_SCAN_MAX_BYTES
    in which the skill's own secret patterns (notify_slack.contains_secret, the one list) find a
    secret -- the zip travels to collaborators (L5b review M3). A [NOTE] line in MANIFEST.txt
    names each, and another names every file too big to be read for secrets;
  * a copy is made in output/custom/.partial/ (scratch: never zipped) and moved into place only
    when it is whole: a copy that fails leaves the earlier copy as it was (C2 errors_ux M4);
  * running finalize again with the same folders replaces their copies (the record names the new
    source); an item from an earlier finalize stays, and is checked against its checksums.
Standard library only (and report_files and notify_slack, themselves standard library only).
"""
import datetime
import fnmatch
import json
import os
import shutil
import stat
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import report_files as rf  # noqa: E402

CUSTOM_DIR = "custom"                  # output/custom
RECORD = "custom.json"                 # output/custom/custom.json
SCHEMA = "custom_work/1"
# kind -> (the finalize flag, what it is: README / AGENTS.md / MANIFEST.txt wording, a folder?)
KINDS = {"report": ("--report", "a report made outside the skill's steps", False),
         "methods": ("--methods", "Methods written outside the skill's steps", False),
         "code": ("--code-dir", "analysis code written outside the session", True),
         "extra": ("--extra-dir", "other work done outside the session", True)}
# Left out of a folder's copy by default (L5b review M3: a .env with an API key went into the
# zip): version control, caches, and the environments a project builds (a Python venv, npm's,
# renv's package library, RStudio's per-user state) -- by folder name, or path for renv/library --
LEFT_OUT_DIRS = frozenset((".git", ".hg", ".svn", "__pycache__", ".ipynb_checkpoints",
                           ".venv", "venv", "env", "node_modules", ".Rproj.user",
                           # credentials and configuration a home folder holds (C2 security H1):
                           # keys, cloud logins, the skill's own ~/.config/ucdavis-proteomics/
                           # (hive.env), Claude Code's conversations
                           ".ssh", ".aws", ".azure", ".config", ".gnupg", ".kube", ".docker",
                           ".claude"))
LEFT_OUT_DIR_PATHS = (("renv", "library"),)            # <...>/renv/library
# ... and files that hold secrets or keys, by name (fnmatch, without regard to case: ID_RSA and
# Credentials.json too). `*token*` takes the skill's own key files, which no `*.token` matched:
# core_submission.KEY_NAME (.coreomics_token) and the STAN service-account secret
# (ht_manifest/fran_deposit TOKEN_CANDIDATES: .pgfarm_token) -- named too, and a test ties them
LEFT_OUT_FILE_RULES = (".DS_Store", ".env", ".env.*", ".Renviron", "*token*", ".coreomics_token",
                       ".pgfarm_token", "*.key", "*.pem", "id_*", "credentials*", ".netrc",
                       ".pgpass", ".git-credentials", "*.p12", "*.pfx", "*.kdbx", ".Rhistory",
                       ".*_history")
SECRET_SCAN_MAX_BYTES = 1 << 20        # a file up to 1 MiB is read for secrets; a bigger one is
                                       # not read (data, not a config or a script) -- and named in
                                       # a [NOTE] as not read (C2 security H1)
# A source inside one of these folders of a session (this one or any other) is refused (C2
# security M2): logs/ holds the conversation and decisions the zip leaves out by path, and
# input/ and output/ are the session's own files already
SESSION_INTERNAL_DIRS = ("logs", "input", "output")
PARTIAL_DIR = ".partial"               # output/custom/.partial/: a copy being made (scratch_files)
MAX_GB = 10.0                          # above this in all, refused (--custom-max-gb raises it)
WARN_GB = 1.0                          # above this, MANIFEST.txt says it is large (it is zipped)
MANIFEST_PART = "Folded-in work"


class Refused(Exception):
    """A source finalize will not fold in; nothing has been written."""


def custom_dir(output_dir):
    return os.path.join(output_dir, CUSTOM_DIR)


def record_path(output_dir):
    return os.path.join(custom_dir(output_dir), RECORD)


def _real(path):
    return os.path.realpath(os.path.abspath(os.path.expanduser(path)))


def _inside(path, root):
    """`path` is `root` or below it (both real paths)."""
    return path == root or path.startswith(root.rstrip(os.sep) + os.sep)


_LEFT_OUT_DIRS_LOWER = frozenset(d.lower() for d in LEFT_OUT_DIRS)


def _left_out_dir(root, d):
    return d.lower() in _LEFT_OUT_DIRS_LOWER or any(
        (os.path.basename(root).lower(), d.lower()) == pair for pair in LEFT_OUT_DIR_PATHS)


def _left_out_file(name):
    return any(fnmatch.fnmatchcase(name.lower(), rule.lower()) for rule in LEFT_OUT_FILE_RULES)


def _holds_secret(path, size):
    """A file up to SECRET_SCAN_MAX_BYTES in which the skill's secret patterns
    (notify_slack.contains_secret: the ONE list, which redacts the Slack posts, the run log and
    the saved conversation) find one. A bigger file is not read: None (named as not read)."""
    if size > SECRET_SCAN_MAX_BYTES:
        return None
    from notify_slack import contains_secret
    with open(path, encoding="utf-8", errors="replace") as fh:
        return contains_secret(fh.read())


def _homes():
    """The home folder(s) of whoever runs finalize: $HOME, and the account's own from the
    password database where there is one (real paths)."""
    homes = {_real("~")}
    try:
        import pwd
        homes.add(_real(pwd.getpwuid(os.getuid()).pw_dir))
    except (ImportError, KeyError):
        pass                                     # no password database (Windows): $HOME only
    return homes


def _session_internal(src, session):
    """(the session, its sub-folder) when `src` lies inside the logs/, input/ or output/ folder
    of this session or of any other (session.is_session_dir, THE test), else None."""
    import session as session_mod
    path = src
    while True:
        up = os.path.dirname(path)
        if up == path:
            return None
        if os.path.basename(path) in SESSION_INTERNAL_DIRS and (
                up == session or session_mod.is_session_dir(up)):
            return up, os.path.basename(path)
        path = up


def _walk(src):
    """(files [(rel, abs, bytes)], symlinks [(rel, target)], left_out count, left_out names,
    not_scanned [rel], unreadable [rel]) of a folder; links are not followed. Left out:
    LEFT_OUT_DIRS / LEFT_OUT_DIR_PATHS folders (named "<rel>/", their files counted),
    LEFT_OUT_FILE_RULES files, and a file holding a secret (named "<rel> (holds a secret)"). A
    file over SECRET_SCAN_MAX_BYTES is copied but not read (not_scanned); one that cannot be read
    is listed (unreadable: plan() refuses)."""
    files, links, left, named, not_scanned, unreadable = [], [], 0, [], [], []
    for root, dirs, names in os.walk(src, followlinks=False):
        keep = []
        for d in sorted(dirs):
            full = os.path.join(root, d)
            if os.path.islink(full):
                links.append((os.path.relpath(full, src), os.readlink(full)))
            elif _left_out_dir(root, d):
                left += sum(len(fs) for _, _, fs in os.walk(full))
                named.append(os.path.relpath(full, src).replace(os.sep, "/") + "/")
            else:
                keep.append(d)
        dirs[:] = keep
        for n in sorted(names):
            full = os.path.join(root, n)
            rel = os.path.relpath(full, src)
            if os.path.islink(full):
                links.append((rel, os.readlink(full)))
            elif _left_out_file(n):
                left += 1
                if n != ".DS_Store":                 # Finder's own file: counted, not named
                    named.append(rel.replace(os.sep, "/"))
            elif os.path.isfile(full):
                size = os.path.getsize(full)
                if not os.access(full, os.R_OK):
                    unreadable.append(rel.replace(os.sep, "/"))
                    continue
                try:
                    secret = _holds_secret(full, size)
                except OSError:
                    unreadable.append(rel.replace(os.sep, "/"))
                    continue
                if secret:
                    left += 1
                    named.append(rel.replace(os.sep, "/") + " (holds a secret)")
                else:
                    files.append((rel, full, size))
                    if secret is None:
                        not_scanned.append(rel.replace(os.sep, "/"))
    return files, links, left, named, not_scanned, unreadable


def plan(session_dir, output_dir, report=None, methods=None, code_dirs=(), extra_dirs=(),
         max_gb=None):
    """What finalize will fold in: [{"kind", "name", "source", "dest" (session-relative), "dir",
    "files", "symlinks", "left_out", "left_out_names", "not_scanned_names", "bytes",
    "report_state"}]. Raises Refused, with every problem, before anything is written."""
    max_gb = MAX_GB if max_gb is None else float(max_gb)
    session = _real(session_dir)
    custom = _real(custom_dir(output_dir))
    homes = _homes()
    wanted = ([("report", report)] if report else []) + ([("methods", methods)] if methods else []) \
        + [("code", d) for d in code_dirs or ()] + [("extra", d) for d in extra_dirs or ()]
    items, problems, names = [], [], {}
    for kind, given in wanted:
        flag, _, is_dir = KINDS[kind]
        src = _real(given)
        name = os.path.basename(src.rstrip(os.sep))
        if not os.path.exists(src):
            problems.append(f"{flag} {given}: not found")
            continue
        if is_dir != os.path.isdir(src):
            problems.append(f"{flag} {given}: " + ("not a folder" if is_dir else "not a file"))
            continue
        if _inside(src, custom):
            problems.append(f"{flag} {given}: it is inside output/{CUSTOM_DIR}/ already (finalize's "
                            f"own copies); give the original")
            continue
        home = next((h for h in sorted(homes) if _inside(h, src)), None)
        if is_dir and home:
            problems.append(f"{flag} {given}: it is " + ("your home folder" if src == home else
                            f"a folder that holds your home folder ({home})") + " -- its keys, "
                            "tokens and configuration would travel in the session zip; give the "
                            "project folder itself")
            continue
        if _inside(session, src):
            problems.append(f"{flag} {given}: it holds the session itself")
            continue
        # folders and single files alike: --methods "$S/logs/decisions.md" is the conversation's
        # neighbour, which the zip leaves out by path (C2 security M2, verify residual)
        internal = _session_internal(src, session)
        if internal:
            owner, sub = internal
            whose = "this session's" if owner == session else f"the session {owner}'s"
            problems.append(f"{flag} {given}: it is inside {whose} {sub}/ folder -- "
                            + ("Core-internal records (the conversation, decisions) the zip "
                               "leaves out; folding them in would ship them" if sub == "logs"
                               else "a session's own files, not work done outside one")
                            + "; give a folder outside the session")
            continue
        if not name:
            problems.append(f"{flag} {given}: no file or folder name to copy it under")
            continue
        if (kind, name) in names:
            problems.append(f"{flag} {given}: the same name ({name}) as {flag} {names[(kind, name)]}"
                            f" -- both would be output/{CUSTOM_DIR}/{kind}/{name}; rename one")
            continue
        names[(kind, name)] = given
        if is_dir:
            files, links, left, left_names, not_scanned, unreadable = _walk(src)
        else:
            files, links, left, left_names = [(name, src, os.path.getsize(src))], [], 0, []
            not_scanned, unreadable, secret = [], [], False
            if not os.access(src, os.R_OK):
                unreadable = [name]
            else:
                try:
                    secret = _holds_secret(src, files[0][2])
                except OSError:
                    unreadable = [name]
            if secret:                           # one file cannot be left out of itself: refused
                problems.append(f"{flag} {given}: the skill's secret patterns find a secret (a "
                                f"key or token) in it -- the session zip travels to "
                                f"collaborators; take it out and run finalize again")
                continue
            elif secret is None:
                not_scanned = [name]
            if kind == "report":
                beside = os.path.join(os.path.dirname(src), rf.SKILL_REPORT_RECORD)
                if rf.skill_report_state(src)[0] == "skill" and os.path.isfile(beside):
                    files.append((rf.SKILL_REPORT_RECORD, beside, os.path.getsize(beside)))
        if unreadable:
            problems.append(f"{flag} {given}: {len(unreadable)} file(s) cannot be read, so the "
                            f"copy would be partial: " + ", ".join(unreadable[:10])
                            + (" ..." if len(unreadable) > 10 else "")
                            + " -- fix their permissions, or move them out of the folder")
            continue
        items.append({"kind": kind, "name": name, "source": src, "dir": is_dir,
                      "dest": "/".join(("output", CUSTOM_DIR, kind, name)),
                      "files": files, "symlinks": links, "left_out": left,
                      "left_out_names": left_names, "not_scanned_names": not_scanned,
                      "bytes": sum(b for _, _, b in files),
                      "report_state": rf.skill_report_state(src) if kind == "report" else None})
    total = sum(i["bytes"] for i in items)
    if total > max_gb * 1024 ** 3:
        problems.append(f"{total / 1024 ** 3:.1f} GB to fold in, above the {max_gb:g} GB limit: "
                        f"the session (and its zip) would carry it. Fold in a smaller folder, or "
                        f"raise the limit with --custom-max-gb <GB>")
    if problems:
        raise Refused("; ".join(problems))
    return items


def report_gaps(items):
    """finalize's gate: one gap per folded-in report the skill did not write (it goes in only
    under --override-finish)."""
    gaps = []
    for i in items:
        state, why = i.get("report_state") or (None, None)
        if i["kind"] == "report" and state not in ("skill", "legacy_skill"):
            gaps.append(f"the custom report {i['source']} (--report): {why} -- it goes in only "
                        f"as a custom report")
    return gaps


def _now():
    return datetime.datetime.now(datetime.timezone.utc).replace(microsecond=0).isoformat()


def load(output_dir):
    """output/custom/custom.json's items ([] when there is none, or it cannot be read)."""
    rec, _ = rf._read_json(record_path(output_dir))
    items = (rec or {}).get("items")
    return [i for i in items if isinstance(i, dict)] if isinstance(items, list) else []


class CopyFailed(Exception):
    """A planned copy that could not be made whole; the earlier copy is left as it was."""


def _copy_file(src, to):
    """One file, never through a link (C2 security L4): a source that is no longer a regular file
    -- swapped for a link since plan() looked -- is not copied."""
    st = os.lstat(src)
    if not stat.S_ISREG(st.st_mode):
        raise CopyFailed(f"{src} is no longer a regular file (a link?) since finalize planned "
                         f"the copy")
    shutil.copy2(src, to, follow_symlinks=False)


def _remove(path):
    if os.path.isdir(path) and not os.path.islink(path):
        shutil.rmtree(path)
    elif os.path.lexists(path):
        os.remove(path)


def apply(session_dir, output_dir, items, override_reason=None):
    """Copy each item to its dest (replacing an earlier copy there) and write the record, which
    keeps the items of earlier finalizes. -> the record's items.

    Every item is first copied whole into output/custom/.partial/<kind>/<name> (PARTIAL_DIR:
    scratch, never zipped) and only then moved into place (C2 errors_ux M4). A copy that fails --
    a file that cannot be read, a full disk -- removes the partial copies and raises: every
    earlier copy, and the record, stay as they were."""
    by_dest = {i.get("dest"): i for i in load(output_dir)}
    partial = os.path.join(custom_dir(output_dir), PARTIAL_DIR)
    staged = []
    try:
        for i in items:
            work = os.path.join(partial, i["kind"], i["name"])
            _remove(work)
            os.makedirs(work)
            sums = {}
            for rel, src, _ in i["files"]:
                to = os.path.join(work, rel)
                os.makedirs(os.path.dirname(to), exist_ok=True)
                _copy_file(src, to)
                sums[rel] = rf._sha256(to)
            staged.append((i, work, sums))
    except BaseException:
        if os.path.lexists(partial):
            _remove(partial)                     # the partial copies only; nothing else changed
        raise
    for i, work, sums in staged:
        dest = os.path.join(session_dir, *i["dest"].split("/"))
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        if i["dir"]:
            old = work + ".old"
            if os.path.lexists(old):
                _remove(old)
            if os.path.lexists(dest):
                os.replace(dest, old)            # the earlier copy aside until the new is in
            try:
                os.replace(work, dest)
            except OSError:
                if os.path.lexists(old):
                    os.replace(old, dest)        # the earlier copy back where it was
                raise
            if os.path.lexists(old):
                _remove(old)
        else:
            for rel in sums:
                os.replace(os.path.join(work, rel), os.path.join(os.path.dirname(dest), rel))
        state = i.get("report_state")
        by_dest[i["dest"]] = {
            "kind": i["kind"], "name": i["name"], "dest": i["dest"], "source": i["source"],
            "dir": i["dir"], "files": sums, "bytes": i["bytes"],
            "symlinks_not_copied": [rel for rel, _ in i["symlinks"]],
            "left_out": i["left_out"], "left_out_names": i.get("left_out_names") or [],
            "not_scanned_names": i.get("not_scanned_names") or [],
            "folded_in_at": _now(),
            "skill_report": (state[0] in ("skill", "legacy_skill")) if state else None,
            "override_reason": (override_reason if state and state[0] not in
                                ("skill", "legacy_skill") else None)}
    if os.path.lexists(partial):
        _remove(partial)                         # what is left: the emptied folders
    if items:
        rf_items = sorted(by_dest.values(), key=lambda x: x["dest"])
        os.makedirs(custom_dir(output_dir), exist_ok=True)
        rf._write_json(record_path(output_dir), {"schema": SCHEMA, "items": rf_items})
    return load(output_dir)


def changed(session_dir, item):
    """The files of a recorded item that are not as they were folded in: [relative name]."""
    dest = os.path.join(session_dir, *str(item.get("dest") or "").split("/"))
    out = []
    for rel, val in sorted((item.get("files") or {}).items()):
        path = os.path.join(dest, rel) if item.get("dir") else os.path.join(os.path.dirname(dest),
                                                                              rel)
        try:
            same = rf._sha256(path) == val
        except OSError:
            same = False
        if not same:
            out.append(rel)
    return out


def describe(item):
    """One line about an item, for MANIFEST.txt, the README and AGENTS.md."""
    what = KINDS.get(item.get("kind"), ("", "work done outside the session", True))[1]
    n = len(item.get("files") or {})
    bits = [f"{what}, from {item.get('source')}"]
    if item.get("dir"):
        bits.append(f"{n} file(s)")
    if item.get("kind") == "report":
        bits.append("the skill's own report" if item.get("skill_report") else
                    "custom: not written by the skill"
                    + (f" (--override-finish: {item['override_reason']})"
                       if item.get("override_reason") else ""))
    return "; ".join(bits)


def manifest_lines(session_dir, items):
    """[(level, part, note)] for MANIFEST.txt: one line per item, then its symlinks, then any
    copy that changed since it was folded in."""
    out = []
    for i in items:
        part = f"{MANIFEST_PART}: {i['dest']}"
        diff = changed(session_dir, i)
        if diff:
            out.append(("SKIPPED", part, f"changed since it was folded in ({len(diff)} file(s): "
                                         + ", ".join(diff[:5]) + (" ..." if len(diff) > 5 else "")
                                         + f"); fold it in again with {KINDS[i['kind']][0]}"))
            continue
        note = describe(i) + "; sha256 of each file in output/custom/custom.json and the bundle"
        if (i.get("bytes") or 0) > WARN_GB * 1024 ** 3:
            note += f"; {i['bytes'] / 1024 ** 3:.1f} GB, in the session zip too"
        if i.get("left_out"):
            note += (f"; {i['left_out']} file(s) left out by the default rule (secrets, keys, "
                     "environments, version control, caches)")
        out.append(("OK", part, note))
        gone = i.get("left_out_names") or []
        if gone:
            out.append(("NOTE", f"{part} (left out)",
                        f"not copied, by the default rule for secrets, keys, environments, "
                        f"version control and caches: " + ", ".join(gone[:20])
                        + (f" and {len(gone) - 20} more" if len(gone) > 20 else "")
                        + " -- they stay where they were; copy one in by hand only if it belongs "
                          "in the record and holds no secret"))
        unread = i.get("not_scanned_names") or []
        if unread:                               # every one named (C2 security H1)
            out.append(("NOTE", f"{part} (not read for secrets)",
                        f"{len(unread)} file(s) over {SECRET_SCAN_MAX_BYTES / (1 << 20):g} MiB "
                        f"copied WITHOUT being read for secrets: " + ", ".join(unread)
                        + " -- check that none holds a key or token before the session zip "
                          "travels"))
        links = i.get("symlinks_not_copied") or []
        if links:
            out.append(("SKIPPED", f"{part} (symlinks)",
                        f"{len(links)} symlink(s) not copied (never followed): "
                        + ", ".join(links[:5]) + (" ..." if len(links) > 5 else "")
                        + " -- copy the file itself into the folder if it belongs in the record"))
    return out


def bundle_merge(output_dir, items):
    """The checksums into the bundle's run_manifest.json (`custom`), again on every finalize (a
    provenance.py re-run writes the file anew). -> (level, note) for MANIFEST.txt."""
    path = os.path.join(output_dir, rf.BUNDLE_RECORD)
    rec, bad = rf._read_json(path)
    if rec is None:
        return "SKIPPED", (bad or f"no output/{rf.BUNDLE_RECORD.replace(os.sep, '/')} (step 10, "
                                  "provenance.py, has not run)") + ": the checksums are in "\
                                  f"output/{CUSTOM_DIR}/{RECORD} only"
    rec["custom"] = {"record": f"output/{CUSTOM_DIR}/{RECORD}", "recorded_by": "session.py finalize",
                     "items": [{k: i.get(k) for k in ("kind", "dest", "source", "files",
                                                      "skill_report")} for i in items]}
    rf._write_json(path, rec)
    return "OK", f"{len(items)} item(s), {sum(len(i.get('files') or {}) for i in items)} file(s)"
