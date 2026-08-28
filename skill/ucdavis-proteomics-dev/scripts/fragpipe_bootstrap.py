#!/usr/bin/env python3
"""
fragpipe_bootstrap.py -- make headless FragPipe actually run on a fresh account.

THE PROBLEM
-----------
Headless FragPipe refuses the Spectral Library Generation module on any account that
has never opened the GUI:

    ERROR - Spectral Library Generation module was not configured correctly.
            Please make sure that Python and FragPipe-SpecLib have been installed.

This is NOT fixable from the command line. It persists with a correct
--config-tools-folder, a correct --config-diann, python with easypqp on PATH, and
--config-python pointed at either the venv or its python binary. (--config-python is
also mis-documented: its help says "the Python directory" but FragPipe execs the path
directly, so handing it a directory gives `Cannot run program ...: error 13`.)

The state FragPipe actually reads is a per-user cache:

    ~/.config/FragPipe/fragpipe/<version>/fragpipe-ui.cache
        fragpipe-config.bin-python=/path/to/python      <-- the key that matters
        fragpipe-config.bin-diann=/path/to/diann-binary

Copying a configured user's cache fixes it instantly; this script writes the same
two keys without needing a GUI, a display, or another user's home directory.

Every student on a fresh cluster account hits this, so it is not an edge case.

Usage
  python3 fragpipe_bootstrap.py --python <python-with-easypqp> \\
      [--diann <diann binary>] [--fragpipe-version 24.0] [--check]

  # UC Davis HIVE (paths already installed):
  python3 fragpipe_bootstrap.py \\
      --python /quobyte/proteomics-grp/fragpipe24/easypqp_venv/bin/python \\
      --diann  /quobyte/proteomics-grp/fragpipe24/fragpipe-24.0/tools/diann/1.8.2_beta_8/linux/diann-1.8.1.8

Then run FragPipe with --config-tools-folder <fragpipe>/tools (the folder holding
MSFragger, IonQuant, diatracer AND speclib/ — not a jars-only directory).
"""
import argparse, os, subprocess, sys

KEY_PY = "fragpipe-config.bin-python"
KEY_DIANN = "fragpipe-config.bin-diann"


def cache_path(version):
    return os.path.join(os.path.expanduser("~"), ".config", "FragPipe",
                        "fragpipe", version, "fragpipe-ui.cache")


def verify_python(py):
    """A python that cannot import easypqp will fail later, less legibly."""
    if not os.path.isfile(py):
        return f"not a file: {py} (pass the python BINARY, not the venv directory)"
    if not os.access(py, os.X_OK):
        return f"not executable: {py}"
    try:
        r = subprocess.run([py, "-c", "import easypqp"], capture_output=True, text=True, timeout=90)
    except Exception as e:
        return f"could not run {py}: {e}"
    if r.returncode != 0:
        return (f"{py} cannot import easypqp. Install it into that interpreter:\n"
                f"    {py} -m pip install easypqp")
    return None


def read_cache(p):
    out = {}
    if os.path.exists(p):
        for line in open(p):
            if "=" in line and not line.lstrip().startswith("#"):
                k, _, v = line.partition("=")
                out[k.strip()] = v.rstrip("\n")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--python", help="python interpreter that can import easypqp")
    ap.add_argument("--diann", help="DIA-NN binary FragPipe should use for quantification")
    ap.add_argument("--fragpipe-version", default="24.0")
    ap.add_argument("--check", action="store_true", help="report current state, write nothing")
    a = ap.parse_args()

    p = cache_path(a.fragpipe_version)
    cur = read_cache(p)

    if a.check or not (a.python or a.diann):
        import json
        print(json.dumps({
            "cache": p,
            "exists": os.path.exists(p),
            KEY_PY: cur.get(KEY_PY),
            KEY_DIANN: cur.get(KEY_DIANN),
            "speclib_ready": bool(cur.get(KEY_PY)),
            "hint": ("Pass --python <interpreter with easypqp> to write it. Without "
                     f"{KEY_PY} set, headless FragPipe rejects the Spectral Library "
                     "Generation module no matter what CLI flags you pass."),
        }, indent=2))
        return

    if a.python:
        err = verify_python(a.python)
        if err:
            sys.exit(f"[fragpipe_bootstrap] {err}")
        cur[KEY_PY] = os.path.abspath(a.python)
    if a.diann:
        if not os.path.isfile(a.diann):
            sys.exit(f"[fragpipe_bootstrap] no such DIA-NN binary: {a.diann}")
        cur[KEY_DIANN] = os.path.abspath(a.diann)

    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        fh.write(f"# FragPipe v{a.fragpipe_version} ui cache\n")
        fh.write("# written by fragpipe_bootstrap.py — headless speclib configuration\n")
        for k in sorted(cur):
            fh.write(f"{k}={cur[k]}\n")

    import json
    print(json.dumps({
        "wrote": p,
        KEY_PY: cur.get(KEY_PY),
        KEY_DIANN: cur.get(KEY_DIANN),
        "next": ("FragPipe should now accept the speclib module headlessly. Run it with "
                 "--config-tools-folder <fragpipe>/tools (the folder that holds MSFragger, "
                 "IonQuant, diatracer AND speclib/)."),
    }, indent=2))


if __name__ == "__main__":
    main()
