"""A fake FragPipe install for tests that generate (and run) a FragPipe search (2.15).

run_search.py takes a FragPipe search's tools folder, Philosopher and bundled DIA-NN from tools.json
`fragpipe_tools` (acquire_tools.sh writes it; INTERFACE_F2.md), and the job checks the Java
FragPipe's launcher runs ($JAVA_HOME/bin/java, else java on PATH) before anything else. This builds
the pieces under a test's folder:

    fp = fake_install(d)                       # {"fragpipe_tools": {...}, "java_home": ...}
    tools_json(d, fragpipe_cmd, fp)            # tools.json with fragpipe + fragpipe_tools
    env["JAVA_HOME"] = fp["java_home"]         # the job's Java (JAVA="1.8.0_292" for an old one)

The fake Philosopher does what Philosopher 5.1.3 does for `workspace --init|--clean --nocheck` and
`database --custom <fasta> --prefix <tag>`: it writes <YYYY-MM-DD>-decoys-<fasta name>.fas in the
current folder, the targets then one `<tag>` decoy per target (FAKE_PHI_HALF=1: only the first).
"""
import json
import os

FAKE_JAVA = """#!/bin/sh
v="${FAKE_JAVA_VERSION:-%s}"
echo "openjdk version \\"$v\\" 2026-01-01" >&2
echo "OpenJDK Runtime Environment (build $v)" >&2
"""

FAKE_PHILOSOPHER = r"""#!/bin/bash
cmd="$1"; shift
case "$cmd" in
  workspace)
    case " $* " in *" --init "*) mkdir -p .meta;; *" --clean "*) rm -rf .meta;; esac
    exit 0;;
  database)
    custom=""; prefix="rev_"
    while [ $# -gt 0 ]; do case "$1" in
      --custom) custom="$2"; shift;; --prefix) prefix="$2"; shift;;
      --contam) echo "contam must never be asked for" >&2; exit 9;; esac; shift; done
    [ -f "$custom" ] || { echo "no such file $custom" >&2; exit 2; }
    out="$(date +%F)-decoys-$(basename "$custom").fas"
    awk -v p="$prefix" -v half="${FAKE_PHI_HALF:-0}" '
      /^>/ { n++; h[n] = substr($0, 2); next } { s[n] = s[n] $0 }
      END { for (i = 1; i <= n; i++) print ">" h[i] "\n" s[i];
            for (i = 1; i <= n; i++) { if (half == "1" && i > 1) break;
              r = ""; for (j = length(s[i]); j > 0; j--) r = r substr(s[i], j, 1);
              print ">" p h[i] "\n" r } }' "$custom" > "$out"
    echo "philosopher: wrote $out"
    exit 0;;
esac
echo "fake philosopher: unknown command $cmd" >&2
exit 1
"""


def _write(path, text, mode=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)
    if mode:
        os.chmod(path, mode)
    return path


def fake_install(d, java="21.0.4", licensed=True):
    """The pieces of a FragPipe install under `d`/fp_install -> {"fragpipe_tools", "java_home",
    "tools_folder", "philosopher"}; `fragpipe_tools` in INTERFACE_F2.md's shape."""
    root = os.path.join(d, "fp_install")
    tools = os.path.join(root, "tools")
    phi = _write(os.path.join(tools, "Philosopher", "philosopher-v5.1.3-RC9"), FAKE_PHILOSOPHER,
                 0o755)
    diann = _write(os.path.join(tools, "diann", "1.8.2_beta_8", "linux", "diann-1.8.1.8"),
                   "#!/bin/sh\nexit 0\n", 0o755)
    if licensed:
        _write(os.path.join(tools, "MSFragger-4.4.1.jar"), "jar")
        _write(os.path.join(tools, "IonQuant-1.11.20.jar"), "jar")
    py = _write(os.path.join(d, "easypqp_venv", "bin", "python"), "#!/bin/sh\nexit 0\n", 0o755)
    jhome = os.path.join(d, "fake_jdk")
    _write(os.path.join(jhome, "bin", "java"), FAKE_JAVA % java, 0o755)
    rec = {"binary": None, "folder": tools if licensed else None,
           "source": "FRAGPIPE_TOOLS_FOLDER",
           "jars": {"msfragger": {"path": os.path.join(tools, "MSFragger-4.4.1.jar"),
                                  "version": "4.4.1"} if licensed else None,
                    "ionquant": {"path": os.path.join(tools, "IonQuant-1.11.20.jar"),
                                 "version": "1.11.20"} if licensed else None,
                    "diatracer": None},
           "philosopher": {"path": phi, "version": "5.1.3-RC9"},
           "diann": {"path": diann, "version": "1.8.2 beta 8"},
           "python": {"path": py, "version": "3.11.0", "easypqp": "0.1.53"},
           "java": {"version": f'openjdk version "{java}"', "major": int(java.split(".")[0]),
                    "ok": True},
           "licensed": licensed, "ready": licensed,
           "why": "a test's fake install" + ("" if licensed else ": no MSFragger/IonQuant jars"),
           "skipped": []}
    return {"fragpipe_tools": rec, "java_home": jhome, "tools_folder": tools, "philosopher": phi}


def tools_json(d, fragpipe_cmd, fp=None, name="tools.json"):
    """<d>/tools.json naming `fragpipe_cmd` and (when `fp` is given) its fragpipe_tools."""
    data = {"fragpipe": fragpipe_cmd, "versions": {"fragpipe": "24.0"}}
    if fp is not None:
        fp["fragpipe_tools"]["binary"] = fragpipe_cmd
        data["fragpipe_tools"] = fp["fragpipe_tools"]
    return _write(os.path.join(d, name), json.dumps(data))
