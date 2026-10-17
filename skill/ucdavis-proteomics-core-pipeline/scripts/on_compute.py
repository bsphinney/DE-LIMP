#!/usr/bin/env python3
"""
on_compute.py -- run one command on a SLURM compute node, never on the cluster's login node.

HPC@UCD's agent skill (ucd-hpc), mandatory on Hive, Farm and Franklin since 2026-10-16: "No
analysis on the login node ... anything past a few seconds of CPU or much I/O -- runs in srun
(blocking) or sbatch, and you wait", with --account, --partition, --time and memory always given.
The searches and DEs are sbatch jobs already; this is the skill's one way to send any OTHER
command to a compute node (a walk of a share, an install, a zip of a session, hashing raw files):

  python3 on_compute.py [--time MIN] [--mem SIZE] [--cpus N] [--name NAME] [--preemptible]
                        [--partition P] [--account A] [--qos Q] [--why TEXT] [--print]
                        -- COMMAND [ARG ...]
  python3 on_compute.py --is-login-node          # exit 0 on a cluster login node, else 1

  * On a cluster login node (no SLURM job, and sbatch on PATH, SLURM at HPC@UCD's fixed CVMFS
    path, or HPCCF's /opt/hpccf/etc/cluster_name.conf): ONE blocking `srun`, the queue written
    out -- --account, --partition, --qos -- and --time (minutes, default 30), --mem (default 4G),
    --ntasks=1, --cpus-per-task (default 1). stdin, stdout and the exit status are the
    command's. The srun line is printed on stderr first.
  * The queue is one that does not preempt: an srun that is preempted is cancelled, not
    requeued, and these commands install, copy or write. genome-center-grp/high when the account
    has it, else another account's `high`; only an account with nothing but `low` gets low (and
    is told so). --preemptible lets run_search.slurm_queue (the skill's queue rule for jobs)
    choose, `low` included. A named --partition/--account/--qos is used as given.
  * Inside a job (SLURM_JOB_ID set) or where there is no SLURM (a laptop): the command runs here,
    unchanged, so one documented line works everywhere.
  * --print: print the srun line on stdout and run nothing.
  * The srun BLOCKS until the job ends (queue wait + the work). An agent runs the call in the
    background, or with a tool timeout of at least --time plus 10 minutes: a killed ssh cancels
    the srun half way.
  * --header LOGDIR (no command): print the top of a job script the skill does not generate (a
    DE, set tests, figures): the same queue, --time/--mem/--cpus, one log per job
    (LOGDIR/NAME_%j.log) and `set -euo pipefail`. Append the commands, then sbatch it. A step
    that must survive a closed laptop, or runs past ~10 minutes, is such a job, not an srun.

From a laptop, `bash hive_exec.sh --job [--time MIN] [--mem SIZE] [--cpus N] '<command>'` runs
this on HIVE. A script that does heavy work calls ensure_compute_node() first, which re-runs it
this way when it was started on a login node.

Exit: the command's status; 2 usage (no command, a size without a unit).
"""
import os
import re
import shlex
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

DEFAULT_MINUTES = 30
DEFAULT_MEM = "4G"
SIZE = re.compile(r"^[1-9][0-9]*[KMGT]$")      # a bare number is MEGABYTES to SLURM (ucd-hpc)


# Where SLURM lives on all three HPC@UCD clusters (ucd-hpc: "a fixed path ... never something to
# search for"); a non-login shell, or one after `module purge`, has no sbatch on PATH. HPCCF's
# own marker of a cluster node. Both overridable for the tests.
SLURM_BIN_DEFAULT = "/cvmfs/hpc.ucdavis.edu/sw/spack/environments/core/view/generic/slurm/bin"
CLUSTER_CONF_DEFAULT = "/opt/hpccf/etc/cluster_name.conf"
# Set in the environment of every srun this starts: whatever runs under it is never sent again
REEXEC_ENV = "PROTEOMICS_ON_COMPUTE"


def _slurm_bin(env):
    return env.get("HPCCF_SLURM_BIN") or SLURM_BIN_DEFAULT


def on_login_node(env=None):
    """A cluster login/submit node: no SLURM job here, and SLURM is (sbatch on PATH, or at the
    fixed CVMFS path) or HPCCF's cluster file is. THE definition: detect_acquisition.
    on_cluster_login_node, setup.sh and ensure_dotnet8.sh (--is-login-node) ask this."""
    env = os.environ if env is None else env
    if env.get("SLURM_JOB_ID") or env.get(REEXEC_ENV):
        return False
    return (shutil.which("sbatch", path=env.get("PATH")) is not None
            or os.access(os.path.join(_slurm_bin(env), "sbatch"), os.X_OK)
            or os.path.isfile(env.get("UCD_HPC_CLUSTER_CONF") or CLUSTER_CONF_DEFAULT))


def _srun_path(env=None):
    env = os.environ if env is None else env
    found = shutil.which("srun", path=env.get("PATH"))
    if found:
        return found
    for d in (_slurm_bin(env), os.path.dirname(shutil.which("sbatch", path=env.get("PATH")) or "")):
        if d and os.access(os.path.join(d, "srun"), os.X_OK):
            return os.path.join(d, "srun")
    return "srun"


def pick_queue(cpus=1, preemptible_ok=False, partition=None, account=None, qos=None):
    """(partition, account, qos) for an srun. Not preemptible unless allowed: run_search.
    slurm_queue moves work to publicgrp/low when the account is at its `high` cap, which suits a
    requeued sbatch job and not an srun (preempted = cancelled, half way through an install or a
    copy). A named queue is used as given (the skill's rule checks it)."""
    import run_search
    if preemptible_ok or partition or account or qos:
        return run_search.slurm_queue(partition, account, qos, peak_cpus=cpus,
                                      preemptible_ok=preemptible_ok)
    steady = [(a, p, q) for a, p, q in run_search._associations()
              if p and p != "low" and not p.startswith("gpu")
              and not (a == "publicgrp" and p == "high" and cpus > 8)]   # the free tier's cap
    if steady:
        steady.sort(key=lambda x: (x[:2] != ("genome-center-grp", "high"), x[1] != "high"))
        a, p, q = steady[0]
        return p, a, (q if q and "," not in q else None)
    p, a, q = run_search.slurm_queue(peak_cpus=cpus, preemptible_ok=False)
    if p == "low":
        sys.stderr.write("[on_compute] this account can use only `low`, which is preemptible: a "
                         "preempted srun is cancelled, not requeued -- run it again if it is.\n")
    return p, a, q


def check_size(mem):
    if not SIZE.match(str(mem)):
        raise ValueError(f"--mem {mem!r} needs a unit (K, M, G or T), e.g. 4G: SLURM reads a bare "
                         f"number as megabytes")
    return mem


def srun_argv(command, minutes=DEFAULT_MINUTES, mem=DEFAULT_MEM, cpus=1, name=None,
              preemptible_ok=False, partition=None, account=None, qos=None):
    """The srun line for `command` (a list), the queue written out. Asks SLURM (sacctmgr) which
    queue this user can submit to (pick_queue)."""
    check_size(mem)
    minutes, cpus = int(minutes), int(cpus)
    if minutes < 1 or cpus < 1:
        raise ValueError("--time and --cpus must be at least 1")
    part, acct, q = pick_queue(cpus, preemptible_ok, partition, account, qos)
    name = name or "pp_" + re.sub(r"[^A-Za-z0-9_.-]", "_", os.path.basename(str(command[0])))[:40]
    argv = [_srun_path()]
    argv += [f"--account={acct}"] if acct else []
    argv += [f"--partition={part}"] if part else []
    argv += [f"--qos={q}"] if q else []
    argv += [f"--time={minutes}", f"--mem={mem}", "--ntasks=1", f"--cpus-per-task={cpus}",
             f"--job-name={name}"]
    return argv + list(command)


def shown(argv):
    """The line as a person would type it: srun by name, not by its path."""
    return " ".join(shlex.quote(a) for a in ["srun"] + list(argv[1:]))


def _exec(argv, srun=False):
    sys.stdout.flush()
    sys.stderr.flush()
    if srun:
        os.environ[REEXEC_ENV] = "1"     # srun passes it on: nothing under it is sent again
    try:
        os.execvp(argv[0], argv)
    except OSError as e:                 # no srun beside sbatch, or no such command: say so, once
        sys.stderr.write(f"on_compute.py: cannot run {argv[0]}: {e.strerror or e}\n")
        sys.exit(127)


def ensure_compute_node(why, minutes=DEFAULT_MINUTES, mem=DEFAULT_MEM, cpus=1,
                        preemptible_ok=False, argv=None, script=None):
    """First thing in a script whose work does not belong on a login node: there, replace this
    process with the same command under srun (never returns). Anywhere else, return. `script`
    (the caller's __file__): only when that script IS this process's program -- a test or
    another script calling its main() in-process is left alone, never re-run under srun."""
    if not on_login_node():
        return
    if script and (not sys.argv or os.path.realpath(sys.argv[0]) != os.path.realpath(script)):
        return
    argv = list(sys.argv if argv is None else argv)
    cmd = [sys.executable or "python3", os.path.abspath(argv[0])] + argv[1:]
    line = srun_argv(cmd, minutes, mem, cpus, preemptible_ok=preemptible_ok,
                     name="pp_" + os.path.splitext(os.path.basename(argv[0]))[0])
    sys.stderr.write(f"[on_compute] {why}: that is work for a compute node, not this login node "
                     f"(HPC@UCD's ucd-hpc rules). Running it there and waiting for it:\n"
                     f"  {shown(line)}\n")
    _exec(line, srun=True)


def _usage(msg):
    sys.stderr.write(f"on_compute.py: {msg}\n"
                     "usage: on_compute.py [--time MIN] [--mem SIZE] [--cpus N] [--name NAME] "
                     "[--preemptible] [--partition P] [--account A] [--qos Q] [--why TEXT] "
                     "[--print] -- COMMAND [ARG ...]\n"
                     "       on_compute.py --header LOGDIR --name NAME [--time MIN] [--mem SIZE] "
                     "[--cpus N] [--no-preempt]\n"
                     "       on_compute.py --is-login-node\n")
    return 2


def sbatch_header(name, log_dir, minutes=DEFAULT_MINUTES, mem=DEFAULT_MEM, cpus=1,
                  preemptible_ok=True, partition=None, account=None, qos=None):
    """The top of a job script the skill does not generate (a DE, set tests, figures): the queue
    from run_search.slurm_queue written out, --time, --mem with its unit, --ntasks=1,
    --cpus-per-task, one log per job (%j) in log_dir -- ucd-hpc's canonical job script -- and
    `set -euo pipefail`. A preemptible queue (low) gets --requeue: a preempted job starts over."""
    check_size(mem)
    minutes, cpus = int(minutes), int(cpus)
    if minutes < 1 or cpus < 1:
        raise ValueError("--time and --cpus must be at least 1")
    if not re.match(r"^[A-Za-z0-9_.-]+$", name or ""):
        raise ValueError("--name: letters, digits, _ . - only")
    log = os.path.join(os.path.abspath(log_dir), f"{name}_%j.log")
    if re.search(r'["\n]', log):
        raise ValueError(f"the log folder holds a quote or a line break, which #SBATCH --output "
                         f"cannot carry: {log_dir!r}")
    import run_search
    part, acct, q = run_search.slurm_queue(partition, account, qos, peak_cpus=cpus,
                                           preemptible_ok=preemptible_ok)
    L = ["#!/bin/bash -l", f"#SBATCH --job-name={name}"]
    L += [f"#SBATCH --account={acct}"] if acct else []
    L += [f"#SBATCH --partition={part}"] if part else []
    L += [f"#SBATCH --qos={q}"] if q else []
    L += [f"#SBATCH --time={minutes}", "#SBATCH --ntasks=1", f"#SBATCH --cpus-per-task={cpus}",
          f"#SBATCH --mem={mem}", f'#SBATCH --output="{log}"', f'#SBATCH --error="{log}"']
    L += ["#SBATCH --requeue"] if part == "low" else []
    L += ["# written by on_compute.py --header (the skill's queue rule; HPC@UCD's ucd-hpc job script)",
          "set -euo pipefail"]
    return L


def main(args):
    i = args.index("--") if "--" in args else len(args)
    opts, cmd = args[:i], args[i + 1:]
    if "-h" in opts or "--help" in opts:
        print(__doc__.strip())
        return 0
    o = {"time": DEFAULT_MINUTES, "mem": DEFAULT_MEM, "cpus": 1, "name": None, "partition": None,
         "account": None, "qos": None, "why": None, "header": None}
    flags = {"print": False, "no-preempt": False, "preemptible": False, "is-login-node": False}
    k = 0
    while k < len(opts):
        f = opts[k]
        if f.startswith("--") and f[2:] in flags:
            flags[f[2:]] = True
            k += 1
        elif f.startswith("--") and f[2:] in o and k + 1 < len(opts):
            o[f[2:]] = opts[k + 1]
            k += 2
        else:
            return _usage(f"unknown or incomplete option {f!r}")
    if flags["is-login-node"]:
        return 0 if on_login_node() else 1
    if o["header"] is not None:
        try:
            print("\n".join(sbatch_header(o["name"] or "pp_job", o["header"], o["time"], o["mem"],
                                           o["cpus"], not flags["no-preempt"], o["partition"],
                                           o["account"], o["qos"])))
        except ValueError as e:
            return _usage(str(e))
        return 0
    if not cmd:
        return _usage("no command: put it after --")
    try:
        check_size(o["mem"])
        if int(o["time"]) < 1 or int(o["cpus"]) < 1:
            raise ValueError("--time and --cpus must be at least 1")
        line = srun_argv(cmd, o["time"], o["mem"], o["cpus"], o["name"],
                         preemptible_ok=flags["preemptible"], partition=o["partition"],
                         account=o["account"], qos=o["qos"]) if (flags["print"] or on_login_node()) \
            else None
    except ValueError as e:
        return _usage(str(e))
    if flags["print"]:
        print(shown(line))
        return 0
    if line is None:                     # inside a job, or no SLURM here: run it in place
        _exec(cmd)
    sys.stderr.write(f"[on_compute] {o['why'] + ': ' if o['why'] else ''}running on a compute node "
                     f"(HPC@UCD's ucd-hpc rules), waiting for it:\n  {shown(line)}\n")
    _exec(line, srun=True)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
