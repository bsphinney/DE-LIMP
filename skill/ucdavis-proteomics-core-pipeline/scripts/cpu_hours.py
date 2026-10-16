#!/usr/bin/env python3
"""
cpu_hours.py -- THE reserved CPU- and GPU-hours of a search: the most it can use.

Summed over every job script the search submits: CPUs per task x time limit in hours x array
tasks (x tasks per job, where a script asks for more than one). Read from the #SBATCH lines the
generator actually wrote, so the figure describes the jobs that will be submitted -- not a
re-derivation that could drift from them. It is an UPPER BOUND: a job that finishes early uses
less, one cut at its time limit uses exactly this.

GPU-hours the same way, from the GPU lines (--gres=gpu[:type]:N, --gpus / -G, --gpus-per-task,
--gpus-per-node): GPUs per job x time limit x array tasks, 0 for a job that asks for none. No route
of the skill asks for a GPU today (2.12), so it is 0 everywhere; it is said anyway, beside the CPU
figure, so a GPU route added later is counted by the same rule.

Why it exists (Brett, 2026-10-15): on the Core board a request under 10,000 CPU- or GPU-hours is
approved automatically and one of 10,000 or more needs a person, so a Claude can work unattended
overnight. Every search states these figures before it is submitted: run_search.py (every --sbatch
route), diann_parallel.py and radiant_parallel.py print `cpu_hours_reserved` and
`gpu_hours_reserved` in their JSON, and each submit.sh carries them in a header comment. Every job
script the skill writes has an #SBATCH --time, so for a SLURM job neither figure is ever null
(tests/test_cpu_hours.py checks every generator).

  python3 cpu_hours.py <job script> [<job script> ...]      # -> JSON

Stdlib only.
"""
import json
import re
import sys

WHAT = ("the most this search can use: CPUs per task (GPUs per job) x time limit (hours) x array "
        "tasks, summed over every job it submits -- an upper bound; jobs that finish early use less")

_CPUS = re.compile(r"^#SBATCH\s+(?:--cpus-per-task[= ]|-c\s*)(\d+)\s*$", re.M)
_TASKS = re.compile(r"^#SBATCH\s+(?:--ntasks[= ]|-n\s*)(\d+)\s*$", re.M)
_TIME = re.compile(r"^#SBATCH\s+(?:--time[= ]|-t\s*)(\S+)\s*$", re.M)
_ARRAY = re.compile(r"^#SBATCH\s+(?:--array[= ]|-a\s*)(\S+)\s*$", re.M)
_GRES = re.compile(r"^#SBATCH\s+--gres[= ](\S+)\s*$", re.M)
_GPUS = re.compile(r"^#SBATCH\s+(?:--gpus[= ]|-G\s*)(\S+)\s*$", re.M)
_GPUS_PER_TASK = re.compile(r"^#SBATCH\s+--gpus-per-task[= ](\S+)\s*$", re.M)
_GPUS_PER_NODE = re.compile(r"^#SBATCH\s+--gpus-per-node[= ](\S+)\s*$", re.M)
_NODES = re.compile(r"^#SBATCH\s+(?:--nodes[= ]|-N\s*)(\d+)(?:-\d+)?\s*$", re.M)


def _count(spec):
    """The count at the end of a GPU spec: 'a100:2' -> 2, '4' -> 4, 'a100' -> 1."""
    last = spec.rsplit(":", 1)[-1]
    return int(last) if last.isdigit() else 1


def gpus_per_job(text, tasks=1):
    """GPUs one job of this script holds: --gpus (the job's total) wins, else --gpus-per-task x
    tasks, else --gpus-per-node or --gres=gpu x nodes. 0 when none is asked for."""
    if _GPUS.search(text):
        return _count(_GPUS.findall(text)[-1])
    if _GPUS_PER_TASK.search(text):
        return _count(_GPUS_PER_TASK.findall(text)[-1]) * tasks
    nodes = int(_NODES.findall(text)[-1]) if _NODES.search(text) else 1
    if _GPUS_PER_NODE.search(text):
        return _count(_GPUS_PER_NODE.findall(text)[-1]) * nodes
    for gres in _GRES.findall(text):
        for item in gres.split(","):
            if item.split(":", 1)[0] == "gpu":
                return _count(item) * nodes
    return 0


def slurm_hours(value):
    """A SLURM time limit in hours: minutes, minutes:seconds, hours:minutes:seconds, days-hours,
    days-hours:minutes or days-hours:minutes:seconds (sbatch --time). ValueError otherwise."""
    v = value.strip().strip('"\'')
    days = 0
    if "-" in v:
        d, v = v.split("-", 1)
        days = int(d)
        parts = [int(x) for x in v.split(":")] + [0, 0]
        h, m, s = parts[0], parts[1], parts[2]
    else:
        parts = [int(x) for x in v.split(":")]
        if len(parts) == 1:
            h, m, s = 0, parts[0], 0
        elif len(parts) == 2:
            h, m, s = 0, parts[0], parts[1]
        elif len(parts) == 3:
            h, m, s = parts
        else:
            raise ValueError(f"not a SLURM time: {value!r}")
    return days * 24 + h + m / 60 + s / 3600


def array_tasks(spec):
    """How many tasks a SLURM --array spec runs: '0-9%5' -> 10, '1,3,5-7' -> 5, '0-10:2' -> 6."""
    n = 0
    for part in spec.split("%", 1)[0].split(","):
        rng, _, step = part.partition(":")
        lo, _, hi = rng.partition("-")
        lo, hi = int(lo), int(hi or lo)
        n += len(range(lo, hi + 1, int(step or 1)))
    return n


def job_reservation(path):
    """{"script", "cpus", "tasks", "hours", "array_tasks", "cpu_hours"} for one job script; hours
    and cpu_hours are None when the script names no time limit (then nothing bounds it)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        text = fh.read()
    cpus = int(_CPUS.findall(text)[-1]) if _CPUS.search(text) else 1      # SLURM's default
    tasks = int(_TASKS.findall(text)[-1]) if _TASKS.search(text) else 1
    hours = slurm_hours(_TIME.findall(text)[-1]) if _TIME.search(text) else None
    arr = array_tasks(_ARRAY.findall(text)[-1]) if _ARRAY.search(text) else 1
    gpus = gpus_per_job(text, tasks)
    return {"script": path, "cpus": cpus, "tasks": tasks, "hours": hours, "array_tasks": arr,
            "gpus": gpus,
            "cpu_hours": None if hours is None else cpus * tasks * hours * arr,
            "gpu_hours": None if hours is None else gpus * hours * arr}


def reservation(paths):
    """{"cpu_hours_reserved": the sum, rounded to 0.1 (None when a job has no time limit),
    "jobs": [job_reservation()], "what": WHAT} for the job scripts a search submits."""
    jobs = [job_reservation(p) for p in paths if p]
    unbounded = [j["script"] for j in jobs if j["cpu_hours"] is None]
    total = None if unbounded else round(sum(j["cpu_hours"] for j in jobs), 1)
    gpu = None if unbounded else round(sum(j["gpu_hours"] for j in jobs), 1)
    rec = {"cpu_hours_reserved": total, "gpu_hours_reserved": gpu, "jobs": jobs, "what": WHAT}
    if unbounded:
        rec["no_time_limit"] = unbounded
    return rec


def submit_comment(rec):
    """The header comment a submit.sh carries."""
    return (f"# cpu_hours_reserved: {rec['cpu_hours_reserved']}  gpu_hours_reserved: "
            f"{rec['gpu_hours_reserved']} -- {WHAT} (cpu_hours.py; {len(rec['jobs'])} job "
            "script(s))")


def main(argv):
    if not argv or argv[0] in ("-h", "--help"):
        print(__doc__)
        sys.exit(0 if argv else 2)
    print(json.dumps(reservation(argv), indent=2))


if __name__ == "__main__":
    main(sys.argv[1:])
