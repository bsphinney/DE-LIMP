#!/usr/bin/env python3
"""
slurm_states.py -- the skill's ONE list of what a Slurm job state means (2.15, C2 review ux H1).

Every reader of `sacct` / `squeue` states takes them from here: watch_run.sh (chain mode and
single-job mode, through `regex`), record_run.py and checkpoint.py. Before 2.15 each kept its own
list, and the watcher's chain mode read PREEMPTED, DEADLINE, BOOT_FAIL and REVOKED as "done, not
failed" while its single-job mode read them as running forever.

Three groups:
  OK      COMPLETED: the job finished and its script exited 0.
  BAD     the job ended and did NOT finish its work: it is done AND failed. CANCELLED is here too
          (a chain with a cancelled step is not a finished search).
  LIVE    the job has not ended: not done, not failed. REQUEUED / REQUEUE_FED / REQUEUE_HOLD and
          SUSPENDED are here: a job Slurm requeued will run again by itself.

Requeue after preemption. `sacct -j <id>` WITHOUT -D/--duplicates returns only the most recent
record of a job id (the sacct man page: "When data for specific jobs are requested with the --jobs
option, sacct returns the most recent job with that number. This behavior can be overridden by
specifying --duplicates"). A job preempted and then requeued therefore reads as its new instance
(PENDING, then RUNNING) -- not done, not failed -- and the skill never passes -D. A PREEMPTED state
in that most recent record means Slurm did NOT requeue it: the job ended there, so it is BAD.

sacct's default State column is 10 characters wide and cuts a longer name to "<first 9>+"
("OUT_OF_MEMORY" -> "OUT_OF_ME+", "SPECIAL_EXIT" -> "SPECIAL_E+"). `regex` therefore matches the
first 9 characters of each name, so a cut name and a full one are read the same.

  python3 slurm_states.py regex bad|live|ok|requeued   # an ERE for grep -E / awk
  python3 slurm_states.py classify <STATE>             # ok | bad | live | unknown
"""
import sys

OK_STATES = frozenset({"COMPLETED"})
BAD_STATES = frozenset({
    "FAILED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "CANCELLED",
    "PREEMPTED",      # preempted and NOT requeued (a requeued job reads as its new instance)
    "DEADLINE",       # killed at its --deadline
    "BOOT_FAIL",      # the node failed to boot
    "REVOKED",        # a federated sibling ran it elsewhere; this copy will never run
    "SPECIAL_EXIT",   # requeued in a held special state: it will not run by itself
})
REQUEUED_STATES = frozenset({"REQUEUED", "REQUEUE_FED", "REQUEUE_HOLD"})
LIVE_STATES = frozenset({
    "PENDING", "RUNNING", "SUSPENDED", "RESIZING", "CONFIGURING", "COMPLETING", "STOPPED",
    "SIGNALING", "STAGE_OUT",
}) | REQUEUED_STATES

GROUPS = {"ok": OK_STATES, "bad": BAD_STATES, "live": LIVE_STATES, "requeued": REQUEUED_STATES}
SACCT_WIDTH = 9          # characters sacct keeps before its "+" in the default State column


def norm(state):
    """'CANCELLED by 123' / 'OUT_OF_ME+' -> the bare name sacct meant ('CANCELLED', 'OUT_OF_ME')."""
    return (state or "").strip().split()[0].rstrip("+") if (state or "").strip() else ""


def classify(state):
    """ok | bad | live | unknown, for one state as sacct or squeue printed it (cut or not)."""
    s = norm(state)
    if not s:
        return "unknown"
    for name, group in (("ok", OK_STATES), ("bad", BAD_STATES), ("live", LIVE_STATES)):
        if s in group or any(g.startswith(s) and len(s) >= SACCT_WIDTH for g in group) \
                or any(s.startswith(g) for g in group):
            return name
    return "unknown"


def regex(group):
    """An ERE matching any state of the group, full or as sacct cuts it (first 9 characters)."""
    names = sorted({g[:SACCT_WIDTH] for g in GROUPS[group]})
    return "|".join(names)


def main(argv):
    if len(argv) == 2 and argv[0] == "regex" and argv[1] in GROUPS:
        print(regex(argv[1]))
        return 0
    if len(argv) == 2 and argv[0] == "classify":
        print(classify(argv[1]))
        return 0
    sys.stderr.write(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
