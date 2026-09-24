"""The ONE environment for a test that EXECUTES a generated job script, or runs finalize.

Every search job this skill writes ends with the job-end hook (notify_slack.wrap_job_script):
it logs the run (record_run.py), stages a finished search for FRAN (fran_deposit.py stage) and
posts to the Core's Slack channel. finalize does the run log and the post too. The hook does
nothing outside SLURM -- but a suite run inside a HIVE allocation HAS a $SLURM_JOB_ID, and then a
test's fake search would be logged in the Core's real run log, staged into FRAN's real drop
directory and announced in the Core channel. Isolation must not depend on where the suite runs,
so every such test builds its environment here:

    from job_env import job_env
    env = job_env(tmpdir)                      # from os.environ
    env = job_env(tmpdir, base=env, FOO="1")   # from another env, plus extras

tests/test_job_env_guard.py fails any test module that runs a generated job script without it.
"""
import os

# Dropped from the base environment before the switches below are set, so nothing inherited --
# a real webhook, a HIVE login, a job id, a FRAN or run-log override -- can leak through.
_DROPPED_PREFIXES = ("SLURM_", "SKILL_SLACK", "FRAN_", "RECORD_RUN", "SKILL_RUNS_DIR", "HIVE_")


def job_env(tmpdir, base=None, **extra):
    """A copy of `base` (default os.environ) with every job-end side effect switched off:
    no Slack post, no FRAN staging, no run log, no HIVE login, no SLURM job -- and the run log
    and FRAN drop directory pointed inside `tmpdir` in case anything ignores its switch."""
    env = {k: v for k, v in (os.environ if base is None else base).items()
           if not k.startswith(_DROPPED_PREFIXES)}
    env.update(SKILL_SLACK="0", FRAN_DEPOSIT="off", FRAN_HEALTH="off", RECORD_RUN="off",
               SKILL_RUNS_DIR=os.path.join(tmpdir, "skill_runs"),
               FRAN_DROP_DIR=os.path.join(tmpdir, "fran_drop"),
               HIVE_ENV_FILE="/nonexistent/hive.env")
    env.update({k: str(v) for k, v in extra.items()})
    return env
