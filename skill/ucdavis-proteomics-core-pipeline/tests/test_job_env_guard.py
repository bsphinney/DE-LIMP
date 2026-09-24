#!/usr/bin/env python3
"""
Every generated search job ends with the job-end hook: run log, FRAN staging, Slack post
(notify_slack.wrap_job_script). A test that EXECUTES such a script with an inherited environment
is isolated only while $SLURM_JOB_ID is unset -- run the suite inside a HIVE allocation and its
fake searches reach the Core's real run log, FRAN's drop directory and the Core channel.
tests/job_env.py is the one environment that switches all of that off. This meta-test fails any
test module that runs a script with bash and does not use it.

A call counts as "runs a generated script" when it is subprocess.run/Popen/call/check_call/
check_output(["bash" | "sh", <script>, ...]) and <script> is neither a flag ("-c", "-n"), nor a
skill script referenced through SCRIPTS or an UPPER_CASE constant (those are the skill's own
scripts, not generated jobs). Such a call must pass env=, and its module must import job_env.
A script that is generated but is not a search job (e.g. prepare_upload.sbatch) says so with a
`# job_env: not a search job ...` comment on the line above the call.
"""
import ast
import glob
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from job_env import job_env  # noqa: E402

CALLS = {"run", "Popen", "call", "check_call", "check_output"}
EXEMPT = "job_env: not a search job"


def _is_skill_script(node):
    if isinstance(node, ast.Name) and node.id.isupper():
        return True                                    # SCRIPT, ENSURE, HIVE_PATH, ...
    if isinstance(node, ast.Call) and getattr(node.func, "attr", "") == "join" and node.args:
        first = node.args[0]
        return isinstance(first, ast.Name) and first.id == "SCRIPTS"
    return False


def script_calls(path):
    """[(line, has_env, exempt)] for every bash/sh call in `path` that runs a script."""
    with open(path) as fh:
        src = fh.read()
    lines = src.splitlines()
    found = []
    for node in ast.walk(ast.parse(src)):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr in CALLS and isinstance(node.func.value, ast.Name)
                and node.func.value.id == "subprocess" and node.args
                and isinstance(node.args[0], ast.List) and len(node.args[0].elts) >= 2):
            continue
        prog, target = node.args[0].elts[:2]
        if not (isinstance(prog, ast.Constant) and prog.value in ("bash", "sh")):
            continue
        if isinstance(target, ast.Constant) and str(target.value).startswith("-"):
            continue                                   # bash -n / bash -c
        if _is_skill_script(target):
            continue
        near = "\n".join(lines[max(0, node.lineno - 2):node.lineno])
        found.append((node.lineno, any(k.arg == "env" for k in node.keywords), EXEMPT in near))
    return found, ("from job_env import job_env" in src or "import job_env" in src)


class JobEnvGuard(unittest.TestCase):
    def test_every_test_that_runs_a_job_script_uses_job_env(self):
        seen, bad = 0, []
        for path in sorted(glob.glob(os.path.join(HERE, "test_*.py"))):
            calls, imports = script_calls(path)
            for line, has_env, exempt in calls:
                if exempt:
                    continue
                seen += 1
                name = f"{os.path.basename(path)}:{line}"
                if not has_env:
                    bad.append(f"{name} runs a script without env= (use job_env)")
                elif not imports:
                    bad.append(f"{name} runs a script but the module never imports job_env")
        self.assertEqual(bad, [], "\n".join(bad))
        # The scan must actually be finding them: at the time of writing, 14 call sites in 8
        # modules run generated scripts. A scan that matches nothing would pass vacuously.
        self.assertGreaterEqual(seen, 12, "the guard found too few script-running calls")

    def test_the_guard_catches_an_inherited_environment(self):
        with tempfile.TemporaryDirectory() as d:
            bad = os.path.join(d, "test_bad.py")
            with open(bad, "w") as fh:
                fh.write("import os, subprocess\n"
                         "def t(out):\n"
                         "    subprocess.run(['bash', os.path.join(out, 'step5_report.sbatch')])\n"
                         "    subprocess.run(['bash', '-n', 'x.sh'])\n")
            calls, imports = script_calls(bad)
        self.assertEqual(calls, [(3, False, False)])
        self.assertFalse(imports)

    def test_job_env_switches_every_side_effect_off(self):
        base = {"SLURM_JOB_ID": "1", "SLURM_ARRAY_JOB_ID": "2", "SLURM_ARRAY_TASK_ID": "3",
                "SKILL_SLACK_WEBHOOK": "https://hooks.slack.com/services/T/B/x",
                "FRAN_DROP_DIR": "/quobyte/proteomics-grp/fran/incoming", "RECORD_RUN": "on",
                "HIVE_USER": "someone", "PATH": "/usr/bin"}
        env = job_env("/tmp/t", base=base, EXTRA="1")
        for k in ("SLURM_JOB_ID", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID",
                  "SKILL_SLACK_WEBHOOK", "HIVE_USER"):
            self.assertNotIn(k, env)
        self.assertEqual((env["SKILL_SLACK"], env["FRAN_DEPOSIT"], env["RECORD_RUN"]),
                         ("0", "off", "off"))
        # never the live FRAN database, even on HIVE as an account that can read the token
        self.assertEqual((env["FRAN_HEALTH"], env["FRAN_CORPUS_QUERY"]), ("off", "off"))
        self.assertEqual(env["SKILL_RUNS_DIR"], "/tmp/t/skill_runs")
        self.assertEqual(env["FRAN_DROP_DIR"], "/tmp/t/fran_drop")
        self.assertEqual(env["HIVE_ENV_FILE"], "/nonexistent/hive.env")
        self.assertEqual((env["PATH"], env["EXTRA"]), ("/usr/bin", "1"))


if __name__ == "__main__":
    unittest.main()
