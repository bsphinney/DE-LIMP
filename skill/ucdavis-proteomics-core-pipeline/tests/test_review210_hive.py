#!/usr/bin/env python3
"""
Independent pre-release review of skill 2.10.0 (HIVE operations and the DIA-NN chain): one failing
test per confirmed issue. Written to drop into skill/ucdavis-proteomics-core-pipeline/tests/; run
from elsewhere with REVIEW_SKILL_DIR=<that skill folder>.

Evidence for the first two (HIVE, read-only sacct, 2026-10-01; every chain task since 2026-07-01):
  * TIME: 6,182 step-2 tasks ran at 16 CPUs; 414 took > 60 min, 84 > 120 min, and 30 of the 91
    arrays with >= 8 files had a task > 60 min. 2.10 sizes those arrays to 8 CPUs on
    genome-center-grp/high (1.96x slower per file, the fixer's own measurement) and keeps 2 h.
  * MEMORY: MaxRSS on HIVE includes page cache (Slurm cgroup/v2 reads memory.current;
    JobAcctGatherParams=NoOverMemoryKill, no NoFileCache), so the lower bound of a task's own
    memory is MaxRSS - MaxDiskRead - MaxDiskWrite. That bound is > 32 GB for 65 step-2 tasks
    (22 arrays) and 9 step-1b jobs, up to 61 GB (24265766_0: a blank-like injection, 1.2 GB .d,
    the standard 1.7 GB human library, MaxRSS 64.0 GB at its 64 GB limit, 2.8 GB read).
Everything SLURM is stubbed; nothing contacts a cluster, Slack, FRAN or the run log.
"""
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import time
import unittest
import zipfile

SKILL = os.environ.get("REVIEW_SKILL_DIR") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(SKILL, "scripts")
TESTS = os.path.join(SKILL, "tests")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, TESTS)
from job_env import job_env  # noqa: E402

SACCTMGR = """#!/bin/sh
case "$*" in
  *"show assoc"*)
    echo 'genome-center-grp|high|genome-center-grp-high-qos'
    echo 'publicgrp|low|publicgrp-low-qos' ;;
  *"show qos"*genome-center-grp-high-qos*)
    echo 'genome-center-grp-high-qos|cpu=64,gres/gpu=0,mem=1T|' ;;
  *"show qos"*publicgrp-low-qos*)
    echo 'publicgrp-low-qos||' ;;
esac
"""


def _exe(path, text):
    with open(path, "w") as fh:
        fh.write(text)
    os.chmod(path, 0o755)


def _read(path):
    with open(path) as fh:
        return fh.read()


def _header(path, key):
    m = re.search(rf"^#SBATCH --{key}=(\S+)$", _read(path), re.M)
    return m.group(1) if m else None


def _hours(t):
    parts = [int(x) for x in t.split("-")[-1].split(":")]
    days = int(t.split("-")[0]) if "-" in t else 0
    return days * 24 + parts[0] + parts[1] / 60


class ChainBase(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.d = self._tmp.name
        self.bin = os.path.join(self.d, "bin")
        os.makedirs(self.bin)
        _exe(os.path.join(self.bin, "sacctmgr"), SACCTMGR)
        _exe(os.path.join(self.bin, "sinfo"), "#!/bin/sh\necho '0/5000/0/5000'\n")
        _exe(os.path.join(self.bin, "squeue"), "#!/bin/sh\necho 0\n")
        self.env = job_env(self.d, PATH=self.bin + os.pathsep + "/usr/bin:/bin",
                           HOME=self.d, USER="someone")

    def tearDown(self):
        self._tmp.cleanup()

    def chain(self, *extra, n=29, window=True, folder="", fasta_name="db.fasta", ok=True):
        base = os.path.join(self.d, folder) if folder else self.d
        raws = []
        for i in range(n):
            raws.append(os.path.join(base, "raw", f"f{i:02d}.d"))
            os.makedirs(raws[-1], exist_ok=True)
        cfg = os.path.join(base, "diann.cfg")
        with open(cfg, "w") as fh:
            fh.write("--qvalue 0.01\n--mass-acc 15\n--mass-acc-ms1 15\n"
                     + ("--window 7\n" if window else ""))
        fasta = os.path.join(base, fasta_name)
        with open(fasta, "w") as fh:
            fh.write(">sp|P1|X\nPEPTIDEK\n")
        out = os.path.join(base, "search out" if folder else "out")
        diann = os.path.join(self.d, "diann-linux")
        if not os.path.exists(diann):
            _exe(diann, "#!/bin/sh\nexit 0\n")
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "diann_parallel.py"),
                            "--diann", diann, "--raw", *raws, "--fasta", fasta, "--out", out,
                            "--cfg", cfg, *([] if not window else ["--no-probe-window"]), *extra],
                           cwd=base, capture_output=True, text=True, env=self.env, timeout=120)
        if ok:
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        return out, r, fasta


class ArrayTimeLimit(ChainBase):
    """HIGH. 2.10 halves the CPUs of a big array on genome-center-grp/high (16 -> 8, the
    documented `run_search.py --threads 16`) but keeps TIME_PER_FILE_HOURS = 2. A file that took
    > 61 min at 16 CPUs -- 6.7% of step-2 tasks on HIVE since July -- now hits TIMEOUT, and one
    timed-out task leaves steps 3-5 DependencyNeverSatisfied: the whole chain is lost."""

    def test_the_wall_clock_per_file_keeps_up_with_fewer_cpus(self):
        out, r, _ = self.chain("--threads-max", "16")
        for step in ("step2_firstpass", "step4_finalpass"):
            path = os.path.join(out, f"{step}.sbatch")
            cpus, hours = int(_header(path, "cpus-per-task")), _hours(_header(path, "time"))
            self.assertEqual(cpus, 8, step)                       # what 2.10 sizes it to
            # 2.9.1 gave each file 16 CPUs x 2 h; fewer CPUs must not shrink that budget
            self.assertGreaterEqual(cpus * hours, 16 * 2,
                                    f"{step}: {cpus} CPUs x {hours} h < 2.9.1's 16 x 2 h")


class ArrayMemory(ChainBase):
    """HIGH. MEM_PER_FILE_GB = 32 for steps 1b and 2 is below the non-cache memory real Core runs
    used on HIVE (lower bound 54-61 GB for four step-2 tasks on 2026-09-29..10-01; > 32 GB for 65
    step-2 tasks and 9 step-1b jobs since July). Each would be OUT_OF_MEMORY and take the chain
    with it. 64 GB costs no concurrency on genome-center-grp/high: 8 tasks x 64 GB = 512 GB, under
    the 1 TB per-user cap."""

    def test_steps_1b_and_2_keep_the_memory_real_runs_used(self):
        out, r, _ = self.chain("--threads-max", "16", window=False)
        for step in ("step1b_window", "step2_firstpass"):
            mem = _header(os.path.join(out, f"{step}.sbatch"), "mem")
            self.assertGreaterEqual(int(mem.rstrip("G")), 64, f"{step}: --mem={mem}")


class PathsWithSpaces(ChainBase):
    """MEDIUM (since before 2.9.1; 2.10's --beside makes it likelier). Core service folders hold
    spaces (2,557 directories within 3 levels of Data/lab/service, e.g. `<Lab Name>
    09182015`). The chain quotes the raw files but splices --out, --fasta, --lib and --temp
    unquoted: step 2's `sed -n ... <out>/file_list.txt` splits on the space, prints `no file for
    task 0` and the array fails; step 1 hands DIA-NN half a FASTA path. Either refuse such a path
    at generation or quote it; never write a chain that cannot run."""

    def test_a_folder_with_a_space_runs_or_is_refused(self):
        out, r, fasta = self.chain(n=3, folder="PROT_0000 Example lab", fasta_name="db v1.fasta",
                                   ok=False)
        if r.returncode != 0:                                     # refused: fine, if nothing ran
            self.assertFalse(os.path.exists(os.path.join(out, "submit.sh")))
            return
        seen = os.path.join(self.d, "argv.json")
        _exe(os.path.join(self.d, "diann-linux"), f"""#!{sys.executable}
import json, os, sys
a = sys.argv[1:]
json.dump(a, open({seen!r}, "w"))
t, f = a[a.index("--temp") + 1], a[a.index("--f") + 1]
open(os.path.join(t, os.path.splitext(os.path.basename(f))[0] + ".quant"), "w").write("q")
""")
        env = dict(self.env, SLURM_ARRAY_TASK_ID="0")
        env.pop("SLURM_JOB_ID", None)
        s2 = os.path.join(out, "step2_firstpass.sbatch")
        p = subprocess.run(["bash", s2], cwd=out, capture_output=True, text=True, env=env,
                           timeout=120)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
        argv = json.load(open(seen))
        self.assertEqual(argv[argv.index("--fasta") + 1], fasta)


class RetryRunsOnce(unittest.TestCase):
    """MEDIUM-LOW. node_fault.py retry has no lock: two retries started together (two sessions on
    the shared brettsp account, or a resumed session beside the original) both read the old
    jobs.txt, both resubmit the failed step and every dependant, and two copies of steps 3-5
    then write the same search folder (and post and stage twice). flock does not lock across
    HIVE's login nodes on /quobyte, so the lock has to be an mkdir in the search folder."""

    SACCT = r"""#!/usr/bin/env python3
import json, os, sys
st = json.load(open(os.environ["NF_STATE"]))
for row in st.get(sys.argv[sys.argv.index("-j") + 1], []):
    print("|".join(row))
"""
    SBATCH = r"""#!/usr/bin/env python3
import json, os, sys, time
time.sleep(1.5)
with open(os.environ["NF_LOG"], "a") as fh:
    fh.write(json.dumps(sys.argv[1:]) + "\n")
print(os.getpid())
"""

    def test_two_retries_at_once_submit_the_step_once(self):
        with tempfile.TemporaryDirectory() as d:
            b = os.path.join(d, "bin")
            os.makedirs(b)
            _exe(os.path.join(b, "sacct"), self.SACCT)
            _exe(os.path.join(b, "sbatch"), self.SBATCH)
            _exe(os.path.join(b, "scancel"), "#!/bin/sh\nexit 0\n")
            out = os.path.join(d, "s", "output", "search")
            os.makedirs(out)
            with open(os.path.join(out, "submit.sh"), "w") as fh:
                fh.write(f'#!/bin/bash\ncd "{out}"\n'
                         "jid1=$(sbatch --parsable step1_libpred.sbatch)\n"
                         "jid2=$(sbatch --parsable --dependency=afterok:$jid1 step2_firstpass.sbatch)\n"
                         "jid3=$(sbatch --parsable --dependency=afterok:$jid2 step3_assembly.sbatch)\n"
                         f'printf "%s\\n" $jid1 $jid2 $jid3 > "{out}/jobs.txt"\n')
            with open(os.path.join(out, "jobs.txt"), "w") as fh:
                fh.write("101\n102\n103\n")
            state = os.path.join(d, "state.json")
            with open(state, "w") as fh:
                json.dump({"101": [["101", "COMPLETED", "0:0", "n0"]],
                           "102": [["102", "FAILED", "75:0", "hive-dc-7-4-50"]],
                           "103": [["103", "PENDING", "0:0", "None assigned"]]}, fh)
            log = os.path.join(d, "sbatch.jsonl")
            env = job_env(d, PATH=b + os.pathsep + os.environ["PATH"], NF_STATE=state, NF_LOG=log)
            argv = [sys.executable, os.path.join(SCRIPTS, "node_fault.py"), "retry", "--out", out]
            procs = [subprocess.Popen(argv, env=env, stdout=subprocess.PIPE,
                                      stderr=subprocess.PIPE, text=True) for _ in range(2)]
            for p in procs:
                p.communicate(timeout=120)
            with open(log) as fh:
                calls = [json.loads(ln) for ln in fh]
        step2 = [c for c in calls if c[-1] == "step2_firstpass.sbatch"]
        self.assertEqual(len(step2), 1, f"step 2 resubmitted {len(step2)} times: {calls}")


class RefinalizeDigest(unittest.TestCase):
    """LOW. The re-finalize digest (record_run.analysis_digest) includes the over-cap zip's size
    text, so finalize run again only to add a podcast -- a few MB more in the zip -- reads as
    'analysis re-finalized, changed': a second master-log entry, activity row and Slack post."""

    def test_the_zip_growing_is_not_a_changed_analysis(self):
        import record_run

        def rec(size):
            why = (f"{size} GB (without .quant, the predicted library and FASTAs) is over the "
                   "2.0 GB cap -- the record points at it instead")
            return {"name": "2026-10-01_X", "prot": None,
                    "zip_copy": {"copied": False, "reason": why},
                    "analysis": {"zip": {"exists": True, "path": "/s/2026-10-01_X.zip",
                                         "reason": why},
                                 "de": {"significant_per_contrast": {"B-A": 12}}}}
        self.assertEqual(record_run.analysis_digest(rec("2.4")),
                         record_run.analysis_digest(rec("2.5")))


class LimsXlsxHostile(unittest.TestCase):
    """LOW. The stdlib .xlsx reader on hand-made workbooks."""

    RUNS = [f"20260101_TEST_60spd_S{i:02d}_A{i}_1_{100 + i}" for i in range(1, 11)]

    def _map(self, d, x):
        out = os.path.join(d, "conditions.csv")
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "collect_conditions.py"),
                            "--map", out, "--runs", ",".join(self.RUNS), "--from-file", x],
                           capture_output=True, text=True, timeout=60)
        return p, out

    def _patch(self, x, name, fn):
        with zipfile.ZipFile(x) as z:
            parts = {n: z.read(n) for n in z.namelist()}
        parts[name] = fn(parts[name])
        with zipfile.ZipFile(x, "w") as z:
            for n, v in parts.items():
                z.writestr(n, v)

    def test_a_bad_shared_string_index_is_a_clean_error(self):
        from test_conditions_lims_xlsx import make_xlsx
        with tempfile.TemporaryDirectory() as d:
            x = os.path.join(d, "s.xlsx")
            make_xlsx(x, [("samples", [["sample", "group"], ["S01", "A"]], False)])
            self._patch(x, "xl/worksheets/sheet1.xml",
                        lambda b: re.sub(rb'(<c r="B2" t="s"><v>)\d+', rb"\g<1>99", b))
            p, _ = self._map(d, x)
        self.assertNotIn("Traceback", p.stderr, p.stderr)
        self.assertIn("save it from Excel", p.stdout + p.stderr)

    def test_merged_group_cells_are_not_an_empty_group(self):
        from test_conditions_lims_xlsx import make_xlsx
        rows = [["sample", "group"]] + [[f"S{i + 1:02d}", "Ctrl" if i == 0 else
                                         "Trt" if i == 5 else None] for i in range(10)]
        with tempfile.TemporaryDirectory() as d:
            x = os.path.join(d, "merged.xlsx")
            make_xlsx(x, [("samples", rows, False)])
            self._patch(x, "xl/worksheets/sheet1.xml", lambda b: b.replace(
                b"</sheetData></worksheet>", b'</sheetData><mergeCells count="2"><mergeCell '
                b'ref="B2:B6"/><mergeCell ref="B7:B11"/></mergeCells></worksheet>'))
            p, out = self._map(d, x)
            res = json.loads(p.stdout)
        # Excel shows Ctrl on rows 2-6 and Trt on 7-11: read it so, or name the runs left
        # without a group -- never a group called "".
        self.assertNotIn("", res["groups"], res["groups"])

    def test_samples_on_a_later_sheet_are_pointed_to(self):
        from test_conditions_lims_xlsx import make_xlsx
        with tempfile.TemporaryDirectory() as d:
            x = os.path.join(d, "s.xlsx")
            make_xlsx(x, [("Instructions", [["Fill in the Samples tab"]], False),
                          ("Samples", [["sample", "group"]] + [
                              [f"S{i + 1:02d}", "A" if i < 5 else "B"] for i in range(10)],
                           False)])
            p, _ = self._map(d, x)
        # the first visible sheet is read; when it has no sample/group columns the message must
        # say which sheet was read and that the others exist (--sheet)
        self.assertIn("Samples", p.stdout + p.stderr)
        self.assertIn("--sheet", p.stdout + p.stderr)


if __name__ == "__main__":
    unittest.main()
