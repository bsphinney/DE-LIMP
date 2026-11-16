"""Two readers of "how much memory does my SLURM allocation hold" must agree (cleanup audit E).

run_search.slurm_allocation_mb is the definition. podcast_local.job_resources is stdlib-only and
runs inside a container with no sibling modules, so it cannot import it and keeps its own read.
This ties the two: the same SLURM variables give the same memory.
"""
import os
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(os.path.dirname(HERE), "scripts"))
import podcast_local  # noqa: E402
import run_search  # noqa: E402


class SameAllocationSameMemory(unittest.TestCase):
    def check(self, env):
        mb = run_search.slurm_allocation_mb(dict(env))
        _, gb = podcast_local.job_resources(dict(env))
        self.assertEqual(gb, None if mb is None else mb / 1024.0, env)

    def test_mem_per_node(self):
        self.check({"SLURM_JOB_ID": "1", "SLURM_MEM_PER_NODE": "65536", "SLURM_CPUS_ON_NODE": "8"})

    def test_mem_per_cpu_times_the_cpus_on_the_node(self):
        self.check({"SLURM_JOB_ID": "1", "SLURM_MEM_PER_CPU": "4096", "SLURM_CPUS_ON_NODE": "8"})

    def test_mem_per_node_wins_over_mem_per_cpu(self):
        self.check({"SLURM_JOB_ID": "1", "SLURM_MEM_PER_NODE": "32768",
                    "SLURM_MEM_PER_CPU": "4096", "SLURM_CPUS_ON_NODE": "8"})

    def test_a_job_that_states_no_memory_is_unknown_to_both(self):
        self.check({"SLURM_JOB_ID": "1", "SLURM_CPUS_ON_NODE": "8"})


if __name__ == "__main__":
    unittest.main()
