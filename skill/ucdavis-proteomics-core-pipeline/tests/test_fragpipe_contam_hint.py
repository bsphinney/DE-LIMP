#!/usr/bin/env python3
"""
The FragPipe `contam_` refusal names the database FragPipe can actually search (S9).

A FragPipe search on a database holding Philosopher's `contam_` entries is refused
(run_search.refuse_fragpipe_contam_database), and the refusal says how to build the skill's own
database instead. Its hint ended at fetch_fasta.py's search.fasta and pointed the workflow at it:
a target-only FASTA, which MSFragger refuses ("No decoys found in the FASTA file",
references/fragpipe-diatracer.md). 2.14.2's hint built the target+decoy `.fas` by hand; since 2.15
(review H3) the search's own job adds the decoys to a target-only FASTA with FragPipe's
Philosopher -- without --contam, which would bring the `contam_` entries back -- so the hint
searches the target FASTA as it is and leaves database.db-path to run_search.py.

Stdlib only; a three-entry FASTA.
"""
import os
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import run_search  # noqa: E402

DB = (">sp|P04406|G3P_HUMAN GAPDH OS=Homo sapiens OX=9606 GN=GAPDH PE=1 SV=3\nMGKVKVGVNGFGRIGR\n"
      ">contam_sp|P00761|TRYP_PIG Trypsin OS=Sus scrofa OX=9823 PE=1 SV=1\nIVGGYTCAANSIPYQV\n"
      ">rev_sp|P04406|G3P_HUMAN decoy\nRGIRGFGNVGVKVKGM\n")


class Hint(unittest.TestCase):
    def test_the_hint_searches_the_target_fasta_and_the_job_adds_decoys(self):
        with tempfile.TemporaryDirectory() as d:
            db = os.path.join(d, "fragpipe.fas")
            with open(db, "w") as fh:
                fh.write(DB)
            wf = os.path.join(d, "fp.workflow")
            with open(wf, "w") as fh:
                fh.write(f"database.decoy-tag=rev_\ndatabase.db-path={db}\n")
            with self.assertRaises(SystemExit) as cm:
                run_search.refuse_fragpipe_contam_database(wf, db)
            msg = str(cm.exception)
        self.assertIn("fetch_fasta.py fetch", msg)
        self.assertIn("--fasta ./search.fasta", msg)          # the target FASTA, as it is
        self.assertIn("job adds the decoys", msg)
        self.assertIn("database --custom", msg)
        self.assertIn("--prefix rev_", msg)
        self.assertIn("without --contam", msg)
        self.assertIn("database.db-path", msg)
        self.assertIn("No decoys found", msg)


if __name__ == "__main__":
    unittest.main()
