#!/usr/bin/env python3
"""No real SLURM job number in the published skill (public-repo hygiene, 2026-10-21).

The repository is public. A HIVE job number cited as evidence ("HIVE srun job <n>", "sbatch <n>",
"jobs <n>–82", a log path under slurm/job<n>) says nothing a reader can check and points at the
Core's cluster accounting, so the prose says "measured on HIVE" / "a test job" instead, and a test
that needs a job number as DATA (a fake sacct line, a log name, SLURM_JOB_ID) uses the reserved
fake range 9900001, 9900002, ...

What is caught: one of job, jobs, jobid, job_id, "job id", srun, sbatch or slurm/job (any case),
then up to 8 characters that are not digits, then a number of 7 to 9 digits.

What is NOT a job number, and so allowed:
  * the reserved fakes, 9900000-9900999;
  * a number that is one digit repeated (1111111, 99999999) or the ascending run 1234567,
    12345678, 123456789 (obvious placeholders);
  * a date, YYYYMMDD with a year 2000-2099 and a real month and day (20261015).

The patterns and the test strings are built from parts so this file does not match itself.
"""
import datetime
import os
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)

KEYWORDS = ("jo" + "b_id", "jo" + "bid", "jo" + "b id", "jo" + "bs", "jo" + "b",
            "sr" + "un", "sba" + "tch", "slurm/" + "jo" + "b")
JOB_NUMBER = re.compile("(?:" + "|".join(KEYWORDS) + r")[^\d\n]{0,8}(?<!\d)(\d{7,9})(?!\d)",
                        re.IGNORECASE)
FAKE_LOW, FAKE_HIGH = 99 * 10 ** 5, 99 * 10 ** 5 + 999          # 9900000-9900999
ASCENDING = {"1234567", "12345678", "123456789"}
TEXT = (".py", ".md", ".R", ".sh", ".json", ".txt", ".yaml", ".tsv", ".html", ".cfg", ".csv")


def is_date(n):
    if len(n) != 8 or not n.startswith("20"):
        return False
    try:
        datetime.date(int(n[:4]), int(n[4:6]), int(n[6:]))
    except ValueError:
        return False
    return True


def allowed(n):
    return (FAKE_LOW <= int(n) <= FAKE_HIGH or len(set(n)) == 1 or n in ASCENDING
            or is_date(n))


def job_numbers(line):
    """The job numbers in a line that are not on the allowed list."""
    return [m.group(1) for m in JOB_NUMBER.finditer(line) if not allowed(m.group(1))]


class NoJobIds(unittest.TestCase):
    def test_the_pattern_catches_the_real_shapes(self):
        for text in ("(HIVE " + "sr" + "un jo" + "b " + "2352" + "8991)",
                     "**Measured** (HIVE, 2026-10-01, " + "jo" + "bs " + "2427" + "6579–82)",
                     "/var/spool/" + "slurm/jo" + "b" + "2417" + "6598/slurm_script",
                     "(review " + "sr" + "un " + "2351" + "2013)",
                     "Submitted batch " + "jo" + "b " + "2399" + "0862",
                     "Jo" + "bId=" + "2417" + "6598",
                     '"jo' + 'b_id": "' + "2417" + '6598"'):
            with self.subTest(text=text):
                self.assertTrue(job_numbers(text), text)

    def test_the_allowed_numbers_pass(self):
        for text in ("Submitted batch " + "jo" + "b " + "99" + "00005",
                     "slurm/jo" + "b" + "99" + "00999/slurm_script",
                     "--" + "jo" + "bs " + "99" + "00009," + "99" + "00010",
                     "jo" + "b " + "1111" + "1111",
                     "jo" + "b " + "1234" + "567",
                     "the " + "jo" + "b of " + "2026" + "1015",
                     "measured on HIVE in a test " + "jo" + "b"):
            with self.subTest(text=text):
                self.assertEqual(job_numbers(text), [], text)

    def test_a_number_next_to_the_fake_range_or_an_impossible_date_is_caught(self):
        self.assertTrue(job_numbers("jo" + "b " + "99" + "01000"))
        self.assertTrue(job_numbers("jo" + "b " + "2026" + "1332"))       # month 13

    def test_none_in_the_skill(self):
        hits = []
        for dp, dns, fns in os.walk(SKILL):
            dns[:] = [d for d in dns if d not in ("__pycache__", ".pytest_cache")]
            for fn in fns:
                if not fn.endswith(TEXT):
                    continue
                path = os.path.join(dp, fn)
                with open(path, encoding="utf-8", errors="replace") as fh:
                    for i, line in enumerate(fh, 1):
                        for n in job_numbers(line):
                            hits.append(f"{os.path.relpath(path, SKILL)}:{i}: {n}")
        self.assertEqual(hits, [], "a real SLURM job number: in prose say \"measured on HIVE\" or "
                                   "\"a test job\"; as test data use 9900001, 9900002, ...")


if __name__ == "__main__":
    unittest.main()
