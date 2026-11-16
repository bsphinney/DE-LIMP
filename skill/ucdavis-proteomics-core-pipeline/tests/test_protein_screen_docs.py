#!/usr/bin/env python3
"""The protein screen's instructions: the router sends "show me the data" to its reference, the
reference's commands run as written, and its route table names every route the DE can write.

SKILL.md is the router (test_skill_structure.py): it keeps the trigger, the never-from-memory rule and
the step 8g pointer; references/protein-screen.md holds the commands. Hermetic: text files only.
"""
import os
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SK = os.path.dirname(HERE)

ENV = "source ~/.proteomics-pipeline/activate.sh"
SCRIPT = "~/proteomics-pipeline/scripts/make_screen_data.R"


def read(*parts):
    with open(os.path.join(SK, *parts), encoding="utf-8") as fh:
        return fh.read()


def section(md, start, stop):
    a = md.index(start)
    return md[a:md.index(stop, a + len(start))]


def flat(text):
    return " ".join(text.split())


REF = read("references", "protein-screen.md")
EIGHT_G = section(REF, "## 8g.", "## Which analyses part 1 covers")
SHOW_ME = section(REF, "## Show me the data for <protein>", "## What the folder holds")


class TheRouter(unittest.TestCase):
    def test_the_trigger_and_its_rule_are_in_the_router(self):
        md = read("SKILL.md")
        route = flat(section(md, "## Route the request", "## Flow"))
        self.assertIn('"Show me the data for <protein>"', route)
        self.assertRegex(route, r"never answer from memory → protein-screen\.md \"Show me the data\"")

    def test_8g_follows_8f_and_is_skipped_with_no_de(self):
        md = read("SKILL.md")
        head = md[md.index("### 8f. / 8g."):md.index("### 9. Analyze the data")]
        self.assertIn('8g: protein-screen.md "8g."', flat(head))
        self.assertIn("none wanted: skip 8–8d, 8f and 8g", md)
        self.assertIn("and 8g", read("references", "no-de-route.md"))


class StepEightGRunsAsWritten(unittest.TestCase):
    def test_the_job_file_route_comes_before_the_blocking_job(self):
        s = EIGHT_G
        self.assertLess(s.index("--constraint"), s.index("bash scripts/hive_exec.sh --job"))
        self.assertIn("%j", s)  # the header route names it
        self.assertNotIn("mkdir -p <session>/logs", s)
        self.assertRegex(flat(s), r"session\.py init[^.]*logs")

    def test_the_skills_r_by_full_path_with_its_environment(self):
        s = EIGHT_G
        self.assertIn(SCRIPT, s)
        self.assertNotIn("Rscript scripts/make_screen_data.R", s)
        self.assertNotIn("--wrap", s)
        sb = s[s.index("on_compute.py --header"):s.index("bash scripts/hive_exec.sh --job")]
        fb = s[s.index("bash scripts/hive_exec.sh --job"):]
        self.assertIn(ENV, sb)
        self.assertIn(ENV, fb)
        self.assertIn('--cpus "$SLURM_CPUS_PER_TASK"', sb)
        self.assertIn("bash scripts/hive_exec.sh 'sbatch --constraint=<that family> <S>/logs/screen_data.sbatch'", sb)
        self.assertIn("bash scripts/hive_exec.sh 'bash -s -- <S>/logs/screen_data.sbatch' < \"$U\"/lint-jobscript.sh", sb)
        self.assertLess(sb.index("lint-jobscript.sh"), sb.index("sbatch --constraint"))
        self.assertIn("bash scripts/hive_exec.sh --job --time 120 --mem 16G --cpus 8 "
                      "'" + ENV + " && Rscript " + SCRIPT + " --session <S> --cpus 8'", fb)
        self.assertRegex(flat(s), r"(?i)login shell[^.]*PATH")

    def test_a_refusal_is_said_not_worked_round(self):
        self.assertRegex(flat(EIGHT_G), r"(?i)exit 2 = refused.*do not work round it")


class ShowMeTheData(unittest.TestCase):
    def test_never_from_memory_and_never_rank_yourself(self):
        s = flat(SHOW_ME)
        self.assertRegex(s, r"(?i)never describe a protein's biology, and never write a number about it, from memory")
        self.assertRegex(s, r"(?i)quote only what is on the page")
        self.assertRegex(s, r"(?i)never pick, rank or list proteins yourself")
        self.assertRegex(s, r"(?i)quoted as it stands")
        self.assertIn("--top 10 --direction up", s)

    def test_the_builder_runs_on_a_compute_node_for_a_hive_session(self):
        # it hashes the search report (gigabytes): on HIVE a compute node, never the login node
        s = flat(SHOW_ME)
        self.assertRegex(s, r"(?i)never on the login node")
        self.assertIn("bash scripts/hive_exec.sh --job --time 60 --mem 4G 'python3 ~/proteomics-pipeline/scripts/"
                      "make_protein_screen.py --screen-data <S>/output/screen_data", s)
        self.assertIn("bash scripts/hive_exec.sh --get '<S>/output/protein_screens/", s)
        self.assertIn("on_compute.py", s)                # hive_local: the Flow section's split
        self.assertIn("ucd-hpc", s)

    def test_scripts_and_references_named_are_real(self):
        for ref in re.findall(r"references/([A-Za-z0-9_.-]+\.md)", REF):
            self.assertTrue(os.path.isfile(os.path.join(SK, "references", ref)), ref)
        for sc in re.findall(r"scripts/([A-Za-z0-9_.-]+\.(?:py|R|sh))", REF):
            self.assertTrue(os.path.isfile(os.path.join(SK, "scripts", sc)), sc)


class EveryRouteIsInTheTable(unittest.TestCase):
    """What works per route is said, never left out (repo rule 4): every pipeline_id run_de.R and
    build_maxlfq.R can write is a row of the reference's table, and only make_screen_data.R's
    SCREEN_PIPELINES get the whole page."""

    def test_the_table_names_every_pipeline_id(self):
        written = set()
        for rel in ("scripts/run_de.R", "scripts/build_maxlfq.R"):
            written |= set(re.findall(r'pipeline_id\s*=\s*"([a-z_]+)"', read(*rel.split("/"))))
        self.assertTrue({"dpc", "maxlfq", "engine_protein", "peptide_max"} <= written, written)
        table = section(REF, "| route (`pipeline_id`)", "\n\n")
        rows = {m for m in re.findall(r"^\| `([a-z_]+)`", table, re.M)}
        self.assertEqual(written - rows, set(), "a route the DE can write is missing from the table")

    def test_only_the_screens_pipelines_get_the_whole_page(self):
        code = read("scripts", "make_screen_data.R")
        m = re.search(r'^SCREEN_PIPELINES <- (.+)$', code, re.M)
        covered = set(re.findall(r'"([a-z_]+)"', m.group(1)))
        table = section(REF, "| route (`pipeline_id`)", "\n\n")
        whole = {r for r, cell in re.findall(r"^\| `([a-z_]+)` \|[^|]*\| ([^|]*) \|$", table, re.M)
                 if not cell.strip().startswith("refused")}
        self.assertEqual(whole, covered)


if __name__ == "__main__":
    unittest.main()
