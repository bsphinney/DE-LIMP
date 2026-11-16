#!/usr/bin/env python3
"""ONE definition of "Core staff" (2.15 Ruling 95.2).

Before 2.15 the router called anyone whose HIVE account can read /quobyte/proteomics-grp
(`check_access.sh`'s core_member) Core staff -- step 0a's question 3, steps 1-1c, the board,
CoreOmics -- while the scripts' staff-only steps (session.py core-client, the survey line, the
podcast's --core-run, the PRM follow-up) asked the admin-owned staff list (staff.core_staff). A
group member not on the list was asked for a submission and shown the board, then refused
core-client and never shown the survey line. Now Core staff is the list everywhere:
check_access.sh reports it as core_staff (tests/test_hive_remote_transport.py runs that), and
core_member says only that the group folder is readable. Text checks; stdlib only.
"""
import os
import re
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL_DIR = os.path.dirname(HERE)


def read(*parts):
    with open(os.path.join(SKILL_DIR, *parts), encoding="utf-8") as fh:
        return fh.read()


def flat(text):
    return " ".join(text.split())


def docs():
    yield "SKILL.md", read("SKILL.md")
    ref = os.path.join(SKILL_DIR, "references")
    for name in sorted(os.listdir(ref)):
        if name.endswith(".md"):
            yield f"references/{name}", read("references", name)


class CoreStaffIsTheList(unittest.TestCase):
    def test_no_staff_gate_is_group_read(self):
        """`core_member: true/false` never decides who is Core staff in the router or a
        reference (it may still say where the software comes from)."""
        for name, text in docs():
            for m in re.finditer(r"`core_member: (?:true|false)`", text):
                around = flat(text[max(0, m.start() - 200):m.end() + 200])
                self.assertNotRegex(around, r"Core staff only|not Core staff \(`core_member",
                                    f"{name}: {around}")

    def test_the_router_gates_name_core_staff(self):
        s = read("SKILL.md")
        rule12 = flat(s[s.index("12. **Staff-only stays staff-only.**"):s.index("13. **")])
        self.assertIn("Anyone who is not Core staff (`core_staff` not true) never hears of "
                      "CoreOmics", rule12)
        step0a = flat(s[s.index("### 0a."):s.index("### 0. ")])
        self.assertIn("says `core_staff: true` — right after that check", step0a)
        self.assertIn("relay `core_staff_detail`", step0a)
        self.assertIn("Pick the mode from `recommended_mode` + `core_member`", step0a)
        step1 = flat(s[s.index("### 1. "):s.index("### 1a.")])
        self.assertIn("(step 0a: yes to question 2 and `core_staff: true`)", step1)

    def test_access_md_defines_both_once(self):
        a = flat(read("references", "access.md"))
        self.assertIn('"Core staff" = this HIVE login is on the Core staff list (`staff.py`', a)
        self.assertIn("It is the ONE test for every staff-only step", a)
        self.assertIn("is not thereby Core staff", a)
        self.assertIn("run step 0, then `check_access.sh` again, and ask question 3 only then", a)
        self.assertEqual(a.count('"Core staff" = '), 1)
        self.assertEqual(a.count('"Core member" = '), 1)

    def test_board_and_core_submissions_follow(self):
        self.assertIn("(step 0a: yes to question 2 and `core_staff: true`)",
                      flat(read("references", "board.md")))
        core = flat(read("references", "core-submissions.md"))
        self.assertIn("Only Core staff run this (`core_staff: true`)", core)

    def test_the_scripts_staff_check_is_staff_py(self):
        """check_access.sh asks staff.py on HIVE (no second list or trust rule in bash)."""
        ca = read("scripts", "check_access.sh")
        self.assertIn("staff.py check --by", ca)
        self.assertNotIn("core_staff.txt", ca)
        self.assertIn('"core_staff": $CORE_STAFF', ca)


class TheSurveyNoteIsRelayed(unittest.TestCase):
    def test_step_9_relays_a_survey_note(self):
        """make_analysis_html.py / finalize print the survey NOTE to Core staff only, never into
        MANIFEST.txt, so step 12's MANIFEST relay does not cover it (SURV215 Edit 2)."""
        s = read("SKILL.md")
        nine = flat(s[s.index("### 9. Analyze the data"):s.index("### 9d.")])
        self.assertIn('Core staff: board.md §3 row "results ready"; relay a survey NOTE here or '
                      "at step 12 as printed.", nine)


if __name__ == "__main__":
    unittest.main()
