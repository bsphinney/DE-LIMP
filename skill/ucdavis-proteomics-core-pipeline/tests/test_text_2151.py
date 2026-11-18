"""2.15.1 text lane: what the Methods and the docs say about what ran (audit lane 3).

Synthetic records only."""
import json
import os
import sys
import tempfile
import unittest
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import make_methods as mm        # noqa: E402
import ai_agent                  # noqa: E402


def write(path, obj):
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh)
    return path


class RadiantVersion(unittest.TestCase):
    """#6: the pinned 2.3.3 is a container tag; the Methods must not call it Radiant's version."""

    def paragraph(self, prov):
        with tempfile.TemporaryDirectory() as tmp:
            rec = mm.search_record(None, write(os.path.join(tmp, "sp.json"), prov), None)
        return mm.search_paragraph(rec)

    def test_a_container_tag_is_printed_as_a_container_tag(self):
        text = self.paragraph({"engine": "radiant", "version": "2.3.3"})
        self.assertIn("radiant-fulcrum container 2.3.3", text)
        self.assertIn("was not recorded", text)
        self.assertNotIn("processed with Radiant 2.3.3", text)

    def test_the_engines_own_version_is_printed_when_recorded(self):
        text = self.paragraph({"engine": "radiant", "version": "2.3.3",
                               "engine_version": "9.9.9"})
        self.assertIn("processed with Radiant 9.9.9", text)
        self.assertNotIn("2.3.3", text)

    def test_other_engines_are_unchanged(self):
        text = self.paragraph({"engine": "diann", "version": "2.7.0"})
        self.assertIn("processed with DIA-NN 2.7.0", text)


class EngineCitations(unittest.TestCase):
    """#14: every route cites what ran, from one table."""

    def paragraph(self, prov):
        with tempfile.TemporaryDirectory() as tmp:
            rec = mm.search_record(None, write(os.path.join(tmp, "sp.json"), prov), None)
        return mm.search_paragraph(rec)

    def test_diann_is_cited(self):
        self.assertIn("Demichev et al. 2020", self.paragraph({"engine": "diann", "version": "2.7.0"}))

    def test_sage_is_cited_from_the_one_table(self):
        from run_search import CITE
        self.assertIn(CITE["sage"], self.paragraph({"engine": "sage", "version": "0.14.7"}))

    def test_radiant_is_cited(self):
        self.assertIn("Just et al. 2026", self.paragraph({"engine": "radiant", "version": "2.3.3"}))


class AgentName(unittest.TestCase):
    """The AI that ran the skill is named from the record or the environment, never assumed."""

    def test_nothing_known_names_no_product(self):
        self.assertEqual(ai_agent.agent_label(env={}), "An AI assistant")

    def test_claude_code_is_named_when_it_says_so(self):
        self.assertEqual(ai_agent.agent_label(env={"CLAUDECODE": "1"}),
                         "Claude Code (Anthropic's AI coding assistant)")

    def test_a_model_only_when_the_environment_names_it(self):
        label = ai_agent.agent_label(env={"CLAUDECODE": "1", "ANTHROPIC_MODEL": "m-1"})
        self.assertIn("model m-1", label)

    def test_a_session_record_wins(self):
        with tempfile.TemporaryDirectory() as tmp:
            write(os.path.join(tmp, "session.json"), {"ai_agent": "Agent Z"})
            self.assertEqual(ai_agent.agent_label(tmp, env={"CLAUDECODE": "1"}), "Agent Z")


class CoreAcknowledgmentSignal(unittest.TestCase):
    def test_no_record_is_not_core(self):
        a = SimpleNamespace(acquired_at_core=False, session=None, submission=None)
        self.assertFalse(mm.core_acquisition_recorded(a))

    def test_the_flag_is_core(self):
        a = SimpleNamespace(acquired_at_core=True, session=None, submission=None)
        self.assertTrue(mm.core_acquisition_recorded(a))

    def test_a_session_answering_a_coreomics_submission_is_core(self):
        with tempfile.TemporaryDirectory() as tmp:
            os.makedirs(os.path.join(tmp, "input"))
            os.makedirs(os.path.join(tmp, "output"))
            write(os.path.join(tmp, "session.json"), {"coreomics": {"internal_id": "PROT_9001"}})
            a = SimpleNamespace(acquired_at_core=False, session=tmp, submission=None)
            self.assertTrue(mm.core_acquisition_recorded(a))

    def test_an_unreadable_session_fails_closed(self):
        a = SimpleNamespace(acquired_at_core=False, session="/nonexistent/x", submission=None)
        self.assertFalse(mm.core_acquisition_recorded(a))


if __name__ == "__main__":
    unittest.main()
