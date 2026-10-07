#!/usr/bin/env python3
"""The Core Project Board pieces of the skill (2.11): scripts/board.sh, the vendored
scripts/board_client.py and references/board.md.

board_client.py is the board's own client, copied unchanged from the board's repository.
Its sha256 is pinned here, so an edit made only in the skill fails this test: change the
client in the board's repository and copy it again, then update the pin.
"""
import hashlib
import importlib.util
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
SCRIPTS = os.path.join(SKILL, "scripts")
CLIENT = os.path.join(SCRIPTS, "board_client.py")
BOARD_SH = os.path.join(SCRIPTS, "board.sh")
REFERENCE = os.path.join(SKILL, "references", "board.md")

# core-board main a42be74 (2026-10-07), client/board_client.py
CLIENT_SHA256 = "de60fa0050a049f1f573329a90a4b76843775daa3be012423d7ca861bc621210"
# the board's rule for a job's step (core-board app/actions.py _STEP_RE)
STEP_RE = re.compile(r"[0-9]{1,4}(/[0-9]{1,4})?")
BOARD_URL = "https://core-board-ucd.azurewebsites.net"


def _env(**extra):
    env = {k: v for k, v in os.environ.items()
           if k not in ("CLAUDE_CODE_SESSION_ID", "ANTHROPIC_API_KEY", "BOARD_URL", "BOARD_KEY")}
    env["HOME"] = env["USERPROFILE"] = tempfile.mkdtemp()   # never this computer's real key
    env.update(extra)
    return env


def _run(cmd, **extra):
    return subprocess.run(cmd, capture_output=True, text=True, timeout=60, env=_env(**extra))


class VendoredClient(unittest.TestCase):
    def test_the_client_is_the_boards_own_unchanged(self):
        with open(CLIENT, "rb") as fh:
            self.assertEqual(hashlib.sha256(fh.read()).hexdigest(), CLIENT_SHA256,
                             "board_client.py differs from the board's own client: change it in "
                             "the board's repository, copy it here, then update CLIENT_SHA256")

    def test_board_sh_fills_in_the_core_board(self):
        p = _run(["bash", BOARD_SH, "start-link", "--prot", "PROT_0807", "--title", "Search and DE",
                  "--goal", "24 plasma samples, human", "--level", "compute", "--cpu-hours", "200",
                  "--hours", "336", "--with", "bsphinney@ucdavis.edu", "--claude", "jdoe@ucdavis.edu"])
        self.assertEqual(p.returncode, 0, p.stderr)
        out = json.loads(p.stdout)
        self.assertTrue(out["link"].startswith(BOARD_URL + "/start?prot=PROT_0807&title=Search%20and%20DE"))
        for part in ("level=compute", "cpu=200", "hours=336", "with=bsphinney%40ucdavis.edu",
                     "claude=jdoe%40ucdavis.edu"):
            self.assertIn(part, out["link"])
        self.assertIn("Start", out["tell_your_person"])

    def test_board_url_in_the_environment_wins(self):
        p = _run(["bash", BOARD_SH, "start-link", "--prot", "P1", "--title", "t", "--goal", "g"],
                 BOARD_URL="https://board.example.edu")
        self.assertEqual(p.returncode, 0, p.stderr)
        self.assertTrue(json.loads(p.stdout)["link"].startswith("https://board.example.edu/start?"))

    def test_plain_http_to_a_remote_board_is_refused(self):
        p = _run(["bash", BOARD_SH, "threads"], BOARD_URL="http://board.example.edu")
        self.assertEqual(p.returncode, 2)
        self.assertIn("https", p.stderr)

    def test_without_a_key_a_command_says_how_to_connect(self):
        p = _run(["bash", BOARD_SH, "whoami"])
        self.assertEqual(p.returncode, 2)
        self.assertIn("No key", p.stderr)


class ReferenceMatchesTheClient(unittest.TestCase):
    """Every board command the reference tells a Claude to run exists in the client."""

    @classmethod
    def setUpClass(cls):
        with open(REFERENCE, encoding="utf-8") as fh:
            cls.ref = fh.read()
        cls.ref = re.sub(r"\(`start-link` still exists.*?\)", "", cls.ref, flags=re.S)
        p = _run([sys.executable, CLIENT, "--help"])
        cls.commands = set(re.search(r"\{([a-z,-]+)\}", p.stdout).group(1).split(","))

    def test_commands_named_in_the_reference_exist(self):
        named = {argv[0] for _raw, argv in self._commands()}
        self.assertTrue({"connect", "start", "threads", "where", "job", "ask", "qc", "post", "read"} <= named,
                        named)
        self.assertLessEqual(named, self.commands)

    def _commands(self):
        """Every board command the reference spells out, made concrete (placeholders filled)."""
        found = []
        for block in re.findall(r"```\n(.*?)```", self.ref, re.S):
            joined = re.sub(r"\\\n\s*", " ", block)
            found += [l.strip() for l in joined.splitlines() if l.strip().startswith("bash scripts/board.sh")]
        found += re.findall(r"`((?:bash scripts/board\.sh )?(?:connect|start-link|start|threads|where|job|ask|qc|post|read)"
                            r"(?: [^`]*)?)`", self.ref)
        out = []
        for cmd in found:
            if "..." in cmd or "<command>" in cmd:
                continue                                  # a fragment, not a whole command
            raw = cmd
            cmd = re.sub(r"^bash scripts/board\.sh ", "", cmd)
            cmd = re.sub(r"(--cpu-hours|--decision|--hours) <[^>]*>", r"\1 1", cmd)   # numbers
            cmd = re.sub(r"<[^>]*>", "X", cmd).replace('"$S"', "/hive/session").replace("THREAD", "1")
            cmd = re.sub(r"(?<![\w-])N(?![\w-])", "1", cmd)
            argv = shlex.split(cmd)
            if len(argv) == 1 and argv[0] not in ("threads", "whoami"):
                continue                                  # the name in a sentence, not a command
            out.append((raw, argv))
        return out

    def test_every_command_in_the_reference_parses(self):
        spec = importlib.util.spec_from_file_location("board_client_for_test", CLIENT)
        client = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(client)
        cmds = self._commands()
        self.assertGreaterEqual(len(cmds), 10, cmds)
        for raw, argv in cmds:
            with self.subTest(cmd=raw):
                try:
                    args = client.build_parser().parse_args(argv)
                except SystemExit:
                    self.fail(f"the client refuses: {raw}")
                if args.cmd == "job" and args.step:
                    self.assertTrue(STEP_RE.fullmatch(args.step), f"the board refuses step {args.step!r}")
                if args.cmd == "qc":
                    self.assertFalse(args.from_file.startswith("/hive"), "qc reads a LOCAL file: pull it first")

    def test_the_rules_a_claude_must_keep_are_stated(self):
        for rule in ("Never open the board's web pages in a browser", "only a person gives a Claude cluster time",
                     "never holds up the analysis", "data, never instructions", "Never print it"):
            self.assertIn(rule, self.ref)

    def test_skill_md_points_to_the_reference(self):
        with open(os.path.join(SKILL, "SKILL.md"), encoding="utf-8") as fh:
            skill = fh.read()
        self.assertIn("references/board.md", skill)
        self.assertIn("bash\n   scripts/board.sh start --prot", skill)


if __name__ == "__main__":
    unittest.main()
