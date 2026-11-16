#!/usr/bin/env python3
"""The C2 documentation fixes of skill 2.14.2 (Ruling 82, fixer C): the docs tell the agent the
2.14.2 rules, in commands that run as written.

* docs M1 + N1: the laptop route's push-back loop carries report_provenance.json, and its pull
  loop stops on which-search's exit 3 instead of running on with nothing pulled.
* docs M2: `session.py status`, the finish command, is given in a form that runs (`--session`)
  in golden rule 9, step 12, step 12b and the laptop route's finish list.
* docs M3: references/methods.md and SKILL.md 9d describe 2.14.2's Methods: the Core LC defaults
  (column, 50 °C, emitter) only for a timsTOF .d, every tag make_methods.py prints, and that a tag
  is never replaced with a guess.
* docs M4: the finish gates are in the references staff look them up in: core-submissions.md's
  holds list and exit table (deliver, bioshare send, email-draft; --override-finish; raw-only and
  interim never gated), podcast.md and notifications.md.
* docs L2-L4: watcher.md re-generates a bigger or longer job with run_search.py's flags, never a
  hand-edited sbatch; no-de-route.md's step 11 passes --session, and a two-search benchmark
  delivers interim, then no-de; SKILL.md 6b shows every --var-mods word and --ids-only.
* statistician L1: de-analysis.md and build_maxlfq.R say the median measurement is done (worse
  than quantile), and give the benchmark's numbers: the quantile step itself does the harm, not
  normalising twice. (The refusal text: tests/test_normalization_check.py.)
* E2E F6: one rule for a no-DE DDA search, said the same in no-de-route.md and SKILL.md steps 3,
  4 and 7: identifications only (--ids-only) unless the user wants amounts; then it quantifies.
* no-de-route.md step 3 describes the no-DE analysis brief as it is (fixer A, C2 biologist H1 +
  statistician M3): its own body, the search's files, no DE demands, the n = 1 caveat.

Hermetic: SKILL.md and the references are read; the pull block runs in bash against a stub
hive_exec.sh in a temp dir.
"""
import ast
import os
import re
import shutil
import subprocess
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SKILL = os.path.dirname(HERE)
REFS = os.path.join(SKILL, "references")


def read(rel):
    with open(os.path.join(SKILL, rel), encoding="utf-8") as fh:
        return fh.read()


def flat(text):
    return " ".join(text.split())


def between(text, start, end):
    a = text.index(start)
    return text[a:text.index(end, a)]


LAPTOP_DOCS = ("SKILL.md", os.path.join("references", "core-submissions.md"))


# 2.15 (the slim router): where each step's text lives now -- the router names the section
def section(rel, head):
    """references/<rel> from the heading `head` to the next ## heading."""
    text = read(os.path.join("references", rel))
    a = text.index(head)
    nxt = text.find("\n## ", a + 3)
    return text[a:nxt if nxt >= 0 else len(text)]


def no_de_block(step):
    """What SKILL.md's no-DE pointer of step 3 or step 8 reaches. 2.15 fix round (Ruling 93, review
    LOW 2): no-de-route.md's copy of the old router blocks is gone; the rule is said once, in the
    route's step 1, and step 8's branch is the route's step 2."""
    sec = section("no-de-route.md", "## The route, in order")
    if step == 3:
        return sec[sec.index("1. **The search (step 7)."):sec.index("2. **Skip steps 8")]
    return sec[sec.index("2. **Skip steps 8"):]


def core_sync():
    return read(os.path.join("scripts", "core_sync.sh"))


class LaptopRoutePushAndPull(unittest.TestCase):
    """docs M1: the push-back loop left out report_provenance.json, so HIVE took every laptop
    report for one made before 2.14.2. N1: `for f in $(... which-search --paths)` ran zero times
    on exit 3 and the steps after it went on."""

    def test_the_push_back_loop_carries_the_reports_record(self):
        # 2.15: the push is core_sync.sh --put (the slim-down: its lists live there, once)
        put = re.search(r'PUT_OUTPUT="([^"]+)"', core_sync()).group(1).split()
        self.assertIn("output/report_provenance.json", put)
        self.assertIn("SKILL.md", LAPTOP_DOCS)
        self.assertIn("core_sync.sh --get", read("SKILL.md"))
        ref = flat(read(os.path.join("references", "core-submissions.md")))
        for before in ref.split("bash scripts/core_sync.sh --put")[:-1]:
            self.assertIn("report_provenance.json", before[-900:],
                          "the prose before each push says why the record goes back")

    def test_the_pull_loop_stops_on_exit_3(self):
        for rel in LAPTOP_DOCS:
            self.assertNotRegex(read(rel), r"for f in \$\(bash scripts/hive_exec\.sh \"[^\n]*which-search",
                                rel)
        # 2.15: core_sync.sh --get asks which-search and stops (tests/test_core_sync.py runs it)
        sync = core_sync()
        self.assertIn("which-search --dir", sync)
        self.assertIn('echo "STOP: $SEARCH_STOP"', sync)
        self.assertIn('--get "$S/$i" "$L/$(folder "$i")" </dev/null', sync,
                      "ssh inside the read loop must not eat it")

    def run_block(self, rel, which_exit):
        d = tempfile.mkdtemp()
        self.addCleanup(shutil.rmtree, d, True)
        os.makedirs(os.path.join(d, "scripts"))
        log = os.path.join(d, "gets.log")
        with open(os.path.join(d, "scripts", "hive_exec.sh"), "w") as fh:
            fh.write(f"""#!/bin/bash
if [ "$1" = "--get" ]; then
  cat >/dev/null                      # as ssh would: it reads whatever stdin it is given
  printf '%s|%s\\n' "$2" "$3" >> {log!r}; exit 0
fi
case "$*" in *which-search*)
  if [ {which_exit} -ne 0 ]; then echo "none chosen, nothing to pull: two searches" >&2; exit {which_exit}; fi
  printf '%s\\n' "output/search sage/search_provenance.json" "output/search sage/report.tsv" \\
                 "input/wf/workflow.manifest.json";;
esac
""")
        env = dict(os.environ, HOME=d, S="/quobyte/x/My Lab/session")
        p = subprocess.run(["bash", "-c", pull_block(read(rel))], cwd=d, env=env,
                           capture_output=True, text=True, timeout=60)
        gets = open(log).read().splitlines() if os.path.exists(log) else []
        return p, gets, d

    def test_the_pull_block_runs_as_written(self):
        """2.15: the documented pull is core_sync.sh --get; tests/test_core_sync.py runs it for
        this case (a lab folder and a search folder with a space; the search's record, report and
        workflow manifest land in the same names)."""
        want = 'bash scripts/core_sync.sh --get --session "$S" --local ~/core/PROT_0807/session'
        for rel in LAPTOP_DOCS:
            self.assertIn(want, flat(read(rel)), rel)
        return
        for rel in LAPTOP_DOCS:
            p, gets, d = self.run_block(rel, 0)
            self.assertEqual(p.returncode, 0, p.stderr)
            self.assertNotIn("STOP", p.stdout)
            base = os.path.join(d, "core", "PROT_0807", "session")
            self.assertEqual(gets, [
                f"/quobyte/x/My Lab/session/output/search sage/search_provenance.json|{base}/output/search sage/",
                f"/quobyte/x/My Lab/session/output/search sage/report.tsv|{base}/output/search sage/",
                f"/quobyte/x/My Lab/session/input/wf/workflow.manifest.json|{base}/input/wf/"], rel)
            self.assertTrue(os.path.isdir(os.path.join(base, "output", "search sage")))

    def test_no_search_chosen_says_stop_and_pulls_nothing(self):
        """2.15: core_sync.sh --get with none chosen: the rest pulled, STOP, exit 1
        (tests/test_core_sync.py test_none_chosen_pulls_the_rest_and_stops)."""
        self.assertIn("no search is chosen", core_sync())
        self.assertIn("choose one (SKILL.md step 12b), then pull again", core_sync())
        return
        for rel in LAPTOP_DOCS:
            p, gets, _d = self.run_block(rel, 3)
            self.assertIn("STOP (which-search exit 3)", p.stdout, rel)
            self.assertIn("none chosen", p.stderr, rel)
            self.assertEqual(gets, [], rel)


def status_flags():
    """The flags session.py's `status` subcommand takes, read from its argparse."""
    src = read(os.path.join("scripts", "session.py"))
    a = src.index('sub.add_parser("status"')
    block = src[a:src.index("sub.add_parser(", a + 10)]
    return set(re.findall(r'add_argument\("(--[a-z-]+)"', block))


class StatusIsRunnable(unittest.TestCase):
    """docs M2: golden rule 9 said "run `session.py status`" with no flags, and steps 12, 12b and
    the laptop finish list never called it; `status` takes --session where its siblings take
    --dir, so a copied `status --dir` was an argparse error."""

    def test_every_status_command_in_the_docs_uses_a_flag_status_takes(self):
        flags = status_flags()
        self.assertIn("--session", flags)
        texts = {"SKILL.md": read("SKILL.md")}
        for fn in os.listdir(REFS):
            if fn.endswith(".md"):
                texts[os.path.join("references", fn)] = read(os.path.join("references", fn))
        seen = 0
        for rel, text in texts.items():
            for m in re.finditer(r"session\.py status (--[a-z-]+)", text):
                seen += 1
                self.assertIn(m.group(1), flags, f"{rel}: session.py status {m.group(1)}")
        self.assertGreaterEqual(seen, 6)

    def test_golden_rule_9_gives_the_command(self):
        # 2.15: golden rule 16 in the router (one line), its text in outputs.md
        self.assertIn('outputs.md "Golden rule 16"', read("SKILL.md"))
        rule9 = flat(section("outputs.md", "## Golden rule 16"))
        self.assertIn("`python3 scripts/session.py status --session <session>`", rule9)
        self.assertIn('session.py status --session $(printf %q "$S")"`', rule9)
        self.assertIn("exit 3 lists each missing item", rule9)

    def test_step_12_and_12b_run_it(self):
        # 2.15: steps 12 and 12b are outputs.md "12." and core-submissions.md "12b."
        self.assertIn("session.py status --session <session>", flat(read("SKILL.md")))
        s12 = section("outputs.md", "## 12. Finalize the session")
        code = between(s12, "```", "\nfinalize refuses")
        self.assertIn("python3 scripts/session.py status --session <session>", code)
        self.assertIn("Run `session.py status --session <session>` and show its result", flat(s12))
        s12b = flat(section("core-submissions.md", "## 12b. Deliver to the collaborator"))
        self.assertIn("scripts/session.py status --session \"$S\"", s12b)

    def test_the_laptop_finish_lists_run_it_on_hive(self):
        want = 'bash scripts/hive_exec.sh "python3 ~/proteomics-pipeline/scripts/session.py status --session $(printf %q "$S")"'
        skill = between(section("core-submissions.md", "## 1c. The service run"),
                        "**Then on HIVE, finish the session:**",
                        "8. **Deliver, share, draft the email**")     # 2.15: SKILL.md 1c's text
        self.assertIn(want, skill)
        core = between(read(os.path.join("references", "core-submissions.md")),
                       "**Then finish the session on HIVE**", "## 6. `deliver`")
        self.assertIn(want, core)


def methods_tags():
    """Every `*_TAG = "[...]"` and DEF make_methods.py prints, read from its source."""
    src = read(os.path.join("scripts", "make_methods.py"))
    tags = dict(re.findall(r'^([A-Z_]*TAG|DEF) = "(\[[^"]+\])"', src, re.M))
    return tags


class MethodsDocs(unittest.TestCase):
    """docs M3: methods.md still said the PepSep column, 50 °C and the 20 µm emitter are printed
    for any run, named neither of 2.14.2's new tags, and asked only for the facility-default tag
    to be resolved -- so an agent "polishing the prose" could replace [vendor — confirm] with a
    guess."""

    def test_methods_md_names_every_tag_and_never_guess(self):
        tags = methods_tags()
        for name in ("DEF", "NR_TAG", "UNREADABLE_TAG", "VENDOR_TAG", "ACK_OUTSIDE_TAG",
                     "EVOTIP_TAG", "DDA_TAG"):
            self.assertIn(name, tags)
        ref = flat(read(os.path.join("references", "methods.md")))
        for name, tag in tags.items():
            self.assertIn(tag, ref, name)
        self.assertIn("[DDA/DIA — confirm]", ref)
        self.assertIn("never replaced with a guess", ref)
        self.assertNotIn("resolve each `[facility default — confirm]`, and keep", ref)

    def test_the_core_lc_defaults_are_for_a_timstof_d_only(self):
        ref = flat(read(os.path.join("references", "methods.md")))
        self.assertGreaterEqual(ref.lower().count("for a timstof `.d` only"), 3)   # column, 50 °C, emitter
        self.assertIn("any other run", ref)
        s9d = flat(section("methods.md", "## 9d. The Methods section"))    # 2.15
        self.assertGreaterEqual(s9d.lower().count("for a timstof `.d` only"), 2)
        self.assertNotIn("STAN's 50 °C and the 20 µm emitter are given as facility defaults, tagged;",
                         s9d)
        tags = methods_tags()
        for name in ("VENDOR_TAG", "ACK_OUTSIDE_TAG", "NR_TAG", "DEF"):
            self.assertIn(tags[name], s9d, name)
        self.assertIn("never replaced with a guess", s9d)


def code_constant(rel, name):
    """A literal constant assigned at the top level of a script, read with ast (not imported)."""
    tree = ast.parse(read(rel))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", None) == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError(f"{name} not found in {rel}")


# where each finish item lives, as the docs name it (report_files._ITEM_LABEL)
ITEM_FILE = {"report": "Analysis_Report.html", "methods": "output/methods.md",
             "reproducibility": "run_manifest.json", "finalized": "FINALIZED.json"}


def argparse_flags(rel):
    return set(re.findall(r'add_argument\((?:[A-Z_]+|"(--[a-z0-9-]+)")', read(rel))) - {""}


class DocLows(unittest.TestCase):
    def test_watcher_regenerates_never_hand_edits(self):
        """L2: a hand-edited sbatch runs at the new size while search_provenance.json still
        records job_mem: 64G (default) -- rule 1, and golden rule 6 (never hand-write the job)."""
        text = read(os.path.join("references", "watcher.md"))
        rows = {m.group(1): m.group(0) for m in
                re.finditer(r"^\| `(out_of_memory|timeout)` \|.*$", text, re.M)}
        self.assertEqual(set(rows), {"out_of_memory", "timeout"})
        self.assertNotIn("raise sbatch `--mem`", text)
        self.assertIn("never by editing the sbatch by hand", rows["out_of_memory"])
        self.assertIn("never by hand", rows["timeout"])
        self.assertIn("run_search.py --mem GB", rows["out_of_memory"])
        self.assertIn("run_search.py --time HOURS", rows["timeout"])
        flags = argparse_flags(os.path.join("scripts", "run_search.py"))
        for row in rows.values():
            for f in re.findall(r"`(?:run_search\.py )?(--[a-z-]+)", row):
                self.assertIn(f, flags, f)


    def test_step_6b_shows_var_mods_and_ids_only(self):
        """L4: the template showed [--var-mods ox] only: no deam, no none, no --ids-only, and
        nothing said a Sage DDA search now gets Ox(M) by default."""
        s6b = section("parameters.md", "## 6b. Estimating the search parameters")    # 2.15
        line = next(ln for ln in s6b.splitlines() if "[--var-mods " in ln)
        words = re.search(r"\[--var-mods ([^\] ]+)\]", line).group(1).split("|")
        rel = os.path.join("scripts", "estimate_params.py")
        known = set(code_constant(rel, "VAR_MODS")) | {code_constant(rel, "VAR_MODS_NONE")}
        self.assertEqual({w for x in words for w in x.split(",")}, known)
        self.assertIn(code_constant(rel, "IDS_ONLY_FLAG"), line)
        # 2.15 (lane F1): the FragPipe DDA preset shares the DDA default; the doc line is Sage's
        self.assertEqual(code_constant(rel, "VAR_MODS_DEFAULT")[("sage", "DDA")], ["ox"])
        text = flat(s6b)
        self.assertIn("Sage DDA searches Ox(M) by default", text)
        self.assertIn("any other word is refused", text)
        self.assertIn("refused with DIA-NN", text)
        # the F6 rule's exception, and where the flag is described (verifier nit on 4313d09d)
        note = text[text.index("**`--ids-only`**"):text.index("**`--labels` for a SILAC sample**")]
        self.assertIn("unless amounts are wanted", note)
        self.assertIn("`references/no-de-route.md`, step 1", note)
        # 2.15: this is parameters.md itself, so the pointer says "above"
        self.assertIn('"Sage for identifications only (`--ids-only`)" above', note)
        self.assertIn("## Sage for identifications only (`--ids-only`)",
                      read(os.path.join("references", "parameters.md")))


class OneSearchRuleForNoDe(unittest.TestCase):
    """E2E F6: SKILL.md step 7 called the quantifying route "the default", while step 3 and
    no-de-route.md offered --ids-only with "a search that quantifies is fine too"; neither said
    which to choose for a detection question."""

    def test_the_route_states_the_rule_and_when_to_quantify(self):
        ref = flat(read(os.path.join("references", "no-de-route.md")))
        step1 = ref[ref.index("1. **The search (step 7)."):ref.index("2. **Skip steps 8")]
        for want in ("One rule, both DDA engines: a no-DE DDA search is identifications only "
                     "unless the user", "wants amounts", "`--ids-only` at step 4",
                     "`estimate_params.py --ids-only` at step 6b too", "**Quantify instead** (no `--ids-only`)",
                     "abundance compared between runs or samples", "a DE that may follow",
                     "A detection question", "ask once", "`log_decision.py`",
                     "step 4 refuses it otherwise"):
            self.assertIn(want, step1, want)
        self.assertNotIn("A search that quantifies is fine too", ref)

    def test_skill_md_says_the_same_rule(self):
        # 2.15: SKILL.md's step-3 block, step 4 and step 7 live in no-de-route.md,
        # parameters.md "4." and search-engines.md "7. On a cluster"; the router points to them
        self.assertIn("no-de-route.md", read("SKILL.md"))
        s3 = flat(no_de_block(3))
        self.assertIn("a no-DE DDA search is identifications only unless the user", s3)
        self.assertIn("**Quantify instead** (no `--ids-only`) when the question needs amounts", s3)
        ref = read(os.path.join("references", "no-de-route.md"))   # said once (LOW 2)
        self.assertEqual(flat(ref).count("a no-DE DDA search is identifications only"), 1)
        s4 = flat(section("parameters.md", "## 4. Resolving the defaults"))
        self.assertIn("identifications only (`--ids-only`), or quantified because amounts are "
                      "wanted", s4)
        s7 = flat(section("search-engines.md", "## 7. On a cluster"))
        self.assertNotIn("**quantification** (the default)", s7)
        self.assertIn("the one rule is in `references/no-de-route.md`, step 1", s7)
        eng = flat(read(os.path.join("references", "search-engines.md")))
        self.assertNotIn("**Quantification** (the default)", eng)


def pull_block(text):
    """The `if P=$(... which-search --paths) ... fi` block of a doc, dedented."""
    lines = text.splitlines()
    i = next(n for n, ln in enumerate(lines) if ln.strip().startswith("if P=$(")
             and "which-search" in ln)
    j = next(n for n in range(i, len(lines)) if re.search(r"\bfi\s*$", lines[n]))
    pad = len(lines[i]) - len(lines[i].lstrip())
    return "\n".join(ln[pad:] for ln in lines[i:j + 1]) + "\n"


if __name__ == "__main__":
    unittest.main()
