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
  normalising twice. (The refusal text: tests/test_check_before_2142.py.)
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


class LaptopRoutePushAndPull(unittest.TestCase):
    """docs M1: the push-back loop left out report_provenance.json, so HIVE took every laptop
    report for one made before 2.14.2. N1: `for f in $(... which-search --paths)` ran zero times
    on exit 3 and the steps after it went on."""

    def test_the_push_back_loop_carries_the_reports_record(self):
        for rel in LAPTOP_DOCS:
            line = next(ln for ln in read(rel).splitlines()
                        if ln.strip().startswith("for f in AI_Analysis_Report.md"))
            self.assertIn("report_provenance.json", line, rel)
            self.assertIn("report_provenance.json", flat(read(rel)).split(line.strip())[0][-900:],
                          f"{rel}: the prose before the loop says why the record goes back")

    def test_the_pull_loop_stops_on_exit_3(self):
        for rel in LAPTOP_DOCS:
            text = read(rel)
            self.assertNotRegex(text, r"for f in \$\(bash scripts/hive_exec\.sh \"[^\n]*which-search",
                                rel)
            block = pull_block(text)
            self.assertIn("else echo \"STOP (which-search exit $?)", block, rel)
            self.assertIn("</dev/null", block, f"{rel}: ssh inside the read loop must not eat it")

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
        rule9 = flat(between(read("SKILL.md"), "9. **Stay on the standard path",
                             "## Audience: assume nothing is installed"))
        self.assertIn("`python3 scripts/session.py status --session <session>`", rule9)
        self.assertIn('session.py status --session $(printf %q "$S")"`', rule9)
        self.assertIn("exit 3 lists each missing item", rule9)

    def test_step_12_and_12b_run_it(self):
        text = read("SKILL.md")
        s12 = between(text, "\n### 12. Finalize the session", "\n### 12b. ")
        code = between(s12, "```", "\nfinalize refuses")
        self.assertIn("python3 scripts/session.py status --session <session>", code)
        self.assertIn("Run `session.py status --session <session>` and show its result", flat(s12))
        s12b = flat(between(text, "\n### 12b. ", "\n## Keeping the record of the analysis"))
        self.assertIn("scripts/session.py status --session \"$S\"", s12b)

    def test_the_laptop_finish_lists_run_it_on_hive(self):
        want = 'bash scripts/hive_exec.sh "python3 ~/proteomics-pipeline/scripts/session.py status --session $(printf %q "$S")"'
        skill = between(read("SKILL.md"), "**Then on HIVE, finish the session:**",
                        "8. **Deliver, share, draft the email**")
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
        s9d = flat(between(read("SKILL.md"), "\n### 9d. ", "\n### 9e. "))
        self.assertGreaterEqual(s9d.lower().count("for a timstof `.d` only"), 2)
        self.assertNotIn("STAN's 50 °C and the 20 µm emitter are given as facility defaults, tagged;",
                         s9d)
        tags = methods_tags()
        for name in ("VENDOR_TAG", "ACK_OUTSIDE_TAG", "NR_TAG", "DEF"):
            self.assertIn(tags[name], s9d, name)
        self.assertIn("never replaced with a guess", s9d)

    def test_methods_md_says_where_the_search_and_mode_come_from(self):
        ref = flat(read(os.path.join("references", "methods.md")))
        for want in ("`--session <S>`", "session.search_choice", "tagged blanks",
                     "input/acquisition.json", "results.json", "fragger.params",
                     "optimised automatically per run", "--ids-only", "first-pass report"):
            self.assertIn(want, ref, want)
        order = [ref.index(x) for x in ("from `--acquisition`", "the workflow manifest's "
                                        "`acquisition`", "the session's `input/acquisition.json`")]
        self.assertEqual(order, sorted(order))


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


class FinishGateDocs(unittest.TestCase):
    """docs M4: a staff agent whose deliver stopped with "held: the skill's finish steps are not
    done" found nothing in core-submissions.md's holds list or exit table; podcast.md and
    notifications.md said nothing of --override-finish or NOT FINALIZED."""

    RF = os.path.join("scripts", "report_files.py")

    def test_the_holds_list_names_the_gate_its_needs_and_the_way_past(self):
        ref = read(os.path.join("references", "core-submissions.md"))
        holds = flat(between(ref, "**Before writing anything, `deliver` refuses (exit 2)** when:",
                             "**Analysis mode**"))
        self.assertIn("held: the skill's finish steps are not done", holds)
        for item in code_constant(self.RF, "FINISH_NEEDS")["deliver"]:
            self.assertIn(ITEM_FILE[item], holds, item)
        for want in ("report_provenance.json", "stale", "finish_gaps", '--override-finish "<reason>"',
                     "finish_override", "[OVERRIDE]", "UNGATED_DELIVERY_MODES", "NOT FINALIZED",
                     "`bioshare send` and `email-draft`"):
            self.assertIn(want, holds, want)
        for mode in code_constant(self.RF, "UNGATED_DELIVERY_MODES"):
            self.assertIn(f"`--mode {mode}`", holds, mode)
        self.assertIn("never gated", holds)

    def test_the_exit_table_rows_name_it(self):
        ref = read(os.path.join("references", "core-submissions.md"))
        rows = {m.group(1): m.group(0) for m in
                re.finditer(r"^\| 2 \| `(deliver|email-draft|bioshare send)` \|.*$", ref, re.M)}
        self.assertEqual(set(rows), {"deliver", "email-draft", "bioshare send"})
        self.assertIn("held: the skill's finish steps are not done", rows["deliver"])
        for cmd, row in rows.items():
            self.assertIn("--override-finish", row, cmd)
        for cmd in ("email-draft", "bioshare send"):
            self.assertIn("not finalized", rows[cmd], cmd)

    def test_podcast_md_names_its_gate(self):
        text = flat(read(os.path.join("references", "podcast.md")))
        gate = text[text.index("**The finish gate**"):text.index("**Adding it later.**")]
        self.assertIn('FINISH_NEEDS["podcast"]', gate)
        for item in code_constant(self.RF, "FINISH_NEEDS")["podcast"]:
            self.assertIn(ITEM_FILE[item], gate, item)
        for want in ("`check`, `link`, `share` and `job` refuse", "only **warn**",
                     '--override-finish "<reason>"', "finish_override"):
            self.assertIn(want, gate, want)

    def test_notifications_md_names_the_announce_gate(self):
        text = flat(read(os.path.join("references", "notifications.md")))
        self.assertIn('FINISH_NEEDS["announce"]', text)
        self.assertEqual(code_constant(self.RF, "FINISH_NEEDS")["announce"], ("finalized",))
        for want in ("needs the finalize stamp", "stale", '--override-finish "<reason>"',
                     "**NOT FINALIZED**", "not gated"):
            self.assertIn(want, text, want)


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

    def test_no_de_route_step_11_and_two_searches(self):
        """L3: the step-11 shorthand passed neither --session nor --search-out, so the no-DE
        catalogue listed no search folder; the route's own benchmark (two searches) was not told
        that no-de then carries none."""
        ref = flat(read(os.path.join("references", "no-de-route.md")))
        self.assertNotIn("`make_report.py --out OUTPUT_FILES.md` (without", ref)
        step11 = ref[ref.index("**The file catalogue (step 11):**"):ref.index("7. **Finalize")]
        self.assertIn("make_report.py --out output/OUTPUT_FILES.md --session .", step11)
        self.assertIn("(without `--de-dir`)", step11)
        two = ref[ref.index("**Two or more searches**"):]
        self.assertLess(two.index("`--mode interim` first"), two.index("then `--mode no-de`"))
        for want in ("carries **no** search folder", "tagged blanks", "superseded_<name>", "owed"):
            self.assertIn(want, two, want)

    def test_step_6b_shows_var_mods_and_ids_only(self):
        """L4: the template showed [--var-mods ox] only: no deam, no none, no --ids-only, and
        nothing said a Sage DDA search now gets Ox(M) by default."""
        s6b = between(read("SKILL.md"), "\n### 6b. ", "\n### 6c. ")
        line = next(ln for ln in s6b.splitlines() if "[--var-mods " in ln)
        words = re.search(r"\[--var-mods ([^\] ]+)\]", line).group(1).split("|")
        rel = os.path.join("scripts", "estimate_params.py")
        known = set(code_constant(rel, "VAR_MODS")) | {code_constant(rel, "VAR_MODS_NONE")}
        self.assertEqual({w for x in words for w in x.split(",")}, known)
        self.assertIn(code_constant(rel, "IDS_ONLY_FLAG"), line)
        self.assertEqual(code_constant(rel, "VAR_MODS_DEFAULT"), {("sage", "DDA"): ["ox"]})
        text = flat(s6b)
        self.assertIn("Sage DDA searches Ox(M) by default", text)
        self.assertIn("any other word is refused", text)
        self.assertIn("refused with DIA-NN", text)
        # the F6 rule's exception, and where the flag is described (verifier nit on 4313d09d)
        note = text[text.index("**`--ids-only`**"):text.index("**`--labels` for a SILAC sample**")]
        self.assertIn("unless amounts are wanted", note)
        self.assertIn("`references/no-de-route.md`, step 1", note)
        self.assertIn('`references/parameters.md` "Sage for identifications only (`--ids-only`)"', note)
        self.assertIn("## Sage for identifications only (`--ids-only`)",
                      read(os.path.join("references", "parameters.md")))


class StatisticianDocLows(unittest.TestCase):
    """statistician L1: "double normalisation" is not what the benchmark shows (a quantile step
    alone, on un-normalised quantities, gave 46/72/69% wrong calls at 4/6/10 m/z), and median
    centring was measured -- worse than quantile -- not "being measured"."""

    def test_de_analysis_says_what_was_measured(self):
        text = flat(read(os.path.join("references", "de-analysis.md")))
        self.assertNotIn("being measured", text)
        block = text[text.index("**What the benchmark measured**"):]
        for want in ("5.5 / 46 / 72 / 69%", "4.8 / 60 / 76 / 74%", "3.4 / 5.5 / 18 / 30%",
                     "worse than quantile", "not normalising twice", "HYE", "Brett's decision"):
            self.assertIn(want, block, want)

    def test_build_maxlfq_comment_says_median_was_measured(self):
        src = read(os.path.join("scripts", "build_maxlfq.R"))
        self.assertNotIn("is measured on the benchmark first", src)
        self.assertIn("was worse than quantile", src)


class OneSearchRuleForNoDe(unittest.TestCase):
    """E2E F6: SKILL.md step 7 called the quantifying route "the default", while step 3 and
    no-de-route.md offered --ids-only with "a search that quantifies is fine too"; neither said
    which to choose for a detection question."""

    def test_the_route_states_the_rule_and_when_to_quantify(self):
        ref = flat(read(os.path.join("references", "no-de-route.md")))
        step1 = ref[ref.index("1. **The search (step 7)."):ref.index("2. **Skip steps 8")]
        for want in ("identifications only (`estimate_params.py --ids-only`, step 6b) unless the "
                     "user", "wants amounts", "**Quantify instead** (no `--ids-only`)",
                     "abundance compared between runs or samples", "a DE that may follow",
                     "A detection question", "ask once", "`log_decision.py`",
                     "refused with DIA-NN"):
            self.assertIn(want, step1, want)
        self.assertNotIn("A search that quantifies is fine too", ref)

    def test_skill_md_says_the_same_rule(self):
        text = read("SKILL.md")
        s3 = flat(text[text.index("**No differential expression?**"):text.index("\n### 3b. ")])
        self.assertIn("a Sage DDA search is identifications only (`estimate_params.py --ids-only`) "
                      "unless the user", s3)
        self.assertIn("then it quantifies", s3)
        s4 = flat(between(text, "\n### 4. ", "\n### 5. "))
        self.assertIn("identifications only (`--ids-only`), or quantified because amounts are "
                      "wanted", s4)
        s7 = flat(between(text, "\n### 7. ", "\n### 8. "))
        self.assertNotIn("**quantification** (the default)", s7)
        self.assertIn("the one rule is in `references/no-de-route.md`, step 1", s7)
        eng = flat(read(os.path.join("references", "search-engines.md")))
        self.assertNotIn("**Quantification** (the default)", eng)


class NoDeBriefDocs(unittest.TestCase):
    """no-de-route.md said only the brief's title changes for a no-DE session; since fixer A's
    analysis_prompt.no_de_body its body is written for it."""

    def step3(self):
        ref = flat(read(os.path.join("references", "no-de-route.md")))
        return ref[ref.index("3. **The report (step 9):**"):ref.index("4. **Methods (step 9d):**")]

    def test_the_route_describes_the_no_de_brief(self):
        s3 = self.step3()
        self.assertIn('"identification and QC (no differential expression)"', s3)
        self.assertIn("not only its title", s3)
        for want in ("`--session <S>`", "`identifications.tsv`", "`sage_identifications.json`",
                     "`report.parquet`", "`search_provenance.json`",
                     '"identifications only: no quantities"', "first-pass caveat",
                     "**It makes no DE demands:**", "significance rule", "normalisation check",
                     "DE Methods", "differences between them were not tested", "(n = 1)",
                     "run-to-run variation", "`session_docs.no_de_condition_caveat`",
                     '"unique to"'):
            self.assertIn(want, s3, want)

    def test_the_sections_named_are_the_briefs(self):
        """Once the no-DE body is in analysis_prompt.py, every `### ` heading it writes is named
        in step 3, and step 3 names no other."""
        src = read(os.path.join("scripts", "analysis_prompt.py"))
        if "def no_de_body(" not in src:
            self.skipTest("analysis_prompt.no_de_body is not on this branch yet (fixer A)")
        a = src.index("def no_de_body(")
        body = src[a:src.index("\ndef ", a + 10)]
        heads = re.findall(r'w\("### ([^"]+)"\)', body)
        self.assertGreaterEqual(len(heads), 6)
        s3 = self.step3()
        for h in heads:
            self.assertIn(h.replace(" (publication-ready)", ""), s3, h)


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
