#!/usr/bin/env python3
"""`session.py finalize --report/--methods/--code-dir/--extra-dir` (2.14.2, Amendment 2 to Task
13, item 2): work done outside the session is folded in.

The beer re-analysis wrote its species-composition code and report in a second folder on the
laptop, outside the HIVE session, and delivered by hand; the session's record never saw them.
What these tests pin (custom_work.py, session.py finalize):
  * each source is copied to output/custom/<kind>/<name>, with output/custom/custom.json (where
    it came from, every file's sha256), the same checksums in the bundle's run_manifest.json, an
    [OK] line in MANIFEST.txt, and a line in AGENTS.md and the README;
  * a report the skill did not write goes in only under --override-finish, recorded, and the
    README calls it custom; the skill's own report from another folder needs no override;
  * refused before anything is written: a source inside output/custom/, one that holds the
    session, two of one kind with one name, a missing one, a file as a folder, too big;
  * a symlink inside a folder is not copied (MANIFEST.txt lists it); .git/__pycache__ stay out;
  * finalize again with the same folder replaces the copy; without the flag, the earlier items
    stay listed, are put back into a bundle provenance.py wrote anew, and a copy edited since is
    reported.
Stdlib only, no network.
"""
import json
import os
import shutil
import sys
import unittest
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import custom_work as cw                       # noqa: E402
import report_files as rf                      # noqa: E402
import test_finish_gates as tfg                # noqa: E402  (FinalizeBase: a no-DE DIA session)

read, write = tfg.read, tfg.write


class FoldIn(tfg.FinalizeBase):
    def setUp(self):
        super().setUp()
        self.finish()
        mac = os.path.join(self.root, "laptop", "beer_local")     # the second, local folder
        self.code = os.path.join(mac, "species_analysis")
        write(os.path.join(self.code, "species.py"), "print('species')\n")
        write(os.path.join(self.code, "R", "plots.R"), "plot(1)\n")
        write(os.path.join(self.code, ".git", "HEAD"), "ref: refs/heads/main\n")
        write(os.path.join(self.code, "__pycache__", "species.cpython-39.pyc"), "x")
        write(os.path.join(self.code, ".DS_Store"), "x")
        os.symlink(os.path.join(self.code, "species.py"), os.path.join(self.code, "latest.py"))
        self.extra = os.path.join(mac, "figures_v5")
        write(os.path.join(self.extra, "composition.png"), "png")
        self.report = write(os.path.join(mac, "report_v5.html"), "<html>pandoc v5</html>")
        self.methods = write(os.path.join(mac, "methods_v2.md"), "# Methods\n\nBy hand.\n")

    def custom(self, *rel):
        return os.path.join(self.out, "custom", *rel)

    def nothing_written(self, r, words):
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn(words, r.stderr)
        self.assertFalse(os.path.exists(self.custom()), "nothing copied")
        self.assertFalse(os.path.exists(self.p["manifest_txt"]), "nothing written")

    def test_code_and_extra_dirs_and_methods_are_copied_with_checksums_and_listed(self):
        r = self.finalize("--zip", "--code-dir", self.code, "--extra-dir", self.extra,
                          "--methods", self.methods)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        dest = self.custom("code", "species_analysis")
        self.assertEqual(read(os.path.join(dest, "species.py")), "print('species')\n")
        self.assertTrue(os.path.isfile(os.path.join(dest, "R", "plots.R")))
        for left in (".git", "__pycache__", ".DS_Store", "latest.py"):
            self.assertFalse(os.path.lexists(os.path.join(dest, left)), left)
        self.assertTrue(os.path.isfile(self.custom("extra", "figures_v5", "composition.png")))
        self.assertTrue(os.path.isfile(self.custom("methods", "methods_v2.md")))
        rec = json.loads(read(self.custom(cw.RECORD)))
        self.assertEqual(rec["schema"], cw.SCHEMA)
        by = {i["dest"]: i for i in rec["items"]}
        code = by["output/custom/code/species_analysis"]
        self.assertEqual(code["source"], os.path.realpath(self.code))
        self.assertEqual(code["files"], {"species.py": rf._sha256(os.path.join(dest, "species.py")),
                                         os.path.join("R", "plots.R"):
                                             rf._sha256(os.path.join(dest, "R", "plots.R"))})
        self.assertEqual(code["symlinks_not_copied"], ["latest.py"])
        self.assertEqual(code["left_out"], 3)
        bundle = json.loads(read(os.path.join(self.out, "reproducibility", "run_manifest.json")))
        self.assertEqual({i["dest"]: i["files"] for i in bundle["custom"]["items"]},
                         {d: i["files"] for d, i in by.items()})
        man = self.manifest()
        for d in by:
            self.assertRegex(man, r"\[OK\] +Folded-in work: " + d.replace("/", r"/") + " ")
        self.assertRegex(man, r"\[SKIPPED\] Folded-in work: output/custom/code/species_analysis "
                              r"\(symlinks\) +-- 1 symlink\(s\) not copied \(never followed\): "
                              r"latest\.py")
        self.assertRegex(man, r"\[OK\] +Folded-in work in the bundle")
        readme, agents = read(self.p["readme"]), read(os.path.join(self.sd, "AGENTS.md"))
        for doc in (readme, agents):
            self.assertIn("output/custom/code/species_analysis", doc)
            self.assertIn("output/custom/extra/figures_v5", doc)
            self.assertIn("output/custom/methods/methods_v2.md", doc)
            self.assertIn(os.path.realpath(self.code), doc)
        names = zipfile.ZipFile(self.sd + ".zip").namelist()
        self.assertTrue(any(n.endswith("output/custom/code/species_analysis/species.py")
                            for n in names), names)
        st = json.loads(read(self.stamp_path()))
        self.assertEqual(st["outputs"]["output/custom/custom.json"],
                         rf._sha256(self.custom(cw.RECORD)))

    def test_a_custom_report_goes_in_only_under_override_and_is_called_custom(self):
        self.nothing_written(self.finalize("--report", self.report),
                             "not written by make_analysis_html.py")
        why = "the beer species report, made by hand"
        r = self.finalize("--report", self.report, "--override-finish", why)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertEqual(read(self.custom("report", "report_v5.html")), "<html>pandoc v5</html>")
        man = self.manifest()
        self.assertRegex(man, r"\[OVERRIDE\] Finish check \(session\.py finalize\) +-- " + why
                              + r" -- went ahead without: the custom report")
        item = json.loads(read(self.custom(cw.RECORD)))["items"][0]
        self.assertEqual((item["skill_report"], item["override_reason"]), (False, why))
        readme = read(self.p["readme"])
        self.assertIn("output/custom/report/report_v5.html", readme)
        self.assertIn("custom: not written by the skill", readme)
        self.assertIn(why, readme)
        st = json.loads(read(self.stamp_path()))
        self.assertEqual(st["finish_override"]["reason"], why)
        self.assertTrue(any("the custom report" in g for g in st["finish_override"]["missing"]))

    def test_the_skills_own_report_from_another_folder_needs_no_override(self):
        other = os.path.join(self.root, "laptop", "beer_local", "output")
        html = tfg.skill_report(other)
        r = self.finalize("--report", html)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        copy = self.custom("report", "Analysis_Report.html")
        self.assertEqual(rf.skill_report_state(copy)[0], "skill", "its record came with it")
        self.assertIn("the skill's own report", read(self.p["readme"]))

    def test_refused_before_anything_is_written(self):
        os.makedirs(self.custom("code", "x"))
        for args, words in (
                (("--code-dir", self.custom("code", "x")), "inside output/custom/"),
                (("--extra-dir", self.root), "it holds the session itself"),
                (("--extra-dir", self.sd), "it holds the session itself"),
                (("--code-dir", self.code, "--code-dir",
                  shutil.copytree(self.code, os.path.join(self.root, "b", "species_analysis"),
                                  symlinks=True)),
                 "the same name (species_analysis)"),
                (("--code-dir", os.path.join(self.root, "nope")), "not found"),
                (("--code-dir", self.methods), "not a folder"),
                (("--methods", self.code), "not a file"),
                (("--extra-dir", self.extra, "--custom-max-gb", "0.0000000001"),
                 "--custom-max-gb")):
            with self.subTest(args=args):
                shutil.rmtree(self.custom("code"), ignore_errors=True)
                if args[0] == "--code-dir" and "inside" in words:
                    os.makedirs(self.custom("code", "x"))
                r = self.finalize(*args)
                self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
                self.assertIn(words, r.stderr)
                self.assertFalse(os.path.exists(self.p["manifest_txt"]), "nothing written")
                self.assertFalse(os.path.exists(self.custom(cw.RECORD)))

    def test_finalize_again_replaces_the_copy_and_keeps_earlier_items(self):
        self.assertEqual(self.finalize("--code-dir", self.code, "--extra-dir",
                                       self.extra).returncode, 0)
        write(os.path.join(self.code, "species.py"), "print('v2')\n")
        os.remove(os.path.join(self.code, "R", "plots.R"))
        r = self.finalize("--code-dir", self.code)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        dest = self.custom("code", "species_analysis")
        self.assertEqual(read(os.path.join(dest, "species.py")), "print('v2')\n")
        self.assertFalse(os.path.exists(os.path.join(dest, "R", "plots.R")), "a fresh copy")
        by = {i["dest"]: i for i in json.loads(read(self.custom(cw.RECORD)))["items"]}
        self.assertEqual(sorted(by), ["output/custom/code/species_analysis",
                                      "output/custom/extra/figures_v5"])
        self.assertEqual(list(by["output/custom/code/species_analysis"]["files"]), ["species.py"])
        # provenance.py run again (step 10) writes run_manifest.json anew; a plain finalize puts
        # the checksums back, and still lists both items
        bundle(self.out)
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        rm = json.loads(read(os.path.join(self.out, "reproducibility", "run_manifest.json")))
        self.assertEqual(len(rm["custom"]["items"]), 2)
        self.assertIn("output/custom/extra/figures_v5", read(self.p["readme"]))

    def test_a_copy_edited_after_it_was_folded_in_is_reported(self):
        self.assertEqual(self.finalize("--code-dir", self.code).returncode, 0)
        write(self.custom("code", "species_analysis", "species.py"), "edited in place\n")
        r = self.finalize()
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertRegex(self.manifest(), r"\[SKIPPED\] Folded-in work: output/custom/code/"
                                          r"species_analysis +-- changed since it was folded in "
                                          r"\(1 file\(s\): species\.py\)")

    def test_no_bundle_keeps_the_checksums_in_the_record_and_says_so(self):
        # finalize builds a missing bundle (L5b review H1): none here only when it cannot -- no
        # record of the raw files -- and finalize goes ahead under an override
        os.remove(os.path.join(self.out, "reproducibility", "run_manifest.json"))
        os.remove(self.p["raw_list"])
        write(os.path.join(self.out, "methods.md"), "# Methods\n\nWritten at step 9d.\n")
        r = self.finalize("--code-dir", self.code, "--override-finish", "a test of the record")
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        self.assertRegex(self.manifest(), r"\[SKIPPED\] Folded-in work in the bundle +-- no "
                                          r"output/reproducibility/run_manifest\.json")
        self.assertTrue(os.path.isfile(self.custom(cw.RECORD)))


class MethodsGivenButNotMade(tfg.FinalizeBase):
    """L5b review L5: `finalize --methods my.md` folds the file in, but the skill's
    output/methods.md is still required -- and with no raw-file list finalize could not write it,
    and refused ("NOT FINALIZED ... output/methods.md: no raw file list") with no way forward."""

    def test_the_refusal_says_how_to_finish(self):
        os.remove(self.p["raw_list"])
        shutil.rmtree(self.p["search_out"])
        mine = write(os.path.join(self.root, "my_methods.md"), "# Methods\nWritten by hand.\n")
        self.finish()
        r = self.finalize("--methods", mine)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("NOT FINALIZED", r.stderr)
        self.assertTrue(os.path.isfile(os.path.join(self.out, "custom", "methods",
                                                    "my_methods.md")))
        tail = r.stderr.split("NOT FINALIZED", 1)[1]
        self.assertIn(f"cp {mine} {os.path.join(self.out, 'methods.md')}", tail)
        self.assertIn("input/raw_files.txt", tail)
        self.assertIn("--override-finish", tail)


class SecretsLeftOut(tfg.FinalizeBase):
    """L5b review M3: `--code-dir`/`--extra-dir` copied a folder whole -- a `.env` holding an API
    key, a `.venv`, `node_modules` -- into output/custom/, which the zip carries to collaborators.
    By default secrets, keys and environments stay out (custom_work.LEFT_OUT_*, plus any small
    file the skill's secret patterns -- notify_slack.contains_secret -- find a secret in), and a
    [NOTE] in MANIFEST.txt names them."""

    def setUp(self):
        super().setUp()
        self.finish()
        self.code = os.path.join(self.root, "beer_local", "analysis")
        write(os.path.join(self.code, "species.py"), "print(1)\n")
        write(os.path.join(self.code, "renv", "activate.R"), "# renv\n")
        self.secret = "GEMINI_API_KEY=" + "AIza" + "B" * 30       # built here, never a real key
        self.left = {".env": self.secret + "\n", ".env.local": "X=1\n", ".Renviron": "Y=2\n",
                     "hf.token": "abc\n", "server.key": "k\n", "cert.pem": "p\n",
                     "id_ed25519": "k\n", "credentials.json": "{}\n",
                     ".venv/lib/site.py": "x\n", "venv/bin/python": "x\n", "env/x.py": "x\n",
                     "node_modules/a/index.js": "x\n", "renv/library/limma/DESCRIPTION": "x\n",
                     ".Rproj.user/state": "x\n", "__pycache__/a.pyc": "x\n", ".git/HEAD": "x\n",
                     "notes.md": "the key is " + "ghp_" + "C" * 36 + "\n"}
        for rel, text in self.left.items():
            write(os.path.join(self.code, *rel.split("/")), text)

    def test_secrets_and_environments_stay_out_and_are_named(self):
        r = self.finalize("--zip", "--code-dir", self.code)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        dest = os.path.join(self.out, "custom", "code", "analysis")
        got = sorted(os.path.relpath(os.path.join(d, f), dest).replace(os.sep, "/")
                     for d, _, fs in os.walk(dest) for f in fs)
        self.assertEqual(got, ["renv/activate.R", "species.py"])
        names = zipfile.ZipFile(self.sd + ".zip").namelist()
        self.assertFalse([n for n in names if n.endswith("/.env") or "node_modules" in n
                          or "/.venv/" in n], names)
        for n in names:
            if n.endswith((".py", ".R", ".md")) and "/custom/" in n:
                self.assertNotIn("AIza", zipfile.ZipFile(self.sd + ".zip").read(n).decode())
        man = self.manifest()
        note = [ln for ln in man.splitlines() if ln.startswith("[NOTE]")
                and "output/custom/code/analysis" in ln]
        self.assertEqual(len(note), 1, man)
        for word in (".env", ".venv/", "node_modules/", "renv/library/", "id_ed25519",
                     "notes.md (holds a secret)"):
            self.assertIn(word, note[0])
        rec = json.loads(read(os.path.join(self.out, "custom", cw.RECORD)))
        item = rec["items"][0]
        self.assertIn(".env", item["left_out_names"])

    def test_the_plan_says_what_it_leaves_out_before_anything_is_written(self):
        items = cw.plan(self.sd, self.out, code_dirs=[self.code])
        self.assertEqual(sorted(rel for rel, _, _ in items[0]["files"]),
                         ["renv/activate.R".replace("/", os.sep), "species.py"])
        self.assertIn("notes.md (holds a secret)", items[0]["left_out_names"])


class CredentialsNeverTravel(tfg.FinalizeBase):
    """C2 security H1 + M2: a fold-in kept the skill's own key files (.coreomics_token, and the
    STAN service-account secret .pgfarm_token -- no `*.token` rule matched either), a home
    folder's .ssh/.aws/.config/..., .netrc/.pgpass, and ID_RSA/Credentials.json (matching was
    case-sensitive); a file over 1 MiB was never read for secrets and nothing said so; and
    `--extra-dir ~` or `--extra-dir <S>/logs` was taken in, so the zip carried a home folder's
    keys or the Core-internal conversation the zip leaves out by path. Now the names and folders
    stay out (case-insensitively), every unread big file is named in a [NOTE], and a home folder
    (or one holding it) and a session's logs/, input/ or output/ are refused."""

    LEFT = (".coreomics_token", ".pgfarm_token", "hf_token.txt", ".netrc", ".pgpass",
            ".git-credentials", "ID_RSA", "Credentials.json", "client.P12", ".bash_history",
            ".Rhistory", ".ssh/known_hosts", ".ssh/config", ".aws/credentials",
            ".azure/accessTokens.json", ".config/ucdavis-proteomics/hive.env", ".gnupg/pubring",
            ".kube/config", ".docker/config.json", ".claude/settings.json")

    def setUp(self):
        super().setUp()
        self.finish()
        self.code = os.path.join(self.root, "work", "analysis")
        write(os.path.join(self.code, "keep.py"), "print(1)\n")
        for rel in self.LEFT:
            write(os.path.join(self.code, *rel.split("/")), "not-a-real-secret\n")

    def test_the_skills_own_keys_and_credential_folders_stay_out(self):
        items = cw.plan(self.sd, self.out, code_dirs=[self.code])
        self.assertEqual([rel for rel, _, _ in items[0]["files"]], ["keep.py"])
        named = items[0]["left_out_names"]
        for rel in self.LEFT:
            top = rel.split("/")[0]
            self.assertIn(top + "/" if "/" in rel else rel, named, rel)

    def test_the_key_names_are_the_ones_the_skill_writes(self):
        """the left-out rules take the skill's own key files by their one definition each"""
        import core_submission
        import ht_manifest
        names = {core_submission.KEY_NAME}
        names |= {os.path.basename(t) for t in ht_manifest.TOKEN_CANDIDATES}
        for n in sorted(names):
            self.assertTrue(cw._left_out_file(n), n)

    def test_the_finalized_zip_carries_none_of_them(self):
        r = self.finalize("--zip", "--code-dir", self.code)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        names = zipfile.ZipFile(self.sd + ".zip").namelist()
        custom = [n for n in names if "/output/custom/code/" in n]
        self.assertEqual([n.rsplit("/", 1)[1] for n in custom if not n.endswith("/")],
                         ["keep.py"], custom)

    def test_every_file_too_big_to_read_is_named_as_not_read(self):
        big = 1024 * 1024 + 1
        for rel in ("data/big_one.bin", "notebook.ipynb"):
            write(os.path.join(self.code, *rel.split("/")), "x" * big)
        items = cw.plan(self.sd, self.out, code_dirs=[self.code])
        self.assertEqual(sorted(items[0]["not_scanned_names"]),
                         ["data/big_one.bin", "notebook.ipynb"])
        r = self.finalize("--code-dir", self.code)
        self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
        note = [ln for ln in self.manifest().splitlines() if ln.startswith("[NOTE]")
                and "not read for secrets" in ln]
        self.assertEqual(len(note), 1, self.manifest())
        for words in ("data/big_one.bin", "notebook.ipynb", "WITHOUT being read for secrets"):
            self.assertIn(words, note[0])

    def test_a_home_folder_or_one_holding_it_is_refused(self):
        home = os.path.join(self.root, "homes", "someone")
        write(os.path.join(home, ".coreomics_token"), "not-a-real-secret\n")
        env = dict(self.env, HOME=home)
        for given, words in ((home, "your home folder"),
                             (os.path.dirname(home), "holds your home folder")):
            with self.subTest(given=given):
                r = self.finalize("--extra-dir", given, env=env)
                self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
                self.assertIn(words, r.stderr)
                self.assertFalse(os.path.exists(os.path.join(self.out, "custom")))
                self.assertFalse(os.path.exists(self.p["manifest_txt"]), "nothing written")

    def test_a_sessions_logs_input_or_output_folder_is_refused(self):
        logs = os.path.join(self.sd, "logs")
        write(os.path.join(logs, "conversation.md"), "the conversation\n")
        other = os.path.join(self.root, "other_session")
        write(os.path.join(other, "session.json"), "{}\n")
        write(os.path.join(other, "logs", "decisions.md"), "decisions\n")
        write(os.path.join(other, "output", "tables", "x.tsv"), "x\n")
        for given, words in ((logs, "this session's logs/"),
                             (os.path.join(self.sd, "input"), "this session's input/"),
                             (os.path.join(self.out, "search"), "this session's output/"),
                             (os.path.join(other, "logs"), "logs/ folder -- Core-internal"),
                             (os.path.join(other, "output", "tables"), "'s output/")):
            with self.subTest(given=given):
                with self.assertRaises(cw.Refused) as cm:
                    cw.plan(self.sd, self.out, extra_dirs=[given])
                self.assertIn(words, str(cm.exception))
        # single files too (verify residual): --methods "$S/logs/decisions.md", --report from
        # logs/ are the Core-internal records the zip leaves out by path
        write(os.path.join(logs, "decisions.md"), "decisions\n")
        write(os.path.join(logs, "review.html"), "<html>review</html>")
        for kw in ({"methods": os.path.join(logs, "decisions.md")},
                   {"report": os.path.join(logs, "review.html")}):
            with self.subTest(**kw):
                with self.assertRaises(cw.Refused) as cm:
                    cw.plan(self.sd, self.out, **kw)
                self.assertIn("this session's logs/", str(cm.exception))
        r = self.finalize("--methods", os.path.join(logs, "decisions.md"))
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("this session's logs/", r.stderr)
        r = self.finalize("--extra-dir", logs)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("this session's logs/", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "custom")))

    def test_a_methods_file_holding_a_secret_is_refused(self):
        mine = write(os.path.join(self.root, "methods_by_hand.md"),
                     "# Methods\n\nkey " + "ghp_" + "C" * 36 + "\n")       # built here, not real
        with self.assertRaises(cw.Refused) as cm:
            cw.plan(self.sd, self.out, methods=mine)
        self.assertIn("the skill's secret patterns find a secret", str(cm.exception))


class CopiesAreWholeOrNotAtAll(tfg.FinalizeBase):
    """C2 errors_ux M4 (and security L4): apply() deleted the earlier copy, then copied file by
    file into the final folder, so a copy that failed part way (an unreadable file, a full disk)
    left half a folder in output/custom/ -- and in the zip -- with the record still naming the old
    checksums; plan() never refused an unreadable file; and a file swapped for a link after
    plan() was copied through the link. Now an unreadable file is refused before anything is
    written, a copy is made whole in output/custom/.partial/ and then moved into place, and a
    failure leaves the earlier copy and the record as they were."""

    def setUp(self):
        super().setUp()
        self.finish()
        self.code = os.path.join(self.root, "work", "species_analysis")
        write(os.path.join(self.code, "species.py"), "print('v1')\n")
        write(os.path.join(self.code, "R", "plots.R"), "plot(1)\n")
        self.dest = os.path.join(self.out, "custom", "code", "species_analysis")

    def fold_in(self):
        return cw.apply(self.sd, self.out, cw.plan(self.sd, self.out, code_dirs=[self.code]))

    def test_an_unreadable_file_is_refused_before_anything_is_written(self):
        if hasattr(os, "geteuid") and os.geteuid() == 0:
            self.skipTest("root reads any file")
        locked = os.path.join(self.code, "R", "plots.R")
        # the mode is put back here, not by addCleanup: FinalizeBase.tearDown removes the
        # temporary folder first, and a cleanup's chmod then fails on the file it removed
        os.chmod(locked, 0)
        try:
            r = self.finalize("--code-dir", self.code)
        finally:
            os.chmod(locked, 0o644)
        self.assertEqual(r.returncode, 2, r.stdout + r.stderr)
        self.assertIn("cannot be read", r.stderr)
        self.assertIn("R/plots.R", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.out, "custom")))

    def test_a_failed_copy_keeps_the_earlier_copy_and_leaves_no_partial(self):
        from unittest import mock
        before = self.fold_in()
        write(os.path.join(self.code, "species.py"), "print('v2')\n")
        items = cw.plan(self.sd, self.out, code_dirs=[self.code])
        real, calls = shutil.copy2, []

        def full_disk(src, dst, *a, **k):
            calls.append(src)
            if len(calls) == 2:
                raise OSError(28, "No space left on device")
            return real(src, dst, *a, **k)

        with mock.patch.object(cw.shutil, "copy2", side_effect=full_disk):
            with self.assertRaises(OSError):
                cw.apply(self.sd, self.out, items)
        self.assertEqual(read(os.path.join(self.dest, "species.py")), "print('v1')\n")
        self.assertTrue(os.path.isfile(os.path.join(self.dest, "R", "plots.R")))
        self.assertEqual(cw.load(self.out), before)                  # the record as it was
        self.assertFalse(os.path.lexists(os.path.join(self.out, "custom", cw.PARTIAL_DIR)))
        self.assertEqual(cw.changed(self.sd, before[0]), [])         # the copy matches its sums

    def test_a_failed_copy_is_said_and_the_earlier_items_still_listed(self):
        import session
        from unittest import mock
        self.fold_in()
        items = cw.plan(self.sd, self.out, code_dirs=[self.code])
        with mock.patch.object(cw.shutil, "copy2", side_effect=OSError(28, "No space left")):
            lines = session.custom_step(session.paths_for(self.sd), {"custom": items})
        levels = [(lvl, part) for lvl, part, _ in lines]
        self.assertIn(("SKIPPED", cw.MANIFEST_PART), levels)
        self.assertIn(("OK", f"{cw.MANIFEST_PART}: output/custom/code/species_analysis"), levels)
        skipped = [n for lvl, part, n in lines if (lvl, part) == ("SKIPPED", cw.MANIFEST_PART)]
        self.assertIn("earlier copies kept", skipped[0])

    def test_a_file_swapped_for_a_link_after_the_plan_is_not_copied(self):
        secret = write(os.path.join(self.root, "elsewhere", "id_key"), "not-a-real-key\n")
        items = cw.plan(self.sd, self.out, code_dirs=[self.code])
        target = os.path.join(self.code, "R", "plots.R")
        os.remove(target)
        os.symlink(secret, target)
        with self.assertRaises(cw.CopyFailed):
            cw.apply(self.sd, self.out, items)
        self.assertFalse(os.path.exists(self.dest))
        self.assertFalse(os.path.lexists(os.path.join(self.out, "custom", cw.PARTIAL_DIR)))

    def test_a_second_fold_in_replaces_the_copy_whole(self):
        """(pin) the move into place: the new copy, nothing of the old left"""
        self.fold_in()
        os.remove(os.path.join(self.code, "R", "plots.R"))
        write(os.path.join(self.code, "species.py"), "print('v2')\n")
        self.fold_in()
        self.assertEqual(read(os.path.join(self.dest, "species.py")), "print('v2')\n")
        self.assertFalse(os.path.exists(os.path.join(self.dest, "R")))
        self.assertFalse(os.path.lexists(os.path.join(self.out, "custom", cw.PARTIAL_DIR)))


bundle = tfg.bundle

if __name__ == "__main__":
    unittest.main()
