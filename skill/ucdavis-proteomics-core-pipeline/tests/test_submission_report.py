"""
The CoreOmics submission in every report (submission_report.py).

The report goes to the collaborator and gets forwarded, so the record is ALLOWLISTED: no email,
phone, billing or internal field may reach it, however CoreOmics or the submitter spells it.
The fixture is PROT_0756's real record with every contact field stripped and fake people's
names (tests/fixtures/coreomics_prot0756); the canaries below are planted at test time.

Stdlib only, no network.
"""
import copy
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import submission_report as sr  # noqa: E402
import record_run  # noqa: E402
import session as session_mod  # noqa: E402

FIXTURE = os.path.join(HERE, "fixtures", "coreomics_prot0756")
EMAIL_RE = re.compile(r"[\w.+-]+@[\w-]+\.[\w.-]+")
RUNS = {f"LRS{n}": f"08132026__60SPD_DIA-LRS-{n}_S3-A1_1_{23600 + n}" for n in range(96, 126)}

# Planted in every place a contact, billing or internal detail can sit in a CoreOmics record.
CANARIES = ("canary.submitter@example.org", "canary.pi@example.org", "canary.contact@example.org",
            "530-555-0142", "(530) 555-0143", "PO-CANARY-7731", "CANARY-PPMS-REF", "CANARY-TOKEN-abc",
            "CANARY-INTERNAL-NOTE", "CANARY-FUTURE-FIELD", "CANARY-STAFF-ROUTING")


def fixture():
    with open(os.path.join(FIXTURE, "submission.json")) as fh:
        return json.load(fh)


def planted():
    r = fixture()
    r.update(email="canary.submitter@example.org", phone="530-555-0142",
             pi_email="canary.pi@example.org", pi_phone="(530) 555-0143",
             payment={"display": {"Payment Type": "PO", "Payment Info": "PO-CANARY-7731",
                                  "PPMS Order Ref #": "CANARY-PPMS-REF"}},
             import_request="CANARY-TOKEN-abc", notes="CANARY-INTERNAL-NOTE",
             future_field="CANARY-FUTURE-FIELD",
             contacts=[{"first_name": "Robin", "last_name": "Contact",
                        "email": "canary.contact@example.org"}])
    r["pi"] = dict(r["pi"], email="canary.pi@example.org", phone="(530) 555-0143")
    sd = r["submission_data"]
    sd["email_ucd_scientist"] = ["CANARY-STAFF-ROUTING"]
    sd["description"] += " Questions: canary.submitter@example.org or 530-555-0142."
    sd["samples"][0] = dict(sd["samples"][0], internal_notes="CANARY-INTERNAL-NOTE")
    return r


def make_session(tmp, raws=None, conditions=None, fasta=None):
    s = os.path.join(tmp, "2026-09-24_PROT_0756")
    os.makedirs(os.path.join(s, "input"))
    os.makedirs(os.path.join(s, "output"))
    if raws is not None:
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
            fh.write("# raw\n" + "".join(f"/data/JUL26/{r}.d\n" for r in raws))
    if conditions is not None:
        with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["File.Name", "Group"] + (["Batch"] if any(len(r) > 2 for r in conditions) else []))
            for row in conditions:
                w.writerow(row)
    if fasta is not None:
        with open(os.path.join(s, "input", "search.fasta.meta.json"), "w") as fh:
            json.dump(fasta, fh)
    return s


def sheet():
    return {s["unique_id"]: s["condition_name"] for s in fixture()["submission_data"]["samples"]}


def by_age_and_ip(uid):
    f = [x.strip() for x in re.split(r"\s+-\s*|\s*-\s+", sheet()[uid])]
    return f"{f[0]}_{f[1]}"


def attach_fixture(s, obj=None):
    return sr.attach(s, sr.sanitize(obj or fixture()))


def run(script, *args):
    return subprocess.run([sys.executable, os.path.join(SCRIPTS, script)] + list(args),
                          capture_output=True, text=True, timeout=120)


def notes_by_id(rec, session=None):
    return {n["id"]: n["text"] for n in sr.quality_notes(rec, session)}


# ------------------------------------------------------------------------ allowlist --
class TestAllowlist(unittest.TestCase):
    def assertClean(self, text, where):
        for c in CANARIES:
            self.assertNotIn(c, text, f"{c} leaked into {where}")
        self.assertIsNone(EMAIL_RE.search(text), f"an email address leaked into {where}")

    def test_no_contact_billing_or_internal_field_reaches_any_output(self):
        rec = sr.sanitize(planted())
        notes = sr.quality_notes(rec)
        self.assertClean(json.dumps(rec), "the record")
        self.assertClean(sr.render_html(rec, notes), "the HTML section")
        self.assertClean(sr.render_markdown(rec, notes), "the markdown section")
        self.assertIn("[email removed]", rec["description"])
        self.assertIn("[phone removed]", rec["description"])

    def test_nor_the_session_files_or_the_finished_report(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            attach_fixture(s, planted())
            for rel in ("session.json", "input/submission.json", "input/samples.tsv"):
                with open(os.path.join(s, rel)) as fh:
                    self.assertClean(fh.read(), rel)
            with open(os.path.join(s, "output", "AI_Analysis_Report.md"), "w") as fh:
                fh.write("# Analysis\n\ntext\n")
            out = os.path.join(tmp, "report.html")
            p = run("make_analysis_html.py", "--session", s, "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                page = fh.read()
            self.assertClean(page, "Analysis_Report.html")
            self.assertIn("PROT_0756", page)

    def test_the_summary_and_facts_given_by_the_user_are_allowlisted_too(self):
        summary = {"schema": "core_submission/1", "internal_id": "PROT_0756", "id": "0022066cd85f",
                   "submitter": {"name": "Sam Sampleton", "email": "canary.submitter@example.org"},
                   "pi": {"name": "Alex Example", "email": "canary.pi@example.org"},
                   "contacts": [{"email": "canary.contact@example.org"}], "samples": []}
        self.assertClean(json.dumps(sr.sanitize(summary)), "a record from the summary")
        given = sr.sanitize({"internal_id": "756", "pi": {"name": "A", "email": "canary.pi@example.org"},
                             "phone": "530-555-0142"})
        self.assertClean(json.dumps(given), "facts given by the user")
        self.assertEqual((given["internal_id"], given["source"]), ("PROT_0756", sr.SOURCE_USER))

    def test_sanitize_is_idempotent(self):
        rec = sr.sanitize(fixture())
        self.assertEqual(sr.sanitize(copy.deepcopy(rec)), rec)

    def test_submitter_text_is_escaped(self):
        r = fixture()
        r["submission_data"]["samples"][0]["condition_name"] = "<script>alert(1)</script>"
        r["submission_data"]["description"] = "<img src=x onerror=alert(1)>"
        page = sr.render_html(sr.sanitize(r))
        self.assertNotIn("<script>alert", page)
        self.assertNotIn("<img", page)


# ------------------------------------------------------------------- the section --
class TestSection(unittest.TestCase):
    def test_renders_everything_the_lab_told_the_core(self):
        rec = sr.sanitize(fixture())
        page = sr.render_html(rec, sr.quality_notes(rec))
        for expect in ("href='https://ucdavis.coreomics.com/submissions/0022066cd85f'>PROT_0756</a>",
                       "Alex Example, SOM: Example Department, UC Davis", "Sam Sampleton",
                       "2026-07-15", "Mouse", "left blank on the form",
                       "DIA (Quantitative); Affinity Purification (magnetic Beads)",
                       "peptides", "the lab prepared peptides",
                       "50 mM Ammonium Bicarbonate", "Protein G coated magnetic beads",
                       "lets discuss/no idea", "I want the proteomics core to do the data analysis",
                       "Sample sheet (30 samples)", "<code>LRS96</code>Old - JPH3 - Mouse 1",
                       "Data quality notes"):
            self.assertIn(expect, page)

    def test_the_description_is_as_written(self):
        """The report once said "chemically cross-linked"; the record says "cross-linked"."""
        rec = sr.sanitize(fixture())
        for text in (sr.render_html(rec), sr.render_markdown(rec)):
            self.assertIn("Immunopurifications of cross-linked brain tissue", text)
            self.assertNotIn("chemically", text)

    def test_markdown_section(self):
        md = sr.render_markdown(sr.sanitize(fixture()))
        self.assertTrue(md.startswith("## Submission\n"))
        self.assertIn("| Submission | [PROT_0756](https://ucdavis.coreomics.com/submissions/0022066cd85f) |", md)
        self.assertIn("| LRS125 | Young - IgG - Mouse 6 |", md)

    def test_facts_given_by_the_user_are_labelled_so(self):
        rec = sr.sanitize({"internal_id": "PROT_0756", "organism": "mouse", "sample_prep": "lab"})
        self.assertIn("not read from CoreOmics", sr.render_html(rec))
        self.assertIn("source_user", notes_by_id(rec))
        self.assertIn("not in the record", sr.render_html(rec), "missing facts say so")

    def test_report_leads_with_it_and_omits_it_when_there_is_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            with open(os.path.join(s, "output", "AI_Analysis_Report.md"), "w") as fh:
                fh.write("# Analysis\n\ntext\n")
            out = os.path.join(tmp, "a.html")
            p = run("make_analysis_html.py", "--session", s, "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                self.assertNotIn('id="submission"', fh.read())
            attach_fixture(s)
            p = run("make_analysis_html.py", "--session", s, "--out", out)
            with open(out) as fh:
                page = fh.read()
            self.assertIn('<ul><li><a href="#submission">Submission</a></li>', page)
            self.assertLess(page.index('id="submission"'), page.index('id="analysis"'))

    def test_an_explicit_record_needs_no_session(self):
        with tempfile.TemporaryDirectory() as tmp:
            rep = os.path.join(tmp, "r.md")
            with open(rep, "w") as fh:
                fh.write("# A\n\nx\n")
            out = os.path.join(tmp, "a.html")
            p = run("make_analysis_html.py", "--report", rep, "--submission", FIXTURE, "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                self.assertIn("PROT_0756", fh.read())

    def test_an_unreadable_record_is_shown_not_dropped(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            attach_fixture(s)
            os.remove(os.path.join(s, "input", "submission.json"))
            self.assertIn("could not be read", sr.html_section(None, s))


# --------------------------------------------------------------------- storage --
class TestAttach(unittest.TestCase):
    def test_one_record_in_the_session(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            with open(os.path.join(s, "session.json"), "w") as fh:
                json.dump({"qc": False}, fh)
            p = run("submission_report.py", "attach", "--session", s, "--record", FIXTURE)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(os.path.join(s, "session.json")) as fh:
                sj = json.load(fh)
            self.assertFalse(sj["qc"], "other session metadata is kept")
            self.assertEqual(sj["coreomics"]["internal_id"], "PROT_0756")
            self.assertEqual(sj["coreomics"]["record"], "input/submission.json")
            self.assertEqual(sr.load(s)["samples"][0]["unique_id"], "LRS96")
            with open(os.path.join(s, "input", "samples.tsv")) as fh:
                self.assertEqual(len(fh.read().splitlines()), 31)

    def test_another_submission_is_not_silently_replaced(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            attach_fixture(s)
            attach_fixture(s)                                   # the same one again: fine
            with self.assertRaises(sr.RecordError):
                sr.attach(s, sr.sanitize({"internal_id": "PROT_0757"}))
            sr.attach(s, sr.sanitize({"internal_id": "PROT_0757"}), replace=True)
            self.assertEqual(sr.load(s)["internal_id"], "PROT_0757")

    def test_given_facts_warn_about_what_was_not_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            p = run("submission_report.py", "attach", "--session", s, "--given",
                    json.dumps({"internal_id": "PROT_0756", "sample_prep": "lab", "email": "x@example.org"}))
            self.assertEqual(p.returncode, 0, p.stderr)
            out = json.loads(p.stdout)
            self.assertEqual((out["source"], out["prepared_by"]), (sr.SOURCE_USER, "lab"))
            self.assertIn("email", out["warnings"][0])

    def test_record_run_finds_both_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            attach_fixture(s)
            got = record_run.find_prot(None, s, None)
            self.assertEqual((got["prot"], got["id"]), ("PROT_0756", "0022066cd85f"))
            self.assertTrue(got["source"].endswith("session.json"))

    def test_the_session_readme_names_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            p = session_mod.paths_for(s)
            self.assertIsNone(session_mod._submission_line(p))
            attach_fixture(s)
            self.assertEqual(session_mod._submission_line(p),
                             "CoreOmics submission: [PROT_0756](https://ucdavis.coreomics.com/submissions/"
                             "0022066cd85f) (source: CoreOmics)")


# ----------------------------------------------------------------- who prepared --
class TestPreparedBy(unittest.TestCase):
    def test_the_form_answers(self):
        self.assertEqual(sr.prepared_by(sr.sanitize(fixture()))[0], "lab")
        core = sr.sanitize({"sample_prep": "I want the proteomics core to prepare my samples"})
        self.assertEqual(sr.prepared_by(core)[0], "core")
        self.assertEqual(sr.prepared_by(sr.sanitize({"prot_or_pep": "peptides"}))[0], "lab")
        self.assertEqual(sr.prepared_by(sr.sanitize({"prot_or_pep": "Intact Proteins"}))[0], "core")
        self.assertEqual(sr.prepared_by(sr.sanitize({}))[0], None)

    def test_a_contradictory_form_is_not_resolved_by_guessing(self):
        rec = sr.sanitize({"sample_prep": "I have prepped my samples, they are peptides and ready "
                                          "to analyzed by LC-MSMS", "prot_or_pep": "Intact Proteins"})
        self.assertIsNone(sr.prepared_by(rec)[0])
        self.assertIn("prep_unclear", notes_by_id(rec))


class TestMethods(unittest.TestCase):
    def methods(self, s):
        out = os.path.join(s, "output", "methods.md")
        p = run("make_methods.py", "--raw", "/nowhere/a.d", "--instrument", "timsTOF HT",
                "--acquisition", "DIA", "--out", out, "--submission", s)
        self.assertEqual(p.returncode, 0, p.stderr)
        with open(out) as fh:
            text = fh.read()
        return text[text.index("## Sample preparation"):text.index("## Liquid chromatography")]

    def test_peptides_from_the_lab_are_the_labs_preparation(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            attach_fixture(s)
            sec = self.methods(s)
            self.assertIn("prepared by the submitting laboratory", sec)
            self.assertIn("as peptides ready for LC-MS/MS (CoreOmics submission PROT_0756)", sec)
            self.assertNotIn("confirm]", sec, "no Core-side placeholder the Core could never fill")
            self.assertNotIn("chemically", sec)
            self.assertIn("“The samples are in 50 mM Ammonium Bicarbonate", sec)

    def test_core_prepared_samples_keep_a_placeholder_for_the_core(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            r = fixture()
            r["submission_data"].update(sample_prep="I want the proteomics core to prepare my samples",
                                        prot_or_pep="Intact Proteins")
            attach_fixture(s, r)
            sec = self.methods(s)
            self.assertIn("prepared by the UC Davis Proteomics Core", sec)
            self.assertIn("[not recorded — confirm]", sec)

    def test_the_deposit_protocol_uses_it_instead_of_the_to_fill(self):
        import make_deposit
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            attach_fixture(s)
            self.methods(s)
            sample, _data = make_deposit.build_protocols(os.path.join(s, "output", "methods.md"))
            self.assertTrue(sample.startswith("Samples were prepared by the submitting laboratory"))
            self.assertNotIn("TO-FILL", sample)
            self.assertNotIn("own protocol", sample, "the author's note is not protocol text")

    def test_without_a_submission_nothing_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = os.path.join(tmp, "m.md")
            p = run("make_methods.py", "--raw", "/nowhere/a.d", "--instrument", "timsTOF HT", "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                self.assertNotIn("## Sample preparation", fh.read())

    def test_the_analysis_brief_quotes_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            attach_fixture(s)
            de = os.path.join(tmp, "tables")
            os.makedirs(de)
            out = os.path.join(tmp, "PROMPT.md")
            p = run("analysis_prompt.py", "--de-dir", de, "--out", out, "--submission", s)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                brief = fh.read()
            self.assertIn("## The submission (PROT_0756) — quote it, add nothing", brief)
            self.assertIn("The submitting lab prepared the samples and sent peptides", brief)
            self.assertIn("Pairing in the sample sheet", brief)


# ------------------------------------------------------------------------- notes --
class TestNotes(unittest.TestCase):
    def session(self, tmp, raws=tuple(RUNS.values()), groups=by_age_and_ip, batch=None,
                fasta=({"organism": "Mus musculus", "taxid": 10090},)):
        conds = [[RUNS[u], groups(u)] + ([batch(u)] if batch else []) for u in RUNS]
        s = make_session(tmp, raws=list(raws), conditions=conds, fasta=fasta[0] if fasta else None)
        attach_fixture(s)
        return s

    def test_prot0756_pairing_is_surfaced(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp)
            n = notes_by_id(sr.load(s), s)
            self.assertIn("pairing", n)
            for part in ("carries the mouse the sample came from (Mouse 1–6)",
                         "Each mouse gave samples under all 5 of JPH3, JPH4, Kv2.1, RyR, IgG",
                         "within-mouse (paired)", "Old vs Young differs between mice "
                         "(Old: Mouse 1–3; Young: Mouse 4–6)", "10 groups of 3",
                         "has no column for the mouse, so the model treated all 30 samples as independent"):
                self.assertIn(part, n["pairing"])
            self.assertEqual(set(n) & {"conditions_differ", "sheet_ids_without_raw",
                                       "raw_without_sheet_id", "organism_mismatch"}, set())

    def test_a_design_column_that_carries_the_mouse_is_recognised(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, batch=lambda u: "m" + sheet()[u].split("Mouse")[1].strip())
            self.assertIn("carries the mouse as “Batch”", notes_by_id(sr.load(s), s)["pairing"])

    def test_blank_uniprot_is_a_note(self):
        self.assertIn("uniprot_blank", notes_by_id(sr.sanitize(fixture())))

    def test_organism_against_the_search_database(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, fasta=({"organism": "Homo sapiens", "taxid": 9606},))
            self.assertIn("but the search database is Homo sapiens (taxid 9606)",
                          notes_by_id(sr.load(s), s)["organism_mismatch"])
        r = fixture()
        r["submission_data"]["organism"] = ""
        self.assertIn("organism_missing", notes_by_id(sr.sanitize(r)))

    def test_sample_ids_against_the_raw_files(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = [v for k, v in RUNS.items() if k != "LRS125"] + ["08142026__60SPD_DIA-QC-HeLa_S3-A1_1_1"]
            s = self.session(tmp, raws=raws)
            n = notes_by_id(sr.load(s), s)
            self.assertIn("1 of 30 sample IDs on the sheet match no raw file analysed here (LRS125)",
                          n["sheet_ids_without_raw"])
            self.assertIn("08142026__60SPD_DIA-QC-HeLa_S3-A1_1_1.d", n["raw_without_sheet_id"])

    def test_conditions_that_differ_from_the_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, groups=lambda u: by_age_and_ip(u).split("_")[0])
            self.assertIn("group “Old” pools the sheet's Old - IgG, Old - JPH3",
                          notes_by_id(sr.load(s), s)["conditions_differ"])
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, groups=lambda u: "Young_JPH3" if u == "LRS96" else by_age_and_ip(u))
            self.assertIn("“Young_JPH3” pools the sheet's Old - JPH3, Young - JPH3",
                          notes_by_id(sr.load(s), s)["conditions_differ"])

    def test_a_sheet_with_no_groups_says_so(self):
        r = fixture()
        for i, smp in enumerate(r["submission_data"]["samples"]):
            smp["condition_name"] = f"sample {i}"
        self.assertIn("sheet_conditions_unique", notes_by_id(sr.sanitize(r)))

    def test_beads_left_blank_for_an_affinity_experiment(self):
        r = fixture()
        r["submission_data"]["magbead_yes_no"] = ""
        self.assertIn("beads_blank", notes_by_id(sr.sanitize(r)))


if __name__ == "__main__":
    unittest.main()
