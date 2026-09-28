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
import session_docs  # noqa: E402

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


def make_session(tmp, raws=None, conditions=None, fasta=None, extra_col="Batch", prov=None):
    s = os.path.join(tmp, "2026-09-24_PROT_0756")
    os.makedirs(os.path.join(s, "input"))
    os.makedirs(os.path.join(s, "output"))
    if raws is not None:
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
            fh.write("# raw\n" + "".join(f"/data/JUL26/{r}.d\n" for r in raws))
    if conditions is not None:
        with open(os.path.join(s, "input", "conditions.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["File.Name", "Group"] + ([extra_col] if any(len(r) > 2 for r in conditions) else []))
            for row in conditions:
                w.writerow(row)
    if fasta is not None:
        with open(os.path.join(s, "input", "search.fasta.meta.json"), "w") as fh:
            json.dump(fasta, fh)
    if prov is not None:
        os.makedirs(os.path.join(s, "output", "tables"))
        with open(os.path.join(s, "output", "tables", "de_provenance.json"), "w") as fh:
            json.dump(prov, fh)
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

    def test_every_phone_form_is_removed(self):
        for phone in ("(530)555-0100", "+15305550100", "5305550100", "530-5550100", "530.555.0100",
                      "+1 (530) 555-0100", "1 530 555 0100", "+44 20 7946 0958", "+49 30 901820",
                      "tel: 555-0100", "Phone 752-1234"):
            got = sr.clean(f"call {phone} today")
            self.assertEqual(got, "call [phone removed] today", phone)

    def test_quantities_dates_and_names_are_not_phones(self):
        for text in ("dilutions 100-200-3000 ug", "50 mM Ammonium Bicarbonate", "Cat#10004D",
                     "sent 2026-07-15", "Kv2.1", "pH 7.4, 150 mM NaCl", "LRS-100", "30 samples"):
            self.assertEqual(sr.clean(text), text)

    def test_links_are_only_coreomics_and_uniprot_pages(self):
        for url in ("https://x.example/?e=a@b.co", "https://ucdavis.coreomics.com/submissions/x)",
                    "javascript:alert(1)", "http://ucdavis.coreomics.com/submissions/0022066cd85f"):
            self.assertIsNone(sr._url(url), url)
        r = fixture()
        r["submission_data"]["uniprot"] = "https://evil.example/UP000000589"
        self.assertNotIn("href='https://evil", sr.render_html(sr.sanitize(r)))
        r["submission_data"]["uniprot"] = "https://www.uniprot.org/proteomes/UP000000589"
        self.assertIn("href='https://www.uniprot.org/proteomes/UP000000589'", sr.render_html(sr.sanitize(r)))

    def test_attach_sanitizes_whatever_it_is_given(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            sr.attach(s, planted())                      # the RAW record, not sanitized first
            with open(os.path.join(s, "input", "submission.json")) as fh:
                self.assertClean(fh.read(), "input/submission.json")

    def test_odd_shapes_neither_leak_nor_pass_as_blank(self):
        r = fixture()
        r["pi_first_name"] = {"email": "canary.pi@example.org", "first": "Alex"}
        r["submission_data"]["uniprot"] = {"link": "x"}
        rec = sr.sanitize(r)
        self.assertEqual(rec["pi"]["name"], "Example")
        self.assertIsNone(rec["uniprot"], "a non-text answer is missing, not 'left blank'")
        self.assertNotIn("uniprot_blank", notes_by_id(rec))

    def test_a_record_from_the_summary_says_what_the_summary_lacks(self):
        summary = {"schema": "core_submission/1", "internal_id": "PROT_0756", "id": "0022066cd85f",
                   "samples": []}
        page = sr.render_html(sr.sanitize(summary))
        self.assertIn("not in submission_summary.json", page)

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
                fh.write("# Analysis\n\ntext\n\n## Overview\n\nmore\n")
            out = os.path.join(tmp, "a.html")
            p = run("make_analysis_html.py", "--session", s, "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                page = fh.read()
            self.assertNotIn('id="submission"', page)
            self.assertNotIn("<b>Submission</b>", page)
            attach_fixture(s)
            p = run("make_analysis_html.py", "--session", s, "--out", out)
            with open(out) as fh:
                page = fh.read()
            # first in the contents, first section, and the header's Submission fact
            self.assertIn('<ol><li><a href="#submission">Submission</a></li>', page)
            self.assertLess(page.index('id="submission"'), page.index('id="overview"'))
            self.assertIn("<b>Submission</b>PROT_0756", page)
            # styled by report_style.py (the one design system), not by a style of its own
            import report_style
            self.assertIn(".subm dl{", report_style.CSS)
            self.assertEqual(page.count("<style>"), 1)

    def test_a_label_is_not_a_record(self):
        """--submission names a record, never a bare number: a PROT id with no record behind it
        is shown as unreadable, and nothing is put in the header."""
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            with open(os.path.join(s, "output", "AI_Analysis_Report.md"), "w") as fh:
                fh.write("# Analysis\n\n## Overview\n\nmore\n")
            out = os.path.join(tmp, "a.html")
            p = run("make_analysis_html.py", "--session", s, "--submission", "PROT_0001",
                    "--out", out)
            self.assertEqual(p.returncode, 0, p.stderr)
            with open(out) as fh:
                page = fh.read()
            self.assertIn("could not be read", page)
            self.assertIn('class="callout warning"', page)
            self.assertNotIn("<b>Submission</b>", page)

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
            self.assertIsNone(session_docs.submission_line(p))
            attach_fixture(s)
            line = ("CoreOmics submission: [PROT_0756](https://ucdavis.coreomics.com/submissions/"
                    "0022066cd85f) (source: CoreOmics)")
            self.assertEqual(session_docs.submission_line(p), line)
            # the README (and AGENTS.md's study summary) that finalize writes carry it
            written = session_docs.write_docs(s)
            for key in ("readme", "agents"):
                with open(written[key], encoding="utf-8") as fh:
                    self.assertIn(f"- {line}", fh.read(), key)
            os.remove(p["submission_record"])           # announced but gone: said, not dropped
            self.assertIn("could not be read", session_docs.submission_line(p))


# ----------------------------------------------------------------- who prepared --
class TestPreparedBy(unittest.TestCase):
    def test_the_form_answers(self):
        self.assertEqual(sr.prepared_by(sr.sanitize(fixture()))[0], "lab")
        core = sr.sanitize({"sample_prep": "I want the proteomics core to prepare my samples"})
        self.assertEqual(sr.prepared_by(core)[0], "core")
        self.assertEqual(sr.prepared_by(sr.sanitize({"prot_or_pep": "peptides"}))[0], "lab")
        self.assertIsNone(sr.prepared_by(sr.sanitize({"prot_or_pep": "Intact Proteins"}))[0],
                          "intact proteins do not say who extracted them")
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

    def test_peptides_only_when_the_form_says_peptides(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            sr.attach(s, {"internal_id": "PROT_0756", "sample_prep": "lab"})
            sec = self.methods(s)
            self.assertIn("prepared by the submitting laboratory", sec)
            self.assertNotIn("peptides", sec)
            self.assertIn("details given by the user", sec)

    def test_multi_line_buffer_text_stays_out_of_the_pride_protocol(self):
        import make_deposit
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            r = fixture()
            r["submission_data"]["buffer"] = "50 mM ABC\nwashed in RIPA *twice*"
            attach_fixture(s, r)
            sec = self.methods(s)
            note = [ln for ln in sec.splitlines() if "own protocol" in ln]
            self.assertEqual(len(note), 1)
            self.assertIn("\u201c50 mM ABC washed in RIPA twice\u201d", note[0])
            sample, _ = make_deposit.build_protocols(os.path.join(s, "output", "methods.md"))
            self.assertNotIn("RIPA", sample)
            self.assertNotIn("own protocol", sample)

    def test_a_methods_file_from_before_the_attach_is_not_enough(self):
        import make_deposit
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp)
            f = {"p": session_mod.paths_for(s), "srec": {}, "params": None, "search_prov_path": None,
                 "de_prov": None}
            self.assertNotIn("Sample preparation", make_deposit.required_sections(f))
            attach_fixture(s)
            self.assertIn("Sample preparation", make_deposit.required_sections(f))

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
            self.assertIn("already lists these notes, so do not copy them", brief)
            self.assertIn("The submitting lab prepared the samples and sent peptides", brief)
            self.assertIn("Pairing in the sample sheet", brief)


# ------------------------------------------------------------------------- notes --
class TestNotes(unittest.TestCase):
    def session(self, tmp, raws=tuple(RUNS.values()), groups=by_age_and_ip, batch=None,
                fasta=({"organism": "Mus musculus", "taxid": 10090},), extra_col="Batch", prov=None):
        conds = [[RUNS[u], groups(u)] + ([batch(u)] if batch else []) for u in RUNS]
        s = make_session(tmp, raws=list(raws), conditions=conds, fasta=fasta[0] if fasta else None,
                         extra_col=extra_col, prov=prov)
        attach_fixture(s)
        return s

    def test_prot0756_pairing_is_surfaced(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp)
            n = notes_by_id(sr.load(s), s)
            self.assertIn("pairing", n)
            for part in ("carry the mouse each sample came from (Mouse 1–6)",
                         "Each mouse gave samples under all 5 of JPH3, JPH4, Kv2.1, RyR, IgG",
                         "within-mouse (paired)", "Old vs Young differs between mice "
                         "(Old: Mouse 1–3; Young: Mouse 4–6)", "10 groups of 3",
                         "No column the DE reads from input/conditions.csv (Group, Batch, Covariate1, "
                         "Covariate2) identifies the mouse, so as set up it will treat the samples as "
                         "independent"):
                self.assertIn(part, n["pairing"])
            self.assertEqual(set(n) & {"conditions_differ", "sheet_ids_without_raw",
                                       "raw_without_sheet_id", "organism_mismatch"}, set())

    def test_a_design_column_that_carries_the_mouse_is_recognised(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, batch=lambda u: "m" + sheet()[u].split("Mouse")[1].strip())
            self.assertIn("includes the mouse as “Batch”, a fixed effect",
                          notes_by_id(sr.load(s), s)["pairing"])

    def test_restarted_numbering_is_not_read_as_pairing(self):
        """Old Mouse 1-3 and Young Mouse 1-3: are they six mice or three? The sheet cannot say."""
        for fmt in ("{age} - Mouse {n}", "Mouse {n} - {age}"):
            r = fixture()
            r["submission_data"]["samples"] = [
                {"unique_id": f"S{i}", "sample_name": f"S{i}", "condition_name": fmt.format(age=age, n=n)}
                for i, (age, n) in enumerate([(a, n) for a in ("Old", "Young") for n in (1, 2, 3)])]
            n = notes_by_id(sr.sanitize(r))
            self.assertIn("does not say whether Mouse 1 under Old and under Young is the same mouse",
                          n["pairing"], fmt)
            self.assertNotIn("(paired)", n["pairing"])
            self.assertNotIn("not independent", n["pairing"])

    def test_one_sample_per_mouse_is_no_pairing(self):
        r = fixture()
        r["submission_data"]["samples"] = [
            {"unique_id": f"S{n}", "sample_name": "", "condition_name": f"{'Old' if n < 4 else 'Young'} - Mouse {n}"}
            for n in range(1, 7)]
        self.assertNotIn("pairing", notes_by_id(sr.sanitize(r)))

    def test_replicate_indexes_are_replicate_groups(self):
        for fmt in ("{g} - Rep {n}", "{g} {n}", "{g}_{n}"):
            r = fixture()
            r["submission_data"]["samples"] = [
                {"unique_id": f"S{g}{n}", "sample_name": "", "condition_name": fmt.format(g=g, n=n)}
                for g in ("Control", "Treated") for n in (1, 2, 3)]
            rec = sr.sanitize(r)
            n = notes_by_id(rec)
            self.assertNotIn("sheet_conditions_unique", n, fmt)
            self.assertNotIn("pairing", n, fmt)
            self.assertEqual(set(sr.sheet_groups(rec).values()), {"Control", "Treated"}, fmt)

    def test_a_sample_outside_the_pattern_keeps_the_pairing_note(self):
        r = fixture()
        r["submission_data"]["samples"].append({"unique_id": "LRS200", "sample_name": "LRS200",
                                                "condition_name": "Pool"})
        n = notes_by_id(sr.sanitize(r))
        self.assertIn("1 sample(s) do not follow this naming and are not part of this reading (LRS200)",
                      n["pairing"])
        self.assertNotIn("sheet_conditions_unique", n)

    def test_a_reinjected_sample_is_not_left_out(self):
        with tempfile.TemporaryDirectory() as tmp:
            conds = [[RUNS[u], by_age_and_ip(u)] for u in RUNS]
            conds.append(["08202026__60SPD_DIA-LRS-96_rerun_S3-A1_1_23999", by_age_and_ip("LRS96")])
            s = make_session(tmp, conditions=conds)
            attach_fixture(s)
            self.assertNotIn("conditions_differ", notes_by_id(sr.load(s), s))

    def test_a_run_carrying_two_sheet_ids_says_nothing_about_groups(self):
        r = fixture()
        r["submission_data"]["samples"] = [
            {"unique_id": u, "sample_name": "", "condition_name": c}
            for u, c in (("AB12", "Old"), ("CD34", "Young"), ("EF56", "Old"), ("GH78", "Young"))]
        with tempfile.TemporaryDirectory() as tmp:
            s = make_session(tmp, conditions=[["x_AB12_CD34", "Old"], ["x_EF56", "Old"], ["x_GH78", "Young"]])
            attach_fixture(s, r)
            self.assertNotIn("conditions_differ", notes_by_id(sr.load(s), s))

    def test_only_columns_the_de_models_count_as_carrying_the_mouse(self):
        mouse = lambda u: sheet()[u].split("Mouse")[1].strip()        # noqa: E731
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, batch=mouse, extra_col="Mouse")
            self.assertIn("input/conditions.csv carries the mouse as \u201cMouse\u201d: the DE models "
                          "it only when run with --block Mouse", notes_by_id(sr.load(s), s)["pairing"],
                          "run_de.R reads a Mouse column only with --block")
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, prov={"design": "~ 0 + groups", "n_samples": 30,
                                        "block": {"applied": False,
                                                  "note": "no --block: samples modelled as independent"}})
            n = notes_by_id(sr.load(s), s)["pairing"]
            self.assertIn("The design analysed (~ 0 + groups) has no term for the mouse, so it "
                          "treated the 30 samples as independent.", n)
            self.assertNotIn("blocked", n, "{'applied': false} is an unblocked run, not a block")

    # run_de.R's block record (feat/skill-de-block, blocking.R block_record): `block` is an
    # object, `block_column` exists only when the design is blocked.
    BLOCKED_RANDOM = {"design": "~ 0 + groups", "n_samples": 30, "block_column": "Mouse",
                      "block": {"applied": True, "column": "Mouse", "effect": "random",
                                "consensus_correlation": 0.412,
                                "contrast_structure": {"Old_JPH3-Old_IgG": "within",
                                                       "Old_JPH3-Young_JPH3": "between"}}}

    def test_a_random_block_on_the_mouse(self):
        mouse = lambda u: sheet()[u].split("Mouse")[1].strip()        # noqa: E731
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, batch=mouse, extra_col="Mouse", prov=self.BLOCKED_RANDOM)
            n = notes_by_id(sr.load(s), s)["pairing"]
            self.assertIn("blocked on the mouse (\u201cMouse\u201d, a random effect; consensus "
                          "within-mouse correlation 0.41).", n)
            self.assertIn("1 comparison(s) (Old_JPH3-Young_JPH3) set different mice against each "
                          "other: their evidence is the number of mice", n)
            self.assertNotIn("samples as independent", n)

    def test_a_fixed_block_on_the_mouse(self):
        mouse = lambda u: sheet()[u].split("Mouse")[1].strip()        # noqa: E731
        prov = {"design": "~ 0 + groups + Mouse", "n_samples": 30, "block_column": "Mouse",
                "block": {"applied": True, "column": "Mouse", "effect": "fixed",
                          "consensus_correlation": None}}
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, batch=mouse, extra_col="Mouse", prov=prov)
            self.assertIn("modelled the mouse (\u201cMouse\u201d) as a fixed effect (--block), so "
                          "every comparison is made within each mouse",
                          notes_by_id(sr.load(s), s)["pairing"])

    def test_a_record_from_before_block_says_not_recorded(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, prov={"design": "~ 0 + groups", "n_samples": 30})
            n = notes_by_id(sr.load(s), s)["pairing"]
            self.assertIn("Whether the analysis modelled the mouse is not recorded", n)
            self.assertNotIn("treated the 30 samples as independent", n)

    def test_block_given_but_every_contrast_between_mice(self):
        prov = {"design": "~ 0 + groups", "n_samples": 30,
                "block": {"applied": False, "column": "Mouse",
                          "note": "--block Mouse given, but every contrast compares different "
                                  "Mouse levels"}}
        with tempfile.TemporaryDirectory() as tmp:
            s = self.session(tmp, prov=prov)
            n = notes_by_id(sr.load(s), s)["pairing"]
            self.assertIn("treated the 30 samples as independent (--block Mouse given, but every "
                          "contrast compares different Mouse levels).", n)

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
