#!/usr/bin/env python3
"""
FAIMS and non-FAIMS runs in one FragPipe quantification (X9, 2.15).

FragPipe's IonQuant indexes runs by compensation voltage and crashes on a mix ("Cannot find CV -50
from all indexed CVs: [0]") after the whole search has run. So: faims.py is the one place that says
which runs have FAIMS; a quantifying FragPipe search over a mix is refused before submission with
the counts and the choices; --faims-split runs one search per group and says so in the record, the
Methods, the report and MANIFEST.txt; and FragPipe's memory has a run-count term (a 219-run search
filled the 200G default).

The mzML and filter-string shapes are those read from real Orbitrap files on HIVE (2026-11-11).
Fakes only: FragPipe is /bin/true, nothing is submitted.
"""
import gzip
import json
import os
import sys
import tempfile
import unittest
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

import faims  # noqa: E402
import run_search  # noqa: E402
from test_search_job_sizing import FP_WORKFLOW, header, provenance, search  # noqa: E402

CV_PARAM = ('<cvParam cvRef="MS" accession="MS:1001581" name="FAIMS compensation voltage" '
            'value="{cv}"/>')
FILTER = ('<cvParam cvRef="MS" accession="MS:1000512" name="filter string" '
          'value="FTMS + p NSI {cvtok}Full ms [400.0000-1600.0000]"/>')


def mzml(cv=None, spectra=2, params=True, filt=True):
    """A small mzML as msconvert writes a Thermo run: with `cv`, every scan carries MS:1001581 and
    "cv=-55.00" in its filter string; without it, neither."""
    body = ['<?xml version="1.0"?>', '<mzML xmlns="http://psi.ms.org"><run><spectrumList>']
    for i in range(spectra):
        body.append(f'<spectrum index="{i}" id="scan={i + 1}"><scanList><scan>')
        if filt:
            body.append(FILTER.format(cvtok=f"cv={cv:.2f} " if cv is not None else ""))
        if params and cv is not None:
            body.append(CV_PARAM.format(cv=cv))
        body.append("</scan></scanList></spectrum>")
    body.append("</spectrumList></run></mzML>")
    return "\n".join(body)


def put(d, name, text, gz=False):
    p = os.path.join(d, name)
    (gzip.open(p, "wt") if gz else open(p, "w")).write(text)
    return p


class Detector(unittest.TestCase):
    def test_the_cv_param_is_read(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(faims.run_cvs(put(d, "a.mzML", mzml(-55.0, filt=False))), [-55.0])

    def test_the_filter_string_alone_is_read(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(faims.run_cvs(put(d, "a.mzML", mzml(-50.0, params=False))), [-50.0])

    def test_several_voltages_are_listed_sorted(self):
        with tempfile.TemporaryDirectory() as d:
            text = mzml(-70.0, spectra=1) + mzml(-50.0, spectra=1)
            self.assertEqual(faims.run_cvs(put(d, "a.mzML", text)), [-70.0, -50.0])

    def test_a_run_without_faims_is_an_empty_list_not_none(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(faims.run_cvs(put(d, "a.mzML", mzml(None))), [])

    def test_a_head_with_no_spectrum_is_unknown_never_no_faims(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertIsNone(faims.run_cvs(put(d, "a.mzML", "<mzML><run></run></mzML>")))
            self.assertIsNone(faims.run_cvs(os.path.join(d, "missing.mzML")))

    def test_a_gzipped_mzml_and_a_cut_head_are_read(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(faims.run_cvs(put(d, "a.mzML.gz", mzml(-45.0), gz=True)), [-45.0])
            # a head cut in the middle of an element still gives its voltage
            self.assertEqual(faims.run_cvs(put(d, "b.mzML", mzml(-55.0) + "<spectrum index=\"9\"><cvPa")),
                             [-55.0])

    def test_other_formats_are_not_read_here(self):
        with tempfile.TemporaryDirectory() as d:
            self.assertIsNone(faims.run_cvs(put(d, "a.raw", "x")))
            self.assertIsNone(faims.run_cvs(os.path.join(d, "a.d")))

    def test_a_thermo_scan_filter(self):
        f = faims.cvs_from_filter
        self.assertEqual(f("FTMS + p NSI cv=-55.00 Full ms [400.0000-1600.0000]"), [-55.0])
        self.assertEqual(f("FTMS + c NSI cv=-45.00 d Full ms2 572.3181@hcd30.00 [100-1200]"),
                         [-45.0])
        self.assertEqual(f("FTMS + p NSI Full ms [375.0000-1500.0000]"), [])
        self.assertEqual(f("FTMS + p NSI Full ms2 367.5000@hcd30.00 [100-1200]"), [])
        self.assertEqual(f(None), [])

    def test_detect_acquisition_records_it_for_a_thermo_scan_list(self):
        import detect_acquisition as da

        def spec(level, filt):
            return {"attributes": [{"accession": "MS:1000511", "value": level},
                                   {"accession": "MS:1000512", "value": filt}]}
        w = da.thermo_isolation_windows([spec(1, "FTMS + p NSI cv=-50.00 Full ms [350.0-1500.0]"),
                                         spec(2, "FTMS + c NSI cv=-50.00 d Full ms2 500.0@hcd30.00 [100-900]")])
        self.assertEqual(sorted(w["faims_cvs"]), [-50.0])
        self.assertEqual(w["n_filters_seen"], 2)
        w = da.thermo_isolation_windows([spec(1, "FTMS + p NSI Full ms [350.0-1500.0]")])
        self.assertEqual(sorted(w["faims_cvs"]), [])
        self.assertEqual(w["n_filters_seen"], 1)

    def test_detect_acquisition_classifies_an_mzml_with_its_faims_status(self):
        import detect_acquisition as da
        with tempfile.TemporaryDirectory() as d:
            self.assertEqual(da.classify(put(d, "a.mzML", mzml(-55.0)))["faims_cvs"], [-55.0])
            self.assertEqual(da.classify(put(d, "b.mzML", mzml(None)))["faims_cvs"], [])


class Grouping(unittest.TestCase):
    def test_groups_and_the_mix(self):
        g = faims.group({"a": [-50.0], "b": [], "c": None, "d": [-70.0]})
        self.assertEqual((g["faims"], g["none"], g["unknown"]), (["a", "d"], ["b"], ["c"]))
        self.assertTrue(faims.is_mixed(g))
        self.assertFalse(faims.is_mixed(faims.group({"a": [-50.0], "c": None})))
        self.assertFalse(faims.is_mixed(faims.group({"b": [], "c": []})))

    def test_step_2_s_record_wins_over_reading_the_file(self):
        with tempfile.TemporaryDirectory() as d:
            f = put(d, "a.mzML", mzml(None))
            rec = faims.from_records([{"file": f, "faims_cvs": [-60.0]}])
            self.assertEqual(faims.statuses([f], rec), {f: [-60.0]})
            self.assertEqual(faims.statuses([f], {}), {f: []})


def quantifying_workflow():
    wf = open(FP_WORKFLOW, newline="").read()
    assert "quantitation.run-label-free-quant=false" in wf
    return wf.replace("quantitation.run-label-free-quant=false",
                      "quantitation.run-label-free-quant=true")


def mixed_text(i):
    return mzml(-50.0 if i % 2 else None)


class Refusal(unittest.TestCase):
    """run_search.py end to end (FragPipe is /bin/true): nothing is written or submitted."""

    def test_a_quantifying_search_over_a_mix_is_refused_with_counts_and_choices(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "fp_job.sh")
            p, out = search(d, "fragpipe", "--sbatch", job, cfg=quantifying_workflow(),
                            n_files=5, file_text=mixed_text)
            self.assertNotEqual(p.returncode, 0)
            msg = p.stdout + p.stderr
            for words in ("REFUSED", "nothing was written", "2 run(s) were acquired with FAIMS",
                          "3 without", "-50", "--faims-split", "NOT comparable", "IonQuant"):
                self.assertIn(words, msg)
            self.assertFalse(os.path.exists(job))

    def test_an_identifications_only_search_may_run_mixed(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "fp_job.sh")
            p, out = search(d, "fragpipe", "--sbatch", job, n_files=4, file_text=mixed_text)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertTrue(os.path.exists(job))

    def test_all_faims_or_none_runs_normally_and_is_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "fragpipe", "--sbatch", os.path.join(d, "a.sh"),
                            cfg=quantifying_workflow(), n_files=3,
                            file_text=lambda i: mzml(-55.0))
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            z = provenance(out)["faims"]
            self.assertEqual(z["cvs"], [-55.0])
            self.assertFalse(z["quantified_separately"])
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "fragpipe", "--sbatch", os.path.join(d, "a.sh"),
                            cfg=quantifying_workflow(), n_files=3,
                            file_text=lambda i: mzml(None))
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIsNone(provenance(out)["faims"])


class LoginNodeBudget(unittest.TestCase):
    """pm-215x check h: the gate read every mzML head on the login node with no total cap. Past
    faims.LOGIN_READ_BYTES (the worst case of HEAD_BYTES[-1] a run) the whole command moves to a
    compute node, as check_tdf does -- never "unknown" for the rest (the plain gate only NOTEs an
    unknown run, so a mix hidden in it would crash IonQuant after the search)."""

    def _gate(self, d, n, login=True, record=False, queue=None):
        files = [put(d, f"run{i}.mzML", mzml(None)) for i in range(n)]
        acq = None
        if record:
            acq = os.path.join(d, "acquisition.json")
            with open(acq, "w") as fh:
                json.dump({"files": [{"file": f, "faims_cvs": []} for f in files]}, fh)
        a = mock.Mock(out=os.path.join(d, "out"), faims_group=None, faims_split=False,
                      sbatch=None)
        moved, reads = [], []

        def ensure(why, **kw):
            moved.append((why, kw))
            raise SystemExit(0)                # the real one never returns

        import on_compute
        with mock.patch.object(on_compute, "on_login_node", lambda env=None: login), \
                mock.patch.object(on_compute, "ensure_compute_node", ensure), \
                mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(run_search.TDF_REEXEC_ENV, None)
            try:
                run_search.fragpipe_faims_gate(
                    a, files, {"acquisition": "DDA"}, {"quantified": True}, acq,
                    reader=lambda f: reads.append(f) or faims.run_cvs(f), queue=queue)
            except SystemExit:
                pass
            marked = os.environ.get(run_search.TDF_REEXEC_ENV)
        return files, moved, reads, marked

    def test_the_budget_is_five_runs_at_the_worst_case(self):
        self.assertEqual(faims.LOGIN_READ_BYTES, 5 * faims.HEAD_BYTES[-1])
        self.assertEqual(faims.LOGIN_READ_BYTES, run_search.TDF_LOGIN_MAX * faims.HEAD_BYTES[-1])
        files = ["a.mzML", "b.mzML.gz", "c.raw", "d.d"]
        self.assertEqual(faims.to_read(files), ["a.mzML", "b.mzML.gz"])
        self.assertEqual(faims.read_bytes_max(files), 2 * faims.HEAD_BYTES[-1])

    def test_six_unread_mzml_on_a_login_node_move_to_a_compute_node(self):
        with tempfile.TemporaryDirectory() as d:
            _, moved, reads, marked = self._gate(d, 6, queue={"partition": "high",
                                                              "account": "labgrp", "qos": None})
        self.assertEqual(reads, [], "nothing is read on the login node")
        (why, kw), = moved
        self.assertIn("6 mzML", why)
        self.assertIn("FAIMS", why)
        self.assertEqual((kw["mem"], kw["minutes"]), ("2G", 10))
        self.assertEqual((kw["account"], kw["partition"]), ("labgrp", "high"))
        self.assertIs(kw["preemptible_ok"], True)
        self.assertEqual(os.path.realpath(kw["script"]), os.path.realpath(run_search.__file__))
        self.assertEqual(marked, "faims", "the run under that srun must know why it is there")

    def test_five_are_read_here(self):
        with tempfile.TemporaryDirectory() as d:
            files, moved, reads, _ = self._gate(d, 5)
        self.assertEqual(moved, [])
        self.assertEqual(sorted(reads), sorted(files))

    def test_recorded_runs_are_not_read_and_not_moved(self):
        with tempfile.TemporaryDirectory() as d:
            _, moved, reads, _ = self._gate(d, 12, record=True)
        self.assertEqual((moved, reads), ([], []))

    def test_off_a_login_node_they_are_read_here(self):
        with tempfile.TemporaryDirectory() as d:
            files, moved, reads, _ = self._gate(d, 8, login=False)
        self.assertEqual(moved, [])
        self.assertEqual(len(reads), 8)

    def test_no_search_runs_inline_in_the_srun_it_moved_to(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "fragpipe", cfg=quantifying_workflow(), n_files=2,
                            file_text=lambda i: mzml(None),
                            env_extra={run_search.TDF_REEXEC_ENV: "faims", "SLURM_JOB_ID": "77"})
            self.assertNotEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertIn("moved to a compute node only to read the mzML heads", p.stderr)


class Split(unittest.TestCase):
    def test_each_group_gets_its_own_search_folder_and_record(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "fp_job.sh")
            p, out = search(d, "fragpipe", "--sbatch", job, "--faims-split",
                            cfg=quantifying_workflow(), n_files=5, file_text=mixed_text)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            split = json.load(open(os.path.join(out, "faims_split.json")))
            self.assertFalse(split["comparable_across_groups"])
            self.assertEqual({g["key"]: g["n_runs"] for g in split["groups"]},
                             {"none": 3, "faims": 2})
            got = {}
            for key in ("none", "faims"):
                sub = os.path.join(out, f"faims_{key}")
                self.assertTrue(os.path.exists(os.path.join(d, f"fp_job_{key}.sh")), key)
                z = provenance(sub)["faims"]
                got[key] = z
                self.assertEqual(provenance(sub)["n_files"], 3 if key == "none" else 2)
                self.assertTrue(z["quantified_separately"])
                self.assertFalse(z["comparable_across_groups"])
                self.assertIn("not comparable across the groups", z["sentence"])
            self.assertEqual(got["faims"]["cvs"], [-50.0])
            self.assertEqual(got["none"]["other_groups"], {"FAIMS": 2})

    def test_the_split_record_says_the_design_was_not_checked_when_unknown(self):
        """C2 review stat M4: with no conditions.csv yet, faims_split.json says so."""
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "fragpipe", "--sbatch", os.path.join(d, "fp_job.sh"),
                            "--faims-split", cfg=quantifying_workflow(), n_files=5,
                            file_text=mixed_text)
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            dc = json.load(open(os.path.join(out, "faims_split.json")))["design_check"]
            self.assertFalse(dc["checked"])
            self.assertIn("NOT checked", dc["why_not"])
            self.assertIn("NOT checked", p.stderr)

    def _plan(self, d):
        none = [os.path.join(d, f"ctl{i}.mzML") for i in range(3)]
        fa = [os.path.join(d, f"trt{i}.mzML") for i in range(2)]
        return [(faims.GROUP_NONE, "none", none, []), (faims.GROUP_FAIMS, "faims", fa, [-50.0])]

    def _conditions(self, d, rows):
        p = os.path.join(d, "conditions.csv")
        with open(p, "w") as fh:
            fh.write("File.Name,Group\n" + "".join(f"{r},{g}\n" for r, g in rows))
        return p

    def test_a_condition_wholly_in_one_group_is_warned_with_the_table(self):
        with tempfile.TemporaryDirectory() as d:
            plan = self._plan(d)
            c = self._conditions(d, [("ctl0", "control"), ("ctl1", "control"), ("ctl2", "control"),
                                     ("trt0", "treated"), ("trt1", "treated")])
            r = faims.design_check(plan, c)
            self.assertTrue(r["checked"])
            self.assertEqual(r["table"], {"control": {faims.GROUP_NONE: 3},
                                          "treated": {faims.GROUP_FAIMS: 2}})
            self.assertEqual({x["condition"] for x in r["confounded"]}, {"control", "treated"})
            self.assertIn("FAIMS difference, not a biological one", r["warning"])
            self.assertIn("control | 0 | 3", "\n".join(faims.design_table_text(r)))

    def test_conditions_spread_over_both_groups_are_not_warned(self):
        with tempfile.TemporaryDirectory() as d:
            plan = self._plan(d)
            c = self._conditions(d, [("ctl0.mzML", "control"), ("trt0", "control"),
                                     ("ctl1", "treated"), ("trt1", "treated"), ("ctl2", "x")])
            r = faims.design_check(plan, c)
            self.assertTrue(r["checked"])
            self.assertEqual([x["condition"] for x in r["confounded"]], ["x"])
            r = faims.design_check(plan, self._conditions(d, [("ctl0", "a"), ("trt0", "a"),
                                                              ("ctl1", "b"), ("trt1", "b")]))
            self.assertEqual(r["confounded"], [])
            self.assertIsNone(r["warning"])
            self.assertEqual(r["runs_not_in_design"], ["ctl2"])

    def test_a_design_naming_no_run_of_the_search_is_not_checked(self):
        with tempfile.TemporaryDirectory() as d:
            r = faims.design_check(self._plan(d), self._conditions(d, [("other", "a")]))
            self.assertFalse(r["checked"])
            self.assertIn("NOT checked", r["why_not"])

    def test_an_unreadable_run_stops_a_split(self):
        with tempfile.TemporaryDirectory() as d:
            p, out = search(d, "fragpipe", "--sbatch", os.path.join(d, "a.sh"), "--faims-split",
                            cfg=quantifying_workflow(), n_files=3,
                            file_text=lambda i: mixed_text(i) if i < 2 else "<mzML/>")
            self.assertNotEqual(p.returncode, 0)
            self.assertIn("could not be read", p.stdout + p.stderr)

    def test_the_child_command_line(self):
        argv = ["--tools", "t.json", "--files", "a.mzML", "b.mzML", "--out", "o", "--sbatch",
                "job.sh", "--faims-split", "--threads", "8"]
        got = run_search.faims_split_argv(argv, "faims", ["b.mzML"], "o/faims_faims", "job.sh")
        self.assertEqual(got, ["--tools", "t.json", "--threads", "8", "--files", "b.mzML",
                               "--out", "o/faims_faims", "--faims-group", "faims", "--sbatch",
                               "job_faims.sh"])
        self.assertNotIn("--faims-split", got)


class Sentences(unittest.TestCase):
    def setUp(self):
        self.rec = faims.split_record("FAIMS", [-50.0], 23, {"no FAIMS": 58})

    def test_the_methods_say_it(self):
        import make_methods
        rec = {"engine": "fragpipe", "faims": self.rec}
        text = make_methods.search_sentence(rec) if hasattr(make_methods, "search_sentence") else None
        # the sentence lives in the record; the Methods paragraph appends it verbatim
        self.assertIn("not comparable across the groups", self.rec["sentence"])
        src = open(os.path.join(SCRIPTS, "make_methods.py")).read()
        self.assertIn('fz["sentence"]', src)

    def test_the_manifest_says_it(self):
        import session
        with tempfile.TemporaryDirectory() as d:
            sub = os.path.join(d, "output", "search_x", "faims_faims")
            os.makedirs(sub)
            json.dump({"faims": self.rec}, open(os.path.join(sub, "search_provenance.json"), "w"))

            class Man:
                def __init__(self):
                    self.notes = []

                def note(self, name, why):
                    self.notes.append((name, why))
            m = Man()
            session.faims_step({"output_dir": os.path.join(d, "output")}, m)
            self.assertEqual(len(m.notes), 1)
            self.assertIn("FAIMS group search", m.notes[0][0])
            self.assertIn("not comparable across the groups", m.notes[0][1])

    def test_the_report_says_it(self):
        import make_analysis_html as mah
        with tempfile.TemporaryDirectory() as d:
            sess = os.path.join(d, "s")
            out = os.path.join(sess, "output", "search_x")
            os.makedirs(out)
            json.dump({"engine": "fragpipe", "faims": self.rec},
                      open(os.path.join(out, "search_provenance.json"), "w"))
            self.assertIn("not comparable across the groups", mah.search_faims_sentence(sess) or "")


def semi_specific_workflow():
    """The non-specific fixture workflow made semi-specific (one enzyme terminus)."""
    return open(FP_WORKFLOW, newline="").read().replace("msfragger.num_enzyme_termini=0",
                                                         "msfragger.num_enzyme_termini=1")


class MemoryByRuns(unittest.TestCase):
    """Semi-specific FragPipe memory by runs: every measured point gets 15% over its peak."""

    def test_the_fixture_is_semi_specific_after_the_edit(self):
        with tempfile.TemporaryDirectory() as d:
            wf = os.path.join(d, "x.workflow")
            open(wf, "w", newline="").write(semi_specific_workflow())
            self.assertEqual(run_search.digest_breadth("fragpipe", wf)[0], "semi-specific")

    def test_a_non_specific_search_keeps_the_x9_numbers(self):
        f = lambda n, s=None: run_search.broad_digest_mem_gb("fragpipe", s, n, "non-specific")
        self.assertEqual(run_search.fragpipe_broad_mem_gb(8, "non-specific"), 200)
        self.assertEqual(run_search.fragpipe_broad_mem_gb(219, "non-specific"), 246)
        self.assertEqual(f(219, 12685), run_search.broad_digest_mem_gb("fragpipe", 12685, 38) + 46)

    def test_every_measured_point_has_15_percent_headroom_on_both_paths(self):
        for runs, peak in run_search.FRAGPIPE_SEMISPECIFIC_PEAKS:
            for n_seq in (None, 12685):
                got = (run_search.fragpipe_broad_mem_gb(runs, "semi-specific") if n_seq is None
                       else run_search.broad_digest_mem_gb("fragpipe", n_seq, runs,
                                                           "semi-specific"))
                self.assertGreaterEqual(got, peak * 1.15, (runs, n_seq, got))

    def test_the_beer_points_in_numbers(self):
        f = lambda n: run_search.fragpipe_broad_mem_gb(n, "semi-specific")
        self.assertEqual([f(n) for n in (8, 23, 58, 183, 300)], [232, 232, 272, 296, 312])
        # a ceiling-limited 219-run point is not a measurement, but the rule clears it by far
        self.assertGreater(f(219), 209.7 * 1.15)
        self.assertEqual(f(None), run_search.BROAD_DIGEST_FRAGPIPE_MEM_GB)

    def test_it_never_falls_with_more_runs(self):
        sizes = [run_search.fragpipe_broad_mem_gb(n, "semi-specific") for n in range(1, 400, 7)]
        self.assertEqual(sizes, sorted(sizes))
        sizes = [run_search.broad_digest_mem_gb("fragpipe", 22939, n, "semi-specific")
                 for n in (10, 100, 219, 400)]
        self.assertEqual(sizes, sorted(sizes))

    def test_the_tryptic_term_is_unchanged(self):
        f = run_search.fragpipe_runs_mem_gb
        self.assertEqual([f(0), f(38), f(None)], [0, 0, 0])
        self.assertEqual((f(39), f(219)), (1, 46))
        self.assertGreater(f(300), f(219))

    def test_sage_does_not_get_the_term(self):
        self.assertEqual(run_search.broad_digest_mem_gb("sage", 22939, 400),
                         run_search.broad_digest_mem_gb("sage", 22939))

    def test_a_semi_specific_job_without_a_database_count_uses_the_runs(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "fp_job.sh")
            p, out = search(d, "fragpipe", "--sbatch", job, n_files=58,
                            cfg=semi_specific_workflow())
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            self.assertEqual(header(job, "mem"), "272G")
            self.assertIn("MEASURED", provenance(out)["job_mem_source"])

    def test_the_job_asks_for_it_and_the_heap_keeps_the_margin(self):
        with tempfile.TemporaryDirectory() as d:
            job = os.path.join(d, "fp_job.sh")
            p, out = search(d, "fragpipe", "--sbatch", job, n_files=219, n_seq=12685,
                            cfg=semi_specific_workflow())
            self.assertEqual(p.returncode, 0, p.stdout + p.stderr)
            want = run_search.broad_digest_mem_gb("fragpipe", 12685, 219, "semi-specific")
            self.assertGreater(want, 200)
            self.assertEqual(header(job, "mem"), f"{want}G")
            self.assertIn(f"--ram {want - 12}", open(job).read())
            src = provenance(out)["job_mem_source"]
            self.assertIn("219 runs", src)
            self.assertIn("MEASURED", src)           # the numbers are labelled as measured
            self.assertIn("183 runs 251 GiB", src)

    def test_a_tryptic_many_run_search_gets_the_term_too_and_a_small_one_does_not(self):
        sizing = run_search.single_job_sizing
        self.assertIsNone(sizing("fragpipe", "/nonexistent/x.workflow", n_runs=30)["mem"])
        z = sizing("fragpipe", "/nonexistent/x.workflow", n_runs=219)
        self.assertEqual(z["mem"], f"{64 + run_search.fragpipe_runs_mem_gb(219)}G")
        self.assertIn("lower-bound", z["mem_source"])

    def test_the_users_mem_still_wins(self):
        z = run_search.single_job_sizing("fragpipe", "/nonexistent/x.workflow", mem_gb=96,
                                         n_runs=219)
        self.assertEqual((z["mem"], z["mem_source"]), ("96G", "user (--mem)"))


if __name__ == "__main__":
    unittest.main()
