"""
The data quality notes on an HT plate whose runs Core staff renamed when staging (2.15).

Two bugs on one HT plate (96 runs; the instrument's sample names were wrong for 94 of them, and
staff staged corrected names into raw/):

  (i)  The notes listed the raw-file names AS ACQUIRED, unlabelled, while the delivered raw/
       folder holds the corrected names: an agent joining by name mismatched most runs.
  (ii) The sheet check compared the CoreOmics sample list's unique_id values -- row IDs (D9001
       ...) that never appear in an HT plate's run names -- and reported "96 of 96 sample IDs
       ... match no raw file" and "96 of 96 raw files carry no sample ID": a data problem that
       is not one. An HT plate run carries the sheet's sample_name.

Now the check matches each run by the rule locate and the HT manifest use
(core_submission.run_matches: an HT plate run by unique_id in its sample field, else by
sample_name; any other run by unique_id), says which basis it used, and says in one neutral line
when no key can match by design. With a stage record that renames runs (the service project's
.core_submission.json `links`: name = the corrected name staged in raw/, target = the run as
acquired -- the one record of a corrected name, core_submission.raw_names), a run is matched by
its delivered name, and every raw name the notes list is labelled "as acquired" with the
delivered name beside it. An ordinary submission, whose run names carry the unique_id, is
checked exactly as before (pinned first).

Fake records only: PROT_9901 and its samples are made up. Stdlib only, no network.
"""
import json
import os
import shutil
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, HERE)
sys.path.insert(0, SCRIPTS)

import session_docs  # noqa: E402
import submission_report as sr  # noqa: E402
import test_submission_report as tsr  # noqa: E402  (as a module: its TestCase classes stay its own)
from test_submission_report import RUNS, make_session, notes_by_id  # noqa: E402

PROT, HEX = "PROT_9901", "9901abcdef12"
N = 12
STAGE = ".core_submission.json"


def uid(i):
    return f"D{9000 + i}"                       # a CoreOmics row ID: never in a run name


def sname(i):
    return f"XY{100 + i}"


def run(i, name=None, well=None):
    """An HT plate run: <date>_<PROT>_<speed>_<sample name>_S<plate>-<well>_1_<acq#>.d"""
    return f"20261023_{PROT}_100spd_{name or sname(i)}_S2-{well or f'A{i}'}_1_{30000 + i}.d"


def record(names=True):
    return {"schema": sr.SCHEMA, "source": sr.SOURCE_COREOMICS, "internal_id": PROT, "id": HEX,
            "samples": [{"unique_id": uid(i), "sample_name": sname(i) if names else "",
                         "condition_name": "Wing"} for i in range(1, N + 1)]}


def acquired_name(i):
    """As acquired: wells A1..A10 typed with the next well's sample (A10 with A1's); A11, A12 right."""
    return run(i, sname(i % 10 + 1)) if i <= 10 else run(i)


EXTRA_ACQ = run(13, "XY199", "B1")            # a run of no sample on the sheet, renamed too
EXTRA_DEL = run(13, "XY198", "B1")


def plate_session(tmp, raws, rec=None, renames=None):
    """A session under the submission's work dir, as stage leaves it: <work>/PROT_9901/sessions/
    <session>, the stage record in the work dir and the service project. `renames`: {as acquired:
    delivered}, each a basename; every raw file is a staged link (unchanged unless renamed)."""
    data = os.path.join(tmp, "instrument", "oct26")
    work = os.path.join(tmp, "work", "off_campus", "Example_Lab", PROT)
    svc = os.path.join(tmp, "service", "off_campus", "Example_Lab", PROT)
    os.makedirs(os.path.join(work, "sessions"))
    os.makedirs(svc)
    s = make_session(os.path.join(work, "sessions"))
    with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
        fh.write("# raw\n" + "".join(os.path.join(data, r) + "\n" for r in raws))
    sr.attach(s, sr.sanitize(rec or record()))
    if renames is not None:
        links = [{"name": renames.get(r, r), "target": os.path.join(data, r),
                  "link": os.path.join(data, r), "relative": False, "action": "create"}
                 for r in raws]
        stage = {"schema": "core_submission_stage/1", "id": HEX, "internal_id": PROT,
                 "label": PROT, "service_dir": svc, "work_dir": work,
                 "files": [x["target"] for x in links], "links": links}
        for d in (work, svc):
            with open(os.path.join(d, STAGE), "w") as fh:
                json.dump(stage, fh)
    return s


def sheet_notes(s):
    return [n for n in sr.quality_notes(sr.load(s), s) if n["id"].startswith(
        ("sheet_", "raw_", "weak_sample_ids"))]


class TestOrdinarySubmissionUnchanged(unittest.TestCase):
    """Pinned first: a submission whose run names carry the unique_id (PROT_9756's LRS96 ...)."""

    def test_unique_ids_in_the_run_names_are_checked_as_before(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = [v for k, v in RUNS.items() if k != "LRS125"] + ["08142026__60SPD_DIA-QC-HeLa_S3-A1_1_1"]
            s = tsr.TestNotes().session(tmp, raws=raws)
            n = notes_by_id(sr.load(s), s)
            self.assertEqual(n["sheet_ids_without_raw"],
                             "Sample sheet vs raw files: 1 of 30 sample IDs on the sheet match no "
                             "raw file analysed here (LRS125).")
            self.assertEqual(n["raw_without_sheet_id"],
                             "Sample sheet vs raw files: 1 of 30 raw files carry no sample ID from "
                             "the sheet (08142026__60SPD_DIA-QC-HeLa_S3-A1_1_1.d).")
            self.assertFalse({"sheet_match_basis", "sheet_not_compared", "raw_names_corrected"}
                             & set(n))

    def test_every_id_matched_says_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = tsr.TestNotes().session(tmp)
            self.assertEqual([n for n in notes_by_id(sr.load(s), s) if n.startswith(
                ("sheet_ids", "raw_", "sheet_match", "sheet_not"))], [])


    def test_a_plate_token_without_a_plate_position_is_an_ordinary_run(self):
        # review M3: an Astral-style run named <YYYYMMDD>_PROT_<n>_... carries the plate token but
        # no plate position (_S<plate>-<well>_<n>_<acq#>): it is not an HT plate run, so it is
        # matched by unique_id only, never by a sample name, and reported in the ordinary words
        with tempfile.TemporaryDirectory() as tmp:
            raws = [f"20261023_{PROT}_Astral_{sname(i)}_0{i}.raw" for i in range(1, 4)] + [
                f"20261023_{PROT}_Astral_{uid(4)}_04.raw"]
            rec = dict(record(), samples=record()["samples"][:4])
            s = plate_session(tmp, raws, rec=rec)
            n = {x["id"]: x for x in sheet_notes(s)}
            self.assertEqual(n["raw_without_sheet_id"]["text"],
                             "Sample sheet vs raw files: 3 of 4 raw files carry no sample ID from "
                             "the sheet (" + ", ".join(raws[:3]) + ").")
            self.assertEqual(n["raw_without_sheet_id"]["basis"], "unique_id")
            self.assertNotIn("sheet_match_basis", n)


class TestHTPlateMatchedBySampleName(unittest.TestCase):
    def test_row_ids_are_not_reported_as_96_of_96_unmatched(self):
        with tempfile.TemporaryDirectory() as tmp:        # staged, nothing renamed
            s = plate_session(tmp, [run(i) for i in range(1, N + 1)], renames={})
            notes = sheet_notes(s)
            ids = [n["id"] for n in notes]
            self.assertNotIn("sheet_ids_without_raw", ids)
            self.assertNotIn("raw_without_sheet_id", ids)
            self.assertEqual(ids, ["sheet_match_basis"])
            self.assertEqual(notes[0]["basis"], "sample_name")
            self.assertEqual(notes[0]["names"], "as analysed")       # no record renames any
            self.assertEqual(notes[0]["text"], (
                "Sample sheet vs raw files: matched by sample name. HT plate runs carry the "
                "sheet's sample names (sample_name), not its sample IDs (unique_id, such as "
                "D9001); every sample on the sheet has a run analysed here, and every run a "
                "sample."))

    def test_a_missing_run_and_a_stray_run_are_reported_by_sample_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = plate_session(tmp, [run(i) for i in range(1, N)] + [EXTRA_ACQ])
            n = {x["id"]: x for x in sheet_notes(s)}
            self.assertEqual(n["sheet_ids_without_raw"]["text"],
                             "Sample sheet vs raw files: 1 of 12 samples on the sheet match no raw "
                             "file analysed here by sample name (D9012 (XY112)).")
            self.assertEqual(n["raw_without_sheet_id"]["text"],
                             f"Sample sheet vs raw files: 1 of 12 raw files carry no sample name "
                             f"from the sheet ({EXTRA_ACQ}).")
            self.assertEqual({x["basis"] for x in n.values()}, {"sample_name"})
            self.assertIn("matched by sample name", n["sheet_match_basis"]["text"])
            self.assertNotIn("every sample on the sheet has a run", n["sheet_match_basis"]["text"])

    def test_no_key_can_match_by_design_is_one_neutral_line(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = plate_session(tmp, [run(i) for i in range(1, N + 1)], rec=record(names=False))
            notes = sheet_notes(s)
            self.assertEqual([n["id"] for n in notes], ["sheet_not_compared"])
            self.assertEqual(notes[0]["basis"], "none")
            text = notes[0]["text"]
            self.assertEqual(text, (
                "Sample sheet vs raw files: not compared by name. HT plate runs carry the plate "
                "position and a sample name, never the sheet's sample IDs (unique_id, such as "
                "D9001), and the sheet gives no sample names."))
            self.assertNotIn("of 12", text)


def dup_record():
    """The real plate's shape, made up: one sample name on 3 samples (QC1), another on 4 (QC2)."""
    rec = record()
    for i, smp in enumerate(rec["samples"], 1):
        if i <= 3:
            smp["sample_name"] = "QC1"
        elif i <= 7:
            smp["sample_name"] = "QC2"
    return rec


def dup_run(i):
    return run(i, "QC1" if i <= 3 else "QC2" if i <= 7 else sname(i))


class TestSharedSampleNames(unittest.TestCase):
    """review H1: two samples with one sample name give one run two owners. locate's curated_rows
    calls such a run ambiguous (staff choose, --choose); the notes claimed a full match."""

    def test_runs_of_a_shared_name_are_ambiguous_not_a_full_match(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = plate_session(tmp, [dup_run(i) for i in range(1, N + 1)], rec=dup_record())
            n = {x["id"]: x for x in sheet_notes(s)}
            self.assertNotIn("every sample on the sheet has a run", n["sheet_match_basis"]["text"])
            self.assertEqual(n["sheet_ambiguous"]["text"], (
                "Sample sheet vs raw files: 7 of 12 raw files carry a sample name that more than "
                "one sample on the sheet has (7 samples: QC1 on 3, QC2 on 4), so the names cannot "
                "tell which of those samples each of these runs is (" + ", ".join(
                    dup_run(i) for i in range(1, 7)) + " and 1 more). The Core's assignment of "
                "runs to samples (core_submission.py locate --choose) settles it."))
            self.assertEqual(n["sheet_ambiguous"]["basis"], "sample_name")
            self.assertNotIn("sheet_ids_without_raw", n)

    def test_a_sample_with_no_run_of_its_own_is_not_counted_as_matched(self):
        # the reviewer's probe: one QC1 run missing -- the notes still said "every sample ... has
        # a run"; the names cannot say which QC1 sample has none, so it is ambiguous, not matched
        with tempfile.TemporaryDirectory() as tmp:
            s = plate_session(tmp, [dup_run(i) for i in range(1, N + 1) if i != 2],
                              rec=dup_record())
            n = {x["id"]: x for x in sheet_notes(s)}
            self.assertIn("6 of 11 raw files carry a sample name that more than one sample on "
                          "the sheet has (7 samples: QC1 on 3, QC2 on 4)", n["sheet_ambiguous"]["text"])
            self.assertNotIn("every sample on the sheet has a run", n["sheet_match_basis"]["text"])

    def test_two_sample_ids_on_one_ordinary_run_are_ambiguous_too(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = ["01012026_DIA-S10-S11_1_100", "01012026_DIA-S12_1_101"]
            rec = {"schema": sr.SCHEMA, "source": sr.SOURCE_COREOMICS, "internal_id": PROT,
                   "id": HEX, "samples": [{"unique_id": u, "sample_name": "", "condition_name": "c"}
                                          for u in ("S10", "S11", "S12")]}
            s = plate_session(tmp, [r + ".d" for r in raws], rec=rec)
            n = {x["id"]: x for x in sheet_notes(s)}
            self.assertEqual(n["sheet_ambiguous"]["text"], (
                "Sample sheet vs raw files: 1 of 2 raw files carry a sample ID that more than one "
                "sample on the sheet has (2 samples: S10, S11), so the names cannot tell which of "
                "those samples each of these runs is (01012026_DIA-S10-S11_1_100.d). The Core's "
                "assignment of runs to samples (core_submission.py locate --choose) settles it."))
            self.assertEqual(n["sheet_ambiguous"]["basis"], "unique_id")


class TestWordings(unittest.TestCase):
    """review L3: the four wordings that said more than the records do."""

    def test_sample_names_equal_to_the_ids_are_matched_by_sample_id(self):
        with tempfile.TemporaryDirectory() as tmp:
            rec = record()
            for smp in rec["samples"]:
                smp["sample_name"] = smp["unique_id"]
            s = plate_session(tmp, [run(i, uid(i)) for i in range(1, N)] + [EXTRA_ACQ], rec=rec)
            n = {x["id"]: x for x in sheet_notes(s)}
            self.assertNotIn("sheet_match_basis", n)                # an ordinary match: by ID
            self.assertEqual(n["raw_without_sheet_id"]["basis"], "unique_id")
            self.assertIn("1 of 12 raw files carry no sample ID from the sheet",
                          n["raw_without_sheet_id"]["text"])

    def test_when_no_run_carries_a_sample_name_it_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = plate_session(tmp, [run(i, f"ZZ{i}") for i in range(1, 4)])
            n = {x["id"]: x for x in sheet_notes(s)}
            self.assertEqual(n["sheet_match_basis"]["text"], (
                "Sample sheet vs raw files: compared by sample name (HT plate runs carry sample "
                "names, not the sheet's sample IDs, such as D9001), and none of the runs analysed "
                "here carries a sample name from the sheet."))
            self.assertNotIn("carry the sheet's sample names", n["sheet_match_basis"]["text"])

    def test_no_example_id_is_made_up(self):
        with tempfile.TemporaryDirectory() as tmp:
            rec = record(names=False)
            for smp in rec["samples"]:
                smp["unique_id"] = ""
            s = plate_session(tmp, [run(i) for i in range(1, 4)], rec=rec)
            for x in sheet_notes(s):
                self.assertNotIn("such as ?", x["text"])
                self.assertNotIn("?", x["text"])

    def test_a_mixed_assignment_names_both_keys_for_an_ht_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = [run(i) for i in range(1, N + 1)] + [EXTRA_ACQ]
            s = plate_session(tmp, raws)
            with open(os.path.join(s, "input", "locate_decisions.json"), "w") as fh:
                json.dump({"assignments": [{"sample": uid(1), "file": run(1)}]}, fh)
            n = {x["id"]: x for x in sheet_notes(s)}
            self.assertEqual(n["raw_without_sheet_id"]["basis"], "mixed")
            self.assertIn("were assigned to no sample and carry no sample ID or sample name from "
                          "the sheet", n["raw_without_sheet_id"]["text"])


class TestNoStageRecord(unittest.TestCase):
    def test_without_the_stage_record_no_full_match_is_claimed(self):
        # review M2: the names analysed may be the instrument's wrong ones; the notes say so
        with tempfile.TemporaryDirectory() as tmp:
            s = plate_session(tmp, [run(i) for i in range(1, N + 1)])
            n = notes_by_id(sr.load(s), s)
            self.assertIn("stage_record_not_found", n)
            self.assertTrue(n["sheet_match_basis"].endswith("(unique_id, such as D9001)."))


class TestRenamedRunsAsAcquired(unittest.TestCase):
    def renamed(self, tmp):
        raws = [acquired_name(i) for i in range(1, N + 1)] + [EXTRA_ACQ]
        ren = {acquired_name(i): run(i) for i in range(1, 11)}
        ren[EXTRA_ACQ] = EXTRA_DEL
        return plate_session(tmp, raws, renames=ren), raws

    def test_every_raw_name_listed_is_labelled_as_acquired_with_its_delivered_name(self):
        with tempfile.TemporaryDirectory() as tmp:
            s, raws = self.renamed(tmp)
            n = {x["id"]: x for x in sr.quality_notes(sr.load(s), s)}
            self.assertEqual(n["raw_without_sheet_id"]["text"],
                             f"Sample sheet vs raw files: 1 of 13 raw files carry no sample name "
                             f"from the sheet ({EXTRA_ACQ} (as acquired; delivered as "
                             f"{EXTRA_DEL})).")
            self.assertIn(
                "Raw file names: Core staff corrected the names of 11 of the 13 runs analysed here "
                "when staging them, and the delivered raw/ folder holds those runs under the "
                "corrected names. input/raw_files.txt, the search and its tables use the names as "
                "acquired. Raw-file names in these notes are labelled as acquired or as "
                "delivered; join on the delivered name or the plate position, never on a name as "
                "acquired.", n["raw_names_corrected"]["text"])
            self.assertNotIn("beside it", n["raw_names_corrected"]["text"])
            basis = n["sheet_match_basis"]
            self.assertEqual(basis["names"], "delivered")
            self.assertIn("Runs Core staff renamed when staging were matched by the name they "
                          "are delivered under.", basis["text"])

    def test_a_renamed_run_is_the_sample_its_delivered_name_says(self):
        with tempfile.TemporaryDirectory() as tmp:
            s, raws = self.renamed(tmp)
            p = sr.paths_for(s)
            owners = sr.run_assignment(sr.load(s), p, raws)[0]
            self.assertEqual(owners[acquired_name(1)], {uid(1)})     # acquired as XY102: XY101's
            self.assertEqual(owners[acquired_name(10)], {uid(10)})   # acquired as XY101
            self.assertEqual(owners[acquired_name(12)], {uid(12)})   # not renamed
            self.assertEqual(owners[EXTRA_ACQ], set())

    def test_agents_md_carries_the_labels(self):
        with tempfile.TemporaryDirectory() as tmp:
            s, raws = self.renamed(tmp)
            with open(session_docs.write_docs(s)["agents"], encoding="utf-8") as fh:
                agents = fh.read()
            notes = agents.split("## Data quality notes to respect", 1)[1]
            self.assertIn(f"{EXTRA_ACQ} (as acquired; delivered as {EXTRA_DEL})", notes)
            self.assertIn("Raw file names: Core staff corrected the names of 11 of the 13 runs",
                          notes)
            self.assertNotIn("of 12 sample IDs", notes)

    def test_a_name_left_as_it_was_is_still_said_to_be_as_acquired(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = [acquired_name(i) for i in range(1, N + 1)] + [EXTRA_ACQ]
            ren = {acquired_name(i): run(i) for i in range(1, 11)}   # the stray run kept its name
            s = plate_session(tmp, raws, renames=ren)
            n = notes_by_id(sr.load(s), s)
            self.assertIn(f"({EXTRA_ACQ} (as acquired and as delivered))",
                          n["raw_without_sheet_id"])

    def test_without_a_stage_record_or_a_rename_nothing_is_labelled(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = [run(i) for i in range(1, N)] + [EXTRA_ACQ]
            for ren in (None, {}):
                with self.subTest(renames=ren), tempfile.TemporaryDirectory() as t2:
                    s = plate_session(t2, raws, renames=ren)
                    n = notes_by_id(sr.load(s), s)
                    self.assertNotIn("raw_names_corrected", n)
                    self.assertNotIn("as acquired", n["raw_without_sheet_id"])
                    self.assertIn(f"({EXTRA_ACQ}).", n["raw_without_sheet_id"])

    def test_another_submissions_stage_record_is_not_read(self):
        with tempfile.TemporaryDirectory() as tmp:
            s, raws = self.renamed(tmp)
            work = os.path.dirname(os.path.dirname(s))
            with open(os.path.join(work, STAGE)) as fh:
                stage = json.load(fh)
            stage.update(id="0000aaaabbbb", internal_id="PROT_9902")
            with open(os.path.join(work, STAGE), "w") as fh:
                json.dump(stage, fh)
            n = notes_by_id(sr.load(s), s)
            self.assertNotIn("raw_names_corrected", n)
            # ... and the notes say so (meaningful on 2.14.2 too, which never looked: no such note)
            self.assertEqual(n["stage_record_not_found"], (
                "Raw file names: no stage record of PROT_9901 (.core_submission.json, written by "
                "core_submission.py stage) was found in the folders above this session, so whether "
                "Core staff corrected run names when staging is not known here; the names in these "
                "notes are as analysed. The stage record above it is PROT_9902's."))


class TestWhichSideOfTheRename(unittest.TestCase):
    """review H2 (direction), M2 (one lookup, damage reported), L4 (only the runs analysed)."""

    def swap_session(self, tmp, links_on_disk, extra_rename=False):
        """Samples S01..S03 (IDs in the names: an ordinary submission). Staff found S01 and S02
        acquired under each other's names and staged them swapped; the session searched the
        corrected links (staff's list), so input/raw_files.txt holds the DELIVERED names."""
        acq = {"01012026_DIA-S01_1_1.d": "01012026_DIA-S02_1_1.d",      # as acquired -> delivered
               "01012026_DIA-S02_1_1.d": "01012026_DIA-S01_1_1.d",
               "01012026_DIA-S03_1_3.d": "01012026_DIA-S03_1_3.d"}
        if extra_rename:                                   # a run renamed without a swap
            acq["01012026_DIA-X04_1_4.d"] = "01012026_DIA-S04_1_4.d"
        rec = {"schema": sr.SCHEMA, "source": sr.SOURCE_COREOMICS, "internal_id": PROT, "id": HEX,
               "samples": [{"unique_id": u, "sample_name": "", "condition_name": "c"}
                           for u in ("S01", "S02", "S03", "S04")]}
        corrected = os.path.join(tmp, "corrected")
        raws = [os.path.join(corrected, d) for d in acq.values()]
        s = plate_session(tmp, list(acq), rec=rec, renames=acq)
        with open(os.path.join(s, "input", "raw_files.txt"), "w") as fh:
            fh.write("# raw\n" + "".join(r + "\n" for r in raws))
        data = os.path.join(tmp, "instrument", "oct26")
        if links_on_disk:
            os.makedirs(data, exist_ok=True)
            os.makedirs(corrected)
            for a, d in acq.items():
                os.makedirs(os.path.join(data, a))
                os.symlink(os.path.join(data, a), os.path.join(corrected, d))
        return s, raws

    def test_a_swap_read_through_its_links_is_each_runs_own_sample(self):
        with tempfile.TemporaryDirectory() as tmp:
            s, raws = self.swap_session(tmp, links_on_disk=True)
            owners = sr.run_assignment(sr.load(s), sr.paths_for(s), raws)[0]
            by_base = {os.path.basename(k): v for k, v in owners.items()}
            self.assertEqual(by_base["01012026_DIA-S02_1_1.d"], {"S02"})   # never S01
            self.assertEqual(by_base["01012026_DIA-S01_1_1.d"], {"S01"})
            n = notes_by_id(sr.load(s), s)
            self.assertIn("use the corrected (delivered) names", n["raw_names_corrected"])

    def test_a_swap_with_nothing_to_tell_the_side_is_matched_to_no_sample(self):
        # a laptop copy: no files to follow, and every renamed run is in the swap
        with tempfile.TemporaryDirectory() as tmp:
            s, raws = self.swap_session(tmp, links_on_disk=False)
            owners = sr.run_assignment(sr.load(s), sr.paths_for(s), raws)[0]
            by_base = {os.path.basename(k): v for k, v in owners.items()}
            self.assertEqual(by_base["01012026_DIA-S02_1_1.d"], set())
            self.assertEqual(by_base["01012026_DIA-S01_1_1.d"], set())
            n = notes_by_id(sr.load(s), s)
            self.assertIn("2 run(s) carry a name that is one run's name as acquired and another's "
                          "delivered name, and nothing here says which: they are matched to no "
                          "sample.", n["raw_names_corrected"])
            self.assertNotIn("raw_without_sheet_id", n)      # not "carry no sample ID" either
            # review2 L4: S01 and S02 are undetermined, not missing; S04 has no run at all
            self.assertEqual(n["sheet_ids_without_raw"], (
                "Sample sheet vs raw files: 1 of 4 sample IDs on the sheet match no raw file "
                "analysed here (S04)."))
            self.assertIn("2 of 4 samples on the sheet cannot be told from the runs analysed "
                          "here", n["sheet_undetermined"])
            self.assertIn("(S01, S02).", n["sheet_undetermined"])

    def test_the_sessions_other_runs_tell_the_side_of_a_swap(self):
        with tempfile.TemporaryDirectory() as tmp:
            s, raws = self.swap_session(tmp, links_on_disk=False, extra_rename=True)
            owners = sr.run_assignment(sr.load(s), sr.paths_for(s), raws)[0]
            by_base = {os.path.basename(k): v for k, v in owners.items()}
            self.assertEqual(by_base["01012026_DIA-S02_1_1.d"], {"S02"})
            self.assertEqual(by_base["01012026_DIA-S04_1_4.d"], {"S04"})
            smap = sr.run_assignment(sr.load(s), sr.paths_for(s), raws)[2]["smap"]
            self.assertEqual(sr._run_label(raws[-1], smap), "01012026_DIA-S04_1_4.d (as "
                             "delivered; acquired as 01012026_DIA-X04_1_4.d)")

    def test_a_session_of_runs_never_renamed_gets_no_labels(self):
        # review L4: a session searching only the runs staff left alone
        with tempfile.TemporaryDirectory() as tmp:
            ren = {acquired_name(i): run(i) for i in range(1, 11)}
            s = plate_session(tmp, [acquired_name(11), acquired_name(12), EXTRA_ACQ],
                              renames=ren)
            n = {x["id"]: x for x in sr.quality_notes(sr.load(s), s)}
            self.assertNotIn("raw_names_corrected", n)
            self.assertEqual(n["raw_without_sheet_id"]["names"], "as analysed")
            self.assertNotIn("as acquired", n["raw_without_sheet_id"]["text"])
            self.assertNotIn("renamed", n["sheet_match_basis"]["text"])

    def test_a_damaged_stage_record_is_said_and_no_full_match_is_claimed(self):
        with tempfile.TemporaryDirectory() as tmp:
            raws = [acquired_name(i) for i in range(1, N + 1)]
            s = plate_session(tmp, raws, renames={acquired_name(i): run(i) for i in range(1, 11)})
            work = os.path.dirname(os.path.dirname(s))
            with open(os.path.join(work, STAGE), "w") as fh:
                fh.write("{")
            n = notes_by_id(sr.load(s), s)
            self.assertEqual(n["stage_record_unreadable"], (
                "Raw file names: the stage record of PROT_9901 above this session "
                "(.core_submission.json) cannot be read, so whether Core staff corrected run names "
                "when staging is not known here; the names in these notes are as analysed. Repair "
                "it, or stage again (core_submission.py stage --apply)."))
            self.assertNotIn("every sample on the sheet has a run", n["sheet_match_basis"])

    def test_deliver_reports_a_damaged_record_instead_of_skipping_it(self):
        import core_submission as cs
        with tempfile.TemporaryDirectory() as tmp:
            s = plate_session(tmp, [run(1)], renames={})
            work = os.path.dirname(os.path.dirname(s))
            with open(os.path.join(work, STAGE), "w") as fh:
                fh.write("[1]")
            rec, where, foreign = cs.find_service_record({"id": HEX, "internal_id": PROT}, s, None,
                                                         search_tree=False)
            self.assertIsNone(rec)
            self.assertIn(STAGE + cs.STAGE_RECORD_UNREADABLE, where)



class TestRound2NameMap(unittest.TestCase):
    """review2 L6 (the delivery README lists RAW_NAMES.tsv), L7 (a session outside its staged
    project: attach --service-project copies the map, so the notes and deliver agree)."""

    def test_the_delivery_readme_lists_the_name_map(self):
        import core_submission as cs
        summary = {"internal_id": PROT, "id": HEX, "pi": {}}
        for delivered, want in (({"Analysis_Report.html", cs.RAW_NAME_MAP}, True),
                                ({"Analysis_Report.html"}, False)):
            text = cs.build_readme(summary, delivered, [], "analysis")
            self.assertEqual(f"| `{cs.RAW_NAME_MAP}` |" in text, want, text)

    def test_attach_with_the_service_project_writes_the_map_for_a_session_outside_it(self):
        import subprocess
        with tempfile.TemporaryDirectory() as tmp:
            ren = {acquired_name(i): run(i) for i in range(1, 4)}
            s = plate_session(tmp, [acquired_name(i) for i in range(1, 4)], renames=ren)
            outside = os.path.join(tmp, "elsewhere", "session")
            shutil.copytree(s, outside)
            svc = os.path.join(tmp, "service", "off_campus", "Example_Lab", PROT)
            rec = os.path.join(s, "input", "submission.json")
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "submission_report.py"),
                                "attach", "--session", outside, "--record", rec, "--replace"],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertFalse(os.path.exists(os.path.join(outside, "input", "raw_names.tsv")))
            r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "submission_report.py"),
                                "attach", "--session", outside, "--record", rec, "--replace",
                                "--service-project", svc],
                               capture_output=True, text=True, timeout=120)
            self.assertEqual(r.returncode, 0, r.stdout + r.stderr)
            self.assertEqual(json.loads(r.stdout)["paths"]["raw_names"]["level"]
                             if "paths" in json.loads(r.stdout) else "OK", "OK")
            self.assertTrue(os.path.exists(os.path.join(outside, "input", "raw_names.tsv")))
            n = notes_by_id(sr.load(outside), outside)
            self.assertNotIn("stage_record_not_found", n)
            self.assertIn("Core staff corrected the names of 3 of the 3 runs",
                          n["raw_names_corrected"])


class TestTheNameMapTravels(unittest.TestCase):
    """review M4 (a join on something delivered) and H4b (a report built on a laptop from the
    step-9 pull equals one built on HIVE)."""

    def renamed(self, tmp):
        raws = [acquired_name(i) for i in range(1, N + 1)] + [EXTRA_ACQ]
        ren = {acquired_name(i): run(i) for i in range(1, 11)}
        ren[EXTRA_ACQ] = EXTRA_DEL
        return plate_session(tmp, raws, renames=ren)

    def test_write_name_map_writes_every_pair_and_nothing_else(self):
        import core_submission as cs
        with tempfile.TemporaryDirectory() as tmp:
            s = self.renamed(tmp)
            self.assertEqual(sr.write_name_map(s), ("OK", "13 staged runs, 11 renamed by Core "
                                                          "staff when staging"))
            with open(os.path.join(s, "input", "raw_names.tsv"), encoding="utf-8") as fh:
                text = fh.read()
            rows = cs.parse_raw_name_map(text)
            self.assertEqual(rows[0], (acquired_name(1), run(1), "S2-A1"))
            self.assertEqual(rows[-1], (EXTRA_ACQ, EXTRA_DEL, "S2-B1"))
            self.assertEqual(len(rows), 13)
            for bad in ("/", "coreomics", HEX, "service", "@"):
                self.assertNotIn(bad, text)
            # nothing renamed now: the old map goes (it would describe renames no longer there)
            work = os.path.dirname(os.path.dirname(s))
            for d in (work, os.path.join(tmp, "service", "off_campus", "Example_Lab", PROT)):
                with open(os.path.join(d, STAGE)) as fh:
                    st = json.load(fh)
                for ln in st["links"]:
                    ln["name"] = os.path.basename(ln["target"])
                with open(os.path.join(d, STAGE), "w") as fh:
                    json.dump(st, fh)
            self.assertEqual(sr.write_name_map(s), (None, "the stage record renames no run"))
            self.assertFalse(os.path.exists(os.path.join(s, "input", "raw_names.tsv")))

    def test_a_laptop_copy_with_the_map_says_what_the_hive_session_says(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = self.renamed(tmp)
            self.assertEqual(sr.write_name_map(s)[0], "OK")
            hive = sr.quality_notes(sr.load(s), s)
            laptop = os.path.join(tmp, "laptop", "core", "PROT_9901", "session")
            shutil.copytree(s, laptop)                     # no stage record above this copy
            copy = sr.quality_notes(sr.load(laptop), laptop)
            self.assertEqual([(n["id"], n["text"]) for n in copy],
                             [(n["id"], n["text"]) for n in hive])
            self.assertIn("raw_names_corrected", {n["id"] for n in copy})
            self.assertIn("input/raw_names.tsv (delivered as RAW_NAMES.tsv, beside the report) "
                          "lists every pair, with the plate position: join on it.",
                          {n["id"]: n["text"] for n in copy}["raw_names_corrected"])
            os.remove(os.path.join(laptop, "input", "raw_names.tsv"))   # pulled without it
            n = notes_by_id(sr.load(laptop), laptop)
            self.assertIn("stage_record_not_found", n)
            self.assertNotIn("raw_names_corrected", n)

    def test_a_damaged_map_is_said(self):
        with tempfile.TemporaryDirectory() as tmp:
            s = plate_session(tmp, [run(i) for i in range(1, 4)])
            with open(os.path.join(s, "input", "raw_names.tsv"), "w") as fh:
                fh.write("not a table\n")
            n = notes_by_id(sr.load(s), s)
            self.assertIn("input/raw_names.tsv cannot be read either (ValueError: not a raw-name "
                          "table", n["stage_record_not_found"])


import test_core_submission as tcs  # noqa: E402  (its deliver fixture; a module: no tests rerun)


class TestDeliverTheNameMap(tcs.DeliverBase):
    rec_kwargs = {"raw_only": True}

    def test_a_renamed_run_is_delivered_with_the_name_map(self):
        wrong = tcs.raw(self.tmp, "tTOF_HT", "mar25", "wrong_sample.d")
        fixed = os.path.join(self.tmp, "corrected", "S07.d")
        os.makedirs(os.path.dirname(fixed))
        os.symlink(wrong, fixed)
        flist = tcs.write(self.tmp, "files_corrected.txt", "\n".join(self.raw_files + [fixed]) + "\n")
        rc, st, pr = tcs.run(["stage", "--summary", self.summary, "--files", flist, "--apply"],
                             self.env)
        self.assertEqual(rc, 0, pr.stdout + pr.stderr)
        rc, out, pr = tcs.run(["deliver", "--summary", self.summary, "--label", "raw", "--apply"],
                              self.env)
        self.assertEqual(rc, 0, pr.stdout + pr.stderr)
        d = os.path.join(self.share, "PROT_0807_raw")
        import core_submission as cs
        rows = cs.parse_raw_name_map(tcs.read(os.path.join(d, "RAW_NAMES.tsv")))
        self.assertIn(("wrong_sample.d", "S07.d", ""), rows)
        self.assertEqual(len(rows), 3)
        self.assertIn("[OK] RAW_NAMES.tsv -- the raw files' names as acquired and as delivered "
                      "(3 runs, 1 renamed)", tcs.read(os.path.join(d, "MANIFEST.txt")))
        self.assertEqual(out["raw_name_map"], {"file": "RAW_NAMES.tsv", "runs": 3, "renamed": 1})

    def test_nothing_renamed_delivers_no_name_map(self):
        rc, out, pr = tcs.run(["deliver", "--summary", self.summary, "--label", "raw", "--apply"],
                              self.env)
        self.assertEqual(rc, 0, pr.stdout + pr.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.share, "PROT_0807_raw", "RAW_NAMES.tsv")))
        self.assertIsNone(out["raw_name_map"])


if __name__ == "__main__":
    unittest.main()
