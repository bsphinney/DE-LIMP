#!/usr/bin/env python3
"""
The Core's LIMS sample sheet (.xlsx) and its per-replicate condition column.

Staff, 2026-09-25 and 09-28: the conditions file was the LIMS export
PROT_####.samples.<date>.xlsx (internal_id, internal_notes, sample_name, unique_id, condition_name,
amt_2_inject). collect_conditions.py had no .xlsx reader and the pipeline env has no openpyxl, so
it was parsed by hand; --map then rejected the columns; and condition_name names each REPLICATE
(X_mix_1 .. X_mix_5), so taken literally it is ten singleton groups and no DE.

Pinned here: the stdlib .xlsx reader (shared, rich and inline strings, numbers, sparse cells,
hidden sheets, --sheet); the LIMS headers (unique_id / condition_name; the identifier column that
names the runs is used; bookkeeping columns skipped and said); a per-replicate column read as
conditions + replicate numbers only when unambiguous, otherwise ASKED (--replicate-labels).
The workbook is synthetic, built here; no real submission, client or run names.
"""
import csv
import json
import os
import subprocess
import sys
import tempfile
import unittest
import zipfile
from xml.sax.saxutils import escape

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(os.path.dirname(HERE), "scripts", "collect_conditions.py")
sys.path.insert(0, os.path.dirname(SCRIPT))
import collect_conditions as cc  # noqa: E402

RUNS = [f"20260101_TEST_60spd_S{i:02d}_A{i}_1_{100 + i}" for i in range(1, 11)]
LIMS_HEADER = ["internal_id", "internal_notes", "sample_name", "unique_id", "condition_name",
               "amt_2_inject"]
NS = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
RNS = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'


def col(i):
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def make_xlsx(path, sheets, inline=False, rich=None):
    """A minimal .xlsx: sheets = [(name, rows, hidden)], a row a list of cells (None = no cell).
    Text goes to sharedStrings (or inline); numbers are numeric cells. `rich` = {text: [runs]}
    stores that text as a rich-text shared string, with a phonetic guide that must not be read."""
    shared, index = [], {}

    def sid(text):
        if text not in index:
            index[text] = len(shared)
            shared.append(text)
        return index[text]
    parts = {}
    for k, (name, rows, _hidden) in enumerate(sheets, 1):
        xr = []
        for ri, row in enumerate(rows, 1):
            cells = []
            for ci, v in enumerate(row):
                ref = f"{col(ci)}{ri}"
                if v is None:
                    continue
                if isinstance(v, (int, float)):
                    cells.append(f'<c r="{ref}"><v>{v}</v></c>')
                elif inline:
                    cells.append(f'<c r="{ref}" t="inlineStr"><is><t>{escape(v)}</t></is></c>')
                else:
                    cells.append(f'<c r="{ref}" t="s"><v>{sid(v)}</v></c>')
            xr.append(f'<row r="{ri}">{"".join(cells)}</row>')
        parts[f"xl/worksheets/sheet{k}.xml"] = (
            f'<?xml version="1.0" encoding="UTF-8"?><worksheet {NS}><sheetData>'
            + "".join(xr) + "</sheetData></worksheet>")
    sis = []
    for t in shared:
        if rich and t in rich:
            sis.append("<si>" + "".join(f"<r><t>{escape(x)}</t></r>" for x in rich[t])
                       + "<rPh><t>PHONETIC</t></rPh></si>")
        else:
            sis.append(f"<si><t>{escape(t)}</t></si>")
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("[Content_Types].xml", '<?xml version="1.0"?><Types xmlns="http://schemas.'
                   'openxmlformats.org/package/2006/content-types"/>')
        z.writestr("xl/workbook.xml",
                   f'<?xml version="1.0"?><workbook {NS} {RNS}><sheets>'
                   + "".join(f'<sheet name="{escape(n)}" sheetId="{k}" r:id="rId{k}"'
                             + (' state="hidden"' if h else "") + "/>"
                             for k, (n, _r, h) in enumerate(sheets, 1))
                   + "</sheets></workbook>")
        z.writestr("xl/_rels/workbook.xml.rels",
                   '<?xml version="1.0"?><Relationships xmlns="http://schemas.openxmlformats.org/'
                   'package/2006/relationships">'
                   + "".join(f'<Relationship Id="rId{k}" Type="worksheet" '
                             f'Target="worksheets/sheet{k}.xml"/>' for k in range(1, len(sheets) + 1))
                   + "</Relationships>")
        for name, text in parts.items():
            z.writestr(name, text)
        if shared:
            z.writestr("xl/sharedStrings.xml", f'<?xml version="1.0"?><sst {NS}>'
                       + "".join(sis) + "</sst>")


def lims_rows(labels=None):
    labels = labels or ([f"CtrlA_p1_mix_{k}" for k in range(1, 6)]
                        + [f"TrtB_p1_mix_{k}" for k in range(1, 6)])
    return [LIMS_HEADER] + [["PROT_0000", "" if i % 3 else "re-injected", f"client{i + 1}",
                             f"S{i + 1:02d}", lab, 200] for i, lab in enumerate(labels)]


def run_map(tmp, sheet_path=None, extra=(), mapping=None, ok=True):
    out = os.path.join(tmp, "conditions.csv")
    args = [sys.executable, SCRIPT, "--map", out, "--runs", ",".join(RUNS), *extra]
    args += (["--mapping-json", json.dumps(mapping)] if mapping is not None
             else ["--from-file", sheet_path])
    p = subprocess.run(args, capture_output=True, text=True, timeout=60)
    if not ok:
        return p
    assert p.returncode == 0, p.stderr
    with open(out, newline="") as fh:
        return json.loads(p.stdout), list(csv.DictReader(fh))


class XlsxReader(unittest.TestCase):
    def test_strings_numbers_and_sparse_cells(self):
        with tempfile.TemporaryDirectory() as d:
            x = os.path.join(d, "s.xlsx")
            make_xlsx(x, [("samples", [["sample", "group", "amount", "note"],
                                        ["S01", "CtrlA_rich", 200, None],
                                        [None, None, None, None],          # an empty row
                                        ["S02", "B", 1.5, "x"]], False)],
                      rich={"CtrlA_rich": ["CtrlA", "_rich"]})
            headers, rows, info = cc.read_xlsx(x)
        self.assertEqual(headers, ["sample", "group", "amount", "note"])
        self.assertEqual(rows, [{"sample": "S01", "group": "CtrlA_rich", "amount": "200",
                                 "note": ""},
                                {"sample": "S02", "group": "B", "amount": "1.5", "note": "x"}])
        self.assertEqual((info["format"], info["sheet"]), ("xlsx", "samples"))

    def test_inline_strings_hidden_sheets_and_choosing_one(self):
        with tempfile.TemporaryDirectory() as d:
            x = os.path.join(d, "s.xlsx")
            make_xlsx(x, [("lookup", [["k"], ["v"]], True),
                          ("samples", [["sample", "group"], ["S01", "A"]], False),
                          ("notes", [["note"], ["n"]], False)], inline=True)
            _, rows, info = cc.read_xlsx(x)
            self.assertEqual((info["sheet"], rows), ("samples", [{"sample": "S01", "group": "A"}]))
            self.assertEqual(info["other_sheets"], ["lookup", "notes"])
            self.assertEqual(cc.read_xlsx(x, "notes")[1], [{"note": "n"}])
            with self.assertRaises(ValueError):
                cc.read_xlsx(x, "nope")

    def test_merged_cells_read_as_excel_shows_them_and_a_bad_index_is_an_error(self):
        with tempfile.TemporaryDirectory() as d:
            x = os.path.join(d, "m.xlsx")
            rows = [["sample", "group"], ["S01", "Ctrl"], ["S02", None], ["S03", "Trt"],
                    ["S04", None]]
            make_xlsx(x, [("samples", rows, False)])
            with zipfile.ZipFile(x) as z:
                parts = {n: z.read(n) for n in z.namelist()}
            sheet = parts["xl/worksheets/sheet1.xml"]
            parts["xl/worksheets/sheet1.xml"] = sheet.replace(
                b"</sheetData></worksheet>", b'</sheetData><mergeCells count="2"><mergeCell '
                b'ref="B2:B3"/><mergeCell ref="B4:B5"/></mergeCells></worksheet>')
            with zipfile.ZipFile(x, "w") as z:
                for n, v in parts.items():
                    z.writestr(n, v)
            _, got, info = cc.read_xlsx(x)
            self.assertEqual([r["group"] for r in got], ["Ctrl", "Ctrl", "Trt", "Trt"])
            self.assertEqual(info["merged_ranges_filled"], 2)
            parts["xl/worksheets/sheet1.xml"] = sheet.replace(b'<v>3</v>', b'<v>99</v>', 1)
            with zipfile.ZipFile(x, "w") as z:
                for n, v in parts.items():
                    z.writestr(n, v)
            with self.assertRaises(ValueError):
                cc.read_xlsx(x)

    def test_old_xls_and_a_non_zip_are_refused_with_what_to_do(self):
        with tempfile.TemporaryDirectory() as d:
            for name in ("old.xls", "fake.xlsx"):
                path = os.path.join(d, name)
                with open(path, "wb") as fh:
                    fh.write(b"\xd0\xcf\x11\xe0 not a zip")
                p = run_map(d, path, ok=False)
                self.assertNotEqual(p.returncode, 0)
                self.assertIn("save it from Excel as .xlsx", p.stderr + p.stdout)


class LimsSheet(unittest.TestCase):
    def test_the_lims_export_maps_and_its_replicates_collapse_once_confirmed(self):
        with tempfile.TemporaryDirectory() as d:
            x = os.path.join(d, "PROT_0000.samples.2026_01_01__00_00.xlsx")
            make_xlsx(x, [("samples", lims_rows(), False), ("notes", [["n"]], False)])
            asked, as_given = run_map(d, x)
            res, rows = run_map(d, x, extra=["--replicate-labels", "collapse"])
        # auto writes the reading into the PROPOSED csv and asks the user to confirm it
        self.assertTrue(asked["needs_confirmation"])
        prop = asked["ambiguities"]["replicate_labels_to_confirm"]
        self.assertTrue(prop["collapsed"])
        self.assertEqual(prop["reasons"], [])
        self.assertIn("ASK the user", prop["to_confirm"])
        self.assertEqual(prop["conditions"]["TrtB_p1_mix"]["3"], "TrtB_p1_mix_3")
        self.assertEqual(asked["groups"], {"CtrlA_p1_mix": 5, "TrtB_p1_mix": 5})
        self.assertEqual(list(as_given[0]), ["File.Name", "Group", "Label"])
        self.assertEqual(res["groups"], {"CtrlA_p1_mix": 5, "TrtB_p1_mix": 5})
        self.assertFalse(res["needs_confirmation"], json.dumps(res["ambiguities"], indent=1))
        self.assertEqual(res["sample_column"], "unique_id")
        self.assertEqual(res["sample_column_matches"], {"sample_name": 0, "unique_id": 10})
        self.assertEqual(set(res["columns_skipped"]),
                         {"internal_id", "internal_notes", "sample_name", "amt_2_inject"})
        self.assertEqual(res["covariate_columns"], {})
        self.assertEqual(res["source"]["sheet"], "samples")
        self.assertTrue(res["replicate_labels"]["collapsed"])
        self.assertEqual(res["replicate_labels"]["conditions"]["TrtB_p1_mix"]["3"],
                         "TrtB_p1_mix_3")
        self.assertEqual(list(rows[0]), ["File.Name", "Group", "Label"])
        self.assertEqual((rows[2]["Group"], rows[2]["Label"]), ("CtrlA_p1_mix", "CtrlA_p1_mix_3"))
        self.assertEqual(res["replicates"][RUNS[7]], 3)

    def test_a_gap_is_asked_not_guessed(self):
        labels = ([f"CtrlA_mix_{k}" for k in (1, 2, 4, 5, 6)]
                  + [f"TrtB_mix_{k}" for k in range(1, 6)])
        with tempfile.TemporaryDirectory() as d:
            x = os.path.join(d, "s.xlsx")
            make_xlsx(x, [("samples", lims_rows(labels), False)])
            res, rows = run_map(d, x)
            self.assertTrue(res["needs_confirmation"])
            amb = res["ambiguities"]["replicate_labels_ambiguous"]
            self.assertFalse(amb["collapsed"])
            self.assertTrue(any("CtrlA_mix: replicate numbers 1, 2, 4, 5, 6" in r
                                for r in amb["reasons"]), amb["reasons"])
            self.assertIn("ASK the user", amb["to_confirm"])
            self.assertEqual(len(res["groups"]), 10)                 # as given: singletons
            self.assertNotIn("Label", rows[0])
            res, rows = run_map(d, x, extra=["--replicate-labels", "collapse"])
            self.assertEqual(res["groups"], {"CtrlA_mix": 5, "TrtB_mix": 5})
            self.assertNotIn("replicate_labels_ambiguous", res["ambiguities"])
            res, _ = run_map(d, x, extra=["--replicate-labels", "keep"])
            self.assertEqual(len(res["groups"]), 10)
            self.assertNotIn("replicate_labels_ambiguous", res["ambiguities"])

    def test_numbers_that_are_the_conditions_are_asked(self):
        for labels, why in (([f"Day{k}" for k in range(1, 11)], "make 1 condition"),
                            ([f"A_{k}" for k in range(1, 6)] + [f"B{x}" for x in "VWXYZ"],
                             "end in no replicate number"),
                            ([f"A_{k}" for k in range(1, 6)] + [f"a-{k}" for k in range(1, 6)],
                             "differ only in case or punctuation")):
            rec = cc.collapse_replicates(
                [{"sample": f"S{i:02d}", "group": g, "extras": {}} for i, g in
                 enumerate(labels, 1)])[1]
            self.assertFalse(rec["collapsed"], labels)
            self.assertTrue(any(why in r for r in rec["reasons"]), (why, rec["reasons"]))

    def test_forcing_labels_without_numbers_collapses_nothing(self):
        intent = [{"sample": f"S{i:02d}", "group": g, "extras": {}}
                  for i, g in enumerate(["A", "A", "B", "B"], 1)]
        out, rec = cc.collapse_replicates(intent, "collapse")
        self.assertFalse(rec["collapsed"])
        self.assertEqual([it["group"] for it in out], ["A", "A", "B", "B"])

    def test_shared_labels_are_left_alone(self):
        intent = [{"sample": f"S{i:02d}", "group": "A" if i < 6 else "B", "extras": {}}
                  for i in range(1, 11)]
        out, rec = cc.collapse_replicates(intent)
        self.assertIsNone(rec)
        self.assertEqual(out, intent)

    def test_the_mapping_json_route_collapses_too(self):
        """core_submission.py conditions hands CoreOmics condition_name to --map this way."""
        mapping = {"mapping": {r: f"{'CtrlA' if i < 5 else 'TrtB'}_{i % 5 + 1}"
                               for i, r in enumerate(RUNS)}}
        with tempfile.TemporaryDirectory() as d:
            asked, _ = run_map(d, mapping=mapping)
            res, rows = run_map(d, mapping=mapping, extra=["--replicate-labels", "collapse"])
        self.assertTrue(asked["needs_confirmation"])
        self.assertIn("replicate_labels_to_confirm", asked["ambiguities"])
        self.assertEqual(res["groups"], {"CtrlA": 5, "TrtB": 5})
        self.assertEqual(rows[0]["Label"], "CtrlA_1")


class WhichRunIsWhichSample(unittest.TestCase):
    """Matching by whole words, never a bare substring and never the plate position (review of
    2.10): a submitter's tubes labelled A1..A4 matched the runs in the Core's wells A1..A4
    (`_S3-A1_`), and a tie between two identifier columns was settled by the sheet's order."""

    def test_words_not_substrings_and_never_a_plate_position(self):
        runs = ["03122025__60SPD_DIA-KG1_S3-A3_1_24003", "03122025__60SPD_DIA-KG12_S3-A1_1_24001",
                "Ctrl1", "Treated_rep1", "run1"]
        m = cc.match_to_runs
        self.assertEqual(m("A1", runs), [])                       # the well, never the sample
        self.assertEqual(m("S3-A3", runs), [])
        self.assertEqual(m("KG1", runs), [runs[0]])               # not KG12
        self.assertEqual(m("Ctrl", runs), ["Ctrl1"])              # a letters-only prefix: kept
        self.assertEqual(m("Treated", runs), ["Treated_rep1"])
        self.assertEqual(m("/data/x/run1.raw", runs), ["run1"])   # a full file name
        self.assertEqual(m("dia-kg1", runs), [runs[0]])

    def test_a_letters_only_label_naming_several_runs_is_asked_not_assigned(self):
        """'Ctrl' names Ctrl1, Ctrl2 and Ctrl3 (a letters-only word names itself plus digits):
        that is a question -- every one of them, or one sample? -- never the first run, and no
        run gets the group until the user answers (--confirm-multi)."""
        runs = ["Ctrl1", "Ctrl2", "Ctrl3", "Trt1", "Trt2", "Trt3"]
        mapping = {"mapping": {"Ctrl": "control", "Trt1": "treated", "Trt2": "treated",
                               "Trt3": "treated"}}
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "conditions.csv")
            base = [sys.executable, SCRIPT, "--map", out, "--runs", ",".join(runs),
                    "--mapping-json", json.dumps(mapping)]
            res = json.loads(subprocess.run(base, capture_output=True, text=True,
                                            timeout=60).stdout)
            with open(out, newline="") as fh:
                group = {r["File.Name"]: r["Group"] for r in csv.DictReader(fh)}
            yes = json.loads(subprocess.run(base + ["--confirm-multi", "Ctrl"],
                                            capture_output=True, text=True, timeout=60).stdout)
        self.assertTrue(res["needs_confirmation"])
        self.assertEqual(res["ambiguities"]["multi_match_identifiers"]["Ctrl"],
                         ["Ctrl1", "Ctrl2", "Ctrl3"])
        self.assertEqual([q["kind"] for q in res["decisions_required"]],
                         ["label_names_several_runs"])
        self.assertEqual(sorted(res["ambiguities"]["awaiting_decision_runs"]),
                         ["Ctrl1", "Ctrl2", "Ctrl3"])
        self.assertEqual([group[r] for r in ("Ctrl1", "Ctrl2", "Ctrl3")], ["", "", ""])
        self.assertFalse(yes["needs_confirmation"], yes["ambiguities"])
        self.assertEqual(yes["groups"], {"control": 3, "treated": 3})

    def test_a_tied_identifier_column_is_asked(self):
        well = {"KG1": "A3", "KG2": "A1", "KG3": "A4", "KG4": "A2"}   # the Core's plate
        runs = [f"03122025__60SPD_DIA-{u}_{w}x_S3-{w}_1_2400{i}" for i, (u, w) in
                enumerate(well.items())]                               # both columns name them
        rows = [LIMS_HEADER] + [["PROT_0000", "", f"{well[f'KG{i}']}x", f"KG{i}",
                                 "ctrl" if i < 3 else "treat", 200] for i in range(1, 5)]
        with tempfile.TemporaryDirectory() as d:
            x = os.path.join(d, "s.xlsx")
            make_xlsx(x, [("samples", rows, False)])
            out = os.path.join(d, "conditions.csv")
            base = [sys.executable, SCRIPT, "--map", out, "--runs", ",".join(runs),
                    "--from-file", x]
            p = subprocess.run(base, capture_output=True, text=True, timeout=60)
            res = json.loads(p.stdout)
            with open(out, newline="") as fh:
                groups = {r["Group"] for r in csv.DictReader(fh)}
            p2 = subprocess.run(base + ["--sample-column", "unique_id"], capture_output=True,
                                text=True, timeout=60)
            res2 = json.loads(p2.stdout)
        self.assertTrue(res["needs_confirmation"])
        tie = res["ambiguities"]["sample_column_tie"]
        self.assertEqual(tie["columns"], {"sample_name": 4, "unique_id": 4})
        self.assertEqual(groups, {""})                              # no guess written
        self.assertFalse(res2["needs_confirmation"], res2["ambiguities"])
        self.assertEqual(res2["groups"], {"ctrl": 2, "treat": 2})

    def test_numbered_labels_beside_a_numbered_control_are_asked(self):
        """Day1..3 beside Ctl1..3 passes every 'unambiguous' test -- and is a time course."""
        labels = ["Day1", "Day2", "Day3", "Ctl1", "Ctl2", "Ctl3"]
        runs = [f"run{i:02d}" for i in range(len(labels))]
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "conditions.csv")
            p = subprocess.run([sys.executable, SCRIPT, "--map", out, "--runs", ",".join(runs),
                                "--mapping-json", json.dumps({"mapping": dict(zip(runs, labels))})],
                               capture_output=True, text=True, timeout=60)
        res = json.loads(p.stdout)
        self.assertTrue(res["needs_confirmation"])
        self.assertEqual(sorted(res["ambiguities"]["replicate_labels_to_confirm"]["conditions"]),
                         ["Ctl", "Day"])
        self.assertIn("time points", res["ambiguities"]["replicate_labels_to_confirm"]["to_confirm"])


if __name__ == "__main__":
    unittest.main()
