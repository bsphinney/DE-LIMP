"""
The CoreOmics submission workflow: every check here guards a failure that does NOT raise.

A service run attaches a PI's name, a sample sheet and a Bioshare folder to a set of raw
files. When any link in that chain is wrong the pipeline still succeeds -- it searches the
next lab's BN1, files the results under the wrong month, or hands the collaborator a share
full of symlinks they cannot open (PROT_0793). So the assertions are on the exit codes the
orchestrator branches on and on the filesystem state a collaborator will actually see.

Stdlib only, no network, no SSH: temp dirs, env overrides, and an in-process fake CoreOmics.
"""
import contextlib
import csv
import datetime as dt
import http.server
import io
import json
import ntpath
import os
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import unittest
import urllib.parse
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)

import core_submission as cs  # noqa: E402

CS = os.path.join(SCRIPTS, "core_submission.py")
HEX = "99922f5337f8"
TOKEN = "test-token-not-real"
SERVER_SHARE = f"/nfs/lssc0/flinders/proteomics/coreomics/projects/2025/03/{HEX}/share"
IS_ROOT = hasattr(os, "geteuid") and os.geteuid() == 0


# ------------------------------------------------------------------------ fixtures --
def record(internal_id="PROT_0807", sid=HEX, submitted="2025-03-10T15:02:53.863095-07:00",
           samples=None, institution="UC Davis", pi=True, raw_only=False, instrument="timsTOF HT",
           pi_last="Placeholder", pi_first="Quinn"):
    if samples is None:
        samples = [("KG1", "ctrl"), ("KG2", "ctrl"), ("KG13", "treat"), ("KG14", "treat")]
    return {
        "id": sid, "internal_id": internal_id, "url": f"https://coreomics.example/submissions/{sid}",
        "status": "Samples Received", "submitted": submitted,
        "first_name": "Ada", "last_name": "Example", "email": "ada.example@example.org",
        "pi_first_name": pi_first, "pi_last_name": pi_last, "pi_email": "pi.placeholder@example.org",
        "pi": {"institution": {"name": institution}, "department": "Example Biology"} if pi else None,
        "institute": "Example Institute",
        "contacts": [{"first_name": "Robin", "last_name": "Contact", "email": "robin.contact@example.org"}],
        "submission_data": {
            "organism": "mouse", "description": "Example description.", "sample_prep": "In-solution digest",
            "data_analysis": ("I only require raw data and will do my own data analysis" if raw_only
                              else "I want the proteomics core to do the data analysis"),
            "proteomics_type": ["Global proteomics"], "mass_spec_wanted": instrument,
            "samples": [{"unique_id": u, "sample_name": f"sample {u}", "condition_name": c}
                        for u, c in samples]},
    }


def roots(tmp):
    return {"CORE_FLINDERS_ROOT": os.path.join(tmp, "flinders"), "CORE_WORK_ROOT": os.path.join(tmp, "work")}


def child_env(tmp, base_url="http://127.0.0.1:9/server/api", **extra):
    e = dict(os.environ)
    for k in ("http_proxy", "https_proxy", "HTTP_PROXY", "HTTPS_PROXY", "ALL_PROXY", "all_proxy"):
        e.pop(k, None)
    e["NO_PROXY"] = e["no_proxy"] = "127.0.0.1,localhost"
    e.update(roots(tmp))
    e.update(COREOMICS_TOKEN=TOKEN, COREOMICS_BASE_URL=base_url)
    e.update(extra)
    return e


def read(path):
    with open(path) as fh:
        return fh.read()


def load(path):
    with open(path) as fh:
        return json.load(fh)


def tsv_rows(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run(args, env):
    p = subprocess.run([sys.executable, CS] + args, capture_output=True, text=True, env=env, timeout=120)
    try:
        out = json.loads(p.stdout)
    except ValueError:
        out = None
    return p.returncode, out, p


def write_summary(tmp, rec=None, neighbors=None, name="submission_summary.json"):
    with mock.patch.dict(os.environ, roots(tmp)):
        s = cs.build_summary(rec or record(), neighbors or [], [], None, 240, False, [])
    path = os.path.join(tmp, name)
    with open(path, "w") as fh:
        json.dump(s, fh)
    return path, s


def local_share(tmp, year="2025", month="03", sid=HEX):
    return os.path.join(tmp, "flinders", "coreomics", "projects", year, month, sid, "share")


def raw(tmp, inst, month, name, mtime=None):
    p = os.path.join(tmp, "flinders", "Data", "raw_data", inst, month, name)
    if name.endswith(".d"):
        os.makedirs(p)
    else:
        os.makedirs(os.path.dirname(p), exist_ok=True)
        open(p, "w").close()
    if mtime:
        ts = dt.datetime(*mtime).timestamp()
        os.utime(p, (ts, ts))
    return p


def tims(tmp, mmddyyyy, uid, n, month="mar25"):
    return raw(tmp, "tTOF_HT", month, f"{mmddyyyy}__60SPD_DIA-{uid}_S3-A{n}_1_{24000 + n}.d")


def gates(out):
    return {g["gate"]: g for g in out["gates"]}


def write(root, rel, text):
    p = os.path.join(root, rel)
    os.makedirs(os.path.dirname(p), exist_ok=True)
    with open(p, "w") as fh:
        fh.write(text)
    return p


class FakeCoreOmics:
    """A tiny CoreOmics: routes are (method, path) -> fn(query, body) -> (status, json[, headers])."""

    def __init__(self):
        self.requests, self.routes = [], {}
        fake = self

        class H(http.server.BaseHTTPRequestHandler):
            def log_message(self, *a):
                pass

            def _do(self, method):
                u = urllib.parse.urlsplit(self.path)
                q = dict(urllib.parse.parse_qsl(u.query))
                n = int(self.headers.get("Content-Length") or 0)
                body = json.loads(self.rfile.read(n)) if n else None
                fake.requests.append({"method": method, "path": u.path, "query": q, "body": body,
                                      "auth": self.headers.get("Authorization")})
                fn = fake.routes.get((method, u.path))
                res = fn(q, body) if fn else (404, {"detail": "Not found."})
                status, obj = res[0], res[1]
                data = json.dumps(obj).encode()
                self.send_response(status)
                for k, v in (res[2] if len(res) > 2 else {}).items():
                    self.send_header(k, v)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(data)))
                self.end_headers()
                self.wfile.write(data)

            def do_GET(self):
                self._do("GET")

            def do_POST(self):
                self._do("POST")

        self.httpd = http.server.ThreadingHTTPServer(("127.0.0.1", 0), H)
        self.thread = threading.Thread(target=self.httpd.serve_forever, daemon=True)

    @property
    def base(self):
        return f"http://127.0.0.1:{self.httpd.server_address[1]}/server/api"

    def __enter__(self):
        self.thread.start()
        return self

    def __exit__(self, *exc):
        self.httpd.shutdown()
        self.httpd.server_close()

    def posts(self):
        return [r for r in self.requests if r["method"] == "POST"]


SHARES = f"/server/api/plugins/bioshare/submissions/{HEX}/submission_shares/"


# ------------------------------------------------------------------ pure functions --
class TestSubmissionIds(unittest.TestCase):
    def test_every_accepted_spelling_names_the_same_submission(self):
        for s in ("807", "0807", "PROT_0807", "prot-807", "#807", " PROT 807 ", "prot_0807"):
            self.assertEqual(cs.normalize_submission(s), ("internal_id", "PROT_0807"), s)

    def test_a_12_hex_id_is_a_detail_lookup(self):
        self.assertEqual(cs.normalize_submission(HEX), ("id", HEX))
        self.assertEqual(cs.normalize_submission(HEX.upper()), ("id", HEX))

    def test_garbage_is_refused_not_guessed(self):
        for s in ("", "abc", "PROT_", "0", "12345678", HEX[:-1], "807a"):
            with self.assertRaises(ValueError, msg=s):
                cs.normalize_submission(s)


class TestCanonicalProjectDir(unittest.TestCase):
    def test_late_evening_submission_stays_in_its_literal_month(self):
        """coreomics_fs takes the literal YYYY-MM; converting to UTC would say 2026/09."""
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(os.environ, roots(tmp)):
            late = "2026-08-31T23:30:00-07:00"
            self.assertEqual(cs.server_project_dir(HEX, late),
                             f"/nfs/lssc0/flinders/proteomics/coreomics/projects/2026/08/{HEX}")
            self.assertEqual(cs.local_project_dir(HEX, late),
                             os.path.join(tmp, "flinders", "coreomics", "projects", "2026", "08", HEX))
            self.assertEqual(cs.parse_date(late), dt.date(2026, 8, 31))

    def test_share_dir_is_the_server_path_whatever_the_local_root(self):
        """Defect 10: Bioshare only knows HIVE's path. A local root (an SMB mount, a test dir)
        must never leak into the summary's share_dir."""
        with tempfile.TemporaryDirectory() as tmp:
            _, s = write_summary(tmp, record(submitted="2026-09-10T15:02:53.863095-07:00"))
            self.assertEqual(s["share_dir"],
                             f"/nfs/lssc0/flinders/proteomics/coreomics/projects/2026/09/{HEX}/share")
            self.assertNotIn(tmp, s["share_dir"])

    def test_server_paths_use_forward_slashes_under_windows_path_rules(self):
        """Defect 10: on native Windows os.path.join would have produced `\\nfs\\...`."""
        with mock.patch.object(cs.os, "path", ntpath):
            d = cs.server_project_dir(HEX, "2025-03-10T15:02:53-07:00")
        self.assertEqual(d, f"/nfs/lssc0/flinders/proteomics/coreomics/projects/2025/03/{HEX}")

    def test_unusable_submitted_is_an_error(self):
        with self.assertRaises(ValueError):
            cs.server_project_dir(HEX, "September 10")


class TestCampus(unittest.TestCase):
    def test_uc_davis_exactly_is_on_campus(self):
        self.assertEqual(cs.classify_campus(record())[0], "on_campus")
        self.assertEqual(cs.classify_campus(record(institution="University of Washington"))[0], "off_campus")

    def test_null_pi_falls_back_to_institute(self):
        r = record(pi=False)
        r["institute"] = "UC Davis"
        self.assertEqual(cs.classify_campus(r), ("on_campus", "UC Davis", "institute"))
        r["institute"] = "Example Corp"
        self.assertEqual(cs.classify_campus(r)[0], "off_campus")
        r["institute"] = None
        self.assertEqual(cs.classify_campus(r)[2], "missing")


class TestFilenameDates(unittest.TestCase):
    TODAY = dt.date(2026, 9, 16)

    def test_each_naming_convention(self):
        cases = {
            "09012026__60SPD_DIA-KG1_S3-A5_1_24026.d": dt.date(2026, 9, 1),     # tims MMDDYYYY
            "20260827_793_100spd_Hel50_S6-A12_1_24026.d": dt.date(2026, 8, 27),  # HT YYYYMMDD
            "Ex08312026_380_JE21.raw": dt.date(2026, 8, 31),                     # Exploris MMDDYYYY
            "Ex040826_HeL50_1.raw": dt.date(2026, 8, 4),                         # Exploris DDMMYY
            "Ex070225_x.raw": dt.date(2025, 2, 7),
            "FL280826_HeL50.raw": dt.date(2026, 8, 28),                          # Lumos DDMMYY
        }
        for name, want in cases.items():
            self.assertEqual(cs.date_from_name(name, self.TODAY), want, name)

    def test_typos_undated_and_future_names_have_no_filename_date(self):
        for name in ("070622026__60SPD_DIA-X.d", "062420266__60SPD.d", "FLsep26_wa_20260909005054.raw",
                     "FL9sep_waID_3.raw", "12312099__future.d", "13452026__bad.d"):
            self.assertIsNone(cs.date_from_name(name, self.TODAY), name)

    def test_undated_file_falls_back_to_mtime_and_says_so(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = raw(tmp, "Lumos1", "sep26", "FLsep26_wa_20260909005054.raw", mtime=(2025, 4, 2, 12, 0))
            self.assertEqual(cs.acquisition_date(p), (dt.date(2025, 4, 2), "mtime"))


class TestTokens(unittest.TestCase):
    def test_delimited_so_kg1_is_not_kg13(self):
        pat = cs.token_pattern("KG1")
        self.assertTrue(pat.search("09012026__60SPD_DIA-KG1_S3-A5_1_24026.d"))
        self.assertFalse(pat.search("09012026__60SPD_DIA-KG13_S3-A5_1_24026.d"))

    def test_separators_are_equivalent_and_case_is_ignored(self):
        pat = cs.token_pattern("EB_001")
        for name in ("x_EB-001_y.d", "x_EB_001.raw", "x EB 001.d", "x_eb-001_y.d"):
            self.assertTrue(pat.search(name), name)
        self.assertTrue(cs.token_pattern("kg1").search("DIA-KG1_S3.d"))

    def test_weak_ids(self):
        for uid in ("A5", "H10", "B2", "001", "7", ""):
            self.assertIsNotNone(cs.weak_reason(uid), uid)
        for uid in ("KG1", "EB_001", "JE21", "I13"):
            self.assertIsNone(cs.weak_reason(uid), uid)

    def test_timstof_names_expose_only_the_sample_field(self):
        """Defect 3: the `_S3-A1_1_` plate position must not be matchable."""
        self.assertEqual(cs.match_space("09012026__60SPD_DIA-KG1_S3-A1_1_24026.d"), "KG1")
        self.assertEqual(cs.match_space("070622026__60SPD_DIA-EB_001_S3-B2_1_24026.d"), "EB_001")
        self.assertEqual(cs.match_space("Ex08312026_380_JE21.raw"), "Ex08312026_380_JE21.raw")


class TestShareName(unittest.TestCase):
    def test_sanitized_to_the_serializer_regex(self):
        s = {"id": HEX, "internal_id": "PROT_0807",
             "pi": {"last": "O'Neill (Lab)/Ext", "first": "Quinn"}}
        name = cs.share_name(s)
        self.assertRegex(name, r"^[\w\d\s'\"\.!\?\-:,]+$")
        self.assertEqual(name, "O'Neill Lab Ext, Quinn: PROT_0807")

    def test_default_shape(self):
        self.assertEqual(cs.share_name({"id": HEX, "internal_id": "PROT_0807",
                                        "pi": {"last": "Placeholder", "first": "Quinn"}}),
                         "Placeholder, Quinn: PROT_0807")


# ------------------------------------------------------------------------ locate --
class TestLocate(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.env = child_env(self.tmp)
        self.out = os.path.join(self.tmp, "loc")

    def tearDown(self):
        self._tmp.cleanup()

    def clean_tree(self):
        return {u: tims(self.tmp, "03122025", u, i) for i, u in enumerate(("KG1", "KG2", "KG13", "KG14"), 1)}

    def locate(self, summary, *extra):
        return run(["locate", "--summary", summary, "--out", self.out, *extra], self.env)

    def files_txt(self):
        return read(os.path.join(self.out, "files.txt"))

    def test_clean_submission_passes_with_one_file_per_sample(self):
        files = self.clean_tree()
        summary, _ = write_summary(self.tmp)
        rc, out, p = self.locate(summary)
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual(sorted(self.files_txt().split()), sorted(files.values()))
        rows = {r["unique_id"]: r for r in tsv_rows(os.path.join(self.out, "sample_files.tsv"))}
        self.assertEqual(rows["KG1"]["file"], files["KG1"], "KG1 must not pick up KG13's file")
        self.assertTrue(all(r["status"] == "matched" for r in rows.values()))
        self.assertEqual(rows["KG1"]["acquired"], "2025-03-12")

    def test_runs_before_the_submission_are_not_its_files(self):
        self.clean_tree()
        old = tims(self.tmp, "02012025", "KG1", 9, month="feb25")
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 0)
        self.assertNotIn(old, self.files_txt())
        self.assertEqual(gates(out)["out_of_window"]["status"], "INFO")

    def test_a_sample_with_only_old_runs_is_unmatched_and_blocks(self):
        for i, u in enumerate(("KG2", "KG13", "KG14"), 2):
            tims(self.tmp, "03122025", u, i)
        tims(self.tmp, "02012025", "KG1", 1, month="feb25")
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        self.assertEqual(gates(out)["unmatched_samples"]["status"], "FAIL")
        self.assertEqual(len(self.files_txt().split()), 3)     # proposal files are still written
        rc, out, _ = self.locate(summary, "--allow-partial")
        self.assertEqual(rc, 0)
        self.assertEqual(gates(out)["unmatched_samples"]["status"], "WARN")

    def test_a_newer_submission_using_the_label_makes_the_run_ambiguous(self):
        for i, u in enumerate(("KG2", "KG13", "KG14"), 2):
            tims(self.tmp, "03122025", u, i)
        tims(self.tmp, "03252025", "KG1", 1)
        neighbor = {"internal_id": "PROT_0810", "id": "aaaaaaaaaaaa",
                    "submitted": "2025-03-20T09:00:00-07:00", "unique_ids": ["KG1", "KG9"]}
        summary, _ = write_summary(self.tmp, neighbors=[neighbor])
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        g = gates(out)["ambiguous_label"]
        self.assertEqual(g["status"], "FAIL")
        self.assertEqual(g["samples"][0]["files"][0]["contenders"][0]["submission"], "PROT_0810")

    def test_a_same_day_neighbour_makes_the_run_ambiguous(self):
        """Defect 2: `lo < nd` missed a neighbour submitted the same day."""
        self.clean_tree()
        neighbor = {"internal_id": "PROT_0808", "id": "aaaaaaaaaaaa",
                    "submitted": "2025-03-10T09:00:00-07:00", "unique_ids": ["KG2"]}
        summary, _ = write_summary(self.tmp, neighbors=[neighbor])
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        self.assertEqual(gates(out)["ambiguous_label"]["samples"][0]["unique_id"], "KG2")

    def test_an_older_neighbour_is_named_with_its_earlier_unambiguous_runs(self):
        """Defect 2, the PROT_0794/0776 shape: an OLDER submission used the label first and
        already has its own runs. Labels alone cannot decide; the detail gives staff the facts."""
        self.clean_tree()
        tims(self.tmp, "02252025", "KG1", 8, month="feb25")          # the older lab's own run
        neighbor = {"internal_id": "PROT_0776", "id": "aaaaaaaaaaaa",
                    "submitted": "2025-02-20T09:00:00-08:00", "unique_ids": ["KG1"]}
        summary, _ = write_summary(self.tmp, neighbors=[neighbor])
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        contender = gates(out)["ambiguous_label"]["samples"][0]["files"][0]["contenders"][0]
        self.assertEqual(contender["submission"], "PROT_0776")
        self.assertTrue(contender["has_earlier_unambiguous_files"])
        self.assertEqual(contender["earlier_files"], 1)

    def test_accept_ambiguous_takes_the_runs_and_is_recorded(self):
        files = self.clean_tree()
        neighbor = {"internal_id": "PROT_0776", "id": "aaaaaaaaaaaa",
                    "submitted": "2025-02-20T09:00:00-08:00", "unique_ids": ["KG1"]}
        summary, _ = write_summary(self.tmp, neighbors=[neighbor])
        rc, out, _ = self.locate(summary, "--accept-ambiguous")
        self.assertEqual(rc, 0)
        self.assertEqual(gates(out)["ambiguous_label"]["status"], "WARN")
        self.assertIn(files["KG1"], self.files_txt())
        self.assertTrue(load(os.path.join(self.out, "locate.json"))["accepted"]["accept_ambiguous"])

    def test_an_unambiguous_run_wins_over_an_ambiguous_one(self):
        files = self.clean_tree()
        tims(self.tmp, "03252025", "KG1", 7)
        neighbor = {"internal_id": "PROT_0810", "id": "aaaaaaaaaaaa",
                    "submitted": "2025-03-20T09:00:00-07:00", "unique_ids": ["kg1"]}
        summary, _ = write_summary(self.tmp, neighbors=[neighbor])
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 0)
        self.assertIn(files["KG1"], self.files_txt())
        self.assertEqual(gates(out)["ambiguous_files_excluded"]["status"], "WARN")

    def test_max_days_wider_than_the_neighbour_window_is_refused(self):
        """Defect 2: beyond fetch's window nobody checked who else used the label."""
        self.clean_tree()
        summary, _ = write_summary(self.tmp)
        self.assertEqual(self.locate(summary, "--max-days", "300")[0], 2)

    def test_several_files_per_sample_most_recent_wins_with_alternates(self):
        files = self.clean_tree()
        newer = tims(self.tmp, "03142025", "KG2", 8)
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 0)
        g = gates(out)["alternates"]
        self.assertEqual(g["status"], "WARN")
        self.assertEqual(g["samples"][0]["chosen"], newer)
        self.assertEqual(g["samples"][0]["alternates"], [files["KG2"]])

    def test_ht_plate_token_exits_4(self):
        self.clean_tree()
        raw(self.tmp, "tTOF_HT", "mar25", "20250315_807_100spd_Hel50_S6-A1_1_24026.d")
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 4)
        self.assertTrue(out["ht_plate"])
        self.assertIn("ht_manifest.py", out["next_step"])

    def test_a_run_counter_that_equals_the_number_is_not_an_ht_plate(self):
        self.clean_tree()
        raw(self.tmp, "Exploris480", "mar25", "Ex03152025_807_JE21.raw")
        summary, _ = write_summary(self.tmp)
        self.assertEqual(self.locate(summary)[0], 0)

    def test_well_like_id_is_weak_and_blocks(self):
        rec = record(samples=[("A5", "ctrl"), ("KG2", "ctrl")])
        tims(self.tmp, "03122025", "KG2", 2)
        summary, _ = write_summary(self.tmp, rec)
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        self.assertEqual(gates(out)["weak_ids"]["status"], "FAIL")

    def test_a_well_plus_suffix_id_does_not_match_the_plate_position_field(self):
        """Defect 3: A1-1 used to match `_S3-A1_1_` in EVERY timsTOF name and exit 0."""
        foreign = [tims(self.tmp, "03122025", f"OTHERLAB{i}", 1) for i in range(4)]
        tims(self.tmp, "03122025", "KG2", 2)
        summary, _ = write_summary(self.tmp, record(samples=[("A1-1", "ctrl"), ("KG2", "ctrl")]))
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        for f in foreign:
            self.assertNotIn(f, self.files_txt())

    def test_an_id_naming_more_than_ten_runs_is_weak(self):
        """Defect 3: frequency across the WHOLE listing, regardless of date."""
        for i in range(11):
            raw(self.tmp, "Exploris480", "mixed", f"Ex0{(i % 9) + 1}122024_{300 + i}_LIVER.raw")
        raw(self.tmp, "Exploris480", "mar25", "Ex03122025_400_LIVER.raw")
        tims(self.tmp, "03122025", "KG2", 2)
        summary, _ = write_summary(self.tmp, record(samples=[("LIVER", "ctrl"), ("KG2", "ctrl")]))
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        self.assertIn("LIVER", [s["unique_id"] for s in gates(out)["weak_ids"]["samples"]])

    def test_the_longest_id_owns_a_run(self):
        """Defect 4: DH1 must not take DH1-1's run just because DH1 is a delimited token in it."""
        for i, u in enumerate(("KG1", "KG2", "KG13"), 1):
            tims(self.tmp, "03122025", u, i)
        other = tims(self.tmp, "03122025", "DH1-1", 5)
        neighbor = {"internal_id": "PROT_0790", "id": "aaaaaaaaaaaa",
                    "submitted": "2025-01-15T09:00:00-08:00", "unique_ids": ["DH1-1"]}
        rec = record(samples=[("KG1", "a"), ("KG2", "a"), ("KG13", "b"), ("DH1", "b")])
        summary, _ = write_summary(self.tmp, rec, neighbors=[neighbor])
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        self.assertNotIn(other, self.files_txt())
        rows = {r["unique_id"]: r for r in tsv_rows(os.path.join(self.out, "sample_files.tsv"))}
        self.assertEqual(rows["DH1"]["status"], "unmatched")
        self.assertEqual(gates(out)["shadowed_by_longer_id"]["status"], "INFO")

    def test_files_from_maps_back_and_validates_existence(self):
        files = self.clean_tree()
        cur = os.path.join(self.tmp, "curated.txt")
        with open(cur, "w") as fh:
            fh.write("\n".join(files.values()) + "\n")
        summary, _ = write_summary(self.tmp)
        rc, out, p = self.locate(summary, "--files-from", cur)
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual(out["counts"]["matched"], 4)
        with open(cur, "a") as fh:
            fh.write(os.path.join(self.tmp, "gone.d") + "\n")
        rc, out, _ = self.locate(summary, "--files-from", cur)
        self.assertEqual(rc, 2)
        self.assertEqual(gates(out)["paths_exist"]["status"], "FAIL")


# ------------------------------------------------------------------------- stage --
class TestStage(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.env = child_env(self.tmp)
        self.sroot = os.path.join(self.tmp, "flinders", "Data", "lab", "service")
        os.makedirs(os.path.join(self.sroot, "on_campus"))
        os.makedirs(os.path.join(self.sroot, "off_campus"))
        self.files = [tims(self.tmp, "03122025", u, i) for i, u in enumerate(("KG1", "KG2"), 1)]
        self.list = os.path.join(self.tmp, "files.txt")
        with open(self.list, "w") as fh:
            fh.write("\n".join(self.files) + "\n")

    def tearDown(self):
        self._tmp.cleanup()

    def stage(self, summary, *extra):
        return run(["stage", "--summary", summary, "--files", self.list, *extra], self.env)

    def folder(self, *parts):
        p = os.path.join(self.sroot, *parts)
        os.makedirs(p)
        return p

    def test_dry_run_creates_nothing(self):
        self.folder("on_campus", "Placeholder lab")
        summary, _ = write_summary(self.tmp)
        rc, out, p = self.stage(summary)
        self.assertEqual(rc, 0, p.stderr)
        self.assertFalse(out["applied"])
        self.assertEqual(out["would_create"], 2)
        self.assertFalse(os.path.exists(out["service_dir"]))
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "work")))

    def test_apply_creates_relative_links_that_resolve_and_is_idempotent(self):
        group = self.folder("on_campus", "Placeholder lab")
        summary, _ = write_summary(self.tmp)
        rc, out, p = self.stage(summary, "--apply")
        self.assertEqual(rc, 0, p.stderr)
        project = os.path.join(group, "PROT_0807")
        self.assertEqual(out["service_dir"], project)
        for f in self.files:
            link = os.path.join(project, "raw", os.path.basename(f))
            self.assertFalse(os.path.isabs(os.readlink(link)), "absolute links are invisible over SMB")
            self.assertEqual(os.path.realpath(link), os.path.realpath(f))
        self.assertTrue(os.path.isfile(os.path.join(project, "SUBMISSION.md")))
        self.assertEqual(load(os.path.join(project, cs.MARKER))["id"], HEX)
        self.assertEqual(load(os.path.join(project, cs.MARKER))["share_dir"], SERVER_SHARE)
        self.assertTrue(os.path.isdir(os.path.join(self.tmp, "work", "on_campus", "Placeholder lab", "PROT_0807")))
        rc, out, _ = self.stage(summary, "--apply")
        self.assertEqual(rc, 0)
        self.assertEqual((out["links_created"], out["links_existing"]), (0, 2))

    def test_no_matching_group_proposes_a_new_folder(self):
        self.folder("on_campus", "Someone Else")
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.stage(summary)
        self.assertEqual(rc, 0)
        self.assertEqual(out["group_resolution"]["mode"], "new")
        self.assertEqual(out["group_dir"], os.path.join(self.sroot, "on_campus", "Placeholder"))

    def test_several_matching_groups_ask_a_human(self):
        self.folder("on_campus", "Placeholder lab")
        self.folder("on_campus", "Placeholder-old")
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.stage(summary)
        self.assertEqual(rc, 2)
        self.assertEqual(len(out["group_resolution"]["candidates"]), 2)
        rc, out, _ = self.stage(summary, "--service-dir", "Placeholder lab")
        self.assertEqual(rc, 0)

    def test_a_same_surname_folder_naming_someone_else_is_ambiguous(self):
        """Defect 8: PI Ying Wang was filed under `Wang Wei` with exit 0."""
        self.folder("on_campus", "Wang Wei")
        summary, _ = write_summary(self.tmp, record(pi_last="Wang", pi_first="Ying"))
        rc, out, _ = self.stage(summary)
        self.assertEqual(rc, 2)
        self.assertIn("wei", out["group_resolution"]["reason"])

    def test_a_same_surname_folder_naming_the_pi_is_reused(self):
        group = self.folder("on_campus", "Wang Ying")
        summary, _ = write_summary(self.tmp, record(pi_last="Wang", pi_first="Ying"))
        rc, out, _ = self.stage(summary)
        self.assertEqual((rc, out["group_dir"]), (0, group))

    def test_surname_particles_do_not_match_an_institution_folder(self):
        """Defect 8: "de la Cruz" was filed under UC_Santa_Cruz."""
        self.folder("off_campus", "UC_Santa_Cruz")
        summary, _ = write_summary(self.tmp, record(pi_last="de la Cruz", pi_first="Quinn",
                                                    institution="Example University"))
        rc, out, _ = self.stage(summary)
        self.assertEqual(rc, 0)
        self.assertNotIn("UC_Santa_Cruz", out["group_dir"])

    def test_institution_folders_must_match_on_distinctive_words(self):
        """Defect 8: University of Washington -> Washington_State_Univ, UT Southwestern -> Texas_AM."""
        self.folder("off_campus", "Washington_State_Univ")
        self.folder("off_campus", "Texas_AM")
        for inst, wrong in (("University of Washington", "Washington_State_Univ"),
                            ("University of Texas Southwestern Medical Center", "Texas_AM")):
            summary, _ = write_summary(self.tmp, record(institution=inst))
            rc, out, _ = self.stage(summary)
            self.assertEqual(rc, 0, inst)
            self.assertNotIn(wrong, out["group_dir"], inst)

    def test_service_dir_outside_the_service_root_is_refused(self):
        summary, _ = write_summary(self.tmp)
        self.assertEqual(self.stage(summary, "--service-dir", self.tmp)[0], 2)

    def test_a_project_folder_owned_by_another_submission_is_refused(self):
        project = self.folder("on_campus", "Placeholder lab", "PROT_0807")
        with open(os.path.join(project, cs.MARKER), "w") as fh:
            json.dump({"id": "bbbbbbbbbbbb", "internal_id": "PROT_0807"}, fh)
        summary, _ = write_summary(self.tmp)
        self.assertEqual(self.stage(summary, "--apply")[0], 2)

    def test_off_campus_institution_alias_finds_the_institution_folder(self):
        self.folder("off_campus", "UCSF", "Other_lab")
        summary, _ = write_summary(self.tmp, record(institution="UC San Francisco"))
        rc, out, _ = self.stage(summary)
        self.assertEqual(rc, 0)
        self.assertEqual(out["group_dir"], os.path.join(self.sroot, "off_campus", "UCSF", "Placeholder"))

    def test_apply_refuses_files_whose_locate_hard_failed(self):
        """Defect 15: files.txt is written even on a FAIL; nothing else stopped it being staged."""
        self.folder("on_campus", "Placeholder lab")
        summary, _ = write_summary(self.tmp)
        locate_json = os.path.join(self.tmp, "locate.json")
        with open(locate_json, "w") as fh:
            json.dump({"mode": "auto", "hard_fail": True, "gates": [{"gate": "weak_ids", "status": "FAIL"}]}, fh)
        self.assertEqual(self.stage(summary)[0], 0, "a dry run may still show the plan")
        rc, out, _ = self.stage(summary, "--apply")
        self.assertEqual(rc, 2)
        self.assertEqual(out["failing_gates"][0]["gate"], "weak_ids")
        with open(locate_json, "w") as fh:
            json.dump({"mode": "files_from", "hard_fail": True, "gates": []}, fh)
        self.assertEqual(self.stage(summary, "--apply")[0], 0)

    def test_a_hand_written_submission_md_is_kept(self):
        """Defect 16."""
        project = self.folder("on_campus", "Placeholder lab", "PROT_0807")
        write(project, "SUBMISSION.md", "# notes typed by staff\n")
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.stage(summary, "--apply")
        self.assertEqual(rc, 0)
        self.assertEqual(read(os.path.join(project, "SUBMISSION.md")), "# notes typed by staff\n")
        self.assertTrue(read(os.path.join(project, "SUBMISSION.core.md")).startswith(cs.GENERATED_MARK))
        rc, out, _ = self.stage(summary, "--apply")                  # ours is overwritten happily
        self.assertEqual(rc, 0)


# -------------------------------------------------------------------- conditions --
class TestConditions(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.env = child_env(self.tmp)
        self.summary, _ = write_summary(self.tmp)
        self.rawdir = os.path.join(self.tmp, "rawdir")
        os.makedirs(self.rawdir)

    def tearDown(self):
        self._tmp.cleanup()

    def conditions(self, pairs, status="matched"):
        tsv = os.path.join(self.tmp, "sample_files.tsv")
        with open(tsv, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["unique_id", "sample_name", "condition_name", "file", "acquired", "date_source",
                        "status", "alternates", "note"])
            for i, (uid, cond) in enumerate(pairs):
                f = os.path.join(self.rawdir, f"03122025__60SPD_DIA-{uid}_S3-A{i}_1_2400{i}.d")
                os.makedirs(f, exist_ok=True)
                w.writerow([uid, f"sample {uid}", cond, f, "2025-03-12", "filename", status, "", ""])
        out = os.path.join(self.tmp, "conditions.csv")
        rc, js, p = run(["conditions", "--summary", self.summary, "--sample-files", tsv, "--out", out], self.env)
        return rc, js, out, p

    def test_clean_two_groups_writes_the_run_names_collect_conditions_uses(self):
        rc, js, out, p = self.conditions([("KG1", "ctrl"), ("KG2", "ctrl"), ("KG13", "treat"), ("KG14", "treat")])
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        self.assertFalse(js["needs_user_input"])
        with open(out, newline="") as fh:
            rows = list(csv.DictReader(fh))
        listed = json.loads(subprocess.run(
            [sys.executable, os.path.join(SCRIPTS, "collect_conditions.py"), "--list-runs",
             "--from-dir", self.rawdir, "--glob", "*.d"], capture_output=True, text=True).stdout)["runs"]
        self.assertEqual(sorted(r["File.Name"] for r in rows), sorted(listed))
        self.assertEqual(dict(Counter_(r["Group"] for r in rows)), {"ctrl": 2, "treat": 2})

    def test_blank_condition_needs_input(self):
        rc, js, _, _ = self.conditions([("KG1", "ctrl"), ("KG2", ""), ("KG13", "treat"), ("KG14", "treat")])
        self.assertEqual(rc, 2)
        self.assertEqual(js["findings"]["blank_conditions"], ["KG2"])
        self.assertTrue(js["questions"])

    def test_single_condition_needs_input(self):
        rc, js, _, _ = self.conditions([("KG1", "all"), ("KG2", "all"), ("KG13", "all")])
        self.assertEqual(rc, 2)
        self.assertIn("single_condition", js["findings"])

    def test_every_condition_unique_needs_input(self):
        rc, js, _, _ = self.conditions([("KG1", "KG1"), ("KG2", "KG2"), ("KG13", "KG13")])
        self.assertEqual(rc, 2)
        self.assertTrue(js["findings"]["all_unique"])

    def test_singleton_group_needs_input(self):
        rc, js, _, _ = self.conditions([("KG1", "ctrl"), ("KG2", "ctrl"), ("KG13", "treat"),
                                        ("KG14", "treat"), ("KG15", "odd")])
        self.assertEqual(rc, 2)
        self.assertEqual(js["findings"]["singleton_groups"], ["odd"])

    def test_conditions_differing_only_in_case_need_input(self):
        """Defect 9: Control x2, control x2, Treated x2 read as three groups with exit 0."""
        rc, js, _, _ = self.conditions([("KG1", "Control"), ("KG2", "Control"), ("KG13", "control"),
                                        ("KG14", "control"), ("KG15", "Treated"), ("KG16", "Treated")])
        self.assertEqual(rc, 2)
        self.assertEqual(js["findings"]["case_variants"], {"control": ["Control", "control"]})


def Counter_(it):
    d = {}
    for x in it:
        d[x] = d.get(x, 0) + 1
    return d


# ----------------------------------------------------------------------- deliver --
class DeliverBase(unittest.TestCase):
    """A staged submission with a finished session under its work dir."""
    rec_kwargs = {}

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.env = child_env(self.tmp)
        self.summary, self.s = write_summary(self.tmp, record(**self.rec_kwargs))
        sroot = os.path.join(self.tmp, "flinders", "Data", "lab", "service")
        os.makedirs(os.path.join(sroot, "on_campus", "Placeholder lab"))
        self.raw_files = [tims(self.tmp, "03122025", u, i) for i, u in enumerate(("KG1", "KG2"), 1)]
        flist = write(self.tmp, "files.txt", "\n".join(self.raw_files) + "\n")
        rc, st, p = run(["stage", "--summary", self.summary, "--files", flist, "--apply"], self.env)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.work, self.service_project = st["work_dir"], st["service_dir"]
        self.session = os.path.join(self.work, "sessions", "2025-04-01_study")
        o = os.path.join(self.session, "output")
        for d in ("tables", "figures", "reproducibility", "search/xic", "raw_data"):
            os.makedirs(os.path.join(o, d), exist_ok=True)
        os.makedirs(os.path.join(self.session, "logs"))
        write(o, "Analysis_Report.html", "<html>report</html>")
        write(o, "AUDIT.md", "# audit")
        # The DE table is a symlink to a file elsewhere: the delivery must hold the BYTES.
        elsewhere = write(self.tmp, "elsewhere/DE_dpc_ctrl.vs.treat.csv", "Protein,logFC\nP1,1\n")
        os.symlink(elsewhere, os.path.join(o, "tables", "DE_dpc_ctrl.vs.treat.csv"))
        write(o, "tables/reproducibility_log.R", "# R")
        write(o, "figures/pca.png", "png")
        write(o, "figures/tmp/scratch.png", "scratch")
        write(o, "tables/run1.quant", "quant")
        write(o, "search/report.parquet", "parquet")
        write(o, "search/xic/t1_xic.parquet", "xic")
        self.output = o
        self.share = local_share(self.tmp)
        self.delivery = os.path.join(self.share, "PROT_0807_t1")

    def tearDown(self):
        for dp, dns, fns in os.walk(self.tmp):               # undo chmods so cleanup can remove
            for n in dns:
                p = os.path.join(dp, n)
                if not os.path.islink(p):
                    os.chmod(p, 0o755)
        self._tmp.cleanup()

    def deliver(self, *extra, env=None):
        args = ["deliver", "--summary", self.summary, "--session", self.session, "--label", "t1"]
        if "--include-raw" not in extra:
            args += ["--include-raw", "no"]
        return run(args + list(extra), env or self.env)


class TestDeliver(DeliverBase):
    def test_dry_run_copies_nothing(self):
        rc, out, p = self.deliver()
        self.assertEqual(rc, 0, p.stderr)
        self.assertFalse(out["applied"])
        self.assertEqual(out["share_dir"], SERVER_SHARE)
        self.assertFalse(os.path.exists(self.share))

    def test_apply_dereferences_records_skips_and_verifies(self):
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        de = os.path.join(self.delivery, "tables", "DE_dpc_ctrl.vs.treat.csv")
        self.assertTrue(os.path.isfile(de) and not os.path.islink(de))
        self.assertIn("logFC", read(de))
        manifest = read(os.path.join(self.delivery, "MANIFEST.txt"))
        self.assertIn("[OK] Analysis_Report.html", manifest)
        self.assertRegex(manifest, r"\[SKIPPED\] methods\.docx -- ")
        self.assertIn("[SKIPPED] figures/tmp/ -- excluded", manifest)
        self.assertIn("[SKIPPED] tables/run1.quant -- excluded", manifest)
        self.assertFalse(os.path.exists(os.path.join(self.delivery, "search", "xic")))
        self.assertFalse(os.path.exists(os.path.join(self.delivery, cs.DELIVERY_MARKER)))
        self.assertIn("  Analysis_Report.html", read(os.path.join(self.delivery, "checksums.sha256")))
        readme = read(os.path.join(self.delivery, "README.md"))
        self.assertIn("Analysis_Report.html", readme)
        self.assertNotIn("methods.docx", readme, "never describe a file that was not delivered")
        self.assertNotIn("expression matrix", readme)
        self.assertNotIn("@", readme, "no email addresses in a collaborator README")
        self.assertNotIn(self.tmp, readme, "no internal paths in a collaborator README")
        self.assertFalse(os.path.exists(os.path.join(self.share, "raw")))
        dj = load(os.path.join(self.session, "delivery.json"))
        self.assertTrue(dj["verified"])
        self.assertEqual((dj["share_dir"], dj["internal_id"]), (SERVER_SHARE, "PROT_0807"))

    def test_missing_analysis_report_blocks(self):
        os.remove(os.path.join(self.output, "Analysis_Report.html"))
        rc, out, _ = self.deliver("--apply")
        self.assertEqual(rc, 2)
        self.assertFalse(os.path.exists(self.delivery))

    def test_size_guard_writes_a_job_instead_of_copying(self):
        rc, out, _ = self.deliver("--apply", "--max-gb", "0.000000001")
        self.assertEqual(rc, 5)
        job = read(out["deliver_job"])
        for want in ("#SBATCH", "--apply", "--no-size-guard", "--session"):
            self.assertIn(want, job)
        self.assertFalse(os.path.exists(self.delivery))

    def test_raw_links_are_relative_and_the_service_project_links_back(self):
        rc, out, p = self.deliver("--apply", "--include-raw", "yes")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        for f in self.raw_files:
            link = os.path.join(self.share, "raw", os.path.basename(f))
            self.assertFalse(os.path.isabs(os.readlink(link)))
            self.assertEqual(os.path.realpath(link), os.path.realpath(f))
        back = out["results_link"]["link"]
        self.assertFalse(os.path.isabs(os.readlink(back)))
        self.assertEqual(os.path.realpath(back), os.path.realpath(self.delivery))
        self.assertIn(cs.RAW_WHITELIST_WARNING, out["warnings"])
        self.assertTrue(out["raw_whitelist_unverified"])

    # -- defect 1: the session must be this submission's
    def test_a_session_from_elsewhere_is_refused_unless_it_searched_the_staged_files(self):
        other = os.path.join(self.tmp, "somewhere", "session")
        shutil.copytree(self.session, other, symlinks=True)
        write(other, "input/raw_files.txt", "# Raw MS files\n/nfs/other/lab_run.d\n")
        rc, out, _ = run(["deliver", "--summary", self.summary, "--session", other, "--include-raw", "no"], self.env)
        self.assertEqual(rc, 2)
        self.assertIn("does not belong", out["error"])
        write(other, "input/raw_files.txt", "# Raw MS files\n" + "\n".join(self.raw_files) + "\n")
        rc, out, _ = run(["deliver", "--summary", self.summary, "--session", other, "--include-raw", "no"], self.env)
        self.assertEqual(rc, 0)

    def test_a_session_under_another_submissions_work_dir_is_refused(self):
        other = os.path.join(self.tmp, "work", "on_campus", "Other", "PROT_0806")
        write(other, cs.MARKER, json.dumps({"id": "bbbbbbbbbbbb", "internal_id": "PROT_0806"}))
        session = os.path.join(other, "sessions", "s")
        shutil.copytree(self.session, session, symlinks=True)
        rc, out, _ = run(["deliver", "--summary", self.summary, "--session", session, "--include-raw", "no"], self.env)
        self.assertEqual(rc, 2)
        self.assertIn("PROT_0806", out["error"])

    def test_without_a_stage_record_force_is_required(self):
        os.remove(os.path.join(self.service_project, cs.MARKER))
        os.remove(os.path.join(self.work, cs.MARKER))
        rc, out, _ = self.deliver()
        self.assertEqual(rc, 2)
        rc, out, _ = self.deliver("--force")
        self.assertEqual(rc, 0)
        self.assertIn("NOT CHECKED", out["session_ownership"])

    # -- defect 5: never mix deliveries; verify knows the plan
    def test_a_non_empty_delivery_folder_is_refused(self):
        write(self.delivery, "old_results.csv", "stale")
        rc, out, _ = self.deliver("--apply")
        self.assertEqual(rc, 2)
        self.assertIn("--label", out["hint"])

    def test_verify_fails_on_a_file_outside_the_plan(self):
        with mock.patch.dict(os.environ, roots(self.tmp)):
            write(self.delivery, "README.md", "x")
            write(self.delivery, "stale.csv", "x")
            for p in (self.delivery, os.path.join(self.delivery, "README.md"), os.path.join(self.delivery, "stale.csv")):
                os.chmod(p, 0o775 if os.path.isdir(p) else 0o664)
            bad, _ = cs.verify_delivery(self.delivery, self.share, cs.flinders_root(), {"README.md"})
        self.assertEqual(len(bad), 1)
        self.assertIn("stale.csv", bad[0])

    # -- defect 6: symlinks anywhere in the share, and never writing through one
    def test_an_existing_link_elsewhere_in_the_share_is_noted_and_left_untouched(self):
        # PROT_0793: share/search -> /quobyte links are valid (they stopped serving only because of
        # a Bioshare https setting). They must not block a later delivery, and are never "repaired".
        os.makedirs(self.share)
        link = os.path.join(self.share, "search")
        os.symlink("/quobyte/proteomics-grp/nowhere/search", link)
        rc, out, _ = self.deliver("--apply")
        self.assertEqual(rc, 0, out)
        self.assertTrue(out["verified"])
        self.assertTrue(any("search" in w and "left untouched" in w for w in out["warnings"]), out["warnings"])
        self.assertEqual(os.readlink(link), "/quobyte/proteomics-grp/nowhere/search")

    def test_a_symlink_inside_the_delivery_is_a_violation(self):
        with mock.patch.dict(os.environ, roots(self.tmp)):
            os.makedirs(self.delivery)
            os.symlink(self.raw_files[0], os.path.join(self.delivery, "sneaky.d"))
            bad, _ = cs.verify_delivery(self.delivery, self.share, cs.flinders_root())
            self.assertTrue(any("inside this delivery" in b for b in bad), bad)

    def test_a_share_dir_that_is_a_symlink_out_of_the_root_is_refused_before_writing(self):
        outside = os.path.join(self.tmp, "quobyte_elsewhere")
        os.makedirs(outside)
        os.makedirs(os.path.dirname(self.share))
        os.symlink(outside, self.share)
        rc, out, _ = self.deliver("--apply")
        self.assertEqual(rc, 2)
        self.assertEqual(os.listdir(outside), [])

    def test_writes_never_follow_a_symlinked_subfolder(self):
        victim_dir = os.path.join(self.tmp, "victim")
        victim = write(victim_dir, "DE_dpc_ctrl.vs.treat.csv", "ORIGINAL")
        os.makedirs(self.delivery)
        write(self.delivery, cs.DELIVERY_MARKER, json.dumps({"session": self.session, "folder": "PROT_0807_t1",
                                                             "mode": "analysis"}))
        os.symlink(victim_dir, os.path.join(self.delivery, "tables"))
        rc, out, _ = self.deliver("--apply")
        self.assertEqual(rc, 2)
        self.assertEqual(read(victim), "ORIGINAL", "a file outside the share was overwritten")

    def test_verify_flags_an_absolute_raw_link_it_created_but_only_notes_an_existing_one(self):
        with mock.patch.dict(os.environ, roots(self.tmp)):
            os.makedirs(os.path.join(self.share, "raw"))
            os.makedirs(self.delivery)
            os.symlink(self.raw_files[0], os.path.join(self.share, "raw", "abs.d"))
            bad, notes = cs.verify_delivery(self.delivery, self.share, cs.flinders_root(), own_raw={"abs.d"})
            self.assertTrue(any("absolute link" in b for b in bad))
            bad, notes = cs.verify_delivery(self.delivery, self.share, cs.flinders_root())
            self.assertEqual(bad, [])
            self.assertTrue(any("abs.d" in n and "left untouched" in n for n in notes), notes)

    # -- defect 7: an error mid-delivery still leaves a MANIFEST and fails loudly
    def test_a_crash_after_copying_still_writes_the_manifest_and_fails(self):
        stdout = io.StringIO()
        with mock.patch.dict(os.environ, roots(self.tmp)), \
                mock.patch.object(cs, "build_readme", side_effect=RuntimeError("simulated crash")), \
                contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(io.StringIO()):
            rc = cs.main(["deliver", "--summary", self.summary, "--session", self.session, "--label", "t1",
                          "--include-raw", "no", "--apply"])
        self.assertEqual(rc, 2)
        out = json.loads(stdout.getvalue())
        self.assertFalse(out["verified"])
        self.assertIn("simulated crash", out["error"])
        self.assertIn("stopped by an error", read(os.path.join(self.delivery, "MANIFEST.txt")))
        self.assertFalse(load(os.path.join(self.session, "delivery.json"))["verified"])

    @unittest.skipIf(IS_ROOT, "root ignores directory permissions")
    def test_a_failed_results_link_is_recorded_not_a_crash(self):
        os.chmod(self.service_project, 0o555)
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout)
        self.assertEqual(out["results_link"]["state"], "SKIPPED")
        self.assertIn("[SKIPPED] results link", read(os.path.join(self.delivery, "MANIFEST.txt")))

    # -- defect 12 and 13
    def test_delivered_files_are_group_and_world_readable_under_umask_077(self):
        os.chmod(os.path.join(self.output, "AUDIT.md"), 0o600)
        old = os.umask(0o077)
        try:
            rc, out, p = self.deliver("--apply")
        finally:
            os.umask(old)
        self.assertEqual(rc, 0, p.stdout)
        for dp, dns, fns in os.walk(self.delivery):
            for n in fns:
                self.assertEqual(os.stat(os.path.join(dp, n)).st_mode & 0o044, 0o044, n)
            for n in dns:
                self.assertEqual(os.stat(os.path.join(dp, n)).st_mode & 0o055, 0o055, n)

    @unittest.skipIf(IS_ROOT, "root can read anything")
    def test_an_unreadable_subfolder_is_recorded_not_skipped_silently(self):
        locked = os.path.join(self.output, "tables", "locked")
        write(locked, "hidden.csv", "x")
        os.chmod(locked, 0)
        rc, out, _ = self.deliver("--apply")
        manifest = read(os.path.join(self.delivery, "MANIFEST.txt"))
        self.assertIn("[SKIPPED] tables/locked/ -- unreadable directory", manifest)

    # -- live finding C
    def test_readme_claims_only_what_was_delivered(self):
        bare = cs.build_readme(self.s, {"Analysis_Report.html"}, [], "analysis")
        self.assertNotIn("compared", bare)
        self.assertNotIn("differ", bare)
        self.assertNotIn("searched", bare)
        full = cs.build_readme(self.s, {"Analysis_Report.html", "tables/DE_dpc_a.vs.b.csv", "search/report.parquet"},
                               [], "analysis")
        self.assertIn("compared between sample groups", full)
        self.assertIn("searched", full)


class TestRawOnlyDelivery(DeliverBase):
    """Defect 11: 'I only require raw data' gets raw links and a README that says so -- no search."""
    rec_kwargs = {"raw_only": True}

    def test_raw_only_needs_no_session_and_no_report(self):
        rc, out, p = run(["deliver", "--summary", self.summary, "--label", "raw", "--apply"], self.env)
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual(out["mode"], "raw-only")
        delivery = os.path.join(self.share, "PROT_0807_raw")
        readme = read(os.path.join(delivery, "README.md"))
        self.assertIn("raw data only", readme)
        self.assertNotIn("Analysis_Report.html", readme)
        for f in self.raw_files:
            self.assertTrue(os.path.islink(os.path.join(self.share, "raw", os.path.basename(f))))
        self.assertIn(cs.RAW_WHITELIST_WARNING, out["warnings"])
        self.assertTrue(load(os.path.join(self.work, "delivery.json"))["verified"])


# ---------------------------------------------------------------------- bioshare --
class TestBioshare(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.summary, self.s = write_summary(self.tmp)

    def tearDown(self):
        self._tmp.cleanup()

    def bs(self, fake, *args, **env):
        return run(["bioshare", *args, "--summary", self.summary], child_env(self.tmp, fake.base, **env))

    def linked_share(self):
        return {"id": 41, "submission": HEX, "bioshare_id": "bsx", "name": "Placeholder, Quinn: PROT_0807",
                "notes": "", "sub_folder": None, "link_to_path": self.s["share_dir"] + "/",
                "url": "https://bioshare.example/bsx/"}

    def delivery(self, **over):
        d = dict({"verified": True, "internal_id": "PROT_0807", "share_dir": SERVER_SHARE}, **over)
        return write(self.tmp, "delivery.json", json.dumps(d))

    def test_status_marks_the_share_linked_to_this_share_dir(self):
        with FakeCoreOmics() as fake:
            other = dict(self.linked_share(), id=40, link_to_path="/somewhere/else")
            fake.routes[("GET", SHARES)] = lambda q, b: (200, {"results": [other, self.linked_share()],
                                                               "next": None})
            rc, out, p = self.bs(fake, "status")
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual([s["linked"] for s in out["shares"]], [False, True])
        self.assertEqual(out["url"], "https://bioshare.example/bsx/")

    def test_ensure_dry_run_posts_nothing(self):
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [])
            rc, out, _ = self.bs(fake, "ensure")
            self.assertEqual(rc, 0)
            self.assertEqual(fake.posts(), [])
            self.assertEqual(out["would_post"]["json"]["link_to_path"], SERVER_SHARE)

    def test_ensure_apply_posts_the_exact_payload(self):
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [])
            fake.routes[("POST", SHARES)] = lambda q, b: (201, dict(self.linked_share(), **b))
            rc, out, _ = self.bs(fake, "ensure", "--apply")
            self.assertEqual(rc, 0)
            self.assertEqual(fake.posts()[0]["body"], {
                "submission": HEX, "name": "Placeholder, Quinn: PROT_0807",
                "notes": "UC Davis Proteomics Core results for PROT_0807",
                "link_to_path": SERVER_SHARE})
            self.assertEqual(fake.posts()[0]["auth"], f"Token {TOKEN}")
            self.assertEqual(out["url"], "https://bioshare.example/bsx/")

    def test_a_local_or_windows_root_never_reaches_link_to_path(self):
        """Defect 10: link_to_path comes from the summary's server path, not CORE_FLINDERS_ROOT."""
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [])
            rc, out, _ = self.bs(fake, "ensure", CORE_FLINDERS_ROOT="R:\\proteomics")
            self.assertEqual(rc, 0)
            self.assertEqual(out["would_post"]["json"]["link_to_path"], SERVER_SHARE)

    def test_ensure_requires_a_share_record_back(self):
        """Defect 14."""
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [])
            fake.routes[("POST", SHARES)] = lambda q, b: (201, [])
            self.assertEqual(self.bs(fake, "ensure", "--apply")[0], 3)

    def test_send_dry_run_lists_submitter_pi_and_contacts_and_posts_nothing(self):
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [self.linked_share()])
            rc, out, _ = self.bs(fake, "send")
            self.assertEqual(rc, 0)
            self.assertEqual(fake.posts(), [])
            r = out["recipients"]
            self.assertEqual(r["submitter"]["email"], "ada.example@example.org")
            self.assertEqual(r["pi"]["email"], "pi.placeholder@example.org")
            self.assertEqual(r["contacts"], [{"name": "Robin Contact", "email": "robin.contact@example.org"}])

    def test_send_apply_requires_a_verified_delivery_of_this_share(self):
        """Defect 15."""
        share_url = SHARES + "41/share/"
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [self.linked_share()])
            fake.routes[("POST", share_url)] = lambda q, b: (200, {"ok": True})
            self.assertEqual(self.bs(fake, "send", "--apply")[0], 2)
            self.assertEqual(self.bs(fake, "send", "--apply", "--delivery", self.delivery(verified=False))[0], 2)
            self.assertEqual(self.bs(fake, "send", "--apply", "--delivery",
                                     self.delivery(share_dir="/nfs/lssc0/other/share"))[0], 2)
            self.assertEqual(self.bs(fake, "send", "--apply", "--delivery", self.delivery(internal_id="PROT_0806"))[0], 2)
            self.assertEqual(fake.posts(), [])
            self.assertEqual(self.bs(fake, "send", "--apply", "--delivery", self.delivery())[0], 0)
            self.assertEqual(len(fake.posts()), 1)

    def test_send_apply_defaults_to_no_notification_email(self):
        share_url = SHARES + "41/share/"
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [self.linked_share()])
            fake.routes[("POST", share_url)] = lambda q, b: (200, {"ok": True})
            d = self.delivery()
            self.assertEqual(self.bs(fake, "send", "--apply", "--delivery", d)[0], 0)
            self.assertEqual(self.bs(fake, "send", "--apply", "--delivery", d, "--email")[0], 0)
            self.assertEqual([p["body"] for p in fake.posts()], [{"email": False}, {"email": True}])

    def test_a_redirected_post_is_an_error_not_a_success(self):
        """Defect 14: urllib turned POST+302 into a GET of the listing and reported applied:true."""
        share_url = SHARES + "41/share/"
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [self.linked_share()])
            fake.routes[("POST", share_url)] = lambda q, b: (302, {}, {"Location": SHARES})
            rc, out, _ = self.bs(fake, "send", "--apply", "--delivery", self.delivery())
            self.assertEqual(rc, 3)
            self.assertEqual(out["status"], 302)
            self.assertEqual(len([r for r in fake.requests if r["method"] == "GET"]), 1,
                             "the redirect must not be followed")

    def test_send_without_a_linked_share_asks_for_ensure(self):
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [])
            self.assertEqual(self.bs(fake, "send", "--apply", "--delivery", self.delivery())[0], 2)
            self.assertEqual(fake.posts(), [])

    def test_http_400_detail_surfaces_with_exit_3(self):
        with FakeCoreOmics() as fake:
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [])
            fake.routes[("POST", SHARES)] = lambda q, b: (400, {"link_to_path": ["Path not allowed."]})
            rc, out, _ = self.bs(fake, "ensure", "--apply")
        self.assertEqual(rc, 3)
        self.assertEqual(out["status"], 400)
        self.assertIn("Path not allowed.", out["detail"])


# ------------------------------------------------------------------------- fetch --
class TestFetch(unittest.TestCase):
    def serve(self, fake, rec, shares=(200, [])):
        def listing(q, b):
            if q.get("internal_id"):
                hits = [rec] if q["internal_id"] == rec["internal_id"] else []
                return 200, {"count": len(hits), "next": None, "results": hits}
            page = int(q.get("page", 1))
            pages = {
                1: [record(internal_id="PROT_0900", sid="cccccccccccc", submitted="2025-12-01T10:00:00-08:00",
                           samples=[("KG1", "x")]),
                    record(internal_id="PROT_0810", sid="aaaaaaaaaaaa", submitted="2025-04-01T10:00:00-07:00",
                           samples=[("KG1", "x"), ("BN1", "y")]),
                    rec,
                    record(internal_id="PROT_0700", sid="dddddddddddd", submitted="2024-09-01T10:00:00-07:00",
                           samples=[("SG001", "x")])],
                2: [record(internal_id="PROT_0500", sid="eeeeeeeeeeee", submitted="2024-05-01T10:00:00-07:00")],
                3: [record(internal_id="PROT_0400", sid="ffffffffffff", submitted="2024-01-01T10:00:00-08:00")],
            }
            nxt = f"{fake.base}/submissions/?page={page + 1}" if page < 3 else None
            return 200, {"count": 6, "next": nxt, "results": pages[page]}
        fake.routes[("GET", "/server/api/submissions/")] = listing
        fake.routes[("GET", SHARES)] = lambda q, b: shares

    def test_fetch_writes_three_files_and_the_neighbours_in_the_window(self):
        rec = record()
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, rec, shares=(200, [{"id": 1, "link_to_path": "/x", "url": "u"}]))
            out_dir = os.path.join(tmp, "PROT_0807")
            rc, out, p = run(["fetch", "#807", "--out", out_dir], child_env(tmp, fake.base))
            self.assertEqual(rc, 0, p.stderr)
            for f in ("submission.json", "samples.tsv", "submission_summary.json"):
                self.assertTrue(os.path.isfile(os.path.join(out_dir, f)), f)
            s = load(os.path.join(out_dir, "submission_summary.json"))
            self.assertEqual({n["internal_id"] for n in s["neighbors"]}, {"PROT_0810", "PROT_0700"})
            self.assertEqual(s["neighbors"][0]["unique_ids"], ["KG1", "BN1"])
            self.assertFalse(any(r["query"].get("page") == "3" for r in fake.requests),
                             "paging must stop once past the old edge of the window")
            self.assertTrue(all(r["auth"] == f"Token {TOKEN}" for r in fake.requests))
            self.assertEqual(s["organism_as_submitted"]["value"], "mouse")
            self.assertFalse(s["organism_as_submitted"]["confirmed"])
            self.assertEqual(s["conditions"]["groups"], {"ctrl": 2, "treat": 2})
            self.assertEqual(s["campus"], "on_campus")
            self.assertEqual(s["share_dir"], SERVER_SHARE)
            self.assertEqual(s["contacts"][0]["email"], "robin.contact@example.org")
            self.assertEqual(len(s["existing_shares"]), 1)
            self.assertEqual(read(os.path.join(out_dir, "samples.tsv")).split("\n")[0].strip(),
                             "unique_id\tsample_name\tcondition_name")

    def test_a_share_listing_error_is_recorded_not_fatal(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, record(), shares=(403, {"detail": "You do not have permission."}))
            rc, out, _ = run(["fetch", "807", "--out", tmp], child_env(tmp, fake.base))
            self.assertEqual(rc, 0)
            s = load(os.path.join(tmp, "submission_summary.json"))
            self.assertIn("permission", s["existing_shares_error"])

    def test_no_such_submission_exits_2(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, record())
            self.assertEqual(run(["fetch", "808", "--out", tmp], child_env(tmp, fake.base))[0], 2)

    def test_unreachable_coreomics_exits_3(self):
        with tempfile.TemporaryDirectory() as tmp:
            with FakeCoreOmics() as fake:
                base = fake.base                     # a port that is closed once the server stops
            self.assertEqual(run(["fetch", "807", "--out", tmp], child_env(tmp, base))[0], 3)


# ------------------------------------------------------------------- email-draft --
class TestEmailDraft(unittest.TestCase):
    def test_draft_names_the_link_and_only_numbers_it_was_given(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary, _ = write_summary(tmp)
            delivery = write(tmp, "delivery.json", json.dumps(
                {"mode": "analysis", "contents": ["Analysis_Report.html", "tables"], "n_files": 12,
                 "bytes": 5 * 1024 ** 2, "raw_links": {"created": 0}}))
            out = os.path.join(tmp, "EMAIL_DRAFT.md")
            rc, js, p = run(["email-draft", "--summary", summary, "--delivery", delivery,
                             "--share-url", "https://bioshare.example/bsx/", "--out", out], child_env(tmp))
            self.assertEqual(rc, 0, p.stderr)
            self.assertFalse(js["sent"])
            text = read(out)
            self.assertIn("To: ada.example@example.org", text)
            self.assertIn("Cc: pi.placeholder@example.org", text)
            self.assertIn("https://bioshare.example/bsx/", text)
            self.assertIn("Analysis_Report.html", text)
            self.assertIn("12 files, 5.0 MB", text)
            self.assertNotIn("raw folder", text, "no raw links were made, so none are promised")
            self.assertNotIn(tmp, text)

    def test_raw_only_draft_does_not_promise_an_analysis(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary, _ = write_summary(tmp, record(raw_only=True))
            delivery = write(tmp, "delivery.json", json.dumps(
                {"mode": "raw-only", "contents": ["MANIFEST.txt", "README.md"], "raw_links": {"created": 4}}))
            out = os.path.join(tmp, "EMAIL_DRAFT.md")
            rc, _, _ = run(["email-draft", "--summary", summary, "--delivery", delivery,
                            "--share-url", "https://bioshare.example/bsx/", "--out", out], child_env(tmp))
            self.assertEqual(rc, 0)
            text = read(out)
            self.assertIn("raw instrument files", text)
            self.assertNotIn("Analysis_Report.html", text)

    def test_missing_link_leaves_a_visible_placeholder_and_exits_2(self):
        with tempfile.TemporaryDirectory() as tmp:
            summary, _ = write_summary(tmp)
            out = os.path.join(tmp, "EMAIL_DRAFT.md")
            rc, js, _ = run(["email-draft", "--summary", summary, "--out", out], child_env(tmp))
            self.assertEqual(rc, 2)
            self.assertIn("Bioshare link", js["placeholders"])
            self.assertIn("[BIOSHARE LINK", read(out))


# --------------------------------------------------------------------- hive_exec --
class TestHiveExec(unittest.TestCase):
    """Connection reuse keeps rapid calls under HIVE's MaxStartups throttle; --get must not
    report a good copy as failed."""

    def _argv(self, *cmd, **env_extra):
        with tempfile.TemporaryDirectory() as tmp:
            log = os.path.join(tmp, "argv")
            for tool in ("ssh", "rsync"):
                fake = os.path.join(tmp, tool)
                with open(fake, "w") as fh:
                    fh.write(f'#!/bin/sh\necho "{tool}" > "{log}"\nfor a in "$@"; do echo "$a"; done >> "{log}"\n')
                os.chmod(fake, 0o755)
            key = os.path.join(tmp, "key")
            open(key, "w").close()
            env = dict(os.environ, PATH=tmp + os.pathsep + os.environ.get("PATH", ""), HOME=tmp,
                       HIVE_USER="someone", HIVE_KEY=key, HIVE_ENV_FILE=os.path.join(tmp, "none"), **env_extra)
            bash = "/bin/bash" if os.path.exists("/bin/bash") else "bash"   # macOS: bash 3.2
            p = subprocess.run([bash, os.path.join(SCRIPTS, "hive_exec.sh"), *(cmd or ("hostname",))],
                               capture_output=True, text=True, env=env)
            self.assertEqual(p.returncode, 0, p.stderr)
            return read(log).split("\n")

    def test_mux_options_are_on_by_default(self):
        argv = self._argv()
        self.assertIn("ControlMaster=auto", argv)
        self.assertTrue(any(a.endswith("/.ssh/cm-%C") for a in argv), argv)

    def test_hive_ssh_mux_0_turns_it_off(self):
        self.assertNotIn("ControlMaster=auto", self._argv(HIVE_SSH_MUX="0"))

    def test_get_does_not_copy_permissions(self):
        """Live finding A: -a from setgid dirs made macOS exit 23 after a good copy."""
        argv = self._argv("--get", "~/core/x", "./x")
        self.assertEqual(argv[0], "rsync")
        self.assertIn("-rlt", argv)
        self.assertNotIn("-a", argv)
        self.assertTrue(any("ControlMaster=auto" in a for a in argv), "the -e ssh string keeps mux")


# ---------------------------------------------------------------- partial scripts --
class TestPartialScriptsDirectory(unittest.TestCase):
    """Live finding B: a partial copy on HIVE ended in a raw traceback, exit 1."""

    def test_missing_sibling_scripts_exit_3_with_the_fix(self):
        with tempfile.TemporaryDirectory() as tmp:
            lone = os.path.join(tmp, "scripts")
            os.makedirs(lone)
            shutil.copy(CS, lone)
            summary, _ = write_summary(tmp)
            tsv = write(tmp, "sample_files.tsv", "unique_id\tcondition_name\tfile\tstatus\nKG1\tc\t/x/a.d\tmatched\n")
            session = os.path.join(tmp, "session")
            os.makedirs(session)
            for args in (["conditions", "--summary", summary, "--sample-files", tsv, "--out", os.path.join(tmp, "c.csv")],
                         ["deliver", "--summary", summary, "--session", session]):
                p = subprocess.run([sys.executable, os.path.join(lone, "core_submission.py"), *args],
                                   capture_output=True, text=True, env=child_env(tmp))
                self.assertEqual(p.returncode, 3, args[0] + p.stderr)
                self.assertIn("sync the whole scripts/ directory", json.loads(p.stdout)["hint"])


if __name__ == "__main__":
    unittest.main()
