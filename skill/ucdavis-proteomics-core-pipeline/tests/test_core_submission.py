"""
The CoreOmics submission workflow: every check here guards a failure that does NOT raise.

A service run attaches a PI's name, a sample sheet and a Bioshare folder to a set of raw
files. When any link in that chain is wrong the pipeline still succeeds -- it searches the
next lab's BN1, files the results under the wrong month, or hands the collaborator a share
full of symlinks they cannot open (PROT_0793). So the assertions are on the exit codes the
orchestrator branches on and on the filesystem state a collaborator will actually see.

Stdlib only, no network, no SSH: temp dirs, env overrides, and an in-process fake CoreOmics.
"""
import codecs
import contextlib
import csv
import datetime as dt
import http.client
import http.server
import importlib
import io
import json
import ntpath
import os
import re
import shlex
import shutil
import stat
import subprocess
import sys
import tempfile
import threading
import unittest
import urllib.error
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
    # never the Core's real staff list (it exists on HIVE): qc_bracket.py ack falls back to
    # refusing agent-like names only
    e["CORE_STAFF_FILE"] = os.path.join(tmp, "no_core_staff.txt")
    # a laptop, wherever the suite runs (on HIVE /quobyte exists): TestKeyOnHive sets "1"
    e["CORE_ON_HIVE"] = "0"
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
    """A tiny CoreOmics: routes are (method, path) -> fn(query, body) -> (status, json[, headers]);
    bytes in place of the json are sent as they are (an HTML error page)."""

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
                raw = isinstance(obj, bytes)            # a web page, not CoreOmics JSON
                data = obj if raw else json.dumps(obj).encode()
                self.send_response(status)
                for k, v in (res[2] if len(res) > 2 else {}).items():
                    self.send_header(k, v)
                self.send_header("Content-Type", "text/html" if raw else "application/json")
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
        # review of f18d95e: Python's \d matched Arabic-Indic and fullwidth digits, which int()
        # reads, so `attach --given '{"internal_id": "٧٥٦"}'` made PROT_0756
        for s in ("", "abc", "PROT_", "0", "12345678", HEX[:-1], "807a", "٧٥٦", "PROT_٠٧٥٦",
                  "１２３", "prot-８０７", "PROT\u00a0807"):
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

    def test_the_hive_root_is_the_share_tables(self):
        """Rule 3: the Flinders share's HIVE path is hive_shares.tsv's (share_map.py, which
        hive_path.sh also reads) -- core_submission keeps no second copy of it."""
        import share_map
        row = next(r for r in share_map.load_table() if r["share"] == share_map.FLINDERS_SHARE)
        self.assertEqual(cs.server_flinders_root(), row["hive"])
        with mock.patch.object(cs.share_map, "load_table", return_value=[dict(row, hive="/x/fl")]):
            self.assertEqual(cs.server_project_dir(HEX, "2025-03-10T15:02:53-07:00"),
                             f"/x/fl/coreomics/projects/2025/03/{HEX}")
        with mock.patch.object(cs.share_map, "load_table", return_value=[]):
            with self.assertRaises(cs.Stop) as e:
                cs.server_flinders_root()
            self.assertEqual(e.exception.code, cs.EXIT_UNREACHABLE)

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


class TestHtPattern(unittest.TestCase):
    """2.9.1: plate files named `<date>_PROT_<n>_...` were not recognised (0 of 96 found)."""

    def test_the_number_after_the_date_with_or_without_prot(self):
        p = cs.ht_pattern("PROT_0807")
        for name in ("20260930_PROT_0807_x.d", "20260930_PROT0807_x.d", "20260930_0807_x.d",
                     "20260930_807_x.d", "20260930_prot_0807_x.d", "20260930_PROT_807_x.d"):
            with self.subTest(name=name):
                self.assertTrue(p.match(name))

    def test_a_longer_number_is_not_the_submission(self):
        p = cs.ht_pattern("PROT_0807")
        for name in ("20260930_PROT_08070_x.d", "20260930_08070_x.d", "20260930_PROT0807A_x.d",
                     "20260930_PROT__0807_x.d", "20260930_X_0807_x.d"):
            with self.subTest(name=name):
                self.assertFalse(p.match(name))

    def test_a_run_counter_does_not_impersonate_a_submission(self):
        self.assertFalse(cs.ht_pattern("PROT_0380").match("Ex08312026_380_JE21.raw"))
        self.assertFalse(cs.ht_pattern("PROT_0380").match("Ex08312026_PROT_0380_JE21.raw"))


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

    def test_a_separator_between_letters_and_digits_is_optional(self):
        """PROT_0756: the sheet says LRS96, every run says DIA-LRS-96."""
        name = "08132026__60SPD_DIA-LRS-96_S3-B1_1_23630.d"
        self.assertTrue(cs.token_pattern("LRS96").search(cs.match_space(name)))
        self.assertTrue(cs.token_pattern("EB_001").search("x_EB001.raw"))
        self.assertTrue(cs.token_pattern("LRS-96").search("DIA-LRS96_S3"))
        # still delimited, and still no separator invented between two runs of digits
        self.assertFalse(cs.token_pattern("LRS96").search("DIA-LRS-960_S3"))
        self.assertFalse(cs.token_pattern("LRS9").search("DIA-LRS-96_S3"))
        self.assertFalse(cs.token_pattern("SG-001-2").search("x_SG0012.raw"))
        self.assertFalse(cs.token_pattern("KG1").search("x_KG-13.raw"))

    def test_prot0756_runs_each_find_their_own_sample(self):
        """All 30 PROT_0756 runs (names as on the share) match 1:1 with the sheet's ids."""
        runs = {f"LRS{n}": f"08132026__60SPD_DIA-LRS-{n}_S3-A1_1_{23600 + n}.d" for n in range(96, 126)}
        samples = [{"unique_id": u, "sample_name": u, "condition_name": "c"} for u in runs]
        entries = [{"path": "/x/" + f, "name": f, "instrument_folder": "tTOF_HT"} for f in runs.values()]
        rows = cs.match_samples(samples, entries, [], dt.date(2026, 7, 15), dt.date(2027, 3, 12))
        got = {r["unique_id"]: os.path.basename(r["file"]) for r in rows if r["status"] == "matched"}
        self.assertEqual(got, runs)

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

    def test_an_unambiguous_and_an_ambiguous_run_ask_staff_which_is_the_sample(self):
        """It used to take the unambiguous run with a WARN (ambiguous_files_excluded). Which run
        IS the sample is a question about identity, so staff answer it (2.10)."""
        files = self.clean_tree()
        later = tims(self.tmp, "03252025", "KG1", 7)
        neighbor = {"internal_id": "PROT_0810", "id": "aaaaaaaaaaaa",
                    "submitted": "2025-03-20T09:00:00-07:00", "unique_ids": ["kg1"]}
        summary, _ = write_summary(self.tmp, neighbors=[neighbor])
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        self.assertNotIn(files["KG1"], self.files_txt())
        g = gates(out)["choose_run"]
        self.assertEqual(g["status"], "FAIL")
        self.assertIn("--choose <key>=<file>", g["detail"])
        self.assertEqual(g["samples"][0]["kind"], "collision")
        (smp,) = g["samples"]
        self.assertEqual(smp["unique_id"], "KG1")
        self.assertEqual([(c["file"], c["acquired"], c["label_also_used_by"]) for c in smp["candidates"]],
                         [(later, "2025-03-25", ["PROT_0810"]), (files["KG1"], "2025-03-12", [])])
        self.assertNotIn("ambiguous_files_excluded", gates(out))
        rc, out, p = self.locate(summary, "--choose", f"KG1={files['KG1']}")
        self.assertEqual(rc, 0, p.stderr)
        self.assertIn(files["KG1"], self.files_txt())
        self.assertNotIn(later, self.files_txt())
        self.assertEqual(gates(out)["staff_choices"]["samples"][0]["chosen"], files["KG1"])
        rec = load(os.path.join(self.out, "locate.json"))["accepted"]["choices"]
        self.assertEqual(rec, [{"unique_id": "KG1", "sample_key": "KG1", "file": files["KG1"],
                                "source": "--choose"}])
        rows = {r["unique_id"]: r for r in tsv_rows(os.path.join(self.out, "sample_files.tsv"))}
        self.assertEqual(rows["KG1"]["chosen_by"], "--choose")
        self.assertIn("chosen by staff (--choose) from 2 candidate run(s)", rows["KG1"]["note"])
        self.assertEqual(rows["KG2"]["chosen_by"], "")

    def test_max_days_wider_than_the_neighbour_window_is_refused(self):
        """Defect 2: beyond fetch's window nobody checked who else used the label."""
        self.clean_tree()
        summary, _ = write_summary(self.tmp)
        self.assertEqual(self.locate(summary, "--max-days", "300")[0], 2)

    def test_several_runs_per_sample_ask_staff_which_one(self):
        """Re-injections: the most recent used to win with a WARN (alternates). Now staff say,
        and an older run is as good an answer as a newer one."""
        files = self.clean_tree()
        newer = tims(self.tmp, "03142025", "KG2", 8)
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        g = gates(out)["choose_run"]
        self.assertEqual([c["file"] for c in g["samples"][0]["candidates"]], [newer, files["KG2"]])
        self.assertEqual(g["samples"][0]["kind"], "reinjections")
        self.assertIn("ONE decision covers them all", g["detail"])
        rows = {r["unique_id"]: r for r in tsv_rows(os.path.join(self.out, "sample_files.tsv"))}
        self.assertEqual((rows["KG2"]["status"], rows["KG2"]["file"]), ("needs_choice", ""))
        self.assertIn("--choose KG2=<file>", rows["KG2"]["note"])
        self.assertNotIn("alternates", gates(out))
        choices = os.path.join(self.tmp, "choices.txt")
        with open(choices, "w") as fh:
            fh.write(f"# staff, 2026-10-01\nKG2={files['KG2']}\n")
        rc, out, p = self.locate(summary, "--choices", choices)
        self.assertEqual(rc, 0, p.stderr)
        self.assertIn(files["KG2"], self.files_txt())
        self.assertNotIn(newer, self.files_txt())
        self.assertEqual(load(os.path.join(self.out, "locate.json"))["accepted"]["choices"][0]["source"],
                         f"--choices {choices}")

    def test_a_choice_must_be_a_run_carrying_the_samples_id(self):
        files = self.clean_tree()
        tims(self.tmp, "03142025", "KG2", 8)
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.locate(summary, "--choose", f"KG2={files['KG13']}")
        self.assertEqual(rc, 2)
        self.assertIn("not one of the runs whose name carries KG2", out["error"])
        rc, out, _ = self.locate(summary, "--choose", f"KG99={files['KG1']}")
        self.assertEqual(rc, 2)
        self.assertIn("no sample of this submission has that unique_id", out["error"])
        rc, out, _ = self.locate(summary, "--choose", "KG2")
        self.assertIn("expected UNIQUE_ID=FILE", out["error"])

    def test_accepting_an_ambiguous_label_still_asks_which_of_several_runs(self):
        self.clean_tree()
        later = tims(self.tmp, "03252025", "KG1", 7)
        neighbor = {"internal_id": "PROT_0810", "id": "aaaaaaaaaaaa",
                    "submitted": "2025-02-20T09:00:00-08:00", "unique_ids": ["KG1"]}
        summary, _ = write_summary(self.tmp, neighbors=[neighbor])
        rc, out, _ = self.locate(summary, "--accept-ambiguous")
        self.assertEqual(rc, 2)
        self.assertEqual(len(gates(out)["choose_run"]["samples"][0]["candidates"]), 2)
        rc, out, _ = self.locate(summary, "--accept-ambiguous", "--choose", f"KG1={later}")
        self.assertEqual(rc, 0)
        rows = {r["unique_id"]: r for r in load(os.path.join(self.out, "locate.json"))["samples"]}
        self.assertTrue(rows["KG1"]["accepted_ambiguous"])

    def test_two_samples_with_one_label_are_asked_and_answered_by_key(self):
        """The label cannot say which of the two samples a run is -- not even with one run."""
        run1 = tims(self.tmp, "03122025", "KG1", 1)
        tims(self.tmp, "03122025", "KG2", 2)
        rec = record(samples=[("KG1", "ctrl"), ("KG2", "ctrl"), ("kg-1", "treat")])
        summary, _ = write_summary(self.tmp, rec)
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 2)
        g = gates(out)["duplicate_ids"]
        self.assertEqual(g["status"], "FAIL")
        self.assertEqual([(x["key"], x["unique_id"]) for x in g["samples"][0]],
                         [("KG1#1", "KG1"), ("kg-1#2", "kg-1")])
        self.assertNotIn(run1, self.files_txt(), "a shared label assigned its one run")
        self.assertEqual({x["key"] for x in gates(out)["choose_run"]["samples"]}, {"KG1#1", "kg-1#2"})
        rc, out, _ = self.locate(summary, "--choose", f"KG1={run1}")
        self.assertIn("2 samples share that unique_id", out["error"])
        self.assertIn("KG1#1 (sample KG1, ctrl)", out["hint"])
        # KG1#1 was run; kg-1#2 was not: its one candidate is taken, so it is unmatched
        rc, out, _ = self.locate(summary, "--choose", f"KG1#1={run1}")
        self.assertEqual(rc, 2)
        self.assertEqual(gates(out)["duplicate_ids"]["status"], "FAIL")
        rows = {r["sample_key"]: r for r in tsv_rows(os.path.join(self.out, "sample_files.tsv"))}
        self.assertEqual((rows["KG1#1"]["file"], rows["kg-1#2"]["status"]), (run1, "unmatched"))
        rc, out, p = self.locate(summary, "--choose", f"KG1#1={run1}", "--allow-partial")
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual(gates(out)["duplicate_ids"]["status"], "INFO")

    def test_a_weak_id_can_be_answered_by_a_choice(self):
        rec = record(samples=[("A5", "ctrl"), ("KG2", "ctrl")])
        mine = raw(self.tmp, "Exploris480", "mar25", "Ex03122025_400_A5.raw")
        tims(self.tmp, "03122025", "KG2", 2)
        summary, _ = write_summary(self.tmp, rec)
        self.assertEqual(self.locate(summary)[0], 2)
        rc, out, p = self.locate(summary, "--choose", f"A5={mine}")
        self.assertEqual(rc, 0, p.stderr)
        self.assertIn(mine, self.files_txt())

    def test_ht_plate_token_exits_4(self):
        self.clean_tree()
        raw(self.tmp, "tTOF_HT", "mar25", "20250315_807_100spd_Hel50_S6-A1_1_24026.d")
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 4)
        self.assertTrue(out["ht_plate"])
        self.assertIn("ht_manifest.py", out["next_step"])

    def test_an_ht_plate_named_with_prot_exits_4(self):
        """2.9.1: a plate named `<date>_PROT_<n>_...` went to sample matching and found 0 of 96."""
        self.clean_tree()
        for i in range(3):
            raw(self.tmp, "tTOF_HT", "mar25", f"20250315_PROT_0807_100spd_Hel50_S6-A{i + 1}_1_2402{i}.d")
        summary, _ = write_summary(self.tmp)
        rc, out, _ = self.locate(summary)
        self.assertEqual(rc, 4)
        self.assertTrue(out["ht_plate"])
        self.assertEqual(out["n_ht_files"], 3)

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


class TestHtFilesFrom(unittest.TestCase):
    """An HT plate's runs carry the submitter's sample_name, never the CoreOmics unique_id
    (a 96-well plate, 2026-10-01: --files-from with the plate's 96 runs left 96 samples
    unmatched, and the only bypass offered, --allow-partial, means "not run"). --files-from
    now recognises this submission's plate runs the way locate does (ht_pattern) and matches
    them on sample_name, in their sample field, by the same token rules."""

    NAMES = [("DB1", "S-101", "a"), ("DB2", "QCPOOL", "a"), ("DB3", "QCPOOL", "b"),
             ("DB4", "BX7", "b"), ("DB5", "12.5", "a"), ("DB6", "S5", "b")]

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.env = child_env(self.tmp)
        self.out = os.path.join(self.tmp, "loc")
        rec = record(samples=[(u, c) for u, _n, c in self.NAMES])
        for smp, (_u, name, _c) in zip(rec["submission_data"]["samples"], self.NAMES):
            smp["sample_name"] = name
        self.summary, _ = write_summary(self.tmp, rec)
        self.files = {}
        for i, (u, name, _c) in enumerate(self.NAMES[:5]):
            tag = "_rep2" if u == "DB3" else ""
            self.files[u] = raw(self.tmp, "tTOF_HT", "mar25",
                                f"20250315_PROT_0807_100spd_{name}{tag}_S5-A{i + 1}_1_{24100 + i}.d")
        # DB6 ("S5") was not run: every plate run has `_S5-` in its position, not its name.
        self.list = os.path.join(self.tmp, "plate.txt")
        with open(self.list, "w") as fh:
            fh.write("\n".join(self.files.values()) + "\n")

    def tearDown(self):
        self._tmp.cleanup()

    def locate(self, *extra):
        return run(["locate", "--summary", self.summary, "--out", self.out,
                    "--files-from", self.list, *extra], self.env)

    def rows(self):
        return {r["unique_id"]: r for r in tsv_rows(os.path.join(self.out, "sample_files.tsv"))
                if r["unique_id"]}

    def test_the_sample_field_is_between_the_plate_token_and_the_position(self):
        pat = cs.ht_pattern("PROT_0807")
        self.assertEqual(cs.ht_sample_field("20250315_PROT_0807_100spd_S-101_S5-A1_1_24100.d", pat),
                         "_100spd_S-101")
        self.assertEqual(cs.ht_sample_field("20250315_0807_rerun_B7_S6-H12_1_9.d", pat), "_rerun_B7")

    def test_plate_runs_match_on_sample_name_and_a_shared_name_asks_staff(self):
        rc, out, p = self.locate()
        self.assertEqual(rc, 2, p.stderr)
        rows = self.rows()
        for u in ("DB1", "DB4", "DB5"):
            self.assertEqual((rows[u]["status"], rows[u]["file"]), ("matched", self.files[u]), u)
            self.assertIn("matched on sample_name (HT plate run)", rows[u]["note"])
        self.assertEqual(rows["DB6"]["status"], "unmatched", "S5 must not match the plate position")
        g = gates(out)
        self.assertEqual(g["ht_plate_runs"]["status"], "INFO")
        self.assertEqual(g["choose_run"]["status"], "FAIL")
        self.assertEqual({x["unique_id"]: sorted(c["file"] for c in x["candidates"])
                          for x in g["choose_run"]["samples"]},
                         {u: sorted([self.files["DB2"], self.files["DB3"]]) for u in ("DB2", "DB3")})
        self.assertEqual(g["unmatched_samples"]["samples"][0]["unique_id"], "DB6")
        self.assertIn("--choose <unique_id>=<file>", g["unmatched_samples"]["detail"])

    def test_staff_choices_settle_the_shared_name(self):
        rc, out, p = self.locate("--choose", f"DB2={self.files['DB2']}",
                                 "--choose", f"DB3={self.files['DB3']}", "--allow-partial")
        self.assertEqual(rc, 0, p.stderr)
        rows = self.rows()
        self.assertEqual((rows["DB3"]["file"], rows["DB3"]["chosen_by"]), (self.files["DB3"], "--choose"))
        self.assertEqual(out["counts"]["matched"], 5)
        self.assertEqual(gates(out)["staff_choices"]["status"], "INFO")

    def test_a_well_like_label_on_a_plate_run_is_asked_not_taken(self):
        """A submitter's tube labelled like a well (B3): a run carrying it in its sample field may
        still be naming a well, so staff say whose run it is (2.10 review, HIGH 3)."""
        rec = record(samples=[("WL1", "a"), ("WL2", "b")])
        for smp, name in zip(rec["submission_data"]["samples"], ("B3", "Kidney")):
            smp["sample_name"] = name
        summary, _ = write_summary(self.tmp, rec, name="wells.json")
        b3 = raw(self.tmp, "tTOF_HT", "mar25", "20250315_PROT_0807_100spd_B3_S5-C1_1_24400.d")
        kid = raw(self.tmp, "tTOF_HT", "mar25", "20250315_PROT_0807_100spd_Kidney_S5-B3_1_24401.d")
        lst = os.path.join(self.tmp, "wells.txt")
        with open(lst, "w") as fh:
            fh.write(f"{b3}\n{kid}\n")
        rc, out, _ = run(["locate", "--summary", summary, "--out", self.out, "--files-from", lst],
                         self.env)
        self.assertEqual(rc, 2)
        need = {x["unique_id"]: x for x in gates(out)["choose_run"]["samples"]}
        self.assertEqual(need["WL1"]["kind"], "well_like_label")
        self.assertEqual([c["file"] for c in need["WL1"]["candidates"]], [b3])
        rows = {r["unique_id"]: r for r in tsv_rows(os.path.join(self.out, "sample_files.tsv"))
                if r["unique_id"]}
        # Kidney's run sits in the Core's well B3; that never makes it WL1's
        self.assertEqual((rows["WL2"]["status"], rows["WL2"]["file"]), ("matched", kid))
        rc, out, p = run(["locate", "--summary", summary, "--out", self.out, "--files-from", lst,
                          "--choose", f"WL1={b3}"], self.env)
        self.assertEqual(rc, 0, p.stderr)

    def test_a_choice_must_be_in_the_list(self):
        other = raw(self.tmp, "tTOF_HT", "mar25", "20250315_PROT_0807_100spd_QCPOOL_S5-B1_1_24200.d")
        rc, out, _ = self.locate("--choose", f"DB2={other}")
        self.assertEqual(rc, 2)
        self.assertIn("not in the --files-from list", out["error"])

    def test_another_submissions_plate_is_not_matched_on_names(self):
        stray = raw(self.tmp, "tTOF_HT", "mar25", "20250315_PROT_9999_100spd_S-101_S5-C1_1_24300.d")
        with open(self.list, "a") as fh:
            fh.write(stray + "\n")
        rc, out, _ = self.locate()
        loose = {r["file"]: r for r in tsv_rows(os.path.join(self.out, "sample_files.tsv"))
                 if r["status"] == "unassigned"}
        self.assertIn("matches no sample id or sample name", loose[stray]["note"])


class TestReinjections(unittest.TestCase):
    """Several runs of ONE sample are re-injections, and one staff decision covers all of them
    (team lead for Brett, 2026-10-01: a 96-well plate can have many): --reinjections
    ask|latest|all. A label two samples or submissions share is a collision and still needs its
    own --choose, whatever the policy."""

    IDS = [f"RJ{i}" for i in range(1, 10)]

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.env = child_env(self.tmp)
        self.out = os.path.join(self.tmp, "loc")
        rec = record(samples=[(u, "ctrl" if i < 5 else "treat") for i, u in enumerate(self.IDS)])
        self.summary, _ = write_summary(self.tmp, rec)
        # every sample injected twice the same day: the counter orders them
        self.first = {u: tims(self.tmp, "03122025", u, 10 + i) for i, u in enumerate(self.IDS)}
        self.rerun = {u: tims(self.tmp, "03122025", u, 50 + i) for i, u in enumerate(self.IDS)}

    def tearDown(self):
        self._tmp.cleanup()

    def locate(self, *extra):
        return run(["locate", "--summary", self.summary, "--out", self.out, *extra], self.env)

    def rows(self):
        return tsv_rows(os.path.join(self.out, "sample_files.tsv"))

    def files_txt(self):
        return read(os.path.join(self.out, "files.txt")).split()

    def test_ask_is_the_default_and_names_the_one_decision(self):
        rc, out, _ = self.locate()
        self.assertEqual(rc, 2)
        g = gates(out)["choose_run"]
        self.assertEqual(len(g["samples"]), 9)
        self.assertEqual({x["kind"] for x in g["samples"]}, {"reinjections"})
        self.assertIn("--reinjections latest", g["detail"])

    def test_latest_settles_nine_reinjections_in_one_call(self):
        rc, out, p = self.locate("--reinjections", "latest")
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual(sorted(self.files_txt()), sorted(self.rerun.values()))
        g = gates(out)["reinjections"]
        self.assertEqual(g["status"], "INFO")
        self.assertIn("--reinjections latest: 9 sample(s)", g["detail"])
        self.assertIn("9 run(s) set aside", g["detail"])
        self.assertEqual({x["key"] for x in g["samples"]}, set(self.IDS))
        rows = {r["unique_id"]: r for r in self.rows()}
        self.assertEqual(rows["RJ3"]["reinjections"],
                         f"latest (set aside: {os.path.basename(self.first['RJ3'])})")
        self.assertIn("kept " + os.path.basename(self.rerun["RJ3"]), rows["RJ3"]["note"])
        loc = load(os.path.join(self.out, "locate.json"))
        self.assertEqual(loc["accepted"]["reinjections"], "latest")
        rec = {r["unique_id"]: r for r in loc["samples"]}["RJ3"]["reinjections"]
        self.assertEqual((rec["kept"], rec["set_aside"]), ([self.rerun["RJ3"]], [self.first["RJ3"]]))

    def test_all_keeps_every_injection_and_says_so(self):
        rc, out, p = self.locate("--reinjections", "all")
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual(len(self.files_txt()), 18)
        g = gates(out)["reinjections"]
        self.assertIn("technical replicate (18 runs)", g["detail"])
        rows = [r for r in self.rows() if r["unique_id"] == "RJ1"]
        self.assertEqual(len(rows), 2)
        self.assertTrue(all(r["reinjections"] == "all (2 runs kept as technical replicates)" for r in rows))
        # the DE never counts them as independent samples: conditions.csv names each run's
        # sample (Sample), which run_de.R blocks on (2.10 review, MED 7)
        out_csv = os.path.join(self.tmp, "conditions.csv")
        rc, js, p = run(["conditions", "--summary", self.summary, "--sample-files",
                         os.path.join(self.out, "sample_files.tsv"), "--out", out_csv], self.env)
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        self.assertEqual(len(js["technical_replicates"]["samples"]), 9)
        with open(out_csv, newline="") as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual({r["Sample"] for r in rows}, set(self.IDS))
        with open(out_csv + ".decisions.json") as fh:
            self.assertEqual(json.load(fh)["technical_replicates"]["column"], "Sample")

    def test_a_collision_blocks_under_every_policy(self):
        rec = record(samples=[("KG1", "ctrl"), ("kg-1", "treat"), ("KG2", "ctrl")])
        summary, _ = write_summary(self.tmp, rec, name="collide.json")
        tims(self.tmp, "03122025", "KG1", 1)
        tims(self.tmp, "03122025", "KG1", 2)
        tims(self.tmp, "03122025", "KG2", 3)
        for policy in ("latest", "all"):
            rc, out, _ = run(["locate", "--summary", summary, "--out", self.out,
                              "--reinjections", policy], self.env)
            self.assertEqual(rc, 2, policy)
            kinds = {x["key"]: x["kind"] for x in gates(out)["choose_run"]["samples"]}
            self.assertEqual(kinds, {"KG1#1": "collision", "kg-1#2": "collision"}, policy)
        # another submission's label on one of the runs is a collision too, unless accepted
        neighbor = {"internal_id": "PROT_0810", "id": "aaaaaaaaaaaa",
                    "submitted": "2025-03-11T09:00:00-07:00", "unique_ids": ["RJ4"]}
        summary2, _ = write_summary(self.tmp, record(samples=[(u, "c") for u in self.IDS]),
                                    neighbors=[neighbor], name="neighbour.json")
        rc, out, _ = run(["locate", "--summary", summary2, "--out", self.out,
                          "--reinjections", "latest"], self.env)
        self.assertEqual(rc, 2)
        self.assertEqual([x["key"] for x in gates(out)["choose_run"]["samples"]
                          if x["kind"] == "collision"], [])
        self.assertEqual(gates(out)["ambiguous_label"]["status"], "FAIL")

    def test_latest_never_orders_same_day_runs_across_instruments(self):
        """2.10 review: the acquisition counter counts ONE instrument's runs. Same-day runs on two
        instruments (or with no counter) cannot be ordered, so latest asks for that sample."""
        other = raw(self.tmp, "Exploris480", "mar25", "Ex03122025_400_RJ1.raw")
        rc, out, _ = self.locate("--reinjections", "latest")
        self.assertEqual(rc, 2)
        (need,) = gates(out)["choose_run"]["samples"]
        self.assertEqual((need["key"], need["kind"]), ("RJ1", "reinjections"))
        self.assertEqual(sorted(c["file"] for c in need["candidates"]),
                         sorted([self.first["RJ1"], self.rerun["RJ1"], other]))
        row = {r["unique_id"]: r for r in self.rows()}["RJ1"]
        self.assertIn("different instruments (Exploris480, tTOF_HT)", row["note"])
        # every other sample was still decided by the one policy
        self.assertEqual(len(gates(out)["reinjections"]["samples"]), 8)
        # same day, one instrument, no counter in the names: also asked
        for n in (401, 402):
            raw(self.tmp, "Exploris480", "mar25", f"Ex03122025_{n}_RJ2.raw")
        rc, out, _ = self.locate("--reinjections", "latest")
        keys = {x["key"] for x in gates(out)["choose_run"]["samples"]}
        self.assertIn("RJ2", keys)

    def test_dates_order_runs_across_instruments(self):
        newer = raw(self.tmp, "Exploris480", "mar25", "Ex03142025_400_RJ1.raw")
        rc, out, p = self.locate("--reinjections", "latest")
        self.assertEqual(rc, 0, p.stderr)
        self.assertIn(newer, self.files_txt())

    def test_latest_never_keeps_a_run_stan_flags_needs_rerun(self):
        stan = os.path.join(self.tmp, "ht_manifest.json")
        with open(stan, "w") as fh:
            json.dump({"entries": [{"raw_path": self.rerun["RJ5"], "needs_rerun": True},
                                   {"raw_path": self.first["RJ5"], "needs_rerun": False}]}, fh)
        rc, out, p = self.locate("--reinjections", "latest", "--ht-manifest", stan)
        self.assertEqual(rc, 0, p.stderr)
        rec = {r["unique_id"]: r for r in load(os.path.join(self.out, "locate.json"))["samples"]}["RJ5"]
        self.assertEqual(rec["file"], self.first["RJ5"])
        self.assertIn("needs_rerun", rec["reinjections"]["reason"])

    def test_all_with_one_name_in_two_folders_keeps_both_and_flags_them(self):
        twin = tims(self.tmp, "03122025", "RJ1", 10, month="apr25")    # same name, another folder
        rc, out, p = self.locate("--reinjections", "all")
        self.assertEqual(rc, 0, p.stderr)
        self.assertIn(twin, self.files_txt())
        g = gates(out)["repeated_names"]
        self.assertEqual(g["status"], "WARN")
        self.assertIn("run_search.py stops before searching them as given", g["detail"])
        self.assertEqual(sorted(g["names"][0]["paths"]), sorted([self.first["RJ1"], twin]))

    def test_curated_lists_follow_the_same_policy(self):
        lst = os.path.join(self.tmp, "list.txt")
        with open(lst, "w") as fh:
            fh.write("\n".join(list(self.first.values()) + list(self.rerun.values())) + "\n")
        rc, out, _ = self.locate("--files-from", lst)
        self.assertEqual(rc, 2)
        self.assertEqual(len(gates(out)["choose_run"]["samples"]), 9)
        self.assertEqual(self.files_txt(), [])
        rc, out, p = self.locate("--files-from", lst, "--reinjections", "latest")
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual(sorted(self.files_txt()), sorted(self.rerun.values()))
        rc, out, p = self.locate("--files-from", lst, "--reinjections", "all")
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual(len(self.files_txt()), 18)

    def test_the_analysis_report_names_every_sample_the_policy_applied_to(self):
        """attach carries locate's decisions into the session; the report's Data Quality Notes
        (submission_report.quality_notes) name each sample, file names only."""
        rc, _out, p = run(["locate", "--summary", self.summary, "--out", self.tmp,
                           "--reinjections", "latest"], self.env)
        self.assertEqual(rc, 0, p.stderr)
        sess = os.path.join(self.tmp, "session")
        os.makedirs(os.path.join(sess, "input"))
        r = subprocess.run([sys.executable, os.path.join(SCRIPTS, "submission_report.py"), "attach",
                            "--session", sess, "--record", self.tmp], capture_output=True, text=True)
        self.assertEqual(r.returncode, 0, r.stderr + r.stdout)
        dec = load(os.path.join(sess, "input", "locate_decisions.json"))
        self.assertEqual(dec["reinjections_policy"], "latest")
        self.assertEqual(dec["reinjections"]["RJ7"]["set_aside"], [os.path.basename(self.first["RJ7"])])
        sys.path.insert(0, SCRIPTS)
        import submission_report as sr
        rec = sr.load(sess)
        notes = {n["id"]: n["text"] for n in sr.quality_notes(rec, sess)}
        self.assertIn("9 sample(s) had several runs of their own, and staff chose the newest "
                      "injection of each (--reinjections latest)", notes["reinjections"])
        self.assertIn(f"RJ7: kept {os.path.basename(self.rerun['RJ7'])}, set aside "
                      f"{os.path.basename(self.first['RJ7'])}", notes["reinjections"])
        self.assertNotIn(self.tmp, notes["reinjections"], "no paths in a collaborator's report")

    def test_stage_records_the_policy_and_each_sample(self):
        sroot = os.path.join(self.tmp, "flinders", "Data", "lab", "service")
        group = os.path.join(sroot, "on_campus", "Placeholder lab")
        os.makedirs(group)
        os.makedirs(os.path.join(sroot, "off_campus"))
        self.assertEqual(self.locate("--reinjections", "latest")[0], 0)
        rc, out, p = run(["stage", "--summary", self.summary, "--files",
                          os.path.join(self.out, "files.txt"), "--apply"], self.env)
        self.assertEqual(rc, 0, p.stderr)
        marker = load(os.path.join(group, "PROT_0807", cs.MARKER))
        self.assertEqual(marker["locate_accepted"]["reinjections"], "latest")
        self.assertEqual(len(marker["reinjections"]), 9)
        md = read(os.path.join(group, "PROT_0807", "SUBMISSION.md"))
        self.assertIn(f"| matched (re-injections: latest (set aside: "
                      f"{os.path.basename(self.first['RJ2'])})) |", md)


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
        self.assertNotIn("@", read(os.path.join(project, "SUBMISSION.md")), "no email on HIVE")
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

    def test_apply_records_the_staff_choices_the_file_list_rests_on(self):
        group = self.folder("on_campus", "Placeholder lab")
        summary, _ = write_summary(self.tmp)
        choice = {"unique_id": "KG2", "file": self.files[1], "source": "--choose"}
        with open(os.path.join(self.tmp, "locate.json"), "w") as fh:
            json.dump({"mode": "auto", "hard_fail": False, "gates": [],
                       "accepted": {"allow_partial": False, "accept_ambiguous": False,
                                    "choices": [choice]}}, fh)
        with open(os.path.join(self.tmp, "sample_files.tsv"), "w") as fh:
            fh.write("unique_id\tsample_name\tcondition_name\tfile\tacquired\tstatus\tchosen_by\n"
                     f"KG1\ts1\tctrl\t{self.files[0]}\t2025-03-12\tmatched\t\n"
                     f"KG2\ts2\tctrl\t{self.files[1]}\t2025-03-12\tmatched\t--choose\n")
        rc, out, p = self.stage(summary, "--apply")
        self.assertEqual(rc, 0, p.stderr)
        project = os.path.join(group, "PROT_0807")
        self.assertEqual(load(os.path.join(project, cs.MARKER))["locate_accepted"]["choices"], [choice])
        md = read(os.path.join(project, "SUBMISSION.md"))
        self.assertIn("| matched (staff choice: --choose) |", md)

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

    def conditions(self, pairs, status="matched", extra=()):
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
        rc, js, p = run(["conditions", "--summary", self.summary, "--sample-files", tsv, "--out", out,
                         *extra], self.env)
        return rc, js, out, p

    def groups_in(self, out):
        with open(out, newline="") as fh:
            return dict(Counter_(r["Group"] for r in csv.DictReader(fh)))

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

    # Integration 2.10: a condition per REPLICATE (the Core LIMS: <sample>_mix_1 .. _n) is
    # read by collect_conditions.py as conditions + replicates, and its proposal is relayed.
    PER_REPLICATE = [("KG1", "CtrlA_mix_1"), ("KG2", "CtrlA_mix_2"), ("KG13", "TrtB_mix_1"),
                     ("KG14", "TrtB_mix_2")]

    def test_per_replicate_conditions_are_proposed_as_conditions_and_asked(self):
        rc, js, out, p = self.conditions(self.PER_REPLICATE)
        self.assertEqual(rc, 2, p.stderr + p.stdout)
        self.assertEqual(js["findings"]["replicate_labels"]["conditions"],
                         {"CtrlA_mix": {"1": "CtrlA_mix_1", "2": "CtrlA_mix_2"},
                          "TrtB_mix": {"1": "TrtB_mix_1", "2": "TrtB_mix_2"}})
        self.assertNotIn("all_unique", js["findings"])          # one question, not two
        self.assertEqual(js["groups"], {"CtrlA_mix": 2, "TrtB_mix": 2})
        self.assertTrue(any("--replicate-labels collapse" in q for q in js["questions"]))
        self.assertEqual(self.groups_in(out), {"CtrlA_mix": 2, "TrtB_mix": 2})

    def test_the_users_collapse_answer_needs_no_more_input(self):
        rc, js, out, p = self.conditions(self.PER_REPLICATE, extra=("--replicate-labels", "collapse"))
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        self.assertFalse(js["needs_user_input"])
        self.assertEqual(self.groups_in(out), {"CtrlA_mix": 2, "TrtB_mix": 2})

    def test_an_ambiguous_per_replicate_column_is_asked_with_its_reasons(self):
        rc, js, out, p = self.conditions([("KG1", "A_1"), ("KG2", "A_3"), ("KG13", "B_1"),
                                          ("KG14", "B_2")])
        self.assertEqual(rc, 2, p.stderr + p.stdout)
        amb = js["findings"]["replicate_labels_ambiguous"]
        self.assertTrue(any("do not run" in r for r in amb["reasons"]), amb)
        self.assertNotIn("all_unique", js["findings"])
        q = [q for q in js["questions"] if "--replicate-labels" in q]
        self.assertEqual(len(q), 1)
        self.assertIn("do not run", q[0])
        self.assertEqual(self.groups_in(out), {"A_1": 1, "A_3": 1, "B_1": 1, "B_2": 1})

    def test_keep_uses_the_labels_as_given(self):
        rc, js, out, p = self.conditions(self.PER_REPLICATE, extra=("--replicate-labels", "keep"))
        self.assertEqual(rc, 2)
        self.assertTrue(js["findings"]["all_unique"])
        self.assertEqual(len(self.groups_in(out)), 4)

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
NORM_DECIDED = {"normalization_check": {
    "schema": "normalization_check", "schema_version": 1, "status": "decided",
    "decision": {"quantities": "normalised", "by": "the skill: the default for this experiment "
                                                  "type, and the data check agreed"}}}


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
        # qc_bracket.py's record (step 8e), verdict good: a check or concern would hold delivery
        write(self.session, "logs/qc_bracket.json",
              json.dumps(qc_record("good", session=self.session)))
        write(o, "Analysis_Report.html", "<html>report</html>")
        write(o, "AUDIT.md", "# audit")
        # The DE table is a symlink to a file elsewhere: the delivery must hold the BYTES.
        elsewhere = write(self.tmp, "elsewhere/DE_dpc_ctrl.vs.treat.csv", "Protein,logFC\nP1,1\n")
        os.symlink(elsewhere, os.path.join(o, "tables", "DE_dpc_ctrl.vs.treat.csv"))
        write(o, "tables/reproducibility_log.R", "# R")
        # run_de.R's record of a DE run under a decided normalisation check (step 8): without one
        # the delivery is held (normalization_check.delivery_gate)
        write(o, "tables/de_provenance.json", json.dumps(NORM_DECIDED))
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


def qc_record(*verdicts, status="ok", session=None):
    """A minimal qc_bracket.py record (schema 1): one instrument per verdict, covering `session`'s
    raw-file list as it is now (the gate refuses a record for another list)."""
    import qc_bracket
    insts = [{"instrument": f"Instrument {i}", "verdict": v} for i, v in enumerate(verdicts, 1)]
    order = {"good": 0, "check": 1, "no_qc_on_record": 2, "concern": 3}
    files = qc_bracket.session_files(session)[0] if session else []
    return {"staff_only_marker": qc_bracket.STAFF_ONLY_MARKER,      # first, as qc_bracket writes it
            "schema": "qc_bracket/1", "schema_version": 1, "status": status,
            "verdict": max(verdicts, key=order.get) if verdicts and status == "ok" else None,
            "instruments": insts, "stan": {"url": "http://stan.invalid", "error": None},
            "files": {"n": len(files), "list_sha256": qc_bracket.list_sha256(files)},
            "not_checked": []}


class TestDeliverQcGate(DeliverBase):
    """Brett, 2026-10-01: the QC verdict is staff-only, and a check or concern holds the client
    delivery until a staff member records an acknowledgement (who, when, a one-line note)."""

    def set_qc(self, *verdicts, **kw):
        write(self.session, "logs/qc_bracket.json",
              json.dumps(qc_record(*verdicts, session=self.session, **kw)))

    def ack(self, note="reviewed the 2 runs after the failed QC; data stand"):
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "qc_bracket.py"), "ack",
                            "--session", self.session, "--by", "staffexample", "--note", note],
                           capture_output=True, text=True, env=self.env)
        self.assertEqual(p.returncode, 0, p.stdout + p.stderr)

    def test_check_or_concern_holds_the_delivery_dry_run_included(self):
        for verdicts in (("check",), ("good", "concern")):
            self.set_qc(*verdicts)
            rc, out, _ = self.deliver()
            self.assertEqual(rc, 2, verdicts)
            self.assertIn("held for staff review", out["error"])
            self.assertIn("qc_bracket.py ack", out["hint"])
            rc, out, _ = self.deliver("--apply")
            self.assertEqual(rc, 2, verdicts)
            self.assertFalse(os.path.exists(self.delivery))

    def test_an_acknowledgement_lets_it_through_and_a_new_record_needs_a_new_one(self):
        self.set_qc("concern")
        self.ack()
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        # delivery.json: whether and when, never who or why (2.10 safety review)
        self.assertEqual(set(out["qc_gate"]), {"proceed", "record_sha256", "ack_at"})
        with open(os.path.join(self.session, "logs", "qc_gate.json")) as fh:
            self.assertEqual(json.load(fh)["gate"]["ack"]["by"], "staffexample")
        self.set_qc("check")                          # a re-run: another record
        rc, out, _ = self.deliver("--label", "t2")
        self.assertEqual(rc, 2)

    def test_good_and_no_qc_on_record_proceed_and_no_qc_is_noted_staff_side(self):
        self.set_qc("good", "no_qc_on_record")
        rc, out, p = self.deliver()
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertEqual(out["qc_gate"]["no_qc_on_record"], ["Instrument 2"])   # printed to staff
        self.assertFalse([w for w in out["warnings"] if "QC" in w])          # not in delivery.json

    def test_no_record_or_a_check_that_could_not_run_needs_an_acknowledgement(self):
        os.remove(os.path.join(self.session, "logs", "qc_bracket.json"))
        rc, out, _ = self.deliver()
        self.assertEqual(rc, 2)
        self.assertIn("were not checked", out["error"])
        self.set_qc(status="stan_unreachable")
        rc, out, _ = self.deliver()
        self.assertEqual(rc, 2)
        self.assertIn("could not be checked", out["error"])

    def test_a_record_for_another_file_list_is_refused_even_acknowledged(self):
        self.set_qc("check")
        self.ack()
        write(self.session, "input/raw_files.txt", "# Raw MS files\n" + "\n".join(self.raw_files) + "\n")
        rc, out, _ = self.deliver()
        self.assertEqual(rc, 2)
        self.assertIn("re-run step 8e", out["error"])
        self.assertTrue(out["qc_gate"]["needs_rerun"])

    def test_delivery_json_carries_no_verdict_and_is_not_world_readable(self):
        """2.10 safety review: deliver wrote the gate -- verdict, the acknowledging login, the
        note -- into <session>/delivery.json at 0664, in a service tree others can read."""
        old = os.umask(0o002)                                # HIVE's umask
        try:
            self.set_qc("concern")
            self.ack(note="CANARYNOTE reviewed")
            rc, out, p = self.deliver("--apply")
            self.assertEqual(rc, 0, p.stdout + p.stderr)
        finally:
            os.umask(old)
        dj = os.path.join(self.session, "delivery.json")
        text = read(dj)
        for canary in ('"verdict"', "concern", "CANARYNOTE", "staffexample", "reason"):
            self.assertNotIn(canary, text.split('"qc_gate"')[1].split("}")[0], canary)
        self.assertNotIn("CANARYNOTE", text)
        self.assertEqual(os.stat(dj).st_mode & 0o777, 0o640)

    def test_a_qc_record_under_another_name_is_never_delivered(self):
        self.set_qc("concern")
        self.ack()
        write(self.output, "figures/notes.json", read(os.path.join(self.session, "logs",
                                                                   "qc_bracket.json")))
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.delivery, "figures", "notes.json")))
        self.assertIn("staff-only", read(os.path.join(self.delivery, "MANIFEST.txt")))

    def test_the_qc_record_is_never_delivered(self):
        self.set_qc("concern")
        self.ack()
        write(self.output, "tables/qc_bracket.json", "{}")       # a stray copy, wherever it lands
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stdout + p.stderr)
        for dp, _dn, fns in os.walk(self.delivery):
            self.assertFalse([f for f in fns if f.startswith("qc_bracket")], dp)


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
        self.assertNotIn(self.tmp, readme.split("## Where this lives on HIVE")[0],
                         "no local paths outside the HIVE table in a collaborator README")
        self.assertFalse(os.path.exists(os.path.join(self.share, "raw")))
        dj = load(os.path.join(self.session, "delivery.json"))
        self.assertTrue(dj["verified"])
        self.assertEqual((dj["share_dir"], dj["internal_id"]), (SERVER_SHARE, "PROT_0807"))

    def test_the_podcast_audio_and_transcript_go_to_the_collaborator(self):
        rc, out, p = self.deliver("--apply", "--label", "t0")
        self.assertEqual(rc, 0, p.stderr)
        manifest_no = read(os.path.join(self.share, "PROT_0807_t0", "MANIFEST.txt"))
        self.assertIn("[SKIPPED] podcast/ -- none was made (optional)", manifest_no)
        for fn in ("podcast.m4a", "transcript.html", "podcast.json", "podcast_script.md", ".cache/a.wav"):
            write(self.output, "podcast/" + fn, fn)
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        got = sorted(os.listdir(os.path.join(self.delivery, "podcast")))
        self.assertEqual(got, ["podcast.m4a", "transcript.html"], "no cache, script or consent record")
        self.assertIn("`podcast/`", read(os.path.join(self.delivery, "README.md")))

    def test_the_report_with_the_audio_built_in_goes_only_while_current(self):
        # make_podcast.py share: the report with the audio and transcript inside, one file the
        # collaborator can send on. A copy built from another report or episode stays behind.
        mp = cs.make_podcast
        write(self.output, "podcast/podcast.json", json.dumps(
            {"show": "Signal to Noise", "title": "T", "audio": "podcast.m4a",
             "transcript": "transcript.html", "duration_s": 60}))
        for fn in ("podcast.m4a", "transcript.html"):
            write(self.output, "podcast/" + fn, fn)
        man = json.loads(read(os.path.join(self.output, "podcast", "podcast.json")))
        write(self.output, mp.SHARE_NAME, f"<html><!-- podcast:share "
                                          f"{mp.share_keys(self.output, man)} kbps=48 unchecked=1 --></html>")
        rc, out, p = self.deliver("--apply", "--label", "t2")
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        d2 = os.path.join(self.share, "PROT_0807_t2")
        self.assertTrue(os.path.isfile(os.path.join(d2, mp.SHARE_NAME)))
        self.assertIn(f"[OK] {mp.SHARE_NAME}", read(os.path.join(d2, "MANIFEST.txt")))
        readme = read(os.path.join(d2, "README.md"))
        self.assertIn(f"**To send the report on with its audio discussion, send `{mp.SHARE_NAME}`**",
                      readme)
        self.assertIn(f"| [`{mp.SHARE_NAME}`]({mp.SHARE_NAME}) | the same report with that audio", readme)
        out_md = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        rc, js, p = run(["email-draft", "--summary", self.summary, "--delivery",
                         os.path.join(self.session, "delivery.json"), "--out", out_md], self.env)
        self.assertIn(f"{mp.SHARE_NAME}: the report with an AI-generated audio discussion built in "
                      "-- to pass the report on with the audio, send this one file", read(out_md))

        write(self.output, "Analysis_Report.html", "<html>report, regenerated</html>")
        rc, out, p = self.deliver("--apply", "--label", "t3")
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        d3 = os.path.join(self.share, "PROT_0807_t3")
        self.assertFalse(os.path.exists(os.path.join(d3, mp.SHARE_NAME)))
        self.assertIn(f"[SKIPPED] {mp.SHARE_NAME} -- built from another version of the report or "
                      "the episode; `make_podcast.py share output/` builds it (finalize does)",
                      read(os.path.join(d3, "MANIFEST.txt")))
        self.assertNotIn(mp.SHARE_NAME, read(os.path.join(d3, "README.md")))

    def test_the_readme_and_the_email_ask_how_we_did(self):
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        line = ("**How did we do?** Tell us what you thought of the report: [a 5-minute survey]"
                f"({SURVEY}?prot=PROT_0807&src=readme).")
        readme = read(os.path.join(self.delivery, "README.md"))
        self.assertIn("Contact the UC Davis Proteomics Core.\n\n" + line + "\n", readme)
        self.assertIn(f'<a href="{SURVEY}?prot=PROT_0807&amp;src=readme">a 5-minute survey</a>',
                      read(os.path.join(self.delivery, "README.html")))
        out_md = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        run(["email-draft", "--summary", self.summary, "--delivery",
             os.path.join(self.session, "delivery.json"), "--out", out_md], self.env)
        self.assertIn("Let us know if anything is unclear.\n\nHow did we do? Tell us what you "
                      f"thought of the report: a 5-minute survey at {SURVEY}?prot=PROT_0807"
                      "&src=email\n\nBest regards,", read(out_md))

    def test_the_review_records_stay_in_the_session(self):
        # the analysis conversation and the decisions log (save_transcript.py, log_decision.py)
        # are Core-internal: not delivered, and the collaborator's AGENTS.md does not send their
        # AI looking for them (review of 32a0780)
        write(self.session, "logs/conversation/11111111-2222-4333-8444-555555555555.jsonl", "{}")
        write(self.session, "logs/conversation/conversation.md", "# c")
        write(self.session, "logs/decisions.md", "# Decisions log")
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        agents = read(os.path.join(self.delivery, "AGENTS.md"))
        for word in ("Reviewing this analysis", "logs/conversation", "conversation.md",
                     "decisions.md", "transcript"):
            self.assertNotIn(word, agents)
        got = [os.path.relpath(os.path.join(dp, f), self.delivery)
               for dp, _, fs in os.walk(self.delivery) for f in fs]
        self.assertFalse([g for g in got if "conversation" in g or "decisions" in g], got)

    def test_readme_html_and_agents_md_go_to_the_collaborator(self):
        """Brett's standing rule: every collaborator folder has README.html to double-click and
        AGENTS.md for an AI assistant -- with links that work IN THE DELIVERY'S layout."""
        rc, out, p = self.deliver("--apply")
        self.assertEqual(rc, 0, p.stderr + p.stdout)
        page = read(os.path.join(self.delivery, "README.html"))
        hrefs = re.findall(r'<a href="([^"#:]+)"', page)
        self.assertIn("Analysis_Report.html", hrefs)
        for h in hrefs:
            self.assertTrue(os.path.exists(os.path.join(self.delivery, urllib.parse.unquote(h))), h)
        self.assertNotIn("Open `README.html`", page, "the HTML does not point at itself")
        readme = read(os.path.join(self.delivery, "README.md"))
        self.assertIn("**Open `README.html`**", readme)
        self.assertIn("## Where this lives on HIVE", readme)
        self.assertIn("[`AGENTS.md`](AGENTS.md)", readme)
        agents = read(os.path.join(self.delivery, "AGENTS.md"))
        self.assertTrue(agents.startswith("# AGENTS.md"), "the delivery note sits under the title")
        self.assertIn("delivery copy of an analysis session", agents)
        self.assertIn("## Where this lives on HIVE", agents)
        for doc in (page, readme, agents):
            self.assertIsNone(re.search(r"[\w.+-]+@[\w-]+\.[\w.-]+", doc), "an email address leaked")
        manifest = read(os.path.join(self.delivery, "MANIFEST.txt"))
        for name in ("AGENTS.md", "README.html", "README.md"):
            self.assertIn(f"[OK] {name}", manifest)
        sums = read(os.path.join(self.delivery, "checksums.sha256"))
        self.assertIn("  README.html", sums)
        self.assertIn("  AGENTS.md", sums)
        self.assertTrue(load(os.path.join(self.session, "delivery.json"))["verified"])

    def test_an_unreadable_session_skips_agents_md_and_says_why(self):
        docs = importlib.import_module("session_docs")
        with tempfile.TemporaryDirectory() as out_dir, \
                mock.patch.object(docs, "gather", side_effect=ValueError("bad record")):
            delivered = {"Analysis_Report.html"}
            lines = cs.delivery_docs(self.s, self.session, out_dir, delivered, [], "analysis",
                                     out_dir, [])
            self.assertEqual(dict(lines)["AGENTS.md"],
                             "the session's records could not be read: ValueError: bad record")
            self.assertIsNone(dict(lines)["README.html"])
            self.assertFalse(os.path.exists(os.path.join(out_dir, "AGENTS.md")))
            self.assertTrue(os.path.isfile(os.path.join(out_dir, "README.md")))

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
        # its QC check covers its own file list (the gate refuses a record for another one)
        write(other, "logs/qc_bracket.json", json.dumps(qc_record("good", session=other)))
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
        # README.html even here; AGENTS.md needs a session to describe, and says so
        self.assertTrue(os.path.isfile(os.path.join(delivery, "README.html")))
        self.assertIn("[`README.html`](README.html)", readme)
        manifest = read(os.path.join(delivery, "MANIFEST.txt"))
        self.assertIn("[OK] README.html", manifest)
        self.assertIn("[SKIPPED] AGENTS.md -- no analysis session to describe (a raw-only delivery)",
                      manifest)
        # no report was made, so nothing asks what they thought of one
        self.assertNotIn(SURVEY, readme)
        out_md = os.path.join(self.tmp, "EMAIL_DRAFT.md")
        run(["email-draft", "--summary", self.summary, "--delivery",
             os.path.join(self.work, "delivery.json"), "--out", out_md], self.env)
        self.assertNotIn(SURVEY, read(out_md))


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

    def test_only_redacted_files_are_written_for_hive(self):
        """SKILL 1c --puts hive/ to HIVE: never an email, contact or billing field."""
        rec = record()
        rec["payment"] = {"display": {"PPMS Order Ref #": "PPMS-CANARY"}}
        rec["submission_data"]["description"] = "call 530-555-0100 or ada.example@example.org"
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, rec)
            out_dir = os.path.join(tmp, "PROT_0807")
            rc, out, p = run(["fetch", "807", "--out", out_dir], child_env(tmp, fake.base))
            self.assertEqual(rc, 0, p.stderr)
            full = load(os.path.join(out_dir, "submission_summary.json"))
            hive_s = load(out["outputs"]["hive_summary"])
            for f in ("hive_summary", "hive_record"):
                text = read(out["outputs"][f])
                for bad in ("@", "555-0100", "PPMS-CANARY", "Robin"):
                    self.assertNotIn(bad, text, f"{bad} in {f}")
            self.assertEqual(hive_s["share_dir"], full["share_dir"])
            self.assertEqual(hive_s["samples"], full["samples"])
            self.assertEqual(load(out["outputs"]["hive_record"])["schema"], "submission_record/1")
            # and the HIVE-side steps still work from it
            s = cs.load_summary(out["outputs"]["hive_summary"])
            self.assertEqual(s["neighbor_window_days"], full["neighbor_window_days"])

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


SURVEY = "https://feedback-ucd-proteomics.azurewebsites.net/"   # pinned here, once


class TestFeedbackSurvey(unittest.TestCase):
    """The Core's feedback survey (Brett, 2026-09-28): one URL and one link builder
    (core_submission.feedback_url / feedback_line), pre-filled with the PROT number only when it
    is one -- read with the one submission-number pattern, normalize_submission."""

    def test_the_url_prefills_only_a_real_prot_number(self):
        self.assertEqual(cs.FEEDBACK_URL, SURVEY)
        for given in ("PROT_0756", "prot-756", "756", "#756", "PROT_756"):
            with self.subTest(given):
                self.assertEqual(cs.feedback_url("report", given), SURVEY + "?prot=PROT_0756&src=report")
        for given in (None, "", "99922f5337f8", "submission", "PROT_0756&src=web", "PROT_0000",
                      "PROT_123456", "PROT_0756 x"):
            with self.subTest(given):
                self.assertEqual(cs.feedback_url("readme", given), SURVEY + "?src=readme")
        for src in cs.FEEDBACK_SOURCES:
            q = urllib.parse.parse_qs(urllib.parse.urlsplit(cs.feedback_url(src, "807")).query,
                                      strict_parsing=True)
            self.assertEqual(q, {"prot": ["PROT_0807"], "src": [src]})
        with self.assertRaises(ValueError):
            cs.feedback_url("slack", "807")

    def test_link_rewords_only_a_survey_line_that_is_there(self):
        # make_podcast.py link asks about the podcast too once there is one: the same place and
        # PROT number, the one wording; a report without the line, or with the same words but
        # no survey link, is left alone
        for fmt in ("html", "md"):
            with self.subTest(fmt):
                wrap = ("<main>x{}</main>" if fmt == "html" else "# R\n\nText.\n\n{}\n")
                page = wrap.format(cs.feedback_line("readme", "PROT_0807", fmt=fmt))
                new, n = cs.refresh_feedback(page, fmt, podcast=True)
                self.assertEqual((n, new), (1, wrap.format(cs.feedback_line(
                    "readme", "PROT_0807", podcast=True, fmt=fmt))))
                self.assertEqual(cs.refresh_feedback(new, fmt, podcast=True), (new, 0))
        for text, fmt in (("<main>no line</main>", "html"), ("# R\n\nText.\n", "md"),
                          ("**How did we do?** Very well, thanks.\n", "md"),
                          ('<p class="feedback">a <a href="https://example.org">x</a></p>', "html")):
            self.assertEqual(cs.refresh_feedback(text, fmt, podcast=True), (text, 0))

    def test_the_md_twins_line_is_stripped_with_its_rule(self):
        line = cs.feedback_line("report", "PROT_0756")
        body = "# R\n\nText.\n"
        self.assertEqual(cs.strip_feedback(body + "\n---\n\n" + line + "\n"), body + "\n")
        self.assertEqual(cs.strip_feedback(body + line + "\n"), body)
        other = body + "\n---\n\n**How did we do?** Very well.\n"    # not the survey's: kept
        self.assertEqual(cs.strip_feedback(other), other)

    def test_the_sources_the_form_reads(self):
        # review of f18d95e: the shareable file is made to be forwarded (src=share); "podcast"
        # was never used (no address is spoken)
        self.assertEqual(cs.FEEDBACK_SOURCES, ("report", "readme", "email", "share", "web"))
        with self.assertRaises(ValueError):
            cs.feedback_url("podcast")
        html_line = cs.feedback_line("report", "PROT_0756", fmt="html")
        new, n = cs.refresh_feedback(html_line, "html", podcast=False, src="share")
        self.assertEqual((n, new), (1, cs.feedback_line("share", "PROT_0756", fmt="html")))
        # what make_podcast's check hashes: the line gone in either form, so link rewording
        # it never makes a check stale
        body = "<main><p>Body 6,112.</p>{}</main>"
        for podcast in (False, True):
            self.assertEqual(cs.without_feedback(body.format(cs.feedback_line(
                "report", "PROT_0756", podcast, fmt="html"))), body.format(""))
            self.assertEqual(cs.without_feedback("# R\n\nText.\n\n---\n\n" + cs.feedback_line(
                "report", "PROT_0756", podcast) + "\n"), "# R\n\nText.\n\n")

    def test_a_share_state_that_cannot_be_read_is_skipped_not_fatal(self):
        with tempfile.TemporaryDirectory() as tmp:
            with mock.patch.object(cs.make_podcast, "share_state",
                                   side_effect=PermissionError(13, "Permission denied")):
                item = cs._plan_share(tmp)
        self.assertEqual(item, {"name": cs.SHARE_FILE, "status": "SKIPPED", "required": False,
                                "reason": "cannot check it is current (Permission denied)"})

    def test_the_line_in_each_format(self):
        self.assertEqual(cs.feedback_line("report", "PROT_0756"),
                         "**How did we do?** Tell us what you thought of this report: [a 5-minute "
                         f"survey]({SURVEY}?prot=PROT_0756&src=report).")
        self.assertEqual(cs.feedback_line("report", "PROT_0756", podcast=True, fmt="html"),
                         '<p class="feedback"><strong>How did we do?</strong> Tell us what you '
                         f'thought of this report and the podcast: <a href="{SURVEY}?prot=PROT_0756'
                         '&amp;src=report">a 5-minute survey</a></p>')
        self.assertEqual(cs.feedback_line("email", "PROT_0756", fmt="text"),
                         "How did we do? Tell us what you thought of the report: a 5-minute survey "
                         f"at {SURVEY}?prot=PROT_0756&src=email")

    def test_the_star_line_reads_its_address_from_plugin_json(self):
        # Brett (2026-09-29): every report asks for a GitHub star. The address is plugin.json's
        # `repository` (skill_version.plugin_meta, the one reader), never a copy of it
        import skill_version
        with open(os.path.join(SCRIPTS, "..", ".claude-plugin", "plugin.json")) as fh:
            repo = json.load(fh)["repository"]
        self.assertEqual(cs.repository_url(), repo)
        self.assertEqual(cs.star_line(), "**Found this report useful?** Please star [the DE-LIMP "
                         f"repository on GitHub]({repo}) — it helps other labs find these "
                         "tools.")
        self.assertEqual(cs.star_line("html"), '<p class="star"><strong>Found this report '
                         f'useful?</strong> Please star <a href="{repo}">the <span class="name">DE-LIMP'
                         "</span> repository on GitHub</a> &mdash; it helps other labs find these "
                         "tools.</p>")
        moved = {"repository": "https://github.com/someone/Moved/"}
        with mock.patch.object(skill_version, "plugin_meta", return_value=moved):
            self.assertIn("[the Moved repository on GitHub](https://github.com/someone/Moved) ",
                          cs.star_line())
        for meta in ({}, {"repository": ""}, {"repository": 7}, {"repository": cs.GITHUB},
                     {"repository": "https://example.org/x"}):      # no address: no line
            with mock.patch.object(skill_version, "plugin_meta", return_value=meta):
                self.assertIsNone(cs.repository_url(), meta)
                self.assertEqual((cs.star_line(), cs.star_line("html")), ("", ""), meta)

    def test_the_star_line_is_stripped_but_never_reworded(self):
        star, survey = cs.star_line(), cs.feedback_line("report", "PROT_0756")
        body = "# R\n\nText.\n"
        # the .md twin's footer (render_md): a rule, the star line, then a Core run's survey
        for tail in ("\n---\n\n" + star + "\n\n" + survey + "\n", "\n---\n\n" + star + "\n"):
            self.assertEqual(cs.strip_feedback(body + tail).rstrip("\n"), body.rstrip("\n"))
            self.assertEqual(cs.without_feedback(body + tail).rstrip("\n"), body.rstrip("\n"))
        other = body + "\n**Found this report useful?** Very.\n"    # links no GitHub: kept
        self.assertEqual(cs.strip_feedback(other), other)
        page = "<main><p>Body.</p>{}{}</main>"
        self.assertEqual(cs.without_feedback(page.format(cs.star_line("html"), cs.feedback_line(
            "report", "PROT_0756", fmt="html"))), page.format("", ""))
        # link rewording the survey line (and share's src=share) leaves the star line as it is
        for fmt in ("html", "md"):
            text = cs.star_line(fmt) + "\n\n" + cs.feedback_line("report", "PROT_0756", fmt=fmt)
            for src in (None, "share"):
                new, n = cs.refresh_feedback(text, fmt, podcast=True, src=src)
                self.assertEqual(n, 1)
                self.assertTrue(new.startswith(cs.star_line(fmt) + "\n\n"), (fmt, src))
                self.assertEqual(new.count("Found this report useful?"), 1)


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


# ---------------------------------------------------------------------- identify --
LRS_RUNS = [f"/data/JUL26/08132026__60SPD_DIA-LRS-{n}_S3-A1_1_{23600 + n}.d" for n in range(96, 126)]
HEX_0756 = "0022066cd85f"


def rec_0756(**kw):
    return record(internal_id="PROT_0756", sid=HEX_0756, submitted="2026-07-15T16:09:50.151363-07:00",
                  samples=[(f"LRS{n}", "c") for n in range(96, 126)], **kw)


class TestIdentifyByName(unittest.TestCase):
    """SKILL.md step 1: the submission is read from the names or the message, never guessed."""

    def test_prot_number_in_a_folder_name(self):
        self.assertEqual(cs.ids_in("/Data/lab/service/on_campus/Lab/PROT_0756/raw/x.d"),
                         [("internal_id", "PROT_0756")])
        self.assertEqual(cs.ids_in("prot-756_rerun"), [("internal_id", "PROT_0756")])

    def test_12_hex_id_in_a_path(self):
        self.assertEqual(cs.ids_in(f"/coreomics/projects/2026/07/{HEX_0756}/share/a.d"), [("id", HEX_0756)])

    def test_numbers_that_are_not_ids(self):
        # a 12-digit timestamp, an Exploris run counter, a hex-looking word with no digit, PROTEOMICS
        for s in ("FLsep26_wa_202609090050.raw", "Ex08312026_380_JE21.raw", "/x/deadbeefcafe/a.raw",
                  "/Volumes/proteomics/PROTEOMICS/a.d", "08132026__60SPD_DIA-LRS-96_S3-B1_1_23630.d"):
            self.assertEqual(cs.ids_in(s), [], s)
        for s in ("Total_prot_10ug_rep1.raw", "WT_prot1.raw", "Nuclear_Prot_50ug.raw",
                  "/tmp/2a31667e-609f-4d85-85ae-f9773f03ff76/raw/a.d"):
            self.assertEqual(cs.ids_in(s), [], s)
        self.assertEqual(cs.ids_in("the 2nd submission 3 weeks later", typed=True), [])
        self.assertEqual(cs.ids_in("submission 756 please"), [], "file-name rules for a path")
        self.assertEqual(cs.ids_in("please search submission #756", typed=True), [("internal_id", "PROT_0756")])

    def test_cli_named_none_and_ambiguous_without_a_lookup(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = child_env(tmp)
            rc, out, p = run(["identify", "--no-lookup", "/svc/PROT_0756/raw/" + os.path.basename(LRS_RUNS[0])], env)
            self.assertEqual((rc, out["status"], out["submission"]), (0, "named", "PROT_0756"), p.stderr)
            self.assertIn("Is that right?", out["ask"])
            self.assertIn(".submissions_db", out["never_search"])
            rc, out, _ = run(["identify", "--no-lookup"] + LRS_RUNS, env)
            self.assertEqual((rc, out["status"]), (2, "none"))
            # no evidence of Core data: ask whether the Core ran them, not "which submission?"
            # (an outside user answering with any number would become a Core run)
            self.assertEqual(out["ask"], "Were these samples run by the UC Davis Proteomics "
                                         "Core? If so, what is the PROT number?")
            rc, out, _ = run(["identify", "--no-lookup", "/svc/PROT_0756/a.d", "--text", "or is it PROT_0757?"], env)
            self.assertEqual((rc, out["status"]), (2, "ambiguous"))
            self.assertEqual({c["submission"] for c in out["candidates"]}, {"PROT_0756", "PROT_0757"})

    def test_no_token_asks_and_says_how_to_get_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = child_env(tmp, HOME=tmp)
            env.pop("COREOMICS_TOKEN")
            rc, out, _ = run(["identify"] + LRS_RUNS, env)
            self.assertEqual((rc, out["status"]), (3, "needs_token"))
            self.assertIn("~/.coreomics_token", out["token_help"])
            self.assertTrue(out["ask"].startswith("Were these samples run by the UC Davis "
                                                  "Proteomics Core? If so, what is the PROT "
                                                  "number?"), out["ask"])


class TestIdentifyBySampleIds(unittest.TestCase):
    def serve(self, fake, pages):
        def listing(q, b):
            if q.get("internal_id"):
                hits = [r for p in pages for r in p if r["internal_id"] == q["internal_id"]]
                return 200, {"count": len(hits), "next": None, "results": hits}
            page = int(q.get("page", 1))
            nxt = f"{fake.base}/submissions/?page={page + 1}" if page < len(pages) else None
            return 200, {"count": sum(map(len, pages)), "next": nxt, "results": pages[page - 1]}
        fake.routes[("GET", "/server/api/submissions/")] = listing
        for p in pages:
            for r in p:
                fake.routes[("GET", f"/server/api/submissions/{r['id']}/")] = (lambda r: lambda q, b: (200, r))(r)

    def test_prot0756_runs_name_one_submission(self):
        other = record(internal_id="PROT_0760", sid="bbbbbbbbbbb1", submitted="2026-07-20T10:00:00-07:00",
                       samples=[("KG1", "x"), ("A3", "y")])
        old = record(internal_id="PROT_0600", sid="bbbbbbbbbbb2", submitted="2026-02-01T10:00:00-08:00",
                     samples=[("LRS66", "x")])
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, [[other, rec_0756(), old]])
            rc, out, p = run(["identify"] + LRS_RUNS, child_env(tmp, fake.base))
            self.assertEqual((rc, out["status"], out["submission"]), (0, "matched", "PROT_0756"), p.stderr)
            c = out["candidates"][0]
            self.assertEqual((c["files_matched"], c["sample_ids_matched"], c["samples_on_sheet"]), (30, 30, 30))
            self.assertIn("30 of 30", out["ask"])

    def test_a_label_two_submissions_use_is_ambiguous(self):
        """locate's ambiguous_label, in reverse: both were submitted before the runs."""
        twin = record(internal_id="PROT_0755", sid="bbbbbbbbbbb3", submitted="2026-07-10T10:00:00-07:00",
                      samples=[(f"LRS{n}", "c") for n in range(96, 100)])
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, [[rec_0756(), twin]])
            rc, out, _ = run(["identify"] + LRS_RUNS, child_env(tmp, fake.base))
            self.assertEqual((rc, out["status"]), (2, "ambiguous"))
            self.assertEqual({f["file"] for f in out["ambiguous_files"]},
                             {os.path.basename(x) for x in LRS_RUNS[:4]})
            self.assertIn("PROT_0755", out["ask"])

    def test_a_submission_made_after_the_runs_does_not_claim_them(self):
        later = record(internal_id="PROT_0799", sid="bbbbbbbbbbb4", submitted="2026-09-01T10:00:00-07:00",
                       samples=[(f"LRS{n}", "c") for n in range(96, 126)])
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, [[later, rec_0756()]])
            rc, out, _ = run(["identify"] + LRS_RUNS, child_env(tmp, fake.base))
            self.assertEqual((rc, out["submission"]), (0, "PROT_0756"))

    def test_one_generic_id_in_thirty_files_is_not_enough(self):
        qc = record(internal_id="PROT_0801", sid="bbbbbbbbbbb6", submitted="2026-07-01T10:00:00-07:00",
                    samples=[("HeLa", "x"), ("Blank", "y"), ("KG77", "z")])
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, [[qc]])
            names = ["/d/08132026__60SPD_DIA-HeLa_S3-A1_1_1.d"] + LRS_RUNS[:29]
            rc, out, _ = run(["identify"] + names, child_env(tmp, fake.base))
            self.assertEqual((rc, out["status"], out["submission"]), (2, "weak", "PROT_0801"))

    def test_weak_ids_are_no_evidence(self):
        wells = record(internal_id="PROT_0761", sid="bbbbbbbbbbb5", submitted="2026-07-01T10:00:00-07:00",
                       samples=[("A1", "x"), ("001", "y")])
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, [[wells]])
            rc, out, _ = run(["identify", "/d/08132026__60SPD_DIA-A1_S3-A1_1_1.d", "/d/Ex08152026_001.raw"],
                             child_env(tmp, fake.base))
            self.assertEqual((rc, out["status"]), (2, "none"))

    def test_a_prot_number_and_its_hex_id_are_merged_by_lookup(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, [[rec_0756()]])
            rc, out, p = run(["identify", f"/coreomics/projects/2026/07/{HEX_0756}/share/x.d"] + LRS_RUNS
                             + ["--text", "this is PROT_0756"], child_env(tmp, fake.base))
            self.assertEqual((rc, out["status"], out["submission"]), (0, "named", "PROT_0756"), p.stderr)
            self.assertTrue(out["verified"])
            self.assertEqual(out["files_matched"], 30)

    def test_a_named_submission_whose_ids_are_in_no_file_name_is_asked_about(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, [[rec_0756()]])
            rc, out, _ = run(["identify", "/svc/PROT_0756/Ex08152026_1_other.raw"], child_env(tmp, fake.base))
            self.assertEqual((rc, out["status"]), (2, "named_unconfirmed"))
            self.assertIn("none of its 30 sample IDs", out["ask"])

    def test_a_named_submission_that_does_not_exist_is_asked_about(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, [[rec_0756()]])
            rc, out, _ = run(["identify", "/svc/PROT_9999/a.d"], child_env(tmp, fake.base))
            self.assertEqual((rc, out["status"]), (2, "not_found"))

    def test_a_failed_lookup_says_so_and_how_to_fix_the_token(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            fake.routes[("GET", "/server/api/submissions/")] = lambda q, b: (401, {"detail": "Invalid token."})
            rc, out, _ = run(["identify"] + LRS_RUNS, child_env(tmp, fake.base))
            self.assertEqual((rc, out["status"]), (3, "lookup_failed"))
            self.assertIn("Invalid token", out["ask"])
            self.assertIn("~/.coreomics_token", out["token_help"])

    def test_a_folder_is_listed_for_its_raw_file_names(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.serve(fake, [[rec_0756()]])
            folder = os.path.join(tmp, "JUL26")
            for r in LRS_RUNS:
                os.makedirs(os.path.join(folder, os.path.basename(r)))
            rc, out, p = run(["identify", folder], child_env(tmp, fake.base))
            self.assertEqual((rc, out["status"], out["submission"]), (0, "matched", "PROT_0756"), p.stderr)

    def test_a_cut_off_submission_list_never_confirms(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            later = record(internal_id="PROT_0799", sid="bbbbbbbbbbb7", submitted="2026-06-01T10:00:00-07:00",
                           samples=[("ZZ1", "x")])
            self.serve(fake, [[rec_0756()], [later]])
            with mock.patch.object(cs, "NEIGHBOR_MAX_PAGES", 1):
                with mock.patch.dict(os.environ, child_env(tmp, fake.base), clear=True):
                    with contextlib.redirect_stdout(io.StringIO()) as buf:
                        rc = cs.main(["identify"] + LRS_RUNS)
            out = json.loads(buf.getvalue())
            self.assertEqual((rc, out["status"]), (2, "weak"))
            self.assertTrue(out["window"]["truncated"])


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


# ------------------------------------------------------------ the CoreOmics API key --
KEY = "k3yNotReal0a1b2c3d4e5f60718293a4b5c6d7e8"     # a planted key: it must never be printed
WIN_ENV = {"USERPROFILE": "C:\\Users\\msalemi", "HOMEDRIVE": "H:", "HOMEPATH": "\\"}


def key_home(tmp, content=None, name=".coreomics_token", raw=None):
    home = os.path.join(tmp, "home")
    os.makedirs(home, exist_ok=True)
    if content is not None or raw is not None:
        with open(os.path.join(home, name), "wb") as fh:
            fh.write(raw if raw is not None else content.encode())
    return home


def key_env(tmp, base_url="http://127.0.0.1:9/server/api", **extra):
    """child_env whose key comes from the files under HOME, unless COREOMICS_TOKEN is passed."""
    e = child_env(tmp, base_url, HOME=key_home(tmp))
    e.pop("COREOMICS_TOKEN")
    e.update(extra)
    return e


def windows_places(env):
    """key_files() and stray_key_files() as Windows Python computes them (ntpath's expanduser
    reads USERPROFILE and never HOME) -- only the environment named here."""
    with mock.patch.dict(os.environ, env, clear=True), mock.patch.object(cs.os, "name", "nt"), \
            mock.patch.object(cs.os, "path", ntpath):
        return cs.key_files(), cs.stray_key_files(), cs.key_target()


class TestKeyPlaces(unittest.TestCase):
    """Windows: Python's ~ is USERPROFILE, Git Bash's ~ is HOME, and on a domain account
    (Michelle's, "AD3+msalemi") they can differ -- a key saved with Git Bash's ~ must be found."""

    def test_windows_reads_the_profile_folder_then_git_bash_home(self):
        read, _, target = windows_places(dict(WIN_ENV, HOME="H:\\"))
        self.assertEqual(read, ["C:\\Users\\msalemi\\.coreomics_token", "H:\\.coreomics_token"])
        self.assertEqual(target, "C:\\Users\\msalemi\\.coreomics_token", "save where Python reads first")

    def test_windows_home_in_the_profile_folder_is_one_place(self):
        for home in ("C:\\Users\\msalemi", "c:/users/MSALEMI/", "/c/Users/msalemi"):
            read, _, _ = windows_places(dict(WIN_ENV, HOME=home))
            self.assertEqual(read, ["C:\\Users\\msalemi\\.coreomics_token"], home)

    def test_windows_an_msys_spelled_home_is_its_drive(self):
        read, _, _ = windows_places(dict(WIN_ENV, HOME="/h/"))
        self.assertEqual(read[1], "H:\\.coreomics_token")
        self.assertEqual(cs._drive_path("/home/msalemi"), "/home/msalemi", "not a drive letter")

    def test_windows_without_home_reads_the_profile_only(self):
        read, stray, _ = windows_places(WIN_ENV)
        self.assertEqual(read, ["C:\\Users\\msalemi\\.coreomics_token"])
        # the home drive and Notepad's .txt are looked at by check, never read as the key
        self.assertIn("H:\\.coreomics_token", stray)
        self.assertIn("C:\\Users\\msalemi\\.coreomics_token.txt", stray)
        self.assertIn("C:\\Users\\msalemi\\coreomics_token", stray)
        self.assertNotIn("C:\\Users\\msalemi\\.coreomics_token", stray)

    def test_posix_reads_home_and_ignores_userprofile(self):
        with mock.patch.dict(os.environ, {"HOME": "/Users/msalemi", "USERPROFILE": "C:\\Users\\x"},
                             clear=True), mock.patch.object(cs.os, "name", "posix"):
            self.assertEqual(cs.key_files(), ["/Users/msalemi/.coreomics_token"])

    def test_the_save_line_writes_where_python_reads(self):
        self.assertIn('cygpath "$USERPROFILE"', cs.SAVE_KEY_LINE["windows"])
        self.assertIn('"$HOME/.coreomics_token"', cs.SAVE_KEY_LINE["posix"])
        for plat, line in cs.SAVE_KEY_LINE.items():
            self.assertNotIn("\n", line, plat)                        # ONE line (TestSaveLine)
            self.assertTrue(line.startswith("read -rs TOK && "), plat) # silent; stops on EOF
            self.assertIn('"$TOK"', line)                             # never the key itself
            self.assertTrue(line.endswith("; unset TOK"), plat)
            self.assertIn('chmod 600 "$F"', line)
        with mock.patch.object(cs.os, "name", "nt"):
            self.assertIn('cygpath "$USERPROFILE"', cs.move_key_steps("H:\\.coreomics_token"))

    def test_windows_documents_and_onedrive_are_looked_at(self):
        """Notepad's Save As usually starts in Documents -- OneDrive's, on a OneDrive PC."""
        _, stray, _ = windows_places(dict(WIN_ENV, OneDrive="C:\\Users\\msalemi\\OneDrive - UC Davis"))
        for p in ("C:\\Users\\msalemi\\Documents\\.coreomics_token.txt",
                  "C:\\Users\\msalemi\\OneDrive\\Documents\\coreomics_token.txt",
                  "C:\\Users\\msalemi\\OneDrive - UC Davis\\Documents\\.coreomics_token.txt"):
            self.assertIn(p, stray)

    def test_macos_documents_is_never_looked_at(self):
        """~/Documents on a Mac is privacy-guarded (a prompt) and often iCloud: not for a guess."""
        with mock.patch.dict(os.environ, {"HOME": "/Users/msalemi"}, clear=True), \
                mock.patch.object(cs.os, "name", "posix"):
            self.assertFalse([p for p in cs.stray_key_files() if "Documents" in p])

    def test_git_bash_home_key_is_read_when_the_profile_has_none(self):
        with tempfile.TemporaryDirectory() as tmp:
            profile, home = os.path.join(tmp, "profile"), key_home(tmp, KEY)
            os.makedirs(profile)
            with mock.patch.dict(os.environ, {"HOME": home}), \
                    mock.patch.object(cs.os, "name", "nt"), \
                    mock.patch.object(cs.os.path, "expanduser",
                                      lambda p: p.replace("~", profile, 1)):
                os.environ.pop("COREOMICS_TOKEN", None)
                self.assertEqual(cs.find_key(), (KEY, os.path.join(home, ".coreomics_token")))
                with open(os.path.join(profile, ".coreomics_token"), "w") as fh:
                    fh.write("\n")                       # an empty profile file hides nothing
                self.assertEqual(cs.api_token(), KEY)
                with open(os.path.join(profile, ".coreomics_token"), "w") as fh:
                    fh.write("Token " + KEY)              # nor does a botched one
                self.assertEqual(cs.api_token(), KEY)


def _shells():
    """(name, argv) for the shells a staff member pastes into: /bin/bash (3.2 on macOS, no
    bracketed paste), zsh, and a bash 4/5 like Git Bash's when one is installed."""
    out = []
    for name in ("bash", "zsh"):
        for path in dict.fromkeys(filter(None, ["/bin/" + name, shutil.which(name)])):
            if os.path.exists(path):
                v = subprocess.run([path, "-c", 'echo "$BASH_VERSION$ZSH_VERSION"'],
                                   capture_output=True, text=True).stdout.strip()
                out.append((f"{name} {v}", path))
    return out


class TestSaveLine(unittest.TestCase):
    """The save line in real shells, fed as a terminal without bracketed paste feeds it: the
    line, then the key typed at the silent prompt. The key must land in the file and never
    run as a command (two lines did: `read` took the second line, the key ran)."""

    def paste(self, shell, line, key, env):
        p = subprocess.run([shell, "-s"], input=f"{line}\n{key}\n", capture_output=True,
                           text=True, env=env, timeout=30)
        return p

    def check_saved(self, path):
        with open(path) as fh:
            self.assertEqual(fh.read(), KEY + "\n")
        self.assertEqual(stat.S_IMODE(os.stat(path).st_mode), 0o600)

    def test_posix_line(self):
        shells = _shells()
        self.assertTrue(shells)
        for name, sh in shells:
            if "zsh" in name:
                continue                  # zsh -s buffers its script; tested with -c below
            with tempfile.TemporaryDirectory() as tmp:
                env = {"PATH": os.environ["PATH"], "HOME": tmp}
                p = self.paste(sh, cs.SAVE_KEY_LINE["posix"], KEY, env)
                self.assertEqual(p.returncode, 0, name + p.stderr)
                self.assertNotIn(KEY, p.stdout + p.stderr, name)
                self.check_saved(os.path.join(tmp, ".coreomics_token"))
                # the hazard this guards against, shown on the same shell
                old = ("read -rs TOK\nF=\"$HOME/.coreomics_token\"; (umask 077; printf '%s\\n' "
                       "\"$TOK\" > \"$F\"); unset TOK; chmod 600 \"$F\"")
                q = self.paste(sh, old, KEY, env)
                self.assertIn(KEY, q.stdout + q.stderr, name)
        for name, sh in shells:
            if "zsh" in name:
                with tempfile.TemporaryDirectory() as tmp:
                    p = subprocess.run([sh, "-c", cs.SAVE_KEY_LINE["posix"]], input=KEY + "\n",
                                       capture_output=True, text=True, timeout=30,
                                       env={"PATH": os.environ["PATH"], "HOME": tmp})
                    self.assertEqual(p.returncode, 0, name + p.stderr)
                    self.assertNotIn(KEY, p.stdout + p.stderr, name)
                    self.check_saved(os.path.join(tmp, ".coreomics_token"))

    def test_git_bash_line_with_a_space_in_the_profile(self):
        for name, sh in _shells():
            if not name.startswith("bash"):
                continue
            with tempfile.TemporaryDirectory() as tmp:
                profile = os.path.join(tmp, "Mary O Brien")
                os.makedirs(profile)
                stub = "cygpath() { printf '%s\\n' \"$1\"; }; "      # Git Bash's, for a POSIX path
                p = self.paste(sh, stub + cs.SAVE_KEY_LINE["windows"], KEY,
                               {"PATH": os.environ["PATH"], "USERPROFILE": profile})
                self.assertEqual(p.returncode, 0, name + p.stderr)
                self.assertNotIn(KEY, p.stdout + p.stderr, name)
                self.check_saved(os.path.join(profile, ".coreomics_token"))

    def test_nothing_is_written_or_chmodded_when_read_gets_no_key(self):
        """Ctrl-D at the prompt leaves an $F from earlier alone; Enter before the paste leaves
        a saved key alone."""
        for name, sh in _shells():
            with tempfile.TemporaryDirectory() as tmp:
                other = os.path.join(tmp, "other.sh")
                with open(other, "w") as fh:
                    fh.write("x")
                os.chmod(other, 0o755)
                env = {"PATH": os.environ["PATH"], "HOME": tmp}
                p = subprocess.run([sh, "-c", f"F={shlex.quote(other)}; " + cs.SAVE_KEY_LINE["posix"]],
                                   input="", capture_output=True, text=True, timeout=30, env=env)
                self.assertEqual(stat.S_IMODE(os.stat(other).st_mode), 0o755, name + p.stderr)
                self.assertFalse(os.path.exists(os.path.join(tmp, ".coreomics_token")), name)
                saved = write(tmp, ".coreomics_token", KEY + "\n")
                subprocess.run([sh, "-c", cs.SAVE_KEY_LINE["posix"]], input="\n",
                               capture_output=True, text=True, timeout=30, env=env)
                self.assertEqual(read(saved), KEY + "\n", name)


class TestKeyContents(unittest.TestCase):
    def test_utf16_and_bom_files_hold_the_key(self):
        """Windows PowerShell's `>` writes UTF-16; Notepad can add a UTF-8 BOM."""
        with tempfile.TemporaryDirectory() as tmp:
            for raw in ((KEY + "\r\n").encode("utf-16"), codecs.BOM_UTF8 + (KEY + "\r\n").encode()):
                p = os.path.join(tmp, "k")
                with open(p, "wb") as fh:
                    fh.write(raw)
                self.assertEqual(cs.read_key_file(p), KEY)

    def test_a_malformed_key_never_reaches_a_header(self):
        for bad in ("Token " + KEY, KEY[:20] + " " + KEY[20:], KEY[:20] + "\n" + KEY[20:],
                    '"' + KEY + '"', KEY + "\u2019"):
            with tempfile.TemporaryDirectory() as tmp, \
                    mock.patch.dict(os.environ, {"HOME": key_home(tmp, bad)}):
                os.environ.pop("COREOMICS_TOKEN", None)
                with self.assertRaises(cs.KeyProblem) as e:
                    cs.api_token()
                self.assertEqual(e.exception.diagnosis, "key_malformed", repr(bad))
                for text in (e.exception.detail, e.exception.fix, str(e.exception)):
                    self.assertNotIn(KEY[:20], text)

    def test_fetch_with_a_two_line_key_exits_3_without_printing_it(self):
        """http.client's error for a line break in a header is "Invalid header value b'Token
        <key>'": before the shape check, fetch died with a traceback holding the key."""
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            env = key_env(tmp, fake.base)
            with open(os.path.join(env["HOME"], ".coreomics_token"), "w") as fh:
                fh.write(KEY[:20] + "\n" + KEY[20:] + "\n")
            rc, out, p = run(["fetch", "807", "--out", os.path.join(tmp, "o")], env)
            self.assertEqual((rc, out["diagnosis"]), (3, "key_malformed"), p.stderr)
            self.assertNotIn(KEY[:20], p.stdout + p.stderr)
            self.assertNotIn(KEY[20:], p.stdout + p.stderr)
            self.assertEqual(fake.requests, [], "nothing was sent")


class TestCheck(unittest.TestCase):
    """`check`: one diagnosis with the fix -- and never the key, whatever the server says."""

    def check(self, env, *extra):
        rc, out, p = run(["check", "--json", *extra], env)
        for part in (KEY, KEY[:20], KEY[20:]):
            self.assertNotIn(part, p.stdout + p.stderr)
        return rc, out, p

    def lab(self, fake, reply):
        fake.routes[("GET", "/server/api/submissions/")] = lambda q, b: reply

    def assert_save_steps(self, fix):
        for s in ("Profile", '"API Key"', "Create", "read -rs TOK", 'chmod 600 "$F"',
                  "Never paste the key into the chat", "coreomics_key_to_hive.sh", "mode 600"):
            self.assertIn(s, fix)

    def test_ok_asks_one_lab_scoped_page_with_the_key(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.lab(fake, (200, {"count": 812, "next": "x", "results": [record()]}))
            env = key_env(tmp, fake.base)
            key_home(tmp, KEY + "\n")
            rc, out, p = self.check(env)
            self.assertEqual((rc, out["status"], out["fix"]), (0, "ok", ""), p.stderr)
            self.assertEqual(out["key_source"], os.path.join(env["HOME"], ".coreomics_token"))
            (req,) = fake.requests
            self.assertEqual(req["query"], {"lab": "PROTEOMICS", "page_size": "1"})
            self.assertEqual(req["auth"], "Token " + KEY)

    def test_no_key_names_where_it_looked_and_how_to_make_one(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            env = key_env(tmp, fake.base)
            rc, out, _ = self.check(env)
            self.assertEqual((rc, out["status"]), (3, "no_key"))
            self.assertEqual(out["looked_in"], [os.path.join(env["HOME"], ".coreomics_token")])
            self.assertIn(out["looked_in"][0], out["say"])
            self.assert_save_steps(out["fix"])
            self.assertEqual(fake.requests, [], "no key, no call")

    def test_the_plain_output_is_for_a_person(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = subprocess.run([sys.executable, CS, "check"], capture_output=True, text=True,
                               env=key_env(tmp), timeout=60)
            self.assertEqual(p.returncode, 3)
            self.assertTrue(p.stdout.startswith("CoreOmics key check: no_key\n"), p.stdout)
            self.assertIn("What to do:", p.stdout)
            self.assertIn("read -rs TOK", p.stdout)

    def test_empty_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = key_env(tmp)
            key_home(tmp, "\n")
            rc, out, _ = self.check(env)
            self.assertEqual((rc, out["status"]), (3, "key_empty"))
            self.assert_save_steps(out["fix"])

    @unittest.skipIf(IS_ROOT, "root reads a mode-000 file")
    def test_an_unreadable_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = key_env(tmp)
            key_home(tmp, KEY)
            path = os.path.join(env["HOME"], ".coreomics_token")
            os.chmod(path, 0)
            try:
                rc, out, _ = self.check(env)
            finally:
                os.chmod(path, 0o600)
            self.assertEqual((rc, out["status"]), (3, "key_unreadable"))
            self.assertIn(path, out["say"])

    def test_a_key_saved_by_notepad_is_in_the_wrong_place(self):
        with tempfile.TemporaryDirectory() as tmp:
            env = key_env(tmp)
            key_home(tmp, KEY, name=".coreomics_token.txt")
            rc, out, _ = self.check(env)
            self.assertEqual((rc, out["status"]), (3, "key_in_wrong_place"))
            src, dst = (os.path.join(env["HOME"], n) for n in (".coreomics_token.txt", ".coreomics_token"))
            self.assertIn(f"mv {shlex.quote(src)} {shlex.quote(dst)}", out["fix"])
            self.assertIn(src, out["say"])

    def test_rejected_key(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.lab(fake, (403, {"detail": "Invalid token."}))
            env = key_env(tmp, fake.base)
            key_home(tmp, KEY)
            rc, out, _ = self.check(env)
            self.assertEqual((rc, out["status"], out["http_status"]), (3, "key_rejected", 403))
            self.assertIn("Invalid token.", out["say"])
            self.assert_save_steps(out["fix"])

    def test_rejected_key_from_the_variable_says_the_variable_wins(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.lab(fake, (401, {"detail": "Invalid token."}))
            rc, out, _ = self.check(key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
            self.assertEqual((rc, out["status"]), (3, "key_rejected"))
            self.assertIn("COREOMICS_TOKEN is set here", out["fix"])

    def test_a_server_that_echoes_the_key_is_redacted(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.lab(fake, (403, {"detail": f"Invalid token {KEY}."}))
            rc, out, _ = self.check(key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
            self.assertEqual(out["detail"], "Invalid token <key>.")

    def test_no_lab_access_two_ways(self):
        for reply in ((200, {"count": 0, "next": None, "results": []}),
                      (403, {"detail": "You do not have permission to perform this action."})):
            with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
                self.lab(fake, reply)
                rc, out, _ = self.check(key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
                self.assertEqual((rc, out["status"]), (3, "no_lab_access"), reply)
                self.assertIn("brettsp", out["fix"])
                self.assertIn("Proteomics lab", out["fix"])

    def test_a_web_page_401_or_403_is_not_about_the_key_or_the_lab(self):
        """Apache's 403 page says "permission"; a proxy says "Access Denied". Only CoreOmics'
        own JSON detail classifies a 401/403."""
        pages = {403: b"<!DOCTYPE HTML PUBLIC><html><head><title>403 Forbidden</title></head><body><h1>Forbidden</h1><p>You don't have permission to access this resource.</p></body></html>",
                 401: b"<html><head><title>Access Denied</title></head><body>Access Denied</body></html>"}
        for status, page in pages.items():
            with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
                self.lab(fake, (status, page))
                rc, out, _ = self.check(key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
                self.assertEqual((rc, out["status"], out["http_status"]), (3, "unexpected", status))
                self.assertIn(f"HTTP {status}", out["say"])
                self.assertNotIn("<", out["say"] + out["detail"])
                self.assertNotIn("brettsp", out["fix"])           # not an account problem

    def test_another_coreomics_refusal_is_quoted_not_guessed(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.lab(fake, (401, {"detail": "User inactive or deleted."}))
            rc, out, _ = self.check(key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
            self.assertEqual((rc, out["status"]), (3, "unexpected"))
            self.assertIn("User inactive or deleted.", out["say"])
            self.assertIn("brettsp", out["fix"])

    def test_a_redirect_to_another_server_never_carries_the_key(self):
        """urllib copies Authorization to wherever a redirect points: before the fix a 302 to a
        second server handed it the key and check said ok."""
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake, FakeCoreOmics() as other:
            other.routes[("GET", "/server/api/submissions/")] = lambda q, b: (200, {"count": 9})
            away = {"Location": other.base + "/submissions/?lab=PROTEOMICS&page_size=1"}
            self.lab(fake, (302, {}, away))
            env = key_env(tmp, fake.base, COREOMICS_TOKEN=KEY)
            rc, out, _ = self.check(env)
            self.assertEqual((rc, out["status"]), (3, "unexpected"))
            self.assertIn("redirect", out["say"])
            rc, _, p = run(["fetch", "807", "--out", os.path.join(tmp, "o")], env)
            self.assertEqual(rc, 3, p.stderr)
            self.assertEqual(other.requests, [], "the other server was never asked")

    def test_a_redirect_with_a_malformed_port_is_a_refused_redirect(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.lab(fake, (302, {}, {"Location": "http://127.0.0.1:abc/x"}))
            rc, out, _ = self.check(key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
            self.assertEqual((rc, out["status"]), (3, "unexpected"))
            self.assertIn("redirect", out["say"])
            self.assertNotIn("COREOMICS_BASE_URL", out["say"])

    def test_a_redirect_on_coreomics_itself_is_followed(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.lab(fake, (301, {}, {"Location": "/server/api/moved/?lab=PROTEOMICS&page_size=1"}))
            fake.routes[("GET", "/server/api/moved/")] = lambda q, b: (200, {"count": 4})
            rc, out, p = self.check(key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
            self.assertEqual((rc, out["status"]), (0, "ok"), p.stderr)
            self.assertEqual(fake.requests[-1]["auth"], "Token " + KEY)

    def test_a_base_url_that_is_not_a_web_address_is_a_diagnosis(self):
        with tempfile.TemporaryDirectory() as tmp:
            rc, out, p = self.check(key_env(tmp, "ucdavis.coreomics.com/server/api", COREOMICS_TOKEN=KEY))
            self.assertEqual((rc, out["status"]), (3, "unexpected"), p.stderr)
            self.assertIn("COREOMICS_BASE_URL", out["say"] + out["fix"])
            self.assertNotIn("Traceback", p.stderr)

    def test_no_key_arriving_is_not_blamed_on_the_key(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.lab(fake, (403, {"detail": "Authentication credentials were not provided."}))
            rc, out, _ = self.check(key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
            self.assertEqual((rc, out["status"]), (3, "unexpected"))

    def test_down_or_refused_is_unreachable(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            self.lab(fake, (503, {"detail": "Service Unavailable"}))
            rc, out, _ = self.check(key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
            self.assertEqual((rc, out["status"]), (3, "unreachable"))
        with tempfile.TemporaryDirectory() as tmp:           # port 9: connection refused
            rc, out, _ = self.check(key_env(tmp, COREOMICS_TOKEN=KEY))
            self.assertEqual((rc, out["status"]), (3, "unreachable"))
            self.assertIn("VPN", out["fix"])
            self.assertIn("not confirmed", out["fix"])


class _Opener:
    """An injected transport: raises `exc`, or answers `body` with 200."""

    def __init__(self, exc=None, body=b""):
        self.exc, self.body = exc, body

    def open(self, req, timeout=None):
        if self.exc is not None:
            raise self.exc
        return io.BytesIO(self.body)


class TestCheckTransport(unittest.TestCase):
    """DNS, timeout, TLS and a captive portal, in-process with an injected transport."""

    def diagnose(self, opener):
        env = {"HOME": tempfile.gettempdir(), "COREOMICS_TOKEN": KEY,
               "COREOMICS_BASE_URL": "https://coreomics.invalid/server/api"}
        with mock.patch.dict(os.environ, env), mock.patch.object(cs, "_OPENER", opener):
            out = io.StringIO()
            with contextlib.redirect_stdout(out), contextlib.redirect_stderr(io.StringIO()) as err:
                rc = cs.main(["check", "--json"])
        self.assertNotIn(KEY, out.getvalue() + err.getvalue())
        return rc, json.loads(out.getvalue())

    def test_dns_timeout_and_tls_are_unreachable_and_say_which(self):
        import socket
        import ssl
        cases = {"DNS": urllib.error.URLError(socket.gaierror(8, "nodename nor servname provided")),
                 "timed out": urllib.error.URLError(socket.timeout("timed out")),
                 "secure connection": urllib.error.URLError(ssl.SSLCertVerificationError(
                     1, "[SSL: CERTIFICATE_VERIFY_FAILED] certificate verify failed"))}
        for words, exc in cases.items():
            rc, out = self.diagnose(_Opener(exc=exc))
            self.assertEqual((rc, out["status"]), (3, "unreachable"), words)
            self.assertIn(words, out["say"])
        rc, out = self.diagnose(_Opener(exc=socket.timeout("timed out")))   # mid-read
        self.assertEqual(out["status"], "unreachable")

    def test_a_sign_in_page_is_not_coreomics(self):
        rc, out = self.diagnose(_Opener(body=b"<html>Sign in to the guest Wi-Fi</html>"))
        self.assertEqual((rc, out["status"]), (3, "unexpected"))
        self.assertIn("web page", out["say"])

    def test_not_http_at_all_and_anything_else_never_crash(self):
        rc, out = self.diagnose(_Opener(exc=http.client.BadStatusLine("HELLO")))
        self.assertEqual((rc, out["status"]), (3, "unexpected"))
        self.assertIn("not a web server", out["say"])
        rc, out = self.diagnose(_Opener(exc=RuntimeError("boom " + KEY)))      # a bug, reported
        self.assertEqual((rc, out["status"]), (3, "unexpected"))
        self.assertIn("RuntimeError: boom <key>", out["say"])

    def test_ok_in_process(self):
        rc, out = self.diagnose(_Opener(body=b'{"count": 3, "next": null, "results": []}'))
        self.assertEqual((rc, out["status"], out["submissions_visible"]), (0, "ok", 3))


class TestKeyNeverPrinted(unittest.TestCase):
    """api_call blanks the key out of everything a server or a library says, so fetch,
    identify and bioshare cannot print it either."""

    def assert_clean(self, p, *paths):
        self.assertNotIn(KEY, p.stdout + p.stderr)
        for root in paths:
            for d, _, fs in os.walk(root):
                for f in fs:
                    self.assertNotIn(KEY, read(os.path.join(d, f)), f)

    def test_fetch_error_replies_that_echo_the_key(self):
        for reply in ((403, {"detail": f"Invalid token {KEY}."}),
                      (400, {"internal_id": [f"bad value; your token is {KEY}"]}),
                      (200, f"<html>signed in as {KEY}</html>".encode())):
            with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
                fake.routes[("GET", "/server/api/submissions/")] = lambda q, b, r=reply: r
                rc, out, p = run(["fetch", "807", "--out", os.path.join(tmp, "o")],
                                 key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
                self.assertEqual(rc, 3, p.stderr)
                self.assert_clean(p)
                self.assertIn("<key>", p.stdout, reply[0])

    def test_a_record_that_holds_the_key_is_written_without_it(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            rec = record()
            rec["submission_data"]["description"] = f"token {KEY} pasted by mistake"
            fake.routes[("GET", "/server/api/submissions/")] = \
                lambda q, b: (200, {"count": 1, "next": None, "results": [rec]})
            fake.routes[("GET", SHARES)] = lambda q, b: (200, [])
            out_dir = os.path.join(tmp, "o")
            rc, _, p = run(["fetch", "807", "--out", out_dir], key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
            self.assertEqual(rc, 0, p.stderr)
            self.assert_clean(p, out_dir)

    def test_identify_lookup_failure_that_echoes_the_key(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            fake.routes[("GET", "/server/api/submissions/")] = lambda q, b: (403, {"detail": f"Invalid token {KEY}."})
            rc, out, p = run(["identify", "/svc/PROT_0756/a.d"], key_env(tmp, fake.base, COREOMICS_TOKEN=KEY))
            self.assertEqual((rc, out["status"]), (3, "lookup_failed"), p.stderr)
            self.assert_clean(p)


class TestKeyHelpElsewhere(unittest.TestCase):
    def test_fetch_on_a_rejected_key_points_at_check(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            fake.routes[("GET", "/server/api/submissions/")] = lambda q, b: (403, {"detail": "Invalid token."})
            rc, out, p = run(["fetch", "807", "--out", os.path.join(tmp, "o")], child_env(tmp, fake.base))
            self.assertEqual(rc, 3, p.stderr)
            self.assertIn("core_submission.py check", out["hint"])

    def test_fetch_without_a_key_gives_the_steps(self):
        with tempfile.TemporaryDirectory() as tmp:
            rc, out, _ = run(["fetch", "807", "--out", os.path.join(tmp, "o")], key_env(tmp))
            self.assertEqual((rc, out["diagnosis"]), (3, "no_key"))
            self.assertIn("read -rs TOK", out["fix"])

    def test_token_help_is_the_steps_not_ask_someone(self):
        with tempfile.TemporaryDirectory() as tmp, \
                mock.patch.dict(os.environ, {"HOME": tmp, "CORE_ON_HIVE": "0"}):
            h = cs.token_help()
        self.assertNotIn("not documented", h)
        for s in ("~/.coreomics_token", '"API Key"', "read -rs TOK", "on THIS computer", "check",
                  "coreomics_key_to_hive.sh"):
            self.assertIn(s, h)
        self.assertNotIn("HIVE has none", h)


# ------------------------------------------------------- the CoreOmics key ON HIVE --
@unittest.skipIf(os.name == "nt", "POSIX file modes")
class TestKeyOnHive(unittest.TestCase):
    """Michelle (2026-10-08, Windows, no Python) saved her key on her laptop, where nothing could
    run the CoreOmics steps, and HIVE had no key. Now the steps run on HIVE from ~/.coreomics_token
    -- in a home folder that is drwxrwsr-x there, open to every account, so the file's own mode is
    the key's only protection: a key anyone else can read is refused, and never sent."""

    def hive_env(self, tmp, base_url="http://127.0.0.1:9/server/api", **extra):
        return key_env(tmp, base_url, CORE_ON_HIVE="1", **extra)

    def keyfile(self, tmp, mode):
        path = os.path.join(key_home(tmp, KEY + "\n"), ".coreomics_token")
        os.chmod(path, mode)
        return path

    def check(self, env):
        rc, out, p = run(["check", "--json"], env)
        for part in (KEY, KEY[:20], KEY[20:]):
            self.assertNotIn(part, p.stdout + p.stderr)
        return rc, out, p

    def test_a_key_others_can_read_is_refused_and_never_sent(self):
        for mode in (0o644, 0o640, 0o604, 0o660):
            with self.subTest(mode=oct(mode)), tempfile.TemporaryDirectory() as tmp, \
                    FakeCoreOmics() as fake:
                env = self.hive_env(tmp, fake.base)
                path = self.keyfile(tmp, mode)
                rc, out, p = self.check(env)
                self.assertEqual((rc, out["status"]), (3, "key_mode_wrong"), p.stderr)
                self.assertEqual(fake.requests, [], "a readable key is not used")
                self.assertEqual((out["runs_on"], out["key_mode"]), ("hive", f"{mode:03o}"))
                self.assertIn(f"mode {mode:03o}", out["say"])
                self.assertIn(f"chmod 600 {shlex.quote(path)}", out["fix"])
                self.assertIn("hive_exec.sh 'chmod 600 ~/.coreomics_token'", out["fix"])
                self.assertEqual(stat.S_IMODE(os.stat(path).st_mode), mode, "check changes nothing")

    def test_owner_only_is_used_and_check_says_where_and_its_mode(self):
        for mode in (0o600, 0o400):
            with self.subTest(mode=oct(mode)), tempfile.TemporaryDirectory() as tmp, \
                    FakeCoreOmics() as fake:
                fake.routes[("GET", "/server/api/submissions/")] = \
                    lambda q, b: (200, {"count": 3, "next": None, "results": [record()]})
                env = self.hive_env(tmp, fake.base)
                path = self.keyfile(tmp, mode)
                rc, out, p = self.check(env)
                self.assertEqual((rc, out["status"]), (0, "ok"), p.stderr)
                self.assertEqual(out["key_source"], path)
                self.assertEqual((out["runs_on"], out["key_mode"]), ("hive", f"{mode:03o}"))
                self.assertIn(f"{path} on HIVE, mode {mode:03o}", out["say"])
                self.assertEqual(fake.requests[0]["auth"], "Token " + KEY)

    def test_fetch_and_identify_refuse_it_too(self):
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            env = self.hive_env(tmp, fake.base)
            self.keyfile(tmp, 0o644)
            rc, out, p = run(["fetch", "807", "--out", os.path.join(tmp, "o")], env)
            self.assertEqual((rc, out["diagnosis"]), (3, "key_mode_wrong"), p.stderr)
            self.assertEqual(out["error"], "no usable CoreOmics API key on HIVE")
            rc, out, p = run(["identify", "/x/run_KG1.d"], env)
            # identify treats an unusable key as none, and its token_help points at check
            self.assertEqual((rc, out["status"]), (3, "needs_token"), p.stderr)
            self.assertIn("core_submission.py check", out["token_help"])
            self.assertEqual(fake.requests, [])
            self.assertNotIn(KEY, p.stdout + p.stderr)

    def test_a_laptop_is_not_held_to_it(self):
        """Off HIVE the home is the user's own (and Windows' chmod does not reach NTFS)."""
        with tempfile.TemporaryDirectory() as tmp, FakeCoreOmics() as fake:
            fake.routes[("GET", "/server/api/submissions/")] = \
                lambda q, b: (200, {"count": 3, "next": None, "results": [record()]})
            self.keyfile(tmp, 0o644)
            rc, out, p = self.check(key_env(tmp, fake.base))
            self.assertEqual((rc, out["status"], out["runs_on"]), (0, "ok", "this_computer"), p.stderr)
            self.assertIn("on this computer, mode 644", out["say"])
        with tempfile.TemporaryDirectory() as tmp:
            path = self.keyfile(tmp, 0o644)
            with mock.patch.dict(os.environ, {ON_HIVE: "0"}):
                self.assertIsNone(cs.key_mode_problem(path))
            with mock.patch.dict(os.environ, {ON_HIVE: "1"}):
                self.assertIn("mode 644", cs.key_mode_problem(path))

    def test_on_hive_is_quobyte_unless_forced(self):
        with mock.patch.dict(os.environ, {}, clear=False):
            os.environ.pop(ON_HIVE, None)
            with mock.patch.object(cs.os.path, "isdir", lambda p: p == "/quobyte"):
                self.assertTrue(cs.on_hive())
            with mock.patch.object(cs.os.path, "isdir", lambda p: False):
                self.assertFalse(cs.on_hive())
                with mock.patch.dict(os.environ, {ON_HIVE: "1"}):
                    self.assertTrue(cs.on_hive())

    def test_no_key_on_hive_names_the_agent_s_two_routes_not_a_typed_line(self):
        """Nobody types anything on HIVE: the agent sends a laptop key, or the user pastes it at
        coreomics_key_to_hive.sh's hidden prompt in their own Git Bash."""
        with tempfile.TemporaryDirectory() as tmp:
            rc, out, _ = self.check(self.hive_env(tmp))
            self.assertEqual((rc, out["status"]), (3, "no_key"))
            self.assertIn("no CoreOmics API key on HIVE", out["say"])
            for s in ("Profile", '"API Key"', "bash scripts/coreomics_key_to_hive.sh`",
                      "coreomics_key_to_hive.sh --paste", "mode 600", "Never paste the key into the chat"):
                self.assertIn(s, out["fix"])
            self.assertNotIn("read -rs TOK", out["fix"])
            rc, out, _ = run(["identify", "/x/run_KG1.d"], self.hive_env(tmp))
            self.assertEqual(out["status"], "needs_token")
            self.assertIn("on HIVE, where this ran", out["token_help"])

    def test_windows_is_told_the_agent_puts_the_key_on_hive(self):
        """The windows branch of the save steps (a Windows Python): the same one-line save, and
        that in hive_remote with no usable Python the agent moves it to HIVE, mode 600."""
        with mock.patch.dict(os.environ, dict(WIN_ENV, HOME="H:\\"), clear=True), \
                mock.patch.object(cs.os, "name", "nt"), mock.patch.object(cs.os, "path", ntpath):
            steps = cs.save_key_steps()
            self.assertFalse(cs.on_hive())
        self.assertIn(cs.SAVE_KEY_LINE["windows"], steps)
        for s in ("Microsoft Store stub", "bash scripts/coreomics_key_to_hive.sh", "HIVE",
                  "~/.coreomics_token", "mode 600", "never shown"):
            self.assertIn(s, steps)
        self.assertNotIn("HIVE has no CoreOmics key", steps)


ON_HIVE = "CORE_ON_HIVE"


@unittest.skipIf(os.name == "nt", "POSIX file modes")
class TestFetchOnHive(unittest.TestCase):
    """fetch run ON HIVE (a laptop with no Python): the record holds emails, phones and PPMS
    billing fields, and HIVE home folders are open to every account. The full record goes into
    <out>/private/ (700, files 600); <out>/ gets exactly what the laptop route --puts there."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.tmp = self._tmp.name
        self.fake = FakeCoreOmics().__enter__()
        rec = record()
        rec["payment"] = {"display": {"PPMS Order Ref #": "PPMS-CANARY"}}
        self.rec = rec
        TestFetch.serve(self, self.fake, rec)
        self.out_dir = os.path.join(self.tmp, "core", "PROT_0807")

    def tearDown(self):
        self.fake.__exit__(None, None, None)
        self._tmp.cleanup()

    def fetch(self):
        rc, out, p = run(["fetch", "807", "--out", self.out_dir],
                         child_env(self.tmp, self.fake.base, CORE_ON_HIVE="1"))
        self.assertEqual(rc, 0, p.stderr)
        return out

    def test_the_full_record_is_private_and_the_top_holds_only_the_hive_copies(self):
        out = self.fetch()
        priv = os.path.join(self.out_dir, "private")
        self.assertEqual(out["ran"], "on HIVE")
        self.assertEqual(stat.S_IMODE(os.stat(priv).st_mode), 0o700)
        for name in ("submission.json", "submission_summary.json"):
            p = os.path.join(priv, name)
            self.assertEqual(stat.S_IMODE(os.stat(p).st_mode), 0o600, name)
            self.assertIn("ada.example@example.org", read(p), name)
        self.assertIn("PPMS-CANARY", read(os.path.join(priv, "submission.json")))
        self.assertEqual(out["outputs"]["summary"], os.path.join(priv, "submission_summary.json"))
        # what locate / stage / attach read: the HIVE copies, where the laptop route --puts them
        self.assertEqual(out["outputs"]["hive_summary"], os.path.join(self.out_dir, "submission_summary.json"))
        self.assertEqual(out["outputs"]["hive_record"], os.path.join(self.out_dir, "submission.json"))
        for name in ("submission_summary.json", "submission.json"):
            text = read(os.path.join(self.out_dir, name))
            for bad in ("@", "PPMS-CANARY", "Robin"):
                self.assertNotIn(bad, text, f"{bad} in {name}")
        self.assertTrue(load(os.path.join(self.out_dir, "submission_summary.json"))["redacted_for_hive"])
        self.assertFalse(os.path.exists(os.path.join(self.out_dir, "hive")), "no hive/ on HIVE")

    def test_a_second_fetch_keeps_it_private_and_a_planted_link_is_refused(self):
        self.fetch()
        priv = os.path.join(self.out_dir, "private")
        os.chmod(priv, 0o755)
        os.chmod(os.path.join(priv, "submission.json"), 0o644)
        self.fetch()
        self.assertEqual(stat.S_IMODE(os.stat(priv).st_mode), 0o700)
        self.assertEqual(stat.S_IMODE(os.stat(os.path.join(priv, "submission.json")).st_mode), 0o600)
        shutil.rmtree(priv)
        elsewhere = os.path.join(self.tmp, "someone_elses")
        os.makedirs(elsewhere)
        os.symlink(elsewhere, priv)
        rc, out, p = run(["fetch", "807", "--out", self.out_dir],
                         child_env(self.tmp, self.fake.base, CORE_ON_HIVE="1"))
        self.assertEqual(rc, 2, p.stderr)
        self.assertIn("not a folder of this account's", out["error"])
        self.assertEqual(os.listdir(elsewhere), [], "nothing was written through the link")

    def test_email_draft_and_bioshare_find_the_full_summary_from_the_hive_copy(self):
        self.fetch()
        hive_copy = os.path.join(self.out_dir, "submission_summary.json")
        draft = os.path.join(self.out_dir, "EMAIL_DRAFT.md")
        env = child_env(self.tmp, self.fake.base, CORE_ON_HIVE="1")
        rc, js, p = run(["email-draft", "--summary", hive_copy, "--share-url",
                         "https://bioshare.example/bsx/", "--out", draft], env)
        self.assertEqual(rc, 0, p.stderr)
        self.assertIn("To: ada.example@example.org", read(draft))
        self.assertEqual(stat.S_IMODE(os.stat(draft).st_mode), 0o600, "names and emails, on HIVE")
        self.assertIn("using the full summary", p.stderr)
        linked = {"id": 41, "submission": HEX, "link_to_path": SERVER_SHARE + "/",
                  "url": "https://bioshare.example/bsx/"}
        self.fake.routes[("GET", SHARES)] = lambda q, b: (200, [linked])
        rc, js, p = run(["bioshare", "send", "--summary", hive_copy], env)     # a dry run
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual(self.fake.posts(), [])
        rec = js["recipients"]
        self.assertEqual(rec["submitter"]["email"], "ada.example@example.org")
        self.assertEqual(rec["pi"]["email"], "pi.placeholder@example.org")
        self.assertEqual([c["email"] for c in rec["contacts"]], ["robin.contact@example.org"])
        rc, js, p = run(["bioshare", "status", "--summary", hive_copy], env)
        self.assertEqual(rc, 0, p.stderr)
        self.assertEqual(js["url"], "https://bioshare.example/bsx/")


if __name__ == "__main__":
    unittest.main()
