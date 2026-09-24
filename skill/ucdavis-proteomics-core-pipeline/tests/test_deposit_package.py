#!/usr/bin/env python3
"""finalize must leave every session with publication Methods and a repository-deposit package.

A paper needs two things the analysis alone does not give: a Methods section that names the
search engine, its version, the parameters, the database and the DE -- and the raw data in a
public ProteomeXchange repository. session.py finalize now ensures output/methods.md (+ .docx)
and writes output/DATA_SUBMISSION/ (HOW_TO_SUBMIT, a pre-filled SDRF, protocol texts, the upload
list and a SLURM prep script), recording each part [OK]/[SKIPPED] in MANIFEST.txt.

What these tests pin:
  * a complete synthetic DIA session -> methods + package + MANIFEST, all in the zip
  * the SDRF: required columns, one row per raw file, TO-FILL for what only the user knows,
    never a guessed sex/age/disease, and values the released validator accepts
  * missing inputs (no raw list, raw files unreachable, no DE, no conditions, Sage DDA) ->
    [SKIPPED] lines or TO-FILL cells, never a crash
  * the raw checksum script is written, not run -- and refuses to run outside SLURM
  * the DE paragraph never states a fold-change threshold run_de.R does not apply
stdlib only, no network.
"""
import csv
import hashlib
import json
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile
import unittest
import zipfile
from unittest import mock

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
from job_env import job_env  # noqa: E402  (env for running job scripts)
sys.path.insert(0, HERE)

import make_deposit as md                            # noqa: E402
import make_methods as mm                            # noqa: E402
from synthetic_tdf import synthetic_tdf_write_uri   # noqa: E402  (the deliberate writer)

PY = sys.executable
RUNS = ["HeLa_ctrl_01", "HeLa_ctrl_02", "HeLa_trt_01"]
DIANN_CFG = """--qvalue 0.01
--matrices
--fasta-search
--predictor
--reanalyse
--cut K*,R*
--missed-cleavages 1
--min-pep-len 7
--max-pep-len 30
--min-pr-charge 2
--max-pr-charge 4
--min-pr-mz 299
--max-pr-mz 1201
--met-excision
--unimod4
--var-mods 1
--var-mod UniMod:35,15.994915,M
--mass-acc 15
--mass-acc-ms1 15
"""
DE_PROV = {
    "pipeline_id": "dpc", "display_label": "DPC-Quant + limma (limpa)",
    "rollup_method": "DPC-Quant (Detection Probability Curve quantification, dpcCN)",
    "de_engine": "limpa::dpcDE (voomaLmFitWithImputation) -> contrasts.fit -> eBayes",
    "missing_policy": "Missing precursors modelled via the detection probability curve.",
    "citation": "Li M, Cobbold SA, Smyth GK (2025) bioRxiv 10.1101/2025.04.28.651125",
    "method": "dpc", "q_cutoff": 0.01, "q_columns": ["Q.Value", "Global.PG.Q.Value"],
    "q_cutoffs": [0.01, 0.05], "logfc": 1, "logfc_role": "reference_line_only", "adjp": 0.05,
    "design": "~ 0 + Group", "contrasts": ["Treated-Control"], "R_version": "4.5.1",
    "packages": {"limpa": "1.2.0", "limma": "3.64.0"}}


def make_d(path, instrument="timsTOF HT"):
    """A minimal Bruker .d that make_methods.bruker_meta() can read (dia-PASEF)."""
    os.makedirs(path)
    tdf = os.path.join(path, "analysis.tdf")
    con = sqlite3.connect(synthetic_tdf_write_uri(tdf), uri=True)
    con.execute("CREATE TABLE GlobalMetadata (Key TEXT, Value TEXT)")
    con.executemany("INSERT INTO GlobalMetadata VALUES (?,?)",
                    [("InstrumentName", instrument), ("MzAcqRangeLower", "100.0"),
                     ("MzAcqRangeUpper", "1700.0"), ("OneOverK0AcqRangeLower", "0.6"),
                     ("OneOverK0AcqRangeUpper", "1.4")])
    con.execute("CREATE TABLE Frames (Id INTEGER, MsMsType INTEGER, AccumulationTime REAL, "
                "RampTime REAL)")
    con.executemany("INSERT INTO Frames VALUES (?,?,?,?)",
                    [(i, 9 if i % 4 else 0, 100.0, 100.0) for i in range(12)])
    con.commit()
    con.close()
    with open(os.path.join(path, "analysis.tdf_bin"), "wb") as fh:
        fh.write(b"\0" * 4096)


def init_session(root, raw_glob, name="demo"):
    r = subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "init", "--name", name,
                        "--date", "2026-09-24", "--raw", raw_glob, "--base", root],
                       capture_output=True, text=True, check=True)
    return json.loads(r.stdout)["paths"]


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


def dia_session(root, conditions=True, de=True, fasta=True, raw_reachable=True):
    """A finished DIA-NN session: 3 human diaPASEF runs, two groups."""
    raw = os.path.join(root, "raw")
    os.makedirs(raw)
    for r in RUNS:
        make_d(os.path.join(raw, r + ".d"))
    p = init_session(root, os.path.join(raw, "*.d"))
    if not raw_reachable:
        shutil.rmtree(raw)
    if conditions:
        write(p["conditions"], "File.Name,Group\nHeLa_ctrl_01,Control\nHeLa_ctrl_02,Control\n"
                               "HeLa_trt_01,Treated\n")
    if fasta:
        write(p["fasta"], ">sp|P1|A_HUMAN\nPEPTIDEK\n")
        write(p["fasta_meta"], json.dumps({
            "fasta": p["fasta"], "organism": "Homo sapiens", "taxid": 9606,
            "proteome": "UP000005640", "proteome_type": "Reference proteome",
            "content_used": "one_per_gene", "uniprot_release": "2026_03", "n_proteome": 20659,
            "n_contaminants_appended": 381, "contaminant_set": "universal",
            "diann_cont_quant_exclude": "Cont_"}))
    write(os.path.join(p["workflow_dir"], "params.cfg"), DIANN_CFG)
    write(p["workflow_manifest"], json.dumps({
        "acquisition": "DIA", "instruments": ["timsTOF HT"], "organism_taxid": 9606,
        "engine": {"name": "diann", "version": "2.7.0"}}))
    for n in ("report.parquet", "report.log.txt", "report.pg_matrix.tsv"):
        write(os.path.join(p["search_out"], n), "x" * 100)
    write(os.path.join(p["search_out"], "quant_step2", "a.quant"), "q")
    write(p["search_prov"], json.dumps({
        "engine": "diann", "version": "2.7.0",
        "params_file": os.path.join(p["workflow_dir"], "params.cfg")}))
    if de:
        write(os.path.join(p["de_dir"], "de_provenance.json"), json.dumps(DE_PROV))
        write(os.path.join(p["de_dir"], "DE_Treated_vs_Control.csv"), "Protein,logFC\nP1,1\n")
    return p


def finalize(sd, *extra):
    # job_env: finalize logs the run (record_run.py) and posts "analysis complete" to the Core's
    # Slack channel, and a test must never do either (or ssh to HIVE to relay the post).
    # test_slack_notify.py covers both, against fakes.
    return subprocess.run([PY, os.path.join(SCRIPTS, "session.py"), "finalize", "--dir", sd,
                           *extra], capture_output=True, text=True,
                          env=job_env(os.path.dirname(sd)))


def read_sdrf(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))
    return rows[0], rows[1:]


def col(header, rows, name, nth=0):
    idx = [i for i, h in enumerate(header) if h == name][nth]
    return [r[idx] for r in rows]


def manifest_lines(p):
    with open(p["manifest_txt"]) as fh:
        return [ln for ln in fh.read().splitlines() if ln.startswith("[")]


class ZipLeavesQuantOut(unittest.TestCase):
    """DIA-NN's per-run .quant intermediates (~30 MB each on HIVE) live inside the session since
    the single-shot search's --temp moved to <search out>/quant: ~3 GB of zip for 100 runs."""

    def test_quant_files_and_the_quant_dir_stay_out_of_the_zip(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp)
            for rel in ("quant/HeLa_ctrl_01.quant", "quant/HeLa_ctrl_02.quant",
                        "quant/nested/extra.bin", "quant_step4/HeLa_trt_01.quant"):
                write(os.path.join(p["search_out"], rel), "q" * 64)
            write(os.path.join(p["de_dir"], "quant", "notes.txt"), "not DIA-NN's")  # kept
            r = finalize(p["session_dir"], "--zip")
            self.assertEqual(r.returncode, 0, r.stderr)
            res = json.loads(r.stdout)
            names = zipfile.ZipFile(res["zip"]).namelist()
            self.assertFalse([n for n in names if n.endswith(".quant")], names)
            self.assertFalse([n for n in names if "/output/search/quant/" in n], names)
            self.assertTrue(any(n.endswith("output/search/report.parquet") for n in names))
            self.assertTrue(any(n.endswith("output/tables/quant/notes.txt") for n in names))
            label = "DIA-NN .quant intermediates (kept on disk where they are)"
            # 3 in quant/, plus quant_step4/HeLa_trt_01.quant and dia_session's quant_step2/a.quant
            self.assertEqual(res["zip_excluded"][label], 5)
            for rel in ("quant/HeLa_ctrl_01.quant", "quant_step4/HeLa_trt_01.quant"):
                self.assertTrue(os.path.isfile(os.path.join(p["search_out"], rel)))

    def test_predicted_library_stays_out_but_empirical_libraries_go_in(self):
        # 687 MB for one 15-file Lumos session; rebuilt exactly from the FASTA + params.
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp)
            write(os.path.join(p["search_out"], "step1.predicted.speclib"), "p" * 64)
            write(os.path.join(p["search_out"], "diann_lib.predicted.speclib"), "p" * 64)
            write(os.path.join(p["search_out"], "step3_empirical.parquet"), "e" * 64)
            write(os.path.join(p["search_out"], "custom.speclib"), "c" * 64)
            r = finalize(p["session_dir"], "--zip")
            self.assertEqual(r.returncode, 0, r.stderr)
            res = json.loads(r.stdout)
            names = zipfile.ZipFile(res["zip"]).namelist()
            self.assertFalse([n for n in names if n.endswith(".predicted.speclib")], names)
            self.assertTrue(any(n.endswith("output/search/step3_empirical.parquet") for n in names))
            self.assertTrue(any(n.endswith("output/search/custom.speclib") for n in names))
            label = [k for k in res["zip_excluded"] if k.startswith("predicted spectral libraries")]
            self.assertEqual(len(label), 1, res["zip_excluded"])
            self.assertEqual(res["zip_excluded"][label[0]], 2)
            self.assertTrue(os.path.isfile(os.path.join(p["search_out"], "step1.predicted.speclib")))


class FullSession(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.mkdtemp()
        cls.p = dia_session(cls.tmp)
        cls.res = finalize(cls.p["session_dir"], "--zip")
        cls.out = json.loads(cls.res.stdout) if cls.res.returncode == 0 else {}
        cls.pkg = cls.p["deposit_dir"]

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def test_finalize_succeeds_with_nothing_skipped(self):
        self.assertEqual(self.res.returncode, 0, self.res.stderr)
        lines = manifest_lines(self.p)
        skipped = [ln for ln in lines if ln.startswith("[SKIPPED]")]
        # .docx needs pandoc or python-docx, which a bare test interpreter may lack; that is
        # recorded, never silent. Everything else must be [OK].
        self.assertTrue(all("Word" in ln for ln in skipped), skipped)
        for part in ("Publication methods (output/methods.md)", "sdrf.tsv", "protocols.txt",
                     "files_to_upload.tsv", "prepare_upload.sbatch", "HOW_TO_SUBMIT.md",
                     "HOW_TO_SUBMIT.html"):
            self.assertTrue(any(part in ln for ln in lines), part)

    def test_package_files_exist(self):
        for n in ("HOW_TO_SUBMIT.md", "HOW_TO_SUBMIT.html", "sdrf.tsv", "protocols.txt",
                  "files_to_upload.tsv", "prepare_upload.sbatch"):
            self.assertTrue(os.path.getsize(os.path.join(self.pkg, n)) > 0, n)

    def test_methods_cover_search_database_de_and_acknowledgment(self):
        with open(self.p["methods_md"], encoding="utf-8") as fh:
            text = fh.read()
        secs = md.md_sections(text)
        for s in ("Liquid chromatography", "Mass spectrometry", "Sequence database",
                  "Database search", "Differential expression", "Acknowledgments"):
            self.assertIn(s, secs)
        self.assertIn("DIA-NN 2.7.0", secs["Database search"])
        self.assertIn("Trypsin/P", secs["Database search"])
        self.assertIn("Oxidation (M)", secs["Database search"])
        self.assertIn("1% FDR", secs["Database search"])
        self.assertIn("UniProt Homo sapiens reference proteome", secs["Sequence database"])
        self.assertIn("Howard Hughes Medical Institute", secs["Acknowledgments"])
        self.assertIn("[facility default — confirm]", text)       # defaults stay tagged

    def test_de_paragraph_states_no_fold_change_filter(self):
        with open(self.p["methods_md"], encoding="utf-8") as fh:
            de = md.md_sections(fh.read())["Differential expression"]
        self.assertIn("No fold-change filter was applied", de)
        self.assertNotRegex(de, r"\|log2FC\|\s*[≥>]")
        self.assertIn("Global.PG.Q.Value ≤ 0.05", de)

    def test_sdrf_rows_columns_and_to_fill(self):
        header, rows = read_sdrf(os.path.join(self.pkg, "sdrf.tsv"))
        self.assertEqual(len(rows), 3)
        for req in ("source name", "characteristics[organism]", "characteristics[organism part]",
                    "characteristics[biological replicate]", "assay name", "technology type",
                    "comment[proteomics data acquisition method]", "comment[label]",
                    "comment[instrument]", "comment[cleavage agent details]",
                    "comment[fraction identifier]", "comment[technical replicate]",
                    "comment[data file]", "characteristics[disease]", "characteristics[sex]",
                    "characteristics[age]"):
            self.assertIn(req, header)
        self.assertEqual(header[0], "source name")
        self.assertTrue(header[-1].startswith("factor value["))
        self.assertEqual(col(header, rows, "comment[data file]"),
                         [r + ".d" for r in RUNS])
        self.assertEqual(set(col(header, rows, "characteristics[organism]")), {"homo sapiens"})
        for unknown in ("characteristics[organism part]", "characteristics[disease]",
                        "characteristics[cell type]", "characteristics[sex]",
                        "characteristics[age]"):
            self.assertEqual(set(col(header, rows, unknown)), {md.TO_FILL}, unknown)
        self.assertEqual(col(header, rows, "characteristics[biological replicate]"),
                         ["1", "2", "1"])
        self.assertEqual(col(header, rows, header[-1]), ["Control", "Control", "Treated"])
        self.assertEqual(header[-1], f"factor value[{md.TO_FILL}]")
        # the one value the released validator accepts in a dia-acquisition SDRF
        self.assertEqual(set(col(header, rows, "comment[proteomics data acquisition method]")),
                         {"Data-independent acquisition"})
        self.assertEqual(set(col(header, rows, "comment[instrument]")),
                         {"NT=timsTOF HT;AC=MS:1003404"})
        self.assertEqual(set(col(header, rows, "comment[cleavage agent details]")),
                         {"NT=Trypsin/P;AC=MS:1001313"})
        self.assertEqual(set(col(header, rows, "comment[label]")), {"label free sample"})
        mods = {col(header, rows, "comment[modification parameters]", i)[0] for i in (0, 1)}
        self.assertEqual(mods, {"NT=Carbamidomethyl;AC=UNIMOD:4;TA=C;MT=fixed;PP=Anywhere",
                                "NT=Oxidation;AC=UNIMOD:35;TA=M;MT=variable;PP=Anywhere"})
        self.assertEqual(col(header, rows, "comment[sdrf template]", 0)[0], "human v1.1.0")
        self.assertEqual(col(header, rows, "comment[sdrf template]", 1)[0],
                         "dia-acquisition v1.1.0")
        # reserved words are the user's to choose, never written as a stand-in
        flat = "\t".join("\t".join(r) for r in rows)
        self.assertNotIn("not available", flat)
        self.assertNotIn("female", flat)

    def test_upload_list(self):
        with open(os.path.join(self.pkg, "files_to_upload.tsv"), newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        by = {r["upload_name"]: r for r in rows}
        for r in RUNS:
            raw = by[r + ".d.tar.gz"]
            self.assertEqual(raw["pride_file_type"], "RAW")
            self.assertIn("prepare_upload.sbatch", raw["md5"])      # never hashed inline
            self.assertIn("compress", raw["before_upload"])
            self.assertTrue(raw["source_path"].endswith(r + ".d"))
        self.assertEqual(by["report.parquet"]["pride_file_type"], "SEARCH")
        self.assertEqual(by["report.log.txt"]["requirement"], "required")
        self.assertEqual(by["sdrf.tsv"]["pride_file_type"], "EXPERIMENTAL_DESIGN")
        self.assertEqual(by["search.fasta"]["pride_file_type"], "FASTA")
        self.assertEqual(by["report.parquet"]["md5"],
                         hashlib.md5(b"x" * 100).hexdigest())
        self.assertFalse(any(n.endswith(".quant") for n in by), "intermediates are not uploads")

    def test_prep_script_is_written_not_run(self):
        script = os.path.join(self.pkg, "prepare_upload.sbatch")
        self.assertEqual(subprocess.run(["bash", "-n", script]).returncode, 0)
        self.assertFalse(os.path.exists(os.path.join(self.pkg, "upload_staging")))
        with open(script) as fh:
            text = fh.read()
        self.assertIn("#SBATCH --job-name=deposit_prep", text)
        for r in RUNS:
            self.assertIn(f"{r}.d\t{r}.d.tar.gz", text)

    def test_prep_script_refuses_outside_slurm(self):
        env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
        env.pop("RUN_HERE", None)
        # job_env: not a search job (the deposit prep script has no job-end hook)
        r = subprocess.run(["bash", os.path.join(self.pkg, "prepare_upload.sbatch")],
                           capture_output=True, text=True, env=env)
        self.assertEqual(r.returncode, 2)
        self.assertIn("sbatch", r.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.pkg, "upload_staging")))

    def test_zip_carries_methods_package_and_manifest(self):
        self.assertTrue(os.path.isfile(self.out.get("zip", "")))
        base = os.path.basename(self.p["session_dir"])
        names = set(zipfile.ZipFile(self.out["zip"]).namelist())
        for rel in ("MANIFEST.txt", "README.md", "output/methods.md",
                    "output/DATA_SUBMISSION/HOW_TO_SUBMIT.md",
                    "output/DATA_SUBMISSION/HOW_TO_SUBMIT.html",
                    "output/DATA_SUBMISSION/sdrf.tsv", "output/DATA_SUBMISSION/protocols.txt",
                    "output/DATA_SUBMISSION/files_to_upload.tsv",
                    "output/DATA_SUBMISSION/prepare_upload.sbatch"):
            self.assertIn(f"{base}/{rel}", names)

    def test_readme_points_to_the_package(self):
        with open(self.p["readme"]) as fh:
            text = fh.read()
        self.assertIn("DATA_SUBMISSION", text)
        self.assertIn("HOW_TO_SUBMIT.md", text)
        self.assertIn("output/methods.md", text)

    def test_protocols_within_pride_limits(self):
        with open(os.path.join(self.pkg, "protocols.txt"), encoding="utf-8") as fh:
            text = fh.read()
        self.assertIn("SAMPLE PROCESSING PROTOCOL", text)
        self.assertIn("DATA PROCESSING PROTOCOL", text)
        self.assertIn(md.TO_FILL, text)                      # sample prep is the user's
        sample, data = md.build_protocols(self.p["methods_md"])
        for block in (sample, data):
            self.assertTrue(md.PROTOCOL_MIN < len(block) < md.PROTOCOL_MAX, len(block))
        self.assertIn("DIA-NN 2.7.0", data)

    def test_howto_html_is_well_formed(self):
        from html.parser import HTMLParser

        class P(HTMLParser):
            def __init__(self):
                super().__init__()
                self.stack, self.bad = [], 0

            def handle_starttag(self, tag, attrs):
                if tag not in ("meta", "br"):
                    self.stack.append(tag)

            def handle_endtag(self, tag):
                if self.stack and self.stack[-1] == tag:
                    self.stack.pop()
                else:
                    self.bad += 1

        p = P()
        with open(os.path.join(self.pkg, "HOW_TO_SUBMIT.html"), encoding="utf-8") as fh:
            p.feed(fh.read())
        self.assertEqual((p.stack, p.bad), ([], 0))


class RunningThePrepScript(unittest.TestCase):
    """RUN_HERE=1 (a laptop with local raws): archives, staging and PRIDE-format checksums."""

    def test_stage_and_checksum(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp)
            man = md.Manifest()
            md.build(p["session_dir"], man)
            env = {k: v for k, v in os.environ.items() if not k.startswith("SLURM_")}
            env["RUN_HERE"] = "1"
            script = os.path.join(p["deposit_dir"], "prepare_upload.sbatch")
            # job_env: not a search job (the deposit prep script has no job-end hook)
            r = subprocess.run(["bash", script], capture_output=True, text=True, env=env)
            self.assertEqual(r.returncode, 0, r.stderr)
            stage = os.path.join(p["deposit_dir"], "upload_staging")
            with open(os.path.join(stage, "checksum.txt")) as fh:
                sums = dict(ln.rstrip("\n").split("\t") for ln in fh)
            for name, sha in sums.items():
                with open(os.path.join(stage, name), "rb") as fh:
                    self.assertEqual(hashlib.sha1(fh.read()).hexdigest(), sha, name)
            self.assertIn("HeLa_trt_01.d.tar.gz", sums)
            listing = subprocess.run(["tar", "-tzf", os.path.join(stage, "HeLa_trt_01.d.tar.gz")],
                                     capture_output=True, text=True).stdout.split()
            self.assertIn("HeLa_trt_01.d/analysis.tdf", listing)
            # the session zip leaves the staged raw archives out
            r = finalize(p["session_dir"], "--zip")
            names = zipfile.ZipFile(json.loads(r.stdout)["zip"]).namelist()
            self.assertFalse(any("upload_staging" in n for n in names))


class MissingInputs(unittest.TestCase):
    """Parts that cannot be made are [SKIPPED] with a reason, or TO-FILL -- never a crash."""

    def test_no_raw_list_no_conditions(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = init_session(tmp, os.path.join(tmp, "nothing", "*.d"))
            os.remove(p["raw_list"]) if os.path.exists(p["raw_list"]) else None
            r = finalize(p["session_dir"], "--zip")
            self.assertEqual(r.returncode, 0, r.stderr)
            lines = manifest_lines(p)
            text = "\n".join(lines)
            self.assertRegex(text, r"\[SKIPPED\] Publication methods \(output/methods\.md\)\s+"
                                   r"-- .*raw file list")
            self.assertRegex(text, r"\[SKIPPED\] Deposit: sdrf\.tsv.* -- .*nothing to describe")
            self.assertRegex(text, r"\[SKIPPED\] Deposit: raw files in the upload plan")
            self.assertRegex(text, r"\[SKIPPED\] Deposit: protocols\.txt.* -- .*methods")
            with open(p["readme"]) as fh:
                self.assertIn("MANIFEST.txt", fh.read())

    def test_raw_files_unreachable_methods_from_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp, raw_reachable=False)
            r = finalize(p["session_dir"])
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(p["methods_md"], encoding="utf-8") as fh:
                text = fh.read()
            self.assertIn("could not be read from where this was run", text)
            self.assertIn("timsTOF HT", text)
            self.assertIn("[raw file not readable here — confirm]", text)
            with open(os.path.join(p["deposit_dir"], "files_to_upload.tsv"), newline="") as fh:
                raws = [r for r in csv.DictReader(fh, delimiter="\t")
                        if r["pride_file_type"] == "RAW"]
            self.assertEqual(len(raws), 3)
            self.assertTrue(all(r["size_bytes"] == "" for r in raws))
            self.assertIn("not reachable", raws[0]["notes"])

    def test_qc_only_no_de_no_conditions(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp, conditions=False, de=False)
            r = finalize(p["session_dir"])
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(p["methods_md"], encoding="utf-8") as fh:
                secs = md.md_sections(fh.read())
            self.assertNotIn("Differential expression", secs)
            self.assertIn("Database search", secs)
            self.assertFalse(any("lacks" in ln for ln in manifest_lines(p)))
            header, rows = read_sdrf(os.path.join(p["deposit_dir"], "sdrf.tsv"))
            self.assertEqual(set(col(header, rows, "characteristics[biological replicate]")),
                             {md.TO_FILL})
            self.assertEqual(set(col(header, rows, header[-1])), {md.TO_FILL})

    def test_no_fasta_meta_organism_to_fill(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp, fasta=False)
            with open(p["workflow_manifest"]) as fh:
                wf = json.load(fh)
            wf.pop("organism_taxid")
            write(p["workflow_manifest"], json.dumps(wf))
            self.assertEqual(finalize(p["session_dir"]).returncode, 0)
            header, rows = read_sdrf(os.path.join(p["deposit_dir"], "sdrf.tsv"))
            self.assertEqual(set(col(header, rows, "characteristics[organism]")), {md.TO_FILL})
            self.assertNotIn("characteristics[sex]", header)   # not known to be human

    def test_existing_methods_kept_and_draft_written(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp)
            write(p["methods_md"], "# My polished methods\n\n## Liquid chromatography\n\nX\n")
            self.assertEqual(finalize(p["session_dir"]).returncode, 0)
            with open(p["methods_md"]) as fh:
                self.assertIn("My polished methods", fh.read())
            draft = os.path.join(p["output_dir"], "methods_complete_draft.md")
            self.assertTrue(os.path.isfile(draft))
            self.assertTrue(any("lacks: Mass spectrometry" in ln for ln in manifest_lines(p)))

    def test_no_deposit_flag(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp)
            self.assertEqual(finalize(p["session_dir"], "--no-deposit").returncode, 0)
            self.assertTrue(os.path.isfile(p["methods_md"]))
            self.assertFalse(os.path.exists(os.path.join(p["deposit_dir"], "sdrf.tsv")))
            self.assertTrue(any(ln.startswith("[SKIPPED] Deposit package") and "--no-deposit"
                                in ln for ln in manifest_lines(p)))


    def test_human_from_taxid_alone(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp, fasta=False)        # the manifest still records taxid 9606
            self.assertEqual(finalize(p["session_dir"]).returncode, 0)
            header, rows = read_sdrf(os.path.join(p["deposit_dir"], "sdrf.tsv"))
            self.assertEqual(set(col(header, rows, "characteristics[organism]")),
                             {"homo sapiens"})
            self.assertIn("characteristics[sex]", header)


class SageDDA(unittest.TestCase):
    def test_sage_dda_session(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw = os.path.join(tmp, "raw")
            os.makedirs(raw)
            for r in ("Ex_A1", "Ex_B1"):
                write(os.path.join(raw, r + ".raw"), "thermo")
            p = init_session(tmp, os.path.join(raw, "*.raw"))
            write(p["conditions"], "File.Name,Group\nEx_A1,A\nEx_B1,B\n")
            write(os.path.join(p["workflow_dir"], "params.json"), json.dumps({
                "database": {"enzyme": {"missed_cleavages": 2, "min_len": 7, "max_len": 30,
                                        "cleave_at": "KR", "restrict": "P"},
                             "static_mods": {"C": 57.0215},
                             "variable_mods": {"[": [42.0106], "M": [15.9949]},
                             "max_variable_mods": 2},
                "precursor_tol": {"ppm": [-10.0, 10.0]}, "fragment_tol": {"ppm": [-10.0, 10.0]},
                "quant": {"lfq": True}}))
            write(p["workflow_manifest"], json.dumps({
                "acquisition": "DDA", "instruments": ["Orbitrap Exploris 480"],
                "engine": {"name": "sage", "version": "0.14.7"}}))
            write(os.path.join(p["search_out"], "results.sage.tsv"), "psm\n")
            write(os.path.join(p["search_out"], "lfq.tsv"), "lfq\n")
            write(p["search_prov"], json.dumps({"engine": "sage", "version": "0.14.7"}))
            r = finalize(p["session_dir"])
            self.assertEqual(r.returncode, 0, r.stderr)
            header, rows = read_sdrf(os.path.join(p["deposit_dir"], "sdrf.tsv"))
            self.assertEqual(set(col(header, rows, "comment[proteomics data acquisition "
                                                   "method]")), {"Data-dependent acquisition"})
            self.assertEqual(set(col(header, rows, "comment[cleavage agent details]")),
                             {"NT=Trypsin;AC=MS:1001251"})
            self.assertEqual(set(col(header, rows, "comment[instrument]")),
                             {"NT=Orbitrap Exploris 480;AC=MS:1003028"})
            mods = {v for i in range(3)
                    for v in col(header, rows, "comment[modification parameters]", i)}
            self.assertIn("NT=Acetyl;AC=UNIMOD:1;MT=variable;PP=Protein N-term", mods)
            self.assertEqual(col(header, rows, "comment[precursor mass tolerance]")[0], "10 ppm")
            self.assertEqual(col(header, rows, "comment[data file]"), ["Ex_A1.raw", "Ex_B1.raw"])
            self.assertIn("ms-proteomics v1.1.0", col(header, rows, "comment[sdrf template]")[0])
            with open(p["methods_md"], encoding="utf-8") as fh:
                search = md.md_sections(fh.read())["Database search"]
            self.assertIn("Sage 0.14.7", search)
            self.assertIn("±10 ppm", search)
            with open(os.path.join(p["deposit_dir"], "files_to_upload.tsv"), newline="") as fh:
                by = {r["upload_name"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(by["Ex_A1.raw"]["pride_file_type"], "RAW")
            self.assertEqual(by["Ex_A1.raw"]["before_upload"], "")   # a file: no compression
            self.assertEqual(by["results.sage.tsv"]["pride_file_type"], "SEARCH")


class FragPipe(unittest.TestCase):
    """An engine whose parameter file is not parsed: named and versioned, the rest TO-FILL."""

    def test_fragpipe_session(self):
        with tempfile.TemporaryDirectory() as tmp:
            raw = os.path.join(tmp, "raw")
            os.makedirs(raw)
            for r in ("S1", "S2"):
                make_d(os.path.join(raw, r + ".d"), instrument="timsTOF Pro 2")
            p = init_session(tmp, os.path.join(raw, "*.d"))
            write(os.path.join(p["workflow_dir"], "fragpipe.workflow"),
                  "msfragger.search_enzyme_name_1=stricttrypsin\n")
            write(p["workflow_manifest"], json.dumps({
                "acquisition": "DIA", "instruments": ["timsTOF Pro 2"],
                "engine": {"name": "fragpipe", "version": "24.0"}}))
            for rel in ("exp1/psm.tsv", "exp1/protein.tsv", "combined_protein.tsv",
                        "fragpipe-files.fp-manifest"):
                write(os.path.join(p["search_out"], rel), "x\n")
            r = finalize(p["session_dir"])
            self.assertEqual(r.returncode, 0, r.stderr)
            with open(p["methods_md"], encoding="utf-8") as fh:
                search = md.md_sections(fh.read())["Database search"]
            self.assertIn("FragPipe 24.0 [pinned version — confirm it is what ran]", search)
            self.assertIn("not parsed here", search)
            header, rows = read_sdrf(os.path.join(p["deposit_dir"], "sdrf.tsv"))
            for c in ("comment[cleavage agent details]", "comment[label]"):
                self.assertEqual(set(col(header, rows, c)), {md.TO_FILL}, c)
            self.assertNotIn("comment[modification parameters]", header)
            self.assertEqual(set(col(header, rows, "comment[instrument]")),
                             {"NT=timsTOF Pro 2;AC=MS:1003230"})
            with open(os.path.join(p["deposit_dir"], "files_to_upload.tsv"), newline="") as fh:
                by = {r["upload_name"]: r for r in csv.DictReader(fh, delimiter="\t")}
            self.assertEqual(by["exp1_psm.tsv"]["pride_file_type"], "SEARCH")
            self.assertEqual(by["fragpipe.workflow"]["requirement"], "required")
            self.assertEqual(by["fragpipe-files.fp-manifest"]["requirement"], "required")


class Units(unittest.TestCase):
    def test_manifest_mirrors_safe_section(self):
        man = md.Manifest()
        self.assertTrue(man.section("fine", lambda: "a note"))
        self.assertFalse(man.section("clean skip", lambda: (_ for _ in ()).throw(
            md.Skip("because"))))
        self.assertFalse(man.section("crash", lambda: 1 / 0))
        self.assertRegex(man.lines[0], r"^\[OK\]\s+fine\s+\(\d+\.\ds\) -- a note$")
        self.assertRegex(man.lines[1], r"^\[SKIPPED\] clean skip\s+-- because$")
        self.assertRegex(man.lines[2], r"^\[SKIPPED\] crash\s+-- ZeroDivisionError")
        self.assertEqual(man.n_skipped, 2)

    def test_raw_data_never_hashed_inline(self):
        with tempfile.TemporaryDirectory() as tmp:
            p = dia_session(tmp)
            seen = []
            real = md._hash2
            with mock.patch.object(md, "_hash2", side_effect=lambda f: seen.append(f) or real(f)):
                md.build(p["session_dir"], md.Manifest())
            raw_dir = os.path.join(tmp, "raw")
            self.assertTrue(seen)
            self.assertFalse([f for f in seen if f.startswith(raw_dir)], seen)

    def test_large_session_file_not_hashed_inline(self):
        rows = [{"kind": "copy", "source_path": __file__, "upload_name": "x.tsv",
                 "size_bytes": md.HASH_MAX_BYTES + 1, "md5": "", "sha1": ""}]
        md.hash_small(rows)
        self.assertIn("prepare_upload.sbatch", rows[0]["md5"])

    def test_unknown_instrument_and_enzyme_are_to_fill(self):
        self.assertIsNone(md._instrument_term("Some Future Instrument"))
        self.assertEqual(md._instrument_term("Thermo Orbitrap Fusion Lumos"),
                         "NT=Orbitrap Fusion Lumos;AC=MS:1002732")
        m = mm._diann_mod(["UniMod:999,1.0,K"], "variable", "cfg")
        self.assertIsNone(md._mod_term(m))

    def test_safe_name(self):
        self.assertEqual(md.safe_name("HeLa 01 (rep).d"), "HeLa_01__rep_.d")
        self.assertTrue(md.PRIDE_NAME_OK.match(md.safe_name("#x.raw")))

    def test_de_paragraph_without_logfc_role_is_not_a_claim(self):
        text = mm.de_paragraph({"display_label": "X", "adjp": 0.05, "logfc": 1})
        self.assertIn(mm.NOT_RECORDED, text)
        self.assertNotIn("No fold-change filter", text)
        self.assertNotRegex(text, r"\|log2FC\|\s*[≥>]")


if __name__ == "__main__":
    unittest.main()
