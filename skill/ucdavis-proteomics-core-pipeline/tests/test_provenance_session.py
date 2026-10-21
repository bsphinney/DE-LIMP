#!/usr/bin/env python3
"""
provenance.py --session <S>: the bundle's arguments from the session's own records.

SKILL.md step 10 had the agent hand-type a 12-line provenance.py command, so a bundle could
quietly get the wrong inputs. With --session each argument not given explicitly is read from
where the session keeps it: the workflow manifest, params, conditions, FASTA + sidecar, the raw
files, the report and the DE from input/ and output/; the engine, its build and tools.json from
search_provenance.json; the DE's settings from de_provenance.json (never the workflow's
defaults); acquisition / instrument from the workflow manifest or acquisition.json; the taxon
from the FASTA's sidecar; the QC record only when there is one. An explicit argument wins; a
value with no record stays unset and is listed under from_session.not_found -- never guessed. A
--session that is not a readable session exits 2.

A small fake session tree; hermetic (no R run: the DE's sessionInfo.txt is there).
"""
import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text if isinstance(text, str) else json.dumps(text, indent=2))
    return path


class SessionFixture(unittest.TestCase):
    def setUp(self):
        self.tmp = os.path.realpath(tempfile.mkdtemp())
        self.addCleanup(shutil.rmtree, self.tmp, True)
        self.S = S = os.path.join(self.tmp, "2026-10-20_PROT_0807")
        # the raw data, outside the session (init records the paths, never copies)
        self.raws = []
        for n in ("A1", "B1"):
            d = os.path.join(self.tmp, "raw", f"{n}.d")
            write(os.path.join(d, "analysis.tdf"), f"tdf {n}")
            self.raws.append(d)
        self.tools = write(os.path.join(self.tmp, "tools", "tools.json"),
                           {"versions": {"diann": "2.3.0"}, "diann": "/opt/diann/diann-linux"})
        self.setup_json = write(os.path.join(self.tmp, "setup.json"), {})
        for sub in ("output/search", "output/tables", "output/reproducibility", "logs",
                    "input/wf"):
            os.makedirs(os.path.join(S, sub), exist_ok=True)
        self.wf = write(os.path.join(S, "input/wf/workflow.manifest.json"), {
            "id": "diann_orbitrap_dia", "name": "diann (DIA)", "acquisition": "DIA",
            "instruments": ["Orbitrap Astral"], "organism_taxid": 10090,
            "engine": {"name": "diann", "version": "2.3.0"},
            # resolve_defaults.py's DEFAULTS: never a record of what the DE ran
            "de": {"method": "dpc", "q_cutoff": 0.01, "logfc": 1.0, "adjp": 0.05},
            "fasta": {"add_contaminants": True},
            "registry": {"source": "local skill defaults", "defaults_version": "test"}})
        self.params = write(os.path.join(S, "input/wf/params.cfg"), "--qvalue 0.01\n")
        self.conditions = write(os.path.join(S, "input/conditions.csv"),
                                "Run,Group\nA1,A\nB1,B\n")
        self.fasta = write(os.path.join(S, "input/search.fasta"), ">sp|P1|X_HUMAN\nPEPTIDEK\n")
        self.sidecar = write(self.fasta + ".meta.json", {
            "organism": "Homo sapiens", "taxid": 9606, "proteome": "UP000005640",
            "content_used": "one_per_gene", "n_proteome": 1, "contaminant_set": "none"})
        write(os.path.join(S, "input/raw_files.txt"),
              "# Raw MS files used in this analysis\n" + "".join(r + "\n" for r in self.raws))
        write(os.path.join(S, "input/acquisition.json"),
              {"overall": "DIA", "instrument": "Orbitrap Exploris 480",
               "instruments_seen": ["Orbitrap Exploris 480"], "files": []})
        self.report = write(os.path.join(S, "output/search/report.parquet"), "PAR1 fake report")
        self.search_prov = os.path.join(S, "output/search/search_provenance.json")
        write(self.search_prov, {
            "engine": "diann", "version": "2.3.0",
            "engine_version": {"value": "2.3.0", "source": "tools.json versions.diann"},
            "tools": self.tools, "bundle": self.wf, "params_file": self.params,
            "fasta": self.fasta, "files": self.raws, "result": {"report": self.report}})
        self.de = os.path.join(S, "output/tables")
        self.de_prov = os.path.join(self.de, "de_provenance.json")
        write(self.de_prov, {
            "method": "maxlfq", "contrasts": ["B-A"], "q_cutoff": 0.02, "logfc": 0.585,
            "adjp": 0.1, "input": self.report, "metadata": self.conditions,
            "packages": {"limpa": "1.4.0"}})
        write(os.path.join(self.de, "sessionInfo.txt"), "R version 4.5.1\nlimpa 1.4.0\n")
        write(os.path.join(self.de, "DE_B_vs_A.csv"), "Protein,logFC\nP1,1.2\n")
        self.commands = write(os.path.join(S, "logs/commands.log"), "$ run_search.py ...\n")
        self.qc = write(os.path.join(S, "logs/qc_bracket.json"),
                        {"schema": "qc_bracket/1", "schema_version": 1,
                         "checked_at": "2026-10-20T00:00:00Z"})

    def run_prov(self, *args, expect=0):
        p = subprocess.run([sys.executable, os.path.join(SCRIPTS, "provenance.py"),
                            "--setup-json", self.setup_json] + list(args),
                           capture_output=True, text=True, timeout=300, cwd=self.tmp)
        self.assertEqual(p.returncode, expect, p.stdout + p.stderr)
        return p

    def bundle(self, *args):
        p = self.run_prov("--session", self.S, *args)
        res = json.loads(p.stdout)
        with open(os.path.join(res["bundle"], "run_manifest.json"), encoding="utf-8") as fh:
            man = json.load(fh)
        with open(os.path.join(res["bundle"], "MANIFEST.txt"), encoding="utf-8") as fh:
            manifest = fh.read()
        return res, man, manifest, p.stderr


class EveryArgumentFromItsRecord(SessionFixture):
    def test_every_argument_is_filled_from_the_right_file(self):
        res, man, manifest, err = self.bundle()
        fs = res["from_session"]
        self.assertEqual(fs["session"], self.S)
        src = {f: v["from"] for f, v in fs["filled"].items()}
        val = {f: v.get("value") for f, v in fs["filled"].items()}
        self.assertEqual(res["bundle"], os.path.join(self.S, "output", "reproducibility"))
        self.assertEqual(src["--workflow-manifest"], "input/wf/workflow.manifest.json")
        self.assertEqual(src["--params"], "input/wf/params.cfg")
        self.assertEqual(src["--conditions"], "input/conditions.csv")
        self.assertEqual(src["--fasta"], "input/search.fasta")
        self.assertEqual(src["--fasta-info-file"], "input/search.fasta.meta.json")
        self.assertIn("search_provenance.json `files`", src["--raw"])
        self.assertEqual(fs["filled"]["--raw"]["n_files"], 2)
        self.assertEqual(src["--report"], "output/search/report.parquet")
        self.assertEqual(src["--de-dir"], "output/tables")
        self.assertEqual(src["--tools-json"], "output/search/search_provenance.json `tools`")
        self.assertEqual(val["--tools-json"], self.tools)
        self.assertEqual((val["--engine"], val["--engine-version"]), ("diann", "2.3.0"))
        for f in ("--de-method", "--contrasts", "--q-cutoff", "--logfc", "--adjp"):
            self.assertTrue(src[f].startswith("output/tables/de_provenance.json"), (f, src[f]))
        self.assertEqual(src["--acquisition"], "input/wf/workflow.manifest.json `acquisition`")
        self.assertEqual(val["--instrument"], "Orbitrap Astral")
        # the database's taxon (sidecar), not the workflow manifest's
        self.assertEqual((val["--organism-taxid"], src["--organism-taxid"]),
                         ("9606", "input/search.fasta.meta.json `taxid`"))
        self.assertEqual(src["--commands"], "logs/commands.log")
        self.assertEqual(src["--qc-bracket"], "logs/qc_bracket.json")
        self.assertIn("now (UTC)", src["--timestamp"])
        self.assertEqual(fs["not_found"], {})
        self.assertEqual(fs["explicit"], [])

        # ... and the bundle carries them
        self.assertEqual(man["engine"], "diann")
        self.assertEqual(man["engine_version"], "2.3.0")
        self.assertEqual(man["de"], {"method": "maxlfq", "contrasts": "B-A", "q_cutoff": "0.02",
                                     "logfc": "0.585", "adjp": "0.1"})
        self.assertEqual(man["query"], {"acquisition": "DIA", "organism_taxid": "9606",
                                        "instrument": "Orbitrap Astral"})
        self.assertEqual(man["inputs"]["raw"], self.raws)
        self.assertEqual(man["inputs"]["fasta"], self.fasta)
        self.assertEqual(man["inputs"]["fasta_info"]["taxid"], 9606)
        self.assertEqual(man["inputs"]["params"], self.params)
        self.assertEqual(man["instrument_qc"]["record"], self.qc)
        self.assertTrue(man["timestamp"].endswith("Z"))
        for name in ("params.cfg", "conditions.csv", "workflow.manifest.json", "commands.log"):
            self.assertTrue(os.path.isfile(os.path.join(res["bundle"], "inputs", name)), name)
        self.assertIn("[NOTE]    arguments read from session", manifest)
        self.assertIn("arguments read from the session's records", err)

    def test_an_explicit_argument_wins(self):
        out = os.path.join(self.tmp, "elsewhere")
        res, man, _, err = self.bundle("--engine", "sage", "--q-cutoff", "0.05", "--outdir", out,
                                       "--instrument", "timsTOF HT", "--raw", self.raws[0])
        fs = res["from_session"]
        self.assertEqual(res["bundle"], out)
        self.assertEqual(man["engine"], "sage")
        self.assertEqual(man["de"]["q_cutoff"], "0.05")
        self.assertEqual(man["de"]["adjp"], "0.1", "the rest still from the session")
        self.assertEqual(man["query"]["instrument"], "timsTOF HT")
        self.assertEqual(man["inputs"]["raw"], [self.raws[0]])
        for f in ("--engine", "--q-cutoff", "--outdir", "--instrument", "--raw"):
            self.assertIn(f, fs["explicit"])
            self.assertNotIn(f, fs["filled"])
        self.assertIn("given explicitly (these win)", err)

    def test_the_search_s_own_file_wins_over_the_session_s_when_they_differ(self):
        other = write(os.path.join(self.tmp, "db", "other.fasta"), ">sp|P2|Y_MOUSE\nPEPK\n")
        write(other + ".meta.json", {"organism": "Mus musculus", "taxid": 10090,
                                     "content_used": "one_per_gene", "contaminant_set": "none"})
        with open(self.search_prov, encoding="utf-8") as fh:
            sp = json.load(fh)
        sp["fasta"] = other
        write(self.search_prov, sp)
        res, man, manifest, _ = self.bundle()
        fs = res["from_session"]
        self.assertEqual(man["inputs"]["fasta"], other)
        self.assertEqual(fs["filled"]["--fasta-info-file"]["value"], other + ".meta.json")
        self.assertEqual(man["query"]["organism_taxid"], "10090")
        self.assertTrue(any("input/search.fasta is not the file" in n for n in fs["notes"]))
        self.assertIn("[NOTE]    --session: --fasta: the session's input/search.fasta", manifest)


class NothingIsInvented(SessionFixture):
    def test_a_missing_record_is_listed_and_left_unset(self):
        for f in (self.de_prov, self.qc, self.sidecar, self.commands,
                  os.path.join(self.S, "input/acquisition.json")):
            os.remove(f)
        with open(self.search_prov, encoding="utf-8") as fh:
            sp = json.load(fh)
        sp.pop("version"), sp.pop("tools")
        sp["engine_version"] = {"value": None, "source": "unknown -- tools.json says latest"}
        write(self.search_prov, sp)
        with open(self.wf, encoding="utf-8") as fh:
            wf = json.load(fh)
        wf["instruments"], wf["organism_taxid"] = [], None
        write(self.wf, wf)
        res, man, manifest, err = self.bundle()
        nf = res["from_session"]["not_found"]
        for f in ("--de-method", "--contrasts", "--q-cutoff", "--logfc", "--adjp",
                  "--qc-bracket", "--fasta-info-file", "--commands", "--engine-version",
                  "--tools-json", "--instrument", "--organism-taxid"):
            self.assertIn(f, nf)
            self.assertNotIn(f, res["from_session"]["filled"])
        self.assertIn("no output/tables/de_provenance.json", nf["--q-cutoff"])
        self.assertIn("tools.json says latest", nf["--engine-version"])
        # the workflow manifest's DE DEFAULTS (q 0.01, adjp 0.05) are not passed off as the DE's
        self.assertEqual(man["de"], {"method": None, "contrasts": "", "q_cutoff": None,
                                     "logfc": None, "adjp": None})
        self.assertNotIn("engine_version", man)
        self.assertEqual(man["query"]["instrument"], "")
        self.assertIsNone(man["query"]["organism_taxid"])
        self.assertIsNone(man["instrument_qc"], "no QC record: none is referenced")
        self.assertIsNone(man["inputs"]["fasta_info"])
        self.assertIn("NOT FOUND in the session's records", err)
        self.assertIn("[NOTE]    --session: --qc-bracket not found", manifest)

    def test_a_search_only_session_has_no_de(self):
        shutil.rmtree(self.de)
        os.makedirs(self.de)
        res, man, _, _ = self.bundle()
        nf = res["from_session"]["not_found"]
        self.assertIn("no DE has run", nf["--de-dir"])
        self.assertNotIn("--de-dir", res["from_session"]["filled"])
        self.assertIsNone(man["de"]["method"])

    def test_without_session_the_output_is_as_before(self):
        p = self.run_prov("--outdir", os.path.join(self.tmp, "plain"),
                          "--timestamp", "2026-10-01T00:00:00Z")
        self.assertNotIn("from_session", json.loads(p.stdout))
        with open(os.path.join(self.tmp, "plain", "run_manifest.json"), encoding="utf-8") as fh:
            man = json.load(fh)
        self.assertNotIn("engine_version", man)
        self.assertEqual(man["query"]["instrument"], "")


class ABadSessionExits2(SessionFixture):
    def test_no_such_folder(self):
        p = self.run_prov("--session", os.path.join(self.tmp, "nope"), expect=2)
        self.assertIn("no such folder", p.stderr)
        self.assertFalse(os.path.exists(os.path.join(self.tmp, "reproducibility")))

    def test_a_folder_that_is_not_a_session(self):
        os.makedirs(os.path.join(self.tmp, "just_a_folder"))
        p = self.run_prov("--session", os.path.join(self.tmp, "just_a_folder"), expect=2)
        self.assertIn("not an analysis session", p.stderr)

    def test_an_unreadable_session_json(self):
        write(os.path.join(self.S, "session.json"), "{not json")
        p = self.run_prov("--session", self.S, expect=2)
        self.assertIn("session.json cannot be read", p.stderr)
        write(os.path.join(self.S, "session.json"), "[1, 2]")
        p = self.run_prov("--session", self.S, expect=2)
        self.assertIn("is not a JSON object", p.stderr)

    def test_a_readable_session_json_is_fine(self):
        write(os.path.join(self.S, "session.json"),
              {"analysis_for": {"submissions": ["PROT_0807"]}})
        res, _, _, _ = self.bundle()
        self.assertEqual(res["from_session"]["session"], self.S)


if __name__ == "__main__":
    unittest.main()
