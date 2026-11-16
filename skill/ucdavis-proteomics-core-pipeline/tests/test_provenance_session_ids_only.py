#!/usr/bin/env python3
"""provenance.py --session builds the bundle finalize builds (2.15 review H1).

SKILL.md step 10 runs `provenance.py --session <session>`, and finalize keeps the bundle it
finds (report_files.bundle_state "present"). In the 2.15 slim port --session worked its flags
out on its own and never passed the search record or the search's other results, so an
identifications-only search (Sage --ids-only, a FragPipe nonspecific search: no report.parquet)
was bundled with `search_record: null`, no `search_outputs`, and a REPRODUCE.md saying "version
not recorded" and "the search's report is the result" -- a report that does not exist. finalize's
own command (make_deposit.bundle_command) had them. Now both take their flags from
make_deposit.bundle_flags, the one definition. What these tests pin:
  * an identifications-only Sage session and a FragPipe nonspecific one: run_manifest.json names
    the search record, checksums.json fingerprints identifications.tsv in `search_outputs`, and
    REPRODUCE.md says "Identifications only" with the engine's version;
  * for those, a quantified DIA-NN session without DE and a DIA-NN DE session, the bundle
    --session builds is the one finalize's command builds, file for file (bar the time it was
    made and the [NOTE] lines that say where --session read each value);
  * bundle_command now passes --engine-version (the one value only --session used to fill);
  * one replay rule (Ruling 96): the FASTA the search ran on, as its params, with its own sidecar
    and taxon, and a "What ran" section in REPRODUCE.md when the session's own file differs.
provenance.py runs for real against a fake skill home (no conda, no R). Stdlib only.
"""
import json
import os
import re
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
sys.path.insert(0, SCRIPTS)
sys.path.insert(0, HERE)

from job_env import job_env                                      # noqa: E402
from test_reproduce_sh_runs import arg                           # noqa: E402
import test_reproduce_replays_what_ran as rw                     # noqa: E402

PY = sys.executable
read = rw.read
write = rw.write


class Sessions(rw.Base):
    sage_ids_only = rw.FinalizesCommand.session        # (p, cfg, ids): Sage --ids-only, no DE

    def fragpipe_nonspecific(self):
        """A FragPipe search with the Nonspecific-peptidome workflow: identifications only, as
        run_search.py --adapt-only leaves it (identifications.tsv + fragpipe_identifications.json,
        no report.parquet) -> (p, workflow, ids)."""
        import test_deposit_package as tdp
        root = os.path.join(self.d, "s")
        os.makedirs(root)
        p = tdp.dia_session(root, de=False)
        for n in ("report.parquet", "report.log.txt", "report.pg_matrix.tsv"):
            os.remove(os.path.join(p["search_out"], n))
        os.remove(os.path.join(p["workflow_dir"], "params.cfg"))
        given = write(os.path.join(p["workflow_dir"], "fragpipe.workflow"),
                      "msfragger.search_enzyme_name_1=nonspecific\n"
                      "msfragger.num_enzyme_termini=0\nquantitation.run-label-free-quant=false\n")
        search = p["search_out"]
        ids = write(os.path.join(search, "identifications.tsv"), "run\tpeptide\nA\tPEPTIDE\n")
        quant = {"status": "not_applicable", "engine": "fragpipe", "quantified": False,
                 "why": "the workflow ran no label-free quantification"}
        write(os.path.join(search, "fragpipe_identifications.json"), json.dumps(
            {"kept_rule": "Philosopher filter --sequential --prot 1", "quantification": quant}))
        write(p["search_prov"], json.dumps({
            "engine": "fragpipe", "version": "24.0", "params_file": given, "quantified": False,
            "quantification": quant, "result": {"engine": "fragpipe", "quantified": False}}))
        write(p["workflow_manifest"], json.dumps({
            "acquisition": "DDA", "instruments": ["timsTOF HT"],
            "engine": {"name": "fragpipe", "version": "24.0"}}))
        return p, given, ids

    def diann(self, de):
        import test_deposit_package as tdp
        root = os.path.join(self.d, "s")
        os.makedirs(root)
        return tdp.dia_session(root, de=de)

    def diann_ran_on_another_fasta(self):
        """A DIA-NN DE session whose search record names a FASTA other than input/search.fasta
        (a mouse database with its own sidecar) -> (p, that FASTA, its sidecar)."""
        p = self.diann(de=True)
        other = write(os.path.join(self.d, "db", "mouse.fasta"), ">sp|P2|Y_MOUSE\nPEPK\n")
        side = write(other + ".meta.json", json.dumps({
            "fasta": other, "organism": "Mus musculus", "taxid": 10090,
            "content_used": "one_per_gene", "contaminant_set": "none"}))
        rec = json.loads(read(p["search_prov"]))
        rec["fasta"] = other
        write(p["search_prov"], json.dumps(rec))
        return p, other, side

    def env(self):
        return job_env(self.d, PROTEOMICS_PIPELINE_HOME=self.home)

    def run_cmd(self, argv):
        r = subprocess.run(argv, capture_output=True, text=True, timeout=300, env=self.env())
        self.assertEqual(r.returncode, 0, r.stdout[-2000:] + r.stderr[-2000:])
        return r

    def by_session(self, p, outdir=None):
        """provenance.py --session (step 10) -> the bundle folder."""
        more = ["--outdir", outdir] if outdir else []
        r = self.run_cmd([PY, os.path.join(SCRIPTS, "provenance.py"), "--session",
                      p["session_dir"], *more])
        return json.loads(r.stdout)["bundle"]

    def by_finalize(self, p, outdir):
        """finalize's command (make_deposit.bundle_command), as ensure_bundle runs it."""
        import make_deposit
        from unittest import mock
        with mock.patch.dict(os.environ, self.env(), clear=True):
            cmd, missing, _ = make_deposit.bundle_command(make_deposit.gather(p["session_dir"]),
                                                          outdir=outdir)
        self.assertEqual(missing, [])
        self.run_cmd(cmd)
        return outdir, cmd


def bundle_files(d):
    """{relative path: text} of a bundle, with what legitimately differs between two builds of
    the same session taken out: the folder's own path, the time it was made, and the [NOTE]
    lines --session adds (where each value was read from)."""
    out = {}
    for dp, _, fns in os.walk(d):
        for fn in fns:
            path = os.path.join(dp, fn)
            rel = os.path.relpath(path, d)
            with open(path, encoding="utf-8", errors="replace") as fh:
                text = fh.read().replace(d, "<bundle>")
            if rel == "run_manifest.json":
                man = json.loads(text)
                man.pop("timestamp", None)
                text = json.dumps(man, indent=1, sort_keys=True)
            if rel == "MANIFEST.txt":
                text = "\n".join(ln for ln in text.splitlines() if not ln.startswith("[NOTE]"))
            text = re.sub(r"\d{4}-\d\d-\d\dT\d\d:\d\d:\d\d(\.\d+)?Z?", "<time>", text)
            out[rel] = text
    return out


class IdentificationsOnlyBundles(Sessions):
    def check(self, p, ids_path, engine, version):
        bundle = self.by_session(p)
        man = json.loads(read(os.path.join(bundle, "run_manifest.json")))
        self.assertEqual(man["search_record"], p["search_prov"])
        self.assertEqual(man.get("engine_version"), version)
        checks = json.loads(read(os.path.join(bundle, "checksums", "checksums.json")))
        self.assertIn(ids_path, [x["path"] for x in checks.get("search_outputs") or []])
        self.assertEqual(checks["search_record"]["path"], p["search_prov"])
        md = read(os.path.join(bundle, "REPRODUCE.md"))
        self.assertIn("Identifications only: the search did not quantify", md)
        self.assertIn(f"{version} (its record: search_provenance.json)", md)
        self.assertNotIn("the search's report is the result", md)
        self.assertNotIn("version not recorded", md)
        return bundle

    def test_sage_ids_only(self):
        p, _cfg, ids = self.sage_ids_only()
        self.check(p, ids["identifications"], "sage", "0.14.7")

    def test_fragpipe_nonspecific(self):
        p, _given, ids = self.fragpipe_nonspecific()
        self.check(p, ids, "fragpipe", "24.0")


class TheSameBundleAsFinalize(Sessions):
    def same(self, p):
        a = self.by_session(p, outdir=os.path.join(self.d, "by_session"))
        b, cmd = self.by_finalize(p, os.path.join(self.d, "by_finalize"))
        fa, fb = bundle_files(a), bundle_files(b)
        self.assertEqual(sorted(fa), sorted(fb))
        for rel in fa:
            self.assertEqual(fa[rel], fb[rel], rel)
        return cmd

    def test_sage_ids_only(self):
        self.same(self.sage_ids_only()[0])

    def test_fragpipe_nonspecific(self):
        self.same(self.fragpipe_nonspecific()[0])

    def test_diann_without_de(self):
        self.same(self.diann(de=False))

    def test_diann_with_de(self):
        cmd = self.same(self.diann(de=True))
        self.assertEqual(arg(cmd, "--engine-version"), "2.7.0")
        self.assertEqual(arg(cmd, "--de-method"), "dpc")
        self.assertNotIn("--note", cmd, "nothing differs: REPRODUCE.md as before")

    def test_a_search_that_ran_on_another_fasta(self):
        self.same(self.diann_ran_on_another_fasta()[0])


class OneReplayRule(Sessions):
    """Ruling 96: the bundle replays what RAN. --params and the FASTA both take the search
    record's file first; REPRODUCE.md says so when the session's own file is another one."""

    def test_finalize_replays_the_fasta_the_search_ran_on(self):
        import make_deposit
        p, other, side = self.diann_ran_on_another_fasta()
        out = os.path.join(self.d, "by_finalize")
        _b, cmd = self.by_finalize(p, out)
        self.assertEqual(arg(cmd, "--fasta"), other)
        self.assertEqual(arg(cmd, "--fasta-info-file"), side)
        self.assertEqual(arg(cmd, "--organism-taxid"), "10090")
        notes = [cmd[i + 1] for i, w in enumerate(cmd) if w == "--note"]
        self.assertEqual(len(notes), 1, notes)
        self.assertIn(f"the search ran on {other}", notes[0])
        man = json.loads(read(os.path.join(out, "run_manifest.json")))
        self.assertEqual(man["inputs"]["fasta"], other)
        self.assertEqual(man["inputs"]["fasta_info"]["taxid"], 10090)
        md = read(os.path.join(out, "REPRODUCE.md"))
        self.assertIn("## What ran, where the session's own files differ", md)
        self.assertIn("the session's input/search.fasta is another file and is not used", md)

    def test_the_params_the_search_was_given_are_said_too(self):
        p, cfg, _ids = self.sage_ids_only()           # an older params.json beside the cfg
        bundle = self.by_session(p)
        md = read(os.path.join(bundle, "REPRODUCE.md"))
        section = md.split("## What ran, where the session's own files differ", 1)[1]
        self.assertIn(f"Parameters: the search was given {cfg}", section.split("##")[0])
        self.assertIn("input/wf/params.json is another file and is not used", section)


if __name__ == "__main__":
    unittest.main()
