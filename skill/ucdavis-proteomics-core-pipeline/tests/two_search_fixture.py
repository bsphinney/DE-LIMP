"""Sessions for the readers of a session's search record (2.14.2, T44). Not a test module.

two_searches(): the real case a 2.14.1 session left the skill's path over -- two timsTOF DDA
runs, searched twice and never into output/search/: a DIA-NN --dda search in
output/search_dda/ whose second pass failed (its first-pass report is what is there), and a Sage
search in output/search_sage/; one DE, which read the Sage report. The workflow manifest pins
DIA-NN (the search the session set out to run), so a reader that falls back to it names the
wrong engine.

one_search(): the same session as 2.14.1 laid it out, its one search in output/search/ (the pin:
a session like this reads as it always has).

layout(): any other arrangement (T44 review, Ruling 78): searches by folder, each with a record
or none, and the report the DE read. FRAGPIPE_DIA is FragPipe's DIA route, whose report is
<workdir>/dia-quant-output/report.tsv beside no record of its own.

Each is laid out by hand, not by `session.py init` (which makes an empty output/search/).
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPTS = os.path.join(os.path.dirname(HERE), "scripts")
if SCRIPTS not in sys.path:
    sys.path.insert(0, SCRIPTS)

import session  # noqa: E402

RUNS = ("runA_01", "runB_01")
DIANN = {"engine": "diann", "version": "2.7.0"}
SAGE = {"engine": "sage", "version": "0.14.7"}
FRAGPIPE = {"engine": "fragpipe", "version": "24.0"}
DE_RECORD = {"pipeline_id": "maxlfq", "display_label": "MaxLFQ + limma", "method": "maxlfq",
             "contrasts": ["B-A"], "adjp": 0.05, "n_samples": 2}
FRAGPIPE_DIA = os.path.join("dia-quant-output", "report.tsv")


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


def _layout(root, name, raw_list):
    sd = os.path.join(root, name)
    p = session.paths_for(sd)
    for d in (p["input_dir"], p["workflow_dir"], p["de_dir"], p["figures_dir"], p["logs_dir"]):
        os.makedirs(d, exist_ok=True)
    raws = [os.path.join(root, "raw", r + ".d") for r in RUNS]
    for r in raws:
        os.makedirs(r, exist_ok=True)
    if raw_list:
        write(p["raw_list"], "# Raw MS files used in this analysis\n" + "\n".join(raws) + "\n")
    write(p["conditions"], "File.Name,Group\n" + "".join(
        f"{r},{g}\n" for r, g in zip(RUNS, ("A", "B"))))
    write(p["workflow_manifest"], json.dumps({
        "acquisition": "DDA", "instruments": ["timsTOF HT"],
        "engine": {"name": "diann", "version": "2.7.0"}}))
    return p, raws


def _search(folder, engine, raws, report, extra=None):
    """A search folder: its report, and its record (`files`, `result.report`) unless engine is
    None (a search with no record: hand-run, or from an older skill)."""
    write(os.path.join(folder, report), "parquet")
    if engine is not None:
        write(os.path.join(folder, "search_provenance.json"), json.dumps(dict(
            engine, files=raws, n_files=len(raws),
            result={"report": os.path.join(folder, report)}, **(extra or {}))))
    return folder


def _de(p, report):
    write(os.path.join(p["de_dir"], "de_provenance.json"), json.dumps(dict(
        DE_RECORD, **({"input": report} if report else {}))))


def two_searches(root, name="dda_two_searches", de=True, raw_list=True):
    """-> (paths_for dict, {"dda": folder, "sage": folder})."""
    p, raws = _layout(root, name, raw_list)
    dda = _search(os.path.join(p["output_dir"], "search_dda"), DIANN, raws,
                  "report-first-pass.parquet", {"search_mode": "dda"})
    write(os.path.join(dda, "report-first-pass.stats.tsv"), "File.Name\tPrecursors.Identified\n")
    write(os.path.join(dda, "params.resolved.cfg"), "--qvalue 0.01\n--dda\n")
    sage = _search(os.path.join(p["output_dir"], "search_sage"), SAGE, raws, "report.parquet")
    for fn in ("lfq.parquet", "results.sage.parquet", "sage.log"):
        write(os.path.join(sage, fn), fn)
    write(os.path.join(sage, "params.resolved.cfg"), '{"database": {"enzyme": {}}}\n')
    if de:
        _de(p, os.path.join(sage, "report.parquet"))
    return p, {"dda": dda, "sage": sage}


def one_search(root, name="dia_one_search", de=True, raw_list=True):
    """-> (paths_for dict, {"search": output/search}): the 2.14.1 layout."""
    p, raws = _layout(root, name, raw_list)
    out = _search(p["search_out"], DIANN, raws, "report.parquet")
    write(os.path.join(out, "params.resolved.cfg"), "--qvalue 0.01\n")
    if de:
        _de(p, os.path.join(out, "report.parquet"))
    return p, {"search": out}


def layout(root, searches, de_read=None, name="layout", raw_list=True):
    """A session with `searches` {folder under output/: (record dict or None, report path in
    it)}; the DE read `de_read` (folder, report) -- or no DE when None, or a DE record without
    `input` when "no-input". Each search gets a params.resolved.cfg naming its folder.
    -> (paths_for dict, {folder: its path})."""
    p, raws = _layout(root, name, raw_list)
    dirs = {}
    for folder, (rec, report) in searches.items():
        dirs[folder] = _search(os.path.join(p["output_dir"], folder), rec, raws, report)
        write(os.path.join(dirs[folder], "params.resolved.cfg"), f"# {folder}\n")
    if de_read == "no-input":
        _de(p, None)
    elif de_read:
        _de(p, os.path.join(dirs[de_read[0]], de_read[1]))
    return p, dirs
