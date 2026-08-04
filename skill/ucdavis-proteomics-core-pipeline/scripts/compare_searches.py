#!/usr/bin/env python3
"""
compare_searches.py  --  Compare two or more SEARCHES of the same raw files, at the
search level: what each engine identified and how well it quantified it.

This is the missing half of the comparison story. compare_analyses.R compares finished
DE tables (which proteins changed); this compares the searches that fed them (which
proteins were seen at all, in how many runs, how consistently). Use it for a tool
bake-off -- e.g. DIA-NN library-free vs FragPipe/diaTracer on the same .d files.

Input is the DE contract every engine in this skill already produces: a DIA-NN-shaped
report.parquet with Run, Protein.Group, PG.MaxLFQ (+ q-value columns when available).

Reports, per search and per pair:
  - protein groups and runs identified
  - per-run protein counts (depth, and how even it is across runs)
  - data completeness: how many proteins are quantified in every run
  - quantitative precision: CV of each protein across runs (replicate-level only --
    with no group labels a low CV means reproducible measurement, not no biology)
  - overlap between searches (shared / unique protein groups, Jaccard)
  - agreement on shared proteins: Pearson/Spearman of log2 intensity, per run

Usage:
  python3 compare_searches.py --out ./search_comparison \\
    --search "DIA-NN:/path/a/report.parquet" \\
    --search "diaTracer:/path/b/report.parquet" [--q 0.01]
"""
import argparse, json, math, os, statistics, sys
from collections import defaultdict


def read_report(path, q_cut):
    """Read a DIA-NN-shaped report into {run: {protein: intensity}}.

    Filters on the GLOBAL q-value when the column exists. Run-level Q.Value is not
    comparable across tools -- the global one is the defensible cross-tool cutoff --
    so prefer Global.Q.Value / Lib.PG.Q.Value and fall back only when absent.
    """
    if not os.path.exists(path):
        sys.exit(f"No such report: {path}")

    # The FragPipe/diaTracer route emits report.TSV, not parquet — FragPipe 24 bundles
    # DIA-NN 1.8.2 beta 8, which has no parquet writer. Assuming parquet made this tool
    # unable to read the very output it exists to compare against. Accept either.
    is_tsv = path.lower().endswith((".tsv", ".txt", ".csv"))
    if is_tsv:
        import csv as _csv
        delim = "," if path.lower().endswith(".csv") else "\t"
        with open(path, newline="") as fh:
            rdr = _csv.DictReader(fh, delimiter=delim)
            names = rdr.fieldnames or []
            rows = list(rdr)
        cols = {c.lower(): c for c in names}
    else:
        try:
            import pyarrow.parquet as pq
        except ImportError:
            sys.exit("pyarrow is required for parquet reports. It ships with the skill's conda env (setup.sh).")
        tbl = pq.read_table(path)
        cols = {c.lower(): c for c in tbl.column_names}

    def col(*names):
        for n in names:
            if n.lower() in cols:
                return cols[n.lower()]
        return None

    c_run, c_pg = col("Run"), col("Protein.Group")
    c_int = col("PG.MaxLFQ", "PG.Quantity", "Intensity")
    if not all([c_run, c_pg, c_int]):
        sys.exit(f"{path} is not a DE-contract report (need Run, Protein.Group, PG.MaxLFQ).")
    c_q = col("Global.PG.Q.Value", "Global.Q.Value", "Lib.PG.Q.Value", "Q.Value")
    q_basis = c_q or "(none — no q-value column, nothing filtered)"

    if is_tsv:
        runs = [r.get(c_run) for r in rows]
        pgs = [r.get(c_pg) for r in rows]
        ints = [r.get(c_int) for r in rows]
        qs = [r.get(c_q) for r in rows] if c_q else [0.0] * len(rows)
    else:
        runs = tbl.column(c_run).to_pylist()
        pgs = tbl.column(c_pg).to_pylist()
        ints = tbl.column(c_int).to_pylist()
        qs = tbl.column(c_q).to_pylist() if c_q else [0.0] * len(runs)

    data, dropped = defaultdict(dict), 0
    for run, pg, val, q in zip(runs, pgs, ints, qs):
        try:
            q = float(q)
        except (TypeError, ValueError):
            q = 0.0
        if q > q_cut:
            dropped += 1
            continue
        try:
            v = float(val)
        except (TypeError, ValueError):
            continue
        if v > 0 and pg:
            data[str(run)][str(pg)] = v
    return data, q_basis, dropped


def cv_percent(values):
    if len(values) < 2:
        return None
    m = statistics.mean(values)
    if m <= 0:
        return None
    return 100.0 * statistics.stdev(values) / m


def pearson(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    return num / (dx * dy) if dx and dy else None


def spearman(xs, ys):
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        r = [0.0] * len(v)
        i = 0
        while i < len(order):                      # average ties
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    return pearson(rank(xs), rank(ys))


def summarize(label, data):
    runs = sorted(data)
    all_pgs = set()
    for r in runs:
        all_pgs |= set(data[r])
    per_run = {r: len(data[r]) for r in runs}
    complete = sum(1 for p in all_pgs if all(p in data[r] for r in runs))
    cvs = []
    for p in all_pgs:
        vals = [data[r][p] for r in runs if p in data[r]]
        c = cv_percent(vals)
        if c is not None:
            cvs.append(c)
    counts = list(per_run.values())
    return {
        "label": label,
        "runs": len(runs),
        "protein_groups": len(all_pgs),
        "per_run_min": min(counts) if counts else 0,
        "per_run_median": int(statistics.median(counts)) if counts else 0,
        "per_run_max": max(counts) if counts else 0,
        "complete_in_all_runs": complete,
        "completeness_pct": round(100.0 * complete / len(all_pgs), 1) if all_pgs else 0.0,
        "median_cv_pct": round(statistics.median(cvs), 1) if cvs else None,
        "per_run_counts": per_run,
        "_pgs": all_pgs,
    }


def compare_pair(a, b, data_a, data_b):
    pa, pb = a["_pgs"], b["_pgs"]
    shared = pa & pb
    union = pa | pb
    # Agreement on shared proteins, per shared run, on log2 intensity.
    per_run = {}
    for run in sorted(set(data_a) & set(data_b)):
        xs, ys = [], []
        for p in shared:
            va, vb = data_a[run].get(p), data_b[run].get(p)
            if va and vb:
                xs.append(math.log2(va))
                ys.append(math.log2(vb))
        if len(xs) >= 3:
            per_run[run] = {"n": len(xs),
                            "pearson": round(pearson(xs, ys), 4),
                            "spearman": round(spearman(xs, ys), 4)}
    rs = [v["pearson"] for v in per_run.values()]
    return {
        "pair": f"{a['label']} vs {b['label']}",
        "shared": len(shared),
        f"only_{a['label']}": len(pa - pb),
        f"only_{b['label']}": len(pb - pa),
        "jaccard": round(len(shared) / len(union), 4) if union else 0.0,
        "runs_compared": len(per_run),
        "median_pearson_log2": round(statistics.median(rs), 4) if rs else None,
        "per_run": per_run,
        "_shared": shared, "_only_a": pa - pb, "_only_b": pb - pa,
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--search", action="append", default=[], metavar="LABEL:PATH",
                    help='repeat: --search "DIA-NN:/path/report.parquet"')
    ap.add_argument("--out", default="search_comparison")
    ap.add_argument("--q", type=float, default=0.01, help="q-value cutoff (default 0.01)")
    a = ap.parse_args()
    if len(a.search) < 2:
        sys.exit('Need >= 2 --search "Label:/path/report.parquet" arguments.')
    os.makedirs(a.out, exist_ok=True)

    loaded, summaries = [], []
    for spec in a.search:
        if ":" not in spec:
            sys.exit(f'--search must be "Label:/path", got: {spec}')
        label, path = spec.split(":", 1)
        data, q_basis, dropped = read_report(path.strip(), a.q)
        if not data:
            sys.exit(f"{label}: no rows survived the q <= {a.q} filter ({path}).")
        s = summarize(label.strip(), data)
        s["q_basis"], s["rows_dropped_by_q"], s["path"] = q_basis, dropped, path.strip()
        loaded.append(data)
        summaries.append(s)

    pairs = []
    for i in range(len(summaries)):
        for j in range(i + 1, len(summaries)):
            pairs.append(compare_pair(summaries[i], summaries[j], loaded[i], loaded[j]))

    # ---- CSVs -----------------------------------------------------------------
    import csv
    with open(os.path.join(a.out, "search_summary.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["search", "runs", "protein_groups", "per_run_min", "per_run_median",
                    "per_run_max", "complete_in_all_runs", "completeness_pct",
                    "median_cv_pct", "q_basis", "rows_dropped_by_q"])
        for s in summaries:
            w.writerow([s["label"], s["runs"], s["protein_groups"], s["per_run_min"],
                        s["per_run_median"], s["per_run_max"], s["complete_in_all_runs"],
                        s["completeness_pct"], s["median_cv_pct"], s["q_basis"],
                        s["rows_dropped_by_q"]])
    for p in pairs:
        tag = p["pair"].replace(" ", "_").replace("/", "-")
        with open(os.path.join(a.out, f"unique_{tag}.csv"), "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["protein_group", "found_only_in"])
            la, lb = p["pair"].split(" vs ")
            for pg in sorted(p["_only_a"]):
                w.writerow([pg, la])
            for pg in sorted(p["_only_b"]):
                w.writerow([pg, lb])

    # ---- markdown -------------------------------------------------------------
    md = ["# Search comparison", "",
          f"Same raw files, {len(summaries)} searches, q &le; {a.q}.", "",
          "## What each search found", "",
          "| Search | Runs | Protein groups | Per-run (min/median/max) | In every run | Median CV |",
          "|---|---|---|---|---|---|"]
    for s in summaries:
        cv = f"{s['median_cv_pct']}%" if s["median_cv_pct"] is not None else "n/a"
        md.append(f"| {s['label']} | {s['runs']} | {s['protein_groups']} | "
                  f"{s['per_run_min']} / {s['per_run_median']} / {s['per_run_max']} | "
                  f"{s['complete_in_all_runs']} ({s['completeness_pct']}%) | {cv} |")
    md += ["",
           "Depth alone does not decide which search is better: a search can identify more "
           "proteins while quantifying them less reproducibly, or find them in fewer runs. "
           "Read the three columns together — protein groups, how many survive in *every* "
           "run, and the CV.", ""]
    for s in summaries:
        md.append(f"- **{s['label']}** — q filtered on `{s['q_basis']}`, "
                  f"{s['rows_dropped_by_q']} rows dropped. `{s['path']}`")
    md += ["", "## How the searches agree", ""]
    for p in pairs:
        la, lb = p["pair"].split(" vs ")
        md += [f"### {p['pair']}", "",
               f"- Shared protein groups: **{p['shared']}** (Jaccard {p['jaccard']})",
               f"- Only {la}: {p[f'only_{la}']}  ·  Only {lb}: {p[f'only_{lb}']}"]
        if p["median_pearson_log2"] is not None:
            md.append(f"- Agreement on shared proteins, log2 intensity: median Pearson "
                      f"**{p['median_pearson_log2']}** across {p['runs_compared']} runs")
        else:
            md.append("- Not enough shared proteins in shared runs to correlate intensities.")
        md += ["", "Proteins found by only one search are the interesting ones — they are "
               "where the two approaches genuinely differ. They are listed in "
               f"`unique_{p['pair'].replace(' ', '_')}.csv`.", ""]
    md += ["## Caveats", "",
           "- CV here is computed across **all runs**, with no knowledge of experimental "
           "groups. If the runs span real conditions, a high CV may be biology rather than "
           "imprecision. It is only a precision measure when the runs are replicates.",
           "- Protein-group identifiers must be comparable between searches. Different "
           "FASTA files, or different protein-inference rules, shift group membership and "
           "will depress the overlap for reasons that have nothing to do with sensitivity.",
           "- This compares searches, not conclusions. Run `compare_analyses.R` on the DE "
           "output to see whether the differences here actually change which proteins come "
           "out as significant.", ""]
    with open(os.path.join(a.out, "SEARCH_COMPARISON.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")

    for s in summaries:
        s.pop("_pgs", None)
    for p in pairs:
        for k in ("_shared", "_only_a", "_only_b"):
            p.pop(k, None)
    result = {"out": a.out, "q_cutoff": a.q, "searches": summaries, "pairs": pairs,
              "report": os.path.join(a.out, "SEARCH_COMPARISON.md")}
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
