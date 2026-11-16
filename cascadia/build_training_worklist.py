"""Worklist for prepare_training_native.py from FRAN: one DIA-NN search per diaPASEF .d run, with
the corpus split (delimp_run_flags, group-aware: a group never straddles splits).

Picks, per run, the DIA-NN search reporting the most precursors for it. Keeps a run only if its
report.parquet and .d exist and the .d has no stale analysis.tdf-wal (memory:
bruker-tdf-wal-truncation). Diversity: at most --per-group runs from one group_key per split, so a
200-run study can't dominate. Read-only.

    build_training_worklist.py OUT.tsv [--per-group 6] [--train 300] [--val 40] [--test 40] [--seed 0]
"""
import argparse, csv, os, random, sys
sys.path.insert(0, "/quobyte/proteomics-grp/brett/glendon/fran_ingest")
from coreomics_import import _conn

ap = argparse.ArgumentParser()
ap.add_argument("out"); ap.add_argument("--per-group", type=int, default=6)
ap.add_argument("--train", type=int, default=300); ap.add_argument("--val", type=int, default=40)
ap.add_argument("--test", type=int, default=40); ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--always-test", default="08302023_mouse_400ng_dia_100spd_beads_2h_hl_s2-g3_1_718",
                help="run_key(s) that must never be trained on (comma-separated): the benchmark run")
a = ap.parse_args()
conn = _conn(); conn.set_session(readonly=True, autocommit=True); c = conn.cursor()
c.execute("SET statement_timeout = '10min'")
c.execute("""
SELECT DISTINCT ON (f.run_key) f.run_key, f.split, coalesce(f.group_key, f.run_key), s.id::text, s.output_dir,
       rf.hive_path, srf.n_precursors
  FROM delimp_searches s
  JOIN search_raw_files srf ON srf.search_id = s.id
  JOIN raw_files rf ON rf.raw_path = srf.raw_path
  JOIN delimp_run_flags f ON f.run_key = lower(rf.raw_basename)
 WHERE s.search_engine = 'diann' AND rf.hive_path ILIKE '%%.d'
   AND coalesce(rf.acquisition_method, '') ~* 'dia' AND s.output_dir IS NOT NULL
   AND coalesce(srf.n_precursors, 0) >= 2000
 ORDER BY f.run_key, srf.n_precursors DESC NULLS LAST""")
rows = c.fetchall()
print("candidate runs:", len(rows))
held = set(x.strip().lower() for x in a.always_test.split(",") if x.strip())
ok, why = [], {}
for rk, split, grp, sid, od, hp, npre in rows:
    rep = os.path.join(od, "report.parquet")
    reason = None
    if not os.path.isfile(rep): reason = "no_report"
    elif not os.path.isdir(hp): reason = "no_d"
    elif not os.path.isfile(os.path.join(hp, "analysis.tdf")): reason = "no_tdf"
    elif os.path.exists(os.path.join(hp, "analysis.tdf-wal")) and os.path.getsize(os.path.join(hp, "analysis.tdf-wal")) > 0:
        reason = "stale_wal"
    if reason:
        why[reason] = why.get(reason, 0) + 1; continue
    if rk in held: split = "test"
    ok.append((rk, split, grp, sid, rep, hp, npre))
print("usable:", len(ok), "dropped:", why)
rng = random.Random(a.seed)
rng.shuffle(ok)
take = {"train": a.train, "val": a.val, "test": a.test}
per = {}; out = []
# benchmark run(s) first so they are always in test
ok.sort(key=lambda r: r[0] not in held)
for r in ok:
    sp, grp = r[1], r[2]
    if sp not in take: continue
    n_sp = sum(1 for x in out if x[1] == sp)
    if n_sp >= take[sp] and r[0] not in held: continue
    if per.get((sp, grp), 0) >= a.per_group and r[0] not in held: continue
    per[(sp, grp)] = per.get((sp, grp), 0) + 1
    out.append(r)
with open(a.out, "w", newline="") as f:
    w = csv.writer(f, delimiter="\t"); w.writerow(["run_key", "split", "group_key", "search_id", "report", "d", "n_precursors"])
    w.writerows(out)
from collections import Counter
print("worklist:", dict(Counter(r[1] for r in out)), "groups:", {s: len({r[2] for r in out if r[1]==s}) for s in take})
print("benchmark in test:", [r[0] for r in out if r[0] in held])
