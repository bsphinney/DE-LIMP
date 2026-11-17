"""Would splitting each diaPASEF window by ion mobility put a peptide's fragments into the top-N peaks?

For a sample of DIA-NN q<=0.01 precursors in one run, take the DIA frame + window nearest the apex
and measure b/y coverage (1+/2+, an ion counted once, 15 ppm, calibrated m/z) in the top-N peaks of:
  whole   : the window centroided over all its mobility scans (what training/inference use now)
  oracle h: scans within +/-h of the fragment-signal apex (upper bound; needs the apex)
  grid W  : fixed slices of W scans at stride W/2; the slice whose centre is nearest the apex
            (what an IM-sliced inference would see, with the label put on the right slice)
Also fits DIA-NN IM (1/K0) against the apex scan, to see if a label's slice can come from DIA-NN.

    python im_cov.py report.parquet RUN RUN.d RUN.calib.json [n=800]
"""
import json
import os
import sqlite3
import sys
from collections import OrderedDict

import numpy as np
import pyarrow.compute as pc
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import prepare_training_native as P  # noqa: E402
from bruker_native_augment import _tdf_meta, _without_appledouble, tof_to_mz  # noqa: E402

rep, run, d, calp = sys.argv[1:5]
N = int(sys.argv[5]) if len(sys.argv) > 5 else 800
PPM = 15.0
cal2 = json.load(open(calp))["ms2"]

cols = ["Run", "Modified.Sequence", "Precursor.Charge", "Precursor.Mz", "RT", "Q.Value"]
names = pq.read_schema(rep).names
if "IM" in names:
    cols.append("IM")
t = pq.read_table(rep, columns=cols, filters=[("Run", "==", run)])
t = t.filter(pc.less_equal(t["Q.Value"], 0.01))
df = t.to_pandas().drop_duplicates(["Modified.Sequence", "Precursor.Charge"])
df = df.iloc[np.random.default_rng(0).permutation(len(df))]

labels = []
for r in df.itertuples(index=False):
    lab, masses = P.to_proforma(r[1])
    if lab is None or not 2 <= int(r[2]) <= 4:
        continue
    nadd = 42.01057 if lab.startswith("[Acetyl]-") else (43.00581 if lab.startswith("[Carbamyl]-") else 0.0)
    labels.append((float(r[3]), float(r[4]) * 60.0, masses, nadd, float(r[6]) if "IM" in names else np.nan))
    if len(labels) >= N:
        break

with _without_appledouble(d) as dr:
    tdf = os.path.join(dr, "analysis.tdf")
    con = sqlite3.connect(f"file:{tdf}?mode=ro&immutable=1", uri=True)
    fr_rows = con.execute("SELECT f.Id, f.Time, i.WindowGroup FROM Frames f JOIN DiaFrameMsMsInfo i ON i.Frame=f.Id").fetchall()
    win_rows = con.execute("SELECT WindowGroup, ScanNumBegin, ScanNumEnd, IsolationMz, IsolationWidth FROM DiaFrameMsMsWindows").fetchall()
    con.close()
    mz_lo, mz_hi, dns = _tdf_meta(dr)
    import timsrust_pyo3 as ts
    reader = ts.FrameReader(dr)

    fid = np.array([r[0] for r in fr_rows]); ftime = np.array([r[1] for r in fr_rows]); fgrp = np.array([r[2] for r in fr_rows])
    wins = {}
    for g, b, e, mz, w in win_rows:
        wins.setdefault(g, []).append((int(b), int(e), float(mz), float(w)))
    cyc = float(np.median(np.diff(np.unique(ftime))))

    # which read_frame index is tdf Id?  (0- or 1-based)
    probe = reader.read_frame(int(fid[len(fid) // 2]))
    off = 0 if abs(probe.rt - ftime[len(fid) // 2]) < 1e-3 else -1
    probe = reader.read_frame(int(fid[len(fid) // 2]) + off)
    assert abs(probe.rt - ftime[len(fid) // 2]) < 1e-3, (probe.rt, ftime[len(fid) // 2])

    def assign(mz, rt):
        best = None
        sel = np.nonzero(np.abs(ftime - rt) < 2 * cyc)[0]
        for k in sel:
            for (b, e, c, w) in wins[fgrp[k]]:
                if c - w / 2 < mz < c + w / 2:
                    dd = abs(ftime[k] - rt)
                    if best is None or dd < best[0]:
                        best = (dd, int(fid[k]), b, e)
        return best

    grp_frames = {}
    for k in np.argsort(ftime):
        grp_frames.setdefault(fgrp[k], []).append(int(fid[k]))
    pos = {f: (g, i) for g, fl in grp_frames.items() for i, f in enumerate(fl)}

    def neighbours(f, k):
        g, i = pos[f]
        fl = grp_frames[g]
        return [fl[j] for j in range(max(0, i - k), min(len(fl), i + k + 1))]

    jobs = []
    for i, p in enumerate(labels):
        a = assign(p[0], p[1])
        if a is not None:
            jobs.append((a[1], a[2], a[3], i))
    jobs.sort()

    cache = OrderedDict()

    def frame(f):
        if f not in cache:
            fr = reader.read_frame(f + off)
            cache[f] = (np.asarray(fr.scan_offsets, dtype=np.int64), np.asarray(fr.tof_indices, dtype=np.int64),
                        np.asarray(fr.intensities, dtype=np.float64))
            if len(cache) > 64:
                cache.popitem(last=False)
        return cache[f]

    def centroid(tof, w):
        if tof.size == 0:
            return np.empty(0), np.empty(0)
        lo = int(tof.min()) - 2
        c = np.concatenate([np.bincount(tof - lo, weights=w), [0.0, 0.0]])
        left, right = np.roll(c, 1), np.roll(c, -1); left[0] = right[-1] = 0.0
        s = left + c + right
        sl, sr = np.roll(s, 1), np.roll(s, -1); sl[0] = sr[-1] = 0.0
        idx = np.nonzero((s > 0) & (s >= sl) & (s > sr) & (c > 0))[0]
        num = (idx - 1) * left[idx] + idx * c[idx] + (idx + 1) * right[idx]
        mz = tof_to_mz(num / s[idx] + lo, mz_lo, mz_hi, dns)
        return P.correct(mz, cal2), s[idx]

    def cov_top(tof, w, masses, nadd, n):
        mz, it = centroid(tof, w)
        mz, _ = P._top(mz, it, n)
        return P.by_coverage(masses, mz, nadd, PPM), mz.size

    H = (8, 16, 32)
    W = (32, 64)
    res = {k: [] for k in ["whole150", "whole300", "whole500", "whole1000", "npk_whole", "sum3_150", "sum5_150",
                           "sum3_300", "sum3_or8_150"] + [f"or{h}" for h in H] + [f"grid{w}" for w in W]
           + [f"npk_grid{w}" for w in W]}
    imfit = []
    nskip = 0
    for f, b, e, i in jobs:
        mzp, rt, masses, nadd, im = labels[i]
        so, tof, inten = frame(f)
        lo_, hi_ = so[b], so[min(e, len(so) - 1)]
        tof_w, int_w = tof[lo_:hi_], inten[lo_:hi_]
        scan_w = np.repeat(np.arange(b, min(e, len(so) - 1)), np.diff(so[b:min(e, len(so) - 1) + 1]))
        # fragment-signal profile over scans
        theo = np.sort(P.by_ions(masses, nadd).ravel())
        mzr = P.correct(tof_to_mz(tof_w, mz_lo, mz_hi, dns), cal2)
        j = np.clip(np.searchsorted(theo, mzr), 1, len(theo) - 1)
        near = np.minimum(np.abs(theo[j] - mzr), np.abs(theo[j - 1] - mzr))
        hit = near <= mzr * 20e-6
        if hit.sum() < 3:
            nskip += 1
            continue
        prof = np.bincount(scan_w[hit] - b, weights=int_w[hit], minlength=e - b)
        prof = np.convolve(prof, np.ones(11), mode="same")
        apex = b + int(np.argmax(prof))
        imfit.append((im, apex))
        c, n = cov_top(tof_w, int_w, masses, nadd, 150); res["whole150"].append(c); res["npk_whole"].append(n)
        for n in (300, 500, 1000):
            res[f"whole{n}"].append(cov_top(tof_w, int_w, masses, nadd, n)[0])
        for k in (1, 2):
            tt, ww, ss = [], [], []
            for g in neighbours(f, k):
                so2, tof2, int2 = frame(g)
                e2 = min(e, len(so2) - 1)
                tt.append(tof2[so2[b]:so2[e2]]); ww.append(int2[so2[b]:so2[e2]])
                ss.append(np.repeat(np.arange(b, e2), np.diff(so2[b:e2 + 1])))
            tt, ww, ss = np.concatenate(tt), np.concatenate(ww), np.concatenate(ss)
            res[f"sum{2 * k + 1}_150"].append(cov_top(tt, ww, masses, nadd, 150)[0])
            if k == 1:
                res["sum3_300"].append(cov_top(tt, ww, masses, nadd, 300)[0])
                m = np.abs(ss - apex) <= 8
                res["sum3_or8_150"].append(cov_top(tt[m], ww[m], masses, nadd, 150)[0])
        for h in H:
            m = np.abs(scan_w - apex) <= h
            res[f"or{h}"].append(cov_top(tof_w[m], int_w[m], masses, nadd, 150)[0])
        for wd in W:
            starts = np.arange(b, max(b + 1, e - wd + 1), wd // 2)
            s0 = starts[np.argmin(np.abs(starts + wd / 2 - apex))]
            m = (scan_w >= s0) & (scan_w < s0 + wd)
            c, n = cov_top(tof_w[m], int_w[m], masses, nadd, 150)
            res[f"grid{wd}"].append(c); res[f"npk_grid{wd}"].append(n)

print(f"{run}: {len(labels)} labels, {len(jobs)} assigned, {nskip} with <3 fragment hits skipped, "
      f"{len(res['whole150'])} measured; window scan span median "
      f"{np.median([e - b for v in wins.values() for b, e, _, _ in v]):.0f}; cycle {cyc:.2f}s")
for k, v in res.items():
    v = np.asarray(v, dtype=float)
    if k.startswith("npk"):
        print(f"  {k:10s} centroided peaks median {np.median(v):.0f}")
    else:
        print(f"  {k:10s} b/y coverage median {np.median(v):.3f} mean {v.mean():.3f}  >=0.25: {np.mean(v >= 0.25):.1%}")
a = np.array(imfit, dtype=float)
ok = np.isfinite(a[:, 0])
if ok.sum() > 20:
    k, c0 = np.polyfit(a[ok, 0], a[ok, 1], 1)
    resid = a[ok, 1] - (k * a[ok, 0] + c0)
    print(f"  apex_scan = {k:.1f} * IM + {c0:.1f}; residual MAD {1.4826 * np.median(np.abs(resid)):.1f} scans, "
          f"|resid|<=16: {np.mean(np.abs(resid) <= 16):.1%}")
