"""CPU checks for bruker_native_augment on one diaPASEF .d with a DIA-NN report.

Usage: check_bruker_native.py <run.d> <report.parquet> <tmpdir>
1. MS1 m/z accuracy: for DIA-NN precursors (q<=0.01) of this run, the nearest MS1 centroid in the MS1
   frame closest to the precursor's RT; median ppm and share within 20 ppm.
2. Builds the augmented ASF (max_charge 4) and reports record counts, peaks per level, RT offsets.
"""
import os
import sys
import time

import numpy as np
import pyarrow.parquet as pq

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import timsrust_pyo3 as timsrust  # noqa: E402

from bruker_native_augment import (_tdf_meta, augment_spectra_bruker_native,  # noqa: E402
                                   centroid_ms1_frame, tof_to_mz)

d, report, tmp = sys.argv[1:4]
run = os.path.basename(d.rstrip("/"))[:-2]
cols = [c for c in ("Run", "Precursor.Mz", "RT", "Q.Value", "Decoy") if c in pq.read_schema(report).names]
df = pq.read_table(report, columns=cols).to_pandas()
df = df[(df["Run"] == run) & (df["Q.Value"] <= 0.01)]
if "Decoy" in df:
    df = df[df["Decoy"] == 0]
print(f"DIA-NN q<=0.01 precursors in {run}: {len(df)}")

lo, hi, dns = _tdf_meta(d)
conv = lambda t: tof_to_mz(t, lo, hi, dns)  # noqa: E731
t0 = time.time()
frames = timsrust.FrameReader(d).read_ms1_frames()
cents = [centroid_ms1_frame(f.tof_indices, f.intensities, conv) for f in frames]
rts = np.array([f.rt for f in frames])
print(f"MS1: {len(frames)} frames, median {np.median([c[0].size for c in cents]):.0f} centroids/frame, "
      f"{time.time()-t0:.0f}s")
samp = df.sample(min(3000, len(df)), random_state=0)
ppm = []
for mz, rt in zip(samp["Precursor.Mz"], samp["RT"]):
    i = int(np.argmin(np.abs(rts - rt * 60.0)))
    cm = cents[i][0]
    if cm.size:
        j = int(np.argmin(np.abs(cm - mz)))
        ppm.append((cm[j] - mz) / mz * 1e6)
ppm = np.array(ppm)
near = ppm[np.abs(ppm) < 50]
print(f"MS1 m/z vs DIA-NN Precursor.Mz: {near.size}/{ppm.size} within 50 ppm; median {np.median(near):+.1f} ppm, "
      f"IQR {np.percentile(near,25):+.1f}..{np.percentile(near,75):+.1f}; within 20 ppm: {np.mean(np.abs(ppm)<20):.1%}")

t0 = time.time()
asf, win, cyc = augment_spectra_bruker_native(d, tmp, max_charge=4)
n = n1 = n2 = 0
off1, off2 = [], []
with open(asf) as f:
    for line in f:
        if line.startswith("BEGIN IONS"):
            n += 1
        elif "\t" in line:
            p = line.split("\t")
            if p[3].strip() == "1":
                n1 += 1
                if len(off1) < 200000:
                    off1.append(float(p[2]))
            else:
                n2 += 1
                if len(off2) < 200000:
                    off2.append(float(p[2]))
off1, off2 = np.abs(off1), np.abs(off2)
print(f"ASF: {n} records ({os.path.getsize(asf)/1e9:.2f} GB) in {time.time()-t0:.0f}s; window {win} m/z, cycle {cyc:.3f} s")
print(f"  peaks/record: MS2 {n2/n:.0f}, MS1 {n1/n:.0f}")
print(f"  |RT offset| s: MS2 mean {off2.mean():.3f} (nonzero {np.mean(off2>0):.0%}), MS1 mean {off1.mean():.3f} max {off1.max():.3f}")
print("first record:")
with open(asf) as f:
    for k, line in enumerate(f):
        if k > 12:
            break
        print("  " + line.rstrip())
