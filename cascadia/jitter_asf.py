"""Concatenate training ASFs, shifting every record's m/z by one random ppm offset (drift augmentation).

A run's calibration can drift from the instrument's stored calibration by several ppm (median
-2.7, outliers -22/+33 ppm over 339 FRAN timsTOF runs), and nothing in the .d records it. Records
from prepare_training_native.py are calibrated per run, so without augmentation the model only
ever sees ~0 ppm error and nothing teaches it to tolerate drift. Each record gets one offset drawn
uniformly from +/-PPM, applied to MS2 and MS1 peaks alike (drift moves both levels together:
median MS1-MS2 difference 0.9 ppm). PEPMASS is the isolation-window centre and is left alone.

    jitter_asf.py OUT.asf PPM SEED IN1.asf [IN2.asf ...]
"""
import random
import sys

out, ppm, seed, ins = sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4:]
rng = random.Random(seed)
n = 0
with open(out + ".tmp", "w") as o:
    for p in ins:
        with open(p) as f:
            k = 1.0
            for line in f:
                c0 = line[0]
                if c0.isdigit():                       # peak line: mz \t intensity \t rt \t level
                    mz, rest = line.split("\t", 1)
                    o.write(f"{float(mz) * k}\t{rest}")
                    continue
                if line.startswith("BEGIN IONS"):
                    k = 1.0 + rng.uniform(-ppm, ppm) * 1e-6
                    n += 1
                o.write(line)
import os
os.replace(out + ".tmp", out)
print(f"{n} records from {len(ins)} files, +/-{ppm} ppm -> {out}")
