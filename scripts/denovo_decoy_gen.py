#!/usr/bin/env python
"""NovoBoard-style decoy SPECTRA generator for de novo target-decoy FDR.

Tran et al. 2024, Mol Cell Proteomics (doi:10.1016/j.mcpro.2024.100849), "random peak
removal" strategy: for each real MS2 spectrum, remove FRAC of its peaks and add back the
same number of peaks sampled from the GLOBAL pool of all peaks (so each noise peak is a
real fragment ion from some peptide). Precursor m/z + charge are PRESERVED, so Casanovo
sees a realistic-statistics spectrum that corresponds to no real peptide.

VALIDATED PARAMETER (DE-LIMP, 2026-06-03 HeLa entrapment): use FRAC >= 0.8. At FRAC=0.5 the
most abundant real peptides (GAPDH, beta-actin, tubulin) are still reconstructible and leak
into the null with top BLAST scores; at >=0.8 the null is clean. See docs/DENOVO_FDR_VALIDATION.md.

Usage: denovo_decoy_gen.py <input.mzML|.mgf> <out.mgf> [FRAC=0.8] [SEED=1337]
"""
import sys, re
import numpy as np

INPUT = sys.argv[1]
OUTMGF = sys.argv[2]
FRAC = float(sys.argv[3]) if len(sys.argv) > 3 else 0.8
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 1337
rng = np.random.default_rng(SEED)


def read_spectra(path):
    """Yield (scan, precursor_mz, charge, mz[], intensity[]) for MS2 spectra."""
    low = path.lower()
    if low.endswith(".mgf"):
        from pyteomics import mgf
        for i, s in enumerate(mgf.read(path)):
            mz = np.asarray(s["m/z array"], float); it = np.asarray(s["intensity array"], float)
            if len(mz) == 0:
                continue
            p = s.get("params", {})
            pmz = (p.get("pepmass") or ("",))[0]
            ch = ""
            if p.get("charge"):
                try: ch = int(p["charge"][0])
                except Exception: ch = ""
            m = re.search(r"scan=(\d+)", str(p.get("title", "")))
            scan = int(m.group(1)) if m else (p.get("scans") or i)
            yield scan, pmz, ch, mz, it
    else:
        from pyteomics import mzml
        for i, s in enumerate(mzml.read(path)):
            if s.get("ms level") != 2:
                continue
            mz = np.asarray(s["m/z array"], float); it = np.asarray(s["intensity array"], float)
            if len(mz) == 0:
                continue
            m = re.search(r"scan=(\d+)", s.get("id", ""))
            scan = int(m.group(1)) if m else i
            pmz = ""; ch = ""
            try:
                pr = s["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]
                pmz = pr.get("selected ion m/z", ""); ch = pr.get("charge state", "")
            except Exception:
                pass
            yield scan, pmz, ch, mz, it


# Pass 1: collect spectra + global peak pool ((m/z, intensity) PAIRS, so swapped-in
# noise peaks carry a real fragment-ion intensity from the pool — not a flat value).
specs = []
pool_mz = []; pool_int = []
for scan, pmz, ch, mz, it in read_spectra(INPUT):
    specs.append((scan, pmz, ch, mz, it)); pool_mz.append(mz); pool_int.append(it)
pool_mz  = np.concatenate(pool_mz)  if pool_mz  else np.array([])
pool_int = np.concatenate(pool_int) if pool_int else np.array([])
print(f"MS2 spectra: {len(specs)} | global peak pool: {len(pool_mz)} | FRAC={FRAC}", flush=True)

# Pass 2: write decoy MGF
with open(OUTMGF, "w") as f:
    written = 0
    for scan, pmz, ch, mz, it in specs:
        npk = len(mz); n_rm = int(round(FRAC * npk))
        keep = rng.permutation(npk)[n_rm:]
        kmz, kit = mz[keep], it[keep]
        if n_rm > 0 and len(pool_mz):
            sidx = rng.integers(0, len(pool_mz), size=n_rm)
            dmz = np.concatenate([kmz, pool_mz[sidx]])
            dit = np.concatenate([kit, pool_int[sidx]])   # real pooled intensities, not a flat fill
        else:
            dmz, dit = kmz, kit
        order = np.argsort(dmz); dmz, dit = dmz[order], dit[order]
        f.write("BEGIN IONS\n")
        f.write(f"TITLE=decoy_scan={scan}\n")
        if ch != "" and ch is not None:
            f.write(f"CHARGE={int(ch)}+\n")
        if pmz != "" and pmz is not None:
            f.write(f"PEPMASS={float(pmz):.6f}\n")
        f.write(f"SCANS={scan}\n")
        for a, b in zip(dmz, dit):
            f.write(f"{a:.5f} {b:.2f}\n")
        f.write("END IONS\n")
        written += 1
print(f"wrote {written} decoy spectra -> {OUTMGF}", flush=True)
