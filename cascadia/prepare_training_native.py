"""Labelled Cascadia training spectra from a diaPASEF .d, built exactly the way inference builds them.

Why this replaces prepare_training_data.py / prepare_training_data_mobfilter.py (Tracks A and B):
  1. Same input representation as inference (cascadia_native_run.py -> cascadia.augment on the
     bruker_native_augment records). At inference every MS2 window spectrum is a centre: PEPMASS is
     its isolation-window CENTRE (the peptide's own m/z is unknown), the centre scan has RT offset 0,
     the nearest MS1 frame within the time window contributes its top_n peaks (4th-root scaled)
     lying within +/-(half window + 1) of the centre, and MS2 peaks are the top_n, sqrt-scaled.
     Records here are the same, centred on the MS2 scan nearest the DIA-NN apex whose window holds
     the precursor. The decoder conditions on PEPMASS, so training on the true m/z would teach it
     to trust a mass it never gets at inference. The old generators wrote co-added MS2 only, no MS1,
     zero RT offsets.
  2. RT units. timsrust RT is seconds; the old generators treated it as minutes, so Tracks A and B
     trained on void-volume spectra (memory: radiant-generator-rt-units).
  3. True charge. The old writer emitted every precursor at charges 2..max_charge with one label.
  4. Modifications. Labels are DIA-NN Modified.Sequence in the ProForma names of Cascadia's
     MassIVE-KB tokenizer; precursors carrying any other mod are skipped.
  5. m/z calibration drift. The boundary TOF->m/z formula reproduces the instrument's stored
     calibration (Bruker libtimsdata, analysis.tdf MzCalibration) to within +/-3 ppm, but a run
     can drift from that calibration with nothing in the file to record it. Over 339 FRAN timsTOF
     runs the offset at m/z 700 has a median of -2.7 ppm (90% between -7.5 and +2.9), with
     outliers to -22 / +33 ppm; the mouse benchmark run is one of them (-25 ppm, calibrated with
     Tuning Mix the day before, 0.16 ppm SD). A per-run linear correction ppm(mz) = a + b*mz is
     fitted separately for MS2 (labelled b/y ions) and MS1 (labelled precursor m/z), applied to
     every peak, and saved in RUN.calib.json so inference can apply the same one
     (cascadia_native_run.py --calib).
  6. Peak budget. A diaPASEF window spectrum holds ~8,400 peaks; inference keeps 150, which drops
     most of a peptide's fragments (median b/y coverage 0.04 at 10 ppm, vs 0.22 at 1,000 peaks and
     0.46 at 3,000). --variants writes one ASF per top_n from the SAME labelled precursors, so the
     peak budget can be compared by fine-tuning on each.

Evidence gate: a precursor is kept only if >= --min-coverage of its b/y ions (1+/2+, counted once
per ion) are within --gate-ppm in the calibrated nearest scan's top --gate-top peaks. This asks
"are the fragments in the raw data", independent of the variant's truncation.

    python prepare_training_native.py --report report.parquet --run RUN --d RUN.d --out OUTPREFIX
           [--variants 150,1000] [--q 0.01] [--max-per-run 1000] [--min-coverage 0.25] [--seed 0]
           [--calib-only]
Writes OUTPREFIX.top<N>.asf per variant, OUTPREFIX.calib.json and OUTPREFIX.json (counts).
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import re
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bruker_native_augment import read_bruker_records  # noqa: E402

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
log = logging.getLogger("prep_native")

# DIA-NN UniMod -> MassIVE-KB tokenizer residue (cascadia/ptms.json + PeptideTokenizer.from_massivekb)
_RESIDUE_MODS = {("M", 35): "M[Oxidation]", ("C", 4): "C[Carbamidomethyl]",
                 ("N", 7): "N[Deamidated]", ("Q", 7): "Q[Deamidated]"}
_NTERM_MODS = {1: "[Acetyl]-", 5: "[Carbamyl]-"}
_NTERM_MASS = {"[Acetyl]-": 42.01057, "[Carbamyl]-": 43.00581}

_AA = {"G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841, "T": 101.04768,
       "C": 103.00919, "L": 113.08406, "I": 113.08406, "N": 114.04293, "D": 115.02694,
       "Q": 128.05858, "K": 128.09496, "E": 129.04259, "M": 131.04049, "H": 137.05891,
       "F": 147.06841, "R": 156.10111, "Y": 163.06333, "W": 186.07931}
_MOD_MASS = {4: 57.02146, 35: 15.99491, 7: 0.98402, 1: 42.01057, 5: 43.00581}
_PROTON, _WATER = 1.007276, 18.010565


def to_proforma(modseq: str):
    """DIA-NN Modified.Sequence -> (Cascadia label, residue masses) or (None, reason)."""
    s = modseq or ""
    nterm = ""
    m = re.match(r"^((?:\(UniMod:\d+\))+)", s)
    if m:                                    # DIA-NN writes an N-terminal mod before residue 1
        ids = [int(x) for x in re.findall(r"UniMod:(\d+)", m.group(1))]
        if len(ids) != 1 or ids[0] not in _NTERM_MODS:
            return None, "nterm_mod"
        nterm = _NTERM_MODS[ids[0]]
        s = s[m.end():]
    out, masses = [], []
    for aa, mods in re.findall(r"([A-Z])((?:\(UniMod:\d+\))*)", s):
        ids = [int(x) for x in re.findall(r"UniMod:(\d+)", mods)]
        if len(ids) > 1:
            return None, "stacked_mod"
        if ids:
            tok = _RESIDUE_MODS.get((aa, ids[0]))
            if tok is None:
                if ids[0] == 1 and not out and not nterm:   # acetyl written on residue 1
                    nterm = "[Acetyl]-"
                    out.append(aa)
                    masses.append(_AA[aa])          # acetyl mass is added as nterm_add
                    continue
                return None, f"mod_{ids[0]}"
            out.append(tok)
            masses.append(_AA[aa] + _MOD_MASS[ids[0]])
        else:
            if aa not in _AA:
                return None, "residue"
            out.append(aa)
            masses.append(_AA[aa])
    return nterm + "".join(out), masses


def by_ions(masses, nterm_add=0.0):
    """Theoretical b/y m/z, shape (n_ions, 2): columns 1+ and 2+."""
    theo = []
    b = nterm_add
    for i in range(len(masses) - 1):
        b += masses[i]
        theo.append((b + _PROTON, (b + 2 * _PROTON) / 2))
    y = _WATER
    for i in range(len(masses) - 1, 0, -1):
        y += masses[i]
        theo.append((y + _PROTON, (y + 2 * _PROTON) / 2))
    return np.asarray(theo, dtype=float).reshape(-1, 2)


def by_coverage(masses, mzs, nterm_add=0.0, ppm=20.0):
    """Fraction of theoretical b/y ions (1+ or 2+, an ion counted once) present in mzs."""
    if len(masses) < 2 or len(mzs) == 0:
        return 0.0
    s = np.sort(np.asarray(mzs, dtype=float))
    theo = by_ions(masses, nterm_add).ravel()
    idx = np.clip(np.searchsorted(s, theo), 1, len(s) - 1)
    near = np.minimum(np.abs(s[idx] - theo), np.abs(s[idx - 1] - theo))
    return float((near <= theo * ppm * 1e-6).reshape(-1, 2).any(axis=1).mean())


def _signed_ppm(mzs, theo, window_ppm):
    """Signed ppm error of the peak nearest each theoretical m/z, for those within window_ppm."""
    s = np.sort(np.asarray(mzs, dtype=float))
    if s.size < 2:
        return np.empty(0), np.empty(0)
    theo = np.asarray(theo, dtype=float)
    idx = np.clip(np.searchsorted(s, theo), 1, len(s) - 1)
    lo, hi = s[idx - 1], s[idx]
    best = np.where(np.abs(hi - theo) < np.abs(lo - theo), hi, lo)
    e = (best - theo) / theo * 1e6
    ok = np.abs(e) <= window_ppm
    return e[ok], theo[ok]


def fit_ppm(errs, mzs, window=80.0, band=8.0):
    """Linear ppm(mz) = a + b*mz through the mode of the error histogram. Random matches spread
    flat over +/-window; true matches pile up at the systematic offset."""
    errs, mzs = np.asarray(errs), np.asarray(mzs)
    if errs.size < 50:
        return {"a": 0.0, "b": 0.0, "n": int(errs.size), "fitted": False}
    h, edges = np.histogram(errs, bins=np.arange(-window, window + 2, 2))
    mode = float(edges[np.argmax(h)] + 1)
    sel = np.abs(errs - mode) <= band
    a, b = mode, 0.0
    if sel.sum() >= 50 and np.ptp(mzs[sel]) > 100:
        b, a = np.polyfit(mzs[sel], errs[sel], 1)
        resid = errs[sel] - (a + b * mzs[sel])
        keep = np.abs(resid) <= 3 * max(1.0, 1.4826 * np.median(np.abs(resid)))
        b, a = np.polyfit(mzs[sel][keep], errs[sel][keep], 1)
    resid = errs[sel] - (a + b * mzs[sel])
    return {"a": float(a), "b": float(b), "mode": mode, "n": int(sel.sum()), "n_window": int(errs.size),
            "mad_ppm": float(1.4826 * np.median(np.abs(resid))) if sel.any() else None, "fitted": True}


def correct(mz, coef):
    """Apply a fitted ppm(mz) correction: true ~= observed * (1 - ppm(observed)*1e-6)."""
    mz = np.asarray(mz, dtype=np.float64)
    return mz * (1.0 - (coef["a"] + coef["b"] * mz) * 1e-6)


def _top(mz, it, n):
    if n is not None and it.size > n:
        k = np.argpartition(it, -n)[-n:]
        mz, it = mz[k], it[k]
    o = np.argsort(mz)
    return mz[o], it[o]


def _augment_peaks(mz, it, n, root):
    """cascadia.augment's per-scan transform: top-n by intensity, m/z order, intensity**root, /max."""
    mz, it = _top(mz, it, n)
    it = it ** root
    if it.size:
        it = it / it.max()
    return mz, it


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", required=True)
    ap.add_argument("--run", required=True, help="DIA-NN Run name")
    ap.add_argument("--d", required=True, help=".d folder (opened read-only by timsrust)")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--variants", default="150,1000", help="comma list of top_n to write")
    ap.add_argument("--q", type=float, default=0.01)
    ap.add_argument("--max-per-run", type=int, default=1000)
    ap.add_argument("--min-coverage", type=float, default=0.25)
    ap.add_argument("--gate-top", type=int, default=3000)
    ap.add_argument("--own-min-coverage", type=float, default=0.0,
                    help="per variant, keep a record only if this fraction of its b/y ions is in its own top-N")
    ap.add_argument("--gate-ppm", type=float, default=10.0)
    ap.add_argument("--scan-width", type=int, default=1, help="cycles of RT window (augment: (w+1)*cycle)")
    ap.add_argument("--max-pep-length", type=int, default=30)
    ap.add_argument("--max-charge", type=int, default=4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--calib-only", action="store_true", help="fit and write OUT.calib.json only")
    a = ap.parse_args()
    t0 = time.time()
    variants = sorted({int(x) for x in a.variants.split(",") if x.strip()})
    import pyarrow.parquet as pq
    import pyarrow.compute as pc

    cols = ["Run", "Modified.Sequence", "Precursor.Charge", "Precursor.Mz", "RT", "Q.Value", "Decoy"]
    have = set(pq.read_schema(a.report).names)
    t = pq.read_table(a.report, columns=[c for c in cols if c in have], filters=[("Run", "==", a.run)])
    t = t.filter(pc.less_equal(t["Q.Value"], a.q))
    if "Decoy" in have:
        t = t.filter(pc.equal(t["Decoy"], 0))
    df = t.to_pandas()
    stats = {"run": a.run, "report_rows_q": int(len(df)), "skipped": {}, "variants": variants}

    def skip(r, k=1):
        stats["skipped"][r] = stats["skipped"].get(r, 0) + k

    # one label per precursor (modseq, z): its most confident row
    df = df.sort_values("Q.Value").drop_duplicates(["Modified.Sequence", "Precursor.Charge"])
    rng = np.random.default_rng(a.seed)
    labels = []
    for ms, z, mz, rt in zip(df["Modified.Sequence"], df["Precursor.Charge"], df["Precursor.Mz"], df["RT"]):
        z = int(z)
        if not 2 <= z <= a.max_charge:
            skip("charge")
            continue
        lab, masses = to_proforma(ms)
        if lab is None:
            skip(masses)
            continue
        if len(masses) > a.max_pep_length:
            skip("too_long")
            continue
        nadd = next((v for k, v in _NTERM_MASS.items() if lab.startswith(k)), 0.0)
        labels.append((float(mz), float(rt) * 60.0, z, lab, masses, nadd))
    if not labels:
        log.error("no usable labels for run %s", a.run)
        return 2
    order = rng.permutation(len(labels))
    cand = [labels[i] for i in order[: max(a.max_per_run * 3, 3000)]]   # gate drops some; calib wants many
    stats["candidates"] = len(cand)

    recs = read_bruker_records(a.d, top_n=max(a.gate_top, max(variants)))
    ms2 = [r for r in recs if r["ms level"] == 2]
    ms1 = [r for r in recs if r["ms level"] == 1]
    del recs
    t2 = np.array([60 * r["scanList"]["scan"][0]["scan start time"] for r in ms2])
    c2 = np.array([r["precursorList"]["precursor"][0]["isolationWindow"]["isolation window target m/z"] for r in ms2])
    h2 = np.array([r["precursorList"]["precursor"][0]["isolationWindow"]["isolation window upper offset"] for r in ms2])
    t1 = np.array([60 * r["scanList"]["scan"][0]["scan start time"] for r in ms1])
    cycle = float(np.median(np.diff(t1))) if len(t1) > 2 else 1.0
    tw = (a.scan_width + 1) * cycle                     # augment's time_width
    o2 = np.argsort(t2, kind="stable")

    def nearest_ms2(mz, rt):
        lo, hi = np.searchsorted(t2[o2], [rt - tw, rt + tw])
        idx = o2[lo:hi]
        idx = idx[(c2[idx] - h2[idx] < mz) & (mz < c2[idx] + h2[idx])]
        if idx.size == 0:
            return None
        d = np.abs(t2[idx] - rt) + 1e-6 * np.abs(c2[idx] - mz)     # tie: the window centred nearer
        return int(idx[np.argmin(d)])

    def nearest_ms1(t):
        j = int(np.argmin(np.abs(t1 - t)))
        return j if abs(t1[j] - t) < tw else None

    # pair every candidate with its centre scan and MS1 frame
    paired = []
    for p in cand:
        i = nearest_ms2(p[0], p[1])
        if i is None:
            skip("no_ms2")
            continue
        j = nearest_ms1(t2[i])
        if j is None:
            skip("no_ms1")
            continue
        paired.append((p, i, j))

    # --- calibration: MS2 from labelled b/y ions, MS1 from labelled precursor m/z (monoisotopic)
    e2, m2, e1, m1 = [], [], [], []
    for p, i, j in paired:
        mz_g, _ = _top(np.asarray(ms2[i]["m/z array"]), np.asarray(ms2[i]["intensity array"]), a.gate_top)
        e, th = _signed_ppm(mz_g, by_ions(p[4], p[5])[:, 0], 80.0)
        e2.append(e); m2.append(th)
        e, th = _signed_ppm(np.asarray(ms1[j]["m/z array"]), [p[0]], 80.0)
        e1.append(e); m1.append(th)
    cal2 = fit_ppm(np.concatenate(e2), np.concatenate(m2))
    cal1 = fit_ppm(np.concatenate(e1), np.concatenate(m1))
    calib = {"run": a.run, "ms2": cal2, "ms1": cal1, "model": "true = obs * (1 - (a + b*obs)*1e-6)"}
    json.dump(calib, open(a.out + ".calib.json", "w"), indent=1)
    log.info("calibration ms2 %s | ms1 %s", json.dumps(cal2), json.dumps(cal1))
    stats["calib"] = calib
    if a.calib_only:
        return 0

    # --- evidence gate on the calibrated top gate_top peaks of the centre scan
    kept = []
    for p, i, j in paired:
        mz_g, _ = _top(correct(ms2[i]["m/z array"], cal2), np.asarray(ms2[i]["intensity array"]), a.gate_top)
        if by_coverage(p[4], mz_g, p[5], a.gate_ppm) < a.min_coverage:
            skip("low_by_coverage")
            continue
        kept.append((p, i, j))
        if len(kept) >= a.max_per_run:
            break
    stats["gated"] = len(kept)

    # --- write one ASF per variant, same precursors, inference's own record construction.
    # --own-min-coverage additionally drops, per variant, records whose b/y ions are not in that
    # variant's own top-N peaks (the top-3000 gate alone left top-150 records with ~4% coverage).
    per_variant = {}
    for v in variants:
        path = f"{a.out}.top{v}.asf"
        n = 0
        with open(path + ".tmp", "w") as out:
            for p, i, j in kept:
                prec_c, half, t_c = float(c2[i]), float(h2[i]), float(t2[i])
                mz2, it2 = _augment_peaks(correct(ms2[i]["m/z array"], cal2),
                                          np.asarray(ms2[i]["intensity array"], dtype=float), v, 0.5)
                if a.own_min_coverage > 0 and by_coverage(p[4], mz2, p[5], a.gate_ppm) < a.own_min_coverage:
                    continue
                n += 1
                mz1, it1 = _augment_peaks(correct(ms1[j]["m/z array"], cal1),
                                          np.asarray(ms1[j]["intensity array"], dtype=float), v, 0.25)
                d1 = float(t1[j] - t_c)
                out.write(f"BEGIN IONS\nTITLE={a.run}:{n}\nPEPMASS={prec_c}\nCHARGE={p[2]}\nSCAN={n}\n"
                          f"RT={t_c}\nSEQ={p[3]}\n")
                out.writelines(f"{x}\t{y}\t0.0\t2\n" for x, y in zip(mz2, it2))
                keep1 = np.abs(mz1 - prec_c) < half + 1
                out.writelines(f"{x}\t{y}\t{d1}\t1\n" for x, y in zip(mz1[keep1], it1[keep1]))
                out.write("END IONS\n")
        os.replace(path + ".tmp", path)
        per_variant[v] = n
    stats.update(written=len(kept), written_per_variant=per_variant, cycle_s=round(cycle, 4), time_width_s=round(tw, 4),
                 rt_s_min=round(min(t2[i] for _, i, _ in kept), 1) if kept else None,
                 rt_s_max=round(max(t2[i] for _, i, _ in kept), 1) if kept else None,
                 seconds=round(time.time() - t0))
    json.dump(stats, open(a.out + ".json", "w"), indent=1)
    log.info(json.dumps({k: v for k, v in stats.items() if k != "calib"}))
    return 0 if kept else 2


if __name__ == "__main__":
    sys.exit(main())
