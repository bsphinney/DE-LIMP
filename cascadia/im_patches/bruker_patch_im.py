"""
bruker_patch_im.py — IM-enhanced Bruker .d augmentation for Cascadia.

Uses timsrust_pyo3's SpectrumReader.new_with_span_step() to split each
DIA window into ~30 ion mobility slices of 0.05 1/K0 width. Each spectrum
gets a real precursor.im value (IM slice midpoint), enabling the model to
learn ion mobility-dependent fragmentation patterns.

Output: 5-column ASF lines: m/z  intensity  RT_offset  MS_level  1/K0

The 1/K0 for MS2 peaks = the IM slice midpoint from precursor.im
The 1/K0 for MS1 peaks = 0.0 (MS1 frames don't have IM filtering)
"""

import numpy as np
import os
import logging
import sqlite3

logger = logging.getLogger(__name__)


def _read_mz_calibration(d_folder):
    """Read m/z calibration parameters from analysis.tdf SQLite."""
    tdf_path = os.path.join(str(d_folder), "analysis.tdf")
    if not os.path.exists(tdf_path):
        return 100.0, 1700.0, 631025  # defaults

    try:
        conn = sqlite3.connect(f"file:{tdf_path}?mode=ro", uri=True)
        meta = dict(conn.execute("SELECT Key, Value FROM GlobalMetadata").fetchall())
        conn.close()
        mz_lower = float(meta.get('MzAcqRangeLower', 100.0))
        mz_upper = float(meta.get('MzAcqRangeUpper', 1700.0))
        dns = int(meta.get('DigitizerNumSamples', 631025))
        return mz_lower, mz_upper, dns
    except Exception as e:
        logger.warning("Failed to read TDF calibration: %s", e)
        return 100.0, 1700.0, 631025


def _tof_to_mz(tof_indices, mz_lower, mz_upper, dns):
    """Convert tof indices to m/z using the Sage/timsrust formula.

    From timsrust src/converters.rs (Tof2MzConverter):
      tof_intercept = sqrt(mz_min)
      tof_slope = (sqrt(mz_max) - sqrt(mz_min)) / tof_max_index
      mz = (tof_intercept + tof_slope * tof_index) ^ 2
    """
    tof = np.array(tof_indices, dtype=np.float64)
    tof_intercept = np.sqrt(mz_lower)
    tof_slope = (np.sqrt(mz_upper) - tof_intercept) / dns
    return (tof_intercept + tof_slope * tof) ** 2


def _read_ms1_frames(d_folder, mz_lower, mz_upper, dns, top_n=150):
    """Read MS1 frames and convert to (rt_min, mz_array, intensity_array) list.

    Applies the same normalization as augment.py MS1: sqrt(sqrt) + normalize.
    """
    import timsrust_pyo3 as timsrust

    reader = timsrust.FrameReader(str(d_folder))
    ms1_frames = reader.read_ms1_frames()
    logger.info("Read %d MS1 frames from %s", len(ms1_frames), os.path.basename(str(d_folder)))

    ms1_list = []
    for frame in ms1_frames:
        tof = np.array(frame.tof_indices, dtype=np.float64)
        ints = np.array(frame.intensities, dtype=np.float64)
        if len(tof) == 0:
            continue

        mzs = _tof_to_mz(tof, mz_lower, mz_upper, dns)

        # Sort by intensity, take top_n, re-sort by m/z (same as augment.py MS1)
        sorted_idx = np.argsort(ints)[-top_n:]
        ints = ints[sorted_idx]
        mzs = mzs[sorted_idx]
        sorted_mz = np.argsort(mzs)
        ints = ints[sorted_mz]
        mzs = mzs[sorted_mz]

        # Apply sqrt(sqrt) + normalize (same as augment.py MS1 path)
        ints = ints ** 0.5
        ints = ints ** 0.5
        if len(ints) > 0 and np.max(ints) > 0:
            ints = ints / np.max(ints)

        ms1_list.append((frame.rt, mzs, ints))  # rt is in minutes

    return ms1_list


def get_centers_bruker_im(d_folder, im_span=0.05, im_step=0.02):
    """Get precursor centers using IM-sliced spectra.

    Uses SpectrumReader.new_with_span_step() instead of read_all_spectra().
    Each spectrum now has a meaningful precursor.im value.

    Returns:
        f_to_mzrt_to_pep: dict of part -> dict of (mz_bin, rt_bin) -> list of (mz, rt, charge, im)
        max_mz: max m/z bin seen
        window_size: isolation window width
        cycle_time: estimated cycle time in seconds
    """
    import timsrust_pyo3 as timsrust

    f_to_mzrt_to_pep = {}
    max_mz = 0
    num_spectra = 0
    part = 0
    cycle_time = None
    window_size = None

    logger.info("Reading IM-sliced spectra from %s (span=%.2f, step=%.2f)",
                os.path.basename(str(d_folder)), im_span, im_step)
    reader = timsrust.SpectrumReader.new_with_span_step(str(d_folder), im_span, im_step)
    total = len(reader)
    logger.info("Total IM-sliced spectra: %d", total)

    # Estimate cycle time from consecutive spectra with different RTs
    prev_rt = None
    count_for_cycle = 0
    for spec in reader:
        rt_sec = float(spec.precursor.rt) * 60.0
        if prev_rt is not None and abs(rt_sec - prev_rt) > 0.01:
            cycle_time = rt_sec - prev_rt
            break
        prev_rt = rt_sec
        count_for_cycle += 1
        if count_for_cycle > 500:
            break

    # Re-create reader (iterators are consumed)
    reader = timsrust.SpectrumReader.new_with_span_step(str(d_folder), im_span, im_step)

    for spec in reader:
        if spec.precursor.mz is None or float(spec.precursor.mz) == 0.0:
            continue

        rt_sec = float(spec.precursor.rt) * 60.0
        window_center = float(spec.precursor.mz)
        im_value = float(spec.precursor.im) if hasattr(spec.precursor, 'im') and spec.precursor.im is not None else 0.0

        if window_size is None:
            if hasattr(spec, 'isolation_width') and spec.isolation_width:
                window_size = float(spec.isolation_width)
            else:
                window_size = 25.0

        if num_spectra % 50000 == 0:
            part += 1
            f_to_mzrt_to_pep[part] = {}
        num_spectra += 1

        key = (int(window_center / 10), int(rt_sec / 10))
        max_mz = max(max_mz, int(window_center / 10))
        # Store im_value as 4th element in the tuple (was charge=1)
        entry = (window_center, rt_sec, 1, im_value)
        if key in f_to_mzrt_to_pep[part]:
            f_to_mzrt_to_pep[part][key].append(entry)
        else:
            f_to_mzrt_to_pep[part][key] = [entry]

    if cycle_time is None:
        cycle_time = 1.0

    logger.info("Centers: %d IM-sliced MS2 spectra, max_mz=%d, window=%.1f, cycle=%.2fs",
                num_spectra, max_mz, window_size, cycle_time)
    return f_to_mzrt_to_pep, max_mz, window_size, cycle_time


def extract_spectra_bruker_im(d_folder, f_to_mzrt_to_pep, part, top_n, time_width, max_mz,
                               im_span=0.05, im_step=0.02):
    """Extract IM-enhanced augmented spectra from Bruker .d files.

    Uses SpectrumReader.new_with_span_step() for IM-sliced MS2 spectra.
    Each MS2 peak gets the IM slice midpoint as its 5th column value.
    MS1 peaks get 0.0 for IM (MS1 frames have no IM filtering).
    """
    import timsrust_pyo3 as timsrust

    prec_to_spec = {}

    # Read m/z calibration for MS1 frame conversion
    mz_lower, mz_upper, dns = _read_mz_calibration(d_folder)

    # Read MS1 frames with tof-to-mz conversion
    ms1_list = _read_ms1_frames(d_folder, mz_lower, mz_upper, dns, top_n)

    # Index MS1 by RT bin (10-second bins in seconds)
    ms1_by_rt = {}
    for rt_min, mzs, ints in ms1_list:
        rt_sec = rt_min * 60.0
        rt_bin = int(rt_sec / 10)
        # MS1 peaks: IM = 0.0 (no IM filtering for MS1 frames)
        mz_int_im = [(m, i, 0.0) for m, i in zip(mzs.tolist(), ints.tolist())]
        ms1_by_rt.setdefault(rt_bin, []).append((rt_sec, mz_int_im))

    logger.info("MS1 index: %d RT bins from %d frames", len(ms1_by_rt), len(ms1_list))

    # Read IM-sliced MS2 spectra
    reader = timsrust.SpectrumReader.new_with_span_step(str(d_folder), im_span, im_step)
    window_half = 12.5  # default half-width

    for spec in reader:
        if spec.precursor.mz is None or float(spec.precursor.mz) == 0.0:
            continue

        rt_sec = float(spec.precursor.rt) * 60.0
        window_center = float(spec.precursor.mz)
        im_value = float(spec.precursor.im) if hasattr(spec.precursor, 'im') and spec.precursor.im is not None else 0.0

        mzs = np.array(spec.mz_values, dtype=np.float64)
        intensities = np.array(spec.intensities, dtype=np.float64)

        if len(intensities) == 0:
            continue

        # Sort by intensity, take top_n, re-sort by m/z (same as augment.py MS2)
        sorted_idx = np.argsort(intensities)[-top_n:]
        intensities = intensities[sorted_idx]
        mzs = mzs[sorted_idx]
        sorted_mz = np.argsort(mzs)
        intensities = intensities[sorted_mz]
        mzs = mzs[sorted_mz]

        # Apply sqrt + normalize (augment.py MS2: single sqrt)
        intensities = intensities ** 0.5
        if np.max(intensities) > 0:
            intensities = intensities / np.max(intensities)

        # MS2 peaks carry the IM slice midpoint
        mz_int_im = [(m, i, im_value) for m, i in zip(mzs.tolist(), intensities.tolist())]

        lower_offset = window_half
        upper_offset = window_half

        # Match MS2 to precursors in window range
        for scan_rt in range(int(rt_sec / 10) - 1, int(rt_sec / 10) + 1):
            for scan_window in range(
                int((window_center - lower_offset) / 10) - 1,
                int((window_center + upper_offset) / 10) + 1
            ):
                key = (scan_window, scan_rt)
                if key not in f_to_mzrt_to_pep.get(part, {}):
                    continue
                for center_mz, center_rt, charge, center_im in f_to_mzrt_to_pep[part][key]:
                    in_mz = (center_mz > window_center - lower_offset and
                             center_mz < window_center + upper_offset)
                    if in_mz and abs(center_rt - rt_sec) < time_width:
                        prec_key = (center_mz, center_rt, charge, center_im)
                        if prec_key not in prec_to_spec:
                            prec_to_spec[prec_key] = {}
                        if 'scans' not in prec_to_spec[prec_key]:
                            prec_to_spec[prec_key]['scans'] = []
                            prec_to_spec[prec_key]['rts'] = []
                            prec_to_spec[prec_key]['window_width'] = max(lower_offset, upper_offset)
                            prec_to_spec[prec_key]['im'] = center_im
                        prec_to_spec[prec_key]['scans'].append(mz_int_im)
                        prec_to_spec[prec_key]['rts'].append(rt_sec - center_rt)

    # Now match MS1 context to each precursor
    n_with_ms1 = 0
    for prec_key, value in prec_to_spec.items():
        center_mz, center_rt, _, _ = prec_key
        rt_bin = int(center_rt / 10)

        for search_bin in range(rt_bin - 1, rt_bin + 2):
            for ms1_rt, ms1_triples in ms1_by_rt.get(search_bin, []):
                if abs(ms1_rt - center_rt) < time_width:
                    if 'ms1_scans' not in value:
                        value['ms1_scans'] = []
                        value['ms1_rts'] = []
                    value['ms1_scans'].append(ms1_triples)
                    value['ms1_rts'].append(ms1_rt - center_rt)

        if 'ms1_scans' in value and len(value['ms1_scans']) > 0:
            n_with_ms1 += 1

    # Filter out precursors without MS1 context
    prec_filtered = {k: v for k, v in prec_to_spec.items()
                     if 'ms1_scans' in v and len(v['ms1_scans']) > 0}

    logger.info("Precursors: %d total, %d with MS1 context (%.0f%%), returning %d",
                len(prec_to_spec), n_with_ms1,
                100 * n_with_ms1 / max(len(prec_to_spec), 1),
                len(prec_filtered))

    return prec_filtered


def write_asf_im(outfile, prec_to_spec, scan_width=1, max_pep_length=30, max_charge=3):
    """Write 5-column ASF with IM values.

    Format: m/z  intensity  RT_offset  MS_level  1/K0

    MS2 peaks: 1/K0 = IM slice midpoint from the precursor
    MS1 peaks: 1/K0 = 0.0 (no IM filtering)
    """
    out = open(outfile, 'a')
    skipped = 0
    count = 0
    for key, value in prec_to_spec.items():
        if 'ms1_scans' not in value:
            skipped += 1
            continue
        prec, rt, charge_orig, im_val = key
        scans = np.array(value['scans'], dtype=object)
        rts = np.array(value['rts'])
        ms1_scans = np.array(value['ms1_scans'], dtype=object)
        ms1_rts = np.array(value['ms1_rts'])
        window_width = value['window_width']

        abs_rts = [np.abs(x) for x in rts]
        sorted_rt_idxs = np.argsort(abs_rts)[:scan_width]
        rts = rts[sorted_rt_idxs]
        scans = scans[sorted_rt_idxs]

        abs_ms1_rts = [np.abs(x) for x in ms1_rts]
        sorted_ms1_rt_idxs = np.argsort(abs_ms1_rts)[:scan_width]
        ms1_rts = ms1_rts[sorted_ms1_rt_idxs]
        ms1_scans = ms1_scans[sorted_ms1_rt_idxs]

        for charge in range(2, max_charge + 1):
            count += 1

            out.write("BEGIN IONS\n")
            out.write(f"TITLE={count}\n")
            out.write(f"PEPMASS={prec}\n")
            out.write(f"CHARGE={charge}\n")
            out.write(f"SCAN={count}\n")
            out.write(f"RT={rt}\n")
            out.write(f"IM={im_val:.4f}\n")
            out.write(f"SEQ={'K' * max_pep_length}\n")

            # MS2 peaks: 5 columns with IM from the spectrum
            for scan, cur_rt in zip(scans, rts):
                for mz, intensity, peak_im in scan:
                    out.write(f"{mz}\t{intensity}\t{cur_rt}\t2\t{peak_im:.4f}\n")

            # MS1 peaks: 5 columns with IM = 0.0
            for scan, cur_rt in zip(ms1_scans, ms1_rts):
                for mz, intensity, peak_im in scan:
                    if np.abs(mz - prec) < window_width + 1:
                        out.write(f"{mz}\t{intensity}\t{cur_rt}\t1\t{peak_im:.4f}\n")

            out.write("END IONS\n")
    out.close()
    return count
