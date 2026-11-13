"""
Prepare Cascadia training data from DIA-NN searched timsTOF .d files.

Reads DIA-NN report.parquet to get peptide-spectrum matches,
reads .d files via timsrust_pyo3 (MS2 spectra) and SQLite (MS1 frames),
creates labeled ASF files for Cascadia fine-tuning.

v3: Fixed m/z matching for DIA windows (search ±30 Da bins),
    added b/y ion quality filter per Cascadia paper (≥25% theoretical ions),
    correct Sage tof-to-mz formula, batch size 16 per paper.
"""

import argparse
import logging
import os
import sys
import time
import sqlite3
import numpy as np
# b/y ion filter uses manual AA masses (no external dependency)

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(message)s')
logger = logging.getLogger("prepare_training")


def compute_theoretical_ions(peptide, max_charge=2):
    """Compute theoretical b and y ion m/z values for a peptide.

    Returns set of (m/z rounded to 1 decimal) for matching.
    Per Cascadia paper: check charge +1 and +2 b- and y-ions.
    """
    ions = set()
    try:
        aa_masses = {
            'G': 57.02146, 'A': 71.03711, 'V': 99.06841, 'L': 113.08406,
            'I': 113.08406, 'P': 97.05276, 'F': 147.06841, 'W': 186.07931,
            'M': 131.04049, 'S': 87.03203, 'T': 101.04768, 'C': 160.03065,  # +CAM
            'Y': 163.06333, 'H': 137.05891, 'D': 115.02694, 'E': 129.04259,
            'N': 114.04293, 'Q': 128.05858, 'K': 128.09496, 'R': 156.10111,
        }
        proton = 1.007276

        # Calculate cumulative masses for b and y ions
        masses = [aa_masses.get(aa, 0) for aa in peptide]
        total = sum(masses) + 18.01056  # add water

        # b ions: cumulative from N-terminus
        cum = 0
        for i in range(len(masses) - 1):
            cum += masses[i]
            for z in range(1, max_charge + 1):
                ions.add(round((cum + z * proton) / z, 1))

        # y ions: cumulative from C-terminus
        cum = 18.01056  # water
        for i in range(len(masses) - 1, 0, -1):
            cum += masses[i]
            for z in range(1, max_charge + 1):
                ions.add(round((cum + z * proton) / z, 1))
    except Exception:
        pass
    return ions


def check_fragment_coverage(fragments, peptide, min_coverage=0.10):
    """Check if spectrum contains sufficient theoretical b/y ions.

    Based on Cascadia paper (25% for Orbitrap), but relaxed to 10% for
    timsTOF diaPASEF where DIA windows are wider and more chimeric.
    Uses ±0.5 Da tolerance to account for tof-to-mz calibration error.
    """
    theoretical = compute_theoretical_ions(peptide, max_charge=2)
    if not theoretical:
        return True

    # Build set of observed m/z values with ±0.5 Da tolerance
    frag_mzs = set()
    for mz, _ in fragments:
        # Round to 0.5 Da bins for matching
        frag_mzs.add(round(mz * 2) / 2)

    n_matched = 0
    for ion_mz in theoretical:
        ion_bin = round(ion_mz * 2) / 2
        if ion_bin in frag_mzs or (ion_bin - 0.5) in frag_mzs or (ion_bin + 0.5) in frag_mzs:
            n_matched += 1

    coverage = n_matched / len(theoretical) if theoretical else 0
    return coverage >= min_coverage


def load_diann_report(report_path, q_threshold=0.01):
    """Load DIA-NN report and extract high-confidence peptide-spectrum matches."""
    import pyarrow.parquet as pq

    logger.info("Loading DIA-NN report: %s", report_path)
    df = pq.read_table(report_path).to_pandas()
    logger.info("Total rows: %d", len(df))

    if 'Q.Value' in df.columns:
        df = df[df['Q.Value'] <= q_threshold]
        logger.info("After Q.Value <= %s filter: %d rows", q_threshold, len(df))

    required = ['Stripped.Sequence', 'Precursor.Mz', 'RT']
    missing = [c for c in required if c not in df.columns]
    if missing:
        logger.error("Missing required columns: %s", missing)
        sys.exit(1)

    if 'File.Name' in df.columns:
        file_col = 'File.Name'
    elif 'Run' in df.columns:
        file_col = 'Run'
    else:
        file_col = df.columns[0]

    logger.info("Unique peptides: %d", df['Stripped.Sequence'].nunique())
    logger.info("Unique files: %d", df[file_col].nunique())

    return df, file_col


def read_ms1_from_tdf(d_path):
    """Read MS1 frames directly from analysis.tdf SQLite with tof-to-mz calibration.

    Returns list of (rt_minutes, mz_array, intensity_array) tuples.
    """
    tdf_path = os.path.join(d_path, "analysis.tdf")
    if not os.path.exists(tdf_path):
        logger.warning("No analysis.tdf in %s", d_path)
        return []

    conn = sqlite3.connect(tdf_path)

    # Get calibration parameters for tof-to-mz conversion
    # TimsTOF uses: mz = c0 + c1/t + c2/t^2 ... or similar
    # The calibration is stored in the GlobalMetadata table
    try:
        meta = dict(conn.execute(
            "SELECT Key, Value FROM GlobalMetadata"
        ).fetchall())

        # Get MZ calibration model
        # Bruker TDF uses: mz = DigitizerTimebase * (C0 + C1*tof + C2*tof^2 + ...)
        # But the simpler approach: read from the Frames table which has
        # pre-calibrated data via the TDF-SDK
        # Actually, timsrust_pyo3 FrameReader does the conversion internally
        # when using read_ms1_frames — the issue is just that Frame objects
        # expose tof_indices, not mz_values.

        # Alternative: read MS1 data from the PropertyGroup/PasefFrameMsMsInfo tables
        # But simplest: just use the Frames table directly with numpy conversion

    except Exception as e:
        logger.warning("Failed to read TDF metadata: %s", e)
        conn.close()
        return []

    # Read MS1 frame IDs and RTs
    try:
        ms1_data = conn.execute("""
            SELECT f.Id, f.Time, f.NumPeaks
            FROM Frames f
            WHERE f.MsMsType = 0
            ORDER BY f.Time
        """).fetchall()
        logger.info("Found %d MS1 frames in TDF", len(ms1_data))
    except Exception as e:
        logger.warning("Failed to read MS1 frames: %s", e)
        conn.close()
        return []

    conn.close()

    # Use timsrust FrameReader to get MS1 frames with proper calibration
    import timsrust_pyo3 as timsrust
    reader = timsrust.FrameReader(str(d_path))
    ms1_frames = reader.read_ms1_frames()

    # FrameReader returns Frame objects with tof_indices + intensities
    # We need to convert tof_indices to m/z using the calibration
    # Since Frame doesn't have mz_values, we need to read from TDF directly
    # Use the GlobalMetadata calibration coefficients

    # Read calibration from TDF
    conn = sqlite3.connect(tdf_path)
    try:
        # Try to get tof-to-mz conversion via property definitions
        # Bruker stores calibration as: mz = (C1 / (10000 * sqrt(tof/C2) - C3))^2
        # Where C1, C2, C3 are per-frame or global calibration constants
        # These are in the Frames table: TimsCalibration column

        # Simpler approach: use the bin-to-mz lookup if available
        # Actually, the standard approach is:
        # Read from Frames table the MzCalibration columns

        # Get calibration for first frame
        calib_data = conn.execute("""
            SELECT Value FROM GlobalMetadata WHERE Key = 'MzAcqRangeLower'
        """).fetchone()
        mz_lower = float(calib_data[0]) if calib_data else 100.0

        calib_data = conn.execute("""
            SELECT Value FROM GlobalMetadata WHERE Key = 'MzAcqRangeUpper'
        """).fetchone()
        mz_upper = float(calib_data[0]) if calib_data else 1700.0

        # Get the DigitizerNumSamples for tof range
        calib_data = conn.execute("""
            SELECT Value FROM GlobalMetadata WHERE Key = 'DigitizerNumSamples'
        """).fetchone()
        n_tof = int(calib_data[0]) if calib_data else 400000

    except Exception as e:
        logger.warning("Failed to read calibration: %s", e)
        mz_lower, mz_upper, n_tof = 100.0, 1700.0, 400000
    finally:
        conn.close()

    # Convert tof indices to approximate m/z using linear interpolation
    # This is an approximation — real calibration is nonlinear (sqrt)
    # But for MS1 context in training data, approximate is sufficient
    # The model uses MS1 for precursor isotope patterns, not exact mass accuracy

    ms1_list = []
    for frame in ms1_frames:
        rt_min = frame.rt  # already in minutes
        tof_idx = np.array(frame.tof_indices, dtype=np.float64)
        ints = np.array(frame.intensities, dtype=np.float64)

        if len(tof_idx) == 0:
            continue

        # timsrust/Sage formula (from src/converters.rs Tof2MzConverter):
        #   tof_intercept = sqrt(mz_min)
        #   tof_slope = (sqrt(mz_max) - sqrt(mz_min)) / DigitizerNumSamples
        #   mz = (tof_intercept + tof_slope * tof_index)^2
        tof_intercept = np.sqrt(mz_lower)
        tof_slope = (np.sqrt(mz_upper) - tof_intercept) / n_tof
        mz_values = (tof_intercept + tof_slope * tof_idx) ** 2

        ms1_list.append((rt_min, mz_values, ints))

    logger.info("Processed %d MS1 frames with m/z calibration", len(ms1_list))
    return ms1_list


def create_labeled_asf(df, file_col, raw_dir, output_path, max_charge=4):
    """Create labeled ASF file from DIA-NN matches."""
    files = df[file_col].unique()
    logger.info("Processing %d files", len(files))
    total_written = 0

    if os.path.exists(output_path):
        os.remove(output_path)

    for file_name in files:
        file_df = df[df[file_col] == file_name]

        base_name = os.path.basename(file_name)
        if not base_name.endswith('.d'):
            base_name = os.path.splitext(base_name)[0] + '.d'

        d_path = os.path.join(raw_dir, base_name)
        if not os.path.isdir(d_path):
            d_path = os.path.join(raw_dir, os.path.splitext(base_name)[0] + '.d')
        if not os.path.isdir(d_path):
            logger.warning("Raw file not found: %s (tried %s)", base_name, d_path)
            continue

        logger.info("Processing %s (%d peptides)", base_name, len(file_df))

        # Read MS2 spectra via timsrust
        try:
            import timsrust_pyo3 as timsrust
            spectra = timsrust.read_all_spectra(str(d_path))
            logger.info("  Read %d MS2 spectra", len(spectra))
        except Exception as e:
            logger.error("  Failed to read %s: %s", d_path, e)
            continue

        # Read MS1 frames from TDF
        ms1_list = read_ms1_from_tdf(d_path)
        logger.info("  Read %d MS1 frames", len(ms1_list))

        # Build MS1 lookup by RT (in minutes)
        ms1_by_rt_bin = {}
        for rt_min, mzs, ints in ms1_list:
            rt_bin = int(rt_min * 6)  # 10-second bins in minutes
            ms1_by_rt_bin.setdefault(rt_bin, []).append((rt_min, mzs, ints))

        # Index MS2 spectra by RT and m/z for matching
        ms2_by_rt_mz = {}
        for spec in spectra:
            if spec.precursor.mz is None or float(spec.precursor.mz) == 0.0:
                continue
            rt_min = float(spec.precursor.rt)  # timsrust gives minutes
            prec_mz = float(spec.precursor.mz)
            mz_arr = np.array(spec.mz_values, dtype=np.float64)
            int_arr = np.array(spec.intensities, dtype=np.float64)

            # Apply same normalization as augment.py MS2: sqrt + normalize
            if len(int_arr) > 0:
                sorted_idx = np.argsort(int_arr)[-150:]
                int_arr = int_arr[sorted_idx]
                mz_arr = mz_arr[sorted_idx]
                sorted_mz = np.argsort(mz_arr)
                int_arr = int_arr[sorted_mz]
                mz_arr = mz_arr[sorted_mz]
                int_arr = int_arr ** 0.5
                if np.max(int_arr) > 0:
                    int_arr = int_arr / np.max(int_arr)

            key = (int(prec_mz / 10), int(rt_min * 6))  # 10-sec bins
            ms2_by_rt_mz.setdefault(key, []).append({
                'mz': prec_mz, 'rt': rt_min,
                'fragments': list(zip(mz_arr.tolist(), int_arr.tolist()))
            })

        # Get isolation half-width from spectra (diaPASEF windows are typically 12-25 Da)
        iso_half = 13.0  # default half-width
        if spectra and hasattr(spectra[0], 'isolation_width') and spectra[0].isolation_width:
            iso_half = float(spectra[0].isolation_width) / 2.0
        logger.info("  Isolation half-width: %.1f Da", iso_half)

        # Match DIA-NN peptides to spectra
        # Key insight: timsrust returns one spectrum per DIA window (not per precursor).
        # A peptide at m/z 523.7 falls in the window centered at 512.5 (±12.5 Da).
        # We must search neighboring m/z bins to find the containing window.
        prec_to_spec = {}
        n_matched = 0
        n_no_ms2 = 0
        n_no_ms1 = 0
        n_low_quality = 0

        for _, row in file_df.iterrows():
            peptide = row['Stripped.Sequence']
            prec_mz = row['Precursor.Mz']
            rt_val = row['RT']
            rt_min = rt_val if rt_val < 200 else rt_val / 60.0

            # Search neighboring m/z bins (±3 bins = ±30 Da covers any DIA window)
            rt_bin = int(rt_min * 6)
            mz_bin = int(prec_mz / 10)
            best_spec = None
            best_dist = float('inf')

            for mz_offset in range(-3, 4):
                for rt_offset in range(-1, 2):
                    key = (mz_bin + mz_offset, rt_bin + rt_offset)
                    for cand in ms2_by_rt_mz.get(key, []):
                        # Check if peptide m/z falls within this DIA window
                        if abs(cand['mz'] - prec_mz) > iso_half + 1.0:
                            continue
                        dist = abs(cand['rt'] - rt_min) + abs(cand['mz'] - prec_mz) * 0.1
                        if dist < best_dist:
                            best_dist = dist
                            best_spec = cand

            if best_spec is None:
                n_no_ms2 += 1
                continue

            # Quality filter: ≥10% of theoretical b/y ions present
            # (relaxed from paper's 25% for chimeric timsTOF diaPASEF windows)
            if not check_fragment_coverage(best_spec['fragments'], peptide, 0.10):
                n_low_quality += 1
                continue

            prec_key = (best_spec['mz'], best_spec['rt'], 0)
            if prec_key not in prec_to_spec:
                prec_to_spec[prec_key] = {
                    'scans': [], 'rts': [],
                    'ms1_scans': [], 'ms1_rts': [],
                    'window_width': 25.0,
                    'peptide': peptide
                }

            prec_to_spec[prec_key]['scans'].append(best_spec['fragments'])
            prec_to_spec[prec_key]['rts'].append(0.0)

            # Find MS1 context near this RT
            rt_bin = int(best_spec['rt'] * 6)
            ms1_found = False
            for rt_search in range(rt_bin - 1, rt_bin + 2):
                for ms1_rt, ms1_mz, ms1_int in ms1_by_rt_bin.get(rt_search, []):
                    if abs(ms1_rt - best_spec['rt']) < 0.2:  # within 12 seconds
                        # Apply same normalization as augment.py MS1: sqrt(sqrt) + normalize
                        top_n = 150
                        if len(ms1_int) > top_n:
                            top_idx = np.argsort(ms1_int)[-top_n:]
                            ms1_int_sel = ms1_int[top_idx]
                            ms1_mz_sel = ms1_mz[top_idx]
                        else:
                            ms1_int_sel = ms1_int.copy()
                            ms1_mz_sel = ms1_mz.copy()

                        sorted_mz = np.argsort(ms1_mz_sel)
                        ms1_int_sel = ms1_int_sel[sorted_mz]
                        ms1_mz_sel = ms1_mz_sel[sorted_mz]

                        ms1_int_sel = ms1_int_sel ** 0.5
                        ms1_int_sel = ms1_int_sel ** 0.5
                        if len(ms1_int_sel) > 0 and np.max(ms1_int_sel) > 0:
                            ms1_int_sel = ms1_int_sel / np.max(ms1_int_sel)

                        ms1_arr = list(zip(ms1_mz_sel.tolist(), ms1_int_sel.tolist()))
                        prec_to_spec[prec_key]['ms1_scans'].append(ms1_arr)
                        prec_to_spec[prec_key]['ms1_rts'].append(ms1_rt - best_spec['rt'])
                        ms1_found = True
                        break
                if ms1_found:
                    break

            if not ms1_found:
                # Synthetic MS1: create a simple isotope pattern centered on precursor m/z
                # This is an approximation but gives the model precursor mass context
                synthetic_mz = [prec_mz + i * 1.003355 for i in range(4)]
                synthetic_int = [1.0, 0.6, 0.25, 0.08]
                prec_to_spec[prec_key]['ms1_scans'].append(
                    list(zip(synthetic_mz, synthetic_int))
                )
                prec_to_spec[prec_key]['ms1_rts'].append(0.0)
                n_no_ms1 += 1

            n_matched += 1

        logger.info("  Matched: %d, no MS2: %d, low quality: %d, synthetic MS1: %d",
                    n_matched, n_no_ms2, n_low_quality, n_no_ms1)

        # Write labeled ASF
        n_written = write_labeled_asf(output_path, prec_to_spec, max_charge=max_charge)
        total_written += n_written
        logger.info("  Wrote %d labeled spectra from %s", n_written, base_name)

    logger.info("Total: %d labeled spectra written to %s", total_written, output_path)
    return total_written


def write_labeled_asf(outfile, prec_to_spec, scan_width=1, max_pep_length=30, max_charge=3):
    """Write ASF with real peptide labels (for training)."""
    count = 0
    with open(outfile, 'a') as out:
        for key, value in prec_to_spec.items():
            if 'ms1_scans' not in value or len(value['ms1_scans']) == 0:
                continue

            prec, rt, _ = key
            peptide = value.get('peptide', 'K' * max_pep_length)

            if len(peptide) > max_pep_length:
                continue

            scans = np.array(value['scans'], dtype=object)
            rts = np.array(value['rts'])
            ms1_scans = np.array(value['ms1_scans'], dtype=object)
            ms1_rts = np.array(value['ms1_rts'])
            window_width = value['window_width']

            for charge in range(2, max_charge + 1):
                count += 1
                out.write("BEGIN IONS\n")
                out.write(f"TITLE={count}\n")
                out.write(f"PEPMASS={prec}\n")
                out.write(f"CHARGE={charge}\n")
                out.write(f"SCAN={count}\n")
                out.write(f"RT={rt}\n")
                out.write(f"SEQ={peptide}\n")

                for scan, cur_rt in zip(scans, rts):
                    for mz, intensity in scan:
                        out.write(f"{mz}\t{intensity}\t{cur_rt}\t2\n")

                for scan, cur_rt in zip(ms1_scans, ms1_rts):
                    for mz, intensity in scan:
                        if abs(mz - prec) < window_width + 1:
                            out.write(f"{mz}\t{intensity}\t{cur_rt}\t1\n")

                out.write("END IONS\n")

    return count


def main():
    parser = argparse.ArgumentParser(description="Prepare Cascadia training data from DIA-NN results")
    parser.add_argument("--report", required=True, help="Path to DIA-NN report.parquet")
    parser.add_argument("--raw-dir", required=True, help="Directory containing .d files")
    parser.add_argument("--output", required=True, help="Output ASF file path")
    parser.add_argument("--q-threshold", type=float, default=0.01, help="Q-value threshold (default 0.01)")
    parser.add_argument("--max-charge", type=int, default=4, help="Max charge state (default 4)")

    args = parser.parse_args()

    t0 = time.time()
    df, file_col = load_diann_report(args.report, args.q_threshold)
    n = create_labeled_asf(df, file_col, args.raw_dir, args.output, args.max_charge)
    logger.info("Done in %.1f min. Total labeled spectra: %d", (time.time() - t0) / 60, n)


if __name__ == "__main__":
    main()
