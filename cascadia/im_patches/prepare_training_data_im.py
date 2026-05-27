"""
Prepare Cascadia IM-enhanced training data from DIA-NN searched timsTOF .d files.

Key differences from prepare_training_data.py:
1. Uses SpectrumReader.new_with_span_step() for IM-sliced spectra
2. Matches peptides using RT + m/z + IM triple match (DIA-NN report has IM column)
3. Writes 5-column ASF with real 1/K0 values per peak
4. MS1 peaks get IM=0.0, MS2 peaks get the IM slice midpoint

The DIA-NN report.parquet contains an 'IM' column with 1/K0 per peptide.
We use this to select the correct IM slice when matching peptides to spectra.
"""

import argparse
import logging
import os
import sys
import time
import sqlite3
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(message)s')
logger = logging.getLogger("prepare_training_im")


def compute_theoretical_ions(peptide, max_charge=2):
    """Compute theoretical b and y ion m/z values for a peptide."""
    ions = set()
    try:
        aa_masses = {
            'G': 57.02146, 'A': 71.03711, 'V': 99.06841, 'L': 113.08406,
            'I': 113.08406, 'P': 97.05276, 'F': 147.06841, 'W': 186.07931,
            'M': 131.04049, 'S': 87.03203, 'T': 101.04768, 'C': 160.03065,
            'Y': 163.06333, 'H': 137.05891, 'D': 115.02694, 'E': 129.04259,
            'N': 114.04293, 'Q': 128.05858, 'K': 128.09496, 'R': 156.10111,
        }
        proton = 1.007276

        masses = [aa_masses.get(aa, 0) for aa in peptide]
        total = sum(masses) + 18.01056

        cum = 0
        for i in range(len(masses) - 1):
            cum += masses[i]
            for z in range(1, max_charge + 1):
                ions.add(round((cum + z * proton) / z, 1))

        cum = 18.01056
        for i in range(len(masses) - 1, 0, -1):
            cum += masses[i]
            for z in range(1, max_charge + 1):
                ions.add(round((cum + z * proton) / z, 1))
    except Exception:
        pass
    return ions


def check_fragment_coverage(fragments, peptide, min_coverage=0.10):
    """Check if spectrum contains sufficient theoretical b/y ions."""
    theoretical = compute_theoretical_ions(peptide, max_charge=2)
    if not theoretical:
        return True

    frag_mzs = set()
    for mz, _, _ in fragments:  # 3-tuple: (mz, intensity, im)
        frag_mzs.add(round(mz * 2) / 2)

    n_matched = 0
    for ion_mz in theoretical:
        ion_bin = round(ion_mz * 2) / 2
        if ion_bin in frag_mzs or (ion_bin - 0.5) in frag_mzs or (ion_bin + 0.5) in frag_mzs:
            n_matched += 1

    coverage = n_matched / len(theoretical) if theoretical else 0
    return coverage >= min_coverage


def load_diann_report(report_path, q_threshold=0.01):
    """Load DIA-NN report with IM column."""
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

    has_im = 'IM' in df.columns
    if has_im:
        logger.info("IM column found in report (1/K0 values)")
        logger.info("  IM range: %.3f - %.3f", df['IM'].min(), df['IM'].max())
    else:
        logger.warning("No IM column in report -- will skip IM matching")

    if 'File.Name' in df.columns:
        file_col = 'File.Name'
    elif 'Run' in df.columns:
        file_col = 'Run'
    else:
        file_col = df.columns[0]

    logger.info("Unique peptides: %d", df['Stripped.Sequence'].nunique())
    logger.info("Unique files: %d", df[file_col].nunique())

    return df, file_col, has_im


def _read_mz_calibration(d_path):
    """Read m/z calibration from analysis.tdf."""
    tdf_path = os.path.join(d_path, "analysis.tdf")
    if not os.path.exists(tdf_path):
        return 100.0, 1700.0, 631025

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
    """Convert tof indices to m/z using Sage formula."""
    tof = np.array(tof_indices, dtype=np.float64)
    tof_intercept = np.sqrt(mz_lower)
    tof_slope = (np.sqrt(mz_upper) - tof_intercept) / dns
    return (tof_intercept + tof_slope * tof) ** 2


def _read_ms1_frames(d_path, mz_lower, mz_upper, dns, top_n=150):
    """Read MS1 frames with tof-to-mz calibration. Returns (rt_min, mz, int) list."""
    import timsrust_pyo3 as timsrust

    reader = timsrust.FrameReader(str(d_path))
    ms1_frames = reader.read_ms1_frames()

    ms1_list = []
    for frame in ms1_frames:
        tof = np.array(frame.tof_indices, dtype=np.float64)
        ints = np.array(frame.intensities, dtype=np.float64)
        if len(tof) == 0:
            continue

        mzs = _tof_to_mz(tof, mz_lower, mz_upper, dns)

        sorted_idx = np.argsort(ints)[-top_n:]
        ints = ints[sorted_idx]
        mzs = mzs[sorted_idx]
        sorted_mz = np.argsort(mzs)
        ints = ints[sorted_mz]
        mzs = mzs[sorted_mz]

        ints = ints ** 0.5
        ints = ints ** 0.5
        if len(ints) > 0 and np.max(ints) > 0:
            ints = ints / np.max(ints)

        ms1_list.append((frame.rt, mzs, ints))

    return ms1_list


def create_labeled_asf_im(df, file_col, raw_dir, output_path, has_im=True,
                           max_charge=4, im_span=0.05, im_step=0.02,
                           im_tolerance=0.03):
    """Create labeled 5-column ASF file from DIA-NN matches using IM-sliced spectra.

    Args:
        df: DIA-NN report DataFrame
        file_col: column name for file paths
        raw_dir: directory containing .d files
        output_path: output ASF file path
        has_im: whether DIA-NN report has IM column
        max_charge: max charge state for augmentation
        im_span: IM slice width for new_with_span_step
        im_step: IM step size for new_with_span_step
        im_tolerance: tolerance for IM matching (1/K0 units)
    """
    import timsrust_pyo3 as timsrust

    files = df[file_col].unique()
    logger.info("Processing %d files with IM-sliced spectra", len(files))
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

        # Read IM-sliced MS2 spectra
        try:
            reader = timsrust.SpectrumReader.new_with_span_step(str(d_path), im_span, im_step)
            n_spectra = len(reader)
            logger.info("  Read %d IM-sliced MS2 spectra", n_spectra)
        except Exception as e:
            logger.error("  Failed to read %s with IM slicing: %s", d_path, e)
            continue

        # Read MS1 frames
        mz_lower, mz_upper, dns = _read_mz_calibration(d_path)
        ms1_list = _read_ms1_frames(d_path, mz_lower, mz_upper, dns)
        logger.info("  Read %d MS1 frames", len(ms1_list))

        # Build MS1 lookup by RT
        ms1_by_rt_bin = {}
        for rt_min, mzs, ints in ms1_list:
            rt_bin = int(rt_min * 6)
            ms1_by_rt_bin.setdefault(rt_bin, []).append((rt_min, mzs, ints))

        # Index IM-sliced MS2 spectra by RT, m/z, and IM
        ms2_by_rt_mz = {}
        for spec in reader:
            if spec.precursor.mz is None or float(spec.precursor.mz) == 0.0:
                continue
            rt_min = float(spec.precursor.rt)
            prec_mz = float(spec.precursor.mz)
            im_value = float(spec.precursor.im) if hasattr(spec.precursor, 'im') and spec.precursor.im is not None else 0.0
            mz_arr = np.array(spec.mz_values, dtype=np.float64)
            int_arr = np.array(spec.intensities, dtype=np.float64)

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

            key = (int(prec_mz / 10), int(rt_min * 6))
            # Store as 3-tuple peaks: (mz, intensity, im)
            fragments = [(m, i, im_value) for m, i in zip(mz_arr.tolist(), int_arr.tolist())]
            ms2_by_rt_mz.setdefault(key, []).append({
                'mz': prec_mz, 'rt': rt_min, 'im': im_value,
                'fragments': fragments
            })

        # Get isolation half-width
        iso_half = 13.0

        # Match DIA-NN peptides to IM-sliced spectra
        prec_to_spec = {}
        n_matched = 0
        n_no_ms2 = 0
        n_no_ms1 = 0
        n_low_quality = 0
        n_im_matched = 0

        for _, row in file_df.iterrows():
            peptide = row['Stripped.Sequence']
            prec_mz = row['Precursor.Mz']
            rt_val = row['RT']
            rt_min = rt_val if rt_val < 200 else rt_val / 60.0
            pep_im = float(row['IM']) if (has_im and 'IM' in row.index) else None

            rt_bin = int(rt_min * 6)
            mz_bin = int(prec_mz / 10)
            best_spec = None
            best_dist = float('inf')

            for mz_offset in range(-3, 4):
                for rt_offset in range(-1, 2):
                    key = (mz_bin + mz_offset, rt_bin + rt_offset)
                    for cand in ms2_by_rt_mz.get(key, []):
                        if abs(cand['mz'] - prec_mz) > iso_half + 1.0:
                            continue

                        # Compute distance: RT + m/z + IM (if available)
                        dist = abs(cand['rt'] - rt_min) + abs(cand['mz'] - prec_mz) * 0.1
                        if pep_im is not None:
                            im_dist = abs(cand['im'] - pep_im)
                            if im_dist > im_tolerance * 3:
                                # Skip spectra far from the expected IM
                                continue
                            dist += im_dist * 5.0  # Weight IM matching
                            if im_dist <= im_tolerance:
                                n_im_matched += 1

                        if dist < best_dist:
                            best_dist = dist
                            best_spec = cand

            if best_spec is None:
                n_no_ms2 += 1
                continue

            # Quality filter
            if not check_fragment_coverage(best_spec['fragments'], peptide, 0.10):
                n_low_quality += 1
                continue

            prec_key = (best_spec['mz'], best_spec['rt'], 0, best_spec['im'])
            if prec_key not in prec_to_spec:
                prec_to_spec[prec_key] = {
                    'scans': [], 'rts': [],
                    'ms1_scans': [], 'ms1_rts': [],
                    'window_width': 25.0,
                    'peptide': peptide,
                    'im': best_spec['im']
                }

            prec_to_spec[prec_key]['scans'].append(best_spec['fragments'])
            prec_to_spec[prec_key]['rts'].append(0.0)

            # Find MS1 context near this RT
            rt_bin_ms1 = int(best_spec['rt'] * 6)
            ms1_found = False
            for rt_search in range(rt_bin_ms1 - 1, rt_bin_ms1 + 2):
                for ms1_rt, ms1_mz, ms1_int in ms1_by_rt_bin.get(rt_search, []):
                    if abs(ms1_rt - best_spec['rt']) < 0.2:
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

                        # MS1 peaks: IM = 0.0
                        ms1_triples = [(m, i, 0.0) for m, i in
                                       zip(ms1_mz_sel.tolist(), ms1_int_sel.tolist())]
                        prec_to_spec[prec_key]['ms1_scans'].append(ms1_triples)
                        prec_to_spec[prec_key]['ms1_rts'].append(ms1_rt - best_spec['rt'])
                        ms1_found = True
                        break
                if ms1_found:
                    break

            if not ms1_found:
                synthetic_mz = [prec_mz + i * 1.003355 for i in range(4)]
                synthetic_int = [1.0, 0.6, 0.25, 0.08]
                prec_to_spec[prec_key]['ms1_scans'].append(
                    [(m, i, 0.0) for m, i in zip(synthetic_mz, synthetic_int)]
                )
                prec_to_spec[prec_key]['ms1_rts'].append(0.0)
                n_no_ms1 += 1

            n_matched += 1

        logger.info("  Matched: %d, no MS2: %d, low quality: %d, synthetic MS1: %d, IM-matched: %d",
                    n_matched, n_no_ms2, n_low_quality, n_no_ms1, n_im_matched)

        # Write labeled 5-column ASF
        n_written = write_labeled_asf_im(output_path, prec_to_spec, max_charge=max_charge)
        total_written += n_written
        logger.info("  Wrote %d labeled IM-enhanced spectra from %s", n_written, base_name)

    logger.info("Total: %d labeled IM-enhanced spectra written to %s", total_written, output_path)
    return total_written


def write_labeled_asf_im(outfile, prec_to_spec, scan_width=1, max_pep_length=30, max_charge=3):
    """Write 5-column ASF with real peptide labels and IM values (for training)."""
    count = 0
    with open(outfile, 'a') as out:
        for key, value in prec_to_spec.items():
            if 'ms1_scans' not in value or len(value['ms1_scans']) == 0:
                continue

            prec, rt, _, im_val = key
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
                out.write(f"IM={im_val:.4f}\n")
                out.write(f"SEQ={peptide}\n")

                # MS2 peaks: 5 columns (mz, intensity, rt_offset, ms_level, im)
                for scan, cur_rt in zip(scans, rts):
                    for mz, intensity, peak_im in scan:
                        out.write(f"{mz}\t{intensity}\t{cur_rt}\t2\t{peak_im:.4f}\n")

                # MS1 peaks: 5 columns with IM = 0.0
                for scan, cur_rt in zip(ms1_scans, ms1_rts):
                    for mz, intensity, peak_im in scan:
                        if abs(mz - prec) < window_width + 1:
                            out.write(f"{mz}\t{intensity}\t{cur_rt}\t1\t{peak_im:.4f}\n")

                out.write("END IONS\n")

    return count


def main():
    parser = argparse.ArgumentParser(
        description="Prepare IM-enhanced Cascadia training data from DIA-NN results"
    )
    parser.add_argument("--report", required=True, help="Path to DIA-NN report.parquet")
    parser.add_argument("--raw-dir", required=True, help="Directory containing .d files")
    parser.add_argument("--output", required=True, help="Output ASF file path")
    parser.add_argument("--q-threshold", type=float, default=0.01, help="Q-value threshold")
    parser.add_argument("--max-charge", type=int, default=4, help="Max charge state")
    parser.add_argument("--im-span", type=float, default=0.05, help="IM slice width (1/K0)")
    parser.add_argument("--im-step", type=float, default=0.02, help="IM step size (1/K0)")
    parser.add_argument("--im-tolerance", type=float, default=0.03,
                        help="IM matching tolerance (1/K0 units)")

    args = parser.parse_args()

    t0 = time.time()
    df, file_col, has_im = load_diann_report(args.report, args.q_threshold)
    n = create_labeled_asf_im(
        df, file_col, args.raw_dir, args.output,
        has_im=has_im, max_charge=args.max_charge,
        im_span=args.im_span, im_step=args.im_step,
        im_tolerance=args.im_tolerance
    )
    logger.info("Done in %.1f min. Total labeled IM-enhanced spectra: %d",
                (time.time() - t0) / 60, n)


if __name__ == "__main__":
    main()
