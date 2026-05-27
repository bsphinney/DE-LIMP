"""
Prepare Cascadia training data from DIA-NN searched timsTOF .d files.

Reads DIA-NN report.parquet to get peptide-spectrum matches,
then creates labeled ASF files for Cascadia fine-tuning.

Usage:
    python prepare_training_data.py \
        --report /path/to/report.parquet \
        --raw-dir /path/to/raw_files/ \
        --output-dir /path/to/training_data/ \
        --split train  # or val/test

Requirements:
    - DIA-NN report.parquet with Stripped.Sequence, Precursor.Mz, RT columns
    - Raw .d files (reads via timsrust_pyo3 or mzML)
    - Cascadia environment with augment.py
"""

import argparse
import logging
import os
import sys
import time
import numpy as np

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(name)s %(message)s')
logger = logging.getLogger("prepare_training")


def load_diann_report(report_path, q_threshold=0.01):
    """Load DIA-NN report and extract high-confidence peptide-spectrum matches."""
    import pyarrow.parquet as pq

    logger.info("Loading DIA-NN report: %s", report_path)
    df = pq.read_table(report_path).to_pandas()
    logger.info("Total rows: %d", len(df))

    # Filter to high confidence
    if 'Q.Value' in df.columns:
        df = df[df['Q.Value'] <= q_threshold]
        logger.info("After Q.Value <= %s filter: %d rows", q_threshold, len(df))

    # Required columns
    required = ['Stripped.Sequence', 'Precursor.Mz', 'RT']
    missing = [c for c in required if c not in df.columns]
    if missing:
        logger.error("Missing required columns: %s", missing)
        logger.info("Available columns: %s", list(df.columns)[:20])
        sys.exit(1)

    # Get unique peptide-precursor pairs per file
    if 'File.Name' in df.columns:
        file_col = 'File.Name'
    elif 'Run' in df.columns:
        file_col = 'Run'
    else:
        file_col = df.columns[0]

    logger.info("Unique peptides: %d", df['Stripped.Sequence'].nunique())
    logger.info("Unique files: %d", df[file_col].nunique())

    return df, file_col


def create_labeled_asf(df, file_col, raw_dir, output_path, max_charge=4):
    """Create labeled ASF file from DIA-NN matches.

    The ASF format for training includes peptide labels in the SEQ field
    (instead of the dummy 'K' * max_pep_length used in inference).
    """
    from cascadia.augment import write_asf

    # Group by file
    files = df[file_col].unique()
    logger.info("Processing %d files", len(files))

    total_written = 0

    # Remove existing output
    if os.path.exists(output_path):
        os.remove(output_path)

    for file_name in files:
        file_df = df[df[file_col] == file_name]

        # Find the .d file
        base_name = os.path.basename(file_name)
        if not base_name.endswith('.d'):
            base_name = os.path.splitext(base_name)[0] + '.d'

        d_path = os.path.join(raw_dir, base_name)
        if not os.path.isdir(d_path):
            # Try without extension
            d_path = os.path.join(raw_dir, os.path.splitext(base_name)[0] + '.d')
        if not os.path.isdir(d_path):
            logger.warning("Raw file not found: %s (tried %s)", base_name, d_path)
            continue

        logger.info("Processing %s (%d peptides)", base_name, len(file_df))

        # Build precursor-to-spectrum map with labels
        prec_to_spec = {}
        try:
            import timsrust_pyo3 as timsrust
            spectra = timsrust.read_all_spectra(str(d_path))
        except Exception as e:
            logger.error("Failed to read %s: %s", d_path, e)
            continue

        # Index spectra by RT and m/z for matching
        ms1_by_rt = {}
        ms2_by_rt_mz = {}

        for spec in spectra:
            rt_sec = spec.precursor.rt * 60.0
            mz_arr = list(zip(
                [float(m) for m in spec.mz_values],
                [float(i) for i in spec.intensities]
            ))

            if spec.precursor.mz is None or spec.precursor.mz == 0.0:
                rt_key = int(rt_sec / 10)
                ms1_by_rt.setdefault(rt_key, []).append((rt_sec, mz_arr))
            else:
                prec_mz = float(spec.precursor.mz)
                key = (int(prec_mz / 10), int(rt_sec / 10))
                ms2_by_rt_mz.setdefault(key, []).append({
                    'mz': prec_mz, 'rt': rt_sec,
                    'fragments': mz_arr[:150]  # top 150 fragments
                })

        # Match DIA-NN peptides to spectra
        for _, row in file_df.iterrows():
            peptide = row['Stripped.Sequence']
            prec_mz = row['Precursor.Mz']
            rt_sec = row['RT'] * 60.0 if row['RT'] < 200 else row['RT']

            # Find closest MS2 spectrum
            key = (int(prec_mz / 10), int(rt_sec / 10))
            candidates = ms2_by_rt_mz.get(key, [])

            best_spec = None
            best_dist = float('inf')
            for cand in candidates:
                dist = abs(cand['rt'] - rt_sec) + abs(cand['mz'] - prec_mz) * 10
                if dist < best_dist:
                    best_dist = dist
                    best_spec = cand

            if best_spec is None:
                continue

            prec_key = (best_spec['mz'], best_spec['rt'], 0)
            if prec_key not in prec_to_spec:
                prec_to_spec[prec_key] = {
                    'scans': [], 'rts': [],
                    'ms1_scans': [], 'ms1_rts': [],
                    'window_width': 25.0,
                    'peptide': peptide  # LABEL
                }

            prec_to_spec[prec_key]['scans'].append(best_spec['fragments'])
            prec_to_spec[prec_key]['rts'].append(0.0)  # delta RT = 0 for center

            # Add MS1 context
            rt_key = int(best_spec['rt'] / 10)
            for ms1_rt, ms1_frags in ms1_by_rt.get(rt_key, [])[:1]:
                prec_to_spec[prec_key]['ms1_scans'].append(ms1_frags[:150])
                prec_to_spec[prec_key]['ms1_rts'].append(ms1_rt - best_spec['rt'])

        # Write labeled ASF (modified write_asf with real peptide labels)
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

            # Truncate or pad peptide to max_pep_length
            if len(peptide) > max_pep_length:
                continue  # Skip peptides longer than model can handle

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
                out.write(f"SEQ={peptide}\n")  # REAL LABEL

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
