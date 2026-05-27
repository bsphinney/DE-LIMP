#!/usr/bin/env python3
"""
Prepare Cascadia training data using mobility-filtered spectrum extraction.

Same purpose as prepare_training_data.py but uses mobility-filtered extraction
instead of collapsed spectra. Key differences:

1. Extracts per-window mobilogram peaks (Mode B from CASCADIA_MOBILITY_FILTER_ADDENDUM)
2. Matches DIA-NN peptides to mobility peaks using RT + m/z + IM triple match
3. Writes 5-column ASF with real 1/K0 per spectrum
4. Quality filter: 10% b/y ion coverage (same as prepare_training_data.py)

Usage:
    python prepare_training_data_mobfilter.py \
        --report /path/to/report.parquet \
        --raw-dir /path/to/raw_files/ \
        --output /path/to/training.asf \
        --q-threshold 0.01
"""

import argparse
import logging
import os
import sys
import time
import sqlite3
import numpy as np
from pathlib import Path
from scipy.signal import find_peaks

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(message)s")
logger = logging.getLogger("prepare_training_mobfilter")


# ============================================================================
# Fragment quality filter (same as prepare_training_data.py)
# ============================================================================

def compute_theoretical_ions(peptide, max_charge=2):
    """Compute theoretical b and y ion m/z values for a peptide."""
    ions = set()
    try:
        aa_masses = {
            "G": 57.02146, "A": 71.03711, "V": 99.06841, "L": 113.08406,
            "I": 113.08406, "P": 97.05276, "F": 147.06841, "W": 186.07931,
            "M": 131.04049, "S": 87.03203, "T": 101.04768, "C": 160.03065,
            "Y": 163.06333, "H": 137.05891, "D": 115.02694, "E": 129.04259,
            "N": 114.04293, "Q": 128.05858, "K": 128.09496, "R": 156.10111,
        }
        proton = 1.007276
        masses = [aa_masses.get(aa, 0) for aa in peptide]

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


def check_fragment_coverage(mz_array, peptide, min_coverage=0.10):
    """Check if spectrum contains sufficient theoretical b/y ions."""
    theoretical = compute_theoretical_ions(peptide, max_charge=2)
    if not theoretical:
        return True

    frag_mzs = set()
    for mz in mz_array:
        frag_mzs.add(round(float(mz) * 2) / 2)

    n_matched = 0
    for ion_mz in theoretical:
        ion_bin = round(ion_mz * 2) / 2
        if ion_bin in frag_mzs or (ion_bin - 0.5) in frag_mzs or (ion_bin + 0.5) in frag_mzs:
            n_matched += 1

    coverage = n_matched / len(theoretical) if theoretical else 0
    return coverage >= min_coverage


# ============================================================================
# Import from bruker_mobility_filter.py
# ============================================================================

# These are imported at runtime to allow standalone testing
from bruker_mobility_filter import (
    MobilityFilteredExtractor,
    _tof_to_mz,
    _read_mz_calibration,
)


def load_diann_report(report_path, q_threshold=0.01):
    """Load DIA-NN report and extract high-confidence peptide-spectrum matches."""
    import pyarrow.parquet as pq

    logger.info("Loading DIA-NN report: %s", report_path)
    df = pq.read_table(report_path).to_pandas()
    logger.info("Total rows: %d", len(df))

    if "Q.Value" in df.columns:
        df = df[df["Q.Value"] <= q_threshold]
        logger.info("After Q.Value <= %s filter: %d rows", q_threshold, len(df))

    required = ["Stripped.Sequence", "Precursor.Mz", "RT"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        logger.error("Missing required columns: %s", missing)
        sys.exit(1)

    has_im = "IM" in df.columns
    logger.info("IM column available: %s", has_im)

    if "File.Name" in df.columns:
        file_col = "File.Name"
    elif "Run" in df.columns:
        file_col = "Run"
    else:
        file_col = df.columns[0]

    logger.info("Unique peptides: %d", df["Stripped.Sequence"].nunique())
    logger.info("Unique files: %d", df[file_col].nunique())

    return df, file_col


def match_peptides_to_spectra(file_df, spectra, im_tol=0.05, rt_tol_min=0.3):
    """Match DIA-NN peptide IDs to mobility-filtered spectra using RT + m/z + IM.

    Args:
        file_df: DIA-NN report rows for this file
        spectra: List of mobility-filtered spectrum dicts from MobilityFilteredExtractor
        im_tol: 1/K0 tolerance for matching (default 0.05 Vs/cm2)
        rt_tol_min: RT tolerance in minutes (default 0.3 min = 18 sec)

    Returns:
        List of matched (spectrum_dict, peptide, charge) tuples
    """
    has_im = "IM" in file_df.columns

    # Index spectra by (mz_bin, rt_bin) for fast lookup
    spec_index = {}
    for spec in spectra:
        mz_bin = int(spec["precursor_mz"] / 10)
        rt_bin = int(spec["rt"] * 6)  # 10-sec bins
        for mb in range(mz_bin - 3, mz_bin + 4):
            for rb in range(rt_bin - 2, rt_bin + 3):
                spec_index.setdefault((mb, rb), []).append(spec)

    matched = []
    n_no_match = 0
    n_low_quality = 0

    for _, row in file_df.iterrows():
        peptide = row["Stripped.Sequence"]
        prec_mz = row["Precursor.Mz"]
        rt_val = row["RT"]
        rt_min = rt_val if rt_val < 200 else rt_val / 60.0
        charge = int(row["Precursor.Charge"]) if "Precursor.Charge" in row.index else 0
        pep_im = float(row["IM"]) if has_im and not np.isnan(row.get("IM", np.nan)) else None

        mz_bin = int(prec_mz / 10)
        rt_bin = int(rt_min * 6)
        key = (mz_bin, rt_bin)

        # Find best matching spectrum
        best_spec = None
        best_score = float("inf")

        for cand in spec_index.get(key, []):
            # Check m/z: peptide must fall within the isolation window
            iso_half = cand["isolation_width"] / 2.0
            if abs(cand["precursor_mz"] - prec_mz) > iso_half + 1.0:
                continue

            # Check RT
            dt = abs(cand["rt"] - rt_min)
            if dt > rt_tol_min:
                continue

            # Check IM if available
            if pep_im is not None and cand["mobility_filtered"]:
                dim = abs(cand["im"] - pep_im)
                if dim > im_tol:
                    continue
                score = dt + dim * 10  # IM matching weighted heavily
            else:
                score = dt + abs(cand["precursor_mz"] - prec_mz) * 0.01

            if score < best_score:
                best_score = score
                best_spec = cand

        if best_spec is None:
            n_no_match += 1
            continue

        # Quality filter: 10% b/y ion coverage
        if not check_fragment_coverage(best_spec["mz"], peptide, 0.10):
            n_low_quality += 1
            continue

        matched.append((best_spec, peptide, charge))

    return matched, n_no_match, n_low_quality


def write_labeled_asf_5col(outfile, matched_list, ms1_list, max_pep_length=30, max_charge=4):
    """Write labeled 5-column ASF with peptide annotations.

    Args:
        outfile: Output file path (appended to)
        matched_list: List of (spectrum_dict, peptide, charge) tuples
        ms1_list: List of (rt_min, mz_array, int_array) from MS1 extraction
        max_pep_length: Maximum peptide length to include
        max_charge: Max charge to enumerate per peptide

    Returns:
        Number of spectra written
    """
    # Build MS1 lookup by RT
    ms1_by_rt_bin = {}
    for rt_min, mzs, ints in ms1_list:
        rt_bin = int(rt_min * 6)
        ms1_by_rt_bin.setdefault(rt_bin, []).append((rt_min, mzs, ints))

    count = 0
    top_n_ms2 = 150

    with open(outfile, "a") as out:
        for spec, peptide, orig_charge in matched_list:
            if len(peptide) > max_pep_length:
                continue

            mz = spec["mz"].copy()
            intensity = spec["intensity"].copy().astype(np.float64)
            rt = spec["rt"]
            prec_mz = spec["precursor_mz"]
            im = spec["im"]

            # MS2 normalization: top_n, sqrt, normalize
            if len(intensity) > top_n_ms2:
                top_idx = np.argsort(intensity)[-top_n_ms2:]
                intensity = intensity[top_idx]
                mz = mz[top_idx]

            sorted_mz = np.argsort(mz)
            intensity = intensity[sorted_mz]
            mz = mz[sorted_mz]

            intensity = intensity ** 0.5
            if len(intensity) > 0 and np.max(intensity) > 0:
                intensity = intensity / np.max(intensity)

            # Find closest MS1 frame
            rt_bin = int(rt * 6)
            best_ms1 = None
            best_dt = float("inf")
            for search_bin in range(rt_bin - 1, rt_bin + 2):
                for ms1_rt, ms1_mz, ms1_int in ms1_by_rt_bin.get(search_bin, []):
                    dt = abs(ms1_rt - rt)
                    if dt < best_dt:
                        best_dt = dt
                        best_ms1 = (ms1_rt, ms1_mz, ms1_int)

            if best_ms1 is None:
                continue

            # Write for each charge state
            charges = [orig_charge] if orig_charge > 0 else range(2, max_charge + 1)
            for charge in charges:
                count += 1
                out.write("BEGIN IONS\n")
                out.write(f"TITLE={count}\n")
                out.write(f"PEPMASS={prec_mz}\n")
                out.write(f"CHARGE={charge}\n")
                out.write(f"SCAN={count}\n")
                out.write(f"RT={rt}\n")
                out.write(f"SEQ={peptide}\n")

                # MS2 peaks with real 1/K0
                for m, i in zip(mz, intensity):
                    out.write(f"{m:.4f}\t{i:.6f}\t0.0\t2\t{im:.4f}\n")

                # MS1 peaks with IM=0.0
                ms1_rt, ms1_mz, ms1_int = best_ms1
                iso_half = spec.get("isolation_width", 25.0) / 2.0 + 1.0
                rt_offset = ms1_rt - rt
                for m, i in zip(ms1_mz, ms1_int):
                    if abs(m - prec_mz) < iso_half + 50:
                        out.write(f"{m:.4f}\t{i:.6f}\t{rt_offset:.4f}\t1\t0.0000\n")

                out.write("END IONS\n")

    return count


def process_file(d_path, file_df, output_path, extractor_kwargs, max_charge=4):
    """Process a single .d file: extract spectra, match peptides, write ASF.

    Returns (n_written, n_matched, n_no_match, n_low_quality, n_spectra)
    """
    logger.info("Processing %s (%d peptides)", os.path.basename(d_path), len(file_df))

    # Extract mobility-filtered spectra
    try:
        extractor = MobilityFilteredExtractor(d_path, **extractor_kwargs)
        spectra = extractor.extract_filtered_spectra()
        logger.info("  Extracted %d spectra", len(spectra))
    except Exception as e:
        logger.error("  Failed to extract spectra from %s: %s", d_path, e)
        return 0, 0, 0, 0, 0

    if not spectra:
        logger.warning("  No spectra extracted")
        return 0, 0, 0, 0, 0

    # Match peptides to spectra
    matched, n_no_match, n_low_quality = match_peptides_to_spectra(file_df, spectra)
    logger.info(
        "  Matched: %d, no match: %d, low quality: %d",
        len(matched), n_no_match, n_low_quality,
    )

    if not matched:
        return 0, 0, n_no_match, n_low_quality, len(spectra)

    # Read MS1 frames for ASF context
    ms1_list = extractor.extract_ms1_frames()
    logger.info("  Read %d MS1 frames", len(ms1_list))

    # Write labeled ASF
    n_written = write_labeled_asf_5col(
        output_path, matched, ms1_list, max_charge=max_charge
    )
    logger.info("  Wrote %d labeled spectra", n_written)

    return n_written, len(matched), n_no_match, n_low_quality, len(spectra)


def main():
    parser = argparse.ArgumentParser(
        description="Prepare Cascadia training data with mobility-filtered extraction"
    )
    parser.add_argument("--report", required=True, help="Path to DIA-NN report.parquet")
    parser.add_argument("--raw-dir", required=True, help="Directory containing .d files")
    parser.add_argument("--output", required=True, help="Output ASF file path")
    parser.add_argument("--q-threshold", type=float, default=0.01, help="Q-value threshold")
    parser.add_argument("--max-charge", type=int, default=4, help="Max charge state")
    parser.add_argument("--mobility-half-width", type=int, default=3, help="Scans on each side of mobility peak")
    parser.add_argument("--mobilogram-prominence", type=float, default=100.0, help="Min prominence for peak detection")
    parser.add_argument("--min-fragments", type=int, default=5, help="Min fragments per spectrum")
    parser.add_argument("--im-tol", type=float, default=0.05, help="1/K0 tolerance for peptide matching")
    parser.add_argument("--single-file", type=str, default=None, help="Process only this .d file (for array jobs)")

    args = parser.parse_args()

    extractor_kwargs = {
        "mobility_half_width": args.mobility_half_width,
        "mobilogram_prominence": args.mobilogram_prominence,
        "min_fragments": args.min_fragments,
    }

    t0 = time.time()

    # Load DIA-NN report
    df, file_col = load_diann_report(args.report, args.q_threshold)

    # Clear output file
    if os.path.exists(args.output):
        os.remove(args.output)

    # Determine which files to process
    all_files = df[file_col].unique()

    if args.single_file:
        # Array job mode: filter report to only rows matching this file
        single_base = os.path.basename(args.single_file)
        if single_base.endswith(".d"):
            single_base = single_base[:-2]
        matched_files = [f for f in all_files if single_base in os.path.basename(f)]
        if not matched_files:
            logger.error("No report rows match file: %s", args.single_file)
            sys.exit(1)
        files = matched_files
        logger.info("Single-file mode: processing %d matching runs for %s", len(files), single_base)
    else:
        files = all_files

    logger.info("Processing %d files", len(files))

    total_written = 0
    total_matched = 0
    total_no_match = 0
    total_low_quality = 0
    total_spectra = 0

    for file_name in files:
        file_df = df[df[file_col] == file_name]

        if args.single_file:
            # Use the provided path directly
            d_path = args.single_file
        else:
            # Resolve .d path
            base_name = os.path.basename(file_name)
            if not base_name.endswith(".d"):
                base_name = os.path.splitext(base_name)[0] + ".d"

            d_path = os.path.join(args.raw_dir, base_name)
            if not os.path.isdir(d_path):
                d_path2 = os.path.join(args.raw_dir, os.path.splitext(base_name)[0] + ".d")
                if os.path.isdir(d_path2):
                    d_path = d_path2
                else:
                    logger.warning("Raw file not found: %s", base_name)
                    continue

        n_wr, n_ma, n_nm, n_lq, n_sp = process_file(
            d_path, file_df, args.output, extractor_kwargs, args.max_charge
        )
        total_written += n_wr
        total_matched += n_ma
        total_no_match += n_nm
        total_low_quality += n_lq
        total_spectra += n_sp

    elapsed = time.time() - t0
    logger.info(
        "\nDone in %.1f min.\n"
        "  Total spectra extracted:     %d\n"
        "  Total peptides matched:      %d\n"
        "  Total unmatched:             %d\n"
        "  Total low quality:           %d\n"
        "  Total labeled spectra written: %d\n"
        "  Output: %s",
        elapsed / 60, total_spectra, total_matched,
        total_no_match, total_low_quality, total_written, args.output,
    )


if __name__ == "__main__":
    main()
