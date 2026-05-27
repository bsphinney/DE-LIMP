#!/usr/bin/env python3
"""
bruker_mobility_filter.py — Mode B: Window-based mobility-filtered spectrum extraction.

Extracts mobility-resolved pseudo-DDA spectra from timsTOF diaPASEF .d files by:
1. Reading raw frame data via FrameReader (preserving per-scan 1/K0 dimension)
2. Building per-window mobilograms (intensity vs scan number)
3. Detecting peaks in mobilograms (scipy.signal.find_peaks)
4. Extracting fragment spectra from narrow mobility slices around each peak
5. Writing 5-column ASF: (m/z, intensity, RT_offset, MS_level, 1/K0)

This produces dramatically cleaner input for de novo sequencing compared to
mobility-summed extraction (the default timsrust read_all_spectra approach).

Key parameters:
    mobility_half_width: 3 scans (~0.004 Vs/cm2 per scan, 7-scan window captures core peak)
    mobilogram_prominence: 100 intensity units (filters noise peaks)
    min_fragments: 5 (minimum peaks per spectrum for Cascadia)
    mz_merge_tol: 0.01 Da (centroid merging within narrow scan slice)
    max_spectra_per_window: 5 (safety cap per isolation window per frame)

Usage:
    python bruker_mobility_filter.py --d-folder /path/to/sample.d --output /path/to/output.asf
    python bruker_mobility_filter.py --d-folder /path/to/sample.d --output /path/to/output.asf --stats-only
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
logger = logging.getLogger("mobility_filter")


def _read_mz_calibration(d_folder):
    """Read m/z calibration parameters from analysis.tdf SQLite."""
    tdf_path = os.path.join(str(d_folder), "analysis.tdf")
    if not os.path.exists(tdf_path):
        return 100.0, 1700.0, 631025

    try:
        conn = sqlite3.connect(f"file:{tdf_path}?mode=ro", uri=True)
        meta = dict(conn.execute("SELECT Key, Value FROM GlobalMetadata").fetchall())
        conn.close()
        mz_lower = float(meta.get("MzAcqRangeLower", 100.0))
        mz_upper = float(meta.get("MzAcqRangeUpper", 1700.0))
        dns = int(meta.get("DigitizerNumSamples", 631025))
        return mz_lower, mz_upper, dns
    except Exception as e:
        logger.warning("Failed to read TDF calibration: %s", e)
        return 100.0, 1700.0, 631025


def _tof_to_mz(tof_indices, mz_lower, mz_upper, dns):
    """Convert tof indices to m/z using Sage formula (timsrust Tof2MzConverter)."""
    tof = np.asarray(tof_indices, dtype=np.float64)
    tof_intercept = np.sqrt(mz_lower)
    tof_slope = (np.sqrt(mz_upper) - tof_intercept) / dns
    return (tof_intercept + tof_slope * tof) ** 2


def _calibrate_scan_to_im(d_folder):
    """Derive scan-to-1/K0 linear calibration from timsrust Spectrum IM values.

    timsrust correctly calibrates 1/K0 on Spectrum.precursor.im using the
    Bruker SDK calibration. We cross-reference these IM values with the
    scan midpoints of their isolation windows (from Frame quadrupole_settings)
    to derive a linear mapping: im = intercept + slope * scan_number.

    Returns (intercept, slope) such that im = intercept + slope * scan.
    """
    import timsrust_pyo3 as timsrust

    # Read a batch of collapsed spectra to get calibrated IM values
    spectra = timsrust.read_all_spectra(str(d_folder))
    if len(spectra) == 0:
        logger.warning("No spectra found for IM calibration, using defaults")
        return 1.3, -0.00074

    # Read frames to get scan ranges per window
    reader = timsrust.FrameReader(str(d_folder))
    dia_frames = reader.read_dia_frames()

    if len(dia_frames) == 0:
        logger.warning("No DIA frames for IM calibration, using defaults")
        return 1.3, -0.00074

    # Build spectrum lookup by (isolation_mz, rt) for matching
    spec_lookup = {}
    for s in spectra[:5000]:  # sample first 5000 spectra
        if s.precursor.im is None or s.precursor.mz is None:
            continue
        key = (round(s.precursor.mz, 1), round(s.precursor.rt, 3))
        spec_lookup[key] = s.precursor.im

    # Collect (scan_midpoint, im) pairs
    pairs = []
    for f in dia_frames[:500]:  # sample first 500 frames
        qs = f.quadrupole_settings
        for i in range(len(qs.isolation_mz)):
            key = (round(qs.isolation_mz[i], 1), round(f.rt, 3))
            if key in spec_lookup:
                scan_mid = (qs.scan_starts[i] + qs.scan_ends[i]) / 2.0
                pairs.append((scan_mid, spec_lookup[key]))

    if len(pairs) < 10:
        logger.warning("Too few calibration pairs (%d), using defaults", len(pairs))
        return 1.3, -0.00074

    scans = np.array([p[0] for p in pairs])
    ims = np.array([p[1] for p in pairs])
    coeffs = np.polyfit(scans, ims, 1)  # slope, intercept
    slope, intercept = coeffs[0], coeffs[1]

    # Verify quality
    pred = np.polyval(coeffs, scans)
    rmse = np.sqrt(np.mean((ims - pred) ** 2))
    logger.info(
        "IM calibration: 1/K0 = %.6f + %.8f * scan (RMSE=%.6f, n=%d pairs)",
        intercept, slope, rmse, len(pairs),
    )

    return intercept, slope


class MobilityFilteredExtractor:
    """Extract mobility-filtered pseudo-DDA spectra from a timsTOF .d folder."""

    def __init__(
        self,
        d_folder,
        mobility_half_width=3,
        mobilogram_prominence=100.0,
        min_fragments=5,
        mz_merge_tol=0.01,
        max_spectra_per_window=5,
    ):
        self.d_folder = Path(d_folder)
        self.mobility_half_width = mobility_half_width
        self.mobilogram_prominence = mobilogram_prominence
        self.min_fragments = min_fragments
        self.mz_merge_tol = mz_merge_tol
        self.max_spectra_per_window = max_spectra_per_window

        tdf_path = self.d_folder / "analysis.tdf"
        if not tdf_path.exists():
            raise FileNotFoundError(f"No analysis.tdf in {d_folder}")

        # Load m/z calibration
        self.mz_lower, self.mz_upper, self.dns = _read_mz_calibration(d_folder)
        logger.info(
            "m/z calibration: lower=%.2f, upper=%.2f, dns=%d",
            self.mz_lower, self.mz_upper, self.dns,
        )

        # Load scan-to-1/K0 calibration (derived from timsrust Spectrum IM values)
        self.im_intercept, self.im_slope = _calibrate_scan_to_im(d_folder)

    def scan_to_1k0(self, scan_num):
        """Convert scan number to 1/K0 using calibrated linear model."""
        return self.im_intercept + self.im_slope * scan_num

    def _build_mobilogram(self, scan_offsets, intensities, scan_begin, scan_end):
        """Build intensity-vs-scan mobilogram for an isolation window.

        Args:
            scan_offsets: Frame.scan_offsets array (length = num_scans + 1)
            intensities: Frame.intensities flat array
            scan_begin: First scan index of isolation window
            scan_end: Last scan index (exclusive) of isolation window

        Returns:
            1D array of summed intensity per scan in [scan_begin, scan_end).
        """
        n_scans = scan_end - scan_begin
        if n_scans <= 0:
            return np.zeros(0)

        mobilogram = np.zeros(n_scans, dtype=np.float64)
        for offset in range(n_scans):
            scan_idx = scan_begin + offset
            if scan_idx + 1 >= len(scan_offsets):
                break
            s_start = scan_offsets[scan_idx]
            s_end = scan_offsets[scan_idx + 1]
            if s_start < s_end and s_end <= len(intensities):
                mobilogram[offset] = np.sum(intensities[s_start:s_end])

        return mobilogram

    def _extract_scan_spectrum(self, scan_offsets, tof_indices, intensities, scan_lo, scan_hi):
        """Extract and centroid-merge fragment peaks from a narrow scan range.

        Args:
            scan_offsets, tof_indices, intensities: Frame arrays
            scan_lo, scan_hi: Inclusive scan range

        Returns:
            (mz_array, intensity_array) after centroid merging and tof-to-mz conversion.
        """
        all_tof = []
        all_int = []

        for scan_idx in range(scan_lo, scan_hi + 1):
            if scan_idx + 1 >= len(scan_offsets):
                break
            s_start = scan_offsets[scan_idx]
            s_end = scan_offsets[scan_idx + 1]
            if s_start >= s_end:
                continue
            all_tof.append(tof_indices[s_start:s_end])
            all_int.append(intensities[s_start:s_end])

        if not all_tof:
            return np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        tof = np.concatenate(all_tof)
        intensity = np.concatenate(all_int).astype(np.float64)

        if len(tof) == 0:
            return np.array([], dtype=np.float32), np.array([], dtype=np.float32)

        # Convert tof to m/z
        mz = _tof_to_mz(tof, self.mz_lower, self.mz_upper, self.dns)

        # Centroid merge: sort by m/z, merge peaks within tolerance
        order = np.argsort(mz)
        mz = mz[order]
        intensity = intensity[order]

        merged_mz = []
        merged_int = []
        i = 0
        while i < len(mz):
            j = i + 1
            while j < len(mz) and (mz[j] - mz[i]) < self.mz_merge_tol:
                j += 1
            cluster_mz = mz[i:j]
            cluster_int = intensity[i:j]
            total_int = np.sum(cluster_int)
            avg_mz = np.average(cluster_mz, weights=cluster_int)
            merged_mz.append(avg_mz)
            merged_int.append(total_int)
            i = j

        return np.array(merged_mz, dtype=np.float32), np.array(merged_int, dtype=np.float32)

    def extract_filtered_spectra(self):
        """Main extraction method.

        Returns list of dicts:
        {
            'mz': np.ndarray,
            'intensity': np.ndarray,
            'rt': float (minutes),
            'precursor_mz': float,
            'isolation_width': float,
            'im': float (1/K0),
            'charge': int (0 = unknown),
            'mobility_filtered': bool,
        }
        """
        import timsrust_pyo3 as timsrust

        reader = timsrust.FrameReader(str(self.d_folder))
        dia_frames = reader.read_dia_frames()
        logger.info("Read %d DIA frames from %s", len(dia_frames), self.d_folder.name)

        spectra = []
        n_filtered = 0
        n_fallback = 0
        n_skipped_few_frags = 0

        for frame_idx, frame in enumerate(dia_frames):
            if frame_idx % 1000 == 0 and frame_idx > 0:
                logger.info(
                    "  Frame %d/%d: %d spectra so far (%d filtered, %d fallback)",
                    frame_idx, len(dia_frames), len(spectra), n_filtered, n_fallback,
                )

            qs = frame.quadrupole_settings
            scan_offsets = np.array(frame.scan_offsets)
            tof_indices = np.array(frame.tof_indices)
            intensities = np.array(frame.intensities, dtype=np.float64)
            rt_min = frame.rt  # timsrust gives minutes

            n_windows = len(qs.isolation_mz)
            for win_idx in range(n_windows):
                iso_mz = qs.isolation_mz[win_idx]
                iso_width = qs.isolation_width[win_idx]
                scan_begin = qs.scan_starts[win_idx]
                scan_end = qs.scan_ends[win_idx]

                if scan_end <= scan_begin:
                    continue

                # Step 1: Build mobilogram
                mobilogram = self._build_mobilogram(
                    scan_offsets, intensities, scan_begin, scan_end
                )
                if len(mobilogram) == 0 or np.max(mobilogram) == 0:
                    continue

                # Step 2: Find mobility peaks
                peak_indices, properties = find_peaks(
                    mobilogram,
                    prominence=self.mobilogram_prominence,
                    distance=self.mobility_half_width,
                )

                if len(peak_indices) == 0:
                    # No clear mobility peaks -- fallback to summing entire window
                    mz, intensity = self._extract_scan_spectrum(
                        scan_offsets, tof_indices, intensities, scan_begin, scan_end - 1
                    )
                    if len(mz) >= self.min_fragments:
                        im_val = self.scan_to_1k0((scan_begin + scan_end) / 2.0)
                        spectra.append({
                            "mz": mz,
                            "intensity": intensity,
                            "rt": rt_min,
                            "precursor_mz": iso_mz,
                            "isolation_width": iso_width,
                            "im": im_val,
                            "charge": 0,
                            "mobility_filtered": False,
                        })
                        n_fallback += 1
                    else:
                        n_skipped_few_frags += 1
                    continue

                # Cap number of spectra per window
                if len(peak_indices) > self.max_spectra_per_window:
                    top_idx = np.argsort(properties["prominences"])[::-1]
                    peak_indices = peak_indices[top_idx[: self.max_spectra_per_window]]

                # Step 3: Extract filtered spectrum for each mobility peak
                for peak_idx in peak_indices:
                    scan_center = scan_begin + peak_idx
                    scan_lo = max(scan_begin, scan_center - self.mobility_half_width)
                    scan_hi = min(scan_end - 1, scan_center + self.mobility_half_width)

                    mz, intensity = self._extract_scan_spectrum(
                        scan_offsets, tof_indices, intensities, scan_lo, scan_hi
                    )

                    if len(mz) < self.min_fragments:
                        n_skipped_few_frags += 1
                        continue

                    im_val = self.scan_to_1k0(scan_center)

                    spectra.append({
                        "mz": mz,
                        "intensity": intensity,
                        "rt": rt_min,
                        "precursor_mz": iso_mz,
                        "isolation_width": iso_width,
                        "im": im_val,
                        "charge": 0,
                        "mobility_filtered": True,
                    })
                    n_filtered += 1

        logger.info(
            "Extraction complete: %d total spectra (%d mobility-filtered, %d fallback, %d skipped <min_frags)",
            len(spectra), n_filtered, n_fallback, n_skipped_few_frags,
        )

        return spectra

    def extract_ms1_frames(self, top_n=150):
        """Read MS1 frames with tof-to-mz conversion.

        Returns list of (rt_min, mz_array, intensity_array) tuples
        with sqrt(sqrt) normalization matching Cascadia's augment.py MS1 path.
        """
        import timsrust_pyo3 as timsrust

        reader = timsrust.FrameReader(str(self.d_folder))
        ms1_frames = reader.read_ms1_frames()
        logger.info("Read %d MS1 frames from %s", len(ms1_frames), self.d_folder.name)

        ms1_list = []
        for frame in ms1_frames:
            tof = np.array(frame.tof_indices, dtype=np.float64)
            ints = np.array(frame.intensities, dtype=np.float64)
            if len(tof) == 0:
                continue

            mzs = _tof_to_mz(tof, self.mz_lower, self.mz_upper, self.dns)

            # Sort by intensity, take top_n, re-sort by m/z
            if len(ints) > top_n:
                sorted_idx = np.argsort(ints)[-top_n:]
                ints = ints[sorted_idx]
                mzs = mzs[sorted_idx]

            sorted_mz = np.argsort(mzs)
            ints = ints[sorted_mz]
            mzs = mzs[sorted_mz]

            # Apply sqrt(sqrt) + normalize (Cascadia MS1 normalization)
            ints = ints ** 0.5
            ints = ints ** 0.5
            if len(ints) > 0 and np.max(ints) > 0:
                ints = ints / np.max(ints)

            ms1_list.append((frame.rt, mzs.astype(np.float32), ints.astype(np.float32)))

        return ms1_list


def write_asf_5col(output_path, spectra, ms1_list, top_n_ms2=150, scan_width=1):
    """Write 5-column ASF matching Cascadia's IM-aware transformer format.

    Format per spectrum block:
        BEGIN IONS
        TITLE=<scan_idx>
        PEPMASS=<precursor_mz>
        CHARGE=<charge>
        SCAN=<scan_idx>
        RT=<rt_minutes>
        <m/z>\t<intensity>\t<rt_offset>\t<ms_level>\t<1/K0>
        ...
        END IONS

    MS2 peaks: ms_level=2, im=actual 1/K0 at mobility peak center
    MS1 peaks: ms_level=1, im=0.0 (MS1 is not mobility-resolved here)
    """
    # Build MS1 lookup by RT
    ms1_by_rt_bin = {}
    for rt_min, mzs, ints in ms1_list:
        rt_bin = int(rt_min * 6)  # 10-second bins
        ms1_by_rt_bin.setdefault(rt_bin, []).append((rt_min, mzs, ints))

    n_written = 0
    with open(output_path, "w") as out:
        for spec in spectra:
            mz = spec["mz"]
            intensity = spec["intensity"]
            rt = spec["rt"]
            prec_mz = spec["precursor_mz"]
            im = spec["im"]
            charge = spec.get("charge", 0)

            # Apply MS2 normalization: sqrt + normalize (Cascadia MS2 path)
            int_norm = intensity.copy().astype(np.float64)
            if len(int_norm) > top_n_ms2:
                top_idx = np.argsort(int_norm)[-top_n_ms2:]
                int_norm = int_norm[top_idx]
                mz = mz[top_idx]
                sorted_mz = np.argsort(mz)
                int_norm = int_norm[sorted_mz]
                mz = mz[sorted_mz]

            int_norm = int_norm ** 0.5
            if len(int_norm) > 0 and np.max(int_norm) > 0:
                int_norm = int_norm / np.max(int_norm)

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
                continue  # skip spectra without MS1 context

            n_written += 1
            out.write("BEGIN IONS\n")
            out.write(f"TITLE={n_written}\n")
            out.write(f"PEPMASS={prec_mz}\n")
            out.write(f"CHARGE={charge}\n")
            out.write(f"SCAN={n_written}\n")
            out.write(f"RT={rt}\n")

            # MS2 peaks with real 1/K0
            for m, i in zip(mz, int_norm):
                out.write(f"{m:.4f}\t{i:.6f}\t0.0\t2\t{im:.4f}\n")

            # MS1 peaks with IM=0.0
            ms1_rt, ms1_mz, ms1_int = best_ms1
            iso_half = spec.get("isolation_width", 25.0) / 2.0 + 1.0
            rt_offset = ms1_rt - rt
            for m, i in zip(ms1_mz, ms1_int):
                if abs(m - prec_mz) < iso_half + 50:  # wider range for MS1
                    out.write(f"{m:.4f}\t{i:.6f}\t{rt_offset:.4f}\t1\t0.0000\n")

            out.write("END IONS\n")

    return n_written


def print_stats(spectra, output_path=None):
    """Print extraction statistics."""
    if not spectra:
        print("No spectra extracted!")
        return

    n_total = len(spectra)
    n_filtered = sum(1 for s in spectra if s["mobility_filtered"])
    n_fallback = n_total - n_filtered
    frag_counts = [len(s["mz"]) for s in spectra]
    im_values = [s["im"] for s in spectra if s["mobility_filtered"]]
    all_ims = [s["im"] for s in spectra]

    print(f"\n{'='*60}")
    print(f"Mobility-Filtered Extraction Summary")
    print(f"{'='*60}")
    print(f"Total spectra:        {n_total:,}")
    print(f"  Mobility-filtered:  {n_filtered:,} ({100*n_filtered/n_total:.1f}%)")
    print(f"  Fallback (summed):  {n_fallback:,} ({100*n_fallback/n_total:.1f}%)")
    print(f"\nFragment counts:")
    print(f"  Median: {np.median(frag_counts):.0f}")
    print(f"  Mean:   {np.mean(frag_counts):.1f}")
    print(f"  Range:  {min(frag_counts)}-{max(frag_counts)}")
    if im_values:
        print(f"\n1/K0 distribution (filtered only):")
        print(f"  Range:  {min(im_values):.4f} - {max(im_values):.4f}")
        print(f"  Median: {np.median(im_values):.4f}")
    print(f"\n1/K0 distribution (all):")
    print(f"  Range:  {min(all_ims):.4f} - {max(all_ims):.4f}")
    print(f"  Median: {np.median(all_ims):.4f}")

    mem_bytes = sum(s["mz"].nbytes + s["intensity"].nbytes for s in spectra)
    print(f"\nMemory: ~{mem_bytes / 1e6:.1f} MB")

    if output_path:
        print(f"Output: {output_path}")
    print(f"{'='*60}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Mobility-filtered spectrum extraction from timsTOF diaPASEF .d files"
    )
    parser.add_argument("--d-folder", required=True, help="Path to .d folder")
    parser.add_argument("--output", required=True, help="Output ASF file path")
    parser.add_argument("--mobility-half-width", type=int, default=3, help="Scans on each side of mobility peak (default: 3)")
    parser.add_argument("--mobilogram-prominence", type=float, default=100.0, help="Min prominence for peak detection (default: 100)")
    parser.add_argument("--min-fragments", type=int, default=5, help="Min fragment peaks per spectrum (default: 5)")
    parser.add_argument("--mz-merge-tol", type=float, default=0.01, help="m/z merge tolerance in Da (default: 0.01)")
    parser.add_argument("--max-spectra-per-window", type=int, default=5, help="Max spectra per isolation window (default: 5)")
    parser.add_argument("--stats-only", action="store_true", help="Print stats without writing ASF")
    parser.add_argument("--compare-collapsed", action="store_true", help="Also extract collapsed spectra for comparison")

    args = parser.parse_args()

    t0 = time.time()

    # Extract mobility-filtered spectra
    extractor = MobilityFilteredExtractor(
        args.d_folder,
        mobility_half_width=args.mobility_half_width,
        mobilogram_prominence=args.mobilogram_prominence,
        min_fragments=args.min_fragments,
        mz_merge_tol=args.mz_merge_tol,
        max_spectra_per_window=args.max_spectra_per_window,
    )

    spectra = extractor.extract_filtered_spectra()
    print_stats(spectra, args.output)

    if args.compare_collapsed:
        # Compare with collapsed extraction
        import timsrust_pyo3 as timsrust
        collapsed = timsrust.read_all_spectra(str(args.d_folder))
        n_collapsed = len(collapsed)
        print(f"\nComparison with collapsed extraction:")
        print(f"  Collapsed spectra: {n_collapsed:,}")
        print(f"  Mobility-filtered: {len(spectra):,}")
        print(f"  Ratio: {len(spectra)/max(n_collapsed,1):.1f}x")

    if not args.stats_only:
        # Read MS1 frames for ASF context
        ms1_list = extractor.extract_ms1_frames()
        logger.info("Read %d MS1 frames for ASF context", len(ms1_list))

        # Write ASF
        n_written = write_asf_5col(args.output, spectra, ms1_list)
        logger.info("Wrote %d spectra to %s", n_written, args.output)

    elapsed = time.time() - t0
    logger.info("Done in %.1f min", elapsed / 60)


if __name__ == "__main__":
    main()
