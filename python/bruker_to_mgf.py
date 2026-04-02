#!/usr/bin/env python3
"""
bruker_to_mgf.py — Convert Bruker timsTOF .d files to MGF for Casanovo de novo sequencing.

Uses timsrust_pyo3 for native .d reading (fast, no mzML intermediate).
Writes standard MGF format that Casanovo accepts directly.

Usage:
    python bruker_to_mgf.py /path/to/sample.d /path/to/output.mgf
    python bruker_to_mgf.py /path/to/raw_dir/ /path/to/mgf_dir/ --batch

Requirements:
    pip install timsrust_pyo3>=0.4.0
"""

import argparse
import logging
import os
import sys
import time
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


def convert_d_to_mgf(d_folder: str, mgf_path: str, min_peaks: int = 6,
                      min_intensity: float = 0.0) -> dict:
    """
    Convert a single Bruker .d directory to MGF format.

    Args:
        d_folder: Path to .d directory
        mgf_path: Output MGF file path
        min_peaks: Minimum number of fragment peaks to include a spectrum
        min_intensity: Minimum base peak intensity to include

    Returns:
        dict with conversion stats: n_total, n_ms2, n_written, elapsed_sec
    """
    import timsrust_pyo3 as timsrust

    path = Path(d_folder)
    if not path.exists():
        raise FileNotFoundError(f"Bruker .d directory not found: {d_folder}")
    if not (path / "analysis.tdf").exists():
        raise FileNotFoundError(f"Missing analysis.tdf in {d_folder}")

    logger.info("Reading %s ...", path.name)
    t0 = time.time()

    try:
        all_spectra = timsrust.read_all_spectra(str(path))
    except Exception as e:
        raise RuntimeError(f"timsrust failed to read {path.name}: {e}") from e

    n_total = len(all_spectra)
    n_ms2 = 0
    n_written = 0

    os.makedirs(os.path.dirname(os.path.abspath(mgf_path)), exist_ok=True)

    with open(mgf_path, "w") as f:
        for i, spec in enumerate(all_spectra):
            # Skip MS1 spectra (no precursor m/z)
            if spec.precursor.mz is None or spec.precursor.mz == 0.0:
                continue
            n_ms2 += 1

            mz_arr = np.array(spec.mz_values, dtype=np.float64)
            int_arr = np.array(spec.intensities, dtype=np.float64)

            # Skip spectra with too few peaks
            if len(mz_arr) < min_peaks:
                continue

            # Skip low-intensity spectra
            if min_intensity > 0 and np.max(int_arr) < min_intensity:
                continue

            charge = int(spec.precursor.charge) if spec.precursor.charge else 2
            rt_sec = spec.precursor.rt * 60.0  # timsrust RT is in minutes

            f.write("BEGIN IONS\n")
            f.write(f"TITLE={path.stem}:scan={i}\n")
            f.write(f"PEPMASS={spec.precursor.mz:.6f}\n")
            f.write(f"CHARGE={charge}+\n")
            f.write(f"RTINSECONDS={rt_sec:.2f}\n")

            for mz, intensity in zip(mz_arr, int_arr):
                f.write(f"{mz:.5f} {intensity:.1f}\n")

            f.write("END IONS\n\n")
            n_written += 1

    elapsed = time.time() - t0
    logger.info(
        "%s: %d total spectra, %d MS2, %d written to MGF (%.1f sec)",
        path.name, n_total, n_ms2, n_written, elapsed
    )

    return {
        "n_total": n_total,
        "n_ms2": n_ms2,
        "n_written": n_written,
        "elapsed_sec": round(elapsed, 1),
    }


def batch_convert(raw_dir: str, mgf_dir: str, min_peaks: int = 6) -> list:
    """Convert all .d directories in raw_dir to MGF files in mgf_dir."""
    raw_path = Path(raw_dir)
    d_dirs = sorted(raw_path.glob("*.d"))

    if not d_dirs:
        logger.warning("No .d directories found in %s", raw_dir)
        return []

    os.makedirs(mgf_dir, exist_ok=True)
    results = []

    for d_dir in d_dirs:
        mgf_file = os.path.join(mgf_dir, d_dir.stem + ".mgf")
        try:
            stats = convert_d_to_mgf(str(d_dir), mgf_file, min_peaks=min_peaks)
            stats["file"] = d_dir.name
            stats["mgf_path"] = mgf_file
            results.append(stats)
        except Exception as e:
            logger.error("Failed to convert %s: %s", d_dir.name, e)
            results.append({
                "file": d_dir.name,
                "error": str(e),
                "n_written": 0,
            })

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Bruker .d files to MGF for Casanovo"
    )
    parser.add_argument("input", help="Path to .d directory or parent directory (with --batch)")
    parser.add_argument("output", help="Output MGF file path or directory (with --batch)")
    parser.add_argument("--batch", action="store_true",
                        help="Convert all .d dirs in input directory")
    parser.add_argument("--min-peaks", type=int, default=6,
                        help="Minimum fragment peaks per spectrum (default: 6)")
    parser.add_argument("--min-intensity", type=float, default=0.0,
                        help="Minimum base peak intensity (default: 0)")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="Verbose logging")

    args = parser.parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s"
    )

    if args.batch:
        results = batch_convert(args.input, args.output, min_peaks=args.min_peaks)
        total = sum(r.get("n_written", 0) for r in results)
        errors = sum(1 for r in results if "error" in r)
        print(f"\nBatch complete: {len(results)} files, {total} spectra written, {errors} errors")
    else:
        stats = convert_d_to_mgf(
            args.input, args.output,
            min_peaks=args.min_peaks,
            min_intensity=args.min_intensity
        )
        print(f"Converted: {stats['n_written']} spectra in {stats['elapsed_sec']}s")
