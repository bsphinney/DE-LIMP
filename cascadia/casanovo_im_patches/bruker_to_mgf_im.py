"""
bruker_to_mgf_im.py — Convert Bruker ddaPASEF .d files to MGF with ion mobility.

Each ddaPASEF spectrum has a precursor with a known 1/K0 value from the
TIMS separation. This script writes MGF files with an extra ION_MOBILITY
field so that Casanovo's IM-enhanced model can use this information.

Usage:
    python bruker_to_mgf_im.py input.d output.mgf [--top-n 150] [--min-mz 50]

Requires:
    timsrust_pyo3 (pip install timsrust-pyo3)

Output MGF format:
    BEGIN IONS
    TITLE=scan_1
    PEPMASS=500.2500
    CHARGE=2+
    RTINSECONDS=120.5
    ION_MOBILITY=1.0500
    SCANS=1
    200.1234 1000.0
    300.2345 2500.0
    ...
    END IONS
"""

import argparse
import logging
import os
import sys
import sqlite3
from pathlib import Path

import numpy as np

logger = logging.getLogger(__name__)


def _read_tdf_metadata(d_folder):
    """Read metadata from analysis.tdf SQLite database.

    Returns
    -------
    dict with keys: mz_lower, mz_upper, digitizer_samples, instrument
    """
    tdf_path = os.path.join(str(d_folder), "analysis.tdf")
    result = {
        "mz_lower": 100.0,
        "mz_upper": 1700.0,
        "digitizer_samples": 631025,
        "instrument": "unknown",
    }
    if not os.path.exists(tdf_path):
        logger.warning("analysis.tdf not found in %s", d_folder)
        return result

    try:
        conn = sqlite3.connect(f"file:{tdf_path}?mode=ro", uri=True)
        meta = dict(
            conn.execute("SELECT Key, Value FROM GlobalMetadata").fetchall()
        )
        conn.close()
        result["mz_lower"] = float(meta.get("MzAcqRangeLower", 100.0))
        result["mz_upper"] = float(meta.get("MzAcqRangeUpper", 1700.0))
        result["digitizer_samples"] = int(
            meta.get("DigitizerNumSamples", 631025)
        )
        result["instrument"] = meta.get("InstrumentName", "unknown")
    except Exception as e:
        logger.warning("Failed to read TDF metadata: %s", e)

    return result


def convert_bruker_to_mgf_im(
    d_folder,
    output_mgf,
    top_n=150,
    min_mz=50.0,
    max_mz=2500.0,
    min_intensity_frac=0.01,
):
    """
    Convert a Bruker ddaPASEF .d folder to an MGF file with ION_MOBILITY.

    Parameters
    ----------
    d_folder : str or Path
        Path to the .d folder.
    output_mgf : str or Path
        Path for the output MGF file.
    top_n : int
        Keep only the top N most intense peaks per spectrum.
    min_mz : float
        Minimum m/z to include.
    max_mz : float
        Maximum m/z to include.
    min_intensity_frac : float
        Remove peaks below this fraction of the base peak intensity.

    Returns
    -------
    dict
        Summary statistics: n_spectra, n_with_im, instrument, etc.
    """
    try:
        import timsrust_pyo3 as timsrust
    except ImportError:
        logger.error(
            "timsrust_pyo3 is required for Bruker .d conversion. "
            "Install with: pip install timsrust-pyo3"
        )
        raise

    d_folder = str(d_folder)
    meta = _read_tdf_metadata(d_folder)
    logger.info(
        "Converting %s (instrument: %s, m/z range: %.0f-%.0f)",
        os.path.basename(d_folder),
        meta["instrument"],
        meta["mz_lower"],
        meta["mz_upper"],
    )

    # Read all DDA spectra using timsrust
    reader = timsrust.SpectrumReader(d_folder)
    n_spectra = 0
    n_with_im = 0
    n_skipped = 0

    with open(output_mgf, "w") as out:
        for spec in reader:
            # Get precursor information
            precursor = spec.precursor
            if precursor is None or precursor.mz is None:
                n_skipped += 1
                continue

            prec_mz = float(precursor.mz)
            if prec_mz <= 0:
                n_skipped += 1
                continue

            # Precursor charge
            prec_charge = (
                int(precursor.charge) if hasattr(precursor, "charge") and precursor.charge else 0
            )
            if prec_charge <= 0:
                # Default to charge 2 if unknown (common in DDA)
                prec_charge = 2

            # Ion mobility (1/K0)
            im_value = 0.0
            if (
                hasattr(precursor, "im")
                and precursor.im is not None
            ):
                im_value = float(precursor.im)
                if im_value > 0:
                    n_with_im += 1

            # Retention time (timsrust returns seconds or minutes depending on version)
            rt_seconds = float(precursor.rt) if hasattr(precursor, "rt") and precursor.rt else 0.0
            # timsrust typically returns RT in seconds already
            # But if values seem like minutes (< 200), convert
            if rt_seconds > 0 and rt_seconds < 200:
                rt_seconds = rt_seconds * 60.0

            # Get m/z and intensity arrays
            mzs = np.array(spec.mz_values, dtype=np.float64)
            intensities = np.array(spec.intensities, dtype=np.float64)

            if len(mzs) == 0:
                n_skipped += 1
                continue

            # Apply m/z range filter
            mask = (mzs >= min_mz) & (mzs <= max_mz)
            mzs = mzs[mask]
            intensities = intensities[mask]

            if len(mzs) == 0:
                n_skipped += 1
                continue

            # Remove precursor peak region (+/- 2 Da)
            prec_mask = np.abs(mzs - prec_mz) > 2.0
            mzs = mzs[prec_mask]
            intensities = intensities[prec_mask]

            if len(mzs) == 0:
                n_skipped += 1
                continue

            # Filter by minimum relative intensity
            max_int = np.max(intensities)
            if max_int > 0:
                rel_mask = intensities >= min_intensity_frac * max_int
                mzs = mzs[rel_mask]
                intensities = intensities[rel_mask]

            if len(mzs) == 0:
                n_skipped += 1
                continue

            # Keep top N by intensity
            if len(mzs) > top_n:
                top_idx = np.argsort(intensities)[-top_n:]
                mzs = mzs[top_idx]
                intensities = intensities[top_idx]

            # Sort by m/z for output
            sort_idx = np.argsort(mzs)
            mzs = mzs[sort_idx]
            intensities = intensities[sort_idx]

            n_spectra += 1

            # Write MGF entry
            out.write("BEGIN IONS\n")
            out.write(f"TITLE=scan_{n_spectra}\n")
            out.write(f"PEPMASS={prec_mz:.6f}\n")
            out.write(f"CHARGE={prec_charge}+\n")
            out.write(f"RTINSECONDS={rt_seconds:.2f}\n")
            out.write(f"ION_MOBILITY={im_value:.4f}\n")
            out.write(f"SCANS={n_spectra}\n")

            for mz, intensity in zip(mzs, intensities):
                out.write(f"{mz:.6f} {intensity:.1f}\n")

            out.write("END IONS\n\n")

    stats = {
        "n_spectra": n_spectra,
        "n_with_im": n_with_im,
        "n_skipped": n_skipped,
        "pct_with_im": 100 * n_with_im / max(n_spectra, 1),
        "instrument": meta["instrument"],
        "output_file": str(output_mgf),
    }

    logger.info(
        "Wrote %d spectra to %s (%d with IM, %.0f%%; %d skipped)",
        n_spectra,
        os.path.basename(str(output_mgf)),
        n_with_im,
        stats["pct_with_im"],
        n_skipped,
    )

    return stats


def main():
    """Command-line entry point."""
    parser = argparse.ArgumentParser(
        description="Convert Bruker ddaPASEF .d to MGF with ion mobility"
    )
    parser.add_argument("input_d", help="Path to .d folder")
    parser.add_argument("output_mgf", help="Path for output .mgf file")
    parser.add_argument(
        "--top-n",
        type=int,
        default=150,
        help="Keep top N peaks per spectrum (default: 150)",
    )
    parser.add_argument(
        "--min-mz",
        type=float,
        default=50.0,
        help="Minimum m/z (default: 50.0)",
    )
    parser.add_argument(
        "--max-mz",
        type=float,
        default=2500.0,
        help="Maximum m/z (default: 2500.0)",
    )
    parser.add_argument(
        "--min-intensity",
        type=float,
        default=0.01,
        help="Min intensity as fraction of base peak (default: 0.01)",
    )

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    stats = convert_bruker_to_mgf_im(
        args.input_d,
        args.output_mgf,
        top_n=args.top_n,
        min_mz=args.min_mz,
        max_mz=args.max_mz,
        min_intensity_frac=args.min_intensity,
    )

    print(f"\nConversion complete:")
    print(f"  Spectra written: {stats['n_spectra']}")
    print(f"  With ion mobility: {stats['n_with_im']} ({stats['pct_with_im']:.0f}%)")
    print(f"  Skipped: {stats['n_skipped']}")
    print(f"  Instrument: {stats['instrument']}")
    print(f"  Output: {stats['output_file']}")


if __name__ == "__main__":
    main()
