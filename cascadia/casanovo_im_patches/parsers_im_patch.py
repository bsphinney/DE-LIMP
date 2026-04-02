"""
Patch for depthcharge's MgfParser to extract ION_MOBILITY from MGF files.

Usage:
    Import and use IMEnabledMgfParser in place of the standard MgfParser
    when building spectrum indexes from IM-enabled MGF files.

    The ION_MOBILITY values are stored in the parser's `ion_mobilities` list,
    parallel to `precursor_mz`, `precursor_charge`, and `scan_id`.
"""

import logging
import numpy as np
from pathlib import Path
from pyteomics.mgf import MGF

logger = logging.getLogger("casanovo")


class IMEnabledMgfParser:
    """
    Parse MGF files and extract ION_MOBILITY values alongside standard fields.

    This is a standalone parser that does NOT subclass depthcharge's BaseParser
    (to avoid modifying the installed package). Instead, it produces parallel
    arrays that can be used to construct an im_array for IMSpectrumDataset.

    The parser reads the MGF file, extracts standard fields (m/z, intensity,
    PEPMASS, CHARGE), and additionally looks for an ION_MOBILITY field in each
    spectrum's parameters.

    Parameters
    ----------
    ms_data_file : str or Path
        The MGF file to parse.
    valid_charge : set of int, optional
        Only keep spectra with these precursor charges.
    """

    def __init__(self, ms_data_file, valid_charge=None):
        self.path = Path(ms_data_file)
        self.valid_charge = (
            set(valid_charge) if valid_charge is not None else None
        )
        self.ion_mobilities = []
        self.precursor_mz = []
        self.precursor_charge = []
        self.scan_id = []
        self._counter = -1

    def read(self):
        """
        Read the MGF file and extract IM values.

        Returns
        -------
        self
        """
        n_with_im = 0
        n_total = 0

        with MGF(str(self.path)) as reader:
            for spectrum in reader:
                self._counter += 1
                params = spectrum.get("params", {})

                precursor_mz = float(params["pepmass"][0])
                precursor_charge = int(params.get("charge", [0])[0])

                # Skip invalid charges
                if (
                    self.valid_charge is not None
                    and precursor_charge not in self.valid_charge
                ):
                    continue

                # Extract ION_MOBILITY if present
                im_value = 0.0
                for key in ("ion_mobility", "ionmobility", "im"):
                    if key in params:
                        try:
                            im_value = float(params[key])
                            n_with_im += 1
                        except (ValueError, TypeError):
                            pass
                        break

                self.precursor_mz.append(precursor_mz)
                self.precursor_charge.append(precursor_charge)
                self.ion_mobilities.append(im_value)
                self.scan_id.append(self._counter)
                n_total += 1

        self.ion_mobilities = np.array(self.ion_mobilities, dtype=np.float64)
        self.precursor_mz = np.array(self.precursor_mz, dtype=np.float64)
        self.precursor_charge = np.array(
            self.precursor_charge, dtype=np.uint8
        )
        self.scan_id = np.array(self.scan_id)

        logger.info(
            "Parsed %s: %d spectra, %d with ion mobility values (%.0f%%)",
            self.path.name,
            n_total,
            n_with_im,
            100 * n_with_im / max(n_total, 1),
        )
        return self


def extract_im_from_mgf(mgf_paths, valid_charge=None):
    """
    Extract ion mobility values from one or more MGF files.

    This produces a single im_array that parallels the depthcharge
    SpectrumIndex built from the same MGF files (in the same order).

    Parameters
    ----------
    mgf_paths : list of str
        Paths to MGF files, in the same order as passed to SpectrumIndex.
    valid_charge : set of int, optional
        Only include spectra with these precursor charges.

    Returns
    -------
    im_array : numpy.ndarray of shape (n_spectra,)
        The 1/K0 values for each spectrum, aligned with the SpectrumIndex.
    """
    all_im = []
    for path in mgf_paths:
        parser = IMEnabledMgfParser(path, valid_charge=valid_charge)
        parser.read()
        all_im.append(parser.ion_mobilities)

    if len(all_im) == 0:
        return np.array([], dtype=np.float64)

    return np.concatenate(all_im)
