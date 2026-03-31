"""
routing_patch.py — Shows the exact code change needed in cascadia/cascadia.py
to add Bruker .d native support via bruker_augment.py.

This is NOT meant to be executed directly. It documents the patch to apply
to the Cascadia fork's main entry point.

Location in Cascadia source:
    cascadia/cascadia.py — find the data loading section where augment.py
    is called (look for `augment_spectra` or `asf_output_path`).

Replace the existing augmentation call with the format-routing block below.
"""

# =============================================================================
# PATCH: Add to cascadia/cascadia.py
#
# Find the line that calls:
#     from cascadia.augment import augment_spectra
#     n_spectra = augment_spectra(str(input_path), str(asf_output_path), config)
#
# Replace with the following block:
# =============================================================================

from pathlib import Path

# --- BEGIN PATCH ---

# Detect input format and route to appropriate augmentation function
input_path = Path(config.input)

if input_path.suffix == ".d" or (
    input_path.is_dir() and (input_path / "analysis.tdf").exists()
):
    # Bruker timsTOF .d directory — use native loader (no mzML conversion)
    from cascadia.bruker_augment import augment_spectra_bruker

    n_spectra = augment_spectra_bruker(
        str(input_path),
        str(asf_output_path),
        augmentation_width=config.augmentation_width,
    )
    logger.info(
        f"Bruker native loader: {n_spectra} augmented spectra "
        f"from {input_path.name}"
    )
else:
    # mzML or other format — use existing augment.py path (unchanged)
    from cascadia.augment import augment_spectra

    n_spectra = augment_spectra(
        str(input_path), str(asf_output_path), config
    )

# --- END PATCH ---

# =============================================================================
# ADDITIONAL CHANGE: pyproject.toml
#
# Add Bruker as an optional dependency group:
#
#   [project.optional-dependencies]
#   bruker = ["timsrust_pyo3>=0.4.0"]
#
# Install with: pip install cascadia[bruker]
#
# For HPC deployment (conda):
#   conda activate cascadia
#   pip install timsrust_pyo3
# =============================================================================
