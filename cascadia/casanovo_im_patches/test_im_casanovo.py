"""
Test script: Verify Casanovo IM integration produces identical output
to the original model when IM values are zero.

Usage (on HPC):
    cd /quobyte/proteomics-grp/de-limp/casanovo_im
    /quobyte/proteomics-grp/conda_envs/cassonovo_env/bin/python3 test_im_casanovo.py

This script:
1. Loads the original Casanovo v4.2.0 checkpoint (use_ion_mobility=False)
2. Runs encoder on a small synthetic batch
3. Loads the same checkpoint via from_pretrained(use_ion_mobility=True)
4. Runs encoder with IM=0 and verifies identical output
5. Runs with IM=1.05 to verify no crash
6. Tests full beam search decode
7. Tests save/load roundtrip
"""

import sys
import os
import logging
import numpy as np
import torch

# PyTorch 2.6+ changed default weights_only=True, need to allowlist numpy
torch.serialization.add_safe_globals([np.core.multiarray.scalar])

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)

# Path to the pre-trained checkpoint
CHECKPOINT = "/quobyte/proteomics-grp/bioinformatics_programs/casanovo_modles/casanovo_v4_2_0.ckpt"


def create_synthetic_batch(batch_size=4, n_peaks=50, with_im=False, im_value=0.0):
    """Create a synthetic batch for testing."""
    torch.manual_seed(42)
    mzs = torch.rand(batch_size, n_peaks) * 1600 + 200
    mzs, _ = mzs.sort(dim=1)
    intensities = torch.rand(batch_size, n_peaks)
    spectra = torch.stack([mzs, intensities], dim=2)

    precursor_mzs = torch.tensor([500.0, 600.0, 700.0, 800.0][:batch_size])
    precursor_charges = torch.tensor([2, 3, 2, 3][:batch_size])
    precursor_masses = (precursor_mzs - 1.007276) * precursor_charges

    if with_im:
        im_values = torch.full((batch_size,), im_value)
        precursors = torch.vstack(
            [precursor_masses, precursor_charges.float(), precursor_mzs, im_values]
        ).T.float()
    else:
        precursors = torch.vstack(
            [precursor_masses, precursor_charges.float(), precursor_mzs]
        ).T.float()

    return spectra, precursors


def test_checkpoint_loading():
    """Test 1: Verify the checkpoint exists and can be examined."""
    logger.info("=" * 60)
    logger.info("TEST 1: Checkpoint inspection")
    logger.info("=" * 60)

    if not os.path.exists(CHECKPOINT):
        logger.error("Checkpoint not found: %s", CHECKPOINT)
        return False

    checkpoint = torch.load(CHECKPOINT, map_location="cpu", weights_only=False)
    hparams = checkpoint.get("hyper_parameters", {})

    logger.info("Checkpoint hyperparameters:")
    for key, val in sorted(hparams.items()):
        if key != "residues":
            logger.info("  %s = %s", key, val)

    state_dict = checkpoint.get("state_dict", {})
    logger.info("State dict: %d keys", len(state_dict))

    dim_model = hparams.get("dim_model", 512)
    logger.info("dim_model: %d", dim_model)
    logger.info("TEST 1 PASSED")
    return True


def test_original_model():
    """Test 2: Load model without IM and run encoder."""
    logger.info("=" * 60)
    logger.info("TEST 2: Model without IM (backward compatible load)")
    logger.info("=" * 60)

    from casanovo.denovo.model import Spec2Pep

    # Load checkpoint manually and create model without IM
    checkpoint = torch.load(CHECKPOINT, map_location="cpu", weights_only=False)
    hparams = checkpoint.get("hyper_parameters", {})

    # Remove non-constructor params
    constructor_kwargs = dict(hparams)
    for key in ("out_writer", "tb_summarywriter"):
        constructor_kwargs.pop(key, None)

    # Explicitly set use_ion_mobility=False (original behavior)
    constructor_kwargs["use_ion_mobility"] = False

    try:
        model = Spec2Pep(**constructor_kwargs)
        # Load weights
        state_dict = checkpoint.get("state_dict", {})
        model.load_state_dict(state_dict, strict=True)
        model.eval()
        logger.info("Model without IM loaded successfully")
    except Exception as e:
        logger.error("Failed to load model: %s", e)
        import traceback
        traceback.print_exc()
        return None

    spectra, precursors = create_synthetic_batch(batch_size=2, n_peaks=30)

    with torch.no_grad():
        memories, mem_masks = model.encoder(spectra)
        logger.info("Encoder output shape: %s", memories.shape)
        logger.info("Encoder output first 5 values: %s", memories[0, 0, :5].tolist())

    logger.info("TEST 2 PASSED")
    return memories.clone(), mem_masks.clone(), spectra, precursors


def test_im_model_zero_init(orig_memories, orig_masks, spectra, precursors):
    """Test 3: IM model with zero IM produces identical output."""
    logger.info("=" * 60)
    logger.info("TEST 3: IM model zero-init backward compatibility")
    logger.info("=" * 60)

    from casanovo.denovo.model import Spec2Pep

    try:
        im_model = Spec2Pep.from_pretrained(
            CHECKPOINT,
            use_ion_mobility=True,
        )
        im_model.eval()
        logger.info("IM model loaded via from_pretrained()")
    except Exception as e:
        logger.error("Failed to load IM model: %s", e)
        import traceback
        traceback.print_exc()
        return False

    # Check IM components
    has_im_encoder = hasattr(im_model, "im_encoder")
    has_im_projector = hasattr(im_model, "im_projector")
    logger.info("Has im_encoder: %s", has_im_encoder)
    logger.info("Has im_projector: %s", has_im_projector)

    if not (has_im_encoder and has_im_projector):
        logger.error("Missing IM components!")
        return False

    # Log projector weights to verify zero-init
    proj_weight = im_model.im_projector.weight.data
    dim_model = proj_weight.shape[0]
    spectrum_part = proj_weight[:, :dim_model]
    im_part = proj_weight[:, dim_model:]

    identity_error = (spectrum_part - torch.eye(dim_model)).abs().max().item()
    im_max = im_part.abs().max().item()
    logger.info("Spectrum part identity error: %.10f", identity_error)
    logger.info("IM part max value: %.10f", im_max)

    # Run with IM=0 (4-column precursors)
    precursors_with_im = torch.cat(
        [precursors, torch.zeros(precursors.shape[0], 1)], dim=1
    )

    with torch.no_grad():
        im_memories, im_masks = im_model._encode_spectra(
            spectra, precursors_with_im
        )

    # Compare outputs
    max_diff = (orig_memories - im_memories).abs().max().item()
    mean_diff = (orig_memories - im_memories).abs().mean().item()

    logger.info("Max absolute difference: %.10f", max_diff)
    logger.info("Mean absolute difference: %.10f", mean_diff)

    tolerance = 1e-5
    if max_diff < tolerance:
        logger.info(
            "IDENTICAL within tolerance (%.0e) -- backward compatible!",
            tolerance,
        )
        logger.info("TEST 3 PASSED")
        return True
    else:
        logger.warning(
            "Outputs differ by %.6f (tolerance %.0e)", max_diff, tolerance
        )
        logger.warning("TEST 3 FAILED")
        return False


def test_im_model_nonzero():
    """Test 4: IM model with actual IM values doesn't crash."""
    logger.info("=" * 60)
    logger.info("TEST 4: IM model with nonzero IM values")
    logger.info("=" * 60)

    from casanovo.denovo.model import Spec2Pep

    im_model = Spec2Pep.from_pretrained(CHECKPOINT, use_ion_mobility=True)
    im_model.eval()

    spectra, precursors_3 = create_synthetic_batch(batch_size=2, n_peaks=30)

    prec_im0 = torch.cat([precursors_3, torch.zeros(2, 1)], dim=1)
    prec_im1 = torch.cat([precursors_3, torch.full((2, 1), 1.05)], dim=1)

    with torch.no_grad():
        mem_im0, _ = im_model._encode_spectra(spectra, prec_im0)
        mem_im1, _ = im_model._encode_spectra(spectra, prec_im1)

    max_diff = (mem_im0 - mem_im1).abs().max().item()
    logger.info("Diff between IM=0 and IM=1.05: %.10f", max_diff)

    if max_diff < 1e-5:
        logger.info(
            "As expected: zero-initialized IM projector produces identical "
            "output regardless of IM input. Fine-tuning will learn to use IM."
        )
    else:
        logger.info("IM channel is active (diff=%.6f).", max_diff)

    logger.info("TEST 4 PASSED (no crash)")
    return True


def test_full_forward():
    """Test 5: Full forward pass (beam search decode) with IM."""
    logger.info("=" * 60)
    logger.info("TEST 5: Full forward pass with beam search")
    logger.info("=" * 60)

    from casanovo.denovo.model import Spec2Pep

    im_model = Spec2Pep.from_pretrained(
        CHECKPOINT, use_ion_mobility=True, n_beams=1
    )
    im_model.eval()

    spectra, precursors_3 = create_synthetic_batch(batch_size=2, n_peaks=30)
    precursors_4 = torch.cat(
        [precursors_3, torch.full((2, 1), 1.05)], dim=1
    )

    logger.info("Running beam search decode...")
    with torch.no_grad():
        results = im_model(spectra, precursors_4)

    for i, spectrum_results in enumerate(results):
        if len(spectrum_results) > 0:
            score, aa_scores, peptide = spectrum_results[0]
            logger.info(
                "Spectrum %d: peptide=%s, score=%.4f, len=%d",
                i, peptide, score, len(peptide),
            )
        else:
            logger.info("Spectrum %d: no predictions", i)

    logger.info("TEST 5 PASSED")
    return True


def test_save_load_roundtrip():
    """Test 6: Save and reload the IM model checkpoint."""
    logger.info("=" * 60)
    logger.info("TEST 6: Save/load roundtrip")
    logger.info("=" * 60)

    from casanovo.denovo.model import Spec2Pep

    im_model = Spec2Pep.from_pretrained(CHECKPOINT, use_ion_mobility=True)

    import tempfile
    with tempfile.NamedTemporaryFile(suffix=".ckpt", delete=False) as f:
        tmp_path = f.name

    try:
        torch.save(
            {
                "state_dict": im_model.state_dict(),
                "hyper_parameters": dict(im_model.hparams),
            },
            tmp_path,
        )
        logger.info("Saved IM model to %s", tmp_path)

        # Reload using from_pretrained (which handles IM params)
        reloaded = Spec2Pep.from_pretrained(tmp_path, use_ion_mobility=True)
        reloaded.eval()
        logger.info("Reloaded model successfully")
        logger.info("Reloaded use_ion_mobility: %s", reloaded.use_ion_mobility)

        # Verify weights match
        spectra, precursors_3 = create_synthetic_batch(batch_size=1, n_peaks=20)
        precursors_4 = torch.cat(
            [precursors_3, torch.full((1, 1), 0.9)], dim=1
        )

        with torch.no_grad():
            mem1, _ = im_model._encode_spectra(spectra, precursors_4)
            mem2, _ = reloaded._encode_spectra(spectra, precursors_4)

        max_diff = (mem1 - mem2).abs().max().item()
        logger.info("Max diff after reload: %.10f", max_diff)

        # fp32 accumulation through 9 transformer layers can introduce
        # small numerical differences (~1e-6) on save/load roundtrip
        if max_diff < 1e-4:
            logger.info("TEST 6 PASSED (within fp32 tolerance)")
            return True
        else:
            logger.warning("TEST 6 FAILED: weights differ after reload (%.6f)", max_diff)
            return False
    finally:
        os.unlink(tmp_path)


def test_3col_fallback():
    """Test 7: IM model works with 3-column precursors (no IM data)."""
    logger.info("=" * 60)
    logger.info("TEST 7: 3-column precursor fallback")
    logger.info("=" * 60)

    from casanovo.denovo.model import Spec2Pep

    im_model = Spec2Pep.from_pretrained(CHECKPOINT, use_ion_mobility=True)
    im_model.eval()

    spectra, precursors_3 = create_synthetic_batch(
        batch_size=2, n_peaks=30, with_im=False
    )

    # _encode_spectra should work with 3-column precursors (skip IM)
    with torch.no_grad():
        memories, masks = im_model._encode_spectra(spectra, precursors_3)

    logger.info("3-col precursors: encoder output shape %s", memories.shape)

    # Full forward should also work
    with torch.no_grad():
        results = im_model(spectra, precursors_3)

    n_pred = sum(len(r) for r in results)
    logger.info("Predictions with 3-col precursors: %d", n_pred)

    logger.info("TEST 7 PASSED")
    return True


def main():
    """Run all tests."""
    logger.info("Casanovo IM Integration Test Suite")
    logger.info("Checkpoint: %s", CHECKPOINT)
    logger.info("PyTorch: %s", torch.__version__)
    logger.info("Device: %s", "cuda" if torch.cuda.is_available() else "cpu")
    logger.info("")

    results = {}

    # Test 1: Checkpoint inspection
    results["checkpoint"] = test_checkpoint_loading()
    if not results["checkpoint"]:
        logger.error("Cannot proceed without checkpoint")
        return

    # Test 2: Model without IM
    orig_result = test_original_model()
    if orig_result is None:
        logger.error("Model without IM failed to load")
        results["original"] = False
    else:
        results["original"] = True
        orig_memories, orig_masks, spectra, precursors = orig_result

        # Test 3: Zero-init backward compatibility
        results["zero_init"] = test_im_model_zero_init(
            orig_memories, orig_masks, spectra, precursors
        )

    # Test 4: Nonzero IM values
    results["nonzero_im"] = test_im_model_nonzero()

    # Test 5: Full forward pass
    results["full_forward"] = test_full_forward()

    # Test 6: Save/load roundtrip
    results["roundtrip"] = test_save_load_roundtrip()

    # Test 7: 3-column fallback
    results["3col_fallback"] = test_3col_fallback()

    # Summary
    logger.info("")
    logger.info("=" * 60)
    logger.info("TEST SUMMARY")
    logger.info("=" * 60)
    all_passed = True
    for name, passed in results.items():
        status = "PASS" if passed else "FAIL"
        logger.info("  %-20s %s", name, status)
        if not passed:
            all_passed = False

    if all_passed:
        logger.info("\nAll tests passed!")
    else:
        logger.info("\nSome tests failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()
