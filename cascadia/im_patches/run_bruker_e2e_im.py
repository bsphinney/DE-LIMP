"""End-to-end Bruker .d -> Cascadia de novo sequencing with IM-enhanced spectra.

Uses SpectrumReader.new_with_span_step() to produce IM-sliced augmented spectra.
Each MS2 peak carries a real 1/K0 value from the IM slice midpoint.

Usage:
    python run_bruker_e2e_im.py <d_folder> <model_path> <outdir> [score_threshold] [batch_size]

The model must support use_ion_mobility=True (5-column input).
For standard 4-column models, use the original run_bruker_e2e.py instead.
"""
import sys
import os
import logging
import time

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("bruker_e2e_im")

d_folder = sys.argv[1]
model_path = sys.argv[2]
outdir = sys.argv[3]
score_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 0.5
batch_size = int(sys.argv[5]) if len(sys.argv) > 5 else 64
scan_width = 1
max_charge = 3
im_span = 0.05
im_step = 0.02

# Import IM-enhanced Cascadia components
from bruker_patch_im import (
    get_centers_bruker_im,
    extract_spectra_bruker_im,
    write_asf_im,
)

# Step 1: Get centers with IM slicing
logger.info("Step 1: Get IM-sliced centers from .d file: %s", os.path.basename(d_folder))
t0 = time.time()
f_to_mzrt_to_pep, max_mz, window_size, cycle_time = get_centers_bruker_im(
    d_folder, im_span=im_span, im_step=im_step
)
logger.info("Centers done in %.1fs", time.time() - t0)

# Step 2: Extract spectra with IM and write 5-column ASF
time_width = (scan_width + 1) * cycle_time
logger.info("Step 2: Extract IM-enhanced spectra (time_width=%.2fs, cycle_time=%.2fs)",
            time_width, cycle_time)
t0 = time.time()
temp_path = os.path.join(outdir, "temp")
os.makedirs(temp_path, exist_ok=True)
outfile = os.path.join(temp_path, "temp.asf")

if os.path.exists(outfile):
    os.remove(outfile)

total_precursors = 0
for part in sorted(f_to_mzrt_to_pep.keys()):
    prec_to_spec = extract_spectra_bruker_im(
        d_folder, f_to_mzrt_to_pep, part, 150, time_width, max_mz,
        im_span=im_span, im_step=im_step
    )
    n_prec = len(prec_to_spec)
    total_precursors += n_prec
    n_augmented = n_prec * (max_charge - 1)
    logger.info("Part %d: %d precursors -> %d augmented spectra", part, n_prec, n_augmented)
    write_asf_im(outfile, prec_to_spec, scan_width=scan_width,
                 max_pep_length=30, max_charge=max_charge)

asf_size = os.path.getsize(outfile) / 1e6
logger.info("IM-enhanced ASF written in %.1fs: %.1f MB (%d total precursors)",
            time.time() - t0, asf_size, total_precursors)

# Step 3: Verify ASF format (check 5-column lines)
logger.info("Step 3: Verify ASF format")
with open(outfile) as f:
    n_4col = 0
    n_5col = 0
    n_header = 0
    for i, line in enumerate(f):
        parts = line.strip().split('\t')
        if len(parts) == 5 and '.' in parts[0]:
            n_5col += 1
        elif len(parts) == 4 and '.' in parts[0]:
            n_4col += 1
        else:
            n_header += 1
        if i > 10000:
            break
    logger.info("ASF format check (first 10k lines): %d 5-col, %d 4-col, %d header lines",
                n_5col, n_4col, n_header)
    if n_4col > 0:
        logger.warning("Found %d 4-column lines -- expected all 5-column!", n_4col)
    if n_5col > 0:
        logger.info("OK: 5-column IM-enhanced format confirmed")

# Step 4: Count IM value distribution
logger.info("Step 4: IM value distribution in ASF")
im_values = []
with open(outfile) as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 5 and '.' in parts[0]:
            im_val = float(parts[4])
            im_values.append(im_val)
            if len(im_values) >= 100000:
                break

import numpy as np
im_arr = np.array(im_values)
ms2_ims = im_arr[im_arr > 0]
ms1_ims = im_arr[im_arr == 0]
logger.info("IM stats (first 100k peaks): %d MS2 (IM>0, range %.3f-%.3f), %d MS1 (IM=0)",
            len(ms2_ims), ms2_ims.min() if len(ms2_ims) > 0 else 0,
            ms2_ims.max() if len(ms2_ims) > 0 else 0, len(ms1_ims))

# Step 5: Run Cascadia inference
# NOTE: This requires a model trained with use_ion_mobility=True (5-column input).
# The standard Cascadia model expects 4-column ASF.
# For now, log a warning and skip inference if model is 4-column.
logger.info("Step 5: Load model and check IM support")

from cascadia.depthcharge.data.spectrum_datasets import AnnotatedSpectrumDataset
from cascadia.depthcharge.data.preprocessing import scale_to_unit_norm, scale_intensity
from cascadia.depthcharge.tokenizers import PeptideTokenizer
from cascadia.model import AugmentedSpec2Pep
import torch

tokenizer = PeptideTokenizer.from_massivekb(
    reverse=False, replace_isoleucine_with_leucine=True
)

# Check if model supports IM
checkpoint = torch.load(model_path, map_location='cpu')
hparams = checkpoint.get('hyper_parameters', {})
has_im = hparams.get('use_ion_mobility', False)

if not has_im:
    logger.warning("Model does not have use_ion_mobility=True in hyperparameters.")
    logger.warning("Standard 4-column model cannot process 5-column IM-enhanced ASF.")
    logger.warning("To use IM-enhanced inference, train a model with the IM-aware encoder.")
    logger.info("Skipping inference. ASF file is at: %s", outfile)
    logger.info("Pipeline complete (augmentation only). "
                "Train an IM-aware model, then run inference on this ASF.")
else:
    logger.info("Model supports IM! Running inference...")
    index_path = os.path.join(temp_path, "index.hdf5")

    train_dataset = AnnotatedSpectrumDataset(
        tokenizer, outfile, index_path=index_path,
        preprocessing_fn=[scale_intensity(scaling="root"), scale_to_unit_norm]
    )
    train_loader = train_dataset.loader(
        batch_size=batch_size, num_workers=4, pin_memory=True
    )
    logger.info("Dataset: %d spectra loaded (batch_size=%d)", len(train_dataset), batch_size)

    model = AugmentedSpec2Pep.load_from_checkpoint(
        model_path,
        d_model=512, n_layers=9, n_head=8, dim_feedforward=1024,
        dropout=0, rt_width=2, tokenizer=tokenizer, max_charge=10,
        residues="massivekb", max_length=31,
        use_ion_mobility=True
    )

    import pytorch_lightning as pl
    if torch.cuda.is_available():
        logger.info("GPU found: %s", torch.cuda.get_device_name(0))
        device = 'gpu'
    else:
        logger.info("No GPU -- running on CPU (will be slow)")
        device = 'cpu'

    trainer = pl.Trainer(logger=False, enable_progress_bar=True,
                         accelerator=device, devices=1)
    predictions = trainer.predict(model, train_loader)

    # Write SSL output
    logger.info("Step 6: Write SSL results (threshold=%.2f)", score_threshold)
    from cascadia.utils import write_results

    ssl_base = os.path.join(outdir, "bruker_im_test")
    write_results(predictions, ssl_base, os.path.basename(d_folder),
                  window_size, score_threshold, scan_width * cycle_time)

    ssl_file = ssl_base + ".ssl"
    if os.path.exists(ssl_file):
        count = sum(1 for _ in open(ssl_file)) - 1
        logger.info("Done! Wrote %d peptides to %s", count, ssl_file)
    else:
        logger.warning("No SSL file produced!")

# Cleanup temp
import shutil
shutil.rmtree(temp_path, ignore_errors=True)
logger.info("Pipeline complete.")
