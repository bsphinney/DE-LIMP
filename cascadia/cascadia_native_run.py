"""Cascadia inference on a Bruker diaPASEF .d with Cascadia's OWN augmentation (bruker_native_augment).

Usage: cascadia_native_run.py [--top-n N] [--calib RUN.calib.json] <run.d> <tag> <outdir> <tmpdir>
                              name=ckpt [name=ckpt ...]
--top-n / --calib must match what the checkpoint was trained on (prepare_training_native.py
--variants / its RUN.calib.json); the stock checkpoint saw 150 peaks and no calibration.
Writes <outdir>/<name>/<tag>.ssl with every prediction (score threshold 0; filter later).
Same model construction and preprocessing as cascadia.cascadia.sequence(); the ASF and HDF5 index
go in <tmpdir> (use node-local $TMPDIR, never Quobyte).
"""
import gc
import hashlib
import os
import sys
import time
import warnings

import numpy as np
import pytorch_lightning as pl
import torch

warnings.filterwarnings("ignore")
torch.set_float32_matmul_precision("high")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import cascadia.depthcharge.data.spectrum_datasets as _sd  # noqa: E402
from cascadia.depthcharge.data.preprocessing import scale_intensity, scale_to_unit_norm  # noqa: E402
from cascadia.depthcharge.data.spectrum_datasets import AnnotatedSpectrumDataset  # noqa: E402
from cascadia.depthcharge.tokenizers import PeptideTokenizer  # noqa: E402
from cascadia.model import AugmentedSpec2Pep  # noqa: E402
from cascadia.utils import write_results  # noqa: E402

from bruker_native_augment import augment_spectra_bruker_native  # noqa: E402

# depthcharge labels its HDF5 index with a dill hash that fails on numpy internals here; label by name.
_sd._hash_obj = lambda obj: hashlib.sha1(
    repr([getattr(f, "__qualname__", repr(f)) for f in obj]).encode()).hexdigest()

MAX_CHARGE, BATCH, WIDTH = 4, 32, 2   # cascadia sequence defaults (-c 4, -b 32, -w 2)

import argparse as _ap  # noqa: E402
import json as _json  # noqa: E402
_p = _ap.ArgumentParser()
_p.add_argument("--top-n", type=int, default=150)
_p.add_argument("--calib", default=None)
_p.add_argument("rest", nargs="+")
_args = _p.parse_args()
dfile, tag, outdir, tmp = _args.rest[:4]
models = dict(a.split("=", 1) for a in _args.rest[4:])
CALIB = _json.load(open(_args.calib)) if _args.calib else None
print("top_n", _args.top_n, "calib", {k: CALIB[k]["a"] for k in ("ms2", "ms1")} if CALIB else None, flush=True)
os.makedirs(tmp, exist_ok=True)
print("numpy", np.__version__, "torch", torch.__version__, "CUDA", torch.cuda.is_available(),
      torch.cuda.get_device_name(0) if torch.cuda.is_available() else "", flush=True)
tok = PeptideTokenizer.from_massivekb(reverse=False, replace_isoleucine_with_leucine=True)

t0 = time.time()
asf, window_size, cycle_time = augment_spectra_bruker_native(dfile, tmp, top_n=_args.top_n, max_charge=MAX_CHARGE,
                                                               calib=CALIB)
n_rec = sum(1 for line in open(asf) if line.startswith("BEGIN IONS"))
print(f"native augment: {n_rec} records, {os.path.getsize(asf)/1e9:.2f} GB, window {window_size} m/z, "
      f"cycle {cycle_time:.3f} s, {time.time()-t0:.0f}s", flush=True)

for name, ckpt in models.items():
    t0 = time.time()
    os.makedirs(os.path.join(outdir, name), exist_ok=True)
    base = os.path.join(outdir, name, tag)
    ds = AnnotatedSpectrumDataset(tok, asf, index_path=os.path.join(tmp, f"index_{name}.hdf5"),
                                  preprocessing_fn=[scale_intensity(scaling="root"), scale_to_unit_norm])
    loader = ds.loader(batch_size=BATCH, num_workers=0, pin_memory=True)
    model = AugmentedSpec2Pep.load_from_checkpoint(ckpt, d_model=512, n_layers=9, n_head=8, dim_feedforward=1024,
                                                   dropout=0, rt_width=2, tokenizer=tok, max_charge=10)
    tr = pl.Trainer(accelerator="gpu" if torch.cuda.is_available() else "cpu", devices=1,
                    enable_progress_bar=False, logger=False)
    preds = tr.predict(model, dataloaders=loader)
    write_results(preds, base, dfile, window_size, 0.0, WIDTH * cycle_time)
    n = sum(1 for _ in open(base + ".ssl")) - 1 if os.path.exists(base + ".ssl") else 0
    print(f"[{name}] {n} PSMs -> {base}.ssl ({time.time()-t0:.0f}s)", flush=True)
    del model, ds, loader, preds
    gc.collect()
    torch.cuda.empty_cache()
    if n == 0:
        sys.exit(f"[{name}] no output")
