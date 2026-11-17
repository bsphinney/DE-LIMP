"""Teacher-forced validation metrics (CE loss, AA and peptide accuracy) of checkpoints on one ASF.

Answers "did fine-tuning improve on the stock model for these spectra?", which the training log
can't: it starts at the first fine-tuned validation (step 10000) and has no stock baseline.

    python val_only.py VAL.asf name=ckpt [name=ckpt ...]
"""
import hashlib
import os
import sys
import tempfile

import cascadia.depthcharge.data.spectrum_datasets as sd


def _describe(o):          # same picklable-free index hash as cascadia_train_wrap.py
    if isinstance(o, (list, tuple)):
        return [_describe(x) for x in o]
    vals = []
    for c in getattr(o, "__closure__", None) or ():
        try:
            vals.append(repr(c.cell_contents))
        except ValueError:
            vals.append("<empty>")
    return [getattr(o, "__module__", ""), getattr(o, "__qualname__", repr(o)), vals]


sd._hash_obj = lambda obj: hashlib.sha1(repr(_describe(obj)).encode()).hexdigest()

import pytorch_lightning as pl  # noqa: E402
import torch  # noqa: E402
from cascadia.depthcharge.data.preprocessing import scale_intensity, scale_to_unit_norm  # noqa: E402
from cascadia.depthcharge.data.spectrum_datasets import AnnotatedSpectrumDataset  # noqa: E402
from cascadia.depthcharge.tokenizers import PeptideTokenizer  # noqa: E402
from cascadia.model import AugmentedSpec2Pep  # noqa: E402

asf = sys.argv[1]
models = dict(a.split("=", 1) for a in sys.argv[2:])
tok = PeptideTokenizer.from_massivekb(reverse=False, replace_isoleucine_with_leucine=True)
tmp = tempfile.mkdtemp(prefix="valonly_", dir=os.environ.get("TMPDIR", "/tmp"))
ds = AnnotatedSpectrumDataset(tok, asf, index_path=os.path.join(tmp, "val.hdf5"),
                              preprocessing_fn=[scale_intensity(scaling="root"), scale_to_unit_norm])
loader = ds.loader(batch_size=32, num_workers=4, pin_memory=True)
for name, ckpt in models.items():
    model = AugmentedSpec2Pep.load_from_checkpoint(ckpt, d_model=512, n_layers=9, n_head=8, dim_feedforward=1024,
                                                   dropout=0, rt_width=2, tokenizer=tok, max_charge=10)
    tr = pl.Trainer(accelerator="gpu" if torch.cuda.is_available() else "cpu", devices=1,
                    enable_progress_bar=False, logger=False)
    out = tr.validate(model, dataloaders=loader, verbose=False)[0]
    print(f"[{name}] " + "  ".join(f"{k}={v:.4f}" for k, v in sorted(out.items()) if "Frag" not in k), flush=True)
