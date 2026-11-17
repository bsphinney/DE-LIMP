"""Run `cascadia train ...` with a picklable-free index hash.

depthcharge's _hash_obj dill-pickles the preprocessing functions to fingerprint the HDF5 index.
scale_intensity("root") is a closure, so dill pickles it by value, reaches numpy's
_ArrayFunctionDispatcher and fails under numpy >= 1.25 ("Can't pickle <class
'numpy._ArrayFunctionDispatcher'>", 2026-11-16). The hash only checks that a cached index was
built with the same preprocessing, so the functions' qualified names (and closure cell values)
are an adequate fingerprint.

    python3 cascadia_train_wrap.py train train.asf val.asf -m stock.ckpt ...
"""
import hashlib
import sys

import cascadia.depthcharge.data.spectrum_datasets as sd


def _describe(o):
    if isinstance(o, (list, tuple)):
        return [_describe(x) for x in o]
    cells = getattr(o, "__closure__", None) or ()
    vals = []
    for c in cells:
        try:
            vals.append(repr(c.cell_contents))
        except ValueError:
            vals.append("<empty>")
    return [getattr(o, "__module__", ""), getattr(o, "__qualname__", repr(o)), vals]


def _hash_obj(obj):
    return hashlib.sha1(repr(_describe(obj)).encode()).hexdigest()


sd._hash_obj = _hash_obj

from cascadia.cascadia import main  # noqa: E402

sys.argv = ["cascadia"] + sys.argv[1:]
sys.exit(main())
