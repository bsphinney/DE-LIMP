#!/bin/bash
# ============================================================
# Cascadia De Novo Sequencing — HIVE HPC Setup Script
# ============================================================
# Run this on HIVE (login node):
#   bash setup_cascadia.sh
#
# What it does:
#   1. Creates a conda environment with Cascadia + timsrust_pyo3
#   2. Downloads the Cascadia model checkpoint
#   3. Verifies GPU access
#   4. Runs a quick sanity check
#
# Prerequisites:
#   - conda/mamba available (module load conda or miniconda installed)
#   - Internet access for pip/conda downloads
#   - GPU partition accessible (for testing)
# ============================================================

set -e

CASCADIA_ENV="cascadia"
CASCADIA_DIR="$HOME/cascadia"
MODEL_DIR="$CASCADIA_DIR/models"

echo "=== Cascadia Setup for HIVE ==="
echo "Install dir: $CASCADIA_DIR"
echo ""

# --- Step 1: Create conda environment ---
echo "[1/5] Creating conda environment '$CASCADIA_ENV'..."

if conda env list 2>/dev/null | grep -q "^${CASCADIA_ENV} "; then
    echo "  Environment '$CASCADIA_ENV' already exists. Updating..."
    conda activate $CASCADIA_ENV 2>/dev/null || source activate $CASCADIA_ENV
else
    conda create -n $CASCADIA_ENV python=3.10 -y -q
    conda activate $CASCADIA_ENV 2>/dev/null || source activate $CASCADIA_ENV
fi

# --- Step 2: Install Cascadia + dependencies ---
echo "[2/5] Installing Cascadia and dependencies..."

# PyTorch with CUDA (check what CUDA version is available on HIVE)
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu121 -q 2>/dev/null || \
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu118 -q 2>/dev/null || \
pip install torch torchvision -q

# Cascadia from GitHub
pip install git+https://github.com/Noble-Lab/cascadia.git -q

# Bruker native loader (Phase 1)
pip install timsrust_pyo3 -q

# DIAMOND for novel peptide BLAST (Phase 4)
if ! command -v diamond &>/dev/null; then
    echo "  Installing DIAMOND..."
    conda install -c bioconda diamond -y -q 2>/dev/null || \
    echo "  WARNING: Could not install DIAMOND via conda. Try: module load diamond"
fi

echo "  Installed packages:"
pip show cascadia 2>/dev/null | grep -E "^(Name|Version)" || echo "  cascadia: not found (check install)"
python -c "import timsrust_pyo3; print('  timsrust_pyo3:', timsrust_pyo3.__version__)" 2>/dev/null || echo "  timsrust_pyo3: not found"
python -c "import torch; print('  torch:', torch.__version__, '(CUDA:', torch.cuda.is_available(), ')')" 2>/dev/null || echo "  torch: not found"

# --- Step 3: Download model checkpoint ---
echo ""
echo "[3/5] Downloading Cascadia model checkpoint..."

mkdir -p "$MODEL_DIR"

# Cascadia model from Zenodo (Noble Lab)
# Check if already downloaded
CKPT_FILE="$MODEL_DIR/cascadia_default.ckpt"
if [ -f "$CKPT_FILE" ]; then
    echo "  Model already exists: $CKPT_FILE"
else
    echo "  Downloading from Zenodo..."
    # Try cascadia's built-in download first
    python -c "
from cascadia.utils import download_model
download_model('$MODEL_DIR')
print('  Download complete')
" 2>/dev/null || {
        echo "  Built-in download failed. Trying manual download..."
        # Manual Zenodo download (update URL if needed)
        wget -q -O "$CKPT_FILE" \
            "https://zenodo.org/records/10837080/files/cascadia_default.ckpt" 2>/dev/null || \
        curl -sL -o "$CKPT_FILE" \
            "https://zenodo.org/records/10837080/files/cascadia_default.ckpt" 2>/dev/null || \
        echo "  WARNING: Could not download model. Check Cascadia docs for latest URL."
    }
fi

if [ -f "$CKPT_FILE" ]; then
    echo "  Model checkpoint: $CKPT_FILE ($(du -h "$CKPT_FILE" | cut -f1))"
else
    echo "  WARNING: No model checkpoint found. Cascadia won't run without it."
fi

# --- Step 4: Verify installation ---
echo ""
echo "[4/5] Verifying installation..."

python -c "
import sys
errors = []

try:
    import cascadia
    print('  cascadia: OK')
except ImportError as e:
    errors.append(f'cascadia: {e}')
    print(f'  cascadia: FAILED ({e})')

try:
    import timsrust_pyo3
    print('  timsrust_pyo3: OK')
except ImportError as e:
    errors.append(f'timsrust_pyo3: {e}')
    print(f'  timsrust_pyo3: FAILED ({e})')

try:
    import torch
    cuda = torch.cuda.is_available()
    print(f'  torch: OK (CUDA available: {cuda})')
    if not cuda:
        errors.append('CUDA not available on login node (expected — test on GPU node)')
except ImportError as e:
    errors.append(f'torch: {e}')
    print(f'  torch: FAILED ({e})')

if errors:
    print(f'\\n  {len(errors)} issue(s) found — see above')
else:
    print('\\n  All checks passed!')
"

# --- Step 5: Test command ---
echo ""
echo "[5/5] Testing cascadia CLI..."

cascadia --help 2>/dev/null | head -3 || \
python -m cascadia --help 2>/dev/null | head -3 || \
echo "  WARNING: cascadia CLI not found. May need: python -m cascadia sequence ..."

# --- Summary ---
echo ""
echo "=== Setup Complete ==="
echo ""
echo "To use Cascadia on HIVE:"
echo "  conda activate $CASCADIA_ENV"
echo "  cascadia sequence --input sample.d --model $CKPT_FILE --output results.ssl"
echo ""
echo "To test on a GPU node:"
echo "  srun --partition=gpu --gres=gpu:1 --time=00:30:00 --pty bash"
echo "  conda activate $CASCADIA_ENV"
echo "  python -c \"import torch; print(torch.cuda.is_available())\"  # should be True"
echo ""
echo "Model checkpoint: $CKPT_FILE"
echo "Conda environment: $CASCADIA_ENV"
echo ""
echo "Next steps:"
echo "  1. Test on GPU node with a single .d file"
echo "  2. If successful, enable Cascadia in DE-LIMP (feature/cascadia-denovo branch)"
echo ""
