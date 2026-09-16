# DirectML Backend Setup for Geneformer

This document describes how to configure Geneformer to use the **DirectML** backend for GPU acceleration on Windows / WSL2 systems without an NVIDIA GPU. DirectML works with AMD, Intel, and Qualcomm GPUs via the DirectX 12 API.

## Overview

Geneformer's `device.py` module automatically selects the best available compute backend:

1. **DirectML (dml)** — Windows / WSL2, any DirectX 12 GPU
2. **CUDA (cuda)** — NVIDIA GPU or AMD ROCm
3. **MPS (mps)** — Apple Silicon
4. **CPU (cpu)** — Fallback

You can force a specific backend with the `GENEFORMER_DEVICE` environment variable.

## Prerequisites

- **Windows 10/11** with WSL2 (Ubuntu 22.04+ recommended)
- **DirectX 12-compatible GPU** (AMD, Intel, NVIDIA, or Qualcomm)
- **WSL GPU driver** installed on the Windows host

## Installation

### 1. Set up the Python environment

```bash
cd /path/to/Geneformer
python3 -m venv .venv
source .venv/bin/activate
```

### 2. Install PyTorch (CPU-only) and torch-directml

```bash
# Install CPU-only PyTorch (required base for torch-directml)
pip install torch==2.4.1 torchvision==0.19.1 --index-url https://download.pytorch.org/whl/cpu

# Install torch-directml matching version
pip install torch-directml
```

> **Note:** torch-directml 0.2.5.dev240914 requires `torch==2.4.1` and `torchvision==0.19.1`. Check the latest compatible versions on [PyPI](https://pypi.org/project/torch-directml/).

### 3. Download Geneformer code and apply patches

```bash
# Download the codebase and models
./download.sh

# Apply multi-backend patches + DML-aware device.py
./patches/apply_patches.sh
```

### 4. Verify the setup

```python
from geneformer.device import get_device, get_device_obj

print(f"Detected device: {get_device()}")        # Should print: dml
print(f"Device object: {get_device_obj()}")       # Should print: privateuseone:0
```

## Usage

### Running with DirectML (automatic)

By default, Geneformer will detect DirectML and use it:

```bash
source .venv/bin/activate
python analysis/04_baseline.py
```

### Forcing a specific backend

```bash
# Force CPU
GENEFORMER_DEVICE=cpu python analysis/04_baseline.py

# Force DirectML
GENEFORMER_DEVICE=dml python analysis/04_baseline.py
```

## How It Works

### Device detection chain (`geneformer/device.py`)

1. `get_device()` checks if `torch-directml` is available and returns `"dml"`.
2. Otherwise, checks `torch.cuda.is_available()` (NVIDIA CUDA / AMD ROCm).
3. Otherwise, checks Apple's MPS.
4. Falls back to `"cpu"`.

### Patch overview (`patches/geneformer_multibackend.patch`)

The patch replaces hardcoded `"cuda"` device strings across 9 files with calls to `get_device_obj()`, `move_to_device()`, `empty_cache()`, and `manual_seed_all()` from `device.py`.

Specific fixes include:

- `device="cuda"` → `device=get_device_obj()`
- `.to("cuda")` → `.to(get_device_obj())`
- `torch.cuda.empty_cache()` → `empty_cache()`
- `torch.device("cuda" if ... else "cpu")` → `get_device_obj()`
- `df["col"][bool_index]` → `df["col"].iloc[bool_index]` (FutureWarning fix in `tokenizer.py`)

## Troubleshooting

| Symptom | Likely Cause | Solution |
|---------|-------------|----------|
| `RuntimeError: Found no NVIDIA driver` | No CUDA, and DirectML not installed | Install `torch-directml` or set `GENEFORMER_DEVICE=cpu` |
| `ImportError: undefined symbol` | torch-directml version mismatch with PyTorch | Install compatible versions (see Installation step 2) |
| `DML device not found` | WSL GPU driver not installed | Install [GPU driver for WSL](https://learn.microsoft.com/en-us/windows/wsl/tutorials/gpu-compute) on Windows host |
| `privateuseone:0` is slow | First-time JIT compilation; performance improves after warmup | Run a small forward pass first to warm up the GPU |
| `Failed to load CPU gemm_4bit_forward` | bitsandbytes CPU fallback warning | Ignore — it's harmless and falls back to pure PyTorch |

## Version Compatibility

| Component | Version |
|-----------|---------|
| Python | 3.10 – 3.12 |
| PyTorch | 2.4.1 (CPU) |
| torch-directml | 0.2.5.dev240914 |
| torchvision | 0.19.1 (CPU) |
| OS | Ubuntu 22.04+ on WSL2 (Windows 10/11) |