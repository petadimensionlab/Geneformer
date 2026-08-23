"""
Geneformer device detection.

Routes compute selection across the best available accelerator with an explicit
override via the ``GENEFORMER_DEVICE`` environment variable.

Supported device names:

- ``cuda`` : NVIDIA CUDA **or** AMD ROCm. ROCm exposes the CUDA API, so a ROCm
  build reports ``torch.cuda.is_available() == True`` and is driven through the
  same code path.
- ``dml``  : DirectML via ``torch-directml`` (Windows / WSL2).
- ``mps``  : Apple Silicon Metal Performance Shader (Mac MPS).
- ``cpu``  : CPU fallback.

Auto-detection priority is DML > CUDA/ROCm > MPS > CPU. Set ``GENEFORMER_DEVICE`` to
force a specific backend (e.g. ``GENEFORMER_DEVICE=mps`` on a Mac that also has
CUDA available via a remote/egpu setup).
"""
from __future__ import annotations

import os

import torch

_VALID = {"cuda", "dml", "mps", "cpu"}


_DML_CHECKED: bool = False
_DML_AVAILABLE: bool = False


def _dml_available() -> bool:
    global _DML_CHECKED, _DML_AVAILABLE
    if _DML_CHECKED:
        return _DML_AVAILABLE
    _DML_CHECKED = True
    try:
        import torch_directml  # noqa: F401

        _DML_AVAILABLE = True
    except (ImportError, OSError):
        _DML_AVAILABLE = False
    return _DML_AVAILABLE


def _detect() -> str:
    """Return the best available device name without consulting the override."""
    if _dml_available():
        return "dml"
    if torch.cuda.is_available():
        return "cuda"
    mps = getattr(torch.backends, "mps", None)
    if mps is not None and mps.is_available():
        return "mps"
    return "cpu"


def get_device() -> str:
    """Return the resolved device name (``cuda``, ``dml``, ``mps``, or ``cpu``).

    Honors the ``GENEFORMER_DEVICE`` environment variable override.
    """
    override = os.environ.get("GENEFORMER_DEVICE", "").strip().lower()
    if override:
        if override not in _VALID:
            raise ValueError(
                f"GENEFORMER_DEVICE={override!r} is invalid; "
                f"choose from {sorted(_VALID)}"
            )
        return override
    return _detect()


def get_device_obj() -> torch.device:
    """Return a ``torch.device`` for the resolved accelerator."""
    name = get_device()
    if name == "dml":
        import torch_directml

        return torch_directml.device()
    return torch.device(name)


def move_to_device(model: torch.nn.Module) -> torch.nn.Module:
    """Move ``model`` to the resolved accelerator if it is not already there."""
    device = get_device_obj()
    model_device = next(model.parameters()).device
    if model_device.type != device.type:
        model.to(device)
    return model


def empty_cache() -> None:
    """Best-effort accelerator cache clear. No-op on MPS/CPU."""
    name = get_device()
    if name == "cuda" and torch.cuda.is_available():
        torch.cuda.empty_cache()


def manual_seed_all(seed: int) -> None:
    """Seed all available generators."""
    torch.manual_seed(seed)
    name = get_device()
    if name == "cuda" and torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)