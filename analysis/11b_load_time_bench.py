#!/usr/bin/env python
"""Fair load-time benchmark: how long does each way of getting a model onto the
GPU actually take, and what does it cost?

Timing covers the WHOLE operation (from_pretrained plus the device move), which
is the number that matters when a pipeline reloads the model per process.

Variants:
  fp32 (cuda())            fp32 checkpoint, from_pretrained then .cuda()
  fp32 devmap              fp32 checkpoint with device_map={"":0}
  bf16 devmap              fp32 checkpoint, cast to bf16 and placed via device_map
  int8 load-time           geneformer's model_type="Pretrained-Quantized" path
  nf4 load-time            4-bit quantized at load
  bf16_pre                 checkpoint pre-saved in bf16 (see analysis/11_...)
  bf16_pre + dtype arg     same, loading with an explicit torch_dtype
  int8_pre / nf4_pre       checkpoint pre-quantized by bitsandbytes

Usage:
    .venv/bin/python analysis/11b_load_time_bench.py
Environment: GENEFORMER_DIR / GENEFORMER_MODEL, QP_QUANT_DIR (default quantized/),
    QL_REPS (default 5)
"""
from __future__ import annotations

import gc
import json
import os
import statistics
import time
from pathlib import Path

import torch
from transformers import BertForMaskedLM, BitsAndBytesConfig

ROOT = Path(__file__).resolve().parent.parent
MODEL_DIR = Path(os.environ.get("GENEFORMER_DIR", ROOT / "geneformer_hf")) / os.environ.get(
    "GENEFORMER_MODEL", "Geneformer-V2-104M"
)
QUANT_ROOT = Path(os.environ.get("QP_QUANT_DIR", ROOT / "quantized"))
REPS = int(os.environ.get("QL_REPS", 5))
NAME = MODEL_DIR.name


def int8() -> BitsAndBytesConfig:
    return BitsAndBytesConfig(load_in_8bit=True)


def nf4() -> BitsAndBytesConfig:
    return BitsAndBytesConfig(
        load_in_4bit=True,
        bnb_4bit_quant_type="nf4",
        bnb_4bit_use_double_quant=True,
        bnb_4bit_compute_dtype=torch.bfloat16,
    )


def timed_load(path, **kw) -> tuple[float, float, str]:
    """Return (seconds for the whole load, weights MiB, weights dtype)."""
    gc.collect()
    torch.cuda.empty_cache()
    torch.cuda.synchronize()
    t0 = time.time()
    model = BertForMaskedLM.from_pretrained(str(path), output_hidden_states=True, **kw)
    quantized = "quantization_config" in kw or _is_quantized_dir(Path(path))
    if not quantized and "device_map" not in kw:
        model = model.cuda()
    torch.cuda.synchronize()
    dt = time.time() - t0
    weights = sum(p.numel() * p.element_size() for p in model.parameters()) / 2**20
    dtype = str(next(model.parameters()).dtype)
    del model
    return dt, weights, dtype


def _is_quantized_dir(path: Path) -> bool:
    cfg = path / "config.json"
    if not cfg.exists():
        return False
    try:
        return "quantization_config" in json.loads(cfg.read_text())
    except json.JSONDecodeError:
        return False


def main() -> None:
    assert torch.cuda.is_available(), "CUDA unavailable"
    pre = {tag: QUANT_ROOT / f"{NAME}-{tag}" for tag in ("int8", "nf4", "bf16")}
    variants = [
        ("fp32 (cuda())", MODEL_DIR, {}),
        ("fp32 devmap", MODEL_DIR, {"device_map": {"": 0}}),
        ("bf16 devmap", MODEL_DIR, {"device_map": {"": 0}, "torch_dtype": torch.bfloat16}),
        ("int8 load-time", MODEL_DIR, {"quantization_config": int8(), "device_map": {"": 0}}),
        ("nf4 load-time", MODEL_DIR, {"quantization_config": nf4(), "device_map": {"": 0}}),
    ]
    if (pre["bf16"] / "model.safetensors").exists():
        variants += [
            ("bf16_pre", pre["bf16"], {}),
            ("bf16_pre +dtype", pre["bf16"],
             {"device_map": {"": 0}, "torch_dtype": torch.bfloat16}),
        ]
    if (pre["int8"] / "model.safetensors").exists():
        variants.append(("int8_pre", pre["int8"], {}))
    if (pre["nf4"] / "model.safetensors").exists():
        variants.append(("nf4_pre", pre["nf4"], {}))

    _ = torch.zeros(1024, device="cuda")
    for _, p, kw in variants:  # warm up page cache, CUDA context and bnb kernels
        m = BertForMaskedLM.from_pretrained(str(p), **kw)
        del m
    gc.collect()
    torch.cuda.empty_cache()

    print(f"model={NAME} reps={REPS}")
    print(f"{'variant':18s} {'median_s':>9s} {'min_s':>7s} {'weights_MiB':>12s} {'dtype':>16s}")
    rows = []
    for name, path, kw in variants:
        times, w, dtype = [], None, ""
        for _ in range(REPS):
            dt, w, dtype = timed_load(path, **kw)
            times.append(dt)
        med = statistics.median(times)
        rows.append({"variant": name, "median_s": med, "min_s": min(times),
                     "weights_mib": w, "dtype": dtype})
        print(f"{name:18s} {med:9.3f} {min(times):7.3f} {w:12.1f} {dtype:>16s}")

    out = QUANT_ROOT / f"loadtime_bench_{NAME}.json"
    out.write_text(json.dumps(rows, indent=2))
    print(f"\nresults -> {out}")


if __name__ == "__main__":
    main()
