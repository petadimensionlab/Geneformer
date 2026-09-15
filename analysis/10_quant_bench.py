#!/usr/bin/env python
"""Quantization throughput / memory benchmark for Geneformer on the local GPU.

Measures fp32 vs bf16 vs bitsandbytes int8 vs bitsandbytes nf4 for a fixed
synthetic batch, across sequence lengths, and reports peak allocated memory.

Usage:
    .venv/bin/python analysis/10_quant_bench.py
Environment:
    GENEFORMER_DIR (default: geneformer_hf)
    GENEFORMER_MODEL (default: Geneformer-V2-104M)
    QB_BATCH (default: 8), QB_ITERS (default: 3)
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import torch
from transformers import BertForMaskedLM, BitsAndBytesConfig

ROOT = Path(__file__).resolve().parent.parent
MODEL_DIR = Path(os.environ.get("GENEFORMER_DIR", ROOT / "geneformer_hf")) / os.environ.get(
    "GENEFORMER_MODEL", "Geneformer-V2-104M"
)
B = int(os.environ.get("QB_BATCH", 8))
ITERS = int(os.environ.get("QB_ITERS", 3))
LENGTHS = (512, 2048, 4096)


def variants():
    yield "fp32", {}, None
    yield "bf16", {"torch_dtype": torch.bfloat16}, None
    yield "int8", {}, BitsAndBytesConfig(load_in_8bit=True)
    yield (
        "nf4",
        {},
        BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_use_double_quant=True,
            bnb_4bit_compute_dtype=torch.bfloat16,
        ),
    )


def main() -> None:
    assert torch.cuda.is_available(), "CUDA unavailable"
    assert MODEL_DIR.is_dir(), f"model dir not found: {MODEL_DIR}"
    print(f"model={MODEL_DIR.name} device={torch.cuda.get_device_name(0)} "
          f"cap={torch.cuda.get_device_capability(0)} batch={B}")
    torch.manual_seed(0)
    vocab = 20275

    for length in LENGTHS:
        ids = torch.randint(3, vocab, (B, length)).cuda()
        for tag, load_kw, qc in variants():
            kw = dict(load_kw)
            if qc is not None:
                kw["quantization_config"] = qc
                kw["device_map"] = {"": 0}
            try:
                model = BertForMaskedLM.from_pretrained(
                    str(MODEL_DIR), output_hidden_states=True, **kw
                )
                if qc is None:
                    model = model.cuda()
                model.eval()
                torch.cuda.reset_peak_memory_stats()
                with torch.no_grad():
                    model(ids)
                    torch.cuda.synchronize()
                    t0 = time.time()
                    for _ in range(ITERS):
                        model(ids)
                    torch.cuda.synchronize()
                    dt = (time.time() - t0) / ITERS
                peak = torch.cuda.max_memory_allocated() / 2**30
                print(f"L={length:5d} {tag:5s} {dt * 1000:8.1f} ms/iter "
                      f"tok/s={B * length / dt:10.0f}  peak={peak:.2f} GiB")
            except Exception as exc:  # keep benchmarking the remaining variants
                print(f"L={length:5d} {tag:5s} FAILED {type(exc).__name__}: {exc}")
            finally:
                try:
                    del model
                except UnboundLocalError:
                    pass
                torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
