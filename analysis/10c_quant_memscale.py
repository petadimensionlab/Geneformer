#!/usr/bin/env python
"""Memory/throughput scaling vs batch size for Geneformer, to size GPU budgets.

Answers: "which (model, dtype, batch, seq_len) combinations fit in N GiB of GPU
memory, and at what throughput". Peak is measured per configuration with
torch.cuda.reset_peak_memory_stats right before the forward pass.

Usage:
    .venv/bin/python analysis/10c_quant_memscale.py
Environment:
    GENEFORMER_DIR / GENEFORMER_MODEL (default Geneformer-V2-104M)
    QS_BATCHES (default "8,32,64"), QS_LEN (default 4096), QS_ITERS (default 2)
    QS_VARIANTS (default "fp32,bf16,int8,nf4"), QS_LIMIT_GIB (skip configs that
      would exceed this estimated peak; default 60)
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
BATCHES = [int(x) for x in os.environ.get("QS_BATCHES", "8,32,64").split(",")]
LEN = int(os.environ.get("QS_LEN", 4096))
ITERS = int(os.environ.get("QS_ITERS", 2))
WANT = os.environ.get("QS_VARIANTS", "fp32,bf16,int8,nf4").split(",")
LIMIT = float(os.environ.get("QS_LIMIT_GIB", 60)) * 2**30
VOCAB = 20275

VARIANTS = {
    "fp32": ({}, None),
    "bf16": ({"torch_dtype": torch.bfloat16}, None),
    "int8": ({}, BitsAndBytesConfig(load_in_8bit=True)),
    "nf4": (
        {},
        BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_use_double_quant=True,
            bnb_4bit_compute_dtype=torch.bfloat16,
        ),
    ),
}


def load(tag):
    load_kw, qc = VARIANTS[tag]
    kw = dict(load_kw)
    if qc is not None:
        kw["quantization_config"] = qc
        kw["device_map"] = {"": 0}
    model = BertForMaskedLM.from_pretrained(str(MODEL_DIR), output_hidden_states=True, **kw)
    if qc is None:
        model = model.cuda()
    return model.eval()


def main() -> None:
    assert torch.cuda.is_available(), "CUDA unavailable"
    n_params = sum(p.numel() for p in load("fp32").parameters())
    torch.cuda.empty_cache()
    print(f"model={MODEL_DIR.name} params={n_params/1e6:.1f}M len={LEN} iters={ITERS}")
    print(f"{'variant':8s} {'batch':>5s} {'ms/iter':>9s} {'tok/s':>9s} {'cells/s':>8s} {'peak_GiB':>9s}")
    torch.manual_seed(0)
    for tag in WANT:
        if tag not in VARIANTS:
            continue
        for batch in BATCHES:
            ids = torch.randint(3, VOCAB, (batch, LEN)).cuda()
            est = batch * LEN * 768 * 4 * 30  # rough activation estimate, fp32 bytes
            if est > LIMIT:
                print(f"{tag:8s} {batch:5d}  skipped (estimated peak > {LIMIT/2**30:.0f} GiB)")
                del ids
                continue
            try:
                model = load(tag)
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
                print(f"{tag:8s} {batch:5d} {dt*1000:9.1f} {batch*LEN/dt:9.0f} "
                      f"{batch/dt:8.1f} {peak:9.2f}")
            except torch.cuda.OutOfMemoryError as exc:
                print(f"{tag:8s} {batch:5d}  CUDA OOM: {str(exc).splitlines()[0][:80]}")
            finally:
                try:
                    del model
                except UnboundLocalError:
                    pass
                del ids
                torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
