#!/usr/bin/env python
"""Real-workload GPU memory probe: actual tokenized cells at the configured
forward batch size, emulating EmbExtractor / InSilicoPerturber batching.

Padding is per forward batch (to the longest cell in the batch), which is what
the collator does, so this is the number that must fit the GPU budget.

Usage:
    .venv/bin/python analysis/10d_real_workload_mem.py
Environment:
    GENEFORMER_DIR / GENEFORMER_MODEL, ADPD_TISSUE (default AD_blood),
    QW_FBS (default 64), QW_NCELLS (default 512), QW_VARIANTS
"""
from __future__ import annotations

import os
import time
from pathlib import Path

import numpy as np
import torch
from datasets import load_from_disk
from transformers import BertForMaskedLM, BitsAndBytesConfig

ROOT = Path(__file__).resolve().parent.parent
MODEL_DIR = Path(os.environ.get("GENEFORMER_DIR", ROOT / "geneformer_hf")) / os.environ.get(
    "GENEFORMER_MODEL", "Geneformer-V2-104M"
)
TISSUE = os.environ.get("ADPD_TISSUE", "AD_blood")
DATA = ROOT / "input" / TISSUE / "tokenized" / f"{TISSUE}.dataset"
FBS = int(os.environ.get("QW_FBS", 64))
NCELLS = int(os.environ.get("QW_NCELLS", 512))
WANT = os.environ.get("QW_VARIANTS", "fp32,bf16,int8,nf4").split(",")

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


def real_batches():
    ds = load_from_disk(str(DATA))
    ds = ds.select(range(min(NCELLS, len(ds))))
    seqs = [np.asarray(x, dtype=np.int64) for x in ds["input_ids"]]
    batches = []
    for i in range(0, len(seqs), FBS):
        chunk = seqs[i : i + FBS]
        lmax = max(len(s) for s in chunk)
        ids = np.zeros((len(chunk), lmax), dtype=np.int64)
        att = np.zeros_like(ids)
        for j, s in enumerate(chunk):
            ids[j, : len(s)] = s
            att[j, : len(s)] = 1
        batches.append((torch.from_numpy(ids), torch.from_numpy(att)))
    return np.array([len(s) for s in seqs]), batches


def main() -> None:
    assert torch.cuda.is_available(), "CUDA unavailable"
    lens, batches = real_batches()
    total_cells = sum(b[0].shape[0] for b in batches)
    worst = max(b[0].shape[1] for b in batches)
    print(f"model={MODEL_DIR.name} tissue={TISSUE} cells={total_cells} fbs={FBS} "
          f"padded_worst={worst} median_len={int(np.median(lens))}")
    print(f"{'variant':8s} {'s/run':>8s} {'cells/s':>9s} {'peak_GiB':>9s} "
          f"{'act_GiB/cell':>13s}")
    for tag in WANT:
        if tag not in VARIANTS:
            continue
        load_kw, qc = VARIANTS[tag]
        kw = dict(load_kw)
        if qc is not None:
            kw["quantization_config"] = qc
            kw["device_map"] = {"": 0}
        model = BertForMaskedLM.from_pretrained(str(MODEL_DIR), output_hidden_states=True, **kw)
        if qc is None:
            model = model.cuda()
        model = model.eval()
        weights = sum(p.numel() * p.element_size() for p in model.parameters()) / 2**30
        # activation cost of the binding (most padded) batch, per cell
        worst_batch = max((b for b in batches), key=lambda b: b[0].numel())
        cost_cells = worst_batch[0].shape[0]
        try:
            torch.cuda.reset_peak_memory_stats()
            torch.cuda.synchronize()
            t0 = time.time()
            with torch.no_grad():
                # two passes over the first batch to amortize warmup
                model(batches[0][0].cuda(), attention_mask=batches[0][1].cuda())
                for ids, att in batches:
                    model(ids.cuda(), attention_mask=att.cuda())
            torch.cuda.synchronize()
            dt = time.time() - t0
            peak = torch.cuda.max_memory_allocated() / 2**30
            act = max(peak - weights, 0.0) / cost_cells
            print(f"{tag:8s} {dt:8.1f} {total_cells/dt:9.1f} {peak:9.2f} {act:13.4f} "
                  f"(weights={weights:.2f} GiB, binding batch={cost_cells} cells x {worst_batch[0].shape[1]} tok)")
        finally:
            del model
            torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
