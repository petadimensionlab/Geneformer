#!/usr/bin/env python
"""Quantization accuracy harness for Geneformer: MLM loss + embedding agreement.

Compares fp32 (reference) against bf16 / bnb-int8 / bnb-nf4 on REAL tokenized
cells, which is the only meaningful calibration-free accuracy probe for
Geneformer: the vocabulary is fixed gene + rank tokens, so synthetic input
would not represent the token distribution the quantizer sees.

Reports per variant:
  - masked-LM loss and top-1 accuracy on 15% masked positions
  - cell embedding (mean-pooled last hidden state) cosine similarity vs fp32
  - top-100 token overlap of the first cell's embedding ranking vs fp32

These are the cheap gates; downstream probes (04_baseline.py), fine-tuned
classifier metrics, and ISP rank agreement must also pass before adoption.
See docs/quantization/PLAN.md (Phase 2).

Usage:
    .venv/bin/python analysis/10b_quant_accuracy.py
Environment:
    GENEFORMER_DIR / GENEFORMER_MODEL, ADPD_TISSUE (default AD_blood),
    QA_NCELLS (default 256), QA_MAXLEN (default 2048), QA_BATCH (default 8)
"""
from __future__ import annotations

import os
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
NCELLS = int(os.environ.get("QA_NCELLS", 256))
MAXLEN = int(os.environ.get("QA_MAXLEN", 2048))
B = int(os.environ.get("QA_BATCH", 8))
MASK_RATE = 0.15
MASK_TOKEN = 103  # BERT [MASK]; Geneformer specials are 0/1/2, gene tokens >= 3


def build_batches():
    ds = load_from_disk(str(DATA))
    ds = ds.select(range(min(NCELLS, len(ds))))
    seqs = [np.asarray(x[:MAXLEN], dtype=np.int64) for x in ds["input_ids"]]
    out = []
    for i in range(0, len(seqs), B):
        chunk = seqs[i : i + B]
        lmax = max(len(s) for s in chunk)
        ids = np.zeros((len(chunk), lmax), dtype=np.int64)
        att = np.zeros_like(ids)
        for j, s in enumerate(chunk):
            ids[j, : len(s)] = s
            att[j, : len(s)] = 1
        out.append((torch.from_numpy(ids), torch.from_numpy(att)))
    return ds, out


def build_masks(batches):
    """Fixed mask plan so every variant sees identical masked positions."""
    plan = []
    for ids, att in batches:
        p = (torch.rand(ids.shape) < MASK_RATE) & (ids >= 3) & att.bool()
        for r in range(ids.shape[0]):
            if p[r].sum() == 0:
                idx = ((ids[r] >= 3) & att[r].bool()).nonzero().flatten()
                if len(idx):
                    p[r, idx[torch.randint(len(idx), (1,))]] = True
        plan.append(p)
    return plan


def run(tag, batches, plan, device, dtype=None, quant_config=None):
    kw = {}
    if dtype is not None:
        kw["torch_dtype"] = dtype
    if quant_config is not None:
        kw["quantization_config"] = quant_config
        kw["device_map"] = {"": 0}
    model = BertForMaskedLM.from_pretrained(str(MODEL_DIR), output_hidden_states=True, **kw)
    if quant_config is None:
        model = model.to(device)
        if dtype is not None:
            model = model.to(dtype)
    model.eval()

    loss_sum = n = correct = 0
    embs = []
    with torch.no_grad():
        for (ids, att), p in zip(batches, plan):
            ids, att, p = ids.to(device), att.to(device), p.to(device)
            masked = ids.clone()
            masked[p] = MASK_TOKEN
            out = model(input_ids=masked, attention_mask=att, output_hidden_states=True)
            logits, target = out.logits[p], ids[p]
            loss_sum += torch.nn.functional.cross_entropy(logits, target, reduction="sum").item()
            n += target.numel()
            correct += (logits.argmax(-1) == target).sum().item()
            h = out.hidden_states[-1] * att.unsqueeze(-1)
            embs.append((h.sum(1) / att.sum(1, keepdim=True)).float().cpu())
    emb = torch.cat(embs)
    print(f"{tag:6s} loss={loss_sum / n:.4f} maskacc={correct / n:.4f} n_masked={n}")
    del model
    if quant_config is None:
        torch.cuda.empty_cache()
    return emb


def main() -> None:
    assert torch.cuda.is_available(), "CUDA unavailable"
    assert DATA.is_dir(), f"tokenized dataset not found: {DATA} (expected a HF Dataset dir)"
    torch.manual_seed(0)
    np.random.seed(0)
    print(f"model={MODEL_DIR.name} tissue={TISSUE} ncells={NCELLS} maxlen={MAXLEN}")

    ds, batches = build_batches()
    lengths = np.array([len(x) for x in ds["input_ids"]])
    print(f"seq len: min={lengths.min()} median={int(np.median(lengths))} max={lengths.max()}")
    plan = build_masks(batches)

    device = "cuda"
    ref = run("fp32", batches, plan, device)
    results = {"fp32": ref}
    results["bf16"] = run("bf16", batches, plan, device, dtype=torch.bfloat16)
    results["int8"] = run(
        "int8", batches, plan, device, quant_config=BitsAndBytesConfig(load_in_8bit=True)
    )
    results["nf4"] = run(
        "nf4",
        batches,
        plan,
        device,
        quant_config=BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_use_double_quant=True,
            bnb_4bit_compute_dtype=torch.bfloat16,
        ),
    )

    print("\n--- cell embedding agreement vs fp32 (mean-pooled last hidden state) ---")
    print(f"{'variant':8s} {'cos_mean':>10s} {'cos_min':>10s} {'relL2':>10s} {'top100_overlap':>15s}")
    ref_rank = ref[0].argsort(descending=True)[:100].tolist()
    for tag, emb in results.items():
        cos = torch.nn.functional.cosine_similarity(ref, emb, dim=1)
        rel = (torch.norm(emb - ref, dim=1) / torch.norm(ref, dim=1)).mean()
        overlap = len(set(ref_rank) & set(emb[0].argsort(descending=True)[:100].tolist()))
        print(f"{tag:8s} {cos.mean():10.6f} {cos.min():10.6f} {rel:10.5f} {overlap:>10d}/100")
    print("\nGates (PLAN.md Phase 2): cos_mean >= 0.9995, cos_min >= 0.999, dLoss <= +0.02")


if __name__ == "__main__":
    main()
