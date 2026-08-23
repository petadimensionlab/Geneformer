#!/usr/bin/env python
"""Empirically measure peak MPS memory of one Geneformer fine-tune train step.

Runs a single forward+backward through the same CellClassifier model used by
06_finetune.py at the worst-case padded sequence length, for several batch
sizes / gradient-checkpointing settings, and reports peak memory so we can pick
a per-device batch size that fits this machine.

Usage:
    python analysis/measure_mps_batch.py            # sweep, prints recommendation
    python analysis/measure_mps_batch.py 4 4096 0   # single probe: B, seq_len, grad_ckpt
"""
from __future__ import annotations

import json
import os
import subprocess
import sys

SEED = 42


def _probe(batch_size: int, seq_len: int, grad_ckpt: bool) -> dict:
    """Run in a fresh process so MPS allocator state doesn't leak between probes."""
    repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    code = r"""
import json, os, sys, torch
batch_size, seq_len, grad_ckpt = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3] == "1"
repo_root = sys.argv[4]

sys.path.insert(0, os.path.join(repo_root, "geneformer_hf"))
MODEL_DIR = os.path.join(repo_root, "geneformer_hf", "Geneformer-V2-104M")

from transformers.models.bert.modeling_bert import BertForSequenceClassification

torch.manual_seed(42)
model = BertForSequenceClassification.from_pretrained(
    MODEL_DIR, num_labels=20, problem_type="single_label_classification",
    attn_implementation="eager")
# mirror 06_finetune.py: freeze the first 6 layers
for param in model.bert.embeddings.parameters():
    param.requires_grad = False
for layer in model.bert.encoder.layer[:6]:
    for param in layer.parameters():
        param.requires_grad = False
if grad_ckpt:
    model.gradient_checkpointing_enable(
        gradient_checkpointing_kwargs={"use_reentrant": False})
model.to("mps").train()

x = torch.randint(low=0, high=20275, size=(batch_size, seq_len), device="mps")
attn = torch.ones_like(x)
labels = torch.randint(low=0, high=20, size=(batch_size,), device="mps")

base = torch.mps.driver_allocated_memory()
out = model(input_ids=x, attention_mask=attn, labels=labels)
loss = out.loss
loss.backward()
peak = torch.mps.driver_allocated_memory()
step_peak = torch.mps.current_allocated_memory()
print(json.dumps({
    "ok": True,
    "batch_size": batch_size, "seq_len": seq_len, "grad_ckpt": grad_ckpt,
    "base_gib": base / 2**30,
    "peak_gib": peak / 2**30,
    "live_after_backward_gib": step_peak / 2**30,
}))
"""
    venv = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), ".venv")
    py = os.path.join(venv, "bin", "python")
    env = dict(os.environ, WANDB_DISABLED="true")
    r = subprocess.run([py, "-c", code, str(batch_size), str(seq_len),
                        "1" if grad_ckpt else "0", repo_root],
                       capture_output=True, text=True, env=env, timeout=900)
    for line in r.stdout.splitlines():
        line = line.strip()
        if line.startswith("{"):
            return json.loads(line)
    return {"ok": False, "batch_size": batch_size, "seq_len": seq_len,
            "grad_ckpt": grad_ckpt, "stderr": r.stderr.strip().splitlines()[-1]
            if r.stderr.strip() else "no output"}


def main() -> None:
    if len(sys.argv) == 4:
        res = _probe(int(sys.argv[1]), int(sys.argv[2]), sys.argv[3] == "1")
        print(json.dumps(res, indent=2))
        return

    SEQ_LEN = 4096          # worst-case padded length (p100 of PD_BM lengths)
    OPTIMIZER_GIB = 1.7     # AdamW moments + fp32 grads for 104M params
    HEADROOM_GIB = 8        # allocator fragmentation + eval/predict spikes
    RAM_BUDGET_GIB = 100    # keep inside physical RAM (128 GB) to avoid swap

    results = []
    for ckpt in (False, True):
        for b in (8, 6, 4, 3, 2):
            print(f"probing B={b} L={SEQ_LEN} grad_ckpt={ckpt} ...", flush=True)
            res = _probe(b, SEQ_LEN, ckpt)
            print("  ->", json.dumps(res), flush=True)
            results.append(res)
            if res.get("ok"):
                need = res["peak_gib"] + OPTIMIZER_GIB + HEADROOM_GIB
                if need > RAM_BUDGET_GIB and b <= 4:
                    break  # larger B only gets worse; stop descending further? no—ascending
            if res.get("ok") and res["peak_gib"] + OPTIMIZER_GIB + HEADROOM_GIB > RAM_BUDGET_GIB * 1.6:
                break  # hopeless, try next setting
    print("\n=== summary ===")
    for r in results:
        if r.get("ok"):
            tag = "grad-ckpt " if r["grad_ckpt"] else "full-bwd   "
            print(f"{tag}B={r['batch_size']:>2} L={r['seq_len']}  "
                  f"peak {r['peak_gib']:6.1f} GiB  (+~{OPTIMIZER_GIB + HEADROOM_GIB:.0f} optimizer/headroom)")
        else:
            print(f"B={r['batch_size']} FAILED: {r.get('stderr', '')[-120:]}")


if __name__ == "__main__":
    main()
