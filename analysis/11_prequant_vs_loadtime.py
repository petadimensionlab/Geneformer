#!/usr/bin/env python
"""Load-time quantization (Geneformer's `quantize` option) vs a pre-quantized
checkpoint saved by bitsandbytes.

Question this answers: is there any speed or accuracy difference between
  (A) quantizing at load time -> BertForMaskedLM.from_pretrained(fp32_dir,
      quantization_config=BitsAndBytesConfig(...), device_map={"":0}),
      which is what geneformer's model_type="Pretrained-Quantized" does, and
  (B) pre-quantizing once, model.save_pretrained(<dir>), then loading that
      already-quantized directory.

The runtime kernels and the quantized values should be the same, so the
expectation is: identical numbers, same inference speed, but (B) loads faster
and is smaller on disk. This script measures that instead of assuming it, and
checks whether the embeddings come out bit-identical.

Geneformer-specific quirk: save_pretrained() fails on this checkpoint because
`cls.predictions.bias` and `cls.predictions.decoder.bias` share storage
("Some tensors share memory ..."). The script breaks that sharing before
saving (numerically identical, just no longer the same tensor).

Usage:
    .venv/bin/python analysis/11_prequant_vs_loadtime.py
Environment:
    GENEFORMER_DIR / GENEFORMER_MODEL, QP_QUANT_DIR (default quantized/),
    QP_TISSUE (default AD_blood), QP_NCELLS (128), QP_LEN (2048),
    QP_BATCH (8), QP_ITERS (3)
"""
from __future__ import annotations

import gc
import json
import os
import shutil
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
QUANT_ROOT = Path(os.environ.get("QP_QUANT_DIR", ROOT / "quantized"))
TISSUE = os.environ.get("QP_TISSUE", "AD_blood")
DATA = ROOT / "input" / TISSUE / "tokenized" / f"{TISSUE}.dataset"
NCELLS = int(os.environ.get("QP_NCELLS", 128))
LEN = int(os.environ.get("QP_LEN", 2048))
BATCH = int(os.environ.get("QP_BATCH", 8))
ITERS = int(os.environ.get("QP_ITERS", 3))
VOCAB = 20275
MASK_TOKEN = 103


def int8() -> BitsAndBytesConfig:
    return BitsAndBytesConfig(load_in_8bit=True)


def nf4() -> BitsAndBytesConfig:
    return BitsAndBytesConfig(
        load_in_4bit=True,
        bnb_4bit_quant_type="nf4",
        bnb_4bit_use_double_quant=True,
        bnb_4bit_compute_dtype=torch.bfloat16,
    )


def real_batches():
    ds = load_from_disk(str(DATA))
    ds = ds.select(range(min(NCELLS, len(ds))))
    seqs = [np.asarray(x, dtype=np.int64) for x in ds["input_ids"]]
    out = []
    for i in range(0, len(seqs), BATCH):
        chunk = seqs[i : i + BATCH]
        lmax = max(len(s) for s in chunk)
        ids = np.zeros((len(chunk), lmax), dtype=np.int64)
        att = np.zeros_like(ids)
        for j, s in enumerate(chunk):
            ids[j, : len(s)] = s
            att[j, : len(s)] = 1
        out.append((torch.from_numpy(ids), torch.from_numpy(att)))
    return out


def mask_plan(batches):
    plan = []
    for ids, att in batches:
        p = (torch.rand(ids.shape) < 0.15) & (ids >= 3) & att.bool()
        for r in range(ids.shape[0]):
            if p[r].sum() == 0:
                idx = ((ids[r] >= 3) & att[r].bool()).nonzero().flatten()
                if len(idx):
                    p[r, idx[torch.randint(len(idx), (1,))]] = True
        plan.append(p)
    return plan


def load(path, quant_config=None, output_hidden_states=True, device_map=None,
         dtype=None):
    """Load a model, quantizing at load time when requested.

    A checkpoint that was itself saved from a quantized model carries
    `quantization_config` inside its config.json, so transformers quantizes it
    automatically on load -- and `.cuda()`/`.to()` is then forbidden. Both
    cases therefore need `device_map` instead of a device move.

    `device_map` may also be passed for a plain (non-quantized) load, which is
    the fair way to compare load times: `from_pretrained(...)` then `.cuda()`
    copies the weights twice (host -> host -> device), and that dominates.
    """
    kw = {"output_hidden_states": output_hidden_states}
    cfg_path = Path(path) / "config.json"
    already_quantized = False
    if cfg_path.exists():
        try:
            already_quantized = "quantization_config" in json.loads(cfg_path.read_text())
        except json.JSONDecodeError:
            already_quantized = False
    if quant_config is not None:
        kw["quantization_config"] = quant_config
    quantized_load = quant_config is not None or already_quantized
    if quantized_load or device_map is not None:
        kw["device_map"] = device_map or {"": 0}
    if dtype is not None:
        kw["torch_dtype"] = dtype
    gc.collect()
    torch.cuda.empty_cache()
    t0 = time.time()
    model = BertForMaskedLM.from_pretrained(str(path), **kw)
    if not quantized_load and device_map is None:
        model = model.cuda()
        if dtype is not None:
            model = model.to(dtype)
    torch.cuda.synchronize()
    return model.eval(), time.time() - t0


def evaluate(model, batches, plan, device="cuda"):
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
    return loss_sum / n, correct / n, torch.cat(embs)


def throughput(model):
    ids = torch.randint(3, VOCAB, (BATCH, LEN)).cuda()
    att = torch.ones_like(ids)
    torch.cuda.reset_peak_memory_stats()
    with torch.no_grad():
        model(ids, attention_mask=att)
        torch.cuda.synchronize()
        t0 = time.time()
        for _ in range(ITERS):
            model(ids, attention_mask=att)
        torch.cuda.synchronize()
        dt = (time.time() - t0) / ITERS
    return BATCH * LEN / dt, torch.cuda.max_memory_allocated() / 2**30


def dir_size_mib(path: Path) -> float:
    if not path.is_dir():
        return 0.0
    return sum(f.stat().st_size for f in path.rglob("*") if f.is_file()) / 2**20


def weights_mib(model) -> float:
    return sum(p.numel() * p.element_size() for p in model.parameters()) / 2**20


def build_prequantized(qc, out: Path, dtype=None) -> str:
    """Quantize (or cast) once from the fp32 checkpoint and save a loadable artifact."""
    if (out / "config.json").exists() and list(out.glob("*.safetensors")):
        return "already built"
    shutil.rmtree(out, ignore_errors=True)
    # Note: load WITHOUT output_hidden_states here. Loading with it sets
    # output_hidden_states=True on the generation config without also setting
    # return_dict_in_generate=True, and GenerationConfig.validate() then makes
    # save_pretrained() fail with a confusing ValueError. Weights are unaffected.
    model, _ = load(MODEL_DIR, qc, output_hidden_states=False, dtype=dtype)
    # break the shared storage between cls.predictions.bias and
    # cls.predictions.decoder.bias, otherwise safetensors refuses to save
    with torch.no_grad():
        bias = model.cls.predictions.bias.data.clone()
        model.cls.predictions.bias = torch.nn.Parameter(bias)
        model.cls.predictions.decoder.bias = torch.nn.Parameter(bias.clone())
    t0 = time.time()
    model.save_pretrained(str(out))
    dt = time.time() - t0
    del model
    gc.collect()
    torch.cuda.empty_cache()
    return f"saved in {dt:.1f}s -> {dir_size_mib(out):.1f} MiB"


def main() -> None:
    assert torch.cuda.is_available(), "CUDA unavailable"
    torch.manual_seed(0)
    np.random.seed(0)
    QUANT_ROOT.mkdir(parents=True, exist_ok=True)
    batches = real_batches()
    plan = mask_plan(batches)
    name = MODEL_DIR.name
    print(f"model={name} tissue={TISSUE} cells={NCELLS} bench_len={LEN} batch={BATCH}")
    print(f"fp32 checkpoint on disk: {dir_size_mib(MODEL_DIR):.1f} MiB\n")

    # --- warm up the CUDA context and bitsandbytes kernels so that the timed
    #     loads below are not charged for one-time initialization
    _ = torch.zeros(1024, device="cuda")
    for qc in (None, int8(), nf4()):
        m, _ = load(MODEL_DIR, qc)
        del m
        gc.collect()
        torch.cuda.empty_cache()

    pre = {}
    for tag, qc, dtype in (("int8", int8(), None), ("nf4", nf4(), None),
                           ("bf16", None, torch.bfloat16)):
        out = QUANT_ROOT / f"{name}-{tag}"
        note = build_prequantized(qc, out, dtype)
        cfg = json.loads((out / "config.json").read_text())
        print(f"pre-quantized {tag}: {note}; quantization_config in config.json="
              f"{'quantization_config' in cfg}; files={sorted(p.name for p in out.iterdir())}")
        pre[tag] = out
    print()

    results = {}

    def record(tag, path, quant_config, note="", device_map=None, dtype=None):
        model, load_cold = load(path, quant_config, device_map=device_map, dtype=dtype)
        del model
        gc.collect()
        torch.cuda.empty_cache()
        model, load_warm = load(path, quant_config, device_map=device_map, dtype=dtype)
        tok_s, peak = throughput(model)
        loss, acc, emb = evaluate(model, batches, plan)
        w = weights_mib(model)
        results[tag] = {"load_first_s": load_cold, "load_second_s": load_warm,
                        "tok_s": tok_s, "peak_gib": peak, "weights_mib": w,
                        "loss": loss, "mask_acc": acc, "emb": emb}
        print(f"[{tag:14s}] load={load_cold:4.2f}s/{load_warm:4.2f}s "
              f"tok/s={tok_s:8.0f} peak={peak:5.2f}GiB w={w:7.1f}MiB "
              f"loss={loss:.4f} acc={acc:.4f} {note}")
        del model
        gc.collect()
        torch.cuda.empty_cache()

    record("fp32", MODEL_DIR, None, "baseline (host->device double copy)")
    record("fp32_devmap", MODEL_DIR, None, "same weights, device_map load",
           device_map={"": 0})
    record("bf16_loadtime", MODEL_DIR, None, "fp32 checkpoint cast at load",
           device_map={"": 0}, dtype=torch.bfloat16)
    record("int8_loadtime", MODEL_DIR, int8(), "<- geneformer option")
    record("nf4_loadtime", MODEL_DIR, nf4())
    record("bf16_pre", pre["bf16"], None, f"<- {dir_size_mib(pre['bf16']):.1f}MiB on disk")
    record("int8_pre", pre["int8"], None, f"<- {dir_size_mib(pre['int8']):.1f}MiB on disk")
    record("nf4_pre", pre["nf4"], None, f"<- {dir_size_mib(pre['nf4']):.1f}MiB on disk")

    # what happens if the geneformer option is pointed at an already-quantized dir?
    for tag, qc in (("int8", int8()), ("nf4", nf4())):
        try:
            m, _ = load(pre[tag], qc)
            print(f"[check] {tag}_pre + quantization_config again: loaded OK, "
                  f"weights={weights_mib(m):.1f}MiB (fp32 would be "
                  f"{dir_size_mib(MODEL_DIR):.0f}MiB) <- possible double quantization")
            del m
            gc.collect()
            torch.cuda.empty_cache()
        except Exception as exc:
            print(f"[check] {tag}_pre + quantization_config again: "
                  f"FAILED {type(exc).__name__}: {str(exc)[:80]}")

    print("\n--- comparison ---")
    ref = results["fp32"]["emb"]
    for tag, r in results.items():
        cos = torch.nn.functional.cosine_similarity(ref, r["emb"], dim=1)
        print(f"{tag:14s} vs fp32: cos_mean={cos.mean():.6f} cos_min={cos.min():.6f} "
              f"dLoss={r['loss']-results['fp32']['loss']:+.4f}")
    for tag in ("int8", "nf4", "bf16"):
        a, b = results[f"{tag}_loadtime"], results[f"{tag}_pre"]
        print(f"\n{tag}: pre-saved vs converted at load time")
        print(f"  embeddings identical      : {torch.equal(a['emb'], b['emb'])}")
        print(f"  max |embedding difference|: {(a['emb']-b['emb']).abs().max().item():.3e}")
        print(f"  dLoss                     : {b['loss']-a['loss']:+.6f}")
        print(f"  tok/s                     : {a['tok_s']:.0f} -> {b['tok_s']:.0f} "
              f"({b['tok_s']/a['tok_s']:.3f}x)")
        print(f"  load time (2nd)           : {a['load_second_s']:.2f}s -> "
              f"{b['load_second_s']:.2f}s ({b['load_second_s']/a['load_second_s']:.3f}x)")
        print(f"  weights in device memory  : {a['weights_mib']:.1f} -> "
              f"{b['weights_mib']:.1f} MiB")

    out_json = QUANT_ROOT / f"prequant_vs_loadtime_{name}.json"
    out_json.write_text(json.dumps(
        {k: {kk: vv for kk, vv in v.items() if kk != "emb"} for k, v in results.items()},
        indent=2))
    print(f"\nresults -> {out_json}")


if __name__ == "__main__":
    main()
