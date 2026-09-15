#!/usr/bin/env python
"""Training-step GPU memory probe for Geneformer cell classifiers.

06_finetune.py fine-tunes BertForSequenceClassification; that is the memory
hungriest stage of the pipeline (forward + backward + optimizer states), and it
is what decides whether a model fits a given GPU budget. This measures peak
memory for one real training step per (model, dtype, batch, checkpointing)
configuration, optionally with 4-bit QLoRA (LoRA adapters on a frozen nf4 base,
as geneformer's CellClassifier(quantize={...}) path does).

Usage:
    .venv/bin/python analysis/10e_train_mem_probe.py
Environment:
    GENEFORMER_DIR / GENEFORMER_MODEL (default Geneformer-V2-104M)
    QT_BATCHES (default "8,16"), QT_LEN (default 4096), QT_NUM_LABELS (25)
    QT_VARIANTS (default "fp32,bf16,bf16+qc,nf4lora"), QT_LIMIT_GIB (default 100)
      qc = gradient checkpointing
"""
from __future__ import annotations

import os
from pathlib import Path

import torch
from transformers import BertForSequenceClassification, BitsAndBytesConfig

ROOT = Path(__file__).resolve().parent.parent
MODEL_DIR = Path(os.environ.get("GENEFORMER_DIR", ROOT / "geneformer_hf")) / os.environ.get(
    "GENEFORMER_MODEL", "Geneformer-V2-104M"
)
BATCHES = [int(x) for x in os.environ.get("QT_BATCHES", "8,16").split(",")]
LEN = int(os.environ.get("QT_LEN", 4096))
NLAB = int(os.environ.get("QT_NUM_LABELS", 25))
WANT = os.environ.get("QT_VARIANTS", "fp32,bf16,bf16+qc,nf4lora,nf4lora+qc").split(",")
LIMIT = float(os.environ.get("QT_LIMIT_GIB", 100)) * 2**30
VOCAB = 20275


def _patch_quantized_init_guard():
    """Work around transformers 4.46 + bitsandbytes: loading a *headless* quantized
    checkpoint (e.g. BertForSequenceClassification with a fresh classifier head)
    calls PreTrainedModel._initialize_weights via model.apply(), which tries
    weight.data.normal_() on already-quantized uint8 parameters and dies with
    `normal_kernel_cpu not implemented for 'Byte'`. Skip quantized modules.
    """
    from transformers.modeling_utils import PreTrainedModel

    original = PreTrainedModel._initialize_weights

    def guarded(self, module):
        weight = getattr(module, "weight", None)
        if weight is not None and weight.dtype == torch.uint8:
            return None
        return original(self, module)

    PreTrainedModel._initialize_weights = guarded


def build(tag):
    """Return (model, is_bnb_quantized) for a configuration tag."""
    grad_ckpt = tag.endswith("+qc")
    base = tag.replace("+qc", "")
    is_bnb = base == "nf4lora"
    kw = {}
    if base == "fp32":
        pass
    elif base == "bf16":
        kw["torch_dtype"] = torch.bfloat16
    elif base == "nf4lora":
        _patch_quantized_init_guard()
        kw["quantization_config"] = BitsAndBytesConfig(
            load_in_4bit=True,
            bnb_4bit_quant_type="nf4",
            bnb_4bit_use_double_quant=True,
            bnb_4bit_compute_dtype=torch.bfloat16,
        )
        kw["device_map"] = {"": 0}
    else:
        raise ValueError(tag)

    model = BertForSequenceClassification.from_pretrained(
        str(MODEL_DIR), num_labels=NLAB, **kw
    )
    if not is_bnb:
        model = model.cuda()
    model.train()
    if grad_ckpt:
        model.gradient_checkpointing_enable(gradient_checkpointing_kwargs={"use_reentrant": False})
    if base == "nf4lora":
        # geneformer's QuantInt path: r=64, alpha=128, dropout=0.1. It requests
        # task_type="TokenClassification", which peft >= 0.12 rejects (the fallback
        # it uses, TOKEN_CLS, is for token tasks); the cell classifier is a
        # sequence classifier, so SEQ_CLS is the correct task type here.
        from peft import LoraConfig, get_peft_model

        model.enable_input_require_grads()
        last_err = None
        for task_type in ("SEQ_CLS", "TOKEN_CLS"):
            try:
                cfg = LoraConfig(
                    lora_alpha=128, lora_dropout=0.1, r=64, bias="none",
                    task_type=task_type,
                )
                model = get_peft_model(model, cfg)
                break
            except ValueError as exc:
                last_err = exc
        else:
            raise last_err
        model = model.to("cuda")
    return model, is_bnb


def main() -> None:
    assert torch.cuda.is_available(), "CUDA unavailable"
    torch.manual_seed(0)
    print(f"model={MODEL_DIR.name} len={LEN} labels={NLAB}")
    print(f"{'config':14s} {'batch':>5s} {'step_s':>8s} {'peak_GiB':>9s}  notes")
    for tag in WANT:
        for batch in BATCHES:
            est = batch * LEN * 1152 * 4 * 40  # rough, scaled to the biggest model
            if est > LIMIT * 1.5:
                print(f"{tag:14s} {batch:5d}  skipped: estimated peak > {LIMIT/2**30:.0f} GiB "
                      f"(set QT_LIMIT_GIB to force)")
                continue
            try:
                model, quant = build(tag)
                trainable = [p for p in model.parameters() if p.requires_grad]
                opt = torch.optim.AdamW(trainable, lr=1e-4)
                ids = torch.randint(3, VOCAB, (batch, LEN)).cuda()
                att = torch.ones_like(ids)
                labels = torch.randint(0, NLAB, (batch,)).cuda()
                torch.cuda.empty_cache()
                torch.cuda.reset_peak_memory_stats()
                t0 = torch.cuda.Event(enable_timing=True)
                t1 = torch.cuda.Event(enable_timing=True)
                t0.record()
                out = model(input_ids=ids, attention_mask=att, labels=labels)
                out.loss.backward()
                opt.step()
                opt.zero_grad(set_to_none=True)
                t1.record()
                torch.cuda.synchronize()
                peak = torch.cuda.max_memory_allocated() / 2**30
                ntrain = sum(p.numel() for p in trainable) / 1e6
                print(f"{tag:14s} {batch:5d} {t0.elapsed_time(t1)/1000:8.2f} {peak:9.2f}  "
                      f"trainable={ntrain:.1f}M")
            except torch.cuda.OutOfMemoryError as exc:
                print(f"{tag:14s} {batch:5d}      OOM  {str(exc).splitlines()[0][:70]}")
            except Exception as exc:
                print(f"{tag:14s} {batch:5d}  FAILED  {type(exc).__name__}: {str(exc)[:70]}")
            finally:
                model = opt = out = ids = att = labels = None
                torch.cuda.empty_cache()


if __name__ == "__main__":
    main()
