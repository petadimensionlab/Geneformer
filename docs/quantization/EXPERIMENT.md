# Geneformer quantization — consolidated experiment log

**English** · [日本語](EXPERIMENT-jp.md)

All numbers below were measured on the machine used for this work. Nothing here
is estimated unless explicitly labelled as such.

| | |
|---|---|
| Host | NVIDIA GB10 (aarch64, sm_121), **119.63 GiB** unified memory shared with the OS, 20 cores |
| Software | torch 2.13.0+cu130, transformers 4.46.3, bitsandbytes 0.50.1, peft 0.20.0, datasets 4.0.0, Python 3.12.13 |
| Models | `Geneformer-V2-104M` (12 layers, hidden 768, 104.4M params), `Geneformer-V2-316M` (18 layers, hidden 1152, 316.4M params) |
| Vocabulary / max length | 20,275 tokens / 4,096 (identical for both models) |
| Data | tokenized single-cell datasets under `input/<TISSUE>/tokenized/` (6 tissues) |

---

## 1. Conclusions

1. **bf16 is the best trade-off and it is not really a quantization choice:
   ~5x faster than fp32 with no measurable accuracy loss** (embedding cosine
   0.99998, MLM loss difference +0.001).
2. **int8 via bitsandbytes is safe but slow** — 1.6-1.7x faster than fp32,
   which is less than half of bf16's gain. Embedding cosine 0.9997.
3. **4-bit nf4 degrades accuracy in a way that matters**: MLM loss worsens and
   the *ranking* of genes changes (top-100 overlap 87/100). Since whole
   pipelines here produce gene rankings (in silico perturbation), nf4 is not
   usable as-is.
4. **Scaling up cells or genes does not unlock an int8/int4 advantage**
   (measured batch 8 → 256, sequence length 512 → 4,096, both models).
   The speed ordering and the ratios are constant, and at large batches the
   memory advantage of int4 over bf16 is 0.12 GiB out of 30 GiB.
5. **For fine-tuning, quantization is counterproductive**: bf16 + gradient
   checkpointing is both smaller and faster than QLoRA (4-bit + LoRA).
6. **Pre-quantizing weights and saving them does not change inference speed or
   accuracy at all** (bit-identical embeddings), but it makes **loading 2-6x
   slower**. Its only benefit is file size.
7. **Practical recommendation**: use **bf16 everywhere** for inference and
   fine-tuning; use int8 only when bf16 is unavailable or to distribute/shrink
   weights; do not use nf4 for any ranking-producing output.

---

## 2. Inference throughput (synthetic input, batch 8, tokens/s)

| Seq len | 104M fp32 | 104M bf16 | 104M int8 | 104M nf4 | 316M fp32 | 316M bf16 | 316M int8 | 316M nf4 |
|---|---|---|---|---|---|---|---|---|
| 512 | 43,590 | **205,466** | 69,489 | 127,436 | 15,656 | **78,059** | 31,766 | 53,281 |
| 2048 | 35,753 | **174,696** | 60,290 | 117,262 | 13,404 | **64,224** | 26,619 | 45,590 |
| 4096 | 29,672 | **154,269** | 55,814 | 107,731 | 11,407 | **55,407** | 24,800 | 40,729 |

Relative to fp32: **bf16 ≈ 4.9-5.2x**, nf4 ≈ 3.3-3.6x, int8 ≈ 1.6-1.7x.
Ordering is **bf16 > nf4 > int8** at every length and both model sizes.
int8 is slower than nf4 because bitsandbytes computes in fp16 and dequantizes
weights on the fly; this GPU's bf16 throughput is high enough that int8's
smaller weights never pay off.

`analysis/10_quant_bench.py`

## 3. Peak memory (GiB)

**104M, seq len 4096**, peak including weights:

| batch | fp32 | bf16 | int8 | nf4 |
|---|---|---|---|---|
| 8 | 4.21 | 2.13 | 2.04 | 2.00 |
| 32 | 15.57 | 7.81 | 7.72 | 7.69 |
| 64 | 30.72 | 15.38 | 15.30 | 15.26 |

**316M, seq len 4096**: batch 8 → fp32 6.64 / bf16 3.40 / int8 3.00 / nf4 2.87;
batch 32 → fp32 22.36 / bf16 11.20 / int8 10.93 / nf4 10.80.

Per cell at seq len 4096: **fp32 ≈ 0.47 GiB, bf16/int8/nf4 ≈ 0.24 GiB.**
The 2x reduction comes from activations being half precision, not from the
weights — which is why bf16 captures the entire memory benefit.

`analysis/10_quant_bench.py`, `analysis/10c_quant_memscale.py`

## 4. Real-data throughput and memory

AD_blood, 512 cells, `forward_batch_size=64`, longest padded batch 3,226 tokens:

| Model | dtype | seconds | cells/s | peak GiB |
|---|---|---|---|---|
| 104M | fp32 | 49.5 | 10.3 | 24.29 |
| 104M | bf16 | 14.1 | **36.3** | 12.16 |
| 104M | int8 | 32.6 | 15.7 | 12.08 |
| 104M | nf4 | 18.2 | 28.1 | 12.04 |
| 316M | fp32 | 130.4 | 3.9 | 34.53 |
| 316M | bf16 | 35.8 | **14.3** | 17.28 |
| 316M | int8 | 71.1 | 7.2 | 17.02 |
| 316M | nf4 | 45.4 | 11.3 | 16.89 |

**316M in bf16 is 1.4x faster than 104M in fp32 while using less memory** — the
central argument for moving to the larger model.

AD_brain (longest tissue: median 2,770 / max 4,096 tokens), 200 cells,
`forward_batch_size=100`, longest padded batch 4,096:

| Model | dtype | cells/s | peak GiB |
|---|---|---|---|
| 104M | fp32 | 4.2 | **47.77 (at the 48 GB wall)** |
| 104M | bf16 | 15.1 | 23.90 |
| 316M | fp32 | 1.6 | 67.31 |
| 316M | bf16 | 5.8 | 33.67 |

`analysis/10d_real_workload_mem.py`

## 5. Does scale help int8/int4? No.

104M, seq len 2048, throughput (tokens/s) and peak memory:

| batch | fp32 | bf16 | int8 | nf4 |
|---|---|---|---|---|
| 8 | 35,753 (2.32 GiB) | 173,872 (1.18) | 60,190 (1.10) | 116,906 (1.06) |
| 128 | 35,458 (30.72) | 178,157 (15.38) | 58,236 (15.30) | 118,756 (15.26) |
| 256 | 36,055 (61.02) | 175,173 (30.53) | 58,042 (30.45) | 118,829 (**30.41**) |

- Throughput is flat from batch 8 to 256: the GPU is already saturated at
  batch 8, and batching does not buy throughput. Ratios never change.
- At batch 256, **bf16 and nf4 differ by 0.12 GiB** (30.41 vs 30.53 GiB).
  int4's weight saving (0.39 → 0.08 GiB for 104M) is negligible next to
  activations once the batch is large.
- More genes in an ISP run only multiply the number of forward passes, so the
  per-forward ratios carry over unchanged. int8/int4 remain useful only for
  **artifact size** and **CPU inference**.

`analysis/10c_quant_memscale.py`

## 6. Accuracy (real cells, 15% masked positions)

AD_blood, 256 cells, 58,888 masked positions. Embedding = mean-pooled last
hidden state; "top-100" = overlap of the token ranking of the first cell's
embedding against fp32.

**104M**

| dtype | MLM loss | mask acc | cos mean | cos min | top-100 |
|---|---|---|---|---|---|
| fp32 | 4.0919 | 0.1135 | 1.000000 | 1.000000 | 100/100 |
| bf16 | 4.0935 | 0.1135 | 0.999976 | 0.999971 | 100/100 |
| int8 | 4.0823 | 0.1135 | 0.999705 | 0.999635 | 96/100 |
| nf4 | 4.3063 | 0.1053 | 0.981259 | 0.970528 | **87/100** |

**316M** (also shows the pretraining-quality gain: loss 3.64 vs 4.09,
mask accuracy 0.137 vs 0.114)

| dtype | MLM loss | mask acc | cos mean | cos min | top-100 |
|---|---|---|---|---|---|
| fp32 | 3.6410 | 0.1368 | 1.000000 | 1.000000 | 100/100 |
| bf16 | 3.6327 | 0.1368 | 0.999949 | 0.999920 | 98/100 |
| int8 | 3.6348 | 0.1367 | 0.999650 | 0.999541 | 98/100 |
| nf4 | 3.5984 | 0.1386 | 0.978739 | 0.969692 | **87/100** |

Note that nf4's *average* loss is not worse on 316M (3.5984 vs 3.6410) yet its
embeddings are off by 2% and one gene in eight moves in the top-100. **Judge
quantization by ranking agreement, not by a scalar loss.**

`analysis/10b_quant_accuracy.py`

## 7. Fine-tuning: quantization is counterproductive

One training step of a `BertForSequenceClassification` cell classifier,
batch 8, seq len 4096, 25 labels (gradient checkpointing = "qc"):

| Config | 104M peak GiB | 104M step s | 316M peak GiB | 316M step s |
|---|---|---|---|---|
| fp32 | 15.40 | 4.05 | 34.22 | 10.05 |
| bf16 | 8.03 | 1.25 | 17.80 | 2.53 |
| **bf16 + qc** | **1.82** | **0.90** | **3.41** | **2.43** |
| nf4 + LoRA (QLoRA) | 9.83 | 1.60 | 21.56 | 3.29 |
| nf4 + LoRA + qc | 1.95 | 1.49 | 3.39 | 3.63 |

- **QLoRA uses more memory and is slower than a bf16 full fine-tune**
  (9.83 GiB / 1.60 s vs 8.03 GiB / 1.25 s on 104M), because bitsandbytes
  stores packed 4-bit weights plus dequantization buffers and computes in fp16.
- The only lever that actually reduces training memory is **gradient
  checkpointing** (8.03 → 1.82 GiB on 104M).
- `analysis/06_finetune.py` already enables gradient checkpointing; add
  `bf16: True` to its training args for the ~4x speedup.

`analysis/10e_train_mem_probe.py`

## 8. Where quantization helps, stage by stage

| Stage | Uses the model? | Cost driver | Quantization effect | Use |
|---|---|---|---|---|
| tokenization | **no** (numpy/pandas + gene dictionaries only) | CPU, disk I/O, `nproc` | none at all | — |
| embedding extraction | forward only | number of cells | large (speed) | **bf16** |
| in silico perturbation | forward only | cells × genes × timepoints | large (speed), but ranking accuracy is the deliverable | **bf16**; int8 only with a rank-agreement check; nf4 no |
| fine-tuning | forward + backward | training steps | **negative** | **bf16 + gradient checkpointing** |

Details and the full accuracy gate design (7 criteria, including ISP
Spearman ρ ≥ 0.99 and top-20 overlap ≥ 18/20) are in
[STAGES.md](STAGES.md) and [PLAN.md](PLAN.md) (Japanese).

## 9. Pre-quantized weights vs quantizing at load time

Two paths compared:
(A) `from_pretrained(fp32_dir, quantization_config=..., device_map={"":0})` —
what geneformer's `model_type="Pretrained-Quantized"` does;
(B) quantize once, `save_pretrained()`, then load the saved directory.

**Inference: no difference at all.** Embeddings are **bit-identical**
(max |difference| 0.000e+00), MLM loss identical to the last digit, throughput
ratio 0.99-1.01x. 104M int8: 59,277 vs 58,943 tok/s.

**Load time: path (B) is slower.**

| Variant | 104M median load | 316M median load | Disk size (104M / 316M) |
|---|---|---|---|
| fp32 | **3.031 s** | **7.251 s** | 398.2 / 1,206.8 MiB |
| bf16 (cast at load) | **0.168 s** | **0.255 s** | 398.2 / 1,206.8 MiB |
| int8 (quantize at load) | **0.217 s** | **0.384 s** | 398.2 / 1,206.8 MiB |
| nf4 (quantize at load) | **0.230 s** | **0.398 s** | 398.2 / 1,206.8 MiB |
| int8 (pre-quantized) | 0.443 s (2.0x slower) | 2.476 s (**6.4x slower**) | 118.5 / 330.8 MiB |
| nf4 (pre-quantized) | 0.590 s (2.6x slower) | 1.442 s (3.6x slower) | 79.1 / 197.9 MiB |
| bf16 (pre-saved, no dtype arg) | 0.149 s | 0.275 s | 199.2 / 603.5 MiB → **loads as fp32** |

Also note **fp32 loading is 13-28x slower than any converted loading path**,
regardless of `device_map` / `low_cpu_mem_usage` / explicit `torch_dtype=float32`.
Practical rule: **always pass a dtype (or a quantization config) when loading.**

Summary: **pre-quantize for distribution size, never for speed.**
`analysis/11_prequant_vs_loadtime.py`, `analysis/11b_load_time_bench.py`

## 10. Environment facts worth knowing

- Tokenized sequence lengths differ a lot by tissue and drive the memory budget:

  | tissue | cells | median | 95th pct | max |
  |---|---|---|---|---|
  | AD_blood | 48,909 | 1,215 | 2,415 | 3,670 |
  | AD_brain | 41,858 | **2,999** | **4,096** | **4,096** |
  | AD_liver | 37,453 | 877 | 2,592 | 4,096 |
  | AD_smallint | 61,951 | 1,085 | 2,319 | 3,877 |
  | AD_spleen | 76,947 | 840 | 1,602 | 2,848 |
  | PD_spleen | 83,783 | 925 | 1,922 | 3,587 |

- 48 GB GPU budget (worst-case 4,096-token batches, 2 GiB reserved):
  fp32 ≈ 96 cells per forward, bf16/int8/nf4 ≈ 194 cells.
- On GB10, GPU and host memory are the same 119.63 GiB pool. Using more VRAM
  leaves less for the ISP perturbation dataset, so `forward_batch_size` above
  the default buys nothing and costs dataset headroom.
- `Geneformer-V1-10M`, `V2-104M_CLcancer` and `V2-316M` shipped as **git-lfs
  pointer stubs** (133-135 bytes). Fetch the real weights first:
  `./download.sh --model V2-316M`.
- The ISP scripts resolve a **fine-tuned cell classifier** from
  `input/<TISSUE>/runs/*/ksplit1`. Running ISP with 316M therefore requires
  fine-tuning a 316M classifier per tissue first.

## 11. Pitfalls found (with fixes)

1. **Never call `.to()` / `.cuda()` / `.half()` on a quantized model.**
   `perturber_utils.load_model()` already skips its device move when quantizing;
   a checkpoint that carries `quantization_config` in `config.json` must be
   placed with `device_map={"": 0}`.
2. **4-bit is unreachable through `model_type`.** `"Pretrained-Quantized"` and
   `"MTLCellClassifier-Quantized"` map to int8 only. Exposing `quantize` on
   `EmbExtractor` / `InSilicoPerturber` (one optional argument, propagated to
   `load_model`) is the minimal patch.
3. **`save_pretrained()` fails: "Some tensors share memory".**
   `cls.predictions.bias` and `cls.predictions.decoder.bias` share storage;
   safetensors refuses to write it. `clone()` one of them before saving
   (numerically identical). Implemented in `analysis/11_prequant_vs_loadtime.py`.
4. **`save_pretrained()` fails after loading with `output_hidden_states=True`.**
   `GenerationConfig.validate()` rejects the combination
   (`output_hidden_states` without `return_dict_in_generate`). Load without
   `output_hidden_states` when the model is only going to be saved.
5. **A pre-saved bf16 checkpoint loads as fp32** unless `torch_dtype` is passed
   explicitly — 2x memory and 5x slower inference (67,338 → 13,428 tok/s on
   316M). `config.json` alone was not enough; pass the argument.
6. **A 4-bit `BertForSequenceClassification` cannot be loaded as-is**:
   `normal_kernel_cpu not implemented for 'Byte'`. The fresh classifier head
   triggers `_initialize_weights` on already-quantized uint8 weights. Guard it
   (see `_patch_quantized_init_guard()` in `analysis/10e_train_mem_probe.py`).
7. **`LoraConfig(task_type="TokenClassification")` raises on peft ≥ 0.12.**
   Valid names are `TOKEN_CLS` / `SEQ_CLS`; the cell classifier is a sequence
   classifier, so `SEQ_CLS` is correct. geneformer falls back to `TOKEN_CLS`,
   which works but is the wrong task type.
8. **`06_finetune.py` writes `runs/<date>_geneformer_cellClassifier_<PREFIX>_celltype/`
   with no model name.** Pointing it at 316M silently reuses the existing 104M
   checkpoint ("reusing existing checkpoint"). Move `runs/` aside or add the
   model name to `RUN_PREFIX` first.
9. **Do not calibrate a quantization on synthetic input.** The vocabulary is
   fixed gene + rank tokens, so the error depends on the real token
   distribution. Use real tokenized cells, stratified by cell type and length,
   and keep the tissues used for final evaluation out of calibration.
10. **Keep `datasets==4.0.0`** (`>=5` makes `perturb_data`'s `dataset.map` hang),
    and keep `transformers==4.46.3` (5.x breaks `SpecialTokensMixin`).
11. **bf16 dies in the embedding export** — numpy has no bfloat16 dtype, so
    `embs.cpu().numpy()` raises `TypeError: Got unsupported ScalarType
    BFloat16` (`emb_extractor.py:256`, 283, 702, `perturber_utils.py:824`).
    Insert `.float()` before `.cpu().numpy()` (no-op on fp32). Patch + note:
    `patches/bf16/`.
12. **A single tied cell aborts the evaluation of a fine-tune.**
    `evaluation_utils.vote()` returns the string `"tie"` on an exact float tie,
    which then reaches sklearn as a mixed-type label array:
    `ValueError: Mix of label input types (string and number)`. The checkpoint
    is already written, but the run exits non-zero and no metrics/report are
    produced. Break ties numerically (lowest class index). Patch + note:
    `patches/eval_tie/`.

## 12. End-to-end ISP measurement (V2-316M, AD_spleen) — see REPORT-316M-ISP.md

Once the 316M path was wired up, the whole pipeline was measured on AD_spleen
rather than on synthetic batches:

| | bf16 | fp32 | change |
|---|---|---|---|
| ISP canary, 55 genes × 3 timepoints × 200 cells | **46 min 25 s, 34.9 Wh** | ≈1 h 46 m, ≈133 Wh *(extrapolated)* | 2.28x faster, 3.82x less energy |
| ISP matched pair, 24 genes × 3 tp × 100 cells | **17 min 58 s, 12.5 Wh, 41.6 W** | 40 min 57 s, 47.7 Wh, 69.9 W | −56% time, −74% energy, −40% power |
| fine-tune, 1 epoch (3,057 steps) | **1 h 59 m 46 s, 137.5 Wh** | ≈8 h *(extrapolated)* | ~4.0x |

Rank agreement of the matched pair (same 316M checkpoint, only dtype differs):
**Spearman ρ = 0.9839, top-20 overlap 20/20**, mean |Δ| = 8.3e-05 against a mean
|Shift| of 7.5e-04. All 4 sign flips sit in the smallest shifts
(|Shift| < 2.1e-04); above the median |Shift| there are **zero** flips. So bf16
changes only the part of the ranking that carries no decision. The
model change (104M fp32 → 316M bf16) is a different story: ρ = 0.588, top-20
13/20, mean |Δ| = 3.2e-03 (82% of the signal).

Consequence for the gate in `PLAN.md` Phase 2: an absolute ρ ≥ 0.99 and a 98%
sign agreement are **not attainable for `Shift_to_goal_end`**, because the metric
is a difference of two ≈1 cosines and therefore carries only 2-3 significant
digits. Use top-N overlap, ρ ≥ 0.98, and sign agreement restricted to genes
above the measured noise floor (3 × median |Δ|). Details and the full tables:
[REPORT-316M-ISP.md](REPORT-316M-ISP.md).

## 13. Reproduce

```bash
cd ~/workspace/research/Geneformer
.venv/bin/python analysis/10_quant_bench.py           # §2, §3 throughput and memory
.venv/bin/python analysis/10c_quant_memscale.py       # §3, §5 batch scaling
.venv/bin/python analysis/10d_real_workload_mem.py    # §4 real-data memory
.venv/bin/python analysis/10b_quant_accuracy.py       # §6 accuracy
.venv/bin/python analysis/10e_train_mem_probe.py      # §7 training memory
.venv/bin/python analysis/11_prequant_vs_loadtime.py  # §9 pre-quantized vs load-time
.venv/bin/python analysis/11b_load_time_bench.py      # §9 load times

# end-to-end (§13) — profiled ISP on the 316M classifier, fp32 vs bf16, and analysis
.venv/bin/python analysis/13_profile.py --label my-run --out docs/quantization/profiles -- \
  .venv/bin/python analysis/07_ad_spleen_early_isp.py     # with GF_DTYPE=bf16 / ISP_EXPERIMENT=...
.venv/bin/python analysis/14_compare_isp.py --a <run A> --b <run B> --label x
.venv/bin/python analysis/15_isp_timing.py docs/quantization/profiles/<label>.log

# any model / tissue / batch size via environment variables
GENEFORMER_MODEL=Geneformer-V2-316M ADPD_TISSUE=AD_brain QW_FBS=100 \
  .venv/bin/python analysis/10d_real_workload_mem.py
```

Artifacts are written under `quantized/` (gitignored) and raw results as JSON
there as well. Profiled runs write to `docs/quantization/profiles/`.

## 14. Related documents (Japanese)

| Document | Content |
|---|---|
| [PLAN.md](PLAN.md) | Overall quantization plan: goals, phases, acceptance gates |
| [STAGES.md](STAGES.md) | Stage-by-stage impact and recommended settings |
| [PREQUANT-VS-LOADTIME.md](PREQUANT-VS-LOADTIME.md) | The pre-quantization investigation in detail |
| [PLAN-316M-128GB.md](PLAN-316M-128GB.md) | Plan for running V2-316M on a 128 GB machine |
| [REPORT-104M-48GB.md](REPORT-104M-48GB.md) | Report: does V2-104M fit a 48 GB GPU |
| [README.md](README.md) | Index of this documentation set |
