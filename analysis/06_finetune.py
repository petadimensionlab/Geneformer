#!/usr/bin/env python
"""Fine-tune a Geneformer cell classifier on the ADPD PD Blood cohort.

Ports the tutorial's cells 20-23 to our data. Same splits as the frozen-embedding
baseline, so the two are directly comparable:
  train = replicate _1  (12 samples)   <- what the baseline probe was fit on
  eval  = replicate _2  (12 samples)   <- model selection during training only
  test  = replicate _3  (12 samples)   <- identical to the baseline's test set
"""
from __future__ import annotations

import json
import os
import pickle
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import torch
from datasets import load_from_disk
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, f1_score

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import resolve as _resolve_tissue

SEED = 42
np.random.seed(SEED)

# ---------------------------------------------------------------- configuration
# Resolve the tissue workspace from input/<TISSUE>/h5ad/*.h5ad. See
# _resolve_tissue.py (ADPD_TISSUE / ADPD_PREFIX / ADPD_ROOT).
ROOT, PREFIX = _resolve_tissue()
GF_ROOT = Path(os.environ.get(
    "GENEFORMER_DIR",
    Path.home() / "Documents/geneformer-uv-starter/geneformer-workspace/Geneformer",
))
MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
MODEL_DIR = GF_ROOT / MODEL_NAME
if not MODEL_DIR.is_dir():
    raise SystemExit(
        f"Geneformer checkpoint not found: {MODEL_DIR}\n"
        "Set GENEFORMER_DIR (and GENEFORMER_MODEL) or run ./setup.sh."
    )
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"
RUN_DIR = ROOT / "runs"
TABLE_DIR = ROOT / "results/tables"
FIG_DIR = ROOT / "results/figures"
for d in (RUN_DIR, TABLE_DIR, FIG_DIR):
    d.mkdir(parents=True, exist_ok=True)

RUN_PREFIX = f"{PREFIX}_celltype"
# The trained checkpoint path is runs/<date>_geneformer_cellClassifier_<RUN_PREFIX>/ksplit1
# and the model name is not part of it, so pointing this script at a different
# model (e.g. Geneformer-V2-316M) would silently reuse the existing classifier
# ("reusing existing checkpoint") instead of training a new one. Tag non-default
# models — or anything set via FINETUNE_RUN_SUFFIX — so each model gets its own
# directory. The default model keeps the original name for backwards compatibility.
RUN_SUFFIX = os.environ.get(
    "FINETUNE_RUN_SUFFIX",
    "" if MODEL_NAME == "Geneformer-V2-104M" else MODEL_NAME,
)
if RUN_SUFFIX:
    RUN_PREFIX = f"{PREFIX}_celltype_{RUN_SUFFIX}"
    print(f"[config] run prefix -> {RUN_PREFIX}", flush=True)
from geneformer.device import get_device

device = get_device()
print(f"fine-tuning on device: {device}", flush=True)

# ------------------------------------------------------------------- splits
tokens = load_from_disk(str(TOKENIZED))
samples = sorted(set(tokens["individual"]))
by_rep: dict[str, list] = {}
for s in samples:
    rep = s.rsplit("_", 1)[-1]
    if rep in ("1", "2", "3"):
        by_rep.setdefault(rep, []).append(s)
TRAIN, EVAL, TEST = by_rep.get("1", []), by_rep.get("2", []), by_rep.get("3", [])
# Replicate-cohort split: train=_1, eval=_2, test=_3. Tissue-dependent sample
# counts (PD_*: 12/rep; AD_*: 9-10/rep) — only require each split non-empty.
assert TRAIN and EVAL and TEST, (
    f"need non-empty _1/_2/_3 replicate cohorts; got "
    f"train={len(TRAIN)} eval={len(EVAL)} test={len(TEST)}"
)
assert not (set(TRAIN) & set(EVAL) | set(TRAIN) & set(TEST) | set(EVAL) & set(TEST))
print(f"samples  train {len(TRAIN)}  eval {len(EVAL)}  test {len(TEST)}", flush=True)

CLASSES = sorted(set(tokens["celltype"]))
print(f"{len(CLASSES)} classes", flush=True)

# bf16 autocast for training (FINETUNE_BF16=1). Measured ~4x faster per step
# than fp32 for V2-316M with no accuracy loss (docs/quantization/EXPERIMENT.md).
USE_BF16 = os.environ.get("FINETUNE_BF16", "").strip().lower() in ("1", "true", "yes")
if USE_BF16:
    print("[config] bf16 autocast enabled", flush=True)

# Activation memory scales with the batch size; on a 64 GB host (M4 Pro) the
# ~64 GiB peak measured at batch 8 + gradient checkpointing is the entire budget.
# Measured here: batch 4 still swapped (14 GB swap fully used, step time drifting
# 25 s -> 500 s), so batch 2 is the safe default on this machine.
FINETUNE_BATCH_SIZE = int(os.environ.get("FINETUNE_BATCH_SIZE", "2"))
print(f"[config] per_device_train_batch_size={FINETUNE_BATCH_SIZE}", flush=True)

# Cap the optimizer steps. A full epoch over the PD_atlas train split is ~960
# steps at batch 4; FINETUNE_MAX_STEPS>0 truncates it so a usable classifier can
# be produced inside a sane wall-clock budget on this host.
FINETUNE_MAX_STEPS = int(os.environ.get("FINETUNE_MAX_STEPS", "0"))
if FINETUNE_MAX_STEPS > 0:
    print(f"[config] max_steps={FINETUNE_MAX_STEPS}", flush=True)

# ------------------------------------------------------------------ classifier
from geneformer import Classifier

classifier = Classifier(
    classifier="cell",
    cell_state_dict={"state_key": "celltype", "states": CLASSES},
    filter_data=None,
    training_args={
        "num_train_epochs": 1,
        "learning_rate": 5e-5,
        # Batch-size/memory budget measured on this machine (M4 Max 128 GB,
        # analysis/measure_mps_batch.py): one train step at the worst-case
        # padded length of 4096 tokens peaks at ~145 GiB WITHOUT gradient
        # checkpointing (-> MPS OOM), but only ~64 GiB WITH it. Keep batch 8
        # and checkpoint instead of shrinking the effective batch.
        # On a 64 GB machine (M4 Pro) ~64 GiB is the whole budget, so the batch
        # size is overridable: FINETUNE_BATCH_SIZE=4 halves the activation peak.
        "per_device_train_batch_size": FINETUNE_BATCH_SIZE,
        "per_device_eval_batch_size": FINETUNE_BATCH_SIZE * 2,
        "gradient_checkpointing": True,
        "gradient_checkpointing_kwargs": {"use_reentrant": False},
        "seed": SEED,
        # Periodic step checkpoints so an interrupted run (OOM/crash/power)
        # resumes instead of losing the whole epoch. Override with
        # FINETUNE_SAVE_STEPS.
        # NOTE: the stock classifier hardcodes save_strategy="epoch" through
        # def_training_args.update(); patches/checkpoints/classifier.py switches
        # that to setdefault, which is what makes these lines take effect.
        "save_strategy": "steps",
        "save_steps": int(os.environ.get("FINETUNE_SAVE_STEPS", "200")),
        "save_total_limit": 2,
        "logging_steps": 100,
        "report_to": "none",
        # max_steps ends the single epoch at the cap. The in-training eval over
        # the ~3,800-cell eval cohort is skipped when a cap is set; the real
        # held-out evaluation runs afterwards via evaluate_saved_model().
        **({"max_steps": FINETUNE_MAX_STEPS,
            "evaluation_strategy": "no",
            "load_best_model_at_end": False} if FINETUNE_MAX_STEPS > 0 else {}),
        **({"bf16": True} if USE_BF16 else {}),
    },
    max_ncells=None,
    freeze_layers=6,
    num_crossval_splits=1,
    forward_batch_size=8,
    nproc=8,
    model_version="V2",
)

prepared_train = RUN_DIR / f"{RUN_PREFIX}_labeled_train.dataset"
prepared_test = RUN_DIR / f"{RUN_PREFIX}_labeled_test.dataset"
id_class_file = RUN_DIR / f"{RUN_PREFIX}_id_class_dict.pkl"

if not (prepared_train.exists() and prepared_test.exists() and id_class_file.exists()):
    print("preparing data...", flush=True)
    classifier.prepare_data(
        input_data_file=str(TOKENIZED),
        output_directory=str(RUN_DIR),
        output_prefix=RUN_PREFIX,
        split_id_dict={"attr_key": "individual", "train": TRAIN + EVAL, "test": TEST},
    )
else:
    print("reusing prepared datasets", flush=True)

# ------------------------------------------------------------------- training
checkpoints = sorted(RUN_DIR.glob(f"*geneformer_cellClassifier_{RUN_PREFIX}/ksplit1"))


def _valid_checkpoint(path) -> bool:
    """A checkpoint dir is reusable if it holds weights directly, OR holds
    step checkpoints (``checkpoint-<N>/``) that do.

    With ``save_strategy="steps"`` (the new default) an interrupted run leaves
    its weights only inside ``checkpoint-<N>/`` subdirectories, so checking the
    top level alone would classify such a run as incomplete and DELETE hours of
    training. Never delete a directory that contains resumable step
    checkpoints.
    """
    if not path.is_dir():
        return False
    if any(
        path.joinpath(name).exists()
        for name in ("pytorch_model.bin", "model.safetensors", "adapter_model.bin")
    ):
        return True
    for sub in path.glob("checkpoint-*"):
        if any(
            sub.joinpath(name).exists()
            for name in ("pytorch_model.bin", "model.safetensors", "adapter_model.bin")
        ):
            return True
    return False


# Drop checkpoint dirs that are empty/incomplete so training re-runs fresh.
for stale in sorted(RUN_DIR.glob(f"*geneformer_cellClassifier_{RUN_PREFIX}/ksplit*")):
    if not _valid_checkpoint(stale):
        print(f"removing incomplete checkpoint dir: {stale}", flush=True)
        shutil.rmtree(stale, ignore_errors=True)
checkpoints = [c for c in checkpoints if _valid_checkpoint(c)]

if not checkpoints:
    print("fine-tuning (retune, 1 epoch, first 6 layers frozen)...", flush=True)
    classifier.validate(
        model_directory=str(MODEL_DIR),
        prepared_input_data_file=str(prepared_train),
        id_class_dict_file=str(id_class_file),
        output_directory=str(RUN_DIR),
        output_prefix=RUN_PREFIX,
        split_id_dict={"attr_key": "individual", "train": TRAIN, "eval": EVAL},
        n_hyperopt_trials=0,
    )
    checkpoints = sorted(RUN_DIR.glob(f"*geneformer_cellClassifier_{RUN_PREFIX}/ksplit1"))
    checkpoints = [c for c in checkpoints if _valid_checkpoint(c)]
    assert checkpoints, "training finished without a valid (non-empty) ksplit1 checkpoint"
else:
    print("reusing existing checkpoint", flush=True)

TRAINED = checkpoints[-1]
(RUN_DIR / "TRAINED_MODEL_PATH.txt").write_text(str(TRAINED))
print("checkpoint:", TRAINED, flush=True)

# ----------------------------------------------------------------- evaluation
print("evaluating on held-out replicate _3 samples...", flush=True)
classifier.evaluate_saved_model(
    model_directory=str(TRAINED),
    id_class_dict_file=str(id_class_file),
    test_data_file=str(prepared_test),
    output_directory=str(RUN_DIR),
    output_prefix=f"{RUN_PREFIX}_heldout",
    predict=True,
    predict_metadata=["cell_id", "individual"],
)

id_to_class = pickle.load(open(id_class_file, "rb"))
payload = pickle.load(open(RUN_DIR / f"{RUN_PREFIX}_heldout_pred_dict.pkl", "rb"))
y_true = np.array([id_to_class[i] for i in payload["label_ids"]])
y_pred = np.array([id_to_class[i] for i in payload["pred_ids"]])
print(f"predicted {len(y_true)} held-out cells", flush=True)

report = pd.DataFrame(
    classification_report(y_true, y_pred, labels=CLASSES, output_dict=True, zero_division=0)
).T
report.to_csv(TABLE_DIR / "adpd_finetuned_classification_report.csv")
print(report.round(3).to_string(), flush=True)

acc = accuracy_score(y_true, y_pred)
mf1 = f1_score(y_true, y_pred, average="macro")
print(f"fine-tuned: accuracy {acc:.4f}  macro F1 {mf1:.4f}", flush=True)

cm = pd.DataFrame(
    confusion_matrix(y_true, y_pred, labels=CLASSES, normalize="true"),
    index=CLASSES, columns=CLASSES,
)
cm.to_csv(TABLE_DIR / "adpd_finetuned_confusion_matrix_normalized.csv")

plt.figure(figsize=(13, 10.5))
sns.heatmap(cm, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=1,
            annot_kws={"size": 6.5}, cbar_kws={"label": "fraction of true class"})
plt.title((PREFIX + " - fine-tuned Geneformer cell classifier\n"
          "held-out samples (replicate _3), row-normalized"))
plt.xlabel("Predicted cell type")
plt.ylabel("True cell type")
plt.tight_layout()
plt.savefig(FIG_DIR / "adpd_finetuned_confusion_matrix_normalized.png", dpi=180,
            bbox_inches="tight")
plt.close()

pd.DataFrame({"cell_id": payload["prediction_metadata"]["cell_id"],
              "individual": payload["prediction_metadata"]["individual"],
              "true": y_true, "predicted": y_pred}
             ).to_csv(TABLE_DIR / "adpd_finetuned_test_predictions.csv", index=False)

# ------------------------------------------------------------- comparison row
baseline = json.loads((TABLE_DIR / "adpd_baseline_summary.json").read_text())
comparison = pd.DataFrame([
    {"method": "Frozen embeddings + logistic regression",
     "accuracy": baseline["accuracy"], "macro_f1": baseline["macro_f1"]},
    {"method": "Fine-tuned Geneformer classifier",
     "accuracy": float(acc), "macro_f1": float(mf1)},
])
comparison.to_csv(TABLE_DIR / "adpd_model_comparison.csv", index=False)
print(comparison.round(4).to_string(index=False), flush=True)

# Record the ACTUAL training regime. A truncated run (FINETUNE_MAX_STEPS>0) is a
# PARTIAL epoch and must never be reported as "1 epoch": that would hide a
# parameter choice that directly affects accuracy.
if FINETUNE_MAX_STEPS > 0:
    print(f"[WARN] TRAINING TRUNCATED: max_steps={FINETUNE_MAX_STEPS} is a partial "
          "epoch -- accuracy may be below a full-epoch run. Report this with the numbers.",
          flush=True)

(TABLE_DIR / "adpd_finetuned_summary.json").write_text(json.dumps({
    "dataset": f"{PREFIX} (mouse -> human ortholog mapped)",
    "method": f"fine-tuned {MODEL_NAME} cell classifier",
    "model_name": MODEL_NAME,
    "learning_rate": 5e-5, "freeze_layers": 6,
    "per_device_train_batch_size": FINETUNE_BATCH_SIZE,
    "bf16_autocast": bool(USE_BF16),
    "gradient_checkpointing": True,
    "num_train_epochs_requested": 1,
    "max_steps_cap": (FINETUNE_MAX_STEPS or None),
    "training_truncated": bool(FINETUNE_MAX_STEPS > 0),
    "n_train_samples": len(TRAIN), "n_eval_samples": len(EVAL), "n_test_samples": len(TEST),
    "n_test_cells": int(len(y_true)), "n_classes": len(CLASSES),
    "accuracy": float(acc), "macro_f1": float(mf1),
    "baseline_accuracy": baseline["accuracy"], "baseline_macro_f1": baseline["macro_f1"],
    "checkpoint": str(TRAINED),
}, indent=2))
print("=== DONE ===", flush=True)
