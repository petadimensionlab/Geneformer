#!/usr/bin/env python
"""Verify HF Trainer-based fine-tuning works on MPS with Geneformer V2-104M.

Trains a tiny CellClassifier for a couple of steps on a small downsampled
subset to prove the Trainer + Classifier path runs on MPS without CUDA. Uses
the geneformer Classifier (which uses HF Trainer) directly.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ["WANDB_DISABLED"] = "true"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import resolve as _resolve_tissue

ROOT, PREFIX = _resolve_tissue()
GF_ROOT = Path(os.environ["GENEFORMER_DIR"])
MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
MODEL_DIR = GF_ROOT / MODEL_NAME
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"
RUN_DIR = ROOT / "runs"
RUN_DIR.mkdir(parents=True, exist_ok=True)

from geneformer.device import get_device
from geneformer import Classifier

device = get_device()
print(f"device: {device}", flush=True)
assert device == "mps", f"expected mps, got {device}"

from datasets import load_from_disk
tokens = load_from_disk(str(TOKENIZED))

# Use a small subset of individuals to keep the run fast on MPS
samples = sorted(set(tokens["individual"]))
train_samples = samples[:4]
eval_samples = samples[4:8]
test_samples = samples[8:10]
CLASSES = sorted(set(tokens["celltype"]))
print(f"{len(CLASSES)} classes; samples train={len(train_samples)} eval={len(eval_samples)}", flush=True)

classifier = Classifier(
    classifier="cell",
    cell_state_dict={"state_key": "celltype", "states": CLASSES},
    filter_data=None,
    training_args={
        "num_train_epochs": 1,
        "learning_rate": 5e-5,
        "per_device_train_batch_size": 8,
        "per_device_eval_batch_size": 16,
        "seed": 42,
        "save_strategy": "epoch",
        "logging_steps": 10,
        "report_to": "none",
        "max_steps": 8,  # tiny, just to prove MPS training works
    },
    max_ncells=None,
    freeze_layers=6,
    num_crossval_splits=1,
    forward_batch_size=16,
    nproc=8,
    model_version="V2",
)

prepared_train = RUN_DIR / "mpsft_labeled_train.dataset"
prepared_test = RUN_DIR / "mpsft_labeled_test.dataset"
id_class_file = RUN_DIR / "mpsft_id_class_dict.pkl"

classifier.prepare_data(
    input_data_file=str(TOKENIZED),
    output_directory=str(RUN_DIR),
    output_prefix="mpsft",
    split_id_dict={"attr_key": "individual", "train": train_samples + eval_samples, "test": test_samples},
)

classifier.validate(
    model_directory=str(MODEL_DIR),
    prepared_input_data_file=str(prepared_train),
    id_class_dict_file=str(id_class_file),
    output_directory=str(RUN_DIR),
    output_prefix="mpsft",
    split_id_dict={"attr_key": "individual", "train": train_samples, "eval": eval_samples},
    n_hyperopt_trials=0,
)
print("=== MPS FINE-TUNE DONE ===", flush=True)
