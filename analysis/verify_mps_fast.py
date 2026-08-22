#!/usr/bin/env python
"""Fast MPS verification of the frozen-embedding + probe pipeline.

Uses a downsampled subset (max_ncells) so the run completes in minutes on a
Mac MPS instead of the multi-hour full-dataset pass. Proves the tokenizer ->
EmbExtractor -> probe path is MPS-correct end to end.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, f1_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

os.environ["WANDB_DISABLED"] = "true"
SEED = 42
np.random.seed(SEED)

ROOT = Path(os.environ.get("ADPD_ROOT", Path.cwd() / "analysis_ws"))
PREFIX = os.environ.get("ADPD_PREFIX", "PD_smallint")
GF_ROOT = Path(os.environ["GENEFORMER_DIR"])
MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
MODEL_DIR = GF_ROOT / MODEL_NAME
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"
OUT = ROOT / "results" / "mps_verify"
OUT.mkdir(parents=True, exist_ok=True)

LABEL_COLUMNS = ["cell_id", "individual", "celltype", "split"]

from geneformer.device import get_device
from geneformer import EmbExtractor

device = get_device()
print(f"device: {device}", flush=True)
assert device == "mps", f"expected mps, got {device}"

MAX_CELLS = int(os.environ.get("MPS_VERIFY_CELLS", "2000"))

ex = EmbExtractor(
    model_type="Pretrained",
    num_classes=0,
    emb_mode="cell",
    max_ncells=MAX_CELLS,
    emb_layer=-1,
    emb_label=LABEL_COLUMNS,
    labels_to_plot=["celltype"],
    forward_batch_size=16,
    nproc=8,
    model_version="V2",
)
emb = ex.extract_embs(
    model_directory=str(MODEL_DIR),
    input_data_file=str(TOKENIZED),
    output_directory=str(OUT),
    output_prefix="mps_verify_cell_embeddings",
)
print(f"embeddings extracted: {emb.shape}", flush=True)

feat = [c for c in emb.columns if c not in LABEL_COLUMNS]
train_mask = emb["split"].eq("train")
test_mask = emb["split"].eq("test")
print(f"train {train_mask.sum()}, test {test_mask.sum()}", flush=True)

probe = make_pipeline(
    StandardScaler(),
    LogisticRegression(max_iter=1000, class_weight="balanced", random_state=SEED),
)
probe.fit(emb.loc[train_mask, feat], emb.loc[train_mask, "celltype"])
y_true = emb.loc[test_mask, "celltype"].to_numpy()
y_pred = probe.predict(emb.loc[test_mask, feat])
acc = accuracy_score(y_true, y_pred)
mf1 = f1_score(y_true, y_pred, average="macro")
print(f"frozen probe accuracy {acc:.4f}  macro F1 {mf1:.4f}", flush=True)

summary = {
    "device": device,
    "model": MODEL_NAME,
    "n_cells": int(len(emb)),
    "n_train": int(train_mask.sum()),
    "n_test": int(test_mask.sum()),
    "accuracy": float(acc),
    "macro_f1": float(mf1),
}
(OUT / "mps_verify_summary.json").write_text(
    __import__("json").dumps(summary, indent=2)
)
print("=== MPS VERIFY DONE ===", flush=True)
