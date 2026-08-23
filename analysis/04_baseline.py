#!/usr/bin/env python
"""Run the tutorial's frozen-embedding baseline on the converted ADPD data.

Mirrors 02_lung_allograft_classification_tutorial.ipynb section 5: tokenize,
extract frozen pretrained Geneformer cell embeddings, fit a logistic-regression
probe on train samples, evaluate on held-out samples. No fine-tuning.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import resolve as _resolve_tissue

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, f1_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 42
np.random.seed(SEED)
os.environ["WANDB_DISABLED"] = "true"

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
INPUT_DIR = ROOT / "h5ad"
TOKEN_DIR = ROOT / "tokenized"
RESULT_DIR = ROOT / "results"
EMB_DIR = RESULT_DIR / "embeddings"
TABLE_DIR = RESULT_DIR / "tables"
for d in (TOKEN_DIR, EMB_DIR, TABLE_DIR):
    d.mkdir(parents=True, exist_ok=True)

TOKENIZED = TOKEN_DIR / f"{PREFIX}.dataset"
LABEL_COLUMNS = ["cell_id", "individual", "celltype", "split"]

# ------------------------------------------------------------------ tokenize
import anndata as _ad

_h5ad = next((p for p in INPUT_DIR.glob("*.h5ad") if not p.name.startswith(".")), None)
if _h5ad is None:
    raise SystemExit(f"no .h5ad found in {INPUT_DIR}; run 03_map_and_build_h5ad.py first")
_obs_cols = list(_ad.read_h5ad(_h5ad, backed="r").obs.columns)
_required = ["cell_id", "individual", "celltype", "split"]
_missing = [c for c in _required if c not in _obs_cols]
if _missing:
    raise SystemExit(f"{_h5ad.name} is missing required obs columns: {_missing}")
_attr_columns = _required + [c for c in _obs_cols
                             if c not in _required + ["n_counts", "filter_pass"]]
print("  tokenizer will carry:", ", ".join(_attr_columns), flush=True)

from geneformer import TranscriptomeTokenizer

if not TOKENIZED.exists():
    print("tokenizing...", flush=True)
    tk = TranscriptomeTokenizer(
        # carry every label column the .h5ad actually has; extra study-specific
        # columns (disease, timepoint, ...) come along when present
        custom_attr_name_dict={c: c for c in _attr_columns},
        nproc=8,
        chunk_size=512,
        model_version="V2",
        use_h5ad_index=True,
    )
    tk.tokenize_data(
        data_directory=str(INPUT_DIR),
        output_directory=str(TOKEN_DIR),
        output_prefix=PREFIX,
        file_format="h5ad",
    )
else:
    print("tokenized dataset exists:", TOKENIZED, flush=True)

from datasets import load_from_disk

tokens = load_from_disk(str(TOKENIZED))
print(tokens, flush=True)
print(pd.Series(tokens["split"]).value_counts(), flush=True)
print("median tokens/cell:", int(np.median(tokens["length"])), flush=True)

# -------------------------------------------------------- frozen embeddings
from geneformer import EmbExtractor

EMB_FILE = EMB_DIR / "pretrained_cell_embeddings.csv"
if EMB_FILE.exists():
    emb = pd.read_csv(EMB_FILE, index_col=0)
    print("reusing embeddings:", emb.shape, flush=True)
else:
    import torch

    from geneformer.device import get_device

    device = get_device()
    print(f"embedding extraction on device: {device}", flush=True)
    ex = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,
        emb_mode="cell",
        max_ncells=None,
        emb_layer=-1,
        emb_label=LABEL_COLUMNS,
        labels_to_plot=["celltype", "split"],
        forward_batch_size=32,
        nproc=8,
        model_version="V2",
    )
    emb = ex.extract_embs(
        model_directory=str(MODEL_DIR),
        input_data_file=str(TOKENIZED),
        output_directory=str(EMB_DIR),
        output_prefix="pretrained_cell_embeddings",
    )
    print("embeddings:", emb.shape, flush=True)

feat = [c for c in emb.columns if c not in LABEL_COLUMNS]

# ------------------------------------------------------------------- probe
train_mask = emb["split"].eq("train")
test_mask = emb["split"].eq("test")
print(f"train cells {train_mask.sum()}, test cells {test_mask.sum()}", flush=True)

probe = make_pipeline(
    StandardScaler(),
    LogisticRegression(max_iter=2_000, class_weight="balanced", random_state=SEED),
)
probe.fit(emb.loc[train_mask, feat], emb.loc[train_mask, "celltype"])

y_true = emb.loc[test_mask, "celltype"].to_numpy()
y_pred = probe.predict(emb.loc[test_mask, feat])

CLASSES = sorted(pd.unique(emb["celltype"]))
report = pd.DataFrame(
    classification_report(y_true, y_pred, labels=CLASSES, output_dict=True, zero_division=0)
).T
report.to_csv(TABLE_DIR / "adpd_baseline_classification_report.csv")
print(report.round(3).to_string(), flush=True)

cm = pd.DataFrame(
    confusion_matrix(y_true, y_pred, labels=CLASSES, normalize="true"),
    index=CLASSES,
    columns=CLASSES,
)
cm.to_csv(TABLE_DIR / "adpd_baseline_confusion_matrix_normalized.csv")

summary = {
    "dataset": f"{PREFIX} (mouse -> human ortholog mapped)",
    "n_cells_total": int(len(emb)),
    "n_train": int(train_mask.sum()),
    "n_test": int(test_mask.sum()),
    "n_classes": len(CLASSES),
    "accuracy": float(accuracy_score(y_true, y_pred)),
    "macro_f1": float(f1_score(y_true, y_pred, average="macro")),
}
(TABLE_DIR / "adpd_baseline_summary.json").write_text(json.dumps(summary, indent=2))
print(json.dumps(summary, indent=2), flush=True)
print("=== DONE ===", flush=True)
