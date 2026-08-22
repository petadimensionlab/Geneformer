#!/usr/bin/env python
"""Figures for the ADPD baseline, mirroring the tutorial's section 7-8 plots.

Reuses the cached embeddings, so no GPU work: refits the same probe, then draws
the normalized confusion heatmap and the held-out-cell UMAP.
"""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

SEED = 42
np.random.seed(SEED)

ROOT = Path(os.environ.get("ADPD_ROOT", Path.home() / "Documents/ADPD_analysis"))
PREFIX = os.environ.get("ADPD_PREFIX", "dataset")
EMB_FILE = ROOT / "results/embeddings/pretrained_cell_embeddings.csv"
FIG_DIR = ROOT / "results/figures"
TABLE_DIR = ROOT / "results/tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)

LABEL_COLUMNS = ["cell_id", "individual", "celltype", "split"]

print("loading embeddings...", flush=True)
emb = pd.read_csv(EMB_FILE, index_col=0)
feat = [c for c in emb.columns if c not in LABEL_COLUMNS]
print(f"  {emb.shape[0]} cells x {len(feat)} dims", flush=True)

train_mask = emb["split"].eq("train")
test_mask = emb["split"].eq("test")

print("refitting probe (same params as 04)...", flush=True)
probe = make_pipeline(
    StandardScaler(),
    LogisticRegression(max_iter=2_000, class_weight="balanced", random_state=SEED),
)
probe.fit(emb.loc[train_mask, feat], emb.loc[train_mask, "celltype"])

test = emb.loc[test_mask].reset_index(drop=True)
pred = probe.predict(test[feat])
CLASSES = sorted(pd.unique(emb["celltype"]))

pd.DataFrame(
    {"cell_id": test["cell_id"], "individual": test["individual"],
     "true": test["celltype"], "predicted": pred}
).to_csv(TABLE_DIR / "adpd_test_predictions.csv", index=False)
print("saved predictions", flush=True)

# ------------------------------------------------------ confusion heatmap
cm = pd.read_csv(TABLE_DIR / "adpd_baseline_confusion_matrix_normalized.csv", index_col=0)
plt.figure(figsize=(13, 10.5))
sns.heatmap(
    cm, annot=True, fmt=".2f", cmap="Blues", vmin=0, vmax=1,
    annot_kws={"size": 6.5}, cbar_kws={"label": "fraction of true class"},
)
plt.title("ADPD PD Blood — frozen Geneformer embeddings + logistic regression\n"
          "held-out samples (replicate _3), row-normalized")
plt.xlabel("Predicted cell type")
plt.ylabel("True cell type")
plt.tight_layout()
plt.savefig(FIG_DIR / "adpd_baseline_confusion_matrix_normalized.png", dpi=180,
            bbox_inches="tight")
plt.close()
print("wrote confusion heatmap", flush=True)

# ------------------------------------------------------------------ UMAP
print("computing UMAP over held-out cells...", flush=True)
view = ad.AnnData(test[feat].to_numpy(dtype=np.float32))
view.obs = test[LABEL_COLUMNS].copy()
view.obs["prediction"] = pred
view.obs["correct"] = np.where(view.obs["celltype"].to_numpy() == pred, "correct", "misclassified")

sc.pp.neighbors(view, n_neighbors=15, use_rep="X", random_state=SEED)
sc.tl.umap(view, random_state=SEED)

sc.pl.umap(view, color=["celltype", "prediction"], wspace=0.45, show=False)
plt.suptitle(f"Pretrained Geneformer embeddings: {PREFIX} held-out samples", y=1.02)
plt.savefig(FIG_DIR / "adpd_baseline_test_embedding_umap.png", dpi=180, bbox_inches="tight")
plt.close()
print("wrote celltype/prediction UMAP", flush=True)

# where the errors sit -- not in the tutorial, but the interesting part here
sc.pl.umap(view, color=["correct"], palette={"correct": "#d9d9d9", "misclassified": "#d62728"},
           show=False)
plt.suptitle("Misclassified held-out cells", y=1.02)
plt.savefig(FIG_DIR / "adpd_baseline_errors_umap.png", dpi=180, bbox_inches="tight")
plt.close()
print("wrote error UMAP", flush=True)

view.write_h5ad(ROOT / "results/embeddings/adpd_test_umap.h5ad")
print("=== DONE ===", flush=True)
