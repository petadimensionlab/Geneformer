#!/usr/bin/env python
"""Extract frozen Geneformer cell embeddings for a tissue, then stop.

This is the "embeddings only" stage that `05_figures.py` depends on. It mirrors
`04_baseline.py`'s tokenize + EmbExtractor steps but does **not** fit the
logistic probe (05 refits it itself). On a fresh tissue it both tokenizes and
extracts; if embeddings already exist it reuses them and exits.

Output (same path 05 reads):
    <input>/<TISSUE>/results/embeddings/pretrained_cell_embeddings.csv

Run per tissue:
    env ADPD_TISSUE=PD_brain .venv/bin/python analysis/04b_extract_embeddings.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import find_h5ad, resolve

os.environ["WANDB_DISABLED"] = "true"

ROOT, PREFIX = resolve()
GF_ROOT = Path(os.environ["GENEFORMER_DIR"])
MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
MODEL_DIR = GF_ROOT / MODEL_NAME
if not MODEL_DIR.is_dir():
    raise SystemExit(
        f"Geneformer checkpoint not found: {MODEL_DIR}\n"
        "Set GENEFORMER_DIR (and GENEFORMER_MODEL) or run ./download.sh."
    )

TOKEN_DIR = ROOT / "tokenized"
EMB_DIR = ROOT / "results/embeddings"
EMB_DIR.mkdir(parents=True, exist_ok=True)
TOKENIZED = TOKEN_DIR / f"{PREFIX}.dataset"
EMB_FILE = EMB_DIR / "pretrained_cell_embeddings.csv"
LABEL_COLUMNS = ["cell_id", "individual", "celltype", "split"]

if EMB_FILE.exists():
    print("embeddings already exist:", EMB_FILE, flush=True)
    raise SystemExit(0)


def ensure_tokenization():
    if TOKENIZED.exists():
        print("tokenized dataset exists:", TOKENIZED, flush=True)
        return
    import anndata as _ad

    h5ad = find_h5ad(ROOT)
    obs = list(_ad.read_h5ad(h5ad, backed="r").obs.columns)
    missing = [c for c in LABEL_COLUMNS if c not in obs]
    if missing:
        raise SystemExit(f"{h5ad.name} missing required obs cols: {missing}")
    attr = LABEL_COLUMNS + [c for c in obs
                            if c not in LABEL_COLUMNS + ["n_counts", "filter_pass"]]
    from geneformer import TranscriptomeTokenizer

    print("tokenizing...", flush=True)
    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={c: c for c in attr},
        nproc=8,
        chunk_size=512,
        model_version="V2",
        use_h5ad_index=True,
    )
    tk.tokenize_data(
        data_directory=str(ROOT / "h5ad"),
        output_directory=str(TOKEN_DIR),
        output_prefix=PREFIX,
        file_format="h5ad",
    )


ensure_tokenization()

from geneformer import EmbExtractor  # noqa: E402
from geneformer.device import get_device  # noqa: E402

print("embedding extraction on device:", get_device(), flush=True)
ex = EmbExtractor(
    model_type="Pretrained",
    num_classes=0,
    emb_mode="cell",
    max_ncells=None,
    emb_layer=-1,
    emb_label=LABEL_COLUMNS,
    labels_to_plot=["celltype", "split"],
    forward_batch_size=4,
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
print("=== EMBEDDINGS DONE ===", flush=True)