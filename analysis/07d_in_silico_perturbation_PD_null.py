#!/usr/bin/env python
"""In silico perturbation — goal_state_shift with NULL SIGNIFICANCE (for
parkinson PD_brain, Microglia). Based on 07b but with `genes_to_perturb="all"`
so the stats step can compare each gene's shift against random (null) shifts.

Configuration: Microglia, start=tg (PD model), goal=WT, alt=[mo,PF], delete mode.
Runs compactly (small max_ncells) because "all" enumerates every deleted gene
per cell -> many perturbed variants per cell (max_ncells * genes-per-cell).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ["WANDB_DISABLED"] = "true"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import resolve as _resolve_tissue

from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats


def _check_datasets_version() -> None:
    import datasets

    major = int(datasets.__version__.split(".")[0])
    if major >= 5:
        print("[WARN] datasets>=5 hangs perturb_data; pin datasets==4.0.0.", file=sys.stderr, flush=True)


NPROC = int(os.environ.get("IS_NPROC", "1"))
_check_datasets_version()

ROOT, PREFIX = _resolve_tissue()
GF_ROOT = Path(os.environ["GENEFORMER_DIR"])
MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

RUN_DIR = ROOT / "runs"
trained_path_file = RUN_DIR / "TRAINED_MODEL_PATH.txt"
if trained_path_file.exists():
    CELLCLASSIFIER_DIR = trained_path_file.read_text().strip()
    print("using fine-tuned cell classifier:", CELLCLASSIFIER_DIR, flush=True)
else:
    CELLCLASSIFIER_DIR = str(GF_ROOT / MODEL_NAME)

ISP_DIR = ROOT / "results" / "isp"
ISP_DIR.mkdir(parents=True, exist_ok=True)

CELL_STATE_KEY = os.environ.get("ISP_STATE_KEY", "disease")
START_STATE = os.environ.get("ISP_START_STATE", "tg")
GOAL_STATE = os.environ.get("ISP_GOAL_STATE", "WT")
ALT_STATES = [s for s in os.environ.get("ISP_ALT_STATES", "mo,PF").split(",") if s]
TARGET_CELLTYPE = os.environ.get("ISP_CELLTYPE", "Microglia")

cell_states_to_model = {
    "state_key": CELL_STATE_KEY,
    "start_state": START_STATE,
    "goal_state": GOAL_STATE,
    "alt_states": ALT_STATES,
}
filter_data_dict = {"celltype": [TARGET_CELLTYPE]}

from datasets import load_from_disk

tokens = load_from_disk(str(TOKENIZED))
celltypes = sorted(set(tokens["celltype"]))
print(f"{len(celltypes)} cell types available", flush=True)
print(f"states: start={START_STATE} goal={GOAL_STATE} alt={ALT_STATES} | celltype={TARGET_CELLTYPE}", flush=True)

MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "50"))
ISP_ODIR = ISP_DIR / f"{TARGET_CELLTYPE}" / "null_all"
ISP_ODIR.mkdir(parents=True, exist_ok=True)

# step 1: state embeddings (reuse same params; not "all")
embex = EmbExtractor(
    model_type="CellClassifier",
    num_classes=len(celltypes),
    filter_data=filter_data_dict,
    max_ncells=1000,
    emb_layer=0,
    summary_stat="exact_mean",
    forward_batch_size=64,
    model_version="V2",
    nproc=NPROC,
)
state_embs_dict = embex.get_state_embs(
    cell_states_to_model, CELLCLASSIFIER_DIR, str(TOKENIZED), str(ISP_DIR), "isp_state_embs",
)
print("state_embs keys:", list(state_embs_dict.keys()), flush=True)

# step 2: perturb ALL genes (necessary for null comparison)
isp = InSilicoPerturber(
    perturb_type="delete",
    perturb_rank_shift=None,
    genes_to_perturb="all",
    combos=0,
    anchor_gene=None,
    model_type="CellClassifier",
    num_classes=len(celltypes),
    emb_mode="cls",
    cell_emb_style="mean_pool",
    filter_data=filter_data_dict,
    cell_states_to_model=cell_states_to_model,
    state_embs_dict=state_embs_dict,
    max_ncells=MAX_CELLS,
    emb_layer=0,
    forward_batch_size=64,
    model_version="V2",
    nproc=NPROC,
)
print(f"perturb_data (null/all): max_ncells={MAX_CELLS} ...", flush=True)
isp.perturb_data(CELLCLASSIFIER_DIR, str(TOKENIZED), str(ISP_ODIR), "isp_perturbation_all")

# step 3: null comparison stats
ispstats = InSilicoPerturberStats(
    mode="goal_state_shift",
    genes_perturbed="all",
    combos=0,
    anchor_gene=None,
    cell_states_to_model=cell_states_to_model,
    model_version="V2",
)
ispstats.get_stats(str(ISP_ODIR), None, str(ISP_ODIR), "isp_stats_null")
print("=== IS NULL PERTURBATION DONE ===", flush=True)