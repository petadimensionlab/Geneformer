#!/usr/bin/env python
"""In silico perturbation — COMBINATION (double-deletion) mode on MPS for
PD_brain, MG interplay: delete SNCA together with a Parkinson's-risk partner
in the same cell (a "double knockout") and measure the combined goal-state shift.

Context / design:
- Single-gene deletions (07b) delete ONE gene. Here we delete SNCA AND a
  partner simultaneously in each cell and compare the combined tg->WT shift.
- Geneformer constraint: when `genes_to_perturb` is a LIST, `combos` is forced
  to 0 and `anchor_gene` is disabled (with cell_states_to_model set), so the
  whole list is deleted together — i.e. a true double knockout.
- Co-expression limit: the deletion is only evaluated on cells that express
  BOTH genes. Measured for tg Microglia: SNCA+LRRK2 = 3 cells (too few),
  SNCA+PINK1 = 47 cells. So only adequately co-expressed pairs are run;
  pairs below `MIN_COEXP_CELLS` are skipped with a note.

State modelling (same as 07b):
    state_key="disease", start="tg" (PD model), goal="WT" (healthy),
    alt=["mo","PF"], restricted to a target cell type (Microglia default).
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
    """Warn loudly if datasets>=5 (known perturb_data map hang)."""
    import datasets

    major = int(datasets.__version__.split(".")[0])
    if major >= 5:
        print(
            f"[WARN] datasets=={datasets.__version__} is installed. "
            "InSilicoPerturber.perturb_data hangs on dataset.map with datasets>=5. "
            "Pin datasets==4.0.0.",
            file=sys.stderr,
            flush=True,
        )


def estimate_perturb_ram(n_cells: int, n_genes: int, combos: int, seq_len: int = 4096) -> float:
    n_variants = max(n_genes, 1)
    total_variants = n_cells * n_variants
    return (total_variants * seq_len * 10) / 1e9


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
    print("no fine-tuned classifier found; using pretrained model", flush=True)

ISP_DIR = ROOT / "results" / "isp"
ISP_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- states
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
print(f"modeled states: key={CELL_STATE_KEY} start={START_STATE} goal={GOAL_STATE} "
      f"alt={ALT_STATES}", flush=True)
print(f"target celltype: {TARGET_CELLTYPE}", flush=True)

# ---------------------------------------------------------------- step 1: state embs
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
    cell_states_to_model,
    CELLCLASSIFIER_DIR,
    str(TOKENIZED),
    str(ISP_DIR),
    "isp_state_embs",
)
print("state_embs keys:", list(state_embs_dict.keys()), flush=True)

# ---------------------------------------------------------------- combinations
SYMBOL = {
    "ENSG00000145335": "SNCA",
    "ENSG00000188906": "LRRK2",
    "ENSG00000177628": "GBA1",
    "ENSG00000158828": "PINK1",
    "ENSG00000116288": "PARK7",
    "ENSG00000185345": "PRKN",
    "ENSG00000154277": "UCHL1",
    "ENSG00000186868": "MAPT",
}

ANCHOR = "ENSG00000145335"  # SNCA
PARTNERS = [
    "ENSG00000188906",  # LRRK2
    "ENSG00000158828",  # PINK1
]

# skip pairs whose target cell type does not co-express both genes well enough
MIN_COEXP_CELLS = int(os.environ.get("ISP_MIN_COEXP", "5"))

MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "200"))


def _coexp_tg_count(gene_a: str, gene_b: str) -> int:
    """Number of start-state cells of the target cell type expressing both genes."""
    import pickle

    td = pickle.load(open(str(GF_ROOT / "geneformer" / "token_dictionary_gc104M.pkl"), "rb"))
    tok_a, tok_b = td[gene_a], td[gene_b]
    n = 0
    for ex in tokens:
        if ex["celltype"] != TARGET_CELLTYPE or ex[CELL_STATE_KEY] != START_STATE:
            continue
        if {tok_a, tok_b}.issubset(set(ex["input_ids"])):
            n += 1
    return n


for partner in PARTNERS:
    symbol_a, symbol_b = SYMBOL[ANCHOR], SYMBOL[partner]
    pair = [ANCHOR, partner]
    name = f"{symbol_a}+{symbol_b}"

    n_coexp = _coexp_tg_count(ANCHOR, partner)
    print(f"\n===== combo {name}: co-expressing {START_STATE} {TARGET_CELLTYPE} cells = {n_coexp} =====", flush=True)
    if n_coexp < MIN_COEXP_CELLS:
        print(f"[SKIP] {name}: only {n_coexp} co-expressing cells (below {MIN_COEXP_CELLS}).", flush=True)
        continue

    combo_isp_dir = ISP_DIR / f"combo_{name}"
    combo_isp_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"isp_combo_{name}"

    _ram_est = estimate_perturb_ram(MAX_CELLS, len(pair), 0)
    print(
        f"perturb_data: max_ncells={MAX_CELLS}, genes={pair}, "
        f"estimated dataset RAM ~{_ram_est:.2f} GB",
        flush=True,
    )

    isp = InSilicoPerturber(
        perturb_type="delete",
        perturb_rank_shift=None,
        genes_to_perturb=pair,
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

    try:
        isp.perturb_data(
            CELLCLASSIFIER_DIR,
            str(TOKENIZED),
            str(combo_isp_dir),
            prefix,
        )
    except Exception as exc:
        print(f"[SKIP] {name}: no perturbable cells ({type(exc).__name__})", flush=True)
        continue

    ispstats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed=pair,
        combos=0,
        anchor_gene=None,
        cell_states_to_model=cell_states_to_model,
        model_version="V2",
    )
    ispstats.get_stats(
        str(combo_isp_dir),
        None,
        str(combo_isp_dir),
        prefix,
    )
    print(f"===== done {name} =====", flush=True)

print("=== IS COMBINATION PERTURBATION DONE ===", flush=True)