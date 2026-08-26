#!/usr/bin/env python
"""In silico perturbation (deletion mode) on MPS — Parkinson's-disease-specific
configuration for PD_brain, built on the same engine as
`analysis/07_in_silico_perturbation.py`.

Key adaptations vs. 07:
- `07` picks the two most-abundant `celltype`s as start/goal (disease-agnostic).
  Here `state_key="disease"` so the start state is the PD-model group and the
  goal is healthy WT:
      start_state = "tg"   (PARK / PD model group)
      goal_state  = "WT"   (healthy control)
      alt_states  = ["mo", "PF"]
- `filter_data` narrows the analysis to a disease-relevant cell type (default
  Microglia; Astrocytes / Neuroblasts are left as commented alternatives).
- `genes_to_perturb_list` targets Parkinson's-linked genes. PARK7 (a
  protective factor; deleting it pushes embeddings AWAY from WT) is left in a
  commented block for a separate follow-up pass.

Everything else (device resolution, datasets<5 check, RAM estimate, nproc=1,
small max_ncells, combos=0, CLS embedding for V2) is inherited unchanged from 07.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ["WANDB_DISABLED"] = "true"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import resolve as _resolve_tissue

from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats
from geneformer.device import get_device


def warn_if_multiproc(name: str, nproc: int, max_ncells: int, n_genes: int) -> None:
    """Emit a visible warning when nproc>1 so users understand the constraints.

    nproc>1 is only safe for `perturb_data` when datasets<5 AND the materialised
    perturbation dataset stays small. Multi-gene / combos inflate variants, so
    nproc>1 with large max_ncells is an OOM/hang risk.
    """
    if nproc <= 1:
        return
    print(
        f"[WARN] {name}: nproc={nproc} (>1). "
        "This is only safe with datasets<5 (pin datasets==4.0.0; 5.x hangs "
        "InSilicoPerturber.perturb_data's dataset.map). "
        "Workers also add RAM that competes with the MPS forward pass. "
        f"With max_ncells={max_ncells} and {n_genes} perturbed gene(s), keep "
        "max_ncells modest (or reduce genes/combos) to avoid OOM. "
        "nproc=1 is the verified-stable default.",
        file=sys.stderr,
        flush=True,
    )


def _check_datasets_version() -> None:
    """Warn loudly if datasets>=5 (known perturb_data map hang)."""
    import datasets

    major = int(datasets.__version__.split(".")[0])
    if major >= 5:
        print(
            f"[WARN] datasets=={datasets.__version__} is installed. "
            "InSilicoPerturber.perturb_data hangs on dataset.map with datasets>=5. "
            "Pin datasets==4.0.0:  uv pip install datasets==4.0.0",
            file=sys.stderr,
            flush=True,
        )


def estimate_perturb_ram(n_cells: int, n_genes: int, combos: int, seq_len: int = 4096) -> float:
    """Rough RAM (GB) `perturb_data` will need to materialise the perturbation
    dataset. Each cell -> n_variants perturbation rows; combos=0 and N genes
    give n_variants ~ N (delete) per cell; combos>0 is the binomial C(N, combos+1)
    of combinations and explodes fast."""
    if combos > 0:
        import math

        n_variants = math.comb(n_genes, combos + 1) if n_genes >= combos + 1 else 1
    else:
        n_variants = max(n_genes, 1)
    total_variants = n_cells * n_variants
    # input_ids + perturb_index + length + dataset/arrow/pandas overhead, ~8-12 B/token
    return (total_variants * seq_len * 10) / 1e9


NPROC = int(os.environ.get("IS_NPROC", "1"))  # verified-stable default (see docstring)
_check_datasets_version()

ROOT, PREFIX = _resolve_tissue()
GF_ROOT = Path(os.environ["GENEFORMER_DIR"])
MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

# Fine-tuned cell classifier produced by 06_finetune.py. Point this at the
# `TRAINED_MODEL_PATH.txt` target if present, else fall back to the pretrained
# model (InSilicoPerturber can run in pretrained "inference" mode with
# model_type="Pretrained", but the tutorial uses a CellClassifier).
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
# Parkinson-specific: model disease state "tg" (PD model) -> goal "WT" control.
# state_key is the metadata column holding the state labels.
CELL_STATE_KEY = os.environ.get("ISP_STATE_KEY", "disease")
START_STATE = os.environ.get("ISP_START_STATE", "tg")
GOAL_STATE = os.environ.get("ISP_GOAL_STATE", "WT")
ALT_STATES = [s for s in os.environ.get("ISP_ALT_STATES", "mo,PF").split(",") if s]

# Target cell type; Kupffer.cells (liver-resident macrophage) is the PD_liver
# default (immune aSyn-clearance axis; closest liver analog of Microglia).
# Alternatives: "Hepatocytes" (metabolic), "Macrophages".
TARGET_CELLTYPE = os.environ.get("ISP_CELLTYPE", "Kupffer.cells")

cell_states_to_model = {
    "state_key": CELL_STATE_KEY,
    "start_state": START_STATE,
    "goal_state": GOAL_STATE,
    "alt_states": ALT_STATES,
}

# keep only cells of the disease-relevant cell type
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

# --------------------------------------------------------- step 2+3: perturb + stats
# NOTE: InSilicoPerturber expects Ensembl IDs, not gene symbols.
# PD_liver / Kupffer.cells gene set — every gene below is expressed in enough
# tg Kupffer cells to give a robust deletion (counts in ~1781 tg Kupffer cells:
# TYROBP=1462, APOE=1693, CTSD=1293, PARK7=467, PINK1=425, IL1B=532).
# The classic neuronal set (SNCA, LRRK2, GBA1, MAPT) is <30 cells here and is
# therefore omitted (those runs would be skipped as "no perturbable cells").
#
# PARK7/DJ-1 (ENSG00000116288) and PINK1 are protective factors (their deletion
# shifts AWAY from WT); a separate follow-up run can test that sign explicitly.
#
# Each gene is perturbed in its OWN run: `InSilicoPerturber` (combos=0, "all
# genes perturbed together") filters to cells containing EVERY listed gene
# (see pu.filter_data_by_tokens), so a batch empties out when even one gene is
# not expressed in the target cell type. Running one gene at a time isolates
# each deletion effect for goal_state_shift.
genes_to_perturb_list = [
    "ENSG00000011600",  # TYROBP
    "ENSG00000130203",  # APOE
    "ENSG00000117984",  # CTSD
    "ENSG00000116288",  # PARK7 (protective; delete -> away from WT)
    "ENSG00000158828",  # PINK1 (protective; delete -> away from WT)
    "ENSG00000125538",  # IL1B
    # "ENSG00000145335",  # SNCA   (only ~1 tg Kupffer cell -> empty)
    # "ENSG00000188906",  # LRRK2  (only ~19 cells -> empty)
    # "ENSG00000177628",  # GBA1   (0 cells -> empty)
    # "ENSG00000185345",  # PRKN   (79 cells; optional)
    # "ENSG00000186868",  # MAPT   (30 cells; optional)
]

# max_ncells for perturb_data. perturb_data materialises the full perturbation
# dataset in RAM (~n_cells x variants x seq_len); 2000 cells OOM'd on this MPS
# host. Keep this modest.
MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "200"))


def _gene_symbol(ensembl_id: str) -> str:
    """Return a short symbol for an Ensembl ID for output filenames."""
    name_to_id = {
        "ENSG00000011600": "TYROBP",
        "ENSG00000130203": "APOE",
        "ENSG00000117984": "CTSD",
        "ENSG00000145335": "SNCA",
        "ENSG00000188906": "LRRK2",
        "ENSG00000177628": "GBA1",
        "ENSG00000158828": "PINK1",
        "ENSG00000116288": "PARK7",
        "ENSG00000125538": "IL1B",
        "ENSG00000185345": "PRKN",
        "ENSG00000154277": "UCHL1",
        "ENSG00000186868": "MAPT",
    }
    return name_to_id.get(ensembl_id, ensembl_id.split(".")[0])


for gene in genes_to_perturb_list:
    symbol = _gene_symbol(gene)
    # Each gene gets its own sub-directory: `read_dictionaries` in
    # InSilicoPerturberStats globs every `*_raw.pickle` in the given directory,
    # so multiple genes must NOT share a directory or their stats get mixed.
    gene_isp_dir = ISP_DIR / f"{TARGET_CELLTYPE}" / f"gene_{symbol}"
    gene_isp_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"isp_gene_{symbol}"
    print(f"\n===== perturbing {symbol} ({gene}) -> {gene_isp_dir.name} =====", flush=True)

    _ram_est = estimate_perturb_ram(MAX_CELLS, 1, 0)
    print(
        f"perturb_data: max_ncells={MAX_CELLS}, n_genes=1, "
        f"estimated dataset RAM ~{_ram_est:.2f} GB",
        flush=True,
    )

    isp = InSilicoPerturber(
        perturb_type="delete",
        perturb_rank_shift=None,
        genes_to_perturb=[gene],
        combos=0,
        anchor_gene=None,
        model_type="CellClassifier",
        num_classes=len(celltypes),
        emb_mode="cls",  # V2 uses CLS token
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
            str(gene_isp_dir),
            prefix,
        )
    except Exception as exc:  # e.g. gene not expressed in the target cell type
        print(f"[SKIP] {symbol}: no perturbable cells ({type(exc).__name__})", flush=True)
        continue

    # step 3: stats — genes_perturbed must match the gene perturbed above.
    ispstats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed=[gene],
        combos=0,
        anchor_gene=None,
        cell_states_to_model=cell_states_to_model,
        model_version="V2",
    )
    ispstats.get_stats(
        str(gene_isp_dir),
        None,
        str(gene_isp_dir),
        prefix,
    )
    print(f"===== done {symbol} =====", flush=True)

print("=== IS PERTURBATION DONE ===", flush=True)