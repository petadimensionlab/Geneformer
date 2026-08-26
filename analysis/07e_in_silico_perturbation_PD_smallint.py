#!/usr/bin/env python
"""In silico perturbation (deletion mode) on MPS — Parkinson's-disease-specific
configuration for PD_smallint (mouse small intestine), built on the same engine
as `analysis/07_in_silico_perturbation.py`.

Key adaptations vs. 07:
- `07` picks the two most-abundant `celltype`s as start/goal (disease-agnostic).
  Here `state_key="disease"` so the start state is the PD-model group and the
  goal is healthy WT:
      start_state = "tg"   (PARK / PD model group)
      goal_state  = "WT"   (healthy control)
      alt_states  = ["mo", "PF"]
- `filter_data` narrows the analysis to the gut-macrophage / enteric-neuroimmune
  axis. Default TARGET_CELLTYPE is Macrophages (the small-intestine resident
  macrophages, the gut analog of PD_brain Microglia and PD_liver Kupffer.cells).
  Rationale: neuronal PD genes (SNCA, LRRK2, GBA1) are essentially unexpressed
  in small-intestine cells (<=~2% of cells), so they cannot be deleted
  meaningfully; Macrophages robustly express the PD-modified immune /
  autophagy / lysosomal gene set below.
- `genes_to_perturb_list` targets Parkinson's-linked genes with adequate
  expression in tg (PD) Macrophages. APOE/TYROBP/CTSD/ATG7/IL1B are
  disease-driver candidates (deleting them should push embeddings TOWARD WT).
  PINK1/PARK7 are protective factors (deleting them should push AWAY from WT);
  they serve as sign/reverse checks. This set matches PD_liver's verified list,
  enabling cross-tissue comparison.

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

# Target cell type; Macrophages (gut-resident) is the PD_smallint default —
# the small-intestine analog of PD_brain Microglia / PD_liver Kupffer.cells
# (immune aSyn-clearance and enteric-neuroinflammation axis).
# Alternatives: "Epithelial.cells" (n=1800; PARK7/PINK1/APOE/CTSD/ATG7
# expressed), "Goblet.cells_M.cells" (n=594; mucin barrier). SNCA/LRRK2/GBA1
# are unexpressed in small intestine (see header), so deleting them is skipped.
TARGET_CELLTYPE = os.environ.get("ISP_CELLTYPE", "Macrophages")

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
# PD_smallint / Macrophages gene set — every gene below is expressed in enough
# tg Macrophages to give a robust deletion (n_tg=514; % of cells expressing:
# APOE=95.5, TYROBP=70.8, CTSD=17.1, ATG7=17.1, IL1B=15.2, PINK1=8.2,
# PARK7=6.0). The classic neuronal set (SNCA, LRRK2, GBA1) is ~unexpressed in
# small-intestine cells (SNCA<=2.5%, LRRK2<=1.2%) and is therefore omitted
# (those runs would be skipped as "no perturbable cells").
#
# APOE/TYROBP/CTSD/ATG7/IL1B are disease-driver candidates (deletion is
# expected to shift embeddings TOWARD WT). PINK1/PARK7 are protective factors
# (deletion expected to shift AWAY from WT) and act as sign checks. This set
# matches the PD_liver / Kupffer.cells list, enabling cross-tissue comparison.
#
# Each gene is perturbed in its OWN run: `InSilicoPerturber` (combos=0, "all
# genes perturbed together") filters to cells containing EVERY listed gene
# (see pu.filter_data_by_tokens), so a batch empties out when even one gene is
# not expressed in the target cell type. Running one gene at a time isolates
# each deletion effect for goal_state_shift.
genes_to_perturb_list = [
    "ENSG00000130203",  # APOE   (95.5% tg Macrophages)
    "ENSG00000011600",  # TYROBP (70.8%)
    "ENSG00000117984",  # CTSD   (17.1%)
    "ENSG00000197548",  # ATG7   (17.1%)
    "ENSG00000125538",  # IL1B   (15.2%)
    "ENSG00000158828",  # PINK1  (8.2%; protective; delete -> away from WT)
    "ENSG00000116288",  # PARK7  (6.0%; protective; delete -> away from WT)
    # "ENSG00000145335",  # SNCA   (~2.5% max -> empty)
    # "ENSG00000188906",  # LRRK2  (~1.2% max -> empty)
    # "ENSG00000177628",  # GBA1   (0% -> empty)
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
        "ENSG00000197548": "ATG7",
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