#!/usr/bin/env python
"""In silico perturbation (deletion mode) on MPS — Alzheimer's-disease-specific
configuration for AD_BM (bone marrow), built on the same engine as
`analysis/07b_in_silico_perturbation_PD.py` (which is itself a disease-state
adaptation of `analysis/07_in_silico_perturbation.py`).

Scientific goal (AD early-detection axis in bone marrow):
  Model the AD (disease) myeloid state as start_state, WT as goal_state, then
  simulate DELETION of AD- / inflammation-linked genes and ask whether the AD
  cell-state embedding shifts TOWARD the healthy (WT) state
  (`mode="goal_state_shift"`). Genes whose deletion pulls AD cells toward WT
  are candidate early-detection/risk factors for AD.

Key adaptations vs. 07b (PD_brain):
- Cell lineage: two bone-marrow myeloid cell types are each modelled in their
  own run (myeloid axis is the leading peripheral-immune axis implicated in AD):
      celltype A: "Ly6c.high.classical.Monocytes"  (CCR2+ classical monocytes)
      celltype B: "Macrophage"                      (bone-marrow macrophages)
- `state_key="disease"`, start=AD, goal=WT, no alt states.
- `genes_to_perturb_list` are AD/inflammation candidates, each run in its own
  gene->celltype sub-directory so InSilicoPerturberStats never mixes
  (single-gene combos=0 "all perturbed together" semantics).

Compute notes inherited from 07/07b (all apply on this MPS host):
  datasets must be < 5 (5.x hangs perturb_data's dataset.map) -> pin 4.0.0.
  nproc=1 default, max_ncells=200 default (2000 OOM'd), combos=0, CLS emb.

Targets are pre-verified against the AD_BM counts matrix so every celltype
x gene pair has >=40% of AD cells expressing the gene (>= the ~few-cell floor
that fills `filter_data_by_tokens`); a low-expressed gene would otherwise
silence a whole combined batch.
"""
from __future__ import annotations

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
    # input_ids + perturb_index + length + dataset overhead, ~8-12 B/token
    return (total_variants * seq_len * 10) / 1e9


NPROC = int(os.environ.get("IS_NPROC", "1"))  # verified-stable default
_check_datasets_version()

ROOT, PREFIX = _resolve_tissue()
GF_ROOT = Path(os.environ["GENEFORMER_DIR"])
MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

# Fine-tuned cell classifier produced by 06_finetune.py. Point this at the
# `TRAINED_MODEL_PATH.txt` target if present, else fall back to the pretrained
# model (InSilicoPerturber can run in pretrained "inference" mode).
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
# Alzheimer start -> WT goal. state_key is the metadata column holding the
# state labels; AD_BM tokenized data carries both `disease` (AD/WT) and
# `celltype`.
CELL_STATE_KEY = os.environ.get("ISP_STATE_KEY", "disease")
START_STATE = os.environ.get("ISP_START_STATE", "AD")
GOAL_STATE = os.environ.get("ISP_GOAL_STATE", "WT")
ALT_STATES = [s for s in os.environ.get("ISP_ALT_STATES", "").split(",") if s]

# Bone-marrow myeloid cell types to model, one run each. Every gene below was
# pre-validated to be expressed in >=40% of the AD cells of BOTH cell types
# (see module docstring), so no combined batch is expected to come up empty.
TARGET_CELLTYPES = [
    s
    for s in os.environ.get("ISP_CELLTYPES", "Ly6c.high.classical.Monocytes,Macrophage").split(",")
    if s
]

cell_states_to_model = {
    "state_key": CELL_STATE_KEY,
    "start_state": START_STATE,
    "goal_state": GOAL_STATE,
    "alt_states": ALT_STATES,
}

from datasets import load_from_disk

tokens = load_from_disk(str(TOKENIZED))
celltypes = sorted(set(tokens["celltype"]))
print(f"{len(celltypes)} cell types available", flush=True)
print(f"modeled states: key={CELL_STATE_KEY} start={START_STATE} goal={GOAL_STATE} "
      f"alt={ALT_STATES}", flush=True)
print(f"target celltypes: {TARGET_CELLTYPES}", flush=True)

# ---------------------------------------------------------------- step 1: state embs
# State embeddings are per-disease-state (AD vs WT), not per-celltype, so a
# single extraction pass over ALL modeled cells (both cell types) suffices and
# is reused by every per-celltype perturb run.
filter_data_all_states = {
    "disease": [START_STATE, GOAL_STATE] + ALT_STATES,
}

embex = EmbExtractor(
    model_type="CellClassifier",
    num_classes=len(celltypes),
    filter_data=filter_data_all_states,
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

# ---------------------------------------------------------------- step 2+3 per celltype/gene

# NOTE: InSilicoPerturber expects Ensembl IDs, not gene symbols. Harbour below.
# AD / inflammation candidates, each verified (counts matrix) to be expressed in
# >=40% of AD cells of BOTH target cell types:
#   Ly6c.high classical Monocytes / Macrophage:
#   TYROBP 100/96, APOE 58/94, S100A9 99/97, S100A8 92/92, LGALS3 98/82,
#   CCR2 96/50, CSF1R 89/81, IRF8 82/67, SPI1 82/56, CTSD 46/60, CTSB 95/93
# Each gene is perturbed in its OWN run: `InSilicoPerturber` (combos=0, "all
# genes perturbed together") filters to cells containing EVERY listed gene
# (pu.filter_data_by_tokens), so a batch empties out when even one gene is not
# expressed in the target cell type. Running one gene at a time isolates each
# deletion effect for goal_state_shift.
GENES_TO_PERTURB = [
    ("ENSG00000011600", "TYROBP"),  # DAP12 adapter; myeloid immune axis
    ("ENSG00000130203", "APOE"),    # largest AD risk factor; lipid-调制 immune
    ("ENSG00000163220", "S100A9"),  # calprotectin bait; blood AD biomarker
    ("ENSG00000143546", "S100A8"),  # calprotectin beta dimer; AD biomarker
    ("ENSG00000131981", "LGALS3"),  # Galectin-3; disease-associated microglia
    ("ENSG00000121807", "CCR2"),    # monocyte migration / brain infiltration
    ("ENSG00000182578", "CSF1R"),   # macrophage survival / microglia
    ("ENSG00000140968", "IRF8"),    # monocyte/TLN master TF
    ("ENSG00000066336", "SPI1"),    # PU.1; AD microglial TF
    ("ENSG00000117984", "CTSD"),    # Cathepsin D; lysosome / Ab processing
    ("ENSG00000164733", "CTSB"),    # Cathepsin B; lysosome / Ab processing
]

# max_ncells for perturb_data. perturb_data materialises the full perturbation
# dataset in RAM (~n_cells x variants x seq_len); 2000 cells OOM'd on this MPS
# host. Keep this modest.
MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "200"))


for celltype in TARGET_CELLTYPES:
    filter_data_dict = {
        "celltype": [celltype],
        "disease": [START_STATE, GOAL_STATE] + ALT_STATES,
    }
    print(f"\n===== processing celltype: {celltype} =====", flush=True)

    for gene, symbol in GENES_TO_PERTURB:
        # Each gene gets its own sub-directory: `read_dictionaries` in
        # InSilicoPerturberStats globs every `*_raw.pickle` in the given dir,
        # so multiple genes must NOT share a directory or their stats mix.
        gene_isp_dir = ISP_DIR / celltype.replace(" ", "_") / f"gene_{symbol}"
        gene_isp_dir.mkdir(parents=True, exist_ok=True)
        prefix = f"isp_gene_{symbol}"
        print(f"\n  ==== perturbing {symbol} ({gene}) in {celltype} ====", flush=True)

        _ram_est = estimate_perturb_ram(MAX_CELLS, 1, 0)
        print(
            f"  perturb_data: max_ncells={MAX_CELLS}, n_genes=1, "
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
            print(f"  [SKIP] {symbol}: no perturbable cells ({type(exc).__name__}): {exc}",
                  flush=True)
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
        print(f"  ==== done {symbol} ====", flush=True)

print("=== IS PERTURBATION DONE ===", flush=True)