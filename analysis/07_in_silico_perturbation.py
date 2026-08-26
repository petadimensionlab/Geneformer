#!/usr/bin/env python
"""In silico perturbation (deletion mode) on MPS — ported from the Geneformer
`examples/in_silico_perturbation.ipynb` tutorial and adapted for V2-104M.

Key adaptations vs. the notebook:
- V1  -> V2 (V2 uses <cls> token; emb_mode="cls"; cls_present=True).
- The notebook's heart-disease states (dcm/nf/hcm) are replaced by cell-type
  states drawn from the tokenized PD_smallint dataset.
- Device is auto-detected (MPS on this Mac) via geneformer.device.

This mirrors the notebook flow: (1) extract state embeddings, (2) perturb_data
to emit intermediate files, (3) InSilicoPerturberStats to get the final CSV.

--------------------------------------------------------------------------
Compute best practices (MPS host) — learned from repeated hangs / OOMs:
- **datasets must be < 5 (4.x).** `datasets>=5` breaks `dataset.map()` and the
  `InSilicoPerturber.perturb_data` step hangs on the map with *any* nproc.
  Pin `datasets==4.0.0` in the env. This is the true root-cause fix; nothing
  about nproc can compensate for the 5.x map hang.
- **nproc**: `make_group_perturbation_batch` is a nested closure. With
  datasets 4.x it pickles fine under `spawn` (verified), so nproc>1 *can*
  work for `perturb_data`'s CPU-only map. The safe default below is nproc=1
  because (a) workers add memory that competes with the MPS forward pass and
  (b) with 2+ genes / combos>0 the per-cell variant count grows and the map
  output becomes type-inconsistent under multiprocessing. nproc=1 is the
  verified-stable configuration; a warning is printed if you override it.
- **max_ncells / OOM**: `perturb_data` materialises the full perturbation
  dataset in RAM (~ n_cells x n_variants_per_cell x seq_len). Single-gene
  delete keeps n_variants ~ n_cells, but N genes (combos=0) multiply by N and
  `combos>0` is combinatorial. 2000 cells OOM'd on this host. Cap max_ncells
  to keep `n_cells * n_variants * seq_len * ~8 B` well under free RAM.
"""
from __future__ import annotations

import json
import os
import pickle
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
    give n_variants ~ N (delete) / N (overexpress) per cell; combos>0 is the
    binomial C(N, combos+1) of combinations and explodes fast."""
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
# Pick cell-type states that actually exist in the tokenized data.
from datasets import load_from_disk

tokens = load_from_disk(str(TOKENIZED))
celltypes = sorted(set(tokens["celltype"]))
print(f"{len(celltypes)} cell types available", flush=True)
# choose the two most abundant as start/goal, next as alt
counts = {}
for ct in celltypes:
    counts[ct] = sum(1 for x in tokens["celltype"] if x == ct)
order = sorted(counts, key=lambda k: counts[k], reverse=True)
start_state, goal_state = order[0], order[1]
alt_states = order[2:5]
print(f"start={start_state} goal={goal_state} alt={alt_states}", flush=True)

cell_states_to_model = {
    "state_key": "celltype",
    "start_state": start_state,
    "goal_state": goal_state,
    "alt_states": alt_states,
}

# keep only cells belonging to the modeled states
modeled_states = set([start_state, goal_state] + alt_states)
filter_data_dict = {"celltype": list(modeled_states)}

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

# ---------------------------------------------------------------- step 2: perturb

# NOTE: InSilicoPerturber expects Ensembl IDs, not gene symbols, and the gene
# MUST be present in cells of the start_state (or at all), otherwise
# perturb_data filters every cell out and raises "No cells in dataset contain
# genes_to_perturb". A hardcoded gene like SNCA (ENSG00000145335, a brain gene)
# is absent from the modeled immune/skin cell types. So auto-select genes
# present in start_state cells; override the count with IS_N_GENES.
from geneformer import TOKEN_DICTIONARY_FILE

N_GENES = int(os.environ.get("IS_N_GENES", "2"))
with open(TOKEN_DICTIONARY_FILE, "rb") as f:
    token_dictionary = pickle.load(f)
tok2gene = {v: k for k, v in token_dictionary.items() if k.startswith("ENSG")}
gene_cell_count = {}
for ex in tokens:
    if ex["celltype"] != start_state:
        continue
    for t in ex["input_ids"]:
        gene = tok2gene.get(t)
        if gene is not None:
            gene_cell_count[gene] = gene_cell_count.get(gene, 0) + 1
gene_order = sorted(
    gene_cell_count, key=lambda g: gene_cell_count[g], reverse=True
)
genes_to_perturb_list = gene_order[:N_GENES]
if genes_to_perturb_list:
    genes_to_perturb = genes_to_perturb_list.copy()
else:
    genes_to_perturb = "all"
print(
    f"genes_to_perturb (present in {start_state} cells): {genes_to_perturb}",
    flush=True,
)

# max_ncells for perturb_data. perturb_data materialises the full perturbation
# dataset in RAM (~n_cells x variants x seq_len); 2000 cells OOM'd on this MPS
# host. Because ~28-32 GB swap was used by a 2000-cell / nproc>1 run, keep this
# modest. Estimate before running and lower if free RAM is tight.
MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "200"))
_ram_est = estimate_perturb_ram(
    MAX_CELLS, len(genes_to_perturb_list) if genes_to_perturb != "all" else 20, 0
)
print(
    f"perturb_data: max_ncells={MAX_CELLS}, n_genes={len(genes_to_perturb_list) if genes_to_perturb != 'all' else 'all'}, "
    f"estimated dataset RAM ~{_ram_est:.2f} GB",
    flush=True,
)
warn_if_multiproc("InSilicoPerturber", NPROC, MAX_CELLS,
                  len(genes_to_perturb_list) if genes_to_perturb != "all" else 20)

isp = InSilicoPerturber(
    perturb_type="delete",
    perturb_rank_shift=None,
    genes_to_perturb=genes_to_perturb,
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

isp.perturb_data(
    CELLCLASSIFIER_DIR,
    str(TOKENIZED),
    str(ISP_DIR),
    "isp_perturbation",
)

# ---------------------------------------------------------------- step 3: stats
# genes_perturbed must match the genes actually perturbed in step 2.
# Using "all" here re-runs stats over every gene in the vocabulary (extremely
# slow on MPS), so restrict it to the same gene(s) as genes_to_perturb.
ispstats = InSilicoPerturberStats(
    mode="goal_state_shift",
    genes_perturbed=genes_to_perturb,
    combos=0,
    anchor_gene=None,
    cell_states_to_model=cell_states_to_model,
    model_version="V2",
)

ispstats.get_stats(
    str(ISP_DIR),
    None,
    str(ISP_DIR),
    "isp_stats",
)
print("=== IS PERTURBATION DONE ===", flush=True)
