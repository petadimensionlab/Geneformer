#!/usr/bin/env python
"""In silico perturbation (deletion mode) for early AD detection — AD_LN
(lymph node). Small-scale CD4+ T-cell test with AD→WT goal-state shift.

Original (07_in_silico_perturbation_AD_LN.py) was a CPU smoke/test run. This
refactored version follows the unified per-organ scheme:

  state_key    = "disease", start = "AD", goal = "WT", alt = []
  celltype     = CD4.T.cells only (small test) — override with ISP_CELLTYPES
  genes        = CD28 / STAT3 / FOXP3 (Treg / Th17 costimulation axis)
  max_ncells   = IS_MAX_CELLS (default 50 for the small test)

Uses the shared `_isp_common.py` helpers and the canonical output layout
`input/<TISSUE>/results/isp/ad_ln/`. State embeddings computed once and
reused across genes; each gene is deleted individually (combos=0) so
`InSilicoPerturberStats` sees only that gene.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ["WANDB_DISABLED"] = "true"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _isp_common import (
    check_datasets_version,
    combined_csv_path,
    estimate_perturb_ram,
    isp_dir_for,
    load_gene_dicts,
    perturb_one_gene,
    resolve_classifier_dir,
    resolve_experiment,
    scan_pool_presence,
    warn_if_multiproc,
    write_combined,
)
from _resolve_tissue import resolve as _resolve_tissue

from geneformer import EmbExtractor
from geneformer.device import get_device

NPROC = int(os.environ.get("IS_NPROC", "1"))
check_datasets_version()
print("device:", get_device(), flush=True)

ROOT, PREFIX = _resolve_tissue()
GF_ROOT = Path(os.environ["GENEFORMER_DIR"])
MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

CELLCLASSIFIER_DIR = resolve_classifier_dir(ROOT, GF_ROOT, MODEL_NAME)

# ---------------------------------------------------------------- config
DISEASE = "AD"
TISSUE = "ln"
EXPERIMENT = resolve_experiment(DISEASE, TISSUE)
ISP_DIR = isp_dir_for(ROOT, EXPERIMENT)
OUT_CSV = combined_csv_path(ROOT, EXPERIMENT)

# Small-scale test: CD4+ T cells only (override with ISP_CELL_POOL).
POOL_CELLTYPES = [
    s for s in os.environ.get("ISP_CELL_POOL", "CD4.T.cells").split(",") if s
]

# Candidate genes (Ensembl IDs verified in the token dictionary). The original
# test used CD28 / STAT3 / FOXP3 (Treg/Thim costimulation & differentiation).
HYPOTHESIS_GENES = [
    ("ENSG00000178562", "CD28"),   # costimulation
    ("ENSG00000168610", "STAT3"),  # Th17/Treg differentiation
    ("ENSG00000049768", "FOXP3"),  # Treg master TF
]
# presence filter + dedupe on the fixed list (no per-symbol ENSG mapping needed)
genes_to_perturb = [(g, s) for g, s in HYPOTHESIS_GENES]

MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "50"))
EMB_CELLS = int(os.environ.get("IS_EMB_CELLS", "1000"))
warn_if_multiproc("InSilicoPerturber", NPROC, MAX_CELLS, 1)

# ---------------------------------------------------------------- data presence
from datasets import load_from_disk  # noqa: E402

tokens = load_from_disk(str(TOKENIZED))
celltypes = sorted(set(tokens["celltype"]))
NUM_CLASSES = len(celltypes)
print(f"{NUM_CLASSES} cell types in dataset", flush=True)

tok2gene, ensg2name, name2ensg = load_gene_dicts()
present = scan_pool_presence(tokens, tok2gene, POOL_CELLTYPES)
genes_to_perturb = [
    (g, s) for g, s in genes_to_perturb
    if g in present
]
print(f"genes_to_perturb (expressed in pool): {len(genes_to_perturb)}", flush=True)
if not genes_to_perturb:
    sys.exit("No candidate genes expressed in the pool cells.")

# ---------------------------------------------------------------- states
cell_states_to_model = {
    "state_key": "disease",
    "start_state": "AD",
    "goal_state": "WT",
    "alt_states": [],
}
filter_data = {"celltype": POOL_CELLTYPES}

# ---------------------------------------------------------------- step 1: state embs (shared)
embex = EmbExtractor(
    model_type="CellClassifier",
    num_classes=NUM_CLASSES,
    filter_data=filter_data,
    max_ncells=EMB_CELLS,
    emb_layer=0,
    summary_stat="exact_mean",
    forward_batch_size=4,   # small test / CPU-friendly
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

# ---------------------------------------------------------------- per-gene perturb + stats
results = []
for gene, symbol in genes_to_perturb:
    gene_dir = ISP_DIR / symbol
    gene_dir.mkdir(parents=True, exist_ok=True)
    _ram_est = estimate_perturb_ram(MAX_CELLS, 1, 0)
    print(
        f"\n=== gene {symbol} ({gene}): max_ncells={MAX_CELLS}, "
        f"est. RAM ~{_ram_est:.2f} GB ===",
        flush=True,
    )
    shift = perturb_one_gene(
        model_dir=CELLCLASSIFIER_DIR,
        tokenized=str(TOKENIZED),
        out_dir=gene_dir,
        gene_ensg=gene,
        cell_states_to_model=cell_states_to_model,
        filter_data=filter_data,
        state_embs_dict=state_embs_dict,
        model_type="CellClassifier",
        num_classes=NUM_CLASSES,
        max_ncells=MAX_CELLS,
        nproc=NPROC,
    )
    if shift is None:
        print(f"[WARN] {symbol}: no Shift_to_goal_end in stats", flush=True)
    results.append({
        "Timepoint": "test",
        "Gene": symbol,
        "Ensembl_ID": gene,
        "Shift_to_goal_end": shift,
    })

# ---------------------------------------------------------------- combine
write_combined(results, OUT_CSV, rank_tp=None)

print("=== IS PERTURBATION (AD LN) DONE ===", flush=True)