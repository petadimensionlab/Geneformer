#!/usr/bin/env python
"""In silico perturbation (deletion mode) for early AD detection — AD_brain.

Disease-state in silico perturbation, hypothesis-driven gene list, focused on
early brain-resident innate immunity (microglia / BAM).

Design:
  state_key    = "disease"
  start_state  = "AD"
  goal_state   = "WT"
  alt_states   = []
  filter_data  = brain-resident innate cell pool + pooled early timepoints
                 (AD3m/4.5m/6m vs WT3m/4.5m/6m)

  perturb_type = "delete", combos=0, emb_mode="cls", model_version="V2",
  one gene at a time (gene loop), state embeddings computed once.

  A gene whose deletion in early-AD cells shifts the embedding toward early WT
  (positive `Shift_to_goal_end`) is an early-driver / intervention-target
  candidate.

Gene list (human Ensembl IDs; AD_brain h5ad is already ortholog-mapped).
Genes not expressed in the pooled cells are auto-dropped with a warning.

All generic logic lives in `_isp_common.py`; this file only configures the
experiment (override the cell pool with ISP_CELLTYPES, comma-separated).
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
TISSUE = "brain"
EXPERIMENT = resolve_experiment(DISEASE, TISSUE)
ISP_DIR = isp_dir_for(ROOT, EXPERIMENT)
OUT_CSV = combined_csv_path(ROOT, EXPERIMENT)

EARLY_TP = ["AD3m", "AD4p5m", "AD6m", "WT3m", "WT4p5m", "WT6m"]
POOL_CELLTYPES = [  # override with ISP_CELLTYPES (comma-separated) if desired
    "Microglia",
]

# ---------------------------------------------------------------- gene list
# Hypothesis-driven fixed list (human Ensembl ENSG). All verified present in the
# token dictionary. DAM / neuroinflammation + AD GWAS risk + inflammatory cytokines.
AD_RISK = [
    "ENSG00000130203",  # APOE
    "ENSG00000095970",  # TREM2
    "ENSG00000011600",  # TYROBP
    "ENSG00000136717",  # BIN1
    "ENSG00000105383",  # CD33
    "ENSG00000064687",  # ABCA7
    "ENSG00000137642",  # SORL1
    "ENSG00000168918",  # INPP5D
    "ENSG00000197943",  # PLCG2
    "ENSG00000081189",  # MEF2C
]
DAM_MICROGLIA = [
    "ENSG00000095970",  # TREM2
    "ENSG00000011600",  # TYROBP
    "ENSG00000130203",  # APOE
    "ENSG00000182578",  # CSF1R
    "ENSG00000173372",  # C1QA
    "ENSG00000173369",  # C1QB
    "ENSG00000125730",  # C3
    "ENSG00000129226",  # CD68
    "ENSG00000204472",  # AIF1 (IBA1)
    "ENSG00000169896",  # ITGAM (CD11b)
    "ENSG00000140678",  # ITGAX (CD11c)
    "ENSG00000118785",  # SPP1
    "ENSG00000108691",  # CCL2
]
INFLAMMATION = [
    "ENSG00000125538",  # IL1B
    "ENSG00000232810",  # TNF
    "ENSG00000019582",  # CD74
]

genes_to_perturb = list(dict.fromkeys(AD_RISK + DAM_MICROGLIA + INFLAMMATION))
print(f"candidate genes (hypothesis list): {len(genes_to_perturb)}", flush=True)

# ---------------------------------------------------------------- gene presence check
tok2gene, ensg2name, name2ensg = load_gene_dicts()
_token_vocab = set(tok2gene.values())

from datasets import load_from_disk  # noqa: E402

tokens = load_from_disk(str(TOKENIZED))

present = scan_pool_presence(tokens, tok2gene, POOL_CELLTYPES, samples4_filter=EARLY_TP)
keep = [g for g in genes_to_perturb if g in _token_vocab and g in present]
dropped = [g for g in genes_to_perturb if g not in keep]
if dropped:
    print(
        f"[WARN] dropping {len(genes_to_perturb) - len(keep)} genes not expressed "
        f"in pool cells: {dropped}",
        flush=True,
    )
genes_to_perturb = keep
print(f"genes_to_perturb (expressed in pool): {len(genes_to_perturb)}", flush=True)

# num_classes: total celltype classes in the trained classifier
celltypes = sorted(set(tokens["celltype"]))
NUM_CLASSES = len(celltypes)

# ---------------------------------------------------------------- states
cell_states_to_model = {
    "state_key": "disease",
    "start_state": "AD",
    "goal_state": "WT",
    "alt_states": [],
}
filter_data_dict = {
    "celltype": POOL_CELLTYPES,
    "samples4": EARLY_TP,
}

# ---------------------------------------------------------------- step 1: state embs
embex = EmbExtractor(
    model_type="CellClassifier",
    num_classes=NUM_CLASSES,
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

# ---------------------------------------------------------------- step 2+3: per-gene perturb + stats
# SINGLE-gene group perturbation per candidate gene: `genes_to_perturb=[gene]`
# produces a per-gene goal_state_shift whose Shift_to_goal_end is exactly the
# early-detection metric (deleting this gene in early-AD cells -> how far toward
# early-WT). Each single-gene run is fast (one perturbed variant per cell).
MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "200"))
warn_if_multiproc("InSilicoPerturber", NPROC, MAX_CELLS, 1)
print(f"per-gene screening: {len(genes_to_perturb)} genes x max_ncells={MAX_CELLS}", flush=True)

results_rows = []
for gene in genes_to_perturb:
    gene_name = ensg2name.get(gene, gene)
    gene_dir = ISP_DIR / gene_name
    gene_dir.mkdir(parents=True, exist_ok=True)
    _ram_est = estimate_perturb_ram(MAX_CELLS, 1, 0)
    print(
        f"\n=== gene {gene_name} ({gene}): max_ncells={MAX_CELLS}, "
        f"est. RAM ~{_ram_est:.2f} GB ===",
        flush=True,
    )
    shift = perturb_one_gene(
        model_dir=CELLCLASSIFIER_DIR,
        tokenized=str(TOKENIZED),
        out_dir=gene_dir,
        gene_ensg=gene,
        cell_states_to_model=cell_states_to_model,
        filter_data=filter_data_dict,
        state_embs_dict=state_embs_dict,
        model_type="CellClassifier",
        num_classes=NUM_CLASSES,
        max_ncells=MAX_CELLS,
        nproc=NPROC,
    )
    if shift is None:
        print(f"[WARN] {gene_name}: no Shift_to_goal_end in stats", flush=True)
    else:
        print(f"{gene_name}: Shift_to_goal_end = {shift:.4f}", flush=True)
    results_rows.append({
        "Timepoint": "early",
        "Gene": gene_name,
        "Ensembl_ID": gene,
        "Shift_to_goal_end": shift,
    })

# ---------------------------------------------------------------- step 4: aggregate candidates
print(
    "\n=== CANDIDATE EARLY-DETECTION RANKING (positive Shift_to_goal_end = "
    "deleting this gene in early-AD pushes toward early-WT) ===",
    flush=True,
)
write_combined(results_rows, OUT_CSV, rank_tp=None)

print("=== IS PERTURBATION (EARLY AD BRAIN) DONE ===", flush=True)