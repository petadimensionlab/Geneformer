#!/usr/bin/env python
"""In silico perturbation (deletion mode) for early AD detection — AD_smallint.

Disease-state in silico perturbation with a hypothesis-driven gene list.

Design:
  state_key    = "disease"
  start_state  = "AD"
  goal_state   = "WT"
  alt_states   = []

  filter_data  = immune + gut-epithelial cell types, pooled early timepoints
                 (3-6m) — a single screen over the early pool.

  perturb_type = "delete", combos=0, emb_mode="cls", model_version="V2". One
  gene at a time (gene loop); state embeddings computed once and shared.

  A gene whose deletion yields a large positive `Shift_to_goal_end` pushes
  early-AD cells toward early-WT — a candidate early driver / target.

Gene list is a fixed hypothesis list (ENSG IDs verified in the token
dictionary): AD risk (GWAS) + inflammation/immune + gut-barrier (mucins,
tight junctions, antimicrobial peptides, stem/metaplasia). Genes absent from
the pool cells are skipped with a warning.

All generic logic lives in `_isp_common.py`; this file only configures the
experiment.
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
TISSUE = "smallint"
EXPERIMENT = resolve_experiment(DISEASE, TISSUE)
ISP_DIR = isp_dir_for(ROOT, EXPERIMENT)
OUT_CSV = combined_csv_path(ROOT, EXPERIMENT)

EARLY_TP = ["AD3m", "AD4p5m", "AD6m", "WT3m", "WT4p5m", "WT6m"]
POOL_CELLTYPES = [
    "Macrophages", "CD4.T.cells", "CD8.T.cells", "NK_ILC1",
    "IECs", "Paneth.cells", "Goblet.cells_M.cells", "Epithelial.cells",
]

# ---------------------------------------------------------------- gene list
AD_RISK = [
    "ENSG00000130203",  # APOE
    "ENSG00000095970",  # TREM2
    "ENSG00000120885",  # CLU
    "ENSG00000136717",  # BIN1
    "ENSG00000105383",  # CD33
    "ENSG00000064687",  # ABCA7
    "ENSG00000137642",  # SORL1
    "ENSG00000011600",  # TYROBP
    "ENSG00000168918",  # INPP5D
    "ENSG00000197943",  # PLCG2
    "ENSG00000081189",  # MEF2C
]
INFLAMMATION = [
    "ENSG00000173372",  # C1QA
    "ENSG00000173369",  # C1QB
    "ENSG00000129226",  # CD68
    "ENSG00000182578",  # CSF1R
    "ENSG00000204472",  # AIF1
    "ENSG00000168329",  # CX3CR1
    "ENSG00000125538",  # IL1B
    "ENSG00000232810",  # TNF
    "ENSG00000136869",  # TLR4
    "ENSG00000090382",  # LYZ
    "ENSG00000130208",  # APOC1
    "ENSG00000204287",  # HLA-DRA
]
GUT_BARRIER_MUCUS = [
    "ENSG00000198788",  # MUC2
    "ENSG00000185499",  # MUC1
    "ENSG00000169894",  # MUC3A
    "ENSG00000173702",  # MUC13
    "ENSG00000169876",  # MUC17
    "ENSG00000205277",  # MUC12
    "ENSG00000117983",  # MUC5B
    "ENSG00000145113",  # MUC4
]
GUT_BARRIER_TJ = [
    "ENSG00000197822",  # OCLN
    "ENSG00000163347",  # CLDN1
    "ENSG00000165376",  # CLDN2
    "ENSG00000165215",  # CLDN3
    "ENSG00000189143",  # CLDN4
    "ENSG00000181885",  # CLDN7
    "ENSG00000106404",  # CLDN15
    "ENSG00000104067",  # TJP1 (ZO-1)
    "ENSG00000119139",  # TJP2
    "ENSG00000105289",  # TJP3 (ZO-3)
    "ENSG00000158769",  # F11R (JAM-1)
    "ENSG00000154721",  # JAM2
]
GUT_BARRIER_SECRETED = [
    "ENSG00000090382",  # LYZ
    "ENSG00000164825",  # DEFB1
    "ENSG00000176797",  # DEFB103A
    "ENSG00000143954",  # REG3G
    "ENSG00000164822",  # DEFA6
    "ENSG00000164047",  # CAMP
]
GUT_BARRIER_STEM_META = [
    "ENSG00000139292",  # LGR5
    "ENSG00000102837",  # OLFM4
    "ENSG00000118785",  # SPP1
    "ENSG00000106541",  # AGR2
    "ENSG00000160180",  # TFF3
    "ENSG00000016490",  # CLCA1
    "ENSG00000143546",  # S100A8
    "ENSG00000163220",  # S100A9
]

genes_to_perturb = (
    AD_RISK + INFLAMMATION
    + GUT_BARRIER_MUCUS + GUT_BARRIER_TJ
    + GUT_BARRIER_SECRETED + GUT_BARRIER_STEM_META
)
genes_to_perturb = list(dict.fromkeys(genes_to_perturb))  # dedupe (LYZ)
print(f"candidate genes (hypothesis list): {len(genes_to_perturb)}", flush=True)

# ---------------------------------------------------------------- presence check
tok2gene, ensg2name, name2ensg = load_gene_dicts()
_token_vocab = set(tok2gene.values())

from datasets import load_from_disk  # noqa: E402

tokens = load_from_disk(str(TOKENIZED))

present = scan_pool_presence(tokens, tok2gene, POOL_CELLTYPES, samples4_filter=EARLY_TP)
genes_to_perturb = [g for g in genes_to_perturb if g in _token_vocab and g in present]
print(f"genes_to_perturb (expressed in early pool): {len(genes_to_perturb)}", flush=True)
if not genes_to_perturb:
    sys.exit("No candidate genes expressed in the pool cells.")

# optional subset for smoke tests: IS_MAX_GENES limits how many genes run
_max_genes = os.environ.get("IS_MAX_GENES")
if _max_genes:
    genes_to_perturb = genes_to_perturb[: int(_max_genes)]
    print(f"[smoke] limited to {len(genes_to_perturb)} genes", flush=True)

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

# ---------------------------------------------------------------- step 1: state embs (shared)
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

# ---------------------------------------------------------------- per-gene perturbation
MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "300"))
warn_if_multiproc("InSilicoPerturber", NPROC, MAX_CELLS, 1)

results = []
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
    results.append({
        "Timepoint": "early",
        "Gene": gene_name,
        "Ensembl_ID": gene,
        "Shift_to_goal_end": shift,
    })

# ---------------------------------------------------------------- combine results
write_combined(results, OUT_CSV, rank_tp=None)

print("=== IS PERTURBATION (EARLY AD SMALLINT) DONE ===", flush=True)