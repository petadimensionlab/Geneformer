#!/usr/bin/env python
"""In silico perturbation (deletion mode) for early AD detection — AD_brain.

Disease-state in silico perturbation, hypothesis-driven gene list, focused on
early brain-resident innate immunity (microglia + BAM).

Design (mirrors 07_ad_spleen_early_isp but retargeted from gut/immune to brain):
  state_key    = "disease"
  start_state  = "AD"
  goal_state   = "WT"
  alt_states   = []
  filter_data  = Microglia + BAM, early timepoints (3m/4.5m/6m)

  perturb_type = "delete", combos=0, emb_mode="cls" (V2), model_version="V2"

A gene whose deletion in early-AD cells shifts the embedding toward early WT
(positive `Shift_to_goal_end`) is an early-driver / intervention-target
candidate: i.e. knocking it out in AD microglia/BAM pushes them toward the
healthy state -> candidate early-detection marker.

Gene list (human Ensembl IDs; the AD_brain h5ad is already ortholog-mapped to
ENSG). Genes not expressed in the pooled cells are auto-dropped with a warning.

Compute constraints (inherited from 07/07_ad_spleen_early_isp, see their docstrings):
- datasets must be < 5 (4.0.0) or perturb_data's dataset.map hangs.
- nproc=1 (verified-stable); max_ncells kept modest (RAM estimate printed).
- stats `genes_perturbed` restricted to the same list as `genes_to_perturb`.
"""
from __future__ import annotations

import os
import pickle
import sys
from pathlib import Path

os.environ["WANDB_DISABLED"] = "true"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import resolve as _resolve_tissue

from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats
from geneformer import TOKEN_DICTIONARY_FILE
from geneformer.device import get_device


def warn_if_multiproc(name: str, nproc: int, max_ncells: int, n_genes: int) -> None:
    if nproc <= 1:
        return
    print(
        f"[WARN] {name}: nproc={nproc} (>1). Only safe with datasets<5 (pin "
        "datasets==4.0.0; 5.x hangs perturb_data's dataset.map). Workers also add "
        f"RAM competing with the MPS forward pass. With max_ncells={max_ncells} "
        f"and {n_genes} perturbed gene(s), keep max_ncells modest to avoid OOM. "
        "nproc=1 is the verified-stable default.",
        file=sys.stderr,
        flush=True,
    )


def _check_datasets_version() -> None:
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
    if combos > 0:
        import math

        n_variants = math.comb(n_genes, combos + 1) if n_genes >= combos + 1 else 1
    else:
        n_variants = max(n_genes, 1)
    return (n_cells * n_variants * seq_len * 10) / 1e9


NPROC = int(os.environ.get("IS_NPROC", "1"))
_check_datasets_version()
print("device:", get_device(), flush=True)

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

ISP_DIR = ROOT / "results" / "isp" / "early_ad_microglia"
ISP_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- states
# Disease axis: early AD -> early WT (healthy). alt_states left empty so stats
# report only Shift_to_goal_end.
EARLY_TP = ["AD3m", "AD4p5m", "AD6m", "WT3m", "WT4p5m", "WT6m"]
POOL_CELLTYPES = [
    # brain-resident innate immunity — Microglia only (microglia alone run)
    "Microglia",
]

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

# ---------------------------------------------------------------- gene list
# Hypothesis-driven fixed list (human Ensembl ENSG). All verified to be present
# in the token dictionary. Microglia/BAM (DAM / neuroinflammation) + AD GWAS
# risk + inflammatory cytokines.
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

genes_to_perturb = (
    AD_RISK + DAM_MICROGLIA + INFLAMMATION
)
# dedupe (TREM2 / APOE appear in both AD_RISK and DAM_MICROGLIA)
genes_to_perturb = list(dict.fromkeys(genes_to_perturb))
print(f"candidate genes (hypothesis list): {len(genes_to_perturb)}", flush=True)

# ---------------------------------------------------------------- gene presence check
# perturb_data raises "No cells in dataset contain genes_to_perturb" if NO cell
# in the filtered pool expresses any gene. Drop genes absent from the pool cells
# (with a warning) rather than fail. Genes still expressed in a subset survive.
from datasets import load_from_disk

tokens = load_from_disk(str(TOKENIZED))
with open(TOKEN_DICTIONARY_FILE, "rb") as f:
    token_dictionary = pickle.load(f)
tok2gene = {v: k for k, v in token_dictionary.items() if k.startswith("ENSG")}

# genes present in the filtered pool cells
present = set()
for ex in tokens:
    if ex["celltype"] not in POOL_CELLTYPES or ex["samples4"] not in EARLY_TP:
        continue
    for t in ex["input_ids"]:
        g = tok2gene.get(t)
        if g is not None:
            present.add(g)

keep = [g for g in genes_to_perturb if g in present]
dropped = [g for g in genes_to_perturb if g not in present]
if dropped:
    print(
        f"[WARN] dropping {len(dropped)} genes not expressed in pool cells: "
        f"{dropped}",
        flush=True,
    )
genes_to_perturb = keep
print(f"genes_to_perturb (expressed in pool): {len(genes_to_perturb)}", flush=True)

# num_classes: total celltype classes in the trained classifier
celltypes = sorted(set(tokens["celltype"]))
NUM_CLASSES = len(celltypes)

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
# PER-GENE early-detection screening. The Geneformer `"all"` mode (each gene in
# each cell deleted one at a time) is O(cells x genes_expressed) forward passes
# and is impractical here (~10+ hrs for 200 cells). Instead we run a SINGLE-gene
# group perturbation per candidate gene: `genes_to_perturb=[gene]` produces a
# per-gene goal_state_shift whose Shift_to_goal_end is exactly the early-detection
# metric we want (deleting this gene in early-AD cells -> how far toward early-WT).
# Each single-gene run is fast (one perturbed variant per cell).
import pandas as pd

MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "200"))
print(f"per-gene screening: {len(genes_to_perturb)} genes x max_ncells={MAX_CELLS}", flush=True)

results_rows = []
for gene in genes_to_perturb:
    print(f"\n=== gene: {gene} ===", flush=True)
    isp = InSilicoPerturber(
        perturb_type="delete",
        perturb_rank_shift=None,
        genes_to_perturb=[gene],
        combos=0,
        anchor_gene=None,
        model_type="CellClassifier",
        num_classes=NUM_CLASSES,
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
        f"isp_perturb_{gene}",
    )

    ispstats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed=[gene],
        combos=0,
        anchor_gene=None,
        cell_states_to_model=cell_states_to_model,
        model_version="V2",
    )
    ispstats.get_stats(
        str(ISP_DIR),
        None,
        str(ISP_DIR),
        f"isp_stats_{gene}",
    )

    gene_csv = ISP_DIR / f"isp_stats_{gene}.csv"
    if gene_csv.exists():
        gdf = pd.read_csv(gene_csv)
        if "Shift_to_goal_end" in gdf.columns and len(gdf) > 0:
            results_rows.append(
                {
                    "Ensembl_ID": gene,
                    "Shift_to_goal_end": float(gdf["Shift_to_goal_end"].iloc[0]),
                }
            )
            print(f"{gene}: Shift_to_goal_end = {gdf['Shift_to_goal_end'].iloc[0]:.4f}", flush=True)
        else:
            print(f"{gene}: no Shift_to_goal_end in stats (possibly no cells expressed it)", flush=True)
    else:
        print(f"{gene}: stats CSV not found", flush=True)

# ---------------------------------------------------------------- step 4: aggregate candidates
if results_rows:
    cand_df = pd.DataFrame(results_rows).sort_values(
        "Shift_to_goal_end", ascending=False
    )
    cand_out = ISP_DIR / "isp_stats_candidates.csv"
    cand_df.to_csv(cand_out, index=False)
    print("\n=== CANDIDATE EARLY-DETECTION RANKING (positive Shift_to_goal_end = "
          "deleting this gene in early-AD pushes toward early-WT) ===", flush=True)
    print(cand_df.to_string(index=False), flush=True)
    print(f"\nsaved to {cand_out}", flush=True)
else:
    print("no per-gene results collected", flush=True)

print("=== IS PERTURBATION (EARLY AD BRAIN) DONE ===", flush=True)