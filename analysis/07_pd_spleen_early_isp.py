#!/usr/bin/env python
"""In silico perturbation (deletion mode) for early PD detection — PD_spleen.

Disease-state in silico perturbation with a hypothesis-driven gene list to find
genes whose deletion pushes **early-PD(PFF) immune cells toward healthy(WT)** —
i.e. candidate early drivers / protective intervention targets.

Design (confirmed with user; mirrors 07_ad_smallint_early_isp.py):

  state_key    = "disease"
  start_state  = "PF"   (PFF alpha-synuclein fibril injection model = sporadic PD)
  goal_state   = "WT"   (healthy control)
  alt_states   = []     (statistically more robust with single-gene loop)

  filter_data  = immune cell-poole subset + one timepoint at a time:
                 {"celltype": IMMUNE_POOL, "samples4": ["PFF6m","WT6m"]} etc.

  model/variant: Pretrained Geneformer-V2-104M ("Pretrained" mode),
                 emb_mode="cls" (V2 uses <cls>), model_version="V2"

  perturb_type = "delete", combos=0, one gene at a time (gene loop), nproc=1.
  A gene whose deletion yields a large positive `Shift_to_goal_end` pushes
  early-PD immune cells toward healthy — a candidate early driver / target.

  Timepoints are evaluated SEPARATELY (6m / 9m / 12m) and compared, so the
  user can rank genes by how strongly they act at the earliest timepoint.

Per-gene loop rationale (inherited from 07_ad_spleen_early_isp): with combos=0 a *list* of
genes_to_perturb makes perturb_data filter for cells containing ALL listed
tokens simultaneously, which is empty once the pool mixes cell-type-specific
genes. The verified-correct approach is combos=0 with ONE gene per run,
reusing the shared state embeddings across genes.

Compute constraints (inherited from 07/07_ad_spleen_early_isp):
- datasets must be < 5 (4.0.0) or perturb_data's dataset.map hangs.
- nproc=1 (verified-stable); max_ncells set by IS_MAX_CELLS (RAM estimate printed).
- stats `genes_perturbed` restricted to the single gene actually perturbed.
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
    resolve_experiment,
    scan_pool_presence_per_key,
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
MODEL_DIR = GF_ROOT / MODEL_NAME
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

# ---------------------------------------------------------------- config
# Model: Pretrained Geneformer-V2-104M (fine-tuned cell classifier on this
# tissue is not available — 0-byte weights). Native "Pretrained" mode supports
# state-embedding extraction and in silico deletion without a label head.
MODEL_TYPE = "Pretrained"
NUM_CLASSES = 0

DISEASE = "PD"
TISSUE = "spleen"
EXPERIMENT = resolve_experiment(DISEASE, TISSUE)
ISP_DIR = isp_dir_for(ROOT, EXPERIMENT)
OUT_CSV = combined_csv_path(ROOT, EXPERIMENT)

# Immune cell pool (spleen immune compartment; exclude erythroid/megakaryocytic).
IMMUNE_POOL = [
    "Macrophages",
    "Ly6c.high.classical.Monocytes",
    "Ly6c.low.nonclassical.Monocytes",
    "DCs",
    "cDC.1",
    "cDC.2",
    "pDCs",
    "Migratory.DCs",
    "CD4.T.cells",
    "CD8.T.cells",
    "NK.T.cells",
    "GammaDelta.T.cells",
    "Naive.B.cells",
    "Marzinal.zone.B.cells",
    "Plasmablasts_Plasma.cells",
    "NK_ILC1",
]

# Timepoints compared: (label, [samples4 PF/start, samples4 WT/goal])
TIMEPOINTS = [
    ("6m", ["PFF6m", "WT6m"]),
    ("9m", ["PFF9m", "WT9m"]),
    ("12m", ["PFF12m", "WT12m"]),
]

# Hypothesis-driven gene list (symbol -> ENSG via ENSEMBL_DICTIONARY_FILE, then
# restricted to vocab + presence in each time point's pool).
HYPOTHESIS_GENES = [
    # alpha-synuclein / PD core
    "SNCA", "LRRK2", "PARK7", "PINK1", "PRKN", "VPS35", "GBA1", "ATP13A2",
    # lysosome / autophagy
    "LAMP2", "CTSB", "CTSD",
    # neuroinflammation / immune
    "TNF", "IL1B", "IL6", "CXCL8", "CX3CR1", "TYROBP", "TREM2", "TREM1",
    "S100A8", "S100A9", "LYZ", "CD68", "CSF1R", "ITGAM", "ITGAX",
    "HLA-DRA", "FOXP3",
    # complement / chemokines
    "C1QA", "C1QB", "CXCL10", "CCL2", "CCL3", "CD14",
]

MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "300"))
EMB_CELLS = int(os.environ.get("IS_EMB_CELLS", "1000"))
warn_if_multiproc("InSilicoPerturber", NPROC, MAX_CELLS, 1)

# optional smoke subset: limit how many genes run (IS_MAX_GENES) and/or which
# timepoint (IS_TIMEPOINTS, comma-separated labels e.g. "6m")
_smoke_tp = os.environ.get("IS_TIMEPOINTS")
if _smoke_tp:
    wanted = {t.strip() for t in _smoke_tp.split(",")}
    TIMEPOINTS = [tp for tp in TIMEPOINTS if tp[0] in wanted]
    print(f"[smoke] limited timepoints to {[tp[0] for tp in TIMEPOINTS]}", flush=True)

# ---------------------------------------------------------------- lookups
tok2gene, ensg2name, name2ensg = load_gene_dicts()
_token_vocab = set(tok2gene.values())

gene_name2ensg = {}
for sym in HYPOTHESIS_GENES:
    ensg = name2ensg.get(sym)
    if ensg is None:
        print(f"[skip] {sym}: no ENSG in dictionary", flush=True)
        continue
    if ensg not in _token_vocab:
        print(f"[skip] {sym} ({ensg}): not in model vocab", flush=True)
        continue
    gene_name2ensg[sym] = ensg
print(f"candidate genes (in vocab): {len(gene_name2ensg)}", flush=True)

# ---------------------------------------------------------------- data presence
from datasets import load_from_disk  # noqa: E402

tokens = load_from_disk(str(TOKENIZED))
celltypes = sorted(set(tokens["celltype"]))
print(f"{len(celltypes)} cell types in dataset", flush=True)

# pre-scan once: gene -> set of (disease, samples4) expressed (in IMMUNE_POOL
# cells). Presence tracked per disease so we require each gene expressed in
# START_STATE cells (a gene absent from all start cells hangs/errors).
gene_presence = scan_pool_presence_per_key(tokens, tok2gene, IMMUNE_POOL)
print("per-gene presence pre-scan done", flush=True)

# ---------------------------------------------------------------- run per timepoint
results = []  # [{timepoint, gene, ensg, shift}]

for tp, samples in TIMEPOINTS:
    print(f"\n===== TIMEPOINT {tp} (start/goal: {samples}) =====", flush=True)
    start_key = ("PF", samples[0])
    genes_tp = [
        (s, e) for s, e in gene_name2ensg.items()
        if start_key in gene_presence.get(e, set())
    ]
    _max_genes = os.environ.get("IS_MAX_GENES")
    if _max_genes:
        genes_tp = genes_tp[: int(_max_genes)]
        print(f"[smoke] limited to {len(genes_tp)} genes", flush=True)
    print(f"genes expressed in {tp} immune pool: {len(genes_tp)}", flush=True)
    if not genes_tp:
        print("[WARN] no hypothesis genes expressed in this timepoint", flush=True)
        continue

    cell_states_to_model = {
        "state_key": "disease",
        "start_state": "PF",
        "goal_state": "WT",
        "alt_states": [],
    }
    filter_data_dict = {
        "celltype": IMMUNE_POOL,
        "samples4": samples,
    }

    embex = EmbExtractor(
        model_type=MODEL_TYPE,
        num_classes=NUM_CLASSES,
        filter_data=filter_data_dict,
        max_ncells=EMB_CELLS,
        emb_layer=0,
        summary_stat="exact_mean",
        forward_batch_size=64,
        model_version="V2",
        nproc=NPROC,
    )
    state_embs_dict = embex.get_state_embs(
        cell_states_to_model,
        str(MODEL_DIR),
        str(TOKENIZED),
        str(ISP_DIR),
        f"isp_state_embs_{tp}",
    )
    print("state_embs keys:", list(state_embs_dict.keys()), flush=True)

    for gene_name, ensg in genes_tp:
        gene_dir = ISP_DIR / tp / gene_name
        gene_dir.mkdir(parents=True, exist_ok=True)
        _ram_est = estimate_perturb_ram(MAX_CELLS, 1, 0)
        print(
            f"\n  [{tp}/{gene_name}] max_ncells={MAX_CELLS}, est. RAM ~{_ram_est:.2f} GB",
            flush=True,
        )
        shift = perturb_one_gene(
            model_dir=str(MODEL_DIR),
            tokenized=str(TOKENIZED),
            out_dir=gene_dir,
            gene_ensg=ensg,
            cell_states_to_model=cell_states_to_model,
            filter_data=filter_data_dict,
            state_embs_dict=state_embs_dict,
            model_type=MODEL_TYPE,
            num_classes=NUM_CLASSES,
            max_ncells=MAX_CELLS,
            nproc=NPROC,
        )
        if shift is None:
            print(f"[WARN] {gene_name}: no Shift_to_goal_end in stats", flush=True)
        results.append({
            "Timepoint": tp,
            "Gene": gene_name,
            "Ensembl_ID": ensg,
            "Shift_to_goal_end": shift,
        })

# ---------------------------------------------------------------- combine results
write_combined(results, OUT_CSV, rank_tp="6m")

print("=== IS PERTURBATION (EARLY PD) DONE ===", flush=True)