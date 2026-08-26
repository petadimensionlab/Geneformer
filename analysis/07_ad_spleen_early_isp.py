#!/usr/bin/env python
"""In silico perturbation (deletion mode) for early AD detection — AD_spleen.

Disease-state in silico perturbation with a hypothesis-driven gene list to find
genes whose deletion pushes **early-AD spleen immune cells toward healthy(WT)** —
i.e. candidate early drivers / spleen-based biomarker targets.

Design (mirrors 07_ad_blood_early_isp.py):

  state_key    = "disease"
  start_state  = "AD"
  goal_state   = "WT"
  alt_states   = []     (statistically more robust with single-gene loop)

  filter_data  = spleen immune-cell pool + one timepoint at a time:
                 {"celltype": SPL_IMMUNE_POOL, "samples4": ["AD3m","WT3m"]} etc.

  model/variant: fine-tuned Geneformer cell classifier (CellClassifier),
                 emb_mode="cls", model_version="V2"

  perturb_type = "delete", combos=0, one gene at a time (gene loop), nproc=1.
  A gene whose deletion yields a large positive `Shift_to_goal_end` pushes
  early-AD spleen cells toward WT — a candidate early driver / target.

  Timepoints evaluated SEPARATELY (3m / 4.5m / 6m) and compared, so the user
  can rank genes by how strongly they act at the EARLIEST timepoint (3m).

All generic logic (guards, dicts, model resolution, per-gene runner, combined
CSV writer) lives in `_isp_common.py`; this file only configures the
experiment. Compute constraints (datasets<5, nproc=1, modest max_ncells)
also live in `_isp_common`.
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
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

CELLCLASSIFIER_DIR = resolve_classifier_dir(ROOT, GF_ROOT, MODEL_NAME)

# ---------------------------------------------------------------- config
DISEASE = "AD"
TISSUE = "spleen"
EXPERIMENT = resolve_experiment(DISEASE, TISSUE)
ISP_DIR = isp_dir_for(ROOT, EXPERIMENT)
OUT_CSV = combined_csv_path(ROOT, EXPERIMENT)

# Spleen immune-cell pool. Low-signal / non-immune / lineage-noise cells dropped:
# Erythroblasts, Megakaryocytes, Basophils, Fibroblasts, endothelial cells.
SPL_IMMUNE_POOL = [
    "Ly6c.high.classical.Monocytes",
    "Ly6c.low.nonclassical.Monocytes",
    "Macrophages",
    "Neutrophils",
    "DCs",
    "cDC.1",
    "cDC.2",
    "pDCs",
    "Migratory.DCs",
    "CD4.T.cells",
    "CD8.T.cells",
    "Immature.T.cells",
    "Gamma.Delta.T.cells",
    "NKT.cells",
    "NK_ILC1",
    "Naive.Memory.B.cells",
    "Marzinal.zone.B.cells",
    "Plasma.cells",
    "Germinal.Center.B.cells",
]

# Timepoints compared: (label, [samples4 AD/start, samples4 WT/goal]).
TIMEPOINTS = [
    ("3m", ["AD3m", "WT3m"]),
    ("4p5m", ["AD4p5m", "WT4p5m"]),
    ("6m", ["AD6m", "WT6m"]),
]

HYPOTHESIS_GENES = [
    # AD risk (GWAS, immune-enriched)
    "APOE", "TREM2", "TYROBP", "CLU", "BIN1", "CD33", "ABCA7", "SORL1",
    "INPP5D", "PLCG2", "MEF2C", "SPI1",
    # innate immunity / neuroinflammation
    "C1QA", "C1QB", "C1QC", "C3", "CD68", "CSF1R", "AIF1", "ITGAM", "ITGAX",
    "SPP1", "CCL2", "IL1B", "TNF", "TLR4", "NLRP3", "LYZ", "APOC1", "CD74",
    "IRF7", "ISG15", "S100A8", "S100A9", "MS4A7", "THBS1",
    # adaptive
    "CD3E", "CD3D", "CD4", "CD8A", "CD19", "MS4A1", "NKG7", "KLRD1", "PRF1",
    # myeloid / megakaryocyte / platelet marker
    "CSF3R", "GCA", "CEACAM1", "FGR", "ITGA2B", "VWF", "CD34", "XCR1",
    "IL3RA", "CLEC9A",
]

MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "200"))
EMB_CELLS = int(os.environ.get("IS_EMB_CELLS", "1000"))
warn_if_multiproc("InSilicoPerturber", NPROC, MAX_CELLS, 1)

# optional smoke subset: limit how many genes run (IS_MAX_GENES) and/or which
# timepoint (IS_TIMEPOINTS, comma-separated labels e.g. "3m")
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
NUM_CLASSES = len(celltypes)
print(f"{NUM_CLASSES} cell types in dataset", flush=True)

gene_presence = scan_pool_presence_per_key(tokens, tok2gene, SPL_IMMUNE_POOL)
print("per-gene presence pre-scan done", flush=True)

# ---------------------------------------------------------------- run per timepoint
results = []  # [{timepoint, gene, ensg, shift}]

for tp, samples in TIMEPOINTS:
    print(f"\n===== TIMEPOINT {tp} (start/goal: {samples}) =====", flush=True)
    start_key = ("AD", samples[0])
    genes_tp = [
        (s, e) for s, e in gene_name2ensg.items()
        if start_key in gene_presence.get(e, set())
    ]
    _max_genes = os.environ.get("IS_MAX_GENES")
    if _max_genes:
        genes_tp = genes_tp[: int(_max_genes)]
        print(f"[smoke] limited to {len(genes_tp)} genes", flush=True)
    print(f"genes expressed in {tp} AD spleen pool: {len(genes_tp)}", flush=True)
    if not genes_tp:
        print("[WARN] no hypothesis genes expressed in this timepoint", flush=True)
        continue

    cell_states_to_model = {
        "state_key": "disease",
        "start_state": "AD",
        "goal_state": "WT",
        "alt_states": [],
    }
    filter_data_dict = {
        "celltype": SPL_IMMUNE_POOL,
        "samples4": samples,
    }

    embex = EmbExtractor(
        model_type="CellClassifier",
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
        CELLCLASSIFIER_DIR,
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
            model_dir=CELLCLASSIFIER_DIR,
            tokenized=str(TOKENIZED),
            out_dir=gene_dir,
            gene_ensg=ensg,
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
            "Timepoint": tp,
            "Gene": gene_name,
            "Ensembl_ID": ensg,
            "Shift_to_goal_end": shift,
        })

# ---------------------------------------------------------------- combine results
write_combined(results, OUT_CSV, rank_tp="3m")

print("=== IS PERTURBATION (EARLY AD SPLEEN) DONE ===", flush=True)