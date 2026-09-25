#!/usr/bin/env python
"""In silico perturbation (deletion mode) on the multi-region PD atlas.

Geneformer V2-316M, bf16 (GF_DTYPE=bf16), MPS. Same verified engine as
analysis/07_pd_spleen_early_isp.py and analysis/07f_in_silico_perturbation_AD_BM.py.

Source data / paper framing
---------------------------
Prashant et al. 2024, Sci Data 11:1274 (doi:10.1038/s41597-024-04117-y),
"A multi-region single nucleus transcriptomic atlas of Parkinson's disease";
CELLxGENE collection d5d0df8f-4eee-49d8-a221-a288f50a1590
(2,096,155 nuclei, 97 donors, 5 regions). Design points taken from the paper:
  * five regions capture the subcortical -> cortical spread of Braak PD
    pathology: DMNX (dorsal motor nucleus of the vagus) and GPI (globus pallidus
    interna) are affected EARLY, PMC/DLPFC late, PVC largely spared.
  * alpha-synuclein (SNCA) / Lewy bodies are the neuropathological hallmark.
  * the atlas contains NO substantia nigra, so dopaminergic-identity genes
    (TH, SLC6A3, SLC18A2, DDC) are essentially absent here (<2% of cells in
    every cell type) and the presence scan drops them automatically.

Design
------
  state_key   = "disease"
  start_state = "Parkinson disease"   (perturb these)
  goal_state  = "normal"              (measure the shift toward this)
  alt_states  = []
  cell pools  = the region-specific neuronal types of the paper's EARLIEST
                affected regions, each run separately:
                    DMNX_Neu  (719/965 of its cells come from the DMNX)
                    GPI_Neu   (400/567 of its cells come from the GPI)
  perturb_type = "delete", combos=0, ONE gene per run (a *list* with combos=0
                 makes perturb_data require all tokens inside the same cell),
                 nproc=1.
  A gene with a large positive Shift_to_goal_end pushes PD cells toward the
  control state -> candidate PD driver / intervention target.

Compute guards (all apply on this 64 GB M4 Pro):
  * datasets < 5 (4.0.0) or perturb_data's dataset.map hangs.
  * nproc=1 verified stable.
  * small max_ncells (IS_MAX_CELLS, default 300): perturb_data materialises the
    whole perturbation dataset in RAM.
  * forward_batch_size 16: the Metal kernel aborts on very large matmuls
    (64 x 4096 on the MLM head died with "too large for kernel").
"""
from __future__ import annotations

import collections
import json
import os
import pickle
import sys
from pathlib import Path

os.environ["WANDB_DISABLED"] = "true"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _isp_common import (  # noqa: E402
    check_datasets_version,
    estimate_perturb_ram,
    isp_dir_for,
    load_gene_dicts,
    perturb_one_gene,
    resolve_classifier_dir,
    resolve_experiment,
    warn_if_multiproc,
)
from _resolve_tissue import resolve as _resolve_tissue  # noqa: E402

from geneformer import EmbExtractor, TOKEN_DICTIONARY_FILE  # noqa: E402
from geneformer.device import get_device  # noqa: E402

NPROC = int(os.environ.get("IS_NPROC", "1"))
check_datasets_version()
print("device:", get_device(), flush=True)

ROOT, PREFIX = _resolve_tissue()
GF_ROOT = Path(os.environ["GENEFORMER_DIR"])
MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

CELLCLASSIFIER_DIR = resolve_classifier_dir(ROOT, GF_ROOT, MODEL_NAME)
MODEL_TYPE = "CellClassifier"

DISEASE, TISSUE = "PD", "atlas"
EXPERIMENT = resolve_experiment(DISEASE, TISSUE)
# The null run for criteria.md E1 reuses this engine, so the experiment label can
# carry a suffix (IS_EXPERIMENT_SUFFIX=null -> pd_atlas_null). Without it a null
# run would overwrite the hypothesis run's combined CSV and provenance.
_SUFFIX = os.environ.get("IS_EXPERIMENT_SUFFIX", "").strip()
if _SUFFIX:
    EXPERIMENT = f"{EXPERIMENT}_{_SUFFIX}"
ISP_DIR = isp_dir_for(ROOT, EXPERIMENT)
OUT_CSV = ISP_DIR / f"{EXPERIMENT}_early_isp_stats_combined.csv"

POOLS = [p.strip() for p in os.environ.get("IS_POOLS", "DMNX_Neu,GPI_Neu").split(",") if p.strip()]

# The gene list can be swapped without editing this file (IS_GENES=A,B,C). That is
# how the null distribution is produced: same engine, same conditions, genes
# chosen to be unrelated to PD (see docs/isp/records/pd_atlas.md for the candidate list
# and the presence fractions that decided the final set).
_GENES_ENV = [g.strip() for g in os.environ.get("IS_GENES", "").split(",") if g.strip()]
_GENES_DEFAULT = [
    # alpha-synuclein / monogenic PD
    "SNCA", "LRRK2", "GBA1", "PRKN", "PINK1", "PARK7", "VPS35", "ATP13A2",
    "MAPT", "TMEM175", "GCH1", "DNAJC13", "RAB39B", "CHCHD2", "FBXO7",
    # dopaminergic identity (absent in these 5 regions; kept so the scan proves it)
    "TH", "SLC6A3", "SLC18A2", "DDC", "NR4A2", "LMX1B", "KCNJ6", "CALB1",
    # neuronal / synaptic
    "SNAP25", "SYN1", "NEFL", "NEFM", "NEFH", "GRN", "RBFOX3",
    # lysosome / autophagy
    "LAMP2", "CTSB", "CTSD", "TFEB", "SQSTM1", "MAP1LC3B",
    # neuroinflammation / glia
    "TREM2", "TYROBP", "C1QA", "C1QB", "CD68", "CX3CR1", "ITGAM", "ITGAX",
    "IL6", "TNF", "CCL2", "CXCL10", "HLA-DRA", "S100A8", "S100A9", "LYZ", "CSF1R",
]
HYPOTHESIS_GENES = _GENES_ENV if _GENES_ENV else _GENES_DEFAULT

MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "300"))
EMB_CELLS = int(os.environ.get("IS_EMB_CELLS", "1000"))
MIN_FRAC = float(os.environ.get("IS_MIN_FRAC", "0.20"))
FBS = int(os.environ.get("IS_FORWARD_BATCH", "16"))
warn_if_multiproc("InSilicoPerturber", NPROC, MAX_CELLS, 1)

# ---------------------------------------------------------------- dictionaries
tok2gene, ensg2name, name2ensg = load_gene_dicts()
_token_vocab = set(tok2gene.values())

gene_name2ensg: dict[str, str] = {}
for sym in HYPOTHESIS_GENES:
    ensg = name2ensg.get(sym)
    if ensg is None:
        print(f"[skip] {sym}: no ENSG in dictionary", flush=True)
    elif ensg not in _token_vocab:
        print(f"[skip] {sym} ({ensg}): not in V2 vocab", flush=True)
    else:
        gene_name2ensg[sym] = ensg
print(f"candidate genes in V2 vocab: {len(gene_name2ensg)}/{len(HYPOTHESIS_GENES)}", flush=True)

_tokdict = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
gene_token = {sym: _tokdict[ensg] for sym, ensg in gene_name2ensg.items()}
want_tokens = set(gene_token.values())
tok2sym = {v: k for k, v in gene_token.items()}

# ---------------------------------------------------------------- one-pass prescan
from datasets import load_from_disk  # noqa: E402

tokens = load_from_disk(str(TOKENIZED))
celltypes = sorted(set(tokens["celltype"]))
NUM_CLASSES = len(celltypes)
print(f"{len(tokens)} cells, {NUM_CLASSES} cell types", flush=True)

n_cells_by_key: collections.Counter = collections.Counter()
n_expr_by_key: dict[tuple, collections.Counter] = collections.defaultdict(collections.Counter)
for ex in tokens:
    key = (ex["celltype"], ex["disease"])
    n_cells_by_key[key] += 1
    for t in want_tokens.intersection(ex["input_ids"]):
        n_expr_by_key[key][tok2sym[t]] += 1
print("presence pre-scan done", flush=True)

def frac(pool: str, disease: str, sym: str) -> float:
    n = n_cells_by_key.get((pool, disease), 0)
    return (n_expr_by_key[(pool, disease)][sym] / n) if n else 0.0

# ---------------------------------------------------------------- run per pool
results = []
for pool in POOLS:
    print(f"\n===== CELL POOL {pool} =====", flush=True)
    n_pd = n_cells_by_key.get((pool, "Parkinson disease"), 0)
    n_norm = n_cells_by_key.get((pool, "normal"), 0)
    print(f"cells: PD {n_pd} / normal {n_norm}", flush=True)
    if n_pd == 0 or n_norm == 0:
        print(f"[WARN] {pool}: need both PD and normal cells; skipping", flush=True)
        continue

    frac_map = {sym: frac(pool, "Parkinson disease", sym) for sym in gene_name2ensg}
    genes_pool = sorted(((s, gene_name2ensg[s]) for s in frac_map if frac_map[s] >= MIN_FRAC),
                        key=lambda se: -frac_map[se[0]])
    _mg = os.environ.get("IS_MAX_GENES")
    if _mg:
        genes_pool = genes_pool[: int(_mg)]
        print(f"[smoke] limited to {len(genes_pool)} genes", flush=True)
    print(f"genes with presence >= {MIN_FRAC} in PD {pool} cells: {len(genes_pool)}", flush=True)
    for s, _ in genes_pool:
        print(f"    {s}: {frac_map[s]:.2f}", flush=True)
    if not genes_pool:
        print("[WARN] no gene passed the presence filter", flush=True)
        continue

    cell_states_to_model = {
        "state_key": "disease",
        "start_state": "Parkinson disease",
        "goal_state": "normal",
        "alt_states": [],
    }
    filter_data_dict = {"celltype": [pool]}

    embex = EmbExtractor(
        model_type=MODEL_TYPE,
        num_classes=NUM_CLASSES,
        filter_data=filter_data_dict,
        max_ncells=EMB_CELLS,
        emb_layer=0,
        summary_stat="exact_mean",
        forward_batch_size=FBS,
        model_version="V2",
        nproc=NPROC,
    )
    state_embs_dict = embex.get_state_embs(
        cell_states_to_model,
        str(CELLCLASSIFIER_DIR),
        str(TOKENIZED),
        str(ISP_DIR),
        f"isp_state_embs_{pool}",
    )
    print("state_embs keys:", list(state_embs_dict.keys()), flush=True)

    for gene_name, ensg in genes_pool:
        gene_dir = ISP_DIR / pool / gene_name
        gene_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n  [{pool}/{gene_name}] max_ncells={MAX_CELLS}, "
              f"est. RAM ~{estimate_perturb_ram(MAX_CELLS, 1, 0):.2f} GB", flush=True)
        shift = perturb_one_gene(
            model_dir=str(CELLCLASSIFIER_DIR),
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
            forward_batch_size=FBS,
        )
        if shift is None:
            print(f"[WARN] {gene_name}: no Shift_to_goal_end", flush=True)
        results.append({
            "Celltype": pool, "Gene": gene_name, "Ensembl_ID": ensg,
            "Presence_frac": round(frac_map[gene_name], 3),
            "Shift_to_goal_end": shift,
        })

# ---------------------------------------------------------------- combine
import pandas as pd  # noqa: E402

if not results:
    print("[WARN] no results produced", flush=True)
    sys.exit(1)
df = pd.DataFrame(results).sort_values(["Celltype", "Shift_to_goal_end"],
                                       ascending=[True, False])
df.to_csv(OUT_CSV, index=False)

# Provenance: the ISP ranking inherits whatever classifier it loaded, so record
# exactly which model was used and its training regime (a truncated fine-tune is
# a result-affecting choice and must be visible next to the numbers).
_ckpt_summary = Path(CELLCLASSIFIER_DIR) / "adpd_finetuned_summary.json"
provenance = {
    "experiment": EXPERIMENT,
    "model_name": MODEL_NAME,
    "gf_dtype": os.environ.get("GF_DTYPE", "float32"),
    "cellclassifier_dir": str(CELLCLASSIFIER_DIR),
    "classifier_training_regime": (
        json.loads(_ckpt_summary.read_text()) if _ckpt_summary.exists() else "not available"
    ),
    "state_key": "disease",
    "start_state": "Parkinson disease",
    "goal_state": "normal",
    "cell_pools": POOLS,
    "min_presence_frac": MIN_FRAC,
    "max_ncells": MAX_CELLS,
    "nproc": NPROC,
    "forward_batch_size": FBS,
    "genes_run": sorted({r["Gene"] for r in results}),
}
(ISP_DIR / "isp_provenance.json").write_text(json.dumps(provenance, indent=2, default=str))
print(f"provenance -> {ISP_DIR / 'isp_provenance.json'}", flush=True)
print(f"\ncombined stats -> {OUT_CSV}", flush=True)
print(df.to_string(index=False), flush=True)
print("\n=== PD ATLAS IS PERTURBATION DONE ===", flush=True)
