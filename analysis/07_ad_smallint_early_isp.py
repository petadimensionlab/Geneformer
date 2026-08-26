#!/usr/bin/env python
"""In silico perturbation (deletion mode) for early AD detection — AD_smallint.

Disease-state in silico perturbation with a hypothesis-driven gene list.

Design (vs. the dummy 07 run which auto-picked the most abundant *celltype*
states and perturbed a single auto-selected gene):

  state_key    = "disease"
  start_state  = "AD"
  goal_state   = "WT"
  alt_states   = []
  filter_data  = immune + gut-epithelial cell types, early timepoints (3-6m)

  perturb_type = "delete", combos=0, emb_mode="cls" (V2), model_version="V2"

Rationale for per-gene loop: with `combos=0` a *list* of genes_to_perturb makes
perturb_data filter for cells containing ALL listed tokens simultaneously
(`filter_data_by_tokens`: intersection == len(tokens)), which is empty once the
pool mixes cell-type-specific genes (e.g. TREM2 in macrophages + MUC2 in goblet
cells). The verified-correct approach is to run combos=0 with ONE gene at a time
and collect that gene's Shift_to_goal_end. State embeddings are computed once
and shared across genes.

A gene whose deletion yields a large positive `Shift_to_goal_end` pushes early-AD
cells toward early-WT (healthy) — a candidate early driver / intervention target.

Hypothesis-driven fixed gene list (Ensembl IDs verified in the token dictionary):
AD risk (GWAS) + inflammation/immune + gut-barrier (mucins, tight junctions,
antimicrobial peptides, stem/metaplasia). Genes absent from the pool cells are
skipped with a warning.

Compute constraints (inherited from 07, see its docstring):
- datasets must be < 5 (4.0.0) or perturb_data's dataset.map hangs.
- nproc=1 (verified-stable); max_ncells kept modest (RAM estimate printed).
- stats `genes_perturbed` restricted to the single gene actually perturbed.
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
from geneformer import ENSEMBL_DICTIONARY_FILE
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

ISP_DIR = ROOT / "results" / "isp" / "early_ad"
ISP_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------- states
EARLY_TP = ["AD3m", "AD4p5m", "AD6m", "WT3m", "WT4p5m", "WT6m"]
POOL_CELLTYPES = [
    "Macrophages", "CD4.T.cells", "CD8.T.cells", "NK_ILC1",
    "IECs", "Paneth.cells", "Goblet.cells_M.cells", "Epithelial.cells",
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
from datasets import load_from_disk

tokens = load_from_disk(str(TOKENIZED))
with open(TOKEN_DICTIONARY_FILE, "rb") as f:
    token_dictionary = pickle.load(f)
tok2gene = {v: k for k, v in token_dictionary.items() if k.startswith("ENSG")}
with open(ENSEMBL_DICTIONARY_FILE, "rb") as f:
    name_id = pickle.load(f)  # symbol -> ENSG
ensg2name = {v: k for k, v in name_id.items()}

present = set()
for ex in tokens:
    if ex["celltype"] not in POOL_CELLTYPES or ex["samples4"] not in EARLY_TP:
        continue
    for t in ex["input_ids"]:
        g = tok2gene.get(t)
        if g is not None:
            present.add(g)

genes_to_perturb = [g for g in genes_to_perturb if g in present]
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

    isp = InSilicoPerturber(
        perturb_type="delete",
        perturb_rank_shift=None,
        genes_to_perturb=[gene],
        combos=0,
        anchor_gene=None,
        model_type="CellClassifier",
        num_classes=NUM_CLASSES,
        emb_mode="cls",
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
        str(gene_dir),
        "isp_perturbation",
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
        str(gene_dir),
        None,
        str(gene_dir),
        "isp_stats",
    )

# ---------------------------------------------------------------- combine results
import pandas as pd

records = []
for gene in genes_to_perturb:
    gene_name = ensg2name.get(gene, gene)
    csv = ISP_DIR / gene_name / "isp_stats.csv"
    if csv.exists():
        df = pd.read_csv(csv, index_col=0)
        # per-gene stats.csv has a single row: Shift_to_goal_end (and possibly
        # Goal_end_vs_random_pval). Annotate with the gene identity.
        rec = {"Gene": gene_name, "Ensembl_ID": gene}
        for col in df.columns:
            if len(df) == 1:
                rec[col] = df[col].iloc[0]
        records.append(rec)
if records:
    combined = pd.DataFrame(records)
    combined = combined.sort_values("Shift_to_goal_end", ascending=False)
    out = ISP_DIR / "early_ad_isp_stats_combined.csv"
    combined.to_csv(out, index=False)
    print(f"\ncombined stats -> {out}")
    print(combined.to_string(index=False))
else:
    print("[WARN] no per-gene stats.csv produced")

print("=== IS PERTURBATION (EARLY AD) DONE ===", flush=True)
