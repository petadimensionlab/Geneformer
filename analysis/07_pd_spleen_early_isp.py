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
import pickle
import sys
from pathlib import Path

os.environ["WANDB_DISABLED"] = "true"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import resolve as _resolve_tissue

from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats
from geneformer import TOKEN_DICTIONARY_FILE, ENSEMBL_DICTIONARY_FILE
from geneformer.device import get_device


def warn_if_multiproc(name: str, nproc: int, max_ncells: int, n_genes: int) -> None:
    if nproc <= 1:
        return
    print(
        f"[WARN] {name}: nproc={nproc} (>1). Only safe with datasets<5 (pin "
        "datasets==4.0.0; 5.x hangs perturb_data's dataset.map). Workers also add "
        f"RAM competing with the forward pass. With max_ncells={max_ncells} "
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
MODEL_DIR = GF_ROOT / MODEL_NAME
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

# ---------------------------------------------------------------- config
# Model: Pretrained Geneformer-V2-104M (fine-tuned cell classifier on this
# tissue is not available — 0-byte weights). Native "Pretrained" mode supports
# state-embedding extraction and in silico deletion without a label head.
MODEL_TYPE = "Pretrained"
NUM_CLASSES = 0

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

# Hypothesis-driven gene list (shared across timepoints; symbol -> ENSG via
# ENSEMBL_DICTIONARY_FILE, then restricted to vocab + presence in each time
# point's pool). PD core + lysosomal/autophagy + neuroinflammation/immune.
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

ISP_DIR = ROOT / "results" / "isp" / "early_pd"
ISP_DIR.mkdir(parents=True, exist_ok=True)

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
with open(TOKEN_DICTIONARY_FILE, "rb") as f:
    token_dictionary = pickle.load(f)
tok2gene = {v: k for k, v in token_dictionary.items() if k.startswith("ENSG")}
with open(ENSEMBL_DICTIONARY_FILE, "rb") as f:
    name_id = pickle.load(f)  # gene symbol -> ENSG

gene_name2ensg = {}
for sym in HYPOTHESIS_GENES:
    ensg = name_id.get(sym)
    if ensg is None:
        print(f"[skip] {sym}: no ENSG in dictionary", flush=True)
        continue
    if ensg not in token_dictionary:
        print(f"[skip] {sym} ({ensg}): not in model vocab", flush=True)
        continue
    gene_name2ensg[sym] = ensg
print(f"candidate genes (in vocab): {len(gene_name2ensg)}", flush=True)
ensg2name = {v: k for k, v in name_id.items()}

# ---------------------------------------------------------------- data presence
from datasets import load_from_disk  # noqa: E402

tokens = load_from_disk(str(TOKENIZED))
celltypes = sorted(set(tokens["celltype"]))
print(f"{len(celltypes)} cell types in dataset", flush=True)

# pre-scan once: gene -> set of (disease, samples4) expressed (in IMMUNE_POOL
# cells). Presence is tracked per disease so we can require each gene to be
# expressed in START_STATE cells (perturb_data filters start-state cells for
# the perturbed gene; a gene absent from all start cells hangs/errors).
gene_presence = {}   # ensg -> set((disease, samples4))
for ex in tokens:
    if ex["celltype"] not in IMMUNE_POOL:
        continue
    key = (ex["disease"], ex["samples4"])
    for t in ex["input_ids"]:
        g = tok2gene.get(t)
        if g is not None:
            gene_presence.setdefault(g, set()).add(key)
print("per-gene presence pre-scan done", flush=True)

# ---------------------------------------------------------------- run per timepoint
results = []  # [{timepoint, gene, ensg, shift}]

for tp, samples in TIMEPOINTS:
    print(f"\n===== TIMEPOINT {tp} (start/goal: {samples}) =====", flush=True)
    # require expression in START-STATE (PF) cells of this timepoint
    start_key = ("PF", samples[0])
    genes_tp = [
        (s, e) for s, e in gene_name2ensg.items()
        if start_key in gene_presence.get(e, set())
    ]
    # optional smoke subset: IS_MAX_GENES limits how many genes run
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

    # step 1: state embs (computed once per timepoint, shared across genes)
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

    # step 2+3: per-gene perturbation + stats (reuse shared state_embs)
    for gene_name, ensg in genes_tp:
        gene_dir = ISP_DIR / tp / gene_name
        gene_dir.mkdir(parents=True, exist_ok=True)
        _ram_est = estimate_perturb_ram(MAX_CELLS, 1, 0)
        print(
            f"\n  [{tp}/{gene_name}] max_ncells={MAX_CELLS}, est. RAM ~{_ram_est:.2f} GB",
            flush=True,
        )

        isp = InSilicoPerturber(
            perturb_type="delete",
            perturb_rank_shift=None,
            genes_to_perturb=[ensg],
            combos=0,
            anchor_gene=None,
            model_type=MODEL_TYPE,
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
            str(MODEL_DIR),
            str(TOKENIZED),
            str(gene_dir),
            "isp_perturbation",
        )

        ispstats = InSilicoPerturberStats(
            mode="goal_state_shift",
            genes_perturbed=[ensg],
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

        csv = gene_dir / "isp_stats.csv"
        if csv.exists():
            import pandas as pd  # noqa: PLC0415

            df = pd.read_csv(csv, index_col=0)
            shift = None
            if len(df) == 1:
                shift = df["Shift_to_goal_end"].iloc[0] if "Shift_to_goal_end" in df.columns else None
            results.append({
                "Timepoint": tp,
                "Gene": gene_name,
                "Ensembl_ID": ensg,
                "Shift_to_goal_end": shift,
            })

# ---------------------------------------------------------------- combine results
import pandas as pd  # noqa: PLC0415

if results:
    combined = pd.DataFrame(results)
    early = combined[combined["Timepoint"] == "6m"].set_index("Gene")["Shift_to_goal_end"]
    combined["Shift_6m"] = combined["Gene"].map(early)
    combined["Early_rank"] = combined["Gene"].map(early.rank(ascending=False).to_dict())
    combined = combined.sort_values("Shift_6m", ascending=False)
    out = ISP_DIR / "pd_spleen_early_isp_stats_combined.csv"
    combined.to_csv(out, index=False)
    print(f"\ncombined stats -> {out}")
    print(combined.to_string(index=False))
else:
    print("[WARN] no per-gene stats.csv produced")

print("=== IS PERTURBATION (EARLY PD) DONE ===", flush=True)