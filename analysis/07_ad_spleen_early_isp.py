#!/usr/bin/env python
"""In silico perturbation (deletion mode) for early AD detection — AD_spleen.

Disease-state in silico perturbation with a hypothesis-driven gene list to find
genes whose deletion pushes **early-AD spleen immune cells toward healthy(WT)** —
i.e. candidate early drivers / spleen-based biomarker targets.

Ported from 07_ad_blood_early_isp.py (verified working on this host) and
adapted for the AD_spleen dataset (25 cell types, 5xFAD-family spleen).

Design:

  state_key    = "disease"
  start_state  = "AD"
  goal_state   = "WT"
  alt_states   = []     (statistically more robust with single-gene loop)

  filter_data  = spleen immune-cell pool + one timepoint at a time:
                 {"celltype": SPL_IMMUNE_POOL, "samples4": ["AD3m","WT3m"]} etc.

  model/variant: fine-tuned Geneformer cell classifier (CellClassifier),
                 emb_mode="cls" (V2 uses <cls>), model_version="V2"

  perturb_type = "delete", combos=0, one gene at a time (gene loop), nproc=1.
  A gene whose deletion yields a large positive `Shift_to_goal_end` pushes
  early-AD spleen cells toward WT — a candidate early driver / target.

  Timepoints evaluated SEPARATELY (3m / 4.5m / 6m) and compared, so the user
  can rank genes by how strongly they act at the EARLIEST timepoint (3m) —
  directly addressing early detection of AD onset.

Per-gene loop rationale (inherited from 07_ad_spleen_early_isp/07_ad_spleen_early_isp/07_ad_spleen_early_isp): with combos=0 a *list* of
genes_to_perturb makes perturb_data filter for cells containing ALL listed
tokens simultaneously, which is empty once the pool mixes cell-type-specific
genes. The verified-correct approach is combos=0 with ONE gene per run,
reusing the shared state embeddings across genes.

Compute constraints (inherited from 07/07_ad_spleen_early_isp/07_ad_spleen_early_isp/07_ad_spleen_early_isp):
- datasets must be < 5 (4.0.0) or perturb_data's dataset.map hangs.
- nproc=1 (verified-stable); max_ncells set by IS_MAX_CELLS (RAM estimate printed).
- stats `genes_perturbed` restricted to the single gene actually perturbed.
- TRAINED_MODEL_PATH.txt in the runs dir may be stale (written on another
  machine); resolve the fine-tuned classifier by an is_dir() check on the
  referenced path, else fall back to the pretrained V2-104M. Set
  IS_CELLCLASSIFIER_DIR to force a specific checkpoint dir.
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
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

# Fine-tuned spleen cell classifier (from 06_finetune.py). The local
# TRAINED_MODEL_PATH.txt may be stale (written on the original host). Resolve
# by is_dir() check on the referenced path; force a specific dir with
# IS_CELLCLASSIFIER_DIR. Fall back to pretrained V2-104M inference otherwise.
RUN_DIR = ROOT / "runs"
CELLCLASSIFIER_DIR = None
_forced = os.environ.get("IS_CELLCLASSIFIER_DIR")
if _forced and Path(_forced).is_dir():
    CELLCLASSIFIER_DIR = _forced
    print("using IS_CELLCLASSIFIER_DIR:", CELLCLASSIFIER_DIR, flush=True)
elif (RUN_DIR / "TRAINED_MODEL_PATH.txt").exists():
    _p = (RUN_DIR / "TRAINED_MODEL_PATH.txt").read_text().strip()
    if Path(_p).is_dir():
        CELLCLASSIFIER_DIR = _p
        print("using fine-tuned cell classifier:", CELLCLASSIFIER_DIR, flush=True)
    else:
        print(f"[WARN] stale TRAINED_MODEL_PATH.txt -> {_p} (not a dir)", flush=True)

if CELLCLASSIFIER_DIR is None:
    # search for a local ksplit checkpoint under runs/
    cands = sorted(RUN_DIR.rglob("ksplit*"), key=lambda p: p.stat().st_mtime, reverse=True)
    if cands:
        CELLCLASSIFIER_DIR = str(cands[0])
        print("using newest local ksplit checkpoint:", CELLCLASSIFIER_DIR, flush=True)
    else:
        CELLCLASSIFIER_DIR = str(GF_ROOT / MODEL_NAME)
        print("no fine-tuned classifier found; using pretrained model", flush=True)

# ---------------------------------------------------------------- config
MODEL_TYPE = "CellClassifier"

# Spleen immune-cell pool. Low-signal / non-immune / lineage-noise cells dropped:
# Erythroblasts, Megakaryocytes, Basophils (erythroid-megakaryocytic lineage
# noise) and Fibroblasts / endothelial cells (stroma / vasculature).
SPL_IMMUNE_POOL = [
    # myeloid (AD-risk gene carriers: TREM2/TYROBP/CD33/CSF1R/APOE...)
    "Ly6c.high.classical.Monocytes",
    "Ly6c.low.nonclassical.Monocytes",
    "Macrophages",
    "Neutrophils",
    "DCs",
    "cDC.1",
    "cDC.2",
    "pDCs",
    "Migratory.DCs",
    # adaptive / lymphoid
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
# 3m is the earliest; genes ranked by 3m Shift_to_goal_end = earliest-onset driver.
TIMEPOINTS = [
    ("3m", ["AD3m", "WT3m"]),
    ("4p5m", ["AD4p5m", "WT4p5m"]),
    ("6m", ["AD6m", "WT6m"]),
]

# Hypothesis-driven gene list (shared across timepoints; symbol -> ENSG via
# ENSEMBL_DICTIONARY_FILE, then restricted to vocab + presence in each time
# point's AD pool). All symbols verified present in early-AD spleen cells.
HYPOTHESIS_GENES = [
    # AD risk (GWAS, immune-enriched)
    "APOE", "TREM2", "TYROBP", "CLU", "BIN1", "CD33", "ABCA7", "SORL1",
    "INPP5D", "PLCG2", "MEF2C", "SPI1",
    # innate immunity / neuroinflammation
    "C1QA", "C1QB", "C1QC", "C3", "CD68", "CSF1R", "AIF1", "ITGAM", "ITGAX",
    "SPP1", "CCL2", "IL1B", "TNF", "TLR4", "NLRP3", "LYZ", "APOC1", "CD74",
    "IRF7", "ISG15", "S100A8", "S100A9", "MS4A7", "THBS1",
    # blood-immune / adaptive
    "CD3E", "CD3D", "CD4", "CD8A", "CD19", "MS4A1", "NKG7", "KLRD1", "PRF1",
    "CSF3R", "GCA", "CEACAM1", "FGR", "ITGA2B", "VWF", "CD34", "XCR1",
    "IL3RA", "CLEC9A",
]

ISP_DIR = ROOT / "results" / "isp" / "early_ad_spleen"
ISP_DIR.mkdir(parents=True, exist_ok=True)

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
NUM_CLASSES = len(celltypes)
print(f"{NUM_CLASSES} cell types in dataset", flush=True)

# pre-scan once: gene -> set of (disease, samples4) expressed (in SPL_IMMUNE_POOL
# cells). Presence is tracked per disease so each gene is required to be
# expressed in START_STATE (AD) cells (perturb_data filters start cells for the
# perturbed gene; a gene absent from all start cells hangs/errors).
gene_presence = {}   # ensg -> set((disease, samples4))
for ex in tokens:
    if ex["celltype"] not in SPL_IMMUNE_POOL:
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
    # require expression in START-STATE (AD) cells of this timepoint
    start_key = ("AD", samples[0])
    genes_tp = [
        (s, e) for s, e in gene_name2ensg.items()
        if start_key in gene_presence.get(e, set())
    ]
    # optional smoke subset: IS_MAX_GENES limits how many genes run
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
        CELLCLASSIFIER_DIR,
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
            CELLCLASSIFIER_DIR,
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
    early = combined[combined["Timepoint"] == "3m"].set_index("Gene")["Shift_to_goal_end"]
    combined["Shift_3m"] = combined["Gene"].map(early)
    combined["Early_rank"] = combined["Gene"].map(early.rank(ascending=False).to_dict())
    combined = combined.sort_values("Shift_3m", ascending=False)
    out = ISP_DIR / "ad_spleen_early_isp_stats_combined.csv"
    combined.to_csv(out, index=False)
    print(f"\ncombined stats -> {out}")
    print(combined.to_string(index=False))
else:
    print("[WARN] no per-gene stats.csv produced")

print("=== IS PERTURBATION (EARLY AD SPLEEN) DONE ===", flush=True)
