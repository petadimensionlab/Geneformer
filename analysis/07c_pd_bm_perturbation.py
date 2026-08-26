#!/usr/bin/env python
"""In silico perturbation (delete) for early Parkinson's detection in BONE MARROW.

Adapted from 07b_pd_perturbation.py (lymph node) for PD_BM.

Target: bone-marrow monocyte / APC lineages; model the early pathological shift
   start = WT 6m  (healthy baseline)
   goal  = PF 6m  (alpha-synuclein preformed-fibril seeding = early pathology)
   alt   = mo 6m / tgPD 6m

For each target cell type we run the tutorial's three steps:
  1. EmbExtractor.get_state_embs          (state embeddings per cell)
  2. InSilicoPerturber.perturb_data       (simulated single-gene deletion)
  3. InSilicoPerturberStats.get_stats     (goal_state_shift -> per-gene table)

Design choices / ported learnings from 07b:
- MPS host; datasets must be <5 (datasets==4.0.0) or perturb_data's map hangs.
- nproc defaults to 1 (workers compete with the MPS forward pass).
- max_ncells kept modest (200) to bound the materialised perturbation dataset.
- genes are a fixed PD/immune target list, never "all".
- Each gene is deleted INDIVIDUALLY and isolated in its OWN subdirectory so
  InSilicoPerturberStats.read_dictionaries sees only that gene (per-gene shift).
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

os.environ["WANDB_DISABLED"] = "true"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import resolve as _resolve_tissue

from geneformer import EmbExtractor, InSilicoPerturber, InSilicoPerturberStats
from geneformer.device import get_device


# ------------------------------------------------------------------- targets
# PD-early-detection gene targets (Ensembl IDs, verified against
# gene_name_id_dict_gc104M.pkl).  InSilicoPerturber expects Ensembl IDs.
TARGET_GENES = [
    # PD core
    "ENSG00000145335",  # SNCA
    "ENSG00000188906",  # LRRK2
    "ENSG00000177628",  # GBA1
    "ENSG00000069329",  # VPS35
    "ENSG00000116288",  # PARK7
    # immune sensors / antigen presentation
    "ENSG00000170458",  # CD14
    "ENSG00000136869",  # TLR4
    "ENSG00000204287",  # HLA-DRA
    "ENSG00000196126",  # HLA-DRB1
    "ENSG00000196735",  # HLA-DQA1
    "ENSG00000011600",  # TYROBP
    "ENSG00000204472",  # AIF1
    "ENSG00000182578",  # CSF1R
    "ENSG00000173369",  # C1QB
    "ENSG00000136235",  # GPNMB
    # adaptive / cytotoxic
    "ENSG00000198821",  # CD247
    "ENSG00000010610",  # CD4
    "ENSG00000153563",  # CD8A
    "ENSG00000105374",  # NKG7
    "ENSG00000180644",  # PRF1
    "ENSG00000115523",  # GNLY
    # expanded: monocyte/microglia & PD-risk immune genes
    "ENSG00000169896",  # ITGAM (CD11b) - monocyte integrin
    "ENSG00000168329",  # CX3CR1 - monocyte/microglia fractalkine receptor
    "ENSG00000203747",  # FCGR3A (CD16) - Fc receptor
    "ENSG00000095970",  # TREM2 - microglia/APC signaling
    "ENSG00000125538",  # IL1B - pro-inflammatory
    "ENSG00000232810",  # TNF - pro-inflammatory
    "ENSG00000137462",  # TLR2
    "ENSG00000129226",  # CD68
    "ENSG00000090382",  # LYZ (lysozyme)
    "ENSG00000140678",  # ITGAX (CD11c)
    "ENSG00000153208",  # MERTK
    "ENSG00000105329",  # TGFB1
    "ENSG00000136244",  # IL6
    "ENSG00000108691",  # CCL2
    "ENSG00000162711",  # NLRP3 - inflammasome
    "ENSG00000150782",  # IL18
    "ENSG00000118785",  # SPP1 (osteopontin)
    "ENSG00000105383",  # CD33
    "ENSG00000130203",  # APOE - PD risk
    "ENSG00000143226",  # FCGR2A
    "ENSG00000177575",  # CD163
]

# Target cell types (bone-marrow monocyte / APC lineage; early immune sensing).
APC_CELLTYPES = [
    "Ly6c.high.classical.Monocytes",
    "Ly6c.low.nonclassical.Monocytes",
    "Macrophage",
    "cDCs",
    "pDCs",
    "MDP",
]

STATE_KEY = "state_early"  # synthetic 6m-only column: WT / PF / mo / tg
START_STATE = "WT"
GOAL_STATE = "PF"
ALT_STATES = ["mo", "tg"]

_SAMPLES4_TO_STATE = {
    "WT6m": "WT", "PFF6m": "PF", "mono6m": "mo", "tgPD6m": "tg",
}


def _check_datasets_version() -> None:
    """Warn loudly if datasets>=5 (known perturb_data map hang)."""
    import datasets

    major = int(datasets.__version__.split(".")[0])
    if major >= 5:
        print(
            f"[WARN] datasets=={datasets.__version__} installed. "
            "InSilicoPerturber.perturb_data hangs on dataset.map with datasets>=5. "
            "Pin datasets==4.0.0:  uv pip install datasets==4.0.0",
            file=sys.stderr, flush=True,
        )


def estimate_perturb_ram(n_cells: int, n_genes: int, seq_len: int = 4096) -> float:
    """Rough RAM (GB) perturb_data needs (~n_cells * n_variants * seq_len * 10)."""
    n_variants = max(n_genes, 1)
    return (n_cells * n_variants * seq_len * 10) / 1e9


def build_early_dataset(tokenized_path: Path, early_path: Path, nproc: int) -> Path:
    """Create a 6m-only dataset with a synthetic disease-state column.

    Kept separate from the shared PD_BM.dataset so 07_in_silico_perturbation.py
    (which uses ALL cells + celltype states) is unaffected.
    """
    if early_path.exists():
        print(f"reusing early dataset {early_path}", flush=True)
        return early_path

    from datasets import load_from_disk

    print(f"building 6m-early dataset from {tokenized_path}...", flush=True)
    ds = load_from_disk(str(tokenized_path))

    def add_state(row):
        # use a concrete string (never None) so Arrow infers a string column;
        # non-6m rows get "other" and are dropped by the filter below.
        row[STATE_KEY] = _SAMPLES4_TO_STATE.get(row["samples4"], "other")
        return row

    ds = ds.map(add_state, num_proc=nproc, desc="add state_early")
    keep = sorted(set(_SAMPLES4_TO_STATE.values()))
    ds = ds.filter(lambda r: r[STATE_KEY] in keep, num_proc=nproc, desc="keep-6m")
    ds.save_to_disk(str(early_path))
    print(f"saved early dataset: {len(ds)} cells -> {early_path}", flush=True)
    return early_path


def expressed_targets(ds, ct: str) -> dict[str, int]:
    """Count WT cells of `ct` expressing each target gene.

    Returns {Ensembl_ID: n_cells}. Genes with 0 WT cells are dropped (they
    cannot be perturbed -> 'no cells contain genes_to_perturb' error).
    """
    import collections

    import pickle

    from geneformer.in_silico_perturber import InSilicoPerturber as _ISP

    # reuse the token dictionary to map Ensembl -> token
    tok_file = Path(os.environ.get(
        "GENEFORMER_DIR", "")) / "geneformer" / "token_dictionary_gc104M.pkl"
    with open(tok_file, "rb") as f:
        tok = pickle.load(f)
    target_tokens = {g: tok[g] for g in TARGET_GENES if g in tok}

    wt = ds.filter(
        lambda r: r[STATE_KEY] == START_STATE and r["celltype"] == ct
    )
    cnt = collections.Counter()
    for ids in wt["input_ids"]:
        s = set(ids)
        for g, t in target_tokens.items():
            if t in s:
                cnt[g] += 1
    return {g: cnt[g] for g in TARGET_GENES if cnt[g] > 0}


def run_celltype(
    ct: str,
    n_classes: int,
    cellclassifier_dir: str,
    tokenized: str,
    isp_dir: Path,
    max_ncells: int,
    nproc: int,
    summary: Path,
    ds,
) -> None:
    """Run the in silico perturbation for one cell type.

    Each target gene is deleted INDIVIDUALLY (genes_to_perturb=[gene]); a fixed
    *list* would be treated as a group deletion requiring cells to contain all
    genes at once, which is neither feasible nor per-gene.  The individual
    per-gene perturbation gives each gene's own Shift_to_goal_end.
    """
    cell_states_to_model = {
        "state_key": STATE_KEY,
        "start_state": START_STATE,
        "goal_state": GOAL_STATE,
        "alt_states": ALT_STATES,
    }
    filter_data = {"celltype": [ct]}
    # dots in a celltype break InSilicoPerturberStats output naming (get_stats
    # applies Path(prefix).with_suffix(".csv"), which treats ".DCs_..." as the
    # extension) -> use a dot-free tag for all file prefixes.
    ct_tag = ct.replace(".", "_")

    print(f"\n===== cell type: {ct} =====", flush=True)

    # ---------------- step 1: state embeddings (once per cell type)
    embex = EmbExtractor(
        model_type="CellClassifier",
        num_classes=n_classes,
        filter_data=filter_data,
        max_ncells=1000,
        emb_layer=0,
        summary_stat="exact_mean",
        forward_batch_size=20,
        model_version="V2",
        nproc=nproc,
    )
    state_embs = embex.get_state_embs(
        cell_states_to_model,
        cellclassifier_dir,
        tokenized,
        str(isp_dir),
        f"isp_state_embs_{ct_tag}",
    )

    # ---------------- step 2+3: per-gene delete + goal-state-shift stats
    present = expressed_targets(ds, ct)
    print(f"[{ct}] {len(present)}/{len(TARGET_GENES)} target genes expressed "
          f"in WT cells; perturbing individually", flush=True)

    import pandas as pd

    rows = []
    for g in TARGET_GENES:
        if g not in present:
            # cannot delete a gene no WT cell expresses -> record as missing
            rows.append({"Ensembl_ID": g, "Shift_to_goal_end": None})
            continue

        # Per-gene subdirectory is REQUIRED: InSilicoPerturberStats.read_dictionaries
        # reads ALL *_raw.pickle files in the given directory, so if all genes
        # share one directory every gene's stats would see the union and collapse
        # to the same Gene[0] value. Isolating each gene gives its own shift.
        g_dir = isp_dir / f"{ct_tag}" / g
        g_dir.mkdir(parents=True, exist_ok=True)

        isp = InSilicoPerturber(
            perturb_type="delete",
            perturb_rank_shift=None,
            genes_to_perturb=[g],
            combos=0,
            anchor_gene=None,
            model_type="CellClassifier",
            num_classes=n_classes,
            emb_mode="cls",            # V2
            cell_emb_style="mean_pool",
            filter_data=filter_data,
            cell_states_to_model=cell_states_to_model,
            state_embs_dict=state_embs,
            max_ncells=max_ncells,
            emb_layer=0,
            forward_batch_size=20,
            model_version="V2",
            nproc=nproc,
        )
        prefix = f"isp_pert_{ct_tag}_{g}"
        isp.perturb_data(cellclassifier_dir, tokenized, str(g_dir), prefix)

        ispstats = InSilicoPerturberStats(
            mode="goal_state_shift",
            genes_perturbed=[g],
            combos=0,
            anchor_gene=None,
            cell_states_to_model=cell_states_to_model,
            model_version="V2",
        )
        ispstats.get_stats(str(g_dir), None, str(g_dir), f"isp_stats_{ct_tag}_{g}")

        row = pd.read_csv(g_dir / f"isp_stats_{ct_tag}_{g}.csv").iloc[0].to_dict()
        row["Ensembl_ID"] = g          # single-gene stats omit gene identity
        row.pop("Unnamed: 0", None)    # stray index col from the grouped output
        rows.append(row)
        print(f"  {g}: Shift_to_goal_end="
              f"{row.get('Shift_to_goal_end')}", flush=True)

    df = pd.DataFrame(rows)
    if "celltype" not in df.columns:
        df.insert(0, "celltype", ct)
    else:
        df["celltype"] = ct
    header = not summary.exists()
    df.to_csv(summary, mode="a", header=header, index=False)
    print(f"[{ct}] saved {len(df)} gene rows; combined -> {summary}", flush=True)


def main() -> None:
    _check_datasets_version()

    ROOT, PREFIX = _resolve_tissue()
    GF_ROOT = Path(os.environ["GENEFORMER_DIR"])
    MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
    TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"
    EARLY = ROOT / "tokenized" / f"{PREFIX}_early.dataset"

    # fine-tuned cell classifier (or pretrained fallback)
    RUN_DIR = ROOT / "runs"
    trained_path = RUN_DIR / "TRAINED_MODEL_PATH.txt"
    if trained_path.exists():
        CELLCLASSIFIER = trained_path.read_text().strip()
        print("using fine-tuned cell classifier:", CELLCLASSIFIER, flush=True)
    else:
        CELLCLASSIFIER = str(GF_ROOT / MODEL_NAME)
        print("no fine-tuned classifier; using pretrained model", flush=True)

    ISP_DIR = ROOT / "results" / "isp"
    ISP_DIR.mkdir(parents=True, exist_ok=True)

    NPROC = int(os.environ.get("IS_NPROC", "1"))
    MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "200"))

    # optional subset for a validation run, e.g. IS_CTYPES="MDP"
    ctypes = APC_CELLTYPES
    override = os.environ.get("IS_CTYPES")
    if override:
        ctypes = [c.strip() for c in override.split(",") if c.strip()]

    tokenized_early = build_early_dataset(TOKENIZED, EARLY, NPROC)

    from datasets import load_from_disk

    ds_early = load_from_disk(str(tokenized_early))
    n_classes = len(set(ds_early["celltype"]))

    # reset combined summary each run
    summary = ISP_DIR / "pd_bm_apc_early_isp_summary.csv"
    if summary.exists():
        summary.unlink()

    for ct in ctypes:
        run_celltype(ct, n_classes, CELLCLASSIFIER, str(tokenized_early),
                     ISP_DIR, MAX_CELLS, NPROC, summary, ds_early)

    print("=== PD BM EARLY IS PERTURBATION DONE ===", flush=True)


if __name__ == "__main__":
    main()