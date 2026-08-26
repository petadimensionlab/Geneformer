#!/usr/bin/env python
"""In silico perturbation (delete) for EARLY ALZHEIMER detection in LIVER.

Adapted from 07d_pd_blood_perturbation.py (Parkinson, blood) for the AD liver
mouse model (input/AD_liver). Models the EARLY pathological disease shift, NOT
cell-type identity:

    start = WT3m   (age-matched healthy baseline)
    goal  = AD3m   (earliest AD timepoint = the early shift we want to detect)
    alt   = AD6m / AD9m / AD12m   (later progression)

For each target cell type we run the Nature-tutorial three steps:
    1. EmbExtractor.get_state_embs          (state embeddings per cell)
    2. InSilicoPerturber.perturb_data       (simulated single-gene deletion)
    3. InSilicoPerturberStats.get_stats     (goal_state_shift -> per-gene table)

Interpretation: a gene whose DELETION pushes the AD3m cell toward the WT state
(Shift_to_goal_end more negative / toward goal) is a potential early driver /
therapeutic target; genes whose deletion moves the cell away may be protective.

Target cells (liver AD axis): Hepatocytes, Kupffer.cells, Ly6c.high.classical
Monocytes, Endothelial.cells, Macrophages.

Design choices / ported learnings from 07b/07c/07d:
- MPS host; datasets must be <5 (datasets==4.0.0) or perturb_data's map hangs.
- nproc defaults to 1 (workers compete with the MPS forward pass).
- max_ncells kept modest (200) to bound the materialised perturbation dataset.
- genes are a fixed AD/immune target list, never "all".
- Each gene is deleted INDIVIDUALLY and isolated in its OWN subdirectory so
  InSilicoPerturberStats.read_dictionaries sees only that gene (per-gene shift).
- Genes not expressed in the WT3m start-state cells are auto-dropped
  (deleting them would raise "no cells contain genes_to_perturb").

LOCAL CHECKPOINT: TRAINED_MODEL_PATH.txt points to a Linux path that does not
exist on this host. Use the real local fine-tuned CellClassifier:
   input/AD_liver/runs/260823_geneformer_cellClassifier_AD_liver_celltype/ksplit1
Override via CELLCLASSIFIER_DIR env var if the path moves.
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

# ------------------------------------------------------------- targets
# Early-AD liver gene targets (Ensembl IDs, verified against
# gene_name_id_dict_gc104M.pkl).  InSilicoPerturber expects Ensembl IDs.
# AD risk (liver/immune axis) + inflammatory + innate immune sensors.
TARGET_GENES = [
    # AD GWAS risk genes expressed in liver/immune cells
    "ENSG00000130203",  # APOE  - apolipoprotein, AD risk, liver-produced
    "ENSG00000095970",  # TREM2 - AD risk, Kupffer/macrophage signaling
    "ENSG00000120885",  # CLU   - clusterin (apoJ), AD risk, Abeta clearance
    "ENSG00000136717",  # BIN1  - AD risk
    "ENSG00000197601",  # FAR1 - fatty acid desaturase (liver lipid axis)
    # pro-inflammatory cytokines
    "ENSG00000125538",  # IL1B
    "ENSG00000136244",  # IL6
    "ENSG00000232810",  # TNF
    "ENSG00000108691",  # CCL2
    "ENSG00000162711",  # NLRP3 - inflammasome
    # innate immune sensors / phagocytosis
    "ENSG00000105383",  # CD33 - AD risk, innate
    "ENSG00000140678",  # ITGAX (CD11c)
    "ENSG00000170458",  # CD14
    "ENSG00000136869",  # TLR4
    "ENSG00000011600",  # TYROBP
    "ENSG00000182578",  # CSF1R
    "ENSG00000173369",  # C1QB
    "ENSG00000136235",  # GPNMB
    "ENSG00000118785",  # SPP1 (osteopontin)
    "ENSG00000153208",  # MERTK
    "ENSG00000090382",  # LYZ (lysozyme)
]

# Target cell types (liver AD axis).
LIVER_CELLTYPES = [
    "Hepatocytes",
    "Kupffer.cells",
    "Ly6c.high.classical.Monocytes",
    "Endothelial.cells",
    "Macrophages",
]

# disease-state column (built from `disease` + `samples4`), early-detection.
STATE_KEY = "state_early_ad"
START_STATE = "WT"
GOAL_STATE = "AD3m"
ALT_STATES = ["AD6m", "AD9m", "AD12m"]

# sample4 -> state; everything else -> "other" (dropped)
_SAMPLES4_TO_STATE = {
    "WT3m": "WT",
    "AD3m": "AD3m",
    "AD6m": "AD6m",
    "AD9m": "AD9m",
    "AD12m": "AD12m",
}

# Real local fine-tuned CellClassifier.  TRAINED_MODEL_PATH.txt on disk points
# to a Linux path; use the actual local ksplit1 run.
LOCAL_CELLCLASSIFIER = (
    Path(os.environ.get(
        "ADPD_ROOT",
        Path(__file__).resolve().parent.parent / "input" / "AD_liver",
    ))
    / "runs" / "260823_geneformer_cellClassifier_AD_liver_celltype" / "ksplit1"
)


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
    """Create a disease-state dataset from `disease` + `samples4`.

    Only WT3m / AD3m / AD6m / AD9m / AD12m rows are kept; the rest get "other"
    and are dropped.  Kept separate from the shared AD_liver.dataset so
    07_in_silico_perturbation.py (which uses ALL cells + celltype states) is
    unaffected.
    """
    if early_path.exists():
        print(f"reusing early dataset {early_path}", flush=True)
        return early_path

    from datasets import load_from_disk

    print(f"building 3m-early dataset from {tokenized_path}...", flush=True)
    ds = load_from_disk(str(tokenized_path))

    def add_state(row):
        # concrete string (never None) so Arrow infers a string column; rows
        # outside the mapped states get "other" and are dropped below.
        row[STATE_KEY] = _SAMPLES4_TO_STATE.get(row["samples4"], "other")
        return row

    ds = ds.map(add_state, num_proc=nproc, desc="add state_early")
    keep = sorted(set([START_STATE, GOAL_STATE] + ALT_STATES))
    ds = ds.filter(lambda r: r[STATE_KEY] in keep, num_proc=nproc, desc="keep-3m")
    ds.save_to_disk(str(early_path))
    print(f"saved early dataset: {len(ds)} cells -> {early_path}", flush=True)
    return early_path


def expressed_targets(ds, ct: str) -> dict[str, int]:
    """Count WT (start-state) cells of `ct` expressing each target gene.

    Returns {Ensembl_ID: n_cells}. Genes with 0 WT cells are dropped (they
    cannot be perturbed -> 'no cells contain genes_to_perturb' error).
    """
    import collections

    import pickle

    from geneformer.in_silico_perturber import InSilicoPerturber as _ISP

    # reuse the token dictionary to map Ensembl -> token
    tok_file = Path(os.environ.get(
        "GENEFORMER_DIR", "")) / "geneformer" / "token_dictionary_gc104M.pkl"
    if not tok_file.is_file():
        tok_file = Path(__file__).resolve().parent.parent / "geneformer_hf" / "geneformer" / "token_dictionary_gc104M.pkl"
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
    # dots in a celltype break InSilicoPerturberStats output naming -> dot-free tag.
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
    GF_ROOT = Path(os.environ.get("GENEFORMER_DIR") or (Path(__file__).resolve().parent.parent / "geneformer_hf"))
    MODEL_NAME = os.environ.get("GENEFORMER_MODEL", "Geneformer-V2-104M")
    TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"
    EARLY = ROOT / "tokenized" / f"{PREFIX}_early.dataset"

    # fine-tuned cell classifier: prefer the env override, else the local path
    # (TRAINED_MODEL_PATH.txt on disk is a stale Linux path).
    CELLCLASSIFIER = os.environ.get("CELLCLASSIFIER_DIR")
    if not CELLCLASSIFIER:
        CELLCLASSIFIER = str(LOCAL_CELLCLASSIFIER)
    if not Path(CELLCLASSIFIER).is_dir():
        print(
            f"[WARN] {CELLCLASSIFIER} not found; falling back to pretrained "
            f"{GF_ROOT / MODEL_NAME}",
            file=sys.stderr, flush=True,
        )
        CELLCLASSIFIER = str(GF_ROOT / MODEL_NAME)
    print("using cell classifier:", CELLCLASSIFIER, flush=True)

    ISP_DIR = ROOT / "results" / "isp"
    ISP_DIR.mkdir(parents=True, exist_ok=True)

    NPROC = int(os.environ.get("IS_NPROC", "1"))
    MAX_CELLS = int(os.environ.get("IS_MAX_CELLS", "200"))

    # optional subset for a validation run, e.g. IS_CTYPES="Kupffer.cells"
    ctypes = LIVER_CELLTYPES
    override = os.environ.get("IS_CTYPES")
    if override:
        ctypes = [c.strip() for c in override.split(",") if c.strip()]

    tokenized_early = build_early_dataset(TOKENIZED, EARLY, NPROC)

    from datasets import load_from_disk

    ds_early = load_from_disk(str(tokenized_early))
    n_classes = len(set(ds_early["celltype"]))

    # reset combined summary each run
    summary = ISP_DIR / "ad_liver_early_isp_summary.csv"
    if summary.exists():
        summary.unlink()

    for ct in ctypes:
        run_celltype(ct, n_classes, CELLCLASSIFIER, str(tokenized_early),
                     ISP_DIR, MAX_CELLS, NPROC, summary, ds_early)

    print("=== AD LIVER EARLY IS PERTURBATION DONE ===", flush=True)


if __name__ == "__main__":
    main()