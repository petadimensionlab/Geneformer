#!/usr/bin/env python
"""Shared helpers for the per-tissue early-disease in-silico-perturbation
screens (AD/PD detection).

The per-organ scripts (07_*_early_isp.py, 07e_*, 07f_*) each ran on a
different host with slightly different output layouts. This module centralises:

1. **Canonical output layout** (unified across all tissues / hosts):
       input/<TISSUE>/results/isp/<experiment>/
           <experiment>_early_isp_stats_combined.csv   # per-gene results
           isp_state_embs*.pkl                          # shared state embs
           <gene>/   or   <tp>/<gene>/                  # per-gene raw stats
   with `<experiment>` = `<disease>_<tissue>` lowercase, e.g. `ad_spleen`,
   `pd_spleen`, `ad_blood`, `ad_smallint`, `ad_brain`, `ad_liver`, `ad_bm`.

2. **Compute guards** (verified on MPS/GB10 hosts):
     - `datasets` must be <5 (5.x hangs `perturb_data`'s dataset.map).
     - `nproc=1` is the verified-stable default; higher nproc needs datasets<5
       AND a small materialised perturbation dataset.
     - `estimate_perturb_ram()` bounds `max_ncells` to avoid OOM.

3. **Gene / model plumbing**: token+ensembl dict loading, fine-tuned-model
   resolution (with stale `TRAINED_MODEL_PATH.txt` fallback to newest local
   `ksplit*` checkpoint), per-gene presence scanning.

4. **Single-gene deletion runner + combined-CSV writer** so the per-tissue
   scripts only configure states/filters/genes and call these.

Behavior (per-gene loop): a *list* of `genes_to_perturb` under `combos=0`
would make `perturb_data` filter for cells containing ALL listed tokens at
once (empty once the pool mixes cell-type-specific genes). The verified
approach is `combos=0` with ONE gene per run; state embeddings are computed
once and shared across genes. A gene whose deletion gives a large positive
`Shift_to_goal_end` pushes the disease cells toward WT — a candidate early
driver / intervention target.
"""
from __future__ import annotations

import math
import os
import pickle
import sys
from pathlib import Path


def warn_if_multiproc(name: str, nproc: int, max_ncells: int, n_genes: int) -> None:
    """Emit a visible warning when nproc>1 so users understand the constraints."""
    if nproc <= 1:
        return
    print(
        f"[WARN] {name}: nproc={nproc} (>1). "
        "This is only safe with datasets<5 (pin datasets==4.0.0; 5.x hangs "
        "InSilicoPerturber.perturb_data's dataset.map). "
        "Workers also add RAM that competes with the forward pass. "
        f"With max_ncells={max_ncells} and {n_genes} perturbed gene(s), keep "
        "max_ncells modest (or reduce genes/combos) to avoid OOM. "
        "nproc=1 is the verified-stable default.",
        file=sys.stderr,
        flush=True,
    )


def check_datasets_version() -> None:
    """Warn loudly if datasets>=5 (known perturb_data map hang)."""
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
    """Rough RAM (GB) `perturb_data` needs to materialise the perturbation dataset.

    Each cell -> n_variants perturbation rows; combos=0 and N genes give
    n_variants ~ N per cell; combos>0 is the binomial C(N, combos+1) which
    explodes fast.
    """
    if combos > 0:
        n_variants = math.comb(n_genes, combos + 1) if n_genes >= combos + 1 else 1
    else:
        n_variants = max(n_genes, 1)
    total_variants = n_cells * n_variants
    # input_ids + perturb_index + length + dataset/arrow/pandas overhead, ~8-12 B/token
    return (total_variants * seq_len * 10) / 1e9


def resolve_experiment(disease: str, tissue: str) -> str:
    """Canonical experiment label `<disease>_<tissue>` (lowercase, no dots)."""
    return f"{disease.lower()}_{tissue.lower().replace('.', '').replace('-', '_')}"


def isp_dir_for(root: Path, experiment: str) -> Path:
    """Canonical ISP output dir: <ROOT>/results/isp/<experiment>/"""
    d = root / "results" / "isp"
    if experiment:
        d = d / experiment
    d.mkdir(parents=True, exist_ok=True)
    return d


def combined_csv_path(root: Path, experiment: str) -> Path:
    """Canonical combined-results CSV inside the experiment's ISP dir."""
    d = isp_dir_for(root, experiment)
    return d / f"{experiment}_early_isp_stats_combined.csv"


# ---------------------------------------------------------------- models / genes


def load_gene_dicts():
    """Load token+ensembl gene dictionaries; returns (tok2gene, ensg2name, name2ensg).

    - tok2gene:  token -> ENSG (from token_dictionary)
    - ensg2name: ENSG -> symbol (from ensembl/gene-name dictionary)
    - name2ensg: symbol -> ENSG
    """
    from geneformer import ENSEMBL_DICTIONARY_FILE, TOKEN_DICTIONARY_FILE

    with open(TOKEN_DICTIONARY_FILE, "rb") as f:
        token_dictionary = pickle.load(f)
    with open(ENSEMBL_DICTIONARY_FILE, "rb") as f:
        name_id = pickle.load(f)  # symbol -> ENSG
    tok2gene = {v: k for k, v in token_dictionary.items() if k.startswith("ENSG")}
    ensg2name = {v: k for k, v in name_id.items()}
    return tok2gene, ensg2name, name_id


def resolve_classifier_dir(root: Path, gf_root: Path, model_name: str) -> str:
    """Return the fine-tuned cell-classifier dir (or pretrained fallback).

    Resolution order:
      1. `IS_CELLCLASSIFIER_DIR` / `CELLCLASSIFIER_DIR` env override (must exist)
      2. `runs/TRAINED_MODEL_PATH.txt` target — but that file is often stale
         (written on the original host, e.g. an absolute Linux path); skip if
         not a dir.
      3. newest local `runs/**/ksplit*` checkpoint
      4. pretrained model dir (GF_ROOT/MODEL_NAME)
    """
    forced = os.environ.get("IS_CELLCLASSIFIER_DIR") or os.environ.get("CELLCLASSIFIER_DIR")
    if forced and Path(forced).is_dir():
        print("using IS_CELLCLASSIFIER_DIR/CELLCLASSIFIER_DIR:", forced, flush=True)
        return forced

    run_dir = root / "runs"
    trained_file = run_dir / "TRAINED_MODEL_PATH.txt"
    if trained_file.exists():
        p = trained_file.read_text().strip()
        if Path(p).is_dir():
            print("using fine-tuned cell classifier:", p, flush=True)
            return p
        print(f"[WARN] stale TRAINED_MODEL_PATH.txt -> {p} (not a dir)", flush=True)

    if run_dir.is_dir():
        cands = sorted(run_dir.rglob("ksplit*"),
                       key=lambda p: p.stat().st_mtime, reverse=True)
        if cands:
            print("using newest local ksplit checkpoint:", cands[0], flush=True)
            return str(cands[0])

    pretrained = str(gf_root / model_name)
    print("no fine-tuned classifier found; using pretrained model", pretrained, flush=True)
    return pretrained


# ---------------------------------------------------------------- presence scan


def scan_pool_presence(tokens, tok2gene, celltype_pool, samples4_filter=None) -> set:
    """Genes (ENSG) expressed in at least one cell of the given pool.

    `samples4_filter` (list) optionally restricts to specific timepoints.
    """
    present = set()
    for ex in tokens:
        if ex["celltype"] not in celltype_pool:
            continue
        if samples4_filter is not None and ex["samples4"] not in samples4_filter:
            continue
        for t in ex["input_ids"]:
            g = tok2gene.get(t)
            if g is not None:
                present.add(g)
    return present


def scan_pool_presence_per_key(tokens, tok2gene, celltype_pool) -> dict:
    """Gene -> set((disease, samples4)) keys where it is expressed in the pool.

    Used for per-timepoint screens so each gene can be required to be expressed
    in START-STATE cells of that timepoint (a gene absent from all start cells
    hangs `perturb_data`).
    """
    gene_presence: dict[str, set] = {}
    for ex in tokens:
        if ex["celltype"] not in celltype_pool:
            continue
        key = (ex["disease"], ex["samples4"])
        for t in ex["input_ids"]:
            g = tok2gene.get(t)
            if g is not None:
                gene_presence.setdefault(g, set()).add(key)
    return gene_presence


# ---------------------------------------------------------------- single-gene runner


def perturb_one_gene(
    *,
    model_dir: str,
    tokenized: str,
    out_dir: Path,
    gene_ensg: str,
    cell_states_to_model: dict,
    filter_data: dict,
    state_embs_dict: dict,
    model_type: str = "CellClassifier",
    num_classes: int,
    max_ncells: int,
    nproc: int,
    forward_batch_size: int = 64,
    emb_layer: int = 0,
    perturb_prefix: str = "isp_perturbation",
    stats_prefix: str = "isp_stats",
    model_version: str = "V2",
) -> float | None:
    """Delete `gene_ensg` (single-gene combos=0) and return Shift_to_goal_end.

    Writes perturb + stats outputs into `out_dir`. Returns the gene's
    `Shift_to_goal_end` or None if no stats CSV / no value.
    """
    from geneformer import InSilicoPerturber, InSilicoPerturberStats

    isp = InSilicoPerturber(
        perturb_type="delete",
        perturb_rank_shift=None,
        genes_to_perturb=[gene_ensg],
        combos=0,
        anchor_gene=None,
        model_type=model_type,
        num_classes=num_classes,
        emb_mode="cls",          # V2 uses CLS token
        cell_emb_style="mean_pool",
        filter_data=filter_data,
        cell_states_to_model=cell_states_to_model,
        state_embs_dict=state_embs_dict,
        max_ncells=max_ncells,
        emb_layer=emb_layer,
        forward_batch_size=forward_batch_size,
        model_version=model_version,
        nproc=nproc,
    )
    isp.perturb_data(model_dir, tokenized, str(out_dir), perturb_prefix)

    ispstats = InSilicoPerturberStats(
        mode="goal_state_shift",
        genes_perturbed=[gene_ensg],
        combos=0,
        anchor_gene=None,
        cell_states_to_model=cell_states_to_model,
        model_version=model_version,
    )
    ispstats.get_stats(str(out_dir), None, str(out_dir), stats_prefix)

    csv = out_dir / f"{stats_prefix}.csv"
    if not csv.exists():
        return None
    import pandas as pd
    df = pd.read_csv(csv, index_col=0)
    if len(df) == 1 and "Shift_to_goal_end" in df.columns:
        return float(df["Shift_to_goal_end"].iloc[0])
    return None


# ---------------------------------------------------------------- combined writer


def write_combined(results: list[dict], out_csv: Path, rank_tp: str | None) -> Path | None:
    """Write the combined per-gene results CSV, ranked by the earliest timepoint.

    `results` rows: {"Timepoint","Gene","Ensembl_ID","Shift_to_goal_end"}.
    If `rank_tp` given, adds `Shift_<rank_tp>` + `Early_rank` columns.
    """
    if not results:
        print("[WARN] no per-gene stats produced", flush=True)
        return None

    import pandas as pd
    combined = pd.DataFrame(results)
    if rank_tp:
        early = combined[combined["Timepoint"] == rank_tp].set_index("Gene")["Shift_to_goal_end"]
        combined[f"Shift_{rank_tp}"] = combined["Gene"].map(early)
        combined["Early_rank"] = combined["Gene"].map(early.rank(ascending=False).to_dict())
        combined = combined.sort_values(f"Shift_{rank_tp}", ascending=False)
    else:
        combined = combined.sort_values("Shift_to_goal_end", ascending=False)
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out_csv, index=False)
    print(f"\ncombined stats -> {out_csv}", flush=True)
    print(combined.to_string(index=False), flush=True)
    return out_csv