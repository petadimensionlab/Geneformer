#!/usr/bin/env python
"""Visualise the AD_BM in silico-deletion results (goal-state shift AD->WT).

Reads the per-gene, per-celltype InSilicoPerturberStats CSVs produced by
`analysis/07f_in_silico_perturbation_AD_BM.py` and writes to
`input/AD_BM/results/report/`:

  * figures/isp_goal_state_shift_by_celltype.png  — grouped bar per gene+celltype
  * figures/isp_goal_state_shift_heatmap.png      — gene x celltype heatmap
  * figures/isp_goal_state_shift_bar.png          — celltype-pooled ranked bar
  * experiments_mapping.csv                        — experiment <-> file index

The figures/tables are consumed by report.md / report.html.
"""
from __future__ import annotations

import glob
import os
import pickle
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(os.environ.get("ADPD_ROOT", "input/AD_BM"))
ISP = ROOT / "results" / "isp"
REPORT = ROOT / "results" / "report"
FIG = REPORT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

RUN_SCRIPT = "analysis/07f_in_silico_perturbation_AD_BM.py"
DISEASE = "AD (Alzheimer's disease)"
TISSUE = "AD_BM (bone marrow)"

CELLTYPE_ORDER = ["Ly6c.high.classical.Monocytes", "Macrophage"]
CELLTYPE_ALIAS = {
    "Ly6c.high.classical.Monocytes": "Ly6c.hi classical Monocytes",
    "Macrophage": "Macrophage",
}

GENE_META = {
    "TYROBP": ("ENSG00000011600", "TYROBP (DAP12/adapter)"),
    "APOE": ("ENSG00000130203", "APOE (apolipoprotein E)"),
    "S100A9": ("ENSG00000163220", "S100A9 (Calprotectin beta)"),
    "S100A8": ("ENSG00000143546", "S100A8 (Calprotectin alpha)"),
    "LGALS3": ("ENSG00000131981", "LGALS3 (Galectin-3)"),
    "CCR2": ("ENSG00000121807", "CCR2 (monocyte chemotaxis)"),
    "CSF1R": ("ENSG00000182578", "CSF1R (macrophage drive)"),
    "IRF8": ("ENSG00000140968", "IRF8 (myeloid master TF)"),
    "SPI1": ("ENSG00000066336", "SPI1 (PU.1 myeloid TF)"),
    "CTSD": ("ENSG00000117984", "CTSD (Cathepsin D)"),
    "CTSB": ("ENSG00000164733", "CTSB (Cathepsin B)"),
}

# % of AD cells expressing each gene (from AD_BM counts matrix, both cell types)
EXPR_PCT_AD = {
    "TYROBP": {"Ly6c.high.classical.Monocytes": 100, "Macrophage": 96},
    "APOE": {"Ly6c.high.classical.Monocytes": 58, "Macrophage": 94},
    "S100A9": {"Ly6c.high.classical.Monocytes": 99, "Macrophage": 97},
    "S100A8": {"Ly6c.high.classical.Monocytes": 92, "Macrophage": 92},
    "LGALS3": {"Ly6c.high.classical.Monocytes": 98, "Macrophage": 82},
    "CCR2": {"Ly6c.high.classical.Monocytes": 96, "Macrophage": 50},
    "CSF1R": {"Ly6c.high.classical.Monocytes": 89, "Macrophage": 81},
    "IRF8": {"Ly6c.high.classical.Monocytes": 82, "Macrophage": 67},
    "SPI1": {"Ly6c.high.classical.Monocytes": 82, "Macrophage": 56},
    "CTSD": {"Ly6c.high.classical.Monocytes": 46, "Macrophage": 60},
    "CTSB": {"Ly6c.high.classical.Monocytes": 95, "Macrophage": 93},
}

CELLTYPE_N_AD = {"Ly6c.high.classical.Monocytes": 3000, "Macrophage": 1015}


def load_shifts() -> pd.DataFrame:
    """Concatenate each per-gene per-celltype stats CSV into one tidy frame."""
    rows = []
    for ct in CELLTYPE_ORDER:
        for csv_path in sorted(glob.glob(str(ISP / ct / "gene_*" / "isp_gene_*.csv"))):
            symbol = Path(csv_path).parent.name.replace("gene_", "")
            ensembl, label = GENE_META[symbol]
            d = pd.read_csv(csv_path)
            rows.append(
                {
                    "celltype": ct,
                    "celltype_alias": CELLTYPE_ALIAS[ct],
                    "gene": symbol,
                    "gene_label": label,
                    "ensembl": ensembl,
                    "shift": float(d["Shift_to_goal_end"].iloc[0]),
                    "expr_pct_ad": EXPR_PCT_AD[symbol][ct],
                    "csv": str(csv_path),
                }
            )
    df = pd.DataFrame(rows)
    return df


def _perturbed_n(celltype: str, symbol: str) -> int:
    """Number of perturbed start-state (AD) cells from the raw pickle."""
    raw = glob.glob(str(ISP / celltype / ("gene_" + symbol) / "*_raw.pickle"))
    if not raw:
        return 0
    try:
        d = pickle.load(open(raw[0], "rb"))
    except Exception:
        return 0
    start = d.get("AD", {})
    for (_tok, _field), vecs in start.items():
        return len(vecs)
    return 0


def plot_grouped_bar(df: pd.DataFrame) -> None:
    genes = list(GENE_META.keys())
    x = np.arange(len(genes))
    w = 0.38
    fig, ax = plt.subplots(figsize=(11, 5.2))
    for i, ct in enumerate(CELLTYPE_ORDER):
        vals = [df.loc[(df["gene"] == g) & (df["celltype"] == ct), "shift"].iloc[0]
                for g in genes]
        ax.bar(x + (i - 0.5) * w, vals, w, label=CELLTYPE_ALIAS[ct].replace(
            "classical ", "cl."))
        for xi, v in zip(x + (i - 0.5) * w, vals):
            ax.text(xi, v + np.sign(v) * 1e-5, f"{v:.3f}", ha="center",
                    va="bottom" if v >= 0 else "top", fontsize=7)
    ax.axhline(0, color="#333", lw=1)
    ax.set_xticks(x)
    ax.set_xticklabels(genes, rotation=0, fontsize=9)
    ax.set_ylabel("Shift to goal (WT) — gene deletion")
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "In silico gene DELETION in AD bone-marrow myeloid cells\n"
        "goal-state shift toward WT (positive = toward healthy)", fontsize=11)
    fig.tight_layout()
    fig.savefig(FIG / "isp_goal_state_shift_by_celltype.png", bbox_inches="tight")
    plt.close(fig)


def plot_heatmap(df: pd.DataFrame) -> None:
    genes = list(GENE_META.keys())
    mat = np.zeros((len(genes), len(CELLTYPE_ORDER)))
    for r, g in enumerate(genes):
        for c, ct in enumerate(CELLTYPE_ORDER):
            mat[r, c] = df.loc[(df["gene"] == g) & (df["celltype"] == ct), "shift"].iloc[0]
    fig, ax = plt.subplots(figsize=(5.5, 6.5))
    clim = np.max(np.abs(mat))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-clim, vmax=clim, aspect="auto")
    ax.set_xticks(range(len(CELLTYPE_ORDER)))
    ax.set_xticklabels([CELLTYPE_ALIAS[c].replace(" classical", "").replace(" classical ", "")
                        for c in CELLTYPE_ORDER], fontsize=9)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=9)
    for c in range(len(CELLTYPE_ORDER)):
        for r in range(len(genes)):
            ax.text(c, r, f"{mat[r, c]:.3f}", ha="center", va="center", fontsize=8,
                    color="white" if abs(mat[r, c]) > 0.5 * clim else "black")
    ax.set_title("Goal-state shift (AD→WT), gene deletion\nred = toward WT (positive)",
                 fontsize=10)
    fig.colorbar(im, ax=ax, shrink=0.8, label="Shift to goal (WT)")
    fig.tight_layout()
    fig.savefig(FIG / "isp_goal_state_shift_heatmap.png", bbox_inches="tight")
    plt.close(fig)


def plot_ranked_bar(df: pd.DataFrame) -> None:
    # pool the two celltypes, sort descending by shift; annotate best celltype
    order = df.sort_values("shift", ascending=False)
    fig, ax = plt.subplots(figsize=(9.5, 4.6))
    colors = ["#2e7d32" if v > 0 else "#c62828" for v in order["shift"]]
    bars = ax.bar(range(len(order)), order["shift"], color=colors, zorder=3)
    for i, (bar, _r) in enumerate(zip(bars, order.itertuples())):
        v = _r.shift
        ax.text(i, v + np.sign(v) * 1e-5,
                f"{_r.gene}\n{_r.celltype_alias.replace('Ly6c.hi classical ', 'Ly6c.hi ') }",
                ha="center", va="bottom" if v >= 0 else "top", fontsize=8)
    ax.axhline(0, color="#333", lw=1)
    ax.set_xticks([])
    ax.set_ylabel("Shift to goal (WT)")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_title("All 22 gene x celltype deletion runs, ranked\n"
                 "(green = positive → more WT-like; red = push away from WT)",
                 fontsize=11)
    ax.legend(handles=[
        plt.Rectangle((0, 0), 1, 1, color="#2e7d32", label="positive → toward WT"),
        plt.Rectangle((0, 0), 1, 1, color="#c62828", label="negative → away from WT"),
    ], frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "isp_goal_state_shift_bar.png", bbox_inches="tight")
    plt.close(fig)


def write_experiments_csv(df: pd.DataFrame) -> None:
    out = pd.DataFrame({
        "run_script": RUN_SCRIPT,
        "disease": DISEASE,
        "organ_tissue": TISSUE,
        "celltype": df["celltype"],
        "start_state": "AD",
        "goal_state": "WT",
        "alt_states": "",
        "perturb_type": "delete",
        "gene": df["gene"],
        "gene_label": df["gene_label"],
        "ensembl": df["ensembl"],
        "expr_pct_AD_cells": df["expr_pct_ad"],
        "max_ncells": 200,
        "n_perturbed_cells_AD": df.apply(
            lambda r: _perturbed_n(r["celltype"], r["gene"]), axis=1),
        "Shift_to_goal_end": df["shift"],
        "csv": df["csv"],
    })
    out = out.sort_values("Shift_to_goal_end", ascending=False).reset_index(drop=True)
    out.to_csv(REPORT / "experiments_mapping.csv", index=False)
    print("wrote", REPORT / "experiments_mapping.csv")


def main() -> None:
    df = load_shifts()
    if df.empty:
        raise SystemExit("no Shift CSVs under " + str(ISP))
    print(df.sort_values("shift", ascending=False).to_string(index=False))
    plot_grouped_bar(df)
    plot_heatmap(df)
    plot_ranked_bar(df)
    write_experiments_csv(df)
    print("figures ->", FIG)


main()