#!/usr/bin/env python
"""Visualize the early-PD in silico perturbation results.

Reads the per-gene goal_state_shift summary (result/report/..._full6.csv) and
writes PNG figures under result/report/assets/:

  1. isp_heatmap.png       — celltype x gene matrix of Shift_to_goal_end
     (positive = gene deletion pushes WT toward early-PD PF state)
  2. isp_gene_ranking.png  — consensus ranking (mean shift across cell types)
     with % of cell types positive, colored by direction
  3. isp_celltype_volcano.png — per-cell-type top genes (bar) showing the
     strongest early-detection candidate deletions

A compact combined figure + a per-gene summary table is also emitted.
"""
from __future__ import annotations

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

_REPO = Path(__file__).resolve().parent.parent
OUT = _REPO / "result" / "report" / "assets"
SRC = _REPO / "result" / "report"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "figure.dpi": 160,
    "savefig.bbox": "tight",
    "font.size": 9,
})


def _load() -> pd.DataFrame:
    df = pd.read_csv(SRC / "pd_apc_early_isp_summary_full6.csv")
    # prefer gene names when present
    if "Gene_name" not in df.columns:
        import pickle

        with open(os.environ["GENEFORMER_DIR"] + "/geneformer/gene_name_id_dict_gc104M.pkl", "rb") as f:
            rev = {v: k for k, v in pickle.load(f).items()}
        df["Gene_name"] = df["Ensembl_ID"].map(lambda e: rev.get(e, e) if pd.notna(e) else None)
    return df


def main() -> None:
    df = _load()
    df = df[df["Shift_to_goal_end"].notna()].copy()
    df["Gene_name"] = df["Gene_name"].fillna(df["Ensembl_ID"])

    # ---- 1. heatmap: celltype x gene ----
    pivot = df.pivot_table(index="celltype", columns="Gene_name",
                           values="Shift_to_goal_end")
    # order celltypes by their mean magnitude for readability
    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]
    # limit to genes present in >=3 cell types for a readable heatmap
    mask = pivot.notna().sum() >= 3
    hp = pivot.loc[:, mask]
    fig, ax = plt.subplots(figsize=(max(6, hp.shape[1] * 0.38), 5))
    cmap = sns.diverging_palette(10, 240, s=80, l=50, as_cmap=True)
    sns.heatmap(hp, ax=ax, cmap=cmap, center=0,
                linewidths=0.3, linecolor="white",
                cbar_kws={"label": "Shift_to_goal_end\n(deletion of gene)",
                          "shrink": 0.8})
    ax.set_title("In silico deletion: shift of WT cells toward early-PD (PF6m)\n"
                 "positive = deletion pushes toward early-PD phenotype")
    ax.set_ylabel("Cell type")
    ax.set_xlabel("Gene deleted")
    plt.xticks(rotation=90, fontsize=8)
    fig.savefig(OUT / "isp_heatmap.png")
    plt.close(fig)

    # ---- 2. gene ranking across cell types ----
    agg = df.groupby("Gene_name").agg(
        mean_shift=("Shift_to_goal_end", "mean"),
        max_shift=("Shift_to_goal_end", "max"),
        n_celltypes=("celltype", "nunique"),
        n_positive=("Shift_to_goal_end", lambda s: int((s > 0).sum())),
    ).reset_index()
    agg["pct_positive"] = agg["n_positive"] / agg["n_celltypes"] * 100
    agg = agg.sort_values("mean_shift", ascending=False)
    top = agg.head(15).iloc[::-1]

    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    colors = ["#d73027" if v > 0 else "#4575b4" for v in top["mean_shift"]]
    ax.barh(top["Gene_name"], top["mean_shift"], color=colors)
    ax.axvline(0, color="k", lw=0.8)
    ax.set_xlabel("Mean Shift_to_goal_end across APC cell types")
    ax.set_title("Deletion of these genes moves WT cells toward early-PD (PF6m)\n"
                 "top 15 by mean shift; red = toward PD, blue = away")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUT / "isp_gene_ranking.png")
    plt.close(fig)

    # ---- 3. per-celltype volcano ----
    n_ct = df["celltype"].nunique()
    ncol, nrow = 3, 2
    fig, axes = plt.subplots(nrow, ncol, figsize=(14, 9))
    for ax, (ct, sub) in zip(axes.ravel(), df.groupby("celltype")):
        sub = sub.sort_values("Shift_to_goal_end", ascending=False)
        x = sub["Shift_to_goal_end"]
        y = np.arange(len(sub), 0, -1)
        cols = ["#d73027" if v > 0 else "#4575b4" for v in x]
        ax.barh(sub["Gene_name"], x, color=cols)
        ax.axvline(0, color="k", lw=0.7)
        ax.set_title(ct, fontsize=9)
        ax.set_xlabel("Shift_to_goal_end", fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
    fig.suptitle("Per-cell-type in silico deletion → early-PD shift", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / "isp_celltype_volcano.png")
    plt.close(fig)

    # ---- gene ranking table ----
    agg.round(5).to_csv(OUT / "isp_gene_ranking.csv", index=False)
    print("wrote:", sorted(p.name for p in OUT.glob("isp_*")))


if __name__ == "__main__":
    main()