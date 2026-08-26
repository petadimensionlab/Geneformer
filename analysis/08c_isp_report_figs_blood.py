#!/usr/bin/env python
"""Blood: reconstruct authoritative per-gene ISP table + visualize results.

Reads the per-gene InSilicoPerturberStats CSVs from each per-gene subdirectory
under input/PD_blood/results/isp/ and:
  1. Rebuilds the authoritative summary (each gene isolated in its own subdir
     is the CORRECT source — the mode='a' combined CSV has a column-order bug).
     -> result/report_blood/pd_blood_apc_early_isp_summary_full6.csv
  2. Computes a consensus gene ranking across the 8 blood cell types
     -> result/report_blood/pd_blood_apc_early_isp_gene_ranking_full6.csv
  3. Writes PNG figures under result/report_blood/assets/:
       isp_heatmap.png       celltype x gene matrix of Shift_to_goal_end
       isp_gene_ranking.png  consensus ranking (mean shift) top 15
       isp_celltype_volcano.png per-cell-type per-gene bar ranking

positive Shift_to_goal_end = gene deletion moves WT cells toward early-PD (PF6m)
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
REPORT = _REPO / "result" / "report_blood"
OUT = REPORT / "assets"
ISP = _REPO / "input" / "PD_blood" / "results" / "isp"
os.makedirs(OUT, exist_ok=True)

# ct_tag -> true celltype name (tag = celltype with '.' removed)
CT_TAG_TO_NAME = {
    "Ly6c_high_classical_Monocytes": "Ly6c.high.classical.Monocytes",
    "Ly6c_low_nonclassical_Monocytes": "Ly6c.low.nonclassical.Monocytes",
    "Neutrophils": "Neutrophils",
    "cDCs": "cDCs",
    "pDCs": "pDCs",
    "CD8_T_cells": "CD8.T.cells",
    "CD4_T_cells": "CD4.T.cells",
    "NK_ILC1": "NK_ILC1",
}

plt.rcParams.update({
    "figure.dpi": 160,
    "savefig.bbox": "tight",
    "font.size": 9,
})


def _gene_name_map() -> dict[str, str]:
    import pickle

    p = _REPO / "geneformer_hf" / "geneformer" / "gene_name_id_dict_gc104M.pkl"
    with open(p, "rb") as f:
        name2id = pickle.load(f)
    return {v: k for k, v in name2id.items()}


def _load_authoritative() -> pd.DataFrame:
    """Rebuild per-gene/per-celltype table from the isolated per-gene subdirs.

    Each gene lives in its own subdir (isp/<ct_tag>/<Ensembl>/isp_stats_*.csv),
    so each row is that gene's own goal_state_shift — the authoritative source.
    """
    name2id = _gene_name_map()
    rows = []
    for d in sorted(ISP.iterdir()):
        if not d.is_dir() or d.name.startswith("isp_state"):
            continue
        ct = CT_TAG_TO_NAME.get(d.name, d.name)
        for gd in sorted(d.iterdir()):
            if not gd.is_dir():
                continue
            stats = list(gd.glob("isp_stats_*.csv"))
            if not stats:
                continue
            row = pd.read_csv(stats[0]).iloc[0].to_dict()
            row.pop("Unnamed: 0", None)
            row["celltype"] = ct
            row["Ensembl_ID"] = gd.name
            row["Gene_name"] = name2id.get(gd.name, gd.name)
            rows.append(row)
    return pd.DataFrame(rows)


def main() -> None:
    df = _load_authoritative()
    df = df[df["Shift_to_goal_end"].notna()].copy()

    # ---- 1. heatmap: celltype x gene ----
    pivot = df.pivot_table(index="celltype", columns="Gene_name",
                           values="Shift_to_goal_end")
    pivot = pivot.loc[pivot.mean(axis=1).sort_values(ascending=False).index]
    mask = pivot.notna().sum() >= 3
    hp = pivot.loc[:, mask]
    fig, ax = plt.subplots(figsize=(max(6, hp.shape[1] * 0.38), 5))
    cmap = sns.diverging_palette(10, 240, s=80, l=50, as_cmap=True)
    sns.heatmap(hp, ax=ax, cmap=cmap, center=0,
                linewidths=0.3, linecolor="white",
                cbar_kws={"label": "Shift_to_goal_end\n(deletion of gene)",
                          "shrink": 0.8})
    ax.set_title("In silico deletion in blood cells: shift of WT toward early-PD (PF6m)\n"
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
    ax.set_xlabel("Mean Shift_to_goal_end across blood cell types")
    ax.set_title("Deletion of these genes moves WT blood cells toward early-PD\n"
                 "top 15 by mean shift; red = toward PD, blue = away")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUT / "isp_gene_ranking.png")
    plt.close(fig)

    # ---- 3. per-celltype volcano ----
    n_ct = df["celltype"].nunique()
    ncol = 3
    nrow = int(np.ceil(n_ct / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(15, 5 * nrow))
    axes = np.atleast_1d(axes).ravel()
    for ax, (ct, sub) in zip(axes, df.groupby("celltype")):
        sub = sub.sort_values("Shift_to_goal_end", ascending=False)
        x = sub["Shift_to_goal_end"]
        cols = ["#d73027" if v > 0 else "#4575b4" for v in x]
        ax.barh(sub["Gene_name"], x, color=cols)
        ax.axvline(0, color="k", lw=0.7)
        ax.set_title(ct, fontsize=9)
        ax.set_xlabel("Shift_to_goal_end", fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
    for ax in axes[n_ct:]:
        ax.axis("off")
    fig.suptitle("Per-cell-type in silico deletion to early-PD shift", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(OUT / "isp_celltype_volcano.png")
    plt.close(fig)

    # ---- write data ----
    df.sort_values(["celltype", "Gene_name"]).to_csv(
        REPORT / "pd_blood_apc_early_isp_summary_full6.csv", index=False)
    agg.round(5).to_csv(REPORT / "pd_blood_apc_early_isp_gene_ranking_full6.csv",
                        index=False)
    OUT.joinpath("isp_gene_ranking.csv").write_text(agg.round(5).to_csv(index=False))

    print("authoritative rows:", len(df),
          "non-null:", df["Shift_to_goal_end"].notna().sum())
    print("wrote summary + ranking + figures:",
          sorted(p.name for p in OUT.glob("isp_*")))


if __name__ == "__main__":
    main()