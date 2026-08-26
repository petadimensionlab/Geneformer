#!/usr/bin/env python
"""Build PD_spleen report data CSVs + visualizations.

Reads the early-PD perturbation results produced by
analysis/07_pd_spleen_early_isp.py
  input/PD_spleen/results/isp/pd_spleen/pd_spleen_early_isp_stats_combined.csv
(one row per gene per timepoint, columns: Timepoint, Gene, Ensembl_ID,
 Shift_to_goal_end, Shift_6m, Early_rank)

and writes, into result/report/:
  data/pd_spleen_gene_list.csv            — gene symbol, ENSG, category
  data/pd_spleen_early_isp_results.csv    — 1 row/gene: Shift at 6m/9m/12m
  figures/pd_spleen_isp_rank_6m.png       — 6m ranked bar (category colors)
  figures/pd_spleen_isp_timeseries.png    — top genes across 3 timepoints
  figures/pd_spleen_isp_early_specific.png — 6m vs later-timepoint scatter
  figures/pd_spleen_isp_heatmap.png       — gene x timepoint heatmap

Usage:
    .venv/bin/python result/report/visualize_isp_pd_spleen.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPORT = Path(__file__).resolve().parent
DATA = REPORT / "data"
FIG = REPORT / "figures"
REPORT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

SRC = Path(__file__).resolve().parent.parent.parent / "input" / "PD_spleen" / "results" / "isp" / "pd_spleen"
SRC_CSV = SRC / "pd_spleen_early_isp_stats_combined.csv"

# hypothesis gene categories (mirror analysis/07_pd_spleen_early_isp.py)
GENE_CATEGORIES = {
    # alpha-synuclein / PD core
    "SNCA": "PD_core", "LRRK2": "PD_core", "PARK7": "PD_core",
    "PINK1": "PD_core", "PRKN": "PD_core", "VPS35": "PD_core",
    "GBA1": "PD_core", "ATP13A2": "PD_core",
    # lysosome / autophagy
    "LAMP2": "lysosome_autophagy", "CTSB": "lysosome_autophagy",
    "CTSD": "lysosome_autophagy",
    # neuroinflammation / immune
    "TNF": "neuroinflammation", "IL1B": "neuroinflammation",
    "IL6": "neuroinflammation", "CXCL8": "neuroinflammation",
    "CX3CR1": "neuroinflammation", "TYROBP": "neuroinflammation",
    "TREM2": "neuroinflammation", "TREM1": "neuroinflammation",
    "S100A8": "neuroinflammation", "S100A9": "neuroinflammation",
    "LYZ": "neuroinflammation", "CD68": "neuroinflammation",
    "CSF1R": "neuroinflammation", "ITGAM": "neuroinflammation",
    "ITGAX": "neuroinflammation", "HLA-DRA": "neuroinflammation",
    "FOXP3": "neuroinflammation",
    # complement / chemokines
    "C1QA": "complement_chemokine", "C1QB": "complement_chemokine",
    "CXCL10": "complement_chemokine", "CCL2": "complement_chemokine",
    "CCL3": "complement_chemokine", "CD14": "complement_chemokine",
    "GAPDH": "housekeeping",
}

CAT_COLORS = {
    "PD_core": "#c0392b",
    "lysosome_autophagy": "#8e44ad",
    "neuroinflammation": "#e67e22",
    "complement_chemokine": "#2980b9",
    "housekeeping": "#7f8c8d",
}
DEFAULT_COLOR = "#7f8c8d"


def main() -> None:
    df = pd.read_csv(SRC_CSV)
    df["Shift_to_goal_end"] = df["Shift_to_goal_end"].astype(float)

    if "Ensembl_ID" not in df.columns:
        df["Ensembl_ID"] = ""
    # ---- gene list csv (dedupe by gene, keep category) ----
    gene_df = df[["Gene", "Ensembl_ID"]].drop_duplicates()
    gene_df["Category"] = gene_df["Gene"].map(GENE_CATEGORIES).fillna("other")
    gene_df = gene_df.sort_values(["Category", "Gene"])
    gene_df.to_csv(DATA / "pd_spleen_gene_list.csv", index=False)
    print("wrote", DATA / "pd_spleen_gene_list.csv")

    # ---- per-gene per-timepoint wide table ----
    wide = df.pivot_table(
        index="Gene", columns="Timepoint", values="Shift_to_goal_end"
    ).reindex(columns=["6m", "9m", "12m"])
    wide["Ensembl_ID"] = df.drop_duplicates("Gene").set_index("Gene")["Ensembl_ID"]
    wide["Category"] = wide.index.map(GENE_CATEGORIES).fillna("other")
    wide = wide[["Ensembl_ID", "Category", "6m", "9m", "12m"]]
    wide = wide.reset_index().rename(columns={"index": "Gene"})
    # early rank: sort by 6m desc
    wide = wide.sort_values("6m", ascending=False)
    wide.insert(0, "Early_rank", range(1, len(wide) + 1))
    wide.to_csv(DATA / "pd_spleen_early_isp_results.csv", index=False)
    print("wrote", DATA / "pd_spleen_early_isp_results.csv")

    _plot_rank_6m(wide)
    _plot_timeseries(wide)
    _plot_early_specificity(wide)
    _plot_heatmap(wide)

    print("figures written to", FIG)


def _colors(series: pd.Series) -> list[str]:
    return [CAT_COLORS.get(c, DEFAULT_COLOR) for c in series]


def _plot_rank_6m(wide: pd.DataFrame) -> None:
    d = wide.sort_values("6m", ascending=True)
    labels = d["Gene"].tolist()
    vals = d["6m"].astype(float).tolist()
    colors = _colors(d["Category"])
    fig, ax = plt.subplots(figsize=(9, max(5, 0.28 * len(vals))))
    ax.barh(range(len(vals)), vals, color=colors, edgecolor="none")
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="black", linewidth=0.9)
    ax.set_xlabel("Shift_to_goal_end at 6m (positive = deletion pushes early-PD → early-WT)")
    ax.set_title("Early-PD (PFF 6m) in-silico gene deletion — Spleen immune pool", fontsize=12)
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v:+.4f}", va="center", fontsize=7)
    present = sorted(set(d["Category"]))
    handles = [plt.Rectangle((0, 0), 1, 1, color=CAT_COLORS.get(c, DEFAULT_COLOR)) for c in present]
    if handles:
        ax.legend(handles, present, loc="lower right", fontsize=7)
    ax.set_xlim(left=min(vals) - 0.001, right=max(vals) + 0.001)
    plt.tight_layout()
    plt.savefig(FIG / "pd_spleen_isp_rank_6m.png", dpi=150)
    plt.close()
    print("wrote", FIG / "pd_spleen_isp_rank_6m.png")


def _plot_timeseries(wide: pd.DataFrame) -> None:
    top = wide.sort_values("6m", ascending=False).head(8)
    order = top["Gene"].tolist()
    x = ["6m", "9m", "12m"]
    fig, ax = plt.subplots(figsize=(10, 5))
    for g in order:
        row = wide[wide["Gene"] == g]
        vals = row[["6m", "9m", "12m"]].astype(float).values.flatten()
        ax.plot(x, vals, marker="o", label=g, linewidth=1.8)
    ax.axhline(0, color="gray", linewidth=0.8, linestyle=":")
    ax.set_xlabel("Timepoint")
    ax.set_ylabel("Shift_to_goal_end")
    ax.set_title("Top-8 genes by 6m shift — time-course (6m → 12m)", fontsize=12)
    ax.legend(fontsize=8, ncol=2, loc="upper right")
    ax.grid(axis="y", linestyle=":", alpha=0.4)
    plt.tight_layout()
    plt.savefig(FIG / "pd_spleen_isp_timeseries.png", dpi=150)
    plt.close()
    print("wrote", FIG / "pd_spleen_isp_timeseries.png")


def _plot_early_specificity(wide: pd.DataFrame) -> None:
    x = wide["6m"].astype(float)
    y = wide[["9m", "12m"]].astype(float).mean(axis=1)
    labels = wide["Gene"].tolist()
    colors = _colors(wide["Category"])
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(x, y, c=colors, s=50, edgecolor="black", linewidth=0.5)
    lim = max(np.abs(x).max(), np.abs(y).max()) * 1.2
    ax.plot([-lim, lim], [-lim, lim], ls="--", color="gray", lw=1)
    ax.axhline(0, color="lightgray", lw=0.8)
    ax.axvline(0, color="lightgray", lw=0.8)
    for g, xi, yi in zip(labels, x, y):
        if abs(xi) > 5e-5 or abs(yi) > 5e-5:
            ax.annotate(g, (xi, yi), textcoords="offset points", xytext=(5, 4), fontsize=8)
    ax.set_xlabel("Shift at 6m (early)")
    ax.set_ylabel("Shift at 9m/12m (mean)")
    ax.set_title("Early-specificity: genes above diagonal act strongest at 6m", fontsize=12)
    present = sorted(set(wide["Category"]))
    handles = [plt.Rectangle((0, 0), 1, 1, color=CAT_COLORS.get(c, DEFAULT_COLOR)) for c in present]
    ax.legend(handles, present, loc="lower right", fontsize=7)
    plt.tight_layout()
    plt.savefig(FIG / "pd_spleen_isp_early_specific.png", dpi=150)
    plt.close()
    print("wrote", FIG / "pd_spleen_isp_early_specific.png")


def _plot_heatmap(wide: pd.DataFrame) -> None:
    d = wide.sort_values("6m", ascending=False)
    genes = d["Gene"].tolist()
    mat = d[["6m", "9m", "12m"]].astype(float).values
    vmax = np.abs(mat).max()
    fig, ax = plt.subplots(figsize=(7, max(6, 0.3 * len(genes))))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)
    ax.set_xticks(range(3))
    ax.set_xticklabels(["6m", "9m", "12m"], fontsize=10)
    ax.set_title("Per-gene Shift_to_goal_end by timepoint (early-PD → WT)", fontsize=12)
    for i in range(len(genes)):
        for j in range(3):
            ax.text(j, i, f"{mat[i, j]:+.4f}", ha="center", va="center", fontsize=6.5)
    fig.colorbar(im, label="Shift_to_goal_end")
    plt.tight_layout()
    plt.savefig(FIG / "pd_spleen_isp_heatmap.png", dpi=150)
    plt.close()
    print("wrote", FIG / "pd_spleen_isp_heatmap.png")


if __name__ == "__main__":
    main()