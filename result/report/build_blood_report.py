#!/usr/bin/env python
"""Build AD_blood report data CSVs + figures for the in-silico-perturbation
report (early AD detection in blood).

Reads:
  input/AD_blood/results/isp/ad_blood/ad_blood_early_isp_stats_combined.csv
Writes:
  result/report/data/ad_blood_early_isp_results.csv  (all gene x timepoint rows)
  result/report/data/ad_blood_gene_list.csv          (Gene, Ensembl_ID, Category)
  result/report/figures/ad_blood_isp_rank_3m.png     (3m rank bar chart)
  result/report/figures/ad_blood_isp_timeseries.png  (per-gene Shift across time)
  result/report/figures/ad_blood_isp_heatmap.png     (gene x timepoint heatmap)
"""
from __future__ import annotations

import pickle
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]            # repo root
REPORT = Path(__file__).resolve().parent
DATA = REPORT / "data"
FIG = REPORT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

COMBINED = ROOT / "input" / "AD_blood" / "results" / "isp" / "ad_blood" / "ad_blood_early_isp_stats_combined.csv"

CATS = {
    # AD GWAS risk (immune-enriched)
    "APOE": "AD_risk", "TREM2": "AD_risk", "TYROBP": "AD_risk", "CLU": "AD_risk",
    "BIN1": "AD_risk", "CD33": "AD_risk", "ABCA7": "AD_risk", "SORL1": "AD_risk",
    "INPP5D": "AD_risk", "PLCG2": "AD_risk", "MEF2C": "AD_risk", "SPI1": "AD_risk",
    # complement / innate inflammation
    "C1QA": "inflammation_immune", "C1QB": "inflammation_immune", "C1QC": "inflammation_immune",
    "C3": "inflammation_immune", "CD68": "inflammation_immune", "CSF1R": "inflammation_immune",
    "AIF1": "inflammation_immune", "ITGAM": "inflammation_immune", "ITGAX": "inflammation_immune",
    "SPP1": "inflammation_immune", "CCL2": "inflammation_immune", "IL1B": "inflammation_immune",
    "TNF": "inflammation_immune", "TLR4": "inflammation_immune", "NLRP3": "inflammation_immune",
    "LYZ": "inflammation_immune", "APOC1": "inflammation_immune", "CD74": "inflammation_immune",
    "IRF7": "inflammation_immune", "ISG15": "inflammation_immune", "S100A8": "inflammation_immune",
    "S100A9": "inflammation_immune", "MS4A7": "inflammation_immune", "THBS1": "inflammation_immune",
    # adaptive (T / NK / B)
    "CD3E": "adaptive_immune", "CD3D": "adaptive_immune", "CD4": "adaptive_immune",
    "CD8A": "adaptive_immune", "CD19": "adaptive_immune", "MS4A1": "adaptive_immune",
    "NKG7": "adaptive_immune", "KLRD1": "adaptive_immune", "PRF1": "adaptive_immune",
    # myeloid / neutrophil / platelet marker
    "CSF3R": "myeloid_neutrophil", "GCA": "myeloid_neutrophil", "CEACAM1": "myeloid_neutrophil",
    "FGR": "myeloid_neutrophil", "ITGA2B": "myeloid_neutrophil", "VWF": "myeloid_neutrophil",
    "CD34": "myeloid_neutrophil", "XCR1": "myeloid_neutrophil", "IL3RA": "myeloid_neutrophil",
    "CLEC9A": "myeloid_neutrophil",
}

CAT_COLORS = {
    "AD_risk": "#c0392b",
    "inflammation_immune": "#e67e22",
    "adaptive_immune": "#16a085",
    "myeloid_neutrophil": "#2980b9",
}
DEFAULT_COLOR = "#7f8c8d"


def main() -> None:
    df = pd.read_csv(COMBINED)
    df["Category"] = df["Gene"].map(CATS).fillna("other")
    df.to_csv(DATA / "ad_blood_early_isp_results.csv", index=False)

    gene_list = (
        df[["Gene", "Ensembl_ID", "Category"]]
        .drop_duplicates("Gene")
        .sort_values("Category")
    )
    gene_list.to_csv(DATA / "ad_blood_gene_list.csv", index=False)
    print(f"wrote {len(df)} result rows ->", DATA / "ad_blood_early_isp_results.csv")
    print(f"wrote {len(gene_list)} genes ->", DATA / "ad_blood_gene_list.csv")

    # ---------------- figure 1: 3m rank ----------------
    tp3 = df[df["Timepoint"] == "3m"].sort_values("Shift_to_goal_end", ascending=False)
    colors = [CAT_COLORS.get(c, DEFAULT_COLOR) for c in tp3["Category"]]
    fig, ax = plt.subplots(figsize=(9, 13))
    ax.barh(tp3["Gene"], tp3["Shift_to_goal_end"], color=colors)
    ax.invert_yaxis()
    ax.set_xlabel("Shift_to_goal_end (deletion -> WT)")
    ax.set_title("AD_blood - in silico gene deletion, 3m (earliest) Shift_to_goal_end", fontsize=11)
    import matplotlib.patches as mpatches
    handles = [mpatches.Patch(color=col, label=lab) for lab, col in CAT_COLORS.items()]
    ax.legend(handles=handles, loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "ad_blood_isp_rank_3m.png", dpi=150)
    plt.close(fig)
    print("wrote", FIG / "ad_blood_isp_rank_3m.png")

    # ---------------- figure 2: timeseries ----------------
    pivot = df.pivot(index="Gene", columns="Timepoint", values="Shift_to_goal_end")
    pivot = pivot[["3m", "4p5m", "6m", "9m", "12m"]]
    top_genes = tp3.head(10)["Gene"].tolist()
    fig, ax = plt.subplots(figsize=(10, 6))
    for g in pivot.index:
        if g in top_genes:
            ax.plot(pivot.columns, pivot.loc[g], marker="o", lw=2, label=g)
        else:
            ax.plot(pivot.columns, pivot.loc[g], color="#d5dbe3", lw=0.8, alpha=0.6)
    ax.axhline(0, color="#999", lw=0.8, ls="--")
    ax.set_xlabel("timepoint (age at draw)")
    ax.set_ylabel("Shift_to_goal_end")
    ax.set_title("AD_blood — per-gene Shift_to_goal_end across timepoints\n(top10 at 3m highlighted)", fontsize=11)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "ad_blood_isp_timeseries.png", dpi=150)
    plt.close(fig)
    print("wrote", FIG / "ad_blood_isp_timeseries.png")

    # ---------------- figure 3: heatmap ----------------
    hm = pivot.copy()
    fig, ax = plt.subplots(figsize=(8, len(hm) * 0.32 + 2))
    im = ax.imshow(hm.values, cmap="RdBu_r", aspect="auto", vmin=-0.005, vmax=0.015)
    ax.set_xticks(range(len(hm.columns))); ax.set_xticklabels(hm.columns)
    ax.set_yticks(range(len(hm.index))); ax.set_yticklabels(hm.index, fontsize=8)
    ax.set_title("AD_blood - Shift_to_goal_end matrix (gene x timepoint)", fontsize=11)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label("Shift_to_goal_end")
    fig.tight_layout()
    fig.savefig(FIG / "ad_blood_isp_heatmap.png", dpi=150)
    plt.close(fig)
    print("wrote", FIG / "ad_blood_isp_heatmap.png")


if __name__ == "__main__":
    main()