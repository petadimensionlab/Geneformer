#!/usr/bin/env python
"""Visualize in-silico-perturbation results for the early-AD detection reports.

Reads the aggregated per-experiment results CSVs and writes horizontal bar
charts (Shift_to_goal_end, descending) colored by gene category, plus a
small-intestine vs brain comparison. Outputs go to result/report/figures/.

Usage:
    .venv/bin/python result/report/visualize_isp.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPORT = Path(__file__).resolve().parent
DATA = REPORT / "data"
FIG = REPORT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

CAT_COLORS = {
    "AD_risk": "#c0392b",
    "inflammation_immune": "#e67e22",
    "DAM_microglia": "#e67e22",
    "inflammation": "#e67e22",
    "gut_barrier_mucin": "#27ae60",
    "gut_barrier_tight_junction": "#16a085",
    "gut_antimicrobial": "#2980b9",
    "gut_stem_metaplasia": "#8e44ad",
    "other": "#7f8c8d",
}
DEFAULT_COLOR = "#7f8c8d"


def load_results(csv_name: str, gene_list_csv: str) -> pd.DataFrame:
    df = pd.read_csv(DATA / csv_name)
    df["Shift_to_goal_end"] = df["Shift_to_goal_end"].astype(float)
    glist = pd.read_csv(DATA / gene_list_csv)[["Gene", "Category"]]
    df = df.merge(glist, on="Gene", how="left")
    df["Category"] = df["Category"].fillna("other")
    return df.sort_values("Shift_to_goal_end", ascending=True)


def plot_rank(df: pd.DataFrame, title: str, out_path: Path, top_n: int | None = None) -> None:
    plot_df = df.tail(top_n) if top_n else df
    labels = [str(g) for g in plot_df["Gene"]]
    vals = plot_df["Shift_to_goal_end"].tolist()
    colors = [CAT_COLORS.get(c, DEFAULT_COLOR) for c in plot_df["Category"]]

    fig, ax = plt.subplots(figsize=(9, max(4, 0.28 * len(vals))))
    ax.barh(range(len(vals)), vals, color=colors, edgecolor="none")
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("Shift_to_goal_end  (positive = deletion pushes early-AD → early-WT)")
    ax.set_title(title, fontsize=12)
    ax.grid(axis="x", linestyle=":", alpha=0.5)
    for i, v in enumerate(vals):
        ax.text(v, i, f" {v:+.4f}", va="center", fontsize=7)
    present = sorted(set(plot_df["Category"]))
    handles = [plt.Rectangle((0, 0), 1, 1, color=CAT_COLORS.get(c, DEFAULT_COLOR)) for c in present]
    if handles:
        ax.legend(handles, present, loc="lower right", fontsize=7)
    ax.set_xlim(left=min(vals) - 0.001, right=max(vals) + 0.001)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print("wrote", out_path)


def main() -> None:
    # ---- small intestine (49 genes) ----
    sm = load_results("ad_smallint_early_isp_results.csv", "ad_smallint_gene_list.csv")
    plot_rank(
        sm,
        "Early-AD in-silico gene deletion — Small Intestine (49 genes)",
        FIG / "ad_smallint_isp_rank.png",
    )

    # ---- brain (23 genes) ----
    br = load_results("ad_brain_early_isp_results.csv", "ad_brain_gene_list.csv")
    plot_rank(
        br,
        "Early-AD in-silico gene deletion — Brain microglia/BAM (23)",
        FIG / "ad_brain_isp_rank.png",
    )

    # ---- combined comparison (top 12 each) ----
    sm_top = sm.sort_values("Shift_to_goal_end", ascending=False).head(12)
    br_top = br.sort_values("Shift_to_goal_end", ascending=False).head(12)
    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for ax, d, t in [(axes[0], sm_top, "Small intestine top-12"), (axes[1], br_top, "Brain top-12")]:
        vals = d["Shift_to_goal_end"].values
        labels = d["Gene"].tolist()
        colors = [CAT_COLORS.get(c, DEFAULT_COLOR) for c in d["Category"]]
        ax.barh(range(len(vals)), vals, color=colors)
        ax.set_yticks(range(len(vals)))
        ax.set_yticklabels(labels, fontsize=9)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_title(t)
        ax.set_xlabel("Shift_to_goal_end")
        ax.invert_yaxis()
        ax.grid(axis="x", linestyle=":", alpha=0.3)
        for i, v in enumerate(vals):
            ax.text(v, i, f" {v:+.4f}", va="center", fontsize=7)
    plt.tight_layout()
    plt.savefig(FIG / "ad_smallint_vs_brain_top.png", dpi=150)
    plt.close()
    print("wrote", FIG / "ad_smallint_vs_brain_top.png")


if __name__ == "__main__":
    main()