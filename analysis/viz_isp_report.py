#!/usr/bin/env python3
"""Visualize the AD_brain in silico perturbation results.

Reads report/isp_results_consolidated.csv and writes figures to report/figures/:
  1. isp_shift_by_gene.png  — horizontal bar of Shift_to_goal_end per gene,
                              Microglia-only vs Microglia+BAM (two panels).
  2. isp_shift_comparison.png — scatter/paired comparison of the two cell pools.
  3. isp_shift_top_bottom.png  — ranked bar of the Microglia-only shifts with a
                              zero reference line.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent / "input" / "AD_brain"
REPORT = ROOT / "results" / "report"
FIG = REPORT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

# [gene, shift_micro_BAM, shift_micro_only] from isp_results_consolidated.csv
DATA = [
    ("APOE",  -0.000215621,  0.000713050),
    ("TREM2", -0.000215621, -0.000383679),
    ("TYROBP",-0.000236759, -0.000006272),
    ("BIN1",  -0.000236759, -0.000006272),
    ("CD33",  -0.000236759, -0.000006272),
    ("ABCA7", -0.000236759, -0.000006272),
    ("SORL1", -0.000236759, -0.000006272),
    ("INPP5D",-0.000236759, -0.000006272),
    ("PLCG2", -0.000236759, -0.000006272),
    ("MEF2C", -0.000236759, -0.000006272),
    ("CSF1R", -0.000236759, -0.000006272),
    ("C1QA",  -0.000236759, -0.000006272),
    ("C1QB",  -0.000236759, -0.000006272),
    ("C3",    -0.000236759, -0.000006272),
    ("CD68",  -0.000236759, -0.000006272),
    ("AIF1",  -0.000236759, -0.000006272),
    ("ITGAM", -0.000236759, -0.000006272),
    ("ITGAX", -0.000236759, -0.000006272),
    ("SPP1",  -0.000236759, -0.000006272),
    ("CCL2",  -0.000236759, -0.000006272),
    ("IL1B",  -0.000236759, -0.000006272),
    ("TNF",   -0.000236759, -0.000006272),
    ("CD74",  -0.000236759, -0.000006272),
]
genes = [d[0] for d in DATA]
shifts_mb = [d[1] for d in DATA]
shifts_mo = [d[2] for d in DATA]

# highlight APOE / TREM2
def color_for(shift, gene):
    if abs(shift) < 1e-4:
        return "#c9c9c9"
    return "#2b6cb0" if shift > 0 else "#c53030"

plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False})

# ---- 1. two-panel bar by gene ----
fig, axes = plt.subplots(1, 2, figsize=(13, 6), sharey=True)
for ax, data, title in [
    (axes[0], shifts_mb, "Microglia + BAM (E1)"),
    (axes[1], shifts_mo, "Microglia only (E2)"),
]:
    colors = [color_for(v, g) for g, v in zip(genes, data)]
    ax.barh(genes, data, color=colors, edgecolor="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_title(title, fontsize=11)
    ax.invert_yaxis()
axes[0].set_xlabel("Shift_to_goal_end (delete -> toward WT)")
axes[1].set_xlabel("Shift_to_goal_end (delete -> toward WT)")
fig.suptitle("In silico deletion shift toward early-WT in early-AD microglia\n"
             "(positive = deleting gene pushes AD cells toward healthy WT)",
             fontsize=12)
fig.tight_layout()
fig.savefig(FIG / "isp_shift_by_gene.png", dpi=150)
plt.close(fig)
print("wrote isp_shift_by_gene.png")

# ---- 2. paired scatter comparison ----
fig, ax = plt.subplots(figsize=(6, 6))
sc = ax.scatter(shifts_mb, shifts_mo, c=["#d33" if g == "TREM2"
                  else "#2b6cb0" if g == "APOE" else "#555" for g in genes])
lim = max(max(map(abs, shifts_mb)), max(map(abs, shifts_mo))) * 1.15
ax.plot([-lim, lim], [-lim, lim], ls="--", color="gray", lw=1)
for g, x, y in zip(genes, shifts_mb, shifts_mo):
    if abs(x) > 1e-4 or abs(y) > 1e-4:
        ax.annotate(g, (x, y), textcoords="offset points", xytext=(5, 4), fontsize=9)
ax.set_xlabel("Microglia + BAM shift")
ax.set_ylabel("Microglia only shift")
ax.set_title("Per-gene shift comparison (E1 vs E2)")
ax.axhline(0, color="lightgray", lw=0.8); ax.axvline(0, color="lightgray", lw=0.8)
fig.tight_layout()
fig.savefig(FIG / "isp_shift_comparison.png", dpi=150)
plt.close(fig)
print("wrote isp_shift_comparison.png")

# ---- 3. ranked bar (microglia only) ----
order = sorted(range(len(genes)), key=lambda i: shifts_mo[i], reverse=True)
rgenes = [genes[i] for i in order]
rshifts = [shifts_mo[i] for i in order]
fig, ax = plt.subplots(figsize=(8, 6))
colors = [color_for(v, g) for g, v in zip(rgenes, rshifts)]
ax.barh(rgenes, rshifts, color=colors, edgecolor="white", linewidth=0.5)
ax.axvline(0, color="black", linewidth=0.8)
ax.invert_yaxis()
ax.set_xlabel("Shift_to_goal_end (Microglia only)")
ax.set_title("Ranked per-gene shift — Microglia only (E2)")
fig.tight_layout()
fig.savefig(FIG / "isp_shift_top_bottom.png", dpi=150)
plt.close(fig)
print("wrote isp_shift_top_bottom.png")
print("figures written to", FIG)