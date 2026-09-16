#!/usr/bin/env python
"""Figures for the 104M vs 316M ISP comparison report — bilingual (ja/en).

TEMPLATE — a starting point (ひな形), not a finished analysis.
The numbers, literature verdicts and prose baked into this file are a worked
example from the AD_spleen 104M-vs-316M check (2026-09). Re-derive every number
and re-verify every PMID before reusing it on another tissue or model.

Reads docs/quantization/isp-model-comparison/{gene_table.csv,summary.json,adjudication.csv}
and docs/quantization/expression/AD_spleen.json, writes PNGs into
docs/quantization/isp-model-comparison/figures/<lang>/.

Usage: .venv/bin/python templates/template_report_figures.py [--lang ja|en|both]

Japanese uses the Noto Sans CJK JP font; the script fails loudly if it cannot be
resolved, and missing glyphs are turned into errors, so a figure never ships with
tofu boxes.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

D = pathlib.Path("docs/quantization/isp-model-comparison")
CANDIDATES = ["Noto Sans CJK JP", "Noto Sans JP", "IPAexGothic", "IPAGothic", "Droid Sans Fallback"]

T = {
    "ja": {
        "f1_title": "同じ細胞・同じ遺伝子・同じ時点で測った 104M と 316M の遺伝子シフト",
        "f1_p1": "全55遺伝子（外れ値を含む）",
        "f1_p2": "拡大（|shift| < 0.004 = 全体の大半）",
        "f1_x": "104M の Shift（3m、正 = WT状態へ）",
        "f1_y": "316M の Shift（3m）",
        "f1_l1": "符号が一致", "f1_l2": "符号の反転", "f1_l3": "完全一致（y = x）",
        "f1_note": "第2・第4象限 = 2モデルで符号が反転\n（どちらが正しいかはモデル自身からは決まらない）",
        "f2_labels": ["精度差\nfp32 対 bf16\n(同じ100セル)",
                      "細胞サンプリング差\n200 対 100セル\n(同じ316M)",
                      "モデル差\n104M 対 316M"],
        "f2_x": "平均 |差| ÷ 平均 |シフト|（%）= シグナルに対する差の大きさ",
        "f2_title": "差の正体を分離する：モデル差はサンプリング誤差の約5倍、精度差の約7倍",
        "f2_txt": "符号反転",
        "f3_title": "2モデルで向きが逆転し、かつ文献が方向を示している遺伝子（7遺伝子）",
        "f3_ylabel": "{g}  検出{det:.0f}% / 文献={lit}",
        "f3_lit_pos": "欠失で改善", "f3_lit_neg": "欠失で悪化",
        "f3_x": "Shift_to_goal_end（3m、正 = WT状態に近づく）",
        "f3_agree": "文献と一致: {m}",
        "f3_note": "LYZ は104Mで突出（1桁大きい）",
        "f4_p1": "シフトの大きさは発現の広さに強く依存する",
        "f4_p2": "316M でも同じ傾向（ただし上限が低い）",
        "f4_x": "標的細胞での検出率 (%)",
        "f4_y": "|Shift|（対数）",
        "f4_sup": "赤 = 2モデルで符号が反転した遺伝子",
        "f5_x": "AD の細胞数の偏り（%、0 = AD/WT 同数）",
        "f5_title": "プール内の細胞構成：AD と WT でほぼ同一（上流で細胞型ごとに最大3000セルに間引き）",
    },
    "en": {
        "f1_title": "Gene shifts of 104M and 316M measured on the same cells, genes and timepoint",
        "f1_p1": "All 55 genes (outliers included)",
        "f1_p2": "Zoom: |shift| < 0.004 (most of the panel)",
        "f1_x": "Shift in 104M (3m; positive = toward the WT state)",
        "f1_y": "Shift in 316M (3m)",
        "f1_l1": "sign agrees", "f1_l2": "sign reversed", "f1_l3": "perfect agreement (y = x)",
        "f1_note": "Quadrants 2 and 4 = the two models disagree on the sign\n(the models alone cannot say which is right)",
        "f2_labels": ["precision\nfp32 vs bf16\n(same 100 cells)",
                      "cell sampling\n200 vs 100 cells\n(same 316M model)",
                      "model\n104M vs 316M"],
        "f2_x": "mean |difference| ÷ mean |shift| (%) = size of the gap relative to the signal",
        "f2_title": "Separating the causes: the model gap is ~5x the sampling error and ~7x the precision error",
        "f2_txt": "sign flips",
        "f3_title": "Genes where the two models disagree in direction and the literature indicates one (7 genes)",
        "f3_ylabel": "{g}  det {det:.0f}% / lit={lit}",
        "f3_lit_pos": "deletion helps", "f3_lit_neg": "deletion harms",
        "f3_x": "Shift_to_goal_end (3m; positive = closer to the WT state)",
        "f3_agree": "literature agrees: {m}",
        "f3_note": "LYZ dominates in 104M (an order of magnitude larger)",
        "f4_p1": "Shift size depends strongly on expression breadth",
        "f4_p2": "Same trend in 316M (with a lower ceiling)",
        "f4_x": "Detection rate in the target cells (%)",
        "f4_y": "|Shift| (log scale)",
        "f4_sup": "Red = genes whose sign is reversed between the two models",
        "f5_x": "Bias in AD cell counts (%, 0 = equal AD/WT)",
        "f5_title": "Pool composition: nearly identical in AD and WT (upstream caps each cell type at 3,000 cells)",
    },
}


def setup_font() -> None:
    have = {f.name for f in fm.fontManager.ttflist}
    for c in CANDIDATES:
        if c in have:
            plt.rcParams["font.family"] = c
            plt.rcParams["axes.unicode_minus"] = False
            print(f"font: {c}")
            return
    raise SystemExit(f"no CJK font found among {CANDIDATES}; have={sorted(have)[:20]}")


def fig1_scatter(t: pd.DataFrame, L: dict, out: pathlib.Path) -> None:
    from matplotlib.lines import Line2D
    fig, axes = plt.subplots(1, 2, figsize=(13.5, 6.2))
    for ax, (lim, title) in zip(axes, [(0.115, L["f1_p1"]), (0.004, L["f1_p2"])]):
        for _, r in t.iterrows():
            ax.scatter(r.shift_104M, r.shift_316M, s=48,
                       c="#d62728" if r.sign_flip else "#1f77b4", alpha=.85,
                       edgecolor="white", linewidth=.6, zorder=3)
        ax.plot([-lim, lim], [-lim, lim], ls="--", lw=1, c="#999", zorder=1)
        ax.axhline(0, lw=.8, c="#bbb", zorder=2); ax.axvline(0, lw=.8, c="#bbb", zorder=2)
        ax.set_xlim(-lim * 1.05, lim * 1.05); ax.set_ylim(-lim * 1.05, lim * 1.05)
        ax.set_xlabel(L["f1_x"], fontsize=10); ax.set_ylabel(L["f1_y"], fontsize=10)
        ax.set_title(title, fontsize=11); ax.grid(alpha=.22, zorder=0)
    for g, off in {"CD8A": (-18, 6), "LYZ": (10, -2)}.items():
        r = t.loc[g]
        axes[0].annotate(g, (r.shift_104M, r.shift_316M), fontsize=9, fontweight="bold",
                         xytext=off, textcoords="offset points")
    for g, off in {"APOE": (-34, -12), "TYROBP": (12, 4), "GCA": (-40, -13), "LYZ": (14, 2),
                   "XCR1": (-30, 8), "S100A8": (10, 8), "ABCA7": (-36, 6), "PLCG2": (10, -12),
                   "CD34": (12, -14), "CLEC9A": (-40, -14)}.items():
        r = t.loc[g]
        if abs(r.shift_104M) < 0.004 and abs(r.shift_316M) < 0.004:
            axes[1].annotate(g, (r.shift_104M, r.shift_316M), fontsize=9,
                             xytext=off, textcoords="offset points")
    axes[1].text(0.012, -0.0032, L["f1_note"], fontsize=9.5, c="#d62728", ha="right")
    handles = [Line2D([], [], marker="o", ls="", color="#1f77b4", label=L["f1_l1"]),
               Line2D([], [], marker="o", ls="", color="#d62728", label=L["f1_l2"]),
               Line2D([], [], ls="--", color="#999", label=L["f1_l3"])]
    axes[0].legend(handles=handles, loc="lower right", fontsize=9, framealpha=.95)
    fig.suptitle(L["f1_title"], fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, .96))
    fig.savefig(out / "fig1_scatter.png", dpi=150); plt.close(fig)


def fig2_noise(s: dict, L: dict, out: pathlib.Path) -> None:
    rows = s["noise_and_contrasts"]
    pct = [r["delta_as_pct_of_signal"] for r in rows]
    flips = [f"{r['sign_flips']}/{r['n_genes']}" for r in rows]
    rho = [r["spearman"] for r in rows]
    fig, ax = plt.subplots(figsize=(9.6, 4.6))
    cols = ["#7fbf7b", "#f0ad4e", "#d62728"]
    bars = ax.barh(L["f2_labels"][::-1], pct[::-1], color=cols[::-1], height=.55)
    for b, p, f, rr in zip(bars, pct[::-1], flips[::-1], rho[::-1]):
        ax.text(b.get_width() + 1.5, b.get_y() + b.get_height() / 2,
                f"{p:.1f}%   {L['f2_txt']} {f}   ρ={rr:+.3f}", va="center", fontsize=10)
    ax.set_xlim(0, 128)
    ax.set_xlabel(L["f2_x"])
    ax.set_title(L["f2_title"], fontsize=12)
    ax.grid(axis="x", alpha=.25)
    fig.tight_layout(); fig.savefig(out / "fig2_noise_ladder.png", dpi=150); plt.close(fig)


def fig3_decisive(adj: pd.DataFrame, L: dict, out: pathlib.Path) -> None:
    d = adj[(adj.sign_flip) & (adj.lit_verdict.isin(["POS", "NEG"])) &
            (adj.detection_AD >= 0.005)].copy().sort_values("shift_104M")
    small, big = d[d.gene != "LYZ"], d[d.gene == "LYZ"]
    fig = plt.figure(figsize=(12.8, 6.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[len(small), 1.15], hspace=.45,
                          left=.22, right=.78, top=.88, bottom=.09)

    def panel(ax, sub, lim, note):
        y = np.arange(len(sub))
        ax.barh(y + .19, sub.shift_104M, height=.36, color="#1f77b4", label="104M")
        ax.barh(y - .19, sub.shift_316M, height=.36, color="#ff7f0e", label="316M")
        ax.set_yticks(y)
        ax.set_yticklabels([L["f3_ylabel"].format(
            g=g, det=det * 100,
            lit=L["f3_lit_pos"] if v == "POS" else L["f3_lit_neg"])
            for g, det, v in zip(sub.gene, sub.detection_AD, sub.lit_verdict)], fontsize=10)
        ax.axvline(0, c="#444", lw=1.1); ax.set_xlim(-lim, lim)
        for i, (a, b) in enumerate(zip(sub.shift_104M, sub.shift_316M)):
            for v, off, col in ((a, .19, "#1f77b4"), (b, -.19, "#ff7f0e")):
                ax.text(v + (lim * .03 if v >= 0 else -lim * .03), i + off,
                        f"{v:+.4f}" if lim > .01 else f"{v:+.5f}", va="center",
                        ha="left" if v >= 0 else "right", fontsize=8, color=col)
        for i, m in enumerate(sub.better_supported_model):
            ax.text(1.02, i, L["f3_agree"].format(m=m), transform=ax.get_yaxis_transform(),
                    va="center", fontsize=9.5, fontweight="bold",
                    color="#2ca02c" if m.startswith("316M") else "#d62728")
        ax.grid(axis="x", alpha=.25)
        if note:
            ax.text(.985, .06, note, transform=ax.transAxes, ha="right", fontsize=8.5, color="#555")

    panel(fig.add_subplot(gs[0]), small, 0.0015, "")
    panel(fig.add_subplot(gs[1]), big, 0.028, L["f3_note"])
    # keep the legend clear of the per-row "literature agrees" column
    fig.axes[0].legend(loc="lower left", bbox_to_anchor=(1.005, 1.06), fontsize=10)
    fig.axes[0].set_xlabel(L["f3_x"], fontsize=10)
    fig.suptitle(L["f3_title"], fontsize=13, y=.975)
    fig.savefig(out / "fig3_decisive_flips.png", dpi=150); plt.close(fig)


def fig4_detection(t: pd.DataFrame, L: dict, out: pathlib.Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, col, title in [(axes[0], "shift_104M", L["f4_p1"]), (axes[1], "shift_316M", L["f4_p2"])]:
        ax.scatter(t.detection_AD * 100, t[col].abs(), s=45,
                   c=["#d62728" if f else "#1f77b4" for f in t.sign_flip], alpha=.85,
                   edgecolor="white", linewidth=.6)
        ax.set_yscale("log")
        ax.set_xlabel(L["f4_x"]); ax.set_ylabel(f"{L['f4_y']}")
        ax.set_title(title, fontsize=11); ax.grid(alpha=.25, which="both")
    for g in ["CD8A", "LYZ", "S100A8", "MS4A1", "CLEC9A", "TREM2", "CCL2", "SPP1"]:
        r = t.loc[g]
        axes[0].annotate(g, (r.detection_AD * 100, abs(r.shift_104M)), fontsize=8,
                         xytext=(4, 3), textcoords="offset points")
    fig.suptitle(L["f4_sup"], fontsize=12)
    fig.tight_layout(); fig.savefig(out / "fig4_detection_vs_shift.png", dpi=150); plt.close(fig)


def fig5_cells(comp: dict, L: dict, out: pathlib.Path) -> None:
    df = pd.DataFrame(comp).T
    df["AD"] = df["AD"].astype(float); df["WT"] = df["WT"].astype(float)
    df["frac"] = df["AD"] / (df["AD"] + df["WT"])
    df = df.sort_values("frac")
    fig, ax = plt.subplots(figsize=(10, 7))
    y = np.arange(len(df))
    ax.barh(y, (df.frac - .5) * 200,
            color=["#d62728" if abs(f - .5) > .03 else "#1f77b4" for f in df.frac])
    ax.axvline(0, c="#444", lw=1)
    ax.set_yticks(y); ax.set_yticklabels(df.index, fontsize=9)
    ax.set_xlabel(L["f5_x"]); ax.set_title(L["f5_title"], fontsize=11)
    ax.grid(axis="x", alpha=.25)
    for i, (n_ad, n_wt) in enumerate(zip(df.AD, df.WT)):
        ax.text((df.frac.iloc[i] - .5) * 200 + (1.2 if df.frac.iloc[i] >= .5 else -1.2), i,
                f"{int(n_ad)}/{int(n_wt)}", va="center",
                ha="left" if df.frac.iloc[i] >= .5 else "right", fontsize=7.5, color="#333")
    fig.tight_layout(); fig.savefig(out / "fig5_cell_composition.png", dpi=150); plt.close(fig)


def build(lang: str) -> None:
    L = T[lang]
    out = D / "figures" / lang
    out.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(D / "gene_table.csv", index_col=0)
    s = json.loads((D / "summary.json").read_text())
    adj = pd.read_csv(D / "adjudication.csv")
    expr = json.loads(pathlib.Path("docs/quantization/expression/AD_spleen.json").read_text())
    fig1_scatter(t, L, out); fig2_noise(s, L, out); fig3_decisive(adj, L, out)
    fig4_detection(t, L, out); fig5_cells(expr["cell_composition"], L, out)
    print(f"[{lang}] -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--lang", default="both", choices=["ja", "en", "both"])
    a = ap.parse_args()
    setup_font()
    warnings.simplefilter("error", UserWarning)      # missing glyphs and layout warnings surface
    for lang in (["ja", "en"] if a.lang == "both" else [a.lang]):
        build(lang)
    for p in sorted((D / "figures").glob("*/*.png")):
        print(f"  {p}  {p.stat().st_size/1024:.0f} KiB")


if __name__ == "__main__":
    main()
