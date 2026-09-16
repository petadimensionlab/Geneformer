#!/usr/bin/env python
"""Qualitative model comparison for the AD_spleen ISP: 104M vs 316M.

TEMPLATE — a starting point (ひな形), not a finished analysis.
The numbers, literature verdicts and prose baked into this file are a worked
example from the AD_spleen 104M-vs-316M check (2026-09). Re-derive every number
and re-verify every PMID before reusing it on another tissue or model.

Separates the three things that can make two ISP runs differ —
  (1) numerical precision (fp32 vs bf16)   -> measured on an identical-cell pair
  (2) cell sampling (200 vs 100 cells)     -> measured within the 316M model
  (3) the model itself (104M vs 316M)      -> the contrast we care about
— and then applies expression-based interpretability filters, because a gene that
is detected in <1% of the target cells cannot carry a meaningful shift.

Writes docs/quantization/isp-model-comparison/{gene_table.csv,summary.json}
and prints a report to stdout.

Usage: .venv/bin/python templates/template_isp_model_compare.py
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
ISP = ROOT / "input/AD_spleen/results/isp"
OUT = ROOT / "docs/quantization/isp-model-comparison"
EXPR = ROOT / "docs/quantization/expression/AD_spleen.json"

RUNS = {
    "m104_fp32_pool": ISP / "ad_spleen/ad_spleen_early_isp_stats_combined.csv",
    "m316_bf16_canary": ISP / "ad_spleen_316m_bf16/ad_spleen_316m_bf16_early_isp_stats_combined.csv",
    "m316_bf16_sub24": ISP / "ad_spleen_316m_bf16_sub24/ad_spleen_316m_bf16_sub24_early_isp_stats_combined.csv",
    "m316_fp32_sub24": ISP / "ad_spleen_316m_fp32_sub24/ad_spleen_316m_fp32_sub24_early_isp_stats_combined.csv",
}


def load(p: Path, shift_col: str = "Shift_3m") -> pd.Series:
    """gene -> shift at the 3m timepoint (the ranking timepoint)."""
    d = pd.read_csv(p)
    s = d.groupby("Gene")[shift_col].first()
    return s.astype(float)


def compare(a: pd.Series, b: pd.Series, label: str) -> dict:
    j = pd.DataFrame({"a": a, "b": b}).dropna()
    rho = spearmanr(j.a, j.b).statistic
    n = len(j)
    flip = int((np.sign(j.a) != np.sign(j.b)).sum())
    mean_abs = float((j.b - j.a).abs().mean())
    sig = float(j.a.abs().mean())
    top_n = max(5, min(20, n // 3))
    ta, tb = set(j.nlargest(top_n, "a").index), set(j.nlargest(top_n, "b").index)
    return {
        "label": label, "n_genes": n, "spearman": round(float(rho), 4),
        "sign_flips": flip, "mean_abs_delta": mean_abs,
        "mean_abs_signal": sig, "delta_as_pct_of_signal": round(100 * mean_abs / sig, 1),
        "top_overlap": f"{len(ta & tb)}/{top_n}",
        "top_only_a": sorted(ta - tb), "top_only_b": sorted(tb - ta),
    }


def main() -> None:
    runs = {k: load(p) for k, p in RUNS.items() if p.exists()}
    missing = [k for k in RUNS if k not in runs]
    if missing:
        print(f"[warn] missing runs: {missing}")

    expr = json.loads(EXPR.read_text())["genes"]
    det = {g: v["detection_AD"] for g, v in expr.items()}

    t = pd.DataFrame({
        "shift_104M": runs["m104_fp32_pool"],
        "shift_316M": runs["m316_bf16_canary"],
        "shift_316M_sub24": runs.get("m316_bf16_sub24"),
        "shift_316M_fp32_sub24": runs.get("m316_fp32_sub24"),
    })
    t["detection_AD"] = [det.get(g, np.nan) for g in t.index]
    t["rank_104M"] = t.shift_104M.rank(ascending=False).astype("Int64")
    t["rank_316M"] = t.shift_316M.rank(ascending=False).astype("Int64")
    t["rank_delta"] = t.rank_316M - t.rank_104M
    t["sign_104M"] = np.sign(t.shift_104M)
    t["sign_316M"] = np.sign(t.shift_316M)
    t["sign_flip"] = t.sign_104M != t.sign_316M
    t["candidate_104M"] = t.sign_104M > 0
    t["candidate_316M"] = t.sign_316M > 0
    t["candidate_flip"] = t.candidate_104M != t.candidate_316M
    t = t.sort_values("rank_104M")
    OUT.mkdir(parents=True, exist_ok=True)
    t.to_csv(OUT / "gene_table.csv", float_format="%.6g")

    summary = {"noise_and_contrasts": [], "expression_filters": []}
    # (1) precision
    if "m316_bf16_sub24" in runs and "m316_fp32_sub24" in runs:
        summary["noise_and_contrasts"].append(compare(
            runs["m316_fp32_sub24"], runs["m316_bf16_sub24"], "precision: 316M fp32 vs bf16 (same 100 cells)"))
    # (2) cell sampling, same model and dtype
    if "m316_bf16_sub24" in runs:
        summary["noise_and_contrasts"].append(compare(
            runs["m316_bf16_canary"], runs["m316_bf16_sub24"], "cell sampling: 316M bf16 200 vs 100 cells"))
    # (3) model contrast
    summary["noise_and_contrasts"].append(compare(
        runs["m104_fp32_pool"], runs["m316_bf16_canary"], "model: 104M fp32 vs 316M bf16 (200 cells each)"))

    for thr in (0.0, 0.01, 0.03, 0.05, 0.10):
        sub = t[t.detection_AD >= thr]
        if len(sub) < 5:
            continue
        rho = spearmanr(sub.shift_104M, sub.shift_316M).statistic
        summary["expression_filters"].append({
            "min_detection": thr, "n_genes": int(len(sub)),
            "spearman": round(float(rho), 4),
            "sign_flips": int(sub.sign_flip.sum()),
            "candidate_flips": int(sub.candidate_flip.sum()),
            "flip_pos_to_neg": sorted(sub.index[sub.candidate_104M & ~sub.candidate_316M]),
            "flip_neg_to_pos": sorted(sub.index[~sub.candidate_104M & sub.candidate_316M]),
        })

    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))

    print("=== contrasts (ranking timepoint 3m) ===")
    for c in summary["noise_and_contrasts"]:
        print(f"{c['label']:48s} n={c['n_genes']:3d} rho={c['spearman']:+.4f} "
              f"signflips={c['sign_flips']:2d} mean|d|={c['mean_abs_delta']:.3e} "
              f"({c['delta_as_pct_of_signal']:5.1f}% of signal) top={c['top_overlap']}")
    print("\n=== interpretability filter (104M vs 316M) ===")
    for f in summary["expression_filters"]:
        print(f"detection >= {f['min_detection']*100:4.0f}%: n={f['n_genes']:2d} rho={f['spearman']:+.3f} "
              f"sign flips={f['sign_flips']:2d} candidate flips={f['candidate_flips']:2d}")
    print("\n=== candidate-list flips among genes detected in >=3% of target cells ===")
    f = [x for x in summary["expression_filters"] if x["min_detection"] == 0.03][0]
    print("104M candidate -> 316M not:", ", ".join(f["flip_pos_to_neg"]))
    print("316M candidate -> 104M not:", ", ".join(f["flip_neg_to_pos"]))
    print(f"\n-> {OUT}/gene_table.csv, {OUT}/summary.json")


if __name__ == "__main__":
    main()
