#!/usr/bin/env python
"""Compare two in-silico-perturbation runs (accuracy gate G5).

Design in docs/quantization/PLAN.md (Phase 2, criterion G5): a quantized or
lower-precision run is only acceptable if it reproduces the *ranking* of
`Shift_to_goal_end`, so this reports rank agreement per timepoint — Spearman
rho, top-N overlap, sign agreement — not just a scalar loss.

Inputs are either an ISP experiment directory (the script finds the
`*_early_isp_stats_combined.csv` inside) or a CSV path directly.

Usage:
    .venv/bin/python analysis/14_compare_isp.py \
        --a input/AD_spleen/results/isp/ad_spleen_316m_fp32 \
        --b input/AD_spleen/results/isp/ad_spleen_316m_bf16 \
        --label 316m_fp32_vs_bf16 --out docs/quantization/g5

Environment: TOPN (default 20)
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

TOPN = int(__import__("os").environ.get("TOPN", "20"))


def resolve(path: str) -> tuple[pd.DataFrame, Path]:
    p = Path(path)
    if p.is_dir():
        cands = sorted(p.glob("*_early_isp_stats_combined.csv")) or sorted(p.glob("*combined*.csv"))
        if not cands:
            raise SystemExit(f"no combined CSV in {p}")
        p = cands[0]
    df = pd.read_csv(p)
    return df, p


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3:
        return float("nan")
    try:
        from scipy.stats import spearmanr

        return float(spearmanr(a, b).statistic)
    except Exception:
        ra = pd.Series(a).rank().to_numpy()
        rb = pd.Series(b).rank().to_numpy()
        return float(np.corrcoef(ra, rb)[0, 1])


def topn_overlap(df_a: pd.DataFrame, df_b: pd.DataFrame, n: int) -> tuple[int, list[str]]:
    ta = df_a.sort_values("Shift_to_goal_end", ascending=False)["Gene"].head(n).tolist()
    tb = df_b.sort_values("Shift_to_goal_end", ascending=False)["Gene"].head(n).tolist()
    shared = [g for g in ta if g in set(tb)]
    return len(shared), shared


def compare(da: pd.DataFrame, db: pd.DataFrame) -> dict:
    key = ["Timepoint", "Gene"]
    m = da.merge(db, on=key, suffixes=("_a", "_b"))
    if m.empty:
        return {"error": "no shared (Timepoint, Gene) rows"}
    x = m["Shift_to_goal_end_a"].to_numpy(dtype=float)
    y = m["Shift_to_goal_end_b"].to_numpy(dtype=float)
    diff = np.abs(x - y)
    out = {
        "n_rows": int(len(m)),
        "n_genes": int(m["Gene"].nunique()),
        "n_timepoints": int(m["Timepoint"].nunique()),
        "all": {
            "spearman": round(spearman(x, y), 5),
            "pearson": round(float(np.corrcoef(x, y)[0, 1]), 5),
            "sign_agreement": round(float(np.mean(np.sign(x) == np.sign(y))), 4),
            "topn_overlap": topn_overlap(m.assign(Shift_to_goal_end=x), m.assign(Shift_to_goal_end=y), TOPN)[0],
            "topn": TOPN,
            "mean_abs_diff": round(float(diff.mean()), 6),
            "max_abs_diff": round(float(diff.max()), 6),
            "mean_abs_shift": round(float(np.abs(x).mean()), 6),
        },
        "per_timepoint": {},
    }
    for tp, grp in m.groupby("Timepoint"):
        gx = grp["Shift_to_goal_end_a"].to_numpy(dtype=float)
        gy = grp["Shift_to_goal_end_b"].to_numpy(dtype=float)
        sub = grp.assign(Shift_to_goal_end=grp["Shift_to_goal_end_a"])
        ov, genes = topn_overlap(sub, grp.assign(Shift_to_goal_end=gy), TOPN)
        out["per_timepoint"][str(tp)] = {
            "n_genes": int(len(grp)),
            "spearman": round(spearman(gx, gy), 5),
            "sign_agreement": round(float(np.mean(np.sign(gx) == np.sign(gy))), 4),
            "topn_overlap": ov,
            "topn_shared": genes,
            "mean_abs_diff": round(float(np.abs(gx - gy).mean()), 6),
            "max_abs_diff": round(float(np.abs(gx - gy).max()), 6),
        }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, help="reference run (dir or CSV)")
    ap.add_argument("--b", required=True, help="candidate run (dir or CSV)")
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", default="docs/quantization/g5")
    args = ap.parse_args()

    da, pa = resolve(args.a)
    db, pb = resolve(args.b)
    res = {
        "a": {"path": str(pa), "rows": int(len(da))},
        "b": {"path": str(pb), "rows": int(len(db))},
        **compare(da, db),
    }
    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / f"{args.label}.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))

    print(f"A (reference): {pa}  rows={len(da)}")
    print(f"B (candidate): {pb}  rows={len(db)}")
    if "error" in res:
        print("!!", res["error"])
        return
    print(f"\nshared rows={res['n_rows']}  genes={res['n_genes']}  timepoints={res['n_timepoints']}")
    print(f"{'scope':10s} {'n':>4s} {'spearman':>9s} {'pearson':>8s} {'sign_agree':>11s} "
          f"{'top' + str(TOPN):>7s} {'mean|d|':>9s} {'max|d|':>9s}")
    s = res["all"]
    print(f"{'ALL':10s} {res['n_genes']:>4d} {s['spearman']:>9.5f} {s['pearson']:>8.5f} "
          f"{s['sign_agreement']:>11.4f} {s['topn_overlap']:>3d}/{TOPN:<3d} "
          f"{s['mean_abs_diff']:>9.6f} {s['max_abs_diff']:>9.6f}")
    for tp, s in res["per_timepoint"].items():
        print(f"{'tp ' + tp:10s} {s['n_genes']:>4d} {s['spearman']:>9.5f} {'-':>8s} "
              f"{s['sign_agreement']:>11.4f} {s['topn_overlap']:>3d}/{TOPN:<3d} "
              f"{s['mean_abs_diff']:>9.6f} {s['max_abs_diff']:>9.6f}")
    print(f"\nA mean |Shift| = {res['all']['mean_abs_shift']:.6f} "
          f"(scale reference for the absolute differences)")
    print(f"result -> {outdir / (args.label + '.json')}")


if __name__ == "__main__":
    main()
