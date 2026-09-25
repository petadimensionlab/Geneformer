#!/usr/bin/env python
"""Null comparison for the PD atlas ISP run: criteria.md E1, E2 and E3.

Reads the null run (pd_atlas_null) and the hypothesis run (pd_atlas) and reports:

  E1  signal threshold   the 95th percentile of the null distribution of |Shift|.
                         criteria.md 2-1: a gene counts only if it exceeds it.
                         Evidence kind A for this data and these conditions,
                         because the distribution is measured here.
  E2  rankability        effect size between the hypothesis genes and the null
                         (AUC = P(|Shift| of a hypothesis gene > that of a null
                         gene), Mann-Whitney U, mean difference).
  E3  negative controls  ACTB and GAPDH. They FAIL the 20% presence rule in these
                         neuron pools, so they cannot be perturbed here and the
                         check cannot be run -> 未検証, not pass or fail.

Reported with the two words that describe what we know (未検証 / 判定不能).
The measured values are printed beside the criterion; no pass/fail word is used
(style-guide.md A-3).

Usage:
    .venv/bin/python analysis/07j_pd_atlas_null_stats.py
"""
import csv
from pathlib import Path
from statistics import mean, median

BASE = Path(__file__).resolve().parent.parent / "input/PD_atlas/results/isp"
NULL_CSV = BASE / "pd_atlas_null" / "pd_atlas_null_early_isp_stats_combined.csv"
HYP_CSV = BASE / "pd_atlas" / "pd_atlas_early_isp_stats_combined.csv"
# ACTB and GAPDH measured presence in the PD cells of each pool (07i output).
# They are below the 20% rule, so they are absent from both runs by construction.
NEG_CONTROLS_PRESENCE = {"ACTB": 0.181, "GAPDH": 0.158}


def load(path: Path) -> dict[str, dict[str, tuple[float, float]]]:
    pools: dict[str, dict[str, tuple[float, float]]] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            pools.setdefault(row["Celltype"], {})[row["Gene"]] = (
                float(row["Shift_to_goal_end"]), float(row["Presence_frac"]))
    return pools


def iqr(xs: list[float]) -> float:
    s = sorted(xs)
    n = len(s)
    if n < 4:
        return float("nan")
    lo, hi = n // 4, (3 * n) // 4
    return median(s[hi:]) - median(s[:lo])


def pct95(xs: list[float]) -> float:
    """95th percentile, linear interpolation (numpy.percentile default)."""
    s = sorted(xs)
    if not s:
        return float("nan")
    k = 0.95 * (len(s) - 1)
    lo = int(k)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def main() -> int:
    if not NULL_CSV.exists():
        print(f"error: {NULL_CSV} not found -- the null run has not finished")
        return 1
    null = load(NULL_CSV)
    hyp = load(HYP_CSV)

    print("=" * 78)
    print("NULL COMPARISON -- criteria.md E1 / E2 / E3")
    print(f"null run      : {NULL_CSV.name}")
    print(f"hypothesis run: {HYP_CSV.name}")
    print("=" * 78)

    for pool in sorted(null):
        n_genes = null[pool]
        h_genes = hyp.get(pool, {})
        abs_null = [abs(v[0]) for v in n_genes.values()]
        thr = pct95(abs_null)

        print(f"\n##### {pool} #####")
        print(f"null set: n={len(abs_null)}")
        print(f"  |Shift| mean={mean(abs_null):.4e}  median={median(abs_null):.4e}  "
              f"IQR={iqr(abs_null):.4e}  min={min(abs_null):.4e}  max={max(abs_null):.4e}")
        print(f"  signed mean={mean(v[0] for v in n_genes.values()):+.4e}")
        print(f"\n  E1 threshold (95th percentile of |Shift|) = {thr:.4e}")
        print(f"  null genes: " + ", ".join(
            f"{g}({abs(v[0]):.2e})" for g, v in sorted(n_genes.items(), key=lambda kv: -abs(kv[1][0]))))

        if not h_genes:
            print(f"  no hypothesis results for {pool}")
            continue
        above = [g for g, v in h_genes.items() if abs(v[0]) > thr]
        print(f"\n  hypothesis genes above the threshold: {len(above)}/{len(h_genes)}")
        for g, v in sorted(h_genes.items(), key=lambda kv: -abs(kv[1][0])):
            mark = "ABOVE" if abs(v[0]) > thr else "below"
            print(f"    {g:10} {v[0]:+.4e}  |Shift|={abs(v[0]):.4e}  {mark}")

        # E2 effect size: AUC = P(|Shift| hyp > |Shift| null)
        a = [abs(v[0]) for v in h_genes.values()]
        wins = sum((x > y) + 0.5 * (x == y) for x in a for y in abs_null)
        auc = wins / (len(a) * len(abs_null))
        try:
            from scipy.stats import mannwhitneyu
            u, p = mannwhitneyu(a, abs_null, alternative="greater")
            p_txt = f"{p:.3g}"
        except Exception as exc:  # pragma: no cover
            p_txt = f"(scipy unavailable: {exc})"
        print(f"\n  E2 effect size: AUC={auc:.3f}  mean |Shift| hyp={mean(a):.4e} vs "
              f"null={mean(abs_null):.4e}  Mann-Whitney p={p_txt}")

    print("\n" + "=" * 78)
    print("E3 negative controls: " + ", ".join(
        f"{g} presence {p:.1%}" for g, p in NEG_CONTROLS_PRESENCE.items()))
    print("  both are below the 20% presence rule, so neither can be perturbed in")
    print("  these pools -> E3 は未検証 (not measurable). Report the observed")
    print("  presence fractions; do not write a pass or fail verdict.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
