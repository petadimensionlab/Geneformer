#!/usr/bin/env python3
"""Cross-pool concordance for the PD multi-region atlas ISP run (docs/isp/records/pd_atlas.md §6.3).

Question: do the genes that move most in one neuron pool move the same way in the
other? If the rankings disagree, a per-pool "top gene" is not a reproducible
signal and must not be reported as one.

Usage:
    python analysis/07h_pd_atlas_concordance.py [combined_stats.csv]

Default input (gitignored, produced by 07g_pd_atlas_perturbation.py):
    input/PD_atlas/results/isp/pd_atlas/pd_atlas_early_isp_stats_combined.csv
"""
import csv
import sys
from pathlib import Path
from statistics import mean, median

DEFAULT = ("input/PD_atlas/results/isp/pd_atlas/"
           "pd_atlas_early_isp_stats_combined.csv")


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(DEFAULT)
    if not path.exists():
        print(f"error: {path} not found -- run the ISP step first", file=sys.stderr)
        return 1

    pools: dict[str, dict[str, float]] = {}
    with path.open() as fh:
        for row in csv.DictReader(fh):
            pools.setdefault(row["Celltype"], {})[row["Gene"]] = float(
                row["Shift_to_goal_end"])

    print("=== pool summary ===")
    for ct, genes in pools.items():
        shifts = list(genes.values())
        print(f"{ct}: n={len(shifts)}  mean={mean(shifts):+.3e}  "
              f"median={median(shifts):+.3e}  max={max(shifts):+.3e}  "
              f"min={min(shifts):+.3e}  positive={sum(1 for s in shifts if s > 0)}"
              f"/{len(shifts)}")

    names = sorted(pools)
    if len(names) < 2:
        print("\nonly one pool present -- nothing to compare")
        return 0
    a, b = names[0], names[1]
    common = sorted(set(pools[a]) & set(pools[b]))

    print(f"\n=== {a} vs {b}: {len(common)} common genes ===")
    print(f"{'gene':10} {a:>12} {b:>12}  sign")
    agree = 0
    for g in common:
        sa, sb = pools[a][g], pools[b][g]
        same = (sa > 0) == (sb > 0)
        agree += same
        print(f"{g:10} {sa:+.6e} {sb:+.6e}  {'same' if same else 'OPPOSITE'}")

    xa = [pools[a][g] for g in common]
    xb = [pools[b][g] for g in common]

    from scipy.stats import pearsonr, spearmanr
    rho_s, p_s = spearmanr(xa, xb)
    rho_p, p_p = pearsonr(xa, xb)

    print(f"\nsign agreement : {agree}/{len(common)}")
    print(f"Spearman rho   : {rho_s:+.3f}  (p={p_s:.3f})")
    print(f"Pearson r      : {rho_p:+.3f}  (p={p_p:.3f})")

    rank_a = {g: i + 1 for i, g in enumerate(sorted(common, key=lambda g: -pools[a][g]))}
    rank_b = {g: i + 1 for i, g in enumerate(sorted(common, key=lambda g: -pools[b][g]))}
    print(f"\n{'gene':10} {a + ' rank':>12} {b + ' rank':>12}")
    for g in sorted(common, key=lambda g: rank_a[g]):
        print(f"{g:10} {rank_a[g]:>12} {rank_b[g]:>12}")

    print(f"\n{a} top gene: {min(common, key=lambda g: rank_a[g])} "
          f"({b} rank {rank_b[min(common, key=lambda g: rank_a[g])]})")
    print(f"{b} top gene: {min(common, key=lambda g: rank_b[g])} "
          f"({a} rank {rank_a[min(common, key=lambda g: rank_b[g])]})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
