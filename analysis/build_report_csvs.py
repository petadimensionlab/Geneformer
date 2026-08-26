#!/usr/bin/env python3
"""Build consolidated CSV reports for the AD_brain in silico perturbation runs.

Aggregates the per-gene goal_state_shift (Shift_to_goal_end) from the two
experiments into per-experiment tables + a consolidated comparison + a manifest.

Outputs to input/AD_brain/results/report/
"""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent / "input" / "AD_brain"
REPORT = ROOT / "results" / "report"
REPORT.mkdir(parents=True, exist_ok=True)

# (Ensembl, symbol, category) — the 23-gene hypothesis list used in 07_ad_spleen_early_isp
GENES = [
    ("ENSG00000130203", "APOE",   "AD risk / lipid"),
    ("ENSG00000095970", "TREM2",  "AD risk / DAM"),
    ("ENSG00000011600", "TYROBP", "DAM / adapter"),
    ("ENSG00000136717", "BIN1",   "AD risk"),
    ("ENSG00000105383", "CD33",   "AD risk / microglia"),
    ("ENSG00000064687", "ABCA7",  "AD risk / lipid"),
    ("ENSG00000137642", "SORL1",  "AD risk / endocytosis"),
    ("ENSG00000168918", "INPP5D", "AD risk / microglia"),
    ("ENSG00000197943", "PLCG2",  "AD risk / microglia"),
    ("ENSG00000081189", "MEF2C",  "AD risk / TF"),
    ("ENSG00000182578", "CSF1R",  "microglia survival"),
    ("ENSG00000173372", "C1QA",   "complement / DAM"),
    ("ENSG00000173369", "C1QB",   "complement / DAM"),
    ("ENSG00000125730", "C3",     "complement"),
    ("ENSG00000129226", "CD68",   "phagocytosis marker"),
    ("ENSG00000204472", "AIF1",   "microglia marker (IBA1)"),
    ("ENSG00000169896", "ITGAM",  "myeloid (CD11b)"),
    ("ENSG00000140678", "ITGAX",  "myeloid (CD11c)"),
    ("ENSG00000118785", "SPP1",   "DAM / neuroinflamm"),
    ("ENSG00000108691", "CCL2",   "chemokine"),
    ("ENSG00000125538", "IL1B",   "cytokine"),
    ("ENSG00000232810", "TNF",    "cytokine"),
    ("ENSG00000019582", "CD74",   "antigen / microglia"),
]

# experiments: (exp_id, suffix, cell_pool, n_cells, run_minutes)
EXPERIMENTS = [
    ("E1", "early_ad_brain",     "Microglia + BAM", 1000, "~72"),
    ("E2", "early_ad_microglia", "Microglia only",  1000, "~80"),
]


def load_shift(isp_dir: Path, ensg: str) -> float | None:
    csv_path = isp_dir / f"isp_stats_{ensg}.csv"
    if not csv_path.exists():
        return None
    with open(csv_path) as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    if len(lines) < 2:
        return None
    return float(lines[1].split(",")[-1])


def fmt(x: float | None) -> str:
    return "NA" if x is None else f"{x:.6g}"


def main() -> None:
    # --- per-experiment gene tables ---
    shift_by_exp = {suffix: {} for _, suffix, *_ in EXPERIMENTS}
    for exp_id, suffix, pool, n_cells, _ in EXPERIMENTS:
        isp_dir = ROOT / "results" / "isp" / suffix
        out = REPORT / f"isp_results_{suffix}.csv"
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["experiment", "Ensembl", "gene_symbol", "category",
                        "cell_pool", "n_cells", "Shift_to_goal_end"])
            for ensg, gene, cat in GENES:
                s = load_shift(isp_dir, ensg)
                shift_by_exp[suffix][ensg] = s
                w.writerow([exp_id, ensg, gene, cat, pool, n_cells, fmt(s)])
        print(f"wrote {out}")

    # --- consolidated comparison ---
    out_comp = REPORT / "isp_results_consolidated.csv"
    with open(out_comp, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["Ensembl", "gene_symbol", "category",
                    "Shift_Microglia_BAM", "Shift_Microglia_only"])
        for ensg, gene, cat in GENES:
            a = shift_by_exp["early_ad_brain"].get(ensg)
            b = shift_by_exp["early_ad_microglia"].get(ensg)
            w.writerow([ensg, gene, cat, fmt(a), fmt(b)])
    print(f"wrote {out_comp}")

    # --- manifest ---
    manifest = [
        {
            "experiment_id": "exp0",
            "script": "analysis/07_ad_brain_early_isp.py",
            "disease": "AD", "organ": "brain",
            "cell_pool": "Microglia + BAM",
            "timepoints": "early (AD3m/4.5m/6m vs WT3m/4.5m/6m)",
            "state_key": "disease", "start_state": "AD", "goal_state": "WT",
            "perturb_type": "delete", "combos": 0, "emb_mode": "cls",
            "model_version": "V2",
            "model": "fine-tuned CellClassifier (24 celltypes)",
            "n_cells": 200, "n_genes": 23,
            "output_dir": "results/isp/early_ad_brain (overwritten by exp1)",
            "results_csv": "not retained (re-run at n=1000)",
            "run_minutes": "~30",
            "notes": "APOE showed +0.0087 but this was sampling noise (not reproduced at n=1000)",
        },
        {
            "experiment_id": "exp1",
            "script": "analysis/07_ad_brain_early_isp.py",
            "disease": "AD", "organ": "brain",
            "cell_pool": "Microglia + BAM",
            "timepoints": "early (AD3m/4.5m/6m vs WT3m/4.5m/6m)",
            "state_key": "disease", "start_state": "AD", "goal_state": "WT",
            "perturb_type": "delete", "combos": 0, "emb_mode": "cls",
            "model_version": "V2",
            "model": "fine-tuned CellClassifier (24 celltypes)",
            "n_cells": 1000, "n_genes": 23,
            "output_dir": "results/isp/early_ad_brain",
            "results_csv": "report/isp_results_early_ad_brain.csv",
            "run_minutes": "~72",
            "notes": "APOE +0.0087 at n=200 but vanishes at n=1000 -> sampling noise; all genes ~0",
        },
        {
            "experiment_id": "exp2",
            "script": "analysis/07_ad_brain_early_isp.py",
            "disease": "AD", "organ": "brain",
            "cell_pool": "Microglia only",
            "timepoints": "early (AD3m/4.5m/6m vs WT3m/4.5m/6m)",
            "state_key": "disease", "start_state": "AD", "goal_state": "WT",
            "perturb_type": "delete", "combos": 0, "emb_mode": "cls",
            "model_version": "V2",
            "model": "fine-tuned CellClassifier (24 celltypes)",
            "n_cells": 1000, "n_genes": 23,
            "output_dir": "results/isp/early_ad_microglia",
            "results_csv": "report/isp_results_early_ad_microglia.csv",
            "run_minutes": "~80",
            "notes": "APOE only positive shift (+0.0007); TREM2 -0.0004; rest ~0",
        },
    ]
    out_m = REPORT / "experiment_manifest.csv"
    with open(out_m, "w", newline="") as f:
        cols = list(manifest[0].keys())
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(manifest)
    print(f"wrote {out_m}")


if __name__ == "__main__":
    main()