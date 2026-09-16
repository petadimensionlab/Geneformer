#!/usr/bin/env python
"""Observed expression direction of the ISP target genes, per cell type and state.

TEMPLATE — a starting point (ひな形), not a finished analysis.
The numbers, literature verdicts and prose baked into this file are a worked
example from the AD_spleen 104M-vs-316M check (2026-09). Re-derive every number
and re-verify every PMID before reusing it on another tissue or model.

The ISP question is "does deleting gene X move AD cells toward the WT state?", so
the independent, data-driven expectation for the sign is "is X up or down in AD
cells compared with WT?". This script extracts that directly from the tokenized
dataset (Geneformer tokens here are gene-only, so detection is exact).

Outputs a JSON + CSV with, for each target gene and each cell type:
  - detection fraction in AD cells and in WT cells
  - the AD-WT difference and log2 ratio
  - detection fraction per samples4 timepoint

Also reports the cell composition (AD vs WT) of the queried pool, which is the
"cell list" side of the same question.

Usage:
    .venv/bin/python templates/template_expression_by_celltype.py
Environment:
    ADPD_TISSUE (default AD_spleen), IS_CELLTYPES (comma list; default = the
    spleen immune pool used by 07_ad_spleen_early_isp.py), EXPR_OUT (default
    docs/quantization/expression/<TISSUE>.json)
"""
from __future__ import annotations

import json
import os
import pickle
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from datasets import load_from_disk

ROOT = Path(__file__).resolve().parent.parent
GF = Path(os.environ.get("GENEFORMER_DIR", ROOT / "geneformer_hf")) / "geneformer"
TISSUE = os.environ.get("ADPD_TISSUE", "AD_spleen")
PREFIX = TISSUE
DEFAULT_POOL = [
    "Ly6c.high.classical.Monocytes", "Ly6c.low.nonclassical.Monocytes", "Macrophages",
    "Neutrophils", "DCs", "cDC.1", "cDC.2", "pDCs", "Migratory.DCs", "CD4.T.cells",
    "CD8.T.cells", "Immature.T.cells", "Gamma.Delta.T.cells", "NKT.cells",
    "NK_ILC1", "Naive.Memory.B.cells", "Marzinal.zone.B.cells", "Plasma.cells",
    "Germinal.Center.B.cells",
]
POOL = [c.strip() for c in os.environ.get("IS_CELLTYPES", ",".join(DEFAULT_POOL)).split(",") if c.strip()]
HYPOTHESIS_GENES = [
    "APOE", "TREM2", "TYROBP", "CLU", "BIN1", "CD33", "ABCA7", "SORL1",
    "INPP5D", "PLCG2", "MEF2C", "SPI1",
    "C1QA", "C1QB", "C1QC", "C3", "CD68", "CSF1R", "AIF1", "ITGAM", "ITGAX",
    "SPP1", "CCL2", "IL1B", "TNF", "TLR4", "NLRP3", "LYZ", "APOC1", "CD74",
    "IRF7", "ISG15", "S100A8", "S100A9", "MS4A7", "THBS1",
    "CD3E", "CD3D", "CD4", "CD8A", "CD19", "MS4A1", "NKG7", "KLRD1", "PRF1",
    "CSF3R", "GCA", "CEACAM1", "FGR", "ITGA2B", "VWF", "CD34", "XCR1",
    "IL3RA", "CLEC9A",
]


def main() -> None:
    name2id = pickle.load(open(GF / "token_dictionary_gc104M.pkl", "rb"))
    name_id = pickle.load(open(GF / "gene_name_id_dict_gc104M.pkl", "rb"))
    # gene symbol -> ENSG -> token id
    sym2id = {}
    missing = []
    for sym in HYPOTHESIS_GENES:
        ensg = name_id.get(sym)
        if ensg is None or ensg not in name2id:
            missing.append(sym)
            continue
        sym2id[sym] = int(name2id[ensg])
    if missing:
        print(f"[warn] not found in dictionary: {missing}")

    ds = load_from_disk(str(ROOT / "input" / TISSUE / "tokenized" / f"{PREFIX}.dataset"))
    pool = set(POOL)
    target_ids = set(sym2id.values())

    # counters
    cells = Counter()                     # (celltype, disease) -> n
    tpcells = Counter()                   # (celltype, disease, samples4) -> n
    hits = defaultdict(Counter)           # gene -> (celltype, disease) -> n cells positive
    tp_hits = defaultdict(Counter)        # gene -> (disease, samples4) -> n
    tp_cells = Counter()                  # (disease, samples4) -> n

    for ex in ds:
        ct = ex["celltype"]
        if ct not in pool:
            continue
        dis, s4 = ex["disease"], ex["samples4"]
        cells[(ct, dis)] += 1
        tpcells[(ct, dis, s4)] += 1
        tp_cells[(dis, s4)] += 1
        present = target_ids.intersection(ex["input_ids"])
        for t in present:
            hits[t][(ct, dis)] += 1
    id2sym = {v: k for k, v in sym2id.items()}

    out = {"tissue": TISSUE, "pool": POOL, "n_cells_pool": int(sum(cells.values())),
           "cell_composition": {}, "genes": {}}
    for (ct, dis), n in sorted(cells.items()):
        out["cell_composition"].setdefault(ct, {})[dis] = int(n)
    for ct, d in out["cell_composition"].items():
        a, w = d.get("AD", 0), d.get("WT", 0)
        d["AD_minus_WT"] = a - w
        d["AD_fraction"] = round(a / (a + w), 4) if (a + w) else None

    for tid, sym in id2sym.items():
        per_cell = {}
        for (ct, dis), n in cells.items():
            if ct not in out["cell_composition"]:
                continue
            per_cell.setdefault(ct, {})[dis] = {
                "n": int(n),
                "detected": int(hits[tid][(ct, dis)]),
                "frac": round(hits[tid][(ct, dis)] / n, 4) if n else None,
            }
        ad_tot = sum(v["detected"] for ct, v in per_cell.items() for d, v in v.items() if d == "AD")
        ad_n = sum(v["n"] for ct, v in per_cell.items() for d, v in v.items() if d == "AD")
        wt_tot = sum(v["detected"] for ct, v in per_cell.items() for d, v in v.items() if d == "WT")
        wt_n = sum(v["n"] for ct, v in per_cell.items() for d, v in v.items() if d == "WT")
        f_ad = ad_tot / ad_n if ad_n else None
        f_wt = wt_tot / wt_n if wt_n else None
        out["genes"][sym] = {
            "token_id": int(tid),
            "detection_AD": round(f_ad, 4) if f_ad is not None else None,
            "detection_WT": round(f_wt, 4) if f_wt is not None else None,
            "delta_AD_minus_WT": round(f_ad - f_wt, 4) if (f_ad is not None and f_wt is not None) else None,
            "log2_ratio": round(float(np.log2((f_ad + 1e-6) / (f_wt + 1e-6))), 3)
            if (f_ad is not None and f_wt is not None) else None,
            "direction": ("up_in_AD" if (f_ad or 0) > (f_wt or 0)
                          else "down_in_AD" if (f_ad or 0) < (f_wt or 0) else "same"),
            "per_celltype": per_cell,
            "per_timepoint": {
                f"{dis}_{s4}": {
                    "n_cells": int(tp_cells[(dis, s4)]),
                    "detected": int(tp_hits[tid][(dis, s4)]),
                    "frac": round(tp_hits[tid][(dis, s4)] / tp_cells[(dis, s4)], 4)
                    if tp_cells[(dis, s4)] else None,
                }
                for (dis, s4) in sorted(tp_cells)
            },
        }

    outp = Path(os.environ.get("EXPR_OUT", ROOT / "docs/quantization/expression" / f"{TISSUE}.json"))
    outp.parent.mkdir(parents=True, exist_ok=True)
    outp.write_text(json.dumps(out, indent=2))

    print(f"pool cells: AD={sum(v.get('AD',0) for v in out['cell_composition'].values())} "
          f"WT={sum(v.get('WT',0) for v in out['cell_composition'].values())}")
    print(f"\n{'gene':9s} {'det_AD':>7s} {'det_WT':>7s} {'delta':>8s} {'log2':>6s}  direction")
    for sym, g in sorted(out["genes"].items(), key=lambda kv: -(kv[1]["delta_AD_minus_WT"] or 0)):
        print(f"{sym:9s} {g['detection_AD']:7.4f} {g['detection_WT']:7.4f} "
              f"{g['delta_AD_minus_WT']:+8.4f} {g['log2_ratio']:+6.2f}  {g['direction']}")
    print(f"\n-> {outp}")


if __name__ == "__main__":
    main()
