#!/usr/bin/env python
"""Presence scan for the candidate NULL gene set (criteria.md E1, selection step).

criteria.md section 2-1 says: pick 20-30 genes that are expressed in the target
cells but biologically unrelated to the disease, run the same deletion
simulation, and take the 95th percentile of the resulting Shift distribution as
the E1 threshold. This script produces the evidence for WHICH candidates are
expressed, so the final set is chosen from measured presence and not by taste.

It also reports the presence of ACTB and GAPDH, which criteria.md E3 needs as
negative controls, and of the hypothesis (PD) genes, so the null set can be
matched to them on expression instead of being systematically higher or lower.

Run (same environment as 07g):
    ADPD_TISSUE=PD_atlas GENEFORMER_DIR=$PWD/geneformer_hf \
    .venv/bin/python analysis/07i_pd_atlas_null_genes.py

Writes nothing; prints a table.
"""
import collections
import os
import pickle
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
os.environ.setdefault("GENEFORMER_DIR", str(REPO / "geneformer_hf"))
os.environ.setdefault("ADPD_TISSUE", "PD_atlas")
sys.path.insert(0, str(REPO / "analysis"))
sys.path.insert(0, os.environ["GENEFORMER_DIR"])

from _isp_common import load_gene_dicts  # noqa: E402
from _resolve_tissue import resolve as _resolve_tissue  # noqa: E402

from datasets import load_from_disk  # noqa: E402
from geneformer import TOKEN_DICTIONARY_FILE  # noqa: E402

POOLS = [p.strip() for p in os.environ.get("IS_POOLS", "DMNX_Neu,GPI_Neu").split(",") if p.strip()]
MIN_FRAC = float(os.environ.get("IS_MIN_FRAC", "0.20"))

# Candidates: expressed in brain tissue, but not implicated in PD or
# neurodegeneration. Grouped by why they are considered unrelated.
CANDIDATES = [
    # negative controls that E3 requires
    "ACTB", "GAPDH",
    # housekeeping: cytoskeleton, glycolysis, mitochondria, translation
    "TUBA1B", "TUBB", "ACTN4", "PFN1", "CFL1", "MYL6", "VIM",
    "PGK1", "ENO1", "SLC2A1", "ATP1A1", "ATP5F1A",
    "RPLP0", "RPS18", "RPL13A", "EEF1A1", "RPS27A",
    # signalling / chaperone, unrelated to PD specifically
    "CALM1", "YWHAZ", "HSP90AA1",
    # second batch (2026-09-16): general cellular machinery, to see whether a
    # 20-30 gene null set exists at all at the 20% presence rule. snRNA-seq
    # reads nuclear transcripts preferentially, so housekeeping genes such as
    # ACTB and GAPDH sit BELOW 20% here; nuclear and translation machinery is
    # the remaining place to look.
    "NPM1", "HNRNPA1", "SRSF1", "POLR2A", "UBA1", "UBB", "UBC",
    "PSMA3", "PSMB4", "PSMD1", "PSMC4",
    "EEF2", "RPS3", "RPL10A", "RPL23A", "RPS24",
    "CANX", "HSPA8", "HSPA5", "CCT5", "TCP1",
    "ATP6V1A", "ATP6V0C", "SDHA", "UQCRC2", "COX5A", "NDUFA4", "ATP5F1B",
    "GNAI2", "GNB1", "PPP1CA", "PPP2CA", "CSK", "RAN",
    "GOT1", "MDH1", "IDH3A", "PDHA1", "ALDOA", "TPI1",
    "ARPC2", "CAPZA2", "DYNLL1", "GSN", "HPRT1", "TBP",
    # lineage markers of tissues that are not present in these pools
    "KRT8", "KRT18", "EPCAM", "CD3E", "CD19", "PTPRC", "COL1A1",
    "HBB", "ALB", "INS", "AMY1A", "CRYAA", "RHO", "OR2H1", "TAS1R1",
    "GFAP", "AQP4", "S100B",
    # hypothesis genes, for expression matching
    "DNAJC13", "RBFOX3", "LRRK2", "KCNJ6", "TMEM175", "FBXO7", "SNCA",
    "ATP13A2", "NEFH", "SNAP25", "PRKN", "NEFM", "PINK1", "VPS35", "NEFL",
    "MAPT", "CALB1",
]

tok2gene, ensg2name, name2ensg = load_gene_dicts()
vocab = set(tok2gene.values())
resolved = {}
for sym in CANDIDATES:
    ensg = name2ensg.get(sym)
    if ensg and ensg in vocab:
        resolved[sym] = ensg

ROOT, PREFIX = _resolve_tissue()
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"
tokdict = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
gene_token = {s: tokdict[e] for s, e in resolved.items()}
want = set(gene_token.values())
tok2sym = {v: k for k, v in gene_token.items()}

tokens = load_from_disk(str(TOKENIZED))
print(f"{len(tokens)} cells; {len(resolved)}/{len(CANDIDATES)} candidates in the V2 vocab", flush=True)

n_cells = collections.Counter()
n_expr = collections.defaultdict(collections.Counter)
for ex in tokens:
    key = (ex["celltype"], ex["disease"])
    n_cells[key] += 1
    for t in want.intersection(ex["input_ids"]):
        n_expr[key][tok2sym[t]] += 1


def frac(pool, disease, sym):
    n = n_cells.get((pool, disease), 0)
    return (n_expr[(pool, disease)][sym] / n) if n else 0.0


for pool in POOLS:
    n_pd = n_cells.get((pool, "Parkinson disease"), 0)
    n_norm = n_cells.get((pool, "normal"), 0)
    print(f"\n===== {pool}: PD {n_pd} cells / normal {n_norm} cells =====", flush=True)
    print(f"{'gene':10} {'presence_PD':>12} {'passes':>7}", flush=True)
    for sym in sorted(resolved, key=lambda s: -frac(pool, "Parkinson disease", s)):
        f = frac(pool, "Parkinson disease", sym)
        print(f"{sym:10} {f:12.3f} {'yes' if f >= MIN_FRAC else 'no':>7}", flush=True)

print("\n=== candidates that pass the presence filter in ALL pools ===")
passed = [
    s for s in sorted(resolved, key=lambda s: -frac(POOLS[0], "Parkinson disease", s))
    if all(frac(p, "Parkinson disease", s) >= MIN_FRAC for p in POOLS)
]
print(", ".join(passed), flush=True)
print(f"count: {len(passed)}", flush=True)