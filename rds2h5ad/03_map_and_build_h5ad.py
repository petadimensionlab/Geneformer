#!/usr/bin/env python
"""Assemble the R export into a Geneformer-ready .h5ad.

Mouse symbols -> human orthologs (MGI homology classes) -> human Ensembl IDs
(Geneformer's own gene_name_id dict) -> filtered to the V2 token vocabulary.

Usage:
  python 03_map_and_build_h5ad.py <export_dir> <out_h5ad> <geneformer_pkg_dir>
"""
from __future__ import annotations

import json
import os
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

export_dir = Path(sys.argv[1])
out_h5ad = Path(sys.argv[2])
gf_dir = Path(sys.argv[3])
ref_dir = Path(__file__).resolve().parent.parent / "ref"

SEED = 42
np.random.seed(SEED)

# ---------------------------------------------------------------- load export
print("loading export from", export_dir, flush=True)
X = sio.mmread(export_dir / "counts.mtx").T.tocsr()  # -> cells x genes
obs = pd.read_csv(export_dir / "obs.csv", index_col=0)
genes = pd.read_csv(export_dir / "genes.csv")["gene"].astype(str).to_numpy()
print(f"  {X.shape[0]} cells x {X.shape[1]} genes", flush=True)
assert X.shape == (obs.shape[0], genes.shape[0]), "shape mismatch across export files"

# ------------------------------------------------- mouse symbol -> human symbol
print("building mouse->human ortholog map (MGI homology classes)", flush=True)
MGI_URL = "https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt"
mgi_file = ref_dir / "mgi_homology.rpt"
if not mgi_file.exists():
    import urllib.request
    ref_dir.mkdir(parents=True, exist_ok=True)
    print(f"  downloading MGI homology table from {MGI_URL}", flush=True)
    urllib.request.urlretrieve(MGI_URL, mgi_file)

# record which version of a moving reference produced these results
import hashlib
digest = hashlib.sha256(mgi_file.read_bytes()).hexdigest()
expected = (ref_dir / "mgi_homology.rpt.sha256")
if expected.exists() and expected.read_text().strip() != digest:
    print(f"  NOTE: MGI table differs from the recorded version ({digest[:12]}...);"
          " gene retention may not match the reference run.", flush=True)
print(f"  MGI homology table sha256: {digest[:16]}...", flush=True)

hom = pd.read_csv(mgi_file, sep="\t", dtype=str)
hom.columns = [c.strip() for c in hom.columns]
key_col, org_col, sym_col = "DB Class Key", "Common Organism Name", "Symbol"

mouse = hom[hom[org_col].str.contains("mouse", case=False, na=False)]
human = hom[hom[org_col].str.contains("human", case=False, na=False)]

human_by_key = defaultdict(list)
for k, s in zip(human[key_col], human[sym_col]):
    human_by_key[k].append(s)
mouse_by_key = defaultdict(list)
for k, s in zip(mouse[key_col], mouse[sym_col]):
    mouse_by_key[k].append(s)

# Keep only unambiguous 1:1 homology classes. Many-to-many classes (expanded
# paralog families) cannot be mapped without an arbitrary choice, so we drop
# them and report how many rather than guessing.
m2h: dict[str, str] = {}
ambiguous = 0
for k, ms in mouse_by_key.items():
    hs = human_by_key.get(k, [])
    if len(ms) == 1 and len(hs) == 1:
        m2h[ms[0]] = hs[0]
    elif hs:
        ambiguous += len(ms)
print(f"  1:1 ortholog pairs: {len(m2h)};  genes in ambiguous classes: {ambiguous}", flush=True)

# ------------------------------------------- human symbol -> Ensembl -> token
name_id = pickle.load(open(gf_dir / "gene_name_id_dict_gc104M.pkl", "rb"))
token = pickle.load(open(gf_dir / "token_dictionary_gc104M.pkl", "rb"))
vocab = {g for g in token if g.startswith("ENSG")}
print(f"  Geneformer: {len(name_id)} symbol->ENSG, {len(vocab)} ENSG tokens", flush=True)

rows, ensg = [], []
stats = dict(no_ortholog=0, no_ensg=0, not_in_vocab=0, mapped=0)
for i, g in enumerate(genes):
    h = m2h.get(g)
    if h is None:
        stats["no_ortholog"] += 1
        continue
    e = name_id.get(h)
    if e is None:
        stats["no_ensg"] += 1
        continue
    if e not in vocab:
        stats["not_in_vocab"] += 1
        continue
    rows.append(i)
    ensg.append(e)
    stats["mapped"] += 1

print("gene mapping:", json.dumps(stats, indent=2), flush=True)
pct = 100 * stats["mapped"] / len(genes)
print(f"  retained {stats['mapped']}/{len(genes)} genes ({pct:.1f}%)", flush=True)

X = X[:, rows]
ensg = np.array(ensg)

# ------------------------------- collapse mouse paralogs landing on same ENSG
uniq, inverse = np.unique(ensg, return_inverse=True)
if len(uniq) < len(ensg):
    print(f"  collapsing {len(ensg) - len(uniq)} duplicate ENSG columns (summing counts)", flush=True)
    collapse = sp.csr_matrix(
        (np.ones(len(inverse)), (np.arange(len(inverse)), inverse)),
        shape=(len(inverse), len(uniq)),
    )
    X = (X @ collapse).tocsr()
    ensg = uniq
X = X.astype(np.float32)
X.eliminate_zeros()

# ---------------------------------------------------------- donor-level splits
# Replicate suffix (_1/_2/_3) within each condition x age group. Every split
# therefore sees all 12 condition-age groups, and no sample spans two splits.
# Column names come from the export manifest that 02_export.R wrote, so this
# script never assumes a particular tissue's naming. Env vars override.
_meta = json.loads((export_dir / "export_meta.json").read_text())
SAMPLE_COL = os.environ.get("ADPD_SAMPLE_COL", _meta.get("donor_col", "samples2"))
CELLTYPE_COL = os.environ.get("ADPD_CELLTYPE_COL", _meta.get("celltype_col", "final_cell.type"))
for _c, _n in ((SAMPLE_COL, "sample"), (CELLTYPE_COL, "cell type")):
    if _c not in obs.columns:
        raise SystemExit(
            f"{_n} column {_c!r} not found in obs.csv.\n"
            f"Available: {sorted(obs.columns)}\n"
            "Set ADPD_SAMPLE_COL / ADPD_CELLTYPE_COL, or re-run 02_export.R."
        )
print(f"  sample column: {SAMPLE_COL} | cell-type column: {CELLTYPE_COL}", flush=True)

rep = obs[SAMPLE_COL].astype(str).str.rsplit("_", n=1).str[-1]
split_map = {"1": "train", "2": "eval", "3": "test"}
obs = obs.copy()
obs["split"] = rep.map(split_map)
unmapped = obs["split"].isna().sum()
assert unmapped == 0, f"{unmapped} cells have an unrecognised replicate suffix"

leak = obs.groupby(SAMPLE_COL)["split"].nunique()
assert leak.max() == 1, "donor leakage: a sample appears in more than one split"
print("donor leakage check: PASS", flush=True)

# ------------------------------------------------------ tutorial obs contract
obs["cell_id"] = obs.index.astype(str)
obs["individual"] = obs[SAMPLE_COL].astype(str)
obs["celltype"] = obs[CELLTYPE_COL].astype(str)
obs["filter_pass"] = 1
obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel().astype(float)

required = ["cell_id", "individual", "celltype", "split", "n_counts", "filter_pass"]
# optional passengers: kept when present, so downstream analyses can use them
optional = [c for c in ("disease", "samples4", "orig.ident", "condition", "age",
                        "treatment", CELLTYPE_COL + "_narrow") if c in obs.columns]
wanted = required + optional
# tissues differ in which optional annotation columns exist
keep_obs = [c for c in wanted if c in obs.columns]
_missing = [c for c in required if c not in obs.columns]
if _missing:
    raise SystemExit(f"required obs columns missing: {_missing}")
var = pd.DataFrame({"ensembl_id": ensg}, index=ensg)

adata = ad.AnnData(X=X, obs=obs[keep_obs].copy(), var=var)
adata.obs_names = adata.obs["cell_id"].to_numpy()

# empty cells cannot be tokenized
nz = np.asarray(adata.X.sum(axis=1)).ravel() > 0
if (~nz).sum():
    print(f"  dropping {(~nz).sum()} cells with zero counts after gene filtering", flush=True)
    adata = adata[nz].copy()

print(adata, flush=True)
print(adata.obs.groupby(["split"], observed=True).size(), flush=True)
print(pd.crosstab(adata.obs["split"], adata.obs["celltype"]).T, flush=True)

out_h5ad.parent.mkdir(parents=True, exist_ok=True)
adata.write_h5ad(out_h5ad)
print("wrote", out_h5ad, flush=True)

summary = dict(
    export_dir=str(export_dir),
    n_cells=int(adata.n_obs),
    n_genes=int(adata.n_vars),
    gene_mapping=stats,
    retained_pct=round(pct, 2),
    n_samples=int(adata.obs["individual"].nunique()),
    n_celltypes=int(adata.obs["celltype"].nunique()),
    splits={k: int(v) for k, v in adata.obs["split"].value_counts().items()},
)
(out_h5ad.parent / (out_h5ad.stem + "_mapping_summary.json")).write_text(
    json.dumps(summary, indent=2)
)
print("=== DONE ===", flush=True)
