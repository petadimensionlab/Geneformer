# rds2h5ad — Seurat `.rds` → Geneformer `.h5ad` pipeline

Converts a mouse single-cell RNA-seq Seurat object into a Geneformer-ready
`.h5ad`. This folder contains steps 1–3 of the ADPD pipeline
(`analysis/mouse2human.md`): inspect → subsample/export → ortholog mapping +
`.h5ad` assembly. Steps 4–6 (tokenize/probe, figures, fine-tune) need a
torch/GPU environment and are **not** set up here.

---

## 1. Layout

```
rds2h5ad/                        (scripts+docs tracked; data/env/outputs ignored)
├── setup.sh                     builds/verifies the environment (idempotent)
├── run.sh                       ./run.sh <tissue> — runs 01 → 02 → 03
├── 01_inspect.R                 inspect the Seurat object
├── 02_export.R                  subsample + export raw counts (mtx + csv)
├── 03_map_and_build_h5ad.py     mouse→human orthologs, filter, build .h5ad
├── env.sh                       source this for $RSCRIPT/$PYTHON/... (generated)
├── manifest.json                environment provenance (generated)
├── .Rprofile renv.lock renv/    renv project (Seurat 5, R 4.6)
├── pyenv/                       Python 3.12 venv (anndata/numpy/pandas/scipy)
├── rds_data/<tissue>/           input .rds files, one folder per tissue
│   └── PD_LN/*.rds
└── ws/<tissue>/                 outputs per tissue
    ├── logs/01_inspect.log 02_export.log 03_map.log
    ├── export/data/{counts.mtx,genes.csv,cells.csv,obs.csv,export_meta.json}
    └── h5ad/<tissue>.h5ad (+ <tissue>_mapping_summary.json)
```

Supporting locations outside this folder:

| Path | Role |
|---|---|
| `Geneformer/geneformer_hf` | Geneformer checkout (token/gene dictionaries, V2-104M checkpoint) |
| `Geneformer/analysis/ref` | MGI mouse-human homology table + sha256 |
| `Geneformer/ref` | symlink → `analysis/ref` (03 looks for `parent.parent/ref`) |

## 2. Requirements (verified environment)

- macOS arm64 (M4, 128 GB RAM)
- R ≥ 4.4 — this machine: Homebrew **R 4.6.0**. `mouse2human.md` asks for
  R 4.4; R 4.6 was accepted by the project owner.
- Apple clang (Xcode CLT) + `gfortran` (Homebrew `gcc`) — the Seurat stack is
  compiled from source (Posit PPM serves **source only** on macOS, no binaries)
- `uv` (or python3.12) for the venv
- The rds must contain: an **RNA assay with a raw-count `counts` layer**
  (never the SCT default), a sample column with `_1/_2/_3` replicate
  suffixes, a cell-type column, and mouse gene symbols.

## 3. Setup

```bash
cd Geneformer/rds2h5ad
bash setup.sh                 # idempotent; first R-package pass ~5-40 min
```

What it builds (all under this folder unless noted):

| Component | Result |
|---|---|
| R base | system Homebrew R 4.6.0 |
| renv project | Seurat 5.5.1, SeuratObject 5.4.0, Matrix, jsonlite (+131 deps), lockfile at `renv.lock` |
| Python venv | `pyenv/` — CPython 3.12.9, anndata 0.13.2, numpy 2.5.2, pandas 3.0.5, scipy 1.18.0 |
| ref symlink | `Geneformer/ref → analysis/ref` (created if absent) |
| generated files | `env.sh`, `manifest.json`, cwd-independent `.Rprofile` |

`setup.sh` self-verifies at the end (R version, renv activation, Seurat load,
python imports) and exits non-zero on any failure.

## 4. Running the pipeline

```bash
cd Geneformer/rds2h5ad
./run.sh PD_LN
```

- Arg 1 = folder name inside `rds_data/`; the first `.rds` found there is used.
- Optional args: `./run.sh <tissue> <sample_col> <ct_col> [max_per_group] [assay]`
  — defaults `samples2 final_cell.type 200 RNA`. Pass `RNA` explicitly; the
  SCT "counts" layer is sctransform-corrected and would silently corrupt
  tokenization.
- `run.sh` re-reads the 01 log and **aborts before 02** if either requested
  column is not present in the object's meta.data.
- Outputs: `ws/<tissue>/{logs,export,h5ad}`. Every step is cached — re-running
  reuses existing outputs (delete an output to force recomputation).

Adding a new tissue: `mkdir rds_data/<name>`, drop its `.rds` in, run
`./run.sh <name> [columns]`. Column names are object-specific — read
`ws/<name>/logs/01_inspect.log` first (meta.data table + candidate columns).

Manual invocation (without `run.sh`):

```bash
source env.sh
mkdir -p "$ADPD_ROOT/PD_LN"/{logs,export,h5ad}
$RSCRIPT   01_inspect.R "$RDS" | tee "$ADPD_ROOT/PD_LN/logs/01_inspect.log"
$RSCRIPT   02_export.R "$RDS" "$ADPD_ROOT/PD_LN/export/data" samples2 final_cell.type 200 RNA \
  | tee "$ADPD_ROOT/PD_LN/logs/02_export.log"
$PYTHON    03_map_and_build_h5ad.py "$ADPD_ROOT/PD_LN/export/data" \
  "$ADPD_ROOT/PD_LN/h5ad/PD_LN.h5ad" "$GENEFORMER_DIR/geneformer" \
  | tee "$ADPD_ROOT/PD_LN/logs/03_map.log"
```

### Environment variables (`env.sh`)

| Variable | Meaning |
|---|---|
| `ADPD_TOOLS` | this folder |
| `R_BIN` | R install dir (`/opt/homebrew/bin`) |
| `RENV_PROJECT`, `R_PROFILE_USER` | renv activation (any cwd) |
| `RSCRIPT` | Rscript with the renv library active |
| `PYTHON` | venv python |
| `GENEFORMER_DIR` | Geneformer checkout (`../geneformer_hf`) |
| `GENEFORMER_MODEL` | checkpoint dir name (`Geneformer-V2-104M`, steps 4+) |
| `RDS` | auto-detected rds (convenience; `run.sh` resolves per tissue) |
| `ADPD_ROOT` | output base (default `ws/`) |
| `ADPD_PREFIX` | default `PD_LN` |

## 5. Verified run — PD_LN (2026-08-22)

Input: `rds_data/PD_LN/LNPD_soup.corrected_all.integrated_doublet.strictly.removed_scDblFinder_annotated_calculated.rds` (11 GB)

| Step | Result |
|---|---|
| 01 | 195,927 cells, 33,750 features, assays RNA+SCT (default SCT — pass RNA), 36 samples, 23 cell types, mouse symbols, counts integer-like: TRUE |
| 02 | 66,467 cells (828 donor×celltype groups, ≤200/group) × 17,139 genes |
| 03 | 15,071/17,139 genes retained (**87.9%**, expected 85–91%), donor leakage check: PASS, `ws/PD_LN/h5ad/PD_LN.h5ad` (828 MB), balanced train/eval/test crosstab over all 23 types |

## 6. Pitfalls & fixes baked into `setup.sh`

Discovered during build-up; all handled automatically, listed so nobody
"fixes" them back:

1. **PPM serves source only on macOS** — no prebuilt packages for mac; the
   Seurat stack compiles locally (fast on M4, ~5 min after the first failure
   round; first full build ~40 min).
2. **Homebrew R/gcc version drift** — the R bottle baked gcc-15 runtime paths
   into `Makeconf` while Homebrew `gcc` was 16; every C package failed at link
   with `library 'emutls_w' not found`. Fix: `setup.sh` symlinks the expected
   versioned lib dir to the one actually present (inside the Homebrew Cellar).
3. **Broken `mamba` binary** on this machine (conda/mamba mismatch) — `setup.sh`
   health-checks mamba/micromamba/conda before use (moot now: system R is used,
   no conda env).
4. **renv 1.2.4 API changes** — `init(prompt=)` removed and no `renv/bin/`
   shims are generated. Activation is done via `R_PROFILE_USER` +
   `RENV_PROJECT` (set in `env.sh`); `setup.sh` overwrites renv's
   relative-path `.Rprofile` with a cwd-independent one.
5. **`Rscript -e` double-unescaping** (R 4.6) — `-e` applies one extra level
   of string unescaping, breaking regex backreferences. `setup.sh` runs its
   verification from a temp file (file mode parses standardly).
6. **`03` ref path** — it looks for `Path(__file__).parent.parent/ref`; with
   the script in `rds2h5ad/` that resolves to `Geneformer/ref`, hence the
   symlink (keeps the recorded MGI sha256 check working).

## 7. Notes

- `rds2h5ad/` scripts and docs are tracked in git; `rds_data/` (11 GB rds),
  `renv/`, `pyenv/`, `ws/` and generated files (`env.sh`, `manifest.json`,
  `.Rprofile`, `renv.lock`) are git-ignored (workspace root `.gitignore`).
  The `Geneformer/ref` symlink is machine-specific and stays untracked —
  `setup.sh` recreates it.
- History: environment originally built at `~/adpd-tools`, scripts in
  `analysis/`, data in `input/PD_LN/`; all relocated here on 2026-08-22.
- The repo's pre-existing `Geneformer/.venv` (python 3.12 + geneformer +
  torch, from earlier runs) is untouched — it is the starting point for
  steps 4–6, whose scripts (`04_baseline.py`, `05_figures.py`,
  `06_finetune.py`) remain in `analysis/`.
