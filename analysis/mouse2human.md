# ADPD Geneformer pipeline

Classify cell types in **mouse** single-cell RNA-seq using **Geneformer**, a
foundation model pretrained on ~95M *human* cells.

Input is a Seurat `.rds`; output is a held-out classification report, confusion
matrix and UMAPs. An optional fine-tuning stage adapts the model itself.

The pipeline exists because three things stand between a Seurat object and
Geneformer: the file is R-serialized, the genes are mouse symbols, and the model
speaks human Ensembl IDs.

```
your.rds ──R──> counts.mtx ──orthologs──> .h5ad ──tokenize──> dataset
                                                                 │
                                              frozen embeddings ─┴─> probe ─> results
                                                   (optional) fine-tune ────> results
```

---

## 1. Requirements

| | |
|---|---|
| OS | Linux (x86_64 or aarch64) |
| GPU | NVIDIA + CUDA — needed for steps 4 and 6 only |
| RAM | ~80 GB. A 15 GB `.rds` expands to 40–60 GB when loaded |
| Disk | ~50 GB (6 GB model, ~2 GB per tissue of intermediates) |
| Tools | `git`, `git-lfs`, `curl` |

Steps 1–3 and 5 are CPU-only. Steps 1–2 are the RAM-hungry ones.

## 2. Install

```bash
./setup.sh ~/adpd-tools
```

Installs at user level, no root: R 4.4 + Seurat 5 (micromamba), the Geneformer
repo + V2-104M checkpoint (git-lfs), and a Python 3.12 environment with a torch
build matched to your driver. Takes 20–40 minutes; safe to re-run.

Then, before each session:

```bash
export ADPD_TOOLS=~/adpd-tools
export RSCRIPT=$ADPD_TOOLS/renv/bin/Rscript
export PYTHON=$ADPD_TOOLS/pyenv/bin/python
export GENEFORMER_DIR=$ADPD_TOOLS/Geneformer
```

## 3. What your `.rds` must contain

Checked by step 1. The pipeline needs:

| Requirement | Why |
|---|---|
| An **`RNA` assay with a `counts` layer** holding raw integer UMIs | Geneformer's rank-value encoding must not receive pre-normalized data |
| A **sample/animal column**, e.g. `samples2` | splits are per-animal, never per-cell |
| Sample names ending `_1` / `_2` / `_3` | the replicate index drives the train/eval/test split |
| A **cell-type column**, e.g. `final_cell.type` | the classification target |
| Mouse gene symbols as rownames | input to ortholog mapping |

Column *names* are arguments, not assumptions. You pass them to step 2; step 3
reads them back from the `export_meta.json` that step 2 writes, so they flow
through automatically. Override with `ADPD_SAMPLE_COL` / `ADPD_CELLTYPE_COL` if
you need to.

Any additional annotation columns present in your object (`disease`, `age`,
`treatment`, ...) are carried through to the `.h5ad` and into the tokenized
dataset, so later analyses can use them without re-running anything. Absent ones
are skipped silently — nothing is required beyond the table above.

> **The most dangerous failure mode:** if the object's default assay is `SCT`,
> its `counts` layer holds sctransform-corrected, depth-normalized values.
> Tokenizing those corrupts the ranking **without raising an error**. Always
> pass `RNA` explicitly.

## 4. Run

Set up per-tissue variables:

```bash
export ADPD_ROOT=~/Documents/my_tissue      # workspace for this tissue
export ADPD_PREFIX=my_tissue                # short name for output files
export RDS=/path/to/YourTissue.rds
mkdir -p $ADPD_ROOT/{logs,export,h5ad,tokenized,results}
```

### Step 1 — Inspect *(~5 min, CPU)*

```bash
$RSCRIPT scripts/01_inspect.R "$RDS" 2>&1 | tee $ADPD_ROOT/logs/01_inspect.log
```

**Read the log before continuing.** Confirm the assay names, your sample and
cell-type column names, that gene names are mouse symbols, and that counts are
integer-like.

### Step 2 — Export *(~5 min, CPU, high RAM)*

```
02_export.R <rds> <outdir> <sample_col> <celltype_col> <max_per_group> <assay>
```

```bash
$RSCRIPT scripts/02_export.R "$RDS" $ADPD_ROOT/export/data \
    samples2 final_cell.type 200 RNA 2>&1 | tee $ADPD_ROOT/logs/02_export.log
```

`max_per_group` caps cells per sample × cell type. 200 keeps rare types whole
while bounding the matrix; raise it if your classes are small.

### Step 3 — Orthologs + `.h5ad` *(~3 min, CPU)*

```bash
$PYTHON scripts/03_map_and_build_h5ad.py \
    $ADPD_ROOT/export/data $ADPD_ROOT/h5ad/$ADPD_PREFIX.h5ad \
    $GENEFORMER_DIR/geneformer 2>&1 | tee $ADPD_ROOT/logs/03_map.log
```

Maps mouse symbols to human orthologs (MGI 1:1 classes) and filters to
Geneformer's vocabulary. The MGI table downloads automatically to `ref/` if
absent; its SHA-256 is recorded so you can tell whether a later run used the
same reference.

Expect **85–91% gene retention**. The report prints where the losses went.

### Steps 4 — Tokenize, embed, probe *(~1 h, GPU)*

```bash
$PYTHON scripts/04_baseline.py 2>&1 | tee $ADPD_ROOT/logs/04_baseline.log
```

Freezes Geneformer, extracts one 768-d vector per cell, fits a logistic probe on
the train samples and evaluates on held-out samples. Writes
`results/tables/adpd_baseline_*`.

### Step 5 — Figures *(~3 min, CPU)*

```bash
$PYTHON scripts/05_figures.py 2>&1 | tee $ADPD_ROOT/logs/05_figures.log
```

Confusion heatmap, a true-vs-predicted UMAP, and an error map. Reuses the cached
embeddings.

### Step 6 — Fine-tune *(optional, ~2 h, GPU)*

```bash
$PYTHON scripts/06_finetune.py 2>&1 | tee $ADPD_ROOT/logs/06_finetune.log
```

Trains the top 6 transformer layers plus a classification head for one epoch and
writes `adpd_model_comparison.csv` — frozen versus fine-tuned on identical
held-out cells.

**Every step is cached**: re-running reuses existing output. Delete an output to
force recomputation. Run long steps under `tmux`.

## 5. Reading the results

- **Macro F1, not accuracy**, when classes are imbalanced. Accuracy is dominated
  by whatever is abundant.
- **The confusion matrix is the real output.** A coherent off-diagonal block
  means one hard sub-problem; scattered errors mean something is wrong.
- **Your effective sample size is the number of animals per split**, typically
  12 — not the number of cells. Cell counts in the thousands do not license
  narrow confidence intervals.
- **This is not zero-shot.** The pretrained checkpoint has no cell-type head; a
  supervised probe is trained downstream. "Without fine-tuning" ≠ "label-free".

## 6. Known limits

- **Ortholog mapping is lossy and non-random.** ~10–14% of genes have no clean
  1:1 human ortholog, and the losses concentrate in fast-evolving families —
  immune receptors, secretory proteins. Mouse-specific genes (e.g. `Ly6c`) are
  simply absent.
- **Rare classes are unreliable.** Below ~100 test cells, per-class F1 is noise.
  In fine-tuning, plain cross-entropy may abandon such a class entirely
  (F1 = 0.000). Add class weighting if that matters to you — note the frozen
  probe already uses `class_weight="balanced"`, so the two stages differ in this
  respect.
- **One fold.** Rotate which replicate is held out and report a spread before
  drawing conclusions.
- **The labels are usually your own annotations**, so a good score validates the
  conversion rather than discovering biology.

## 7. Layout

```
scripts/01_inspect.R              inspect the object (run first, always)
scripts/02_export.R               subsample + export raw counts
scripts/03_map_and_build_h5ad.py  orthologs -> .h5ad
scripts/04_baseline.py            tokenize -> frozen embeddings -> probe
scripts/05_figures.py             confusion matrix + UMAPs
scripts/06_finetune.py            optional fine-tuning + comparison
ref/mgi_homology.rpt              MGI mouse-human homology (auto-downloaded)
setup.sh                          installs R, Geneformer, Python env
```

Configuration is entirely by environment variable — no paths or column names
are hardcoded:

| Variable | Meaning | Default |
|---|---|---|
| `ADPD_ROOT` | workspace for this tissue | `~/Documents/ADPD_analysis` |
| `ADPD_PREFIX` | short dataset name, used in output filenames | `dataset` |
| `GENEFORMER_DIR` | the Geneformer checkout | set by `setup.sh` |
| `GENEFORMER_MODEL` | checkpoint directory name | `Geneformer-V2-104M` |
| `ADPD_SAMPLE_COL` | override the sample column | from `export_meta.json` |
| `ADPD_CELLTYPE_COL` | override the cell-type column | from `export_meta.json` |

## Verified on

Two tissues of a mouse Parkinson's study (36 animals each, 4 conditions ×
3 timepoints), run end to end on an NVIDIA GB10 / aarch64 / CUDA 13 host:

| | cells | classes | gene retention | frozen | fine-tuned |
|---|---|---|---|---|---|
| PD Blood | 65,013 | 15 | 90.7% | 0.955 / 0.845 | 0.976 / 0.869 |
| PD SmallInt | 64,614 | 25 | 86.6% | 0.867 / 0.826 | see `results/` |

(accuracy / macro F1 on held-out animals)

The two objects differed in metadata columns — one lacked a column the other
had — which is why steps 3 and 4 adapt to whatever your object contains rather
than requiring a fixed schema.

## 8. Troubleshooting

| Symptom | Cause |
|---|---|
| `Geneformer checkpoint not found` | `GENEFORMER_DIR` unset or `setup.sh` not run |
| `KeyError` on a metadata column | column name differs — check the step-1 log |
| `assert leak.max() == 1` fails | a sample appears in two splits; check the replicate suffix |
| Retention below ~80% | gene names may not be mouse symbols — check step 1 |
| `counts integer-like: FALSE` | you are reading a normalized layer; pass `RNA` |
| CUDA OOM | lower `forward_batch_size` in the script |
