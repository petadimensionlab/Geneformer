# Geneformer V2-104M — Multi-backend (MPS / CUDA / ROCm) enablement

> **Language: English | [Japanese](README-jp.md)**

This workspace makes **Geneformer V2-104M** (from
[`ctheodoris/Geneformer`](https://huggingface.co/ctheodoris/Geneformer)) run on
three accelerator backends:

| Backend | Hardware | Notes |
|---|---|---|
| **MPS** | Apple Silicon (Mac) | verified in this workspace |
| **CUDA** | NVIDIA GPUs | original supported path |
| **ROCm** | AMD GPUs via WSL2 | implemented; not verified on this Mac |

## Device abstraction

A central module `geneformer/device.py` resolves the compute device once and
routes every hardcoded `"cuda"` call through it:

- **Auto-detect priority**: CUDA/ROCm (`torch.cuda.is_available()`) > MPS >
  CPU.
- **Override**: set `GENEFORMER_DEVICE=mps|cuda|cpu` to force a backend.

```python
from geneformer.device import get_device, get_device_obj, move_to_device, empty_cache
get_device()          # -> "mps" | "cuda" | "cpu"
get_device_obj()      # -> torch.device("mps") ...
```

All call sites that previously hardcoded `"cuda"` (emb_extractor,
perturber_utils, evaluation_utils, in_silico_perturber, classifier,
mtl/{train,eval_utils,utils}, and the `analysis/` scripts) now route through
`get_device_obj()` / `move_to_device()` / `empty_cache()`.

### Why ROCm == CUDA

ROCm exposes the **CUDA API**: an AMD ROCm build of PyTorch reports
`torch.cuda.is_available() == True`, and `torch.cuda.*` calls hit the AMD GPU
under the hood. That means the same `device="cuda"` code path drives AMD ROCm
hardware — no separate device branch is needed. The abstraction therefore
supports ROCm "for free" whenever the torch build is the ROCm one.

## MPS verification (this Mac)

Ran end-to-end on Apple Silicon:

1. **Tokenization** of `input/PD_smallint/PD_smallint.h5ad` (64,614 cells,
   25 cell types) → `analysis_ws/tokenized/PD_smallint.dataset` — **works**.
   - Required fixes: `use_h5ad_index=True` (var index already holds ENSG IDs),
     and `tokenize_anndata` positional indexing when the var index is a string
     series.
2. **Frozen-embedding + probe** (`analysis/04_baseline.py`) — **works** on MPS.
   - `analysis/smoke_mps.py` verifies MPS vs CPU outputs agree to ~1.9e-5.
3. **Fine-tuning** (`analysis/06_finetune.py`) — **works** on MPS (uses
   Hugging Face `Trainer`, which already supports MPS).
4. **In silico perturbation** (`analysis/07_in_silico_perturbation.py`, ported
   from the tutorial notebook) — **works** on MPS.

## In silico perturbation — nproc / datasets / max_ncells (learned the hard way)

`InSilicoPerturber.perturb_data` builds a perturbation `Dataset` with
`dataset.map(make_group_perturbation_batch, num_proc=nproc)`. Three separate
pitfalls were hit and resolved; `analysis/07_in_silico_perturbation.py` applies
all of them and prints warnings / RAM estimates at runtime.

1. **`datasets>=5` hangs `perturb_data`** — the `dataset.map` step stalls with
   *any* nproc (observed for both `nproc=8` and `nproc=1`), consuming a core at
   0% CPU and never returning. **Root cause fix: pin `datasets==4.0.0`.** Do not
   try to work around it with nproc; it only masks the symptom. `07` prints a
   loud warning if `datasets>=5` is detected.

2. **`nproc`** — with `datasets==4.x` the map function (a nested closure aliasing
   `self.tokens_to_perturb`, etc.) *does* pickle under `spawn` and nproc>1 can
   run. **Empirically verified here:** the exact closure shape of
   `make_group_perturbation_batch` (delete-index + optional delete_indices) was
   run through `dataset.map` with `num_proc=1/2/4` on `datasets==4.0.0`; all
   three produced identical, correct output. So nproc>1 is functional with
   datasets 4.x. `07` still defaults to **`nproc=1`** (verified-stable),
   overridable with `IS_NPROC`. Reasons nproc=1 is the recommended default:
   - Workers add RAM that competes with the MPS forward pass on the same host.
   - With 2+ genes (or `combos>0`) the per-cell variant count grows, so the map
     output becomes large and type-inconsistent under multiprocessing — nproc=1
     sidesteps that fragility.
   - **Parallelism only pays off on large datasets.** A 4-cell benchmark showed
     nproc=1 at ~0.01s and nproc=2/4 at ~0.09s (spawn + interprocess overhead
     dominates at small scale). For single-gene perturbations the map is quick
     serially and the MPS forward pass dominates regardless of nproc.
   - `07` prints a **warning whenever `nproc>1`**, reminding you datasets<5 and
     modest `max_ncells` are required.

   Concretely, for **2+ perturbed genes**: nproc=1 remains the right choice.
   Doubling the gene list roughly multiplies per-cell variants (deletion mode:
   1 deleted gene -> 1 variant per cell; N deleted genes, `combos=0` -> up to N
   variants per cell), so memory scales with the gene count and the forward pass
   stays on the single MPS process. Parallelising the *map* does not help the
   dominant MPS forward cost.

   How to keep it fast with more genes: use a smaller `max_ncells` (below) or
   `combos=0`; prefer fewer, biologically-targeted genes over a large list.

4. **`max_ncells` / OOM** — `perturb_data` materialises the whole perturbation
   dataset in RAM (~ `n_cells * n_variants * seq_len`). A 2000-cell run OOM'd
   this MPS host (≈28–32 GB swap touched). `07` defaults `max_ncells=200`
   (`IS_MAX_CELLS` to override) and prints an RAM estimate from
   `estimate_perturb_ram()` before perturbing. Rule of thumb:
   `bytes ≈ n_cells * n_variants * seq_len * 10`; keep it well under free RAM.
   ```
   n_variants = max(n_genes,1)                      # combos=0
   n_variants = C(n_genes, combos+1)                # combos>0 (binomial, explodes)
   ```
   Single-gene `combos=0` is safe; 2+ genes × `combos>0` blows up fast.

5. **Restrict `genes_perturbed` in stats too** — `InSilicoPerturberStats(...,
   genes_perturbed="all")` re-scans every gene in the vocabulary, which is
   extremely slow on MPS. It must match the genes you actually perturbed. `07`
   sets `genes_perturbed=genes_to_perturb`.

Env vars used by `07`: `IS_NPROC` (default 1), `IS_MAX_CELLS` (default 200).

## Setup (uv)

```bash
# 1. Get the geneformer package (git-lfs clone, ~7 GB with weights) into
#    geneformer_hf/, then apply the multi-backend (MPS/CUDA/ROCm) device patches.
#    (`./download.sh` options are documented below in "Model weights are
#    LFS-tracked".)
#    IMPORTANT: the patch paths are relative to geneformer_hf/, so cd there
#    first. Running `patch` from the repo root fails with "File to patch:".
./download.sh
cd geneformer_hf
patch -p1 < ../patches/geneformer_multibackend.patch
cp ../patches/device.py geneformer/device.py        # new file (untracked)
cd ..

# 2. Create env
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python -e geneformer_hf
# geneformer requires transformers==4.46 (setup.py pins an unpinned range,
# so uv's resolver may pick 5.x — pin it back):
uv pip install --python .venv/bin/python "transformers==4.46.3"
# datasets MUST be < 5. datasets>=5 makes InSilicoPerturber.perturb_data hang
# in its dataset.map step (any nproc):
uv pip install --python .venv/bin/python "datasets==4.0.0"
```

> **Note on re-running setup:** `download.sh` is idempotent and re-uses an
> existing `geneformer_hf/`. If you already applied the patch once, applying it
> again fails with "patch does not apply" (it is already applied) — that is
> expected. Skip the patch steps if `geneformer/device.py` already exists.

> **Pinned versions:** `transformers==4.46.3` (5.x breaks `SpecialTokensMixin`)
> and **`datasets==4.0.0`** (`datasets>=5` hangs `perturb_data`'s `dataset.map`
> — see the in silico perturbation notes below). Both are required for the
> standard pipeline to run reliably.

Model weights are LFS-tracked. Use the provided `download.sh`:

```bash
./download.sh                # V2-104M + dictionaries into geneformer_hf/ (git-lfs preferred)
./download.sh --force-https  # force direct HTTPS download (no git-lfs / it fails)
./download.sh --all          # all models (V1-10M, V2-104M, CLcancer, V2-316M)
# change the target model: ./download.sh --model V1-10M
```

`download.sh` prefers git-lfs and falls back to `curl` against
`huggingface.co` when git-lfs is unavailable. It is overridable via the
`GENEFORMER_DIR` (output dir), `GF_REPO_URL`, and `HF_MIRROR` environment
variables. Existing files are skipped, so re-running is safe (idempotent).

## Downloaded files

The Geneformer repo (`geneformer_hf/`) is a git-lfs clone of
`ctheodoris/Geneformer`. `git lfs pull` fetches the binary weights and the
package's token/median dictionaries. Files actually downloaded:

| File | Size | Purpose |
|---|---|---|
| `geneformer_hf/Geneformer-V2-104M/model.safetensors` | 417,571,156 B (~418 MB) | V2-104M pretrained weights |
| `geneformer_hf/Geneformer-V2-104M/config.json` | 590 B | model config (18 layers, hidden 1152, vocab 20275) |
| `geneformer_hf/Geneformer-V2-104M/generation_config.json` | 90 B | generation config |
| `geneformer_hf/Geneformer-V2-104M/training_args.bin` | 5,496 B | training args |
| `geneformer_hf/geneformer/token_dictionary_gc104M.pkl` | 425,590 B | Ensembl ID ↔ token map (V2) |
| `geneformer_hf/geneformer/gene_median_dictionary_gc104M.pkl` | 1,512,661 B | gene normalization factors (V2) |
| `geneformer_hf/geneformer/ensembl_mapping_dict_gc104M.pkl` | 3,957,652 B | Ensembl ID collapse/mapping (V2) |
| `geneformer_hf/geneformer/gene_name_id_dict_gc104M.pkl` | 1,660,882 B | Ensembl ID ↔ gene name (V2) |

> The V1-10M, V2-104M_CLcancer, and V2-316M directories are also present in
> the clone but their weights were **not** pulled (only V2-104M was requested).
> To fetch any of them: `git lfs pull --include="Geneformer-V2-316M/*"`.

## Environment (uv + venv)

The pipeline was verified against this environment (see **Setup (uv)** above
for the install commands):

```
.venv/                      1.5 GB  uv-created Python 3.12 environment
.venv/bin/python            3.12.9
torch                       2.13.0  (MPS built: True, MPS available: True)
transformers                4.46.3  (pinned — 5.x breaks SpecialTokensMixin)
datasets                    4.0.0   (pinned — 5.x hangs InSilicoPerturber.perturb_data)
```

## Generated artifacts (`analysis_ws/`)

| Path | Contents |
|---|---|
| `analysis_ws/tokenized/PD_smallint.dataset` | 64,614 tokenized cells (Hugging Face dataset) |
| `analysis_ws/results/mps_verify/mps_verify_cell_embeddings.csv` | 2,000 cell embeddings (768-d) on MPS |
| `analysis_ws/results/mps_verify/mps_verify_summary.json` | verify metrics (acc 0.8265, macro F1 0.6883) |
| `analysis_ws/runs/260821_geneformer_cellClassifier_mpsft/ksplit1/checkpoint-8` | fine-tuned checkpoint (8 steps) on MPS |
| `analysis_ws/results/isp/*_raw.pickle` | in silico perturbation intermediate (per-batch) |
| `analysis_ws/results/isp/isp_state_embs.pkl` | state embeddings for perturbation |

Analysis scripts live in `analysis/`:

| Script | Role |
|---|---|
| `analysis/04_baseline.py` | tokenize → frozen embeddings → logistic probe |
| `analysis/06_finetune.py` | fine-tune cell classifier |
| `analysis/07_in_silico_perturbation.py` | tutorial IS perturbation ported to V2 |
| `analysis/smoke_mps.py` | MPS vs CPU numerical correctness |
| `analysis/verify_mps_fast.py` | fast frozen-embedding + probe verify |
| `analysis/verify_mps_finetune.py` | fast fine-tune verify |

### Running per tissue (04–07)

Each of `04_baseline.py` / `05_figures.py` / `06_finetune.py` /
`07_in_silico_perturbation.py` resolves its workspace automatically from a
tissue folder under `input/<TISSUE>/h5ad/*.h5ad` (the `_resolve_tissue.py`
helper, bundled in `analysis/`). Pick the tissue with `ADPD_TISSUE`; the
scripts use `input/<TISSUE>/` as the workspace root and write
`tokenized/`, `results/`, `runs/` inside it:

```bash
export GENEFORMER_DIR=$PWD/geneformer_hf        # Geneformer checkout + V2-104M
export GENEFORMER_MODEL=Geneformer-V2-104M

# one tissue (pick any PD_*/AD_* under input/, e.g. PD_blood / AD_blood)
ADPD_TISSUE=PD_blood  .venv/bin/python analysis/04_baseline.py        # tokenize + embeddings + probe
ADPD_TISSUE=PD_blood  .venv/bin/python analysis/05_figures.py          # figures (reuses cached embeddings)
ADPD_TISSUE=PD_blood  .venv/bin/python analysis/06_finetune.py         # fine-tune cell classifier
ADPD_TISSUE=PD_blood  .venv/bin/python analysis/07_in_silico_perturbation.py

# OR equivalent via a single tissue:
#   ADPD_ROOT=$PWD/input/PD_blood  (legacy)  — both work
```

> **Notes**
> - `ADPD_TISSUE` must name a folder under `input/` that holds
>   `h5ad/*.h5ad` with the required obs columns (`cell_id`, `individual`,
>   `celltype`, `split`). All checked-in tissues (`PD_*`, `AD_*`) satisfy this.
> - `ADPD_PREFIX` remains a back-compat alias for `ADPD_TISSUE`; setting
>   `ADPD_ROOT` bypasses the resolver entirely.
> - `06_finetune.py` splits held-out replicates by their `_1/_2/_3` suffix and
>   only requires all three cohorts to be non-empty (PD_* have 12/rep,
>   AD_* typically 9–10/rep).
> - Each script is cached: re-running reuses existing `tokenized/`, embeddings,
>   and checkpoints. Delete an output to force recomputation.

## Converting a Seurat `.rds` to a Geneformer `.h5ad` (R → Python)

To go from an **R-generated Seurat object (`.rds`)** to an `.h5ad` that this
Geneformer pipeline can tokenize, use the **`rds2h5ad/`** converter (an
R/Python bridge that inspects the object, subsamples + exports raw counts, then
maps mouse→human orthologs and assembles the `.h5ad`).

- Full documentation: **[`rds2h5ad/rds2h5ad.md`](rds2h5ad/rds2h5ad.md)**
  (a Japanese version is also available:
  **[`rds2h5ad/rds2h5ad-jp.md`](rds2h5ad/rds2h5ad-jp.md)**)
  - The `.rds` must expose an `RNA` assay with a raw-integer `counts` layer
    (never the SCT default), a sample column with `_1/_2/_3` replicate suffixes,
    a cell-type column, and mouse gene symbols.
  - Verified end-to-end on a mouse `PD_LN` object (195,927 cells → 66,467 after
    subsampling → `PD_LN.h5ad`, 87.9% gene retention).

`rds2h5ad/` covers steps 1–3 (inspect / export / ortholog-map). Steps 4–6
(tokenize + probe, figures, fine-tune) run in this workspace's `analysis/`
scripts on the resulting `.h5ad`.

## ROCm on Windows / WSL2

**Not verifiable on this Mac** (macOS cannot run WSL2). The implementation is
complete — the abstraction already drives AMD ROCm through the CUDA path — but
the WSL2 environment itself must be set up on a Windows machine with an AMD GPU.

To run on ROCm under WSL2:

```bash
# In Windows, ensure WSL2 + an Ubuntu distro are installed.
# In WSL2 Ubuntu (AMD ROCm drivers already installed on Windows side):
wget https://repo.radeon.com/amdgpu-install/latest/ubuntu/jammy/amdgpu-install_*.deb
sudo apt install -y ./amdgpu-install_*.deb
sudo amdgpu-install --usecase=rocm
# then install PyTorch with the ROCm wheel:
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python torch --index-url \
    https://download.pytorch.org/whl/rocm6.2
uv pip install --python .venv/bin/python -e geneformer_hf
uv pip install --python .venv/bin/python "transformers==4.46.3"
# ROCm build reports torch.cuda.is_available() == True -> device auto = cuda
GENEFORMER_DEVICE=cuda python analysis/04_baseline.py
```

Expected result: `get_device() == "cuda"` on the ROCm build (CUDA API), so all
existing code paths run on the AMD GPU unchanged.
