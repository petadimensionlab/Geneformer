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
     series. The upstream code indexes the pandas `var["ensembl_id_collapsed"]`
     Series with integer positions via `[...]`, which does **label-based**
     lookup on the gene-name index and selects wrong genes. The
     `patches/geneformer_multibackend.patch` fixes this by using
     `.iloc[coding_miRNA_loc]` (positional) in `tokenize_anndata`
     (`norm_factor_vector` and `coding_miRNA_ids`).
   - On Linux, the tokenizer can also fail with
     `OSError: Unable to synchronously open file (file signature not found)`
     when the input dir contains **hidden files** such as macOS AppleDouble
     metadata (`._*.h5ad`, `.DS_Store`, created on SMB/NFS or Mac→Linux copy).
     The patch skips any file whose name starts with `.` in `tokenize_files`,
     and `04_baseline.py` excludes hidden `.h5ad` files when globbing.
2. **Frozen-embedding + probe** (`analysis/04_baseline.py`) — **works** on MPS.
   - `analysis/smoke_mps.py` verifies MPS vs CPU outputs agree to ~1.9e-5.
3. **Fine-tuning** (`analysis/06_finetune.py`) — **works** on MPS (uses
   Hugging Face `Trainer`, which already supports MPS).
   - **Memory budget (M4 Max, 128 GB unified memory):** the collator pads each
     batch to its longest cell (PD_BM: median 1201 tokens, max 4096), and one
     train step at that worst-case length peaks at ~145 GiB with plain
     backward → MPS OOM (`max allowed: 182.78 GiB`) within a few steps.
     `analysis/measure_mps_batch.py` measures this empirically per batch size:

     | config (L=4096 worst case)          | peak MPS memory |
     |-------------------------------------|-----------------|
     | batch 8, no checkpointing           | ~145 GiB (OOM)  |
     | batch 4, no checkpointing           | ~68 GiB         |
     | **batch 8 + gradient checkpointing**| **~64 GiB**     |

Fix applied in `06_finetune.py`: keep
      `per_device_train_batch_size: 8` and enable
      `gradient_checkpointing: True` (`use_reentrant: False`) — same effective
      batch / optimizer dynamics, ~2.3x less activation memory. Runs at
      ~6 s/it (~4x faster than before, which was swapping) with an ETA of a few
      hours for the full epoch.
   - **Incomplete-checkpoint recovery:** if a prior fine-tune run is interrupted
     (e.g. MPS OOM) before `trainer.save_model()` writes weights, the
     `<TISSUE>/runs/<prefix>.../ksplitN/` dir may be left **empty**. `06_finetune.py`
     treats a checkpoint as reusable only if it actually contains saved weights
     (`pytorch_model.bin` / `model.safetensors` / `adapter_model.bin`); any empty
     or incomplete `ksplit*` dir is deleted automatically so a fresh fine-tune
     runs instead of erroring on a missing model file.
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

> **Multi-organ early-disease screens** — separate AD/PD in-silico-perturbation
> runs across several organs (smallint / brain / blood / spleen / liver / bone
> marrow) and hosts are documented in the
> **[In Silico Perturbation Wiki](docs/isp/README.md)** (unified output layout,
> shared `_isp_common.py`, per-organ configs & results).

## Setup (uv)

```bash
# 1. Download the geneformer package + weights into geneformer_hf/ (curl,
#    no git / git-lfs needed), then apply the multi-backend (MPS/CUDA/ROCm)
#    device patches.
./download.sh
cd geneformer_hf
git apply ../patches/geneformer_multibackend.patch   # cwd-independent, works for any user
cp ../patches/device.py geneformer/device.py         # new file (untracked)
cd ..
# If `git apply` fails (e.g. tokenizer.py was already edited by hand), apply
# the fixes directly by copying the ready-made files instead:
#   cp patches/tokenizer.py geneformer_hf/geneformer/tokenizer.py
#   cp patches/device.py    geneformer_hf/geneformer/device.py
# (tokenizer.py carries the .iloc positional-index and hidden-file-skip fixes)

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

Model weights are stored via git-lfs upstream, but `download.sh` does **not**
use git at all — it fetches every file (code, weights, dictionaries) directly
from `huggingface.co` with `curl`:

```bash
./download.sh              # V2-104M + package code + dictionaries
./download.sh --all        # all models (V1-10M, V2-104M, CLcancer, V2-316M)
# change the target model: ./download.sh --model V1-10M
```

Because it uses plain `curl` (not `git lfs`), it works on any host regardless
of whether git-lfs is installed — no LFS smudge hangs, no LFS-pointer files
(a common failure with `git lfs`), and no `safetensors: header too large` /
`pickle.UnpicklingError: invalid load key, 'v'` load errors. Output dir is
overridable via `GENEFORMER_DIR` and `HF_MIRROR`. Existing files are skipped,
so re-running is safe (idempotent).

**Every run ends with a full-file integrity check** (`verify_all`): it verifies
that `setup.py`, `geneformer/*.py` code, `model.safetensors`, `config.json`,
`generation_config.json`, `training_args.bin`, and all `geneformer/*.pkl`
dictionaries exist as **real data** (not LFS pointers) with a minimum size.
`pyproject.toml` is not tracked upstream (Geneformer ships `setup.py` only), so
`download.sh` auto-generates a minimal one after the download — no manual
`cp` step is needed. If anything is missing, `download.sh` fails with a
non-zero exit — so a "success" means the checkout is genuinely usable.

## Downloaded files

The Geneformer repo (`geneformer_hf/`) is populated by `download.sh` via direct
HTTPS download (no git needed). It fetches the package code, the V2-104M model
weights, and the package's token/median dictionaries. Files actually downloaded:

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

> The V1-10M, V2-104M_CLcancer, and V2-316M weights are **not** downloaded by
> default (only V2-104M is). To fetch any of them: `./download.sh --model V2-316M`
> or `./download.sh --all`.

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
| `analysis/04b_extract_embeddings.py` | embeddings only (no probe) — see below |
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
# Run from the repo root (so the relative path .venv/bin/python resolves).
cd /path/to/Geneformer

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

> **Must run from the repo root.** The scripts use the relative path
> `.venv/bin/python`, so running from another directory fails with
> `no such file or directory: .venv/bin/python` (the resolver itself is
> cwd-independent — it locates `input/` from the script path). The scripts'
> `ADPD_TISSUE=... cmd` syntax is POSIX; if your shell is **csh/tcsh** (which
> does not support `VAR=value cmd`), use `env` instead:
> ```bash
> env ADPD_TISSUE=PD_blood .venv/bin/python analysis/05_figures.py
> ```

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

### 04 vs `04b_extract_embeddings.py`

`04_baseline.py` runs **tokenize → extract embeddings → fit logistic probe** in
one go. The embedding step is the slow part (it runs the full 104M model over
every cell), so `04b_extract_embeddings.py` extracts **only the embeddings**
(no probe — `05_figures.py` refits the probe itself) and writes the same file
`05` reads:

```
<input>/<TISSUE>/results/embeddings/pretrained_cell_embeddings.csv
```

If `05_figures.py` fails with
`No such file or directory: .../pretrained_cell_embeddings.csv`, the embeddings
were never produced (04 was not run, or stopped before the embedding step). Run
the embedding stage (idempotent — reuses an existing `tokenized/`):

```bash
env ADPD_TISSUE=<TISSUE> .venv/bin/python analysis/04b_extract_embeddings.py
```

**How to tell whether tokenization finished and only embeddings remain** — check
for each artifact (all created by `04b`/`04`, in `<input>/<TISSUE>/`):

| Artifact | Meaning |
|---|---|
| `tokenized/<PREFIX>.dataset/` (a directory, not empty) | tokenization done |
| `results/embeddings/pretrained_cell_embeddings.csv` | embeddings done (05 is ready) |

- If **only the tokenized dir exists** (embeddings CSV missing), tokenization
  finished but the run stopped during embedding — just run `04b_extract_embeddings.py`;
  it reuses `tokenized/` and only does the embedding.
- If **neither exists**, run `04b_extract_embeddings.py` and it will both
  tokenize and embed.

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
