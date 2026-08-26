#!/usr/bin/env python
"""Rebuild the corrupted tokenized PD_spleen dataset via split-and-parallel
tokenization, then concatenate and verify.

Why (the corruption):
    input/PD_spleen/tokenized/PD_spleen.dataset is unreadable by
    datasets.load_from_disk (ArrowInvalid: "Expected to read 4161536 metadata
    bytes, but only read 4"). 5 arrow files (incl. the main data file) lost a
    4 MB-aligned tail block during Mac->Linux transfer — the trailing IPC
    footer/metadata is gone, so the files cannot be read even with open_stream.
    The source h5ad is intact (verified readable).

Strategy (per user choice = split & parallel):
    1. Split PD_spleen.h5ad (83,783 cells) into N cell-shards (same obs & var).
    2. Tokenize each shard in a SEPARATE worker subprocess -> shard{}.dataset.
    3. Concatenate shard datasets into PD_spleen.dataset.
    4. Verify load_from_disk returns 83,783 rows (matches n_cells).

Run:
    env ADPD_TISSUE=PD_spleen GENEFORMER_DIR=... \
        .venv/bin/python analysis/08_rebuild_tokenized_pd_spleen.py

Env knobs:
    REBUILD_SHARDS      number of shards (default 8)
    REBUILD_NPROC_PER   tokenizer nproc per worker (default 1; shards give the
                        parallelism, keep 1 per worker to avoid nesting)
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

os.environ["WANDB_DISABLED"] = "true"

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _resolve_tissue import find_h5ad, resolve

ROOT, PREFIX = resolve()
N_SHARDS = int(os.environ.get("REBUILD_SHARDS", "8"))
NPROC_PER = int(os.environ.get("REBUILD_NPROC", "1"))

H5AD = find_h5ad(ROOT)
WORK = ROOT / "tokenize_work"
SHARD_DIR = WORK / "shards_h5ad"
DS_DIR = WORK / "shards_dataset"
SHARD_DIR.mkdir(parents=True, exist_ok=True)
DS_DIR.mkdir(parents=True, exist_ok=True)
TOKENIZED = ROOT / "tokenized" / f"{PREFIX}.dataset"

ATTR_COLUMNS = ["cell_id", "individual", "celltype", "split",
                "disease", "samples4", "orig.ident"]


def split_shards() -> list[Path]:
    """Write N cell-shard h5ad files preserving obs + ENSG var index."""
    import anndata as ad

    adata = ad.read_h5ad(str(H5AD), backed="r")
    n = adata.n_obs
    paths = []
    for i in range(N_SHARDS):
        lo = n * i // N_SHARDS
        hi = n * (i + 1) // N_SHARDS
        if lo == hi:
            continue
        sub = adata[lo:hi].to_memory().copy()
        p = SHARD_DIR / f"shard{i:02d}.h5ad"
        sub.write_h5ad(str(p))
        paths.append(p)
        print(f"  wrote {p.name}: {sub.n_obs} cells", flush=True)
    print(f"split {n} cells into {len(paths)} shards", flush=True)
    return paths


def worker(shard_path: Path) -> None:
    """Tokenize one shard h5ad -> shard{nn}.dataset (run in a subprocess)."""
    from geneformer import TranscriptomeTokenizer

    tk = TranscriptomeTokenizer(
        custom_attr_name_dict={c: c for c in ATTR_COLUMNS},
        nproc=NPROC_PER,
        chunk_size=512,
        model_version="V2",
        use_h5ad_index=True,
    )
    tk.tokenize_data(
        data_directory=str(SHARD_DIR),
        output_directory=str(DS_DIR),
        output_prefix=shard_path.stem,
        file_format="h5ad",
        input_identifier=shard_path.stem,
    )
    print(f"done {shard_path.stem}", flush=True)


def run_parallel(shards: list[Path]) -> None:
    """Launch one worker subprocess per shard (bounded concurrency = N_SHARDS)."""
    proc = sys.executable
    script = os.path.abspath(__file__)
    procs = []
    for s in shards:
        env = dict(os.environ)
        env["REBUILD_WORKER_SHARD"] = s.name
        procs.append(subprocess.Popen(
            [proc, script, "--worker"], env=env, cwd=str(Path(__file__).parent),
        ))
    for p in procs:
        rc = p.wait()
        if rc != 0:
            raise SystemExit(f"worker failed rc={rc}")


def assemble_and_verify() -> None:
    from datasets import load_from_disk

    ds_paths = sorted(DS_DIR.glob("shard[0-9][0-9].dataset"))
    if not ds_paths:
        raise SystemExit("no shard datasets produced")
    print(f" concatenating {len(ds_paths)} shard datasets...", flush=True)
    from datasets import concatenate_datasets, load_from_disk

    shards = [load_from_disk(str(p)) for p in ds_paths]
    merged = concatenate_datasets(shards)
    TOKENIZED.parent.mkdir(parents=True, exist_ok=True)
    # remove the corrupt target dir if present
    if TOKENIZED.exists():
        import shutil

        shutil.rmtree(str(TOKENIZED))
    merged.save_to_disk(str(TOKENIZED))
    print("saved:", TOKENIZED, flush=True)

    # verification: re-load through the exact path ISP uses
    check = load_from_disk(str(TOKENIZED))
    print(f"VERIFY: {len(check)} rows", flush=True)
    from collections import Counter

    c = Counter(check["celltype"])
    print("celltype counts:", dict(c), flush=True)
    print("columns:", list(check.features), flush=True)
    n_exp = 83783
    if len(check) != n_exp:
        raise SystemExit(f"VERIFY FAIL: got {len(check)} rows, expected {n_exp}")
    print("=== TOKENIZE REBUILD OK ===", flush=True)


if __name__ == "__main__":
    worker_shard = os.environ.get("REBUILD_WORKER_SHARD")
    if worker_shard:
        worker(SHARD_DIR / worker_shard)
        sys.exit(0)

    print(f"rebuilding tokenized {PREFIX} (n_shards={N_SHARDS})", flush=True)
    shards = split_shards()
    run_parallel(shards)
    assemble_and_verify()