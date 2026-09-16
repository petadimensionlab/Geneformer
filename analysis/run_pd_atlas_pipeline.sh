#!/usr/bin/env bash
# Full PD multi-region atlas pipeline (Geneformer V2-316M, bf16, MPS).
#
#   [0] 04b_extract_embeddings.py  — tokenize + frozen cell embeddings (started separately;
#                                    this script waits for it)
#   [1] 04_baseline.py             — frozen-embedding logistic probe (reuses the CSV)
#   [2] 06_finetune.py             — fine-tuned cell-type classifier
#   [3] 07g_pd_atlas_perturbation.py — in silico gene deletion, PD -> normal shift
#
# All stages are idempotent: existing tokenized data / embeddings / checkpoints are reused.
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."

export GENEFORMER_DIR="$PWD/geneformer_hf"
export GENEFORMER_MODEL="${GENEFORMER_MODEL:-Geneformer-V2-316M}"
export ADPD_TISSUE=PD_atlas
export GF_DTYPE=bf16                 # perturber_utils.load_model -> bfloat16
export FINETUNE_BF16=1
# 64 GB host: batch 4 still exhausted swap (step time drifted 25 s -> 500 s),
# so batch 2 + a step cap is the workable setting on this machine.
export FINETUNE_BATCH_SIZE="${FINETUNE_BATCH_SIZE:-2}"
export FINETUNE_MAX_STEPS="${FINETUNE_MAX_STEPS:-600}"
export IS_NPROC=1
export IS_MAX_CELLS="${IS_MAX_CELLS:-300}"
export IS_FORWARD_BATCH=16

EMB="input/PD_atlas/results/embeddings/pretrained_cell_embeddings.csv"

echo "=== [0/3] waiting for the embedding extraction to finish ==="
while pgrep -f "04b_extract_embeddings" >/dev/null; do sleep 30; done
if [ ! -f "$EMB" ]; then
  echo "ERROR: $EMB was not produced — embedding stage failed"; exit 1
fi
echo "embeddings present: $(( $(wc -l < "$EMB") - 1 )) cells"

echo; echo "=== [1/3] frozen-embedding baseline probe (04_baseline.py) ==="
.venv/bin/python analysis/04_baseline.py || { echo "FAILED: 04_baseline"; exit 1; }

echo; echo "=== [2/3] fine-tune cell classifier (06_finetune.py) ==="
.venv/bin/python analysis/06_finetune.py || { echo "FAILED: 06_finetune"; exit 1; }

echo; echo "=== [3/3] in silico perturbation (07g_pd_atlas_perturbation.py) ==="
.venv/bin/python analysis/07g_pd_atlas_perturbation.py || { echo "FAILED: 07g"; exit 1; }

echo; echo "=== PD ATLAS PIPELINE COMPLETE ==="
