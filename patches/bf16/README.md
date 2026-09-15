# bf16 support for embedding extraction / ISP (patches)

## Why

`GF_DTYPE=bf16` (see `analysis/_isp_common.py`) casts every model that
`perturber_utils.load_model()` returns to bfloat16. bf16 is the recommended
inference dtype on this hardware: ~5x the fp32 throughput with an embedding
cosine of 0.99995-0.99998 (docs/quantization/EXPERIMENT.md).

Upstream geneformer then fails, because the embedding export paths hand a
bfloat16 tensor straight to numpy, and numpy has no bfloat16 dtype:

```
TypeError: Got unsupported ScalarType BFloat16
  geneformer/emb_extractor.py:256  pd.DataFrame(embs.cpu().numpy())
  geneformer/emb_extractor.py:283  .cpu().numpy()
  geneformer/emb_extractor.py:702  pd.DataFrame(embs.cpu().numpy()).T
  geneformer/perturber_utils.py:824 torch.cat(tensor_list).cpu().numpy()
```

## Fix

Insert `.float()` before `.cpu().numpy()` at the four sites above. On an fp32
tensor this is a no-op; on bf16 it upcasts. Embeddings are exported as float32
either way, so results for fp32 runs are unchanged bit-for-bit.

Applied in place here, since `geneformer_hf/` is a clone with local patches
(see `patches/geneformer_multibackend.patch`). Ready-made copies of the two
touched files live next to this note:

```bash
cp patches/bf16/emb_extractor.py    geneformer_hf/geneformer/emb_extractor.py
cp patches/bf16/perturber_utils.py  geneformer_hf/geneformer/perturber_utils.py
```

Note: `patches/bf16/*.py` also contain the earlier multi-backend device
patches, because those files carry both changes. Copying them over a fresh
`download.sh` checkout therefore applies the device patches too — which is
what the workspace does anyway.

## Verifying

```bash
# 1 gene, 1 timepoint smoke test on the spleen classifier
GENEFORMER_DIR=geneformer_hf ADPD_TISSUE=AD_spleen GF_DTYPE=bf16 \
IS_MAX_GENES=1 IS_TIMEPOINTS=3m ISP_EXPERIMENT=smoke_bf16 \
.venv/bin/python analysis/07_ad_spleen_early_isp.py
```

Without the fix this dies with the TypeError above; with it the run produces
`*_early_isp_stats_combined.csv` normally.
