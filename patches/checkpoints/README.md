# Fine-tuning のチェックポイント（中断に強い学習）— patches

## Why

Geneformer の `Classifier` は Trainer の保存設定を **`save_strategy="epoch"` /
`save_total_limit=1` にハードコード**していました（旧 `classifier.py` の
`def_training_args.update({"save_strategy": "epoch", "save_total_limit": 1})`）。

そのため **エポックの途中で落ちる（OOM・クラッシュ・スリープ・強制終了）と、
それまでの全ステップが失われます**。実測でこれを踏みました:

- V2-316M・MPS・batch 4 の 1 エポック = **6 時間 52 分**
- 学習自体は完走したが、直後の評価ステップがデータ都合で例外終了
- `save_strategy="epoch"` だったため、**完了した 6.9 時間分の重みがディスクに無い**
  （同時に、後続の自動化が「失敗」と誤判定して `rm -rf` したため完全に消失）

## Fix

`classifier.py` の 2 つの学習経路（`train_classifier` のハイパラ探索経路と、
`validate` 経由の通常経路）の両方で:

1. **保存既定値をユーザー指定で上書きできるようにする**（`update` → `setdefault`）
2. 既定を **`save_strategy="steps"` / `save_total_limit=2`** にする
   （`save_steps` は `logging_steps` から自動決定、`FINETUNE_SAVE_STEPS` で上書き可）
3. `load_best_model_at_end=True` のときは HF の要求どおり
   `eval_strategy` も `"steps"` に揃える（戦略不一致で例外になるため）
4. `trainer.train(resume_from_checkpoint=...)` で
   **最新の `checkpoint-*` から自動再開**する（`transformers.trainer_utils.get_last_checkpoint`）

これにより、中断しても最後のチェックポイント以降だけをやり直せば済みます。

## 適用方法

```bash
cp patches/checkpoints/classifier.py geneformer_hf/geneformer/classifier.py
```

## 注意: `download.sh` はパッチを上書きする

`./download.sh --model V2-316M` は `geneformer/*.py` も取得し直すため、
**適用済みのローカルパッチ（`patches/bf16` の bf16 修正や multi-backend の device 修正）が
元に戻ります**。実際にこれを踏み、BF16 の ISP が
`TypeError: Got unsupported ScalarType BFloat16` で落ちました。

fresh checkout の後は必ず次を実行してください:

```bash
cp patches/bf16/emb_extractor.py    geneformer_hf/geneformer/emb_extractor.py
cp patches/bf16/perturber_utils.py  geneformer_hf/geneformer/perturber_utils.py
cp patches/checkpoints/classifier.py geneformer_hf/geneformer/classifier.py
```

## 検証方法

```bash
# 短い run で checkpoint-* が step 単位で出ることを確認
FINETUNE_SAVE_STEPS=20 FINETUNE_BF16=1 GENEFORMER_MODEL=Geneformer-V2-316M \
ADPD_TISSUE=<TISSUE> .venv/bin/python analysis/06_finetune.py
ls input/<TISSUE>/runs/*cellClassifier*/ksplit1/   # checkpoint-20/, checkpoint-40/, ... が出る

# 中断 → 再開: 途中で Ctrl-C して再実行すると "Resuming fine-tuning from checkpoint" が出る
```
