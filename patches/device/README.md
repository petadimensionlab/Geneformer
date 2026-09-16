# device 修正（CUDA 直書きの漏れ）— patches

> `patches/geneformer_multibackend.patch` は device 抽象化を入れるパッチですが、
> **`evaluation_utils.py` と `in_silico_perturber.py` の 6 箇所を漏らしています。**
> CUDA の無い環境（macOS / MPS）では、この 6 箇所が実行時に例外になります。

## Why — 実測した症状

### 1. fine-tune の評価ステップが必ず落ちる

`06_finetune.py` は学習完走後に `classifier.validate()` の中で評価します。そこが
`evaluation_utils.py` の CUDA 直書きに当たります:

```
File "geneformer/classifier.py", line 837, in validate
    result = self.evaluate_model(
File "geneformer/classifier.py", line 1343, in evaluate_model
    y_pred, y_true, logits_list, predict_metadata_all = eu.classifier_predict(
File "geneformer/evaluation_utils.py", line 130, in classifier_predict
    input_ids=input_data_batch.to("cuda"),
AssertionError: Torch not compiled with CUDA enabled
```

実害: **学習は 600 ステップ完走してチェックポイントも保存済みなのに、評価で落ちて
`adpd_finetuned_summary.json` / `adpd_model_comparison.csv` / 混同行列が一切書かれない。**
（本ワークスペースでは V2-316M の 2 時間 2 分の学習のあと、ここで停止しました。）

### 2. in silico perturbation も落ちる

ISP は `cell_states_to_model` を指定した時点で `mean_nonpadding_embs` に
`device="cuda"` のテンソルを渡します。つまり **state 埋め込みを計算する最初の段階で
同じ AssertionError** になります。評価だけ直しても ISP では確実に落ちます。

## What — 変更は 6 箇所 + import 2 行だけ

| ファイル | 行（変更前） | 変更 |
|---|---|---|
| `evaluation_utils.py` | 130–132 | `input_data_batch/attn_msk_batch/label_batch.to("cuda")` → `.to(get_device())` |
| `evaluation_utils.py` | (import) | `from .device import get_device` を追加 |
| `in_silico_perturber.py` | 688 | `torch.tensor(minibatch["length"], device="cuda")` → `device=get_device()` |
| `in_silico_perturber.py` | 693 | `torch.tensor(perturbation_batch["length"], device="cuda")` → `device=get_device()` |
| `in_silico_perturber.py` | 748 | `torch.tensor(nonpadding_lens, device="cuda")` → `device=get_device()` |
| `in_silico_perturber.py` | (import) | `from .device import get_device` を追加 |

`get_device()` は `geneformer/device.py`（`patches/device.py`）が返す
`"cuda" | "mps" | "cpu"` です。CUDA 環境では従来と同一の挙動になります。

## 適用方法

```bash
cd geneformer_hf
# 方法 A: patch で当てる（推奨・レビュー可能・行番号付き）
git apply -p1 ../patches/device/geneformer_cuda_device.patch
# 方法 B: 用意済みファイルをコピー（download.sh 後の復旧に便利）
cp ../patches/device/evaluation_utils.py      geneformer/
cp ../patches/device/in_silico_perturber.py   geneformer/
```

> ⚠️ `./download.sh` は `geneformer/*.py` を取り直すため、**この修正は毎回外れます。**
> download.sh 実行後は `patches/bf16/`, `patches/checkpoints/`, `patches/device/` を
> まとめて再適用してください。

## 検証した内容

- patch のベースラインは **HF から `download.sh` が取得する素のファイル**
  （このリポジトリのパッチ前の状態）にしてあります
- `patch -p1 --dry-run` → **exit 0**、`git apply --check` → **クリーン**
- 生成した patch を素のコピーに適用した結果が、稼働中の `geneformer_hf/` の
  該当行と一致することを確認済み
- この修正後、`06_finetune.py` の評価が完走し
  `adpd_finetuned_summary.json`（accuracy / macro F1）が出力されることを実測

## 残っている既知のギャップ（未修正）

`geneformer/mtl/*.py` には

```python
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
```

が 3 箇所残っています。これは例外にはならず **CPU にフォールバックする**だけなので、
今回の経路（`analysis/` の ISP / fine-tune）では影響しません。ただし MTL を MPS で
使うと遅くなります。MTL を使う場合は `get_device_obj()` に置き換えてください
（未検証のため、このパッチには含めていません）。
