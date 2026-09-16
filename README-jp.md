# Geneformer V2-104M — マルチバックエンド (MPS / CUDA / ROCm) 対応

> [English](README.md) · **日本語**

このワークスペースでは、[`ctheodoris/Geneformer`](https://huggingface.co/ctheodoris/Geneformer) の
**Geneformer V2-104M** を、3 つのアクセラレータで動作させます:

| バックエンド | ハードウェア | 備考 |
|---|---|---|
| **MPS** | Apple Silicon (Mac) | 本ワークスペースで検証済み |
| **CUDA** | NVIDIA GPU | 元来の対応パス |
| **ROCm** | AMD GPU (WSL2 経由) | 実装済み。この Mac では未検証 |

## デバイス抽象化

中央モジュール `geneformer/device.py` が計算デバイスを一度だけ解決し、
ハードコードされた `"cuda"` 呼び出しをすべてここに通します:

- **自動検出の優先度**: CUDA/ROCm (`torch.cuda.is_available()`) > MPS > CPU。
- **上書き**: `GENEFORMER_DEVICE=mps|cuda|cpu` でバックエンドを強制。

```python
from geneformer.device import get_device, get_device_obj, move_to_device, empty_cache
get_device()          # -> "mps" | "cuda" | "cpu"
get_device_obj()      # -> torch.device("mps") ...
```

以前 `"cuda"` を直書きしていた呼び出し箇所すべて(emb_extractor、
perturber_utils、evaluation_utils、in_silico_perturber、classifier、
mtl/{train,eval_utils,utils}、`analysis/` のスクリプト)は
`get_device_obj()` / `move_to_device()` / `empty_cache()` 経由に変更済みです。

### ROCm が CUDA と等しい理由

ROCm は **CUDA API を公開**しています。AMD 用 ROCm ビルドの PyTorch は
`torch.cuda.is_available() == True` を返し、`torch.cuda.*` 呼び出しが AMD GPU
を駆動します。つまり、同じ `device="cuda"` コードパスが AMD ROCm ハードウェア
をそのまま駆動するため、別ブランチは不要です。したがって torch が ROCm ビルド
なら、この抽象化は ROCm を「そのまま」サポートします。

## MPS 検証(この Mac)

Apple Silicon 上でエンドツーエンド実行:

1. **トークナイゼーション** `input/PD_smallint/PD_smallint.h5ad`
   (64,614 cells、25 cell types)→ `analysis_ws/tokenized/PD_smallint.dataset`
   — **動作**。
- 必要だった修正: `use_h5ad_index=True`(var index がすでに ENSG ID を持つ)、
      および var index が文字列 Series の場合の `tokenize_anndata` 位置インデックス。
      upstream のコードは pandas の `var["ensembl_id_collapsed"]` Series に整数位置を
      `[...]` で渡すため、gene 名 index に対する**ラベル検索**になり、誤った遺伝子を
      選びます。`patches/geneformer_multibackend.patch` は `tokenize_anndata` 内の
      `norm_factor_vector` と `coding_miRNA_ids` で `.iloc[coding_miRNA_loc]`(位置検索)
      を使うよう修正します。
- Linux では、入力ディレクトリに **隠しファイル**(macOS の AppleDouble メタ
    ファイル `._*.h5ad` や `.DS_Store`。SMB/NFS や Mac→Linux コピーで混入)があると、
     `OSError: Unable to synchronously open file (file signature not found)`
     でトークナイズが失敗することがあります。patch は `tokenize_files` で
     `.` で始まるファイルをスキップし、`04_baseline.py` も glob 時に
     隠し `.h5ad` を除外します。
2. **Frozen-embedding + probe**(`analysis/04_baseline.py`)— MPS で**動作**。
   - `analysis/smoke_mps.py` が MPS vs CPU の出力一致を ~1.9e-5 で検証。
3. **Fine-tuning**(`analysis/06_finetune.py`)— MPS で**動作**(HF `Trainer` 使用。
   こちらは MPS をサポート済み)。
   - **メモリ予算(M4 Max・ユニファイドメモリ 128 GB)**: collator はバッチ内の
     最長セルに合わせてパディングするため(PD_BM: 中央値 1201 トークン、最大
     4096)、最悪長の 1 train step は素の backward で ~145 GiB に達し、数ステップ
     で MPS OOM(`max allowed: 182.78 GiB`)になります。`analysis/measure_mps_batch.py`
     がバッチサイズごとに実測します:

     | 設定(L=4096 最悪ケース)              | MPS ピークメモリ |
     |---------------------------------------|------------------|
     | batch 8、checkpointing 無し           | ~145 GiB(OOM)   |
     | batch 4、checkpointing 無し           | ~68 GiB          |
     | **batch 8 + gradient checkpointing**  | **~64 GiB**      |

     `06_finetune.py` への対応:`per_device_train_batch_size: 8` は維持しつつ
     `gradient_checkpointing: True`(`use_reentrant: False`)を有効化 —
     実効バッチ / オプティマイザ挙動は変えずに活性化メモリを約 1/2.3 化。
     スワップ解消もあり ~6 s/it(従来比 ~4 倍高速)で、1 epoch は数時間で完了。
   - **不完全チェックポイントの復旧**: 以前の fine-tune が`trainer.save_model()`
     が重みを書き出す前に中断された場合(MPS OOM など)、
     `<TISSUE>/runs/<prefix>.../ksplitN/` ディレクトリは**空のまま**残ることがあります。
     `06_finetune.py` は実際に重み(`pytorch_model.bin` / `model.safetensors` /
     `adapter_model.bin`)を含む場合のみチェックポイントを再利用可能とみなし、
     空・不完全な `ksplit*` ディレクトリは自動削除して、モデルファイル欠落エラーになる
     代わりに新規 fine-tune を実行します。
   - **MPS の INT_MAX 制約（V2-316M で顕在化）** — MPS の attention は
     `batch × heads × seq²` 要素のスコアテンソルを確保します。316M は 18 ヘッドなので
     系列長 4096 では **batch 8 で 2.42e9 要素**となり `INT_MAX (2^31 = 2.147e9)` を
     超えて `RuntimeError: MPSGraph does not support tensor dims larger than INT_MAX`
     で**最初の backward で即落ち**します（104M は 12 ヘッドで 1.61e9 のため batch 8 でも通る）。
     対策は次のいずれかです。
     - **batch を 4 以下にする**（1.21e9。実測で安定。本ワークスペースの 316M 解析はこれで実行）
     - 系列長を切り詰める（batch 8 なら 3072 で 1.36e9、2048 で 6.0e8）

     `gradient_checkpointing` は活性化メモリを減らしますが、この制約には効きません。
     同じ上限は ISP と `EmbExtractor` の `forward_batch_size` にも当てはまるので、
     316M を MPS で回すときは `forward_batch_size` も 4 に落としてください。
   - **チェックポイント（中断に強い学習）** — 長時間の fine-tune の前に
     `cp patches/checkpoints/classifier.py geneformer_hf/geneformer/classifier.py`
     を適用してください。上流は `save_strategy="epoch"` / `save_total_limit=1` を
     ハードコードしていたため、**エポック途中で落ちると完了済みの全ステップが失われます**
     （V2-316M の MPS 1 エポックは約 7 時間）。パッチ版の `Classifier` は
     ステップ単位のチェックポイントを既定にし（ユーザー指定が優先）、
     `trainer.train(resume_from_checkpoint=...)` で最新の `checkpoint-*` から自動再開します。
     詳細: `patches/checkpoints/README.md`。
   - **`download.sh` はローカルパッチを打ち消します** — `geneformer/*.py` を取り直すため、
     fresh checkout の直後は bf16 / device パッチが外れています（実際にこれで
     `GF_DTYPE=bf16` の ISP が `TypeError: Got unsupported ScalarType BFloat16` で落ちました）。
     `./download.sh` の後は必ず再適用してください:
     `cp patches/bf16/emb_extractor.py patches/bf16/perturber_utils.py
     patches/checkpoints/classifier.py patches/device/evaluation_utils.py
     patches/device/in_silico_perturber.py geneformer_hf/geneformer/`。
     （`patches/device/` は CUDA 直書きの漏れ 6 箇所の修正。当てないと **fine-tune の
     評価と ISP が `AssertionError: Torch not compiled with CUDA enabled` で落ちます** —
     学習は完走するのに指標が一切書かれないので気づきにくい。詳細は
     [`patches/device/README.md`](patches/device/README.md)）
4. **In silico perturbation**(`analysis/07_in_silico_perturbation.py`、チュートリアル
   ノートブック移植)— MPS で**動作**。

## In silico perturbation — nproc / datasets / max_ncells(苦労して得た知見)

`InSilicoPerturber.perturb_data` は `dataset.map(make_group_perturbation_batch,
num_proc=nproc)` で perturbation `Dataset` を構築します。3 つの別々の落とし穴に
遭遇し解決済みで、`analysis/07_in_silico_perturbation.py` はすべてを適用し、
実行時に warning / RAM 見積もりを出力します。

1. **`datasets>=5` で `perturb_data` がハングする** — `dataset.map` ステップが
   *どの nproc でも*(`nproc=8` と `nproc=1` の両方で確認)CPU 0% のまま戻らない。
   **根本対策: `datasets==4.0.0` に固定。** nproc での回避は症状を隠すだけ。
   `07` は `datasets>=5` を検出すると大きな warning を出力します。

2. **`nproc`** — `datasets==4.x` では map 関数(ネストしたクロージャで
   `self.tokens_to_perturb` 等を参照)は `spawn` で pickle でき、nproc>1 も動作可。
   **ここで実証済み:** `make_group_perturbation_batch` と同じクロージャ形状を
   `datasets==4.0.0` で `num_proc=1/2/4` に渡して `dataset.map` を実行し、
   3 つとも同一の正しい結果を返しました。つまり datasets 4.x なら nproc>1 は
   機能します。`07` はそれでもデフォルト **`nproc=1`**(検証済み安定)、
   `IS_NPROC` で上書き可。nproc=1 を推奨する理由:
   - ワーカーが RAM を消費し、同一ホストの MPS forward pass と競合。
   - 2 遺伝子以上(または `combos>0`)だとセルあたりバリアント数が増え、map 出力が
     大きくなり multiprocessing で型が不安定になりやすい — nproc=1 はそれを回避。
   - **並列化が効くのは大規模データのみ。** 4 セルのベンチでは nproc=1 が ~0.01s、
     nproc=2/4 が ~0.09s(spawn + プロセス間通信のオーバーヘッドが小規模で支配的)。
     単一遺伝子の perturbation では map は直列でも十分速く、MPS forward が常に支配的。
   - `07` は **`nproc>1` のとき必ず warning** を出し、datasets<5 と控えめな
     `max_ncells` が必須であることを再通知。

   具体的に **2 遺伝子以上**の場合も nproc=1 が正解です。遺伝子リストを 2 倍にすると
   セルあたりのバリアントが概ね倍増(削除モード: 1 遺伝子削除→セルあたり 1 バリアント、
   N 遺伝子削除 `combos=0`→最大 N バリアント/セル)するため、メモリは遺伝子数に
   比例し、forward pass は単一 MPS プロセスのまま。map の並列化は支配的な MPS
   forward コストには効きません。

   遺伝子を増やしつつ速く保つには: `max_ncells` を小さく(次項)、`combos=0`、
   大きなリストより生物学的に狙った少数遺伝子。

4. **`max_ncells` / OOM** — `perturb_data` は perturbation データセット全体を
   RAM にマテリアライズ(~ `n_cells * n_variants * seq_len`)。2000 セルの実行で
   この MPS ホストが OOM しました(スワップ約28–32 GB 使用)。`07` はデフォルト
   `max_ncells=200`(`IS_MAX_CELLS` で上書き)、perturb 前に
   `estimate_perturb_ram()` による RAM 見積もりを表示します。目安:
   `バイト数 ≈ n_cells * n_variants * seq_len * 10`、空き RAM より十分小さく。
   ```
   n_variants = max(n_genes,1)                      # combos=0
   n_variants = C(n_genes, combos+1)                # combos>0(二項係数で爆発)
   ```
   単一遺伝子 `combos=0` は安全。2 遺伝子以上 × `combos>0` は急速に爆発します。

5. **stats 側の `genes_perturbed` も絞る** — `InSilicoPerturberStats(...、
   genes_perturbed="all")` は語彙の全遺伝子を再走査し、MPS では極端に遅い。
   実際に perturb した遺伝子に一致させる必要があります。`07` は
   `genes_perturbed=genes_to_perturb` に設定。

`07` が使う環境変数: `IS_NPROC`(既定 1)、`IS_MAX_CELLS`(既制 200)。

> **複数臓器の早期疾患スクリーニング** — 複数臓器(小腸/脳/血液/脾臓/肝臓/骨髄)・
> 複数ホストで実施した AD/PD の in silico perturbation は、
> **[In Silico Perturbation Wiki](docs/isp/README.md)** にまとめています
> (出力先の統一規則・共通 `_isp_common.py`・各臓器の設定と結果)。

## モデルの量子化 — 実測した知見

<a id="quant-ja"></a>
**日本語** · [English](README.md#quant-en)

**推論・埋め込み抽出・in silico perturbation・微調整**のそれぞれで、
bf16 / bitsandbytes の int8 / 4bit(nf4) がどう違うか。さらに、モデルサイズ
（V2-104M と V2-316M）と規模（batch 8〜256、系列長 512〜4,096）を
変えるとどうなるかを実測しました。

測定値・方法・落とし穴の一覧は
**[EXPERIMENT-jp.md](docs/quantization/EXPERIMENT-jp.md)**（日本語）/
**[EXPERIMENT.md](docs/quantization/EXPERIMENT.md)**（English）にまとめています。

GB10（統合メモリ 119.63 GiB）での主な結果:

- **bf16 は fp32 の約5倍速く、精度の低下は測定できないレベル**です
  （埋め込み cos 0.99998）。埋め込み抽出・ISP・微調整のすべてで最適解です。
- **int8 は安全ですが bf16 より遅い**です（fp32 の1.6〜1.7倍）。
  bf16 が使えない環境と、配布物のサイズ削減のときだけ使います。
- **4bit(nf4) は遺伝子の順位を壊します**（top-100 の一致が 87/100）。
  perturbation の成果物が遺伝子ランキングである以上、使えません。
- **セル数や遺伝子数を増やしても int4/int8 の利点は出ません** —
  batch 8〜256 で速さはほぼ一定、batch 256 での bf16 と nf4 のメモリ差は 0.12 GiB。
  int4/int8 が有利なのは配布物のサイズと CPU 推論だけです。
- **微調整では量子化は逆効果**です。bf16 + 勾配チェックポイントが
  QLoRA よりメモリも速度も有利です（104M・batch 8・系列長 4096 で 1.82 GiB 対 9.83 GiB）。
- **事前量子化は配布用には有効ですが、速度のためには損**です —
  推論は同一（埋め込みはビット単位で一致）なのにロードは 2〜6倍遅くなります。
  bf16 で事前保存したチェックポイントは `torch_dtype` を明示しないと
  fp32 として読み込まれます。
- **V2-316M の分類器を作って AD_spleen で端から端まで実測**しました —
  同一条件ペア（24遺伝子 × 3時点）で bf16 は ISP の所要時間を **2.28倍**、
  GPU エネルギーを **3.82倍**削減（fp32 40分57秒 / 47.7 Wh 対 **bf16 17分58秒 /
  12.5 Wh**）、そのうえで **top-20 の遺伝子一致は 20/20**。55遺伝子の canary は 46分、
  微調整は 2時間・137 Wh でした。詳細は
  [REPORT-316M-ISP-jp.md](docs/quantization/REPORT-316M-ISP-jp.md)。
- **316M 分類器は細胞種分類を改善しません。** AD_spleen では 104M と引き分けです
  （26,116 held-out セルで accuracy 0.9320 対 0.9312、macro F1 0.8890 対 0.8867）。
  一方で ISP の遺伝子順位は大きく変わるため、大きいモデルへの移行は生物学的な
  妥当性で正当化が必要です。**bf16 にはその必要がありません。**

| ドキュメント | 言語 | 内容 |
|---|---|---|
| [EXPERIMENT-jp.md](docs/quantization/EXPERIMENT-jp.md) | JP | **最初にこれ。** 実験のまとめ（全測定値・方法・段階ごとの指針・落とし穴一覧） |
| [REPORT-316M-ISP-jp.md](docs/quantization/REPORT-316M-ISP-jp.md) | JP | V2-316M 分類器 + AD_spleen の ISP canary、fp32 対 bf16 の順位判定、時間・エネルギー・負荷の比較 |
| [EXPERIMENT.md](docs/quantization/EXPERIMENT.md) | EN | 同上（英語版） |
| [PLAN.md](docs/quantization/PLAN.md) | JP | 量子化の全体計画（目的の整理・フェーズ・合格基準） |
| [STAGES.md](docs/quantization/STAGES.md) | JP | 段階ごとの影響（tokenization / embedding / ISP / fine-tuning） |
| [PREQUANT-VS-LOADTIME.md](docs/quantization/PREQUANT-VS-LOADTIME.md) | JP | 事前量子化とロード時量子化の比較 |
| [PLAN-316M-128GB.md](docs/quantization/PLAN-316M-128GB.md) | JP | V2-316M を 128GB のマシンで動かす計画 |
| [REPORT-104M-48GB.md](docs/quantization/REPORT-104M-48GB.md) | JP | 104M は 48GB の GPU で動くかの検証レポート |
| [docs/quantization/README.md](docs/quantization/README.md) | EN/JP | 上記の索引（バイリンガル） |

測定スクリプトは `analysis/10*_quant*.py`、`analysis/10e_train_mem_probe.py`、
`analysis/11*_prequant*.py` / `analysis/11b_load_time_bench.py` にあります。
JSON の生データと事前量子化アーティファクトは `quantized/` に出力されます（git 管理外）。

なお **`Geneformer-V2-316M` の重みは既定では入っていません**
（`model.safetensors` は 135 バイトの LFS ポインタ）。316M を使う前に
`./download.sh --model V2-316M` で実体を取得してください。

## bf16 で in silico perturbation を回す

<a id="isp-bf16-ja"></a>
**日本語** · [English](README.md#isp-bf16-en)

`GF_DTYPE=bf16` を付けると、geneformer が読み込むモデル（`perturber_utils.load_model`）が
すべて bfloat16 で動きます。AD_spleen の実測では **所要時間 2.28倍・GPUエネルギー 3.82倍削減・
平均電力 40%減・最高温度 7°C 低下**で、**上位20遺伝子は fp32 と 20/20 一致**しました
（[docs/quantization/REPORT-316M-ISP-jp.md](docs/quantization/REPORT-316M-ISP-jp.md)）。

### 1. 初回だけ必要な準備

numpy は bfloat16 を扱えないため、パッチを当てないと埋め込みの書き出しで
`TypeError: Got unsupported ScalarType BFloat16` で落ちます。

```bash
cp patches/bf16/emb_extractor.py    geneformer_hf/geneformer/emb_extractor.py
cp patches/bf16/perturber_utils.py  geneformer_hf/geneformer/perturber_utils.py
```

### 2. 実行

通常の ISP との違いは `GF_DTYPE=bf16` の1行だけです。

```bash
GENEFORMER_DIR=geneformer_hf \
GENEFORMER_MODEL=Geneformer-V2-316M \
ADPD_TISSUE=AD_spleen \
GF_DTYPE=bf16 \
ISP_EXPERIMENT=ad_spleen_316m_bf16 \
IS_CELLCLASSIFIER_DIR=$PWD/input/AD_spleen/runs/260915_geneformer_cellClassifier_AD_spleen_celltype_Geneformer-V2-316M/ksplit1 \
.venv/bin/python analysis/07_ad_spleen_early_isp.py
```

組織によってスクリプト名だけが変わります（環境変数の並びは同じ）:
AD_blood → `07_ad_blood_early_isp.py`、AD_brain → `07_ad_brain_early_isp.py`、
AD_smallint → `07_ad_smallint_early_isp.py`、PD_spleen → `07_pd_spleen_early_isp.py`、
AD_liver → `07e_ad_liver_perturbation.py`。

### 3. 必ず設定する4つ

| 変数 | 意味 | 注意 |
|---|---|---|
| `GF_DTYPE=bf16` | モデルを bfloat16 で動かす | 未設定 / `none` なら従来の fp32 |
| `ISP_EXPERIMENT=<名前>` | 出力フォルダ名 | **付けないと既存の fp32 結果を上書きします** |
| `IS_CELLCLASSIFIER_DIR=<path>` | 使う微調整済み分類器 | **明示推奨** — 自動解決は最新の `runs/**/ksplit*` を選ぶため、104M 実行でも 316M 分類器を掴みます |
| `GENEFORMER_MODEL` | 事前学習モデル | `Geneformer-V2-316M` または `Geneformer-V2-104M` |

分類器のパスは `ls -d input/<TISSUE>/runs/*/ksplit1` で確認できます。

### 4. 先にスモークテスト（1〜2分）

```bash
IS_MAX_GENES=1 IS_TIMEPOINTS=3m IS_MAX_CELLS=50 \
GENEFORMER_DIR=geneformer_hf GENEFORMER_MODEL=Geneformer-V2-316M ADPD_TISSUE=AD_spleen \
GF_DTYPE=bf16 ISP_EXPERIMENT=smoke_bf16 \
IS_CELLCLASSIFIER_DIR=$PWD/input/AD_spleen/runs/<run>/ksplit1 \
.venv/bin/python analysis/07_ad_spleen_early_isp.py
```

### 5. bf16 が効いているか確認する

ログに次の行が必ず出ます。出ていなければ fp32 のままで、結果は変わりません。

```
[config] GF_DTYPE -> bfloat16 (geneformer load_model)
```

### 6. `GF_DTYPE` が効くスクリプト

対応（`_isp_common` を読み込んでいるもの）: `07_ad_spleen_early_isp.py`、
`07_ad_blood_early_isp.py`、`07_ad_brain_early_isp.py`、
`07_ad_smallint_early_isp.py`、`07_ad_ln_early_isp.py`、
`07_pd_spleen_early_isp.py`、`07e_ad_liver_perturbation.py`、
`07f_in_silico_perturbation_AD_BM.py`。

古い系統（`07_in_silico_perturbation.py`、`07b_*`、`07c_*`、`07d_*`、
`07e_in_silico_perturbation_PD_smallint.py`）は無視します。
bf16 を使いたい場合は、import の後に2行足してください。

```python
from _isp_common import apply_dtype_override
apply_dtype_override()
```

### 7. 結果の検証と計測

```bash
# fp32 の結果との順位一致（スピアマン・top-20・符号一致）
.venv/bin/python analysis/14_compare_isp.py \
  --a input/AD_spleen/results/isp/<fp32の実験名> \
  --b input/AD_spleen/results/isp/<bf16の実験名> \
  --label fp32_vs_bf16 --out docs/quantization/g5

# 所要時間・電力・エネルギー・温度＋ログからの遺伝子別所要時間
.venv/bin/python analysis/13_profile.py --label isp-bf16 \
  --out docs/quantization/profiles -- .venv/bin/python analysis/07_ad_spleen_early_isp.py
.venv/bin/python analysis/15_isp_timing.py docs/quantization/profiles/isp-bf16.log
```

### 8. つまずきやすい点

1. **量子化（`model_type="Pretrained-Quantized"`）と `GF_DTYPE=bf16` は併用できません。**
   量子化済みモデルには `.to(dtype)` を呼べないためです。bf16 を使うなら量子化は不要です。
2. **`IS_MAX_CELLS` を上げるほど bf16 の利得が伸びます。** 遺伝子×時点の1ユニットあたり
   約12秒は dtype と無関係な固定費（データセット構築・stats・pickle）で、そこは速くなりません。
3. **`datasets==4.0.0` と `IS_NPROC=1` を維持**してください（上の ISP セクション参照）。
4. 期待値: 所要時間 **2.28倍**、エネルギー **3.82倍**。上位20遺伝子は fp32 と同一で、
   差は `|Shift| ≈ 2e-04` 以下（ノイズフロア内）に収まります。

## Setup (uv)

```bash
# 1. Geneformer パッケージ + 重みを geneformer_hf/ に取得(curl のみ・
#    git / git-lfs 不要)、マルチバックエンド(MPS/CUDA/ROCm)の device パッチ適用。
./download.sh
cd geneformer_hf
git apply ../patches/geneformer_multibackend.patch   # cwd 非依存。どのユーザーでも動作
cp ../patches/device.py geneformer/device.py        # 新規ファイル(未追跡)
cd ..
# `git apply` が失敗する場合(例: tokenizer.py を手で編集済み)は、用意済みの
# 修正済みファイルを直接コピーして適用できます:
#   cp patches/tokenizer.py geneformer_hf/geneformer/tokenizer.py
#   cp patches/device.py    geneformer_hf/geneformer/device.py
# (tokenizer.py には .iloc 位置インデックス修正と隠しファイルスキップが含まれます)

# 2. 環境作成
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python -e geneformer_hf
# geneformer は transformers==4.46 を要求(setup.py は範囲指定のため、
# uv のリゾルバが 5.x を選ぶ可能性がある — 固定し直す):
uv pip install --python .venv/bin/python "transformers==4.46.3"
# datasets は必ず < 5。datasets>=5 だと InSilicoPerturber.perturb_data が
# dataset.map でハングする(任意 nproc):
uv pip install --python .venv/bin/python "datasets==4.0.0"
```

> **再実行時の注意:** `download.sh` は冪等で、既存の `geneformer_hf/` を
> 再利用します。一度 patch を適用済みなら、再度適用すると
> "patch does not apply"(適用済み)で失敗します — これは正常です。
> `geneformer/device.py` が既に存在する場合は patch 手順をスキップしてください。

> **固定バージョン:** `transformers==4.46.3`(5.x は `SpecialTokensMixin` が壊れる)
> と **`datasets==4.0.0`**(`datasets>=5` は `perturb_data` の `dataset.map` が
> ハング — 上記 IS perturbation の注を参照)。どちらも標準パイプラインを
> 確実に動かすために必須です。

モデル重みの管理は upstream では LFS ですが、`download.sh` は git を全く
使いません — 全ファイル(コード、重み、辞書)を `huggingface.co` から `curl` で
直接取得します:

```bash
./download.sh              # V2-104M + パッケージコード + 辞書
./download.sh --all        # 全モデル (V1-10M, V2-104M, CLcancer, V2-316M)
# 対象変更: ./download.sh --model V1-10M
```

git-lfs を使わず素の `curl` で取得するため、git-lfs のインストール有無に
かかわらずどのホストでも動きます — LFS smudge のハングや LFS ポインタ残存
(git lfs の典型的な失敗)、`safetensors: header too large` /
`pickle.UnpicklingError: invalid load key, 'v'` の読み込みエラーも起きません。
環境変数 `GENEFORMER_DIR`(出力先)、`HF_MIRROR` で上書き可。既存ファイルは
スキップされるので再実行も安全(冪等)。

**実行のたびに全ファイルの整合性チェック(`verify_all`)を行います。**
`setup.py` と `geneformer/*.py` のコード、`model.safetensors`、`config.json`、
`generation_config.json`、`training_args.bin`、`geneformer/*.pkl` の辞書が
**LFS ポインタではなく実データ**(最小サイズ以上)として存在するかを確認します。
`pyproject.toml` は upstream が未管理(Geneformer は `setup.py` のみ)のため、
`download.sh` がダウンロード後に最小構成を自動生成します — 手動の `cp` は
不要です。不足があれば非ゼロで失敗します — つまり「成功」は本当に使える
状態であることを保証します。

## ダウンロードされたファイル

Geneformer リポジトリ(`geneformer_hf/`)は `download.sh` が HTTPS で直接
展開します(git 不要)。パッケージコード、V2-104M の重み、トークン/メディアン
辞書を取得します。実際にダウンロードされるファイル:

| ファイル | サイズ | 用途 |
|---|---|---|
| `geneformer_hf/Geneformer-V2-104M/model.safetensors` | 417,571,156 B (~418 MB) | V2-104M 学習済み重み |
| `geneformer_hf/Geneformer-V2-104M/config.json` | 590 B | モデル設定 (18層, hidden 1152, vocab 20275) |
| `geneformer_hf/Geneformer-V2-104M/generation_config.json` | 90 B | 生成設定 |
| `geneformer_hf/Geneformer-V2-104M/training_args.bin` | 5,496 B | 学習引数 |
| `geneformer_hf/geneformer/token_dictionary_gc104M.pkl` | 425,590 B | Ensembl ID ↔ token 辞書(V2) |
| `geneformer_hf/geneformer/gene_median_dictionary_gc104M.pkl` | 1,512,661 B | 遺伝子正規化係数(V2) |
| `geneformer_hf/geneformer/ensembl_mapping_dict_gc104M.pkl` | 3,957,652 B | Ensembl ID 折りたたみ/マッピング(V2) |
| `geneformer_hf/geneformer/gene_name_id_dict_gc104M.pkl` | 1,660,882 B | Ensembl ID ↔ 遺伝子名(V2) |

> V1-10M, V2-104M_CLcancer, V2-316M の重みは既定では**取得しません**
> (V2-104M のみ)。取得するには `./download.sh --model V2-316M` または
> `./download.sh --all`。

## 環境 (uv + venv)

パイプラインはこの環境で検証済みです(インストール手順は上記 **Setup (uv)** を
参照):

```
.venv/                      1.5 GB  uv 作成の Python 3.12 環境
.venv/bin/python            3.12.9
torch                       2.13.0  (MPS built: True, MPS available: True)
transformers                4.46.3  (固定 — 5.x は SpecialTokensMixin が壊れる)
datasets                    4.0.0   (固定 — 5.x は InSilicoPerturber.perturb_data がハング)
```

## 生成された成果物(`analysis_ws/`)

| パス | 内容 |
|---|---|
| `analysis_ws/tokenized/PD_smallint.dataset` | 64,614 細胞のトークナイズ(Hugging Face データセット) |
| `analysis_ws/results/mps_verify/mps_verify_cell_embeddings.csv` | MPS 上の 2,000 細胞埋め込み(768次元) |
| `analysis_ws/results/mps_verify/mps_verify_summary.json` | 検証メトリクス(acc 0.8265, macro F1 0.6883) |
| `analysis_ws/runs/260821_geneformer_cellClassifier_mpsft/ksplit1/checkpoint-8` | MPS 上でのファインチューン済みチェックポイント(8 ステップ) |
| `analysis_ws/results/isp/*_raw.pickle` | in silico perturbation の中間(バッチ別) |
| `analysis_ws/results/isp/isp_state_embs.pkl` | perturbation 用状態埋め込み |

分析スクリプトは `analysis/`:

| スクリプト | 役割 |
|---|---|
| `analysis/04_baseline.py` | tokenize → frozen 埋め込み → ロジスティック回帰 probe |
| `analysis/04b_extract_embeddings.py` | 埋め込みのみ(probe なし)— 下記参照 |
| `analysis/06_finetune.py` | セル分類器の微調整 |
| `analysis/07_in_silico_perturbation.py` | チュートリアル IS perturbation を V2 用に移植 |
| `analysis/smoke_mps.py` | MPS vs CPU の数値一致検証 |
| `analysis/verify_mps_fast.py` | frozen 埋め込み + probe の高速検証 |
| `analysis/verify_mps_finetune.py` | 高速微調整検証 |

### ティッシュごとの実行(04–07)

`04_baseline.py` / `05_figures.py` / `06_finetune.py` /
`07_in_silico_perturbation.py` は、`input/<TISSUE>/h5ad/*.h5ad` のティッシュ
フォルダからワークスペースを自動解決します(`analysis/` に同梱の
`_resolve_tissue.py` ヘルパーを使用)。`ADPD_TISSUE` でティッシュを選ぶと、
`input/<TISSUE>/` を作業ルートとし、`tokenized/`, `results/`, `runs/` を
書き出します:

```bash
# リポジトリ root から実行(相対パス .venv/bin/python を解決するため)
cd /path/to/Geneformer

export GENEFORMER_DIR=$PWD/geneformer_hf        # Geneformer checkout + V2-104M
export GENEFORMER_MODEL=Geneformer-V2-104M

# 1 ティッシュ(input/ 配下の任意の PD_*/AD_*、例: PD_blood / AD_blood)
ADPD_TISSUE=PD_blood  .venv/bin/python analysis/04_baseline.py        # tokenize + 埋め込み + probe
ADPD_TISSUE=PD_blood  .venv/bin/python analysis/05_figures.py          # 図(cached 埋め込みを再利用)
ADPD_TISSUE=PD_blood  .venv/bin/python analysis/06_finetune.py         # セル分類器の微調整
ADPD_TISSUE=PD_blood  .venv/bin/python analysis/07_in_silico_perturbation.py

# または ADPD_ROOT で直接指定も可(従来方式):
#   ADPD_ROOT=$PWD/input/PD_blood
```

> **リポジトリ root から実行すること。** スクリプトは相対パス
> `.venv/bin/python` を使うため、他のディレクトリから実行すると
> `no such file or directory: .venv/bin/python` で失敗します(パス解決自体は
> cwd 非依存で、スクリプト位置から `input/` を特定します)。スクリプトの
> `ADPD_TISSUE=... cmd` 構文は POSIX です。**csh/tcsh** 系シェル(`VAR=value cmd`
> 非対応)の場合は `env` を使ってください:
> ```bash
> env ADPD_TISSUE=PD_blood .venv/bin/python analysis/05_figures.py
> ```

> **補足**
> - `ADPD_TISSUE` は `input/` 下で `h5ad/*.h5ad` を持つフォルダ名を指定。
>   必要な obs 列(`cell_id`, `individual`, `celltype`, `split`)があること。
>   含まれる全ティッシュ(`PD_*`, `AD_*`)がこの条件を満たします。
> - `ADPD_PREFIX` は `ADPD_TISSUE` の後方互換エイリアス。`ADPD_ROOT` を
>   設定するとリゾルバを完全にバイパスします。
> - `06_finetune.py` は `_1/_2/_3` サフィックスでリプリケートを分割し、
>   3 群がすべて空でないことのみを要求します(PD_* は 12/rep、AD_* は
>   通常 9–10/rep)。
> - 各スクリプトはキャッシュされます。再実行は既存の `tokenized/`,
>   埋め込み、チェックポイントを再利用します。再計算は出力を削除して。

### 04 と `04b_extract_embeddings.py` の関係

`04_baseline.py` は **tokenize → 埋め込み抽出 → ロジスティック回帰 probe** を
一括で実行します。埋め込み抽出が最も重く(全細胞に対して 104M モデルを実行)、
`04b_extract_embeddings.py` は **埋め込みのみ** を抽出します(probe はなし —
`05_figures.py` が自身で probe を再実行)。05 が読むのと同じファイルを書き出します:

```
<input>/<TISSUE>/results/embeddings/pretrained_cell_embeddings.csv
```

`05_figures.py` が
`No such file or directory: .../pretrained_cell_embeddings.csv` で失敗したら、
埋め込みが生成されていません(04 未実行、または埋め込み前に停止)。埋め込み
ステージを実行します(冪等 — 既存の `tokenized/` があれば再利用):

```bash
env ADPD_TISSUE=<TISSUE> .venv/bin/python analysis/04b_extract_embeddings.py
```

**tokenize 完了・embedding のみ残っているかの確認方法** — 各成果物(いずれも
`04b`/`04` が `<input>/<TISSUE>/` に作成)を確認:

| 成果物 | 意味 |
|---|---|
| `tokenized/<PREFIX>.dataset/`(ディレクトリで空でない) | tokenize 完了 |
| `results/embeddings/pretrained_cell_embeddings.csv` | 埋め込み完了(05 実行可) |

- **tokenized ディレクトリのみ存在**し埋め込み CSV が無い = tokenize は終わって
  いるが、埋め込み中に停止 → `04b_extract_embeddings.py` を実行(既存の
  `tokenized/` を再利用し、埋め込みのみ実施)。
- **どちらも存在しない** = `04b_extract_embeddings.py` が tokenize + 埋め込みの
  両方を行います。

## Seurat `.rds` から Geneformer 用 `.h5ad` への変換(R → Python)

**R で生成された Seurat オブジェクト(`.rds`)** を、この Geneformer パイプラインで
トークナイズ可能な `.h5ad` にするには、**`rds2h5ad/`** コンバータを使ってください
(R/Python ブリッジで、オブジェクト検査 → サブサンプル+生 counts エクスポート →
マウス→ヒト オルソログ変換 → `.h5ad` 組み立て)。

- 詳細ドキュメント:
  - English: **[`rds2h5ad/rds2h5ad.md`](rds2h5ad/rds2h5ad.md)**
  - 日本語版: **[`rds2h5ad/rds2h5ad-jp.md`](rds2h5ad/rds2h5ad-jp.md)**
  - `.rds` に必要なもの: 生整数の `counts` 層を持つ `RNA` アッセイ(SCT デフォルトは不可)、
    `_1/_2/_3` のレプリケートサフィックス付きサンプル列、細胞型列、マウス遺伝子シンボル。
  - マウス `PD_LN` オブジェクト(195,927 細胞 → サブサンプル後 66,467 → `PD_LN.h5ad`、
    遺伝子保持率 87.9%)でエンドツーエンド検証済み。

`rds2h5ad/` は step 1–3(検査/エクスポート/オルソログ変換)をカバーします。step 4–6
(tokenize+probe、図、fine-tune)は、その `.h5ad` をこのワークスペースの
`analysis/` スクリプトで実行します。

## DirectML (Windows / WSL2 + AMD GPU)

Geneformer は **DirectML** (`torch-directml`) 経由で Windows または WSL2 上の
AMD GPU でも動作します。`GENEFORMER_DEVICE=dml` で有効化します。

### 検証済み構成

| 項目 | 値 |
|---|---|
| GPU | AMD Radeon(TM) 8060S Graphics |
| ドライバー | 32.0.22032.6002 (2025-12-15) |
| 専用 GPU メモリ | 4.0 GB |
| 共有 GPU メモリ | 43.6 GB |
| GPU メモリ (合計) | 47.6 GB |
| `forward_batch_size` | **4** (下記参照) |

### バッチサイズ

埋め込み抽出スクリプト (`04_baseline.py`, `04b_extract_embeddings.py`) は
`forward_batch_size=4` を使用します。これは AMD Radeon 8060S (共有 GPU メモリ
43.6 GB) 向けに調整されています。

データセットは系列長の降順でソートされる（最長の細胞が最初に来る）ため、
最初のバッチには常に最長の系列（最大 3,861 トークン）が含まれます。3,861 トークン
の 4 細胞バッチは `(4, 3861, 768)` の隠れ状態テンソル（約 45 MB）に加えて、
中間のアテンション行列（O(n²) メモリ）を生成します。8 細胞ではアテンション行列が
2 倍になり、この GPU では **segfault**（DML の OOM がクリーンなエラーにならず）
が発生していました。`forward_batch_size=4` が安全な上限です。

GPU の共有メモリが**より少ない**場合はさらに減らしてください（例:
`forward_batch_size=2`）。**より多い**場合は増やせますが、必ず最長系列で
テストしてください。

### DML セットアップ

```bash
# torch-directml のインストール (詳細は docs/directml_setup.md)
uv pip install --python .venv/bin/python torch-directml
# DML バックエンドで実行
GENEFORMER_DEVICE=dml .venv/bin/python analysis/04_baseline.py
```

## ROCm on Windows / WSL2

**この Mac では検証不可**(macOS は WSL2 を実行できない)。実装は完了しています
— 抽象層は CUDA パス経由で AMD ROCm を駆動できる — ただし WSL2 環境自体は
AMD GPU のある Windows マシン上でセットアップする必要があります。

ROCm / WSL2 で動かすには:

```bash
# Windows 側で WSL2 + Ubuntu ディストリを導入。
# WSL2 Ubuntu 内で(AMD ROCm ドライバは Windows 側に導入済み):
wget https://repo.radeon.com/amdgpu-install/latest/ubuntu/jammy/amdgpu-install_*.deb
sudo apt install -y ./amdgpu-install_*.deb
sudo amdgpu-install --usecase=rocm
# PyTorch を ROCm ホイールでインストール:
uv venv .venv --python 3.12
uv pip install --python .venv/bin/python torch --index-url \
    https://download.pytorch.org/whl/rocm6.2
uv pip install --python .venv/bin/python -e geneformer_hf
uv pip install --python .venv/bin/python "transformers==4.46.3"
# ROCm ビルドは torch.cuda.is_available() == True を返す → device 自動 = cuda
GENEFORMER_DEVICE=cuda python analysis/04_baseline.py
```

期待結果: ROCm ビルドで `get_device() == "cuda"`(CUDA API)となり、既存の全
コードパスが AMD GPU 上でそのまま動きます。