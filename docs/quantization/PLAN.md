# Geneformer 重み量子化の計画

対象モデル: `Geneformer-V2-104M`（12層・hidden 768・語彙 20,275）
実行マシン: **NVIDIA GB10**（aarch64 / sm_121 / メモリ 119.6 GiB / 20 cores）
ソフト: torch 2.13.0+cu130、bitsandbytes 0.50.1、transformers 4.46.3

---

## この計画の結論（先に要点）

1. このマシンでは **bf16 が一番速く、しかも精度はほぼ無損失**です。
   - 速度: fp32 の **約5倍**
   - 精度: fp32 と比べて MLM loss の差は +0.002、埋め込みのコサイン類似度 0.99998
2. **int8 量子化は安全**ですが、速さは bf16 の半分以下です（fp32 の約1.7倍）。
   - 精度: 埋め込みコサイン 0.9997、MLM loss の差 −0.009
3. **4bit (nf4) は精度が落ちます**。今回の測定では採用できません。
   - 埋め込みコサイン 0.981、遺伝子ランキング top-100 の一致が 87/100
   - 「平均の誤差は小さいが、遺伝子の順位が入れ替わる」タイプの劣化なので、
     in silico perturbation の結論が変わり得ます。
4. したがって **「速くしたいだけ」なら量子化は不要で、bf16 に切り替えるのが最適解**です。
   量子化が本当に効くのは次の場合です。
   - より大きいモデル（V2-316M）を載せたい
   - 配布物・保存サイズを小さくしたい、CPU やエッジで動かしたい
   - QLoRA で微調整したい
   - FP8（Blackwell の fp8 演算）を使いたい

---

## 1. 現状（コードを読んで確認したこと）

- **量子化のしくみはすでに実装済み**です。
  `geneformer_hf/geneformer/perturber_utils.py` の `load_model()` が
  `quantize=False / True / dict(bnb_config=..., peft_config=...)` を受け取ります。
- `EmbExtractor`（埋め込み抽出）と `InSilicoPerturber`（in silico perturbation）は、
  どちらも内部で `load_model()` を呼びます。ですから
  **`model_type="Pretrained-Quantized"` を渡すだけで int8 推論になります**。コード変更は不要です。
- 逆に **4bit (nf4) は外から指定できません**（`model_type` 経由だと int8 に固定）。
  4bit を使うには、`EmbExtractor` と `InSilicoPerturber` に `quantize` 引数を足す小改造が必要です。
- 比較に使える既存の成果物があります。
  - トークナイズ済みデータ: `input/<TISSUE>/tokenized/<PREFIX>.dataset`（6組織）
  - 微調整済みセル分類器: `input/<TISSUE>/runs/*/ksplit1`（6組織）
  - ISP の結果: `isp_result/report_*`、`result/report_*`
- **注意**: `Geneformer-V1-10M` / `V2-104M_CLcancer` / `V2-316M` の
  `model.safetensors` は git-lfs のポインタ（133〜135バイト）で、重みの実体がありません。
  使う前に `./download.sh --model V2-316M` などで取得してください。

---

## 2. 実測ベンチマーク

再現コマンドは付録にあります。スクリプトは `analysis/10_quant_bench.py` と
`analysis/10b_quant_accuracy.py` です。

### 2.1 速さ（batch=8、人工データ、単位は tokens/秒）

| 系列長 | fp32 | bf16 | int8 | nf4 |
|---|---|---|---|---|
| 512 | 43,590 | **205,466** | 69,489 | 127,436 |
| 2048 | 35,753 | **174,696** | 60,290 | 117,262 |
| 4096 | 29,672 | **154,269** | 55,814 | 107,731 |

読み方: **bf16 は fp32 の約5倍、nf4 は約3.4倍、int8 は約1.7倍**です。
速さの順番は **bf16 > nf4 > int8** になります。
int8 が遅いのは、計算を fp16 で行い、そのたびに重みを戻す処理が入るためです。
この GPU は bf16 の演算が非常に速いので、int8 の利点が出ません。

### 2.2 ピークメモリ（batch=8、系列長 2048）

| fp32 | bf16 | int8 | nf4 |
|---|---|---|---|
| 2.32 GiB | 1.18 GiB | 1.10 GiB | 1.06 GiB |

104M のモデルでは、メモリは足かせになりません。
メモリが問題になるのは、もっと大きいモデルを大量バッチ・長い系列で回すときです。

### 2.3 精度（実データ: AD_blood の 256 セル、マスク位置 58,855 箇所）

| 方式 | MLM loss | マスク予測 正解率 | 埋め込み cos 平均 | cos 最小 | 遺伝子 top-100 一致 |
|---|---|---|---|---|---|
| fp32 | 4.0869 | 0.1133 | 1.000000 | 1.000000 | 100/100 |
| bf16 | 4.0876 | 0.1136 | 0.999976 | 0.999971 | 100/100 |
| **int8** | 4.0776 | 0.1133 | 0.999704 | 0.999649 | 96/100 |
| nf4 | 4.2988 | 0.1052 | 0.981277 | 0.970528 | 87/100 |

- bf16 は実質的に無損失です。
- int8 も安全な範囲です。
- **nf4 は劣化がはっきり出ます**。しかも「平均の誤差」ではなく「遺伝子の順位」に効いてきます。

---

## 3. 進め方（フェーズ）

### フェーズ0: 目的の確定と基準の固定（0.5日）
1. 目的を決めます（速くしたい / 大きいモデル / サイズ削減 / QLoRA / FP8）。
2. **fp32 の出力を「基準」として保存**します。比較の土台になります。
   - 埋め込み: `input/*/results/embeddings/`
   - ISP の結果: `input/*/results/isp/`
   - logistic regression プローブの正解率: `analysis/04_baseline.py`
   - 微調整済みセル分類器の評価値: `input/*/runs/*/ksplit1`

成果物: `docs/quantization/BASELINE.md`

### フェーズ1: bf16 と int8 の適用（0.5〜1日）
1. **bf16 化（最優先）**: `torch_dtype=torch.bfloat16` を指定します。
   量子化したモデルには `.to(dtype)` を使ってはいけません。
   いっぽう bf16 は通常のモデルなので、今までのコードがそのまま動きます。
2. **int8**: `model_type="Pretrained-Quantized"` を渡すだけです。
   `perturber_utils` 側が面倒を見てくれます。
3. **nf4 はこの段階では使いません**（2.3 の劣化が確認済みのため）。

### フェーズ2: 精度の判定（1〜1.5日）— ここが本番です
`analysis/10b_quant_accuracy.py` を拡張し、複数の組織で測ります
（AD_blood / AD_spleen / PD_spleen / AD_liver）。

| # | 指標 | 合格ライン | なぜこの基準か |
|---|---|---|---|
| G1 | MLM loss | fp32 との差が +0.02 以内 | int8 は −0.009、nf4 は +0.212 なので、この線で分けられます |
| G2 | マスク予測の正解率 | 差が −0.002 以内 | 同上 |
| G3 | セル埋め込みの cos | 平均 0.9995 以上・最小 0.999 以上 | int8 は合格、nf4 は不合格 |
| G4 | プローブ（celltype 分類） | 正解率・macro-F1 の差が 0.01 以内 | 下流タスクでの実害を見ます |
| G5 | **ISP の順位一致**（`Shift_to_goal_end`） | スピアマン相関 0.99 以上、top-20 の一致 18/20 以上、符号一致 98% 以上 | **最重要**。生物学的な結論を守るため |
| G6 | 微調整済み分類器の評価値 | macro-F1 の差が 0.01 以内 | 結論に直結します |
| G7 | 少数派の最悪ケース | 希少な celltype や短いセル（100トークン未満）でも cos 0.99 以上 | 平均では隠れる劣化を探すため |

**nf4 は G3 と G5 で落ちる見込み**です。4bit を使いたい場合はフェーズ3に進みます。

成果物: `docs/quantization/ACCURACY.md`

### フェーズ3: キャリブレーション型・静的量子化（条件付き・1〜3日）
フェーズ2で nf4 が落ち、それでも 4bit 級のサイズと速度が必要なときだけ実施します。

1. **キャリブレーション用データの設計**（Geneformer 特有の注意点）
   - 語彙は遺伝子IDと「値の順位トークン」で固定されています。量子化誤差は
     トークンの出現頻度に強く依存するので、**ランダムな人工データでは意味がありません。実セルを使ってください**。
   - 使い方: `input/<TISSUE>/tokenized/*.dataset` から celltype と系列長で層化して 512〜1024 セル。
   - **リーク防止**: 最終評価に使う AD_blood / AD_spleen をキャリブレーションに使わないこと。
     PD_spleen / AD_liver / AD_smallint を使います。
2. 手法の候補（有利な順）
   - **FP8 W8A8**（Blackwell は fp8 演算ユニットを積んでいるので一番期待できます。
     `llmcompressor` + `compressed-tensors`。aarch64 向けホイールの有無を確認）
   - **HQQ**（キャリブレーション不要で実装が軽い）
   - **GPTQ / AWQ（W4A16）**（`gptqmodel` / `autoawq` の aarch64 対応を確認）
   - **torchao**（int8 / int4 / fp8 と `torch.compile` の融合）
3. 出力形式を決めます（HF の `quantization_config` か、`compressed-tensors` 形式か）。
4. **保護する層**: `LayerNorm`、埋め込み、分類ヘッド、pooler は量子化から外します。
   Geneformer は埋め込みの平均をセルの表現に使うので、埋め込みテーブルの誤差が最も効きます。

### フェーズ4: QLoRA による微調整（目的が「QLoRA で微調整」の場合・1〜2日）
- `analysis/06_finetune.py` の `CellClassifier(quantize={...})` を使います（実装済み）。
  `bnb_config` に 4bit-nf4 か 8bit、`peft_config` に
  `LoraConfig(r=64, alpha=128, dropout=0.1, task_type=TokenClassification)` を指定します。
- 注意点
  - `_ensure_complete_checkpoint()` は `adapter_model.bin` がある run を再利用可能とみなします。
    LoRA の run も同じ扱いになるので、`ksplit1` の名前が衝突しないようにしてください。
  - 勾配チェックポイントと 4bit の併用には `enable_input_require_grads()` が必要です
    （`load_model` 内で処理済み）。
  - **このマシンではメモリ削減の効果はほとんどありません**（119 GiB あり、fp32 でも 2.3 GiB）。
    利点は速さと、MPS マシンと同じ制約で回せる再現性です。
- 比較対象は既存の fp32 微調整済み `ksplit1`（6組織）で、macro-F1 の差 0.01 以内を合格とします。

### フェーズ5: 統合と展開（1日）
1. 環境変数で切り替えられるようにします（`device.py` の `GENEFORMER_DEVICE` と同じやり方）。
   `GF_QUANT=none|bf16|int8|nf4|fp8` を `analysis/_isp_common.py` のモデル生成部分で読みます。
2. 4bit を使う場合だけ、wrapper に `quantize` を公開するパッチを当てます。
   - `EmbExtractor.__init__` と `InSilicoPerturber.__init__` に `quantize=None` を足し、
     `load_model(..., quantize=self.quantize)` に渡します。
   - 既定は `None`（＝今までと同じ動き）にして、既存の結果を壊しません。
   - パッチは `patches/geneformer_quantization.patch` として管理します。
3. **カナリア実行**: 1組織（AD_spleen）の ISP を量子化モデルで回し、
   `isp_result/report_AD_spleen` の既存結果とフェーズ2の G5 を再確認します。
4. ドキュメントを更新します（`README-jp.md`、`docs/quantization/RESULTS.md`）。

---

## 4. つまずきやすい点

1. **量子化したモデルに `.to()` / `.half()` / `.cuda()` を呼ばないでください。**
   `load_model` は量子化時にデバイス移動をスキップする作りになっています
   （`perturber_utils.py` の196行目付近）。
2. **`device_map={"":0}` が必須**です（bnb のロード時）。
3. **`datasets==4.0.0` を維持してください**（5.x では `perturb_data` が止まります）。量子化の検証でも同じ制約です。
4. ISP の時間は「perturbation データセットの構築」と「forward」に分かれます。
   量子化で速くなるのは **forward だけ**です。`max_ncells` を絞るなどの既存の運用知識はそのまま有効です。
5. **nf4 の劣化は平均ではなく順位に出ます**（top-100 の一致が 87/100）。
   遺伝子ランキングが成果物である以上、MLM loss だけで合否を決めてはいけません。
6. **LFS ポインタに注意**。実体のないモデルを読むと `file signature not found` などで失敗します。
7. **再現性**: 量子化カーネルは fp32 とビット単位では一致しません。
   比較するときは必ず「同じ seed・同じバッチ構成」で fp32 と量子化を並べて測ってください。

---

## 5. スケジュールの目安

| フェーズ | 内容 | 見積り |
|---|---|---|
| フェーズ0 | 目的の確定・基準の固定 | 0.5日 |
| フェーズ1 | bf16 / int8 の適用 | 0.5〜1日 |
| フェーズ2 | 精度判定（4組織・7指標） | 1〜1.5日 |
| フェーズ3 | 静的量子化（条件付き） | 1〜3日 |
| フェーズ4 | QLoRA（その目的の場合） | 1〜2日 |
| フェーズ5 | 統合・カナリア・文書化 | 1日 |

**最小構成（おすすめ）**: フェーズ0 → 1（bf16 と int8）→ 2 → 5。
3〜4日で「速さは bf16 で5倍、配布や保存は int8 で精度0.99」という結論まで到達できます。
nf4 と FP8 は、フェーズ2の G5（ISP の順位一致）を通ったときだけ採用してください。

---

## 付録: 再現コマンド

```bash
cd ~/workspace/research/Geneformer
.venv/bin/python analysis/10_quant_bench.py         # 2.1 と 2.2（速さ・メモリ）
.venv/bin/python analysis/10b_quant_accuracy.py     # 2.3（精度・埋め込み一致）
.venv/bin/python analysis/10c_quant_memscale.py     # バッチサイズごとのピークメモリ
.venv/bin/python analysis/10d_real_workload_mem.py  # 実データでのピークメモリ
.venv/bin/python analysis/10e_train_mem_probe.py    # 微調整時のピークメモリ

# 既存のしくみで int8 の ISP を回す例（フェーズ1の実運用形）
GENEFORMER_DIR=geneformer_hf ADPD_TISSUE=AD_spleen GF_QUANT=int8 \
  .venv/bin/python analysis/07_ad_spleen_early_isp.py
```
