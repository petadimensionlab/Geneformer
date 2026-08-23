# Geneformer V2-104M — マルチバックエンド (MPS / CUDA / ROCm) 対応

> **🇬🇧 English version → [`README.md`](README.md)**
> **日本語 (Japanese) | [English](README.md)**

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

## Setup (uv)

```bash
# 1. Geneformer パッケージ + 重みを geneformer_hf/ に取得(curl のみ・
#    git / git-lfs 不要)、マルチバックエンド(MPS/CUDA/ROCm)の device パッチ適用。
./download.sh
cd geneformer_hf
git apply ../patches/geneformer_multibackend.patch   # cwd 非依存。どのユーザーでも動作
cp ../patches/device.py geneformer/device.py        # 新規ファイル(未追跡)
cd ..

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