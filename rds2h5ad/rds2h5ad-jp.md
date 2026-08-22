# rds2h5ad — Seurat `.rds` → Geneformer `.h5ad` パイプライン

マウス single-cell RNA-seq の Seurat オブジェクトを、Geneformer 対応の
`.h5ad` へ変換します。このフォルダは ADPD パイプライン
(`analysis/mouse2human.md`) の step 1–3 を一式備えます:
inspect → サブサンプル/エクスポート → 直系種(オルソログ)マッピング +
`.h5ad` 組み立て。step 4–6(tokenize/probe、図、fine-tune)は
torch/GPU 環境が必要であり、**ここでは構築していません**。

---

## 1. 構成

```
rds2h5ad/                        (スクリプト+ドキュメントは git 管理、データ/環境/出力は無視)
├── setup.sh                     環境の構築/検証(冪等・再実行安全)
├── run.sh                       ./run.sh <tissue> — 01 → 02 → 03 を実行
├── 01_inspect.R                 Seurat オブジェクトの検査
├── 02_export.R                  サブサンプル + 生 counts のエクスポート (mtx + csv)
├── 03_map_and_build_h5ad.py     マウス→ヒトオルソログ変換、フィルタ、.h5ad 生成
├── env.sh                       $RSCRIPT/$PYTHON/... を出すために source(生成物)
├── manifest.json                環境の出所情報(生成物)
├── .Rprofile renv.lock renv/    renv プロジェクト(Seurat 5、R 4.6)
├── pyenv/                       Python 3.12 venv(anndata/numpy/pandas/scipy)
├── rds_data/<tissue>/           入力 .rds(ティッシュごとに 1 フォルダ)
│   └── PD_LN/*.rds
└── ws/<tissue>/                 ティッシュごとの出力
    ├── logs/01_inspect.log 02_export.log 03_map.log
    ├── export/data/{counts.mtx,genes.csv,cells.csv,obs.csv,export_meta.json}
    └── h5ad/<tissue>.h5ad (+ <tissue>_mapping_summary.json)
```

このフォルダ外の関連パス:

| パス | 役割 |
|---|---|
| `Geneformer/geneformer_hf` | Geneformer checkout(token/gene 辞書、V2-104M チェックポイント) |
| `Geneformer/analysis/ref` | MGI マウス-ヒト相同性テーブル + sha256 |
| `Geneformer/ref` | `analysis/ref` へのシンボリックリンク(03 は `parent.parent/ref` を見る) |

## 2. 前提(検証済み環境)

- macOS arm64(M4、RAM 128 GB)
- R ≥ 4.4 — 本機は Homebrew **R 4.6.0**。`mouse2human.md` は R 4.4 を指定して
  いますが、R 4.6 を project owner が承認済み。
- Apple clang(Xcode CLT)+ `gfortran`(Homebrew `gcc`)— Seurat スタックは
  ソースからコンパイルします(Posit PPM は macOS に**ソースのみ**提供、
  バイナリなし)
- venv 用 `uv`(または python3.12)
- rds に必要なもの: **生 counts を持つ RNA アッセイ**(SCT デフォルトは不可)、
  `_1/_2/_3` のレプリケートサフィックスを持つサンプル列、細胞型列、
  マウス遺伝子シンボル

## 3. セットアップ

```bash
cd Geneformer/rds2h5ad
bash setup.sh                 # 冪等。R パッケージ初回導入は ~5-40 分
```

構築内容(特に明記しない限りすべてこのフォルダ内):

| 要素 | 結果 |
|---|---|
| R base | システム Homebrew R 4.6.0 |
| renv プロジェクト | Seurat 5.5.1, SeuratObject 5.4.0, Matrix, jsonlite(+依存 131)。lockfile は `renv.lock` |
| Python venv | `pyenv/` — CPython 3.12.9, anndata 0.13.2, numpy 2.5.2, pandas 3.0.5, scipy 1.18.0 |
| ref シンボリックリンク | `Geneformer/ref → analysis/ref`(なければ作成) |
| 生成ファイル | `env.sh`, `manifest.json`, cwd 非依存の `.Rprofile` |

`setup.sh` は最後に自己検証(R バージョン、renv アクティベーション、
Seurat ロード、python import)を行い、失敗時は非ゼロで終了します。

## 4. パイプラインの実行

```bash
cd Geneformer/rds2h5ad
./run.sh PD_LN
```

- 引数 1 = `rds_data/` 内のフォルダ名。そこで最初に見つかった `.rds` を使用。
- オプション引数: `./run.sh <tissue> <sample_col> <ct_col> [max_per_group] [assay]`
  — デフォルトは `samples2 final_cell.type 200 RNA`。アッセイは必ず `RNA`
  を指定すること(SCT の "counts" レイヤーは sctransform 補正済みで、
  トークナイズを**エラーなく壊す**)。
- `run.sh` は 01 のログを再読み込みし、指定した列のいずれかがオブジェクトの
  meta.data に無ければ **02 の前に中止**します(11GB 再読込の回避)。
- 出力: `ws/<tissue>/{logs,export,h5ad}`。全ステップはキャッシュされ、
  再実行は既存出力を再利用します(再計算は出力を削除して)。

新規ティッシュの追加: `mkdir rds_data/<名前>` に `.rds` を入れて
`./run.sh <名前> [列名...]`。列名はオブジェクトごとに異なるので、まず
`ws/<名前>/logs/01_inspect.log` の meta.data テーブルと候補列を確認すること。

`run.sh` を使わない手動実行:

```bash
source env.sh
mkdir -p "$ADPD_ROOT/PD_LN"/{logs,export,h5ad}
$RSCRIPT   01_inspect.R "$RDS" | tee "$ADPD_ROOT/PD_LN/logs/01_inspect.log"
$RSCRIPT   02_export.R "$RDS" "$ADPD_ROOT/PD_LN/export/data" samples2 final_cell.type 200 RNA \
  | tee "$ADPD_ROOT/PD_LN/logs/02_export.log"
$PYTHON    03_map_and_build_h5ad.py "$ADPD_ROOT/PD_LN/export/data" \
  "$ADPD_ROOT/PD_LN/h5ad/PD_LN.h5ad" "$GENEFORMER_DIR/geneformer" \
  | tee "$ADPD_ROOT/PD_LN/logs/03_map.log"
```

### 環境変数(`env.sh`)

| 変数 | 意味 |
|---|---|
| `ADPD_TOOLS` | このフォルダ |
| `R_BIN` | R インストールディレクトリ(`/opt/homebrew/bin`) |
| `RENV_PROJECT`, `R_PROFILE_USER` | renv アクティベーション(任意の cwd で有効) |
| `RSCRIPT` | renv ライブラリ有効の Rscript |
| `PYTHON` | venv の python |
| `GENEFORMER_DIR` | Geneformer checkout(`../geneformer_hf`) |
| `GENEFORMER_MODEL` | チェックポイント名(`Geneformer-V2-104M`、step 4+) |
| `RDS` | 自動検出された rds(便宜用。`run.sh` はティッシュ単位で解決) |
| `ADPD_ROOT` | 出力ベース(デフォルト `ws/`) |
| `ADPD_PREFIX` | デフォルト `PD_LN` |

## 5. 検証済み実行 — PD_LN(2026-08-22)

入力: `rds_data/PD_LN/LNPD_soup.corrected_all.integrated_doublet.strictly.removed_scDblFinder_annotated_calculated.rds`(11 GB)

| Step | 結果 |
|---|---|
| 01 | 195,927 cells, 33,750 features, アッセイ RNA+SCT(デフォルトは SCT なので RNA を渡す), サンプル 36, 細胞型 23 種, マウスシンボル, counts integer-like: TRUE |
| 02 | 66,467 cells(828 の donor×celltype 群、群あたり ≤200)× 17,139 genes |
| 03 | 15,071/17,139 遺伝子保持(**87.9%**,期待値 85–91% 内), donor leakage check: PASS, `ws/PD_LN/h5ad/PD_LN.h5ad`(828 MB), 23 種すべてで train/eval/test のバランスの良いクロス集計 |

## 6. 構築中に発見した落とし穴と対策(`setup.sh` に組み込み済み)

すべて自動処理済み。「直して」戻さないための記録です:

1. **PPM は macOS にソースのみ提供** — mac 向けプリビルドパッケージが無く、
   Seurat スタックはローカルコンパイル(M4 では高速。初回フルビルド ~40 分)。
2. **Homebrew R/gcc のバージョンずれ** — R bottle は Makeconf に gcc-15
   ランタイムパスを焼き込み済みだが Homebrew `gcc` は 16 で、全 C
   パッケージがリンク時に `library 'emutls_w' not found` で失敗。
   対策: `setup.sh` が期待されるバージョン付きライブラリディレクトリを実在する
   ディレクトリへシンボリックリンク(Homebrew Cellar 内)。
3. **本機の `mamba` バイナリが壊れている**(conda/mamba の不一致)—
   `setup.sh` は使用前に mamba/micromamba/conda を健全性チェック
   (現在はシステム R を使うため conda env なしで無関係)。
4. **renv 1.2.4 の API 変更** — `init(prompt=)` が削除され、`renv/bin/`
   シムも生成されない。アクティベーションは `R_PROFILE_USER` +
   `RENV_PROJECT`(`env.sh` で設定)で行い、`setup.sh` は renv 標準の
   相対パス `.Rprofile` を cwd 非依存版で上書き。
5. **`Rscript -e` の二重アンエスケープ**(R 4.6)— `-e` は文字列アン
   エスケープを 1 回余分に処理し、正規表現のバックリファレンスを壊す。
   `setup.sh` の検証は一時ファイル実行方式(ファイルモードは標準的に解析)。
6. **`03` の ref パス** — `Path(__file__).parent.parent/ref` を見るため、
   スクリプトが `rds2h5ad/` にあると `Geneformer/ref` を解決する。
   上記のシンボリックリンクで記録済み MGI sha256 チェックを維持。

## 7. 補足

- `rds2h5ad/` のスクリプトとドキュメントは git 管理。一方 `rds_data/`
  (11 GB rds)、`renv/`、`pyenv/`、`ws/`、生成ファイル(`env.sh`,
  `manifest.json`, `.Rprofile`, `renv.lock`)は git 無視(workspace 直下の
  `.gitignore`)。`Geneformer/ref` のシンボリックリンクはマシン依存のため
  未追跡のまま(`setup.sh` が再作成)。
- 沿革: 環境は当初 `~/adpd-tools` に構築、スクリプトは `analysis/`、
  データは `input/PD_LN/` にありました。2026-08-22 にすべてへこのフォルダへ
  移設。
- リポジトリ直下にある既存の `Geneformer/.venv`(python 3.12 + geneformer +
  torch、以前のランの残り)はそのままです。step 4–6 の出発点に利用可能で、
  そのスクリプト(`04_baseline.py`, `05_figures.py`, `06_finetune.py`)は
  `analysis/` に残っています。
