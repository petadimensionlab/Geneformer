# 再現手順（How to run）

各臓器の in silico perturbation を再実行する手順です。

## 前提
- Geneformer リポジトリ（`geneformer_hf/`）が `./download.sh` で取得済み
- `.venv` に依存解決済み、`datasets==4.0.0`
- `input/<TISSUE>/` にトークナイズ済みデータセット（`tokenized/<PREFIX>.dataset`）と
  ファインチューン済みセル分類器（`runs/**/ksplit*`）があること

## 共通実行コマンド

```bash
# 例: <TISSUE を設定>（例: 脾臓の AD）
GENEFORMER_DIR=geneformer_hf ADPD_TISSUE=AD_spleen \
  .venv/bin/python analysis/07_ad_spleen_early_isp.py
```

- `ADPD_TISSUE` は `input/<TISSUE>` を指定（`AD_spleen`, `AD_blood`, `AD_brain`,
  `AD_smallint`, `PD_spleen`, `AD_liver`, `AD_BM` 等）。
- 出力は統一規則で `input/<TISSUE>/results/isp/<experiment>/` に書かれます。

## 各実験の実行

| 実験 | コマンド（ADPD_TISSUE=...） |
|---|---|
| ad_smallint | `AD_smallint` |
| ad_brain | `AD_brain` |
| ad_blood | `AD_blood` |
| ad_spleen | `AD_spleen` |
| pd_spleen | `PD_spleen` |
| ad_liver | `AD_liver`（`07e_ad_liver_perturbation.py`） |
| ad_bm | `AD_BM`（`07f_in_silico_perturbation_AD_BM.py`） |

> liver / bm は独自スクリプト（`07e_` / `07f_`）を使用。実行コマンドは各スクリプトの docstring を参照。

## smoke 実行（短時間で動作確認）

```bash
# 3m タイムポイントのみ・先頭 3 遺伝子で実施
IS_TIMEPOINTS=3m IS_MAX_GENES=3 \
  GENEFORMER_DIR=geneformer_hf ADPD_TISSUE=AD_spleen \
  .venv/bin/python analysis/07_ad_spleen_early_isp.py
```

## 出力の確認

実行後、`input/<TISSUE>/results/isp/<experiment>/<experiment>_early_isp_stats_combined.csv`
に全遺伝子の `Shift_to_goal_end` が出力されます。

## レポート生成

結果の可視化・HTML レポート・officecli 用 PPT プロンプトは `result/report/` にあります
（各 build、render、visualize スクリプト）。一覧は[リポジトリの README](../../README.md) にあります。

文章としてのレポート（DOCX、日本語と英語）が必要な場合は、`templates/` の雛形を使います。

```bash
.venv/bin/python templates/template_report_figures.py --lang both
.venv/bin/python templates/template_build_report_docx.py --lang both
```

`templates/README_templates_ja.md` に、入力（前段の出力）と実行順、日本語 PDF への変換方法を書いてあります。
節の並べ方は[レポートの標準構成](report_template.md)、書き方の決まりは[遵守項目](checklist.md) の C-1（レポートを書くときに守ること）を参照してください。

## 計算上の注意（必読）

1. **`datasets==4.0.0` に固定**（`>=5` は `perturb_data` がハング）。
2. **`nproc=1`** が安定（`IS_NPROC` で上書き可能だが RAM 競合に注意）。
3. `max_ncells`（`IS_MAX_CELLS`）を上げすぎない（OOM 防止）。RAM 見積りがプリントされます。