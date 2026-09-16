# モデル比較・文献照合テンプレート集

Geneformer の 2 モデル（または 2 条件）の出力を比較し、対象の遺伝子と細胞型について
文献と突き合わせて、読めるレポート（DOCX/PDF、日英）にまとめるための再利用用テンプレート。

**これらは雛形（TEMPLATE）であり、完成した解析ではありません。**
ファイル内に埋め込まれている遺伝子判定・レポート本文・文献検索式は、
`AD_spleen` の 104M vs 316M 実行（2026-09）の**記入例**です。
別の組織・モデル・時点に使うときは、数値を再計算し、PMID を再検証してください。

| テンプレート | 役割 | 主な出力 |
|---|---|---|
| `template_expression_by_celltype.py` | 細胞種・病態別の遺伝子検出率（発現の広さ） | `docs/quantization/expression/<TISSUE>.json` |
| `template_isp_model_compare.py` | 数値精度の差、細胞サンプリングの差、モデルの差を分離し、符号の反転を数える | `isp-model-comparison/{gene_table.csv,summary.json}` |
| `template_literature_check.py` | Europe PMC 検索（遺伝子ごとの検索式を内蔵） | `isp-model-comparison/literature/*.json` |
| `template_literature_adjudicate.py` | 文献判定の統合と PMID の機械照合 | `isp-model-comparison/adjudication.{csv,json}` |
| `template_report_figures.py` | 図の生成（日英） | `isp-model-comparison/figures/<lang>/*.png` |
| `template_build_report_docx.py` | レポートの生成（DOCX、日英） | `isp-model-comparison/REPORT-*.docx` |

## 実行順

```bash
.venv/bin/python templates/template_expression_by_celltype.py
.venv/bin/python templates/template_isp_model_compare.py
.venv/bin/python templates/template_literature_check.py --batch all
.venv/bin/python templates/template_literature_adjudicate.py
.venv/bin/python templates/template_report_figures.py --lang both
.venv/bin/python templates/template_build_report_docx.py --lang both
```

PDF が必要な場合（日本語フォントは Noto Sans CJK JP を使用）:

```bash
cd docs/quantization/isp-model-comparison
soffice --headless --norestore --convert-to pdf \
    REPORT-ISP-104M-vs-316M.docx REPORT-ISP-104M-vs-316M-jp.docx
```

## 前提と注意

- **出力はコミットしません。** `docs/quantization/isp-model-comparison/` と
  `docs/quantization/expression/` は `.gitignore` 済み（テンプレートのみ公開）。
- **`template_isp_model_compare.py` の前段**として ISP 本体の実行が必要です
  （`analysis/07_*_early_isp.py` など）。テンプレートは ISP の出力 CSV を読みます。
- **`template_literature_check.py` は Europe PMC REST を使います。** ページ本文の
  スクレイピングは行いません（`web_extract` のバックエンド未設定でも動くようにするため）。
- **PMID は必ず引き直して検証**します。存在確認だけでなく、主張した遺伝子と
  取得した題名が対応するかまで見てください（対応しない例が実際にありました）。
- **`template_build_report_docx.py` の文献番号は 2 パス方式**で、初出順の連番になります
  （キー基準の固定番号だと 1,2,3,5,8… と飛んで読者に壊れて見えます）。
- 環境変数: `ADPD_TISSUE`（既定 AD_spleen）、`GENEFORMER_DIR`、`GENEFORMER_MODEL`、`GF_DTYPE`。
- 日本語の図・文書には CJK フォント（Noto Sans CJK JP）が必要です。フォントが
  解決できない場合、スクリプトは黙って□を出さずにエラーで停止します。
