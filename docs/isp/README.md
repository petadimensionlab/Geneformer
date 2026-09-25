# ISP のルールブック（索引）

> **このページ群は「ISP（in silico perturbation）解析の規則」を置く場所です。**
> 規則は**個別の実行に依存しない形**で書きます。実行ごとの結果・解釈・実測値は
> **[実行の記録](#2-実行の記録規則ではない)**に置きます。両者を混ぜないでください。

> ⛔ **作業を始める前に、必ず [記述ルール](style-guide.md) を読んでください。**
> 用語、表記、判定の語彙の唯一の基準です。読む順序は次のとおりです。
>
> | 順 | 文書 | 何が書いてあるか |
> |---|---|---|
> | **1** | **[記述ルール](style-guide.md)** | **A データ分析、B 技術・コード、C 共通（表現と表記）。書き方の唯一の基準** |
> | 2 | [評価基準](criteria.md) | 数値の基準（E1〜E6）と根拠の種類（A・B・C） |
> | 3 | [用語と手法の定義](glossary.md) | 用語の意味 |
> | 4 | [方法と共通設定](methods.md) | 操作と統計量の定義、null 分布の数式 |
> | 5 | [遵守項目](checklist.md) | 実行時の行動規範 |
> | 6 | [レポートの標準構成](report_template.md) | レポートの節立て |
> | 7 | [図表の視覚設計ルール](design.md) | スライドとポスターのレイアウト・文字サイズ・密度 |
> | 8 | [再現手順](how_to_run.md) | 実行の手順 |
> | 9 | [執筆者名簿とページ担当](owners.md) | 誰が編集するか（**2026-09-25 に AG-GH 1 名へ一本化**） |
> | 10 | [競合の記録](conflicts.md) | 競合と裁定の記録（規則を変えた理由の記録でもあります） |
>
> 記述が食い違ったときは、勝手に上書きせず [競合の記録](conflicts.md) に追記してください。

---

## 1. 規則（このページ群）

| 文書 | 何が書いてあるか |
|---|---|
| [記述ルール](style-guide.md) | 用語、表現、数値の表記、根拠の種類、判定の語彙、リポジトリ操作 |
| [評価基準](criteria.md) | 何をもって成功��するか（E1〜E6）、基準の作り方、実行前後の手順 |
| [用語と手法の定義](glossary.md) | トークン化・埋め込み・`Shift_to_goal_end` の式・ノイズ床・状態分離 |
| [方法と共通設定](methods.md) | 操作と統計量の定義、null 分布の数式、共通パラメータ |
| [遵守項目](checklist.md) | ���行時の行動規範（解析前・解釈・文献・文書・運用） |
| [レポートの標準構成](report_template.md) | レポートの節立て（背景 → 目的 → 問い → プロセス → 結果 → 解釈） |
| [視覚設計ルール](design.md) | スライドとポスターの視覚設計の規則（出典つき） |
| [再現手順](how_to_run.md) | 実行の手順 |
| [執筆者名簿](owners.md) | 誰が編集するか、規則の置き場所 |
| [競合の記録](conflicts.md) | 競合の経緯と裁定、旧節番号の対応 |

## 2. 実行の記録（規則ではない）

| 記録 | 何が書いてあるか |
|---|---|
| [結果のまとめ](records/results.md) | 実行ごとの結果と解釈、**基準の根拠になった実測**、限界 |
| [PD の現状](records/pd.md) | PD の実施状況と残タスク |
| [PD 多領域アトラス](records/pd_atlas.md) | ヒト多領域アトラスの実施記録（null 比較を含む） |
| [肝臓](records/ad_liver.md)、[骨髄](records/ad_bm.md) | 個別の実行の記録 |
| [文献の要約](records/ad_morabito2021.md) | 文献のセクション別要約 |
| [視覚設計の適用記録](records/design_record.md) | ポスターに規則を当てたときの実測（版ごとの検証） |

**規則と記録の分け方**

- 規則ページには**実行に依存しないこと**を書きます。**実測値は書きません**（実行ごとに変わるため）。
- 記録ページには**実行ごとの測定値・判定・出所**を書きます。規則を書き換えません。
- 迷ったら、[執筆者名簿](owners.md) §4 の「内容ごとの正典」に従います。

---

## 3. 出力先の命名規則（規則）

各マシンで出力先が異なっていたため、次の規則に統一しました。**新しい実行はこの形にしてください。**

```
input/<TISSUE>/results/isp/<experiment>/
    <experiment>_early_isp_stats_combined.csv   # 全遺伝子の Shift_to_goal_end
    isp_state_embs*.pkl                          # 共有状態 embedding
    <timepoint>/<Gene>/isp_stats.csv             # 各遺伝子の raw 統計
```

`<experiment>` は **`<disease>_<tissue>` の小文字**です（例: `<疾患>_<臓器>` の形）。

**対応スクリプト** `analysis/_isp_common.py` の `resolve_experiment()` / `isp_dir_for()` /
`combined_csv_path()` ���この規則を返します。

---

## 4. 共通リファクタ（`analysis/_isp_common.py`）

各実行のスクリプトに重複していた処理を集約しました。

| 関数 | 役割 |
|---|---|
| `warn_if_multiproc()` | nproc>1 の警告（`datasets<5` 前提） |
| `check_datasets_version()` | `datasets>=5` で `perturb_data` がハングする警告 |
| `estimate_perturb_ram()` | `max_ncells` の RAM 見積り |
| `load_gene_dicts()` | token→ENSG / ENSG→symbol / symbol→ENSG 辞書 |
| `resolve_classifier_dir()` | 微調整済みモデルの解決（stale パス / 最新 ksplit fallback） |
| `scan_pool_presence()` / `_per_key()` | 対象プールでの遺伝子発現チェック |
| `perturb_one_gene()` | 1 遺伝子削除と `goal_state_shift` 統計（`Shift_to_goal_end` を返す） |
| `write_combined()` | 全遺伝子結果の結合 CSV 出力（タイムポイント順位つき） |

---

## 5. 実行時の共通設定（規則）

- **1 遺伝子ずつ**削除します（`combos=0`。リストを渡すと「全遺伝子を同時に持つ細胞」条件になり細胞が空になります）
- **`perturb_type = delete`**、`emb_mode="cls"`、`model_version="V2"`
- 状態 embedding はタイムポイントごとに 1 回計算し、全遺伝子で共有します
- **`nproc=1`**（安定）、**`datasets==4.0.0`**（`>=5` は `perturb_data` がハングします）

詳細と環境変数は [方法と共通設定](methods.md) と [再現手順](how_to_run.md) にあります。

---

## 6. 成果物の置き場所（規則）

統一後の結果 CSV（実計算値）は **Git にコミットしません**。ローカルの
`input/<TISSUE>/results/isp/<experiment>/<experiment>_early_isp_stats_combined.csv` にあります（`.gitignore` 対象）。

レポート用に整形したデータは `result/report/data/`、図は `result/report/figures/` にあります
（いずれも計算成果物のため Git 対象外）。**実行スクリプト・このページ群・レポートテンプレートは Git 追跡対象です。**
