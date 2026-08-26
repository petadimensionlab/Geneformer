# In Silico Perturbation Wiki — 早期疾患（AD/PD）発症因子の仮想遺伝子削除スクリーニング

> この Wiki は、Geneformer による **in silico perturbation（遺伝子仮想削除）** を複数臓器・複数 PC で実施した際の **スクリプト・出力先・実施結果** を整理したものです。
> 出力先を統一した命名規則・共通リファクタ（`analysis/_isp_common.py`）・各臓器の設定と結果をまとめています。

---

## 1. 概要

**目的**: アルツハイマー病（AD）/ パーキンソン病（PD）の **早期発見・早期発症因子** を、単細胞トランスクリプトームから探索する。

**手法**: Geneformer V2-104M の埋め込み表現を用いて、**疾患状態の細胞から 1 遺伝子を仮想的に削除（delete）** し、その embedding が **健常（WT）状態へどれだけ近づくか**（`Shift_to_goal_end`）を評価。

- `Shift_to_goal_end > 0`（正）: 削除で疾患細胞が WT に近づく = **早期駆動因子 / 介入標的候補**
- `Shift_to_goal_end < 0`（負）: 削除で WT から離れる = **健常状態の維持因子**

### 対象臓器と実施 PC（3 台で並行実施）

| 実験 | 疾患 | 臓器 | スクリプト | モデル | 実施ホスト |
|---|---|---|---|---|---|
| `ad_smallint` | AD | 小腸（gut-brain axis） | `07_ad_smallint_early_isp.py` | Fine-tuned CellClassifier | PC-A |
| `ad_brain` | AD | 脳（ミクログリア） | `07_ad_brain_early_isp.py` | Fine-tuned CellClassifier | PC-A |
| `ad_blood` | AD | 血液（末梢免疫） | `07_ad_blood_early_isp.py` | Fine-tuned CellClassifier | PC-B |
| `ad_spleen` | AD | 脾臓（免疫） | `07_ad_spleen_early_isp.py` | Fine-tuned CellClassifier | PC-B |
| `ad_ln` | AD | リンパ節（CD4+ T テスト） | `07_ad_ln_early_isp.py` | Fine-tuned CellClassifier | ローカル |
| `pd_spleen` | PD | 脾臓（免疫） | `07_pd_spleen_early_isp.py` | **Pretrained** V2-104M | PC-C |

> **PD の完成度は AD より低い**（データ=脾臓のみ・fine-tune 済みモデル未完成=Pretrained 実行）。
> 詳細は **[PD の現状・残タスク（pd_spleen）](pd.md)** を参照。
| `ad_liver` | AD | 肝臓 | `07e_ad_liver_perturbation.py` | Fine-tuned CellClassifier | PC-C |
| `ad_bm` | AD | 骨髄 | `07f_in_silico_perturbation_AD_BM.py` | Fine-tuned CellClassifier | PC-C |

> **AD_LN（リンパ節）**: `07_ad_ln_early_isp.py` は CD4+ T 細胞の小規模テスト（CD28/STAT3/FOXP3）。
> `input/AD_LN` の h5ad/tokenized/fine-tune モデルは未配置（スクリプトのみ）のため、実行にはデータの配置が必要。
> **LN は AD/PD ともスクリプト（AD は本スクリプト）はあるが、データ・レポートは未整備**（PD_LN は `07b_pd_perturbation.py` が命名上の対象、実データなし）。

> 各 PC で作成された script 名・出力先が異なっていたため、本 Wiki と併せて **出力先を統一**し、共通部分を `_isp_common.py` に集約しました（後述）。

---

## 2. 出力先の統一（命名規則）

3 台の PC で出力先（`input/<TISSUE>/results/isp/` 配下）がバラバラだったため、**以下の共通規則**に統一しました。

```
input/<TISSUE>/results/isp/<experiment>/
    <experiment>_early_isp_stats_combined.csv   # 全遺伝子の Shift_to_goal_end
    isp_state_embs*.pkl                          # 共有状態 embedding
    <timepoint>/<Gene>/isp_stats.csv             # 各遺伝子の raw 統計
```

`<experiment>` は **`<disease>_<tissue>` の小文字**（`ad_spleen`, `pd_spleen`, `ad_blood`, `ad_smallint`, `ad_brain`, `ad_liver`, `ad_bm`）。

| 旧パス（PC 別） | 新パス（統一） |
|---|---|
| `results/isp/early_ad` | `results/isp/ad_smallint` |
| `results/isp/early_ad_microglia` | `results/isp/ad_brain` |
| `results/isp/early_ad_blood` | `results/isp/ad_blood` |
| `results/isp/early_ad_spleen` | `results/isp/ad_spleen` |
| `results/isp/early_pd` | `results/isp/pd_spleen` |
| `results/isp/` (直置き) | `results/isp/ad_liver` |
| `results/isp/` (直置き) | `results/isp/ad_bm` |

**対応スクリプト** `analysis/_isp_common.py` の `resolve_experiment()` / `isp_dir_for()` / `combined_csv_path()` がこの規則を返します。

---

## 3. 共通リファクタ（`analysis/_isp_common.py`）

各臓器スクリプトに重複していた以下を集約しました。

| 関数 | 役割 |
|---|---|
| `warn_if_multiproc()` | nproc>1 の警告（datasets<5 前提） |
| `check_datasets_version()` | datasets>=5 で `perturb_data` ハング警告 |
| `estimate_perturb_ram()` | `max_ncells` の RAM 見積り |
| `load_gene_dicts()` | token→ENSG / ENSG→symbol / symbol→ENSG 辞書 |
| `resolve_classifier_dir()` | ファインチューン済みモデル解決（stale パス / 最新 ksplit fallback） |
| `scan_pool_presence()` / `_per_key()` | 対象プールでの遺伝子発現チェック |
| `perturb_one_gene()` | 1 遺伝子削除 + `goal_state_shift` 統計（`Shift_to_goal_end` 返却） |
| `write_combined()` | 全遺伝子結果の結合 CSV 出力（early タイムポイント順位付き） |

---

## 4. メソッド共通点

- **state_key = `disease`**, start = `AD`（または `PF`）, goal = `WT`
- **perturb_type = `delete`**, `combos=0`, `emb_mode="cls"`, `model_version="V2"`
- **1 遺伝子ずつ**削除（combos=0 でリストを渡すと全遺伝子同時保持条件になり細胞が空になるため）
- 状態 embedding はタイムポイントごとに 1 回計算し、全遺伝子で共有
- `nproc=1`（veri-stable）, `datasets==4.0.0`（>=5 は `perturb_data` がハング）

---

## 5. ページ一覧

- [方法・共通設定（詳細）](methods.md)
- [実行結果のまとめ](results.md)
- [肝臓（ad_liver）](ad_liver.md)
- [骨髄（ad_bm）](ad_bm.md)
- [PD の現状・残タスク（pd_spleen）](pd.md)
- [再現手順（How to run）](how_to_run.md)

各臓器のスクリプト設定（細胞プール・遺伝子リスト・タイムポイント）と結果トップは各ページに記載しています。

---

## 6. 結果データの場所

統一後の結果 CSV（実計算値）は Git にコミットせず、ローカルの
`input/<TISSUE>/results/isp/<experiment>/<experiment>_early_isp_stats_combined.csv`
にあります（`.gitignore` 対象）。

レポート用に整形したデータは `result/report/data/` に、図は `result/report/figures/` にあります（いずれも計算成果物のため Git 対象外）。実行スクリプト・Wiki・レポートテンプレートは Git 追跡対象です。