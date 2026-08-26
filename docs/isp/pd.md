# PD（パーキンソン病）in silico perturbation — 現状と整理

> PD の in silico perturbation は **AD に比べて完成度が大幅に低く**、リポジトリ上に
> **データ・モデル・スクリプト・結果のすべてが「1 臓器分」しか揃っていません**。
> ここでは現状を正確に記録し、AD と揃えるための残タスクを整理します。

---

## 1. 現状サマリ（PD vs AD）

| 項目 | AD | **PD** |
|---|---|---|
| データ臓器 | blood, brain, smallint, spleen, liver, BM（+LN 未取得） | **spleen のみ**（script は LN/BM/blood/brain/smallint 分あり） |
| h5ad / tokenized | 6 臓器完備 | PD_spleen のみ（h5ad 667MB, 83,783 細胞, 24 細胞型, 36 個体） |
| fine-tune 済み分類器 | 実重みあり | **未完成（0 バイト）** |
| ISP 結果 | 複数臓器 | PD_spleen 1 件のみ（85 レコード） |
| 専用スクリプト | 5 + liver/bm | `07_pd_spleen_early_isp.py` のみ |
| レポート | 各臓器 .md/.html | `report_pd_spleen.*` のみ |

---

## 2. 判明している不足（具体的に）

### 2-1. データ臓器が spleen のみ
- `input/` 配下の PD は **`PD_spleen` のみ**。
- ユーザ提示の全臓器リスト（blood / brain / smallint / spleen / BM / liver / LN）は
  **AD 用**であり、PD は元から「脾臓 1 臓器」しか入力されていない。

### 2-2. fine-tune 済み分類器が未完成（0 バイト）
- `input/PD_spleen/runs/260823_geneformer_cellClassifier_PD_spleen_celltype/` は **8–12 KB**。
  - `ksplit1/config.json`, `ksplit1/model.safetensors` が **すべて 0 バイト**
  - `TRAINED_MODEL_PATH.txt`, `*_pred_dict.pkl` 等も 0 バイト
- そのため `07_pd_spleen_early_isp.py` は **`model_type="Pretrained"`（未ファインチューン）** で実行。
  → AD 系の `CellClassifier`（fine-tuned）と非対称で、組織特異シグナルの抽出精度に差がある。
- `results/tables/adpd_*`（fine-tuned 評価）も **空**（0 バイト）＝ fine-tune 結果が存在しない。

### 2-3. PD スクリプトの存在状況（2026-08-26 更新）

> **更新**: `analysis/07b_in_silico_perturbation_PD.py`（PD_brain）と
> `analysis/07e_in_silico_perturbation_PD_smallint.py`（PD_smallint）が
> ローカルに実在することを確認し、リポジトリへ**取り込み commit** 済み。
> あわせて `analysis/07c_in_silico_perturbation_PD_combos.py`（PD_brain の
> SNCA×partner 削除）と `analysis/07d_in_silico_perturbation_PD_null.py`
> （PD_brain の null 有意性）も取り込み済み。

| スクリプト | 対象 | 存在（更新後） |
|---|---|---|
| `07b_in_silico_perturbation_PD.py` | PD_brain | ✅ 取込済 |
| `07e_in_silico_perturbation_PD_smallint.py` | PD_smallint | ✅ 取込済 |
| `07c_in_silico_perturbation_PD_combos.py` | PD_brain（combos） | ✅ 取込済 |
| `07d_in_silico_perturbation_PD_null.py` | PD_brain（null） | ✅ 取込済 |
| `07_pd_spleen_early_isp.py` | PD_spleen | ✅ 追跡済 |
| `07d_pd_blood_perturbation.py` | PD_blood | ✅ 追跡済 |
| `07c_pd_bm_perturbation.py` | PD_BM | ✅ 追跡済 |
| `07b_pd_perturbation.py` | PD_LN | ✅ 追跡済 |

> 注: 全リスト **LN（リンパ節）** は、AD / PD いずれも **データ・レポートが存在しない**。
> - **AD_LN**: スクリプト `07_ad_ln_early_isp.py`（CD4+ T 小規模テスト）は **2026-08-26 追跡済**。データ（`input/AD_LN`）は未配置。
> - **PD_LN**: `07b_pd_perturbation.py` が命名上「LN」を対象にしているが、実データ（`input/PD_LN` / h5ad）は
>   リポジトリ上に確認できない。

### 2-4. PD 用の Wiki ページが無い
- 現在 PD の記述は `docs/isp/README.md`（1 行）と `results.md`（1 節）のみ。
- 本ページ（`docs/isp/pd.md`）を新設し、現状と残タスクを記録。

---

## 3. 現在までに実施済みの PD 解析

### PD_spleen（`07_pd_spleen_early_isp.py`, **Pretrained**）
- **state**: `disease` PF(α-syn PFF 注入 = 孤発性 PD) → WT
- **タイムポイント**: 6m / 9m / 12m を個別評価、6m でランキング
- **細胞プール**: 脾臓免疫 16 種
- **遺伝子**: 34（α-syn/PD core + リソソーム + 神経炎症/免疫 + 補体）
- **結果**: 85 レコード（`results/isp/pd_spleen/pd_spleen_early_isp_stats_combined.csv`）
- **実行時間**: 約 30 分（2026-08-25）
- **上位（6m Shift）**: S100A8 / S100A9 / FOXP3 / C1QB / CD14 / CXCL10

---

## 4. 残タスク（AD と揃えるために）

### 優先度 A（書類整理、データ不要）
- [x] Wiki に `docs/isp/pd.md` を新設し現状を記録
- [x] `experiments.csv` に PD の完成度・欠落を注記
- [x] 他 PD スクリプト（PD_brain / PD_smallint / combos / null）をローカルから取り込み push

### 優先度 B（データ構築、実行/別 PC 必要）
- [ ] PD_spleen の **fine-tune 再実行**（`06_finetune.py`）→ 実モデル重み生成
- [ ] 他 PD 臓器（blood / brain / LN）データの取得・構築
- [ ] fine-tune 済みモデルで `07_*` を `CellClassifier` で再実行（AD と対称に）

### 注意（解釈）
- PD は現状 **Pretrained ベース**。fine-tune 完了までは、AD の fine-tuned 結果と直接比較しないこと。
- 脾臓単一臓器・PFF モデルはヒト PD の全容を再現しない。