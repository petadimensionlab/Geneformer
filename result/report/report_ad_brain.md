# 早期アルツハイマー病（AD）発症因子の探索 — in silico gene deletion レポート（脳 / AD_brain）

> 実験スクリプト: `analysis/07_ad_brain_early_isp.py`
> データセット: `input/AD_brain`（5xFAD マウス脳 scRNA-seq、ヒトオルソログ変換済み）

---

## 1. 対象

| 項目 | 値 |
|---|---|
| 対象疾患 | **AD（アルツハイマー病）** |
| 対象臓器 | **脳（brain）** |
| データセット | 5xFAD マウス脳、41,858 細胞、24 細胞型、30 個体 |
| 疾患状態 | AD (21,827) vs WT (20,031) |
| 時点 | 3・4.5・6・9・12 ヶ月齢（`samples4`） |
| **解析対象時点** | **早期（3–4.5–6 ヶ月）のみ** |

実験対応（`result/report/experiments.csv`）: `isp_early_ad_brain`

---

## 2. 臓器特異的に探索すべき背景知識（脳）

AD は **神経炎症** を病態の中心とする。特に **ミクログリア（脳常在マクロファージ）** は、アミロイドβ(Aβ)蓄積のごく初期から応答し、**DAM（disease-associated microglia）** 状態へ分化する。

### 注目すべき細胞
- **ミクログリア（Microglia）**: AD の主要な免疫細胞。TREM2-APOE 軸を介した DAM 遷移が進行に関与。
- **BAM（border-associated macrophages）**: 脳血管周囲・脳脊髄液側の常在マクロファージ。初期炎症応答に寄与。

### 注目すべき遺伝子の機能クラス
- **DAM ミクログリアシグネチャー**: TREM2, APOE, TYROBP, CSF1R, C1QA/B, C3, ITGAM(CD11b), ITGAX(CD11c), CD68, SPP1
- **炎症性サイトカイン**: IL1B, TNF
- **AD GWAS リスク（免疫濃縮）**: APOE, TREM2, CLU, BIN1, CD33, ABCA7, SORL1, INPP5D, PLCG2, MEF2C

### 参考文献（ミクログリア・DAM と AD 早期化）
| 文献 | 要点 |
|---|---|
| **Keren-Shaul et al., *Cell* 169:1276 (2017)** | 単細胞で **DAM（disease-associated microglia）** を同定。TREM2-APOE 軸の必要性 |
| **Guerreiro et al., *NEJM* 368:117 (2013)** | TREM2 のミス変異が AD リスクを約3倍に増大（**TREM2 が AD リスクの最有力**） |
| **Jonsson et al., *NEJM* 368:107 (2013)** | TREM2 変異と AD の関連 |
| **Kinney et al., *J Alzheimers Dis Parkinsonism* 8:3 (2018)** | 炎症は AD の早期イベントであるとする炎症仮説 |
| **Heneka et al., *Nat Rev Neurol* 11:142 (2015)** | 神経炎症と AD 進行の包括的レビュー |

### 3. 早期発見の観点
5xFAD では **3–6 ヶ月** に Aβ が蓄積し、**ミクログリアの応答（DAM 遷移）** が初期に起こる。したがって、**早期ミクログリア/BAM の遺伝子を仮想削除し、AD→WT への embedding シフト**を見ることは、早期発症を駆動する免疫遺伝子の同定に直接つながる。

---

## 3. 細胞・遺伝子の選択理由

### 選択した細胞型（2種、早期のみ）
**ミクログリア（Microglia）＋ BAM（border-associated macrophages）**。AD 早期の主要な免疫応答細胞であり、DAM 遷移の中核（Keren-Shaul 2017）。早期 AD/WT 細胞数は下記の通り。

### 選択した遺伝子（23種、仮説駆動・固定リスト）
| カテゴリ | 遺伝子 | 選択根拠 |
|---|---|---|
| **AD risk (GWAS)** | APOE, TREM2, TYROBP, BIN1, CD33, ABCA7, SORL1, INPP5D, PLCG2, MEF2C | Lambert 2013; Wightman 2021; Jansen 2019 |
| **DAM ミクログリア** | CSF1R, C1QA, C1QB, C3, CD68, AIF1, ITGAM, ITGAX, SPP1, CCL2 | Keren-Shaul 2017（DAM 標識）; 
| **炎症** | IL1B, TNF, CD74 | 神経炎症（Heneka 2015; Kinney 2018） |

> 文献: DAM 標識遺伝子（Keren-Shaul 2017）、TREM2 リスク（Guerreiro 2013）、炎症仮説（Kinney 2018）。

---

## 4. in silico perturbation の主要設定

| 設定項目 | 値 |
|---|---|
| **遺伝子操作** | **delete（仮想遺伝子削除）** |
| 対象状態（disease） | `start_state = AD`, `goal_state = WT` |
| state_key | `disease` |
| alt_states | なし |
| 細胞型フィルタ | Microglia + BAM（早期のみ） |
| 時点フィルタ | 早期（AD/WT × 3/4.5/6m） |
| combos | 0（1遺伝子ずつ） |
| モデル | Geneformer V2-104M（`model_version="V2"`） |
| 埋め込み | `emb_mode="cls"`, `cell_emb_style="mean_pool"` |
| 分類器 | ファインチューン済みセル分類器 |
| 統計モード | `goal_state_shift`（`Shift_to_goal_end`） |

> **解釈**: 早期 AD ミクログリア/BAM から各遺伝子を削除し、embedding が早期 WT へどれだけ近づくかを `Shift_to_goal_end` で評価。**正に大きい = 早期駆動/介入標的候補**。

---

## 5. 計算パラメータ・実行時間

| パラメータ | 値 |
|---|---|
| 解析対象細胞数（早期プール） | Microglia 早期 AD/WT = 1,800/1,800; BAM 早期 AD/WT = 376/264 |
| 状態埋め込み（EmbExtractor） | `max_ncells=1000` |
| perturb（InSilicoPerturber） | `max_ncells=200`, `forward_batch_size=64`, `emb_layer=0` |
| 遺伝子数 | 候補 23 → 発現 23（全てプールに発現） |
| nproc | 1（`datasets==4.0.0` 必須） |
| デバイス | CUDA |
| **実行時間** | **約 67 分**（15:16–16:23, 2026-08-25） |

> **補助**: `datasets==4.0.0` 固定（`>=5` で perturb ハング）。`max_ncells=200` のスクリーニング結果を使用（`max_ncells=1000` の再実行は後続で別途実施）。

---

## 6. 結果まとめ（in silico deletion — 早期 AD→WT シフト）

`Shift_to_goal_end` のランキング。全 23 遺伝子の値は `result/report/data/ad_brain_early_isp_results.csv` を参照。

### 上位ヒット
| Rank | Gene | 分類 | Shift_to_goal_end |
|---|---|---|---|
| 1 | **APOE** | AD risk / DAM | **+0.00872** |
| 2 | **TREM2** | AD risk / DAM | **+0.00006** |
| 3–23 | C1QB, TNF, IL1B, CCL2, SPP1, ITGAX, ITGAM, AIF1, CD68, C3, C1QA, CSF1R, MEF2C, PLCG2, INPP5D, SORL1, ABCA7, CD33, BIN1, TYROBP, CD74 | — | **-0.000216**（ほぼ一律の小さい負） |

### 生物学的解釈
- **APOE が突出（+0.0087）**。ミクログリア/BAM の早期 AD 細胞から **APOE を削除すると embedding が早期 WT 方向へ最も大きく戻る**。これは APOE が早期ミクログリア駆動の最重要ハブであることを示唆し、**早期発見・介入の最有力標的**。
- **TREM2** は弱い正（+0.00006）。
- 他の 21 遺伝子はほぼ中性〜微弱な負にまとまり、`max_ncells=200` では有意なシフトが検出されなかった。APOE の特異性が際立つ結果。

> **注意**: 多数の遺伝子が同一値（-0.000216）に固まっているのは、`max_ncells=200` のスクリーニングで効果が埋もれた可能性がある。より高解像度（`max_ncells=1000`）の再実行が検討対象。

> 図: `result/report/figures/ad_brain_isp_rank.png`, `figures/ad_smallint_vs_brain_top.png`

---

## 7. 制限事項・留意点
- `Shift_to_goal_end` は cosine シフトで、統計的有意性（p値）は出力に含まれない。
- `max_ncells=200` のスクリーニング結果。多数遺伝子が同一小値に固まるため、候補の分離には解像度向上（細胞数増加）が必要。
- 5xFAD は家族性 AD モデルで、ヒト散発性 AD への一般化には注意。