# 実行結果のまとめ（In Silico Perturbation）

> 各実験の `Shift_to_goal_end` 上位（正 = 疾患細胞が WT に近づく = 早期駆動因子/介入標的候補）をまとめます。
> 全結果は `input/<TISSUE>/results/isp/<experiment>/<experiment>_early_isp_stats_combined.csv`。

---

## 一覧

| 実験 | 疾患 | 臓器 | 遺伝子数 | タイムポイント | モデル | 実行時間(目安) |
|---|---|---|---|---|---|---|
| ad_smallint | AD | 小腸 | 49 | early(3-6m pooled) | Fine-tuned | ~45分 |
| ad_brain | AD | 脳(ミクログリア) | 23 | early(3-6m pooled) | Fine-tuned | ~67分 |
| ad_blood | AD | 血液 | 55 | 3m/4.5m/6m/9m/12m | Fine-tuned | ~1h35m |
| ad_spleen | AD | 脾臓 | 55 | 3m/4.5m/6m | Fine-tuned | ~50分 |
| pd_spleen | PD | 脾臓 | 34 | 6m/9m/12m | **Fine-tuned** | ~1h(finetune)+ISP |
| ad_liver | AD | 肝臓 | 21 | 3m→後期 | Fine-tuned | - |
| ad_bm | AD | 骨髄 | 11 | AD→WT | Fine-tuned | - |

> 実測値（Shift_to_goal_end）はローカル成果物のため数値は下表のとおり（生成時点）。最新値は各 `*_combined.csv` を参照。

---

## AD（アルツハイマー病）

### 小腸 `ad_smallint` — 腸-脳軸（gut-brain axis）

**細胞プール**: Macrophages, CD4/CD8 T, NK_ILC1, IECs, Paneth, Goblet, Epithelial（免疫 + バリア細胞）
**遺伝子**: AD risk + 炎症 + 粘液（MUC）+ タイトジャンクション + 抗菌ペプチド + 幹細胞

| Gene | Shift_to_goal_end |
|---|---|
| **CLDN1**（タイトジャンクション） | **+0.0093** |
| MUC5B（粘液） | +0.0038 |
| TREM2 | +0.0029 |
| SORL1 | +0.0023 |
| MUC2 | +0.0018 |
| TYROBP | +0.0016 |

**解釈**: 腸バリア整合性（CLDN1）が早期 AD で最も disease-driver。粘液（MUC）と免疫/ミクログリア（TREM2）も正。

### 脳 `ad_brain`（ミクログリア / BAM）

| プール | 遺伝子数 | 上位（Shift） |
|---|---|---|
| Microglia only | 23 | **APOE +0.0007**（他 ~0） |
| Microglia + BAM | 23 | 全遺伝子 ~0（APOE も微小） |

**解釈**: ミクログリアのみでは **APOE** が唯一の正シフト（早期ミクログリアの最有力ターゲット）。BAM を含めるとシグナルが薄まる（max_ncells=200 の限界）。

### 血液 `ad_blood`（末梢免疫）

| Gene | 3m Shift |
|---|---|
| **CD3E**（TCR コア） | **+0.0142** |
| **CD4**（ヘルパーT） | **+0.0138** |
| CD74（MHC-II） | +0.0053 |
| S100A8（警報因子） | +0.0050 |
| C1QC | +0.0049 |
| CD19 | +0.0041 |

**解釈**: **T 細胞軸（CD3E/CD4）が最有力**。血液では APOE ほぼ中性（脳と対照）。CD8A は強く負（-0.11）= CD8 細胞毒性は WT 維持に寄与。

### 脾臓 `ad_spleen`（免疫）

| Gene | 3m Shift |
|---|---|
| **CLEC9A**（cDC1） | **+0.0066** |
| **ITGA2B**（血小板） | +0.0062 |
| CD19（B） | +0.0034 |
| CD4 | +0.0034 |
| CLU | +0.0029 |
| CD33 | +0.0014 |

**解釈**: cDC1（クロス提示）・血小板（ITGA2B）・B 細胞（CD19）軸が最早期に disease を駆動。APOE 中性。CD8A 強く負。

### 肝臓 `ad_liver`

> 早期 AD 肝臓の免疫/炎症軸（Kupffer, 単球, マクロファージ, 肝細胞）。WT 3m baseline → AD3m goal、AD6m/9m/12m を alt として 進行方向のシフトを評価。

### 骨髄 `ad_bm`

> 骨髄骨髄系（Ly6c.high Monocytes / Macrophage）。11 個の AD/炎症候補遺伝子を 1 遺伝子ずつ削除。

---

## PD（パーキンソン病）

### 脾臓 `pd_spleen`（免疫）— **Fine-tuned CellClassifier（2026-08-26 再実行）**

| Gene | 6m Shift |
|---|---|
| **S100A8**（カルプロテクチン） | **+0.0356** |
| **S100A9** | +0.0221 |
| **LYZ** | +0.0118 |
| ITGAX | +0.0028 |
| C1QA | +0.0020 |
| TREM1 | +0.0019 |
| IL6 | +0.0016 |
| ITGAM | +0.0015 |
| SNCA | +0.0015 |

**解釈**: PD 脾臓では **S100A8/S100A9（警報因子）** と **LYZ（ライソザイム）** が最有力。α-synuclein 注入（PFF）モデルで、神経炎症の末梢プロキシ。SNCA 自体も正（+0.0015）で、α-syn 軸が脾臓免疫の疾病状態に関与。

> **2026-08-26 更新**: `06_finetune.py` で fine-tune 完了（accuracy 0.9149 / macro F1 0.9079）後、`07_pd_spleen_early_isp.py` を **Pretrained → CellClassifier** に変更して再実行。S100A8/S100A9 は Pretrained 時（+0.0009/+0.0006）から大幅に増強（+0.036/+0.022）。Pretrained 時の FOXP3/C1QB/CD14 上位は fine-tuned ではランク外（モデル差）。詳細は **[pd.md](pd.md)** を参照。

---

## 補足：血液 vs 脾臓（AD）の比較

| Gene | 血中 3m | 脾臓 3m | 解釈 |
|---|---|---|---|
| CD3E | +0.0142 | +0.0007 | 血中では T 細胞が最有力 |
| CLEC9A | +0.0022 | +0.0066 | 両方で正 → cDC1 一貫 |
| ITGA2B | +0.0028 | +0.0062 | 血小板系一貫 |
| S100A8 | +0.0050 | -0.0134 | 符号が逆（血=正, 脾=負）|
| CD8A | -0.112 | -0.098 | 両方で強く負 |

## 制限事項
- `Shift_to_goal_end` は cosine シフト（統計的有意性 p 値なし）
- `max_ncells=200-300` のスクリーニング（上位候補の分離にはより多細胞で再実行が有効）
- 5xFAD（家族性AD）/ PFF（散発性PD）モデル由来。ヒト散発性疾患への一般化には注意