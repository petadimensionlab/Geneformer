# 実行結果のまとめ（In Silico Perturbation）

> Shift の値は有効数字 3 桁の指数表記で書いています（例: +6.60e-03）。割合は小数 1 桁です。

> 各実験の `Shift_to_goal_end` 上位（正 = 疾患細胞が WT に近づく = 早期駆動因子/介入標的候補）をまとめます。
> 全結果は `input/<TISSUE>/results/isp/<experiment>/<experiment>_early_isp_stats_combined.csv`。

---

## 一覧

| 実験 | 疾患 | 臓器 | 遺伝子数 | タイムポイント | モデル | 実行時間(目安) |
|---|---|---|---|---|---|---|
| ad_smallint | AD | 小腸 | 49 | early(3-6m pooled) | 微調整済み | ~45分 |
| ad_brain | AD | 脳(ミクログリア) | 23 | early(3-6m pooled) | 微調整済み | ~67分 |
| ad_blood | AD | 血液 | 55 | 3m/4.5m/6m/9m/12m | 微調整済み | ~1h35m |
| ad_spleen | AD | 脾臓 | 55 | 3m/4.5m/6m | 微調整済み | ~50分 |
| pd_spleen | PD | 脾臓 | 34 | 6m/9m/12m | **微調整済み** | ~1h(finetune)+ISP |
| ad_liver | AD | 肝臓 | 21 | 3m→後期 | 微調整済み | - |
| ad_bm | AD | 骨髄 | 11 | AD→WT | 微調整済み | - |

> 実測値（Shift_to_goal_end）はローカル成果物のため数値は下表のとおり（生成時点）。最新値は各 `*_combined.csv` を参照。

---

## AD（アルツハイマー病）

### 小腸 `ad_smallint` — 腸-脳軸（gut-brain axis）

**細胞プール**: Macrophages, CD4/CD8 T, NK_ILC1, IECs, Paneth, Goblet, Epithelial（免疫 + バリア細胞）
**遺伝子**: AD risk + 炎症 + 粘液（MUC）+ タイトジャンクション + 抗菌ペプチド + 幹細胞

| Gene | Shift_to_goal_end |
|---|---|
| **CLDN1**（タイトジャンクション） | **+9.30e-03** |
| MUC5B（粘液） | +3.80e-03 |
| TREM2 | +2.90e-03 |
| SORL1 | +2.30e-03 |
| MUC2 | +1.80e-03 |
| TYROBP | +1.60e-03 |

**解釈**: 腸バリア整合性（CLDN1）が早期 AD で最も disease-driver。粘液（MUC）と免疫/ミクログリア（TREM2）も正。

### 脳 `ad_brain`（ミクログリア / BAM）

| プール | 遺伝子数 | 上位（Shift） |
|---|---|---|
| Microglia only | 23 | **APOE +7.00e-04**（他 ~0） |
| Microglia + BAM | 23 | 全遺伝子 ~0（APOE も微小） |

**解釈**: ミクログリアのみでは **APOE** が唯一の正の Shift（早期ミクログリアの最有力ターゲット）。BAM を含めるとシグナルが薄まる（max_ncells=200 の限界）。

### 血液 `ad_blood`（末梢免疫）

| Gene | 3m Shift |
|---|---|
| **CD3E**（TCR コア） | **+1.42e-02** |
| **CD4**（ヘルパーT） | **+1.38e-02** |
| CD74（MHC-II） | +5.30e-03 |
| S100A8（警報因子） | +5.00e-03 |
| C1QC | +4.90e-03 |
| CD19 | +4.10e-03 |

**解釈**: **T 細胞軸（CD3E/CD4）が最有力**。血液では APOE ほぼ中性（脳と対照）。CD8A は強く負（-0.11）= CD8 細胞毒性は WT 維持に寄与。

### 脾臓 `ad_spleen`（免疫）

| Gene | 3m Shift |
|---|---|
| **CLEC9A**（cDC1） | **+6.60e-03** |
| **ITGA2B**（血小板） | +6.20e-03 |
| CD19（B） | +3.40e-03 |
| CD4 | +3.40e-03 |
| CLU | +2.90e-03 |
| CD33 | +1.40e-03 |

**解釈**: cDC1（クロス提示）・血小板（ITGA2B）・B 細胞（CD19）軸が最早期に disease を駆動。APOE 中性。CD8A 強く負。

> **モデル依存について（2026-09 追試）**: 上の表は 104M の結果です。同じデータで 316M に替えると順位は動きます。
> 55 遺伝子のうち 18 個で Shift の符号が反転し、差の大きさは計算誤差の約 7 倍、細胞の選び方を変えた差の約 5 倍でした。
> 上の主要 6 遺伝子（CLEC9A、ITGA2B、CD19、CD4、CLU、CD33）は両モデルで正のままですが、値は小さくなり順位も入れ替わります
> （例: CD19 は 3 位 → 16 位、CLEC9A は両方とも 1 位）。CD8A と APOE は両モデルで負のままです。
> 一方、TYROBP・GCA・LYZ は符号が反転します。したがって上の順位は暫定として扱い、
> 標的候補は「両モデルで同じ向きが出た遺伝子」に限ってください（[遵守項目](checklist.md) B 章）。

### 肝臓 `ad_liver`

> 早期 AD 肝臓の免疫/炎症軸（Kupffer, 単球, マクロファージ, 肝細胞）。WT 3m baseline → AD3m goal、AD6m/9m/12m を alt として 進行方向のシフトを評価。

### 骨髄 `ad_bm`

> 骨髄骨髄系（Ly6c.high Monocytes / Macrophage）。11 個の AD/炎症候補遺伝子を 1 遺伝子ずつ削除。

---

## PD（パーキンソン病）

### 脾臓 `pd_spleen`（免疫）— **微調整済み CellClassifier（2026-08-26 再実行）**

| Gene | 6m Shift |
|---|---|
| **S100A8**（カルプロテクチン） | **+3.56e-02** |
| **S100A9** | +2.21e-02 |
| **LYZ** | +1.18e-02 |
| ITGAX | +2.80e-03 |
| C1QA | +2.00e-03 |
| TREM1 | +1.90e-03 |
| IL6 | +1.60e-03 |
| ITGAM | +1.50e-03 |
| SNCA | +1.50e-03 |

**解釈**: PD 脾臓では **S100A8/S100A9（警報因子）** と **LYZ（ライソザイム）** が最有力。α-synuclein 注入（PFF）モデルで、神経炎症の末梢プロキシ。SNCA 自体も正（+1.50e-03）で、α-syn 軸が脾臓免疫の疾病状態に関与。

> **2026-08-26 更新**: `06_finetune.py` で微調整が完了（accuracy 0.9149 / macro F1 0.9079）後、`07_pd_spleen_early_isp.py` を **Pretrained → CellClassifier** に変更して再実行。S100A8/S100A9 は Pretrained 時（+9.00e-04/+6.00e-04）から大幅に増強（+3.60e-02/+2.20e-02）。Pretrained 時の FOXP3/C1QB/CD14 上位は微調整後ではランク外（モデル差）。詳細は **[pd.md](pd.md)** を参照。

---

## 補足：血液 vs 脾臓（AD）の比較

| Gene | 血中 3m | 脾臓 3m | 解釈 |
|---|---|---|---|
| CD3E | +1.42e-02 | +7.00e-04 | 血中では T 細胞が最有力 |
| CLEC9A | +2.20e-03 | +6.60e-03 | 両方で正 → cDC1 一貫 |
| ITGA2B | +2.80e-03 | +6.20e-03 | 血小板系一貫 |
| S100A8 | +5.00e-03 | -1.34e-02 | 符号が反転（血=正、脾=負）|
| CD8A | -0.112 | -9.80e-02 | 両方で強く負 |

## 制限事項
- `Shift_to_goal_end` はコサイン類似度の差（統計的有意性の p 値は出せない）
- `max_ncells=200-300` のスクリーニング（上位候補の分離にはより多細胞で再実行が有効）
- 5xFAD（家族性AD）/ PFF（散発性PD）モデル由来。ヒト散発性疾患への一般化には注意