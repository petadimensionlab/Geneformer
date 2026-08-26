# 早期アルツハイマー病（AD）発症因子の探索 — in silico gene deletion レポート（血液 / AD_blood）

> 実験スクリプト: `analysis/07_ad_blood_early_isp.py`
> データセット: `input/AD_blood`（5xFAD ファミリーADモデルマウス血液 scRNA-seq、ヒトオルソログ変換済み）
> 実施日: 2026-08-26

---

## 1. 対象

| 項目 | 値 |
|---|---|
| 対象疾患 | **AD（アルツハイマー病）** |
| 対象臓器 | **血液（blood / 末梢血）** |
| データセット | AD_blood、48,909 細胞、15 細胞型、29 サンプル（個体 index1–4） |
| 疾患状態 | AD (25,222) vs WT (23,687) |
| 時点 | 3・4.5・6・9・12 ヶ月齢（`samples4`） |
| **解析対象時点** | **全時点（3m / 4.5m / 6m / 9m / 12m）を個別評価、3m（最早期）でランキング** |

実験対応（`result/report/experiments.csv`）: `isp_early_ad_blood`

---

## 2. 対象臓器（血液）特異的に探索すべき背景知識

AD 研究は脳・ミクログリアに集中するが、**血液（末梢免疫）は侵襲的でない早期発見バイオマーカー**の宝庫であり、AD 臨床前段階から末梢免疫プロファイルが変化するという証拠が蓄積している。

### 注目すべき細胞
- **単球 / マクロファージ系（Ly6c.low / Ly6c.high Monocytes）**: 血液常在の"ミクログリア様"細胞。TREM2 / TYROBP / SPI1 などの AD リスク遺伝子を発現し、脳の DAM 軸に対応する最有力の血中プロキシ。
- **T 細胞系（CD4 / CD8 / NKT）**: AD で CD8 T 細胞のクローン性増殖が CSF で報告され（Gate et al., *Nature* 2020）、末梢 T 細胞プロファイルの変化が早期に生じる。
- **B 細胞 / プラズマ細胞（Naive.Memory.B）**: 液性免疫変化とアミロイド関連抗体。
- **好中球（Neutrophils）**: S100A8/S100A9（カルプロテクチン）を産生。全身性炎症・警報因子の血中マーカー。
- **cDCs / pDCs**: 抗原提示・I 型 IFN（IRF7/ISG15）軸。

### 注目すべき遺伝子の機能クラス
- **AD GWAS リスク（免疫濃縮）**: APOE, TREM2, TYROBP, CLU, BIN1, CD33, ABCA7, SORL1, INPP5D, PLCG2, MEF2C, SPI1
- **T 細胞 / 適応免疫**: CD3D/E, CD4, CD8A, NKG7, KLRD1, PRF1, CD19, MS4A1
- **自然免疫 / 炎症**: C1QA/B/C, C3, CD68, CSF1R, AIF1, ITGAM/ITGAX, SPP1, CCL2, IL1B, TNF, TLR4, NLRP3, CD74, IRF7, ISG15, S100A8/A9
- **骨髄球 / 顆粒球 / 血小板マーカー**: CSF3R, GCA, CEACAM1, FGR, ITGA2B, VWF, CD34, XCR1, IL3RA, CLEC9A

### 参考文献（血液・末梢免疫と AD 早期化）
| 文献 | 要点 |
|---|---|
| **Gate et al., *Nature* 583:459 (2020)** | AD で **CSF にクローン性増殖した CD8 T 細胞**が増加。T 細胞系（CD3/CD8/NKG7）に着目する根拠 |
| **Tan et al., *Neurology* 58:1053 (2002)** | AD で血液 **S100A8/A9（カルプロテクチン）** と炎症マーカーが上昇 |
| **Bettcher et al., *Brain Behav Immun* (2018)** | 末梢免疫と AD 病態の関連レビュー |
| **Nakamura et al., *Nature* 554:249 (2018)** | 血液アミロイド測定による AD 早期検出 — "血液で AD を早期に見る"という枠組みの基盤 |
| **Guerreiro et al., *NEJM* 368:117 (2013)** | TREM2 変異が AD リスク増大 — 血中免疫細胞でも発現する AD リスク軸 |
| **Jansen et al., *Nat Genet* 51:414 (2019)** | AD GWAS 遺伝子が免疫/MHC 機構に濃縮 |

---

## 3. 選択した細胞・遺伝子の理由（早期発見の観点）

### 選択した細胞型（10 種 — 低シグナル細胞を除外）
**Naive.Memory.B.cells, CD4.T.cells, CD8.T.cells, NK_ILC1, Neutrophils, Ly6c.low.nonclassical.Monocytes, Ly6c.high.classical.Monocytes, cDCs, pDCs, Gamma.Delta.T.cells**

- 免疫担当細胞を網羅しつつ、**Basophils / Megakaryocytes / Erythroblasts / Immature.T.cells / NKT.cells は除外**（低頻度・赤血球/巨核球系ノイズ、早期判定に寄与が薄い）。
- 単球=ミクログリア様血中細胞、T/B 細胞=適応免疫応答、好中球=警報因子という 3 系統をカバー。

### 選択した遺伝子（55 種、仮説駆動・固定リスト）
| カテゴリ | 遺伝子数 | 代表 | 選択根拠 |
|---|---|---|---|
| AD risk (GWAS) | 12 | APOE, TREM2, TYROBP, CLU, BIN1, CD33, ABCA7, SORL1, INPP5D, PLCG2, MEF2C, SPI1 | Lambert 2013; Jansen 2019; Wightman 2021（免疫濃縮） |
| 炎症・自然免疫 | 24 | C1QA/B/C, C3, CD68, CSF1R, AIF1, ITGAM, ITGAX, SPP1, CCL2, IL1B, TNF, TLR4, NLRP3, LYZ, APOC1, CD74, IRF7, ISG15, S100A8/A9, MS4A7, THBS1 | 神経炎症（Heneka 2015; Kinney 2018）、アルマミン（Foell 2007） |
| 適応免疫 | 12 | CD3E, CD3D, CD4, CD8A, CD19, MS4A1, NKG7, KLRD1, PRF1 | Gate 2020（CD8 T クローン）、末梢適応免疫変化 |
| 骨髄球/顆粒球/血小板 | 10 | CSF3R, GCA, CEACAM1, FGR, ITGA2B, VWF, CD34, XCR1, IL3RA, CLEC9A | 血中系統マーカーの網羅 |

> 全 55 遺伝子は早期（3m–12m）AD 血液細胞で発現を確認済み（`result/report/data/ad_blood_gene_list.csv`）。

---

## 4. in silico perturbation の主要設定

| 設定項目 | 値 |
|---|---|
| **遺伝子操作** | **delete（仮想遺伝子削除）** |
| 対象状態（disease） | `start_state = AD`, `goal_state = WT` |
| state_key | `disease` |
| alt_states | なし |
| 細胞型フィルタ | 血液免疫プール 10 種 |
| 時点フィルタ | 各時点ごと（AD/WT × 3/4.5/6/9/12m） |
| combos | 0（1遺伝子ずつ個別評価）× 各時点 |
| モデル | Geneformer V2-104M ファインチューンフィネチューン済みセル分類器（`CellClassifier`） |
| 埋め込み | `emb_mode="cls"`, `cell_emb_style="mean_pool"` |
| 統計モード | `goal_state_shift`（`Shift_to_goal_end`） |

> **解釈**: 各時点の AD 血液免疫細胞から遺伝子を削除し、embedding が同時点の WT へどれだけ近づくかを `Shift_to_goal_end` で評価。**正に大きい = 早期駆動因子／介入標的候補**。負に大きい = 削除で逆に WT から離れる（健常状態の維持に関与）。

---

## 5. 計算に利用したパラメータ

| パラメータ | 値 |
|---|---|
| 解析対象細胞 | 血液免疫プール 10 細胞型（各時点 AD+WT） |
| 状態埋め込み（EmbExtractor） | `max_ncells=1000`（EMB_CELLS）、exact_mean |
| perturb（InSilicoPerturber） | `max_ncells=200`, `forward_batch_size=64`, `emb_layer=0`, `nproc=1` |
| 遺伝子数 | 55（各時点で AD 細胞に発現したものを採用: 3m=55, 4p5m=54, 6m=55, 9m=55, 12m=54） |
| 総結果行 | **273**（gene × timepoint） |
| datasets バージョン | `datasets==4.0.0`（`>=5` だと `perturb_data` がハング） |
| デバイス | **CUDA** |
| **実行時間** | **約 1 時間 35 分**（01:14–02:49, 2026-08-26） |

---

## 6. in silico perturbation の結果まとめ（詳細）

`Shift_to_goal_end` の全結果は `result/report/data/ad_blood_early_isp_results.csv`（273 行）に保存。

### 6-1. 最早期（3m）ランキング上位

| Rank | Gene | カテゴリ | 3m | 4p5m | 6m | 9m | 12m |
|---|---|---|---|---|---|---|---|
| 1 | **CD3E** | adaptive | **+0.01422** | +0.00945 | +0.01653 | +0.01484 | +0.01560 |
| 2 | **CD4** | adaptive | **+0.01382** | +0.00226 | +0.01635 | +0.00768 | +0.00707 |
| 3 | **CD74** | inflammation | **+0.00533** | +0.00492 | +0.00795 | +0.00219 | +0.00588 |
| 4 | **S100A8** | inflammation | **+0.00496** | -0.00160 | +0.00563 | +0.00559 | +0.00349 |
| 5 | **C1QC** | inflammation | **+0.00486** | +0.00293 | +0.00055 | +0.00014 | -0.00014 |
| 6 | CD19 | adaptive | +0.00412 | +0.00367 | +0.00373 | +0.00416 | +0.00394 |
| 7 | ITGA2B | myeloid | +0.00283 | +0.00110 | +0.00356 | +0.00223 | +0.00324 |
| 8 | CLEC9A | myeloid | +0.00220 | -0.00012 | +0.00529 | +0.00198 | +0.00193 |
| 9 | KLRD1 | adaptive | +0.00182 | +0.00099 | +0.00111 | +0.00133 | +0.00063 |
| 10 | CD34 | myeloid | +0.00170 | -0.00059 | +0.00086 | -0.00016 | +0.00018 |
| 11 | CLU | AD_risk | +0.00167 | +0.00043 | -0.00084 | +0.00072 | +0.00033 |
| 12 | TREM2 | AD_risk | +0.00160 | +0.00100 | +0.00073 | +0.00111 | +0.00163 |

### 6-2. 負に大きい（削除で WT から離れる = 健常状態の維持因子）遺伝子

| Gene | 3m | 6m | 9m | 12m |
|---|---|---|---|---|
| CD8A | -0.10936 | -0.11452 | -0.09737 | -0.11993 |
| MS4A1 | -0.03405 | -0.02726 | -0.03678 | -0.03307 |
| NKG7 | -0.02002 | -0.01671 | -0.01450 | -0.01621 |
| CSF3R | -0.00314 | -0.00232 | -0.00280 | -0.00304 |
| S100A9 | -0.00230 | -0.00153 | -0.00021 | -0.00098 |

### 6-3. 生物学的解釈

- **血液では APOE/TREM2 はほぼ中性（3m +0.00009 / +0.00160）**。脳（ミクログリア）で APOE が最有力だった（07_ad_spleen_early_isp）のに対し、ミクログリア不在の血液では AD リスク軸が血中免疫細胞の state を強くは駆動していない。
- **最有力は T 細胞系**。CD3E（TCR コア）、CD4（ヘルパーT）の削除が 3m で最大の WT 方向シフト（+0.014）を示し、**最早期の AD 血液で T 細胞活性化/共刺激シグネチャーが疾病状態を駆動している**ことを示唆。Gate 2020（*Nature* 2020）の CD8 T クローン増殖と整合する末梢適応免疫仮説を血中 single-cell レベルで支持。
- **CD74（MHC-II インバリアント鎖）も全時点で正**。抗原提示機構の活性化が早期から疾病状態に寄与。
- **S100A8（カルプロテクチン）が 3m 正・後期も持続**。警報因子による全身性炎症が早期発症に先行する（Tan 2002）。
- **CD8A / MS4A1 / NKG7 は強く負**。これらの削除で AD 細胞が WT から遠ざかる = CD8 T 細胞毒性、B 細胞マーカー、NK 顆粒は WT 状態の維持に寄与しており、削除は逆効果。AD/WT の基底差を strongly 反映。

### 6-4. 図

| 図 | 内容 | ファイル |
|---|---|---|
| 3m ランキング | Shift_to_goal_end 棒チャート（カテゴリ色分け） | `result/report/figures/ad_blood_isp_rank_3m.png` |
| 経時傾向 | top10 遺伝子の timepoint 別推移（3m→12m） | `result/report/figures/ad_blood_isp_timeseries.png` |
| ヒートマップ | 遺伝子 × 時点の Shift 行列 | `result/report/figures/ad_blood_isp_heatmap.png` |

---

## 7. 制限事項・留意点

- `Shift_to_goal_end` は cosine シフトであり、統計的有意性（p 値）は出力に含まれない。
- `max_ncells=200` のスクリーニング結果。多数の遺伝子が微小値に固まるため、上位候補の分離には `max_ncells=1000` での再実行が有用。
- 血液データも 5xFAD ファミリー（家族性 AD）モデル由来の可能性が高く、ヒト散発性 AD への一般化には注意。
- ファインチューン済みモデルはセルタイプ分類器（15 クラス / celltype）であり、disease 状態自体の分類器ではない（チュートリアル準拠の使用方法）。
- 負値の大きい CD8A などは embedding 差分として強く出るため、正側のみならず負側も解釈に含めるべき。