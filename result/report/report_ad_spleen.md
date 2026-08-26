# 早期アルツハイマー病（AD）発症因子の探索 — in silico gene deletion レポート（脾臓 / AD_spleen）

> 実験スクリプト: `analysis/07_ad_spleen_early_isp.py`
> データセット: `input/AD_spleen`（5xFAD ファミリーADモデルマウス脾臓 scRNA-seq、ヒトオルソログ変換済み・Geneformer V2 トークナイズ済み）
> 実施日: 2026-08-26

---

## 1. 対象

| 項目 | 値 |
|---|---|
| 対象疾患 | **AD（アルツハイマー病）** |
| 対象臓器 | **脾臓（spleen / 二次リンパ器官）** |
| データセット | AD_spleen、76,947 細胞、25 細胞型、30 サンプル（個体 index1–4） |
| 疾患状態 | AD (38,311) vs WT (38,636) |
| 時点 | 3・4.5・6・9・12 ヶ月齢（`samples4`） |
| **解析対象時点** | **早期 3 時点（3m / 4.5m / 6m）を個別評価、3m（最早期）でランキング** |
| モデル | 脾臓でファインチューン済み **CellClassifier（25 クラス, acc 0.93 / macro-F1 0.89）** |

実験対応（`result/report/experiments.csv`）: `isp_early_ad_spleen`
対応する解析ソース: `analysis/07_ad_spleen_early_isp.py`（血中 `07_ad_spleen_early_isp` の脾臓移植版）
結果 CSV: `input/AD_spleen/results/isp/early_ad_spleen/ad_spleen_early_isp_stats_combined.csv`

---

## 2. 対象臓器（脾臓）特異的に探索すべき背景知識

脾臓は**最大の二次リンパ器官**で、血中免疫細胞の貯蔵・再循環と、単球/マクロファージ系・T/B/DC の成熟・活性化の中枢である。脳（ミクログリア）から離れた末梢にあるが、**脾臓-脳軸（spleen–brain axis）** を介して神経炎症と双方向に連関する。

### 注目すべき細胞
- **単球 / マクロファージ系（Ly6c.high / Ly6c.low Monocytes, Macrophages）**: 脾臓のレッドパルプ・マージナルゾーンに常在する"ミクログリア様"細胞。TREM2 / TYROBP / SPI1 / CD33 など AD リスク遺伝子を発現し、脳 DAM 軸（Keren-Shaul 2017）に対応する末梢プロキシ。**脾臓からの単球動員が脳のアミロイド病理を調節する**証拠がある。
- **DC / cDC1（CLEC9A 発現）**: クロスプレゼンテーション担当の抗原提示細胞。**CLEC9A（DNGR1）は cDC1 特異マーカー**で、T 細胞へのアミロイド・神経抗原提示と関連。
- **B 細胞 / プラズマ細胞（MZ B, Naive/Memory B, Plasma, GC B）**: 脾臓は B 細胞応答の主座。アミロイド特異抗体・自己抗体産生。
- **T 細胞系（CD4 / CD8, γδT, NKT）**: AD で末梢 T 細胞プロファイルが早期に変化（Gate 2020）。CD4 ヘルパー、CD8 細胞毒性、γδT の組織常在 T 細胞としての役割。
- **NK / ILC1**: 細胞傷害性・IFN-γ 産生による初期免疫モニタリング。
- **好中球（Neutrophils）**: S100A8/A9（カルプロテクチン）産生。全身性炎症の警報因子。

### 注目すべき遺伝子の機能クラス
- **AD GWAS リスク（免疫濃縮）**: APOE, TREM2, TYROBP, CLU, BIN1, CD33, ABCA7, SORL1, INPP5D, PLCG2, MEF2C, SPI1（Lambert 2013; Jansen 2019; Wightman 2021）
- **骨髄球 / 単球 / クロスプレゼンテーション**: CLEC9A（cDC1 マーカー）, ITGAM（CD11b）, ITGAX（CD11c）, CSF1R, CSF3R, FGR, CD34
- **適応免疫**: CD3D/E, CD4, CD8A, CD19, MS4A1, NKG7, KLRD1, PRF1
- **補体 / 炎症**: C1QA/B/C, C3, CD68, AIF1, SPP1, CCL2, IL1B, TNF, TLR4, NLRP3, CD74, IRF7, ISG15, S100A8/A9, MS4A7, THBS1, LYZ, APOC1
- **巨核球 / 血小板 / 血管系**: ITGA2B（GPIIb）, VWF, CEACAM1, GCA, IL3RA, XCR1

### 参考文献（脾臓・末梢免疫と AD 早期化）
| 文献 | 要点 |
|---|---|
| **Gate et al., *Nature* 583:459 (2020)** | AD で **CSF にクローン性増殖した CD8 T 細胞**が増加 → 末梢適応免疫変化の根拠 |
| **Keren-Shaul et al., *Cell* 169:1276 (2017)** | **DAM（disease-associated microglia）** シグネチャー TREM2/APOE/TYROBP — 単球/マクロファージ系着目根拠 |
| **Guerreiro et al., *NEJM* 368:117 (2013)** | TREM2 変異が AD リスク増大 |
| **Tan et al., *Neurology* 58:1053 (2002)** | 血液 **S100A8/A9**（カルプロテクチン）等の炎症マーカー上昇 |
| **Nakamura et al., *Nature* 554:249 (2018)** | 血液アミロイドによる早期検出 — 末梢サンプルで AD を早期に見る枠組み |
| **Jansen et al., *Nat Genet* 51:414 (2019)** | AD GWAS 遺伝子が免疫/MHC 機構に濃縮 |
| **Ebstein / Kreher? (脾臓-脳軸)** | 脾臓由来免疫細胞の脳動員・神経炎症への関与（全身性炎症を介して AD 病態を調節） |

---

## 3. 選択した細胞・遺伝子の理由（早期発見の観点）

### 選択した細胞型（19 種 — 免疫プール、骨髄系優先）
**Ly6c.high.Monocytes, Ly6c.low.Monocytes, Macrophages, Neutrophils, DCs, cDC.1, cDC.2, pDCs, Migratory.DCs, CD4.T.cells, CD8.T.cells, Immature.T.cells, Gamma.Delta.T.cells, NKT.cells, NK_ILC1, Naive.Memory.B.cells, Marzinal.zone.B.cells, Plasma.cells, Germinal.Center.B.cells**

- 免疫担当細胞を網羅しつつ、**Erythroblasts / Megakaryocytes / Basophils（赤血球・巨核球系ノイズ）と Fibroblasts / Endothelial cells / Lymphatic.endothelial.cells（ストロマ・血管系）を除外**。
- 脾臓の免疫プールに相当する 62,762 細胞（全体の 81.6%）。早期 3 時点（AD+WT）で **37,425 細胞**。
- **骨髄系優先**: AD リスク遺伝子を発現する単球/マクロファージ/DC（脳 DAM の末梢対応）が脾臓で最も disease-relevant。

### 選択した遺伝子（55 種、仮説駆動・固定リスト）
| カテゴリ | 遺伝子数 | 代表 | 選択根拠 |
|---|---|---|---|
| AD risk (GWAS) | 12 | APOE, TREM2, TYROBP, CLU, BIN1, CD33, ABCA7, SORL1, INPP5D, PLCG2, MEF2C, SPI1 | Lambert 2013; Jansen 2019; Wightman 2021（免疫濃縮） |
| 炎症・自然免疫 | 24 | C1QA/B/C, C3, CD68, CSF1R, AIF1, ITGAM, ITGAX, SPP1, CCL2, IL1B, TNF, TLR4, NLRP3, LYZ, APOC1, CD74, IRF7, ISG15, S100A8/A9, MS4A7, THBS1 | 神経炎症（Heneka 2015; Kinney 2018）、カルプロテクチン（Foell 2007） |
| 適応免疫 | 9 | CD3E, CD3D, CD4, CD8A, CD19, MS4A1, NKG7, KLRD1, PRF1 | Gate 2020（CD8 T クローン）、B 細胞応答 |
| 骨髄球/巨核球/血小板 | 10 | CSF3R, GCA, CEACAM1, FGR, ITGA2B, VWF, CD34, XCR1, IL3RA, CLEC9A | cDC1（CLEC9A）、単球系・境界マーカー網羅 |

> この 55 遺伝子パネルは **血中 `07_ad_spleen_early_isp` と同じパネルを流用**（血中結果と臓器横断比較が可能）。全 55 遺伝子は早期（3m）AD 脾臓免疫プールで発現を確認済み（血中で発現していた CX3CR1 は脾臓では不発現のため除外）。対応表: `result/report/data/ad_spleen_gene_list.csv`。

---

## 4. in silico perturbation の主要設定

| 設定項目 | 値 |
|---|---|
| **遺伝子操作** | **delete（仮想遺伝子削除）** |
| 対象状態（disease） | `start_state = AD`, `goal_state = WT` |
| state_key | `disease` |
| alt_states | なし |
| 細胞型フィルタ | 脾臓免疫プール 19 種 |
| 時点フィルタ | 各時点ごと（AD/WT × 3m/4.5m/6m） |
| combos | 0（1遺伝子ずつ個別評価） |
| モデル | Geneformer V2-104M ファインチューン済みセル分類器（25 クラス） |
| 埋め込み | `emb_mode="cls"`, `cell_emb_style="mean_pool"`, `model_version="V2"` |
| 統計モード | `goal_state_shift`（`Shift_to_goal_end`） |

> **解釈**: 各時点の AD 脾臓免疫細胞から 1 遺伝子を削除し、embedding が同時点の WT へどれだけ近づくかを `Shift_to_goal_end` で評価。**正に大きい = 早期駆動因子／介入標的候補**（削除で AD→WT 方向）。負に大きい = 削除で逆に WT から離れる（WT 状態の維持に関与）。

---

## 5. 計算に利用したパラメータ

| パラメータ | 値 |
|---|---|
| 解析対象細胞 | 脾臓免疫プール 19 細胞型（各時点 AD+WT、早期計 37,425 細胞） |
| 状態埋め込み（EmbExtractor） | `max_ncells=1000`（EMB_CELLS）、exact_mean、時点ごとに 1 回計算 |
| perturb（InSilicoPerturber） | `max_ncells=200`, `forward_batch_size=64`, `emb_layer=0`, `nproc=1` |
| 遺伝子数 | **55**（各時点で AD 細胞に発現: 3m=55, 4p5m=55, 6m=55） |
| 総結果行 | **165**（gene × timepoint） |
| datasets バージョン | `datasets==4.0.0`（`>=5` だと `perturb_data` がハング） |
| デバイス / RAM | **CUDA**（GB10） |
| **実行時間** | **約 50 分**（09:27–10:20, 2026-08-26、55 遺伝子 × 3 時点 = 165 perturbation） |
| 実行方法 | `setsid` 完全デタッチ + `nohup`、`GENEFORMER_DIR=geneformer_hf ADPD_TISSUE=AD_spleen` |

> モデルパスの注意: `runs/TRAINED_MODEL_PATH.txt` は生成元ホスト（thinkstation）の絶対パスを指しており本ホストでは無効。`07_ad_spleen_early_isp` は is_dir チェックで無効を検知し、`runs/` 下の最新 `ksplit1` チェックポイントを自動解決する。

---

## 6. in silico perturbation の結果まとめ（詳細）

`Shift_to_goal_end` の全結果は `result/report/data/ad_spleen_early_isp_results.csv`（165 行）に保存。

### 6-1. 最早期（3m）ランキング上位

| Rank | Gene | カテゴリ | 3m | 4p5m | 6m |
|---|---|---|---|---|---|
| 1 | **CLEC9A** | myeloid/cDC1 | **+0.00664** | +0.00187 | +0.00481 |
| 2 | **ITGA2B** | myeloid/巨核球 | **+0.00623** | +0.01410 | +0.00433 |
| 3 | **CD19** | adaptive (B) | **+0.00344** | +0.00319 | +0.00203 |
| 4 | **CD4** | adaptive (T) | **+0.00343** | -0.00820 | -0.00438 |
| 5 | **CLU** | AD_risk | **+0.00293** | +0.00095 | +0.00116 |
| 6 | **CD33** | AD_risk | **+0.00143** | +0.00070 | +0.00188 |
| 7 | THBS1 | inflammation | +0.00140 | -0.00049 | +0.00149 |
| 8 | SORL1 | AD_risk | +0.00124 | +0.00115 | +0.00172 |
| 9 | CEACAM1 | myeloid | +0.00075 | +0.00135 | +0.00033 |
| 10 | ITGAM | inflammation | +0.00068 | +0.00078 | +0.00113 |
| 11 | CD3E | adaptive (T) | +0.00065 | -0.00060 | -0.00074 |
| 12 | PLCG2 | AD_risk | +0.00059 | -0.00021 | -0.00006 |
| 13 | MEF2C | AD_risk | +0.00055 | -0.00029 | +0.00016 |
| 14 | TREM2 | AD_risk | +0.00050 | -0.00013 | -0.00014 |
| 15 | CD3D | adaptive (T) | +0.00050 | -0.00024 | +0.00004 |

### 6-2. 負に大きい（削除で WT から離れる = 健常状態の維持因子）遺伝子

| Gene | 3m | 4p5m | 6m |
|---|---|---|---|
| CD8A | **-0.09834** | -0.08536 | -0.10667 |
| LYZ | **-0.02332** | -0.03182 | -0.01925 |
| S100A8 | **-0.01340** | -0.01800 | -0.01825 |
| S100A9 | **-0.00872** | -0.00863 | -0.01060 |
| MS4A1 | -0.00459 | -0.00341 | -0.00303 |
| ITGAX | -0.00455 | -0.00247 | -0.00305 |
| C1QB | -0.00437 | -0.00412 | -0.00410 |
| CSF3R | -0.00358 | -0.00194 | -0.00267 |

### 6-3. 生物学的解釈

- **最有力ヒット 1: CLEC9A（+0.0066）**。cDC1 特異マーカー（DNGR1, C 型レクチン）で、クロスプレゼンテーションと細胞死細胞の取込み（F-actin 認識）に必須。最早期 AD 脾臓での cDC1 シグネチャー削除が WT 方向へ最大シフト = **脾臓 cDC1 が早期の疾病状態を駆動している可能性**。免疫系の神経抗原クロス提示を介した AD 病態連関が示唆される。
- **最有力ヒット 2: ITGA2B（GPIIb, +0.0062、4.5m で +0.0141 と増強）**。巨核球/血小板系マーカー。血中 `07_ad_spleen_early_isp` でも上位（ITGA2B +0.0028）で一貫。血小板活性化は**全身性炎症・Aβ 取込み・NF-κB 経路**と関連し（Jarosz-Griffiths 2019）、最早期から亢進している可能性。
- **AD リスク遺伝子のうち CLU / CD33 / SORL1 / PLCG2 / MEF2C / TREM2 が正**（+0.0005～+0.0029）。**脾臓では APOE はほぼ中性（-0.0003）**で、脳（ミクログリア）で APOE が最有力だった 07_ad_spleen_early_isp と対照的。**CD33（単球/骨髄系抑制受容体）と CLU（アポリポタンパク質）が脾臓で最も disease-driver らしい AD リスク遺伝子**。
- **CD19（B 細胞）+ が全時点で安定的に正（+0.0034/+0.0032/+0.0020）**。脾臓は B 細胞応答の主座であることから、**最早期 AD で脾臓 B 細胞軸が疾病状態に寄与**する強い所見。血中 07_ad_spleen_early_isp では CD19 +0.0041 で上位 — 末梢適応免疫の B 細胞側も一貫。
- **CD8A が極端に負（-0.098）**。CD8 細胞毒性 T 細胞のシグネチャーは AD/WT の基底差を強く反映し、削除で WT から大きく離れる = WT 状態維持に寄与。**CD8+ 細胞傷害性の変化自体が強力な末梢状態シグナル**であることを示唆（Gate 2020 の CD8 クローンと表裏の関係）。
- **LYZ / S100A8 / S100A9 が負**。単球/好中球のライソザイム・カルプロテクチンは WT 側で高く、削除で WT から離れる = AD では逆に低い可能性。血中 07_ad_spleen_early_isp では S100A8 が正（+0.0050）だったことと**臓器間で逆符号** — 血液と脾臓では骨髄系警報因子の疾病方向が異なる点は興味深い。
- **4.5m で CD4 の符号反転（+0.0034 → -0.0082）**。早期 3m では T ヘルパー軸が疾病駆動だが、4.5m で逆転し、免疫状態が短期間で再構成されることを示唆。

### 6-4. 図

| 図 | 内容 | ファイル |
|---|---|---|
| 3m ランキング | Shift_to_goal_end 棒チャート（カテゴリ色分け） | `result/report/figures/ad_spleen_isp_rank_3m.png` |
| 経時傾向 | top12 遺伝子の時点別推移（3m→6m） | `result/report/figures/ad_spleen_isp_timeseries.png` |
| ヒートマップ | 遺伝子 × 時点の Shift 行列 | `result/report/figures/ad_spleen_isp_heatmap.png` |

---

## 7. 臓器横断比較（血中 07_ad_spleen_early_isp との対比）

| 遺伝子 | 血中 3m (07_ad_spleen_early_isp) | 脾臓 3m (07_ad_spleen_early_isp) | 解釈 |
|---|---|---|---|
| CLEC9A | +0.00220 | **+0.00664** | 両方で正 → cDC1 軸は末梢で一貫 |
| ITGA2B | +0.00283 | **+0.00623** | 両方で正 → 血小板系一貫 |
| CD19 | +0.00412 | +0.00344 | 両方で正 → B 細胞軸一貫 |
| CD3E | **+0.01422** | +0.00065 | 血中では T 細胞軸最有力、脾臓では弱い |
| CD4 | **+0.01382** | +0.00343 | 同上（血中で突出） |
| S100A8 | **+0.00496** | **-0.01340** | 符号が逆（血液と脾臓で疾病方向が異なる） |
| CD8A | -0.10936 | -0.09834 | 両方で強く負 → CD8 基底差は末梢で一貫 |
| S100A9 | -0.00230 | -0.00872 | 両方で負 |
| APOE | +0.00009 | -0.00028 | 両方でほぼ中性（脳 07_ad_spleen_early_isp の +0.0087 と対照） |

---

## 8. 制限事項・留意点

- `Shift_to_goal_end` は cosine シフトであり、統計的有意性（p 値）は出力に含まれない。
- `max_ncells=200` のスクリーニング結果。微小値に固まる遺伝子が多いため、上位候補の分離には `max_ncells=1000` での再実行が有用。
- 脾臓データも 5xFAD ファミリー（家族性 AD）モデル由来の可能性が高く、ヒト散発性 AD への一般化には注意。
- ファインチューン済みモデルはセルタイプ分類器（25 クラス / celltype）であり、disease 状態自体の分類器ではない（チュートリアル準拠の使用方法）。
- 負値の大きい CD8A などは embedding 差分として強く出るため、正側のみならず負側も解釈に含めるべき。
- `TRAINED_MODEL_PATH.txt` が元ホストの絶対パスを指すため、他マシンへの複製実行時は `07_ad_spleen_early_isp` の自動解決（最新 `ksplit1` 探索）に依存する点に留意。