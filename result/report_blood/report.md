# Geneformer を利用したパーキンソン病（PD）早期発見を目指す In Silico Perturbation 解析レポート（血液・Blood）

- **対象疾患**：パーキンソン病（Parkinson's Disease, PD）
- **対象臓器**：血液（peripheral blood, PD_blood）
- **実施日**：2026-08-25（検証 23:13 〜、本番 23:57 〜）
- **使用モデル**：Geneformer V2-104M（ファインチューン済みセル分類器，`input/PD_blood/runs/260823_geneformer_cellClassifier_PD_blood_celltype/ksplit1`）
- **実行環境**：Apple Silicon MPS（macOS）、`datasets==4.0.0`、`nproc=1`

---

## 1. はじめに：目的

パーキンソン病は従来「中枢ドパミン神経の変性疾患」として捉えられるが、近年**末梢免疫応答・炎症の関与**が発症・進行に重要であることが示されている。特に**血液**は臨床で最も採取が容易な生体材料であり、**血液中の免疫細胞の変化は早期バイオマーカー**として大きな期待がある。

本解析では、血液の免疫細胞（単球・好中球・樹状細胞・T 細胞・NK 細胞など）を対象に、Geneformer の **in silico perturbation（遺伝子の仮想 deletion）** を実施し、**健常（WT6m）→ 早期 PD（PFF6m）という早期病理への遷移を駆動し得る細胞型 × 遺伝子**を同定することを目的とした。

### 目的の位置づけ
- 発症・早期病理に関与する免疫ドライバー遺伝子を、実験的に遺伝子を操作する前に in silico で探索する。
- 血液という「非侵襲的に採取できる材料」の視点から、**血液・末梢で検出できる早期バイオマーカー**の候補を得る。

---

## 2. 対象臓器と特異的探索方針（背景知識）

### なぜ血液の免疫細胞に注目するのが妥当か
1. **血液は早期病変の「窓口」**：PD の病理である α-synuclein（SNCA）は中枢だけでなく血液の免疫細胞でも発現・検出される。PFF（preformed fibrils）が自然免疫（単球・好中球・樹状細胞）を活性化し、その活性化シグナルが血液で観測できる。
2. **単球の再分布と炎症偏倚**：PD では**炎症性単球（Ly6c^high classical）の増加・機能偏倚**が複数研究で報告されており、血液単球の炎症シグネチャは早期バイオマーカー候補。
3. **GWAS リスク遺伝子の免疫偏重**：PD の感受性遺伝子（`LRRK2`、`GBA1`、`VPS35`、`PARK7`）はリソソーム/免疫シグナルに関与し、特に**単球・ミクログリア**で高発現する。血液の骨髄系免疫細胞はミクログリアと共通の活性化機構を持つ。

### 注目すべき細胞（血液・単球/APC 系 + 適応免疫系）
- **Ly6c.high.classical.Monocytes（炎症性単球）**：`CD14`、`TLR4`、`TREM2`、`TYROBP`、`IL1B` を高発現。**PD 早期炎症の主役**。
- **Ly6c.low.nonclassical.Monocytes（非古典的単球）**：内皮パトロール・組織修復。炎症-抗炎症バランスに関与。
- **Neutrophils（好中球）**：`LYZ`、`TLR2`、`NLRP3`、`IL1B` 高発現。PD での好中球活性化・NET の関与が報告。
- **cDCs / pDCs（樹状細胞）**：α-synuclein を取り込み T 細胞へ抗原提示。pDC はⅠ型 IFN 応答で早期免疫のセンサー。
- **CD8.T / CD4.T / NK_ILC1**：適応免疫・細胞毒性。**α-syn 特異的な T 細胞応答は PD 早期に検出され得る**バイオマーカー候補。

### 注目すべき遺伝子（PD コア + 免疫センサー）
- **PD コア**：`SNCA`、`LRRK2`、`VPS35`、`PARK7`（※ `GBA1` は血液では全細胞型未発現 → 除外）
- **免疫センサー / 抗原提示**：`CD14`、`TLR4`、`TLR2`、`TYROBP`、`AIF1`、`CSF1R`、`C1QB`、`GPNMB`（※ `HLA-DRA/HLA-DRB1` は血液未発現 → 除外）
- **適応 / 細胞毒性**：`CD247`、`CD4`、`CD8A`、`NKG7`、`PRF1`、`GNLY`
- **追加（PD / 単球・ミクログリア）**：`ITGAM`、`TREM2`、`IL1B`、`TNF`、`CD68`、`LYZ`、`ITGAX`、`MERTK`、`TGFB1`、`IL6`、`CCL2`、`NLRP3`、`IL18`、`SPP1`、`CD33`、`APOE`、`CLEC7A`、`FCGR2A`、`CD163`（※ `CX3CR1`、`FCGR3A` は血液未発現 → 除外）

### 選択理由（文献に基づく根拠）
- **α-synuclein（SNCA）**：in silico perturbation の方法論は Theodoris et al., *Nature Protocols*（2026）および *Nature*（2023）で確立。PD の本体遺伝子。
- **LRRK2 / GBA1 / VPS35 / PARK7**：PD の主要リスク遺伝子。単球・ミクログリアのリソソーム・オートファジー・炎症応答を制御（Tagliafierro ほか；Panicker ほか）。
- **TREM2 / C1QB / TYROBP / APOE / SPP1**：ミクログリア・APC の「病態応答（DAM）」シグネチャ。PD ではミクログリア活性化が早期から観察され、血液の単球系にも共通。
- **CD14 / TLR4 / IL1B**：α-synuclein を関連分子パターン（PAMP）として認識する経路。TLR4 が α-syn を認識し NF-κB を誘導（Codolo et al., *PLoS One* 2013）。
- **補体 C1QB と単球炎症上昇**：PD で上昇。血液単球炎症シグネチャ（IL1B、NLRP3、TNF 等）は早期検出バイオマーカーとして報告。
- これらのセットは **GWAS・単一細胞 RNA-seq の免疫マップ**で PD と結びつくことが報告された遺伝子を中心に選択。

---

## 3. 実験実施の概要と使用スクリプト

### 3.1 実装スクリプト
| スクリプト | 役割 |
|---|---|
| `analysis/07d_pd_blood_perturbation.py` | PD 早期向け in silico perturbation 実行（血液 8 型 × 37 遺伝子の個別 deletion） |
| `analysis/08c_isp_report_figs_blood.py` | 結果の可視化＋正しい結果テーブルの再構築（ヒートマップ・ランキング・細胞型別） |
| `analysis/06_finetune.py` | セル分類器のファインチューニング（`07d` の状態抽出に使用） |
| `analysis/07b_pd_perturbation.py`・`07c_pd_bm_perturbation.py` | リンパ節・骨髄版（比較用・他組織で実施済み） |
| `analysis/09c_render_report_blood.py` | `report.md` → `report.html` への描画 |

> 前処理：`07d` 内で `tokenized/PD_blood_early.dataset` を自動作成（6m 時点のみ抽出し、合成状態列 `state_early`=WT/PF/mo/tg を追加）。

### 3.2 各実験とファイル対応
実施実験と出力ファイルの対応は **`result/report_blood/experiment_log.csv`** に記載。

---

## 4. In Silico Perturbation の主要設定

- **perturbation の種類**：`delete`（仮想遺伝子削除）
  - `perturb_rank_shift=None`、`combos=0`（単一遺伝子ずつ）
  - 過剰発現（`overexpress`）は今回対象外。早期発見マーカー探索を delete で実施。
- **状態定義**：
  - `state_key = "state_early"`（6m のみ・合成列）
  - `start_state = "WT"`（健常 6m）
  - `goal_state  = "PF"`（α-syn fibril（PFF）播種・早期 PD 6m）
  - `alt_states = ["mo", "tg"]`（monomer 陰性対照・transgenic PD）
- **統計**：`InSilicoPerturberStats(mode="goal_state_shift")` — 各遺伝子 deletion が WT embedding を **goal（PF6m）へどれだけ近づけるか**を cos シフトで評価。
- **重要な仕様**：
  - `genes_to_perturb` にリストを渡すと **グループ削除**になり、細胞が全遺伝子を同時発現する必要があり、かつ stats が 1 行平均に潰れる → **1 遺伝子ずつ個別削除**に変更。
  - 各遺伝子の raw pickle を**個別サブディレクトリ**に分離（`read_dictionaries` が指定ディレクトリの全ファイルを読むため）。

---

## 5. 計算パラメータ

| パラメータ | 値 |
|---|---|
| 対象細胞型 | Ly6c.high.classical.Monocytes, Ly6c.low.nonclassical.Monocytes, Neutrophils, cDCs, pDCs, CD8.T.cells, CD4.T.cells, NK_ILC1（8 型） |
| 選定遺伝子数 | **37**（発現遺伝子は各型 13〜32 個、計 **207** 個の細胞×遺伝子を個別削除） |
| `max_ncells` | **200**（検証 run は 50） |
| `nproc` | 1（MPS でのメモリ競合回避） |
| `forward_batch_size` | 20 |
| `combos` | 0 |
| `perturb_type` | delete |
| データセット | `PD_blood.dataset`（65,013 cells）→ `PD_blood_early.dataset`（6m のみ 20,741 cells） |

### 選定遺伝子リスト（37、Ensembl ID・Geneformer 辞書 gc104M で照合済み・全37件 OK）
| グループ | 遺伝子 |
|---|---|
| PD コア | SNCA, LRRK2, VPS35, PARK7 |
| 免疫センサー / 抗原提示 | CD14, TLR4, TLR2, TYROBP, AIF1, CSF1R, C1QB, GPNMB |
| 適応 / 細胞毒性 | CD247, CD4, CD8A, NKG7, PRF1, GNLY |
| 追加（PD / 単核・ミクログリア） | ITGAM, TREM2, IL1B, TNF, CD68, LYZ, ITGAX, MERTK, TGFB1, IL6, CCL2, NLRP3, IL18, SPP1, CD33, APOE, CLEC7A, FCGR2A, CD163 |

> ※ 元の 42 遺伝子リストから、血液で未発現の GBA1, HLA-DRA, HLA-DRB1, CX3CR1, FCGR3A の 5 つを除外し **37** に削減。

---

## 6. 計算時間

| 実験 | 対象 | max_ncells | 時間（目安） |
|---|---|---|---|
| ISP-blood-trial-1 | cDCs（検証、37 遺伝子） | 50 | ~1 分 |
| ISP-blood-full-8-37 | 血液 8 型 × 37 遺伝子 | 200 | **~28 分**（実測 23:57 起動） |

> 単一遺伝子の deletion を遺伝子数分ループするため、forward の総数は 細胞×遺伝子。血液 8 型は 27+31+25+27+13+25+27+27 = **207 遺伝子の個別 perturb + stats** を約 28 分で実行（状態抽出・データセット構築を含む）。

---

## 7. In Silico Perturbation の結果まとめ

### 7.1 全体像
- 血液 8 細胞型 × 37 遺伝子のうち、発現遺伝子 **207 個**で goal_state_shift を算出。
- **80 個は正（deletion が WT → PF6m へ近づける）、87 個は負（遠ざける）**。
- 発現遺伝子数（型別）：Ly6c.high=32, Ly6c.low=31, Neutrophils=25, cDCs=27, pDCs=13, CD8.T=25, CD4.T=27, NK_ILC1=27。
- 全 207 個が遺伝子固有の値（個別サブディレクトリ修正後）。

### 7.2 コンセンサスランキング（平均シフト・全8型）
| 遺伝子 | mean_shift | max_shift | 型数 | 正の型割合 |
|---|---|---|---|---|
| **CD8A** | +0.00519 | +0.13453 | 6 | 67% |
| **AIF1** | +0.00376 | +0.01079 | 6 | 67% |
| **IL6** | +0.00268 | +0.00268 | 1 | 100% |
| **ITGAM** | +0.00148 | +0.00727 | 7 | 57% |
| **SNCA** | +0.00138 | +0.00249 | 6 | 83% |
| **CSF1R** | +0.00119 | +0.00999 | 7 | 57% |
| **PRF1** | +0.00083 | +0.00308 | 5 | 60% |
| **FCGR2A** | +0.00067 | +0.00859 | 5 | 29% |
| **LRRK2** | +0.00064 | +0.00504 | 7 | 14% |
| **TLR4** | +0.00055 | +0.00353 | 7 | 57% |
| **NKG7** | +0.00051 | +0.00755 | 7 | 43% |
| **CD68** | +0.00022 | +0.00221 | 8 | 75% |
| **TREM2** | +0.00010 | +0.00082 | 2 | 50% |
| **IL1B** | +0.00007 | +0.00056 | 8 | 63% |
| （下位） **TYROBP** | −0.00887 | +0.02084 | 8 | 25% |
| （下位） **CD33** | −0.00293 | +0.00422 | 6 | 50% |
| （下位） **CD4** | −0.00286 | +0.05512 | 7 | 43% |
| （下位） **C1QB** | −0.01509 | +0.00026 | 2 | 50% |

詳細は `result/report_blood/pd_blood_apc_early_isp_gene_ranking_full6.csv`。

### 7.3 細胞型別の主要遺伝子（deletion で早期 PD 方向・正のシフト）
| 細胞型 | 正（早期 PD 方向） | 負（保護的） |
|---|---|---|
| CD4.T.cells | CD8A(+0.135), NKG7(+0.008), AIF1(+0.008), LRRK2(+0.005) | CD4(−0.043) |
| CD8.T.cells | CD4(+0.055), TYROBP(+0.004), ITGAM(+0.003) | CD8A(−0.105) |
| Ly6c.high.Monocytes | NKG7(+0.005), SNCA(+0.001) | LYZ(−0.016), TYROBP(−0.009) |
| Ly6c.low.Monocytes | AIF1(+0.011), SNCA(+0.002), CD8A(+0.002) | C1QB(−0.030), TYROBP(−0.006) |
| NK_ILC1 | CD4(+0.003), IL6(+0.003) | TYROBP(−0.020) |
| Neutrophils | CD8A, AIF1, CD33（微正） | TYROBP(−0.001) |
| cDCs | LYZ(+0.023), TYROBP(+0.021), CSF1R(+0.010) | CD4(−0.035), CD33(−0.021) |
| pDCs | CD33(+0.004), ITGAX(+0.001), APOE(+0.001) | TYROBP(−0.059) |

### 7.4 解釈と科学的示唆
- **CD8A の deletion が CD4.T で最大の正シフト（+0.135）**：CD8 T のシグナル欠失が CD4 T を早期 PD 側へ強く動かす。**T 細胞免疫応答が早期 PD 病理の駆動因子**である可能性を示唆。
- **SNCA（α-synuclein）deletion が 6 型中 5 型で正（83%）**：α-syn の欠失が複数の血液免疫細胞を早期 PD 側へ寄せる＝**α-syn 発現自体が免疫恒常の維持に寄与**している可能性。
- **TYROBP が cell 型依存的に反転**：cDCs では正（+0.021）、Ly6c.low・NK・pDCs では負（−0.006 〜 −0.059）。**同一遺伝子でも細胞型によって早期 PD への寄与が反転**することを示す。
- **LRRK2**：CD4.T で正シフト（+0.005）。LRRK2 キナーゼ阻害（PD 免疫創薬の最前線）と免疫系の関連で注目。
- **TGFB1 / VPS35 / C1QB は主に負（保護的）**：これらの欠失は早期 PD から遠ざける＝炎症を支える遺伝子群。

### 7.5 結果の制約と注意
- 単一遺伝子の stats は**平均 shift のみで有意性（p 値/FDR）を伴わない**。有意性はチュートリアル通り `genes_to_perturb="all"` + null 分布比較で得られる（MPS では時間増大）。
- シフト絶対値は ~1e-2 以下と小さい。**ランキングの方向を参考**にし、絶対値の解釈は慎重に行うべき。
- ファインチューニング済み分類器は**細胞型**ベース。疾病状態（WT/PD）を直接学習していないため疾病シグナルが埋もれる可能性。
- 同じパイプラインをリンパ節（PD_LN）・骨髄（PD_BM）でも実施済み。**組織間比較で血液特異的シグナルを分離**できる。

### 7.6 可視化

**コンセンサスランキング**（8 型の平均シフト・上位 15 遺伝子。赤=欠失で早期 PD 方向・青=遠ざかる）

![gene_ranking](./assets/isp_gene_ranking.png)

**ヒートマップ**（細胞型×遺伝子の Shift_to_goal_end 行列。正=欠失で早期 PD 方向）

![heatmap](./assets/isp_heatmap.png)

**細胞型別ランキング**（各細胞型で欠失が早期 PD へ寄る遺伝子）

![celltype_volcano](./assets/isp_celltype_volcano.png)

> 図は `analysis/08c_isp_report_figs_blood.py` で `result/report_blood/assets/` に生成。元データは `result/report_blood/pd_blood_apc_early_isp_summary_full6.csv`・`pd_blood_apc_early_isp_gene_ranking_full6.csv`。

---

## 8. 追加・補足情報

### 8.1 主要な実装上の発見（プロジェクト知識）
1. **`datasets>=5` で `perturb_data` の `dataset.map` がハング** → `datasets==4.0.0` に固定。
2. **`genes_to_perturb` リスト＝グループ削除**、細胞が全遺伝子を同時に発現する必要 → 1 遺伝子ずつ削除。
3. **`InSilicoPerturberStats.read_dictionaries` が指定ディレクトリの全 `*_raw.pickle` を読む** → 各遺伝子を個別サブディレクトリへ。
4. **結合 CSV の列順バグ**：`to_csv(mode='a')` を跨ぐと列順が不整合になり行が欠落する → **per-gene サブディレクトリの `isp_stats_*.csv` を再走査して再構築**する方式が正（本レポートはこの方式で 207 行すべて取得）。
5. 細胞型名のドット（`Ly6c.high.classical.Monocytes` 等）が `Path(prefix).with_suffix(".csv")` でファイル名を壊す → prefix をドット無しに正規化。

### 8.2 次ステップ候補
- 有意性検定（null 比較）を加える。
- 疾病状態（WT/PF）で fine-tune したモデルで再実行。
- 全 GWAS-PD リスク遺伝子へ拡張。
- リンパ節・骨髄・血液の結果を統合した**組織横断比較**。

---

*Generated for internal research reporting. Data: `result/report_blood/pd_blood_apc_early_isp_summary_full6.csv`（全 8 型・207 件）、`result/report_blood/assets/isp_gene_ranking.csv`。*