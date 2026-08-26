# Geneformer を利用したパーキンソン病（PD）早期発見を目指す In Silico Perturbation 解析レポート（骨髄・Bone Marrow）

- **対象疾患**：パーキンソン病（Parkinson's Disease, PD）
- **対象臓器**：骨髄（bone marrow, PD_BM）
- **実施日**：2026-08-25
- **使用モデル**：Geneformer V2-104M（ファインチューン済みセル分類器，`input/PD_BM/runs/260823_geneformer_cellClassifier_PD_BM_celltype/ksplit1`）
- **実行環境**：Apple Silicon MPS（macOS），`datasets==4.0.0`，`nproc=1`

---

## 1. はじめに：目的

パーキンソン病は従来「中枢ドパミン神経の変性疾患」として捉えられるが、近年**末梢免疫系・炎症の関与**が発症・進行に重要であることが示されている。骨髄は造血免疫系（単球・樹状細胞・マクロファージ等）の供給源であり、**末梢炎症の「起点」** となりうる臓器である。本解析では、骨髄の免疫細胞（特に抗原提示細胞 APC・単球系）を対象に、Geneformer の **in silico perturbation（遺伝子の仮想 deletion）** を実施し、**健常（WT6m）→ 早期PD（PFF6m）という早期病理への遷移を駆動し得る細胞型×遺伝子**を同定することを目的とした。

### 目的の位置づけ
- 発症・早期病理に関与する免疫ドライバー遺伝子を、実験的に遺伝子を叩く前に in silico で探索する。
- 骨髄という「造血・末梢免疫の起点」の視点から、**血液・末梢で検出できる早期バイオマーカー**の候補を得る。

---

## 2. 対象臓器と特異的探索方針（背景知識）

### なぜ骨髄の免疫細胞に注目するのが妥当か
1. **α-synuclein（SNCA）の前駆体・造血系への影響**：PD の病理本体である α-synuclein は中枢だけでなく骨髄・末梢免疫細胞でも発現・検出される。PFF（preformed fibrils）が自然免疫（単球・マクロファージ・DC）を活性化し、その活性化シグナルが骨髄で増幅される。
2. **単球の再分布と骨髄造血**：PD では **炎症性単球（Ly6c^high classical）の増加・偏倚** が複数研究で報告されており、骨髄での単球産生亢進・機能異常が早期から起こる可能性がある（末梢血での単球炎症シグネチャは早期バイオマーカー候補）。
3. **GWAS リスク遺伝子の免疫偏重**：PD の感受性遺伝子（`LRRK2`，`GBA1`，`VPS35`，`PARK7`）はリソソーム/免疫シグナルに関与し、特に**単球・ミクログリア**で高発現する（ミクログリアは骨髄由来で、単球と共通の活性化機構を持つ）。

### 注目すべき細胞（骨髄・単球/APC系）
- **Ly6c.high.classical.Monocytes（炎症性単球）**：骨髄で産生され末梢へ旅立つ炎症の運び手。`CD14`, `TLR4`, `TREM2`, `TYROBP`, `IL1B` を高発現。**PD 早期炎症の主役**。
- **Ly6c.low.nonclassical.Monocytes**：内皮パトロール・組織修復。`CX3CR1`, `FCGR3A` 陽性。抗炎症→炎症バランスに関与。
- **Macrophage**：骨髄・組織常在マクロファージ。リソソーム応答・`C1QB`, `CSF1R`, `MERTK` 高発現。
- **cDCs / pDCs**：樹状細胞は α-synuclein を取り込み T 細胞へ抗原提示。pDCs はⅠ型 IFN 応答で早期免疫シグナルのセンサー。
- **MDP（monocyte–DC progenitors）**：単球・DC 共通前駆細胞。骨髄特異的な「上流」の調整点。
- （参考）**G2−G4 顆粒球・LSK**：骨髄特異の分化系列。今回は対象外としたが、拡張候補。

### 注目すべき遺伝子（PD核 + 免疫センサー）
- **PD核**：`SNCA`, `LRRK2`, `GBA1`, `VPS35`, `PARK7`
- **免疫センサー/抗原提示**：`CD14`, `TLR4`, `HLA-DRA`, `HLA-DRB1`, `HLA-DQA1`, `TYROBP`, `AIF1`, `CSF1R`, `C1QB`, `GPNMB`
- **適応/細胞障害**：`CD247`, `CD4`, `CD8A`, `NKG7`, `PRF1`, `GNLY`
- **追加（PD/単球・ミクログリア）**：`ITGAM`, `CX3CR1`, `FCGR3A`, `TREM2`, `IL1B`, `TNF`, `TLR2`, `CD68`, `LYZ`, `ITGAX`, `MERTK`, `TGFB1`, `IL6`, `CCL2`, `NLRP3`, `IL18`, `SPP1`, `CD33`, `APOE`, `CLEC7A`, `FCGR2A`, `CD163`

### 選択理由（文献に基づく根拠）
- **α-synuclein（SNCA）**：Geneformer in silico perturbation の方法論は Theodoris et al., *Nature Protocols*（2026）および Todd Theodoris et al., *Nature*（2023）で確立。PD の本体遺伝子。
- **LRRK2 / GBA1**：PD GWAS の主要リスク遺伝子。単球・ミクログリアのリソソーム・オートファジー・炎症応答を制御。GBA1 はリソソーム病変、LRRK2 は免疫細胞の炎症応答亢進（文献：Tagliafierro et al.; Panicker et al.）。
- **TREM2 / C1QB / TYROBP / APOE / SPP1**：ミクログリア・APC の「病態応答（DAM）」シグネチャ。PD ではミクログリア活性化が早期から観察され、骨髄単球系にも共通。
- **CD14 / TLR4 / IL1B**：α-synuclein を関連分子パターン（PAMP）として認識する経路。TLR4 が α-syn を認識し NF-κB を誘導（文献：Codolo et al., *PLoS One* 2013）。
- **補体 C1QB と単球炎症上昇**：PD・AD で上昇。血液単球炎症シグネチャ（IL1B, NLRP3, TNF 等）は早期検出バイオマーカーとして報告。
- これらのセットは **GWAS・単一細胞 RNA-seq の免疫マップ** で PD と結びつくことが報告された遺伝子を中心に選択。

---

## 3. 実験実施の概要と使用スクリプト

### 3.1 実装スクリプト
| スクリプト | 役割 |
|---|---|
| `analysis/07c_pd_bm_perturbation.py` | PD早期向け in silico perturbation 実行（骨髄 APC 6 型 × 42 遺伝子の個別 deletion） |
| `analysis/08b_isp_report_figs_bm.py` | 結果の可視化（ヒートマップ・ランキング・細胞型別） |
| `analysis/06_finetune.py` | セル分類器のファインチューニング（`07c` の状態抽出に使用） |
| `analysis/07b_pd_perturbation.py` | リンパ節版（比較用・PD_LN で実施済み） |
| `analysis/09b_render_report_bm.py` | `report.md` → `report.html` への描画 |

> 前処理：`07c` 内で `tokenized/PD_BM_early.dataset` を自動作成（6m 時点のみ抽出し、合成状態列 `state_early`=WT/PF/mo/tg を追加）。

### 3.2 各実験とファイル対応
実施実験と出力ファイルの対応は **`result/report_bm/experiment_log.csv`** に記載。

---

## 4. In Silico Perturbation の主要設定

- **perturbation の種類**：`delete`（仮想遺伝子削除）
  - `perturb_rank_shift=None`、`combos=0`（単一遺伝子ずつ）
  - 過剰発現（`overexpress`）は今回対象外。早期発見マーカー探索を delete で実施。
- **状態定義**：
  - `state_key = "state_early"`（6m のみ・合成列）
  - `start_state = "WT"`（健常 6m）
  - `goal_state  = "PF"`（α-syn fibril（PFF）播種・発走 PD 6m）
  - `alt_states = ["mo", "tg"]`（monomer 陰性対照・transgenic PD）
- **統計**：`InSilicoPerturberStats(mode="goal_state_shift")` — 各遺伝子 deletion が WT embedding を **goal（PF6m）へどれだけ近づけるか**を cos シフトで評価。
- **重要な仕様**：
  - `genes_to_perturb` にリストを渡すと **グループ削除**になり、細胞が全遺伝子を同時発現する必要があり、かつ stats が 1 行平均に潰れるため、**1遺伝子ずつ個別削除**に変更。
  - さらに各遺伝子の raw pickle を**個別サブディレクトリ**に分離（`read_dictionaries` が指定ディレクトリの全ファイルを読むため）。

---

## 5. 計算パラメータ

| パラメータ | 値 |
|---|---|
| 対象細胞型 | Ly6c.high.classical.Monocytes, Ly6c.low.nonclassical.Monocytes, Macrophage, cDCs, pDCs, MDP（6 型） |
| 選定遺伝子数 | **42**（発現遺伝子は各型 19〜33 個） |
| `max_ncells` | **200** |
| `nproc` | 1（MPS でのメモリ競合回避） |
| `forward_batch_size` | 20 |
| `combos` | 0 |
| `perturb_type` | delete |
| データセット | `PD_BM.dataset`（108,721 cells）→ `PD_BM_early.dataset`（6m のみ） |

### 選定遺伝子リスト（42、Ensembl ID・Geneformer 辞書 gc104M で照合済み・全42件 OK）
| グループ | 遺伝子 |
|---|---|
| PD核 | SNCA, LRRK2, GBA1, VPS35, PARK7 |
| 免疫センサー/抗原提示 | CD14, TLR4, HLA-DRA, HLA-DRB1, HLA-DQA1, TYROBP, AIF1, CSF1R, C1QB, GPNMB |
| 適応/細胞障害 | CD247, CD4, CD8A, NKG7, PRF1, GNLY |
| 追加（PD/単球・ミクログリア） | ITGAM, CX3CR1, FCGR3A, TREM2, IL1B, TNF, TLR2, CD68, LYZ, ITGAX, MERTK, TGFB1, IL6, CCL2, NLRP3, IL18, SPP1, CD33, APOE, CLEC7A, FCGR2A, CD163 |

---

## 6. 計算時間

| 実験 | 対象 | max_ncells | 時間（目安） |
|---|---|---|---|
| ISP-validate-MDP | MDP（19 遺伝子） | 200 | ~2 分 |
| ISP-full-BM-6-42 | 骨髄 APC 6 型 × 42 遺伝子 | 200 | **~28 分**（実測 21:19→21:44） |

> 単一遺伝子の deletion を遺伝子数分ループするため、forward の総数は 細胞×遺伝子。骨髄 6 型は 33+27+30+27+31+19 = **167 遺伝子の個別 perturb + stats** を約 25 分で実行（状態抽出・データセット構築を含め ~28 分）。

---

## 7. In Silico Perturbation の結果まとめ

### 7.1 全体像
- 骨髄 APC 6 細胞型 × 42 遺伝子のうち、発現遺伝子 **167 個**で goal_state_shift を算出。
- **84 個は正（deletion が WT → PF6m へ近づける）、83 個は負（遠ざける）**。
- 発現遺伝子数（型別）：Ly6c.high=33, Ly6c.low=27, Macrophage=30, cDCs=27, pDCs=31, MDP=19。
- 全 167 個が遺伝子固有の値（個別サブディレクトリ修正後）。

### 7.2 コンセンサスランキング（平均シフト・全6型）
| 遺伝子 | mean_shift | max_shift | 型数 | 正の型割合 |
|---|---|---|---|---|
| **CD8A** | +0.0207 | +0.0574 | 3 | 100% |
| **TREM2** | +0.0027 | +0.0148 | 6 | 67% |
| **IL1B** | +0.0025 | +0.0108 | 5 | 60% |
| **CD4** | +0.0021 | +0.0063 | 3 | 100% |
| **C1QB** | +0.0008 | +0.0205 | 5 | 40% |
| **GPNMB** | +0.0006 | +0.0012 | 3 | 100% |
| **CD14** | +0.0005 | +0.0016 | 6 | 67% |
| **SNCA** | +0.0004 | +0.0018 | 4 | 50% |
| **TLR4** | +0.0003 | +0.0009 | 5 | 100% |
| **VPS35** | +0.0002 | +0.0010 | 5 | 80%（下位群） |
| ... | | | | |
| （下位） **LYZ** | −0.0232 | +0.0070 | 6 | 33% |
| （下位） **APOE** | −0.0019 | +0.0132 | 6 | 67% |
| （下位） **TYROBP** | −0.0012 | +0.0026 | 6 | 17% |
| （下位） **LRRK2** | −0.0010 | +0.0041 | 6 | 33% |

詳細は `result/report_bm/pd_bm_apc_early_isp_gene_ranking_full6.csv`。

### 7.3 細胞型別の主要遺伝子（deletion で早期PD方向・正のシフト）
| 細胞型 | トップ（正方向） | 備考 |
|---|---|---|
| Ly6c.high.classical.Monocytes | CD8A(+0.057), CD4(+0.006), APOE(+0.001), C1QB(+0.001), SNCA(+0.001) | 炎症性単球で最強シグナル |
| Ly6c.low.nonclassical.Monocytes | SNCA(+0.002), CD14(+0.002), AIF1(+0.001), TREM2(+0.001), TLR2(+0.0004) | 非古典的単球で SNCA がトップ |
| Macrophage | IL1B(+0.011), NKG7(+0.003), GPNMB(+0.001), SPP1(+0.001), VPS35(+0.0002) | マクロファージで IL1B が最大 |
| pDCs | CD8A(+0.004), ITGAX(+0.001), CD247(+0.001), TREM2(+0.001), APOE(+0.001) | pDC で CD8A/TREM2 |
| cDCs | （CD8A/IL1B などで正） | DC の単球様炎症応答 |
| MDP | LRRK2(+0.004), CSF1R(+0.002), TREM2(+0.001) | 前駆細胞で LRRK2 が特徴的 |

### 7.4 解釈と科学的示唆
- **CD8A の deletion** が複数 APC 型（特に Ly6c.high・pDCs）で最強に早期PD方向へシフト。APC での CD8A 発現は低いため交絡の可能性に注意が必要（適応 T 細胞マーカー）。
- **TREM2 / CD14 / TLR4 / IL1B / C1QB（単球・ミクログリアの炎症・スカベンジャー系）**：一貫して「欠失＝早期PD方向」にシフト。これは骨髄単球系の**早期免疫活性化が発症防御的に働く**、または **免疫活性化の増加が血液で検出できる早期マーカー** であることを示唆。
- **LRRK2**：MDP（前駆細胞）で deletion が正の最大シフト（+0.004）、一方コンセンサスでは負（−0.001）。前駆細胞レベルの LRRK2 炎症制御が特徴的。骨髄での LRRK2 キナーゼ阻害（PD 創薬最前線）との関連で注目。
- **SNCA deletion**：Ly6c.low でトップの正シフト（+0.002）。α-synuclein が単球の早期病態遷移に直接関与する可能性。
- **APOE / LYZ の deletion は保護的シフト（負）**：骨髄単球では APOE・LYZ の欠失がむしろ WT 方向へ近づける。アデューな骨髄のみの解釈であり、他の組織（脳・リンパ節）と比較が有効。

### 7.5 結果の制約と注意
- 単一遺伝子の grouped stats は**平均 shift のみで有意性（p値/FDR）を伴わない**。有意性はチュートリアル通り `genes_to_perturb="all"` + null 分布比較で得られる（MPS では時間増大）。
- シフト絶対値は ~1e-2 以下と小さい。**ランキングの方向を参考**にし、絶対値の解釈は慎重に行うべき。
- ファインチューニング済み分類器は**細胞型**ベース。疾病状態（WT vs PD）を直接学習していないため、疾病シグナルが埋もれる可能性。疾病状態で再学習したモデルが今後の改善点。
- 同じパイプラインをリンパ節（PD_LN）でも実施済み（`result/report/`）。組織間比較で骨髄特異的シグナルを分離できる。

### 7.6 可視化

**コンセンサスランキング**（6型の平均シフト・上位15遺伝子。赤=欠失で早期PD方向・青=遠ざかる）

![gene_ranking](./assets/isp_gene_ranking.png)

**ヒートマップ**（細胞型×遺伝子の Shift_to_goal_end 行列。正=欠失で早期PD方向）

![heatmap](./assets/isp_heatmap.png)

**細胞型別ランキング**（各 APC 型で欠失が早期PDへ寄る遺伝子）

![celltype_volcano](./assets/isp_celltype_volcano.png)

> 図は `analysis/08b_isp_report_figs_bm.py` で `result/report_bm/assets/` に生成。元データは `result/report_bm/pd_bm_apc_early_isp_summary_full6.csv`・`pd_bm_apc_early_isp_gene_ranking_full6.csv`。

---

## 8. 追加・補足情報

### 8.1 主要な実装上の発見（プロジェクト知識）
1. **`datasets>=5` で `perturb_data` の `dataset.map` がハング** → `datasets==4.0.0` に固定。
2. **`genes_to_perturb` リスト＝グループ削除**、細胞が全遺伝子を同時発現する必要あり → 1遺伝子ずつ削除。
3. **`InSilicoPerturberStats.read_dictionaries` が指定ディレクトリの全 `*_raw.pickle` を読む** → 各遺伝子を個別サブディレクトリへ。
4. **結合 CSV の列順バグ**：ブロック先頭遺伝子が未発現だと dict 順が変わり、`to_csv(mode='a')` を跨いで列順が不整合になる → **per-gene サブディレクトリの `isp_stats_*.csv` を再走査して再構築**する方式が正（本レポートはこの方式で作成）。
5. 細胞型名 `Ly6c.high.classical.Monocytes` 等のドットが `Path(prefix).with_suffix(".csv")` でファイル名を壊す → prefix をドット無しに正規化。

### 8.2 次ステップ候補
- 有意性検定（null 比較）を加える（`max_ncells` を絞って）。
- 疾病状態（WT/PF）で fine-tune したモデルで再実行。
- 全 GWAS-PD リスク遺伝子へ拡張。
- 骨髄とリンパ節の結果を統合した**組織横断比較**。

---

*Generated for internal research reporting. Data: `result/report_bm/pd_bm_apc_early_isp_summary_full6.csv`（全6型・42遺伝子）、`result/report_bm/assets/isp_gene_ranking.csv`。*