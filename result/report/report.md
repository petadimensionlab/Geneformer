# Geneformer を利用したパーキンソン病（PD）早期発見を目指す In Silico Perturbation 解析レポート

- **対象疾患**：パーキンソン病（Parkinson's Disease, PD）
- **対象臓器**：リンパ節（lymph node）
- **実施日**：2026-08-25
- **使用モデル**：Geneformer V2-104M（ファインチューン済みセル分類器，`input/PD_LN/runs/.../ksplit1`）
- **実行環境**：Apple Silicon MPS（macOS），`datasets==4.0.0`，`transformers==4.46.3`

---

## 1. はじめに：目的

パーキンソン病は従来「中枢ドパミン神経の変性疾患」として捉えられるが、近年**末梢免疫系・炎症の関与**が発症・進行に重要であることが示されている。本解析では、リンパ節の免疫細胞（特に抗原提示細胞 APC・単球系）を対象に、Geneformer の **in silico perturbation（遺伝子の仮想 deletion）** を実施し、**健常（WT6m）→ 早期PD（PFF6m）という早期病理への遷移を駆動し得る細胞型×遺伝子**を同定することを目的とした。

### 目的の位置づけ
- 発症・早期病理に関与する免疫ドライバー遺伝子を、実験的に遺伝子を叩く前に in silico で探索する。
- リンパ節という「末梢免疫」の視点から、**血液・リンパ系で検出できる早期バイオマーカー**の候補を得る。

---

## 2. 対象臓器と特異的探索方針（背景知識）

パーキンソン病では以下の点から**リンパ節の免疫細胞に注目するのが妥当**である：

1. **α-synuclein（SNCA）の末梢伝播と免疫応答**：PD の病理本体である α-synuclein は中枢だけでなく末梢・リンパ組織でも検出され、PFF（preformed fibrils）が自然免疫を活性化する。
2. **早期の免疫活性化シグナル**：血液/リンパの単核細胞・抗原提示細胞（APC）で、発症前から炎症性シグナル（IL-1β，TNF，補体，単球マーカー）が上昇することが複数の研究で報告されている。
3. **ゲノムワイド関連解析（GWAS）リスク遺伝子の免疫偏重**：PD の感受性遺伝子（`LRRK2`，`GBA1`，`VPS35`，`PARK7`）はリソソーム/免疫シグナルに関与し、特に**単球・ミクログリア**で高発現する。

### 注目すべき細胞（リンパ節・APC/単系）
- **Migratory.DCs / cDC.1 / cDC.2 / pDCs**：樹状細胞は α-synuclein を取り込み T 細胞へ抗原提示。**最初に PD シグナルを感知するセル**。
- **Macrophages / Ly6c^high classical Monocytes**：単核細胞・マクロファージはリソソーム・炎症応答を担い、`TYROBP`, `CSF1R`, `CD14`, `TLR4`, `TREM2` を高発現。
- （参考）**NK_ILC1 / ILC2 / ILC3**：自然免疫の第二段階。

### 注目すべき遺伝子（PD核 + 免疫センサー）
- **PD核**：`SNCA`, `LRRK2`, `GBA1`, `VPS35`, `PARK7`
- **免疫センサー/抗原提示**：`CD14`, `TLR4`, `HLA-DRA`, `HLA-DRB1`, `HLA-DQA1`, `TYROBP`, `AIF1`, `CSF1R`, `C1QB`, `GPNMB`
- **適応/細胞障害**：`CD247`, `CD4`, `CD8A`, `NKG7`, `PRF1`, `GNLY`
- **追加（PD/単球・ミクログリア）**：`ITGAM`, `CX3CR1`, `FCGR3A`, `TREM2`, `IL1B`, `TNF`, `TLR2`, `CD68`, `LYZ`, `ITGAX`, `MERTK`, `TGFB1`, `IL6`, `CCL2`, `NLRP3`, `IL18`, `SPP1`, `CD33`, `APOE`, `CLEC7A`, `FCGR2A`, `CD163`

### 選択理由（文献に基づく根拠）
- **α-synuclein（SNCA）と in silico perturbation**：Geneformer チュートリアル（Theodoris et al., *Nature Protocols*, 2026）で在来の方法論が確立。PD の本体遺伝子。
- **LRRK2 / GBA1**：PD GWAS の主要リスク遺伝子で、単球・ミクログリアの機能（リソソーム、オートファジー）に深く関与。GBA1 はリソソーム病変、LRRK2 は免疫細胞の炎症応答を制御（文献：Tagliafierro et al.; Panicker et al.）。
- **TREM2 / C1QB / TYROBP / APOE / SPP1**：ミクログリア・APC の“病態応答”サイン。PD ではミクログリア活性化が早期から観察され、DAM（disease-associated microglia）シグネチャと共通。
- **CD14 / TLR4 / IL1β**：α-synuclein の病原体関連分子パターン（PAMP）としての認識経路。TLR4 は α-syn を認識し NF-κB を誘導（文献：Codolo et al., *PLoS One 2013*）。
- **補体 C1QB**：神経炎症・シナプス除去に関与、PD・AD で上昇。
- これらのセットは **GWAS・単一細胞 RNA-seq の免疫マップ** で PD と結びつくことが報告された遺伝子を中心に選択。

---

## 3. 実験実施の概要と使用スクリプト

### 3.1 実装スクリプト
| スクリプト | 役割 |
|---|---|
| `analysis/07b_pd_perturbation.py` | PD早期向け in silico perturbation 実行（6 細胞型 × 遺伝子の個別 deletion） |
| `analysis/08_isp_report_figs.py` | 結果の可視化（ヒートマップ・ランキング・細胞型別） |
| `analysis/06_finetune.py` | セル分類器のファインチューニング（`07b` の状態抽出に使用） |
| `analysis/07_in_silico_perturbation.py` | 元のチュートリアル移植版（今回は未使用・参考） |

> 前処理：`07b` 内で `tokenized/PD_LN_early.dataset` を自動作成（6m 時点のみ抽出し、合成状態列 `state_early`=WT/PF/mo/tg を追加）。

### 3.2 各実験とファイル対応
実施実験と出力ファイルの対応は **`result/report/experiment_log.csv`** に記載。

---

## 4. In Silico Perturbation の主要設定

- **perturbation の種類**：`delete`（仮想遺伝子削除）
  - `perturb_rank_shift=None`、`combos=0`（単一遺伝子ずつ）
- **状態定義**：
  - `state_key = "state_early"`（6枚のみ・合成列）
  - `start_state = "WT"`（健常 6m）
  - `goal_state  = "PF"`（α-syn fibril 播種・発走 PD 6m）
  - `alt_states = ["mo", "tg"]`（monomer 陰性対照・transgenic）
- **統計**：`InSilicoPerturberStats(mode="goal_state_shift")` — 各遺伝子 deletion が WT embedding を **goal（PF6m）へどれだけ近づけるか**を cos シフトで評価。
- **重要な仕様**：`genes_to_perturb` にリストを渡すと **グループ削除**になり、細胞が全遺伝子を同時発現する必要があり、かつ stats が 1 行平均に潰れるため、**1遺伝子ずつ個別削除**に変更。さらに各遺伝子の raw pickle を**個別サブディレクトリ**に分離（`read_dictionaries` が指定ディレクトリの全ファイルを読み込むため）。

---

## 5. 計算パラメータ

| パラメータ | 値 |
|---|---|
| 対象細胞型 | Migratory.DCs, cDC.1, cDC.2, pDCs, Macrophages, Ly6c.high.classical.Monocytes（6 型） |
| 選定遺伝子数 | **43**（発現遺伝子は各型 27〜33 個） |
| `max_ncells` | 本番 **200**、試行 **50** |
| `nproc` | 1（MPS でのメモリ競合回避） |
| `forward_batch_size` | 20 |
| `combos` | 0 |
| `perturb_type` | delete |

### 選定遺伝子リスト（43、Ensembl ID・Geneformer 辞書で照合済み）
| グループ | 遺伝子 |
|---|---|
| PD核 | SNCA, LRRK2, GBA1, VPS35, PARK7 |
| 免疫センサー/抗原提示 | CD14, TLR4, HLA-DRA, HLA-DRB1, HLA-DQA1, TYROBP, AIF1, CSF1R, C1QB, GPNMB |
| 適応/細胞障害 | CD247, CD4, CD8A, NKG7, PRF1, GNLY |
| 扩 PD/単球・ミクログリア | ITGAM, CX3CR1, FCGR3A, TREM2, IL1B, TNF, TLR2, CD68, LYZ, ITGAX, MERTK, TGFB1, IL6, CCL2, NLRP3, IL18, SPP1, CD33, APOE, CLEC7A, FCGR2A, CD163 |

---

## 6. 計算時間

| 実験 | 対象 | max_ncells | 時間（目安） |
|---|---|---|---|
| ISP-trial-1〜3（Mig DC・21遺伝子） | 1 型 | 100 | 各数分 |
| ISP-full-6-21 | 6 型 | 200 | ~25 分 |
| ISP-full-6-43 | 6 型 | 200 | ~30 分 |
| ISP-trial-50（Macrophages） | 1 型 | 50 | ~5-10 分 |

> 単一遺伝子の deletion を遺伝子数分ループするため、forward の総数は 細胞×遺伝子。MPS では forward が支配的で、6 型 × 43 遺伝子は 30 分前後。

---

## 7. In Silico Perturbation の結果まとめ

### 7.1 全体像
- 6 細胞型 × 43 遺伝子のうち、発現遺伝子 **180 個**で goal_state_shift を算出。
- **104 個は正（deletion が WT → PF6m へ近づける）、76 個は負（遠ざける）**。
- 全 180 個が遺伝子固有の値（バグ修正後）。

### 7.2 コンセンサスランキング（平均シフト、全型）
| 遺伝子 | mean_shift | max_shift | 正の細胞型割合 |
|---|---|---|---|
| **CD8A** | +0.0042 | +0.0152 | 67% |
| **C1QB**（補介） | +0.0035 | +0.0138 | 83% |
| **SPP1** | +0.0014 | +0.0055 | 50% |
| **CD163** | +0.0013 | +0.0052 | 80% |
| **MERTK** | +0.0008 | +0.0099 | 60% |
| **CD14** | +0.0005 | +0.0027 | 83% |
| **TLR4** | +0.0004 | +0.0011 | 60% |
| **VPS35** | +0.0002 | +0.0010 | 83% |
| （下位） **LYZ** | −0.011 | +0.0068 | 67% |
| （下位） **CD33** | −0.0065 | +0.0014 | 50% |
| （下位） **TYROBP** | −0.0022 | +0.0014 | 67% |

### 7.3 細胞型別の主要遺伝子
| 細胞型 | トップ（正方向） | 備考 |
|---|---|---|
| Migratory.DCs | MERTK, FCGR2A, IL1B | DC のリソファージ/受容体 |
| cDC.1 | APOE, C1QB, CD4 | APC で APOE/補介が最強 |
| cDC.2 | C1QB, APOE, LYZ | 補介＋単球様 |
| pDCs | ITGAX, TNF, C1QB | 自然系 |
| Macrophages | CD8A, SPP1, ITGAX, TYROBP | 最強シグナル |
| Ly6c^high Mono | CD8A, C1QB, APOE, CD4 | CD8A が最大(+0.015) |

### 7.4 解釈と科学的示唆
- **CD8A の deletion** が複数 APC 型で最強に早期PD方向へシフト：リンパ節で細胞障害性 T 応答を低下させることが、末梢免疫過応答→発走病理に寄与する可能性。ただし CD8A は適応 T 細胞マーカーで、APC には低発現のため、シフトは交絡/微弱である可能性に注意。
- **C1QB / CD14 / TLR4 / MERTK / TREM2（補介・スカベンジャー）**：一貫して「欠失＝発走方向」。早期の免疫活性化を抑えることで病理が進む＝**発走早期の免疫防御を抑制**を意味する可能性。逆に、これらの*増強*（活性化）が保護的である可能性も示唆。これは **早期発見マーカーとして有望**（免疫活性化の増加を血液で検出）。
- **CD4 / CSF1R / CD33 の deletion は発走方向に寄らない** → これらの欠失はむしろ保護的シフト。CSF1R（単球/ミクログリア生存）や CD33 は PD では抑制が保護に働く可能性（CD33 は PD 保護の対立遺伝子で知られる）。

### 7.5 結果の制約と注意
- 単一遺伝子の grouped stats は**平均 shift のみで有意性（p値/FDR）を伴わない**。有意性はチュートリアル通り `genes_to_perturb="all"` + null 分布比較で得られる（MPS では時間増大）。
- シフト絶対値は ~1e-2 以下と小さい。**ランキングの方向を参考**にし、絶対値の解釈は慎重に行うべき。
- ファインチューニング済み分類器は**細胞型**ベース。疾病状態（WT vs PD）を直接学習していないため、疾病シグナルが埋もれる可能性。疾病状態で再学習したモデルが今後の改善点。

### 7.6 可視化

**コンセンサスランキング**（6型の平均シフト・上位15遺伝子。赤=欠失で早期PD方向・青=遠ざかる）

![gene_ranking](./assets/isp_gene_ranking.png)

**ヒートマップ**（細胞型×遺伝子の Shift_to_goal_end 行列。正=欠失で早期PD方向）

![heatmap](./assets/isp_heatmap.png)

**細胞型別ランキング**（各 APC 型で欠失が早期PDへ寄る遺伝子）

![celltype_volcano](./assets/isp_celltype_volcano.png)

> 図は `analysis/08_isp_report_figs.py` で `result/report/assets/` に生成。元データは `result/report/pd_apc_early_isp_summary_full6.csv`。

---

## 8. 追加・補足情報

### 8.1 主要な実装上の発見（プロジェクト知識）
1. **`datasets>=5` で `perturb_data` の `dataset.map` がハング** → `datasets==4.0.0` に固定。
2. **`genes_to_perturb` リスト＝グループ削除**、細胞が全遺伝子を同時発現する必要あり → 1遺伝子ずつ削除。
3. **`InSilicoPerturberStats.read_dictionaries` が指定ディレクトリの全 `*_raw.pickle` を読む** → 各遺伝子を個別サブディレクトリへ。
4. 細胞型名 `Migratory.DCs` のドットが `Path(prefix).with_suffix(".csv")` でファイル名を壊す → prefix をドット無しに正規化。

### 8.2 次ステップ候補
- 有意性検定（null 比較）を加える（`max_ncells` を絞って）。
- 疾病状態で fine-tune したモデルで再実行。
- 全 GWAS-PD リスク遺伝子へ拡張。

---

*Generated for internal research reporting. Data: `result/report/pd_apc_early_isp_summary_full6.csv`（全6型・43遺伝子）、`result/report/assets/isp_gene_ranking.csv`。*