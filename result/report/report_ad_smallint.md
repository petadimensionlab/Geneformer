# 早期アルツハイマー病（AD）発症因子の探索 — in silico gene deletion レポート（小腸 / AD_smallint）

> 実験スクリプト: `analysis/07_ad_smallint_early_isp.py`
> データセット: `input/AD_smallint`（5xFAD マウス小腸 scRNA-seq、ヒトオルソログ変換済み）

---

## 1. 対象

| 項目 | 値 |
|---|---|
| 対象疾患 | **AD（アルツハイマー病）** |
| 対象臓器 | **小腸（small intestine）** |
| データセット | 5xFAD マウス小腸、61,951 細胞、25 細胞型、30 個体 |
| 疾患状態 | AD (31,071) vs WT (30,880) |
| 時点 | 3・4.5・6・9・12 ヶ月齢（`samples4`） |
| **解析対象時点** | **早期（3–4.5–6 ヶ月）のみ** |

実験対応（`result/report/experiments.csv`）: `isp_early_ad_smallint`

---

## 2. 臓器特異的に探索すべき背景知識（小腸）

AD 研究の多くは脳に注目するが、**腸-脳軸（gut–brain axis）** が AD の病態・早期発見に重要な役割を持つという証拠が蓄積している。小腸は病原体・内毒素に対する第一のバリアであり、全身の炎症・免疫状態に直結する。

### 注目すべき細胞
- **腸上皮 / バリア細胞（IECs・Paneth・Goblet・Epithelial）**: 粘液層とタイトジャンクションを形成し、腸管透過性を制御。AD では透過性亢進→内毒素(リポ多糖)の血中移行→全身性炎症が示唆される。
- **組織常在マクロファージ・単球（Macrophages・Ly6c Monocytes）**: 炎症性サイトカイン産生、感染・損傷応答。
- **T 細胞（CD4/CD8）・NK（NK_ILC1）**: 全身性免疫応答の中心。AD 早期に末梢免疫プロファイルが変化するという報告がある。

### 注目すべき遺伝子の機能クラス
- **粘液層（Mucin）**: MUC2, MUC1, MUC3A, MUC13, MUC5B, MUC4 — 腸バリアの物理防御
- **タイトジャンクション（TJ）**: CLDN, OCLN, TJP1/2/3, F11R, JAM2 — 透過性制御
- **抗菌ペプチド**: LYZ, DEFB1, REG3G, DEFA6, CAMP
- **幹細胞/再生/化生**: LGR5, OLFM4, SPP1, TFF3, AGR2
- **免疫・炎症**: IL1B, TNF, TLR4, S100A8/9, CD68

### 参考文献（腸-脳軸と AD）
| 文献 | 要点 |
|---|---|
| **Vogt et al., *Cell* 170:1253 (2017)** | AD 患者・モデルで腸内細菌組成変化を報告。腸-脳軸の生物学的基盤 |
| **Kimura et al., *eLife* 10:e64683 (2021)** | AD 前駆期に腸内細菌由来の代謝物変化が先行する可能性 |
| **Kowalski & Mulak 2019** | 腸-脳軸と神経変性・透過性変化のレビュー |
| **Cryan et al., *Physiol Rev* 99:1877 (2019)** | 微生物-腸-脳軸の生理的役割の包括的レビュー |

### 3. AD のリスク遺伝子と免疫シグネチャー
AD の GWAS リスク遺伝子は **免疫・ミクログリア機能** に濃縮される（Jansen et al. *Nat Genet* 2019; Wightman et al. 2021）。その多くは小腸でも発現するため、早期発見マーカーの手がかりになる。

---

## 3. 細胞・遺伝子の選択理由

### 選択した細胞型（8種、早期のみ）
小腸内の **免疫細胞（Macrophages, CD4.T.cells, CD8.T.cells, NK_ILC1）＋ バリア・上皮細胞（IECs, Paneth.cells, Goblet.cells_M.cells, Epithelial.cells）** をプール。

- 腸バリアと免疫が AD の早期全身性炎症で重要とされるため。
- 各細胞型に早期 AD/WT が十分な細胞数（後述）あることを確認して選択。

### 選択した遺伝子（49種、仮説駆動・固定リスト）
| カテゴリ | 遺伝子 | 選択根拠（文献） |
|---|---|---|
| **AD risk (GWAS)** | APOE, TREM2, CLU, BIN1, CD33, ABCA7, SORL1, TYROBP, INPP5D, PLCG2, MEF2C | Lambert et al. 2013; Jansen 2019; Wightman 2021; Kunkle 2019 |
| **炎症/免疫** | C1QA, C1QB, CD68, CSF1R, AIF1, CX3CR1, IL1B, TNF, TLR4, LYZ, APOC1, HLA-DRA | 神経炎症・免疫応答（Heneka 2015; 補体） |
| **粘液バリア** | MUC2, MUC1, MUC3A, MUC13, MUC17, MUC12, MUC5B, MUC4 | 腸バリア（Bergstrom 2017; Johansson 2011） |
| **タイトジャンクション** | OCLN, CLDN1-15, TJP1-3, F11R, JAM2 | 透過性（Ulluwishewa 2011） |
| **抗菌/分泌** | LYZ, DEFB1, DEFB103A, REG3G, DEFA6, CAMP | 抗菌ペプチド（Bevins 2004） |
| **幹細胞/化生** | LGR5, OLFM4, SPP1, TFF3, CLCA1, S100A8, S100A9 | 上皮再生・化生（Barker 2010; 5xFAD 腸で S100A8/9 上昇報告） |

> **文献**: 上記は代表例。GWAS（Lambert 2013; Wightman 2021）、免疫炎症（Kinney 2018; Heneka 2015）、腸バリア（Bergstrom 2017）を背景とした仮説駆動選択である。

---

## 4. in silico perturbation の主要設定

| 設定項目 | 値 |
|---|---|
| **遺伝子操作** | **delete（仮想遺伝子削除）** |
| 対象状態（disease） | `start_state = AD`, `goal_state = WT` |
| state_key | `disease` |
| alt_states | なし |
| 細胞型フィルタ | 免疫＋上皮 8 細胞型（早期のみ） |
| 時点フィルタ | 早期（AD/WT × 3/4.5/6m） |
| combos | 0（1遺伝子ずつ） |
| モデル | Geneformer V2-104M（`model_version="V2"`） |
| 埋め込み | `emb_mode="cls"`, `cell_emb_style="mean_pool"` |
| 分類器 | ファインチューン済みセル分類器（`TRAINED_MODEL_PATH`） |
| 統計モード | `goal_state_shift`（`Shift_to_goal_end`） |

> **解釈**: 各遺伝子を早期 AD 細胞から削除し、embedding が早期 WT（健常）へどれだけ近づくかを `Shift_to_goal_end`（cosine シフト）で評価。**正に大きい = その遺伝子の削除が AD を WT へ戻す = 早期駆動/介入標的候補**。

---

## 5. 計算パラメータ・実行時間

| パラメータ | 値 |
|---|---|
| 解析対象細胞数（早期プール） | **17,949 細胞** |
| 状態埋め込み（EmbExtractor） | `max_ncells=1000` |
| perturb（InSilicoPerturber） | `max_ncells=300`, `forward_batch_size=64`, `emb_layer=0` |
| 遺伝子数 | 候補 56 → 発現 49（early プールに発現あり） |
| nproc | 1（`datasets==4.0.0` 必須、多プロセスは不安定） |
| デバイス | CUDA |
| **実行時間** | **約 45 分**（15:38–16:23, 2026-08-25） |

> **補助**: `datasets>=5` だと `InSilicoPerturber.perturb_data` がハングするため `datasets==4.0.0` に固定。`max_ncells` は RAM 見積もり（`estimate_perturb_ram`）で調整。プール細胞に発現しない 7 遺伝子（CX3CR1, HLA-DRA, MUC3A, MUC17, MUC12, DEFB103A, DEFA6）は warning 付きで除外。

---

## 6. 結果まとめ（in silico deletion — 早期 AD→WT シフト）

`Shift_to_goal_end` の降順ランキング（正 = 削除で早期 AD が WT へ近づく）。全 49 遺伝子の値は `result/report/data/ad_smallint_early_isp_results.csv` を参照。

### 上位ヒット（正のシフト = 早期駆動/介入標的候補）
| Rank | Gene | 分類 | Shift_to_goal_end |
|---|---|---|---|
| 1 | **CLDN1** | タイトジャンクション | **+0.00926** |
| 2 | **MUC5B** | 粘液バリア | **+0.00377** |
| 3 | **TREM2** | 免疫/ミクログリア | **+0.00291** |
| 4 | **SORL1** | AD risk | **+0.00234** |
| 5 | **MUC2** | 粘液バリア | **+0.00184** |
| 6 | **TYROBP** | 免疫/ミクログリア | **+0.00157** |
| 7 | **C1QA** | 補体 | +0.00087 |
| 8 | **MUC4** | 粘液バリア | +0.00054 |
| 9 | **TNF** | 炎症 | +0.00049 |
| 10 | **MEF2C** | AD risk | +0.00046 |

### 中位（ほぼ中性、-0.001〜+0.0004）
CSF1R, JAM2, S100A9, CLU, BIN1, IL1B, PLCG2, AGR2, TJP1, INPP5D, CLDN3, ABCA7, DEFB1, OCLN, TLR4, TJP2, F11R, TFF3, CLCA1, AIF1, LGR5, CLDN4, REG3G, CLDN15, SPP1, CAMP, TJP3, CD68, CD33 …

### 下位（負 = 削除で WT から遠ざかる）
| Gene | 分類 | Shift_to_goal_end |
|---|---|---|
| OLFM4 | 幹細胞 | -0.01029 |
| LYZ | 抗菌 | -0.00841 |
| APOC1 | 免疫/リポ蛋白 | -0.00763 |
| C1QB | 補体 | -0.00593 |
| MUC13 | 粘液 | -0.00211 |
| CLDN2 | タイトジャンクション | -0.00170 |

### 生物学的解釈
- **腸バリア遺伝子（CLDN1・MUC5B・MUC2・MUC4）と免疫/ミクログリア遺伝子（TREM2・TYROBP・C1QA・TNF）** を早期 AD 細胞で除去すると **WT 方向へ戻る**ことが示唆された。特に **CLDN1** が最上位で、小腸の腸バリア整合性が早期 AD の重要な駆動・標的部位である可能性を示す。
- **下位（負）** の遺伝子（OLFM4, LYZ, APOC1, C1QB）は、むしろこれらが存在することで AD 細胞が WT に近い状態を保っている（＝削除すると AD 側へシフト）ことを示し、保護的または AD 特異的な発現状態に関与する可能性。

> 図: `result/report/figures/ad_smallint_isp_rank.png`, `figures/ad_smallint_vs_brain_top.png`

---

## 7. 制限事項・留意点
- `Shift_to_goal_end` は cosine シフトで、**効果の統計的有意性（p値）は本実装の stats 出力に含まれない**。順位付けはシフトの大きさに基づく。
- 複数細胞型のプールにより、細胞型特異的信号は平均化される（個別解析はより明瞭に区別可能）。
- 5xFAD マウス（ヒトオルソログ変換）は AD の家族性モデルであり、ヒト散発性 AD の一般化には検討が必要。