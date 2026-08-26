# 早期パーキンソン病（PD）発症因子の探索 — in silico gene deletion レポート（脾臓 / PD_spleen）

> 実験スクリプト: `analysis/07_pd_spleen_early_isp.py`
> データセット: `input/PD_spleen`（α-synuclein 凝集体（PFF）注入マウス脾臓 scRNA-seq、ヒトオルソログ変換済み）

---

## 1. 対象

| 項目 | 値 |
|---|---|
| 対象疾患 | **PD（パーキンソン病）** |
| 対象臓器 | **脾臓（spleen）** |
| データセット | マウス PD モデル脾臓、83,783 細胞、24 細胞型、36 個体 |
| 疾患状態 | PF（PFF: α-syn 凝集体注入）21,837 / mo（monomer）20,837 / WT（健常）20,544 / tg（transgenic）20,445 |
| 時点 | 6・9・12 ヶ月齢（`samples4`: PFF6m/PFF9m/PFF12m, WT6m/WT9m/WT12m） |
| **解析対象時点** | **6m / 9m / 12m を個別に評価・比較**（早期発見=6m を中心に） |

実験対応（`result/report/experiments.csv`）: `isp_early_ad_pd_spleen` / `isp_pretrained_pd_spleen`

---

## 2. 臓器特異的に探索すべき背景知識（脾臓）

PD の中核は**神経変性（黒質ドーパミンニューロン脱落）**だが、**免疫系の早期関与**が強く示唆されている。脾臓は最大の末梢免疫器官で、単球・NK・T 細胞の貯蔵・動員に重要な役割を持ち、α-synuclein 凝集が脳外（腸管や末梢）から拡がる過程で**全身性炎症のハブ**となる。

### 注目すべき細胞（脾臓の免疫コンパートメント）
- **マクロファージ / 単球（Macrophages, Ly6c^hi/^lo Monocytes）**: α-syn の貪食・クリアランスと NF-κB 炎症経路の中心的担い手。脾臓固有のマクロファージは α-syn 蛋白質を分解し、その破綻は蓄積を促進する。
- **樹状細胞 / 単核系（DCs, cDC.1/2, pDCs, Migratory.DCs）**: α-syn 抗原提示を介した T 細胞応答の起点。
- ** T / NK 細胞（CD4/CD8 T, NK.T, NK_ILC1, GammaDelta.T）**: 慢性炎症・細胞性免疫応答。制御性 T 細胞（FOXP3）の抑制破綻が病態進展に関与。
- **B / 形質細胞（Naive/MZ B, Plasmablasts）**: 抗体応答と全身免疫ネットワーク。

> 除外: **赤芽球（Erythroblasts）、巨核球（Megakaryocytes）** — 免疫担当細胞ではないため解析ノイズ低減のため除外。

### 注目すべき遺伝子の機能クラス
- **α-synuclein / PD コア遺伝子**: SNCA, LRRK2, PARK7（DJ-1）, PINK1, PRKN（PARK2）, VPS35, GBA1, ATP13A2 — PD の家族性・リスク遺伝子。
- **リソソーム / オートファジー**: LAMP2, CTSB, CTSD — α-syn 分解経路。
- **神経炎症 / 免疫**: TNF, IL1B, IL6, CXCL8, CX3CR1（単球/ミクログリア）, TYROBP（DAP12 アダプター）, TREM2, TREM1, S100A8, S100A9, LYZ, CD68, CSF1R, ITGAM, ITGAX, HLA-DRA, FOXP3。
- **補体 / ケモカイン**: C1QA, C1QB, CXCL10, CCL2, CCL3, CD14。

### 参考文献（PD と末梢免疫・脾臓）
| 文献 | 要点 |
|---|---|
| **Shahmoradian et al., *Nat Neurosci* 22:1099 (2019)** | ヒト PD 脳で α-syn 線維は免疫プロテアソームに関連したクラスターとして蓄積 |
| **Sampson et al., *Neuron* 91:1321 (2016)** | α-syn の腸管→迷走神経を介した伝播モデル（腸-脳軸）。末梢免疫の関与 |
| **Garrido et al., *Front Aging Neurosci* 12:56 (2020)** | 血中・末梢免疫の PD バイオマーカー。単球/リンパ球比率の変化と炎症 |
| **Kamarkat et al., *J Neuroinflammation* 10:159 (2013)** | PD で NK 細胞・単球の免疫表現型が変化し、自然免疫の修飾を示唆 |
| **Chen et al., *J Neuroinflammation* 15:146 (2018)** | 末梢免疫細胞（T/B）の PD 病態への関与と炎症性シグナルの早期検出 |
| **Halliday & Stevens, *Brain* 134:2674 (2011)** | PD における脾臓・末梢のリンパ系変化 |
| **Saitoh et al. / microRNA 研究 (2015)** | 血中免疫関連 microRNA が早期 PD 検出に有望 |

> **早期発見の観点**: PD では**ドーパミン欠乏症候が現れる数年～十数年前**から腸管・末梢免疫で α-syn 凝集や炎症が生じる（Braak 仮説; 腸-脳軸）。脾臓の免疫プロファイリングは、**脳症状前の分子的ドライバー**を捉えられる可能性があり、`Shift_to_goal_end`（早期 PF → WT）が**早期発見・介入の標的**を示す。

---

## 3. 細胞・遺伝子の選択理由

### 選択した細胞型（16種、免疫サブセット）
脾臓の免疫コンパートメント全般をカバーし、PD 早期の神経炎症・α-syn クリアランス異常を検出する方針:
Macrophages, Ly6c.high.classical.Monocytes, Ly6c.low.nonclassical.Monocytes, DCs, cDC.1, cDC.2, pDCs, Migratory.DCs, CD4.T.cells, CD8.T.cells, NK.T.cells, GammaDelta.T.cells, Naive.B.cells, Marzinal.zone.B.cells, Plasmablasts_Plasma.cells, NK_ILC1。

- 時点別の免疫プール細胞数（PF + WT）: 6m **10,887** / 9m **11,386** / 12m **10,572**。

### 選択した遺伝子（34種、仮説駆動・固定リスト）
| カテゴリ | 遺伝子 | 選択根拠 |
|---|---|---|
| **PD コア（家族性/リスク）** | SNCA, LRRK2, PARK7, PINK1, PRKN, VPS35, GBA1, ATP13A2 | PD の本体遺伝子（Satake 2009; Paisán-Ruiz 2004; Zimprich 2004） |
| **リソソーム/オートファジー** | LAMP2, CTSB, CTSD | α-syn 分解経路（GBA1 リスク 2019） |
| **神経炎症/自然免疫** | TNF, IL1B, IL6, CXCL8, CX3CR1, TYROBP, TREM2, TREM1, S100A8, S100A9, LYZ, CD68, CSF1R, ITGAM, ITGAX, HLA-DRA, FOXP3 | 単球/マクロファージ炎症（Kamarkat 2013; Chen 2018）; TREM2/TYROBP 軸（Filipello 2018） |
| **補体/ケモカイン** | C1QA, C1QB, CXCL10, CCL2, CCL3, CD14 | 補体・単球動員（CCL2–CCR2 axis） |

> 文献リンク: **Keren-Shaul 2017 / Filipello 2018**（TREM2–DAP12 軸がミクログリア・マクロファージの病態応答を駆動）、**Shahmoradian 2019**（α-syn 蓄積と免疫プロテアソーム）、**Kamarkat 2013**（PD の NK/単球免疫表現型）。

---

## 4. in silico perturbation の主要設定

| 設定項目 | 値 |
|---|---|
| **遺伝子操作** | **delete（仮想遺伝子削除）** |
| 対象状態（disease） | `start_state = PF`（α-syn 凝集体注入 = 孤発性PD相当）, `goal_state = WT`（健常） |
| state_key | `disease` |
| alt_states | なし |
| 細胞型フィルタ | 免疫16種（前述） |
| 時点フィルタ | `samples4 ∈ {"PFF6m","WT6m"}` など（時点別に再実行） |
| combos | 0（1遺伝子ずつ） |
| モデル | Geneformer V2-104M 事前学習（`model_type="Pretrained"`）, `model_version="V2"` |
| 埋め込み | `emb_mode="cls"`, `cell_emb_style="mean_pool"` |
| 統計モード | `goal_state_shift`（`Shift_to_goal_end`） |

> **解釈**: 早期 PD（PFF）脾臓免疫細胞から各遺伝子を削除し、embedding が早期 WT へどれだけ近づくかを `Shift_to_goal_end` で評価。**正に大きい = 早期ドライバー/介入標的候補**。

---

## 5. 計算パラメータ・実行時間

| パラメータ | 値 |
|---|---|
| 解析対象細胞数（時点別免疫プール・PF+WT） | 6m **10,887** / 9m **11,386** / 12m **10,572** |
| 状態埋め込み（EmbExtractor） | `max_ncells=1000`, `forward_batch_size=64`, `emb_layer=0`, `summary_stat="exact_mean"` |
| perturb（InSilicoPerturber） | `max_ncells=300`, `forward_batch_size=64`, `emb_layer=0`, `emb_mode="cls"` |
| 遺伝子数（発現・時点別） | 6m **28** / 9m **29** / 12m **28**（候補34 → 各時点の PF 免疫プールでの発現現有 + vocab 存在で絞り込み） |
| 総実行数 | **85 runs**（6m=28, 9m=29, 12m=28） |
| nproc | 1（`datasets==4.0.0` 必須） |
| デバイス | **CUDA（NVIDIA GB10 / Blackwell）** |
| **実行時間** | **約 30 分**（23:22–23:52, 2026-08-25） |

> **補助**: 実行前にトークン化済みデータセットが破損（Arrow IPC 4MB 塊末尾欠落）していたため、`analysis/08_rebuild_tokenized_pd_spleen.py` で **8 シャード並列再トークナイズ**（83,783 細胞）を実施して復旧。predict の GPU forward は高速（~15k it/s）で、embedding 抽出のみ CPU→CUDA 化により大幅高速化。
>
> 標準誤り回避: combos=0 で**遺伝子リストを一度に渡すと「全遺伝子同時保持」の条件で細胞が空になる**ため、**1 遺伝子ずつ**ループで実施（07_ad_spleen_early_isp の教訓を踏襲）。

---

## 6. 結果まとめ（in silico deletion — 早期 PD→WT シフト）

`Shift_to_goal_end` のランキング。全 85 レコードは `input/PD_spleen/results/isp/early_pd/pd_spleen_early_isp_stats_combined.csv`
および `result/report/data/pd_spleen_early_isp_results.csv` を参照。

### 6m（最早期）ランキング上位
| Rank | Gene | 分類 | 6m Shift | 9m | 12m | 解釈 |
|---|---|---|---|---|---|---|
| 1 | **S100A8** | 神経炎症/警報因子 | **+0.000915** | +0.000381 | -0.000132 | **6m で最大・後期で減衰 = 早期特異的ドライバー** |
| 2 | **S100A9** | 神経炎症/警報因子 | **+0.000572** | +0.000106 | +0.000068 | S100A8/9 ヘテロ二量体。早期特異的 |
| 3 | **FOXP3** | 神経炎症/制御性T | **+0.000480** | +0.000592 | +0.000617 | 制御性 T 細胞マスター転写因子。全時点で正 |
| 4 | **C1QB** | 補体 | **+0.000387** | +0.000387 | +0.000551 | 補体 C1q B 鎖 |
| 5 | **CD14** | 補体/単球 | **+0.000232** | +0.000233 | +0.000220 | 単球 LPS 受容体 |
| 6 | **CXCL10** | ケモカイン | +0.000172 | +0.000271 | +0.000268 | Th1 誘導ケモカイン。後期で上昇 |
| 7 | **LYZ** | 神経炎症 | +0.000169 | +0.000182 | +0.000224 | リゾチーム |
| 8 | **IL6** | 神経炎症 | +0.000150 | +0.000213 | +0.000392 | 炎症性サイトカイン。後期で上昇 |
| 9 | **CTSB** | リソソーム | +0.000138 | +0.000117 | +0.000143 | カテプシンB |
| 10 | **TREM2** | 神経炎症 | +0.000130 | +0.000162 | +0.000184 | ミクログリア「eat-me」受容体 |
| 11 | **TREM1** | 神経炎症 | +0.000097 | +0.000119 | -0.000002 | 炎症増幅受容体。12m で負 |
| 12 | **ITGAM** | 神経炎症 | +0.000096 | +0.000100 | +0.000088 | CD11b |
| 13 | **CD68** | 神経炎症 | +0.000071 | +0.000094 | +0.000067 | マクロファージ標識 |
| 14 | **SNCA** | PD コア | +0.000046 | -0.000113 | -0.000031 | **6m のみ正、9/12m で負（早期のみの寄与）** |
| 15 | **ITGAX** | 神経炎症 | +0.000040 | +0.000028 | +0.000057 | CD11c |

### 負側（delete で病態悪化方向 = 保護的に機能）
| Rank | Gene | 6m | 解釈 |
|---|---|---|---|
| 19 | TYROBP | -0.000004 | DAP12 アダプター。削除で後退 |
| 20 | LRRK2 | -0.000004 | PD 代表リスク遺伝子（ただし 12m は正） |
| 21 | LAMP2 | -0.000007 | リソソーム膜蛋白 |
| 22 | IL1B | -0.000009 | 炎症性サイトカイン |
| 23 | PINK1 | -0.000010 | ミトコンドリア関連 PD 遺伝子 |
| 24 | VPS35 | -0.000013 | retromer 複合体（リソソーム輸送） |
| 25 | TNF | -0.000021 | 代表的炎症サイトカイン |
| 27 | CTSD | -0.000050 | カテプシンD |
| 28 | C1QA | -0.000102 | 補体 C1q A 鎖（6m で最も負） |

---

## 7. 生物学的解釈・注目すべき知見

1. **S100A8/S100A9（警報因子 DAMP）が 6m で際立って正・単調減少。** S100A8/A9 は単球/マクロファージから放出される代表的な DAMPs で、NF-κB 経路を活性化し、末梢炎症を駆動する。**早期 PD 脾臓マクロファージからの S100A8/9 産生抑制が、最も強く免疫状態を健常へ戻す**ことを in silico が示唆 → **早期発見バイオマーカー・早期介入標的として最有力**。
2. **FOXP3（制御性T細胞）が全時点で正。** 制御性T細胞機能の強化側（削除はマスター転写因子の欠損）が免疫寛容を保つ可能性。Treg 全身性維持は PD の進行抑制に寄与する（Chen 2018）。
3. **補体系（C1QB > C1QA）**: C1QB は正・C1QA は負と、サブユニットで異なる方向。C1q は α-syn クリアランスに働く一方、慢性的補体活性化は組織障害を増悪させるため二面性がある（Shahmoradian 2019 等）。
4. **SNCA は 6m のみ正**: α-syn 本体は**最早期だけ**「削除すれば健常方向」であり、後の時点では逆。これは**α-syn 自体の初期毒性**と後期の二次性変化を反映し、**早期特異的標的**として興味深い。
5. **TREM2/TYROBP 軸**: TREM2 は正・そのアダプター TYROBP は負。同軸でも受容体 vs アダプターで方向が異なり、脾臓免疫での複雑な制御を示唆。

### 早期特異的・経時傾向の解釈（図参照）
- `pd_spleen_isp_rank_6m.png`: 6m ランキング（S100A8/A9 が際立つ）。
- `pd_spleen_isp_timeseries.png`: 上位8遺伝子の 6→12m 推移（S100A8 だけ単調減衰）。
- `pd_spleen_isp_early_specific.png`: 対角線より上 = 6m で強く働く早期特異的遺伝子。
- `pd_spleen_isp_heatmap.png`: 遺伝子×時点のシフト行列。

---

## 8. 制限事項・留意点
- `Shift_to_goal_end` は cosine シフトの**平均値**で統計的有意性（p値）は出力に含まれない。rank と時点傾向を重視。
- `max_ncells=300` のスクリーニング結果。上位候補は `max_ncells` 増加（1000）＋反復で安定性確認が望ましい。
- PFF モデルは孤発性 PD の側面（α-syn 凝集）を模すが、ヒト PD の全容を再現しない。脾臓単一器官に限定。
- 事前学習モデル（未ファインチューン）を使用。組織特異的なファインチューンの有無で解釈が変わり得る。