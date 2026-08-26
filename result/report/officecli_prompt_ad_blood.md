# officecli（PPT 生成）用プロンプト — AD_blood

このファイルは、`result/report/report_ad_blood.md`（早期アルツハイマー病の in silico gene deletion 解析レポート・血液）を **PowerPoint スライド** に変換するための officecli への指示書です。

---

## 使用方法（officecli コマンド）

```bash
# スライド作成ガイドを読み込み
officecli load_skill pptx

# 新規デッキ作成
officecli create report_ad_blood.pptx

# スライド検証・保存（必須）
officecli view report_ad_blood.pptx issues
officecli save report_ad_blood.pptx
```

> 詳細なスキーマは `officecli help pptx shape` で確認。画像は
> `result/report/figures/ad_blood_isp_rank_3m.png`
> `result/report/figures/ad_blood_isp_timeseries.png`
> `result/report/figures/ad_blood_isp_heatmap.png` を使用。

---

## PPT レイアウト設計

1. **タイトルスライド** — 疾患（AD）、臓器（血液）、手法（in silico gene deletion, Geneformer V2-104M）
2. **対象・データセット** — AD_blood、48,909 細胞、15 細胞型、AD vs WT、時点 3/4.5/6/9/12m
3. **背景知識（血液特異的）** — 末梢免疫と AD 早期化、Gate 2020（CD8 T クローン）、Tan 2002（S100A8/9）、Nakamura 2018（血液 AD 検出）
4. **細胞・遺伝子の選択理由** — 血液免疫 10 種（低シグナル除外）、55 遺伝子（GWAS/炎症/適応免疫/骨髄球）、文献
5. **ISP 設定・計算パラメータ** — delete, combos=0, AD→WT, Fine-tuned CellClassifier V2-104M, emb=cls, max_ncells=200, 実行約1時間35分
6. **結果ランキング（図）** — `figures/ad_blood_isp_rank_3m.png` を大きく表示（CD3E +0.0142, CD4 +0.0138）
7. **経時傾向（図）** — `figures/ad_blood_isp_timeseries.png` と `ad_blood_isp_heatmap.png`
8. **結果表** — 上位12遺伝子と Shift_to_goal_end（3m/4p5m/6m/9m/12m）
9. **生物学的解釈・考察** — T 細胞軸（CD3E/CD4/CD74）が最有力、血液では APOE 中性、CD8A/MS4A1/NKG7 は負
10. **制限事項** — 統計的有意性なし、max_ncells=200、5xFAD モデルの一般化注意

---

## 推奨プロンプト（officecli に渡す文言）

```
result/report/report_ad_blood.md の内容を基に、10 枚の PowerPoint を作成して下さい。
各スライド:
- タイトル: 早期アルツハイマー病発症因子の in silico gene deletion — 血液（AD_blood）
- 2枚目: 対象データセット（AD_blood, 48,909細胞, 15細胞型, AD 25,222 vs WT 23,687, 時点3/4.5/6/9/12m）
- 3枚目: 背景知識（末梢免疫と AD 早期化）と注目細胞（単球/マクロファージ系, T細胞, B細胞, 好中球, DC）・遺伝子クラス、文献（Gate2020, Tan2002, Nakamura2018, Guerreiro2013）
- 4枚目: 細胞(10種, 低シグナル除外)・遺伝子(55)選択理由とカテゴリ（AD GWAS 12, 炎症24, 適応免疫12, 骨髄球10）
- 5枚目: ISP設定(delete, combos=0, AD→WT, Fine-tuned CellClassifier V2-104M, emb_mode=cls, max_ncells=200, forward_batch=64, nproc=1)と実行時間約1時間35分
- 6枚目: 結果ランキング図（figures/ad_blood_isp_rank_3m.png）を挿入（CD3E +0.0142, CD4 +0.0138, CD74 +0.0053, S100A8 +0.0050, C1QC +0.0049, CD19 +0.0041）
- 7枚目: 経時傾向図（figures/ad_blood_isp_timeseries.png, ad_blood_isp_heatmap.png）— T 細胞系が全時点で安定して正
- 8枚目: 上位12遺伝子表（Gene, カテゴリ, 3m/4p5m/6m/9m/12m Shift）を掲載
- 9枚目: 生物学的解釈（T細胞活性化/共刺激シグナルが最早期ドライバー, 血液ではAPOE中性, CD8A/MS4A1/NKG7削除は逆効果=WT維持因子）
- 10枚目: 制限事項（統計的有意性なし, max_ncells=200, 5xFADは家族性AD, fine-tunedモデルはcelltype分類器）
配色はクリーンで視認性の高いもの（ブルー#1f6feb系アクセント）。テキストは箇条書きを基本とし、過剰な装飾は避ける。
```
