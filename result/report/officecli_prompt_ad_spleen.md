# officecli（PPT 生成）用プロンプト — AD_spleen

このファイルは、`result/report/report_ad_spleen.md`（早期アルツハイマー病の in silico gene deletion 解析レポート・脾臓）を **PowerPoint スライド** に変換するための officecli への指示書です。

---

## 使用方法（officecli コマンド）

```bash
# スライド作成ガイドを読み込み
officecli load_skill pptx

# 新規デッキ作成
officecli create report_ad_spleen.pptx

# スライド検証・保存（必須）
officecli view report_ad_spleen.pptx issues
officecli save report_ad_spleen.pptx
```

> 詳細なスキーマは `officecli help pptx shape` で確認。画像は
> `result/report/figures/ad_spleen_isp_rank_3m.png`
> `result/report/figures/ad_spleen_isp_timeseries.png`
> `result/report/figures/ad_spleen_isp_heatmap.png` を使用。

---

## PPT レイアウト設計

1. **タイトルスライド** — 疾患（AD）、臓器（脾臓）、手法（in silico gene deletion, Geneformer V2-104M ファインチューン済み）
2. **対象・データセット** — AD_spleen、76,947 細胞、25 細胞型、AD 38,311 vs WT 38,636、早期時点 3/4.5/6m
3. **背景知識（脾臓特異的）** — 脾臓-脳軸（spleen–brain axis）、骨髄系（単球/マクロファージ/DC/cDC1）が AD リスク遺伝子を発現、Gate 2020・Keren-Shaul 2017
4. **細胞・遺伝子の選択理由** — 脾臓免疫プール 19 種（骨髄系優先・ストロマ/血球系ノイズ除外）、55 遺伝子（血中 07_ad_spleen_early_isp と同じパネル）、文献
5. **ISP 設定・計算パラメータ** — delete, combos=0, AD→WT, Fine-tuned CellClassifier V2-104M（25クラス, acc0.93）, emb=cls, max_ncells=200, 実行約50分
6. **結果ランキング（図）** — `figures/ad_spleen_isp_rank_3m.png` を大きく表示（CLEC9A +0.0066, ITGA2B +0.0062, CD19 +0.0034）
7. **経時傾向（図）** — `figures/ad_spleen_isp_timeseries.png` と `ad_spleen_isp_heatmap.png`
8. **結果表** — 上位12遺伝子と Shift_to_goal_end（3m/4p5m/6m）
9. **生物学的解釈・考察** — cDC1（CLEC9A）・血小板（ITGA2B）・B細胞（CD19）軸が最有力、APOE は中性、CD8A は強く負
10. **臓器横断比較・制限事項** — 血中 07_ad_spleen_early_isp との対比（CLEC9A/ITGA2B/CD19 一貫、S100A8 逆符号）、統計的有意性なし・max_ncells=200

---

## 推奨プロンプト（officecli に渡す文言）

```
result/report/report_ad_spleen.md の内容を基に、10 枚の PowerPoint を作成して下さい。
各スライド:
- 1枚目: タイトル「早期アルツハイマー病発症因子の in silico gene deletion — 脾臓（AD_spleen）」
- 2枚目: 対象データセット（AD_spleen, 76,947細胞, 25細胞型, AD 38,311 vs WT 38,636, 早期時点3m/4.5m/6m, 5xFADファミリー）
- 3枚目: 背景知識（脾臓=最大の二次リンパ器官・spleen-brain axis）と注目細胞（単球/マクロファージ系, DC/cDC1, T細胞, B細胞, 好中球）・遺伝子クラス、文献（Gate2020, Keren-Shaul2017, Guerreiro2013, Nakamura2018）
- 4枚目: 細胞(19種, 骨髄系優先, ストロマ/赤血球・巨核球系ノイズ除外)・遺伝子(55, 血中07_ad_spleen_early_ispと同一パネル)選択理由とカテゴリ（AD GWAS 12, 炎症24, 適応免疫9, 骨髄球10）
- 5枚目: ISP設定(delete, combos=0, AD→WT, Fine-tuned CellClassifier V2-104M 25クラス, emb_mode=cls, max_ncells=200, forward_batch=64, nproc=1)と実行時間約50分
- 6枚目: 結果ランキング図（figures/ad_spleen_isp_rank_3m.png）を挿入（CLEC9A +0.0066, ITGA2B +0.0062, CD19 +0.0034, CD4 +0.0034, CLU +0.0029）
- 7枚目: 経時傾向図（figures/ad_spleen_isp_timeseries.png, ad_spleen_isp_heatmap.png）— CLEC9A/ITGA2B/CD19 が安定して正
- 8枚目: 上位12遺伝子表（Gene, カテゴリ, 3m/4p5m/6m Shift）を掲載
- 9枚目: 生物学的解釈（cDC1・血小板・B細胞軸が最早期ドライバー候補, CLU/CD33/SORL1/TREM2が正, APOEは中性, CD8A/LYZは削除逆効果=WT維持因子, 4.5mでCD4反転）
- 10枚目: 臓器横断比較（血中07_ad_spleen_early_isp: CLEC9A/ITGA2B/CD19両方正, CD3E/CD4は血中で突出, S100A8は逆符号, CD8A両方強く負）と制限事項（統計的有意性なし, max_ncells=200, 5xFADは家族性AD, fine-tunedはcelltype分類器）
配色はクリーンで視認性の高いもの（ブルー#1f6feb系アクセント vs 血中07_ad_spleen_early_ispと区別できるエメラルド系#16a085をタイトルアクセントにするなど）。テキストは箇条書きを基本とし、過剰な装飾は避ける。
```

> 図のパスは officecli 実行時、`figures/` が `result/report/figures/` を指すように調整してください（officecli の CWD に応じて `result/report/figures/ad_spleen_isp_rank_3m.png` 等に置換可）。