# officecli（PPT 生成）用プロンプト — PD_spleen

このファイルは、`result/report/report_pd_spleen.md`（早期パーキンソン病の in silico gene deletion 解析レポート）を **PowerPoint スライド** に変換するための officecli への指示書です。

---

## 使用方法（officecli コマンド）

```bash
# geneformer_hf の skill を読み込み、スライド作成ガイドを取得
officecli load_skill pptx

# 例: 新規デッキ作成
officecli create report_pd_spleen.pptx

# スライド追加例（タイトルスライド）
officecli add report_pd_spleen.pptx /slide[1] --type title --prop text="早期パーキンソン病発症因子の in silico gene deletion（脾臓）"

# スライドの検証・保存（必須）
officecli view report_pd_spleen.pptx issues
officecli save report_pd_spleen.pptx
```

> 詳細なスキーマは `officecli help pptx shape` で確認。画像は
> `result/report/figures/pd_spleen_isp_rank_6m.png` / `pd_spleen_isp_timeseries.png`
> / `pd_spleen_isp_heatmap.png` / `pd_spleen_isp_early_specific.png` を使用。

---

## PPT レイアウト設計

1. **タイトルスライド** — 疾患（PD）、臓器（脾臓）、手法（in silico gene deletion, Geneformer V2-104M）
2. **対象・データセット** — PFF α-syn 凝集体注入マウス脾臓、83,783 細胞、36 個体、PF vs WT、時点 6/9/12m
3. **背景知識（脾臓特異的）** — 腸-脳軸・末梢免疫と PD 早期、単球/マクロファージ/T/NK/B に注目
4. **細胞・遺伝子の選択理由** — 免疫16種、34遺伝子（PDコア/リソソーム/神経炎症/補体）、文献
5. **ISP 設定・計算パラメータ** — delete, combos=0, PF→WT, Pretrained V2-104M, max_ncells=300, 実行約30分
6. **結果ランキング（図）** — `figures/pd_spleen_isp_rank_6m.png` を大きく表示（S100A8/A9, FOXP3, C1QB）
7. **経時傾向（図）** — `figures/pd_spleen_isp_timeseries.png` と `pd_spleen_isp_early_specific.png`（早期特異性）
8. **結果表** — 上位遺伝子と Shift_to_goal_end（6m/9m/12m）
9. **生物学的解釈・考察** — S100A8/A9 警報因子の早期特異的シフト、FOXP3/Treg、補体の二面性
10. **制限事項** — 統計的有意性なし、max_ncells=300、PFF モデルの一般化注意

---

## 推奨プロンプト（officecli に渡す文言）

```
result/report/report_pd_spleen.md の内容を基に、10 枚の PowerPoint を作成して下さい。
各スライド:
- タイトル: 早期パーキンソン病発症因子の in silico gene deletion — 脾臓（PD_spleen）
- 2枚目: 対象データセット（PFF α-syn 凝集体注入マウス脾臓, 83,783細胞, PF vs WT, 時点6/9/12m）
- 3枚目: 背景知識（腸-脳軸・末梢免疫）と注目細胞（マクロファージ/単球, DC, T, NK, B）・遺伝子クラス
- 4枚目: 細胞(16種)・遺伝子(34)選択理由と文献（Shahmoradian2019, Sampson2016, Kamarkat2013, Chen2018）
- 5枚目: ISP設定(delete, combos=0, PF→WT, Pretrained V2-104M, emb_mode=cls, max_ncells=300, forward_batch=64)と実行時間約30分
- 6枚目: 結果ランキング図（figures/pd_spleen_isp_rank_6m.png）を挿入（S100A8 +0.000915, S100A9 +0.000572, FOXP3 +0.000480, C1QB +0.000387, CD14 +0.000232）
- 7枚目: 経時傾向図（figures/pd_spleen_isp_timeseries.png, pd_spleen_isp_early_specific.png）— S100A8 が 6m で最大・後期減衰
- 8枚目: 上位遺伝子表（Gene, Category, 6m/9m/12m Shift）を掲載
- 9枚目: 生物学的解釈（S100A8/A9 警報因子が早期特異的ドライバー, FOXP3/Treg が全時点で正, 補体C1QB/C1QAの方向差）
- 10枚目: 制限事項（統計的有意性なし, max_ncells=300, PFFは孤発性PD側面のみ, 事前学習モデル）
配色はクリーンで視認性の高いもの（紫#8e44ad系アクセント）。テキストは箇条書きを基本とし、過剰な装飾は避ける。
検証（view ... issues）と保存（save）を必ず行うこと。
```

---

## 図の参照

- `result/report/figures/pd_spleen_isp_rank_6m.png` — 6m（最早期）28遺伝子ランキング
- `result/report/figures/pd_spleen_isp_timeseries.png` — 上位8遺伝子の 6→12m 推移
- `result/report/figures/pd_spleen_isp_early_specific.png` — 6m vs 後期平均（早期特異性）
- `result/report/figures/pd_spleen_isp_heatmap.png` — 遺伝子×時点ヒートマップ
