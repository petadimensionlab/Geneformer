# officecli（PPT 生成）用プロンプト

このファイルは、`result/report/report_ad_smallint.md` と `result/report/report_ad_brain.md`
（早期 AD の in silico gene deletion 解析レポート）を **PowerPoint スライド** に変換するための
officecli への指示書です。

---

## 使用方法（officecli コマンド）

```bash
# geneformer_hf の skill を読み込み、スライド作成ガイドを取得
officecli load_skill pptx

# 例: 新規デッキ作成（小腸レポート）
officecli create report_ad_smallint.pptx

# スライド追加例（タイトルスライド）
officecli add report_ad_smallint.pptx /slide[1] --type title --prop text="早期AD発症因子のin silico gene deletion（小腸）"

# スライドの検証・保存（必須）
officecli view report_ad_smallint.pptx issues
officecli save report_ad_smallint.pptx
```

> 詳細なスキーマは `officecli help pptx shape` で確認。画像は
> `result/report/figures/ad_smallint_isp_rank.png` / `ad_brain_isp_rank.png` を使用。

---

## PPT レイアウト設計（両レポート共通）

1. **タイトルスライド** — 疾患（AD）、臓器、手法（in silico gene deletion, Geneformer V2-104M）
2. **対象・データセット** — 5xFAD、細胞数、時点（早期3–6m）、AD vs WT
3. **背景知識（臓器特異的）** — 腸-脳軸 / ミクログリア-DAM と文献
4. **細胞・遺伝子の選択理由** — 細胞型、49/23遺伝子、カテゴリ、文献
5. **ISP 設定・計算パラメータ** — delete, combos=0, AD→WT, max_ncells, 実行時間
6. **結果ランキング（図）** — `figures/ad_*_isp_rank.png` を大きく表示
7. **結果表** — 上位遺伝子と Shift_to_goal_end
8. **生物学的解釈・考察** — 主要ヒット（小腸: CLDN1 等 / 脳: APOE）
9. **制限事項** — 統計的有意性なし、5xFAD モデル、解像度

---

## 推奨プロンプト（officecli に渡す文言）

```
result/report/report_ad_smallint.md の内容を基に、8-9 枚の PowerPoint を作成して下さい。
各スライド:
- タイトル: 早期AD発症因子のin silico gene deletion — 小腸（AD_smallint）
- 2枚目: 対象データセット（5xFAD小腸, 61,951細胞, 早期3-6m, AD vs WT）
- 3枚目: 背景知識（腸-脳軸）と注目細胞・遺伝子
- 4枚目: 細胞(8種)・遺伝子(49)選択理由と文献（Vogt2017, Kimura2021, Lambert2013, Wightman2021等）
- 5枚目: ISP設定(delete, combos=0, AD→WT, V2-104M, max_ncells=300, forward_batch=64)と実行時間45分
- 6枚目: 結果ランキング図（figures/ad_smallint_isp_rank.png）を挿入
- 7枚目: 上位遺伝子表（CLDN1, MUC5B, TREM2, SORL1, MUC2, TYROBP と Shift_to_goal_end）
- 8枚目: 生物学的解釈（腸バリア+免疫が早期標的）
- 9枚目: 制限事項（統計的有意性なし, 5xFADは家族性モデル）
配色はクリーンで視認性の高いもの。テキストは箇条書きを基本とし、過剰な装飾は避ける。
検証（view ... issues）と保存（save）を必ず行うこと。
```

（脳レポート `report_ad_brain.md` も同様に、ミクログリア/BAM・APOE を強調した内容で作成）

---

## 図の参照

- `result/report/figures/ad_smallint_isp_rank.png` — 小腸 49遺伝子ランキング
- `result/report/figures/ad_brain_isp_rank.png` — 脳 23遺伝子ランキング
- `result/report/figures/ad_smallint_vs_brain_top.png` — 小腸 vs 脳 上位比較
