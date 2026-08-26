# officecli 用プロンプト：in silico perturbation レポート（骨髄・BM）を PPT 化

このプロンプトを Claude Code / opencode に渡し、`officecli`（pptx）で
`result/report_bm/report.md` の内容を要約したレポート用 PPT を作成してください。

## 入力資料
- 本文: `result/report_bm/report.md`
- 図: `result/report_bm/assets/isp_heatmap.png`, `isp_gene_ranking.png`, `isp_celltype_volcano.png`
- 表: `result/report_bm/assets/isp_gene_ranking.csv`
- 実施一覧: `result/report_bm/experiment_log.csv`

## 対象（officecli スキル）
`load_skill pptx` を先に読み、`officecli help pptx` でスキーマを確認してから
以下の内容で作成する。

## 作成すべき構成（10 スライド目安）
1. **タイトル**「Geneformer in silico perturbation によるパーキンソン病早期検出」
   - 副題: 骨髄 APC/単球 6 型 × 42 遺伝子、健常(WT6m)→早期PD(PFF6m)
2. **背景と目的**
   - 末梢免疫（骨髄造血系）が PD 発症・早期病理に関与。単球/APC が α-syn を感知。
   - in silico deletion で早期病理遷移を駆動し得る細胞×遺伝子を探索。
3. **対象臓器と探索方針**
   - 骨髄の Ly6c.high/Ly6c.low Monocytes / Macrophage / cDCs / pDCs / MDP。
   - PD 核(SNCA/LRRK2/GBA1) + 免疫センサー(CD14/TLR4/TREM2/C1QB)に注目。
   - 骨髄特異性：単球産生の起点・MDP 前駆細胞の調整点。
4. **設定とパラメータ**
   - perturb_type=delete, combos=0, goal_state_shift（WT→PF6m）。
   - 6 型 × 42 遺伝子（発現 167 コンボ）、max_ncells=200, nproc=1, datasets=4.0.0。
5. **遺伝子リスト（表）**
   - 42 遺伝子を PD核 / 免疫センサー / 適応 / 追加 に分類した表。
6. **方法論の重要な仕様（実装ノウハウ）**
   - リスト→グループ削除になるため 1 遺伝子ずつ削除。
   - read_dictionaries が全 raw を読むため、各遺伝子を個別サブディレクトリへ。
   - 結合 CSV の列順不整合は per-gene stats 再走査で復元。
7. **結果：ヒートマップ**（`isp_heatmap.png`）
8. **結果：コンセンサスランキング**（`isp_gene_ranking.png` + 上位表）
   - CD8A / TREM2 / IL1B / CD4 / C1QB / CD14 / SNCA / TLR4 等が「欠失＝早期PD方向」。
9. **結果：細胞型別**（`isp_celltype_volcano.png`）
   - Ly6c.high で CD8A(+0.057)、Macrophage で IL1B(+0.011)、MDP で LRRK2(+0.004)。
10. **結論と示唆**
    - TREM2/CD14/TLR4/IL1B の単球炎症・スカベンジャー系が早期マーカーとして有望。
    - MDP での LRRK2 deletion シグナル（前駆レベル炎症制御）が骨髄特異的。
    - APOE/LYZ の欠失はむしろ保護的シフト。制約（p 値なし・絶対値小）に注意。
11. **実施一覧（表）** experiment_log.csv

## デザイン要件
- シンプルで視認性が高い、無地ライト背景、アクセント色は青(#2563EB)。
- 各スライドに見出し＋要点 bullet。画像は元の縦横比を保つ。
- 表は列ヘッダを強調し、数値は右寄せ。
- 最後に参考文献/補足スライドを1枚。

## 検証（納品ゲート）
- `officecli validate <file>` → エラーなし。
- `view <file> issues` → overflow / 配置問題なし。
- `view <file> screenshot --page N` で全ページを目視確認（overlap・はみ出し・
  見にくい色合いの有無を判定的に確認）し、問題があれば位置/サイズを調整して再確認。
- `save <file>` でディスクへ確実に書き出す。