# officecli 用プロンプト：in silico perturbation レポートを PPT 化

このプロンプトを Claude Code / opencode に渡し、`officecli`（pptx）で
`result/report/report.md` の内容を要約したレポート用 PPT を作成してください。

## 入力資料
- 本文: `result/report/report.md`
- 図: `result/report/assets/isp_heatmap.png`, `isp_gene_ranking.png`, `isp_celltype_volcano.png`
- 表: `result/report/assets/isp_gene_ranking.csv`
- 実施一覧: `result/report/experiment_log.csv`

## 対象（officecli スキル）
`load_skill pptx` を先に読み、`officecli help pptx` でスキーマを確認してから
以下の内容で作成する。

## 作成すべき構成（10 スライド目安）
1. **タイトル**「Geneformer in silico perturbation によるパーキンソン病早期検出」
   - 副題: リンパ節 APC/単球 6 型 × 43 遺伝子、健常(WT6m)→早期PD(PFF6m)
2. **背景と目的**
   - 末梢免疫が PD 発症・早期病理に関与。APC が α-syn を感知し T 細胞応答を誘導。
   - in silico deletion で早期病理遷移を駆動し得る細胞×遺伝子を探索。
3. **対象臓器と探索方針**
   - リンパ節の Migratory.DCs / cDC.1/2 / pDCs / Macrophages / Ly6c^high Monocytes。
   - PD 核(SNCA/LRRK2/GBA1) + 免疫センサー(CD14/TLR4/TREM2/C1QB)に注目。
4. **設定とパラメータ**
   - perturb_type=delete, combos=0, goal_state_shift。
   - 6 型 × 43 遺伝子、max_ncells=200（試行50）、nproc=1, datasets=4.0.0。
5. **遺伝子リスト（表）**
   - 43 遺伝子を PD核 / 免疫センサー / 適応 / 追加 に分類した表。
6. **方法論の重要な仕様（実装ノウハウ）**
   - リスト→グループ削除になるため 1 遺伝子ずつ削除。
   - read_dictionaries が全 raw を読むため、各遺伝子を個別サブディレクトリへ。
7. **結果：ヒートマップ**（`isp_heatmap.png`）
8. **結果：コンセンサスランキング**（`isp_gene_ranking.png` + 上位表）
9. **結果：細胞型別**（`isp_celltype_volcano.png`）
10. **結論と示唆**
   - CD8A / C1QB / CD14 / TLR4 / VPS35 が「欠失＝早期PD方向」で有望。
   - CD4 / CSF1R / CD33 の欠失はむしろ保護的シフト。
   - 制約（有意性 p値なし、絶対値小）と次ステップ。
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