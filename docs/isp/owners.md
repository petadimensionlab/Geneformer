# 執筆者名簿とページ担当

このページは、同じ Wiki を複数のエージェントが同時に編集して競合が繰り返されたため、**誰が何を書くか**を決めたものです。
用語と表記の決まりは[記述ルール](style-guide.md)、数値の基準は[評価基準](criteria.md)、これまでの競合は[競合の記録](conflicts.md)を参照してください。

## 1. 執筆者の名簿

git の著者メールで識別できます（`git log --pretty='%h %ae %s'` で確認できます）。呼び名は commit メッセージでも使います。

| 呼び名 | マシン | git の著者メール | 役割 |
|---|---|---|---|
| **AG-GB10** | このマシン（ホスト名 `spark-c083`、NVIDIA GB10） | `petadimensionlab@gmail.com` | 文書の索引、チェックリスト、競合の記録、結果ページ、方法ページ、テンプレート、レポート生成 |
| **AG-MAC** | Mac mini M4 Pro | `petadimensionlab@PetadimensionlabMacMiniM4Pro-3.local` | PD アトラス、用語集、レポート標準構成、記述ルール |
| **AG-GH** | GitHub の noreply メールを使うマシン（ホスト名は不明） | `petadimensionlab@users.noreply.github.com` | 評価基準、解析コード（`analysis/`）、パッチ（`patches/`） |
| **Nakaoka さん** | MacBook Pro | `nakaoka@MacBook-Pro-2025Masaharu-Nagayama.local` | 指示と最終判断 |

## 2. ページ担当

**編集してよいのは主担当だけ**です。他の執筆者は、直したい点を見つけたら commit メッセージに「提案」として書き、主担当に渡します。

| ページ | 主担当 | 決めた理由 |
|---|---|---|
| [記述ルール](style-guide.md) | AG-MAC | これまでの変更が最も多いため |
| [評価基準](criteria.md) | AG-GH | null 分布とブートストラップによる基準の作り直しを担当したため |
| [用語と手法の定義](glossary.md) | AG-MAC | 用語の定義を最も多く書いているため |
| [レポート標準構成](report_template.md) | AG-MAC | 節立ての提案と修正を担当したため |
| [遵守項目](checklist.md) | AG-GB10 | 新規作成と書き直しを担当したため |
| [競合の記録](conflicts.md) | AG-GB10 | 新規作成を担当したため |
| [結果](results.md) | AG-GB10 | 脾臓・肝臓・骨髄の結果の記載を担当したため |
| [方法と共通設定](methods.md) | AG-GB10 | 解釈の節の追記を担当したため |
| [再現手順](how_to_run.md) | AG-GB10 | 実行手順とレポート生成の記載を担当したため |
| [PD の現状](pd.md) | AG-GB10 | 実行状況の更新を担当したため |
| [PD 多領域アトラス](pd_atlas.md) | AG-MAC | 記載の書き直しを担当したため |
| [肝臓](ad_liver.md)、[骨髄](ad_bm.md) | AG-GB10 | 結果の記載を担当したため |
| [索引](README.md) | AG-GB10 | ページ一覧の維持を担当したため |
| `templates/`（テンプレート集） | AG-GB10 | レポート生成の雛形を担当したため |
| `analysis/`、`patches/` | AG-GH と AG-MAC | 実行したマシンの担当とする（`analysis/06_finetune.py` は AG-GH、`analysis/04b_extract_embeddings.py` は AG-MAC） |

## 3. 運用ルール

1. **編集の前に `git fetch` をします。** 他の執筆者が先に進んでいる場合が多いためです。
2. **主担当以外は、そのページを直接書き換えません。** 直したい点は commit メッセージ、または主担当への依頼として書きます。
3. **用語や基準を変えるときは、順番を守ります。** 先に[記述ルール](style-guide.md)の統一表を直し、次に[評価基準](criteria.md)の記号と各ページの参照を直し、最後にレポートを再生成します。
4. **数値には根拠の種類（A は実測、B は借用、C は未確立）を添えます。** 合否の線に使えるのは A だけです。
5. **判定の言葉は「合格、不合格、未検証、判定不能」の4語だけを使います。**
6. **push が拒否されたら、取り込んでから push し直します。** 競合が出た場合は、[競合の記録](conflicts.md)の方針に従います。
7. **commit メッセージには呼び名を書きます。** 例: `docs(isp): AG-MAC の提案を反映`。誰の変更かが後から分かるようにするためです。

## 4. 共有のしかた

このページを編集したら、[索引（README）](README.md)のページ一覧から辿れるようにしておきます。
他の執筆者は `git pull` でこのページを読み、担当を確認してください。
担当を変える場合は、変更の理由を添えてこのページを直してください。
