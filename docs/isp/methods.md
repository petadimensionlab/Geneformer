# 方法・共通設定（詳細）

このページでは、各臓器の in silico perturbation に共通する **設定** を解説します。

## 1. 基本フロー（チュートリアル準拠）

Geneformer 公式チュートリアル（`examples/in_silico_perturbation.ipynb`）に従い、以下 3 ステップで実施：

1. **状態 embedding 抽出** — `EmbExtractor.get_state_embs()` で、疾患（AD/PF）と健常（WT）の状態 embedding を計算（タイムポイントごとに 1 回、全遺伝子で共有）。
2. **仮想遺伝子削除** — `InSilicoPerturber.perturb_data()` で 1 遺伝子ずつ `delete` した変異細胞 embedding を生成。
3. **統計** — `InSilicoPerturberStats.get_stats()`（`mode="goal_state_shift"`）で `Shift_to_goal_end` を算出。

## 2. 共通パラメータ

| パラメータ | 値 | 説明 |
|---|---|---|
| `perturb_type` | `delete` | 仮想遺伝子削除 |
| `combos` | `0` | 1 遺伝子ずつ個別評価 |
| `emb_mode` | `cls` | V2 は CLS トークン使用 |
| `cell_emb_style` | `mean_pool` | |
| `model_version` | `"V2"` | Geneformer V2 |
| `state_key` | `"disease"` | 状態列 |
| `start_state` | `AD` / `PF` | 疾患状態 |
| `goal_state` | `WT` | 健常状態 |
| `alt_states` | `[]` | 通常は空 |
| `emb_layer` | `0` | |
| `summary_stat` | `"exact_mean"` | |
| `forward_batch_size` | `64`（liver は 20） | |
| `max_ncells` | `200-300`（`IS_MAX_CELLS`） | perturb 時の細胞数上限 |
| `emb_max_ncells` | `1000`（`IS_EMB_CELLS`） | 状態 embedding 細胞数 |

### 計算制約（3 台の PC で検証済み）

- **`datasets` は 4.x 必須**。`datasets>=5` だと `InSilicoPerturber.perturb_data` の `dataset.map()` がハングします。
  `uv pip install datasets==4.0.0`
- **`nproc=1` が検証済みの安定デフォルト**。`nproc>1` は datasets<5 かつ小さなデータセットでのみ安全（フォワードパスと RAM 競合）。
- **`max_ncells` は控えめに**。`perturb_data` は全変種データセットを RAM に実体化（~n_cells × variants × seq_len）。2000 細胞で MPS ホストは OOM。
- `estimate_perturb_ram()` で RAM 見積りを表示。

## 3. 1 遺伝子ずつ削除する理由

`combos=0` の状態で遺伝子リストを `genes_to_perturb` に渡すと、
`perturb_data` が「全リストの遺伝子を同時に持つ細胞」でフィルタします
（`filter_data_by_tokens`: intersection == len(tokens)）。細胞型特異的遺伝子が混ざると
この条件は空になります。正しい方法は **combos=0 で 1 遺伝子ずつ** 実行し、
各遺伝子の `Shift_to_goal_end` を個別に集めることです。

## 4. 解釈

- **`Shift_to_goal_end > 0`（正・大きい）**: この遺伝子を削除すると疾患細胞が WT に近づく
  → **早期発症の駆動因子 / 介入標的候補**。
- **`Shift_to_goal_end < 0`（負・大きい）**: 削除で WT から遠ざかる = この遺伝子は健常状態の
  維持に寄与。削除は逆効果（CD8A など基底差を反映）。
- 絶対値が小さい遺伝子は微小値に固まるため、`max_ncells` を上げて再スクリーニングが有用。

## 5. 環境変数（smoke / 制御）

| 変数 | デフォルト | 説明 |
|---|---|---|
| `IS_NPROC` | `1` | 並列数 |
| `IS_MAX_CELLS` | `200-300` | perturb 細胞数 |
| `IS_EMB_CELLS` | `1000` | embedding 細胞数 |
| `IS_MAX_GENES` | (なし) | 遺伝子数制限（smoke） |
| `IS_TIMEPOINTS` | (全) | タイムポイント制限（smoke） |
| `CELLCLASSIFIER_DIR` / `IS_CELLCLASSIFIER_DIR` | 自動 | モデルディレクトリ強制 |
| `GENEFORMER_DIR` | - | geneformer_hf の場所 |
| `ADPD_TISSUE` | 自動 | `input/<TISSUE>` |

## 6. モデル解決

ファインチューン済みセル分類器は `runs/TRAINED_MODEL_PATH.txt` を参照しますが、
**他ホストで生成された絶対パス**（例：thinkstation3 の `/home/thinkstation3/...`）を指すことがあり
ローカルでは無効です。`_isp_common.resolve_classifier_dir()` は以下の順で解決：

1. `CELLCLASSIFIER_DIR` / `IS_CELLCLASSIFIER_DIR`（env 指定）
2. `runs/TRAINED_MODEL_PATH.txt` のパスが存在する場合
3. ローカル `runs/**/ksplit*` の最新チェックポイント
4. プリトレーン済み `Geneformer-V2-104M`

これにより **どの PC で実行しても同じ出力**が得られます。