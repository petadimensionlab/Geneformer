# 骨髄（ad_bm）

早期アルツハイマー病（AD）の **骨髄** における in silico gene deletion。

## スクリプト
`analysis/07f_in_silico_perturbation_AD_BM.py`（3 台目 PC で作成、`_isp_common.py` を利用）

## 設計
- **state_key = `disease`**, **start = `AD`**, **goal = `WT`**, alt = なし
- 対象細胞型（骨髄骨髄系）: **Ly6c.high.classical.Monocytes**, **Macrophage**
- 細胞型ごとに 1 回の state embedding を共有し、遺伝子は 1 遺伝子ずつ削除

## 遺伝子リスト（11 個）
TYROBP, APOE, S100A9, S100A8, LGALS3, CCR2, CSF1R, IRF8, SPI1, CTSD, CTSB

> 各遺伝子は counts 行列で **AD 細胞の ≥40% が発現**することを確認済み
> （低発現遺伝子は `filter_data_by_tokens` でバッチ全体を空にするため）。

## 出力先（統一後）
```
input/AD_BM/results/isp/ad_bm/<celltype>/gene_<symbol>/
```

## 実行
```bash
ADPD_TISSUE=AD_BM GENEFORMER_DIR=geneformer_hf \
  .venv/bin/python analysis/07f_in_silico_perturbation_AD_BM.py
```

## 背景
骨髄骨髄系（Ly6c.high 古典単球 / マクロファージ）は AD の末梢免疫軸の先導役。
CCR2+ 単球の脳浸潤、TREM2/CSF1R ミクログリア生存、S100A8/A9 カルプロテクチン
（血中バイオマーカー）の骨髄起源に着目。

## 制御（env）
| 変数 | デフォルト | 説明 |
|---|---|---|
| `ISP_STATE_KEY` | `disease` | |
| `ISP_START_STATE` | `AD` | |
| `ISP_GOAL_STATE` | `WT` | |
| `ISP_CELLTYPES` | `Ly6c.high.classical.Monocytes,Macrophage` | 対象細胞型 |
| `IS_MAX_CELLS` | `200` | |