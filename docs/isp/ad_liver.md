# 肝臓（ad_liver）

早期アルツハイマー病（AD）の **肝臓** における in silico gene deletion。

## スクリプト
`analysis/07e_ad_liver_perturbation.py`（3 台目 PC で作成、`_isp_common.py` を利用）

## 設計
- **state_key = `state_early_ad`**（`disease` + `samples4` から構築）
- **start = `WT`**（3m の健常）、**goal = `AD3m`**（最早期の病理変化）
  - **注意**: この実験は start と goal が他の実験と逆です（他の実験は start = 疾患、goal = WT）。
    したがって `Shift_to_goal_end` の符号の読み方も他の実験と逆になります。
- **alt = `AD6m/AD9m/AD12m`**（進行方向）
- 対象細胞: Hepatocytes, Kupffer.cells, Ly6c.high Monocytes, Endothelial.cells, Macrophages

## 遺伝子リスト（21 個）
APOE, TREM2, CLU, BIN1, FAR1, IL1B, IL6, TNF, CCL2, NLRP3, CD33, ITGAX,
CD14, TLR4, TYROBP, CSF1R, C1QB, GPNMB, SPP1, MERTK, LYZ

## 出力先（統一後）
```
input/AD_liver/results/isp/ad_liver/ad_liver_early_isp_summary.csv
```

## 解釈
WT 3m から AD3m への**最早期の変化**を goal にする設計。肝臓の免疫（Kupffer/単球/
マクロファージ）と肝細胞の脂質（APOE/FAR1）・炎症軸が早期 AD の病態形成に関与するとの仮説に基づく。

## 実行
```bash
ADPD_TISSUE=AD_liver GENEFORMER_DIR=geneformer_hf \
  .venv/bin/python analysis/07e_ad_liver_perturbation.py
```
- `CELLCLASSIFIER_DIR` env でモデルを強制可能（`runs/260823_.../ksplit1` がローカル既定）
- `IS_CTYPES` で細胞型限定（例: `IS_CTYPES="Kupffer.cells"`）