# PD multi-region atlas（`pd_atlas`）— ISP 実施記録

> **状態（2026-09-16 完了）**: データ取得・トークナイズ・埋め込み・ベースライン・fine-tune・
> **ISP（DMNX_Neu 16 / GPI_Neu 15 遺伝子）まで完走**。
> **ただし本実行から遺伝子を結論として引用しないでください** — §6.3 のとおり 2 領域間で
> 順位が再現せず（Spearman ρ = +0.165, p = 0.573）、§7 の留保が全て残っています。
> fine-tune は **部分 epoch（31%）** です。

`docs/isp/pd.md`（PD_spleen 単一臓器）とは別系統です。こちらは
**公開ヒト snRNA-seq アトラス（5 領域・97 ドナー）** を入力にした初の PD 解析です。

---

## 1. 入力データ（元論文）

| 項目 | 値 |
|---|---|
| 論文 | Prashant NM, Fullard JF, … Roussos P. **A multi-region single nucleus transcriptomic atlas of Parkinson's disease.** *Sci Data* **11**, 1274 (2024) |
| DOI | [10.1038/s41597-024-04117-y](https://doi.org/10.1038/s41597-024-04117-y) |
| 提供元 | CELLxGENE Discover collection `d5d0df8f-4eee-49d8-a221-a288f50a1590`（AMP PD Knowledge Platform と同一データ） |
| 規模 | **2,096,155 nuclei × 17,266 genes / 97 ドナー（75 PD / 25 control）/ 5 領域** |
| 5 領域 | DMNX（迷走神経背側核）, GPI（淡蒼球内節）, PMC（一次運動野）, PFC（前頭前野）, PVC（一次視覚野） |
| 取得物 | h5ad 30,547,659,019 bytes（**サイズ完全一致で検証**） |

**論文の設計をどう使ったか**

- 5 領域が **Braak PD 病理の皮質下 → 皮質方向の広がり**を捉える。**DMNX と GPI が早期障害**、PMC/PFC が後期、PVC はほぼ温存 → **ISP の細胞プールに DMNX_Neu / GPI_Neu を選定**
- **α-シヌクレイン（SNCA）/ Lewy 体**が神経病理の中心 → **遺伝子リストの主軸**

> ⚠️ **この 5 領域に黒質（substantia nigra）は含まれません。** したがって
> dopaminergic identity 遺伝子（TH, SLC6A3, SLC18A2, DDC）は**全細胞型で <2% の細胞にしか
> 発現していません**（実測）。存在確認スキャンで自動的に除外され、解析対象になりません。
> 「PD なのにドパミン神経を見ていない」という批判はこのデータセットの構造的制約です。

## 2. 作業サブセットとその根拠（推測ではなく実測で決定）

フルアトラス（2.1M 細胞）は 64 GiB 機では非現実的なため、**11,323 細胞**に絞りました。

| 決定 | 根拠（実測） |
|---|---|
| 11,323 細胞（`derived_class2` あたり最大 1,000、ドナー均等） | **V2-316M の MPS スループット実測**: bf16 batch 8 で 7.23 ktok/s → 50,050 細胞では埋め込みだけで 5.5 時間。実用的な規模に再設計 |
| `_1`/`_2`/`_3` の 3 ドナー群に分割 | `06_finetune.py` が要求するレプリケート分割。**同一ドナーが 2 群に入らない** |
| 分割比 | train 3,840 / eval 3,758 / test 3,725 細胞（12 細胞型） |

## 3. 実行環境とパイプライン

`macminim4pro`（Apple M4 Pro, **64 GiB** unified, MPS）/ `.venv`（Python 3.12.13,
torch 2.14.0, **transformers 4.46.3**, **datasets 4.0.0**）。
モデルは **Geneformer-V2-316M**（316.4M params, hidden 1152, 18 layers, vocab 20275, max_pos 4096）を
`./download.sh --model V2-316M` で新規取得、**bf16**（`GF_DTYPE=bf16`）で実行。

```bash
cd ~/workspace/research/Geneformer
# 入力 h5ad（CELLxGENE 30GB → 11,323 細胞のパイプライン形式）
#   input/PD_atlas/h5ad/PD_atlas.h5ad   obs: cell_id, individual, celltype, split
bash analysis/run_pd_atlas_pipeline.sh      # tokenize → embed → baseline → fine-tune → ISP
```

| 段階 | 実測 |
|---|---|
| トークナイズ | 11,323 細胞、median 2,767 tokens（max 4,096）、557 MiB |
| 凍結埋め込み抽出 | 11,323 × 1152、**1 時間 10 分**（bf16, `forward_batch_size=4`） |
| ベースライン probe | 数分（埋め込みを再利用） |
| fine-tune（600 ステップ） | `train_runtime` **2 時間 2 分**（batch 2 + gradient checkpointing） |
| ISP | DMNX_Neu 16 遺伝子 ≈ 80 分 / GPI_Neu 15 遺伝子 実行中 |

## 4. fine-tune — **部分 epoch（31%）であることの明示**

| 項目 | 値 |
|---|---|
| モデル | Geneformer-V2-316M（6/18 層を凍結） |
| `per_device_train_batch_size` | **2** |
| 1 epoch のステップ数 | 1,920（train 3,840 細胞 / batch 2） |
| **実行したステップ** | **600 = 1 epoch の 31%** |
| `train_loss` | 0.2495（0.4999 → 0.1837 → 0.0735 と低下） |
| チェックポイント | `input/PD_atlas/runs/260916_geneformer_cellClassifier_PD_atlas_celltype_Geneformer-V2-316M/ksplit1` |

**なぜ打ち切ったか**: batch 4 でスワップを使い切り（swap 12.8/14.3 GB、Metal 割当が
スワップ落ち、ステップ時間 25 秒 → 500 秒、1 epoch 換算 40 時間超）、batch 2 に落として
**パイプラインを最後まで到達させるため**ステップ上限を設定しました。フル epoch は実測
**10.37 時間**です（`FINETUNE_MAX_STEPS=0` で実行可）。

> **この値は `adpd_finetuned_summary.json` に自動記録されます**
> （`model_name` / `per_device_train_batch_size` / `max_steps_cap` / `steps_for_full_epoch` /
> `training_truncated` / `epoch_fraction_trained` / `baseline_accuracy`）。
> 部分 epoch を「1 epoch」と報告しないための仕組みです。

## 5. 結果（分類器）

test cohort（`_3`、3,725 細胞、12 クラス）:

| 手法 | accuracy | macro F1 |
|---|---:|---:|
| 凍結埋め込み + ロジスティック回帰 | 0.9812 | 0.9758 |
| **微調整 V2-316M（600/1920, ~31%）** | **0.9783** | **0.9701** |

**微調整版はベースラインをわずかに下回りました。** 上位 10 クラスは F1 ≥ 0.97 ですが、
**ISP の対象クラスが最も弱い**という重大な偏りがあります:

| クラス | precision | recall | F1 | support |
|---|---:|---:|---:|---:|
| **GPI_Neu** | 0.992 | **0.663** | 0.795 | 178 |
| **DMNX_Neu** | **0.854** | 0.978 | 0.912 | 316 |
| 他 10 クラス | — | — | ≥ 0.97 | — |

**GPI_Neu は 3 分の 1 を取りこぼし、DMNX_Neu は誤検出が多い。** ISP の状態埋め込みは
この 2 クラス上で計算されるため、**遺伝子順位はこの分類精度に依存**します。

## 6. In silico perturbation

**設計**（`analysis/07g_pd_atlas_perturbation.py`、`07_pd_spleen_early_isp.py` と同じエンジン）

```
state_key    = "disease"
start_state  = "Parkinson disease"      # ここを削除
goal_state   = "normal"                 # この方向への変位を測る
cell pools   = DMNX_Neu, GPI_Neu        # 各プール独立に実行
perturb      = "delete", combos=0, 1 遺伝子ずつ, nproc=1, max_ncells=300
```

遺伝子は実行前に**プール内 PD 細胞での検出率 ≥20%** を確認（低発現遺伝子は
combined バッチを空にしてしまうため）。結果は
`input/PD_atlas/results/isp/pd_atlas/pd_atlas_early_isp_stats_combined.csv`、
使用モデルと学習レジームは同ディレクトリの `isp_provenance.json` に記録されます。

### DMNX_Neu（16 遺伝子、完了）— `Shift_to_goal_end`

| 遺伝子 | Shift | 検出率 | | 遺伝子 | Shift | 検出率 |
|---|---:|---:|---|---|---:|---:|
| **DNAJC13** | **+0.002132** | 0.60 | | ATP13A2 | +0.000644 | 0.62 |
| **RBFOX3** | **+0.001572** | 0.42 | | NEFH | +0.000583 | 0.25 |
| **LRRK2** | **+0.001477** | 0.28 | | SNAP25 | +0.000571 | 0.34 |
| **KCNJ6** | **+0.001243** | 0.48 | | PRKN | +0.000315 | 0.37 |
| TMEM175 | +0.000916 | 0.37 | | NEFM | +0.000231 | 0.33 |
| FBXO7 | +0.000803 | 0.29 | | PINK1 | +0.000039 | 0.34 |
| SNCA | +0.000685 | 0.29 | | VPS35 | −0.000013 | 0.36 |
| — | — | — | | NEFL | −0.000112 | 0.36 |
| — | — | — | | **MAPT** | **−0.000482** | 0.57 |

解釈（`docs/isp/README.md` の定義）: 正 = 削除で PD 細胞が健常状態へ近づく
= **早期駆動因子 / 介入標的候補**。

### GPI_Neu（15 遺伝子、完了）— `Shift_to_goal_end`

| 遺伝子 | Shift | 検出率 | | 遺伝子 | Shift | 検出率 |
|---|---:|---:|---|---|---:|---:|
| **LRRK2** | **+0.000700** | 0.50 | | TMEM175 | −0.000049 | 0.38 |
| **SNCA** | **+0.000469** | 0.55 | | RBFOX3 | −0.000087 | 0.75 |
| ATP13A2 | +0.000394 | 0.51 | | MAPT | −0.000237 | 0.59 |
| PRKN | +0.000296 | 0.42 | | KCNJ6 | −0.000467 | 0.42 |
| SNAP25 | +0.000153 | 0.38 | | CALB1 | −0.000475 | 0.23 |
| NEFL | +0.000067 | 0.20 | | PINK1 | −0.000820 | 0.32 |
| DNAJC13 | +0.000021 | 0.55 | | VPS35 | −0.001236 | 0.32 |
| — | — | — | | **FBXO7** | **−0.001502** | 0.21 |

全 15 遺伝子で `|Shift| ≤ 1.5e-3`（最大は最下位 FBXO7 の負方向）で、DMNX_Neu の
最上位（DNAJC13 +2.1e-3）に届きません。**8/15 が負**で、プール全体の移動方向が
DMNX_Neu（13/16 が正）と逆です。

### 6.3 2 領域間の一致度（本実行の最も重要な結果）

共通して解析した **14 遺伝子**で両プールの順位を比較:

| 指標 | 値 | 判定 |
|---|---|---|
| 符号一致 | **8 / 14** | 偶然（7/14）と区別できない |
| Spearman ρ | **+0.165** (p = 0.573) | **有意でない** |
| Pearson r | +0.257 (p = 0.376) | **有意でない** |

順位の食い違い:

| 遺伝子 | DMNX_Neu 順位 | GPI_Neu 順位 | 備考 |
|---|---:|---:|---|
| LRRK2 | 3 | **1** | 両方で上位なのはこれのみ |
| DNAJC13 | **1** | 7 | DMNX 最上位が GPI では中位 |
| RBFOX3 | 2 | 9 | |
| KCNJ6 | 4 | 11 | |
| **FBXO7** | 6 | **14（最下位）** | **符号も逆**（+8.0e-4 → −1.5e-3） |

**結論: 本実行の遺伝子順位は領域間で再現しません。** 単一プールの上位遺伝子を
PD の知見として報告することはできないため、§7 の留保が解消するまで
**遺伝子名を結論に書かない**方針とします。

## 7. 解釈の留保（重要）

1. **2 領域間で再現しない（最も重い留保）。** 共通 14 遺伝子の順位相関は
   **Spearman ρ = +0.165 (p = 0.573)**、符号一致 8/14 で、**偶然と区別できません**（§6.3）。
   一方の領域で上位の遺伝子が他方では下位になり（DNAJC13: DMNX 1位 → GPI 7位）、
   **FBXO7 は DMNX 6位 → GPI 最下位かつ符号が反転**します。単一プールの上位遺伝子を
   知見として報告できません。
2. **絶対値が小さい。** DMNX_Neu の `Shift` は最大 2.1e-3 で、リポジトリが bf16 の
   ノイズフロアとして挙げる **|Shift| ≈ 2e-4** の 10 倍程度にすぎません
   （GPI_Neu は最大 7.0e-4 ＝ 3.5 倍）。**単独では生物学的結論を出せません。**
3. **有意性検定は未実施。** 順位を主張するには
   `analysis/07d_in_silico_perturbation_PD_null.py` 相当の**並べ替え null 比較**が必要です。
4. **分類器の偏り（§5）。** ISP の 2 プールは分類器が最も苦手なクラスであり、
   GPI_Neu の recall は 0.663 です。さらに **fine-tune は 31% epoch** です。
5. **プール間で全体の向きすら違う。** DMNX_Neu は 13/16 が正（平均 +6.6e-4）、
   GPI_Neu は 7/15 が正（平均 −1.8e-4）。「PD → 健常」方向の移動量が領域で
   一貫していません。
6. **バッチサイズは結果に影響しません。** ISP の padding は `mean_nonpadding_embs` で
   除外されるため、`forward_batch_size` は細胞埋め込みを変えません（実行時の制約のみ）。

## 8. 次にやるべきこと

ISP は完走しましたが §7 の留保が残るため、**遺伝子順位を主張する前に**次が必要です。

- [ ] **null 比較**（順位の有意性）— **最優先**。`analysis/07d_in_silico_perturbation_PD_null.py` 相当
- [ ] **フル epoch で再学習**（`FINETUNE_MAX_STEPS=0`、約 10.4 時間）。`patches/checkpoints/classifier.py`
      適用後は**中断してもチェックポイントから自動再開**できます。31% epoch の分類器のままでは
      プール埋め込み自体が信用できません
- [ ] **分類器の弱クラス対策**（クラス重み、GPI_Neu の支援サンプリング）。
      GPI_Neu の recall 0.663 は §6.3 の非再現性の一因でありうる
- [ ] 他領域への拡張（PVC を対照に PMC/PFC と比較）— 3 領域目があれば再現性の判定力が上がる

## 9. この解析で踏んだ罠（再発防止）

| 症状 | 原因と対処 |
|---|---|
| `AssertionError: Torch not compiled with CUDA enabled`（fine-tune 評価と ISP） | `patches/geneformer_multibackend.patch` が `evaluation_utils.py` / `in_silico_perturber.py` の CUDA 直書き 6 箇所を漏らしていた → **`patches/device/`** に修正を追加（学習は完走するのに指標が書かれないため発見が遅れる） |
| fine-tune が 164 秒/ステップまで劣化 | batch 4 でスワップ枯渇 → **batch 2** + gradient checkpointing |
| Metal カーネル abort | `batch × heads × seq²` が INT_MAX 超過 → **`forward_batch_size=4`**（MPS + 316M） |
| `TypeError: Got unsupported ScalarType BFloat16` | bf16 埋め出しの numpy 変換 → **`patches/bf16/`** |
| トークナイズの KeyError | `tokenize_anndata` の位置インデックス → **`patches/tokenizer.py`** |
| 部分 epoch が「1 epoch」と記録される | サマリに実レジームを記録するよう `06_finetune.py` を修正 |

> `./download.sh` は `geneformer/*.py` を取り直すため**上記パッチは全て外れます**。
> 実行後は必ず `patches/apply_patches.sh`（または README の `cp` 手順）で再適用してください。
