# Geneformer quantization docs / Geneformer 量子化ドキュメント

<a id="en"></a>
**English** · [日本語](#ja)

Everything in this folder is based on measurements on the NVIDIA GB10 host
(119.63 GiB unified memory, aarch64, torch 2.13.0+cu130, bitsandbytes 0.50.1).
Start with the consolidated experiment log, then drill into the topic documents.

| Document | Language | Content |
|---|---|---|
| [EXPERIMENT.md](EXPERIMENT.md) | EN | **Start here.** Consolidated experiment log: all measured numbers, methodology, stage-by-stage guidance, pitfall catalogue |
| [EXPERIMENT-jp.md](EXPERIMENT-jp.md) | JP | 同上（日本語版） |
| [REPORT-316M-ISP.md](REPORT-316M-ISP.md) | EN | V2-316M classifier + AD_spleen ISP canary, G5 (fp32 vs bf16) and the time/energy/load comparison |
| [REPORT-316M-ISP-jp.md](REPORT-316M-ISP-jp.md) | JP | 同上（日本語版） |
| [PLAN.md](PLAN.md) | JP | Overall quantization plan: goals, phases, acceptance gates |
| [STAGES.md](STAGES.md) | JP | Impact per pipeline stage (tokenization / embedding / ISP / fine-tuning) |
| [PREQUANT-VS-LOADTIME.md](PREQUANT-VS-LOADTIME.md) | JP | Pre-quantized weights vs quantizing at load time |
| [PLAN-316M-128GB.md](PLAN-316M-128GB.md) | JP | Plan for running V2-316M on a 128 GB machine |
| [REPORT-104M-48GB.md](REPORT-104M-48GB.md) | JP | Report: does V2-104M fit on a 48 GB GPU |

Raw evidence: `profiles/` (per-run wall time / power / energy / temperature
profiles, with timestamped logs) and `g5/` (rank-agreement JSON).

## The short version

1. **bf16 is the answer for every stage** — ~5x faster than fp32, embedding
   cosine 0.99998 (no measurable accuracy loss).
2. **int8 is safe but slower than bf16** (1.6-1.7x vs fp32); useful only where
   bf16 is unavailable or for artefact size.
3. **4-bit nf4 changes gene rankings** (top-100 overlap 87/100) — unusable where
   the output is a gene ranking, i.e. in silico perturbation.
4. **Scaling cells or genes does not make int4/int8 win** — throughput is flat
   from batch 8 to 256, and bf16 vs nf4 memory differs by 0.12 GiB at batch 256.
5. **Quantization hurts fine-tuning.** Use bf16 + gradient checkpointing
   (104M: 1.82 GiB/step) instead of QLoRA (9.83 GiB/step).
6. **Pre-quantizing for distribution is fine, for speed it is a loss** — same
   inference speed and bit-identical embeddings, but 2-6x slower loading.

## Measurement scripts (`analysis/`)

| Script | What it measures |
|---|---|
| `10_quant_bench.py` | throughput and peak memory per dtype / sequence length |
| `10b_quant_accuracy.py` | MLM loss, mask accuracy, embedding agreement |
| `10c_quant_memscale.py` | memory and throughput vs batch size |
| `10d_real_workload_mem.py` | peak memory on real tokenized cells |
| `10e_train_mem_probe.py` | one fine-tuning step (including the QLoRA comparison) |
| `11_prequant_vs_loadtime.py` | pre-quantized vs load-time quantization (speed, accuracy, load time) |
| `11b_load_time_bench.py` | repeated load-time benchmark per variant |
| `13_profile.py` | run any command under a load/energy profile (wall time, GPU %, power, Wh, CPU, RSS, °C) |
| `14_compare_isp.py` | rank agreement between two ISP runs (Spearman, top-N, sign agreement) |
| `15_isp_timing.py` | per-gene × timepoint cost from a profiled ISP log |

All scripts take the model, tissue and batch sizes from environment variables —
see each docstring.

---

## 日本語

<a id="ja"></a>
**日本語** · [English](#en)

このフォルダの内容はすべて NVIDIA GB10（統合メモリ 119.63 GiB、aarch64、
torch 2.13.0+cu130、bitsandbytes 0.50.1）での実測に基づいています。
まず実験のまとめを読み、必要に応じて各トピックのドキュメントに進んでください。

| ドキュメント | 言語 | 内容 |
|---|---|---|
| [EXPERIMENT-jp.md](EXPERIMENT-jp.md) | JP | **最初にこれ。** 実験のまとめ（全測定値・方法・段階ごとの指針・落とし穴一覧） |
| [EXPERIMENT.md](EXPERIMENT.md) | EN | 同上（英語版） |
| [REPORT-316M-ISP-jp.md](REPORT-316M-ISP-jp.md) | JP | V2-316M 分類器 + AD_spleen の ISP canary、G5（fp32 対 bf16）、時間・エネルギー・負荷の比較 |
| [REPORT-316M-ISP.md](REPORT-316M-ISP.md) | EN | 同上（英語版） |
| [PLAN.md](PLAN.md) | JP | 量子化の全体計画（目的の整理・フェーズ・基準） |
| [STAGES.md](STAGES.md) | JP | 段階ごとの影響（tokenization / embedding / ISP / fine-tuning） |
| [PREQUANT-VS-LOADTIME.md](PREQUANT-VS-LOADTIME.md) | JP | 事前量子化とロード時量子化の比較 |
| [PLAN-316M-128GB.md](PLAN-316M-128GB.md) | JP | V2-316M を 128GB のマシンで動かす計画 |
| [REPORT-104M-48GB.md](REPORT-104M-48GB.md) | JP | 104M は 48GB の GPU で動くかの検証レポート |

生の証跡: `profiles/`（実行ごとの所要時間・電力・エネルギー・温度のプロファイルと時刻付きログ）、
`g5/`（順位一致の JSON）。

### 結論だけ知りたい場合

1. **どの段階でも bf16 が正解**です — fp32 の約5倍速く、埋め込み cos 0.99998
   （精度低下は測定できないレベル）。
2. **int8 は安全ですが bf16 より遅い**です（fp32 の1.6〜1.7倍）。
   bf16 が使えない環境か、配布物のサイズ削減のときだけ使います。
3. **4bit(nf4) は遺伝子の順位を変えます**（top-100の一致 87/100）。
   遺伝子ランキングが成果物になる in silico perturbation では使えません。
4. **セル数や遺伝子数を増やしても int4/int8 は有利になりません** —
   batch 8〜256 で速さはほぼ一定、batch 256 での bf16 と nf4 のメモリ差は 0.12 GiB。
5. **微調整では量子化が逆効果**です。QLoRA（9.83 GiB/step）ではなく
   bf16 + 勾配チェックポイント（104M で 1.82 GiB/step）を使ってください。
6. **事前量子化は配布用には有効ですが、速度のためには損**です —
   推論速度と埋め込みは同一（ビット単位で一致）なのに、ロードは 2〜6倍遅くなります。

### 測定スクリプト（`analysis/`）

| スクリプト | 測るもの |
|---|---|
| `10_quant_bench.py` | 方式・系列長ごとの速さとピークメモリ |
| `10b_quant_accuracy.py` | MLM loss・マスク正解率・埋め込み一致 |
| `10c_quant_memscale.py` | batch サイズを変えたときのメモリとスループット |
| `10d_real_workload_mem.py` | 実データ（トークナイズ済みセル）でのピークメモリ |
| `10e_train_mem_probe.py` | 微調整 1 ステップ（QLoRA の比較つき） |
| `11_prequant_vs_loadtime.py` | 事前量子化とロード時量子化の比較（速度・精度・ロード時間） |
| `11b_load_time_bench.py` | 方式ごとのロード時間の反復測定 |
| `13_profile.py` | 任意のコマンドを負荷計測つきで実行（所要時間・GPU%・電力・Wh・CPU・RSS・温度） |
| `14_compare_isp.py` | ISP 2実行の順位一致（スピアマン・top-N・符号一致） |
| `15_isp_timing.py` | 計測ログから遺伝子×時点の単価を抽出 |

モデル・組織・batch はすべて環境変数で切り替えられます。各スクリプトの docstring を参照してください。
