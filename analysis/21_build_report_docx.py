#!/usr/bin/env python
"""Build the readable comparison report (DOCX) in Japanese and English.

TEMPLATE — a starting point (ひな形), not a finished analysis.
The numbers, literature verdicts and prose baked into this file are a worked
example from the AD_spleen 104M-vs-316M check (2026-09). Re-derive every number
and re-verify every PMID before reusing it on another tissue or model.

Usage: .venv/bin/python analysis/21_build_report_docx.py [--lang ja|en|both]

Outputs (docs/quantization/isp-model-comparison/):
    REPORT-ISP-104M-vs-316M.docx      (+ .pdf)  English
    REPORT-ISP-104M-vs-316M-jp.docx   (+ .pdf)  Japanese

Figures come from figures/<lang>/ (see 20_report_figures.py). Every literature
claim carries a numbered reference whose PMID/DOI was verified against Europe
PMC; the reference numbering is assigned by order of first citation in a
discarded first pass, so the printed list has no gaps.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re

import pandas as pd
from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml.ns import qn
from docx.shared import Cm, Pt, RGBColor

D = pathlib.Path("docs/quantization/isp-model-comparison")
LIT = D / "literature"
FIGROOT = D / "figures"
FONT = "Noto Sans CJK JP"

OUTNAME = {"en": "REPORT-ISP-104M-vs-316M.docx", "ja": "REPORT-ISP-104M-vs-316M-jp.docx"}

# ----------------------------------------------------------------- references
REFS: list[tuple[str, str]] = [
    ("28612290", "Haure-Mirande JV et al. (2017) Deficiency of TYROBP, an adapter protein for TREM2 and CR3 receptors, is neuroprotective in a mouse model of early Alzheimer's pathology. Acta Neuropathol. PMID 28612290"),
    ("30283031", "Audrain M et al. (2019) Deficiency of TYROBP in a tauopathy mouse model reduces C1q and normalizes clinical phenotype… Mol Psychiatry. PMID 30283031"),
    ("40931065", "ABCA7 variants impact phosphatidylcholine and mitochondria in neurons. (2025) Nature. PMID 40931065"),
    ("41125877", "Transcriptional and epigenetic targets of MEF2C in human microglia… (2025) Nat Immunol. PMID 41125877"),
    ("41620758", "Protective PLCG2 variants associate with a delayed onset of Alzheimer's disease… (2026) Alzheimers Res Ther. PMID 41620758"),
    ("41928507", "Deletion of SPI1 in microglia exacerbates amyloid pathology by impairing microglial response… (2026) Neuron. PMID 41928507"),
    ("41821106", "PU.1: the conductor of Alzheimer's symphony. (2026) Mol Neurodegener. PMID 41821106"),
    ("37949676", "Bone Marrow-Derived GCA+ Immune Cells Drive Alzheimer's Disease Progression. (2023) Adv Sci. PMID 37949676"),
    ("41152873", "Soluble CSF1R promotes microglial activation and amyloid clearance in Alzheimer's disease. (2025) J Neuroinflammation. PMID 41152873"),
    ("39555667", "Microglial Lyzl4 Facilitates β-Amyloid Clearance in Alzheimer's Disease. (2025) Adv Sci. PMID 39555667"),
    ("42676482", "Direct toll-like receptor 4 inhibition in Alzheimer's disease models: a focused scoping review. (2026) Front Aging Neurosci. PMID 42676482"),
    ("40544243", "Edaravone-Dexborneol slows down pathological progression and cognitive decline via inhibiting S100A9. (2025) Alzheimers Res Ther. PMID 40544243"),
    ("28566429", "Complement C3 deficiency protects against neurodegeneration in aged plaque-rich APP/PS1 mice. (2017) Sci Transl Med. PMID 28566429"),
    ("18562603", "Complement C3 deficiency leads to accelerated amyloid beta plaque deposition and neurodegeneration. (2008) J Neurosci. PMID 18562603"),
    ("23623698", "Griciuc A et al. (2013) Alzheimer's disease risk gene CD33 inhibits microglial uptake of amyloid beta. Neuron. PMID 23623698"),
    ("41789102", "TREM2 and microglial immunity in Alzheimer's disease: mechanisms, genetics, and therapeutic opportunities. (2026) Front Immunol. PMID 41789102"),
    ("41491073", "NLRP3 inflammasome and Alzheimer's disease: bridging inflammation and neurodegeneration. (2026) Inflammopharmacology. PMID 41491073"),
    ("36730200", "Definition of the contribution of an Osteopontin-producing CD11c+ microglial subset to Alzheimer's disease. (2023) PNAS. PMID 36730200"),
    ("42444329", "Young Adult Microglial Deletion of C1q Reduces Engulfment of Synapses… (2026) Glia. PMID 42444329"),
    ("40275327", "Apolipoprotein E in Alzheimer's disease: molecular insights and therapeutic opportunities. (2025) Mol Neurodegener. PMID 40275327"),
    ("42033879", "Adversarial validation of in silico perturbation profiles… (2026) Comput Biol Chem. PMID 42033879"),
    ("10.64898/2026.08.04.732812", "A confound-diagnostic toolkit for in silico perturbation with single-cell foundation models. (2026) preprint, doi:10.64898/2026.08.04.732812"),
    ("42321169", "Croese T et al. (2026) Disruption of the brain-spleen axis impairs monocyte-microglia communication and accelerates disease progression. Nat Commun. PMID 42321169"),
    ("41778859", "The spleen-brain axis in Alzheimer's disease and related dementias. (2026) Alzheimers Dement. PMID 41778859"),
    ("42711429", "Bone marrow myelopoiesis dysfunction in Alzheimer's disease limits monocyte homing to the brain. (2026) Nat Neurosci. PMID 42711429"),
    ("42469217", "CD33 and clusterin interact biophysically and genetically to modulate Alzheimer risk. (2026) Nat Commun. PMID 42469217"),
    ("28100745", "Jay TR et al. (2017) Disease Progression-Dependent Effects of TREM2 Deficiency in a Mouse Model of Alzheimer's Disease. J Neurosci 37:637-647. PMID 28100745"),
    ("40269681", "Benchmarking foundation cell models for post-perturbation RNA-seq prediction. (2025) BMC Genomics. PMID 40269681"),
    ("42482180", "Systematic evaluation of single-cell foundation model interpretability: attention-derived edge scores add no incremental value. (2026) BMC Genomics. PMID 42482180"),
    ("10.1101/2025.05.11.653338", "Evaluating Foundation Models for In-Silico Perturbation. (2025) bioRxiv preprint, doi:10.1101/2025.05.11.653338"),
]


def _extend_refs() -> None:
    """Add any PMID cited by the adjudication table that is not curated above,
    using the metadata Europe PMC returned (never hand-transcribed)."""
    verified = {}
    for f in ("tier1-verified-pmids.json", "subagent-verified-pmids.json"):
        fp = LIT / f
        if fp.exists():
            verified.update(json.loads(fp.read_text()))
    adjp = D / "adjudication.csv"
    if not adjp.exists():
        return
    col = pd.read_csv(adjp, dtype={"lit_pmid": str}).lit_pmid.fillna("")
    have = {k for k, _ in REFS}
    for pm in sorted({re.sub(r"\.0$", "", str(p)) for p in col if str(p).strip()} - have):
        m = verified.get(pm)
        if m:
            REFS.append((pm, f"{m['title']} ({m['year']}) {m['journal']}. PMID {pm}"))


_extend_refs()
BY_KEY = dict(REFS)
USED: set[int] = set()
MODE = "collect"
ORDER: list[str] = []
NUM: dict[str, int] = {}


def key_of(pmid) -> str:
    if pmid is None:
        return ""
    s = str(pmid).strip()
    if s in ("", "nan", "NaN"):
        return ""
    try:
        return str(int(float(s)))
    except ValueError:
        return s


def cite(*keys: str) -> str:
    ids = []
    for k in keys:
        if k not in BY_KEY:
            raise SystemExit(f"unknown reference key: {k}")
        if MODE == "collect":
            if k not in ORDER:
                ORDER.append(k)
        else:
            ids.append(f"[{NUM[k]}]")
            USED.add(NUM[k])
    return "".join(ids)


# ----------------------------------------------------------------- strings
S: dict[str, dict] = {
    "ja": {
        "kicker": "Geneformer in-silico perturbation (ISP)",
        "title": "V2-104M と V2-316M の比較、および文献との照合",
        "subtitle": "AD_spleen（脾臓免疫細胞）／55遺伝子 × 3時点／文献 {n}件（全件 PMID・DOI を機械照合）",
        "meta": "作成日 2026-09-15　リポジトリ petadimensionlab/Geneformer　データ AD_spleen（19細胞集団、AD 31,371セル・WT 31,391セル）",
        "s1": "1. 要約", "s2": "2. 背景と目的", "s3": "3. 方法", "s4": "4. 結果",
        "s5": "5. 考察", "s6": "6. 再現手順", "s7": "7. 引用文献",
        "s31": "3.1 データと計算条件", "s32": "3.2 3種類の差を分離する",
        "s33": "3.3 発現量による解釈可能性フィルタ", "s34": "3.4 文献調査と証拠の階層",
        "s41": "4.1 モデル差は測定ノイズではない", "s42": "4.2 遺伝子の3分の1で符号が逆転する",
        "s43": "4.3 候補遺伝子リストが入れ替わる", "s44": "4.4 文献との照合：316M の向きが一致する例が多い",
        "s45": "4.5 両モデルがそろって文献と食い違う遺伝子群", "s46": "4.6 シフトの大きさは発現の広さに依存する",
        "s47": "4.7 細胞集団リストの文献照合", "s51": "5.1 言えること", "s52": "5.2 限界", "s53": "5.3 推奨",
        "sum_intro": "知りたかったのは「モデルを 104M から 316M に替えると、遺伝子欠失の結果は質的に変わるのか」です。"
                     "結論は「変わる。しかも測定ノイズではなく、遺伝子の3分の1で向き（正負）が逆転する」でした。"
                     "どちらが正しいかは、片方のモデルでは決められません。文献と突き合わせると 316M の向きが一致する例が多く、"
                     "ただし両モデルがそろって文献と食い違う遺伝子群も残ります。",
        "sum_b1": "① モデル差は本物です。同じ細胞・同じ遺伝子・同じ時点で測った 104M と 316M の差は、シグナルの {model:.1f}%"
                  "（平均 |差| ÷ 平均 |シフト|）。これは bf16 化の差（{prec:.1f}%）の約7倍、"
                  "細胞サンプリングの差（{samp:.1f}%）の約5倍です。",
        "sum_b2": "② 55遺伝子のうち {flips} 個でシフトの符号が逆転しました（{pct:.0f}%）。順位の相関（スピアマン ρ）は {rho:+.3f} しかありません。",
        "sum_b3": "③ 発現量で絞っても消えません。標的細胞での検出率が 3% 以上の 44 遺伝子に限っても 13 個が逆転します。"
                  "つまり「あまり発現していない遺伝子のノイズ」ではありません。",
        "sum_b4": "④ 文献が方向を示せる遺伝子は {ndec} 個。うち両モデルが向きを争っているのは {ndis} 個で、"
                  "316M が文献に一致 {n316} 個、104M が一致 {n104} 個でした（LYZ の 1 個は証拠が弱く、LYZ 本体の論文ではありません）。",
        "sum_b5": "⑤ 逆に、両モデルがそろって文献と食い違う群が {ncontra} 個あります（TNF・IL1B・C1QA/B/C・S100A8/A9）。"
                  "分泌性の炎症遺伝子に偏っており、この手法の系統的な弱点を示しています。",
        "sum_b6": "⑥ 実務的な推奨は「両モデルで同じ向きが出た遺伝子だけを報告する」。片方だけの結果は、"
                  "モデル選択で入れ替わるので候補として持ち出さないことです。マーカー遺伝子（CD8A 等）は除外します。",
        "bg1": "in-silico perturbation（ISP、計算機上の遺伝子操作）は、細胞の遺伝子発現プロファイルから"
               "「この遺伝子を消したら細胞の状態がどう動くか」を予測する Geneformer の機能です。"
               "本研究では AD（アルツハイマー病）モデルマウスの脾臓から取った免疫細胞に、"
               "55 個の仮説遺伝子それぞれを 1 つずつ欠失させ、その細胞状態が健常（WT）の状態にどれだけ近づくかを測っています。",
        "bg2": "指標は Shift_to_goal_end です。ある遺伝子を欠失させたときの埋め込みベクトルと WT 状態の"
               "コサイン類似度から、元の状態と WT のコサイン類似度を引いた値です。"
               "正の値は「その遺伝子を消すと WT 状態に近づく」＝その遺伝子は病気を進める側（治療標的候補）、"
               "負の値は「消すと WT から遠ざかる」＝その遺伝子は健常状態を支えている、と読みます。"
               "この符号の規約が本報告のすべての議論の前提です。",
        "bg3": "ところが我々は、モデルを 104M から 316M（パラメータ約3倍）に替えると"
               "同じデータでも遺伝子の順位が大きく変わることを観測しました。本報告の目的は、"
               "(a) その差が測定ノイズなのか本物なのかを定量し、(b) 対象とした 55 遺伝子と 19 細胞集団について文献を調べ、"
               "(c) どちらのモデルの向きが公表された生物学と整合するかを判定することです。",
        "m_t_h": ["項目", "内容"],
        "m_rows": [
            ["組織・細胞", "AD_spleen（脾臓）19 集団の免疫細胞（単球・マクロファージ・好中球・DC 4種・T 4種・NKT・NK/ILC1・B 4種・形質細胞ほか）"],
            ["細胞数", "AD 31,371 / WT 31,391（プール内）。上流で細胞種ごとに最大 3,000 セルに間引き"],
            ["時点", "3m・4p5m・6m（3・4.5・6 か月齢）"],
            ["遺伝子", "仮説 55 遺伝子（AD リスク遺伝子、自然免疫、補体、適応免疫、系統マーカー）"],
            ["比較した 2 モデル", "V2-104M（fp32、200セル）と V2-316M（bf16、200セル）"],
            ["指標", "Shift_to_goal_end（3m 時点の Shift_3m で順位付け）。正 = WT 状態に近づく"],
            ["実行環境", "NVIDIA GB10（aarch64, sm_121）、torch 2.13.0+cu130、transformers 4.46.3"],
        ],
        "m32": "2 つの ISP 実行の差には、少なくとも 3 つの原因が混ざります。(1) 数値精度（fp32 と bf16）、"
               "(2) 細胞のサンプリング（200 セルと 100 セル）、(3) モデルそのもの（104M と 316M）。"
               "これを切り分けるため、同一 checkpoint・同一細胞・同一遺伝子で dtype だけ変えたペアと、"
               "同一モデルで細胞数だけ変えたペアを用意し、それぞれの差を測りました。",
        "m33": "標的細胞でほとんど発現していない遺伝子は、欠失させても動かせる細胞がごく僅かで、"
               "シフト値は意味を持ちません。トークン化データから遺伝子ごとの検出率"
               "（その遺伝子を持つ細胞の割合）を直接計算し、検出率で結果を層別しました。",
        "m34a": "文献は Europe PMC の API で検索し、要旨を取得しました。検索は遺伝子ごとに"
                "「その遺伝子を欠失・阻害すると病態が改善するか悪化するか」を狙って組み立てています。"
                "得られた PMID は 1 件ずつ API で引き直し、題名・雑誌・年が一致することを機械的に確認しました"
                "（本報告で引用する文献はすべて確認済み）。",
        "m34b": "証拠は 2 層に分けています。第 1 層は本調査で直接取得して要旨を読んだもの。第 2 層は並列に走らせた"
                "文献調査の結果で、PMID を機械照合したもの。第 1 層が沈黙している遺伝子に限って第 2 層を使い、"
                "どちらも無い場合は「証拠不十分」として推測で埋めていません。",
        "m34c": "第 2 層の引用は、PMID の存在確認に加えて「主張された遺伝子と取得した題名が対応するか」も確認しました。"
                "対応しなかった例（C1QC と IL1B）は、C1q は複合体レベルの知見に置き換え、IL-1 は引用せず「保留」としています。"
                "推測で埋めないことを優先しました。",
        "r_th": ["比較", "遺伝子数", "ρ（順位相関）", "符号の逆転", "平均|差| ÷ 平均|シフト|"],
        "r41a": "モデル差の平均 |差| はシグナルの {model:.1f}% で、順位相関は {rho:+.3f}。一方、精度差は {prec:.1f}%"
                "（ρ={prho:+.3f}）、サンプリング差は {samp:.1f}%（ρ={srho:+.3f}）です。"
                "つまり 104M と 316M の不一致は、測定誤差の範囲では説明できません。",
        "r41b": "副次的に分かったこととして、細胞サンプリングの差（17〜18%）は bf16 化の差（12%）より大きいです。"
                "bf16 を使うことに起因する不確かさは、細胞の選び方による不確かさより小さい、と言えます"
                "（bf16 化の判断自体は別報告のとおり妥当でした）。",
        "fig1c": "図1　差の正体を分離した結果。モデル差（赤）はサンプリング差（橙）の約5倍、精度差（緑）の約7倍。",
        "fig2c": "図2　同じ条件で測った 104M と 316M のシフト。赤は符号が逆転した遺伝子。右は拡大図（全体の大半はここに入る）。",
        "r42a": "{flips} 個の遺伝子で 104M と 316M のシフトの符号が逆になりました。内訳は「104M が正・316M が負」が {p2n} 個、"
                "「104M が負・316M が正」が {n2p} 個で、どちらか一方に偏っていません。"
                "つまり「316M の方が保守的」といった系統的な傾向ではなく、遺伝子ごとにどちらが正しいかを争っている状態です。",
        "r42th": ["遺伝子", "検出率", "104M の Shift", "順位", "316M の Shift", "順位", "文献"],
        "r42note": "（順位は 55 遺伝子中の順位で、1 が最も「WT に近づく」側。文献欄は後述 4.4 の判定。）",
        "r43a": "実務上いちばん困るのは、標的候補として挙げる遺伝子リストが入れ替わることです。"
                "標的細胞での検出率が 3% 以上で、かつシフトが正（＝消すと WT に近づく）という基準で拾うと、次のようになります。",
        "r43th": ["遺伝子", "104M の判定", "316M の判定", "検出率", "文献の向き", "Shift 104M / 316M"],
        "r43cand104": "標的候補（104M のみ）", "r43cand316": "標的候補（316M のみ）", "r43na": "—",
        "r43b": "104M だけが「候補」とする遺伝子には ABCA7・MEF2C・PLCG2・SPI1 という AD の代表的リスク遺伝子が含まれます。"
                "逆に 316M だけが「候補」とするのは TYROBP・GCA・CSF1R・LYZ・TLR4 です。",
        "r44a": "文献が「欠失で改善」か「悪化」かを明確に決められる遺伝子は {ndec} 個ありました。"
                "このうち両モデルが向きを争っている {ndis} 個では、316M が文献に一致したのが {n316} 個、"
                "104M が一致したのが {n104} 個です。",
        "r44th": ["遺伝子", "104M", "316M", "文献の向き（根拠）", "文献と一致"],
        "r44sens": "判定の感度について。7 個のうち 2 個は、文献自身が条件に依存しています。TREM2 は疾患段階に依存し、"
                   "早期のアミロイドでは欠損が病態を軽減し、後期では悪化させます{cite}。本データは 3〜6 か月齢の早期なので、"
                   "早期の読み方を採れば TREM2 の符号は反転し得ます。TLR4 は指標に依存し、炎症指標では阻害が有益とする報告が"
                   "優勢ですが、アミロイド処理では不一致です。仮に TREM2 を判定不能、TLR4 を改善方向と数えると、"
                   "316M が 6 個・104M が 1 個（母数 7）になります。個々の判定は揺れますが、"
                   "「316M の向きの方が文献に合う」という結論の向きは変わりません。",
        "fig3c": "図3　2モデルで向きが逆転し、文献が方向を示している 7 遺伝子。上段は 6 遺伝子（横軸 ±0.0015）、"
                 "下段は突出する LYZ を別スケールで示す。",
        "bullets44": [
            "ABCA7：機能喪失変異が AD リスクを大きく上げる＝保護的遺伝子。欠失は悪化方向と予測される{cite:40931065}。316M は負（一致）、104M は正（不一致）。",
            "MEF2C：マイクログリアの炎症を抑える因子で、低下が疾患に寄与する{cite:41125877}。316M が一致。",
            "PLCG2：保護的な P522R は「軽度亢進型（hypermorph）」で、酵素活性が上がる方向が保護的{cite:41620758}。欠失は逆方向なので悪化と予測。316M が一致。",
            "TREM2：機能喪失変異は免疫監視を破綻させアミロイド病理を悪化させる{cite:41789102}。316M が一致（ただし TREM2 の検出率は 0.6% と低く、断定はできません）。",
            "TYROBP：欠損が APP/PSEN1 マウスで神経保護的{cite:28612290}、タウオパチーモデルでも臨床表現型が正常化する{cite:30283031}。欠失は改善方向＝正。316M が一致。TYROBP は検出率 43% と本データで最も広く発現しており、この一致は信頼度が高いです。",
            "GCA（グランカルシン）：骨髄由来の GCA 陽性免疫細胞が AD 進行を駆動し、組み換え GCA の投与でアミロイド斑と認知が悪化する{cite:37949676}。末梢免疫細胞を扱った知見で、脾臓プールに直結します。316M が一致。",
            "LYZ（リゾチーム）：リゾチーム様酵素 Lyzl4 は Aβ 除去を促進する{cite:39555667}ため欠失は悪化方向と予測され、これだけ 104M が一致します。ただし Lyzl4 は LYZ そのものではなく、証拠としては弱いものです。104M の値（−0.0233）は突出して大きく、CD8A に次ぐ外れ値で、系統マーカー削除のアーティファクトを疑う余地があります。",
        ],
        "r45a": "{n} 個の遺伝子では、どちらのモデルも文献と逆を向きました。内容を見ると TNF・IL1B・C1QA/B/C・S100A8/A9 と、"
                "分泌性の炎症メディエーターと補体に偏っています。",
        "r45th": ["遺伝子", "検出率", "104M", "316M", "文献の向き（根拠）"],
        "r45b": "これらは「欠失させると病態が改善する」という文献（TNF 欠失で Aβ 産生が減る、C1q 阻害でシナプス貪食が減る、"
                "S100A9 阻害でプラークが減る）に対し、ISP ではどちらのモデルも負（＝WT から遠ざかる）を返しています。"
                "分泌性の遺伝子を消すと「活性化状態を失った」ように埋め込みが動き、「病気から遠ざかる」として現れる可能性があります。"
                "この手法をそのまま治療標的探索に使う際の系統的な限界と考えられます。",
        "fig4c": "図4　検出率（横軸）とシフトの大きさ（縦軸、対数）。発現の広い遺伝子ほどシフトが大きい。赤は符号が逆転した遺伝子。",
        "r46a": "シフトの大きさは検出率と強く連動します。CD8A（検出 9%）は 104M で −0.098 と他を圧倒しますが、"
                "CD8A は細胞の系統を決めるマーカーであり、「CD8A を消すと細胞の identity が崩れる」という効果が、"
                "「WT 状態から遠ざかる」として現れていると解釈するのが自然です。"
                "in-silico 欠失の敵対的検証では、発現量マッチのヌルと比べて有意な効果がほとんど残らないという報告もあります{cite:42033879,10.64898/2026.08.04.732812}。"
                "マーカー遺伝子の削除結果は候補リストから除外することを勧めます。",
        "r46b": "手法そのものへの批判も文献にあります。単細胞基盤モデルの摂動応答予測を正面から評価した研究では、"
                "単純な平均・線形のベースラインが基盤モデルを上回ること{cite:40269681}、"
                "注意由来のペアワイズ情報が追加の価値を持たないこと{cite:42482180}、"
                "ゼロショットの Geneformer が評価対象の摂動の半数程度でしか状態分離の基準を満たさないこと{cite:10.1101/2025.05.11.653338}が"
                "報告されています。本報告の慎重な結論（両モデル一致のみを候補とする）は、この批判を前提にしています。",
        "r47a": "対象とした 19 集団について、AD での変化を文献で確認しました。研究間の不一致が多く、確度は集団によって大きく異なります。",
        "r47th": ["細胞集団", "AD の細胞数", "WT の細胞数", "文献での報告", "確度"],
        "r47b": "重要な制約があります。このデータセットは上流の前処理で細胞種ごとに最大 3,000 セルに間引かれており、"
                "CD4 T・CD8 T・ナイーブ B・辺縁帯 B は AD/WT ともちょうど 3,000 セルです。したがって「AD で B 細胞が増える」"
                "といった文献の主張を、このデータの細胞数比で検証することは原理的にできません。文献の変化方向は、"
                "ISP の結果を解釈する背景としてのみ使い、比率の一致・不一致は判定に含めていません（図5）。",
        "fig5c": "図5　プール内の細胞構成。AD と WT でほぼ同一（間引きのため細胞数比の比較は不可）。",
        "r47c": "脾臓と脳をつなぐ軸そのものは、AD の病態に関与することが報告されています。脾臓神経の遮断は単球–マイクログリア間の"
                "情報伝達を壊し病態を加速させ{cite:42321169}、脾臓–脳軸は免疫と代謝の制御を通じて疾患に関わります{cite:41778859}。"
                "また AD では骨髄の造血が乱れ、単球の中枢への移動が損なわれます{cite:42711429}。"
                "本解析の対象（脾臓の免疫細胞）は、この軸の末梢側に位置づけられます。",
        "d51": [
            "モデル選択は結論を変えます。104M と 316M の不一致は測定ノイズの 5〜7 倍で、55 遺伝子中 18 個（33%）で符号が逆転します。片方のモデルだけで標的候補を出すのは危険です。",
            "文献が方向を決められる遺伝子に限ると、316M の向きが一致する例が多い（7 対 1）。特に TYROBP（検出 43%）と GCA（骨髄由来免疫細胞が AD を駆動）は、末梢免疫という文脈に直結する知見で 316M が支持されます。",
            "一方で、これは「316M を採用せよ」という結論にはなりません。判定できた遺伝子が 7 個しかなく、しかも逆の例（LYZ、あるいはヒト遺伝学とマウス欠失実験が対立する SPI1）が残ります。",
            "補体・炎症性メディエーター（C1q、TNF、IL1B、S100A8/A9）では、両モデルがそろって文献と食い違います。これはモデル選択の問題ではなく、手法（単一遺伝子のトークン削除 → 埋め込み距離）の限界である可能性が高いです。",
        ],
        "d52": [
            "目標状態の意味：指標が測っているのは「脾臓免疫細胞の WT 状態への近さ」であり、「AD 病理の軽さ」ではありません。APOE が良い例で、APOE は AD 最大のリスク遺伝子ですが、両モデルとも負（WT から遠ざかる）を返しました。脳アミロイドの知見と脾臓免疫の状態は、同じ方向を向くとは限りません。",
            "発現交絡：シフトの大きさは検出率に強く依存します（図4）。発現量マッチのヌルに対する検証が必須で、本解析では未実施です。",
            "マーカー遺伝子：CD8A・CD4・CD19 などの系統マーカーは、削除すると identity が崩れる方向に埋め込みが動くため、治療標的として解釈できません。",
            "検出力：文献が方向を示せる遺伝子は 55 個中 17 個、そのうちモデル間で争いがあり発現も十分なのは 7 個だけです。7 対 1 という数字はその 7 個の中での話で、統計的な一般化はできません。",
            "両モデル一致は独立な再現ではありません：104M と 316M はトークン化と rank エンコードを共有しているため、発現量に由来する効果は両方で一致して現れます。したがって「両モデルで一致」は必要条件ではあっても、十分条件ではありません。独立な検証には CRISPRi/Perturb-seq など別モダリティが必要です。",
            "比較の交絡：104M と 316M の比較には dtype（fp32 と bf16）も同時に変わっています。ただし精度差はサンプリング差より小さく、順位への影響は限定的と判断しました。",
        ],
        "d53": [
            "両モデルで同じ向きのシフトが出た遺伝子だけを「候補」として報告する。",
            "片方のモデルだけで得た候補は、モデル選択で入れ替わるので単独では提示しない。",
            "系統マーカー（CD8A、CD4、CD19、NKG7 など）は標的候補から除外し、集団レベルの効果として別に扱う。",
            "標的細胞での検出率（例：5% 未満）を明記し、低発現遺伝子の順位は解釈しない。",
            "補体・分泌性メディエーターの結果は、この手法の弱点として扱い、文献の因果知見（阻害実験）を優先する。",
            "指標の意味を常に添える：測っているのは WT 免疫状態への近さで、病理の改善ではない。",
        ],
        "s6a": "本報告の数値・図・表は、すべてリポジトリ内のスクリプトで再生成できます。",
        "s6th": ["スクリプト", "内容"],
        "s6rows": [
            ["analysis/16_gene_expression_by_celltype.py", "細胞種・病態別の遺伝子検出率（発現の広さ）"],
            ["analysis/17_isp_model_compare.py", "3種類の差の分離、符号反転、候補リストの比較"],
            ["analysis/18_literature_check.py", "Europe PMC からの文献取得（遺伝子ごとの検索式を内蔵）"],
            ["analysis/19_adjudicate.py", "文献判定の統合と PMID の機械照合"],
            ["analysis/20_report_figures.py", "図1〜図5の生成（日英）"],
            ["analysis/21_build_report_docx.py", "本報告（DOCX、日英）の生成"],
        ],
        "s6b": "生データは docs/quantization/isp-model-comparison/ 以下（gene_table.csv、summary.json、adjudication.csv、literature/、figures/）にあります。",
        "s7a": "PMID・DOI は Europe PMC API で 1 件ずつ照合し、題名・雑誌・年が一致することを確認済みです。被引用の無い文献は掲載していません。",
        "pos": "欠失で改善", "neg": "欠失で悪化", "mixed": "相反", "none": "—",
        "mixed_long": "相反（判定不能）", "none_long": "因果的証拠なし",
        "figref": "図",
        "lit_hdr": "文献の向き",
    },
    "en": {
        "kicker": "Geneformer in-silico perturbation (ISP)",
        "title": "V2-104M versus V2-316M: comparison and literature check",
        "subtitle": "AD_spleen (splenic immune cells) / 55 genes x 3 timepoints / {n} references (all PMID/DOI machine-verified)",
        "meta": "Prepared 2026-09-15  Repo petadimensionlab/Geneformer  Data AD_spleen (19 cell populations; 31,371 AD and 31,391 WT cells)",
        "s1": "1. Summary", "s2": "2. Background and objective", "s3": "3. Methods", "s4": "4. Results",
        "s5": "5. Discussion", "s6": "6. How to reproduce", "s7": "7. References",
        "s31": "3.1 Data and compute conditions", "s32": "3.2 Separating the three sources of difference",
        "s33": "3.3 Interpretability filter on expression", "s34": "3.4 Literature search and evidence tiers",
        "s41": "4.1 The model gap is not measurement noise", "s42": "4.2 One gene in three flips sign",
        "s43": "4.3 The candidate list is not stable between models", "s44": "4.4 Literature check: 316M agrees in most contested cases",
        "s45": "4.5 Genes where both models contradict the literature", "s46": "4.6 Shift size tracks expression breadth",
        "s47": "4.7 Literature check on the cell-population list", "s51": "5.1 What we can say", "s52": "5.2 Limitations",
        "s53": "5.3 Recommendations",
        "sum_intro": "The question was whether swapping the model from 104M to 316M changes the outcome of gene-deletion "
                     "perturbation qualitatively. The answer is yes, and not through measurement noise: in one gene in three "
                     "the direction of the shift flips sign. A single model cannot decide which direction is right. Checked "
                     "against the literature, 316M's direction agrees in most of the contested cases, but a group of genes "
                     "remains where both models contradict the published evidence.",
        "sum_b1": "(1) The model gap is real. On the same cells, genes and timepoint the difference between 104M and 316M is "
                  "{model:.1f}% of the signal (mean |difference| divided by mean |shift|). That is about 7x the bf16 rounding "
                  "difference ({prec:.1f}%) and about 5x the cell-sampling difference ({samp:.1f}%).",
        "sum_b2": "(2) In {flips} of the 55 genes the sign of the shift reverses ({pct:.0f}%). The rank correlation (Spearman rho) is only {rho:+.3f}.",
        "sum_b3": "(3) Filtering by expression does not remove it. Restricting to the 44 genes detected in at least 3% of the "
                  "target cells still leaves 13 sign reversals, so this is not noise from barely expressed genes.",
        "sum_b4": "(4) The literature fixes a direction for {ndec} genes. In {ndis} of them the two models disagree, and 316M "
                  "agrees with the literature in {n316} cases against {n104} for 104M (the single 104M case, LYZ, rests on weak "
                  "evidence about a related gene, not LYZ itself).",
        "sum_b5": "(5) Conversely, in {ncontra} genes both models contradict the literature (TNF, IL1B, C1QA/B/C, S100A8/A9). "
                  "They are concentrated among secreted inflammatory mediators and complement, which points to a systematic "
                  "weakness of the assay rather than a model-choice problem.",
        "sum_b6": "(6) The operational recommendation is to report only genes whose direction agrees across both models. A "
                  "candidate found with one model alone moves with the model choice, so it should not be presented on its own. "
                  "Lineage markers (CD8A and similar) must be excluded.",
        "bg1": "In-silico perturbation (ISP) is Geneformer's function for predicting how a cell's state moves when a gene is "
               "removed, from the cell's own gene-expression profile. Here we delete each of 55 candidate genes, one at a time, "
               "in splenic immune cells taken from a mouse model of Alzheimer's disease (AD), and measure how much closer the "
               "cell state moves to the healthy (WT) state.",
        "bg2": "The metric is Shift_to_goal_end: the cosine similarity between the perturbed embedding and the WT state, minus "
               "the cosine similarity between the original embedding and the WT state. A positive value means deleting the gene "
               "moves the cell toward the WT state, i.e. the gene pushes the disease forward and is a candidate therapeutic "
               "target; a negative value means deleting it moves the cell away from WT, i.e. the gene supports the healthy "
               "state. This sign convention underpins every statement in this report.",
        "bg3": "We observed that moving from 104M to 316M (about 3x the parameters) substantially reorders the genes on the "
               "same data. This report therefore (a) quantifies whether that difference is measurement noise or real, "
               "(b) examines the literature for the 55 target genes and the 19 target cell populations, and (c) judges which "
               "model's direction is consistent with the published biology.",
        "m_t_h": ["Item", "Detail"],
        "m_rows": [
            ["Tissue / cells", "AD_spleen, 19 splenic immune populations (monocytes, macrophages, neutrophils, 4 DC types, 4 T-cell types, NKT, NK/ILC1, 4 B-cell types, plasma cells)"],
            ["Cell counts", "31,371 AD / 31,391 WT in the pool; upstream preprocessing caps each cell type at 3,000 cells"],
            ["Timepoints", "3m, 4p5m, 6m (3, 4.5 and 6 months of age)"],
            ["Genes", "55 hypothesis genes (AD risk genes, innate immunity, complement, adaptive immunity, lineage markers)"],
            ["Models compared", "V2-104M (fp32, 200 cells) and V2-316M (bf16, 200 cells)"],
            ["Metric", "Shift_to_goal_end, ranked by the 3m timepoint (Shift_3m); positive = toward the WT state"],
            ["Environment", "NVIDIA GB10 (aarch64, sm_121), torch 2.13.0+cu130, transformers 4.46.3"],
        ],
        "m32": "Any difference between two ISP runs mixes at least three sources: (1) numerical precision (fp32 vs bf16), "
               "(2) which cells entered the run (200 vs 100), and (3) the model itself (104M vs 316M). To separate them we "
               "measured a pair that differs only in dtype on identical cells and genes, and a pair that differs only in the "
               "number of cells within one model.",
        "m33": "A gene that is barely expressed in the target cells can only move a handful of cells when deleted, so its shift "
               "value carries no meaning. We computed the detection rate of every gene directly from the tokenised data (the "
               "fraction of cells containing it) and stratified all results by that rate.",
        "m34a": "The literature was searched through the Europe PMC API and abstracts were retrieved. Each query was written to "
                "target the perturbation question, i.e. whether deleting or inhibiting the gene improves or worsens the "
                "pathology. Every PMID used was then re-fetched one by one and its title, journal and year checked "
                "mechanically; all references cited here passed that check.",
        "m34b": "Evidence is kept in two tiers. Tier 1 is what this study retrieved and read directly. Tier 2 comes from parallel "
                "delegated literature sweeps with their PMIDs re-verified. Tier 2 is used only where tier 1 is silent, and where "
                "neither exists the verdict is reported as insufficient evidence rather than guessed.",
        "m34c": "For tier-2 citations we also checked that the retrieved title matches the claimed subject, not just that the "
                "PMID exists. Two such citations failed (C1QC and IL1B): the C1QC claim was replaced with complex-level C1q "
                "evidence, and the IL-1 claim is left unresolved instead of being cited.",
        "r_th": ["Comparison", "Genes", "rho (rank corr.)", "Sign flips", "mean|diff| / mean|shift|"],
        "r41a": "For the model gap the mean absolute difference is {model:.1f}% of the signal and the rank correlation is "
                "{rho:+.3f}. The precision difference is {prec:.1f}% (rho={prho:+.3f}) and the cell-sampling difference is "
                "{samp:.1f}% (rho={srho:+.3f}). The disagreement between 104M and 316M therefore cannot be explained as "
                "measurement error.",
        "r41b": "One secondary finding: the cell-sampling difference (17-18%) is larger than the bf16 difference (12%). The "
                "uncertainty introduced by using bf16 is smaller than the uncertainty introduced by which cells are selected, "
                "which supports the earlier conclusion that bf16 is the right default.",
        "fig1c": "Figure 1  Separating the sources of the gap. The model difference (red) is about 5x the sampling difference "
                 "(orange) and about 7x the precision difference (green).",
        "fig2c": "Figure 2  Shifts measured under identical conditions in 104M and 316M. Red points flip sign between models; "
                 "the right panel zooms in on the bulk of the data.",
        "r42a": "In {flips} genes the sign of the shift reverses between 104M and 316M: {p2n} cases go from positive (104M) to "
                "negative (316M), and {n2p} go the other way. The reversals are not biased in either direction, so this is not "
                "a systematic tendency such as 316M being more conservative; the two models disagree per gene.",
        "r42th": ["Gene", "Detection", "Shift 104M", "Rank", "Shift 316M", "Rank", "Literature"],
        "r42note": "(Rank is out of 55 genes, 1 being closest to the WT state. The literature column reports the verdict of section 4.4.)",
        "r43a": "The practical problem is that the list of candidate targets changes. Taking genes detected in at least 3% of "
                "the target cells with a positive shift (deletion moves the cells toward WT) gives the following.",
        "r43th": ["Gene", "Verdict in 104M", "Verdict in 316M", "Detection", "Literature direction", "Shift 104M / 316M"],
        "r43cand104": "candidate (104M only)", "r43cand316": "candidate (316M only)", "r43na": "-",
        "r43b": "The genes that only 104M calls candidates include ABCA7, MEF2C, PLCG2 and SPI1, all canonical AD risk genes. "
                "The genes only 316M calls candidates are TYROBP, GCA, CSF1R, LYZ and TLR4.",
        "r44a": "The literature fixes a direction for {ndec} genes. In the {ndis} where the two models disagree, 316M matches "
                "the literature in {n316} cases and 104M in {n104}.",
        "r44th": ["Gene", "104M", "316M", "Literature direction (source)", "Agrees"],
        "r44sens": "Sensitivity of these verdicts. Two of the seven are conditional in the literature itself. TREM2 depends on "
                   "disease stage: deficiency reduces pathology early in amyloid deposition and worsens it later{cite}. Our data "
                   "are early (3-6 months), so an early-stage reading could invert the TREM2 sign. TLR4 depends on the endpoint: "
                   "inflammation endpoints favour benefit from inhibition, whereas amyloid handling is inconsistent. Treating "
                   "TREM2 as unresolved and TLR4 as beneficial gives 6 for 316M and 1 for 104M out of 7. Individual verdicts move, "
                   "but the direction of the conclusion, that 316M's sign matches the literature more often, does not.",
        "fig3c": "Figure 3  The seven genes where the models disagree in direction and the literature indicates one. The upper "
                 "panel holds six genes (x-axis plus/minus 0.0015); LYZ is shown separately because it dominates in 104M.",
        "bullets44": [
            "ABCA7: loss-of-function variants substantially increase AD risk, so the gene is protective and deletion is predicted to worsen the state{cite:40931065}. 316M is negative (agrees); 104M is positive (does not).",
            "MEF2C: MEF2C restrains microglial inflammation and its loss contributes to disease{cite:41125877}. 316M agrees.",
            "PLCG2: the protective P522R allele is a mildly hyperactive (hypermorph) form, so higher enzyme activity is protective{cite:41620758}; deletion acts in the opposite direction and should worsen the state. 316M agrees.",
            "TREM2: loss-of-function variants break immune surveillance and worsen amyloid pathology{cite:41789102}. 316M agrees, but TREM2 is detected in only 0.6% of the target cells, so this cannot be conclusive.",
            "TYROBP: deficiency is neuroprotective in APP/PSEN1 mice{cite:28612290} and normalises the clinical phenotype in a tauopathy model{cite:30283031}, so deletion is predicted to help (positive). 316M agrees. TYROBP is detected in 43% of the target cells, the broadest expression in this panel, so this agreement carries the most weight.",
            "GCA (grancalcin): bone-marrow-derived GCA-positive immune cells drive AD progression, and recombinant GCA worsens plaques and cognition{cite:37949676}. This is peripheral-immune biology and maps directly onto a splenic pool. 316M agrees.",
            "LYZ (lysozyme): the lysozyme-like enzyme Lyzl4 promotes amyloid-beta clearance{cite:39555667}, so deletion is predicted to be harmful and 104M is the model that agrees. The evidence concerns Lyzl4 rather than LYZ itself and is therefore weak, and the 104M value (-0.0233) is an outlier second only to CD8A, which raises the possibility of a lineage-artefact.",
        ],
        "r45a": "In {n} genes both models point opposite to the literature. They are TNF, IL1B, C1QA/B/C and S100A8/A9, i.e. "
                "concentrated among secreted inflammatory mediators and complement.",
        "r45th": ["Gene", "Detection", "104M", "316M", "Literature direction (source)"],
        "r45b": "The literature says that deleting these genes improves the pathology (TNF deletion lowers amyloid-beta "
                "production, C1q inhibition reduces synaptic engulfment, S100A9 inhibition reduces plaque), yet both models "
                "return a negative shift, away from the WT state. Removing a secreted gene may displace the embedding as if the "
                "cell had lost an activation state, which then reads as moving away from disease. This looks like a systematic "
                "limitation of using the assay directly for target discovery.",
        "fig4c": "Figure 4  Detection rate (x-axis) against the size of the shift (y-axis, log scale). Broadly expressed genes "
                 "shift more. Red points flip sign between models.",
        "r46a": "Shift size is tightly coupled to detection rate. CD8A (detected in 9% of target cells) reaches -0.098 in 104M, "
                "dwarfing everything else, but CD8A defines the T-cell lineage; the natural reading is that deleting it "
                "displaces the cell's identity, which surfaces as movement away from the WT state. Adversarial validation of "
                "in-silico deletion reports that almost no effect survives comparison with expression-matched nulls{cite:42033879,10.64898/2026.08.04.732812}, "
                "so deletion results for lineage markers should be excluded from the candidate list.",
        "r46b": "The method itself is also criticised in the literature. Evaluations of perturbation-response prediction with "
                "single-cell foundation models report that simple mean or linear baselines beat the foundation models{cite:40269681}, "
                "that attention-derived pairwise information adds no incremental value{cite:42482180}, and that zero-shot "
                "Geneformer meets the state-separation criterion for only about half of the perturbations evaluated{cite:10.1101/2025.05.11.653338}. "
                "The conservative conclusion of this report, to accept only genes where both models agree, assumes that criticism.",
        "r47a": "For the 19 target populations we checked the literature for reported changes in AD. Between-study disagreement "
                "is common and confidence varies widely by population.",
        "r47th": ["Population", "AD cells", "WT cells", "Reported in the literature", "Confidence"],
        "r47b": "One important constraint: upstream preprocessing caps each cell type at 3,000 cells, so CD4 T, CD8 T, naive B "
                "and marginal-zone B cells are exactly 3,000 in both AD and WT. Claims such as increased B-cell abundance in AD "
                "therefore cannot be tested with the cell counts in this dataset. The reported directions were used only as "
                "background for interpreting the ISP results; agreement in proportions was not part of the adjudication (Figure 5).",
        "fig5c": "Figure 5  Composition of the pool. Nearly identical in AD and WT; the upstream cap makes count ratios uninformative.",
        "r47c": "The spleen-brain axis itself is implicated in AD. Denervating the spleen breaks monocyte-microglia communication "
                "and accelerates disease{cite:42321169}, and the axis acts on immune and metabolic regulation{cite:41778859}. "
                "AD is also accompanied by disturbed bone-marrow myelopoiesis that limits monocyte homing to the brain{cite:42711429}. "
                "The cells analysed here sit on the peripheral side of that axis.",
        "d51": [
            "The choice of model changes the conclusion. The 104M/316M gap is 5-7x the measurement noise and flips the sign in 18 of 55 genes (33%), so deriving targets from one model alone is unsafe.",
            "Among the genes where the literature fixes a direction, 316M agrees markedly more often (7 to 1). TYROBP (43% detection) and GCA (bone-marrow-derived immune cells driving AD) are the strongest cases because the supporting evidence is about peripheral immune cells, the same compartment as the assay.",
            "This is not an argument to adopt 316M. Only seven genes could be adjudicated, and counterexamples remain, including LYZ and SPI1, where human genetics and direct deletion experiments disagree with each other.",
            "For complement and inflammatory mediators (C1q, TNF, IL1B, S100A8/A9) both models contradict the literature. That is most likely a limitation of the method, token deletion of a single gene followed by an embedding-distance readout, rather than a model-selection problem.",
        ],
        "d52": [
            "What the goal state means: the metric measures proximity to the WT state of splenic immune cells, not the severity of AD pathology. APOE illustrates this. It is the strongest AD risk gene, yet both models return a negative shift. Brain amyloid biology and the splenic immune state need not point the same way.",
            "Expression confounding: shift size depends strongly on detection rate (Figure 4). Validation against expression-matched nulls is essential and has not been done here.",
            "Marker genes: deleting lineage markers such as CD8A, CD4 or CD19 displaces the embedding along the lineage's own axis, so those shifts cannot be read as therapeutic targets.",
            "Statistical power: the literature fixes a direction for 17 of the 55 genes, and only 7 of those are both contested between the models and adequately expressed. The 7-to-1 count describes those seven genes and does not generalise.",
            "Agreement between the two models is not independent replication: 104M and 316M share their tokenisation and rank encoding, so effects driven by expression appear in both. Agreement is necessary but not sufficient, and independent confirmation requires another modality such as CRISPRi or Perturb-seq.",
            "Confounding in the comparison: dtype also changes between the 104M (fp32) and 316M (bf16) runs. The precision difference is smaller than the sampling difference, so its effect on the ranking was judged limited.",
        ],
        "d53": [
            "Report as candidates only the genes whose shift direction agrees across both models.",
            "Do not present single-model candidates on their own: they move with the model choice.",
            "Exclude lineage markers (CD8A, CD4, CD19, NKG7 and similar) from the target list and track population-level effects separately.",
            "State the detection rate in the target cells (for instance, below 5%) and do not interpret ranks of barely expressed genes.",
            "Treat complement and secreted-mediator results as a known weakness of the assay and prefer causal evidence from inhibition experiments in the literature.",
            "Always state what the metric measures: proximity to the WT immune state, not improvement of pathology.",
        ],
        "s6a": "Every number, figure and table in this report can be regenerated with the scripts in the repository.",
        "s6th": ["Script", "Purpose"],
        "s6rows": [
            ["analysis/16_gene_expression_by_celltype.py", "Per-cell-type, per-condition gene detection rates (expression breadth)"],
            ["analysis/17_isp_model_compare.py", "Separation of the three sources of difference, sign flips, candidate-list comparison"],
            ["analysis/18_literature_check.py", "Europe PMC retrieval with the per-gene queries built in"],
            ["analysis/19_adjudicate.py", "Merges the literature verdicts and machine-verifies the PMIDs"],
            ["analysis/20_report_figures.py", "Generates Figures 1-5 (Japanese and English)"],
            ["analysis/21_build_report_docx.py", "Generates this report in DOCX (Japanese and English)"],
        ],
        "s6b": "Raw data live under docs/quantization/isp-model-comparison/ (gene_table.csv, summary.json, adjudication.csv, literature/, figures/).",
        "s7a": "Every PMID and DOI was checked one by one against the Europe PMC API and its title, journal and year confirmed. "
               "Uncited references are not listed.",
        "pos": "deletion helps", "neg": "deletion harms", "mixed": "mixed", "none": "-",
        "mixed_long": "mixed (unresolved)", "none_long": "no causal evidence",
        "figref": "Figure ",
        "lit_hdr": "Literature",
    },
}

LIT_NOTE_EN = {
    "C1QA": "Microglial C1q deletion reduces synaptic engulfment and improves cognition",
    "C1QB": "C1q complex level: inhibition or deletion reduces synaptic engulfment and early synapse loss (Science 2016)",
    "C1QC": "Same C1q-complex evidence as C1QB; no subunit-specific data",
    "IL1B": "IL-1 blockade is reported to help, but no primary source could be pinned down here (left unresolved)",
    "S100A8": "Same direction as its obligate heterodimer partner S100A9 (calprotectin)",
    "S100A9": "Inhibiting S100A9 improves cognition and reduces amyloid-beta plaque in APP/PS1",
}

POP = {
    "ja": {
        "Ly6c.high.classical.Monocytes": ("脾臓・脳で増加。脾臓が供給源となり CCL2-CCR2 で中枢へ動員", "中〜高"),
        "Ly6c.low.nonclassical.Monocytes": ("ヒト血で古典型から中間型へのシフトとして増加", "低〜中"),
        "Macrophages": ("脾臓では数より表現型の変化が主体", "低"),
        "Neutrophils": ("脾臓で増加（3xTg-AD では明確、5xFAD は脾腫なし＝系統依存）。ヒト血はメタ解析で増加", "中〜高"),
        "DCs": ("ヒト血で mDC 減少の報告。方向は研究間で相反", "低"),
        "cDC.1": ("機能的な関与の報告はあるが、数の変化は未確立", "低"),
        "cDC.2": ("AD の脾臓・血での測定が見つからない", "なし"),
        "pDCs": ("アミロイド陽性者で増加との報告と、変化なしとの報告が対立", "低"),
        "Migratory.DCs": ("測定した研究が見つからない", "なし"),
        "CD4.T.cells": ("変化は報告されるが方向は相反", "中"),
        "CD8.T.cells": ("ヒト血で比率低下・疲弊増加。脾臓では 5xFAD で変化なし、3xTg-AD で減少（脾臓特異的変化は弱い）", "中"),
        "Immature.T.cells": ("該当集団を定量した研究が見つからない", "なし"),
        "Gamma.Delta.T.cells": ("報告が見つからない", "なし"),
        "NKT.cells": ("方向相反、研究の質も低い", "低"),
        "NK_ILC1": ("AD 血で減少＋機能変化（2コホート一致、メタ解析は有意差なし）。NK と ILC1 の分離は文献側でも未整理", "中"),
        "Naive.Memory.B.cells": ("アミロイド陽性健常者で B 細胞増加、別コホートでは減少（相反）", "相反"),
        "Marzinal.zone.B.cells": ("定量した研究が見つからない", "なし"),
        "Plasma.cells": ("マウス AD で増加の報告（脾臓は間接的）", "低"),
        "Germinal.Center.B.cells": ("直接の測定が見つからない", "なし"),
    },
    "en": {
        "Ly6c.high.classical.Monocytes": ("Expand in spleen and brain; the spleen is the reservoir, mobilised to the CNS via CCL2-CCR2", "moderate-high"),
        "Ly6c.low.nonclassical.Monocytes": ("Reported to expand in human blood as a shift from classical to intermediate", "low-moderate"),
        "Macrophages": ("In spleen mainly a phenotypic shift rather than a change in number", "low"),
        "Neutrophils": ("Expand in spleen (clear in 3xTg-AD, no splenomegaly in 5xFAD, so strain-dependent); human blood elevated in meta-analysis", "moderate-high"),
        "DCs": ("Myeloid DCs reported decreased in human blood; direction conflicting between studies", "low"),
        "cDC.1": ("Functional role reported, but no established change in abundance", "low"),
        "cDC.2": ("No AD-specific spleen or blood measurement found", "none"),
        "pDCs": ("Increased in amyloid-positive participants in one cohort, unchanged in a larger one: directly conflicting", "low"),
        "Migratory.DCs": ("No study quantifying this population found", "none"),
        "CD4.T.cells": ("Changes reported but the direction conflicts", "moderate"),
        "CD8.T.cells": ("Lower percentage and higher exhaustion in human blood; in spleen unchanged in 5xFAD and decreased in 3xTg-AD, so the spleen-specific change is weak", "moderate"),
        "Immature.T.cells": ("No study quantifying this population found", "none"),
        "Gamma.Delta.T.cells": ("No report found", "none"),
        "NKT.cells": ("Directions conflict and study quality is low", "low"),
        "NK_ILC1": ("Decreased with a functional shift in AD blood (two cohorts agree, one meta-analysis is null); NK versus ILC1 is not cleanly resolved in the literature either", "moderate"),
        "Naive.Memory.B.cells": ("B cells increased in amyloid-positive cognitively normal participants, decreased in another cohort: conflicting", "conflicting"),
        "Marzinal.zone.B.cells": ("No study quantifying this population found", "none"),
        "Plasma.cells": ("Reported to increase in mouse AD (spleen evidence indirect)", "low"),
        "Germinal.Center.B.cells": ("No direct measurement found", "none"),
    },
}


# ----------------------------------------------------------------- doc helpers
def setup(doc: Document) -> None:
    st = doc.styles["Normal"]
    st.font.name = FONT
    st.font.size = Pt(10.5)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Cm(2.0)
        s.left_margin = s.right_margin = Cm(2.2)


def h(doc: Document, text: str, level: int = 1) -> None:
    p = doc.add_heading(text, level=level)
    for r in p.runs:
        r.font.name = FONT
        r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        r.font.color.rgb = RGBColor(0x1F, 0x2A, 0x44)


def para(doc: Document, text: str, size: float = 10.5, align=None, space_after: int = 6):
    p = doc.add_paragraph()
    r = p.add_run(text)
    r.font.size = Pt(size)
    r.font.name = FONT
    r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    if align:
        p.alignment = align
    p.paragraph_format.space_after = Pt(space_after)
    return p


def bullets(doc: Document, items: list[str]) -> None:
    for it in items:
        p = doc.add_paragraph(style="List Bullet")
        r = p.add_run(it)
        r.font.size = Pt(10.5)
        r.font.name = FONT
        r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        p.paragraph_format.space_after = Pt(3)


def table(doc: Document, headers: list[str], rows: list[list[str]],
          widths: list[float] | None = None, font: float = 8.5) -> None:
    t = doc.add_table(rows=1, cols=len(headers))
    t.style = "Light Grid Accent 1"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, htxt in enumerate(headers):
        c = t.rows[0].cells[i]
        c.text = ""
        r = c.paragraphs[0].add_run(htxt)
        r.bold = True
        r.font.size = Pt(font + .5)
        r.font.name = FONT
        r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    for row in rows:
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            r = cells[i].paragraphs[0].add_run(str(v))
            r.font.size = Pt(font)
            r.font.name = FONT
            r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    if widths:
        for r_ in t.rows:
            for i, w in enumerate(widths):
                r_.cells[i].width = Cm(w)


def figure(doc: Document, lang: str, name: str, caption: str, width: float = 16.0) -> None:
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(FIGROOT / lang / name), width=Cm(width))
    c = doc.add_paragraph()
    c.alignment = WD_ALIGN_PARAGRAPH.CENTER
    r = c.add_run(caption)
    r.font.size = Pt(8.5)
    r.italic = True
    r.font.name = FONT
    r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    c.paragraph_format.space_after = Pt(10)


def render_bullet(tpl: str) -> str:
    """Bullets store citations as {cite:<key>[,<key>]} placeholders so that the
    citation order collected in pass 1 matches the order printed in pass 2."""
    def repl(m):
        keys = m.group(1)
        if not keys:
            raise SystemExit("bare {cite} without a reference key in a bullet")
        return cite(*keys.split(","))
    return re.sub(r"\{cite(?::([^}]+))?\}", repl, tpl)


# ----------------------------------------------------------------- main build
def build(lang: str) -> Document:
    L = S[lang]
    t = pd.read_csv(D / "gene_table.csv", index_col=0)
    summ = json.loads((D / "summary.json").read_text())
    adj = pd.read_csv(D / "adjudication.csv")
    expr = json.loads(pathlib.Path("docs/quantization/expression/AD_spleen.json").read_text())
    c = summ["noise_and_contrasts"]
    prec, samp, modl = c[0], c[1], c[2]
    flips = adj[adj.sign_flip]
    dec = adj[(adj.lit_verdict.isin(["POS", "NEG"])) & (adj.detection_AD >= 0.005)]
    dis = dec[dec.sign_flip]
    contra = dec[dec.better_supported_model == "neither"]
    n316 = int(dis.better_supported_model.str.startswith("316M").sum())
    n104 = int((dis.better_supported_model == "104M").sum())
    fmt = dict(model=modl["delta_as_pct_of_signal"], prec=prec["delta_as_pct_of_signal"],
               samp=samp["delta_as_pct_of_signal"], rho=modl["spearman"],
               prho=prec["spearman"], srho=samp["spearman"], flips=len(flips),
               pct=len(flips) / len(adj) * 100, ndec=len(dec), ndis=len(dis),
               n316=n316, n104=n104, ncontra=len(contra),
               p2n=int((flips.shift_104M > 0).sum()), n2p=int((flips.shift_104M < 0).sum()))

    doc = Document()
    setup(doc)

    para(doc, L["kicker"], size=11, align=WD_ALIGN_PARAGRAPH.CENTER)
    h(doc, L["title"], 0)
    # count what is actually cited (known after pass 1), not the size of the pool
    para(doc, L["subtitle"].format(n=len(NUM) if NUM else len(REFS)), size=12,
         align=WD_ALIGN_PARAGRAPH.CENTER)
    para(doc, L["meta"], size=9, align=WD_ALIGN_PARAGRAPH.CENTER, space_after=16)

    # 1 summary
    h(doc, L["s1"], 1)
    para(doc, L["sum_intro"])
    bullets(doc, [L[k].format(**fmt) for k in ("sum_b1", "sum_b2", "sum_b3", "sum_b4", "sum_b5", "sum_b6")])

    # 2 background
    h(doc, L["s2"], 1)
    for k in ("bg1", "bg2", "bg3"):
        para(doc, L[k])

    # 3 methods
    h(doc, L["s3"], 1)
    h(doc, L["s31"], 2)
    table(doc, L["m_t_h"], L["m_rows"], widths=[3.6, 12.4], font=9)
    h(doc, L["s32"], 2); para(doc, L["m32"])
    h(doc, L["s33"], 2); para(doc, L["m33"])
    h(doc, L["s34"], 2)
    for k in ("m34a", "m34b", "m34c"):
        para(doc, L[k])

    # 4 results
    h(doc, L["s4"], 1)
    h(doc, L["s41"], 2)
    table(doc, L["r_th"], [[x["label"], x["n_genes"], f"{x['spearman']:+.3f}",
                            f"{x['sign_flips']}/{x['n_genes']}", f"{x['delta_as_pct_of_signal']:.1f}%"]
                           for x in (prec, samp, modl)], widths=[6.0, 2.0, 2.4, 2.4, 3.6], font=9)
    para(doc, "")
    figure(doc, lang, "fig2_noise_ladder.png", L["fig1c"], width=15.0)
    para(doc, L["r41a"].format(**fmt))
    para(doc, L["r41b"])
    figure(doc, lang, "fig1_scatter.png", L["fig2c"], width=16.0)

    h(doc, L["s42"], 2)
    para(doc, L["r42a"].format(**fmt))
    table(doc, L["r42th"], [[r.gene, f"{r.detection_AD*100:.1f}%", f"{r.shift_104M:+.5f}", int(r.rank_104M),
                             f"{r.shift_316M:+.5f}", int(r.rank_316M),
                             {"POS": L["pos"], "NEG": L["neg"], "MIXED": L["mixed"], "NONE": L["none"]}[r.lit_verdict]]
                            for _, r in flips.sort_values("detection_AD", ascending=False).iterrows()],
          widths=[2.2, 1.8, 3.0, 1.4, 3.0, 1.4, 2.6], font=8)
    para(doc, L["r42note"], size=8.5)

    h(doc, L["s43"], 2)
    para(doc, L["r43a"])
    f3 = [x for x in summ["expression_filters"] if abs(x["min_detection"] - 0.03) < 1e-9][0]
    ab = adj.set_index("gene")
    lit_short = {"POS": L["pos"], "NEG": L["neg"], "MIXED": L["mixed_long"], "NONE": L["none_long"]}
    rows43 = []
    # flip_pos_to_neg = candidate in 104M but not in 316M, so those labels belong in
    # the 104M / 316M columns respectively (and the reverse for flip_neg_to_pos).
    for g in f3["flip_pos_to_neg"]:
        r = ab.loc[g]
        rows43.append([g, L["r43cand104"], L["r43na"], f"{r.detection_AD*100:.0f}%",
                       lit_short[r.lit_verdict], f"{r.shift_104M:+.5f} / {r.shift_316M:+.5f}"])
    for g in f3["flip_neg_to_pos"]:
        r = ab.loc[g]
        rows43.append([g, L["r43na"], L["r43cand316"], f"{r.detection_AD*100:.0f}%",
                       lit_short[r.lit_verdict], f"{r.shift_104M:+.5f} / {r.shift_316M:+.5f}"])
    table(doc, L["r43th"], rows43, widths=[2.0, 3.2, 3.2, 1.6, 2.6, 3.0], font=8)
    para(doc, L["r43b"])

    h(doc, L["s44"], 2)
    para(doc, L["r44a"].format(**fmt))
    rows44 = []
    for _, r in dis.iterrows():
        src = (L["pos"] if r.lit_verdict == "POS" else L["neg"])
        if key_of(r.lit_pmid):
            src += " " + cite(key_of(r.lit_pmid))
        rows44.append([r.gene, f"{r.shift_104M:+.5f}", f"{r.shift_316M:+.5f}", src,
                       r.better_supported_model.replace("both", "both" if lang == "en" else "両方")])
    table(doc, L["r44th"], rows44, widths=[2.2, 2.4, 2.4, 5.4, 2.4], font=8.5)
    para(doc, L["r44sens"].format(cite=cite("28100745")))
    para(doc, "")
    figure(doc, lang, "fig3_decisive_flips.png", L["fig3c"], width=16.0)
    bullets(doc, [render_bullet(b) for b in L["bullets44"]])

    h(doc, L["s45"], 2)
    para(doc, L["r45a"].format(n=len(contra)))
    table(doc, L["r45th"], [[r.gene, f"{r.detection_AD*100:.1f}%", f"{r.shift_104M:+.5f}", f"{r.shift_316M:+.5f}",
                             ((LIT_NOTE_EN.get(r.gene) or r.lit_note)[:70] if lang == "en" else r.lit_note[:60])
                             + (cite(key_of(r.lit_pmid)) if key_of(r.lit_pmid) else "")]
                            for _, r in contra.iterrows()],
          widths=[2.2, 1.8, 2.4, 2.4, 6.0], font=8)
    para(doc, L["r45b"])

    h(doc, L["s46"], 2)
    figure(doc, lang, "fig4_detection_vs_shift.png", L["fig4c"], width=16.0)
    para(doc, render_bullet(L["r46a"]))
    para(doc, render_bullet(L["r46b"]))

    h(doc, L["s47"], 2)
    para(doc, L["r47a"])
    pop = POP[lang]
    table(doc, L["r47th"], [[p, int(v["AD"]), int(v["WT"]), pop[p][0], pop[p][1]]
                            for p, v in sorted(expr["cell_composition"].items(),
                                               key=lambda kv: kv[1]["AD"] / (kv[1]["AD"] + kv[1]["WT"]))],
          widths=[4.6, 2.2, 2.2, 6.0, 1.6], font=8)
    para(doc, L["r47b"])
    figure(doc, lang, "fig5_cell_composition.png", L["fig5c"], width=15.0)
    para(doc, render_bullet(L["r47c"]))

    # 5 discussion
    h(doc, L["s5"], 1)
    h(doc, L["s51"], 2); bullets(doc, L["d51"])
    h(doc, L["s52"], 2); bullets(doc, L["d52"])
    h(doc, L["s53"], 2); bullets(doc, L["d53"])

    # 6 reproduce
    h(doc, L["s6"], 1)
    para(doc, L["s6a"], space_after=4)
    table(doc, L["s6th"], L["s6rows"], widths=[7.6, 8.4], font=8.5)
    para(doc, L["s6b"], size=9.5)

    # 7 references
    doc.add_page_break()
    h(doc, L["s7"], 1)
    para(doc, L["s7a"], size=9.5)
    for key in ORDER:
        p = doc.add_paragraph()
        r = p.add_run(f"[{NUM.get(key, 0)}] {BY_KEY[key]}")
        r.font.size = Pt(9)
        r.font.name = FONT
        r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        p.paragraph_format.space_after = Pt(3)
    return doc


def main() -> None:
    global MODE, NUM, ORDER
    ap = argparse.ArgumentParser()
    ap.add_argument("--lang", default="both", choices=["ja", "en", "both"])
    a = ap.parse_args()
    for lang in (["en", "ja"] if a.lang == "both" else [a.lang]):
        MODE, ORDER, NUM = "collect", [], {}
        USED.clear()
        build(lang)                                        # pass 1: citation order
        NUM = {k: i + 1 for i, k in enumerate(ORDER)}
        MODE = "write"
        doc = build(lang)                                  # pass 2: consecutive numbers
        out = D / OUTNAME[lang]
        doc.save(out)
        print(f"[{lang}] wrote {out}  ({out.stat().st_size/1024:.0f} KiB)  refs: {len(ORDER)}")


if __name__ == "__main__":
    main()
