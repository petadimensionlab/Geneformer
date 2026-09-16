#!/usr/bin/env python
"""Build the compact comparison report (DOCX) in Japanese and English.

TEMPLATE — a starting point (ひな形), not a finished analysis.
The numbers and literature verdicts baked into this file are a worked example
from the AD_spleen 104M-vs-316M check (2026-09). Re-derive every number and
re-verify every PMID before reusing it on another tissue or model.

Usage: .venv/bin/python templates/template_build_report_docx.py [--lang ja|en|both]

Outputs (docs/quantization/isp-model-comparison/):
    REPORT-ISP-104M-vs-316M.docx      (+ .pdf)  English
    REPORT-ISP-104M-vs-316M-jp.docx   (+ .pdf)  Japanese

Layout is deliberately short: aim -> what was done -> numbers -> literature ->
discussion -> terms -> references. Background history is not repeated here.
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

REFS = [
    ("40931065", "ABCA7 variants impact phosphatidylcholine and mitochondria in neurons. (2025) Nature. PMID 40931065"),
    ("41125877", "Transcriptional and epigenetic targets of MEF2C in human microglia. (2025) Nat Immunol. PMID 41125877"),
    ("41620758", "Protective PLCG2 variants associate with a delayed onset of Alzheimer's disease. (2026) Alzheimers Res Ther. PMID 41620758"),
    ("41789102", "TREM2 and microglial immunity in Alzheimer's disease. (2026) Front Immunol. PMID 41789102"),
    ("28612290", "Deficiency of TYROBP is neuroprotective in a mouse model of early Alzheimer's pathology. (2017) Acta Neuropathol. PMID 28612290"),
    ("30283031", "Deficiency of TYROBP in a tauopathy mouse model reduces C1q and normalizes clinical phenotype. (2019) Mol Psychiatry. PMID 30283031"),
    ("37949676", "Bone Marrow-Derived GCA+ Immune Cells Drive Alzheimer's Disease Progression. (2023) Adv Sci. PMID 37949676"),
    ("39555667", "Microglial Lyzl4 Facilitates β-Amyloid Clearance in Alzheimer's Disease. (2025) Adv Sci. PMID 39555667"),
    ("28100745", "Disease Progression-Dependent Effects of TREM2 Deficiency in a Mouse Model of Alzheimer's Disease. (2017) J Neurosci. PMID 28100745"),
    ("42676482", "Direct toll-like receptor 4 inhibition in Alzheimer's disease models: a scoping review. (2026) Front Aging Neurosci. PMID 42676482"),
    ("41152873", "Soluble CSF1R promotes microglial activation and amyloid clearance in Alzheimer's disease. (2025) J Neuroinflammation. PMID 41152873"),
    ("40544243", "Edaravone-Dexborneol slows cognitive decline via inhibiting S100A9. (2025) Alzheimers Res Ther. PMID 40544243"),
    ("42444329", "Young Adult Microglial Deletion of C1q Reduces Engulfment of Synapses. (2026) Glia. PMID 42444329"),
    ("27033548", "Complement and microglia mediate early synapse loss in Alzheimer mouse models. (2016) Science. PMID 27033548"),
    ("28826177", "Genetic Deletion of TNF-α Attenuates Amyloid-β Production in the 5XFAD Model. (2017) J Alzheimers Dis. PMID 28826177"),
    ("42033879", "Adversarial validation of in silico perturbation profiles. (2026) Comput Biol Chem. PMID 42033879"),
    ("10.64898/2026.08.04.732812", "A confound-diagnostic toolkit for in silico perturbation with single-cell foundation models. (2026) preprint, doi:10.64898/2026.08.04.732812"),
    ("40269681", "Benchmarking foundation cell models for post-perturbation RNA-seq prediction. (2025) BMC Genomics. PMID 40269681"),
    ("42321169", "Disruption of the brain-spleen axis impairs monocyte-microglia communication. (2026) Nat Commun. PMID 42321169"),
    ("41778859", "The spleen-brain axis in Alzheimer's disease and related dementias. (2026) Alzheimers Dement. PMID 41778859"),
    ("42711429", "Bone marrow myelopoiesis dysfunction in Alzheimer's disease limits monocyte homing. (2026) Nat Neurosci. PMID 42711429"),
]
BY_KEY = dict(REFS)


def _extend_refs() -> None:
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
    return "".join(ids)


# ----------------------------------------------------------------- content
S = {
    "ja": {
        "title": "脾臓データ（AD_spleen）の遺伝子削除シミュレーション: V2-104M と V2-316M の比較（文献照合つき）",
        "sub": "脾臓の免疫細胞 19 細胞型・55 個の遺伝子×3 時点・文献 {n} 件（PMID 検証済み）",
        "h1": "1. 狙い",
        "aim": "先に前提を説明します。AD_spleen は脾臓の単一細胞データ（AD と健常の比較）を指します。"
               "ISP は in silico perturbation の略で、計算機上で遺伝子を1つ削除する解析です。"
               "V2-104M と V2-316M は Geneformer の事前学習モデルで、数字はパラメータ数（1億400万と3億1600万）を表します。"
               "fp32 は 32 ビットの浮動小数点、bf16 は bfloat16（16 ビット）の計算精度です。"
               "確かめたかったことは2つあります。1つ目は、使うモデルを V2-104M から V2-316M に替えると、"
               "遺伝子を削除したときの結果の向き（改善か悪化か）と順位が変わるのかどうかです。2つ目は、変わる場合に、"
               "どちらのモデルが示す向きが公表された文献と合っているのかです。この2点だけを調べました。",
        "h2": "2. やったこと",
        "did": [
            "同じ細胞、同じ遺伝子、同じ時点で、104M と 316M の結果を突き合わせました"
            "（55 個の遺伝子、3 か月齢・4.5 か月齢・6 か月齢）。",
            "差の原因を切り分けるために、数値精度だけを替えた比較と、細胞数だけを替えた比較を行いました。"
            "モデルを替えた比較では、手元にある実行が数値精度と細胞数も違うため、それらを固定できていません。",
            "対象の細胞でその遺伝子がどのくらい発現しているかを計算し、そもそも削除しても動かせない遺伝子は"
            "判定から外しました。",
            "遺伝子ごとに Europe PMC（文献データベース）で文献を検索し、得られた PMID（文献の識別番号）を1件ずつ引き直して、"
            "題名と内容が主張と合っているかを確認しました。",
            "文献が示す向き（削除すると病態が改善するのか、悪化するのか）と、2つのモデルが示した向きを"
            "突き合わせました。",
        ],
        "h3": "3. 結果",
        "t1h": ["変えた条件", "遺伝子数", "順位相関（スピアマン）", "符号が反転した数", "差の大きさ（%）"],
        "t1note": "「差の大きさ」は、条件を変えたときの平均の差を、平均のシフトで割った値（%）です。"
                  "大きいほど、その条件が結果に与える影響が大きいことを意味します。"
                  "「順位相関」はスピアマンの順位相関係数で、1.0 に近いほど順位が同じことを表します。"
                  "「符号が反転した数」は、その行で比べた2つの結果が正と負で食い違った遺伝子の数です。"
                  "分母は「遺伝子数」の列の値です。"
                  "1行目と2行目は、比べる条件以外を固定して測りました。"
                  "3行目は、手元にある V2-104M と V2-316M の実行が数値精度と細胞数も違うため、他の条件を固定できていません。"
                  "したがって 84.9% はモデル差の上限として読んでください。",
        "res": [
            "モデルを替えたときの差は、シグナルの 84.9% でした（精度と細胞数も同時に違うため、上限の値です）。これは、計算精度だけを替えたときの差"
            "（12.3%）の約7倍、細胞の選び方だけを変えたときの差（17.5%）の約5倍にあたります。"
            "**つまりこの違いは測定のばらつきではなく、モデルの選択が結果を左右しています。**",
            "55 個の遺伝子のうち 18 個（33%）で、2つのモデルが正と負で逆の結論を出しました（符号の反転）。"
            "対象の細胞でよく発現している遺伝子（検出率 3% 以上、44 個）だけに絞っても、13 個で逆のままでした。",
            "文献が向きを決められる遺伝子のうち、2つのモデルが争っているのは 7 個でした。"
            "この 7 個では、316M が文献と一致したのが 6 個、104M が一致したのが 1 個（LYZ）でした。",
        ],
        "h4": "4. 文献との照合（争いのある 7 遺伝子）",
        "t2h": ["遺伝子", "検出率", "104M", "316M", "文献の向き（出典）", "文献と一致"],
        "t2note": "Shift の値は有効数字 3 桁の指数表記で書いています（例: -5.50e-04）。"
                  "「文献の向き」は、その遺伝子を削除すると病態が改善するのか、悪化するのかを示します。"
                  "「文献と一致」は、そのモデルが示した向きが文献と合っているかを示します。",
        "h5": "5. 両モデルが文献と食い違う 7 遺伝子",
        "t3h": ["遺伝子", "検出率", "104M", "316M", "文献の向き（出典）"],
        "t3note": "Shift の値は表2と同じく有効数字 3 桁の指数表記です。"
                  "7 個すべてについて、文献は「削除すると改善する」と報告しています。ところが2つのモデルは"
                  "どちらも負の値を返しました。負は「健常な状態から遠ざかる」という意味です。"
                  "この 7 個は分泌性の炎症メディエーターと補体に偏っており、モデルの選択ではなく"
                  "手法そのものの弱点だと考えられます。",
        "h6": "6. 考察と推奨",
        "disc": [
            "**モデルを替えると結論が変わります。** したがって、片方のモデルだけで標的候補を出すのは危険です。",
            "**文献との一致は 316M の方が多い結果でした**（7 個中 6 個が 316M、1 個が 104M）。"
            "ただし判定できたのは 7 個だけですので、これを「316M が正しい」ことの証明として扱うことはできません。",
            "**候補として報告するのは、2つのモデルで同じ向きが出た遺伝子だけにしてください。**"
            "片方のモデルだけに出た候補は、モデルを替えると入れ替わります。",
            "**系統マーカー（CD8A など）は候補から外してください。** これらを削除すると、細胞の種類そのものが"
            "崩れる方向に結果が動くため、治療の標的として読むことはできません。",
            "**発現が少ない遺伝子の順位は解釈しないでください。** 摂動できる細胞が 500 個に満たない遺伝子が 9 個ありました"
            "（TREM2、SPP1、CCL2 など。TREM2 は 176 個）。値が安定しないため、これらの順位は結果として扱いません。",
        ],
        "h7": "7. 用語（この文書で使うラベル）",
        "terms": [
            ("Shift（Shift_to_goal_end）", "遺伝子を1つ消したとき、細胞の状態が健常（WT）の状態にどれだけ近づいたかを表す値。正 = 近づく（＝その遺伝子は病気を進める側）、負 = 遠ざかる。"),
            ("符号の反転", "同じ遺伝子で V2-104M と V2-316M の Shift の正負が食い違うこと。順位が動くことより深刻。"),
            ("モデル差・精度差・サンプリング差", "条件を1つだけ変えた比較の呼び名。モデル差 = 104M 対 316M、精度差 = fp32 対 bf16、サンプリング差 = 200 細胞 対 100 細胞。"),
            ("検出率", "標的細胞のうち、その遺伝子を持つ細胞の割合。低いと動かせる細胞が少なく、Shift の意味が薄れる。"),
            ("文献の向き", "その遺伝子を欠失・阻害すると病態が改善するか悪化するか。改善 → Shift 正を期待、悪化 → 負を期待。"),
            ("E4", "摂動できる細胞数についての項目。目安は 500 個以上だが、これは実測 2 例の間を取った線であり、合否を決める基準ではない。少ない遺伝子は順位を解釈しない。"),
            ("文献の証拠の層", "第1層 = 本調査で直接取得して要旨を読んだもの。第2層 = 並列の文献調査の結果を PMID で機械照合したもの。第2層は第1層が黙っている遺伝子にのみ使用。"),
        ],
        "h8": "8. 参考文献",
        "ref_note": "PMID・DOI は Europe PMC で1件ずつ照合し、題名・雑誌・年が一致することを確認済み。",
        "pos": "改善（正を期待）", "neg": "悪化（負を期待）", "mixed": "相反", "none": "—",
        "cap_t1": "表1　条件を1つずつ変えた比較です。モデルを替えたときの差だけが突出しています。",
        "cap_t2": "表2　両モデルが争い、文献が向きを示している 7 遺伝子。",
        "cap_f1": "図1　差の正体を分けた結果です。モデルを替えたときの差は、細胞の選び方を変えたときの差の約5倍、"
                   "計算精度を変えたときの差の約7倍あります。",
        "cap_f2": "図2　2つのモデルで符号が反転した 7 個の遺伝子です。LYZ だけ値が 1 桁大きいため、"
                   "下段は横軸の縮尺を変えています。",
    },
    "en": {
        "title": "AD_spleen ISP: 104M versus 316M, with the literature check",
        "sub": "19 splenic immune populations / 55 genes x 3 timepoints / {n} references (PMID-verified)",
        "h1": "1. Aim",
        "aim": "First, the setting. AD_spleen is a single-cell spleen dataset comparing AD with healthy controls, and ISP stands for in silico perturbation: deleting one gene on the computer. V2-104M and V2-316M are pretrained Geneformer models; the numbers are parameter counts (104 million, 316 million). fp32 is 32-bit floating point and bf16 is bfloat16, a 16-bit format. The report asks whether replacing the model (V2-104M -> V2-316M) changes the outcome of gene deletion (direction and ranking), and if so, "
               "which model's direction agrees with the published literature? Nothing else was tested.",
        "h2": "2. What was done",
        "did": [
            "Compared the 104M and 316M outputs on the same cells, genes and timepoints (55 genes, 3/4.5/6 months).",
            "Split that difference by running a precision-only comparison and a cell-count-only comparison. In the model comparison the runs in hand differ in precision and cell count as well, so those were not held fixed.",
            "Computed per-gene detection rates in the target cells and dropped genes that cannot move anyway.",
            "Searched Europe PMC, the literature database, for each gene and re-fetched every PMID (the paper identifier) to verify it.",
            "Matched the literature direction (does deletion help or harm) against each model's direction.",
        ],
        "h3": "3. Results",
        "t1h": ["Condition varied", "Genes", "Rank corr. (Spearman)", "Sign flips", "Size of the gap (%)"],
        "t1note": "Size of the gap = mean |difference| / mean |shift| (%). Larger means that factor matters more. "
                  "Rank corr. is Spearman's rank correlation: 1.0 means the two rankings are identical. "
                  "Sign flips counts the genes where the two results disagree on the sign; the denominator is the value in the Genes column. "
                  "Rows 1 and 2 hold every other condition fixed. Row 3 does not: the V2-104M and V2-316M runs in hand also differ in "
                  "precision and cell count, so read 84.9% as an upper bound on the model effect.",
        "res": [
            "The model gap is 84.9% (an upper bound: precision and cell count differ too): about 7x the precision difference (12.3%) and about 5x the cell-sampling "
            "difference (17.5%). **This is not noise; the model choice drives the result.**",
            "18 of 55 genes (33%) flip sign. Restricting to well-expressed genes (detected in >=3% of target cells, 44 genes) still leaves 13 flips.",
            "Where the literature fixes a direction, the two models disagree in 7 genes: 316M agrees with the literature in 6, 104M in 1 (LYZ).",
        ],
        "h4": "4. Literature check: the 7 contested genes",
        "t2h": ["Gene", "Detection", "104M", "316M", "Literature direction (source)", "Agrees"],
        "t2note": "Shifts are given in scientific notation with three significant digits (for example -5.50e-04). "
                  "Literature direction = whether deleting the gene improves or worsens the pathology. "
                  "Agrees = whether that model's sign matches it.",
        "h5": "5. The 7 genes where both models contradict the literature",
        "t3h": ["Gene", "Detection", "104M", "316M", "Literature direction (source)"],
        "t3note": "The literature says deletion helps in every case, yet both models return a negative shift (away from WT). "
                  "They are concentrated in secreted inflammatory mediators and complement, which points to a weakness of the assay.",
        "h6": "6. Discussion and recommendations",
        "disc": [
            "**Changing the model changes the conclusion.** Targets derived from one model alone are unsafe.",
            "**316M fits the literature better** (7 to 1), but only seven genes could be adjudicated, so this is not proof that 316M is correct.",
            "**Report as candidates only genes whose direction agrees across both models.** Single-model candidates move with the model choice.",
            "**Exclude lineage markers (CD8A and similar).** Deleting them displaces the cell's identity, which cannot be read as a therapeutic effect.",
            "**Do not interpret ranks of barely expressed genes.** Nine genes can be perturbed in fewer than 500 cells (TREM2 has 176, SPP1 has 31, CCL2 has 13). Their ranks are not treated as results because the values are unstable.",
        ],
        "h7": "7. Terms used in this document",
        "terms": [
            ("Shift (Shift_to_goal_end)", "How much closer a cell moves to the healthy (WT) state when one gene is deleted. Positive = closer (the gene pushes the disease), negative = further away."),
            ("Sign flip", "The two models disagree on the sign of the Shift for the same gene. More serious than a mere rank change."),
            ("Model / precision / sampling gap", "Names for the one-condition-at-a-time comparisons: model = 104M vs 316M, precision = fp32 vs bf16, sampling = 200 vs 100 cells."),
            ("Detection rate", "Share of target cells carrying that gene. Low values mean few cells can be perturbed, so the Shift carries little meaning."),
            ("Literature direction", "Whether deleting or inhibiting the gene improves or worsens the pathology. Improves -> expect a positive Shift; worsens -> expect negative."),
            ("E4", "The item on how many cells can be perturbed. The 500-cell guide sits between two measured cases and is not a pass/fail criterion. Genes below it are not ranked."),
            ("Evidence tiers", "Tier 1 = retrieved and read in this study. Tier 2 = a parallel sweep with PMIDs machine-verified. Tier 2 is used only where tier 1 is silent."),
        ],
        "h8": "8. References",
        "ref_note": "Every PMID and DOI was checked against Europe PMC and its title, journal and year confirmed.",
        "pos": "helps (expect +)", "neg": "harms (expect -)", "mixed": "mixed", "none": "-",
        "cap_t1": "Table 1  One condition varied at a time. Only the model gap stands out.",
        "cap_t2": "Table 2  The seven genes where the models disagree and the literature indicates a direction.",
        "cap_f1": "Figure 1  Sources of the gap. The model gap is about 5x the sampling gap and 7x the precision gap.",
        "cap_f2": "Figure 2  The seven sign-flipped genes. The two panels use different x-scales (LYZ is an order of magnitude larger).",
    },
}

POP = {
    "ja": {
        "Ly6c.high.classical.Monocytes": ("脾臓・脳で増加、脾臓が供給源", "中〜高"),
        "Ly6c.low.nonclassical.Monocytes": ("ヒト血で増加（中間型へのシフト）", "低〜中"),
        "Macrophages": ("脾臓では表現型の変化が主体", "低"),
        "Neutrophils": ("脾臓で増加（系統依存）、ヒト血でも増加", "中〜高"),
        "DCs": ("ヒト血で mDC 減少、研究間で相反", "低"),
        "cDC.1": ("機能の報告あり、数の変化は未確立", "低"),
        "cDC.2": ("AD での測定なし", "なし"),
        "pDCs": ("増加報告と変化なし報告が対立", "低"),
        "Migratory.DCs": ("測定なし", "なし"),
        "CD4.T.cells": ("変化は報告されるが方向は相反", "中"),
        "CD8.T.cells": ("ヒト血で減少・疲弊増加、脾臓では弱い", "中"),
        "Immature.T.cells": ("測定なし", "なし"),
        "Gamma.Delta.T.cells": ("報告なし", "なし"),
        "NKT.cells": ("方向相反、研究の質も低い", "低"),
        "NK_ILC1": ("AD 血で減少＋機能変化", "中"),
        "Naive.Memory.B.cells": ("増加報告と減少報告が対立", "相反"),
        "Marzinal.zone.B.cells": ("測定なし", "なし"),
        "Plasma.cells": ("マウス AD で増加（脾臓は間接）", "低"),
        "Germinal.Center.B.cells": ("測定なし", "なし"),
    },
    "en": {
        "Ly6c.high.classical.Monocytes": ("Expand in spleen and brain; spleen is the reservoir", "moderate-high"),
        "Ly6c.low.nonclassical.Monocytes": ("Expand in human blood (shift to intermediate)", "low-moderate"),
        "Macrophages": ("Mainly a phenotypic shift in spleen", "low"),
        "Neutrophils": ("Expand in spleen (strain-dependent) and in human blood", "moderate-high"),
        "DCs": ("Myeloid DCs decreased in human blood; studies conflict", "low"),
        "cDC.1": ("Functional role reported; no abundance change established", "low"),
        "cDC.2": ("No measurement in AD", "none"),
        "pDCs": ("Increased in one cohort, unchanged in another", "low"),
        "Migratory.DCs": ("No measurement", "none"),
        "CD4.T.cells": ("Changes reported, direction conflicts", "moderate"),
        "CD8.T.cells": ("Lower and more exhausted in human blood; weak in spleen", "moderate"),
        "Immature.T.cells": ("No measurement", "none"),
        "Gamma.Delta.T.cells": ("Not reported", "none"),
        "NKT.cells": ("Directions conflict, low study quality", "low"),
        "NK_ILC1": ("Decreased with functional change in AD blood", "moderate"),
        "Naive.Memory.B.cells": ("Increase and decrease both reported", "conflicting"),
        "Marzinal.zone.B.cells": ("No measurement", "none"),
        "Plasma.cells": ("Increase in mouse AD (spleen indirect)", "low"),
        "Germinal.Center.B.cells": ("No measurement", "none"),
    },
}


# ----------------------------------------------------------------- helpers
def setup(doc: Document) -> None:
    st = doc.styles["Normal"]
    st.font.name = FONT
    st.font.size = Pt(10.5)
    st.element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    for s in doc.sections:
        s.top_margin = s.bottom_margin = Cm(2.0)
        s.left_margin = s.right_margin = Cm(2.2)


def h(doc, text, level=1):
    p = doc.add_heading(text, level=level)
    for r in p.runs:
        r.font.name = FONT
        r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        r.font.color.rgb = RGBColor(0x1F, 0x2A, 0x44)


def para(doc, text, size=10.5, space_after=6):
    p = doc.add_paragraph()
    bold = text.startswith("**") and "**" in text[2:]
    while "**" in text:
        head, _, rest = text.partition("**")
        if head:
            r = p.add_run(head)
            r.font.size = Pt(size)
            r.font.name = FONT
            r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        mid, _, text = rest.partition("**")
        if mid:
            r = p.add_run(mid)
            r.bold = True
            r.font.size = Pt(size)
            r.font.name = FONT
            r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    if text:
        r = p.add_run(text)
        r.font.size = Pt(size)
        r.font.name = FONT
        r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    p.paragraph_format.space_after = Pt(space_after)
    return p


def bullets(doc, items, size=10.5):
    for it in items:
        p = doc.add_paragraph(style="List Bullet")
        while "**" in it:
            head, _, rest = it.partition("**")
            if head:
                r = p.add_run(head)
                r.font.size = Pt(size)
                r.font.name = FONT
                r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
            mid, _, it = rest.partition("**")
            if mid:
                r = p.add_run(mid)
                r.bold = True
                r.font.size = Pt(size)
                r.font.name = FONT
                r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        if it:
            r = p.add_run(it)
            r.font.size = Pt(size)
            r.font.name = FONT
            r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
        p.paragraph_format.space_after = Pt(3)


def table(doc, headers, rows, widths=None, font=8.5):
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


def caption(doc, text, size=8.5):
    c = doc.add_paragraph()
    r = c.add_run(text)
    r.font.size = Pt(size)
    r.italic = True
    r.font.name = FONT
    r._element.rPr.rFonts.set(qn("w:eastAsia"), FONT)
    c.paragraph_format.space_after = Pt(9)


def figure(doc, lang, name, cap, width=15.0):
    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.add_run().add_picture(str(FIGROOT / lang / name), width=Cm(width))
    caption(doc, cap)


def render(tpl: str) -> str:
    def repl(m):
        keys = m.group(1)
        if not keys:
            raise SystemExit("bare {cite} in a bullet")
        return cite(*keys.split(","))
    return re.sub(r"\{cite(?::([^}]+))?\}", repl, tpl)


# ----------------------------------------------------------------- build
def build(lang: str) -> Document:
    L = S[lang]
    adj = pd.read_csv(D / "adjudication.csv")
    summ = json.loads((D / "summary.json").read_text())
    expr = json.loads(pathlib.Path("docs/quantization/expression/AD_spleen.json").read_text())
    prec, samp, modl = summ["noise_and_contrasts"]
    flips = adj[adj.sign_flip]
    dis = adj[(adj.lit_verdict.isin(["POS", "NEG"])) & (adj.detection_AD >= 0.005) & adj.sign_flip]
    contra = adj[(adj.lit_verdict == "POS") & (adj.detection_AD >= 0.005)
                 & (~adj.sign_flip) & (adj.shift_104M < 0) & (adj.shift_316M < 0)]
    LIT = {"POS": L["pos"], "NEG": L["neg"], "MIXED": L["mixed"], "NONE": L["none"]}

    doc = Document()
    setup(doc)
    h(doc, L["title"], 0)
    para(doc, L["sub"].format(n=len(NUM) if NUM else len(REFS)), size=9)

    h(doc, L["h1"], 1)
    para(doc, L["aim"])

    h(doc, L["h2"], 1)
    bullets(doc, L["did"])

    h(doc, L["h3"], 1)
    table(doc, L["t1h"], [[x["label"], x["n_genes"], f"{x['spearman']:+.3f}",
                           f"{x['sign_flips']}/{x['n_genes']}", f"{x['delta_as_pct_of_signal']:.1f}%"]
                          for x in (prec, samp, modl)], widths=[6.4, 1.8, 2.2, 2.6, 2.8], font=9)
    caption(doc, L["cap_t1"])
    bullets(doc, [L["res"][0]])
    para(doc, L["t1note"], size=8.5)
    figure(doc, lang, "fig2_noise_ladder.png", L["cap_f1"], width=14.5)
    bullets(doc, L["res"][1:])

    h(doc, L["h4"], 1)
    para(doc, L["t2note"], size=9)
    rows = []
    for _, r in dis.sort_values("detection_AD", ascending=False).iterrows():
        src = (L["pos"] if r.lit_verdict == "POS" else L["neg"])
        if key_of(r.lit_pmid):
            src += " " + cite(key_of(r.lit_pmid))
        rows.append([r.gene, f"{r.detection_AD*100:.1f}%", f"{r.shift_104M:+.2e}", f"{r.shift_316M:+.2e}",
                     src, r.better_supported_model])
    table(doc, L["t2h"], rows, widths=[2.0, 1.7, 2.4, 2.4, 5.4, 2.1], font=8.5)
    caption(doc, L["cap_t2"])
    figure(doc, lang, "fig3_decisive_flips.png", L["cap_f2"], width=15.0)

    h(doc, L["h5"], 1)
    rows = []
    for _, r in contra.iterrows():
        src = (L["pos"] if r.lit_verdict == "POS" else L["neg"])
        if key_of(r.lit_pmid):
            src += " " + cite(key_of(r.lit_pmid))
        rows.append([r.gene, f"{r.detection_AD*100:.1f}%", f"{r.shift_104M:+.2e}", f"{r.shift_316M:+.2e}", src])
    table(doc, L["t3h"], rows, widths=[2.2, 1.8, 2.4, 2.4, 6.2], font=8.5)
    para(doc, L["t3note"], size=9)

    h(doc, L["h6"], 1)
    bullets(doc, L["disc"])

    h(doc, L["h7"], 1)
    table(doc, ["用語" if lang == "ja" else "Term", "説明" if lang == "ja" else "Meaning"],
          [[k, v] for k, v in L["terms"]], widths=[4.4, 11.6], font=8.5)

    # cell populations, compact
    h(doc, "付録. 対象とした 19 細胞集団の文献状況" if lang == "ja"
         else "Appendix. The 19 target cell populations in the literature", 1)
    pop = POP[lang]
    table(doc, (["細胞集団", "AD", "WT", "文献での報告", "確度"] if lang == "ja"
                else ["Population", "AD", "WT", "Reported in the literature", "Confidence"]),
          [[p, int(v["AD"]), int(v["WT"]), pop[p][0], pop[p][1]]
           for p, v in sorted(expr["cell_composition"].items(),
                              key=lambda kv: kv[1]["AD"] / (kv[1]["AD"] + kv[1]["WT"]))],
          widths=[4.4, 1.7, 1.7, 6.6, 1.6], font=8)
    note = ("注意: 上流の前処理で細胞型ごとに最大 3,000 セルに間引かれているため、"
               "このデータの細胞数比で文献の主張（例: AD で B 細胞が増える）は検証できません。"
               "脾臓と脳をつなぐ軸自体は AD に関与することが報告されています{cite:42321169}。AD では骨髄の造血が乱れ、"
               "単球の中枢への移動が損なわれます{cite:42711429}。"
               if lang == "ja" else
               "Note: upstream preprocessing caps each cell type at 3,000 cells, so claims such as increased B cells in AD "
               "cannot be tested with these counts. The spleen-brain axis itself is implicated in AD{cite:42321169}, and disturbed "
               "bone-marrow myelopoiesis limits monocyte homing to the brain{cite:42711429}.")
    para(doc, render(note), size=8.5)

    doc.add_page_break()
    h(doc, L["h8"], 1)
    para(doc, L["ref_note"], size=9)
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
        build(lang)
        NUM = {k: i + 1 for i, k in enumerate(ORDER)}
        MODE = "write"
        doc = build(lang)
        out = D / OUTNAME[lang]
        doc.save(out)
        print(f"[{lang}] wrote {out}  ({out.stat().st_size/1024:.0f} KiB)  refs: {len(ORDER)} -> {sorted(NUM.values())[:3]}...")


if __name__ == "__main__":
    main()
