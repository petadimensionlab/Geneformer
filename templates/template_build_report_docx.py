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
        "title": "AD_spleen ISP: 104M と 316M の比較（文献照合つき）",
        "sub": "脾臓免疫細胞 19 集団・55 遺伝子×3 時点・文献 {n} 件（PMID 検証済み）",
        "h1": "1. 狙い",
        "aim": "モデルを 104M から 316M に替えると、遺伝子欠失（ISP）の結果は質的に変わるのか。"
               "変わるなら、どちらの向きが公表された文献と合うのか。これだけを確かめた。",
        "h2": "2. やったこと",
        "did": [
            "同じ細胞・同じ遺伝子・同じ時点で 104M と 316M の結果を突き合わせた（55 遺伝子、3m/4.5m/6m）。",
            "その差が計算誤差や細胞の選び方によるものかを、条件を1つずつ変えた比較で切り分けた。",
            "標的細胞での遺伝子の検出率を計算し、そもそも動かせない遺伝子を判定から外した。",
            "各遺伝子について Europe PMC で文献を検索し、得た PMID を1件ずつ引き直して検証した。",
            "文献の向き（欠失で改善か悪化か）と 2 モデルの向きを突き合わせた。",
        ],
        "h3": "3. 結果",
        "t1h": ["条件を1つだけ変えた比較", "遺伝子数", "順位相関", "符号が逆になった数", "差の大きさ"],
        "t1note": "差の大きさ = 平均 |差| ÷ 平均 |シフト|（%）。数値が大きいほど、その要素の影響が大きい。"
                  "符号が逆 = 同じ遺伝子で 104M と 316M の正負が逆。",
        "res": [
            "モデル差は 84.9%。計算精度差（12.3%）の約7倍、細胞サンプリング差（17.5%）の約5倍。"
            "**ノイズではなく、モデル選択が結果を左右している。**",
            "55 遺伝子中 18 個（33%）で符号が逆になった。発現が十分な遺伝子（検出率3%以上、44個）に絞っても 13 個が逆転する。",
            "文献が向きを決められる遺伝子のうち、2 モデルが争っているのは 7 個。316M が文献に一致 6 個、104M が 1 個（LYZ）。",
        ],
        "h4": "4. 文献との照合（争いのある 7 遺伝子）",
        "t2h": ["遺伝子", "検出率", "104M", "316M", "文献の向き（出典）", "文献と一致"],
        "t2note": "「文献の向き」= その遺伝子を欠失させると病態が改善するか悪化するか。"
                  "「文献と一致」= モデルの向きがそれと合っているか。",
        "h5": "5. 両モデルが文献と食い違う 7 遺伝子",
        "t3h": ["遺伝子", "検出率", "104M", "316M", "文献の向き（出典）"],
        "t3note": "文献はいずれも「欠失で改善」。ところが両モデルとも負（＝WT から遠ざかる）を返した。"
                  "分泌性の炎症メディエーターと補体に偏っており、手法側の弱点と考えられる。",
        "h6": "6. 考察と推奨",
        "disc": [
            "**モデルを替えると結論が変わる。** 片方のモデルだけで標的候補を出すのは危険。",
            "**317M の方が文献に合う**（7 対 1）。ただし判定できたのは 7 遺伝子だけなので、これを「316M が正しい」の証明とは扱わない。",
            "**両モデルで同じ向きの遺伝子だけを候補として報告する。** 片方だけの結果はモデル選択で入れ替わる。",
            "**系統マーカー（CD8A など）は候補から外す。** 削除すると細胞の種類そのものが崩れる方向に動くため、治療標的として読めない。",
            "**検出率の低い遺伝子は順位を解釈しない。** 摂動できる細胞が 500 未満の遺伝子が 9 個あり（TREM2・SPP1・CCL2 など）、チームの合格基準 C4 では不合格になる。",
        ],
        "h7": "7. 用語（この文書で使うラベル）",
        "terms": [
            ("Shift（Shift_to_goal_end）", "遺伝子を1つ消したとき、細胞の状態が健常（WT）の状態にどれだけ近づいたかを表す値。正 = 近づく（＝その遺伝子は病気を進める側）、負 = 遠ざかる。"),
            ("符号が逆（反転）", "同じ遺伝子で 104M と 316M の Shift の正負が食い違うこと。順位が動くことより深刻。"),
            ("モデル差・精度差・サンプリング差", "条件を1つだけ変えた比較の呼び名。モデル差 = 104M 対 316M、精度差 = fp32 対 bf16、サンプリング差 = 200 細胞 対 100 細胞。"),
            ("検出率", "標的細胞のうち、その遺伝子を持つ細胞の割合。低いと動かせる細胞が少なく、Shift の意味が薄れる。"),
            ("文献の向き", "その遺伝子を欠失・阻害すると病態が改善するか悪化するか。改善 → Shift 正を期待、悪化 → 負を期待。"),
            ("C4", "チームが事前に定めた合格基準の1つ。摂動できる細胞が 500 以上であること。未達の結果は解釈しない。"),
            ("文献の証拠の層", "第1層 = 本調査で直接取得して要旨を読んだもの。第2層 = 並列の文献調査の結果を PMID で機械照合したもの。第2層は第1層が黙っている遺伝子にのみ使用。"),
        ],
        "h8": "8. 参考文献",
        "ref_note": "PMID・DOI は Europe PMC で1件ずつ照合し、題名・雑誌・年が一致することを確認済み。",
        "pos": "改善（正を期待）", "neg": "悪化（負を期待）", "mixed": "相反", "none": "—",
        "cap_t1": "表1　条件を1つずつ変えた比較。モデル差（赤）だけが突出している。",
        "cap_t2": "表2　両モデルが争い、文献が向きを示している 7 遺伝子。",
        "cap_f1": "図1　正体を分けた結果。モデル差はサンプリング差の約5倍、精度差の約7倍。",
        "cap_f2": "図2　符号が逆になった 7 遺伝子。上下で横軸の縮尺が違う（LYZ が 1 桁大きいため）。",
    },
    "en": {
        "title": "AD_spleen ISP: 104M versus 316M, with the literature check",
        "sub": "19 splenic immune populations / 55 genes x 3 timepoints / {n} references (PMID-verified)",
        "h1": "1. Aim",
        "aim": "Does replacing the model (104M -> 316M) change the outcome of gene deletion (ISP) qualitatively, and if so, "
               "which model's direction agrees with the published literature? Nothing else was tested.",
        "h2": "2. What was done",
        "did": [
            "Compared the 104M and 316M outputs on the same cells, genes and timepoints (55 genes, 3/4.5/6 months).",
            "Split that difference into precision, cell sampling and model by varying one condition at a time.",
            "Computed per-gene detection rates in the target cells and dropped genes that cannot move anyway.",
            "Searched Europe PMC per gene and re-fetched every PMID to verify it.",
            "Matched the literature direction (does deletion help or harm) against each model's direction.",
        ],
        "h3": "3. Results",
        "t1h": ["One condition varied", "Genes", "Rank corr.", "Sign flips", "Size of the gap"],
        "t1note": "Size of the gap = mean |difference| / mean |shift| (%). Larger means that factor matters more. "
                  "Sign flip = the two models disagree on the sign for the same gene.",
        "res": [
            "The model gap is 84.9%: about 7x the precision difference (12.3%) and about 5x the cell-sampling "
            "difference (17.5%). **This is not noise; the model choice drives the result.**",
            "18 of 55 genes (33%) flip sign. Restricting to well-expressed genes (detected in >=3% of target cells, 44 genes) still leaves 13 flips.",
            "Where the literature fixes a direction, the two models disagree in 7 genes: 316M agrees with the literature in 6, 104M in 1 (LYZ).",
        ],
        "h4": "4. Literature check: the 7 contested genes",
        "t2h": ["Gene", "Detection", "104M", "316M", "Literature direction (source)", "Agrees"],
        "t2note": "Literature direction = whether deleting the gene improves or worsens the pathology. "
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
            "**Do not interpret ranks of barely expressed genes.** Nine genes can be perturbed in fewer than 500 cells (TREM2, SPP1, CCL2 and others), which fails the team's C4 criterion.",
        ],
        "h7": "7. Terms used in this document",
        "terms": [
            ("Shift (Shift_to_goal_end)", "How much closer a cell moves to the healthy (WT) state when one gene is deleted. Positive = closer (the gene pushes the disease), negative = further away."),
            ("Sign flip", "The two models disagree on the sign of the Shift for the same gene. More serious than a mere rank change."),
            ("Model / precision / sampling gap", "Names for the one-condition-at-a-time comparisons: model = 104M vs 316M, precision = fp32 vs bf16, sampling = 200 vs 100 cells."),
            ("Detection rate", "Share of target cells carrying that gene. Low values mean few cells can be perturbed, so the Shift carries little meaning."),
            ("Literature direction", "Whether deleting or inhibiting the gene improves or worsens the pathology. Improves -> expect a positive Shift; worsens -> expect negative."),
            ("C4", "One of the team's pre-agreed pass criteria: at least 500 perturbable cells. Results below it are not interpreted."),
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
        rows.append([r.gene, f"{r.detection_AD*100:.0f}%", f"{r.shift_104M:+.5f}", f"{r.shift_316M:+.5f}",
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
        rows.append([r.gene, f"{r.detection_AD*100:.1f}%", f"{r.shift_104M:+.5f}", f"{r.shift_316M:+.5f}", src])
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
    note = ("注意: 上流の前処理で細胞種ごとに最大 3,000 セルに間引かれているため、"
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
