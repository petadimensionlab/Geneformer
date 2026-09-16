#!/usr/bin/env python
"""Merge the 104M/316M ISP comparison with the literature verdict per gene.

TEMPLATE — a starting point (ひな形), not a finished analysis.
The numbers, literature verdicts and prose baked into this file are a worked
example from the AD_spleen 104M-vs-316M check (2026-09). Re-derive every number
and re-verify every PMID before reusing it on another tissue or model.

Inputs (all produced earlier in this task):
  docs/quantization/isp-model-comparison/gene_table.csv      (17_isp_model_compare.py)
  docs/quantization/expression/AD_spleen.json                (16_gene_expression_by_celltype.py)
  docs/quantization/isp-model-comparison/literature/*.json   (18_literature_check.py + subagents)
  docs/quantization/isp-model-comparison/literature/verified-pmids.json

Two tiers of literature evidence are kept apart on purpose:
  tier 1 (mine)      — PMIDs I retrieved and read the abstract of myself; quotation-ready
  tier 2 (subagent)  — PMIDs taken from the delegated literature sweep and re-verified
                       against Europe PMC by ID; used only where tier 1 is silent
Anything unverified is reported as "no evidence" rather than guessed.

Output: adjudication.csv + adjudication.json (+ a printed summary).
"""
from __future__ import annotations

import json
import pathlib
import re
import urllib.request

import numpy as np
import pandas as pd

D = pathlib.Path("docs/quantization/isp-model-comparison")
LIT = D / "literature"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

# ---------------------------------------------------------------- tier 1 verdicts
# direction: what deleting the gene is expected to do to the (WT-ward) shift.
#   POS = deletion is expected to move cells toward the healthy state
#   NEG = deletion is expected to be harmful
#   MIXED = published evidence contradicts itself
#   NONE = no causal evidence (marker / lineage only)
T1: dict[str, tuple[str, str, str]] = {
    "ABCA7": ("NEG", "40931065", "機能喪失変異が AD リスクを大きく上昇させる（ABCA7 は保護的）"),
    "MEF2C": ("NEG", "41125877", "MEF2C はマイクログリア炎症の抑制因子で、低下が疾患に寄与する"),
    "PLCG2": ("NEG", "41620758", "P522R は軽度亢進型（hypermorph）で保護的。欠失は逆方向"),
    "SPI1": ("MIXED", "41928507", "Spi1 欠失はアミロイド病理を悪化（Neuron 2026）vs PU.1 低値が保護的というヒト遺伝学"),
    "TYROBP": ("POS", "28612290", "Tyrobp 欠損は APP/PSEN1 マウスで神経保護的"),
    "CSF1R": ("MIXED", "41152873", "可溶性 CSF1R はアミロイド除去を促進（阻害は有害）vs 持続阻害はプラーク形成を抑制"),
    "GCA": ("POS", "37949676", "GCA+ 骨髄由来免疫細胞が AD 進行を駆動し、GCA 投与でプラークが悪化"),
    "LYZ": ("NEG", "39555667", "Lyzl4（リゾチーム様酵素）は Aβ 除去を促進。ただし LYZ 本体の証拠ではない"),
    "TLR4": ("MIXED", "42676482", "TLR4 阻害は改善と悪化の双方が報告されている（scoping review）"),
    "S100A9": ("POS", "40544243", "S100A9 阻害で認知低下と Aβ プラークが改善（APP/PS1）"),
    "S100A8": ("POS", "40544243", "S100A9 とのヘテロ二量体（カルプロテクチン）として同方向"),
    "C3": ("MIXED", "28566429", "C3 欠損は保護（Sci Transl Med 2017, PMID 28566429）vs 欠損がプラークを加速（J Neurosci 2008, PMID 18562603）"),
    "CD33": ("POS", "23623698", "CD33 は Aβ 取り込みを阻害する。欠失で取り込み増、アミロイド減少"),
    "TREM2": ("NEG", "41789102", "機能喪失変異は免疫監視を破綻させアミロイド病理を悪化させる"),
    "NLRP3": ("POS", "41491073", "阻害または遺伝的欠失でアミロイド負荷と認知障害が改善"),
    "SPP1": ("POS", "36730200", "OPN の遺伝的除去・抗 OPN 抗体でプラークと認知が改善（5xFAD）"),
    "C1QA": ("POS", "42444329", "マイクログリア C1q 欠失はシナプス貪食を減らし認知を改善"),
    "CLU": ("MIXED", "", "KO の方向について一致した一次証拠を得られず"),
    "APOE": ("MIXED", "40275327", "AD 最大のリスク遺伝子だが、末梢免疫細胞の状態に対する欠失方向の一次証拠が得られず"),
    "C1QB": ("POS", "27033548", "C1q 複合体レベル：阻害・欠失でシナプス貪食と早期シナプス損失が減る（Science 2016）"),
    "C1QC": ("POS", "27033548", "C1QB と同じ C1q 複合体の知見。サブユニット別の証拠ではない"),
    "IL1B": ("POS", "", "IL-1 遮断は改善方向とされるが、本調査で一次文献を確定できず（判定は保留）"),
    "CD8A": ("NONE", "42033879", "系統マーカー。削除の効果は系統同一性アーティファクトと解釈すべき"),
    "CCL2": ("NONE", "", "標的細胞での検出率 0.04%。シグナルを担えない"),
    "CD34": ("NONE", "", "標的細胞での検出率 0.1%。シグナルを担えない"),
}
# lineaeage / marker genes with no causal role -> tier 1 NONE
for g in ["CD3D", "CD3E", "CD4", "CD19", "MS4A1", "NKG7", "KLRD1", "PRF1", "CD34",
          "ITGA2B", "XCR1", "IL3RA", "CLEC9A", "VWF", "CSF3R", "CD68", "AIF1",
          "ITGAM", "ITGAX", "MS4A7", "CEACAM1", "FGR", "CD74", "IRF7", "ISG15",
          "THBS1", "APOC1"]:
    T1.setdefault(g, ("NONE", "", "因果的な AD 証拠なし（マーカー／系統遺伝子）"))

EXTRA_SRC = {  # tier-1 sources that carry a second conflicting paper
    "C3": ["28566429", "18562603"],
    "SPI1": ["41928507", "41821106"],
}


def verify_pmids(pmids: list[str]) -> dict[str, dict]:
    """Confirm each PMID exists in Europe PMC and capture its real title."""
    out: dict[str, dict] = {}
    for p in pmids:
        url = (f"{API}?query=EXT_ID:{p}%20AND%20SRC:MED&format=json&resultType=core")
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                d = json.load(r)
            res = d.get("resultList", {}).get("result", [])
            if not res:
                continue
            rec = res[0]
            ji = rec.get("journalInfo", {}) or {}
            out[str(rec.get("pmid"))] = {
                "title": re.sub(r"<[^>]+>", "", rec.get("title") or ""),
                "journal": (ji.get("journal") or {}).get("title"),
                "year": rec.get("pubYear"), "doi": rec.get("doi"),
                "url": f"https://europepmc.org/article/MED/{rec.get('pmid')}"}
        except Exception as e:
            print(f"[warn] PMID {p} lookup failed: {e}")
    return out


def main() -> None:
    t = pd.read_csv(D / "gene_table.csv", index_col=0)
    sub_in = json.loads((LIT / "subagent-innate-genes.json").read_text())
    sub_li = json.loads((LIT / "subagent-lineage-artifacts.json").read_text())

    # ---- collect every PMID mentioned by the subagents and verify them in bulk
    mentioned: dict[str, list[str]] = {}
    for rec in sub_in["genes"]:
        blob = " ".join(str(rec.get(k, "")) for k in ("ad_direction", "perturbation_evidence", "uncertainty"))
        mentioned[rec["gene"]] = re.findall(r"PMID\s*:?\s*(\d{6,9})", blob)
    for rec in sub_li["genes"]:
        mentioned[rec["gene"]] = re.findall(r"PMID\s*:?\s*(\d{6,9})", json.dumps(rec))
    all_pm = sorted({p for v in mentioned.values() for p in v})
    verified = verify_pmids(all_pm)
    print(f"subagent PMIDs mentioned={len(all_pm)}  verified={len(verified)}  "
          f"unverified={sorted(set(all_pm) - set(verified))[:8]}")
    (LIT / "subagent-verified-pmids.json").write_text(json.dumps(verified, indent=2, ensure_ascii=False))

    t1_pm = sorted({p for _, p, _ in T1.values() if p} | {p for v in EXTRA_SRC.values() for p in v})
    t1_verified = verify_pmids([p for p in t1_pm if p not in verified])
    t1_verified.update({p: v for p, v in verified.items() if p in t1_pm})
    print(f"tier-1 PMIDs={len(t1_pm)}  verified={len(t1_verified)}")
    (LIT / "tier1-verified-pmids.json").write_text(json.dumps(t1_verified, indent=2, ensure_ascii=False))

    sub_pred = {r["gene"]: r for r in sub_in["genes"]}
    sub_int = {r["gene"]: r for r in sub_li["genes"]}

    rows = []
    for g, r in t.iterrows():
        verdict, pmid, note = T1.get(g, ("", "", ""))
        tier, prov = 1, "本調査で取得・要旨を確認"
        if not verdict:                                    # fall back to the subagent tier
            if g in sub_pred:
                verdict = {"POSITIVE": "POS", "NEGATIVE": "NEG", "MIXED": "MIXED",
                           "no causal evidence (marker)": "NONE"}.get(
                    sub_pred[g]["predicted_shift"].split(" (")[0].strip(), "MIXED")
                note = (sub_pred[g].get("perturbation_evidence") or "")[:220]
                tier, prov = 2, "委託した文献調査（PMID を機械照合）"
                pm = [p for p in mentioned.get(g, []) if p in verified]
                pmid = pm[0] if pm else ""
            elif g in sub_int:
                verdict, tier, prov = "NONE", 2, "委託した文献調査（マーカー判定）"
                note = (sub_int[g].get("interpretation") or "")[:200]
            else:
                verdict, note = "NONE", "文献判定なし"
        match = ""
        if verdict == "POS":
            a, b = r.sign_104M > 0, r.sign_316M > 0
        elif verdict == "NEG":
            a, b = r.sign_104M < 0, r.sign_316M < 0
        else:
            a = b = False
        if verdict in ("POS", "NEG"):
            match = "both" if (a and b) else ("104M" if a else ("316M" if b else "neither"))
        elif verdict == "NONE":
            match = "（判定対象外）"
        else:
            match = "（判定不能）"
        if verdict in ("POS", "NEG") and (r.detection_AD or 0) < 0.005:
            match += "（発現が極小）"
        rows.append({
            "gene": g, "detection_AD": r.detection_AD,
            "shift_104M": r.shift_104M, "rank_104M": int(r.rank_104M),
            "shift_316M": r.shift_316M, "rank_316M": int(r.rank_316M),
            "rank_delta": int(r.rank_delta), "sign_flip": bool(r.sign_flip),
            "candidate_104M": bool(r.candidate_104M), "candidate_316M": bool(r.candidate_316M),
            "lit_verdict": verdict, "lit_note": note, "lit_pmid": str(pmid) if pmid else "", "lit_tier": tier,
            "lit_provenance": prov, "better_supported_model": match,
        })
    adj = pd.DataFrame(rows).sort_values(["sign_flip", "rank_104M"], ascending=[False, True])
    adj.to_csv(D / "adjudication.csv", index=False, float_format="%.6g")
    (D / "adjudication.json").write_text(json.dumps(rows, indent=2, ensure_ascii=False))

    # ---------------- headline numbers
    flips = adj[adj.sign_flip]
    dec = adj[adj.lit_verdict.isin(["POS", "NEG"])]
    dec_exp = dec[dec.detection_AD >= 0.005]          # expressed well enough to carry a signal
    dis = dec_exp[dec_exp.sign_flip]                  # models disagree AND the literature is clear
    print("\n=== headline ===")
    print(f"genes compared: {len(adj)}   sign flips: {len(flips)}")
    print(f"clear literature direction: {len(dec)}  (of which detection>=0.5%: {len(dec_exp)})")
    for lbl, sub in [("all clear-direction genes", dec_exp), ("genes where the models disagree", dis)]:
        if len(sub):
            print(f"  {lbl:34s} n={len(sub):2d}  both={sum(sub.better_supported_model=='both')} "
                  f"only104M={sum(sub.better_supported_model=='104M')} "
                  f"only316M={sum(sub.better_supported_model.str.startswith('316M'))} "
                  f"neither={sum(sub.better_supported_model=='neither')}")
    print("\n=== the decisive set (models disagree, literature clear, gene expressed) ===")
    for _, r in dis.iterrows():
        print(f"  {r.gene:9s} det={r.detection_AD*100:5.1f}%  104M={r.shift_104M:+.5f} (rank {r.rank_104M:2d})  "
              f"316M={r.shift_316M:+.5f} (rank {r.rank_316M:2d})  lit={r.lit_verdict} "
              f"(PMID {r.lit_pmid}) -> {r.better_supported_model}")
    print(f"\n=== where BOTH models contradict the literature ===")
    for _, r in dec_exp[dec_exp.better_supported_model == 'neither'].iterrows():
        print(f"  {r.gene:9s} det={r.detection_AD*100:5.1f}%  104M={r.shift_104M:+.5f}  316M={r.shift_316M:+.5f}  "
              f"lit={r.lit_verdict} (PMID {r.lit_pmid}) — {r.lit_note[:90]}")
    print(f"\n-> {D}/adjudication.csv")


if __name__ == "__main__":
    main()
