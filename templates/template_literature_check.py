#!/usr/bin/env python
"""Literature check for the ISP target genes, via the Europe PMC REST API.

TEMPLATE — a starting point (ひな形), not a finished analysis.
The numbers, literature verdicts and prose baked into this file are a worked
example from the AD_spleen 104M-vs-316M check (2026-09). Re-derive every number
and re-verify every PMID before reusing it on another tissue or model.

For each query it stores the top records (title, journal, year, PMID, DOI,
abstract) under docs/quantization/isp-model-comparison/literature/ and prints
the abstract sentences that carry direction-of-effect information, so the
"does loss of function help or hurt?" verdict in the comparison report can be
traced to a specific paper.

Usage:
    .venv/bin/python templates/template_literature_check.py [--batch 1|2|3|all]
"""
from __future__ import annotations

import argparse
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs/quantization/isp-model-comparison/literature"
API = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"

QUERIES: dict[str, dict[str, str]] = {
    # batch 1 — genes where 104M calls the gene a target candidate and 316M does not
    "ABCA7": {"q": 'TITLE:"ABCA7" AND (loss-of-function OR knockout) AND Alzheimer', "why": "104M +, 316M -"},
    "MEF2C": {"q": 'TITLE:"MEF2C" AND microglia AND (inflammatory OR Alzheimer)', "why": "104M +, 316M -"},
    "PLCG2": {"q": 'TITLE:"PLCG2" AND (protective OR variant) AND Alzheimer', "why": "104M +, 316M -"},
    "SPI1": {"q": 'TITLE:"PU.1" AND Alzheimer AND expression', "why": "104M +, 316M -"},
    "CD68": {"q": 'CD68 AND (Alzheimer OR amyloid) AND (knockout OR deficiency)', "why": "104M +, 316M -"},
    "FGR": {"q": 'Fgr AND (Alzheimer OR amyloid OR microglia)', "why": "104M +, 316M -"},
    "CEACAM1": {"q": 'CEACAM1 AND (Alzheimer OR amyloid OR neuroinflammation)', "why": "104M +, 316M -"},
    "CD3D": {"q": 'CD3D AND (Alzheimer OR amyloid) AND (T cell OR lymphocyte)', "why": "104M +, 316M -"},
    # batch 2 — genes where 316M calls the gene a target candidate and 104M does not
    "TYROBP": {"q": 'TITLE:"Deficiency of TYROBP"', "why": "316M +, 104M -"},
    "CSF1R": {"q": 'TITLE:"CSF1R" AND (inhibitor OR depletion) AND (plaque OR Alzheimer)', "why": "316M +, 104M -"},
    "GCA": {"q": 'grancalcin AND (Alzheimer OR "bone marrow" OR inflammation)', "why": "316M +, 104M -"},
    "LYZ": {"q": '(lysozyme OR Lyzl4) AND amyloid AND (clearance OR microglia)', "why": "316M +, 104M -"},
    "TLR4": {"q": 'TLR4 AND (deficiency OR knockout) AND (amyloid OR Alzheimer)', "why": "316M +, 104M -"},
    # batch 3 — the largest disagreements in magnitude and the genes both models flag
    "APOE": {"q": 'apoE AND (deficiency OR knockout) AND fibrillar amyloid deposition', "why": "both -, literature?"},
    "S100A9": {"q": 'TITLE:"S100A9" AND (knockout OR knockdown) AND (memory OR amyloid)', "why": "both -"},
    "C3": {"q": 'TITLE:"C3" AND complement AND Alzheimer AND (deficien OR knockout)', "why": "both -"},
    "CD33": {"q": 'TITLE:"CD33" AND Alzheimer AND (loss of function OR protective OR knockout)', "why": "both +"},
    "TREM2": {"q": 'TITLE:"TREM2" AND (loss of function OR deficiency) AND amyloid AND seeding', "why": "both ~0/-"},
    "NLRP3": {"q": 'TITLE:"NLRP3" AND (inhibition OR deletion OR knockout) AND Alzheimer', "why": "both +"},
    "SPP1": {"q": 'osteopontin AND (ablation OR knockout) AND (plaque OR 5XFAD)', "why": "316M +, 104M ~0"},
    "CD74": {"q": 'CD74 AND (MIF OR invariant chain) AND (Alzheimer OR amyloid) AND (deficiency OR knockout)', "why": "both -"},
}
BATCH = {"1": list(QUERIES)[:8], "2": list(QUERIES)[8:13], "3": list(QUERIES)[13:]}
KW = re.compile(r"(protect|neuroprotect|reduc|attenuat|ameliorat|improv|worsen|aggravat|increase|decreas|"
                r"loss.of.function|gain.of.function|knockout|knock-out|ablat|deficien|delay|risk|clearance|"
                r"no effect|did not)", re.I)


def fetch(query: str, n: int = 3) -> list[dict]:
    url = f"{API}?query={urllib.parse.quote(query)}&format=json&resultType=core&pageSize={n}"
    with urllib.request.urlopen(url, timeout=40) as r:
        d = json.load(r)
    out = []
    for rec in d.get("resultList", {}).get("result", []):
        ji = rec.get("journalInfo", {}) or {}
        out.append({
            "title": rec.get("title"),
            "journal": (ji.get("journal") or {}).get("title"),
            "year": rec.get("pubYear"), "pmid": rec.get("pmid"), "doi": rec.get("doi"),
            "abstract": re.sub(r"<[^>]+>", "", rec.get("abstractText") or ""),
            "url": f"https://europepmc.org/article/MED/{rec.get('pmid')}" if rec.get("pmid") else None,
        })
    return out


def sentences(abstract: str, k: int = 3) -> list[str]:
    s = [re.sub(r"\s+", " ", x).strip() for x in re.split(r"(?<=[.!?])\s+", abstract)]
    return [x for x in s if 50 < len(x) < 320 and KW.search(x)][:k]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--batch", default="all")
    a = ap.parse_args()
    genes = list(QUERIES) if a.batch == "all" else BATCH[a.batch]
    OUT.mkdir(parents=True, exist_ok=True)
    for g in genes:
        spec = QUERIES[g]
        try:
            recs = fetch(spec["q"])
        except Exception as e:  # network hiccup: keep going
            print(f"##### {g} [{spec['why']}] QUERY FAILED: {e}")
            continue
        (OUT / f"{g}.json").write_text(json.dumps({"gene": g, "query": spec["q"], "records": recs}, indent=2))
        print(f"##### {g}  [{spec['why']}]  ({len(recs)} records)")
        for r in recs[:2]:
            print(f"  - {r['title']} | {r['journal']} {r['year']} | PMID {r['pmid']} | {r['doi']}")
            for s in sentences(r["abstract"], 2):
                print(f"      · {s}")
        print()
        time.sleep(0.4)


if __name__ == "__main__":
    main()
