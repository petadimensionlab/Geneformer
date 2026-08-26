#!/usr/bin/env python
"""Render the PD_spleen tabbed HTML report.

Reads:
  result/report/report_pd_spleen.md
  result/report/data/pd_spleen_early_isp_results.csv
  result/report/data/pd_spleen_gene_list.csv
  result/report/experiments.csv
  result/report/templates/template_pd_spleen.html

Writes:
  result/report/report_pd_spleen.html  (standalone, open in any browser)
"""
from __future__ import annotations

import html
from pathlib import Path

import markdown
import pandas as pd

REPORT = Path(__file__).resolve().parent
DATA = REPORT / "data"
TEMPLATE = REPORT / "templates" / "template_pd_spleen.html"
MD = REPORT / "report_pd_spleen.md"
OUT = REPORT / "report_pd_spleen.html"


def experiments_table() -> str:
    df = pd.read_csv(REPORT / "experiments.csv")
    cols = ["experiment_id", "target_disease", "organ", "script", "n_genes", "max_ncells", "status"]
    rows = []
    for _, r in df.iterrows():
        rows.append(
            "<tr>"
            f"<td><code>{html.escape(str(r['experiment_id']))}</code></td>"
            f"<td>{html.escape(str(r['target_disease']))}</td>"
            f"<td>{html.escape(str(r['organ']))}</td>"
            f"<td><code>{html.escape(str(r['script']))}</code></td>"
            f"<td>{r['n_genes']}</td>"
            f"<td>{r['max_ncells']}</td>"
            f"<td>{html.escape(str(r['status']))}</td></tr>"
        )
    head = "".join(f"<th>{c}</th>" for c in cols)
    return f"<table><tr>{head}</tr>{''.join(rows)}</table>"


def md_section(md: str, header: str) -> str:
    """Return the markdown body under the given ## header, converted to html."""
    lines = md.splitlines()
    out, on = [], False
    for ln in lines:
        if ln.startswith("## ") and header in ln:
            on = True
            continue
        if on and ln.startswith("## "):
            break
        if on:
            out.append(ln)
    return markdown.markdown("\n".join(out), extensions=["tables", "fenced_code"])


def str2float(x) -> float:
    try:
        return float(x)
    except (TypeError, ValueError):
        return 0.0


def main() -> None:
    md_text = MD.read_text(encoding="utf-8")
    cat_df = pd.read_csv(DATA / "pd_spleen_gene_list.csv")
    df = pd.read_csv(DATA / "pd_spleen_early_isp_results.csv")

    df6 = df.sort_values("6m", ascending=False).reset_index(drop=True)
    df9 = df.sort_values("9m", ascending=False).reset_index(drop=True)
    df12 = df.sort_values("12m", ascending=False).reset_index(drop=True)

    # per-timepoint tables (Rank ordering uses that timepoint's column)
    def rank_table(d, col):
        rows = []
        for i, r in enumerate(d.iterrows(), 1):
            _, r = r
            v = str2float(r[col])
            cls = "pos" if v >= 0 else "neg"
            rows.append(
                f"<tr><td>{i}</td><td>{html.escape(str(r['Gene']))}</td>"
                f"<td><code>{html.escape(str(r['Ensembl_ID']))}</code></td>"
                f"<td>{html.escape(str(r['Category']))}</td>"
                f"<td class='pos'>{str2float(r['6m']):+.5f}</td>"
                f"<td class='{'pos' if str2float(r['9m'])>=0 else 'neg'}'>{str2float(r['9m']):+.5f}</td>"
                f"<td class='{'pos' if str2float(r['12m'])>=0 else 'neg'}'>{str2float(r['12m']):+.5f}</td></tr>"
            )
        head = ("<tr><th>Rank</th><th>Gene</th><th>Ensembl_ID</th><th>Category</th>"
                "<th>6m</th><th>9m</th><th>12m</th></tr>")
        return f"<table>{head}{''.join(rows)}</table>"

    # background/rationale from MD sections 2 & 3
    background_html = md_section(md_text, "臓器特異的に探索すべき背景知識")
    if "注目すべき" not in background_html:
        # fallback: sections 2,3,4 combined
        background_html = (md_section(md_text, "臓器特異的に探索すべき背景知識")
                           + md_section(md_text, "細胞・遺伝子の選択理由"))
    rationale_html = md_section(md_text, "細胞・遺伝子の選択理由")
    settings_html = md_section(md_text, "in silico perturbation の主要設定") + md_section(md_text, "計算パラメータ")
    interpretation_html = md_section(md_text, "生物学的解釈")
    results_summary_html = md_section(md_text, "結果まとめ")

    tpl = TEMPLATE.read_text(encoding="utf-8")

    top6 = df6.head(6)["Gene"].tolist()
    top_line = ", ".join(top6)

    ctx = {
        "{{title}}": "PD_spleen — 早期パーキンソン病発症因子の in silico gene deletion レポート",
        "{{subtitle}}": "Geneformer V2-104M 仮想遺伝子削除（delete）· 脾臓（免疫プール）· disease軸 PF→WT（α-syn PFF モデル）",
        "{{organ}}": "脾臓 (PD_spleen)",
        "{{kpi_line}}": "細胞: 時点別免疫プール <b>10,887 / 11,386 / 10,572</b>（PF+WT） · 遺伝子: <b>34</b>候補 → <b>28–29</b>/時点発現 · 実行: <b>約30分</b>",
        "{{top_line}}": top_line,
        "{{rank_img}}": "figures/pd_spleen_isp_rank_6m.png",
        "{{table_6m}}": rank_table(df6, "6m"),
        "{{table_9m}}": rank_table(df9, "9m"),
        "{{table_12m}}": rank_table(df12, "12m"),
        "{{background_html}}": background_html,
        "{{rationale_html}}": rationale_html,
        "{{settings_html}}": settings_html,
        "{{interpretation_html}}": interpretation_html,
        "{{results_summary_html}}": results_summary_html,
        "{{experiments_table}}": experiments_table(),
        "{{date}}": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"),
    }
    for k, v in ctx.items():
        tpl = tpl.replace(k, v)

    OUT.write_text(tpl, encoding="utf-8")
    print("wrote", OUT)


if __name__ == "__main__":
    main()