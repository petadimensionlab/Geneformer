#!/usr/bin/env python
"""Build an interactive, tabbed HTML report for the AD_BM ISA analysis.

Merges the three user-custom tabs:
  1. lab report   — rendered markdown (report.md -> report_template.html)
  2. figures      —— generated ISP images
  3. experiments  ——  csv   manifest table (with +/- colouring)

Usage:
    ADPD_ROOT=input/AD_BM .venv/bin/python analysis/08h_make_html_report_AD_BM.py
"""
from __future__ import annotations

import os
import re
from pathlib import Path

import markdown
import pandas as pd
from jinja2 import Template

ROOT = Path(os.environ.get("ADPD_ROOT", "input/AD_BM"))
REPORT = ROOT / "results" / "report"
MD = REPORT / "report.md"
TPL = REPORT / "report_template.html"
OUT = REPORT / "report.html"
CSV = REPORT / "experiments_mapping.csv"

FIG_G1 = "figures/isp_goal_state_shift_by_celltype.png"
FIG_G2 = "figures/isp_goal_state_shift_heatmap.png"
FIG_G3 = "figures/isp_goal_state_shift_bar.png"

MD_RENDER = markdown.Markdown(extensions=["tables", "fenced_code", "attr_list", "nl2br"])


def split_sections(md_text: str) -> tuple[str, list[dict]]:
    lines = md_text.splitlines()
    title, cur_label, cur = "", None, []
    sections = []
    for ln in lines:
        if ln.startswith("# "):
            title = ln[2:].strip()
        elif ln.startswith("## "):
            if cur_label is not None:
                sections.append({"label": cur_label, "body": "\n".join(cur)})
            cur_label, cur = ln[3:].strip(), []
        elif ln.startswith("> **⚠"):
            # keep corruption callouts with the current section, but hoist them
            sections[-1]["callout"] = ln
        else:
            cur.append(ln)
    if cur_label is not None:
        sections.append({"label": cur_label, "body": "\n".join(cur)})
    return title, sections


def figures_tab() -> str:
    return f"""
    <div class="figures">
      <figure><img src="{FIG_G1}" alt="grouped"><figcaption>図1 ｜ 各遺伝子 × 各celltype の goal_state_shift（deletion、AD→WT）。緑=正（WT方向）、赤=負。</figcaption></figure>
      <figure><img src="{FIG_G2}" alt="heatmap"><figcaption>図2 ｜ gene × celltype のヒートマップ。APOE（Macrophage）が突出した正。</figcaption></figure>
      <figure><img src="{FIG_G3}" alt="ranked"><figcaption>図3 ｜ 全22ランを shift 順に並べた棒グラフ（正=WT 方向へ）。</figcaption></figure>
    </div>"""


def csv_tab() -> str:
    if not CSV.exists():
        return "<p>experiments_mapping.csv がありません。</p>"
    df = pd.read_csv(CSV)
    rows = []

    def fmt(v):
        return f"{float(v):+.6f}"

    for r in df.to_dict("records"):
        cls = "pos" if float(r["Shift_to_goal_end"]) > 0 else "neg"
        rows.append(
            "<tr>"
            f"<td>{r['run_script']}</td>"
            f"<td>{r['celltype']}</td>"
            f"<td>{r['start_state']}&rarr;{r['goal_state']}</td>"
            f"<td>{r['gene']}<br><code>{r['ensembl']}</code></td>"
            f"<td>{int(r['expr_pct_AD_cells'])}</td>"
            f"<td>{int(r['n_perturbed_cells_AD'])}</td>"
            f"<td class='{cls}'>{fmt(r['Shift_to_goal_end'])}</td>"
            f"<td><code>{Path(r['csv']).name}</code></td>"
            "</tr>"
        )
    col = ["script", "celltype", "state", "gene<br>Ensembl", "AD expr%", "AD used", "adj shift", "CSV"]
    table = "".join(f"<th>{c}</th>" for c in col)
    return f"""
    <p>各実験（run）と出力ファイル・結果の対応。 <a href='experiments_mapping.csv' download>CSV をダウンロード</a></p>
    <div class='table-wrap'><table>
      <thead><tr>{table}</tr></thead>
      <tbody>{''.join(rows)}</tbody>
    </table></div>"""


def main() -> None:
    md_text = MD.read_text(encoding="utf-8")
    title, sections = split_sections(md_text)
    tabs = [{"label": s["label"], "html": MD_RENDER.convert(s["body"])}
            for s in sections]
    for t in tabs:
        MD_RENDER.reset()

    tabs.append({"label": "Figures", "html": figures_tab()})
    tabs.append({"label": "Experiment CSV", "html": csv_tab()})

    subtitle = ("AD（Alzheimer's disease）／ 骨髄（AD_BM，mouse） ／ 2 cells × 11 genes = 22 runs ／ "
                "delete / goal_state_shift（AD → WT） ／ 実行: analysis/07f_in_silico_perturbation_AD_BM.py "
                "・可視化: analysis/08g_isp_visualization_AD_BM.py")

    tpl = Template(TPL.read_text(encoding="utf-8"))
    OUT.write_text(tpl.render(title=title, subtitle=subtitle, tabs=tabs), encoding="utf-8")
    print("wrote", OUT, f"(~{len(OUT.read_text())}B, {len(tabs)} tabs)")


main()