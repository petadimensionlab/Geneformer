#!/usr/bin/env python3
"""Build the tabbed HTML report from the AD_brain ISP markdown + CSV results.

Reads:
  report/ISP_report_AD_brain.md
  report/isp_results_consolidated.csv
  report/isp_results_early_ad_brain.csv
  report/isp_results_early_ad_microglia.csv
  report/experiment_manifest.csv
  report/templates/report_template.html
  report/figures/*.png

Writes:
  report/ISP_report_AD_brain.html  (standalone, open in any browser)
"""
from __future__ import annotations

import csv
import html
import re
from pathlib import Path

import markdown
from jinja2 import Template

ROOT = Path(__file__).resolve().parent.parent / "input" / "AD_brain"
REPORT = ROOT / "results" / "report"
FIG = REPORT / "figures"
TEMPLATE = REPORT / "templates" / "report_template.html"

MD = REPORT / "ISP_report_AD_brain.md"


def md_to_html(text: str) -> str:
    return markdown.markdown(text, extensions=["tables", "fenced_code", "toc"])


def split_sections(md_text: str) -> tuple[dict[str, list[str]], list[str]]:
    """Split the markdown into top-level `## N. Title` sections."""
    lines = md_text.splitlines()
    sections: dict[str, list[str]] = {}
    order: list[str] = []
    cur = None
    for ln in lines:
        m = re.match(r"^##\s+(\d+)(?:\\.|\.)\s*(.*)$", ln)
        if m:
            cur = m.group(1)
            order.append(cur)
            sections[cur] = [f"## {m.group(1)}. {m.group(2)}"]
        elif cur is not None:
            sections[cur].append(ln)
    return {k: "\n".join(v) for k, v in sections.items()}, order


def csv_to_html(path: Path, color_key: str | None = None) -> str:
    with open(path) as f:
        rows = list(csv.reader(f))
    if not rows:
        return "<p>no data</p>"
    head = rows[0]
    out = ["<table><thead><tr>"]
    out += [f"<th>{html.escape(c)}</th>" for c in head]
    out.append("</tr></thead><tbody>")
    for r in rows[1:]:
        out.append("<tr>")
        for i, cell in enumerate(r):
            cls = ""
            if color_key and head[i] == color_key:
                try:
                    v = float(cell)
                    cls = " pos" if v > 1e-5 else (" neg" if v < -1e-5 else " zero")
                except ValueError:
                    pass
            out.append(f'<td class="{cls.strip()}">{html.escape(cell)}</td>')
        out.append("</tr>")
    out.append("</tbody></table>")
    return "".join(out)


def main() -> None:
    md_text = MD.read_text()
    md = markdown.Markdown(extensions=["tables", "fenced_code", "toc"])
    full_html = md.convert(md_text)
    sections, order = split_sections(md_text)

    # overview: keep sections 1-4 and 7 (context), conclusions: section 6 core
    keep_ids = ["1", "2", "3", "4", "7"]
    overview_blocks = []
    for sid in keep_ids:
        if sid in sections:
            overview_blocks.append(md_to_html(sections[sid]))

    conclusions_html = md_to_html(sections.get("6", ""))

    results_html = csv_to_html(REPORT / "isp_results_consolidated.csv",
                               color_key="Shift_Microglia_only")
    results_html += "<h3>E1: Microglia + BAM</h3>"
    results_html += csv_to_html(REPORT / "isp_results_early_ad_brain.csv",
                                color_key="Shift_to_goal_end")
    results_html += "<h3>E2: Microglia only</h3>"
    results_html += csv_to_html(REPORT / "isp_results_early_ad_microglia.csv",
                                color_key="Shift_to_goal_end")

    manifest_html = csv_to_html(REPORT / "experiment_manifest.csv")

    figures = [
        {"path": "figures/isp_shift_by_gene.png",
         "caption": "遺伝子別 Shift_to_goal_end（左: Microglia+BAM, 右: Microglia 単独）"},
        {"path": "figures/isp_shift_comparison.png",
         "caption": "細胞プール間の shift 比較（E1 vs E2）"},
        {"path": "figures/isp_shift_top_bottom.png",
         "caption": "Microglia 単独でのランキング（正 = AD→WT 方向）"},
    ]

    tpl = Template(TEMPLATE.read_text())
    page = tpl.render(
        page_title="Geneformer in silico perturbation — AD brain 早期検出レポート",
        page_subtitle="AD (brain) · ミクログリア/骨髄系の仮想的遺伝子欠失が早期 WT 方向へ遷移させるか",
        kpi_disease="AD", kpi_organ="brain", kpi_genes="23", kpi_cells="1,000",
        overview_blocks=overview_blocks,
        conclusions=conclusions_html,
        results_html=results_html,
        figures=figures,
        manifest_html=manifest_html,
        full_md=full_html,
        generated_at="2026-08-25",
    )
    out = REPORT / "ISP_report_AD_brain.html"
    out.write_text(page)
    print("wrote", out)


if __name__ == "__main__":
    main()