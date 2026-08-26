#!/usr/bin/env python
"""Render per-organ tabbed HTML reports by injecting content into
result/report/template.html.

Outputs (independent per-organ files):
  result/report/report_ad_smallint.html
  result/report/report_ad_brain.html
  result/report/report_ad_blood.html
  result/report/report_ad_spleen.html

Usage:
    .venv/bin/python result/report/render_html.py
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REPORT = Path(__file__).resolve().parent
DATA = REPORT / "data"


def results_table(df: pd.DataFrame) -> str:
    rows = []
    for _, r in df.iterrows():
        v = float(r["Shift_to_goal_end"])
        cls = "pos" if v >= 0 else "neg"
        rows.append(
            f"<tr><td>{r['Gene']}</td>"
            f"<td><code>{r['Ensembl_ID']}</code></td>"
            f"<td>{r.get('Category', '')}</td>"
            f"<td class='{cls}'>{v:+.5f}</td></tr>"
        )
    head = "<tr><th>Gene</th><th>Ensembl_ID</th><th>カテゴリ</th><th>Shift_to_goal_end</th></tr>"
    return f"<table>{head}{''.join(rows)}</table>"


def experiments_table() -> str:
    df = pd.read_csv(REPORT / "experiments.csv")
    cols = ["experiment_id", "organ", "script", "celltypes", "n_genes", "max_ncells", "status"]
    rows = []
    for _, r in df.iterrows():
        rows.append(
            "<tr>"
            f"<td><code>{r['experiment_id']}</code></td>"
            f"<td>{r['organ']}</td>"
            f"<td><code>{r['script']}</code></td>"
            f"<td>{r['celltypes']}</td>"
            f"<td>{r['n_genes']}</td>"
            f"<td>{r['max_ncells']}</td>"
            f"<td>{r['status']}</td>"
            "</tr>"
        )
    head = "".join(f"<th>{c}</th>" for c in cols)
    return f"<table><tr>{head}</tr>{''.join(rows)}</table>"


def html(p) -> str:
    return f"<p>{p}</p>"


def smallint_context() -> dict:
    return dict(
        title="小腸（AD_smallint）— 早期 AD 発症因子の in silico gene deletion レポート",
        subtitle="Geneformer V2-104M 仮想遺伝子削除（delete）· 小腸（腸-脳軸）· disease軸 AD→WT",
        organ="小腸 (AD_smallint)",
        dataset="AD_smallint",
        kpi_line="細胞: <b>17,949</b>（早期プール） · 遺伝子: <b>49</b> · 実行: <b>45分</b>",
        top_line="<b>CLDN1</b> (+0.0093), MUC5B (+0.0038), TREM2 (+0.0029), SORL1 (+0.0023), MUC2 (+0.0018)",
        rank_img="figures/ad_smallint_isp_rank.png",
        background_html=(
            "<p>AD は脳疾患だが、<b>腸-脳軸（gut-brain axis）</b>が早期病態に関与。小腸は病原体・内毒素の第一バリアで、透過性亢進→全身性炎症が示唆される（Vogt 2017; Kimura 2021; Cryan 2019）。</p>"
            "<ul>"
            "<li><b>注目細胞</b>: 腸上皮/バリア（IECs, Paneth, Goblet, Epithelial）, マクロファージ, T細胞, NK</li>"
            "<li><b>注目遺伝子クラス</b>: 粘液（MUC2/5B）, タイトジャンクション（CLDN/OCLN/TJP）, 抗菌ペプチド, 免疫・炎症</li>"
            "</ul>"
        ),
        rationale_html=(
            "<ul>"
            "<li>免疫（Macrophages, CD4/CD8 T, NK）+ バリア（IECs, Paneth, Goblet, Epithelial）8種をプール</li>"
            "<li>49遺伝子（仮説駆動）: AD risk（APOE/TREM2/CLU/BIN1/CD33/ABCA7/SORL1/TYROBP/INPP5D/PLCG2/MEF2C）、炎症、粘液、TJ、抗菌、幹細胞</li>"
            "<li>文献: Lambert 2013; Wightman 2021（GWAS）、Bergstrom 2017（粘液）、Ulluwishewa 2011（TJ）</li>"
            "</ul>"
        ),
        settings_html=(
            "<p><b>delete</b> · combos=0 · disease（AD→WT） · V2-104M · emb=cls</p>"
            "<table><tr><th>パラメータ</th><th>値</th></tr>"
            "<tr><td>解析細胞（早期）</td><td>17,949</td></tr>"
            "<tr><td>max_ncells</td><td>300</td></tr>"
            "<tr><td>遺伝子数</td><td>49（発現）</td></tr>"
            "<tr><td>forward_batch_size</td><td>64</td></tr>"
            "<tr><td>実行時間</td><td>45分</td></tr></table>"
        ),
        interpretation_html=(
            "<p><b>CLDN1</b>（タイトジャンクション）が最上位（+0.0093）。粘液バリア（MUC5B, MUC2, MUC4）と免疫/ミクログリア遺伝子（TREM2, TYROBP, C1QA, TNF）も正で、<b>早期 AD 小腸のバリア整合性と免疫応答が重要標的</b>。</p>"
        ),
    )


def brain_context() -> dict:
    return dict(
        title="脳（AD_brain）— 早期 AD 発症因子の in silico gene deletion レポート",
        subtitle="Geneformer V2-104M 仮想遺伝子削除（delete）· 脳（ミクログリア/BAM）· disease軸 AD→WT",
        organ="脳 (AD_brain)",
        dataset="AD_brain",
        kpi_line="細胞: ミクログリア 早期 AD/WT 1,800/1,800, BAM 376/264 · 遺伝子: <b>23</b> · 実行: <b>67分</b>",
        top_line="<b>APOE</b> (+0.0087), TREM2 (+0.00006)",
        rank_img="figures/ad_brain_isp_rank.png",
        background_html=(
            "<p>AD は<b>神経炎症</b>が中心。ミクログリアは Aβ 蓄積のごく初期から <b>DAM</b> へ分化（Keren-Shaul 2017）。TREM2 は AD リスク最有力（Guerreiro 2013）。</p>"
            "<ul>"
            "<li><b>注目細胞</b>: ミクログリア（Microglia）, BAM（border-associated macrophages）</li>"
            "<li><b>注目遺伝子クラス</b>: DAM シグネチャー（TREM2/APOE/TYROBP/CSF1R/C1QA/B/C3/ITGAM/ITGAX/CD68/SPP1）、炎症（IL1B/TNF）、AD GWAS リスク</li>"
            "</ul>"
        ),
        rationale_html=(
            "<ul>"
            "<li>ミクログリア + BAM（脳常在免疫）を早期時点で選択 — DAM 遷移の中核（Keren-Shaul 2017）</li>"
            "<li>23遺伝子: AD risk（APOE/TREM2/TYROBP/BIN1/CD33/ABCA7/SORL1/INPP5D/PLCG2/MEF2C）+ DAM（CSF1R/C1QA/B/C3/CD68/AIF1/ITGAM/ITGAX/SPP1/CCL2）+ 炎症（IL1B/TNF/CD74）</li>"
            "<li>文献: Guerreiro 2013（TREM2 リスク）; Keren-Shaul 2017（DAM）; Kinney 2018（炎症仮説）</li>"
            "</ul>"
        ),
        settings_html=(
            "<p>細胞: <b>ミクログリア + BAM</b>（早期） · 23遺伝子（AD risk + DAM + 炎症） · delete · V2-104M · max_ncells=200</p>"
            "<table><tr><th>パラメータ</th><th>値</th></tr>"
            "<tr><td>max_ncells</td><td>200</td></tr>"
            "<tr><td>遺伝子数</td><td>23（全て発現）</td></tr>"
            "<tr><td>forward_batch_size</td><td>64</td></tr>"
            "<tr><td>実行時間</td><td>67分</td></tr></table>"
        ),
        interpretation_html=(
            "<p><b>APOE が突出（+0.0087）</b>。早期ミクログリア/BAM からの APOE 削除が embedding を早期 WT に最も戻す＝<b>早期発見・介入の最有力標的</b>。他の21遺伝子は ~0（max_ncells=200 の限界）。</p>"
        ),
    )


def blood_context() -> dict:
    return dict(
        title="血液（AD_blood）— 早期 AD 発症因子の in silico gene deletion レポート",
        subtitle="Geneformer V2-104M（ファインチューン済みセル分類器）仮想遺伝子削除（delete）· 血液（末梢免疫）· disease軸 AD→WT",
        organ="血液 (AD_blood)",
        dataset="AD_blood",
        kpi_line="細胞: <b>48,909</b>（うち免疫プール 10 細胞型） · 遺伝子: <b>55</b>（全時点 273 結果） · 実行: <b>約1時間35分</b>",
        top_line="<b>CD3E</b> (+0.0142), CD4 (+0.0138), CD74 (+0.0053), S100A8 (+0.0050), C1QC (+0.0049), CD19 (+0.0041)",
        rank_img="figures/ad_blood_isp_rank_3m.png",
        extra_figures=(
            "<figure><img src='figures/ad_blood_isp_timeseries.png' alt='timeseries'>"
            "<figcaption>Shift_to_goal_end の時点別推移（3m→12m, top10@3m 強調）</figcaption></figure>"
            "<figure><img src='figures/ad_blood_isp_heatmap.png' alt='heatmap'>"
            "<figcaption>遺伝子 × 時点の Shift_to_goal_end 行列</figcaption></figure>"
        ),
        results_table2=_blood_tp_table(),
        background_html=(
            "<p>AD は脳疾患だが、<b>血液（末梢免疫）は侵襲性の低い早期発見バイオマーカーの宝庫</b>。臨床前段階から末梢免疫プロファイルが変化する（Gate 2020; Tan 2002）。</p>"
            "<ul>"
            "<li><b>注目細胞</b>: 単球（Ly6c.low/high Monocytes）= ミクログリア様血中細胞, T 細胞（CD4/CD8）, B 細胞, 好中球, DC（c/p）</li>"
            "<li><b>注目遺伝子クラス</b>: AD GWAS リスク（APOE/TREM2/TYROBP…）, T 細胞/適応免疫（CD3/CD4/CD8/NKG7/CD19）, 炎症/警報因子（S100A8/A9, CD74, C1Q, IL1B/TNF）, 骨髄球/顆粒球/血小板</li>"
            "</ul>"
        ),
        rationale_html=(
            "<ul>"
            "<li>血液免疫 10 種を選択（低シグナル: Basophils/Megakaryocytes/Erythroblasts/Immature.T/NKT を除外）</li>"
            "<li>55 遺伝子: AD GWAS 12 + 炎症/自然免疫 24 + 適応免疫 12 + 骨髄球/血小板 10（全時点の AD 細胞で発現確認）</li>"
            "<li>文献: Gate 2020（CD8 T クローン, Nature）; Tan 2002（S100A8/9）; Nakamura 2018（血液 AD 検出）; Guerreiro 2013（TREM2）; Jansen 2019（GWAS 免疫濃縮）</li>"
            "</ul>"
        ),
        settings_html=(
            "<p><b>delete</b> · combos=0 · disease（AD→WT） · ファインチューン済み CellClassifier（V2-104M） · emb=cls · 各時点個別評価（3m/4p5m/6m/9m/12m）</p>"
            "<table><tr><th>パラメータ</th><th>値</th></tr>"
            "<tr><td>免疫プール細胞型</td><td>10 種（CD4/CD8 T, NK_ILC1, B, 好中球, 単球×2, cDC, pDC, γδT）</td></tr>"
            "<tr><td>max_ncells</td><td>200</td></tr>"
            "<tr><td>遺伝子数</td><td>55（時点別発現: 55/54/55/55/54）</td></tr>"
            "<tr><td>forward_batch_size</td><td>64</td></tr>"
            "<tr><td>nproc / datasets</td><td>1 / ==4.0.0</td></tr>"
            "<tr><td>実行時間</td><td>約1時間35分（CUDA）</td></tr></table>"
        ),
        interpretation_html=(
            "<p><b>T 細胞軸（CD3E +0.0142, CD4 +0.0138）が最有力</b>。最早期（3m）の AD 血液で T 細胞活性化/共刺激シグネチャーの削除が WT 方向へ最大のシフト＝<b>早期発見・介入の最有力標的</b>。CD74（MHC-II）・S100A8（警報因子）も全時点で正。</p>"
            "<p>血液では <b>APOE はほぼ中性</b>（脳で最有力だった 07_ad_spleen_early_isp と対照的、ミクログリア不在のため）。<b>CD8A / MS4A1 / NKG7 は強く負</b>＝削除すると WT から離れる（CD8 細胞毒性・B 細胞・NK 顆粒は健常状態の維持に寄与）。</p>"
        ),
    )


def _blood_tp_table() -> str:
    df = pd.read_csv(DATA / "ad_blood_early_isp_results.csv")
    tp3 = df[df["Timepoint"] == "3m"].sort_values("Shift_to_goal_end", ascending=False).head(12)
    rows = []
    for _, r in tp3.iterrows():
        cells = []
        for tt in ["4p5m", "6m", "9m", "12m"]:
            v = df[(df["Gene"] == r["Gene"]) & (df["Timepoint"] == tt)]["Shift_to_goal_end"]
            val = float(v.iloc[0]) if len(v) else float("nan")
            cls = "pos" if val >= 0 else "neg"
            cells.append(f"<td class='{cls}'>{val:+.5f}</td>")
        cls3 = "pos" if float(r["Shift_to_goal_end"]) >= 0 else "neg"
        rows.append(
            f"<tr><td>{r['Gene']}</td><td>{r['Category']}</td>"
            f"<td class='{cls3}'>{r['Shift_to_goal_end']:+.5f}</td>{''.join(cells)}</tr>"
        )
    head = ("<tr><th>Gene</th><th>カテゴリ</th><th>3m</th>"
            "<th>4.5m</th><th>6m</th><th>9m</th><th>12m</th></tr>")
    return (
        "<h3>上位12遺伝子の時点別 Shift_to_goal_end</h3>"
        f"<table>{head}{''.join(rows)}</table>"
    )


def _spleen_tp_table() -> str:
    df = pd.read_csv(DATA / "ad_spleen_early_isp_results.csv")
    tp3 = df[df["Timepoint"] == "3m"].sort_values("Shift_to_goal_end", ascending=False).head(12)
    rows = []
    for _, r in tp3.iterrows():
        cells = []
        for tt in ["4p5m", "6m"]:
            v = df[(df["Gene"] == r["Gene"]) & (df["Timepoint"] == tt)]["Shift_to_goal_end"]
            val = float(v.iloc[0]) if len(v) else float("nan")
            cls = "pos" if val >= 0 else "neg"
            cells.append(f"<td class='{cls}'>{val:+.5f}</td>")
        cls3 = "pos" if float(r["Shift_to_goal_end"]) >= 0 else "neg"
        rows.append(
            f"<tr><td>{r['Gene']}</td><td>{r['Category']}</td>"
            f"<td class='{cls3}'>{r['Shift_to_goal_end']:+.5f}</td>{''.join(cells)}</tr>"
        )
    head = ("<tr><th>Gene</th><th>カテゴリ</th><th>3m</th><th>4.5m</th><th>6m</th></tr>")
    return (
        "<h3>上位12遺伝子の時点別 Shift_to_goal_end</h3>"
        f"<table>{head}{''.join(rows)}</table>"
    )


def spleen_context() -> dict:
    return dict(
        title="脾臓（AD_spleen）— 早期 AD 発症因子の in silico gene deletion レポート",
        subtitle="Geneformer V2-104M（ファインチューン済みセル分類器）仮想遺伝子削除（delete）· 脾臓（二次リンパ器官・免疫）· disease軸 AD→WT",
        organ="脾臓 (AD_spleen)",
        dataset="AD_spleen",
        kpi_line="細胞: <b>76,947</b>（うち免疫プール 19 種類, 62,762 / 早期 37,425） · 遺伝子: <b>55</b>（全時点 165 結果） · 実行: <b>約50分</b>",
        top_line="<b>CLEC9A</b> (+0.0066), ITGA2B (+0.0062), CD19 (+0.0034), CD4 (+0.0034), CLU (+0.0029), CD33 (+0.0014)",
        rank_img="figures/ad_spleen_isp_rank_3m.png",
        extra_figures=(
            "<figure><img src='figures/ad_spleen_isp_timeseries.png' alt='timeseries'>"
            "<figcaption>Shift_to_goal_end の時点別推移（3m→6m, top12@3m 強調）</figcaption></figure>"
            "<figure><img src='figures/ad_spleen_isp_heatmap.png' alt='heatmap'>"
            "<figcaption>遺伝子 × 時点の Shift_to_goal_end 行列</figcaption></figure>"
        ),
        results_table2=_spleen_tp_table(),
        background_html=(
            "<p>脾臓は<b>最大の二次リンパ器官</b>で、血中免疫細胞の貯蔵・再循環と単球/マクロファージ系・T/B/DC の成熟・活性化の中枢。脳（ミクログリア）から離れた末梢にあるが、<b>脾臓-脳軸（spleen-brain axis）</b>を介して神経炎症と双方向に連関する。</p>"
            "<ul>"
            "<li><b>注目細胞</b>: 単球/マクロファージ系（Ly6c.high/low・Macrophages）, DC / cDC1（CLEC9A）, T 細胞（CD4/CD8）, B 細胞（MZ・Naive/Memory・Plasma）, 好中球, NK/ILC1</li>"
            "<li><b>注目遺伝子クラス</b>: AD GWAS リスク（APOE/TREM2/TYROBP/CD33…）, cDC1 クロス提示（CLEC9A/ITGAM/ITGAX）, 適応免疫（CD3/CD4/CD8/B）, 補体/炎症（C1Q/S100A8/A9/IL1B/TNF）, 血小板系（ITGA2B/VWF）</li>"
            "</ul>"
        ),
        rationale_html=(
            "<ul>"
            "<li>脾臓免疫 19 種を選択（骨髄系=単球/マクロファージ/DC/cDC1 を優先、ストロマ・血管・赤血球/巨核球系ノイズ: Erythroblasts/Megakaryocytes/Basophils/Fibroblasts/endothelial を除外）</li>"
            "<li>55 遺伝子は<b>血中 07_ad_spleen_early_isp と同一パネル</b>（臓器横断比較が可能）。全遺伝子が 3m AD 脾臓免疫プールで発現（CX3CR1 のみ脾臓では不発現のため除外）</li>"
            "<li>文献: Gate 2020（CD8 T クローン, Nature）; Keren-Shaul 2017（DAM, TREM2/APOE/TYROBP）; Guerreiro 2013（TREM2）; Nakamura 2018（血液 AD 検出）; Jansen 2019（GWAS 免疫濃縮）</li>"
            "</ul>"
        ),
        settings_html=(
            "<p><b>delete</b> · combos=0 · disease（AD→WT） · ファインチューン済み CellClassifier（V2-104M, 25クラス, acc 0.93） · emb=cls · 各時点個別評価（3m/4p5m/6m）</p>"
            "<table><tr><th>パラメータ</th><th>値</th></tr>"
            "<tr><td>免疫プール細胞型</td><td>19 種（骨髄系 + リンパ系、ストロマ/血球系ノイズ除外）</td></tr>"
            "<tr><td>max_ncells</td><td>200（perturb）/ 1000（状態埋め込み）</td></tr>"
            "<tr><td>遺伝子数</td><td>55（各時点で AD 細胞に発現: 3m=55, 4.5m=55, 6m=55）</td></tr>"
            "<tr><td>forward_batch_size</td><td>64</td></tr>"
            "<tr><td>nproc / datasets</td><td>1 / ==4.0.0</td></tr>"
            "<tr><td>実行時間</td><td>約50分（CUDA）</td></tr></table>"
        ),
        interpretation_html=(
            "<p><b>cDC1 軸: CLEC9A (+0.0066)</b> が最有力。cDC1 特異マーカー（DNGR1）でクロスプレゼンテーション・細胞死取込みに必須。最早期 AD 脾臓での cDC1 削除が WT 方向へ最大シフト＝<b>脾臓 cDC1 が早期の疾病状態を駆動</b>。</p>"
            "<p><b>血小板系 ITGA2B (+0.0062, 4.5m +0.0141)</b> と <b>B 細胞 CD19 (+0.0034, 全時点で正)</b> も安定的に上位（血中 07_ad_spleen_early_isp で一貫）。AD リスク遺伝子では <b>CLU/CD33/SORL1/TREM2 が正</b>、脾臓では<b> APOE は中性</b>（脳 07_ad_spleen_early_isp の APOE +0.0087 と対照）。</p>"
            "<p><b>CD8A (-0.098) は極端に負</b>＝CD8 細胞毒性シグネチャーは WT 状態の維持に寄与。LYZ/S100A8/S100A9 も負（血中の S100A8 正と逆符号・臓器差）。4.5m で CD4 符号反転（+0.0034→-0.0082）＝免疫状態の短期再構成を示唆。</p>"
        ),
    )


def render(tpl: str, ctx: dict, df: pd.DataFrame) -> str:
    rep = {
        "{{title}}": ctx["title"],
        "{{subtitle}}": ctx["subtitle"],
        "{{organ}}": ctx["organ"],
        "{{dataset}}": ctx["dataset"],
        "{{kpi_line}}": ctx["kpi_line"],
        "{{top_line}}": ctx["top_line"],
        "{{rank_img}}": ctx["rank_img"],
        "{{background_html}}": ctx["background_html"],
        "{{rationale_html}}": ctx["rationale_html"],
        "{{settings_html}}": ctx["settings_html"],
        "{{interpretation_html}}": ctx["interpretation_html"],
        "{{results_table}}": results_table(df),
        "{{results_table2}}": ctx.get("results_table2", ""),
        "{{extra_figures}}": ctx.get("extra_figures", ""),
        "{{experiments_table}}": experiments_table(),
        "{{date}}": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M"),
    }
    for k, v in rep.items():
        tpl = tpl.replace(k, v)
    return tpl


def main() -> None:
    tpl = (REPORT / "template.html").read_text(encoding="utf-8")

    sm = pd.read_csv(DATA / "ad_smallint_early_isp_results.csv").merge(
        pd.read_csv(DATA / "ad_smallint_gene_list.csv")[["Gene", "Category"]],
        on="Gene",
        how="left",
    ).sort_values("Shift_to_goal_end", ascending=False)
    (REPORT / "report_ad_smallint.html").write_text(
        render(tpl, smallint_context(), sm), encoding="utf-8"
    )
    print("wrote", REPORT / "report_ad_smallint.html")

    br = pd.read_csv(DATA / "ad_brain_early_isp_results.csv").merge(
        pd.read_csv(DATA / "ad_brain_gene_list.csv")[["Gene", "Category"]],
        on="Gene",
        how="left",
    ).sort_values("Shift_to_goal_end", ascending=False)
    (REPORT / "report_ad_brain.html").write_text(
        render(tpl, brain_context(), br), encoding="utf-8"
    )
    print("wrote", REPORT / "report_ad_brain.html")

    bl = pd.read_csv(DATA / "ad_blood_early_isp_results.csv").merge(
        pd.read_csv(DATA / "ad_blood_gene_list.csv")[["Gene", "Category"]],
        on="Gene",
        how="left",
    )
    bl3 = bl[bl["Timepoint"] == "3m"].sort_values("Shift_to_goal_end", ascending=False)
    (REPORT / "report_ad_blood.html").write_text(
        render(tpl, blood_context(), bl3), encoding="utf-8"
    )
    print("wrote", REPORT / "report_ad_blood.html")

    sp = pd.read_csv(DATA / "ad_spleen_early_isp_results.csv")
    sp3 = sp[sp["Timepoint"] == "3m"].sort_values("Shift_to_goal_end", ascending=False)
    (REPORT / "report_ad_spleen.html").write_text(
        render(tpl, spleen_context(), sp3), encoding="utf-8"
    )
    print("wrote", REPORT / "report_ad_spleen.html")


if __name__ == "__main__":
    main()