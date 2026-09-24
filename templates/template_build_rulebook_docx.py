#!/usr/bin/env python
"""テンプレート: ISP ルールブック（docs/isp の規範ページ）を1冊の文書にまとめます。

このファイルはテンプレートです。自分のリポジトリに合わせて SOURCES を書き換えてください。

使い方:
    python templates/template_build_rulebook_docx.py
    soffice --headless --convert-to pdf rulebook/RULEBOOK-ISP-ja.docx

出力:
    rulebook/RULEBOOK-ISP-ja.md    （原稿。差分を取るにはこちら）
    rulebook/RULEBOOK-ISP-ja.docx  （Word。以降 PDF に変換）
    rulebook/RULEBOOK-ISP-ja.pdf   （soffice で生成）

対応する Markdown: 見出し（#〜####）、表（パイプ）、箇条書き（- と入れ子）、
チェックボックス（- [ ]）、番号付きリスト、コードブロック（```）、引用（>）、
インラインの **強調** と `コード`、水平線（---）。
"""
from __future__ import annotations

import re
import subprocess
from datetime import date
from pathlib import Path

from docx import Document
from docx.enum.text import WD_BREAK
from docx.shared import Pt, RGBColor

ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "rulebook"

# まとめる順番（規範 → 基準 → 遵守項目 → 定義 → 手順 → 記録）
SOURCES = [
    ("docs/isp/README.md", "はじめに（索引）"),
    ("docs/isp/style-guide.md", "記述ルール（唯一の基準）"),
    ("docs/isp/criteria.md", "評価基準（E1〜E6）"),
    ("docs/isp/checklist.md", "遵守項目（チェックリスト）"),
    ("docs/isp/methods.md", "方法・共通設定（用語定義と null 分布の定義）"),
    ("docs/isp/glossary.md", "用語と手法の定義"),
    ("docs/isp/design.md", "図表・スライド・ポスターの視覚設計ルール"),
    ("docs/isp/report_template.md", "レポートの標準構成"),
    ("docs/isp/owners.md", "執筆者名簿とページ担当"),
    ("docs/isp/requests/README.md", "受信箱の説明"),
    ("docs/isp/conflicts.md", "競合と裁定の記録"),
    ("docs/isp/how_to_run.md", "付録: 再現手順"),
]

BOLD = re.compile(r"\*\*(.+?)\*\*")
CODE = re.compile(r"`([^`]+)`")
LINK = re.compile(r"\[([^\]]+)\]\(([^)]+)\)")
BULLET = re.compile(r"^(\s*)[-*]\s+(.*)$")
CHECKBOX = re.compile(r"^\[( |x)\]\s*(.*)$")
NUMBERED = re.compile(r"^(\s*)(\d+)\.\s+(.*)$")


def git_commit() -> str:
    try:
        out = subprocess.run(["git", "-C", str(ROOT), "log", "-1", "--format=%h %ad", "--date=short"],
                             capture_output=True, text=True, check=True)
        return out.stdout.strip()
    except Exception:
        return "unknown"


def add_runs(par, text: str, size: float | None = None, mono: bool = False) -> None:
    """Write text into a paragraph, honouring **bold**, `code` and [link](url)."""
    pos = 0
    tokens: list[tuple[int, int, str]] = []
    for rx, kind in ((BOLD, "bold"), (CODE, "code"), (LINK, "link")):
        for m in rx.finditer(text):
            tokens.append((m.start(), m.end(), kind))
    tokens.sort()
    for start, end, kind in tokens:
        if start < pos:
            continue
        if start > pos:
            r = par.add_run(text[pos:start])
            if size:
                r.font.size = Pt(size)
            if mono:
                r.font.name = "Consolas"
        if kind == "bold":
            r = par.add_run(BOLD.match(text, start).group(1))
            r.bold = True
        elif kind == "code":
            r = par.add_run(CODE.match(text, start).group(1))
            r.font.name = "Consolas"
        else:
            m = LINK.match(text, start)
            r = par.add_run(f"{m.group(1)}（{m.group(2)}）")
        if size:
            r.font.size = Pt(size)
        pos = end
    if pos < len(text):
        r = par.add_run(text[pos:])
        if size:
            r.font.size = Pt(size)
        if mono:
            r.font.name = "Consolas"


def render_table(doc: Document, rows: list[list[str]]) -> None:
    cols = max(len(r) for r in rows)
    t = doc.add_table(rows=0, cols=cols)
    t.style = "Light Grid Accent 1"
    for row in rows:
        cells = t.add_row().cells
        for i in range(cols):
            text = row[i] if i < len(row) else ""
            cells[i].text = ""
            par = cells[i].paragraphs[0]
            add_runs(par, text, size=8)
            for r in par.runs:
                r.font.size = Pt(8)


def render_markdown(doc: Document, text: str) -> int:
    lines = text.splitlines()
    i = 0
    tables = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # code fence
        if stripped.startswith("```"):
            i += 1
            buf = []
            while i < len(lines) and not lines[i].strip().startswith("```"):
                buf.append(lines[i])
                i += 1
            i += 1
            for code_line in buf:
                par = doc.add_paragraph()
                par.paragraph_format.space_after = Pt(0)
                r = par.add_run(code_line if code_line.strip() else " ")
                r.font.name = "Consolas"
                r.font.size = Pt(7.5)
            continue

        # table
        if stripped.startswith("|") and i + 1 < len(lines) and set(lines[i + 1].strip()) <= set("|-: "):
            rows = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                cells = [c.strip() for c in lines[i].strip().strip("|").split("|")]
                if not set("".join(cells)) <= set("-: "):
                    rows.append(cells)
                i += 1
            render_table(doc, rows)
            tables += 1
            continue

        if not stripped:
            i += 1
            continue

        if stripped.startswith("---") or stripped.startswith("***"):
            i += 1
            continue

        m = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if m:
            level = len(m.group(1))
            par = doc.add_paragraph()
            par.paragraph_format.space_before = Pt(10 if level > 1 else 14)
            par.paragraph_format.space_after = Pt(4)
            r = par.add_run(re.sub(r"[*`]", "", m.group(2)))
            r.bold = True
            r.font.size = Pt(16 - 1.4 * level)
            if level <= 2:
                r.font.color.rgb = RGBColor(0x1A, 0x36, 0x5D)
            i += 1
            continue

        if stripped.startswith(">"):
            par = doc.add_paragraph()
            par.paragraph_format.left_indent = Pt(18)
            add_runs(par, stripped.lstrip("> ").rstrip(), size=9.5)
            for r in par.runs:
                r.italic = True
                r.font.size = Pt(9.5)
            i += 1
            continue

        b = BULLET.match(line)
        if b:
            indent, rest = b.group(1), b.group(2)
            cb = CHECKBOX.match(rest)
            par = doc.add_paragraph()
            par.paragraph_format.left_indent = Pt(12 + 12 * (len(indent) // 2))
            par.paragraph_format.space_after = Pt(2)
            if cb:
                mark = "☑ " if cb.group(1) == "x" else "☐ "
                add_runs(par, mark + cb.group(2), size=9.5)
            else:
                add_runs(par, "・" + rest, size=9.5)
            i += 1
            continue

        n = NUMBERED.match(line)
        if n:
            par = doc.add_paragraph()
            par.paragraph_format.left_indent = Pt(12 + 12 * (len(n.group(1)) // 2))
            par.paragraph_format.space_after = Pt(2)
            add_runs(par, f"{n.group(2)}. {n.group(3)}", size=9.5)
            i += 1
            continue

        par = doc.add_paragraph()
        par.paragraph_format.space_after = Pt(4)
        add_runs(par, stripped, size=9.5)
        i += 1
    return tables


def main() -> None:
    OUT_DIR.mkdir(exist_ok=True)
    doc = Document()
    style = doc.styles["Normal"]
    style.font.name = "Noto Sans CJK JP"
    style.font.size = Pt(9.5)

    # cover
    t = doc.add_paragraph()
    r = t.add_run("ISP ルールブック（Geneformer in silico perturbation）")
    r.bold = True
    r.font.size = Pt(20)
    head = doc.add_paragraph()
    head.add_run(f"まとめ日: {date.today().isoformat()}　リポジトリのコミット: {git_commit()}")
    doc.add_paragraph()
    doc.add_paragraph("収録したページ（この順で掲載しています）:")
    for i, (rel, title) in enumerate(SOURCES, 1):
        par = doc.add_paragraph()
        par.paragraph_format.space_after = Pt(2)
        add_runs(par, f"{i}. {title}（`{rel}`）", size=9.5)
    doc.add_paragraph()
    note = doc.add_paragraph()
    add_runs(note, "この文書は `templates/template_build_rulebook_docx.py` が "
                   "`docs/isp/` の Markdown から自動生成します。"
                   "内容の正は各ページ（とくに記述ルールと評価基準）です。"
                   "文書とページが食い違う場合はページが正です。", size=9)

    combined = [f"# ISP ルールブック（Geneformer in silico perturbation）\n",
                f"まとめ日: {date.today().isoformat()}　リポジトリのコミット: {git_commit()}\n"]

    tables = 0
    for n, (rel, title) in enumerate(SOURCES):
        path = ROOT / rel
        if not path.exists():
            print(f"[skip] {rel} not found")
            continue
        body = path.read_text()
        combined.append(f"\n\n<!-- ===== {rel} ===== -->\n\n" + body)
        brk = doc.add_paragraph()
        brk.add_run().add_break(WD_BREAK.PAGE)
        banner = doc.add_paragraph()
        rb = banner.add_run(f"【{n + 1}】{title}")
        rb.bold = True
        rb.font.size = Pt(15)
        rb.font.color.rgb = RGBColor(0x8B, 0x1A, 0x1A)
        src = doc.add_paragraph()
        add_runs(src, f"出典: `{rel}`", size=8.5)
        tables += render_markdown(doc, body)
        print(f"[ok] {rel}: {len(body.splitlines())} lines")

    (OUT_DIR / "RULEBOOK-ISP-ja.md").write_text("".join(combined))
    docx_path = OUT_DIR / "RULEBOOK-ISP-ja.docx"
    doc.save(str(docx_path))
    print(f"wrote {docx_path} ({tables} tables)")


if __name__ == "__main__":
    main()
