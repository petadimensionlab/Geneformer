# Model comparison and literature check — templates

Reusable templates for comparing two Geneformer models (or two conditions),
checking the target genes and cell populations against the published literature,
and writing it up as a readable report (DOCX/PDF, Japanese and English).

**These are TEMPLATES, not a finished analysis.** The gene verdicts, report prose
and literature queries embedded in them are a worked example from the `AD_spleen`
104M vs 316M run (2026-09). Re-derive every number and re-verify every PMID
before reusing them on another tissue, model or timepoint.

| Template | Purpose | Main output |
|---|---|---|
| `template_expression_by_celltype.py` | gene detection rates per cell type and condition (expression breadth) | `docs/quantization/expression/<TISSUE>.json` |
| `template_isp_model_compare.py` | separates precision / sampling / model gaps and counts sign flips | `isp-model-comparison/{gene_table.csv,summary.json}` |
| `template_literature_check.py` | Europe PMC retrieval with per-gene queries built in | `isp-model-comparison/literature/*.json` |
| `template_literature_adjudicate.py` | merges verdicts and machine-verifies the PMIDs | `isp-model-comparison/adjudication.{csv,json}` |
| `template_report_figures.py` | generates the figures (ja/en) | `isp-model-comparison/figures/<lang>/*.png` |
| `template_build_report_docx.py` | builds the report (DOCX, ja/en) | `isp-model-comparison/REPORT-*.docx` |

## Run order

```bash
.venv/bin/python templates/template_expression_by_celltype.py
.venv/bin/python templates/template_isp_model_compare.py
.venv/bin/python templates/template_literature_check.py --batch all
.venv/bin/python templates/template_literature_adjudicate.py
.venv/bin/python templates/template_report_figures.py --lang both
.venv/bin/python templates/template_build_report_docx.py --lang both
```

For PDFs (Japanese needs the Noto Sans CJK JP font):

```bash
cd docs/quantization/isp-model-comparison
soffice --headless --norestore --convert-to pdf \
    REPORT-ISP-104M-vs-316M.docx REPORT-ISP-104M-vs-316M-jp.docx
```

## Prerequisites and caveats

- **Outputs are not committed.** `docs/quantization/isp-model-comparison/` and
  `docs/quantization/expression/` are gitignored; only the templates are public.
- **Run the ISP itself first** (`analysis/07_*_early_isp.py` and friends); these
  templates read the ISP output CSVs.
- **`template_literature_check.py` uses the Europe PMC REST API** rather than page
  scraping, so it works even when the `web_extract` backend is unconfigured.
- **Always re-fetch every PMID.** Existence is not enough: check that the
  retrieved title matches the claimed subject (mismatches did occur).
- **`template_build_report_docx.py` numbers references in two passes**, so the
  printed list runs 1..N in order of first citation. Fixed key-based numbering
  leaves visible gaps (1, 2, 3, 5, 8, ...) that look broken to a reader.
- Environment variables: `ADPD_TISSUE` (default AD_spleen), `GENEFORMER_DIR`,
  `GENEFORMER_MODEL`, `GF_DTYPE`.
- Japanese figures and documents need a CJK font (Noto Sans CJK JP). If the font
  cannot be resolved the scripts fail loudly instead of shipping empty boxes.
