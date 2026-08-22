#!/usr/bin/env python
"""Resolve tissue workspace + prefix for the analysis scripts (04-07).

Input layout (this repo):
    input/<TISSUE>/h5ad/<something>.h5ad

The analysis scripts each need:
    ROOT    — workspace that holds h5ad/ (input), tokenized/, results/, runs/
    PREFIX  — short dataset name used in output filenames

Resolution precedence:
  1. ADPD_ROOT  — if set, use it verbatim (legacy; the workspace root).
  2. ADPD_TISSUE — tissue dir inside `input/` (e.g. PD_blood, AD_smallint).
     The `.h5ad` is searched as input/<TISSUE>/h5ad/*.h5ad.
  3. ADPD_PREFIX — used as the tissue name too (back-compat shim); if a matching
     `input/<PREFIX>` dir exists, behaves like ADPD_TISSUE.
  4. default — scan `input/` and pick the single PD_*/AD_* dir; error if multiple.

PREFIX (the output filename stem) is the tissue name unless ADPD_PREFIX is set.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

DEFAULT_INPUT = Path(__file__).resolve().parent.parent / "input"  # repo/input


def _tissue_dirs(base: Path) -> list[str]:
    if not base.is_dir():
        return []
    return sorted(
        d.name for d in base.iterdir() if d.is_dir() and d.name.startswith(("PD_", "AD_"))
    )


def resolve() -> tuple[Path, str]:
    root_env = os.environ.get("ADPD_ROOT")
    if root_env:
        root = Path(root_env)
        prefix = os.environ.get("ADPD_PREFIX", root.name)
        return root, prefix

    base = Path(os.environ.get("ADPD_INPUT", DEFAULT_INPUT))
    tissue: str | None = os.environ.get("ADPD_TISSUE") or os.environ.get("ADPD_PREFIX")

    if not tissue:
        dirs = _tissue_dirs(base)
        if len(dirs) == 1:
            tissue = dirs[0]
        elif len(dirs) > 1:
            sys.exit(
                "Multiple tissue dirs found; set ADPD_TISSUE (e.g. AD_blood) or "
                f"ADPD_ROOT (e.g. {base}/PD_blood): {dirs}"
            )
        elif dirs == []:
            sys.exit(
                f"no PD_*/AD_* dir under {base}; set ADPD_TISSUE with the h5ad folder "
                f"created under input/<TISSUE>/h5ad/*.h5ad"
            )

    assert tissue is not None, "tissue unexpectedly None"
    prefix: str = os.environ.get("ADPD_PREFIX", tissue)

    root = base / tissue
    if not (root / "h5ad").is_dir():
        sys.exit(f"no h5ad dir at {root / 'h5ad'}; expected input/<TISSUE>/h5ad/*.h5ad")
    return root, prefix


def find_h5ad(root: Path) -> Path:
    matches = list((root / "h5ad").glob("*.h5ad"))
    if not matches:
        sys.exit(f"no .h5ad under {root / 'h5ad'}")
    return matches[0]