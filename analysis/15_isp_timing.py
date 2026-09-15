#!/usr/bin/env python
"""Per-unit cost breakdown of an ISP run from its profiled log.

`analysis/13_profile.py` prefixes every output line with `[+  123.4s]`. The ISP
script prints one line per (timepoint, gene) unit as `[3m/APOE] max_ncells=...`,
so consecutive markers give the wall time of each gene × timepoint unit. That is
the number to compare between dtypes/models, independent of startup cost.

Usage:
    .venv/bin/python analysis/15_isp_timing.py docs/quantization/profiles/isp-316m-bf16-full.log \
        [--label 316m-bf16] [--out docs/quantization/profiles/isp-timings.json] [more logs...]
"""
from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path

STAMP = re.compile(r"^\[\+\s*([0-9.]+)s\]")
UNIT = re.compile(r"\[(3m|4p5m|6m|9m|12m)/([A-Za-z0-9_.-]+)\]")


def parse(path: Path) -> dict:
    rows: list[tuple[float, str, str]] = []
    for line in path.read_text(errors="replace").splitlines():
        m_stamp = STAMP.match(line)
        m_unit = UNIT.search(line)
        if m_stamp and m_unit:
            rows.append((float(m_stamp.group(1)), m_unit.group(1), m_unit.group(2)))
    # deduplicate: the marker line is printed once per unit
    seen: set[tuple[str, str]] = set()
    uniq: list[tuple[float, str, str]] = []
    for t, tp, gene in rows:
        if (tp, gene) in seen:
            continue
        seen.add((tp, gene))
        uniq.append((t, tp, gene))

    if len(uniq) < 2:
        return {"units": len(uniq), "error": "not enough unit markers to time"}

    deltas = [uniq[i + 1][0] - uniq[i][0] for i in range(len(uniq) - 1)]
    first_start = uniq[0][0]
    out = {
        "log": str(path),
        "units": len(uniq),
        "startup_s": round(first_start, 1),
        "first_unit_s": round(uniq[1][0] - uniq[0][0], 2),
        "per_unit_s": {
            "mean": round(statistics.mean(deltas), 2),
            "median": round(statistics.median(deltas), 2),
            "min": round(min(deltas), 2),
            "max": round(max(deltas), 2),
            "stdev": round(statistics.pstdev(deltas), 2),
        },
        "span_first_to_last_s": round(uniq[-1][0] - uniq[0][0], 1),
    }
    by_tp: dict[str, list[float]] = {}
    for i, (_, tp, _gene) in enumerate(uniq[:-1]):
        by_tp.setdefault(tp, []).append(deltas[i])
    out["per_timepoint_unit_s"] = {
        tp: {"n": len(v), "mean": round(statistics.mean(v), 2)} for tp, v in by_tp.items()
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("logs", nargs="+")
    ap.add_argument("--label", action="append", default=None,
                    help="repeat once per log; defaults to the log stem")
    ap.add_argument("--out", default="docs/quantization/profiles/isp-timings.json")
    args = ap.parse_args()

    labels = args.label or [Path(p).stem for p in args.logs]
    result = {}
    for label, p in zip(labels, args.logs):
        r = parse(Path(p))
        result[label] = r
        if "error" in r:
            print(f"{label:24s} {r['error']} (units={r['units']})")
            continue
        pu = r["per_unit_s"]
        print(f"{label:24s} units={r['units']:4d} startup={r['startup_s']:7.1f}s "
              f"per-unit mean={pu['mean']:7.2f}s median={pu['median']:7.2f}s "
              f"min={pu['min']:6.2f} max={pu['max']:6.2f} span={r['span_first_to_last_s']:7.1f}s")
        print(f"{'':24s} per timepoint: " +
              ", ".join(f"{tp}={v['mean']:.2f}s(n={v['n']})"
                        for tp, v in sorted(r["per_timepoint_unit_s"].items())))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    prev = json.loads(out.read_text()) if out.exists() else {}
    prev.update(result)
    out.write_text(json.dumps(prev, indent=2))
    print(f"\nresult -> {out}")


if __name__ == "__main__":
    main()
