#!/usr/bin/env python
"""Run a command under a load/energy profile and summarise it.

Samples GPU utilisation/power/clock/temperature (nvidia-smi), CPU utilisation
(/proc/stat), load average, system memory and the process-tree RSS at a fixed
interval, while streaming the child's output with an elapsed-time prefix into a
log file (so per-step cost can be derived afterwards).

Outputs, next to --out:
  <label>.json   metrics: wall time, energy, mean/max utilisation, peak RSS ...
  <label>.csv    raw samples
  <label>.log    child output, each line prefixed with "[+seconds]"

Usage:
    .venv/bin/python analysis/13_profile.py --label name --out results/profiles \
        --  <command> [args...]
Environment:
    PROF_INTERVAL (default 2.0 s), PROF_GPU_INDEX (default 0)
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path

INTERVAL = float(os.environ.get("PROF_INTERVAL", "2.0"))
GPU = os.environ.get("PROF_GPU_INDEX", "0")


def _read_cpu_total() -> tuple[int, int]:
    """Return (busy_jiffies, total_jiffies) from /proc/stat."""
    with open("/proc/stat") as fh:
        parts = fh.readline().split()[1:]
    vals = [int(x) for x in parts]
    idle = vals[3] + (vals[4] if len(vals) > 4 else 0)
    return sum(vals) - idle, sum(vals)


def _gpu_sample() -> dict:
    query = ("utilization.gpu", "power.draw", "clocks.current.sm", "temperature.gpu")
    try:
        out = subprocess.run(
            ["nvidia-smi", "-i", GPU, f"--query-gpu={','.join(query)}",
             "--format=csv,noheader,nounits"],
            capture_output=True, text=True, timeout=10, check=False,
        ).stdout.strip()
        if not out:
            return {}
        vals = [v.strip() for v in out.splitlines()[0].split(",")]
        return {
            "gpu_util_pct": float(vals[0]),
            "gpu_power_w": float(vals[1]) if vals[1] not in ("N/A", "[N/A]") else None,
            "gpu_sm_mhz": float(vals[2]) if vals[2] not in ("N/A", "[N/A]") else None,
            "gpu_temp_c": float(vals[3]) if vals[3] not in ("N/A", "[N/A]") else None,
        }
    except Exception:
        return {}


def _mem_available_gib() -> float:
    with open("/proc/meminfo") as fh:
        for line in fh:
            if line.startswith("MemAvailable:"):
                return int(line.split()[1]) / 2**20
    return float("nan")


def _tree_rss_gib(pgid: int) -> float:
    """RSS of every process in the child's process group (GiB)."""
    total_kb = 0
    for entry in Path("/proc").iterdir():
        if not entry.name.isdigit():
            continue
        try:
            with open(entry / "stat") as fh:
                fields = fh.read().rsplit(")", 1)[-1].split()
            if int(fields[2]) != pgid:  # field 5 overall = pgrp
                continue
            with open(entry / "status") as fh:
                for line in fh:
                    if line.startswith("VmRSS:"):
                        total_kb += int(line.split()[1])
                        break
        except (OSError, IndexError, ValueError):
            continue
    return total_kb / 2**20


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("command", nargs=argparse.REMAINDER)
    args = ap.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command:
        ap.error("no command given (use: -- <command>)")

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    label = args.label
    json_path = outdir / f"{label}.json"
    csv_path = outdir / f"{label}.csv"
    log_path = outdir / f"{label}.log"

    samples: list[dict] = []
    stop = threading.Event()
    log_fh = open(log_path, "w")

    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
        start_new_session=True,  # own process group so RSS accounting works
    )
    pgid = os.getpgid(proc.pid)
    t0 = time.time()

    def reader() -> None:
        assert proc.stdout is not None
        for line in proc.stdout:
            stamp = f"[+{time.time() - t0:7.1f}s] "
            log_fh.write(stamp + line)
            log_fh.flush()
            sys.stdout.write(stamp + line)
            sys.stdout.flush()

    def sampler() -> None:
        prev_busy, prev_total = _read_cpu_total()
        while not stop.wait(INTERVAL):
            busy, total = _read_cpu_total()
            dcpu = max(total - prev_total, 1)
            cpu_pct = 100.0 * (busy - prev_busy) / dcpu
            prev_busy, prev_total = busy, total
            row = {
                "t": round(time.time() - t0, 2),
                "cpu_pct": round(cpu_pct, 2),
                "loadavg1": os.getloadavg()[0],
                "mem_avail_gib": round(_mem_available_gib(), 2),
                "rss_gib": round(_tree_rss_gib(pgid), 3),
                **_gpu_sample(),
            }
            samples.append(row)

    th_read = threading.Thread(target=reader, daemon=True)
    th_samp = threading.Thread(target=sampler, daemon=True)
    th_read.start()
    th_samp.start()

    def forward(signum, _frame):  # keep Ctrl-C / TERM behaviour predictable
        proc.send_signal(signum)

    for sig in (signal.SIGINT, signal.SIGTERM):
        signal.signal(sig, forward)

    try:
        rc = proc.wait()
    finally:
        stop.set()
        th_samp.join(timeout=INTERVAL * 3)
        th_read.join(timeout=30)
        log_fh.close()
        wall = time.time() - t0

    # ---- summarise
    def col(key: str) -> list[float]:
        return [s[key] for s in samples if isinstance(s.get(key), (int, float))]

    powers = col("gpu_power_w")
    energy_wh = None
    if powers and len(samples) > 1:
        # trapezoidal integration over the sampled window
        energy_wh = sum(powers) / len(powers) * (wall / 3600.0)
    utils = col("gpu_util_pct")
    cpus = col("cpu_pct")
    rss = col("rss_gib")
    mem = col("mem_avail_gib")
    temps = col("gpu_temp_c")

    summary = {
        "label": label,
        "command": command,
        "wall_s": round(wall, 2),
        "wall_hms": time.strftime("%H:%M:%S", time.gmtime(wall)),
        "exit_code": rc,
        "samples": len(samples),
        "gpu_util_pct": {
            "mean": round(sum(utils) / len(utils), 2) if utils else None,
            "max": max(utils) if utils else None,
            "active_fraction": (round(sum(1 for u in utils if u > 10) / len(utils), 3)
                                if utils else None),
        },
        "gpu_power_w": {
            "mean": round(sum(powers) / len(powers), 2) if powers else None,
            "max": max(powers) if powers else None,
            "min": min(powers) if powers else None,
        },
        "gpu_energy_wh": round(energy_wh, 3) if energy_wh is not None else None,
        "gpu_temp_c_max": max(temps) if temps else None,
        "cpu_pct": {
            "mean": round(sum(cpus) / len(cpus), 2) if cpus else None,
            "max": max(cpus) if cpus else None,
        },
        "rss_peak_gib": max(rss) if rss else None,
        "mem_available_min_gib": min(mem) if mem else None,
        "interval_s": INTERVAL,
    }
    json_path.write_text(json.dumps(summary, indent=2))
    if samples:
        with open(csv_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=sorted({k for s in samples for k in s}))
            writer.writeheader()
            writer.writerows(samples)

    print("\n=== profile: " + label + " ===")
    print(json.dumps(summary, indent=2))
    print(f"\nlog  -> {log_path}\ncsv  -> {csv_path}\njson -> {json_path}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
