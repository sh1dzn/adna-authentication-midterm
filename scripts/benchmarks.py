#!/usr/bin/env python3
"""Aggregate Snakemake benchmark files into one table: step, sample, wall time, peak RSS.

Usage: python scripts/benchmarks.py benchmarks > results/summary/runtime.tsv
"""
import csv
import sys
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else "benchmarks")
print("step\tsample\twall_seconds\tmax_rss_mb\tcpu_seconds")
for f in sorted(root.glob("*/*.tsv")):
    with open(f) as fh:
        row = next(csv.DictReader(fh, delimiter="\t"), None)
    if row:
        # Snakemake records NA for memory/CPU when a job finishes before its
        # sampling interval; report that honestly instead of inventing a value.
        def num(key, fmt):
            try:
                return format(float(row[key]), fmt)
            except ValueError:
                return "NA"
        print(f"{f.parent.name}\t{f.stem}\t{num('s', '.1f')}\t{num('max_rss', '.0f')}\t{num('cpu_time', '.1f')}")
