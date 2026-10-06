#!/usr/bin/env python3
"""Fragment length distributions, and the quantitative case against long reads (task 5).

Why this script exists
----------------------
Task 5 asks us to show that long-read sequencing is essentially inapplicable to
ancient DNA, and to support that argument with numbers rather than assertion.
The argument has two legs, and this script supplies the first and frames the
second.

Leg 1 - there is no long molecule to read.
    Post-mortem depurination cuts the phosphodiester backbone continuously, so
    fragment length decays roughly exponentially with time since death. Ancient
    libraries are typically 40-70 bp. A platform whose advantage only appears
    above several kilobases has nothing to work with: you cannot read a 10 kb
    molecule out of a sample whose longest molecule is 150 bp. We measure the
    distribution directly and report what fraction of fragments would clear the
    length thresholds at which ONT and PacBio become worthwhile.

Leg 2 - even if a long molecule existed, the error rate would bury the signal.
    Authentication depends on measuring a C->T rate of a few per cent at
    fragment ends. Raw ONT single-pass error is ~5% and PacBio CLR ~10-15%,
    both an order of magnitude above the signal we are trying to detect, and
    both producing substitution errors that are not strand-symmetric in a way
    that cleanly separates from deamination. PacBio HiFi reaches ~0.1% but only
    by sequencing one molecule many times, which needs a template long enough
    to circularise - exactly what leg 1 rules out. We report the signal-to-error
    ratio so the report can state this as a measured quantity.

Note on measuring aDNA fragment length
--------------------------------------
For a collapsed single-end library the sequenced read length IS the molecule
length, because the molecule was shorter than the read cycle count and the
sequencer ran off the end of it into adapter. That is why the adapter trimming
step must collapse overlapping pairs rather than merely trim: after collapsing,
read length is a direct physical measurement of the ancient fragment.

Reads that reach the full cycle length are right-censored - the molecule was at
least that long but we cannot say how much longer. We report the censored
fraction rather than silently treating it as a true length.

Usage
-----
    python scripts/fragment_lengths.py in.bam --out results/lengths/SAMPLE --min-mapq 25
"""
from __future__ import annotations

import argparse
import json
import statistics
from collections import Counter
from pathlib import Path

import pysam

# Thresholds below which each platform loses its reason to exist. These are the
# lengths at which the technology's advantage over Illumina actually appears,
# not the vendor's absolute minimum.
PLATFORM_THRESHOLDS = {
    "ONT_useful_1kb": 1_000,
    "ONT_typical_10kb": 10_000,
    "PacBio_HiFi_15kb": 15_000,
    "PacBio_CLR_20kb": 20_000,
}

# Published single-pass substitution error rates, used for the signal-to-error
# comparison. Cited in the report rather than measured here.
PLATFORM_ERROR_RATES = {
    "Illumina_NextSeq": 0.001,
    "ONT_R10_simplex": 0.05,
    "PacBio_CLR": 0.12,
    "PacBio_HiFi": 0.001,
}


def collect(bam_path: str, min_mapq: int) -> Counter:
    bam = pysam.AlignmentFile(bam_path, "rb")
    lengths: Counter = Counter()
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        lengths[read.query_length] += 1
    bam.close()
    return lengths


def describe(lengths: Counter, ct_delta: float | None = None) -> dict:
    """Summarise the distribution and answer the long-read question with it."""
    if not lengths:
        return {"n": 0}

    expanded = []
    for length, count in lengths.items():
        expanded.extend([length] * count)
    expanded.sort()
    n = len(expanded)

    def pct(p: float) -> int:
        return expanded[min(int(n * p), n - 1)]

    longest = expanded[-1]
    # Right-censoring: reads sitting at the maximum observed length may be
    # truncated by the read cycle count rather than by the molecule ending.
    censored = lengths[longest] / n

    fractions_above = {
        name: sum(c for l, c in lengths.items() if l >= thr) / n
        for name, thr in PLATFORM_THRESHOLDS.items()
    }

    out = {
        "n": n,
        "mean": statistics.mean(expanded),
        "median": statistics.median(expanded),
        "stdev": statistics.stdev(expanded) if n > 1 else 0.0,
        "min": expanded[0],
        "max": longest,
        "p01": pct(0.01), "p25": pct(0.25), "p50": pct(0.50),
        "p75": pct(0.75), "p95": pct(0.95), "p99": pct(0.99),
        "fraction_at_max_length_right_censored": censored,
        "fraction_above_platform_threshold": fractions_above,
        "histogram": dict(sorted(lengths.items())),
    }

    # Leg 2, stated as a ratio the report can quote directly.
    if ct_delta is not None:
        out["signal_to_error"] = {
            "damage_signal_ct_delta": ct_delta,
            "ratio_vs_platform_error": {
                name: (ct_delta / err if err else None)
                for name, err in PLATFORM_ERROR_RATES.items()
            },
            "note": (
                "A ratio below ~1 means the platform's own substitution error is "
                "larger than the damage signal being measured, so the terminal "
                "C>T spike cannot be separated from miscalled bases."
            ),
        }
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bam")
    p.add_argument("--out", required=True)
    p.add_argument("--min-mapq", type=int, default=25)
    p.add_argument("--ct-delta", type=float, default=None,
                   help="terminal C>T excess from damage_profile.py, for the signal-to-error ratio")
    a = p.parse_args()

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    lengths = collect(a.bam, a.min_mapq)
    result = describe(lengths, a.ct_delta)
    result["bam"] = a.bam

    with open(f"{a.out}.json", "w") as fh:
        json.dump(result, fh, indent=2)
    with open(f"{a.out}.histogram.tsv", "w") as fh:
        fh.write("length\tcount\n")
        for length, count in sorted(lengths.items()):
            fh.write(f"{length}\t{count}\n")

    if result["n"] == 0:
        print(f"{a.bam}: no reads passed filters")
        return

    print(f"{a.bam}: n={result['n']} mean={result['mean']:.1f} "
          f"median={result['median']} max={result['max']}")
    for name, frac in result["fraction_above_platform_threshold"].items():
        print(f"  fraction >= {name}: {frac:.2e}")


if __name__ == "__main__":
    main()
