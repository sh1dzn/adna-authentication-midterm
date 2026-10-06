#!/usr/bin/env python3
"""Reference bias under damage-aware vs standard alignment (task 6).

The problem
-----------
An aligner decides whether a read belongs to a locus by counting how badly it
disagrees with the reference. Ancient reads disagree more than modern ones for
two reasons that have nothing to do with where they came from: they are short,
so there is less sequence to anchor them, and they carry deamination
substitutions, so they mismatch even where they align correctly.

The consequence is systematic rather than random. A read carrying the
non-reference allele at a variant site already has one mismatch; add a
deamination event and it crosses the aligner's threshold and is discarded,
while the read carrying the reference allele survives. The surviving pileup is
therefore enriched for the reference allele, and the genotypes called from it
are pulled toward the reference. This is reference bias, and in ancient DNA it
is strong enough to distort population-genetic conclusions.

Why the aDNA settings look reckless and are not
-----------------------------------------------
    bwa aln -l 1024 -n 0.01 -o 2

`-l 1024` pushes the seed length past any realistic read length, which disables
seeding. Seeding requires the first 32 bases to match near-perfectly, and in an
ancient read those are exactly the bases most likely to carry deamination, so
seeding preferentially rejects the genuinely ancient molecules. `-n 0.01`
loosens the permitted edit distance so that a short read carrying damage can
still be placed. Both changes trade specificity for the ability to see ancient
molecules at all; the mapping-quality filter applied afterwards is what buys
the specificity back.

What this script measures
-------------------------
It compares two alignments of the *same* trimmed reads - one permissive, one
with bwa's defaults - and reports:

  * how many reads each setting places, and the mapping-rate difference;
  * the per-read mismatch distribution under each;
  * the terminal damage rate among reads placed ONLY by the permissive run.

The third number is the argument. If the reads that strict alignment throws
away are *more* damaged than the reads it keeps, then strict alignment is
discarding ancient molecules in preference to modern ones, and the alignment
step is itself selecting against the evidence the project is trying to measure.

Usage
-----
    python scripts/reference_bias.py --permissive perm.bam --strict strict.bam \
        --ref ref.fa --out results/refbias/SAMPLE
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pysam

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"


def read_stats(bam_path: str, ref_path: str, min_mapq: int) -> tuple[dict, dict]:
    """Return (per-read-name summary, aggregate stats) for one alignment."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref = pysam.FastaFile(ref_path)

    per_read: dict[str, dict] = {}
    mismatch_hist: Counter = Counter()
    mapq_hist: Counter = Counter()
    n_total = n_mapped = n_pass = 0

    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        n_total += 1
        if read.is_unmapped:
            continue
        n_mapped += 1
        mapq_hist[read.mapping_quality] += 1
        if read.mapping_quality < min_mapq:
            continue
        n_pass += 1

        pairs = []
        start = read.reference_start
        segment = ref.fetch(read.reference_name, start, read.reference_end).upper()
        for qpos, rpos in read.get_aligned_pairs(matches_only=True):
            rb = read.query_sequence[qpos].upper()
            fb = segment[rpos - start]
            if rb in BASES and fb in BASES:
                pairs.append((rb, fb))
        if read.is_reverse:
            pairs = [(a.translate(COMPLEMENT), b.translate(COMPLEMENT))
                     for a, b in reversed(pairs)]

        mismatches = sum(1 for rb, fb in pairs if rb != fb)
        mismatch_hist[mismatches] += 1
        has_ct = any(fb == "C" and rb == "T" for rb, fb in pairs[:2])
        has_ga = any(fb == "G" and rb == "A" for rb, fb in pairs[-2:])
        # Damage-SPECIFIC bookkeeping. Counting "reads with a terminal C>T" is
        # confounded: reads rescued by permissive alignment are rescued because
        # they carry extra mismatches of ANY kind, so that count rises even in a
        # modern library (it did, 123x, in our modern control). What separates
        # damage from divergence is the SPECTRUM of terminal mismatches: damage
        # makes them overwhelmingly C>T (5') and G>A (3'), divergence does not.
        term = [(rb, fb) for rb, fb in pairs[:2] if rb != fb] + \
               [(rb, fb) for rb, fb in pairs[-2:] if rb != fb]
        n_term = len(term)
        n_dmg = (sum(1 for rb, fb in pairs[:2] if fb == "C" and rb == "T")
                 + sum(1 for rb, fb in pairs[-2:] if fb == "G" and rb == "A"))

        per_read[read.query_name] = {
            "mapq": read.mapping_quality,
            "mismatches": mismatches,
            "length": read.query_length,
            "damaged": bool(has_ct or has_ga),
            "term_mm": n_term,
            "term_dmg": n_dmg,
        }

    bam.close()
    ref.close()

    total_mm = sum(k * v for k, v in mismatch_hist.items())
    n_mm_reads = sum(mismatch_hist.values())
    agg = {
        "reads_seen": n_total,
        "reads_mapped": n_mapped,
        "reads_passing_mapq": n_pass,
        "mapping_rate": (n_mapped / n_total) if n_total else 0.0,
        "pass_rate": (n_pass / n_total) if n_total else 0.0,
        "mean_mismatches": (total_mm / n_mm_reads) if n_mm_reads else 0.0,
        "mismatch_histogram": dict(sorted(mismatch_hist.items())),
        "damaged_fraction": (
            sum(1 for r in per_read.values() if r["damaged"]) / len(per_read)
            if per_read else 0.0
        ),
    }
    return per_read, agg


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--permissive", required=True)
    p.add_argument("--strict", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--min-mapq", type=int, default=25)
    a = p.parse_args()

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    perm_reads, perm_agg = read_stats(a.permissive, a.ref, a.min_mapq)
    strict_reads, strict_agg = read_stats(a.strict, a.ref, a.min_mapq)

    only_perm = set(perm_reads) - set(strict_reads)
    shared = set(perm_reads) & set(strict_reads)

    def damaged_frac(names: set[str], table: dict) -> float:
        if not names:
            return 0.0
        return sum(1 for n in names if table[n]["damaged"]) / len(names)

    def mean_len(names: set[str], table: dict) -> float:
        if not names:
            return 0.0
        return sum(table[n]["length"] for n in names) / len(names)

    def damage_share(names: set[str], table: dict) -> tuple[float | None, int]:
        """Share of terminal mismatches that are damage-type, and how many there were."""
        tot = sum(table[n]["term_mm"] for n in names)
        dmg = sum(table[n]["term_dmg"] for n in names)
        return ((dmg / tot) if tot else None), tot

    share_only, n_only = damage_share(only_perm, perm_reads)
    share_shared, n_shared = damage_share(shared, perm_reads)

    d_only = damaged_frac(only_perm, perm_reads)
    d_shared = damaged_frac(shared, perm_reads)

    result = {
        "permissive_bam": a.permissive,
        "strict_bam": a.strict,
        "min_mapq": a.min_mapq,
        "permissive": perm_agg,
        "strict": strict_agg,
        "reads_recovered_by_permissive_only": len(only_perm),
        "reads_shared": len(shared),
        "extra_yield_fraction": (
            (perm_agg["reads_passing_mapq"] - strict_agg["reads_passing_mapq"])
            / strict_agg["reads_passing_mapq"]
            if strict_agg["reads_passing_mapq"] else None
        ),
        "damaged_fraction_permissive_only": d_only,
        "damaged_fraction_shared": d_shared,
        "damage_enrichment_in_recovered_reads": (d_only / d_shared) if d_shared else None,
        "terminal_mismatches_recovered": n_only,
        "terminal_mismatches_shared": n_shared,
        "damage_type_share_of_terminal_mismatches_recovered": share_only,
        "damage_type_share_of_terminal_mismatches_shared": share_shared,
        "mean_length_permissive_only": mean_len(only_perm, perm_reads),
        "mean_length_shared": mean_len(shared, perm_reads),
        "interpretation": (
            "damage_enrichment_in_recovered_reads is CONFOUNDED and should not be read "
            "as evidence of damage: it is large even for a modern library, because reads "
            "rescued by permissive alignment carry extra mismatches of every kind. The "
            "damage-specific quantity is damage_type_share_of_terminal_mismatches_*: a "
            "high share among recovered reads means strict alignment is discarding "
            "deaminated molecules; a share near the modern control's means it is "
            "discarding divergent ones."
        ),
    }

    with open(f"{a.out}.json", "w") as fh:
        json.dump(result, fh, indent=2)

    print(f"permissive passing={perm_agg['reads_passing_mapq']}  "
          f"strict passing={strict_agg['reads_passing_mapq']}")
    print(f"  recovered only by permissive: {len(only_perm)}")
    print(f"  damaged fraction: recovered={d_only:.4f} shared={d_shared:.4f} "
          f"enrichment={result['damage_enrichment_in_recovered_reads'] or float('nan'):.2f}x")
    print(f"  damage-type share of terminal mismatches: recovered="
          f"{share_only if share_only is not None else float('nan'):.3f} (n={n_only}) "
          f"shared={share_shared if share_shared is not None else float('nan'):.3f} (n={n_shared})")
    print(f"  mean length: recovered={result['mean_length_permissive_only']:.1f} "
          f"shared={result['mean_length_shared']:.1f}")


if __name__ == "__main__":
    main()
