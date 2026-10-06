#!/usr/bin/env python3
"""Measure post-mortem deamination damage at fragment ends (task 7).

The biology we are measuring
---------------------------
After death, cytosine in DNA loses its amino group (hydrolytic deamination) and
becomes uracil. A polymerase reads uracil as thymine, so the sequenced molecule
carries C->T where the original had C.

This happens everywhere in the molecule, but it *accumulates at fragment ends*
for a mechanical reason: the ends of an ancient fragment are frayed into short
single-stranded overhangs, and deamination is roughly two orders of magnitude
faster in single-stranded DNA than in the double-stranded interior. So the
signal we look for is not "more C->T than expected" but specifically "a C->T
rate that decays from the fragment terminus inwards".

For a conventional double-stranded library the signature is asymmetric:

    C->T rising at the 5' end   and   G->A rising at the 3' end

because the 3' G->A is the same deamination event observed on the complementary
strand after the second-strand synthesis fills in the overhang.

What this script reports
------------------------
For every position i counting in from each end, the misincorporation rate

    rate_C>T(i) = (# reads with ref C at i and read T) / (# reads with ref C at i)

together with all twelve substitution types, so that the C->T term can be
compared against a background of substitutions that damage does not cause.
A library with genuine ancient molecules shows a C->T terminal spike against a
flat background; a modern library shows a flat line at the sequencing error
rate; and a UDG-treated ancient library shows the flat line too, because the
enzyme excised the uracil before sequencing.

Orientation
-----------
pysam reports a reverse-strand alignment already reverse-complemented onto the
forward reference, which would put the molecule's 5' end on the right. We undo
that, so index 0 is always the 5' end of the original molecule.

Usage
-----
    python scripts/damage_profile.py in.bam ref.fa --out results/damage/SAMPLE \
        --length 25 --min-mapq 25
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import pysam

BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
# The two substitutions deamination produces, and a control that it does not.
DAMAGE_KEYS = ("C>T", "G>A")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def aligned_pairs(read: pysam.AlignedSegment, ref: pysam.FastaFile,
                  min_baseq: int = 0) -> list[tuple[str, str]]:
    """Return [(read_base, ref_base), ...] in the ORIGINAL molecule's orientation.

    Insertions and deletions are skipped: a position only contributes if the
    read has a base aligned against a reference base. Soft-clipped bases are
    excluded by pysam from the aligned pairs, which is what we want, since a
    clipped base has no reference to be compared against.

    Base-quality filtering, and why it does not erase the signal
    -----------------------------------------------------------
    The trimming step deliberately leaves the 5' terminus un-quality-trimmed
    (--preserve5p), because that is where deamination lives and trimming it
    would delete the evidence. The cost is that genuinely low-quality terminal
    bases survive too, and those inflate *every* substitution type at position
    1 - including A->G, which deamination cannot produce.

    Dropping low-quality bases separates the two. A deamination event is not a
    sequencing error: the uracil is read as a real, confidently-called T, so it
    carries a high base quality. Sequencer noise at the read terminus does not.
    Filtering on base quality therefore removes the artefact while leaving the
    damage, which is exactly what the A->G control is there to verify.

    Both numerator and denominator are filtered, since a base excluded from the
    mismatch count must also leave the count of opportunities.
    """
    quals = read.query_qualities
    # One reference fetch per read, not per base. Profiling showed a per-base
    # ref.fetch() dominated runtime (cost scaled with total bases, not reads).
    start = read.reference_start
    segment = ref.fetch(read.reference_name, start, read.reference_end).upper()
    pairs: list[tuple[str, str]] = []
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        if min_baseq and quals is not None and quals[qpos] < min_baseq:
            continue
        read_base = read.query_sequence[qpos].upper()
        ref_base = segment[rpos - start]
        if read_base not in BASES or ref_base not in BASES:
            continue
        pairs.append((read_base, ref_base))

    if read.is_reverse:
        # Put the molecule back in its sequenced orientation so that index 0 is
        # the 5' terminus rather than whichever end happened to map leftmost.
        pairs = [
            (rb.translate(COMPLEMENT), fb.translate(COMPLEMENT))
            for rb, fb in reversed(pairs)
        ]
    return pairs


def profile(bam_path: str, ref_path: str, length: int, min_mapq: int,
            min_baseq: int = 0) -> dict:
    """Accumulate per-position substitution counts from both fragment ends."""
    ref = pysam.FastaFile(ref_path)
    bam = pysam.AlignmentFile(bam_path, "rb")

    # numerator[end][i][("C","T")] and denominator[end][i]["C"]
    subs = {"5p": [Counter() for _ in range(length)], "3p": [Counter() for _ in range(length)]}
    refcount = {"5p": [Counter() for _ in range(length)], "3p": [Counter() for _ in range(length)]}

    n_reads = 0
    n_skipped = 0
    lengths: Counter = Counter()

    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            n_skipped += 1
            continue
        if read.mapping_quality < min_mapq:
            n_skipped += 1
            continue
        pairs = aligned_pairs(read, ref, min_baseq)
        if not pairs:
            n_skipped += 1
            continue

        n_reads += 1
        lengths[read.query_length] += 1

        # 5' end: walk inwards from the start of the molecule.
        for i, (read_base, ref_base) in enumerate(pairs[:length]):
            refcount["5p"][i][ref_base] += 1
            if read_base != ref_base:
                subs["5p"][i][f"{ref_base}>{read_base}"] += 1

        # 3' end: walk inwards from the far terminus.
        for i, (read_base, ref_base) in enumerate(reversed(pairs[-length:])):
            refcount["3p"][i][ref_base] += 1
            if read_base != ref_base:
                subs["3p"][i][f"{ref_base}>{read_base}"] += 1

    bam.close()
    ref.close()

    def rates(end: str) -> list[dict]:
        out = []
        for i in range(length):
            row: dict = {"position": i + 1}
            for frm in BASES:
                denom = refcount[end][i][frm]
                for to in BASES:
                    if frm == to:
                        continue
                    key = f"{frm}>{to}"
                    row[key] = (subs[end][i][key] / denom) if denom else 0.0
            row["_depth"] = sum(refcount[end][i][b] for b in BASES)
            out.append(row)
        return out

    return {
        "bam": bam_path,
        "reads_used": n_reads,
        "reads_skipped": n_skipped,
        "min_mapq": min_mapq,
        "min_baseq": min_baseq,
        "positions": length,
        "profile_5p": rates("5p"),
        "profile_3p": rates("3p"),
        "length_histogram": dict(sorted(lengths.items())),
    }


def terminal_summary(result: dict) -> dict:
    """Condense the profile into the few numbers the report actually argues from.

    `delta` is the quantity that separates ancient from modern: the terminal
    C->T rate minus the interior rate. Damage decays within roughly ten bases,
    so positions 10+ serve as the sample's own internal background, which is
    more honest than comparing against a theoretical error rate.
    """
    def end_stats(prof: list[dict], key: str) -> dict:
        terminal = prof[0][key] if prof else 0.0
        interior = [row[key] for row in prof[10:]] or [0.0]
        background = sum(interior) / len(interior)
        return {
            "terminal": terminal,
            "interior_background": background,
            "delta": terminal - background,
            "first_five": [prof[i][key] for i in range(min(5, len(prof)))],
        }

    return {
        "ct_5p": end_stats(result["profile_5p"], "C>T"),
        "ga_3p": end_stats(result["profile_3p"], "G>A"),
        # A control substitution: deamination cannot produce it, so a terminal
        # rise here would mean our estimator is picking up an artefact instead.
        "ag_5p_control": end_stats(result["profile_5p"], "A>G"),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bam")
    p.add_argument("ref")
    p.add_argument("--out", required=True, help="output prefix (directory is created)")
    p.add_argument("--length", type=int, default=25, help="positions to profile from each end")
    p.add_argument("--min-mapq", type=int, default=25)
    p.add_argument("--min-baseq", type=int, default=20,
                   help="drop aligned bases below this Phred score; see aligned_pairs()")
    a = p.parse_args()

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    result = profile(a.bam, a.ref, a.length, a.min_mapq, a.min_baseq)
    result["summary"] = terminal_summary(result)

    with open(f"{a.out}.json", "w") as fh:
        json.dump(result, fh, indent=2)

    # A flat table, so the profile can be plotted or eyeballed without parsing JSON.
    with open(f"{a.out}.misincorporation.tsv", "w") as fh:
        keys = [f"{f}>{t}" for f in BASES for t in BASES if f != t]
        fh.write("end\tposition\tdepth\t" + "\t".join(keys) + "\n")
        for end, prof in (("5p", result["profile_5p"]), ("3p", result["profile_3p"])):
            for row in prof:
                fh.write(
                    f"{end}\t{row['position']}\t{row['_depth']}\t"
                    + "\t".join(f"{row[k]:.6f}" for k in keys)
                    + "\n"
                )

    s = result["summary"]
    print(f"{a.bam}: {result['reads_used']} reads used, {result['reads_skipped']} skipped")
    print(f"  5' C>T terminal={s['ct_5p']['terminal']:.4f} "
          f"background={s['ct_5p']['interior_background']:.4f} "
          f"delta={s['ct_5p']['delta']:+.4f}")
    print(f"  3' G>A terminal={s['ga_3p']['terminal']:.4f} "
          f"background={s['ga_3p']['interior_background']:.4f} "
          f"delta={s['ga_3p']['delta']:+.4f}")
    print(f"  5' A>G control delta={s['ag_5p_control']['delta']:+.4f} (should be ~0)")


if __name__ == "__main__":
    main()
