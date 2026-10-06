#!/usr/bin/env python3
"""Endogenous content and contamination fraction (task 8).

What "contamination" means here
-------------------------------
An ancient sample carries at least three populations of DNA:

  endogenous   - molecules from the organism we care about, genuinely ancient
  environmental- soil and gut bacteria that colonised the remains after burial
  handling     - modern DNA from excavators, museum staff and the lab

Only the first is evidence. The other two are usually the large majority: a
plague victim's tooth is mostly soil bacteria, and a well-preserved ancient
human library is often under 1% endogenous.

Two independent estimators, because each fails differently
----------------------------------------------------------
1. Competitive mapping. Map the library against the target genome and a panel
   of plausible contaminants simultaneously, and let reads go to whichever
   reference they fit best. This catches environmental and handling DNA, but it
   can only see contaminants that are in the panel, and it cannot distinguish a
   modern molecule of the target species from an ancient one.

2. Damage partitioning. Split the reads that map to the target into those
   carrying a terminal C->T and those that do not, then ask what fraction of
   the library is damaged. A modern contaminant of the *same species* is
   invisible to competitive mapping but shows up here as a deficit of damage.

Estimator 2 is the one that answers the question estimator 1 cannot, and it is
also the one that depends on our damage model being right - which is why the
pipeline runs it against a modern negative control.

Interpreting the damage-based fraction
--------------------------------------
If genuinely ancient molecules carry terminal damage at rate `d_ref` and modern
molecules carry none, then an observed terminal rate `r` implies an ancient
fraction of roughly `r / d_ref`. The weakness is that `d_ref` is not universal:
it depends on the sample's age, burial temperature and the library protocol
(UDG treatment removes the signal entirely). We therefore report the raw
observed rate alongside the inferred fraction, and we state the `d_ref` used
rather than hiding it inside a point estimate.

Usage
-----
    python scripts/authenticate.py --bam results/bam/S.bam --ref data/ref/r.fa \
        --flagstat results/bam/S.flagstat --trim-json results/trim/S.json \
        --out results/auth/S [--d-ref 0.30]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pysam

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
BASES = "ACGT"


def read_terminal_damage(read: pysam.AlignedSegment, ref: pysam.FastaFile,
                         window: int) -> tuple[bool, bool]:
    """Does this read carry a C->T at its 5' end or a G->A at its 3' end?

    Returns (has_5p_ct, has_3p_ga) considering the first/last `window` bases of
    the molecule in its original orientation.
    """
    pairs = []
    start = read.reference_start
    segment = ref.fetch(read.reference_name, start, read.reference_end).upper()
    for qpos, rpos in read.get_aligned_pairs(matches_only=True):
        rb = read.query_sequence[qpos].upper()
        fb = segment[rpos - start]
        if rb in BASES and fb in BASES:
            pairs.append((rb, fb))
    if not pairs:
        return False, False

    if read.is_reverse:
        pairs = [(a.translate(COMPLEMENT), b.translate(COMPLEMENT))
                 for a, b in reversed(pairs)]

    has_ct = any(fb == "C" and rb == "T" for rb, fb in pairs[:window])
    has_ga = any(fb == "G" and rb == "A" for rb, fb in pairs[-window:])
    return has_ct, has_ga


def damage_partition(bam_path: str, ref_path: str, min_mapq: int,
                     window: int, d_ref: float) -> dict:
    bam = pysam.AlignmentFile(bam_path, "rb")
    ref = pysam.FastaFile(ref_path)

    n = n_ct = n_ga = n_either = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        ct, ga = read_terminal_damage(read, ref, window)
        n += 1
        n_ct += ct
        n_ga += ga
        n_either += (ct or ga)

    bam.close()
    ref.close()

    if n == 0:
        return {"reads_assessed": 0}

    frac_either = n_either / n
    return {
        "reads_assessed": n,
        "window_bases": window,
        "fraction_with_5p_CT": n_ct / n,
        "fraction_with_3p_GA": n_ga / n,
        "fraction_with_any_terminal_damage": frac_either,
        "d_ref_assumed": d_ref,
        "inferred_ancient_fraction": min(frac_either / d_ref, 1.0) if d_ref else None,
        "inferred_modern_contamination": (
            max(0.0, 1.0 - min(frac_either / d_ref, 1.0)) if d_ref else None
        ),
        "caveat": (
            "A UDG-treated library has no terminal uracil left to detect, so a low "
            "damaged fraction here means 'damage removed in the lab' and NOT "
            "'modern contamination'. Read this number together with the library "
            "protocol recorded in the sample sheet."
        ),
    }


def parse_flagstat(path: str | None) -> dict:
    """Endogenous content straight from samtools flagstat."""
    if not path or not Path(path).exists():
        return {}
    total = mapped = dups = 0
    for line in Path(path).read_text().splitlines():
        if " in total " in line:
            total = int(line.split()[0])
        elif " mapped (" in line and "primary" not in line:
            mapped = int(line.split()[0])
        elif line.split()[1:3] == ["+", "0"] and "duplicates" in line and "primary" not in line:
            dups = int(line.split()[0])
    return {
        "reads_total": total,
        "reads_mapped": mapped,
        "endogenous_fraction": (mapped / total) if total else 0.0,
        "duplicates": dups,
        "duplicate_fraction": (dups / mapped) if mapped else 0.0,
    }


def parse_trim(path: str | None, reads_fastq: str | None, layout: str) -> dict:
    """Trimming statistics, including the paired-end MERGE RATE.

    Merge rate is reported but is NOT used as evidence of antiquity. Mates
    merge whenever the molecule is shorter than twice the read length, so the
    rate depends on the sequencing configuration as much as on the sample: a
    modern 300 bp insert read with 2x270 bp reads still merges. In this project
    the modern control merges at 86%, close to the ancient samples' 92%, so the
    discriminating quantity is the merged fragment LENGTH, not the merge rate.
    """
    out: dict = {}
    if path and Path(path).exists():
        try:
            d = json.loads(Path(path).read_text())
            out["reads_input"] = d["summary"]["input"]["reads"]
            out["reads_output_all"] = d["summary"]["output"]["reads"]
            out["input_mean_length"] = d["summary"]["input"].get("mean_length")
        except (json.JSONDecodeError, KeyError, TypeError):
            pass
    if reads_fastq and Path(reads_fastq).exists():
        import gzip
        with gzip.open(reads_fastq, "rt") as fh:
            out["reads_after_trim"] = sum(1 for _ in fh) // 4
    if layout == "PAIRED" and out.get("reads_input") and out.get("reads_after_trim") is not None:
        pairs = out["reads_input"] / 2
        out["read_pairs"] = pairs
        out["merge_rate"] = out["reads_after_trim"] / pairs
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--bam", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--flagstat")
    p.add_argument("--trim-json")
    p.add_argument("--reads-fastq")
    p.add_argument("--layout", default="SINGLE")
    p.add_argument("--min-mapq", type=int, default=25)
    p.add_argument("--window", type=int, default=2,
                   help="bases from each terminus counted as 'terminal'")
    p.add_argument("--d-ref", type=float, default=0.30,
                   help="terminal damage rate expected of a fully ancient, non-UDG library")
    a = p.parse_args()

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    result = {
        "bam": a.bam,
        "endogenous": parse_flagstat(a.flagstat),
        "trimming": parse_trim(a.trim_json, a.reads_fastq, a.layout),
        "damage_partition": damage_partition(a.bam, a.ref, a.min_mapq, a.window, a.d_ref),
    }

    with open(f"{a.out}.json", "w") as fh:
        json.dump(result, fh, indent=2)

    e = result["endogenous"]
    dp = result["damage_partition"]
    if e:
        print(f"{a.bam}: endogenous={e['endogenous_fraction']:.4f} "
              f"({e['reads_mapped']}/{e['reads_total']})")
    if dp.get("reads_assessed"):
        print(f"  damaged reads={dp['fraction_with_any_terminal_damage']:.4f} "
              f"-> ancient fraction ~{dp['inferred_ancient_fraction']:.3f} "
              f"(d_ref={dp['d_ref_assumed']})")


if __name__ == "__main__":
    main()
