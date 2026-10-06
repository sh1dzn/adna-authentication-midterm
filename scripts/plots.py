#!/usr/bin/env python3
"""Figures for the report: damage profiles and fragment length distributions.

Two figures carry the argument.

damage_profiles.png
    The classic aDNA plot: misincorporation rate against distance from the
    fragment terminus, 5' C->T on the left and 3' G->A on the right, one line
    per sample. The dashed line is the A->G control. A reader should be able to
    see the ancient samples separate from the controls without reading a table,
    and should be able to see the control lines staying flat underneath.

fragment_lengths.png
    Length distributions on a log x-axis with the long-read platform
    thresholds marked. The figure's job is to make the task 5 negative result
    visible: the entire distribution sits two orders of magnitude to the left of
    where long-read sequencing starts being useful.

Usage
-----
    python scripts/plots.py --samples config/samples.tsv --out results/figures
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Colour by what we expect, so the controls are visually distinct from the
# samples that should show damage.
PROTOCOL_STYLE = {
    "undeclared": {"color": "#762a83", "ls": "-", "lw": 1.8},
    "partial_udg": {"color": "#ef8a62", "ls": "-", "lw": 2.0},
    "full_udg": {"color": "#2166ac", "ls": "--", "lw": 1.6},
    "modern": {"color": "#4d4d4d", "ls": ":", "lw": 1.6},
}

PLATFORM_LINES = [
    (1_000, "ONT useful\n1 kb"),
    (10_000, "ONT typical\n10 kb"),
    (15_000, "PacBio HiFi\n15 kb"),
]


def read_samples(path: str) -> dict:
    rows = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            rows[f[0]] = {"species": f[3], "protocol": f[6], "expect": f[7]}
    return rows


_UNDECLARED_COLOURS = ["#762a83", "#1b7837", "#d6604d", "#4393c3", "#8c510a", "#35978f"]
_seen: dict[str, int] = {}


def style_for(protocol: str, name: str = "") -> dict:
    """Controls keep one fixed style; undeclared samples each get their own colour."""
    st = dict(PROTOCOL_STYLE.get(protocol, {"color": "#999999", "ls": "-", "lw": 1.5}))
    if protocol == "undeclared":
        idx = _seen.setdefault(name, len(_seen))
        st["color"] = _UNDECLARED_COLOURS[idx % len(_UNDECLARED_COLOURS)]
    return st


def plot_damage(samples: dict, out_path: Path) -> None:
    """Two rows: full scale (the strong signal) and zoomed (everything else).

    One strongly damaged library flattens every other line against zero on a
    shared axis, so the controls and the weak samples would be unreadable. The
    bottom row rescales to 0-3% and is where the negative controls are judged.
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for row, ymax in enumerate((None, 0.03)):
        ax5, ax3 = axes[row]
        for name, meta in samples.items():
            path = Path(f"results/damage/{name}.json")
            if not path.exists():
                continue
            d = json.loads(path.read_text())
            st = style_for(meta["protocol"], name)
            label = f"{name} ({meta['protocol']})"
            p5, p3 = d["profile_5p"], d["profile_3p"]
            ax5.plot([r["position"] for r in p5], [r["C>T"] for r in p5],
                     label=label if row == 0 else None, **st)
            ax3.plot([r["position"] for r in p3], [r["G>A"] for r in p3], **st)
        for ax in (ax5, ax3):
            ax.axhline(0.02, color="black", lw=0.7, ls="--", alpha=0.5)
            ax.grid(alpha=0.25, lw=0.5)
            if ymax:
                ax.set_ylim(-0.001, ymax)
        ax5.set_ylabel("misincorporation rate" + (" (zoom)" if ymax else ""))
    axes[0][0].set_title("5' terminus: C$\\rightarrow$T")
    axes[0][1].set_title("3' terminus: G$\\rightarrow$A")
    axes[0][0].legend(fontsize=8, frameon=False)
    axes[1][0].annotate("authentication threshold (0.02)", xy=(10, 0.0215), fontsize=8, alpha=0.7)
    for ax in axes[1]:
        ax.set_xlabel("distance from fragment terminus (bp)")
    fig.suptitle("Post-mortem deamination accumulates at fragment ends", fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_lengths(samples: dict, out_path: Path) -> None:
    fig, (ax, axlog) = plt.subplots(1, 2, figsize=(12, 5))

    for name, meta in samples.items():
        path = Path(f"results/lengths/{name}.json")
        if not path.exists():
            continue
        d = json.loads(path.read_text())
        hist = {int(k): v for k, v in d.get("histogram", {}).items()}
        if not hist:
            continue
        total = sum(hist.values())
        xs = sorted(hist)
        ys = [hist[x] / total for x in xs]
        st = style_for(meta["protocol"], name)
        label = f"{name} (median {d.get('median')} bp)"
        ax.plot(xs, ys, label=label, **st)
        axlog.plot(xs, ys, label=label, **st)

    ax.set_xlabel("fragment length (bp)")
    ax.set_ylabel("fraction of reads")
    ax.set_title("Fragment length distribution")
    ax.grid(alpha=0.25, lw=0.5)
    ax.legend(fontsize=8, frameon=False)

    axlog.set_xscale("log")
    axlog.set_xlim(20, 50_000)
    axlog.set_xlabel("fragment length (bp, log scale)")
    axlog.set_title("Why long reads cannot be used")
    axlog.grid(alpha=0.25, lw=0.5)
    for x, text in PLATFORM_LINES:
        axlog.axvline(x, color="black", lw=0.8, ls="--", alpha=0.6)
        axlog.annotate(text, xy=(x, axlog.get_ylim()[1] * 0.75), fontsize=7,
                       rotation=90, va="top", ha="right", alpha=0.8)

    fig.suptitle("Ancient fragments sit two orders of magnitude below the long-read regime",
                 fontsize=13)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--samples", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    samples = read_samples(a.samples)

    plot_damage(samples, out / "damage_profiles.png")
    plot_lengths(samples, out / "fragment_lengths.png")
    print(f"wrote {out}/damage_profiles.png and {out}/fragment_lengths.png")


if __name__ == "__main__":
    main()
