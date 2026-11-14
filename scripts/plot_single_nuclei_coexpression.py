"""Dot plot of EFHC1, CACNA1E and FOXJ1 across selected human brain cell types.

Reads the tables written by hpa_single_nuclei_coexpression.py.
Usage: python plot_single_nuclei_coexpression.py <coexpression_dir> <out_prefix>
"""
import csv
import sys
from pathlib import Path

import ajur_figure_style
import matplotlib.pyplot as plt

# Categorical slots 1-3 of the reference palette (validated all-pairs, light mode).
SERIES = [("EFHC1", "#2a78d6", "o"), ("CACNA1E", "#eb6834", "s"),
          ("FOXJ1", "#1baf7a", "^")]
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
FLOOR = 0.1

ROWS = [  # (label, source, key)
    ("Ependymal cells", "type", "ependymal cells"),
    ("Choroid plexus epithelial cells", "type", "choroid plexus epithelial cells"),
    ("Thalamic PVALB+ GABAergic (c-13)", "cluster", "c-13"),
    ("Thalamic projection neurons", "type", "thalamic projection neurons"),
    ("Dopaminergic neurons", "type", "dopaminergic"),
    ("Forebrain excitatory neurons", "type", "forebrain excitatory"),
    ("MGE-derived interneurons", "type", "mge derived interneuron"),
]


def main(d, out_prefix):
    d = Path(d)
    types = {r["cell_type"]: r for r in csv.DictReader(
        open(d / "cell_type_nCPM.tsv", encoding="utf-8"), delimiter="\t")}
    clusters = {r["cluster"]: r for r in csv.DictReader(
        open(d / "cluster_nCPM_selected.tsv", encoding="utf-8"), delimiter="\t")
        if r["group"] == "thalamic_reticular_like"}

    ajur_figure_style.apply()
    fig, ax = plt.subplots(figsize=(6.0, 3.2), dpi=300)
    n = len(ROWS)
    for i, (label, src, key) in enumerate(ROWS):
        rec = types[key] if src == "type" else clusters[key]
        y = n - 1 - i
        for j, (gene, color, marker) in enumerate(SERIES):
            v = max(float(rec[gene.split()[0]]), FLOOR)
            ax.scatter(v, y + (j - 1) * 0.22, s=34, color=color, marker=marker,
                       edgecolors="white", linewidths=0.8, zorder=3,
                       label=gene if i == 0 else None)
    ax.axhline(n - 2.5, color=MUTED, lw=0.6, ls=(0, (2, 2)))
    ax.text(1.02, n - 1.5, "ciliated", transform=ax.get_yaxis_transform(),
            color=MUTED, va="center")
    ax.text(1.02, (n - 3) / 2, "neuronal", transform=ax.get_yaxis_transform(),
            color=MUTED, va="center")
    ax.set_xscale("log")
    ax.set_xlim(0.08, 2000)
    ax.set_yticks(range(n))
    ax.set_yticklabels([r[0] for r in reversed(ROWS)], color=INK)
    ax.set_xlabel("Expression in adult human brain nuclei (nCPM, log scale)", color=MUTED)
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors=MUTED)
    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.45, 1.13), ncol=3, frameon=False,
                    handletextpad=0.3, columnspacing=1.4)
    for t in leg.get_texts():  # gene symbols in italics
        t.set_fontstyle("italic")
    # Values below 0.1 nCPM are drawn at 0.1; data source and this note go in the caption.
    fig.tight_layout()
    for ext in ("png", "pdf"):
        # No creation date in the PDF, so that a rerun gives an identical file
        fig.savefig(f"{out_prefix}.{ext}", bbox_inches="tight", facecolor="white",
                    metadata={"CreationDate": None} if ext == "pdf" else None)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
