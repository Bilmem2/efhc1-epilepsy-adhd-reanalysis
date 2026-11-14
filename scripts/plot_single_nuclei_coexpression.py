"""Dot plot of EFHC1, CACNA1E and FOXJ1 in human brain cell types (Figure 3).

Panel a: adult brain, from the tables written by hpa_single_nuclei_coexpression.py.
Panel b (optional): first-trimester neurons, from fetal_neurons_pseudobulk.tsv
written by fetal_brain_single_cell.py.

Usage: python plot_single_nuclei_coexpression.py <coexpression_dir> <out_prefix> [<developmental_dir>]
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

ROWS_ADULT = [  # (label, source, key)
    ("Ependymal cells", "type", "ependymal cells"),
    ("Choroid plexus epithelial cells", "type", "choroid plexus epithelial cells"),
    (r"Thalamic $\mathit{PVALB}$+ GABAergic (c-13)", "cluster", "c-13"),
    ("Thalamic projection neurons", "type", "thalamic projection neurons"),
    ("Dopaminergic neurons", "type", "dopaminergic"),
    ("Forebrain excitatory neurons", "type", "forebrain excitatory"),
    ("MGE-derived interneurons", "type", "mge derived interneuron"),
]
ROWS_FETAL = [  # (label, group in fetal_neurons_pseudobulk.tsv)
    ("Thalamic neurons", "thalamus, all neurons"),
    (r"Thalamic $\mathit{GAD2}$+ neurons", "thalamus, GAD2+ neurons"),
    (r"Midbrain $\mathit{TH}$+ neurons", "midbrain, TH+ neurons"),
    ("Other midbrain neurons", "midbrain, other neurons"),
    ("Cortical neurons", "cortex, all neurons"),
]


def dot_panel(ax, records, labels):
    n = len(records)
    for i, rec in enumerate(records):
        y = n - 1 - i
        for j, (gene, color, marker) in enumerate(SERIES):
            v = max(float(rec[gene]), FLOOR)
            ax.scatter(v, y + (j - 1) * 0.22, s=34, color=color, marker=marker,
                       edgecolors="white", linewidths=0.8, zorder=3,
                       label=gene if i == 0 else None)
    ax.set_xscale("log")
    ax.set_xlim(0.08, 2000)
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(labels)), color=INK)
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    for s in ["top", "right", "left"]:
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors=MUTED)


def main(d, out_prefix, dev_dir=None):
    d = Path(d)
    types = {r["cell_type"]: r for r in csv.DictReader(
        open(d / "cell_type_nCPM.tsv", encoding="utf-8"), delimiter="\t")}
    clusters = {r["cluster"]: r for r in csv.DictReader(
        open(d / "cluster_nCPM_selected.tsv", encoding="utf-8"), delimiter="\t")
        if r["group"] == "thalamic_reticular_like"}
    adult = [types[k] if src == "type" else clusters[k] for _, src, k in ROWS_ADULT]

    ajur_figure_style.apply()
    if dev_dir:
        fetal = {r["group"]: r for r in csv.DictReader(
            open(Path(dev_dir) / "fetal_neurons_pseudobulk.tsv", encoding="utf-8"), delimiter="\t")
            if r["age"] == "all ages"}
        fetal_recs = [{g: fetal[k][f"{g}_per_million"] for g, _, _ in SERIES} for _, k in ROWS_FETAL]
        fig, (ax, bx) = plt.subplots(2, 1, figsize=(6.0, 4.9), dpi=300, height_ratios=[7, 5])
    else:
        fig, ax = plt.subplots(figsize=(6.0, 3.2), dpi=300)

    dot_panel(ax, adult, [r[0] for r in ROWS_ADULT])
    n = len(ROWS_ADULT)
    ax.axhline(n - 2.5, color=MUTED, lw=0.6, ls=(0, (2, 2)))
    ax.text(1.02, n - 1.5, "ciliated", transform=ax.get_yaxis_transform(),
            color=MUTED, va="center")
    ax.text(1.02, (n - 3) / 2, "neuronal", transform=ax.get_yaxis_transform(),
            color=MUTED, va="center")
    ax.set_xlabel("Adult human brain nuclei (nCPM, log scale)", color=MUTED)
    leg = ax.legend(loc="upper center", bbox_to_anchor=(0.45, 1.16 if dev_dir else 1.13), ncol=3,
                    frameon=False, handletextpad=0.3, columnspacing=1.4)
    for t in leg.get_texts():  # gene symbols in italics
        t.set_fontstyle("italic")

    if dev_dir:
        dot_panel(bx, fetal_recs, [r[0] for r in ROWS_FETAL])
        bx.set_xlabel("First-trimester human brain neurons (counts per million, log scale)",
                      color=MUTED)
        for a, letter in ((ax, "a"), (bx, "b")):
            a.text(-0.42, 1.02, letter, transform=a.transAxes, fontweight="bold",
                   va="bottom", color=INK)
        fig.tight_layout(h_pad=1.6)
    else:
        fig.tight_layout()
    # Values below 0.1 are drawn at 0.1; data sources and this note go in the caption.
    for ext in ("png", "pdf"):
        # No creation date in the PDF, so that a rerun gives an identical file
        fig.savefig(f"{out_prefix}.{ext}", bbox_inches="tight", facecolor="white",
                    metadata={"CreationDate": None} if ext == "pdf" else None)


if __name__ == "__main__":
    if len(sys.argv) not in (3, 4):
        sys.exit(__doc__)
    main(*sys.argv[1:])
