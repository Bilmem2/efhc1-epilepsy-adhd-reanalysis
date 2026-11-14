"""Attrition chart of the variant-prioritization funnel (manuscript Figure 1).

Reads funnel_summary.json written by filter_opencravat_export.py.
Usage: python plot_prioritization_funnel.py <funnel_summary.json> <out_prefix>
"""
import json
import sys

import ajur_figure_style
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter  # noqa: E402

BLUE, BLUE_DARK = "#2a78d6", "#1c5aa6"
INK, MUTED, GRID = "#0b0b0b", "#52514e", "#e4e3df"
MANUAL_FINAL = 6  # candidates remaining after manual review (Table 2)


def main(summary_path, out_prefix):
    s = json.load(open(summary_path, encoding="utf-8"))
    steps = [
        ("Annotated variant alleles", s["input_rows"]),
        ("Protein-altering or splice-site", s["protein_altering_or_splice"]),
        ("Not benign in ClinVar", s["not_benign_in_ClinVar"]),
        ("Pass quality filters", s["after_QC"]),
        (r"Heterozygous, frequency $\leq$ 0.1%", s["model_A_dominant_het_AF_le_0.001"]),
        ("Gene with known disease (HPO)", s["model_A_disease_gene"]),
        ("Neurodevelopmental phenotype", s["model_A_disease_gene_neuro_terms"]),
        ("Dominant or X-linked inheritance", s["model_A_neuro_and_inheritance_fits"]),
        ("Unresolved after manual review", MANUAL_FINAL),
    ]
    ajur_figure_style.apply()
    fig, ax = plt.subplots(figsize=(6.0, 3.0), dpi=300)
    n = len(steps)
    for i, (label, count) in enumerate(steps):
        y = n - 1 - i
        color = BLUE_DARK if i == n - 1 else BLUE
        ax.barh(y, count, height=0.62, color=color, left=1, zorder=2)
        ax.text((count + 1) * 1.22, y, f"{count:,}", va="center", ha="left", color=INK,
                fontweight="bold" if i == n - 1 else "normal")
    # The step at which EFHC1 p.Glu604Asp was removed and the recessive model are
    # described in the caption.
    ax.set_xscale("log")
    ax.set_xlim(1, 3e5)
    ax.xaxis.set_major_locator(FixedLocator([1, 10, 100, 1000, 10000, 100000]))
    ax.set_xticklabels(["1", "10", "100", "1,000", "10,000", "100,000"])
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(range(n))
    ax.set_yticklabels([lab for lab, _ in reversed(steps)], color=INK)
    ax.set_xlabel("Variants remaining (log scale)", color=MUTED)
    ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", colors=MUTED)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        # No creation date in the PDF, so that a rerun gives an identical file
        fig.savefig(f"{out_prefix}.{ext}", bbox_inches="tight", facecolor="white",
                    metadata={"CreationDate": None} if ext == "pdf" else None)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
