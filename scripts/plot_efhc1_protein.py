"""EFHC1 protein map and EF-hand close-up (manuscript Figure 2).

Panel a: domain map of EFHC1 (UniProt Q5JVL4) with historical JME-associated
variants (UniProt "in EJM1" natural variants) and the proband's p.Glu604Asp.
Panel b: sequence of residues 577-609 aligned to the canonical 12-residue
EF-hand loop motif (putative loop 587-598), the positions corresponding to the
canonical calcium-ligand positions (1, 3, 5, 7, 9, 12) and the helices flanking
it in the AlphaFold model. Glu604 lies on the following helix and does not
correspond to a canonical ligand position. This is a sequence/motif comparison
on a predicted model, not an experimental map of calcium coordination.

Usage: python plot_efhc1_protein.py <uniprot_Q5JVL4.json> <out_prefix>
"""
import json
import sys

import ajur_figure_style
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle  # noqa: E402

BLUE, BLUE_LIGHT, ORANGE = "#2a78d6", "#cfe0f6", "#eb6834"
INK, MUTED, LINE, GREY = "#0b0b0b", "#52514e", "#b9b8b2", "#8a8984"
LOOP_START = 587           # first residue of the putative loop (motif alignment)
LIGAND_POSITIONS = (1, 3, 5, 7, 9, 12)
HELIX_IN = (579, 586)      # helical backbone angles in the AlphaFold model
HELIX_OUT = (596, 606)
VARIANT = 604
WINDOW = (577, 609)


def main(uniprot_path, out_prefix):
    u = json.load(open(uniprot_path, encoding="utf-8"))
    seq = u["sequence"]["value"]
    length = len(seq)
    domains, jme = [], []
    for f in u["features"]:
        s, e = f["location"]["start"]["value"], f["location"]["end"]["value"]
        if f["type"] == "Domain":
            domains.append((f["description"], s, e))
        if f["type"] == "Natural variant" and "in EJM1" in f.get("description", ""):
            jme.append(s)
    jme = sorted(set(jme))

    ajur_figure_style.apply()
    fig = plt.figure(figsize=(7.0, 3.1), dpi=300)
    ax = fig.add_axes([0.04, 0.55, 0.92, 0.40])
    bx = fig.add_axes([0.04, 0.04, 0.92, 0.42])

    # ---- panel a: domain map
    ax.add_patch(Rectangle((1, -0.08), length - 1, 0.16, color=LINE, lw=0))
    for name, s, e in domains:
        color = ORANGE if name == "EF-hand" else BLUE
        ax.add_patch(Rectangle((s, -0.28), e - s, 0.56, color=color, lw=0))
        label = name.replace("DM10 ", "DM10-")
        ax.text((s + e) / 2, -0.62, label, ha="center", va="top", color=INK)
    for p in jme:
        ax.plot([p, p], [0.3, 0.78], color=GREY, lw=0.7)
        ax.plot(p, 0.82, "o", ms=3.2, color=GREY, mec="white", mew=0.4)
    ax.plot([VARIANT, VARIANT], [0.3, 1.1], color=ORANGE, lw=1.1)
    ax.plot(VARIANT, 1.15, "o", ms=5, color=ORANGE, mec="white", mew=0.6)
    ax.text(VARIANT - 8, 1.15, "p.Glu604Asp (proband)", ha="right", va="center",
            color=ORANGE, fontweight="bold")
    ax.text(jme[0], 1.15, "Historical JME-associated variants", ha="left", va="center",
            color=GREY)
    for t in (1, 100, 200, 300, 400, 500, 600, length):
        ax.text(t, -1.05, str(t), ha="center", va="top", color=MUTED)
    ax.set_xlim(-10, length + 10)
    ax.set_ylim(-1.35, 1.4)
    ax.axis("off")
    ax.text(-10, 1.4, "a", fontweight="bold", va="top")

    # ---- panel b: EF-hand close-up
    start, end = WINDOW
    w = 1.0
    for i, pos in enumerate(range(start, end + 1)):
        x = i * w
        loop_pos = pos - LOOP_START + 1
        face, edge, text = "white", LINE, INK
        if 1 <= loop_pos <= 12 and loop_pos in LIGAND_POSITIONS:
            face, edge = BLUE_LIGHT, BLUE
        if pos == VARIANT:
            face, edge, text = ORANGE, ORANGE, "white"
        bx.add_patch(Rectangle((x + 0.06, 0), w - 0.12, 1, facecolor=face,
                               edgecolor=edge, lw=0.8))
        if loop_pos == 12:
            bx.add_patch(Rectangle((x + 0.06, 0), w - 0.12, 1, facecolor="none",
                                   edgecolor=BLUE, lw=0.8, hatch="////"))
        bx.text(x + w / 2, 0.5, seq[pos - 1], ha="center", va="center", color=text,
                fontweight="bold")
        if 1 <= loop_pos <= 12:
            bx.text(x + w / 2, -0.32, str(loop_pos), ha="center", va="center",
                    color=BLUE if loop_pos in LIGAND_POSITIONS else MUTED)
        if pos % 5 == 0:
            bx.text(x + w / 2, 1.32, str(pos), ha="center", va="center", color=MUTED)

    def span(a, b, y, label, color):
        x0, x1 = (a - start) * w + 0.1, (b - start + 1) * w - 0.1
        bx.plot([x0, x1], [y, y], color=color, lw=2.2, solid_capstyle="round")
        bx.text((x0 + x1) / 2, y + 0.3, label, ha="center", va="center", color=color)

    # Short labels only; shading, hatching and the AlphaFold source are explained in the caption.
    span(*HELIX_IN, 1.85, "helix", GREY)
    span(*HELIX_OUT, 1.85, "helix", GREY)
    span(LOOP_START, LOOP_START + 11, -0.75, "", BLUE)
    lx = (LOOP_START - start + 6) * w
    bx.text(lx, -1.15, "putative EF-hand loop", ha="center", va="center", color=BLUE)
    vx = (VARIANT - start + 0.5) * w
    bx.text(vx, -0.32, "Glu604", ha="center", va="center", color=ORANGE, fontweight="bold")
    bx.set_xlim(-0.3, (end - start + 1) * w + 0.3)
    bx.set_ylim(-1.45, 2.4)
    bx.axis("off")
    bx.text(-0.3, 2.4, "b", fontweight="bold", va="top")

    for ext in ("pdf", "png"):
        # No creation date in the PDF, so that a rerun gives an identical file
        fig.savefig(f"{out_prefix}.{ext}", bbox_inches="tight", facecolor="white",
                    metadata={"CreationDate": None} if ext == "pdf" else None)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
