"""Shared matplotlib style for the figures: Garamond, 8 pt, final print size.

The AJUR style guide requires Garamond 8 pt for text inside figures. Figures
are drawn at their final width (text width 7.0 in) and included without
scaling, so 8 pt in the figure is 8 pt on the page.

Font, in order of preference:
  1. TrueType versions of garamondx (URW Garamond No. 8, "NewG8"), the font of
     the official AJUR LaTeX template, made by convert_garamondx_ttf.py and
     placed in the directory given by the environment variable FIGURE_FONT_DIR
     (default: fonts_local/garamondx in the repository root). The fonts are not
     distributed with this repository.
  2. Any Garamond family installed on the system (e.g. EB Garamond).
  3. DejaVu Serif, which ships with matplotlib.
The figures are otherwise identical; only the typeface changes.
"""
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
from matplotlib import font_manager  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

TEXT_WIDTH_IN = 7.0
REPO_ROOT = Path(__file__).resolve().parent.parent
GARAMONDX_FILES = ("NewG8-Regular.ttf", "NewG8-Italic.ttf", "NewG8-Bold.ttf",
                   "NewG8-BoldItalic.ttf")
SYSTEM_GARAMONDS = ("EB Garamond", "Garamond", "Adobe Garamond Pro", "URW Garamond")


def _font_family():
    font_dir = Path(os.environ.get("FIGURE_FONT_DIR", REPO_ROOT / "fonts_local" / "garamondx"))
    files = [font_dir / f for f in GARAMONDX_FILES]
    if all(f.exists() for f in files):
        for f in files:
            font_manager.fontManager.addfont(str(f))
        return "NewG8"
    installed = {f.name for f in font_manager.fontManager.ttflist}
    for family in SYSTEM_GARAMONDS:
        if family in installed:
            return family
    return "DejaVu Serif"


def apply():
    family = _font_family()
    plt.rcParams.update({
        "font.family": family,
        "font.size": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        # Text in $...$ uses the same family; symbols it lacks fall back to
        # Computer Modern, as in the LaTeX body text.
        "mathtext.fontset": "custom",
        "mathtext.rm": family,
        "mathtext.it": f"{family}:italic",
        "mathtext.bf": f"{family}:bold",
        "mathtext.fallback": "cm",
        "pdf.fonttype": 42,
    })
    return family
