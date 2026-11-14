"""Convert the garamondx Type 1 fonts (URW Garamond No. 8, "NewG8") to TrueType.

The manuscript is typeset with garamondx, the font of the official AJUR LaTeX
template. Matplotlib cannot read Type 1 (.pfb) fonts, so the figures would
otherwise use a different Garamond. This script converts the four text faces
(Regular, Italic, Bold, Bold Italic) to TrueType files (cubic outlines
converted to quadratic, max error 1 unit) with a Unicode cmap and the AFM
kerning pairs, for use by ajur_figure_style.py only. TrueType is used so that
matplotlib embeds the fonts in PDFs as ordinary TrueType (pdf.fonttype 42).

The fonts are licensed under the Aladdin Free Public License; the converted
files are for local figure production and are not redistributed.

Usage:
    python scripts/convert_garamondx_ttf.py <dir with NewG8-*.pfb and .afm> fonts_local/garamondx
"""
import sys
from pathlib import Path

from fontTools import afmLib, agl
from fontTools.fontBuilder import FontBuilder
from fontTools.pens.cu2quPen import Cu2QuPen
from fontTools.pens.ttGlyphPen import TTGlyphPen
from fontTools.t1Lib import T1Font
from fontTools.ttLib import newTable
from fontTools.ttLib.tables._k_e_r_n import KernTable_format_0

FACES = {  # file stem: (style name, weight class, italic, bold)
    "NewG8-Reg": ("Regular", 400, False, False),
    "NewG8-Ita": ("Italic", 400, True, False),
    "NewG8-Bol": ("Bold", 700, False, True),
    "NewG8-BolIta": ("Bold Italic", 700, True, True),
}
FAMILY = "NewG8"


def unicode_of(name):
    u = agl.toUnicode(name)
    return ord(u) if len(u) == 1 else None


def convert(pfb, afm, out_dir, style, weight, italic, bold):
    t1 = T1Font(str(pfb))
    t1.parse()
    glyphs = t1.getGlyphSet()
    order = [".notdef"] + sorted(n for n in glyphs.keys() if n != ".notdef")

    tt_glyphs, widths, cmap = {}, {}, {}
    for name in order:
        g = glyphs[name]
        tt_pen = TTGlyphPen(None)
        # PostScript contours run counter-clockwise; TrueType wants clockwise.
        g.draw(Cu2QuPen(tt_pen, max_err=1.0, reverse_direction=True))
        tt_glyphs[name] = tt_pen.glyph()
        widths[name] = round(g.width)
        u = unicode_of(name)
        if u is not None and u not in cmap:
            cmap[u] = name
    # Unicode minus and en-dash fallbacks used by matplotlib labels
    if 0x2212 not in cmap and "minus" in glyphs:
        cmap[0x2212] = "minus"

    info = t1.font["FontInfo"]
    ps_name = f"{FAMILY}-{style.replace(' ', '')}"
    fb = FontBuilder(1000, isTTF=True)
    fb.setupGlyphOrder(order)
    fb.setupCharacterMap(cmap)
    fb.setupGlyf(tt_glyphs)
    # left side bearings from the glyph bounds
    from fontTools.pens.boundsPen import BoundsPen
    metrics = {}
    for name in order:
        bp = BoundsPen(glyphs)
        glyphs[name].draw(bp)
        lsb = round(bp.bounds[0]) if bp.bounds else 0
        metrics[name] = (widths[name], lsb)
    fb.setupHorizontalMetrics(metrics)
    a = afmLib.AFM(str(afm))
    ascent = int(getattr(a, "Ascender", 750))
    descent = int(getattr(a, "Descender", -250))
    fb.setupHorizontalHeader(ascent=ascent, descent=descent)
    fb.setupNameTable({"familyName": FAMILY, "styleName": style,
                       "psName": ps_name, "fullName": f"{FAMILY} {style}",
                       "copyright": str(info.get("Notice", ""))})
    fs_sel = (0x01 if italic else 0) | (0x20 if bold else 0) or 0x40
    fb.setupOS2(usWeightClass=weight, fsSelection=fs_sel,
                sTypoAscender=ascent, sTypoDescender=descent,
                usWinAscent=ascent, usWinDescent=-descent)
    fb.setupPost(isFixedPitch=0, italicAngle=float(info.get("ItalicAngle", 0)))
    fb.font["head"].macStyle = (1 if bold else 0) | (2 if italic else 0)

    kern = newTable("kern")
    kern.version = 0
    sub = KernTable_format_0()
    sub.version, sub.coverage, sub.format = 0, 1, 0
    sub.kernTable = {pair: int(v) for pair, v in a._kerning.items()
                     if pair[0] in glyphs and pair[1] in glyphs}
    kern.kernTables = [sub]
    fb.font["kern"] = kern

    out = Path(out_dir) / f"{ps_name}.ttf"
    fb.save(str(out))
    return out, len(cmap), len(sub.kernTable)


def main(src, out_dir):
    src, out_dir = Path(src), Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    for stem, (style, weight, italic, bold) in FACES.items():
        out, n_cmap, n_kern = convert(src / f"{stem}.pfb", src / f"{stem}.afm", out_dir,
                                      style, weight, italic, bold)
        print(f"{out.name}: {n_cmap} mapped characters, {n_kern} kerning pairs")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
