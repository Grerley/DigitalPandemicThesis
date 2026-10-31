"""Black-and-white, journal-style rendering of the thesis figures.

Activated by the environment variable DPM_BW=1 (run_all.py calls install() at import). The figure code
of the analysis scripts is left unchanged: the five palette colours are mapped, artist by artist, to
black lines distinguished by line style, grey fills distinguished by shade and hatching, and greyscale
colour maps. Figures are re-sized to the printed text width, typeset in Computer Modern to match the
LaTeX text, and saved as vector PDF (plus a PNG preview) in the directory named by DPM_BW_DIR.
"""
from __future__ import annotations

import os, pathlib, re

import matplotlib
from matplotlib import colors as mcolors

# palette of run_all.py -> (line colour, line style, marker face, fill grey, hatch)
PALETTE = {
    "#1f2933": ("black", "-", "black", "0.20", None),   # INK
    "#2f6db5": ("black", "--", "white", "0.45", "///"),  # C1
    "#c2571a": ("black", "-.", "0.55", "0.70", ".."),   # C2
    "#4b8b3b": ("black", ":", "white", "0.85", "xx"),   # C3
    "#9aa5b1": ("0.55", None, "0.55", "0.80", None),    # GREY
}
# a plain circle marker is replaced by a shape that identifies the series
MARKERS = {"#1f2933": "o", "#2f6db5": "s", "#c2571a": "^", "#4b8b3b": "D", "#9aa5b1": "v"}
# matplotlib default cycle, in case a script relies on it
DEFAULT_CYCLE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd"]
CYCLE_MAP = dict(zip(DEFAULT_CYCLE, ["#1f2933", "#2f6db5", "#c2571a", "#4b8b3b", "#9aa5b1"]))
TEXT_WIDTH_IN = 6.0
LW_SCALE = 0.7
GREYS_LIGHT = mcolors.LinearSegmentedColormap.from_list("greys_light", ["#ffffff", "#8c8c8c"])
LINE_FUNCS = ("plot", "step", "semilogx", "semilogy", "loglog", "errorbar")
GREEK = {"α": r"\alpha", "β": r"\beta", "γ": r"\gamma", "δ": r"\delta", "η": r"\eta", "σ": r"\sigma",
         "ν": r"\nu", "λ": r"\lambda", "χ": r"\chi", "ξ": r"\xi", "φ": r"\phi", "ψ": r"\psi", "Δ": r"\Delta",
         "ε": r"\varepsilon", "κ": r"\kappa", "π": r"\pi", "ρ": r"\rho", "θ": r"\theta", "μ": r"\mu"}
SUB = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")


def _key(c):
    try:
        h = mcolors.to_hex(c).lower()
    except (ValueError, TypeError):
        return None
    return CYCLE_MAP.get(h, h)


def _tex(s: str) -> str:
    """Turn Unicode Greek letters, subscripts and a few symbols into mathtext, outside existing $...$."""
    if not isinstance(s, str) or "$" in s or not s:
        return s
    s = re.sub(r"R₀", r"$\\mathcal{R}_0$", s)
    s = re.sub(r"R_net", r"$\\mathcal{R}_{\\mathrm{net}}$", s)

    def gr(m):
        g, sub = m.group(1), m.group(2)
        body = GREEK[g] + (("_{" + sub.translate(SUB) + "}") if sub else "")
        return "$" + body + "$"
    s = re.sub("([" + "".join(GREEK) + "])([₀-₉]+)?", gr, s)
    s = s.replace("²", "$^2$").replace("→", r"$\rightarrow$").replace("≥", r"$\geq$").replace("≤", r"$\leq$")
    s = s.replace("×", r"$\times$").replace("±", r"$\pm$")
    s = re.sub(r"\$\s*\$", "", s) if "$$" in s else s
    return s


def _style_line(line):
    k = _key(line.get_color())
    if k in PALETTE:
        col, ls, mfc, _, _ = PALETTE[k]
        line.set_color(col)
        if ls is not None and line.get_linestyle() in ("-", "solid"):
            line.set_linestyle(ls)
        if line.get_marker() == "o":
            line.set_marker(MARKERS[k])
        if line.get_marker() not in (None, "None", "", " "):
            line.set_markerfacecolor(mfc); line.set_markeredgecolor("black" if col == "black" else col)


def _style_patch(p, hatch_ok=True):
    fk = _key(p.get_facecolor()); ek = _key(p.get_edgecolor())
    if fk in PALETTE and p.get_fill():
        _, _, _, grey, hatch = PALETTE[fk]
        p.set_facecolor(grey)
        if hatch and hatch_ok:
            p.set_hatch(hatch); p.set_edgecolor("black"); p.set_linewidth(0.4)
    if ek in PALETTE:
        col, ls, _, _, _ = PALETTE[ek]
        p.set_edgecolor(col)
        if not p.get_fill() and ls is not None and p.get_linestyle() in ("-", "solid"):
            p.set_linestyle(ls)


def _legend_rect(fig):
    """Reserve room for figure-level legends above or below the axes when re-laying out."""
    if not fig.legends:
        return (0, 0, 1, 1)
    fig.canvas.draw()
    inv = fig.transFigure.inverted()
    lo, hi = 0.0, 1.0
    for lg in fig.legends:
        bb = lg.get_window_extent().transformed(inv)
        if bb.y1 < 0.5:
            lo = max(lo, bb.y1 + 0.01)
        else:
            hi = min(hi, bb.y0 - 0.01)
    return (0, lo, 1, hi)


def install(plt):
    import matplotlib.axes as maxes
    import matplotlib.figure as mfig
    import matplotlib.text as mtext
    from cycler import cycler
    outdir = pathlib.Path(os.environ.get("DPM_BW_DIR", "figures_bw")); outdir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "font.family": ["cmr10", "DejaVu Serif"], "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
        "font.size": 8.5, "axes.titlesize": 8.5, "axes.labelsize": 8.5, "legend.fontsize": 7.5,
        "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
        "axes.spines.top": True, "axes.spines.right": True, "axes.linewidth": 0.6, "axes.edgecolor": "black",
        "axes.grid": True, "grid.color": "0.88", "grid.linewidth": 0.4, "grid.linestyle": "-",
        "xtick.direction": "out", "ytick.direction": "out", "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "lines.linewidth": 1.1, "lines.markersize": 4, "patch.linewidth": 0.5,
        "legend.frameon": False, "image.cmap": "Greys", "pdf.fonttype": 42, "savefig.dpi": 300,
        "axes.prop_cycle": cycler(color=["#1f2933", "#2f6db5", "#c2571a", "#4b8b3b", "#9aa5b1"]),
        "hatch.linewidth": 0.4, "axes.titleweight": "normal",
    })

    def wrap(name, post):
        orig = getattr(maxes.Axes, name)

        def f(self, *a, **k):
            if name in ("contourf", "imshow", "pcolormesh"):
                k["cmap"] = GREYS_LIGHT     # light greys, so that black lines drawn on top stay legible
            explicit = ("ls" in k or "linestyle" in k or
                        (isinstance(k.get("fmt"), str) and re.search(r"--|-\.|:|-", k["fmt"])) or
                        any(isinstance(x, str) and re.search(r"--|-\.|:|-", x) for x in a))
            out = orig(self, *a, **k)
            try:
                if name in LINE_FUNCS:
                    post(out, explicit)
                else:
                    post(out)
            except Exception:
                pass
            return out
        setattr(maxes.Axes, name, f)

    def lines_post(out, explicit=False):
        for ln in (out if isinstance(out, list) else [out]):
            ls = ln.get_linestyle()
            _style_line(ln)
            if explicit:
                ln.set_linestyle(ls)      # a line style chosen in the script is kept
            ln.set_linewidth(ln.get_linewidth() * LW_SCALE)

    def hist_post(out):
        patches = out[2]
        seq = patches if isinstance(patches, (list, tuple)) else [patches]
        for grp in seq:
            for p in (grp if hasattr(grp, "__iter__") else [grp]):
                _style_patch(p, hatch_ok=False); p.set_hatch(None)
                if p.get_fill():
                    p.set_edgecolor("none")

    def bars_post(out):
        for p in out.patches:
            _style_patch(p)

    def coll_post(out):
        fc = out.get_facecolor()
        if len(fc):
            k = _key(fc[0])
            if k in PALETTE:
                a = fc[0][3]
                g = float(PALETTE[k][3])
                out.set_facecolor((g, g, g, max(a, 0.25)))
        ec = out.get_edgecolor()
        if len(ec):
            k = _key(ec[0])
            if k in PALETTE:
                out.set_edgecolor(PALETTE[k][0])

    def text_post(out):
        k = _key(out.get_color())
        if k in PALETTE:
            out.set_color("black" if PALETTE[k][0] == "black" else PALETTE[k][0])

    for nm in ("plot", "step", "semilogx", "semilogy", "loglog"):
        wrap(nm, lines_post)
    def errorbar_post(out, explicit=False):
        data, caps, bars = out.lines
        k = _key(data.get_color()) if data is not None else None
        col = PALETTE[k][0] if k in PALETTE else "black"
        if data is not None:
            ls = data.get_linestyle(); _style_line(data)
            if explicit:
                data.set_linestyle(ls)
            data.set_linewidth(data.get_linewidth() * LW_SCALE)
        for c in caps:
            c.set_color(col); c.set_markeredgecolor(col); c.set_markerfacecolor(col)
        for b in bars:
            b.set_color(col); b.set_linewidth(0.8)
    wrap("errorbar", errorbar_post)

    def refline_post(ln):
        k = _key(ln.get_color())
        if k in PALETTE:
            ln.set_color(PALETTE[k][0])
    wrap("axhline", refline_post); wrap("axvline", refline_post)
    def field_post(out):
        try:
            out.axes.grid(False)
        except Exception:
            pass
    for nm in ("contourf", "imshow", "pcolormesh"):
        wrap(nm, field_post)
    wrap("hist", hist_post)
    wrap("bar", bars_post); wrap("barh", bars_post)
    wrap("fill_between", coll_post); wrap("fill_betweenx", coll_post)
    wrap("scatter", coll_post)
    wrap("text", text_post); wrap("annotate", text_post)
    wrap("axhspan", lambda p: _style_patch(p, hatch_ok=False)); wrap("axvspan", lambda p: _style_patch(p, hatch_ok=False))

    orig_set_text = mtext.Text.set_text
    mtext.Text.set_text = lambda self, s: orig_set_text(self, _tex(s) if isinstance(s, str) else s)

    orig_save = mfig.Figure.savefig

    def savefig(self, fname, *a, **k):
        stem = pathlib.Path(str(fname)).stem
        w, h = self.get_size_inches()
        if w > TEXT_WIDTH_IN:
            self.set_size_inches(TEXT_WIDTH_IN, h * TEXT_WIDTH_IN / w)
        try:
            if not self.get_constrained_layout():
                self.tight_layout(rect=_legend_rect(self))
        except Exception:
            pass
        k.pop("dpi", None)
        orig_save(self, outdir / f"{stem}.pdf", *a, **k)
        orig_save(self, outdir / f"{stem}.png", *a, dpi=200, **k)
    mfig.Figure.savefig = savefig
