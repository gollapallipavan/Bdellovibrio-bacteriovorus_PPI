"""Shared plotting helpers for the PPI network figures.

Publication defaults: Times New Roman (falls back to a serif face if Times New
Roman is not installed), white background, 600-dpi TIFF + PDF output.

To use Times New Roman, drop the TTF files (times.ttf, timesbd.ttf, timesi.ttf,
timesbi.ttf) into figures/fonts/ or install them system-wide; otherwise the
figures render in the default serif face.
"""
import os
import io
import csv
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Circle, FancyArrowPatch
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(ROOT, "outputs")
os.makedirs(OUT, exist_ok=True)

# ---- fonts -----------------------------------------------------------------
_fontdir = os.path.join(HERE, "fonts")
if os.path.isdir(_fontdir):
    for fn in os.listdir(_fontdir):
        if fn.lower().endswith(".ttf"):
            try:
                font_manager.fontManager.addfont(os.path.join(_fontdir, fn))
            except Exception:
                pass
_have_times = any("Times New Roman" == f.name for f in font_manager.fontManager.ttflist)
FAMILY = "Times New Roman" if _have_times else "serif"
if not _have_times:
    print("[plot_utils] Times New Roman not found; using the default serif face. "
          "Add the TTFs to figures/fonts/ to match the manuscript exactly.")

plt.rcParams.update({
    "font.family": FAMILY, "mathtext.fontset": "cm",
    "font.size": 11, "axes.labelsize": 12, "axes.titlesize": 12,
    "xtick.labelsize": 11, "ytick.labelsize": 11, "legend.fontsize": 10.5,
    "axes.linewidth": 1.0, "figure.facecolor": "white", "axes.facecolor": "white",
    "savefig.facecolor": "white", "pdf.fonttype": 42, "ps.fonttype": 42,
})

# ---- palette ---------------------------------------------------------------
BB_COL, ST_COL = "#1e7b34", "#b03a2e"
CATEGORY_COLOR = {
    "Recognition": "#3498db", "Motility": "#3498db", "Hydrolytic": "#e67e22",
    "Regulatory": "#f1c40f", "Chaperone": "#5b2c6f", "Cell_division": "#5b2c6f",
    "T6SS": "#a569bd", "Structural_OMP": "#ec7063", "Virulence": "#922b21",
    "Cell_wall": "#922b21", "Biofilm": "#1a5276", "AMR": "#839192",
}
EDGE_STYLE = {
    "cross_organism": dict(color="#34495e", lw=1.5, label="Cross-organism (Bb → St)"),
    "intra_predator": dict(color="#1e8449", lw=1.3, label="Intra-predator"),
    "intra_prey": dict(color="#cb4335", lw=1.3, label="Intra-prey"),
}


# ---- data loading ----------------------------------------------------------
def load_nodes():
    with open(os.path.join(DATA, "ppi_nodes.csv"), newline="") as f:
        return list(csv.DictReader(f))


def load_edges():
    with open(os.path.join(DATA, "ppi_edges.csv"), newline="") as f:
        return list(csv.DictReader(f))


def build_graph():
    import networkx as nx
    G = nx.DiGraph()
    for r in load_nodes():
        G.add_node(r["node_id"], **r)
    for r in load_edges():
        G.add_edge(r["source"], r["target"], **r)
    return G


# ---- output ----------------------------------------------------------------
def save(fig, name, dpi=600):
    """Save a 600-dpi RGB TIFF (LZW) and a vector PDF to outputs/."""
    buf = io.BytesIO()
    fig.savefig(buf, dpi=dpi, format="png", facecolor="white", bbox_inches="tight", pad_inches=0.08)
    buf.seek(0)
    im = Image.open(buf).convert("RGB")
    im.save(os.path.join(OUT, name + ".tif"), dpi=(dpi, dpi), compression="tiff_lzw")
    fig.savefig(os.path.join(OUT, name + ".pdf"), facecolor="white", bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print(f"  wrote outputs/{name}.tif and .pdf  ({im.size[0]}x{im.size[1]} px)")


# ---- geometry helpers ------------------------------------------------------
def _lum(hexc):
    r, g, b = matplotlib.colors.to_rgb(hexc)
    return 0.299 * r + 0.587 * g + 0.114 * b


def text_color(hexc):
    return "black" if _lum(hexc) > 0.55 else "white"


def text_width_pt(texts, fontsize, weight="bold"):
    fig = plt.figure(figsize=(4, 4), dpi=72)
    r = fig.canvas.get_renderer()
    out = {}
    for t in texts:
        tt = fig.text(0, 0, t, fontsize=fontsize, fontweight=weight)
        out[t] = tt.get_window_extent(renderer=r).width
        tt.remove()
    plt.close(fig)
    return out


def arc_points(p1, p2, rad, n=40):
    p1 = np.asarray(p1, float); p2 = np.asarray(p2, float); d = p2 - p1
    c = (p1 + p2) / 2 + rad * np.array([d[1], -d[0]])
    t = np.linspace(0.0, 1.0, n)[:, None]
    return (1 - t) ** 2 * p1 + 2 * (1 - t) * t * c + t ** 2 * p2


def edge_rads(edges, rad_single, rad_pair, pos=None, bulge=0.2):
    es = set(edges); out = {}
    for e in edges:
        if (e[1], e[0]) in es and pos is not None:
            L = np.linalg.norm(np.subtract(pos[e[1]], pos[e[0]]))
            out[e] = min(rad_pair, 2 * bulge / max(L, 1e-6))
        elif (e[1], e[0]) in es:
            out[e] = rad_pair
        else:
            out[e] = rad_single
    return out


def draw_network(ax, pos, edges, edge_kind, radius, node_color, fontsize=10.5, rads=None):
    """edge_kind(u,v) -> one of EDGE_STYLE keys."""
    for (u, v) in edges:
        st = EDGE_STYLE[edge_kind(u, v)]
        rad = rads[(u, v)] if rads is not None else 0.0
        p1, p2 = np.array(pos[u], float), np.array(pos[v], float)
        d = p2 - p1; L = np.linalg.norm(d); ang = np.arctan(2 * rad) if rad else 0.0

        def rot(vec, a):
            c, s = np.cos(a), np.sin(a)
            return np.array([c * vec[0] - s * vec[1], s * vec[0] + c * vec[1]])
        u_dir = d / L
        s = p1 + rot(u_dir, ang) * radius[u]
        e = p2 + rot(-u_dir, -ang) * (radius[v] + 0.01)
        ax.add_patch(FancyArrowPatch(
            s, e, connectionstyle=f"arc3,rad={rad}",
            arrowstyle="-|>,head_length=5,head_width=2.6", color=st["color"],
            lw=st["lw"], alpha=0.9, shrinkA=0, shrinkB=0, zorder=1, mutation_scale=1.6))
    for n, (x, y) in pos.items():
        c = node_color[n]
        ax.add_patch(Circle((x, y), radius[n], facecolor=c, edgecolor="black", lw=1.2, zorder=3))
        ax.text(x, y, n, ha="center", va="center", fontsize=fontsize, fontweight="bold",
                color=text_color(c), zorder=4)
