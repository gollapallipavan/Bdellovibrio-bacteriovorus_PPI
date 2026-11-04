"""Full 25-protein predator-prey PPI network.
Node layout is read from data/network_layout.json (a force-refined spring layout;
delete that file to recompute a fresh layout)."""
import os, json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import plot_utils as pu

G = pu.build_graph()
NODES = list(G.nodes()); EDGES = list(G.edges())
BB = {n for n in G if G.nodes[n]["role"] == "predator"}
deg = dict(G.degree())
node_color = {n: pu.CATEGORY_COLOR[G.nodes[n]["functional_category"]] for n in NODES}

# legend groups (colour -> label)
CAT_LABELS = [
    ("#3498db", "Bb: Recognition / Motility"), ("#e67e22", "Bb: Hydrolytic enzymes"),
    ("#f1c40f", "Bb: Regulatory"), ("#5b2c6f", "Bb: Chaperone / Division"), ("#a569bd", "Bb: T6SS"),
    ("#ec7063", "St: Structural OMPs"), ("#922b21", "St: Virulence"),
    ("#1a5276", "St: Biofilm"), ("#839192", "St: AMR"),
]

layout = os.path.join(pu.DATA, "network_layout.json")
if os.path.exists(layout):
    pos = {k: tuple(v) for k, v in json.load(open(layout)).items() if k in NODES}
else:  # recompute with a simple spring layout as a fallback
    import networkx as nx
    pos = nx.spring_layout(G, seed=42, k=2.0)

FS = 10.5
W = pu.text_width_pt(NODES, FS)
XMIN, XMAX, FIGW = -1.2, 10.9, 7.2
unit_pt = FIGW * 72 / (XMAX - XMIN)
radius = {n: max(0.26 + 0.025 * deg[n], (W[n] / unit_pt) / 2 + 0.07) for n in NODES}
RADS = pu.edge_rads(EDGES, 0.06, 0.2, pos, bulge=0.2)


def kind(u, v):
    if u in BB and v in BB: return "intra_predator"
    if u not in BB and v not in BB: return "intra_prey"
    return "cross_organism"


YMIN, YMAX = -3.0, 10.1
FIGH = FIGW * (YMAX - YMIN) / (XMAX - XMIN) + 1.35
fig = plt.figure(figsize=(FIGW, FIGH))
axh = FIGW * (YMAX - YMIN) / (XMAX - XMIN) / FIGH
ax = fig.add_axes([0, 0, 1, axh])
ax.set_xlim(XMIN, XMAX); ax.set_ylim(YMIN, YMAX); ax.set_aspect("equal"); ax.axis("off")
pu.draw_network(ax, pos, EDGES, kind, radius, node_color, fontsize=FS, rads=RADS)
handles = [Patch(facecolor=c, edgecolor="black", label=l) for c, l in CAT_LABELS]
handles += [Line2D([0], [0], color=pu.EDGE_STYLE[k]["color"], lw=2.2, label=pu.EDGE_STYLE[k]["label"])
            for k in ("cross_organism", "intra_predator", "intra_prey")]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=True,
           fontsize=10.5, title="Node / Edge Category", title_fontproperties={"weight": "bold", "size": 11.5},
           handlelength=1.8, columnspacing=1.2, borderpad=0.5, edgecolor="#888888")
pu.save(fig, "fig_network")
