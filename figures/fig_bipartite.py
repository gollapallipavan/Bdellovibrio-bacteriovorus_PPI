"""Bipartite cross-organism map: predator -> prey edges coloured by the
reconstructed interaction mechanism (from data/ppi_edges.csv).

Carbohydrate/lipid/nuclease mechanisms that had been assigned to PROTEIN targets
are shown as 'Predicted (unverified mechanism)'."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle, Ellipse, FancyArrowPatch
import plot_utils as pu

G = pu.build_graph()
BB = {n for n in G if G.nodes[n]["role"] == "predator"}
EMECH = {(u, v): G.edges[u, v]["mechanism_reconstructed"]
         for u, v in G.edges() if (u in BB) != (v in BB)}

MECH = {
    "Physical binding": ("#1f4e79", "-"), "Mechanical contact": ("#2e86de", "-"),
    "Motility mediated": ("#48c9e0", (0, (5, 2))), "T6SS effector delivery": ("#8e44ad", "-"),
    "Proteolytic degradation": ("#922b21", "-"), "Proteolytic disruption": ("#c0392b", (0, (5, 2))),
    "Predicted (unverified mechanism)": ("#7f8c8d", (0, (2, 2))),
}
LEFT = ["VgrG", "PilB", "OapA", "FliF", "Bd3507", "Bd3499", "Bd3279", "Bd2352", "Bd1972", "Bd0934", "Bd0537"]
RIGHT = ["SseB", "SipA", "OmpF", "OmpC", "OmpA", "MurA", "FimC", "BcsA", "AmpC"]
LEFT = [n for n in LEFT if n in {u for u, _ in EMECH}]
RIGHT = [n for n in RIGHT if n in {v for _, v in EMECH}]
NODE_COLOR = {n: pu.CATEGORY_COLOR[G.nodes[n]["functional_category"]] for n in G}

fig = plt.figure(figsize=(7.2, 8.6))
ax = fig.add_axes([0, 0.17, 1, 0.83])
ax.set_xlim(0, 10); ax.set_ylim(-0.2, 11.6); ax.axis("off")
xl, xr = 2.3, 7.7
yl = {n: 10.0 - i for i, n in enumerate(LEFT)}
yr = {n: 10.0 - i * 1.12 for i, n in enumerate(RIGHT)}
inc, out = {}, {}
for v in RIGHT:
    srcs = sorted([u for (u, w) in EMECH if w == v], key=lambda u: -yl[u]); k = len(srcs)
    for i, u in enumerate(srcs): inc[(u, v)] = (0.22 - 0.44 * i / (k - 1)) if k > 1 else 0.0
for u in LEFT:
    tg = sorted([w for (x, w) in EMECH if x == u], key=lambda w: -yr[w]); k = len(tg)
    for i, w in enumerate(tg): out[(u, w)] = (0.18 - 0.36 * i / (k - 1)) if k > 1 else 0.0
for (u, v), m in EMECH.items():
    c, ls = MECH[m]
    ax.add_patch(FancyArrowPatch((xl + 0.36, yl[u] + out[(u, v)]), (xr - 0.34, yr[v] + inc[(u, v)]),
                 connectionstyle="arc3,rad=0.06", arrowstyle="-|>,head_length=6,head_width=3",
                 color=c, lw=1.8, linestyle=ls, shrinkA=0, shrinkB=0, zorder=2))
for n in LEFT:
    ax.add_patch(Ellipse((xl, yl[n]), 0.72, 0.6, facecolor=NODE_COLOR[n], edgecolor="black", lw=1.2, zorder=3))
    ax.text(xl - 0.5, yl[n], n, ha="right", va="center", fontsize=12, fontweight="bold", color="#14532d",
            bbox=dict(boxstyle="round,pad=0.25", fc="#eaf6ee", ec="#1e7b34", lw=0.8))
for n in RIGHT:
    ax.add_patch(Rectangle((xr - 0.32, yr[n] - 0.29), 0.64, 0.58, facecolor=NODE_COLOR[n], edgecolor="black", lw=1.2, zorder=3))
    ax.text(xr + 0.5, yr[n], n, ha="left", va="center", fontsize=12, fontweight="bold", color="#7b241c",
            bbox=dict(boxstyle="round,pad=0.25", fc="#fbeceb", ec="#b03a2e", lw=0.8))
ax.text(xl - 0.6, 11.15, "B. bacteriovorus\n(predator proteins)", ha="center", va="center", fontsize=12.5,
        fontweight="bold", color="white", style="italic", bbox=dict(boxstyle="round,pad=0.4", fc=pu.BB_COL, ec="black", lw=0.8))
ax.text(xr + 0.6, 11.15, "S. Typhimurium\n(prey proteins)", ha="center", va="center", fontsize=12.5,
        fontweight="bold", color="white", style="italic", bbox=dict(boxstyle="round,pad=0.4", fc=pu.ST_COL, ec="black", lw=0.8))
h = [Line2D([0], [0], color=c, lw=2.6, ls=ls, label=m) for m, (c, ls) in MECH.items()]
fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=10.5, frameon=True,
           title="Interaction mechanism (carbohydrate/lipid/nuclease labels on protein\ntargets shown as predicted — see text)",
           title_fontproperties={"weight": "bold", "size": 10.5}, handlelength=2.6, columnspacing=1.6, edgecolor="#888888")
pu.save(fig, "fig_bipartite")
