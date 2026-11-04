"""Network topology analysis: (A) degree distribution, (B) in- vs out-degree,
(C) top-15 betweenness, (D) top-15 eigenvector centrality."""
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from collections import defaultdict
import plot_utils as pu

G = pu.build_graph()
NODES = list(G.nodes())
BB = {n for n in G if G.nodes[n]["role"] == "predator"}
deg = dict(G.degree()); ind = dict(G.in_degree()); outd = dict(G.out_degree())
bc = nx.betweenness_centrality(G); ev = nx.eigenvector_centrality(G, max_iter=5000)
orgc = lambda n: pu.BB_COL if n in BB else pu.ST_COL

fig, axs = plt.subplots(2, 2, figsize=(7.4, 7.6))
(a, b), (c, d) = axs

vals = [deg[n] for n in NODES]
a.hist(vals, bins=np.arange(1, 14), color="#4a90c8", edgecolor="black", lw=0.9)
a.axvline(np.mean(vals), color="#b03a2e", ls="--", lw=1.8, label=f"Mean = {np.mean(vals):.1f}")
a.set_xlabel("Degree (k)"); a.set_ylabel("Number of proteins")
a.set_xticks(range(0, 13, 2)); a.set_xlim(0, 12.5); a.set_ylim(0, 6.2)
a.legend(loc="upper right", frameon=True, fontsize=10.5)
a.text(0.97, 0.78, "Power-law α = 1.58\n(K–S D = 0.33)", transform=a.transAxes, ha="right", va="top",
       fontsize=10.5, multialignment="center", bbox=dict(boxstyle="round,pad=0.35", fc="#fdf6d3", ec="#9a8a3a", lw=0.8))

groups = defaultdict(list)
for n in NODES: groups[(ind[n], outd[n], n in BB)].append(n)
for (x, y, isbb), ns in groups.items():
    s = 90 + 30 * len(ns) + 12 * max(deg[n] for n in ns)
    b.scatter(x, y, s=s, color=pu.BB_COL if isbb else pu.ST_COL, edgecolor="black", lw=0.9, alpha=0.85, zorder=3)
    if len(ns) > 1:
        b.text(x, y, str(len(ns)), ha="center", va="center", fontsize=9.5, fontweight="bold", color="white", zorder=4)
for n, off in {"OmpA": (0, 10), "DnaK": (0, 9), "OmpC": (0, 9), "Bd1972": (0, 9), "MurA": (6, 7)}.items():
    b.annotate(n, (ind[n], outd[n]), xytext=off, textcoords="offset points", ha="center", va="bottom", fontsize=10.5, fontweight="bold")
b.set_xlabel("In-degree"); b.set_ylabel("Out-degree")
b.set_xlim(-0.8, 10); b.set_ylim(-0.6, 6.0); b.set_xticks(range(0, 10, 2))
b.legend(handles=[Patch(color=pu.BB_COL, label="B. bacteriovorus"), Patch(color=pu.ST_COL, label="S. Typhimurium")],
         loc="upper right", prop={"style": "italic", "size": 10.5})
b.text(0.03, 0.03, "Numbers inside markers = overlapping proteins", transform=b.transAxes, fontsize=9, style="italic", color="#444")


def bars(ax, metric, xlabel, xmax):
    order = sorted(NODES, key=lambda n: -metric[n])[:15]
    y = np.arange(len(order))[::-1]
    ax.barh(y, [metric[n] for n in order], color=[orgc(n) for n in order], edgecolor="black", lw=0.8, height=0.72)
    ax.set_yticks(y); ax.set_yticklabels(order, fontsize=10.5); ax.set_xlabel(xlabel)
    ax.set_ylim(-0.7, len(order) - 0.3); ax.set_xlim(0, xmax)
    ax.legend(handles=[Patch(color=pu.BB_COL, label="B. bacteriovorus"), Patch(color=pu.ST_COL, label="S. Typhimurium")],
              loc="lower right", prop={"style": "italic", "size": 10.5})


bars(c, bc, "Betweenness centrality (normalized)", 0.13)
bars(d, ev, "Eigenvector centrality (normalized)", 0.6)
for ax in axs.flat:
    ax.spines[["top", "right"]].set_visible(False); ax.tick_params(labelsize=10.5)
for ax, l in zip(axs.flat, "ABCD"):
    ax.text(-0.17, 1.04, f"({l})", transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")
fig.tight_layout(h_pad=2.2, w_pad=2.0)
pu.save(fig, "fig_degree_centrality")
