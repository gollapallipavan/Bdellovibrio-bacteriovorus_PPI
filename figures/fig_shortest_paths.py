"""Shortest-path-length heatmap: directed predator -> prey hop counts."""
import numpy as np
import networkx as nx
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
import plot_utils as pu

G = pu.build_graph()
PRED = [n for n in G if G.nodes[n]["role"] == "predator"]
PREY = sorted(n for n in G if G.nodes[n]["role"] == "prey")
PRED = sorted(PRED)
sp = dict(nx.all_pairs_shortest_path_length(G))
M = np.full((len(PRED), len(PREY)), np.nan)
for i, u in enumerate(PRED):
    for j, v in enumerate(PREY):
        d = sp.get(u, {}).get(v)
        if d is not None:
            M[i, j] = d

cmap = mpl.cm.RdYlGn_r.copy(); cmap.set_bad("#d9d9d9")
norm = Normalize(vmin=1, vmax=5)
fig = plt.figure(figsize=(7.4, 8.8))
ax = fig.add_axes([0.14, 0.11, 0.74, 0.845])
im = ax.imshow(np.ma.masked_invalid(M), cmap=cmap, norm=norm, aspect="auto")
for i in range(len(PRED)):
    for j in range(len(PREY)):
        d = M[i, j]
        if np.isnan(d):
            ax.text(j, i, "∞", ha="center", va="center", fontsize=11, fontweight="bold")
        else:
            n = norm(min(d, 5))
            ax.text(j, i, str(int(d)), ha="center", va="center", fontsize=11,
                    color="white" if (n <= 0.2 or n >= 0.80) else "black")
ax.set_xticks(range(len(PREY))); ax.set_xticklabels(PREY, rotation=45, ha="right", rotation_mode="anchor", fontsize=11.5, fontweight="bold")
ax.set_yticks(range(len(PRED))); ax.set_yticklabels(PRED, fontsize=11.5, fontweight="bold")
for t in ax.get_xticklabels(): t.set_color(pu.ST_COL)
for t in ax.get_yticklabels(): t.set_color(pu.BB_COL)
ax.set_xticks(np.arange(-0.5, len(PREY), 1), minor=True); ax.set_yticks(np.arange(-0.5, len(PRED), 1), minor=True)
ax.grid(which="minor", color="white", lw=1.4); ax.tick_params(which="minor", length=0)
ax.set_xlabel("S. Typhimurium proteins (prey)", fontsize=13, labelpad=8, fontstyle="italic")
ax.set_ylabel("B. bacteriovorus proteins (predator)", fontsize=13, labelpad=8, fontstyle="italic")
cax = fig.add_axes([0.895, 0.11, 0.028, 0.845])
cb = fig.colorbar(im, cax=cax); cb.set_label("Shortest path length (hops)", fontsize=12.5, labelpad=8)
cb.set_ticks([1, 2, 3, 4, 5]); cb.ax.tick_params(labelsize=11)
ax.legend(handles=[Patch(facecolor="#d9d9d9", edgecolor="black", label="∞  unreachable")],
          loc="lower left", bbox_to_anchor=(0.0, 1.005), fontsize=11, frameon=True, edgecolor="#888888", handlelength=1.4)
pu.save(fig, "fig_shortest_paths")
