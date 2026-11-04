"""Degree-distribution power-law validation (CCDF, log-log), 25-protein network.
Panel A: restricted fit (x_min=4); Panel B: full-range fit (x_min=1)."""
import numpy as np, warnings
warnings.filterwarnings("ignore")
import powerlaw
from matplotlib.lines import Line2D
import plot_utils as pu
import matplotlib.pyplot as plt

G = pu.build_graph()
deg = np.array([G.degree(n) for n in G])

STYLE = {"power_law": ("#d62728", "-"), "lognormal": ("#1f5fb4", (0, (6, 3))),
         "exponential": ("#2ca02c", (0, (7, 2, 1.5, 2))), "stretched_exponential": ("#7b4fa6", (0, (1.5, 2)))}
LABEL = {"power_law": "Power law", "lognormal": "Log-normal",
         "exponential": "Exponential", "stretched_exponential": "Stretched exponential"}


def emp_ccdf(d, xmin):
    d = np.sort(d[d >= xmin]); ks = np.unique(d)
    return ks, np.array([(d >= k).sum() / len(d) for k in ks])


fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.6, 3.7))


def panel(ax, xmin, dists, letter):
    fit = powerlaw.Fit(deg, xmin=xmin, discrete=True, verbose=False)
    ks, y = emp_ccdf(deg, xmin)
    ax.scatter(ks, y, s=55, color="black", zorder=5)
    for nm in dists:
        c, ls = STYLE[nm]
        getattr(fit, nm).plot_ccdf(ax=ax, color=c, linestyle=ls, lw=2.0)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("Degree (k)", fontsize=12.5); ax.set_ylabel("P(K ≥ k)", fontsize=12.5)
    ax.tick_params(labelsize=10.5, which="both"); ax.spines[["top", "right"]].set_visible(False)
    hs = [Line2D([0], [0], marker="o", color="none", markerfacecolor="black", markersize=8, label="Empirical data")]
    for nm in dists:
        c, ls = STYLE[nm]
        lab = f"Power law (x$_{{min}}$={xmin}, α={fit.power_law.alpha:.2f})" if nm == "power_law" else LABEL[nm]
        hs.append(Line2D([0], [0], color=c, ls=ls, lw=2.0, label=lab))
    ax.legend(handles=hs, loc="lower left", fontsize=9, frameon=True, edgecolor="#999999", handlelength=2.4)
    ax.text(-0.17, 1.03, f"{letter})", transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")


panel(axA, 4, ["power_law", "lognormal", "exponential", "stretched_exponential"], "A")
panel(axB, 1, ["power_law", "lognormal", "exponential"], "B")
fig.tight_layout(w_pad=3.0)
pu.save(fig, "fig_powerlaw_ccdf")
