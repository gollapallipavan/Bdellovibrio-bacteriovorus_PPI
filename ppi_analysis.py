#!/usr/bin/env python3
"""
Topology analysis of the Bdellovibrio bacteriovorus - Salmonella Typhimurium
protein-protein interaction (PPI) network.

Builds the directed network from data/ppi_nodes.csv and data/ppi_edges.csv and
reproduces every value reported in Tables 1-3 of the manuscript:
  - Table 1: per-protein centrality metrics          -> outputs/table1_centrality.csv
  - Table 2: global topology properties              -> outputs/table2_global.csv
  - Table 3: power-law vs alternative-distribution fits -> outputs/table3_powerlaw.csv
  - community detection (greedy modularity)

Network: 25 proteins, 52 directed edges. The non-protein species LPS and DNA/RNA
that were present in an earlier version are NOT part of this protein-protein
interaction network.

Usage:
    python ppi_analysis.py
Requires: networkx>=3, numpy, powerlaw (optional; skipped if not installed).
"""
import os
import csv
import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)


def load_graph():
    """Build a directed graph with node attributes from the CSV files."""
    G = nx.DiGraph()
    with open(os.path.join(DATA, "ppi_nodes.csv"), newline="") as f:
        for r in csv.DictReader(f):
            G.add_node(r["node_id"], **r)
    with open(os.path.join(DATA, "ppi_edges.csv"), newline="") as f:
        for r in csv.DictReader(f):
            G.add_edge(r["source"], r["target"], **r)
    return G


def node_table(G):
    bc = nx.betweenness_centrality(G)
    ev = nx.eigenvector_centrality(G, max_iter=5000)
    rows = []
    for n in G.nodes():
        rows.append(dict(
            protein=n, organism=G.nodes[n]["organism"],
            category=G.nodes[n]["functional_category"],
            degree=G.degree(n), in_degree=G.in_degree(n), out_degree=G.out_degree(n),
            betweenness=round(bc[n], 4), eigenvector=round(ev[n], 4),
            community=G.nodes[n]["community"]))
    rows.sort(key=lambda r: -r["betweenness"])
    with open(os.path.join(OUT, "table1_centrality.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)
    return rows


def global_properties(G):
    BB = [n for n in G if G.nodes[n]["role"] == "predator"]
    ST = [n for n in G if G.nodes[n]["role"] == "prey"]
    E = list(G.edges())
    cross = [e for e in E if (e[0] in BB) != (e[1] in BB)]
    ibb = [e for e in E if e[0] in BB and e[1] in BB]
    ist = [e for e in E if e[0] in ST and e[1] in ST]
    U = G.to_undirected()
    deg = np.array([d for _, d in G.degree()])
    comm = nx.community.greedy_modularity_communities(U)
    props = [
        ("Total nodes (proteins)", G.number_of_nodes()),
        ("Total directed edges", G.number_of_edges()),
        ("Cross-organism edges (Bb->St)", f"{len(cross)} ({100*len(cross)/len(E):.1f}%)"),
        ("Intra-predator edges", f"{len(ibb)} ({100*len(ibb)/len(E):.1f}%)"),
        ("Intra-prey edges", f"{len(ist)} ({100*len(ist)/len(E):.1f}%)"),
        ("Network density", round(nx.density(G), 3)),
        ("Average degree", round(deg.mean(), 2)),
        ("Degree SD", round(deg.std(ddof=1), 2)),
        ("Max degree (hub)", f"{deg.max()} ({max(G.degree(), key=lambda x: x[1])[0]})"),
        ("Average clustering coefficient", round(nx.average_clustering(U), 3)),
        ("Average shortest path length", round(nx.average_shortest_path_length(U), 2)),
        ("Number of communities", len(comm)),
        ("Modularity Q", round(nx.community.modularity(U, comm), 3)),
    ]
    with open(os.path.join(OUT, "table2_global.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["property", "value"]); w.writerows(props)
    return props, comm, deg


def powerlaw_table(deg):
    try:
        import powerlaw
    except ImportError:
        print("\n[powerlaw not installed -> skipping Table 3; `pip install powerlaw`]")
        return
    import warnings
    warnings.filterwarnings("ignore")
    rows = []
    for xmin in (1, 4):
        fit = powerlaw.Fit(deg, xmin=xmin, discrete=True, verbose=False)
        n = int((deg >= xmin).sum())
        for alt in ("lognormal", "exponential", "truncated_power_law", "stretched_exponential"):
            R, p = fit.distribution_compare("power_law", alt, normalized_ratio=True)
            rows.append(dict(fit_range=f"xmin={xmin} (n={n})",
                             alpha=round(fit.power_law.alpha, 2), KS_D=round(fit.power_law.D, 3),
                             alternative=alt, R=round(R, 2), p=round(p, 3)))
    with open(os.path.join(OUT, "table3_powerlaw.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    return rows


def main():
    G = load_graph()
    print(f"Loaded network: {G.number_of_nodes()} proteins, {G.number_of_edges()} directed edges")
    node_table(G)
    props, comm, deg = global_properties(G)
    print("\n=== Global topology (Table 2) ===")
    for k, v in props:
        print(f"  {k:34s} {v}")
    print("\n=== Communities (greedy modularity) ===")
    for i, c in enumerate(sorted(comm, key=len, reverse=True), 1):
        print(f"  Community {i} (n={len(c)}): {', '.join(sorted(c))}")
    pl = powerlaw_table(deg)
    if pl:
        print("\n=== Power-law fits (Table 3) ===")
        for r in pl:
            print(f"  {r['fit_range']:16s} alpha={r['alpha']:.2f} D={r['KS_D']:.3f} "
                  f"vs {r['alternative']:22s} R={r['R']:+.2f} p={r['p']:.3f}")
    print(f"\nOutput tables written to: {OUT}/")


if __name__ == "__main__":
    main()
