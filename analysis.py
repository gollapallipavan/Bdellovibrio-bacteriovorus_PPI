"""
Systems-level PPI network analysis for the Bdellovibrio bacteriovorus HD100 -
Salmonella Typhimurium ATCC 14028 predicted interactome.

Reads   : data/ppi_nodes.csv, data/ppi_edges.csv
Writes  : outputs/metrics.json, outputs/powerlaw.json,
          outputs/table1_centrality.csv, outputs/table2_global.csv,
          outputs/table3_powerlaw.csv

Centralities are computed on the UNDIRECTED graph (standard for interactome
analysis; avoids directed-sink degeneracy); in/out-degree are retained from the
directed graph. No interaction is fabricated: every edge in ppi_edges.csv carries
a mechanism, an edge class and an evidence tier (literature vs functional inference).
"""
import csv, json, os, warnings
warnings.filterwarnings("ignore")
import numpy as np
import networkx as nx

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data")
OUT  = os.path.join(HERE, "outputs")
os.makedirs(OUT, exist_ok=True)


def load():
    nodes = list(csv.DictReader(open(os.path.join(DATA, "ppi_nodes.csv"))))
    edges = list(csv.DictReader(open(os.path.join(DATA, "ppi_edges.csv"))))
    G = nx.DiGraph()
    G.add_nodes_from(n["node_id"] for n in nodes)
    for e in edges:
        G.add_edge(e["source"], e["target"], **e)
    # validate: every endpoint is a declared node
    ids = {n["node_id"] for n in nodes}
    bad = [(e["source"], e["target"]) for e in edges
           if e["source"] not in ids or e["target"] not in ids]
    assert not bad, f"edge endpoints not in node list: {bad}"
    return nodes, edges, G


def centrality_table(nodes, G):
    UG = G.to_undirected()
    deg = dict(UG.degree()); ind = dict(G.in_degree()); outd = dict(G.out_degree())
    bc = nx.betweenness_centrality(UG)
    ev = nx.eigenvector_centrality(UG, max_iter=10000)
    comms = list(nx.algorithms.community.greedy_modularity_communities(UG))
    comm = {}
    for i, c in enumerate(sorted(comms, key=lambda s: -len(s)), 1):
        for n in c:
            comm[n] = i
    return deg, ind, outd, bc, ev, comm, comms


def power_law(G):
    import powerlaw
    UG = G.to_undirected()
    degs = np.array(sorted([d for _, d in UG.degree()]))
    res = {}
    for xmin in (1, 2):
        fit = powerlaw.Fit(degs, xmin=xmin, discrete=True, verbose=False)
        row = {"xmin": xmin, "n": int((degs >= xmin).sum()),
               "alpha": round(fit.power_law.alpha, 3), "D": round(fit.power_law.D, 3),
               "cmp": {}}
        for alt in ("lognormal", "exponential", "truncated_power_law", "stretched_exponential"):
            R, p = fit.distribution_compare("power_law", alt, normalized_ratio=True)
            row["cmp"][alt] = (round(R, 2), round(p, 4))
        res[f"xmin{xmin}"] = row
    # naive OLS log-log fit (reported for completeness only)
    from collections import Counter
    c = Counter(degs); ks = sorted(c)
    x = np.log10(ks); y = np.log10([c[k] for k in ks])
    A = np.vstack([x, np.ones_like(x)]).T
    m, b = np.linalg.lstsq(A, y, rcond=None)[0]
    yp = A @ [m, b]
    r2 = 1 - np.sum((y - yp) ** 2) / np.sum((y - y.mean()) ** 2)
    res["ols"] = {"slope": round(m, 3), "r2": round(r2, 3)}
    return res


def shortest_paths(nodes, G):
    BB = [n["node_id"] for n in nodes if n["role"] == "predator"]
    ST = [n["node_id"] for n in nodes if n["role"] == "prey"]
    hops = []
    pairs = [(u, v) for u in BB for v in ST]
    for u, v in pairs:
        if nx.has_path(G, u, v):
            hops.append(nx.shortest_path_length(G, u, v))
    hops = np.array(hops)
    return {"pairs": len(pairs), "reachable": int(len(hops)),
            "reach_pct": round(100 * len(hops) / len(pairs), 1),
            "within3": int((hops <= 3).sum()),
            "within3_pct_reach": round(100 * (hops <= 3).sum() / len(hops), 1),
            "mean_hops": round(float(hops.mean()), 2),
            "sd_hops": round(float(hops.std(ddof=1)), 2)}


def main():
    nodes, edges, G = load()
    UG = G.to_undirected()
    deg, ind, outd, bc, ev, comm, comms = centrality_table(nodes, G)
    BB = set(n["node_id"] for n in nodes if n["role"] == "predator")
    ST = set(n["node_id"] for n in nodes if n["role"] == "prey")
    degs = [deg[n] for n in G]
    comps = list(nx.connected_components(UG)); giant = max(comps, key=len)
    cross = sum(1 for e in edges if e["edge_class"] == "cross")
    ipred = sum(1 for e in edges if e["edge_class"] == "intra_predator")
    iprey = sum(1 for e in edges if e["edge_class"] == "intra_prey")
    lit = sum(1 for e in edges if e["evidence_tier"] == "literature")
    E = G.number_of_edges()
    summary = {
        "N": G.number_of_nodes(), "E": E, "density": round(nx.density(UG), 3),
        "mean_deg": round(float(np.mean(degs)), 2), "sd_deg": round(float(np.std(degs, ddof=1)), 2),
        "max_deg": max(degs), "hub": sorted([n for n in G if deg[n] == max(degs)]),
        "clustering": round(nx.average_clustering(UG), 3),
        "avg_path": round(nx.average_shortest_path_length(UG.subgraph(giant)), 2),
        "Q": round(nx.algorithms.community.modularity(UG, comms), 3),
        "ncomm": len(comms), "components": len(comps), "giant": len(giant),
        "cross": cross, "ipred": ipred, "iprey": iprey, "lit": lit, "inf": E - lit,
        "cross_pct": round(100 * cross / E, 1), "ipred_pct": round(100 * ipred / E, 1),
        "iprey_pct": round(100 * iprey / E, 1),
        "lit_pct": round(100 * lit / E, 1), "inf_pct": round(100 * (E - lit) / E, 1),
    }
    pl = power_law(G); sp = shortest_paths(nodes, G)

    json.dump({"summary": summary, "deg": deg, "ind": ind, "outd": outd,
               "bc": {k: round(v, 4) for k, v in bc.items()},
               "ev": {k: round(v, 4) for k, v in ev.items()},
               "comm": comm, "communities": [sorted(c) for c in sorted(comms, key=lambda s: -len(s))],
               "BB": sorted(BB), "ST": sorted(ST)},
              open(os.path.join(OUT, "metrics.json"), "w"), indent=1)
    json.dump({**pl, "shortest_paths": sp}, open(os.path.join(OUT, "powerlaw.json"), "w"), indent=1)

    # Table 1 - centrality
    with open(os.path.join(OUT, "table1_centrality.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["protein", "organism", "degree", "in_degree",
                                       "out_degree", "betweenness", "eigenvector", "community"])
        for n in sorted(G, key=lambda n: -bc[n]):
            org = "B. bacteriovorus" if n in BB else "S. Typhimurium"
            w.writerow([n, org, deg[n], ind[n], outd[n], f"{bc[n]:.3f}", f"{ev[n]:.3f}", comm[n]])
    # Table 2 - global
    with open(os.path.join(OUT, "table2_global.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["property", "value"])
        for k, v in summary.items():
            w.writerow([k, v if not isinstance(v, list) else ", ".join(v)])
    # Table 3 - power-law
    with open(os.path.join(OUT, "table3_powerlaw.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(["fit_range", "statistic", "value"])
        for tag, lab in (("xmin1", "full (xmin=1)"), ("xmin2", "restricted (xmin=2)")):
            r = pl[tag]
            w.writerow([lab, "n", r["n"]]); w.writerow([lab, "alpha", r["alpha"]])
            w.writerow([lab, "KS_D", r["D"]])
            for alt, (R, p) in r["cmp"].items():
                w.writerow([lab, f"LLR_vs_{alt} (R;p)", f"{R:+.2f}; {p}"])
        w.writerow(["full", "OLS_slope; R2", f"{pl['ols']['slope']}; {pl['ols']['r2']}"])

    print("== GLOBAL SUMMARY ==")
    print(json.dumps(summary, indent=1))
    print("== POWER LAW ==")
    for tag in ("xmin1", "xmin2"):
        print(tag, pl[tag]["alpha"], "D", pl[tag]["D"], pl[tag]["cmp"])
    print("== PREDATOR->PREY SHORTEST PATHS ==", sp)
    print("\nwrote outputs/ (metrics.json, powerlaw.json, table1-3 CSVs)")


if __name__ == "__main__":
    main()
