#!/usr/bin/env python3
"""
network_analysis.py

Reproduces the predator-prey PPI network construction and analysis described in:
"Systems-Level Interactome Analysis of Bdellovibrio bacteriovorus Reveals Surface
Recognition Targets and Lytic Effectors Against Salmonella Typhimurium"

Inputs
------
ppi_nodes.csv : node attribute table (protein_id, organism, functional_category,
                subcellular_localization, molecular_weight_kDa, uniprot_accession)
ppi_edges.csv : curated interaction table (source, target, interaction_type,
                confidence_score, evidence_type, directionality, literature_reference)

Outputs
-------
outputs/table1_centrality_metrics.csv   Degree/betweenness/eigenvector centrality
                                         and community assignment per protein (Table 1)
outputs/table2_global_properties.csv    Global network topology summary (Table 2)
outputs/table3_community_detection.csv  Community membership summary (Table 3)
outputs/fig7_full_network.png           Full predator-prey network (Fig. 7)
outputs/fig8_degree_centrality.png      Degree distribution & centrality panels (Fig. 8)
outputs/fig9_bipartite_subnetwork.png   Cross-organism bipartite subnetwork (Fig. 9)
outputs/fig10_community_detection.png   Community detection visualization (Fig. 10)
outputs/fig11_confidence_evidence.png   Confidence score / evidence type panels (Fig. 11)
outputs/fig12_shortest_paths.png        Shortest-path heatmap, predator to prey (Fig. 12)

Requirements
------------
pip install networkx==3.6.1 matplotlib pandas numpy scipy python-louvain

Usage
-----
python network_analysis.py --nodes ppi_nodes.csv --edges ppi_edges.csv --outdir outputs
"""

import argparse
import os

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd


# --------------------------------------------------------------------------- #
# 1. Data loading and graph construction
# --------------------------------------------------------------------------- #

def load_network(nodes_csv: str, edges_csv: str) -> nx.DiGraph:
    """Build a weighted directed graph from the curated node/edge CSVs.

    Symmetric ('undirected') interactions are bidirectionally encoded as two
    directed edges, matching the manuscript's stated methodology (50 unique
    interactions -> 66 directed edges after bidirectional expansion).
    """
    nodes = pd.read_csv(nodes_csv)
    edges = pd.read_csv(edges_csv)

    G = nx.DiGraph()
    for _, row in nodes.iterrows():
        G.add_node(row["protein_id"], **row.drop("protein_id").to_dict())

    for _, row in edges.iterrows():
        attrs = row.drop(["source", "target"]).to_dict()
        G.add_edge(row["source"], row["target"], weight=row["confidence_score"], **attrs)
        if str(row["directionality"]).lower() == "undirected":
            G.add_edge(row["target"], row["source"], weight=row["confidence_score"], **attrs)

    return G


# --------------------------------------------------------------------------- #
# 2. Centrality analysis (Table 1)
# --------------------------------------------------------------------------- #

def compute_centrality_table(G: nx.DiGraph, communities: dict) -> pd.DataFrame:
    deg = dict(G.degree())
    indeg = dict(G.in_degree())
    outdeg = dict(G.out_degree())
    bc = nx.betweenness_centrality(G, weight="weight", normalized=True)
    try:
        ec = nx.eigenvector_centrality(G, weight="weight", max_iter=1000, tol=1e-6)
    except nx.PowerIterationFailedConvergence:
        ec = {n: 0.0 for n in G.nodes()}

    rows = []
    for n, data in G.nodes(data=True):
        rows.append({
            "Protein": n,
            "Organism": data.get("organism", ""),
            "Category": data.get("functional_category", ""),
            "Deg": deg[n],
            "In": indeg[n],
            "Out": outdeg[n],
            "BC": round(bc[n], 3),
            "EC": round(ec[n], 3),
            "Comm": communities.get(n, -1),
        })
    df = pd.DataFrame(rows).sort_values("BC", ascending=False).reset_index(drop=True)
    return df


# --------------------------------------------------------------------------- #
# 3. Global network properties (Table 2)
# --------------------------------------------------------------------------- #

def compute_global_properties(G: nx.DiGraph, communities: dict, modularity_q: float,
                               edges_df: pd.DataFrame) -> pd.DataFrame:
    UG = G.to_undirected()
    density = nx.density(G)
    avg_degree = np.mean([d for _, d in G.degree()])
    clustering = nx.average_clustering(UG)
    try:
        avg_path_len = nx.average_shortest_path_length(UG)
    except nx.NetworkXError:
        largest_cc = max(nx.connected_components(UG), key=len)
        avg_path_len = nx.average_shortest_path_length(UG.subgraph(largest_cc))

    deg_seq = sorted([d for _, d in G.degree() if d > 0])
    gamma, r2 = fit_power_law(deg_seq)

    orgs = nx.get_node_attributes(G, "organism")
    cross = sum(1 for u, v in G.edges() if orgs.get(u) != orgs.get(v))
    predator_org = "Bdellovibrio bacteriovorus"
    intra_pred = sum(1 for u, v in G.edges() if orgs.get(u) == predator_org and orgs.get(v) == predator_org)
    intra_prey = sum(1 for u, v in G.edges() if orgs.get(u) != predator_org and orgs.get(v) != predator_org
                      and orgs.get(u) == orgs.get(v))

    mean_conf = edges_df["confidence_score"].mean()
    pct_exp = (edges_df["evidence_type"] == "Experimental").mean() * 100

    rows = [
        ("Total Nodes", G.number_of_nodes()),
        ("Total Directed Edges", G.number_of_edges()),
        ("Cross-organism edges (Bb->St)", cross),
        ("Intra-predator edges", intra_pred),
        ("Intra-prey edges", intra_prey),
        ("Network Density", round(density, 3)),
        ("Average Degree", round(avg_degree, 2)),
        ("Average Clustering Coefficient", round(clustering, 3)),
        ("Average Shortest Path Length", round(avg_path_len, 2)),
        ("Power-law exponent (gamma)", round(gamma, 2) if gamma is not None else "NA"),
        ("Power-law R2", round(r2, 2) if r2 is not None else "NA"),
        ("Modularity Q", round(modularity_q, 3)),
        ("Number of Communities", len(set(communities.values()))),
        ("Mean Confidence Score", round(mean_conf, 3)),
        ("% Experimental Evidence", round(pct_exp, 1)),
    ]
    return pd.DataFrame(rows, columns=["Network Property", "Value"])


def fit_power_law(degree_sequence):
    """Simple log-log linear fit of the degree distribution tail."""
    if len(degree_sequence) < 3:
        return None, None
    values, counts = np.unique(degree_sequence, return_counts=True)
    mask = values > 0
    log_v = np.log(values[mask])
    log_c = np.log(counts[mask])
    if len(log_v) < 2:
        return None, None
    slope, intercept = np.polyfit(log_v, log_c, 1)
    pred = slope * log_v + intercept
    ss_res = np.sum((log_c - pred) ** 2)
    ss_tot = np.sum((log_c - np.mean(log_c)) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return -slope, r2


# --------------------------------------------------------------------------- #
# 4. Community detection (Table 3, Fig. 10)
# --------------------------------------------------------------------------- #

def detect_communities(G: nx.DiGraph):
    UG = G.to_undirected()
    comm_generator = nx.community.greedy_modularity_communities(UG, weight="weight")
    communities = {}
    for i, comm in enumerate(comm_generator, start=1):
        for n in comm:
            communities[n] = i
    q = nx.community.modularity(UG, comm_generator, weight="weight")
    return communities, q, list(comm_generator)


def community_summary_table(comm_list, G: nx.DiGraph) -> pd.DataFrame:
    rows = []
    for i, comm in enumerate(comm_list, start=1):
        members = sorted(comm)
        types = pd.Series([G.nodes[n].get("functional_category", "") for n in members])
        rows.append({
            "Community": i,
            "Size": len(members),
            "Member Proteins": ", ".join(members),
            "Dominant Functional Categories": ", ".join(types.value_counts().index[:3]),
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# 5. Cross-organism bipartite subgraph (Fig. 9) and shortest paths (Fig. 12)
# --------------------------------------------------------------------------- #

def cross_organism_subgraph(G: nx.DiGraph, predator_org="Bdellovibrio bacteriovorus"):
    orgs = nx.get_node_attributes(G, "organism")
    cross_edges = [(u, v, d) for u, v, d in G.edges(data=True)
                   if orgs.get(u) == predator_org and orgs.get(v) != predator_org]
    H = nx.DiGraph()
    H.add_edges_from(cross_edges)
    return H


def shortest_path_matrix(G: nx.DiGraph, predator_org="Bdellovibrio bacteriovorus"):
    orgs = nx.get_node_attributes(G, "organism")
    predators = sorted([n for n, o in orgs.items() if o == predator_org])
    prey = sorted([n for n, o in orgs.items() if o != predator_org])
    lengths = dict(nx.all_pairs_shortest_path_length(G))
    matrix = pd.DataFrame(index=predators, columns=prey, dtype=float)
    for p in predators:
        for q in prey:
            matrix.loc[p, q] = lengths.get(p, {}).get(q, np.inf)
    return matrix


# --------------------------------------------------------------------------- #
# 6. Plotting
# --------------------------------------------------------------------------- #

def plot_full_network(G, communities, outpath):
    pos = nx.spring_layout(G, k=2.8, seed=42)
    orgs = nx.get_node_attributes(G, "organism")
    colors = ["#2e7d32" if o == "Bdellovibrio bacteriovorus" else "#c62828" for o in orgs.values()]
    sizes = [80 + 40 * G.degree(n) for n in G.nodes()]
    plt.figure(figsize=(9, 7))
    nx.draw_networkx_edges(G, pos, alpha=0.25, arrows=True, width=0.6)
    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=sizes, alpha=0.9)
    nx.draw_networkx_labels(G, pos, font_size=6)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def plot_degree_and_centrality(G, table1, outpath):
    fig, axes = plt.subplots(2, 2, figsize=(10, 8))

    degrees = [d for _, d in G.degree()]
    axes[0, 0].hist(degrees, bins=range(1, max(degrees) + 2), color="#455a64", edgecolor="white")
    axes[0, 0].set_title("A) Degree distribution")
    axes[0, 0].set_xlabel("Degree")
    axes[0, 0].set_ylabel("Count")

    orgs = nx.get_node_attributes(G, "organism")
    indeg = dict(G.in_degree())
    outdeg = dict(G.out_degree())
    for n in G.nodes():
        c = "#2e7d32" if orgs.get(n) == "Bdellovibrio bacteriovorus" else "#c62828"
        axes[0, 1].scatter(indeg[n], outdeg[n], color=c, s=40, alpha=0.8)
    axes[0, 1].set_title("B) In-degree vs out-degree")
    axes[0, 1].set_xlabel("In-degree")
    axes[0, 1].set_ylabel("Out-degree")

    top_bc = table1.sort_values("BC", ascending=False).head(15)
    axes[1, 0].barh(top_bc["Protein"], top_bc["BC"], color="#5c6bc0")
    axes[1, 0].invert_yaxis()
    axes[1, 0].set_title("C) Top-15 betweenness centrality")

    top_ec = table1.sort_values("EC", ascending=False).head(15)
    axes[1, 1].barh(top_ec["Protein"], top_ec["EC"], color="#26a69a")
    axes[1, 1].invert_yaxis()
    axes[1, 1].set_title("D) Top-15 eigenvector centrality")

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def plot_bipartite_subnetwork(H, outpath):
    predators = sorted({u for u, _ in H.edges()})
    prey = sorted({v for _, v in H.edges()})
    pos = {}
    for i, n in enumerate(predators):
        pos[n] = (0, -i)
    for i, n in enumerate(prey):
        pos[n] = (1, -i)
    plt.figure(figsize=(6, 9))
    nx.draw_networkx_nodes(H, pos, nodelist=predators, node_color="#2e7d32", node_size=300)
    nx.draw_networkx_nodes(H, pos, nodelist=prey, node_color="#c62828", node_size=300)
    nx.draw_networkx_edges(H, pos, arrows=True, alpha=0.5)
    nx.draw_networkx_labels(H, pos, font_size=7)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def plot_community_detection(G, communities, outpath):
    pos = nx.spring_layout(G, k=2.8, seed=42)
    cmap = plt.get_cmap("tab10")
    colors = [cmap((communities.get(n, 0) - 1) % 10) for n in G.nodes()]
    sizes = [80 + 40 * G.degree(n) for n in G.nodes()]
    plt.figure(figsize=(8, 7))
    nx.draw_networkx_edges(G, pos, alpha=0.2, width=0.6)
    nx.draw_networkx_nodes(G, pos, node_color=colors, node_size=sizes)
    nx.draw_networkx_labels(G, pos, font_size=6)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def plot_confidence_evidence(edges_df, outpath):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))

    axes[0].hist(edges_df["confidence_score"], bins=10, color="#5c6bc0", edgecolor="white")
    axes[0].axvline(edges_df["confidence_score"].mean(), color="red", linestyle="--", label="mean")
    axes[0].set_title("A) Confidence score distribution")
    axes[0].legend()

    evidence_counts = edges_df["evidence_type"].value_counts()
    axes[1].pie(evidence_counts, labels=evidence_counts.index, autopct="%1.0f%%",
                colors=["#26a69a", "#ef5350"])
    axes[1].set_title("B) Evidence type")

    type_counts = edges_df["interaction_type"].value_counts().head(8)
    axes[2].barh(type_counts.index, type_counts.values, color="#8d6e63")
    axes[2].invert_yaxis()
    axes[2].set_title("C) Interaction type frequency")

    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


def plot_shortest_path_heatmap(matrix, outpath):
    plt.figure(figsize=(8, 6))
    display_matrix = matrix.replace(np.inf, np.nan)
    plt.imshow(display_matrix.astype(float), cmap="RdYlGn_r", aspect="auto")
    plt.colorbar(label="Shortest path length (hops)")
    plt.xticks(range(len(matrix.columns)), matrix.columns, rotation=90, fontsize=7)
    plt.yticks(range(len(matrix.index)), matrix.index, fontsize=7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()


# --------------------------------------------------------------------------- #
# 7. Main pipeline
# --------------------------------------------------------------------------- #

def main():
    parser = argparse.ArgumentParser(description="B. bacteriovorus - S. Typhimurium PPI network analysis")
    parser.add_argument("--nodes", default="ppi_nodes.csv")
    parser.add_argument("--edges", default="ppi_edges.csv")
    parser.add_argument("--outdir", default="outputs")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    G = load_network(args.nodes, args.edges)
    edges_df = pd.read_csv(args.edges)

    communities, modularity_q, comm_list = detect_communities(G)

    table1 = compute_centrality_table(G, communities)
    table1.to_csv(os.path.join(args.outdir, "table1_centrality_metrics.csv"), index=False)

    table2 = compute_global_properties(G, communities, modularity_q, edges_df)
    table2.to_csv(os.path.join(args.outdir, "table2_global_properties.csv"), index=False)

    table3 = community_summary_table(comm_list, G)
    table3.to_csv(os.path.join(args.outdir, "table3_community_detection.csv"), index=False)

    plot_full_network(G, communities, os.path.join(args.outdir, "fig7_full_network.png"))
    plot_degree_and_centrality(G, table1, os.path.join(args.outdir, "fig8_degree_centrality.png"))

    H = cross_organism_subgraph(G)
    plot_bipartite_subnetwork(H, os.path.join(args.outdir, "fig9_bipartite_subnetwork.png"))

    plot_community_detection(G, communities, os.path.join(args.outdir, "fig10_community_detection.png"))
    plot_confidence_evidence(edges_df, os.path.join(args.outdir, "fig11_confidence_evidence.png"))

    sp_matrix = shortest_path_matrix(G)
    sp_matrix.to_csv(os.path.join(args.outdir, "shortest_path_matrix.csv"))
    plot_shortest_path_heatmap(sp_matrix, os.path.join(args.outdir, "fig12_shortest_paths.png"))

    print(f"Nodes: {G.number_of_nodes()}, Directed edges: {G.number_of_edges()}")
    print(f"Modularity Q = {modularity_q:.3f}, Communities = {len(comm_list)}")
    print(f"Outputs written to: {os.path.abspath(args.outdir)}")


if __name__ == "__main__":
    main()
