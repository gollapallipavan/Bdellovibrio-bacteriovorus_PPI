"""
build_network.py
================
PPI network CONSTRUCTION for the Bdellovibrio bacteriovorus HD100 - Salmonella
Typhimurium ATCC 14028 predicted interactome.

This single script defines the complete network - every node (with its UniProt
accession) and every interaction (with its mechanism, edge class, evidence tier
and citation) - and builds the directed graph. It is the authoritative, auditable
record of how the network was assembled.

HOW THE NETWORK WAS BUILT
-------------------------
* Nodes: the 25 proteins supplied for the study (15 B. bacteriovorus HD100 predator
  proteins and 10 S. Typhimurium ATCC 14028 prey proteins), each mapped to a UniProt
  accession. No node is invented.
* Edges: no experimentally determined cross-organism interactome exists for this
  predator-prey pair, and public interaction databases hold no inter-species edges,
  so interactions were NOT downloaded from a database. Each candidate interaction was
  assembled by hand from (a) documented B. bacteriovorus predation biology and
  (b) the UniProt/KEGG functional annotation of the two proteins.
* Evidence tiering (no fabrication): every edge is tagged
      "literature" - a specific peer-reviewed study supports the mechanism, or
      "inference"  - the edge is predicted from protein function alone (a hypothesis).
  Only the 5 literature-tier edges are presented as experimentally grounded; the
  other 33 are explicitly labelled functional-inference hypotheses to be tested.
* The graph is DIRECTED (predator -> prey for cross-organism edges; functional
  direction for intra-organism edges). No bidirectional expansion is applied.

OUTPUTS (written under ../data and ../outputs)
----------------------------------------------
    data/ppi_nodes.csv      data/ppi_edges.csv      outputs/network.graphml

Run:  python build_network.py
"""
import csv, os
import networkx as nx

# ---------------------------------------------------------------------------
# 1. NODES  (protein, organism, role, uniprot_id, function, reference)
# ---------------------------------------------------------------------------
NODES = [
    ("DnaK", "B. bacteriovorus HD100", "predator", "Q6MNF8", "Chaperone / stress response / ATP binding and hydrolysis", "Rendulic et al. 2004"),
    ("Bd1972", "B. bacteriovorus HD100", "predator", "Q6MLN5", "ABC-type dipeptide/oligopeptide/nickel transport system", "Rendulic et al. 2004"),
    ("Bd0867", "B. bacteriovorus HD100", "predator", "Q6MPI6", "Fimbrial assembly / secretory protein", "Rendulic et al. 2004"),
    ("Bd1509", "B. bacteriovorus HD100", "predator", "Q6MMW4", "ATP hydrolysis / pilus binding", "Rendulic et al. 2004"),
    ("Bd3459", "B. bacteriovorus HD100", "predator", "Q6MHT0", "D-alanyl-D-alanine carboxypeptidase / proteolysis", "Lerner et al. 2012"),
    ("Bd0934", "B. bacteriovorus HD100", "predator", "Q6MPC3", "Hydrolase / secreted nuclease", "Rendulic et al. 2004"),
    ("Bd0109", "B. bacteriovorus HD100", "predator", "Q6MRH5", "Cell wall-associated protein (WapA-like)", "Rendulic et al. 2004"),
    ("Bd0108", "B. bacteriovorus HD100", "predator", "Q6MRH6", "Hit locus orf4", "Rendulic et al. 2004"),
    ("Bd0537", "B. bacteriovorus HD100", "predator", "Q6MQD5", "Negative regulator of flagellin synthesis", "Rendulic et al. 2004"),
    ("Bd2352", "B. bacteriovorus HD100", "predator", "Q6MKN6", "Lipoprotein", "Rendulic et al. 2004"),
    ("Bd3279", "B. bacteriovorus HD100", "predator", "Q6MI90", "Polysaccharide deacetylase", "Lambert et al. 2016"),
    ("Bd3507", "B. bacteriovorus HD100", "predator", "Q6MHN3", "Secreted nuclease / hydrolase", "Rendulic et al. 2004"),
    ("FliF", "B. bacteriovorus HD100", "predator", "Q6MHY0", "Flagellar M-ring protein", "Rendulic et al. 2004"),
    ("FtsZ", "B. bacteriovorus HD100", "predator", "Q6MIG9", "Cell division protein FtsZ", "Rendulic et al. 2004"),
    ("Bd0875", "B. bacteriovorus HD100", "predator", "Q6MPH9", "VWFA domain-containing protein", "Rendulic et al. 2004"),
    ("OmpA", "S. Typhimurium 14028", "prey", "A0A0F6AZN3", "Outer membrane protein A / porin activity", "UniProt"),
    ("MurA", "S. Typhimurium 14028", "prey", "A0A0F6B793", "UDP-N-acetylglucosamine 1-carboxyvinyltransferase / peptidoglycan biosynthesis", "UniProt"),
    ("BcsA", "S. Typhimurium 14028", "prey", "A0A0F6B898", "Cellulose synthase catalytic subunit / biofilm", "UniProt"),
    ("SipA", "S. Typhimurium 14028", "prey", "A0A0F6B5V0", "Secreted effector protein (T3SS)", "UniProt"),
    ("OmpC", "S. Typhimurium 14028", "prey", "A0A0F6B3Z7", "Porin activity", "UniProt"),
    ("RfaL", "S. Typhimurium 14028", "prey", "A0A0F6B8L0", "O-antigen ligase", "UniProt"),
    ("OmpF", "S. Typhimurium 14028", "prey", "A0A0F6AZF0", "Porin activity", "UniProt"),
    ("YcfR", "S. Typhimurium 14028", "prey", "A0A0F6B051", "Outer membrane protein / biofilm", "UniProt"),
    ("SseB", "S. Typhimurium 14028", "prey", "A0A0F6B4V0", "Enhanced serine sensitivity protein SseB (T3SS translocon)", "UniProt"),
    ("FimC", "S. Typhimurium 14028", "prey", "A0A0F6AY32", "Periplasmic fimbrial chaperone", "UniProt"),
]

# ---------------------------------------------------------------------------
# 2. EDGES  (source, target, edge_class, mechanism, evidence_tier, citation)
#    edge_class in {cross, intra_predator, intra_prey}
#    evidence_tier in {literature, inference}
# ---------------------------------------------------------------------------
EDGES = [
    ("Bd3279", "MurA", "cross", "Deacetylation of prey peptidoglycan during predation", "literature", "Lambert et al. 2016 (Sci Rep 6:26010)"),
    ("Bd3459", "MurA", "cross", "DD-carboxypeptidase remodelling of prey peptidoglycan", "literature", "Lerner et al. 2012 (PLoS Pathog 8:e1002524)"),
    ("Bd0867", "OmpF", "cross", "Fimbrial/secretory attachment to prey porin recognition receptor", "inference", "Mun et al. 2022 (Microbiol Spectr 10:e03094-22)"),
    ("Bd1509", "OmpF", "cross", "Pilus-mediated attachment to prey porin receptor", "inference", "Mun et al. 2022 (Microbiol Spectr 10:e03094-22)"),
    ("Bd0875", "OmpA", "cross", "VWFA-domain adhesion to prey outer-membrane protein", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("Bd0875", "OmpF", "cross", "VWFA-domain adhesion to prey porin", "inference", "Mun et al. 2022 (Microbiol Spectr 10:e03094-22)"),
    ("Bd0867", "OmpC", "cross", "Fimbrial recognition of prey porin", "inference", "Mun et al. 2022 (Microbiol Spectr 10:e03094-22)"),
    ("Bd1509", "OmpA", "cross", "Pilus attachment to prey outer-membrane protein", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("Bd3279", "BcsA", "cross", "Polysaccharide deacetylase acting on prey exopolysaccharide/biofilm", "inference", "Lambert et al. 2016 (Sci Rep 6:26010)"),
    ("Bd0109", "MurA", "cross", "WapA-like cell-wall-associated engagement of prey peptidoglycan", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("Bd0875", "RfaL", "cross", "Adhesion engaging prey O-antigen/LPS surface", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("Bd0108", "OmpF", "cross", "Hit-locus/pilus-linked recognition of prey porin", "inference", "Capeness et al. 2013 (PLoS One 8:e79759); Mun et al. 2022"),
    ("Bd0108", "Bd0109", "intra_predator", "Hit-locus pair regulating type IVa pilus extrusion/retraction", "literature", "Capeness et al. 2013 (PLoS One 8:e79759)"),
    ("Bd1509", "Bd0867", "intra_predator", "Pilus-assembly ATPase with fimbrial/secretory protein", "inference", "Capeness et al. 2013 (PLoS One 8:e79759)"),
    ("Bd0108", "Bd1509", "intra_predator", "Hit-locus control of pilus extrusion/retraction", "inference", "Capeness et al. 2013 (PLoS One 8:e79759)"),
    ("Bd0109", "Bd0867", "intra_predator", "Pilus-associated cell-wall protein with fimbrial protein", "inference", "Capeness et al. 2013 (PLoS One 8:e79759)"),
    ("Bd0934", "Bd3507", "intra_predator", "Co-regulated sequentially released secreted nucleases", "literature", "Bukowska-Faniband et al. 2020 (J Bacteriol 202:e00150-20)"),
    ("Bd0537", "FliF", "intra_predator", "Flagellin-synthesis regulator with flagellar M-ring", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("FtsZ", "Bd3459", "intra_predator", "Cell-division peptidoglycan remodelling", "inference", "Lerner et al. 2012 (PLoS Pathog 8:e1002524)"),
    ("DnaK", "FtsZ", "intra_predator", "Chaperone assisting division-protein folding", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("DnaK", "FliF", "intra_predator", "Chaperone assisting flagellar-protein folding", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("DnaK", "Bd1509", "intra_predator", "Chaperone assisting pilus ATPase", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("DnaK", "Bd1972", "intra_predator", "Chaperone assisting ABC-transporter assembly", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("Bd0875", "Bd0867", "intra_predator", "VWFA adhesion protein with fimbrial assembly protein", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("Bd2352", "Bd0875", "intra_predator", "Surface lipoprotein with adhesion protein", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("Bd1972", "Bd3459", "intra_predator", "Peptide transport with peptidoglycan carboxypeptidase", "inference", "Rendulic et al. 2004 (Science 303:689)"),
    ("OmpA", "OmpC", "intra_prey", "Outer-membrane porin family co-occurrence", "inference", "Mun et al. 2022 (Microbiol Spectr 10:e03094-22)"),
    ("OmpC", "OmpF", "intra_prey", "EnvZ-OmpR co-regulated porins", "literature", "Mun et al. 2022 (Microbiol Spectr 10:e03094-22)"),
    ("OmpA", "OmpF", "intra_prey", "Outer-membrane porin family co-occurrence", "inference", "Mun et al. 2022 (Microbiol Spectr 10:e03094-22)"),
    ("BcsA", "FimC", "intra_prey", "Biofilm cellulose matrix with fimbrial chaperone", "inference", "UniProt functional annotation"),
    ("BcsA", "YcfR", "intra_prey", "Cellulose biosynthesis with biofilm-associated OMP", "inference", "UniProt functional annotation"),
    ("FimC", "YcfR", "intra_prey", "Fimbrial/biofilm module", "inference", "UniProt functional annotation"),
    ("RfaL", "OmpA", "intra_prey", "O-antigen ligase with outer-membrane protein biogenesis", "inference", "UniProt functional annotation"),
    ("RfaL", "OmpC", "intra_prey", "LPS/outer-membrane biogenesis", "inference", "UniProt functional annotation"),
    ("SipA", "SseB", "intra_prey", "T3SS effector with translocon component", "inference", "UniProt functional annotation"),
    ("MurA", "RfaL", "intra_prey", "Cell-envelope biosynthesis (peptidoglycan with LPS)", "inference", "UniProt functional annotation"),
    ("FimC", "OmpA", "intra_prey", "Fimbrial chaperone with outer-membrane protein", "inference", "UniProt functional annotation"),
    ("MurA", "BcsA", "intra_prey", "Envelope polysaccharide synthesis (peptidoglycan with cellulose)", "inference", "UniProt functional annotation"),
]

# ---------------------------------------------------------------------------
# 3. CONSTRUCT the directed graph
# ---------------------------------------------------------------------------
def build_graph():
    G = nx.DiGraph()
    for pid, org, role, uni, func, ref in NODES:
        G.add_node(pid, organism=org, role=role, uniprot_id=uni,
                   function=func, reference=ref)
    for s, t, cls, mech, tier, cite in EDGES:
        assert s in G and t in G, f"edge endpoint not a declared node: {s}->{t}"
        G.add_edge(s, t, edge_class=cls, mechanism=mech,
                   evidence_tier=tier, citation=cite)
    return G


def export(G, root=None):
    root = root or os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    data = os.path.join(root, "data"); out = os.path.join(root, "outputs")
    os.makedirs(data, exist_ok=True); os.makedirs(out, exist_ok=True)
    with open(os.path.join(data, "ppi_nodes.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n"); w.writerow(["node_id", "organism", "role",
                                       "uniprot_id", "function", "reference"])
        w.writerows(NODES)
    with open(os.path.join(data, "ppi_edges.csv"), "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n"); w.writerow(["source", "target", "edge_class",
                                       "mechanism", "evidence_tier", "citation"])
        w.writerows(EDGES)
    nx.write_graphml(G, os.path.join(out, "network.graphml"))


def summary(G):
    BB = [n for n, d in G.nodes(data=True) if d["role"] == "predator"]
    ST = [n for n, d in G.nodes(data=True) if d["role"] == "prey"]
    ec = lambda k: sum(1 for _, _, d in G.edges(data=True) if d["edge_class"] == k)
    cross, ipred, iprey = ec("cross"), ec("intra_predator"), ec("intra_prey")
    lit = sum(1 for _, _, d in G.edges(data=True) if d["evidence_tier"] == "literature")
    E = G.number_of_edges()
    print(f"Nodes: {G.number_of_nodes()}  (predator {len(BB)}, prey {len(ST)})")
    print(f"Edges: {E}  (cross {cross}, intra-predator {ipred}, intra-prey {iprey})")
    print(f"Evidence: literature {lit} ({100*lit/E:.1f}%), "
          f"inference {E-lit} ({100*(E-lit)/E:.1f}%)")
    print("Literature-documented edges:")
    for s, t, d in G.edges(data=True):
        if d["evidence_tier"] == "literature":
            print(f"   {s:8}-> {t:6}  {d['mechanism']}  [{d['citation']}]")


if __name__ == "__main__":
    G = build_graph()
    summary(G)
    export(G)
    print("\nwrote data/ppi_nodes.csv, data/ppi_edges.csv, outputs/network.graphml")
