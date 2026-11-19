# Bdellovibrio bacteriovorus HD100 – Salmonella Typhimurium ATCC 14028 predicted PPI network

Data and fully reproducible analysis code for the systems-level interactome in:

> *In Vitro and In Vivo Predation Activity of* Bdellovibrio bacteriovorus *Against*
> Salmonella *Typhimurium and Systems-Level Interactome Analysis*

The network is a **directed graph of 25 proteins (15 *B. bacteriovorus* HD100
predator + 10 *S.* Typhimurium ATCC 14028 prey) and 38 evidence-annotated
interactions**. Every node carries its UniProt accession and every edge carries a
mechanism, an edge class (cross-organism / intra-predator / intra-prey) and an
**evidence tier**:

| Tier | Meaning | Count |
|------|---------|-------|
| `literature` | mechanism supported by a specific peer-reviewed study | 5 (13.2%) |
| `inference`  | predicted from protein function alone — a hypothesis | 33 (86.8%) |

> **No interaction is fabricated.** The five literature-tier edges are the only
> ones presented as experimentally grounded; all others are explicitly labelled
> functional-inference hypotheses to be tested experimentally.

## Repository layout

```
.
├── data/
│   ├── ppi_nodes.csv          # 25 proteins: organism, role, UniProt ID, function, reference
│   ├── ppi_edges.csv          # 38 edges: class, mechanism, evidence_tier, citation
│   └── network_layout.json    # cached node coordinates (collision-resolved) for Fig 6/9
├── src/
│   ├── analysis.py            # topology, centrality, communities, power-law, shortest paths
│   └── plotting.py            # shared styling + loaders (Times New Roman with serif fallback)
├── figures/
│   ├── fig06_network.py       # full predator–prey network (+ shared layout/draw)
│   ├── fig07_topology.py      # degree dist, in/out-degree, betweenness, eigenvector
│   ├── fig08_bipartite.py     # cross-organism subnetwork by evidence tier
│   ├── fig09_community.py     # 5-community structure
│   ├── fig10_edgeclass_evidence.py  # edge-class + evidence-tier composition
│   ├── fig11_shortestpath.py  # predator→prey shortest-path heatmap
│   └── fig12_powerlaw_ccdf.py # degree-distribution CCDF vs alternatives
├── outputs/                   # generated: metrics.json, powerlaw.json, table1–3 CSVs, figures
├── requirements.txt
└── README.md
```

## Install

```bash
python -m venv venv && source venv/bin/activate
pip install -r requirements.txt
```
Python ≥ 3.9. Dependencies: `networkx`, `numpy`, `matplotlib`, `powerlaw`, `pillow`.

## Reproduce

```bash
python src/analysis.py            # writes outputs/metrics.json, powerlaw.json, table1–3.csv
for f in figures/fig*.py; do python "$f"; done   # writes figures/*.jpg (600 dpi)
```

### Expected global metrics (printed by `src/analysis.py`)
25 nodes, 38 directed edges, density 0.127, mean degree 3.04 (SD 1.67),
max degree 6 (OmpA, OmpF), clustering 0.225, giant component = 21 nodes
(3 components), average shortest path 2.70, modularity **Q = 0.483**, 5 communities.
Edge classes: cross-organism 12 (31.6%), intra-predator 14 (36.8%), intra-prey 12 (31.6%).
Power law (maximum-likelihood, Clauset framework): α = 1.71 (x_min = 1, K–S D = 0.25);
**the power law is statistically rejected against all alternative heavy-tailed
distributions** (log-normal, exponential, truncated power law, stretched exponential).

## Edge evidence sources

The five literature-tier interactions cite:

- **Rendulic et al. 2004** *Science* 303:689 — HD100 genome / protein annotation.
- **Lerner et al. 2012** *PLoS Pathog* 8(2):e1002524 — Bd3459 DD-carboxypeptidase remodelling of prey peptidoglycan.
- **Lambert et al. 2016** *Sci Rep* 6:26010 — Bd3279 deacetylation of prey peptidoglycan.
- **Capeness et al. 2013** *PLoS ONE* 8(11):e79759 — Bd0108/Bd0109 hit-locus / type IVa pilus.
- **Bukowska-Faniband et al. 2020** *J Bacteriol* 202(18):e00150-20 — Bd0934/Bd3507 secreted nucleases.
- **Mun et al. 2022** *Microbiol Spectr* 10(6):e03094-22 — prey OmpF as a *Bdellovibrio* recognition receptor.

The full per-edge mechanism and citation are in `data/ppi_edges.csv`.

## Fonts

Figures use **Times New Roman** when available and fall back to the default serif
face otherwise. To match the manuscript, place `times*.ttf` in `figures/fonts/`.

## Data availability statement (for the manuscript)

> The curated node and edge datasets (`ppi_nodes.csv`, `ppi_edges.csv`) and the
> analysis scripts underlying the PPI network analysis are available at
> https://github.com/<user>/bdellovibrio-salmonella-ppi
> (archived at https://doi.org/10.5281/zenodo.XXXXXXXX).

## License

MIT (see `LICENSE`).
