# Bdellovibrio bacteriovorus – Salmonella Typhimurium PPI network

Data and reproducible analysis code for the protein–protein interaction (PPI)
network in:

> *In Vitro and In Vivo Predation Activity of* Bdellovibrio bacteriovorus *Against*
> Salmonella *Typhimurium and Systems-Level Interactome Analysis*

The network is a directed, weighted **protein–protein interaction network of 25
proteins and 52 directed edges** (15 *B. bacteriovorus* HD100-derived predator
proteins and 10 *S.* Typhimurium prey proteins). The non-protein species LPS (a
glycolipid) and DNA/RNA (nucleic acids), considered in an earlier draft, are
**not** part of this PPI network.

> **Important:** the curated interactions are HD100-derived, computationally
> generated hypotheses. None are experimentally validated in the specific
> LBD5–*S.* Typhimurium system. Mechanism labels assigned to protein targets that
> are inconsistent with the target's identity (e.g. carbohydrate/lipid/nuclease
> activities on protein nodes) are retained only as *predicted (unverified)*
> associations.

## Repository layout

```
.
├── ppi_analysis.py            # builds the network; reproduces Tables 1–3 + communities + power-law
├── data/
│   ├── ppi_nodes.csv          # 25 proteins: organism, role, category, function, centralities, community
│   ├── ppi_edges.csv          # 52 directed edges: edge_class + reconstructed mechanism
│   └── network_layout.json    # node coordinates used by fig_network.py
├── figures/
│   ├── plot_utils.py          # shared styling + helpers (Times New Roman with serif fallback)
│   ├── fig_network.py         # full predator–prey network
│   ├── fig_degree_centrality.py   # degree distribution + in/out scatter + centrality bars
│   ├── fig_bipartite.py       # cross-organism map coloured by interaction mechanism
│   ├── fig_shortest_paths.py  # predator→prey shortest-path heatmap
│   └── fig_powerlaw_ccdf.py   # power-law vs alternative-distribution CCDF fits
├── outputs/                   # generated tables (CSV) and figures (TIFF + PDF)
├── requirements.txt
└── README.md
```

## Requirements

```
pip install -r requirements.txt
```

* Python ≥ 3.9, `networkx` ≥ 3, `numpy`, `matplotlib`, `pillow`
* `powerlaw` (for Table 3 / the CCDF figure; the analysis otherwise skips it)

## Reproduce the analysis

```bash
python ppi_analysis.py          # writes outputs/table1_centrality.csv, table2_global.csv, table3_powerlaw.csv
cd figures
python fig_network.py
python fig_degree_centrality.py
python fig_bipartite.py
python fig_shortest_paths.py
python fig_powerlaw_ccdf.py      # requires `powerlaw`
```

Expected global metrics (printed by `ppi_analysis.py`): 25 nodes, 52 edges,
density 0.087, mean degree 4.16, max degree 11 (OmpA), clustering 0.192,
average path length 2.91, 5 communities, modularity Q = 0.424. Power-law
(maximum-likelihood, Clauset framework): α = 1.58 (x_min = 1, K–S D = 0.33),
α = 3.00 (x_min = 4, 14 nodes); the power law is rejected against all alternative
heavy-tailed distributions at x_min = 1.

## Fonts

The figures use **Times New Roman** when available and fall back to the default
serif face otherwise (the manuscript figures use Times New Roman). Times New
Roman is a licensed font and is **not** bundled here; to match the manuscript,
place `times.ttf`, `timesbd.ttf`, `timesi.ttf`, `timesbi.ttf` in `figures/fonts/`.

## Notes / limitations

* **Edge mechanisms** (`mechanism_reconstructed` in `ppi_edges.csv`) are provided
  for the cross-organism edges only and are reconstructed from the manuscript
  figures; biochemically inconsistent labels are set to
  `Predicted (unverified mechanism)`.
* **Confidence scores, evidence provenance, and interaction-type frequencies**
  (the "interaction-quality" figure in the manuscript) depend on the authors'
  full edge-level curation file and are **not** reproduced here; add the
  `confidence_score`, `evidence_type` and `literature_reference` columns to
  `ppi_edges.csv` to regenerate them.

## Data availability statement (for the manuscript)

> The curated node and edge datasets (`ppi_nodes.csv`, `ppi_edges.csv`) and the
> analysis scripts underlying the PPI network analysis are available at
> https://doi.org/10.5281/zenodo.XXXXXXXX  *(replace with the DOI assigned on
> deposition)*.
