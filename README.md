# B. bacteriovorus - S. Typhimurium PPI Network: Data and Code

This repository accompanies "Systems-Level Interactome Analysis of *Bdellovibrio
bacteriovorus* Reveals Surface Recognition Targets and Lytic Effectors Against
*Salmonella* Typhimurium." It contains the curated node/edge tables and the
NetworkX analysis pipeline used to construct the predator-prey protein-protein
interaction (PPI) network and generate the centrality, community-detection, and
shortest-path analyses reported in the manuscript.

## Files

| File | Description |
|---|---|
| `ppi_nodes.csv` | 27 protein nodes (15 *B. bacteriovorus*, 12 *S.* Typhimurium) with organism, functional category, and subcellular localization. |
| `ppi_edges.csv` | 50 curated interactions (source, target, interaction type, confidence score, evidence type, directionality, literature reference); expands to 66 directed edges after bidirectional encoding of symmetric interactions. |
| `network_analysis.py` | NetworkX pipeline: graph construction, centrality analysis (Table 1), global topology metrics (Table 2), greedy-modularity community detection (Table 3, Fig. 10), the cross-organism bipartite subnetwork (Fig. 9), confidence/evidence summaries (Fig. 11), and the predator-to-prey shortest-path heatmap (Fig. 12). |
| `requirements.txt` | Python dependencies. |

## Usage

```bash
pip install -r requirements.txt
python network_analysis.py --nodes ppi_nodes.csv --edges ppi_edges.csv --outdir outputs
```

Outputs (tables as CSV, figures as PNG) are written to `outputs/`.

## ⚠️ Important note on data provenance

The manuscript text describes the network construction methodology and reports
summary values (27 nodes, 66 directed edges, density = 0.094, modularity Q =
0.381, mean confidence = 0.874, etc.) and a top-10 highest-confidence interaction
table, but does not enumerate the complete original 50-interaction edge list
in the main text.

**`ppi_edges.csv` in this repository is a reconstruction**, built from:
- The 10 fully-specified interactions in the manuscript's top-confidence table (exact source, target, type, confidence score, and reference).
- The community memberships, degree profiles, and interaction descriptions given in the Results and Discussion sections (e.g., which predator enzymes target which prey proteins, and which relationships are symmetric such as chaperone-client or protein-complex pairs).

Running the pipeline on this reconstruction reproduces the manuscript's node/edge
counts exactly (27 nodes, 66 directed edges) and closely approximates its
summary statistics (e.g., density = 0.094 and average degree = 4.89 match
exactly, since these follow directly from node/edge counts; modularity Q ≈ 0.42
vs. the manuscript's 0.381), but individual confidence scores, evidence-type
labels, and literature references for interactions beyond the original top-10
are estimates, not the authors' original source values.

**Before public deposition, the corresponding author should replace
`ppi_edges.csv` (and, if available, `ppi_nodes.csv` molecular weight and
UniProt accession fields, currently marked `TBD`) with the original files used
to generate the manuscript's figures and tables.** This reconstruction is
provided as a working, reproducible starting point and a template for the
repository structure — not as a substitute for the authors' primary dataset.
