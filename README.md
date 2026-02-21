# Functional Enrichment & GSEA App

A Streamlit-based interactive application for **KEGG pathway enrichment** (ORA & GSEA) and **COG functional category analysis** of bacterial RNA-seq data. Inspired by [pathsnake](https://github.com/RickGelhausen/pathsnake).

![Python](https://img.shields.io/badge/Python-3.10-blue)
![Streamlit](https://img.shields.io/badge/Streamlit-1.28+-red)
![License](https://img.shields.io/badge/License-MIT-green)

---

## Features

- **DESeq2 Results Parsing** — Loads tab/comma-separated DESeq2 output files with GFF-style annotations
- **KEGG Over-Representation Analysis (ORA)** — Fisher's exact test with Benjamini-Hochberg correction for up/down/all DE genes
- **Gene Set Enrichment Analysis (GSEA)** — Prerank-based GSEA using [GSEApy](https://github.com/zqfang/GSEApy)
- **COG Category Analysis** — COG enrichment and distribution (inferred from product descriptions, fetched from KEGG, or uploaded)
- **Interactive Plots** — Volcano, MA, dot plots, bar charts, lollipop, waterfall, treemap, heatmap, and more (Plotly)
- **Cross-strain ID Mapping** — Maps gene names/symbols across KEGG organism databases (e.g., RefSeq locus tags → KEGG gene IDs)
- **Downloadable Results** — Export all analysis tables as TSV

## Input Format

The app accepts DESeq2 results files (TSV or CSV) containing:

| Column | Description |
|--------|-------------|
| `log2FoldChange` | Log2 fold change |
| `padj` | Adjusted p-value |
| `baseMean` | Mean expression |
| `Attributes` | GFF-style annotation string with `locus_tag=`, `gene=`, `product=`, etc. |

An example file is included: `deseq_comp_InSPI2_vs_LSP_with_annotation_and_countings.csv`

## Quick Start

### 1. Create the conda environment

```bash
conda env create -f environment.yml
```

### 2. Activate and run

```bash
conda activate enrichment_app
streamlit run app.py --server.port 8501
```

Or use the launcher script:

```bash
bash run_app.sh
```

### 3. Open in browser

Navigate to [http://localhost:8501](http://localhost:8501)

## Usage

1. **Upload or load** your DESeq2 results file in the sidebar (or use the built-in example)
2. **Set the KEGG organism code** (e.g., `sey` for *Salmonella enterica* SL1344, `eco` for *E. coli* K-12)
3. **Adjust parameters** — p-value cutoff, log2FC cutoff, min gene set size, GSEA permutations
4. **Click "Run Enrichment Analysis"**
5. **Explore results** across five tabs:
   - 🌋 **Overview** — Volcano & MA plots, DE gene summary
   - 🧬 **KEGG Pathways** — ORA dot/bar/lollipop/treemap plots
   - 📈 **GSEA** — NES dot/bar/waterfall plots
   - 🏷️ **COG Analysis** — Category distribution, enrichment, heatmap
   - 📊 **Summary** — Combined KEGG + COG overview

## Project Structure

```
functional_enrichment_app/
├── app.py                  # Main Streamlit application
├── analysis/
│   ├── __init__.py
│   ├── data_parser.py      # DESeq2 file parsing & gene ID extraction
│   ├── kegg_analysis.py    # KEGG API queries, ORA, GSEA
│   ├── cog_analysis.py     # COG category mapping & enrichment
│   └── plotting.py         # All Plotly visualization functions
├── environment.yml         # Conda environment specification
├── run_app.sh              # Launch script
└── README.md
```

## KEGG Organism Codes

Common codes for Salmonella strains:

| Code | Strain |
|------|--------|
| `sey` | *S. enterica* serovar Typhimurium **SL1344** |
| `stm` | *S. enterica* serovar Typhimurium **LT2** |
| `seb` | *S. enterica* serovar Typhimurium **ST4/74** |
| `seo` | *S. enterica* serovar Typhimurium **14028S** |
| `eco` | *E. coli* K-12 MG1655 |

Find your organism code: [KEGG Organism List](https://www.genome.jp/kegg/catalog/org_list.html)

## Dependencies

- Python 3.10
- Streamlit ≥ 1.28
- Pandas, NumPy, SciPy
- Plotly ≥ 5.15
- GSEApy ≥ 1.0
- statsmodels ≥ 0.14
- requests

## Acknowledgments

- Analysis approach inspired by [pathsnake](https://github.com/RickGelhausen/pathsnake)
- KEGG pathway data from [KEGG REST API](https://rest.kegg.jp)
- GSEA powered by [GSEApy](https://github.com/zqfang/GSEApy)
