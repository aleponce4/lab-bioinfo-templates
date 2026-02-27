# Jonsson Lab Bioinformatics Templates

A centralized collection of reusable analysis templates for virology and genomics research.
Browse the live gallery at: **https://jonsson-lab.github.io/lab-bioinfo-templates**

## Templates

| # | Template | Language | Key Analysis |
|---|----------|----------|-------------|
| 01 | Plaque Assay + Violin Plots | R | One-way ANOVA, Tukey HSD, log10 titer |
| 02 | Image Infection + Dose-Response | R | 4PL regression, EC50/CC50 |
| 03 | Multi-Condition ANOVA | R | Welch ANOVA, Holm correction |
| 04 | MagPix / Luminex Multiplex | R | Two-way ANOVA, emmeans, FDR |
| 05 | Phylo-Geographic Analysis | R | Mantel test, scatter pies, maps |
| 06 | GO Enrichment Analysis | R | clusterProfiler, KEGG/Reactome |
| 07 | RNA-seq / DESeq2 | R | Volcano, bubble, PCA, triangle plots |
| 08 | VCF Mutation Analysis | Python | Entropy, mutation frequency, heatmaps |

## Quick Start

### R Templates (01–07)

```r
# 1. Install all R packages once
source("install_packages.R")

# 2. Open a template and edit the USER CONFIGURATION block at the top
# 3. Knit / render
```

### Python Template (08)

```bash
conda env create -f environment.yml
conda activate jonsson-bioinfo
jupyter notebook templates/08_vcf-mutation-analysis/template.ipynb
```

## How Templates Work

Every template has a clearly marked configuration section at the very top:

```r
## ── USER CONFIGURATION ──────────────────────────────────────────────────────
DATA_FILE   <- "data/your_data.csv"   # <-- change this
GROUP_COL   <- "Treatment"
# ... (all parameters documented inline)
## ────────────────────────────────────────────────────────────────────────────
```

**Only edit the config block.** The analysis code below uses those variables and should
not need modification for standard use cases.

## Synthetic Example Data

Each template ships with a `data/simulate_data.R` (or `.py`) script that generates
a small, realistic synthetic dataset. These are the datasets used by the Quarto gallery
site so every template renders a meaningful example plot.

To regenerate example data:
```r
source("templates/01_plaque-assay-violin/data/simulate_data.R")
```

## Plot Style

All templates preserve the **publication-quality plot style** established in the lab:
- `theme_classic()` or `theme_bw()` base themes
- Arial font family throughout
- Okabe-Ito colorblind-safe palettes where applicable
- TIFF / SVG / PDF output at 300–600 dpi
- Consistent axis tick lengths and line widths

Do **not** modify the `theme()` calls unless intentionally changing the style for a
specific figure.

## Contributing

1. Fork the repo
2. Copy an existing template folder as a starting point
3. Follow the conversion rules in the template header comments
4. Add a row to `index.qmd`'s gallery table
5. Open a pull request

## License

MIT — free to use and adapt for academic research.
