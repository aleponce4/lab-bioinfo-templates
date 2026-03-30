# install_packages.R
# One-shot installer for all R packages used across all templates.
# Run once before using any template: source("install_packages.R")

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ── CRAN packages ─────────────────────────────────────────────────────────────
cran_pkgs <- c(
  # Core tidyverse
  "tidyverse", "readr", "readxl", "writexl", "janitor", "stringr",
  # Plotting
  "ggplot2", "ggpubr", "ggridges", "ggrepel", "ggbeeswarm", "ggh4x", "hexbin",
  "viridis", "RColorBrewer", "scales", "patchwork", "svglite", "ragg",
  "scatterpie", "ggspatial", "pheatmap",
  # Statistics
  "car", "broom", "rstatix", "emmeans", "coin", "drc", "vegan",
  # Spatial / maps
  "sf", "rnaturalearth", "rnaturalearthdata", "gtools", "ape",
  # Misc
  "knitr", "DT", "gt", "here", "forcats", "geosphere", "survival"
)

to_install <- cran_pkgs[!cran_pkgs %in% installed.packages()[, "Package"]]
if (length(to_install)) {
  message("Installing CRAN packages: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
} else {
  message("All CRAN packages already installed.")
}

# ── Bioconductor packages ─────────────────────────────────────────────────────
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Use Posit Public Package Manager for pre-built Bioconductor binaries (CI speed)
options(BioC_mirror = "https://packagemanager.posit.co/bioconductor")

bioc_pkgs <- c(
  "Biostrings",          # 05_phylo-geographic
  "DESeq2",              # 07_rnaseq-deseq2
  "SummarizedExperiment",
  "AnnotationDbi",
  "org.Mm.eg.db",        # mouse (default); swap to org.Hs.eg.db for human
  "org.Hs.eg.db",        # human
  "clusterProfiler",     # 06_go-enrichment
  "enrichplot",
  "ggtree",              # transitive dep of enrichplot; explicit pin ensures ggplot2-compatible version
  "GOplot",
  "ggVennDiagram",
  "ReactomePA"
)

bioc_to_install <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
if (length(bioc_to_install)) {
  message("Installing Bioconductor packages: ", paste(bioc_to_install, collapse = ", "))
  BiocManager::install(bioc_to_install, ask = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

message("\nDone! All packages installed.")
