# install_packages.R
# One-shot installer for all R packages used across all templates.
# Run once before using any template: source("install_packages.R")

# Repos configured by r-lib/actions/setup-r (use-public-rspm: true) or CRAN fallback

# ── CRAN packages ─────────────────────────────────────────────────────────────
# Keep this list to what the templates actually use — `Rscript scripts/validate_packages.R`
# reports anything referenced but undeclared, and anything declared but unreferenced.
cran_pkgs <- c(
  # Core tidyverse
  "tidyverse", "readr", "readxl", "writexl", "janitor", "stringr",
  # Plotting
  "ggplot2", "ggpubr", "ggridges", "ggbeeswarm",
  "viridis", "RColorBrewer", "scales", "patchwork", "svglite", "ragg",
  "scatterpie", "ggspatial", "pheatmap",
  "gprofiler2", "igraph",
  # Statistics
  "car", "broom", "rstatix", "emmeans", "vegan", "mclust", "minpack.lm",
  "ashr",                  # 07_rnaseq-deseq2: lfcShrink(type = "ashr")
  # Spatial / maps
  "sf", "rnaturalearth", "rnaturalearthdata", "gtools", "ape",
  # Misc
  "knitr", "DT", "gt", "forcats", "geosphere"
)

# Rscript sets no CRAN mirror, so install.packages() fails non-interactively without this.
if (is.null(getOption("repos")[["CRAN"]]) || getOption("repos")[["CRAN"]] %in% c("@CRAN@", NA)) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
}

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

# Bioconductor mirror set by r-lib/actions/setup-r or defaults

bioc_pkgs <- c(
  "Biostrings",          # 05_phylo-geographic
  "DESeq2",              # 07_rnaseq-deseq2
  "SummarizedExperiment",
  "AnnotationDbi",
  "org.Mm.eg.db",        # mouse (default); swap to org.Hs.eg.db for human
  "org.Hs.eg.db",        # human
  "clusterProfiler",     # 06_go-enrichment
  "enrichplot",
  "impute",                # WGCNA Bioconductor dep (not on CRAN)
  "preprocessCore",        # WGCNA Bioconductor dep
  "GO.db",                 # WGCNA Bioconductor dep
  "WGCNA"                  # 09_wgcna (installed after its Bioc deps)
  # NOTE: decoupleR / OmnipathR / CARNIVAL are deliberately NOT listed. Template 10
  # is a self-contained toy demo built on hardcoded CollecTRI-style and
  # CARNIVAL-style networks and never calls those packages, so installing them
  # would add three heavy dependencies for zero executed code.
)

bioc_to_install <- bioc_pkgs[!bioc_pkgs %in% installed.packages()[, "Package"]]
if (length(bioc_to_install)) {
  message("Installing Bioconductor packages: ", paste(bioc_to_install, collapse = ", "))
  BiocManager::install(bioc_to_install, ask = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

# install.packages()/BiocManager::install() only warn on a failed build, so verify
# rather than reporting success unconditionally.
still_missing <- setdiff(c(cran_pkgs, bioc_pkgs), rownames(installed.packages()))
if (length(still_missing)) {
  stop("Failed to install: ", paste(still_missing, collapse = ", "),
       "\nCheck the build log above for missing system libraries.")
}

message("\nDone! All packages installed.")
