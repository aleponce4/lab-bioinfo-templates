# simulate_data.R
# Generates a synthetic DEG table for template 06 (GO Enrichment).
# Output: data/deg_results.csv
# Run: Rscript data/simulate_data.R

set.seed(42)

# Ensure working directory is set to template directory
if (!file.exists("template.qmd")) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    setwd(dirname(script_dir))
  }
}
dir.create("data", showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({
  library(dplyr)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

comparisons <- c("Condition_A", "Condition_B", "Condition_C")

# These pathway-focused symbols are valid mouse genes that map through org.Mm.eg.db
# and are coherent enough to yield visible GO and KEGG enrichment in the demo.
spike_sets <- list(
  Condition_A = c(
    "Stat1", "Stat2", "Irf7", "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1",
    "Oas1a", "Oas2", "Rsad2", "Ddx58", "Ifih1", "Usp18", "Bst2"
  ),
  Condition_B = c(
    "Mki67", "Top2a", "Cdk1", "Ccnb1", "Ccnb2", "Aurkb", "Bub1", "Ube2c",
    "Cenpf", "Nusap1", "Tpx2", "Pttg1", "Cdc20", "Plk1", "Kif11"
  ),
  Condition_C = c(
    "Il6", "Cxcl10", "Ccl2", "Ccl5", "Tnf", "Nfkbia", "Socs3", "Icam1",
    "Jun", "Fos", "Stat3", "Ptgs2", "Cxcl1", "Il1b", "Tlr2"
  )
)

# One direction per pathway, fixed here so the demo tells a consistent story:
# an induced antiviral/interferon programme, a suppressed proliferation
# programme, and an induced inflammatory programme.
pathway_direction <- list(
  Condition_A =  1,   # interferon-stimulated genes: up
  Condition_B = -1,   # cell-cycle / mitotic genes: down
  Condition_C =  1    # NF-kB / inflammatory genes: up
)
stopifnot(setequal(names(pathway_direction), comparisons))

required_genes <- unique(unlist(spike_sets, use.names = FALSE))

all_symbols <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL")
valid_symbols <- unique(all_symbols[!is.na(all_symbols) & all_symbols != ""])
valid_symbols <- valid_symbols[grepl("^[A-Za-z][A-Za-z0-9.-]+$", valid_symbols)]

background_pool <- setdiff(valid_symbols, required_genes)

# The template uses the full tested gene list as the enrichment BACKGROUND
# (universe), and applies MIN_GENE_SET/MAX_GENE_SET to gene sets after they are
# intersected with it. A few hundred genes is therefore too small to be a
# realistic demo: nearly every GO term and KEGG pathway falls below the minimum
# set size and enrichment returns nothing. A real DESeq2 result carries
# thousands of tested genes, so the demo does too.
N_BACKGROUND_GENES <- 2000
background_genes <- sample(background_pool, N_BACKGROUND_GENES)
gene_names <- c(required_genes, background_genes)
n_genes <- length(gene_names)

rows <- list()
for (comp in comparisons) {
  lfc <- rnorm(n_genes, 0, 0.9)
  pval <- runif(n_genes, 0.2, 1)

  # Add a modest background of unrelated DE genes so the examples do not look synthetic.
  bg_sig_idx <- sample(seq_along(gene_names), round(0.02 * n_genes))
  bg_sign <- sample(c(-1, 1), length(bg_sig_idx), replace = TRUE)
  lfc[bg_sig_idx] <- rnorm(length(bg_sig_idx), mean = 1.9 * bg_sign, sd = 0.25)
  pval[bg_sig_idx] <- runif(length(bg_sig_idx), 1e-3, 0.03)

  # Force a coherent pathway signal for each comparison.
  # The sign is drawn ONCE PER PATHWAY, not once per gene. Drawing it per gene
  # makes each "coherent pathway signal" move half up and half down at the same
  # time, which is what the interferon panel used to show.
  spike_idx  <- match(spike_sets[[comp]], gene_names)
  spike_sign <- pathway_direction[[comp]]
  lfc[spike_idx] <- rnorm(length(spike_idx), mean = 2.8 * spike_sign, sd = 0.3)
  pval[spike_idx] <- runif(length(spike_idx), 1e-10, 1e-5)

  rows[[comp]] <- data.frame(
    Gene = gene_names,
    Comparison = comp,
    log2FoldChange = round(lfc, 4),
    padj = p.adjust(pval, method = "BH")
  )
}

df <- bind_rows(rows)
write.csv(df, "data/deg_results.csv", row.names = FALSE)
message("Wrote data/deg_results.csv (", nrow(df), " rows)")
