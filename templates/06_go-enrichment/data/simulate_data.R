# simulate_data.R
# Generates a synthetic DEG table for template 06 (GO Enrichment).
# Output: data/deg_results.csv
# Run: Rscript data/simulate_data.R

set.seed(42)
library(dplyr)

# ── Parameters ────────────────────────────────────────────────────────────────
n_genes      <- 800
comparisons  <- c("Condition_A", "Condition_B", "Condition_C")
gene_names   <- paste0("Gene", seq_len(n_genes))

# ── Simulate ──────────────────────────────────────────────────────────────────
rows <- list()
for (comp in comparisons) {
  lfc  <- rnorm(n_genes, 0, 1.5)
  pval <- runif(n_genes, 0, 1)
  # Make ~15% significant
  sig_idx <- sample(n_genes, round(n_genes * 0.15))
  lfc[sig_idx]  <- rnorm(length(sig_idx), sign(rnorm(length(sig_idx))) * 2.5, 0.5)
  pval[sig_idx] <- runif(length(sig_idx), 0, 0.01)

  rows[[comp]] <- data.frame(
    Gene           = gene_names,
    Comparison     = comp,
    log2FoldChange = round(lfc, 4),
    padj           = round(pmin(pval * n_genes, 1), 4)   # simple Bonferroni-like
  )
}

df <- do.call(rbind, rows)
write.csv(df, "data/deg_results.csv", row.names = FALSE)
message("Wrote data/deg_results.csv  (", nrow(df), " rows)")
