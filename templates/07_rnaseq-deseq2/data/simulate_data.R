# simulate_data.R
# Creates a synthetic DESeqDataSet (dds_template.RData) for template 07.
# Run from template folder: Rscript data/simulate_data.R

set.seed(42)

library(DESeq2)

# ── Parameters ────────────────────────────────────────────────────────────────
n_genes      <- 500
n_samples    <- 12   # 4 groups × 3 replicates
group_labels <- rep(c("Control", "Group_A", "Group_B", "Group_C"), each = 3)
sample_names <- paste0(group_labels, "_rep", rep(1:3, 4))
gene_names   <- paste0("Gene", seq_len(n_genes))

# ── Simulate count matrix ──────────────────────────────────────────────────────
base_means <- round(runif(n_genes, 50, 2000))
counts <- matrix(
  sapply(base_means, function(m) rnbinom(n_samples, mu = m, size = 10)),
  nrow = n_genes, ncol = n_samples,
  dimnames = list(gene_names, sample_names)
)

# Introduce some DE: 15% of genes upregulated in Group_A/B/C vs Control
de_idx <- sample(n_genes, round(n_genes * 0.15))
for (grp_col in which(group_labels != "Control")) {
  fold <- sample(c(3, 5, 8, 0.3, 0.2), length(de_idx), replace = TRUE)
  counts[de_idx, grp_col] <- round(counts[de_idx, grp_col] * fold)
}

# ── Build DESeqDataSet ────────────────────────────────────────────────────────
col_data <- data.frame(
  group   = factor(group_labels, levels = c("Control", "Group_A", "Group_B", "Group_C")),
  timepoint = rep(c("0h", "24h", "48h", "72h"), each = 3),
  row.names = sample_names
)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = col_data,
  design    = ~ group
)
dds <- DESeq(dds)

save(dds, file = "data/dds_template.RData")
message("Wrote data/dds_template.RData")
