# simulate_data.R
# Generates synthetic data for template 09 (WGCNA Co-Expression):
#   - data/dds_wgcna_template.RData  (DESeqDataSet, post DESeq(), ready to load)
# Run: Rscript data/simulate_data.R

set.seed(42)

library(DESeq2)
library(SummarizedExperiment)

# ── Design ──────────────────────────────────────────────────────────────────
# 18 samples: 3 DPI (one / three / five) × 2 conditions (virus / pbs) × 3 reps
# 600 host genes arranged in 5 co-expression modules (~80 genes each) + 100 background
# 1 viral gene (SIM_VIRUS) — high in virus samples
#
# Co-expression structure:
#   Module 1 (turquoise, 100 genes): strongly up in virus, rises with DPI
#   Module 2 (blue,      100 genes): down in virus, anti-correlated with DPI
#   Module 3 (red,        80 genes): rises with DPI regardless of treatment (time effect)
#   Module 4 (green,      80 genes): slight up in virus, stable over time
#   Module 5 (yellow,     80 genes): noisy, weak signal
#   Background (grey,    160 genes): no consistent pattern

n_samples <- 18
dpi_levels  <- c("one", "three", "five")
trt_levels  <- c("virus", "pbs")
n_reps      <- 3

sample_meta <- expand.grid(
  dpi       = dpi_levels,
  treatment = trt_levels,
  rep       = seq_len(n_reps),
  stringsAsFactors = FALSE
)
sample_meta$group <- paste0(sample_meta$treatment, sample_meta$dpi)
sample_meta$rep   <- NULL
rownames(sample_meta) <- paste0("S", seq_len(nrow(sample_meta)))

dpi_num <- c(one = 1, three = 3, five = 5)
sample_meta$dpi_n     <- dpi_num[sample_meta$dpi]
sample_meta$is_virus  <- as.integer(sample_meta$treatment == "virus")

# ── Eigengene templates (one value per sample) ──────────────────────────────
eig1 <- with(sample_meta,  is_virus * (dpi_n / 5))        # up in virus, rises with DPI
eig2 <- with(sample_meta, -is_virus * (dpi_n / 5))        # down in virus
eig3 <- with(sample_meta,  dpi_n / 5)                      # time effect only
eig4 <- with(sample_meta,  0.5 * is_virus)                 # modest virus effect
eig5 <- rnorm(n_samples, 0, 0.3)                           # weak/noisy

# ── Build count matrix ───────────────────────────────────────────────────────
make_module_counts <- function(eigengene, n_genes, base_mean = 200, noise_sd = 0.3) {
  loadings <- runif(n_genes, 0.5, 1.5)
  mat <- sapply(seq_len(n_samples), function(i) {
    mu <- base_mean * exp(loadings * eigengene[i] + rnorm(n_genes, 0, noise_sd))
    rpois(n_genes, lambda = mu)
  })
  rownames(mat) <- paste0(substitute(eigengene), "_G", seq_len(n_genes))
  mat
}

# Give modules distinct gene name prefixes
make_mod <- function(prefix, eigengene, n_genes, base_mean = 200) {
  loadings <- runif(n_genes, 0.5, 1.5)
  mat <- sapply(seq_len(n_samples), function(i) {
    mu <- base_mean * exp(loadings * eigengene[i] + rnorm(n_genes, 0, 0.25))
    rpois(n_genes, lambda = pmax(mu, 1))
  })
  rownames(mat) <- paste0(prefix, "_G", seq_len(n_genes))
  mat
}

m1  <- make_mod("Turq",   eig1, 100, 800)
m2  <- make_mod("Blue",   eig2, 100, 700)
m3  <- make_mod("Red",    eig3,  80, 600)
m4  <- make_mod("Green",  eig4,  80, 500)
m5  <- make_mod("Yellow", eig5,  80, 400)

# Background genes: moderate counts so VST > 5 after normalization
bg_mat <- matrix(
  rpois(160 * n_samples, lambda = sample(300:800, 160, replace = TRUE)),
  nrow = 160, ncol = n_samples
)
rownames(bg_mat) <- paste0("BG_G", seq_len(160))

# Viral gene: high counts in virus, near-zero in PBS
virus_counts <- with(sample_meta,
  rpois(n_samples, lambda = ifelse(treatment == "virus", 500 + 300 * (dpi_n / 5), 2))
)
virus_mat <- matrix(virus_counts, nrow = 1)
rownames(virus_mat) <- "SIM_VIRUS"

counts_mat <- rbind(m1, m2, m3, m4, m5, bg_mat, virus_mat)
colnames(counts_mat) <- rownames(sample_meta)

# ── Build DESeqDataSet ────────────────────────────────────────────────────────
col_data <- sample_meta[, c("treatment", "dpi", "group")]
col_data$group <- factor(col_data$group)

dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = col_data,
  design    = ~ group
)

# Run DESeq so object is analysis-ready
dds <- DESeq(dds, quiet = TRUE)

# Save
dir.create("data", showWarnings = FALSE)
save(dds, file = "data/dds_wgcna_template.RData")
message("Wrote data/dds_wgcna_template.RData  (",
        nrow(dds), " genes × ", ncol(dds), " samples)")
