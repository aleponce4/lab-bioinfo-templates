# simulate_data.R
# Builds realistic synthetic inputs for template 07 (RNA-seq / DESeq2):
# - Realistic count distribution (3,536 genes spanning 10 - 50,000 counts)
# - Realistic negative binomial mean-dispersion trend matching DESeq2 parametric model
# - Classic volcano distribution with up-regulated, down-regulated, and non-DE background genes
# Outputs:
# 1. demo_gene_counts.csv / demo_transcript_counts.csv
# 2. demo_metadata.csv / demo_metadata.xlsx
# 3. demo_gene_se.rds / demo_transcript_se.rds

# Ensure working directory is set to template directory
if (!file.exists("template.qmd")) {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[1])))
    setwd(dirname(script_dir))
  } else if (dir.exists("templates/07_rnaseq-deseq2")) {
    setwd("templates/07_rnaseq-deseq2")
  }
}
dir.create("data", showWarnings = FALSE, recursive = TRUE)

set.seed(42)


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(writexl)
  library(SummarizedExperiment)
  library(S4Vectors)
})

# Curated biological pathway gene sets
gene_sets <- list(
  Interferon_Response = c("Stat1", "Stat2", "Irf7", "Isg15", "Ifit1", "Ifit2", "Ifit3", "Mx1", "Oas1a", "Oas2", "Irf9", "Usp18", "B2m"),
  Cell_Cycle = c("Mki67", "Top2a", "Cdk1", "Ccnb1", "Ccnb2", "Aurkb", "Ube2c", "Birc5", "Plk1", "Cdc20", "Cdk2"),
  Inflammation = c("Il6", "Ccl2", "Cxcl10", "Nfkbia", "Socs3", "Jun", "Fos", "Icam1", "Tnf", "Il1b", "Cxcl9", "Ccl5")
)

pathway_genes <- unique(unlist(gene_sets, use.names = FALSE))
n_bg <- 3500
bg_genes <- paste0("Gene_", sprintf("%04d", seq_len(n_bg)))
all_genes <- c(pathway_genes, bg_genes)
n_genes <- length(all_genes)

# Samples metadata
metadata <- tidyr::expand_grid(
  tissue = c("Brain", "Lung"),
  pair_id = paste0("Pair_", 1:4),
  group = c("Control", "Treatment")
) %>%
  arrange(tissue, pair_id, group) %>%
  mutate(
    sample_id = paste(tissue, pair_id, group, sep = "_"),
    alternate_id = paste0("s", sprintf("%02d", row_number())),
    batch = rep(c("Batch_1", "Batch_2"), length.out = n()),
    row_id = sample_id
  )

sample_ids <- metadata$sample_id

# 1. Base gene means: log-normal distribution (dynamic range: ~10 to 50,000 counts)
log_means <- rnorm(n_genes, mean = 4.8, sd = 1.8)
base_means <- exp(log_means)
base_means <- pmax(base_means, 10)
names(base_means) <- all_genes

# 2. Mean-dispersion relationship (DESeq2 parametric curve: alpha = alpha0 + alpha1 / mean)
alpha0 <- 0.03
alpha1 <- 3.5
log_disp_noise <- rnorm(n_genes, mean = 0, sd = 0.25)
dispersions <- (alpha0 + alpha1 / base_means) * exp(log_disp_noise)
dispersions <- pmax(dispersions, 0.005)
names(dispersions) <- all_genes

# 3. Treatment Fold Changes (log2 scale)
l2fc <- numeric(n_genes)
names(l2fc) <- all_genes

# Pathway genes
l2fc[gene_sets$Interferon_Response] <- runif(length(gene_sets$Interferon_Response), 1.8, 4.0)
l2fc[gene_sets$Inflammation]         <- runif(length(gene_sets$Inflammation), 1.5, 3.8)
l2fc[gene_sets$Cell_Cycle]          <- runif(length(gene_sets$Cell_Cycle), -2.8, -1.0)

# Background genes DE assignment (~11% up, ~11% down, ~78% unchanged)
# NOTE: no second set.seed() here. Re-seeding mid-script restarts the RNG stream
# and makes everything after this point independent of the set.seed(42) above,
# which quietly breaks the "one seed reproduces the whole file" contract.
bg_indices <- seq(length(pathway_genes) + 1, n_genes)
n_bg_de <- round(length(bg_indices) * 0.11)
up_bg <- sample(bg_indices, n_bg_de)
down_bg <- sample(setdiff(bg_indices, up_bg), n_bg_de)

l2fc[up_bg]   <- runif(length(up_bg), 0.8, 3.2)
l2fc[down_bg] <- runif(length(down_bg), -3.2, -0.8)

# Unchanged genes get realistic background biological noise
null_bg <- setdiff(bg_indices, c(up_bg, down_bg))
l2fc[null_bg] <- rnorm(length(null_bg), mean = 0, sd = 0.12)

# Pair/Tissue effects per sample
pair_effects <- setNames(exp(rnorm(4, 0, 0.08)), paste0("Pair_", 1:4))
tissue_effects <- setNames(exp(rnorm(2, 0, 0.10)), c("Brain", "Lung"))

# 4. Generate Count Matrix
gene_counts <- matrix(0L, nrow = n_genes, ncol = nrow(metadata),
                      dimnames = list(all_genes, sample_ids))

for (i in seq_len(nrow(metadata))) {
  row <- metadata[i, ]
  
  sample_l2fc <- if (row$group == "Treatment") l2fc else numeric(n_genes)
  
  if (row$tissue == "Brain" && row$group == "Treatment") {
    sample_l2fc[gene_sets$Interferon_Response] <- sample_l2fc[gene_sets$Interferon_Response] + 0.5
  }
  
  sample_means <- base_means * (2^sample_l2fc) * pair_effects[row$pair_id] * tissue_effects[row$tissue]
  sizes <- 1 / dispersions
  gene_counts[, i] <- rnbinom(n_genes, mu = sample_means, size = sizes)
}

storage.mode(gene_counts) <- "integer"

# 5. Simulate viral features (transcript counts)
viral_features <- c("Pathogen_genome", "Pathogen_gRNA", "Pathogen_sgRNA")
viral_lengths <- c(Pathogen_genome = 11500, Pathogen_gRNA = 11500, Pathogen_sgRNA = 3900)

transcript_counts <- sapply(seq_len(nrow(metadata)), function(i) {
  sample_row <- metadata[i, ]
  base_signal <- if (sample_row$group == "Treatment") {
    if (sample_row$tissue == "Brain") 4000 else 1800
  } else {
    25
  }

  c(
    Pathogen_genome = rnbinom(1, mu = base_signal, size = 12),
    Pathogen_gRNA = rnbinom(1, mu = base_signal * 0.9, size = 12),
    Pathogen_sgRNA = rnbinom(1, mu = base_signal * 1.7, size = 12)
  )
})
colnames(transcript_counts) <- sample_ids
rownames(transcript_counts) <- viral_features
storage.mode(transcript_counts) <- "integer"

# 6. Save output files
gene_count_tbl <- tibble::tibble(Gene = rownames(gene_counts)) %>%
  bind_cols(as.data.frame(gene_counts, check.names = FALSE))
transcript_count_tbl <- tibble::tibble(Transcript = rownames(transcript_counts)) %>%
  bind_cols(as.data.frame(transcript_counts, check.names = FALSE))

write.csv(gene_count_tbl, "data/demo_gene_counts.csv", row.names = FALSE)
write.csv(transcript_count_tbl, "data/demo_transcript_counts.csv", row.names = FALSE)
write.csv(metadata %>% select(-row_id), "data/demo_metadata.csv", row.names = FALSE)
writexl::write_xlsx(list(metadata = metadata %>% select(-row_id)), "data/demo_metadata.xlsx")

metadata_se <- as.data.frame(metadata %>% select(-row_id), stringsAsFactors = FALSE)
rownames(metadata_se) <- metadata_se$sample_id

gene_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = gene_counts),
  colData = S4Vectors::DataFrame(metadata_se)
)

transcript_se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(counts = transcript_counts),
  rowData = S4Vectors::DataFrame(feature_length_bp = viral_lengths[rownames(transcript_counts)]),
  colData = S4Vectors::DataFrame(metadata_se)
)

saveRDS(gene_se, "data/demo_gene_se.rds")
saveRDS(transcript_se, "data/demo_transcript_se.rds")

message("Wrote demo_gene_counts.csv, demo_transcript_counts.csv, demo_metadata.csv/xlsx, demo_gene_se.rds, and demo_transcript_se.rds (", n_genes, " genes)")
