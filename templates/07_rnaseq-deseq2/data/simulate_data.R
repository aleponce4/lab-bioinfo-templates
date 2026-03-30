# simulate_data.R
# Builds synthetic inputs for template 07 in both supported modes:
# 1. count matrix + metadata
# 2. SummarizedExperiment RDS files

set.seed(42)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(writexl)
  library(SummarizedExperiment)
  library(S4Vectors)
})

gene_sets <- list(
  Interferon_Response = c("Stat1", "Stat2", "Irf7", "Isg15", "Ifit1", "Ifit3", "Mx1", "Oas1a"),
  Cell_Cycle = c("Mki67", "Top2a", "Cdk1", "Ccnb1", "Ccnb2", "Aurkb", "Ube2c", "Birc5"),
  Inflammation = c("Il6", "Ccl2", "Cxcl10", "Nfkbia", "Socs3", "Jun", "Fos", "Icam1")
)

required_genes <- unique(unlist(gene_sets, use.names = FALSE))
background_genes <- paste0("Gene_", sprintf("%03d", seq_len(260)))
gene_ids <- c(required_genes, background_genes)

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

pair_effects <- setNames(runif(length(unique(metadata$pair_id)), 0.9, 1.15), unique(metadata$pair_id))

simulate_gene_means <- function(sample_row) {
  mu <- rgamma(length(gene_ids), shape = 2.5, scale = 60)
  names(mu) <- gene_ids

  if (sample_row$tissue == "Brain") {
    mu[gene_sets$Interferon_Response] <- mu[gene_sets$Interferon_Response] * 1.4
    mu[gene_sets$Cell_Cycle] <- mu[gene_sets$Cell_Cycle] * 0.9
  } else {
    mu[gene_sets$Inflammation] <- mu[gene_sets$Inflammation] * 1.5
    mu[gene_sets$Cell_Cycle] <- mu[gene_sets$Cell_Cycle] * 1.2
  }

  if (sample_row$group == "Treatment") {
    if (sample_row$tissue == "Brain") {
      mu[gene_sets$Interferon_Response] <- mu[gene_sets$Interferon_Response] * 4.5
      mu[gene_sets$Inflammation] <- mu[gene_sets$Inflammation] * 2.0
    } else {
      mu[gene_sets$Inflammation] <- mu[gene_sets$Inflammation] * 4.0
      mu[gene_sets$Cell_Cycle] <- mu[gene_sets$Cell_Cycle] * 2.5
    }
  }

  mu * pair_effects[[sample_row$pair_id]]
}

gene_counts <- sapply(seq_len(nrow(metadata)), function(i) {
  mu <- simulate_gene_means(metadata[i, ])
  rnbinom(length(mu), mu = mu, size = 18)
})
colnames(gene_counts) <- sample_ids
rownames(gene_counts) <- gene_ids
storage.mode(gene_counts) <- "integer"

viral_features <- c("VEEV_genome", "VEEV_49S", "VEEV_26S")
viral_lengths <- c(VEEV_genome = 11446, VEEV_49S = 11446, VEEV_26S = 3880)

transcript_counts <- sapply(seq_len(nrow(metadata)), function(i) {
  sample_row <- metadata[i, ]
  base_signal <- if (sample_row$group == "Treatment") {
    if (sample_row$tissue == "Brain") 4000 else 1800
  } else {
    25
  }

  c(
    VEEV_genome = rnbinom(1, mu = base_signal, size = 12),
    VEEV_49S = rnbinom(1, mu = base_signal * 0.9, size = 12),
    VEEV_26S = rnbinom(1, mu = base_signal * 1.7, size = 12)
  )
})
colnames(transcript_counts) <- sample_ids
rownames(transcript_counts) <- viral_features
storage.mode(transcript_counts) <- "integer"

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

message("Wrote demo_gene_counts.csv, demo_transcript_counts.csv, demo_metadata.csv/xlsx, demo_gene_se.rds, and demo_transcript_se.rds")
