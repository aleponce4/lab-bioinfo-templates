# simulate_data.R
# Generates small synthetic DESeq2-style inputs for template 10.
# Outputs:
#   data/example_deg_results.csv
#   data/example_gene_modules.csv
# Run from the template folder: Rscript data/simulate_data.R

set.seed(101)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

dir.create("data", showWarnings = FALSE, recursive = TRUE)

comparisons <- c("Treatment_24h_vs_Control", "Treatment_48h_vs_Control", "Recovery_96h_vs_Control")

modules <- tribble(
  ~module, ~gene,
  "Antiviral_Response", "Ifit1",
  "Antiviral_Response", "Isg15",
  "Antiviral_Response", "Mx1",
  "Antiviral_Response", "Oas1",
  "Antiviral_Response", "Cxcl10",
  "Antiviral_Response", "Irf7",
  "Inflammatory_Signaling", "Tnf",
  "Inflammatory_Signaling", "Ccl2",
  "Inflammatory_Signaling", "Ptgs2",
  "Inflammatory_Signaling", "Icam1",
  "Inflammatory_Signaling", "Vcam1",
  "Vascular_Stress", "Vegfa",
  "Vascular_Stress", "Nos2",
  "Vascular_Stress", "Angpt2",
  "Vascular_Stress", "Thbd",
  "ECM_Remodeling", "Mmp3",
  "ECM_Remodeling", "Mmp9",
  "ECM_Remodeling", "Col4a1",
  "ECM_Remodeling", "Spp1",
  "Leukocyte_Infiltration", "Itgam",
  "Leukocyte_Infiltration", "Lyz2",
  "Leukocyte_Infiltration", "S100a9"
)

background <- c(
  "Gapdh", "Actb", "Rpl13a", "Hprt", "B2m", "Stat1", "Irf1", "Irf7",
  "Nfkb1", "Rela", "Jun", "Fos", "Hif1a", "Spi1", "Ets1", "Klf2",
  "Mapk1", "Ikbkb", "Jak1", "Pik3ca", "Tbk1", "Src", "Mapk14",
  "Ddx58", "Ifih1", "Socs3", "Cxcl1", "Cxcl2", "Cxcl9", "Ccl5"
)

genes <- unique(c(modules$gene, background))

spikes <- list(
  Treatment_24h_vs_Control = c("Ifit1", "Isg15", "Mx1", "Oas1", "Cxcl10", "Irf7"),
  Treatment_48h_vs_Control = c("Tnf", "Ccl2", "Ptgs2", "Icam1", "Vcam1", "Nos2", "Vegfa", "Mmp3", "Mmp9", "Itgam", "Lyz2", "S100a9"),
  Recovery_96h_vs_Control = c("Mmp3", "Mmp9", "Col4a1", "Spp1", "Angpt2", "Thbd")
)

rows <- lapply(comparisons, function(comp) {
  n <- length(genes)
  lfc <- rnorm(n, 0, 0.45)
  p <- runif(n, 0.20, 0.99)

  idx <- match(spikes[[comp]], genes)
  idx <- idx[!is.na(idx)]
  lfc[idx] <- rnorm(length(idx), 2.2, 0.35)
  p[idx] <- runif(length(idx), 1e-8, 5e-4)

  # Include a small number of down-regulated genes so the measured color scale is visible.
  down_genes <- if (comp == "Recovery_96h_vs_Control") c("Ifit1", "Cxcl10", "Icam1") else c("Klf2", "Thbd")
  didx <- match(down_genes, genes)
  didx <- didx[!is.na(didx)]
  lfc[didx] <- rnorm(length(didx), -1.8, 0.25)
  p[didx] <- runif(length(didx), 1e-5, 1e-3)

  tibble(
    Gene = genes,
    Comparison = comp,
    log2FoldChange = round(lfc, 4),
    padj = signif(p.adjust(p, method = "BH"), 4),
    stat = round(lfc / runif(n, 0.18, 0.45), 4)
  )
}) %>% bind_rows()

readr::write_csv(rows, "data/example_deg_results.csv")
readr::write_csv(modules, "data/example_gene_modules.csv")

message("Wrote data/example_deg_results.csv (", nrow(rows), " rows)")
message("Wrote data/example_gene_modules.csv (", nrow(modules), " rows)")
