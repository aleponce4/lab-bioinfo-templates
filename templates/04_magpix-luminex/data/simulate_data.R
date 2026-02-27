# simulate_data.R
# Generates synthetic MagPix/Luminex multiplex cytokine data for template 04.
# Output: data/luminex_data.csv
# Run from template folder: Rscript data/simulate_data.R

set.seed(42)

# ── Parameters ────────────────────────────────────────────────────────────────
analytes  <- c("CCL2", "CCL5", "IL-6", "IL-10", "TNF-alpha",
               "IFN-gamma", "CXCL10", "IL-1beta")
factor1   <- c("PBS", "Disease")    # two levels of Factor 1
factor2   <- c("Vehicle", "Drug")   # two levels of Factor 2
n_per_group <- 5
lod       <- 2.0  # pg/mL, limit of detection

# True mean concentrations on log10 scale (PBS+Vehicle as baseline)
means <- list(
  CCL2      = c(2.5, 2.6, 3.2, 3.1),
  CCL5      = c(2.2, 2.3, 2.8, 2.5),
  `IL-6`    = c(1.8, 1.9, 3.5, 3.0),
  `IL-10`   = c(1.5, 1.6, 2.4, 2.0),
  `TNF-alpha` = c(1.2, 1.3, 2.8, 2.4),
  `IFN-gamma` = c(1.0, 1.1, 3.0, 2.5),
  CXCL10    = c(2.0, 2.1, 3.8, 3.2),
  `IL-1beta`= c(1.3, 1.4, 2.6, 2.1)
)
sd_val <- 0.25

# ── Generate data ─────────────────────────────────────────────────────────────
samples <- expand.grid(
  F1 = factor1,
  F2 = factor2,
  rep = seq_len(n_per_group)
)
samples$Sample <- paste0(samples$F1, "_", samples$F2, "_", samples$rep)

group_order <- c("PBS+Vehicle", "PBS+Drug", "Disease+Vehicle", "Disease+Drug")

rows <- list()
for (i in seq_len(nrow(samples))) {
  grp_idx <- which(group_order == paste0(samples$F1[i], "+", samples$F2[i]))
  row <- data.frame(
    Location    = paste0("Row", ceiling(i / 4)),
    Sample      = samples$Sample[i],
    Original_ID = paste0("ID_", i),
    `Total Events` = sample(80:120, 1),
    check.names = FALSE
  )
  for (an in analytes) {
    true_val <- 10^rnorm(1, means[[an]][grp_idx], sd_val)
    if (true_val < lod) {
      row[[an]] <- paste0("< ", lod)  # censored
    } else {
      row[[an]] <- round(true_val, 2)
    }
  }
  rows[[i]] <- row
}

df <- do.call(rbind, rows)
write.csv(df, "data/luminex_data.csv", row.names = FALSE)
message("Wrote data/luminex_data.csv  (", nrow(df), " rows, ",
        length(analytes), " analytes)")
