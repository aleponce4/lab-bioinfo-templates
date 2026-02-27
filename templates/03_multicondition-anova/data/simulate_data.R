# simulate_data.R
# Generates synthetic multi-condition infection data for template 03.
# Output: data/condition_A.csv, data/condition_B.csv, data/condition_C.csv
# Run from the template folder: Rscript data/simulate_data.R

set.seed(42)

out_dir <- "data/"   # script is always run from the template directory

# ── Parameters ────────────────────────────────────────────────────────────────
cell_types <- c("CellType_A", "CellType_B", "CellType_C", "CellType_D")
n_reps     <- 3       # technical replicates per cell type
conditions <- c("Condition_A", "Condition_B", "Condition_C")

# Mean log10 titer per condition × cell type (CellType_A is the "reference")
means <- list(
  Condition_A = c(CellType_A = 6.5, CellType_B = 4.2, CellType_C = 3.8, CellType_D = 5.1),
  Condition_B = c(CellType_A = 5.8, CellType_B = 3.5, CellType_C = 4.1, CellType_D = 4.7),
  Condition_C = c(CellType_A = 6.0, CellType_B = 3.9, CellType_C = 3.3, CellType_D = 5.3)
)
sd_val <- 0.35   # within-group SD on log10 scale

# ── Generate and save CSVs ────────────────────────────────────────────────────
for (cond in conditions) {
  m <- means[[cond]]
  # Wide format: each column = cell type, rows = replicates
  mat <- sapply(cell_types, function(ct) {
    vals <- 10^rnorm(n_reps, mean = m[ct], sd = sd_val)
    vals
  })
  # First replicate in plain column; subsequent as "CellType.1", "CellType.2"
  df <- as.data.frame(mat)
  # Add extra replicate columns (mimicking the original file structure)
  df2 <- df
  colnames(df2) <- cell_types
  # Append second replicate block as ".1" columns (original format convention)
  df3 <- df2
  colnames(df3) <- paste0(cell_types, ".1")
  out <- cbind(df2, df3)
  write.csv(out, file = paste0(out_dir, cond, ".csv"), row.names = FALSE)
  message("Wrote data/", cond, ".csv")
}
