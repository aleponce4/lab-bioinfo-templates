# simulate_data.R
# Generates synthetic multi-condition infection data for template 03.
# Output: data/condition_A.csv, data/condition_B.csv, data/condition_C.csv
# Run from the template folder: Rscript data/simulate_data.R

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
out_dir <- "data/"


# ── Parameters ────────────────────────────────────────────────────────────────
cell_types <- c("CellType_A", "CellType_B", "CellType_C", "CellType_D")
n_reps     <- 3       # replicates per cell type, per experimental block
n_blocks   <- 2       # independent experimental blocks written as ".1" columns
conditions <- c("Condition_A", "Condition_B", "Condition_C")

# Mean log10 titer per condition × cell type (CellType_A is the "reference")
means <- list(
  Condition_A = c(CellType_A = 6.5, CellType_B = 4.2, CellType_C = 3.8, CellType_D = 5.1),
  Condition_B = c(CellType_A = 5.8, CellType_B = 3.5, CellType_C = 4.1, CellType_D = 4.7),
  Condition_C = c(CellType_A = 6.0, CellType_B = 3.9, CellType_C = 3.3, CellType_D = 5.3)
)
sd_val <- 0.35   # within-group SD on log10 scale

# ── Generate and save CSVs ────────────────────────────────────────────────────
# Each experimental block is drawn INDEPENDENTLY. (An earlier version copied the
# first block into the ".1" columns, which doubled the apparent n without adding
# any information and made every p-value far too small.)
draw_block <- function(m) {
  as.data.frame(sapply(cell_types, function(ct) {
    10^rnorm(n_reps, mean = m[[ct]], sd = sd_val)
  }))
}

for (cond in conditions) {
  m <- means[[cond]]

  # Block 1 -> plain column names; block 2 -> ".1" suffix, block 3 -> ".2", ...
  blocks <- lapply(seq_len(n_blocks), function(b) {
    blk <- draw_block(m)
    colnames(blk) <- if (b == 1) cell_types else paste0(cell_types, ".", b - 1)
    blk
  })
  out <- do.call(cbind, blocks)

  # Guard against the copy-paste bug returning: no two replicate blocks may be
  # identical, and no value may be non-positive (the template log-transforms).
  if (n_blocks > 1) {
    for (b in 2:n_blocks) {
      stopifnot(!isTRUE(all.equal(unname(as.matrix(blocks[[1]])),
                                  unname(as.matrix(blocks[[b]])))))
    }
  }
  stopifnot(all(out > 0), !any(is.na(out)))

  write.csv(out, file = paste0(out_dir, cond, ".csv"), row.names = FALSE)
  message("Wrote data/", cond, ".csv  (", nrow(out), " rows x ", ncol(out), " cols)")
}
