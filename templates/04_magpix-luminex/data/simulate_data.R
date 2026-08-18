# simulate_data.R
# Generates synthetic MagPix/Luminex multiplex cytokine data for template 04.
# Output: data/luminex_data.csv
# Run from template folder: Rscript data/simulate_data.R

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

analytes  <- c("CCL2", "CCL5", "IL-6", "IL-10", "TNF-alpha",
               "IFN-gamma", "CXCL10", "IL-1beta")
factor1   <- c("PBS", "Disease")    # two levels of Factor 1
factor2   <- c("Vehicle", "Drug")   # two levels of Factor 2
n_per_group <- 5

# ── Assay dynamic range, per analyte ─────────────────────────────────────────
# Real Luminex kits have a separate standard curve per analyte, so the limit of
# detection (LoD) and the upper limit of quantification (ULoQ) differ per bead
# region. The values below are chosen so the demo data actually EXERCISES both
# kinds of censoring that the template handles:
#   - low-abundance analytes (TNF-alpha, IFN-gamma, IL-1beta) fall below the LoD
#     in the unstimulated PBS groups  -> "< LoD"  -> LoD/sqrt(2) substitution
#   - high-abundance analytes (CCL2, IL-6, CXCL10) saturate the top standard in
#     the Disease groups                -> "> ULoQ" -> ULoQ substitution
# Without this, the censoring-rate QC table renders as all zeros and the
# template's headline LoD/sqrt(2) step is never executed.
lod_map <- c(
  "CCL2" = 2.0, "CCL5" = 2.0, "IL-6" = 2.0, "IL-10" = 2.0,
  "TNF-alpha" = 10.0, "IFN-gamma" = 8.0, "CXCL10" = 2.0, "IL-1beta" = 12.0
)
uloq_map <- c(
  "CCL2" = 2500, "CCL5" = 10000, "IL-6" = 3500, "IL-10" = 10000,
  "TNF-alpha" = 10000, "IFN-gamma" = 10000, "CXCL10" = 7500, "IL-1beta" = 10000
)

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

rows        <- list()
n_below_lod <- 0L
n_above_ul  <- 0L

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
    if (true_val < lod_map[[an]]) {
      row[[an]] <- paste0("< ", lod_map[[an]])   # left-censored
      n_below_lod <- n_below_lod + 1L
    } else if (true_val > uloq_map[[an]]) {
      row[[an]] <- paste0("> ", uloq_map[[an]])  # right-censored
      n_above_ul <- n_above_ul + 1L
    } else {
      row[[an]] <- round(true_val, 2)
    }
  }
  rows[[i]] <- row
}

df <- do.call(rbind, rows)

# ── Non-biological wells ─────────────────────────────────────────────────────
# Every real Luminex export carries standard-curve, blank and background rows.
# They are included here so the template's non-sample filter has something to
# remove and the reported drop count is non-zero.
qc_wells <- data.frame(
  Location    = paste0("Row", nrow(df) %/% 4 + seq_len(4)),
  Sample      = c("Standard1", "Standard7", "Blank", "Background0"),
  Original_ID = c("STD1", "STD7", "BLK", "BKG"),
  `Total Events` = sample(80:120, 4, replace = TRUE),
  check.names = FALSE
)
for (an in analytes) {
  qc_wells[[an]] <- c(
    round(uloq_map[[an]] * 0.9, 2),        # top standard
    round(lod_map[[an]] * 1.2, 2),         # bottom standard
    paste0("< ", lod_map[[an]]),           # blank
    paste0("< ", lod_map[[an]])            # background
  )
}
df <- rbind(df, qc_wells)

write.csv(df, "data/luminex_data.csv", row.names = FALSE)
message("Wrote data/luminex_data.csv  (", nrow(df), " rows incl. ",
        nrow(qc_wells), " non-sample QC wells, ", length(analytes), " analytes)")
message("  left-censored  (< LoD):  ", n_below_lod, " sample measurements")
message("  right-censored (> ULoQ): ", n_above_ul, " sample measurements")
stopifnot(n_below_lod > 0, n_above_ul > 0)
