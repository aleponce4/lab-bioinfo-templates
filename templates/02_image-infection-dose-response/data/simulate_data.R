# simulate_data.R
# Generates fully synthetic, high-content-screening-style image-based infection data
# for template 02. Column names follow the generic per-cell export schema produced by
# most HCS platforms (Area/Volume, Sphericity, DAPI, Virus CH2). All values are
# invented by this script from a fixed seed; no experimental data is involved.
# Outputs: data/nucleus.csv  (DAPI channel)   data/virus.csv   (virus channel)
# Run from template folder: Rscript data/simulate_data.R
#
# The cell-level features are deliberately given the structure the template's QC
# panels claim to gate on:
#   * Sphericity is bimodal — round single nuclei near 0.90, plus a doublet /
#     debris population near 0.55 that is also roughly twice the area. This is
#     what the "Doublet & Shape Discrimination" panel separates.
#   * Integrated DNA (Area x DAPI) is bimodal within the singlets — a G0/G1 peak
#     and a G2/M peak at ~2N, with S-phase cells in between and a small
#     sub-G1 (apoptotic) tail. This is what the "Cell Cycle" panel shows.
#   * Toxicity follows a real 4PL with a CC50 inside the tested dose range, so
#     the template's CC50 fit and selectivity index are actually estimable.

set.seed(123)

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

# ── Ground-truth pharmacology ────────────────────────────────────────────────
ic50      <- 8                          # true antiviral EC50 (µM)
hill      <- 0.8                        # antiviral hill slope
cc50      <- 40                         # true cytotoxic CC50 (µM)
hill_tox  <- 2.5                        # cytotoxicity hill slope
# => true selectivity index SI = CC50 / EC50 = 5. CC50 sits inside the 0.01-50
#    µM dose range, so viability crosses 50% and the 4PL fit is identifiable.

mock_rows <- c("C3","C4","C5","C6")     # 4 mock wells
ab_rows   <- c("C11")                   # 1 mock+antibody well (for gate)
tox_rows  <- paste0("D", 1:12)          # toxicity control wells (12 doses)

# Fraction of objects that are doublets (two touching nuclei) and clumps
# (three or more merged nuclei). Doublets are separated by the Sphericity gate;
# clumps are what the AREA_MAX ceiling is for.
doublet_frac <- 0.09
clump_frac   <- 0.02
# Cell-cycle composition among true singlets.
g1_frac  <- 0.60
s_frac   <- 0.22
g2m_frac <- 0.15
subg1_frac <- 1 - (g1_frac + s_frac + g2m_frac)   # apoptotic / pyknotic tail

viability_frac <- function(dose) {
  # 4PL: 1 at dose 0, falling to ~0 at high dose, half-maximal at cc50.
  1 / (1 + (dose / cc50)^hill_tox)
}

all_wells <- list()

# ── Helper: vectorised per-class draw (base R only, no package dependencies) ──
# by_class(labels, G1 = <vector>, S = <vector>, ...) picks, for each element of
# `labels`, the value from the vector named after that label.
by_class <- function(labels, ...) {
  opts <- list(...)
  stopifnot(all(labels %in% names(opts)))
  out <- numeric(length(labels))
  for (nm in names(opts)) {
    idx <- labels == nm
    out[idx] <- opts[[nm]][idx]
  }
  out
}

# ── Helper: generate synthetic cell-level data for one well ──────────────────
# Area, Sphericity and DAPI are drawn jointly per cell-class so that the QC
# panels show genuine structure instead of three independent clouds.
make_cells <- function(well, treatment, conc_uM, n_cells,
                       infected_pct, area_g1 = 230, dapi_mean = 1100) {
  n_cells <- max(1L, as.integer(n_cells))

  # Assign each object to a class.
  singlet_frac <- 1 - doublet_frac - clump_frac
  class_lab <- sample(
    c("doublet", "clump", "subG1", "G1", "S", "G2M"),
    size    = n_cells,
    replace = TRUE,
    prob    = c(doublet_frac,
                clump_frac,
                singlet_frac * subg1_frac,
                singlet_frac * g1_frac,
                singlet_frac * s_frac,
                singlet_frac * g2m_frac)
  )

  # DNA content relative to G1 (drives both nuclear area and DAPI intensity).
  dna_rel <- by_class(
    class_lab,
    doublet = runif(n_cells, 1.85, 2.25),   # two touching nuclei
    clump   = runif(n_cells, 3.20, 5.50),   # three or more merged nuclei
    subG1   = runif(n_cells, 0.35, 0.75),   # apoptotic / pyknotic
    G1      = rnorm(n_cells, 1.00, 0.055),
    S       = runif(n_cells, 1.15, 1.80),
    G2M     = rnorm(n_cells, 2.00, 0.075)
  )

  # Nuclear area scales with DNA content (sub-linearly), plus segmentation noise.
  area <- pmax(15, area_g1 * dna_rel^0.85 * rlnorm(n_cells, 0, 0.11))

  # Sphericity: round singlets high; doublets, clumps and debris markedly lower.
  # This is the bimodality the doublet-discrimination panel gates on.
  sph_mu <- by_class(class_lab,
                     doublet = rep(0.55,  n_cells),
                     clump   = rep(0.44,  n_cells),
                     subG1   = rep(0.78,  n_cells),
                     G1      = rep(0.90,  n_cells),
                     S       = rep(0.90,  n_cells),
                     G2M     = rep(0.90,  n_cells))
  sph_sd <- by_class(class_lab,
                     doublet = rep(0.075, n_cells),
                     clump   = rep(0.090, n_cells),
                     subG1   = rep(0.045, n_cells),
                     G1      = rep(0.045, n_cells),
                     S       = rep(0.045, n_cells),
                     G2M     = rep(0.045, n_cells))
  sphericity <- pmax(0.05, pmin(0.99, rnorm(n_cells, sph_mu, sph_sd)))

  # Mean DAPI intensity: chromatin is more condensed in mitotic and pyknotic
  # nuclei, so mean intensity is NOT flat across the cycle. Integrated intensity
  # (Area x DAPI, computed by the template) then tracks DNA content.
  dapi_scale <- ifelse(class_lab == "subG1", 1.35,
                ifelse(class_lab == "G2M",   1.15, 1.00))
  dapi_scale <- as.numeric(dapi_scale)
  dapi <- pmax(200, dapi_mean * dna_rel^0.20 * dapi_scale * rlnorm(n_cells, 0, 0.16))

  # Virus reporter channel: an infected subpopulation on top of background.
  n_inf     <- round(n_cells * infected_pct)
  virus_inf <- c(rnorm(n_inf, 4800, 1400), rnorm(n_cells - n_inf, 320, 120))
  virus_inf <- pmax(40, virus_inf)[seq_len(n_cells)]

  data.frame(
    WellName     = well,
    FieldIndex   = 1,
    ObjectNumber = seq_len(n_cells),
    dapi_raw     = dapi,
    area_raw     = area,
    sphericity   = sphericity,
    virus_raw    = virus_inf,
    treatment    = treatment,
    conc_uM      = conc_uM
  )
}

# ── Mock wells ────────────────────────────────────────────────────────────────
for (w in mock_rows) {
  n_c <- round(rnorm(1, 400, 25))
  all_wells[[w]] <- make_cells(w, "Mock", NA, n_c, 0.00)
}
all_wells[["C11"]] <- make_cells("C11", "Mock+Ab", NA, 400, 0.00)

# ── Virus only (dose = 0): wells C7, C8 — matches PLATE_MAP ─────────────────
for (r in 1:2) {
  w <- paste0("C", 7 + r - 1)
  n_c <- round(rnorm(1, 400, 30))
  vo_pct <- max(0.55, min(0.75, 0.65 + rnorm(1, 0, 0.045)))
  all_wells[[w]] <- make_cells(w, "Virus Only", 0, n_c, vo_pct)
}

# ── Virus + Compound dose series ──────────────────────────────────────────────
# 12 concentrations × 2 replicates = 24 wells (rows E-F, cols 1-12)
# Adds realistic well-to-well biological replicate noise (sd = 0.045) and background floor (~0.02)
doses <- c(0.01, 0.1, 0.5, 1, 2, 4, 8, 12, 16, 20, 30, 50)
vc_layout <- data.frame(
  conc = rep(doses, 2),
  well = c(paste0("E", 1:12), paste0("F", 1:12)),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(vc_layout))) {
  conc <- vc_layout$conc[i]
  w    <- vc_layout$well[i]
  expected_pct <- 0.65 / (1 + (conc / ic50)^hill)
  # Well-to-well replicate variance + background floor
  well_pct <- max(0.02, min(0.72, expected_pct + rnorm(1, 0, 0.045)))
  # The compound is cytotoxic in the infected arm too, on the same CC50 curve.
  viab <- max(0.02, viability_frac(conc) + rnorm(1, 0, 0.03))
  n_c  <- round(rnorm(1, 400, 35) * viab)
  all_wells[[w]] <- make_cells(w, "Virus + Compound", conc, n_c, well_pct)
}

# ── Toxicity controls (drug only, no virus) — the CC50 arm ───────────────────
tox_doses <- doses   # match virus+compound dose series
for (i in seq_along(tox_rows)) {
  viab   <- max(0.02, viability_frac(tox_doses[i]) + rnorm(1, 0, 0.03))
  n_surv <- round(400 * viab)
  all_wells[[tox_rows[i]]] <- make_cells(tox_rows[i], "Toxicity Control",
                                          tox_doses[i], n_surv, 0.00)
}

df <- do.call(rbind, all_wells)

# ── Sanity checks on the generated ground truth ──────────────────────────────
tox_counts <- vapply(tox_rows, function(w) nrow(all_wells[[w]]), integer(1))
mock_count <- mean(vapply(mock_rows, function(w) nrow(all_wells[[w]]), integer(1)))
tox_viab   <- tox_counts / mock_count
stopifnot(
  # Viability must actually cross 50% inside the dose range, or no CC50 exists.
  max(tox_viab) > 0.5, min(tox_viab) < 0.5,
  # Sphericity must be bimodal enough for a gate to separate the two modes.
  diff(range(quantile(df$sphericity, c(0.02, 0.98)))) > 0.25,
  # Clumps must exist above the AREA_MAX ceiling, or that knob is a no-op.
  mean(df$area_raw > 700) > 0.005,
  all(df$area_raw > 0), all(df$dapi_raw > 0), all(df$virus_raw > 0)
)

# Write split nucleus / virus files (generic two-channel per-cell HCS export layout)
nucleus <- data.frame(
  WellName                       = df$WellName,
  FieldIndex                     = df$FieldIndex,
  ObjectNumber                   = df$ObjectNumber,
  `(nucleus) MeanIntensity CH1`   = df$dapi_raw,
  `(nucleus) Area`                = df$area_raw,
  `(nucleus) Sphericity`          = df$sphericity,
  check.names = FALSE
)
virus <- data.frame(
  WellName                       = df$WellName,
  FieldIndex                     = df$FieldIndex,
  ObjectNumber                   = df$ObjectNumber,
  `(Virus) MeanIntensity CH2`    = df$virus_raw,
  check.names = FALSE
)

write.csv(nucleus, "data/nucleus.csv",  row.names = FALSE)
write.csv(virus,   "data/virus.csv",    row.names = FALSE)
message("Wrote data/nucleus.csv and data/virus.csv  (", nrow(df), " cells)")
message("Ground truth: EC50 = ", ic50, " uM, CC50 = ", cc50,
        " uM, SI = ", round(cc50 / ic50, 2))
