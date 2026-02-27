# simulate_data.R
# Generates synthetic CQ1-style image-based infection data for template 02.
# Outputs: data/nucleus.csv  (DAPI channel)   data/virus.csv   (virus channel)
# Run from template folder: Rscript data/simulate_data.R

set.seed(123)

# ── Parameters ────────────────────────────────────────────────────────────────
doses_uM  <- c(0, 4, 8, 12, 16, 20)   # treatment concentrations (0 = virus only)
n_wells   <- 2                          # replicate wells per dose
n_cells   <- 400                        # cells per well

ic50      <- 8                          # true IC50 (µM)
hill      <- 0.8                        # hill slope
mock_rows <- c("C3","C4","C5","C6")     # 4 mock wells
ab_rows   <- c("C11")                   # 1 mock+antibody well (for gate)
tox_rows  <- paste0("D", 3:8)           # toxicity control wells

all_wells <- list()

# ── Helper: generate cell-level data for one well ────────────────────────────
make_cells <- function(well, treatment, conc_uM, n_cells,
                       infected_pct, area_mean = 150, dapi_mean = 3000) {
  area   <- pmax(40, rnorm(n_cells, area_mean, 30))
  dapi   <- pmax(500, rnorm(n_cells, dapi_mean, 600))
  n_inf  <- round(n_cells * infected_pct)
  virus_inf <- c(rnorm(n_inf, 4000, 800), rnorm(n_cells - n_inf, 250, 80))
  virus_inf <- pmax(50, virus_inf)
  data.frame(
    WellName     = well,
    FieldIndex   = 1,
    ObjectNumber = seq_len(n_cells),
    dapi_raw     = dapi,
    area_raw     = area,
    virus_raw    = virus_inf,
    treatment    = treatment,
    conc_uM      = conc_uM
  )
}

# ── Mock wells ────────────────────────────────────────────────────────────────
for (w in mock_rows) {
  all_wells[[w]] <- make_cells(w, "Mock", NA, n_cells, 0.00)
}
all_wells[["C11"]] <- make_cells("C11", "Mock+Ab", NA, n_cells, 0.00)

# ── Virus only + dose series ──────────────────────────────────────────────────
for (d in seq_along(doses_uM)) {
  conc <- doses_uM[d]
  # 4PL: infection fraction = 1 / (1 + (conc/IC50)^hill)
  inf_pct <- ifelse(conc == 0, 0.65,
                    0.65 / (1 + (conc / ic50)^hill))
  for (r in seq_len(n_wells)) {
    w <- paste0("E", r + 2, "_", conc, "uM")
    if (conc == 0) w <- paste0("C", 7 + r - 1)
    all_wells[[w]] <- make_cells(w, ifelse(conc == 0, "Virus Only", "Virus + Compound"),
                                 conc, n_cells, inf_pct)
  }
}

# ── Toxicity controls ─────────────────────────────────────────────────────────
tox_doses <- c(20, 16, 12, 8, 4, 1)
for (i in seq_along(tox_rows)) {
  viab <- max(0.60, 1 - 0.01 * tox_doses[i])   # mild toxicity
  n_surv <- round(n_cells * viab)
  all_wells[[tox_rows[i]]] <- make_cells(tox_rows[i], "Toxicity Control",
                                          tox_doses[i], n_surv, 0.00)
}

df <- do.call(rbind, all_wells)

# Write split nucleus / virus files (mimicking CQ1 export format)
nucleus <- data.frame(
  WellName     = df$WellName,
  FieldIndex   = df$FieldIndex,
  ObjectNumber = df$ObjectNumber,
  `(nucleus) MeanIntensity CH1` = df$dapi_raw,
  `(nucleus) Area`              = df$area_raw,
  check.names = FALSE
)
virus <- data.frame(
  WellName     = df$WellName,
  FieldIndex   = df$FieldIndex,
  ObjectNumber = df$ObjectNumber,
  `(Virus) MeanIntensity CH2`  = df$virus_raw,
  check.names = FALSE
)

write.csv(nucleus, "data/nucleus.csv",  row.names = FALSE)
write.csv(virus,   "data/virus.csv",    row.names = FALSE)
message("Wrote data/nucleus.csv and data/virus.csv  (", nrow(df), " cells)")
