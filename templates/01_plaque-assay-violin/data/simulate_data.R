# simulate_data.R
# Generates synthetic plaque assay data for template 01.
# Output: data/titer_data.csv
# Run from template folder: Rscript data/simulate_data.R

set.seed(42)

groups <- c(
  "Mock + Compound A",
  "Mock + Compound B",
  "Vehicle",
  "Drug X",
  "Compound A",
  "Compound B",
  "Compound A + Drug X",
  "Compound B + Drug X"
)

# Mean log10 titer per group (mock groups near 0, infected groups high)
mean_log <- c(0.5, 0.5, 6.2, 5.9, 5.8, 5.7, 5.0, 4.8)
sd_log   <- 0.4
n_per_group <- 6

rows <- list()
mouse_id <- 1
for (i in seq_along(groups)) {
  for (j in seq(n_per_group)) {
    log_t <- rnorm(1, mean_log[i], sd_log)
    titer <- max(0, round(10^log_t))
    # generate two replicate plaque counts around the titer
    r1 <- max(0, round(rnorm(1, titer / 100, titer / 300)))
    r2 <- max(0, round(rnorm(1, titer / 100, titer / 300)))
    rows[[length(rows) + 1]] <- data.frame(
      `Mouse ID`      = mouse_id,
      Group           = paste0("Group ", i),
      Virus           = ifelse(i <= 2, "Mock PBS", "Pathogen X"),
      `# plaques counted R1` = r1,
      `# plaques counted R2` = r2,
      `Average plaques from duplicate wells in a 12-well plate` = (r1 + r2) / 2,
      `Plated 10e-1 to 10e-6 dilution` = "10e-2",
      `Total Volume (ml)` = 1.0,
      `pfu/mL right-side lung homogenate` = titer,
      check.names = FALSE
    )
    mouse_id <- mouse_id + 1
  }
}

df <- do.call(rbind, rows)
write.csv(df, "data/titer_data.csv", row.names = FALSE)
message("Wrote data/titer_data.csv  (", nrow(df), " rows)")
