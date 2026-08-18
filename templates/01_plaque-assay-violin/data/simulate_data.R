# simulate_data.R
# Generates synthetic plaque assay data for template 01.
# Output: data/titer_data.csv
# Run from template folder: Rscript data/simulate_data.R
#
# IMPORTANT — single source of truth for the group semantics
# --------------------------------------------------------
# GROUP_DEFS below is the ONE place where "Group N" is tied to a treatment
# label and to that treatment's expected mean log10 titer. The template's
# GROUP_MAP must use the same "Group N" -> label mapping; the template asserts
# this with stopifnot() after loading, so the two cannot silently drift.

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

# ── Group definitions: label + mean log10 titer, keyed by CSV "Group" value ────
# Mock arms sit near zero, Vehicle (untreated infected control) is the highest,
# single agents reduce titer modestly, and the combination arms are the lowest.
GROUP_DEFS <- data.frame(
  group     = paste0("Group ", 1:8),
  label     = c(
    "Mock + Compound A",
    "Mock + Compound B",
    "Compound A",
    "Compound B",
    "Drug X",
    "Compound A + Drug X",
    "Compound B + Drug X",
    "Vehicle"
  ),
  mean_log  = c(0.5, 0.5, 5.8, 5.7, 5.9, 5.0, 4.8, 6.2),
  stringsAsFactors = FALSE
)
rownames(GROUP_DEFS) <- GROUP_DEFS$group

# Named vector used for the draws; names are the CSV "Group N" values.
GROUP_MEAN_LOG <- setNames(GROUP_DEFS$mean_log, GROUP_DEFS$group)

# Mock arms are identified by their label, not by row position.
is_mock <- setNames(grepl("^Mock", GROUP_DEFS$label), GROUP_DEFS$group)

stopifnot(
  !any(duplicated(GROUP_DEFS$group)),
  !any(duplicated(GROUP_DEFS$label)),
  # Vehicle must be the highest-titer arm; mocks must be the lowest.
  GROUP_MEAN_LOG[["Group 8"]] == max(GROUP_MEAN_LOG),
  all(GROUP_MEAN_LOG[is_mock] < min(GROUP_MEAN_LOG[!is_mock]))
)

sd_log      <- 0.4
n_per_group <- 12

rows <- list()
mouse_id <- 1
for (g in names(GROUP_MEAN_LOG)) {
  for (j in seq_len(n_per_group)) {
    log_t <- rnorm(1, GROUP_MEAN_LOG[[g]], sd_log)
    titer <- max(0, round(10^log_t))
    # generate two replicate plaque counts around the titer
    r1 <- max(0, round(rnorm(1, titer / 100, titer / 300)))
    r2 <- max(0, round(rnorm(1, titer / 100, titer / 300)))
    rows[[length(rows) + 1]] <- data.frame(
      `Mouse ID`      = mouse_id,
      Group           = g,
      Virus           = ifelse(is_mock[[g]], "Mock PBS", "Pathogen X"),
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
