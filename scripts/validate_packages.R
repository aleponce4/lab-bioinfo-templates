# validate_packages.R
# Scans all template Rmd files and simulate_data scripts for library() / require()
# calls and verifies every referenced package is installed.
# Exits non-zero if any are missing, so CI fails fast before the render phase.
# Usage: Rscript scripts/validate_packages.R

all_files <- c(
  list.files("templates", pattern = "template\\.Rmd$",
             full.names = TRUE, recursive = TRUE),
  list.files("templates", pattern = "simulate_data\\.R$",
             full.names = TRUE, recursive = TRUE)
)

pkgs <- character(0)
for (f in all_files) {
  lines   <- readLines(f, warn = FALSE)
  matches <- regmatches(
    lines,
    gregexpr("(?:library|require)\\(([A-Za-z0-9._]+)\\)", lines, perl = TRUE)
  )
  for (m in unlist(matches)) {
    pkg <- sub(".*\\(([A-Za-z0-9._]+)\\).*", "\\1", m, perl = TRUE)
    pkgs <- c(pkgs, pkg)
  }
}

pkgs    <- sort(unique(pkgs))
inst    <- rownames(installed.packages())
missing <- pkgs[!pkgs %in% inst]

if (length(missing) > 0) {
  cat("MISSING packages (imported in templates but not installed — check install_packages.R and system deps):\n",
      paste0("  - ", missing, collapse = "\n"), "\n", sep = "")
  cat("\nIf these packages ARE in install_packages.R, their installation failed silently.\n")
  cat("Check the 'Install R packages' step log above for errors.\n")
  quit(status = 1)
}

cat(sprintf(
  "OK — all %d packages imported across templates are installed.\n",
  length(pkgs)
))
