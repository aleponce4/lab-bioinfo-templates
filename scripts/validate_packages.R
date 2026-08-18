# validate_packages.R
# Scans every R file in the repo for package references and checks two things:
#   1. every referenced package is DECLARED in install_packages.R, and
#   2. every referenced package is actually INSTALLED.
#
# (1) is the check that matters. The previous version only did (2), which passes
# whenever a package happens to arrive as somebody else's transitive dependency — so
# packages could be used for months without ever being declared.
#
# Usage: Rscript scripts/validate_packages.R

if (!dir.exists("templates")) {
  stop("Run this from the repository root.", call. = FALSE)
}

# ── Which files to scan ───────────────────────────────────────────────────────
# .qmd is included so index.qmd and the gallery wrappers are covered too; Python
# templates simply yield no matches.
all_files <- c(
  list.files("templates", pattern = "\\.(R|Rmd|qmd)$", full.names = TRUE, recursive = TRUE),
  list.files("scripts", pattern = "\\.R$", full.names = TRUE),
  "index.qmd"
)
all_files <- unique(all_files[file.exists(all_files)])

# ── Extract package references ────────────────────────────────────────────────
# Matches library(x), require(x), requireNamespace("x") and x::fn — the last is how
# most of the Bioconductor usage in these templates is actually written.
PATTERNS <- c(
  "(?:library|require)\\(\\s*[\"']?([A-Za-z][A-Za-z0-9._]*)[\"']?\\s*[,)]",
  "requireNamespace\\(\\s*[\"']([A-Za-z][A-Za-z0-9._]*)[\"']",
  "([A-Za-z][A-Za-z0-9._]*):::?[A-Za-z._]"
)

pkgs <- character(0)
for (f in all_files) {
  lines <- readLines(f, warn = FALSE)
  # Drop comments so commented-out code and prose don't register.
  lines <- sub("#.*$", "", lines)
  for (pattern in PATTERNS) {
    m <- regmatches(lines, regexec(pattern, lines, perl = TRUE))
    hits <- vapply(m, function(x) if (length(x) >= 2) x[2] else NA_character_, character(1))
    pkgs <- c(pkgs, hits[!is.na(hits)])
  }
}

# Base and recommended packages ship with R and are never declared.
base_pkgs <- rownames(installed.packages(priority = c("base", "recommended")))
pkgs <- sort(unique(setdiff(pkgs, c(base_pkgs, "base", "stats", "utils", "grid", "tools",
                                    "methods", "graphics", "grDevices", "parallel"))))

# ── Compare against install_packages.R ────────────────────────────────────────
# Optional at runtime: every use site is behind a requireNamespace() guard and the
# template degrades gracefully without them. Not required to be declared or installed.
OPTIONAL <- c("gurobi", "Rcplex", "lpSolve", "CARNIVAL")

# Pulled in as a hard dependency of a package that IS declared.
TRANSITIVE <- c("S4Vectors", "IRanges", "GenomicRanges", "Biobase", "BiocGenerics")

# Members of tidyverse, added to `declared` below. Tracked separately so they are not
# then reported as "declared but unused" — they were never listed by hand.
TIDYVERSE_MEMBERS <- c("dplyr", "tidyr", "tibble", "purrr", "readr",
                       "stringr", "forcats", "ggplot2", "lubridate")

declared <- character(0)
if (file.exists("install_packages.R")) {
  src <- paste(readLines("install_packages.R", warn = FALSE), collapse = "\n")
  # Read only the cran_pkgs <- c(...) / bioc_pkgs <- c(...) vectors, so unrelated
  # strings elsewhere in the installer aren't mistaken for package names.
  for (v in c("cran_pkgs", "bioc_pkgs")) {
    block <- regmatches(src, regexpr(paste0(v, "\\s*<-\\s*c\\((?:[^()]|\\([^()]*\\))*\\)"),
                                     src, perl = TRUE))
    if (length(block)) {
      quoted <- regmatches(block, gregexpr("\"[A-Za-z][A-Za-z0-9._]*\"", block))
      declared <- c(declared, gsub("\"", "", unlist(quoted)))
    }
  }
}
# tidyverse pulls in its own members; treat them as declared.
if ("tidyverse" %in% declared) {
  declared <- c(declared, TIDYVERSE_MEMBERS)
}
declared <- unique(declared)

undeclared <- setdiff(pkgs, c(declared, OPTIONAL, TRANSITIVE))
inst <- rownames(installed.packages())
missing <- setdiff(pkgs, c(inst, OPTIONAL))

status <- 0

if (length(undeclared)) {
  cat("UNDECLARED — used in a template but not listed in install_packages.R:\n")
  cat(paste0("  - ", undeclared, collapse = "\n"), "\n\n", sep = "")
  status <- 1
}

if (length(missing)) {
  cat("NOT INSTALLED — referenced by a template but absent from this library:\n")
  cat(paste0("  - ", missing, collapse = "\n"), "\n", sep = "")
  cat("\nIf these ARE in install_packages.R, their installation failed silently;\n")
  cat("check the install step log for missing system libraries.\n\n")
  status <- 1
}

unused <- setdiff(declared, c(pkgs, OPTIONAL, TIDYVERSE_MEMBERS, "BiocManager", "tidyverse"))
if (length(unused)) {
  cat("note: declared in install_packages.R but never referenced:\n")
  cat(paste0("  - ", unused, collapse = "\n"), "\n\n", sep = "")
}

if (status == 0) {
  cat(sprintf("OK — all %d referenced packages are declared and installed.\n", length(pkgs)))
}

quit(status = status)
