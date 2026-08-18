#!/usr/bin/env Rscript
# Parses every R chunk in the repo without loading a single package.
#
# This is deliberately dependency-free (base R only) so CI can run it on a bare
# `setup-r` step: a syntax error in a template.Rmd then fails in seconds instead of
# forty minutes into the render job.

# Run from the repo root.
if (!dir.exists("templates")) {
  stop("Run this from the repository root (no templates/ directory here).", call. = FALSE)
}

# Pull the bodies of ```{r ...} chunks out of an .Rmd/.qmd file.
extract_r_chunks <- function(path) {
  lines <- readLines(path, warn = FALSE)
  opens <- grep("^\\s*```+\\s*\\{[rR][ ,}]", lines)
  fences <- grep("^\\s*```+\\s*$", lines)

  unlist(lapply(opens, function(start) {
    end <- fences[fences > start]
    if (!length(end)) return(character(0))
    # An empty chunk body (e.g. ```{r child="template.Rmd"}) has nothing between the
    # fences; seq() would count backwards and pick up the fences themselves.
    if (end[1] <= start + 1L) return(character(0))
    lines[seq(start + 1L, end[1] - 1L)]
  }))
}

targets <- c(
  Sys.glob("templates/*/template.Rmd"),
  Sys.glob("templates/*/template.qmd"),
  Sys.glob("templates/*/data/simulate_data.R"),
  Sys.glob("scripts/*.R"),
  "index.qmd",
  "install_packages.R"
)
targets <- targets[file.exists(targets)]

failures <- character(0)

for (f in targets) {
  code <- if (grepl("\\.(Rmd|qmd)$", f)) extract_r_chunks(f) else readLines(f, warn = FALSE)

  # A .qmd holding Python chunks yields no R code; nothing to check.
  if (!length(code) || !any(nzchar(trimws(code)))) {
    cat(sprintf("  skip  %-55s (no R chunks)\n", f))
    next
  }

  res <- tryCatch({
    parse(text = paste(code, collapse = "\n"))
    NULL
  }, error = function(e) conditionMessage(e))

  if (is.null(res)) {
    cat(sprintf("  ok    %s\n", f))
  } else {
    cat(sprintf("  FAIL  %s\n        %s\n", f, res))
    failures <- c(failures, f)
  }
}

if (length(failures)) {
  stop("R syntax errors in: ", paste(failures, collapse = ", "), call. = FALSE)
}

cat(sprintf("\nAll %d file(s) parsed cleanly.\n", length(targets)))
