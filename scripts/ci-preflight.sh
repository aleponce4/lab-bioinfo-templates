#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO_ROOT"

readonly PYTHON_TEMPLATE_DIR="templates/08_vcf-mutation-analysis"
KEEP_GOING=0
RENDER_ONLY=0
SELECTED_PAGE=""

PAGES=(
  "index.qmd"
  "templates/01_plaque-assay-violin/template.qmd"
  "templates/02_image-infection-dose-response/template.qmd"
  "templates/03_multicondition-anova/template.qmd"
  "templates/04_magpix-luminex/template.qmd"
  "templates/05_phylo-geographic/template.qmd"
  "templates/06_go-enrichment/template.qmd"
  "templates/07_rnaseq-deseq2/template.qmd"
  "templates/08_vcf-mutation-analysis/template.qmd"
)

EXPECTED_OUTPUTS=(
  "docs/index.html"
  "docs/templates/01_plaque-assay-violin/template.html"
  "docs/templates/02_image-infection-dose-response/template.html"
  "docs/templates/03_multicondition-anova/template.html"
  "docs/templates/04_magpix-luminex/template.html"
  "docs/templates/05_phylo-geographic/template.html"
  "docs/templates/06_go-enrichment/template.html"
  "docs/templates/07_rnaseq-deseq2/template.html"
  "docs/templates/08_vcf-mutation-analysis/template.html"
)

usage() {
  cat <<'EOF'
Usage: bash scripts/ci-preflight.sh [--render-only] [--keep-going] [--page <path>]

Runs the same page-by-page Quarto render flow used in CI.
Without arguments it also regenerates synthetic example data first.
Use --keep-going to continue past page render failures and report them all at the end.
Use --page to render just one page from the configured site page list.
EOF
}

require_cmd() {
  local cmd=$1

  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Missing required command: $cmd" >&2
    exit 1
  fi
}

# ── Static preflight: catches obvious regressions without running anything ────
static_preflight() {
  local errors=0

  echo "==> Running static preflight checks"

  # 1. No runtime package installs in simulate_data scripts or template Rmd files.
  #    install_packages.R is the single source of truth for dependencies.
  for f in templates/*/data/simulate_data.R templates/*/template.Rmd; do
    [[ -f "$f" ]] || continue
    if grep -qE "BiocManager::install|install\.packages" "$f"; then
      echo "  ERROR: runtime package install found in $f" >&2
      errors=$((errors + 1))
    fi
  done

  # 2. Every simulate_data.R must produce the data file(s) its template.Rmd reads.
  #    We detect the DATA_FILE / DDS_FILE assignment and check the path exists
  #    relative to the template directory (after data generation has run).
  #    Skip this sub-check when --render-only is set (data not yet generated).
  if [[ $RENDER_ONLY -eq 0 ]]; then
    for rmd in templates/*/template.Rmd; do
      [[ -f "$rmd" ]] || continue
      tmpl_dir=$(dirname "$rmd")
      while IFS= read -r assignment; do
        data_path=$(echo "$assignment" | sed 's/.*<- *"\(.*\)".*/\1/')
        full_path="$tmpl_dir/$data_path"
        if [[ ! -f "$full_path" ]]; then
          echo "  ERROR: $rmd declares $data_path but $full_path does not exist" >&2
          errors=$((errors + 1))
        fi
      done < <(grep -E '^\s*(DATA_FILE|DDS_FILE)\s*<-\s*"' "$rmd" || true)
    done
  fi

  if [[ $errors -gt 0 ]]; then
    echo "Static preflight failed with $errors error(s). Fix before rendering." >&2
    return 1
  fi

  echo "==> Static preflight passed"
}

prepare_output_dirs() {
  rm -rf docs

  for output in $(expected_outputs); do
    mkdir -p "$(dirname "$output")"
  done
}

selected_pages() {
  if [[ -n "$SELECTED_PAGE" ]]; then
    printf '%s\n' "$SELECTED_PAGE"
    return
  fi

  printf '%s\n' "${PAGES[@]}"
}

expected_output_for_page() {
  case "$1" in
    index.qmd)
      printf '%s\n' "docs/index.html"
      ;;
    templates/*/template.qmd)
      printf '%s\n' "docs/${1%.qmd}.html"
      ;;
    *)
      echo "Unknown page: $1" >&2
      exit 2
      ;;
  esac
}

expected_outputs() {
  local page

  while IFS= read -r page; do
    expected_output_for_page "$page"
  done < <(selected_pages)
}

generate_r_data() {
  local script template_dir

  shopt -s nullglob

  for script in templates/*/data/simulate_data.R; do
    template_dir=$(dirname "$(dirname "$script")")
    echo "==> Generating R data: $script (cwd: $template_dir)"
    (cd "$template_dir" && Rscript data/simulate_data.R)
  done

  shopt -u nullglob
}

generate_python_data() {
  if ! command -v conda >/dev/null 2>&1; then
    echo "Warning: conda not found; skipping Python synthetic data generation." >&2
    return
  fi

  if ! conda run -n jonsson-bioinfo python -c "import sys" >/dev/null 2>&1; then
    echo "Warning: conda environment 'jonsson-bioinfo' is unavailable; skipping Python synthetic data generation." >&2
    return
  fi

  echo "==> Generating Python data: $PYTHON_TEMPLATE_DIR/data/simulate_data.py"
  (
    cd "$PYTHON_TEMPLATE_DIR"
    conda run -n jonsson-bioinfo python data/simulate_data.py
  )
}

render_pages() {
  local failures=()
  local page

  require_cmd quarto
  require_cmd Rscript

  prepare_output_dirs

  while IFS= read -r page; do
    echo "==> Rendering $page"

    if quarto render "$page"; then
      continue
    fi

    if [[ $KEEP_GOING -eq 1 ]]; then
      failures+=("$page")
      echo "Render failed for $page; continuing because --keep-going is enabled." >&2
      continue
    fi

    return 1
  done < <(selected_pages)

  if [[ ${#failures[@]} -gt 0 ]]; then
    echo "Render failures:" >&2
    printf ' - %s\n' "${failures[@]}" >&2
    return 1
  fi
}

validate_manifest() {
  local missing=0
  local output

  while IFS= read -r output; do
    if [[ ! -f "$output" ]]; then
      echo "Missing expected output: $output" >&2
      missing=1
    fi
  done < <(expected_outputs)

  if [[ $missing -ne 0 ]]; then
    exit 1
  fi

  echo "Render manifest verified."
}

main() {
  while [[ $# -gt 0 ]]; do
    case "$1" in
      --render-only)
        RENDER_ONLY=1
        shift
        ;;
      --keep-going)
        KEEP_GOING=1
        shift
        ;;
      --page)
        [[ $# -ge 2 ]] || { echo "--page requires a value" >&2; exit 2; }
        SELECTED_PAGE=$2
        shift 2
        ;;
      --help|-h)
        usage
        exit 0
        ;;
      *)
        usage >&2
        exit 2
        ;;
    esac
  done

  if [[ -n "$SELECTED_PAGE" ]] && [[ ! " ${PAGES[*]} " =~ [[:space:]]${SELECTED_PAGE//\//\\/}[[:space:]] ]]; then
    echo "--page must be one of the configured site pages." >&2
    exit 2
  fi

  static_preflight || exit 1

  if [[ $RENDER_ONLY -eq 0 ]]; then
    require_cmd Rscript
    generate_r_data
    generate_python_data
  fi

  if command -v quarto >/dev/null 2>&1; then
    render_pages
    validate_manifest
  else
    echo "Quarto not found; skipping render and manifest validation." >&2
  fi
}

main "$@"
