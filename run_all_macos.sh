#!/bin/bash
# run_all_macos.sh -- runs the notebooks on macOS or Linux.
#
#   ./run_all_macos.sh setup        register the Jupyter kernel (install the R packages
#                                   first: Rscript R/00_install.R)
#   ./run_all_macos.sh datasets     notebook 1: GEO series -> per-dataset tables, Table 1
#   ./run_all_macos.sh analysis     notebook 2: signature, enrichment, figures
#   ./run_all_macos.sh validation   notebook 3: external validation on GSE267650
#   ./run_all_macos.sh all          notebooks 1, 2 and 3, in this order
#
# Jupyter is taken from $JUPYTER if set, otherwise from the PATH. The locale is fixed to
# en_US.UTF-8 because the sort order of gene symbols depends on it.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export RENAL_IRI_ROOT="$ROOT"
export R_LIBS_USER="$ROOT/env/rlib"
export JUPYTER_DATA_DIR="$ROOT/env/jupyter"
export JUPYTER_PATH="$ROOT/env/py/share/jupyter"
export JUPYTER_CONFIG_DIR="$ROOT/env/jupyter_config"   # independent of any personal Jupyter configuration
mkdir -p "$JUPYTER_CONFIG_DIR"
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

JUPYTER="${JUPYTER:-$(command -v jupyter || true)}"
[ -n "$JUPYTER" ] || { echo "jupyter not found: set JUPYTER=/path/to/jupyter" >&2; exit 2; }
export PATH="$(dirname "$JUPYTER"):$PATH"               # R/01_register_kernel.R calls jupyter

run_nb () {
  echo "== notebook/$1 =="
  cd "$ROOT"
  "$JUPYTER" nbconvert --to notebook --execute --inplace \
      --ExecutePreprocessor.timeout=36000 \
      --ExecutePreprocessor.kernel_name=ir-renaliri \
      "notebook/$1"
}

case "${1:-}" in
  setup)      Rscript "$ROOT/R/01_register_kernel.R" ;;
  datasets)   run_nb 01_datasets.ipynb ;;
  analysis)   run_nb 02_analysis.ipynb ;;
  validation) run_nb 03_validation.ipynb ;;
  all)        run_nb 01_datasets.ipynb; run_nb 02_analysis.ipynb; run_nb 03_validation.ipynb ;;
  *)          echo "usage: $0 {setup|datasets|analysis|validation|all}" >&2; exit 2 ;;
esac
