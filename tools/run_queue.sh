#!/usr/bin/env bash
set -euo pipefail

export IO=1
export RUN=1
export CLEAN=0
export DEBUG=0
export RDC=0
export PROGRESS=1
export PROGRESS_HZ=2
export WARMUP=100
export REAL=float

stencils=(
  "D2Q9"
)

grids=(
  32
  64
  128
  256
  512
  1024
  2048
  4096
)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

for stencil in "${stencils[@]}"; do
  for grid in "${grids[@]}"; do
    ts="$(date +%Y%m%d_%H%M%S)"

    export STENCIL="${stencil}"
    export GRID="${grid}"
    export RUN_ID="${ts}_${stencil}_G${GRID}"

    echo "================================================="
    echo "[CASE] STENCIL=${STENCIL}  GRID=${GRID}  RUN_ID=${RUN_ID}"
    echo "================================================="

    bash "${ROOT_DIR}/compile.sh"

    echo
  done
done

echo "ALL CASES DONE."