#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

: "${DRY_RUN:=0}"

parse_grid() {
  local grid="$1"
  if [[ "${grid}" =~ ^([0-9]+)$ ]]; then
    GRID_NX="${BASH_REMATCH[1]}"
    GRID_NY="${BASH_REMATCH[1]}"
  elif [[ "${grid}" =~ ^([0-9]+)[xX]([0-9]+)$ ]]; then
    GRID_NX="${BASH_REMATCH[1]}"
    GRID_NY="${BASH_REMATCH[2]}"
  else
    echo "Error: invalid grid '${grid}'. Expected N or NXxNY, e.g. 256 or 256x256" >&2
    exit 1
  fi
}

# (Opcional) Ajustes gerais de execução:
export IO=1
export RUN=1
export CLEAN=0
export DEBUG=0
export RDC=0
export PROGRESS=0
export PROGRESS_HZ=2
export WARMUP=100
export REAL=float   # ou double se quiser

# Tamanhos de grid. Use "256" para 256x256, ou "512x256" para retangular.
grids=(
  "32"
  "64"
  "128"
)

cases=(
  "D2Q9  3200"
)

for grid in "${grids[@]}"; do
  parse_grid "${grid}"

  for c in "${cases[@]}"; do
    set -- $c
    stencil="$1"
    re="$2"

    ts="$(date +%Y%m%d_%H%M%S)"
    run_id="${ts}_${stencil}_${GRID_NX}x${GRID_NY}_RE${re}"

    echo "================================================="
    echo "[CASE] STENCIL=${stencil}  RE=${re}  GRID=${GRID_NX}x${GRID_NY}  RUN_ID=${run_id}"
    echo "================================================="

    cmd=(bash "${ROOT_DIR}/compile.sh" \
      --stencil "${stencil}" \
      --re "${re}" \
      --grid "${GRID_NX}x${GRID_NY}" \
      --run_id "${run_id}")

    if [[ "${DRY_RUN}" == "1" ]]; then
      printf '[DRY_RUN]'
      printf ' %q' "${cmd[@]}"
      printf '\n'
    else
      "${cmd[@]}"
    fi

    echo
  done
done

echo "ALL CASES DONE."
