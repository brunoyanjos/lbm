#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

: "${DRY_RUN:=0}"
: "${T_FINAL:=2000}"

show_help() {
  cat <<EOF
Usage:
  bash tools/run_queue.sh [options]

Options:
  --stencil VALUE[,VALUE...]  Stencil(s) to run (default: D2Q9)
  --grid VALUE[,VALUE...]     Grid size(s), N or NXxNY (default: 512)
  --re VALUE[,VALUE...]       Reynolds number(s) (default: 100,400,1000,3200,5000,7500,10000)
  --t_final VALUE             Final t* passed to compile.sh --t_star_end (default: ${T_FINAL})
  --dry-run                   Print commands without running them
  -h, --help                  Show this help
EOF
}

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

append_csv_values() {
  local target_name="$1"
  local csv="$2"
  local value
  IFS=',' read -ra values <<< "${csv}"
  for value in "${values[@]}"; do
    [[ -n "${value}" ]] || continue
    eval "${target_name}+=(\"\${value}\")"
  done
}

# (Opcional) Ajustes gerais de execução:
export IO=1
export RUN=1
export CLEAN=0
export DEBUG=0
export RDC=0
export PROGRESS=1
export PROGRESS_HZ=2
export WARMUP=100
export REAL=float   # ou double se quiser

stencils=()
grids=()
res=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stencil)
      append_csv_values stencils "$2"
      shift 2
      ;;
    --grid)
      append_csv_values grids "$2"
      shift 2
      ;;
    --re)
      append_csv_values res "$2"
      shift 2
      ;;
    --t_final|--t-final|--t_star_end)
      T_FINAL="$2"
      shift 2
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    -h|--help)
      show_help
      exit 0
      ;;
    *)
      echo "Error: unknown argument '$1'" >&2
      show_help >&2
      exit 1
      ;;
  esac
done

[[ "${T_FINAL}" =~ ^[0-9]+$ ]] || {
  echo "Error: t_final must be numeric: '${T_FINAL}'" >&2
  exit 1
}

if ((${#stencils[@]} == 0)); then
  stencils=("D2Q9")
fi

# Tamanhos de grid. Use "256" para 256x256, ou "512x256" para retangular.
if ((${#grids[@]} == 0)); then
  grids=("512")
fi

if ((${#res[@]} == 0)); then
  res=("100" "400" "1000" "3200" "5000" "7500" "10000")
fi

for grid in "${grids[@]}"; do
  parse_grid "${grid}"

  for stencil in "${stencils[@]}"; do
    for re in "${res[@]}"; do

      ts="$(date +%Y%m%d_%H%M%S)"
      run_id="${ts}_${stencil}_${GRID_NX}x${GRID_NY}_RE${re}_T${T_FINAL}"

      echo "================================================="
      echo "[CASE] STENCIL=${stencil}  RE=${re}  GRID=${GRID_NX}x${GRID_NY}  T_FINAL=${T_FINAL}  RUN_ID=${run_id}"
      echo "================================================="

      cmd=(bash "${ROOT_DIR}/compile.sh" \
        --stencil "${stencil}" \
        --re "${re}" \
        --grid "${GRID_NX}x${GRID_NY}" \
        --t_star_end "${T_FINAL}" \
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
done

echo "ALL CASES DONE."
