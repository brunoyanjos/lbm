#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

: "${DRY_RUN:=0}"
: "${T_FINAL:=2000}"
: "${NX_MULT:=4}"

show_help() {
  cat <<EOF
Usage:
  bash tools/run_queue.sh [options]

Options:
  --stencil VALUE[,VALUE...]  Stencil(s) to run (default: D2Q9)
  --grid VALUE[,VALUE...]     Base NY size(s) (default: 512)
  --ny VALUE[,VALUE...]       Alias for --grid
  --nx_mult VALUE             NX multiplier, NX = VALUE * NY (default: ${NX_MULT})
  --re VALUE[,VALUE...]       Reynolds number(s) (default: 100,400,1000,3200,5000,7500,10000)
  --device VALUE[,VALUE...]   CUDA device id(s), assigned round-robin (default: 0)
  --t_final VALUE             Final t* passed to compile.sh --t_star_end (default: ${T_FINAL})
  --dry-run                   Print commands without running them
  -h, --help                  Show this help
EOF
}

parse_grid() {
  local grid="$1"
  if [[ "${grid}" =~ ^([0-9]+)$ ]]; then
    GRID_NY="${BASH_REMATCH[1]}"
    GRID_NX=$((GRID_NY * NX_MULT))
  else
    echo "Error: invalid grid '${grid}'. Expected base NY only, e.g. 256" >&2
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
devices=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stencil)
      append_csv_values stencils "$2"
      shift 2
      ;;
    --grid|--ny)
      append_csv_values grids "$2"
      shift 2
      ;;
    --nx_mult|--nx-mult)
      NX_MULT="$2"
      shift 2
      ;;
    --re)
      append_csv_values res "$2"
      shift 2
      ;;
    --device|--devices)
      append_csv_values devices "$2"
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

[[ "${NX_MULT}" =~ ^[0-9]+$ ]] || {
  echo "Error: nx_mult must be numeric: '${NX_MULT}'" >&2
  exit 1
}
((NX_MULT > 0)) || {
  echo "Error: nx_mult must be positive: '${NX_MULT}'" >&2
  exit 1
}

if ((${#stencils[@]} == 0)); then
  stencils=("D2Q9")
fi

# Tamanho base. O domínio final é (NX_MULT * NY) x NY.
if ((${#grids[@]} == 0)); then
  grids=("512")
fi

if ((${#res[@]} == 0)); then
  res=("100" "400" "1000" "3200" "5000" "7500" "10000")
fi

if ((${#devices[@]} == 0)); then
  devices=("0")
fi

for device in "${devices[@]}"; do
  [[ "${device}" =~ ^[0-9]+$ ]] || {
    echo "Error: device must be numeric: '${device}'" >&2
    exit 1
  }
done

case_index=0
for grid in "${grids[@]}"; do
  parse_grid "${grid}"

  for stencil in "${stencils[@]}"; do
    for re in "${res[@]}"; do
      device="${devices[$((case_index % ${#devices[@]}))]}"
      ((case_index += 1))

      ts="$(date +%Y%m%d_%H%M%S)"
      run_id="${ts}_${stencil}_${GRID_NX}x${GRID_NY}_RE${re}_T${T_FINAL}_GPU${device}"

      echo "================================================="
      echo "[CASE] STENCIL=${stencil}  RE=${re}  GRID=${GRID_NX}x${GRID_NY}  NX_MULT=${NX_MULT}  T_FINAL=${T_FINAL}  DEVICE=${device}  RUN_ID=${run_id}"
      echo "================================================="

      cmd=(bash "${ROOT_DIR}/compile.sh" \
        --stencil "${stencil}" \
        --re "${re}" \
        --ny "${GRID_NY}" \
        --nx_mult "${NX_MULT}" \
        --device "${device}" \
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
