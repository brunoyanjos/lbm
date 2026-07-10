#!/usr/bin/env bash
set -euo pipefail

# =====================================================
# Config (overridable via environment variables)
# =====================================================

: "${ARCHES:=}"
: "${BUILD_ROOT:=build}"
: "${CLEAN:=0}"
: "${DEBUG:=0}"
: "${DEVICE:=0}"
: "${EXEC_NAME:=sim}"
: "${GRID:=}"
: "${IO:=1}"
: "${MAXRREGCOUNT:=0}"
: "${NY:=128}"
: "${NX_MULT:=4}"
: "${OUT_ROOT:=runs}"
: "${PROGRESS:=1}"
: "${PROGRESS_HZ:=2}"
: "${PTXAS_VERBOSE:=0}"
: "${RDC:=0}"
: "${REAL:=float}"
: "${RUN:=1}"
: "${RE:=}"
: "${RUN_ID:=}"
: "${STENCIL:=D2Q9}"
: "${T_STAR_END:=1000}"
: "${VERBOSE:=0}"
: "${VTI_INTERVAL:=0}"
: "${WARMUP:=100}"

# =====================================================
# Helpers
# =====================================================

die() { echo "Error: $*" >&2; exit 1; }

detect_arches_from_system() {
  local caps
  caps="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | tr -d '.' | awk '/^[0-9]+$/ {print}')" || true
  if [[ -z "${caps}" ]]; then echo ""; return 0; fi
  echo "${caps}" | awk '!seen[$0]++' | sort -n | tr '\n' ' ' | sed 's/[[:space:]]*$//'
}

make_gencodes() {
  local arches="$1"
  GENCODES=()
  for a in ${arches}; do
    [[ "${a}" =~ ^[0-9]+$ ]] || die "ARCHES contains non-numeric entry: '${a}'"
    GENCODES+=(-gencode "arch=compute_${a},code=sm_${a}")
  done
  local lowest
  lowest="$(echo "${arches}" | tr ' ' '\n' | sort -n | head -n 1)"
  GENCODES+=(-gencode "arch=compute_${lowest},code=compute_${lowest}")
}

parse_grid() {
  local grid="$1"
  if [[ "${grid}" =~ ^([0-9]+)$ ]]; then
    GRID_NY="${BASH_REMATCH[1]}"
    GRID_NX=$((GRID_NY * NX_MULT))
  else
    die "Invalid grid '${grid}'. Expected base NY only, e.g. 256"
  fi
  ((GRID_NX > 0 && GRID_NY > 0)) || die "Grid dimensions must be positive: '${grid}'"
}

# =====================================================
# Parse CLI args
# =====================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stencil)        STENCIL="$2"; shift 2 ;;
    --real)           REAL="$2"; shift 2 ;;
    --re)             RE="$2"; shift 2 ;;
    --run)            RUN="$2"; shift 2 ;;
    --clean)          CLEAN="$2"; shift 2 ;;
    --debug)          DEBUG="$2"; shift 2 ;;
    --device)         DEVICE="$2"; shift 2 ;;
    --rdc)            RDC="$2"; shift 2 ;;
    --arches)         ARCHES="$2"; shift 2 ;;
    --exec)           EXEC_NAME="$2"; shift 2 ;;
    --grid)           GRID="$2"; shift 2 ;;
    --ny)             NY="$2"; shift 2 ;;
    --nx_mult|--nx-mult) NX_MULT="$2"; shift 2 ;;
    --build_root)     BUILD_ROOT="$2"; shift 2 ;;
    --out_root)       OUT_ROOT="$2"; shift 2 ;;
    --run_id)         RUN_ID="$2"; shift 2 ;;
    --io)             IO="$2"; shift 2 ;;
    --warmup)         WARMUP="$2"; shift 2 ;;
    --t_star_end)     T_STAR_END="$2"; shift 2 ;;
    --vti_interval)   VTI_INTERVAL="$2"; shift 2 ;;
    --verbose)        VERBOSE="$2"; shift 2 ;;
    --progress)       PROGRESS="$2"; shift 2 ;;
    --progress_hz)    PROGRESS_HZ="$2"; shift 2 ;;
    --ptxas_verbose)  PTXAS_VERBOSE="$2"; shift 2 ;;
    --maxrregcount)   MAXRREGCOUNT="$2"; shift 2 ;;
    -h|--help)
      cat <<EOF
Usage:
  STENCIL=D2Q9 REAL=float bash compile.sh
  bash compile.sh --stencil D2Q9 --real float --run 1
  bash compile.sh --ny 256
  bash compile.sh --grid 256
  bash compile.sh --ny 256 --nx_mult 4
  bash compile.sh --device 0
  bash compile.sh --t_star_end 1000
  bash compile.sh --vti_interval 10000
EOF
      exit 0
      ;;
    *)
      die "Unknown argument: '$1'"
      ;;
  esac
done

# =====================================================
# Validate
# =====================================================

[[ "${T_STAR_END}" =~ ^[0-9]+$ ]] || die "T_STAR_END must be numeric: '${T_STAR_END}'"

case "${STENCIL}" in D2Q9|D2V17|D2V37) ;; *) die "Unknown STENCIL='${STENCIL}'" ;; esac
case "${REAL}" in float|double) ;; *) die "Unknown REAL='${REAL}'" ;; esac
[[ "${NY}" =~ ^[0-9]+$ ]] || die "NY must be numeric: '${NY}'"
[[ "${NX_MULT}" =~ ^[0-9]+$ ]] || die "NX_MULT must be numeric: '${NX_MULT}'"
((NX_MULT > 0)) || die "NX_MULT must be positive: '${NX_MULT}'"
[[ "${DEVICE}" =~ ^[0-9]+$ ]] || die "DEVICE must be numeric: '${DEVICE}'"

if [[ -n "${GRID}" ]]; then
  parse_grid "${GRID}"
else
  parse_grid "${NY}"
fi

# =====================================================
# ARCHES
# =====================================================

if [[ -z "${ARCHES}" ]]; then ARCHES="$(detect_arches_from_system)"; fi
if [[ -z "${ARCHES}" ]]; then ARCHES="80 86 89 90"; fi
make_gencodes "${ARCHES}"

# =====================================================
# Build dirs
# =====================================================

CFG_TAG="${STENCIL}_${REAL}_${GRID_NX}x${GRID_NY}"
BUILD_DIR="${BUILD_ROOT}/${CFG_TAG}"
OBJ_DIR="${BUILD_DIR}/obj"
BIN_PATH="${BUILD_DIR}/${EXEC_NAME}"

echo "ARCHES: ${ARCHES}"
echo "STENCIL=${STENCIL} REAL=${REAL} GRID=${GRID_NX}x${GRID_NY} NX_MULT=${NX_MULT} T_STAR_END=${T_STAR_END}"
echo "DEVICE=${DEVICE}"
echo "DEBUG=${DEBUG} RDC=${RDC} CLEAN=${CLEAN}"
echo "BUILD_DIR: ${BUILD_DIR}"

if [[ "${CLEAN}" == "1" ]]; then
  echo "[CLEAN] rm -rf '${BUILD_DIR}'"
  rm -rf "${BUILD_DIR}"
fi

mkdir -p "${OBJ_DIR}"
rm -f "${BIN_PATH}"

# =====================================================
# Flags
# =====================================================

NVCCFLAGS=(-std=c++17 -lineinfo --restrict -Isrc)
NVCCFLAGS+=(-Xcompiler -Wall -Xcompiler -Wextra)

if [[ "${DEBUG}" == "1" ]]; then
  NVCCFLAGS+=(-O0 -g -G)
else
  NVCCFLAGS+=(-O3)
fi

if [[ "${RDC}" == "1" ]]; then
  NVCCFLAGS+=(-rdc=true)
else
  NVCCFLAGS+=(-rdc=false)
fi

case "${STENCIL}" in
  D2Q9)
    NVCCFLAGS+=(-DLBM_STENCIL_D2Q9)
    ;;
  D2V17)
    NVCCFLAGS+=(-DLBM_STENCIL_D2V17)
    ;;
  D2V37)
    NVCCFLAGS+=(-DLBM_STENCIL_D2V37)
    ;;
  *)
    die "Unknown STENCIL='${STENCIL}'"
    ;;
esac

NVCCFLAGS+=(-DLBM_NX="${GRID_NX}" -DLBM_NY="${GRID_NY}")

if [[ -n "${RE}" ]]; then
  NVCCFLAGS+=(-DLBM_RE="${RE}")
fi

NVCCFLAGS+=(-DLBM_T_STAR_END="${T_STAR_END}")

if [[ "${REAL}" == "double" ]]; then
  NVCCFLAGS+=(-DREAL_T_IS_DOUBLE)
fi

if [[ "${PTXAS_VERBOSE}" == "1" ]]; then
  NVCCFLAGS+=(-Xptxas -v)
fi

if [[ "${MAXRREGCOUNT}" != "0" ]]; then
  NVCCFLAGS+=(-maxrregcount="${MAXRREGCOUNT}")
fi

NVCCFLAGS+=("${GENCODES[@]}")

# =====================================================
# Sources
# =====================================================

mapfile -d '' CU_FILES < <(find src -type f -name "*.cu" -print0)
TOTAL=${#CU_FILES[@]}
echo "Sources: ${TOTAL}"

# =====================================================
# Compile
# =====================================================

OBJ_FILES=()
COUNT=0

for cu in "${CU_FILES[@]}"; do
  ((++COUNT))

  rel="${cu#./}"
  objbase="${rel//\//_}"
  obj="${OBJ_DIR}/${objbase%.cu}.o"
  OBJ_FILES+=("${obj}")

  printf "\r[NVCC %2d/%2d] %-70s" "${COUNT}" "${TOTAL}" "${cu}"

  if [[ "${RDC}" == "1" ]]; then
    nvcc "${NVCCFLAGS[@]}" -dc "${cu}" -o "${obj}"
  else
    nvcc "${NVCCFLAGS[@]}" -c "${cu}" -o "${obj}"
  fi
done

printf '\r\033[K'

# =====================================================
# Link
# =====================================================

echo "[LINK] ${BIN_PATH}"
LINKFLAGS=()
if [[ "${RDC}" == "1" ]]; then LINKFLAGS+=(-lcudadevrt); fi
nvcc "${NVCCFLAGS[@]}" "${OBJ_FILES[@]}" "${LINKFLAGS[@]}" -o "${BIN_PATH}"

echo "✔ Build successful: ${BIN_PATH}"

# =====================================================
# Run
# =====================================================

if [[ "${RUN}" == "1" ]]; then
  if [[ -z "${RUN_ID}" ]]; then
    CURRENT_RUN_ID="$(date +%Y%m%d_%H%M%S)_${STENCIL}_${REAL}_${GRID_NX}x${GRID_NY}${RE:+_RE${RE}}"
  else
    CURRENT_RUN_ID="${RUN_ID}"
  fi

  OUT_DIR="${OUT_ROOT}/${CURRENT_RUN_ID}"
  mkdir -p "${OUT_DIR}/vtk" "${OUT_DIR}/logs" "${OUT_DIR}/outputs"

  echo "[RUN] out_dir=${OUT_DIR}"

  CURRENT_PROGRESS="${PROGRESS}"
  if [[ "${CURRENT_PROGRESS}" == "1" && ! -t 2 ]]; then
    CURRENT_PROGRESS=0
  fi

  "${BIN_PATH}" \
    --out "${OUT_DIR}" \
    --device "${DEVICE}" \
    --io "${IO}" \
    --warmup "${WARMUP}" \
    --vti_interval "${VTI_INTERVAL}" \
    --verbose "${VERBOSE}" \
    --progress "${CURRENT_PROGRESS}" \
    --progress_hz "${PROGRESS_HZ}" \
    > "${OUT_DIR}/logs/stdout.txt"

  echo "✔ Run finished. Outputs in: ${OUT_DIR}"
fi
