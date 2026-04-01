#!/usr/bin/env bash
set -euo pipefail

# =====================================================
# Config (overridable via environment variables)
# =====================================================

: "${ARCHES:=}"
: "${BUILD_ROOT:=build}"
: "${CLEAN:=0}"
: "${DEBUG:=0}"
: "${EXEC_NAME:=sim}"
: "${IO:=1}"
: "${MAXRREGCOUNT:=0}"
: "${OUT_ROOT:=runs}"
: "${PROGRESS:=1}"
: "${PROGRESS_HZ:=2}"
: "${PTXAS_VERBOSE:=0}"
: "${RDC:=0}"
: "${REAL:=float}"
: "${RUN:=1}"
: "${GRID:=512}"
: "${RE:=}"
: "${RUN_ID:=}"
: "${VERBOSE:=0}"
: "${WARMUP:=100}"

# =====================================================
# Helpers
# =====================================================

die() { echo "Error: $*" >&2; exit 1; }

detect_arches_from_system() {
  local caps
  caps="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | tr -d '.[:space:]')" || true
  if [[ -z "${caps}" ]]; then echo ""; return 0; fi
  echo "${caps}" | awk '!seen[$0]++' | sort -n | tr '\n' ' ' | sed 's/[[:space:]]*$//'
}

make_gencodes() {
  local arches="$1"
  local -a gencodes=()
  for a in ${arches}; do
    [[ "${a}" =~ ^[0-9]+$ ]] || die "ARCHES contains non-numeric entry: '${a}'"
    gencodes+=(-gencode "arch=compute_${a},code=sm_${a}")
  done
  local lowest
  lowest="$(echo "${arches}" | tr ' ' '\n' | sort -n | head -n 1)"
  gencodes+=(-gencode "arch=compute_${lowest},code=compute_${lowest}")
  printf '%s\n' "${gencodes[@]}"
}

# =====================================================
# Parse CLI args
# =====================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    --real)           REAL="$2"; shift 2 ;;
    --re)             RE="$2"; shift 2 ;;
    --run)            RUN="$2"; shift 2 ;;
    --clean)          CLEAN="$2"; shift 2 ;;
    --debug)          DEBUG="$2"; shift 2 ;;
    --rdc)            RDC="$2"; shift 2 ;;
    --arches)         ARCHES="$2"; shift 2 ;;
    --exec)           EXEC_NAME="$2"; shift 2 ;;
    --build_root)     BUILD_ROOT="$2"; shift 2 ;;
    --out_root)       OUT_ROOT="$2"; shift 2 ;;
    --reg_order)      REG_ORDER="$2"; shift 2 ;;
    --run_id)         RUN_ID="$2"; shift 2 ;;
    --io)             IO="$2"; shift 2 ;;
    --warmup)         WARMUP="$2"; shift 2 ;;
    --verbose)        VERBOSE="$2"; shift 2 ;;
    --progress)       PROGRESS="$2"; shift 2 ;;
    --progress_hz)    PROGRESS_HZ="$2"; shift 2 ;;
    --ptxas_verbose)  PTXAS_VERBOSE="$2"; shift 2 ;;
    --maxrregcount)   MAXRREGCOUNT="$2"; shift 2 ;;
    --grid)           GRID="$2"; shift 2 ;;
    -h|--help)
      cat <<EOF
Usage:
  REAL=float bash compile.sh
  bash compile.sh --real float --run 1
  bash compile.sh --grid 512
  GRID=1024 bash compile.sh
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

case "${REAL}" in
  float|double) ;;
  *)
    die "Unknown REAL='${REAL}'"
    ;;
esac

case "${GRID}" in
  ''|*[!0-9]*)
    die "GRID must be a positive integer, got '${GRID}'"
    ;;
  *)
    ;;
esac

if [[ "${GRID}" -le 0 ]]; then
  die "GRID must be > 0"
fi

# =====================================================
# ARCHES
# =====================================================

if [[ -z "${ARCHES}" ]]; then ARCHES="$(detect_arches_from_system)"; fi
if [[ -z "${ARCHES}" ]]; then ARCHES="80 86 89 90"; fi
mapfile -t GENCODES < <(make_gencodes "${ARCHES}")

# =====================================================
# Build dirs
# =====================================================

CFG_TAG="D2Q9_G${GRID}_${REAL}"
BUILD_DIR="${BUILD_ROOT}/${CFG_TAG}"
OBJ_DIR="${BUILD_DIR}/obj"
BIN_PATH="${BUILD_DIR}/${EXEC_NAME}"

echo "ARCHES: ${ARCHES}"
echo "STENCIL=D2Q9 REAL=${REAL}"
echo "GRID=${GRID}"
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
NVCCFLAGS+=(-DLBM_GRID="${GRID}")

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

if [[ -n "${RE}" ]]; then
  NVCCFLAGS+=(-DLBM_RE="${RE}")
fi

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
if [[ "${RDC}" == "1" ]]; then
  LINKFLAGS+=(-lcudadevrt)
fi
nvcc "${NVCCFLAGS[@]}" "${OBJ_FILES[@]}" "${LINKFLAGS[@]}" -o "${BIN_PATH}"

echo "✔ Build successful: ${BIN_PATH}"

# =====================================================
# Run
# =====================================================

if [[ "${RUN}" == "1" ]]; then
  if [[ -z "${RUN_ID}" ]]; then
    RUN_ID="$(date +%Y%m%d_%H%M%S)_D2Q9_G${GRID}_annul_rbc"
  fi

  OUT_DIR="${OUT_ROOT}/${RUN_ID}"
  mkdir -p "${OUT_DIR}/vtk" "${OUT_DIR}/logs" "${OUT_DIR}/outputs"

  echo "[RUN] out_dir=${OUT_DIR}"

  "${BIN_PATH}" \
    --out "${OUT_DIR}" \
    --io "${IO}" \
    --warmup "${WARMUP}" \
    --verbose "${VERBOSE}" \
    --progress "${PROGRESS}" \
    --progress_hz "${PROGRESS_HZ}" \
    > "${OUT_DIR}/logs/stdout.txt"

  echo "✔ Run finished. Outputs in: ${OUT_DIR}"
fi