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
: "${GEOM:=square_cavity}"
: "${IO:=1}"
: "${MAXRREGCOUNT:=0}"
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
# Parse CLI args (optional)
# =====================================================
# Mantém compatível com o uso por env vars, mas permite:
#   ./compile.sh --geom annul --stencil D2Q9 --real double --clean 1 --run 0 ...
while [[ $# -gt 0 ]]; do
  case "$1" in
    --geom|-g)        GEOM="$2"; shift 2 ;;
    --stencil)        STENCIL="$2"; shift 2 ;;
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
    --run_id)         RUN_ID="$2"; shift 2 ;;
    --io)             IO="$2"; shift 2 ;;
    --warmup)         WARMUP="$2"; shift 2 ;;
    --verbose)        VERBOSE="$2"; shift 2 ;;
    --progress)       PROGRESS="$2"; shift 2 ;;
    --progress_hz)    PROGRESS_HZ="$2"; shift 2 ;;
    --ptxas_verbose)  PTXAS_VERBOSE="$2"; shift 2 ;;
    --maxrregcount)   MAXRREGCOUNT="$2"; shift 2 ;;
    -h|--help)
      cat <<EOF
Usage:
  GEOM=square_cavity STENCIL=D2Q9 REAL=float ./compile.sh
  ./compile.sh --geom annul --stencil D2Q9 --real float --run 1

Valid:
  --geom: annul | channel | couette | jet | square_cavity
  --stencil: D2Q9 | D2V17
  --real: float | double
EOF
      exit 0
      ;;
    *)
      die "Unknown argument: '$1' (use --help)"
      ;;
  esac
done

# =====================================================
# Validate
# =====================================================

case "${STENCIL}" in D2Q9|D2V17) ;; *) die "Unknown STENCIL='${STENCIL}'" ;; esac
case "${REAL}" in float|double) ;; *) die "Unknown REAL='${REAL}'" ;; esac
case "${GEOM}" in annul|channel|couette|jet|square_cavity) ;; *) die "Unknown GEOM='${GEOM}'" ;; esac

# =====================================================
# GEOM -> macro
# =====================================================

GEOM_DEF=""
case "${GEOM}" in
  annul)         GEOM_DEF="-DLBM_GEOM_ANNUL" ;;
  channel)       GEOM_DEF="-DLBM_GEOM_CHANNEL" ;;
  couette)       GEOM_DEF="-DLBM_GEOM_COUETTE" ;;
  jet)           GEOM_DEF="-DLBM_GEOM_JET" ;;
  square_cavity) GEOM_DEF="-DLBM_GEOM_SQUARE_CAVITY" ;;
esac

# =====================================================
# ARCHES
# =====================================================

if [[ -z "${ARCHES}" ]]; then ARCHES="$(detect_arches_from_system)"; fi
if [[ -z "${ARCHES}" ]]; then ARCHES="80 86 89 90"; fi
mapfile -t GENCODES < <(make_gencodes "${ARCHES}")

# =====================================================
# Build dirs
# =====================================================

CFG_TAG="${STENCIL}_${REAL}_${GEOM}"
BUILD_DIR="${BUILD_ROOT}/${CFG_TAG}"
OBJ_DIR="${BUILD_DIR}/obj"
BIN_PATH="${BUILD_DIR}/${EXEC_NAME}"

echo "ARCHES: ${ARCHES}"
echo "STENCIL=${STENCIL} REAL=${REAL} GEOM=${GEOM}"
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

# stencil macro
if [[ "${STENCIL}" == "D2Q9" ]]; then
  NVCCFLAGS+=(-DLBM_STENCIL_D2Q9)
else
  NVCCFLAGS+=(-DLBM_STENCIL_D2V17)
fi

# geometry macro (NOVO)
NVCCFLAGS+=("${GEOM_DEF}")

# RE macro (se você usa isso em compile-time)
if [[ -n "${RE}" ]]; then
  NVCCFLAGS+=(-DLBM_RE="${RE}")
fi

# real macro
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
echo "Sources: ${#CU_FILES[@]}"

# =====================================================
# Compile
# =====================================================

OBJ_FILES=()

for cu in "${CU_FILES[@]}"; do
  rel="${cu#./}"
  objbase="${rel//\//_}"
  obj="${OBJ_DIR}/${objbase%.cu}.o"
  OBJ_FILES+=("${obj}")

  echo "[NVCC] ${cu}"
  if [[ "${RDC}" == "1" ]]; then
    nvcc "${NVCCFLAGS[@]}" -dc "${cu}" -o "${obj}"
  else
    nvcc "${NVCCFLAGS[@]}" -c "${cu}" -o "${obj}"
  fi
done

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
    RUN_ID="$(date +%Y%m%d_%H%M%S)_${STENCIL}_${REAL}_${GEOM}${RE:+_RE${RE}}"
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