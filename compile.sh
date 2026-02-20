#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Config (overridable)
# -------------------------------
: "${RDC:=0}"
: "${DEBUG:=0}"
: "${RUN:=1}"
: "${ARCHES:=}"

: "${STENCIL:=D2Q9}"     
: "${REAL:=float}"       
: "${CLEAN:=0}"          

: "${BUILD_ROOT:=build}"
: "${EXEC_NAME:=sim}"

: "${OUT_ROOT:=runs}"
: "${RUN_ID:=}"        

# -------------------------------
# Helpers
# -------------------------------
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
  gencodes+=(-gencode "arch=compute_${lowest},code=compute_${lowest}")  # PTX fallback
  printf '%s\n' "${gencodes[@]}"
}

# -------------------------------
# Depfile helpers (incremental rebuild with headers)
# -------------------------------
needs_rebuild() {
  local obj="$1"
  local dep="$2"
  local src="$3"

  [[ -f "$obj" ]] || return 0
  [[ -f "$dep" ]] || return 0
  [[ "$src" -nt "$obj" ]] && return 0

  local d
  while read -r d; do
    [[ -z "$d" ]] && continue
    [[ -f "$d" ]] || continue
    [[ "$d" -nt "$obj" ]] && return 0
  done < <(
    tr ' \\' '\n\n' < "$dep" | sed -e 's/:$//' -e '1d' -e '/^$/d'
  )

  return 1
}

# -------------------------------
# Validate
# -------------------------------
case "${STENCIL}" in D2Q9|D2V17) ;; *) die "Unknown STENCIL='${STENCIL}'" ;; esac
case "${REAL}" in float|double) ;; *) die "Unknown REAL='${REAL}'" ;; esac

# -------------------------------
# ARCHES
# -------------------------------
if [[ -z "${ARCHES}" ]]; then ARCHES="$(detect_arches_from_system)"; fi
if [[ -z "${ARCHES}" ]]; then ARCHES="80 86 89 90"; fi
mapfile -t GENCODES < <(make_gencodes "${ARCHES}")

# -------------------------------
# Build dirs
# -------------------------------
CFG_TAG="${STENCIL}_${REAL}"
BUILD_DIR="${BUILD_ROOT}/${CFG_TAG}"
OBJ_DIR="${BUILD_DIR}/obj"
BIN_PATH="${BUILD_DIR}/${EXEC_NAME}"

echo "ARCHES: ${ARCHES}"
echo "RDC=${RDC} DEBUG=${DEBUG} RUN=${RUN} STENCIL=${STENCIL} REAL=${REAL} CLEAN=${CLEAN}"
echo "BUILD_DIR: ${BUILD_DIR}"

if [[ "${CLEAN}" == "1" ]]; then
  echo "[CLEAN] rm -rf '${BUILD_DIR}'"
  rm -rf "${BUILD_DIR}"
fi

mkdir -p "${OBJ_DIR}"
rm -f "${BIN_PATH}"

# -------------------------------
# Flags
# -------------------------------
NVCCFLAGS=(-std=c++17 -lineinfo --restrict -Isrc)
NVCCFLAGS+=(-Xcompiler -Wall -Xcompiler -Wextra)

if [[ "${DEBUG}" == "1" ]]; then NVCCFLAGS+=(-O0 -g -G); else NVCCFLAGS+=(-O3); fi
if [[ "${RDC}" == "1" ]]; then NVCCFLAGS+=(-rdc=true); else NVCCFLAGS+=(-rdc=false); fi

if [[ "${STENCIL}" == "D2Q9" ]]; then NVCCFLAGS+=(-DLBM_STENCIL_D2Q9); else NVCCFLAGS+=(-DLBM_STENCIL_D2V17); fi
if [[ "${REAL}" == "double" ]]; then NVCCFLAGS+=(-DREAL_T_IS_DOUBLE); fi

NVCCFLAGS+=("${GENCODES[@]}")

# -------------------------------
# Sources (main agora está em src/)
# -------------------------------
mapfile -d '' CU_FILES < <(find src -type f -name "*.cu" -print0)
echo "Sources: ${#CU_FILES[@]}"

# -------------------------------
# Compile
# -------------------------------
OBJ_FILES=()
for cu in "${CU_FILES[@]}"; do
  rel="${cu#./}"
  objbase="${rel//\//_}"
  obj="${OBJ_DIR}/${objbase%.cu}.o"
  dep="${obj%.o}.d"
  OBJ_FILES+=("${obj}")

  if needs_rebuild "${obj}" "${dep}" "${cu}"; then
    echo "[NVCC] ${cu} -> ${obj}"
    DEPFLAGS=(-MMD -MP -MF "${dep}" -MT "${obj}")

    if [[ "${RDC}" == "1" ]]; then
      nvcc "${NVCCFLAGS[@]}" "${DEPFLAGS[@]}" -dc "${cu}" -o "${obj}"
    else
      nvcc "${NVCCFLAGS[@]}" "${DEPFLAGS[@]}" -c  "${cu}" -o "${obj}"
    fi
  fi
done

# -------------------------------
# Link
# -------------------------------
echo "[LINK] ${BIN_PATH}"
LINKFLAGS=()
if [[ "${RDC}" == "1" ]]; then LINKFLAGS+=(-lcudadevrt); fi
nvcc "${NVCCFLAGS[@]}" "${OBJ_FILES[@]}" "${LINKFLAGS[@]}" -o "${BIN_PATH}"
echo "✔ Build successful: ${BIN_PATH}"

# -------------------------------
# Run
# -------------------------------
if [[ "${RUN}" == "1" ]]; then
  if [[ -z "${RUN_ID}" ]]; then
    RUN_ID="$(date +%Y%m%d_%H%M%S)"
  fi

  OUT_DIR="${OUT_ROOT}/${RUN_ID}"
  mkdir -p "${OUT_DIR}/vtk" "${OUT_DIR}/logs"

  echo "[RUN] out_dir=${OUT_DIR}"

  "${BIN_PATH}" --out "${OUT_DIR}" \
    > "${OUT_DIR}/logs/stdout.txt"

  echo "✔ Run finished. Outputs in: ${OUT_DIR}"
fi