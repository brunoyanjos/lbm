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
: "${NX:=128}"
: "${NY:=128}"
: "${OUT_ROOT:=runs}"
: "${PROGRESS:=1}"
: "${PROGRESS_HZ:=2}"
: "${PTXAS_VERBOSE:=0}"
: "${RDC:=0}"
: "${REAL:=float}"
: "${RECURRENCE:=0}"
: "${REG_ORDER:=2}"
: "${RUN:=1}"
: "${SYMBOLIC_BOUNDARY:=0}"
: "${RE:=}"
: "${RESTART:=0}"
: "${CHECKPOINT_RUN_ID:=}"
: "${CHECKPOINT_DIR:=}"
: "${CHECKPOINT_T_STAR:=}"
: "${CHECKPOINT_FILE:=}"
: "${CHECKPOINT_N_STEPS:=}"
: "${CHECKPOINT_SAVE_INTERVAL:=}"
: "${CHECKPOINT_VTI_SAVE_INTERVAL:=}"
: "${CHECKPOINT_STEP:=}"
: "${RUN_ID:=}"
: "${STENCIL:=D2Q9}"
: "${T_STAR_END:=1000}"
: "${AVG_START_T_STAR:=0}"
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
    GRID_NX="${BASH_REMATCH[1]}"
    GRID_NY="${BASH_REMATCH[1]}"
  elif [[ "${grid}" =~ ^([0-9]+)[xX]([0-9]+)$ ]]; then
    GRID_NX="${BASH_REMATCH[1]}"
    GRID_NY="${BASH_REMATCH[2]}"
  else
    die "Invalid grid '${grid}'. Expected N or NXxNY, e.g. 256 or 256x256"
  fi
  ((GRID_NX > 0 && GRID_NY > 0)) || die "Grid dimensions must be positive: '${grid}'"
}

import_checkpoint_build_config() {
  local checkpoint_path="$1"
  local checkpoint_t_star="$2"
  local parsed=()

  mapfile -t parsed < <(python3 - "${checkpoint_path}" "${checkpoint_t_star}" <<'PY'
import pathlib
import struct
import sys

path = pathlib.Path(sys.argv[1])
requested_t_star = sys.argv[2]

fmt = "<8sIII" + ("i" * 9) + ("q" * 7) + ("d" * 4) + "I16s"

def read_header(candidate):
    with candidate.open("rb") as f:
        data = f.read(struct.calcsize(fmt))

    if len(data) != struct.calcsize(fmt):
        raise SystemExit(f"Could not read checkpoint header: {candidate}")

    values = struct.unpack(fmt, data)
    magic = values[0]
    if magic != b"LBMCHK1\0":
        raise SystemExit(f"Invalid checkpoint magic in: {candidate}")

    version, header_bytes, endian = values[1:4]
    if version != 2:
        raise SystemExit(f"Unsupported checkpoint version {version} in: {candidate}")
    if endian != 0x01020304:
        raise SystemExit(f"Unsupported checkpoint endian marker in: {candidate}")

    return values

if path.is_dir():
    candidates = sorted(path.glob("checkpoint*.bin"))
    if not candidates:
        raise SystemExit(f"No checkpoint*.bin files found in: {path}")

    if requested_t_star:
        requested = int(requested_t_star)
        selected = None
        selected_values = None
        for candidate in candidates:
            values = read_header(candidate)
            step, save_interval = values[16], values[18]
            t_star = step // save_interval
            if t_star == requested:
                selected = candidate
                selected_values = values
                break

        if selected is None:
            raise SystemExit(f"No checkpoint found for t_star={requested} in: {path}")

        path = selected
        values = selected_values
    else:
        path = candidates[-1]
        values = read_header(path)
elif path.is_file():
    values = read_header(path)
else:
    raise SystemExit(f"Checkpoint path does not exist: {path}")

nx, ny, q, real_bytes, real_is_double, _cur, _field_count, reg_order, recurrence = values[4:13]
node_count, field_bytes, payload_bytes, step, n_steps, save_interval, vti_save_interval = values[13:20]
re, _u_lid, _tau, _omega = values[20:24]
stencil_name = values[25].split(b"\0", 1)[0].decode("ascii")
checkpoint_t_star = step // save_interval

real = "double" if real_is_double else "float"

print(f"CHECKPOINT_FILE={path}")
print(f"STENCIL={stencil_name}")
print(f"REAL={real}")
print(f"NX={nx}")
print(f"NY={ny}")
print(f"RE={re:.17g}")
print(f"REG_ORDER={reg_order}")
print(f"RECURRENCE={recurrence}")
print(f"CHECKPOINT_N_STEPS={n_steps}")
print(f"CHECKPOINT_SAVE_INTERVAL={save_interval}")
print(f"CHECKPOINT_VTI_SAVE_INTERVAL={vti_save_interval}")
print(f"CHECKPOINT_STEP={step}")
print(f"CHECKPOINT_T_STAR={checkpoint_t_star}")
print(f"CHECKPOINT_Q={q}")
print(f"CHECKPOINT_REAL_BYTES={real_bytes}")
print(f"CHECKPOINT_NODE_COUNT={node_count}")
print(f"CHECKPOINT_FIELD_BYTES={field_bytes}")
print(f"CHECKPOINT_PAYLOAD_BYTES={payload_bytes}")
PY
  )

  for kv in "${parsed[@]}"; do
    local key="${kv%%=*}"
    local value="${kv#*=}"
    case "${key}" in
      CHECKPOINT_FILE) CHECKPOINT_FILE="${value}" ;;
      STENCIL) STENCIL="${value}" ;;
      REAL) REAL="${value}" ;;
      NX) NX="${value}" ;;
      NY) NY="${value}" ;;
      RE) RE="${value}" ;;
      REG_ORDER) REG_ORDER="${value}" ;;
      RECURRENCE) RECURRENCE="${value}" ;;
      CHECKPOINT_N_STEPS) CHECKPOINT_N_STEPS="${value}" ;;
      CHECKPOINT_SAVE_INTERVAL) CHECKPOINT_SAVE_INTERVAL="${value}" ;;
      CHECKPOINT_VTI_SAVE_INTERVAL) CHECKPOINT_VTI_SAVE_INTERVAL="${value}" ;;
      CHECKPOINT_STEP) CHECKPOINT_STEP="${value}" ;;
      CHECKPOINT_T_STAR) CHECKPOINT_T_STAR="${value}" ;;
      CHECKPOINT_Q) CHECKPOINT_Q="${value}" ;;
      CHECKPOINT_REAL_BYTES) CHECKPOINT_REAL_BYTES="${value}" ;;
      CHECKPOINT_NODE_COUNT) CHECKPOINT_NODE_COUNT="${value}" ;;
      CHECKPOINT_FIELD_BYTES) CHECKPOINT_FIELD_BYTES="${value}" ;;
      CHECKPOINT_PAYLOAD_BYTES) CHECKPOINT_PAYLOAD_BYTES="${value}" ;;
    esac
  done
}

# =====================================================
# Parse CLI args
# =====================================================

while [[ $# -gt 0 ]]; do
  case "$1" in
    --stencil)        STENCIL="$2"; shift 2 ;;
    --real)           REAL="$2"; shift 2 ;;
    --re)             RE="$2"; shift 2 ;;
    --restart)        RESTART="$2"; shift 2 ;;
    --checkpoint_run_id) CHECKPOINT_RUN_ID="$2"; shift 2 ;;
    --checkpoint_dir) CHECKPOINT_DIR="$2"; shift 2 ;;
    --checkpoint_t_star) CHECKPOINT_T_STAR="$2"; shift 2 ;;
    --run)            RUN="$2"; shift 2 ;;
    --clean)          CLEAN="$2"; shift 2 ;;
    --debug)          DEBUG="$2"; shift 2 ;;
    --device)         DEVICE="$2"; shift 2 ;;
    --rdc)            RDC="$2"; shift 2 ;;
    --arches)         ARCHES="$2"; shift 2 ;;
    --exec)           EXEC_NAME="$2"; shift 2 ;;
    --grid)           GRID="$2"; shift 2 ;;
    --nx)             NX="$2"; shift 2 ;;
    --ny)             NY="$2"; shift 2 ;;
    --build_root)     BUILD_ROOT="$2"; shift 2 ;;
    --out_root)       OUT_ROOT="$2"; shift 2 ;;
    --reg_order)      REG_ORDER="$2"; shift 2 ;;
    --recurrence)     RECURRENCE="$2"; shift 2 ;;
    --symbolic_boundary|--symbolic-boundary) SYMBOLIC_BOUNDARY="$2"; shift 2 ;;
    --run_id)         RUN_ID="$2"; shift 2 ;;
    --io)             IO="$2"; shift 2 ;;
    --warmup)         WARMUP="$2"; shift 2 ;;
    --t_star_end)     T_STAR_END="$2"; shift 2 ;;
    --avg_start_t_star) AVG_START_T_STAR="$2"; shift 2 ;;
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
  bash compile.sh --grid 256
  bash compile.sh --device 0
  bash compile.sh --grid 256x128
  bash compile.sh --t_star_end 1000
  bash compile.sh --avg_start_t_star 500
  bash compile.sh --vti_interval 10000
  bash compile.sh --reg_order 3 --recurrence 1
  bash compile.sh --symbolic_boundary 1
  bash compile.sh --restart 1 --checkpoint_run_id 20260525_120000_D2Q9_float_128x128 --checkpoint_t_star 900 --t_star_end 2000

Regularization:
  --reg_order 2|3              Regularization order (default: 2)
  --recurrence 0|1             Use recurrence for third-order moments (default: 0)
  --symbolic_boundary 0|1      Use symbolic boundary solver path (default: 0)
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

case "${RESTART}" in 0|1) ;; *) die "RESTART must be 0 or 1: '${RESTART}'" ;; esac
[[ "${T_STAR_END}" =~ ^[0-9]+$ ]] || die "T_STAR_END must be numeric: '${T_STAR_END}'"
[[ "${AVG_START_T_STAR}" =~ ^[0-9]+$ ]] || die "AVG_START_T_STAR must be numeric: '${AVG_START_T_STAR}'"
if [[ -n "${CHECKPOINT_T_STAR}" ]]; then
  [[ "${CHECKPOINT_T_STAR}" =~ ^[0-9]+$ ]] || die "CHECKPOINT_T_STAR must be numeric: '${CHECKPOINT_T_STAR}'"
fi

if [[ "${RESTART}" == "1" ]]; then
  if [[ -z "${CHECKPOINT_RUN_ID}" && -z "${CHECKPOINT_DIR}" ]]; then
    die "Restart requires --checkpoint_run_id <run_id> or --checkpoint_dir <path>"
  fi

  if [[ -z "${CHECKPOINT_DIR}" ]]; then
    CHECKPOINT_DIR="${OUT_ROOT}/${CHECKPOINT_RUN_ID}/checkpoints"
  fi

  import_checkpoint_build_config "${CHECKPOINT_DIR}" "${CHECKPOINT_T_STAR}"
  [[ -n "${CHECKPOINT_FILE}" ]] || die "Could not import checkpoint build config from '${CHECKPOINT_DIR}'"
  GRID="${NX}x${NY}"
fi

case "${STENCIL}" in D2Q9|D2V17|D2V37) ;; *) die "Unknown STENCIL='${STENCIL}'" ;; esac
case "${REAL}" in float|double) ;; *) die "Unknown REAL='${REAL}'" ;; esac
case "${REG_ORDER}" in 2|3) ;; *) die "REG_ORDER must be 2 or 3: '${REG_ORDER}'" ;; esac
case "${RECURRENCE}" in 0|1) ;; *) die "RECURRENCE must be 0 or 1: '${RECURRENCE}'" ;; esac
case "${SYMBOLIC_BOUNDARY}" in 0|1) ;; *) die "SYMBOLIC_BOUNDARY must be 0 or 1: '${SYMBOLIC_BOUNDARY}'" ;; esac
if [[ "${RECURRENCE}" == "1" && "${REG_ORDER}" != "3" ]]; then
  die "RECURRENCE=1 requires REG_ORDER=3"
fi
[[ "${NX}" =~ ^[0-9]+$ ]] || die "NX must be numeric: '${NX}'"
[[ "${NY}" =~ ^[0-9]+$ ]] || die "NY must be numeric: '${NY}'"
[[ "${DEVICE}" =~ ^[0-9]+$ ]] || die "DEVICE must be numeric: '${DEVICE}'"

if [[ -n "${GRID}" ]]; then
  parse_grid "${GRID}"
else
  parse_grid "${NX}x${NY}"
fi

if [[ "${RESTART}" == "1" ]]; then
  CHECKPOINT_N_STEPS=$((T_STAR_END * CHECKPOINT_SAVE_INTERVAL))
  if (( CHECKPOINT_N_STEPS <= CHECKPOINT_STEP )); then
    die "t_star_end=${T_STAR_END} must be greater than checkpoint_t_star=${CHECKPOINT_T_STAR}"
  fi
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

REC_TAG=""
if [[ "${RECURRENCE}" == "1" ]]; then
  REC_TAG="_rec"
fi

SYM_TAG=""
if [[ "${SYMBOLIC_BOUNDARY}" == "1" ]]; then
  SYM_TAG="_symbc"
fi

CFG_TAG="${STENCIL}_${REAL}_${GRID_NX}x${GRID_NY}_reg${REG_ORDER}${REC_TAG}${SYM_TAG}"
BUILD_DIR="${BUILD_ROOT}/${CFG_TAG}"
OBJ_DIR="${BUILD_DIR}/obj"
BIN_PATH="${BUILD_DIR}/${EXEC_NAME}"

echo "ARCHES: ${ARCHES}"
echo "STENCIL=${STENCIL} REAL=${REAL} GRID=${GRID_NX}x${GRID_NY} T_STAR_END=${T_STAR_END} AVG_START_T_STAR=${AVG_START_T_STAR}"
echo "REG_ORDER=${REG_ORDER} RECURRENCE=${RECURRENCE} SYMBOLIC_BOUNDARY=${SYMBOLIC_BOUNDARY}"
echo "DEVICE=${DEVICE}"
if [[ "${RESTART}" == "1" ]]; then
  echo "RESTART=1 CHECKPOINT_RUN_ID=${CHECKPOINT_RUN_ID:-<direct-dir>} CHECKPOINT_DIR=${CHECKPOINT_DIR}"
  echo "CHECKPOINT_FILE=${CHECKPOINT_FILE} STEP=${CHECKPOINT_STEP} T_STAR=${CHECKPOINT_T_STAR}"
  echo "CHECKPOINT_COMPILE RE=${RE} N_STEPS=${CHECKPOINT_N_STEPS} SAVE_INTERVAL=${CHECKPOINT_SAVE_INTERVAL} VTI_SAVE_INTERVAL=${CHECKPOINT_VTI_SAVE_INTERVAL}"
fi
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
NVCCFLAGS+=(-DLBM_REG_ORDER="${REG_ORDER}" -DLBM_USE_RECURRENCE="${RECURRENCE}")
NVCCFLAGS+=(-DLBM_USE_SYMBOLIC_BOUNDARY="${SYMBOLIC_BOUNDARY}")

if [[ -n "${RE}" ]]; then
  NVCCFLAGS+=(-DLBM_RE="${RE}")
fi

NVCCFLAGS+=(-DLBM_T_STAR_END="${T_STAR_END}")
NVCCFLAGS+=(-DLBM_AVG_START_T_STAR="${AVG_START_T_STAR}")

if [[ -n "${CHECKPOINT_N_STEPS}" ]]; then
  NVCCFLAGS+=(-DLBM_N_STEPS="${CHECKPOINT_N_STEPS}")
fi

if [[ -n "${CHECKPOINT_SAVE_INTERVAL}" ]]; then
  NVCCFLAGS+=(-DLBM_SAVE_INTERVAL="${CHECKPOINT_SAVE_INTERVAL}")
fi

if [[ -n "${CHECKPOINT_VTI_SAVE_INTERVAL}" ]]; then
  NVCCFLAGS+=(-DLBM_VTI_SAVE_INTERVAL="${CHECKPOINT_VTI_SAVE_INTERVAL}")
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
if [[ "${RDC}" == "1" ]]; then LINKFLAGS+=(-lcudadevrt); fi
nvcc "${NVCCFLAGS[@]}" "${OBJ_FILES[@]}" "${LINKFLAGS[@]}" -o "${BIN_PATH}"

echo "✔ Build successful: ${BIN_PATH}"

# =====================================================
# Run
# =====================================================

if [[ "${RUN}" == "1" ]]; then
  if [[ -z "${RUN_ID}" ]]; then
    if [[ "${RESTART}" == "1" ]]; then
      CURRENT_RUN_ID="$(date +%Y%m%d_%H%M%S)_restart_${CHECKPOINT_RUN_ID:-checkpoint}"
    else
      CURRENT_RUN_ID="$(date +%Y%m%d_%H%M%S)_${STENCIL}_${REAL}_${GRID_NX}x${GRID_NY}_reg${REG_ORDER}${REC_TAG}${SYM_TAG}${RE:+_RE${RE}}"
    fi
  else
    CURRENT_RUN_ID="${RUN_ID}"
  fi

  OUT_DIR="${OUT_ROOT}/${CURRENT_RUN_ID}"
  mkdir -p "${OUT_DIR}/vtk" "${OUT_DIR}/logs" "${OUT_DIR}/outputs" "${OUT_DIR}/checkpoints"

  echo "[RUN] out_dir=${OUT_DIR}"

  CURRENT_PROGRESS="${PROGRESS}"
  if [[ "${CURRENT_PROGRESS}" == "1" && ! -t 2 ]]; then
    CURRENT_PROGRESS=0
  fi

  "${BIN_PATH}" \
    --out "${OUT_DIR}" \
    --device "${DEVICE}" \
    --restart "${RESTART}" \
    --checkpoint_run_id "${CHECKPOINT_RUN_ID}" \
    --checkpoint_dir "${CHECKPOINT_FILE:-${CHECKPOINT_DIR}}" \
    --io "${IO}" \
    --warmup "${WARMUP}" \
    --vti_interval "${VTI_INTERVAL}" \
    --verbose "${VERBOSE}" \
    --progress "${CURRENT_PROGRESS}" \
    --progress_hz "${PROGRESS_HZ}" \
    > "${OUT_DIR}/logs/stdout.txt"

  echo "✔ Run finished. Outputs in: ${OUT_DIR}"
fi
