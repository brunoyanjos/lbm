#!/usr/bin/env bash
set -euo pipefail

# =====================================================
# Exec defaults (override via env)
# =====================================================
: "${IO:=1}"
: "${RUN:=1}"
: "${CLEAN:=0}"
: "${DEBUG:=0}"
: "${RDC:=0}"
: "${PROGRESS:=1}"
: "${PROGRESS_HZ:=2}"
: "${WARMUP:=100}"
: "${REAL:=float}"
: "${STENCIL:=D2Q9}"
: "${ARCHES:=}"          # opcional: "86 89" etc
: "${BUILD_ROOT:=build}" # opcional
: "${OUT_ROOT:=runs}"    # opcional

# =====================================================
# Root discovery (script inside tools/)
# =====================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# =====================================================
# One of each case (all D2Q9)
# Pick Re values that make sense for each geometry in your solver.
# You can adjust freely.
# =====================================================
cases=(
  "annul         25"
  "channel       200"
  "jet           100"
)

echo "ROOT_DIR: ${ROOT_DIR}"
echo "STENCIL=${STENCIL} REAL=${REAL}"
echo "IO=${IO} WARMUP=${WARMUP} PROGRESS=${PROGRESS} PROGRESS_HZ=${PROGRESS_HZ}"
echo

for c in "${cases[@]}"; do
  set -- $c
  geom="$1"
  re="$2"

  ts="$(date +%Y%m%d_%H%M%S)"
  run_id="${ts}_${STENCIL}_${REAL}_${geom}_RE${re}"

  echo "================================================="
  echo "[CASE] GEOM=${geom}  STENCIL=${STENCIL}  REAL=${REAL}  RE=${re}"
  echo "       RUN_ID=${run_id}"
  echo "================================================="

  bash "${ROOT_DIR}/compile.sh" \
    --geom "${geom}" \
    --stencil "${STENCIL}" \
    --real "${REAL}" \
    --re "${re}" \
    --run "${RUN}" \
    --clean "${CLEAN}" \
    --debug "${DEBUG}" \
    --rdc "${RDC}" \
    --io "${IO}" \
    --warmup "${WARMUP}" \
    --progress "${PROGRESS}" \
    --progress_hz "${PROGRESS_HZ}" \
    --build_root "${BUILD_ROOT}" \
    --out_root "${OUT_ROOT}" \
    ${ARCHES:+--arches "${ARCHES}"} \
    --run_id "${run_id}"

  echo
done

echo "ALL CASES DONE."