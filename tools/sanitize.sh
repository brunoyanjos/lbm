#!/usr/bin/env bash
set -euo pipefail

: "${STENCIL:=D2V17}"
: "${REAL:=float}"
: "${GEOM:=square_cavity}"
: "${BUILD_ROOT:=build}"
: "${EXEC_NAME:=sim}"

: "${OUT_ROOT:=runs}"
: "${RUN_ID:=}"
: "${ARGS:=--io 0 --progress 0 --verbose 0}"

CFG_TAG="${STENCIL}_${REAL}_${GEOM}"
BIN="${BUILD_ROOT}/${CFG_TAG}/${EXEC_NAME}"

[[ -x "${BIN}" ]] || { echo "Binary not found/executable: ${BIN}" >&2; exit 1; }

if [[ -z "${RUN_ID}" ]]; then
  RUN_ID="$(date +%Y%m%d_%H%M%S)_sanitize_${STENCIL}_${REAL}_${GEOM}"
fi

OUT_DIR="${OUT_ROOT}/${RUN_ID}"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

common=(
  --print-limit 1000
  --check-exit-code no
  --check-api-memory-access yes
  --track-unused-memory
)

run_tool () {
  local tool="$1"
  local log="${LOG_DIR}/${tool}.txt"

  echo "[SANITIZE] tool=${tool}"
  echo "  bin:  ${BIN}"
  echo "  args: ${ARGS}"

  compute-sanitizer "${common[@]}" --tool="${tool}" "${BIN}" ${ARGS} > "${log}" 2>&1
  local code=$?

  echo "  exit=${code}"
  echo "  -> log: ${log}"
  return ${code}
}

run_tool memcheck
run_tool initcheck

set +e
run_tool racecheck
echo "[WARN] racecheck exit=$?"
run_tool synccheck
echo "[WARN] synccheck exit=$?"
set -e

echo "✔ Sanitizers finished. Logs in: ${LOG_DIR}"