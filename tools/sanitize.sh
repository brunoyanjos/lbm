#!/usr/bin/env bash
set -euo pipefail

# -------------------------------
# Config (overridable)
# -------------------------------
: "${STENCIL:=D2V17}"     # D2Q9 | D2V17
: "${REAL:=float}"        # float | double
: "${BUILD_ROOT:=build}"
: "${EXEC_NAME:=sim}"

: "${OUT_ROOT:=runs}"
: "${RUN_ID:=}"           # se vazio, usa timestamp
: "${ARGS:=--io 0 --progress 0 --verbose 0}"  # defaults p/ log limpo

# -------------------------------
# Paths
# -------------------------------
CFG_TAG="${STENCIL}_${REAL}"
BIN="${BUILD_ROOT}/${CFG_TAG}/${EXEC_NAME}"

[[ -x "${BIN}" ]] || { echo "Binary not found/executable: ${BIN}" >&2; exit 1; }

if [[ -z "${RUN_ID}" ]]; then
  RUN_ID="$(date +%Y%m%d_%H%M%S)"
fi

OUT_DIR="${OUT_ROOT}/${RUN_ID}"
LOG_DIR="${OUT_DIR}/logs"
mkdir -p "${LOG_DIR}"

# -------------------------------
# Sanitizer common flags
# -------------------------------
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
  compute-sanitizer "${common[@]}" --tool="${tool}" "${BIN}" ${ARGS} \
    > "${log}" 2>&1
  echo "  -> log: ${log}"
}

# -------------------------------
# Run (sequential)
# -------------------------------
run_tool memcheck
run_tool initcheck

# No WSL/WDDM, racecheck/synccheck podem falhar. Rodar best-effort.
set +e
run_tool racecheck
echo "[WARN] racecheck exit=$?"
run_tool synccheck
echo "[WARN] synccheck exit=$?"
set -e

echo "✔ Sanitizers finished. Logs in: ${LOG_DIR}"