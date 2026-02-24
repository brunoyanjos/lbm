#!/usr/bin/env bash
set -euo pipefail

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

# Se quiser evitar que um caso sobrescreva build do outro, seu compile.sh já separa por STENCIL/REAL.
# Só não misture REAL no meio sem querer.

cases=(
  "D2Q9  3200"
  "D2Q9 10000"
  "D2V17 3200"
  "D2V17 10000"
)

for c in "${cases[@]}"; do
  set -- $c
  stencil="$1"
  re="$2"

  ts="$(date +%Y%m%d_%H%M%S)"
  export STENCIL="$stencil"
  export RE="$re"
  export RUN_ID="${ts}_${stencil}_RE${re}"

  echo "================================================="
  echo "[CASE] STENCIL=${STENCIL}  RE=${RE}  RUN_ID=${RUN_ID}"
  echo "================================================="

  ./compile.sh

  echo
done

echo "ALL CASES DONE."