#!/usr/bin/env bash
set -euo pipefail

flow_bin="$1"
consumer_deck="$2"
np="${3:-1}"

# Ensure no matching FLUX artifacts exist for the consumer deck base name.
rm -f FLUX_USE_MISSING_SMOKE.FLUX
rm -f FLUX_USE_MISSING_SMOKE.FLUX????

cp "$consumer_deck" FLUX_USE_MISSING_SMOKE.DATA
if (( np > 1 )); then
    run_cmd=(mpirun -np "$np" "$flow_bin" FLUX_USE_MISSING_SMOKE.DATA)
else
    run_cmd=("$flow_bin" FLUX_USE_MISSING_SMOKE.DATA)
fi

if "${run_cmd[@]}" > useflux_missing_smoke.log 2>&1; then
    echo "Expected USEFLUX consumer run to fail when no matching FLUX file exists" >&2
    exit 1
fi

grep -q "USEFLUX requested but no FLUX file found" useflux_missing_smoke.log
