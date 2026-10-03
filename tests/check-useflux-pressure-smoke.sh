#!/usr/bin/env bash
set -euo pipefail

flow_bin="$1"
consumer_deck="$2"

if [[ ! -f FLUX_DUMP_PRESSURE_SMOKE.FLUX ]]; then
    echo "Expected FLUX_DUMP_PRESSURE_SMOKE.FLUX in test result directory" >&2
    exit 1
fi

cp "$consumer_deck" FLUX_USE_PRESSURE_SMOKE.DATA
"$flow_bin" FLUX_USE_PRESSURE_SMOKE.DATA > useflux_pressure_smoke.log 2>&1

grep -q "End of simulation" useflux_pressure_smoke.log
