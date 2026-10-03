#!/usr/bin/env bash
set -euo pipefail

flow_bin="$1"
consumer_deck="$2"

if [[ ! -f FLUX_DUMP_USEFLUX_BOTH_SMOKE.FLUX ]]; then
    echo "Expected FLUX_DUMP_USEFLUX_BOTH_SMOKE.FLUX in test result directory" >&2
    exit 1
fi

cp FLUX_DUMP_USEFLUX_BOTH_SMOKE.FLUX FLUX_USE_BOTH_SMOKE.FLUX
cp "$consumer_deck" FLUX_USE_BOTH_SMOKE.DATA
"$flow_bin" FLUX_USE_BOTH_SMOKE.DATA > useflux_both_smoke.log 2>&1

grep -q "Total number of active cells: 100" useflux_both_smoke.log
grep -q "End of simulation" useflux_both_smoke.log
