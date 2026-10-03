#!/usr/bin/env bash
set -euo pipefail

flow_bin="$1"
consumer_deck="$2"

if [[ ! -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE.FLUX ]]; then
    echo "Expected FLUX_DUMP_USEFLUX_SUFFIX_SMOKE.FLUX in test result directory" >&2
    exit 1
fi

# Exercise generalized .FLUXdddd fallback by staging only a .FLUX0002 file.
cp FLUX_DUMP_USEFLUX_SUFFIX_SMOKE.FLUX FLUX_USE_SUFFIX_SMOKE.FLUX0002
cp "$consumer_deck" FLUX_USE_SUFFIX_SMOKE.DATA
"$flow_bin" FLUX_USE_SUFFIX_SMOKE.DATA > useflux_suffix2_smoke.log 2>&1

grep -q "Total number of active cells: 100" useflux_suffix2_smoke.log
grep -q "End of simulation" useflux_suffix2_smoke.log
