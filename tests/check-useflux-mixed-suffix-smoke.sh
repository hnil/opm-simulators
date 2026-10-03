#!/usr/bin/env bash
set -euo pipefail

flow_bin="$1"
consumer_deck="$2"
np="${3:-1}"

if [[ ! -f FLUX_DUMP_USEFLUX_SUFFIX_SMOKE.FLUX ]]; then
    echo "Expected FLUX_DUMP_USEFLUX_SUFFIX_SMOKE.FLUX in test result directory" >&2
    exit 1
fi

# Stage mixed candidates:
# - valid lowest numeric suffix with real payload (.FLUX0002)
# - valid higher numeric suffix with intentionally bad payload (.FLUX0007)
# - invalid non-numeric suffix (.FLUXABCD)
# The run should still succeed by selecting the lowest valid numeric suffix.
rm -f FLUX_USE_MIXED_SUFFIX_SMOKE.FLUX
rm -f FLUX_USE_MIXED_SUFFIX_SMOKE.FLUX????
rm -f FLUX_USE_MIXED_SUFFIX_SMOKE.FLUXABCD
cp FLUX_DUMP_USEFLUX_SUFFIX_SMOKE.FLUX FLUX_USE_MIXED_SUFFIX_SMOKE.FLUX0002
printf "not-a-valid-flux-file\n" > FLUX_USE_MIXED_SUFFIX_SMOKE.FLUX0007
cp FLUX_DUMP_USEFLUX_SUFFIX_SMOKE.FLUX FLUX_USE_MIXED_SUFFIX_SMOKE.FLUXABCD

cp "$consumer_deck" FLUX_USE_MIXED_SUFFIX_SMOKE.DATA
if (( np > 1 )); then
    run_cmd=(mpirun -np "$np" "$flow_bin" FLUX_USE_MIXED_SUFFIX_SMOKE.DATA)
else
    run_cmd=("$flow_bin" FLUX_USE_MIXED_SUFFIX_SMOKE.DATA)
fi

"${run_cmd[@]}" > useflux_mixed_suffix_smoke.log 2>&1

grep -q "Total number of active cells: 100" useflux_mixed_suffix_smoke.log
grep -q "End of simulation" useflux_mixed_suffix_smoke.log
