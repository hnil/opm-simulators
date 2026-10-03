#!/usr/bin/env bash
set -euo pipefail

flow_bin="$1"
mutate_bin="$2"
consumer_deck="$3"
mutation="$4"
np="${5:-1}"

if [[ ! -f FLUX_DUMP_USEFLUX_FLUX_SMOKE.FLUX ]]; then
    echo "Expected FLUX_DUMP_USEFLUX_FLUX_SMOKE.FLUX in test result directory" >&2
    exit 1
fi

case "$mutation" in
    nan|inf|ninf)
        ;;
    *)
        echo "Unsupported mutation '$mutation' (expected nan|inf|ninf)" >&2
        exit 1
        ;;
esac

case_name="FLUX_USE_INVALID_RATE_${mutation^^}_SMOKE"
flux_file="${case_name}.FLUX"
deck_file="${case_name}.DATA"
log_file="useflux_invalid_rate_${mutation}_smoke.log"

rm -f "$flux_file" "$deck_file" "$log_file"
"$mutate_bin" FLUX_DUMP_USEFLUX_FLUX_SMOKE.FLUX "$flux_file" "$mutation"
cp "$consumer_deck" "$deck_file"

# The canonical consumer deck now uses the USEFLUX basename argument, so retarget
# the copied deck at the mutated FLUX file for this negative test.
perl -0pi -e "s/FLUX_DUMP_USEFLUX_FLUX_SMOKE/${case_name}/g" "$deck_file"

if (( np > 1 )); then
    run_cmd=(mpirun -np "$np" "$flow_bin" "$deck_file")
else
    run_cmd=("$flow_bin" "$deck_file")
fi

if "${run_cmd[@]}" > "$log_file" 2>&1; then
    echo "Expected USEFLUX consumer run to fail with $mutation FLUX rate" >&2
    exit 1
fi

grep -q "Invalid USEFLUX FLUX-mode component mass rate" "$log_file"
