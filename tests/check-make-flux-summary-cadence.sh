#!/usr/bin/env bash
set -euo pipefail

# The boundary can only be described where the parent wrote a restart, but its
# summary is written at every time step, and a reduced run takes the rates of
# the wells outside the sector, and everything its UDQs and group controls
# read, from that. make_flux used to embed the summary at the restart cadence
# too, taking the value at each such report step -- which for a rate is the
# rate over the last time step before it, not over the report step.
#
# This parent restarts at every third report step while its summary has a
# sample a day, and its group target, and so the rates, move on every step.
# The embedded summary has to have every one of those samples, independently of
# the two boundary records; thinned on request, a sample's rate has to be the
# average over what it covers; and a reduced run stepping over several samples
# at once has to add up to the parent's production over its step. The check on
# all three is the cumulative of P2, which lies outside the sector and so is
# known to the reduced run only through the embedded summary.

flow_bin="$1"
make_flux_bin="$2"
inspect_bin="$3"
summary_bin="$4"
parent_deck="$5"
consumer_deck="$6"

parent_base="$(basename "$parent_deck" .DATA)"
consumer_base="$(basename "$consumer_deck" .DATA)"

cp "$parent_deck" "${parent_base}.DATA"
"$flow_bin" "${parent_base}.DATA" --output-dir="$(pwd)" \
            --solver-max-time-step-in-days=1 > cadence_parent.log 2>&1
grep -q "End of simulation" cadence_parent.log

printf 'FLUXNUM\n  2*0 3*1 6*0 /\n' > cadence_sector.grdecl

make_flux() {
    local output="$1"
    shift
    "$make_flux_bin" --parent="$parent_base" --fluxnum=cadence_sector.grdecl \
                     --regions=1 --output="$output" "$@" > "${output}.log" 2>&1 || {
        echo "check-make-flux-summary-cadence: make_flux failed for $output" >&2
        grep "make_flux:" "${output}.log" >&2 || true
        exit 1
    }
}

make_flux ALL.FLUX
make_flux THIN.FLUX --summary-min-interval=1000

parent_steps="$("$summary_bin" "${parent_base}.SMSPEC" TIME | awk '$1 ~ /^[0-9.]+$/' | wc -l)"
parent_total="$("$summary_bin" "${parent_base}.SMSPEC" WOPT:P2 | awk 'NF { v = $1 } END { print v }')"

samples() {
    "$inspect_bin" "$1" --summary=WOPR:P2 | awk '$1 !~ /^#/'
}

# 1. Every summary step, whatever the restart cadence. The boundary still has
#    its two records.
if ! grep -q "2 report steps" ALL.FLUX.log; then
    echo "check-make-flux-summary-cadence: expected two boundary records" >&2
    grep "Wrote" ALL.FLUX.log >&2 || true
    exit 1
fi

all_samples="$(samples ALL.FLUX | wc -l)"
if [[ "$all_samples" -ne "$parent_steps" ]]; then
    echo "check-make-flux-summary-cadence: ${all_samples} summary samples embedded," \
         "the parent wrote ${parent_steps}" >&2
    exit 1
fi

# 2. Thinned to the report steps, each rate is the average over its sample's
#    interval, so integrating the samples gives back the parent's cumulative.
thin_samples="$(samples THIN.FLUX | wc -l)"
if [[ "$thin_samples" -ge "$all_samples" ]]; then
    echo "check-make-flux-summary-cadence: --summary-min-interval kept ${thin_samples}" \
         "of ${all_samples} samples" >&2
    exit 1
fi

for file in ALL.FLUX THIN.FLUX; do
    integral="$(samples "$file" | awk '{ total += $2 * ($1 - t); t = $1 } END { print total }')"
    awk -v got="$integral" -v want="$parent_total" -v f="$file" 'BEGIN {
        if ((got - want) ^ 2 > (1.0e-4 * want) ^ 2) {
            printf "check-make-flux-summary-cadence: the rates in %s integrate to %s, the parent produced %s\n", f, got, want > "/dev/stderr";
            exit 1;
        }
    }'
done

# 3. The reduced run steps over several samples at a time. What it reports for
#    P2 is accumulated from the rates it was handed step by step, so it only
#    comes out right if each is the average over the step.
for file in ALL.FLUX THIN.FLUX; do
    dir="use_${file%.FLUX}"
    rm -rf "$dir"
    mkdir -p "$dir"
    cp "$consumer_deck" "$dir/${consumer_base}.DATA"
    cp "$file" "$dir/${consumer_base}.FLUX"

    (cd "$dir" && "$flow_bin" "${consumer_base}.DATA" --output-dir=. > use.log 2>&1) || {
        echo "check-make-flux-summary-cadence: reduced run failed with $file" >&2
        exit 1
    }

    total="$("$summary_bin" "$dir/${consumer_base}.SMSPEC" WOPT:P2 | awk 'NF { v = $1 } END { print v }')"
    awk -v got="$total" -v want="$parent_total" -v f="$file" 'BEGIN {
        if ((got - want) ^ 2 > (1.0e-4 * want) ^ 2) {
            printf "check-make-flux-summary-cadence: with %s the reduced run has WOPT:P2 = %s, the parent %s\n", f, got, want > "/dev/stderr";
            exit 1;
        }
    }'
done

echo "check-make-flux-summary-cadence: ok"
