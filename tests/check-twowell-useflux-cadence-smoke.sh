#!/usr/bin/env bash
set -euo pipefail

# Sector boundary data is written per time step by default, and once per report
# step when --flux-boundary-report-steps-only is given.
#
# This checks that the two cadences produce the expected number of records, that
# both still reproduce the parent, and that the throttle bounds the record count
# without dropping the report step boundaries.

flow_bin="$1"
compare_bin="$2"
producer_deck="$3"
consumer_deck="$4"

producer_base="$(basename "$producer_deck" .DATA)"
consumer_base="$(basename "$consumer_deck" .DATA)"

cp "$consumer_deck" "${consumer_base}.DATA"

# Number of boundary records a consumer run reports having read.
record_count() {
    sed -n 's/.*USEFLUX: \([0-9][0-9]*\) boundary record(s).*/\1/p' "$1" | head -1
}

run_pair() {
    local tag="$1"
    shift

    rm -f "${producer_base}.FLUX" "${consumer_base}.FLUX"
    "$flow_bin" "$producer_deck" --output-dir="$(pwd)" "$@" > "cadence_dump_${tag}.log" 2>&1
    cp "${producer_base}.FLUX" "${consumer_base}.FLUX"
    "$flow_bin" "${consumer_base}.DATA" --output-dir="$(pwd)" > "cadence_use_${tag}.log" 2>&1

    grep -q "End of simulation" "cadence_use_${tag}.log"
}

run_pair default
n_default="$(record_count cadence_use_default.log)"

run_pair reportonly --flux-boundary-report-steps-only=true
n_report="$(record_count cadence_use_reportonly.log)"

run_pair dense --flux-boundary-min-interval-between-samples=0
n_dense="$(record_count cadence_use_dense.log)"

echo "boundary records: default=${n_default} reportonly=${n_report} dense=${n_dense}"

if [[ -z "$n_default" || -z "$n_report" || -z "$n_dense" ]]; then
    echo "Could not determine the boundary record count from the consumer logs" >&2
    exit 1
fi

# Per time step must produce at least as many records as per report step, and
# for this deck strictly more, since the first report step is split into
# several time steps.
if (( n_default <= n_report )); then
    echo "Expected more boundary records per time step (${n_default}) than per report step (${n_report})" >&2
    exit 1
fi

# Removing the throttle cannot reduce the number of records.
if (( n_dense < n_default )); then
    echo "Expected at least as many records without a throttle (${n_dense}) as with one (${n_default})" >&2
    exit 1
fi

# The report-step cadence must be reported as such, and the default must not be.
grep -q "per report step" cadence_use_reportonly.log
grep -q "per time step" cadence_use_default.log

# Both cadences have to reproduce the parent. The well inside the sector and the
# field aggregates are the quantities the sector is supposed to get right.
for tag in default reportonly; do
    case "$tag" in
        default)    opts=() ;;
        reportonly) opts=(--flux-boundary-report-steps-only=true) ;;
    esac

    run_pair "$tag" "${opts[@]}"

    for keyword in WOPR:P1 WBHP:P1 FOPR FOPT; do
        "$compare_bin" -d -t SMRY -k "$keyword" "$producer_base" "$consumer_base" 2e-2 2e-2
    done
done
