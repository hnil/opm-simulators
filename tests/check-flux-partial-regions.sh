#!/usr/bin/env bash
set -euo pipefail

# A reduced run answers its ACTIONX conditions from the vectors embedded in the
# FLUX file. Which vectors those are is decided by keyword, and the expansion
# from a keyword to the objects it covers used to stop at wells and groups:
# nothing in a schedule says which regions of which region set a run reports
# on, so region, segment and node quantities were dropped on the floor.
#
# The condition here reads RPR 3 REC. Region 3 of the REC set lies wholly
# outside the sector, so the reduced run has no cells to compute it from and no
# way to arrive at the value other than being handed the parent's. Without it
# the run aborts with
#
#   Summary vector RPR__REC:3 is unknown
#
# The parent's own SUMMARY section never asks for region 3 either. opm-common
# computes it regardless, because an ACTIONX has to be answerable whether or
# not the deck asked for the vector in writing, and that is the value the FLUX
# file has to carry.

flow_bin="$1"
producer_deck="$2"
consumer_deck="$3"
num_ranks="$4"

producer_base="$(basename "$producer_deck" .DATA)"
consumer_base="$(basename "$consumer_deck" .DATA)"

cp "$producer_deck" "${producer_base}.DATA"
"$flow_bin" "${producer_base}.DATA" --output-dir="$(pwd)" > partial_regions_dump.log 2>&1
grep -q "End of simulation" partial_regions_dump.log

# The vector has to reach the file at all. Checking here as well as through the
# consumer says which half is at fault when this breaks.
if ! strings "${producer_base}.FLUX" | tr -s ' ' '\n' | grep -qx "RPR__REC:3"; then
    echo "check-flux-partial-regions: the producer did not embed RPR__REC:3," \
         "which its ACTIONX condition names outright" >&2
    echo "embedded region vectors:" >&2
    strings "${producer_base}.FLUX" | tr -s ' ' '\n' | grep -E "^RPR" >&2 || true
    exit 1
fi

# While we are here: the depth of each exterior cell has to be recorded too.
# A reduced run imposes the file's pressures at the boundary face, and they
# were measured at the centre of the cell on the far side; without the depth it
# cannot carry them from the one to the other, and a dipping boundary then
# drives flow through a sector that should be standing still.
if ! strings "${producer_base}.FLUX" | grep -qw "FLXEXDP"; then
    echo "check-flux-partial-regions: the producer wrote no FLXEXDP array, so a" \
         "reduced run has no way to tell where the recorded pressures were taken" >&2
    exit 1
fi

run_consumer() {
    local dir="$1"
    shift

    rm -rf "$dir"
    mkdir -p "$dir"
    cp "$consumer_deck" "$dir/${consumer_base}.DATA"
    cp "${producer_base}.FLUX" "$dir/${consumer_base}.FLUX"

    (cd "$dir" && "$@" "$flow_bin" "${consumer_base}.DATA" --output-dir=. \
        > consumer.log 2>&1) || {
        echo "check-flux-partial-regions: consumer failed on $dir" >&2
        grep -iE "exception|unknown" "$dir/consumer.log" >&2 || true
        return 1
    }

    grep -q "End of simulation" "$dir/consumer.log"
}

run_consumer serial

if [[ "$num_ranks" -gt 1 ]]; then
    # The same vector is wanted on every rank, so this is not a serial-only
    # concern even though nothing about it is rank dependent.
    run_consumer "parallel_${num_ranks}" mpirun -np "$num_ranks" --oversubscribe
fi

echo "check-flux-partial-regions: ok"
