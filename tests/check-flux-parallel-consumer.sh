#!/usr/bin/env bash
set -euo pipefail

# A USEFLUX run must give the same answer however many ranks it uses.
#
# The boundary faces in a FLUX file are described against the region's own cell
# numbering, which a serial reduced run happens to reproduce because it
# activates the region's cells in the same order. A parallel rank does not: it
# holds an arbitrary subset of the region under its own numbering, so the
# mapping has to be looked up in the grid rather than assumed. Getting that
# wrong took the interior cell of every face from the wrong place, which the
# consumer noticed only because the indices ran off the end of its arrays.
#
# Comparing against the serial answer covers the case where they would not have
# run off the end and the run would simply have imposed the boundary on the
# wrong cells.

flow_bin="$1"
summary_bin="$2"
producer_deck="$3"
consumer_deck="$4"
num_ranks="$5"
shift 5
keys=("$@")

producer_base="$(basename "$producer_deck" .DATA)"
consumer_base="$(basename "$consumer_deck" .DATA)"

# One producer run supplies the boundary for both consumer runs.
cp "$producer_deck" "${producer_base}.DATA"
"$flow_bin" "${producer_base}.DATA" --output-dir="$(pwd)" > consumer_dump.log 2>&1
grep -q "End of simulation" consumer_dump.log

run_consumer() {
    local dir="$1"
    shift

    rm -rf "$dir"
    mkdir -p "$dir"
    cp "$consumer_deck" "$dir/${consumer_base}.DATA"
    cp "${producer_base}.FLUX" "$dir/${consumer_base}.FLUX"

    ( cd "$dir" && "$@" "$flow_bin" "${consumer_base}.DATA" --output-dir=. ) \
        > "consumer_${dir}.log" 2>&1

    grep -q "End of simulation" "consumer_${dir}.log"
}

run_consumer serial
run_consumer parallel mpirun -np "$num_ranks" --oversubscribe

last_value() {
    "$summary_bin" "$1/${consumer_base}.SMSPEC" "$2" 2>/dev/null \
        | awk 'NF { value = $1 } END { print value }'
}

failed=0
for key in "${keys[@]}"; do
    serial_value="$(last_value serial "$key")"
    parallel_value="$(last_value parallel "$key")"

    if [[ -z "$serial_value" || -z "$parallel_value" ]]; then
        echo "check-flux-parallel-consumer: $key is missing from one of the runs" >&2
        failed=1
        continue
    fi

    # Loose enough for the linear solver to order its work differently, far
    # tighter than imposing the boundary on the wrong cells would be.
    if ! awk -v s="$serial_value" -v p="$parallel_value" -v k="$key" 'BEGIN {
             scale = (s < 0 ? -s : s);
             if (scale < 1.0) { scale = 1.0; }
             diff = (s - p < 0 ? p - s : s - p);
             if (diff / scale > 1.0e-4) {
                 printf "  %s: serial %s vs %d-rank %s\n", k, s, ranks, p > "/dev/stderr";
                 exit 1;
             }
         }' ranks="$num_ranks"; then
        failed=1
    fi
done

if [[ "$failed" -ne 0 ]]; then
    echo "check-flux-parallel-consumer: the ${num_ranks}-rank run disagrees with the serial one." \
         "Each rank holds part of the region, so the region's cells have to be" \
         "resolved against that rank's own grid." >&2
    exit 1
fi
