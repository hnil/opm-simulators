#!/usr/bin/env bash
set -euo pipefail

# A FLUX file must not depend on how many ranks produced it.
#
# The dumpers used to live on the IO rank alone, so it reported on the cells of
# its own partition and wrote zeros for the rest: zero transmissibility, zero
# exterior pressure, zero relative permeability, and region sums covering only
# the part of the model that happened to land on rank 0. None of that is
# visible from inside a parallel run, and the resulting file is not obviously
# broken either -- a consumer applies the zero transmissibilities and simply
# runs a sealed sector.
#
# Comparing against the serial file catches all of it at once and needs no
# tolerance: the same deck on the same build gives the same boundary, so the
# two files should agree exactly.

flow_bin="$1"
compare_bin="$2"
producer_deck="$3"
num_ranks="$4"

parent_base="$(basename "$producer_deck" .DATA)"

cp "$producer_deck" "${parent_base}.DATA"

rm -f "${parent_base}.FLUX"
"$flow_bin" "${parent_base}.DATA" --output-dir="$(pwd)" > parallel_serial_run.log 2>&1
grep -q "End of simulation" parallel_serial_run.log

if [[ ! -f "${parent_base}.FLUX" ]]; then
    echo "check-flux-parallel-consistency: the serial run wrote no FLUX file" >&2
    exit 1
fi
mv "${parent_base}.FLUX" serial.FLUX

rm -f "${parent_base}.FLUX"
mpirun -np "$num_ranks" --oversubscribe "$flow_bin" "${parent_base}.DATA" \
       --output-dir="$(pwd)" > parallel_mpi_run.log 2>&1

if [[ ! -f "${parent_base}.FLUX" ]]; then
    echo "check-flux-parallel-consistency: the ${num_ranks}-rank run wrote no FLUX file" >&2
    exit 1
fi

# Guard against the comparison passing because the file is empty. Both must
# describe the same non-trivial boundary, so insist there is one.
faces="$(strings serial.FLUX | grep -c FACECELL || true)"
if [[ "$faces" -lt 1 ]]; then
    echo "check-flux-parallel-consistency: the serial file describes no boundary faces," \
         "so this comparison could not detect a difference" >&2
    exit 1
fi

if ! "$compare_bin" serial.FLUX "${parent_base}.FLUX"; then
    echo "check-flux-parallel-consistency: the ${num_ranks}-rank file differs from the serial one." \
         "A boundary record is gathered from every rank, so each has to reach" \
         "the collection and contribute the cells it owns." >&2
    exit 1
fi
