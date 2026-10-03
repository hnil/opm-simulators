#!/usr/bin/env bash
set -euo pipefail

# make_flux has to pair each restart step with the summary step of the same
# time. Doing that by position assumes the parent wrote a restart at every
# report step, which a deck is under no obligation to do.
#
# This parent has four report steps and restarts at 2 and 4 only, so the two
# sequences differ in length -- which used to be refused outright -- and the
# same position means a different time in each. Reading the summary by position
# would put report step 1 against restart step 2, and report step 2 against
# restart step 4.

make_flux_bin="$1"
inspect_flux_bin="$2"
parent_deck="$3"
sector_grdecl="$4"

parent_base="$(basename "$parent_deck" .DATA)"

if [[ ! -f "${parent_base}.SMSPEC" ]]; then
    echo "expected ${parent_base}.SMSPEC from the parent run" >&2
    exit 1
fi

cp "$parent_deck" "${parent_base}.DATA"

if ! "$make_flux_bin" --parent="$parent_base" \
                      --fluxnum="$sector_grdecl" \
                      --regions=1 \
                      --output=SPARSE.FLUX > sparse.log 2>&1; then
    echo "check-make-flux-sparse-restart: make_flux refused a parent that wrote" \
         "restarts at every second report step" >&2
    cat sparse.log >&2
    exit 1
fi

# One record per restart step, not per report step.
if ! grep -q "2 report steps" sparse.log; then
    echo "check-make-flux-sparse-restart: expected two records, one per restart step" >&2
    cat sparse.log >&2
    exit 1
fi

# The times decide it. Reading the summary by position would give 1 and 2;
# pairing by report step number gives 2 and 4, which is when the parent
# actually wrote those restarts. The initial state comes first, at 0, as it does
# in a DUMPFLUX file.
times="$("$inspect_flux_bin" SPARSE.FLUX 1 1 2 | awk '$1 ~ /^[0-9]+$/ { print $1 }' | tr '\n' ' ')"

if [[ "$times" != "0 2 4 " ]]; then
    echo "check-make-flux-sparse-restart: records are at times [$times], expected [0 2 4 ]" >&2
    echo "A pairing by position would give [0 1 2 ]." >&2
    exit 1
fi

echo "check-make-flux-sparse-restart: ok"
