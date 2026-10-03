#!/usr/bin/env bash
set -euo pipefail

# What a UDQ or an ACTIONX says it needs is a bare KEYWORD. A UDQ over WBHP
# reports "WBHP"; the wells it applies to are held separately, and the parent's
# summary holds WBHP:P1 and WBHP:P2. Comparing the two directly matches
# nothing, so every keyword naming an object is called missing and the tool
# refuses to run over what is in fact a complete summary.
#
# The same comparison decides which of the parent's vectors are worth carrying.
# Region and segment quantities suffer twice over: their objects cannot be
# enumerated from the schedule at all, so an expansion built from the schedule
# never contains them and they are dropped even when the parent reported them.
#
# This deck refers to WBHP from a UDQ and to GEFF and RPR from an ACTIONX, and
# its summary holds all three against real objects.

make_flux_bin="$1"
parent_deck="$2"
sector_grdecl="$3"

parent_base="$(basename "$parent_deck" .DATA)"

if [[ ! -f "${parent_base}.SMSPEC" ]]; then
    echo "expected ${parent_base}.SMSPEC from the parent run" >&2
    exit 1
fi

cp "$parent_deck" "${parent_base}.DATA"

if ! "$make_flux_bin" --parent="$parent_base" \
                      --fluxnum="$sector_grdecl" \
                      --regions=1 \
                      --output=UDQKEYS.FLUX > udqkeys.log 2>&1; then
    echo "check-make-flux-udq-keys: make_flux refused the parent summary" >&2
    cat udqkeys.log >&2
    exit 1
fi

# The keys live in a C032 array, space padded and several to a line, so split
# on whitespace before matching them whole.
carried_keys() {
    strings UDQKEYS.FLUX | tr -s ' ' '\n'
}

# Every vector the deck's expressions refer to has to be carried, named against
# the object it belongs to. RPR is the one that matters most: nothing but the
# parent's summary says which regions it reported on.
for key in "WBHP:P1" "WBHP:P2" "GEFF:G1" "RPR:1" "RPR:2"; do
    if ! carried_keys | grep -qx "$key"; then
        echo "check-make-flux-udq-keys: $key is referred to by the deck and present in the" \
             "parent summary, but was not carried into the FLUX file" >&2
        echo "carried:" >&2
        carried_keys | grep -E "^(WBHP|GEFF|RPR|FPR)" >&2 || true
        exit 1
    fi
done

# And nothing else: the point of selecting at all is to carry what a reduced
# run can use rather than the whole of the parent's summary.
if carried_keys | grep -qx "FPR"; then
    echo "check-make-flux-udq-keys: FPR is not referred to by the deck and should not" \
         "have been carried" >&2
    exit 1
fi
