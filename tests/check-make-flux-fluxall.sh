#!/usr/bin/env bash
set -euo pipefail

# A sector run still evaluates the parent deck's UDQ and ACTIONX expressions,
# but it no longer computes what they refer to, so the parent has to have
# reported those vectors. Working out which ones by hand means reading every
# expression in the deck and then naming the objects each could match; FLUXALL
# is how a deck says "report whatever my own expressions will ask for".
#
# This deck's SUMMARY section names nothing at all -- it holds FLUXALL and
# nothing else -- while a UDQ refers to WBHP and an ACTIONX to GEFF and RPR.
# So if FLUXALL works, make_flux has everything it needs; if it does not,
# make_flux refuses to run.

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
                      --output=FLUXALL.FLUX > fluxall.log 2>&1; then
    echo "check-make-flux-fluxall: make_flux refused a summary built by FLUXALL," \
         "which is the one thing FLUXALL exists to prevent" >&2
    cat fluxall.log >&2
    exit 1
fi

# The keys live in a C032 array, space padded and several to a line, so split
# on whitespace before matching them whole.
carried_keys() {
    strings FLUXALL.FLUX | tr -s ' ' '\n'
}

# FLUXALL has to reach the object as well as the keyword: naming every well for
# WBHP, every group for GEFF, and -- the one that takes the most work to get
# right -- every region of the FIPNUM set for RPR.
#
# The two user defined quantities have to be carried in their own right as
# well. WUFOUR because the ACTIONX compares against it, and WUTWICE because
# WUFOUR is defined in terms of it: a reduced run can only recompute a UDQ
# whose inputs it still has.
for key in "WBHP:P1" "WBHP:P2" "GEFF:G1" "RPR:1" "RPR:2" "WUTWICE:P1" "WUTWICE:P2" "WUFOUR:P1" "WUFOUR:P2"; do
    if ! carried_keys | grep -qx "$key"; then
        echo "check-make-flux-fluxall: $key is referred to by the deck's expressions," \
             "but FLUXALL did not put it in the parent summary" >&2
        echo "carried:" >&2
        carried_keys | grep -E "^(WBHP|GEFF|RPR|WUTWICE|WUFOUR)" >&2 || true
        exit 1
    fi
done

echo "check-make-flux-fluxall: ok"
