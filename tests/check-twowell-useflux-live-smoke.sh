#!/usr/bin/env bash
set -euo pipefail

# Live-oil / wet-gas sector case. One producer inside the FLUX region and one
# outside, no group control, so the only thing keeping the sector on track is
# the boundary data itself.
#
# This is the regression guard for the component split of the boundary mass
# rate. The oil phase here carries a sixth of its mass as dissolved gas and the
# gas phase carries vaporised oil, so a producer that writes PHASE masses and
# calls them components puts a large part of the inflow into the wrong
# conservation equation. That showed up as a 16% error in WBHP:P1, against the
# 0.5% that remains once the split is done properly.
#
# The producer deck also leaves FLORES out of RPTRST on purpose, so a
# regression in the automatic FLORES activation would zero the boundary rates
# and strand the sector as a closed region, which this check would catch too.

flow_bin="$1"
compare_bin="$2"
producer_deck="$3"
consumer_deck="$4"

producer_base="$(basename "$producer_deck" .DATA)"
consumer_base="$(basename "$consumer_deck" .DATA)"

"$flow_bin" "$producer_deck" --output-dir="$(pwd)" > live_dump.log 2>&1
cp "${producer_base}.FLUX" "${consumer_base}.FLUX"
cp "$consumer_deck" "${consumer_base}.DATA"
"$flow_bin" "${consumer_base}.DATA" --output-dir="$(pwd)" > live_use.log 2>&1

grep -q "End of simulation" live_use.log

# The producer must not have fallen back to a closed boundary.
if grep -q "no FLORES data is available" live_dump.log; then
    echo "check-twowell-useflux-live-smoke: producer did not enable FLORES" >&2
    exit 1
fi

# Local pressure inside the sector. This is the quantity the component split
# moves: booking dissolved gas as oil leans out the oil, shrinks its formation
# volume factor and drops the pressure.
for keyword in WBHP:P1 FPR; do
    "$compare_bin" -d -t SMRY -k "$keyword" "$producer_base" "$consumer_base" 2e-2 2e-2
done

# The well inside the sector is on a rate target it can always meet, so its
# cumulative oil must match the parent exactly. The produced gas is the
# sensitive one: it is what goes missing when gas mass is booked as oil.
"$compare_bin" -d -t SMRY -k FOPT "$producer_base" "$consumer_base" 1e-3 1e-3
"$compare_bin" -d -t SMRY -k FGPT "$producer_base" "$consumer_base" 1e-2 1e-2
