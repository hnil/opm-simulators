#!/usr/bin/env bash
set -euo pipefail

# Sector whose boundary runs along a non-neighbour connection.
#
# An ordinary boundary face is imposed through the grid face it sits on, found
# from the interior cell and a direction. An NNC has no such face in a reduced
# run: the cell on the far side is outside the sector and absent from the grid,
# and the connection went with it. Those faces carry no direction either, so
# for a long time nothing ever looked them up and everything the parent
# recorded across one of them was dropped in silence.
#
# Here that connection carries about three quarters of the oil entering the
# sector, which makes the failure impossible to miss: the sector ends the run
# 46% below the parent on both WBHP:P1 and FPR, against the 4% that remains
# once the flux is applied. The tolerances below sit an order of magnitude
# inside that gap.
#
# The remaining 4% is the usual drift of a prescribed-rate boundary, not a
# defect in the NNC handling. Measured against how far the pressure actually
# moves over the run it is the same 6 to 7% the axis-aligned live-oil case
# shows, and the mass arriving over the boundary closes against the sector's
# own balance to within 0.05%.

flow_bin="$1"
compare_bin="$2"
producer_deck="$3"
consumer_deck="$4"

producer_base="$(basename "$producer_deck" .DATA)"
consumer_base="$(basename "$consumer_deck" .DATA)"

"$flow_bin" "$producer_deck" --output-dir="$(pwd)" > nnc_dump.log 2>&1
cp "${producer_base}.FLUX" "${consumer_base}.FLUX"
cp "$consumer_deck" "${consumer_base}.DATA"
"$flow_bin" "${consumer_base}.DATA" --output-dir="$(pwd)" > nnc_use.log 2>&1

grep -q "End of simulation" nnc_use.log

# Guard against the test passing for the wrong reason. The whole point is that
# the connection cannot be made in the reduced grid, so the flux across it has
# to arrive through the FLUX file. If a change ever left the connection intact
# here, the sector would get the flow through its own grid and this check would
# no longer be testing anything.
if ! grep -q "NNC between active and inactive cells" nnc_use.log; then
    echo "check-twowell-useflux-nnc-smoke: the boundary NNC was not severed in" \
         "the reduced grid, so this case no longer exercises the FLUX path" >&2
    exit 1
fi

# Pressure inside the sector. This is what collapses when the NNC flux is
# discarded, because it is most of the inflow.
for keyword in WBHP:P1 FPR; do
    "$compare_bin" -d -t SMRY -k "$keyword" "$producer_base" "$consumer_base" 8e-2 8e-2
done

# The well inside the sector holds its oil target either way, so the oil
# cumulative is only a sanity check. The gas is the telling one: starved of the
# boundary inflow the sector drops below bubble point and liberates gas it
# should not have, which shows up as roughly 19% too much produced gas.
"$compare_bin" -d -t SMRY -k FOPT "$producer_base" "$consumer_base" 1e-3 1e-3
"$compare_bin" -d -t SMRY -k FGPT "$producer_base" "$consumer_base" 5e-2 5e-2
