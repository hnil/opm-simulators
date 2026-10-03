#!/usr/bin/env bash
set -euo pipefail

# A fault with throw joins each cell on one side to whichever cells face it on
# the other, a layer or several up or down. Those connections are not in the
# deck: the grid works them out from the corner-point geometry. Where such a
# fault is the sector boundary, every one of them is a boundary face, and a
# FLUX file that leaves them out seals the fault in the reduced run however
# open it is in the full one.
#
# Here the fault has a full layer of throw, so there is no Cartesian
# connection across it at all and the sector meets the rest of the model only
# through those grid-generated NNCs. The producer outside draws its side down
# by two hundred bar over ninety days; a sector that has lost the connections
# stays at 300 bar, one that has them follows the parent.
#
# This is a pressure-mode sector, so the NNC faces have to be driven from the
# exterior pressure like any other boundary face rather than from a recorded
# mass.

flow_bin="$1"
summary_bin="$2"
inspect_bin="$3"
producer_deck="$4"
consumer_deck="$5"
make_flux_bin="$6"

parent_base="$(basename "$producer_deck" .DATA)"
sector_base="$(basename "$consumer_deck" .DATA)"

cp "$consumer_deck" "${sector_base}.DATA"
"$flow_bin" "${sector_base}.DATA" --output-dir="$(pwd)" > throw_use.log 2>&1
grep -q "End of simulation" throw_use.log

# The connections have to be in the file at all. inspect_flux names them NNC.
"$inspect_bin" "${parent_base}.FLUX" 3 1 2 2> inspect_throw.log > /dev/null || true
if ! grep -qE "^ +NNC.* exterior cell .*\(2,1,3\)" inspect_throw.log; then
    echo "check-flux-fault-throw: the producer recorded no NNC face from (3,1,2)" \
         "to (2,1,3) across the fault" >&2
    cat inspect_throw.log >&2
    exit 1
fi

# The offline route has to find them too, and with the parent's
# transmissibility. It has no grid to work them out from, only the parent's
# output, and its EGRID is the one place that lists them all.
cp "$producer_deck" "${parent_base}.DATA"
printf 'FLUXNUM\n  2*0 2*1\n  2*0 2*1\n  2*0 2*1\n  2*0 2*1 /\n' > throw_sector.grdecl
"$make_flux_bin" --parent="$parent_base" --fluxnum=throw_sector.grdecl --regions=1 \
                 --output=THROW_OFFLINE.FLUX > throw_make_flux.log 2>&1

"$inspect_bin" THROW_OFFLINE.FLUX 3 1 2 2> inspect_throw_offline.log > /dev/null || true
live_nnc="$(grep -E '^ +NNC' inspect_throw.log)"
offline_nnc="$(grep -E '^ +NNC' inspect_throw_offline.log || true)"
if [[ "$live_nnc" != "$offline_nnc" ]]; then
    echo "check-flux-fault-throw: make_flux does not describe the NNC faces as DUMPFLUX does" >&2
    echo "DUMPFLUX:"  >&2; echo "$live_nnc" >&2
    echo "make_flux:" >&2; echo "$offline_nnc" >&2
    exit 1
fi

value_at() {
    "$summary_bin" "$1.SMSPEC" TIME "$2" 2>/dev/null \
        | awk -v t="$3" 'NF == 2 && ($1 - t) ^ 2 < 1.0e-8 { print $2 }'
}

for cell in 3,1,2 4,1,4; do
    for day in 30 60 90; do
        parent="$(value_at "$parent_base" "BPR:${cell}" "$day")"
        sector="$(value_at "$sector_base" "BPR:${cell}" "$day")"

        awk -v p="$parent" -v s="$sector" -v c="$cell" -v d="$day" 'BEGIN {
            if (p == "" || s == "") {
                printf "check-flux-fault-throw: no BPR:%s at day %s\n", c, d > "/dev/stderr";
                exit 1;
            }
            # The parent must have moved, or a sealed sector would pass too.
            if (p > 250.0) {
                printf "check-flux-fault-throw: parent BPR:%s is still %s bar at day %s\n", c, p, d > "/dev/stderr";
                exit 1;
            }
            if ((s - p) ^ 2 > 0.25) {
                printf "check-flux-fault-throw: BPR:%s at day %s is %s bar in the sector, %s in the parent\n", c, d, s, p > "/dev/stderr";
                exit 1;
            }
        }'
    done
done

echo "check-flux-fault-throw: ok"
