#!/usr/bin/env bash
set -euo pipefail

# A defaulted THPRES entry is not a number anyone wrote down. It is the largest
# initial potential difference found anywhere along that region boundary, so
# arriving at it means looking at the whole boundary. A sector holds part of
# one, and working the entry out from a sector therefore gives a number that is
# too small -- on a real field, 0.01 bar where the parent had 0.81, which let
# flow through faces the parent held shut and left the sector drifting away
# from its parent by a third of a bar.
#
# The producer records its table and the consumer adopts it. This checks both
# halves: that the arrays reach the file, and that the consumer's own report of
# its thresholds is the parent's table rather than the one it would have
# derived.
#
# The decks are built so the two genuinely differ: the region boundary runs the
# full depth of the grid with the pressure jump growing a bar per layer, while
# the sector holds only the top three layers. Derived locally the entry comes
# out around 3 bar; the parent's is around 10.
#
# The same file is also built offline by make_flux from the producer's restart,
# and each file is given both to a consumer on the parent grid with everything
# outside the sector inactive and to one on a grid cut down to the sector. The
# offline route has to find the table in the restart, where it is stored in
# output units rather than SI, and the cut-out grid has to find the parent's
# cells at their place in a grid of its own; either going wrong changes the
# thresholds or the pressures the sector reaches.

flow_bin="$1"
producer_deck="$2"
consumer_deck="$3"
make_flux_bin="$4"
summary_bin="$5"
cutout_deck="$6"

producer_base="$(basename "$producer_deck" .DATA)"

cp "$producer_deck" "${producer_base}.DATA"
"$flow_bin" "${producer_base}.DATA" --output-dir="$(pwd)" > thpres_dump.log 2>&1
grep -q "End of simulation" thpres_dump.log

for array in FLXTHPR FLXEQLN; do
    if ! strings "${producer_base}.FLUX" | grep -qw "${array}"; then
        echo "check-flux-threshold-pressure: the producer wrote no ${array} array," \
             "so a reduced run cannot be told what the thresholds were" >&2
        exit 1
    fi
done

"$make_flux_bin" --parent="$PWD/${producer_base}" \
    --mapping-inline='box 4 1 1 7 1 3' \
    --output="$PWD/${producer_base}.MAKE.FLUX" > thpres_make_flux.log 2>&1 || {
    echo "check-flux-threshold-pressure: make_flux failed" >&2
    cat thpres_make_flux.log >&2
    exit 1
}

# The threshold table as the run itself reports it, one "from to value" triple
# per line. Comparing the report rather than the raw array keeps the check on
# the number the simulator will actually use.
thresholds() {
    grep -E '^ +[0-9]+ +[0-9]+ +[0-9.eE+-]+ +BARSA *$' "$1" \
        | awk '{print $1, $2, $3}'
}

last_value() {
    "$summary_bin" "$1.SMSPEC" "$2" 2>/dev/null | awk 'NF { value = $1 } END { print value }'
}

parent_thresholds="$(thresholds "${producer_base}.PRT")"

if [[ -z "$parent_thresholds" ]]; then
    echo "check-flux-threshold-pressure: the producer reported no threshold pressures," \
         "so there is nothing for this test to compare" >&2
    exit 1
fi

# Parent cells and where the cut-out grid has them, i shifted by the three
# columns cut away.
parent_cells=(4,1,1 6,1,2 7,1,3)
cutout_cells=(1,1,1 3,1,2 4,1,3)

run_consumer() {
    local label="$1" deck="$2" flux_file="$3"
    shift 3
    local cells=("$@")
    local base dir
    base="$(basename "$deck" .DATA)"
    dir="sector_${label}"

    rm -rf "$dir"
    mkdir -p "$dir"
    cp "$deck" "${dir}/${base}.DATA"
    # The consumer's USEFLUX names the producer's base, so the file takes it.
    cp "$flux_file" "${dir}/${producer_base}.FLUX"

    (cd "$dir" && "$flow_bin" "${base}.DATA" --output-dir=. > consumer.log 2>&1) || {
        echo "check-flux-threshold-pressure: ${label} consumer failed" >&2
        grep -iE "exception|error" "${dir}/consumer.log" >&2 || true
        exit 1
    }
    grep -q "End of simulation" "${dir}/consumer.log"

    if grep -q "have no matching grid" "${dir}/consumer.log" "${dir}/${base}.PRT"; then
        echo "check-flux-threshold-pressure: ${label} consumer could not place" \
             "the boundary on its grid" >&2
        grep -h "have no matching grid" "${dir}/${base}.PRT" >&2 || true
        exit 1
    fi

    # The adopted table is logged after the locally derived one, so take the
    # last report in the file.
    local sector_thresholds
    sector_thresholds="$(thresholds "${dir}/${base}.PRT" | tail -n "$(wc -l <<< "$parent_thresholds")")"
    if [[ "$sector_thresholds" != "$parent_thresholds" ]]; then
        echo "check-flux-threshold-pressure: the ${label} consumer did not adopt" \
             "the parent's threshold pressures" >&2
        echo "parent:" >&2
        echo "$parent_thresholds" >&2
        echo "sector:" >&2
        echo "$sector_thresholds" >&2
        exit 1
    fi

    local n parent sector
    for n in "${!parent_cells[@]}"; do
        parent="$(last_value "$producer_base" "BPR:${parent_cells[$n]}")"
        sector="$(last_value "${dir}/${base}" "BPR:${cells[$n]}")"
        awk -v p="$parent" -v s="$sector" -v c="${parent_cells[$n]}" -v l="$label" 'BEGIN {
            if (p == "" || s == "" || (p - s) ^ 2 > 1.0e-4) {
                printf "check-flux-threshold-pressure: BPR:%s is %s bar in the %s consumer, %s in the parent\n",
                       c, s, l, p > "/dev/stderr";
                exit 1;
            }
        }'
    done
}

run_consumer dumpflux "$consumer_deck" "${producer_base}.FLUX" "${parent_cells[@]}"
run_consumer dumpflux_cutout "$cutout_deck" "${producer_base}.FLUX" "${cutout_cells[@]}"
run_consumer make_flux "$consumer_deck" "${producer_base}.MAKE.FLUX" "${parent_cells[@]}"
run_consumer make_flux_cutout "$cutout_deck" "${producer_base}.MAKE.FLUX" "${cutout_cells[@]}"

echo "check-flux-threshold-pressure: ok"
