#!/usr/bin/env bash
set -euo pipefail

# A boundary record holds the parent's state at the end of the interval it
# covers. Held over the whole interval, that state reaches the sector early,
# and by as much as the interval is long: a file built from a parent that wrote
# restarts once a year imposes each year's end-of-year pressure from January.
# On Drogon that showed as a sector well's bottom-hole pressure 3 bar off for
# the two years the parent's prediction wrote no restart.
#
# This parent restarts at every third report step, so the file has two records
# of fifteen days each after the initial state, and the pressure falls about 9
# bar over each. By default a reduced run interpolates between records, from
# the initial state over the first one, and has to follow the parent to within
# a bar or so at every step. Asked not to, it holds each record's end state and
# is 7 bar or more off early in each interval, which is what the default has to
# improve on.

flow_bin="$1"
make_flux_bin="$2"
summary_bin="$3"
parent_deck="$4"
consumer_deck="$5"

parent_base="$(basename "$parent_deck" .DATA)"
consumer_base="$(basename "$consumer_deck" .DATA)"

cp "$parent_deck" "${parent_base}.DATA"
"$flow_bin" "${parent_base}.DATA" --output-dir="$(pwd)" \
            --solver-max-time-step-in-days=1 > interp_parent.log 2>&1
grep -q "End of simulation" interp_parent.log

printf 'FLUXNUM\n  2*0 3*1 6*0 /\n' > interp_sector.grdecl
"$make_flux_bin" --parent="$parent_base" --fluxnum=interp_sector.grdecl \
                 --regions=1 --output="${consumer_base}.FLUX" > interp_make_flux.log 2>&1 || {
    echo "check-flux-boundary-interpolation: make_flux failed" >&2
    cat interp_make_flux.log >&2
    exit 1
}

# WBHP:P1 at the given days, one "day value" pair per line. P1 is inside the
# sector, so its pressure is the sector's.
bhp_at() {
    "$summary_bin" "$1.SMSPEC" TIME WBHP:P1 2>/dev/null \
        | awk -v days="5 10 20 25" 'BEGIN { n = split(days, d, " "); for (i = 1; i <= n; ++i) want[d[i]] = 1 }
                                     $1 ~ /^[0-9.]+$/ && (($1 + 0) in want) { print $1 + 0, $2 }'
}

run_consumer() {
    local dir="$1"
    shift
    rm -rf "$dir"
    mkdir -p "$dir"
    cp "$consumer_deck" "${dir}/${consumer_base}.DATA"
    cp "${consumer_base}.FLUX" "${dir}/"
    (cd "$dir" && "$flow_bin" "${consumer_base}.DATA" --output-dir=. "$@" > use.log 2>&1) || {
        echo "check-flux-boundary-interpolation: reduced run in ${dir} failed" >&2
        exit 1
    }
}

run_consumer default
run_consumer held --flux-boundary-interpolate=false

parent="$(bhp_at "$parent_base")"

# Largest deviation from the parent over the chosen days, or -1 when the two
# runs do not both report all of them.
worst() {
    paste -d ' ' <(echo "$parent") <(bhp_at "$1") \
        | awk 'BEGIN { w = -1 } NF != 4 || $1 != $3 { w = -1; exit }
               { d = $2 - $4; if (d < 0) d = -d; if (d > w) w = d }
               END { print (NR == 4) ? w : -1 }'
}

default_worst="$(worst "default/${consumer_base}")"
held_worst="$(worst "held/${consumer_base}")"

awk -v got="$default_worst" 'BEGIN { exit !(got >= 0 && got < 1.5) }' || {
    echo "check-flux-boundary-interpolation: interpolating, WBHP:P1 is up to" \
         "${default_worst} bar off the parent" >&2
    exit 1
}

awk -v got="$held_worst" 'BEGIN { exit !(got > 5.0) }' || {
    echo "check-flux-boundary-interpolation: holding each record, WBHP:P1 is at" \
         "most ${held_worst} bar off the parent, so this test no longer shows" \
         "what interpolating is for" >&2
    exit 1
}

echo "check-flux-boundary-interpolation: ok (worst ${default_worst} bar interpolated," \
     "${held_worst} bar held)"
