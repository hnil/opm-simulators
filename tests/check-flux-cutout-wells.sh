#!/usr/bin/env bash
set -euo pipefail

# A sector cut out of the full model as a grid of its own usually keeps only
# its own wells in its SCHEDULE section. The wells outside it still produce and
# inject in the parent, still count towards every group and field total, still
# take their share of a group target, and a UDQ over them still has a value.
# The FLUX file carries their rates and says which group each is in and whether
# it produces or injects, and the reduced run adds them to its schedule without
# connections, so that all of that carries on as in the parent.
#
# This parent has P1 in the sector, and P2 in the same group outside it with an
# efficiency factor, and I1 outside it as the only well of a group the sector
# deck does not have. The sector deck lists P1 alone. G1's target is shared
# between P1 and P2, so P1's own rate is only right if P2 takes its share. The
# UDQs read P2's rate, bottom-hole pressure and block average pressure, sum
# over all wells, and read the field pressure, the pressure of a FIPNUM region
# that straddles the sector's edge and the oil rate of one wholly outside it,
# which the sector deck's own FIPNUM does not even have. None of those can be
# computed from the sector, and each has to come out as in the parent.

flow_bin="$1"
make_flux_bin="$2"
summary_bin="$3"
parent_deck="$4"
sector_deck="$5"

parent_base="$(basename "$parent_deck" .DATA)"
sector_base="$(basename "$sector_deck" .DATA)"

rm -rf sector
mkdir -p sector

cp "$parent_deck" "${parent_base}.DATA"

"$make_flux_bin" --parent="$PWD/${parent_base}" \
                 --mapping-inline='box 3 1 1 5 1 1' \
                 --output="$PWD/sector/${parent_base}.FLUX" > cutout_make_flux.log 2>&1 || {
    echo "check-flux-cutout-wells: make_flux failed" >&2
    cat cutout_make_flux.log >&2
    exit 1
}

cp "$sector_deck" "sector/${sector_base}.DATA"
(cd sector && "$flow_bin" "${sector_base}.DATA" --output-dir=. > sector.log 2>&1) || {
    echo "check-flux-cutout-wells: the sector run failed" >&2
    grep -iE "error|problem" -A4 sector/sector.log >&2 || true
    exit 1
}

if ! grep -q "USEFLUX: 2 well(s) of the parent run are not in this deck" "sector/${sector_base}.PRT"; then
    echo "check-flux-cutout-wells: the sector run did not add P2 and I1" >&2
    grep "USEFLUX" "sector/${sector_base}.PRT" >&2 || true
    exit 1
fi

vectors=(FOPR FOPT FWIR FWIT FVPR FVPT GOPR:G1 GOPT:G1 GWIR:WI
         WOPR:P1 WOPR:P2 WWIR:I1 FU_P2 FU_WSUM FU_BHP2
         FU_BP4P2 FU_RPR1 FU_ROPR2 FU_FPR)

values() {
    "$summary_bin" -r -n "$1.SMSPEC" "${vectors[@]}"
}

parent_values="$(values "$parent_base")"
sector_values="$(values "sector/${sector_base}")"

paste <(echo "$parent_values") <(echo "$sector_values") \
    | awk -v n="${#vectors[@]}" -v names="${vectors[*]}" '
        BEGIN { split(names, name, " "); bad = 0; rows = 0 }
        NF == 0 { next }
        {
            ++rows;
            if (NF != 2 * n) { printf "check-flux-cutout-wells: row %d has %d columns\n", rows, NF > "/dev/stderr"; bad = 1; exit }
            for (i = 1; i <= n; ++i) {
                p = $i; s = $(i + n);
                scale = (p < 0 ? -p : p); if (scale < 1) scale = 1;
                d = p - s; if (d < 0) d = -d;
                if (d > 1.0e-4 * scale) {
                    printf "check-flux-cutout-wells: %s at report step %d is %s in the sector, %s in the parent\n",
                           name[i], rows, s, p > "/dev/stderr";
                    bad = 1;
                }
            }
        }
        END { if (rows == 0) { print "check-flux-cutout-wells: no report steps to compare" > "/dev/stderr"; bad = 1 } exit bad }'

echo "check-flux-cutout-wells: ok"
