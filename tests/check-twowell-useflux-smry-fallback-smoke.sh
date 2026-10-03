#!/usr/bin/env bash
set -euo pipefail

# The parent run's summary vectors normally travel embedded in the .FLUX file.
# When they are absent - for instance when the .FLUX file was produced by an
# older run or by a tool that does not embed them - the sector run has to fall
# back to the parent run's own summary output.
#
# This checks that the fallback is taken, that the run says so, and that the
# field and group aggregates still reproduce the parent.

flow_bin="$1"
mutate_bin="$2"
compare_bin="$3"
producer_deck="$4"
consumer_deck="$5"

producer_base="$(basename "$producer_deck" .DATA)"
case_name="TWOWELL_USEFLUX_SMRYFALLBACK_1D"

rm -f "${case_name}".* smryfallback_dump.log smryfallback_use.log

# Parent run: writes the .FLUX file and its own summary output.
"$flow_bin" "$producer_deck" --output-dir="$(pwd)" > smryfallback_dump.log 2>&1

if [[ ! -f "${producer_base}.FLUX" ]]; then
    echo "Expected ${producer_base}.FLUX from the parent run" >&2
    exit 1
fi

if [[ ! -f "${producer_base}.SMSPEC" && ! -f "${producer_base}.ESMRY" ]]; then
    echo "Expected a parent summary file for ${producer_base}" >&2
    exit 1
fi

# Remove the embedded summary vectors, leaving the parent's summary file as the
# only remaining source.
"$mutate_bin" "${producer_base}.FLUX" "${producer_base}.FLUX.stripped" strip-summary
mv "${producer_base}.FLUX.stripped" "${producer_base}.FLUX"

# Point the consumer at the parent basename, so that both the .FLUX file and the
# fallback summary file resolve to the parent run rather than to the sector's
# own output.
cp "$consumer_deck" "${case_name}.DATA"
perl -0pi -e "s/USEFLUX\s*\n\s*1\*\s*\//USEFLUX\n  '${producer_base}' \//" "${case_name}.DATA"

if ! grep -q "'${producer_base}'" "${case_name}.DATA"; then
    echo "Failed to retarget USEFLUX in ${case_name}.DATA at ${producer_base}" >&2
    exit 1
fi

"$flow_bin" "${case_name}.DATA" --output-dir="$(pwd)" > smryfallback_use.log 2>&1

grep -q "End of simulation" smryfallback_use.log

# The fallback must be taken, and reported.
grep -q "USEFLUX: summary vectors loaded from parent summary file" smryfallback_use.log

# ... and the embedded source must NOT be claimed, since it was stripped.
if grep -q "USEFLUX: summary vectors loaded from flux file" smryfallback_use.log; then
    echo "Sector run used the embedded summary vectors even though they were stripped" >&2
    exit 1
fi

# Parity has to survive the fallback: the well inside the sector, and the field
# and group aggregates that depend on the well outside it.
#
# A parent summary file carries every vector, including the cumulatives, which
# makes this the case that catches seeding a cumulative instead of letting it
# accumulate -- that double counts production and shows up in FOPT/WOPT long
# before it shows up in any rate.
for keyword in WOPR:P1 FOPR FOPT GOPR:G1 GOPT:G1 WOPT:P1 WOPT:P2; do
    "$compare_bin" -d -t SMRY -k "$keyword" "$producer_base" "$case_name" 2e-2 2e-2
done
