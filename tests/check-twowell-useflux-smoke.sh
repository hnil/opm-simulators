#!/usr/bin/env bash
set -euo pipefail

# Two-well sector case: one producer inside the FLUX region and one outside,
# sharing a group target. The sector run has to drop the well outside the
# region and still count what that well produced in the parent run against the
# group, handing the well inside only what is left. The target is driven by a
# UDQ in some of the decks and is a plain number in others; with a plain number
# nothing but the group control itself can take the absent well into account.

flow_bin="$1"
compare_bin="$2"
producer_deck="$3"
consumer_deck="$4"

producer_base="$(basename "$producer_deck" .DATA)"
consumer_base="$(basename "$consumer_deck" .DATA)"

"$flow_bin" "$producer_deck" --output-dir="$(pwd)" > twowell_dump.log 2>&1
cp "${producer_base}.FLUX" "${consumer_base}.FLUX"
cp "$consumer_deck" "${consumer_base}.DATA"
"$flow_bin" "${consumer_base}.DATA" --output-dir="$(pwd)" > twowell_use.log 2>&1

grep -q "End of simulation" twowell_use.log

# The sector must report where it read the parent's summary vectors from.
grep -q "USEFLUX: summary vectors loaded from flux file" twowell_use.log

# Compare on report steps only. The FLUX payload carries the parent's summary
# values once per report step, while the UDQ in this deck advances on every
# time step, so the sector can only reproduce the parent at report step
# boundaries - which is where the exchanged boundary data is defined.
for keyword in WOPR:P1 WBHP:P1; do
    "$compare_bin" -d -t SMRY -k "$keyword" "$producer_base" "$consumer_base" 2e-2 2e-2
done

# Field and group aggregates must match the parent even though P2 lies outside
# the sector: its rates are injected into the summary aggregation from the
# parent run. WOPR:P2 is reported straight from the parent series.
#
# The cumulatives are the sharper check of the two. They are built by
# accumulation rather than assignment, so a cumulative that was seeded instead
# of accumulated would show up here as double counting even when the rates
# still agree.
for keyword in FOPR FOPT GOPR:G1 GOPT:G1 WOPR:P2 WOPT:P1 WOPT:P2; do
    "$compare_bin" -d -t SMRY -k "$keyword" "$producer_base" "$consumer_base" 2e-2 2e-2
done
