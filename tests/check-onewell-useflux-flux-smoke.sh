#!/usr/bin/env bash
set -euo pipefail

flow_bin="$1"
compare_bin="$2"
producer_deck="$3"
consumer_deck="$4"

producer_base="$(basename "$producer_deck" .DATA)"
consumer_base="$(basename "$consumer_deck" .DATA)"

"$flow_bin" "$producer_deck" --output-dir="$(pwd)" > onewell_dump_flux.log 2>&1
cp "${producer_base}.FLUX" "${consumer_base}.FLUX"
cp "$consumer_deck" "${consumer_base}.DATA"
"$flow_bin" "${consumer_base}.DATA" --output-dir="$(pwd)" > onewell_use_flux.log 2>&1

for keyword in WOPR:P1 WBHP:P1 GOPR:G1 WWPR:P1 WGPR:P1; do
    "$compare_bin" -t SMRY -k "$keyword" "$producer_base" "$consumer_base" 1e-3 1e-3
done

grep -q "End of simulation" onewell_use_flux.log