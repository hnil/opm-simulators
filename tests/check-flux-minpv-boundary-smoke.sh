#!/usr/bin/env bash
set -euo pipefail

# A cell the simulator drops for having too small a pore volume is not a cell
# the sector can exchange anything with, so the face towards it must not be in
# the FLUX file. The deck's own ACTNUM still calls that cell active, so
# anything deciding the boundary from the deck rather than from the grid the
# simulator built will write a face here.
#
# The sector is the single cell (1,1,2) of a 10x1x10 grid. Of its six faces:
#   I-  runs off the edge of the model
#   J-  and J+ likewise, the grid being one cell thick
#   I+  reaches (2,1,2), which MINPV removes
#   K-  reaches (1,1,1) and K+ reaches (1,1,3), both of which survive
# so exactly two faces should be described, and neither may name (2,1,2).

make_flux_bin="$1"
inspect_flux_bin="$2"
compare_bin="$3"
parent_deck="$4"
sector_grdecl="$5"

parent_base="$(basename "$parent_deck" .DATA)"

if [[ ! -f "${parent_base}.FLUX" ]]; then
    echo "expected ${parent_base}.FLUX from the parent DUMPFLUX run" >&2
    exit 1
fi

cp "$parent_deck" "${parent_base}.DATA"

# The cell MINPV removes, as inspect_flux reports exterior cells: global index
# 11, printed alongside its cartesian position.
removed_cell_ijk="(2,1,2)"

check_boundary() {
    local label="$1"
    local flux_file="$2"

    # inspect_flux puts the face geometry on stderr and the time series on
    # stdout, so the geometry is what gets examined here.
    "$inspect_flux_bin" "$flux_file" 1 1 2 2> "inspect_${label}.log" > "series_${label}.txt"

    local faces
    faces="$(grep -c "exterior cell" "inspect_${label}.log" || true)"

    if [[ "$faces" -ne 2 ]]; then
        echo "check-flux-minpv-boundary-smoke: $label describes $faces faces on the sector cell, expected 2" >&2
        cat "inspect_${label}.log" >&2
        exit 1
    fi

    if grep -q -- "$removed_cell_ijk" "inspect_${label}.log"; then
        echo "check-flux-minpv-boundary-smoke: $label has a face towards $removed_cell_ijk," \
             "which MINPV removed from the grid" >&2
        cat "inspect_${label}.log" >&2
        exit 1
    fi

    for direction in "K+" "K-"; do
        if ! grep -qE "^ +${direction//+/\\+} " "inspect_${label}.log"; then
            echo "check-flux-minpv-boundary-smoke: $label is missing the $direction face" >&2
            cat "inspect_${label}.log" >&2
            exit 1
        fi
    done

    # The time series must carry one row per record on top of the two headers.
    local rows
    rows="$(grep -vc '^#' "series_${label}.txt" || true)"
    if [[ "$rows" -lt 1 ]]; then
        echo "check-flux-minpv-boundary-smoke: $label produced no time series rows" >&2
        exit 1
    fi
}

# 1. What the simulator wrote for itself.
check_boundary dumpflux "${parent_base}.FLUX"

# 2. What the offline route builds from the same run's restart. It resolves the
#    active set from the EGRID, so it has to reach the same conclusion.
"$make_flux_bin" --parent="$parent_base" \
                 --fluxnum="$sector_grdecl" \
                 --regions=1 \
                 --output=MINPV_FROM_FLUXNUM.FLUX > make_flux_minpv.log 2>&1

check_boundary makeflux MINPV_FROM_FLUXNUM.FLUX

# 3. And the two must agree on everything, not merely on the face count.
"$compare_bin" "${parent_base}.FLUX" MINPV_FROM_FLUXNUM.FLUX
