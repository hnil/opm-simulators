#!/usr/bin/env bash
set -euo pipefail

# The --fluxnum route: build a sector boundary file from a parent run's restart
# and summary alone, with the region given as a FLUXNUM array rather than by
# running the parent with DUMPFLUX.
#
# The first check is the one that matters. The parent deck's own FLUXNUM selects
# a single cell, and the supplied grdecl selects the same one, so the file built
# from the restart must match the file the parent wrote for itself. Anything the
# route gets wrong about geometry, transmissibility, exterior state or PVT
# regions shows up as a mismatch.

make_flux_bin="$1"
compare_bin="$2"
shape_bin="$3"
parent_deck="$4"
one_sector="$5"
two_sectors="$6"

parent_base="$(basename "$parent_deck" .DATA)"

if [[ ! -f "${parent_base}.FLUX" ]]; then
    echo "expected ${parent_base}.FLUX from the parent DUMPFLUX run" >&2
    exit 1
fi

# make_flux resolves the case from the deck, which the driver leaves behind in
# the source directory rather than here.
cp "$parent_deck" "${parent_base}.DATA"

# 1. Same region, built without DUMPFLUX, must equal the parent's own file.
"$make_flux_bin" --parent="$parent_base" \
                 --fluxnum="$one_sector" \
                 --regions=1 \
                 --output=FROM_FLUXNUM.FLUX > fluxnum_one.log 2>&1

"$compare_bin" "${parent_base}.FLUX" FROM_FLUXNUM.FLUX

# The route must carry the exterior PVT region and the external region sums;
# without the latter a RESV well in the sector gets a different target than it
# had in the full model.
for array in FLXPVTN FLXRCON; do
    if ! grep -q "$array" <(strings FROM_FLUXNUM.FLUX); then
        echo "check-make-flux-fluxnum-smoke: $array missing from the generated file" >&2
        exit 1
    fi
done

# 2. Several regions in one run write one file each, using the numbered names
#    the consumer's input selection already understands.
"$make_flux_bin" --parent="$parent_base" \
                 --fluxnum="$two_sectors" \
                 --regions=1,2 \
                 --output=MULTI.FLUX > fluxnum_multi.log 2>&1

for suffix in 0001 0002; do
    if [[ ! -f "MULTI.FLUX${suffix}" ]]; then
        echo "check-make-flux-fluxnum-smoke: expected MULTI.FLUX${suffix}" >&2
        exit 1
    fi
done

# The parent grid is 10x1x10. Region 2 is the single cell at i 7, k 7, which is
# the shape the check tool can express; it confirms the header box, mode and
# face count of a file the route built on its own.
"$shape_bin" MULTI.FLUX0002 10 1 10 7 1 7 pressure

# Region 1 spans four cells, so check it structurally instead: it must be a
# readable FLUX file carrying the same boundary arrays.
for array in FLUXHEAD LOC2GLOB FACECELL FACEDIR FLUXTRAN FLXPRES; do
    if ! grep -q "$array" <(strings MULTI.FLUX0001); then
        echo "check-make-flux-fluxnum-smoke: $array missing from MULTI.FLUX0001" >&2
        exit 1
    fi
done

# A plain --output name must survive when only one region is asked for.
if [[ -f "FROM_FLUXNUM.FLUX0001" ]]; then
    echo "check-make-flux-fluxnum-smoke: a single region should keep the plain name" >&2
    exit 1
fi

# 3. Guard rails. Each of these must be refused rather than silently doing
#    something surprising.
expect_failure() {
    local description="$1"
    shift

    if "$@" > guard.log 2>&1; then
        echo "check-make-flux-fluxnum-smoke: expected failure for $description" >&2
        exit 1
    fi
}

expect_failure "--fluxnum together with --mapping-inline" \
    "$make_flux_bin" --parent="$parent_base" --fluxnum="$one_sector" \
                     --regions=1 --mapping-inline="global 10" --output=GUARD.FLUX

expect_failure "--fluxnum without --regions" \
    "$make_flux_bin" --parent="$parent_base" --fluxnum="$one_sector" --output=GUARD.FLUX

expect_failure "flux mode on the fluxnum route" \
    "$make_flux_bin" --parent="$parent_base" --fluxnum="$one_sector" \
                     --regions=1 --mode=flux --output=GUARD.FLUX

expect_failure "a region the array does not hold" \
    "$make_flux_bin" --parent="$parent_base" --fluxnum="$one_sector" \
                     --regions=7 --output=GUARD.FLUX

expect_failure "--regions without --fluxnum" \
    "$make_flux_bin" --parent="$parent_base" --mapping-inline="global 10" \
                     --regions=1 --output=GUARD.FLUX
