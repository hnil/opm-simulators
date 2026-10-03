#!/bin/bash

# A deck with FLUXNUM but without DUMPFLUX is an ordinary run. The well in it
# straddles two FLUXNUM regions, which a DUMPFLUX run would reject, and no
# FLUX file may be written.

set -euo pipefail

base="FLUX_FLUXNUM_ONLY_SMOKE"

if ! grep -q "End of simulation" "${base}.PRT"; then
    echo "The run with FLUXNUM but no DUMPFLUX did not complete" >&2
    exit 1
fi

if compgen -G "${base}.FLUX*" > /dev/null; then
    echo "A FLUX file was written although the deck has no DUMPFLUX" >&2
    ls "${base}".FLUX* >&2
    exit 1
fi
