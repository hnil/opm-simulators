#!/bin/bash

set -euo pipefail

flux_file="FLUX_DUMP_INVALID_TYPE_SMOKE.FLUX"
prt_file="FLUX_DUMP_INVALID_TYPE_SMOKE.PRT"

"$1" "${flux_file}" 10 1 10 10 1 5 FLUX
grep -q "Unrecognized FLUXTYPE BC_TYPE value 'BANANA', defaulting to FLUX output mode" "${prt_file}"
