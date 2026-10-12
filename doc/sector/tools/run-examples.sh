#!/bin/bash
# Make the runs of the worked examples on the Drogon model.
#
# Usage: run-examples.sh <drogon0 directory> <scratch directory>
#
# The drogon0 directory is the one holding DUMP-SECTOR.DATA,
# DUMP-SECTOR-R2.DATA and sector/job/SECTOR-UF.DATA. Nothing is written
# there: the decks are copied to the scratch directory, and the include
# folders linked. FLOW and MAKE_FLUX may be set to the executables to use,
# FLOW_ARGS to extra arguments such as --threads-per-process.
set -euo pipefail

src="$(cd "$1" && pwd)"
run="$2"
flow="${FLOW:-flow}"
make_flux="${MAKE_FLUX:-make_flux}"
read -r -a flow_args <<< "${FLOW_ARGS:-}"
box='box 10 37 1 41 64 32'

mkdir -p "${run}/sector/job"
run="$(cd "${run}" && pwd)"
ln -sfn "${src}/include" "${run}/include"
ln -sfn "${src}/sector/job/include" "${run}/sector/job/include"
cp "${src}"/DUMP-SECTOR.DATA "${src}"/DUMP-SECTOR-R2.DATA "${src}"/USE-SECTOR.DATA "${run}/"
cp "${src}"/sector/job/SECTOR-UF.DATA "${run}/sector/job/"

cd "${run}"

# Flux-mode variants of example 1.
sed "s/^  'PRESSURE' \//  'FLUX' \//" DUMP-SECTOR.DATA > DUMP-SECTOR-FLUX.DATA
sed -e "s/^  'PRESSURE' \//  'FLUX' \//" -e "s/^  'DUMP-SECTOR' \//  'DUMP-SECTOR-FLUX' \//" \
    USE-SECTOR.DATA > USE-SECTOR-FLUX.DATA

# The cut-out sector on the other FLUX files of examples 2 and 3.
for v in R2:TEST-R2 DF:DUMP-SECTOR UDQ:TEST-UDQ; do
    sed "s#^  '../../TEST' /#  '../../${v#*:}' /#" sector/job/SECTOR-UF.DATA \
        > "sector/job/SECTOR-UF-${v%%:*}.DATA"
done

# Example 3: two UDQs ahead of the prediction schedule, reported in SUMMARY.
cp DUMP-SECTOR-R2.DATA DUMP-SECTOR-UDQ.DATA
python3 - <<'PY'
import re

udq = """-- Oil rate of the producers outside and inside the sector
UDQ
  DEFINE FUOUTOPR  WOPR 'A1' + WOPR 'A2' /
  UNITS  FUOUTOPR  SM3/DAY /
  DEFINE FUSECOPR  WOPR 'A3' + WOPR 'A4' /
  UNITS  FUSECOPR  SM3/DAY /
/

"""

def add_udqs(deck):
    s = open(deck).read()
    m = re.search(r"INCLUDE\s*\n\s*'include/drogon_pred\.sch'", s)
    s = s[:m.start()] + udq + s[m.start():]
    m = re.search(r"^SUMMARY\s*\n(?:.*\n)*?ALL\s*\n", s, re.M)
    s = s[:m.end()] + "FUOUTOPR\nFUSECOPR\n" + s[m.end():]
    open(deck, "w").write(s)

add_udqs("DUMP-SECTOR-UDQ.DATA")
add_udqs("sector/job/SECTOR-UF-UDQ.DATA")
PY

run_flow() {
    "${flow}" "$1.DATA" --enable-esmry=true "${flow_args[@]}" > "$1.log" 2>&1
}

for deck in DUMP-SECTOR DUMP-SECTOR-FLUX DUMP-SECTOR-R2 DUMP-SECTOR-UDQ; do
    run_flow "${deck}"
done

"${make_flux}" --parent=DUMP-SECTOR --mapping-inline="${box}" --output=TEST.FLUX > make_flux.log 2>&1
"${make_flux}" --parent=DUMP-SECTOR-R2 --mapping-inline="${box}" --output=TEST-R2.FLUX > make_flux_r2.log 2>&1
"${make_flux}" --parent=DUMP-SECTOR-UDQ --mapping-inline="${box}" --output=TEST-UDQ.FLUX > make_flux_udq.log 2>&1

for deck in USE-SECTOR USE-SECTOR-FLUX; do
    run_flow "${deck}"
done

cd sector/job
for deck in SECTOR-UF SECTOR-UF-R2 SECTOR-UF-DF SECTOR-UF-UDQ; do
    run_flow "${deck}"
done

echo "All runs done in ${run}"
