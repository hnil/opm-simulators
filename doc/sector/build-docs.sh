#!/bin/bash
# Build the sector modelling guide as a single HTML page and a PDF, and the
# presentation as a PowerPoint file.
#
# Usage: build-docs.sh [output directory]
#
# PANDOC and PDF_ENGINE may be set to the executables to use, PYTHON to the
# interpreter for the presentation. The PDF is skipped, with a message, when
# no PDF engine is found, and the presentation when python-pptx is missing.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
out="${1:-${here}/build}"
pandoc="${PANDOC:-pandoc}"
pdf_engine="${PDF_ENGINE:-xelatex}"

mkdir -p "${out}"
cd "${here}"

chapters=(01-introduction.md 02-examples.md 03-technical.md 04-limitations.md)
common=(--from=markdown --metadata-file=metadata.yaml
        --lua-filter=tools/chapter-links.lua
        --resource-path="${here}"
        --toc --toc-depth=2 --number-sections)

# Pandoc 2.19 renamed --self-contained.
if "${pandoc}" --help | grep -q -- '--embed-resources'; then
    embed=(--embed-resources --standalone)
else
    embed=(--self-contained)
fi

"${pandoc}" "${common[@]}" --to=html5 "${embed[@]}" \
    --mathml --css=tools/style.css \
    --output="${out}/sector-modelling.html" "${chapters[@]}"
echo "Wrote ${out}/sector-modelling.html"

if command -v "${pdf_engine}" > /dev/null 2>&1; then
    "${pandoc}" "${common[@]}" --to=latex --pdf-engine="${pdf_engine}" \
        --variable=tables \
        --output="${out}/sector-modelling.pdf" "${chapters[@]}"
    echo "Wrote ${out}/sector-modelling.pdf"
else
    echo "No PDF engine '${pdf_engine}' found, PDF not built" >&2
fi

python="${PYTHON:-python3}"
if "${python}" -c 'import pptx, PIL' > /dev/null 2>&1; then
    "${python}" tools/make-slides.py "${out}/sector-modelling.pptx"
else
    echo "python-pptx or Pillow not found, presentation not built" >&2
fi
