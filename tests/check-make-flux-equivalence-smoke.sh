#!/usr/bin/env bash
set -euo pipefail

make_flux_bin="$1"
compare_bin="$2"
parent_deck="$3"
parent_root="$4"
mode="$5"
mapping_text="$6"
mapping_text="${mapping_text#\'}"
mapping_text="${mapping_text%\'}"

expected_flux_file=""
compare_args=()
no_summary=0

for extra in "${@:7}"; do
    case "$extra" in
        --expected=*)
            expected_flux_file="${extra#--expected=}"
            ;;
        --compare-arg=*)
            compare_args+=("${extra#--compare-arg=}")
            ;;
        --no-summary)
            no_summary=1
            ;;
        *)
            if [[ -z "$expected_flux_file" ]]; then
                expected_flux_file="$extra"
            else
                echo "Unexpected argument '$extra'" >&2
                exit 1
            fi
            ;;
    esac
done

if [[ -n "$expected_flux_file" && ! -f "$expected_flux_file" ]]; then
    expected_flux_file="$PWD/$expected_flux_file"
fi

if [[ -z "$expected_flux_file" ]]; then
    if [[ -f "${parent_root}.FLUX" ]]; then
        expected_flux_file="$PWD/${parent_root}.FLUX"
    else
        shopt -s nullglob
        candidates=("$PWD/${parent_root}.FLUX"[0-9][0-9][0-9][0-9])
        shopt -u nullglob
        if (( ${#candidates[@]} == 0 )); then
            echo "Expected ${parent_root}.FLUX or ${parent_root}.FLUXdddd in test result directory" >&2
            exit 1
        fi
        expected_flux_file="${candidates[0]}"
    fi
fi

if [[ ! -f "$expected_flux_file" ]]; then
    echo "Expected FLUX reference file '$expected_flux_file'" >&2
    exit 1
fi

cp "$parent_deck" "${parent_root}.DATA"

mapping_args=()
if [[ "$mapping_text" == inline:* ]]; then
    mapping_inline="${mapping_text#inline:}"
    mapping_args+=("--mapping-inline=$mapping_inline")
else
    printf '%s\n' "$mapping_text" > "${parent_root}.MAP"
    mapping_args+=("--mapping=$PWD/${parent_root}.MAP")
fi

summary_args=()
if (( no_summary == 1 )); then
    summary_args+=("--no-summary")
fi

"$make_flux_bin" \
    --parent="$PWD/${parent_root}" \
    "${mapping_args[@]}" \
    "${summary_args[@]}" \
    --output="$PWD/${parent_root}.MAKE.FLUX" \
    --mode="$mode" \
    > "make_flux_${mode}_equivalence.log" 2>&1

"$compare_bin" "$expected_flux_file" "${parent_root}.MAKE.FLUX" "${compare_args[@]}"
