#!/usr/bin/env python3
# Copyright 2026 Equinor ASA.
# This file is part of the Open Porous Media project (OPM), licensed under the
# GNU General Public License version 3 or later.
"""Write a two-component compositional deck from a black-oil deck.

The components are stock-tank oil (STO) and surface gas (GAS). flow_comp then
takes oil and gas properties from the deck's PVTO/PVDG tables (black-oil PVT
mode of GenericOilGasWaterFluidSystem). The EOS keywords written here are
placeholders: they are read, but not used for the phase behaviour.

Usage: bo_to_comp_deck.py IN.DATA OUT_BO.DATA OUT_COMP.DATA

OUT_BO is the black-oil reference: IN without DRSDT, since an equilibrium
flash always lets free gas re-dissolve.

Handles decks of the SPE1 kind: one equilibration region, RSVD, oil and gas
surface densities in DENSITY, FIELD or METRIC units.
"""
import re
import sys

MW_OIL = 200.0  # arbitrary; the results do not depend on it


def strip_comment(line):
    return line.split('--', 1)[0]


def keyword_records(lines, name):
    """Data lines of the first occurrence of keyword `name`, until the slash
    that closes it (a lone '/' or the last record of a single-record keyword)."""
    out, inside = [], False
    for line in lines:
        s = strip_comment(line).strip()
        if not inside:
            if s == name:
                inside = True
            continue
        if not s:
            continue
        if re.match(r'^[A-Z][A-Z0-9_]{0,7}$', s):
            break
        out.append(s)
    return out


def numbers(records):
    vals = []
    for r in records:
        for tok in r.replace('/', ' / ').split():
            if tok == '/':
                vals.append('/')
            else:
                vals.append(float(tok))
    return vals


def gas_injectors(lines):
    recs = keyword_records(lines, 'WCONINJE')
    return [r.split()[0].strip("'") for r in recs
            if r != '/' and len(r.split()) > 1 and r.split()[1].strip("'") == 'GAS']


def main(src, out_bo, out_comp):
    lines = open(src).read().splitlines()
    field = any(strip_comment(l).strip() == 'FIELD' for l in lines)

    rho_o, _rho_w, rho_g = [v for v in numbers(keyword_records(lines, 'DENSITY')) if v != '/'][:3]
    # Gas molar mass from the ideal gas at standard conditions.
    if field:
        vm_std = 379.48  # scf/lb-mol at 60 F, 14.696 psia
        rs_scale = 1000.0  # Mscf/stb -> scf/stb
        oil_vol = 5.6146  # ft3/stb
    else:
        vm_std = 23.645  # sm3/kmol at 15.56 C, 1.01325 bar
        rs_scale = 1.0
        oil_vol = 1.0
    mw_gas = rho_g * vm_std

    rsvd = [v for v in numbers(keyword_records(lines, 'RSVD')) if v != '/']
    depth_rs = list(zip(rsvd[0::2], rsvd[1::2]))
    zmfvd = []
    for depth, rs in depth_rs:
        n_gas = rs * rs_scale * rho_g / mw_gas
        n_oil = oil_vol * rho_o / MW_OIL
        zg = n_gas / (n_gas + n_oil)
        zmfvd.append((depth, 1.0 - zg, zg))

    if field:
        stcond, rtemp = '60.0 14.696', '200.0'
        tcrit, pcrit = '1300.0 343.0', '250.0 667.0'
    else:
        stcond, rtemp = '15.56 1.01325', '93.0'
        tcrit, pcrit = '722.0 190.6', '17.2 46.0'

    # black-oil reference: drop DRSDT and its record
    bo, skip = [], False
    for line in lines:
        s = strip_comment(line).strip()
        if s == 'DRSDT':
            skip = True
            continue
        if skip:
            if s:
                skip = False
            continue
        bo.append(line)
    open(out_bo, 'w').write('\n'.join(bo) + '\n')

    comp, stream_done = [], False
    for line in bo:
        s = strip_comment(line).strip()
        if s == 'GRID':
            comp += ['COMPS', ' 2 /', '']
        if s == 'SOLUTION':
            comp += ['-- placeholder EOS: phase behaviour comes from PVTO/PVDG',
                     'EOS', ' PR /', 'CNAMES', " 'STO' 'GAS' /",
                     'MW', f' {MW_OIL} {mw_gas:.6f} /', 'TCRIT', f' {tcrit} /',
                     'PCRIT', f' {pcrit} /', 'ACF', ' 0.6 0.011 /',
                     'ZCRIT', ' 0.25 0.288 /', 'BIC', ' 0.0 /',
                     'STCOND', f' {stcond} /', 'ZMFVD']
            comp += [f' {d} {zo:.10f} {zg:.10f}' for d, zo, zg in zmfvd]
            comp += [' /', '']
        if s == 'SUMMARY':
            comp += ['RTEMP', f' {rtemp} /', '']
        if s == 'TSTEP' and not stream_done:
            comp += ['WELLSTRE', " 'SEPGAS' 0.0 1.0 /", '/', 'WINJGAS']
            comp += [f" '{w}' STREAM SEPGAS /" for w in gas_injectors(bo)]
            comp += ['/', '']
            stream_done = True
        comp.append(line)
    open(out_comp, 'w').write('\n'.join(comp) + '\n')
    print(f'MW gas {mw_gas:.4f}, ZMFVD {zmfvd}')


if __name__ == '__main__':
    main(*sys.argv[1:4])
