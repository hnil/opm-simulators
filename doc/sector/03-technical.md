# Technical documentation {#technical-documentation}

This chapter describes how sector modelling is implemented, for developers
and for users who need to know exactly what a reduced run does. File and
class names refer to opm-common and opm-simulators. In the table below,
`flux/` stands for `opm/simulators/flow/flux/` in opm-simulators.

## Overview of the implementation

| Component | Where | Role |
|:--|:--|:--|
| `FluxFile` | opm-common, `opm/io/eclipse` | Reads and writes the FLUX file format. |
| `IOConfig` | opm-common | Holds DUMPFLUX, USEFLUX, FLUXREG and FLUXTYPE. |
| `SummaryConfig` | opm-common | Expands FLUXALL; adds FPRFLUX to reduced runs. |
| `Schedule::insertKeywords` | opm-common | Adds the parent's wells to a reduced run's schedule. |
| `FluxRegions` | `flux/` | Cells and boundary faces of each FLUXNUM region. |
| `FluxDumper` | `flux/` | Collects and appends one region's boundary data. |
| `FluxBoundary` | `flux/` | Maps a FLUX file onto the reduced run's grid. |
| `FluxActivation` | `flux/` | Narrows a USEFLUX deck's ACTNUM to the selected region. |
| `ParentSummary` | `flux/` | Time series of the parent's summary vectors. |
| `FluxSummaryKeys` | `flux/` | Which summary vectors go into a FLUX file. |
| `FluxParentWells` | `flux/` | Describes the parent's wells, and adds missing ones to a schedule. |
| `EclWriter` | opm-simulators | The DUMPFLUX producer, and parent values in the reduced run's summary. |
| `FlowProblem`, `FlowProblemBlackoil` | opm-simulators | The USEFLUX boundary conditions and parent coupling. |
| `make_flux`, `inspect_flux` | opm-simulators, `examples` | Offline producer and inspection tool. |

## Sectors and boundary faces

A sector is the set of active cells with a given positive FLUXNUM value.
`FluxRegions::extract` finds, for every region, its cells and the faces
between a cell in the region and an active cell outside it. Three kinds of
face are found:

* Faces between Cartesian neighbours, identified by the interior cell and
  the direction of the face, I−, I+, J−, J+, K− or K+.
* Non-neighbour connections (NNCs) from the deck, from pinch-outs and from
  the grid processing, such as faults with throw, where one cell is inside the
  region and the other outside.
* Faces towards inactive cells are not boundary faces, since nothing flows
  through them in the parent either.

Each face records the global (Cartesian) index of the interior and the
exterior cell. A DUMPFLUX run requires every well to be completed within a
single region, and stops with a list of the offending wells otherwise, since
a reduced run would only see the part of a straddling well that lies inside
its sector.

## The producer: DUMPFLUX

`EclWriter` sets up one `FluxDumper` per region when the run has DUMPFLUX
(`IOConfig::getDumpFlux()`) and a FLUXNUM array, and from then on collects
the boundary data on every rank for the cells that rank owns. The data is
reduced onto the I/O rank, so a FLUX file written in parallel is the same as
one written serially.

### Rates

At the end of every converged time step, `sampleFluxDumperRates()` reads the
phase volumetric fluxes (FLORES) of every boundary face from the linearizer,
which holds the converged fluxes of the step. They are turned into component
masses with the state of the upstream cell: the oil component is
$(q_o b_o + R_v q_g b_g)\,\rho_o^{ref}$, the gas component
$(q_g b_g + R_s q_o b_o)\,\rho_g^{ref}$ and water $q_w b_w \rho_w^{ref}$,
where $b$ is the inverse formation volume factor of the upstream cell. The
sign convention is positive into the sector. The same masses are also read
directly from the residual's component fluxes (FLOWS), and the run warns if
the two disagree, since the file would then not describe what was
transported.

The masses are accumulated over the steps between two records and written as
time averages, so that a reduced run imposing them over the record's
interval moves exactly the mass the parent moved.

DUMPFLUX makes the run compute FLORES and FLOWS internally whether or not
RPTRST asks for them, including for NNCs, without writing them to the restart
file.

### Exterior state

For a pressure boundary, `updateFluxDumpers_()` records for every face the
state of the exterior cell at the end of the record's interval:

* pressure of the oil phase, or of the first active phase without oil,
* water and gas saturation, Rs, Rv and temperature,
* relative permeability and capillary pressure of every active phase, read
  from the converged state, so that the cell's saturation region, scaled end
  points and hysteresis state are all included.

The static part of the file records, per face, the transmissibility the
parent had across the face, the PVT region, the equilibration region and the
depth of the exterior cell.

### Other data

* **Pore-volume weighted sums outside the sector.** For every record, sums of
  pressure, temperature, Rs, Rv, Rsw, Rvw, pore volume and salt
  concentration over all cells outside the region, weighted by hydrocarbon
  pore volume and by total pore volume. They are stored as sums rather than
  averages, so that the reduced run can add its own cells to them.
* **Threshold pressures.** The full table of threshold pressures between
  equilibration regions, as the parent arrived at it after equilibration.
* **Parent wells.** For every report step: well names, the group of each
  well, its type, status and efficiency factor, the group tree and group
  efficiency factors.
* **Parent summary.** The summary vectors chosen by `fluxSummaryKeys()`,
  sampled from the run's summary state after every time step. Rate vectors
  are averaged over the interval since the previous sample.

### Cadence

Boundary records and summary samples are written at every time step, but no
closer together than one day by default, and always at the end of a report
step. Records are appended to the file as they are formed, so a long run
does not hold its boundary data in memory and a FLUX file can be read while
the parent is still running. The cadence is controlled by the [run-time
parameters](#run-time-parameters).

## The FLUX file format

A FLUX file is an ECL-format keyword file, binary or formatted (`.FFLUX`),
which can be read with any tool that reads restart files, for instance with
`opm.io.ecl.EclFile` in Python, and converted between the two forms with
`convertECL`. It has a static section, written once, followed by any number
of record blocks and summary blocks in the order they were formed. All values
are in SI units; times are in seconds from the start of the run.

### Static section

| Array | Type | Size | Contents |
|:--|:--|:--|:--|
| `FLUXHEAD` | INTE | 22 | Header, see below. |
| `FLUXNAMS` | C0nn | names + keys | Parent case name and region number, followed by the summary keys. |
| `FLUXNCNT` | INTE | 2 | Number of names, and number of summary keys, in FLUXNAMS. |
| `LOC2GLOB` | INTE | cells | Global (Cartesian, zero-based) index of every cell of the region. |
| `FACECELL` | INTE | faces | Interior cell of every face, as an index into LOC2GLOB. |
| `FACEDIR` | INTE | faces | Face direction: 1 I+, 2 I−, 4 J+, 8 J−, 16 K+, 32 K−, 0 for an NNC. |
| `FACEGLNB` | INTE | faces | Global index of the exterior cell. |
| `FLUXTRAN` | DOUB | faces | Transmissibility across the face, in the parent. |
| `FLXPVTN` | INTE | faces | PVT region of the exterior cell. |
| `FLXEXDP` | DOUB | faces | Depth of the centre of the exterior cell. |
| `FLXEQLN` | INTE | faces | Equilibration region of the exterior cell, zero-based. |
| `FLXTHPR` | DOUB | n² | Threshold pressures; entry `r1*n + r2` is for flow from region r1 to r2. Optional. |
| `FLXMINT` | DOUB | 1 | Minimum interval between boundary records the producer enforced. |
| `FLXPSTRT` … `FLXGEFAC` | | | Parent wells, see below. Optional. |
| `SMRYMINT` | DOUB | 1 | Minimum interval between summary samples. Present when summary keys are. |

The header `FLUXHEAD` holds, in order: format version (currently 7), parent
grid dimensions NX, NY, NZ, the one-based corner I1, J1, K1 and the
dimensions of the region's bounding box, number of cells, number of boundary
faces, number of records (zero, derived on read), number of phases, a
temperature flag, the mode (1 flux, 2 pressure, 3 both), the sampling (0
instant, 1 averaged), the phase mask (1 oil, 2 water, 4 gas), the number of
summary keys, and three fields derived on read.

### Parent wells

| Array | Size | Contents |
|:--|:--|:--|
| `FLXPSTRT` | steps | Start time of every report step. |
| `FLXWELLS` | wells | Well names. |
| `FLXGRUPS` | groups | Group names. |
| `FLXWGRP` | steps × wells | Group of each well, one-based into FLXGRUPS, 0 if not defined yet. |
| `FLXWTYPE` | steps × wells | 1 producer, 2 water, 3 gas, 4 oil injector. |
| `FLXWSTAT` | steps × wells | 1 open, 2 stop, 3 shut, 4 auto. |
| `FLXWEFAC` | steps × wells | Well efficiency factor. |
| `FLXWDREF` | wells | Reference depth of the bottom-hole pressure. |
| `FLXGPAR` | steps × groups | Parent group, one-based into FLXGRUPS, 0 for FIELD. |
| `FLXGEFAC` | steps × groups | Group efficiency factor. |

### Record blocks

Each block holds one or more boundary records and starts with `FLXSEQ`,
the block number and the number of records in it. Per-face arrays are
face-major: all phases, or all values, of the first face come first.

| Array | Size per record | Contents |
|:--|:--|:--|
| `FLXSEQ` | (2 per block) | Block number and number of records. |
| `FLXSTEP`, `FLXSIM` | 1 | Report step and simulator step. |
| `FLXTIME`, `FLXDT` | 1 | Start time and length of the interval the record covers. |
| `FLXMASS` | faces × phases | Component masses per second, positive into the sector. Flux mode. |
| `FLXPRES` | faces | Exterior pressure. Pressure mode. |
| `FLXSATW`, `FLXSATG` | faces | Exterior water and gas saturation. Pressure mode. |
| `FLXRS`, `FLXRV`, `FLXTEMP` | faces | Exterior Rs, Rv and temperature. Pressure mode, optional. |
| `FLXKR`, `FLXPC` | faces × phases | Exterior relative permeability and capillary pressure, phases in the order oil, water, gas. Pressure mode, optional. |
| `FLXRCON` | 16 | Pore-volume weighted sums outside the region. |

A rate stored in a record is the average over the record's interval, and a
state is the state at its end.

### Summary blocks

| Array | Size | Contents |
|:--|:--|:--|
| `SMRYSEQ` | 2 | Block number and number of samples. |
| `SMRYTIME` | samples | End time of the interval each sample covers. |
| `SMRYVALS` | samples × keys | Values in the parent's output units, sample-major. Rates are interval averages. |

`inspect_flux` prints the contents of a file. For one cell it gives the
boundary data of every face of the cell as a time series, with the face
geometry on standard error:

```
$ inspect_flux DUMP-SECTOR.FLUX 10 50 5 --start-date=2018-01-01
inspect_flux: DUMP-SECTOR.FLUX
  parent case      DUMP-SECTOR
  region           1
  parent grid      46x73x32
  records          157
  faces on cell    1 of 2427 in the region
  I-     exterior cell     15694       (9,50,5)  transmissibility 4.664e-13 m3.Pa.s  pvt region 1
#       DATE      I-:PRES      I-:SWAT      I-:SGAS        I-:RS        I-:RV  ...
  2018-01-01    3.015e+07       0.6966            0        140.8    0.0002351  ...
```

and with `--summary` it lists embedded summary vectors.

## The consumer: USEFLUX

### Setting up the reduced grid

USEFLUX is processed before the schedule is built. `applyUseFluxActnum()`
sets ACTNUM to zero for every cell outside the regions FLUXREG selects
(intersected with the deck's own ACTNUM), so that the grid, the
transmissibilities and the well connections are all formed for the sector
alone. A cut-out sector deck covers just the sector, and its FLUXNUM is one
everywhere.

`FlowProblem::readFluxBoundary_()` then reads the FLUX file and checks it.
The file is `<name>.FLUX` in the deck's directory, where `<name>` is the
USEFLUX item or the deck's own base name; if no such file exists, the
lowest-numbered `<name>.FLUXnnnn` is used. The parent's cell indices are
mapped onto the reduced grid, which must have the dimensions of either the
parent grid or the region's bounding box recorded in the header. Every
boundary face is matched to a face of the reduced grid; faces that cannot be
matched are reported, since whatever crossed them would be lost.

### Choosing the record

Before every time step, `refreshFluxBoundaryRecord_()` selects the record
covering the start of the step. In pressure mode, the boundary state is
interpolated linearly between that record and the previous one, evaluated at
the end of the step, since the step is implicit. Where a step ends on a
record's time, the recorded state is used as it is. Rates are averages over
their interval and are not interpolated.

### Flux boundary

In flux mode, the boundary condition of a face is the recorded mass rate of
each component, divided by the face area, with the sign changed to the
outward normal that the discretisation uses. Since the rates are component
masses, nothing about the fluid outside the sector is needed to impose them.

### Pressure boundary

In pressure mode, every face the file describes is a free boundary, and the
flux across it is the two-point flux the parent would have computed:

$$
q_\alpha = T\,\lambda_\alpha^{up}\,\big(p_\alpha^{ext} - p_\alpha^{int}
           + \bar\rho_\alpha g\,(z^{int} - z^{face})\big)
$$

with these ingredients:

* **Transmissibility.** The parent's transmissibility across the face,
  `FLUXTRAN`, replaces the default transmissibility of the outer boundary.
* **Exterior state.** `fluxExteriorFluidState_()` builds the fluid state on
  the far side from the record: pressure, saturations, Rs, Rv and
  temperature, the phase pressures from the recorded capillary pressures, and
  densities, formation volume factors and viscosities evaluated with the
  exterior cell's PVT region.
* **Depth.** The recorded pressure is that of the exterior cell at its
  centre. It is carried hydrostatically, with the exterior fluid's densities,
  to the face where the boundary condition is imposed. Without this, a
  sector in a dipping reservoir sees an error of the hydrostatic head between
  the two depths.
* **Upstream mobility.** When fluid flows into the sector, its relative
  permeability is the exterior cell's recorded value,
  `boundaryRelativePermeabilities()`, not this cell's saturation functions
  evaluated at the exterior saturations.
* **Upstream direction.** The upstream side is decided from the
  gravity-corrected potential difference, consistently with the flux itself.
* **Threshold pressures.** If the deck has threshold pressures, the
  threshold between the interior cell's and the exterior cell's
  equilibration regions is applied across the face as across an interior
  face. The reduced run takes the parent's table of threshold pressures
  instead of computing its own, since a defaulted threshold pressure is the
  largest initial potential difference along a region boundary, of which the
  sector only holds part.

### NNC faces

A boundary face reached through a non-neighbour connection has no grid face
in the reduced run, since the exterior cell and the connection are gone. Its
flux is added as a source term in the interior cell instead: in flux mode the
recorded masses, in pressure mode the two-point flux above, computed in
`addFluxNncPressureSource_()` with the recorded transmissibility and the
cells' centre depths.

## Parent coupling in the reduced run

### Parent summary

The parent's summary vectors are read from the FLUX file into a
`ParentSummary` series. If the file holds none, the parent's own summary
file next to it (`.ESMRY`, `.SMSPEC` or `.FSMSPEC`) is used, at the
resolution it was written with. Values are looked up by time: rates are held
constant over the parent's sample intervals, where they are interval
averages, and all other vectors are interpolated linearly.

### Wells outside the sector

A well is outside the sector if none of its connections are in the reduced
grid. Such a well is never solved for. Instead:

* **Summary output.** Its rates over the step, and its bottom-hole and
  tubing-head pressure at the end of the step, are taken from the parent
  summary and given to the summary evaluation as if the well had been
  simulated (`injectParentOnlyWellData_()`). Field and group vectors are then
  formed over all wells in the usual way, and cumulatives accumulate from the
  rates.
* **Group control.** Its surface and reservoir rates are handed to the group
  state helper as fixed rates (`GroupStateHelper::setFixedWellRates()`). A
  well with fixed rates counts towards its groups' production and injection
  and comes off group targets before the rest is shared out, but is never
  given a share of a target itself. This is done at the start of every time
  step, with the parent's average rates over the step.

### Adding wells missing from a cut-out deck

A cut-out sector deck usually keeps only its own wells. Before the schedule
is completed, `addFluxParentWells()` compares the deck's wells with the
parent wells in the FLUX file and, for every well the deck does not define,
generates at each report step where something changes:

* `GRUPTREE` and `GEFAC` for groups the deck lacks,
* `WELSPECS` placing the well in its parent group, with its head at (1, 1),
* `WCONPROD` or `WCONINJE` with the parent's status and a zero target,
* `WEFAC` with the parent's efficiency factor.

The keywords are inserted with `Schedule::insertKeywords()`, which rebuilds
the schedule once from the first affected report step, so that the deck's
own keywords that refer to an added well are handled with the well in place.
The controls are placeholders: the added wells have no connections, and
their rates come from the parent summary as for any well outside the sector.
While the schedule is being read, the checks that would reject a UDQ over a
well that does not exist yet, or more wells and groups than WELLDIMS allows,
are turned into warnings.

### Summary vectors and UDQs

At the start of every time step, the parent's values of all embedded
non-cumulative vectors are written into the reduced run's summary state
(`seedParentSummaryState_()`), so that UDAs, group controls and ACTIONX
conditions see them. Cumulatives are not seeded: they are built up from the
rates by the normal accumulation, which would otherwise count the parent's
production twice.

The summary evaluation then overwrites every vector the reduced run reports
with what it computes itself. For quantities that lie outside the sector,
that is wrong: a well with no connections, a region the sector holds only
part of, a block that is inactive. `applyParentOnlySummaryValues_()` runs
between the summary evaluation and the UDQ evaluation, and puts the parent's
values back for

* well, connection and segment vectors of wells outside the sector, except
  rates and pressures already supplied as described above,
* field and group vectors that are not sums of well rates, such as pressures,
* region vectors,
* block vectors of cells that are not active in the reduced run.

The UDQs are then evaluated from a summary state that holds the parent's
values for everything the sector cannot compute, and the reduced run's own
values for everything it can. UDQs that refer to regions the sector does not
hold are accepted, the region sets being extended with the regions the
parent reported (`FIPRegionStatistics::raiseMaximumRegionID()`). A UDQ used
as a group target keeps its start-of-step value, as in the parent, and is
the same on every rank.

### Field pressure and rate conversion

The external pore-volume weighted sums of the record covering the step are
added to the reduced run's own sums in two places: in the rate converter that
turns reservoir voidage rates into surface rates for wells on RESV control,
and in the FPRFLUX summary vector. Both therefore average over the whole
parent model, while still reflecting changes inside the sector.

## FLUXALL and summary embedding

The summary vectors embedded in a FLUX file are chosen by
`fluxSummaryKeys()`:

* well rates and cumulatives at surface and reservoir conditions, for
  production and injection of every active phase, expanded over all wells,
* every keyword referred to by a UDQ or an ACTIONX condition, expanded over
  the wells and groups the schedule defines,
* region, block, connection, segment, aquifer and network node vectors that
  the UDQs or ACTIONX conditions refer to and the parent reports,
* the UDQs that other UDQs or ACTIONX conditions refer to.

A vector can only be embedded if the parent reports it. FLUXALL in the
parent's SUMMARY section makes it report exactly this list: UDQ and ACTIONX
expressions now report the vectors they name in full, such as `BPR 10 10 3`
or `SOFR 'P1' 12`, and FLUXALL expands the rest over the objects they can
apply to. Whatever FLUXALL cannot expand is listed in the print file.

## make_flux

make_flux builds the same file from the parent's output:

* **Region.** From a FLUXNUM array in a grdecl file (`--fluxnum` and
  `--regions`), or from a mapping (`--mapping` or `--mapping-inline`) of
  boxes, cells or index ranges.
* **Records.** One per restart step of the parent. Pressure, saturations,
  Rs and Rv are read from the restart. The exterior relative permeabilities
  and capillary pressures are evaluated with the simulator's material law
  from the parent's saturation functions and scaled end points; the
  hysteresis state is not in the restart, and the drainage curves are used.
* **Flux mode.** Reads the FLOWS arrays of the restart (component surface
  volumes per face, including NNCs) and multiplies them by the reference
  density. It is not available with `--fluxnum`.
* **Threshold pressures** from the THRESHPR restart array. If the parent has
  THPRES but wrote no THRESHPR, which is the case for restarts written with
  `NORST=1`, make_flux warns, and the reduced run computes its own.
* **Parent wells and summary** as in DUMPFLUX, from the deck and the summary
  file. `--summary-min-interval` thins the summary, averaging rates over
  what each sample covers, and `--smry-define` supplies a vector the parent
  did not write.

## Run-time parameters {#run-time-parameters}

`--flux-boundary-min-interval-between-samples` (default 1 day)
:   Minimum time between boundary records written by DUMPFLUX. 0 writes
    every time step.

`--flux-boundary-report-steps-only` (default false)
:   Write boundary records at report steps only.

`--flux-summary-min-interval-between-samples` (default 1 day)
:   Minimum time between embedded summary samples. 0 samples every time
    step.

`--flux-boundary-interpolate` (default true)
:   In a reduced run, interpolate the boundary state in time between
    records.

Report step boundaries are always written, whatever the interval.

## Tests

The tests are in `fluxTests.cmake`, with the decks in `tests/flux` and
`tests/data/sector`. Unit tests cover the file format, region extraction,
the dumper, the boundary reader and the parent summary. The integration
tests run a parent and a reduced run and compare them, for both boundary
types, in serial and in parallel, with NNC and fault boundaries, threshold
pressures, interpolation in time, wells outside the sector under group
control, UDQs over regions outside the sector, and make_flux against
DUMPFLUX.
