# Introduction {#introduction}

Sector modelling lets a simulation run on part of a reservoir model, a
*sector*, while the rest of the model is represented by what crosses the
sector's boundary. The data describing that boundary is taken from a run of
the full model, called the *parent* run here, and stored in a FLUX file. A
run on the sector alone, the *reduced* run, reads the FLUX file and uses it as
boundary conditions. As long as the sector is left as it was, the reduced run
reproduces the parent's results for the sector to within a few tenths of a
percent, at a fraction of the cost. Changes made inside the sector, such as a
new well, a different completion or modified rock properties, can then be
studied without simulating the whole field again.

OPM Flow implements the keywords that Eclipse uses for this, DUMPFLUX,
USEFLUX, FLUXNUM, FLUXREG and FLUXTYPE, and extends them in a few ways that
are described below. This chapter gives an overview of the features. The
[worked examples](02-examples.md#worked-examples) show them on the Drogon
model, the [technical documentation](03-technical.md#technical-documentation)
describes how they work, and the [limitations](04-limitations.md#limitations)
section lists what is not supported.

## Workflow

A sector study has three parts.

1. **Define the sector.** The sector is one region of the FLUXNUM array,
   which is defined over the whole grid of the parent model, typically with
   EQUALS, a BOX or an INCLUDE file. FLUXNUM is an ordinary region array:
   defining it does not by itself make anything happen.

2. **Produce the boundary data.** There are two ways to obtain a FLUX file:

   * **DUMPFLUX.** The parent run is made with DUMPFLUX in its GRID section.
     While it runs it writes, for every FLUXNUM region, what crosses the
     region's boundary at every time step. This gives the most accurate
     boundary, but has to be decided before the parent run is made.
   * **make_flux.** A command-line tool that builds the FLUX file afterwards,
     from the output of an ordinary run: its deck, grid, INIT, restart and
     summary files. No rerun of the full model is needed, but the boundary
     is only described at the restart steps the parent wrote.

3. **Run the sector.** The reduced run has USEFLUX in its GRID section,
   naming the FLUX file, and FLUXREG selecting the region. It can either be
   the parent deck itself, where everything outside the selected region is
   simply made inactive, or a separate deck with a grid of its own that has
   been cut out of the parent model.

![Workflow of a sector study. A parent run with DUMPFLUX, or make_flux applied
to the output of an ordinary parent run, produces a FLUX file. A reduced run
with USEFLUX reads it, either on the parent deck or on a cut-out sector
deck.](images/workflow.png){width=95%}

## Keywords

| Keyword | Section | Purpose |
|:--|:--|:--|
| FLUXNUM | GRID | Region array defining the sectors, over the whole grid. Zero means no sector. |
| DUMPFLUX | GRID | Write a FLUX file for each FLUXNUM region during the run. |
| FLUXTYPE | GRID | `'FLUX'`, `'PRESSURE'` or `'BOTH'`: which kind of boundary data a DUMPFLUX run writes. Default `'FLUX'`. |
| USEFLUX | GRID | Run on one sector only, with boundary conditions from a FLUX file. The optional item names the file; by default it is the deck's own base name. |
| FLUXREG | GRID | Which FLUXNUM region the USEFLUX run covers. |
| FLUXALL | SUMMARY | In the parent run: report the summary vectors a reduced run needs. |
| FPRFLUX | SUMMARY | In a reduced run: field pressure averaged over the whole parent model. Added automatically to every USEFLUX run. |

A DUMPFLUX run writes `CASE.FLUX` when FLUXNUM has a single region, and
`CASE.FLUX0001`, `CASE.FLUX0002` and so on, one per region, when it has
several. A USEFLUX run reads one of those files. FLUXNUM alone, without
DUMPFLUX, is just a region array, as it is when MULTREGT uses it.

## Pressure and flux boundaries

The boundary of the sector can be driven in two ways. FLUXTYPE in the parent
deck decides which, and the reduced run follows whatever the file contains.

* **Pressure boundary** (`FLUXTYPE 'PRESSURE'`). The FLUX file records the
  state of the cell on the far side of every boundary face: its pressure,
  saturations, dissolved gas-oil ratio Rs, vaporised oil-gas ratio Rv,
  relative permeabilities and capillary pressures. The reduced run imposes
  that state on the far side of the face, and computes the flow across it
  from the transmissibility the parent had there. The flow therefore responds
  to what happens inside the sector: if a well is added near the boundary,
  more fluid is drawn in from outside. This is the mode to use for studying
  changes inside the sector.

* **Flux boundary** (`FLUXTYPE 'FLUX'`). The FLUX file records the mass of
  each component (oil, gas, water) that crossed every boundary face. The
  reduced run imposes those rates as they are. A flux boundary reproduces the
  parent very closely, since the same amount of fluid enters and leaves as in
  the parent, but it does not respond to changes inside the sector.

* **Both** (`FLUXTYPE 'BOTH'`). Both kinds of data are written to the file.
  A reduced run on such a file imposes the rates.

## Boundary accuracy

Several details of the boundary matter for how closely a reduced run follows
its parent, and are handled automatically:

* Boundary data is written at every time step of the parent, not only at
  report steps, and interpolated in time between records by the reduced run.
* Rates are time averages over the parent's step, taken from the converged
  solution, so the boundary moves exactly the mass the parent moved across
  it.
* The pressure on the far side is carried hydrostatically from the centre of
  the exterior cell to the face where the boundary condition is imposed.
* The fluid entering the sector is evaluated with the PVT region, relative
  permeabilities and capillary pressures of the cell it comes from, which
  need not be those of the interior cell.
* Faces that are connected through faults with throw, pinch-outs or other
  non-neighbour connections are part of the boundary. Cells that are
  inactive in the parent are not part of the sector, and a face towards one
  is a no-flow boundary, as it was in the parent.
* Threshold pressures (THPRES) across the boundary are honoured. Because a
  defaulted threshold pressure depends on the equilibration of the whole
  field, the reduced run uses the threshold pressures the parent arrived at.

## Wells, groups, UDQs and summary vectors outside the sector

A reduced run only simulates the wells inside its sector, but the deck's
group controls, user defined quantities (UDQs) and ACTIONX conditions are
normally written for the whole field. A FLUX file therefore also carries what
the reduced run needs about the rest of the field:

* **Summary vectors of the parent.** Rates and cumulatives of all wells, plus
  every summary vector the deck's UDQs and ACTIONX conditions refer to, are
  embedded at every time step of the parent.
* **Wells and groups of the parent.** The name, group, type, status and
  efficiency factor of every well, and the group tree, at every report step.

With these, a reduced run behaves as the parent did towards the rest of the
field:

* Wells outside the sector take part in group control at the rates they had
  in the parent. They count towards group and field totals and are
  subtracted from group targets, so the wells inside the sector are given the
  share of a group target they had in the parent.
* Field and group summary vectors, such as FOPR and FOPT, cover the whole
  field and are directly comparable with the parent's.
* UDQs, UDAs and ACTIONX conditions are evaluated with the parent's values
  for every quantity the sector cannot compute itself, such as the rates of
  wells outside it or the pressure of a region it does not contain.
* Rate conversion for wells on reservoir volume control, and the FPRFLUX
  summary vector, average over the whole parent model rather than over the
  sector alone.

When the reduced run uses the parent deck, the wells outside the sector are
already in its schedule; they simply have no connections in the sector. A
cut-out sector deck usually keeps only its own wells, and the others are then
added to its schedule from the FLUX file, without connections.

## FLUXALL

Which summary vectors a parent has to report for the reduced run depends on
the deck: on the wells, and on what the UDQs and ACTIONX conditions refer
to. FLUXALL in the parent's SUMMARY section asks for exactly those vectors,
so that they do not have to be listed by hand. It is required for make_flux,
which can only embed what the parent wrote to its summary file, and
recommended for DUMPFLUX, which embeds the same vectors from the run's own
summary state.

## make_flux

make_flux builds a FLUX file from the output of a parent run that was made
without DUMPFLUX. It reads the parent's deck and its EGRID, INIT, restart and
summary files. The sector can be selected with a FLUXNUM array in a separate
file, or with a box or a list of cells given on the command line. The FLUX
file contains a boundary record at every restart step of the parent, and the
embedded summary at every step of the parent's summary file. In pressure
mode, the exterior relative permeabilities and capillary pressures are
evaluated from the parent's saturation functions. Flux mode needs the
FLORES and FLOWS restart arrays.

Since a restart is a snapshot, the boundary of a file built by make_flux is
only as detailed as the parent's restart output. Requesting a restart at
every report step, as shown in the [worked
examples](02-examples.md#example-2-make_flux-and-a-cut-out-sector-deck),
brings a reduced run close to what a DUMPFLUX file gives.

## Other features

* **Parallel runs.** Both DUMPFLUX and USEFLUX runs can be run in parallel.
  A FLUX file written in parallel is the same as one written serially, and a
  reduced run gives the same results whatever the number of processes.
* **inspect_flux.** A command-line tool that prints the boundary data a FLUX
  file holds for one cell, or the embedded summary vectors, as a time series.
* **Run-time parameters** control how often DUMPFLUX writes boundary records
  and summary samples, and whether the reduced run interpolates the boundary
  state in time. They are listed in the [technical
  documentation](03-technical.md#run-time-parameters).
