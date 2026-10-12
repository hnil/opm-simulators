---
title: Sector modelling in OPM Flow
subtitle: DUMPFLUX, USEFLUX and make_flux
---

# Sector modelling in OPM Flow {#top}

This guide describes how to run a simulation on part of a reservoir model,
with boundary conditions taken from a run of the whole model.

1. [Introduction](01-introduction.md#introduction): features, keywords and
   workflow.
2. [Worked examples](02-examples.md#worked-examples) on the Drogon model:
   DUMPFLUX and USEFLUX on the parent deck, make_flux and a cut-out sector
   deck, and UDQs over wells outside the sector.
3. [Technical documentation](03-technical.md#technical-documentation):
   boundary conditions, the FLUX file format, wells, groups and UDQs outside
   the sector.
4. [Limitations](04-limitations.md#limitations).

The markdown files in this folder are the source. `make sector_docs` in the
opm-simulators build directory, or `build-docs.sh` here, turns them into a
single HTML page and a PDF, and makes a PowerPoint presentation for users
with `tools/make-slides.py`.
