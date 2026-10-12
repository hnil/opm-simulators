# Limitations {#limitations}

* **Black-oil models only.** DUMPFLUX and USEFLUX work with the black-oil
  simulators. A compositional run that would write a FLUX file stops with an
  error. Water-borne gas (Rsw) and vaporised water (Rvw), solvent, polymer,
  brine, energy and tracers are not carried across the boundary: the fluid
  entering the sector has no solvent, polymer, salt or tracer, and its
  temperature is the recorded one only for thermal runs.

* **One sector per reduced run.** A FLUX file describes one FLUXNUM region,
  and a reduced run covers one region. FLUXREG accepts a list of regions,
  but the boundary data only fits if together they make up the region the
  FLUX file describes.

* **Wells must lie within one region.** A DUMPFLUX run stops if a well has
  connections in more than one FLUXNUM region. A reduced run would only see
  the part of such a well that lies in its sector.

* **No feedback to the outside.** The boundary data is the parent's. Changes
  made inside the sector do not change the state outside it, so a change
  that would have affected the rest of the field, and through it the flow
  back into the sector, is only accounted for by the pressure boundary's
  response to the local pressure difference. The larger the change, the
  further from the boundary it should be.

* **Wells outside the sector follow the parent.** Their rates are the
  parent's, whatever happens inside the sector. Group controls are applied
  to the wells inside the sector with the outside wells' parent rates held
  fixed. ACTIONX conditions and UDQs see the parent's values for everything
  the sector cannot compute; an ACTIONX that would have acted on a well
  outside the sector changes that well in the schedule, but not its rates.

* **Same grid geometry.** A reduced run must use either the parent grid, or
  a grid of the dimensions of the region's bounding box with the same cells.
  Refining or coarsening the sector, or a sector grid that is not the
  bounding box, is not supported.

* **make_flux is limited by the parent's output.** Boundary records exist
  only where the parent wrote a restart, and the boundary is interpolated in
  time between them. Exterior relative permeabilities and capillary
  pressures are evaluated from saturation functions, without hysteresis.
  Flux mode needs FLOWS in the restart files and a mapping rather than a
  FLUXNUM file. Only summary vectors the parent reported can be embedded,
  which in practice means the parent needs FLUXALL.

* **Restarts and local grid refinements.** Restarted runs and models with
  local grid refinements have not been tested with DUMPFLUX or USEFLUX.

* **File size.** A DUMPFLUX file grows with the number of boundary faces and
  time steps. Drogon's 2427 faces over 157 time steps take 34 MB in pressure
  mode. `--flux-boundary-min-interval-between-samples` reduces the number of
  records.
