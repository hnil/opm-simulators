/*
  Copyright 2026 Equinor ASA.

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef OPM_FACILITY_COUNTERS_HEADER_INCLUDED
#define OPM_FACILITY_COUNTERS_HEADER_INCLUDED

namespace Opm {

// What the facility solve spends, whoever drives it. Diagnostic, serial only;
// logged with OPM_FACILITY_CHECK.
struct FacilityCounters
{
    long outer_iterations = 0;     // control + network passes before a linearisation
    long network_sub_iterations = 0;
    long well_solves = 0, well_solves_failed = 0;
    long well_linearisations = 0;  // one per inner well iteration

    static FacilityCounters& get()
    {
        static FacilityCounters counters;
        return counters;
    }
};

} // namespace Opm

#endif
