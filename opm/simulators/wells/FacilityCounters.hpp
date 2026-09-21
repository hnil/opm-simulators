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

#include <map>
#include <string>

namespace Opm {

// What the facility solve spends, whoever drives it. Diagnostic, serial only;
// logged with OPM_FACILITY_CHECK.
struct FacilityCounters
{
    long outer_iterations = 0;     // control + network passes before a linearisation
    long network_sub_iterations = 0;
    long well_solves = 0, well_solves_failed = 0;
    long well_linearisations = 0;  // one per inner well iteration
    long ipr_assemblies = 0;       // of those, the implicit IPR's (not a solve iteration)
    /// OPM_WELL_SOLVE_STATS=1: solves split by what changed since the well's previous solve.
    /// first = first in this Newton iteration; same = same control and target; changed = not.
    long first_solves = 0, first_lins = 0, same_solves = 0, same_lins = 0, same_cold = 0;
    long changed_solves = 0, changed_lins = 0;
    struct LastSolve { double time = -1; int iteration = -1; int cmode = -1; double target = 0; };
    std::map<std::string, LastSolve> last;

    static FacilityCounters& get()
    {
        static FacilityCounters counters;
        return counters;
    }
};

} // namespace Opm

#endif
