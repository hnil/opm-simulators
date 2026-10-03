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

#ifndef OPM_FLUX_PARENT_WELLS_HPP
#define OPM_FLUX_PARENT_WELLS_HPP

#include <opm/io/eclipse/FluxFile.hpp>

#include <string>
#include <vector>

namespace Opm {

class ErrorGuard;
class ParseContext;
class Schedule;

/// The producing run's wells and groups, for the FLUX file.
///
/// Shared by DUMPFLUX and make_flux, so that a file built either way tells a
/// reduced run the same thing about the wells it does not simulate.
EclIO::FluxFile::ParentWells describeFluxParentWells(const Schedule& schedule);

/// Give a reduced run's schedule the producer's wells it does not define.
///
/// A sector cut out of the full model normally keeps only its own wells, and
/// perhaps only the groups they are in. The wells outside it still count
/// towards every group and field quantity, every group target and every UDQ,
/// and the FLUX file carries their rates, but a well the schedule does not
/// have cannot take part in any of that. Each missing well is therefore added
/// here, in its parent group, as a producer or injector as in the parent and
/// with the parent's status and efficiency factor, but without connections.
/// From then on it is exactly what a well outside the sector is in a reduced
/// run on the parent deck: its rates come from the parent's summary, in the
/// summary output and in group control alike. Groups the reduced run lacks
/// are added likewise.
///
/// Must be called before the SummaryConfig is built from the schedule, so
/// that the added wells and groups are reported too, and on a schedule that
/// kept its keywords, since inserting keywords replays the rest of it. The
/// replay is held to \p parseContext, as the original reading was.
///
/// \return Names of the wells added.
std::vector<std::string>
addFluxParentWells(Schedule& schedule,
                   const EclIO::FluxFile::ParentWells& parent,
                   const ParseContext& parseContext,
                   ErrorGuard& errors);

} // namespace Opm

#endif // OPM_FLUX_PARENT_WELLS_HPP
