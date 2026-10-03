/*
  Copyright 2026 Equinor ASA

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

#include <opm/simulators/flow/flux/FluxActivation.hpp>

#include <opm/common/ErrorMacros.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/IOConfig/IOConfig.hpp>
#include <opm/simulators/flow/flux/FluxRegions.hpp>

namespace Opm {

bool applyUseFluxActnum(EclipseState& eclipseState)
{
    const auto& ioConfig = eclipseState.getIOConfig();
    if (!ioConfig.getUseFlux()) {
        return false;
    }

    const auto& fieldProps = eclipseState.globalFieldProps();

    if (!fieldProps.has_int("FLUXNUM")) {
        OPM_THROW(std::invalid_argument,
                  "applyUseFluxActnum(): USEFLUX requires FLUXNUM in the sector deck");
    }

    const auto regionValues = fieldProps.get_global_int("FLUXNUM");

    const auto& selectedRegions = ioConfig.getUseFluxRegions();
    auto actnum = selectedRegions.empty()
        ? FluxRegions::buildActnum(regionValues, FluxRegions::uniqueSelectedRegion(regionValues))
        : FluxRegions::buildActnum(regionValues, selectedRegions);

    const auto& currentActnum = eclipseState.globalFieldProps().actnumRaw();
    for (std::size_t index = 0; index < actnum.size(); ++index) {
        actnum[index] = actnum[index] && currentActnum[index];
    }

    eclipseState.reset_actnum(actnum);
    return true;
}

} // namespace Opm
