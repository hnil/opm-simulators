/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2017 IRIS
  Copyright 2019 Norce

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

#ifndef OPM_WELLINTERFACE_INDICES_HEADER_INCLUDED
#define OPM_WELLINTERFACE_INDICES_HEADER_INCLUDED

#include <opm/simulators/wells/WellInterfaceFluidSystem.hpp>
#include <opm/simulators/wells/WellInterfaceEval.hpp>

namespace Opm
{

template<class FluidSystem, class Indices, class Scalar>
class WellInterfaceIndices : public WellInterfaceFluidSystem<FluidSystem>
                           , public WellInterfaceEval<FluidSystem>
{
public:
    using WellInterfaceFluidSystem<FluidSystem>::Gas;
    using WellInterfaceFluidSystem<FluidSystem>::Oil;
    using WellInterfaceFluidSystem<FluidSystem>::Water;

protected:
    WellInterfaceIndices(const Well& well,
                         const ParallelWellInfo& parallel_well_info,
                         const int time_step,
                         const typename WellInterfaceFluidSystem<FluidSystem>::RateConverterType& rate_converter,
                         const int pvtRegionIdx,
                         const int num_components,
                         const int num_phases,
                         const int index_of_well,
                         const int first_perf_index,
                         const std::vector<PerforationData>& perf_data);

    int flowPhaseToEbosCompIdx( const int phaseIdx ) const;
    int ebosCompIdxToFlowCompIdx( const unsigned compIdx ) const;

    double scalingFactor(const int phaseIdx) const;
};

}

#endif // OPM_WELLINTERFACE_INDICES_HEADER_INCLUDED
