/*
  Copyright 2020, NORCE AS

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
#include "config.h"
#include <opm/simulators/flow/Main.hpp>
#include <opm/grid/polyhedralgrid.hh>
#include <ebos/eclpolyhedralgridvanguard.hh>
#include <ebos/equil/initstateequil_impl.hh>
#include <opm/models/blackoil/blackoillocalresidualtpfa.hh>
#include <opm/models/discretization/common/tpfalinearizer.hh>

namespace Opm {
namespace Properties {
    namespace TTag {
        struct EclFlowProblemPoly {
            using InheritsFrom = std::tuple<EclFlowProblem>;
        };
    }

    template<class TypeTag>
    struct Linearizer<TypeTag, TTag::EclFlowProblemPoly> { using type = TpfaLinearizer<TypeTag>; };

    template<class TypeTag>
    struct LocalResidual<TypeTag, TTag::EclFlowProblemPoly> { using type = BlackOilLocalResidualTPFA<TypeTag>; };

    template<class TypeTag>
    struct EnableDiffusion<TypeTag, TTag::EclFlowProblemPoly> { static constexpr bool value = false; };

    template<class TypeTag>
    struct Grid<TypeTag, TTag::EclFlowProblemPoly> {
        using type = Dune::PolyhedralGrid<3, 3>;
    };
    template<class TypeTag>
    struct EquilGrid<TypeTag, TTag::EclFlowProblemPoly> {
        //using type = Dune::CpGrid;
        using type = GetPropType<TypeTag, Properties::Grid>;
    };

    template<class TypeTag>
    struct Vanguard<TypeTag, TTag::EclFlowProblemPoly> {
        using type = Opm::EclPolyhedralGridVanguard<TypeTag>;
    };
}
}
int main(int argc, char** argv)
{
    using TypeTag = Opm::Properties::TTag::EclFlowProblemPoly;
    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    auto ret = mainObject->runStatic<TypeTag>();
    // Destruct mainObject as the destructor calls MPI_Finalize!
    mainObject.reset();
    return ret;
}
