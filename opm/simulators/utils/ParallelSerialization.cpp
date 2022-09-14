/*
  Copyright 2020 Equinor AS.

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

#include <config.h>

#include <opm/simulators/utils/ParallelSerialization.hpp>

#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/Schedule/Schedule.hpp>
#include <opm/input/eclipse/Schedule/Action/ASTNode.hpp>
#include <opm/input/eclipse/Schedule/MSW/SICD.hpp>
#include <opm/input/eclipse/Schedule/MSW/Valve.hpp>
#include <opm/input/eclipse/Schedule/Action/State.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQState.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQASTNode.hpp>
#include <opm/input/eclipse/Schedule/UDQ/UDQConfig.hpp>
#include <opm/input/eclipse/Schedule/Well/WList.hpp>
#include <opm/input/eclipse/Schedule/Well/WListManager.hpp>
#include <opm/input/eclipse/Schedule/Well/WellTestState.hpp>
#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/input/eclipse/EclipseState/Grid/TransMult.hpp>

#include <ebos/eclmpiserializer.hh>

#include <dune/common/parallel/mpihelper.hh>

namespace Opm {


void eclStateBroadcast(Parallel::Communication comm, EclipseState& eclState, Schedule& schedule,
                       SummaryConfig& summaryConfig,
                       UDQState& udqState,
                       Action::State& actionState,
                       WellTestState&  wtestState)
{
    Opm::EclMpiSerializer ser(comm);
    ser.broadcast(eclState);
    ser.broadcast(schedule);
    ser.broadcast(summaryConfig);
    ser.broadcast(udqState);
    ser.broadcast(actionState);
    ser.broadcast(wtestState);
}

template <class T>
void eclBroadcast(Parallel::Communication comm, T& data)
{
    Opm::EclMpiSerializer ser(comm);
    ser.broadcast(data);
}


template void eclBroadcast<TransMult>(Parallel::Communication, TransMult&);
template void eclBroadcast<Schedule>(Parallel::Communication, Schedule&);

}
