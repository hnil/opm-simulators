/*
  Copyright 2022 Equinor ASA.

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

#ifndef OPM_GASLIFT_COMMON_HEADER_INCLUDED
#define OPM_GASLIFT_COMMON_HEADER_INCLUDED

#include <opm/simulators/utils/DeferredLogger.hpp>
#include <opm/simulators/wells/WellState.hpp>

#include <string>
#include <fmt/format.h>

namespace Opm
{
class GasLiftCommon
{
public:
    virtual ~GasLiftCommon() = default;

protected:
    GasLiftCommon(
        WellState &well_state,
        DeferredLogger &deferred_logger,
        bool debug
    );
    int debugUpdateGlobalCounter_() const;
    virtual void displayDebugMessage_(const std::string& msg) const = 0;

    WellState &well_state_;
    DeferredLogger &deferred_logger_;
    bool debug;
};

} // namespace Opm

#endif // OPM_GASLIFT_COMMON_INCLUDED
