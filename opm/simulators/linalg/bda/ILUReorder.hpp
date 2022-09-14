/*
  Copyright 2020 Equinor ASA

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

#ifndef ILUREORDER_HEADER_INCLUDED
#define ILUREORDER_HEADER_INCLUDED

namespace Opm
{
namespace Accelerator
{
    // Level Scheduling respects the dependencies in the original matrix, and behaves like Dune and cusparse
    // Graph Coloring is more aggresive and is likely to increase the number of linearizations and linear iterations to converge significantly, but can still be faster on GPU because it results in more parallelism
        
    enum class ILUReorder {
        LEVEL_SCHEDULING,
        GRAPH_COLORING,
        NONE
    };

} // namespace Accelerator
} // namespace Opm

#endif
