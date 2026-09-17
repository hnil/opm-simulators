/*
  Copyright 2026 SINTEF Digital

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

#ifndef OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED
#define OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED

#include <opm/simulators/wells/ProdGroupTreeNode.hpp>

#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm {

class DeferredLogger;
class GuideRate;
class SummaryState;
template<typename Scalar, typename IndexTraits> class BlackoilWellModelGeneric;

} // namespace Opm

namespace Opm::ProdGroupTreeBalancer {

/// Type alias for the tree map.
template<class Scalar>
using Tree = std::map<std::string, ProdGroupTreeNode<Scalar>>;

/// Top-level entry point: build tree, balance it, validate, and apply.
/// All internal functions are implementation details not exposed through this interface.
///
/// \param[in]    wellModel     Well model (schedule, well/group state, guide rate, RESV coefficients)
/// \param[in]    summaryState  Summary state
/// \param[in]    reportStep    Current report step index
/// \param[in]    tol           Convergence tolerance
/// \param[in]    limits        Globally gathered well limits (from prepareWellsForBalancing_*)
/// \param[in]    logger        Deferred logger
/// \param[in]    assignTargets Also write the allocation as each well's group target
///                             (setTargets); otherwise legacy derives the targets from the modes
/// \param[out]   decidedWells  If given, the wells the tree assigned a control to
/// \return       true if the result passed checkTreeValidity
template<class Scalar, typename IndexTraits>
bool runGroupTreeBalancer(BlackoilWellModelGeneric<Scalar, IndexTraits>& wellModel,
                          const SummaryState& summaryState,
                          int reportStep,
                          Scalar tol,
                          const std::unordered_map<std::string, std::pair<int, Scalar>>& limits,
                          DeferredLogger& logger,
                          bool assignTargets = false,
                          std::vector<std::string>* decidedWells = nullptr,
                          bool writeRates = true);

/// Balance a tree built by hand (unit tests): the same top-down pass
/// buildTree() ends with, then the algorithm and the validity check.
template<class Scalar>
bool balanceTreeForTesting(Tree<Scalar>& tree,
                           const GuideRate& guideRate,
                           Scalar tol,
                           DeferredLogger& logger,
                           bool assignTargets = false);

} // namespace Opm::ProdGroupTreeBalancer

#endif // OPM_PROD_GROUP_TREE_BALANCER_HEADER_INCLUDED
