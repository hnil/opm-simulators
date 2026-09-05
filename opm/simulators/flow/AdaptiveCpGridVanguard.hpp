// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 * \copydoc Opm::AdaptiveCpGridVanguard
 */
#ifndef OPM_ADAPTIVE_CPGRID_VANGUARD_HPP
#define OPM_ADAPTIVE_CPGRID_VANGUARD_HPP

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/models/utils/parametersystem.hpp>

#include <opm/simulators/flow/AdaptiveLgr.hpp>
#include <opm/simulators/flow/WellZoneLgr.hpp>
#include <opm/simulators/flow/CpGridVanguard.hpp>

// A grid backend may require a refinement builder to be registered before a
// post-construction refinement request.  opm-gridrefined does, for a deck
// without CARFIN: it only retains the corner-point description when the deck
// declares LGRs, so its own fallback has nothing to build from.  Upstream
// opm-grid refines from the grid it already holds and has no such notion.
#if __has_include(<opm/grid/cpgrid/refinement/RefinementBuilder.hpp>)
#define OPM_ADAPTIVE_HAVE_REFINEMENT_BUILDER 1
#include <opm/grid/cpgrid/refinement/ConformingBlockBuilder.hpp>
#include <opm/grid/cpgrid/refinement/RefinementBuilder.hpp>
#endif

#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Opm {

/*!
 * \brief A CpGridVanguard that additionally performs adaptive local refinement
 *        on a deck that has no CARFIN, driven by the --adaptive-lgr parameter.
 *
 * This is the AdaptiveCpGrid demonstration as a *separate Vanguard* (used only by
 * the flow_blackoil_adaptive executable), so the standard CpGridVanguard and the
 * standard flow binaries are completely unchanged. After the coarse grid is
 * built, addLgrs() refines the requested region by reusing the exact same static
 * refinement machinery (CpGrid::addLgrsUpdateLeafView with the conforming
 * builder), so the resulting grid -- and the simulation -- match the equivalent
 * static-CARFIN run. Serial demonstration; the builder is constructed from the
 * input grid's (sanitized-at-construction) COORD/ZCORN.
 */
template <class TypeTag>
class AdaptiveCpGridVanguard : public CpGridVanguard<TypeTag>
{
    using Base = CpGridVanguard<TypeTag>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;

public:
    explicit AdaptiveCpGridVanguard(Simulator& simulator)
        : Base(simulator)
    {}

    static void registerParameters()
    {
        Base::registerParameters();
        Parameters::Register<Parameters::AdaptiveLgr>(
            "Adaptive local grid refinement on a deck with no CARFIN: one or more "
            "CARFIN-style boxes 'I1 I2 J1 J2 K1 K2 NX NY NZ' (1-based, "
            "';'-separated). Empty (default) => behaves like flow_blackoil.");
        Parameters::Register<Parameters::WellRefine>(
            "Nested refinement rings around wells, ';'-separated zones "
            "PATTERN:FX,FY,FZ:rings=W1[,W2,...][:layers=all][:klayers=N]. Each ring "
            "is refined FX x FY x FZ relative to the ring outside it; the K range "
            "covers the perforated layers (plus N above and below) or every layer.");
    }

    //! \brief The refinement spec to apply; overridable so a dynamic driver
    //! can substitute a new mark set on a mid-run rebuild.
    virtual std::string adaptiveLgrSpec() const
    { return Parameters::Get<Parameters::AdaptiveLgr>(); }

    //! \brief Add deck LGRs (base) then any --adaptive-lgr refinement.
    void addLgrs()
    {
        // Honour a deck CARFIN if present (an adaptive deck normally has none).
        Base::addLgrs();

        const std::string spec = this->adaptiveLgrSpec();
        const std::string wellSpec = Parameters::Get<Parameters::WellRefine>();
        if (spec.empty() && wellSpec.empty()) {
            return;
        }
        const auto& cartDims = this->grid_->logicalCartesianSize();

        // Explicit boxes and well zones go into one request list: the builder
        // refines once and checks the boxes against each other.
        std::vector<Opm::Refinement::BlockRefinement> requests;
        for (const auto& b : parseAdaptiveLgrSpec(spec)) {
            Opm::Refinement::BlockRefinement req;
            req.name = b.name;
            req.cellsPerDim = b.cellsPerDim;
            req.startIJK = b.startIJK;
            req.endIJK = b.endIJK;
            requests.push_back(std::move(req));
        }
        const auto zoneRequests = wellZoneRefinements(parseWellZoneSpec(wellSpec),
                                                     this->schedule(), cartDims);
        requests.insert(requests.end(), zoneRequests.begin(), zoneRequests.end());
        if (requests.empty()) {
            return;
        }

        const auto& inputGrid = this->eclState().getInputGrid();

        std::vector<LgrCellBox> cellBoxes;
        for (const auto& r : requests) {
            if (r.parentGridName == "GLOBAL") {
                cellBoxes.push_back({ r.startIJK, r.endIJK });
            }
            OpmLog::info(fmt::format("  refinement {} on {}: cells [{},{}]x[{},{}]x[{},{}] by {}x{}x{}",
                                     r.name, r.parentGridName,
                                     r.startIJK[0] + 1, r.endIJK[0], r.startIJK[1] + 1, r.endIJK[1],
                                     r.startIJK[2] + 1, r.endIJK[2],
                                     r.cellsPerDim[0], r.cellsPerDim[1], r.cellsPerDim[2]));
        }
        refuseDeckConnectionsInsideBoxes(this->eclState(), cartDims, cellBoxes);
        OpmLog::info("\nAdaptive refinement: refining "
                     + std::to_string(requests.size())
                     + " box(es) on the coarse grid post-construction");

        // Refine through the grid's own entry point, exactly as a deck CARFIN
        // does.
#if OPM_ADAPTIVE_HAVE_REFINEMENT_BUILDER
        // This backend needs a builder supplied for a deck with no CARFIN,
        // since it retains the corner-point description only when the deck
        // declares LGRs.  Build one from the input grid for the duration of
        // the call and restore whatever was registered before.
        std::vector<int> actnum = inputGrid.getACTNUM();
        auto previous = Opm::Refinement::setBuilder(
            std::make_unique<Opm::Refinement::ConformingBlockBuilder>(
                inputGrid.getNXYZ(), inputGrid.getCOORD(), inputGrid.getZCORN(),
                std::move(actnum)));
        try {
            this->grid_->addLgrsUpdateLeafView(requests);
        }
        catch (...) {
            Opm::Refinement::setBuilder(std::move(previous));
            throw;
        }
        Opm::Refinement::setBuilder(std::move(previous));
#else
        this->grid_->addLgrsUpdateLeafView(requests);
#endif

        // Same post-refinement bookkeeping the base addLgrs() does for a deck
        // CARFIN: refresh the leaf view, rebuild the stale Cartesian->compressed
        // map (so wells map to the right leaf cell), refresh depths/thickness.
        this->updateGridView_();
        this->updateCartesianToCompressedMapping_();
        this->updateCellDepths_();
        this->updateCellThickness_();

        // Wells in the refined region (plan S6b): COMPDAT wells carry no
        // trajectory to replay, so synthesize an equivalent one from their
        // connection cells (coarse input-grid geometry -- the polyline is
        // purely geometric) and then re-derive every trajectory well's
        // connections against the refined leaf grid.
        if (this->grid_->maxLevel() > 0) {
            this->schedule().synthesizeWellTrajectories(
                [&inputGrid](std::size_t globalIdx)
                { return inputGrid.getCellCenter(globalIdx); },
                [&inputGrid](std::size_t globalIdx)
                { return inputGrid.getCellDims(globalIdx); });

            this->recomputeWellTrajectoriesInLgr_();
        }
    }
};

} // namespace Opm

#endif // OPM_ADAPTIVE_CPGRID_VANGUARD_HPP
