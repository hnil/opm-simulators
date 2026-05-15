/*
  Copyright 2024, SINTEF Digital

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
#ifndef FLOW_COMP_HPP
#define FLOW_COMP_HPP

#include <memory>

#include <opm/material/constraintsolvers/PTFlash.hpp>
#include <opm/material/fluidsystems/GenericOilGasWaterFluidSystem.hpp>

#include <opm/models/discretization/common/baseauxiliarymodule.hh>
#include <opm/models/nonlinear/newtonmethod.hh>
#include <opm/models/ptflash/flashmodel.hh>

#include <opm/simulators/flow/Main.hpp>
#include <opm/simulators/flow/FlowProblemComp.hpp>
#include <opm/simulators/flow/FlowProblemCompProperties.hpp>

#include <opm/simulators/linalg/parallelbicgstabbackend.hh>

#include <flowexperimental/comp/EmptyModel.hpp>
#include <flowexperimental/comp/wells/CompWellModel.hpp>

// // the current code use eclnewtonmethod adding other conditions to proceed_ should do the trick for KA
// // adding linearshe sould be chaning the update_ function in the same class with condition that the error is reduced.
// the trick is to be able to recalculate the residual from here.
// unsure where the timestepping is done from suggestedtime??
// suggestTimeStep is taken from newton solver in problem.limitTimestep
namespace Opm {

template<int numComp, bool EnableWater>
int dispatchFlowComp(int argc, char** argv);

template <class TypeTag>
class FlowCompNewtonMethod : public Opm::NewtonMethod<TypeTag>
{
public:
    using ParentType = Opm::NewtonMethod<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    enum { pressure0Idx = Indices::pressure0Idx };
    enum { z0Idx = Indices::z0Idx };
    enum { numComponents = getPropValue<TypeTag, Properties::NumComponents>() };

    static constexpr bool waterEnabled = Indices::waterEnabled;

    explicit FlowCompNewtonMethod(Simulator& simulator)
        : ParentType(simulator)
    {}

    using ParentType::preSolve_;
    using ParentType::update_;

protected:
    friend ParentType;
    friend NewtonMethod<TypeTag>;

    /*!
     * \copydoc FvBaseNewtonMethod::updatePrimaryVariables_
     * Apply update limiters to prevent divergence after control changes
     */
    void updatePrimaryVariables_(unsigned /* globalDofIdx */,
                                 PrimaryVariables& nextValue,
                                 const PrimaryVariables& currentValue,
                                 const EqVector& update,
                                 const EqVector& /* currentResidual */)
    {
        // normal Newton-Raphson update
        nextValue = currentValue;
        nextValue -= update;

        ////
        // Pressure updates
        ////
        // limit pressure reference change relative to the total value per iteration
        constexpr Scalar max_percent_change = 0.2;
        constexpr Scalar upper_bound = 1. + max_percent_change;
        constexpr Scalar lower_bound = 1. - max_percent_change;
        nextValue[pressure0Idx] = std::clamp(nextValue[pressure0Idx],
                                             currentValue[pressure0Idx] * lower_bound,
                                             currentValue[pressure0Idx] * upper_bound);

        ////
        // z updates
        ////
        // restrict update
        Scalar maxDeltaZ = 0.0;  // in update vector
        Scalar sumDeltaZ = 0.0; // changes in last component (not in update vector)
        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
            maxDeltaZ = std::max(std::abs(update[z0Idx + compIdx]), maxDeltaZ);
            sumDeltaZ += update[z0Idx + compIdx];
        }
        maxDeltaZ = std::max(std::abs(sumDeltaZ), maxDeltaZ);

        // if max. update is above limit, restrict that one to limit and adjust the rest
        // accordingly (s.t. last comp. update is sum of the changes in update vector)
        constexpr Scalar deltaz_limit = 0.2;
        if (maxDeltaZ > deltaz_limit) {
            const Scalar alpha = deltaz_limit / maxDeltaZ;
            for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
                nextValue[z0Idx + compIdx] = currentValue[z0Idx + compIdx] - alpha * update[z0Idx + compIdx];
            }
        }

        // ensure that z-values are less than tol or more than 1-tol
        constexpr Scalar tol = 1e-8;
        for (unsigned compIdx = 0; compIdx < numComponents - 1; ++compIdx) {
           nextValue[z0Idx + compIdx] = std::clamp(nextValue[z0Idx + compIdx], tol, 1-tol);
        }

        if constexpr (waterEnabled) {
            // limit change in water saturation
            constexpr Scalar dSwMax = 0.2;
            if (update[Indices::water0Idx] > dSwMax) {
                nextValue[Indices::water0Idx] = currentValue[Indices::water0Idx] - dSwMax;
            }
        }
    }
};

}

namespace Opm::Properties {
namespace TTag {

template<int NumComp, bool EnableWater>
struct FlowCompProblem {
   using InheritsFrom = std::tuple<FlowBaseProblemComp, FlashModel>;
};

}

template<class TypeTag, int NumComp, bool EnableWater>
struct SparseMatrixAdapter<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    enum { numEq = getPropValue<TypeTag, Properties::NumEq>() };
    using Block = MatrixBlock<Scalar, numEq, numEq>;

public:
    using type = typename Linear::IstlSparseMatrixAdapter<Block>;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct NewtonMethod<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>>
{
    using type = Opm::FlowCompNewtonMethod<TypeTag>;
};

#if 0
template<class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::FlowCompProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager = ::Opm::EclThermalLawManager<Scalar, FluidSystem>;

    using type = typename EclThermalLawManager::SolidEnergyLaw;
};

// Set the material law for thermal conduction
template<class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::FlowCompProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    using EclThermalLawManager = ::Opm::EclThermalLawManager<Scalar, FluidSystem>;

    using type = typename EclThermalLawManager::ThermalConductionLaw;
};


template <class TypeTag>
struct SpatialDiscretizationSplice<TypeTag, TTag::FlowCompProblem>
{
    using type = TTag::EcfvDiscretization;
};

template <class TypeTag>
struct LocalLinearizerSplice<TypeTag, TTag::FlowCompProblem>
{
    using type = TTag::AutoDiffLocalLinearizer;
};
#endif

// Set the problem property
template <class TypeTag, int NumComp, bool EnableWater>
struct Problem<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>>
{
    using type = FlowProblemComp<TypeTag>;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct AquiferModel<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    using type = EmptyModel<TypeTag>;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct WellModel<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    using type = CompWellModel<TypeTag>;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct TracerModel<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    using type = EmptyModel<TypeTag>;
};


template <class TypeTag, int NumComp, bool EnableWater>
struct FlashSolver<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;

public:
    using type = Opm::PTFlash<Scalar, FluidSystem>;
};


template <class TypeTag, class MyTypeTag>
struct NumComp { using type = UndefinedProperty; };

// TODO: this is unfortunate, have to check why we need to hard-code it
template <class TypeTag, int NumComp_, bool EnableWater_>
struct NumComp<TypeTag, TTag::FlowCompProblem<NumComp_, EnableWater_>> {
    static constexpr int value = NumComp_;
};

template <class TypeTag, class MyTypeTag>
struct EnableDummyWater { using type = UndefinedProperty; };

template <class TypeTag, int NumComp_, bool EnableWater_>
struct EnableDummyWater<TypeTag, TTag::FlowCompProblem<NumComp_, EnableWater_>> {
    static constexpr bool value = EnableWater_;
};
#if 0
struct Temperature { using type = UndefinedProperty; };

 template <class TypeTag>
 struct Temperature<TypeTag, TTag::FlowCompProblem> {
     using type = GetPropType<TypeTag, Scalar>;
     static constexpr type value = 423.25;
 };
#endif

template <class TypeTag, int NumComp_, bool EnableWater_>
struct FluidSystem<TypeTag, TTag::FlowCompProblem<NumComp_, EnableWater_>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    static constexpr int num_comp = getPropValue<TypeTag, Properties::NumComp>();
    static constexpr bool enable_water = getPropValue<TypeTag, Properties::EnableDummyWater>();

public:
    using type = Opm::GenericOilGasWaterFluidSystem<Scalar, num_comp, enable_water>;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableMech<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableDisgasInWater<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct Stencil<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;

public:
    using type = EcfvStencil<Scalar, GridView>;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableApiTracking<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableSaltPrecipitation<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnablePolymerMW<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnablePolymer<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableDispersion<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableBrine<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableVapwat<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

template<class TypeTag, int NumComp, bool EnableWater>
struct EnableSolvent<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableFoam<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableExtbo<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};
template<class TypeTag, int NumComp, bool EnableWater>
struct EnableBioeffects<TypeTag, TTag::FlowCompProblem<NumComp, EnableWater>> {
    static constexpr bool value = false;
};

// disable thermal flux boundaries by default
#if 0
template<class TypeTag>
struct EnableThermalFluxBoundaries<TypeTag, TTag::FlowCompProblem> {
    static constexpr bool value = false;
};
#endif

} // namespace Opm::Properties

namespace Opm {

template<int numComp, bool EnableWater>
int dispatchFlowComp(int argc, char** argv)
{
    using TypeTag = Properties::TTag::FlowCompProblem<numComp, EnableWater>;

    auto mainObject = std::make_unique<Opm::Main>(argc, argv);
    const auto ret = mainObject->runStatic<TypeTag>();
    mainObject.reset();
    return ret;
}

} // namespace Opm

#endif
