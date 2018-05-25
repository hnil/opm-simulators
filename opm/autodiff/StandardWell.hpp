/*
  Copyright 2017 SINTEF Digital, Mathematics and Cybernetics.
  Copyright 2017 Statoil ASA.
  Copyright 2016 - 2017 IRIS AS.

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


#ifndef OPM_STANDARDWELL_HEADER_INCLUDED
#define OPM_STANDARDWELL_HEADER_INCLUDED


#include <opm/autodiff/WellInterface.hpp>
#include <opm/autodiff/ISTLSolver.hpp>
#include <opm/autodiff/RateConverter.hpp>
#include <dune/istl/matrixmarket.hh>
#include <opm/autodiff/ISTLSolver.hpp>

namespace Opm
{

    template<typename TypeTag>
    class StandardWell: public WellInterface<TypeTag>
    {

    public:
        typedef WellInterface<TypeTag> Base;
        // TODO: some functions working with AD variables handles only with values (double) without
        // dealing with derivatives. It can be beneficial to make functions can work with either AD or scalar value.
        // And also, it can also be beneficial to make these functions hanle different types of AD variables.
        using typename Base::Simulator;
        using typename Base::WellState;
        using typename Base::IntensiveQuantities;
        using typename Base::FluidSystem;
        using typename Base::MaterialLaw;
        using typename Base::ModelParameters;
        using typename Base::Indices;
        using typename Base::PolymerModule;
        using typename Base::RateConverterType;

        using Base::numEq;

        // the positions of the primary variables for StandardWell
        // there are three primary variables, the second and the third ones are F_w and F_g
        // the first one can be total rate (G_t) or bhp, based on the control

        static const bool gasoil = numEq == 2 && (Indices::compositionSwitchIdx >= 0);
        static const int XvarWell = 0;
        static const int WFrac = gasoil? -1000: 1;
        static const int GFrac = gasoil? 1: 2;
        static const int SFrac = 3;


        using typename Base::Scalar;
        using typename Base::ConvergenceReport;


        using Base::has_solvent;
        using Base::has_polymer;
        using Base::name;
        using Base::Water;
        using Base::Oil;
        using Base::Gas;

        // TODO: with flow_ebos，for a 2P deck, // TODO: for the 2p deck, numEq will be 3, a dummy phase is already added from the reservoir side.
        // it will cause problem here without processing the dummy phase.
        //static const int numWellEq = GET_PROP_VALUE(TypeTag, EnablePolymer)? numEq-1 : numEq; // number of wellEq is only numEq - 1 for polymer
        static const int numWellEq =  numEq;
        static const int numAdjoint = GET_PROP_VALUE(TypeTag, numAdjoint);
        static const int control_index=numEq + numWellEq;
        using typename Base::Mat;
        using typename Base::BVector;
        using typename Base::Eval;

        // sparsity pattern for the matrices
        //[A C^T    [x       =  [ res
        // B  D ]   x_well]      res_well]

        // the vector type for the res_well and x_well
        typedef Dune::FieldVector<Scalar, numWellEq> VectorBlockWellType;
        typedef Dune::BlockVector<VectorBlockWellType> BVectorWell;

        // vector type or cntrl need for adjoint
        typedef Dune::FieldVector<Scalar, 1> VectorBlockWellCtrlType;
        typedef Dune::BlockVector<VectorBlockWellCtrlType> BVectorWellCtrl;

        // vector type for res need for adjoint
        //typedef Dune::FieldVector<Scalar, numEq >        VectorBlockResType;
        //typedef Dune::BlockVector<VectorBlockResType>      BVectorRes;

        // the matrix type for the diagonal matrix D
        typedef Dune::FieldMatrix<Scalar, numWellEq, numWellEq > DiagMatrixBlockWellType;
        typedef Dune::FieldMatrix<Scalar, numWellEq, 1 > DiagMatrixBlockWellAdjointType;

        typedef Dune::BCRSMatrix <DiagMatrixBlockWellType> DiagMatWell;
        typedef Dune::BCRSMatrix <DiagMatrixBlockWellAdjointType> DiagMatWellCtrl;

        // the matrix type for the non-diagonal matrix B and C^T
        typedef Dune::FieldMatrix<Scalar, numWellEq, numEq>  OffDiagMatrixBlockWellType;
        typedef Dune::BCRSMatrix<OffDiagMatrixBlockWellType> OffDiagMatWell;

        // for adjoint
        typedef Dune::FieldMatrix<Scalar, 1, numEq>  OffDiagMatrixBlockWellAdjointType;
        typedef Dune::BCRSMatrix<OffDiagMatrixBlockWellAdjointType> OffDiagMatWellCtrl;

        // added extra space in derivative to have control derivatives
        typedef DenseAd::Evaluation<double, /*size=*/numEq + numWellEq+numAdjoint> EvalWell;

        using Base::contiSolventEqIdx;
        using Base::contiPolymerEqIdx;


        StandardWell(const Well* well, const int time_step, const Wells* wells,
                     const ModelParameters& param,
                     const RateConverterType& rate_converter,
                     const int pvtRegionIdx,
                     const int num_components);

        virtual void init(const PhaseUsage* phase_usage_arg,
                          const std::vector<double>& depth_arg,
                          const double gravity_arg,
                          const int num_cells);


        virtual void initPrimaryVariablesEvaluation() const;

        virtual void assembleWellEq(Simulator& ebosSimulator,
                                    const double dt,
                                    WellState& well_state,
                                    bool only_wells);

        /// updating the well state based the control mode specified with current
        // TODO: later will check wheter we need current
        virtual void updateWellStateWithTarget(WellState& well_state) const;

        /// check whether the well equations get converged for this well
        virtual ConvergenceReport getWellConvergence(const std::vector<double>& B_avg) const;

        /// Ax = Ax - C D^-1 B x
        virtual void apply(const BVector& x, BVector& Ax) const;
        /// r = r - C D^-1 Rw
        virtual void apply(BVector& r) const;

        /// Ax = Atx - Bt Dt^-1 C x
        virtual void applyt(const BVector& x, BVector& Ax) const;
        /// r = r - Bt Dt^-1 Rw
        virtual void applyt(BVector& r) const;

        // adjoint right hand side of well equations
        // this may at a later point depend on the adjont vectors for prevois step of
        // reservoir and well equations
        void rhsAdjointWell();//(const BVectorWell& lamda_w);

        // add the contributions of the well to the righ has side of the reservoir
        // adjoint equations
        void rhsAdjointRes(BVector& adjRes) const;

        // compute objective derivative contributions used for forming the right hand side
        // and calculating the objective
        void computeObj(Simulator& ebosSimulator,
                        const double dt);


        // update derivative contribution from this well
        void objectDerivative(const BVector& lam_r ,const BVectorWell& lam_w);

        // get the results NB only valid after compute objective Derivative
        void addAdjointResult(AdjointResults& adjres) const;
        // print object function
        void printObjective(std::ostream& os) const;
        // recover adjoint variables for wells and update well_state
        virtual void recoverWellAdjointAndUpdateAdjointState(const BVector& x,
                                                             WellState& well_state);


        /// using the solution x to recover the solution xw for wells and applying
        /// xw to update Well State
        virtual void recoverWellSolutionAndUpdateWellState(const BVector& x,
                                                           WellState& well_state) const;



        /// computing the well potentials for group control
        virtual void computeWellPotentials(const Simulator& ebosSimulator,
                                           const WellState& well_state,
                                           std::vector<double>& well_potentials) /* const */;

        virtual void updatePrimaryVariables(const WellState& well_state) const;

        virtual void solveEqAndUpdateWellState(WellState& well_state);

        virtual void calculateExplicitQuantities(const Simulator& ebosSimulator,
                                                 const WellState& well_state); // should be const?
        void printMatrixes() const{
            std::cout << "duneB " << std::endl;
            Dune::writeMatrixMarket(duneB_, std::cout);
            std::cout << std::endl;
            std::cout << "duneC " << std::endl;
            Dune::writeMatrixMarket(duneC_, std::cout);
            std::cout << std::endl;
            std::cout << "duneD " << std::endl;
            // diagonal matrix for the well
            Dune::writeMatrixMarket(duneD_, std::cout);
            std::cout << "invDuneD " << std::endl;
            // diagonal matrix for the well
            Dune::writeMatrixMarket(invDuneD_, std::cout);
            //std::cout << std::endl;
            // for adjoint
            std::cout << "duneCA " << std::endl;
            Dune::writeMatrixMarket(duneCA_, std::cout);
            std::cout << std::endl;
            //OffDiagMatWellAdjoint duneCA_;
            std::cout << "duneDA " << std::endl;
            Dune::writeMatrixMarket(duneDA_, std::cout);
            std::cout << std::endl;
            std::cout << "adjWell_ " << std::endl;
            Dune::writeMatrixMarket(adjWell_, std::cout);
            std::cout << "objder_adjres_ " << std::endl;
            Dune::writeMatrixMarket(objder_adjres_, std::cout);
            std::cout << "objder_adjwell_ " << std::endl;
            Dune::writeMatrixMarket(objder_adjwell_, std::cout);
            std::cout << "objder_adjctrl_ " << std::endl;
            Dune::writeMatrixMarket(objder_adjctrl_, std::cout);
            std::cout << "adjont_variables " << std::endl;
            Dune::writeMatrixMarket(adjoint_variables_, std::cout);
        }
        BVectorWell getResWell() const {return resWell_;}
        virtual void  addWellContributions(Mat& mat) const;

        /// \brief Wether the Jacobian will also have well contributions in it.
        virtual bool jacobianContainsWellContributions() const
        {
            return param_.matrix_add_well_contributions_;
        }
    protected:

        // protected functions from the Base class
        using Base::getAllowCrossFlow;
        using Base::phaseUsage;
        using Base::flowPhaseToEbosCompIdx;
        using Base::ebosCompIdxToFlowCompIdx;
        using Base::wsolvent;
        using Base::wpolymer;
        using Base::wellHasTHPConstraints;
        using Base::mostStrictBhpFromBhpLimits;
        using Base::scalingFactor;

        // protected member variables from the Base class
        using Base::vfp_properties_;
        using Base::gravity_;
        using Base::param_;
        using Base::well_efficiency_factor_;
        using Base::first_perf_;
        using Base::ref_depth_;
        using Base::perf_depth_;
        using Base::well_cells_;
        using Base::number_of_perforations_;
        using Base::number_of_phases_;
        using Base::saturation_table_number_;
        using Base::comp_frac_;
        using Base::well_index_;
        using Base::index_of_well_;
        using Base::well_controls_;
        using Base::well_type_;
        using Base::num_components_;

        using Base::perf_rep_radius_;
        using Base::perf_length_;
        using Base::bore_diameters_;

        // densities of the fluid in each perforation
        std::vector<double> perf_densities_;
        // pressure drop between different perforations
        std::vector<double> perf_pressure_diffs_;

        // residuals of the well equations
        BVectorWell resWell_;

        // adjoint rhs of the well equations
        BVectorWell adjWell_;

        // two off-diagonal matrices
        OffDiagMatWell duneB_;
        OffDiagMatWell duneC_;

        // diagonal matrix for the well
        DiagMatWell duneD_;// not striktly neeed
        DiagMatWell invDuneD_;

        // for adjoint
        OffDiagMatWellCtrl duneCA_;
        //OffDiagMatWellAdjoint duneCA_;
        DiagMatWellCtrl duneDA_;

        // quatities forobjective function
        Scalar objval_;
        mutable BVectorWellCtrl objder_;
        // well, cells, res primary var
        mutable BVector  objder_adjres_;
        // ... well, well_primary variables
        mutable BVectorWell  objder_adjwell_;
        // ... well, control variables
        mutable BVectorWellCtrl objder_adjctrl_;



        // several vector used in the matrix calculation
        mutable BVectorWell Bx_;
        mutable BVectorWell invDrw_;

        // several vector used in the matrix calculation of transpose solve
        mutable BVectorWell Ctx_;
        mutable BVectorWell invDtadj_;

        // the values for the primary varibles
        // based on different solutioin strategies, the wells can have different primary variables
        mutable std::vector<double> primary_variables_;

        // adjoint variables for well
        mutable BVectorWell adjoint_variables_;
        //mutable std::vector<double> adjoint_variables_;

        // the Evaluation for the well primary variables, which contain derivativles and are used in AD calculation
        mutable std::vector<EvalWell> primary_variables_evaluation_;

        // the saturations in the well bore under surface conditions at the beginning of the time step
        std::vector<double> F0_;

        // TODO: this function should be moved to the base class.
        // while it faces chanllenges for MSWell later, since the calculation of bhp
        // based on THP is never implemented for MSWell yet.
        EvalWell getBhp() const;

        // TODO: it is also possible to be moved to the base class.
        EvalWell getQs(const int comp_idx) const;

        EvalWell wellVolumeFractionScaled(const int phase) const;

        EvalWell wellVolumeFraction(const unsigned compIdx) const;

        EvalWell wellSurfaceVolumeFraction(const int phase) const;

        EvalWell extendEval(const Eval& in) const;

        bool crossFlowAllowed(const Simulator& ebosSimulator) const;

        // xw = inv(D)*(rw - C*x)
        void recoverSolutionWell(const BVector& x, BVectorWell& xw) const;

        // recover well adoint variables
        void recoverAdjointWell(const BVector& x, BVectorWell& xw) const;

        // updating the well_state based on well solution dwells
        void updateWellState(const BVectorWell& dwells,
                             WellState& well_state) const;

        // updating the well_state based on well solution dwells
        void updateAdjointState(const BVectorWell& dwells, WellState& well_state) const;

        // calculate the properties for the well connections
        // to calulate the pressure difference between well connections.
        void computePropertiesForWellConnectionPressures(const Simulator& ebosSimulator,
                                                         const WellState& well_state,
                                                         std::vector<double>& b_perf,
                                                         std::vector<double>& rsmax_perf,
                                                         std::vector<double>& rvmax_perf,
                                                         std::vector<double>& surf_dens_perf) const;

        // TODO: not total sure whether it is a good idea to put this function here
        // the major reason to put here is to avoid the usage of Wells struct
        void computeConnectionDensities(const std::vector<double>& perfComponentRates,
                                        const std::vector<double>& b_perf,
                                        const std::vector<double>& rsmax_perf,
                                        const std::vector<double>& rvmax_perf,
                                        const std::vector<double>& surf_dens_perf);

        void computeConnectionPressureDelta();

        void computeWellConnectionDensitesPressures(const WellState& well_state,
                                                    const std::vector<double>& b_perf,
                                                    const std::vector<double>& rsmax_perf,
                                                    const std::vector<double>& rvmax_perf,
                                                    const std::vector<double>& surf_dens_perf);

        // computing the accumulation term for later use in well mass equations
        void computeAccumWell();

        void computeWellConnectionPressures(const Simulator& ebosSimulator,
                                                    const WellState& well_state);

        // TODO: to check whether all the paramters are required
        void computePerfRate(const IntensiveQuantities& intQuants,
                             const std::vector<EvalWell>& mob_perfcells_dense,
                             const double Tw, const EvalWell& bhp, const double& cdp,
                             const bool& allow_cf, std::vector<EvalWell>& cq_s,
                             double& perf_dis_gas_rate, double& perf_vap_oil_rate) const;

        // TODO: maybe we should provide a light version of computePerfRate, which does not include the
        // calculation of the derivatives
        void computeWellRatesWithBhp(const Simulator& ebosSimulator,
                                     const EvalWell& bhp,
                                     std::vector<double>& well_flux) const;

        std::vector<double> computeWellPotentialWithTHP(const Simulator& ebosSimulator,
                                                        const double initial_bhp, // bhp from BHP constraints
                                                        const std::vector<double>& initial_potential) const;

        template <class ValueType>
        ValueType calculateBhpFromThp(const std::vector<ValueType>& rates, const int control_index) const;

        double calculateThpFromBhp(const std::vector<double>& rates, const int control_index, const double bhp) const;

        // get the mobility for specific perforation
        void getMobility(const Simulator& ebosSimulator,
                         const int perf,
                         std::vector<EvalWell>& mob) const;

        void updateWaterMobilityWithPolymer(const Simulator& ebos_simulator,
                                            const int perf,
                                            std::vector<EvalWell>& mob_water) const;
    };

}

#include "StandardWell_impl.hpp"

#endif // OPM_STANDARDWELL_HEADER_INCLUDED
