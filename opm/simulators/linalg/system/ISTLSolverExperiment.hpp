#pragma once

#include <opm/simulators/linalg/ISTLSolver.hpp>
#include "WellMatrixMerger.hpp"
namespace Opm
{
        namespace SystemSolver
        {
                const int numResDofs = 3;
                const int numWellDofs = 4;

                // Define matrix and vector types
                // using RRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numResDofs>>;
                using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
                // using RWtype = Dune::FieldMatrix<double, numResDofs, numWellDofs>;
                using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
                using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
                using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
                using RVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
                using WVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
                // Define system matrix and vector types
                using SystemMatrix = Dune::MultiTypeBlockMatrix<
                    Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                    Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>>;
                using SystemVector = Dune::MultiTypeBlockVector<RVector, WVector>;
                using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;

                void getSystemSolver(std::unique_ptr<Dune::InverseOperator<SystemVector, SystemVector>>& systemsolver,
                        std::unique_ptr<Dune::InverseOperator<SystemVector, SystemVector>>& syspreconditioner,
                        std::unique_ptr< Dune::AssembledLinearOperator<SystemMatrix, SystemVector, SystemVector>>& system_operator,
                        std::shared_ptr< Dune::ScalarProduct<SystemVector>>& scalarproduct,
                        const SystemComm& systemcomm,
                        const SystemMatrix &S,
                        const std::function<RVector()> &weightCalculator,
                        int pressureIndex, const Opm::PropertyTree &prm);
                void getSystemSolver(Dune::InverseOperator<SystemVector, SystemVector> &systemsolver,
                        const SystemMatrix &S,
                Dune::InverseOperatorResult solveSystem(const SystemMatrix &S, SystemVector &x, const SystemVector &b,  
                        const std::function<RVector()> &weightCalculator,
                        int pressureIndex, const Opm::PropertyTree &prm, const Comm &comm);
                Dune::InverseOperatorResult solveSystem(const SystemMatrix &S, SystemVector &x, const SystemVector &b,  
                        const std::function<RVector()> &weightCalculator,
                        int pressureIndex, const Opm::PropertyTree &prm);        
               
        } // namespace SystemSolver
        template <class TypeTag>
        class ISTLSolverExperiment : public ISTLSolver<TypeTag>
        {

        protected:
                using GridView = GetPropType<TypeTag, Properties::GridView>;
                using Scalar = GetPropType<TypeTag, Properties::Scalar>;
                using SparseMatrixAdapter = GetPropType<TypeTag, Properties::SparseMatrixAdapter>;
                using Vector = GetPropType<TypeTag, Properties::GlobalEqVector>;
                using Indices = GetPropType<TypeTag, Properties::Indices>;
                using WellModel = GetPropType<TypeTag, Properties::WellModel>;
                using Simulator = GetPropType<TypeTag, Properties::Simulator>;
                using Matrix = typename SparseMatrixAdapter::IstlMatrix;
                using ThreadManager = GetPropType<TypeTag, Properties::ThreadManager>;
                using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
                using AbstractSolverType = Dune::InverseOperator<Vector, Vector>;
                using AbstractOperatorType = Dune::AssembledLinearOperator<Matrix, Vector, Vector>;
                using AbstractPreconditionerType = Dune::PreconditionerWithUpdate<Vector, Vector>;
                using WellModelOperator = WellModelAsLinearOperator<WellModel, Vector, Vector>;
                using ElementMapper = GetPropType<TypeTag, Properties::ElementMapper>;
                using ElementChunksType = ElementChunks<GridView, Dune::Partitions::All>;

                constexpr static std::size_t pressureIndex = GetPropType<TypeTag, Properties::Indices>::pressureSwitchIdx;

                using SystemComm = Dune::MultiCommunicator<const Dune::OwnerOverlapCopyCommunication<int, int>&,const WellComm&>;
                const int numResDofs = 3;
                const int numWellDofs = 4;
                // using Indices = GetPropType<TypeTag, Properties::Indices>;
                // static constexpr int numResDofs = Indices::numEq;
                // static constexpr int numWellDofs = numWellDofs + 1;
                //  Define matrix and vector types
                // using RRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numResDofs>>;
                using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
                // using RWtype = Dune::FieldMatrix<double, numResDofs, numWellDofs>;
                using RWMatrix
                    = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
                using WRMatrix
                    = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
                using WWMatrix
                    = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
                using RVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
                using WVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
                // Define system matrix and vector types
                using SystemMatrix
                    = Dune::MultiTypeBlockMatrix<Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                                                 Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>>;
                using SystemVector = Dune::MultiTypeBlockVector<RVector, WVector>;


                enum
                {
                        enablePolymerMolarWeight = getPropValue<TypeTag, Properties::EnablePolymerMW>()
                };
                constexpr static bool isIncompatibleWithCprw = enablePolymerMolarWeight;

#if HAVE_MPI
                using CommunicationType = Dune::OwnerOverlapCopyCommunication<int, int>;
#else
                using CommunicationType = Dune::Communication<int>;
#endif
                using Parent = ISTLSolver<TypeTag>;

        public:
                ISTLSolverExperiment(const Simulator &simulator,
                                     const FlowLinearSolverParameters &parameters,
                                     bool forceSerial = false)
                    : Parent(simulator, parameters, forceSerial)
                {
                }
                explicit ISTLSolverExperiment(const Simulator &simulator)
                    : Parent(simulator)
                {
                }
                void prepare(const Matrix& M, Vector& b) override
                {               
                        OPM_TIMEBLOCK(istlSolverPrepare);
                        try {
                        Parent::initPrepare(M,b);
                        const auto& prm = this->prm_[this->activeSolverNum_];
                        bool solve_system = prm.get("use_system_solver", false);
                        if(!solve_system){
                                Parent::prepareFlexibleSolver();
                        }else{
                             systemsolver_ = nullptr; // reset system solver
                             systemsolver_   
                        }        
                        } OPM_CATCH_AND_RETHROW_AS_CRITICAL_ERROR("This is likely due to a faulty linear solver JSON specification. Check for errors related to missing nodes.");
                        
                }
                void setupSystemSolver()
                {
                        OPM_TIMEBLOCK(ISTLSolverExperiment_solve);
                        // Here we could add experimental features before or after calling the base class solve.
                        
                                
                                // get reservoir matrix
                                const int numResDofs = 3;
                                const int numWellDofs = 4;
                                // using Indices = GetPropType<TypeTag, Properties::Indices>;
                                // static constexpr int numResDofs = Indices::numEq;
                                // static constexpr int numWellDofs = numWellDofs + 1;
                                //  Define matrix and vector types
                                // using RRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numResDofs>>;
                                using RRMatrix = Dune::BCRSMatrix<Opm::MatrixBlock<double, numResDofs, numResDofs>>;
                                // using RWtype = Dune::FieldMatrix<double, numResDofs, numWellDofs>;
                                using RWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numResDofs, numWellDofs>>;
                                using WRMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numResDofs>>;
                                using WWMatrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, numWellDofs, numWellDofs>>;
                                using RVector = Dune::BlockVector<Dune::FieldVector<double, numResDofs>>;
                                using WVector = Dune::BlockVector<Dune::FieldVector<double, numWellDofs>>;
                                // Define system matrix and vector types
                                using SystemMatrix = Dune::MultiTypeBlockMatrix<
                                    Dune::MultiTypeBlockVector<RRMatrix, RWMatrix>,
                                    Dune::MultiTypeBlockVector<WRMatrix, WWMatrix>>;
                                using SystemVector = Dune::MultiTypeBlockVector<RVector, WVector>;

                                // Define helper constants for accessing MultiTypeBlockMatrix elements
                                constexpr auto _0 = Dune::Indices::_0;
                                constexpr auto _1 = Dune::Indices::_1;
                                A_copy_ = *Parent::matrix_;
                                Opm::WellMatrixMerger merger(A_copy.N());
                                WRMatrix B1;
                                RWMatrix C1;
                                WWMatrix D1;
                                WVector w_res;
                                std::vector<int> cells1 = {}; // Cells for well 1
                                std::vector<WRMatrix> b_matrices;
                                std::vector<RWMatrix> c_matrices;
                                std::vector<WWMatrix> d_matrices;
                                std::vector<std::vector<int>> wcells;
                                std::vector<WVector> residual;

                                this->simulator_.problem().wellModel().addBCDMatrix(b_matrices, c_matrices, d_matrices, wcells, residual);
                                for (size_t i = 0; i < b_matrices.size(); ++i)
                                {
                                        merger.addWell(b_matrices[i], c_matrices[i], d_matrices[i], wcells[i], static_cast<int>(i), "Well" + std::to_string(i + 1));
                                }
                                //   merger.addWell(B1, C1, D1, cells1, 0, "Well1");

                                merger.finalize();
                                


                                // Get the merged matrices
                                const auto &mergedB = merger.getMergedB();
                                const auto &mergedC = merger.getMergedC();
                                const auto &mergedD = merger.getMergedD();
                                
                                C_copy_ = mergedC;
                                B_copy_ = mergedB;
                                D_copy_ = mergedD;
                                S_[_0][_0] = A_copy;
                                S_[_0][_1] = C_copy;
                                S_[_1][_0] = B_copy;
                                S_[_1][_1] = D_copy;

                                //not used now
                                for (size_t i = 0; i < residual.size(); ++i)
                                {
                                        for (size_t j = 0; j < residual[i].size(); ++j)
                                        {
                                                w_res[well_dof] += residual[i][j];
                                                well_dof++;
                                        }
                                }
                                //std::cout << "Well residual norm: " << w_res.two_norm2() << std::endl;
                                //w_res = 0.0;// this should be applied to reservoir at this point
                                const auto& prm_system = prm.get_child("system_solver");
                                weightCalculator_ = this->getWeightsCalculator(prm_system.get_child("preconditioner.reservoir_solver"), this->getMatrix(), pressureIndex);
                                #if HAVE_MPI
                                        using Comm = Dune::OwnerOverlapCopyCommunication<int, int>;
                                #endif
                                        
                                const Comm& comm = *(this->comm_.get());        
                                if(this->comm_->communicator().size() > 1)
                                {
                                  //std::cout << "Parallel run with system solver..." << this->comm_->communicator().rank() <<  std::endl;
                                        const auto result = SystemSolver::solveSystem(S, x_s, r_s, weightCalculator, pressureIndex, prm_system, comm);
                                        this->iterations_ = result.iterations;
                                }else
                                {
                                   const auto result = SystemSolver::solveSystem(S, x_s, r_s, weightCalculator, pressureIndex, prm_system);        // Serial run
                                   this->iterations_ = result.iterations;
                                }
                }                
        bool solve(Vector &x) override{
                        const auto& prm = this->prm_[this->activeSolverNum_];
                        bool solve_system = prm.get("use_system_solver", false);
                        if (solve_system)
                        {
                                RVector x_r(A_copy_.N());
                                x_r = 0.0;
                                WVector x_w(D_copy_.N());
                                x_w = 0.0;
                                RVector r_res = *Parent::rhs_;
                                WVector w_res(D_copy.N());
                                w_res = x_w;
                                SystemVector x_s{x_r, x_w};
                                SystemVector r_s{r_res, w_res};
                                systemsolver_->apply(x_s, r_s);      
                                x = x_s[_0];         
                        }else
                        {
                                // Call the base class solve method
                                return Parent::solve(x);
                        }
                        return true;
                }

        private:
                RRMatrix A_copy_;
                RWMatrix C_copy_;
                WRMatrix B_copy_;
                WWMatrix D_copy_;
                SystemMatrix S_;
                std::function<RVector()> weightCalculator_;
                std::unique_ptr<Dune::InverseOperator<SystemVector, SystemVector>> systemsolver_;
                WellComm seqComm_;
                SystemComm systemComm_;//(comm,seqComm);
                //std::shared_ptr<Dune::Preconditioner<SystemVector, SystemVector>> precond_;
                std::shared_ptr<Dune::Preconditioner<SystemVector, SystemVector>> sysprecond_;//, SystemComm, SystemPreconditionerParallel> sysprecond_;//(precond, systemComm);
                std::shared_ptr< Dune::ScalarProduct<SystemVector> > scalarproduct_;
                std::shared_ptr< Dune::AssembledLinearOperator<SystemMatrix, SystemVector, SystemVector> > system_operator_;
        };

}
