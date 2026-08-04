/*
  Copyright 2019 SINTEF Digital, Mathematics and Cybernetics.

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

#ifndef OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED
#define OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED

#include <dune/common/fvector.hh>

#include <opm/grid/utility/ElementChunks.hpp>
#include <opm/simulators/utils/DeferredLoggingErrorHelpers.hpp>
#include <opm/material/common/MathToolbox.hpp>
#include <opm/models/parallel/threadmanager.hpp>
#include <algorithm>
#include <cmath>

#if HAVE_CUDA
#if USE_HIP
#include <opm/simulators/linalg/gpuistl_hip/detail/cpr_amg_operations.hpp>
#else
#include <opm/simulators/linalg/gpuistl/detail/cpr_amg_operations.hpp>
#endif
#endif


namespace Opm
{

namespace Details
{
    template <class DenseMatrix>
    DenseMatrix transposeDenseMatrix(const DenseMatrix& M)
    {
        DenseMatrix tmp;
        for (int i = 0; i < M.rows; ++i)
            for (int j = 0; j < M.cols; ++j)
                tmp[j][i] = M[i][j];

        return tmp;
    }
} // namespace Details

namespace Amg
{
    /*!
     * \brief The quasi-IMPES weight of one matrix row.
     *
     * Built from the row's diagonal block alone, so it asks nothing of the grid.
     */
    template <class VectorBlockType, class Matrix>
    VectorBlockType quasiImpesWeightForRow(const Matrix& A,
                                           const int rowIdx,
                                           const int pressureVarIndex,
                                           const bool transpose)
    {
        using MatrixBlockType = typename Matrix::block_type;

        VectorBlockType rhs(0.0);
        rhs[pressureVarIndex] = 1.0;

        MatrixBlockType diag_block(0.0);
        const auto row_it = A.begin() + rowIdx;
        const auto endj = (*row_it).end();
        for (auto j = (*row_it).begin(); j != endj; ++j) {
            if (row_it.index() == j.index()) {
                diag_block = (*j);
                break;
            }
        }

        VectorBlockType bweights;
        if (transpose) {
            diag_block.solve(bweights, rhs);
        } else {
            MatrixBlockType diag_block_transpose = Details::transposeDenseMatrix(diag_block);
            diag_block_transpose.solve(bweights, rhs);
        }

        const double abs_max =
            *std::ranges::max_element(bweights,
                                      [](double a, double b)
                                      { return std::fabs(a) < std::fabs(b); });
        bweights /= std::fabs(abs_max);

        return bweights;
    }

    /*!
     * \brief Give the auxiliary degrees of freedom a CPR weight.
     *
     * The true-IMPES weights are built by walking the grid, so they leave the auxiliary
     * degrees of freedom -- which have no element -- untouched.  That is not a small
     * inaccuracy: the weight multiplies the whole row on its way into the coarse pressure
     * system, so a zero weight deletes the row and makes that system singular.
     *
     * They get the quasi-IMPES weight instead, which needs only the assembled diagonal
     * block.  Mixing the two is a compromise on scaling, not on correctness; an auxiliary
     * cell's true-IMPES weight needs the storage term evaluated without an element
     * context, which is the same thing the TPFA linearizer does and is the natural
     * follow-up.
     */
    template <class Matrix, class Vector>
    void getAuxiliaryDofWeights(const Matrix& matrix,
                                const int firstAuxiliaryRow,
                                const int pressureVarIndex,
                                const bool transpose,
                                Vector& weights)
    {
        using VectorBlockType = typename Vector::block_type;

        for (int rowIdx = firstAuxiliaryRow; rowIdx < static_cast<int>(matrix.N()); ++rowIdx) {
            weights[rowIdx] = quasiImpesWeightForRow<VectorBlockType>(matrix, rowIdx,
                                                                     pressureVarIndex, transpose);
        }
    }

    template <class Matrix, class Vector>
    void getQuasiImpesWeights(const Matrix& matrix,
                              const int pressureVarIndex,
                              const bool transpose,
                              Vector& weights,
                              [[maybe_unused]] bool enable_thread_parallel)
    {
        using VectorBlockType = typename Vector::block_type;
        using MatrixBlockType = typename Matrix::block_type;
        const Matrix& A = matrix;

        VectorBlockType rhs(0.0);
        rhs[pressureVarIndex] = 1.0;

        // Declare variables outside the loop to avoid repetitive allocation
        MatrixBlockType diag_block;
        VectorBlockType bweights;
        MatrixBlockType diag_block_transpose;

        // Use OpenMP to parallelize over matrix rows (runtime controlled via if clause)
#ifdef _OPENMP
#pragma omp parallel for private(diag_block, bweights, diag_block_transpose) if(enable_thread_parallel)
#endif
        for (int row_idx = 0; row_idx < static_cast<int>(A.N()); ++row_idx) {
            diag_block = MatrixBlockType(0.0);
            // Find diagonal block for this row
            const auto row_it = A.begin() + row_idx;
            const auto endj = (*row_it).end();
            for (auto j = (*row_it).begin(); j != endj; ++j) {
                if (row_it.index() == j.index()) {
                    diag_block = (*j);
                    break;
                }
            }
            if (transpose) {
                diag_block.solve(bweights, rhs);
            } else {
                diag_block_transpose = Details::transposeDenseMatrix(diag_block);
                diag_block_transpose.solve(bweights, rhs);
            }

            const double abs_max =
                *std::ranges::max_element(bweights,
                                          [](double a, double b)
                                          { return std::fabs(a) < std::fabs(b); });
            bweights /= std::fabs(abs_max);
            weights[row_idx] = bweights;
        }
    }

    template <class Matrix, class Vector>
    Vector getQuasiImpesWeights(const Matrix& matrix,
                                const int pressureVarIndex,
                                const bool transpose,
                                bool enable_thread_parallel)
    {
        Vector weights(matrix.N());
        getQuasiImpesWeights(matrix, pressureVarIndex, transpose, weights, enable_thread_parallel);
        return weights;
    }

#if HAVE_CUDA
    template <typename T>
    std::vector<int> precomputeDiagonalIndices(const gpuistl::GpuSparseMatrixWrapper<T>& matrix) {
        std::vector<int> diagonalIndices(matrix.N(), -1);
        const auto rowIndices = matrix.getRowIndices().asStdVector();
        const auto colIndices = matrix.getColumnIndices().asStdVector();

        for (auto row = 0; row < Opm::gpuistl::detail::to_int(matrix.N()); ++row) {
            for (auto i = rowIndices[row]; i < rowIndices[row+1]; ++i) {
                if (colIndices[i] == row) {
                    diagonalIndices[row] = i;
                    break;
                }
            }
        }
        return diagonalIndices;
    }

    // GPU version that delegates to the GPU implementation
    template <typename T, bool transpose>
    void getQuasiImpesWeights(const gpuistl::GpuSparseMatrixWrapper<T>& matrix,
                             const int pressureVarIndex,
                             gpuistl::GpuVector<T>& weights,
                             const gpuistl::GpuVector<int>& diagonalIndices)
    {
        gpuistl::detail::getQuasiImpesWeights<T, transpose>(matrix, pressureVarIndex, weights, diagonalIndices);
    }

    template <typename T, bool transpose>
    gpuistl::GpuVector<T> getQuasiImpesWeights(const gpuistl::GpuSparseMatrixWrapper<T>& matrix,
                                              const int pressureVarIndex,
                                              const gpuistl::GpuVector<int>& diagonalIndices)
    {
        gpuistl::GpuVector<T> weights(matrix.N() * matrix.blockSize());
        getQuasiImpesWeights<T, transpose>(matrix, pressureVarIndex, weights, diagonalIndices);
        return weights;
    }
#endif

    template<class Vector, class ElementContext, class Model, class ElementChunksType>
    void getTrueImpesWeights(int pressureVarIndex, Vector& weights,
                             const ElementContext& elemCtx,
                             const Model& model,
                             const ElementChunksType& element_chunks,
                             [[maybe_unused]] bool enable_thread_parallel)
    {
        using VectorBlockType = typename Vector::block_type;
        using Matrix = typename std::decay_t<decltype(model.linearizer().jacobian())>;
        using MatrixBlockType = typename Matrix::MatrixBlock;
        constexpr int numEq = VectorBlockType::size();
        using Evaluation = typename std::decay_t<decltype(model.localLinearizer(ThreadManager::threadId()).localResidual().residual(0))>
            ::block_type;

        using LocalResidual = std::decay_t<
            decltype(model.localLinearizer(ThreadManager::threadId()).localResidual())>;

        VectorBlockType rhs(0.0);
        rhs[pressureVarIndex] = 1.0;

        // Declare variables outside the loop to avoid repetitive allocation
        MatrixBlockType block;
        VectorBlockType bweights;
        MatrixBlockType block_transpose;
        Dune::FieldVector<Evaluation, numEq> storage;

        // Turn one degree of freedom's accumulation term into its weight.  Only the
        // storage derivatives and the volume the storage is scaled by enter it.
        const auto weightFromStorage = [pressureVarIndex, &rhs, &block_transpose]
            (const auto& stor, const auto storage_scale, auto& bw)
        {
            const double pressure_scale = 50e5;

            // Build the transposed matrix directly to avoid separate transpose step
            for (int ii = 0; ii < numEq; ++ii) {
                for (int jj = 0; jj < numEq; ++jj) {
                    block_transpose[jj][ii] = stor[ii].derivative(jj)/storage_scale;
                    if (jj == pressureVarIndex) {
                        block_transpose[jj][ii] *= pressure_scale;
                    }
                }
            }
            block_transpose.solve(bw, rhs);

            const double abs_max =
                *std::ranges::max_element(bw,
                                          [](double a, double b)
                                          { return std::fabs(a) < std::fabs(b); });
            // probably a scaling which could give approximately total compressibility would be better
            bw /=  std::fabs(abs_max); // given normal densities this scales weights to about 1.
        };

        // Walk the degrees of freedom when that is possible: the accumulation term of a
        // degree of freedom is a function of its own intensive quantities and its own
        // volume, so no element is needed to form it -- and a degree of freedom appended
        // after the grid has no element to be reached through.
        //
        // It takes two things the element walk gets from the ElementContext: intensive
        // quantities by index, which is what the cache is, and a storage term computed
        // from them alone, which the TPFA local residual has and the element-by-element
        // one does not.  Those are also precisely the conditions under which a degree of
        // freedom outside the grid can exist at all -- the TPFA linearizer is the one
        // that assembles their equations -- so when either is missing there is nothing
        // beyond the grid to miss, and the element walk below is complete.
        if constexpr (requires (Dune::FieldVector<Evaluation, numEq>& s, const Model& m)
                      {
                          LocalResidual::template computeStorage<Evaluation>
                              (s, m.intensiveQuantities(0u, 0u));
                          m.intensiveQuantityCacheEnabled();
                          m.dofTotalVolume(0);
                      })
        {
            if (model.intensiveQuantityCacheEnabled()) {
                const auto numDof = static_cast<int>(model.numTotalDof());
                const auto dt = elemCtx.simulator().timeStepSize();

                OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for private(bweights, block_transpose, storage) if(enable_thread_parallel)
#endif
                const auto& matrix = model.linearizer().jacobian().istlMatrix();

                for (int dofIdx = 0; dofIdx < numDof; ++dofIdx) {
                    const auto& intQuants =
                        model.intensiveQuantities(static_cast<unsigned>(dofIdx), /*timeIdx=*/0);

                    const auto scvVolume =
                        model.dofTotalVolume(dofIdx) * intQuants.extrusionFactor();

                    // The true-IMPES weight is the accumulation term scaled by the volume
                    // it accumulates in, so a degree of freedom occupying no volume has
                    // no such weight: dividing by it gives NaN, and the whole linear
                    // solve with it.  A reserved but not yet occupied auxiliary degree of
                    // freedom is exactly that -- it carries an identity row until
                    // something claims it -- and the weight its own row implies is what is
                    // wanted for it.  Grid cells always have a volume, so this never fires
                    // where there are no auxiliary degrees of freedom.
                    if (!(scvVolume > 0.0)) {
                        weights[dofIdx] = quasiImpesWeightForRow<VectorBlockType>
                            (matrix, dofIdx, pressureVarIndex, /*transpose=*/false);
                        continue;
                    }

                    LocalResidual::template computeStorage<Evaluation>(storage, intQuants);

                    weightFromStorage(storage, scvVolume / dt, bweights);
                    weights[dofIdx] = bweights;
                }
                OPM_END_PARALLEL_TRY_CATCH("getTrueImpesWeights() failed: ",
                                           elemCtx.simulator().vanguard().grid().comm());

                return;
            }
        }

        OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for private(block, bweights, block_transpose, storage) if(enable_thread_parallel)
#endif
        for (const auto& chunk : element_chunks) {
            const std::size_t thread_id = ThreadManager::threadId();
            ElementContext localElemCtx(elemCtx.simulator());

            for (const auto& elem : chunk) {
                localElemCtx.updatePrimaryStencil(elem);
                localElemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                model.localLinearizer(thread_id).localResidual().computeStorage(storage, localElemCtx, /*spaceIdx=*/0, /*timeIdx=*/0);

                auto extrusionFactor = localElemCtx.intensiveQuantities(0, /*timeIdx=*/0).extrusionFactor();
                auto scvVolume = localElemCtx.stencil(/*timeIdx=*/0).subControlVolume(0).volume() * extrusionFactor;

                weightFromStorage(storage,
                                  scvVolume / localElemCtx.simulator().timeStepSize(),
                                  bweights);

                const auto index = localElemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                weights[index] = bweights;
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("getTrueImpesWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }

    template <class Vector, class ElementContext, class Model, class ElementChunksType>
    void getTrueImpesWeightsAnalytic([[maybe_unused]] int pressureVarIndex,
                                     Vector& weights,
                                     const ElementContext& elemCtx,
                                     const Model& model,
                                     const ElementChunksType& element_chunks,
                                     [[maybe_unused]] bool enable_thread_parallel)
    {
        // The sequential residual is a linear combination of the
        // mass balance residuals, with coefficients equal to (for
        // water, oil, gas):
        //    1/bw,
        //    (1/bo - rs/bg)/(1-rs*rv)
        //    (1/bg - rv/bo)/(1-rs*rv)
        // These coefficients must be applied for both the residual and
        // Jacobian.
        using FluidSystem = typename Model::FluidSystem;
        using LhsEval = double;

        using PrimaryVariables = typename Model::PrimaryVariables;
        using VectorBlockType = typename Vector::block_type;
        using Evaluation =
            typename std::decay_t<decltype(model.localLinearizer(ThreadManager::threadId()).localResidual().residual(0))>::block_type;
        using Toolbox = MathToolbox<Evaluation>;

        const auto& solution = model.solution(/*timeIdx*/ 0);
        VectorBlockType bweights;

        // The weight of one degree of freedom.  It is a function of that degree of
        // freedom's fluid state and primary variables and of nothing else -- no
        // neighbours, no geometry -- so it can be had by index as readily as by element.
        const auto weightOf = [&solution]
            (const auto& intQuants, const auto index, auto& bw)
        {
            const auto& fs = intQuants.fluidState();

            if (FluidSystem::phaseIsActive(FluidSystem::waterPhaseIdx)) {
                const unsigned activeCompIdx = FluidSystem::canonicalToActiveCompIdx(
                    FluidSystem::solventComponentIndex(FluidSystem::waterPhaseIdx));
                bw[activeCompIdx]
                    = Toolbox::template decay<LhsEval>(1 / fs.invB(FluidSystem::waterPhaseIdx));
            }

            double denominator = 1.0;
            double rs = Toolbox::template decay<double>(fs.Rs());
            double rv = Toolbox::template decay<double>(fs.Rv());
            const auto& priVars = solution[index];
            if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rv) {
                rs = 0.0;
            }
            if (priVars.primaryVarsMeaningGas() == PrimaryVariables::GasMeaning::Rs) {
                rv = 0.0;
            }
            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)
                && FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                denominator = Toolbox::template decay<LhsEval>(1 - rs * rv);
            }

            if (FluidSystem::phaseIsActive(FluidSystem::oilPhaseIdx)) {
                const unsigned activeCompIdx = FluidSystem::canonicalToActiveCompIdx(
                    FluidSystem::solventComponentIndex(FluidSystem::oilPhaseIdx));
                bw[activeCompIdx] = Toolbox::template decay<LhsEval>(
                    (1 / fs.invB(FluidSystem::oilPhaseIdx) - rs / fs.invB(FluidSystem::gasPhaseIdx))
                    / denominator);
            }
            if (FluidSystem::phaseIsActive(FluidSystem::gasPhaseIdx)) {
                const unsigned activeCompIdx = FluidSystem::canonicalToActiveCompIdx(
                    FluidSystem::solventComponentIndex(FluidSystem::gasPhaseIdx));
                bw[activeCompIdx] = Toolbox::template decay<LhsEval>(
                    (1 / fs.invB(FluidSystem::gasPhaseIdx) - rv / fs.invB(FluidSystem::oilPhaseIdx))
                    / denominator);
            }
        };

        // Walking the degrees of freedom is what this wants to do: it needs no element,
        // and the auxiliary degrees of freedom appended after the grid have none to be
        // reached through.  It is possible exactly when their intensive quantities can be
        // had by index, which is also the condition under which they exist at all -- so
        // when the cache is off there is nothing beyond the grid to miss, and the element
        // walk below is complete.
        if (model.intensiveQuantityCacheEnabled()) {
            const auto numDof = static_cast<int>(model.numTotalDof());
            const auto& matrix = model.linearizer().jacobian().istlMatrix();

#ifdef _OPENMP
#pragma omp parallel for private(bweights) if(enable_thread_parallel)
#endif
            for (int dofIdx = 0; dofIdx < numDof; ++dofIdx) {
                // A degree of freedom occupying no volume holds no fluid, so there is no
                // fluid state to take a weight from -- and the cached one, which nothing
                // has had reason to keep current, is not required to mean anything.  A
                // reserved but unoccupied auxiliary degree of freedom is exactly that; it
                // carries an identity row until something claims it, and the weight that
                // row implies is both defined and the one that belongs to it.  A grid
                // cell always has a volume, so this cannot fire where there are no
                // auxiliary degrees of freedom.
                if (!(model.dofTotalVolume(dofIdx) > 0.0)) {
                    weights[dofIdx] = quasiImpesWeightForRow<VectorBlockType>
                        (matrix, dofIdx, pressureVarIndex, /*transpose=*/false);
                    continue;
                }

                weightOf(model.intensiveQuantities(static_cast<unsigned>(dofIdx), /*timeIdx=*/0),
                         dofIdx, bweights);
                weights[dofIdx] = bweights;
            }

            return;
        }

        // Use OpenMP to parallelize over element chunks (runtime controlled via if clause)
        OPM_BEGIN_PARALLEL_TRY_CATCH();
#ifdef _OPENMP
#pragma omp parallel for private(bweights) if(enable_thread_parallel)
#endif
        for (const auto& chunk : element_chunks) {

            // Each thread gets a unique copy of elemCtx
            ElementContext localElemCtx(elemCtx.simulator());

            for (const auto& elem : chunk) {
                localElemCtx.updatePrimaryStencil(elem);
                localElemCtx.updatePrimaryIntensiveQuantities(/*timeIdx=*/0);

                const auto index = localElemCtx.globalSpaceIndex(/*spaceIdx=*/0, /*timeIdx=*/0);
                weightOf(localElemCtx.intensiveQuantities(/*spaceIdx=*/0, /*timeIdx=*/0),
                         index, bweights);
                weights[index] = bweights;
            }
        }
        OPM_END_PARALLEL_TRY_CATCH("getTrueImpesAnalyticWeights() failed: ", elemCtx.simulator().vanguard().grid().comm());
    }
} // namespace Amg

} // namespace Opm

#endif // OPM_GET_QUASI_IMPES_WEIGHTS_HEADER_INCLUDED
