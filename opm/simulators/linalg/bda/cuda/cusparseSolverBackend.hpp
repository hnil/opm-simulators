/*
  Copyright 2019 Equinor ASA

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

#ifndef OPM_CUSPARSESOLVER_BACKEND_HEADER_INCLUDED
#define OPM_CUSPARSESOLVER_BACKEND_HEADER_INCLUDED


#include "cublas_v2.h"
#include "cusparse_v2.h"

#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>

namespace Opm
{
namespace Accelerator
{

/// This class implements a cusparse-based ilu0-bicgstab solver on GPU
template <unsigned int block_size>
class cusparseSolverBackend : public BdaSolver<block_size> {

    typedef BdaSolver<block_size> Base;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;
    using Base::deviceID;
    using Base::maxit;
    using Base::tolerance;
    using Base::initialized;

private:

    cublasHandle_t cublasHandle;
    cusparseHandle_t cusparseHandle;
    cudaStream_t stream;
    cusparseMatDescr_t descr_B, descr_M, descr_L, descr_U;
    bsrilu02Info_t info_M;
    bsrsv2Info_t info_L, info_U;
    // b: bsr matrix, m: preconditioner
    double *d_bVals, *d_mVals;
    int *d_bCols, *d_mCols;
    int *d_bRows, *d_mRows;
    double *d_x, *d_b, *d_r, *d_rw, *d_p;     // vectors, used during linear solve
    double *d_pw, *d_s, *d_t, *d_v;
    void *d_buffer;
    double *vals_contiguous;                  // only used if COPY_ROW_BY_ROW is true in cusparseSolverBackend.cpp

    bool analysis_done = false;


    /// Solve linear system using ilu0-bicgstab
    /// \param[in] wellContribs   contains all WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void gpu_pbicgstab(WellContributions& wellContribs, BdaResult& res);

    /// Initialize GPU and allocate memory
    /// \param[in] N                number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] nnz              number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] dim              size of block
    void initialize(int N, int nnz, int dim);

    /// Clean memory
    void finalize();

    /// Copy linear system to GPU
    /// \param[in] vals        array of nonzeroes, each block is stored row-wise, contains nnz values
    /// \param[in] rows        array of rowPointers, contains N/dim+1 values
    /// \param[in] cols        array of columnIndices, contains nnzb values
    /// \param[in] b           input vector, contains N values
    void copy_system_to_gpu(double *vals, int *rows, int *cols, double *b);

    // Update linear system on GPU, don't copy rowpointers and colindices, they stay the same
    /// \param[in] vals        array of nonzeroes, each block is stored row-wise, contains nnz values
    /// \param[in] rows        array of rowPointers, contains N/dim+1 values, only used if COPY_ROW_BY_ROW is true
    /// \param[in] b           input vector, contains N values
    void update_system_on_gpu(double *vals, int *rows, double *b);

    /// Reset preconditioner on GPU, ilu0-decomposition is done inplace by cusparse
    void reset_prec_on_gpu();

    /// Analyse sparsity pattern to extract parallelism
    /// \return true iff analysis was successful
    bool analyse_matrix();

    /// Perform ilu0-decomposition
    /// \return true iff decomposition was successful
    bool create_preconditioner();

    /// Solve linear system
    /// \param[in] wellContribs   contains all WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void solve_system(WellContributions& wellContribs, BdaResult &res);

public:


    /// Construct a cusparseSolver
    /// \param[in] linear_solver_verbosity    verbosity of cusparseSolver
    /// \param[in] maxit                      maximum number of iterations for cusparseSolver
    /// \param[in] tolerance                  required relative tolerance for cusparseSolver
    /// \param[in] deviceID                   the device to be used
    cusparseSolverBackend(int linear_solver_verbosity, int maxit, double tolerance, unsigned int deviceID);

    /// Destroy a cusparseSolver, and free memory
    ~cusparseSolverBackend();

    /// Solve linear system, A*x = b, matrix A must be in blocked-CSR format
    /// \param[in] N              number of rows, divide by dim to get number of blockrows
    /// \param[in] nnz            number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] dim            size of block
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] rows           array of rowPointers, contains N/dim+1 values
    /// \param[in] cols           array of columnIndices, contains nnz values
    /// \param[in] b              input vector, contains N values
    /// \param[in] wellContribs   contains all WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    /// \return                   status code
    SolverStatus solve_system(int N, int nnz, int dim, double *vals, int *rows, int *cols, double *b, WellContributions& wellContribs, BdaResult &res) override;

    /// Get resulting vector x after linear solve, also includes post processing if necessary
    /// \param[inout] x        resulting x vector, caller must guarantee that x points to a valid array
    void get_result(double *x) override;

}; // end class cusparseSolverBackend

} // namespace Accelerator
} // namespace Opm

#endif

