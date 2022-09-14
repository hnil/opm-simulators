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

#ifndef OPM_OPENCLSOLVER_BACKEND_HEADER_INCLUDED
#define OPM_OPENCLSOLVER_BACKEND_HEADER_INCLUDED

#include <opm/simulators/linalg/bda/opencl/opencl.hpp>
#include <opm/simulators/linalg/bda/opencl/openclKernels.hpp>
#include <opm/simulators/linalg/bda/BdaResult.hpp>
#include <opm/simulators/linalg/bda/BdaSolver.hpp>
#include <opm/simulators/linalg/bda/ILUReorder.hpp>
#include <opm/simulators/linalg/bda/WellContributions.hpp>
#include <opm/simulators/linalg/bda/opencl/BILU0.hpp>

#include <opm/simulators/linalg/bda/opencl/Preconditioner.hpp>

#include <tuple>

namespace Opm
{
namespace Accelerator
{

template <unsigned int block_size>
class CPR;

/// This class implements a opencl-based ilu0-bicgstab solver on GPU
template <unsigned int block_size>
class openclSolverBackend : public BdaSolver<block_size>
{
    typedef BdaSolver<block_size> Base;

    using Base::N;
    using Base::Nb;
    using Base::nnz;
    using Base::nnzb;
    using Base::verbosity;
    using Base::platformID;
    using Base::deviceID;
    using Base::maxit;
    using Base::tolerance;
    using Base::initialized;

private:
    double *rb = nullptr;                 // reordered b vector, if the matrix is reordered, rb is newly allocated, otherwise it just points to b
    std::vector<double> vals_contiguous;  // only used if COPY_ROW_BY_ROW is true in openclSolverBackend.cpp

    // OpenCL variables must be reusable, they are initialized in initialize()
    cl::Buffer d_Avals, d_Acols, d_Arows;        // (reordered) matrix in BSR format on GPU
    cl::Buffer d_x, d_b, d_rb, d_r, d_rw, d_p;   // vectors, used during linear solve
    cl::Buffer d_pw, d_s, d_t, d_v;              // vectors, used during linear solve
    cl::Buffer d_tmp;                            // used as tmp GPU buffer for dot() and norm()
    cl::Buffer d_toOrder;                        // only used when reordering is used

    std::vector<cl::Device> devices;

    std::unique_ptr<Preconditioner<block_size> > prec;
                                                                  // can perform blocked ILU0 and AMG on pressure component
    bool is_root;                                                 // allow for nested solvers, the root solver is called by BdaBridge
    int *toOrder = nullptr, *fromOrder = nullptr;                 // BILU0 reorders rows of the matrix via these mappings
    bool analysis_done = false;
    std::unique_ptr<BlockedMatrix> mat = nullptr;                 // original matrix
    BlockedMatrix *rmat = nullptr;                                // reordered matrix (or original if no reordering), used for spmv
    ILUReorder opencl_ilu_reorder;                                // reordering strategy
    std::vector<cl::Event> events;
    cl_int err;

    /// Divide A by B, and round up: return (int)ceil(A/B)
    /// \param[in] A    dividend
    /// \param[in] B    divisor
    /// \return         rounded division result
    unsigned int ceilDivision(const unsigned int A, const unsigned int B);

    /// Calculate dot product between in1 and in2, partial sums are stored in out, which are summed on CPU
    /// \param[in] in1           input vector 1
    /// \param[in] in2           input vector 2
    /// \param[out] out          output vector containing partial sums
    /// \return                  dot product
    double dot_w(cl::Buffer in1, cl::Buffer in2, cl::Buffer out);

    /// Calculate the norm of in, partial sums are stored in out, which are summed on the CPU
    /// Equal to Dune::DenseVector::two_norm()
    /// \param[in] in          input vector
    /// \param[out] out        output vector containing partial sums
    /// \return                norm
    double norm_w(cl::Buffer in, cl::Buffer out);

    /// Perform axpy: out += a * in
    /// \param[in] in         input vector
    /// \param[in] a          scalar value to multiply input vector
    /// \param[inout] out     output vector
    void axpy_w(cl::Buffer in, const double a, cl::Buffer out);

    /// Perform scale: vec *= a
    /// \param[inout] vec     vector to scale
    /// \param[in] a          scalar value to multiply vector
    void scale_w(cl::Buffer vec, const double a);

    /// Custom function that combines scale, axpy and add functions in bicgstab
    /// p = (p - omega * v) * beta + r
    /// \param[inout] p      output vector
    /// \param[in] v         input vector
    /// \param[in] r         input vector
    /// \param[in] omega     scalar value
    /// \param[in] beta      scalar value
    void custom_w(cl::Buffer p, cl::Buffer v, cl::Buffer r, const double omega, const double beta);

    /// Sparse matrix-vector multiply, spmv
    /// b = A * x
    /// Matrix A, must be in BCRS format
    /// \param[in] vals     nnzs of matrix A
    /// \param[in] cols     columnindices of matrix A
    /// \param[in] rows     rowpointers of matrix A
    /// \param[in] x        input vector
    /// \param[out] b       output vector
    void spmv_blocked_w(cl::Buffer vals, cl::Buffer cols, cl::Buffer rows, cl::Buffer x, cl::Buffer b);

    /// Solve linear system using ilu0-bicgstab
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void gpu_pbicgstab(WellContributions& wellContribs, BdaResult& res);

    /// Initialize GPU and allocate memory
    /// \param[in] N              number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] nnz            number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] dim            size of block
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] rows           array of rowPointers, contains N/dim+1 values
    /// \param[in] cols           array of columnIndices, contains nnz values
    void initialize(int N, int nnz, int dim, double *vals, int *rows, int *cols);

    /// Clean memory
    void finalize();

    /// Copy linear system to GPU
    void copy_system_to_gpu();

    /// Reorder the linear system so it corresponds with the coloring
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] b              input vectors, contains N values
    /// \param[out] wellContribs  WellContributions, to set reordering
    void update_system(double *vals, double *b, WellContributions &wellContribs);

    /// Update linear system on GPU, don't copy rowpointers and colindices, they stay the same
    void update_system_on_gpu();

    /// Analyze sparsity pattern to extract parallelism
    /// \return true iff analysis was successful
    bool analyze_matrix();

    /// Perform ilu0-decomposition
    /// \return true iff decomposition was successful
    bool create_preconditioner();

    /// Solve linear system
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    void solve_system(WellContributions &wellContribs, BdaResult &res);

public:
    std::shared_ptr<cl::Context> context;
    std::shared_ptr<cl::CommandQueue> queue;

    /// Construct a openclSolver
    /// \param[in] linear_solver_verbosity    verbosity of openclSolver
    /// \param[in] maxit                      maximum number of iterations for openclSolver
    /// \param[in] tolerance                  required relative tolerance for openclSolver
    /// \param[in] platformID                 the OpenCL platform to be used
    /// \param[in] deviceID                   the device to be used
    /// \param[in] opencl_ilu_reorder         select either level_scheduling or graph_coloring, see BILU0.hpp for explanation
    /// \param[in] linsolver                  indicating the preconditioner, equal to the --linsolver cmdline argument
    ///                                       only ilu0 and cpr_quasiimpes are supported
    openclSolverBackend(int linear_solver_verbosity, int maxit, double tolerance, unsigned int platformID, unsigned int deviceID,
        ILUReorder opencl_ilu_reorder, std::string linsolver);

    /// For the CPR coarse solver
    openclSolverBackend(int linear_solver_verbosity, int maxit, double tolerance, ILUReorder opencl_ilu_reorder);

    /// Destroy a openclSolver, and free memory
    ~openclSolverBackend();

    /// Solve linear system, A*x = b, matrix A must be in blocked-CSR format
    /// \param[in] N              number of rows, divide by dim to get number of blockrows
    /// \param[in] nnz            number of nonzeroes, divide by dim*dim to get number of blocks
    /// \param[in] nnz_prec       number of nonzeroes of matrix for ILU0, divide by dim*dim to get number of blocks
    /// \param[in] dim            size of block
    /// \param[in] vals           array of nonzeroes, each block is stored row-wise and contiguous, contains nnz values
    /// \param[in] rows           array of rowPointers, contains N/dim+1 values
    /// \param[in] cols           array of columnIndices, contains nnz values
    /// \param[in] vals_prec      array of nonzeroes for preconditioner
    /// \param[in] rows_prec      array of rowPointers for preconditioner
    /// \param[in] cols_prec      array of columnIndices for preconditioner
    /// \param[in] b              input vector, contains N values
    /// \param[in] wellContribs   WellContributions, to apply them separately, instead of adding them to matrix A
    /// \param[inout] res         summary of solver result
    /// \return                   status code
    SolverStatus solve_system(int N, int nnz, int dim, double *vals, int *rows, int *cols, double *b, WellContributions& wellContribs, BdaResult &res) override;

    /// Solve scalar linear system, for example a coarse system of an AMG preconditioner
    /// Data is already on the GPU
    // SolverStatus solve_system(BdaResult &res);

    /// Get result after linear solve, and peform postprocessing if necessary
    /// \param[inout] x          resulting x vector, caller must guarantee that x points to a valid array
    void get_result(double *x) override;

    /// Set OpenCL objects
    /// This class either creates them based on platformID and deviceID or receives them through this function
    /// \param[in] context   the opencl context to be used
    /// \param[in] queue     the opencl queue to be used
    void setOpencl(std::shared_ptr<cl::Context>& context, std::shared_ptr<cl::CommandQueue>& queue);

}; // end class openclSolverBackend

} // namespace Accelerator
} // namespace Opm

#endif


