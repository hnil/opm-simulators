// Standalone micro-benchmark for the (OpenMP) Hypre BoomerAMG preconditioner,
// representing the CPR pressure stage. Builds a scalar (1x1) 3D 7-point Poisson
// matrix and times HyprePreconditioner::apply() at the current OMP_NUM_THREADS.
// Dedicated solver benchmark (no flow). Usage: OMP_NUM_THREADS=N bench_hypre [nx ny nz reps]

#include <config.h>

// Dune ISTL headers first: HyprePreconditioner.hpp pulls in <_hypre_utilities.h>
// which #defines DATA (clashing with dune matrixmarket's enum), so anything that
// transitively includes matrixmarket.hh must be processed before it.
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/paamg/pinfo.hh>

#include <opm/simulators/linalg/PropertyTree.hpp>
#include <opm/simulators/linalg/HyprePreconditioner.hpp>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>

#if HAVE_OPENMP
#include <omp.h>
#endif

namespace {

using Block = Dune::FieldMatrix<double, 1, 1>;
using Matrix = Dune::BCRSMatrix<Block>;
using Vector = Dune::BlockVector<Dune::FieldVector<double, 1>>;

struct Grid { int nx, ny, nz; std::size_t N() const { return std::size_t(nx)*ny*nz; }
    std::size_t id(int i,int j,int k) const { return std::size_t(k)*ny*nx + std::size_t(j)*nx + i; } };

Matrix buildPoisson(const Grid& g)
{
    const std::size_t N = g.N();
    const int off[6][3] = {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
    std::vector<std::vector<std::size_t>> cols(N);
    std::size_t nnz = 0;
    for (int k=0;k<g.nz;++k) for (int j=0;j<g.ny;++j) for (int i=0;i<g.nx;++i) {
        const std::size_t r = g.id(i,j,k);
        cols[r].push_back(r);
        for (auto& o : off) { int ni=i+o[0],nj=j+o[1],nk=k+o[2];
            if (ni>=0&&ni<g.nx&&nj>=0&&nj<g.ny&&nk>=0&&nk<g.nz) cols[r].push_back(g.id(ni,nj,nk)); }
        nnz += cols[r].size();
    }
    Matrix A(N, N, nnz, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row)
        for (auto c : cols[row.index()]) row.insert(c);
    for (std::size_t r=0;r<N;++r) {
        const int ncon = int(cols[r].size())-1;
        for (auto c : cols[r]) A[r][c] = (c==r) ? double(ncon) + 0.1 : -1.0;
    }
    return A;
}

} // namespace

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    HYPRE_Initialize();   // required by Hypre 3.x before any HYPRE call (flow does this in Main.cpp)
    Grid g{40,40,40};
    int reps = 100;
    if (argc>=4){ g.nx=std::atoi(argv[1]); g.ny=std::atoi(argv[2]); g.nz=std::atoi(argv[3]); }
    if (argc>=5) reps=std::atoi(argv[4]);
    int threads = 1;
#if HAVE_OPENMP
    threads = omp_get_max_threads();
#endif

    Matrix A = buildPoisson(g);
    Vector d(A.N()), v(A.N());
    for (std::size_t i=0;i<A.N();++i){ d[i]=1.0; v[i]=0.0; }

    Opm::PropertyTree prm;
    prm.put("use_gpu", false);
    prm.put("print_level", 0);
    prm.put("max_iter", 1);
    prm.put("tolerance", 0.0);
    prm.put("strong_threshold", 0.5);
    prm.put("agg_trunc_factor", 0.3);
    prm.put("interp_type", 6);
    prm.put("max_levels", 15);
    prm.put("coarsen_type", 10);
    prm.put("relax_type", 18);
    prm.put("agg_num_levels", 1);
    prm.put("agg_interp_type", 4);

    Dune::Amg::SequentialInformation comm;
    auto t0 = std::chrono::steady_clock::now();
    Opm::linalg::HyprePreconditioner<Matrix, Vector, Vector, Dune::Amg::SequentialInformation> prec(A, prm, comm);
    auto t1 = std::chrono::steady_clock::now();

    for (int r=0;r<3;++r){ v=0.0; prec.apply(v, d); }   // warmup
    auto a0 = std::chrono::steady_clock::now();
    for (int r=0;r<reps;++r){ v=0.0; prec.apply(v, d); }
    auto a1 = std::chrono::steady_clock::now();

    double setup = std::chrono::duration<double>(t1-t0).count();
    double apply_s = std::chrono::duration<double>(a1-a0).count();
    std::printf("grid=%dx%dx%d N=%zu threads=%d reps=%d  hypre_setup=%.4f s  apply_total=%.4f s  per_apply=%.6f s\n",
                g.nx,g.ny,g.nz, A.N(), threads, reps, setup, apply_s, apply_s/reps);
    MPI_Finalize();
    return 0;
}
