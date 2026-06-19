// Standalone micro-benchmark for the threaded DILU preconditioners.
// Builds a representative 3D 7-point block matrix (block size 3, black-oil-like),
// then times preconditioner apply() at the current OMP_NUM_THREADS.
// No flow, no MPI: a dedicated solver benchmark per the project workflow.
//
// Usage: OMP_NUM_THREADS=N bench_dilu [nx ny nz reps]

#include <config.h>

#include <opm/simulators/linalg/DILU.hpp>
#include <opm/simulators/linalg/DILU2.hpp>

#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/bvector.hh>

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <vector>

#if HAVE_OPENMP
#include <omp.h>
#endif

namespace {

constexpr int bz = 3;
using Block = Dune::FieldMatrix<double, bz, bz>;
using Matrix = Dune::BCRSMatrix<Block>;
using Vector = Dune::BlockVector<Dune::FieldVector<double, bz>>;

// 7-point stencil index helpers on an nx*ny*nz grid (natural/lexicographic order).
struct Grid {
    int nx, ny, nz;
    std::size_t N() const { return std::size_t(nx) * ny * nz; }
    std::size_t id(int i, int j, int k) const { return std::size_t(k) * ny * nx + std::size_t(j) * nx + i; }
};

std::vector<std::array<int,3>> neighborOffsets()
{
    return {{-1,0,0},{1,0,0},{0,-1,0},{0,1,0},{0,0,-1},{0,0,1}};
}

Matrix build7pt(const Grid& g)
{
    const std::size_t N = g.N();
    // Count nonzeroes.
    std::size_t nnz = 0;
    std::vector<std::vector<std::size_t>> cols(N);
    for (int k = 0; k < g.nz; ++k)
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                const std::size_t r = g.id(i,j,k);
                cols[r].push_back(r);
                for (auto o : neighborOffsets()) {
                    int ni=i+o[0], nj=j+o[1], nk=k+o[2];
                    if (ni>=0&&ni<g.nx&&nj>=0&&nj<g.ny&&nk>=0&&nk<g.nz)
                        cols[r].push_back(g.id(ni,nj,nk));
                }
                nnz += cols[r].size();
            }

    Matrix A(N, N, nnz, Matrix::row_wise);
    for (auto row = A.createbegin(); row != A.createend(); ++row) {
        for (auto c : cols[row.index()]) {
            row.insert(c);
        }
    }
    // Fill: strongly diagonally dominant so DILU is well defined.
    for (std::size_t r = 0; r < N; ++r) {
        const int ncon = int(cols[r].size()) - 1;
        for (auto c : cols[r]) {
            Block b(0.0);
            if (c == r) {
                for (int d = 0; d < bz; ++d) b[d][d] = 2.0 * ncon + 4.0;
                // small off-diagonal block coupling
                b[0][1] = 0.3; b[1][0] = 0.2; b[1][2] = 0.1; b[2][1] = 0.15;
            } else {
                for (int d = 0; d < bz; ++d) b[d][d] = -1.0;
            }
            A[r][c] = b;
        }
    }
    return A;
}

template <class Prec>
double benchApply(Prec& prec, const Vector& d, int reps)
{
    Vector v(d.size());
    v = 0.0;
    // warmup
    for (int r = 0; r < 3; ++r) prec.apply(v, d);
    auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < reps; ++r) prec.apply(v, d);
    auto t1 = std::chrono::steady_clock::now();
    return std::chrono::duration<double>(t1 - t0).count();
}

} // namespace

int main(int argc, char** argv)
{
    Grid g{40, 40, 40};
    int reps = 200;
    if (argc >= 4) { g.nx = std::atoi(argv[1]); g.ny = std::atoi(argv[2]); g.nz = std::atoi(argv[3]); }
    if (argc >= 5) reps = std::atoi(argv[4]);

    int threads = 1;
#if HAVE_OPENMP
    threads = omp_get_max_threads();
#endif

    Matrix A = build7pt(g);
    Vector d(A.N());
    for (std::size_t i = 0; i < A.N(); ++i)
        for (int c = 0; c < bz; ++c) d[i][c] = 1.0 + 0.001 * double((i + c) % 97);

    std::printf("grid=%dx%dx%d N=%zu nnz=%zu blocks  threads=%d reps=%d\n",
                g.nx, g.ny, g.nz, A.N(), A.nonzeroes(), threads, reps);

    {
        auto t0 = std::chrono::steady_clock::now();
        Dune::MultithreadDILU<Matrix, Vector, Vector> prec(A);
        auto t1 = std::chrono::steady_clock::now();
        double apply_s = benchApply(prec, d, reps);
        std::printf("  DILU  (wavefront)  construct=%.4f s  apply_total=%.4f s  per_apply=%.6f s\n",
                    std::chrono::duration<double>(t1 - t0).count(), apply_s, apply_s / reps);
    }
    {
        auto t0 = std::chrono::steady_clock::now();
        Dune::MultithreadDILU2<Matrix, Vector, Vector> prec(A);
        auto t1 = std::chrono::steady_clock::now();
        double apply_s = benchApply(prec, d, reps);
        std::printf("  DILU2 (multicolor) construct=%.4f s  apply_total=%.4f s  per_apply=%.6f s  colors=%zu\n",
                    std::chrono::duration<double>(t1 - t0).count(), apply_s, apply_s / reps, prec.numColors());
    }
    return 0;
}
