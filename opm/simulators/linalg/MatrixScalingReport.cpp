/*
  Copyright 2026 SINTEF Digital, Mathematics and Cybernetics.

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

#include "config.h"

#include <opm/simulators/linalg/MatrixScalingReport.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>
#include <opm/simulators/linalg/matrixblock.hh>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <limits>

namespace {

// Determinant histogram covers 1e-60 .. 1e+20 in decades.
constexpr int detBinMin = -60;
constexpr int detBinMax = 20;
constexpr int numDetBins = detBinMax - detBinMin + 1;

//! LU with partial pivoting on a copy. Own implementation rather than
//! FieldMatrix::determinant() so that a singular block cannot throw.
template<class Scalar>
Scalar denseDeterminant(std::vector<Scalar> m, const int n)
{
    Scalar det{1};
    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(m[i * n + k]) > std::abs(m[piv * n + k])) {
                piv = i;
            }
        }
        if (m[piv * n + k] == Scalar{0}) {
            return Scalar{0};
        }
        if (piv != k) {
            for (int j = 0; j < n; ++j) {
                std::swap(m[k * n + j], m[piv * n + j]);
            }
            det = -det;
        }
        det *= m[k * n + k];
        for (int i = k + 1; i < n; ++i) {
            const Scalar f = m[i * n + k] / m[k * n + k];
            for (int j = k; j < n; ++j) {
                m[i * n + j] -= f * m[k * n + j];
            }
        }
    }
    return det;
}

} // anonymous namespace

namespace Opm::detail {

int scalingReportLevel()
{
    static const int level = []() {
        const char* env = std::getenv("OPM_SCALING_REPORT");
        return env == nullptr ? 0 : std::atoi(env);
    }();
    return level;
}

template<class Scalar>
void ruizEquilibrate(const std::vector<Scalar>& matrix,
                     const int n,
                     std::vector<Scalar>& rowScale,
                     std::vector<Scalar>& colScale,
                     const bool powerOfTwo,
                     const int maxIter,
                     const Scalar tol)
{
    std::vector<Scalar> m(matrix);
    rowScale.assign(n, Scalar{1});
    colScale.assign(n, Scalar{1});

    for (int iter = 0; iter < maxIter; ++iter) {
        Scalar deviation{0};

        for (int i = 0; i < n; ++i) {
            Scalar mx{0};
            for (int j = 0; j < n; ++j) {
                mx = std::max(mx, std::abs(m[i * n + j]));
            }
            if (!(mx > Scalar{0})) {
                continue;
            }
            deviation = std::max(deviation, std::abs(Scalar{1} - mx));
            const Scalar r = std::sqrt(mx);
            rowScale[i] /= r;
            for (int j = 0; j < n; ++j) {
                m[i * n + j] /= r;
            }
        }

        for (int j = 0; j < n; ++j) {
            Scalar mx{0};
            for (int i = 0; i < n; ++i) {
                mx = std::max(mx, std::abs(m[i * n + j]));
            }
            if (!(mx > Scalar{0})) {
                continue;
            }
            deviation = std::max(deviation, std::abs(Scalar{1} - mx));
            const Scalar c = std::sqrt(mx);
            colScale[j] /= c;
            for (int i = 0; i < n; ++i) {
                m[i * n + j] /= c;
            }
        }

        if (deviation < tol) {
            break;
        }
    }

    // Powers of two keep the scaling exact in IEEE-754, so an equilibrated
    // system can be compared bit-for-bit against the unscaled one.
    if (powerOfTwo) {
        const auto snap = [](const Scalar s) {
            return s > Scalar{0} ? std::exp2(std::round(std::log2(s))) : Scalar{1};
        };
        std::transform(rowScale.begin(), rowScale.end(), rowScale.begin(), snap);
        std::transform(colScale.begin(), colScale.end(), colScale.begin(), snap);
    }
}

template<class Scalar>
std::vector<Scalar> singularValues(const std::vector<Scalar>& matrix, const int n)
{
    std::vector<Scalar> b(matrix);
    const Scalar eps = std::numeric_limits<Scalar>::epsilon();

    for (int sweep = 0; sweep < 30; ++sweep) {
        Scalar off{0};
        for (int p = 0; p < n - 1; ++p) {
            for (int q = p + 1; q < n; ++q) {
                Scalar alpha{0}, beta{0}, gamma{0};
                for (int i = 0; i < n; ++i) {
                    alpha += b[i * n + p] * b[i * n + p];
                    beta += b[i * n + q] * b[i * n + q];
                    gamma += b[i * n + p] * b[i * n + q];
                }
                const Scalar norm = std::sqrt(alpha * beta);
                if (gamma == Scalar{0} || !(std::abs(gamma) > eps * norm)) {
                    continue;
                }
                off = std::max(off, std::abs(gamma) / norm);

                const Scalar zeta = (beta - alpha) / (Scalar{2} * gamma);
                const Scalar t = (zeta >= Scalar{0} ? Scalar{1} : Scalar{-1})
                    / (std::abs(zeta) + std::sqrt(Scalar{1} + zeta * zeta));
                const Scalar c = Scalar{1} / std::sqrt(Scalar{1} + t * t);
                const Scalar s = c * t;

                for (int i = 0; i < n; ++i) {
                    const Scalar bp = b[i * n + p];
                    const Scalar bq = b[i * n + q];
                    b[i * n + p] = c * bp - s * bq;
                    b[i * n + q] = s * bp + c * bq;
                }
            }
        }
        if (!(off > Scalar{1e-14})) {
            break;
        }
    }

    std::vector<Scalar> sv(n);
    for (int j = 0; j < n; ++j) {
        Scalar s{0};
        for (int i = 0; i < n; ++i) {
            s += b[i * n + j] * b[i * n + j];
        }
        sv[j] = std::sqrt(s);
    }
    std::sort(sv.begin(), sv.end(), std::greater<Scalar>{});
    return sv;
}

template<class Scalar>
Scalar conditionNumber(const std::vector<Scalar>& matrix, const int n)
{
    const auto sv = singularValues(matrix, n);
    if (sv.empty() || !(sv.back() > Scalar{0})) {
        return std::numeric_limits<Scalar>::infinity();
    }
    return sv.front() / sv.back();
}

template<class Block>
MatrixScalingStats<typename Block::field_type>
computeMatrixScalingStats(const Dune::BCRSMatrix<Block>& matrix)
{
    using Scalar = typename Block::field_type;
    constexpr int bs = Block::rows;

    MatrixScalingStats<Scalar> stats;
    stats.blockSize = bs;
    stats.colMax.assign(bs, Scalar{0});
    stats.colMin.assign(bs, std::numeric_limits<Scalar>::max());
    stats.rowMax.assign(bs, Scalar{0});
    stats.rowMin.assign(bs, std::numeric_limits<Scalar>::max());
    stats.typeMatrix.assign(bs * bs, Scalar{0});
    stats.detHistogram.assign(numDetBins, 0);
    stats.detMin = std::numeric_limits<Scalar>::max();
    stats.detMax = Scalar{0};

    std::vector<Scalar> diag(bs * bs);

    // log10 of the per-cell ideal column scale, one sample per diagonal block.
    std::vector<std::vector<Scalar>> cellColScale(bs);
    std::vector<Scalar> cellRow, cellCol;

    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            ++stats.numBlocks;
            const auto& block = *col;

            for (int i = 0; i < bs; ++i) {
                for (int j = 0; j < bs; ++j) {
                    const Scalar v = std::abs(block[i][j]);
                    if (!(v > Scalar{0})) {
                        continue;
                    }
                    stats.colMax[j] = std::max(stats.colMax[j], v);
                    stats.colMin[j] = std::min(stats.colMin[j], v);
                    stats.rowMax[i] = std::max(stats.rowMax[i], v);
                    stats.rowMin[i] = std::min(stats.rowMin[i], v);
                    stats.typeMatrix[i * bs + j] = std::max(stats.typeMatrix[i * bs + j], v);
                }
            }

            if (col.index() == row.index()) {
                for (int i = 0; i < bs; ++i) {
                    for (int j = 0; j < bs; ++j) {
                        diag[i * bs + j] = block[i][j];
                    }
                }
                const Scalar det = std::abs(denseDeterminant(diag, bs));
                stats.detMin = std::min(stats.detMin, det);
                stats.detMax = std::max(stats.detMax, det);
                const int bin = (det > Scalar{0})
                    ? static_cast<int>(std::floor(std::log10(det))) - detBinMin
                    : 0;
                ++stats.detHistogram[static_cast<std::size_t>(std::clamp(bin, 0, numDetBins - 1))];

                // Equilibrate this cell on its own; the spread of the resulting
                // column scales across cells is what decides static vs per-cell.
                // Not power-of-two here - we want the ideal value, not a snapped one.
                ruizEquilibrate(diag, bs, cellRow, cellCol, /*powerOfTwo=*/false);
                for (int j = 0; j < bs; ++j) {
                    if (cellCol[j] > Scalar{0}) {
                        cellColScale[j].push_back(std::log10(cellCol[j]));
                    }
                }
            }
        }
    }

    stats.condBefore = conditionNumber(stats.typeMatrix, bs);

    ruizEquilibrate(stats.typeMatrix, bs, stats.rowScale, stats.colScale);

    std::vector<Scalar> equilibrated(stats.typeMatrix);
    for (int i = 0; i < bs; ++i) {
        for (int j = 0; j < bs; ++j) {
            equilibrated[i * bs + j] *= stats.rowScale[i] * stats.colScale[j];
        }
    }
    stats.condAfter = conditionNumber(equilibrated, bs);

    stats.cellColScaleP10.assign(bs, Scalar{0});
    stats.cellColScaleP50.assign(bs, Scalar{0});
    stats.cellColScaleP90.assign(bs, Scalar{0});
    for (int j = 0; j < bs; ++j) {
        auto& s = cellColScale[j];
        if (s.empty()) {
            continue;
        }
        std::sort(s.begin(), s.end());
        const auto pick = [&s](const double q) {
            const auto k = static_cast<std::size_t>(q * static_cast<double>(s.size() - 1));
            return s[k];
        };
        stats.cellColScaleP10[j] = pick(0.10);
        stats.cellColScaleP50[j] = pick(0.50);
        stats.cellColScaleP90[j] = pick(0.90);
    }

    return stats;
}

template<class Scalar>
std::string formatMatrixScalingStats(const MatrixScalingStats<Scalar>& stats,
                                     std::string_view tag)
{
    const int n = stats.blockSize;
    const auto decades = [](const Scalar hi, const Scalar lo) {
        return (hi > Scalar{0} && lo > Scalar{0} && lo <= hi)
            ? std::log10(hi / lo) : Scalar{0};
    };

    std::string out = fmt::format("scaling-report [{}] blocks={} blockSize={}\n",
                                  tag, stats.numBlocks, n);

    out += "  idx |     col max |     col min | decades |     row max |     row min | decades\n";
    for (int k = 0; k < n; ++k) {
        out += fmt::format("  {:3d} | {:11.4e} | {:11.4e} | {:7.2f} | {:11.4e} | {:11.4e} | {:7.2f}\n",
                           k,
                           stats.colMax[k], stats.colMin[k], decades(stats.colMax[k], stats.colMin[k]),
                           stats.rowMax[k], stats.rowMin[k], decades(stats.rowMax[k], stats.rowMin[k]));
    }

    out += fmt::format("  type-matrix cond {:.4e} -> {:.4e} after Ruiz\n",
                       stats.condBefore, stats.condAfter);

    out += "  row scale:";
    for (const auto s : stats.rowScale) {
        out += fmt::format(" {:.4e}", s);
    }
    out += "\n  col scale:";
    for (const auto s : stats.colScale) {
        out += fmt::format(" {:.4e}", s);
    }

    out += "\n  per-cell col scale log10 p10/p50/p90:";
    for (int k = 0; k < n; ++k) {
        out += fmt::format(" [{}] {:.2f}/{:.2f}/{:.2f}", k,
                           stats.cellColScaleP10[k],
                           stats.cellColScaleP50[k],
                           stats.cellColScaleP90[k]);
    }

    out += fmt::format("\n  |det| diag blocks: min {:.4e} max {:.4e}\n",
                       stats.detMin, stats.detMax);
    out += "  |det| decades:";
    for (std::size_t b = 0; b < stats.detHistogram.size(); ++b) {
        if (stats.detHistogram[b] > 0) {
            out += fmt::format(" 1e{}:{}", detBinMin + static_cast<int>(b), stats.detHistogram[b]);
        }
    }

    return out;
}

template<class Scalar>
void reportDynamicBlock(const Dune::DynamicMatrix<Scalar>& block, std::string_view tag)
{
    if (scalingReportLevel() <= 0) {
        return;
    }

    const int n = static_cast<int>(block.N());
    if (n <= 0 || static_cast<int>(block.M()) != n) {
        return;
    }

    std::vector<Scalar> m(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            m[i * n + j] = block[i][j];
        }
    }

    std::vector<Scalar> rowScale, colScale;
    ruizEquilibrate(m, n, rowScale, colScale);

    std::vector<Scalar> equilibrated(m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            equilibrated[i * n + j] *= rowScale[i] * colScale[j];
        }
    }

    std::string out = fmt::format("scaling-report [{}] dyn {}x{} cond {:.4e} -> {:.4e}"
                                  " |det| {:.4e} -> {:.4e}\n",
                                  tag, n, n,
                                  conditionNumber(m, n), conditionNumber(equilibrated, n),
                                  std::abs(denseDeterminant(m, n)),
                                  std::abs(denseDeterminant(equilibrated, n)));

    out += "  colmax:";
    for (int j = 0; j < n; ++j) {
        Scalar mx{0};
        for (int i = 0; i < n; ++i) {
            mx = std::max(mx, std::abs(m[i * n + j]));
        }
        out += fmt::format(" {:.4e}", mx);
    }
    out += "\n  rowmax:";
    for (int i = 0; i < n; ++i) {
        Scalar mx{0};
        for (int j = 0; j < n; ++j) {
            mx = std::max(mx, std::abs(m[i * n + j]));
        }
        out += fmt::format(" {:.4e}", mx);
    }
    out += "\n  colscale:";
    for (const auto s : colScale) {
        out += fmt::format(" {:.4e}", s);
    }
    OpmLog::info(out);
}

int scalingEquilibrateMode()
{
    static const int mode = []() {
        const char* env = std::getenv("OPM_SCALING_EQUILIBRATE");
        return env == nullptr ? 0 : std::atoi(env);
    }();
    return mode;
}

template<class Block>
void computeEquilibration(const Dune::BCRSMatrix<Block>& matrix,
                          std::vector<typename Block::field_type>& rowScale,
                          std::vector<typename Block::field_type>& colScale)
{
    using Scalar = typename Block::field_type;
    constexpr int bs = Block::rows;

    // Only the type matrix, not the per-cell work computeMatrixScalingStats does.
    std::vector<Scalar> typeMatrix(bs * bs, Scalar{0});
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            for (int i = 0; i < bs; ++i) {
                for (int j = 0; j < bs; ++j) {
                    const Scalar v = std::abs((*col)[i][j]);
                    typeMatrix[i * bs + j] = std::max(typeMatrix[i * bs + j], v);
                }
            }
        }
    }

    ruizEquilibrate(typeMatrix, bs, rowScale, colScale);

    const int mode = scalingEquilibrateMode();
    if (mode == 2) {
        colScale.assign(bs, Scalar{1});
    }
    else if (mode == 3) {
        rowScale.assign(bs, Scalar{1});
    }
}

template<class Block>
void reportMatrixScaling(const Dune::BCRSMatrix<Block>& matrix,
                         std::string_view tag)
{
    if (scalingReportLevel() <= 0) {
        return;
    }
    const auto stats = computeMatrixScalingStats(matrix);
    OpmLog::info(formatMatrixScalingStats(stats, tag));
}

template<class Block>
void reportDenseBlock(const Block& block, std::string_view tag)
{
    if (scalingReportLevel() <= 0) {
        return;
    }

    using Scalar = typename Block::field_type;
    constexpr int bs = Block::rows;

    std::vector<Scalar> m(bs * bs);
    for (int i = 0; i < bs; ++i) {
        for (int j = 0; j < bs; ++j) {
            m[i * bs + j] = block[i][j];
        }
    }

    std::vector<Scalar> rowScale, colScale;
    ruizEquilibrate(m, bs, rowScale, colScale);

    std::vector<Scalar> equilibrated(m);
    for (int i = 0; i < bs; ++i) {
        for (int j = 0; j < bs; ++j) {
            equilibrated[i * bs + j] *= rowScale[i] * colScale[j];
        }
    }

    std::string out = fmt::format("scaling-report [{}] dense {}x{} cond {:.4e} -> {:.4e}"
                                  " |det| {:.4e} -> {:.4e}\n",
                                  tag, bs, bs,
                                  conditionNumber(m, bs), conditionNumber(equilibrated, bs),
                                  std::abs(denseDeterminant(m, bs)),
                                  std::abs(denseDeterminant(equilibrated, bs)));
    out += "  col scale:";
    for (const auto s : colScale) {
        out += fmt::format(" {:.4e}", s);
    }
    OpmLog::info(out);
}

#define INSTANTIATE_SCALAR(T)                                                   \
    template void ruizEquilibrate(const std::vector<T>&, int,                   \
                                  std::vector<T>&, std::vector<T>&,             \
                                  bool, int, T);                                \
    template std::vector<T> singularValues(const std::vector<T>&, int);         \
    template T conditionNumber(const std::vector<T>&, int);                     \
    template std::string formatMatrixScalingStats(const MatrixScalingStats<T>&, \
                                                  std::string_view);            \
    template void reportDynamicBlock(const Dune::DynamicMatrix<T>&, std::string_view);

// Both block types occur: the reservoir matrix uses Opm::MatrixBlock, the well
// D matrices use plain Dune::FieldMatrix.
// BLOCKTPL is passed unapplied: the commas of its template arguments would
// otherwise be read as macro argument separators.
#define INSTANTIATE_FOR_BLOCK(BLOCKTPL, T, bs)                                    \
    template MatrixScalingStats<T>                                                \
    computeMatrixScalingStats(const Dune::BCRSMatrix<BLOCKTPL<T,bs,bs>>&);        \
    template void reportMatrixScaling(const Dune::BCRSMatrix<BLOCKTPL<T,bs,bs>>&, \
                                      std::string_view);                          \
    template void reportDenseBlock(const BLOCKTPL<T,bs,bs>&, std::string_view);   \
    template void computeEquilibration(const Dune::BCRSMatrix<BLOCKTPL<T,bs,bs>>&, \
                                       std::vector<T>&, std::vector<T>&);

#define INSTANTIATE_BLOCK(T, bs)                      \
    INSTANTIATE_FOR_BLOCK(Dune::FieldMatrix, T, bs)   \
    INSTANTIATE_FOR_BLOCK(Opm::MatrixBlock, T, bs)

#define INSTANTIATE_ALL(T)  \
    INSTANTIATE_SCALAR(T)   \
    INSTANTIATE_BLOCK(T, 1) \
    INSTANTIATE_BLOCK(T, 2) \
    INSTANTIATE_BLOCK(T, 3) \
    INSTANTIATE_BLOCK(T, 4) \
    INSTANTIATE_BLOCK(T, 5) \
    INSTANTIATE_BLOCK(T, 6) \
    INSTANTIATE_BLOCK(T, 7)

INSTANTIATE_ALL(double)

#if FLOW_INSTANTIATE_FLOAT
INSTANTIATE_ALL(float)
#endif

} // namespace Opm::detail
