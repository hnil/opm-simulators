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
#ifndef OPM_MATRIX_SCALING_REPORT_HPP
#define OPM_MATRIX_SCALING_REPORT_HPP

#include <dune/common/dynmatrix.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>

#include <string>
#include <string_view>
#include <vector>

namespace Opm::detail {

//! \brief Per-variable and per-equation magnitude summary of a block matrix.
//! \details Diagnostic only. All vectors are indexed by the block index, so
//!          colMax[j] is the largest |A[.][.][i][j]| over every block entry.
template<class Scalar>
struct MatrixScalingStats
{
    int blockSize{};
    std::size_t numBlocks{};

    std::vector<Scalar> colMax, colMin;   //!< over the variable index
    std::vector<Scalar> rowMax, rowMin;   //!< over the equation index

    //! \brief blockSize*blockSize, row-major: max |A[i][j]| over all blocks.
    std::vector<Scalar> typeMatrix;

    //! \brief Ruiz diagonals equilibrating typeMatrix, and the resulting condition numbers.
    std::vector<Scalar> rowScale, colScale;
    Scalar condBefore{}, condAfter{};

    //! \brief Diagonal-block determinants, log10 histogram from 1e-60 upwards in decades.
    std::vector<std::size_t> detHistogram;
    Scalar detMin{}, detMax{};

    //! \brief Per-cell ideal column scale: 10th/50th/90th percentile of
    //!        log10(colScale) over all diagonal blocks, equilibrated one by one.
    //! \details p90-p10 is the number that decides static vs per-cell scaling.
    std::vector<Scalar> cellColScaleP10, cellColScaleP50, cellColScaleP90;
};

//! \brief Value of OPM_SCALING_REPORT, or 0 when unset. Diagnostic gate.
int scalingReportLevel();

//! \brief Two-sided Ruiz equilibration of a dense row-major n*n matrix.
//! \details Fills rowScale/colScale so that diag(rowScale)*M*diag(colScale) has
//!          unit infinity-norm rows and columns. Scales are rounded to powers of
//!          two when powerOfTwo, which keeps the scaling exact in IEEE-754.
template<class Scalar>
void ruizEquilibrate(const std::vector<Scalar>& matrix,
                     int n,
                     std::vector<Scalar>& rowScale,
                     std::vector<Scalar>& colScale,
                     bool powerOfTwo = true,
                     int maxIter = 20,
                     Scalar tol = Scalar(1e-3));

//! \brief Singular values of a dense row-major n*n matrix, descending.
//! \details One-sided Jacobi; n is at most numEq here, so cost is irrelevant.
template<class Scalar>
std::vector<Scalar> singularValues(const std::vector<Scalar>& matrix, int n);

//! \brief 2-norm condition number of a dense row-major n*n matrix.
template<class Scalar>
Scalar conditionNumber(const std::vector<Scalar>& matrix, int n);

//! \details Templated on the block, not on (Scalar, bs), so that it deduces for
//!          both Dune::FieldMatrix and Opm::MatrixBlock blocks.
template<class Block>
MatrixScalingStats<typename Block::field_type>
computeMatrixScalingStats(const Dune::BCRSMatrix<Block>& matrix);

//! \brief Format a stats struct as a multi-line human-readable report.
template<class Scalar>
std::string formatMatrixScalingStats(const MatrixScalingStats<Scalar>& stats,
                                     std::string_view tag);

//! \brief Compute and log a report. No-op unless OPM_SCALING_REPORT is set.
template<class Block>
void reportMatrixScaling(const Dune::BCRSMatrix<Block>& matrix,
                         std::string_view tag);

//! \brief Same, for a single dense block - used for the well D matrix.
template<class Block>
void reportDenseBlock(const Block& block, std::string_view tag);

//! \brief Report a runtime-sized dense block - the well D matrix is a DynamicMatrix.
//! \details Prints per-column and per-row magnitudes so the bhp column can be
//!          compared against the rate columns, plus cond and |det| before and
//!          after equilibration.
template<class Scalar>
void reportDynamicBlock(const Dune::DynamicMatrix<Scalar>& block, std::string_view tag);

//! \brief Value of OPM_SCALING_EQUILIBRATE: 0 off, 1 two-sided, 2 rows only, 3 columns only.
int scalingEquilibrateMode();

//! \brief Power-of-two Ruiz scales for the aggregated numEq x numEq type matrix.
//! \details One scale per variable and per equation for the whole matrix, i.e.
//!          what a static primary-variable/equation scaling could achieve.
//!          Powers of two so that applying and undoing is exact in IEEE-754.
template<class Block>
void computeEquilibration(const Dune::BCRSMatrix<Block>& matrix,
                          std::vector<typename Block::field_type>& rowScale,
                          std::vector<typename Block::field_type>& colScale);

//! \brief A[i][j] *= rowScale[i]*colScale[j] for every block. Inline: trivial, and
//!        avoids instantiating for every block type in use.
template<class Matrix, class Scalar>
void applyEquilibration(Matrix& matrix,
                        const std::vector<Scalar>& rowScale,
                        const std::vector<Scalar>& colScale)
{
    for (auto row = matrix.begin(); row != matrix.end(); ++row) {
        for (auto col = row->begin(); col != row->end(); ++col) {
            for (std::size_t i = 0; i < rowScale.size(); ++i) {
                for (std::size_t j = 0; j < colScale.size(); ++j) {
                    (*col)[i][j] *= rowScale[i] * colScale[j];
                }
            }
        }
    }
}

//! \brief v[.][i] *= scale[i] for every block of a block vector.
template<class Vector, class Scalar>
void scaleBlockVector(Vector& v, const std::vector<Scalar>& scale)
{
    for (auto& blk : v) {
        for (std::size_t i = 0; i < scale.size(); ++i) {
            blk[i] *= scale[i];
        }
    }
}

} // namespace Opm::detail

#endif // OPM_MATRIX_SCALING_REPORT_HPP
