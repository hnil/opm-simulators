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

#ifndef OPM_NEUMAIER_SUM_HPP
#define OPM_NEUMAIER_SUM_HPP

#include <cmath>

namespace Opm {

/// Compensated (Neumaier/Kahan-Babuska) accumulator.
///
/// Quantities that feed convergence decisions -- the signed mass-balance
/// residual sum in particular -- suffer catastrophic cancellation, so their
/// low-order digits depend on summation order and on compiler instruction
/// selection (FMA contraction, vectorization). Since accept/reject decisions
/// amplify such round-off differences into diverging solution trajectories,
/// these sums are accumulated with first-order error compensation, which
/// makes the result independent of those effects for all practical purposes.
///
/// The compensated add costs a handful of flops per element and is branch-free
/// except for one fabs comparison; the loops using it are a negligible part of
/// a nonlinear iteration.
template <typename Scalar>
class NeumaierSum
{
public:
    NeumaierSum() = default;

    explicit NeumaierSum(const Scalar value)
        : sum_(value)
    {}

    void add(const Scalar value)
    {
        const Scalar t = sum_ + value;
        if (std::fabs(sum_) >= std::fabs(value)) {
            comp_ += (sum_ - t) + value;
        } else {
            comp_ += (value - t) + sum_;
        }
        sum_ = t;
    }

    NeumaierSum& operator+=(const Scalar value)
    {
        this->add(value);
        return *this;
    }

    /// Final compensated value.
    Scalar value() const
    {
        return sum_ + comp_;
    }

private:
    Scalar sum_{0};
    Scalar comp_{0};
};

} // namespace Opm

#endif // OPM_NEUMAIER_SUM_HPP
