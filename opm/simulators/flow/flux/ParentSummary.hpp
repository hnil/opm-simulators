/*
  Copyright 2026 Equinor ASA.

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

#ifndef OPM_PARENT_SUMMARY_HPP
#define OPM_PARENT_SUMMARY_HPP

#include <opm/input/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#include <opm/io/eclipse/FluxFile.hpp>

#include <filesystem>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace Opm {

class SummaryState;

/// Time series of a parent (full field) run's summary vectors.
///
/// A reduced USEFLUX run uses this to seed its SummaryState at every time
/// step, so that quantities it cannot compute itself -- most importantly the
/// rates of wells that fall outside the sector -- are still available to UDQ,
/// UDA and group control through the standard code paths.
///
/// The series is queried by time rather than by report step, because the
/// reduced run's time stepping will generally not line up with the parent's.
/// Interpolation is type dependent:
///
///   - Rate (and the discrete Mode/Count) quantities are held PIECEWISE
///     CONSTANT.  Sample times are interval END times, so the value stored at
///     t[i] applies throughout (t[i-1], t[i]].  When the sample series comes
///     from a .FLUX file those rate entries are time averages over exactly
///     that interval, which makes the hold integral preserving: it reproduces
///     the parent's production over the interval exactly.
///
///   - Everything else -- Total, Pressure, Ratio, ProdIndex -- is interpolated
///     LINEARLY, which is exact for cumulatives.
///
/// Queries outside the sampled range clamp to the first or last sample.
class ParentSummary
{
public:
    /// Build from the summary samples embedded in a .FLUX file.
    ///
    /// Returns nothing when the file carries no summary data.  Times in the
    /// flux file are already in seconds.
    static std::optional<ParentSummary>
    fromFluxFile(const EclIO::FluxFile::Data& data, const std::string& source);

    /// Build from a parent run's summary output.
    ///
    /// Supports .SMSPEC/.FSMSPEC (via EclIO::ESmry) and .ESMRY (via
    /// EclIO::ExtESmry).  The TIME vector of a summary file is expressed in
    /// DAYS and is converted to seconds here, so that both sources present the
    /// same time base.  Returns nothing when the file cannot be read.
    static std::optional<ParentSummary>
    fromSummaryFile(const std::filesystem::path& path);

    /// Whether \p key is present in the series.
    bool has(const std::string& key) const;

    /// Value of \p key at \p time, in seconds since the simulation start.
    ///
    /// Returns NaN when \p key is not present.
    double valueAt(const std::string& key, double time) const;

    /// Value of \p key over the interval (\p start, \p end], in seconds.
    ///
    /// For a rate this is the time average over the interval, taken across
    /// however many samples it spans. The parent's samples need not line up
    /// with the steps asking for them: a reduced run may step over several of
    /// them at once, and the rate at the end of the step then describes only
    /// the last sliver of it, where the average is what the parent actually
    /// produced over the step and what a cumulative built from it has to add
    /// up to. Everything else, and an empty interval, is valueAt(\p end).
    ///
    /// Intervals, or parts of them, where the value is NaN are left out of the
    /// average; NaN when nothing is left.
    double valueOver(const std::string& key, double start, double end) const;

    const std::vector<std::string>& keys() const { return this->keys_; }
    const std::vector<double>& times() const { return this->times_; }

    /// Human readable description of where the data came from, for logging.
    const std::string& sourceDescription() const { return this->source_; }

    bool empty() const { return this->keys_.empty() || this->times_.empty(); }

    /// Whether values of \p type are held piecewise constant rather than
    /// interpolated linearly.
    static bool isHeldConstant(SummaryConfigNode::Type type);

    /// Summary type of the keyword part of a composite key such as "WOPR:P1".
    static SummaryConfigNode::Type keyType(const std::string& key);

    /// Put \p value into \p state under \p key, a composite key such as
    /// "WOPR:P1" or "SOFR:P1:3", the way the summary evaluation stores it.
    static void assign(SummaryState& state, const std::string& key, double value);

    /// Region set, as FIPNUM or FIPxxx, and region ID of a region key such as
    /// "RPR__REC:6". Nothing for any other key.
    static std::optional<std::pair<std::string, int>> regionOf(const std::string& key);

private:
    struct Series
    {
        std::vector<double> values;
        SummaryConfigNode::Type type = SummaryConfigNode::Type::Undefined;
    };

    std::vector<std::string> keys_;

    /// Sample times in seconds, non-decreasing.
    std::vector<double> times_;

    std::unordered_map<std::string, Series> series_;
    std::string source_;
};

} // namespace Opm

#endif // OPM_PARENT_SUMMARY_HPP
