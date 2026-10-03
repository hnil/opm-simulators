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

#include <opm/simulators/flow/flux/ParentSummary.hpp>

#include <opm/common/OpmLog/OpmLog.hpp>

#include <opm/input/eclipse/Schedule/SummaryState.hpp>
#include <opm/input/eclipse/Units/Units.hpp>

#include <opm/io/eclipse/SummaryNode.hpp>

#include <opm/io/eclipse/ESmry.hpp>
#include <opm/io/eclipse/ExtESmry.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iterator>
#include <limits>
#include <regex>
#include <string>
#include <utility>

namespace {

/// Lower-case file extension of \p path.
std::string extensionOf(const std::filesystem::path& path)
{
    auto ext = path.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(),
                   [](const unsigned char c) { return std::tolower(c); });
    return ext;
}

/// Copy a summary reader's vectors into a key -> values map.
///
/// \tparam Reader Either Opm::EclIO::ESmry or Opm::EclIO::ExtESmry.  The two
///   have compatible but not identical interfaces (ExtESmry's accessors are
///   non-const), so this is a template rather than a function taking a base.
template <class Reader>
bool loadFromReader(Reader&                                              reader,
                    std::vector<std::string>&                            keys,
                    std::vector<double>&                                 times,
                    std::unordered_map<std::string, std::vector<double>>& values)
{
    if (!reader.hasKey("TIME")) {
        return false;
    }

    const auto& timeVector = reader.get("TIME");
    if (timeVector.empty()) {
        return false;
    }

    times.clear();
    times.reserve(timeVector.size());
    for (const auto t : timeVector) {
        // Summary files express TIME in days; the rest of this class works in
        // seconds so that the flux-file and summary-file sources agree.
        times.push_back(Opm::unit::convert::from(static_cast<double>(t),
                                                 Opm::unit::day));
    }

    keys.clear();
    values.clear();

    for (const auto& key : reader.keywordList()) {
        const auto& vec = reader.get(key);
        if (vec.size() != timeVector.size()) {
            continue;
        }

        auto& dst = values[key];
        dst.assign(vec.begin(), vec.end());
        keys.push_back(key);
    }

    return !keys.empty();
}

} // namespace

namespace Opm {

bool ParentSummary::isHeldConstant(const SummaryConfigNode::Type type)
{
    // A rate is a property of an interval, not of an instant, so blending two
    // neighbouring rates would invent a value the parent never produced.  Mode
    // and Count are discrete, so interpolating them is meaningless.
    return (type == SummaryConfigNode::Type::Rate)
        || (type == SummaryConfigNode::Type::Mode)
        || (type == SummaryConfigNode::Type::Count);
}

SummaryConfigNode::Type ParentSummary::keyType(const std::string& key)
{
    const auto colon = key.find(':');
    const auto keyword = (colon == std::string::npos)
        ? key
        : key.substr(0, colon);

    if (keyword.empty()) {
        return SummaryConfigNode::Type::Undefined;
    }

    return parseKeywordType(keyword);
}

void ParentSummary::assign(SummaryState& state, const std::string& key, const double value)
{
    const auto firstColon = key.find(':');
    const auto secondColon = (firstColon == std::string::npos)
        ? std::string::npos
        : key.find(':', firstColon + 1);

    const auto keyword = key.substr(0, firstColon);

    switch (EclIO::SummaryNode::category_from_keyword(keyword)) {
    case EclIO::SummaryNode::Category::Well:
        if (firstColon != std::string::npos) {
            state.update_well_var(key.substr(firstColon + 1), keyword, value);
        }
        break;

    case EclIO::SummaryNode::Category::Group:
        if (firstColon != std::string::npos) {
            state.update_group_var(key.substr(firstColon + 1), keyword, value);
        }
        break;

    case EclIO::SummaryNode::Category::Connection:
    case EclIO::SummaryNode::Category::Completion:
        if (secondColon != std::string::npos) {
            state.update_conn_var(key.substr(firstColon + 1, secondColon - firstColon - 1),
                                  keyword,
                                  static_cast<std::size_t>(std::stoul(key.substr(secondColon + 1))),
                                  value);
        }
        break;

    case EclIO::SummaryNode::Category::Segment:
        if (secondColon != std::string::npos) {
            state.update_segment_var(key.substr(firstColon + 1, secondColon - firstColon - 1),
                                     keyword,
                                     static_cast<std::size_t>(std::stoul(key.substr(secondColon + 1))),
                                     value);
        }
        break;

    case EclIO::SummaryNode::Category::Region:
        // Held apart from other values, keyed by region set and number, so
        // that a UDQ can find them.
        if (const auto region = regionOf(key); region.has_value()) {
            state.update_region_var(region->first, keyword, region->second, value);
        }
        else {
            state.set(key, value);
        }
        break;

    default:
        state.set(key, value);
        break;
    }
}

std::optional<std::pair<std::string, int>> ParentSummary::regionOf(const std::string& key)
{
    // RPR, RPR__REC, RGPR_EQL: the quantity, then optionally the region set's
    // name padded to eight characters with underscores.
    static const auto pattern = std::regex { R"((R[A-Z]{2,4})_{0,2}([A-Z0-9]{1,3})?:([0-9]+))" };

    auto pieces = std::smatch {};
    if (!std::regex_match(key, pieces, pattern)) {
        return std::nullopt;
    }

    const auto set = pieces[2].matched ? "FIP" + pieces[2].str() : std::string { "FIPNUM" };
    return std::make_pair(set, std::stoi(pieces[3].str()));
}

std::optional<ParentSummary>
ParentSummary::fromFluxFile(const EclIO::FluxFile::Data& data, const std::string& source)
{
    if (data.summaryKeys.empty() || data.summarySamples.empty()) {
        return std::nullopt;
    }

    ParentSummary summary;
    summary.keys_ = data.summaryKeys;
    summary.source_ = source;

    summary.times_.reserve(data.summarySamples.size());
    for (const auto& sample : data.summarySamples) {
        summary.times_.push_back(sample.time);
    }

    for (std::size_t k = 0; k < data.summaryKeys.size(); ++k) {
        auto& series = summary.series_[data.summaryKeys[k]];
        series.type = keyType(data.summaryKeys[k]);
        series.values.reserve(data.summarySamples.size());

        for (const auto& sample : data.summarySamples) {
            series.values.push_back(sample.values[k]);
        }
    }

    return summary;
}

std::optional<ParentSummary>
ParentSummary::fromSummaryFile(const std::filesystem::path& path)
{
    if (!std::filesystem::exists(path)) {
        return std::nullopt;
    }

    auto keys = std::vector<std::string>{};
    auto times = std::vector<double>{};
    auto values = std::unordered_map<std::string, std::vector<double>>{};

    const auto ext = extensionOf(path);

    try {
        if (ext == ".esmry") {
            EclIO::ExtESmry reader(path.string());
            if (!loadFromReader(reader, keys, times, values)) {
                return std::nullopt;
            }
        }
        else {
            EclIO::ESmry reader(path.string());
            if (!loadFromReader(reader, keys, times, values)) {
                return std::nullopt;
            }
        }
    }
    catch (const std::exception& e) {
        OpmLog::warning(fmt::format("Failed to read parent summary file '{}': {}",
                                    path.string(), e.what()));
        return std::nullopt;
    }

    ParentSummary summary;
    summary.keys_ = std::move(keys);
    summary.times_ = std::move(times);
    summary.source_ = path.string();

    for (auto& [key, vec] : values) {
        auto& series = summary.series_[key];
        series.type = keyType(key);
        series.values = std::move(vec);
    }

    return summary;
}

bool ParentSummary::has(const std::string& key) const
{
    return this->series_.find(key) != this->series_.end();
}

double ParentSummary::valueAt(const std::string& key, const double time) const
{
    const auto it = this->series_.find(key);
    if ((it == this->series_.end()) || this->times_.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const auto& series = it->second;
    if (series.values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // First sample at or after 'time'.
    const auto pos = std::lower_bound(this->times_.begin(), this->times_.end(), time);

    if (pos == this->times_.begin()) {
        return series.values.front();
    }

    if (pos == this->times_.end()) {
        return series.values.back();
    }

    const auto hi = static_cast<std::size_t>(std::distance(this->times_.begin(), pos));
    const auto lo = hi - 1;

    if (isHeldConstant(series.type)) {
        // Sample times are interval end times, so the value stored at times_[hi]
        // is the one that applies throughout (times_[lo], times_[hi]].
        return series.values[hi];
    }

    const auto t0 = this->times_[lo];
    const auto t1 = this->times_[hi];
    if (!(t1 > t0)) {
        return series.values[hi];
    }

    const auto weight = (time - t0) / (t1 - t0);
    return series.values[lo] + weight * (series.values[hi] - series.values[lo]);
}

double ParentSummary::valueOver(const std::string& key, const double start, const double end) const
{
    const auto it = this->series_.find(key);
    if ((it == this->series_.end())
        || (it->second.type != SummaryConfigNode::Type::Rate)
        || !(end > start))
    {
        return this->valueAt(key, end);
    }

    const auto& values = it->second.values;
    const auto& times = this->times_;
    if (values.empty() || times.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Sample k holds for (times[k-1], times[k]]; the first one reaches back
    // indefinitely and the last one forward, as in valueAt().
    double integral = 0.0;
    double covered = 0.0;
    const auto add = [&](const double value, const double from, const double to)
    {
        const auto length = std::min(to, end) - std::max(from, start);
        if ((length > 0.0) && std::isfinite(value)) {
            integral += value * length;
            covered += length;
        }
    };

    const auto lowest = std::numeric_limits<double>::lowest();
    const auto highest = std::numeric_limits<double>::max();

    for (std::size_t k = 0; k < values.size() && k < times.size(); ++k) {
        const auto from = (k == 0) ? lowest : times[k - 1];
        if (from >= end) {
            break;
        }
        add(values[k], from, times[k]);
    }
    add(values.back(), times.back(), highest);

    return (covered > 0.0)
        ? integral / covered
        : std::numeric_limits<double>::quiet_NaN();
}

} // namespace Opm
