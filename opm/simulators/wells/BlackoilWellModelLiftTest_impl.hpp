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

#ifndef OPM_BLACKOILWELLMODEL_LIFT_TEST_IMPL_HEADER_INCLUDED
#define OPM_BLACKOILWELLMODEL_LIFT_TEST_IMPL_HEADER_INCLUDED

// Does a producer lift at a wellhead pressure, and at what rate? One answer for every caller
// (the controller's status decisions, the facility check, the reference), always determined.

#include <opm/common/OpmLog/OpmLog.hpp>

#include <fmt/format.h>

#include <algorithm>
#include <cstdlib>
#include <optional>
#include <vector>

namespace Opm {

template<typename TypeTag>
typename BlackoilWellModel<TypeTag>::WellInterfacePtr
BlackoilWellModel<TypeTag>::
makeProbe_(const std::string& name, DeferredLogger& deferred_logger) const
{
    const int step = simulator_.episodeIndex();
    auto* self = const_cast<BlackoilWellModel<TypeTag>*>(this);
    auto probe = self->createWellForWellTest(name, step, deferred_logger);
    const auto live = std::find_if(well_container_.begin(), well_container_.end(),
                                   [&name](const auto& wp) { return wp->name() == name; });
    probe->init(depth_, gravity_, B_avg_, true);
    probe->setWellEfficiencyFactor(live != well_container_.end() ? (*live)->wellEfficiencyFactor() : Scalar(1));
    probe->setVFPProperties(this->vfp_properties_.get());
    probe->setGuideRate(&this->guideRate_);
    if (probe->isVFPActive(deferred_logger)) {
        probe->setPrevSurfaceRates(self->wellState(), this->prevWellState());
    }
    self->network_.initializeWell(*probe);
    probe->calculateExplicitQuantities(simulator_, this->groupStateHelper());
    return probe;
}

template<typename TypeTag>
std::vector<typename BlackoilWellModel<TypeTag>::Scalar>
BlackoilWellModel<TypeTag>::
inflowAtBhp_(const WellInterface<TypeTag>& well, const Scalar bhp, DeferredLogger& deferred_logger) const
{
    std::vector<Scalar> r(this->phaseUsage().numActivePhases(), Scalar(0));
    well.computeWellRatesWithBhp(simulator_, bhp, r, deferred_logger);
    return r;
}

template<typename TypeTag>
typename BlackoilWellModel<TypeTag>::Scalar
BlackoilWellModel<TypeTag>::
thpMarginAt_(const WellInterface<TypeTag>& well, const Scalar bhp, const Scalar thp, const Scalar alq,
             std::vector<Scalar>& rates, DeferredLogger& deferred_logger) const
{
    rates = inflowAtBhp_(well, bhp, deferred_logger);
    std::vector<Scalar> rv = rates;
    well.adaptRatesForVFP(rv);
    return WellBhpThpCalculator<Scalar, IndexTraits>(well).thpMargin(
        [&rv](const Scalar) { return rv; }, bhp, this->summaryState(), well.refDensity(), alq, thp);
}

template<typename TypeTag>
typename BlackoilWellModel<TypeTag>::LiftAnswer
BlackoilWellModel<TypeTag>::
liftByScan_(const WellInterface<TypeTag>& well, const Scalar thp, const Scalar alq, DeferredLogger& deferred_logger) const
{
    // The margin bhp - tubing(q(bhp)) over bhp, with the inflow at fixed bhp (unique). Lifts where it is
    // positive somewhere the well produces; the stable crossing is the lower-bhp sign change.
    const int oil = this->phaseUsage().canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
    LiftAnswer a;
    a.thp = thp;
    a.determined = true;
    a.by_scan = true;
    auto produces = [oil](const std::vector<Scalar>& r) { return -r[oil] > Scalar(0); };
    const Scalar b_lo = WellBhpThpCalculator<Scalar, IndexTraits>(well).mostStrictBhpFromBhpLimits(this->summaryState());
    std::vector<Scalar> r_lo;
    const Scalar m_lo = thpMarginAt_(well, b_lo, thp, alq, r_lo, deferred_logger);
    ++a.solves;
    if (!produces(r_lo)) {
        return a;
    }
    if (m_lo > Scalar(0)) {
        a.lifts = a.stable = true;
        a.bhp = b_lo;
        a.flux = r_lo;
        return a;
    }
    Scalar hi = b_lo, step = 10.0e5;
    std::vector<Scalar> r;
    for (int k = 0; k < 8; ++k, step *= 2) {
        hi += step;
        r = inflowAtBhp_(well, hi, deferred_logger);
        ++a.solves;
        if (!produces(r)) {
            break;
        }
    }
    Scalar best_b = b_lo, best_m = m_lo;
    constexpr int npts = 12;
    for (int i = 1; i < npts; ++i) {
        const Scalar b = b_lo + (hi - b_lo) * i / npts;
        const Scalar m = thpMarginAt_(well, b, thp, alq, r, deferred_logger);
        ++a.solves;
        if (produces(r) && m > best_m) {
            best_m = m;
            best_b = b;
        }
    }
    if (best_m <= Scalar(0)) {
        return a;
    }
    Scalar lo = b_lo, up = best_b;
    std::vector<Scalar> r_up = inflowAtBhp_(well, up, deferred_logger);
    for (int k = 0; k < 14 && up - lo > 0.01e5; ++k) {
        const Scalar mid = 0.5 * (lo + up);
        const Scalar m = thpMarginAt_(well, mid, thp, alq, r, deferred_logger);
        ++a.solves;
        if (m > Scalar(0)) { up = mid; r_up = r; } else { lo = mid; }
    }
    a.lifts = a.stable = true;
    a.bhp = up;
    a.flux = r_up;
    return a;
}

template<typename TypeTag>
auto
BlackoilWellModel<TypeTag>::
liftFractions_(const std::string& name) -> LiftFractions
{
    if (name == "explicit") { return LiftFractions::Explicit; }
    if (name == "implicit") { return LiftFractions::Implicit; }
    return LiftFractions::Inherited;
}

template<typename TypeTag>
typename BlackoilWellModel<TypeTag>::LiftAnswer
BlackoilWellModel<TypeTag>::
liftTest_(const WellInterface<TypeTag>& live, const Scalar thp, const Scalar alq, DeferredLogger& deferred_logger,
          const bool probe_stopped, const LiftFractions fractions) const
{
    // The fractions are the well's flag for the test's duration; the well is left as it was.
    struct FlagGuard {
        WellInterface<TypeTag>* w; bool was; bool armed;
        FlagGuard(const WellInterface<TypeTag>& well, const LiftFractions f)
            : w(const_cast<WellInterface<TypeTag>*>(&well)), was(well.vfpExplicitFlag()), armed(f != LiftFractions::Inherited)
        { if (armed) { w->setVfpExplicitFlag(f == LiftFractions::Explicit); } }
        ~FlagGuard() { if (armed) { w->setVfpExplicitFlag(was); } }
    } flag_guard(live, fractions);
    // The well's own equations under thp control, started from its last flowing point. A settled crossing is kept
    // only on the stable branch (positive margin just above its bhp); an unsettled or unstable answer is decided
    // by the scan and, where it lifts, the equations are asked again from the scan's crossing. A stopped well is
    // asked on a probe only where asked for: on Norne the probe's answers revive wells that then die (cliffs
    // 22 against 5, Newton 3217 against 1116). The inflow one bar above the crossing goes with the answer.
    using ThpPoint = typename WellInterface<TypeTag>::ThpPoint;
    // OPM_LIFTTEST_COMPARE=1: a stopped well is asked both ways and the two answers logged.
    static const bool compare = std::getenv("OPM_LIFTTEST_COMPARE") != nullptr;
    static thread_local bool comparing = false;    // the two inner calls must not compare again
    if (compare && !comparing && live.wellIsStopped() && !probe_stopped) {
        comparing = true;
        const auto own = liftTest_(live, thp, alq, deferred_logger, false, fractions);
        auto probe = makeProbe_(live.name(), deferred_logger);
        const auto pr = liftTest_(*probe, thp, alq, deferred_logger, false, fractions);
        comparing = false;
        const int oil = this->phaseUsage().canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
        const Scalar b_cmp = pr.lifts ? pr.bhp : own.lifts ? own.bhp
                           : WellBhpThpCalculator<Scalar, IndexTraits>(live).mostStrictBhpFromBhpLimits(this->summaryState()) + 20.0e5;
        const auto q_own = inflowAtBhp_(live, b_cmp, deferred_logger);
        const auto q_pr = inflowAtBhp_(*probe, b_cmp, deferred_logger);
        auto say = [&](const LiftAnswer& a) {
            return fmt::format("det {} lifts {} stable {} scan {} bhp {:.2f} oil {:.1f}", a.determined ? 1 : 0, a.lifts ? 1 : 0,
                               a.stable ? 1 : 0, a.by_scan ? 1 : 0, a.bhp / 1e5, a.lifts ? -a.flux[oil] * 86400.0 : 0.0);
        };
        OpmLog::debug(fmt::format("LIFTCMP step={} {} thp {:.2f} alq {:.0f} | own: {} | probe: {} | inflow at {:.2f} bar: own {:.1f} probe {:.1f} sm3/d"
                                  " | rho own {:.1f} probe {:.1f} | explicit own {} probe {}",
                                          simulator_.episodeIndex(), live.name(), thp / 1e5, alq * 86400.0, say(own), say(pr),
                                          b_cmp / 1e5, -q_own[oil] * 86400.0, -q_pr[oil] * 86400.0,
                                          live.refDensity(), probe->refDensity(), live.useVfpExplicit() ? 1 : 0,
                                          probe->useVfpExplicit() ? 1 : 0));
        return own;
    }
    WellInterfacePtr probe;
    if (probe_stopped && live.wellIsStopped()) {
        probe = makeProbe_(live.name(), deferred_logger);
        if (fractions != LiftFractions::Inherited) { probe->setVfpExplicitFlag(fractions == LiftFractions::Explicit); }
    }
    const WellInterface<TypeTag>& well = probe ? *probe : live;
    auto finish = [&](LiftAnswer a) {
        if (a.lifts) {
            a.flux_up = inflowAtBhp_(well, a.bhp + 1.0e5, deferred_logger);
            ++a.solves;
        }
        return a;
    };
    const int oil = this->phaseUsage().canonicalToActivePhaseIdx(IndexTraits::oilPhaseIdx);
    std::optional<ThpPoint> start;
    if (const auto fp = this->controller_flowing_point_.find(well.name()); fp != this->controller_flowing_point_.end()) {
        start.emplace();
        start->thp = fp->second[0];
        start->bhp = fp->second[1];
        start->flux.assign(fp->second.begin() + 2, fp->second.end());
        start->lifts = true;
    }
    LiftAnswer a;
    a.thp = thp;
    auto settled = [&](const ThpPoint& pt) {
        a.determined = true;
        a.lifts = pt.lifts;
        a.bhp = pt.bhp;
        a.flux = pt.flux;
        a.rates = pt.rates;
        if (a.lifts) {
            std::vector<Scalar> r;
            const Scalar m = thpMarginAt_(well, pt.bhp + 0.5e5, thp, alq, r, deferred_logger);
            ++a.solves;
            a.stable = -r[oil] > Scalar(0) && m > Scalar(0);
        }
    };
    auto curve = well.sampleThpCurve(simulator_, this->groupStateHelper(), alq, {thp}, a.solves, start ? &*start : nullptr);
    if (!curve.empty() && curve.front().determined) {
        settled(curve.front());
        if (!a.lifts || a.stable) {
            return finish(a);
        }
        ++a.unstable;
    }
    LiftAnswer sc = liftByScan_(well, thp, alq, deferred_logger);
    sc.solves += a.solves;
    sc.unstable = a.unstable;
    if (!sc.lifts) {
        return finish(sc);
    }
    ThpPoint from;
    from.thp = thp;
    from.bhp = sc.bhp;
    from.flux = sc.flux;
    from.lifts = true;
    curve = well.sampleThpCurve(simulator_, this->groupStateHelper(), alq, {thp}, sc.solves, &from);
    if (!curve.empty() && curve.front().determined && curve.front().lifts) {
        LiftAnswer again;
        again.thp = thp;
        again.solves = sc.solves;
        again.unstable = sc.unstable;
        again.by_scan = true;
        again.determined = again.lifts = true;
        again.bhp = curve.front().bhp;
        again.flux = curve.front().flux;
        again.rates = curve.front().rates;
        std::vector<Scalar> r;
        const Scalar m = thpMarginAt_(well, again.bhp + 0.5e5, thp, alq, r, deferred_logger);
        ++again.solves;
        again.stable = -r[oil] > Scalar(0) && m > Scalar(0);
        if (again.stable) {
            return finish(again);
        }
    }
    sc.rates = sc.flux;
    well.adaptRatesForVFP(sc.rates);
    return finish(sc);
}

} // namespace Opm

#endif
