#pragma once

#include <limits>
#include <cmath>
#include <optional>

#include <basics.h>
#include <solve.h>

namespace scurvy::impl {
    inline void maybe_update_best(const problem_t &prob, const std::optional<periods_t> &periods, solution_type_t type, std::optional<solution_t> &best_sol, double &best_time) {
        if(periods.has_value() && periods->time() < best_time) {
            best_sol.emplace(solution_t { prob, *periods, type });
            best_time = periods->time();
        }
    }

    inline std::optional<solution_t> ncv_ca_cd(const problem_t &prob) {
        log("%s\n", __func__);

        auto [V, A, D, J, L, v0, vf] = prob;

        auto a = A * (A/D + 1);
        auto b = 1/(J*D) * (A + D) * (A*D - 2*(A*A) + 2*v0*J);
        auto c = -2*L - 1/D * (v0 + vf - (A*A)/J) * (vf - v0 + ((A*A) - (D*D))/J);

        auto roots = solve_quadratic(a, b, c);
        auto best_time = std::numeric_limits<double>::max();
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x : roots) {
            auto vp = v0 - (A*A)/J + A*x;
            auto x_bar = ((D*D) - vf*J + J*vp) / (D*J);
            auto periods = get_periods(prob, x, 0, x_bar, vp, false, true, true);
            maybe_update_best(prob, periods, solution_type_t::NCV_CA_CD, best_sol, best_time);
        }

        return best_sol;
    }

    inline std::optional<solution_t> ncv_nca_cd(const problem_t &prob) {
        log("%s\n", __func__);

        auto [V, A, D, J, L, v0, vf] = prob;
        auto a = J*J / (16*D);
        auto b = 0.25*J;
        auto c = 0.25*(2 * (J*v0) / D + D);
        auto d = 2*v0;
        auto e = -2*L + 1/D * (v0 + vf) * (v0 - vf + D*D/J);

        auto roots = solve_poly(a, b, c, d, e);

        auto best_time = std::numeric_limits<double>::max();
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x : roots) {
            if(std::isnan(x)) {
                continue;
            }

            auto vp = v0 + 0.25*J*x*x;
            auto x_bar = (D*D - vf*J + J*vp) / (D*J);
            auto periods = get_periods(prob, x, 0, x_bar, vp, false, false, true);
            maybe_update_best(prob, periods, solution_type_t::NCV_NCA_CD, best_sol, best_time);
        }

        return best_sol;
    }

    inline std::optional<solution_t> ncv_ca_ncd(const problem_t &prob) {
        log("%s\n", __func__);

        auto [V, A, D, J, L, v0, vf] = prob;
        auto a = J*J/(16*A);
        auto b = 0.25*J;
        auto c = 0.25*(2 * (J*vf)/A + A);
        auto d = 2*vf;
        auto e = -2*L + 1/A * (vf + v0) * (vf - v0 + A*A/J);

        auto roots = solve_poly(a, b, c, d, e);

        auto best_time = std::numeric_limits<double>::max();
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x_bar : roots) {
            if(std::isnan(x_bar)) {
                continue;
            }

            auto vp = vf + 0.25*J*x_bar*x_bar;
            auto x = (A*A + J*vp - J*v0) / (A*J);
            auto periods = get_periods(prob, x, 0, x_bar, vp, false, true, false);
            maybe_update_best(prob, periods, solution_type_t::NCV_CA_NCD, best_sol, best_time);
        }

        return best_sol;
    }

    inline std::optional<solution_t> ncv_nca_ncd(const problem_t &prob) {
        log("%s\n", __func__);

        auto [V, A, D, J, L, v0, vf] = prob;
        auto a = 0.25*(vf - v0)*J;
        auto b = J*L;
        auto c = -std::pow(vf - v0, 2);
        auto d = 8*v0*L;
        auto e = -4*(L*L + 1.0/J * std::pow(v0 + vf, 2) * (vf - v0));

        auto roots = solve_poly(a, b, c, d, e);

        auto best_time = std::numeric_limits<double>::max();
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x : roots) {
            if(std::isnan(x)) {
                continue;
            }

            auto vp = v0 + 0.25*J*(x*x);
            // use complex sqrt or just clamp input
            auto x_bar = (2 * std::sqrt(std::max(0.0, vp - vf))) / std::sqrt(J);
            auto periods = get_periods(prob, x, 0, x_bar, vp, false, false, false);
            maybe_update_best(prob, periods, solution_type_t::NCV_NCA_NCD, best_sol, best_time);
        }

        return best_sol;
    }
}