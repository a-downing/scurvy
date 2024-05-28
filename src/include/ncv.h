#pragma once

#include <limits>
#include <cmath>
#include <optional>

#include <basics.h>
#include <solve.h>

namespace scurvy::impl {
    inline void update_best(const std::optional<periods_t> &periods, std::optional<solution_t> &best_sol, double &best_time, const solution_t &better_sol) {
        if(periods.has_value() && periods->time() < best_time) {
            best_sol = better_sol;
            best_time = periods->time();
        }
    }

    inline std::optional<solution_t> ncv_ca_cd(const problem_t &prob) {
        if(DEBUG) {
            std::println(__func__);
        }

        auto [V, A, D, J, L, v_0, v_f] = prob;

        auto a = A * (A/D + 1);
        auto b = 1/(J*D) * (A + D) * (A*D - 2*(A*A) + 2*v_0*J);
        auto c = -2*L - 1/D * (v_0 + v_f - (A*A)/J) * (v_f - v_0 + ((A*A) - (D*D))/J);

        auto roots = solve_quadratic(a, b, c);
        auto best_time = std::numeric_limits<double>::max();
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x : roots) {
            auto v_p = v_0 - (A*A)/J + A*x;
            auto x_bar = ((D*D) - v_f*J + J*v_p) / (D*J);
            //auto x_bar = -((v_0 + v_p)*x - 2*L)/(v_f + v_p);
            auto periods = get_periods(prob, x, 0, x_bar, false, true, true, v_p);
            update_best(periods, best_sol, best_time, { prob, *periods, solution_type_t::NCV_CA_CD });
        }

        return best_sol;
    }

    inline std::optional<solution_t> ncv_nca_cd(const problem_t &prob) {
        if(DEBUG) {
            std::println(__func__);
        }

        auto [V, A, D, J, L, v_0, v_f] = prob;
        auto a = J*J / (16*D);
        auto b = 0.25*J;
        auto c = 0.25*(2 * (J*v_0) / D + D);
        auto d = 2*v_0;
        auto e = -2*L + 1/D * (v_0 + v_f) * (v_0 - v_f + D*D/J);

        auto roots = solve_poly(a, b, c, d, e);

        auto best_time = std::numeric_limits<double>::max();
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x : roots) {
            if(std::isnan(x)) {
                continue;
            }

            auto v_p = v_0 + 0.25*J*x*x;
            auto x_bar = (D*D - v_f*J + J*v_p) / (D*J);
            //auto x_bar = -((v_0 + v_p)*x - 2*L)/(v_f + v_p);
            auto periods = get_periods(prob, x, 0, x_bar, false, false, true, v_p);
            update_best(periods, best_sol, best_time, { prob, *periods, solution_type_t::NCV_NCA_CD });
        }

        return best_sol;
    }

    inline std::optional<solution_t> ncv_ca_ncd(const problem_t &prob) {
        if(DEBUG) {
            std::println(__func__);
        }

        auto [V, A, D, J, L, v_0, v_f] = prob;
        auto a = J*J/(16*A);
        auto b = 0.25*J;
        auto c = 0.25*(2 * (J*v_f)/A + A);
        auto d = 2*v_f;
        auto e = -2*L + 1/A * (v_f + v_0) * (v_f - v_0 + A*A/J);

        auto roots = solve_poly(a, b, c, d, e);

        auto best_time = std::numeric_limits<double>::max();
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x_bar : roots) {
            if(std::isnan(x_bar)) {
                continue;
            }

            auto v_p = v_f + 0.25*J*x_bar*x_bar;
            auto x = (A*A + J*v_p - J*v_0) / (A*J);
            //auto x = -((v_f + v_p)*x_bar - 2*L)/(v_0 + v_p);
            auto periods = get_periods(prob, x, 0, x_bar, false, true, false, v_p);
            update_best(periods, best_sol, best_time, { prob, *periods, solution_type_t::NCV_CA_NCD });
        }

        return best_sol;
    }

    inline std::optional<solution_t> ncv_nca_ncd(const problem_t &prob) {
        if(DEBUG) {
            std::println(__func__);
        }

        auto [V, A, D, J, L, v_0, v_f] = prob;
        auto a = 0.25*(v_f - v_0)*J;
        auto b = J*L;
        auto c = -std::pow(v_f - v_0, 2);
        auto d = 8*v_0*L;
        auto e = -4*(L*L + 1.0/J * std::pow(v_0 + v_f, 2) * (v_f - v_0));

        auto roots = solve_poly(a, b, c, d, e);

        auto best_time = std::numeric_limits<double>::max();
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x : roots) {
            if(std::isnan(x)) {
                continue;
            }

            auto v_p = v_0 + 0.25*J*(x*x);

            auto x_bar = (2 * std::sqrt(v_p - v_f)) / std::sqrt(J);
            //auto x_bar = -((v_0 + v_p)*x - 2*L)/(v_f + v_p);

            auto periods = get_periods(prob, x, 0, x_bar, false, false, false, v_p);
            update_best(periods, best_sol, best_time, { prob, *periods, solution_type_t::NCV_NCA_NCD });
        }

        return best_sol;
    }
}