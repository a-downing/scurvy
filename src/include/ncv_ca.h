#pragma once

#include <cmath>
#include <optional>
#include <array>

#include <basics.h>

namespace scurvy::impl {
    inline std::array<double, 2> ncv_ca_x_roots(const scurvy::problem_t &prob) {
        auto [V, A, D, J, L, v0, vf] = prob;
        auto a = 1.0;
        auto b = -A/J + 2*v0/A;
        auto c = -2*L/A;
        return solve_quadratic(a, b, c);
    }

    inline std::optional<solution_t> ncv_ca(const scurvy::problem_t &_prob) {
        log(solution_type_t::NCV_CA, "%s\n", __func__);

        auto prob = _prob;
        prob.J *= 0.99;

        const auto [V, A, D, J, L, v0, vf] = prob;
        auto x_roots = ncv_ca_x_roots(prob);

        auto best_dv = 0.0;
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x : x_roots) {
            auto vp = v0 - A*A/J + A*x;
            auto dist = (prob.v0 + vp)/2 * x;

            log(solution_type_t::NCV_CA, "\n");
            log(solution_type_t::NCV_CA, "    x: %.17g\n", x, x);
            log(solution_type_t::NCV_CA, "    v0: %.17g\n", v0);
            log(solution_type_t::NCV_CA, "    vp: %.17g\n", vp);
            log(solution_type_t::NCV_CA, "    vf: %.17g\n", vf);
            log(solution_type_t::NCV_CA, "    dist: %.17g, L: %.17g, err: %.17g\n", dist, L, dist - L);

            if(approx_lt(x, 0)) {
                log(solution_type_t::NCV_CA, "    bad: negative x\n");
                continue;
            }

            if(!is_close(dist, prob.L)) {
                log(solution_type_t::NCV_CA, "    bad: distance\n");
                continue;
            }

            if(approx_lt(vp, 0.0) && prob.afp() || approx_lt(-vp, 0.0) && !prob.afp()) {
                log(solution_type_t::NCV_CA, "    bad: negative velocity\n");
                continue;
            }

            // if(approx_gt(vp , prob.vf) && prob.afp() || approx_lt(-vp, -prob.vf) && !prob.afp()) {
            //     log(solution_type_t::NCV_CA, "    bad: vf over/undershoot: err: %g\n", vp - prob.vf);
            //     continue;
            // }

            auto T1 = A/J;
            auto T2 = x - 2*T1;
            auto T3 = T1;

            if(approx_lt(T2, 0.0)) {
                log(solution_type_t::NCV_CA, "    bad: T2 < 0 (%g)\n", T2);
                continue;
            }

            log(solution_type_t::NCV_NCA, "    success\n");

            auto dv = vp - prob.v0;

            if(dv > best_dv) {
                log(solution_type_t::NCV_NCA, "    best so far\n");
                best_dv = dv;
                best_sol.emplace(solution_t { prob, { T1, T2, T3, 0, 0, 0, 0 }, solution_type_t::NCV_CA });
            }
        }

        return best_sol;
    }
}