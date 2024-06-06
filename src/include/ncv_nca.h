#pragma once

#include <cmath>
#include <optional>
#include <array>

#include <maths.h>
#include <basics.h>
#include <solve.h>

namespace scurvy::impl {
    inline std::array<double, 3> ncv_nca_x_roots(const scurvy::problem_t &prob) {
        auto [V, A, D, J, L, v0, vf] = prob;
        auto a = 1.0;
        auto b = 0.0;
        auto c = 8.0*v0/J;
        auto d = -8.0*L/J;

        return solve_cubic(a, b, c, d);
    }

    inline std::optional<solution_t> ncv_nca(const scurvy::problem_t &_prob) {
        log(solution_type_t::NCV_NCA, "%s\n", __func__);

        auto prob = _prob;
        prob.J *= 0.99;

        const auto [V, A, D, J, L, v0, vf] = prob;
        auto x_roots = ncv_nca_x_roots(prob);

        auto best_dv = 0.0;
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto x : x_roots) {
            auto vp = v0 + 0.25*J*(x*x);
            auto dist = (v0 + vp)/2 * x;

            log(solution_type_t::NCV_NCA, "\n");
            log(solution_type_t::NCV_NCA, "    x: %.17g\n", x);
            log(solution_type_t::NCV_NCA, "    v0: %.17g\n", v0);
            log(solution_type_t::NCV_NCA, "    vp: %.17g\n", vp);
            log(solution_type_t::NCV_NCA, "    vf: %.17g\n", vf);
            log(solution_type_t::NCV_NCA, "    J * 0.5*x: %.17g\n", J * 0.5*x);
            log(solution_type_t::NCV_NCA, "    dist: %.17g, L: %.17g, err: %.17g\n", dist, L, dist - L);

            if(approx_lt(x, 0.0)) {
                log(solution_type_t::NCV_NCA, "    bad: negative x\n");
                continue;
            }

            if(!is_close(dist, prob.L)) {
                log(solution_type_t::NCV_NCA, "    bad: distance\n");
                continue;
            }

            if(approx_gt(J * 0.5*std::abs(x), A)) {
                log(solution_type_t::NCV_NCA, "    bad: over acc\n");
                continue;
            }

            if(approx_lt(vp, 0.0) && prob.afp() || approx_lt(-vp, 0.0) && !prob.afp()) {
                log(solution_type_t::NCV_NCA, "    bad: negative velocity\n");
                continue;
            }

            // if(approx_gt(vp , vf) && prob.afp() || approx_lt(-vp, -vf) && !prob.afp()) {
            //     log(solution_type_t::NCV_NCA, "    bad: vf over/undershoot: err: %g\n", vp - prob.vf);
            //     continue;
            // }

            auto T1 = 0.5*x;
            auto T2 = 0.0;
            auto T3 = T1;

            auto dv = vp - prob.v0;

            log(solution_type_t::NCV_NCA, "    success\n");

            if(dv > best_dv) {
                log(solution_type_t::NCV_NCA, "    best so far\n");
                best_dv = dv;
                best_sol.emplace(solution_t { prob, { T1, T2, T3, 0, 0, 0, 0 }, solution_type_t::NCV_NCA });
            }
        }

        return best_sol;
    }
}