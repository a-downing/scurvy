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

    inline std::array<std::optional<solution_t>, 2> ncv_ca(const scurvy::problem_t &prob) {
        auto result = std::array<std::optional<solution_t>, 2>();
        log(solution_type_t::NCV_CA, "%s\n", __func__);

        // hack, there are precision issues when this case slightly overshoots vf and the other solutions are left with only a tiny time for the deceleration phase
        //auto prob = _prob;
        //prob.J *= 1.0 - 1e-2;

        const auto [V, A, D, J, L, v0, vf] = prob;
        auto x_roots = ncv_ca_x_roots(prob);

        for(size_t i = 0; i < x_roots.size(); i++) {
            auto x = x_roots[i];
            auto vp = v0 - A*A/J + A*x;
            auto dist = (prob.v0 + vp)/2 * x;

            if(i > 0) {
                log(solution_type_t::NCV_CA, "\n");
            }

            log(solution_type_t::NCV_CA, "    x: %g\n", x, x);
            log(solution_type_t::NCV_CA, "    v0: %g\n", v0);
            log(solution_type_t::NCV_CA, "    vp: %g\n", vp);
            log(solution_type_t::NCV_CA, "    vf: %g\n", vf);
            log(solution_type_t::NCV_CA, "    dist: %g, L: %g, err: %g\n", dist, L, dist - L);

            if(near_zero(x)) {
                x = 0;
            }

            if(near_zero(vp)) {
                vp = 0;
            }

            if(is_close(vp, vf)) {
                vp = vf;
            }

            if(x < 0) {
                log(solution_type_t::NCV_CA, "    bad: negative x\n");
                continue;
            }

            if(!is_close(dist, prob.L)) {
                log(solution_type_t::NCV_CA, "    bad: distance\n");
                continue;
            }

            if(vp < 0 && prob.afp() || -vp < 0 && !prob.afp()) {
                log(solution_type_t::NCV_CA, "    bad: negative velocity\n");
                continue;
            }

            if(vp > prob.vf && prob.afp() || -vp < -prob.vf && !prob.afp()) {
                log(solution_type_t::NCV_CA, "    bad: vf over/undershoot: err: %g\n", vp - prob.vf);
                continue;
            }

            auto T1 = A/J;
            auto T2 = x - 2*T1;
            auto T3 = T1;

            if(T2 < 0) {
                log(solution_type_t::NCV_CA, "    bad: T2 < 0 (%g)\n", T2);
                continue;
            }

            result[i].emplace(solution_t { prob, { T1, T2, T3, 0, 0, 0, 0 }, solution_type_t::NCV_CA });
            log(solution_type_t::NCV_NCA, "    success\n");
        }

        return result;
    }
}