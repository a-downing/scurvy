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

    inline std::array<std::optional<solution_t>, 3> ncv_nca(const scurvy::problem_t &prob) {
        auto result = std::array<std::optional<solution_t>, 3>();
        log(solution_type_t::NCV_NCA, "%s\n", __func__);

        // hack, there are precision issues when this case slightly overshoots vf and the other solutions are left with only a tiny time for the deceleration phase
        //auto prob = _prob;
        //prob.J *= 1.0 - 1e-2;

        const auto [V, A, D, J, L, v0, vf] = prob;
        auto x_roots = ncv_nca_x_roots(prob);

        for(size_t i = 0; i < x_roots.size(); i++) {
            auto x = x_roots[i];
            auto vp = v0 + 0.25*J*(x*x);
            auto dist = (v0 + vp)/2 * x;

            if(i > 0) {
                log(solution_type_t::NCV_NCA, "\n");
            }

            log(solution_type_t::NCV_NCA, "    x: %g\n", x);
            log(solution_type_t::NCV_NCA, "    v0: %g\n", v0);
            log(solution_type_t::NCV_NCA, "    vp: %g\n", vp);
            log(solution_type_t::NCV_NCA, "    vf: %g\n", vf);
            log(solution_type_t::NCV_NCA, "    J * 0.5*x: %g\n", J * 0.5*x);
            log(solution_type_t::NCV_NCA, "    dist: %g, L: %g, err: %g\n", dist, L, dist - L);

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
                log(solution_type_t::NCV_NCA, "    bad: negative x\n");
                continue;
            }

            if(!is_close(dist, prob.L)) {
                log(solution_type_t::NCV_NCA, "    bad: distance\n");
                continue;
            }

            if(std::abs(J * 0.5*x) > A) {
                log(solution_type_t::NCV_NCA, "    bad: over acc\n");
                continue;
            }

            if(vp < 0 && prob.afp() || -vp < 0 && !prob.afp()) {
                log(solution_type_t::NCV_NCA, "    bad: negative velocity\n");
                continue;
            }

            if(vp > vf && prob.afp() || -vp < -vf && !prob.afp()) {
                log(solution_type_t::NCV_NCA, "    bad: vf over/undershoot: err: %g\n", vp - prob.vf);
                continue;
            }

            auto T1 = 0.5*x;
            auto T2 = 0.0;
            auto T3 = T1;

            log(solution_type_t::NCV_NCA, "    T1: %g, T2: %g, T3: %g\n", T1, T2, T3);

            result[i].emplace(solution_t { prob, { T1, T2, T3, 0, 0, 0, 0 }, solution_type_t::NCV_NCA });
            log(solution_type_t::NCV_NCA, "    success\n");
        }

        return result;
    }
}