#pragma once

#include <cmath>
#include <optional>
#include <array>

#include <basics.h>

namespace scurvy::impl {
    inline std::array<std::complex<double>, 2> ncv_ca_x_roots(const scurvy::problem_t &prob) {
        auto [V, A, D, J, L, v_0, v_f] = prob;
        auto a = 1.0;
        auto b = -A/J + 2*v_0/A;
        auto c = -2*L/A;
        return solve_quadratic(a+0i, b+0i, c+0i);
    }

    inline std::optional<solution_t> ncv_ca(const scurvy::problem_t &_prob) {
        if(DEBUG) {
            std::println(__func__);
        }

        // hack, there are precision issues when this case slightly overshoots v_f and the other solutions are left with only a tiny time for the deceleration phase
        auto prob = _prob;
        prob.J *= 1.0 - 1e-2;

        auto [V, A, D, J, L, v_0, v_f] = prob;
        auto x_roots = ncv_ca_x_roots(prob);

        for(auto xc : x_roots) {
            auto x = xc.real();

            auto vp = v_0 - A*A/J + A*x;

            if(DEBUG) {
                std::println("x: {}", x);
                std::println("vp: {}", vp);
            }

            if(vp < 0 && prob.afp() || -vp < 0 && !prob.afp()) {
                continue;
            }

            if(vp > prob.vf && prob.afp() || -vp < -prob.vf && !prob.afp()) {
                continue;
            }

            auto T1 = A/J;
            auto T2 = x - 2*T1;
            auto T3 = T1;

            if(DEBUG) {
                std::println("T1: {}, T2: {}, T3: {}", T1, T2, T3);
            }

            if(T2 < 0) {
                continue;
            }

            auto L = (prob.v0 + vp)/2 * x;

            if(DEBUG) {
                std::println("dist: {}", L);
            }

            if(!is_close(L, prob.L)) {
                continue;
            }

            return solution_t { prob, { T1, T2, T3, 0, 0, 0, 0 }, solution_type_t::NCV_CA };
        }

        return std::nullopt;
    }
}