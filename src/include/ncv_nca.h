#pragma once

#include <cmath>
#include <optional>
#include <array>

#include <maths.h>
#include <solve.h>

namespace scurvy::impl {
    inline std::array<std::complex<double>, 3> ncv_nca_x_roots(const scurvy::problem_t &prob) {
        auto [V, A, D, J, L, v_0, v_f] = prob;
        auto a = 1.0;
        auto b = 0.0;
        auto c = 8*v_0/J;

        auto ra = 27*L*L + 32*std::pow(v_0, 3)/J;
        auto sqrt3 = std::sqrt(3.0);

        auto z = std::pow(J, 3)*(sqrt3*std::sqrt(ra+0i)/J + 9*L/J);

        if(DEBUG) {
            std::fprintf(stderr, "a: %g\n", a);
            std::fprintf(stderr, "b: %g\n", b);
            std::fprintf(stderr, "c: %g\n", c);
            std::fprintf(stderr, "ra: %g\n", ra);
            std::fprintf(stderr, "z: %g + %gi\n", z.real(), z.imag());
        }

        if(z == 0.0) {
            return { NAN_CD, NAN_CD, NAN_CD };
        }
        
        auto d = -4.0/9.0*sqrt3*std::sqrt(ra+0i)/J - 4*L/J + 128.0/3.0*std::pow(v_0, 3)/z;

        return solve_cubic(a, b, c, d);
    }

    inline std::optional<solution_t> ncv_nca(const scurvy::problem_t &_prob) {
        if(DEBUG) {
            std::printf("%s\n", __func__);
        }

        // hack, there are precision issues when this case slightly overshoots v_f and the other solutions are left with only a tiny time for the deceleration phase
        auto prob = _prob;
        prob.J *= 1.0 - 1e-2;

        const auto [V, A, D, J, L, v0, vf] = prob;
        auto x_roots = ncv_nca_x_roots(prob);

        for(auto xc : x_roots) {
            auto x = xc.real();
            auto vp = v0 + 0.25*J*(x*x);

            if(DEBUG) {
                std::fprintf(stderr, "x: %g\n", x);
                std::fprintf(stderr, "vp: %g\n", vp);
            }

            if(vp < 0 && prob.afp() || -vp < 0 && !prob.afp()) {
                continue;
            }

            if(vp > vf && prob.afp() || -vp < -vf && !prob.afp()) {
                continue;
            }

            if(J * 0.5*x > A) {
                continue;
            }

            auto T1 = 0.5*x;
            auto T2 = 0.0;
            auto T3 = T1;

            auto l = (v0 + vp)/2 * x;

            if(DEBUG) {
                std::fprintf(stderr, "T1: %g, T2: %g, T3: %g\n", T1, T2, T3);
                std::fprintf(stderr, "l: %g, L: %g, err: %g\n", l, L, l - L);
            }

            if(!is_close(l, prob.L)) {
                continue;
            }

            return solution_t { prob, { T1, T2, T3, 0, 0, 0, 0 }, solution_type_t::NCV_NCA };
        }

        return std::nullopt;;
    }
}