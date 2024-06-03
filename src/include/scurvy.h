#pragma once

#include <optional>
#include <array>

#include <basics.h>
#include <cv.h>
#include <ncv.h>
#include <ncv_nca.h>
#include <ncv_ca.h>

namespace scurvy {
    inline std::optional<solution_t> solve(const problem_t &_prob) {
        auto p = _prob.regularized();
        auto sol = impl::ncv_nca(p.is_acc() ? p : p.as_dfp());
        
        if(sol.has_value()) {
            return sol;
        }

        sol = impl::ncv_ca(p.is_acc() ? p : p.as_dfp());

        if(sol.has_value()) {
            return sol;
        }

        if(p.dfp_optimal()) {
            p = p.as_dfp();
        }

        std::array sols = {
            impl::ncv_ca_cd(p),
            impl::ncv_nca_ncd(p),
            impl::ncv_nca_cd(p),
            impl::ncv_ca_ncd(p),
            impl::cv_ca_cd(p),
            impl::cv_ca_ncd(p),
            impl::cv_nca_cd(p),
            impl::cv_nca_ncd(p),
        };

        auto best_time = std::numeric_limits<double>::max();
        std::optional<solution_t> best_sol = std::nullopt;

        for(auto sol : sols) {
            if(!sol.has_value()) {
                continue;
            }

            if(sol->periods.time() < best_time) {
                best_time = sol->periods.time();
                best_sol = sol;
            }
        }

        if(best_sol.has_value()) {
            return *best_sol;
        }

        std::fprintf(stderr, "no solutions found\n");
        return std::nullopt;
    }

    inline std::optional<std::vector<solution_t>> solve_path(std::vector<problem_t> &probs) {
        std::vector<solution_t> solutions(probs.size());
        int backtrack = 1;
        int iterations = 0;

        for(auto &prob : probs) {
            if(&prob == &probs.front()) {
                prob.vf = prob.V;
            } else if(&prob == &probs.back()) {
                prob.v0 = prob.V;
            } else {
                prob.v0 = prob.V;
                prob.vf = prob.V;
            }
        }

        for(int i = 0; i < probs.size(); i++) {
            iterations++;
            auto &prob = probs[i];
            auto next = i < probs.size() - 1 ? &probs[i + 1] : nullptr;

            prob = prob.regularized();

            if(next) {
                *next = next->regularized();
                prob.vf = std::min(next->v0, prob.V);
            }

            auto sol = solve(prob);

            if(!sol.has_value()) {
                std::fprintf(stderr, "solve_path: solve() failed\n");
                prob.print();
                return std::nullopt;
            }

            if(next) {
                if(std::abs(sol->vf()) > next->v0 && !impl::is_close(std::abs(sol->vf()), next->v0)) {
                    auto err = std::abs(sol->vf()) - next->v0;

                    std::fprintf(stderr, "solve_path: overshoot: i: %d, sol->vf(): %.17g -> next->v0: %.17g\n", i, sol->vf(), next->v0);

                    // not efficient
                    prob.v0 *= 0.99;
                    i = i - backtrack - 1;
                    continue;
                }

                if(std::abs(sol->vf()) < next->v0 && !impl::is_close(std::abs(sol->vf()), next->v0)) {
                    std::fprintf(stderr, "solve_path: undershoot: i: %d, sol->vf(): %g -> next->v0: %g\n", i, sol->vf(), next->v0);
                    next->v0 = std::abs(sol->vf());
                }
            } else {
                if(std::abs(sol->vf()) > prob.vf && !impl::is_close(std::abs(sol->vf()), prob.vf)) {
                    auto err = std::abs(sol->vf()) - prob.vf;

                    std::fprintf(stderr, "solve_path: end overshoot: i: %d, %g -> %g\n", i, std::abs(sol->vf()), prob.vf);

                    // not efficient
                    prob.v0 *= 0.99;
                    i = i - backtrack - 1;
                    continue;
                }
            }

            solutions[i] = *sol;
        }

        std::fprintf(stderr, "solved %d segments in %d iterations\n", probs.size(), iterations);

        return solutions;
    }
}