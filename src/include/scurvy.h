#pragma once

#include <optional>
#include <array>

#include <basics.h>
#include <cv.h>
#include <ncv.h>
#include <ncv_nca.h>
#include <ncv_ca.h>

namespace scurvy {
    inline std::optional<solution_t> solve(const problem_t &prob) {
        auto sol = impl::ncv_nca(prob.is_acc() ? prob : prob.as_dfp());
        
        if(sol.has_value()) {
            return sol;
        }

        sol = impl::ncv_ca(prob.is_acc() ? prob : prob.as_dfp());

        if(sol.has_value()) {
            return sol;
        }

        auto p = prob;

        if(prob.dfp_optimal()) {
            p = prob.as_dfp();
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
        std::vector<solution_t> solutions;

        for(int i = 0; i < probs.size(); i++) {
            auto prob = probs[i];
            auto next = i < probs.size() - 1 ? std::optional(&probs[i + 1]) : std::nullopt;

            std::fprintf(stderr, "solve_path: %d, v0: %g, vf: %g\n", i, prob.v0, prob.vf);

            if(next) {
                prob.vf = std::min(prob.V, (*next)->v0) * 0.99; // temporary hack
            }

            auto sol = solve(prob);

            if(!sol.has_value()) {
                std::fprintf(stderr, "solve_path: solve() failed\n");
                return std::nullopt;
            }

            if(next) {
                if(std::abs(sol->vf()) > (*next)->v0) {
                    std::fprintf(stderr, "solve_path: overshot, %g -> %g\n", std::abs(sol->vf()), (*next)->v0);

                    sol = solve(prob.inverse());

                    if(sol) {
                        prob.v0 = sol->vf();
                        // temporary hack
                        i = -1;
                        solutions.clear();
                        continue;
                    }

                    return std::nullopt;
                }

                if(std::abs(sol->vf()) < (*next)->v0) {
                    std::fprintf(stderr, "solve_path: undershot: %g -> %g (%g) \n", (*next)->v0, std::abs(sol->vf()), sol->prob.vf);
                    (*next)->v0 = std::abs(sol->vf());
                }
            } else {
                // TODO: do tolerance stuff the right way
                if(std::abs(sol->vf()) > prob.vf + impl::ABSTOL) {
                    std::fprintf(stderr, "solve_path: end overshot, %g -> %g\n", std::abs(sol->vf()), (*next)->vf);

                    sol = solve(prob.inverse());

                    if(sol) {
                        prob.v0 = sol->vf();
                        // temporary hack
                        i = -1;
                        solutions.clear();
                        continue;
                    }
                }
            }

            solutions.push_back(*sol);
        }

        return solutions;
    }
}