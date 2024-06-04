#pragma once

#include <optional>
#include <array>

#include <basics.h>
#include <cv.h>
#include <ncv.h>
#include <ncv_nca.h>
#include <ncv_ca.h>

namespace scurvy
{
    inline std::optional<solution_t> solve_half(const problem_t &_prob)
    {
        auto p = _prob.regularized();
        auto best_time = 0.0;
        std::optional<solution_t> best_sol = std::nullopt;

        auto sols_nca = impl::ncv_nca(p.is_acc() ? p : p.as_dfp());

        for(auto sol : sols_nca) {
            if(sol) {
                if(sol->periods.time() > best_time) {
                    best_time = sol->periods.time();
                    best_sol = sol;
                }
            }
        }

        auto sols_ca = impl::ncv_ca(p.is_acc() ? p : p.as_dfp());

        for(auto sol : sols_ca) {
            if(sol) {
                if(sol->periods.time() > best_time) {
                    best_time = sol->periods.time();
                    best_sol = sol;
                }
            }
        }

        if(!best_sol) {
            log("solve_half: no solutions found\n");
        }

        return best_sol;
    }

    inline std::optional<solution_t> solve_full(const problem_t &_prob) {
        auto p = _prob.regularized();

        if(p.dfp_optimal() && p.afp()) {
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

        if(!best_sol) {
            log("solve_full: no solutions found\n");
        }

        return best_sol;
    }

    inline std::optional<solution_t> solve(const problem_t &prob) {
        auto sol = solve_full(prob);

        if(sol) {
            return sol;
        }

        return solve_half(prob);
    }

    inline std::optional<std::vector<solution_t>> solve_path(std::vector<problem_t> &probs) {
        std::vector<solution_t> solutions(probs.size());

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

        int backtrack = 1;
        int iterations = 0;
        bool backtracking = false;

        for(int i = 0; i < probs.size(); i++) {
            iterations++;
            auto &prob = probs[i];
            auto next = i < probs.size() - 1 ? &probs[i + 1] : nullptr;

            log("solve_path: iteration: %d, backtracking: %d\n", i, backtracking);

            prob = prob.regularized();

            if(next) {
                *next = next->regularized();
                prob.vf = std::min(next->v0, prob.V);
            }

            auto sol = solve(prob);

            if(!sol.has_value()) {
                log("solve_path: solve() failed\n");
                return std::nullopt;
            }

            if(next) {
                if(std::abs(sol->vf()) > next->v0 && !impl::is_close(std::abs(sol->vf()), next->v0)) {
                    backtracking = true;
                    log("solve_path: overshoot: i: %d, %s: sol->prob.v0: %g, sol->vf(): %.17g(%.17g) -> next->v0: %.17g\n", i, sol->type_name(), sol->prob.v0, sol->vf(), sol->prob.vf, next->v0);

                    auto inv = solve(prob.as_dfp().as_inverse());

                    if(inv) {
                        log("USING: %s: inv->vf(): %g\n", inv->type_name(), inv->vf());
                        prob.v0 = inv->vf();
                        i = i - backtrack - 1;
                        continue;
                    }

                    log("solve_path: failed to resolve overshoot\n");
                    return std::nullopt;
                }

                if(std::abs(sol->vf()) < next->v0 && !impl::is_close(std::abs(sol->vf()), next->v0)) {
                    next->v0 = std::abs(sol->vf());
                }
            } else {
                if(std::abs(sol->vf()) > prob.vf && !impl::is_close(std::abs(sol->vf()), prob.vf)) {
                    backtracking = true;
                    log("solve_path: end overshoot: i: %d, %s: sol->prob.v0: %g, sol->vf(): %.17g(%.17g) -> prob.vf: %.17g\n", i, sol->type_name(), sol->prob.v0, sol->vf(), sol->prob.vf, prob.vf);

                    auto inv = solve_half(prob.as_dfp().as_inverse());

                    if(inv) {
                        log("USING: %s: inv->vf(): %g\n", inv->type_name(), inv->vf());
                        prob.v0 = inv->vf();
                        i = i - backtrack - 1;
                        continue;
                    }

                    log("solve_path: failed to resolve end overshoot\n");
                    return std::nullopt;
                }
            }

            backtracking = false;
            solutions[i] = *sol;
        }

        //log("solved %d segments in %d iterations\n", probs.size(), iterations);

        return solutions;
    }
}