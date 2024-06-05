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
    inline std::optional<solution_t> solve_half(const problem_t &prob)
    {
        const auto _prob = !prob.is_acc() && prob.afp() ? prob.as_dfp() : prob;
        auto best_time = 0.0;
        std::optional<solution_t> best_sol = std::nullopt;

        auto sols_nca = impl::ncv_nca(_prob);

        for(auto sol : sols_nca) {
            if(sol) {
                if(sol->periods.time() > best_time) {
                    best_time = sol->periods.time();
                    best_sol = sol;
                }
            }
        }

        auto sols_ca = impl::ncv_ca(_prob);

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

    inline std::optional<solution_t> solve_full(const problem_t &prob) {
        const auto _prob = prob.dfp_optimal() && prob.afp() ? prob.as_dfp() : prob;

        const std::array sols = {
            impl::ncv_ca_cd(_prob),
            impl::ncv_nca_ncd(_prob),
            impl::ncv_nca_cd(_prob),
            impl::ncv_ca_ncd(_prob),
            impl::cv_ca_cd(_prob),
            impl::cv_ca_ncd(_prob),
            impl::cv_nca_cd(_prob),
            impl::cv_nca_ncd(_prob),
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
        if(impl::is_close(prob.v0, prob.V) && impl::is_close(prob.vf, prob.V)) {
            return std::make_optional(solution_t { prob, { 0, 0, 0, std::abs(prob.L) / prob.V, 0, 0, 0 }, solution_type_t::CV });
        }

        const auto _prob = prob.regularized();

        auto sol = solve_full(_prob);

        if(sol) {
            return sol;
        }

        return solve_half(_prob);
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

        int iterations = 0;

        for(size_t i = 0; i < probs.size(); i++) {
            iterations++;
            auto &prob = probs[i];
            auto next = i < probs.size() - 1 ? &probs[i + 1] : nullptr;

            log("solve_path: iteration: %d\n", i);

            if(next) {
                prob.vf = std::min(next->v0, prob.V);
            }

            auto sol = solve(prob);

            if(!sol.has_value()) {
                log("solve_path: solve() failed\n");
                return std::nullopt;
            }

            if(next && std::abs(sol->vf()) < next->v0 && !impl::is_close(std::abs(sol->vf()), next->v0)) {
                next->v0 = std::abs(sol->vf());
            }

            if(std::abs(sol->vf()) > prob.vf && !impl::is_close(std::abs(sol->vf()), prob.vf)) {
                log("solve_path: overshoot: i: %d, %s: sol->prob.v0: %g, sol->vf(): %g -> prob.vf: %g\n", i, sol->type_name(), sol->prob.v0, sol->vf(), sol->prob.vf, prob.vf);

                if(i == 0) {
                    log("first segment unsolvable and nowhere to backtrack\n");
                    return std::nullopt;
                }

                auto inv = solve(prob.as_dfp().as_inverse());

                if(inv) {
                    log("USING: %s: inv->vf(): %g\n", inv->type_name(), inv->vf());
                    prob.v0 = inv->vf();
                    i -= 2;
                    continue;
                }

                log("solve_path: failed to resolve end overshoot\n");
                return std::nullopt;
            }

            solutions[i] = *sol;
        }

        log("solved %zu segments in %d iterations\n", probs.size(), iterations);

        return solutions;
    }
}