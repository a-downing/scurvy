#pragma once

#include <optional>
#include <array>
#include <cassert>

#include <basics.h>
#include <cv.h>
#include <ncv.h>
#include <ncv_nca.h>
#include <ncv_ca.h>

namespace scurvy
{
    inline std::optional<solution_t> solve_half(const problem_t &prob, bool auto_transform_dfp = true) {
        const auto _prob = !prob.is_acc() && prob.afp() && auto_transform_dfp ? prob.as_dfp() : prob;

        auto best_dv = -1.0;
        std::optional<solution_t> best_sol = std::nullopt;

        auto sols_nca = impl::ncv_nca(_prob);
        auto sols_ca = impl::ncv_ca(_prob);

        // rank the best solution by the highest change in velocity
        if(sols_nca && sols_nca->vp() - sols_nca->prob.v0 > best_dv) {
            assert(impl::is_close(sols_nca->vf(), sols_nca->vp()));
            best_dv = sols_nca->vp() - sols_nca->prob.v0;
            best_sol = sols_nca;
        }

        if(sols_ca && sols_ca->vp() - sols_ca->prob.v0 > best_dv) {
            assert(impl::is_close(sols_ca->vf(), sols_ca->vp()));
            best_sol = sols_ca;
        }

        if(!best_sol) {
            log("solve_half: no solutions found\n");
        }

        return best_sol;
    }

    inline std::optional<solution_t> solve_full(const problem_t &prob, bool auto_transform_dfp = true) {
        const auto _prob = prob.dfp_optimal() && prob.afp() && auto_transform_dfp ? prob.as_dfp() : prob;

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

    inline std::optional<solution_t> solve(const problem_t &prob, bool auto_transform_dfp = true) {
        if(impl::is_close(prob.v0, prob.V) && impl::is_close(prob.vf, prob.V)) {
            return std::make_optional(solution_t { prob, { 0, 0, 0, std::abs(prob.L) / prob.V, 0, 0, 0 }, solution_type_t::CV });
        }

        const auto _prob = prob.regularized();

        auto sol_half = solve_half(_prob, auto_transform_dfp);

        if(sol_half && impl::is_close(sol_half->vf(), sol_half->prob.vf)) {
            solve_full(_prob, auto_transform_dfp);
            log("solve: using sol_half, has exact solution: err: %g\n", sol_half->vf() - sol_half->prob.vf);
            return sol_half;
        }

        auto sol_full = solve_full(_prob, auto_transform_dfp);

        if(sol_full) {
            log("solve: using solve_full\n");
            return sol_full;
        }

        if(!sol_half) {
            log("solve: no solution\n");
        }

        return sol_half;
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
        size_t i_backtrack_last = 0;

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

            if(DEBUG) {
                sol->prob.print();
            }

            if(next && impl::approx_lt(sol->vf(true), next->v0, 1e-8, 1e-8)) {
                log("solve_path: undershoot: i: %d, %s: prob.vf: %.17g, next->v0: %.17g, sol->vf(true): %.17g, err: %g\n", i, sol->type_name(), prob.vf, next->v0, sol->vf(true), sol->vf(true) - next->v0);
                next->v0 = sol->vf(true);
            }

            log("sol->prob.v0: %.17g\n", sol->prob.v0);
            log("sol->vf(): %.17g\n", sol->vf());
            log("sol->vf(true): %.17g\n", sol->vf(true));
            log("sol->vp(): %.17g\n", sol->vp());
            log("sol->vp(true): %.17g\n", sol->vp(true));
            log("sol->prob.vf: %.17g\n", sol->prob.vf);
            log("prob.vf: %.17g\n", prob.vf);

            if(impl::approx_gt(sol->vf(true), prob.vf, 1e-8, 1e-8)) {
                if(i_backtrack_last == i) {
                    log("stuck in infinite loop backtracking\n");
                    return std::nullopt;
                }

                i_backtrack_last = i;

                log("solve_path: overshoot: i: %d, %s: prob.v0: %.17g, prob.vf: %.17g, sol->vf(true): %.17g, err: %g\n", i, sol->type_name(), prob.v0, prob.vf, sol->vf(true), sol->vf(true) - prob.vf);

                if(i == 0) {
                    log("first segment unsolvable and nowhere to backtrack\n");
                    return std::nullopt;
                }

                auto inv = solve(prob.as_dfp().as_inverse(), false);

                if(inv) {
                    log("USING: %s: inv->vf(): %g\n", inv->type_name(), inv->vf());
                    assert(inv->vf() > 0);
                    prob.v0 = inv->vf();
                    i -= 2;
                    continue;
                }

                log("solve_path: failed to resolve end overshoot\n");
                return std::nullopt;
            }

            i_backtrack_last = 0;
            solutions[i] = *sol;
        }

        log("solved %zu segments in %d iterations\n", probs.size(), iterations);

        return solutions;
    }
}
