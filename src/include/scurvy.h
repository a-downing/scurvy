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
    inline std::optional<solution_t> solve_half(const problem_t &prob, bool auto_transform_dfp = true)
    {
        const auto _prob = !prob.is_acc() && prob.afp() && auto_transform_dfp ? prob.as_dfp() : prob;

        auto best_dv = 0.0;
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
            best_dv = sols_ca->vp() - sols_ca->prob.v0;
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

        if(sol_half && impl::is_close(sol_half->vp(), sol_half->prob.vf)) {
            log("solve: using sol_half, has exact solution\n");
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

            if(next && sol->vf(true) < next->v0 && !impl::is_close(sol->vf(true), next->v0)) {
                next->v0 = std::abs(sol->vf());
            }

            if(sol->vf(true) > prob.vf && !impl::is_close(sol->vf(true), prob.vf)) {
                if(i_backtrack_last == i) {
                    log("stuck in infinite loop backtracking\n");
                    return std::nullopt;
                }

                i_backtrack_last = i;

                log("solve_path: overshoot: i: %d, %s: sol->prob.v0: %g, sol->vf(): %g -> sol->prob.vf: %g, err: %.17g\n", i, sol->type_name(), sol->prob.v0, sol->vf(), sol->prob.vf, sol->vf() - sol->prob.vf);

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
