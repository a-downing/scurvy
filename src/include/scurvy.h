#pragma once

#include <optional>
#include <expected>
#include <array>

#include <basics.h>
#include <cv.h>
#include <ncv.h>
#include <ncv_nca.h>
#include <ncv_ca.h>

namespace scurvy {
    inline std::expected<solution_t, const char *> solve(const problem_t &prob) {
        auto sol = impl::ncv_nca(prob.is_acc() ? prob : prob.as_dfp());
        
        if(sol.has_value()) {
            return *sol;
        }

        sol = impl::ncv_ca(prob.is_acc() ? prob : prob.as_dfp());

        if(sol.has_value()) {
            return *sol;
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

        return std::unexpected("no solutions found");
    }
}