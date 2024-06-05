#pragma once

#include <optional>

#include <basics.h>

namespace scurvy::impl {
    inline std::optional<solution_t> cv_ca_cd(const problem_t &prob) {
        log("%s\n", __func__);

        auto [V, A, D, J, L, v0, vf] = prob;

        auto x = (V - v0)/A + A/J;
        auto x_bar = (V - vf)/D + D/J;
        auto x_hat = calc_x_hat(V, L, v0, vf, x, x_bar);

        auto regions = get_periods(prob, x, x_hat, x_bar, V, true, true, true);

        if(!regions.has_value()) {
            return std::nullopt;
        }

        return std::make_optional<solution_t>({ prob, regions.value(), solution_type_t::CV_CA_CD });
    }

    inline std::optional<solution_t> cv_nca_cd(const problem_t &prob) {
        log("%s\n", __func__);

        auto [V, A, D, J, L, v0, vf] = prob;

        auto x = 2*sqrt((V - v0)/J);
        auto x_bar = (V - vf)/D + D/J;
        auto x_hat = calc_x_hat(V, L, v0, vf, x, x_bar);

        auto regions = get_periods(prob, x, x_hat, x_bar, V, true, false, true);

        if(!regions.has_value()) {
            return std::nullopt;
        }
        
        return std::make_optional<solution_t>({ prob, regions.value(), solution_type_t::CV_NCA_CD });
    }

    inline std::optional<solution_t> cv_ca_ncd(const problem_t &prob) {
        log("%s\n", __func__);

        auto [V, A, D, J, L, v0, vf] = prob;

        auto x = (V - v0)/A + A/J;
        auto x_bar = 2*sqrt((V - vf)/J);
        auto x_hat = calc_x_hat(V, L, v0, vf, x, x_bar);

        auto regions = get_periods(prob, x, x_hat, x_bar, V, true, true, false);

        if(!regions.has_value()) {
            return std::nullopt;
        }

        return std::make_optional<solution_t>({ prob, regions.value(), solution_type_t::CV_CA_NCD });
    }

    inline std::optional<solution_t> cv_nca_ncd(const problem_t &prob) {
        log("%s\n", __func__);

        auto [V, A, D, J, L, v0, vf] = prob;
        
        auto x = (2*sqrt(V - v0)) / sqrt(J);
        auto x_bar = (2*sqrt(V - vf)) / sqrt(J);
        auto x_hat = calc_x_hat(V, L, v0, vf, x, x_bar);

        auto regions = get_periods(prob, x, x_hat, x_bar, V, true, false, false);

        if(!regions.has_value()) {
            return std::nullopt;
        }

        return std::make_optional<solution_t>({ prob, regions.value(), solution_type_t::CV_NCA_NCD });
    }
}