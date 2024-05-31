#pragma once

#include <optional>

#include <basics.h>

namespace scurvy::impl {
    inline std::optional<solution_t> cv_ca_cd(const problem_t &prob) {
        if(DEBUG) {
            std::printf("%s\n", __func__);
        }

        auto [V, A, D, J, L, v_0, v_f] = prob;

        auto x = (V - v_0)/A + A/J;
        auto x_bar = (V - v_f)/D + D/J;
        auto x_hat = calc_x_hat(V, L, v_0, v_f, x, x_bar);

        auto regions = get_periods(prob, x, x_hat, x_bar, true, true, true, V);

        if(!regions.has_value()) {
            return std::nullopt;
        }

        return solution_t { prob, regions.value(), solution_type_t::CV_CA_CD};
    }

    inline std::optional<solution_t> cv_nca_cd(const problem_t &prob) {
        if(DEBUG) {
            std::printf("%s\n", __func__);
        }

        auto [V, A, D, J, L, v_0, v_f] = prob;

        auto x = 2*sqrt((V - v_0)/J);
        auto x_bar = (V - v_f)/D + D/J;
        auto x_hat = calc_x_hat(V, L, v_0, v_f, x, x_bar);

        auto regions = get_periods(prob, x, x_hat, x_bar, true, false, true, V);

        if(!regions.has_value()) {
            return std::nullopt;
        }
        
        return solution_t { prob, regions.value(), solution_type_t::CV_NCA_CD};
    }

    inline std::optional<solution_t> cv_ca_ncd(const problem_t &prob) {
        if(DEBUG) {
            std::printf("%s\n", __func__);
        }

        auto [V, A, D, J, L, v_0, v_f] = prob;

        auto x = (V - v_0)/A + A/J;
        auto x_bar = 2*sqrt((V - v_f)/J);
        auto x_hat = calc_x_hat(V, L, v_0, v_f, x, x_bar);

        auto regions = get_periods(prob, x, x_hat, x_bar, true, true, false, V);

        if(!regions.has_value()) {
            return std::nullopt;
        }
        
        return solution_t { prob, regions.value(), solution_type_t::CV_CA_NCD};
    }

    inline std::optional<solution_t> cv_nca_ncd(const problem_t &prob) {
        if(DEBUG) {
            std::printf("%s\n", __func__);
        }

        auto [V, A, D, J, L, v_0, v_f] = prob;
        
        auto x = (2*sqrt(V - v_0)) / sqrt(J);
        auto x_bar = (2*sqrt(V - v_f)) / sqrt(J);
        auto x_hat = calc_x_hat(V, L, v_0, v_f, x, x_bar);

        auto regions = get_periods(prob, x, x_hat, x_bar, true, false, false, V);

        if(!regions.has_value()) {
            return std::nullopt;
        }
        
        return solution_t { prob, regions.value(), solution_type_t::CV_NCA_NCD};
    }
}