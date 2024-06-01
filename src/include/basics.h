#pragma once

#include <cstdio>
#include <complex>
#include <utility>

#include <Eigen/Eigen>
#include <Eigen/Dense>

#include <maths.h>

namespace scurvy {
    constexpr bool DEBUG = false;

    enum class solution_type_t {
        NCV_CA,
        NCV_NCA_CD,
        NCV_NCA_NCD,
        CV_CA_CD,
        NCV_CA_NCD,
        NCV_CA_CD,
        NCV_NCA,
        CV_CA_NCD,
        CV_NCA_CD,
        CV_NCA_NCD
    };

    inline const char* solution_type_to_string(const solution_type_t type) {
        switch (type) {
            case solution_type_t::NCV_CA: return "NCV_CA";
            case solution_type_t::NCV_NCA_CD: return "NCV_NCA_CD";
            case solution_type_t::NCV_NCA_NCD: return "NCV_NCA_NCD";
            case solution_type_t::CV_CA_CD: return "CV_CA_CD";
            case solution_type_t::NCV_CA_NCD: return "NCV_CA_NCD";
            case solution_type_t::NCV_CA_CD: return "NCV_CA_CD";
            case solution_type_t::NCV_NCA: return "NCV_NCA";
            case solution_type_t::CV_CA_NCD: return "CV_CA_NCD";
            case solution_type_t::CV_NCA_CD: return "CV_NCA_CD";
            case solution_type_t::CV_NCA_NCD: return "CV_NCA_NCD";
            default: return "Unknown solution_type_t";
        }
    }
    
    struct problem_t {
        double V;
        double A;
        double D;
        double J;
        double L;
        double v0;
        double vf;

        problem_t(double V, double A, double D, double J, double L, double v_0, double v_f): V(V), A(A), D(D), J(J), L(L), v0(v_0), vf(v_f) {

        }

        problem_t() = default;

        bool is_acc() const {
            return vf >= v0;
        }

        problem_t as_dfp() const {
            return { V, D, A, J, -L, -v0, -vf }; // A and D swapped too
        }

        problem_t inverse() const {
            return { V, D, A, J, L, vf, v0 };
        }

        bool afp() const {
            return L >= 0;
        }

        double min_dv_time() const {
            if(v0 <= vf) {
                if (vf - v0 <= A*A / J) {
                    return 2 * sqrt((vf - v0) / J);
                } else {
                    return (vf - v0)/A + A/J;
                }
            } else {
                if (v0 - vf <= D*D / J) {
                    return 2 * sqrt((v0 - vf) / J);
                } else {
                    return (v0 - vf)/D + D/J;
                }
            }
        }

        bool dfp_optimal() const {
            auto Lm = (v0 + vf)/2 * min_dv_time();
            return Lm > L;
        }

        void print() const {
            std::fprintf(stderr, "Problem:\n");
            std::fprintf(stderr, "V = %.17g;\n", V);
            std::fprintf(stderr, "A = %.17g;\n", A);
            std::fprintf(stderr, "D = %.17g;\n", D);
            std::fprintf(stderr, "J = %.17g;\n", J);
            std::fprintf(stderr, "L = %.17g;\n", L);
            std::fprintf(stderr, "v0 = %.17g;\n", v0);
            std::fprintf(stderr, "vf = %.17g;\n", vf);
        }
    };

    struct periods_t {
        double T1, T2, T3, T4, T5, T6, T7;

        double time() const {
            return T1 + T2 + T3 + T4 + T5 + T6 + T7;
        }

        double acc_time() const {
            return T1 + T2 + T3;
        }

        double cv_time() const {
            return T4;
        }

        double dec_time() const {
            return T5 + T6 + T7;
        }

        double distance(const problem_t &prob, double v_p) const {
            return 0.5*(prob.v0 + v_p)*acc_time() + 0.5*(v_p + prob.vf)*dec_time() + v_p*cv_time();
        }

        double vf(const problem_t &prob) const {
            auto a = !impl::near_zero(T2) ? prob.A : 0.5*prob.J*acc_time();
            auto d = !impl::near_zero(T6) ? prob.D : 0.5*prob.J*dec_time();

            auto v0 = prob.v0;
            auto v1 = v0 + 0.5*prob.J * T1*T1;
            auto v2 = v1 + a*T2;
            auto v3 = v2 + a*T3 - 0.5*prob.J * T3*T3;
            auto v4 = v3;
            auto v5 = v4 - 0.5 * prob.J * T5*T5;
            auto v6 = v5 - d*T6;
            auto v7 = v6 - d*T7 + 0.5*prob.J * T7*T7;

            return v7;
        }

        double vt(const problem_t &prob, double t) const {
            auto a = !impl::near_zero(T2) ? prob.A : 0.5*prob.J*acc_time();
            auto d = !impl::near_zero(T6) ? prob.D : 0.5*prob.J*dec_time();

            auto v0 = prob.v0;

            auto done = 0.0;
            if(t <= T1) {
                return v0 + 0.5*prob.J * t*t;
            }
            
            auto v1 = v0 + 0.5*prob.J * T1*T1;

            done += T1;
            if(t <= done+T2) {
                t -= done;
                return v1 + a*t;
            }

            auto v2 = v1 + a*T2;

            done += T2;
            if(t <= done+T3) {
                t -= done;
                return v2 + a*t - 0.5*prob.J * t*t;
            }

            auto v3 = v2 + a*T3 - 0.5*prob.J * T3*T3;

            done += T3;
            if(t <= done+T4) {
                return v3;
            }

            auto v4 = v3;

            done += T4;
            if(t <= done+T5) {
                t -= done;
                return v4 - 0.5 * prob.J * t*t;
            }

            auto v5 = v4 - 0.5 * prob.J * T5*T5;

            done += T5;
            if(t <= done+T6) {
                t -= done;
                return v5 - d*t;
            }

            auto v6 = v5 - d*T6;

            done += T6;
            if(t <= done+T7) {
                t -= done;
                return v6 - d*t + 0.5*prob.J * t*t;
            }

            auto v7 = v6 - d*T7 + 0.5*prob.J * T7*T7;
            return v7;
        }

        void print() const {
            printf("T1: %.17g\n", T1);
            printf("T2: %.17g\n", T2);
            printf("T3: %.17g\n", T3);
            printf("T4: %.17g\n", T4);
            printf("T5: %.17g\n", T5);
            printf("T6: %.17g\n", T6);
            printf("T7: %.17g\n", T7);
        }
    };

    struct solution_t {
        problem_t prob;
        periods_t periods;
        solution_type_t type;

        double vf() const {
            return periods.vf(prob);
        }

        double vt(double t) const {
            return periods.vt(prob, t);
        }

        double vp() const {
            if(impl::near_zero(periods.T2)) {
                return prob.v0 + 0.25 * prob.J * std::pow(periods.acc_time(), 2);
            }

            return prob.v0 - std::pow(prob.A, 2)/prob.J + prob.A*periods.acc_time();
        }

        double distance() const {
            return periods.distance(prob, vp());
        }

        const char *type_name() const {
            return solution_type_to_string(type);
        }

        bool cv_case() const {
            return type == solution_type_t::CV_CA_CD || type == solution_type_t::CV_CA_NCD || type == solution_type_t::CV_NCA_CD || type == solution_type_t::CV_NCA_NCD;
        }
    };
}

namespace scurvy::impl {

    inline double calc_x_hat(const double V, const double L, const double v_0, const double v_f, const double x, const double x_bar) {
        return (2*L - (v_0 + V)*x - (V + v_f)*x_bar) / (2*V);
    }

    inline periods_t calc_periods(const double x, const double x_hat, const double x_bar, const double A, const double D, const double J) {
        auto T1 = A/J, T5 = D/J;

        return {
            A/J,
            x - 2*T1,
            T1,
            x_hat,
            T5,
            x_bar - 2*T5,
            T5,
        };
    }

    inline std::optional<periods_t> get_periods(const problem_t &prob, const double x, const double x_hat, const double x_bar, const bool cv, const bool ca, const bool cd, const double v_p) {
        auto [V, A, D, J, L, v_0, v_f] = prob;

        if(DEBUG) {
            std::fprintf(stderr, "x: %g, x_hat: %g, x_bar: %g\n", x, x_hat, x_bar);
        }

        if(x < 0 || x_hat < 0 || x_bar < 0) {
            return std::nullopt;
        }

        if(v_p > V && L >= 0 || -v_p > V && L < 0) {
            return std::nullopt;
        }

        if(v_p < 0 && L >= 0 || -v_p < 0 && L < 0) {
            return std::nullopt;
        }

        auto a = ca ? A : 0.5*J*x;
        auto d = cd ? D : 0.5*J*x_bar;
        auto periods = calc_periods(x, x_hat, x_bar, a, d, J);

        if(periods.T1 * J > (A * (1+RELTOL)) || periods.T5 * J > (D * (1+RELTOL))) {
            return std::nullopt;
        }

        if(cv && periods.T4 < -ABSTOL || !cv && !near_zero(periods.T4)) {
            return std::nullopt;
        }

        if(ca && periods.T2 < -ABSTOL || !ca && !near_zero(periods.T2)) {
            return std::nullopt;
        }

        if(cd && periods.T6 < -ABSTOL || !cd && !near_zero(periods.T6)) {
            return std::nullopt;
        }

        if(DEBUG) {
            std::fprintf(stderr, "l: %g, L: %g, err: %g\n", periods.distance(prob, v_p), L, periods.distance(prob, v_p) - L);
        }

        if(!is_close(periods.distance(prob, v_p), L, RELTOL_DIST, ABSTOL_DIST)) {
            return std::nullopt;
        }

        return periods;
    }
}