#pragma once

#include <cstdio>
#include <complex>

#include <maths.h>

namespace scurvy {
    constexpr bool DEBUG = false;

    enum class solution_type_t {
        CV,
        NCV_CA,
        NCV_NCA,
        NCV_NCA_CD,
        NCV_NCA_NCD,
        NCV_CA_NCD,
        NCV_CA_CD,
        CV_CA_CD,
        CV_CA_NCD,
        CV_NCA_CD,
        CV_NCA_NCD
    };

    inline const char* solution_type_to_string(const solution_type_t type) {
        switch (type) {
            case solution_type_t::CV: return "CV";
            case solution_type_t::NCV_CA: return "NCV_CA";
            case solution_type_t::NCV_NCA: return "NCV_NCA";
            case solution_type_t::NCV_NCA_CD: return "NCV_NCA_CD";
            case solution_type_t::NCV_NCA_NCD: return "NCV_NCA_NCD";
            case solution_type_t::NCV_CA_NCD: return "NCV_CA_NCD";
            case solution_type_t::NCV_CA_CD: return "NCV_CA_CD";
            case solution_type_t::CV_CA_CD: return "CV_CA_CD";
            case solution_type_t::CV_CA_NCD: return "CV_CA_NCD";
            case solution_type_t::CV_NCA_CD: return "CV_NCA_CD";
            case solution_type_t::CV_NCA_NCD: return "CV_NCA_NCD";
            default: return "Unknown solution_type_t";
        }
    }

    inline void log([[maybe_unused]] solution_type_t type, const char *format, ...) {
        if(!DEBUG) {
            return;
        }

        std::va_list args;
        va_start (args, format);
        std::vfprintf(stderr, format, args);
        va_end (args);
    }

    inline void log(const char *format, ...) {
        if(!DEBUG) {
            return;
        }

        std::va_list args;
        va_start (args, format);
        std::vfprintf(stderr, format, args);
        va_end (args);
    }
    
    struct problem_t {
        double V;
        double A;
        double D;
        double J;
        double L;
        double v0;
        double vf;

        problem_t(double V, double A, double D, double J, double L, double v0, double vf): V(V), A(A), D(D), J(J), L(L), v0(v0), vf(vf) {

        }

        problem_t() = default;

        bool is_acc() const {
            return vf >= v0;
        }

        problem_t as_dfp() const {
            return { V, D, A, J, -L, -v0, -vf };
        }

        problem_t as_inverse() const {
            return { V, A, D, J, -L, -vf, -v0 };
        }

        problem_t regularized() const {
            auto _v0 = impl::near_zero(v0) ? 0.0 : v0;
            auto _vf = impl::near_zero(vf) ? 0.0 : vf;

            if(impl::is_close(_v0, _vf)) {
                _vf = _v0;
            }

            return { V, A, D, J, L, _v0, _vf };
        }

        bool afp() const {
            return L >= 0;
        }

        double min_dv_time() const {
            if(v0 <= vf) {
                if (vf - v0 <= A*A / J) {
                    return 2 * sqrt((vf - v0) / J);
                }

                return (vf - v0)/A + A/J;
            }

            if (v0 - vf <= D*D / J) {
                return 2 * sqrt((v0 - vf) / J);
            }

            return (v0 - vf)/D + D/J;
        }

        bool dfp_optimal() const {
            auto Lm = (v0 + vf)/2 * min_dv_time();
            log("dfp_optimal: %d, Lm: %g, L: %g, diff: %g\n", Lm > L, Lm, L, Lm - L);
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

        bool validate() const {
            return T1 > 0-impl::ABSTOL && T2 > 0-impl::ABSTOL && T3 > 0-impl::ABSTOL && T4 > 0-impl::ABSTOL && T5 > 0-impl::ABSTOL && T6 > 0-impl::ABSTOL && T7 > 0-impl::ABSTOL;
        }

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

        double vp(const problem_t &prob) const {
            if(impl::near_zero(T2)) {
                return prob.v0 + 0.25 * prob.J * std::pow(acc_time(), 2);
            }

            return prob.v0 - std::pow(prob.A, 2)/prob.J + prob.A*acc_time();
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
            std::fprintf(stderr, "T1: %.17g\n", T1);
            std::fprintf(stderr, "T2: %.17g\n", T2);
            std::fprintf(stderr, "T3: %.17g\n", T3);
            std::fprintf(stderr, "T4: %.17g\n", T4);
            std::fprintf(stderr, "T5: %.17g\n", T5);
            std::fprintf(stderr, "T6: %.17g\n", T6);
            std::fprintf(stderr, "T7: %.17g\n", T7);
        }
    };

    struct solution_t {
        problem_t prob;
        periods_t periods;
        solution_type_t type;

        double vf(bool actual = false) const {
            return actual ? (prob.afp() ? periods.vf(prob) : -periods.vf(prob)) : periods.vf(prob);
        }

        double vt(double t) const {
            return periods.vt(prob, t);
        }

        double vp(bool actual = false) const {
            return actual ? (prob.afp() ? periods.vp(prob) : -periods.vp(prob)) : periods.vp(prob);
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
    inline double calc_x_hat(const double V, const double L, const double v0, const double vf, const double x, const double x_bar) {
        return (2*L - (v0 + V)*x - (V + vf)*x_bar) / (2*V);
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

    inline std::optional<periods_t> get_periods(const problem_t &prob,  double x,  double x_hat,  double x_bar, double vp, const bool cv, const bool ca, const bool cd) {
        auto [V, A, D, J, L, v0, vf] = prob;

        auto dist = 0.5*(prob.v0 + vp)*x + 0.5*(vp + prob.vf)*x_bar + vp*x_hat;

        log("\n");
        log("    x: %.17g, x_hat: %.17g, x_bar: %.17g\n", x, x_hat, x_bar);
        log("    v0: %.17g\n", prob.v0);
        log("    vp: %.17g\n", vp);
        log("    vf: %.17g\n", prob.vf);
        log("    dist: %.17g, L: %.17g, err: %.17g\n", dist, L, dist - L);

        if(!is_close(dist, L, RELTOL_DIST, ABSTOL_DIST)) {
            log("    bad: distance\n");
            return std::nullopt;
        }

        if(approx_lt(x, 0.0) || approx_lt(x_hat, 0.0) || approx_lt(x_bar, 0.0)) {
            log("    bad: negative time\n");
            return std::nullopt;
        }

        if(approx_gt(vp, V) && prob.afp() || approx_gt(-vp, V) && !prob.afp()) {
            log("    bad: vp > V\n");
            return std::nullopt;
        }

        if(approx_lt(vp, 0.0) && prob.afp() || approx_lt(-vp, 0.0) && !prob.afp()) {
            log("    bad: negative velocity\n");
            return std::nullopt;
        }

        auto a = ca ? A : 0.5*J*x;
        auto d = cd ? D : 0.5*J*x_bar;
        auto periods = calc_periods(x, x_hat, x_bar, a, d, J);

        // if(approx_gt(periods.vp(prob), V) && prob.afp() || approx_gt(-periods.vp(prob), V) && !prob.afp()) {
        //     log("    bad: vp > V\n");
        //     return std::nullopt;
        // }
        //
        // if(approx_lt(vp, 0.0) && prob.afp() || approx_lt(-vp, 0.0) && !prob.afp()) {
        //     log("    bad: negative velocity\n");
        //     return std::nullopt;
        // }

        if(!is_close(periods.vf(prob), prob.vf)) {
            log("    bad: vf(): %.17g, err: %g\n", periods.vf(prob), periods.vf(prob) - prob.vf);
            return std::nullopt;
        }

        if(approx_gt(periods.T1 * J, A)) {
            log("    bad: over acc\n");
            return std::nullopt;
        }

        if(approx_gt(periods.T5 * J, D)) {
            log("    bad: over dec\n");
            return std::nullopt;
        }

        if(!periods.validate()) {
            log("    bad: negative time period\n");
        }

        if(cv && approx_lt(periods.T4, 0.0) || !cv && !near_zero(periods.T4)) {
            log("    bad: T4 mismatch: %g\n", periods.T4);
            return std::nullopt;
        }

        if(ca && approx_lt(periods.T2, 0.0) || !ca && !near_zero(periods.T2)) {
            log("    bad: T2 mismatch: %g\n", periods.T2);
            return std::nullopt;
        }

        if(cd && approx_lt(periods.T6, 0.0) || !cd && !near_zero(periods.T6)) {
            log("    bad: T6 mismatch: %g\n", periods.T6);
            return std::nullopt;
        }

        log("    success\n");

        return periods;
    }
}