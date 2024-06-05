#include <cmath>
#include <array>

#include <Eigen/Eigen>
#include <Eigen/Dense>

#include <maths.h>

namespace scurvy::impl {
    using namespace std::complex_literals;

    std::array<double, 2> solve_quadratic(double a, double b, double c) {
        auto q = -0.5 * (b + std::copysign(std::sqrt(diff_of_products(b, b, 4.0*a, c)), b));
        return { q / a, c / q };
    }

    std::array<std::complex<double>, 2> solve_quadratic(std::complex<double> a, std::complex<double> b, std::complex<double> c) {
        auto cs = Eigen::VectorXcd(3);
        cs << a, b, c;
        cs /= a;

        auto m = Eigen::MatrixXcd({
            {0.0+0i, -cs[2]},
            {1.0+0i, -cs[1]},
        });

        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(m, false);
        
        if (solver.info() == Eigen::Success) {
            auto values = solver.eigenvalues();
            return { values[0], values[1] };
        }
        
        //assert("Eigenvalue computation did not succeed." == nullptr);
        return { NAN_CD, NAN_CD };
    }

    std::array<std::complex<double>, 3> solve_cubic(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d) {
        if(a == 0.0) {
            auto rs = solve_quadratic(b, c, d);
            return { NAN_CD, rs[0], rs[1] };
        }

        auto cs = Eigen::VectorXcd(4);
        cs << a, b, c, d;
        cs /= a;

        auto m = Eigen::MatrixXcd({
            {0.0+0i, 0.0+0i, -cs[3]},
            {1.0+0i, 0.0+0i, -cs[2]},
            {0.0+0i, 1.0+0i, -cs[1]},
        });

        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> solver(m, false);
        
        if (solver.info() == Eigen::Success) {
            auto values = solver.eigenvalues();
            return { values[0], values[1], values[2] };
        }
        
        //assert("Eigenvalue computation did not succeed." == nullptr);
        return { NAN_CD, NAN_CD, NAN_CD };
    }

    std::array<double, 3> solve_cubic(double a, double b, double c, double d) {
        if(a == 0.0) {
            auto rs = solve_quadratic(b, c, d);
            return { NAN_D, rs[0], rs[1] };
        }

        auto cs = Eigen::VectorXd(4);
        cs << a, b, c, d;
        cs /= a;

        auto m = Eigen::MatrixXd({
            {0.0, 0.0, -cs[3]},
            {1.0, 0.0, -cs[2]},
            {0.0, 1.0, -cs[1]},
        });

        Eigen::EigenSolver<Eigen::MatrixXd> solver(m, false);
        
        if (solver.info() == Eigen::Success) {
            auto values = solver.eigenvalues();
            return { values[0].real(), values[1].real(), values[2].real() };
        }
        
        //assert("Eigenvalue computation did not succeed." == nullptr);
        return { NAN_D, NAN_D, NAN_D };
    }

    std::array<double, 4> solve_quartic(double a, double b, double c, double d, double e) {
        auto cs = Eigen::VectorXd(5);
        cs << a, b, c, d, e;
        cs /= a;

        auto m = Eigen::MatrixXd({
            {0, 0, 0, -cs[4]},
            {1, 0, 0, -cs[3]},
            {0, 1, 0, -cs[2]},
            {0, 0, 1, -cs[1]},
        });

        Eigen::EigenSolver<Eigen::MatrixXd> solver(m, false);
        
        if (solver.info() == Eigen::Success) {
            auto values = solver.eigenvalues();
            return { values[0].real(), values[1].real(), values[2].real(), values[3].real() };
        }
        
        //assert("Eigenvalue computation did not succeed." == nullptr);
        return { NAN_D, NAN_D, NAN_D, NAN_D };
    }

    std::array<double, 4> solve_poly(double a, double b, double c, double d, double e) {
        if(a == 0 && b == 0) {
            auto rs = solve_quadratic(c, d, e);
            return { NAN_D, NAN_D, rs[0], rs[1] };
        } if(a == 0) {
            auto rs = solve_cubic(b, c, d, e);
            return { NAN_D, rs[0], rs[1], rs[2] };
        }

        return solve_quartic(a, b, c, d, e);
    }
}