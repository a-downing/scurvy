#pragma once

#include <array>
#include <complex>

namespace scurvy::impl {
    std::array<double, 2> solve_quadratic(double a, double b, double c);
    std::array<std::complex<double>, 2> solve_quadratic(std::complex<double> a, std::complex<double> b, std::complex<double> c);
    std::array<std::complex<double>, 3> solve_cubic(std::complex<double> a, std::complex<double> b, std::complex<double> c, std::complex<double> d);
    std::array<std::complex<double>, 3> solve_cubic(double a, double b, double c, double d);
    std::array<double, 4> solve_quartic(double a, double b, double c, double d, double e);
    std::array<double, 4> solve_poly(double a, double b, double c, double d, double e);
}