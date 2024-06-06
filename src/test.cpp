#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cstdarg>

#include <random>
#include <string>
#include <unordered_map>
#include <chrono>

#include <scurvy.h>

void fail([[maybe_unused]] const scurvy::problem_t &prob, const char *format, ...) {
    prob.print();
    std::va_list args;
    va_start (args, format);
    std::vfprintf(stderr, format, args);
    va_end (args);
    std::exit(1);
}

int simulate(const scurvy::solution_t &sol) {
    char str[1024];
    std::snprintf(str, sizeof(str), "python ../simulate.py %g %g %g %g %g %g %g %g %g", sol.periods.T1, sol.periods.T2, sol.periods.T3, sol.periods.T4, sol.periods.T5, sol.periods.T6, sol.periods.T7, sol.prob.J, sol.prob.v0);
    return system(str);
}

void print_scurve(const scurvy::solution_t &sol, int res = 1000) {
    for(int i = 0; i <= res; i++) {
        auto t = sol.periods.time() / res * i;
        std::printf("%g %g\n", t, std::abs(sol.vt(t)));
    }

    std::printf("\n");
}

int main() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis1(0.01, 100.0/60);
    std::uniform_real_distribution<> dis2(0.001, 10.0);

    uint64_t num_problems = 0;
    auto start = std::chrono::high_resolution_clock::now();

    for(;;) {
        auto V = dis1(gen);
        auto A = dis1(gen);
        auto D = dis1(gen);
        auto J = dis1(gen);
        auto v0 = dis1(gen);
        auto vf = dis1(gen);
        auto L = dis2(gen);

        if(v0 > V || vf > V) {
            continue;
        }

        auto prob = scurvy::problem_t(V, A, D, J, L, v0, vf);
        num_problems++;

        auto sol_nvc_nca = scurvy::impl::ncv_nca(prob);
        auto sol_nvc_ca = scurvy::impl::ncv_ca(prob);

        if(num_problems % 1000000 == 0) {
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = now - start;
            auto per = duration.count() / num_problems / 1e-9;
            std::printf("solved %lu random problems (%.1f ns/problem)\n", num_problems, per);
        }
    }
}
