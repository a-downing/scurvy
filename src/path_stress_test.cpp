#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cstdarg>
#include <random>
#include <string>

#include <scurvy.h>

void fail(const char *format, ...) {
    std::va_list args;
    va_start (args, format);
    std::vfprintf(stderr, format, args);
    va_end (args);
    std::exit(1);
}

int main() {
    constexpr int SEGS = 100;
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis1(0.01, 100.0/60);
    std::uniform_real_distribution<> dis2(0.001, 10.0);
    std::vector<scurvy::problem_t> path;

    for(int i = 0; i < SEGS; i++) {
        auto V = dis1(gen);
        auto A = dis1(gen);
        auto D = dis1(gen);
        auto J = dis1(gen);
        auto L = dis2(gen);

        V = 100.0/60;
        A = 40.0/60;
        J = 20.0/60;
        L = 1.0;
        D = A;

        // v0 and vf of interior segments will be solved by scurvy::solve_path
        auto prob = scurvy::problem_t(V, A, D, J, L, 0.0, 0.0);
        path.push_back(prob);
    }

    // the first v0 and last vf are respected
    path.front().v0 = 0.0;
    path.back().vf = 0.0;

    auto sols = scurvy::solve_path(path);

    if(!sols) {
        fail("no solutions\n");
    }

    auto t_start = 0.0;

    for(auto sol : *sols) {
        constexpr int spaces = 10000;

        for(int i = 0; i <= spaces; i++) {
            auto t = sol.periods.time() / spaces * i;
            std::printf("%g %g\n", t_start + t, std::abs(sol.vt(t)));
        }

        std::printf("\n");
        t_start += sol.periods.time();
    }

    return 0;
}
