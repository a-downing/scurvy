#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <cstdarg>
#include <random>
#include <chrono>

#include <scurvy.h>

void fail(const std::vector<scurvy::problem_t> &path, const char *format, ...) {
    // for(size_t i = 0; i < path.size(); i++) {
    //     auto prob = path[i];
    //
    //     std::fprintf(stderr, "path[%zu] = scurvy::problem_t {\n", i);
    //     std::fprintf(stderr, "    %.17g,\n", prob.V);
    //     std::fprintf(stderr, "    %.17g,\n", prob.J);
    //     std::fprintf(stderr, "    %.17g,\n", prob.A);
    //     std::fprintf(stderr, "    %.17g,\n", prob.D);
    //     std::fprintf(stderr, "    %.17g,\n", prob.L);
    //     std::fprintf(stderr, "    %.17g,\n", prob.v0);
    //     std::fprintf(stderr, "    %.17g\n", prob.vf);
    //     std::fprintf(stderr, "};\n");
    // }

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
    std::uniform_real_distribution<> dis1(0.001, 100.0/60);
    std::uniform_real_distribution<> dis2(0.000001, 10.0);
    std::vector<scurvy::problem_t> path;

    for(size_t i = 0;; i++) {
        path.clear();
        
        for(int i = 0; i < SEGS; i++) {
            auto V = dis1(gen);
            auto A = dis1(gen);
            auto D = dis1(gen);
            auto J = dis1(gen);
            auto L = dis2(gen);

            //V = 100.0/60;
            // A = 40.0/60;
            // D = 40.0/60;
            // J = 20.0/60;
            //L = 1.0;

            // v0 and vf of interior segments will be solved by scurvy::solve_path
            auto prob = scurvy::problem_t(V, A, D, J, L, 0.0, 0.0);
            path.push_back(prob);
        }

        // the first v0 and last vf are respected
        path.front().v0 = 0.0;
        path.back().vf = 0.0;

        // save the unmodified path to reproduce a test case on failure
        auto path_copy = path;

        auto start = std::chrono::high_resolution_clock::now();
        auto sols = scurvy::solve_path(path);
        std::chrono::duration<double> duration = std::chrono::high_resolution_clock::now() - start;

        std::fprintf(stderr, "%zu: solved %zu problems in %gs\n", i, path.size(), duration.count());

        if(!sols) {
            fail(path_copy, "no solutions\n");
        }

        continue;

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

        //return 0;
    }
}
