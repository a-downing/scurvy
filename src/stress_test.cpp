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
    std::unordered_map<scurvy::solution_type_t, uint64_t> stats;
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

        // TODO: there is a precision issue when v0 and vf are extremely close, but not equal
        // TODO: still need to figure out the best way to deal with this
        if(scurvy::impl::is_close(v0, vf, 1e-5, 1e-5)) {
            vf = v0;
        }

        auto prob = scurvy::problem_t(V, A, D, J, L, v0, vf);

        num_problems++;

        auto sol = scurvy::solve(prob);

        /*
        scurvy::solution_t sol contains a potentially modified version of
        scurvy::problem_t prob that was passed to scurvy::solve(). It's modified
        in one or two ways. If the solution is a deceleration first solution
        then L, v0, and vf will be negative and A and D are swapped. This is because
        this transformation allows the same algorithm to be used for acceleration first
        solutions and deceleration first solutions. The actual distance, initial velocity,
        and final velocity are the absolute values of L, v0, and vf. The other way
        that prob might be modified is if there is not enough distance to reach
        the final velocity from the initial velocity. In these cases the max jerk
        constraint J is perturbed down by 1%. This is to avoid a potential corner
        case that results in a loss of precision. For this reason it's important
        to use the modified problem in the solution instead of the original problem
        passed to scurvy::solve(). The absolute value of vf will still be the value
        of the original problem even if it can't be reached. To get the actual
        final velocity use scurvy::solution_t::vf(). This will also be negative if
        this is a deceleration first solution. This is pretty low level and not that
        ergonomic to use. It was my intention that the resulting solution_t would be
        transformed into something more convinient to use.
        */

        if(!sol.has_value()) {
            fail(prob, "no solution\n");
        }

        /*
        for(int i = 0; i <= 1000; i++) {
            auto t = sol->periods.time() / 1000 * i;
            std::printf("%g %g\n", t, sol->vt(t));
        }
        return 0;*/

        if(!scurvy::impl::is_close(sol->distance(), sol->prob.L, scurvy::impl::RELTOL_DIST, scurvy::impl::ABSTOL_DIST)) {
            fail(sol->prob, "%s: wrong distance: %g vs %g, err: %g\n", sol->type_name(), sol->distance(), sol->prob.L, sol->distance() - sol->prob.L);
        }

        if(sol->type == scurvy::solution_type_t::NCV_CA || sol->type == scurvy::solution_type_t::NCV_NCA) {
            if(sol->vf() > sol->prob.vf && sol->prob.afp() || -sol->vf() < -sol->prob.vf && !sol->prob.afp()) {
                // not actually a failure, the solution is impossible
                //fail(sol->prob, "%s: wrong final velocity: %g vs %g, err: %g\n", sol->type_name(), sol->vf(), prob.vf, sol->vf() - prob.vf);
            }
        } else {
            if(!scurvy::impl::is_close(sol->vf(), sol->prob.vf)) {
                fail(sol->prob, "%s: wrong final velocity: %g vs %g, err: %g\n", sol->type_name(), sol->vf(), sol->prob.vf, sol->vf() - sol->prob.vf);
            }
        }

        if(sol->periods.T2 < -scurvy::impl::ABSTOL || sol->periods.T4 < -scurvy::impl::ABSTOL || sol->periods.T6 < -scurvy::impl::ABSTOL) {
            sol->periods.print();
            fail(sol->prob, "%s: bad time period\n", sol->type_name());
        }

        if(sol->cv_case()) {
            if(!scurvy::impl::is_close(sol->vp(), sol->prob.V)) {
                auto err = sol->vp() - sol->prob.V;
                fail(sol->prob, "%s: peak velocity for constant velocity case should be V: %g vs %g, err: %g\n", sol->type_name(), sol->vp(), sol->prob.V, err);
            }
        } else {
            if(sol->vp() > sol->prob.V) {
                auto err = sol->vp() - sol->prob.V;
                fail(sol->prob, "%s: peak velocity over V: %g vs %g, err: %g\n", sol->type_name(), sol->vp(), sol->prob.V, err);
            }
        }

        if(num_problems % 1000000 == 0) {
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = now - start;
            auto per = duration.count() / num_problems / 1e-9;
            std::printf("solved %lu random problems (%.1f ns/problem)\n", num_problems, per);

            for(auto [type, num_solutions] : stats) {
                std::printf("%s: %.2f\n", scurvy::solution_type_to_string(type), double(num_solutions)/num_problems*100);
            }
        }

        stats[sol->type] += 1;
    }
}
