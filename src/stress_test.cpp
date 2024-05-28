#include <cstdint>
#include <cstdlib>
#include <cmath>
#include <random>
#include <string>
#include <print>
#include <format>
#include <unordered_map>
#include <chrono>

#include <scurvy.h>

void display(const scurvy::solution_t &sol) {
    auto cmd = std::format("python ../simulate.py {} {} {} {} {} {} {} {} {}", sol.periods.T1, sol.periods.T2, sol.periods.T3, sol.periods.T4, sol.periods.T5, sol.periods.T6, sol.periods.T7, sol.prob.J, sol.prob.v0);
    system(cmd.c_str());
}

void fail(const scurvy::problem_t &prob, std::string message) {
    prob.print();
    std::println("fail: {}", message);
    std::exit(1);
}

int main() {
    std::random_device rd;
    std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis1(0.01, 100.0/60);
    std::uniform_real_distribution<> dis2(0.001, 100.0);

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

        auto prob = scurvy::problem_t(V, A, D, J, L, v0, vf);

        num_problems++;

        auto sol = scurvy::solve(prob);

        if(!sol.has_value()) {
            fail(prob, "no solution");
        }

        if(!scurvy::impl::is_close(sol->distance(), sol->prob.L, scurvy::impl::RELTOL_DIST, scurvy::impl::ABSTOL_DIST)) {
            fail(sol->prob, std::format("{}: wrong distance: {} vs {}, err: {}", sol->type_name(), sol->distance(), sol->prob.L, sol->distance() - sol->prob.L));
        }

        if(sol->type == scurvy::solution_type_t::NCV_CA || sol->type == scurvy::solution_type_t::NCV_NCA) {
            if(sol->vf() > sol->prob.vf && sol->prob.afp() || -sol->vf() < -sol->prob.vf && !sol->prob.afp()) {
                fail(sol->prob, std::format("{}: wrong final velocity: {} vs {}, err: {}", sol->type_name(), sol->vf(), prob.vf, sol->vf() - prob.vf));
            }
        } else {
            if(!scurvy::impl::is_close(sol->vf(), prob.vf)) {
                fail(sol->prob, std::format("{}: wrong final velocity: {} vs {}, err: {}", sol->type_name(), sol->vf(), prob.vf, sol->vf() - prob.vf));
            }
        }

        if(sol->periods.T2 < -scurvy::impl::ABSTOL || sol->periods.T4 < -scurvy::impl::ABSTOL || sol->periods.T6 < -scurvy::impl::ABSTOL) {
            sol->periods.print();
            fail(sol->prob, std::format("{}: bad time period", sol->type_name()));
        }

        if(sol->cv_case()) {
            if(!scurvy::impl::is_close(sol->vp(), sol->prob.V)) {
                auto err = sol->vp() - sol->prob.V;
                fail(sol->prob, std::format("{}: peak velocity for constant velocity case should be V: {} vs {}, err: {}", sol->type_name(), sol->vp(), sol->prob.V, err));
            }
        } else {
            if(sol->vp() > sol->prob.V) {
                auto err = sol->vp() - sol->prob.V;
                fail(sol->prob, std::format("{}: peak velocity over V: {} vs {}, err: {}", sol->type_name(), sol->vp(), sol->prob.V, err));
            }
        }

        if(num_problems % 1000000 == 0) {
            auto now = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> duration = now - start;
            auto per = duration.count() / num_problems / 1e-9;
            std::println("solved {} random problems ({:.1f} ns/problem)", num_problems, per);

            for(auto [type, num_solutions] : stats) {
                std::println("{}: {:.2f}%", scurvy::solution_type_to_string(type), double(num_solutions)/num_problems*100);
            }
        }

        stats[sol->type] += 1;
        //display(*sol);
    }

    return 0;
}
