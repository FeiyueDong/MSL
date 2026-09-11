#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

#include "ode.hpp"
#include "test_utils.hpp"

using namespace msl;

int test_ode() {
    auto exponential = [](double, const ode::state &y) {
        return ode::state{y[0]};
    };
    std::vector<double> y0{1.0};

    auto euler_result = ode::euler(exponential, y0, 0.0, 1.0, 1e-3);
    EXPECT_TRUE(euler_result.success());
    EXPECT_EQ(euler_result.status, ode::ode_status::success);
    EXPECT_NEAR(euler_result.final_time(), 1.0, 1e-12);
    EXPECT_NEAR(euler_result.final_state()[0], std::exp(1.0), 2e-3);

    auto heun_result = ode::heun(exponential, y0, 0.0, 1.0, 1e-2);
    EXPECT_TRUE(heun_result);
    EXPECT_NEAR(heun_result.final_state()[0], std::exp(1.0), 5e-5);

    auto rk4_result = ode::rk4(exponential, y0, 0.0, 1.0, 0.1);
    EXPECT_TRUE(rk4_result.success());
    EXPECT_NEAR(rk4_result.final_state()[0], std::exp(1.0), 3e-6);
    EXPECT_NEAR(rk4_result.value()[0], std::exp(1.0), 3e-6);

    auto rk4_final = ode::rk4_final(exponential, y0, 0.0, 1.0, 0.1);
    EXPECT_NEAR(rk4_final[0], std::exp(1.0), 3e-6);

    std::vector<double> t_eval{0.0, 0.25, 0.5, 1.0};
    auto euler_eval = ode::euler_eval(exponential, y0, t_eval, 0.1);
    EXPECT_TRUE(euler_eval.success());
    EXPECT_EQ(euler_eval.t.size(), t_eval.size());
    EXPECT_EQ(euler_eval.y.size(), t_eval.size());
    EXPECT_NEAR(euler_eval.t[1], 0.25, 1e-12);
    EXPECT_NEAR(euler_eval.y.back()[0], std::exp(1.0), 2e-1);

    auto heun_eval = ode::heun_eval(exponential, y0, t_eval, 0.1);
    EXPECT_TRUE(heun_eval.success());
    EXPECT_NEAR(heun_eval.t[2], 0.5, 1e-12);
    EXPECT_NEAR(heun_eval.y.back()[0], std::exp(1.0), 5e-3);

    auto rk4_eval = ode::rk4_eval(exponential, y0, t_eval, 0.1);
    EXPECT_TRUE(rk4_eval.success());
    EXPECT_NEAR(rk4_eval.final_time(), 1.0, 1e-12);
    EXPECT_NEAR(rk4_eval.final_state()[0], std::exp(1.0), 3e-6);

    ode::ode_options adaptive_options;
    adaptive_options.rtol = 1e-8;
    adaptive_options.atol = 1e-10;
    auto ode45_result = ode::ode45(exponential, y0, 0.0, 1.0, adaptive_options);
    EXPECT_TRUE(ode45_result.success());
    EXPECT_NEAR(ode45_result.final_time(), 1.0, 1e-12);
    EXPECT_NEAR(ode45_result.final_state()[0], std::exp(1.0), 1e-7);
    EXPECT_TRUE(ode45_result.steps > 0);

    auto ode45_final =
        ode::ode45_final(exponential, y0, 0.0, 1.0, adaptive_options);
    EXPECT_NEAR(ode45_final[0], std::exp(1.0), 1e-7);

    auto ode45_eval =
        ode::ode45_eval(exponential, y0, t_eval, adaptive_options);
    EXPECT_TRUE(ode45_eval.success());
    EXPECT_EQ(ode45_eval.t.size(), t_eval.size());
    EXPECT_NEAR(ode45_eval.t[1], 0.25, 1e-12);
    EXPECT_NEAR(ode45_eval.final_state()[0], std::exp(1.0), 1e-7);

    auto oscillator = [](double, const ode::state &y) {
        return ode::state{y[1], -y[0]};
    };
    std::vector<double> osc0{1.0, 0.0};
    auto osc_result = ode::rk4(oscillator, osc0, 0.0, std::acos(-1.0), 0.01);
    EXPECT_TRUE(osc_result.success());
    EXPECT_NEAR(osc_result.final_state()[0], -1.0, 1e-8);
    EXPECT_NEAR(osc_result.final_state()[1], 0.0, 1e-8);

    ode::ode_options limited;
    limited.max_steps = 2;
    auto limited_result = ode::rk4(exponential, y0, 0.0, 1.0, 0.1, limited);
    EXPECT_EQ(limited_result.status, ode::ode_status::max_steps);
    bool value_threw = false;
    try {
        (void)limited_result.value();
    } catch (const std::runtime_error &) {
        value_threw = true;
    }
    EXPECT_TRUE(value_threw);

    auto non_finite = ode::euler(
        [](double, const ode::state &) {
            return ode::state{std::numeric_limits<double>::infinity()};
        },
        y0,
        0.0,
        1.0,
        0.1);
    EXPECT_EQ(non_finite.status, ode::ode_status::non_finite_value);

    bool threw = false;
    try {
        (void)ode::rk4(exponential, y0, 1.0, 0.0, 0.1);
    } catch (const std::invalid_argument &) {
        threw = true;
    }
    EXPECT_TRUE(threw);

    bool eval_threw = false;
    try {
        std::vector<double> bad_t_eval{0.0, 0.5, 0.5};
        (void)ode::rk4_eval(exponential, y0, bad_t_eval, 0.1);
    } catch (const std::invalid_argument &) {
        eval_threw = true;
    }
    EXPECT_TRUE(eval_threw);

    auto compare_dir = msl::test::result_dir("ode");
    std::ofstream exp_file(compare_dir / "ode_exp_eval.txt");
    exp_file << std::setprecision(17);
    for (size_t i = 0; i < t_eval.size(); ++i) {
        exp_file << t_eval[i] << " " << euler_eval.y[i][0] << " "
                 << heun_eval.y[i][0] << " " << rk4_eval.y[i][0] << " "
                 << ode45_eval.y[i][0] << "\n";
    }

    std::vector<double> osc_eval_t{0.0, std::acos(-1.0) / 2.0, std::acos(-1.0)};
    auto osc_eval = ode::rk4_eval(oscillator, osc0, osc_eval_t, 0.01);
    std::ofstream osc_file(compare_dir / "ode_oscillator_rk4.txt");
    osc_file << std::setprecision(17);
    for (size_t i = 0; i < osc_eval.t.size(); ++i) {
        osc_file << osc_eval.t[i] << " " << osc_eval.y[i][0] << " "
                 << osc_eval.y[i][1] << "\n";
    }

    return msl::test::summary();
}

int main() { return test_ode(); }
