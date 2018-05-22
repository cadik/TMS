// Pallas Solver
// Copyright 2015. All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: ryan.latture@gmail.com (Ryan Latture)

#include "gtest/gtest.h"
#include "rapidjson/document.h"
#include "rapidjson/error/en.h"

#include "pallas/basinhopping.h"
#include "pallas/brute.h"
#include "pallas/differential_evolution.h"
#include "pallas/simulated_annealing.h"
#include "pallas/internal/test_functions.h"

namespace pallas {

    class Rosenbrock : public pallas::GradientCostFunction {
    public:
        virtual ~Rosenbrock() {}

        virtual bool Evaluate(const double* parameters,
                              double* cost,
                              double* gradient) const {
            const double x = parameters[0];
            const double y = parameters[1];

            cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
            if (gradient != NULL) {
                gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
                gradient[1] = 200.0 * (y - x * x);
            }
            return true;
        }

        virtual int NumParameters() const { return 2; }
    };

    TEST(RandomNumberGenerator, IntRespectsBounds)
    {
        int min = 0;
        int max = 1;
        internal::RandomNumberGenerator<int> r(min, max);
        std::vector<int> v(20);
        for (std::vector<int>::iterator it = v.begin(); it != v.end(); ++it)
        {
            *it = r();
        }

        for (size_t i = 0; i < v.size(); i++)
        {
            EXPECT_TRUE(0 == v[i] || 1==v[i]);
        }
    }

    TEST(RandomNumberGenerator, DoubleRespectsBounds)
    {
        double min = 2.0;
        double max = 2.1;
        internal::RandomNumberGenerator<double> r(min, max);
        std::vector<double> v(20);
        for (std::vector<double>::iterator it = v.begin(); it != v.end(); ++it)
        {
            *it = r();
        }

        for (size_t i = 0; i < v.size(); i++)
        {
            EXPECT_TRUE(v[i] < max);
            EXPECT_TRUE(v[i] > min);
        }

    }

    TEST(Shuffler, ShufflesVector)
    {
        const unsigned int population_size = 10;
        Eigen::VectorXi indices(population_size);
        indices.setLinSpaced(population_size, 0, population_size);

        Eigen::VectorXi indices_orig = indices;

        internal::Shuffler shuffler(population_size);

        shuffler.Shuffle(indices.data());

        unsigned int num_unchanged = 0;
        for (size_t i = 0; i < indices.size(); i++)
        {
            if(indices[i] == indices_orig[i]) {
                ++num_unchanged;
            }
        }

        EXPECT_TRUE(num_unchanged != population_size) << " Shuffler did not rearrange input data.";

    }

    TEST(Metropolis, TestAcceptReject) {
        internal::Metropolis m;
        double costNew = 1.0;
        double costOld = 2.0;
        double costHigh = 10000000000.0;
        bool accept = m.accept_reject(costNew, costOld);
        bool reject = m.accept_reject(costHigh, costOld);

        EXPECT_TRUE(accept);
        EXPECT_FALSE(reject);
    }

    TEST(Basinhopping, SolvesRosenbrockWithDefaultOptions) {
        const double expected_tolerance = 1e-9;
        double parameters[2] = {-1.2, 0.0};

        pallas::Basinhopping::Options options;
        pallas::Basinhopping::Summary summary;
        pallas::GradientProblem problem(new Rosenbrock());
        pallas::Solve(options, problem, parameters, &summary);

        EXPECT_EQ(TerminationType::CONVERGENCE, summary.termination_type);
        EXPECT_NEAR(1.0, parameters[0], expected_tolerance);
        EXPECT_NEAR(1.0, parameters[1], expected_tolerance);
    }

    TEST(Basinhopping, SavesHistoryOutput) {
        double parameters[2] = {-1.2, 0.0};

        pallas::Basinhopping::Options options;
        options.max_iterations = 2;
        options.history_save_frequency = 1;
        pallas::Basinhopping::Summary summary;
        pallas::GradientProblem problem(new Rosenbrock());
        pallas::Solve(options, problem, parameters, &summary);

        rapidjson::StringBuffer sb;
        HistoryWriter writer(sb);
        rapidjson::Document d;
        dump(summary.history, writer);
        d.Parse(sb.GetString());
        EXPECT_TRUE(!d.Parse(sb.GetString()).HasParseError()) << "Error parsing dumped history data: " << rapidjson::GetParseError_En(d.GetParseError());

        std::vector<std::string> expected_members = {"iteration_number", "stagnant_iterations", "current_solution", "best_cost", "best_solution"};
        for (auto i = 0; i < d.Size(); ++i) {
            for (auto& member: expected_members)
                EXPECT_TRUE(d[i].HasMember(member.c_str())) << "History output missing member: " << member ;
        }

        double history_best_cost = d[d.Size() - 1]["best_cost"].GetDouble();
        double history_best_x = d[d.Size() - 1]["best_solution"].GetArray()[0].GetDouble();
        double history_best_y = d[d.Size() - 1]["best_solution"].GetArray()[1].GetDouble();

        EXPECT_EQ(summary.history.size(), summary.num_iterations);
        EXPECT_DOUBLE_EQ(summary.final_cost, history_best_cost);
        EXPECT_DOUBLE_EQ(parameters[0], history_best_x);
        EXPECT_DOUBLE_EQ(parameters[1], history_best_y);
    }

    TEST(SimulatedAnnealing, SolvesRosenbrockWithoutPolishing) {
        const double expected_tolerance = 0.1;
        double parameters[2] = {-1.2, 0.0};

        Vector upper_bounds(2);
        upper_bounds.setConstant(5);
        Vector lower_bounds = -upper_bounds;

        scoped_ptr<StepFunction> step_function (new BoundedStepFunction(0.1,
                                                                        upper_bounds.data(),
                                                                        lower_bounds.data(),
                                                                        upper_bounds.size()));

        pallas::SimulatedAnnealing::Options options;
        options.set_step_function(step_function);
        options.dwell_iterations = 10000;
        options.max_iterations = 1000;
        options.max_stagnant_iterations = 10000;
        options.cooling_schedule_options.initial_temperature = 1000;
        options.minimum_cost = 0.001;
        pallas::SimulatedAnnealing::Summary summary;
        pallas::GradientProblem problem(new Rosenbrock());
        pallas::Solve(options, problem, parameters, &summary);

        EXPECT_TRUE(summary.termination_type == TerminationType::USER_SUCCESS
                    || summary.termination_type == TerminationType::CONVERGENCE);
        EXPECT_NEAR(1.0, parameters[0], expected_tolerance);
        EXPECT_NEAR(1.0, parameters[1], expected_tolerance);
    }

    TEST(SimulatedAnnealing, SolvesRosenbrockWithPolishing) {
        const double expected_tolerance = 1e-8;
        double parameters[2] = {-1.2, 0.0};

        Vector upper_bounds(2);
        upper_bounds.setConstant(5);
        Vector lower_bounds = -upper_bounds;

        scoped_ptr<StepFunction> step_function(new BoundedStepFunction(0.1,
                                                                       upper_bounds.data(),
                                                                       lower_bounds.data(),
                                                                       upper_bounds.size()));

        pallas::SimulatedAnnealing::Options options;
        options.set_step_function(step_function);
        options.dwell_iterations = 1000;
        options.max_iterations = 1000;
        options.max_stagnant_iterations = 1000;
        options.cooling_schedule_options.initial_temperature = 1000;
        options.minimum_cost = 0.0001;
        options.polish_output = true;
        pallas::SimulatedAnnealing::Summary summary;
        pallas::GradientProblem problem(new Rosenbrock());
        pallas::Solve(options, problem, parameters, &summary);

        EXPECT_TRUE(summary.termination_type == TerminationType::USER_SUCCESS
                    || summary.termination_type == TerminationType::CONVERGENCE);
        EXPECT_NEAR(1.0, parameters[0], expected_tolerance);
        EXPECT_NEAR(1.0, parameters[1], expected_tolerance);
    }

    TEST(SimulatedAnnealing, SavesHistoryOutput) {
        double parameters[2] = {-1.2, 0.0};

        pallas::SimulatedAnnealing::Options options;
        options.dwell_iterations = 1;
        options.max_iterations = 2;
        options.history_save_frequency = 1;
        pallas::SimulatedAnnealing::Summary summary;
        pallas::GradientProblem problem(new Rosenbrock());
        pallas::Solve(options, problem, parameters, &summary);

        rapidjson::StringBuffer sb;
        HistoryWriter writer(sb);
        rapidjson::Document d;
        dump(summary.history, writer);
        EXPECT_FALSE(d.Parse(sb.GetString()).HasParseError()) << "Error parsing dumped history data: " << rapidjson::GetParseError_En(d.GetParseError());

        std::vector<std::string> expected_members = {"iteration_number", "stagnant_iterations", "temperature", "current_solution", "best_cost", "best_solution"};
        for (auto i = 0; i < d.Size(); ++i) {
            for (auto& member: expected_members)
                EXPECT_TRUE(d[i].HasMember(member.c_str())) << "History output missing member: " << member ;
        }

        double history_best_cost = d[d.Size() - 1]["best_cost"].GetDouble();
        double history_best_x = d[d.Size() - 1]["best_solution"].GetArray()[0].GetDouble();
        double history_best_y = d[d.Size() - 1]["best_solution"].GetArray()[1].GetDouble();

        EXPECT_EQ(summary.history.size(), summary.num_iterations);
        EXPECT_DOUBLE_EQ(summary.final_cost, history_best_cost);
        EXPECT_DOUBLE_EQ(parameters[0], history_best_x);
        EXPECT_DOUBLE_EQ(parameters[1], history_best_y);
    }

    TEST(Brute, SolvesRosenbrockWithDefaults) {
        const double expected_tolerance = 1e-8;

        pallas::GradientProblem problem(new Rosenbrock());

        pallas::Brute brute;
        pallas::Brute::Options options;
        pallas::Brute::Summary summary;

        std::vector<Brute::ParameterRange> ranges = {Brute::ParameterRange(-3.0, 3.0, 7),
                                                     Brute::ParameterRange(-3.0, 3.0, 7)};

        Vector parameters(2);
        brute.Solve(options, problem, ranges, parameters.data(), &summary);
        EXPECT_EQ(TerminationType::USER_SUCCESS, summary.termination_type);
        EXPECT_NEAR(1.0, parameters[0], expected_tolerance);
        EXPECT_NEAR(1.0, parameters[1], expected_tolerance);
    }

    TEST(Brute, SavesHistoryOutput) {
        const double expected_tolerance = 1e-8;

        pallas::GradientProblem problem(new Rosenbrock());

        pallas::Brute brute;
        pallas::Brute::Options options;
        options.history_save_frequency = 1;
        pallas::Brute::Summary summary;

        std::vector<Brute::ParameterRange> ranges = {Brute::ParameterRange(-3.0, 3.0, 7),
                                                     Brute::ParameterRange(-3.0, 3.0, 7)};

        Vector parameters(2);
        brute.Solve(options, problem, ranges, parameters.data(), &summary);

        rapidjson::StringBuffer sb;
        HistoryWriter writer(sb);
        rapidjson::Document d;
        dump(summary.history, writer);
        d.Parse(sb.GetString());
        EXPECT_TRUE(!d.Parse(sb.GetString()).HasParseError()) << "Error parsing dumped history data: " << rapidjson::GetParseError_En(d.GetParseError());

        std::vector<std::string> expected_members = {"iteration_number", "current_solution", "best_cost", "best_solution"};
        for (auto i = 0; i < d.Size(); ++i) {
            for (auto& member: expected_members)
                EXPECT_TRUE(d[i].HasMember(member.c_str())) << "History output missing member: " << member;
        }

        double history_best_x = d[d.Size() - 1]["best_solution"].GetArray()[0].GetDouble();
        double history_best_y = d[d.Size() - 1]["best_solution"].GetArray()[1].GetDouble();

        EXPECT_EQ(summary.num_iterations, summary.history.size());
        EXPECT_DOUBLE_EQ(parameters[0], history_best_x);
        EXPECT_DOUBLE_EQ(parameters[1], history_best_y);
    }

    TEST(DifferentialEvolution, SolvesRosenbrockCustomBoundsNoPolish) {
        const double expected_tolerance = 1e-4;
        double parameters[2] = {-1.2, 0.0};

        Vector upper(2);
        upper << 10.0, 10.0;

        Vector lower(2);
        lower << -10.0, -10.0;

        pallas::DifferentialEvolution::Options options;
        options.upper_bounds = upper;
        options.lower_bounds = lower;
        pallas::DifferentialEvolution::Summary summary;
        pallas::GradientProblem problem(new Rosenbrock());
        pallas::Solve(options, problem, parameters, &summary);

        EXPECT_EQ(TerminationType::CONVERGENCE, summary.termination_type);
        EXPECT_NEAR(1.0, parameters[0], expected_tolerance);
        EXPECT_NEAR(1.0, parameters[1], expected_tolerance);
    }

    TEST(DifferentialEvolution, SolvesRosenbrockCustomBoundsWithPolish) {
        const double expected_tolerance = 1e-9;
        double parameters[2] = {-1.2, 0.0};

        Vector upper(2);
        upper << 10.0, 10.0;

        Vector lower(2);
        lower << -10.0, -10.0;

        pallas::DifferentialEvolution::Options options;
        options.upper_bounds = upper;
        options.lower_bounds = lower;
        options.polish_output = true;
        pallas::DifferentialEvolution::Summary summary;
        pallas::GradientProblem problem(new Rosenbrock());
        pallas::Solve(options, problem, parameters, &summary);

        EXPECT_EQ(TerminationType::CONVERGENCE, summary.termination_type);
        EXPECT_NEAR(1.0, parameters[0], expected_tolerance);
        EXPECT_NEAR(1.0, parameters[1], expected_tolerance);
    }

    TEST(DifferentialEvolution, SavesHistoryOutput) {
        double parameters[2] = {-1.2, 0.0};

        Vector upper(2);
        upper << 10.0, 10.0;

        Vector lower(2);
        lower << -10.0, -10.0;

        pallas::DifferentialEvolution::Options options;
        options.upper_bounds = upper;
        options.lower_bounds = lower;
        options.max_iterations = 2;
        options.history_save_frequency = 1;
        pallas::DifferentialEvolution::Summary summary;
        pallas::GradientProblem problem(new Rosenbrock());
        pallas::Solve(options, problem, parameters, &summary);

        rapidjson::StringBuffer sb;
        HistoryWriter writer(sb);
        rapidjson::Document d;
        dump(summary.history, writer);
        d.Parse(sb.GetString());
        EXPECT_TRUE(!d.Parse(sb.GetString()).HasParseError()) << "Error parsing dumped history data: " << rapidjson::GetParseError_En(d.GetParseError());

        std::vector<std::string> expected_members = {"iteration_number", "population", "best_cost", "best_solution"};
        for (auto i = 0; i < d.Size(); ++i) {
            for (auto& member: expected_members)
                EXPECT_TRUE(d[i].HasMember(member.c_str())) << "History output missing member: " << member;
        }

        double history_best_x = d[d.Size() - 1]["best_solution"].GetArray()[0].GetDouble();
        double history_best_y = d[d.Size() - 1]["best_solution"].GetArray()[1].GetDouble();

        EXPECT_EQ(summary.num_iterations + 1, summary.history.size());
        EXPECT_DOUBLE_EQ(parameters[0], history_best_x);
        EXPECT_DOUBLE_EQ(parameters[1], history_best_y);
    }
} // namespace pallas

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}