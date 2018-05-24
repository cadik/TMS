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

#include "Eigen/Core"
#include "glog/logging.h"

#include "pallas/basinhopping.h"
#include "pallas/internal/solver_utils.h"
#include "pallas/internal/stringprintf.h"
#include "pallas/internal/wall_time.h"

namespace pallas {

    using std::string;

    using ceres::LineSearchDirectionTypeToString;
    using ceres::TerminationTypeToString;

    using pallas::internal::StringAppendF;
    using pallas::internal::StringPrintf;
    using pallas::internal::WallTimeInSeconds;

    namespace {

        bool Evaluate(const GradientProblem &problem,
                      Vector &x,
                      internal::State *state,
                      std::string *message) {
            if (!problem.Evaluate(x.data(),
                                  &(state->cost),
                                  state->gradient.data())) {
                *message = "Gradient evaluation failed.";
                return false;
            }
            Vector negative_gradient = -state->gradient;
            Vector projected_gradient_step(x.size());
            if (!problem.Plus(x.data(),
                              negative_gradient.data(),
                              projected_gradient_step.data())) {
                *message = "projected_gradient_step = Plus(x, -gradient) failed.";
                return false;
            }

            state->gradient_squared_norm = (x - projected_gradient_step).squaredNorm();
            state->gradient_max_norm =
                    (x - projected_gradient_step).lpNorm<Eigen::Infinity>();
            return true;
        }

    }  // namespace

    Basinhopping::Summary::Summary()
            : termination_type(TerminationType::FAILURE),
              message("pallas::Basinhopping was not called."),
              initial_cost(-1.0),
              final_cost(-1.0),
              num_parameters(0),
              num_iterations(0),
              total_time_in_seconds(0.0),
              local_minimization_time_in_seconds(0.0),
              step_time_in_seconds(0.0),
              cost_evaluation_time_in_seconds(0.0) {

    };

    std::string Basinhopping::Summary::BriefReport() const {
        return StringPrintf(
                "Pallas basinhopping report: "
                "iterations: %d, "
                "initial cost: %e, "
                "final cost: %e, "
                "termination: %s\n",
                num_iterations,
                initial_cost,
                final_cost,
                TerminationTypeToString(termination_type));
    };

    string Basinhopping::Summary::FullReport() const {

        string report = string("\nSolver Summary\n\n");

        StringAppendF(&report, "Parameters          %25d\n", num_parameters);

        string line_search_direction_string = LineSearchDirectionTypeToString(line_search_direction_type);

        StringAppendF(&report, "Line search direction     %19s\n",
                      line_search_direction_string.c_str());

        StringAppendF(&report, "\nCost:\n");
        StringAppendF(&report, "  Initial        %28e\n", initial_cost);
        if (termination_type != TerminationType::FAILURE &&
            termination_type != TerminationType::USER_FAILURE) {
            StringAppendF(&report, "  Final          %28e\n", final_cost);
            StringAppendF(&report, "  Change         %28e\n",
                          initial_cost - final_cost);
        }

        StringAppendF(&report, "\nMinimizer iterations         %16d\n",
                      num_iterations);

        StringAppendF(&report, "\nTime (in seconds):\n");

        StringAppendF(&report, "  Cost evaluation     %23.4f\n",
                      cost_evaluation_time_in_seconds);

        StringAppendF(&report, "  Local minimization   %22.4f\n",
                      local_minimization_time_in_seconds);

        StringAppendF(&report, "  Step function   %27.4f\n",
                      step_time_in_seconds);

        StringAppendF(&report, "  Total               %23.4f\n\n",
                      total_time_in_seconds);

        StringAppendF(&report, "Termination: %2s (%s)\n",
                      TerminationTypeToString(termination_type), message.c_str());
        return report;
    };

    void Basinhopping::Solve(const Basinhopping::Options &options,
                             const GradientProblem &problem,
                             double *parameters,
                             Basinhopping::Summary *global_summary) {

        double start_time = WallTimeInSeconds();
        double t1;

        global_summary->line_search_direction_type =
                options.local_minimizer_options.line_search_direction_type;

        bool is_not_silent = !options.is_silent;
        bool new_global_min = false;
        bool accept;
        num_iterations_ = 0;
        num_stagnant_iterations_ = 0;

        const unsigned int num_parameters = static_cast<unsigned int>(problem.NumParameters());
        VectorRef x(parameters, num_parameters);

        global_summary->num_parameters = num_parameters;

        current_state_ = internal::State(num_parameters);
        current_state_.x = x;

        GradientLocalMinimizer::Summary local_summary;

        // evaluate problem with initial parameters
        if (!Evaluate(problem, current_state_.x, &current_state_, &global_summary->message)) {
            global_summary->termination_type = TerminationType::FAILURE;
            global_summary->message = "Initial cost and jacobian evaluation failed. "
                                              "More details: " + global_summary->message;
            LOG_IF(WARNING, is_not_silent) << "Terminating: " << global_summary->message;
            prepare_final_summary_(global_summary, local_summary);
            return;
        }

        global_summary->initial_cost = current_state_.cost;

        t1 = WallTimeInSeconds();
        // minimize problem with initial parameters
        GradientLocalMinimizer local_minimizer;
        local_minimizer.Solve(options.local_minimizer_options,
                              problem,
                              current_state_.x.data(),
                              &local_summary);

        global_summary->local_minimization_time_in_seconds += WallTimeInSeconds() - t1;

        // check the minimization exited without failure and update state
        // variables if so.
        if (local_summary.termination_type == TerminationType::FAILURE ||
            local_summary.termination_type == TerminationType::USER_FAILURE) {
            global_summary->termination_type = TerminationType::FAILURE;
            global_summary->message = "Initial local minimization iteration failed."
                                              "More details: " + local_summary.message;
            LOG_IF(WARNING, is_not_silent) << "Terminating: " << global_summary->message;
            prepare_final_summary_(global_summary, local_summary);
            return;
        }

        t1 = WallTimeInSeconds();
        if (!Evaluate(problem, current_state_.x, &current_state_, &global_summary->message)) {
            global_summary->termination_type = TerminationType::FAILURE;
            global_summary->message = "Initial cost and jacobian evaluation failed. "
                                              "More details: " + global_summary->message;
            LOG_IF(WARNING, is_not_silent) << "Terminating: " << global_summary->message;
            prepare_final_summary_(global_summary, local_summary);
            return;
        }
        global_summary->cost_evaluation_time_in_seconds += WallTimeInSeconds() - t1;

        ++num_iterations_;
        // initialize values of remaining state variables with current state
        candidate_state_ = current_state_;
        global_minimum_state_ = current_state_;

        if (options.history_save_frequency > 0)
            global_summary->history.push_back(HistoryOutput(num_iterations_, num_stagnant_iterations_, current_state_.x, global_minimum_state_.cost, global_minimum_state_.x));

        // check that initial minimization didn't satisfy termination conditions
        // before entering main loop
        if (check_for_termination_(options, &global_summary->message, &global_summary->termination_type)) {
            prepare_final_summary_(global_summary, local_summary);
            if (internal::IsSolutionUsable(global_summary))
                x = global_minimum_state_.x;
            return;
        }

        // start main loop
        while (true) {
            t1 = WallTimeInSeconds();
            options.step_function->Step(candidate_state_.x.data(), num_parameters);
            global_summary->step_time_in_seconds += WallTimeInSeconds() - t1;

            t1 = WallTimeInSeconds();
            GradientLocalMinimizer::Summary local_summary;
            local_minimizer.Solve(options.local_minimizer_options,
                                  problem,
                                  candidate_state_.x.data(),
                                  &local_summary);
            global_summary->local_minimization_time_in_seconds += WallTimeInSeconds() - t1;

            t1 = WallTimeInSeconds();
            if (!Evaluate(problem, candidate_state_.x, &candidate_state_, &global_summary->message)) {
                global_summary->termination_type = TerminationType::FAILURE;
                global_summary->message = "Cost and jacobian evaluation failed. "
                                                  "More details: " + global_summary->message;
                LOG_IF(WARNING, is_not_silent) << "Terminating: " << global_summary->message;
                prepare_final_summary_(global_summary, local_summary);
                return;
            }
            global_summary->cost_evaluation_time_in_seconds += WallTimeInSeconds() - t1;

            accept = metropolis_(candidate_state_.cost, current_state_.cost);

            if (accept) {
                current_state_ = candidate_state_;
                new_global_min = global_minimum_state_.update(current_state_);
            }

            if (new_global_min) {
                num_stagnant_iterations_ = 0;
            } else {
                ++num_stagnant_iterations_;
            }

            ++num_iterations_;

            if (options.history_save_frequency > 0 && num_iterations_ % options.history_save_frequency == 0)
                global_summary->history.push_back(HistoryOutput(num_iterations_, num_stagnant_iterations_, current_state_.x, global_minimum_state_.cost, global_minimum_state_.x));

            if (check_for_termination_(options, &global_summary->message, &global_summary->termination_type)) {
                prepare_final_summary_(global_summary, local_summary);
                if (internal::IsSolutionUsable(global_summary))
                    x = global_minimum_state_.x;

                if (options.history_save_frequency > 0 && num_iterations_ % options.history_save_frequency != 0)
                    global_summary->history.push_back(HistoryOutput(num_iterations_, num_stagnant_iterations_, current_state_.x, global_minimum_state_.cost, global_minimum_state_.x));
                global_summary->total_time_in_seconds = WallTimeInSeconds() - start_time;
                return;
            }
        }

    };

    bool Basinhopping::check_for_termination_(const Basinhopping::Options &options,
                                              std::string *message,
                                              TerminationType *termination_type) {
        if (global_minimum_state_.cost < options.minimum_cost) {
            *message = "Prescribed minimum cost reached.";
            *termination_type = TerminationType::USER_SUCCESS;
            return true;
        } else if (num_iterations_ >= options.max_iterations) {
            *message = "Maximum number of iterations reached.";
            *termination_type = TerminationType::NO_CONVERGENCE;
            return true;
        } else if (num_stagnant_iterations_ >= options.max_stagnant_iterations) {
            *message = "Maximum number of stagnant iterations reached.";
            *termination_type = TerminationType::CONVERGENCE;
            return true;
        }

        return false;
    };

    void Basinhopping::prepare_final_summary_(Basinhopping::Summary *global_summary,
                                              const GradientLocalMinimizer::Summary &local_summary) {
        global_summary->final_cost = global_minimum_state_.cost;
        global_summary->num_iterations = num_iterations_;
        global_summary->local_minimization_summary = local_summary;
    }

    void Solve(const Basinhopping::Options& options,
               const GradientProblem& problem,
               double* parameters,
               Basinhopping::Summary* summary) {
        Basinhopping solver;
        solver.Solve(options, problem, parameters, summary);
    }

    void dump(const Basinhopping::HistoryOutput &h, HistoryWriter& writer) {
        writer.StartObject();
        writer.String("iteration_number");
        writer.Uint(h.iteration_number);

        writer.String("stagnant_iterations");
        writer.Uint(h.stagnant_iterations);

        writer.String("current_solution");
        writer.StartArray();
        for (auto i = 0; i < h.current_solution.size(); ++i)  {
            writer.Double(h.current_solution[i]);
        }
        writer.EndArray();

        writer.String("best_cost");
        writer.Double(h.best_cost);

        writer.String("best_solution");
        writer.StartArray();
        for (auto i = 0; i < h.best_solution.size(); ++i)  {
            writer.Double(h.best_solution[i]);
        }
        writer.EndArray();
        writer.EndObject();
    }
} // namespace pallas
