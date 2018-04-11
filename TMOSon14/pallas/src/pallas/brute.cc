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

#include <cfloat>
#include "pallas/brute.h"
#include "pallas/internal/state.h"
#include "pallas/internal/solver_utils.h"
#include "pallas/internal/stringprintf.h"
#include "pallas/internal/wall_time.h"

namespace pallas {

    using std::string;

    using pallas::internal::StringAppendF;
    using pallas::internal::StringPrintf;
    using pallas::internal::WallTimeInSeconds;

    namespace {

        bool Evaluate(const GradientProblem &problem,
                      Vector &x,
                      internal::State *state,
                      std::string *message) {
            if(!problem.Evaluate(x.data(), &(state->cost), state->gradient.data())) {
                *message = "Problem evaluation failed.";
                return false;
            }
            return true;
        }

    } // namespace

    Brute::Summary::Summary()
            : termination_type(TerminationType::FAILURE),
              message("pallas::Brute was not called."),
              final_cost(-1.0),
              num_parameters(0),
              num_iterations(0),
              total_time_in_seconds(0.0),
              local_minimization_time_in_seconds(0.0),
              cost_evaluation_time_in_seconds(0.0),
              permutation_build_time_in_seconds(0.0),
              was_polished(false){

    };

    std::string Brute::Summary::BriefReport() const {
        return StringPrintf(
                "Pallas brute report: "
                        "iterations: %d, "
                        "final cost: %e, "
                        "termination: %s\n",
                num_iterations,
                final_cost,
                TerminationTypeToString(termination_type));
    };

    std::string Brute::Summary::FullReport() const {

        string report = string("\nSolver Summary\n\n");

        StringAppendF(&report, "Parameters          %25d\n", num_parameters);

        StringAppendF(&report, "\n");

        if (termination_type != TerminationType::FAILURE &&
            termination_type != TerminationType::USER_FAILURE) {
            StringAppendF(&report, "Final cost          %25e\n", final_cost);
        }

        StringAppendF(&report, "\nTotal iterations         %20d\n",
                      num_iterations);

        StringAppendF(&report, "\nTime (in seconds):\n");

        StringAppendF(&report, "  Calculate permutations     %16.4f",
                      permutation_build_time_in_seconds);

        StringAppendF(&report, "\n  Cost evaluation     %23.4f\n",
                      cost_evaluation_time_in_seconds);

        if (was_polished) {
            StringAppendF(&report, "  Local minimization   %22.4f\n",
                          local_minimization_time_in_seconds);
        }

        StringAppendF(&report, "  Total               %23.4f\n\n",
                      total_time_in_seconds);

        StringAppendF(&report, "Termination: %2s (%s)\n",
                      TerminationTypeToString(termination_type), message.c_str());
        return report;
    };

    void Brute::Solve(const Brute::Options options,
                      const GradientProblem &problem,
                      const std::vector<Brute::ParameterRange> &parameter_ranges,
                      double* parameters,
                      Brute::Summary *global_summary) {
        double start_time = WallTimeInSeconds();
        double t1;

        bool is_not_silent = !options.is_silent;

        const unsigned int num_parameters = static_cast<unsigned int>(problem.NumParameters());

        VectorRef x(parameters, num_parameters);

        t1 = WallTimeInSeconds();
        std::vector<Vector> expanded_ranges = expand_parameter_ranges_(parameter_ranges);

        std::vector<Vector> permutations = build_permutations_(expanded_ranges);
        global_summary->permutation_build_time_in_seconds = WallTimeInSeconds() - t1;

        global_summary->num_parameters = num_parameters;
        global_summary->num_iterations = static_cast<unsigned int>(permutations.size());

        internal::State current_state(num_parameters);
        internal::State global_minimum_state(num_parameters);
        global_minimum_state.cost = DBL_MAX;

        t1 = WallTimeInSeconds();
        for (unsigned int i = 0; i < permutations.size(); ++i) {
            current_state.x = permutations[i];
            if (!Evaluate(problem, current_state.x, &current_state, &global_summary->message)) {
                global_summary->termination_type = TerminationType::FAILURE;
                global_summary->message = "Initial cost and jacobian evaluation failed. "
                                                  "More details: " + global_summary->message;
                LOG_IF(WARNING, is_not_silent) << "Terminating: " << global_summary->message;
                return;
            }
            global_minimum_state.update(current_state);

            if (options.history_save_frequency > 0 && i % options.history_save_frequency == 0)
                global_summary->history.push_back(HistoryOutput(i, current_state.x, global_minimum_state.cost, global_minimum_state.x));
        }
        global_summary->cost_evaluation_time_in_seconds = WallTimeInSeconds() - t1;

        if(options.polish_output) {
            t1 = WallTimeInSeconds();
            GradientLocalMinimizer local_minimizer;
            local_minimizer.Solve(options.local_minimizer_options,
                                  problem,
                                  global_minimum_state.x.data(),
                                  &global_summary->local_minimization_summary);
            global_summary->local_minimization_time_in_seconds = WallTimeInSeconds() - t1;

            t1 = WallTimeInSeconds();
            if (!Evaluate(problem, global_minimum_state.x, &global_minimum_state, &global_summary->message)) {
                global_summary->termination_type = TerminationType::FAILURE;
                global_summary->message = "Cost evaluation of global minimum state failed after polishing step "
                                                  "More details: " + global_summary->message;
                LOG_IF(WARNING, is_not_silent) << "Terminating: " << global_summary->message;
                return;
            }
            global_summary->cost_evaluation_time_in_seconds += WallTimeInSeconds() - t1;
            global_summary->was_polished = true;
        }

        global_summary->message = "Specified search of parameter space successfully completed.";
        global_summary->termination_type = TerminationType::USER_SUCCESS;
        global_summary->final_cost = global_minimum_state.cost;

        x = global_minimum_state.x;
        global_summary->total_time_in_seconds = WallTimeInSeconds() - start_time;
    };

    std::vector<Vector> Brute::expand_parameter_ranges_(
            const std::vector <Brute::ParameterRange> &parameter_ranges) {

        std::vector<Vector> expanded_ranges(parameter_ranges.size());
        double start, stop;
        int size;

        for(unsigned int i = 0; i < expanded_ranges.size(); ++i) {
            Vector v;
            start = parameter_ranges[i].start;
            stop = parameter_ranges[i].stop;
            size = parameter_ranges[i].size;
            expanded_ranges[i] = v.setLinSpaced(size, start, stop);
        }

        return expanded_ranges;
    }

    std::vector<Vector> Brute::build_permutations_(
            const std::vector<Vector> &expanded_ranges) {

        const unsigned int num_parameters = static_cast<unsigned int>(expanded_ranges.size());

        Eigen::VectorXi cyclic_counter(num_parameters);
        cyclic_counter.setZero();

        unsigned int num_permutations = 1;

        for (unsigned int i = 0; i < expanded_ranges.size(); ++i) {
            num_permutations *= expanded_ranges[i].size();
        }

        std::vector<Vector> output;
        output.resize(num_permutations);
        for (unsigned int i = 0; i < num_permutations; ++i) {
            output[i].resize(num_parameters);
        }

        unsigned int col_idx = num_parameters - 1;

        Vector current_permutation(num_parameters);
        for (unsigned int i = 0; i < num_parameters; ++i) {
            current_permutation[i] = expanded_ranges[i][0];
        }
        output[0] = current_permutation;

        bool completed_cycle = false;
        for (unsigned int i = 1; i < num_permutations; ++i) {
            ++cyclic_counter[col_idx];
            cyclic_counter[col_idx] %= expanded_ranges[col_idx].size();

            current_permutation[col_idx] = expanded_ranges[col_idx][cyclic_counter[col_idx]];

            if (cyclic_counter[col_idx] == 0)
                completed_cycle = true;

            while (completed_cycle && i < num_permutations) {
                for (unsigned int j = col_idx - 1; j>=0; --j) {
                    ++cyclic_counter[j];
                    cyclic_counter[j] %= expanded_ranges[j].size();

                    current_permutation[j] = expanded_ranges[j][cyclic_counter[j]];
                    if (cyclic_counter[j] != 0) {
                        completed_cycle = false;
                        break;
                    }
                }
            }
            output[i] = current_permutation;
        }
        return output;
    };

    void Solve(const Brute::Options options,
               const GradientProblem& problem,
               const std::vector<Brute::ParameterRange> &parameter_ranges,
               double* parameters,
               Brute::Summary *global_summary) {
        Brute solver;
        solver.Solve(options, problem, parameter_ranges, parameters, global_summary);
    }

    void dump(const Brute::HistoryOutput &h, HistoryWriter& writer) {
        writer.StartObject();
        writer.String("iteration_number");
        writer.Uint(h.iteration_number);

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
