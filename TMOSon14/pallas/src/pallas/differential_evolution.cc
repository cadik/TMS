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

#include "pallas/differential_evolution.h"
#include "pallas/internal/solver_utils.h"
#include "pallas/internal/stringprintf.h"
#include "pallas/internal/wall_time.h"

#include "glog/logging.h"


namespace pallas {
    using std::string;

    using ceres::TerminationTypeToString;

    using pallas::internal::StringAppendF;
    using pallas::internal::StringPrintf;
    using pallas::internal::WallTimeInSeconds;


    namespace {

        bool Evaluate(const GradientProblem &problem,
                      Vector &x,
                      double *cost,
                      string *message) {
            if (!problem.Evaluate(x.data(),
                                  cost,
                                  NULL)) {
                *message = "Problem evaluation failed";
                return false;
            }
            return true;
        }
    } // namespace

    DifferentialEvolution::Summary::Summary()
            : termination_type(TerminationType::FAILURE),
              message("pallas::DifferentialEvolution was not called."),
              final_cost(-1.0),
              num_parameters(0),
              num_iterations(0),
              total_time_in_seconds(0.0),
              local_minimization_time_in_seconds(0.0),
              cost_evaluation_time_in_seconds(0.0),
              was_polished(false){

    };

    string DifferentialEvolution::Summary::BriefReport() const {
        return StringPrintf(
                "Pallas differential evolution report: "
                        "iterations: %d, "
                        "final cost: %e, "
                        "termination: %s\n",
                num_iterations,
                final_cost,
                TerminationTypeToString(termination_type));
    };

    string DifferentialEvolution::Summary::FullReport() const {

        string report = string("\nSolver Summary\n\n");

        StringAppendF(&report, "Parameters          %25d\n", num_parameters);

        string mutation_strategy_string = MutationStrategyTypeToString(mutation_strategy);

        StringAppendF(&report, "Mutation strategy     %23s\n",
                      mutation_strategy_string.c_str());

        if (termination_type != TerminationType::FAILURE &&
            termination_type != TerminationType::USER_FAILURE) {
            StringAppendF(&report, "\nFinal cost          %25e\n", final_cost);
        }

        StringAppendF(&report, "\nMinimizer iterations         %16d\n",
                      num_iterations);

        StringAppendF(&report, "\nTime (in seconds):\n");

        StringAppendF(&report, "  Cost evaluation     %23.4f\n",
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

    void DifferentialEvolution::Solve(const DifferentialEvolution::Options& options,
                                      const GradientProblem& problem,
                                      double* parameters,
                                      DifferentialEvolution::Summary* global_summary) {
        double start_time = WallTimeInSeconds();
        double t1;
        bool is_not_silent = !options.is_silent;

        num_parameters_ = static_cast<unsigned int>(problem.NumParameters());
        global_summary->num_parameters = num_parameters_;

        init_member_variables_(options);

        VectorRef x(parameters, num_parameters_);

        GradientLocalMinimizer::Summary local_summary;

        global_summary->mutation_strategy = options.mutation_strategy;
        global_summary->crossover_strategy = options.crossover_strategy;

        std::vector<Vector> scaled_population;
        if (options.history_save_frequency > 0) {
            scaled_population.resize(population_size_);
            for (unsigned int i = 0; i < population_size_; ++i) {
                scaled_population[i].resize(num_parameters_);
            }
        }

        Vector scaled_trial(num_parameters_);
        t1 = WallTimeInSeconds();
        for(unsigned int i = 0; i < population_energies_.size(); ++i) {
            scale_parameters_(population_[i], scaled_trial);
            if (!Evaluate(problem, scaled_trial, &(population_energies_[i]), &global_summary->message)) {
                global_summary->termination_type = TerminationType::FAILURE;
                global_summary->message = "Initial cost evaluation failed. "
                                                  "More details: " + global_summary->message;
                LOG_IF(WARNING, is_not_silent) << "Terminating: " << global_summary->message;
                prepare_final_summary_(global_summary, local_summary);
                return;
            }
        }
        global_summary->cost_evaluation_time_in_seconds += WallTimeInSeconds() - t1;

        pallas::Vector::Index min_idx;

        global_minimum_state_.cost = population_energies_.minCoeff(&min_idx);
        global_minimum_state_.x = population_[min_idx];

        if (options.history_save_frequency > 0) {
            for (size_t i = 0; i < population_.size(); ++i)
                scale_parameters_(population_[i], scaled_population[i]);
            global_summary->history.push_back(HistoryOutput(num_iterations_, scaled_population, global_minimum_state_.cost, global_minimum_state_.x));
        }

        if(check_for_termination_(options, &global_summary->message, &global_summary->termination_type)) {
            prepare_final_summary_(global_summary, local_summary);
            if (internal::IsSolutionUsable(global_summary)||internal::IsSolutionUsable(local_summary))
                x = global_minimum_state_.x;
            return;
        }

        // move fittest individual to first slot
        std::swap(population_[min_idx], population_[0]);
        std::swap(population_energies_[min_idx], population_energies_[0]);

        internal::State current_state(num_parameters_);

        double trial_energy;
        Vector trial(num_parameters_);
        while (true) {
            scale_ = (*random_dither_)();

            for (unsigned int i = 0; i < population_size_; ++i) {
                trial = population_[i];
                mutate_(trial, i);
                ensure_constraint_(trial);
                scale_parameters_(trial, scaled_trial);
                t1 = WallTimeInSeconds();
                if (!Evaluate(problem, scaled_trial, &trial_energy, &global_summary->message)) {
                    global_summary->termination_type = TerminationType::FAILURE;
                    global_summary->message = "Cost evaluation failed. "
                                                      "More details: " + global_summary->message;
                    LOG_IF(WARNING, is_not_silent) << "Terminating: " << global_summary->message;
                    prepare_final_summary_(global_summary, local_summary);
                    return;
                }
                global_summary->cost_evaluation_time_in_seconds += WallTimeInSeconds() - t1;

                if (trial_energy < population_energies_[i]) {
                    population_[i] = trial;
                    population_energies_[i] = trial_energy;

                    if (trial_energy < population_energies_[0]) {
                        population_[0] = trial;
                        population_energies_[0] = trial_energy;
                    }
                }
            }

            scale_parameters_(population_[0], current_state.x);
            current_state.cost = population_energies_[0];

            ++num_iterations_;
            global_minimum_state_.update(current_state);

            update_std_dev_();

            if (options.history_save_frequency > 0 && num_iterations_ % options.history_save_frequency == 0) {
                for (size_t i = 0; i < population_.size(); ++i)
                    scale_parameters_(population_[i], scaled_population[i]);
                global_summary->history.push_back(HistoryOutput(num_iterations_, scaled_population, global_minimum_state_.cost, global_minimum_state_.x));
            }

            if(check_for_termination_(options, &global_summary->message, &global_summary->termination_type)) {
                if (options.polish_output) {
                    t1 = WallTimeInSeconds();
                    GradientLocalMinimizer local_minimizer;
                    local_minimizer.Solve(options.local_minimizer_options,
                                          problem,
                                          global_minimum_state_.x.data(),
                                          &local_summary);
                    global_summary->was_polished = true;
                    global_summary->local_minimization_time_in_seconds = WallTimeInSeconds() - t1;
                }

                t1 = WallTimeInSeconds();
                if (!Evaluate(problem, global_minimum_state_.x, &global_minimum_state_.cost, &global_summary->message)) {
                    global_summary->termination_type = TerminationType::FAILURE;
                    global_summary->message = "Cost evaluation of global minimum state failed after polishing step "
                                                      "More details: " + global_summary->message;
                    LOG_IF(WARNING, is_not_silent) << "Terminating: " << global_summary->message;
                }
                global_summary->cost_evaluation_time_in_seconds += WallTimeInSeconds() - t1;
                prepare_final_summary_(global_summary, local_summary);
                if (internal::IsSolutionUsable(global_summary)||internal::IsSolutionUsable(local_summary))
                    x = global_minimum_state_.x;

                if (options.history_save_frequency > 0 && num_iterations_ % options.history_save_frequency != 0) {
                    for (size_t i = 0; i < population_.size(); ++i)
                        scale_parameters_(population_[i], scaled_population[i]);
                    global_summary->history.push_back(HistoryOutput(num_iterations_, scaled_population, global_minimum_state_.cost, global_minimum_state_.x));
                }
                global_summary->total_time_in_seconds = WallTimeInSeconds() - start_time;
                return;
            }
        }
    };

    void DifferentialEvolution::init_member_variables_(const DifferentialEvolution::Options &options) {
        fractional_std_dev_ = DBL_MAX;

        population_size_ = options.population_size;
        init_mutation_strategy_(options.mutation_strategy);
        init_crossover_strategy_(options.crossover_strategy, options.crossover_probability);
        init_random_dither_(options.dither);

        scoped_ptr<internal::RandomNumberGenerator<double>> tmp_rng(new internal::RandomNumberGenerator<double>());
        swap(random_number_, tmp_rng);

        scoped_ptr<internal::Shuffler> tmp_shuffler(new internal::Shuffler(population_size_));
        swap(shuffler_, tmp_shuffler);

        CHECK(options.upper_bounds.size() == num_parameters_) << "Upper bounds of size" << options.upper_bounds.size()
                                                              << " does not have the length as the number of parameters.";
        upper_bounds_ = options.upper_bounds;

        CHECK(options.lower_bounds.size() == num_parameters_) << "Lower bounds of size" << options.lower_bounds.size()
                                                              << " does not have the length as the number of parameters.";
        lower_bounds_ = options.lower_bounds;

        scale_ = (*random_number_)();

        scale_arg1_ = upper_bounds_ + lower_bounds_;
        scale_arg1_ *= 0.5;
        scale_arg2_ = (upper_bounds_ - lower_bounds_).cwiseAbs();

        population_.resize(population_size_);
        for (unsigned int i = 0; i < population_size_; ++i) {
            population_[i].resize(num_parameters_);
        }
        population_idx_.resize(population_size_);
        population_idx_.setLinSpaced(population_size_, 0, population_size_ - 1);
        shuffler_->Shuffle(population_idx_.data(), 100);

        init_population_(options.population_initialization);

        population_energies_.resize(population_size_);
        population_energies_.setConstant(DBL_MAX);

        num_iterations_ = 0;
    }

    void DifferentialEvolution::init_mutation_strategy_(MutationStrategyType type) {
        scoped_ptr<internal::MutationStrategy> tmp_mutation_strategy(
                internal::MutationStrategy::Create(type));
        swap(mutation_strategy_, tmp_mutation_strategy);
    };

    void DifferentialEvolution::init_crossover_strategy_(CrossoverStrategyType type,
                                                         double crossover_probability) {
        scoped_ptr<internal::CrossoverStrategy> tmp_crossover_strategy(
                internal::CrossoverStrategy::Create(type, crossover_probability, num_parameters_));
        swap(crossover_strategy_, tmp_crossover_strategy);
    };

    void DifferentialEvolution::init_random_dither_(const Vector2d &dither) {
        scoped_ptr<internal::RandomNumberGenerator<double>> tmp_rd(
                new internal::RandomNumberGenerator<double>(dither[0], dither[1]));
        swap(random_dither_, tmp_rd);
    };

    void DifferentialEvolution::init_population_(PopulationInitializationType type) {

        Vector v(num_parameters_);

        if (type == LATIN_HYPERCUBE) {
            double segsize = 1.0 / population_size_;
            std::vector<Vector> rdrange(population_size_);
            for (unsigned int i = 0; i < population_size_; ++i) {
                rdrange[i].resize(num_parameters_);
            }
            Vector arange(population_size_);
            arange.setLinSpaced(population_size_, 0.0, 1.0);
            Vector arange_vector(num_parameters_);

            for (unsigned int i = 0; i < population_size_; ++i) {
                for (unsigned int j = 0; j < num_parameters_; ++j) {
                    v[j] = (*random_number_)();
                }
                rdrange[i] = (segsize * v);
                arange_vector.setConstant(arange[i]);
                rdrange[i] += arange_vector;
                for (unsigned int j = 0; j < num_parameters_; ++j) {
                    if (rdrange[i][j] > 1.0) {
                        rdrange[i][j] = 1.0;
                    }
                    else if (rdrange[i][j] < 0.0) {
                        rdrange[i][j] = 0.0;
                    }
                }
            }

            int idx;
            for (unsigned int j = 0; j < num_parameters_; ++j) {
                shuffler_->Shuffle(population_idx_.data());
                for (unsigned int i = 0; i < population_size_; ++i) {
                    idx = population_idx_[i];
                    population_[i][j] = rdrange[idx][j];
                }
            }
        }
        else if (type == RANDOM) {
            for (unsigned int i = 0; i < population_size_; ++i) {
                for (unsigned int j = 0; j < v.size(); ++j) {
                    v[j] = (*random_number_)();
                }
                population_[i] = v;
            }
        }
        else
            LOG(ERROR) << "Unknown population initialization type: " << PopulationInitializationTypeToString(type);
    };


    void DifferentialEvolution::scale_parameters_(const Vector& trial_in, Vector& trial_out) {
        for (unsigned int i = 0; i < num_parameters_; ++i) {
            trial_out[i] = scale_arg1_[i] + (trial_in[i] - 0.5) * scale_arg2_[i];
        }
    };

    void DifferentialEvolution::unscale_parameters_(const Vector& trial_in, Vector& trial_out) {
        for (unsigned int i = 0; i < num_parameters_; ++i) {
            trial_out[i] = (trial_in[i] - scale_arg1_[i]) / scale_arg2_[i] + 0.5;
        }
    };

    void DifferentialEvolution::ensure_constraint_(Vector& trial) {
        for (unsigned int i = 0; i < trial.size(); ++i) {
            if(trial[i] < 0.0 || trial[i] > 1.0) {
                trial[i] = (*random_number_)();
            }
        }
    };

    void DifferentialEvolution::update_std_dev_() {
        double mean = population_energies_.mean();
        double temp, variance = 0.0;
        for(unsigned int i = 0; i < population_size_; ++i) {
            temp = (population_energies_[i] - mean);
            variance += temp * temp;
        }
        variance /= (population_size_ - 1);
        fractional_std_dev_ = std::sqrt(variance) / std::abs(mean);
    }

    void DifferentialEvolution::mutate_(Vector& candidate, unsigned int idx) {
        shuffler_->Shuffle(population_idx_.data());

        crossover_strategy_->Crossover(candidate, mutation_strategy_->get_bprime(population_,
                                                                                 idx,
                                                                                 population_idx_.data(),
                                                                                 scale_));
    };

    bool DifferentialEvolution::check_for_termination_(const DifferentialEvolution::Options& options,
                                string *message,
                                TerminationType * termination_type) {

        if (global_minimum_state_.cost < options.minimum_cost) {
            *message = "Prescribed minimum cost reached.";
            *termination_type = TerminationType::USER_SUCCESS;
            return true;
        } else if (num_iterations_ >= options.max_iterations) {
            *message = "Maximum number of iterations reached.";
            *termination_type = TerminationType::NO_CONVERGENCE;
            return true;
        } else if (fractional_std_dev_ < options.tolerance) {
            *message = "Fractional standard deviation of population less than specified tolerance.";
            *termination_type = TerminationType::CONVERGENCE;
            return true;
        }
        return false;
    };

    void DifferentialEvolution::prepare_final_summary_(DifferentialEvolution::Summary *global_summary,
                                const GradientLocalMinimizer::Summary &local_summary) {
        global_summary->final_cost = global_minimum_state_.cost;
        global_summary->num_iterations = num_iterations_;
        global_summary->local_minimization_summary = local_summary;
    };

    void Solve(const DifferentialEvolution::Options& options,
               const GradientProblem& problem,
               double* parameters,
               DifferentialEvolution::Summary* summary) {
        DifferentialEvolution solver;
        solver.Solve(options, problem, parameters, summary);
    };

    void dump(const DifferentialEvolution::HistoryOutput &h, HistoryWriter& writer) {
        writer.StartObject();
        writer.String("iteration_number");
        writer.Uint(h.iteration_number);

        writer.String("population");
        writer.StartArray();
        for (auto& individual : h.population) {
            writer.StartArray();
            for (auto i = 0; i < individual.size(); ++i)  {
                writer.Double(individual[i]);
            }
            writer.EndArray();
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
