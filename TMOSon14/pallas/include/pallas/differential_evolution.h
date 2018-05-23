/*!
* \file differential_evolution.h
*
* \author Ryan Latture
* \date 5-18-15
*
* This file contains a C++ implementation of the differential algorithm due to Storn 
* and Price: Differential Evolution - a Simple and Efficient Heuristic for Global Optimization over Continuous Spaces.
* 
* This code is based on the SciPy implementation of differential evolution found in scipy.optimize and relies on
* the local optimization algorithms in the Google Ceres project.
*/

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

#ifndef PALLAS_DIFFERENTIAL_EVOLUTION_H
#define PALLAS_DIFFERENTIAL_EVOLUTION_H

#include <cfloat>

#include "pallas/history_concept.h"
#include "pallas/types.h"
#include "pallas/internal/crossover_strategy.h"
#include "pallas/internal/mutation_strategy.h"
#include "pallas/internal/shuffler.h"
#include "pallas/internal/state.h"

namespace pallas {

    /**
     * @brief Minimizes an objective function by continuously evolving a population of candidate solutions.
     * @details <B>Example</B>
     * @code
    #include "glog/logging.h"

    // Each solver is defined in its own header file.
    // include the solver you wish you use:
    #include "pallas/differential_evolution.h"

    // define a problem you wish to solve by inheriting
    // from the pallas::GradientCostFunction interface
    // and implementing the Evaluate and NumParameters methods.
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

    int main(int argc, char** argv) {
        google::InitGoogleLogging(argv[0]);

        // define the starting point for the optimization
        double parameters[2] = {-1.2, 0.0};

        // set up global optimizer options only initialization
        // is need to accept the default options
        pallas::DifferentialEvolution::Options options;

        // differential evolution requires that upper
        // and lower bounds are placed on the parameters
        // in the form of pallas::Vector which is a
        // just a typedef of an Eigen::VectorXd.
        pallas::Vector upper_bounds(2);
        upper_bounds << 5, 5;

        pallas::Vector lower_bounds(2);
        lower_bounds << -5, -5;

        options.upper_bounds = upper_bounds;
        options.lower_bounds = lower_bounds;

        // initialize a summary object to hold the
        // optimization details
        pallas::DifferentialEvolution::Summary summary;

        // create a problem from your cost function
        pallas::GradientProblem problem(new Rosenbrock());

        // solve the problem and store the optimal position
        // in parameters and the optimization details in
        // the summary
        pallas::Solve(options, problem, parameters, &summary);

        std::cout << summary.FullReport() << std::endl;
        std::cout << "Global minimum found at:" << std::endl;
        std::cout << "\tx: " << parameters[0] << "\ty: " << parameters[1] << std::endl;

        return 0;
    }
     * @endcode
     */
    class DifferentialEvolution {
    public:
        /**
         * Configurable options for modifying the default behaviour of the differential evolution algorithm.
         */
        struct Options {
            /**
             * @brief Default constructor
             * @details This tries to set up reasonable defaults for optimization. It is highly
             * recommended that the user read and overwrite the defaults based on the cost function.
             */
            Options() {
                local_minimizer_options = GradientLocalMinimizer::Options();
                mutation_strategy = BEST_1;
                crossover_strategy = BINOMIAL;
                population_initialization = LATIN_HYPERCUBE;
                max_iterations = 1000;
                population_size = 15;
                tolerance = 0.01;
                minimum_cost = -DBL_MAX;
                dither << 0.5, 1.0;
                crossover_probability = 0.7;
                is_silent = true;
                polish_output = false;
                history_save_frequency = 0;
            };

            /**
             * Contains any changes to the default options for the local minimization algorithm.
             * See the documentation for ceres::GradientProblemSolver::Options for relevant options
             */
            GradientLocalMinimizer::Options local_minimizer_options;

            /**
             * The mutation strategy to use. Should be one of `BEST_1`, `RAND_1`, `RAND_TO_BEST_1`, `BEST_2`, or `RAND_2`.
             */
            MutationStrategyType mutation_strategy;

            /**
             * The crossover strategy to use. Should be one of `BINOMIAL` or `EXPONENTIAL`.
             */
            CrossoverStrategyType crossover_strategy;

            /**
             * Specify how the population initialization is performed. 
             * Should be one of: `LATIN_HYPERCUBE` or `RANDOM`. The default is `LATIN_HYPERCUBE`. 
             * Latin Hypercube sampling tries to maximize coverage of the available parameter space. 
             * `RANDOM` initializes the population randomly - this has the drawback that clustering 
             * can occur, preventing the whole of parameter space being covered.
             */
            PopulationInitializationType population_initialization;

            /**
             * The maximum number of times the entire population is evolved.
             */
            unsigned int max_iterations;


            /**
             * The total number of individuals in the population.
             */
            unsigned int population_size;

            /**
             * User specified minimum cost. Minimization will halt if `minimum_cost` is reached.
             */
            double minimum_cost;

            /**
             * It should be in the range [0, 2]. 
             * Dithering randomly changes the mutation constant on a generation by generation basis. 
             * The mutation constant for that generation is taken from U[min, max). Dithering can help speed convergence significantly. 
             * Increasing the mutation constant increases the search radius, but will slow down convergence.
             */
            Vector2d dither;

            /**
             * Upper bounds for variables.
             */
            Vector upper_bounds;

            /**
             * Lower bounds for variables
             */
            Vector lower_bounds;

            /**
             * The recombination constant, should be in the range [0, 1]. 
             * Increasing this value allows a larger number of mutants to 
             * progress into the next generation, but at the risk of population stability.
             */
            double crossover_probability;

            /**
             * When the mean of the population energies, multiplied by tolerance, 
             * divided by the standard deviation of the population energies is greater 
             * than 1 the solving process terminates: `convergence = mean(pop) * tol / stdev(pop) > 1`.
             */
            double tolerance;

            /**
             * Whether the global minimum found through differential evolution should be subjected to a local minimization "polishing" step before returning the result.
             */
            bool polish_output;

            /**
             * Whether to log failure information relating the to global optimization algorithm using glog.
             */
            bool is_silent;

            /**
             * Frequency to save the state of the system. Values will be appended to a `HistorySeries` contained
             * in the optimization summary. Default is 0. If 0 then history is not saved. Otherwise, the state of
             * the system will be appended to the series when `i % history_save_frequency == 0`. If there are a
             * large number of iterations this can lead to a lot of data being stored in memory.
             */
            unsigned int history_save_frequency;
        };

        /**
         * @brief Contains a summary of the optimization.
         * @details This struct contains the result of the optimization and has convenience methods for printing reports of a completed optimization.
         */
        struct Summary{
            /**
            * @brief Default constructor
            */
            Summary();

            std::string BriefReport() const;/**<A brief one line description of the state of the solver after termination.*/

            std::string FullReport() const;/**<A full multi-line description of the state of the solver after termination.*/

            TerminationType termination_type;/**<Reason optimization was terminated*/

            /**
             * The mutation strategy to use. Should be one of `BEST_1`, `RAND_1`, `RAND_TO_BEST_1`, `BEST_2`, or `RAND_2`.
             */
            MutationStrategyType mutation_strategy;

            /**
             * The crossover strategy to use. Should be one of `BINOMIAL` or `EXPONENTIAL`.
             */
            CrossoverStrategyType crossover_strategy;

            std::string message;/**<Message describing why the solver terminated.*/

            double final_cost;/**<Cost of the problem (value of the objective function) after the optimization.*/

            GradientLocalMinimizer::Summary local_minimization_summary;/**<Summary from the local minimization polishing step (if performed).*/

            unsigned int num_parameters;/**<Number of parameters in the problem.*/

            unsigned int num_iterations;/**<Number of times the population was evolved*/

            double total_time_in_seconds;/**<total time elapsed in global minimizer*/

            double local_minimization_time_in_seconds;/**<time spent in local minimizer (if polishing was done).*/

            double cost_evaluation_time_in_seconds;/**<time spent evaluating cost function (outside local minimization)*/

            bool was_polished;/**<whether global minimum was polished after differential evolution completed*/

            HistorySeries history;/**<History of the system saved on the interval specified by the `history_save_frequency` option.*/
        };

        /**
         * @brief Stores information about the state of the system for at a given iteration number
         */
        struct HistoryOutput {

            /**
             * @brief Constructor
             *
             * @param iteration_number unsigned int. The number of global optimization iterations that have elapsed.
             * @param population std::vector<Vector>. Candidate solutions for the current iteration.
             * @param best_cost double. Cost associated with the best individual found at any iteration thus far during optimization.
             * @param best_solution Vector. Best solution found at any iteration thus far during optimization.
             */
            HistoryOutput(unsigned int iteration_number,
                          const std::vector<Vector> &population,
                          double best_cost,
                          const Vector &best_solution)
                    : iteration_number(iteration_number),
                      population(population),
                      best_cost(best_cost),
                      best_solution(best_solution) {}
            unsigned int iteration_number;/**<The number of global optimization iterations that have elapsed.*/
            std::vector<Vector> population;/**<Candidate solutions for the current iteration.*/
            double best_cost;/**<Cost associated with the best individual found at any iteration thus far during optimization.*/
            Vector best_solution;/**<Best solution found at any iteration thus far during optimization.*/
        };

        /**
         * @brief Default constructor
         */
        DifferentialEvolution() {};

        /**
         * @brief Minimizes the specified gradient problem.
         * @details The specified options are used to setup a differential evolution instance which
         * is then used to minimize the GradientProblem. The optimal solution is stored
         * in `parameters` and a summary of the global optimization can be found in `summary`.
         * 
         * @param options pallas::DifferentialEvolution::Options. Options used to configure the optimization.
         * @param problem pallas::GradientProblem. The problem to optimize.
         * @param parameters double*. The starting point for further optimization.
         * @param summary DifferentialEvolution::Summary*. Summary instance to store the optimization details.
         */        
        void Solve(const DifferentialEvolution::Options& options,
                   const GradientProblem& problem,
                   double* parameters,
                   DifferentialEvolution::Summary* global_summary);
    private:
        /**
         * @brief Initializes the state of member variables based on the optimizer options.
         */
        void init_member_variables_(const DifferentialEvolution::Options& options);

        /**
         * @brief Creates the mutation strategy member function.
         */
        void init_mutation_strategy_(MutationStrategyType type);

        /**
         * @brief Creates the crossover strategy member function.
         */
        void init_crossover_strategy_(CrossoverStrategyType type,
                                      double crossover_probability);

        /**
         * @brief Creates the random dither member function.
         */        
        void init_random_dither_(const Vector2d& dither);

        /**
         * @brief Initializes the population using either Latin Hypercube samples or randomly
         */
        void init_population_(PopulationInitializationType type);

        /**
         * @brief Scales the parameters from local to global coordinates.
         * @details Within the differential algorithm, all variables are scaled to lie on the range
         * `[0,1]`. This function returns the parameters back to their original scale based on the 
         * upper and lower ranges provided for the functions.
         * 
         * @param trial_in pallas::Vector. The vector in local space.
         * @param trial_out pallas::Vector. The vector to write the position in global space to.
         */
        void scale_parameters_(const Vector& trial_in, Vector& trial_out);

        /**
         * @brief Scales the parameters from global to local coordinates.
         * @details Within the differential algorithm, all variables are scaled to lie on the range
         * `[0,1]`. This function converts parameters back to fall between 0 and 1 based on the 
         * upper and lower ranges provided for the functions.
         * 
         * @param trial_in pallas::Vector. The vector in global space.
         * @param trial_out pallas::Vector. The vector to write the position in local space to.
         */  
        void unscale_parameters_(const Vector& trial_in, Vector& trial_out);

        /**
         * @brief Ensures the trail candidate's parameters fall within the bounds.
         * @details If a parameters falls outside of the bounds a random number is chosen to replace it.
         * 
         * @param trial pallas::Vector&. Trial candidate to enforce boundary constraints on.
         */
        void ensure_constraint_(Vector& trial);

        /**
         * @brief Update the fractional standard deviation of the population.
         * @details When the fraction standard deviation of the population falls below 
         * the tolerance specified in the DifferentialEvolution::Options variable
         * optimization is terminated.
         */
        void update_std_dev_();

        /**
         * @brief Produce a new candidate from the `ith` member of the population. 
         * 
         * @param candidate pallas::Vector&. Member of population to mutate.
         * @param int unsigned int. Index candidate appears in the population.
         */
        void mutate_(Vector& candidate, unsigned int idx);

        /**
         * @brief Checks to see if any termination conditions were met.
         * 
         * @param options pallas::DifferentialEvolution::Options. Options used to configure the optimization.
         * @param message std::string*. If a termination condition is met, a message describing the satisfied condition is stored in the variable.
         * @param termination_type pallas::TerminationType*. This
         * @return Returns `true` if a termination condition was meet, `false` otherwise. 
         */
        bool check_for_termination_(const DifferentialEvolution::Options& options,
                                    std::string *message,
                                    TerminationType * termination_type);

        /**
         * @brief Updates the global summary before exiting the differential evolution algorithm.
         */
        void prepare_final_summary_(DifferentialEvolution::Summary *global_summary,
                                    const GradientLocalMinimizer::Summary &local_summary);

        scoped_ptr<internal::MutationStrategy> mutation_strategy_;
        scoped_ptr<internal::CrossoverStrategy> crossover_strategy_;
        Vector upper_bounds_;
        Vector lower_bounds_;
        scoped_ptr< internal::RandomNumberGenerator< double > > random_number_;
        scoped_ptr< internal::RandomNumberGenerator< double > > random_dither_;
        scoped_ptr< internal::Shuffler > shuffler_;
        double scale_; /**<Random mutation constant for each iteration generated by `random_dither_`*/
        double fractional_std_dev_;/**<The fractional standard deviation of the population. Controls convergence.*/
        Vector scale_arg1_;/**<Precomputed parameter to scale between global and local parameter space.*/
        Vector scale_arg2_;/**<Precomputed parameter to scale between global and local parameter space.*/
        Vector population_energies_;/**<Cost associated with each member of the current population.*/
        std::vector<Vector> population_;/**<All solutions currently being evolved.*/
        Eigen::VectorXi population_idx_;/**<Randomly shuffled indices of the population used for mutation.*/
        internal::State global_minimum_state_;/**The best solution found during any stage of the optimization. It is this value that is returned when minimization concludes.*/
        unsigned int num_parameters_;/**<The number of degrees of freedom each candidate solution possesses.*/
        unsigned int population_size_;/**<The number total number of individuals currently being evolved.*/
        unsigned int num_iterations_;/**<The number of differential evolution iterations the global optimizer has performed.*/
    };

    /**
     * @brief Helper function which avoids going through the interface of the pallas::DifferentialEvolution class.
     * @details The specified options are used to setup a differential evolution instance which
     * is then used to minimize the GradientProblem. The optimal solution is stored
     * in `parameters` and a summary of the global optimization can be found in `summary`.
     * 
     * @param options pallas::DifferentialEvolution::Options. Options used to configure the optimization.
     * @param problem pallas::GradientProblem. The problem to optimize.
     * @param parameters double*. The starting point for further optimization.
     * @param summary DifferentialEvolution::Summary*. Summary instance to store the optimization details.
     */
    void Solve(const DifferentialEvolution::Options& options,
               const GradientProblem& problem,
               double* parameters,
               DifferentialEvolution::Summary* summary);

    /**
     * @brief Dumps the system state contained in the history output into the stream contained by the writer.
     *
     * @param h DifferentialEvolution::HistoryOutput. State of the system for a specific iteration.
     * @param writer HistoryWriter. Object responsible for writing the history output to a stream.
     */
    void dump(const DifferentialEvolution::HistoryOutput &h, HistoryWriter& writer);

} // namespace pallas

#endif // PALLAS_DIFFERENTIAL_EVOLUTION_H
