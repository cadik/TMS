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

#ifndef PALLAS_SIMULATED_ANNEALING_H
#define PALLAS_SIMULATED_ANNEALING_H

#include <cfloat>

#include "pallas/cooling_schedule.h"
#include "pallas/history_concept.h"
#include "pallas/step_function.h"
#include "pallas/types.h"
#include "pallas/internal/state.h"
#include "pallas/internal/metropolis.h"
#include "pallas/scoped_ptr.h"


namespace pallas {
    /**
     * @brief Minimizes a function using simulated annealing.
     * @details Uses simulated annealing, a random algorithm that uses no 
     * derivative information from the function being optimized. Other names 
     * for this family of approaches include: “Monte Carlo”, “Metropolis”, “Metropolis-Hastings”, etc.\n
     * <B>Example</B>
     * @code
     #include "glog/logging.h"
     
     // Each solver is defined in its own header file.
     // include the solver you wish you use:
     #include "pallas/simulated_annealing.h"
     
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
         pallas::SimulatedAnnealing::Options options;
     
         // Increase the number of iterations
         // SA often requires many iterations to get
         // with 1 to 2 significant figures of the
         // optimal solution
         options.max_iterations = 1000;
         options.dwell_iterations = 1000;
         options.max_stagnant_iterations = 1000;
     
         // set a high initial temperature to allow
         // the SA algorithm to fully explore the
         // parameter space
         options.cooling_schedule_options.initial_temperature = 1000;
     
         // quit the optimization of cost gets within
         // 3 significant figures of global minimum
         options.minimum_cost = 0.001;
     
         // define custom step function which will bound the
         // randomized candidate solution in order to limit
         // the search and speed up convergence
         double upper_bounds [2] = {5, 5};
         double lower_bounds [2] = {-5, -5};
     
         unsigned int num_parameters = 2;
     
         double step_size = 0.1;
     
         pallas::scoped_ptr<pallas::StepFunction> step_function (new pallas::BoundedStepFunction(step_size,
                                                                                                 upper_bounds,
                                                                                                 lower_bounds,
                                                                                                 num_parameters));
         options.set_step_function(step_function);
     
         // initialize a summary object to hold the
         // optimization details
         pallas::SimulatedAnnealing::Summary summary;
     
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
    class SimulatedAnnealing {
    public:
        /**
         * Configurable options for modifying the default behaviour of the simulated annealing algorithm.
         */
        struct Options {
            Options () {
                cooling_schedule_options = CoolingSchedule::Options();
                local_minimizer_options = GradientLocalMinimizer::Options();
                scoped_ptr<StepFunction> default_step(new DefaultStepFunction(1.0));
                step_function.swap(default_step);
                max_iterations = 1000;
                max_stagnant_iterations = 100;
                dwell_iterations = 20;
                minimum_cost = -DBL_MAX;
                is_silent = true;
                polish_output = false;
                history_save_frequency = 0;
            };

            /**
             * @brief Convenience function for changing the default step function.
             * @details This function simply swaps the scoped_ptr to the user 
             * defined `StepFunction` with the scoped_ptr to the `StepFunction`
             * held within the `pallas::SimulatedAnnealing::Options` struct.
             * 
             * @param user_step_function pallas::scoped_ptr<pallas::StepFunction>. This function generates randomized candidate solutions based on the current position.
             */
            void set_step_function(scoped_ptr<StepFunction>& user_step_function) {
                swap(user_step_function, step_function);
            };

            /**
             * The annealing schedule to use. Must be one of ‘FAST’, ‘CAUCHY’ or ‘BOLTZMANN’. see pallas::CoolingSchedule.
             */
            CoolingSchedule::Options cooling_schedule_options;

            /**
             * Contains any changes to the default options for the local minimization algorithm.
             * See the documentation for ceres::GradientProblemSolver::Options for relevant options
             */
            GradientLocalMinimizer::Options local_minimizer_options;

            /**
             * Function that produces randomized candidate solutions.
             */
            scoped_ptr<StepFunction> step_function;

            /**
             * Maximum number of simulated annealing iterations.
             */
            unsigned int max_iterations;

            /**
             * Maximum number sequential simulated annealing iterations allowed without finding a new global minimum.
             */
            unsigned int max_stagnant_iterations;

            /**
             * The number candidate solutions to generate before decreasing the temperature.
             */
            unsigned int dwell_iterations;

            /**
             * User specified minimum cost. Minimization will halt if `minimum_cost` is reached.
             */
            double minimum_cost;

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
        struct Summary {
            /**
             * @brief Default constructor
             */
            Summary();

            std::string BriefReport() const;/**<A brief one line description of the state of the solver after termination.*/

            std::string FullReport() const;/**<A full multi-line description of the state of the solver after termination.*/

            TerminationType termination_type;/**<Reason optimization was terminated*/

            /**
             * The annealing schedule to used. Must be one of ‘FAST’, ‘CAUCHY’ or ‘BOLTZMANN’. see pallas::CoolingSchedule.
             */
            CoolingScheduleType cooling_schedule;

            std::string message;/**<Message describing why the solver terminated.*/

            double initial_cost;/**<Cost of the problem (value of the objective function) before the optimization.*/

            double final_cost;/**<Cost of the problem (value of the objective function) after the optimization.*/

            GradientLocalMinimizer::Summary local_minimization_summary;/**<Summary from the local minimization polishing step (if performed).*/

            unsigned int num_parameters;/**<Number of parameters in the problem.*/

            unsigned int num_iterations;/**<Number of iterations*/

            double total_time_in_seconds;/**<total time elapsed in global minimizer*/

            double local_minimization_time_in_seconds;/**<time spent in local minimizer (if polishing was done).*/

            double step_time_in_seconds;/**<time spent calling step function*/

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
             * @param stagnant_iterations unsigned int. The number of iterations that have elapsed without finding a new global minimum.
             * @param temperature double. Current temperature of the system.
             * @param current_solution Vector. Candidate solution vector for the current iteration.
             * @param best_cost double. Cost associated with the best solution found at any iteration thus far during optimization.
             * @param best_solution Vector. Best solution found at any iteration thus far during optimization.
             */
            HistoryOutput(unsigned int iteration_number,
                          unsigned int stagnant_iterations,
                          double temperature,
                          const Vector &current_solution,
                          double best_cost,
                          const Vector &best_solution)
                    : iteration_number(iteration_number),
                      stagnant_iterations(stagnant_iterations),
                      temperature(temperature),
                      current_solution(current_solution),
                      best_cost(best_cost),
                      best_solution(best_solution) {}
            unsigned int iteration_number;/**<The number of global optimization iterations that have elapsed.*/
            unsigned int stagnant_iterations;/**<The number of iterations that have elapsed without finding a new global minimum.*/
            double temperature;/**<Current temperature of the system.*/
            Vector current_solution;/**<Candidate solution for the current iteration.*/
            double best_cost;/**<Cost associated with the best solution found at any iteration thus far during optimization.*/
            Vector best_solution;/**<Best solution found at any iteration thus far during optimization.*/
        };

        /**
         * @brief Default constructor
         */
        SimulatedAnnealing() {};

        /**
         * @brief Minimizes the specified gradient problem.
         * @details The specified options are used to setup a simulated annealing instance which
         * is then used to minimize the GradientProblem. The optimal solution is stored
         * in `parameters` and a summary of the global optimization can be found in `summary`.
         * 
         * @param options pallas::SimulatedAnnealing::Options. Options used to configure the optimization.
         * @param problem pallas::GradientProblem. The problem to optimize.
         * @param parameters double*. The starting point for further optimization.
         * @param summary SimulatedAnnealing::Summary*. Summary instance to store the optimization details.
         */
        void Solve(const SimulatedAnnealing::Options& options,
                   const GradientProblem& problem,
                   double* parameters,
                   SimulatedAnnealing::Summary* global_summary);

    private:
        /**
         * @brief Checks to see if any termination conditions were met.
         * 
         * @param options pallas::SimulatedAnnealing::Options. Options used to configure the optimization.
         * @param message std::string*. If a termination condition is met, a message describing the satisfied condition is stored in the variable.
         * @param termination_type pallas::TerminationType*. This
         * @return Returns `true` if a termination condition was meet, `false` otherwise. 
         */
        bool check_for_termination_(const SimulatedAnnealing::Options& options,
                                    std::string *message,
                                    TerminationType * termination_type);

        /**
         * @brief Updates the global summary before exiting the simulated annealing algorithm.
         */
        void prepare_final_summary_(SimulatedAnnealing::Summary *global_summary,
                                    const GradientLocalMinimizer::Summary &local_summary);

        scoped_ptr<CoolingSchedule> cooling_schedule_;/**<Responsible for updating the temperature of the system. Higher temperatures make accepting a worse candidate solution more likely.*/

        internal::Metropolis metropolis_;/**<Determines whether to accept a higher cost candidate solution*/
        internal::State current_state_;/**<The current state of the optimization*/
        internal::State candidate_state_;/**<A randomized candidate solution that is then minimized using a local minimization algorithm and compared to the current state.*/
        internal::State global_minimum_state_;/**The best solution found during any stage of the optimization. It is this value that is returned when minimization concludes.*/

        unsigned int num_iterations_;/**<The number of local optimization iterations the global optimizer has performed.*/
        unsigned int num_stagnant_iterations_;/**<The number of iterations that have elapsed without finding a new global minimum.*/
    };

    /**
     * @brief Helper function that avoids going through the interface of the pallas::SimulatedAnnealing class.
     * @details The specified options are used to setup a simulated instance which
     * is then used to minimize the GradientProblem. The optimal solution is stored
     * in `parameters` and a summary of the global optimization can be found in `summary`.
     * 
     * @param options pallas::SimulatedAnnealing::Options. Options used to configure the optimization.
     * @param problem pallas::GradientProblem. The problem to optimize.
     * @param parameters double*. The starting point for further optimization.
     * @param summary SimulatedAnnealing::Summary*. Summary instance to store the optimization details.
     */
    void Solve(const SimulatedAnnealing::Options& options,
               const GradientProblem& problem,
               double* parameters,
               SimulatedAnnealing::Summary* summary);

    /**
     * @brief Dumps the system state contained in the history output into the stream contained by the writer.
     *
     * @param h SimulatedAnnealing::HistoryOutput. State of the system for a specific iteration.
     * @param writer HistoryWriter. Object responsible for writing the history output to a stream.
     */
    void dump(const SimulatedAnnealing::HistoryOutput &h, HistoryWriter& writer);

} // namespace pallas

#endif // PALLAS_SIMULATED_ANNEALING_H