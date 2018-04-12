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
       

        /* if (parameters[0] < 0.5) {
            return false;
        }*/
        for (int i = 0; i < 500000; i++) {
            cost[0] += parameters[i]*i;
        }
       /* if (gradient != NULL) {
            gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
            gradient[1] = 200.0 * (y - x * x);
        }*/
        return true;
    }

    virtual int NumParameters() const { return 500000; }
};

int main(int argc, char** argv) {
    google::InitGoogleLogging(argv[0]);

    // define the starting point for the optimization
    double parameters[500000];

    for (int i = 0; i < 500000; i++) {
        parameters[i] = 1.0;
    }

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
    double upper_bounds [250000];
    std::fill_n(upper_bounds, 250000, 5);
    double lower_bounds [250000];
    std::fill_n(lower_bounds, 250000, -5);

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