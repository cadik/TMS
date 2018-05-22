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

#ifndef PALLAS_STEP_FUNCTOR_H
#define PALLAS_STEP_FUNCTOR_H

#include "pallas/types.h"
#include "pallas/internal/random_number_generator.h"
#include "pallas/scoped_ptr.h"

namespace pallas {

    /**
     * @brief Interface to produce randomized candidate solutions.
     * @details Candidate solutions are produced by randomly perturbing degrees of freedom
     * by a prescribed amount (and often limited by a maximum step size and/or upper and lower bounds).
     */
    class StepFunction {
    public:
        /**
         * @brief Destructor
         */
        virtual ~StepFunction() {}

        /**
         * @brief Produces a randomized candidate solution by modifying the input variable `x`.
         * @details Must override in derived class to implement custom step functions.
         * 
         * @param x double*. The current state of the parameters. Should be modified in place to achieve new state.
         * @param num_parameters. unsigned int. The number of parameters to modify in `x`.
         */
        virtual void Step(double* x, unsigned int num_parameters) = 0;
    };

    
    /**
     * @brief Simple candidate generator that modifies the input by a random amount between +/- `step_size`.
     * 
     * @param step_size double. Maximum magnitude a degree of freedom can be changed from its current position.
     */
    class DefaultStepFunction : public StepFunction {
    public:
        /**
         * @brief Constructor
         * 
         * @param step_size double. Maximum magnitude a degree of freedom can be changed from its current position.
         */
        DefaultStepFunction(double step_size);

        /**
         * @brief Modifies `x` in place to a new random position.
         * @details The maximum difference between the old and new position for
         * each degree of freedom is determined by the `step_size`.
         * 
         * @param x double*. The current state of the parameters. Modified in place to achieve new state.
         * @param num_parameters. unsigned int. The number of parameters to modify in `x`.
         */
        void Step(double* x, unsigned int num_parameters);

    private:
        /**
         * @brief Generates random numbers between +/- `step_size`.
         */
        scoped_ptr<internal::RandomNumberGenerator<double>> random_number_;
    };

    /**
     * @brief A new candidate solution is generated between upper and lower bounds. Each degree of freedom in the
     *  candidate solution is at most `step_size` from the current solution.
     * 
     * @param step_size double. Maximum magnitude a degree of freedom can be changed from its current position.
     * @param upper_bounds double*. Pointer to the upper bounds on the candidate's variables.
     * @param lower_bounds double*. Pointer to the lower bounds on the candidate's variables.
     * @param num_parameters unsigned int. The number of parameters to modify in `x`.
     */
    class BoundedStepFunction : public StepFunction {

    public:
        /**
         * @brief Constructor
         * 
         * @param step_size double. Maximum magnitude a degree of freedom can be changed from its current position.
         * @param upper_bounds double*. Pointer to the upper bounds on the candidate's variables.
         * @param lower_bounds double*. Pointer to the lower bounds on the candidate's variables.
         * @param num_parameters unsigned int. The number of parameters to modify in `x`. Should also be the number
         * of constraints held within the upper and lower bounds.
         */
        BoundedStepFunction(double step_size,
                            const double* upper_bounds,
                            const double* lower_bounds,
                            unsigned int num_parameters);

        /**
         * @brief Modifies `x` in place to a new random position.
         * @details The maximum difference between the old and new position for
         * each degree of freedom is determined by the `step_size`. If any degree
         * of freedom exceeds the upper or lower bounds, it is adjusted to fall
         * within the boundary.
         * 
         * @param x double*. The current state of the parameters. Modified in place to achieve new state.
         * @param num_parameters. unsigned int. The number of parameters to modify in `x`.
         */
        void Step(double* x, unsigned int num_parameters);

    private:
        /**
         * @brief Generates random numbers between +/- `step_size`.
         */
        scoped_ptr<internal::RandomNumberGenerator<double>> random_number_;

        Vector upper_bounds_;/**<Upper bounds on candidate solutions*/
        Vector lower_bounds_;/**<Lower bounds on candidate solutions*/
    };
} // namespace pallas

#endif //PALLAS_STEP_FUNCTOR_H
