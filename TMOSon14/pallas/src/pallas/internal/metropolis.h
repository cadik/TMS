/*!
* \file Metropolis.h
*
* \author Ryan Latture
* \date 8-17-14
*
* Metropolis acceptance/rejection criterion.
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

#ifndef PALLAS_METROPOLIS_H
#define PALLAS_METROPOLIS_H

#include <math.h>

#include "pallas/internal/random_number_generator.h"

namespace pallas {
    namespace internal {
        /**
        * @brief Implements a probabilistic acceptance criterion for candidate solutions.
        */
        class Metropolis {
            double beta;
            RandomNumberGenerator<double> *random_num;/**<generates random numbers to determine if candidate solution is accepted.*/

        public:

            Metropolis() {
                beta = 1.0;
                random_num = new RandomNumberGenerator<double>();
            };

            /**Constructor
            * @param T double. "Temperature" --higher temperatures will result in higher cost functions being accepted with higher probability.
            */
            Metropolis(double T) {
                beta = 1.0 / T;
                random_num = new RandomNumberGenerator<double>();
            };

            /**Default destructor*/
            ~Metropolis() {
                delete random_num;
            };

            /**
            * @brief Returns a bool whether to accept the candidate solution.
            *
            * Accepts the candidate solution based on the function:
            * /code
            * 	double w = std::min(1.0, exp(-1.0*(cost_new - cost_old)*beta));
            *	double r = (*random_num)();
            *	return w >= r;
            * /endcode
            *
            * @param cost_new double. cost of the candidate solution.
            * @param cost_old cost of the current solution.
            *
            * @return <B>Bool</B> whether the candidate solution was accepted.
            */
            bool accept_reject(double cost_new, double cost_old) {
                double w = std::min(1.0, exp(-1.0 * (cost_new - cost_old) * beta));
                double r = (*random_num)();
                return w >= r;
            };

            /** @brief Calls the accept_reject function
            * Accepts the candidate solution based on the function:
            * /code
            * 	double w = std::min(1.0, exp(-1.0*(cost_new - cost_old)*beta));
            *	double r = (*random_num)();
            *	return w >= r;
            * /endcode
            * @param cost_new double. cost of the candidate solution.
            * @param cost_old cost of the current solution.
            *
            * @return <B>Bool</B> whether the candidate solution was accepted.
            */
            inline bool operator()(double cost_new, double cost_old) {
                return accept_reject(cost_new, cost_old);
            };
        };
    } // namespace internal
} // namespace pallas

#endif //METROPOLIS_H