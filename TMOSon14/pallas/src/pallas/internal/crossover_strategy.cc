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

#include "pallas/internal/crossover_strategy.h"

namespace pallas {
    namespace internal {

        CrossoverStrategy* CrossoverStrategy::Create(CrossoverStrategyType type,
                                                     double crossover_probability,
                                                     unsigned int num_parameters) {
            switch (type) {
                case BINOMIAL:
                    return new BinomialCrossover(crossover_probability, num_parameters);
                case EXPONENTIAL:
                    return new ExponentialCrossover(crossover_probability, num_parameters);
                default:
                    LOG(ERROR) << "Unknown crossover strategy type: " << CrossoverStrategyTypeToString(type);
                    return NULL;
            }
        };

        void CrossoverStrategy::set_crossover_probability(double crossover_probability) {
            crossover_probability_ = crossover_probability;
        };

        BinomialCrossover::BinomialCrossover(double crossover_probability, unsigned int num_parameters) {
            scoped_ptr <RandomNumberGenerator<double>> tmp_rng_dbl(new RandomNumberGenerator<double>());
            swap(random_double_, tmp_rng_dbl);

            scoped_ptr <RandomNumberGenerator<unsigned int>> tmp_rng_uint(
                    new RandomNumberGenerator <unsigned int>(0, num_parameters - 1));
            swap(random_uint_, tmp_rng_uint);

            crossover_probability_ = crossover_probability;
        };

        void BinomialCrossover::Crossover(Vector &trial, const Vector &bprime) {
            for (unsigned int i = 0; i < trial.size() - 1; ++i) {
                if ((*random_double_)() < crossover_probability_) {
                    trial[i] = bprime[i];
                }
            }

            // one random entry in trial will always have a crossover
            const unsigned int fill_point = (*random_uint_)();
            trial[fill_point] = bprime[fill_point];
        };

        ExponentialCrossover::ExponentialCrossover(double crossover_probability, unsigned int num_parameters) {
            scoped_ptr <RandomNumberGenerator<double>> tmp_rng_dbl(new RandomNumberGenerator<double>());
            swap(random_double_, tmp_rng_dbl);


            scoped_ptr <RandomNumberGenerator<unsigned int>> tmp_rng_uint(
                    new RandomNumberGenerator <unsigned int>(0, num_parameters - 1));
            swap(random_uint_, tmp_rng_uint);

            crossover_probability_ = crossover_probability;
        };

        void ExponentialCrossover::Crossover(Vector &trial, const Vector &bprime) {
            unsigned int i = 0;
            const unsigned int num_parameters = trial.size();
            unsigned int fill_point = (*random_uint_)();

            while (i < num_parameters && (*random_double_)() < crossover_probability_) {
                trial[fill_point] = bprime[fill_point];
                fill_point = (fill_point + 1) % num_parameters;
                ++i;
            }
        };

    } // namespace internal
} // namespace pallas