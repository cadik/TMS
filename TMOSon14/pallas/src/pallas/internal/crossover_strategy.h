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

#ifndef PALLAS_CROSSOVER_STRATEGY_H
#define PALLAS_CROSSOVER_STRATEGY_H

#include "pallas/types.h"
#include "pallas/scoped_ptr.h"
#include "pallas/internal/random_number_generator.h"

namespace pallas {
    namespace internal {

        class CrossoverStrategy {
        public:
            static CrossoverStrategy* Create(CrossoverStrategyType type,
                                             double crossover_probability,
                                             unsigned int num_parameters);

            virtual ~CrossoverStrategy() { };

            virtual void Crossover(Vector &trial, const Vector &bprime) = 0;

            void set_crossover_probability(double crossover_probability);

        protected:
            scoped_ptr <RandomNumberGenerator<double>> random_double_;
            scoped_ptr <RandomNumberGenerator<unsigned int>> random_uint_;
            double crossover_probability_;
        };

        class BinomialCrossover : public CrossoverStrategy {
        public:
            BinomialCrossover(double crossover_probability, unsigned int num_parameters);

            void Crossover(Vector &trial, const Vector &bprime);
        };

        class ExponentialCrossover : public CrossoverStrategy {
        public:
            ExponentialCrossover(double crossover_probability, unsigned int num_parameters);

            void Crossover(Vector &trial, const Vector &bprime);
        };

    } // namespace internal
} // namespace pallas

#endif // PALLAS_CROSSOVER_STRATEGY_H
