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

#ifndef PALLAS_MUTATION_STRATEGY_H
#define PALLAS_MUTATION_STRATEGY_H

#include "pallas/types.h"

namespace pallas {
    namespace internal {

        class MutationStrategy {
        public:
            static MutationStrategy* Create(MutationStrategyType type);

            virtual ~MutationStrategy() {};

            virtual Vector get_bprime(const std::vector<Vector> &population,
                                      unsigned int candidate,
                                      const int* samples,
                                      double scale) = 0;

            unsigned int NumSamples() const;

        protected:
            unsigned int num_samples_;
        };

        class MutateBest1 : public MutationStrategy {
        public:
            MutateBest1();

            Vector get_bprime(const std::vector<Vector> &population,
                              unsigned int candidate,
                              const int* samples,
                              double scale);
        };

        class MutateRand1 : public MutationStrategy {
        public:
            MutateRand1();

            Vector get_bprime(const std::vector<Vector> &population,
                              unsigned int candidate,
                              const int* samples,
                              double scale);
        };

        class MutateRandToBest1 : public MutationStrategy {
        public:
            MutateRandToBest1();

            Vector get_bprime(const std::vector<Vector> &population,
                              unsigned int candidate,
                              const int* samples,
                              double scale);

        };

        class MutateBest2 : public MutationStrategy {
        public:
            MutateBest2();

            Vector get_bprime(const std::vector<Vector> &population,
                              unsigned int candidate,
                              const int* samples,
                              double scale);

        };

        class MutateRand2 : public MutationStrategy {
        public:
            MutateRand2();

            Vector get_bprime(const std::vector<Vector> &population,
                              unsigned int candidate,
                              const int* samples,
                              double scale);

        };

    } // namespace internal
} // namespace pallas

#endif // PALLAS_MUTATION_STATEGY_H
