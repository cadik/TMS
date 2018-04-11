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

#include "pallas/internal/mutation_strategy.h"
#include "glog/logging.h"

namespace pallas {
    namespace internal {

        MutationStrategy* MutationStrategy::Create(MutationStrategyType type) {
            switch(type) {
                case BEST_1:
                    return new MutateBest1();
                case RAND_1:
                    return new MutateRand1();
                case RAND_TO_BEST_1:
                    return new MutateRandToBest1();
                case BEST_2:
                    return new MutateBest2();
                case RAND_2:
                    return new MutateRand2();
                default:
                    LOG(ERROR) << "Unknown mutation strategy type: " << MutationStrategyTypeToString(type);
                    return NULL;
            }
        };

        unsigned int MutationStrategy::NumSamples() const {
            return num_samples_;
        };

        MutateBest1::MutateBest1() {
            num_samples_ = 2;
        };

        Vector MutateBest1::get_bprime(const std::vector<Vector> &population,
                                       unsigned int candidate,
                                       const int* samples,
                                       double scale) {
            const int r0 = samples[0];
            const int r1 = samples[1];

            return population[0] + scale * (population[r0] - population[r1]);
        };

        MutateRand1::MutateRand1() {
            num_samples_ = 3;
        };

        Vector MutateRand1::get_bprime(const std::vector<Vector> &population,
                                       unsigned int candidate,
                                       const int* samples,
                                       double scale) {
            const int r0 = samples[0];
            const int r1 = samples[1];
            const int r2 = samples[2];

            return population[r0] + scale * (population[r1] - population[r2]);
        };

        MutateRandToBest1::MutateRandToBest1() {
            num_samples_ = 2;
        };

        Vector MutateRandToBest1::get_bprime(const std::vector<Vector> &population,
                                             unsigned int candidate,
                                             const int* samples,
                                             double scale) {
            const int r0 = samples[0];
            const int r1 = samples[1];
            Vector bprime = population[candidate];
            bprime += scale * (population[0] - bprime);
            return scale * (population[r0] - population[r1]);
        };

        MutateBest2::MutateBest2() {
            num_samples_ = 4;
        };

        Vector MutateBest2::get_bprime(const std::vector<Vector> &population,
                                       unsigned int candidate,
                                       const int* samples,
                                       double scale) {
            const int r0 = samples[0];
            const int r1 = samples[1];
            const int r2 = samples[2];
            const int r3 = samples[3];

            return population[0] + scale * (population[r0] + population[r1] - population[r2] - population[r3]);
        };

        MutateRand2::MutateRand2() {
            num_samples_ = 5;
        };

        Vector MutateRand2::get_bprime(const std::vector<Vector> &population,
                                       unsigned int candidate,
                                       const int* samples,
                                       double scale) {
            const int r0 = samples[0];
            const int r1 = samples[1];
            const int r2 = samples[2];
            const int r3 = samples[3];
            const int r4 = samples[4];

            return population[r0] + scale * (population[r1] + population[r2] - population[r3] - population[r4]);
        };

    } // namespace internal
} // namespace pallas