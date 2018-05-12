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

#ifndef PALLAS_INTERNAL_SHUFFLER_H
#define PALLAS_INTERNAL_SHUFFLER_H

#include "pallas/internal/random_number_generator.h"
#include "pallas/scoped_ptr.h"


namespace pallas {
    namespace internal {
        class Shuffler {
        public:
            Shuffler(unsigned int population_size)
                    : population_size_(population_size) {
                scoped_ptr<RandomNumberGenerator<unsigned int>> tmp_random_idx(
                        new RandomNumberGenerator<unsigned int>(0, population_size_ - 1));
                swap(random_idx_, tmp_random_idx);
            };

            template<typename T>
            void Shuffle(T* data, unsigned int num_times = 1) {
                for (unsigned int i = 0; i < num_times; ++i) {
                    for (unsigned int j = population_size_ - 1; j > 0; --j) {
                        std::swap(data[j], data[(*random_idx_)()]);
                    }
                }
            };

        private:
            const unsigned int population_size_;
            scoped_ptr<RandomNumberGenerator<unsigned int>> random_idx_;

        };

    } // namespace internal
} // namespace pallas

#endif //PALLAS_RANDOM_NUMBER_GENERATOR_H