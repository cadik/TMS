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

#include "pallas/step_function.h"

namespace pallas {

    DefaultStepFunction::DefaultStepFunction(double step_size) :
            random_number_(new internal::RandomNumberGenerator<double>(-step_size, step_size)) {

    };

    void DefaultStepFunction::Step(double* x, unsigned int num_parameters) {
        for(unsigned int i =0; i < num_parameters; ++i) {
            x[i] += (*random_number_)();
        }
    };

    BoundedStepFunction::BoundedStepFunction(double step_size,
                                             const double *upper_bounds,
                                             const double *lower_bounds,
                                             unsigned int num_parameters)
            : random_number_(new internal::RandomNumberGenerator<double>(-step_size, step_size)) {
        ConstVectorRef upper_ref(upper_bounds, num_parameters);
        ConstVectorRef lower_ref(lower_bounds, num_parameters);

        upper_bounds_ = upper_ref;
        lower_bounds_ = lower_ref;
    };

    void BoundedStepFunction::Step(double *x, unsigned int num_parameters) {

        double candidate_x;
        for(unsigned int i =0; i < num_parameters; ++i) {
            candidate_x = x[i] + (*random_number_)();

            if (candidate_x > upper_bounds_(i))
                x[i] = 2.0 * upper_bounds_(i) - candidate_x;
            else if (candidate_x < lower_bounds_(i))
                x[i] = 2.0 * lower_bounds_(i) - candidate_x;
            else
                x[i] = candidate_x;
        }
    };

} // namespace pallas
