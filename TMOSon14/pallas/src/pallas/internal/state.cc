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

#include "pallas/internal/state.h"

namespace pallas {
    namespace internal {

        State::State(unsigned int num_parameters)
                : cost(0),
                  x(num_parameters),
                  gradient(num_parameters),
                  gradient_squared_norm(0.0),
                  gradient_max_norm(0.0),
                  directional_derivative(0.0),
                  tolerance(1e-10) {

        };

        State::State(const State& s) {
                clone_(s);
        }

        bool State::update(const State& s) {
            if (s.cost + tolerance < cost) {
                clone_(s);
                return true;
            }

            return false;
        }

        void State::clone_(const State& s) {
            cost = s.cost;
            x = s.x;
            gradient = s.gradient;
            gradient_squared_norm = s.gradient_squared_norm;
            gradient_max_norm = s.gradient_max_norm;
            directional_derivative = s.directional_derivative;
        }

    }
}
