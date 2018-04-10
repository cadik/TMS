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

#include <cfloat>
#include <cmath>
#include "glog/logging.h"
#include "pallas/cooling_schedule.h"

namespace pallas {

    CoolingSchedule* CoolingSchedule::Create(const CoolingSchedule::Options& options) {
        
        switch(options.type) {
            case FAST:
                return new FastCooling(options);

            case CAUCHY:
                return new CauchyCooling(options);

            case BOLTZMANN:
                return new BoltzmannCooling(options);

            default:
                LOG(ERROR) << "Unknown cooling schedule type: " << options.type;
                return NULL;
        }
    }

    double CoolingSchedule::get_temperature() const {
        return temperature;
    };

    double CoolingSchedule::get_initial_temperature() const {
        return initial_temperature;
    };

    void CoolingSchedule::set_temperature(double T) {
        temperature = T;
    };

    void CoolingSchedule::calc_start_temperature(const GradientProblem &problem,
                                                 internal::State& best_state,
                                                 StepFunction* step_function) {
        double cost_max = DBL_MIN;
        double cost_min = DBL_MAX;
        const unsigned int num_parameters = static_cast<unsigned int>(best_state.x.size());
        internal::State state = best_state; 

        for(unsigned int i = 0; i < 100; ++i) {
            step_function->Step(state.x.data(), num_parameters);
            if(!problem.Evaluate(state.x.data(), &state.cost, state.gradient.data())) {
                LOG(ERROR) << "Problem evaluation failed at:\n" << state.x << "\n in CoolingSchedule::calc_start_temperature.";
            }

            if (state.cost > cost_max) cost_max = state.cost;
            if (state.cost < cost_min) {
                cost_min = state.cost;
                best_state.update(state);
            }
        } 

        initial_temperature = (cost_max - cost_min) * 1.2;
        temperature = initial_temperature;
    };

    FastCooling::FastCooling(const CoolingSchedule::Options& options) {
        temperature = options.initial_temperature;
        initial_temperature = options.initial_temperature;
        fast_m_param = options.fast_m_param;
        fast_n_param = options.fast_n_param;
        fast_quench_param = options.fast_quench_param;
        boltzmann_constant = options.boltzmann_constant;
        fast_c_param = fast_m_param * std::exp(-fast_n_param * fast_quench_param);
    }

    void FastCooling::update_temperature() {
        temperature = std::pow(-fast_c_param * boltzmann_constant, fast_quench_param);
        temperature = initial_temperature * std::exp(temperature);
        ++boltzmann_constant;
    };

    CauchyCooling::CauchyCooling(const CoolingSchedule::Options& options) {
        temperature = options.initial_temperature;
        initial_temperature = options.initial_temperature;
        boltzmann_constant = options.boltzmann_constant;

    };

    void CauchyCooling::update_temperature() {
        temperature = initial_temperature / (1.0 + boltzmann_constant);;
        ++boltzmann_constant;
    };

    BoltzmannCooling::BoltzmannCooling(const CoolingSchedule::Options& options) {
        temperature = options.initial_temperature;
        initial_temperature = options.initial_temperature;
        boltzmann_constant = options.boltzmann_constant;
    };

    void BoltzmannCooling::update_temperature() {
        ++boltzmann_constant;
        temperature = initial_temperature / std::log(1.0 + boltzmann_constant);
    };

} // namespace pallas
