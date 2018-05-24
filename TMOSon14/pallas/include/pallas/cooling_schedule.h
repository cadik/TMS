/*!
* \file cooling_schedule.h
*
* \author Ryan Latture
* \date 5-18-15
*
* This file contains a C++ implementation of a cooling schedules used in the pallas::SimulatedAnnealing algorithm.
* This code is based on the SciPy implementation of simulated annealing found in scipy.optimize.
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

#ifndef PALLAS_COOLING_SCHEDULE_H
#define PALLAS_COOLING_SCHEDULE_H

#include "pallas/types.h"
#include "pallas/internal/state.h"
#include "pallas/step_function.h"


namespace pallas {

    /**
     * @brief Determines the rate of cooling within pallas::SimulatedAnnealing.
     */
    class CoolingSchedule {
    public:
        /**
         * Configurable options for modifying the default behaviour of the cooling schedule.
         */
        struct Options {
            /**
             * @brief Default constructor
             * @details Allows the type of cooling schedule and its respective parameters to be set as well as the initial temperature to start the simulation.
             */
            Options() {
                type = BOLTZMANN;
                initial_temperature = -1.0;
                final_temperature = 1.0e-12;
                boltzmann_constant = 1.0;
                fast_quench_param = 1.0;
                fast_m_param = 1.0;
                fast_n_param = 1.0;
            }

            CoolingScheduleType type;/**<Algorithm to use to decrease temperature. Can be one of `BOLTZMANN`, `CAUCHY`, or `FAST`.*/
            double initial_temperature;/**<Temperature to start the minimization at. Higher temperatures mean that it is more likely that algorithm accepts worse candidate solutions relative to its current state.*/
            double final_temperature;/**<The final temperature to cool to. Optimization terminates if the temperature reaches this value.*/
            double boltzmann_constant;/**<Boltzmann constant in the probabilistic acceptance criteria (increase for less stringent criteria at each temperature).*/
            double fast_quench_param;/**<Parameter to alter the fast simulated annealing schedule. See pallas::FastCooling.*/
            double fast_m_param;/**<Parameter to alter the fast simulated annealing schedule. See pallas::FastCooling.*/
            double fast_n_param;/**<Parameter to alter the fast simulated annealing schedule. See pallas::FastCooling.*/

        };

        /**
         * @brief Creates a pointer to a cooling schedule.
         * @details Based on the options the proper cooling schedule, with correct parameters, is returned from this function.
         * 
         * @param options pallas::CoolingSchedule::Options. Configured options to specify the cooling schedule behaviour.
         * @return pallas::CoolingSchedule*. A pointer to a cooling schedule.
         */
        static CoolingSchedule* Create(const Options& options);

        /**
         * @brief Default destructor
         */
        virtual ~CoolingSchedule() {};

        /**
         * @brief Interface to update the temperature member variable.
         * @details Should override in derived cooling schedules.
         */
        virtual void update_temperature() = 0;

        /**
         * @brief Returns the current temperature stored withing the cooling schedule.
         */
        double get_temperature() const;
        /**
         * @brief Returns the initial temperature to start the minimization at.
         */
        double get_initial_temperature() const;

        /**
         * @brief Sets the current temperature.
         * 
         * @param T double. The current temperature to set the cooling schedule at.
         */
        void set_temperature(double T);

        /**
         * @brief Estimates an appropriate starting temperature.
         * @details If no starting temperature is supplied by the user, this method attempts to 
         * approximate an acceptable starting temperature by sampling random points about the starting point
         * in order to gauge the magnitude of the cost function.
         * 
         * @param problem pallas::GradientProblem. The problem to optimize.
         * @param state pallas::internal::State. The details of the current state of the minimization algorithm.
         * @param step_function pallas::StepFunction. This function produces randomized candidate solutions.
         */
        void calc_start_temperature(const GradientProblem& problem,
                                    internal::State& state,
                                    StepFunction* step_function);


    protected:
        double temperature;/**<Current temperature of the cooling schedule.*/
        double initial_temperature;/**Initial temperature the cooling schedule should begin at.*/
        double boltzmann_constant;/**<Boltzmann constant in the probabilistic acceptance criteria (increase for less stringent criteria at each temperature).*/

    };

    /**
     * @brief Fast cooling schedule updates.
     * @details Updates the temperature according to the formula:
     * \code
     T = pow(-c * k, quench);
     T = T0 * exp(T);
     ++k;
     * \endcode
     * where `c = m * exp(-n * quench)` and `k = boltzmann_constant`.
     * @param options Configurable options for modifying the default behaviour of the cooling schedule. Can be used to the m, n, k, and quench parameters for the cooling schedule.
     */
    class FastCooling : public CoolingSchedule {
    public:
        /**
         * @brief Default constructor
         * @param options 
         */
        FastCooling(const CoolingSchedule::Options& options);

        /**
         * @brief Updates the `temperature` member variable according to the cooling schedule method.
         */
        void update_temperature();

    private:
        double fast_m_param;/**<See class description for details. Increase for faster cooling.*/
        double fast_n_param;/**<See class description for details. Increase for slower cooling.*/
        double fast_c_param;/**<See class description for details. Calculated from m, n, and quench parameters.*/
        double fast_quench_param;/**<See class description for details. Increase for faster cooling.*/
    };

    /**
     * @brief Cauchy cooling schedule.
     *
     * Temperature updates are done using:
     * \code
     T = T0 / (1.0 + k);
     k++;
     * \endcode
     * where `k= boltzmann_constant`.
     */
    class CauchyCooling : public CoolingSchedule {
    public:
        CauchyCooling(const CoolingSchedule::Options& options);

        /**
         * @brief Updates the `temperature` member variable according to the cooling schedule method.
         */
        void update_temperature();
    };

    /**
     * @brief Boltzmann cooling schedule.
     *
     * Temperature updates are done using:
     * \code
     k++;
     T = T0 / log(k + 1.0);
     * \endcode
     * where `k = boltzmann_constant`.
     */
    class BoltzmannCooling : public CoolingSchedule {
    public:
        /**
         * Default constructor
         */
        BoltzmannCooling(const CoolingSchedule::Options& options);

        /**
         * @brief Updates the `temperature` member variable according to the cooling schedule method.
         */
        void update_temperature();
    };

} // namespace pallas

#endif // PALLAS_COOLING_SCHEDULE_H
