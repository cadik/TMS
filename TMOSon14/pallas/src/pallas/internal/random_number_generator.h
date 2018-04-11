/*!
* \file random_number_generator.h
*
* \author Ryan Latture
* \date 8-17-14
*
* Class that generates random numbers
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

#ifndef PALLAS_RANDOM_NUMBER_GENERATOR_H
#define PALLAS_RANDOM_NUMBER_GENERATOR_H

#include <random>
#include <type_traits>

namespace pallas {
    namespace internal {
        /**
        * @brief Generates random numbers.
        * @details The engine specified by the typedef Engine is used to generate random numbers with a distribution defined by the typedef Distribution.
        *
        */
        template <class T, class Enable = void> class RandomNumberGenerator;

        template<class T>
        class RandomNumberGenerator<T, typename std::enable_if<std::is_integral<T>::value >::type> {
            typedef std::mt19937 Engine;
            typedef std::uniform_int_distribution<T> Distribution;

        private:
            Engine *_eng;/**<Pointer to an Engine*/
            Distribution *_dist;/**<Pointer to a Distribution*/
            std::random_device rd;/**<Random device used to seed the Engine*/

        public:
            /**
            * @brief Constructor with min and max values specified.
            * @details Produces numbers on the range `[minVal, maxVal)` seeded by a random device.
            *
            * @param minVal T. The minimum value that the random number can output (inclusive).
            * @param maxVal T. The random number generator will output values lower than maxVal (exclusive).
            */
            RandomNumberGenerator(T minVal, T maxVal) {
                std::random_device rd;
                _eng = new Engine(rd());
                _dist = new Distribution(minVal, maxVal);
            };

            /**
            * @brief Constructor with min and max values specified and explicit integer seed.
            * @details Produces numbers on the range `[minVal, maxVal)` seeded by the specified integer seed.
            * Use explicit seed if predictable outputs from the random number generator are needed.
            *
            * @param minVal T. The minimum value that the random number can output (inclusive).
            * @param maxVal T. The random number generator will output values lower than maxVal (exclusive).
            * @param seed int. A seeding value. An entire state sequence is generated from this value using
            * a linear random generator.
            */
            RandomNumberGenerator(T minVal, T maxVal, int seed) {
                std::random_device rd;
                _eng = new Engine(seed);
                _dist = new Distribution(minVal, maxVal);
            };

            /**
            * @brief Default constructor
            * @details constructs a random number generator on the interval [0.0,1.0) seeded with a random device.
            */
            RandomNumberGenerator() {
                std::random_device rd;
                _eng = new Engine(rd());
                _dist = new Distribution(0, 1);
            };

            /**
            * @brief Constructor with default range '[0, 1)' and explicit integer seed.
            * @details Produces numbers on the range `[0, 1)` seeded by the specified integer seed.
            * Use explicit seed if predictable outputs from the random number generator are needed.
            *
            * @param seed int. A seeding value. An entire state sequence is generated from this value using
            * a linear random generator.
            */
            RandomNumberGenerator(int seed) {
                std::random_device rd;
                _eng = new Engine(seed);
                _dist = new Distribution(0, 1);
            };

            ~RandomNumberGenerator() {
                delete _dist;
                delete _eng;
            };

            /**
            * @brief Returns a random number
            *
            * @return Returns a number between the minimum and maximum values specified when constructing the distribution.
            */
            inline T operator()() {
                return (*_dist)(*_eng);
            };
        };

        template<class T>
        class RandomNumberGenerator<T, typename std::enable_if<std::is_floating_point<T>::value >::type> {
            typedef std::mt19937 Engine;
            typedef std::uniform_real_distribution<T> Distribution;

        private:
            Engine *_eng;/**<Pointer to an Engine*/
            Distribution *_dist;/**<Pointer to a Distribution*/
            std::random_device rd;/**<Random device used to seed the Engine*/

        public:
            /**
            * @brief Constructor with min and max values specified.
            * @details Produces numbers on the range `[minVal, maxVal)` seeded by a random device.
            *
            * @param minVal T. The minimum value that the random number can output (inclusive).
            * @param maxVal T. The random number generator will output values lower than maxVal (exclusive).
            */
            RandomNumberGenerator(T minVal, T maxVal) {
                std::random_device rd;
                _eng = new Engine(rd());
                _dist = new Distribution(minVal, maxVal);
            };

            /**
            * @brief Constructor with min and max values specified and explicit integer seed.
            * @details Produces numbers on the range `[minVal, maxVal)` seeded by the specified integer seed.
            * Use explicit seed if predictable outputs from the random number generator are needed.
            *
            * @param minVal T. The minimum value that the random number can output (inclusive).
            * @param maxVal T. The random number generator will output values lower than maxVal (exclusive).
            * @param seed int. A seeding value. An entire state sequence is generated from this value using
            * a linear random generator.
            */
            RandomNumberGenerator(T minVal, T maxVal, int seed) {
                std::random_device rd;
                _eng = new Engine(seed);
                _dist = new Distribution(minVal, maxVal);
            };

            /**
            * @brief Default constructor
            * @details constructs a random number generator on the interval [0.0,1.0) seeded with a random device.
            */
            RandomNumberGenerator() {
                std::random_device rd;
                _eng = new Engine(rd());
                _dist = new Distribution(0, 1);
            };

            /**
            * @brief Constructor with default range '[0, 1)' and explicit integer seed.
            * @details Produces numbers on the range `[0, 1)` seeded by the specified integer seed.
            * Use explicit seed if predictable outputs from the random number generator are needed.
            *
            * @param seed int. A seeding value. An entire state sequence is generated from this value using
            * a linear random generator.
            */
            RandomNumberGenerator(int seed) {
                std::random_device rd;
                _eng = new Engine(seed);
                _dist = new Distribution(0, 1);
            };

            ~RandomNumberGenerator() {
                delete _dist;
                delete _eng;
            };

            /**
            * @brief Returns a random number
            *
            * @return Returns a number between the minimum and maximum values specified when constructing the distribution.
            */
            inline T operator()() {
                return (*_dist)(*_eng);
            };
        };

    } // namespace internal
} // namespace pallas

#endif //PALLAS_RANDOM_NUMBER_GENERATOR_H