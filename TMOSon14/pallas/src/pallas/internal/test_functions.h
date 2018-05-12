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

#ifndef PALLAS_TEST_FUNCTIONS_H
#define PALLAS_TEST_FUNCTIONS_H

#include "pallas/types.h"

#include <cmath>

#define PI           3.14159265358979323846

namespace pallas {
    namespace internal {
        using std::cos;
        using std::sin;
        using std::sqrt;
        using std::exp;
        using std::abs;

        static const double pow_LUT[2] = { 1.0, -1.0 };

        /**
         * @brief Schwefel's function
         * @details Schwefel's function is deceptive in that the global minimum is geometrically distant,
         * over the parameter space, from the next best local minima. Therefore, the search algorithms are
         * potentially prone to convergence in the wrong direction.
         *
         * <B>Global minimum:</B>
         * f(x)=0; x(i)=420.9687, i=1:n
         *
         * <B>Bounds:</B>
         * -500 < x<sub>i</sub> < 500
         *
         * @param x column_vector.
         * @return <B>cost</B>	double.
         */
        double schwefel(const Vector &x)
        {
            unsigned int num = x.size();
            double s = 0.0;
            for (unsigned int i = 0; i < num; ++i)
            {
                s += -x(i) * sin(sqrt(abs(x(i))));
            }
            return 418.9829 * num + s;
        }

        /**
         * @brief Rosenbrock's function
         * @details Rosenbrock's valley is a classic optimization problem, also known as Banana function.
         * The global optimum is inside a long, narrow, parabolic shaped flat valley. To find the valley is trivial;
         * however, convergence to the global optimum is difficult and hence this problem has been repeatedly
         * used in assess the performance of optimization algorithms.
         *
         * <B>Global minimum:</B>
         * f(x)=0; x(i)=1, i=1:n.
         *
         * <B>Bounds:</B>
         * -2.048 < x<sub>i</sub> < 2.048
         *
         * @param x column_vector.
         * @return <B>cost</B>	double.
         */
        double rosenbrock(const Vector &x)
        {
            unsigned int num = x.size();
            double s = 0.0;
            for (unsigned int i = 0; i < (num - 1); ++i)
            {
                s += (x(i+1) - x(i)*x(i))*(x(i+1) - x(i)*x(i)) + (1-x(i)) * (1-x(i));
            }
            return 100.0 * s;
        }

        /**
         * @brief Easom's function
         * @details The Easom function is a unimodal test function, where the global minimum has
         * a small area relative to the search space. The function was inverted for minimization.
         *
         * <B>Global minimum:</B>
         * f(x)=0.0; x(i)=&pi, i=1:n.
         *
         * <B>Bounds:</B>
         * -2&pi < x<sub>i</sub> < 2&pi
         *
         * @param x column_vector.
         * @return <B>cost</B>	double.
         */
        double easom(const Vector &x)
        {
            unsigned int num = x.size();
            double sm = 0.0;
            double tm = 1.0;
            for (unsigned int i = 0; i < num; ++i)
            {
                sm += (x(i) - PI)*(x(i) - PI);
                tm *= cos(x(i))*cos(x(i));
            }
            return -(pow_LUT[num % 2]) * tm * exp(-sm) + 1.0;
        }

        /**
         * @brief Xin-She Yang 2 test objective function
         * @details The Xin-She Yang 2 function is a very difficult multimodal minimization problem. Gradient
         * Solvers often become trapped within the extreme fluctuations at the edges of the bounds without converging
         * on the global minimimum at the origin.
         *
         * <B>Global minimum:</B>
         * f(x)=0; x(i)=0, i=1:n.
         *
         * <B>Bounds:</B>
         * -2&pi < x<sub>i</sub> < 2&pi
         *
         * @param x	column_vector.
         * @return <B>cost</B>	double.
         */
        double xsy02(const Vector &x)
        {
            unsigned int num = x.size();
            double sm1 = 0.0;
            double sm2 = 0.0;
            for (unsigned int i = 0; i < num; ++i)
            {
                sm1 += abs(x(i));
                sm2 += sin(x(i)*x(i));
            }
            return sm1 / exp(sm2);
        }
    }
}

#endif // PALLAS_TEST_FUNCTIONS_H