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

#include <cctype>
#include <string>
#include "pallas/types.h"

namespace pallas {

    using std::string;

#define CASESTR(x) case x: return #x
#define STRENUM(x) if (value == #x) { *type = x; return true;}

    static void UpperCase(string* input) {
        std::transform(input->begin(), input->end(), input->begin(), ::toupper);
    }

    const char* CoolingScheduleTypeToString(CoolingScheduleType type) {
        switch (type) {
            CASESTR(BOLTZMANN);
            CASESTR(CAUCHY);
            CASESTR(FAST);
            default:
                return "UNKNOWN";
        }
    }

    bool StringToCoolingScheduleType(string value, CoolingScheduleType* type) {
        UpperCase(&value);
        STRENUM(BOLTZMANN);
        STRENUM(CAUCHY);
        STRENUM(FAST);
        return false;
    }

    const char* MutationStrategyTypeToString(MutationStrategyType type) {
        switch (type) {
            CASESTR(BEST_1);
            CASESTR(RAND_1);
            CASESTR(RAND_TO_BEST_1);
            CASESTR(BEST_2);
            CASESTR(RAND_2);
            default:
                return "UNKNOWN";
        }
    }

    bool StringToMutationStrategyType(std::string value, MutationStrategyType* type) {
        UpperCase(&value);
        STRENUM(BEST_1);
        STRENUM(RAND_1);
        STRENUM(RAND_TO_BEST_1);
        STRENUM(BEST_2);
        STRENUM(RAND_2);
        return false;
    }

    const char* CrossoverStrategyTypeToString(CrossoverStrategyType type) {
        switch (type) {
            CASESTR(BINOMIAL);
            CASESTR(EXPONENTIAL);
            default:
                return "UNKNOWN";
        }
    }

    bool StringToCrossoverStrategyType(std::string value, CrossoverStrategyType* type) {
        UpperCase(&value);
        STRENUM(BINOMIAL);
        STRENUM(EXPONENTIAL);
        return false;
    }

    const char* PopulationInitializationTypeToString(PopulationInitializationType type) {
        switch (type) {
            CASESTR(LATIN_HYPERCUBE);
            CASESTR(RANDOM);
            default:
                return "UNKNOWN";
        }
    }

    bool StringToPopulationInitializationType(std::string value, PopulationInitializationType* type) {
        UpperCase(&value);
        STRENUM(LATIN_HYPERCUBE);
        STRENUM(RANDOM);
        return false;
    }

} // namespace pallas
